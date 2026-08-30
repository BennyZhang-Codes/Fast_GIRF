function report = check_HighOrderSequence(seqFile, manifestFile)
% Validate a generated Fast High-Order GSTF Pulseq sequence.
%
% This checker parses the written .seq file directly. It verifies the final
% scanner sequence rather than relying on variables in the generator.
%
% Usage:
%   report = check_HighOrderSequence('fast_gstf_HighOrder_*.seq');
%   report = check_HighOrderSequence(seqFile,'acq_manifest_HighOrder.mat');
%
% Checks include:
%   - number of TRs / ADCs and total duration
%   - x/y/z ordering
%   - acquisition order: triangle -> PE -> polarity -> slice
%   - physical slice positions
%   - actual 2D PE indices and ordering
%   - positive/negative probing-gradient polarity
%   - probing-gradient rise time / amplitude consistency and start time
%   - ADC sample count, dwell and seven ADC start times
%   - optional acquisition-by-acquisition manifest agreement
%
% Copyright (c) 2026 Jinyuan Zhang

if nargin < 1 || isempty(seqFile)
    files = dir('fast_gstf_HighOrder_*.seq');
    if isempty(files)
        error('No fast_gstf_HighOrder_*.seq file found in the current folder.');
    end
    [~,idx] = max([files.datenum]);
    seqFile = fullfile(files(idx).folder,files(idx).name);
end
if nargin < 2
    manifestFile = [];
end

if ~exist(seqFile,'file')
    error('Sequence file does not exist: %s',seqFile);
end

fprintf('==============================================================\n');
fprintf('Fast High-Order GSTF sequence validation\n');
fprintf('File: %s\n',seqFile);
fprintf('==============================================================\n');

seqText = fileread(seqFile);
defs = readDefinitions(seqText);

%% Read protocol definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TR = getDefinitionNumbers(defs,'TR');
blockRaster = getDefinitionNumbers(defs,'BlockDurationRaster');
totalDurationDefinition = getDefinitionNumbers(defs,'TotalDuration');
nPE = getDefinitionNumbers(defs,'HighOrder_nPE');
numPEAcq = getDefinitionNumbers(defs,'HighOrder_numPEAcq');
PEFOV = getDefinitionNumbers(defs,'HighOrder_PEFOV');
PEMode = getDefinitionString(defs,'HighOrder_PEMode');
sliceOffsets = getDefinitionNumbers(defs,'HighOrder_sliceOffsets_m');
nRepetition = getDefinitionNumbers(defs,'nRepetition');
probeRise = getDefinitionNumbers(defs,'GSTF_rise');
probeAmp_mT = getDefinitionNumbers(defs,'GSTF_amp_mT');
probeDelays = getDefinitionNumbers(defs,'GSTF_delays');
probeScheme = getDefinitionNumbers(defs,'GSTF_scheme_type');
adcSamplesDefinition = getDefinitionNumbers(defs,'adcSamples');
adcDwellDefinition = getDefinitionNumbers(defs,'adcDwell');
triangleStartDefinition = getDefinitionNumbers(defs,'TriangleStartTime');

numSlices = length(sliceOffsets);
numTriangles = length(probeRise);
numAxes = 3;
numPolarities = 2;

if numTriangles ~= 7
    error('Expected seven Fast-GIRF probing measurements, found %d.',numTriangles);
end
if length(probeAmp_mT) ~= numTriangles || length(probeScheme) ~= numTriangles || ...
        length(probeDelays) ~= numTriangles
    error('GSTF probe definitions have inconsistent lengths.');
end

expectedOrderString = 'axis-repetition-triangle-PE-polarity-slice';
if isfield(defs,'HighOrder_AcquisitionOrder')
    orderString = getDefinitionString(defs,'HighOrder_AcquisitionOrder');
    if ~strcmp(orderString,expectedOrderString)
        error('Unsupported High-Order acquisition-order definition: %s',orderString);
    end
else
    orderString = 'not defined (actual order will still be checked)';
end

expectedTRs = numAxes*nRepetition*numTriangles*numPolarities*numPEAcq*numSlices;

%% Parse Pulseq event tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blocks = readNumericSection(seqText,'BLOCKS',8);
trapEvents = readNumericSection(seqText,'TRAP',6);
adcEvents = readNumericSection(seqText,'ADC',9);
rfEvents = readNumericSection(seqText,'RF',11);

if isempty(blocks) || isempty(trapEvents) || isempty(adcEvents) || isempty(rfEvents)
    error('Could not parse one or more required Pulseq sections.');
end

numADCBlocks = sum(blocks(:,7) ~= 0);
if numADCBlocks ~= expectedTRs
    error('ADC count mismatch: expected %d, found %d.',expectedTRs,numADCBlocks);
end

actualDuration = sum(blocks(:,2))*blockRaster;
expectedDuration = expectedTRs*TR;
tolTime = max(blockRaster,1e-9);
if abs(actualDuration-expectedDuration) > tolTime
    error('Sequence duration does not equal expectedTRs x TR.');
end
if abs(actualDuration-totalDurationDefinition) > tolTime
    error('Pulseq TotalDuration definition does not match block durations.');
end

%% Split block list into individual TRs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TRRaster = round(TR/blockRaster);
TRStart = zeros(expectedTRs,1);
TREnd = zeros(expectedTRs,1);
trCounter = 0;
currentStart = 1;
currentDuration = 0;

for blockIndex=1:size(blocks,1)
    currentDuration = currentDuration + blocks(blockIndex,2);
    if currentDuration == TRRaster
        trCounter = trCounter + 1;
        if trCounter > expectedTRs
            error('More TRs were found than expected.');
        end
        TRStart(trCounter) = currentStart;
        TREnd(trCounter) = blockIndex;
        currentStart = blockIndex + 1;
        currentDuration = 0;
    elseif currentDuration > TRRaster
        error('Block durations exceed the requested TR near block %d.',blockIndex);
    end
end
if currentDuration ~= 0 || trCounter ~= expectedTRs
    error('Could not partition the complete sequence into %d TRs.',expectedTRs);
end

%% Expected PE mask/order %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~mod(nPE,2)
    error('High-Order nPE must be odd.');
end
peIndices = -(nPE-1)/2:(nPE-1)/2;
peAreasExpected = peIndices/PEFOV;
[PE1Grid,PE2Grid] = ndgrid(1:nPE,1:nPE);
if strcmpi(PEMode,'elliptical')
    radius = (nPE-1)/2;
    if radius == 0
        peMask = true(nPE,nPE);
    else
        peMask = (peIndices(PE1Grid)/radius).^2 + ...
                 (peIndices(PE2Grid)/radius).^2 <= 1 + eps;
    end
elseif strcmpi(PEMode,'full')
    peMask = true(nPE,nPE);
else
    error('Unsupported HighOrder_PEMode: %s',PEMode);
end
[validPE1,validPE2] = find(peMask);
validPhaseEncodes = [validPE1,validPE2];
if size(validPhaseEncodes,1) ~= numPEAcq
    error('HighOrder_numPEAcq does not match nPE/PEMode.');
end

%% Inspect every TR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adcStart_us = zeros(expectedTRs,1);
axisCode = zeros(expectedTRs,1);       % 1=x, 2=y, 3=z
rfID = zeros(expectedTRs,1);
peArea = zeros(expectedTRs,2);
testAmplitude_HzPerM = zeros(expectedTRs,1);
testRise_us = zeros(expectedTRs,1);
testStart_us = zeros(expectedTRs,1);
adcSamplesActual = zeros(expectedTRs,1);
adcDwell_s = zeros(expectedTRs,1);

TRsPerPE = numPolarities*numSlices;
TRsPerTriangle = numPEAcq*TRsPerPE;
TRsPerRepetition = numTriangles*TRsPerTriangle;
TRsPerAxis = nRepetition*TRsPerRepetition;

for tr=1:expectedTRs
    firstBlock = TRStart(tr);
    lastBlock = TREnd(tr);
    localBlocks = blocks(firstBlock:lastBlock,:);

    adcLocal = find(localBlocks(:,7) ~= 0);
    if numel(adcLocal) ~= 1
        error('TR %d contains %d ADC events; exactly one is required.',tr,numel(adcLocal));
    end
    adcBlock = firstBlock + adcLocal - 1;
    adcEventID = blocks(adcBlock,7);
    adcRow = find(adcEvents(:,1) == adcEventID,1);
    if isempty(adcRow)
        error('ADC event ID %d referenced by TR %d was not found.',adcEventID,tr);
    end
    adcSamplesActual(tr) = adcEvents(adcRow,2);
    adcDwell_s(tr) = adcEvents(adcRow,3)*1e-9; % Pulseq ADC dwell is stored in ns
    adcDelay_us = adcEvents(adcRow,4);
    adcStart_us(tr) = sum(blocks(firstBlock:adcBlock-1,2))*blockRaster*1e6 + adcDelay_us;

    rfLocal = find(localBlocks(:,3) ~= 0);
    if numel(rfLocal) ~= 1
        error('TR %d contains %d RF events; exactly one is required.',tr,numel(rfLocal));
    end
    rfBlock = firstBlock + rfLocal - 1;
    rfID(tr) = blocks(rfBlock,3);

    sliceGradientEvents = blocks(rfBlock,4:6);
    drivenAxis = find(sliceGradientEvents ~= 0);
    if numel(drivenAxis) ~= 1
        error('Could not determine slice-selection axis in TR %d.',tr);
    end
    axisCode(tr) = drivenAxis;

    spatialBlock = rfBlock + 1;
    if spatialBlock > lastBlock
        error('Spatial-encoding block is missing in TR %d.',tr);
    end
    spatialEvents = blocks(spatialBlock,4:6);
    if drivenAxis == 1
        peAxes = [2,3];
    elseif drivenAxis == 2
        peAxes = [1,3];
    else
        peAxes = [2,1];
    end
    peArea(tr,1) = getTrapArea(spatialEvents(peAxes(1)),trapEvents);
    peArea(tr,2) = getTrapArea(spatialEvents(peAxes(2)),trapEvents);

    % Triangle is an outer acquisition dimension for both supported loop
    % orders, so it can be determined independently of PE/polarity order.
    trWithinAxis = mod(tr-1,TRsPerAxis)+1;
    trWithinRep = mod(trWithinAxis-1,TRsPerRepetition)+1;
    triangleIndex = floor((trWithinRep-1)/TRsPerTriangle)+1;

    if probeScheme(triangleIndex) == 1
        testBlock = adcBlock;
    else
        testBlock = firstBlock;
    end
    testEventID = blocks(testBlock,3+drivenAxis);
    testRow = find(trapEvents(:,1) == testEventID,1);
    if isempty(testRow)
        error('Could not identify probing gradient in TR %d.',tr);
    end
    testAmplitude_HzPerM(tr) = trapEvents(testRow,2);
    testRise_us(tr) = trapEvents(testRow,3);
    testStart_us(tr) = sum(blocks(firstBlock:testBlock-1,2))*blockRaster*1e6 + trapEvents(testRow,6);
end

if any(adcSamplesActual ~= adcSamplesDefinition)
    error('One or more ADC events do not match the adcSamples definition.');
end
if max(abs(adcDwell_s-adcDwellDefinition)) > 1e-12
    error('One or more ADC dwell times do not match the adcDwell definition.');
end
if max(abs(testStart_us-triangleStartDefinition*1e6)) > 0.5
    error('The probing gradient does not start at TriangleStartTime in every TR.');
end

%% Recover physical slice index and PE index from actual events %%%%%%%%%%%
uniqueRF = unique(rfID,'stable');
if numel(uniqueRF) ~= numSlices
    error('Expected %d slice RF events, found %d.',numSlices,numel(uniqueRF));
end
actualSliceIndex = zeros(expectedTRs,1);
for sl=1:numSlices
    actualSliceIndex(rfID == uniqueRF(sl)) = sl;
end

firstRFBlock = TRStart(1) + find(blocks(TRStart(1):TREnd(1),3) ~= 0,1) - 1;
sliceGradientID = blocks(firstRFBlock,4);
sliceGradientRow = find(trapEvents(:,1) == sliceGradientID,1);
if isempty(sliceGradientRow)
    error('Could not identify slice-selection gradient event.');
end
sliceGradient_HzPerM = trapEvents(sliceGradientRow,2);
actualSliceOffsets = zeros(1,numSlices);
for i=1:numSlices
    row = find(rfEvents(:,1) == uniqueRF(i),1);
    if isempty(row)
        error('RF event ID %d was not found.',uniqueRF(i));
    end
    actualSliceOffsets(i) = rfEvents(row,10)/sliceGradient_HzPerM;
end
if max(abs(actualSliceOffsets-sliceOffsets)) > 1e-6
    error('RF ordering/offsets do not reproduce HighOrder_sliceOffsets_m.');
end

actualPE1Index = zeros(expectedTRs,1);
actualPE2Index = zeros(expectedTRs,1);
for tr=1:expectedTRs
    [d1,i1] = min(abs(peAreasExpected-peArea(tr,1)));
    [d2,i2] = min(abs(peAreasExpected-peArea(tr,2)));
    if d1 > 1e-4 || d2 > 1e-4
        error('TR %d contains a PE area outside the requested PE grid.',tr);
    end
    actualPE1Index(tr) = i1;
    actualPE2Index(tr) = i2;
end
actualPolarity = sign(testAmplitude_HzPerM);

%% Verify exact acquisition ordering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expectedAxisCode = zeros(expectedTRs,1);
expectedRep = zeros(expectedTRs,1);
expectedTri = zeros(expectedTRs,1);
expectedPE1 = zeros(expectedTRs,1);
expectedPE2 = zeros(expectedTRs,1);
expectedPolarity = zeros(expectedTRs,1);
expectedSlice = zeros(expectedTRs,1);

row = 0;
for ax=1:numAxes
    for rep=1:nRepetition
        for tri=1:numTriangles
            for pe=1:numPEAcq
                for pol=[1,-1]
                    for sl=1:numSlices
                        row = row + 1;
                        expectedAxisCode(row) = ax;
                        expectedRep(row) = rep;
                        expectedTri(row) = tri;
                        expectedPE1(row) = validPhaseEncodes(pe,1);
                        expectedPE2(row) = validPhaseEncodes(pe,2);
                        expectedPolarity(row) = pol;
                        expectedSlice(row) = sl;
                    end
                end
            end
        end
    end
end

checkVector(axisCode,expectedAxisCode,'gradient axis');
checkVector(actualPE1Index,expectedPE1,'PE1 order');
checkVector(actualPE2Index,expectedPE2,'PE2 order');
checkVector(actualPolarity,expectedPolarity,'polarity order');
checkVector(actualSliceIndex,expectedSlice,'slice order');

%% Check probing amplitudes/rise times and +/- symmetry %%%%%%%%%%%%%%%%%%%
ampArray = reshape(testAmplitude_HzPerM, ...
    [numSlices,numPolarities,numPEAcq,numTriangles,nRepetition,numAxes]);
riseArray = reshape(testRise_us, ...
    [numSlices,numPolarities,numPEAcq,numTriangles,nRepetition,numAxes]);

inferredGamma = zeros(1,numTriangles);
for triangle=1:numTriangles
    positiveAmp = ampArray(:,1,:,triangle,:,:);
    negativeAmp = ampArray(:,2,:,triangle,:,:);
    if any(positiveAmp(:) <= 0) || any(negativeAmp(:) >= 0)
        error('Positive/negative polarity is incorrect for triangle %d.',triangle);
    end
    if max(abs(abs(positiveAmp(:))-abs(negativeAmp(:)))) > 1e-6*max(abs(positiveAmp(:)))
        error('Positive/negative amplitudes differ for triangle %d.',triangle);
    end
    currentRise = riseArray(:,:,:,triangle,:,:);
    if max(abs(currentRise(:)-probeRise(triangle)*1e6)) > 0.5
        error('Probing-gradient rise time is incorrect for triangle %d.',triangle);
    end
    inferredGamma(triangle) = mean(abs(positiveAmp(:))) / (probeAmp_mT(triangle)*1e-3);
end
if (max(inferredGamma)-min(inferredGamma))/mean(inferredGamma) > 1e-4
    error('Probing-gradient amplitudes are inconsistent with GSTF_amp_mT.');
end

%% Check seven ADC start times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adcArray = reshape(adcStart_us, ...
    [numSlices,numPolarities,numPEAcq,numTriangles,nRepetition,numAxes]);
adcStartByTriangle_us = zeros(1,numTriangles);
for triangle=1:numTriangles
    currentADC = adcArray(:,:,:,triangle,:,:);
    if max(currentADC(:))-min(currentADC(:)) > 0.5
        error('ADC start time varies across spatial encodings for triangle %d.',triangle);
    end
    adcStartByTriangle_us(triangle) = currentADC(1);
end

if isfield(defs,'TimeShiftADC_us')
    shiftDefinition = getDefinitionNumbers(defs,'TimeShiftADC_us');
    if length(shiftDefinition) ~= numTriangles || ...
            max(abs(shiftDefinition-adcStartByTriangle_us)) > 0.5
        error('TimeShiftADC_us is inconsistent with actual ADC block timing.');
    end
    shiftDefinitionMode = 'physical time [us]';
else
    shiftDefinitionMode = 'not available';
end
if isfield(defs,'TimeShiftADC_index')
    shiftIndex = getDefinitionNumbers(defs,'TimeShiftADC_index');
    if length(shiftIndex) ~= numTriangles || ...
            max(abs((shiftIndex-1)-adcStartByTriangle_us)) > 0.5
        error('TimeShiftADC_index is inconsistent with the 1-us nominal input grid.');
    end
end

%% Optional acquisition-by-acquisition manifest check %%%%%%%%%%%%%%%%%%%%%
manifestMatched = false;
if ~isempty(manifestFile)
    if ischar(manifestFile) || isstring(manifestFile)
        manifestData = load(manifestFile);
        if ~isfield(manifestData,'acq_manifest')
            error('Manifest file does not contain acq_manifest.');
        end
        acq_manifest = manifestData.acq_manifest;
    elseif istable(manifestFile)
        acq_manifest = manifestFile;
    else
        error('manifestFile must be a path or table.');
    end

    if height(acq_manifest) ~= expectedTRs
        error('Manifest row count does not match the sequence.');
    end
    checkVector(lower(string(acq_manifest.Axis)), ...
        [repmat("x",TRsPerAxis,1);repmat("y",TRsPerAxis,1);repmat("z",TRsPerAxis,1)],'manifest Axis');
    checkVector(double(acq_manifest.Repetition),expectedRep,'manifest Repetition');
    checkVector(double(acq_manifest.Triangle),expectedTri,'manifest Triangle');
    checkVector(double(acq_manifest.PE1Index),actualPE1Index,'manifest PE1Index');
    checkVector(double(acq_manifest.PE2Index),actualPE2Index,'manifest PE2Index');
    checkVector(double(acq_manifest.Polarity),actualPolarity,'manifest Polarity');
    checkVector(double(acq_manifest.SliceIndex),actualSliceIndex,'manifest SliceIndex');

    if ismember('Scheme',acq_manifest.Properties.VariableNames)
        expectedScheme = probeScheme(expectedTri).';
        checkVector(double(acq_manifest.Scheme),expectedScheme,'manifest Scheme');
    end
    if ismember('Delay_s',acq_manifest.Properties.VariableNames)
        expectedDelay = probeDelays(expectedTri).';
        if max(abs(double(acq_manifest.Delay_s(:))-expectedDelay(:))) > 1e-9
            error('Manifest Delay_s does not match the sequence.');
        end
    end
    if ismember('SlicePosition_m',acq_manifest.Properties.VariableNames)
        expectedPosition = sliceOffsets(actualSliceIndex).';
        if max(abs(double(acq_manifest.SlicePosition_m(:))-expectedPosition(:))) > 1e-9
            error('Manifest SlicePosition_m does not match the sequence.');
        end
    end
    manifestMatched = true;
end

%% Report %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('TRs / ADCs             : %d / %d\n',expectedTRs,numADCBlocks);
fprintf('Total duration         : %.3f s (%.2f min)\n',actualDuration,actualDuration/60);
fprintf('Acquisition order      : %s\n',orderString);
fprintf('Spatial encoding       : %d x %d %s, %d acquired PE points\n', ...
    nPE,nPE,PEMode,numPEAcq);
fprintf('Slice positions [mm]   :'); fprintf(' %.1f',actualSliceOffsets*1e3); fprintf('\n');
fprintf('TRs per gradient axis  : %d\n',TRsPerAxis);
fprintf('ADC start times [us]   :'); fprintf(' %.0f',adcStartByTriangle_us); fprintf('\n');
fprintf('Inferred gamma [MHz/T] : %.6f\n',mean(inferredGamma)/1e6);
fprintf('TimeShiftADC definition: %s\n',shiftDefinitionMode);
if manifestMatched
    fprintf('Manifest pairing       : acquisition-by-acquisition match\n');
end
fprintf('--------------------------------------------------------------\n');
fprintf('High-Order sequence validation PASSED.\n');
fprintf('==============================================================\n');

report = struct();
report.seqFile = seqFile;
report.passed = true;
report.numTRs = expectedTRs;
report.numADCs = numADCBlocks;
report.totalDuration_s = actualDuration;
report.scanTime_min = actualDuration/60;
report.nPE = nPE;
report.PEMode = PEMode;
report.numPEAcq = numPEAcq;
report.sliceOffsets_m = actualSliceOffsets;
report.adcStartByTriangle_us = adcStartByTriangle_us;
report.inferredGamma_HzPerT = mean(inferredGamma);
report.TimeShiftADCDefinitionMode = shiftDefinitionMode;
report.manifestMatched = manifestMatched;
report.actualPE1Index = actualPE1Index;
report.actualPE2Index = actualPE2Index;
report.actualPolarity = actualPolarity;
report.actualSliceIndex = actualSliceIndex;

end

%% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function defs = readDefinitions(seqText)
section = getSectionText(seqText,'DEFINITIONS');
lines = regexp(section,'\r?\n','split');
defs = struct();
for i=1:length(lines)
    line = strtrim(lines{i});
    if isempty(line) || startsWith(line,'#'), continue; end
    parts = regexp(line,'\s+','split');
    if length(parts) >= 2
        defs.(parts{1}) = strjoin(parts(2:end),' ');
    end
end
end

function values = getDefinitionNumbers(defs,name)
if ~isfield(defs,name), error('Required Pulseq definition is missing: %s',name); end
values = sscanf(defs.(name),'%f').';
end

function value = getDefinitionString(defs,name)
if ~isfield(defs,name), error('Required Pulseq definition is missing: %s',name); end
value = strtrim(defs.(name));
end

function rows = readNumericSection(seqText,sectionName,numColumns)
section = getSectionText(seqText,sectionName);
lines = regexp(section,'\r?\n','split');
rows = zeros(0,numColumns);
for i=1:length(lines)
    line = strtrim(lines{i});
    if isempty(line) || startsWith(line,'#'), continue; end
    values = sscanf(line,'%f').';
    if length(values) >= numColumns
        rows(end+1,:) = values(1:numColumns); %#ok<AGROW>
    end
end
end

function section = getSectionText(seqText,sectionName)
header = ['[',sectionName,']'];
startIndex = strfind(seqText,header);
if isempty(startIndex), error('Pulseq section not found: %s',header); end
startIndex = startIndex(1)+length(header);
remaining = seqText(startIndex:end);
nextHeader = regexp(remaining,'\r?\n\[','once');
if isempty(nextHeader), section = remaining; else, section = remaining(1:nextHeader-1); end
end

function area = getTrapArea(eventID,trapEvents)
if eventID == 0, area = 0; return; end
row = find(trapEvents(:,1) == eventID,1);
if isempty(row), error('Gradient event ID %d is not a trapezoid event.',eventID); end
amplitude = trapEvents(row,2);
rise = trapEvents(row,3)*1e-6;
flat = trapEvents(row,4)*1e-6;
fall = trapEvents(row,5)*1e-6;
area = amplitude*(0.5*rise+flat+0.5*fall);
end

function checkVector(actual,expected,name)
actual = actual(:); expected = expected(:);
if numel(actual) ~= numel(expected), error('%s length mismatch.',name); end
idx = find(actual ~= expected,1);
if ~isempty(idx), error('%s mismatch at acquisition %d.',name,idx); end
end
