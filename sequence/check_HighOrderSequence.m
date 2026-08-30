function report = check_HighOrderSequence(seqFile)
% Validate a generated Fast High-Order GSTF Pulseq sequence.
%
% This checker parses the written .seq file directly. It therefore verifies
% the final sequence file that will be transferred to the scanner rather
% than relying on variables left in the sequence-generation workspace.
%
% Usage:
%   report = check_HighOrderSequence('fast_gstf_HighOrder_7T_pe5_elliptical_tr500_fa90_rep1.seq');
%
% The following are checked:
%   - number of TRs and ADC events
%   - total sequence duration
%   - x/y/z acquisition ordering
%   - four-slice RF frequency offsets / physical slice positions
%   - acquired 2D PE pattern and PE areas
%   - positive/negative probing-gradient polarity
%   - probing-gradient rise times and relative amplitudes
%   - ADC start time for all seven Fast-GIRF probing measurements
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

if ~exist(seqFile,'file')
    error('Sequence file does not exist: %s',seqFile);
end

fprintf('==============================================================\n');
fprintf('Fast High-Order GSTF sequence validation\n');
fprintf('File: %s\n',seqFile);
fprintf('==============================================================\n');

seqText = fileread(seqFile);
defs = readDefinitions(seqText);

%% Read High-Order protocol definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
probeScheme = getDefinitionNumbers(defs,'GSTF_scheme_type');

numSlices = length(sliceOffsets);
numTriangles = length(probeRise);
numAxes = 3;
numPolarities = 2;

if numTriangles ~= 7
    error('Expected seven Fast-GIRF probing measurements, found %d.',numTriangles);
end
if length(probeAmp_mT) ~= numTriangles || length(probeScheme) ~= numTriangles
    error('GSTF probe definitions have inconsistent lengths.');
end

expectedTRs = numAxes * nRepetition * numTriangles * ...
    numPolarities * numPEAcq * numSlices;

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

actualDuration = sum(blocks(:,2)) * blockRaster;
expectedDuration = expectedTRs * TR;
tolTime = max(blockRaster,1e-9);

if abs(actualDuration-expectedDuration) > tolTime
    error('Sequence duration does not equal expectedTRs x TR.');
end
if abs(actualDuration-totalDurationDefinition) > tolTime
    error('Pulseq TotalDuration definition does not match the block durations.');
end

%% Split the block list into individual TRs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TRRaster = round(TR/blockRaster);
TRStart = zeros(expectedTRs,1);
TREnd = zeros(expectedTRs,1);
trCounter = 0;
currentStart = 1;
currentDuration = 0;

for blockIndex = 1:size(blocks,1)
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

%% Inspect every TR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adcStart_us = zeros(expectedTRs,1);
axisCode = zeros(expectedTRs,1);       % 1=x, 2=y, 3=z
rfID = zeros(expectedTRs,1);
peArea = zeros(expectedTRs,2);
testAmplitude_HzPerM = zeros(expectedTRs,1);
testRise_us = zeros(expectedTRs,1);

TRsPerPolarity = numPEAcq*numSlices;
TRsPerTriangle = numPolarities*TRsPerPolarity;
TRsPerRepetition = numTriangles*TRsPerTriangle;
TRsPerAxis = nRepetition*TRsPerRepetition;

for tr = 1:expectedTRs
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
        error('Could not determine the slice-selection axis in TR %d.',tr);
    end
    axisCode(tr) = drivenAxis;

    spatialBlock = rfBlock + 1;
    if spatialBlock > lastBlock
        error('Spatial-encoding block is missing in TR %d.',tr);
    end
    spatialEvents = blocks(spatialBlock,4:6);

    if drivenAxis == 1
        peAxes = [2,3];       % x slice: PE1=y, PE2=z
    elseif drivenAxis == 2
        peAxes = [1,3];       % y slice: PE1=x, PE2=z
    else
        peAxes = [2,1];       % z slice: PE1=y, PE2=x
    end

    peArea(tr,1) = getTrapArea(spatialEvents(peAxes(1)),trapEvents);
    peArea(tr,2) = getTrapArea(spatialEvents(peAxes(2)),trapEvents);

    trWithinAxis = mod(tr-1,TRsPerAxis) + 1;
    trWithinRep = mod(trWithinAxis-1,TRsPerRepetition) + 1;
    triangleIndex = floor((trWithinRep-1)/TRsPerTriangle) + 1;

    if probeScheme(triangleIndex) == 1
        testBlock = adcBlock;
    else
        testBlock = firstBlock;
    end

    testEventID = blocks(testBlock,3+drivenAxis);
    testRow = find(trapEvents(:,1) == testEventID,1);
    if isempty(testRow)
        error('Could not identify the probing gradient in TR %d.',tr);
    end
    testAmplitude_HzPerM(tr) = trapEvents(testRow,2);
    testRise_us(tr) = trapEvents(testRow,3);
end

%% Check x/y/z ordering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expectedAxisCode = [ones(TRsPerAxis,1); 2*ones(TRsPerAxis,1); 3*ones(TRsPerAxis,1)];
if ~isequal(axisCode,expectedAxisCode)
    error('The sequence does not follow the expected x -> y -> z axis ordering.');
end

%% Check the four physical slice offsets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uniqueRF = unique(rfID,'stable');
if numel(uniqueRF) ~= numSlices
    error('Expected %d slice RF events, found %d.',numSlices,numel(uniqueRF));
end

% The RF-block slice gradient is the same physical gradient event for all
% axes. Its Hz/m amplitude together with the RF frequency offset gives the
% physical slice position directly.
firstRFBlock = TRStart(1) + find(blocks(TRStart(1):TREnd(1),3) ~= 0,1) - 1;
sliceGradientID = blocks(firstRFBlock,4);
sliceGradientRow = find(trapEvents(:,1) == sliceGradientID,1);
if isempty(sliceGradientRow)
    error('Could not identify the slice-selection gradient event.');
end
sliceGradient_HzPerM = trapEvents(sliceGradientRow,2);

actualSliceOffsets = zeros(1,numSlices);
for i = 1:numSlices
    row = find(rfEvents(:,1) == uniqueRF(i),1);
    if isempty(row)
        error('RF event ID %d was not found in the RF table.',uniqueRF(i));
    end
    rfFrequency_Hz = rfEvents(row,10);
    actualSliceOffsets(i) = rfFrequency_Hz/sliceGradient_HzPerM;
end

if max(abs(sort(actualSliceOffsets)-sort(sliceOffsets))) > 1e-6
    error('RF frequency offsets do not reproduce the requested physical slice positions.');
end

%% Check the acquired 2D PE pattern %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peIndices = -(nPE-1)/2:(nPE-1)/2;
peAreasExpected = peIndices/PEFOV;
[PE1Grid,PE2Grid] = ndgrid(1:nPE,1:nPE);

if strcmpi(PEMode,'elliptical')
    radius = (nPE-1)/2;
    peMask = (peIndices(PE1Grid)/radius).^2 + ...
             (peIndices(PE2Grid)/radius).^2 <= 1 + eps;
elseif strcmpi(PEMode,'full')
    peMask = true(nPE,nPE);
else
    error('Unsupported HighOrder_PEMode in the sequence: %s',PEMode);
end

expectedPE = [peAreasExpected(PE1Grid(peMask)).', peAreasExpected(PE2Grid(peMask)).'];
expectedPE = sortrows(expectedPE,[1 2]);

actualPEByAxis = cell(1,3);
for ax = 1:3
    currentPE = peArea(axisCode == ax,:);
    currentPE = round(currentPE*1e6)/1e6;
    currentPE = unique(currentPE,'rows');
    currentPE = sortrows(currentPE,[1 2]);
    actualPEByAxis{ax} = currentPE;

    if size(currentPE,1) ~= numPEAcq || size(expectedPE,1) ~= numPEAcq
        error('PE point count mismatch for axis %d.',ax);
    end
    if max(abs(currentPE(:)-expectedPE(:))) > 1e-4
        error('The actual 2D PE areas do not match the requested %s mask for axis %d.',PEMode,ax);
    end
end

%% Check probing gradient polarity, rise time and relative amplitudes %%%%%
ampArray = reshape(testAmplitude_HzPerM, ...
    [numSlices,numPEAcq,numPolarities,numTriangles,nRepetition,numAxes]);
riseArray = reshape(testRise_us, ...
    [numSlices,numPEAcq,numPolarities,numTriangles,nRepetition,numAxes]);

inferredGamma = zeros(1,numTriangles);
for triangle = 1:numTriangles
    positiveAmp = ampArray(:,:,1,triangle,:,:);
    negativeAmp = ampArray(:,:,2,triangle,:,:);

    if any(positiveAmp(:) <= 0) || any(negativeAmp(:) >= 0)
        error('Positive/negative polarity ordering is incorrect for triangle %d.',triangle);
    end

    if max(abs(abs(positiveAmp(:))-abs(negativeAmp(:)))) > 1e-6*max(abs(positiveAmp(:)))
        error('Positive and negative probing-gradient amplitudes differ for triangle %d.',triangle);
    end

    currentRise = riseArray(:,:,:,triangle,:,:);
    if max(abs(currentRise(:)-probeRise(triangle)*1e6)) > 0.5
        error('Probing-gradient rise time is incorrect for triangle %d.',triangle);
    end

    inferredGamma(triangle) = mean(abs(positiveAmp(:))) / (probeAmp_mT(triangle)*1e-3);
end

if (max(inferredGamma)-min(inferredGamma))/mean(inferredGamma) > 1e-4
    error('The probing-gradient amplitudes are not mutually consistent with GSTF_amp_mT.');
end

%% Check the seven ADC start times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
adcArray = reshape(adcStart_us, ...
    [numSlices,numPEAcq,numPolarities,numTriangles,nRepetition,numAxes]);
adcStartByTriangle_us = zeros(1,numTriangles);

for triangle = 1:numTriangles
    currentADC = adcArray(:,:,:,triangle,:,:);
    if max(currentADC(:))-min(currentADC(:)) > 0.5
        error('ADC start time is not constant across spatial encodings for triangle %d.',triangle);
    end
    adcStartByTriangle_us(triangle) = currentADC(1);
end

shiftDefinition = [];
if isfield(defs,'TimeShiftADC_us')
    shiftDefinition = getDefinitionNumbers(defs,'TimeShiftADC_us');
end

shiftDefinitionMode = 'not available';
if ~isempty(shiftDefinition)
    if length(shiftDefinition) ~= numTriangles
        error('TimeShiftADC_us definition has the wrong number of entries.');
    end

    deltaShift = shiftDefinition-adcStartByTriangle_us;
    if max(abs(deltaShift)) < 0.5
        shiftDefinitionMode = 'physical time [us]';
    elseif max(abs(deltaShift-1)) < 0.5
        % Legacy High-Order generator stored the 1-based 1-us sample index
        % under the TimeShiftADC_us name. The physical timing is index-1.
        shiftDefinitionMode = 'legacy 1-us sample index (physical time = value - 1 us)';
        warning(['TimeShiftADC_us contains MATLAB 1-based sample indices rather than ' ...
            'physical microseconds. The sequence timing itself is correct.']);
    else
        error('TimeShiftADC_us is inconsistent with the actual ADC block timing.');
    end
end

%% Report %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axisNames = {'x','y','z'};

fprintf('TRs / ADCs             : %d / %d\n',expectedTRs,numADCBlocks);
fprintf('Total duration         : %.3f s (%.2f min)\n',actualDuration,actualDuration/60);
fprintf('Spatial encoding       : %d x %d %s, %d acquired PE points\n', ...
    nPE,nPE,PEMode,numPEAcq);
fprintf('Slice positions [mm]   :');
fprintf(' %.1f',actualSliceOffsets*1e3);
fprintf('\n');
fprintf('TRs per gradient axis  : %d\n',TRsPerAxis);
fprintf('ADC start times [us]   :');
fprintf(' %.0f',adcStartByTriangle_us);
fprintf('\n');
fprintf('Inferred gamma [MHz/T] : %.6f\n',mean(inferredGamma)/1e6);
fprintf('TimeShiftADC definition: %s\n',shiftDefinitionMode);

for ax = 1:3
    fprintf('%s-axis PE points       : %d\n',axisNames{ax},size(actualPEByAxis{ax},1));
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
report.actualPEByAxis = actualPEByAxis;

end

%% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function defs = readDefinitions(seqText)
section = getSectionText(seqText,'DEFINITIONS');
lines = regexp(section,'\r?\n','split');
defs = struct();

for i = 1:length(lines)
    line = strtrim(lines{i});
    if isempty(line) || startsWith(line,'#')
        continue;
    end
    parts = regexp(line,'\s+','split');
    if length(parts) < 2
        continue;
    end
    defs.(parts{1}) = strjoin(parts(2:end),' ');
end
end

function values = getDefinitionNumbers(defs,name)
if ~isfield(defs,name)
    error('Required Pulseq definition is missing: %s',name);
end
values = sscanf(defs.(name),'%f').';
end

function value = getDefinitionString(defs,name)
if ~isfield(defs,name)
    error('Required Pulseq definition is missing: %s',name);
end
value = strtrim(defs.(name));
end

function rows = readNumericSection(seqText,sectionName,numColumns)
section = getSectionText(seqText,sectionName);
lines = regexp(section,'\r?\n','split');
rows = zeros(0,numColumns);

for i = 1:length(lines)
    line = strtrim(lines{i});
    if isempty(line) || startsWith(line,'#')
        continue;
    end
    values = sscanf(line,'%f').';
    if length(values) >= numColumns
        rows(end+1,:) = values(1:numColumns); %#ok<AGROW>
    end
end
end

function section = getSectionText(seqText,sectionName)
header = ['[',sectionName,']'];
startIndex = strfind(seqText,header);
if isempty(startIndex)
    error('Pulseq section not found: %s',header);
end
startIndex = startIndex(1) + length(header);
remaining = seqText(startIndex:end);
nextHeader = regexp(remaining,'\r?\n\[','once');

if isempty(nextHeader)
    section = remaining;
else
    section = remaining(1:nextHeader-1);
end
end

function area = getTrapArea(eventID,trapEvents)
if eventID == 0
    area = 0;
    return;
end

row = find(trapEvents(:,1) == eventID,1);
if isempty(row)
    error('Gradient event ID %d is not a trapezoid event.',eventID);
end

amplitude = trapEvents(row,2);  % Hz/m
rise = trapEvents(row,3)*1e-6;
flat = trapEvents(row,4)*1e-6;
fall = trapEvents(row,5)*1e-6;
area = amplitude*(0.5*rise + flat + 0.5*fall); % 1/m
end
