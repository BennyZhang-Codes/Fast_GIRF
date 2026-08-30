function validateHighOrderManifest(acq_manifest, Protocol)
% Validate a High-Order acquisition manifest against the sequence protocol.
%
% The manifest is not only checked for its total number of rows. Every row
% is checked against the expected acquisition order:
%
%   axis -> repetition -> triangle -> PE -> polarity -> slice
%
% This prevents a manifest from another protocol with the same total ADC
% count from silently being used to sort the Twix data.

if ischar(acq_manifest) || isstring(acq_manifest)
    data = load(acq_manifest);
    if ~isfield(data,'acq_manifest')
        error('Manifest file does not contain acq_manifest.');
    end
    acq_manifest = data.acq_manifest;
end

if ~istable(acq_manifest)
    error('acq_manifest must be a table or a path to a manifest .mat file.');
end

requiredVariables = {'Axis','Repetition','Triangle','Polarity', ...
    'PE1Index','PE2Index','SliceIndex'};
for i=1:length(requiredVariables)
    if ~ismember(requiredVariables{i},acq_manifest.Properties.VariableNames)
        error('Manifest is missing required variable: %s',requiredVariables{i});
    end
end

if height(acq_manifest) ~= Protocol.expectedAcquisitions
    error('Manifest contains %d rows; sequence protocol expects %d acquisitions.', ...
        height(acq_manifest),Protocol.expectedAcquisitions);
end

%% Expected PE pattern %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~mod(Protocol.nPE,2)
    error('High-Order PE matrix size must be odd.');
end

peIndices = -(Protocol.nPE-1)/2:(Protocol.nPE-1)/2;
[PE1Grid,PE2Grid] = ndgrid(1:Protocol.nPE,1:Protocol.nPE);

if strcmpi(Protocol.PEMode,'elliptical')
    radius = (Protocol.nPE-1)/2;
    if radius == 0
        peMask = true(Protocol.nPE,Protocol.nPE);
    else
        peMask = (peIndices(PE1Grid)/radius).^2 + ...
                 (peIndices(PE2Grid)/radius).^2 <= 1 + eps;
    end
elseif strcmpi(Protocol.PEMode,'full')
    peMask = true(Protocol.nPE,Protocol.nPE);
else
    error('Unsupported High-Order PE mode: %s',Protocol.PEMode);
end

[validPE1,validPE2] = find(peMask);
validPhaseEncodes = [validPE1,validPE2];

if size(validPhaseEncodes,1) ~= Protocol.numPEAcq
    error('Sequence numPEAcq does not match the PE mask implied by nPE/PEMode.');
end

%% Build expected acquisition columns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numRows = Protocol.expectedAcquisitions;
expectedAxis = strings(numRows,1);
expectedRep = zeros(numRows,1);
expectedTri = zeros(numRows,1);
expectedPol = zeros(numRows,1);
expectedPE1 = zeros(numRows,1);
expectedPE2 = zeros(numRows,1);
expectedSlice = zeros(numRows,1);
expectedScheme = zeros(numRows,1);
expectedDelay = zeros(numRows,1);
expectedSlicePosition = zeros(numRows,1);

axes_to_test = {'x','y','z'};
row = 0;
for ax_idx=1:length(axes_to_test)
    for rep=1:Protocol.nRepetition
        for tri=1:Protocol.numTriangles
            for pe=1:Protocol.numPEAcq
                for pol=[1,-1]
                    for sl=1:Protocol.numSlices
                        row = row + 1;
                        expectedAxis(row) = axes_to_test{ax_idx};
                        expectedRep(row) = rep;
                        expectedTri(row) = tri;
                        expectedPol(row) = pol;
                        expectedPE1(row) = validPhaseEncodes(pe,1);
                        expectedPE2(row) = validPhaseEncodes(pe,2);
                        expectedSlice(row) = sl;
                        expectedScheme(row) = Protocol.probeScheme(tri);
                        expectedDelay(row) = Protocol.probeDelays_s(tri);
                        expectedSlicePosition(row) = Protocol.sliceOffsets_m(sl);
                    end
                end
            end
        end
    end
end

if row ~= numRows
    error('Internal manifest validation error: expected row count was not reached.');
end

%% Compare required columns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
manifestAxis = string(acq_manifest.Axis);
manifestAxis = lower(strtrim(manifestAxis));

checkExact(manifestAxis,expectedAxis,'Axis');
checkExact(double(acq_manifest.Repetition),expectedRep,'Repetition');
checkExact(double(acq_manifest.Triangle),expectedTri,'Triangle');
checkExact(double(acq_manifest.Polarity),expectedPol,'Polarity');
checkExact(double(acq_manifest.PE1Index),expectedPE1,'PE1Index');
checkExact(double(acq_manifest.PE2Index),expectedPE2,'PE2Index');
checkExact(double(acq_manifest.SliceIndex),expectedSlice,'SliceIndex');

%% Compare optional metadata when present %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember('Scheme',acq_manifest.Properties.VariableNames)
    checkExact(double(acq_manifest.Scheme),expectedScheme,'Scheme');
end
if ismember('Delay_s',acq_manifest.Properties.VariableNames)
    checkNumeric(double(acq_manifest.Delay_s),expectedDelay,1e-9,'Delay_s');
end
if ismember('SlicePosition_m',acq_manifest.Properties.VariableNames)
    checkNumeric(double(acq_manifest.SlicePosition_m),expectedSlicePosition,1e-9,'SlicePosition_m');
end
if ismember('PE1Value',acq_manifest.Properties.VariableNames)
    expectedValue = peIndices(expectedPE1).';
    checkExact(double(acq_manifest.PE1Value),expectedValue,'PE1Value');
end
if ismember('PE2Value',acq_manifest.Properties.VariableNames)
    expectedValue = peIndices(expectedPE2).';
    checkExact(double(acq_manifest.PE2Value),expectedValue,'PE2Value');
end

fprintf('High-Order acquisition manifest validation passed: %d rows.\n',numRows);

end

%% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkExact(actual,expected,name)
actual = actual(:);
expected = expected(:);
if numel(actual) ~= numel(expected)
    error('Manifest %s has the wrong number of entries.',name);
end
idx = find(actual ~= expected,1);
if ~isempty(idx)
    error('Manifest %s mismatch at acquisition %d.',name,idx);
end
end

function checkNumeric(actual,expected,tolerance,name)
actual = actual(:);
expected = expected(:);
if numel(actual) ~= numel(expected)
    error('Manifest %s has the wrong number of entries.',name);
end
idx = find(~isfinite(actual) | abs(actual-expected)>tolerance,1);
if ~isempty(idx)
    error('Manifest %s mismatch at acquisition %d.',name,idx);
end
end
