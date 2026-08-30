function Protocol = readHighOrderSequenceDefinitions(seqFile)
% Read High-Order Fast_GIRF protocol parameters from a written Pulseq file.
%
% Arguments:
%   seqFile: generated fast_gstf_HighOrder_*.seq file
%
% Returns:
%   Protocol: structure containing the acquisition parameters required for
%             High-Order raw-data preparation.
%
% Reading these values from the actual .seq file avoids manually repeating
% dwell time, PE FOV, slice positions, PE matrix size, and repetition count
% during raw-data processing.

if ~exist(seqFile,'file')
    error('Sequence file does not exist: %s',seqFile);
end

seqText = fileread(seqFile);
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

requiredDefinitions = {'TR','adcDwell','adcSamples', ...
    'HighOrder_nPE','HighOrder_PEFOV','HighOrder_PEMode', ...
    'HighOrder_numPEAcq','HighOrder_sliceOffsets_m','nRepetition', ...
    'GSTF_rise','GSTF_amp_mT','GSTF_delays','GSTF_scheme_type'};

for i = 1:length(requiredDefinitions)
    if ~isfield(defs,requiredDefinitions{i})
        error('Required High-Order sequence definition is missing: %s', ...
            requiredDefinitions{i});
    end
end

Protocol = struct();
Protocol.seqFile = seqFile;
Protocol.TR = getNumbers(defs.TR);
Protocol.adcDwell = getNumbers(defs.adcDwell);
Protocol.adcSamples = round(getNumbers(defs.adcSamples));
Protocol.nPE = round(getNumbers(defs.HighOrder_nPE));
Protocol.PEFOV_m = getNumbers(defs.HighOrder_PEFOV);
Protocol.PEFOV_mm = Protocol.PEFOV_m*1e3;
Protocol.PEMode = strtrim(defs.HighOrder_PEMode);
Protocol.numPEAcq = round(getNumbers(defs.HighOrder_numPEAcq));
Protocol.sliceOffsets_m = getNumbers(defs.HighOrder_sliceOffsets_m);
Protocol.sliceOffsets_mm = Protocol.sliceOffsets_m*1e3;
Protocol.nRepetition = round(getNumbers(defs.nRepetition));

Protocol.probeRise_s = getNumbers(defs.GSTF_rise);
Protocol.probeAmp_mT = getNumbers(defs.GSTF_amp_mT);
Protocol.probeDelays_s = getNumbers(defs.GSTF_delays);
Protocol.probeScheme = round(getNumbers(defs.GSTF_scheme_type));

Protocol.numSlices = length(Protocol.sliceOffsets_m);
Protocol.numTriangles = length(Protocol.probeRise_s);
Protocol.numPolarities = 2;
Protocol.numAxes = 3;

if Protocol.numTriangles < 1 || ...
        length(Protocol.probeAmp_mT) ~= Protocol.numTriangles || ...
        length(Protocol.probeDelays_s) ~= Protocol.numTriangles || ...
        length(Protocol.probeScheme) ~= Protocol.numTriangles
    error('GSTF probe definitions in the sequence have inconsistent lengths.');
end

Protocol.expectedAcquisitions = Protocol.numAxes * ...
    Protocol.nRepetition * Protocol.numTriangles * ...
    Protocol.numPolarities * Protocol.numPEAcq * Protocol.numSlices;

if isfield(defs,'TotalDuration')
    Protocol.totalDuration_s = getNumbers(defs.TotalDuration);
else
    Protocol.totalDuration_s = Protocol.expectedAcquisitions*Protocol.TR;
end

if isfield(defs,'TimeShiftADC_us')
    Protocol.TimeShiftADC_us = getNumbers(defs.TimeShiftADC_us);
else
    Protocol.TimeShiftADC_us = [];
end
if isfield(defs,'TimeShiftADC_index')
    Protocol.TimeShiftADC_index = round(getNumbers(defs.TimeShiftADC_index));
else
    Protocol.TimeShiftADC_index = [];
end

if isfield(defs,'HighOrder_AcquisitionOrder')
    Protocol.AcquisitionOrder = strtrim(defs.HighOrder_AcquisitionOrder);
else
    Protocol.AcquisitionOrder = '';
end
if isfield(defs,'ScannerType')
    Protocol.ScannerType = strtrim(defs.ScannerType);
else
    Protocol.ScannerType = '';
end

fprintf('Read High-Order protocol from sequence:\n');
fprintf('    nPE / mode       : %d / %s\n',Protocol.nPE,Protocol.PEMode);
fprintf('    acquired PE      : %d\n',Protocol.numPEAcq);
fprintf('    PE FOV           : %.1f mm\n',Protocol.PEFOV_mm);
fprintf('    slices           :');
fprintf(' %.1f',Protocol.sliceOffsets_mm);
fprintf(' mm\n');
fprintf('    ADC              : %d samples, %.3f us dwell\n', ...
    Protocol.adcSamples,Protocol.adcDwell*1e6);
fprintf('    probes           : %d\n',Protocol.numTriangles);
fprintf('    expected ADCs    : %d\n',Protocol.expectedAcquisitions);

end

%% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function values = getNumbers(valueString)
values = sscanf(valueString,'%f').';
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
