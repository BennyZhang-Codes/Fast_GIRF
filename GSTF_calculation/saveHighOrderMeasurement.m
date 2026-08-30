function saveHighOrderMeasurement(outputFile, kspace, dwelltime, FOV, axisName, PosSlices)
% Save spatial High-Order GIRF data in the format expected by Fast_GIRF.
%
% Arguments:
%   outputFile: path of the measurement .mat file
%   kspace:     sorted High-Order data from sortHighOrderAcquisitions()
%   dwelltime:  ADC dwell time [s]
%   FOV:        phase-encoding FOV [mm]
%   axisName:   test-gradient/slice-selection axis: 'x', 'y', or 'z'
%   PosSlices:  slice positions relative to isocenter [mm]
%
% The saved variables match calculateOutputGradient.m:
% dwelltime, FOV, orientation, numSlices, PosSlices, and kspace.

if ~ismember(axisName,{'x','y','z'})
    error('axisName must be ''x'', ''y'', or ''z''.');
end

if strcmp(axisName,'x')
    orientation = 'dSag';
elseif strcmp(axisName,'y')
    orientation = 'dCor';
else
    orientation = 'dTra';
end

PosSlices = PosSlices(:).';
numSlices = length(PosSlices);

if size(kspace,5) ~= numSlices
    error('Number of k-space slices does not match PosSlices.');
end
if size(kspace,3) ~= size(kspace,4)
    error('Fast_GIRF requires the two PE dimensions to have equal size.');
end

save(outputFile,'dwelltime','FOV','orientation','numSlices','PosSlices','kspace','-v7.3');
fprintf('Saved High-Order GIRF measurement: %s\n',outputFile);

end
