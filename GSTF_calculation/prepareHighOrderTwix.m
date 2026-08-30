function outputFiles = prepareHighOrderTwix(datFile, manifestFile, outputPath, dwelltime, FOV, PosSlices, nPE, firstAcquisition)
% Prepare x/y/z High-Order Fast_GIRF measurement files from Siemens Twix.
%
% Arguments:
%   datFile:          Siemens raw-data .dat file
%   manifestFile:     path to acq_manifest_HighOrder.mat
%   outputPath:       folder for measurement_H_fast_HighOrder_*.mat
%   dwelltime:        ADC dwell time [s]
%   FOV:              phase-encoding FOV [mm]
%   PosSlices:        slice positions relative to isocenter [mm]
%   nPE:              full PE matrix size; use [] to infer from manifest
%   firstAcquisition: optional first image ADC if Twix contains extra image
%                     acquisitions. Leave empty when counts match exactly.
%
% Returns:
%   outputFiles: structure containing x/y/z measurement paths
%
% This is a convenience wrapper around readHighOrderTwix() and
% prepareHighOrderMeasurements(). Keeping the two underlying steps separate
% makes raw-data extraction and manifest sorting independently testable.

if nargin < 8
    firstAcquisition = [];
end
if nargin < 7
    nPE = [];
end
if nargin < 6
    error(['prepareHighOrderTwix requires datFile, manifestFile, outputPath, ' ...
        'dwelltime, FOV, and PosSlices.']);
end

%% Read sequential image ADCs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rawFID,~] = readHighOrderTwix(datFile,manifestFile,firstAcquisition);

%% Sort and save Fast_GIRF measurement files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputFiles = prepareHighOrderMeasurements(rawFID,manifestFile,outputPath, ...
    dwelltime,FOV,PosSlices,nPE);

end
