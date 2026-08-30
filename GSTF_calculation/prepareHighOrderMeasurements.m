function outputFiles = prepareHighOrderMeasurements(rawFID, manifestFile, outputPath, dwelltime, FOV, PosSlices, nPE)
% Prepare x/y/z High-Order Fast_GIRF measurement files from sequential ADCs.
%
% Arguments:
%   rawFID:       sequential image ADC data [numROP, coils, acquisitions]
%   manifestFile: acquisition manifest table or path to acq_manifest_HighOrder.mat
%   outputPath:   folder in which measurement_H_fast_HighOrder_*.mat are saved
%   dwelltime:    ADC dwell time [s]
%   FOV:          phase-encoding FOV [mm]
%   PosSlices:    slice positions relative to isocenter [mm]
%   nPE:          full PE matrix size. If empty, infer it from the manifest.
%
% Returns:
%   outputFiles: structure containing the saved x/y/z measurement paths
%
% This helper intentionally starts from sequential image ADCs. Scanner-
% specific Siemens Twix extraction is kept separate from the Fast_GIRF
% spatial sorting step.

if nargin < 7
    nPE = [];
end
if nargin < 6
    error(['prepareHighOrderMeasurements requires rawFID, manifestFile, ' ...
        'outputPath, dwelltime, FOV, and PosSlices.']);
end

if isempty(outputPath)
    outputPath = '.';
end
if ~exist(outputPath,'dir')
    mkdir(outputPath);
end

%% Load acquisition manifest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if istable(manifestFile)
    acq_manifest = manifestFile;
elseif ischar(manifestFile) || isstring(manifestFile)
    manifestData = load(manifestFile);
    if ~isfield(manifestData,'acq_manifest')
        error('Manifest file does not contain acq_manifest.');
    end
    acq_manifest = manifestData.acq_manifest;
else
    error('manifestFile must be a table or path to a manifest .mat file.');
end

if isempty(nPE)
    nPE = max([acq_manifest.PE1Index; acq_manifest.PE2Index]);
end

PosSlices = PosSlices(:).';
numSlices = length(PosSlices);
axes_to_prepare = {'x','y','z'};
outputFiles = struct();

%% Sort and save each physical gradient axis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ax_idx = 1:length(axes_to_prepare)
    ax = axes_to_prepare{ax_idx};

    disp(['Prepare High-Order measurement for ',ax,'-axis...'])
    [kspace,~] = sortHighOrderAcquisitions( ...
        rawFID,acq_manifest,ax,nPE,numSlices);

    outputFile = fullfile(outputPath, ...
        ['measurement_H_fast_HighOrder_',ax,'.mat']);

    saveHighOrderMeasurement(outputFile,kspace,dwelltime,FOV,ax,PosSlices);
    outputFiles.(ax) = outputFile;
end

disp('High-Order measurement preparation finished.')

end
