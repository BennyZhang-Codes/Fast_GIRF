function [ outputFiles, Protocol ] = prepareHighOrderTwixFromSequence(datFile, manifestFile, seqFile, outputPath, firstAcquisition)
% Prepare High-Order Fast_GIRF Twix data using the exact acquisition files.
%
% Arguments:
%   datFile:          Siemens raw-data .dat file
%   manifestFile:     matching acq_manifest_HighOrder.mat
%   seqFile:          exact fast_gstf_HighOrder_*.seq used for acquisition
%   outputPath:       folder for measurement_H_fast_HighOrder_*.mat
%   firstAcquisition: optional first image ADC when Twix contains extra
%                     image acquisitions. Leave empty when counts match.
%
% Returns:
%   outputFiles: x/y/z measurement file paths
%   Protocol:    protocol parameters read from the .seq file
%
% The .seq file, acquisition manifest, nominal input_H_fast_HighOrder.mat,
% and Twix stream are cross-checked before any measurement files are saved.

if nargin < 5
    firstAcquisition = [];
end
if nargin < 4 || isempty(outputPath)
    outputPath = '.';
end
if ~exist(outputPath,'dir')
    mkdir(outputPath);
end

%% Read protocol from the exact Pulseq file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Protocol = readHighOrderSequenceDefinitions(seqFile);

%% Load and strictly validate acquisition manifest %%%%%%%%%%%%%%%%%%%%%%%%
manifestData = load(manifestFile);
if ~isfield(manifestData,'acq_manifest')
    error('Manifest file does not contain acq_manifest.');
end
acq_manifest = manifestData.acq_manifest;
validateHighOrderManifest(acq_manifest,Protocol);

%% Validate the nominal input generated with this sequence %%%%%%%%%%%%%%%%
[seqPath,~,~] = fileparts(seqFile);
if isempty(seqPath)
    seqPath = '.';
end
inputSource = fullfile(seqPath,'input_H_fast_HighOrder.mat');
if ~exist(inputSource,'file')
    error(['The nominal input file was not found next to the sequence: %s. ' ...
        'Keep the .seq, manifest, and input_H_fast_HighOrder.mat generated ' ...
        'in the same sequence-generation run.'],inputSource);
end

inputData = load(inputSource);
requiredInputVariables = {'grad_input','lengthADC','shift'};
for i=1:length(requiredInputVariables)
    if ~isfield(inputData,requiredInputVariables{i})
        error('Nominal input file is missing variable: %s',requiredInputVariables{i});
    end
end

if size(inputData.grad_input,2) ~= Protocol.numTriangles
    error('Nominal grad_input contains %d probes; sequence contains %d.', ...
        size(inputData.grad_input,2),Protocol.numTriangles);
end
if numel(inputData.shift) ~= Protocol.numTriangles
    error('Nominal input shift array does not match the number of sequence probes.');
end
expectedLengthADC_us = round(Protocol.adcSamples*Protocol.adcDwell/1e-6);
if double(inputData.lengthADC) ~= expectedLengthADC_us
    error('Nominal input lengthADC (%d us) does not match the sequence ADC duration (%d us).', ...
        double(inputData.lengthADC),expectedLengthADC_us);
end

% Keep the exact nominal input together with the sorted measurement files,
% which is where main_H_fast_HighOrder.m expects it by default.
inputDestination = fullfile(outputPath,'input_H_fast_HighOrder.mat');
if ~strcmp(inputSource,inputDestination)
    [copyOK,copyMessage] = copyfile(inputSource,inputDestination,'f');
    if ~copyOK
        error('Could not copy nominal input file to processing folder: %s',copyMessage);
    end
end

%% Read sequential image ADCs from Twix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rawFID,~] = readHighOrderTwix(datFile,acq_manifest,firstAcquisition);

if size(rawFID,1) ~= Protocol.adcSamples
    error(['Twix ADC sample count (%d) does not match the sequence ' ...
        'definition (%d). Check flagRemoveOS and confirm that the correct ' ...
        '.dat and .seq files are paired.'], ...
        size(rawFID,1),Protocol.adcSamples);
end

if size(rawFID,3) ~= Protocol.expectedAcquisitions
    error('Unexpected number of High-Order ADCs after Twix extraction.');
end

%% Manifest-based sorting and Fast_GIRF measurement files %%%%%%%%%%%%%%%%%
outputFiles = prepareHighOrderMeasurements(rawFID,acq_manifest, ...
    outputPath,Protocol.adcDwell,Protocol.PEFOV_mm, ...
    Protocol.sliceOffsets_mm,Protocol.nPE);
outputFiles.input = inputDestination;

fprintf('High-Order Twix preparation completed using the exact acquisition protocol.\n');

end
