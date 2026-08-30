function [ outputFiles, Protocol ] = prepareHighOrderTwixFromSequence(datFile, manifestFile, seqFile, outputPath, firstAcquisition)
% Prepare High-Order Fast_GIRF Twix data using parameters from the .seq file.
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
% This wrapper is preferred over manually entering dwell time, FOV, slice
% positions, and nPE because those quantities are read from the exact
% sequence file used on the scanner.

if nargin < 5
    firstAcquisition = [];
end
if nargin < 4 || isempty(outputPath)
    outputPath = '.';
end

%% Read protocol from the exact Pulseq file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Protocol = readHighOrderSequenceDefinitions(seqFile);

%% Check manifest length against sequence definitions %%%%%%%%%%%%%%%%%%%%%
manifestData = load(manifestFile);
if ~isfield(manifestData,'acq_manifest')
    error('Manifest file does not contain acq_manifest.');
end

numManifestAcquisitions = height(manifestData.acq_manifest);
if numManifestAcquisitions ~= Protocol.expectedAcquisitions
    error(['Manifest/sequence mismatch: sequence expects %d ADCs but ' ...
        'manifest contains %d rows. Do not mix files generated from ' ...
        'different High-Order protocols.'], ...
        Protocol.expectedAcquisitions,numManifestAcquisitions);
end

%% Read sequential image ADCs from Twix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rawFID,~] = readHighOrderTwix(datFile,manifestData.acq_manifest,firstAcquisition);

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
outputFiles = prepareHighOrderMeasurements(rawFID,manifestData.acq_manifest, ...
    outputPath,Protocol.adcDwell,Protocol.PEFOV_mm, ...
    Protocol.sliceOffsets_mm,Protocol.nPE);

fprintf('High-Order Twix preparation completed using the exact .seq protocol.\n');

end
