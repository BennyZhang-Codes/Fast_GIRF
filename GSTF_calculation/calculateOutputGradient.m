function [ out_signals, dwelltime ] = calculateOutputGradient( file, numRepPerGIRF, numGIRF, numTriang, numDelays, calcChannels, singleCoil, skipCoils)
% This function calculates the measured output gradient time courses
% from the phases of the recorded FID signals
% of the measurements with positive and negative gradient triangles
% Arguments:
%   file:          path to the .mat file containing the raw data
%   numRepPerGIRF: number of measurement repetitions to be averaged
%   numGIRF:       if numRepetitions/numRepPerGIRF > 1, which averaged dataset should be evaluated
%   numTriang:     number of triangles measured per delay
%   numDelays:     number of delays used between excitation and test gradient
%   calcChannels:  number of basis functions to be used for the field expansion
%   singleCoil:    number of the coil element to be evaluated if a combination is not wanted
%   skipCoils:     which coil elements to leave out of the calculation
%
% Returns:
%   out_signals: measured gradient waveforms
%   dwelltime:   dwell time in seconds

% Copyright (c) 2022 Hannah Scholten

%% Read raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawdata = load(file);
requiredVariables = {'dwelltime','FOV','orientation','numSlices','PosSlices','kspace'};
for i=1:length(requiredVariables)
    if ~isfield(rawdata,requiredVariables{i})
        error('Measurement file is missing required variable: %s',requiredVariables{i});
    end
end

dwelltime = rawdata.dwelltime;
FOV = rawdata.FOV;
orientation = rawdata.orientation;
numSlices = rawdata.numSlices;
PosSlices = rawdata.PosSlices;
kspace = rawdata.kspace;

% dimensions: [read-out-points, coils, PE1, PE2, slices,
%              polarity/triangles, cardiac(1), ADCs, repetitions]
weigh_equal = 0;

%% Prepare and validate raw data array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(kspace)==6
    % MATLAB removes trailing singleton dimensions. Explicitly restore the
    % Fast_GIRF layout so size(kspace,7:9) remains meaningful.
    kspace = reshape(kspace,[size(kspace,1),size(kspace,2),size(kspace,3), ...
        size(kspace,4),size(kspace,5),size(kspace,6),1,1,1]);
end

if size(kspace,3) ~= size(kspace,4)
    error('The two High-Order PE dimensions must have equal size.');
end
if size(kspace,5) ~= numSlices || length(PosSlices) ~= numSlices
    error('k-space slice dimension, numSlices, and PosSlices are inconsistent.');
end
if size(kspace,6) ~= 2*numTriang*numDelays
    error(['Measurement contains %d polarity/triangle entries, but processing ' ...
        'expects 2*numTriang*numDelays = %d.'],size(kspace,6),2*numTriang*numDelays);
end
if numRepPerGIRF < 1 || round(numRepPerGIRF) ~= numRepPerGIRF
    error('numRepPerGIRF must be a positive integer.');
end

numRep = size(kspace,9);
if mod(numRep,numRepPerGIRF) ~= 0
    error('Number of repetitions (%d) is not divisible by numRepPerGIRF (%d).', ...
        numRep,numRepPerGIRF);
end
numIter = numRep/numRepPerGIRF;
if numGIRF < 1 || numGIRF > numIter || round(numGIRF) ~= numGIRF
    error('numGIRF=%d is outside the available averaged datasets 1:%d.',numGIRF,numIter);
end

numROP = size(kspace,1);
numPE = size(kspace,3);
numADC = size(kspace,8);

%% Fourier-transform k-space data along the two PE directions %%%%%%%%%%%%%
kspace = fft_1D(kspace,3);
kspace = fft_1D(kspace,4);

%% Sort out undesired coil elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for coil=size(kspace,2):-1:1
    if ismember(coil,skipCoils)
        kspace(:,coil,:,:,:,:,:,:,:) = [];
    end
end
if size(kspace,2) < 1
    error('All receive-coil elements were removed.');
end
if singleCoil > size(kspace,2)
    error('singleCoil exceeds the number of remaining receive channels.');
end

%% Get coil-combined magnitude and phase data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[diff_phase,magnitude] = combineCoils(kspace,dwelltime,weigh_equal,singleCoil,numRepPerGIRF);
% [numROP, numPE, numPE, numSlices, triangles, ADCs, repetitions]
clearvars kspace;

%% Average magnitude data over numRepPerGIRF repetitions %%%%%%%%%%%%%%%%%%
magnitude_avg = zeros(numROP,numPE,numPE,numSlices,2*numTriang*numDelays,numADC,numIter);
for i=1:numIter
    magnitude_avg(:,:,:,:,:,:,i) = mean( ...
        magnitude(:,:,:,:,:,:,(i-1)*numRepPerGIRF+1:i*numRepPerGIRF),7);
end
magnitude_avg = magnitude_avg(:,:,:,:,:,:,numGIRF);
clearvars magnitude;

%% Separate positive and negative probing gradients %%%%%%%%%%%%%%%%%%%%%%%
diff_phase_plus = diff_phase(:,:,:,:,1:2:end-1,:,:);
diff_phase_minus = diff_phase(:,:,:,:,2:2:end,:,:);
clearvars diff_phase;

magnitude_plus = magnitude_avg(:,:,:,:,1:2:end-1,:);
clearvars magnitude_avg;

%% Calculate polarity difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_diff_phase = (diff_phase_plus-diff_phase_minus)/2;
clearvars diff_phase_plus diff_phase_minus;

% Average the phase derivative only after polarity subtraction.
delta_diff_phase_avg = zeros(numROP,numPE,numPE,numSlices,numTriang*numDelays,numADC,numIter);
for i=1:numIter
    delta_diff_phase_avg(:,:,:,:,:,:,i) = mean( ...
        delta_diff_phase(:,:,:,:,:,:,(i-1)*numRepPerGIRF+1:i*numRepPerGIRF),7);
end
clearvars delta_diff_phase;
delta_diff_phase_avg = delta_diff_phase_avg(:,:,:,:,:,:,numGIRF);

delta_diff_phase_avg = permute(delta_diff_phase_avg,[2,3,4,1,5,6,7]);
delta_diff_phase_avg = reshape(delta_diff_phase_avg, ...
    [numPE*numPE*numSlices,numROP,numTriang*numDelays,numADC]);

%% Get virtual-voxel positions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
positions = createPositionArray(orientation,numPE,numSlices,FOV,PosSlices);
positions = reshape(positions,[numPE*numPE*numSlices,3]);

%% Prepare magnitude and finite-data QC for voxel selection %%%%%%%%%%%%%%%
magnitude_plus = permute(magnitude_plus,[2,3,4,1,5,6]);
magnitude_plus = reshape(magnitude_plus, ...
    [numPE*numPE*numSlices,numROP,numTriang*numDelays,numADC]);

magnitude_gui = mean(magnitude_plus,3);
magnitude_gui = mean(magnitude_gui,4);
magnitude_gui = reshape(magnitude_gui,[numPE*numPE*numSlices,numROP]);

field2D = reshape(delta_diff_phase_avg,size(delta_diff_phase_avg,1),[]);
finiteVoxelMask = all(isfinite(magnitude_gui),2) & all(isfinite(field2D),2);
if numPE > 1
    fprintf('Finite High-Order voxels before selection: %d / %d\n', ...
        sum(finiteVoxelMask),length(finiteVoxelMask));
end

%% Sort out unusable voxels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validVoxels = true(size(positions,1),1);

if numPE > 1
    [selectionPath,selectionName,~] = fileparts(file);
    defaultSelectionFile = fullfile(selectionPath,[selectionName,'_voxelSelection.mat']);
    useSavedSelection = false;

    if exist(defaultSelectionFile,'file')
        savedSelection = load(defaultSelectionFile);
        savedPositionsOK = true;

        if isfield(savedSelection,'selection') && isfield(savedSelection.selection,'validVoxels')
            savedMask = logical(savedSelection.selection.validVoxels(:));
            if isfield(savedSelection.selection,'positions')
                savedPositions = savedSelection.selection.positions;
                savedPositionsOK = isequal(size(savedPositions),size(positions)) && ...
                    max(abs(savedPositions(:)-positions(:))) < 1e-9;
            end
        elseif isfield(savedSelection,'validVoxels')
            savedMask = logical(savedSelection.validVoxels(:));
        else
            savedMask = [];
        end

        savedFiniteOK = numel(savedMask)==numel(finiteVoxelMask) && ...
            ~any(savedMask & ~finiteVoxelMask);

        if numel(savedMask)==size(positions,1) && savedPositionsOK && savedFiniteOK
            savedA = createProbingMatrix(positions(savedMask,:),calcChannels);
            if sum(savedMask)>=calcChannels && size(savedA,2)==calcChannels && rank(savedA)==calcChannels
                validVoxels = savedMask;
                useSavedSelection = true;
                disp(['Loaded saved voxel selection: ',defaultSelectionFile]);
            else
                disp('Saved voxel selection is rank deficient for the requested SH expansion; opening GUI.');
            end
        else
            disp('Saved voxel selection does not match current geometry/data validity; opening GUI.');
        end
    end

    if ~useSavedSelection
        [validVoxels,voxelSelection] = selectVoxelsGUI(positions, ...
            magnitude_gui,delta_diff_phase_avg,numPE,numSlices,calcChannels, ...
            defaultSelectionFile,finiteVoxelMask);

        selection = voxelSelection; %#ok<NASGU>
        save(defaultSelectionFile,'selection','validVoxels');
        disp(['Saved voxel selection: ',defaultSelectionFile]);
    end
else
    % Legacy thin-slice behavior: retain the original 60% magnitude rule,
    % while preventing NaNs/Infs from contaminating the threshold.
    idx1 = min(10,numROP);
    idx2 = min(50,numROP);
    magnitude_metric = squeeze(mean(magnitude_plus(:,idx1:idx2,1,1),2));
    magnitude_metric(~isfinite(magnitude_metric)) = 0;
    max_mag = max(magnitude_metric);
    if size(positions,1)>2 && max_mag>0
        validVoxels = magnitude_metric >= max_mag*0.6;
    end
end

if ~any(validVoxels)
    error('No valid voxels remain for the field fit.');
end
positions = positions(validVoxels,:);
delta_diff_phase_avg = delta_diff_phase_avg(validVoxels,:,:,:);
delta_diff_phase_avg = reshape(delta_diff_phase_avg, ...
    [size(delta_diff_phase_avg,1),numROP*numTriang*numDelays*numADC]);
clearvars magnitude_plus magnitude_gui field2D;

%% Get probing matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(positions,1)<4
    probingMatrix = zeros(size(positions,1),calcChannels);
    for slice=1:size(positions,1)
        probingMatrix(slice,1) = 1;
        if calcChannels>1
            if strcmp(orientation,'dTra')
                slicePosition = positions(slice,3);
            elseif strcmp(orientation,'dSag')
                slicePosition = positions(slice,1);
            elseif strcmp(orientation,'dCor')
                slicePosition = positions(slice,2);
            else
                error('Unknown slice orientation: %s',orientation);
            end
            probingMatrix(slice,2) = slicePosition;
            if calcChannels>2
                probingMatrix(slice,3) = 2*slicePosition*slicePosition;
            end
        end
    end
else
    probingMatrix = createProbingMatrix(positions,calcChannels);
end

if size(probingMatrix,2) ~= calcChannels
    error('Probing matrix returned %d channels; requested %d.',size(probingMatrix,2),calcChannels);
end

%% Check spatial conditioning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rankA = rank(probingMatrix);
columnNorm = sqrt(sum(abs(probingMatrix).^2,1));
columnNorm(columnNorm==0) = 1;
probingMatrixNormalized = probingMatrix./columnNorm;

if rankA < calcChannels
    error(['Selected voxels are rank deficient for the requested field expansion. ' ...
        'rank(A)=%d, requested channels=%d.'],rankA,calcChannels);
end
condA = cond(probingMatrixNormalized);
disp(['Spatial probing matrix: ',num2str(size(probingMatrix,1)), ...
    ' voxels, rank ',num2str(rankA),'/',num2str(calcChannels), ...
    ', normalized cond(A) = ',num2str(condA)]);
clearvars probingMatrixNormalized columnNorm;

%% Calculate output field coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = 267.513*10^6; % rad/s/T
out_signals = 1/gamma*(probingMatrix\delta_diff_phase_avg);
out_signals = reshape(out_signals,[calcChannels,numROP,numTriang*numDelays,numADC]);
out_signals = permute(out_signals,[1,4,2,3]);
out_signals(~isfinite(out_signals)) = 0;
% [channels, ADC readouts, read-out-points, test gradients]

end
