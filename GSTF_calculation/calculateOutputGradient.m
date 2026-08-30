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
%   skipCoils:     which coil elements to leave out of the calculation (if they have too low signal, for example)
%
% Returns:
%   out_signals: measured gradient waveforms
%   dwelltime:   dwell time in seconds

% Copyright (c) 2022 Hannah Scholten

%% Read raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawdata = load(file);

dwelltime = rawdata.dwelltime;      % dwelltime in seconds
FOV = rawdata.FOV;                  % FOV width and height in mm
orientation = rawdata.orientation;  % slice orientation
numSlices = rawdata.numSlices;      % number of slices
PosSlices = rawdata.PosSlices;      % slice positions (order should match the raw data)

kspace = rawdata.kspace;            % FID raw data
% dimensions: [read-out-points, coils, PE-steps(PE-dir), PE-steps(RO-dir), slices, triangles, Cardiac-Phases(1), ADC readouts, repetitions(if>1)]
% If phase encoding was applied in the two directions perpendicular to the
% slice selection direction to enable the determination of higher order
% terms, the respective data should be stored along the dimensions
% "PE-steps(PE-dir)" and "PE-steps(RO-dir)" in the raw data array.
% The following code assumes that the number of PE-steps in both directions
% is equal.
% Since every triangle is played out twice (every second time with inverted
% sign), size(kspace,6) = 2*numTriang*numDelays.

weigh_equal = 0;                    % parameter for the coil weighting

%% Prepare raw data array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (length(size(kspace))==6)
    % In case of just 1 ADC and 1 measurement: add extra dimensions, so the dimensions in the following are still correct
    kspace=repmat(kspace,1,1,1,1,1,1,1,1,numRepPerGIRF);
end
clearvars rawdata;

%% Determine some parameters for easier array reshaping later on %%%%%%%%%%
numRep = size(kspace,9);                % number of repetitions
numIter = floor(numRep/numRepPerGIRF);  % number of loop-counts
numROP = size(kspace,1);                % number of Read Out Points
numPE = size(kspace,3);                 % number of phase encode steps
numADC = size(kspace,8);                % number of ADCs per TR

%% Fourier-transform k-space data along the two phase-encoding directions %
kspace = fft_1D(kspace, 3);
kspace = fft_1D(kspace, 4);

%% Sort out undesired coil elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for coil=size(kspace,2):-1:1
    if ismember(coil, skipCoils)
        kspace(:,coil,:,:,:,:,:,:,:) = [];
    end
end

%% Get coil-combined magnitude and phase data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[diff_phase, magnitude] = combineCoils(kspace, dwelltime, weigh_equal, singleCoil, numRepPerGIRF);
% dimensions: [numROP, numPE, numPE, numSlices, triangles, ADCs, repetitions]
clearvars kspace;

%% Average magnitude data over numRepPerGIRF repetitions %%%%%%%%%%%%%%%%%%
% This is not done for the phase data because the slope of the measured phase changes over time,
% so we have to subtract the reference data before averaging over numRepPerGIRF
magnitude_avg = zeros(numROP, numPE, numPE, numSlices, 2*numTriang*numDelays, numADC, numIter);
for i=1:1:numIter
    magnitude_avg(:,:,:,:,:,:,i) = mean(magnitude(:,:,:,:,:,:,(i-1)*numRepPerGIRF+1:i*numRepPerGIRF),7);
end
magnitude_avg = magnitude_avg(:,:,:,:,:,:,numGIRF);
clearvars magnitude;

%% Separate triangles with positive and negative sign %%%%%%%%%%%%%%%%%%%%%
diff_phase_plus = diff_phase(:,:,:,:,1:2:end-1,:,:); % [numROP, numPE, numPE, numSlices, numTriang, numADC, numRep]
diff_phase_minus = diff_phase(:,:,:,:,2:2:end,:,:);
clearvars diff_phase;

magnitude_plus = magnitude_avg(:,:,:,:,1:2:end-1,:); % [numROP, numPE, numPE, numSlices, numTriang, numADC]
clearvars magnitude_avg;

%% Calculate difference between positive and negative triangles %%%%%%%%%%%
delta_diff_phase = (diff_phase_plus - diff_phase_minus)/2;
clearvars diff_phase_plus diff_phase_minus;

% Average over numRepPerGIRF repetitions
delta_diff_phase_avg = zeros(numROP, numPE, numPE, numSlices, numTriang*numDelays, numADC, numIter);
for i=1:1:numIter
    delta_diff_phase_avg(:,:,:,:,:,:,i) = mean(delta_diff_phase(:,:,:,:,:,:,(i-1)*numRepPerGIRF+1:i*numRepPerGIRF),7);
end
clearvars delta_diff_phase;

% Take the specified iteration for the GIRF
delta_diff_phase_avg = delta_diff_phase_avg(:,:,:,:,:,:,numGIRF);

delta_diff_phase_avg = permute(delta_diff_phase_avg, [2,3,4,1,5,6,7]);
% [numPE, numPE, numSlices, numROP, numTriang*numDelays, numADC, 1]
delta_diff_phase_avg = reshape(delta_diff_phase_avg, ...
    [numPE*numPE*numSlices, numROP, numTriang*numDelays, numADC]);

%% Get the positions of the measured voxels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
positions = createPositionArray(orientation, numPE, numSlices, FOV, PosSlices);
positions = reshape(positions, [numPE*numPE*numSlices, 3]);

%% Prepare magnitude data for voxel selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
magnitude_plus = permute(magnitude_plus, [2,3,4,1,5,6]);
magnitude_plus = reshape(magnitude_plus, ...
    [numPE*numPE*numSlices, numROP, numTriang*numDelays, numADC]);

% Mean magnitude time course across the positive-polarity test gradients.
% This is used for display in the High-Order voxel-selection GUI.
magnitude_gui = mean(magnitude_plus,3);
magnitude_gui = mean(magnitude_gui,4);
magnitude_gui = reshape(magnitude_gui,[numPE*numPE*numSlices,numROP]);

%% Sort out unusable voxels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validVoxels = true(size(positions,1),1);

if numPE > 1
    % Spatially resolved High-Order acquisition. Reuse a saved selection
    % automatically when it is compatible and supports the requested SH
    % expansion; otherwise open the interactive voxel-selection GUI.
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

        if numel(savedMask) == size(positions,1) && savedPositionsOK
            savedA = createProbingMatrix(positions(savedMask,:),calcChannels);
            if sum(savedMask) >= calcChannels && rank(savedA) == calcChannels
                validVoxels = savedMask;
                useSavedSelection = true;
                disp(['Loaded saved voxel selection: ',defaultSelectionFile]);
            else
                disp('Saved voxel selection is rank deficient for the requested SH expansion; opening GUI.');
            end
        else
            disp('Saved voxel selection does not match the current spatial geometry; opening GUI.');
        end
    end

    if ~useSavedSelection
        [validVoxels, voxelSelection] = selectVoxelsGUI(positions, ...
            magnitude_gui,delta_diff_phase_avg,numPE,numSlices,calcChannels, ...
            defaultSelectionFile);

        % Save the confirmed selection automatically so that subsequent
        % processing of the same dataset is reproducible and non-interactive.
        selection = voxelSelection; %#ok<NASGU>
        save(defaultSelectionFile,'selection','validVoxels');
        disp(['Saved voxel selection: ',defaultSelectionFile]);
    end
else
    % Legacy thin-slice behavior. Keep the original 60% magnitude
    % threshold unchanged for backwards compatibility.
    idx1 = min(10,numROP);
    idx2 = min(50,numROP);
    magnitude_metric = squeeze(mean(magnitude_plus(:,idx1:idx2,1,1),2));
    max_mag = max(magnitude_metric);

    if size(positions,1)>2 && max_mag>0
        validVoxels = magnitude_metric >= max_mag*0.6;
    end
end

positions = positions(validVoxels,:);
delta_diff_phase_avg = delta_diff_phase_avg(validVoxels,:,:,:);
delta_diff_phase_avg = reshape(delta_diff_phase_avg, ...
    [size(delta_diff_phase_avg,1), numROP*numTriang*numDelays*numADC]);

clearvars magnitude_plus magnitude_gui;

%% Get the probing matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to translate voxel positions into spherical harmonics
if size(positions,1)<4
    probingMatrix = zeros(size(positions,1), calcChannels);
    for slice=1:1:size(positions,1)
        probingMatrix(slice,1) = 1;
        if calcChannels>1
            if strcmp(orientation,'dTra')
                slicePosition = positions(slice,3);
            elseif strcmp(orientation,'dSag')
                slicePosition = positions(slice,1);
            elseif strcmp(orientation,'dCor')
                slicePosition = positions(slice,2);
            end
            probingMatrix(slice,2) = slicePosition;
            if calcChannels>2
                probingMatrix(slice,3) = 2*slicePosition*slicePosition;
            end
        end
    end
else
    probingMatrix = createProbingMatrix(positions,calcChannels); % [numValidVoxels, 4/9/16]
end

%% Check spatial conditioning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rankA = rank(probingMatrix);
columnNorm = sqrt(sum(abs(probingMatrix).^2,1));
columnNorm(columnNorm == 0) = 1;
probingMatrixNormalized = probingMatrix ./ columnNorm;

if rankA < calcChannels
    error(['The selected voxels are rank deficient for the requested field expansion. ' ...
        'rank(A) = %d, requested channels = %d.'],rankA,calcChannels);
end

condA = cond(probingMatrixNormalized);
disp(['Spatial probing matrix: ',num2str(size(probingMatrix,1)), ...
    ' voxels, rank ',num2str(rankA),'/',num2str(calcChannels), ...
    ', normalized cond(A) = ',num2str(condA)]);
clearvars probingMatrixNormalized columnNorm;

%% Calculate the output signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = 267.513*10^6; % rad/s/T
out_signals = 1/gamma * (probingMatrix\delta_diff_phase_avg); % [channels, numTimePoints]
out_signals = reshape(out_signals, [calcChannels, numROP, numTriang*numDelays, numADC]);
out_signals = permute(out_signals, [1,4,2,3]);
out_signals(isnan(out_signals)) = 0;
% Dimensions of out_signals: [channels, ADC readouts, read-out-points, test gradients]

end