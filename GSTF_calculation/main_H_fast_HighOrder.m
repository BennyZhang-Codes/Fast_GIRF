% Copyright (c) 2022 Hannah Scholten
% High-order spatial extension for Fast_GIRF

clear all;
% Change current directory to that of this .m file
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);

t_tic_1 = tic; % For measuring the time

%% Select data files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = '../triangle_measurements/';     % Specify the path to the folder where the measurement data are stored

meas_name = 'measurement_H_fast_HighOrder_x.mat';  % x-axis
% meas_name = 'measurement_H_fast_HighOrder_y.mat';  % y-axis
% meas_name = 'measurement_H_fast_HighOrder_z.mat';  % z-axis
files = dir([path,meas_name]);
measfiles = {files.name};

ax_name = 'x';  % Specify the axis corresponding to the measurement
% ax_name = 'y';
% ax_name = 'z';

input_path = path;
input_name = 'input_H_fast_HighOrder.mat';
files = dir([input_path,input_name]);
inputfiles = {files.name};

% IMPORTANT: Check that the measurement files and the input files correspond to each other!!!
for i=1:1:length(measfiles)
    disp(['measurement file #',num2str(i),': ',measfiles{i}])
    disp(['      input file #',num2str(i),': ',inputfiles{i}])
    disp(' ')
end

if isempty(measfiles)
    error('No High-Order measurement file found. Check path and meas_name.');
end
if isempty(inputfiles)
    error('No High-Order input file found. Check input_path and input_name.');
end

doSaveGIRFs = 1;
path2save = './results/';
name2save = ['HighOrder_',meas_name(13:end)];

%% Set high-order field expansion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 channels  = B0 + first order
% 9 channels  = B0 + first + second order
% 16 channels = B0 + first + second + third order
calcChannels = 16;

if ~ismember(calcChannels,[4,9,16])
    error('calcChannels must be 4, 9, or 16 for spatial High-Order processing.');
end

% The self-term depends on the physical test-gradient axis because the
% spatial basis is always ordered as [B0, X, Y, Z, ...].
if strcmp(ax_name,'x')
    selfTerm = 2;
elseif strcmp(ax_name,'y')
    selfTerm = 3;
elseif strcmp(ax_name,'z')
    selfTerm = 4;
else
    error('Unknown ax_name. Use x, y, or z.');
end
term2plot = selfTerm;

%% Set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prep = preparation;
% The preparation object will hold the input and measured output data and
% associated parameters necessary for preparing the data for the GIRF calculation.
prep.numRepPerGIRF = 1;         % number of acquired repetitions/averages used for one GIRF
prep.numGIRF = 1;               % which averaged dataset should be evaluated
prep.numADC = 1;                % number of readouts per TR to evaluate
prep.numTriang = 1;             % number of triangles per delay
prep.numTriangs4fft = 2;        % number of triangles used for H_FFT
prep.numDelays = 7;             % number of Fast-GIRF delays
prep.skipTriangles = cell(1,length(measfiles));
prep.singleCoil = 0;
prep.skipCoils = [];
prep.resamp_factors = [10];     % time resolution for matrix calculation in us
prep.cut_time = 80;             % points cut at the beginning/end of each measurement period
prep.VarPreMeas = 0;
prep.calcChannels = calcChannels;

% Some more parameters needed for the GIRF-calculation
corr_delay = 0.95e-4;       % compensate phase error introduced by matrix construction/resampling
lambda = sqrt(1e-7);        % Tikhonov regularization weighting factor
A = 600;                    % frequency-dependent regularization growth factor
omega_0 = 2500;             % cutoff between LF and HF GSTF [Hz]

%% Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Prepare High-Order input and output data...')
% For spatial data (numPE > 1), calculateOutputGradient() reuses a saved
% voxel selection when available and compatible. Otherwise it opens the
% voxel-selection GUI and automatically saves the confirmed selection.
prep.prepare_Data(path,measfiles,input_path,inputfiles);
dts_out = prep.dts_resamp;
prep.findMax_tshift();
% Dimensions:
% prep.input{k}: [numGRTTimePoints, numTriang]
% prep.output{k}: [calcChannels, numADC, numGRTTimePoints, numTriang]
% prep.input4fft: [numADCTimePoints, numTriang]
% prep.output4fft: [calcChannels, numADCTimePoints, numTriang]

% Exclude the final part of scheme-i measurements because the same time
% window is covered by the measurement with the 1-ms delayed excitation.
prep.discardData{1}(:,1,1000:end,1:2) = 0;
disp('    High-Order data prepared.')

%% Calculate H_FFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculate H_FFT...')
H_fft = GIRF_fft(1e-6);
H_fft.calcH_fft(prep.input4fft,prep.output4fft,prep.skipTriangles{1})
disp('    H_FFT calculated.')

%% Put input and output data into matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Prepare matrices for time domain calculation...')
t_tic_2 = tic;

numInputChannels = 1;
lengthADC = size(prep.output{1},3);
dt_in = dts_out(1);

lengthH = floor((prep.findMax_tshift() - prep.t_shift{1}(1,1))/dt_in + ...
    lengthADC - 990/prep.resamp_factors(1));
if ~mod(lengthH,2)
    lengthH = lengthH - 1;
end

inputMatrix_LF = prep.calcInputMatrix(numInputChannels,lengthH,[],[]);
outputMatrix_LF = prep.calcOutputMatrix([],[]);

lengthH_HF = floor(lengthADC - 990/prep.resamp_factors(1)) - 200;
if ~mod(lengthH_HF,2)
    lengthH_HF = lengthH_HF - 1;
end
inputMatrix_HF = prep.calcInputMatrix_firstMeas(numInputChannels,lengthH_HF,0,prep.numTriangs4fft/prep.numTriang);
outputMatrix_HF = prep.calcOutputMatrix_firstMeas(0,prep.numTriangs4fft/prep.numTriang);

disp(['    size(inputMatrix_LF) = ',num2str(size(inputMatrix_LF))])
disp(['    size(outputMatrix_LF) = ',num2str(size(outputMatrix_LF))])
disp(['    size(inputMatrix_HF) = ',num2str(size(inputMatrix_HF))])
disp(['    size(outputMatrix_HF) = ',num2str(size(outputMatrix_HF))])
disp('    Matrices prepared.')

%% Calculate GIRFs in the time domain with the matrix inversion method %%%%
disp('Calculate High-Order H_LF and H_HF in time domain...')
H_LF = GIRF_matrix(dt_in,lengthH);
H_HF = GIRF_matrix(dt_in,lengthH_HF);

H_LF.calcH_matrix_Tikhonov(inputMatrix_LF,outputMatrix_LF,lambda,0);
alpha_array = alpha_func(H_HF.f_axis,omega_0,H_HF.f_axis(end),A);
H_HF.calcH_matrix_Tikhonov_freqWeight(inputMatrix_HF,outputMatrix_HF,lambda,alpha_array);

%% Phase correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the physical self-term as the phase reference. For the full spatial
% basis this is X=2, Y=3, or Z=4 depending on the driven gradient axis.
% Keep the original GIRF_matrix.m unchanged and apply the High-Order
% reference-channel extension through a dedicated helper.
phaseAtZero_fft = angle(H_fft.gstf(floor(size(H_fft.gstf,1)/2)+3,selfTerm));
correctHighOrderGSTFPhase(H_LF,phaseAtZero_fft,corr_delay,selfTerm);
correctHighOrderGSTFPhase(H_HF,phaseAtZero_fft,corr_delay,selfTerm);

disp('    High-Order H in time domain calculated.')

%% Dwell time compensation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Perform dwelltime compensation...')
H_fft.dwelltime_compensation(prep.dwelltime);
H_LF.dwelltime_compensation(prep.dwelltime);
H_HF.dwelltime_compensation(prep.dwelltime);
disp('    Dwelltime compensation done.');

%% Combine H_LF and H_HF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Combine High-Order GSTFs...')
H_combined = GIRF_combined();
H_combined.combineGSTFs_cutoffFreq(H_LF.gstf,H_LF.f_axis,H_HF.gstf, ...
    H_HF.f_axis,omega_0,H_LF.fieldOffsets,corr_delay);
disp('    High-Order GSTFs combined.')
t_toc_2 = toc(t_tic_2);

%% Save GIRFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if doSaveGIRFs
    save([path2save,name2save],'H_combined','lambda','alpha_array','H_LF', ...
        'H_HF','H_fft','dts_out','corr_delay','omega_0','calcChannels','selfTerm');
    disp('Saved High-Order results.');
else
    disp('High-Order results were not saved.');
end
t_toc_1 = toc(t_tic_1);

%% Plot GSTFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plotting...')

% Keep the original Fast_GIRF comparison plot for the physical self-term.
plotter = GIRF_plotter();
plotter.plot_GSTFs(H_fft,['H^f^a^s^t_F_F_T_,_',ax_name], ...
    H_LF,['H^f^a^s^t_L_F_,_',ax_name], ...
    H_HF,['H^f^a^s^t_H_F_,_',ax_name], ...
    H_combined,['H^f^a^s^t_,_',ax_name],4,term2plot,ax_name);

% Plot all available spatial terms grouped by SH order. Different spatial
% orders have different GSTF units and are therefore shown separately.
plotHighOrderGSTFs(H_combined,ax_name,5);

disp('    main_H_fast_HighOrder.m finished.')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alpha = alpha_func(omega,omega_0,omega_max,A)
    alpha = zeros(size(omega));
    for i=1:size(omega,2)
        if abs(omega(i))<omega_0
            alpha(i) = 0;
        else
            alpha(i) = A*sqrt((abs(omega(i))-omega_0)/(omega_max-omega_0));
        end
    end
end