% Copyright (c) 2022 Hannah Scholten
% High-order spatial extension for Fast_GIRF

clear all;
% Change current directory to that of this .m file
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name); %#ok<ASGLU>
cd(pathstr);

t_tic_1 = tic; % For measuring the time

%% Select data files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = '../triangle_measurements/';

meas_name = 'measurement_H_fast_HighOrder_x.mat';
% meas_name = 'measurement_H_fast_HighOrder_y.mat';
% meas_name = 'measurement_H_fast_HighOrder_z.mat';
files = dir(fullfile(path,meas_name));
measfiles = {files.name};

% Derive the physical input axis from the measurement file name instead of
% maintaining a second manually edited ax_name variable.
axisToken = regexp(meas_name,'_([xyz])\.mat$','tokens','once');
if isempty(axisToken)
    error('Could not determine x/y/z axis from measurement file name: %s',meas_name);
end
ax_name = axisToken{1};

input_path = path;
input_name = 'input_H_fast_HighOrder.mat';
files = dir(fullfile(input_path,input_name));
inputfiles = {files.name};

if isempty(measfiles)
    error('No High-Order measurement file found. Check path and meas_name.');
end
if length(measfiles) ~= 1
    error('Expected exactly one High-Order measurement file for the selected axis.');
end
if isempty(inputfiles)
    error(['No input_H_fast_HighOrder.mat found in the measurement folder. ' ...
        'Use prepareHighOrderTwixFromSequence() to copy the exact nominal input.']);
end
if length(inputfiles) ~= 1
    error('Expected exactly one High-Order nominal input file.');
end

% Verify that the measurement metadata agrees with the axis encoded in its
% file name before any GIRF calculation is performed.
measurementMeta = load(fullfile(path,measfiles{1}),'orientation');
if ~isfield(measurementMeta,'orientation')
    error('Measurement file does not contain orientation metadata.');
end
if strcmp(ax_name,'x')
    expectedOrientation = 'dSag';
elseif strcmp(ax_name,'y')
    expectedOrientation = 'dCor';
else
    expectedOrientation = 'dTra';
end
if ~strcmp(measurementMeta.orientation,expectedOrientation)
    error(['Measurement/axis mismatch: %s implies %s, but the file contains %s.'], ...
        ax_name,expectedOrientation,measurementMeta.orientation);
end

fprintf('measurement file: %s\n',measfiles{1});
fprintf('input file      : %s\n\n',inputfiles{1});

doSaveGIRFs = 1;
path2save = './results/';
if ~exist(path2save,'dir')
    mkdir(path2save);
end
name2save = ['HighOrder_GSTF_',ax_name,'.mat'];

%% Set high-order field expansion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 channels  = B0 + first order
% 9 channels  = B0 + first + second order
% 16 channels = B0 + first + second + third order
calcChannels = 16;

if ~ismember(calcChannels,[4,9,16])
    error('calcChannels must be 4, 9, or 16 for spatial High-Order processing.');
end

% Spatial basis is always [B0, X, Y, Z, ...].
if strcmp(ax_name,'x')
    selfTerm = 2;
elseif strcmp(ax_name,'y')
    selfTerm = 3;
else
    selfTerm = 4;
end
term2plot = selfTerm;

%% Set parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prep = preparation;
prep.numRepPerGIRF = 1;
prep.numGIRF = 1;
prep.numADC = 1;
prep.numTriang = 1;
prep.numTriangs4fft = 2;
prep.numDelays = 7;
prep.skipTriangles = cell(1,length(measfiles));
prep.singleCoil = 0;
prep.skipCoils = [];
prep.resamp_factors = [10];
prep.cut_time = 80;
prep.VarPreMeas = 0;
prep.calcChannels = calcChannels;

% Fast_GIRF regularization parameters retained from the original workflow.
corr_delay = 0.95e-4;
lambda = sqrt(1e-7);
A = 600;
omega_0 = 2500;

%% Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Prepare High-Order input and output data...')
prep.prepare_Data(path,measfiles,input_path,inputfiles);
dts_out = prep.dts_resamp;
prep.findMax_tshift();

% Exclude the final part of scheme-i measurements because the same time
% window is covered by the measurement with 1-ms delayed excitation.
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

lengthH = floor((prep.findMax_tshift()-prep.t_shift{1}(1,1))/dt_in + ...
    lengthADC-990/prep.resamp_factors(1));
if ~mod(lengthH,2)
    lengthH = lengthH-1;
end

inputMatrix_LF = prep.calcInputMatrix(numInputChannels,lengthH,[],[]);
outputMatrix_LF = prep.calcOutputMatrix([],[]);

lengthH_HF = floor(lengthADC-990/prep.resamp_factors(1))-200;
if ~mod(lengthH_HF,2)
    lengthH_HF = lengthH_HF-1;
end
inputMatrix_HF = prep.calcInputMatrix_firstMeas(numInputChannels,lengthH_HF,0,prep.numTriangs4fft/prep.numTriang);
outputMatrix_HF = prep.calcOutputMatrix_firstMeas(0,prep.numTriangs4fft/prep.numTriang);

disp(['    size(inputMatrix_LF) = ',num2str(size(inputMatrix_LF))])
disp(['    size(outputMatrix_LF) = ',num2str(size(outputMatrix_LF))])
disp(['    size(inputMatrix_HF) = ',num2str(size(inputMatrix_HF))])
disp(['    size(outputMatrix_HF) = ',num2str(size(outputMatrix_HF))])
disp('    Matrices prepared.')

%% Calculate GIRFs in time domain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculate High-Order H_LF and H_HF in time domain...')
H_LF = GIRF_matrix(dt_in,lengthH);
H_HF = GIRF_matrix(dt_in,lengthH_HF);

H_LF.calcH_matrix_Tikhonov(inputMatrix_LF,outputMatrix_LF,lambda,0);
alpha_array = alpha_func(H_HF.f_axis,omega_0,H_HF.f_axis(end),A);
H_HF.calcH_matrix_Tikhonov_freqWeight(inputMatrix_HF,outputMatrix_HF,lambda,alpha_array);

%% Phase correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retain the original Fast_GIRF phase-reference frequency-bin convention,
% but use the physical self-term (X/Y/Z) as the reference channel.
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
    save(fullfile(path2save,name2save),'H_combined','lambda','alpha_array','H_LF', ...
        'H_HF','H_fft','dts_out','corr_delay','omega_0','calcChannels','selfTerm','ax_name');
    disp('Saved High-Order results.');
else
    disp('High-Order results were not saved.');
end
t_toc_1 = toc(t_tic_1); %#ok<NASGU>

%% Plot GSTFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plotting...')
plotter = GIRF_plotter();
plotter.plot_GSTFs(H_fft,['H^f^a^s^t_F_F_T_,_',ax_name], ...
    H_LF,['H^f^a^s^t_L_F_,_',ax_name], ...
    H_HF,['H^f^a^s^t_H_F_,_',ax_name], ...
    H_combined,['H^f^a^s^t_,_',ax_name],4,term2plot,ax_name);
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
