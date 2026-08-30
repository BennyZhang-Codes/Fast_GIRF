% Run all scanner-independent High-Order Fast_GIRF regression tests.

clear; clc;

mfile_name = mfilename('fullpath');
[pathstr,~,~] = fileparts(mfile_name);
repoPath = fullfile(pathstr,'..');

addpath(fullfile(repoPath,'GSTF_calculation'));
addpath(pathstr);

%% High-Order spatial basis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('==============================================================')
disp('Test 1/3: High-Order spatial spherical-harmonic recovery')
disp('==============================================================')
run(fullfile(pathstr,'test_HighOrderSpatial.m'));

%% Manifest-based acquisition sorting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('==============================================================')
disp('Test 2/3: High-Order acquisition manifest sorting')
disp('==============================================================')
run(fullfile(pathstr,'test_HighOrderSorting.m'));

%% High-Order phase correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('==============================================================')
disp('Test 3/3: High-Order self-term phase correction')
disp('==============================================================')
run(fullfile(pathstr,'test_HighOrderPhaseCorrection.m'));

%% Finished %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('==============================================================')
disp('All High-Order Fast_GIRF scanner-independent tests passed.')
disp('==============================================================')
