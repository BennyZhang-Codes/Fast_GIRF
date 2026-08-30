% Regression test for High-Order self-term phase correction.
% This test does not require scanner data.

clc;

mfile_name = mfilename('fullpath');
[pathstr,~,~] = fileparts(mfile_name);
addpath(fullfile(pathstr,'..','GSTF_calculation'));

%% Selected self-term reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 10e-6;
lengthH = 101;
H = GIRF_matrix(dt,lengthH);

numChannels = 4;
H.gstf = ones(lengthH,numChannels);

% Add a pi offset to the Y channel so that channel 3 must be used as the
% phase reference. This distinguishes the High-Order helper from the legacy
% channel-2 reference used by GIRF_matrix.correct_GSTF_phase().
H.gstf(:,3) = -1;

corr_delay = 0;
referenceChannel = 3;
phaseIndex = floor(size(H.gstf,1)/2)+3;
phaseAtZero = angle(H.gstf(phaseIndex,referenceChannel));

correctHighOrderGSTFPhase(H,phaseAtZero,corr_delay,referenceChannel);

correctedPhase = angle(H.gstf(phaseIndex,referenceChannel));
phaseError = angle(exp(1i*(correctedPhase-phaseAtZero)));
assert(abs(phaseError) < 1e-12, ...
    'High-Order GSTF phase correction did not preserve the selected reference channel.');

%% Phase-wrap boundary regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H2 = GIRF_matrix(dt,lengthH);
H2.gstf = ones(lengthH,numChannels);

% These phases are physically only 0.02 rad apart, although direct
% subtraction would be close to 2*pi. The parity must therefore not switch.
wrappedCurrent = pi-0.01;
wrappedReference = -pi+0.01;
H2.gstf(:,referenceChannel) = exp(1i*wrappedCurrent);

linear_array = (1:1:lengthH)*pi;
expected = H2.gstf.*exp(1i*linear_array).';
correctHighOrderGSTFPhase(H2,wrappedReference,0,referenceChannel);

relError = norm(H2.gstf-expected,'fro') / norm(expected,'fro');
assert(relError < 1e-12, ...
    'Phase-wrap boundary incorrectly changed the High-Order phase-correction parity.');

disp('High-Order phase correction regression test passed.');
