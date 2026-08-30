% Synthetic regression test for the High-Order spatial field expansion.
% This test does not require scanner data.

clc;

mfile_name = mfilename('fullpath');
[pathstr,~,~] = fileparts(mfile_name);
addpath(fullfile(pathstr,'..','GSTF_calculation'));

%% Test parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FOV = 200;                         % mm
numPE = 5;
numSlices = 4;
PosSlices = [-34 -17 17 34];     % mm
calcChannels = 16;
orientations = {'dSag','dCor','dTra'};

% Scale the synthetic coefficients approximately by spatial order so that
% the contributions of different SH orders have comparable magnitude.
coeffScale = [1, ...
    repmat(10,1,3), ...
    repmat(100,1,5), ...
    repmat(1000,1,7)].';
coeffTrue = coeffScale .* (1:calcChannels).';

%% Run all three slice orientations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(orientations)
    orientation = orientations{i};

    positions = createPositionArray(orientation,numPE,numSlices,FOV,PosSlices);
    positions = reshape(positions,[numPE*numPE*numSlices,3]);
    A = createProbingMatrix(positions,calcChannels);

    rankA = rank(A);
    assert(rankA == calcChannels, ...
        'High-Order probing matrix is rank deficient for %s.',orientation);

    columnNorm = sqrt(sum(abs(A).^2,1));
    columnNorm(columnNorm == 0) = 1;
    condA = cond(A ./ columnNorm);
    assert(condA < 10, ...
        'Default High-Order spatial geometry is unexpectedly ill-conditioned for %s.',orientation);

    fieldSynthetic = A * coeffTrue;
    coeffEstimate = A \ fieldSynthetic;
    relError = norm(coeffEstimate-coeffTrue) / norm(coeffTrue);

    fprintf('%s: rank(A) = %d/%d, normalized cond(A) = %.3g, relative error = %.3e\n', ...
        orientation,rankA,calcChannels,condA,relError);

    assert(relError < 1e-8, ...
        'High-Order coefficient recovery failed for %s.',orientation);
end

disp('High-Order spatial regression test passed for dSag, dCor, and dTra.');
