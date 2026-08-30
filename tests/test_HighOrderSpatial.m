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
calcChannelsList = [4 9 16];
orientations = {'dSag','dCor','dTra'};

%% Run all three slice orientations and field-expansion orders %%%%%%%%%%%%
for i=1:length(orientations)
    orientation = orientations{i};

    positions = createPositionArray(orientation,numPE,numSlices,FOV,PosSlices);
    positions = reshape(positions,[numPE*numPE*numSlices,3]);

    for c=1:length(calcChannelsList)
        calcChannels = calcChannelsList(c);
        A = createProbingMatrix(positions,calcChannels);

        assert(size(A,2) == calcChannels, ...
            'Probing matrix has %d columns although calcChannels=%d.',size(A,2),calcChannels);

        rankA = rank(A);
        assert(rankA == calcChannels, ...
            'High-Order probing matrix is rank deficient for %s, %d channels.', ...
            orientation,calcChannels);

        columnNorm = sqrt(sum(abs(A).^2,1));
        columnNorm(columnNorm == 0) = 1;
        condA = cond(A ./ columnNorm);
        assert(condA < 10, ...
            'Default spatial geometry is unexpectedly ill-conditioned for %s, %d channels.', ...
            orientation,calcChannels);

        % Scale the synthetic coefficients approximately by spatial order
        % so that contributions of different orders have comparable size.
        coeffScale16 = [1, ...
            repmat(10,1,3), ...
            repmat(100,1,5), ...
            repmat(1000,1,7)].';
        coeffTrue = coeffScale16(1:calcChannels) .* (1:calcChannels).';

        fieldSynthetic = A * coeffTrue;
        coeffEstimate = A \ fieldSynthetic;
        relError = norm(coeffEstimate-coeffTrue) / norm(coeffTrue);

        fprintf('%s, %2d channels: rank(A) = %d/%d, normalized cond(A) = %.3g, relative error = %.3e\n', ...
            orientation,calcChannels,rankA,calcChannels,condA,relError);

        assert(relError < 1e-8, ...
            'Coefficient recovery failed for %s, %d channels.',orientation,calcChannels);
    end
end

disp('High-Order spatial regression test passed for 4, 9, and 16 channels in dSag, dCor, and dTra.');
