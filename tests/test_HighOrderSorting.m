% Synthetic regression test for High-Order acquisition sorting.
% This test does not require scanner data.

clc;

mfile_name = mfilename('fullpath');
[pathstr,~,~] = fileparts(mfile_name);
addpath(fullfile(pathstr,'..','GSTF_calculation'));

%% Test parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numROP = 8;
numCoils = 2;
nPE = 5;
numSlices = 4;
numTriangles = 7;
numRep = 2;
axisName = 'x';

peIndices = -(nPE-1)/2:(nPE-1)/2;
[PE1Grid,PE2Grid] = ndgrid(1:nPE,1:nPE);
radius = (nPE-1)/2;
peMask = (peIndices(PE1Grid)/radius).^2 + ...
         (peIndices(PE2Grid)/radius).^2 <= 1 + eps;
[validPE1,validPE2] = find(peMask);
validPhaseEncodes = [validPE1,validPE2];
numPEAcq = size(validPhaseEncodes,1);

%% Build a synthetic manifest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numAcq = numRep * numTriangles * 2 * numPEAcq * numSlices;
Axis = repmat({axisName},numAcq,1);
Repetition = zeros(numAcq,1);
Triangle = zeros(numAcq,1);
Polarity = zeros(numAcq,1);
PE1Index = zeros(numAcq,1);
PE2Index = zeros(numAcq,1);
SliceIndex = zeros(numAcq,1);

row = 0;
for rep = 1:numRep
    for tri = 1:numTriangles
        for pol = [1,-1]
            for pe = 1:numPEAcq
                for sl = 1:numSlices
                    row = row + 1;
                    Repetition(row) = rep;
                    Triangle(row) = tri;
                    Polarity(row) = pol;
                    PE1Index(row) = validPhaseEncodes(pe,1);
                    PE2Index(row) = validPhaseEncodes(pe,2);
                    SliceIndex(row) = sl;
                end
            end
        end
    end
end

acq_manifest = table(Axis,Repetition,Triangle,Polarity, ...
    PE1Index,PE2Index,SliceIndex);

%% Encode the acquisition index into the raw data %%%%%%%%%%%%%%%%%%%%%%%%%
% Each acquisition has a unique value. After sorting, this allows every
% destination index to be checked exactly.
rawFID = zeros(numROP,numCoils,numAcq,'single');
for i=1:numAcq
    rawFID(:,:,i) = single(i);
end

[kspace,axisManifest] = sortHighOrderAcquisitions( ...
    rawFID,acq_manifest,axisName,nPE,numSlices);

assert(height(axisManifest) == numAcq, ...
    'Axis manifest has an unexpected number of acquisitions.');
assert(size(kspace,3) == nPE && size(kspace,4) == nPE, ...
    'Sorted PE matrix has the wrong dimensions.');
assert(size(kspace,5) == numSlices, ...
    'Sorted data has the wrong number of slices.');
assert(size(kspace,6) == 2*numTriangles, ...
    'Sorted polarity/triangle dimension has the wrong size.');
assert(size(kspace,9) == numRep, ...
    'Sorted repetition dimension has the wrong size.');

%% Verify every acquired location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numAcq
    pe1 = PE1Index(i);
    pe2 = PE2Index(i);
    sl = SliceIndex(i);
    tri = Triangle(i);
    rep = Repetition(i);

    if Polarity(i) > 0
        triPol = 2*tri - 1;
    else
        triPol = 2*tri;
    end

    value = kspace(1,1,pe1,pe2,sl,triPol,1,1,rep);
    assert(value == single(i), ...
        'Acquisition %d was placed at the wrong k-space index.',i);
end

%% Verify unacquired elliptical PE locations remain zero %%%%%%%%%%%%%%%%%
for pe1=1:nPE
    for pe2=1:nPE
        if ~peMask(pe1,pe2)
            value = kspace(1,1,pe1,pe2,1,1,1,1,1);
            assert(value == 0, ...
                'Unacquired elliptical PE location (%d,%d) is not zero.',pe1,pe2);
        end
    end
end

fprintf(['High-Order sorting regression passed: %d acquisitions, %d PE points, ' ...
    '%d slices, %d repetitions.\n'],numAcq,numPEAcq,numSlices,numRep);
