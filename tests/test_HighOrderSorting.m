% Synthetic regression test for High-Order acquisition sorting and manifest order.
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

%% Build a synthetic manifest in the current sequence order %%%%%%%%%%%%%%%
% triangle -> PE -> polarity -> slice
numAcq = numRep*numTriangles*2*numPEAcq*numSlices;
Axis = repmat({axisName},numAcq,1);
Repetition = zeros(numAcq,1);
Triangle = zeros(numAcq,1);
Scheme = zeros(numAcq,1);
Delay_s = zeros(numAcq,1);
Polarity = zeros(numAcq,1);
PE1Index = zeros(numAcq,1);
PE2Index = zeros(numAcq,1);
PE1Value = zeros(numAcq,1);
PE2Value = zeros(numAcq,1);
SliceIndex = zeros(numAcq,1);
SlicePosition_m = zeros(numAcq,1);

probeScheme = [1 1 2 2 2 2 2];
probeDelays = [0 0 1e-3 39e-3 77e-3 115e-3 153e-3];
sliceOffsets_m = [-34 -17 17 34]*1e-3;

row = 0;
for rep=1:numRep
    for tri=1:numTriangles
        for pe=1:numPEAcq
            for pol=[1,-1]
                for sl=1:numSlices
                    row = row+1;
                    Repetition(row) = rep;
                    Triangle(row) = tri;
                    Scheme(row) = probeScheme(tri);
                    Delay_s(row) = probeDelays(tri);
                    Polarity(row) = pol;
                    PE1Index(row) = validPhaseEncodes(pe,1);
                    PE2Index(row) = validPhaseEncodes(pe,2);
                    PE1Value(row) = peIndices(PE1Index(row));
                    PE2Value(row) = peIndices(PE2Index(row));
                    SliceIndex(row) = sl;
                    SlicePosition_m(row) = sliceOffsets_m(sl);
                end
            end
        end
    end
end

acq_manifest = table(Axis,Repetition,Triangle,Scheme,Delay_s,Polarity, ...
    PE1Index,PE2Index,PE1Value,PE2Value,SliceIndex,SlicePosition_m);

Protocol = struct();
Protocol.nPE = nPE;
Protocol.PEMode = 'elliptical';
Protocol.numPEAcq = numPEAcq;
Protocol.nRepetition = numRep;
Protocol.numTriangles = numTriangles;
Protocol.numSlices = numSlices;
Protocol.numAxes = 1;
Protocol.numPolarities = 2;
Protocol.probeScheme = probeScheme;
Protocol.probeDelays_s = probeDelays;
Protocol.sliceOffsets_m = sliceOffsets_m;
Protocol.expectedAcquisitions = numAcq;

% validateHighOrderManifest assumes the complete x/y/z acquisition. For a
% sorter-only single-axis test, make a temporary three-axis manifest below.
Axis3 = [repmat({'x'},numAcq,1);repmat({'y'},numAcq,1);repmat({'z'},numAcq,1)];
manifest3 = [acq_manifest;acq_manifest;acq_manifest];
manifest3.Axis = Axis3;
Protocol.numAxes = 3;
Protocol.expectedAcquisitions = 3*numAcq;
validateHighOrderManifest(manifest3,Protocol);

%% Encode acquisition index into raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawFID = zeros(numROP,numCoils,numAcq,'single');
for i=1:numAcq
    rawFID(:,:,i) = single(i);
end

[kspace,axisManifest] = sortHighOrderAcquisitions( ...
    rawFID,acq_manifest,axisName,nPE,numSlices);

assert(height(axisManifest)==numAcq,'Axis manifest has unexpected size.');
assert(size(kspace,3)==nPE && size(kspace,4)==nPE,'Sorted PE matrix has wrong dimensions.');
assert(size(kspace,5)==numSlices,'Sorted data has wrong number of slices.');
assert(size(kspace,6)==2*numTriangles,'Sorted polarity/triangle dimension has wrong size.');
assert(size(kspace,9)==numRep,'Sorted repetition dimension has wrong size.');

%% Verify every acquired location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numAcq
    pe1 = PE1Index(i);
    pe2 = PE2Index(i);
    sl = SliceIndex(i);
    tri = Triangle(i);
    rep = Repetition(i);
    if Polarity(i)>0
        triPol = 2*tri-1;
    else
        triPol = 2*tri;
    end
    value = kspace(1,1,pe1,pe2,sl,triPol,1,1,rep);
    assert(value==single(i),'Acquisition %d was placed at wrong k-space index.',i);
end

%% Verify unacquired elliptical PE locations remain zero %%%%%%%%%%%%%%%%%
for pe1=1:nPE
    for pe2=1:nPE
        if ~peMask(pe1,pe2)
            value = kspace(1,1,pe1,pe2,1,1,1,1,1);
            assert(value==0,'Unacquired PE location (%d,%d) is not zero.',pe1,pe2);
        end
    end
end

%% Legacy polarity-before-PE ordering must be rejected %%%%%%%%%%%%%%%%%%%%
legacy = acq_manifest;
legacyRow = 0;
for rep=1:numRep
    for tri=1:numTriangles
        for pol=[1,-1]
            for pe=1:numPEAcq
                for sl=1:numSlices
                    legacyRow = legacyRow+1;
                    legacy.Repetition(legacyRow) = rep;
                    legacy.Triangle(legacyRow) = tri;
                    legacy.Scheme(legacyRow) = probeScheme(tri);
                    legacy.Delay_s(legacyRow) = probeDelays(tri);
                    legacy.Polarity(legacyRow) = pol;
                    legacy.PE1Index(legacyRow) = validPhaseEncodes(pe,1);
                    legacy.PE2Index(legacyRow) = validPhaseEncodes(pe,2);
                    legacy.PE1Value(legacyRow) = peIndices(legacy.PE1Index(legacyRow));
                    legacy.PE2Value(legacyRow) = peIndices(legacy.PE2Index(legacyRow));
                    legacy.SliceIndex(legacyRow) = sl;
                    legacy.SlicePosition_m(legacyRow) = sliceOffsets_m(sl);
                end
            end
        end
    end
end
legacy3 = [legacy;legacy;legacy];
legacy3.Axis = Axis3;
legacyRejected = false;
try
    validateHighOrderManifest(legacy3,Protocol);
catch
    legacyRejected = true;
end
assert(legacyRejected,'Legacy polarity-before-PE manifest was not rejected.');

fprintf(['High-Order sorting regression passed: %d acquisitions/axis, %d PE points, ' ...
    '%d slices, %d repetitions.\n'],numAcq,numPEAcq,numSlices,numRep);
