% =========================================================================
% Pulseq implementation for Fast High-Order GSTF measurement at 7T
% Based on: Scholten et al., MRM 2023 (delayed-excitation thin-slice)
% Spatial encoding based on the phase-encoded thin-slice approach for
% higher-order field monitoring.
% =========================================================================
clc; clear; close all;
%%
addpath(genpath(''));
addpath(genpath('check' ));
addpath(genpath('../pulseq'));


sys = mr.opts('MaxGrad',70, 'GradUnit', 'mT/m', ...
    'MaxSlew', 200, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 100e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6, 'adcRasterTime',100e-9,...
    'rfRasterTime',1e-6,'gradRasterTime',10e-6, ...
    'B0', 6.98, 'adcSamplesDivisor', 4);

sys_soft = mr.opts('MaxGrad',40, 'GradUnit', 'mT/m', ...
    'MaxSlew', 50, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 100e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6, 'adcRasterTime',100e-9,...
    'rfRasterTime',1e-6,'gradRasterTime',10e-6, ...
    'B0', 6.98, 'adcSamplesDivisor', 4);

seq = mr.Sequence(sys);

RF    = struct();
Grad  = struct();
ADC   = struct();
Delay = struct();
Label = struct();
%% ==================== Core Parameter Settings ====================

Setup.ScannerType      = 'Terra-XR';  % for pns_check

% A shorter TR and one repetition are used as the default development
% protocol because the addition of 2D phase encoding substantially
% increases the number of TRs. These values can be changed as needed.
Setup.TR = 500e-3;
Setup.nRepetition = 1;

Setup.sliceThickness = 3e-3;                       % Thin slice thickness 3 mm
Setup.sliceOffsets = [-34, -17, 17, 34] * 1e-3;  % Four offset slices

Setup.PEFOV = 200e-3;             % Phase-encoding FOV [m]
Setup.nPE = 5;                    % Odd number of PE steps per direction
Setup.PEMode = 'elliptical';      % 'elliptical' or 'full'

Setup.adcDwell = 9.8e-6;          % Receiver dwell time
Setup.adcDuration = 40.18e-3;     % ADC readout duration
Setup.adcSamples = round(Setup.adcDuration / Setup.adcDwell);

Setup.RFflipEx   = 90;
Setup.RFDuration = 2e-3;
Setup.RFTBP      = 4;

Setup.TriangleStartTime = 5.5e-3;

%% Parameters for the test gradients (calculated based on 180 T/m/s slew rate)
% Index 1-2: Short triangles for high-frequency probing (Scheme i)
% Index 3-7: Long triangles for low-frequency & long-term eddy current probing (Scheme ii)
tr_rise     = [50e-6, 80e-6, 140e-6, 310e-6, 310e-6, 310e-6, 310e-6];
tr_amp_mT   = [9.0,   14.4,  25.2,   55.8,   55.8,   55.8,   55.8];
tr_amp      = tr_amp_mT * 1e-3 * sys.gamma; % Convert to Pulseq standard unit Hz/m
delays      = [0,     0,     1e-3,   39e-3,  77e-3,  115e-3, 153e-3];
scheme_type = [1,     1,     2,      2,      2,      2,      2]; % 1: excitation first, 2: delayed excitation

%% Phase-encoding pattern %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~mod(Setup.nPE,2)
    error('Setup.nPE must be odd.');
end

peIndices = -(Setup.nPE-1)/2:(Setup.nPE-1)/2;
peAreas = peIndices / Setup.PEFOV; % Pulseq gradient area [1/m]
[PE1Grid,PE2Grid] = ndgrid(1:Setup.nPE,1:Setup.nPE);

if strcmpi(Setup.PEMode,'elliptical')
    radius = (Setup.nPE-1)/2;
    if radius == 0
        peMask = true(Setup.nPE,Setup.nPE);
    else
        peMask = (peIndices(PE1Grid)/radius).^2 + ...
                 (peIndices(PE2Grid)/radius).^2 <= 1 + eps;
    end
elseif strcmpi(Setup.PEMode,'full')
    peMask = true(Setup.nPE,Setup.nPE);
else
    error('Unknown Setup.PEMode. Use ''elliptical'' or ''full''.');
end

[validPE1,validPE2] = find(peMask);
validPhaseEncodes = [validPE1, validPE2];
numPEAcq = size(validPhaseEncodes,1);

fprintf('High-Order GIRF spatial encoding: %d x %d %s PE, %d acquired PE points.\n', ...
    Setup.nPE,Setup.nPE,Setup.PEMode,numPEAcq);

%%
Actual = Setup;

adc = mr.makeAdc(Actual.adcSamples, 'Dwell', Actual.adcDwell, 'system', sys);
[RF.rf_base, Grad.gz_base, Grad.gzReph_base] = mr.makeSincPulse(Actual.RFflipEx/180*pi, 'system', sys_soft, ...
    'Duration', Actual.RFDuration, 'SliceThickness', Actual.sliceThickness, ...
    'apodization', 0.5, 'timeBwProduct', Actual.RFTBP, 'use', 'excitation');

% Ensure that the slice-refocusing block is long enough to contain the
% largest phase-encoding gradient without changing the excitation timing
% more than necessary.
maxPEArea = max(abs(peAreas));
if maxPEArea > 0
    peTemp = mr.makeTrapezoid('x','system',sys_soft,'area',maxPEArea);
    peDuration = max(mr.calcDuration(Grad.gzReph_base),mr.calcDuration(peTemp));
else
    peDuration = mr.calcDuration(Grad.gzReph_base);
end

if peDuration > mr.calcDuration(Grad.gzReph_base) + eps
    gzRephArea = Grad.gzReph_base.area;
    Grad.gzReph_base = mr.makeTrapezoid('z','system',sys_soft, ...
        'area',gzRephArea,'duration',peDuration);
end

DurationExcitation = mr.calcDuration(Grad.gz_base) + mr.calcDuration(Grad.gzReph_base);
Delay.delayTriangleScheme1 = round((Actual.TriangleStartTime - DurationExcitation) / sys.gradRasterTime) * sys.gradRasterTime;
Delay.delayTriangleScheme2 = round(Actual.TriangleStartTime / sys.gradRasterTime) * sys.gradRasterTime;

if (Delay.delayTriangleScheme1 < -eps(0))
    error(['Total time (%f ms) of excitation/spatial-encoding blocks is ' ...
        'longer than desired TriangleStartTime (%f ms)!'], ...
        1e3*DurationExcitation,1e3*Actual.TriangleStartTime);
end

Delay.delayMeasShift = mr.makeDelay(1);
Delay.delayTRFill = mr.makeDelay(1);
MinTRActual = 0;
axes_to_test = {'x', 'y', 'z'};

numSlices = length(Actual.sliceOffsets);
TotalTRs = Actual.nRepetition * length(axes_to_test) * 7 * 2 * numPEAcq * numSlices;
TRCounter = 0;

TriangleBlock = zeros(1,TotalTRs);
ReadoutBlock  = zeros(1,TotalTRs);
TRStartBlock  = zeros(1,TotalTRs);
NominalTRIndex = zeros(1,7);

% Acquisition manifest. This does not change the scanner sequence; it is
% saved alongside the sequence so that raw-data sorting can reproduce the
% exact acquisition ordering.
Manifest = repmat(struct('Axis','','Repetition',0,'Triangle',0,'Scheme',0, ...
    'Delay_s',0,'Polarity',0,'PE1Index',0,'PE2Index',0,'PE1Value',0, ...
    'PE2Value',0,'SliceIndex',0,'SlicePosition_m',0),TotalTRs,1);

for ax_idx = 1:length(axes_to_test)
    ax = axes_to_test{ax_idx};

    Grad.GEx = Grad.gz_base;
    Grad.GEx.channel = ax;

    Grad.GExRe = Grad.gzReph_base;
    Grad.GExRe.channel = ax;

    % Match the two PE dimensions to createPositionArray.m exactly:
    %   dSag (slice x): PE1 = y, PE2 = z
    %   dCor (slice y): PE1 = x, PE2 = z
    %   dTra (slice z): PE1 = y, PE2 = x
    % Keeping this convention avoids any coordinate swap in the existing
    % Fast_GIRF spatial reconstruction code.
    if strcmp(ax,'x')
        pe_channels = {'y','z'};
    elseif strcmp(ax,'y')
        pe_channels = {'x','z'};
    else
        pe_channels = {'y','x'};
    end

    % Create RF pulses for all offset slice positions.
    RF.rfs = cell(1,numSlices);
    for sl_idx = 1:numSlices
        RF.rfs{sl_idx} = RF.rf_base;
        RF.rfs{sl_idx}.freqOffset = Grad.GEx.amplitude * Actual.sliceOffsets(sl_idx);
        RF.rfs{sl_idx}.phaseOffset = -2*pi*RF.rfs{sl_idx}.freqOffset*mr.calcRfCenter(RF.rfs{sl_idx});
    end

    for irep = 1:Actual.nRepetition
        % 7 triangles x 2 polarities x spatial PE points x offset slices
        for t_idx = 1:7
            for pol = [1, -1]
                for pe_idx = 1:numPEAcq
                    pe1_idx = validPhaseEncodes(pe_idx,1);
                    pe2_idx = validPhaseEncodes(pe_idx,2);
                    peArea1 = peAreas(pe1_idx);
                    peArea2 = peAreas(pe2_idx);

                    for sl_idx = 1:numSlices
                        TRCounter = TRCounter + 1;
                        TRStartBlock(TRCounter) = length(seq.blockEvents)+1;
                        TimeInTR = 0;

                        test_grad = mr.makeTrapezoid(ax,'system',sys, ...
                            'amplitude',tr_amp(t_idx)*pol, ...
                            'riseTime',tr_rise(t_idx),'fallTime',tr_rise(t_idx), ...
                            'flatTime',0,'delay',0);

                        Delay.delayMeasShift.delay = round((delays(t_idx)-2*tr_rise(t_idx)) / sys.gradRasterTime) * sys.gradRasterTime;

                        if scheme_type(t_idx) == 1
                            test_grad.delay = Delay.delayTriangleScheme1;
                            seq.addBlock(RF.rfs{sl_idx},Grad.GEx); TimeInTR = TimeInTR + seq.blockDurations(end);
                            addSpatialEncodingBlock(seq,Grad.GExRe,pe_channels,peArea1,peArea2,peDuration,sys_soft);
                            TimeInTR = TimeInTR + seq.blockDurations(end);
                            seq.addBlock(test_grad,adc); TimeInTR = TimeInTR + seq.blockDurations(end);
                            TriangleBlock(TRCounter) = length(seq.blockEvents);
                        else
                            test_grad.delay = Delay.delayTriangleScheme2;
                            seq.addBlock(test_grad); TimeInTR = TimeInTR + seq.blockDurations(end);
                            TriangleBlock(TRCounter) = length(seq.blockEvents);
                            seq.addBlock(Delay.delayMeasShift); TimeInTR = TimeInTR + seq.blockDurations(end);
                            seq.addBlock(RF.rfs{sl_idx},Grad.GEx); TimeInTR = TimeInTR + seq.blockDurations(end);
                            addSpatialEncodingBlock(seq,Grad.GExRe,pe_channels,peArea1,peArea2,peDuration,sys_soft);
                            TimeInTR = TimeInTR + seq.blockDurations(end);
                            seq.addBlock(adc); TimeInTR = TimeInTR + seq.blockDurations(end);
                        end
                        ReadoutBlock(TRCounter) = length(seq.blockEvents);

                        if ax_idx == 1 && irep == 1 && pol == 1 && pe_idx == 1 && sl_idx == 1
                            NominalTRIndex(t_idx) = TRCounter;
                        end

                        Manifest(TRCounter).Axis = ax;
                        Manifest(TRCounter).Repetition = irep;
                        Manifest(TRCounter).Triangle = t_idx;
                        Manifest(TRCounter).Scheme = scheme_type(t_idx);
                        Manifest(TRCounter).Delay_s = delays(t_idx);
                        Manifest(TRCounter).Polarity = pol;
                        Manifest(TRCounter).PE1Index = pe1_idx;
                        Manifest(TRCounter).PE2Index = pe2_idx;
                        Manifest(TRCounter).PE1Value = peIndices(pe1_idx);
                        Manifest(TRCounter).PE2Value = peIndices(pe2_idx);
                        Manifest(TRCounter).SliceIndex = sl_idx;
                        Manifest(TRCounter).SlicePosition_m = Actual.sliceOffsets(sl_idx);

                        MinTRActual = max(MinTRActual,TimeInTR);
                        TRFill = Actual.TR - TimeInTR;
                        TRFill = round(TRFill / sys.gradRasterTime) * sys.gradRasterTime;

                        if (TRFill < -eps(0))
                            error(['Total time (%f ms) of blocks within current TR (#%d) is ' ...
                                'longer than desired TR (%f ms)!'], ...
                                1e3*TimeInTR,TRCounter,1e3*Actual.TR);
                        end

                        Delay.delayTRFill.delay = TRFill;
                        seq.addBlock(Delay.delayTRFill);
                        TimeInTR = TimeInTR + seq.blockDurations(end);
                    end
                end
            end
        end
    end
end
outpath = '';

%% Generate nominal input gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waveform = seq.waveforms_and_times();
BlockEndTimes   = cumsum(seq.blockDurations);
BlockStartTimes = BlockEndTimes - seq.blockDurations;

window = 210e-3;
dt_1us = 1e-6;
dt     = sys.gradRasterTime;
grad_nominal = cell(1,7);
TimeShiftADC = zeros(1,7);

for itri = 1:7
    idx = NominalTRIndex(itri);
    if idx == 0
        error('Could not determine nominal TR index for triangle %d.',itri);
    end

    tADCStart = BlockStartTimes(ReadoutBlock(idx)) - BlockStartTimes(TRStartBlock(idx)) + adc.delay;
    TimeShiftADC(itri) = round(tADCStart / sys.gradRasterTime) * sys.gradRasterTime;

    t_s = BlockStartTimes(TRStartBlock(idx));
    t_e = t_s + window;
    t_win = t_s:dt_1us:t_e;

    t_e = BlockEndTimes(TriangleBlock(idx));
    t_s = BlockStartTimes(TriangleBlock(idx));
    t_grad = t_s:dt:t_e;

    grad_shape = zeros(length(t_win),1);
    if size(waveform{1},2) > 0
        g = interp1(waveform{1}(1,:),waveform{1}(2,:),t_grad,'linear',0);
        grad_shape(:,1) = interp1(t_grad,g,t_win,'linear',0);
    end
    grad_shape = mr.convert(grad_shape,'Hz/m','mT/m','gamma',sys.gamma);
    grad_nominal{itri} = reshape(grad_shape,[1,size(grad_shape)]);
end

t_nominal = 0:dt_1us:window;
grad_nominal = cell2mat(grad_nominal')';

figure; hold on;
for i = 1:7
    plot(t_nominal*1e3,grad_nominal(:,i))
end
xlabel('Time (ms)');
ylabel('Gradient (mT/m)');
title('Nominal Fast-GIRF probing gradients');

grad_input = single(grad_nominal);
lengthADC = round(adc.numSamples * adc.dwell / dt_1us);
shift = round(TimeShiftADC / dt_1us) + 1;

save(strcat(outpath,'input_H_fast_HighOrder.mat'),'grad_input','lengthADC','shift');

%% Save acquisition manifest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acq_manifest = struct2table(Manifest);
save(strcat(outpath,'acq_manifest_HighOrder.mat'),'acq_manifest','validPhaseEncodes','peMask','peAreas','peIndices');
writetable(acq_manifest,strcat(outpath,'acq_manifest_HighOrder.csv'));

%% Sequence checks and definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[seq] = check_Timing(seq);
Actual.ScannerType = 'Terra-XJ';
[seq] = check_PNS(seq,Actual);

seq.setDefinition('FOV'                  , [Actual.PEFOV, Actual.PEFOV, Actual.PEFOV]);
seq.setDefinition('BaseResolution'       , Actual.adcSamples);
seq.setDefinition('MatrixSize'           , [Actual.adcSamples, Actual.nPE, Actual.nPE]);

seq.setDefinition('TR'                   , Actual.TR);
seq.setDefinition('TE'                   , 5e-3);
seq.setDefinition('Excit_FlipAngle'      , Actual.RFflipEx);
seq.setDefinition('Excit_Duration'       , Actual.RFDuration);
seq.setDefinition('Excit_TBP'            , Actual.RFTBP);

seq.setDefinition('GSTF_rise'            , tr_rise);
seq.setDefinition('GSTF_amp_mT'          , tr_amp_mT);
seq.setDefinition('GSTF_delays'          , delays);
seq.setDefinition('GSTF_scheme_type'     , scheme_type);

seq.setDefinition('adcDwell'             , Actual.adcDwell);
seq.setDefinition('adcDuration'          , Actual.adcDuration);
seq.setDefinition('adcSamples'           , Actual.adcSamples);

seq.setDefinition('TriangleStartTime'    , Actual.TriangleStartTime);
seq.setDefinition('TimeShiftADC_us'      , shift);
seq.setDefinition('nRepetition'          , Actual.nRepetition);

seq.setDefinition('HighOrder_nPE'        , Actual.nPE);
seq.setDefinition('HighOrder_PEFOV'      , Actual.PEFOV);
seq.setDefinition('HighOrder_PEMode'     , Actual.PEMode);
seq.setDefinition('HighOrder_numPEAcq'   , numPEAcq);
seq.setDefinition('HighOrder_sliceOffsets_m', Actual.sliceOffsets);

seq.setDefinition('Developer'            , 'Jinyuan Zhang');
seq.setDefinition('Name'                 , 'gstf_7T_scholten_HighOrder');

seq.checkTiming();
seq.plot()

seqname = sprintf('fast_gstf_HighOrder_7T_pe%d_%s_tr%s_fa%s_rep%d', ...
    Actual.nPE,Actual.PEMode,num2str(Actual.TR*1e3), ...
    num2str(Actual.RFflipEx),Actual.nRepetition);
seq.write(strcat(outpath,seqname,'.seq'));

fprintf('Minimum required TR: %.2f ms\n',1e3*MinTRActual);
fprintf('Total number of TRs: %d\n',TotalTRs);
fprintf('Estimated scan time: %.2f min\n',TotalTRs*Actual.TR/60);

% %%
% rep = seq.testReport;
% fprintf([rep{:}]);

%% Local function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addSpatialEncodingBlock(seq,GExRe,pe_channels,peArea1,peArea2,peDuration,sys_soft)
% Add slice refocusing and the two orthogonal phase-encoding gradients in
% one block. Zero-area PE events are omitted.

events = {GExRe};

if abs(peArea1) > eps
    GPE1 = mr.makeTrapezoid(pe_channels{1},'system',sys_soft, ...
        'area',peArea1,'duration',peDuration);
    events{end+1} = GPE1;
end

if abs(peArea2) > eps
    GPE2 = mr.makeTrapezoid(pe_channels{2},'system',sys_soft, ...
        'area',peArea2,'duration',peDuration);
    events{end+1} = GPE2;
end

seq.addBlock(events{:});
end
