% =========================================================================
% Pulseq implementation for Fast GSTF measurement at 7T
% Based on: Scholten et al., MRM 2023 (Delayed-excitation thin-slice)
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

Setup.TR = 2000e-3;            % Repetition time 500 ms
Setup.sliceThickness = 3e-3;  % Thin slice thickness 3 mm
Setup.sliceOffset = 16.5e-3;  % Offset excitation distance +/- 16.5 mm
Setup.adcDwell = 9.8e-6;      % Receiver bandwidth approx 200 kHz
Setup.adcDuration = 40.18e-3;    % ADC readout duration 40 ms
Setup.adcSamples = round(Setup.adcDuration / Setup.adcDwell); % Number of ADC samples

Setup.RFflipEx   = 90;
Setup.RFDuration = 2e-3;
Setup.RFTBP      = 4;

Setup.TriangleStartTime = 5.5e-3;

Setup.nRepetition = 8;
%% Parameters for the test gradients (calculated based on 180 T/m/s slew rate)
% Index 1-2: Short triangles for high-frequency probing (Scheme i)
% Index 3-7: Long triangles for low-frequency & long-term eddy current probing (Scheme ii)
tr_rise     = [50e-6, 80e-6, 140e-6, 310e-6, 310e-6, 310e-6, 310e-6];
tr_amp_mT   = [9.0,   14.4,  25.2,   55.8,   55.8,   55.8,   55.8];
tr_amp      = tr_amp_mT * 1e-3 * sys.gamma; % Convert to Pulseq standard unit Hz/m
delays      = [0,     0,     1e-3,   39e-3,  77e-3,  115e-3, 153e-3];
scheme_type = [1,     1,     2,      2,      2,      2,      2]; % 1: Excitation first (Scheme i), 2: Delayed excitation (Scheme ii)

%%
Actual = Setup;

adc = mr.makeAdc(Actual.adcSamples, 'Dwell', Actual.adcDwell, 'system', sys);
[RF.rf_base, Grad.gz_base, Grad.gzReph_base] = mr.makeSincPulse(Actual.RFflipEx/180*pi, 'system', sys_soft, ...
    'Duration', Actual.RFDuration, 'SliceThickness', Actual.sliceThickness, ...
    'apodization', 0.5, 'timeBwProduct', Actual ...
    .RFTBP, 'use', 'excitation');

DurationExcitation = mr.calcDuration(Grad.gz_base) + mr.calcDuration(Grad.gzReph_base);
Delay.delayTriangleScheme1 = round((Actual.TriangleStartTime - DurationExcitation) / sys.gradRasterTime) * sys.gradRasterTime;
Delay.delayTriangleScheme2 = round(Actual.TriangleStartTime / sys.gradRasterTime) * sys.gradRasterTime;

if (Delay.delayTriangleScheme1 < -eps(0)) % Sanity check
    error(['Total time (%f ms) of excitation blocks is ' ...
    'longer than desired TriangleStartTime (%f ms)!'], 1e3*DurationExcitation, 1e3*Actual.TriangleStartTime) ;
end

Delay.delayMeasShift = mr.makeDelay(1);
Delay.delayTRFill = mr.makeDelay(1) ;
MinTRActual = 0 ;
axes_to_test = {'x', 'y', 'z'};

TotalTRs = Actual.nRepetition * length(axes_to_test) * 7 * 2 * 2;
TRCounter = 0;

TriangleBlock = zeros(1, TotalTRs);
ReadoutBlock  = zeros(1, TotalTRs);
TRStartBlock  = zeros(1, TotalTRs);


for ax_idx = 1:length(axes_to_test)
    ax = axes_to_test{ax_idx};
    
    Grad.GEx = Grad.gz_base;
    Grad.GEx.channel = ax;
    
    Grad.GExRe = Grad.gzReph_base;
    Grad.GExRe.channel = ax;
    
    % Create two RF pulses for different offset positions
    RF.rf_pos = RF.rf_base;
    RF.rf_pos.freqOffset = Grad.GEx.amplitude * Actual.sliceOffset;
    RF.rf_pos.phaseOffset = -2 * pi * RF.rf_pos.freqOffset * mr.calcRfCenter(RF.rf_pos);
    
    RF.rf_neg = RF.rf_base;
    RF.rf_neg.freqOffset = Grad.GEx.amplitude * (-Actual.sliceOffset);
    RF.rf_neg.phaseOffset = -2 * pi * RF.rf_neg.freqOffset * mr.calcRfCenter(RF.rf_neg);
    RF.rfs = {RF.rf_pos, RF.rf_neg};
    
    for irep = 1:Actual.nRepetition
        % Start sequence inner loop (7 triangles x 2 polarities x 2 offset positions)
        for t_idx = 1:7
            for pol = [1, -1]
                for sl_idx = 1:2
                    TRCounter = TRCounter + 1;
    
                    TRStartBlock(TRCounter) = length(seq.blockEvents)+1;
    
                    TimeInTR = 0 ; % [s]
                    
                    % Generate current probing triangular gradient (flat top time = 0)
                    test_grad = mr.makeTrapezoid(ax, 'system', sys, ...
                        'amplitude', tr_amp(t_idx) * pol, ...
                        'riseTime', tr_rise(t_idx), 'fallTime', tr_rise(t_idx), 'flatTime', 0, 'delay', 0);
    
                    Delay.delayMeasShift.delay = round((delays(t_idx)-2*tr_rise(t_idx)) / sys.gradRasterTime) * sys.gradRasterTime;
    
                    if scheme_type(t_idx) == 1
                        test_grad.delay = Delay.delayTriangleScheme1;
                        seq.addBlock(RF.rfs{sl_idx}, Grad.GEx); TimeInTR = TimeInTR + seq.blockDurations(end);
                        seq.addBlock(Grad.GExRe);               TimeInTR = TimeInTR + seq.blockDurations(end);
                        seq.addBlock(test_grad, adc);           TimeInTR = TimeInTR + seq.blockDurations(end);
                        TriangleBlock(TRCounter) = length(seq.blockEvents);
                    else
                        test_grad.delay = Delay.delayTriangleScheme2;
                        seq.addBlock(test_grad);                TimeInTR = TimeInTR + seq.blockDurations(end);
                        TriangleBlock(TRCounter) = length(seq.blockEvents);
                        seq.addBlock(Delay.delayMeasShift);     TimeInTR = TimeInTR + seq.blockDurations(end);
                        seq.addBlock(RF.rfs{sl_idx}, Grad.GEx); TimeInTR = TimeInTR + seq.blockDurations(end);
                        seq.addBlock(Grad.GExRe);               TimeInTR = TimeInTR + seq.blockDurations(end);
                        seq.addBlock(adc);                      TimeInTR = TimeInTR + seq.blockDurations(end);
                    end
                    ReadoutBlock(TRCounter) = length(seq.blockEvents);
    
                    MinTRActual = max(MinTRActual, TimeInTR) ;
    
                    TRFill = Actual.TR - TimeInTR ;
                    TRFill = round(TRFill / sys.gradRasterTime) * sys.gradRasterTime;
    
                    % Sanity check
                    if (TRFill < -eps(0))
                        error(['Total time (%f ms) of blocks within current TR (#%d) is ' ...
                        'longer than desired TR (%f ms)!'], 1e3*TimeInTR, TRCounter, 1e3*Actual.TR) ;
                    end
    
                    Delay.delayTRFill.delay = TRFill ; % update delay of eTRFill
                    seq.addBlock(Delay.delayTRFill)  ;  % Add delay to the sequence
                    TimeInTR = TimeInTR + seq.blockDurations(end) ;
    
                end
            end
        end
    end
end
outpath = '';

%%
waveform = seq.waveforms_and_times();
BlockEndTimes   = cumsum(seq.blockDurations);
BlockStartTimes = BlockEndTimes - seq.blockDurations;

window = 210e-3;
dt_1us = 1e-6;
dt     = sys.gradRasterTime; % [s]
grad_nominal = cell(1, 7);
TimeShiftADC = zeros(1, 7);

idxs = 1:4:28;
for itri = 1:7
    idx = idxs(itri);

    tADCStart = BlockStartTimes(ReadoutBlock(idx)) - BlockStartTimes(TRStartBlock(idx)) + adc.delay;
    TimeShiftADC(itri) = round(tADCStart / sys.gradRasterTime) * sys.gradRasterTime;
    % disp(t_ADC - t_TR)

    % t_e = BlockEndTimes(TRStartBlock(idx+1)-1);
    t_s = BlockStartTimes(TRStartBlock(idx));
    t_e = t_s + window;
    t_win = t_s:dt_1us:t_e; % time of the TR

    t_e = BlockEndTimes(TriangleBlock(idx));
    t_s = BlockStartTimes(TriangleBlock(idx));
    t_grad = t_s:dt:t_e; % time of the readout gradient

    grad_shape = zeros(length(t_win), 1);
    for i = 1:1
        if size(waveform{i},2) > 0
            g = interp1(waveform{i}(1,:), waveform{i}(2,:), t_grad, 'linear', 0);
            grad_shape(:,i) = interp1(t_grad, g, t_win, 'linear', 0);
        end
    end
    grad_shape = mr.convert(grad_shape, 'Hz/m', 'mT/m', 'gamma', sys.gamma);
    grad_nominal{itri} = reshape(grad_shape, [1, size(grad_shape)]);
end

t_nominal = 0:dt_1us:window;
grad_nominal = cell2mat(grad_nominal')';

figure; hold on;
for i = 1:7
    plot(t_nominal*1e3, grad_nominal(:, i))
end

%
grad_input = single(grad_nominal);
lengthADC = round(adc.numSamples * adc.dwell / dt_1us);
shift = round(TimeShiftADC / dt_1us) + 1;

save(strcat(outpath, 'input_H_fast.mat'), 'grad_input', 'lengthADC', 'shift');

%%
[seq] = check_Timing(seq);
Actual.ScannerType = 'Terra-XJ';
[seq] = check_PNS(seq, Actual);

seq.setDefinition('FOV'                  , [200, 200, 200] * 1e-3      );
seq.setDefinition('BaseResolution'       , Actual.adcSamples           ); 
seq.setDefinition('MatrixSize'           , [Actual.adcSamples, 1, 1]   );

seq.setDefinition('TR'                   , Actual.TR                   );
seq.setDefinition('TE'                   , 5e-3                        );
seq.setDefinition('Excit_FlipAngle'      , Actual.RFflipEx             );
seq.setDefinition('Excit_Duration'       , Actual.RFDuration           );
seq.setDefinition('Excit_TBP'            , Actual.RFTBP                );
seq.setDefinition('Excit_TBP'            , Actual.RFTBP                );

seq.setDefinition('GSTF_rise'            , tr_rise                     );
seq.setDefinition('GSTF_amp_mT'          , tr_amp_mT                   );
seq.setDefinition('GSTF_delays'          , delays                      );
seq.setDefinition('GSTF_scheme_type'     , scheme_type                 );

seq.setDefinition('adcDwell'             , Setup.adcDwell              );
seq.setDefinition('adcDuration'          , Setup.adcDuration           );
seq.setDefinition('adcSamples'           , Setup.adcSamples            );

seq.setDefinition('TriangleStartTime'    , Actual.TriangleStartTime    );
seq.setDefinition('TimeShiftADC_us'      , shift                       );


seq.setDefinition('nRepetition'          , Actual.nRepetition          );

seq.setDefinition('Developer'            , 'Jinyuan Zhang'             );
seq.setDefinition('Name'                 , 'gstf_7T_scholten'          );

seq.checkTiming();
seq.plot()


seqname = sprintf('fast_gstf_7T_tr%s_fa%s_rep%d', num2str(Actual.TR*1e3), num2str(Actual.RFflipEx), Actual.nRepetition);

seq.write(strcat(outpath, seqname,'.seq'));

% %%
% rep = seq.testReport;
% fprintf([rep{:}]);
