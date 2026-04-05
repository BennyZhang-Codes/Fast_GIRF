% =========================================================================
% Pulseq implementation for Trapezoid Validation at 7T
% =========================================================================
clc; clear; close all;
addpath(genpath('pulseq')); % 请根据实际情况修改路径

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


Setup.TR = 500e-3;
Setup.nRepetition = 8;
Setup.sliceThickness = 3e-3;
Setup.sliceOffset = 16.5e-3;
Setup.adcDwell = 9.8e-6;
Setup.adcDuration = 40.18e-3;
Setup.adcSamples = round(Setup.adcDuration / Setup.adcDwell);

Setup.RFflipEx   = 90;  
Setup.RFDuration = 2e-3;
Setup.RFTBP      = 4;

Setup.delayAfterEx  = 2e-3;
Setup.delayTestGrad = 1e-3;

Trapz.amp  = 5.4;    % mT/m
Trapz.rise = 40e-6;  % 40 us
Trapz.flat = 300e-6; % 300 us
Trapz.fall = 40e-6;  % 40 us

adc = mr.makeAdc(Setup.adcSamples, 'Dwell', Setup.adcDwell, 'system', sys);
[RF_base, Gz_base, GzReph_base] = mr.makeSincPulse(Setup.RFflipEx/180*pi, 'system', sys_soft, ...
    'Duration', Setup.RFDuration, 'SliceThickness', Setup.sliceThickness, ...
    'apodization', 0.5, 'timeBwProduct', Setup.RFTBP, 'use', 'excitation');

%%
Actual = Setup;
DelayAfterEx = mr.makeDelay(Setup.delayAfterEx);
DelayTRFill = mr.makeDelay(1);
axes_to_test = {'x', 'y', 'z'};

for ax_idx = 1:length(axes_to_test)
    ax = axes_to_test{ax_idx};
    
    GEx = Gz_base;       GEx.channel = ax;
    GExRe = GzReph_base; GExRe.channel = ax;
    
    rf_pos = RF_base;
    rf_pos.freqOffset = GEx.amplitude * Setup.sliceOffset;
    rf_pos.phaseOffset = -2 * pi * rf_pos.freqOffset * mr.calcRfCenter(rf_pos);
    
    rf_neg = RF_base;
    rf_neg.freqOffset = GEx.amplitude * (-Actual.sliceOffset);
    rf_neg.phaseOffset = -2 * pi * rf_neg.freqOffset * mr.calcRfCenter(rf_neg);
    rfs = {rf_pos, rf_neg};
    
    test_grad = mr.makeTrapezoid(ax, 'system', sys, ...
        'amplitude', Trapz.amp * 1e-3 * sys.gamma, ...
        'riseTime', Trapz.rise, 'flatTime', Trapz.flat, 'fallTime', Trapz.fall, 'delay', 0);
    
    for rep = 1:Actual.nRepetition
        for pol = [1, -1]
            for sl_idx = 1:2
                TimeInTR = 0;
                
                test_grad_pol = test_grad;
                test_grad_pol.amplitude = test_grad.amplitude * pol;
                test_grad_pol.delay = Actual.delayTestGrad;
                
                seq.addBlock(rfs{sl_idx}, GEx);   TimeInTR = TimeInTR + seq.blockDurations(end);
                seq.addBlock(GExRe);              TimeInTR = TimeInTR + seq.blockDurations(end);
                seq.addBlock(DelayAfterEx);       TimeInTR = TimeInTR + seq.blockDurations(end);
                seq.addBlock(test_grad_pol, adc); TimeInTR = TimeInTR + seq.blockDurations(end);
                
                TRFill = Actual.TR - TimeInTR;
                TRFill = round(TRFill / sys.gradRasterTime) * sys.gradRasterTime;
                DelayTRFill.delay = TRFill;
                seq.addBlock(DelayTRFill);
            end
        end
    end
end

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

seq.setDefinition('adcDwell'             , Setup.adcDwell              );
seq.setDefinition('adcDuration'          , Setup.adcDuration           );
seq.setDefinition('adcSamples'           , Setup.adcSamples            );

seq.setDefinition('delayAfterEx'         , Actual.delayAfterEx         );
seq.setDefinition('delayTestGrad'        , Actual.delayTestGrad        );

seq.setDefinition('nRepetition'          , Actual.nRepetition          );

seq.setDefinition('Developer'            , 'Jinyuan Zhang'             );
seq.setDefinition('Name'                 , 'gstf_7T_scholten'          );

seq.checkTiming();
seq.plot()

outpath = '';
seqname = sprintf('val_trapz_7T_tr%d_fa%s_rep%d', round(Actual.TR*1e3), num2str(Actual.RFflipEx), Actual.nRepetition);
seq.write([seqname, '.seq']);


%%
%% =========================================================================
% 生成 input_trapz.mat 供 Fast-GIRF 预测对比使用 (带精准延迟补偿)
% =========================================================================
disp('Generating nominal waveform...');

dt_1us = 1e-6;  
window = Actual.adcDuration;
t_nominal = 0:dt_1us:window;
grad_input = zeros(length(t_nominal), 1);

relative_delay_s = Actual.delayTestGrad - adc.delay; 
delay_pts = round(relative_delay_s / dt_1us); 


idx_up = (1:41) + delay_pts; 
grad_input(idx_up) = linspace(0, Trapz.amp, 41);


idx_flat = (42:341) + delay_pts; 
grad_input(idx_flat) = Trapz.amp;


idx_down = (342:382) + delay_pts; 
grad_input(idx_down) = linspace(Trapz.amp, 0, 41);


grad_input = single(grad_input);

lengthADC = round(Actual.adcSamples * Actual.adcDwell / dt_1us);

shift = 1; 

save('input_trapz.mat', 'grad_input', 'lengthADC', 'shift');
disp(['Saved input_trapz.mat successfully! Gradient starts at point: ', num2str(delay_pts + 1)]);