% =========================================================================
% Pulseq implementation for Fast GSTF measurement at 7T
% Based on: Scholten et al., MRM 2023 (Delayed-excitation thin-slice)
% =========================================================================
clc; clear; close all;

addpath(genpath(''));
addpath(genpath('check' ));
addpath(genpath('./pulseq'));


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
Setup.adcDwell = 5.0e-6;      % Receiver bandwidth approx 200 kHz
Setup.adcDuration = 40e-3;    % ADC readout duration 40 ms
Setup.adcSamples = round(Setup.adcDuration / Setup.adcDwell); % Number of ADC samples

Setup.RFflipEx   = 90;
Setup.RFDuration = 2e-3;
Setup.RFTBP      = 3;

Setup.DelayAfterEx = 2e-3;

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

Delay.delayTRFill = mr.makeDelay(1) ;
MinTRActual = 0 ;
axes_to_test = {'x', 'y', 'z'};
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
    
    % Start sequence inner loop (7 triangles x 2 polarities x 2 offset positions)
    for t_idx = 1:7
        for pol = [1, -1]
            for sl_idx = 1:2
                TimeInTR = 0 ; % [s]
                
                % Generate current probing triangular gradient (flat top time = 0)
                test_grad = mr.makeTrapezoid(ax, 'system', sys, ...
                    'amplitude', tr_amp(t_idx) * pol, ...
                    'riseTime', tr_rise(t_idx), 'fallTime', tr_rise(t_idx), 'flatTime', 0, 'delay', Actual.DelayAfterEx);

                delay_blk = mr.makeDelay(round(delays(t_idx) / sys.gradRasterTime) * sys.gradRasterTime);

                if scheme_type(t_idx) == 1
                    seq.addBlock(RF.rfs{sl_idx}, Grad.GEx); TimeInTR = TimeInTR + seq.blockDurations(end);
                    seq.addBlock(Grad.GExRe);               TimeInTR = TimeInTR + seq.blockDurations(end);
                    seq.addBlock(test_grad, adc);           TimeInTR = TimeInTR + seq.blockDurations(end);
                else
                    seq.addBlock(test_grad);                TimeInTR = TimeInTR + seq.blockDurations(end);
                    seq.addBlock(delay_blk);                TimeInTR = TimeInTR + seq.blockDurations(end);
                    seq.addBlock(RF.rfs{sl_idx}, Grad.GEx); TimeInTR = TimeInTR + seq.blockDurations(end);
                    seq.addBlock(Grad.GExRe);               TimeInTR = TimeInTR + seq.blockDurations(end);
                    seq.addBlock(adc);                      TimeInTR = TimeInTR + seq.blockDurations(end);
                end

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
% seq.plot()
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

seq.setDefinition('Developer'            , 'Jinyuan Zhang'             );
seq.setDefinition('Name'                 , 'gstf_7T_scholten'          );

seq.checkTiming();
seq.plot()

seq.write('fast_gstf_7T_scholten.seq');

% %%
% rep = seq.testReport;
% fprintf([rep{:}]);
