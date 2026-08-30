Fast High-Order GIRF: pre-scan validation
========================================

1. Generate the High-Order sequence

   sequence/write_GSTF_Scholten2023_HighOrder.m

2. Validate the written .seq file directly

   report = check_HighOrderSequence( ...
       '../sequence/fast_gstf_HighOrder_7T_pe5_elliptical_tr500_fa90_rep1.seq');

   The checker verifies:
     - total TR and ADC count
     - total duration
     - x/y/z ordering
     - four physical slice positions
     - 2D PE pattern
     - +/- probing-gradient polarity
     - probing-gradient rise time/amplitude consistency
     - seven ADC start times

3. Keep these three files together for the scan and processing

   fast_gstf_HighOrder_*.seq
   input_H_fast_HighOrder.mat
   acq_manifest_HighOrder.mat

   Do not regenerate the manifest after the scanner acquisition. The
   manifest must correspond to the exact .seq file used on the scanner.

4. After acquisition, prepare the Twix data using the exact .seq file

   [outputFiles, Protocol] = prepareHighOrderTwixFromSequence( ...
       'meas_MIDxxxx.dat', ...
       '../sequence/acq_manifest_HighOrder.mat', ...
       '../sequence/fast_gstf_HighOrder_7T_pe5_elliptical_tr500_fa90_rep1.seq', ...
       '../triangle_measurements/');

   This reads ADC dwell time, nPE, PE FOV and slice positions directly from
   the .seq Definitions. No manual duplication of these protocol values is
   required.

5. Process one physical gradient axis

   GSTF_calculation/main_H_fast_HighOrder.m

   Set meas_name / ax_name to x, y or z and choose calcChannels = 4, 9 or 16.
