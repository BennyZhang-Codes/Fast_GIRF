# High-Order Fast_GIRF development

This document describes the current development workflow on the `dev_HighOrder` branch. The original Fast_GIRF scripts remain available for first-order/two-slice measurements.

## Goal

Extend the Scholten et al. fast GSTF measurement with spatial encoding for field expansions up to third order while retaining the original delayed-excitation and LF/HF GSTF calculation framework.

The spatial basis remains the existing Fast_GIRF/Skope-compatible 16-term real solid-harmonic basis implemented in `GSTF_calculation/createProbingMatrix.m`.

## Acquisition

Generate the High-Order Pulseq sequence with:

```matlab
sequence/write_GSTF_Scholten2023_HighOrder.m
```

The development defaults are:

- four offset slices: `[-34 -17 17 34] mm`
- `5 x 5` PE matrix
- elliptical PE mask (13 acquired PE points)
- TR = 500 ms
- one repetition
- the original seven Fast_GIRF probing gradients and delayed-excitation timings
- positive and negative gradient polarity
- x, y, and z gradient axes in one sequence

The PE axes are matched to the existing `createPositionArray.m` convention:

- x / sagittal: PE1 = y, PE2 = z
- y / coronal: PE1 = x, PE2 = z
- z / transverse: PE1 = y, PE2 = x

The sequence also writes:

- `input_H_fast_HighOrder.mat`
- `acq_manifest_HighOrder.mat`
- `acq_manifest_HighOrder.csv`

The manifest records axis, repetition, test gradient, scheme, delay, polarity, PE indices, and slice position for every ADC.

## Siemens Twix extraction

The High-Order processing uses the acquisition manifest as the source of truth for spatial and gradient ordering. Siemens image ADCs should therefore first be read in acquisition order.

`mapVBVD` must be on the MATLAB path. Then use:

```matlab
[rawFID, twix] = readHighOrderTwix( ...
    'meas_MIDxxxx.dat', ...
    'acq_manifest_HighOrder.mat');
```

`readHighOrderTwix.m` uses:

```matlab
twix.image.unsorted()
```

and returns:

```matlab
rawFID  % [numROP, coils, acquisitions]
```

Readout oversampling and averaging are deliberately disabled so that the acquired ADC stream is kept unchanged before manifest-based sorting.

If the Twix image stream contains extra image ADCs, specify the first High-Order GIRF acquisition explicitly:

```matlab
rawFID = readHighOrderTwix(datFile, manifestFile, firstAcquisition);
```

## Raw-data sorting

Sequential ADCs are sorted into the existing Fast_GIRF k-space layout with:

```matlab
[kspace_x, manifest_x] = sortHighOrderAcquisitions( ...
    rawFID, 'acq_manifest_HighOrder.mat', 'x', 5, 4);
```

The sorter fills the complete PE matrix and leaves unacquired elliptical PE positions at zero.

Save the Fast_GIRF measurement file with:

```matlab
saveHighOrderMeasurement( ...
    'measurement_H_fast_HighOrder_x.mat', ...
    kspace_x, 9.8e-6, 200, 'x', [-34 -17 17 34]);
```

For all three physical axes at once, use:

```matlab
outputFiles = prepareHighOrderMeasurements( ...
    rawFID, ...
    'acq_manifest_HighOrder.mat', ...
    '../triangle_measurements/', ...
    9.8e-6, ...
    200, ...
    [-34 -17 17 34], ...
    5);
```

The complete Siemens `.dat` to Fast_GIRF measurement-file workflow can be run with:

```matlab
outputFiles = prepareHighOrderTwix( ...
    'meas_MIDxxxx.dat', ...
    'acq_manifest_HighOrder.mat', ...
    '../triangle_measurements/', ...
    9.8e-6, ...
    200, ...
    [-34 -17 17 34], ...
    5);
```

This creates:

- `measurement_H_fast_HighOrder_x.mat`
- `measurement_H_fast_HighOrder_y.mat`
- `measurement_H_fast_HighOrder_z.mat`

## Processing

Run:

```matlab
GSTF_calculation/main_H_fast_HighOrder.m
```

Set `meas_name` and `ax_name` to the axis to be processed. The default field expansion is:

```matlab
calcChannels = 16;
```

Possible values are:

- 4: B0 + first order
- 9: B0 + first + second order
- 16: B0 + first + second + third order

The same spatial dataset can therefore be processed progressively with 4, 9, and 16 channels for regression and conditioning checks.

## Voxel selection

For spatial data (`numPE > 1`), `calculateOutputGradient.m` performs voxel selection before the spherical-harmonic fit.

The GUI provides:

- automatic magnitude-based voxel selection
- manual voxel add/remove
- FID magnitude display
- polarity-difference phase-derivative display
- rank of the selected probing matrix
- column-normalized condition number
- load/save of voxel selections

After confirmation, the selection is automatically saved next to the measurement file as:

```text
measurement_H_fast_HighOrder_x_voxelSelection.mat
```

On subsequent processing of the same dataset, the saved mask is loaded automatically when it matches the spatial grid and remains full rank for the requested `calcChannels`. Otherwise the GUI is opened again.

The original `numPE = 1` path retains the legacy Fast_GIRF 60% magnitude threshold.

## Spatial conditioning

Before the spherical-harmonic fit, the selected probing matrix is checked for:

```matlab
rank(A)
cond(A_normalized)
```

The column-normalized condition number is used only as a geometry/QC metric. The actual field fit remains the original unweighted least-squares/backslash solution:

```matlab
probingMatrix \ delta_diff_phase_avg
```

## Phase correction

`GIRF_matrix.correct_GSTF_phase()` accepts an optional reference channel. The default remains channel 2 for backwards compatibility. High-Order processing uses the physical self-term:

- x input: channel 2 (X)
- y input: channel 3 (Y)
- z input: channel 4 (Z)

## High-Order plotting

The original `GIRF_plotter.m` is kept unchanged and is still used for the self-term LF/HF comparison.

`plotHighOrderGSTFs.m` additionally plots all available terms grouped by spatial order:

- B0
- first order
- second order
- third order

The orders are separated because their GSTF units differ for a first-order input gradient:

- B0: m
- first order: dimensionless
- second order: 1/m
- third order: 1/m^2

## Regression tests

Run:

```matlab
tests/test_HighOrderSpatial.m
```

This test constructs a known 16-channel synthetic field on four slices and a `5 x 5` spatial grid, then verifies full rank and coefficient recovery for sagittal, coronal, and transverse orientations.

Run:

```matlab
tests/test_HighOrderSorting.m
```

This test assigns a unique value to every synthetic sequential ADC and verifies manifest-based placement into PE1, PE2, slice, triangle, polarity, and repetition dimensions. It also verifies that unacquired elliptical PE locations remain zero.

## Current limitations

- The current SH fit remains the original unweighted least-squares/backslash fit.
- The current LF/HF regularization is unchanged from the original Fast_GIRF implementation.
- The automatic voxel-selection threshold currently follows the original magnitude-based criterion and will be evaluated with real 7 T data before adding more complex quality metrics.
- CVP, EPI, and spiral validation are planned after the spatial acquisition/processing path is validated with real data.
