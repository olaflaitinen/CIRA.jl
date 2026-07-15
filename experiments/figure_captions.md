# Figure captions

Descriptions are placed in captions, not inside the figures, following the target
journal style. Figures are provided as high-resolution PNG (400 dpi) and vector PDF.
File names use the journal numbering convention (Fig1, Fig2, Fig3).

## Fig1

Computational scaling of the CIRA reference pipeline. Analysis runtime as a
function of the number of parameters q for a family of linear catenary models,
plotted on logarithmic axes. Markers are the fastest of three repetitions; the
dashed line is a power-law fit whose exponent is reported in the legend. The
near-quadratic scaling reflects the singular value decomposition of the local
sensitivity matrix.

## Fig2

Structural and practical identifiability. (a) Relative standard errors of the two
SIR parameters at five percent measurement noise, confirming that the globally
identifiable model yields finite, well-constrained estimates. (b) Fisher-information
condition number for the seven structurally non-identifiable benchmarks, comparing
the original parameterization with the reduced parameterization returned by CIRA;
the reduction removes the near-singular directions and lowers the condition number
by many orders of magnitude. (c) Growth of the SIR relative standard errors as a
function of measurement noise.

## Fig3

Robustness of the identifiability verdict to the integration grid. (a) Recovered
numerical rank for every benchmark across integration grids ranging from 0.75x to
3x the baseline density; cells are annotated with the rank and shaded green when it
matches the reference value. (b) Smallest identifiable singular value, normalized by
the largest singular value, as a function of grid density. All identifiable
directions remain far above the numerical null space (shaded band), so the verdict
is unchanged across grids.
