# Simulation study for volume-outcome analyses

<!-- badges: start -->
<!-- badges: end -->

This R code provides tools for volume-outcome analyses
as they are proposed in the methodological paper and to ascertain the
reproducibility of the simulation study conducted in the paper.

## simulation\_generate.R

This file contains code to simulate and save the data to be used in the
analysis. One may freely adjust all parameters introduced in the paper.

## simulation\_analyze.R

A volume-outcome-analysis is conducted on the data generated
through `simulation_generate.R`. This includes coefficient estimates,
p-values for the volume effect and odds ratios with respect to the
volume effect. All results can be saved.

## simulation\_plot.R

This file contains code for constructing and exporting a number of plots
based on the results obtained in `simulation_analyze.R`. It covers the
estimated (versus underlying) volume effects and corresponding p-values,
random effects standard deviations and odds ratios with respect to the
volume effects.

Successively executing `simulation_generate.R`,
`simulation_analyze.R` and `simulation_plot.R` produces the figures as
included in the paper.

## License

MIT
