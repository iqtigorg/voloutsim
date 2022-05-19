# Simulation study for volume-outcome analyses

<!-- badges: start -->
<!-- badges: end -->

This R code provides tools for volume-outcome analyses
as they are proposed in the paper
[*Modelling volume-outcome relationships in health care*](https://arxiv.org/abs/2203.12927v1)
and to ascertain the reproducibility of the simulation study
conducted in that paper.

## Files

The simulation is contained in four R scripts.
Successively executing `simulation_generate.R`,
`simulation_analyze.R` and `simulation_plot.R` produces the numbers and figures as
included in the paper.

1. `simulation_generate.R`:
    This file contains code to simulate data to be used in the analysis. One may
    freely adjust all parameters introduced in the paper.

2. `simulation_analyze.R`:
    A volume-outcome-analysis is conducted on the data generated
    through `simulation_generate.R`. This includes coefficient estimates,
    p-values for the volume effect and odds ratios with respect to the
    volume effect.

3. `simulation_plot.R`:
    This file contains code for constructing and exporting a number of plots
    based on the results obtained in `simulation_analyze.R`. It covers the
    estimated (versus underlying) volume effects and corresponding p-values,
    random effects standard deviations and odds ratios with respect to the
    volume effects.

The full simulation study with all combinations of parameter values as in the
paper takes around two to three days on a server with 16 cores at 2.9 GHz.
For a single combination of the parameter values, generating and analyzing
the data can be done within a few minutes on a normal desktop computer.

## System requirements

The code assumes a recent version of [R](https://www.r-project.org/) and makes
use of the libraries `mgcv`, `dplyr`, `tidyr`, `ggplot2`, `scales`,
`extrafont` and `sessioninfo`.
The simulation results in the paper were obtained with R-4.0.2 on a
Windows server (x64).

## License

MIT
