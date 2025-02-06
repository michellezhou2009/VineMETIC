# VineTMIC

This repository provides the R code for the simulation studies and data analysis in the manuscript "Modeling times to multiple events under informative censoring with $C$-vine copula" by Xinyuan Chen, Yiwei Li, and Qian M. Zhou.

- "runSimu1.R" and "executeSimu1.R" for the simulation study I (Section 4.1 of the manuscript)
- "runSimu2.R" and "executeSimu2.R" for the simulation study II (Section 4.2 of the manuscript)
- "runDataExample.R" for data example: the analysis of bivariate nonterminal event times (Section 5.2 of the manuscript) and trivariate nonterminal events (Section 5.3 of the manuscript)

The above three R files require some or all of the following R files: 

- "helpers.R", "helpers_bvic.R", "helpers_trivic.R", and "cvine_trivic.R" for the proposed stagewise estimation and inference procedures;

- "helpers_Li.R" for implementing the method by Li, D., Hu, X. J., and Wang, R. (2023), “Evaluating association between two event times with observations subject to informative censoring,” Journal of the American Statistical Association, 118, 1282–1294;

- "fitSPT.R" for fitting a semi-parametric transformation model;

- "DataGenFuns.R" for generating simulated data for simulation study I.



