# InformativeEmptiness
Code to run simulations for "Informatively empty clusters with application to multigenerational studies." 

## Simulations
#### simulation_poisson_slopes.R; simulation_ZIP_slopes.R 
Generates data from shared random effects models, with cluster sizes from Poisson/Zero-Inflated-Poisson distributions, specified conditionally and marginally. Fits conditional models via an outcome-only (naive) GLMM, joint conditional model excluding empty clusters, and full joint conditional model. Fits marginal models via GEE-exchangeable, IEE (independence estimating equations), WEE (inverse cluster-size weighted GEE), and joint marginalized model.

Default is RR=3,000 datasets with KK=4,000 clusters each. To test simulation and speed up process, set RR=100 and KK=1,000.

#### simulation_poisson.R; simulation_poisson.R
Same as above but for shared random intercepts models.
#### make_results_tables.R
Take results from simulations and generate summary tables and plots.
#### simulate_cluster_size_dist.R
Generate a population of 1,000,000 clusters in each of the simulation set-ups used above, and summarize cluster size distributions and outcome prevalences. Reports tables.

## Functions
Generic functions to generate data and fit models via proposed joint marginalized model.
#### JMMICSpack
Rcpp package containing functions used in fitting joint marginalized models (and joint conditional models).
#### JMMICS.R
Wrapper functions for fitting joint marginalized models (and joint conditional models).
#### genICS.R
Generate data with informative cluster size via Poisson or ZIP model, via shared random intercepts or random slopes. Specifies either a conditionally-specified or marginally-specified model.

## Data Analysis
Code for analysis of NHSII data (study of ADHD/Diethylstilbestrol). Note: data not publicly available, so code will not run.
#### analysis_run.R
Run data analysis. Calls "analysis_clean_dat.R" and "analysis_EDA.R". Fits conditional models via an outcome-only (naive) GLMM, joint conditional model excluding empty clusters, and full joint conditional model. Fits marginal models via GEE-exchangeable, IEE (independence estimating equations), WEE (inverse cluster-size weighted GEE), and joint marginalized model.
#### analysis_clean_dat.R
Clean NHSII data and format for model fitting.
#### analysis_EDA.R
Generate descriptive tables and plots for G1 and G2 level.
#### analysis_make_results_tables.R
Take results from "analysis_run.R" and report tables of results with estimates and confidence intervals.

