# Data generation
## data.generation.R

Generates 200 datasets following specifications outlined in the Simulation Data Generation section. In order to create the 4 data generating mechanisms under which 200 sets of replicated datasets are simulated, we may toggle the number of instrumental/prognostic/noise variables which are included in the covariate matrix output, though the core data generating mechanism is always the same.

# Run code
## designfunction.R

Included function performs the PS estimation in a MLE or Bayesian manner. Requires treatment (vector), covariates (matrix or data.frame where column number = number of covariates). User may specify how many draws from the posterior distribution of PS they wish to output (K). 

Outputs estimated PS in a n by K matrix format, where n is the sample size of the inputted data and K is 1 if PS estimation is performed in a MLE manner.

## analysisfunction.R 

Included function (analysisfxn) performs one of 5 (user-specified) implementations + analyses 

Outputs estimated ATEs, estimated SEs, and implementation/analysis details such as number of matches, proportion of treated in each strata etc. Some outputs are implementation specific. 

## BPSAfxn.R

Included function performs BPSA (calls designfunction.R and analysisfunction.R) for all 5 implementations, assembles outputs for all implementations.

Outputs bpsa.result, a large matrix with row names:
 
- ATE (estimate of delta)
- between.design.SE (SD of ATE estimates across K iterations)
- within.design.SE (average of K SE estimates)
- within.design.sandwich.SE (average of K SE estimates, calculated with a robust sandwich estimator)
- method (implementation method)
- rubinse (marginal SE of delta)
- interval.lower (lower bound for a 95% credible interval)
- inteval.upper (upper bound for a 95% credible interval)
- percused (percentage of data used, only relevant for matching implementations)
- marginal.SD.nu (marginal SD of nu, averaged over within-observation standard deviation of nu across K draws of the PS)
- marginal SD.nu.treated (marginal SD of nu for treated observations)
- marginal.SD.nu.control (marginal SD of nu for control observations)
- ncontrol (number of control observations in matched set, only relevant for matching algorithms)
- ntreat (number of treated observations in matched set, only relevant for matching algorithms)

## PSAfxn.R
Performs PSA (calls designfunction.R and analysisfunction.R) for all 5 implementations, assembles outputs for all implementations. Used for primary simulation.

Outputs psa.result, a matrix with row names:

- ATE (estimated treatment effect)
- se (SE estimated via asymptotic methods)
- method (implementation method)
- percused (percentage of data used, only relevant for matching implementations)
- ncontrol (number of control observations in matched set, only relevant for matching algorithms)
- ntreat (number of treated observations in matched set, only relevant for matching algorithms)
- cov.balance (measure covariate balance, calculated via standardized difference in covariates between treated and control, averaged over covariates)

# Code for application

## app.data.R

List object containing "treatment" (vector length 22,723), "covariates" 22,723 by 17 matrix of covariates, "outcome" (vector length 22,723)

## data.manip.R

Reads in original source_oriented_data_noMEDICARE.csv, performs a complete case analysis, restricts regions to Northeast, Southeast and Industrial Midwest, dichotomizes the exposure, removes some covariates and re-organizes relevant data into a list app.data (written as app.data.R).

## run.app.R

Calls on PSAfxn.R and BPSAfxn.R, performs both PSA and BPSA on app.data.R

# Data for application

## source_oriented_data_noMEDICARE.csv

Data from Cummisky et. al. stripped of Medicare data, with simulated values for the outcome. 

# Figure code

## p1_sim_graphs.R

Code to create Figures 1-3, Table B1

## appfig.R

Code to create Figures 4 and 5

