# Introduction 

This is a repository containing the R code for simulation studies conducted in the project "Spatial Dependence Modeling". 
There are two set of R codes, namely sim1 and sim2, corresponding to each simulations study, respectively. 

## sim1

utility functions are included in commonf.cpp; 

likelihood function for estimation and inference purposes are included in loglik.cpp and ts_inference_fct.R; 

dataGen_fct.R: functions for simulating data sets;

sim_script.R: the running script for generating results presented in Table 2 and Web Table 1. It displays the flow of a simulation run: 1) specifying the parameter settings; 2) solving for unknown parameters; 3) generating a dataset; 4) implement the two-stage estimation procedure and 5) conduct the inference on parameters of interest. 


## sim2

prepare_param_sim.R: functions for calculating parameters in different configurations in the 2nd set of simulation studies; 

dataGen_fcts.R: functions for simulating data sets;

TwoStage_est_fct.R, TwoStage_inference_fct.R: functions for estimation and inference purposes;

sim_script.R: the running script for generating results presented in Table 3. Similarly, it presents the flow of a simulation run: 1) specifying the parameter settings; 2) solving for unknown parameters; 3) generating a dataset; 4) implement the two-stage estimation procedure and 5) conduct the inference on parameters of interest. 


## PseudoScoreTest.R

Functions for conducting pseudo-score tests based on working independence composite likelihoods. 

## sample_data.txt

A toy data set showing the simulated data structure: 

- id: subject id;

- type: joint type, j=1, ..., J;

- joint: k=1, ..., $K_j$;

- x: a binary covariate representing a demographic feature such as gender; 

- g: a binary genetic marker;

- l: the left end of the censored interval; 

- r: the right end of the censored interval; 

- truez: true values of the (partially) latent susceptibility indicators. 



