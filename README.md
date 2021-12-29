RSV Transmission model

# THIS CODE HAS BEEN SUPERSEEDED
The fitting part of the code has an R implementation: [dchodge/rsv_model_R_fit](https://github.com/dchodge/rsv_model_R_fit)
The intevention part of the code has an R implementaiton here: [dchodge/rsv_model_R_custom](https://github.com/dchodge/rsv_model_R_custom)

# David Hodgson's RSV transmission model 

This is the c++ code is used to determine the cost-effectiveness of different intervention programmes agaisnt RSV infection in England. The model comprises three main parts:
1) A code to calibrate the dynamic transmission model to RSV positive samples using a parallel tempering MCMC method
2) Codes to evaluate the impact of the intervention prgorammes at preventing outcomes, including symptomatic infections, hospitalised cases, deaths, GP consultations, and bed days.
3) A cost-effectiveness analysis to convert the outcomes averted into QALY loss and costs.

## LICENSE
Copyright (C) 2019  David Hodgson

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Installation
A c++14 compliant complier is required

## Dependencies 
A header folder should be made and the following header-only packages should be placed there:
- Ascent https://github.com/AnyarInc/Ascent
- Eigen http://eigen.tuxfamily.org/index.php?title=Main_Page

Boost is also required, and the boost_random package should be linked. 

# Outline of .cpp and .h files in src folder
## Main.cpp
This section outlines the five main functions;

FUNC1: parallel_temp_mcmc -> Run the MCMC parallel tempering algorithm
FUNC2:"posterior_inc -> Determine incidence from posterior samples
FUNC3:"find_optimal_week" -> Determine the optimal month for seasonal programmes to begin
FUNC5a:"intervention_p" -> Run the simualtions for the 14 intervention progammes
FUNC5b:"intervention_p_SA" -> Run the simualtions for the sensitivtiy analysis

I comment out the functions I don't want and run then build and run. Probably very easy to make an interface for this. 

## pre.h
A few type defs, probability distributions functions, and global functions are defined here 

## epmgp.h
This is just a c++ versions of a script outlined here >>> https://github.com/cunni/epmgp. It just helps evaluated a point x in a multivariate truncated normal distribution with support: a<x<b, mean: mu, and covariance matrix Sigma. 

## mcmc.h
Code for the adaptive parallel tempering mcmc methods >> https://arxiv.org/abs/1205.1076. amh namespace refers to metropolis hasting on each of the chains, mhp namespace refers to parameters relevant to all the chains

## model.h
Priors, ODES, likelihood, data import for the calibration dynamic transmission model. all data is in the inputs folder.

## outcomes.h
Code used to determine the probabilities of each outcome.

## interventions.h
Code outlining the characteristics of each of the intervention programmes, generating the dynamic calendar, and the ODEs for the augmented model 

## writing.h
Formating the vectors into a exportable format (see outputs folder)

## cea.h
code for converting the outcomes averted into QALY and costs

