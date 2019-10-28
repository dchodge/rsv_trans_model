//
//  main.cpp
//  flu_forecast
//
//  Created by David Hodgson on 13/06/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

// Necessary includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <chrono>
#include <ctime>
#include <ctime>
#include <map>

// Preprocessor derivatives
#define PI 3.14159265   // Global definition of pi
#define Cal5            // Which detection model to use (calX for detection model X)
#define DynamicMat      // Which maternal immunity model to use (comment out -> static immunity)
//#define Minutes

using namespace std;

// Non-boost dependencies (in header files)
#include <math.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>

// Boost dependencies
#include <boost/random.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

// Prepackaged header files (look in header folder)
#include "Eigen/Dense"
#include "ascent/Ascent.h"
using namespace asc;
using namespace Eigen;

// Model header files
#include "pre.h"            // Soem preliminary stuff
#include "epmgp.h"          // Contains code for evalauting from a truncated multivariate normal distribution
#include "model.h"          // Description of model to be calibrated
#include "mcmc.h"           // Outlines the mcmc methods for calibration
#include "interventions.h"  // Intervention model code and scenarios
#include "outcomes.h"       // Calculate the outcomes
#include "cea.h"            // Cost-effectiveness analysis code
#include "writing.h"        // Code for writing out the results


// FUNC1:"parallel_temp_mcmc -> Run the MCMC parallel tempering algorithm
// Arg:"paramFit" -> Parameters to calibrate
// Arg:"Nburn" -> Burn-in
// Arg:"Nrun" -> Total length of run
// Arg:"Nthin_all" -> Thinning across all parameters (to view)
// Arg:"Nthin_pos" -> Thnning after burn-in (gives posterior)
// Arg:"Nchain" -> Number of chains used
// Arg:"A" -> Identifying character in output
void parallel_temp_mcmc(const std::vector<std::string>& paramFit, int Nburn, int Nrun, int Nthin_all, int Nthin_pos, int Nchain, char A)
{
    //Define the vector of mcmc_states and par_state
    std::vector<amh::amh_state_t> mcmc_state;
    std::vector<param::param_state_t> pars_r;
    mhp::mhp_state_t mhp_state = mhp::pre_initialize(Nchain, Nrun, Nthin_all);
    
    // Initialise the mcmc class spaces
    for (int j = 0; j < Nchain; j++)
    {
        amh::amh_state_t mcmc_state_t = mhp::run_MH::initialize_mcmc(paramFit, pars_r, Nburn, Nrun, Nthin_all, Nthin_pos);
        mcmc_state_t.T = mhp_state.T_full(0, j);
        mcmc_state.push_back(mcmc_state_t);
    }
    mhp_state = mhp::post_initialize(mhp_state, mcmc_state);
    
    for (int i = 1; i < Nrun ; i++)
    {
        // Swap chains adjacent chains
        mhp::update_swap(mhp_state, mcmc_state, i);

        // Find new position of each of the Markov chains using Metropolis hastings algorithm
        mhp::run_MH::update_mcmc(mhp_state, mcmc_state, pars_r);
        
        // Update the temperature ladders
        mhp::T_edit(mhp_state, mcmc_state);
        
        // Save at thinning points
        if(i%Nthin_all==0)
            mhp::save(mhp_state, mcmc_state);
        
    }
    
    // Write data for the 4 coldest chains
    amh::write_mcmc (mcmc_state[0], pars_r[0], 'A');
    amh::write_mcmc (mcmc_state[1], pars_r[1], 'B');
    amh::write_mcmc (mcmc_state[2], pars_r[2], 'C');
    amh::write_mcmc (mcmc_state[3], pars_r[3], 'D');

    // Write data on overall mcmc chains
    mhp::write_mhp_state (mhp_state, A);
}

// FUNC2:"posterior_inc -> Determine incidence from posterior samples
// Arg:"paramFit" -> parameters to calibrate
// Arg:"mcmc_state" -> class of the current mcmc state
// Arg:"pars" -> class of the current parameter values
// Arg:"seed" -> seed vector
void posterior_inc(std::vector<std::string> paramFit, amh::amh_state_t& mcmc_state, param::param_state_t& pars, num_vec seed)
{
    // WRITE THE POS AND SAMPLE SOLUTIONS (1000 SAMPLES FROM POSTERIOR)
    sim_ouput::write_inc (pars, mcmc_state, 'Z', seed);
    
    // WRITE THE R0 VALUES (not working, not needed)
    //sim_ouput::write_Rs(pars, mcmc_state, 100);
    
    // WRITE THE FORCE OF INFECTION (FOI) AND PR (proportion of new born with are protected)
    int s = 450; // seed value
    sim_ouput::write_foi(pars, mcmc_state, s);
    sim_ouput::write_pR(pars, mcmc_state, s);
    
}

// FUNC3:"find_optimal_week" -> Determine the optimal month for seasonal programmes to begin
// Arg:"pars" -> class of the current parameter values
// Arg:"mcmc_state" -> class of the current mcmc state
// Arg:"inter_data" -> structure of the intervention data
// Arg:"seed" -> seed vector
void find_optimal_week(param::param_state_t& pars, amh::amh_state_t& mcmc_state, cal::inter_data_t inter_data, num_vec seed)
{
    asc::Recorder recorder; // record the values (in ascent package)
    num_vec ind_s = {0, 2, 3, 4, 5, 6, 8, 10, 11, 12, 13, 14};
    str_vec sea_name = {"NONE","MAB_VHR_S","MAB_HR_S","MAB_HR_S+","MAB_ALL_S","MAB_ALL_S+","MAT_S","VAC_INF_S","VAC_2_4_S","VAC_5_9_S","VAC_5_14_S","VAC_75_S","VAC_65_S"};
    int s = 0;
    num_vec tot_incp(12,0);
    
    for (int w = 0; w < 52; w = w + 4)
    {
        for (int i = 0; i < 12 ; i++)
        {
            cea_state_t cea_state;

            num_vec up_take_base = cal::gen_daily(inter_data.uprate[ind_s[i]], w);
            inter_data.rate = up_take_base;
            inter_data.start_w[ind_s[i]] = w;
            inter_data.end_w[ind_s[i]] = w+21;

            cal::Calendar_full cal(inter_data, s, ind_s[i]);
            
            sim_ouput::int_post(pars, mcmc_state, cal, inter_data, cea_state, seed[s], ind_s[i], false);

            tot_incp[i] = cea_state.Q;
            cout << "Starting week number: " << w << ". Programme: " << sea_name[i] << ". Total QALY loss: " << tot_incp[i] << endl;
        }
        recorder({(double)w});
        recorder.add(tot_incp);
    }
    recorder.csv(get_ll::dout + "soln/" + "find_optimal", sea_name);
}

// FUNC5a:"intervention_p" -> Run the simualtions for the 14 intervention progammes
// Arg:"paramFit" -> parameters to calibrate
// Arg:"mcmc_state" -> class of the current mcmc state
// Arg:"pars" -> class of the current parameter values
// Arg:"inter_data" -> structure of the intervention data
// Arg:"seed" -> seed vector
void intervention_p(std::vector<std::string> paramFit, amh::amh_state_t& mcmc_state, param::param_state_t& pars, cal::inter_data_t inter_data, num_vec seed)
{
    // Intervention programme characteristics
    for (int i = 0; i < 16; i++)
    {
        cout << i << endl;
        sim_ouput::write_interventions(pars, mcmc_state, inter_data, seed, i);
    }
}

// FUNC5b:"intervention_p_SA" -> Run the simualtions for the sensitivtiy analysis
// Arg:"paramFit" -> parameters to calibrate
// Arg:"mcmc_state" -> class of the current mcmc state
// Arg:"pars" -> class of the current parameter values
// Arg:"inter_data" -> structure of the intervention data
// Arg:"seed" -> seed vector
void intervention_p_SA(std::vector<std::string> paramFit, amh::amh_state_t& mcmc_state, param::param_state_t& pars, cal::inter_data_t inter_data, num_vec seed)
{
    // Intervention programme characteristics
    num_vec ind_s1 = {0, 1, 2, 3, 4, 5, 6};
    num_vec ind_s2 = {0, 1, 2, 3, 4, 5, 6};
    num_vec ind_s3 = {0, 1, 2, 7, 8};

    inter_data.prog_no =  {"NONE_L_PRO", "PAL_VHR_S_L_PRO", "MAB_VHR_S_L_PRO", "MAB_HR_S_L_PRO", "MAB_HR_S+_L_PRO",  "MAB_ALL_S_L_PRO",  "MAB_ALL_S+_L_PRO"};
    inter_data.om_mab = 1.0/150.0;
    for (int i = 0; i < 7; i++)
    {
        cout << ind_s1[i] << endl;
        sim_ouput::write_interventions(pars, mcmc_state, inter_data, seed, ind_s1[i]);
    }
    
    inter_data.prog_no =   {"NONE_H_PRO", "PAL_VHR_S_H_PRO", "MAB_VHR_S_H_PRO", "MAB_HR_S_H_PRO", "MAB_HR_S+_H_PRO",  "MAB_ALL_S_H_PRO",  "MAB_ALL_S+_H_PRO"};
    inter_data.om_mab = 1.0/365.0;
    for (int i = 0; i < 7; i++)
    {
        cout << ind_s1[i] << endl;
        sim_ouput::write_interventions(pars, mcmc_state, inter_data, seed, ind_s2[i]);
    }

    inter_data.prog_no = {"NONE_H_MAT",   "PAL_VHR_S_H_MAT", "MAB_VHR_S_H_MAT", "MAB_HR_S_H_MAT", "MAB_HR_S+_H_MAT",  "MAB_ALL_S_H_MAT",  "MAB_ALL_S+_H_MAT",  "MAT_S_H_MAT",  "MAT_A_H_MAT"};
    inter_data.om_mab = 1.0/250.0;
    inter_data.xi_b = 0.75;

    for (int i = 0; i < 5; i++)
        sim_ouput::write_interventions(pars, mcmc_state, inter_data, seed, ind_s3[i]);

}

int main(int argc, const char * argv[]) {
    
    // Var:"paramFitA -> List of parameters the mcmc model is to fit. Parameters not in this list are fixed.
    std::vector<std::string> paramFitA = {"xi", "si", "g0", "g1", "g2", "om",
        "pA1", "pA2", "pA3", "pA4", "alpha_i", "d1", "d2", "d3",
        "phi", "qp", "qc", "b1",
        "c5ep1", "c5ep2",
        "ep5", "ep6",
        "seed1", "seed2", "psi"
    };
    
/**********************************/
/**      1. CALIBRATE MODEL         **/
/**********************************/
    // FUNC1:"parallel_temp_mcmc -> Run the MCMC parallel tempering algorithm
    // DATA REQUIRED: EVERYTHING IN THE INPUTS FOLDER (see model.h for importing of thesis data)
    //parallel_temp_mcmc(paramFitA, 25000, 50000, 100, 20, 12, 'A');
    
/**********************************/
/**      2. EVALUATE MODEL FIT        **/
/**********************************/
    // IMPORT POSTERIORS FROM CALIBRATION IN MCMC
    param::param_state_t pars(paramFitA);
    amh::amh_state_t mcmc_state(25000, 50000, 100, 20);
    mcmc_state = amh::initialize(pars, 25000, 50000, 100, 20);
    vector2D const post_in = get_ll::get_2d(get_ll::fileposterior_in, mcmc_state.NK2, pars.dim_cal);
    for (int d = 0; d < pars.dim_cal; d++)
        for (int s = 0; s < mcmc_state.NK2; s++)
            mcmc_state.posterior(s,d) = post_in[s][d];
    
    // SEED SET OF SAMPLES FROM THE POSTERIOR DISTRIBUTIONS
    num_vec seed;
    for (int s = 0; s < 10; s++)
        seed.push_back(uniform_dist_disc(0, mcmc_state.NK2-1, 'r'));
    
    // DETERMINE INCIDENCE SAMPLES FROM POSTERIOR, R0/REFF VALUES AND INCIDENCE
    // FUNC2:"posterior_inc -> Determine incidence from posterior samples
    // DATA REQUIRED: EVERYTHING IN THE INPUTS FOLDER (see model.h for importing of thesis data) AND THE POSTERIOR_A FILE
    //posterior_inc(paramFitA, mcmc_state, pars, seed);
    
/**********************************/
/**     3. FIND OPTIMAL WEEK            **/
/**********************************/
    // DATA REQUIRED: EVERYTHING IN THE INPUTS FOLDER (see model.h for importing of thesis data), EVERYTHING IN THE inter class AND THE POSTERIOR_A FILE
    cal::inter_data_t inter_data;
    inter_data.get_eff(inter_data, seed.size());
    
    // FUNC3:"find_optimal_week" -> Determine the optimal month for seasonal programmes to begin
    // DATA REQUIRED: EVERYTHING IN THE INPUTS FOLDER (see model.h for importing of thesis data) AND THE POSTERIOR_A FILE
    //find_optimal_week(pars, mcmc_state, inter_data, seed);
    
/**********************************/
/**     4. COMPARE CONSISTENCY (other file)         **/
/**********************************/
    
    // FUNC4a:"consistency_checks" -> XXXX
    // FUNC4b:"write_foi_mat" -> XXXX
    
    //sim_ouput::consistency_checks(pars, mcmc_state, inter_data);
    //sim_ouput::write_foi_mat(pars, mcmc_state, inter_data, 100);
    
/**********************************/
/**  5. EVALUATE INTERVENTION PROGRAMMES        **/
/**********************************/
    // FUNC5a:"intervention_p" -> Run the simualtions for the 14 intervention progammes
    // FUNC5b:"intervention_p_SA" -> Run the simualtions for the sensitivtiy analysis
    // DATA REQUIRED: EVERYTHING IN THE INPUTS FOLDER (see model.h for importing of thesis data), EVERYTHING IN THE inter class AND THE POSTERIOR_A FILE
    
    intervention_p(paramFitA, mcmc_state, pars, inter_data, seed);
    //intervention_p_SA(paramFitA, mcmc_state, pars, inter_data, seed);
}
