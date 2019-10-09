//
//  main.cpp
//  flu_forecast
//
//  Created by David Hodgson on 13/06/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <chrono>
#include <ctime>

#define PI 3.14159265
#define Cal5
#define DynamicMat
//#define Minutes

using namespace std;

#include <math.h>
#include <cmath>

#include <boost/random.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

#include <algorithm>    // std::max
#include <vector>
#include <random>

#include "Eigen/Dense"
#include "ascent/Ascent.h"

using namespace asc;
using namespace Eigen;

#include "pre.h"
#include "epmgp.h"
#include "model.h"
#include "mcmc.h"
#include "interventions.h"
#include "outcomes.h"
#include "cea.h"
#include "writing.h"


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
        mcmc_state_t.adap_covar = true;
        mcmc_state_t.T = mhp_state.T_full(0, j);
        mcmc_state.push_back(mcmc_state_t);
    }
    mhp_state = mhp::post_initialize(mhp_state, mcmc_state);
    
    for (int i = 1; i < Nrun ; i++)
    {
        ///auto wcts_2 = std::chrono::system_clock::now();
        // Swap chains adjacent chains
        mhp::update_swap(mhp_state, mcmc_state, i);

        // Find new position of each of the Markov chains using Metropolis hastings algorithm
        mhp::run_MH::update_mcmc(mhp_state, mcmc_state, pars_r);
        
        // Update the temperature ladders
        mhp::T_edit(mhp_state, mcmc_state);
        
        //
        if(i%Nthin_all==0)
            mhp::save(mhp_state, mcmc_state);
        
    }

    // Write data for the 4 coldest chains
    amh::write_mcmc (mcmc_state[0], pars_r[0], 'A');
    amh::write_mcmc (mcmc_state[1], pars_r[1], 'B');
    amh::write_mcmc (mcmc_state[2], pars_r[2], 'C');
    amh::write_mcmc (mcmc_state[3], pars_r[3], 'D');

    mhp::write_mhp_state (mhp_state, A);
}


void posterior_inc(std::vector<std::string> paramFit, amh::amh_state_t& mcmc_state, param::param_state_t& pars,  num_vec seed)
{
    // WRITE THE POS AND SAMPLE SOLUTIONS (1000 SAMPLES FROM POSTERIOR)
    sim_ouput::write_inc (pars, mcmc_state, 'Z', seed);
    
    // WRITE THE R0 VALUES (N SAMPLES)
    //sim_ouput::write_Rs(pars, mcmc_state, 100);
    
    // WRITE THE FOI AND PR (s= SEED VALUE)
    int s = 450;
    sim_ouput::write_foi(pars, mcmc_state, s);
    //sim_ouput::write_pR(pars, mcmc_state, s);
    
}

void find_optimal_week(param::param_state_t& pars, amh::amh_state_t& mcmc_state, num_vec seed, int time_hor, double disc)
{
    asc::Recorder recorder;
    VectorXd eff_pal = cal::get_eff(seed.size(), 1);
    VectorXd eff_mab =  cal::get_eff(seed.size(), 2);
    VectorXd eff_vac =  cal::get_eff(seed.size(), 3);
    VectorXd eff_mat =  cal::get_eff(seed.size(), 4);

    str_vec prog_no =   {"P_", "P1_", "P2_", "P3_", "P4_" , "P5_", "P7_", "P9_", "P10_", "P11_", "P12_", "P13_"};
    str_vec prog_no_out = {"P_", "P1_", "P2_", "P3_", "P4_" , "P5_", "P7_", "P9_", "P10_", "P11_", "P12_", "P13_"};
    str_vec prog_name = {"Base_", "mABHR_","mABHR_","mAB_" ,"mAB_" , "matS_","infS_","Pre"    ,"Sch1_" ,"Sch2_"  ,"Eld1_"   ,"Eld2_"};
    
    //Intervention calendar generation for each programme
    str_vec cal_type = {"None", "Mhr","Mhr_p","Mlr","Mlr_p", "mat", "LAV_inf", "LAV_ald", "LAV_ald", "LAV_ald", "LAV_ald", "LAV_ald"};
    vector2D t_group = {cal::G_base, cal::G_0mo, cal::G_0mo, cal::G_0mo, cal::G_0mo, cal::G_par, cal::G_2mo, cal::G_2_4, cal::G_5_10, cal::G_5_14, cal::G_75_, cal::G_65_};
    vector2D uprate = {{0},{0},{0},{0},{0},{0},{0}, cal::up_week_2t3, cal::up_week_u65, cal::up_week_u65, cal::up_week_o65, cal::up_week_o65};
    num_vec cov =     {0.0,0.9,0.9,0.9,0.9,0.6,0.9, 0.45, 0.6, 0.6, 0.7, 0.7};
    num_vec cov_c =   {0.0,0.0,0.0,0.0,0.0,0.6,0.0,0.0,0.0,0.0,0.0,0.0};
    num_vec c_ad =  {0.0, 11.5, 11.5, 11.5, 11.5, 9, 9, 9,9, 9, 9, 9, 9, 9};
    num_vec Pal_ind = {false,false,false,false,false,true,true, true, true, true, true, true};
    num_vec inc;
    
    num_vec S_tot(NoAgeG, 0); num_vec H_tot(NoAgeG, 0); num_vec D_tot(NoAgeG, 0);
    num_vec GP_tot(NoAgeG, 0); num_vec BD_tot(NoAgeG, 0);
    vector2D doses(2,vector<double>(time_hor*52,0));
    double Q, CP, CT;
    
    int s = 0;
    num_vec tot_incp(12,0);
    for (int w = 0; w < 52; w = w + 4)
    {
        for (int iN = 0; iN < 12 ; iN++)
        {
            Q = CP = CT = 0;
            for (int a = 0; a < NoAgeG; a++)
            {
                S_tot[a] = 0; H_tot[a] = 0; D_tot[a] = 0;
                GP_tot[a] = 0; BD_tot[a] = 0;
            }
            for (int t = 0; t < time_hor*52; t++)
            {
                doses[0][t] = 0; doses[1][t] = 0;
            }
            
            num_vec up_take_base = cal::gen_daily(uprate[iN], w);
            cal::Calendar_full cal(t_group[iN], cov[iN], up_take_base, w, w+21, Pal_ind[iN], cal_type[iN], s, eff_pal, eff_mab, eff_vac, eff_mat);
            
            num_vec inciall = sim_ouput::int_post(pars, mcmc_state, cal, false, seed[s], cov_c[iN], S_tot, H_tot, D_tot, GP_tot, BD_tot, doses, Q, CP, CT, c_ad[iN], time_hor, disc, 1.0/250.0, 1);

            tot_incp[iN] = Q;
            cout << tot_incp[iN] << endl;
        }
        recorder({(double)w});
        recorder.add(tot_incp);
    }
    recorder.csv(get_ll::dout + "soln/" + "find_optimal" , prog_no_out);
}


void intervention_p(std::vector<std::string> paramFit, amh::amh_state_t& mcmc_state, param::param_state_t& pars,  num_vec seed, int time_hor, double disc)
{
    // Intervention programme identifiers
    //{{{26}}, {{26}}, {{4}}, {{46}}, {{4}}, {{13}}, {{17}}, {{21}}, {{21}}}
    // Optimal fit is XXXXX //
    VectorXd eff_pal = cal::get_eff(seed.size(), 1);
    VectorXd eff_mab =  cal::get_eff(seed.size(), 2);
    VectorXd eff_vac =  cal::get_eff(seed.size(), 3);
    VectorXd eff_mat =  cal::get_eff(seed.size(), 4);
    
    for (int i = 0; i < 16; i++)
    {
        cout << i << endl;
        sim_ouput::write_interventions(pars, mcmc_state, cal::prog_no[i], cal::prog_name[i], seed, cal::cov_c[i], cal::c_ad[i], time_hor, eff_pal, eff_mab, eff_vac, eff_mat, i, disc, 1.0/250.0, 1);
    }
}

void intervention_p_SA(std::vector<std::string> paramFit, amh::amh_state_t& mcmc_state, param::param_state_t& pars,  num_vec seed, int time_hor, double disc)
{
    // Intervention programme identifiers
    VectorXd eff_pal = cal::get_eff(seed.size(), 1);
    VectorXd eff_mab =  cal::get_eff(seed.size(), 2);
    VectorXd eff_vac =  cal::get_eff(seed.size(), 3);
    VectorXd eff_mat =  cal::get_eff(seed.size(), 4);
    
    for (int i = 0; i < 7; i++)
    {
        cout << i << endl;
        sim_ouput::write_interventions(pars, mcmc_state, cal::prog_no_SA_1[i], cal::prog_name[i], seed, cal::cov_c[i], cal::c_ad[i], time_hor, eff_pal, eff_mab, eff_vac, eff_mat, i, disc, 1.0/150.0, 1);
        sim_ouput::write_interventions(pars, mcmc_state, cal::prog_no_SA_2[i], cal::prog_name[i], seed, cal::cov_c[i], cal::c_ad[i], time_hor, eff_pal, eff_mab, eff_vac, eff_mat, i, disc, 1.0/365.0, 1);
    }

    for (int i = 0; i < 2; i++)
        sim_ouput::write_interventions(pars, mcmc_state, cal::prog_no_SA_3[i], cal::prog_name[i], seed, cal::cov_c[i], cal::c_ad[i], time_hor, eff_pal, eff_mab, eff_vac, eff_mat, i, disc, 1.0/250.0, 0.75);

    for (int i = 7; i < 9; i++)
        sim_ouput::write_interventions(pars, mcmc_state, cal::prog_no_SA_3[i], cal::prog_name[i], seed, cal::cov_c[i], cal::c_ad[i], time_hor, eff_pal, eff_mab, eff_vac, eff_mat, i, disc, 1.0/250.0, 0.75);

}

void test_random(int s)
{
    PRNG_s rng(s);
    boost::random::weibull_distribution<> Ev(11.8975, 0.7317);
    cout << Ev(rng) << endl;
    cout << Ev(rng) << endl;
}

int main(int argc, const char * argv[]) {
        
    // arguments (vector string listing the parameters to be fitted,
        // number of steps of burn in,
        //  number of steps overall,
        //  thinning of whole chain,
       // thinning of posterior,
        // number of parallel chains,
        // year data is fitted to,
        // Char for output )
    std::vector<std::string> paramFitA = {"xi", "si", "g0", "g1", "g2", "om",
        "pA1", "pA2", "pA3", "pA4", "alpha_i", "d1", "d2", "d3",
        "phi", "qp", "qc", "b1",
        "c5ep1", "c5ep2",
        "ep5", "ep6",
        "seed1", "seed2", "psi"
    };
    
/**********************************/
/**      CALIBRATE MODEL         **/
/**********************************/
    //parallel_temp_mcmc(paramFitA, 25000, 50000, 100, 20, 12, 'A');
/**********************************/
/**      IMPORT POSTERIOR FOR REST        **/
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
    int time_hor = 10;
    double disc = 0.035;
/**********************************/
/**      EVALUATE MODEL FIT        **/
/**********************************/
    // DETERMINE INCIDENCE SAMPLES FROM POSTERIOR, R0/REFF VALUES AND INCIDENCE
   // posterior_inc(paramFitA, mcmc_state, pars, seed);
   // output_fitted_dist();
    
/**********************************/
/** FIND OPTIMAL WEEK            **/
/**********************************/
    // DETERMINE THE OPTIMAL WEEEK FOR WHICH TO START THE TREATMENT
    find_optimal_week(pars, mcmc_state, seed, time_hor, disc);
    
/**********************************/
/** COMPARE CONSISTENCT         **/
/**********************************/
    
    //sim_ouput::consistency_checks(pars, mcmc_state, 1, 1.0/250.0, 1.0);
    //sim_ouput::write_foi_mat(pars, mcmc_state, 100);
    
/**********************************/
/**  EVALUATE INTERVENTION PROGRAMMES        **/
/**********************************/
    // RUN THE INTERVENTION PROGRAMMES

    //intervention_p(paramFitA, mcmc_state, pars, seed, time_hor, disc);
    //intervention_p_SA(paramFitA, mcmc_state, pars, seed, time_hor, disc);
}
