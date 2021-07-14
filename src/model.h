//
//  param.h
//  flu_forecast
//
//  Created by David Hodgson on 17/06/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

#ifndef model_h
#define model_h


using namespace boost::math;

namespace param
{
    // Everything in this strcture outlines the five differen detection model structures (outlined in Appendix)
    struct param_state_t
    {
        str_vec pA_vec = {"pA1", "pA1", "pA1", "pA1", "pA1", "pA1", "pA1", "pA1", "pA1", "pA1", "pA1", "pA1", "pA2", "pA2", "pA2", "pA2", "pA3", "pA3", "pA4", "pA4", "pA4", "pA4", "pA4", "pA4", "pA4"
        };
        
        str_vec c1m = {
            "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1", "c1ep1",
            "ep5", "ep5", "ep5", "ep5", "ep5", "ep5",
            "ep6", "ep6", "ep6"
        };
        
        str_vec c2m = {
            "c2ep1", "c2ep1", "c2ep1", "c2ep2", "c2ep2", "c2ep2", "c2ep2", "c2ep2", "c2ep2", "c2ep2", "c2ep2", "c2ep2", "c2ep2", "c2ep2", "c2ep2", "c2ep2",
            "ep5", "ep5", "ep5", "ep5", "ep5", "ep5",
            "ep6", "ep6", "ep6"
        };
        
        str_vec c3m = {
            "c3ep1", "c3ep1", "c3ep1", "c3ep2", "c3ep2", "c3ep2", "c3ep2", "c3ep2", "c3ep3", "c3ep3", "c3ep3", "c3ep3", "c3ep3", "c3ep3", "c3ep3", "c3ep3",
            "ep5", "ep5", "ep5", "ep5", "ep5", "ep5"
            "ep6", "ep6", "ep6"
        };
        
        str_vec c4m = {
            "c4ep1", "c4ep1", "c4ep1", "c4ep2", "c4ep2", "c4ep2", "c4ep3", "c4ep3", "c4ep3", "c4ep3", "c4ep3", "c4ep3", "c4ep4", "c4ep4", "c4ep4", "c4ep4",
            "ep5", "ep5", "ep5", "ep5", "ep5", "ep5",
            "ep6", "ep6", "ep6"
        };
        
        str_vec c5mA = {
            "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1", "c5ep1",
            "ep5", "ep5", "ep5", "ep5", "ep5", "ep5",
            "ep6", "ep6", "ep6"
        };
        
        str_vec c5mB = {
            "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2", "c5ep2"
        };
        
        // vector names of names of parameters
        str_vec paramAll =   {"xi", "si", "g0", "g1", "g2", "om",
            "pA1", "pA2", "pA3", "pA4", "rho", "alpha_i", "d1", "d2", "d3",
            "phi", "qp", "qc", "b1",
            "c1ep1",
            "c2ep1", "c2ep2",
            "c3ep1", "c3ep2", "c3ep3",
            "c4ep1", "c4ep2", "c4ep3", "c4ep4",
            "c5ep1", "c5ep2",
            "ep5", "ep6",
            "seed1", "seed2","psi"
        };
        
        // Lower bound for the 25 parameters (for the truncated normal, values cannot be less than these values )
        num_vec Low_b_all =  {  14,      2,      2,      0.5,      0.5,     60,    // xi, si, g0, g1, g2, om
            0,     0,     0,    0,      0.0,      0.05,      0.0,      0.0,     0.0,         // pA1, pA2, pA3, pA4, rho, alpha_i, d1, d2, d3
            0,     0,     0,    0,                                          // phi, qp, qc, b1
            0,                                                              // c1ep1
            0,     0,                                                       // c2ep1, c2ep2,
            0,     0,     0,                                                // c3ep1, c3ep2, c3ep3,
            0,     0,     0,     0,                                         // c4ep1, c4ep2, c4ep3,, c4ep4,
            -10,     -10,                                                       // c5ep1, c5ep2
            0,     0,                                                       // ep5, ep6
            0,     0,  0                                                      // seed1, seed2
        };
        
        // Upper bound for the 25 parameter for the truncated normal, values cannot exceed these values )
        num_vec Upp_b_all =  {    180,     8,     20,    1,      1,      365,    // xi, si, g0, g1, g2, om
            1,     1,     1,      1,     1,     0.95,    1,      1,      1,         // pA1, pA2, pA3, pA4, rho, alpha_i, d1, d2, d3
            1,     0.1,     1,    2,                                                // phi,  qp, qc, b1
            1,                                                              // c1ep1
            1,     1,                                                       // c2ep1, c2ep2,
            1,     1,     1,                                                // c3ep1, c3ep2, c3ep3,
            1,     1,     1,     1,                                         // c4ep1, c4ep2, c4ep3,, c4ep4,
            10,     10,                                                     // c5ep1, c5ep2
            1.0,     1.0,                                                    // ep5, ep6
            1,     1,    1                                                  // seed1, seed2
        };
        
        str_vec paramFit;

        std::map<std::string, double> mapofParam;
        std::map<std::string, double> mapofParamUall;
        std::map<std::string, double> mapofParamLall;
        int dim_full = 36;
        
/** Useful matrices for the support */
        VectorXd support_range;
        int dim_cal;
        bool OOS;
        int year;
    
        param_state_t()
        {
            for (int i = 0; i < dim_full; i++)
            {
                mapofParam.insert(std::make_pair(paramAll[i], 0));
                mapofParamUall.insert(std::make_pair(paramAll[i], Upp_b_all[i]));
                mapofParamLall.insert(std::make_pair(paramAll[i], Low_b_all[i]));
            }
            OOS = false;
        }

        param_state_t(std::vector<std::string> paramFit_t): paramFit(paramFit_t)
        {
            for (int i = 0; i < dim_full; i++)
            {
                mapofParam.insert(std::make_pair(paramAll[i], 0));
                mapofParamUall.insert(std::make_pair(paramAll[i], Upp_b_all[i]));
                mapofParamLall.insert(std::make_pair(paramAll[i], Low_b_all[i]));
            }
            OOS = false;
            dim_cal = paramFit.size();
            support_range = Eigen::VectorXd::Zero(dim_cal);
            for (int i = 0; i < dim_cal; i++)
            {
                support_range(i) =  mapofParamUall.find(paramFit[i])->second -  mapofParamLall.find(paramFit[i])->second;
            }
            
        }
        
        // FUNC: lprior_dist ->  get values for the prior out
        // ARG: par -> parameter for which prior is needed
        // ARG: pars -> structure of parameter values
        // ARG: sample -> false: evaluate at point in prior, true: sample a point from dist

        double lprior_dist(string par, param::param_state_t& pars, bool sample)
        {
            if (par == "xi")
            {
                if (sample){boost::random::uniform_real_distribution<double> U_xi(14.0,180.0); return U_xi(rng);}
                else{       boost::math::uniform_distribution<double> U_xi(14.0,180.0); return pdf(U_xi, pars.mapofParam.find(par)->second);}
            }
            else if (par == "si")
            {
                if (sample){boost::random::gamma_distribution<double> G_si(7.111, 0.563); return G_si(rng);}
                else{       boost::math::gamma_distribution<double> G_si(7.111, 0.563); return pdf(G_si, pars.mapofParam.find(par)->second);}
            }
            else if (par == "g0")
            {
                if (sample){boost::random::weibull_distribution<double> W_g0(4.137, 8.303); return W_g0(rng);}
                else{       boost::math::weibull_distribution<double> W_g0(4.137, 8.303); return pdf(W_g0, pars.mapofParam.find(par)->second);}
            }
            else if (par == "g1")
            {
                if (sample){boost::random::weibull_distribution<double> W_g1(34.224, 0.879); return W_g1(rng);}
                else{       boost::math::weibull_distribution<double> W_g1(34.224, 0.879); return pdf(W_g1, pars.mapofParam.find(par)->second);}
            }
            else if (par == "g2")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_g2(-0.561, 0.163); return LN_g2(rng);}
                else{       boost::math::lognormal_distribution<double> LN_g2(-0.561, 0.163); return pdf(LN_g2, pars.mapofParam.find(par)->second);}
            }
            else if (par == "om")
            {
                if (sample){boost::random::normal_distribution<double> N_om(200.0, 35.0); return N_om(rng);}
                else{       boost::math::normal_distribution<double> N_om(200.0, 35.0); return pdf(N_om, pars.mapofParam.find(par)->second);}
            }
            else if (par == "pA1")
            {
                if (sample){boost::random::beta_distribution<double> B_pA1(3.003, 29.997); return B_pA1(rng);}
                else{       boost::math::beta_distribution<double> B_pA1(3.003, 29.997); return pdf(B_pA1, pars.mapofParam.find(par)->second);}
            }
            else if (par == "pA2")
            {
                if (sample){boost::random::beta_distribution<double> B_pA2(8.996, 43.004); return B_pA2(rng);}
                else{       boost::math::beta_distribution<double> B_pA2(8.996, 43.004); return pdf(B_pA2, pars.mapofParam.find(par)->second);}
            }
            else if (par == "pA3")
            {
                if (sample){boost::random::beta_distribution<double> B_pA3(38.033, 34.967); return B_pA3(rng);}
                else{       boost::math::beta_distribution<double> B_pA3(38.033, 34.967); return pdf(B_pA3, pars.mapofParam.find(par)->second);}
            }
            else if (par == "pA4")
            {
                if (sample){boost::random::beta_distribution<double> B_pA4(35.955, 11.045); return B_pA4(rng);}
                else{       boost::math::beta_distribution<double> B_pA4(35.955, 11.045); return pdf(B_pA4, pars.mapofParam.find(par)->second);}
            }
            else if (par == "d1")
            {
                if (sample){boost::random::beta_distribution<double> B_d1(35.583, 11.41); return B_d1(rng);}
                else{       boost::math::beta_distribution<double> B_d1(35.583, 11.41); return pdf(B_d1, pars.mapofParam.find(par)->second);}
            }
            else if (par == "d2")
            {
                if (sample){boost::random::beta_distribution<double> B_d2(23.209, 2.790); return B_d2(rng);}
                else{       boost::math::beta_distribution<double> B_d2(23.209, 2.7909); return pdf(B_d2, pars.mapofParam.find(par)->second);}
            }
            else if (par == "d3")
            {
                if (sample){boost::random::beta_distribution<double> B_d3(8.352, 10.648); return B_d3(rng);}
                else{       boost::math::beta_distribution<double> B_d3(8.352, 10.648); return pdf(B_d3, pars.mapofParam.find(par)->second);}
            }
            else if (par == "phi")
            {
                if (sample){boost::random::uniform_real_distribution<double> U_phi(0.35,0.55); return U_phi(rng);}
                else{       boost::math::uniform_distribution<double> U_phi(0.0, 1.0); return pdf(U_phi, pars.mapofParam.find(par)->second);}
            }
#ifdef Minutes
            else if (par == "qp")
            {
                if (sample){boost::random::uniform_real_distribution<double> U_R0(0.0, 0.001); return U_R0(rng);}
                else{       boost::math::uniform_distribution<double> U_R0(0.0, 0.001); return pdf(U_R0, pars.mapofParam.find(par)->second);}
            }
#else
            else if (par == "qp")
            {
                if (sample){boost::random::uniform_real_distribution<double> U_R0(0.0, 0.1); return U_R0(rng);}
                else{       boost::math::uniform_distribution<double> U_R0(0.0, 0.1 ); return pdf(U_R0, pars.mapofParam.find(par)->second);}
            }
#endif
            else if (par == "alpha_i")
            {
                if (sample){boost::random::normal_distribution<double> N_alp(0.5, 0.1); return N_alp(rng);}
                else{       boost::math::normal_distribution<double> N_alp(0.5, 0.1); return pdf(N_alp, pars.mapofParam.find(par)->second);}
            }
            else if (par == "b1")
            {
                if (sample){boost::random::uniform_real_distribution<double> U_b1(0.0, 1.0); return U_b1(rng);}
                else{       boost::math::uniform_distribution<double> U_b1(0.0, 10.0); return pdf(U_b1, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c1ep1")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_ep1(-4.92843, 0.227847); return LN_ep1(rng);}
                else{       boost::math::lognormal_distribution<double> LN_ep1(-4.92843, 0.227847); return pdf(LN_ep1, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c2ep1")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_ep1(-4.26992, 0.0969545); return LN_ep1(rng);}
                else{       boost::math::lognormal_distribution<double> LN_ep1(-4.26992, 0.0969545); return pdf(LN_ep1, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c2ep2")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_ep2(-5.70866, 0.290209); return LN_ep2(rng);}
                else{       boost::math::lognormal_distribution<double> LN_ep2(-5.70866, 0.290209); return pdf(LN_ep2, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c3ep1")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_ep1(-4.269, 0.0969); return LN_ep1(rng);}
                else{       boost::math::lognormal_distribution<double> LN_ep1(-4.269, 0.0969); return pdf(LN_ep1, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c3ep2")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_ep2(-5.2330, 0.18086); return LN_ep2(rng);}
                else{       boost::math::lognormal_distribution<double> LN_ep2(-5.2330, 0.18086); return pdf(LN_ep2, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c3ep3")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_ep3(-6.2009, 0.3319); return LN_ep3(rng);}
                else{       boost::math::lognormal_distribution<double> LN_ep3(-6.2009, 0.3319); return pdf(LN_ep3, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c4ep1")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_ep1(-4.26891, 0.00344986); return LN_ep1(rng);}
                else{       boost::math::lognormal_distribution<double> LN_ep1(-4.26891, 0.00344986); return pdf(LN_ep1, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c4ep2")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_ep2(-5.02743, 0.144197); return LN_ep2(rng);}
                else{       boost::math::lognormal_distribution<double> LN_ep2(-5.02743, 0.144197); return pdf(LN_ep2, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c4ep3")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_ep3(-5.78595, 0.270815); return LN_ep3(rng);}
                else{       boost::math::lognormal_distribution<double> LN_ep3(-5.78595, 0.270815); return pdf(LN_ep3, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c4ep4")
            {
                if (sample){boost::random::lognormal_distribution<double> LN_ep4(-6.54447, 0.376344); return LN_ep4(rng);}
                else{       boost::math::lognormal_distribution<double> LN_ep4(-6.54447, 0.376344); return pdf(LN_ep4, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c5ep1")
            {
                if (sample){boost::random::normal_distribution<double> LN_ep1(-3.97331, 0.135741); return LN_ep1(rng);}
                else{       boost::math::normal_distribution<double> LN_ep1(-3.97331, 0.135741); return pdf(LN_ep1, pars.mapofParam.find(par)->second);}
            }
            else if (par == "c5ep2")
            {
                if (sample){boost::random::normal_distribution<double> LN_ep2(-0.17945, 0.0413348); return LN_ep2(rng);}
                else{       boost::math::normal_distribution<double> LN_ep2(-0.17945, 0.0413348); return pdf(LN_ep2, pars.mapofParam.find(par)->second);}
            }
            else if (par == "ep5")
            {
                if (sample){boost::random::gamma_distribution<double> G_ep5(35.0678,0.00000261628); return G_ep5(rng);}
                else{       boost::math::gamma_distribution<double> G_ep5(35.0678,0.00000261628); return pdf(G_ep5, pars.mapofParam.find(par)->second);}
            }
            else if (par == "ep6")
            {
                if (sample){boost::random::gamma_distribution<double> G_ep6(59.2461,0.00000228079); return G_ep6(rng);}
                else{       boost::math::gamma_distribution<double> G_ep6(59.2461,0.00000228079); return pdf(G_ep6, pars.mapofParam.find(par)->second);}
            }
            else if (par == "seed1")
            {
                if (sample){boost::random::uniform_real_distribution<double> U_s1(0, 1); return U_s1(rng);}
                else{       boost::math::uniform_distribution<double> U_s1(0, 1); return pdf(U_s1, pars.mapofParam.find(par)->second);}
            }
            else if (par == "psi")
            {
                if (sample){boost::random::uniform_real_distribution<double> U_s1(0, 1); return U_s1(rng);}
                else{       boost::math::uniform_distribution<double> U_s1(0, 1); return pdf(U_s1, pars.mapofParam.find(par)->second);}
            }
            else
            {
                if (sample){boost::random::uniform_real_distribution<double> U_o(0.0, 1.0); return U_o(rng);}
                else{       boost::math::uniform_distribution<double> U_o(0.0, 1.0); return pdf(U_o, pars.mapofParam.find(par)->second);}
            }
        }
    };
    
    // FUNC: init_prior_sample ->  get the initial prior samples
    // ARG: pars -> structure of parameter values
    void init_prior_sample(param_state_t& pars)
    {
        for (int i = 0; i < pars.dim_cal; i++)
            pars.mapofParam[pars.paramFit[i]] = pars.lprior_dist(pars.paramFit[i], pars, true);
    }
    
    // FUNC: init_prior_sample ->  get the log prior for metropolis hastings
    // ARG: pars -> structure of parameter values
    double prior_dist_get(param_state_t& pars)
    {
        double pp = 0;
        for (int i = 0; i < pars.dim_cal; i++)
            pp += log(pars.lprior_dist(pars.paramFit[i], pars, false));

        return pp;
    }
    
    // FUNC: init_prior_sample ->  exact a vector of the prior values from the pars structure
    // ARG: pars -> structure of parameter values
    VectorXd parameterfitonly(param_state_t& pars)
    {
        VectorXd fitonly(pars.dim_cal);
        for (int i = 0; i < pars.dim_cal; i++)
            fitonly(i) = pars.mapofParam[pars.paramFit[i]];
        
        return fitonly;
    }
    
    // FUNC: update -> update the pars structure with the new point
    // ARG: pars -> structure of parameter values
    // ARG: prop_par -> proposed point from the metropolis hastings

    void update(param_state_t& pars, VectorXd prop_par)
    {
        for (int i = 0; i < pars.dim_cal; i++)
            pars.mapofParam[pars.paramFit[i]] = prop_par(i);
    }
}

/// Describe method for obtaining the likelihood given the parameters

namespace get_ll
{
    /// Describe method for obtaining the likelihood given the parameters
    /** ///////////////////////////////////////////////////////////////////////////////////
     ////////////////////////////// 1. Import stuff //////////////////////////////
     //////////////////////////////////////////////////////////////////////////////////// **/
    string const filepath = "/Users/davidhodgson/documents/research/rsv/packages/rsv_trans_model/"; // change to a relative path 
    
    string const din = filepath + "inputs/";
    // Demographic files
    const string in_demo = din + "demo/";
    string FileDemo = in_demo + "pop_age_17_conv.txt"; //The source of the demography of the UK age groups
    string FileBDrate = in_demo + "bd_rate.txt"; //The source of the demography of the UK age groups
    string FileAgerate = in_demo +  "age_rate.txt"; //The source of the demography of the UK age groups
    string FilePropmat = in_demo + "prop_par.txt"; //The source of the demography of the UK age groups
    string FileMatdist = in_demo + "prop_mat.txt"; //The source of the demography of the UK age groups
    string FilePropRisk = in_demo + "proprisk.txt"; //The source of the demography of the UK age groups
    string FilePalAge = in_demo +  "palpro_age.txt";

    // Surveillance data
    const string in_rdms = din + "rdms/";
    string FileRDMSRSVpos = in_rdms + "RSV_RSDM_pos_trim.txt";
    string FileRDMSScalarTest = in_rdms +  "scalar.txt";
    
    // Contact matrices
    const string in_cnt_mat = din + "cnt_mat/";
#ifdef Minutes
    string filecontact_phys = in_cnt_mat + "cntMinPA.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_phys_home = in_cnt_mat + "cntMinPAH.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_pwp_p = in_cnt_mat + "cntMinPpwp.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_pwn_p = in_cnt_mat + "cntMinPpwn.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_nwp_p = in_cnt_mat + "cntMinPnwp.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_nwn_p = in_cnt_mat + "cntMinPnwn.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_conv = in_cnt_mat + "cntMinCA.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_conv_home = in_cnt_mat + "cntMinCAH.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_pwp_c = in_cnt_mat + "cntMinCpwp.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_pwn_c = in_cnt_mat + "cntMinCpwn.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_nwp_c = in_cnt_mat + "cntMinCnwp.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_nwn_c = in_cnt_mat + "cntMinCnwn.txt"; //The source of the contact data from AJ STudy (<1yrs)
#else
    string filecontact_phys = in_cnt_mat + "cntPA.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_phys_home = in_cnt_mat + "cntPAH.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_pwp_p = in_cnt_mat + "cntPpwp.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_pwn_p = in_cnt_mat + "cntPpwn.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_nwp_p = in_cnt_mat + "cntPnwp.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_nwn_p = in_cnt_mat + "cntPnwn.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_conv = in_cnt_mat + "cntCA.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_conv_home = in_cnt_mat + "cntCAH.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_pwp_c = in_cnt_mat + "cntCpwp.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_pwn_c = in_cnt_mat + "cntCpwn.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_nwp_c = in_cnt_mat + "cntCnwp.txt"; //The source of the contact data from AJ STudy (<1yrs)
    string filecontact_nwn_c = in_cnt_mat + "cntCnwn.txt"; //The source of the contact data from AJ STudy (<1yrs)
#endif
    
    
    //Out Files //
    string const dout = filepath + "outputs/";
    string FilePar = dout + "data/par.txt";      //The list of parameters used in the solution
    string FileData = dout + "data/pos_samples.txt";      //Output of the reconstructed calibration data
    
    
    // Posterior samples
#ifdef Minutes
    string model_t = dout + "cal/m2d5min/"; //
#else
    string model_t = dout + "cal/m2d5cnt/"; //
#endif
    
    string fileposterior_in = model_t + "posterior_A.txt";
    //string fileposterior_in = dout + calout_name + "mc/posterior_A_edit.txt";
    string fileannualinc = dout + "soln/annual_sample_Z.txt";

    
    // functions to help with the import
    namespace import
    {
        
        num_vec get_1d(string filename, int x)
        {
            ifstream temp(filename.c_str());
            num_vec cnt_1D;
            for (int j = 0; j < x; j++)
            {
                double tmp;
                temp >> tmp;
                cnt_1D.push_back(tmp);
            }
            
            return cnt_1D;
        }
        
        
        vector2D get_2d(string filename, int x, int y)
        {
            ifstream temp(filename.c_str());
            vector2D cnt_2D;
            num_vec cnt_1D;
            for (int i = 0; i < x; i++)
            {
                cnt_1D.clear();
                for (int j = 0; j < y; j++)
                {
                    double tmp;
                    temp >> tmp;
                    cnt_1D.push_back(tmp);
                }
                cnt_2D.push_back(cnt_1D);
            }
            return cnt_2D;
        }
        
        vector3D get_3d(string filename, int x, int y, int z)
        {
            ifstream temp(filename.c_str());
            vector3D cnt_3D;
            vector2D cnt_2D;
            num_vec cnt_1D;
            for (int i = 0; i < x; i++)
            {
                cnt_2D.clear();
                for (int j = 0; j < y; j++)
                {
                    cnt_1D.clear();
                    for (int k = 0; k < z; k++)
                    {
                        double tmp;
                        temp >> tmp;
                        cnt_1D.push_back(tmp);
                    }
                    cnt_2D.push_back(cnt_1D);
                }
                cnt_3D.push_back(cnt_2D);
            }
            return cnt_3D;
        }
        
        num_vec get_p(num_vec popsize, double tot_pop)
        {
            num_vec pop_p;
            for (int a = 0; a < NoAgeG; a++)
            {
                pop_p.push_back(popsize[a]/tot_pop);
            }
            return pop_p;
        }
        
        num_vec get_pop_mat(num_vec popsize, num_vec p_mat, bool ind)
        {
            num_vec popsize_alt;
            if (ind)
            {
                for (int a = 0; a < NoAgeG; a++)
                    popsize_alt.push_back(popsize[a]*p_mat[a]);
            }
            else
            {
                for (int a = 0; a < NoAgeG; a++)
                    popsize_alt.push_back(popsize[a]*(1-p_mat[a]));
            }
            return popsize_alt;
        }
        
        // All the data needed for the model
        vector2D const d_w = get_2d(FileRDMSRSVpos, 364, 25);
        vector2D const Scalar_test = get_2d(FileRDMSScalarTest, 25, 7);
        
        vector2D const cnt_matrix_p = get_2d(filecontact_phys, 25, 25);
        vector2D const cnt_matrix_p_h = get_2d(filecontact_phys_home, 25, 25);
        vector2D const pwp_p = get_2d(filecontact_pwp_p, 25, 25);
        vector2D const pwn_p = get_2d(filecontact_pwn_p, 25, 25);
        vector2D const nwp_p = get_2d(filecontact_nwp_p, 25, 25);
        vector2D const nwn_p = get_2d(filecontact_nwn_p, 25, 25);
        
        vector2D const cnt_matrix_c = get_2d(filecontact_conv, 25, 25);
        vector2D const cnt_matrix_c_h = get_2d(filecontact_conv_home, 25, 25);
        vector2D const pwp_c = get_2d(filecontact_pwp_c, 25, 25);
        vector2D const pwn_c = get_2d(filecontact_pwn_c, 25, 25);
        vector2D const nwp_c = get_2d(filecontact_nwp_c, 25, 25);
        vector2D const nwn_c = get_2d(filecontact_nwn_c, 25, 25);
        
        vector2D const risk = get_2d(FilePropRisk, 3, 25);
        num_vec const pVHR = risk[0];
        num_vec const pHR = risk[1];
        num_vec const pLR = risk[2];
        
        num_vec const eta = get_1d(FileAgerate, 26);
        num_vec const popsize = get_1d(FileDemo, 25);
        num_vec const br = get_1d(FileBDrate, 2);
        num_vec const p_mat = get_1d(FilePropmat, 25);
        num_vec const u_p = get_1d(FileMatdist, 3);
        double const u18p = u_p[0];  double const u19p = u_p[1]; double const u20p = u_p[2];
        num_vec popsize_m = get_pop_mat(popsize, p_mat, true);
        num_vec popsize_n = get_pop_mat(popsize, p_mat, false);

        const double muB = br[0];
        
        double tot_pop = accumulate(popsize.begin(),popsize.end(),0);
        num_vec pop_p = get_p(popsize, tot_pop);

        vector2D ann_inc = get_2d(fileannualinc, 25,  2);
        
                
    }
    using namespace import;

    // FUNC: prop_init_ex- > function to get the initial proportion of each age group in each exposure group
    num_vec prop_init_ex(double l, double a1, double a2)
    {
        num_vec prop(4);
        prop[0] = abs(poisson_cdf(l,a2,0)-poisson_cdf(l,a1,0))/((a2-a1)*l);
        prop[1] = abs(poisson_cdf(l,a2,1)-poisson_cdf(l,a1,1))/((a2-a1)*l);
        prop[2] = abs(poisson_cdf(l,a2,2)-poisson_cdf(l,a1,2))/((a2-a1)*l);
        prop[3] = 1 - (prop[2]+prop[1]+prop[0]);
        return prop;
    }
    
    // FUNC: initial_M - > function to get the initial proportion of each age group in each MATERNAL group
    num_vec initial_M(const param::param_state_t& pars)
    {
        double xi = 1.0/(pars.mapofParam.find("xi")->second);
        boost::math::exponential_distribution <> exp ( xi );
        
        num_vec init_con;
        for (int i = 0; i< NoAgeG; i++)
        {
            double init_con_temp = (cdf(exp, 365*ageGroupP[i+1])-cdf(exp, 365*ageGroupP[i]))/((365*ageGroupP[i+1]-365*ageGroupP[i])*xi);
            init_con.push_back(init_con_temp*popsize[i]);
        }
        return init_con;
    }
    
    // FUNC: init_cond - > function to get the initial conditions for the epidemiological state variables

    num_vec init_cond(const param::param_state_t& pars)
    {
        
        num_vec age_size_v = popsize;
        num_vec comp_size_init;
        num_vec init_con_m = initial_M(pars);
        VectorXd pA(NoAgeG);

        for (int a = 0; a < NoAgeG; a++)
            pA(a) = pars.mapofParam.find(pars.pA_vec[a])->second;
        
        double init_inf = pars.mapofParam.find("seed1")->second;
        double init_rec = pars.mapofParam.find("seed2")->second;
        double init_exp_time = 0.5; //1/(assumed duration of infection cycle)
        
        
        double pI1; double pI2; double pI3; double pI4;
        double a1, a2;
        double tot = 0;
        double si = 1.0/(pars.mapofParam.find("si")->second);
        double g0 = 1.0/(pars.mapofParam.find("g0")->second);
        double g1 = 1.0/((pars.mapofParam.find("g0")->second)*(pars.mapofParam.find("g1")->second));
        double g2 = 1.0/((pars.mapofParam.find("g0")->second)*(pars.mapofParam.find("g1")->second)*(pars.mapofParam.find("g2")->second));
        double d1 = pars.mapofParam.find("d1")->second;
        double d2 = (pars.mapofParam.find("d1")->second)*(pars.mapofParam.find("d2")->second);
        double d3 = (pars.mapofParam.find("d1")->second)*(pars.mapofParam.find("d2")->second)*(pars.mapofParam.find("d3")->second);
        
        /** Get ratio of number of each person in each exposure group per age group **/
        for(int i = 0; i < NoAgeG; i++)
        {
            if (i < 24)
            {
                a1 = ageGroup[i]; a2 = ageGroup[i+1];
            }
            
            else
            {
                a1 = ageGroup[i]; a2 = 90;
            }
            num_vec prop_init = prop_init_ex(init_exp_time, a1, a2);
            pI1 = prop_init[0]; pI2 = prop_init[1]; pI3 = prop_init[2]; pI4 = prop_init[3];
            double age_size = age_size_v[i] - init_con_m[i];
        

            comp_size_init.push_back(init_con_m[i]); // Number  in M group
            
            comp_size_init.push_back(pI1*age_size*(1.0 - init_inf)*(1.0-init_rec));      //Sus
            comp_size_init.push_back(pI1*age_size*init_inf*si/(si+g0));       //Exp
            comp_size_init.push_back(pI1*age_size*init_inf*g0/(si+g0)*pA(i));
            comp_size_init.push_back(pI1*age_size*init_inf*g0/(si+g0)*(1-pA(i)));   //Inf S
            comp_size_init.push_back(pI1*age_size*(1.0 - init_inf)*init_rec);      //Rec
            
            comp_size_init.push_back(pI2*age_size*(1.0 - d1*init_inf)*(1.0-init_rec));      //Sus
            comp_size_init.push_back(pI2*age_size*d1*init_inf*si/(si+g1));       //Exp
            comp_size_init.push_back(pI2*age_size*d1*init_inf*g1/(si+g1)*pA(i));
            comp_size_init.push_back(pI2*age_size*d1*init_inf*g1/(si+g1)*(1-pA(i)));   //Inf S
            comp_size_init.push_back(pI2*age_size*(1.0 - d1*init_inf)*init_rec);      //Rec
            
            comp_size_init.push_back(pI3*age_size*(1.0 - d2*init_inf)*(1.0-init_rec));      //S/Sus
            comp_size_init.push_back(pI3*age_size*d2*init_inf*si/(si+g2));       //Exp
            comp_size_init.push_back(pI3*age_size*d2*init_inf*g2/(si+g2)*pA(i));
            comp_size_init.push_back(pI3*age_size*d2*init_inf*g2/(si+g2)*(1-pA(i)));   //Inf S
            comp_size_init.push_back(pI3*age_size*(1.0 - d2*init_inf)*init_rec);      //ec
            
            comp_size_init.push_back(pI4*age_size*(1.0 - d3*init_inf)*(1.0-init_rec));        //Sus
            comp_size_init.push_back(pI4*age_size*d3*init_inf*si/(si+g2));       //Exp
            comp_size_init.push_back(pI4*age_size*d3*init_inf*g2/(si+g2)*pA(i));
            comp_size_init.push_back(pI4*age_size*d3*init_inf*g2/(si+g2)*(1-pA(i)));   //Inf S
            comp_size_init.push_back(pI4*age_size*(1.0 - d3*init_inf)*init_rec);     //Rec
            
            comp_size_init.push_back(0.0);     //Rec
            comp_size_init.push_back(0.0);     //Rec
            
        }
        for (int i = 0; i < 23*25; i++)
            tot += comp_size_init[i];
        
        // cout << "comp_size_init: " << tot << endl;
        return comp_size_init;
    }
    

    struct ODE_dynamics
    {
        double t_start, t_burn, t_end, dt;
        bool get_incidence;
        
        ODE_dynamics(const double t_start_i, const double t_burn_i, const double t_end_i, const double dt_i): t_start(t_start_i), t_burn(t_burn_i), t_end(t_end_i) , dt(dt_i) {}
    };
    
    // Class dscribing the system of ODEs for calibrating
    class ODE_desc
    {
        param::param_state_t  pars;
        
    public:
        double pA[NoAgeG];
        double dsI[NoAgeG];
        
        double p_vul;
        
        /** Single infection parameters **/
        double xi = 1.0/(pars.mapofParam.find("xi")->second);
        double si = 1.0/(pars.mapofParam.find("si")->second);
        double ga0 = 1.0/((pars.mapofParam.find("g0")->second));
        double ga1 = 1.0/((pars.mapofParam.find("g0")->second)*(pars.mapofParam.find("g1")->second));
        double ga2 = 1.0/((pars.mapofParam.find("g0")->second)*(pars.mapofParam.find("g1")->second)*(pars.mapofParam.find("g2")->second));
        double ga3 = ga2;
        double om = 1.0/(pars.mapofParam.find("om")->second);
        /** Long-term infection parameters **/
        //double rho = pars.mapofParam.find("rho")->second;
        double rho = 1.0;
        double alpha_i = pars.mapofParam.find("alpha_i")->second;
        double d1 = (pars.mapofParam.find("d1")->second);
        double d2 = (pars.mapofParam.find("d1")->second)*(pars.mapofParam.find("d2")->second);
        double d3 = (pars.mapofParam.find("d1")->second)*(pars.mapofParam.find("d2")->second)*(pars.mapofParam.find("d3")->second);
        double a1 = 1.0, a2 = 1.0, a3 = 1.0;
        /** Transmission  parameters **/
        double phi =(pars.mapofParam.find("phi")->second);
        double qp = (pars.mapofParam.find("qp")->second);
        double qc = (pars.mapofParam.find("qc")->second);
        double b1 = (pars.mapofParam.find("b1")->second);
        double psi = (pars.mapofParam.find("psi")->second);

        double M, S0, S1, S2, S3, E0, E1, E2, E3, A0, A1, A2, A3, I0, I1, I2, I3, R1, R2, R3, N;
        double I[NoAgeG];
        double N_tot[NoAgeG];
        double dx;
        double x_tot, I_temp;
        double beta;
        double t1;
        double Nall1;
        double Nall2;
        double Nall2_vec[NoAgeG];
        double sum_wcb;
        vector2D pI; num_vec pI_a;
        num_vec d, g;
        
        ODE_desc(const param::param_state_t& pars_t): pars(pars_t) {
            
            d = {1, d1, d2, d3};
            g = {ga0, ga1, ga2, ga3};

            for (int a = 0; a < NoAgeG; a++)
                pA[a] = pars.mapofParam.find(pars.pA_vec[a])->second;
            
            for (int a = 0; a < NoAgeG ; a++)
            {
                pI_a = prop_init_ex(0.5, ageGroupP[a], ageGroupP[a+1]);
                pI.push_back(pI_a);
            }
            p_vul = 0.0;
        }
        
    /*    double find_REff(ODE_desc &ODE_i)
        {
            MatrixXd T =  MatrixXd::Zero(NoAgeG*12,NoAgeG*12);
            MatrixXd Sigma =  MatrixXd::Zero(NoAgeG*12,NoAgeG*12);
            MatrixXd K =  MatrixXd::Zero(NoAgeG*12,NoAgeG*12);
            
            for (int a1 = 0; a1 < NoAgeG; a1++)
                for (int a2 = 0; a2 < NoAgeG; a2++)
                    for (int i1 = 0; i1 < 4; i1++)
                        for (int i2 = 0; i2 < 4; i2++){
                            T(a1*12 + i1*3, a2*12 + i2*3  + 1) = ODE_i.alpha_i*ODE_i.d[i2]*cnt_matrix[a1][a2]*pI[a1][i1]*pI[a2][i2];
                            T(a1*12 + i1*3, a2*12 + i2*3  + 2) = ODE_i.d[i2]*cnt_matrix[a1][a2]*pI[a2][i1]*pI[a2][i2];
                        }


            for (int a1 = 0; a1 < NoAgeG; a1++)
            {
                for (int i1 = 0; i1 < 4; i1++)
                {
                    Sigma(a1*12 + i1*3, a1*12 + i1*3 + 0) = (si - eta[a1+1])*pI[a1][i1];
                    Sigma(a1*12 + i1*3 + 1, a1*12 + i1*3 + 0) = pA[a1]*si*pI[a1][i1];
                    Sigma(a1*12 + i1*3 + 2, a1*12 + i1*3 + 0) = (1-pA[a1])*si*pI[a1][i1];
                    Sigma(a1*12 + i1*3 + 1, a1*12 + i1*3 + 1) = (ODE_i.g[i1]*rho - eta[a1+1])*pI[a1][i1];
                    Sigma(a1*12 + i1*3 + 2, a1*12 + i1*3 + 2) = (ODE_i.g[i1] - eta[a1+1])*pI[a1][i1];
                }
            }

            for (int a1 = 0; a1 < NoAgeG-1; a1++)
                for (int i1 = 0; i1 < 4; i1++)
                {
                    Sigma((a1+1)*12 + i1*3, a1*12 + i1*3 + 0) = eta[a1+1]*pI[a1][i1];
                    Sigma((a1+1)*12 + i1*3 + 1, a1*12 + i1*3 + 1) = eta[a1+1]*pI[a1][i1];
                    Sigma((a1+1)*12 + i1*3 + 2, a1*12 + i1*3 + 2) = eta[a1+1]*pI[a1][i1];
                }
            
            MatrixXd SigmaInv = Sigma.householderQr().solve(MatrixXd::Identity(NoAgeG*12, NoAgeG*12));
            K = -T*SigmaInv;
            double evs = K.eigenvalues().real().maxCoeff();
            double find_REff = ODE_i.b0*evs;
            return find_REff;
        }
        
        double find_R0(ODE_desc &ODE_i)
        {
            MatrixXd T =  MatrixXd::Zero(NoAgeG*12,NoAgeG*12);
            MatrixXd Sigma =  MatrixXd::Zero(NoAgeG*12,NoAgeG*12);
            MatrixXd K =  MatrixXd::Zero(NoAgeG*12,NoAgeG*12);
            
            for (int a1 = 0; a1 < NoAgeG; a1++)
                for (int a2 = 0; a2 < NoAgeG; a2++)
                    for (int i1 = 0; i1 < 4; i1++)
                        for (int i2 = 0; i2 < 4; i2++){
                            T(a1*12 + i1*3, a2*12 + i2*3  + 1) = ODE_i.alpha_i*cnt_matrix[a1][a2]*pI[a1][i1]*pI[a2][i2];
                            T(a1*12 + i1*3, a2*12 + i2*3  + 2) = cnt_matrix[a1][a2]*pI[a2][i1]*pI[a2][i2];
                        }
            
            
            for (int a1 = 0; a1 < NoAgeG; a1++)
            {
                for (int i1 = 0; i1 < 4; i1++)
                {
                    Sigma(a1*12 + i1*3, a1*12 + i1*3 + 0) = (si - eta[a1+1])*pI[a1][i1];
                    Sigma(a1*12 + i1*3 + 1, a1*12 + i1*3 + 0) = pA[a1]*si*pI[a1][i1];
                    Sigma(a1*12 + i1*3 + 2, a1*12 + i1*3 + 0) = (1-pA[a1])*si*pI[a1][i1];
                    Sigma(a1*12 + i1*3 + 1, a1*12 + i1*3 + 1) = (ODE_i.g[i1]*rho - eta[a1+1])*pI[a1][i1];
                    Sigma(a1*12 + i1*3 + 2, a1*12 + i1*3 + 2) = (ODE_i.g[i1] - eta[a1+1])*pI[a1][i1];
                }
            }
            
            for (int a1 = 0; a1 < NoAgeG-1; a1++)
                for (int i1 = 0; i1 < 4; i1++)
                {
                    Sigma((a1+1)*12 + i1*3, a1*12 + i1*3 + 0) = eta[a1+1]*pI[a1][i1];
                    Sigma((a1+1)*12 + i1*3 + 1, a1*12 + i1*3 + 1) = eta[a1+1]*pI[a1][i1];
                    Sigma((a1+1)*12 + i1*3 + 2, a1*12 + i1*3 + 2) = eta[a1+1]*pI[a1][i1];
                }
            
            MatrixXd SigmaInv = Sigma.householderQr().solve(MatrixXd::Identity(NoAgeG*12, NoAgeG*12));
            K = -T*SigmaInv;
            double evs = K.eigenvalues().real().maxCoeff();
            double find_R0= ODE_i.b0*evs;
            return find_R0;
        }*/
        
        void operator() (  num_vec &x , num_vec &dxdt , const double  t  )
        {
            t1 = (int)t%365;
           // beta = (1 + b1*(1 + exp(-((t1/365 - phi))*((t1/365 - phi))/(2*psi*psi))));
            //beta = (1 + b1*cos((t/365.0-phi)*2*PI));
#ifdef DynamicMat
            num_vec num_vec_wcba;
            sum_wcb = 0.0;
            
            for(int j = 18; j < 21 ; j++)
            {
                double CB2_temp = (x[1 + j*23] + x[6 + j*23] + x[11 + j*23] + x[16 + j*23] + x[2 + j*23] + x[7 + j*23] + x[12 + j*23] + x[17 + j*23])/(double)popsize[j];
                num_vec_wcba.push_back(CB2_temp/3.0);
            }
            for(int i=0; i < boost::size(num_vec_wcba) ;i++){sum_wcb = sum_wcb + num_vec_wcba[i];}
            double p_vul = sum_wcb;
#else
            
#endif
            
            for(int j = 0; j < NoAgeG; j++)
            {
                Nall1 = 0;
                for (int k = 0; k < 21 ;k++)
                    Nall1 += x[k + 23*j];
                
                I_temp = 0;
                for (int k = 0; k < NoAgeG ; k++)
                {
                    I_temp += (x[3 + 23*k]*alpha_i+x[4+23*k]+a1*(x[8+23*k]*alpha_i+x[9+23*k])+a2*(x[13+ 23*k]*alpha_i+x[14+23*k])+a3*(x[18+23*k]*alpha_i+x[19+23*k]))*(qp*(cnt_matrix_p[j][k]+qc*cnt_matrix_c[j][k]))/(popsize[k]);
                }
                
                int pj = max(23*(j-1),0);
                //cout << pj << endl;
                double muBv;
                if (j == 0)
                    muBv = muB;
                else
                    muBv = 0;
                // These differentail equations appear at the start of the appendix
                dxdt[0 + 23*j] = (1.0-p_vul)*muBv - x[0 + 23*j]*xi                                   - x[0+23*j]*eta[j+1] + x[0+pj]*eta[j];
                
                dxdt[1 + 23*j] =  p_vul*muBv + x[0 + 23*j]*xi       - x[1 + 23*j]*I_temp*beta        - x[1+23*j]*eta[j+1] + x[1+pj]*eta[j];
                
                dxdt[2 + 23*j] = x[1 + 23*j]*I_temp*beta               - x[2+23*j]*si                - x[2+23*j]*eta[j+1] + x[2+pj]*eta[j];
                dxdt[3 + 23*j] = x[2 + 23*j]*si*pA[j]                - x[3+23*j]*ga0*rho             - x[3+23*j]*eta[j+1] + x[3+pj]*eta[j];
                dxdt[4 + 23*j] = x[2 + 23*j]*si*(1.0-pA[j])            - x[4+23*j]*ga0               - x[4+23*j]*eta[j+1] + x[4+pj]*eta[j];
                dxdt[5 + 23*j] = x[4 + 23*j]*ga0 + x[3+23*j]*ga0*rho        - x[5+23*j]*om           - x[5+23*j]*eta[j+1] + x[5+pj]*eta[j];
                
                dxdt[6 + 23*j] = x[5+23*j]*om                       - d1*x[6+23*j]*I_temp*beta       - x[6+23*j]*eta[j+1] + x[6+pj]*eta[j];
                dxdt[7 + 23*j] = d1*x[6+23*j]*I_temp*beta            - x[7+23*j]*si                  - x[7+23*j]*eta[j+1] + x[7+pj]*eta[j];
                dxdt[8 + 23*j] = x[7+23*j]*si*pA[j]                - x[8+23*j]*ga1*rho               - x[8+23*j]*eta[j+1] + x[8+pj]*eta[j];
                dxdt[9 + 23*j] = x[7+23*j]*si*(1.0-pA[j])            - x[9+23*j]*ga1                 - x[9+23*j]*eta[j+1] + x[9+pj]*eta[j];
                dxdt[10 + 23*j] = x[9+23*j]*ga1 + x[8+23*j]*ga1*rho       - x[10+23*j]*om            - x[10+23*j]*eta[j+1] + x[10+pj]*eta[j];
                
                dxdt[11 + 23*j] = x[10+23*j]*om                     - d2*x[11+23*j]*I_temp*beta      - x[11+23*j]*eta[j+1] + x[11+pj]*eta[j];
                dxdt[12 + 23*j] = d2*x[11+23*j]*I_temp*beta           - x[12+23*j]*si                - x[12+23*j]*eta[j+1] + x[12+pj]*eta[j];
                dxdt[13 + 23*j] = x[12+23*j]*si*pA[j]              - x[13+23*j]*ga2*rho              - x[13+23*j]*eta[j+1] + x[13+pj]*eta[j];
                dxdt[14 + 23*j] = x[12+23*j]*si*(1.0-pA[j])           - x[14+23*j]*ga2               - x[14+23*j]*eta[j+1] + x[14+pj]*eta[j];
                dxdt[15 + 23*j] = x[14+23*j]*ga2 + x[13+23*j]*ga2*rho       - x[15+23*j]*om          - x[15+23*j]*eta[j+1] + x[15+pj]*eta[j];
                
                dxdt[16 + 23*j] = x[15+23*j]*om + x[20+23*j]*om     - d3*x[16+23*j]*I_temp*beta      - x[16+23*j]*eta[j+1] + x[16+pj]*eta[j];
                dxdt[17 + 23*j] = d3*x[16+23*j]*I_temp*beta           - x[17+23*j]*si                - x[17+23*j]*eta[j+1] + x[17+pj]*eta[j];
                dxdt[18 + 23*j] = x[17+23*j]*si*pA[j]               - x[18+23*j]*ga3*rho             - x[18+23*j]*eta[j+1] + x[18+pj]*eta[j];
                dxdt[19 + 23*j] = x[17+23*j]*si*(1.0-pA[j])           - x[19+23*j]*ga3               - x[19+23*j]*eta[j+1] + x[19+pj]*eta[j];
                dxdt[20 + 23*j] = x[19+23*j]*ga3 + x[18+23*j]*ga3*rho       - x[20+23*j]*om          - x[20+23*j]*eta[j+1] + x[20+pj]*eta[j];
                
                dxdt[21 + 23*j] = 0;
                
                dxdt[22 + 23*j] = si*(x[2 + 23*j] + x[7 + 23*j] + x[12 + 23*j] + x[17 + 23*j]);
            }
        }
    };

    // llikelihood_func -> function to calculate the likelihood (at a time t)
    double llikelihood_func(VectorXd& inc, const param::param_state_t& pars, int t_w, VectorXd& ep_t){
        
        double ll = 0;
        for (int a = 0; a < NoAgeG; a++)
        {
            double N = inc(a);
            double x = d_w[t_w][a]/Scalar_test[a][t_w/52];
            
            if (x > N)
                return log(0);
            
            // log binomial -> use stirlings approximation for log (x take N)
            ll += stirl(N) - stirl(x) - stirl(N-x) + x*log(ep_t(a)) + (N-x)*log(1-ep_t(a));
        }
        return ll;
    }
    
    
   // get_lprior -> function for log prior
    double get_lprior(double ll, param::param_state_t &pars)
    {
        if (std::isinf(ll) || std::isnan(ll))
            return 0;
        else
           return param::prior_dist_get(pars);
    }
    
    // check_stability -> check the system of ODEs is working (no negative values)
    inline double check_stability(num_vec x, double t)
    {
        long double X_w = 0;
        for (int j = 0; j < NoAgeG*23; j++)
        {
            if (x[j] < 0)
            {
                cout << "Negative value at pos: " << j << ". At time: " << t <<  '\n';
                return log(0);
            }
            else
                X_w += x[j];
        }
        
        if (std::isinf(X_w) || std::isnan(X_w))
            return log(0);
        else
            return 1.0;
    }
    
    
    inline double check_incidence(VectorXd inc_tot)
    {
        double pl = 0;
        for (int a = 0; a < 25; a++)
        {
            if (inc_tot(a) > popsize[a]*0.8)
            {
                return log(0);
            }
        }

        
        return pl;
    }
    
    // get_likelihood -> Functionto get the likehood over all time for a set of parameter values pars
    double get_likelihood(param::param_state_t& pars) //(7 years)
    {

        ODE_dynamics ode_par(0, 365, 365*3+1, 1.0);
        
        EulerT<state_t> integrator;
        SystemT<state_t, system_t> System;
        asc::Recorder recorder;
        
        VectorXd inc(NoAgeG);
        VectorXd inc_tot(NoAgeG);
        inc_tot = VectorXd::Zero(NoAgeG);
        
        ODE_desc ODE_desc_inst(pars);
        num_vec x0 = get_ll::init_cond(pars);

        double t = ode_par.t_start;
        int t_d = 0;
        int t_w = 0;
        double ll = 0;
        double val;
        
        VectorXd ep_t(NoAgeG);
#ifdef Cal1
        for (int a = 0; a < NoAgeG; a++)
            ep_t(a) = pars.mapofParam.find(pars.c1m[a])->second;
#endif
#ifdef Cal2
        for (int a = 0; a < NoAgeG; a++)
            ep_t(a) = pars.mapofParam.find(pars.c2m[a])->second;
#endif
#ifdef Cal3
        for (int a = 0; a < NoAgeG; a++)
            ep_t(a) = pars.mapofParam.find(pars.c3m[a])->second;
#endif
#ifdef Cal4
        for (int a = 0; a < NoAgeG; a++)
            ep_t(a) = pars.mapofParam.find(pars.c4m[a])->second;
#endif
#ifdef Cal5
        for (int a = 0; a < 16; a++)
            ep_t(a) = exp(pars.mapofParam.find(pars.c5mA[a])->second + a*(pars.mapofParam.find(pars.c5mB[a])->second));
        
        for (int a = 16; a < 25; a++)
        {
            ep_t(a) = pars.mapofParam.find(pars.c5mA[a])->second;
        }
#endif
        
        while (t < ode_par.t_end)
        {
            // Burn in section if required
            integrator(ODE_desc_inst, x0, t, 1.0);
            
            val = check_stability(x0, t);
            if (val < 0.0)
            {
                return log(0);
            }
            if (t > ode_par.t_burn)
            {
                if (t_d == 0)
                {
                    for (int a = 0; a < NoAgeG; a++)
                        x0[22 + 23*a] = 0.0; //Incidence at t_d = 0;
                }
                if (t_d%7 == 0 && t_d > 0)
                {
                    for (int a = 0; a < NoAgeG; a++)
                    {
                        inc(a) = x0[22 + 23*a]; //Incidence at t_d = 7;
                        inc_tot(a) += inc(a);
                        x0[22 + 23*a] = 0.0;
                    }
                    
                    ll += llikelihood_func(inc, pars, t_w, ep_t);
                    if (std::isinf(ll))
                    {
                        return log(0);
                    }
                    t_w ++;
                }
                if (t_d%365 == 0 && t_d > 0)
                {
                    ll += check_incidence(inc_tot);
                    if (std::isinf(ll))
                    {
                        return log(0);
                    }
                    inc_tot = VectorXd::Zero(NoAgeG);
                }
                t_d++;
            }
        }
        return ll;
    }
    
    //  I think everything from here is redundant.... 
    num_vec init_cond_par(const param::param_state_t& pars)
    {
        
        num_vec age_size_v = get_ll::popsize;
        num_vec comp_size_init;
        num_vec init_con_m = get_ll::initial_M(pars);
        VectorXd pA(NoAgeG);
        
        for (int a = 0; a < NoAgeG; a++)
            pA(a) = pars.mapofParam.find(pars.pA_vec[a])->second;
        
        double init_inf = pars.mapofParam.find("seed1")->second;
        double init_rec = pars.mapofParam.find("seed2")->second;
        double init_exp_time = 0.5; //1/(assumed duration of infection cycle)
        
        
        double pI1; double pI2; double pI3; double pI4;
        double a1, a2;
        double tot = 0;
        /** Get ratio of number of each person in each exposure group per age group **/
        for(int i = 0; i < NoAgeG; i++)
        {
            double age_size, matp_size;
            if (i < 24)
            {
                a1 = ageGroup[i]; a2 = ageGroup[i+1];
            }
            
            else
            {
                a1 = ageGroup[i]; a2 = 90;
            }
            num_vec prop_init = get_ll::prop_init_ex(init_exp_time, a1, a2);
            pI1 = prop_init[0]; pI2 = prop_init[1]; pI3 = prop_init[2]; pI4 = prop_init[3];
            
            double si = 1.0/(pars.mapofParam.find("si")->second);
            double g0 = 1.0/(pars.mapofParam.find("g0")->second);
            double g1 = 1.0/((pars.mapofParam.find("g0")->second)*(pars.mapofParam.find("g1")->second));
            double g2 = 1.0/((pars.mapofParam.find("g0")->second)*(pars.mapofParam.find("g1")->second)*(pars.mapofParam.find("g2")->second));
            double d1 = pars.mapofParam.find("d1")->second;
            double d2 = (pars.mapofParam.find("d1")->second)*(pars.mapofParam.find("d2")->second);
            double d3 = (pars.mapofParam.find("d1")->second)*(pars.mapofParam.find("d2")->second)*(pars.mapofParam.find("d3")->second);
                    
            for (int s = 0; s < 2; s++)
            {
                if (s==0)
                {
                    age_size = p_mat[i]*(age_size_v[i] - init_con_m[i]);
                    matp_size = p_mat[i]*init_con_m[i];
                }
                else if (s==1)
                {
                    age_size = (1-p_mat[i])*(age_size_v[i] - init_con_m[i]);
                    matp_size = (1-p_mat[i])*init_con_m[i];
                }
                
                comp_size_init.push_back(matp_size); // Number  in M group
                
                comp_size_init.push_back(pI1*age_size*(1.0 - init_inf)*(1.0-init_rec));      //Sus
                comp_size_init.push_back(pI1*age_size*init_inf*si/(si+g0));       //Exp
                comp_size_init.push_back(pI1*age_size*init_inf*g0/(si+g0)*pA(i));
                comp_size_init.push_back(pI1*age_size*init_inf*g0/(si+g0)*(1-pA(i)));   //Inf S
                comp_size_init.push_back(pI1*age_size*(1.0 - init_inf)*init_rec);      //Rec
                
                comp_size_init.push_back(pI2*age_size*(1.0 - d1*init_inf)*(1.0-init_rec));      //Sus
                comp_size_init.push_back(pI2*age_size*d1*init_inf*si/(si+g1));       //Exp
                comp_size_init.push_back(pI2*age_size*d1*init_inf*g1/(si+g1)*pA(i));
                comp_size_init.push_back(pI2*age_size*d1*init_inf*g1/(si+g1)*(1-pA(i)));   //Inf S
                comp_size_init.push_back(pI2*age_size*(1.0 - d1*init_inf)*init_rec);      //Rec
                
                comp_size_init.push_back(pI3*age_size*(1.0 - d2*init_inf)*(1.0-init_rec));      //S/Sus
                comp_size_init.push_back(pI3*age_size*d2*init_inf*si/(si+g2));       //Exp
                comp_size_init.push_back(pI3*age_size*d2*init_inf*g2/(si+g2)*pA(i));
                comp_size_init.push_back(pI3*age_size*d2*init_inf*g2/(si+g2)*(1-pA(i)));   //Inf S
                comp_size_init.push_back(pI3*age_size*(1.0 - d2*init_inf)*init_rec);      //ec
                
                comp_size_init.push_back(pI4*age_size*(1.0 - d3*init_inf)*(1.0-init_rec));        //Sus
                comp_size_init.push_back(pI4*age_size*d3*init_inf*si/(si+g2));       //Exp
                comp_size_init.push_back(pI4*age_size*d3*init_inf*g2/(si+g2)*pA(i));
                comp_size_init.push_back(pI4*age_size*d3*init_inf*g2/(si+g2)*(1-pA(i)));   //Inf S
                comp_size_init.push_back(pI4*age_size*(1.0 - d3*init_inf)*init_rec);     //Rec
                
                comp_size_init.push_back(0.0);     //Rec
                comp_size_init.push_back(0.0);     //Rec
            }
        }
        //for (int i = 0; i < 23*25; i++)
        //tot += comp_size_init[i];
        
        // cout << "comp_size_init: " << tot << endl;
        return comp_size_init;
    }
    
    
    class ODE_desc_par
    {
        param::param_state_t  pars;
        
    public:
        double pA[NoAgeG];
        double dsI[NoAgeG];
        
        double p_vul;
        
        ODE_desc_par(const param::param_state_t& pars_t): pars(pars_t) {
            
            for (int a = 0; a < NoAgeG; a++)
                pA[a] = pars.mapofParam.find(pars.pA_vec[a])->second;
            
            p_vul = 0.0;
        }
        
        /** Single infection parameters **/
        double xi = 1.0/(pars.mapofParam.find("xi")->second);
        double si = 1.0/(pars.mapofParam.find("si")->second);
        double ga0 = 1.0/((pars.mapofParam.find("g0")->second));
        double ga1 = 1.0/((pars.mapofParam.find("g0")->second)*(pars.mapofParam.find("g1")->second));
        double ga2 = 1.0/((pars.mapofParam.find("g0")->second)*(pars.mapofParam.find("g1")->second)*(pars.mapofParam.find("g2")->second));
        double ga3 = ga2;
        double om = 1.0/(pars.mapofParam.find("om")->second);
        /** Long-term infection parameters **/
       // double rho = pars.mapofParam.find("rho")->second;
        double rho = 1.0;
        double alpha_i = pars.mapofParam.find("alpha_i")->second;
        double d1 = (pars.mapofParam.find("d1")->second);
        double d2 = (pars.mapofParam.find("d1")->second)*(pars.mapofParam.find("d2")->second);
        double d3 = (pars.mapofParam.find("d1")->second)*(pars.mapofParam.find("d2")->second)*(pars.mapofParam.find("d3")->second);
        double a1 = 1.0, a2 = 1.0, a3 = 1.0;
        /** Transmission  parameters **/
        double phi =(pars.mapofParam.find("phi")->second);
        double qp = (pars.mapofParam.find("qp")->second);
        double qc = (pars.mapofParam.find("qc")->second);
        double b1 = (pars.mapofParam.find("b1")->second);
        double psi = (pars.mapofParam.find("psi")->second);
        
        
        double M, S0, S1, S2, S3, E0, E1, E2, E3, A0, A1, A2, A3, I0, I1, I2, I3, R0, R1, R2, R3, N;
        double I[NoAgeG];
        double N_tot[NoAgeG];
        double dx;
        double x_tot, I_temp, I_temp_n, I_temp_m;
        double beta;
        double t1;
        double Nall1;
        double Nall2;
        double Nall2_vec[NoAgeG];
        double sum_wcb;
        long double prop_mg, pop_n, pop_m, catsize;
        int ag, mg;
        /** Description of ODE equations **/
        void operator() (  num_vec &x , num_vec &dxdt , const double  t  )
        {
            //beta = (1 + b1*cos((t/365.0-phi)*2*PI));
            t1 = (int)t%365;
            beta = (1 + b1*(1 + exp(-((t1/365 - phi))*((t1/365 - phi))/(2*psi*psi))));
            ag = 46; mg = 23;
#ifdef DynamicMat
            num_vec num_vec_wcba;
            sum_wcb = 0.0;
            
            for(int j = 18; j < 21 ; j++)
            {
                double CB2_temp = (x[1+mg+j*ag] + x[6+mg+j*ag] + x[11+mg+j*ag] + x[16+mg+j*ag] + x[2+mg+j*ag] + x[7+mg+j*ag] + x[12+mg+j*ag] + x[17+mg+j*ag])/((double)get_ll::popsize_n[j]);
                num_vec_wcba.push_back(CB2_temp*(1-p_mat[j])/(3-p_mat[18]-p_mat[19]-p_mat[20]));
            }
            for(int i=0; i < boost::size(num_vec_wcba) ;i++){sum_wcb = sum_wcb + num_vec_wcba[i];}
            double p_vul = sum_wcb;
#else
            
#endif
            for(int j = 0; j < NoAgeG; j++)
            {
                I_temp_n = I_temp_m = 0;
                
                for (int k = 0; k < NoAgeG ; k++)
                {
                    if (get_ll::popsize_n[k] < 1)
                        pop_n = 0;
                    else
                        pop_n = 1.0/get_ll::popsize_n[k];
                    
                    if (get_ll::popsize_m[k] < 1)
                        pop_m = 0;
                    else
                        pop_m = 1.0/get_ll::popsize_m[k];
                    
                    
                    I_temp_m += (x[3+0*mg+ag*k]*alpha_i+x[4+0*mg+ag*k]+a1*(x[8+0*mg+ag*k]*alpha_i+x[9+0*mg+ag*k])+a2*(x[13+0*mg+ag*k]*alpha_i+x[14+0*mg+ag*k])+a3*(x[18+0*mg+ag*k]*alpha_i+x[19+0*mg+ag*k]))*(qp*(pwp_p[j][k]+qc*pwp_c[j][k]))*pop_m + (x[3+1*mg+ag*k]*alpha_i+x[4+1*mg+ag*k]+a1*(x[8+1*mg+ag*k]*alpha_i+x[9+1*mg+ag*k])+a2*(x[13+1*mg+ag*k]*alpha_i+x[14+1*mg+ag*k])+a3*(x[18+1*mg+ag*k]*alpha_i+x[19+1*mg+ag*k]))*(qp*(pwn_p[j][k]+qc*pwn_c[j][k]))*pop_n;
                    
                    I_temp_n += (x[3+0*mg+ag*k]*alpha_i+x[4+0*mg+ag*k]+a1*(x[8+0*mg+ag*k]*alpha_i+x[9+0*mg+ag*k])+a2*(x[13+0*mg+ag*k]*alpha_i+x[14+0*mg+ag*k])+a3*(x[18+0*mg+ag*k]*alpha_i+x[19+0*mg+ag*k]))*(qp*(nwp_p[j][k]+qc*nwp_c[j][k]))*pop_m + (x[3+1*mg+ag*k]*alpha_i+x[4+1*mg+ag*k]+a1*(x[8+1*mg+ag*k]*alpha_i+x[9+1*mg+ag*k])+a2*(x[13+1*mg+ag*k]*alpha_i+x[14+1*mg+ag*k])+a3*(x[18+1*mg+ag*k]*alpha_i+x[19+1*mg+ag*k]))*(qp*(nwn_p[j][k]+qc*nwn_c[j][k]))*pop_n;
                }
                
                for (int s = 0; s < 2; s++)
                {
                    if (s == 0)
                    {
                        prop_mg = p_mat[j];
                        I_temp = I_temp_m;
                        catsize = popsize_m[j];
                    }
                    else if (s == 1)
                    {
                        prop_mg = 1-p_mat[j];
                        I_temp = I_temp_n;
                        catsize = popsize_n[j];
                    }
                    
                    int pj = max(ag*(j-1),0);
                    //cout << pj << endl;
                    double muBv;
                    if (j == 0 && s == 1)
                        muBv = muB;
                    else
                        muBv = 0;
                    
                    dxdt[0+s*mg+j*ag] = (1.0-p_vul)*muBv - x[0+s*mg+j*ag]*xi                                     - x[0+s*mg+j*ag]*eta[j+1] + (x[0+pj]+x[0+mg+pj])*eta[j]*prop_mg;
                    
                    dxdt[1+s*mg+j*ag] =  p_vul*muBv + x[0+s*mg+j*ag]*xi   - x[1+s*mg+j*ag]*I_temp*beta           - x[1+s*mg+j*ag]*eta[j+1] + (x[1+pj]+x[1+mg+pj])*eta[j]*prop_mg;
                    
                    dxdt[2+s*mg+j*ag] = x[1+s*mg+j*ag]*I_temp*beta               - x[2+s*mg+j*ag]*si             - x[2+s*mg+j*ag]*eta[j+1] + (x[2+pj]+x[2+mg+pj])*eta[j]*prop_mg;
                    dxdt[3+s*mg+j*ag] = x[2+s*mg+j*ag]*si*pA[j]                - x[3+s*mg+j*ag]*ga0*rho          - x[3+s*mg+j*ag]*eta[j+1] + (x[3+pj]+x[3+mg+pj])*eta[j]*prop_mg;
                    dxdt[4+s*mg+j*ag] = x[2+s*mg+j*ag]*si*(1.0-pA[j])            - x[4+s*mg+j*ag]*ga0            - x[4+s*mg+j*ag]*eta[j+1] + (x[4+pj]+x[4+mg+pj])*eta[j]*prop_mg;
                    dxdt[5+s*mg+j*ag] = x[4+s*mg+j*ag]*ga0 + x[3+s*mg+j*ag]*ga0*rho        - x[5+s*mg+j*ag]*om   - x[5+s*mg+j*ag]*eta[j+1] + (x[5+pj]+x[5+mg+pj])*eta[j]*prop_mg;
                    
                    dxdt[6+s*mg+j*ag] = x[5+s*mg+j*ag]*om                      - d1*x[6+s*mg+j*ag]*I_temp*beta   - x[6+s*mg+j*ag]*eta[j+1] + (x[6+pj]+x[6+mg+pj])*eta[j]*prop_mg;
                    dxdt[7+s*mg+j*ag] = d1*x[6+s*mg+j*ag]*I_temp*beta            - x[7+s*mg+j*ag]*si             - x[7+s*mg+j*ag]*eta[j+1] + (x[7+pj]+x[7+mg+pj])*eta[j]*prop_mg;
                    dxdt[8+s*mg+j*ag] = x[7+s*mg+j*ag]*si*pA[j]                - x[8+s*mg+j*ag]*ga1*rho          - x[8+s*mg+j*ag]*eta[j+1] + (x[8+pj]+x[8+mg+pj])*eta[j]*prop_mg;
                    dxdt[9+s*mg+j*ag] = x[7+s*mg+j*ag]*si*(1.0-pA[j])            - x[9+s*mg+j*ag]*ga1            - x[9+s*mg+j*ag]*eta[j+1] + (x[9+pj]+x[9+mg+pj])*eta[j]*prop_mg;
                    dxdt[10+s*mg+j*ag] = x[9+s*mg+j*ag]*ga1 + x[8+s*mg+j*ag]*ga1*rho       - x[10+s*mg+j*ag]*om  - x[10+s*mg+j*ag]*eta[j+1] +(x[10+pj]+x[10+mg+pj])*eta[j]*prop_mg;
                    
                    dxdt[11+s*mg+j*ag] = x[10+s*mg+j*ag]*om                     - d2*x[11+s*mg+j*ag]*I_temp*beta - x[11+s*mg+j*ag]*eta[j+1] +(x[11+pj]+x[11+mg+pj])*eta[j]*prop_mg;
                    dxdt[12+s*mg+j*ag] = d2*x[11+s*mg+j*ag]*I_temp*beta           - x[12+s*mg+j*ag]*si           - x[12+s*mg+j*ag]*eta[j+1] +(x[12+pj]+x[12+mg+pj])*eta[j]*prop_mg;
                    dxdt[13+s*mg+j*ag] = x[12+s*mg+j*ag]*si*pA[j]              - x[13+s*mg+j*ag]*ga2*rho         - x[13+s*mg+j*ag]*eta[j+1] +(x[13+pj]+x[13+mg+pj])*eta[j]*prop_mg;
                    dxdt[14+s*mg+j*ag] = x[12+s*mg+j*ag]*si*(1.0-pA[j])           - x[14+s*mg+j*ag]*ga2          - x[14+s*mg+j*ag]*eta[j+1] +(x[14+pj]+x[14+mg+pj])*eta[j]*prop_mg;
                    dxdt[15+s*mg+j*ag] = x[14+s*mg+j*ag]*ga2 + x[13+s*mg+j*ag]*ga2*rho       - x[15+s*mg+j*ag]*om- x[15+s*mg+j*ag]*eta[j+1] +(x[15+pj]+x[15+mg+pj])*eta[j]*prop_mg;
                    
                    dxdt[16+s*mg+j*ag] = x[15+s*mg+j*ag]*om + x[20+s*mg+j*ag]*om - d3*x[16+s*mg+j*ag]*I_temp*beta- x[16+s*mg+j*ag]*eta[j+1] +(x[16+pj]+x[16+mg+pj])*eta[j]*prop_mg;
                    dxdt[17+s*mg+j*ag] = d3*x[16+s*mg+j*ag]*I_temp*beta           - x[17+s*mg+j*ag]*si           - x[17+s*mg+j*ag]*eta[j+1] +(x[17+pj]+x[17+mg+pj])*eta[j]*prop_mg;
                    dxdt[18+s*mg+j*ag] = x[17+s*mg+j*ag]*si*pA[j]                 - x[18+s*mg+j*ag]*ga3*rho      - x[18+s*mg+j*ag]*eta[j+1] +(x[18+pj]+x[18+mg+pj])*eta[j]*prop_mg;
                    dxdt[19+s*mg+j*ag] = x[17+s*mg+j*ag]*si*(1.0-pA[j])           - x[19+s*mg+j*ag]*ga3          - x[19+s*mg+j*ag]*eta[j+1] +(x[19+pj]+x[19+mg+pj])*eta[j]*prop_mg;
                    dxdt[20+s*mg+j*ag] = x[19+s*mg+j*ag]*ga3 + x[18+s*mg+j*ag]*ga3*rho - x[20+s*mg+j*ag]*om      - x[20+s*mg+j*ag]*eta[j+1] +(x[20+pj]+x[20+mg+pj])*eta[j]*prop_mg;
                    
                    dxdt[21+s*mg+j*ag] = x[21+s*mg+j*ag]*eta[j+1] - x[21+s*mg+j*ag]*eta[j+1];
                    
                    dxdt[22+s*mg+j*ag] = si*(x[2+s*mg+j*ag] + x[7+s*mg+j*ag] + x[12+s*mg+j*ag] + x[17+s*mg+j*ag]);
                }
            }
        }
    };
    
    
    vector2D get_foi(param::param_state_t& pars, num_vec x, double t)
    {
        double ag = 46; double mg = 23;
        double alpha_i = pars.mapofParam.find("alpha_i")->second;
        double qp = (pars.mapofParam.find("qp")->second);
        double qc = (pars.mapofParam.find("qc")->second);
        double a1 = 1; double a2 = 1; double a3 = 1;
        double I_temp_nwn, I_temp_nwp, I_temp_pwn, I_temp_pwp ;
        double pop_n, pop_m;
        num_vec I_pwp; num_vec I_nwp; num_vec I_pwn; num_vec I_nwn;
        vector2D I;
        
        for(int a = 0; a < NoAgeG; a++)
        {
            I_temp_nwn = I_temp_nwp = I_temp_pwn = I_temp_pwp = 0;
            for (int k = 0; k < NoAgeG ; k++)
            {
                if (get_ll::popsize_n[k] < 1)
                    pop_n = 0;
                else
                    pop_n = 1.0/get_ll::popsize_n[k];
                
                if (get_ll::popsize_m[k] < 1)
                    pop_m = 0;
                else
                    pop_m = 1.0/get_ll::popsize_m[k];
                
                I_temp_nwp = (x[3+0*mg+ag*k]*alpha_i+x[4+0*mg+ag*k]+a1*(x[8+0*mg+ag*k]*alpha_i+x[9+0*mg+ag*k])+a2*(x[13+0*mg+ag*k]*alpha_i+x[14+0*mg+ag*k])+a3*(x[18+0*mg+ag*k]*alpha_i+x[19+0*mg+ag*k]))*(qp*(nwp_p[a][k]+qc*nwp_c[a][k]))*pop_m;
                
                I_temp_nwn = (x[3+1*mg+ag*k]*alpha_i+x[4+1*mg+ag*k]+a1*(x[8+1*mg+ag*k]*alpha_i+x[9+1*mg+ag*k])+a2*(x[13+1*mg+ag*k]*alpha_i+x[14+1*mg+ag*k])+a3*(x[18+1*mg+ag*k]*alpha_i+x[19+1*mg+ag*k]))*(qp*(nwn_p[a][k]+qc*nwn_c[a][k]))*pop_n;
                
                I_temp_pwp = (x[3+0*mg+ag*k]*alpha_i+x[4+0*mg+ag*k]+a1*(x[8+0*mg+ag*k]*alpha_i+x[9+0*mg+ag*k])+a2*(x[13+0*mg+ag*k]*alpha_i+x[14+0*mg+ag*k])+a3*(x[18+0*mg+ag*k]*alpha_i+x[19+0*mg+ag*k]))*(qp*(pwp_p[a][k]+qc*pwp_c[a][k]))*pop_m;
                
                I_temp_pwn = (x[3+1*mg+ag*k]*alpha_i+x[4+1*mg+ag*k]+a1*(x[8+1*mg+ag*k]*alpha_i+x[9+1*mg+ag*k])+a2*(x[13+1*mg+ag*k]*alpha_i+x[14+1*mg+ag*k])+a3*(x[18+1*mg+ag*k]*alpha_i+x[19+1*mg+ag*k]))*(qp*(pwn_p[a][k]+qc*pwn_c[a][k]))*pop_n;
                
                I_nwp.push_back(I_temp_nwp);
                I_nwn.push_back(I_temp_nwn);
                I_pwp.push_back(I_temp_pwp);
                I_pwn.push_back(I_temp_pwn);
            }
        }
        I.push_back(I_nwn);I.push_back(I_nwp);
        I.push_back(I_pwn);I.push_back(I_pwp);
        return I;
    }
    
    double get_pR(param::param_state_t& pars, num_vec x, double t)
    {
        double sum_wcb; num_vec num_vec_wcba;
        double ag = 46; double mg = 23;
        for(int j = 18; j < 21 ; j++)
        {
            double CB2_temp = (x[1+mg+j*ag] + x[6+mg+j*ag] + x[11+mg+j*ag] + x[16+mg+j*ag] + x[2+mg+j*ag] + x[7+mg+j*ag] + x[12+mg+j*ag] + x[17+mg+j*ag])/((double)get_ll::popsize_n[j]);
            num_vec_wcba.push_back(CB2_temp/3.0);
        }
        for(int i = 0; i < boost::size(num_vec_wcba) ;i++){sum_wcb = sum_wcb + num_vec_wcba[i];}
        return sum_wcb;
    }
};
using namespace get_ll;


#endif /* model_h */
