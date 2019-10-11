//
//  cea.h
//  calibration
//
//  Created by David Hodgson on 08/08/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

#ifndef cea_h
#define cea_h

struct cea_state_t
{
    // Empty vectors
    num_vec S_tot; // total number of symptomatic infections
    num_vec H_tot; // total number of hospitalised cases
    num_vec D_tot; // total number of deaths
    num_vec GP_tot;// total number of GP consultations
    num_vec BD_tot;// total number of hospital bed days

    num_vec doses_pal; // total number of doses of palivizumab
    num_vec doses_pro; // total number of doses of prophylactics
    
    num_vec inci; // incidence in each age group, maternal group and risk group at equilibirum year
    num_vec inci_pri; // incidence of primary incidence in each year group

    double Q, CP, CT; // Total QALY loss, cost of purchasing, cost of administration
    
    int time_hor;   // time horizon
    double disc;    // discounting
    
    // initial values in constructor
    cea_state_t()
    {
        time_hor = 10;
        disc = 0.035;

        for (int i = 0; i < NoAgeG; i++)
        {
            S_tot.push_back(0); H_tot.push_back(0); D_tot.push_back(0);
            GP_tot.push_back(0); BD_tot.push_back(0);
            inci_pri.push_back(0);
        }
        for (int i = 0; i < time_hor*52; i++)
        {
            doses_pal.push_back(0); doses_pro.push_back(0);
        }
        for (int i = 0; i < NoAgeG*9; i++)
        {
            inci.push_back(0);
        }
        
        Q = CP = CT = 0;
    }
    
    // clear all the values
    void cea_state_clear(cea_state_t& cea_state)
    {
        for (int i = 0; i < NoAgeG; i++)
        {
            cea_state.S_tot[i] = 0; cea_state.H_tot[i] = 0; cea_state.D_tot[i] = 0;
            cea_state.GP_tot[i] = 0; cea_state.BD_tot[i] = 0; cea_state.inci_pri[i] = 0;
        }
        for (int i = 0; i < time_hor*52; i++)
        {
            cea_state.doses_pal[i] = 0; cea_state.doses_pro[i] = 0;
        }
        for (int i = 0; i < time_hor*52; i++)
        {
            cea_state.inci[i] = 0;
        }
        
        cea_state.Q = cea_state.CP = cea_state.CT = 0;
    }
};


// get_QALY -> convert the cases averted to QALY loss
double get_QALY(num_vec inci, amh::amh_state_t& mcmc_state, num_vec& S_tot, num_vec& H_tot, num_vec& D_tot, int s)
{
    num_vec S_a; num_vec H_a; num_vec D_a;

    PRNG_s rng(s);
    
    S_a = get_S(inci, mcmc_state, s);
    H_a = get_H(inci, s);
    D_a = get_D(inci, s);
    
    // Keep track of totals
    for (int a = 0; a < NoAgeG; a++)
    {
        S_tot[a] += S_a[a]; H_tot[a] += H_a[a]; D_tot[a] += D_a[a];
    }
    
    double tot_Q = 0;
    double val = 0;
    
    for (int a = 0; a < 16; a++){
        val = Q_S1(rng);
        tot_Q += S_a[a]*val;
    }
    
    for (int a = 16; a < NoAgeG; a++){
        val = Q_S2(rng);
        tot_Q += S_a[a]*val;
    }
    
    for (int a = 0; a < 16; a++){
        val = Q_H1(rng);
        tot_Q += H_a[a]*val;
    }
    
    for (int a = 16; a < NoAgeG; a++){
        val = Q_H2(rng);
        tot_Q += H_a[a]*val;
    }
    
    for (int a = 0; a < NoAgeG; a++){
        val = Q_D[a](rng);
        tot_Q += D_a[a]*val;
    }
    return tot_Q;
}


// get_QALY -> convert the cases averted to QALY loss
double get_QALY_2(num_vec inci, amh::amh_state_t& mcmc_state, num_vec& S_tot, num_vec& GP_tot, num_vec& H_tot, num_vec& D_tot, int s)
{
    num_vec S_a; num_vec GP_a; num_vec H_a; num_vec D_a;
    
    PRNG_s rng(s);
    
    S_a = get_S(inci, mcmc_state, s);
    GP_a = get_GP(inci, s);
    H_a = get_H(inci, s);
    D_a = get_D(inci, s);
    
    // Keep track of totals
    for (int a = 0; a < NoAgeG; a++)
    {
        S_tot[a] += S_a[a]; GP_tot[a] += GP_a[a]; H_tot[a] += H_a[a]; D_tot[a] += D_a[a];
    }
    
    double tot_Q = 0;
    double val = 0;
    
    for (int a = 0; a < 16; a++){
        val = Q_NHS1(rng);
        tot_Q += (S_a[a] - GP_a[a] - H_a[a])*val;
    }
    
    for (int a = 16; a < NoAgeG; a++){
        val = Q_NHS2(rng);
        tot_Q += (S_a[a] - GP_a[a] - H_a[a])*val;
    }
    
    for (int a = 0; a < 16; a++){
        val = Q_HS1(rng);
        tot_Q += (GP_a[a] + H_a[a])*val;
    }
    
    for (int a = 16; a < NoAgeG; a++){
        val = Q_HS2(rng);
        tot_Q += (GP_a[a] + H_a[a])*val;
    }
    
    for (int a = 0; a < NoAgeG; a++){
        val = Q_D[a](rng);
        tot_Q += D_a[a]*val;
    }
    return tot_Q;
}

// get_CostT -> convert the cases averted to cost
double get_CostT(num_vec inci, num_vec& GP_tot, num_vec& BD_tot, int s)
{
    num_vec GP_a; num_vec BD_a;

    PRNG_s rng(s);
    GP_a = get_GP(inci, s);
    BD_a = get_BD(inci, s);
    
    // Keep track of totals
    for (int a = 0; a < NoAgeG; a++)
    {
        GP_tot[a] += GP_a[a]; BD_tot[a] += BD_a[a];
    }
    double tot_CT = 0;
    double GP_c = 37.4;
    
    for (int a = 0; a < NoAgeG; a++)
        tot_CT += GP_a[a]*GP_c;
    
    for (int a = 0; a < 16; a++)
        tot_CT += BD_a[a]*C_H1(rng);
    
    for (int a = 16; a < NoAgeG; a++)
        tot_CT += BD_a[a]*C_H1(rng);
    
    return tot_CT;
}

// get_CostT -> determine the cost of administratoin and purchasing the prophylactic
double get_CostP(double no_pal, double no_vac, double p_ad, cea_state_t& cea_state, int t_w, int s)
{
    cea_state.doses_pal[t_w] = no_pal;
    cea_state.doses_pro[t_w] = no_vac;
    return no_pal*57.5 + no_vac*p_ad;
}


namespace CEA{
     //Import the relevant files (PXCP, PXCT, PXQ, PXDose)
    class CSVReader
    {
        std::string fileName;
        std::string delimeter;
    
    public:
        CSVReader(std::string filename, string delm = ",") : fileName(filename), delimeter(delm)
        {}
    vector2D getData;
    
    };
};
#endif /* cea_h */
