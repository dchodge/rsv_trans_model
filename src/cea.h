//
//  cea.h
//  calibration
//
//  Created by David Hodgson on 08/08/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

#ifndef cea_h
#define cea_h


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
        S_tot[a] += S_a[a];H_tot[a] += H_a[a];D_tot[a] += D_a[a];
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

double get_CostP(double no_pal, double no_vac, double p_ad, vector2D& doses, int t_w, int s)
{
    doses[0][t_w] = no_pal;
    doses[1][t_w] = no_vac;
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
    
   /* vector2D CSVReader::getData()
    {
        std::ifstream file(FileName);
        vector2D dataList;
        std::string line = "";
        
        while (getline(file, line))
        {
            num_vec vec;
            boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
            dataList.push_back(vec);
        }
        file.close();
        
        return dataList;
    }
    string PBQ_f = get_ll::dout + get_ll::calout_name + "inter/P_/qaly.csv"; //The source of the contact data from AJ STudy (<1yrs)
    CSVReader reader(PBQ_f);
    vector2D dataList = reader.getData();
    for (num_vec vec : dataList)
    {
        for(double data : vec)
        {
            cout << data << ",";
        }
        cout << endl;
    }*/
    
    
   // vector2D const PBQ = get_2d(PBQ_f, 2, 52*10);

    
    
   /* double get_incQALY(double )
    {
        
    }
    
    double get_ppd_pal(double )
    {
        
    }

    double get_ppd_int(double )
    {
        
    }*/

};
#endif /* cea_h */
