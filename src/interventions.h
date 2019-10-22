//
//  interventions.h
//  calibration
//
//  Created by David Hodgson on 12/07/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

#ifndef interventions_h
#define interventions_h


namespace manip_contact {
    
    /**********************************************************************/
    /** CwX Contacts for persons who are cocooned ****/
    /**********************************************************************/
    
    vector2D get_cwn(double prop_c, char s)
    {
        vector2D cwn_e(NoAgeG,vector<double>(NoAgeG,0));
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++)        //infants with
                for (int j = 0; j < NoAgeG; j++)        //infants
                    cwn_e[i][j] = (nwn_p[i][j]);
            
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwn_e[i][j] = (nwn_p[i][j])*(1-prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++)        //infants with
                for (int j = 0; j < NoAgeG; j++)        //infants
                    cwn_e[i][j] = (pwn_p[i][j]);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwn_e[i][j] = (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(1-prop_c);

        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++)        //infants with
                for (int j = 0; j < NoAgeG; j++)        //infants
                    cwn_e[i][j] = (nwn_c[i][j]);
            
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwn_e[i][j] = (nwn_c[i][j])*(1-prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++)        //infants with
                for (int j = 0; j < NoAgeG; j++)        //infants
                    cwn_e[i][j] = (pwn_c[i][j]);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwn_e[i][j] = (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return cwn_e;
    }
    
    vector2D get_cwp(double prop_c, char s)
    {
        vector2D cwp_e(NoAgeG,vector<double>(NoAgeG,0));
        
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwp_e[i][j] = (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(p_mat[j])*(1-prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwp_e[i][j] = (pwp_p[i][j])*(1-prop_c); // parental w. infant plus non-parental with infant who are cocooned
    
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwp_e[i][j] = (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(p_mat[j])*(1-prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwp_e[i][j] = (pwp_c[i][j])*(1-prop_c); // parental w. infant plus non-parental with infant who are cocooned
        }
        else
            cout << "Error cont" << endl;
        
        return cwp_e;
    }
    
    vector2D get_cwc(double prop_c, char s)
    {
        vector2D cwc_e(NoAgeG,vector<double>(NoAgeG,0));

        if (s == 'p')
        {
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwc_e[i][j] = (nwn_p[i][j])*(prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwc_e[i][j] = cnt_matrix_p_h[i][j] + (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(prop_c);
            
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwc_e[i][j] = cnt_matrix_p_h[i][j]/2.0 + (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(p_mat[j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwc_e[i][j] = (pwp_p[i][j])*(prop_c); // parental w. infant plus non-parental with infant who are cocooned
    
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwc_e[i][j] = (nwn_c[i][j])*(prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwc_e[i][j] = cnt_matrix_c_h[i][j] + (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(prop_c);
            
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwc_e[i][j] = cnt_matrix_c_h[i][j]/2.0 + (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(p_mat[j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwc_e[i][j] = (pwp_c[i][j])*(prop_c); // parental w. infant plus non-parental with infant who are cocooned

        }
        else
            cout << "Error cont" << endl;
        return cwc_e;
    }
    
    /**********************************************************************/
    /** PwX Contacts for persons who are mothers but not cocooned ****/
    /**********************************************************************/
    
    
    vector2D get_pwn(double prop_c, char s)
    {
        vector2D pwn_e(NoAgeG,vector<double>(NoAgeG,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < NoAgeG; j++) // all
                    pwn_e[i][j] = pwn_p[i][j]; // all contacts
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwn_e[i][j] = cnt_matrix_p_h[i][j] + (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(1-prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < NoAgeG; j++) // all
                    pwn_e[i][j] = pwn_c[i][j]; // all contacts
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwn_e[i][j] = cnt_matrix_c_h[i][j] + (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return pwn_e;
    }
    
    vector2D get_pwp(double prop_c, char s)
    {
        vector2D pwp_e(NoAgeG,vector<double>(NoAgeG,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // all
                    pwp_e[i][j] = pwp_p[i][j]*(1-prop_c); // all contacts
        }
        else if (s == 'c')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // all
                    pwp_e[i][j] = pwp_c[i][j]*(1-prop_c); // all contacts
        }
        else
            cout << "Error cont" << endl;
        
        return pwp_e;
    }
    
    vector2D get_pwc(double prop_c, char s)
    {
        // Non cocooned parent -> cocooned infant (family and outside family)
        vector2D pwc_e(NoAgeG,vector<double>(NoAgeG,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwc_e[i][j] = (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    pwc_e[i][j] = pwn_p[i][j]*prop_c; //
        }
        else if (s == 'c')
        {
            // Parent -> parents (all)
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwc_e[i][j] = (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    pwc_e[i][j] = pwn_c[i][j]*prop_c; //
        }
        else
            cout << "Error cont" << endl;
        
        return pwc_e;
    }
    
    
 /**********************************************************************/
/** NwX Contacts for persons who are neither cocooned nor mothers ****/
 /**********************************************************************/
    
    vector2D get_nwn(double prop_c, char s)
    {
        // cnt_matrix_p cnt_matrix_p_h
        // (cHij[[i]][[j]]/2.0 ) + (cij[[i]][[j]] - cHij[[i]][[j]])*(1 - proppar[[j]])
        vector2D nwn_e(NoAgeG,vector<double>(NoAgeG,0));
        if (s == 'p')
        {
            for (int i = 0; i < NoAgeG; i++) // population with.
                for (int j = 0; j < NoAgeG; j++) //     infants
                    nwn_e[i][j] = nwn_p[i][j]; //number of contacts made by population with non-cocooned prop_ci
            
            for (int i = 0; i < NoAgeG; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwn_e[i][j] = nwn_p[i][j]*(1-prop_c); //number
        }
        else if (s == 'c')
        {
            for (int i = 0; i < NoAgeG; i++) // population with.
                for (int j = 0; j < NoAgeG; j++) //     infants
                    nwn_e[i][j] = nwn_c[i][j]; //number of contacts made by population with non-cocooned prop_ci
            
            for (int i = 0; i < NoAgeG; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwn_e[i][j] = nwn_c[i][j]*(1-prop_c); //numbe

        }
        else
            cout << "Error cont" << endl;

        return nwn_e;
    }
    
    vector2D get_nwp(double prop_c, char s)
    {
        // Population with non-cocooned parents)
        vector2D nwp_e(NoAgeG,vector<double>(NoAgeG,0));
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = (cnt_matrix_p_h[i][j]/2.0) + (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*p_mat[j]*(1-prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = cnt_matrix_p[i][j]*p_mat[j]*(1-prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = (cnt_matrix_c_h[i][j]/2.0) + (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*p_mat[j]*(1-prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = cnt_matrix_c[i][j]*p_mat[j]*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return nwp_e;
    }
    
    vector2D get_nwc(double prop_c, char s)
    {
        
        vector2D nwc_e(NoAgeG,vector<double>(NoAgeG,0));
        if (s == 'p')
        {
            for (int i = 0; i < NoAgeG; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwc_e[i][j] = nwn_p[i][j]*(prop_c); //numbe
            
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = (cnt_matrix_p[i][j] - cnt_matrix_p_h[i][j])*p_mat[j]*(prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = cnt_matrix_p[i][j]*p_mat[j]*(prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 0; i < NoAgeG; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwc_e[i][j] = nwn_c[i][j]*(prop_c); //numbe
            
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = (cnt_matrix_c[i][j] - cnt_matrix_c_h[i][j])*p_mat[j]*(prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = cnt_matrix_c[i][j]*p_mat[j]*(prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return nwc_e;
    }
}






namespace cal
{
    // Indicator function for target groups with age (25 age groups)
    num_vec G_base = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    num_vec G_0mo = {1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    num_vec G_plus = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    num_vec G_2mo = {0,0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    num_vec G_pal = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    num_vec G_2_4 =  {0,0,0,0,0,0,0,0,0,0,0,0,0, 1.0, 1.0, 1.0,           0,0,0,0,0,0,0,0,0};
    num_vec G_5_10 = {0,0,0,0,0,0,0,0,0,0,0,0,0, 0.0, 0.0, 0.0, 1.0,        0,0,0,0,0,0,0,0};
    num_vec G_5_14 = {0,0,0,0,0,0,0,0,0,0,0,0,0, 0.0, 0.0, 0.0, 1.0, 1.0, 0,  0,  0,  0,0,0,0};
    num_vec G_par =  {0,0,0,0,0,0,0,0,0,0,0,0,0, 0,   0,   0,   0,   0,   1.0,1.0,1.0,0,0,0,0};
    num_vec G_65_ = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 1.0, 1.0};
    num_vec G_75_ = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 1.0};
    
    // Proportion of total coverage achieved by month i for
    num_vec up_week_o65 = {0., 0.046, 0.177, 0.311, 0.444, 0.64, 0.737, 0.802, 0.856, 0.894, 0.920, 0.941, 0.953, 0.967, 0.978, 0.980, 0.986, 0.988, 0.988, 0.996, 1.}; //over 65
    num_vec up_week_u65 = {0., 0.031, 0.075, 0.16, 0.274, 0.387, 0.491, 0.579, 0.664, 0.726, 0.783, 0.830, 0.868, 0.906, 0.931, 0.953, 0.953, 0.965, 0.981, 0.987, 1.}; //under 65
    num_vec up_week_2t3 = {0., 0.000, 0.000, 0.00, 0.040, 0.120, 0.227, 0.369, 0.511, 0.629, 0.719, 0.809, 0.868, 0.910, 0.941, 0.965, 0.969, 0.976, 0.990, 0.993, 1.}; //2-3 year olds
    num_vec up_week_preg ={0., 0.041, 0.095, 0.192, 0.32, 0.451, 0.568, 0.664, 0.734, 0.794, 0.843, 0.881, 0.917, 0.934, 0.944, 0.967, 0.974, 0.973, 0.983, 0.990, 1.}; //pregnant women
    
    struct inter_data_t
    {
        // Each vector has 14 entries corresponding to the 14 intervention programme of the order:
        // None, PAL-VHR-S, MAB-VHR-S, MAB-HR-S, MAB-HR-A, MAB-ALL-S, MAB-ALL-A, MAT-S, MAT-A, VAC-INF-S, VAC-INF-A, VAC-2-4-S,VAC-5-9-S, VAC-5-14-S, VAC-75-S, VAC-75-S
        str_vec prog_no; // Folder names for output-files
        
        str_vec prog_name;       // Name of programme (Redundant?)
        str_vec cal_type;        // Type of calendar used
        vector2D t_group;        // Target age group
        num_vec c_ad;            // Cost of prophylatics (palivizumab cost = 57.5)
        
        vector2D uprate;         //  Uptake
        num_vec cov;             // Coverage
        num_vec cov_c;           // Coverage of maternal group (proportion in group c)
        num_vec start_w;         // Start week (0=Jul 1st)
        num_vec end_w;           // End week (0=Jul 1st)
        num_vec Pal_ind;         // Indicator for Palivizumab
        
        // efficacy holders
        VectorXd eff_pal;     // Palivizumab
        VectorXd eff_mab;     // long-acting monoclonal antibodies
        VectorXd eff_vac;     // vaccines (elderly + infants)
        VectorXd eff_mat;     // maternal
        VectorXd vac_delay;   // delay between vacciantion and immunity (due to build up of antibodies)

        
        num_vec rate;
        double om_mab; double xi_b;
        
        inter_data_t()
        {
            prog_no =   {"NONE",   "PAL_VHR_S", "MAB_VHR_S", "MAB_HR_S", "MAB_HR_S+",  "MAB_ALL_S",  "MAB_ALL_S+",  "MAT_S",  "MAT_A",  "VAC_INF_S",  "VAC_INF_A",  "VAC_2_4_S","VAC_5_9_S", "VAC_5_14_S", "VAC_75_S", "VAC_65_S"};
            prog_name = {"Base_","Pal_","mABVHR_","mABHR_","mABHR_","mAB_","mAB_","matS_","matA_","infS_"  ,"infA_"  ,"Pre"    ,"Sch1_" ,"Sch2_"  ,"Eld1_"   ,"Eld2_"};
            cal_type = {"None", "None", "Mvhr", "Mhr", "Mhr_p", "Mlr", "Mlr_p", "mat", "mat", "LAV_inf", "LAV_inf", "LAV_ald", "LAV_ald", "LAV_ald", "LAV_ald", "LAV_ald"};
            t_group = {cal::G_base, cal::G_base, cal::G_base, cal::G_0mo, cal::G_0mo, cal::G_0mo, cal::G_0mo, cal::G_par, cal::G_par, cal::G_2mo, cal::G_2mo, cal::G_2_4, cal::G_5_10, cal::G_5_14, cal::G_75_, cal::G_65_};
            c_ad =   {0, 0, 11.5,11.5, 11.5, 11.5, 11.5,       9, 9, 9, 9, 9, 9, 9, 9, 9};
            uprate = {{0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, cal::up_week_2t3, cal::up_week_u65, cal::up_week_u65, cal::up_week_o65, cal::up_week_o65};
            cov = {0.0, 0.9, 0.9,0.9, 0.9, 0.9, 0.9, 0.6, 0.6, 0.9, 0.9, 0.45, 0.6, 0.6, 0.7, 0.7};
            cov_c = {0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            start_w = {0.0,  16, 16, 12,  12,    12,   12,       0,    0,  4,    0,    16,    16,    16,    16,     16};
            end_w   = {0.0,  37, 37, 12+21, 12+21, 12+21, 12+21, 0+21, 52, 4+21, 52,   16+21, 8+21, 8+21, 16+21, 16+21};
            Pal_ind = {false, true, false, false, false, false, false, true, true, true, true, true, true, true, true, true, true};
            
            om_mab = 1.0/250.0;
            xi_b = 1.0;
        }
        
        void get_eff(inter_data_t& inter_data, int ls)
        {
            // Distributions assocaited with efficacy
            boost::random::weibull_distribution<> Pal_E(12.4311, 0.772251);
            boost::random::weibull_distribution<> Mab_E(11.8975, 0.7317);
            boost::random::weibull_distribution<> Vac_E(31.4637, 0.844676);
            boost::random::weibull_distribution<> Mat_E(3.32676, 0.461343);
            
            inter_data.eff_pal = VectorXd::Zero(ls); inter_data.eff_mab = VectorXd::Zero(ls);
            inter_data.eff_vac = VectorXd::Zero(ls); inter_data.eff_mat = VectorXd::Zero(ls);

            // Get samples of efficacies from distributions
            for (int s = 0; s < ls; s++)
            {
                inter_data.eff_pal(s) = Pal_E(rng);
                inter_data.eff_mab(s) = Mab_E(rng);
                inter_data.eff_vac(s) = Vac_E(rng);
                inter_data.eff_mat(s) = Mat_E(rng);
            }
            
        }
    };
    
    // Turns weekly uptake ratein daily uptake rates
    num_vec gen_daily(num_vec up_week_raw, int week_s)
    {
        num_vec up_week(52,0); num_vec up_day(365,0);
        
        int pos_week, pos_ist;
        for (int k = week_s; k < week_s + 21; k++)
        {
            pos_week = (k)%52;
            pos_ist = (k-week_s);
            up_week[pos_week] = up_week_raw[pos_ist];
        }
        
        int pos_day;
        
        for (int i = 0 ; i < 365; i++)
        {
            pos_day = (i+week_s*7)%365;
            up_day[pos_day] = (up_week[(i/7+1+week_s)%52] - up_week[(i/7+week_s)%52])/7.0;
            if (up_day[pos_day] < 0)
                up_day[pos_day] = 0;
        }
        
        return up_day;
    }
    
    //Dating is 0 = 1st Jul, 1 = 1st Aug, 2 = 1st Sep, 3 = 1st Oct, 4 = 1st Nov, 5 = 1st Dec, 6 = 1st Jan, 7 = 1st Feb, 8 = 1st Mar, 9 = 1st Apr, 10 = 1st May, 11 = 1st Jun.
    struct Calendar_full
    {
        bool ind_OOB;
        double cov;
        int t_start, t_end, dt;
        num_vec t_group, rate;
        bool static_ind;
        char pal_ind;
        string proph_ind;
        MatrixXd init = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_pal = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_mAB_VHR = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_mAB_HR = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_mAB_LR = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_LAV_HR = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_LAV_LR = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_mat_LR = MatrixXd::Zero(365, NoAgeG);
        
        MatrixXd cal_pal_dose = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_mAB_VHR_dose = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_mAB_HR_dose = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_mAB_LR_dose = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_LAV_HR_dose = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_LAV_LR_dose = MatrixXd::Zero(365, NoAgeG);
        MatrixXd cal_mat_LR_dose = MatrixXd::Zero(365, NoAgeG);
        int s;

        Calendar_full(inter_data_t& inter_data, int s, int iN)
        {
            // Calendar type
            string proph_ind = inter_data.cal_type[iN];
            // Palivizumab used?
            bool pal_ind = inter_data.Pal_ind[iN];

            // Properties of calendar
            double t_start = inter_data.start_w[iN];
            double t_end = inter_data.end_w[iN];
            double cov = inter_data.cov[iN];
            num_vec t_group = inter_data.t_group[iN];
            num_vec rate = inter_data.rate;
            
            // Efficacy
            double eff_pal = inter_data.eff_pal[s];
            double eff_mab = inter_data.eff_mab[s];
            double eff_vac = inter_data.eff_vac[s];
            double eff_mat = inter_data.eff_mat[s];

            ind_OOB = false;
            if (pal_ind == true)
            {
                cal_pal = PAL_calendar(init, 'p', eff_pal);
                cal_pal_dose = PAL_calendar(init, 'd', eff_pal);
            }
            
            if (proph_ind == "Mvhr")
            {
                cal_mAB_VHR = mAB_vhr_calendar(init, t_start, t_end, 'p', eff_mab);
                cal_mAB_VHR_dose = mAB_vhr_calendar(init, t_start, t_end, 'd', eff_mab);
            }
            else if (proph_ind == "Mhr")
            {
                cal_mAB_VHR = mAB_vhr_calendar(init, t_start, t_end, 'p', eff_mab);
                cal_mAB_HR = mAB_calendar(init, t_group, cov, t_start, t_end, 'p', eff_mab);
                cal_mAB_VHR_dose = mAB_vhr_calendar(init, t_start, t_end, 'd', eff_mab);
                cal_mAB_HR_dose = mAB_calendar(init, t_group, cov, t_start, t_end, 'd', eff_mab);
            }
            else if (proph_ind == "Mhr_p")
            {
                cal_mAB_VHR = mAB_vhr_calendar(init, t_start, t_end, 'p', eff_mab);
                cal_mAB_HR = mAB_calendar_plus(init, t_group, cov, t_start, t_end, 'p', eff_mab);
                cal_mAB_VHR_dose = mAB_vhr_calendar(init, t_start, t_end, 'd', eff_mab);
                cal_mAB_HR_dose = mAB_calendar_plus(init, t_group, cov, t_start, t_end, 'd', eff_mab);
            }
            else if (proph_ind == "Mlr")
            {
                cal_mAB_VHR = mAB_vhr_calendar(init, t_start, t_end, 'p', eff_mab);
                cal_mAB_HR = mAB_calendar(init, t_group, cov, t_start, t_end, 'p', eff_mab);
                cal_mAB_LR = mAB_calendar(init, t_group, cov, t_start, t_end, 'p', eff_mab);
                cal_mAB_VHR_dose = mAB_vhr_calendar(init, t_start, t_end, 'd', eff_mab);
                cal_mAB_HR_dose = mAB_calendar(init, t_group, cov, t_start, t_end, 'd', eff_mab);
                cal_mAB_LR_dose = mAB_calendar(init, t_group, cov, t_start, t_end, 'd', eff_mab);
            }
            else if (proph_ind == "Mlr_p")
            {
                cal_mAB_VHR = mAB_vhr_calendar(init, t_start, t_end, 'p', eff_mab);
                cal_mAB_HR = mAB_calendar_plus(init, t_group, cov, t_start, t_end, 'p', eff_mab);
                cal_mAB_LR = mAB_calendar_plus(init, t_group, cov, t_start, t_end, 'p', eff_mab);
                cal_mAB_VHR_dose = mAB_vhr_calendar(init, t_start, t_end, 'd', eff_mab);
                cal_mAB_HR_dose = mAB_calendar_plus(init, t_group, cov, t_start, t_end, 'd', eff_mab);
                cal_mAB_LR_dose = mAB_calendar_plus(init, t_group, cov, t_start, t_end, 'd', eff_mab);
            }
            else if (proph_ind == "LAV_inf")
            {
                cal_LAV_HR = LAV_calendar(init, t_group, cov, rate, t_start, t_end, 2, 'p', eff_vac);
                cal_LAV_LR = LAV_calendar(init, t_group, cov, rate, t_start, t_end, 2, 'p', eff_vac);
                cal_LAV_HR_dose = LAV_calendar(init, t_group, cov, rate, t_start, t_end, 2, 'd', eff_vac);
                cal_LAV_LR_dose = LAV_calendar(init, t_group, cov, rate, t_start, t_end, 2, 'd', eff_vac);
            }
            else if (proph_ind == "LAV_ald")
            {
                cal_LAV_LR = LAV_calendar(init, t_group, cov, rate, t_start, t_end, 1, 'p', eff_vac);
                cal_LAV_LR_dose = LAV_calendar(init, t_group, cov, rate, t_start, t_end, 1, 'd', eff_vac);
            }
            else if (proph_ind == "mat")
            {
                cal_mat_LR = mat_calendar(init, t_group, cov, t_start, t_end, 'p', eff_mat, eff_vac);
                cal_mat_LR_dose = mat_calendar(init, t_group, cov, t_start, t_end, 'd', eff_mat, eff_vac);
            }
            else
                cout << "No cal" << endl;
        }
        
        MatrixXd PAL_calendar(MatrixXd calendar, char t, double eff_pal)
        {
            PRNG_s rng(0);
            MatrixXd sero = MatrixXd::Zero(365, NoAgeG);
            MatrixXd outmat = MatrixXd::Zero(365, NoAgeG);

            for (int i = 15*7 ; i < 36*7+1; i++)
                for (int j = 0 ; j < NoAgeG; j++)
                    calendar(i%365,j) = cal::G_pal[j]*0.9/30.0;
            
            for (int i = 15*7 ; i < 36*7+1; i++)
                for (int j = 0 ; j < NoAgeG; j++)
                    sero(i%365,j) = calendar(i%365,j)*eff_pal;
            
            if (t == 'p')
                outmat = sero;
            else
                outmat = calendar;
            
            return outmat;
        }
        
        MatrixXd mAB_vhr_calendar(MatrixXd calendar, int t_start, int t_end, char t, double eff_mab)
        {
            PRNG_s rng(0);
            MatrixXd sero = MatrixXd::Zero(365, NoAgeG);
            MatrixXd outmat = MatrixXd::Zero(365, NoAgeG);

            for (int i = 15*7 ; i < 36*7+1; i++)
                for (int j = 0 ; j < NoAgeG; j++)
                    calendar(i%365,j) = cal::G_pal[j]*0.9/30.0;
            
            for (int i = 15*7 ; i < 36*7+1; i++)
                for (int j = 0 ; j < NoAgeG; j++)
                    sero(i%365,j) = calendar(i%365,j)*eff_mab;
            
            if (t == 'p')
                outmat = sero;
            else
                outmat = calendar;
            
            return outmat;
        }
        
        MatrixXd LAV_calendar(MatrixXd calendar, num_vec t_group, double cov, num_vec rate, int t_start, int t_end, int up_type, char t, double eff_vac)
        {
            MatrixXd sero = MatrixXd::Zero(365, NoAgeG);
            MatrixXd outmat = MatrixXd::Zero(365, NoAgeG);

            boost::math::weibull_distribution<> W(2.42, 12.87); // delay between vaccination and immunity
            
            if (up_type == 1)
            {
                for (int i = 0; i < 365; i++)
                    for (int j = 0 ; j < NoAgeG; j++)
                        calendar(i%365,j) = rate[i]*t_group[j]*cov;
            }
            else
            {
                for (int i = t_start*7 ; i < t_end*7+1; i++)
                    for (int j = 0 ; j < NoAgeG; j++)
                        calendar(i%365,j) = t_group[j]*cov/30.0;
            }
            
            for(int i = 0; i < 365; i++)
                for(int j = 0; j < NoAgeG; j++)
                    for(int k = 1; k < 30; k++)
                        sero((i + k)%365, j) += pdf(W,k)*calendar((i)%365, j)*eff_vac;
            
            if (t == 'p')
                outmat = sero;
            else
                outmat = calendar;
            return outmat;
        }
        
        MatrixXd mAB_calendar(MatrixXd calendar, num_vec t_group, double cov, int t_start, int t_end, char t, double eff_mab)
        {
            MatrixXd sero = MatrixXd::Zero(365, NoAgeG);
            MatrixXd outmat = MatrixXd::Zero(365, NoAgeG);
            
            for (int i = t_start*7 ; i < t_end*7+1; i++)
                for (int j = 0 ; j < NoAgeG; j++)
                    calendar(i%365,j) = t_group[j]*cov/30.0;
            
            for (int i = t_start*7 ; i < t_end*7+1; i++)
                for (int j = 0 ; j < NoAgeG; j++)
                    sero(i%365,j) = calendar(i%365,j)*eff_mab;
            
            
            if (t == 'p')
                outmat = sero;
            else
                outmat = calendar;
            
            return outmat;
        }
        
        MatrixXd mAB_calendar_plus(MatrixXd calendar, num_vec t_group, double cov, int t_start, int t_end, char t, double eff_mab)
        {
            MatrixXd sero = MatrixXd::Zero(365, NoAgeG);
            MatrixXd outmat = MatrixXd::Zero(365, NoAgeG);
            
            for (int i = t_start*7 ; i < t_end*7+1; i++)
                for (int j = 0 ; j < NoAgeG; j++)
                    calendar(i%365,j) = t_group[j]*cov/30.0;
            
            
            for (int i = t_start*7 ; i < (t_start+4)*7+1; i++)
                for (int j = 0 ; j < NoAgeG; j++)
                    calendar(i%365,j) = cal::G_plus[j]*cov/30.0;
            
            for (int i = t_start*7 ; i < t_end*7+1; i++)
                for (int j = 0 ; j < NoAgeG; j++)
                    sero(i%365,j) = calendar(i%365,j)*eff_mab;
            
            if (t == 'p')
                outmat = sero;
            else
                outmat = calendar;
            
            return outmat;
        }
        
        MatrixXd mat_calendar(MatrixXd calendar, num_vec t_group, double cov, int t_start, int t_end, char t, double eff_mat, double eff_vac)
        {
            // vaccination can occur between 28-32 weeks gestational age for pregnant women
            MatrixXd sero = MatrixXd::Zero(365, NoAgeG);
            MatrixXd outmat = MatrixXd::Zero(365, NoAgeG);

            boost::math::weibull_distribution<> W(2.42, 12.87); // delay between vaccination and immunity

            for (int i = t_start*7 ; i < t_end*7+1; i++)
                for (int j = 0 ; j < NoAgeG; j++)
                    calendar(i%365,j) = t_group[j]/365.0;
            
            boost::math::uniform_distribution<> U(8*7, 12*7); // delay between vaccination of mother and birth of neonate

            for(int i = 0; i < 365; i++)
            {
                for(int j = 0; j < NoAgeG; j++)
                {
                    for(int k = 1; k < 31; k++)
                    {
                        sero((i + k)%365, j) += pdf(W,k)*calendar((i)%365, j)*eff_vac;
                    }
                }
                for(int k = 8*7; k < 12*7; k++)
                    sero((i + k)%365, 0) += pdf(U,k)*calendar((i)%365, 19)*eff_mat*365;
                
            }
            if (t == 'p')
                outmat = sero;
            else
                outmat = calendar;
            
            return outmat;
        }
    };
}


namespace interventions
{
    VectorXd no_doses(num_vec cov, num_vec start_w, num_vec end_w)
    {
        double pal_num = 2596;
        double prop_hr = 0.038;
        VectorXd no_dose;
        no_dose = VectorXd::Zero(14);
        no_dose(0) = pal_num*cov[1];                                                             // Pal
        no_dose(1) = popsize[0]*12*cov[2]*prop_hr*(end_w[2]-start_w[2])/52.0 + pal_num*cov[1];                    // HR, S
        no_dose(2) = popsize[0]*12*cov[3]*prop_hr*(end_w[3]-start_w[3])/52.0 + pal_num*cov[1];                    // HR, A
        no_dose(3) = popsize[0]*12*cov[4]*(end_w[4]-start_w[4])/52.0 + pal_num*cov[1];                            // HR, S
        no_dose(4) = popsize[0]*12*cov[5]*(end_w[5]-start_w[5])/52.0 + pal_num*cov[1];                            // HR, S
        no_dose(5) = popsize[0]*12*cov[6]*(end_w[6]-start_w[6])/52.0 + pal_num*cov[1];                            // mat, S
        no_dose(6) = popsize[0]*12*cov[7]*(end_w[7]-start_w[7])/52.0 + pal_num*cov[1];                            // mat, A
        no_dose(7) = popsize[0]*12*cov[8]*(end_w[8]-start_w[8])/52.0 + pal_num*cov[1];                            // inf, s
        no_dose(8) = popsize[0]*12*cov[9]*(end_w[9]-start_w[9])/52.0 + pal_num*cov[1];                            // inf, A
        no_dose(9) = (popsize[13]+popsize[14]+popsize[15])*cov[10] + pal_num*cov[1]; // 2-4
        no_dose(10) = popsize[16]*cov[11] + pal_num*cov[1];                           // 5-9
        no_dose(11) = (popsize[16]+popsize[17])*cov[12] + pal_num*cov[1];             // 5-14
        no_dose(12) = popsize[24]*cov[13] + pal_num*cov[1];                           // 75+
        no_dose(13) = (popsize[24]+popsize[23]) + pal_num*cov[1];             // 65+
        
        string outvacno = get_ll::dout +  "inter/";
        string File_vacno = outvacno + "vac_no.txt";      //Output MC
        ofstream File_vacno_out(File_vacno.c_str());
        File_vacno_out << no_dose;
        File_vacno_out.close();
        
        return no_dose;
    }
    
    num_vec initial_cond(const param::param_state_t& pars, double cov_c)
    {
        num_vec age_size_v = get_ll::popsize;
        num_vec comp_size_init;
        num_vec init_con_m = get_ll::initial_M(pars);
        VectorXd pA(NoAgeG);
        
        for (int a = 0; a < NoAgeG; a++)
            pA(a) = pars.mapofParam.find(pars.pA_vec[a])->second;
        
        double init_inf = pars.mapofParam.find("seed1")->second;
        double init_rec = pars.mapofParam.find("seed2")->second;
        double init_exp_time = 0.5;
        double s_prop, r_prop;
        double pI1; double pI2; double pI3; double pI4;
        
        double si = 1.0/(pars.mapofParam.find("si")->second);
        double g0 = 1.0/(pars.mapofParam.find("g0")->second);
        double g1 = 1.0/((pars.mapofParam.find("g0")->second)*(pars.mapofParam.find("g1")->second));
        double g2 = 1.0/((pars.mapofParam.find("g0")->second)*(pars.mapofParam.find("g1")->second)*(pars.mapofParam.find("g2")->second));
        double d1 = pars.mapofParam.find("d1")->second;
        double d2 = (pars.mapofParam.find("d1")->second)*(pars.mapofParam.find("d2")->second);
        double d3 = (pars.mapofParam.find("d1")->second)*(pars.mapofParam.find("d2")->second)*(pars.mapofParam.find("d3")->second);

        double a1, a2;
        /** Get ratio of number of each person in each exposure group per age group **/
        for(int a = 0; a < NoAgeG; a++)
        {
            if (a < 24)
            {
                a1 = ageGroup[a]; a2 = ageGroup[a+1];
            }
            
            else
            {
                a1 = ageGroup[a]; a2 = 90;
            }
            num_vec prop_init = get_ll::prop_init_ex(init_exp_time, a1, a2);
            pI1 = prop_init[0]; pI2 = prop_init[1]; pI3 = prop_init[2]; pI4 = prop_init[3];
            /** Get ratio of number of number of each person in each compoartment group **/
            double age_size = age_size_v[a] - init_con_m[a];
            
            for (int s = 0; s < 6; s++)
            {
                if (a < 12)
                {
                    if (s == 0 || s == 3)
                        s_prop = (0.0);
                    else if (s == 1 || s == 4)
                        s_prop = (cov_c);
                    else if (s == 2 || s == 5)
                        s_prop = (1-cov_c);
                    else
                        cout << "OOPS" << '\n';
                }
                else
                {
                    if (s == 0 || s == 3)
                        s_prop = p_mat[a]*(1-cov_c);
                    else if (s == 1 || s == 4)
                        s_prop = p_mat[a]*(cov_c);
                    else if (s == 2 || s == 5)
                        s_prop = 1-p_mat[a];
                    
                    else
                        cout << "OOPS" << '\n';
                }
                for (int r = 0; r < 3; r++)
                {
                    if (r == 0)
                        r_prop = pVHR[a];
                    else if (r == 1)
                        r_prop = pHR[a];
                    else if (r == 2)
                        r_prop = pLR[a];
                    else
                        cout << "ERROR" << endl;
                    
                    comp_size_init.push_back(r_prop*s_prop*init_con_m[a]); // Number  in M group
                    
                    comp_size_init.push_back(r_prop*s_prop*pI1*age_size*(1.0 - init_inf)*(1.0-init_rec));      //Sus
                    comp_size_init.push_back(r_prop*s_prop*pI1*age_size*init_inf*si/(si+g0));       //Exp
                    comp_size_init.push_back(r_prop*s_prop*pI1*age_size*init_inf*g0/(si+g0)*pA(a));
                    comp_size_init.push_back(r_prop*s_prop*pI1*age_size*init_inf*g0/(si+g0)*(1-pA(a)));   //Inf S
                    comp_size_init.push_back(r_prop*s_prop*pI1*age_size*(1.0 - init_inf)*init_rec);      //Rec
                    
                    comp_size_init.push_back(r_prop*s_prop*pI2*age_size*(1.0 - d1*init_inf)*(1.0-init_rec));      //Sus
                    comp_size_init.push_back(r_prop*s_prop*pI2*age_size*d1*init_inf*si/(si+g1));       //Exp
                    comp_size_init.push_back(r_prop*s_prop*pI2*age_size*d1*init_inf*g1/(si+g1)*pA(a));
                    comp_size_init.push_back(r_prop*s_prop*pI2*age_size*d1*init_inf*g1/(si+g1)*(1-pA(a)));   //Inf S
                    comp_size_init.push_back(r_prop*s_prop*pI2*age_size*(1.0 - d1*init_inf)*init_rec);      //Rec
                    
                    comp_size_init.push_back(r_prop*s_prop*pI3*age_size*(1.0 - d2*init_inf)*(1.0-init_rec));      //S/Sus
                    comp_size_init.push_back(r_prop*s_prop*pI3*age_size*d2*init_inf*si/(si+g2));       //Exp
                    comp_size_init.push_back(r_prop*s_prop*pI3*age_size*d2*init_inf*g2/(si+g2)*pA(a));
                    comp_size_init.push_back(r_prop*s_prop*pI3*age_size*d2*init_inf*g2/(si+g2)*(1-pA(a)));   //Inf S
                    comp_size_init.push_back(r_prop*s_prop*pI3*age_size*(1.0 - d2*init_inf)*init_rec);      //ec
                    
                    comp_size_init.push_back(r_prop*s_prop*pI4*age_size*(1.0 - d3*init_inf)*(1.0-init_rec));        //Sus
                    comp_size_init.push_back(r_prop*s_prop*pI4*age_size*d3*init_inf*si/(si+g2));       //Exp
                    comp_size_init.push_back(r_prop*s_prop*pI4*age_size*d3*init_inf*g2/(si+g2)*pA(a));
                    comp_size_init.push_back(r_prop*s_prop*pI4*age_size*d3*init_inf*g2/(si+g2)*(1-pA(a)));   //Inf S
                    comp_size_init.push_back(r_prop*s_prop*pI4*age_size*(1.0 - d3*init_inf)*init_rec);     //Rec
                    
                    comp_size_init.push_back(0);     //Rec
                    comp_size_init.push_back(0);     //Rec
                    comp_size_init.push_back(0);     //Rec

                }
            }
            for (int j = 0; j < 23; j++)
                comp_size_init.push_back(0.0);      //inf exposure 1
        }
        return comp_size_init;
        
    }
    
    
    class ODE_Desc_inter
    {
        param::param_state_t pars;
        MatrixXd vac_cal_mat;
        MatrixXd vac_cal_pal;
        MatrixXd vac_cal_init;
        MatrixXd cal_PVM;
        cal::Calendar_full cal_full;
        double cov_c;
        bool direct, PVM;
        double om_mab;
        double xi_boost;


    public:
        double pA[NoAgeG];
        double p_vul;
        
        ODE_Desc_inter(const param::param_state_t& pars_t, cal::Calendar_full& cal_full_i, bool direct_t, double cov_c_t, double om_mab_t, double xi_boost_t): pars(pars_t), cal_full(cal_full_i), direct(direct_t), cov_c(cov_c_t), om_mab(om_mab_t), xi_boost(xi_boost_t) {
            for (int a = 0; a < NoAgeG; a++)
                pA[a] = pars.mapofParam.find(pars.pA_vec[a])->second;
            
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
        double beta;
        num_vec VHR_g = {0.653012, 0.291842, 1., 1., 1., 0.0321774, 1., 1., 0,0,0,0,0,0,0,0,0,0,0};
        num_vec LR_g = {1.00126, 1.00168, 1., 1., 1., 1.00067, 1., 1., 1.00002, 1, 1, 1, 1, 1, 1, 1, 1};
        
        double M, S0, S1, S2, S3, SN, E0, E1, E2, E3, A0, A1, A2, A3, I0, I1, I2, I3, R0, R1, R2, R3, M_, S0_, S1_, S2_, S3_, E0_, E1_, E2_, E3_, A0_, A1_, A2_, A3_, I0_, I1_, I2_, I3_, R0_, R1_, R2_, R3_, N;
        double Mv, S0v, S1v, S2v, S3v, E0v, E1v, E2v, E3v, A0v, A1v, A2v, A3v, I0v, I1v, I2v, I3v, R0v, R1v, R2v, R3v, Mv_, S0v_, S1v_, S2v_, S3v_, E0v_, E1v_, E2v_, E3v_, A0v_, A1v_, A2v_, A3v_, I0v_, I1v_, I2v_, I3v_, R0v_, R1v_, R2v_, R3v_;
        
        double I_n[NoAgeG];     double I_c[NoAgeG];     double I_p[NoAgeG];
        double I_n_v[NoAgeG];   double I_c_v[NoAgeG];   double I_p_v[NoAgeG];
        
        double N_tot[NoAgeG];
        double N_tot_n[NoAgeG];         double N_tot_c[NoAgeG];         double N_tot_p[NoAgeG];
        double N_tot_n_v[NoAgeG];       double N_tot_c_v[NoAgeG];       double N_tot_p_v[NoAgeG];
        double N_tot_n_v_inv, N_tot_p_v_inv, N_tot_c_v_inv, N_tot_n_inv, N_tot_p_inv, N_tot_c_inv;
        
        double prop_n[NoAgeG];  double prop_c[NoAgeG];  double prop_p[NoAgeG];
        double prop_nv[NoAgeG]; double prop_cv[NoAgeG]; double prop_pv[NoAgeG];
        double prop_empty[NoAgeG];
        
        double PS[24];
        
        double I_temp_p, I_temp_c, I_temp_n, I_temp_p_v, I_temp_c_v, I_temp_n_v;
        double x_tot, x_tot_1, x_tot_2;
        //Calendar
        MatrixXd cal_pal = cal_full.cal_pal;
        MatrixXd cal_mAB_VHR = cal_full.cal_mAB_VHR;
        MatrixXd cal_mAB_HR = cal_full.cal_mAB_HR;
        MatrixXd cal_mAB_LR = cal_full.cal_mAB_LR;
        MatrixXd cal_LAV_HR = cal_full.cal_LAV_HR;
        MatrixXd cal_LAV_LR = cal_full.cal_LAV_LR;
        MatrixXd cal_mat_LR = cal_full.cal_mat_LR;
        
        MatrixXd cal_pal_dose = cal_full.cal_pal_dose;
        MatrixXd cal_mAB_VHR_dose = cal_full.cal_mAB_VHR_dose;
        MatrixXd cal_mAB_HR_dose = cal_full.cal_mAB_HR_dose;
        MatrixXd cal_mAB_LR_dose = cal_full.cal_mAB_LR_dose;
        MatrixXd cal_LAV_HR_dose = cal_full.cal_LAV_HR_dose;
        MatrixXd cal_LAV_LR_dose = cal_full.cal_LAV_LR_dose;
        MatrixXd vac_cal_dose = cal_full.cal_mat_LR_dose;
        double xi_b;
        
        /** Description of ODE equations **/
        void operator() (  num_vec &x , num_vec &dxdt , const double  t  )
        {
            int ag = 455;
            int sg = 72;
            int rg = 24;
            /*/////////////*/
            /* CALENDARS /*/
            /*/////////////*/
            MatrixXd vac_cal_vhr; MatrixXd vac_cal;
            if (t < 365)
            {
                cal_pal = MatrixXd::Zero(365, NoAgeG);
                vac_cal = MatrixXd::Zero(365, NoAgeG);
                cal_mAB_VHR = MatrixXd::Zero(365, NoAgeG);
                cal_mAB_HR = MatrixXd::Zero(365, NoAgeG);
                cal_mAB_LR = MatrixXd::Zero(365, NoAgeG);
                cal_LAV_HR = MatrixXd::Zero(365, NoAgeG);
                cal_LAV_LR = MatrixXd::Zero(365, NoAgeG);

            }
            else{
                cal_pal = cal_full.cal_pal;
                vac_cal = cal_full.cal_mat_LR;
                cal_mAB_VHR = cal_full.cal_mAB_VHR;
                cal_mAB_HR = cal_full.cal_mAB_HR;
                cal_mAB_LR = cal_full.cal_mAB_LR;
                cal_LAV_HR = cal_full.cal_LAV_HR;
                cal_LAV_LR = cal_full.cal_LAV_LR;
            }
            
            /*/////////////*/
            /* POPULATIONS /*/
            /*/////////////*/
            //Find sum of every group
            double tot_coc = 0, tot_coc_v = 0, tot_inf = 0, tot_inf_v = 0;
            
            for(int a = 0; a < NoAgeG; a++)
            {
                double r_prop = 0;
                double tot_temp_n = 0, tot_temp_c = 0, tot_temp_p = 0, tot_temp_nv = 0, tot_temp_cv = 0, tot_temp_pv = 0;
                double tot_temp_vhr = 0, tot_temp_vhrv = 0, tot_temp_hr = 0, tot_temp_hrv = 0, tot_temp_lr = 0, tot_temp_lrv = 0;

                for (int j = 0; j < 24; j++){
                    for (int r = 0; r < 3; r++){
                        tot_temp_p +=  x[a*ag+0*sg+r*rg+j];        tot_temp_c += x[a*ag+1*sg+r*rg+j];            tot_temp_n += x[a*ag+2*sg+r*rg+j];
                        tot_temp_pv += x[a*ag+3*sg+r*rg+j];       tot_temp_cv += x[a*ag+4*sg+r*rg+j];           tot_temp_nv += x[a*ag+5*sg+r*rg+j];
                    }
                }
                
                for (int j = 0; j < 24; j++){
                    for (int s = 0; s < 3; s++){
                        tot_temp_vhr +=  x[a*ag+s*sg+0*rg+j];        tot_temp_hr += x[a*ag+s*sg+1*rg+j];            tot_temp_lr += x[a*ag+s*sg+2*rg+j];
                        tot_temp_vhrv += x[a*ag+(s+3)*sg+0*rg+j];       tot_temp_hrv += x[a*ag+(s+3)*sg+1*rg+j];           tot_temp_lrv += x[a*ag+(s+3)*sg+2*rg+j];
                    }
                }
                
                if (a < 12){
                    tot_coc += tot_temp_c;                  tot_coc_v += tot_temp_cv;
                    tot_inf += tot_temp_c + tot_temp_n;     tot_inf_v += tot_temp_cv + tot_temp_nv;
                }
            
                N_tot_n[a] = tot_temp_n;    N_tot_c[a] = tot_temp_c;    N_tot_p[a] = tot_temp_p;
                N_tot_n_v[a] = tot_temp_nv; N_tot_c_v[a] = tot_temp_cv; N_tot_p_v[a] = tot_temp_pv;

            }
            
            /*/////////////*/
            /* 2. Force of Infection /*/
            /*/////////////*/
            double phi_c = cov_c;
            double phi_c_v = cov_c;
            
            //  2.1 Get contact matrices
            vector2D cnt_matrix_cwc_p = manip_contact::get_cwc(phi_c,'p');
            vector2D cnt_matrix_cwp_p = manip_contact::get_cwp(phi_c,'p');
            vector2D cnt_matrix_cwn_p = manip_contact::get_cwn(phi_c,'p');
            vector2D cnt_matrix_pwc_p = manip_contact::get_pwc(phi_c,'p');
            vector2D cnt_matrix_pwp_p = manip_contact::get_pwp(phi_c,'p');
            vector2D cnt_matrix_pwn_p = manip_contact::get_pwn(phi_c,'p');
            vector2D cnt_matrix_nwc_p = manip_contact::get_nwc(phi_c,'p');
            vector2D cnt_matrix_nwp_p = manip_contact::get_nwp(phi_c,'p');
            vector2D cnt_matrix_nwn_p = manip_contact::get_nwn(phi_c,'p');
            
            vector2D cnt_matrix_cwc_c = manip_contact::get_cwc(phi_c,'c');
            vector2D cnt_matrix_cwp_c = manip_contact::get_cwp(phi_c,'c');
            vector2D cnt_matrix_cwn_c = manip_contact::get_cwn(phi_c,'c');
            vector2D cnt_matrix_pwc_c = manip_contact::get_pwc(phi_c,'c');
            vector2D cnt_matrix_pwp_c = manip_contact::get_pwp(phi_c,'c');
            vector2D cnt_matrix_pwn_c = manip_contact::get_pwn(phi_c,'c');
            vector2D cnt_matrix_nwc_c = manip_contact::get_nwc(phi_c,'c');
            vector2D cnt_matrix_nwp_c = manip_contact::get_nwp(phi_c,'c');
            vector2D cnt_matrix_nwn_c = manip_contact::get_nwn(phi_c,'c');

            vector2D cnt_matrix_cwc_p_v = manip_contact::get_cwc(phi_c_v,'p');
            vector2D cnt_matrix_cwp_p_v = manip_contact::get_cwp(phi_c_v,'p');
            vector2D cnt_matrix_cwn_p_v = manip_contact::get_cwn(phi_c_v,'p');
            vector2D cnt_matrix_pwc_p_v = manip_contact::get_pwc(phi_c_v,'p');
            vector2D cnt_matrix_pwp_p_v = manip_contact::get_pwp(phi_c_v,'p');
            vector2D cnt_matrix_pwn_p_v = manip_contact::get_pwn(phi_c_v,'p');
            vector2D cnt_matrix_nwc_p_v = manip_contact::get_nwc(phi_c_v,'p');
            vector2D cnt_matrix_nwp_p_v = manip_contact::get_nwp(phi_c_v,'p');
            vector2D cnt_matrix_nwn_p_v = manip_contact::get_nwn(phi_c_v,'p');
            
            vector2D cnt_matrix_cwc_c_v = manip_contact::get_cwc(phi_c_v,'c');
            vector2D cnt_matrix_cwp_c_v = manip_contact::get_cwp(phi_c_v,'c');
            vector2D cnt_matrix_cwn_c_v = manip_contact::get_cwn(phi_c_v,'c');
            vector2D cnt_matrix_pwc_c_v = manip_contact::get_pwc(phi_c_v,'c');
            vector2D cnt_matrix_pwp_c_v = manip_contact::get_pwp(phi_c_v,'c');
            vector2D cnt_matrix_pwn_c_v = manip_contact::get_pwn(phi_c_v,'c');
            vector2D cnt_matrix_nwc_c_v = manip_contact::get_nwc(phi_c_v,'c');
            vector2D cnt_matrix_nwp_c_v = manip_contact::get_nwp(phi_c_v,'c');
            vector2D cnt_matrix_nwn_c_v = manip_contact::get_nwn(phi_c_v,'c');
            
            //  2.2 Seasonal forcing
            //beta = (1 + b1*cos((t/365.0-phi)*2*PI));
            double t1 = (int)t%365;
            beta = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))));

            /*/////////////*/
            /*  3. DYNAMIC MATERNAL PROTECTION /*/
            /*/////////////*/
           // num_vec pos_wcb_ag = {18,21};
            double sum_wcb = 0;
            double sum_wcb_v = 0;
            double CB2_temp, CB2_temp_v;
            CB2_temp = CB2_temp_v = 0;
            for(int a = 18; a < 21 ; a++)
            {
                for (int r = 0; r < 3; r++)
                {
                    CB2_temp += (x[a*ag+2*sg+r*rg+1] + x[a*ag+2*sg+r*rg+6] + x[a*ag+2*sg+r*rg+11] + x[a*ag+2*sg+r*rg+16] + x[a*ag+2*sg+r*rg+2] + x[a*ag+2*sg+r*rg+7] + x[a*ag+2*sg+r*rg+12] + x[a*ag+2*sg+r*rg+17]);
                    CB2_temp_v += (x[a*ag+5*sg+r*rg+1] + x[a*ag+5*sg+r*rg+6] + x[a*ag+5*sg+r*rg+11] + x[a*ag+5*sg+r*rg+16] + x[a*ag+5*sg+r*rg+2] + x[a*ag+5*sg+r*rg+7] + x[a*ag+5*sg+r*rg+12] + x[a*ag+5*sg+r*rg+17]);
                }
            }
            sum_wcb = CB2_temp/((double)N_tot_n[18]+(double)N_tot_n[19]+(double)N_tot_n[20]);
            sum_wcb_v = (CB2_temp_v)/((double)N_tot_n[18]+(double)N_tot_n[19]+(double)N_tot_n[20]);
            
            double tot_C = N_tot_c[18] + N_tot_c[19] + N_tot_c[20];
            double tot_C_inv;
            if (tot_C < 1.0)
                tot_C_inv = 0;
            else
                tot_C_inv = 1.0/tot_C;
            /*/////////////*/
            /* 4. ODES /*/
            /*/////////////*/
            //4.1 Age groups
            double protectpal;
            double protectmabs;
            double protectLAV;
            double protectmat;

            for(int a = 0; a < NoAgeG; a++)
            {
                double r_prop = 0;
                protectpal = protectmabs = protectLAV = protectmat = 0;
                I_temp_c = 0.0;    I_temp_p = 0.0; I_temp_n = 0.0;
                for (int k = 0; k < NoAgeG; k++)
                {
                    
                    if (N_tot_n[k] < 0.1){N_tot_n_inv = 0;}else{N_tot_n_inv=1.0/N_tot_n[k];}
                    if (N_tot_c[k] < 0.1){N_tot_c_inv = 0;}else{N_tot_c_inv=1.0/N_tot_c[k];}
                    if (N_tot_p[k] < 0.1){N_tot_p_inv = 0;}else{N_tot_p_inv=1.0/N_tot_p[k];}
                    
                    for (int r = 0; r < 3; r++)
                    {
                        I_temp_p += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_pwp_p[a][k]+qc*cnt_matrix_pwp_c[a][k]))*N_tot_p_inv +(x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_pwc_p[a][k]+qc*cnt_matrix_pwc_c[a][k]))*N_tot_c_inv +(x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_pwn_p[a][k]+qc*cnt_matrix_pwn_c[a][k]))*N_tot_n_inv;
                        
                        I_temp_c += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_cwp_p[a][k]+qc*cnt_matrix_cwp_c[a][k]))*N_tot_p_inv +(x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_cwc_p[a][k]+qc*cnt_matrix_cwc_c[a][k]))*N_tot_c_inv +(x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_cwn_p[a][k]+qc*cnt_matrix_cwn_c[a][k]))*N_tot_n_inv;
                        
                        I_temp_n += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_nwp_p[a][k]+qc*cnt_matrix_nwp_c[a][k]))*N_tot_p_inv +(x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_nwc_p[a][k]+qc*cnt_matrix_nwc_c[a][k]))*N_tot_c_inv +(x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_nwn_p[a][k]+qc*cnt_matrix_nwn_c[a][k]))*N_tot_n_inv;
                    }

                }

                
                I_temp_c_v = 0.0;  I_temp_p_v = 0.0; I_temp_n_v = 0.0;
                for (int k = 0; k < NoAgeG ; k++)
                {
                    if (N_tot_n_v[k] < 1){N_tot_n_v_inv = 0;}else{N_tot_n_v_inv=1.0/N_tot_n_v[k];}
                    if (N_tot_c_v[k] < 1){N_tot_c_v_inv = 0;}else{N_tot_c_v_inv=1.0/N_tot_c_v[k];}
                    if (N_tot_p_v[k] < 1){N_tot_p_v_inv = 0;}else{N_tot_p_v_inv=1.0/N_tot_p_v[k];}
                    
                    for (int r = 0; r < 3; r++)
                    {
                        I_temp_p_v += (x[k*ag+3*sg+r*rg+3]*alpha_i+x[k*ag+3*sg+r*rg+4]+a1*(x[k*ag+3*sg+r*rg+8]*alpha_i+x[k*ag+3*sg+r*rg+9])+a2*(x[k*ag+3*sg+r*rg+13]*alpha_i+x[k*ag+3*sg+r*rg+14])+a3*(x[k*ag+3*sg+r*rg+18]*alpha_i+x[k*ag+3*sg+r*rg+19]))*(qp*(cnt_matrix_pwp_p_v[a][k]+qc*cnt_matrix_pwp_c_v[a][k]))*N_tot_p_v_inv + (x[k*ag+4*sg+r*rg+3]*alpha_i+x[k*ag+4*sg+r*rg+4]+a1*(x[k*ag+4*sg+r*rg+8]*alpha_i+x[k*ag+4*sg+r*rg+9])+a2*(x[k*ag+4*sg+r*rg+13]*alpha_i+x[k*ag+4*sg+r*rg+14])+a3*(x[k*ag+4*sg+r*rg+18]*alpha_i+x[k*ag+4*sg+r*rg+19]))*(qp*(cnt_matrix_pwc_p_v[a][k]+qc*cnt_matrix_pwc_c_v[a][k]))*N_tot_c_v_inv  + (x[k*ag+5*sg+r*rg+3]*alpha_i+x[k*ag+5*sg+r*rg+4]+a1*(x[k*ag+5*sg+r*rg+8]*alpha_i+x[k*ag+5*sg+r*rg+9])+a2*(x[k*ag+5*sg+r*rg+13]*alpha_i+x[k*ag+5*sg+r*rg+14])+a3*(x[k*ag+5*sg+r*rg+18]*alpha_i+x[k*ag+5*sg+r*rg+19]))*(qp*(cnt_matrix_pwn_p_v[a][k]+qc*cnt_matrix_pwn_c_v[a][k]))*N_tot_n_v_inv;
                        
                        I_temp_c_v += (x[k*ag+3*sg+r*rg+3]*alpha_i+x[k*ag+3*sg+r*rg+4]+a1*(x[k*ag+3*sg+r*rg+8]*alpha_i+x[k*ag+3*sg+r*rg+9])+a2*(x[k*ag+3*sg+r*rg+13]*alpha_i+x[k*ag+3*sg+r*rg+14])+a3*(x[k*ag+3*sg+r*rg+18]*alpha_i+x[k*ag+3*sg+r*rg+19]))*(qp*(cnt_matrix_cwp_p_v[a][k]+qc*cnt_matrix_cwp_c_v[a][k]))*N_tot_p_v_inv +(x[k*ag+4*sg+r*rg+3]*alpha_i+x[k*ag+4*sg+r*rg+4]+a1*(x[k*ag+4*sg+r*rg+8]*alpha_i+x[k*ag+4*sg+r*rg+9])+a2*(x[k*ag+4*sg+r*rg+13]*alpha_i+x[k*ag+4*sg+r*rg+14])+a3*(x[k*ag+4*sg+r*rg+18]*alpha_i+x[k*ag+4*sg+r*rg+19]))*(qp*(cnt_matrix_cwc_p_v[a][k]+qc*cnt_matrix_cwc_c_v[a][k]))*N_tot_c_v_inv +(x[k*ag+5*sg+r*rg+3]*alpha_i+x[k*ag+5*sg+r*rg+4]+a1*(x[k*ag+5*sg+r*rg+8]*alpha_i+x[k*ag+5*sg+r*rg+9])+a2*(x[k*ag+5*sg+r*rg+13]*alpha_i+x[k*ag+5*sg+r*rg+14])+a3*(x[k*ag+5*sg+r*rg+18]*alpha_i+x[k*ag+5*sg+r*rg+19]))*(qp*(cnt_matrix_cwn_p_v[a][k]+qc*cnt_matrix_cwn_c_v[a][k]))*N_tot_n_v_inv;
                        
                        I_temp_n_v += (x[k*ag+3*sg+r*rg+3]*alpha_i+x[k*ag+3*sg+r*rg+4]+a1*(x[k*ag+3*sg+r*rg+8]*alpha_i+x[k*ag+3*sg+r*rg+9])+a2*(x[k*ag+3*sg+r*rg+13]*alpha_i+x[k*ag+3*sg+r*rg+14])+a3*(x[k*ag+3*sg+r*rg+18]*alpha_i+x[k*ag+3*sg+r*rg+19]))*(qp*(cnt_matrix_nwp_p_v[a][k]+qc*cnt_matrix_nwp_c_v[a][k]))*N_tot_p_v_inv + (x[k*ag+4*sg+r*rg+3]*alpha_i+x[k*ag+4*sg+r*rg+4]+a1*(x[k*ag+4*sg+r*rg+8]*alpha_i+x[k*ag+4*sg+r*rg+9])+a2*(x[k*ag+4*sg+r*rg+13]*alpha_i+x[k*ag+4*sg+r*rg+14])+a3*(x[k*ag+4*sg+r*rg+18]*alpha_i+x[k*ag+4*sg+r*rg+19]))*(qp*(cnt_matrix_nwc_p_v[a][k]+qc*cnt_matrix_nwc_c_v[a][k]))*N_tot_c_v_inv + (x[k*ag+5*sg+r*rg+3]*alpha_i+x[k*ag+5*sg+r*rg+4]+a1*(x[k*ag+5*sg+r*rg+8]*alpha_i+x[k*ag+5*sg+r*rg+9])+a2*(x[k*ag+5*sg+r*rg+13]*alpha_i+x[k*ag+5*sg+r*rg+14])+a3*(x[k*ag+5*sg+r*rg+18]*alpha_i+x[k*ag+5*sg+r*rg+19]))*(qp*(cnt_matrix_nwn_p_v[a][k]+qc*cnt_matrix_nwn_c_v[a][k]))*N_tot_n_v_inv;
                    }
                }

                if (direct)
                {
                    I_temp_p_v = I_temp_p; I_temp_c_v = I_temp_c; I_temp_n_v = I_temp_n;
                }
                
                int pj = max(ag*(a-1),0);
                int cj = max(ag*(a),0);
                double muBp, muBc, muBn, muBpv, muBcv, muBnv, mu_mat, cl, Icp, Icc, up, kpc, kpd, pro, pro_v_in, pro_v_out, rp;
                double p_vulp, p_vulc, p_vuln, p_vulpv, p_vulcv, p_vulnv;
                double xi_bp, xi_bc, xi_bn, xi_bpv, xi_bcv, xi_bnv;

                double ej1  = eta[a+1];
                double ej   = eta[a];
                double lossP, lossMS0, lossMS1;
                
                // Birth rate into each social group
                if (a == 0)
                {
                    p_vulp = 0; p_vulc = sum_wcb; p_vuln = sum_wcb; p_vulpv = 0; p_vulcv = sum_wcb_v; p_vulnv = sum_wcb_v;
                    
                    muBp = 0;    muBc = muB*(cov_c);   muBn = muB*(1-cov_c);
                    muBpv = 0;  muBcv = muB*(cov_c);  muBnv = muB*(1-cov_c);
                    xi_bp = 1;    xi_bc = 1;   xi_bn = 1;
                    xi_bpv = 1;   xi_bcv = xi_boost;  xi_bnv = 1;
                }
                else
                {
                    p_vulp = 0; p_vulc = 0; p_vuln = 0; p_vulpv = 0; p_vulcv = 0; p_vulnv = 0;
                    muBp = 0; muBc = 0; muBn = 0;
                    muBpv = 0; muBcv = 0; muBnv = 0;
                    xi_bp = 1;    xi_bc = 1;   xi_bn = 1;
                    xi_bpv = 1;   xi_bcv = 1;  xi_bnv = 1;
                }
                
                //4.2 Social groups
                for(int s = 0; s < 6; s++)
                {
                    
                    kpd = kpc = pro = pro_v_in = pro_v_out = 0;
                    cl = 0;     Icp = 1;    Icc = 0;    up = 0;
                    
                    double u, In, mu, vac_c_o, vac_c_i;
                    int kp;
                    if (a < 12)
                    {
                        if (s == 0)
                        {   kp = 0*sg;       u = 0;             In = I_temp_p;   mu = muBp;  mu_mat = 0; p_vul = p_vulp; xi_b = xi_bp;
                        }
                        else if (s == 1)
                        {   kp = 1*sg;       u = (cov_c);       In = I_temp_c;   mu = muBc;  mu_mat = 0; p_vul = p_vulc; xi_b = xi_bc;
                        }
                        else if (s == 2)
                        {   kp = 2*sg;       u = (1-cov_c);     In = I_temp_n;   mu = muBn;  mu_mat = 0; p_vul = p_vuln; xi_b = xi_bn;
                        }
                        else if (s == 3)
                        {   kp = 3*sg;       u = 0;             In = I_temp_p_v; mu = muBpv; mu_mat = 0; p_vul = p_vulpv; xi_b = xi_bpv;
                        }
                        else if (s == 4)
                        {   kp = 4*sg;       u = (cov_c);       In = I_temp_c_v; mu = muBcv*(1-vac_cal(t1,0)); mu_mat = muBcv*vac_cal(t1,0); p_vul = p_vulcv;  xi_b = xi_bcv;
                        }
                        else
                        {   kp = 5*sg;       u = (1-cov_c);     In = I_temp_n_v; mu = muBnv; mu_mat = 0; p_vul = p_vulnv;  xi_b = xi_bnv;
                        }
                    }
                    else
                    {
                        if (s == 0)
                        {   kp = 0*sg;       u = p_mat[a]*(1-cov_c);        In = I_temp_p;   mu = mu_mat = 0;  p_vul = p_vulp;
                        }
                        else if (s == 1)
                        {   kp = 1*sg;       u = p_mat[a]*(cov_c);          In = I_temp_c;   mu = mu_mat = 0;  p_vul = p_vulc;
                        }
                        else if (s == 2)
                        {   kp = 2*sg;       u = (1-p_mat[a]);              In = I_temp_n;   mu = mu_mat = 0;  p_vul = p_vuln;
                        }
                        else if (s == 3)
                        {   kp = 3*sg;       u = p_mat[a]*(1-cov_c);        In = I_temp_p_v; mu = mu_mat = 0; p_vul = p_vulpv;
                        }
                        else if (s == 4)
                        {   kp = 4*sg;       u = p_mat[a]*(cov_c);          In = I_temp_c_v; mu = mu_mat = 0; p_vul = p_vulcv;
                        }
                        else
                        {   kp = 5*sg;       u = (1-p_mat[a]);             In = I_temp_n_v; mu = mu_mat = 0; p_vul = p_vulnv;
                        }
                    }
                    //4.3 Risk groups
                    for (int r = 0; r < 3; r++)
                    {
                        double PST = 0;
                        x_tot = x_tot_1 = x_tot_2 = 0;
                        if (r==0){
                            rp = pVHR[a];
                        }
                        else if (r==1){
                            rp = pHR[a];}
                        else if (r==2){
                            rp = pLR[a];}
                        else
                            cout << "ERROR" << endl;
                        
                        if (r == 2)
                        {
                            if (a==18 || a==19 || a==20)
                            {
                                if (a==18)
                                    up = u18p;
                                else if (a==19)
                                    up = u19p;
                                else if (a==20)
                                    up = u20p;
                                else
                                    cout << "ERROR" << endl;
                            }
                        }
                    
                        if (s < 3)
                        {
                            if (a < 12)
                            {
                                /*if (a > 0)
                                    rp = 1;
                                
                                if (r == 0)
                                    r_prop = VHR_g[a];
                                else if (r == 2)
                                    r_prop = LR_g[a];
                                else
                                    r_prop = 1;
                                
                                for (int i = 0; i < 24; i++)
                                    PS[i] = (x[pj+0*sg+r*rg+i]+x[pj+1*sg+r*rg+i]+x[pj+2*sg+r*rg+i])*r_prop;*/
                                for (int i = 0; i < 21; i++)
                                    PS[i] = (x[pj+0*sg+0*rg+i]+x[pj+0*sg+1*rg+i]+x[pj+0*sg+2*rg+i] + x[pj+1*sg+0*rg+i]+x[pj+1*sg+1*rg+i]+x[pj+1*sg+2*rg+i] + x[pj+2*sg+0*rg+i] + x[pj+2*sg+1*rg+i]+x[pj+2*sg+2*rg+i]);
                                
                                for (int i = 21; i < 24; i++)
                                    PS[i] = x[pj+0*sg+r*rg+i] + x[pj+1*sg+r*rg+i] + x[pj+2*sg+r*rg+i];
                            }
                            else
                            {
                                for (int i = 0; i < 24; i++)
                                    PS[i] = (x[pj+0*sg+0*rg+i]+x[pj+0*sg+1*rg+i]+x[pj+0*sg+2*rg+i] + x[pj+1*sg+0*rg+i]+x[pj+1*sg+1*rg+i]+x[pj+1*sg+2*rg+i] + x[pj+2*sg+0*rg+i] + x[pj+2*sg+1*rg+i]+x[pj+2*sg+2*rg+i]);
                            }
                        }
                        else
                        {
                            if (a < 12)
                            {
                                /*if (a > 0)
                                    rp = 1;
                                
                                if (r == 3)
                                    r_prop = VHR_g[a];
                                else if (r == 5)
                                    r_prop = LR_g[a];
                                else
                                    r_prop = 1;
                                
                                for (int i = 0; i < 24; i++)
                                    PS[i] = (x[pj+3*sg+r*rg+i]+x[pj+4*sg+r*rg+i]+x[pj+5*sg+r*rg+i])*r_prop;
                                 */
                                for (int i = 0; i < 21; i++)
                                    PS[i] = x[pj+3*sg+0*rg+i] + x[pj+3*sg+1*rg+i] + x[pj+3*sg+2*rg+i] + x[pj+4*sg+0*rg+i] + x[pj+4*sg+1*rg+i] + x[pj+4*sg+2*rg+i] + x[pj+5*sg+0*rg+i] + x[pj+5*sg+1*rg+i] + x[pj+5*sg+2*rg+i];
                                
                                for (int i = 21; i < 24; i++)
                                    PS[i] = x[pj+3*sg+r*rg+i] + x[pj+4*sg+r*rg+i] + x[pj+5*sg+r*rg+i];
       
                            }
                            else
                            {
                                for (int i = 0; i < 24; i++)
                                    PS[i] = x[pj+3*sg+0*rg+i] + x[pj+3*sg+1*rg+i] + x[pj+3*sg+2*rg+i] + x[pj+4*sg+0*rg+i] + x[pj+4*sg+1*rg+i] + x[pj+4*sg+2*rg+i] + x[pj+5*sg+0*rg+i] + x[pj+5*sg+1*rg+i] + x[pj+5*sg+2*rg+i];
                            }
                        }
                        
                        double cpmu, cpo, cMmu, cMo, cLo, cMa, cpmu_dose, cpo_dose, cMmu_dose, cMo_dose, cLo_dose, cMa_dose;

                        int p = a*ag + s*sg + r*rg;
                        int q = a*ag + kpc + r*rg;
                        int o = 0;
                        if (s>3)
                            o = a*ag + (s-3)*sg + r*rg;
                        
                        if (s < 3)
                        {
                            cpmu = 0; cpmu_dose = 0;
                            cpo = 0; cpo_dose = 0;
                            cMmu = 0; cMmu_dose = 0;
                            cMo = 0; cMo_dose = 0;
                            cLo = 0; cLo_dose = 0;
                            cMa = 0; cMa_dose = 0;
                            SN = 1;
                            lossMS1 = lossMS0 = lossP = 0;

                        }
                        else
                        {
                                                        //if (a > 4){
                            lossP = x[a*ag + s*sg + r*rg + 21]*(1.0/60.0);
                            lossMS0 = x[a*ag + s*sg + r*rg + 22]*(om_mab);
                            lossMS1 = x[a*ag + s*sg + r*rg + 23]*(om_mab);
                            
                            cpmu = 0; cpmu_dose = 0;
                            cpo = 0; cpo_dose = 0;
                            cMmu = 0; cMmu_dose = 0;
                            cMo = 0; cMo_dose = 0;
                            cLo = 0; cLo_dose = 0;
                            cMa = 0; cMa_dose = 0;

                            if (a == 0){
                                if (r == 0){
                                    cpo = cal_pal(t1, 0);
                                    cMo = cal_mAB_VHR(t1, 0);
                                    cpo_dose = cal_pal_dose(t1, 0);
                                    cMo_dose = cal_mAB_VHR_dose(t1, 0);}
                                else if (r == 1){
                                    cMo = cal_mAB_HR(t1, 0);
                                    cMo_dose = cal_mAB_HR_dose(t1, 0);
                                }
                                else{
                                    cMo = cal_mAB_LR(t1, 0);
                                    cMo_dose = cal_mAB_LR_dose(t1, 0);
                                }
                            }
                            else if (a > 0)
                            {
                                if (r == 0){
                                    cpo = cal_pal(t1, a);
                                    cpo_dose = cal_pal_dose(t1, a);
                                    cMo = cal_mAB_VHR(t1, a);
                                    cMo_dose = cal_mAB_VHR_dose(t1, a);
                                    S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];
                                    SN = 1;
                                }
                                else if (r == 1){
                                    cMo = cal_mAB_HR(t1, a);
                                    cLo = cal_LAV_HR(t1, a);
                                    cMo_dose = cal_mAB_HR_dose(t1, a);
                                    cLo_dose = cal_LAV_HR_dose(t1, a);
                                    S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];
                                    SN = 1;
                                }
                                else
                                {
                                    cMo = cal_mAB_LR(t1, a);
                                    cMo_dose = cal_mAB_LR_dose(t1, a);
                                    if (s==4)
                                    {
                                        cMa = vac_cal(t1, a);
                                        cMa_dose = vac_cal_dose(t1, a);
                                        S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];
                                    }
                                    else
                                    {
                                        cMa = 0;
                                        cMa_dose = 0;
                                        SN = 1;
                                        cLo = cal_LAV_LR(t1, a);
                                        cLo_dose = cal_LAV_LR_dose(t1, a);
                                        S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];
                                    }
                                }
                            }
                            else
                            {
                                
                            }
                        }
                        
                        for (int i = 0; i < 21; i++)
                        {
                            x_tot += x[o+i];
                            PST += PS[i];
                            if (PST < 1)
                                PST = 1;
                        }
                        for (int i = 0; i < 2; i++)
                            x_tot_1 += x[o+i];
                        
                        for (int i = 2; i < 6; i++)
                            x_tot_2 += x[o+i];
                        
                        
                        protectpal += cpo_dose*x_tot;
                        protectmabs += cMo_dose*x_tot;
                        protectLAV += cLo_dose*x_tot;
                        protectmat += cMa_dose*x_tot;
                        
                        dxdt[p+0] = (1.0-p_vul)*mu*rp + mu_mat*rp - x[p+0]*xi*xi_b - (x[p+0])*ej1 + PS[0]*ej*rp*u - x[o+0]*(cMo) - x[o+0]*(cpo);
                        
                        dxdt[p+1] = p_vul*mu*rp  + x[p+0]*xi*xi_b + lossP + lossMS0 - x[p+1]*In*beta - (x[p+1])*ej1 + PS[1]*ej*rp*u - x[o+1]*(cMo) - x[o+1]*(cpo) - cLo*S0 - cMa*S0;
                        dxdt[p+2] = x[p+1]*In*beta                 - x[p+2]*si          - (x[p+2])*ej1 + PS[2]*ej*rp*u - x[o+2]*(cMo) - x[o+2]*(cpo);
                        dxdt[p+3] = x[p+2]*si*pA[a]                - x[p+3]*ga0*rho     - (x[p+3])*ej1 + PS[3]*ej*rp*u - x[o+3]*(cMo) - x[o+3]*(cpo);
                        dxdt[p+4] = x[p+2]*si*(1.0-pA[a])          - x[p+4]*ga0         - (x[p+4])*ej1 + PS[4]*ej*rp*u - x[o+4]*(cMo) - x[o+4]*(cpo);
                        dxdt[p+5] = x[p+4]*ga0 + x[p+3]*ga0*rho    - x[p+5]*om          - (x[p+5])*ej1 + PS[5]*ej*rp*u - x[o+5]*(cMo) - x[o+5]*(cpo) + cLo*S0 + cMa*S0;
                        
                        dxdt[p+6] = x[p+5]*om                      - d1*x[p+6]*In*beta + lossMS1  - (x[p+6])*ej1 + PS[6]*ej*rp*u - x[o+6]*(cMo) - x[o+6]*(cpo) - cLo*S1 - cMa*S1;
                        dxdt[p+7] = d1*x[p+6]*In*beta              - x[p+7]*si          - (x[p+7])*ej1 + PS[7]*ej*rp*u - x[o+7]*0 - x[o+7]*(cpo);
                        dxdt[p+8] = x[p+7]*si*pA[a]                - x[p+8]*ga1*rho     - (x[p+8])*ej1 + PS[8]*ej*rp*u - x[o+8]*0 - x[o+8]*(cpo);
                        dxdt[p+9] = x[p+7]*si*(1.0-pA[a])          - x[p+9]*ga1         - (x[p+9])*ej1 + PS[9]*ej*rp*u - x[o+9]*0 - x[o+9]*(cpo);
                        dxdt[p+10] = x[p+9]*ga1 + x[p+8]*ga1*rho   - x[p+10]*om         - (x[p+10])*ej1 + PS[10]*ej*rp*u - x[o+10]*0 - x[o+10]*(cpo) + cLo*S1 + cMa*S1;
                        
                        dxdt[p+11] = x[p+10]*om                    - d2*x[p+11]*In*beta - (x[p+11])*ej1 + PS[11]*ej*rp*u - x[o+11]*0 - x[o+11]*(cpo) - cLo*S2 - cMa*S2;
                        dxdt[p+12] = d2*x[p+11]*In*beta            - x[p+12]*si         - (x[p+12])*ej1 + PS[12]*ej*rp*u - x[o+12]*0 - x[o+12]*(cpo);
                        dxdt[p+13] = x[p+12]*si*pA[a]              - x[p+13]*ga2*rho    - (x[p+13])*ej1 + PS[13]*ej*rp*u - x[o+13]*0 - x[o+13]*(cpo);
                        dxdt[p+14] = x[p+12]*si*(1.0-pA[a])        - x[p+14]*ga2        - (x[p+14])*ej1 + PS[14]*ej*rp*u - x[o+14]*0 - x[o+14]*(cpo);
                        dxdt[p+15] = x[p+14]*ga2 + x[p+13]*ga2*rho - x[p+15]*om         - (x[p+15])*ej1 + PS[15]*ej*rp*u - x[o+15]*0 - x[o+15]*(cpo) + cLo*S2 + cMa*S2;
                        
                        dxdt[p+16] = x[p+15]*om + x[p+20]*om       - d3*x[p+16]*In*beta - (x[p+16])*ej1 + PS[16]*ej*rp*u - x[o+16]*0 - x[o+16]*(cpo) - cLo*S3 - cMa*S3;
                        dxdt[p+17] = d3*x[p+16]*In*beta            - x[p+17]*si         - (x[p+17])*ej1 + PS[17]*ej*rp*u - x[o+17]*0 - x[o+17]*(cpo);
                        dxdt[p+18] = x[p+17]*si*pA[a]              - x[p+18]*ga3*rho    - (x[p+18])*ej1 + PS[18]*ej*rp*u - x[o+18]*0 - x[o+18]*(cpo);
                        dxdt[p+19] = x[p+17]*si*(1.0-pA[a])        - x[p+19]*ga3        - (x[p+19])*ej1 + PS[19]*ej*rp*u - x[o+19]*0 - x[o+19]*(cpo);
                        dxdt[p+20] = x[p+19]*ga3 + x[p+18]*ga3*rho - x[p+20]*om         - (x[p+20])*ej1 + PS[20]*ej*rp*u - x[o+20]*0 - x[o+20]*(cpo) + cLo*S3 + cMa*S3;
                        
                        // Vaccination groups
                        dxdt[p+21] = x_tot*cpo - lossP - x[p+21]*ej1 + PS[21]*ej*u*rp;  //pal
                        dxdt[p+22] = x_tot_1*cMo - lossMS0 - x[p+22]*ej1 + PS[22]*ej*u*rp; // + x[(a-1)*ag + s*sg + r*rg +22]*ej;
                        dxdt[p+23] = x_tot_2*cMo - lossMS1 - x[p+23]*ej1 + PS[23]*ej*u*rp; // + x[(a-1)*ag +

                       // dxdt[p+21] = mu*rp*cpmu + cpo*PST     - x[p+21]*ej1 + x[p+21]*ej1;
                       // dxdt[p+22] = mu*rp*cMmu + PST2*ej*u*rp*cMo - loss    - x[p+22]*ej1 + PS[22]*ej*u*rp*(1-cpo)*(1-cMo) ; // + x[(a-1)*ag + s*sg + r*rg +22]*ej;
                    }
                }
                // Vaccine groups
                dxdt[a*ag + 6*sg + 0] = si*(x[cj+3*sg+0*rg+2] + x[cj+4*sg+0*rg+2] + x[cj+5*sg+0*rg+2] + x[cj+3*sg+1*rg+2] + x[cj+4*sg+1*rg+2] + x[cj+5*sg+1*rg+2] + x[cj+3*sg+2*rg+2] + x[cj+4*sg+2*rg+2] + x[cj+5*sg+2*rg+2]);
                dxdt[a*ag + 6*sg + 1] = protectpal;
                dxdt[a*ag + 6*sg + 2] = protectmabs + protectLAV + protectmat;
               // cout << protectmabs << endl;
                /*dxdt[a*ag + 6*sg + 2] = protectmabs + vac_cal(t1, a)*(x[a*ag+4*sg+2*rg+1]+x[a*ag+4*sg+2*rg+6]+x[a*ag+4*sg+2*rg+11]+x[a*ag+4*sg+2*rg+16])
                + (cal_LAV_HR(t1, a))*((x[a*ag+0*sg+1*rg+1]+x[a*ag+0*sg+1*rg+6]+x[a*ag+0*sg+1*rg+11]+x[a*ag+0*sg+1*rg+16]) + (x[a*ag+1*sg+1*rg+1]+x[a*ag+1*sg+1*rg+6]+x[a*ag+1*sg+1*rg+11]+x[a*ag+1*sg+1*rg+16]))
                + (cal_LAV_LR(t1, a))*((x[a*ag+0*sg+2*rg+1]+x[a*ag+0*sg+2*rg+6]+x[a*ag+0*sg+2*rg+11]+x[a*ag+0*sg+2*rg+16]) + (x[a*ag+1*sg+2*rg+1]+x[a*ag+1*sg+2*rg+6]+x[a*ag+1*sg+2*rg+11]+x[a*ag+1*sg+2*rg+16]));
                */
                dxdt[a*ag + 6*sg + 3] = 0;
                
                dxdt[a*ag + 6*sg + 4] = 0;
                dxdt[a*ag + 6*sg + 5] = 0;
                dxdt[a*ag + 6*sg + 6] = 0;
                dxdt[a*ag + 6*sg + 7] = 0;
                
                // Monitoring Parent
                dxdt[a*ag + 6*sg + 8] =  si*(x[cj+3*sg+0*rg+2]+x[cj+3*sg+0*rg+7]+x[cj+3*sg+0*rg+12]+x[cj+3*sg+0*rg+17]); //VHR
                dxdt[a*ag + 6*sg + 9] =  si*(x[cj+3*sg+1*rg+2]+x[cj+3*sg+1*rg+7]+x[cj+3*sg+1*rg+12]+x[cj+3*sg+1*rg+17]); //HR
                dxdt[a*ag + 6*sg + 10] = si*(x[cj+3*sg+2*rg+2]+x[cj+3*sg+2*rg+7]+x[cj+3*sg+2*rg+12]+x[cj+3*sg+2*rg+17]); //LR
                
                // Monitoring Cocoon
                dxdt[a*ag + 6*sg + 11] = si*(x[cj+4*sg+0*rg+2]+x[cj+4*sg+0*rg+7]+x[cj+4*sg+0*rg+12]+x[cj+4*sg+0*rg+17]); //VHR
                dxdt[a*ag + 6*sg + 12] = si*(x[cj+4*sg+1*rg+2]+x[cj+4*sg+1*rg+7]+x[cj+4*sg+1*rg+12]+x[cj+4*sg+1*rg+17]); //HR
                dxdt[a*ag + 6*sg + 13] = si*(x[cj+4*sg+2*rg+2]+x[cj+4*sg+2*rg+7]+x[cj+4*sg+2*rg+12]+x[cj+4*sg+2*rg+17]); //LR
                
                // Monitoring Neither
                dxdt[a*ag + 6*sg + 14] = si*(x[cj+5*sg+0*rg+2]+x[cj+5*sg+0*rg+7]+x[cj+5*sg+0*rg+12]+x[cj+5*sg+0*rg+17]); //VHR
                dxdt[a*ag + 6*sg + 15] = si*(x[cj+5*sg+1*rg+2]+x[cj+5*sg+1*rg+7]+x[cj+5*sg+1*rg+12]+x[cj+5*sg+1*rg+17]); //HR
                dxdt[a*ag + 6*sg + 16] = si*(x[cj+5*sg+2*rg+2]+x[cj+5*sg+2*rg+7]+x[cj+5*sg+2*rg+12]+x[cj+5*sg+2*rg+17]); //LR
                
                x[a*ag + 6*sg + 17] = N_tot_n[a];
                x[a*ag + 6*sg + 18] = N_tot_c[a];
                x[a*ag + 6*sg + 19] = N_tot_p[a];
                
                x[a*ag + 6*sg + 20] = N_tot_n_v[a];
                x[a*ag + 6*sg + 21] = N_tot_c_v[a];
                x[a*ag + 6*sg + 22] = N_tot_p_v[a];
            }
        }
    };
    
    vector2D get_foi_mat1(param::param_state_t& pars, num_vec x, double t)
    {
        int ag = 455;
        int sg = 72;
        int rg = 24;
        long double alpha_i = pars.mapofParam.find("alpha_i")->second;
        long double qp = (pars.mapofParam.find("qp")->second);
        long double qc = (pars.mapofParam.find("qc")->second);
        long double a1 = 1; double a2 = 1; double a3 = 1;
        long double N_tot[NoAgeG];
        long double N_tot_n[NoAgeG];         double N_tot_c[NoAgeG];         double N_tot_p[NoAgeG];
        long double N_tot_n_inv, N_tot_p_inv, N_tot_c_inv;
        
        long double I_temp_nwn, I_temp_nwp, I_temp_nwc, I_temp_pwn, I_temp_pwp, I_temp_pwc, I_temp_cwc, I_temp_cwp, I_temp_cwn;
        long double pop_n, pop_m;
        num_vec I_pwp; num_vec I_pwc; num_vec I_pwn;
        num_vec I_cwp; num_vec I_cwc; num_vec I_cwn;
        num_vec I_nwp; num_vec I_nwc; num_vec I_nwn;
       
        vector2D I;
        long double tot_coc = 0;
        long double tot_inf = 0;
        
        for(int a = 0; a < NoAgeG; a++)
        {
            long double tot_temp_n = 0, tot_temp_c = 0, tot_temp_p = 0, tot_temp_nv = 0, tot_temp_cv = 0, tot_temp_pv = 0;
            long double tot_temp_vhr = 0, tot_temp_vhrv = 0, tot_temp_hr = 0, tot_temp_hrv = 0, tot_temp_lr = 0, tot_temp_lrv = 0;
            
            for (int j = 0; j < 24; j++){
                for (int r = 0; r < 3; r++){
                    tot_temp_p +=  x[a*ag+0*sg+r*rg+j];        tot_temp_c += x[a*ag+1*sg+r*rg+j];            tot_temp_n += x[a*ag+2*sg+r*rg+j];
                    tot_temp_pv += x[a*ag+3*sg+r*rg+j];       tot_temp_cv += x[a*ag+4*sg+r*rg+j];           tot_temp_nv += x[a*ag+5*sg+r*rg+j];
                }
            }
            
            for (int j = 0; j < 24; j++){
                for (int s = 0; s < 3; s++){
                    tot_temp_vhr +=  x[a*ag+s*sg+0*rg+j];        tot_temp_hr += x[a*ag+s*sg+1*rg+j];            tot_temp_lr += x[a*ag+s*sg+2*rg+j];
                    tot_temp_vhrv += x[a*ag+(s+3)*sg+0*rg+j];       tot_temp_hrv += x[a*ag+(s+3)*sg+1*rg+j];           tot_temp_lrv += x[a*ag+(s+3)*sg+2*rg+j];
                }
            }
            // Find proportions
            /*prop_n[a] = tot_temp_n/(popsize[a]);    prop_c[a] = tot_temp_c/(popsize[a]);   prop_p[a] = tot_temp_p/(popsize[a]);
             prop_nv[a] = tot_temp_nv/(popsize[a]); prop_cv[a] = tot_temp_cv/(popsize[a]); prop_pv[a] = tot_temp_pv/(popsize[a]);*/
            
            if (a < 12){
                tot_coc += tot_temp_c;
                tot_inf += tot_temp_c + tot_temp_n;
            }
            N_tot_n[a] = tot_temp_n;    N_tot_c[a] = tot_temp_c;    N_tot_p[a] = tot_temp_p;
        }
        
        /*/////////////*/
        /* 2. Force of Infection /*/
        /*/////////////*/
        long double phi_c = tot_coc/tot_inf;
        
        //  2.1 Get contact matrices
        vector2D cnt_matrix_cwc_p = manip_contact::get_cwc(phi_c,'p');
        vector2D cnt_matrix_cwp_p = manip_contact::get_cwp(phi_c,'p');
        vector2D cnt_matrix_cwn_p = manip_contact::get_cwn(phi_c,'p');
        vector2D cnt_matrix_pwc_p = manip_contact::get_pwc(phi_c,'p');
        vector2D cnt_matrix_pwp_p = manip_contact::get_pwp(phi_c,'p');
        vector2D cnt_matrix_pwn_p = manip_contact::get_pwn(phi_c,'p');
        vector2D cnt_matrix_nwc_p = manip_contact::get_nwc(phi_c,'p');
        vector2D cnt_matrix_nwp_p = manip_contact::get_nwp(phi_c,'p');
        vector2D cnt_matrix_nwn_p = manip_contact::get_nwn(phi_c,'p');
        
        vector2D cnt_matrix_cwc_c = manip_contact::get_cwc(phi_c,'c');
        vector2D cnt_matrix_cwp_c = manip_contact::get_cwp(phi_c,'c');
        vector2D cnt_matrix_cwn_c = manip_contact::get_cwn(phi_c,'c');
        vector2D cnt_matrix_pwc_c = manip_contact::get_pwc(phi_c,'c');
        vector2D cnt_matrix_pwp_c = manip_contact::get_pwp(phi_c,'c');
        vector2D cnt_matrix_pwn_c = manip_contact::get_pwn(phi_c,'c');
        vector2D cnt_matrix_nwc_c = manip_contact::get_nwc(phi_c,'c');
        vector2D cnt_matrix_nwp_c = manip_contact::get_nwp(phi_c,'c');
        vector2D cnt_matrix_nwn_c = manip_contact::get_nwn(phi_c,'c');
        
        //  2.2 Seasonal forcing
        
        for(int a = 0; a < NoAgeG; a++)
        {
            for (int k = 0; k < NoAgeG; k++)
            {
                I_temp_pwp = 0.0;    I_temp_pwc = 0.0; I_temp_pwn = 0.0;
                I_temp_cwp = 0.0;    I_temp_cwc = 0.0; I_temp_cwn = 0.0;
                I_temp_nwp = 0.0;    I_temp_nwc = 0.0; I_temp_nwn = 0.0;

                if (N_tot_n[k] < 0.1){N_tot_n_inv = 0;}else{N_tot_n_inv=1.0/N_tot_n[k];}
                if (N_tot_c[k] < 0.1){N_tot_c_inv = 0;}else{N_tot_c_inv=1.0/N_tot_c[k];}
                if (N_tot_p[k] < 0.1){N_tot_p_inv = 0;}else{N_tot_p_inv=1.0/N_tot_p[k];}
                
                for (int r = 0; r < 3; r++)
                {
                    I_temp_pwp += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_pwp_p[a][k]+qc*cnt_matrix_pwp_c[a][k]))*N_tot_p_inv;
                    
                    I_temp_pwc += (x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_pwc_p[a][k]+qc*cnt_matrix_pwc_c[a][k]))*N_tot_c_inv;
                    
                    I_temp_pwn += (x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_pwn_p[a][k]+qc*cnt_matrix_pwn_c[a][k]))*N_tot_n_inv;
                    
                    I_temp_cwp += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_cwp_p[a][k]+qc*cnt_matrix_cwp_c[a][k]))*N_tot_p_inv;
                    
                    I_temp_cwc += (x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_cwc_p[a][k]+qc*cnt_matrix_cwc_c[a][k]))*N_tot_c_inv;
                    
                    I_temp_cwn += (x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_cwn_p[a][k]+qc*cnt_matrix_cwn_c[a][k]))*N_tot_n_inv;
                    
                    I_temp_nwp += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_nwp_p[a][k]+qc*cnt_matrix_nwp_c[a][k]))*N_tot_p_inv;
                    
                    I_temp_nwc += (x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_nwc_p[a][k]+qc*cnt_matrix_nwc_c[a][k]))*N_tot_c_inv;
                    
                    I_temp_nwn += (x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_nwn_p[a][k]+qc*cnt_matrix_nwn_c[a][k]))*N_tot_n_inv;
                }
                I_pwp.push_back(I_temp_pwp); I_pwc.push_back(I_temp_pwc); I_pwn.push_back(I_temp_pwn);
                I_cwp.push_back(I_temp_cwp); I_cwc.push_back(I_temp_cwc); I_cwn.push_back(I_temp_cwn);
                I_nwp.push_back(I_temp_nwp); I_nwc.push_back(I_temp_nwc); I_nwn.push_back(I_temp_nwn);
            }
        }
        I.push_back(I_pwp);I.push_back(I_pwc);I.push_back(I_pwn);
        I.push_back(I_cwp);I.push_back(I_cwc);I.push_back(I_cwn);
        I.push_back(I_nwp);I.push_back(I_nwc);I.push_back(I_nwn);

        return I;
    }
}


#endif /* interventions_h */
