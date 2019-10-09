//
//  writing.h
//  flu_forecast
//
//  Created by David Hodgson on 17/06/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

#ifndef writing_h
#define writing_h

namespace sim_ouput
{
    /**********************************/
    /**    WRTIE RESULTS OF CALIBRATION  **/
    /**********************************/
    num_vec inc_annual_post(param::param_state_t& pars, amh::amh_state_t& mcmc_state, char F)
    {
        VectorXd sample_post;
        EulerT<state_t> get_newInf;
        SystemT<state_t, system_t> System;
        asc::Recorder recorder;
        get_ll::ODE_dynamics ode_par(0, 365*2, 365*3, 1.0);
        int i = uniform_dist_disc(0, mcmc_state.NK2-1, 'r');
        int t_w = 0;
        num_vec inci(NoAgeG,0);
        VectorXd inc_tot(NoAgeG);
        inc_tot = VectorXd::Zero(NoAgeG);
        
        sample_post = mcmc_state.posterior.row(i).transpose();
        update(pars, sample_post);
        
        double t = ode_par.t_start;
        int t_d = 0;
        
        get_ll::ODE_desc ODE_desc_inst(pars);
        num_vec x0 = get_ll::init_cond(pars);
        
        while (t < ode_par.t_end)
        {
            get_newInf(ODE_desc_inst, x0, t, 1);
            if (t > ode_par.t_burn)
            {
                if (t_d == 0)
                {
                    for (int a = 0; a < NoAgeG; a++)
                        x0[22 + 23*a] = 0.0; //Incidence at t_d = 0;
                }
                t_d++;
            }
        }
        for (int a = 0; a < NoAgeG; a++)
        {
            inci[a] = (x0[22 + 23*a]);

            inc_tot(a) = inci[a];
        }
        
        return inci;
    }
    
    vector2D inc_soln_post(param::param_state_t& pars, amh::amh_state_t& mcmc_state, char F)
    {
        VectorXd sample_post;
        EulerT<state_t> get_newInf;
        SystemT<state_t, system_t> System;
        asc::Recorder recorder;
        get_ll::ODE_dynamics ode_par(0, 365, 365*8, 1.0);
        int i = uniform_dist_disc(0, mcmc_state.NK2-1, 'r');
        int t_w = 0;
        vector2D inci(52*7+1,vector<double>(NoAgeG,0));
        
        sample_post = mcmc_state.posterior.row(i).transpose();
        update(pars, sample_post);
        
        double t = ode_par.t_start;
        int t_d = 0;
        
        get_ll::ODE_desc ODE_desc_inst(pars);
        num_vec x0 = get_ll::init_cond(pars);
        
        while (t < ode_par.t_end)
        {
            get_newInf(ODE_desc_inst, x0, t, 1);
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
                        inci[t_w][a] = x0[22 + 23*a]; //Incidence at t_d = 7;
                        x0[22 + 23*a] = 0.0;
                    }
                    t_w ++;
                }
                t_d++;
            }
        }
        return inci;
    }
    
    vector2D pos_soln_post(param::param_state_t& pars, amh::amh_state_t& mcmc_state, char F)
    {
        VectorXd sample_post;
        EulerT<state_t> get_newInf;
        SystemT<state_t, system_t> System;
        asc::Recorder recorder;
        
        get_ll::ODE_dynamics ode_par(0, 365, 365*8, 1.0);
        VectorXd ep_t(NoAgeG);
        int i = uniform_dist_disc(0, mcmc_state.NK2-1, 'r');
        int t_w = 0;
        vector2D inci(52*7+1,vector<double>(NoAgeG,0));
        
        sample_post = mcmc_state.posterior.row(i).transpose();
        update(pars, sample_post);
        
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
            ep_t(a) = pars.mapofParam.find(pars.c5mA[a])->second ;
#endif
        
        double t = ode_par.t_start;
        int t_d = 0;
        
        get_ll::ODE_desc ODE_desc_inst(pars);
        num_vec x0 = get_ll::init_cond(pars);
        
        while (t < ode_par.t_end)
        {
            get_newInf(ODE_desc_inst, x0, t, 1);
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
                        inci[t_w][a] = x0[22 + 23*a]*ep_t(a); //Incidence at t_d = 7;
                        x0[22 + 23*a] = 0.0;
                    }
                    t_w ++;
                }
                t_d++;
            }
        }
        return inci;
        //recorder.csv(dout + calout_name + "soln/" + "possoln_sample_" + F, {"x"});
    }
    
    void write_inc(param::param_state_t& pars, amh::amh_state_t& mcmc_state, char F, num_vec seed)
    {
        //ode_soln_post_par(pars, mcmc_state, F, seed[0]);
        
        asc::Recorder recorder_inci;
        asc::Recorder recorder_pos;
        asc::Recorder recorder_annual;
        
        int no_samples = 1000;
        int no_weeks = 7*52+1;
        vector3D x_inci_t(no_samples,vector<vector<double>>(no_weeks,vector<double>(NoAgeG,0)));
        vector3D x_inci(NoAgeG,vector<vector<double>>(no_weeks,vector<double>(no_samples,0)));
        vector3D x_pos_t(no_samples,vector<vector<double>>(no_weeks,vector<double>(NoAgeG,0)));
        vector3D x_pos(NoAgeG,vector<vector<double>>(no_weeks,vector<double>(no_samples,0)));
        vector2D x_annual_t(no_samples,vector<double>(NoAgeG,0));
        vector2D x_annual(NoAgeG,vector<double>(no_samples,0));

        num_vec filtered_inci;
        num_vec filtered_pos;
        num_vec filtered_annual;

        
        for (int s = 0; s < no_samples; s++)
        {
            x_annual_t[s] = inc_annual_post(pars, mcmc_state, F);
            x_inci_t[s] = inc_soln_post(pars, mcmc_state, F);
            x_pos_t[s]  = pos_soln_post(pars, mcmc_state, F);
        }
        
        for (int a = 0; a < NoAgeG; a++)
            for (int t = 0; t < 7*52; t++)
                for (int s = 0; s < no_samples; s++)
                {
                    x_inci[a][t][s] = x_inci_t[s][t][a];
                    x_pos[a][t][s] = x_pos_t[s][t][a];
                }
        
        for (int a = 0; a < NoAgeG; a++)
            for (int s = 0; s < no_samples; s++)
            {
                x_annual[a][s] = x_annual_t[s][a];
            }
        
        string outmcyr = dout + calout_name + "soln/" + "annual_sample_" + F + ".txt";
        //string FileMC_chain_M_full = outmcyr + "MC_chain_M_full_" + F + ".txt";      //Output MC
        ofstream annual_out(outmcyr.c_str());
        for (int a = 0; a < NoAgeG; a++)
        {
            
            std::sort (x_annual[a].begin(), x_annual[a].end());
            annual_out << x_annual[a][5]*(1-pHR[a]) << endl;
            annual_out << x_annual[a][5]*pHR[a] << endl;
        }
        annual_out.close();
        
        for (int a = 0; a < NoAgeG; a++)
        {
            for (int t = 0; t < 7*52; t++)
            {
                filtered_inci.clear();
                std::sort (x_inci[a][t].begin(), x_inci[a][t].end());
                filtered_inci.push_back(x_inci[a][t][25]);
                filtered_inci.push_back(x_inci[a][t][250]);
                filtered_inci.push_back(x_inci[a][t][500]);
                filtered_inci.push_back(x_inci[a][t][750]);
                filtered_inci.push_back(x_inci[a][t][975]);
                
                recorder_inci({(double)a});
                recorder_inci.add(t);
                recorder_inci.add(filtered_inci);
            }
        }
        
        recorder_inci.csv(get_ll::dout + get_ll::calout_name + "soln/" + "inci_sample_" + F, {"Age group", "t", "0.025", "0.25","0.5", "0.75", "0.975"});
        
        for (int a = 0; a < NoAgeG; a++)
        {
            for (int t = 0; t < 7*52; t++)
            {
                filtered_pos.clear();
                std::sort (x_pos[a][t].begin(), x_pos[a][t].end());
                filtered_pos.push_back(x_pos[a][t][25]);
                filtered_pos.push_back(x_pos[a][t][250]);
                filtered_pos.push_back(x_pos[a][t][500]);
                filtered_pos.push_back(x_pos[a][t][750]);
                filtered_pos.push_back(x_pos[a][t][975]);
                
                recorder_pos({(double)t});
                recorder_pos.add(filtered_pos);
            }
        }
        recorder_pos.csv(get_ll::dout + get_ll::calout_name + "soln/" + "pos_sample_" + F, {"Age group", "t", "0.025", "0.25","0.5", "0.75", "0.975"});
    
    }
    
   /* void write_Rs(param::param_state_t& pars, amh::amh_state_t& mcmc_state, int no_samples)
    {
        asc::Recorder recorder;
        VectorXd sample_post;
        VectorXd R0_v;
        VectorXd REff_v;
        R0_v = VectorXd::Zero(no_samples);
        REff_v = VectorXd::Zero(no_samples);
        
        for (int i = 0; i < no_samples; i++)
        {
            int s = uniform_dist_disc(0, mcmc_state.NK2-1, 'r');
            sample_post = mcmc_state.posterior.row(s).transpose();
            update(pars, sample_post);
            
            get_ll::ODE_desc ODE_i(pars);
            R0_v(i) = ODE_i.find_R0(ODE_i);
            REff_v(i) = ODE_i.find_REff(ODE_i);
            
            recorder({(double)s});
            recorder.add(R0_v(i));
            recorder.add(REff_v(i));
        }
        recorder.csv(get_ll::dout + get_ll::calout_name + "soln/" + "R0_sample" , {"s","R0","REff"});
    }
    */
    
    void write_foi(param::param_state_t& pars, amh::amh_state_t& mcmc_state, int i)
    {
        VectorXd sample_post;
        asc::Recorder recorder;
        asc::Euler integrator;

        get_ll::ODE_dynamics ode_par(0.0, 365.0*2, 365*3, 1.0);
        
        sample_post = mcmc_state.posterior.row(i).transpose();
        update(pars, sample_post);
        ODE_desc_par ODE_desc_inst(pars);
        num_vec x0 = init_cond_par(pars);
        
        num_vec inci_tot;
        vector2D foi(4,num_vec(25*25,0));
        num_vec foi_sum_pwp(25*25,0);
        num_vec foi_sum_nwp(25*25,0);
        num_vec foi_sum_pwn(25*25,0);
        num_vec foi_sum_nwn(25*25,0);
        
        double t = ode_par.t_start;
        while (t < ode_par.t_end)
        {
            integrator(ODE_desc_inst, x0, t, 1.0);
            
            if (t > ode_par.t_burn)
            {
                foi = get_foi(pars, x0, t);
                for (int i = 0; i < 3; i++)
                    for (int a = 0; a < 25*25; a++)
                    {
                        foi_sum_nwn[a] += foi[0][a];
                        foi_sum_nwp[a] += foi[1][a];
                        foi_sum_pwn[a] += foi[2][a];
                        foi_sum_pwp[a] += foi[3][a];
                    }
            }
        }
        recorder({(double)0});
        recorder.add(foi_sum_nwn);
        recorder({(double)1});
        recorder.add(foi_sum_nwp);
        recorder({(double)2});
        recorder.add(foi_sum_pwn);
        recorder({(double)3});
        recorder.add(foi_sum_pwp);
        recorder.csv(get_ll::dout + get_ll::calout_name + "soln/" + "foi_" + "All", {"Foi", "value"});
    }
    
    void write_foi_mat(param::param_state_t& pars, amh::amh_state_t& mcmc_state, int i, double om_mab, double xi_b)
    {
        VectorXd sample_post;
        asc::Recorder recorder1, recorder2;
        asc::Euler integrator;
        
        VectorXd eff_pal = cal::get_eff(10, 1);
        VectorXd eff_mab =  cal::get_eff(10, 2);
        VectorXd eff_vac =  cal::get_eff(10, 3);
        VectorXd eff_mat =  cal::get_eff(10, 4);
        
        get_ll::ODE_dynamics ode_par(0.0, 0.0, 365*3, 1.0);
        
        sample_post = mcmc_state.posterior.row(i).transpose();
        update(pars, sample_post);
        
        cal::Calendar_full cal(cal::G_base, 0.0, {0}, 0, 0, true, "None", i, eff_pal, eff_mab, eff_vac, eff_mat);
        interventions::ODE_Desc_inter ODE_desc_int_inst(pars, cal, false, 0.6, om_mab, xi_b);
        num_vec x01 = interventions::initial_cond(pars, 0.6);
        
        cal::Calendar_full cal_mat(cal::G_par, 0.6, {0}, 48, 48+21, true, "mat", i, eff_pal, eff_mab, eff_vac, eff_mat);
        interventions::ODE_Desc_inter ODE_desc_int_inst_mat(pars, cal_mat, false, 0.6, om_mab, xi_b);
        num_vec x02 = interventions::initial_cond(pars, 0.6);
        
        num_vec inci_tot;
        vector2D foi(9,num_vec(25*25,0));

        double t1 = ode_par.t_start;
        while (t1 < ode_par.t_end)
        {
            integrator(ODE_desc_int_inst, x01, t1, 1.0);
            
            if (t1 > ode_par.t_burn)
            {
                foi = interventions::get_foi_mat1(pars, x01, t1);
                recorder1({(double)t1});
                recorder1.add(foi[0]);
                recorder1.add(foi[1]);
                recorder1.add(foi[2]);
                recorder1.add(foi[3]);
                recorder1.add(foi[4]);
                recorder1.add(foi[5]);
                recorder1.add(foi[6]);
                recorder1.add(foi[7]);
                recorder1.add(foi[8]);
            }
        }
        
        double t2 = ode_par.t_start;
        while (t2 < ode_par.t_end)
        {
            integrator(ODE_desc_int_inst_mat, x02, t2, 1.0);
            
            if (t2 > ode_par.t_burn)
            {
                foi = interventions::get_foi_mat1(pars, x02, t2);
                recorder2({(double)t2});
                recorder2.add(foi[0]);
                recorder2.add(foi[1]);
                recorder2.add(foi[2]);
                recorder2.add(foi[3]);
                recorder2.add(foi[4]);
                recorder2.add(foi[5]);
                recorder2.add(foi[6]);
                recorder2.add(foi[7]);
                recorder2.add(foi[8]);
            }
        }
        recorder1.csv(get_ll::dout + get_ll::calout_name + "soln/" + "foi_mat_" + "1", {"Foi", "value"});
        recorder2.csv(get_ll::dout + get_ll::calout_name + "soln/" + "foi_mat_" + "2", {"Foi", "value"});

    }
    
    void write_pR(param::param_state_t& pars, amh::amh_state_t& mcmc_state, int i)
    {
        VectorXd sample_post;
        asc::Recorder recorder;
        asc::Euler integrator;
        
        get_ll::ODE_dynamics ode_par(0.0, 365.0*2, 365*3, 1.0);
        
        sample_post = mcmc_state.posterior.row(i).transpose();
        update(pars, sample_post);
        ODE_desc_par ODE_desc_inst(pars);
        num_vec x0 = init_cond_par(pars);
        
        num_vec inci_tot;
        num_vec pR(365,0);
        int t_d = 0;
        double t = ode_par.t_start;
        while (t < ode_par.t_end)
        {
            integrator(ODE_desc_inst, x0, t, 1.0);
            
            if (t > ode_par.t_burn)
            {
                pR[t_d] = get_pR(pars, x0, t_d);
                t_d++;
            }
        }
        recorder({(double)0});
        recorder.add(pR);
        recorder.csv(get_ll::dout + get_ll::calout_name + "soln/" + "pR", { "s", "week_no"});
    }
    
    
    void consistency_checks(param::param_state_t& pars, amh::amh_state_t& mcmc_state, int s, double om_mab, double xi_b)
    {
        VectorXd eff_pal = cal::get_eff(10, 1);
        VectorXd eff_mab =  cal::get_eff(10, 2);
        VectorXd eff_vac =  cal::get_eff(10, 3);
        VectorXd eff_mat =  cal::get_eff(10, 4);

        //cal::Calendar_full cal(cal::G_base, 0.0, {0}, 0, 0, false, "None", s, eff_pal, eff_mab, eff_vac, eff_mat);
        
        cal::Calendar_full cal_mab(cal::G_0mo, 0.0, {0}, 0, 0, true, "None", s, eff_pal, eff_mab, eff_vac, eff_mat);
        cal::Calendar_full cal_mab_2(cal::G_0mo, 0.9, {0}, 12, 12+21, false, "Mlr_p", s, eff_pal, eff_mab, eff_vac, eff_mat);

        EulerT<state_t> get_newInf;
        SystemT<state_t, system_t> System;
        asc::Recorder recorder_1; asc::Recorder recorder_2; asc::Recorder recorder_3; asc::Recorder recorder_4;

        get_ll::ODE_dynamics ode_par(0, 365*2, 365*3, 1.0);
        
        VectorXd sample_post;
        sample_post = mcmc_state.posterior.row(s).transpose();
        update(pars, sample_post);
        
        double t1 = ode_par.t_start; double t2 = ode_par.t_start; double t3 = ode_par.t_start;  double t4 = ode_par.t_start;
        
        get_ll::ODE_desc ODE_desc_inst(pars);
        get_ll::ODE_desc_par ODE_desc_par_inst(pars);
        interventions::ODE_Desc_inter ODE_desc_int_inst_m1(pars, cal_mab, false, 0.0, om_mab, xi_b);
        interventions::ODE_Desc_inter ODE_desc_int_inst_m2(pars, cal_mab_2, false, 0.0, om_mab, xi_b);
        //interventions::ODE_Desc_inter ODE_desc_int_inst(pars, cal, false, 0.6, om_mab, xi_b);
        //interventions::ODE_Desc_inter ODE_desc_int_inst_mat(pars, cal_mat, false, 0.6, om_mab, xi_b);

        num_vec x01 = get_ll::init_cond(pars);
        num_vec x02 = init_cond_par(pars);
        num_vec x03 = interventions::initial_cond(pars, 0.0);
        num_vec x04 = interventions::initial_cond(pars, 0.0);
        
        while (t1 < ode_par.t_end)
        {
            recorder_1({(double)t1});
            recorder_1.add(x01);
            get_newInf(ODE_desc_inst, x01, t1, 1);
        }
        while (t2 < ode_par.t_end)
        {
            recorder_2({(double)t2});
            recorder_2.add(x02);
            get_newInf(ODE_desc_par_inst, x02, t2, 1);
        }
        while (t3 < ode_par.t_end)
        {
            recorder_3({(double)t3});
            recorder_3.add(x03);
            get_newInf(ODE_desc_int_inst_m1, x03, t3, 1);
        }
        while (t4 < ode_par.t_end)
        {
            recorder_4({(double)t4});
            recorder_4.add(x04);
            get_newInf(ODE_desc_int_inst_m2, x04, t4, 1);
        }
        
        recorder_1.csv(get_ll::dout + get_ll::calout_name + "soln/" + "full_1", { "t", "x0"});
        recorder_2.csv(get_ll::dout + get_ll::calout_name + "soln/" + "full_2", { "t", "x0"});
        recorder_3.csv(get_ll::dout + get_ll::calout_name + "soln/" + "full_3", { "t", "x0"});
        recorder_4.csv(get_ll::dout + get_ll::calout_name + "soln/" + "full_3_mat", { "t", "x0"});
    }

    
    /**********************************/
    /**    WRTIE INTERVENTION PROGRAMMES  **/
     /**********************************/
    
    num_vec int_post(param::param_state_t& pars, amh::amh_state_t& mcmc_state, cal::Calendar_full cal, bool dir_ind, int s, double cov_c, num_vec& S_tot, num_vec& H_tot, num_vec& D_tot, num_vec& GP_tot, num_vec& BD_tot, vector2D& doses, double& Q, double& CT, double& CP, double p_cost, int time_hor, double disc, double om_mab, double xi_b)
    {
        VectorXd sample_post;
        asc::Euler integrator;
        get_ll::ODE_dynamics ode_par(0.0, 365.0*1, 365.0*4, 1.0);
        
        sample_post = mcmc_state.posterior.row(s).transpose();
        update(pars, sample_post);
        
        interventions::ODE_Desc_inter ODE_desc_inst(pars, cal, dir_ind, cov_c, om_mab, xi_b);
        num_vec x0 = interventions::initial_cond(pars, cov_c);
        
        num_vec inci_tot;
        num_vec inci_temp;
        double r = disc;
        double t = ode_par.t_start;
        int t_d = 0;
        vector2D x_inc(52,vector<double>(NoAgeG*9,0));
        vector2D x_pri(52,vector<double>(NoAgeG,0));
        double x_pal = 0;
        double x_vac = 0;

        num_vec x_pal_inc(52,0);
        num_vec x_vac_inc(52,0);
        
        double x_inc_annual, x_pri_annual;
        int t_w = 0;
   
    // Run model until equilibrium
        while (t < ode_par.t_end)
        {
            integrator(ODE_desc_inst, x0, t, 1.0);
            
            if (t > ode_par.t_burn)
            {
                if (t_d == 0 || t_d%7 == 0)
                {
                    if (t_d > 0)
                    {
                        inci_temp.clear();
                        // First get the incidence
                        for (int a = 0; a < NoAgeG; a++)
                            for (int j = 0; j < 9; j++)
                            {
                                x_inc[(t_w)%52][9*a + j] = x0[455*a + 72*6 + 8 + j];
                                inci_temp.push_back(x0[455*a + 72*6 + 8 + j]);
                            }
                        // Primary incidence get the incidence
                        for (int a = 0; a < NoAgeG; a++)
                            x_pri[(t_w)%52][a] = (x0[455*a + 72*6 + 0]);
                        // Get number of Palivizumab
                        for (int a = 0; a < NoAgeG; a++)
                            x_pal += (x0[455*a + 72*6 + 1]);
                        // Get number of Non-palivizumab
                        for (int a = 0; a < NoAgeG; a++)
                            x_vac += (x0[455*a + 72*6 + 2]);
   
                        x_pal_inc[t_w%52] = x_pal;
                        x_vac_inc[t_w%52] = x_vac;

                        Q += get_QALY_2(inci_temp, mcmc_state, S_tot, GP_tot, H_tot, D_tot, s)*exp(-t_w*r/52.0);
                        CT += get_CostT(inci_temp, GP_tot, BD_tot, s)*exp(-t_w*r/52.0);
                        CP += get_CostP(x_pal, x_vac, p_cost, doses, t_w, s)*exp(-t_w*r/52.0);
                        t_w++;
                    }
                    for (int a = 0; a < NoAgeG; a++)
                    {
                        x_vac = 0; x_pal = 0;
                        x0[455*a + 72*6 + 0] = 0;
                        x0[455*a + 72*6 + 1] = 0;
                        x0[455*a + 72*6 + 2] = 0;
                        for (int j = 0; j < 9; j++)
                        {
                            x0[455*a + 72*6 + 8 + j] = 0;
                        }
                    }
                }
                t_d++;
            }
        }
   // After equilibirum just use last year with an increasing t_w;
        for (int t_w = 4*52; t_w  < (time_hor+1)*52; t_w++)
        {
            inci_temp.clear();
            for (int a = 0; a < NoAgeG; a++)
                for (int j = 0; j < 9; j++)
                    inci_temp.push_back(x_inc[(t_w)%52][9*a + j]);
            
            if (t_w == (time_hor)*52)
            {
                for (int a = 0; a < NoAgeG; a++)
                {
                    S_tot[a] = 0; H_tot[a] = 0; D_tot[a] = 0;
                    GP_tot[a] = 0; BD_tot[a] = 0;
                }
            }
            Q += get_QALY_2(inci_temp, mcmc_state, S_tot, GP_tot, H_tot, D_tot, s)*exp(-t_w*r/52.0);
            CT += get_CostT(inci_temp, GP_tot, BD_tot, s)*exp(-t_w*r/52.0);
            CP += get_CostP(x_pal_inc[t_w%52], x_vac_inc[t_w%52], p_cost, doses, t_w, s)*exp(-t_w*r/52.0);
        }
  // Find the total annual incidences for the equilibirum year
        for (int a = 0; a < NoAgeG; a++)
        {
            for (int j = 0; j < 9; j++)
            {
                x_inc_annual = 0;
                for (int t_t = 0; t_t < 52 ; t_t++)
                {
                    x_inc_annual += x_inc[t_t][9*a + j];
                }
                inci_tot.push_back(x_inc_annual);
            }
        }
        for (int a = 0; a < NoAgeG; a++)
        {
            x_pri_annual = 0;
            for (int t_t = 0; t_t < 52; t_t++)
                x_pri_annual += x_pri[t_t][a];

            inci_tot.push_back(x_pri_annual);
        }
        return inci_tot;
    }
        
    void write_interventions(param::param_state_t& pars, amh::amh_state_t& mcmc_state, string prog_no, string prog_name, num_vec seed, double cov_c, double p_cost, int time_hor, VectorXd eff_pal, VectorXd mono_pal, VectorXd vac_pal, VectorXd mat_pal, int iN, double disc, double om_mab, double xi_b)
    {
        asc::Recorder record_inc, record_inc_pri, record_s, record_h, record_d, record_gp, record_bd, record_dose, record_q, record_cp, record_ct;
        asc::Recorder record_inc_d, record_inc_d_pri, record_s_d, record_h_d, record_d_d, record_gp_d, record_bd_d, record_dose_d, record_q_d, record_cp_d, record_ct_d;
        
        str_vec col_name = {"seed", "<1 mo", "1 mo", "2 mo", "3 mo", "4 mo", "5 mo", "6 mo", "7 mo", "8 mo", "9 mo", "10 mo", "11 mo", "1 yr", "2 yr", "3 yr", "4 yr", "5-9yrs", "10-14 yrs", "15-24 yrs", "25-34 yrs", "35-44 yrs", "44-54 yrs", "55-64yrs", "65-74yrs", "75yrs+"};
        num_vec inciall;
        num_vec inci(25*9,0);
        num_vec incipri(25,0);
        num_vec S_tot(NoAgeG, 0); num_vec H_tot(NoAgeG, 0); num_vec D_tot(NoAgeG, 0);
        num_vec GP_tot(NoAgeG, 0); num_vec BD_tot(NoAgeG, 0);
        vector2D doses(2,vector<double>(time_hor*52,0));
        double Q, CP, CT;
        for (int s = 0; s < seed.size(); s++)
        {
            num_vec up_take_base = cal::gen_daily(cal::uprate[iN], cal::start_w[iN]);
            cal::Calendar_full cal(cal::t_group[iN], cal::cov[iN], up_take_base, cal::start_w[iN], cal::end_w[iN], cal::Pal_ind[iN], cal::cal_type[iN], s, eff_pal, mono_pal, vac_pal, mat_pal);
            
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
            
            inciall = int_post(pars, mcmc_state, cal, false, seed[s], cov_c, S_tot, H_tot, D_tot, GP_tot, BD_tot, doses, Q, CP, CT, p_cost, time_hor, disc, om_mab, xi_b);
            for (int i = 0; i < 9*25; i++)
                inci[i] = inciall[i];
            
            for (int i = 0; i < 25; i++)
                incipri[i] = inciall[i + 25*9];
            
            record_inc({(double)seed[s]}); record_inc.add(inci);
            record_inc_pri({(double)seed[s]}); record_inc_pri.add(incipri);
            record_s({(double)seed[s]}); record_s.add(S_tot);
            record_h({(double)seed[s]}); record_h.add(H_tot);
            record_d({(double)seed[s]}); record_d.add(D_tot);
            record_gp({(double)seed[s]}); record_gp.add(GP_tot);
            record_bd({(double)seed[s]}); record_bd.add(BD_tot);
            record_dose({(double)seed[s]}); record_dose.add(doses[0]);
            record_dose({(double)seed[s]}); record_dose.add(doses[1]);
            record_q({(double)seed[s]}); record_q.add(Q);
            record_cp({(double)seed[s]}); record_cp.add(CP);
            record_ct({(double)seed[s]}); record_ct.add(CT);
        }
        
        record_inc.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "inc", col_name);
        record_inc_pri.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "incpri", col_name);
        record_s.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "symp", col_name);
        record_h.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "hosp", col_name);
        record_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "death", col_name);
        record_gp.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "gp", col_name);
        record_bd.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "bd", col_name);
        record_dose.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "no_doses",  {"seed", "Total numver of doses"});
        record_q.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "qaly", {"seed", "Total QALY"});
        record_cp.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "cost_pro",  {"seed", "Total Cost T"});
        record_ct.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "cost_tre",  {"seed", "Total Cost T"});

        for (int s = 0; s < seed.size(); s++)
        {
            num_vec up_take_base = cal::gen_daily(cal::uprate[iN], cal::start_w[iN]);
            cal::Calendar_full cal(cal::t_group[iN], cal::cov[iN], up_take_base, cal::start_w[iN], cal::end_w[iN], cal::Pal_ind[iN], cal::cal_type[iN], s, eff_pal, mono_pal, vac_pal, mat_pal );
            

            Q = CP = CT = 0;
            for (int a = 0; a < NoAgeG; a++)
            {
                S_tot[a] = 0; H_tot[a] = 0; D_tot[a] = 0;
                GP_tot[a] = 0; BD_tot[a] = 0;
            }
            
            inciall = int_post(pars, mcmc_state, cal, true, seed[s], cov_c, S_tot, H_tot, D_tot, GP_tot, BD_tot, doses, Q, CP, CT, p_cost, time_hor, disc, om_mab, xi_b);
            for (int i = 0; i < 9*25; i++)
                inci[i] = inciall[i];
            
            for (int i = 0; i < 25; i++)
                incipri[i] = inciall[i + 25*9];
            
            record_inc_d({(double)seed[s]}); record_inc_d.add(inci);
            record_inc_d_pri({(double)seed[s]}); record_inc_d_pri.add(incipri);
            
            record_s_d({(double)seed[s]}); record_s_d.add(S_tot);
            record_h_d({(double)seed[s]}); record_h_d.add(H_tot);
            record_d_d({(double)seed[s]}); record_d_d.add(D_tot);
            record_gp_d({(double)seed[s]}); record_gp_d.add(GP_tot);
            record_bd_d({(double)seed[s]}); record_bd_d.add(BD_tot);
            record_dose_d({(double)seed[s]}); record_dose_d.add(doses[0]);
            record_dose_d({(double)seed[s]}); record_dose_d.add(doses[1]);

            record_q_d({(double)seed[s]}); record_q_d.add(Q);
            record_cp_d({(double)seed[s]}); record_cp_d.add(CP);
            record_ct_d({(double)seed[s]}); record_ct_d.add(CT);
        }
        
        record_inc_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "inc_d", col_name);
        record_inc_d_pri.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "inc_pri_d", col_name);
        record_s_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "symp_d", col_name);
        record_h_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "hosp_d", col_name);
        record_d_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "death_d", col_name);
        record_gp_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "gp_d", col_name);
        record_bd_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "bd_d", col_name);
        record_dose_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "no_doses_d",  {"seed", "Total numver of doses"});
        record_q_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "qaly_d", {"seed", "Total QALY"});
        record_cp_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "cost_pro_d", {"seed", "Total Cost T"});
        record_ct_d.csv(get_ll::dout + get_ll::calout_name + "inter/" + prog_no + "/" + "cost_tre_d", {"seed", "Total Cost T"});
    }
}



#endif /* writing_h */
