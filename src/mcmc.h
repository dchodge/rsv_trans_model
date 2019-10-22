//
//  mcmc.h
//  flu_forecast
//
//  Created by David Hodgson on 17/06/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

#ifndef amh_h
#define amh_h


// Structure for monitoring the variables associated with the chains

/////
/// amh -> adaptive metropolis hasting. Everything in the amh namespace is used to find a new Markov chain element , including estimating the likelihood value
////

namespace amh
{
    struct amh_state_t  // a structure which contains all variables assocaited with an mcmc chain
    {
        //! Keep track of means of the mcmc samples
        Eigen::VectorXd means_parameters, prop_parameters, M_wg, M_v, prop_lprior_mwg, prop_llikelihood_mwg, curr_lprior_mwg, curr_llikelihood_mwg, my_acceptance_rate_mwg;
        MatrixXd emp_cov_matrix_fix, init_id_matrix, emp_cov_matrix, emp_cov_matrix1, emp_cov_matrix2, emp_cov_matrix3;
        
        // Important parameterrs in the run
        MatrixXd explore, result, posterior;
        VectorXd llikelihoods, llikelihoods_post, accept_ratios, accept_ratios_post, mean_parameters, curr_parameters;
        double curr_llikelihood, curr_lprior;
        
        bool OOB;           // indicator function to determine if sampled point in proposal is out of bounds
        bool OOB_indi;      // indicator function to determine if sampled point in proposal is out of bounds

        double adaptive_scaling;
        size_t no_accepted;     // number of accepted points in the markov chain
        size_t no_adaptive;     // number of adaptive points accepted in the markov chain
        bool adaptive_step;     // indicator function to detemrine if the step is adaptive

        double prop_llikelihood;    // holder for the proposed log-likelihood
        double prop_lprior;         // holder for the proposed log-prior
        double prop_hnr;                    // holder for the random probability samples to compare with the acceptance probability
        double my_acceptance_rate;          // holder for the acceptance probability
        double T;
        int L;
        double Ms;
        double lambda;
        int i, B;
        
        /// Stuff that is not changeable
        int k1, k2, Mburn, N0, N1, NK1, NK2, nk_1, nk_2;         //thinning
        int Nburn, Nrun, Nthin_all, Nthin_pos;
        double cov_scal_b;
        
        amh_state_t()
        {

        }
        
        amh_state_t(const int Nburn_t, const int Nrun_t, const int Nthinall_t, const int Nthinpos_t): Nburn(Nburn_t), Nrun(Nrun_t), Nthin_all(Nthinall_t), Nthin_pos(Nthinpos_t)
        {
            no_accepted = 0;
            no_adaptive = 0;
            adaptive_step = false;
            
            i = 0;
            N0 = 100;
            k1 = Nthin_all;
            k2 = Nthin_pos;
            Mburn = Nburn;
            N1 = Nrun;
            NK1 = Nrun/Nthin_all;
            NK2 = (Nrun-Nburn)/Nthin_pos - 1;
            cov_scal_b = 1.0;
        }
    };
    
    
    amh_state_t initialize(param::param_state_t& pars, int Nburn, int Nrun, int Nthin_all, int Nthin_pos)
    {
        amh_state_t state(Nburn, Nrun, Nthin_all, Nthin_pos);
        state.OOB = false;
        
        state.prop_parameters = Eigen::VectorXd::Zero(pars.dim_cal); state.curr_parameters = Eigen::VectorXd::Zero(pars.dim_cal); state.mean_parameters = Eigen::VectorXd::Zero(pars.dim_cal);
        
        state.explore = Eigen::MatrixXd::Zero(Nrun, pars.dim_cal);
        state.result = Eigen::MatrixXd::Zero(Nrun, pars.dim_cal);
        state.llikelihoods = Eigen::VectorXd::Zero(Nrun);
        state.accept_ratios = Eigen::VectorXd::Zero(Nrun);
        
        state.posterior = Eigen::MatrixXd::Zero(state.NK2, pars.dim_cal);
        state.llikelihoods_post = Eigen::VectorXd::Zero(state.NK2);
        state.accept_ratios_post = Eigen::VectorXd::Zero(state.NK2);
        
        state.M_wg = Eigen::VectorXd::Zero(pars.dim_cal);
        
        MatrixXd init_cov_matrix = MatrixXd::Zero(pars.dim_cal, pars.dim_cal);
        
        state.emp_cov_matrix_fix =
        MatrixXd::Zero(pars.dim_cal,pars.dim_cal);
        
        state.init_id_matrix = MatrixXd::Identity(pars.dim_cal,pars.dim_cal);
        
        state.emp_cov_matrix =
        MatrixXd::Zero(pars.dim_cal, pars.dim_cal);
        
        for(int i = 0; i < pars.dim_cal; i++)
        {
            state.emp_cov_matrix(i,i) = pars.support_range(i)/1000.0;
            state.emp_cov_matrix_fix(i,i) = pars.support_range(i)/1000.0;
        }
        
        state.lambda = 0.0;
        state.Ms =  0.0;
        
        state.nk_1 = 1;
        state.nk_2 = 0;
        
        return state;
    }
    
    // gain factors
    inline double gain_factor1(size_t n)
    {
        return 1.0/pow(n + 1.0, 1.0);
    }
    
    inline double gain_factor2(size_t n)
    {
        return 1.0/pow(n + 1.0, 0.5);
    }
    
    // update_cov -> update the adaptive covariance matrix according to the new point in markov chain
    void update_cov( amh_state_t& state)
    {
        double gf2 = gain_factor2(state.i);
        /*update of the variance-covariance matrix and the mean vector*/
        if (state.i==1)
        {
            state.mean_parameters = state.curr_parameters;
        }
        else
        {
            state.mean_parameters = state.mean_parameters + gf2*(state.curr_parameters-state.mean_parameters);
            state.emp_cov_matrix = state.emp_cov_matrix + gf2*((state.curr_parameters-state.mean_parameters)*((state.curr_parameters-state.mean_parameters).transpose())) - gf2*state.emp_cov_matrix;
        }
    }
    
    // accepted -> stochastic approximation to chain the scaling factor in the covariance matrix
    void accepted( amh_state_t& state, bool accepted, bool adaptive_step_size)
    {
        if (state.adaptive_step)
            ++state.no_adaptive;
        
        if (accepted)
            ++state.no_accepted;
        
        if(state.i < state.N0 || state.prop_hnr < 0.05)
        {
            state.lambda = state.lambda + gain_factor1(state.i)*(state.my_acceptance_rate - 0.234);
            if (state.lambda > state.cov_scal_b)
                state.lambda  = state.cov_scal_b;
        }
        else if (adaptive_step_size == true)
        {
            if (state.adaptive_step)
            {
                state.Ms = state.Ms + gain_factor1(state.no_adaptive)*(state.my_acceptance_rate - 0.234);
                if (state.Ms > state.cov_scal_b)
                    state.Ms  = state.cov_scal_b;
            }
        }
    }
    
    // everything in this namespace helps sample frmo the proposal distributions
    namespace get_proposal
    {
        // check_support -> check the sampled point is in the support of the truncated multivariate normal distribution
        void check_support(param::param_state_t& pars, VectorXd chain_pos)
        {

            /** Single infection parameters, must be positive **/
            for (int i = 0; i < pars.dim_cal; i++)
            {
                if (chain_pos(i) < pars.mapofParamLall.find(pars.paramFit[i])->second || chain_pos(i) > pars.mapofParamUall.find(pars.paramFit[i])->second){pars.OOS = true; return;}
            }
            pars.OOS = false;
        }
        
        // check_support -> sample a point from a random normal
        Eigen::VectorXd random_normal(amh::amh_state_t &state, param::param_state_t& pars)
        {
            VectorXd trans_temp(pars.dim_cal);
            int o = 0;
            
            EigenMultivariateNormal normX_solver(state.curr_parameters, exp(state.lambda)*state.emp_cov_matrix_fix, true, engine());
            trans_temp = normX_solver.samples(1);
            check_support(pars, trans_temp);
            while(pars.OOS)
            {
                o++;
                trans_temp = normX_solver.samples(1);
                check_support(pars, trans_temp);
                if (o == 100000)
                    cout << "Taking a long time to sample a point, consider changing the support of the parameters and/or transformaing the parameters" << '\n';
                if (o == 1000000)
                    cout << "Taking extremely long time to sample a point, you must change the support of the parameters and/or transform the parameters" << '\n';
            }
            return trans_temp;
        }
        
        // covar_adap -> sample a point from the adaptive covariance matrix
        Eigen::VectorXd covar_adap(amh::amh_state_t &state, param::param_state_t& pars)
        {
            VectorXd trans_temp(pars.dim_cal);
            int o = 0;
            
            EigenMultivariateNormal normX_solver(state.curr_parameters, exp(state.Ms)*state.emp_cov_matrix, true, engine());
            trans_temp = normX_solver.samples(1);
            check_support(pars, trans_temp);
            while(pars.OOS)
            {
                o++;
                trans_temp = normX_solver.samples(1);
                check_support(pars, trans_temp);
                if (o == 100000)
                {
                    EigenMultivariateNormal normX_solver(state.curr_parameters, exp(state.Ms)*state.emp_cov_matrix, true, engine());
                    cout << "Taking a long time to sample a point, consider changing the support of the parameters and/or transformaing the parameters" << '\n';
                }
                if (o == 1000000)
                    cout << "Taking extremely long time to sample a point, you must change the support of the parameters and/or transform the parameters" << '\n';
            }
            return trans_temp;
        }
        
        // sherlock -> using the sherlock proposal:
        void sherlock(amh::amh_state_t &state, param::param_state_t& pars)
        {
            VectorXd prop_temp;
            state.prop_hnr = uniform_dist(0, 1, 'r');
            
            if (state.i < state.N0){
                state.adaptive_step = false;
                state.prop_parameters = get_proposal::random_normal(state, pars);
            }
            else if (state.prop_hnr > 0.05){
                state.adaptive_step = true;
                state.prop_parameters = get_proposal::covar_adap(state, pars);
            }
            else {
                state.adaptive_step = false;
                state.prop_parameters = get_proposal::random_normal(state, pars);
            }
            
            param::update(pars, state.prop_parameters);
        }
    }
    
    // write_mcmc ->  write imporant variables for a chain
    void write_mcmc (amh::amh_state_t& mcmc_state,  param::param_state_t& pars, char F)
    {
        
        string FileMC_chain_explore = get_ll::model_t + "MC_chain_explore_" + F + ".txt";      //Output MC
        ofstream FileMC_chain_explore_out(FileMC_chain_explore.c_str());
        FileMC_chain_explore_out << mcmc_state.explore;
        FileMC_chain_explore_out.close();
        
        string FileMC_chain = get_ll::model_t + "MC_chain_" + F + ".txt";      //Output MC
        ofstream FileMC_chain_out(FileMC_chain.c_str());
        FileMC_chain_out << mcmc_state.result;
        FileMC_chain_out.close();
        
        string FileMC_chain_like = get_ll::model_t + "MC_chain_like_" + F + ".txt";      //Output MC
        ofstream FileMC_chain_like_out(FileMC_chain_like.c_str());
        FileMC_chain_like_out << mcmc_state.llikelihoods;
        FileMC_chain_like_out.close();
        
        string FileMC_acceptance = get_ll::model_t + "MC_acceptance_" + F + ".txt";      //Output MC
        ofstream FileMC_accpt(FileMC_acceptance.c_str());
        FileMC_accpt << mcmc_state.accept_ratios;
        FileMC_accpt.close();
        
        string FileMC_posterior = get_ll::model_t + "posterior_" + F + ".txt";      //Output MC
        ofstream FileMC_posterior_out(FileMC_posterior.c_str());
        FileMC_posterior_out << mcmc_state.posterior;
        FileMC_posterior_out.close();
        
        string FileMC_chain_like_con = get_ll::model_t + "MC_chain_like_post_" + F + ".txt";      //Output MC
        ofstream FileMC_chain_like_out_con(FileMC_chain_like_con.c_str());
        FileMC_chain_like_out_con << mcmc_state.llikelihoods_post;
        FileMC_chain_like_out_con.close();
        
        string FileMC_acceptance_con = get_ll::model_t + "MC_acceptance_post_" + F + ".txt";      //Output MC
        ofstream FileMC_accpt_con(FileMC_acceptance_con.c_str());
        FileMC_accpt_con << mcmc_state.accept_ratios_post;
        FileMC_accpt_con.close();
        
        
        /* Write out chain files */
        
        string FileMC_covar = get_ll::model_t + "covar_" + F + ".txt";      //Output MC
        ofstream FileMC_chain_covar(FileMC_covar.c_str());
        FileMC_chain_covar << mcmc_state.emp_cov_matrix_fix;
        FileMC_chain_covar << endl;
        FileMC_chain_covar << mcmc_state.emp_cov_matrix;
        FileMC_chain_covar << endl;;
        FileMC_chain_covar << mcmc_state.mean_parameters;
        FileMC_chain_covar.close();
    }
};

namespace mhp
{
    struct mhp_state_t // structure for parameter values which are used across all the chains
    {
        VectorXd S, no_accep, no_attemp, means_parameters, mean_parameters1, mean_parameters2, mean_parameters3;
        MatrixXd A, A_full, M_full, T_full, S_full, emp_cov_matrix, emp_cov_matrix1, emp_cov_matrix2, emp_cov_matrix3;
        num_vec T_ladder;
        int N, NK1;
        double t,  i, nk; //no of partition points
        
        mhp_state_t()
        {
            
            t = 0;
            nk = 1;
            i = 0;
            
        }
    };
    
    // initialise this first
    mhp_state_t pre_initialize(int no_Chain, int Nrun, int Nthin_all)
    {
        mhp_state_t mhp_state;
        mhp_state.N = no_Chain;
        mhp_state.NK1 = Nrun/Nthin_all;
        mhp_state.S = VectorXd::Zero( mhp_state.N - 1 );
        mhp_state.no_accep = VectorXd::Zero( mhp_state.N );
        mhp_state.no_attemp = VectorXd::Zero( mhp_state.N );
        mhp_state.A = VectorXd::Zero( mhp_state.N - 1 );
        
        mhp_state.A_full = MatrixXd::Zero( mhp_state.NK1, mhp_state.N - 1 );
        mhp_state.S_full = MatrixXd::Zero( mhp_state.NK1, mhp_state.N - 1 );
        mhp_state.T_full = MatrixXd::Zero( mhp_state.NK1, mhp_state.N  );
        mhp_state.M_full = MatrixXd::Zero( mhp_state.NK1, mhp_state.N  );
        
       // mhp_state.means_parameters = VectorXd::Zero( pars.dim_cal );
       // mhp_state.emp_cov_matrix = MatrixXd::Zero( pars.dim_cal, pars.dim_cal );

        //initialist temperature ladder
        for (int i = 0; i < mhp_state.N; i++)
            mhp_state.T_ladder.push_back(pow(10, 7.0*(i+1-1.0)/(mhp_state.N-1.0)));

        for (int i = 0; i < mhp_state.N-1 ;i++)
            mhp_state.A(i) = 0.234;
        
        for (int i = 0; i < mhp_state.N-1 ;i++)
        {
            mhp_state.S(i) = log(mhp_state.T_ladder[i+1] - mhp_state.T_ladder[i]);
        }
        
        mhp_state.T_full(0,0) = 1.0;
        for (int i = 0; i < mhp_state.N-1 ;i++)
        {
            mhp_state.T_full(0,i + 1) = mhp_state.T_full(0,i) + exp(mhp_state.S(i));
        }
        
        return mhp_state;
        
    }
    
    // initialise this second
    mhp_state_t post_initialize(mhp::mhp_state_t &mhp_state, std::vector<amh::amh_state_t> &mcmc_state)
    {

        for (int i = 0; i < mhp_state.N; i++)
        {
            mhp_state.M_full(0,i) = mcmc_state[i].Ms;
        }
        for(int i = 0; i < mhp_state.N-1; i++)
        {
            mhp_state.A_full(0,i) = mhp_state.A(i);
            mhp_state.S_full(0,i) = mhp_state.S(i);
        }
        
        return mhp_state;
    }
    
    // gain factor
    inline double gain_factor1(size_t n)
    {
        return 1.0/pow(n + 1.0, 1.0);
    }
    
    double gain_factor(mhp_state_t &mhp_state)
    {
        return 1.0/(pow(mhp_state.t + 1.0, 0.5));
    }
    
    // S_edit -> update S value
    void S_edit(mhp_state_t &mhp_state, int p)
    {
        mhp_state.S(p) = mhp_state.S(p) + gain_factor(mhp_state)*(mhp_state.A(p) - 0.234);
    }
    
    // T_edit -> update T value
    void T_edit(mhp_state_t &mhp_state, std::vector<amh::amh_state_t> &mcmc_state)
    {
        for (int i = 0; i < mhp_state.N-2; i++)
            mcmc_state[i+1].T = mcmc_state[i].T + exp(mhp_state.S(i));
    }
    
    // ratio_swap -> calculate the acceptance to swap chains
    inline double ratio_swap(amh::amh_state_t &state_mcmc_i, amh::amh_state_t &state_mcmc_j)
    {
        
        return min(1.0, exp((state_mcmc_j.curr_llikelihood + state_mcmc_j.curr_lprior - state_mcmc_i.curr_llikelihood - state_mcmc_i.curr_lprior)*(1.0/state_mcmc_i.T - 1.0/state_mcmc_j.T)) );
    }
    
    // save -> save elements of the mhp
    void save(mhp_state_t &mhp_state, std::vector<amh::amh_state_t> &mcmc_state)
    {
        for (int i = 0; i < mhp_state.N; i++)
        {
            mhp_state.T_full(mhp_state.nk,i) = mcmc_state[i].T;
            mhp_state.M_full(mhp_state.nk,i) = mcmc_state[i].Ms;
        }
        
        for(int i = 0; i < mhp_state.N-1; i++)
        {
            mhp_state.A_full(mhp_state.nk,i) = mhp_state.A(i);
            mhp_state.S_full(mhp_state.nk,i) = mhp_state.S(i);
        }
        mhp_state.nk++;
    }
    
    // update_swap -> swap adjacents chains according to the swapping probability
    void update_swap(mhp::mhp_state_t &mhp_state, std::vector<amh::amh_state_t> &mcmc_state, int i)
    {
        for (int j = 0; j < mhp_state.N; j++ )
        {
            int p = uniform_dist_disc(0, mhp_state.N - 2, 'r');
            int q = p + 1;
            mhp_state.no_attemp(p) ++;
            double r = mhp::ratio_swap(mcmc_state[p], mcmc_state[q]);
            
            if (uniform_dist(0,1,'r') < r)
            {
                mhp_state.t ++;
                mhp_state.no_accep(p) ++;
                
                std::swap(mcmc_state[p].curr_parameters, mcmc_state[q].curr_parameters);
                std::swap(mcmc_state[p].curr_lprior, mcmc_state[q].curr_lprior);
                std::swap(mcmc_state[p].curr_llikelihood, mcmc_state[q].curr_llikelihood);
            }
            
            mhp_state.A(p) = min(1.0, exp((mcmc_state[p+1].curr_llikelihood + mcmc_state[p+1].curr_lprior - mcmc_state[p].curr_llikelihood - mcmc_state[p].curr_lprior)*(1.0/mcmc_state[p].T - 1.0/mcmc_state[p+1].T)) );
            mhp::S_edit(mhp_state, p);
        }
    }
    
    // everything in here is used during the
    namespace run_MH
    {
        //  initialize_mcmc -> initialise the mcmc algorithm
        amh::amh_state_t initialize_mcmc(std::vector<std::string> paramFit_t, std::vector<param::param_state_t> &pars_r, int Nburn, int Nrun, int Nthin_all, int Nthin_pos)
        {
            
            //get initial sample
            param::param_state_t pars(paramFit_t);
            amh::amh_state_t mcmc_state = amh::initialize(pars, Nburn, Nrun, Nthin_all, Nthin_pos);

 /***** FUNCTIONS CALLING THE MODEL.H FILE */
            param::init_prior_sample(pars);
            mcmc_state.curr_parameters = param::parameterfitonly(pars);
            amh::get_proposal::check_support(pars, mcmc_state.curr_parameters);
            mcmc_state.curr_llikelihood = get_ll::get_likelihood(pars);

            //If initial parameters give inifinite likelihood then try again by sampling new parameters
            while (std::isinf(mcmc_state.curr_llikelihood) ==  true || std::isnan(mcmc_state.curr_llikelihood) ==  true || pars.OOS == true)
            {
/***** FUNCTIONS CALLING THE MODEL.H FILE */
                param::init_prior_sample(pars);
                mcmc_state.curr_parameters = param::parameterfitonly(pars);
                amh::get_proposal::check_support(pars, mcmc_state.curr_parameters);

                //check likelihood value
                mcmc_state.curr_llikelihood = get_ll::get_likelihood(pars);
                cout << "mcmc_state.curr_llikelihood: " << mcmc_state.curr_llikelihood << '\n';
                cout << "pars.OOS: " << pars.OOS << endl;
            }
            
/***** FUNCTIONS CALLING THE MODEL.H FILE */
            mcmc_state.curr_lprior = param::prior_dist_get(pars);
            pars_r.push_back(pars);
            
            //output values proposal state and the first entry of all the parameters
            cout << "mcmc_state.curr_llikelihood: " << mcmc_state.curr_llikelihood << '\n';
            cout << "mcmc_state.curr_lprior: " << mcmc_state.curr_lprior << '\n';
            cout << "mcmc_state.curr_parameters: " << mcmc_state.curr_parameters << '\n';
            
            //initial proposal state and the first entry of all the parameters
            mcmc_state.llikelihoods(0) = mcmc_state.curr_llikelihood;
            mcmc_state.result.row(0) = mcmc_state.curr_parameters;
            mcmc_state.explore.row(0) = mcmc_state.curr_parameters;
            
            return mcmc_state;
        }
        
        //  get_pop_llikelihood -> get the likelihood for each of the newly sampled points for each markov chain
        void get_pop_llikelihood(mhp::mhp_state_t &mhp_state, std::vector<amh::amh_state_t> &mcmc_state, std::vector<param::param_state_t> &pars)
        {
            //Sample point from proposal distribution, q() (see proposal.h)
            for (int i = 0; i < mhp_state.N; i++)
                amh::get_proposal::sherlock(mcmc_state[i], pars[i]);
            
            // Determine likelihood and prior values from model.h
            for (int i = 0; i < mhp_state.N; i++)
            {
/***** FUNCTIONS CALLING THE MODEL.H FILE */
                mcmc_state[i].prop_llikelihood = get_ll::get_likelihood(pars[i]);
                mcmc_state[i].prop_lprior = get_ll::get_lprior(mcmc_state[i].prop_llikelihood, pars[i]);
            }
        }
        
        // update_mcmc -> function to determine the correction factor due to asymetry arising from truncation
        double correction_mvn(VectorXd curr, VectorXd prop, MatrixXd Sigma, param::param_state_t& pars )
        {
            mvn::epmgp_t epmgp;
            epmgp = epmgp.initialise(pars.dim_cal);
            epmgp.K = Sigma;
            
            for (int i = 0; i < pars.dim_cal; i++) {
                epmgp.Ub(i) = pars.Upp_b_all[i]; epmgp.Lb(i) = pars.Low_b_all[i];
            }
            
            epmgp.m = curr; double retcur = epmgp_axis(epmgp);
            epmgp.m = prop; double retprop = epmgp_axis(epmgp);
            
            return retcur - retprop;
        }
        
        // update_mcmc -> function to calculate the acceptance probability of metropolis hasting algorithm
        void get_acceptance(mhp::mhp_state_t &mhp_state, std::vector<amh::amh_state_t> &mcmc_state, std::vector<param::param_state_t> &pars)
        {
            for (int i = 0; i < mhp_state.N; i++)
            {
                if(std::isinf(mcmc_state[i].prop_llikelihood) || std::isnan(mcmc_state[i].prop_llikelihood)  )
                    mcmc_state[i].my_acceptance_rate = -1;
                else
                {
                    if (mcmc_state[i].OOB_indi){
                        mcmc_state[i].my_acceptance_rate = 0.0; mcmc_state[i].OOB_indi = false;
                    }
                    else{
                        double correction = 0;
                        if (mcmc_state[i].i < mcmc_state[i].N0)
                        {
                            correction = correction_mvn(mcmc_state[i].curr_parameters, mcmc_state[i].prop_parameters, exp(mcmc_state[i].lambda)*mcmc_state[i].emp_cov_matrix_fix, pars[i] );
                            if (std::isinf(correction) ==  true || std::isnan(correction) ==  true)
                            {
                                mcmc_state[i].my_acceptance_rate = -1;
                                cout << "Correction funny" << endl;
                            }
                            else
                                mcmc_state[i].my_acceptance_rate = min(1.0, exp((mcmc_state[i].prop_llikelihood + mcmc_state[i].prop_lprior  - mcmc_state[i].curr_llikelihood - mcmc_state[i].curr_lprior + correction )/mcmc_state[i].T ));

                            
                        }
                        else if (mcmc_state[i].prop_hnr > 0.05)
                        {
                            correction =  correction_mvn(mcmc_state[i].curr_parameters, mcmc_state[i].prop_parameters, exp(mcmc_state[i].Ms)*mcmc_state[i].emp_cov_matrix_fix, pars[i] );
                            if (std::isinf(correction) ==  true || std::isnan(correction) ==  true)
                            {
                                mcmc_state[i].my_acceptance_rate = -1;
                                cout << "Correction funny" << endl;
                            }
                            else
                                mcmc_state[i].my_acceptance_rate = min(1.0, exp((mcmc_state[i].prop_llikelihood + mcmc_state[i].prop_lprior  - mcmc_state[i].curr_llikelihood - mcmc_state[i].curr_lprior + correction  )/mcmc_state[i].T ));
                        }
                        else
                        {
                            correction = correction_mvn(mcmc_state[i].curr_parameters, mcmc_state[i].prop_parameters, exp(mcmc_state[i].lambda)*mcmc_state[i].emp_cov_matrix_fix, pars[i] );
                            if (std::isinf(correction) ==  true || std::isnan(correction) ==  true)
                            {
                                mcmc_state[i].my_acceptance_rate = -1;
                                cout << "Correction funny" << endl;
                            }
                            else
                                mcmc_state[i].my_acceptance_rate = min(1.0, exp((mcmc_state[i].prop_llikelihood + mcmc_state[i].prop_lprior  - mcmc_state[i].curr_llikelihood - mcmc_state[i].curr_lprior + correction)/mcmc_state[i].T ));
                        }
                    }
                }
            }
        }
        
        // update_mcmc ->
        void update_mcmc (mhp::mhp_state_t &mhp_state, std::vector<amh::amh_state_t> &mcmc_state, std::vector<param::param_state_t> &pars)
        {
            // Update the chain number, adaptive covariance matrix.
            for (int i = 0; i < mhp_state.N; i++)
            {
                mcmc_state[i].i++;
                amh::update_cov(mcmc_state[i]);
            }
            
            run_MH::get_pop_llikelihood(mhp_state, mcmc_state, pars);
            run_MH::get_acceptance(mhp_state, mcmc_state, pars);
            
            for (int i = 0; i < mhp_state.N; i++)
            {
                if(mcmc_state[i].i%mcmc_state[i].k1==0)
                    mcmc_state[i].explore.row(mcmc_state[i].nk_1) = mcmc_state[i].prop_parameters;
                
                if(uniform_dist(0,1,'r') < mcmc_state[i].my_acceptance_rate)
                {
                    /** Update the acceptance rate */
                    amh::accepted(mcmc_state[i], true, true);
                    mcmc_state[i].curr_parameters = mcmc_state[i].prop_parameters;
                    mcmc_state[i].curr_llikelihood = mcmc_state[i].prop_llikelihood;
                    mcmc_state[i].curr_lprior = mcmc_state[i].prop_lprior;
                }
                else //if reject
                {
                    amh::accepted(mcmc_state[i], false, true);
                }
                
                //Save thinnned chain places (full model)
                if(mcmc_state[i].i%mcmc_state[i].k1==0)
                {
                    mcmc_state[i].llikelihoods(mcmc_state[i].nk_1) = mcmc_state[i].curr_llikelihood +  mcmc_state[i].curr_lprior;
                    mcmc_state[i].result.row(mcmc_state[i].nk_1) = mcmc_state[i].curr_parameters;
                    mcmc_state[i].accept_ratios(mcmc_state[i].nk_1) = double(mcmc_state[i].no_accepted)/mcmc_state[i].i;
                    mcmc_state[i].nk_1 ++;
                }
                
                //Save thinnned chains convergence (converged model model)
                if(mcmc_state[i].i%mcmc_state[i].k2==0 && mcmc_state[i].i > mcmc_state[i].Mburn)
                {
                    mcmc_state[i].llikelihoods_post(mcmc_state[i].nk_2) = mcmc_state[i].curr_llikelihood +  mcmc_state[i].curr_lprior;
                    mcmc_state[i].posterior.row(mcmc_state[i].nk_2) = mcmc_state[i].curr_parameters;
                    mcmc_state[i].accept_ratios_post(mcmc_state[i].nk_2) = double(mcmc_state[i].no_accepted)/mcmc_state[i].i;
                    
                    if (mcmc_state[i].i == mcmc_state[i].Mburn + mcmc_state[i].k2)
                        mcmc_state[i].mean_parameters = mcmc_state[i].curr_parameters;
                    else
                        mcmc_state[i].mean_parameters = mcmc_state[i].mean_parameters + 1.0/(mcmc_state[i].nk_2)*(mcmc_state[i].curr_parameters - mcmc_state[i].mean_parameters);
                    
                    mcmc_state[i].nk_2 ++;
                }
                
                // Output some informations
                if(mcmc_state[i].i%(mcmc_state[i].k1) == 0)
                {
                    cout << "chain cout: " << mcmc_state[i].i << "." <<'\t';
                    cout << "# accep ->" << double(mcmc_state[i].no_accepted) << "." <<'\t';
                    cout << "curr_like -> " << mcmc_state[i].curr_llikelihood << "." <<'\t';
                    cout << "curr_prior -> " << mcmc_state[i].curr_lprior << "." <<'\t';
                    cout << "T -> " << mcmc_state[i].T << "." <<'\t';
                    cout << "Ms -> " << mcmc_state[i].Ms << "." << '\t' << "curr -> " << mcmc_state[i].curr_parameters << '\n';
                    cout << "__________________________________________________________________________________________" << '\n';
                }
            }
        }
    }
    
    // write_mhp_state -> function to output information from the mhp state
    void write_mhp_state(mhp_state_t &mhp_state, char F)
    {
        string FileMC_chain_M_full = get_ll::model_t + "MC_chain_M_full_" + F + ".txt";      //Output MC
        ofstream FileMC_chain_M_full_out(FileMC_chain_M_full.c_str());
        FileMC_chain_M_full_out << mhp_state.M_full;
        FileMC_chain_M_full_out.close();
        
        string FileMC_chain_T_full = get_ll::model_t + "MC_chain_T_full_" + F + ".txt";      //Output MC
        ofstream FileMC_chain_T_full_out(FileMC_chain_T_full.c_str());
        FileMC_chain_T_full_out << mhp_state.T_full;
        FileMC_chain_T_full_out.close();
        
        string FileMC_chain_A_full = get_ll::model_t + "MC_chain_A_full_" + F + ".txt";      //Output MC
        ofstream FileMC_chain_A_full_out(FileMC_chain_A_full.c_str());
        FileMC_chain_A_full_out << mhp_state.A_full;
        FileMC_chain_A_full_out.close();
        
        string FileMC_chain_S_full = get_ll::model_t + "MC_chain_S_full_" + F + ".txt";      //Output MC
        ofstream FileMC_chain_S_full_out(FileMC_chain_S_full.c_str());
        FileMC_chain_S_full_out << mhp_state.S_full;
        FileMC_chain_S_full_out.close();
        
    }
};

#endif /* amh_h */
