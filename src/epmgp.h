
//  Created by David Hodgson on 21/06/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

#ifndef epmgp_h
#define epmgp_h

// Everything in this class is just for quick sampling from a multivariate normal distribution

namespace Eigen {
    namespace internal {
        //template<typename Scalar>
        struct scalar_normal_dist_op
        {
            //static boost::mt19937 rng;    // The uniform pseudo-random algorithm
            mutable boost::normal_distribution<double> norm;  // The gaussian combinator.
            
            EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)
            
            template<typename Index>
            inline const double operator() (Index, Index = 0) const { return norm(rng); }
            inline void seed(const uint64_t &s) { rng.seed(s); }
        };
        
        // template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;
        
        template<>
        struct functor_traits<scalar_normal_dist_op >
        { enum { Cost = 50 * NumTraits<double>::MulCost, PacketAccess = false, IsRepeatable = false }; };
    } // end namespace internal
    
    
    /**
     Find the eigen-decomposition of the covariance matrix
     and then store it for sampling from a multi-variate normal
     */
    //template<typename double>
    class EigenMultivariateNormal
    {
        MatrixXd _covar;
        MatrixXd _transform;
        VectorXd _mean;
        //internal::scalar_normal_dist_op<double> randN; // Gaussian functor
        internal::scalar_normal_dist_op randN; // Gaussian functor
        bool _use_cholesky;
        SelfAdjointEigenSolver<MatrixXd > _eigenSolver; // drawback: this creates a useless eigenSolver when using Cholesky decomposition, but it yields access to eigenvalues and vectors
        
    public:
        EigenMultivariateNormal(const VectorXd& mean, const MatrixXd& covar,
                                bool use_cholesky=false, const uint64_t &seed=boost::mt19937::default_seed)
        :_use_cholesky(use_cholesky)
        {
            randN.seed(seed);
            setMean(mean);
            setCovar(covar);
        }
        
        void setMean(const VectorXd& mean) { _mean = mean; }
        void setCovar(const MatrixXd& covar)
        {
            _covar = covar;
            
            // Assuming that we'll be using this repeatedly,
            // compute the transformation matrix that will
            // be applied to unit-variance independent normals
            
            if (_use_cholesky)
            {
                Eigen::LLT<Eigen::MatrixXd > cholSolver(_covar);  //Decompose covar in the cholesky decomposition
                // We can only use the cholesky decomposition if
                // the covariance matrix is symmetric, pos-definite.
                // But a covariance matrix might be pos-semi-definite.
                // In that case, we'll go to an EigenSolver
                if (cholSolver.info()==Eigen::Success)
                {
                    // Use cholesky solver
                    _transform = cholSolver.matrixL(); //Retrieve matrix L from decomposition
                }
                else
                {
                    _eigenSolver = SelfAdjointEigenSolver<MatrixXd >(_covar);  //computes Eigenvalues and vectors of self adjoint matrices
                    _transform = _eigenSolver.eigenvectors()*_eigenSolver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
                }
            }
            else
            {
                _eigenSolver = SelfAdjointEigenSolver<MatrixXd >(_covar);  //computes Eigenvalues and vectors of self adjoint matrices
                _transform = _eigenSolver.eigenvectors()*_eigenSolver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal(); //find max then square root of eigen values x times vectors and take diganoal
            }
        }
        
        /// Draw nn samples from the gaussian and return them
        /// as columns in a Dynamic by nn matrix
        Matrix<double,Dynamic,-1> samples(int nn)
        {
            return (_transform * Matrix<double,Dynamic,-1>::NullaryExpr(_covar.rows(),nn,randN)).colwise() + _mean;
        }
    }; // end class EigenMultivariateNormal
    
} // end namespace Eigen






using namespace Eigen;

// Everything in this class is just for evaluating a point in a multivaraite truncated normal distribution

namespace mvn
{
    
    struct epmgp_t
    {
        // Define stuff
        int dim_cal;
        MatrixXd K; //K, the positive semi-definite covariance matrix of the Gaussian
        VectorXd Lb, Ub, m; //the mean of the multivariate Gaussian, the lower/upper bound of the hyper-rectangle
        
        VectorXd muHat, sigmaHat, logZHat, ZHat; //Inputs/outputs of truncNormMoments
        
        VectorXd mu; // the mean of the truncated Gaussian (calculated by EPGCD)
        MatrixXd Sigma; //the covariance of the truncated Gaussian (calculated by EPGCD)
        
        /********************************/
        /**  convergence criteria for stopping EPGCD **/
        /********************************/
        double epsCon = 1e-5;
        
        epmgp_t()
        {
            
        }
        
        mvn::epmgp_t initialise(int dim)
        {
            mvn::epmgp_t epmgp_s;
            epmgp_s.dim_cal = dim;
            epmgp_s.Lb = VectorXd::Zero(epmgp_s.dim_cal);
            epmgp_s.Ub = VectorXd::Zero(epmgp_s.dim_cal);
            epmgp_s.m = VectorXd::Zero(epmgp_s.dim_cal);
            
            epmgp_s.muHat = VectorXd::Zero(epmgp_s.dim_cal);
            epmgp_s.sigmaHat = VectorXd::Zero(epmgp_s.dim_cal);
            epmgp_s.logZHat = VectorXd::Zero(epmgp_s.dim_cal);
            epmgp_s.ZHat = VectorXd::Zero(epmgp_s.dim_cal);
            
            epmgp_s.mu = VectorXd::Zero(epmgp_s.dim_cal);
            
            epmgp_s.K = MatrixXd::Zero(epmgp_s.dim_cal,epmgp_s.dim_cal);
            
            epmgp_s.Sigma = MatrixXd::Zero(epmgp_s.dim_cal,epmgp_s.dim_cal);
            return epmgp_s;
        }
        
    };
    
    inline int sign(double val) {
        return (0 < val) - (val < 0);
    }
    
    inline double erfcx(double b){
        if (b > 25)
            return 1.0/(sqrt(PI)*b);
        else
            return exp(pow(b,2))*erfc(b);
    }
    
    void truncNormMoments(mvn::epmgp_t&  epmgp_s)
    {
        // Variable used in this space only
        double lb, ub, mu, sigma;
        double a,b;
        double logZhatOtherTail, logZhat, Zhat, muHat, sigmaHat;
        double meanConst, varConst;
        double maxab, minab;
        
        
        /********************************/
        /** Main loop **/
        /********************************/
        /**Calculate for each input univariate Gaussian */
        for (int i = 0; i < epmgp_s.dim_cal; i++)
        {
            lb = epmgp_s.Lb(i);  ub = epmgp_s.Ub(i);
            mu = epmgp_s.muHat(i);  sigma = epmgp_s.sigmaHat(i);
            
            if (lb > ub)
                cout << "Bounds are incorrect" << endl;
            
            /********************************/
            /** Establish bounds  **/
            /********************************/
            a = (lb - mu)/sqrt(2*sigma); b = (ub - mu)/sqrt(2*sigma);
            

            if (std::isinf(a) && std::isinf(b)) /** Problem case A*/
            {
                cout << "Error. Problem case A" << endl;
                if (sign(a) == sign(b))
                {
                    epmgp_s.logZHat(i) = -INFINITY; epmgp_s.ZHat(i) = 0;
                    epmgp_s.muHat(i) = a; epmgp_s.sigmaHat(i) = 0;
                    continue;
                }
                else
                { 
                    epmgp_s.logZHat(i) = 0; epmgp_s.ZHat(i) = 1;
                    epmgp_s.muHat(i) = mu; epmgp_s.sigmaHat(i) = sigma;
                    continue;
                }
            }
            else if (a > b)             /** Problem case B*/
            {
                cout << "Error. Problem case B" << endl;
                epmgp_s.logZHat(i) = -INFINITY; epmgp_s.ZHat(i) = 0;
                epmgp_s.muHat(i) = mu; epmgp_s.sigmaHat(i) = 0;
                continue;
            }
            else if (a==-INFINITY)  /** Real case A when lower limit infinite */
            {
                if (b > 25)
                {
                    //logZhatOtherTail = log(0.5) + pow(b,2) + log(erfc(b)) - pow(b,2);
                    logZhatOtherTail = log(0.5) + log(erfc(b));
                    logZhat = log1p(-exp(logZhatOtherTail));
                }
                else
                    logZhat = log(0.5) + log(erfcx(-b)) - pow(b,2);
                
                meanConst = -2./erfcx(-b);
                varConst = -2./erfcx(-b)*(ub + mu);
            }
            else if (b==INFINITY) /** Real case B when lower limit infinite */
            {
                if (a < -25)
                {
                    //logZhatOtherTail = log(0.5) + log(erfcx(-a)) - pow(a,2);
                    logZhatOtherTail = log(0.5) + log(erfc(-a));
                    logZhat = log1p(-exp(logZhatOtherTail));
                }
                else
                    logZhat = log(0.5) + log(erfcx(a)) - pow(a,2);
                
                meanConst = 2./erfcx(a);
                varConst = 2./erfcx(a)*(lb + mu);
            }
            else    /** Real case C */
            {
                if (sign(a) == sign(b))
                {
                    maxab = fmax(abs(a),abs(b));
                    minab = fmin(abs(a),abs(b));
                    logZhat = log(0.5) - pow(minab,2) + log( abs( exp(-(pow(maxab,2)-pow(minab,2)))*erfcx(maxab) - erfcx(minab)) );
                    
                    meanConst = 2*sign(a)*(1./((erfcx(abs(a)) - exp(pow(a,2)-pow(b,2))*erfcx(abs(b)))) - 1/((exp(pow(b,2)-pow(a,2))*erfcx(abs(a)) - erfcx(abs(b)))));
                    varConst =  2*sign(a)*((lb+mu)/((erfcx(abs(a)) - exp(pow(a,2)-pow(b,2))*erfcx(abs(b)))) - (ub+mu)/((exp(pow(b,2)-pow(a,2))*erfcx(abs(a)) - erfcx(abs(b)))));
                }
                else
                {
                    if (abs(b) > abs(a) || abs(b) == abs(a))
                    {
                        if (a >= -25) /** Do things normally*/
                        {

                            logZhat = log(0.5) - pow(a,2) + log( erfcx(a) - exp(-(pow(b,2) - pow(a,2)))*erfcx(b));
                            meanConst = 2*(1.0/((erfcx(a) - exp(pow(a,2)-pow(b,2))*erfcx(b))) - 1.0/((exp(pow(b,2)-pow(a,2))*erfcx(a) - erfcx(b))));
                            varConst = 2*((lb+mu)/((erfcx(a) - exp(pow(a,2)-pow(b,2))*erfcx(b))) - (ub+mu)/((exp(pow(b,2)-pow(a,2))*erfcx(a) - erfcx(b))));
                        }
                        else /** a very small and  will be unstable */
                        {
                            logZhat = log(0.5) + log( 2 - erfc(b) - erfc(-a) ); ///
                            meanConst = 2*(1.0/(exp(pow(a,2))*(erfc(a) - erfc(b))) - 1.0/(exp(pow(b,2))*(2 - erfc(b))));
                            //meanConst = 2*( - 1.0/(exp(pow(b,2))*2 - erfcx(b)));
                            //varConst = 2*((lb+mu)/((erfcx(a) - exp(pow(a,2)-pow(b,2))*erfcx(b))) - (ub+mu)/(exp(pow(b,2))*2 - erfcx(b)));
                            varConst = 2*((lb+mu)/((exp(pow(a,2))*(erfc(a) - erfc(b)))) - (ub+mu)/(exp(pow(b,2))*(2 - erfc(b))));
                            //varConst = 2*(- (ub+mu)/(exp(pow(b,2))*2 - erfcx(b)));
 
                        }
                    }
                    else
                    {
                        if (b <= 25) /*Do things normally*/
                        {
                            logZhat = log(0.5) - pow(b,2) + log( erfcx(-b) - exp(-(pow(a,2) - pow(b,2)))*erfcx(-a));
                            meanConst = -2*(1/((erfcx(-a) - exp(pow(a,2)-pow(b,2))*erfcx(-b))) - 1/((exp(pow(b,2)-pow(a,2))*erfcx(-a) - erfcx(-b))));
                            varConst = -2*((lb+mu)/((erfcx(-a) - exp(pow(a,2)-pow(b,2))*erfcx(-b))) - (ub+mu)/((exp(pow(b,2)-pow(a,2))*erfcx(-a) - erfcx(-b))));
                        }
                        else /*b very small and calculations will be unstable */
                        {
                            //logZhat = log(0.5) + log( 2 - exp(-pow(a,2))*erfcx(-a) - exp(-pow(b,2))*erfcx(b) );
                            logZhat = log(0.5) + log( 2 - erfc(-a) - erfc(b) );
                            meanConst = -2*(1/((erfc(-a) - 2)*exp(pow(a,2))) - 1/(exp(pow(b,2))*(erfc(-a) - erfc(-b))));
                            //meanConst = -2*(1/(erfcx(-a) - exp(pow(a,2))*2));
                            //varConst = -2*((lb + mu)/(erfcx(-a) - exp(pow(a,2))*2) - (ub + mu)/(exp(pow(b,2)-pow(a,2))*erfcx(-a) - erfcx(-b)));
                            varConst = -2*((lb + mu)/(exp(pow(a,2))*(erfc(-a) -  2)) - (ub + mu)/(exp(pow(b,2))*(erfc(-a) - erfc(-b))));
                            //varConst = -2*((lb + mu)/(erfcx(-a) - exp(pow(a,2))*2));
                        }
                    }
                }
            }
            
            Zhat = exp(logZhat);
            muHat = mu + meanConst*sqrt(sigma/(2*PI));
            sigmaHat = sigma + varConst*sqrt(sigma/(2*PI)) + pow(mu,2) - pow(muHat,2);
            
            epmgp_s.logZHat(i) = logZhat;        epmgp_s.ZHat(i) = Zhat;
            epmgp_s.muHat(i) = muHat;            epmgp_s.sigmaHat(i) = sigmaHat;
        }
    }
    
    double epmgp_axis(epmgp_t& epmgp_s)
    {
        // Initialise containers
        VectorXd deltatauSite, nuSite, nuCavity;
        VectorXd tauSite, tauCavity;
        VectorXd KinvM, muLast;
        double logZEP = 0;
        bool converged;
        double k, lZ1, lZ2, lZ3, lZ4;
        VectorXd Sigmadiag;
        VectorXd ones;
        MatrixXd temp1;
        VectorXd temp2;
        
        MatrixXd SsiteHalf;
        MatrixXd Id(epmgp_s.dim_cal, epmgp_s.dim_cal);
        MatrixXd cholmat;
        MatrixXd L;
        MatrixXd V;
        
        Id.setIdentity();
        
        ones = VectorXd::Ones(epmgp_s.dim_cal);
        
        nuSite = VectorXd::Zero(epmgp_s.dim_cal);       tauSite = VectorXd::Zero(epmgp_s.dim_cal);
        nuCavity = VectorXd::Zero(epmgp_s.dim_cal);     tauCavity = VectorXd::Zero(epmgp_s.dim_cal);
        KinvM = VectorXd::Zero(epmgp_s.dim_cal);        muLast = VectorXd::Ones(epmgp_s.dim_cal);
        deltatauSite = VectorXd::Zero(epmgp_s.dim_cal);
        SsiteHalf = MatrixXd::Zero(epmgp_s.dim_cal,epmgp_s.dim_cal);
        temp1 = MatrixXd::Zero(epmgp_s.dim_cal,epmgp_s.dim_cal);
        temp2 = VectorXd::Zero(epmgp_s.dim_cal);
        
        /********************************/
        /**  convergence criteria for stopping EPGCD **/
        /********************************/
        
        
        /********************************/
        /** initialize the distribution q(x)**/
        /********************************/
        epmgp_s.Sigma = epmgp_s.K;
        epmgp_s.mu = (epmgp_s.Lb + epmgp_s.Ub)/2.0; // works fine
        for (int i = 0; i < epmgp_s.dim_cal ; i++)
        {
            if (std::isinf(epmgp_s.mu(i)))
                epmgp_s.mu(i) = sign(epmgp_s.mu(i))*100;
        }
        
        /********************************/
        /** Precalculate a few quantities **/
        /********************************/
        KinvM = epmgp_s.K.ldlt().solve(epmgp_s.m);  // Works for identity matrix
        logZEP = INFINITY;
        muLast = -INFINITY*muLast; // works
        converged = 0;
        k = 1;
        
        /********************************/
        /** MAIN axisEPMGP ALGORITHM LOOP **/
        /********************************/
        while(!converged)
        {
            /********************************/
            /** make the cavity distribution **/
            /********************************/
            Sigmadiag = epmgp_s.Sigma.diagonal(); // makes vector of diagonal
            tauCavity = ones.cwiseQuotient(Sigmadiag) - tauSite; // work
            nuCavity = epmgp_s.mu.cwiseQuotient(Sigmadiag) - nuSite; // work
            
            /********************************/
            /** Compute moments using truncated normals **/
            /********************************/
            epmgp_s.muHat = nuCavity.cwiseQuotient(tauCavity); // work
            epmgp_s.sigmaHat = ones.cwiseQuotient(tauCavity); // work
            truncNormMoments(epmgp_s);
            
            /********************************/
            /** Update the site parameters **/
            /********************************/
            deltatauSite = ones.cwiseQuotient(epmgp_s.sigmaHat) - tauCavity - tauSite;
            tauSite = tauSite + deltatauSite;
            nuSite = epmgp_s.muHat.cwiseQuotient(epmgp_s.sigmaHat) - nuCavity;
            
            /// tauSite is probably nonnegative. If a result is negative, it is either numerical precision (then we change to 0) or an error in the algorithm
            for (int i = 0; i < epmgp_s.dim_cal; i++)
            {
                if (tauSite(i) < 0 && tauSite(i) > -1e-6)
                {
                    tauSite(i) = 0;
                }
                else if  (tauSite(i) < 0)
                {
                    tauSite(i) = 0;
                    cout << "Error the code is wrong please check moments" << endl;
                }
            }
            
            /********************************/
            /** Having now iterated through all sites, update q(x) (Sigma and mu) **/
            /********************************/
            SsiteHalf = tauSite.cwiseSqrt().asDiagonal(); //works fine
            cholmat = Id + SsiteHalf*epmgp_s.K*SsiteHalf;
            Eigen::LLT<Eigen::MatrixXd > cholSolver(cholmat);
            L = cholSolver.matrixL().transpose();
            V = (L.transpose()).colPivHouseholderQr().solve(SsiteHalf*epmgp_s.K);
            epmgp_s.Sigma = epmgp_s.K - V.transpose()*V;

            epmgp_s.mu = epmgp_s.Sigma*(nuSite + KinvM);
            
            /********************************/
            /** Check convergence criteria **/
            /********************************/
            if ((muLast-epmgp_s.mu).norm() < epmgp_s.epsCon)
                converged = 1;
            else
                muLast = epmgp_s.mu;
            k++;
        }
        /********************************/
        /** Finished algorithm main loop **/
        /********************************/

        /********************************/
        /** Compute logZEP, the EPGCD cumulative density, from q(x) **/
        /********************************/
        if (logZEP != -INFINITY) // ensure that no 0 moment was found, in which case this calculation is not needed
        {
            // find lZ1
            lZ1 = 0;
            temp1 = ones + tauSite.cwiseQuotient(tauCavity);
            temp2 = L.diagonal();
            for (int i = 0; i < epmgp_s.dim_cal; i++)
                lZ1 += 0.5*log(temp1(i)) - log(temp2(i));
            
            // find lZ2
            temp1 = (ones.cwiseQuotient(tauCavity + tauSite)).asDiagonal();
            lZ2 = 0.5*(nuSite - tauSite.cwiseProduct(epmgp_s.m)).transpose()*(epmgp_s.Sigma - temp1)*(nuSite - tauSite.cwiseProduct(epmgp_s.m));
            
            // find lZ3
            temp1 = (tauSite + tauCavity).asDiagonal();
            temp2 = (tauSite.cwiseProduct(nuCavity)).cwiseQuotient(tauCavity) - 2*nuSite;
            lZ3 = 0.5*nuCavity.transpose()*(temp1.colPivHouseholderQr().solve(temp2));
            
            // find lZ4
            temp1 = (tauSite + tauCavity).asDiagonal();
            temp2 = (tauSite.cwiseProduct(epmgp_s.m) - 2*nuSite);
            lZ4 = - 0.5*(tauCavity.cwiseProduct(epmgp_s.m)).transpose()*(temp1.colPivHouseholderQr().solve(temp2));
            logZEP = lZ1 + lZ2 + lZ3 + lZ4 + epmgp_s.logZHat.sum();
        }
        return logZEP;
        
    }
    
}

#endif /* epmgp_h */
