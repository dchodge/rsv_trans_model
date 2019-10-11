//
//  pre.h
//  flu_forecast
//
//  Created by David Hodgson on 13/06/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

#ifndef pre_h
#define pre_h

/** ///////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////// 1. Preamble //////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////// **/

using namespace std;
using namespace boost::math;

// Random seeding stuff
std::random_device dev;
std::mt19937 engine(dev());
typedef boost::mt19937 PRNG_s;
PRNG_s rng(engine());

// Important type defs
typedef Array<bool,Dynamic,1> ArrayXb;      // General-purpose boolean vector (eigen)
typedef vector< double > num_vec;           // General-purpose numerical vector
typedef std::vector<std::string> str_vec;   // General-purpose string vector
typedef vector<vector<double> > vector2D;   // General-purpose 2-d numerical vector
typedef vector<vector<vector<double> > > vector3D; // General-purpose 3-d numerical vector
typedef vector<vector<vector<vector<double> > > > vector4D; // General-purpose 4-d numerical vector

// Age stratification
double const ageGroup_nv[] = {0.0, 1.0/12.0, 2.0/12.0, 3.0/12.0, 4.0/12.0, 5.0/12.0, 6.0/12.0, 7.0/12.0, 8.0/12.0, 9.0/12.0, 10.0/12.0, 11.0/12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};    // Age stratification in non-vector
const int NoAgeG = (sizeof(ageGroup_nv)/sizeof(*ageGroup_nv));    //Derive number of age groups from non-vector
const num_vec ageGroup(ageGroup_nv, ageGroup_nv + NoAgeG);        // Convert to vector
num_vec ageGroup_d(ageGroup_nv, ageGroup_nv + NoAgeG);            // Derive number of age groups from vector
double ageGroupP[] = {0.0, 1.0/12.0, 2.0/12.0, 3.0/12.0, 4.0/12.0, 5.0/12.0, 6.0/12.0, 7.0/12.0, 8.0/12.0, 9.0/12.0, 10.0/12.0, 11.0/12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 90.0};   //Age groups with upper limit

using namespace Eigen;

/** ///////////////////////////////////////////////////////////////////////////////////
 ///////// 2. Probability distributions + stirling approximation //////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////// **/

namespace distribution_def
{
    // binomial distribution
    double binomial_dist(int n, double p, char r, int x = 0)
    {
        if(p < 0) {cout << "Distribution OOB: probability value in binomial is < 0" << endl;}
        if (r =='r')
        {boost::random::binomial_distribution<> B(n,p); return B(rng);}
        else if ( r =='m')
        {boost::math::binomial_distribution<> B(n,p); return mean(B);}
        else
        {boost::math::binomial_distribution<> B(n,p);return pdf(B,x);}
    }
    
    // wiebull distribution
    double weibull_dist(double k, double t, char r, double x = 0)
    {
        if (r =='r')
        {boost::random::weibull_distribution<> W(k,t); return W(rng);}
        else if ( r =='m')
        {boost::math::weibull_distribution<> W(k,t); return mean(W);}
        else
        {boost::math::weibull_distribution<> W(k,t); return pdf(W,x);}
    }
    
    // poisson distribution
    double poisson_dist(double lambda, char r, int x = 0)
    {
        if(lambda < 0) {cout << "Distribution OOB: probability value in binomial is < 0" << endl;}
        if (r =='r')
        {boost::random::poisson_distribution<> P(lambda); return P(rng);}
        else if ( r =='m')
        {boost::math::poisson_distribution<> P(lambda); return mean(P);}
        else
        {boost::math::poisson_distribution<> P(lambda);return pdf(P,x);}
    }
    
    // cdf of poisson distribution
    double poisson_cdf(double l, double a, double x)
    {
        if( l == 0.0 || a == 0.0)
        {
            boost::math::poisson_distribution<> p(0.000001); return cdf(p,x);
        }
        else
        {
            boost::math::poisson_distribution<> p(l*a); return cdf(p,x);
        }
    }
    
    // normal distribution
    double normal_dist(double mu1, double sig1, char r, double x = 0)
    {
        if(sig1 < 0) {cout << "Distribution OOB: standard deviation value in normal is < 0" << endl;}
        if (r =='r')
        {boost::random::normal_distribution<> n(mu1,sig1); return n(rng);}
        else if ( r =='m')
        {boost::math::normal_distribution<> n(mu1,sig1); return mean(n);}
        else if ( r =='c')
        {boost::math::normal_distribution<> n(mu1,sig1);return cdf(n,x);}
        else
        {boost::math::normal_distribution<> n(mu1,sig1);return pdf(n,x);}
    }
    
    // lognormal distribution
    double lognormal_dist(double mu1, double sig1, char r, double x = 0)
    {
        if (r =='r')
        {boost::random::lognormal_distribution<> n(mu1,sig1); return n(rng);}
        else if ( r =='m')
        {boost::math::lognormal_distribution<> n(mu1,sig1); return mean(n);}
        else
        {boost::math::lognormal_distribution<> n(mu1,sig1);return pdf(n,x);}
    }
    
    //  beta distribution
    double beta_dist(double alpha, double beta, char r, double x = 0)
    {
        if(alpha < 0) {cout << "Distribution OOB: alpha value in beta is < 0" << endl;}
        else if (beta < 0) {cout << "Distribution OOB: beta value in beta is < 0" << endl;}
        if (r =='r')
        {boost::random::beta_distribution<> b(alpha,beta); return b(rng);}
        else if ( r =='m')
        {boost::math::beta_distribution<> b(alpha,beta); return mean(b);}
        else
        {boost::math::beta_distribution<> b(alpha,beta); return pdf(b,x);}
    }
    
    // gamma distribution
    double gamma_dist(double k, double theta, char r, double x = 0)
    {
        if(k < 0) {cout << "Distribution OOB: k value in gamma is < 0" << endl;}
        else if (theta < 0) {cout << "Distribution OOB: theta value in gamma is < 0" << endl;}
        
        if (r =='r')
        {boost::random::gamma_distribution<> g(k,theta); return g(rng);}
        else if ( r =='m')
        {boost::math::gamma_distribution<> g(k,theta); return mean(g);}
        else
        {boost::math::gamma_distribution<> g(k,theta); return pdf(g,x);}
    }
    
    // inverse gamma distribution
    double gamma_dist_inv(double k, double theta, char r, double x = 0)
    {
        if(k < 0) {cout << "Distribution OOB: k value in gamma is < 0" << endl;}
        else if (theta < 0) {cout << "Distribution OOB: theta value in gamma is < 0" << endl;}
        
        if ( r =='m')
        {boost::math::inverse_gamma_distribution<> g(k,theta); return mean(g);}
        else
        {boost::math::inverse_gamma_distribution<> g(k,theta); return pdf(g,x);}
    }
    
    // uniform distribution (continuous)
    double uniform_dist(double a, double b, char r, double x = 0)
    {
        if (r =='r')
        {boost::random::uniform_real_distribution<> u(a,b); return u(rng);}
        else if ( r =='m')
        {boost::math::uniform_distribution<> u(a,b); return mean(u);}
        else
        {boost::math::uniform_distribution<> u(a,b); return pdf(u,x);}
    }
    
    // uniform distribution (discrete)
    double uniform_dist_disc(int a, int b, char r, double x = 0)
    {
        if (r =='r')
        {boost::random::uniform_int_distribution<> u(a,b); return u(rng);}
        else
            return 0;
    }
    
    // inverse uniform distribution
    double uniform_dist_inv(double a, double b, char r, double x = 0)
    {
        if ( r =='m')
        {
            return (log(1.0/a) - log(1.0/b))/(b-a);
        }
        else
        {
            if (x < 1.0/b ||  x > 1.0/a)
                return 0;
            else
                return 1.0/((b-a)*x*x);
        }
    }
    
    // striling approximation
    long double stirl(double n)
    {
        double x = n + 1;
        return (x - 0.5)*log(x) - x + 0.5*log(2*PI) + 1.0/(12*x) - 1.0/(360.0*pow(x,3)); // https://en.wikipedia.org/wiki/Stirling%27s_approximation#Speed_of_convergence_and_error_estimates
    }
   
};

using namespace distribution_def;


#endif /* pre_h */
