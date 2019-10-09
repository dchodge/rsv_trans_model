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

// Define some important type defs
typedef Array<bool,Dynamic,1> ArrayXb;
typedef vector< double > num_vec;       //General-purpose numerical vector
typedef std::vector<std::string> str_vec;
typedef vector<vector<double> > vector2D;
typedef vector<vector<vector<double> > > vector3D;
typedef vector<vector<vector<vector<double> > > > vector4D;
typedef boost::mt19937 PRNG_s;
PRNG_s rng(engine()); //Generate non-static random numbers (pick different numbers from prior distribution each run)
//PRNG_s rng; //Generate non-static random numbers (pick different numbers from prior distribution each run)


double const ageGroup_nv[] = {0.0, 1.0/12.0, 2.0/12.0, 3.0/12.0, 4.0/12.0, 5.0/12.0, 6.0/12.0, 7.0/12.0, 8.0/12.0, 9.0/12.0, 10.0/12.0, 11.0/12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};  //Lower age group limits (upper limit not needed)
//{0.0, 0.5, 1.0, 5.0, 15.0, 45.0, 65.0};
//{0.0, 1.0/12.0, 2.0/12.0, 3.0/12.0, 4.0/12.0, 5.0/12.0, 6.0/12.0, 7.0/12.0, 8.0/12.0, 9.0/12.0, 10.0/12.0, 11.0/12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 15.0, 45.0, 65.0};  //Lower age group limits (upper limit not needed)

const int NoAgeG = (sizeof(ageGroup_nv)/sizeof(*ageGroup_nv));     //Derive number of age groups
const num_vec ageGroup(ageGroup_nv, ageGroup_nv + NoAgeG);                //Convert to vector
num_vec ageGroup_d(ageGroup_nv, ageGroup_nv + NoAgeG);
const int no_mo = 12; //Number of age groups <1
double ageGroupP[] = {0.0, 1.0/12.0, 2.0/12.0, 3.0/12.0, 4.0/12.0, 5.0/12.0, 6.0/12.0, 7.0/12.0, 8.0/12.0, 9.0/12.0, 10.0/12.0, 11.0/12.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 90.0};

using namespace Eigen;

/** ///////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////// 5. Distributions //////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////// **/


namespace distribution_def
{
    /**Hypergemoetric function 2F1*/
    double Hyper2F1(double a, double b, double c, double z) //Method 4b >> https://www.math.ucla.edu/~mason/research/pearson_final.pdf
    {
        double C = 1; double S = C;
        for (double j=0; j<100; j++)
        {
            C = C*(a+j)*(b+j)/(c+j)*z/(j+1);
            S = S + C;
        }
        return S;
    }
    
    /**Normal distribution*/
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
    
    double weibull_dist(double k, double t, char r, double x = 0)
    {
        if (r =='r')
        {boost::random::weibull_distribution<> W(k,t); return W(rng);}
        else if ( r =='m')
        {boost::math::weibull_distribution<> W(k,t); return mean(W);}
        else
        {boost::math::weibull_distribution<> W(k,t); return pdf(W,x);}
    }
    
    
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
    
    boost::random::normal_distribution<> sn(0,1);
    
    /**Log-Normal distribution*/
    double lognormal_dist(double mu1, double sig1, char r, double x = 0)
    {
        if (r =='r')
        {boost::random::lognormal_distribution<> n(mu1,sig1); return n(rng);}
        else if ( r =='m')
        {boost::math::lognormal_distribution<> n(mu1,sig1); return mean(n);}
        else
        {boost::math::lognormal_distribution<> n(mu1,sig1);return pdf(n,x);}
    }
    
    /**Beta distribution*/
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
    
    /**Gamma distribution*/
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
    
    double gamma_dist_inv(double k, double theta, char r, double x = 0)
    {
        if(k < 0) {cout << "Distribution OOB: k value in gamma is < 0" << endl;}
        else if (theta < 0) {cout << "Distribution OOB: theta value in gamma is < 0" << endl;}
        
        if ( r =='m')
        {boost::math::inverse_gamma_distribution<> g(k,theta); return mean(g);}
        else
        {boost::math::inverse_gamma_distribution<> g(k,theta); return pdf(g,x);}
    }
    
    /**Uniform distribution*/
    double uniform_dist(double a, double b, char r, double x = 0)
    {
        if (r =='r')
        {boost::random::uniform_real_distribution<> u(a,b); return u(rng);}
        else if ( r =='m')
        {boost::math::uniform_distribution<> u(a,b); return mean(u);}
        else
        {boost::math::uniform_distribution<> u(a,b); return pdf(u,x);}
    }
    
    double uniform_dist_disc(int a, int b, char r, double x = 0)
    {
        if (r =='r')
        {boost::random::uniform_int_distribution<> u(a,b); return u(rng);}
        else
            return 0;
    }
    
    
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
    
    /** Define special poisson distribution */
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
    
    long double stirl(double n)
    {
        double x = n + 1;
        return (x - 0.5)*log(x) - x + 0.5*log(2*PI) + 1.0/(12*x) - 1.0/(360.0*pow(x,3)); // https://en.wikipedia.org/wiki/Stirling%27s_approximation#Speed_of_convergence_and_error_estimates
    }
    
    double prior_on_incidence(double inc_tot_a, int a, num_vec popsize)
    {
        double ll = 0;
        if (a == 0)
            ll = log(gamma_dist(95.26, 0.00523172, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==1)
            ll = log(gamma_dist(109.191, 0.00449698, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==2)
            ll = log(gamma_dist(117.614, 0.00414962, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==3)
            ll = log(gamma_dist(123.156, 0.00394015, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==4)
            ll = log(gamma_dist(134.022, 0.00357981, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==5)
            ll = log(gamma_dist(139.718, 0.00341238, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==6)
            ll = log(gamma_dist(153.665, 0.00307141, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==7)
            ll = log(gamma_dist(156.759,  0.00298442, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==8)
            ll = log(gamma_dist(171.802, 0.0027138, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==9)
            ll = log(gamma_dist(168.809, 0.00273133, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==10)
            ll =  log(gamma_dist(189.384, 0.00242842, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==11)
            ll = log(gamma_dist(192.554, 0.00237422, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==12)
            ll = log(gamma_dist(287.976, 0.00152691, 'p', inc_tot_a/popsize[a]))/12.0;
        else if (a==13)
            ll = log(gamma_dist(429.08,  0.000948251, 'p', inc_tot_a/popsize[a]));
        else if (a==14)
            ll = log(gamma_dist(324.545, 0.0011647, 'p', inc_tot_a/popsize[a]));
        else if (a==15)
            ll = log(gamma_dist(196.136, 0.0017685, 'p', inc_tot_a/popsize[a]));
        else if (a==16)
            ll = log(gamma_dist(50.1489,  0.00561644, 'p', inc_tot_a/popsize[a]));
        else if (a==17)
            ll = log(gamma_dist(18.3679, 0.0126687, 'p', inc_tot_a/popsize[a]))*5;
        else if (a==18)
            ll = log(weibull_dist(4.28855, 0.240087, 'p', inc_tot_a/popsize[a]))*5;
        else if (a==19)
            ll = log(gamma_dist(16.0045, 0.0138897, 'p', inc_tot_a/popsize[a]))*10;
        else if (a==20)
            ll = log(weibull_dist(4.48801, 0.240849, 'p', inc_tot_a/popsize[a]))*10;
        else if (a==21)
            ll = log(weibull_dist(4.31245, 0.241728, 'p', inc_tot_a/popsize[a]))*10;
        else if (a==22)
            ll = log(weibull_dist(4.39624,0.243956, 'p', inc_tot_a/popsize[a]))*10;
        else if (a==23)
            ll = log(weibull_dist(4.50322, 0.239956, 'p', inc_tot_a/popsize[a]))*10;
        else if (a==24)
            ll = log(weibull_dist(4.30876, 0.240868, 'p', inc_tot_a/popsize[a]))*10;
        else
            cout << "Error" << endl;
        
        return ll;
    }
    
};

using namespace distribution_def;


#endif /* pre_h */
