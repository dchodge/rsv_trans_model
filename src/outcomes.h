//
//  outcomes.h
//  calibration
//
//  Created by David Hodgson on 16/07/2019.
//  Copyright © 2019 David Hodgson. All rights reserved.
//

#ifndef outcomes_h
#define outcomes_h

double ep = 1e-6;
num_vec bisection_Gamma(num_vec ci_to_fit)
{
    double mu = ci_to_fit[0];
    double ll = ci_to_fit[1];
    double ul = ci_to_fit[2];
    num_vec out;
    double a = 1e-6, b = 10000;
    double x = a;
    while ((b-a) >= ep)
    {
        x = (a+b)/2.0;
        boost::math::gamma_distribution<> Fx(x,mu/x); boost::math::gamma_distribution<> Fa(a,mu/a);
        
        if ((cdf(Fx, ul) - cdf(Fx, ll) - 0.95)*(cdf(Fa, ul) - cdf(Fa, ll) - 0.95) < 0.0)
            b = x;
        else
            a = x;
    }
    boost::math::gamma_distribution<> Fx(x,mu/x);
    double err = (abs(cdf(Fx, ll) - 0.025) + abs(cdf(Fx, ul) - 0.975));
    out.push_back(1.0);  out.push_back(x); out.push_back(mu/x); out.push_back(err);
    return out;
}

num_vec bisection_LN(num_vec ci_to_fit)
{
    double mu = ci_to_fit[0];
    double ll = ci_to_fit[1];
    double ul = ci_to_fit[2];
    num_vec out;
    double a = 1e-6, b = 10000;
    double x = a;
    while ((b-a) >= ep)
    {
        x = (a+b)/2;
        boost::math::lognormal_distribution<> Fx(log(mu),x); boost::math::lognormal_distribution<> Fa(log(mu),a);
        
        if ((cdf(Fx, ul) - cdf(Fx, ll) - 0.95)*(cdf(Fa, ul) - cdf(Fa, ll) - 0.95) < 0.0)
            b = x;
        else
            a = x;
    }
    boost::math::lognormal_distribution<> Fx(log(mu),x);
    double err = (abs(cdf(Fx, ll) - 0.025) + abs(cdf(Fx, ul) - 0.975));
    out.push_back(2.0);  out.push_back(log(mu)); out.push_back(x); out.push_back(err);
    return out;
}

num_vec bisection_Weibull(num_vec ci_to_fit)
{
    double mu = ci_to_fit[0];
    double ll = ci_to_fit[1];
    double ul = ci_to_fit[2];
    num_vec out;
    double a = 1e-2, b = 10000;
    double x = a;
    while ((b-a) >= ep)
    {
        x = (a+b)/2;
        boost::math::weibull_distribution<> Fx(x, mu/tgamma(1 + 1/x)); boost::math::weibull_distribution<> Fa(a, mu/tgamma(1 + 1.0/a));
        if ((cdf(Fx, ul) - cdf(Fx, ll) - 0.95)*(cdf(Fa, ul) - cdf(Fa, ll) - 0.95) < 0.0)
            b = x;
        else
            a = x;
    }
    boost::math::weibull_distribution<> Fx(x, mu/tgamma(1 + 1/x));
    double err = (abs(cdf(Fx, ll) - 0.025) + abs(cdf(Fx, ul) - 0.975));
    
    out.push_back(3.0);  out.push_back(x); out.push_back(mu/tgamma(1 + 1/x)); out.push_back(err);
    return out;
}

string fit_CI_to_dist(string name_dist, num_vec data_to_fit)
{
    num_vec errlist;
    string out;
    num_vec c1 = bisection_Gamma(data_to_fit);
    num_vec c2 = bisection_LN(data_to_fit);
    num_vec c3 = bisection_Weibull(data_to_fit);
    errlist.push_back(c1[3]); errlist.push_back(c2[3]); errlist.push_back(c3[3]);
    
    double minElementIndex = std::min_element(errlist.begin(),errlist.end()) - errlist.begin();
    if (minElementIndex == 0)
    {
        string par1 = to_string(c1[1]);
        string par2 = to_string(c1[2]);
        out = "boost::random::gamma_distribution<> " + name_dist + "(" + par1 + "," + par2 + ");";
    }
    else if (minElementIndex == 1)
    {
        string par1 = to_string(c2[1]);
        string par2 = to_string(c2[2]);
        out = "boost::random::lognormal_distribution<> " + name_dist + "(" + par1 + "," + par2 + ");";
    }
    else if (minElementIndex == 2)
    {
        string par1 = to_string(c3[1]);
        string par2 = to_string(c3[2]);
        out = "boost::random::weibull_distribution<> " + name_dist + "(" + par1 + "," + par2 + ");";
    }
    else
        cout << "Error" << endl;
    return out ;
}

namespace nr_data
{
    vector2D H_inc = {{1683.4, 1626.4, 1817.35}, {3965.3, 3795.25, 4096.4}, {3024.8, 2919.35, 3202.45}, {2206.85, 2090.95, 2328.45}, {1912.35, 1800.25, 2065.3}, {1577.95, 1484.85, 1701.45}, {1291.05, 1205.55, 1426.9}, {1064.95, 998.45, 1210.3}, {865.45, 794.2, 998.45}, {732.45, 663.1, 852.15}, {551., 501.6, 658.35}, {468.35, 398.05, 572.85}, {2563, 2153, 3188}, {2563, 2153, 3188}, {2563, 2153, 3188}, {2563, 2153, 3188}, {135.853, 101.89, 169.816}, {135.853, 101.89, 169.816}, {271.706, 203.78, 339.633}, {271.706, 203.78, 339.633}, {271.706, 203.78, 339.633}, {2037.79, 1494.38, 2445.35}, {2037.79, 1494.38, 2445.35}, {4822.78, 3559.35, 6079.42}, {19561.7, 14476.3, 24647.1}};
    
    vector2D D_inc = {{1.40381, 1.34721, 1.46042}, {1.40381, 1.34721, 1.46042}, {1.40381, 1.34721, 1.46042}, {1.40381, 1.34721, 1.46042}, {1.40381, 1.34721, 1.46042}, {1.40381, 1.34721, 1.46042}, {0.362275, 0.362275, 0.384917}, {0.362275, 0.362275, 0.384917}, {0.362275, 0.362275,0.384917}, {0.362275, 0.362275, 0.384917}, {0.362275, 0.362275, 0.384917}, {0.362275, 0.362275, 0.384917}, {4.3473, 4.3473, 4.619}, {4.3473, 4.3473, 4.619}, {4.3473, 4.3473, 4.619}, {4.3473,  4.3473, 4.619}, {14.63, 1., 174.}, {14.63, 1., 43.5}, {15.3333, 1, 96}, {15.3333, 1, 96}, {15.3333, 1, 96.}, {291, 17, 564}, {291, 17, 564}, {291, 17, 564}, {2806, 848, 4873}};
    
    vector2D GP_inc = {{2490.64, 2479.32, 2496.3}, {2490.64, 2479.32, 2496.3}, {2490.64, 2479.32, 2496.3}, {2490.64, 2479.32, 2496.3}, {2490.64, 2479.32, 2496.3}, {2490.64, 2479.32, 2496.3}, {283.027, 254.724, 311.33}, {283.027, 254.724, 311.33}, {283.027, 254.724,  311.33}, {283.027, 254.724, 311.33}, {283.027, 254.724, 311.33}, {283.027, 254.724, 311.33}, {3396.33, 3056.69,  3735.96}, {3396.33, 3056.69, 3735.96}, {3396.33, 3056.69, 3735.96}, {3396.33, 3056.69, 3735.96}, {37835.1, 21464.8, 54205.3}, {37835.1, 21464.8, 54205.3}, {45986.2, 29412.2, 57737.5}, {45986.2, 29412.2, 57737.5}, {45986.2, 29412.2, 57737.5}, {90002.6, 63035.8, 104743.}, {90002.6, 63035.8, 104743.}, {118328., 85519.5, 138434.}, {169644., 121208., 196241.}};
    
    vector2D BD_inc = {{8317.91, 7969.52, 9071.57}, {10055.9, 9714.63, 10358.5}, {6895.91, 6722.9, 7203.22}, {4538.55, 4375.81, 4723.41}, {3799.11, 3539.99, 3991.08}, {2831.36, 2695.48, 3042.29}, {2227.01, 2091.13, 2424.51}, {2054.79, 1962.36, 2233.33}, {1621.87, 1538.92, 1805.94}, {1389.61, 1292.44, 1557.09}, {1191.32, 1125.75, 1326.41}, {825.55, 729.96, 956.69}, {5126, 4306, 6376}, {5126, 4306, 6376}, {5126, 4306, 6376}, {5126, 4306, 6376}, {271.706, 203.78, 339.633}, {271.706, 203.78, 339.633}, {815.118, 611.339, 1018.9}, {815.118, 611.339, 1018.9}, {815.118, 611.339, 1018.9}, {6113.38, 4483.15, 7336.06}, {6113.38, 4483.15, 7336.06}, {14468.3, 10678., 18238.3}, {58685.1, 43428.9, 73941.4}};
    
    string risklab = "R1_A";
    
    vector< string > get_labels(string prefix)
    {
        vector< string > str;
        for (int i = 0; i < 25; i++)
        {
            string pre = prefix + risklab + to_string(i);
            str.push_back(pre);
        }
        return str;
    }
    vector< string > labels_H = get_labels("H_");
    vector< string > labels_D = get_labels("D_");
    vector< string > labels_GP = get_labels("GP_");
    vector< string > labels_BD = get_labels("BD_");
}

namespace hr_data
{
    vector2D H_inc = {{88.6, 85.6, 95.65}, {208.7, 199.75, 215.6}, {159.2, 153.65, 168.55}, {116.15, 110.05, 122.55}, {100.65, 94.75, 108.7}, {83.05, 78.15, 89.55}, {67.95, 63.45, 75.1}, {56.05, 52.55, 63.7}, {45.55, 41.8, 52.55}, {38.55, 34.9, 44.85}, {29., 26.4, 34.65}, {24.65, 20.95, 30.15}};
    
    vector2D BD_inc = {{2211.09, 2118.48, 2411.43}, {2673.09, 2582.37, 2753.52}, {1833.09, 1787.1, 1914.78}, {1206.45, 1163.19, 1255.59}, {1009.89, 941.01, 1060.92}, {752.64, 716.52, 808.71}, {591.99, 555.87, 644.49}, {546.21, 521.64, 593.67}, {431.13, 409.08, 480.06}, {369.39, 343.56, 413.91}, {316.68, 299.25, 352.59}, {219.45, 194.04, 254.31}};
    
    string risklab = "R2_A";
    
    vector< string > get_labels(string prefix)
    {
        vector< string > str;
        for (int i = 0; i < 25; i++)
        {
            string pre = prefix + risklab + to_string(i);
            str.push_back(pre);
        }
        return str;
    }
    
    vector< string > labels_H = get_labels("H_");
    vector< string > labels_BD = get_labels("BD_");
}

namespace qaly_data
{
    vector2D S;
    vector2D H;
    vector2D D;

}

void output_fitted_dist()
{
    string outmcyr = dout + "soln/" + "risks";
    //string FileMC_chain_M_full = outmcyr + "MC_chain_M_full_" + F + ".txt";      //Output MC
    ofstream risks_out(outmcyr.c_str());
    num_vec inc;
    risks_out << "// Risks for hospital admission no-at-risk persons" << endl;
    cout << "Risks for hospital admission no-at-risk persons" << endl;
    for (int a = 0; a < NoAgeG; a++)
    {
        inc.clear();
        inc.push_back(nr_data::H_inc[a][0]/ann_inc[a][0]);
        inc.push_back(nr_data::H_inc[a][1]/ann_inc[a][0]);
        inc.push_back(nr_data::H_inc[a][2]/ann_inc[a][0]);
        risks_out << fit_CI_to_dist(nr_data::labels_H[a], inc) << endl;
        cout << inc[0] << endl;
    }
    risks_out << "// Risks for deaths no-at-risk persons" << endl;
    cout << "Risks for deaths no-at-risk persons" << endl;

    for (int a = 0; a < NoAgeG; a++)
    {
        inc.clear();
        inc.push_back(nr_data::D_inc[a][0]/ann_inc[a][0]);
        inc.push_back(nr_data::D_inc[a][1]/ann_inc[a][0]);
        inc.push_back(nr_data::D_inc[a][2]/ann_inc[a][0]);
        risks_out << fit_CI_to_dist(nr_data::labels_D[a], inc) << endl;
        cout << inc[0] << endl;

    }
    risks_out << "// Risks for GP consultations no-at-risk persons" << endl;
    cout << "Risks for GP consultations no-at-risk persons" << endl;

    for (int a = 0; a < NoAgeG; a++)
    {
        inc.clear();
        inc.push_back(nr_data::GP_inc[a][0]/ann_inc[a][0]);
        inc.push_back(nr_data::GP_inc[a][1]/ann_inc[a][0]);
        inc.push_back(nr_data::GP_inc[a][2]/ann_inc[a][0]);
        risks_out << fit_CI_to_dist(nr_data::labels_GP[a], inc) << endl;
        cout << inc[0] << endl;

    }
    risks_out << "// Risks for number of bed days no-at-risk persons" << endl;
    cout << " Risks for number of bed days no-at-risk persons" << endl;

    for (int a = 0; a < NoAgeG; a++)
    {
        inc.clear();
        inc.push_back(nr_data::BD_inc[a][0]/ann_inc[a][0]);
        inc.push_back(nr_data::BD_inc[a][1]/ann_inc[a][0]);
        inc.push_back(nr_data::BD_inc[a][2]/ann_inc[a][0]);
        risks_out << fit_CI_to_dist(nr_data::labels_BD[a], inc) << endl;
        cout << inc[0] << endl;

    }
    
    risks_out << "// Risks for hospital admission high-risk persons" << endl;
    cout << " Risks for hospital admission high-risk persons" << endl;
    for (int a = 0; a < 12; a++)
    {
        inc.clear();
        inc.push_back(hr_data::H_inc[a][0]/ann_inc[a][1]);
        inc.push_back(hr_data::H_inc[a][1]/ann_inc[a][1]);
        inc.push_back(hr_data::H_inc[a][2]/ann_inc[a][1]);
        risks_out << fit_CI_to_dist(hr_data::labels_H[a], inc) << endl;
        cout << inc[0] << endl;

    }
    
    risks_out << "// Risks for number of bed days high-risk persons" << endl;
    cout << "// Risks for number of bed days high-risk persons" << endl;

    for (int a = 0; a < 12; a++)
    {
        inc.clear();
        inc.push_back(hr_data::BD_inc[a][0]/ann_inc[a][1]);
        inc.push_back(hr_data::BD_inc[a][1]/ann_inc[a][1]);
        inc.push_back(hr_data::BD_inc[a][2]/ann_inc[a][1]);
        risks_out << fit_CI_to_dist(hr_data::labels_BD[a], inc) << endl;
        cout << inc[0] << endl;

    }
    risks_out.close();

}
namespace out_dist
{
    // Risks for hospital admission no-at-risk persons
    boost::random::lognormal_distribution<> H_R1_A0(-3.115625,0.020926);
    boost::random::weibull_distribution<> H_R1_A1(62.227553,0.094955);
    boost::random::lognormal_distribution<> H_R1_A2(-2.639371,0.021130);
    boost::random::lognormal_distribution<> H_R1_A3(-2.956136,0.027445);
    boost::random::lognormal_distribution<> H_R1_A4(-3.103970,0.034059);
    boost::random::lognormal_distribution<> H_R1_A5(-3.300267,0.033973);
    boost::random::lognormal_distribution<> H_R1_A6(-3.501204,0.040090);
    boost::random::lognormal_distribution<> H_R1_A7(-3.683768,0.039078);
    boost::random::lognormal_distribution<> H_R1_A8(-3.874510,0.051394);
    boost::random::lognormal_distribution<> H_R1_A9(-4.018271,0.058688);
    boost::random::lognormal_distribution<> H_R1_A10(-4.274218,0.056815);
    boost::random::lognormal_distribution<> H_R1_A11(-4.405675,0.090839);
    boost::random::lognormal_distribution<> H_R1_A12(-5.332327,0.097721);
    boost::random::lognormal_distribution<> H_R1_A13(-5.102881,0.097721);
    boost::random::lognormal_distribution<> H_R1_A14(-5.000518,0.097721);
    boost::random::lognormal_distribution<> H_R1_A15(-4.935537,0.097721);
    boost::random::gamma_distribution<> H_R1_A16(60.924023,0.000001);
    boost::random::gamma_distribution<> H_R1_A17(60.924023,0.000001);
    boost::random::gamma_distribution<> H_R1_A18(60.923051,0.000001);
    boost::random::gamma_distribution<> H_R1_A19(60.923051,0.000002);
    boost::random::gamma_distribution<> H_R1_A20(60.923051,0.000001);
    boost::random::weibull_distribution<> H_R1_A21(10.009425,0.000849);
    boost::random::weibull_distribution<> H_R1_A22(10.009425,0.001267);
    boost::random::gamma_distribution<> H_R1_A23(55.761551,0.000075);
    boost::random::gamma_distribution<> H_R1_A24(56.300550,0.000270);
    // Risks for deaths no-at-risk persons
    boost::random::gamma_distribution<> D_R1_A0(2362.125146,0.00000001);
    boost::random::gamma_distribution<> D_R1_A1(2362.125146,0.00000001);
    boost::random::gamma_distribution<> D_R1_A2(2362.125146,0.00000001);
    boost::random::gamma_distribution<> D_R1_A3(2362.125146,0.00000001);
    boost::random::gamma_distribution<> D_R1_A4(2362.125146,0.00000001);
    boost::random::gamma_distribution<> D_R1_A5(2362.125146,0.00000001);
    boost::random::weibull_distribution<> D_R1_A6(9999.999999,0.000008);
    boost::random::weibull_distribution<> D_R1_A7(9999.999999,0.000009);
    boost::random::weibull_distribution<> D_R1_A8(9999.999999,0.000009);
    boost::random::weibull_distribution<> D_R1_A9(9999.999999,0.000009);
    boost::random::weibull_distribution<> D_R1_A10(9999.999999,0.000009);
    boost::random::weibull_distribution<> D_R1_A11(9999.999999,0.000009);
    boost::random::weibull_distribution<> D_R1_A12(9999.999999,0.000008);
    boost::random::weibull_distribution<> D_R1_A13(9999.999999,0.000010);
    boost::random::weibull_distribution<> D_R1_A14(9999.999999,0.000011);
    boost::random::weibull_distribution<> D_R1_A15(9999.999999,0.000012);
    boost::random::lognormal_distribution<> D_R1_A16(-11.908821,1.312034);
    boost::random::weibull_distribution<> D_R1_A17(1.309704,0.000007);
    boost::random::lognormal_distribution<> D_R1_A18(-12.275800,1.077579);
    boost::random::lognormal_distribution<> D_R1_A19(-12.169540,1.077579);
    boost::random::lognormal_distribution<> D_R1_A20(-12.245614,1.077579);
    boost::random::weibull_distribution<> D_R1_A21(2.058097,0.000130);
    boost::random::weibull_distribution<> D_R1_A22(2.058097,0.000194);
    boost::random::weibull_distribution<> D_R1_A23(2.058097,0.000286);
    boost::random::weibull_distribution<> D_R1_A24(2.908690,0.002443);
    // Risks for GP consultations no-at-risk persons
    boost::random::weibull_distribution<> GP_R1_A0(788.190729,0.065666);
    boost::random::weibull_distribution<> GP_R1_A1(788.190729,0.059147);
    boost::random::weibull_distribution<> GP_R1_A2(788.190729,0.058839);
    boost::random::weibull_distribution<> GP_R1_A3(788.190729,0.058752);
    boost::random::weibull_distribution<> GP_R1_A4(788.190729,0.058482);
    boost::random::weibull_distribution<> GP_R1_A5(788.190729,0.058244);
    boost::random::gamma_distribution<> GP_R1_A6(383.596616,0.000017);
    boost::random::gamma_distribution<> GP_R1_A7(383.596616,0.000017);
    boost::random::gamma_distribution<> GP_R1_A8(383.596616,0.000018);
    boost::random::gamma_distribution<> GP_R1_A9(383.596616,0.000018);
    boost::random::gamma_distribution<> GP_R1_A10(383.596616,0.000019);
    boost::random::gamma_distribution<> GP_R1_A11(383.596616,0.000019);
    boost::random::gamma_distribution<> GP_R1_A12(383.600556,0.000017);
    boost::random::gamma_distribution<> GP_R1_A13(383.600556,0.000021);
    boost::random::gamma_distribution<> GP_R1_A14(383.600556,0.000023);
    boost::random::gamma_distribution<> GP_R1_A15(383.600556,0.000025);
    boost::random::weibull_distribution<> GP_R1_A16(5.159418,0.018924);
    boost::random::weibull_distribution<> GP_R1_A17(5.159418,0.019388);
    boost::random::weibull_distribution<> GP_R1_A18(7.620541,0.014886);
    boost::random::weibull_distribution<> GP_R1_A19(7.620541,0.016555);
    boost::random::weibull_distribution<> GP_R1_A20(7.620541,0.015342);
    boost::random::weibull_distribution<> GP_R1_A21(11.119634,0.037348);
    boost::random::weibull_distribution<> GP_R1_A22(11.119634,0.055732);
    boost::random::weibull_distribution<> GP_R1_A23(11.023472,0.107717);
    boost::random::weibull_distribution<> GP_R1_A24(11.630569,0.137592);
    // Risks for number of bed days no-at-risk persons
    boost::random::lognormal_distribution<> BD_R1_A0(-1.518030,0.025948);
    boost::random::gamma_distribution<> BD_R1_A1(3804.312368,0.000063);
    boost::random::lognormal_distribution<> BD_R1_A2(-1.815288,0.015253);
    boost::random::lognormal_distribution<> BD_R1_A3(-2.235095,0.019425);
    boost::random::weibull_distribution<> BD_R1_A4(39.946103,0.090392);
    boost::random::lognormal_distribution<> BD_R1_A5(-2.715636,0.028779);
    boost::random::lognormal_distribution<> BD_R1_A6(-2.956000,0.036121);
    boost::random::lognormal_distribution<> BD_R1_A7(-3.026522,0.027760);
    boost::random::lognormal_distribution<> BD_R1_A8(-3.246425,0.031848);
    boost::random::lognormal_distribution<> BD_R1_A9(-3.377888,0.043004);
    boost::random::lognormal_distribution<> BD_R1_A10(-3.503136,0.034244);
    boost::random::lognormal_distribution<> BD_R1_A11(-3.838841,0.067924);
    boost::random::lognormal_distribution<> BD_R1_A12(-4.639179,0.097721);
    boost::random::lognormal_distribution<> BD_R1_A13(-4.409734,0.097721);
    boost::random::lognormal_distribution<> BD_R1_A14(-4.307371,0.097721);
    boost::random::lognormal_distribution<> BD_R1_A15(-4.242390,0.097721);
    boost::random::gamma_distribution<> BD_R1_A16(60.923051,0.000002);
    boost::random::gamma_distribution<> BD_R1_A17(60.923051,0.000002);
    boost::random::gamma_distribution<> BD_R1_A18(60.922448,0.000004);
    boost::random::gamma_distribution<> BD_R1_A19(60.922448,0.000005);
    boost::random::gamma_distribution<> BD_R1_A20(60.922448,0.000004);
    boost::random::weibull_distribution<> BD_R1_A21(10.009442,0.002548);
    boost::random::weibull_distribution<> BD_R1_A22(10.009442,0.003802);
    boost::random::gamma_distribution<> BD_R1_A23(55.759807,0.000226);
    boost::random::gamma_distribution<> BD_R1_A24(56.300150,0.000809);
    // Risks for hospital admission high-risk persons
    boost::random::lognormal_distribution<> H_R2_A0(-2.828635,0.020926);
    boost::random::weibull_distribution<> H_R2_A1(62.227553,0.126519);
    boost::random::lognormal_distribution<> H_R2_A2(-2.352380,0.021130);
    boost::random::lognormal_distribution<> H_R2_A3(-2.669146,0.027445);
    boost::random::lognormal_distribution<> H_R2_A4(-2.816981,0.034059);
    boost::random::lognormal_distribution<> H_R2_A5(-3.013278,0.033973);
    boost::random::lognormal_distribution<> H_R2_A6(-3.214214,0.040090);
    boost::random::lognormal_distribution<> H_R2_A7(-3.396781,0.039078);
    boost::random::lognormal_distribution<> H_R2_A8(-3.587518,0.051394);
    boost::random::lognormal_distribution<> H_R2_A9(-3.731282,0.058688);
    boost::random::lognormal_distribution<> H_R2_A10(-3.987226,0.056815);
    boost::random::lognormal_distribution<> H_R2_A11(-4.118684,0.090839);
    // Risks for number of bed days high-risk persons
    boost::random::lognormal_distribution<> BD_R2_A0(0.388474,0.025948);
    boost::random::gamma_distribution<> BD_R2_A1(3804.674598,0.000422);
    boost::random::lognormal_distribution<> BD_R2_A2(0.091217,0.015253);
    boost::random::lognormal_distribution<> BD_R2_A3(-0.328591,0.019425);
    boost::random::weibull_distribution<> BD_R2_A4(39.946103,0.608296);
    boost::random::lognormal_distribution<> BD_R2_A5(-0.809133,0.028779);
    boost::random::lognormal_distribution<> BD_R2_A6(-1.049497,0.036121);
    boost::random::lognormal_distribution<> BD_R2_A7(-1.120021,0.027760);
    boost::random::lognormal_distribution<> BD_R2_A8(-1.339918,0.031848);
    boost::random::lognormal_distribution<> BD_R2_A9(-1.471385,0.043004);
    boost::random::lognormal_distribution<> BD_R2_A10(-1.596630,0.034244);
    boost::random::lognormal_distribution<> BD_R2_A11(-1.932337,0.067924);


    // Distributions for H (for very-high risk groups)
    boost::random::gamma_distribution<> H_R3_A0(1707.23, 0.000218917); //0-2 months, VH risk
    boost::random::gamma_distribution<> H_R3_A3(1097.73, 0.000218973); //3-5 months, VH risk
    boost::random::gamma_distribution<> H_R3_A6(905.012, 0.000156833); //6-8 months, VH risk
    
    // Distributions for BD (for very-high risk groups)
    boost::random::gamma_distribution<> BD_R3_A0(91.7501,0.0254); //0-2 months, VH risk
    boost::random::gamma_distribution<> BD_R3_A3(87.2681,0.01718); //3-5 months, VH risk
    boost::random::gamma_distribution<> BD_R3_A6(87.4392,0.01014); //6-8 months, VH risk
    
    // Distributions for D (for very-high risk groups)
    boost::random::gamma_distribution<> D_R3_A0(8.31507, 0.00166481); //0-2 months, VH risk
    boost::random::gamma_distribution<> D_R3_A3(8.18185, 0.00108933); //3-5 months, VH risk
    boost::random::gamma_distribution<> D_R3_A6(8.27439, 0.00063814); //6-8 months, VH risk

    // Distributions for the QALY loss due to each infection
    boost::random::gamma_distribution<> Q_S1(1.57764, 0.0014807); // <5years, symptomatic
    boost::random::gamma_distribution<> Q_S2(1.41075, 0.0010264); // >5 years, symptomatic
    boost::random::gamma_distribution<> Q_H1(2.0017, 0.00204726); // <5years, hospitalised
    boost::random::gamma_distribution<> Q_H2(1.60289, 0.00186539); // >5 years, hospitalised

    boost::random::gamma_distribution<> Q_NHS1(1.6578, 0.0018241); // <5years, symptomatic
    boost::random::gamma_distribution<> Q_NHS2(1.36973, 0.0011265); // >5 years, symptomatic
    boost::random::gamma_distribution<> Q_HS1(1.7927, 0.00213254); // <5years, hospitalised
    boost::random::lognormal_distribution<> Q_HS2(-6.23993, 0.933905); // >5 years, hospitalised
    
         /// https://onlinelibrary.wiley.com/doi/full/10.1002/hec.1823 for the Quality adjusted life year per age
    boost::random::normal_distribution<> Q_D1(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D2(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D3(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D4(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D5(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D6(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D7(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D8(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D9(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D10(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D11(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D12(23.32, 23.32*0.1);
    boost::random::normal_distribution<> Q_D13(23.2186, 23.2186*0.1);
    boost::random::normal_distribution<> Q_D14(23.1135, 23.1135*0.1);
    boost::random::normal_distribution<> Q_D15(23.0048, 23.0048*0.1);
    boost::random::normal_distribution<> Q_D16(22.8921, 22.8921*0.1);
    boost::random::normal_distribution<> Q_D17(22.5251, 22.5251*0.1);
    boost::random::normal_distribution<> Q_D18(21.8287, 21.8287*0.1);
    boost::random::normal_distribution<> Q_D19(20.7512, 20.7512*0.1);
    boost::random::normal_distribution<> Q_D20(19.3996, 19.3996*0.1);
    boost::random::normal_distribution<> Q_D21(17.616, 17.616*0.1);
    boost::random::normal_distribution<> Q_D22(15.2231, 15.2231*0.1);
    boost::random::normal_distribution<> Q_D23(12.038, 12.038*0.1);
    boost::random::normal_distribution<> Q_D24(7.63045, 7.63045*0.1);
    boost::random::normal_distribution<> Q_D25(3.04609, 3.04609*0.1);

    vector< boost::random::normal_distribution<double>> Q_D = {Q_D1, Q_D2, Q_D3, Q_D4, Q_D5, Q_D6, Q_D7, Q_D8, Q_D9, Q_D10, Q_D11, Q_D12, Q_D13, Q_D14, Q_D15, Q_D16, Q_D17, Q_D18, Q_D19, Q_D20, Q_D21, Q_D22, Q_D23, Q_D24, Q_D25};
    
    //Costs
    boost::random::normal_distribution<> C_H1(725.293, 4.12643);
    boost::random::normal_distribution<> C_H2(425.242, 5.27808);
}

using namespace out_dist;


num_vec get_GP(num_vec annual_inc, int s)
{

    PRNG_s rng(s);
    num_vec gp_tot(NoAgeG,0);
    // Covert to risk and age only
    num_vec vhr(NoAgeG,0);
    num_vec hr(NoAgeG,0);
    num_vec nr(NoAgeG,0);
    for (int a = 0; a < NoAgeG; a++)
        for (int s = 0; s < 3; s++){
            vhr[a] += annual_inc[9*a+3*s+0]; hr[a] += annual_inc[9*a+3*s+1]; nr[a] += annual_inc[9*a+3*s+2];
        }
    
    // Calculate the number of GP consultations
    gp_tot[0] = vhr[0]*GP_R1_A0(rng) + hr[0]*GP_R1_A0(rng) + nr[0]*GP_R1_A0(rng);      // 0 months
    gp_tot[1] = vhr[1]*GP_R1_A1(rng) + hr[1]*GP_R1_A1(rng) + nr[1]*GP_R1_A1(rng);      // 1 months
    gp_tot[2] = vhr[2]*GP_R1_A2(rng) + hr[2]*GP_R1_A2(rng) + nr[2]*GP_R1_A2(rng);      // 2 months
    gp_tot[3] = vhr[3]*GP_R1_A3(rng) + hr[3]*GP_R1_A3(rng) + nr[3]*GP_R1_A3(rng);      // 3 months
    gp_tot[4] = vhr[4]*GP_R1_A4(rng) + hr[4]*GP_R1_A4(rng) + nr[4]*GP_R1_A4(rng);      // 4 months
    gp_tot[5] = vhr[5]*GP_R1_A5(rng) + hr[5]*GP_R1_A5(rng) + nr[5]*GP_R1_A5(rng);      // 5 months
    gp_tot[6] = vhr[6]*GP_R1_A6(rng) + hr[6]*GP_R1_A6(rng) + nr[6]*GP_R1_A6(rng);      // 6 months
    gp_tot[7] = vhr[7]*GP_R1_A7(rng) + hr[7]*GP_R1_A7(rng) + nr[7]*GP_R1_A7(rng);      // 7 months
    gp_tot[8] = vhr[8]*GP_R1_A8(rng) + hr[8]*GP_R1_A8(rng) + nr[8]*GP_R1_A8(rng);      // 8 months
    gp_tot[9] = vhr[9]*0             + hr[9]*GP_R1_A9(rng) + nr[9]*GP_R1_A9(rng);     // 9 months
    gp_tot[10] = vhr[10]*0           + hr[10]*GP_R1_A10(rng) + nr[10]*GP_R1_A10(rng);      // 10 months
    gp_tot[11] = vhr[11]*0           + hr[11]*GP_R1_A11(rng) + nr[11]*GP_R1_A11(rng);      // 11 months
    gp_tot[12] = vhr[12]*0           + hr[12]*0              + nr[12]*GP_R1_A12(rng);      // 1 yr
    gp_tot[13] = vhr[13]*0           + hr[13]*0              + nr[13]*GP_R1_A13(rng);      // 2 yr
    gp_tot[14] = vhr[14]*0           + hr[14]*0              + nr[14]*GP_R1_A14(rng);      // 3 yr
    gp_tot[15] = vhr[15]*0           + hr[15]*0              + nr[15]*GP_R1_A15(rng);      // 4 yr
    gp_tot[16] = vhr[16]*0           + hr[16]*0              + nr[16]*GP_R1_A16(rng);      // 5-9 yr
    gp_tot[17] = vhr[17]*0           + hr[17]*0              + nr[17]*GP_R1_A17(rng);      // 10-14 yr
    gp_tot[18] = vhr[18]*0           + hr[18]*0              + nr[18]*GP_R1_A18(rng);      // 15-24 yr
    gp_tot[19] = vhr[19]*0           + hr[19]*0              + nr[19]*GP_R1_A19(rng);      // 25-34 yr
    gp_tot[20] = vhr[20]*0           + hr[20]*0              + nr[20]*GP_R1_A20(rng);      // 35-44 yr
    gp_tot[21] = vhr[21]*0           + hr[21]*0              + nr[21]*GP_R1_A21(rng);      // 45-54 yr
    gp_tot[22] = vhr[22]*0           + hr[22]*0              + nr[22]*GP_R1_A22(rng);      // 55-64 yr
    gp_tot[23] = vhr[23]*0           + hr[23]*0              + nr[23]*GP_R1_A23(rng);      // 65-74 yr
    gp_tot[24] = vhr[24]*0           + hr[24]*0              + nr[24]*GP_R1_A24(rng);      // 75+ yr

    return gp_tot;
}

num_vec get_BD(num_vec annual_inc, int s)
{

    PRNG_s rng(s);

    num_vec bd_tot(NoAgeG,0);
    // Covert to risk and age only
    num_vec vhr(NoAgeG,0);
    num_vec hr(NoAgeG,0);
    num_vec nr(NoAgeG,0);
    for (int a = 0; a < NoAgeG; a++)
        for (int s = 0; s < 3; s++){
            vhr[a] += annual_inc[9*a+3*s+0]; hr[a] += annual_inc[9*a+3*s+1]; nr[a] += annual_inc[9*a+3*s+2];
        }


    // Calculate the number of GP consultations
    bd_tot[0] = vhr[0]*BD_R3_A0(rng)   + hr[0]*BD_R2_A0(rng)      + nr[0]*BD_R1_A0(rng);      // 0 months
    bd_tot[1] = vhr[1]*BD_R3_A0(rng)   + hr[1]*BD_R2_A1(rng)      + nr[1]*BD_R1_A1(rng);      // 1 months
    bd_tot[2] = vhr[2]*BD_R3_A0(rng)   + hr[2]*BD_R2_A2(rng)      + nr[2]*BD_R1_A2(rng);      // 2 months
    bd_tot[3] = vhr[3]*BD_R3_A3(rng)   + hr[3]*BD_R2_A3(rng)      + nr[3]*BD_R1_A3(rng);      // 3 months
    bd_tot[4] = vhr[4]*BD_R3_A3(rng)   + hr[4]*BD_R2_A4(rng)       + nr[4]*BD_R1_A4(rng);      // 4 months
    bd_tot[5] = vhr[5]*BD_R3_A3(rng)   + hr[5]*BD_R2_A5(rng)      + nr[5]*BD_R1_A5(rng);      // 5 months
    bd_tot[6] = vhr[6]*BD_R3_A6(rng)   + hr[6]*BD_R2_A6(rng)      + nr[6]*BD_R1_A6(rng);      // 6 months
    bd_tot[7] = vhr[7]*BD_R3_A6(rng)   + hr[7]*BD_R2_A7(rng)      + nr[7]*BD_R1_A7(rng);      // 7 months
    bd_tot[8] = vhr[8]*BD_R3_A6(rng)   + hr[8]*BD_R2_A8(rng)      + nr[8]*BD_R1_A8(rng);      // 8 months
    bd_tot[9] = vhr[9]*0                + hr[9]*BD_R2_A9(rng)     + nr[9]*BD_R1_A9(rng);     // 9 months
    bd_tot[10] = vhr[10]*0              + hr[10]*BD_R2_A10(rng)    + nr[10]*BD_R1_A10(rng);      // 10 months
    bd_tot[11] = vhr[11]*0              + hr[11]*BD_R2_A11(rng)    + nr[11]*BD_R1_A11(rng);      // 11 months
    bd_tot[12] = vhr[12]*0              + hr[12]*0                  + nr[12]*BD_R1_A12(rng);      // 1 yr
    bd_tot[13] = vhr[13]*0              + hr[13]*0                  + nr[13]*BD_R1_A13(rng);      // 2 yr
    bd_tot[14] = vhr[14]*0              + hr[14]*0                  + nr[14]*BD_R1_A14(rng);      // 3 yr
    bd_tot[15] = vhr[15]*0              + hr[15]*0                  + nr[15]*BD_R1_A15(rng);      // 4 yr
    bd_tot[16] = vhr[16]*0              + hr[16]*0                  + nr[16]*BD_R1_A16(rng);      // 5-9 yr
    bd_tot[17] = vhr[17]*0              + hr[17]*0                  + nr[17]*BD_R1_A17(rng);      // 10-14 yr
    bd_tot[18] = vhr[18]*0              + hr[18]*0                  + nr[18]*BD_R1_A18(rng);      // 15-24 yr
    bd_tot[19] = vhr[19]*0              + hr[19]*0                  + nr[19]*BD_R1_A19(rng);      // 25-34 yr
    bd_tot[20] = vhr[20]*0              + hr[20]*0                  + nr[20]*BD_R1_A20(rng);      // 35-44 yr
    bd_tot[21] = vhr[21]*0              + hr[21]*0                  + nr[21]*BD_R1_A21(rng);      // 45-54 yr
    bd_tot[22] = vhr[22]*0              + hr[22]*0                  + nr[22]*BD_R1_A22(rng);      // 55-64 yr
    bd_tot[23] = vhr[23]*0              + hr[23]*0                  + nr[23]*BD_R1_A23(rng);      // 65-74 yr
    bd_tot[24] = vhr[24]*0              + hr[24]*0                  + nr[24]*BD_R1_A24(rng);      // 75+ yr
    
    return bd_tot;
}

num_vec get_S(num_vec annual_inc, amh::amh_state_t& mcmc_state, int s1)
{
    int s;

    s = s1;
    
    num_vec s_tot(NoAgeG,0);
    // Covert to risk and age only
    num_vec vhr(NoAgeG,0);
    num_vec hr(NoAgeG,0);
    num_vec nr(NoAgeG,0);
    for (int a = 0; a < NoAgeG; a++)
        for (int p = 0; p < 3; p++){
            vhr[a] += annual_inc[9*a+3*p+0]; hr[a] += annual_inc[9*a+3*p+1]; nr[a] += annual_inc[9*a+3*p+2];
        }
    double pA1, pA2, pA3, pA4;
    

    pA1 = 1 - mcmc_state.posterior(s,6);
    pA2 = 1 - mcmc_state.posterior(s,7);
    pA3 = 1 - mcmc_state.posterior(s,8);
    pA4 = 1 - mcmc_state.posterior(s,9);

    // Calculate the number of GP consultations
    s_tot[0] = vhr[0]*pA1   + hr[0]*pA1      + nr[0]*pA1;     // 0 months
    s_tot[1] = vhr[1]*pA1   + hr[1]*pA1      + nr[1]*pA1;       // 1 months
    s_tot[2] = vhr[2]*pA1   + hr[2]*pA1     + nr[2]*pA1;       // 2 months
    s_tot[3] = vhr[3]*pA1   + hr[3]*pA1     + nr[3]*pA1;       // 3 months
    s_tot[4] = vhr[4]*pA1   + hr[4]*pA1      + nr[4]*pA1;       // 4 months
    s_tot[5] = vhr[5]*pA1   + hr[5]*pA1     + nr[5]*pA1;       // 5 months
    s_tot[6] = vhr[6]*pA1   + hr[6]*pA1    + nr[6]*pA1;       // 6 months
    s_tot[7] = vhr[7]*pA1   + hr[7]*pA1    + nr[7]*pA1;       // 7 months
    s_tot[8] = vhr[8]*pA1   + hr[8]*pA1     + nr[8]*pA1;      // 8 months
    s_tot[9] = vhr[9]*0     + hr[9]*pA1     + nr[9]*pA1;      // 9 months
    s_tot[10] = vhr[10]*0   + hr[10]*pA1    + nr[10]*pA1;       // 10 months
    s_tot[11] = vhr[11]*0   + hr[11]*pA1   + nr[11]*pA1;      // 11 months
    s_tot[12] = vhr[12]*0   + hr[12]*0              + nr[12]*pA2;     // 1 yr
    s_tot[13] = vhr[13]*0   + hr[13]*0              + nr[13]*pA2;        // 2 yr
    s_tot[14] = vhr[14]*0              + hr[14]*0              + nr[14]*pA2;       // 3 yr
    s_tot[15] = vhr[15]*0              + hr[15]*0              + nr[15]*pA2;        // 4 yr
    s_tot[16] = vhr[16]*0              + hr[16]*0              + nr[16]*pA3;       // 5-9 yr
    s_tot[17] = vhr[17]*0              + hr[17]*0              + nr[17]*pA3;        // 10-14 yr
    s_tot[18] = vhr[18]*0              + hr[18]*0              + nr[18]*pA4;         // 15-24 yr
    s_tot[19] = vhr[19]*0              + hr[19]*0              + nr[19]*pA4;     // 25-34 yr
    s_tot[20] = vhr[20]*0              + hr[20]*0              + nr[20]*pA4;       // 35-44 yr
    s_tot[21] = vhr[21]*0              + hr[21]*0              + nr[21]*pA4;       // 45-54 yr
    s_tot[22] = vhr[22]*0              + hr[22]*0              + nr[22]*pA4;       // 55-64 yr
    s_tot[23] = vhr[23]*0              + hr[23]*0              + nr[23]*pA4;      // 65-74 yr
    s_tot[24] = vhr[24]*0              + hr[24]*0              + nr[24]*pA4;     // 75+ yr
    
    return s_tot;
}

num_vec get_H(num_vec annual_inc, int s)
{

    
    PRNG_s rng(s);

    num_vec h_tot(NoAgeG,0);
    // Covert to risk and age only
    num_vec vhr(NoAgeG,0);
    num_vec hr(NoAgeG,0);
    num_vec nr(NoAgeG,0);
    for (int a = 0; a < NoAgeG; a++)
        for (int s = 0; s < 3; s++){
            vhr[a] += annual_inc[9*a+3*s+0]; hr[a] += annual_inc[9*a+3*s+1]; nr[a] += annual_inc[9*a+3*s+2];
        }

    // Calculate the number of GP consultations
    h_tot[0] = vhr[0]*H_R3_A0(rng)   + hr[0]*H_R2_A0(rng)      + nr[0]*H_R1_A0(rng);      // 0 months
    h_tot[1] = vhr[1]*H_R3_A0(rng)   + hr[1]*H_R2_A1(rng)      + nr[1]*H_R1_A1(rng);      // 1 months
    h_tot[2] = vhr[2]*H_R3_A0(rng)   + hr[2]*H_R2_A2(rng)      + nr[2]*H_R1_A2(rng);      // 2 months
    h_tot[3] = vhr[3]*H_R3_A3(rng)   + hr[3]*H_R2_A3(rng)      + nr[3]*H_R1_A3(rng);      // 3 months
    h_tot[4] = vhr[4]*H_R3_A3(rng)   + hr[4]*H_R2_A4(rng)       + nr[4]*H_R1_A4(rng);      // 4 months
    h_tot[5] = vhr[5]*H_R3_A3(rng)   + hr[5]*H_R2_A5(rng)      + nr[5]*H_R1_A5(rng);      // 5 months
    h_tot[6] = vhr[6]*H_R3_A6(rng)   + hr[6]*H_R2_A6(rng)      + nr[6]*H_R1_A6(rng);      // 6 months
    h_tot[7] = vhr[7]*H_R3_A6(rng)   + hr[7]*H_R2_A7(rng)      + nr[7]*H_R1_A7(rng);      // 7 months
    h_tot[8] = vhr[8]*H_R3_A6(rng)   + hr[8]*H_R2_A8(rng)      + nr[8]*H_R1_A8(rng);      // 8 months
    h_tot[9] = vhr[9]*0                + hr[9]*H_R2_A9(rng)     + nr[9]*H_R1_A9(rng);     // 9 months
    h_tot[10] = vhr[10]*0              + hr[10]*H_R2_A10(rng)    + nr[10]*H_R1_A10(rng);      // 10 months
    h_tot[11] = vhr[11]*0              + hr[11]*H_R2_A11(rng)    + nr[11]*H_R1_A11(rng);      // 11 months
    h_tot[12] = vhr[12]*0              + hr[12]*0                  + nr[12]*H_R1_A12(rng);      // 1 yr
    h_tot[13] = vhr[13]*0              + hr[13]*0                  + nr[13]*H_R1_A13(rng);      // 2 yr
    h_tot[14] = vhr[14]*0              + hr[14]*0                  + nr[14]*H_R1_A14(rng);      // 3 yr
    h_tot[15] = vhr[15]*0              + hr[15]*0                  + nr[15]*H_R1_A15(rng);      // 4 yr
    h_tot[16] = vhr[16]*0              + hr[16]*0                  + nr[16]*H_R1_A16(rng);      // 5-9 yr
    h_tot[17] = vhr[17]*0              + hr[17]*0                  + nr[17]*H_R1_A17(rng);      // 10-14 yr
    h_tot[18] = vhr[18]*0              + hr[18]*0                  + nr[18]*H_R1_A18(rng);      // 15-24 yr
    h_tot[19] = vhr[19]*0              + hr[19]*0                  + nr[19]*H_R1_A19(rng);      // 25-34 yr
    h_tot[20] = vhr[20]*0              + hr[20]*0                  + nr[20]*H_R1_A20(rng);      // 35-44 yr
    h_tot[21] = vhr[21]*0              + hr[21]*0                  + nr[21]*H_R1_A21(rng);      // 45-54 yr
    h_tot[22] = vhr[22]*0              + hr[22]*0                  + nr[22]*H_R1_A22(rng);      // 55-64 yr
    h_tot[23] = vhr[23]*0              + hr[23]*0                  + nr[23]*H_R1_A23(rng);      // 65-74 yr
    h_tot[24] = vhr[24]*0              + hr[24]*0                  + nr[24]*H_R1_A24(rng);      // 75+ yr
    
    return h_tot;
}

num_vec get_D(num_vec annual_inc, int s)
{

    PRNG_s rng(s);

    num_vec d_tot(NoAgeG,0);
    // Covert to risk and age only
    num_vec vhr(NoAgeG,0);
    num_vec hr(NoAgeG,0);
    num_vec nr(NoAgeG,0);
    for (int a = 0; a < NoAgeG; a++)
        for (int s = 0; s < 3; s++){
            vhr[a] += annual_inc[9*a+3*s+0]; hr[a] += annual_inc[9*a+3*s+1]; nr[a] += annual_inc[9*a+3*s+2];
        }
    
    // Calculate the number of GP consultations
    d_tot[0] = vhr[0]*D_R3_A0(rng)   + hr[0]*D_R1_A0(rng)      + nr[0]*D_R1_A0(rng);      // 0 months
    d_tot[1] = vhr[1]*D_R3_A0(rng)   + hr[1]*D_R1_A1(rng)      + nr[1]*D_R1_A1(rng);      // 1 months
    d_tot[2] = vhr[2]*D_R3_A0(rng)   + hr[2]*D_R1_A2(rng)      + nr[2]*D_R1_A2(rng);      // 2 months
    d_tot[3] = vhr[3]*D_R3_A3(rng)   + hr[3]*D_R1_A3(rng)      + nr[3]*D_R1_A3(rng);      // 3 months
    d_tot[4] = vhr[4]*D_R3_A3(rng)   + hr[4]*D_R1_A4(rng)       + nr[4]*D_R1_A4(rng);      // 4 months
    d_tot[5] = vhr[5]*D_R3_A3(rng)   + hr[5]*D_R1_A5(rng)      + nr[5]*D_R1_A5(rng);      // 5 months
    d_tot[6] = vhr[6]*D_R3_A6(rng)   + hr[6]*D_R1_A6(rng)      + nr[6]*D_R1_A6(rng);      // 6 months
    d_tot[7] = vhr[7]*D_R3_A6(rng)   + hr[7]*D_R1_A7(rng)      + nr[7]*D_R1_A7(rng);      // 7 months
    d_tot[8] = vhr[8]*D_R3_A6(rng)   + hr[8]*D_R1_A8(rng)      + nr[8]*D_R1_A8(rng);      // 8 months
    d_tot[9] = vhr[9]*0                + hr[9]*D_R1_A9(rng)     + nr[9]*D_R1_A9(rng);     // 9 months
    d_tot[10] = vhr[10]*0              + hr[10]*D_R1_A10(rng)    + nr[10]*D_R1_A10(rng);      // 10 months
    d_tot[11] = vhr[11]*0              + hr[11]*D_R1_A11(rng)    + nr[11]*D_R1_A11(rng);      // 11 months
    d_tot[12] = vhr[12]*0              + hr[12]*0                  + nr[12]*D_R1_A12(rng);      // 1 yr
    d_tot[13] = vhr[13]*0              + hr[13]*0                  + nr[13]*D_R1_A13(rng);      // 2 yr
    d_tot[14] = vhr[14]*0              + hr[14]*0                  + nr[14]*D_R1_A14(rng);      // 3 yr
    d_tot[15] = vhr[15]*0              + hr[15]*0                  + nr[15]*D_R1_A15(rng);      // 4 yr
    d_tot[16] = vhr[16]*0              + hr[16]*0                  + nr[16]*D_R1_A16(rng);      // 5-9 yr
    d_tot[17] = vhr[17]*0              + hr[17]*0                  + nr[17]*D_R1_A17(rng);      // 10-14 yr
    d_tot[18] = vhr[18]*0              + hr[18]*0                  + nr[18]*D_R1_A18(rng);      // 15-24 yr
    d_tot[19] = vhr[19]*0              + hr[19]*0                  + nr[19]*D_R1_A19(rng);      // 25-34 yr
    d_tot[20] = vhr[20]*0              + hr[20]*0                  + nr[20]*D_R1_A20(rng);      // 35-44 yr
    d_tot[21] = vhr[21]*0              + hr[21]*0                  + nr[21]*D_R1_A21(rng);      // 45-54 yr
    d_tot[22] = vhr[22]*0              + hr[22]*0                  + nr[22]*D_R1_A22(rng);      // 55-64 yr
    d_tot[23] = vhr[23]*0              + hr[23]*0                  + nr[23]*D_R1_A23(rng);      // 65-74 yr
    d_tot[24] = vhr[24]*0              + hr[24]*0                  + nr[24]*D_R1_A24(rng);      // 75+ yr
    
    return d_tot;
}


#endif /* outcomes_h */
