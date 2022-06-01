/*!
    Program for Exercise set 12 Q.5
*/
#include <armadillo>
#include <iostream>
#include <sstream>
#include <cmath>
#include <ctime>

using namespace arma;

int main(void)
{
    clock_t ct = clock();
    ///-- set random seed
    arma::arma_rng::set_seed_random();
    ///-- load data
    mat dat;
    dat.load("fitinput.txt");
    const vec x = dat.col(0);
    const vec y = dat.col(1);
    ///-- @function dXda -> dX^2/da = sum{(ax^H - y)x^H}
    auto dXda = [x, y](const double aa, const double HH) -> double {
        return arma::sum((aa * arma::pow(x, HH) - y) % arma::pow(x, HH));
    };
    ///-- @function dXdH -> dX^2/dH = sum{(ax^H - y)* ax^H * ln(x)}
    auto dXdH = [x, y](const double aa, const double HH) -> double {
        return aa * arma::sum((aa * arma::pow(x, HH) - y) % arma::pow(x, HH) % arma::log(x));
    };
    ///-- @param M_MCMC, @param M_HMC -> number of steps
    const int M_MCMC = 1000;
    const int M_HMC = 10;
    ///-- @param a and @param H -> the proposed parameters
    double a = 1.0, H = 0.0;
    ///-- @param a_best and @param H_best -> accepted parameters
    double a_best = 1.0, H_best = 0.0;
    ///-- @param kesiA = da/dt
    ///-- @param kesiH = dH/dt
    double kesiA, kesiH;
    kesiA = arma::randn();
    kesiH = arma::randn();
    ///-- @param Xold and @param Kold are "potential" and "kinetic" energy corresponding to a_best and H_best
    ///-- @param Xnew and @param Knew are "potential" and "kinetic" energy corresponding to proposed a and H
    double Xold, Xnew;
    double Kold, Knew;
    Xold = arma::sum(arma::pow(a * arma::pow(x, H) - y, 2)) / 2.0;
    Kold = (kesiA * kesiA + kesiH * kesiH) / 2.0;
    ///--- loop MCMC method
    for (int i = 0; i < M_MCMC; ++i)
    {
        ///-- HMC step forward with Verlet algorithm
        for (int j = 1; j <= M_HMC; ++j)
        {
            double aprev = a, Hprev = H;
            double kesiAprev = kesiA, kesiHprev = kesiH;
            ///-- @param eps -> t
            double eps = 1e-05 * j;
            a = aprev + kesiAprev * eps - 0.5 * eps * eps * dXda(aprev, Hprev);
            H = Hprev + kesiHprev * eps - 0.5 * eps * eps * dXdH(aprev, Hprev);
            kesiA = kesiAprev - 0.5 * eps * (dXda(a, H) + dXda(aprev, Hprev));
            kesiH = kesiHprev - 0.5 * eps * (dXdH(a, H) + dXdH(aprev, Hprev));
        }
        //-- compute X for the proposed a and H
        Xnew = arma::sum(arma::pow(a * arma::pow(x, H) - y, 2)) / 2.0;
        //-- compute K for the new kesiA and kesiH
        Knew = (kesiA * kesiA + kesiH * kesiH) / 2.0;
        double DH = Xnew + Knew - Xold - Kold;
        if (exp(-DH) >= randu())
        {
            a_best = a;
            H_best = H;
            Xold = Xnew;
            Kold = Knew;
        }
        else
        {
            ///-- if the step is not accepted return a and H to its original values
            a = a_best;
            H = H_best;
        }
        ///-- arma::randn() = N(0,1)
        kesiA = arma::randn();
        kesiH = arma::randn();
    }
    //-- output and plot results
    std::cout.setf(std::ios_base::scientific);
    std::cout << "a_best = "<<a_best << "\tH_best = " << H_best 
                << "\terror = " << exp(-Xold / x.size()) << std::endl;
    std::stringstream gnuplot;
    gnuplot << "gnuplot -e \" set terminal png enhanced;"
            << "set xlabel 'x'; set ylabel 'y'; set grid;"
            << "set title 'HMC approximation: a = " << a_best << ", H = " << H_best << "';"
            << "set key left; set output 'HMCFit-res.png';"
            << "p 'fitinput.txt' u 1:2 w p pt 6 ps 1.0 lc 'red' t 'data',";
    gnuplot.setf(std::ios_base::scientific);
    gnuplot << "'' u 1:(" << a_best << "* \\$1 **" << H_best << ") w l lc 'blue' t 'ax^H'\"";
    if (system(gnuplot.str().c_str()))
    {
        std::cout << "error plot!" << std::endl;
        std::cout << gnuplot.str() << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Program 'HMCFit' finished in " <<
             (double)(clock() - ct)/CLOCKS_PER_SEC << " sec" << std::endl;
    return 0;
}
