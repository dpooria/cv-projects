/*!
    Exercise set 12 Q.4
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
    mat dat;
    dat.load("fitinput.txt");
    const vec x = dat.col(0);
    const vec y = dat.col(1);
    ///-- current step
    double a_best = 1.0, H_best = 0.0, Xcur;
    ///-- suppose central limit theorem applies here (Normal distribution)
    Xcur = sum(pow(a_best * pow(x, H_best) - y, 2)) / 2.0;
    ///-- for the next step(proposed step)
    double Xpr;
    double a, H;
    arma::arma_rng::set_seed_random();
    const int M_MCMC = 100000;
    for (int i = 0; i < M_MCMC; ++i)
    {
        ///-- a and H considered to be normally distributed(al_ij = randn())
        a = a_best + randn();
        H = H_best + randn();
        Xpr = sum(pow(a * pow(x, H) - y, 2)) / 2.0;
        ///-- @param r uniform random number 0 <= r <= 1
        double r = randu();
        ///-- Metropolis algorithm: p_i/p_j > 1 ? (p_i/p_j = exp((X^2 - X_next^2)/2sd) > r)
        ///-- is next step accepted?
        if (exp((Xcur - Xpr) / 2.0) > r)
        {
            a_best = a;
            H_best = H;
            Xcur = Xpr;
        }
    }
    ///-- output results
    std::cout.setf(std::ios_base::scientific);
    std::cout << "a_best = " << a_best << "\t"
              << "H_best = " << H_best << "\t"
              << "error = " << exp(Xcur / x.size()) << std::endl;
    std::stringstream gnuplot;
    gnuplot << "gnuplot -e \""
            << "set terminal png; set key left;"
            << "set grid; set title 'MCMC approximation: a = "
            << a_best
            << ", H = " << H_best << "';";
    gnuplot.setf(std::ios_base::scientific);
    gnuplot << "set xlabel 'x'; set ylabel 'y';"
            << "set output 'MCMCFit-res.png';"
            << "p 'fitinput.txt' u 1:2 w p pt 6 ps 1.0 lc 'red' t 'data', '' u 1:("
            << a_best << "*(\\$1**" << H_best
            << ")) w l lc 'blue' t 'ax^H';\"";
    if (system(gnuplot.str().c_str()))
    {
        std::cout << "Plot error!" << std::endl;
        std::cout << gnuplot.str() << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Program 'MCMCFit' finished in " <<
             (double)(clock() - ct)/CLOCKS_PER_SEC << " sec" << std::endl;
    
    return 0;
}