/*!
*    Exercise set 15  
*    You can compile and run the program with main15.sh script
*    if not please first compile and run 'setFitConditions.cpp'
*    if the result wasn't satisfying try running the program again     
*/
/*!-- 
*   'fitConditions.h' Contains the constants produced 
*   with 'setFitConditions.cpp' (from the file 'fitConditions.txt' by default)
*   We needed to do this because bitset class length needs to be constexpr(the 
*   length should be declared before the compilation) 
*/
#include "fitConditions.h"
#include <armadillo>
#include <bitset>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>

typedef unsigned int uint;

template <typename GenomeVec>
GenomeVec RouletteWheel(const arma::vec &L, const GenomeVec &prvgen);

template <typename GenomeVec>
void cross_over(const double pc, GenomeVec &g);

template <typename GenomeVec>
void Mutation(const double pm, GenomeVec &g);

int main(int argc, char *argv[])
{
    const std::clock_t ct = std::clock();
    ///-- set random seed
    arma::arma_rng::set_seed_random();
    ///-- import data
    arma::mat dat;
    if (!dat.load("fitinput.txt"))
    {
        std::cout << "error! 'fitinput.txt' not found" << std::endl;
        return -1;
    }
    const arma::vec x = dat.col(0);
    const arma::vec y = dat.col(1);
    ///-- error function X^2 = sum((y - expected_y)^2) = sum((ax^H - expected_y)^2)
    auto ferror = [x, y](const arma::vec &a, const arma::vec &H) -> arma::vec {
        const unsigned M = (unsigned)a.size();
        arma::vec X(M);
        for (uint i = 0; i < M; ++i)
        {
            X[i] = arma::sum(arma::pow(a[i] * arma::pow(x, H[i]) - y, 2));
        }
        return X;
    };
    ///-- M-> the number of individuals in each generation M <= Na and NH preference M = 100
    const uint M = Na >= 100 && NH >= 100 ? 100 : std::min(Na, NH);
    ///-- initialize a and H randomly amin <= a <= amax, Hmin <= H <= Hmax
    arma::vec a = (amax - amin) * arma::randu<arma::vec>(M) + amin;
    arma::vec H = (Hmax - Hmin) * arma::randu<arma::vec>(M) + Hmin;
    ///-- likelihood of a and H according to error L(a,H) = exp(-(ax^H -y)^2/2sd^2)
    ///-- here because the standard deviations of vector x is too large sd considered sd = 500
    arma::vec L = arma::exp(-ferror(a, H) / 1000.0);
    ///-- currentLsum -> L~
    double currentLsum = arma::sum(L);
    L /= currentLsum;
    ///-- create a generation by coding each individual a and H (a vector of 'genomes')
    std::vector<std::bitset<na>> currentAgeneration(M);
    std::vector<std::bitset<nH>> currentHgeneration(M);
    ///-- coding -> b_a(i) = bin{round [(2^n - 1) (a_i - amin)/(amax - amin)]}
    ///--           b_H(i) = bin{round [(2^n - 1) (H_i - Hmin)/(Hmax - Hmin)]}
    for (uint i = 0; i < M; ++i)
    {
        currentAgeneration[i] = std::lround((pow(2, na) - 1) * (a[i] - amin) / (amax - amin));
        currentHgeneration[i] = std::lround((pow(2, nH) - 1) * (H[i] - Hmin) / (Hmax - Hmin));
    }
    double Convergence = 999.0;
    ///-- D -> Delta -> a small number
    const double D = 1e-06;
    // double pc = 1.0, pm = 0.1;
    ///-- pc(t) -> probability of cross_over (large)
    ///-- I supposed that the probabilities of cross_over and
    ///-- mutation decay exponentially as the "generations" pass
    auto pc = [](double t) -> double {
        return 1.0 * exp(-t / 1e04);
    };
    ///-- pm(t) -> probability of mutation (small)
    ///-- mutation's probability decays faster
    auto pm = [](double t) -> double {
        return 0.1 * exp(-t / 1e03);
    };
    ///-- generation counter
    uint t = 0;
    while (Convergence > D)
    {
        std::vector<std::bitset<na>> nextAgeneration;
        std::vector<std::bitset<nH>> nextHgeneration;
        ///--select the next "generation" with Roulette Wheel selection
        nextAgeneration = RouletteWheel(L, currentAgeneration);
        nextHgeneration = RouletteWheel(L, currentHgeneration);
        ///-- if the convergence is taking too long reduce the probabiliy of cross_over
        // if (t == 1000)
        // {
        //     pc = 0.7;
        // }
        ///-- cross over the generation with the probability pc(t)
        cross_over(pc(t), nextAgeneration);
        cross_over(pc(t), nextHgeneration);
        ///-- if the convergence is taking too long reduce the probabiliy of mutation
        // if (t == 1000)
        // {
        //     pm = 1e-03;
        // }
        ///-- Mutate "genomes" of the "generation" with the porbiliy pm(t)
        Mutation(pm(t), nextAgeneration);
        Mutation(pm(t), nextHgeneration);
        ///-- decoding -> a_i = amin + (amax - amin)bin^-1(ba_i)/(2^n -1)
        ///               H_i = Hmin + (Hmax - Hmin)bin^-1(bH_i)/(2^n -1)
        for (uint i = 0; i < M; ++i)
        {
            a[i] = nextAgeneration[i].to_ulong() * (amax - amin) / (pow(2, na) - 1.0) + amin;
            H[i] = nextHgeneration[i].to_ulong() * (Hmax - Hmin) / (pow(2, nH) - 1.0) + Hmin;
        }
        ///-- save the current "genomes" for the next "generation"
        currentAgeneration = nextAgeneration;
        currentHgeneration = nextHgeneration;
        ///-- calculate the new "generation"'s likelihood
        L = arma::exp(-ferror(a, H) / 1000.0);
        double nextLsum = arma::sum(L);
        L /= nextLsum;
        ///-- Convergence = | L~(t) - L~(t+1)|
        Convergence = std::abs(nextLsum - currentLsum);
        currentLsum = nextLsum;
        t++;
    }
    ///-- the best indexes for a and H are those which have maximum likelihood
    const arma::uvec best = arma::find(L == L.max());
    ///-- it is possible that multiple pairs of a and H have the same likelihood so we average over them
    const double a_best = arma::mean(a(best)), H_best = arma::mean(H(best));
    const arma::vec yres = a_best * arma::pow(x, H_best);
    ///-- output the results
    std::cout.setf(std::ios_base::scientific);
    std::cout << "GA converged in " << t << "'th generation" << std::endl;
    std::cout << "a_best = " << a_best << "\tGenes: " << currentAgeneration[best[0]]
              << "\nH_best = " << H_best << "\tGenes: " << currentHgeneration[best[0]]
              << "\nL_max = " << L.max() << std::endl;
    ///-- save the results
    std::ofstream fres("GAResults.txt");
    if (!fres)
    {
        std::cout << "error! opening file" << std::endl;
        return -1;
    }
    fres.setf(std::ios_base::scientific);
    for (uint i = 0; i < M; ++i)
    {
        fres << x[i] << "\t" << y[i] << "\t" << yres[i] << "\n";
    }
    fres.close();
    ///-- plot the results
    std::stringstream gp;
    gp << "gnuplot -e \" set terminal png enhanced; set grid;"
       << "set output 'GAResults.png'; set key left;"
       << "set xlabel 'x'; set ylabel 'y';"
       << "set title 'Simple Genetic Algorithm: a = " << a_best
       << " H = " << H_best << "';"
       << "p 'GAResults.txt' u 1:2 w p ps 0.5 pt 6 lc 'red' t 'data'"
       << ", '' u 1:3 w l lc 'blue' t 'ax^H';\"";
    if (system(gp.str().c_str()))
    {
        std::cout << "error plot!" << std::endl;
        std::cout << gp.str() << std::endl;
        return -2;
    }
    std::cout << "Program '" << argv[0] << "' finished in "
              << (double)(std::clock() - ct) / CLOCKS_PER_SEC << " sec" << std::endl;
    return 0;
}
///============================================================================
///-- GenomeVec -> a vector of genomes (like a matrix of 0 and 1's)
///--               this class must have a size() method and an overloadede [] operator
///-- @param L -> Likelihood of each individual
///-- @param prevgen -> previous generation's genomes
///-- returns the selected next generation with Roulette wheel selection
template <typename GenomeVec>
GenomeVec RouletteWheel(const arma::vec &L, const GenomeVec &prvgen)
{
    const unsigned M = (unsigned)L.size();
    GenomeVec nxtgen(M);
    for (uint i = 0; i < M - 1; ++i)
    {
        bool selected = false;
        double X = arma::randu();
        double sum1 = L[i];
        double sum2 = 0.0;
        ///-- role the wheel
        for (uint l = 0; l < M - 1; ++l)
        {
            if (sum2 <= X && X < sum1)
            {
                nxtgen[i] = prvgen[l];
                selected = true;
                ///-- if the individual is selected go to the next one
                break;
            }
            sum2 = sum1;
            sum1 += L[l + 1];
        }
        ///-- if no one is selected, select the last one
        if (!selected)
        {
            nxtgen[i] = prvgen[M - 1];
        }
    }
    return nxtgen;
}
///=============================================================
///-- takes a generation and cross over it with the probability pc
template <typename GenomeVec>
void cross_over(const double pc, GenomeVec &g)
{
    const uint M = (uint)g.size();
    const uint n = (uint)g[0].size();
    ///-- (M,2) = M(M-1)/2
    for (uint l = 0; l < (M * (M - 1)) / 2; ++l)
    {
        ///-- with a probability pc the y1 to y2'th genomes of individual X1
        ///-- will cross_over with the y1 to y2'th genomes of individual X2
        if (pc >= arma::randu())
        {
            ///-- choose two random individuals
            uint X1 = arma::randi(arma::distr_param(0, M - 1));
            uint X2 = arma::randi(arma::distr_param(0, M - 1));
        ///-- the same individual can't cross_over with itself
        l:
            if (X2 == X1)
            {
                X2 = arma::randi(arma::distr_param(0, M - 1));
                goto l;
            }
            uint y1 = arma::randi(arma::distr_param(0, n - 1));
            uint y2 = arma::randi(arma::distr_param(0, n - 1));
            for (uint k = std::min(y1, y2); k <= std::max(y1, y2); ++k)
            {
                uint aux = g[X1][k];
                g[X1][k] = g[X2][k];
                g[X2][k] = aux;
            }
        }
    }
}
///=====================================================================
///-- Mutates random individuals of a generation with the probabiliy pm
template <typename GenomeVec>
void Mutation(const double pm, GenomeVec &g)
{
    const uint M = (uint)g.size();
    const uint n = (uint)g[0].size();
    for (uint l = 0; l < M; ++l)
    {
        if (pm >= arma::randu())
        {
            ///-- select a random genome of the l'th individual and flip it
            uint X = arma::randi(arma::distr_param(0, n - 1));
            g[l][X] = !g[l][X];
        }
    }
}
