/*!
    Program of set3 to compare the results with MCMC and HMC
*/
#include <armadillo>
#include <iostream>
#include <sstream>
#include <ctime>

using namespace arma;

int main(int argc, char **argv)
{
    clock_t t = clock();
    //----Read the data-------
    mat dat;
    if(!dat.load("fitinput.txt"))
    {
        std::cout << "load file error!" << std::endl;
        exit(EXIT_FAILURE);
    }
    const vec x = dat.col(0);
    const vec y = dat.col(1);
    //--initialize a, H and X
    double a, H, X;
    double atmp = 0.0, Htmp = 0.0, Xtmp;
    X = arma::sum(arma::pow(y - atmp * arma::pow(x, Htmp), 2));
    //--fit data by searching the phase space
    //--search phase space for -10<a<10 && -3<H<3
    for (int i = 0; i < 200; ++i)
    {
        atmp += 0.005;
        Htmp = 0.0;
        for (int j = 0; j < 200; ++j)
        {
            Htmp += 0.005;
            Xtmp = arma::sum(arma::pow(y - atmp * arma::pow(x, Htmp), 2));
            //--find minimum of X and therefor values of a_best & H_best
            if (Xtmp < X)
            {
                X = Xtmp;
                a = atmp;
                H = Htmp;
            }
        }
    }
        //-- output results
        std::iostream::sync_with_stdio(false);
        std::cout.setf(std::ios_base::scientific);
        std::cout << "a_best = " << a << "\t"
                  << "H_best = " << H << "\t, error ="
                  << exp(-X / x.size()) << std::endl;
        std::stringstream gnuplot;
        gnuplot << "gnuplot -e \" set terminal png enhanced;"
                << "set grid; set key left; set xlabel 'x';"
                << "set ylabel 'y'; set output 'Set3Fit-res.png';"
                << "set title 'searching phase space: a = " << a
                << " , H = " << H << "';";
        gnuplot.setf(std::ios_base::scientific);
        gnuplot << "p 'fitinput.txt' u 1:2 w p pt 6 ps 1.0 lc 'red'"
                << " t 'data' , '' u 1:(" << a << " * \\$1**" << H
                << ") w l lc 'blue' t 'ax^H';\"";
        if(system(gnuplot.str().c_str()))
        {
            std::cout << "plot error!" << std::endl;
            std::cout << gnuplot.str() << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Program 'Set3Fit' finished in " << (double)(clock() - t)/CLOCKS_PER_SEC 
                    << " sec" << std::endl;
        return 0;
    }
