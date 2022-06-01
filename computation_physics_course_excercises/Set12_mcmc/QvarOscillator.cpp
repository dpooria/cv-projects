/*!
    Program for Exersice set 12 Q. 3
*/

#include <random>
#include <fstream>
#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;

int main(void)
{
    clock_t ct = clock();
    /*!
    @param E energy in each state(alpha)
    @param x positions
    @param r1, r2 random numbers between 0 and 1
    @param E0 the ground state energy (E_min)
    @param al_best parameter corresponding to the ground state energy
    @param m # of positions
    */
    iostream::sync_with_stdio(false);
    double E, x, r1, r2;
    double al_best = -1.0, E0 = 999.0;
    const int m = 10000;
    ///-- open a file to save the result
    ofstream fres("QvarOscillator-res.txt");
    if (!fres.is_open())
    {
        cout << "opening file error!" << endl;
        exit(EXIT_FAILURE);
    }
    fres.setf(ios_base::scientific);
    default_random_engine g;
    g.seed((unsigned)time(0));
    uniform_real_distribution<double> d;
    /// @param al the parameter alpha supposed to be between 0.1, 2.0 and d(al) = 0.01
    for (double al = 0.1; al <= 2.0; al += 0.01)
    {
        E = 0.0;
        for (int i = 0; i < m; ++i)
        {
            ///-- p(al,x) ~ exp(-2*al*x^2) so x is generated with Box-Muller method
            r1 = d(g);
            r2 = d(g);
            x = sqrt(-2.0 * log(1 - r1)) * cos(2.0 * M_PI * r2);
            ///-- in Box-Muller method sd = 1 but here sd_x = 1/sqrt(4*al) => x = sd_x * x_boxMuller
            x /= sqrt(4.0*al);
            ///-- E = sum{al * (1/2 - 2*al^2)x^2}/m
            E += al + x * x * (0.5 - 2.0 * al * al);
        }
        E /= m;
        fres << al << "\t" << E << "\n";
        ///-- finding the minimum energy and the correspondence alpha
        if (E < E0)
        {
            E0 = E;
            al_best = al;
        }
    }
    fres.close();
    cout << "al_best = " << al_best << "\tE_0 = " << E0 << endl;
    ///-- plot the result
    const char *pltfrmt = R"E(gnuplot -e "
    set terminal png enhanced; set grid;
    unset key; set xlabel '{/Symbol a}';
    set ylabel 'E({/Symbol a})'; 
    set output 'QvarOscillator-res.png'; 
    set title '{/Symbol a}_{best} = %5.3lf, E_0 = %5.3lf';
    p 'QvarOscillator-res.txt' u 1:2 w l lc 'blue', 
        '' u 1:2 w p pt 6 ps 0.5 lc 'red';" )E";
    char gnuplot[500];
    snprintf(gnuplot, 500, pltfrmt, al_best, E0);
    if (system(gnuplot))
    {
        cout << "plot error!" << endl;
        cout << gnuplot << endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Program 'QvarOscillator' finished in " <<
             (double)(clock() - ct)/CLOCKS_PER_SEC << " sec" << std::endl;
    return 0;
}