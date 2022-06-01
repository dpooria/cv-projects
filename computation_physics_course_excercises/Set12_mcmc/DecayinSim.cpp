/*!
    Program for Exersice set 12 Q. 1
*/

#include <functional>
#include <fstream>
#include <iostream>
#include <random>
#include <ctime>
#include <cmath>

using namespace std;
///-- @param dt1 time step for case A
const double dt1 = 0.1;
///-- @param dt2 time step for case B ( p(t) = l2 * dt2 / t )
const double dt2 = 1.0;
///-- @param N0 number of particles
const int N0 = 1000;
///-- @param l1 decay constant (case A)
const double l1 = 1.1e-03;
///-- @param l2 decay constant of
const double l2 = 0.8;
///-- @function pA -> decay probability function for case A
const double pA(const double t)
{
    return l1 * dt1;
}
///-- @function pA -> decay probability function for case B
const double pB(const double t)
{
    return l2 * dt2 / t;
}
/*! @function DecaySim simulates decaying
        @param N0 -> initial number of particles
        @param dt -> time step
        @param p -> decay probability function
        @param Nt -> theory function
        @param fname -> file name to save the results in
*/
void DecaySim(const int N0, const double dt, const double (*p)(double),
              function<double(double)> Nt, const char *fname)
{

    ofstream fres(fname);
    default_random_engine g;
    uniform_real_distribution<double> d;
    if (!fres.is_open())
    {
        cout << "error opening file!" << endl;
        return;
    }
    fres.setf(ios_base::scientific);
    ///-- for t = 0 number of particle = N0
    int N = N0;
    double t = 0.0;
    fres << t << "\t" << N << "\t" << N0 << endl;
    g.seed((unsigned)time(0));
    ///-- repeat  the simulation until all of particles has been decayed
    while (N > 0)
    {
        /// step forward
        t += dt;
        ///-- @param dN -> number of particles that will decay in this time interval
        int dN = 0;
        double pa = p(t);
        for (int i = 1; i <= N; ++i)
        {
            ///-- d(g) -> random number 0 <= R <= 1
            if (d(g) <= pa)
            {
                dN++;
            }
        }
        N -= dN;
        fres << t << "\t" << N << "\t" << Nt(t) << "\n";
    }
    fres.close();
}

int main(void)
{
    clock_t ct = clock();
    DecaySim(
        N0, dt2, pB, [=](double t) { return N0 * pow(t, -l2); }, "DecayinSim-resB.txt");
    DecaySim(
        N0, dt1, pA, [=](double t) { return N0 * exp(-l1 * t); }, "DecayinSim-resA.txt");
    const char *pltfrmt = R"E( gnuplot -e "
    set terminal png enhanced; set grid;
    set xlabel 't'; set ylabel 'N';
    set output 'DecayinSim-resA.png';
    set title 
    'p(t) = {/Symbol l}{/Symbol D}t, N0 = %d, {/Symbol l} = %4.2lf, {/Symbol D}t=%.2lf';
    p 'DecayinSim-resA.txt' u 1:2 w l lc 'blue' t 'simulation result',
        '' u 1:3 w l lc 'red' t 'N_0e^{-{/Symbol l}t}';
    set output 'DecayinSim-resB.png';
    set title 
    'p(t) = {/Symbol l}{/Symbol D}t/t, N0 = %d, {/Symbol l} = %4.2lf, {/Symbol D}t = %.2lf';
    p 'DecayinSim-resB.txt' every ::0::1000 u 1:3 w l lc 'red' t 'N_0t^{-{/Symbol l}}',
     '' every ::0::1000 u 1:2 w l lw 1.5 lc 'blue' t 'simulation result';"
    )E";
    char gnuplot[800];
    snprintf(gnuplot, 800, pltfrmt, N0, l1, dt1, N0, l2, dt2);
    if (system(gnuplot))
    {
        cout << "plot error!" << endl;
        cout << gnuplot << endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Program 'DecayinSim' finished in " <<
             (double)(clock() - ct)/CLOCKS_PER_SEC << " sec" << std::endl;
    return 0;
}