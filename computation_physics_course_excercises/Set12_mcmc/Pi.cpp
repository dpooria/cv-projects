/*!
    Program for set 12 Q.2
*/
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
#include <ctime>
#include <cmath>

using namespace std;

double calPi(const int N)
{
    default_random_engine g;
    g.seed((unsigned)time(0));
    uniform_real_distribution<double> d;
    int N0 = 0;
    for (int i = 1; i <= N; ++i)
    {
        double x = d(g);
        double y = d(g);
        if (y <= sqrt(1.0 - x * x))
        {
            N0++;
        }
    }
    double pi = 4.0 * N0 / (N + 0.0);
    return pi;
}

int main(void)
{
    clock_t ct = clock();
    iostream::sync_with_stdio(false);
    ofstream fres("Pi-result.txt");
    if(!fres.is_open())
    {
        cout << "Error opening file!" << endl;
    }
    fres.setf(ios_base::scientific);
    for(int N = 100; N <= (int)1e7; N *= 2)
    {
        double pi = calPi(N);
        fres << log(N)/log(10.) << "\t" << pi << "\n";
    }
    fres.close();
    stringstream gnuplot;
    gnuplot << "gnuplot -e \" set terminal png enhanced;"
            << "set xlabel 'log_{10}N'; set ylabel '{/Symbol p}(N)';"
            << "set output 'Pi-result.png'; set grid; set key off;"
            << "p 'Pi-result.txt' u 1:2 w p ps 1.5 pt 6 lc 'red',"
            << M_PI << " w l lc 'blue';\"";

    if(system(gnuplot.str().c_str()))
    {
        cout << "plot error!" << endl;
        cout << gnuplot.str() << endl;
    }
    std::cout << "Program 'Pi' finished in " <<
             (double)(clock() - ct)/CLOCKS_PER_SEC << " sec" << std::endl;
    return 0;
}