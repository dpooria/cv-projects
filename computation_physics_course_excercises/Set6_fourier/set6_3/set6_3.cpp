#include <armadillo>
#include <complex>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
using namespace std;
using namespace arma;

//-- calculates the correlation function of vector x
vec correlcal(const vec &x)
{
    size_t n = x.size();
    vec c = zeros<vec>(n);
    for (uint t = 0; t < n; t++)
    {
        for (auto i = 0; i < n - t; i++)
        {
            c(t) += x(i) * x(i + t);
        }
        c(t) /= n - t;
    }
    return c;
}
int plotthis(string fname1, string title1 = "x(t)", string xlabel = "t",
             string ylabel = "x",
             string fname2 = "",
             string title2 = "y(t)")
{
    string gnuplot = "gnuplot -e \"set terminal jpeg enhanced;\
    set xlabel '" + xlabel +
                     "'; set ylabel '" + ylabel + "';\
    set style line 2 lc rgb 'blue' pt 5 ps 0.8; set grid;\
    set style line 1 lc rgb 'black' pt 6 ps 0.8; set grid;p '" +
                     fname1 + ".txt' w p ls 1 t '" + title1 + "'";
    if (fname2.empty())
    {
        gnuplot = gnuplot + ";\" >'" + fname1 + ".jpg'";
    }
    else
    {
        gnuplot = gnuplot + ", '" + fname2 + ".txt' w p ls 2 t '" + title2 + "';\"\
        > '" + fname2 +
                  "and.jpg'";
    }
    return system(gnuplot.c_str());
}

void set6_3_ABC(int N)
{
    string dir = "res" + to_string(N);
    string mkdir = "mkdir -p " + dir;
    if (system(mkdir.c_str()))
    {
        cout << "mkdir : error!\n";
    }
    vec t = linspace(0, 5, N);
    //-- create a random sery
    vec randat = randu<vec>(N);
    vec A = zeros<vec>(10);
    double s;
    for (int i = 0; i < 10; i++)
    {
        s = 0;
        for (int j = 0; j < N; j++)
        {
            A(i) += randat(j) * sin(t(j) * M_PI / (i + 1.0));
            s += powl(sin(t(j) * M_PI / (i + 1.0)), 2);
        }
        A(i) /= s;
    }
    //-- create the sinusoidal time sery
    vec x = zeros<vec>(N);
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < 10; i++)
        {
            x(j) += A(i) * sin(t(j) * M_PI / (i + 1.0));
        }
    }
    randat.save(dir + "/randat.txt", arma::raw_ascii);
    x.save(dir + "/x(sinusoidalSery).txt", arma::raw_ascii);
    plotthis(dir + "/randat", "random data", "t/dt",
             "x(t)", dir + "/x(sinusoidalSery)", "x");
    cx_vec cx_psx = fft(correlcal(x));
    cx_psx.save(dir + "/PSofX.txt", arma::raw_ascii);
    vec psx = imag(cx_psx);
    psx.save(dir + "/imagPSofX.txt", arma::raw_ascii);
    plotthis(dir + "/imagPSofX", "Power Spectrum of X", "index",
             "x(t)");
}

int main(int argc, char **argv)
{
    //--set6_3_A:C
    set6_3_ABC(100);
    //-- f = 10Hz
    set6_3_ABC(50);
    //-- f = 15Hz
    set6_3_ABC(75);

    //--set6_3_D:
    mat dat;
    dat.load("sunspot.txt");
    vec sunspot = dat.col(1);
    vec pssunspot = real(fft(correlcal(sunspot)));
    pssunspot = log(pssunspot);
    pssunspot.save("LogPSofSunspot.txt", raw_ascii);
    plotthis("LogPSofSunspot", "power spectrum of sunspot.txt", "t/dt", "PS");
    sunspot.save("sunspotcol2.txt", raw_ascii);
    plotthis("sunspotcol2", "sunspot.txt");
    return 0;
}