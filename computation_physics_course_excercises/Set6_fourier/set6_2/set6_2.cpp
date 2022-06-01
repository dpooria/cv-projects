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
int plotthis(string fname1, string title1, string xlabel,
             string ylabel,
             string fname2 = "",
             string title2 = "y(t)")
{
    string gnuplot = "gnuplot -e \"set terminal jpeg enhanced;\
    set xlabel '" + xlabel +
                     "'; set ylabel '" + ylabel + "';\
    set style line 1 lc rgb 'blue' pt 5 ps 0.8; set grid;\
    set style line 2 lc rgb 'black' pt 6 ps 0.8; set grid;p '" +
                     fname1 + ".txt' w p ls 1 t '" + title1 + "'";
    if (fname2.empty())
    {
        gnuplot = gnuplot + ";\" >'" + fname1 + ".jpg'";
    }
    else
    {
        gnuplot = gnuplot + ", '" + fname2 + ".txt' w p ls 2 t '" + title2 + "';\"\
        > '" + fname1 +
                  "-" + fname2 + ".jpg'";
    }
    return system(gnuplot.c_str());
}

int main()
{
    int N = 100;
    //-- create a random data
    vec randat = randu<vec>(N);
    //-- T_i => period time scale
    vec T(10);
    for (int i = 0; i < 10; i++)
    {
        T(i) = 2 * N / (i + 1.0);
    }
    //-- A_i => amplitude
    vec A = zeros<vec>(10);
    for (int i = 0; i < 10; i++)
    {
        for (int t = 0; t < N; t++)
        {
            A(i) += randat(t) * sin(M_PI * (i + 1) * t / (N + 0.0));
        }
        A(i) *= 2.0 / N;
    }
    //-- x => time sery
    vec x(N);
    for (int t = 0; t < N; t++)
    {
        x(t) = sum(A % arma::sin(2 * M_PI * t / T));
    }
    randat.save("randomData.txt", raw_ascii);
    x.save("x.txt", raw_ascii);
    int response = plotthis("x", "x(t)", "t", "x(t)", "randomData", "random Data");

    //-- set6_2_B: power spectrum is the fourier transform of correlation function
    cx_vec cxpsx = fft(correlcal(x));
    cxpsx.save("powerspectrumX.txt", raw_ascii);

    vec ipsx = imag(cxpsx);
    ipsx.save("imagXpowerSpecturm.txt", arma::raw_ascii);
    response += plotthis("imagXpowerSpecturm", "Power Spectrum of x(t)",
                         "index", "power spectrum");

    //-- f => fourier transform of x(t)
    cx_vec f = arma::fft(x);
    double meanf = mean(imag(f));

    //-- filter function
    cx_vec lpf(N);
    cx_vec hpf(N);
    for (int i = 0; i < N; i++)
    {
        /*-- if fft is bigger than its mean average -> lpf */
        if (f(i).imag() > meanf)
        {
            lpf(i) = f(i);
            hpf(i) = cx_double(0, 0);
        }
        else
        {
            lpf(i) = cx_double(0, 0);
            hpf(i) = f(i);
        }
    }

    //-- inverse fast fourier transform for deriving the new data
    vec xlpf = real(ifft(lpf));
    vec xhpf = real(ifft(hpf));
    xhpf.save("highpasseddata.txt", arma::raw_ascii);
    xlpf.save("lowpasseddata.txt", arma::raw_ascii);
    response += plotthis("highpasseddata", "Data after High Pass Filter",
                         "t", "x(t)", "x", "x(t)");
    response += plotthis("lowpasseddata", "Data after Low Pass Filter",
                         "t", "x(t)", "x", "x(t)");
    if (response != 0)
    {
        cout << "gnuplot: There was an error!\n";
    }
    return 0;
}