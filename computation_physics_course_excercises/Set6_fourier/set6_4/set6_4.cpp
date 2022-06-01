#include <armadillo>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace arma;
using namespace std;

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

mat pdfcal(const vec &x)
{
    double xmin = x.min();
    double xmax = x.max();
    uint N = x.size();
    double dx = (xmax - xmin) / N;
    int kmin = static_cast<int>(xmin / dx), kmax = static_cast<int>(xmax / dx);

    std::map<int, uint> n;
    for (int k = kmin; k <= kmax; ++k)
    {
        n.insert({k, 0});
    }
    for (uint i = 0; i < N; ++i)
    {
        n.at(static_cast<int>(x.at(i) / dx))++;
    }
    vec p(N), xc(N);
    int i = 0;
    for (auto nit : n)
    {
        xc(i) = nit.first * dx;
        p(i) = nit.second / (N * dx);
        i++;
    }
    return join_horiz(xc, p);
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

vec generateGaussi(const vec &R1, const vec &R2, double al)
{
    int N = R1.size();
    vec C(N);
    vec t = linspace(0, 5, N);
    for (int i = 0; i < N; i++)
    {
        C(i) = (powl(t(i), 2.0) + 1.0) * exp(-al * powl(t(i), 2.0));
    }
    C.save("expectedCorrelation" + to_string(int(al)) + ".txt", raw_ascii);
    cx_vec S = fft(C);
    S *= (2.0 * M_PI) / (N * (t(1) - t(0)));
    vector<cx_double> X(N);

    for (uint i = 0; i < N; i++)
    {
        X.push_back(sqrt(S(i) / 2.0) * cx_double(R1(i), R2(i)));
    }
    ofstream snorm;

    vec x = real(ifft(conv_to<cx_vec>::from(X)));
    snorm.open("expectedPDF" + to_string(int(al)) + ".txt", ios::out);
    for (auto sit : S)
    {
        snorm << norm(sit) << "\n";
    }
    snorm.close();
    return x;
}

int main()
{
    int N = 100;
    //-- Generate random numbers
    vec R1 = randu<vec>(N);
    vec R2 = randu<vec>(N);
    //-- alpha = 1
    vec x1 = generateGaussi(R1, R2, 1);
    x1.save("x1.txt", arma::raw_ascii);
    plotthis("x1");
    mat p1 = pdfcal(x1);
    p1.save("pdfx1.txt", arma::raw_ascii);
    plotthis("pdfx1", "PDF of x1", "x", "p(x)", "expectedPDF1", "Expected PDF Function of x1");
    vec c1 = correlcal(x1);
    c1.save("correlationx1.txt", arma::raw_ascii);
    plotthis("correlationx1", "Correlation of x1", "{/Symbol t}", "C({/Symbol t})",
             "expectedCorrelation1", "Expected Correlation");
    //-- alpha = 2
    vec x2 = generateGaussi(R1, R2, 2);
    x2.save("x2.txt", arma::raw_ascii);
    plotthis("x2");
    mat p2 = pdfcal(x1);
    p2.save("pdfx2.txt", arma::raw_ascii);
    plotthis("pdfx2", "PDF of x2", "x", "p(x)", "expectedPDF2", "Expected PDF Function of x2");
    vec c2 = correlcal(x1);
    c2.save("correlationx2.txt", arma::raw_ascii);
    plotthis("correlationx2", "Correlation of x2", "{/Symbol t}", "C({/Symbol t})",
             "expectedCorrelation2", "Expected Correlation for x2");

    return 0;
}