/* Exercise set 8 part 2 */
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>

using namespace std;
typedef vector<double> vec;

vec solveExpilicit(double (*df)(double), double f0, double dx, int N)
{
    vec f(N, 0.0);
    f[0] = f0;
    for (int i = 0; i < N - 1; ++i)
    {
        f[i + 1] = f[i] + dx * df(f[i]);
    }
    return f;
}

int singlePlot(string name, string title, string res)
{
    string gnuplot =
        "gnuplot -e \"set terminal jpeg enhanced;set xlabel 'x';set ylabel 'f(x)';set grid;\
        set style line 1 lc rgb 'blue' pt 6 ps 0.8;\
        p '" +
        name + "' u 1:2 w lp ls 1 t '" + title + "';\" > '" + res + "'";
    return system(gnuplot.c_str());
}

int plot(string name1, string title1, string name2, string title2, string res)
{
    string gnuplot =
        "gnuplot -e \"set terminal jpeg enhanced;set xlabel 'x';set ylabel 'f(x)';set grid;\
        set style line 1 lc rgb 'red' pt 5 ps 0.8;\
        set style line 2 lc rgb 'blue' pt 6 ps 0.8;\
        p '" +
        name1 + "' u 1:2 w lp ls 1 t '" + title1 + "','" +
        name2 + "' u 1:2 w lp ls 2 t '" + title2 + "';\" > " + res;
    return system(gnuplot.c_str());
}
int main(int argc, char *argv[])
{
    const int N{100};
    const double dx{0.5}, f0{1.0}, x0{1.0};
    auto dfA = [](double ff) {
        return (double)powl(ff, 2.0);
    };
    auto dfB = [](double ff) {
        return -ff;
    };
    //-- Solution of equation part A with explicit method
    vec fAE = solveExpilicit(dfA, f0, dx, N);
    //-- Solution of part B with explicit method
    vec fBE = solveExpilicit(dfB, f0, dx, N);
    //-- Solution of equation part A with implicit method(refer to readme8.pdf)
    /* Since Delta = 1 - 4dxf(x) will be below zero for some x, f considered to be a complex vector;
    there is to possible solutions according to readme8.pdf one with 1 + sqrt(Delta) -> fAI1 and other with
    1 - sqrt(Delta) -> fAI2. 
    */
    vector<complex<double>> fAI1(N, complex<double>(0.0, 0.0));
    fAI1[0] = complex<double>(f0, 0.0);
    vector<complex<double>> fAI2(N, complex<double>(0.0, 0.0));
    fAI2[0] = complex<double>(f0, 0.0);
    complex<double> Delta(0.0, 0.0);
    for (int i = 0; i < N - 1; ++i)
    {
        Delta = 1.0 - 4.0 * dx * fAI1[i];
        fAI1[i + 1] = (1.0 + sqrt(Delta)) / (2.0 * dx);
        Delta = 1.0 - 4.0 * dx * fAI2[i];
        fAI2[i + 1] = (1.0 - sqrt(Delta)) / (2.0 * dx);
    }

    //-- Solution of equation part B with implicit method(refer to readme8.pdf)
    vec fBI(N, 0.0);
    fBI[0] = f0;
    for (int i = 0; i < N - 1; ++i)
    {
        fBI[i + 1] = fBI[i] / (1.0 + dx);
    }
    //-- Save results
    ofstream AE, AI1, AI2, AI3, BE, BI, AI4;
    AE.open("eAExplicit.txt");
    //-- eAImplicit1.txt -> |fAI1|
    AI1.open("eAImplicit1.txt");
    //-- eAImplicit4.txt -> real(fAI1)
    AI4.open("eAImplicit4.txt");
    //-- eAImplicit2.txt -> real(fAI2)
    AI2.open("eAImplicit2.txt");
    //-- eAImplicit3.txt -> imag(fAI1)
    AI3.open("eAImplicit3.txt");
    BE.open("eBExplicit.txt");
    BI.open("eBImplicit.txt");
    double x = x0;
    for (int i = 0; i < N; ++i)
    {
        AE << x << "\t" << fAE[i] << "\n";
        AI1 << x << "\t" << abs(fAI1[i]) << "\n";
        AI4 << x << "\t" << abs(fAI1[i].real()) << "\n";
        AI2 << x << "\t" << fAI2[i].real() << "\n";
        AI3 << x << "\t" << fAI1[i].imag() << "\n";
        BE << x << "\t" << fBE[i] << "\n";
        BI << x << "\t" << fBI[i] << "\n";
        x += dx;
    }
    AE.close();
    AI1.close();
    AI2.close();
    BE.close();
    BI.close();
    AI3.close();
    AI4.close();
    int respond{0};
    respond += plot("eAExplicit.txt", "Explicit", "eAImplicit1.txt", "Implicit", "resultA1.jpg");
    respond += plot("eBExplicit.txt", "Explicit", "eBImplicit.txt", "Implicit", "resultB.jpg");
    respond += singlePlot("eAExplicit.txt", "part A with Explicit method", "resA-Explicit.jpg");
    respond += singlePlot("eAImplicit1.txt", "Implicit |fAI1|", "resA-Implicit1.jpg");
    respond += singlePlot("eAImplicit2.txt", "Implicit Real(fAI2)", "resA-Implicit2.jpg");
    respond += singlePlot("eAImplicit3.txt", "Implicit Imag(fAI1)", "resA-Implicit3.jpg");
    respond += singlePlot("eAImplicit4.txt", "Implicit Real(fAI1)", "resA-Implicit4.jpg");

    if (respond)
    {
        cerr << "Plot Error!" << endl;
    }
    return 0;
}
