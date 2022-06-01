#include <matplotlibcpp.h>
#include <armadillo>
#include <map>
#include <vector>
#include <iostream>
#include <math.h>

using namespace std;
using namespace arma;
namespace plt = matplotlibcpp;
typedef vector<double> stdvect;

double gwf(const double &, const double &);
void pdfcal(const vec &dat, double &s, string name);
int main(int argc, char **argv)
{
    vec dat;
    dat.load("marks.txt");
    double s1 = 2, s2 = 0.2;
    string res1 = "s1", res2 = "s2";
    pdfcal(dat, s1, res1);
    pdfcal(dat, s2, res2);
    return 0;
}
//--Gaussian window function
double gwf(const double &y, const double &s)
{
    return exp(-powl(y, 2) / (2 * powl(s, 2))); /*/sqrt(M_PI*8);*/
}

void pdfcal(const vec &dat, double &s, string name)
{
    double dx = 0.1;
    ivec k = arma::conv_to<ivec>::from(dat / dx);
    int kmin = k.min(), kmax = k.max();
    const uint N = k.size();
    //--map for calculating pdf by using top-hat window function
    map<int, uint> n;
    for (int i = kmin; i <= kmax; ++i)
    {
        n.insert({i, 0});
    }
    //--map for normalize the n with gaussian window function
    map<int, double> gn;
    for (int i = kmin; i <= kmax; ++i)
    {
        gn.insert({i, 0});
    }

    for (uint i = 0; i < N; ++i)
    {
        n.at(k[i])++;
    }
    int nmin = 0, nmax = 0;
    for (int i = kmin; i <= kmax; ++i)
    {
        nmin = min(n.at(i), uint(i - kmin));
        nmax = min(n.at(i), uint(kmax - i));
        for (int j = i - nmin; j <= i + nmax; ++j)
        {
            gn.at(i) += gwf(i * dx - j * dx, s) * n.at(j) / (N * dx);
        }
    }
    stdvect x, p;
    x.resize(gn.size());
    p.resize(gn.size());
    ofstream res;
    res.open(name + ".txt", ios::out);
    for (auto gnit : gn)
    {
        x.push_back(gnit.first * dx);
        p.push_back(gnit.second / (N * dx));
        res << gnit.first * dx << "\t" << gnit.second / (N * dx) << endl;
    }
    res.close();
    plt::named_plot("$\\sigma$ = "+to_string(_Float32(s)),x, p, "b.");
    plt::legend();
    plt::plot(x, p, "r:");
    plt::xlabel("x");
    plt::ylabel("p(x)");
    plt::grid("on");
    plt::save(name + ".png");
    plt::clf();
}
