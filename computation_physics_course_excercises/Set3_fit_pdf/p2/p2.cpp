#include <matplotlibcpp.h>
#include <armadillo>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <time.h>

using namespace std;
using namespace arma;
namespace plt = matplotlibcpp;
typedef vector<double> stdvect;
typedef Col<double> armavect;

void slowpdf(const armavect &ku, const armavect &kg);
void fastpdf(const armavect &ku, const armavect &kg);

int main(int argc, char **argv)
{
    arma::arma_rng::set_seed_random();
    //--uniformly distributed random
    armavect ku = randu(1000, 1);
    //--normally distributed random
    armavect kg = randn(1000, 1);
    //--change the range of 0<ku<1 to range of kg
    ku = (kg.max() - kg.min()) * ku + kg.min();
    clock_t t = clock();
    fastpdf(ku, kg);
    cout << "t(fastpdf)\t" << (double)(clock() - t) / CLOCKS_PER_SEC << endl;
    t = clock();
    slowpdf(ku, kg);
    cout << "t(slowpdf)\t" << (double)(clock() - t) / CLOCKS_PER_SEC << endl;
    return 0;
}
//--fast pdf calculation
void fastpdf(const armavect &datu, const armavect &datg)
{
    const uint N = datg.size();
    double dxg = (datg.max() - datg.min()) / N;
    double dxu = (datu.max() - datu.min()) / N;
    ivec ku = conv_to<ivec>::from(datu / dxu);
    ivec kg = conv_to<ivec>::from(datg / dxg);
    double kminu = ku.min(), kmaxu = ku.max();
    double kming = kg.min(), kmaxg = kg.max();
    map<int, uint> ng, nu;
    for (int k = kming; k <= kmaxg; ++k)
    {
        ng.insert({k, 0});
     }
    for (uint i = 0; i < kg.size(); ++i)
    {
        ng.at(kg[i])++;
    }
    fstream resultg;
    resultg.open("gaussian.txt", ios::out);
    for (auto nitg:ng)
    {
        resultg << nitg.first * dxg << "\t" << nitg.second / (N * dxg) << endl;
    }
    resultg.close();

    for (int k = kminu; k <= kmaxu; ++k)
    {
        nu.insert({k, 0});
    }
    for (uint i = 0; i < ku.size(); ++i)
    {
        nu.at(ku[i])++;
    }
    fstream resultu;
    resultu.open("uniform.txt", ios::out);
    for (auto nitu : nu)
    {
        resultu << nitu.first * dxu << "\t" << nitu.second / (N * dxu) << endl;
    }
    resultu.close();
    string gnuplot= "gnuplot -e 'set terminal jpeg; set xlabel \"x\";\
        set ylabel \"p(x)\";set grid; p \"uniform.txt\" w p ls 1, \
        \"gaussian.txt\" w p ls 2;'> randompdf-fast.jpg";
    if(system(gnuplot.c_str())){
        cout << "error in plot!" << endl;
    }
}

//--slow pdf calculation
void slowpdf(const armavect &ku, const armavect &kg)
{
    double xming = kg.min(), xmaxg = kg.max();
    double xminu = ku.min(), xmaxu = ku.max();
    const uint N = kg.size();
    const double dxg = (xmaxg - xming) / N, dxu = (xmaxu - xminu) / N;
    stdvect sxg, sxu;
    sxg.push_back(xming);
    sxu.push_back(xminu);
    uint i = 0;
    do
    {
        sxg.push_back(sxg[i] + dxg);
        ++i;
    } while (sxg[i] < xmaxg);
    i = 0;
    do
    {
        sxu.push_back(sxu[i] + dxu);
        ++i;
    } while (sxu[i] < xmaxu);
    armavect ng = arma::zeros(sxg.size(), 1);
    armavect nu = arma::zeros(sxu.size(), 1);
    //--calculate pdf with simple algorithm
    for (i = 0; i < sxg.size(); ++i)
    {
        for (uint j = 0; j < kg.size(); ++j)
        {
            if (kg[j] >= sxg[i] && kg[j] < sxg[i + 1])
            {
                ++ng[i];
            }
        }
    }
    for (i = 0; i < sxu.size(); ++i)
    {
        for (uint j = 0; j < ku.size(); ++j)
        {
            if (ku[j] >= sxu[i] && ku[j] < sxu[i + 1])
            {
                ++nu[i];
            }
        }
    }
    plt::named_plot("Uniform", sxu, arma::conv_to<stdvect>::from(nu / N), "b*");
    plt::named_plot("Gaussian", sxg, arma::conv_to<stdvect>::from(ng / N), "r.");
    plt::legend();
    plt::xlabel("x");
    plt::ylabel("p(x)");
    plt::save("randompdf-slow.jpg");
    plt::close();
}
