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

double gwf(const double &);
int main(int argc, char **argv)
{
    double dx = 0.2;
    vec dat;
    dat.load("marks.txt");
    ivec k = arma::conv_to<ivec>::from(dat/dx);
    int kmin = k.min(), kmax = k.max();
    const uint N = k.size();
    //--create map to store the number of members of each class
    map<int, double> n;
    for (int i = kmin; i <= kmax; ++i)
    {
        n.insert({i, 0});
    }
    
    for (uint i = 0; i < N; ++i)
    {
        n.at(k[i])+=gwf((k[i]*dx-dat[i])/dx);
    }
   
    stdvect x,p;
    x.resize(n.size());
    p.resize(n.size());
    ofstream res;
    res.open("marksPDF.txt",ios::out);
    for(auto nit:n){
        x.push_back(nit.first*dx);
        p.push_back(nit.second/(N*dx));
        res << nit.first*dx << "\t" << nit.second/(N*dx) << endl; 
    }
    res.close();
    plt::plot(x,p,"b.");
    plt::plot(x,p,"r:");
    plt::xlabel("x");
    plt::ylabel("p(x)");
    plt::grid("on");
    plt::save("marks GWF PDF.png");
    plt::clf();
    return 0;
}
//--Gaussian window function
double gwf(const double &y)
{
    return exp(-powl(y, 2) / 8.0);
}