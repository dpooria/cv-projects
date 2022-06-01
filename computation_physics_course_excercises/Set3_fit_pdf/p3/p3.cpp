#include <matplotlibcpp.h>
#include <armadillo>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <time.h>

using namespace arma;
using namespace std;
typedef vector<double> stdvect;
namespace plt = matplotlibcpp;

void fpdfcal(const vec &data, const double &dx, string &name);
vector<stdvect> spdfcal(const vec &data, const double &dx);
int main(int argc, char **argv)
{
    vec dat;
    dat.load("data.txt");
    uint N = dat.size();
    double dx1 = 0.1, dx2 = 0.01, dx3 = 0.001;
    string res1 = "dx1", res2 = "dx2", res3 = "dx3";
    clock_t t = clock();
    fpdfcal(dat, dx1, res1);
    fpdfcal(dat, dx2, res2);
    fpdfcal(dat, dx3, res3);
    cout << "t(fpdfcal)\t" << ((double)(clock() - t))/CLOCKS_PER_SEC << "sec" << endl;
    
    //--this part has been commented out intentionally due to its long runtime(~20sec)
    // t = clock();
    // vector<stdvect> spdf1=spdfcal(dat, dx1);
    // vector<stdvect> spdf2=spdfcal(dat, dx2);
    // vector<stdvect> spdf3=spdfcal(dat, dx3);
    // cout << "t(spdfcal)\t" << double((clock() - t) / CLOCKS_PER_SEC) << " sec" << endl;
    // plt::named_plot("$\\Delta$x=0.001",spdf3[0],spdf3[1],"b.");
    // plt::named_plot("$\\Delta$x=0.01",spdf2[0],spdf2[1],"r.");
    // plt::named_plot("$\\Delta$x=0.1",spdf1[0],spdf1[1],"k*");
    // plt::xlabel("x");
    // plt::ylabel("p(x)");
    // plt::legend();
    // plt::save("all.jpg");
    cout << "---------done!---------" << endl;
    return 0;
}
//--fast pdf calculation
void fpdfcal(const vec &data, const double &dx, string &name)
{
    ivec k = conv_to<ivec>::from(data / dx);
    int kmin = k.min(), kmax = k.max(), N = k.size();
    //--create a map with classes as key and number of members as value
    map<int, uint> n;
    //--initialize map with zeros
    for (int i = kmin; i <= kmax; ++i)
    {
        n.insert({i, 0});
    }
    //--loop trough data and count the number of members of each class
    for (uint i = 0; i < N; ++i)
    {
        n.at(k[i])++;
    }
    //--write x and p(x) to file name.txt
    fstream result;
    result.open(name + ".txt", ios::out);
    for (auto nit = n.begin(); nit != n.end(); ++nit)
    {
        result << (*nit).first * dx << "\t" << (*nit).second / (N * dx) << endl;
    }
    result.close();

    string gnuplot = "gnuplot -e 'set terminal jpeg; set xlabel \"x\";";
    gnuplot += "set ylabel \"p(x)\";set grid; p \"" + name + ".txt\" w p'>" + name + ".jpg";
    int response = system(gnuplot.c_str());
    if (response)
    {
        cout << "error!" << endl;
    }
}
//--slow pdf calculation
vector<stdvect> spdfcal(const vec &data, const double &dx)
{
    double xmin = data.min(), xmax = data.max();
    uint N = data.size();
    stdvect x,n;
    x.push_back(xmin);
    n.push_back(0);
    uint i = 0;
    do
    {
        x.push_back(x[i] + dx);
        n.push_back(0);
        for (int j = 0; j < N; ++j)
        {
            if (data[j] >= x[i] && data[j] < x[i + 1])
            {
                ++n[i];
            }
        }
        n[i] /= N * dx;
        ++i;
    } while (x[i] < xmax);
    
    return vector<stdvect>({x,n});
}
