#include <matplotlibcpp.h>
#include <armadillo>
#include <time.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;
namespace plt = matplotlibcpp;
typedef vector<double> stdvect;
typedef arma::Col<double> armavect;

stdvect fitsph(const armavect &, const armavect &);

int main(int argc, char **argv)
{
    //----Read the data-------
    ifstream fitinput;
    fitinput.open("fitinput.txt");
    double xtmp, ytmp;
    stdvect x, y;
    while (fitinput >> xtmp && fitinput >> ytmp)
    {
        x.push_back(xtmp);
        y.push_back(ytmp);
    }
    fitinput.close();

    //--convert vector<double> to arma::Col<double> vector
    armavect xa = arma::conv_to<armavect>::from(x);
    armavect ya = arma::conv_to<armavect>::from(y);

    plt::named_plot("data", x, y, "b.");
    plt::named_plot("a$x^H$", x, fitsph(xa, ya), "r-");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::legend();
    plt::grid(true);
    plt::save("fitinput.png");

    //--calcualte error of a and H
    double a, H;
    fstream sph;
    sph.open("sph.txt", ios::in);
    sph >> a >> H;
    sph.close();
    int N = xa.size();
    double meanx = arma::mean(xa), meany = arma::mean(ya);
    double ex, ey, ea, eH;
    ex = arma::sum(arma::pow(xa - meanx, 2)) / (N * (N - 1.0));
    ey = arma::sum(arma::pow(ya - meany, 2)) / (N * (N - 1.0));
    eH = ex * arma::sum(arma::pow(arma::log(ya / a) / (xa % arma::pow(arma::log(xa), 2)), 2)) +
         ey * arma::sum(arma::pow(1 / (arma::log(xa) % ya), 2));
    eH /= N;
    armavect diffay = arma::pow(1 / arma::pow(xa, H), 2);
    ea = ex * pow(H, 2) * arma::sum(arma::pow(ya / arma::pow(xa, H + 1), 2)) +
         ey * arma::sum(diffay);
    ea /= N;
    cout << "a_best = " << a << "\t"
         << "H_best = " << H << endl;
    cout << "error x and y:\nex = " << sqrt(ex) << "\t"
         << "ey = " << sqrt(ey) << endl;
    cout << "error a and H:\nea= " << sqrt(ea) << "\t"
         << "eH= " << sqrt(eH) << endl;
    sph.open("sph.txt", ios::out | ios::app);
    sph << sqrt(ea) << "\t" << sqrt(eH) << endl;
    sph.close();
    return 0;
}

//--fit data by searching the phase space
stdvect fitsph(const armavect &x, const armavect &y)
{
    //--initialize a, H and X
    double a, H, X;
    double atmp = -10.0, Htmp = -3.0, Xtmp;
    X = arma::sum(arma::pow(y - atmp * arma::pow(x, Htmp), 2));
    //--search phase space for -10<a<10 && -3<H<3
    clock_t t = clock();
    for (int i = 0; i < 200; ++i)
    {
        atmp += 0.1;
        Htmp = -3.0;
        for (int j = 0; j < 600; ++j)
        {
            Htmp += 0.01;
            Xtmp = arma::sum(arma::pow(y - atmp * arma::pow(x, Htmp), 2));
            //--find minimum of X and therefor values of a_best & H_best
            if (Xtmp < X)
            {
                X = Xtmp;
                a = atmp;
                H = Htmp;
            }
        }
    }
    cout << "t(sph) = " << double(clock() - t) / CLOCKS_PER_SEC << " sec" << endl;

    //--save a_best & H_best in sph.txt
    ofstream sph;
    sph.open("sph.txt", ios::out);
    sph << a << "\t" << H << endl;
    sph.close();
    //--reurn yt = ax^H
    return arma::conv_to<stdvect>::from(a * arma::pow(x, H));
}
