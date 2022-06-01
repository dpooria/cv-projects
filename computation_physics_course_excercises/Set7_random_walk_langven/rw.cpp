#include <armadillo>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

int plot(const string fname, string title, int n)
{
    string gnuplot = "gnuplot -e \"set terminal jpeg enhanced; set xlabel 't';\
    set style line 2 lc rgb 'red' lt 1 lw 0.5 pt 7 ps 1; set grid;\
    set title '" +
                     title + "';\
                     p '" +
                     fname + ".txt' u 1:" + to_string(n) + " w lp ls 2;\" >'" +
                     title + ".jpg'";
    return system(gnuplot.c_str());
}

//-- random walk
void RW(arma::vec (*arand)(int), const string name)
{
    //-- N => time, M => ensembl
    const int N{100}, M{100};
    const double dt{0.1};
    double t{0.0};
    arma::vec xmean(N), var(N);
    xmean = var = arma::zeros<arma::vec>(N);
    //-- x0 => initial position
    arma::vec x0 = arma::zeros<arma::vec>(M);
    ofstream res;
    res.open(name + ".txt");
    //-- loop through time
    for (int i = 0; i < N; ++i)
    {
        //-- get steps randomly according to passed function
        arma::vec s = arand(M);
        arma::vec x(M);
        //-- step forward\backward according to s and previous position in each ensembl
        x = x0 + s;
        //-- save current position in x0
        x0 = x;
        xmean(i) = arma::mean<arma::vec>(x);
        //-- calculate variance
        for (int j = 0; j < M; ++j)
        {
            var(i) += powl(x(j) - xmean(i), 2.0);
        }
        var(i) /= M;
        res << t << "\t" << xmean(i) << "\t" << var(i) << "\n";
        t += dt;
    }
    res.close();
}

int main(int argc, char *argv[])
{
    srand((unsigned)time(0));
    //-- uniform random data in range [0,1] with armadillo function randu
    auto lrandU1 = [](int M) -> arma::vec {
        return arma::randu<arma::vec>(M);
    };
    //-- uniform random data in range [-1,1] with armadillo function randu
    auto lrandU2 = [](int M) -> arma::vec {
        return 2 * arma::randu<arma::vec>(M) - 1;
    };
    //-- uniform random data in range [-1,1] with c function rand()
    auto lrandU3 = [](int M) -> arma::vec {
        arma::vec x(M);
        for (int i = 0; i < M; ++i)
        {
            x(i) = 2 * (rand() / (RAND_MAX + 0.0)) - 1.0;
        }
        return x;
    };
    //-- normal random data with armadillo function randn
    auto lrandG1 = [](int M) -> arma::vec {
        return arma::randn<arma::vec>(M);
    };
    //-- normal random data with Box-Muller method
    auto lrandG2 = [](int M) -> arma::vec {
        double r1, r2;
        arma::vec x(M);
        for (int i = 0; i < M; ++i)
        {
            r1 = rand() / (RAND_MAX + 0.0);
            r2 = rand() / (RAND_MAX + 0.0);
            x(i) = sqrt(-2 * log(1 - r1)) * cos(2 * M_PI * r2);
        }
        return x;
    };
    RW(lrandU1, "RWuniform1");
    RW(lrandU2, "RWuniform2");
    RW(lrandU3, "RWuniform3");
    RW(lrandG1, "RWgaussian1");
    RW(lrandG2, "RWgaussian2");
    plot("RWuniform1", "rwUmean1", 2);
    plot("RWuniform1", "rwUVar1", 3);
    plot("RWuniform2", "rwUmean2", 2);
    plot("RWuniform2", "rwUVar2", 3);
    plot("RWuniform3", "rwUmean3", 2);
    plot("RWuniform3", "rwUVar3", 3);
    plot("RWgaussian1", "rwGmean1", 2);
    plot("RWgaussian1", "rwGVar1", 3);
    plot("RWgaussian2", "rwGmean2", 2);
    plot("RWgaussian2", "rwGVar2", 3);
    return 0;
}
