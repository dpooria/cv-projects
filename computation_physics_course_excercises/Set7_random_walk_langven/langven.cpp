#include <armadillo>
#include <cmath>
#include <ctime>

using namespace std;

typedef arma::vec avec;

avec pdf(const avec &x)
{
    int N = x.size();
    double xmax = x.max(), xmin = x.min();
    double dx = (xmax - xmin) / N;
    int kmin = static_cast<int>(xmin / dx), kmax = static_cast<int>(xmax / dx);
    avec p = arma::zeros<avec>(N);
    for (int i = 0; i < N; ++i)
    {
        int k = static_cast<int>(x(i) / dx);
        p(k - kmin)++;
    }
    return p / (N * dx);
}

arma::mat jointpdf(const avec &x)
{
    int N = x.size();
    double xmin = x.min(), xmax = x.max();
    double dx = (xmax - xmin) / N;
    int kmin = static_cast<int>(xmin / dx);
    int kmax = static_cast<int>(xmax / dx);
    arma::mat p = arma::zeros<arma::mat>(N, N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N - i; ++j)
        {
            int k1 = static_cast<int>(x(j) / dx) - kmin;
            int k2 = static_cast<int>(x(j + i) / dx) - kmin;
            p(k1, k2)++;
        }
    }
    return p / (N * dx);
}

avec correlation(const avec &x)
{
    int N = x.size();
    avec c = arma::zeros<avec>(N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N - i; ++j)
        {
            c(i) += x(i) * x(i + j);
        }
        c(i) /= (N - i);
    }
    return c;
}

int plot(string name, string ylabel)
{
    string gnuplot = "gnuplot -e \"set terminal jpeg enhanced; set xlabel '{/Symbol t}/dt';\
    set ylabel '" + ylabel +
                     "'; set style line 1 lc rgb 'red' pt 5 ps 0.8;\
    set style line 2 lc rgb 'blue' pt 6 ps 0.8;set grid;p '" +
                     name + ".txt' w p ls 1 t 'Simulation Results','" +
                     name + "T.txt' w p ls 2 t 'expected';\">'" + name + ".jpg'";
    return system(gnuplot.c_str());
}

int main(int argc, char *argv[])
{
    const int N{1000};
    const double g{0.5}, k{0.1}, dt{0.01};
    avec eta(N), x(N), v(N);
    x(0) = v(0) = 0.0;
    srand((unsigned)time(0));
    for (int i = 1; i < N; ++i)
    {
        double r1 = rand() / (RAND_MAX + 0.0);
        double r2 = rand() / (RAND_MAX + 0.0);
        //-- generate eta using Box-Muller method
        eta(i) = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * 2 * g;
        v(i) = v(i - 1) - k * v(i - 1) * dt + eta(i) * sqrt(dt);
        x(i) = x(i - 1) + v(i) * dt;
    }
    cout << "\tSimulation Results\t"
         << "Expected Values\n";
    cout << "___________________________________________________\n";
    cout << "<v(t)> : \t" << arma::mean(v) << "\t\t" << 0.0 << "\n";
    cout << "<v^2(t)> :\t" << arma::mean(v % v)
         << "\t\t" << (g / k) * (1 - exp(-2 * k * 10)) << "\n";
    cout << "<x(t)> : \t" << arma::mean(x)
         << "\t\t" << 0.0 << "\n";
    cout << "<x^2(t)> :\t" << arma::mean(x % x)
         << "\t\t" << 2 * g / powl(k, 2) << "\n";
    avec tau = arma::linspace(0, 10, N);
    avec vcorrel = correlation(v);
    avec vcorrelT = (g / k) * arma::exp(-k * tau);
    vcorrel.save("vcorrel.txt", arma::raw_ascii);
    vcorrelT.save("vcorrelT.txt", arma::raw_ascii);
    plot("vcorrel", "<v(t)v(t+{/Symbol t})>");
    avec xcorrel = correlation(x);
    avec xcorrelT = vcorrelT * 50;
    xcorrel.save("xcorrel.txt", arma::raw_ascii);
    xcorrelT.save("xcorrelT.txt", arma::raw_ascii);
    plot("xcorrel", "<x(t)x(t+{/Symbol t})>");
    avec vpdf = pdf(v);
    avec vpdfT = arma::exp(-k * arma::pow(v, 2) / (2 * g));
    vpdf.save("vpdf.txt", arma::raw_ascii);
    vpdfT.save("vpdfT.txt", arma::raw_ascii);
    plot("vpdf", "p(v)");
    arma::mat vjointpdf = jointpdf(v);
    vjointpdf.save("vjointpdf.txt", arma::raw_ascii);
    avec Delta = arma::zeros<avec>(N);
    for (int t = 0; t < N; ++t)
    {
        for (int i = 0; i < N - t; ++i)
        {
            Delta(t) += abs(vjointpdf(i + t, i) - vpdf(i) * vpdf(i + t));
        }
    }
    Delta.save("DeltaJointPdf.txt", arma::raw_ascii);
    plot("DeltaJointPdf", "{/Symbol D}({/Symbol t})");
    return 0;
}
