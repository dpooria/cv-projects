/* Exercise set 8 part 1 */
#include <armadillo>
#include <iostream>
#include <cstring>

using namespace arma;

vec cd3(const vec &x, const vec &f);
vec cd5(const vec &x, const vec &f);
vec cd7(const vec &x, const vec &f);
vec cd9(const vec &x, const vec &f);

int main(int argc, char *argv[])
{
    mat data;
    if (!data.load("dataprofile.txt"))
    {
        std::cerr << "dataprofile.txt not found\n";
        return -1;
    }
    vec x = data.col(0);
    vec f = data.col(1);
    vec cdf3 = cd3(x, f);
    cdf3.save("cdf3.txt", arma::raw_ascii);
    vec cdf5 = cd5(x, f);
    cdf3.save("cdf5.txt", arma::raw_ascii);
    vec cdf7 = cd7(x, f);
    cdf7.save("cdf7.txt", arma::raw_ascii);
    vec cdf9 = cd9(x, f);
    cdf9.save("cdf9.txt", arma::raw_ascii);
    char gnuplot[] =
        "gnuplot -e \"set terminal jpeg enhanced;set xlabel 'x/dx';set ylabel 'df/dx';set grid;\
     p 'cdf3.txt' w p t '3-Pts','cdf5.txt' w l t '5-Pts','cdf7.txt' w p t '7-Pts','cdf9.txt' \
     w l t '9-Pts';\" > 'cdf.jpg'";
    if (system(gnuplot))
    {
        cerr << "Plot Error!\n";
    }
    return 0;
}
//-- 3-Point central differential
vec cd3(const vec &x, const vec &f)
{
    int N = x.size();
    double dx = x[1] - x[0];
    vec cdf(N, fill::zeros);
    for (int i = 1; i < N - 1; ++i)
    {
        cdf(i) = (-f(i - 1) + f(i + 1)) / (2 * dx);
    }
    return cdf;
}

//-- 5-Point central differential
vec cd5(const vec &x, const vec &f)
{
    int N = x.size();
    double dx = x[1] - x[0];
    vec cdf(N, fill::zeros);
    for (int i = 2; i < N - 2; ++i)
    {
        cdf(i) = (f(i - 2) - 8 * f(i - 1) + 8 * f(i + 1) - f(i + 2)) / (12 * dx);
    }
    return cdf;
}
//-- 7-Point central differential(refer to readme.pdf)
vec cd7(const vec &x, const vec &f)
{
    int N = x.size();
    double dx = x[1] - x[0];
    vec cdf(N, fill::zeros);
    for (int i = 3; i < N - 3; ++i)
    {
        cdf(i) = (-f(i - 3) + 9 * f(i - 2) - 45 * f(i - 1) +
                  45 * f(i + 1) - 9 * f(i - 2) + f(i + 3)) /
                 (60 * dx);
    }
    return cdf;
}

//-- 9-Point central differential(refer to readme.pdf)
vec cd9(const vec &x, const vec &f)
{
    int N = x.size();
    double dx = x[1] - x[0];
    vec cdf(N, fill::none);
    for (int i = 4; i < N - 4; ++i)
    {
        cdf(i) = (3 * f(i - 4) - 32 * f(i - 3) + 168 * f(i - 2) - 672 * f(i - 1) +
                  672 * f(i + 1) - 168 * f(i - 2) + 32 * f(i + 3) - 3 * f(i + 4)) /
                 (840 * dx);
    }
    return cdf;
}