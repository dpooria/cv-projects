/* Set9_2 refer to readme9.pdf */
#include <armadillo>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

int main(int argc, char *argv[])
{
    clock_t et = clock();
    const int N{10000000};
    //-- b -> infty
    const double b{1e04}, beta{0.2}, m{0.01};
    const double al = beta * m / 2.0, A = pow(beta * m / (2.0 * M_PI), 1.5);
    auto f = [al, A](double rr) -> double {
        return (4.0 / 3.0) * M_PI * A * powl(rr, 4.0) * exp(-al * powl(rr, 2.0));
    };
    //-- findout the maximum of the function
    arma::vec r = arma::linspace(0, b, N);
    arma::vec F = (4.0 / 3.0) * M_PI * A * arma::pow(r, 4.0) % arma::exp(-al * arma::pow(r, 2.0));
    const double fMax = F.max();
    srand((unsigned)time(0));
    unsigned int Ns{0};
    double x, y;
    for (int i = 0; i < N; ++i)
    {
        //-- 0 <= x <= b
        x = (double)rand() / RAND_MAX * b;
        //-- 0 <= y <= fMax
        y = (double)rand() / RAND_MAX * fMax;
        if (y <= f(x))
        {
            Ns++;
        }
    }
    //-- fMax * b -> area of the rectangular
    double ans = ((Ns + 0.0) / N) * (fMax * b);
    std::ofstream res;
    res.open("IntegralResult.txt");
    std::cout << "for beta = " << beta << " and m = " << m << std::endl;
    std::cout << "<vz^2> =\n\tMonte Carlo Approximation: " << ans << std::endl;
    std::cout << "\tTheory: " << 0.5 * A * pow(M_PI, 1.5) * pow(al, -2.5) << std::endl;
    res << "for beta = " << beta << " and m = " << m << std::endl;
    res << "<vz^2> =\nMonte Carlo Approximation: " << ans << std::endl;
    res << "Theory: " << 0.5 * A * pow(M_PI, 1.5) * pow(al, -2.5) << std::endl;
    res.close();
    std::cout << "-Program 'Integral' finished in "
              << (double)(clock() - et) / CLOCKS_PER_SEC << " sec" << std::endl;
    return 0;
}