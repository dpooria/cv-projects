/*!
    Exercise set 14: please first take a look at 'readme14.pdf'. you can build and run this program by running 'main14.sh'
*/
#include <armadillo>
#include <functional>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace arma;

///-- x and dx and N defined globaly so every function could use'em
// const double dx = 0.25;
const double dx = 0.025;
///-- 0 <= x <= 1
const vec x = regspace(0, dx, 1);
const int N = x.size();

const vec u(double xl, double xu, int i);
const double intui(int i);
const double intuiuj(int i, int j);
const double intduiduj(int i, int j);
const vec FEM(std::function<double(int, int)> integralL, std::function<double(int)> integralfu);

int main(void)
{
    ///-- eqution: Ly = f(x)
    ///-- in both parts B and C : L = -d^2/dx^2 + 1
    ///-- according to part A : integralL(i, j) =
    ///-- integral_{xmin}^{xmax}(L{sum_{ij}(al_i * u_i (x) * u_j(x))}) =
    ///--          integral{du_i/dx * du_j/dx dx} + integral{u_i(x) * u_j(x)dx}
    const auto integralL = [](int i, int j) {
        return intduiduj(i, j) + intuiuj(i, j);
    };
    ///=== part B: f(x) = 1 ->
    ///-- integralfu(i) = integral_{xmin}^{xmax}(f(x) * u_i(x) dx) = integral(u_i dx)
    const vec y1 = FEM(integralL, intui);
    ///-- anaylitical solution for part B
    const vec y1t = (-1.0 / (1.0 + exp(1))) * (exp(-x + 1) + exp(x)) + 1.0;
    ///=== part C:
    ///-- f(x) = sin(x) -> integralfu(i) = integral_{xmin}^{xmax}(f(x) * u_i(x) dx) = sum(integral_{xl}^{xu}(f(x) * u_i(x) dx)})
    ///-- integral_{xl}^{xu}(f(x) * u_i(x) dx)} = integral{sin(x) * a_i * (x - b_i)/h dx} = (a_i/h)*[integral{xsin(x)dx} - b_i*integral{sin(x)}
    ///-- = (a_i/h) * [-xcos(x) + sin(x) + b_i * cos(x)]|_{xu}^{xl}
    const auto integralFU2 = [](int i) -> double {
        double s = 0.0;
        for (int k = 1; k < N; ++k)
        {
            double xL = x[k - 1];
            double xU = x[k];
            vec r = u(xL, xU, i);
            s += r[0] * (-xU * cos(xU) + sin(xU) + r[1] * cos(xU));
            s -= r[0] * (-xL * cos(xL) + sin(xL) + r[1] * cos(xL));
        }
        return s / dx;
    };
    const vec y2 = FEM(integralL, integralFU2);
    ///-- anaylitical solution for part C
    const vec y2t = (-0.5 * sin(1.0) / sinh(1.0)) * sinh(x) + 0.5 * sin(x);
    //=========================================================================
    ///-- save and plot the results
    std::ofstream fres1("FEMResultF1.txt");
    std::ofstream fres2("FEMResultFsin.txt");
    if (!fres1 || !fres2)
    {
        std::cout << "error opening file!" << std::endl;
        return -1;
    }
    fres1.setf(std::ios_base::scientific);
    fres2.setf(std::ios_base::scientific);

    for (int i = 0; i < N; ++i)
    {
        fres1 << x[i] << "\t" << y1[i] << "\t" << y1t[i] << "\n";
        fres2 << x[i] << "\t" << y2[i] << "\t" << y2t[i] << "\n";
    }
    fres1.close();
    fres2.close();
    const char *gnuplot = R"END(gnuplot -e "
    set terminal png enhanced; set grid;
    set xlabel 'x'; set ylabel '{/Symbol F}(x)';
    set output 'FEMResultF1.png'; set key bottom center;
    set title 'equation:  -{/Symbol F}\"(x) + {/Symbol F}(x) = 1';
    p 'FEMResultF1.txt' u 1:2 w p pt 6 ps 1.0 lc 'red' t 'FEM results',
        '' u 1:3 w l lc 'blue' t 'Analytical solution';
    set output 'FEMResultFsin.png';
    set title 'equation:  -{/Symbol F}\"(x) + {/Symbol F}(x) = sin(x)';
    p 'FEMResultFsin.txt' u 1:2 w p pt 6 ps 1.0 lc 'red' t 'FEM results',
        '' u 1:3 w l lc 'blue' t 'Analytical solution';"
    )END";
    if (system(gnuplot))
    {
        std::cout << "plot error!" << std::endl;
        std::cout << gnuplot << std::endl;
    }

    return 0;
}
//==========================================================================
///--   u_i(x) = (x_{i+1} - x)/dx if x_i <= x <= x_{i+1}
///--   u_i(x) = (x - x_{i-1})/dx if x_{i-1} <= x <= x_i
///--   else  u_i(x) = 0
///-- @function u -> in all the functions instead of above notation I use:u_i(x) = a(x - b)/dx
///--                 a = -1 & b = x_{i+1} if x_i <= x <= x_{i+1}
///--                 a = 1 & b = x_{i-1} if x_{i-1} <= x <= x_i
///-- therefore this function instead of taking one point x takes an interval [xl, xu]
///-- and returns a and b -> a = -1 & b = x_{i+1} if x_i <= xl && xu <= x_{i+1}
///--                        a = 1 & b = x_{i-1} if x_{i-1} <= xl && xu <= x_i
///-- this form is more suitable for integration as you can see in 'readme14.pdf'
const vec u(double xl, double xu, int i)
{
    double a, b;
    if (xl >= x[i] && xu <= x[i + 1])
    {
        a = -1;
        b = x[i + 1];
    }
    else if (xl >= x[i - 1] && xu <= x[i])
    {
        a = 1;
        b = x[i - 1];
    }
    else
    {
        a = 0;
        b = 0;
    }
    return {a, b};
}
//================================================================================
///-- @function intui -> returns integral_{0}^{1}{u_i(x)dx}
///-- u_i(x) = a(x-b)/h => int_{xl}^{xu}{u_i(x) dx} = a(xu^2/2 - b*xu)/h - a(xl^2/2 - b*xl)/h
///-- integral_{0}^{1}{u_i(x)dx} = sum(int_{xl}^{xu}{u_i(x) dx}) for all x's
const double intui(int i)
{
    double s = 0.0;
    for (int k = 1; k < N; ++k)
    {
        double xU = x[k];
        double xL = x[k - 1];
        vec r = u(xL, xU, i);
        s += r[0] * ((xU * xU - xL * xL) / 2.0 - r[1] * (xU - xL));
    }
    return s / dx;
}
//=================================================================================
///-- @function intuiuj -> returns integral_{0}^{1}{u_i(x)u_j(x)dx}
///-- u_i * u_j = a_i*a_j*(x^2 - (b_i + b_j)x + b_i * b_j)
///-- int_{xl}^{xu}{u_i*u_j dx} = a_i*a_j*(x^3/3 - (b_i+b_j)x^2/2 + b_i*b_j*x)|_{xl}^{xu}
///-- integral_{0}^{1}{u_i(x)u_j(x)dx} = sum(int_{xl}^{xu}{u_i*u_j dx}) for all x's
const double intuiuj(int i, int j)
{
    double s = 0.0;
    const auto ll = [](double xx, double a, double b, double c) -> double {
        return a * (pow(xx, 3) / 3.0 - b * xx * xx / 2.0 + c * xx);
    };
    for (int k = 1; k < N; ++k)
    {
        double xU = x[k];
        double xL = x[k - 1];
        vec ri = u(xL, xU, i);
        vec rj = u(xL, xU, j);
        s += ll(xU, ri[0] * rj[0], ri[1] + rj[1], ri[1] * rj[1]);
        s -= ll(xL, ri[0] * rj[0], ri[1] + rj[1], ri[1] * rj[1]);
    }
    return s / (dx * dx);
}
//==================================================================================
///-- @function intduiduj -> returns integral_{xmin}^{xmax}(du_i/dx*du_j/dx dx)
///-- du_i/dx = a_i/dx => du_i/dx * du_j/dx = a_i * a_j /dx^2 ->
///-- int_{xl}^{xu}{du_i/dx * du_j/dx dx} = a_i * a_j * (xu - xl)/dx^2
///-- integral_{xmin}^{xmax}(du_i/dx*du_j/dx dx) = sum(int_{xl}^{xu}{du_i/dx * du_j/dx dx}) for all x's
const double intduiduj(int i, int j)
{
    double s = 0.0;
    for (int k = 1; k < N; ++k)
    {
        double xU = x[k];
        double xL = x[k - 1];
        vec ri = u(xL, xU, i);
        vec rj = u(xL, xU, j);
        s += ri[0] * rj[0] * (xU - xL);
    }
    return s / (dx * dx);
}
//===================================================================================
///-- @function FEM -> returns the solution (y = sum_i(al_i * u_i(x))) according to given functions and with all boundary conditions = 0
///-- for hypothetical equation  Ly = f(x) (where L is an operator) ->
///-- integralL(i, j) = integral_{xmin}^{xmax}(L{sum_{ij}(al_i * u_i (x) * u_j(x))})
///-- integralfu(i) = integral_{xmin}^{xmax}(f(x) * u_i(x))
const vec FEM(std::function<double(int, const int)> integralL,
              std::function<double(int)> integralfu)
{
    vec F = zeros<vec>(N - 2);
    for (int i = 0; i < N - 2; ++i)
    {
        F[i] = integralfu(i + 1);
    }
    mat A = zeros<mat>(N - 2, N - 2);
    for (int i = 0; i < N - 2; ++i)
    {
        for (int j = 0; j < N - 2; ++j)
        {
            A(i, j) = integralL(i + 1, j + 1);
        }
    }
    const vec al = solve(A, F);
    vec y = zeros<vec>(N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N - 2; ++j)
        {
            vec r = u(x[i], x[i], j + 1);
            y[i] += al[j] * r[0] * (x[i] - r[1]) / dx;
        }
    }
    return y;
}
