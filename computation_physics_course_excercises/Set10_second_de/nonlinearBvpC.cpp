/* Exercise set 10 problem 3 part C 
    explained in 'readme10.pdf'
*/

#include <armadillo>
#include <iostream>
#include <cmath>

using namespace arma;

const double f(double t, double y)
{
    double f = -std::pow(y, 2.0);
    f += std::pow(t, -2.5) * (9.0 + 30 * t + 105 * t * t) / 16.0;
    f += pow(t, 3.0) * pow(1.0 - t, 4.0);
    return f;
}

const double g(double y, double z)
{
    return -2 * y * z;
}

int main()
{
    const int N{100};
    const double a = 0.01, b = 1.0;
    const double eps{1e-04};
    const vec t = linspace<vec>(a, b, N);
    const double dt = t[1] - t[0];
    const double alpha{0.0}, beta{0.0};
    const double alphap{0.0}, betap{0.0};
    double T = (beta - alpha) / (b - a);
    double Tp = (betap - alphap) / (b - a);
    double Tol{100.0}, Tolp{100.0};
    double py, qy, ry;
    double prevPy, prevQy, prevRy;
    double pz, qz, rz;
    double prevPz, prevQz, prevRz;
    vec y = zeros<vec>(N);
    vec z = zeros<vec>(N);
    while (Tol > eps || Tolp > eps)
    {
        y[0] = alpha;
        z[0] = 0.0;
        prevPz = 0.0;
        prevPy = alphap;
        prevQz = 1.0;
        prevQy = T;
        prevRz = 0.0;
        prevRy = Tp;
        for (int i = 0; i < N - 1; ++i)
        {
            //-- solve equation y'''' = f(t,y) with Euler method
            y[i + 1] = y[i] + dt * prevPy;
            py = prevPy + dt * prevQy;
            qy = prevQy + dt * prevRy;
            ry = prevRy + dt * f(t[i], y[i]);
            //-- solve equation z'''' = g(y,z) with Euler method
            z[i + 1] = z[i] + dt * prevPz;
            pz = prevPz + dt * prevQz;
            qz = prevQz + dt * prevRz;
            rz = prevRz + dt * g(y[i], z[i]);
            //-- update
            prevPy = py;
            prevPz = pz;
            prevQy = qy;
            prevQz = qz;
            prevRz = rz;
            prevRy = ry;
        }
        T += (beta - y[N - 1]) / z[N - 1];
        Tp += (betap - py) / pz;
        Tol = std::abs(beta - y[N - 1]);
        Tolp = std::abs(betap - py);
    }
    //-- exact solution
    vec yE = pow(t, 1.5) % pow(1.0 - t, 2.0);
    //-- save results
    std::ofstream fres;
    fres.open("nonlinearBvpC-res.txt");
    for (int i = 0; i < N; ++i)
    {
        fres << t[i] << "\t" << y[i] << "\t" << yE[i] << "\n";
    }
    fres.close();
    //-- plot'em
    char gnuplotC[] = R"END( gnuplot -e "
    set terminal png; set grid;
    set style line 1 lc rgb 'blue';
    set style line 2 lc rgb 'red' pt 6 ps 0.5;
    set xlabel 't'; set ylabel 'y(t)';
    set output 'nonlinearBvpC-res.png';
    plot 'nonlinearBvpC-res.txt' u 1:3 w l ls 1 t 'exact solution',
    '' u 1:2 w p ls 2 t 'numerical results';"
    )END";
    if (system(gnuplotC))
    {
        std::cout << "plot error!" << std::endl;
        std::cout << gnuplotC << std::endl;
    }
    return 0;
}