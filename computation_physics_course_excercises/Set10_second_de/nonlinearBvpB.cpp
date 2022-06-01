/* Exercise set 10 problem 3 part B  
    explained in 'readme10.pdf'
*/
#include <armadillo>
#include <iostream>
#include <cmath>

using namespace arma;

const double f(double y, double py, double qy)
{
    return -(y * qy) + (py * py) - 1.0;
}

const double g(double y, double py, double qy, double z, double pz, double qz)
{
    return -(qy * z) + (2.0 * py * pz) - y * qz;
}

int main()
{
    const int N{1000};
    const double a = 0.0, b = 1.0;
    const vec t = linspace<vec>(a, b, N);
    const double dt = t[1] - t[0];
    const double alpha{0.0}, beta{0.0};
    const double alphap{0.0};
    double T = (beta - alpha) / (b - a);
    vec y = zeros<vec>(N);
    vec z = zeros<vec>(N);
    double py, qy;
    double pz, qz;
    double prevPy, prevQy;
    double prevPz, prevQz;
    double Tol = 1000.0;
    const double eps{1e-04};
    while (Tol > eps)
    {
        y[0] = alpha;
        z[0] = 0.0;
        prevPz = 0.0;
        prevPy = alphap;
        prevQz = 1.0;
        prevQy = T;
        for (int i = 0; i < N - 1; ++i)
        {
            //-- solve equation y''' = fy with Euler Method
            y[i + 1] = y[i] + prevPy * dt;
            py = prevPy + dt * prevQy;
            qy = prevQy + dt * f(y[i], prevPy, prevQy);
            //-- solve equation z''' = gz with Euler Method
            z[i + 1] = z[i] + prevPz * dt;
            pz = prevPz + dt * prevQz;
            qz = prevQz + dt * g(y[i], prevPy, prevQy, z[i], prevPz, prevQz);
            //-- update
            prevPy = py;
            prevPz = pz;
            prevQy = qy;
            prevQz = qz;
        }
        T += (beta - y[N - 1]) / z[N - 1];
        Tol = std::abs(beta - y[N - 1]);
    }
    std::ofstream fres;
    fres.open("nonlinearBvpB-res.txt");
    for(int i = 0; i < N; ++i)
    {
        fres << t[i] << "\t" << y[i] << "\n";
    }
    fres.close();
    char gnuplotB[] = R"END( gnuplot -e "
    set terminal png; set grid; unset key;
    set style line 1 lc rgb 'blue' pt 6 ps 0.3;
    set xlabel 't'; set ylabel 'y(t)';
    set output 'nonlinearBvpB-res.png';
    plot 'nonlinearBvpB-res.txt' w p ls 1;"
    )END";
    if(system(gnuplotB))
    {
        std::cout << "plot error!" << std::endl;
        std::cout << gnuplotB << std::endl;
    }
    return 0;
}