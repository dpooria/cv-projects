/* Exercise set 10 problem 3 part A  
    explained in 'readme10.pdf'
*/
#include <armadillo>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace arma;

int main()
{

    const int N{100};
    const double eps{1e-03};
    vec t = arma::linspace(1.0, 2.0, N);
    const double dt = t[1] - t[0];
    const double alpha = 2.0, beta = 2.5;
    vec y = zeros<vec>(N);
    auto f = [](double tt, double yy) -> double {
        return 2.0 * pow(yy, 3.0) - 6.0 * yy - 2.0 * pow(tt, 3.0);
    };
    auto dfdyz = [](double yy, double zz) -> double {
        return 6.0 * (pow(yy, 2.0) - 1.0) * zz;
    };
    double Tol = 1000.0;
    double k1, k2, k3, k4, m1, m2, m3, m4;
    double kz1, kz2, kz3, kz4, mz1, mz2, mz3, mz4;
    double T = (beta - alpha) / N; //-- T = (beta-alpha)/(t[N-1]-t[0]) -> doesn't converge!
    double p, pz;
    while (Tol > eps)
    {
        y[0] = alpha;
        vec z = zeros<vec>(N);
        z[0] = 0.0;
        p = T;
        pz = 1.0;
        for (int i = 0; i < N - 1; ++i)
        {
            //-- solve y''=f(t,y) equation with RKF45 method
            k1 = p;
            m1 = f(t[i], y[i]);
            k2 = p + dt * m1 / 2.0;
            m2 = f(t[i] + dt / 2.0, y[i] + dt * k1 / 2.0);
            k3 = p + dt * m2 / 2.0;
            m3 = f(t[i] + dt / 2.0, y[i] + dt * k2 / 2.0);
            k4 = p + dt * m3;
            m4 = f(t[i] + dt, y[i] + dt * k3);
            y[i + 1] = y[i] + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            p += (dt / 6.0) * (m1 + 2.0 * m2 + 2.0 * m3 + m4);
            //-- solve z''=dfdyz(y,z) equation with RKF45 method
            kz1 = pz;
            mz1 = dfdyz(y[i], z[i]);
            kz2 = pz + dt * mz1 / 2.0;
            mz2 = dfdyz(y[i] + dt * k1 / 2.0, z[i] + dt * kz1 / 2.0);
            kz3 = pz + dt * mz2 / 2.0;
            mz3 = dfdyz(y[i] + dt * k2 / 2.0, z[i] + dt * kz2 / 2.0);
            kz4 = pz + dt * mz3;
            mz4 = dfdyz(y[i] + dt * k3, z[i] + dt * kz3);
            z[i + 1] = z[i] + (dt / 6.0) * (kz1 + 2.0 * kz2 + 2.0 * kz3 + kz4);
            pz += (dt / 6.0) * (mz1 + 2.0 * mz2 + 2.0 * mz3 + mz4);
        }
        T += (beta - y[N - 1]) / z[N - 1];
        Tol = std::abs(beta - y[N - 1]);
    }
    vec yE = t + 1.0 / t;
    std::ofstream fres;
    fres.open("nonlinearBvpA-res.txt");
    for (int i = 0; i < N; ++i)
    {
        fres << t[i] << "\t" << y[i] << "\t" << yE[i] << "\n";
    }
    fres.close();
    char gnuplotA[] = R"END( gnuplot -e "
    set terminal png enhanced; set grid;
    set style line 1 lc rgb 'blue';
    set style line 2 lc rgb 'red' pt 1 ps 0.3;
    set xlabel 't'; set ylabel 'y(t)';
    set output 'nonlinearBvpA-res.png';
    p 'nonlinearBvpA-res.txt' u 1:3 w l ls 1 t 'exact solution',
        '' u 1:2 w p ls 2 t 'numerical results';"
    )END";
    if (system(gnuplotA))
    {
        std::cout << "plot error!" << std::endl;
        std::cout << gnuplotA << std::endl;
    }
    return 0;
}