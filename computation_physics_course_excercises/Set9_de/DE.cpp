/* Set9_3 and Set9_4 explained in readme.pdf*/
#include <armadillo>
#include <iostream>
#include <cmath>
#include <ctime>
#include <string>

int plot(std::string name1, std::string t1, std::string name2, std::string t2, std::string res)
{
    std::string gnuplot =
        "gnuplot -e \"set terminal png enhanced; set output '" +
        res + "';set style line 1 lc rgb 'red' pt 6 ps 0.8;\
        set style line 2 lc rgb 'blue';\
        set xlabel 't/dt'; set ylabel 'y(t)'; set grid; \
        p '" +
        name1 + "' w p ls 1 t '" +
        t1 + "','" + name2 + "' w l ls 2 t '" + t2 + "';\"";
    return system(gnuplot.c_str());
}

int main(int argc, char *argv[])
{
    clock_t et = clock();
    const double A{1.0}, w{0.5}, a = 2.0 * w;
    const int N{1000};
    auto f = [A, w, a](double tt, double yy, double pp) -> double {
        return cos(w * tt) - a * pp - pow(w, 2.0) * yy;
    };
    arma::vec t = arma::linspace(0.0, 10.0, N);
    const double dt = t[1] - t[0];
    arma::vec yE(N), yRF4(N), yFD(N);
    yE[0] = yRF4[0] = yFD[0] = A;
    double p = 0.0;
    //-- Euler Method
    for (int i = 0; i < N - 1; ++i)
    {
        yE[i + 1] = yE[i] + p * dt;
        p += f(t[i], yE[i], p);
    }
    //-- RF4 Method
    double k1, k2, k3, k4, m1, m2, m3, m4;
    p = 0.0;
    for (int i = 0; i < N - 1; ++i)
    {
        k1 = p;
        m1 = f(t[i], yRF4[i], p);
        k2 = p + m1 / 2.0;
        m2 = f(t[i] + dt / 2.0, yRF4[i] + k1 / 2.0, k2);
        k3 = p + m2 / 2.0;
        m3 = f(t[i] + dt / 2.0, yRF4[i] + k2 / 2.0, k3);
        k4 = p + m3;
        m4 = f(t[i] + dt, yRF4[i] + k3, k4);
        yRF4[i + 1] = yRF4[i] + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        p += (dt / 6.0) * (m1 + 2 * m2 + 2 * m3 + m4);
    }

    //-- Theory
    arma::vec yt = (A + ((1.0 / a) - A) * t) % arma::exp(-a * t / 2.0) + arma::sin(w * t) / (a * w);

    //-- Finite Difference Method
    yFD[N - 1] = yt[N - 1];
    arma::vec C(N - 2);
    arma::mat D = arma::zeros<arma::mat>(N - 2, N - 2);
    D.diag(0) = 2.0 + pow(dt, 2.0) * -pow(w, 2.0) * arma::ones<arma::vec>(N - 2);
    D.diag(1) = ((dt / 2.0) * -a - 1.0) * arma::ones<arma::vec>(N - 3);
    D.diag(-1) = ((dt / 2.0) * a - 1.0) * arma::ones<arma::vec>(N - 3);
    C = -pow(dt, 2.0) * arma::cos(w * t.rows(1, N - 2));
    C[0] += (-a * dt / 2.0 + 1.0) * yFD[0];
    C[N - 3] += (a * dt / 2.0 + 1.0) * yFD[N - 1];
    yFD.rows(1, N - 2) = arma::solve(D, C);
    //-- save results
    yFD.save("DE-FiniteDifference.txt", arma::raw_ascii);
    yE.save("DE-Euler.txt", arma::raw_ascii);
    yRF4.save("DE-RF4.txt", arma::raw_ascii);
    yt.save("DE-Theory.txt", arma::raw_ascii);
    int respond = plot("DE-Theory.txt", "Theory", "DE-Euler.txt", "Euler", "DE-Euler.png");
    respond += plot("DE-Theory.txt", "Theory", "DE-RF4.txt", "RF4", "DE-RF4.png");
    respond += plot("DE-Theory.txt", "Theory", "DE-FiniteDifference.txt",
                    "Finite Difference", "DE-FiniteDifference.png");
    respond += plot("DE-Euler.txt", "Euler", "DE-FiniteDifference.txt",
                    "Finite Difference", "DE-EulerVsFD.png");
    respond += plot("DE-RF4.txt", "RF4", "DE-FiniteDifference.txt",
                    "Finite Difference", "DE-RF4VsFD.png");

    if (respond)
    {
        std::cerr << "plot error!" << std::endl;
    }
    std::cout << "-Program 'DE' finished in "
              << (double)(clock() - et) / CLOCKS_PER_SEC << " sec" << std::endl;

    return 0;
}