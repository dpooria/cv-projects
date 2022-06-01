/* Set9 part 1 */
#include <armadillo>
#include <fstream>
#include <iostream>
#include <ctime>

using namespace arma;

int main(int argc, char *argv[])
{
    clock_t et = clock();
    const int K{197}, Nx{1000}, Nt{1000};
    const vec x = arma::linspace(0.0, 10.0, Nx);
    const double dx = x[1] - x[0];
    const vec t = arma::linspace(0.0, 10.0, Nt);
    const double dt = t[1] - t[0];
    mat TE = zeros<mat>(Nt, Nx);
    mat TIm = zeros<mat>(Nt, Nx);
    TE.row(0) = arma::exp(-x / 5.0).t();
    TIm.row(0) = arma::exp(-x / 5.0).t();
    TE.col(0) = 100.0 * arma::ones<vec>(Nt);
    TIm.col(0) = 100.0 * arma::ones<vec>(Nt);
    TE.col(Nx - 1) = 25.0 * arma::ones<vec>(Nt);
    TIm.col(Nx - 1) = 25.0 * arma::ones<vec>(Nt);
    double r = dt / (K * dx);
    //-- Explicit
    for (int n = 0; n < Nt - 1; ++n)
    {
        for (int j = 1; j < Nx - 1; ++j)
        {
            TE(n + 1, j) = (1 - 2 * r) * TE(n, j) + r * TE(n, j + 1) + r * TE(n, j - 1);
        }
    }
    TE.save("TP-Exp.csv", arma::csv_ascii);
    //-- Implicit
    mat A = zeros<mat>(3, 3);
    A.diag(0) = 1.0 + 2.0 * r * arma::ones<vec>(3);
    A.diag(-1) = -r * arma::ones<vec>(2);
    A.diag(1) = -r * arma::ones<vec>(2);
    double em2, ep2;
    for (int n = 0; n < Nt - 1; ++n)
    {
        for (int j = 2; j < Nx - 2; ++j)
        {
            em2 = r * TIm(n + 1, j - 2);
            ep2 = r * TIm(n + 1, j + 2);
            vec u = TIm(n, span(j - 1, j + 1)).t();
            u[0] += em2;
            u[2] += ep2;
            TIm(n + 1, span(j - 1, j + 1)) = solve(A, u).t();
        }
    }
    TIm.save("TP-Imp.csv", arma::csv_ascii);
    //-- Insulated endpoints
    mat TIns = zeros<mat>(Nt, Nx);
    TIns.row(0) = arma::exp(-x / 5.0).t();
    vec dT = -(1.0 / 5.0) * arma::exp(-x / 5.0);

    for (int n = 0; n < Nt - 1; ++n)
    {
        dT[0] = 0.0;
        dT[Nx - 1] = 0.0;
        for (int j = 0; j < Nx - 1; ++j)
        {
            TIns(n + 1, j) = TIns(n, j) + (K * dt / dx) * (dT(j + 1) - dT(j));
        }
        TIns(n + 1, Nx - 1) = 25.0;
        for (int j = 0; j < Nx - 2; ++j)
        {
            dT(j + 1) = (TIns(n + 1, j + 1) - TIns(n + 1, j)) / dx;
        }
    }
    TIns.save("TP-Ins.csv", arma::csv_ascii);
    //-- save and plot
    std::ofstream fileTEx;
    fileTEx.open("TP-Explicit.txt");
    std::ofstream fileTIm;
    fileTIm.open("TP-Implicit.txt");
    std::ofstream fileTIns;
    fileTIns.open("TP-InsulatorEnd.txt");

    for (int n = 0; n < Nt; ++n)
    {
        for (int j = 0; j < Nx; ++j)
        {
            fileTEx << n * dt << "\t" << j * dx << "\t" << TE(n, j) << "\n";
            fileTIm << n * dt << "\t" << j * dx << "\t" << TIm(n, j) << "\n";
            fileTIns << n * dt << "\t" << j * dx << "\t" << TIns(n, j) << "\n";
        }
    }
    fileTEx.close();
    fileTIm.close();
    fileTIns.close();
    char gnuplotE[] =
        "gnuplot -e \"set terminal png;\
        set hidden3d;set contour base;set dgrid3d;\
        set key off; set xlabel 't';\
        set ylabel 'x'; set zlabel 'ln(T(t,x))';set grid;\
        set output 'TP-Explicit.png';\
        splot 'TP-Explicit.txt' u 1:2:3 w l;\
        set output 'TP-InsulatorEnd.png';\
        splot 'TP-InsulatorEnd.txt' u 1:2:(\\$3 > 1e20  ? 1e20 : (\\$3 < -1e20 ? -1e20:\\$3)) w l;\
        set output 'TP-Implicit.png';\
        splot 'TP-Implicit.txt' u 1:2:3 w l;\"";

    if (std::system(gnuplotE))
    {
        std::cerr << "plot error!" << std::endl;
    }
    std::cout << "-Program 'TP' finished in "
              << (double)(clock() - et) / CLOCKS_PER_SEC << " sec" << std::endl;
    return 0;
}
