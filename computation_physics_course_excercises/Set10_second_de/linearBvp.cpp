/* Exercise set10 problem 2 */

#include <armadillo>
#include <iostream>
#include <cmath>

using namespace arma;

int main()
{
    const int N = 100;
    const vec t = arma::linspace(0.0, 1.0, N);
    const double dt = t[1] - t[0];
    vec y = zeros<vec>(N);
    mat D = zeros<mat>(N - 2, N - 2);
    vec C = zeros<vec>(N - 2);
    //-- boundary conditions
    y[0] = 1.0;
    y[N - 1] = 3.0;
    //-- finite difference method
    C[0] = (dt + 1.0) * y[0];
    C[N - 3] = (dt + 1) * y[N - 1];
    D.diag(0) = (2.0 - pow(dt, 2.0)) * ones<vec>(N-2);
    D.diag(1) = (-dt - 1.0) * ones<vec>(N-3);
    D.diag(-1) = (dt - 1.0) * ones<vec>(N-3);
    y.rows(1,N-2) = solve(D,C);
    y.save("linearBvp-yFDM.txt", raw_ascii);
    //-- exact solution
    vec yt = ((3.0*exp(1.0) -1.0)*t + 1.0)%exp(-t);
    yt.save("linearBvp-yA.txt", raw_ascii);
    //-- plot results
    char gnuplot[] = R"end(
    gnuplot -e "set terminal png enhanced;
    set grid; set xlabel 't/dt'; set ylabel 'y(t)';
    set style line 1 lc rgb 'red' pt 6 ps 0.4;
    set style line 2 lc rgb 'blue';
    set output 'linearBvp.png';
    plot 'linearBvp-yA.txt' w l ls 2 t 'exact Solution',
    'linearBvp-yFDM.txt' w p ls 1 t 'FDM result';"
    )end";

    if(system(gnuplot))
    {
        std::cout << "error plot!" << std::endl;
        std::cout << gnuplot << std::endl;
    }
    return 0;
}
