/*!
     Exercise set 13
     You can bulid and run the program with main13.sh 
*/
#include <armadillo>
#include <random>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

using namespace arma;

///-- generates an N x N  matrix(spin) with random element either 1(spin up) or -1(spin down)
Mat<int> GenerateS(const int N);
///-- sum over the nearest neighbors of the given point
int sumNN(const Mat<int> &s, int i, int j);
///-- Ising Hamiltonian in absence of external magnetic field
int HI(const Mat<int> &s);
///-- do the Ising simulation according to the external magnetic field
void IsingSim(const double Bext, const char *fname, const int N_atoms);
int main(void)
{
    IsingSim(0.0, "IsingBext0.txt", 400);
    IsingSim(1.0, "IsingBext1.txt", 400);
    ///-- plot the results
    ///-- file format of Ising function is like: T <H> <M> X C_v
    const char *gnuplot = R"END(gnuplot -e "
        set terminal png enhanced; set grid;
        set xlabel 'T'; unset key;
        set style line 1 lc rgb 'blue' pt 6 ps 1.0;
        set title 'B_{ext} = 0';
        set output 'IsingBext0Energy.png';
        set ylabel '<H>'; 
        p 'IsingBext0.txt' u 1:2 w lp ls 1;
        set output 'IsingBext0Magnetization.png';
        set ylabel '<M>'; 
        p 'IsingBext0.txt' u 1:3 w lp ls 1;
        set output 'IsingBext0X.png';
        set ylabel '{/Symbol K}'; 
        p 'IsingBext0.txt' u 1:4 w lp ls 1;
        set output 'IsingBext0C_v.png';
        set ylabel 'C_v'; 
        p 'IsingBext0.txt' u 1:5 w lp ls 1;set title 'B_{ext} = 0';
        set title 'B_{ext} = 1';
        set output 'IsingBext1Energy.png';
        set ylabel '<H>'; 
        p 'IsingBext1.txt' u 1:2 w lp ls 1;
        set output 'IsingBext1Magnetization.png';
        set ylabel '<M>'; 
        p 'IsingBext1.txt' u 1:3 w lp ls 1;
        set output 'IsingBext1X.png';
        set ylabel '{/Symbol K}'; 
        p 'IsingBext1.txt' u 1:4 w lp ls 1;
        set output 'IsingBext1C_v.png';
        set ylabel 'C_v'; 
        p 'IsingBext1.txt' u 1:5 w lp ls 1;"
    )END";
    if (system(gnuplot))
    {
        std::cout << "error plot!" << std::endl;
        std::cout << gnuplot << std::endl;
    }
    return 0;
}

///-- generates an N x N  matrix(spin) with random element either 1(spin up) or -1(spin down)
Mat<int> GenerateS(const int N)
{
    ///-- initialize the random generators
    std::random_device rd;
    std::default_random_engine g(rd());
    ///-- d -> returns either 0 or 1
    std::uniform_int_distribution<int> d(0, 1);
    Mat<int> s = zeros<Mat<int>>(N, N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            ///-- s(i,j) = 1 or -1
            s(i, j) = 2 * d(g) - 1;
        }
    }
    return s;
}
///-- sum over the nearest neighbors of the given point
int sumNN(const Mat<int> &s, int i, int j)
{
    const int n = s.n_rows;
    ///-- @param shorz -> sum of spin of the first horizontal neighbors
    int shorz = 0;
    ///-- @param svert -> sum of spin of the first vertical neighbors
    int svert = 0;
    ///-- shorz = ?
    if (i == 0)
    {
        shorz = s(i + 1, j);
    }
    else if (i == n - 1)
    {
        shorz = s(i - 1, j);
    }
    else
    {
        shorz = s(i - 1, j) + s(i + 1, j);
    }
    ///-- svert = ?
    if (j == 0)
    {
        svert = s(i, j + 1);
    }
    else if (j == n - 1)
    {
        svert = s(i, j - 1);
    }
    else
    {
        svert = s(i, j + 1) + s(i, j - 1);
    }
    return shorz + svert;
}

///-- Ising Hamiltonian without external magnetic field -> J*sum_<ij>{S_i * S_j}
int HI(const Mat<int> &s, const double J)
{
    const int n = s.n_rows;
    int sNN = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            sNN += s(i, j) * sumNN(s, i, j);
        }
    }
    return -J * sNN;
}

///-- Ising model function
void IsingSim(const double Bext, const char *fname, const int N_atoms)
{
    /// open file to save the results
    std::ofstream fres(fname);
    if (!fres.is_open())
    {
        std::cout << "error opening file!" << std::endl;
        return;
    }
    fres.setf(std::ios_base::scientific);

    ///-- initial conditions
    const int N = sqrt(N_atoms);
    const int M_MCMC = (int)1e06;
    const double dT{0.1}, J{1.0};
    arma::Mat<int> sOld = GenerateS(N);
    ///-- M_old = sum_i{S_i}
    double Mold = sum(sum(sOld));
    ///-- Hold = -J*sum_<ij>{S_i*S_j} - B_ext * sum_i{S_i}
    double Hold = HI(sOld, J) - Bext * Mold;
    ///-- 1 <= T <= 4
    vec T = regspace(1.0, dT, 4.0);

    /// initialize the random generator
    std::random_device rd;
    std::mt19937 g(rd());
    /// @param d -> random int between 0 and N-1(for indexes)
    std::uniform_int_distribution<int> d(0, N - 1);
    /// @param R -> random real number between 0 and 1(for metropolis algorithm)
    std::uniform_real_distribution<double> R(0.0, 1.0);

    /// loop over temperature
    // for(int t = 0; t < T.size(); t++)
    for (int t = T.size() - 1; t != 0; --t)
    {
        vec M = zeros<vec>(M_MCMC);
        vec H = zeros<vec>(M_MCMC);
        ///-- loop over MCMC
        for (int v = 1; v <= M_MCMC; ++v)
        {
            ///-- select a spin randomly
            int x = d(g);
            int y = d(g);
            ///-- save the current spin
            int soldOld = sOld(x, y);
            ///-- dE = 2J * sum_<ij> * S_i + 2 * B_ext * S_i
            double dE = 2.0 * sOld(x, y) * (J * sumNN(sOld, x, y) + Bext);
            ///-- accept the step with Metropolis algorithm
            if (exp(-dE / T[t]) >= R(g))
            {
                Hold += dE;
                sOld(x, y) = -sOld(x, y);
            }
            H[v] = Hold;
            M[v] = std::abs(Mold - sOld(x, y) + soldOld);
        }
        /// fres -> T  <H>  <M>  X  C_v
        fres << T[t] << "\t" << mean(H) << "\t" << mean(M) << "\t" << var(M) / T[t] << "\t" << var(H) / (T[t] * T[t]) << "\n";
    }
    fres.close();
}
