/* Exercise set 10 problem 1
*/
#include <armadillo>
#include <cmath>
#include <ctime>
#include <fstream>

using namespace arma;

int main()
{
    const int m{300}, n{300};
    mat phi = arma::zeros<mat>(n, m);
    //-- imply boundary conditions
    
    for (int j = 0; j < m; ++j)
    {
        phi(0, j) = static_cast<double>(j * j);
        phi(n - 1, j) = 0.0;
    }

    for (int i = 0; i < n; ++i)
    {
        phi(i, 0) = static_cast<double>(i);
        phi(i, m - 1) = 1.0;
    }
    phi(n-1,0) = 0.0;

    //-- initialize other elements by random numbers 
    srand((unsigned)time(0));
    for (int i = 1; i < n - 1; ++i)
    {
        for (int j = 1; j < m - 1; ++j)
        {
            phi(i, j) = (double)rand() / RAND_MAX;
        }
    }

    //-- choose elements randomly and replace them by average of their neighbors
    srand((unsigned)time(0));
    for (int k = 0; k < (int)1e06; ++k)
    {
        int i = floor((n - 2 - 1) * (double)rand() / RAND_MAX + 1);
        int j = floor((m - 2 - 1) * (double)rand() / RAND_MAX + 1);
        if(i == 0 || j == 0 || i == n-1 || j == m - 1 )
        {
            std::cout << "somthing's wrong\n";
        }
        phi(i, j) = (phi(i + 1, j) + phi(i - 1, j) + phi(i, j + 1) + phi(i, j - 1)) / 4.0;
    }
    
    //-- save the results
    std::ofstream fres;
    fres.open("laplace-phi.txt");
    for(int i = 0; i < m; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            fres << i << "\t" << j << "\t" << phi(i,j) << "\n";
        }
    }
    fres.close();

    //-- plot'em
    char gnuplot[] = R"END(
    gnuplot -e "set terminal png enhanced;
    set grid; set xlabel 'x'; set ylabel 'y'; set zlabel '{/Symbol F}(x,y)';
    set dgrid3d 30,30;set hidden3d;
    set output 'laplace-phi.png';
    splot 'laplace-phi.txt' u 1:2:3 w l ls 1 t '{/Symbol F}(x,y)';
    set output 'laplace-Logphi.png';
    splot 'laplace-phi.txt' u 1:2:(\$3==0 ? 0:(log(\$3))) w l t 'log({/Symbol F}(x,y))';"
    )END";
    if (system(gnuplot))
    {
        std::cout << "plot error!" << std::endl;
        std::cout << gnuplot << std::endl;
    }
    return 0;
}