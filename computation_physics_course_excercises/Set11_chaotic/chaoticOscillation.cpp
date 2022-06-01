/* Exercise set 11 Q. 2 */
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

//-- solve 2nd order differential equation with Euler method
void solveEuler(double (*F)(double, double, double), const int N,
                const double *t, double *theta, double *pth, const double f)
{
    const double dt = t[1] - t[0];
    for (int i = 0; i < N - 1; ++i)
    {
        theta[i + 1] = theta[i] + dt * pth[i];
        pth[i + 1] = pth[i] + dt * F(theta[i], pth[i], f);
    }
}

double F(double th, double p, double f)
{
    static const double w0{1.0}, w{0.666}, al{0.2};
    return -pow(w0, 2) * sin(th) - al * p + f * cos(w);
}

int main()
{
    const int N{1000};
    const double dt = 100.0 / (N - 1.0);
    double *t = new double[N]{0.0};
    double *theta = new double[N]{0.0};
    double *pth = new double[N]{0.0};
    for (int i = 0; i < N - 1; ++i)
    {
        t[i + 1] = t[i] + dt;
    }
    // -- theta[0] = 5 degrees
    theta[0] = 5.0 * M_PI / 180.0;
    //-- d(theta)[0]/dt = 0.5 rps
    pth[0] = M_PI;
    solveEuler(F, N, t, theta, pth, 0.52);
    FILE *fres = fopen("chaoticOscillation-phaseDiagram1.txt", "w");
    for (int i = 0; i < N; ++i)
    {
        fprintf(fres, "%5.3f\t%5.3f\t%5.3f\n", t[i], theta[i], pth[i]);
    }
    fclose(fres);
    // -- theta[0] = 0
    theta[0] = 0.0;
    //-- d(theta)[0]/dt = 0.01 rps
    pth[0] = 0.02 * M_PI;
    solveEuler(F, N, t, theta, pth, 0.52);
    fres = fopen("chaoticOscillation-phaseDiagram2.txt", "w");
    for (int i = 0; i < N; ++i)
    {
        fprintf(fres, "%5.3f\t%5.3f\t%5.3f\n", t[i], theta[i], pth[i]);
    }
    fclose(fres);
    //---------------------------------------------------------
    fres = fopen("chaoticOscillation-dthetadtf.txt", "w");
    //-- 0 <= f <= 1
    double f = 0.0;
    for (int j = 0; j <= 1000; ++j)
    {
        //-- set initial conditions
        //-- theta[0] = 0 degrees
        theta[0] = 0.0;
        //-- d(theta)[0]/dt = 0.01 rps
        pth[0] = 0.02 * M_PI;
        //-- solve the equation for each f
        solveEuler(F, N, t, theta, pth, f);
        //-- average over d(theta)/dt
        double s = (accumulate(&pth[0], &pth[N], 0.0)) / N;
        //-- save the results
        fprintf(fres, "%5.3f\t%5.3f\n", f, s);
        f += 0.002;
    }
    fclose(fres);

    delete[] t;
    delete[] theta;
    delete[] pth;

    //-- plot the results
    char gnuplot[] = R"END(gnuplot -e "
    set terminal png enhanced; set grid; unset key;
    set style line 1 lc rgb 'red' ps 0.4 pt 6;
    set style line 2 lc rgb 'blue';
    set xlabel '{/Symbol q}'; set ylabel 'd{/Symbol q}/dt';
    set output 'chaoticOscillation-phaseDiagram1.png';
    set title '{/Symbol q}_0 = 5^o & d{/Symbol q}(0)/dt = 0.5 rps';
    p 'chaoticOscillation-phaseDiagram1.txt' u 2:3 w p ls 1, '' u 2:3 w l ls 2;;
    set output 'chaoticOscillation-phaseDiagram2.png';
    set title '{/Symbol q}_0 = 0^o & d{/Symbol q}(0)/dt = 0.01 rps';
    p 'chaoticOscillation-phaseDiagram2.txt' u 2:3 w p ls 1, '' u 2:3 w l ls 2;
    set output 'chaoticOscillation-theta.png';
    set xlabel 't'; set ylabel '{/Symbol q}(t)';
    p 'chaoticOscillation-phaseDiagram2.txt' u 1:2 w p ls 1, '' u 1:2 w l ls 2;
    set output 'chaoticOscillation-dtheta.png';
    set xlabel 't'; set ylabel 'd{/Symbol q}/dt';
    p 'chaoticOscillation-phaseDiagram2.txt' u 1:3 w p ls 1, '' u 1:3 w l ls 2;
    set output 'chaoticOscillation-dthetadtf.png';
    set xlabel 'f';set ylabel '|d{/Symbol q}/dt|';
    p 'chaoticOscillation-dthetadtf.txt' w p ls 1, '' w l ls 2;;"
    )END";

    if (system(gnuplot))
    {
        perror("plot error!\n");
        puts(gnuplot);
    }
    return 0;
}