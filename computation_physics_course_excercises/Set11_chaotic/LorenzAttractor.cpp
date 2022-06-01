/* Exercise set 11 Q.3 */

#include <cstdio>
#include <cmath>

using namespace std;

int main()
{

    auto fx = [](double xx, double yy, double zz) -> double {
        return 10.0 * (yy - xx);
    };
    auto fy = [](double xx, double yy, double zz) -> double {
        return -xx * zz + 28.0 * xx - yy;
    };
    auto fz = [](double xx, double yy, double zz) -> double {
        return xx * yy - (8.0 / 3.0) * zz;
    };

    const int N{1000};
    double *x = new double[N]{0.0};
    double *y = new double[N]{0.0};
    double *z = new double[N]{0.0};
    //-- suppose 0 <= t <= 10
    const double dt = 10.0 / (N - 1.0);
    //-- set initial conditions (arbitrary)
    x[0] = -1.0;
    y[0] = 1.0;
    z[0] = 0.0;
    //-- solve the coupled equations with Euler method
    for (int i = 0; i < N - 1; ++i)
    {
        x[i + 1] = x[i] + dt * fx(x[i], y[i], z[i]);
        y[i + 1] = y[i] + dt * fy(x[i], y[i], z[i]);
        z[i + 1] = z[i] + dt * fz(x[i], y[i], z[i]);
    }
    //-- save results
    double px, py, pz;
    FILE *fres = fopen("LorenzAttractor-res.txt", "w");
    for (int i = 0; i < N; ++i)
    {
        px = fx(x[i], y[i], z[i]);
        py = fy(x[i], y[i], z[i]);
        pz = fz(x[i], y[i], z[i]);
        //-- write x dx/dt , y dy/dt and z dz/dt
        fprintf(fres, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",
                x[i], px, y[i], py, z[i], pz);
    }
    fclose(fres);

    delete[] x;
    delete[] y;
    delete[] z;
    //-- plot phase diagram of x, y and z
    const char gnuplot[] = R"end(gnuplot -e "
    set terminal png enhanced; set grid; unset key;
    set style line 1 lc rgb 'blue';
    set style line 2 lc rgb 'red' pt 6 ps 0.5;
    set output 'LorenzAttractor-phaseDiagramX.png';
    set xlabel 'x'; set ylabel 'dx/dt';
    plot 'LorenzAttractor-res.txt' u 1:2 w l ls 1, '' u 1:2 w p ls 2;
    set output 'LorenzAttractor-phaseDiagramY.png';
    set xlabel 'y'; set ylabel 'dy/dt';
    plot 'LorenzAttractor-res.txt' u 3:4 w l ls 1, '' u 3:4 w p ls 2;
    set output 'LorenzAttractor-phaseDiagramZ.png';
    set xlabel 'z'; set ylabel 'dz/dt';
    plot 'LorenzAttractor-res.txt' u 5:6 w l ls 1, '' u 5:6 w p ls 2;"
    )end";
    if (system(gnuplot))
    {
        perror("error plot!\n");
        puts(gnuplot);
    }
    return 0;
}