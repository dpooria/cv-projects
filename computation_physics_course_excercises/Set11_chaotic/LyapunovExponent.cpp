/* Exercise set 11 Q.5 */
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

//-- function to compute Lyapunov Exponent of the Logistic Map
const void lyapunovLogisticMap(const int n, const char *fname, const double mumin = 0.0)
{
    const int N = 1000;
    //-- initialize mu
    double mu = mumin;
    const double dmu = (4.0 - mumin) / (N - 1);
    //-- lyapunov exponent
    double l;
    srand((unsigned)time(NULL));
    FILE *f = fopen(fname, "w");
    for (int i = 0; i < N; ++i)
    {
        double x = (double)rand() / (RAND_MAX + 0.0);
        l = 0.0;
        for (int j = 0; j < n; ++j)
        {
            x = x * mu * (1.0 - x);
            l += log(abs(mu - 2.0 * mu * x));
        }
        l /= n;
        fprintf(f, "%5.3f\t%5.3f\n", mu, l);
        mu += dmu;
    }
    fclose(f);
}

int main()
{
    lyapunovLogisticMap(100, "LyapunovExponent-logisticmap-n100.txt");
    lyapunovLogisticMap(1000, "LyapunovExponent-logisticmap1-n1000.txt");
    lyapunovLogisticMap(1000, "LyapunovExponent-logisticmap2-n1000.txt", 2.5);

    //-- calculate Lyapunov exponent of Chaotic Oscillator(Q2)
    FILE *fCO = fopen("chaoticOscillation-phaseDiagram1.txt", "r");
    FILE *fCO2 = fopen("chaoticOscillation-phaseDiagram2.txt", "r");
    if (fCO == NULL || fCO2 == NULL)
    {
        perror("error!: couldn't open chaotic oscillator file\n");
        return -1;
    }
    float *tmp = new float;
    double *dxCo = new double;
    double s = 0.0;
    int nco = 0;
    while (fscanf(fCO, "%f\t%f\t%lf\n", tmp, tmp, dxCo) != EOF)
    {
        nco++;
        s += *dxCo;
    }
    fclose(fCO);
    puts("Lyapunov Exponent of Chaotic Oscillator -> Q2 (theta0 = 5 & w0 = 0.5)");
    printf("\t\tLE = %lf\n", s / nco);
    s = 0.0;
    nco = 0;
    while (fscanf(fCO2, "%f\t%f\t%lf\n", tmp, tmp, dxCo) != EOF)
    {
        nco++;
        s += *dxCo;
    }
    delete dxCo;
    fclose(fCO2);
    puts("Lyapunov Exponent of Chaotic Oscillator -> Q2 (theta0 = 0 & w0 = 0.01) ");
    printf("\t\tLE = %lf\n", s / nco);
    //-- calculate Lyapunov exponent of Lorenz attractor
    FILE *fLa = fopen("LorenzAttractor-res.txt", "r");
    if (fLa == NULL)
    {
        perror("error!:\tcouldn't open Lorenz Attractor file.\n");
        return -1;
    }
    int nLa = 0;
    double *dx = new double;
    double *dy = new double;
    double *dz = new double;
    double sx{0.0}, sy{0.0}, sz{0.0};
    //-- format of file "LorenzAttractor-res.txt" is like x dx y dy z dz
    while (fscanf(fLa, "%f\t%lf\t%f\t%lf\t%f\t%lf\n", tmp, dx, tmp, dy, tmp, dz) != EOF)
    {
        nLa++;
        sx += *dx;
        sy += *dy;
        sz += *dz;
    }
    fclose(fLa);
    puts("Lyapunov Exponent of Lorenz Attractor(Q3)");
    printf("\t\tLE_x = %f\n\t\tLE_y = %f\n\t\tLE_z = %f\n", sx / nLa, sy / nLa, sz / nLa);
    delete tmp;
    delete dx;
    delete dy;
    delete dz;
    //-- plot the results
    const char gnuplot[] = R"END(gnuplot -e "
    set terminal png enhanced; unset key;set grid;
    set xlabel '{/Symbol m}'; set ylabel '{/Symbol l}';
    set style line 1 lc rgb 'blue' pt 6 ps 0.3;
    set style line 2 lc rgb 'red' pt 6 ps 0.3;
    set output 'LyapunovExponent-logisticmap-n100.png';
    set title 'Lyapunov Exponent of Logistic map, n = 100';
    p 'LyapunovExponent-logisticmap-n100.txt' w l ls 1;
    set output 'LyapunovExponent-logisticmap1-n1000.png';
    set title 'Lyapunov Exponent of Logistic map, n = 1000';
    p 'LyapunovExponent-logisticmap1-n1000.txt' w l ls 1;
    set output 'LyapunovExponent-logisticmap2-n1000.png';
    p 'LyapunovExponent-logisticmap2-n1000.txt' w l ls 1;"
    )END";
    if (system(gnuplot))
    {
        perror("error!\n");
    }
    return 0;
}