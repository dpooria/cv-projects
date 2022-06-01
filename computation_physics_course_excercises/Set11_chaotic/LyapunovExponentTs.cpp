/* Exercise set 11, Q.5 */
#include <cstdio>
#include <cmath>

using namespace std;

//-- Lyapunov exponent of a given time sery
const void LyapunovTs(const double *x, const int N, const char *fname)
{
    double l;
    static const int n = 100;
    FILE *f = fopen(fname, "w");
    for (int t = 1; t < 100; ++t)
    {
        l = 0.0;
        double d0, dn;
        for (int i = 0; i < N - n - t; ++i)
        {
            d0 = abs(x[i + t] - x[i]);
            dn = abs(x[i + n + t] - x[i + n]);
            if (d0 != 0.0 && dn != 0.0)
            {
                l += (1.0 / n) * (log(dn / d0));
            }
        }
        fprintf(f, "%d\t%.10lf\n", t, l / N);
    }
    fclose(f);
}

//-- loads a file and saves it in an array
double *loadF(const char *fname, int *N)
{
    FILE *f = fopen(fname, "r");
    if (f == NULL)
    {
        printf("couldn't open file %s\n", fname);
        return nullptr;
    }
    int n = 0;
    float *tmp = new float;
    while (fscanf(f, "%f\t%f\n", tmp, tmp) != EOF)
    {
        n++;
    }
    *N = n;
    double *x = new double[n];
    rewind(f);
    for (int i = 0; i < n; ++i)
    {
        fscanf(f, "%f\t%lf\n", tmp, &x[i]);
    }
    delete tmp;
    return x;
}

bool LyapunovCal(const char *fin, const char *fout)
{
    int N = 0;
    double *x = loadF(fin, &N);
    if (x == nullptr)
    {
        return false;
    }
    LyapunovTs(x, N, fout);
    delete[] x;
    return true;
}

int main()
{
    LyapunovCal("chaotic_data/0.200.txt", "LyapunovEponentTs-0.200.txt");
    LyapunovCal("chaotic_data/0.900.txt", "LyapunovEponentTs-0.900.txt");
    LyapunovCal("chaotic_data/logistic_map4.txt", "LyapunovEponentTs-logistic_map4.txt");
    LyapunovCal("chaotic_data/position.txt", "LyapunovExponentTs-position.txt");
    LyapunovCal("chaotic_data/position2.txt", "LyapunovEponentTs-position2.txt");
    LyapunovCal("chaotic_data/velocity.txt", "LyapunovEponentTs-velocity.txt");
    //---------------------------------------------------------------------
    const char gnuplot[] = R"END(gnuplot -e "
    set terminal png enhanced; set grid;unset key;
    set xlabel '{/Symbol t}'; set ylabel '{/Symbol l}({/Symbol t})';
    set style line 1 lc rgb 'blue';
    set output 'LyapunovEponentTs-0.200.png';
    set title 'LE of "0.200.txt"';
    p 'LyapunovEponentTs-0.200.txt' w l ls 1;
    set output 'LyapunovEponentTs-0.900.png';
    set title 'LE of "0.900.txt"';
    p 'LyapunovEponentTs-0.900.txt' w l ls 1;
    set output 'LyapunovEponentTs-logistic_map4.png';
    set title 'LE of "logistic\\_map4.txt"';
    p 'LyapunovEponentTs-logistic_map4.txt' w l ls 1;
    set output 'LyapunovExponentTs-position.png';
    set title 'LE of "position.txt"';
    p 'LyapunovExponentTs-position.txt' w l ls 1;
    set output 'LyapunovExponentTs-position2.png';
    set title 'LE of "position2.txt"';
    p 'LyapunovEponentTs-position2.txt' w l ls 1;
    set output 'LyapunovExponentTs-velocity.png';
    set title 'LE of "velocity.txt"';
    p 'LyapunovEponentTs-velocity.txt' w l ls 1;"
    )END";
    if (system(gnuplot))
    {
        puts("error!: plot");
        puts(gnuplot);
        return -1;
    }
    return 0;
}