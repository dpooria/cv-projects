/* Exercise set 11 Q. 1 */
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

int const bifurcation(const double muMin, const double muMax, const char *res)
{
    double x = 0.0;
    const double dmu = (muMax - muMin) / 100;
    srand((unsigned)time(NULL));
    char fn[100];
    snprintf(fn, 100, "%s.txt", res);
    FILE *fres = fopen(fn, "w");
    for (double mu = muMin; mu <= muMax; mu += dmu)
    {
        for (int j = 0; j < 100; ++j)
        {
            x = (double)rand() / (RAND_MAX + 0.0);
            for (int i = 1; i <= 1000; ++i)
            {
                x = x * mu * (1.0 - x);
            }
            fprintf(fres, "%5.3f\t%5.3f\n", mu, x);
        }
    }
    fclose(fres);
    static const char plotfrmt[] = "gnuplot -e \"\
    set terminal png enhanced; unset key;\
    set title 'bifurcation diagram';\
    set style line 1 lc rgb 'blue' pt 7 ps 0.5;\
    set style line 2 lc rgb 'red' pt 6 ps 0.5;\
    set xlabel '{/Symbol m}'; set ylabel 'x';\
    set output '%s.png'; set grid;\
    p '%s' w l ls 1, '' w p ls 2;\"";

    char gnuplot[500];
    snprintf(gnuplot, 500, plotfrmt, res, fn);
    return system(gnuplot);
}

int const cycles(const int nmax, const double mu, const char *res, const char *title, bool withLine = true)
{
    double x;
    srand((unsigned)time(NULL));
    char fn[100];
    snprintf(fn, 100, "%s.txt", res);
    FILE *fres = fopen(fn, "w");
    for (int n = 1; n <= nmax; ++n)
    {
        x = (double)rand() / (RAND_MAX + 0.0);
        for (int i = 1; i <= n; ++i)
        {
            x = x * mu * (1.0 - x);
        }
        fprintf(fres, "%d\t%5.3f\n", n, x);
    }

    fclose(fres);
    static const char plotfrmt[] = "gnuplot -e \"\
    set terminal png enhanced; unset key;\
    set title '%s';\
    set style line 1 lc rgb 'blue' pt 7 ps 0.5;\
    set style line 2 lc rgb 'red' pt 6 ps 0.5;\
    set xlabel 'n'; set ylabel 'x_n';\
    set output '%s.png'; set grid;\
    p '%s' w p ls 1 %s";
    char gnuplot[500];
    if (withLine)
    {
        snprintf(gnuplot, 500, plotfrmt, title, res, fn, ", '' w l ls 2;\"");
    }
    else
    {
        snprintf(gnuplot, 500, plotfrmt, title, res, fn, ";\"");
    }

    return system(gnuplot);
}

int main(int argc, char *argv[])
{
    bifurcation(0.0, 5.0, "logisticmap-bifurcationdiagram");
    //-- one_cycle 1.0 <= mu <= 3.0
    cycles(100, 2.4, "logisticmap-one_cycle", "one\\_cycle {/Symbol m}= 2.4");
    cycles(1000, 2.4, "logisticmap-one_cycle-n1000", "one\\_cycle {/Symbol m}= 2.4", false);
    //-- twoycle 3.0 <= mu <= 3.45
    cycles(100, 3.2, "logisticmap-two_cycle", "two\\_cycle {/Symbol m}= 3.2");
    cycles(1000, 3.2, "logisticmap-two_cycle-n1000", "two\\_cycle {/Symbol m}=3.2", false);
    //-- foucycle 3.45 < mu < 3.56
    cycles(100, 3.52, "logisticmap-four_cycle", "four\\_cycle {/Symbol m}= 3.52");
    cycles(1000, 3.52, "logisticmap-four_cycle-n1000", "four\\_cycle {/Symbol m}= 3.52", false);
    //-- chaotic regime: mu >= 3.56
    cycles(100, 3.64, "logisticmap-chaotic_regime", "chaotic\\_regime {/Symbol m}= 3.56");
    cycles(1000, 3.64, "logisticmap-chaotic_regime-n1000", "chaotic\\_regime {/Symbol m}= 3.56", false);

    return 0;
}