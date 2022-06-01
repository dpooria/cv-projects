
#include <armadillo>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace arma;
using namespace std;

#define loop(i, b) for (int i = 0; i < b; ++i)

void msolve(double b, double k, ivec &I, ivec &S, ivec &R)
{
    const int N = 101;
    const int Nt = 25;
    double dt = 1.0;
    I[0] = 1.0;
    S[0] = 100.0;
    R[0] = 0.0;
    loop(i, Nt - 1)
    {
        double Ip = (b * dt / N) * S[i] * I[i] - k * dt * I[i];
        double Sp = (-b * dt / N) * S[i] * I[i];
        double Rp = dt * k * I[i];
        I[i + 1] = floor(Ip) + I[i];
        S[i + 1] = floor(Sp) + S[i];
        R[i + 1] = floor(Rp) + R[i];
    }
}

int main()
{
    const int N = 101;
    const int Nt = 25;
    vec t = arma::linspace(0.0, 500.0, Nt);
    double dt = 1.0;
    double b = 2.0, k = 0.05;
    ivec S = zeros<ivec>(Nt);
    ivec I = zeros<ivec>(Nt);
    ivec R = zeros<ivec>(Nt);
    I[0] = 1.0;
    S[0] = 100.0;
    R[0] = 0.0;
    loop(i, Nt - 1)
    {
        double Ip = (b * dt / N) * S[i] * I[i] - k * dt * I[i];
        double Sp = (-b * dt / N) * S[i] * I[i];
        double Rp = dt * k * I[i];
        I[i + 1] = floor(Ip) + I[i];
        S[i + 1] = floor(Sp) + S[i];
        R[i + 1] = floor(Rp) + R[i];
    }
    ofstream fI, fR, fS;
    fI.open("I.txt");
    fR.open("R.txt");
    fS.open("S.txt");
    loop(i, Nt)
    {
        fI << t[i] << "\t" << I[i] << "\n";
        fR << t[i] << "\t" << R[i] << "\n";
        fS << t[i] << "\t" << S[i] << "\n";
    }
    fI.close();
    fR.close();
    fS.close();
    char gnuplot[] = R"end(
    gnuplot -e "
    set terminal png enhanced;
    set grid;
    set style line 1 lc rgb 'blue' pt 6 ps 0.5;
    set style line 2 lc rgb 'red' pt 7 ps 0.5;
    set output 'I.png';
    p 'I.txt' w p ls 1;
    set output 'R.png';
    p 'R.txt' w p ls 1;
    set output 'S.png';
    p 'S.txt' w p ls 1;"
    
    )end";
    if (system(gnuplot))
    {
        cout << "error plot!" << endl;
    }
    mat dataI, dataS, dataR;
    dataI.load("Infected.txt");
    dataS.load("Susceptible.txt");
    dataR.load("Recovered.txt");
    t = dataI.col(0);
    double k1 = 0.0, b1 = 0.0;
    loop(i, 10)
    {
        b1 += 0.1;
        k1 = 0.0;
        loop(j, 10)
        {
            msolve(b1, k1, I, S, R);
            k1 += 0.01;
            string gnuplot1 = "gnuplot -e \"\
                set terminal png enhanced;\
                set grid;\
                set style line 1 lc rgb 'blue' pt 6 ps 0.5;\
            set style line 2 lc rgb 'red' pt 7 ps 0.5;\
            set output 'I" + to_string(i) +
                              ".png';\
            p 'I" + to_string(i) +
                              ".txt' w p ls 1, 'Infected.txt' w p ls 2;\
            set output 'R" + to_string(i) +
                              ".png';\
            p 'R" + to_string(i) +
                              ".txt' w p ls 1, 'Recovered.txt' w p ls 2;\
            set output 'S" + to_string(i) +
                              ".png';\
            p 'S" + to_string(i) +
                              ".txt' w p ls 1, 'Susceptible.txt' w p ls 2;\"";
            // if (system(gnuplot1.c_str()))
            // {
            //     cout << "error plot!" << endl;
            // }
        }
    }
    

    return 0;
}