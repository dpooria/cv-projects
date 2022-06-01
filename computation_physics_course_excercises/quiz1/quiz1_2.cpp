#include <armadillo>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace arma;

int main(){
    arma::vec x;
    x.load("IRANCovid19.txt");
    int N = x.size();
    double n = 0, l = 0, X;
    double ntmp = -10, ltmp = -10, Xtmp;
    vec t = linspace(0,N,N-1);
    t *= 0.1;
    X = arma::sum(arma::pow(x - ntmp * arma::exp(ltmp*t),2));
    for(int i = 0; i < 1000; i++){
        ntmp = -10;
        for(int j = 0; j<1000; j++){
            Xtmp = arma::sum(arma::pow(x - ntmp * arma::exp(ltmp*t),2));
            
 if (Xtmp < X)
            {
                X = Xtmp;
                n = ntmp;
                l = ltmp;
            }
            ntmp += 0.02;
        }
        ltmp += 0.02;
 }
 cout << n << "\t" << l;
    return 0;
}