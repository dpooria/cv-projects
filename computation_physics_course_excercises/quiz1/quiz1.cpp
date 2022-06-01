
#include <vector>
#include <numeric>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
    vector<double> x;
    double atmp;
    while(cin >> atmp){
        x.push_back(atmp);
    }
    size_t N = x.size();
    double s = accumulate(x.begin(), x.end(), 0.0);
    double mean = s / N;
    double va = 0.0;
    ofstream resc;
    resc.open("output10.txt", ios::out);
    for (auto xit : x)
    {
        va += powl(xit - mean, 2);
        resc << xit*10 << "\n";
    }
    resc.close();
    va /= N - 1.0;
    double sm = sqrt(va / N);
    ofstream resB;
    resB.open("output.txt", ios::out);
    resB << mean << "\t" << sm << "\n";
    resB.close();
    return 0;
}