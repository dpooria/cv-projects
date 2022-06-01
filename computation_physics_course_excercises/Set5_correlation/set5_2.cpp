/* unweighted correlation function with different levels has been calculated for
    datasets that weighted correlation function of them was plotted in set5_1;
    I wasn't able to calculate 0 level and -2sd level for 'data.txt' since Npeak
    would be about 3E05 I decided to just calculte 2sd level for 'data.txt'.
*/
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <time.h>

using namespace std;

void calunweightedCorrelation(vector<double> &data, string fname, double s)
{
    uint N = data.size();
    vector<double> peak;
    vector<int> peakpos;
    //-- Just peaks that are higher than s has been picked
    if (data[0] - data[1] > 0 && data[0] > s)
    {
        peak.push_back(data[0]);
        peakpos.push_back(0);
    }

    int Npeak = 0;
    for (int i = 1; i < N - 1; i++)
    {
        if (data[i] > data[i + 1] && data[i] > data[i - 1])
        {
            if (data[i] > s)
            {
                peak.push_back(data[i]);
                peakpos.push_back(i);
                Npeak++;
            }
        }
    }
    //--save peaks and their position
    ofstream fpeak;
    fpeak.open("res2/peak/peaksInLevel" + fname + ".txt", ios::out);
    for (int i = 0; i < Npeak; i++)
    {
        fpeak << peakpos.at(i) << "\t" << peak.at(i) << endl;
    }
    fpeak.close();
    //-- find the maximum and minimum of peakpos therefore maximum of R
    int maxindex = *max_element(peakpos.begin(), peakpos.end());
    int minindex = *min_element(peakpos.begin(), peakpos.end());
    double dR = double(maxindex - minindex) / (Npeak + 0.0);
    int R;
    map<int, int> p;
    for (int i = 0; i < Npeak - 1; ++i)
    {
        for (int j = i + 1; j < Npeak; ++j)
        {
            R = static_cast<int>(round((peakpos[j] - peakpos[i]) / dR));
            //--find the iterator proportional to R
            auto it = p.find(R);
            //-- if it wasn't there insert it
            if (it == p.end())
            {
                p.insert({R, 1});
            }
            //-- if it was add to its value
            else
            {
                (*it).second++;
            }
        }
    }
    int Rmin = 0;
    int Rmax = static_cast<int>((maxindex - minindex) / dR);
    for (R = Rmin; R <= Rmax; R++)
    {
        auto it = p.find(R);
        //-- insert another possible values of R into p
        if (it == p.end())
        {
            p.insert({R, 0});
        }
    }
    //-- save results and plot'em
    ofstream fres;
    fname = "unweighted-Correlation-Level-" + fname;
    fres.open("res2/un-weightedCorrelation/" + fname + ".txt", ios::out);
    for (auto pit : p)
    {
        fres << pit.first * dR << "\t" << pit.second / (dR * Npeak * Npeak / (2.0 * N)) - 1.0 << "\n";
    }
    fres.close();
    string linestyle = "set style line 2 lc rgb \"blue\" lt 2 lw 0.3 pt 7 ps 0.5;";
    string lop;
    if (s > 0)
        lop = "lp";
    else
        lop = "p";   
    string gnuplot = "gnuplot -e 'set terminal jpeg enhanced; set xlabel \"R\"; \
    set ylabel \"{/Symbol Y}(R)\";set title \"" +
                     fname +
                     "\"; set grid;unset key;" + linestyle +
                     "p \"res2/un-weightedCorrelation/" +
                     fname + ".txt\" w "+lop+" ls 2;'>res2/un-weightedCorrelation/" + fname + ".jpg";
    if (system(gnuplot.c_str()))
    {
        cout << "gnuplot: error!\n";
        cout << gnuplot << endl;
    }
}

int main()
{
    if (system("mkdir -p res2 res2/peak res2/un-weightedCorrelation"))
    {
        cout << "mkdir: error!\n";
        return -1;
    }
    cout << "set5_2: running...\n";
    ifstream fdata;
    clock_t t = clock();
    vector<string> dataset = {"dataset2", "dataset25", "dataset55",
                              "dataset75", "dataset97", "data"};
    for (string f : dataset)
    {
        //-- import data
        fdata.open("data/" + f + ".txt", ios::in);
        if (f == "data")
            fdata.open("data.txt", ios::in);
        if (!fdata.is_open())
        {
            cout << "error!: you should first run set5_1.o!\n";
            return -1;
        }
        vector<double> data;
        double xtmp, meandat = 0;
        while (fdata >> xtmp)
        {
            data.push_back(xtmp);
            meandat += xtmp;
        }
        fdata.close();
        uint N = data.size();
        meandat /= N;
        double sd = 0;
        for (auto &d : data)
        {
            d -= meandat;
            sd += powl(d, 2.0);
        }
        sd /= N;

        calunweightedCorrelation(data, "2SD-for-" + f, 2 * sd);
        if (f != "data")
        {
            calunweightedCorrelation(data, "minus2SD-for-" + f, -2 * sd);
            calunweightedCorrelation(data, "zero-for-" + f, 0.0);
        }
    }
    cout << "t(set5_2) = " << (double)(clock() - t) / CLOCKS_PER_SEC << " sec" << endl;

    return 0;
}
