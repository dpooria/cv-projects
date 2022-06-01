#include <armadillo>
#include <string>
#include <fstream>
#include <iostream>
#include <time.h>

using namespace std;
using namespace arma;

//---Function for Exercise 5 - 1 part A
void set5_1_A(const mat &dat, uint dataset)
{

    //-- dataset => number of datasets
    const uint n = static_cast<uint>(dat.size() / dataset);
    //-- split data into ${dataset} parts by reshaping it to a matrix
    //-- with columns as each dataset and rows as values of data
    mat data(dat.memptr(), n, dataset);
    //-- initialize c to zeros
    mat c(n, n);
    c = zeros<mat>(n, n);
    //-- loop through each dataset
    for (int k = 0; k < dataset; ++k)
    {
        //-- copy value of each dataset to x
        vec x = data.col(k);
        //--loop through c columns
        for (uint i = 0; i < n; ++i)
        {
            //-- loop through c rows
            for (uint j = 0; j < n; ++j)
            {
                c(i, j) += x(i) * x(j);
            }
        }
    }
    //-- average c over number of dataset
    c /= dataset;
    //-- save and plot results
    string filename = "cA-" + to_string(dataset);
    ofstream f;
    f.open("res1/A/" + filename + ".txt", ios::out);
    f << c << endl;
    f.close();
    string gnuplot = "gnuplot -e 'set terminal jpeg;\
    set xrange [0:" + to_string(n) +
                     "];set yrange [0:" + to_string(n) + "];\
    set xlabel \"row\"; set ylabel \"col\";\
    set title \"C(i,j) with " +
                     to_string(dataset) + " dataset\";\
    p \"res1/A/" + filename +
                     ".txt\" matrix w image;' > 'res1/A/" + filename + ".jpg'";
    if (system(gnuplot.c_str()))
    {
        cout << dataset << ": gnuplot: error!" << endl;
    }
}
//-- Function for Exercise 5 - 1 part B
void set5_1_B(const mat &dat)
{
    //-- n => number of data in each dataset
    uint n = static_cast<uint>(dat.size() / 100);
    //-- split data into 100 parts by reshaping data into a matrix with
    //-- columns as each dataset
    mat data(dat.memptr(), n, 100);
    mat c(n, 100);
    c = zeros<mat>(n, 100);
    //-- loop through datasets
    ofstream fcorrelation, fdata;
    for (uint j = 0; j < 100; ++j)
    {
        //-- create a pointer to each column of data and c
        auto ptrc = c.begin_col(j);
        auto x0 = data.begin_col(j);
        double *x;
        //-- loop through tau
        for (uint t = 0; t < n; t++)
        {
            //-- loop through the dataset
            for (x = x0; x < x0 + n - t; x++)
            {
                *ptrc += (*x) * (*(x + t));
            }
            *ptrc /= n - t;
            ptrc++;
        }
        fcorrelation.open("res1/B/correlations/correlation" +
                              to_string(j + 1) + ".txt",
                          ios::out);
        fdata.open("data/dataset" + to_string(j + 1) + ".txt", ios::out);
        vec cv = c.col(j);
        vec datv = data.col(j);
        fcorrelation << cv << endl;
        fdata << datv << endl;
        fcorrelation.close();
        fdata.close();
    }

    vec tau(n);
    tau = linspace(0, n - 1, n);
    ivec dataset = {2, 25, 55, 75, 97};
    mat cdat(n, 2);
    for (auto it : dataset)
    {
        string fname = "Correlation for dataset " + to_string(it);
        fcorrelation.open("res1/B/" + fname + ".txt", ios::out);
        cdat = join_horiz(tau, c.col(it));
        fcorrelation << cdat << endl;
        fcorrelation.close();
        string gnuplot = "gnuplot -e 'set terminal jpeg;\
        set xlabel \"t\"; set ylabel \"C_i(t)\"; set title \"" +
                         fname + "\";\
        set grid; unset key; set style line 2 lc rgb \"black\" pt 7 ps 0.5; \
        p \"res1/B/" +
                         fname + ".txt\" w p ls 2;' > " +
                         "\"res1/B/" + fname + ".jpg\"";
        if (system(gnuplot.c_str()))
        {
            cout << fname << ": gnuplot: error!" << endl;
            cout << gnuplot << endl;
        }
    }
}

int main()
{
    vec data;
    if (!data.load("data.txt"))
    {
        cout << "error! : 'data.txt' file not found!\n";
        return -1;
    }
    if (system("mkdir -p res1 res1/A res1/B res1/B/correlations data"))
        cout << "error!\n";
    cout << "set5_1_A: running...\n";
    /*---The exercise wants us to split data into 100 datasets but
        for part A, C(i,j) would be too large(about 1GB file) and therefore
        computations for 100 datasets is too heavy since C size is N/100*N/100
        and N ~ 1E06 so for part A I decided to split data into 1000 and 10000
        parts so size of C would be about 1E06*1E06 and 1E04*1E04 */
    clock_t t = clock();
    set5_1_A(data, 1E03);
    set5_1_A(data, 1E04);
    cout << "set5_1_B: running...\n";
    set5_1_B(data);
    cout << "t(set5_1) = " << (double)(clock() - t) / CLOCKS_PER_SEC << " sec" << endl;
    return 0;
}
