#include "Dataset.h"
#include "set6.h"

bool splitdata(const vector<double> &data)
{
    int N = data.size(), n = static_cast<int>(N / 100);
    ofstream fds;
    fds.open("data/data1.txt", ios::out);
    int i = 0, j = 1;
    for (auto dit : data)
    {
        i++;
        fds << dit << "\n";
        if (i >= n)
        {
            i = 0;
            j++;
            fds.close();
            fds.open("data/data" + to_string(j) + ".txt", ios::out);
            if (!fds.is_open())
                return false;
        }
    }
    return true;
}
void plotpdferrorbar(string fname)
{
    string gnuplot = "gnuplot -e \"set terminal jpeg enhanced;\
                        set style line 2 lc rgb 'blue' pt 3 ps 0.1;\
                        set style line 3 lc rgb 'red' pt 3 ps 0.3;\
                        set xlabel 'x'; set ylabel 'p(x)';\
                        p [*:*] [0:*] 'A/" +
                     fname + ".txt' u 1:2:(\\$2-\\$3) w yerr ls 2 t 'errorbar',\
                         '' u 1:2 w p ls 3 t 'PDF';\"> 'B/" +
                     fname +
                     ".jpg'";
    if (system(gnuplot.c_str()))
    {
        cout << "error gnuplot!\n";
        cout << gnuplot.c_str();
    }
}
int main()
{
    clock_t t = clock();
    if (system("mkdir -p data A B C"))
    {
        cout << "mkdir error!\n";
        return -1;
    }
    //-- import data
    ifstream fdata;
    fdata.open("data.txt", ios::in);
    if (!fdata.is_open())
    {
        cout << "error: Could not find 'data.txt'!\n";
        return -1;
    }
    vector<double> data;
    double dtmp;
    while (fdata >> dtmp)
    {
        data.push_back(dtmp);
    }
    fdata.close();
    //-- split data into 100 parts
    if (!splitdata(data))
    {
        cout << "error: splitdata!\n";
    }

    for (short i = 1; i <= 100; ++i)
    {
        Dataset dataset(i);
        dataset.savepdf("A/pdfdata" + to_string(i) + ".txt");
    }
    for (uint i : {2, 25, 55, 75, 97})
    {
        plotpdferrorbar("pdfdata" + to_string(i));
        Dataset dataset(i);
        dataset.initjointPDF();
        dataset.initDelta();
        // dataset.savejointpdf("C/jointpdfData" + to_string(i) + ".txt");
        string f = "C/DeltaData" + to_string(i);
        dataset.savedelta(f + ".txt");
        string gnuplot = "gnuplot -e \"set terminal jpeg enhanced; \
        set xlabel '{/Symbol t}'; set ylabel '{/Symbol D}({/Symbol t})'; p '" +
                         f + ".txt' w p pt 5 ps 0.5 t\
        '{/Symbol D}_{" + to_string(i) +
                         "}({/Symbol t}) = " +
                         "|P(x(t+{/Symbol t}), x(t)) − P(x(t+{/Symbol t}))P(x(t))|';\
                         \" > '" +
                         f + ".jpg'";
        if (system(gnuplot.c_str()))
        {
            cout << "error gnuplot!\n";
        }
    }
    cout << "program set6_1 finished in " << double(clock() - t) / CLOCKS_PER_SEC;
    cout << " sec" << endl;

    return 0;
}
