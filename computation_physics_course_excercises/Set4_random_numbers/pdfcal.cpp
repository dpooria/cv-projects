/*---Calculates the probability density of given data(const vector<double>& x) 
	and plots and saves it in ${dirname}/${name}pdf.png 
	and ${dirname}/${name}pdf.txt it also plots 
	the passed function(double (*w)(double))
	to compare it to the obtained PDF -------*/
#include "set4.h"
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

void pdfcal(const vector<double> &x, double (*w)(double),
            const string &dirname, const string &name)
{
    double xmin = *min_element(x.begin(), x.end());
    double xmax = *max_element(x.begin(), x.end());
    uint N = x.size();
    double dx = (xmax - xmin) / N;
    int kmin = static_cast<int>(xmin / dx), kmax = static_cast<int>(xmax / dx);

    map<int, uint> n;
    for (int k = kmin; k <= kmax; ++k)
    {
        n.insert({k, 0});
    }
    for (uint i = 0; i < N; ++i)
    {
        n.at(static_cast<int>(x.at(i) / dx))++;
    }
    vector<double> p(N), pt(N), xc(N);
    for (auto nit : n)
    {
        xc.push_back(nit.first * dx);
        p.push_back(nit.second / (N * dx));
        //--expected PDF (Gaussian PDF)
        pt.push_back(w(nit.first * dx));
    }
    plt::named_plot("PDF of " + name, xc, p, "b.");
    plt::named_plot("expected pdf function", xc, pt, "r:");
    plt::legend();
    plt::xlabel(name);
    plt::ylabel("p(" + name + ")");
    plt::grid(true);
    plt::save(dirname + "/" + name + "pdf.png");
    plt::clf();
    ofstream res;
    res.open(dirname + "/" + name + "pdf.txt", ios::out);
    for (uint i = 0; i < N; ++i)
    {
        res << x.at(i) << "\t" << n.at(static_cast<int>(x.at(i) / dx)) / (N * dx) << endl;
    }
    res.close();
}
