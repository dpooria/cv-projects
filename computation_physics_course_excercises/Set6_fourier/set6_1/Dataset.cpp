#ifndef DATA_SETS_SIZE
#define DATA_SETS_SIZE 10485
#endif
#include "Dataset.h"

Dataset::Dataset(uint name)
{
    m_name = name;
    ifstream fdata;
    fdata.open("data/data" + to_string(m_name) + ".txt", ios::in);
    if (!fdata.is_open())
    {
        cout << "error making dataset object: data file not found!\n";
        *this = Dataset();
    }
    else
    {
        double xtmp;
        while (fdata >> xtmp)
        {
            m_data.push_back(xtmp);
        }
        fdata.close();
        m_N = m_data.size();
        m_pdf = pdfcal(m_data);
    }
}
uint Dataset::getname()
{
    return m_name;
}
vector<double> Dataset::getdata()
{
    return m_data;
}
vector<vector<double>> Dataset::getpdf()
{
    return m_pdf;
}
vector<vector<double>> Dataset::getjointpdf()
{
    return m_jointpdf;
}
double Dataset::geterrordata()
{
    return errorcal(m_data);
}
double Dataset::geterrorpdf()
{
    return errorcal(m_pdf[1]);
}
size_t Dataset::getsize()
{
    return m_N;
}
void Dataset::initjointPDF()
{
    m_jointpdf = jointpdfcal(m_data);
}
void Dataset::initDelta()
{
    m_delta.resize(2);
    size_t n = m_pdf[1].size();
    for (size_t t = 0; t < n; t++)
    {
        m_delta[0].push_back(t);
        m_delta[1].push_back(0.0);
        for (size_t i = 0; i < n - t; i++)
        {
            m_delta[1].at(t) += abs(m_jointpdf.at(i + t).at(i) -
                                    m_pdf[1].at(i) * m_pdf[1].at(i + t));
        }
        m_delta[1].at(t) /= n;
    }
}
double errorcal(const vector<double> &x)
{
    int n = x.size();
    double mean = accumulate(x.begin(), x.end(), 0.0) / double(n);
    double smd = 0;
    for (auto xit : x)
    {
        smd += powl(xit - mean, 2);
    }
    smd = sqrt(smd / (n * (n - 1)));
    return smd;
}

void Dataset::savepdf(string name)
{
    double perr = geterrorpdf();
    ofstream res;
    res.open(name, ios::out);
    for (size_t j = 0; j < m_pdf[1].size(); ++j)
    {
        res << m_pdf[0].at(j) << "\t" << m_pdf[1].at(j) << "\t" << perr << "\n";
    }
    res.close();
}

void Dataset::savejointpdf(string name)
{
    savethis(m_jointpdf, name);
}

void Dataset::savedelta(string name)
{
    savethis(m_delta, name);
}

vector<vector<double>> pdfcal(const vector<double> &x)
{
    double xmin = *min_element(x.begin(), x.end());
    double xmax = *max_element(x.begin(), x.end());
    uint N = x.size();
    double dx = (xmax - xmin) / N;
    int kmin = static_cast<int>(xmin / dx), kmax = static_cast<int>(xmax / dx);
    map<int, uint> n;
    for (auto xit : x)
    {
        int k = static_cast<int>(xit / dx);
        auto nit = n.find(k);
        if (nit == n.end())
        {
            n.insert({k, 1});
        }
        else
        {
            (*nit).second++;
        }
    }
    for (int k = kmin; k <= kmax; ++k)
    {
        auto nit = n.find(k);
        if (nit == n.end())
        {
            n.insert({k, 0});
        }
    }
    vector<vector<double>> p(2);
    for (auto &nit : n)
    {
        p[0].push_back(nit.first * dx);
        p[1].push_back(nit.second / (N * dx));
    }
    return p;
}
vector<vector<double>> jointpdfcal(const vector<double> &x)
{
    size_t N = x.size();
    double xmin = *min_element(x.begin(), x.end());
    double xmax = *max_element(x.begin(), x.end());
    double dx = (xmax - xmin) / N;
    int kmin = static_cast<int>(xmin / dx);
    int kmax = static_cast<int>(xmax / dx);
    vector<vector<double>> n(kmax - kmin + 1, vector<double>(kmax - kmin + 1));
    for (auto &rows : n)
    {
        for (auto &cols : rows)
        {
            cols = 0.0;
        }
    }
    double w = (1 / (N * dx));
    for (size_t t = 0; t < N; t++)
    {
        for (size_t i = 0; i < N - t; i++)
        {
            int k1 = static_cast<int>(x.at(i) / dx) - kmin;
            int k2 = static_cast<int>(x.at(i + t) / dx) - kmin;
            n.at(k1).at(k2) += w;
        }
    }
    return n;
}

void savethis(const vector<vector<double>> &dat, string fname)
{
    ofstream res;
    res.open(fname, ios::out);
    size_t rows = dat.size();
    size_t cols = dat.at(0).size();
    for (size_t i = 0; i < cols; i++)
    {
        for (size_t j = 0; j < rows; j++)
        {
            res << dat.at(j).at(i) << "  ";
        }
        res << "\n";
    }
    res.close();
}