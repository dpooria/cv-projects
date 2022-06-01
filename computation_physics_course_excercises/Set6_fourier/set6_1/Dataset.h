#ifndef DATASET_H
#define DATASET_H
#include "set6.h"
double errorcal(const vector<double> &);
vector<vector<double>> pdfcal(const vector<double> &);
vector<vector<double>> jointpdfcal(const vector<double> &);
void savethis(const vector<vector<double>>&, string);
class Dataset
{
private:
    size_t m_N;
    uint m_name;
    vector<double> m_data;
    vector<vector<double>> m_pdf, m_jointpdf, m_delta;

public:
    Dataset() {}
    Dataset(uint);
    uint getname();
    vector<double> getdata();
    vector<vector<double>> getpdf();
    vector<vector<double>> getjointpdf();
    double geterrordata();
    double geterrorpdf();
    size_t getsize();
    void setname(uint);
    void initjointPDF();
    void initDelta();
    void savepdf(string);
    void savejointpdf(string);
    void savedelta(string);
};
#endif