#include "set4.h"
//--given p(x)
double w(double x)
{
    return sin(powl(x, 2.0) / 100.0) + 1.0 / cos(powl(x, 3.0) / 100.0) + powl(x, -3.0);
}
int main()
{
    srand(unsigned(time(0)));
    //-- M = max(p(x))
    double x, y, M = w(5.0);
    vector<double> xv;
    //-- Von-Neumann Method
    for (int i = 0; i < 1E04; ++i)
    {
        //-- 1 < x < 5
        x = (5.0 - 1.0) * rand() / (RAND_MAX + 0.0) + 1.0;
        //-- 0 < y < w(5)
        y = rand() / (RAND_MAX + 0.0) * M;
        if (y <= w(x))
        {
            xv.push_back(x);
        }
    }
    string name("x"), dirname("set4_2-RES");
    if (system(string("mkdir -p " + dirname).c_str()))
    {
        cout << "error making directory!" << endl;
    }
    ofstream result;
    result.open(dirname + "/" + "GeneratedRandomData.txt");
    for (auto xvit : xv)
    {
        result << xvit << endl;
    }
    pdfcal(xv, w, dirname, name);
    return 0;
}
