#include "set4.h"
//-- Gaussian Function
double w(double x)
{
    return exp(-powl(x, 2.0) / 2.0);
}
int main()
{
    //--set the seed with time
    srand((unsigned)time(0));
    double r1, r2;
    vector<double> x(1000), y(1000);
    for (uint i = 0; i < 1000; ++i)
    {
        r1 = rand() / (RAND_MAX + 0.0); //==> 0<r1<1
        r2 = rand() / (RAND_MAX + 0.0); //==> 0<r2<1
        //--convert x and y using Box_Muller 
        x.at(i) = sqrt(-2 * log(1 - r1)) * cos(2 * M_PI * r2);
        y.at(i) = sqrt(-2 * log(1 - r1)) * sin(2 * M_PI * r2);
    }
    string namex("x"), namey("y"), dirname("set4_1-RES");
    if(system(string("mkdir -p " + dirname).c_str())){
        cout << "error making directory!" << endl;
    }
    //--save x and y 
    ofstream result;
    result.open(dirname + "/GeneratedRandomData.txt");
    result << "x:\n";
    for (auto xit : x)
    {
        result << xit << endl;
    }
    result << "y:\n";
    for (auto yit : y)
    {
        result << yit << endl;
    }
    result.close();
    //--calculate pdf of x and y
    pdfcal(x, w, dirname, namex);
    pdfcal(y, w, dirname, namey);
    return 0;
}
