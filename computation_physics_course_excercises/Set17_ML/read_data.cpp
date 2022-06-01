#include <armadillo>
#include <iostream>
#include <fstream>

extern "C" void read_hoda(void *xptr, void *yptr, int n, int m, int c)
{
    arma::Mat<int> x((int *)xptr, m, n);
    arma::inplace_trans(x);
    arma::Mat<int> y((int *)yptr, c, n);
    arma::inplace_trans(y);
    const char *fnameX = "img.bin";
    const char *fnameY = "target.bin";
    x.save(fnameX);
    y.save(fnameY);
    std::ofstream header("data_info.h");
    if (!header.is_open())
    {
        std::cout << "error opening file!" << std::endl;
        return;
    }
    header << "///generated with read_data.cpp" << std::endl;
    header << "const char* fnameX = \"" << fnameX << "\";" << std::endl;
    header << "const char* fnameY = \"" << fnameY << "\";" << std::endl;
    header << "constexpr int n = " << n << ";" << std::endl;
    header << "constexpr int m = " << m << ";" << std::endl;
    header << "constexpr int c = " << c << ";" << std::endl;
    header.close();
}

extern "C" void read_boston(const void *Xptr, const void *yptr, int n, int m)
{
    arma::mat X((double *)Xptr, m, n);
    arma::vec y((double *)yptr, n);
    arma::inplace_trans(X);
    const char *fnameBostonFeat = "X_boston.bin";
    const char *fnameBostonPrice = "y_boston.bin";
    X.save(fnameBostonFeat);
    y.save(fnameBostonPrice);
    std::ofstream header;
    header.open("data_info.h", std::ios_base::app);
    header << "//-- Boston" << std::endl;
    header << "const char* fnameBostonFeat = \"" << fnameBostonFeat << "\";" << std::endl;
    header << "const char* fnameBostonPrice = \"" << fnameBostonPrice << "\";" << std::endl;
    header.close();
}
