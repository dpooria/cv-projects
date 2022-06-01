#include <map>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <time.h>
using namespace std;
void pdfcal(const vector<double> &x, double (*w)(double x),
            const string &dirname, const string &name);
