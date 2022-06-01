
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

using namespace std;

int main(int argc, char** argv){
    
    if(argv[1]==nullptr){
	cout << "Please enter the folder path" << endl;
	return 0;
}

    ifstream data;
    // the directory name will be passed as an argument to this program
    string dirname(argv[1]);
    //opening the data file in that directory
    data.open(dirname+"/data");
    //defining a vector to put data samples in it
    vector<double> dat_v;
    double temp;
    //read from data and put it in the vector datv
    while(data >> temp){
        dat_v.push_back(temp);
    }
    data.close();
    // Calculate and Write the Mean of data
    double dat_sum = 0; //sum of data elements
    for(auto& v : dat_v){
        dat_sum += v;
        }
    int n = dat_v.size(); //# of elements in data file
    double dat_mean = dat_sum/n;
    
    //writing mean to data.${j}/mean file
    ofstream fmean; 
    fmean.open(dirname+"/mean");
    fmean << dat_mean;
    fmean.close();
    
    double dat_squareSum = 0;//sum((x-mean)^2)
    for(auto& v: dat_v){
        dat_squareSum += pow(v-dat_mean,2);
    }
    double dat_variance = dat_squareSum/(n-1);
    double dat_sd = pow(dat_variance,2); //standard deviation
    double dat_msd = dat_sd/sqrt(n);    // mean standard deviation
    ofstream fsd;
    fsd.open(dirname+"/standard deviation");
    fsd << dat_sd;
    fsd.close();
    ofstream fmsd;
    fmsd.open(dirname+"/mean standard deviation");
    fmsd << dat_msd;
    fmsd.close();
    
    return 0;
}
