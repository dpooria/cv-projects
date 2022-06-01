#include "iostream"
#include "string"

using namespace std;
int main(int argc, char** argv){
	cout << string("cpp::Hello from ")+argv[1]+string(" !") << endl;
	return 0;
}
