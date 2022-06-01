#! /bin/bash
#path to the python2.7 library in my system is /usr/include/python2.7 it might be different in yours
g++ -c pdfcal.cpp -I/usr/include/python2.7 -lpython2.7
g++ -c set4_1.cpp
g++ -c set4_2.cpp
g++ pdfcal.o set4_1.o -I/usr/include/python2.7 -lpython2.7 -o set4_1
g++ pdfcal.o set4_2.o -I/usr/include/python2.7 -lpython2.7 -o set4_2
./set4_1
./set4_2
echo "---done!---" 
