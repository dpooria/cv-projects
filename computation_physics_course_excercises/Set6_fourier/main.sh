#! /bin/bash
armadillopath=$(whereis armadillo | awk '{print $2}')
openblaspath=$(whereis openblas | awk '{print $2}')
if [ ! -z "$armadillopath" ] && [ ! -z "$openblaspath" ]; then
cd set6_1
echo "_Compiling set6_1.cpp ..."
g++  -c Dataset.cpp set6_1.cpp -std=c++11 -O2
echo "_Running set6_1.o ..." 
g++  Dataset.o set6_1.o -o set6_1 -std=c++11 -O2 
./set6_1
cd ..
cd set6_2
echo "_Compiling set6_2.cpp ..."
g++ set6_2.cpp -o set6_2.o -O2 -std=c++11 -L${armadillopath} -larmadillo -lopenblas
echo "_Running set6_2.o ..."
./set6_2.o
cd ..
cd set6_3
echo "_Compiling set6_3.cpp ..."
g++ set6_3.cpp -o set6_3.o -O2 -std=c++11 -L${armadillopath} -larmadillo -lopenblas
echo "_Running set6_3.o ..."
./set6_3.o
cd ..
cd set6_4
echo "_Compiling set6_4.cpp ..."
g++ set6_4.cpp -o set6_4.o -O2 -std=c++11 -L${armadillopath} -larmadillo -lopenblas
echo "_Running set6_4.o ..."
./set6_4.o
echo "---dOnE!---"
fi

