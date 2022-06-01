#! /bin/bash
armadillopath=$(whereis armadillo | awk '{print $2}')
openblaspath=$(whereis openblas | awk '{print $2}')
if [ ! -z "$armadillopath" ] && [ ! -z "$openblaspath" ]; then
echo "Compiling cdf.cpp"
g++ cdf.cpp -o cdf.o -std=c++11 -O2 -L$armadillopath -larmadillo -lopenblas
echo "Running cdf.o"
./cdf.o
echo "Compiling diffEquation.cpp"
g++ diffEquation.cpp -o diffEquation.o -std=c++11 -O2
echo "Running diffEquation.o"
./diffEquation.o
echo "______DONE!_______"
else
echo "error! couldn't find armadillo or openblas library"
fi
