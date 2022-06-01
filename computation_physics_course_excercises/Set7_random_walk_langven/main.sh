#! /bin/bash
armadillopath=$(whereis armadillo | awk '{print $2}')
openblaspath=$(whereis openblas | awk '{print $2}')
if [ ! -z "$armadillopath" ] && [ ! -z "$openblaspath" ]; then
echo "Compiling rw.cpp"
g++ rw.cpp -o rw.o -std=c++11 -O2 -L$armadillopath -larmadillo -lopenblas
echo "Running rw.o"
./rw.o
echo "Compiling langven.cpp"
g++ langven.cpp -o langven.o -std=c++11 -O2 -L$armadillopath -larmadillo -lopenblas
echo "Running langven.o"
./langven.o
echo "______DONE!_______"
else
echo "error! couldn't find armadillo or openblas library"
fi
