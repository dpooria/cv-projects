#! /bin/bash
#-- find the path of armadillo and openblas libraries
armadillopath=$(whereis armadillo | awk '{print $2}')
openblaspath=$(whereis openblas | awk '{print $2}')
#-- checks if the required libraries installed
if [ ! -z "$armadillopath" ] && [ ! -z "$openblaspath" ]; then
echo 'main.sh: Compliling set5_1.cpp...'
g++ set5_1.cpp -o set5_1.o -std=c++11 -O2 -L${armadillopath} -larmadillo -lopenblas
echo 'main.sh: Running set5_1.o...'
./set5_1.o
echo 'main.sh: Compliling set5_2.cpp...'
g++ set5_2.cpp -o set5_2.o -std=c++11 -O2
echo 'main.sh: Running set5_2.o...'
./set5_2.o
echo '--------------------Done!------------------';
else  
echo "main.sh: error!: armadillo or openblas library not found";
fi
