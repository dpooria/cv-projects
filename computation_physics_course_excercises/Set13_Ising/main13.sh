#! /bin/bash

armapath=$(whereis armadillo | awk '{print $2}')
blaspath=$(whereis openblas | awk '{print $2}')

if [ -z "$armapath" ] || [ -z "$blaspath" ]; then
	echo "Couldn't find armadillo or openblas libraries!"
	echo "Please insert the armadillo path:"
	read armapath
	echo "Please insert the openblas path:"
	read blaspath
fi

echo "Compiling 'Ising.cpp'"
if g++ Ising.cpp -o Ising -std=c++11 -O2 -Wall \
	-L$armapath -larmadillo -L$balspath -lopenblas; then
	echo "Running 'Ising'..."
	./Ising
	echo "Done!"
fi
 
