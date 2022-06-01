#! /bin/bash

armapath=$(whereis armadillo | awk '{print $2}')
blaspath=$(whereis libopenblas | awk '{print $2}')

if [ -z "$armapath" ] || [ -z "$blaspath" ]; then
	echo "Couldn't find armadillo or openblas libraries!"
	echo "Please insert the armadillo path:"
	read armapath
	echo "Please insert the openblas path:"
	read blaspath
fi
# a_min = 0.0, a_max = 1.0, da = 0.005, H_min = 0.0, H_max = 1.0, dH = 0.005
echo -e "0.0\n1.0\n0.005\n0.0\n1.0\n0.005\n" >'fitConditions.txt'
echo "Compiling 'setFitConditioins.cpp'"
if g++ setFitConditions.cpp -o setFitConditions -std=c++11 -Wall; then
	echo "Running 'setFitConditions'..."
	if ./setFitConditions 'fitConditions.txt' == 0; then
		echo "Compiling 'GA.cpp'"
		if g++ GA.cpp -o GA -L$armapath -larmadillo -L$blaspath -lopenblas -O2 -std=c++11 -Wall; then
			echo "Running 'GA'..."
			if ./GA == 0; then
				echo "Done!"
			fi
		fi
	fi

fi
