#! /bin/bash

armapath=$(whereis armadillo | awk '{print $2}')
obpath=$(whereis openblas | awk '{print $2}')

if [ -z "$armapath" ] || [ -z "$obpath" ]; then
	echo "Couldn't find armadillo or openblas libraries path"
	echo "try insert the paths manually"
	echo "path to armadillo:"
	read armapath
	echo "path to openblas:"
	read obpath
fi

echo "$armapath $obpath"
g++ -c -Wall -fpic read_data.cpp
g++ -shared -o libread_data.so read_data.o -L$armapath -larmadillo -Wall

echo "Running 'read_data.py'..."
if python3 read_data.py; then
	echo "Compiling ml.cpp..."
	g++ -c ml.cpp svm.cpp -O2 -std=c++11
	g++ ml.o svm.o -o ml -L$armapath -larmadillo -L$obpath -lopenblas -O2 -std=c++11
	echo "Running ml"
	if ./ml == 0; then
		echo "Done!"
	fi
fi
