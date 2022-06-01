#! /bin/bash

for p in "DecayinSim" "Pi" "QvarOscillator"; do
	echo "Compiling '${p}.cpp'..."
	g++ ${p}.cpp -o ${p} -Wall -lm -std=c++11 -O2
done
# works if armadillo and openblas libraries are in gcc's default path
let respond=0
for p in "Set3Fit" "MCMCFit" "HMCFit"; do
	echo "Compiling '${p}.cpp'..."
	if ! g++ ${p}.cpp -o ${p} \
		-Wall -lm -larmadillo -lopenblas -std=c++11 -O2;
	then
		echo "Error Compiling '${p}.cpp"
		let respond=-1
	fi
done	
if [ $respond!=-1 ]; then
	for p in "DecayinSim" "Pi" "QvarOscillator" "Set3Fit" "MCMCFit" "HMCFit"; do
		echo "Running '$p'..."
		./$p
	done
	echo "-----------------DONE!----------------"
fi
