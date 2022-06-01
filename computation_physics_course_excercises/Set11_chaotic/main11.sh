#! /bin/bash

for l in "logisticmap" "chaoticOscillation" "LorenzAttractor" \
	"LyapunovExponent" "LyapunovExponentTs"; do
	if g++ ${l}.cpp -o ${l}.o -std=c++11 -Wall; then
		echo Running ${l}
		./${l}.o
	else
		echo error compiling ${l}
	fi
done
echo "_________________DONE!_____________________"
