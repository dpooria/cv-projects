#-- Exercise Set 9
#! /bin/bash
armadillopath=$(whereis armadillo | awk '{print $2}')
openblaspath=$(whereis openblas | awk '{print $2}')
if [ ! -z "$armadillopath" ] && [ ! -z "$openblaspath" ]; then
for p in "TP" "Integral" "DE"
do
echo "-Compiling '${p}.cpp'"
g++ ${p}.cpp -o ${p} -std=c++11 -O2 \
-L${armadillopath} -larmadillo -L${openblaspath} -lopenblas
echo "-Running '${p}'"
./${p}
done
echo "______DONE!_______"
else
echo "-error!: couldn't find armadillo or openblas library"
fi
