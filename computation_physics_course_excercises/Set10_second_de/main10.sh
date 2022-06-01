#-- Exercise Set 10
#! /bin/bash
parma=$(whereis armadillo | awk '{print $2}')
poblas=$(whereis openblas | awk '{print $2}')
if [ ! -z "$parma" ] && [ ! -z "$poblas" ]; then
    for p in "laplace" "linearBvp" "nonlinearBvpA" "nonlinearBvpB" "nonlinearBvpC"; do
        echo "-Compiling '${p}.cpp'"
        if g++ ${p}.cpp -o ${p} -std=c++11 -O2 -Wall \
            -L${parma} -larmadillo -L${poblas} -lopenblas; then
            echo "-Running '${p}'"
            ./${p}
        else
            echo "error compiling ${p}"
        fi
    done
    echo "______DONE!_______"
else
    echo "-error!: couldn't find armadillo or openblas library"
fi
