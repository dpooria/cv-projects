#!/bin/bash
if [[ -f mean ]]; then
rm mean;
fi
if [[ -f 'mean standard deviation' ]]; then
rm 'mean standard deviation';
fi
if [[ -f 'standard deviation' ]]; then
rm 'standard deviation';
fi

for ((i=1; i<=100;i++)); 
do
#get the varibles from each folder and put it in a single file
printf '\n'$i' ' | cat - dat/data.$i/mean >> mean
printf '\n'$i' ' | cat - dat/data.$i/'standard deviation' >> 'standard deviation'
printf '\n'$i' ' | cat - dat/data.$i/'mean standard deviation' >> 'mean standard deviation';
done
#plot the varibles using gnuplot package
gnuplot -e 'set terminal jpeg; set xlabel "data"; set ylabel "mean";set grid; p "mean" w l' > mean.jpg
gnuplot -e 'set terminal jpeg; set xlabel "data"; set ylabel "standard deviation";set grid; p "standard deviation" w l' > 'standard deviation.jpg'
gnuplot -e 'set terminal jpeg; set xlabel "data"; set ylabel "mean standard deviation";set grid; p "mean standard deviation" w l' > 'mean standard deviation.jpg'


