#!/bin/bash

if [[ -d dat ]]; #checks if directory dat exists
then
rm -r dat #remove it
fi
mkdir dat
if [[ ! -f program ]]; #checks if binary file of the c++ program doesn't exist
then
g++ main.cpp -o program; #creates it
fi
n=$(cat data.txt | wc -l) #assigns number of lines in data.txt to the variable n
let i=0 #for controlling the number of lines that have been read
let j=1 #for changing the file name
filename='dat/data1'
while read -r l; do #reads data.txt line by line and puts it in variable l
let i++ #counting the lines
if ((i>((n/100)))); #after every n/100 lines that've been read we change the output file
then
mkdir dat/data.$j
mv dat/data$j dat/data.$j/data
./program dat/data.$j #execute program to calculate mean, sd, msd by passing the directory path
echo $j%
let i=0 #reset i to count again
let j++
filename=dat/data${j}
fi
echo $l >> $filename; #write the data into the file
done < data.txt #giving the data.txt as input
mkdir dat/data.100
mv dat/data100 dat/data.100/data
./program 'dat/data.100';
./plot.sh; #executing plot.sh bash file for plotting the data
echo '100% done!'

