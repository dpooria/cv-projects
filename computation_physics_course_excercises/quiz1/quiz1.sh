#! /bin/bash
g++ quiz1.cpp -o quiz1.o -O2 -std=c++11
./quiz1.o < data.txt

let n=48
let i=0;
while read -r l; do 
let i++
if ((i>n)); 
then
break;
fi
mkdir -p $l
cp output.txt $l
cp output10.txt $l
done < list.txt

