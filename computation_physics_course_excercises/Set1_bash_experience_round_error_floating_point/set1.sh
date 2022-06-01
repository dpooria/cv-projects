#! /bin/bash
clear
if [[ -d result ]]; then
	rm -r result
fi
mkdir result
cp list_arrange result
cp sunspot.txt result
cp prog.cpp result
cd result
g++ prog.cpp -o prog
let cn=48 #$(cat -s list_arrange| wc -l)-1
let spn=$(cat -s sunspot.txt| wc -l)
let i=0

while read -r c; do
let i++
if ((i>cn)); 
then
break;
fi
mkdir "d$i"
echo '' > "d$i/$c"
./prog "$c"
done <list_arrange

let k=1
let m=0
fname="d1/ARGENTINA"
while read -r d; do
let m++
echo $d >> $fname
if ((m>((spn/cn)))); then
let k++ 
let m=0
f=$(ls "d$k")
fname="d$k/$f"
let per=k*100/cn
echo "--------$per%------------"
fi
done < sunspot.txt
echo "-------Done!--------"	
