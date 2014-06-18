#!/bin/sh

# $x input file

# Get tag
TAG=`echo $1 | cut -d'_' -f 1`
echo Cooking tag \"$TAG\"
echo Output file \"cook.dat\"
rm aux.dat
touch aux.dat

for arg in "$@"
do
	threads=`echo $arg | cut -d'_' -f 3`
	threads=`echo $threads | cut -d'.' -f 1`
	echo Cooking $arg, $threads threads...
	cat $arg | awk -v THR=$threads '{sum[$1]+=$2; count[$1]++;} 
					END {
						for (x in sum) 
							print x " " THR " " sum[x]/count[x];
					}' >> aux.dat	
done

sort -s -n -k 2,2 aux.dat > cook.dat
rm aux.dat
