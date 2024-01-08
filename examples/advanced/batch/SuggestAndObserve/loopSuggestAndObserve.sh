#!/bin/bash

rm -f cache.txt cache[1-9]*.txt x*.txt f*.txt param[1-9]*.txt

i=0
until [ $i -gt 10 ]
do
    ./suggest.exe cache$i.txt param$i.txt > x$i.txt
    nbPoints=`wc -l x$i.txt | awk '{print $1}'`
    echo "Iteration $i generated $nbPoints points"
    if [ 0 -eq $nbPoints ]
    then
        break
    fi
    echo "Evaluate $nbPoints points"
    ./bbr.exe x$i.txt > f$i.txt
    j=$((i+1))
    echo "Observe $nbPoints points"
    ./observe.exe x$i.txt f$i.txt cache$i.txt cache$j.txt param$i.txt > param$j.txt
    i=$j
done
echo "Normal termination"
