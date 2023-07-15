#!/bin/bash

num=1
counter=1

for f in *.sam 
do
    mv $f $num"_"$f
    counter=($counter + 1)
if ($counter == 8000)
then
    counter=1
    num=($num + 1)
fi
done