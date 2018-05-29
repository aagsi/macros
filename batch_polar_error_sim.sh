#!/bin/bash
MAX=150
for ((i=20; i <= MAX ; i+=5)) ; do
for df in $(seq -0.9 0.01 0.9)
do


total=`echo $i + $df | bc`


#echo sbatch sim $total $i

sbatch sim $total $i

done
done
