#!/bin/bash
MAX=150
for ((i=20; i <= MAX ; i+=5)) ; do
for df in $(seq 0 1 100)
do

#echo sbatch sim $i $df 

sbatch sim $i $df

done
done
