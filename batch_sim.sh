#!/bin/bash

for t in $(seq 20 10 150)
do


#echo $t

  sbatch sim $t

done


