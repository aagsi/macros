#!/bin/bash

for t in $(seq 20 5 150)
do

  sbatch sim $t

done


