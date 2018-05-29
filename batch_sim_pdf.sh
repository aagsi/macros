#!/bin/bash

for t in $(seq 100 1 122)
do

  sbatch sim_pdf $t

done


