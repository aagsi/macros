#!/bin/bash


for dtheta in $(seq -0.02 0.001 0.02)
do
        for dphi in $(seq -0.01 0.0005 0.01)
        do


#for dtheta in $(seq 0 1 1)
#do
#        for dphi in $(seq 0 1 1)
#        do




#echo $dtheta $dphi

 sbatch beam_correction $dtheta $dphi

        done
done

