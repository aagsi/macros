#!/bin/bash


for xstep in $(seq 60 0.5 74 )
do
        for ystep in $(seq 13 0.1 17)
        do


#for dtheta in $(seq 0 1 1)
#do
#        for dphi in $(seq 0 1 1)
#        do




#echo $xstep $ystep

 sbatch prism_lens_step_scan $xstep $ystep

        done
done

