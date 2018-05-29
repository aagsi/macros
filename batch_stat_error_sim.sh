#!/bin/bash
recoAngle_proton[0]=0.817927
recoAngle_proton[1]=0.816989
recoAngle_proton[2]=0.816539
recoAngle_proton[3]=0.815916
recoAngle_proton[4]=0.814639
recoAngle_proton[5]=0.813995
recoAngle_proton[6]=0.813616
recoAngle_proton[7]=0.813517
recoAngle_proton[8]=0.813906
recoAngle_proton[9]=0.814171
recoAngle_proton[10]=0.814349
recoAngle_proton[11]=0.815740
recoAngle_proton[12]=0.816496
recoAngle_proton[13]=0.816786

recoAngle_pi[0]=0.825962
recoAngle_pi[1]=0.824957
recoAngle_pi[2]=0.824360
recoAngle_pi[3]=0.824033
recoAngle_pi[4]=0.822690
recoAngle_pi[5]=0.822330
recoAngle_pi[6]=0.821902
recoAngle_pi[7]=0.822180
recoAngle_pi[8]=0.822104
recoAngle_pi[9]=0.822211
recoAngle_pi[10]=0.822613
recoAngle_pi[11]=0.823945
recoAngle_pi[12]=0.824355
recoAngle_pi[13]=0.824832

COUNTER=0

MAX=150
for ((i=20; i <= MAX ; i+=5)) ; do
for df in $(seq 0 1 100)
do

echo $i ${recoAngle_proton[$COUNTER]} ${recoAngle_pi[$COUNTER]} $df

#sbatch $i ${recoAngle_proton[$COUNTER]} ${recoAngle_pi[$COUNTER]} $df

done
COUNTER=$[$COUNTER +1]
done

