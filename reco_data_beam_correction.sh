#!/bin/bash

DIRr="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332"
DIRd="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/bardata/332/data"
DIRl="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/lut"

nevents=500000
chAngleCut=0.04
timeCut=4
recoAngle_p=0.81676
recoAngle_pi=0.824869

#proton cuts

recoAngle20_proton=0.818334
timeCut20_proton=1.842077
chAngleCut20_proton=0.045394
recoAngle30_proton=0.816051
timeCut30_proton=1.910087
chAngleCut30_proton=0.045984
recoAngle40_proton=0.816986
timeCut40_proton=1.889148
chAngleCut40_proton=0.050538
recoAngle50_proton=0.817087
timeCut50_proton=1.892745
chAngleCut50_proton=0.043855
recoAngle60_proton=0.815137
timeCut60_proton=1.895385
chAngleCut60_proton=0.047770
recoAngle70_proton=0.813078
timeCut70_proton=1.960066
chAngleCut70_proton=0.046759
recoAngle80_proton=0.812915
timeCut80_proton=2.091160
chAngleCut80_proton=0.047968
recoAngle90_proton=0.810881
timeCut90_proton=1.828611
chAngleCut90_proton=0.049158
recoAngle100_proton=0.811012
timeCut100_proton=1.395309
chAngleCut100_proton=0.038529
recoAngle110_proton=0.813848
timeCut110_proton=1.384700
chAngleCut110_proton=0.042495
recoAngle120_proton=0.814565
timeCut120_proton=1.377825
chAngleCut120_proton=0.047388
recoAngle130_proton=0.817052
timeCut130_proton=1.407649
chAngleCut130_proton=0.046720
recoAngle140_proton=0.816796
timeCut140_proton=1.441139
chAngleCut140_proton=0.050637
recoAngle150_proton=0.815613
timeCut150_proton=1.527406
chAngleCut150_proton=0.045013

# create spr
./prtdirc -s 2 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/30_test_p_data.root" -recop 1 -recopi 0  -e $nevents -v 2 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/30_test_pi_data.root" -recop 0 -recopi 1  -e $nevents -v 2 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton -t1 0 -t2 0

# create spr
#./prtdirc -s 2 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/test_open1.root" -recop 1 -recopi 0  -e 500 -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_pro$  -openchcorr 1
#./prtdirc -s 2 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/test_open2.root" -recop 0 -recopi 1  -e 500 -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_pr$ -openchcorr 2



