 #!/bin/bash

DIRs="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/sim/332"
DIRr="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/sim/332"
DIRl="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/lut"
DIRsep="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/prtdirc/build"

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



# proton
./prtdirc -s 2 -i "${DIRs}/new_theta_20_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/20_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_30_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/30_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_40_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/40_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_50_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/50_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_60_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/60_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_70_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/70_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_80_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/80_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_90_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/90_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_100_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/100_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_110_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/110_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_120_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/120_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_130_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/130_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_140_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/140_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_150_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/150_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton


# pion
./prtdirc -s 2 -i "${DIRs}/new_theta_20_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/20_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_30_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/30_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_40_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/40_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_50_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/50_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_60_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/60_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_70_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/70_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_80_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/80_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_90_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/90_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_100_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/100_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_110_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/110_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_120_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/120_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_130_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/130_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_140_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/140_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_150_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/150_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton

