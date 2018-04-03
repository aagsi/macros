 #!/bin/bash



DIRs="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/sim/332"
DIRr="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/sim/332"
DIRl="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/lut"
DIRsep="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/prtdirc/build"

nevents=0
chAngleCut=0.04
timeCut=4
recoAngle_p=0.81676
recoAngle_pi=0.824869

# bar_3sph_grease

#./prtdirc -s 2 -i "${DIRs}/theta_20_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/20_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_30_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/30_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_40_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/40_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_50_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/50_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_60_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/60_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_70_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/70_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_80_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/80_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_90_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/90_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_100_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/100_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_110_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/110_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_120_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/120_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_130_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/130_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_140_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/140_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p
#./prtdirc -s 2 -i "${DIRs}/theta_150_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/150_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_p

#./prtdirc -s 2 -i "${DIRs}/theta_20_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/20_sph_pi_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_30_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/30_sph_pi_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_40_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/40_sph_pi_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_50_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/50_sph_pi_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_60_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/60_sph_pi_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_70_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/70_sph_pi_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_80_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/80_sph_pi_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_90_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/90_sph_pi_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_100_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/100_sph_proton_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_110_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/110_sph_proton_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_120_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/120_sph_proton_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_130_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/130_sph_proton_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_140_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/140_sph_proton_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi
#./prtdirc -s 2 -i "${DIRs}/theta_150_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/150_sph_proton_sim.root" -e $nevents -v 3 -tr $timeCut -chcut $chAngleCut -recoangle $recoAngle_pi




###################################################################################################################################################################################################

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




./prtdirc -s 2 -i "${DIRs}/new_theta_20_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/20_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_30_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/30_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_40_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/40_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_50_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/50_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_60_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/60_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_70_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/70_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_80_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/80_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_90_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/90_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_100_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/100_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_110_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/110_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_120_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/120_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_130_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/130_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_140_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/140_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_150_3lsph_proton_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/150_sph_proton_sim.root" -e $nevents -v 2 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton





./prtdirc -s 2 -i "${DIRs}/new_theta_20_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/20_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_30_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/30_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_40_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/40_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_50_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/50_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_60_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/60_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_70_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/70_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_80_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/80_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_90_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/90_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_100_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/100_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_110_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/110_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_120_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/120_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_130_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/130_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_140_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/140_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton
./prtdirc -s 2 -i "${DIRs}/new_theta_150_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/150_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton



############################################################################################################################################################################################################

# pion cuts

recoAngle20_pi=0.826920
timeCut20_pi=1.863995
chAngleCut20_pi=0.047589
recoAngle30_pi=0.824533
timeCut30_pi=1.901147
chAngleCut30_pi=0.049349
recoAngle40_pi=0.825388
timeCut40_pi=1.877434
chAngleCut40_pi=0.054559
recoAngle50_pi=0.824601
timeCut50_pi=1.880113
chAngleCut50_pi=0.044268
recoAngle60_pi=0.823733
timeCut60_pi=1.886543
chAngleCut60_pi=0.049438
recoAngle70_pi=0.820671
timeCut70_pi=1.931819
chAngleCut70_pi=0.044584
recoAngle80_pi=0.820424
timeCut80_pi=2.088451
chAngleCut80_pi=0.050023
recoAngle90_pi=0.819132
timeCut90_pi=1.794568
chAngleCut90_pi=0.048559
recoAngle100_pi=0.820097
timeCut100_pi=1.395242
chAngleCut100_pi=0.045446
recoAngle110_pi=0.821585
timeCut110_pi=1.371660
chAngleCut110_pi=0.040614
recoAngle120_pi=0.822756
timeCut120_pi=1.385644
chAngleCut120_pi=0.051380
recoAngle130_pi=0.824813
timeCut130_pi=1.403864
chAngleCut130_pi=0.045468
recoAngle140_pi=0.825208
timeCut140_pi=1.440179
chAngleCut140_pi=0.049217
recoAngle150_pi=0.823934
timeCut150_pi=1.521143
chAngleCut150_pi=0.048570

#./prtdirc -s 2 -i "${DIRs}/theta_20_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/20_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut20_pi -chcut $chAngleCut20_pi -recoangle $recoAngle20_pi
#./prtdirc -s 2 -i "${DIRs}/theta_30_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/30_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut30_pi -chcut $chAngleCut30_pi -recoangle $recoAngle30_pi
#./prtdirc -s 2 -i "${DIRs}/theta_40_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/40_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut40_pi -chcut $chAngleCut40_pi -recoangle $recoAngle40_pi
#./prtdirc -s 2 -i "${DIRs}/theta_50_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/50_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut50_pi -chcut $chAngleCut50_pi -recoangle $recoAngle50_pi
#./prtdirc -s 2 -i "${DIRs}/theta_60_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/60_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut60_pi -chcut $chAngleCut60_pi -recoangle $recoAngle60_pi
#./prtdirc -s 2 -i "${DIRs}/theta_70_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/70_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut70_pi -chcut $chAngleCut70_pi -recoangle $recoAngle70_pi
#./prtdirc -s 2 -i "${DIRs}/theta_80_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/80_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut80_pi -chcut $chAngleCut80_pi -recoangle $recoAngle80_pi
#./prtdirc -s 2 -i "${DIRs}/theta_90_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/90_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut90_pi -chcut $chAngleCut90_pi -recoangle $recoAngle90_pi
#./prtdirc -s 2 -i "${DIRs}/theta_100_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/100_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut100_pi -chcut $chAngleCut100_pi -recoangle $recoAngle100_pi
#./prtdirc -s 2 -i "${DIRs}/theta_110_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/110_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut110_pi -chcut $chAngleCut110_pi -recoangle $recoAngle110_pi
#./prtdirc -s 2 -i "${DIRs}/theta_120_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/120_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut120_pi -chcut $chAngleCut120_pi -recoangle $recoAngle120_pi
#./prtdirc -s 2 -i "${DIRs}/theta_130_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/130_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut130_pi -chcut $chAngleCut130_pi -recoangle $recoAngle130_pi
#./prtdirc -s 2 -i "${DIRs}/theta_140_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/140_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut140_pi -chcut $chAngleCut140_pi -recoangle $recoAngle140_pi
#./prtdirc -s 2 -i "${DIRs}/theta_150_3lsph_pi_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/150_sph_pi_sim.root" -e $nevents -v 2 -tr $timeCut150_pi -chcut $chAngleCut150_pi -recoangle $recoAngle150_pi

###############################################################################################################################################################################################################################
# Separation power:


#./prtdirc -s 3 -i "${DIRs}/theta_20_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_20_sph_sim.root" -e $nevents -v 3 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton
#./prtdirc -s 3 -i "${DIRs}/theta_30_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_30_sph_sim.root" -e $nevents -v 3 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton
#./prtdirc -s 3 -i "${DIRs}/theta_40_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_40_sph_sim.root" -e $nevents -v 3 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton
#./prtdirc -s 3 -i "${DIRs}/theta_50_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_50_sph_sim.root" -e $nevents -v 3 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton
#./prtdirc -s 3 -i "${DIRs}/theta_60_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_60_sph_sim.root" -e $nevents -v 3 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton
#./prtdirc -s 3 -i "${DIRs}/theta_70_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_70_sph_sim.root" -e $nevents -v 3 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton
#./prtdirc -s 3 -i "${DIRs}/theta_80_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_80_sph_sim.root" -e $nevents -v 3 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton
#./prtdirc -s 3 -i "${DIRs}/theta_90_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_90_sph_sim.root" -e $nevents -v 3 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton
#./prtdirc -s 3 -i "${DIRs}/theta_100_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_100_sph_sim.root" -e $nevents -v 3 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton
#./prtdirc -s 3 -i "${DIRs}/theta_110_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_110_sph_sim.root" -e $nevents -v 3 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton
#./prtdirc -s 3 -i "${DIRs}/theta_120_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_120_sph_sim.root" -e $nevents -v 3 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton
#./prtdirc -s 3 -i "${DIRs}/theta_130_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_130_sph_sim.root" -e $nevents -v 3 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton
#./prtdirc -s 3 -i "${DIRs}/theta_140_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_140_sph_sim.root" -e $nevents -v 3 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton
#./prtdirc -s 3 -i "${DIRs}/theta_150_3lsph_all_sim.root" -u "${DIRl}/lut_db_332_2018_cs_avr.root" -o "${DIRr}/separation_150_sph_sim.root" -e $nevents -v 3 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton


