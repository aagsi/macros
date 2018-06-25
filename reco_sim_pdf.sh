#!/bin/bash

DIRs="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/sim/332"
DIRr="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/sim/332"
DIRl="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/lut"

nevents=0

chAngleCut=0.04
timeCut=4
recoAngle_p=0.81676
recoAngle_pi=0.824869


# new cuts:
#p
recoAngle20_proton=0.817927
recoAngle30_proton=0.816989
recoAngle40_proton=0.816539
recoAngle50_proton=0.815916
recoAngle60_proton=0.814639
recoAngle70_proton=0.813995
recoAngle80_proton=0.813616
recoAngle90_proton=0.813517
recoAngle100_proton=0.813906
recoAngle110_proton=0.814171
recoAngle120_proton=0.814349
recoAngle130_proton=0.815740
recoAngle140_proton=0.816496
recoAngle150_proton=0.816786

timeCut20_proton=1.862127
timeCut30_proton=1.910943
timeCut40_proton=1.878375
timeCut50_proton=1.895868
timeCut60_proton=1.894080
timeCut70_proton=1.966618
timeCut80_proton=2.132993
timeCut90_proton=1.823477
timeCut100_proton=1.402801
timeCut110_proton=1.383924
timeCut120_proton=1.385515
timeCut130_proton=1.412302
timeCut140_proton=1.437529
timeCut150_proton=1.529106

chAngleCut20_proton=0.038775
chAngleCut30_proton=0.041868
chAngleCut40_proton=0.042588
chAngleCut50_proton=0.041037
chAngleCut60_proton=0.045299
chAngleCut70_proton=0.046581
chAngleCut80_proton=0.045421
chAngleCut90_proton=0.054653
chAngleCut100_proton=0.047778
chAngleCut110_proton=0.046473
chAngleCut120_proton=0.044891
chAngleCut130_proton=0.043477
chAngleCut140_proton=0.041921
chAngleCut150_proton=0.042944
##########################
#pi
recoAngle20_pi=0.825962
recoAngle30_pi=0.824957
recoAngle40_pi=0.824360
recoAngle50_pi=0.824033
recoAngle60_pi=0.822690
recoAngle70_pi=0.822330
recoAngle80_pi=0.821902
recoAngle90_pi=0.822180
recoAngle100_pi=0.822104
recoAngle110_pi=0.822211
recoAngle120_pi=0.822613
recoAngle130_pi=0.823945
recoAngle140_pi=0.824355
recoAngle150_pi=0.824832

timeCut20_pi=1.862249
timeCut30_pi=1.898656
timeCut40_pi=1.863541
timeCut50_pi=1.868133
timeCut60_pi=1.878101
timeCut70_pi=1.942607
timeCut80_pi=2.100297
timeCut90_pi=1.815946
timeCut100_pi=1.399277
timeCut110_pi=1.376950
timeCut120_pi=1.381105
timeCut130_pi=1.400709
timeCut140_pi=1.437260
timeCut150_pi=1.524456

chAngleCut20_pi=0.038741
chAngleCut30_pi=0.040075
chAngleCut40_pi=0.040655
chAngleCut50_pi=0.041165
chAngleCut60_pi=0.045630
chAngleCut70_pi=0.046156
chAngleCut80_pi=0.045905
chAngleCut90_pi=0.053783
chAngleCut100_pi=0.046921
chAngleCut110_pi=0.044334
chAngleCut120_pi=0.044868
chAngleCut130_pi=0.042434
chAngleCut140_pi=0.040151
chAngleCut150_pi=0.040652

# Crear PDF to test number of events which generate PDF
: <<'END'
./prtdirc -s 2 -i "${DIRs}/forPDF_1_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_1_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_1_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_1_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_2_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_2_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_2_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_2_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_3_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_3_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_3_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_3_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_4_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_4_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_5_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_5_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_5_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_5_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_6_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_6_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_6_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_6_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_7_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_7_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_7_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_7_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_8_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_8_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_8_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_8_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_9_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_9_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_9_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_9_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_10_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_10_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_10_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_10_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_11_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_11_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_11_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_11_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_12_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_12_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_12_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_12_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_13_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_13_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_13_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_13_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_14_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_14_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_14_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_14_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_15_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_15_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_15_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_15_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_16_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_16_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_16_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_16_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_17_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_17_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_17_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_17_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_18_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_18_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_18_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_18_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_19_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_19_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_19_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_19_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_20_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_20_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_20_theta_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_20_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
END
######
# Test different PDF with Different events number:
: <<'END'
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_1_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 1
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_2_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 2
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_3_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 3
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_4_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 4
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_5_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 5
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_6_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 6
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_7_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 7
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_8_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 8
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_9_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 9
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_10_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 10
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_11_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 11
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_12_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 12
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_13_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 13
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_14_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 14
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_15_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 15
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_16_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 16
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_17_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 17
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_18_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 18
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_19_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 19
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/pdf_neventTest_20_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 20
END

########################################################
#creat PDF
: <<'END'
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_20_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_20_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0
#: <<'END'
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_30_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_30_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_30_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_30_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_40_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_40_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_40_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_40_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_50_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_50_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_50_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_50_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_60_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_60_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_60_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_60_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_70_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_70_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_70_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_70_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_80_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_80_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_80_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_80_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_90_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_90_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_90_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_90_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_100_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_100_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_100_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_100_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_110_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_110_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_110_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_110_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_120_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_120_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_120_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_120_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_130_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_130_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_130_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_130_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_140_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_140_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_140_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_140_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton  -recop 1 -recopi 0
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_150_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_150_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton  -recop 0 -recopi 1
./prtdirc -s 2 -i "${DIRs}/forPDF_theta_150_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_sim_4BarRefl_150_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton  -recop 1 -recopi 0
END

: <<'END'
# separation power
./prtdirc -s 3 -i "${DIRs}/new_20_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_20_sph.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -openchcorr 0
#: <<'END'
./prtdirc -s 3 -i "${DIRs}/new_30_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_30_sph.root"  -e $nevents -v 2 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_40_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_40_sph.root"  -e $nevents -v 2 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_50_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_50_sph.root"  -e $nevents -v 2 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_60_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_60_sph.root"  -e $nevents -v 2 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_70_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_70_sph.root"  -e $nevents -v 2 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_80_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_80_sph.root"  -e $nevents -v 2 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_90_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_90_sph.root"  -e $nevents -v 2 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_100_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_100_sph.root"  -e $nevents -v 2 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_110_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_110_sph.root"  -e $nevents -v 2 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_120_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_120_sph.root"  -e $nevents -v 2 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_130_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_130_sph.root"  -e $nevents -v 2 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_140_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_140_sph.root"  -e $nevents -v 2 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton  -openchcorr 0
./prtdirc -s 3 -i "${DIRs}/new_150_3lsph_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sep_std_sim_150_sph.root"  -e $nevents -v 2 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton  -openchcorr 0
END




#spr
#./prtdirc -s 2 -i "${DIRs}/new_theta_20_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_20_sph_pi_sim.root" -e 0 -v 2  -gpdf 0 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_pi  -recop 0 -recopi 1 -openchcorr 2
./prtdirc -s 2 -i "${DIRs}/new_theta_20_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_20_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_30_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_30_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_30_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_30_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_40_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_40_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_40_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_40_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_50_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_50_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_50_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_50_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_60_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_60_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_60_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_60_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_70_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_70_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_70_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_70_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_80_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_80_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_80_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_80_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_90_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_90_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_90_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_90_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_100_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_100_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_100_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_100_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_110_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_110_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_110_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_110_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_120_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_120_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_pi -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_120_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_120_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_130_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_130_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_130_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_130_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_140_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_140_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_140_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_140_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton  -recop 1 -recopi 0 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_150_3lsph_pi_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_150_sph_pi_sim.root"  -e 0 -v 2 -gpdf 0 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_pi  -recop 0 -recopi 1 -openchcorr 2
#./prtdirc -s 2 -i "${DIRs}/new_theta_150_3lsph_proton_sim.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/spr_wt_150_sph_p_sim.root" -e 0 -v 2 -gpdf 0 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton  -recop 1 -recopi 0 -openchcorr 2





