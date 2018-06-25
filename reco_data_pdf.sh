#!/bin/bash

DIRr="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332"
DIRd="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/bardata/332/data"
DIRl="/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/lut"

nevents=500000
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


timeCut20_proton=5
#timeCut20_proton=1.862127
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

chAngleCut20_proton=5
#chAngleCut20_proton=0.038775
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

# create cherenkove PDF
#./prtdirc -s 2 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -recop 1 -recopi 0 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton -t1 0 -t2 0
#./prtdirc -s 2 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton -t1 0 -t2 0
echo before comment
: <<'END'
./prtdirc -s 2 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_30_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_30_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_40C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_40_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_40C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_40_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_50C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_50_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_50C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_50_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_60C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_60_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_60C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_60_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_70C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_70_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_70C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_70_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_80C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_80_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_80C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_80_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_90C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_90_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_90C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_90_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_100C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_100_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton -t1 0 -t2 0 
./prtdirc -s 2 -i "${DIRd}/beam_s332_100C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_100_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_110C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_110_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_110C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_110_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_120C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_120_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_120C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_120_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_130C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_130_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_130C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_130_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_140C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_140_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_140C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_140_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_150C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_150_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_150C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_150_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton -t1 0 -t2 0

END
echo after comment


# separation power
#./prtdirc -s 3 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_20_sph_data.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton -t1 0 -t2 0
echo start separation
: <<'END'
./prtdirc -s 3 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_30_sph_data.root"  -e $nevents -v 2 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_40C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_40_sph_data.root"  -e $nevents -v 2 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_50C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_50_sph_data.root"  -e $nevents -v 2 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_60C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_60_sph_data.root"  -e $nevents -v 2 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_70C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_70_sph_data.root"  -e $nevents -v 2 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_80C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_80_sph_data.root"  -e $nevents -v 2 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_90C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_90_sph_data.root"  -e $nevents -v 2 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_100C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_100_sph_data.root"  -e $nevents -v 2 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_110C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_110_sph_data.root"  -e $nevents -v 2 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_120C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_120_sph_data.root"  -e $nevents -v 2 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_130C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_130_sph_data.root"  -e $nevents -v 2 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_140C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_140_sph_data.root"  -e $nevents -v 2 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_150C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_150_sph_data.root"  -e $nevents -v 2 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton -t1 0 -t2 0 -gpdf 2

END

# sep corrected
./prtdirc -s 3 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_20_sph_data.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_30_sph_data.root"  -e $nevents -v 2 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_40C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_40_sph_data.root"  -e $nevents -v 2 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_50C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_50_sph_data.root"  -e $nevents -v 2 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_60C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_60_sph_data.root"  -e $nevents -v 2 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_70C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_70_sph_data.root"  -e $nevents -v 2 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_80C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_80_sph_data.root"  -e $nevents -v 2 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_90C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_90_sph_data.root"  -e $nevents -v 2 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_100C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_100_sph_data.root"  -e $nevents -v 2 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_110C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_110_sph_data.root"  -e $nevents -v 2 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_120C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_120_sph_data.root"  -e $nevents -v 2 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_130C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_130_sph_data.root"  -e $nevents -v 2 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_140C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_140_sph_data.root"  -e $nevents -v 2 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_150C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_150_sph_data.root"  -e $nevents -v 2 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3

# sep std w/o correction 
#./prtdirc -s 3 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_20_sph_data.root"  -e $nevents -v 2 -tr $timeCut20_proton -chcut $chAngleCut20_proton -recoangle $recoAngle20_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_30_sph_data.root"  -e $nevents -v 2 -tr $timeCut30_proton -chcut $chAngleCut30_proton -recoangle $recoAngle30_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_40C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_40_sph_data.root"  -e $nevents -v 2 -tr $timeCut40_proton -chcut $chAngleCut40_proton -recoangle $recoAngle40_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_50C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_50_sph_data.root"  -e $nevents -v 2 -tr $timeCut50_proton -chcut $chAngleCut50_proton -recoangle $recoAngle50_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_60C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_60_sph_data.root"  -e $nevents -v 2 -tr $timeCut60_proton -chcut $chAngleCut60_proton -recoangle $recoAngle60_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_70C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_70_sph_data.root"  -e $nevents -v 2 -tr $timeCut70_proton -chcut $chAngleCut70_proton -recoangle $recoAngle70_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_80C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_80_sph_data.root"  -e $nevents -v 2 -tr $timeCut80_proton -chcut $chAngleCut80_proton -recoangle $recoAngle80_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_90C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_90_sph_data.root"  -e $nevents -v 2 -tr $timeCut90_proton -chcut $chAngleCut90_proton -recoangle $recoAngle90_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_100C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_100_sph_data.root"  -e $nevents -v 2 -tr $timeCut100_proton -chcut $chAngleCut100_proton -recoangle $recoAngle100_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_110C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_110_sph_data.root"  -e $nevents -v 2 -tr $timeCut110_proton -chcut $chAngleCut110_proton -recoangle $recoAngle110_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_120C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_120_sph_data.root"  -e $nevents -v 2 -tr $timeCut120_proton -chcut $chAngleCut120_proton -recoangle $recoAngle120_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_130C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_130_sph_data.root"  -e $nevents -v 2 -tr $timeCut130_proton -chcut $chAngleCut130_proton -recoangle $recoAngle130_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_140C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_140_sph_data.root"  -e $nevents -v 2 -tr $timeCut140_proton -chcut $chAngleCut140_proton -recoangle $recoAngle140_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_150C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_150_sph_data.root"  -e $nevents -v 2 -tr $timeCut150_proton -chcut $chAngleCut150_proton -recoangle $recoAngle150_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3



