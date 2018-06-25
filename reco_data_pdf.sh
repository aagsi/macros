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




# create cherenkove PDF
#./prtdirc -s 2 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_20_sph_p_data.root" -e 0 -v 2  -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle20_proton -t1 0 -t2 0
#./prtdirc -s 2 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_20_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle20_proton -t1 0 -t2 0
echo before comment
: <<'END'
./prtdirc -s 2 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_30_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle30_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_30_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle30_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_40C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_40_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle40_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_40C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_40_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle40_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_50C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_50_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle50_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_50C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_50_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle50_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_60C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_60_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle60_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_60C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_60_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle60_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_70C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_70_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle70_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_70C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_70_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle70_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_80C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_80_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle80_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_80C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_80_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle80_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_90C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_90_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle90_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_90C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_90_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle90_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_100C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_100_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle100_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_100C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_100_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle100_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_110C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_110_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle110_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_110C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_110_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle110_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_120C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_120_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle120_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_120C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_120_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle120_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_130C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_130_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle130_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_130C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_130_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle130_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_140C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_140_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle140_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_140C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_140_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle140_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_150C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_150_sph_p_data.root"  -e 0 -v 2 -gpdf 1 -recop 1 -recopi 0 -tr 5.0-chcut 5.0 -recoangle $recoAngle150_proton -t1 0 -t2 0
./prtdirc -s 2 -i "${DIRd}/beam_s332_150C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/histo_150_sph_pi_data.root" -e 0 -v 2 -gpdf 1 -recop 0 -recopi 1 -tr 5.0-chcut 5.0 -recoangle $recoAngle150_proton -t1 0 -t2 0

END
echo after comment


# separation power
#./prtdirc -s 3 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_20_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle20_proton -t1 0 -t2 0
echo start separation
#: <<'END'
./prtdirc -s 3 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_30_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle30_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_40C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_40_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle40_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_50C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_50_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle50_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_60C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_60_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle60_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_70C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_70_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle70_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_80C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_80_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle80_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_90C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_90_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle90_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_100C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_100_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle100_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_110C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_110_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle110_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_120C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_120_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle120_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_130C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_130_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle130_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_140C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_140_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle140_proton -t1 0 -t2 0 -gpdf 2
./prtdirc -s 3 -i "${DIRd}/beam_s332_150C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepPDF_150_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle150_proton -t1 0 -t2 0 -gpdf 2

#END

# sep corrected
#./prtdirc -s 3 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_20_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle20_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_30_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle30_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_40C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_40_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle40_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_50C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_50_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle50_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_60C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_60_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle60_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_70C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_70_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle70_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_80C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_80_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle80_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_90C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_90_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle90_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_100C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_100_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle100_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_110C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_110_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle110_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_120C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_120_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle120_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_130C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_130_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle130_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_140C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_140_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle140_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_150C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_corrected_150_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle150_proton -t1 99 -t2 0 -openchcorr 1 -gpdf 3

# sep std w/o correction 
#./prtdirc -s 3 -i "${DIRd}/beam_s332_20C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_20_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle20_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_30C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_30_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle30_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_40C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_40_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle40_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_50C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_50_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle50_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_60C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_60_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle60_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_70C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_70_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle70_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_80C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_80_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle80_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_90C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_90_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle90_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_100C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_100_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle100_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_110C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_110_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle110_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_120C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_120_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle120_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_130C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_130_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle130_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_140C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_140_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle140_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3
#./prtdirc -s 3 -i "${DIRd}/beam_s332_150C.root" -u "${DIRl}/lut_opt_332_2018_cs_avr.root" -o "${DIRr}/sepSTD_150_sph_data.root"  -e $nevents -v 2 -tr 5.0-chcut 5.0 -recoangle $recoAngle150_proton -t1 0 -t2 0 -openchcorr 0 -gpdf 3



