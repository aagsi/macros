 //2018  ok// ----------w-------------------------------
// PrtLutReco.cpp
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtLutReco.h"
#include "PrtManager.h"
#include "PrtLutNode.h"
#include "PrtTrackInfo.h"
#include "PrtPhotonInfo.h"
#include "PrtAmbiguityInfo.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRotation.h"
#include "TGraph.h"
#include <TVirtualFitter.h>
#include <TArc.h>
#include <TLegend.h>
#define prt__sim
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"
#include "TGraph.h"
#include "TMultiGraph.h"
using std::cout;
using std::endl;

TH1F*  fHist0 = new TH1F("timediff",";t_{calc}-t_{measured} [ns];entries [#]", 500,-10,10);
TH1F*  fHist0i = new TH1F("timediffi",";t_{calc}-t_{measured} [ns];entries [#]", 500,-10,10);
TH1F*  fHist0i_bg = new TH1F("timediffi_bg",";t_{calc}-t_{measured} [ns];entries [#]", 500,-10,10);
//TH1F*  fHist1 = new TH1F("time1",";measured time [ns];entries [#]",   1000,-500,500);  // test
//TH1F*  fHist2 = new TH1F("time2",";calculated time [ns];entries [#]", 1000,-500,500); //test

TH1F*  fHist1 = new TH1F("time1",";measured time [ns];entries [#]",   500,0,50);  // test
TH1F*  fHist2 = new TH1F("time2",";calculated time [ns];entries [#]", 500,0,50); //test

TH2F*  fHist3 = new TH2F("time3",";calculated time [ns];measured time [ns]", 500,0,80, 500,0,40);
TH2F*  fHist4 = new TH2F("time4",";#theta_{c}sin(#varphi_{c});#theta_{c}cos(#varphi_{c}", 100,-1,1, 100,-1,1);
TH2F*  fHist5 = new TH2F("time5",";#theta_{c}sin(#varphi_{c});#theta_{c}cos(#varphi_{c}", 100,-1,1, 100,-1,1);
TH1F*  falpha = new TH1F("alpha",";|#alpha|[degree];entries [#]",   1000,-360,360);
TH1F*  falphai = new TH1F("alphai",";|#alphaf| [degree];entries [#]",   1000,-360,360);
TH1F*  fHistPhotonEnergy = new TH1F("fHistPhotonEnergy",";|#alpha|[degree];entries [#]", 200, 0, 8);

TH1F*  fnHits = new TH1F("fnHits",";number of photons per track [#];entries [#]",   100,0,200);
TH1F*  fnHits_p = new TH1F("fnHits_p",";number of photons per proton track [#];entries [#]",   100,0,200);
TH1F*  fnHits_p_good = new TH1F("fnHits_p_good",";number of Good photons per proton track [#];entries [#]",   100,0,200);

TH1F*  nHits_dac = new TH1F("nHits_dac",";number of photons solutions per proton track [#];entries [#]",   100,0,200);
TH1F*  nHits_dac_syscut_p = new TH1F("nHits_dac_syscut_p",";number of photons per proton track after sys cut [#];entries [#]",   100,0,200);
TH1F*  fnHits_true_sim = new TH1F("fnHits_true_sim",";number of photons per track [#];entries [#]",   100,0,200);

TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",200,-600,600);
TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",200,-600,600);

//TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",1000,-1000,1000);
//TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",1000,-1000,1000);


TH1F*  test_hist = new TH1F("test_hist",";t_{calc}-t_{measured} [ns];entries [#]", 500,-1000,1000);

TF1 *gF1 = new TF1("gaus0","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
TF1 *gF2= new TF1("gaus0","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
Int_t gg_i(0), gg_ind(0);
TGraph gg_gr;
PrtLutNode *fLutNode[5000];
TH1F*  fHistMcp[12]; // changed
TH1F*  fHistMcp_same_path[12]; // changed
TH1F*  fHistCh[960], *fHistCh_read_p[960], *fHistCh_read_pi[960], *hist_nph_wtc_p, *hist_nph_wtc_pi;
TGraph *fHistCh_graph_p[960], *fHistCh_graph_pi[960];

// tof
Double_t TOT_tof1, TOT_tof2, TOT_trg1, TOT_trg2, TOT_trgmzH, TOT_trgmzV, delta_tof2tof1;
//tof2tof1
Int_t bin_d_tof2tof1=100; //400;
Double_t x_d_tof2tof1=30.7; // 30;//25;//65
Double_t y_d_tof2tof1=34; // 35;//85;//75



TAxis *axis_data[960], *axis_data_pi[960];
Int_t bmin_data[960],  bmax_data[960], bmin_data_pi[960], bmax_data_pi[960];
Double_t integral_data[960], integral_data_pi[960];
Double_t xmin_data = 0.6;
Double_t xmax_data = 1.0;

/*
 // LE
 TH1F * htof1_le =       new TH1F("htof1_le","",400, -500, 500);
 TH1F * htof2_le =       new TH1F("htof2_le","",300, -275, -125); // 150
 TH1F * htrg1_le =       new TH1F("htrg1_le","",100, -160, -110); //50
 TH1F * htrg2_le =       new TH1F("htrg2_le","",400, -500, 500);
 TH1F * htrgmzH_le =      new TH1F("htrgmzH_le","",400, -500, 500);
 TH1F * htrgmzV_le =      new TH1F("htrgmzV_le","",400, -500, 500);
 // TOT histo
 TH1F * htof1_tot =      new TH1F("htof1_tot","",1000, 110, 150);
 TH1F * htof2_tot =      new TH1F("htof2_tot","",1000, 110, 150);
 TH1F * htrg1_tot =      new TH1F("htrg1_tot","",100, 110, 125);
 TH1F * htrg2_tot =      new TH1F("htrg2_tot","",1000, 110, 150);
 TH1F * htrgmzH_tot =     new TH1F("htrgmzH_tot","",1000, 110, 150);
 TH1F * htrgmzV_tot =     new TH1F("htrgmzV_tot","",1000, 110, 150);
 */

// delta
TH1F *hdelta_tof2tof1= new TH1F("hdelta_tof2tof1","",bin_d_tof2tof1, x_d_tof2tof1, y_d_tof2tof1);
TH1F *hdelta_tof2tof1_isproton= new TH1F("hdelta_tof2tof1_isproton","",bin_d_tof2tof1, x_d_tof2tof1, y_d_tof2tof1);


TH1F *histo_photon_ambiguity_wo= new TH1F("histo_photon_ambiguity_wo","",100, 0, 100);
TH1F *histo_photon_ambiguity_wt= new TH1F("histo_photon_ambiguity_wt","",100, 0, 100);
TH1F *histo_photon_ambiguity_wtc= new TH1F("histo_photon_ambiguity_wtc","",100, 0, 100);



/*
 // hodo
 const int nHodoPixel=16;
 TH2F * hHodo= new TH2F("hHodo","",16,0,16,16,0,16); // 9 sep 2017
 TH1F * hHodoPixelH[nHodoPixel],hHodoPixelV[nHodoPixel] ;
 */
TH1D * hodoV = new TH1D("1","vertical fibers",16,0,16);
TH1D * hodoH = new TH1D("2","horizontal fibers",16,0,16);
TH2D * hodoF = new TH2D("3","4",16,0,16,16,0,16);

TH1D * hodoV_select = new TH1D("1_select","vertical fibers selection",16,0,16);
TH1D * hodoH_select = new TH1D("2_select","horizontal fibers selection ",16,0,16);
TH2D *hodo_afterCut = new TH2D("5","6",16,0,16,16,0,16);

TH2D *hodo_multi_withmedHfiber = new TH2D("7","8",16,0,16,16,0,16);
TH2D *hodo_multi_withmedVfiber = new TH2D("9","10",16,0,16,16,0,16);

// 2D histo
// LE TOT histo
TH2F *htof1_le_tot = new TH2F("htof1_le_tot","",500, -300, -50, 100, 10, 60);
TH2F *htof2_le_tot = new TH2F("htof2_le_tot","",500, -300, -50, 100, 10, 60);
TH2F *htrg1_le_tot = new TH2F("htrg1_le_tot","",500, -300, -50, 200, 105, 126);
TH2F *htrg2_le_tot = new TH2F("htrg2_le_tot","",500, -300, -50, 200, 125, 140);
TH2F *htrgmzH_le_tot = new TH2F("htrgmzH_le_tot","",500,-300, -50, 200, 128, 145);
TH2F *htrgmzV_le_tot = new TH2F("htrgmzV_le_tot","",500,-300, -50, 200, 133, 150);

// 1D histo
TH1F * countmulti_tof1 =        new TH1F("countmulti_tof1","",7, 0, 7);
TH1F * countmulti_tof2 =        new TH1F("countmulti_tof2","",7, 0, 7);
TH1F * countmulti_trg1 =        new TH1F("countmulti_trg1","",7, 0, 7);
TH1F * countmulti_trg2 =        new TH1F("countmulti_trg2","",7, 0, 7);
TH1F * countmulti_trgmzH =       new TH1F("countmulti_trgmzH","",7, 0, 7);
TH1F * countmulti_trgmzV =       new TH1F("countmulti_trgmzV","",7, 0, 7);

TH1F * countmulti_hodoV =        new TH1F("countmulti_hodoV","",7, 0, 7);
TH1F * countmulti_hodoH =        new TH1F("countmulti_hodoH","",7, 0, 7);

// -----   Default constructor   -------------------------------------------
PrtLutReco::PrtLutReco(TString infile, TString lutfile, Int_t verbose) {
    fVerbose = verbose;
    fChain = new TChain("data");
    fChain->Add(infile);
    fChain->SetBranchAddress("PrtEvent", &fEvent);
    fFile = new TFile(lutfile);
    fTree=(TTree *) fFile->Get("prtlut") ;
    fLut = new TClonesArray("PrtLutNode");
    fTree->SetBranchAddress("LUT",&fLut);
    fTree->GetEntry(0);
    fHist = new TH1F("fHist",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80     // warning
    fHist_copy = new TH1F("fHist_copy",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80
    fHist_correction = new TH1F("fHist_correction",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80
    
    fHist_same_path = new TH1F("fHist_same_path",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80
    fHist_same_path_wotc = new TH1F("fHist_same_path_wotc",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80
    
    
    fHist_bg = new TH1F("fHist_bg",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80
    fHistPi = new TH1F("chrenkov_angle_hist_Pi",  "chrenkov angle pi;#theta_{C} [rad];entries [#]", 80,0.6,1); //150
    fHisti = new TH1F("chrenkov_angle_histi","chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150
    
    
    fFit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
    fSpect = new TSpectrum(10);
    if(infile.Contains("beam_")) {
        TString fileid(infile);
        fileid.Remove(0,fileid.Last('/')+1);
        fileid.Remove(fileid.Last('.')-1);
        prt_data_info = getDataInfo(fileid);
        TString opath(infile);
        opath.Remove(opath.Last('/'));
        if(infile.Contains("C.root")) {
            prt_savepath = opath+Form("/%dr/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());
        } else {
            prt_savepath = opath+Form("/%ds/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());
        }
    } else prt_savepath="data/sim";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    for(Int_t i=0; i<5000; i++) {
        fLutNode[i] = (PrtLutNode*) fLut->At(i);
    }
    cout << "-I- PrtLutReco: Intialization successfull" << endl;
    for(Int_t i=0; i<prt_nmcp; i++) {
        fHistMcp[i] = new TH1F(Form("fHistMcp_%d",i),Form("fHistMcp_%d;#theta_{C} [rad];entries [#]",i), 80,0.6,1); //150
    }
    
    for(Int_t i=0; i<prt_nmcp; i++) {
        fHistMcp_same_path[i] = new TH1F(Form("fHistMcp_same_path_%d",i),Form("fHistMcp_same_path_%d;#theta_{C} [rad];entries [#]",i), 80,0.6,1); //150
    }
    
    for(Int_t i=0; i<960; i++) {
        fHistCh[i] = new TH1F(Form("fHistCh_%d",i),Form("fHistCh_%d;#theta_{C} [rad];entries [#]",i), 2000,0.6,1); //150
    }
    
    
    /*
     // hodo
     for(Int_t p=0; p<nHodoPixel; p++){
     hHodoPixelH[p]  = new TH1F(Form("hHodoPixelH_%d",p),Form("pixelH %d", p),400,0,25);
     hHodoPixelV[p]  = new TH1F(Form("hHodoPixelV_%d",p),Form("pixelV %d", p),400,0,25);
     }
     */
    gr_pi= new TGraph();
    gr_p= new TGraph();
    mg = new TMultiGraph();
}

// -----   Destructor   ----------------------------------------------------
PrtLutReco::~PrtLutReco() {
    
}

Int_t mcpdata[12][65]; // changed
Int_t cluster[12][65]; // changed
Int_t lneighbours[65];
Int_t lsize(0);
Int_t getneighbours(Int_t m, Int_t p) {
    for(Int_t i=0; i<65; i++) if(p==lneighbours[i]) return -1;
    lneighbours[lsize]=p;
    lsize++;
    for(Int_t t=0; t<65; t++) {
        if(mcpdata[m][t]) {
            for(Int_t i=0; i<65; i++) if(t==lneighbours[i]) continue;
            if((t==p-1 && p%8!=0) || (t==p+1 && p%8!=7) ||
               (t==p+8 && p<57) || (t==p-8 && p>8)) getneighbours(m,t);
        }
    }
    return lsize;
}

void getclusters() {
    for(Int_t m=0; m<prt_nmcp; m++) {
        for(Int_t p=0; p<65; p++) {
            if(mcpdata[m][p])  cluster[m][p] = getneighbours(m,p);
            lsize=0;
            for(Int_t i=0; i<65; i++) lneighbours[i]=0;
        }
    }
}

//-------------- Loop over tracks ------------------------------------------
void PrtLutReco::Run(Int_t start, Int_t end) {
    
    gr_pi->SetPoint(1,0.6025,0.0064464);
    gr_pi->SetPoint(2,0.6075,0.00700676);
    gr_pi->SetPoint(3,0.6125,0.00614543);
    gr_pi->SetPoint(4,0.6175,0.00712542);
    gr_pi->SetPoint(5,0.6225,0.00800538);
    gr_pi->SetPoint(6,0.6275,0.00721376);
    gr_pi->SetPoint(7,0.6325,0.00853717);
    gr_pi->SetPoint(8,0.6375,0.00892778);
    gr_pi->SetPoint(9,0.6425,0.00941106);
    gr_pi->SetPoint(10,0.6475,0.00948425);
    gr_pi->SetPoint(11,0.6525,0.00931406);
    gr_pi->SetPoint(12,0.6575,0.010495);
    gr_pi->SetPoint(13,0.6625,0.009709);
    gr_pi->SetPoint(14,0.6675,0.00837867);
    gr_pi->SetPoint(15,0.6725,0.0106539);
    gr_pi->SetPoint(16,0.6775,0.0102642);
    gr_pi->SetPoint(17,0.6825,0.0099212);
    gr_pi->SetPoint(18,0.6875,0.010896);
    gr_pi->SetPoint(19,0.6925,0.0132535);
    gr_pi->SetPoint(20,0.6975,0.0134423);
    gr_pi->SetPoint(21,0.7025,0.0102248);
    gr_pi->SetPoint(22,0.7075,0.011624);
    gr_pi->SetPoint(23,0.7125,0.0117742);
    gr_pi->SetPoint(24,0.7175,0.0123671);
    gr_pi->SetPoint(25,0.7225,0.0106206);
    gr_pi->SetPoint(26,0.7275,0.0118717);
    gr_pi->SetPoint(27,0.7325,0.0127308);
    gr_pi->SetPoint(28,0.7375,0.0112108);
    gr_pi->SetPoint(29,0.7425,0.0125052);
    gr_pi->SetPoint(30,0.7475,0.0133769);
    gr_pi->SetPoint(31,0.7525,0.0119002);
    gr_pi->SetPoint(32,0.7575,0.0144098);
    gr_pi->SetPoint(33,0.7625,0.012639);
    gr_pi->SetPoint(34,0.7675,0.0142283);
    gr_pi->SetPoint(35,0.7725,0.0134882);
    gr_pi->SetPoint(36,0.7775,0.0144453);
    gr_pi->SetPoint(37,0.7825,0.0141313);
    gr_pi->SetPoint(38,0.7875,0.01544);
    gr_pi->SetPoint(39,0.7925,0.0174524);
    gr_pi->SetPoint(40,0.7975,0.0150364);
    gr_pi->SetPoint(41,0.8025,0.0142699);
    gr_pi->SetPoint(42,0.8075,0.0144128);
    gr_pi->SetPoint(43,0.8125,0.019305);
    gr_pi->SetPoint(44,0.8175,0.0274632);
    gr_pi->SetPoint(45,0.8225,0.0316876);
    gr_pi->SetPoint(46,0.8275,0.0391829);
    gr_pi->SetPoint(47,0.8325,0.0268998);
    gr_pi->SetPoint(48,0.8375,0.0244085);
    gr_pi->SetPoint(49,0.8425,0.0166811);
    gr_pi->SetPoint(50,0.8475,0.0158432);
    gr_pi->SetPoint(51,0.8525,0.015679);
    gr_pi->SetPoint(52,0.8575,0.015123);
    gr_pi->SetPoint(53,0.8625,0.0141205);
    gr_pi->SetPoint(54,0.8675,0.014214);
    gr_pi->SetPoint(55,0.8725,0.014544);
    gr_pi->SetPoint(56,0.8775,0.0152295);
    gr_pi->SetPoint(57,0.8825,0.014638);
    gr_pi->SetPoint(58,0.8875,0.0118751);
    gr_pi->SetPoint(59,0.8925,0.011546);
    gr_pi->SetPoint(60,0.8975,0.010731);
    gr_pi->SetPoint(61,0.9025,0.0114815);
    gr_pi->SetPoint(62,0.9075,0.0105145);
    gr_pi->SetPoint(63,0.9125,0.0113845);
    gr_pi->SetPoint(64,0.9175,0.010003);
    gr_pi->SetPoint(65,0.9225,0.0121583);
    gr_pi->SetPoint(66,0.9275,0.0110632);
    gr_pi->SetPoint(67,0.9325,0.0116508);
    gr_pi->SetPoint(68,0.9375,0.0116768);
    gr_pi->SetPoint(69,0.9425,0.0123541);
    gr_pi->SetPoint(70,0.9475,0.0126208);
    gr_pi->SetPoint(71,0.9525,0.0123918);
    gr_pi->SetPoint(72,0.9575,0.0125156);
    gr_pi->SetPoint(73,0.9625,0.00817254);
    gr_pi->SetPoint(74,0.9675,0.00723238);
    gr_pi->SetPoint(75,0.9725,0.00602547);
    gr_pi->SetPoint(76,0.9775,0.00545904);
    gr_pi->SetPoint(77,0.9825,0.00580938);
    gr_pi->SetPoint(78,0.9875,0.00485017);
    gr_pi->SetPoint(79,0.9925,0.0046081);
    gr_pi->SetPoint(80,0.9975,0.00406548);
    gr_p->SetPoint(1,0.6025,0.00663638);
    gr_p->SetPoint(2,0.6075,0.00767694);
    gr_p->SetPoint(3,0.6125,0.00663401);
    gr_p->SetPoint(4,0.6175,0.00727819);
    gr_p->SetPoint(5,0.6225,0.00830735);
    gr_p->SetPoint(6,0.6275,0.00720793);
    gr_p->SetPoint(7,0.6325,0.00808756);
    gr_p->SetPoint(8,0.6375,0.00962893);
    gr_p->SetPoint(9,0.6425,0.00980979);
    gr_p->SetPoint(10,0.6475,0.00981027);
    gr_p->SetPoint(11,0.6525,0.00961041);
    gr_p->SetPoint(12,0.6575,0.0105019);
    gr_p->SetPoint(13,0.6625,0.0101041);
    gr_p->SetPoint(14,0.6675,0.00929094);
    gr_p->SetPoint(15,0.6725,0.0108067);
    gr_p->SetPoint(16,0.6775,0.0113412);
    gr_p->SetPoint(17,0.6825,0.0104426);
    gr_p->SetPoint(18,0.6875,0.0116232);
    gr_p->SetPoint(19,0.6925,0.0131384);
    gr_p->SetPoint(20,0.6975,0.0133079);
    gr_p->SetPoint(21,0.7025,0.0106382);
    gr_p->SetPoint(22,0.7075,0.0108451);
    gr_p->SetPoint(23,0.7125,0.0109121);
    gr_p->SetPoint(24,0.7175,0.0120775);
    gr_p->SetPoint(25,0.7225,0.0101971);
    gr_p->SetPoint(26,0.7275,0.0104782);
    gr_p->SetPoint(27,0.7325,0.0113431);
    gr_p->SetPoint(28,0.7375,0.0112828);
    gr_p->SetPoint(29,0.7425,0.0136397);
    gr_p->SetPoint(30,0.7475,0.0142222);
    gr_p->SetPoint(31,0.7525,0.0127454);
    gr_p->SetPoint(32,0.7575,0.0148459);
    gr_p->SetPoint(33,0.7625,0.0127895);
    gr_p->SetPoint(34,0.7675,0.0152319);
    gr_p->SetPoint(35,0.7725,0.014356);
    gr_p->SetPoint(36,0.7775,0.0145868);
    gr_p->SetPoint(37,0.7825,0.014113);
    gr_p->SetPoint(38,0.7875,0.0155499);
    gr_p->SetPoint(39,0.7925,0.0172176);
    gr_p->SetPoint(40,0.7975,0.0156302);
    gr_p->SetPoint(41,0.8025,0.0188824);
    gr_p->SetPoint(42,0.8075,0.0235316);
    gr_p->SetPoint(43,0.8125,0.0304338);
    gr_p->SetPoint(44,0.8175,0.0368063);
    gr_p->SetPoint(45,0.8225,0.0310192);
    gr_p->SetPoint(46,0.8275,0.028619);
    gr_p->SetPoint(47,0.8325,0.0182543);
    gr_p->SetPoint(48,0.8375,0.0181437);
    gr_p->SetPoint(49,0.8425,0.0142255);
    gr_p->SetPoint(50,0.8475,0.0138021);
    gr_p->SetPoint(51,0.8525,0.0133815);
    gr_p->SetPoint(52,0.8575,0.0133193);
    gr_p->SetPoint(53,0.8625,0.013146);
    gr_p->SetPoint(54,0.8675,0.0126661);
    gr_p->SetPoint(55,0.8725,0.0132804);
    gr_p->SetPoint(56,0.8775,0.013138);
    gr_p->SetPoint(57,0.8825,0.0139312);
    gr_p->SetPoint(58,0.8875,0.0116089);
    gr_p->SetPoint(59,0.8925,0.0112049);
    gr_p->SetPoint(60,0.8975,0.0104976);
    gr_p->SetPoint(61,0.9025,0.0119189);
    gr_p->SetPoint(62,0.9075,0.0111167);
    gr_p->SetPoint(63,0.9125,0.0117907);
    gr_p->SetPoint(64,0.9175,0.0107127);
    gr_p->SetPoint(65,0.9225,0.0124335);
    gr_p->SetPoint(66,0.9275,0.0118183);
    gr_p->SetPoint(67,0.9325,0.0128375);
    gr_p->SetPoint(68,0.9375,0.0133606);
    gr_p->SetPoint(69,0.9425,0.0137546);
    gr_p->SetPoint(70,0.9475,0.0126153);
    gr_p->SetPoint(71,0.9525,0.0111565);
    gr_p->SetPoint(72,0.9575,0.00988622);
    gr_p->SetPoint(73,0.9625,0.00658607);
    gr_p->SetPoint(74,0.9675,0.00646834);
    gr_p->SetPoint(75,0.9725,0.00558539);
    gr_p->SetPoint(76,0.9775,0.00502808);
    gr_p->SetPoint(77,0.9825,0.00540785);
    gr_p->SetPoint(78,0.9875,0.00471952);
    gr_p->SetPoint(79,0.9925,0.00488757);
    gr_p->SetPoint(80,0.9975,0.00407725);
    
    
    TVector3 dird, dir, momInBar(0,0,1),posInBar,cz;
    Double_t mom, cangle,spr,tangle,likelihood(0),boxPhi,weight,evtime,bartime, lenz,dirz,luttheta, barHitTime, hitTime, photonEnergy, alphAngle;
    Int_t  tofPid(0),distPid(0),likePid(0),pdgcode, evpointcount=0;
    Bool_t reflected = kFALSE;
    gStyle->SetOptFit(111);
    TVector3 fnX1 = TVector3 (1,0,0);
    TVector3 fnY1 = TVector3( 0,1,0);
    bool testTrRes = false;
    Double_t angdiv,dtheta,dtphi,prtangle, counter(0);
    
    /////////////////////
    /////////////////////
    /////////////////////
    
    Double_t theta(0),phi(0), trr(0),  nph(0),par1(0), par2(0), par3(0), par4(0), par5(0), par6(0), test1(0), test2(0), test3(0),separation(0),recoAngle(0), chAngleCut(0), timeRes(0);
    Double_t recoP(0), recoPi(0), gPDF(0), openChCorr(0), method_type(-1);
    Int_t solution_number_approach_selection(0),solution_number(0);
    
    chAngleCut = PrtManager::Instance()->GetchAngleCut();
    recoAngle = PrtManager::Instance()->GetrecoAngle();
    timeRes = PrtManager::Instance()->GetTimeRes();
    recoP = PrtManager::Instance()->GetrecoP();
    recoPi = PrtManager::Instance()->GetrecoPi();
    gPDF = PrtManager::Instance()->GetgPDF();
    openChCorr = PrtManager::Instance()->GetopenChCorr();
    method_type = PrtManager::Instance()->GetRunType();
    
    cout<<"@@@@@@@@@@@@ gPDF="<< gPDF << endl;
    cout<<"@@@@@@@@@@@@ recoP="<< recoP << endl;
    cout<<"@@@@@@@@@@@@ recoPi="<< recoPi << endl;
    cout<<"@@@@@@@@@@@@ method_type="<< method_type << endl;
    
    
    TFile *ffile_data_p, *ffile_data_pi;
    Double_t prtangle_pdf, pdf_nph_p, pdf_nph_pi;
    TString cherenkov_data_p_path, cherenkov_data_pi_path;
    
    if (method_type == 3) {
        fChain->GetEntry(20);
        prtangle_pdf = fEvent->GetAngle();
        //cherenkov_data_p_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/ambiguit_pdf/histo_%g_sph_p_data_cherenkovPDF.root", prtangle_pdf);
        //cherenkov_data_pi_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/ambiguit_pdf/histo_%g_sph_pi_data_cherenkovPDF.root", prtangle_pdf);
        //cherenkov_data_p_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/ambiguit_pdf/histo_2BarRefl_%g_sph_p_data_cherenkovPDF.root", prtangle_pdf);
        //cherenkov_data_pi_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/ambiguit_pdf/histo_2BarRefl_%g_sph_pi_data_cherenkovPDF.root", prtangle_pdf);
	cherenkov_data_p_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/ambiguit_pdf/histo_4BarRefl_%g_sph_p_data_cherenkovPDF.root", prtangle_pdf);
        cherenkov_data_pi_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/ambiguit_pdf/histo_4BarRefl_%g_sph_pi_data_cherenkovPDF.root", prtangle_pdf);
	//cherenkov_data_p_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/histo_%d_sph_p_data_cherenkovPDF.root", 40);
        //cherenkov_data_pi_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/histo_%d_sph_pi_data_cherenkovPDF.root", 40);
	cout<<"cherenkov_data_p_path= " <<cherenkov_data_p_path<<endl;
        cout<<"cherenkov_data_pi_path= " <<cherenkov_data_pi_path<<endl;
        ffile_data_p  = new TFile(cherenkov_data_p_path, "READ");
        ffile_data_pi  = new TFile(cherenkov_data_pi_path, "READ");
        for(Int_t pix=0; pix<960; pix++) {
            //fHistCh_graph_p[pix] =new TGraph(   (TH1F*)ffile_data_p->Get(Form("fHistCh_%d",pix))   );
            //i[pix] =new TGraph(   (TH1F*)ffile_data_pi->Get(Form("fHistCh_%d",pix))   );
            fHistCh_read_p[pix] = (TH1F*)ffile_data_p->Get(Form("fHistCh_%d",pix));
            fHistCh_read_pi[pix] = (TH1F*)ffile_data_pi->Get(Form("fHistCh_%d",pix));
        }
        hist_nph_wtc_p=(TH1F*)ffile_data_p->Get("fnHits_p_good");
        hist_nph_wtc_pi=(TH1F*)ffile_data_pi->Get("fnHits_p_good");
        //Double_t pdf_nph_p=42556;
        //Double_t pdf_nph_pi=45334;
        pdf_nph_p=hist_nph_wtc_p->GetEntries();
        pdf_nph_pi=hist_nph_wtc_pi->GetEntries();
        cout<<"@@@@@@@@@@@@ pdf_nph_p="<< pdf_nph_p << endl;
        cout<<"@@@@@@@@@@@@ pdf_nph_pi="<< pdf_nph_pi << endl;
        for(Int_t pix=0; pix<960; pix++) {
            axis_data[pix] = fHistCh_read_p[pix]->GetXaxis();
            if (prtangle_pdf==90)xmin_data  = 0.9;
            if (prtangle_pdf==90)xmax_data  = 1.0;
            bmin_data[pix] = axis_data[pix]->FindBin(xmin_data);
            bmax_data[pix] = axis_data[pix]->FindBin(xmax_data);
            integral_data[pix] = fHistCh_read_p[pix]->Integral(bmin_data[pix],bmax_data[pix]);
            //fHistCh_read_p[pix]->Scale(1/integral_data[pix]);
            fHistCh_read_p[pix]->Scale(1/pdf_nph_p);
            axis_data_pi[pix] = fHistCh_read_pi[pix]->GetXaxis();
            if (prtangle_pdf==90)xmin_data = 0.9;
            if (prtangle_pdf==90)xmax_data = 1.0;
            bmin_data_pi[pix] = axis_data_pi[pix]->FindBin(xmin_data);
            bmax_data_pi[pix] = axis_data_pi[pix]->FindBin(xmax_data);
            integral_data_pi[pix] = fHistCh_read_pi[pix]->Integral(bmin_data_pi[pix],bmax_data_pi[pix]);
            //fHistCh_read_pi[pix]->Scale(1/integral_data_pi[pix]);
            fHistCh_read_pi[pix]->Scale(1/pdf_nph_pi);
            fHistCh_read_pi[pix]->SetLineColor(kRed);
            fHistCh_graph_p[pix] =new TGraph(fHistCh_read_p[pix]);
            fHistCh_graph_pi[pix] =new TGraph(fHistCh_read_pi[pix]);
        }
    }
    TString outFile;
    if(method_type==3) {
        outFile = PrtManager::Instance()->GetOutName()+"_separation.root";
    } else if(gPDF==1) {
        outFile = PrtManager::Instance()->GetOutName()+"_cherenkovPDF.root";
    } else {
        outFile = PrtManager::Instance()->GetOutName()+"_spr.root";
    }

    //TString outFile =PrtManager::Instance()->GetOutName()+"_spr.root" ;
    
    
    
    
    
    
    Double_t minChangle(0);
    Double_t maxChangle(1);
    Double_t rad = TMath::Pi()/180.;
    Double_t criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
    prt_setRootPalette(1);
    prt_createMap();
    prt_initDigi();
    TFile file(outFile,"recreate");
    TTree tree("dirc","SPR");
    tree.Branch("mom", &mom,"mom/D");
    tree.Branch("tofPid", &tofPid,"tofPid/I");
    tree.Branch("distPid", &distPid,"distPid/I");
    tree.Branch("likePid", &likePid,"likePid/I");
    tree.Branch("spr", &spr,"spr/D");
    tree.Branch("trr", &trr,"trr/D");
    tree.Branch("nph",&nph,"nph/D");
    tree.Branch("cangle",&cangle,"cangle/D");
    tree.Branch("likelihood",&likelihood,"par3/D");
    tree.Branch("separation",&separation,"separation/D");
    tree.Branch("par5",&par5,"par5/D");
    tree.Branch("par6",&par6,"par6/D");
    tree.Branch("test1",&test1,"test1/D");
    tree.Branch("test2",&test2,"test2/D");
    tree.Branch("test3",&test3,"test3/D");
    tree.Branch("theta",&theta,"theta/D");
    tree.Branch("phi",&phi,"phi/D");
    
    tree.Branch("chAngleCut",&chAngleCut,"chAngleCut/D");
    tree.Branch("recoAngle",&recoAngle,"recoAngle/D");
    tree.Branch("timeRes",&timeRes,"timeRes/D");
    tree.Branch("end",&end,"end/I");
    tree.Branch("solution_number_approach_selection",&end,"solution_number_approach_selection/I");
    tree.Branch("solution_number",&end,"solution_number/I");
    
    
    
    
    
    test1 = PrtManager::Instance()->GetTest1();
    test2 = PrtManager::Instance()->GetTest2();
    test3 = PrtManager::Instance()->GetTest3();
    Int_t radiator = PrtManager::Instance()->GetRadiator();
    fMethod = PrtManager::Instance()->GetRunType();
    Int_t nEvents = fChain->GetEntries();
    if(end==0) end = nEvents;
    if (gPDF==1)start = 500001;
    
    cout<<"@@@@@@@@@@@@ gPDF="<< gPDF << "    start= "<< start << endl;
    
    std::cout<<"Run started for ["<<start<<","<<end <<"]"<<std::endl;
    Int_t nsHits(0),nsEvents(0),studyId(0), nHits(0), ninfit(1);
    if(start<0) {
        ninfit=abs(start);
        start=0;
    }
    for (Int_t ievent=start; ievent<end; ievent++) { //&& ievent<end
        fChain->GetEntry(ievent);
        nHits = fEvent->GetHitSize();
        nHits_dac->Fill(nHits);
        if(ievent%1000==0) std::cout<<"Event # "<< ievent << " has "<< nHits <<" hits"<<std::endl;
        //if (ievent > 1) if (nHits < 40 || nHits > 50) continue;
        if (ievent >1 && nHits < 17) continue;
        if(ievent-start==0) {
            tree.SetTitle(fEvent->PrintInfo());
            prtangle = fEvent->GetAngle(); //changed 2017
            // here
            std::cout<<"@@@@@@@@@@@@@@@@@@@"<<" prtangle= "<<prtangle<<std::endl;
            studyId = fEvent->GetGeometry();
            if(studyId==152 || studyId==153 || studyId==161 || studyId==162 || studyId==171 || studyId==172 || studyId==173 || studyId==175 || studyId==176 || studyId==177 || studyId==178) {
                radiator=2;
            }
            mom=fEvent->GetMomentum().Mag(); // changed 2017
            //mom=7.0; // here
            std::cout<<"No Problem1  "<< " mom ="<< mom <<std::endl; // changed
            Double_t beam_corr(0);
            //if(studyId==151) beam_corr = 0.0045; //125 deg //!
            // if(studyId==151) beam_corr = 0.001; // 20 deg
            // if(studyId==151) beam_corr = -0.003; // 25 deg
            //beam_corr = 0.002; // 125 deg s160
            if(fEvent->GetType()==0) {
                momInBar.RotateY(TMath::Pi()-prtangle*rad-test1);// 0  test1
                momInBar.RotateX(test2);//0 test2
                
                //momInBar.RotateY(TMath::Pi()-prtangle*rad+0.0);
                //momInBar.RotateX(0.0);
                // -+, --, +-, ++
                
            } else {
                momInBar.RotateY(TMath::Pi()-prtangle*rad);
            }
            // momInBar = fEvent->GetMomentum().Unit();
            if(fVerbose==3) {
                cz = momInBar.Unit();
                cz = TVector3(-cz.X(),cz.Y(),cz.Z());
            }
        }
        //if(nHits<5) continue;  // changed
        //std::cout<<"@@@@@@@@"<<" prtangle= "<<prtangle<<std::endl; // changed
        Double_t momentum=fEvent->GetMomentum().Mag(); // changed 2017
        //Double_t momentum=7.0; // here
        if( fEvent->GetType()==1) momentum /= 1000;
        tofPid=fEvent->GetParticle();
        if(tofPid==212) tofPid=211;
        //std::cout<<"$$$$$$$$$$$$$$$$$$$ up $$$$$$ tofPid "<<tofPid <<std::endl;
        
        Int_t pdg[]= {11,13,211,321,2212};
        Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
        Double_t angle1(0), angle2(0),sum1(0),sum2(0), sigma(0.009),range(5*sigma),noise(0.3);
        fAngleP = acos(sqrt(momentum*momentum+ mass[4]*mass[4])/momentum/1.4738)-0.00; //1.4738 = 370 = 3.35
        fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00; //-0.0014 for 160 25deg
        //std::cout<<" ####### ###### chrencove Angle  "<< " proton ="<< fAngleP <<std::endl; // changed
        //std::cout<<" ####### ###### chrencove Angle  "<< " pi  ="<< fAnglePi <<std::endl;   // changed
        gF1->SetParameter(0,1);
        gF2->SetParameter(0,1);
        
        // move the model
        
        gF1->SetParameter(1,fAngleP);
        gF2->SetParameter(1,fAnglePi);
        //std::cout<<"No Problem2  " <<std::endl;
        /*
         if(prtangle == 20)gF1->SetParameter(1,0.818334);
         if(prtangle == 30)gF1->SetParameter(1,0.816051);
         if(prtangle == 40)gF1->SetParameter(1,0.816986);
         if(prtangle == 50)gF1->SetParameter(1,0.817087);
         if(prtangle == 60)gF1->SetParameter(1,0.815137);
         if(prtangle == 70)gF1->SetParameter(1,0.813078);
         if(prtangle == 80)gF1->SetParameter(1,0.812915);
         if(prtangle == 90)gF1->SetParameter(1,0.810881);
         if(prtangle == 100)gF1->SetParameter(1,0.811012);
         if(prtangle == 110)gF1->SetParameter(1,0.813848);
         if(prtangle == 120)gF1->SetParameter(1,0.814565);
         if(prtangle == 130)gF1->SetParameter(1,0.817052);
         if(prtangle == 140)gF1->SetParameter(1,0.816796);
         if(prtangle == 150)gF1->SetParameter(1,0.815613);
         
         
         if(prtangle == 20)gF2->SetParameter(1,0.826920);
         if(prtangle == 30)gF2->SetParameter(1,0.824533);
         if(prtangle == 40)gF2->SetParameter(1,0.825388);
         if(prtangle == 50)gF2->SetParameter(1,0.824601);
         if(prtangle == 60)gF2->SetParameter(1,0.823733);
         if(prtangle == 70)gF2->SetParameter(1,0.820671);
         if(prtangle == 80)gF2->SetParameter(1,0.820424);
         if(prtangle == 90)gF2->SetParameter(1,0.819132);
         if(prtangle == 100)gF2->SetParameter(1,0.820097);
         if(prtangle == 110)gF2->SetParameter(1,0.821585);
         if(prtangle == 120)gF2->SetParameter(1,0.822756);
         if(prtangle == 130)gF2->SetParameter(1,0.824813);
         if(prtangle == 140)gF2->SetParameter(1,0.825208);
         if(prtangle == 150)gF2->SetParameter(1,0.823934);
         */
        
        gF1->SetParameter(2,sigma);
        gF2->SetParameter(2,sigma);
        //////////////////////////////
        // select proton sim
        if(fMethod==2 && tofPid!=2212 && fEvent->GetType()==1&& recoP==1 && recoPi==0) continue; //2212, 211 will remove pion , proton respectively in simulation
        if(fMethod==2 && tofPid!=211 && fEvent->GetType()==1&& recoP==1 && recoPi==0) goto goon;
        
        //select pi sim
        if(fMethod==2 && tofPid!=211 && fEvent->GetType()==1&& recoP==0 && recoPi==1) continue; //2212, 211 will remove pion , proton respectively in simulation
        if(fMethod==2 && tofPid!=2212 && fEvent->GetType()==1&& recoP==0 && recoPi==1) goto goon;
        
        if(fMethod==3 && fEvent->GetType()==1) goto goon;
        
        if(fEvent->GetType()==0) { // Pions and protons will enter the loop the filltering will be applyied from the tof system
            Bool_t tof1_bool (false);
            Bool_t tof2_bool (false);
            Bool_t trg1_bool (false);
            Bool_t trg2_bool = false;
            Bool_t trgmzH_bool = false;
            Bool_t trgmzV_bool = false;
            Bool_t multi_bool = false;
            Bool_t multi_bool_tof1 = false;
            Bool_t multi_bool_tof2 = false;
            Bool_t multi_bool_trg1 = false;
            Bool_t multi_bool_trg2 = false;
            Bool_t multi_bool_trgmzH = false;
            Bool_t multi_bool_trgmzV = false;
            Bool_t bad_bool_hodoH = false;
            Bool_t bad_bool_hodoV = false;
            
            
            //Vertical
            Bool_t hodo_bool_1349= false;
            Bool_t hodo_bool_1350= false;
            Bool_t hodo_bool_1351= false;
            Bool_t hodo_bool_1352= false;
            
            //Horizental
            Bool_t hodo_bool_1367= false;
            Bool_t hodo_bool_1368= false;
            Bool_t hodo_bool_1369= false;
            Bool_t hodo_bool_1370= false;
            Bool_t hodo_bool_1371= false;
            
            Double_t LE_tof1=0;
            Double_t LE_tof2=0;
            Double_t LE_trg1=0;
            Double_t LE_trg2=0;
            Double_t LE_trgmzH=0;
            Double_t LE_trgmzV=0;
            
            TOT_tof1=0;
            TOT_tof2=0;
            TOT_trg1=0;
            TOT_trg2=0;
            TOT_trgmzH=0;
            TOT_trgmzV=0;
            delta_tof2tof1=0;
            
            Int_t count_tof1 =0;
            Int_t count_tof2 =0;
            Int_t count_trg1 =0;
            Int_t count_trg2 =0;
            Int_t count_trgmzH =0;
            Int_t count_trgmzV =0;
            
            Int_t count_hodoH_1368 =0;
            Int_t count_hodoH_1369 =0;
            Int_t count_hodoH_1370 =0;
            
            
            Int_t count_hodoV_1350 =0;
            Int_t count_hodoV_1351 =0;
            Int_t count_hodoV_1352 =0;
            
            
            
            
            for(Int_t h=0; h<fEvent->GetHitSize(); h++) {
                fHit = fEvent->GetHit(h);
                Int_t gch=fHit.GetChannel();
                //std::cout<<"starting hits loop "<<std::endl;
                //std::cout<<"ievent"<<ievent<< " 	ihits= "<< h <<std::endl;
                // TOF 1
                if (gch == 1392) { //1392
                    LE_tof1 = fHit.GetLeadTime();
                    TOT_tof1= fHit.GetTotTime();
                    tof1_bool = true;
                    ++count_tof1;
                }
                if (count_tof1 > 1) {
                    multi_bool = true;
                    multi_bool_tof1= true;
                }
                // TOF 2
                if (gch == 1398) { //1398
                    LE_tof2 = fHit.GetLeadTime();
                    TOT_tof2= fHit.GetTotTime();
                    tof2_bool = true;
                    ++count_tof2;
                }
                if (count_tof2 > 1) {
                    multi_bool = true;
                    multi_bool_tof2= true;
                }
                // Trg1
                if (gch == 816) {
                    LE_trg1 = fHit.GetLeadTime();
                    TOT_trg1= fHit.GetTotTime();
                    trg1_bool = true;
                    ++count_trg1;
                }
                if (count_trg1 > 1) {
                    multi_bool = true;
                    multi_bool_trg1= true;
                }
                // trg2
                if (gch == 817) {
                    LE_trg2 = fHit.GetLeadTime();
                    TOT_trg2= fHit.GetTotTime();
                    trg2_bool = true;
                    ++count_trg2;
                }
                if (count_trg2 > 1) {
                    // multi_bool = true;
                    multi_bool_trg2= true;
                }
                // TrgmzH
                if (gch == 818) {
                    LE_trgmzH = fHit.GetLeadTime();
                    TOT_trgmzH= fHit.GetTotTime();
                    trgmzH_bool = true;
                    ++count_trgmzH;
                }
                if (count_trgmzH > 1) {
                    //  multi_bool = true;
                    multi_bool_trgmzH= true;
                }
                // TrgmzV
                if (gch == 819) {
                    LE_trgmzV = fHit.GetLeadTime();
                    TOT_trgmzV= fHit.GetTotTime();
                    trgmzV_bool = true;
                    ++count_trgmzV;
                }
                if (count_trgmzV > 1) {
                    //  multi_bool = true;
                    multi_bool_trgmzV= true;
                }
                
                if(0<=(gch-1344) && (gch-1344)<16) { // hodoscope V all
                    hodoV->Fill(gch-1344);
                }
                
                if(16<=(gch-1344) && (gch-1344)<32) { // hodoscope H all
                    hodoH->Fill(gch-1344-16);
                }
                
                
                // select bad fibers
                if (gch==1344||gch==1345||gch==1346||gch==1347||gch==1348||gch==1349||gch==1353||gch==1354||gch==1355||gch==1356||gch==1357||gch==1358||gch==1359) bad_bool_hodoV=true;
                if (gch==1360||gch==1361||gch==1362||gch==1363||gch==1364||gch==1365||gch==1366||gch==1367||gch==1371||gch==1372||gch==1373||gch==1374||gch==1375) bad_bool_hodoH=true;
                // fill good fibers
                
                // fill good fibers
                if (gch==1350 && bad_bool_hodoV==false && bad_bool_hodoH== false) { // Vertical
                    hodo_bool_1350=true;
                    for(Int_t h4=0; h4<fEvent->GetHitSize(); h4++) {
                        fHit4 = fEvent->GetHit(h4);
                        Int_t gch4=fHit4.GetChannel();
                        if (gch4==1344 ) ++count_hodoV_1350;
                        if (gch4==1345 ) ++count_hodoV_1350;
                        if (gch4==1346 ) ++count_hodoV_1350;
                        if (gch4==1347 ) ++count_hodoV_1350;
                        if (gch4==1348 ) ++count_hodoV_1350;
                        if (gch4==1349 ) ++count_hodoV_1350;
                        if (gch4==1351 ) ++count_hodoV_1350;
                        if (gch4==1352 ) ++count_hodoV_1350;
                        if (gch4==1353 ) ++count_hodoV_1350;
                        if (gch4==1354 ) ++count_hodoV_1350;
                        if (gch4==1355 ) ++count_hodoV_1350;
                        if (gch4==1356 ) ++count_hodoV_1350;
                        if (gch4==1357 ) ++count_hodoV_1350;
                        if (gch4==1358 ) ++count_hodoV_1350;
                        if (gch4==1359 ) ++count_hodoV_1350;
                        if (count_hodoV_1350==0) {
                            if (0<=(gch4-1344) && (gch4-1344)<16)  hodoV_select->Fill(gch4-1344);
                        }
                    }
                }
                if (gch==1351 && bad_bool_hodoV==false && bad_bool_hodoH== false) { // Vertical
                    hodo_bool_1351=true;
                    for(Int_t h2=0; h2<fEvent->GetHitSize(); h2++) {
                        fHit2 = fEvent->GetHit(h2);
                        Int_t gch2=fHit2.GetChannel();
                        if (gch2==1344 ) ++count_hodoV_1351;
                        if (gch2==1345 ) ++count_hodoV_1351;
                        if (gch2==1346 ) ++count_hodoV_1351;
                        if (gch2==1347 ) ++count_hodoV_1351;
                        if (gch2==1348 ) ++count_hodoV_1351;
                        if (gch2==1349 ) ++count_hodoV_1351;
                        if (gch2==1350 ) ++count_hodoV_1351;
                        if (gch2==1352 ) ++count_hodoV_1351;
                        if (gch2==1353 ) ++count_hodoV_1351;
                        if (gch2==1354 ) ++count_hodoV_1351;
                        if (gch2==1355 ) ++count_hodoV_1351;
                        if (gch2==1356 ) ++count_hodoV_1351;
                        if (gch2==1357 ) ++count_hodoV_1351;
                        if (gch2==1358 ) ++count_hodoV_1351;
                        if (gch2==1359 ) ++count_hodoV_1351;
                        if (count_hodoV_1351==0) {
                            if (0<=(gch2-1344) && (gch2-1344)<16)  hodoV_select->Fill(gch2-1344); // fill vertical fiber
                        }
                    }
                    
                }
                if (gch==1352 && bad_bool_hodoV==false && bad_bool_hodoH== false) {  // Vertical
                    hodo_bool_1352=true;
                    for(Int_t h5=0; h5<fEvent->GetHitSize(); h5++) {
                        fHit5 = fEvent->GetHit(h5);
                        Int_t gch5=fHit5.GetChannel();
                        if (gch5==1344 ) ++count_hodoV_1352;
                        if (gch5==1345 ) ++count_hodoV_1352;
                        if (gch5==1346 ) ++count_hodoV_1352;
                        if (gch5==1347 ) ++count_hodoV_1352;
                        if (gch5==1348 ) ++count_hodoV_1352;
                        if (gch5==1349 ) ++count_hodoV_1352;
                        if (gch5==1350 ) ++count_hodoV_1352;
                        if (gch5==1351 ) ++count_hodoV_1352;
                        if (gch5==1353 ) ++count_hodoV_1352;
                        if (gch5==1354 ) ++count_hodoV_1352;
                        if (gch5==1355 ) ++count_hodoV_1352;
                        if (gch5==1356 ) ++count_hodoV_1352;
                        if (gch5==1357 ) ++count_hodoV_1352;
                        if (gch5==1358 ) ++count_hodoV_1352;
                        if (gch5==1359 ) ++count_hodoV_1352;
                        if (count_hodoV_1352==0) {
                            if (0<=(gch5-1344) && (gch5-1344)<16)  hodoV_select->Fill(gch5-1344); // fill vertical fiber
                        }
                    }
                }
                /////////////////////////////////////////////////////
                if (gch==1368 && bad_bool_hodoV==false && bad_bool_hodoH== false) { // Horizental
                    hodo_bool_1368=true;
                    for(Int_t h6=0; h6<fEvent->GetHitSize(); h6++) {
                        fHit6 = fEvent->GetHit(h6);
                        Int_t gch6=fHit6.GetChannel();
                        if (gch6==1360 ) ++count_hodoH_1368;
                        if (gch6==1361 ) ++count_hodoH_1368;
                        if (gch6==1362 ) ++count_hodoH_1368;
                        if (gch6==1363 ) ++count_hodoH_1368;
                        if (gch6==1364 ) ++count_hodoH_1368;
                        if (gch6==1365 ) ++count_hodoH_1368;
                        if (gch6==1366 ) ++count_hodoH_1368;
                        if (gch6==1367 ) ++count_hodoH_1368;
                        if (gch6==1369 ) ++count_hodoH_1368;
                        if (gch6==1370 ) ++count_hodoH_1368;
                        if (gch6==1371 ) ++count_hodoH_1368;
                        if (gch6==1372 ) ++count_hodoH_1368;
                        if (gch6==1373 ) ++count_hodoH_1368;
                        if (gch6==1374 ) ++count_hodoH_1368;
                        if (count_hodoH_1368==0) {
                            if (16<=(gch6-1344) && (gch6-1344)<32) hodoH_select->Fill(gch6-1344-16);
                        }
                    }
                }
                if (gch==1369 && bad_bool_hodoV==false && bad_bool_hodoH== false) { // Horizental
                    hodo_bool_1369=true;
                    for(Int_t h3=0; h3<fEvent->GetHitSize(); h3++) {
                        fHit3 = fEvent->GetHit(h3);
                        Int_t gch3=fHit3.GetChannel();
                        if (gch3==1360 ) ++count_hodoH_1369;
                        if (gch3==1361 ) ++count_hodoH_1369;
                        if (gch3==1362 ) ++count_hodoH_1369;
                        if (gch3==1363 ) ++count_hodoH_1369;
                        if (gch3==1364 ) ++count_hodoH_1369;
                        if (gch3==1365 ) ++count_hodoH_1369;
                        if (gch3==1366 ) ++count_hodoH_1369;
                        if (gch3==1367 ) ++count_hodoH_1369;
                        if (gch3==1368 ) ++count_hodoH_1369;
                        if (gch3==1370 ) ++count_hodoH_1369;
                        if (gch3==1371 ) ++count_hodoH_1369;
                        if (gch3==1372 ) ++count_hodoH_1369;
                        if (gch3==1373 ) ++count_hodoH_1369;
                        if (gch3==1374 ) ++count_hodoH_1369;
                        if (count_hodoH_1369==0) {
                            if (16<=(gch3-1344) && (gch3-1344)<32) hodoH_select->Fill(gch3-1344-16);
                        }
                    }
                    
                }
                if (gch==1370 && bad_bool_hodoV==false && bad_bool_hodoH== false) { // Horizental
                    hodo_bool_1370=true;
                    for(Int_t h7=0; h7<fEvent->GetHitSize(); h7++) {
                        fHit7 = fEvent->GetHit(h7);
                        Int_t gch7=fHit7.GetChannel();
                        if (gch7==1360 ) ++count_hodoH_1370;
                        if (gch7==1361 ) ++count_hodoH_1370;
                        if (gch7==1362 ) ++count_hodoH_1370;
                        if (gch7==1363 ) ++count_hodoH_1370;
                        if (gch7==1364 ) ++count_hodoH_1370;
                        if (gch7==1365 ) ++count_hodoH_1370;
                        if (gch7==1366 ) ++count_hodoH_1370;
                        if (gch7==1367 ) ++count_hodoH_1370;
                        if (gch7==1368 ) ++count_hodoH_1370;
                        if (gch7==1369 ) ++count_hodoH_1370;
                        if (gch7==1371 ) ++count_hodoH_1370;
                        if (gch7==1372 ) ++count_hodoH_1370;
                        if (gch7==1373 ) ++count_hodoH_1370;
                        if (gch7==1374 ) ++count_hodoH_1370;
                        if (count_hodoH_1370==0) {
                            if (16<=(gch7-1344) && (gch7-1344)<32) hodoH_select->Fill(gch7-1344-16);
                        }
                    }
                }
                
            }// end of hit loop
            countmulti_tof1->Fill(count_tof1);
            countmulti_tof2->Fill(count_tof2);
            countmulti_trg1->Fill(count_trg1);
            countmulti_trg2->Fill(count_trg2);
            countmulti_trgmzH->Fill(count_trgmzH);
            countmulti_trgmzV->Fill(count_trgmzV);
            ///////
            countmulti_hodoV->Fill(count_hodoV_1351);
            countmulti_hodoH->Fill(count_hodoH_1369);
            // hodo all
            {
                for(int h=0; h<16; h++)
                    for(int v=0; v<16; v++)
                        hodoF->Fill(v,(15-h),hodoH->GetBinContent(h)+hodoV->GetBinContent(v));
            } // end Fill
            
            // hodo Vertical
            {
                for(int h=0; h<16; h++)
                    for(int v=0; v<16; v++)
                        hodo_multi_withmedVfiber->Fill(v,(15-h),hodoV_select->GetBinContent(v));
            }
            // hodo Horizental
            {
                for(int h=0; h<16; h++)
                    for(int v=0; v<16; v++)
                        hodo_multi_withmedHfiber->Fill(v,(15-h),hodoH_select->GetBinContent(h));
            }
            
            
            //////////////////////////////////////////////////////////////////////////////////////////
            delta_tof2tof1=LE_tof2-LE_tof1;
            
            //htof1_le->Fill(LE_tof1);
            //htof1_tot->Fill(TOT_tof1);
            htof1_le_tot->Fill(LE_tof1,TOT_tof1);
            //htof2_le->Fill(LE_tof2);
            //htof2_tot->Fill(TOT_tof2);
            htof2_le_tot->Fill(LE_tof2,TOT_tof2);
            //htrg1_le->Fill(LE_trg1);
            //htrg1_tot->Fill(TOT_trg1);
            htrg1_le_tot->Fill(LE_trg1, TOT_trg1);
            //htrg2_le->Fill(LE_trg2);
            //htrg2_tot->Fill(TOT_trg2);
            htrg2_le_tot->Fill(LE_trg2, TOT_trg2);
            //htrgmzV_le->Fill(LE_trgmzV);
            //htrgmzV_tot->Fill(TOT_trgmzV);
            htrgmzV_le_tot->Fill(LE_trgmzV, TOT_trgmzV);
            //htrgmzH_le->Fill(LE_trgmzH);
            //htrgmzH_tot->Fill(TOT_trgmzH);
            htrgmzH_le_tot->Fill(LE_trgmzH, TOT_trgmzH);
            
            hdelta_tof2tof1->Fill(delta_tof2tof1);
            
            if(/*hodo_bool_1351==true  &&  hodo_bool_1369==true &&
                count_hodoV_1351==0 && count_hodoH_1369==0 &&*/  // tight cut
               
               
               bad_bool_hodoV==false && bad_bool_hodoH== false &&
               (hodo_bool_1350==true || hodo_bool_1351==true|| hodo_bool_1352==true) && (hodo_bool_1368==true || hodo_bool_1369==true || hodo_bool_1370==true ) &&
               count_hodoV_1350 ==0 && count_hodoV_1351==0 && count_hodoV_1352==0 && count_hodoH_1368==0 && count_hodoH_1369==0 && count_hodoH_1370==0 &&
               
               
               
               
               tof1_bool == true && tof2_bool == true &&
               trg1_bool==true && trg2_bool==true && trgmzH_bool==true && trgmzV_bool==true &&
               
               LE_trg1>-144 && LE_trg1 < -132 &&
               TOT_trg1>118.4 && TOT_trg1<119.8 &&
               
               LE_tof1>-282 && LE_tof1<-270 &&
               TOT_tof1>45 && TOT_tof1<50 &&
               
               LE_tof2>-250 && LE_tof2 < -238 &&
               TOT_tof2>44 && TOT_tof2<50 &&
               
               LE_trg2> -125 && LE_trg2<-110 &&
               LE_trg1> -143 && LE_trg1<-130 &&
               
               multi_bool ==false ) { // TOF1, TOF2, TRG1
                
                // All good candidates
                //hdelta_tof2tof1->Fill(delta_tof2tof1);
                //goto goon;
                
                if (fMethod==2 && delta_tof2tof1 > 32.4 && delta_tof2tof1< 32.9 && recoP==1 && recoPi==0) goto goon; //proton Candidate
                if (fMethod==2 && delta_tof2tof1 > 31.5 && delta_tof2tof1< 32.0 && recoP==0 && recoPi==1 ) goto goon; //pi Candidate
                //std::cout<<"No Problem  methos =  " <<fMethod<<std::endl;
                if (fMethod==3 && fEvent->GetType()==0 && ((delta_tof2tof1 > 32.4 && delta_tof2tof1< 32.9) || (delta_tof2tof1 > 31.5 && delta_tof2tof1< 32.0)) ) goto goon; //with tof cuts
                //if (fMethod==3 && fEvent->GetType()==0 ) goto goon; //with tof cuts
            }
            continue; // outside the event loop
        }// if data condition
        
    goon:
        //std::cout<<"No Problem  goon  " <<std::endl;
        hdelta_tof2tof1_isproton->Fill(delta_tof2tof1);
        nHits_dac_syscut_p->Fill(nHits);
        // hodo after cuts
        {
            for(int h=0; h<16; h++)
                for(int v=0; v<16; v++)
                    hodo_afterCut->Fill(v,(15-h),hodoH_select->GetBinContent(h)+hodoV_select->GetBinContent(v));
        } // end FillF()
        
        //   //clusters search
        //   for(Int_t h=0; h<nHits; h++) {
        //     Int_t mid=fEvent->GetHit(h).GetMcpId();
        //     Int_t pid=fEvent->GetHit(h).GetPixelId()-1;
        //     mcpdata[mid][pid]=1;
        //   }
        //   getclusters();
        
        //   Int_t bad(0);
        //   for(Int_t j=0; j<prt_nmcp; j++){
        //     for(Int_t i=0; i<65; i++){
        //   	if(cluster[j][i]>7){
        // 	  bad=cluster[j][i];
        // 	  goto lllll;
        // 	}
        //     }
        //   }
        
        // lllll:
        //   std::cout<<"bad  "<<bad <<std::endl;
        
        //   if(bad){
        //     for(Int_t j=0; j<prt_nmcp; j++){
        //   	for(Int_t i=0; i<65; i++){
        //   	  mcpdata[j][i]=0;
        //   	  cluster[j][i]=0;
        //   	}
        //     }
        
        //     continue;
        //   }
        
        Int_t nGoodPhotons=0;
        Int_t nGoodPhotons_TimeCutOnly=0;
        Int_t nCandidat=0;
        Int_t nGoodPhotons_true_sim=0;
        //std::cout<<"######### Hit loop will start"<<std::endl;
        for(Int_t h=0; h<nHits; h++) {
            fHit = fEvent->GetHit(h);
            hitTime = fHit.GetLeadTime();
            if(fEvent->GetType()!=0) hitTime+=fRand.Gaus(0,0.25); // time resol. in case it was not simulated 200 ps 250
            photonEnergy= fHit.GetMomentum().Mag()*1E6;
            if(fEvent->GetType()==0) {// shift hit time for data only
                if(prtangle==20) hitTime=hitTime-0.189155;
                if(prtangle==30) hitTime=hitTime+0.030525;
                if(prtangle==40) hitTime=hitTime+0.196482;
                if(prtangle==50) hitTime=hitTime+0.149906;
                if(prtangle==60) hitTime=hitTime+0.081199;
                if(prtangle==70) hitTime=hitTime-0.102829;
                if(prtangle==80) hitTime=hitTime-0.228011;
                if(prtangle==90) hitTime=hitTime+0.038770;
                if(prtangle==100) hitTime=hitTime-0.182761;
                if(prtangle==110) hitTime=hitTime-0.189065;
                if(prtangle==120) hitTime=hitTime+0.050843;
                if(prtangle==130) hitTime=hitTime+0.222726;
                if(prtangle==140) hitTime=hitTime+0.235055;
                if(prtangle==150) hitTime=hitTime+0.019390;
            }
            //======================================== dynamic cuts for sim and data
            if(fEvent->GetType()==1 || fEvent->GetType()==0) {
                
                Double_t cut1(0);
                {   //time cuts
                    if(prtangle==20) {
                        if(hitTime<9.0 || hitTime>25 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==25 ) {
                        if(hitTime<9.5 || hitTime>30 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==30 ) {
                        if(hitTime<9.5 || hitTime>35 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==35 ) {
                        if(hitTime<9.5 || hitTime>35 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==40 ) {
                        if(hitTime<10 || hitTime>40 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==45 ) {
                        if(hitTime<10 || hitTime>43 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==50 ) {
                        if(hitTime<10 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==55 ) {
                        if(hitTime<10 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==60 ) {
                        if(hitTime<10.5 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==65 ) {
                        if(hitTime<11 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==70 ) {
                        if(hitTime<11 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==75 ) {
                        if(hitTime<11 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle==80 ) {
                        if(hitTime<11 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if (prtangle==85.0) {
                        if(hitTime< 3.0 || hitTime>50.0 ) continue;
                        //if(hitTime< 13.5 && hitTime>10.0 ) continue;
                        if(hitTime<=12)   reflected = kFALSE;
                        if(hitTime>12) reflected = kTRUE;
                    }
                    else if (prtangle==90) {
                        if(hitTime< 3.0 || hitTime>50.0 ) continue;
                        //if(hitTime< 14.0 && hitTime>10.0 ) continue;
                        if(hitTime<=12)   reflected = kFALSE;
                        if(hitTime>12) reflected = kTRUE;
                    }
                    else if(prtangle>94) {
                        if(hitTime<3 || hitTime>25) continue;
                        reflected = kFALSE;
                    }
                }
            }
            //======================================== End of dynamic cuts
            // Double_t radiatorL = (radiator==2)? 1224.9 : 1250; //plate : bar
            Double_t radiatorL = (radiator==2)? 1224.9 : 1200; //plate : bar // changed
            if( fEvent->GetType()==1) {
                lenz = radiatorL/2.-fHit.GetPosition().Z();
            } else {
                lenz = fHit.GetPosition().Z();
                // old calculaions to calculate photon path length on the bar
                if (false) {
                    //Double_t z =  fEvent->GetBeamZ()+25;
                    Double_t z =  fEvent->GetBeamZ(); // changed
                    Double_t yrot= 89.3+17.9/2;
                    Double_t xrot= 146.3;
                    //Double_t z =  600.5; // changed
                    //lenz = z-1/tan(prtangle*rad)*(122+(z-96)/tan((135-0.5*prtangle)*rad)); // changed comminted
                    // Double_t b = 122*tan(0.5*((prtangle-90)*rad));
                    // Double_t lenz = (z-96+b)/cos((prtangle-90)*rad)+b+96;
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //lenz = z-1/tan(prtangle*rad)*(122+(z-96)/tan((135-0.5*prtangle)*rad)); // 2015 configuration
                    //lenz = z-1/tan(prtangle*rad)*(81.5+(z-97)/tan((135-0.5*prtangle)*rad));       //2016
                    lenz = z-1/tan(prtangle*rad)*(yrot+(z-xrot)/tan((135-0.5*prtangle)*rad));       //2017
                    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // 81.5 vertical distance between the top of the bar to the rotation angle
                    // 97 distance between the end of the bar (prism side) and the norm between rotation point and the top of the bar
                    // 135 constant                                                                                                                                                                                                   $
                    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    // new calculations give the asme resuilts
                    Double_t alpha= 90-prtangle;
                    Double_t x = yrot* tan(alpha/2*rad);
                    Double_t zo = (z-xrot-(1+cos(alpha*rad))*x)/cos(alpha*rad); //2017
                    Double_t zfinal = zo + xrot;
                    //std::cout<<"@@@@@@@@@@@"<<" lenZ Lee= "<<lenz<<std::endl;
                    //std::cout<<"@@@@@@@@@@@"<<" lenZ Ahmed= "<<zfinal<<std::endl;
                    // path correction
                    Double_t d_correction = (zo+x)*cos(prtangle*rad); // original should be
                    //std::cout<<"@@@@@@@@@@@"<<" distance correction = "<<d_correction<<std::endl;
                    // Constants
                    Double_t beam_momentum= 7000; // MeV
                    Double_t nano = 0.000000001; // 10e-9
                    Double_t m_proton = 938.28;
                    Double_t m_pi = 139.570;
                    Double_t c = 299792458; // speed of light
                    Double_t v_proton = c* sqrt(1-(m_proton*m_proton) /(beam_momentum * beam_momentum+ m_proton*m_proton)); // 2016
                    // different way to calculate the velocity
                    //Double_t mom_check(7000), mass_check (938.272046);
                    //Double_t v_proton = ((mom_check/sqrt(mass_check*mass_check+mom_check*mom_check)*299792458));
                    Double_t v_pi = c* sqrt(1-( m_pi*m_pi) /(beam_momentum * beam_momentum+ m_pi*m_pi)); // 2016
                    Double_t  t_correction_proton =  d_correction / (v_proton*1000*nano); // 1.55228 is the calculated d_correction at polar angle 25 (calibration satge)
                    Double_t  t_correction_aligment =  (32 +yrot) / (v_proton*1000*nano); // 32 mm distance between interaction point at bar and prependicular line to the reotation point
                    std::cout << std::fixed;
                    //std::cout << std::setprecision(10);
                    if(counter==0)std::cout<<"@@@@@@@@@@@"<<" prtangle= "<<prtangle<<" d_correction_proton= "<<d_correction<< "	lenZ lee= "<< lenz<<"  lenZ Ahmed= "<<zfinal << " chAngleCut =" <<chAngleCut<<std::endl;
                    // third method for calculating lenz and d_correction
                    /*Double_t lenz_3rd_method = xrot + (447-yrot*sin((90-prtangle)*rad))/sin(prtangle*rad);
                     Double_t d_correction_3rd= lenz_3rd_method * cos(prtangle*rad);
                     if(counter==0)std::cout<<"###########"<<" prtangle= "<<prtangle<<" d_correction_proton_3rd= "<<d_correction_3rd<< "     lenz_3rd_method= "<< lenz_3rd_method<<std::endl;
                     */
                    // method 3
                    //Double_t zo_test= (447-xrot-yrot*sin((90-prtangle)*rad))/sin(prtangle*rad);
                    //Double_t lenz_test = xrot + zo_test;
                    //Double_t d_correction_test= sqrt(zo_test *zo_test- (447-xrot-yrot*sin((90-prtangle)*rad))*(447-xrot-*sin((90-prtangle)*rad)) )+yrot-yrot*cos((90-prtangle)*rad);
                    //Double_t t_correction_proton_test =  d_correction_test / (v_proton*1000*nano); // - at theta < 90  + at theta > 90 or visversa
                    //if(counter==0)std::cout<<"#############"<<" prtangle= "<<prtangle<<" d_correction_proton_test= "<<d_correction_test<< "     lenZ tes= "<< lenz_test<<std::endl;
                    //if(counter==0)std::cout<<"@@@@@@@@@@@$$$$$$$$"<<" v_proton_check= "<<v_proton_check<< "     v_proton= "<< v_proton<<std::endl;
                    ++counter;
                    /*  Double_t t_test1= 2277.6 / (v_proton*1000*nano);
                     Double_t t_test2= fHit.GetLeadTime()+ 2;
                     //std::cout<<"@@@@@@@@@@@"<<" proton from TOF1 to DIRC = "<<t_test1<<"	calib=  "<<t_test2<<std::endl;
                     if (t_test2>0){
                     count= count+1;
                     allcalib = allcalib + t_test2 ;
                     avcalib =  allcalib/count ;
                     std::cout<<"@@@@@@@@@@@"<<" avcalib= "<<avcalib<<"	count= "<<count<<std::endl;
                     }
                     */
                    //std::cout<<"@@@@@@@@@@@"<<" lenz= "<<lenz<<std::endl;
                    //hitTime = fHit.GetLeadTime() +195 - t_correction_proton ;
                    Double_t shift = -2;
                    /*if(prtangle == 20) shift = shift -1.37006;
                     if(prtangle == 25) shift = shift -1.06119;
                     if(prtangle == 30) shift = shift -0.773135;
                     //if(prtangle == 35) shift = shift ;
                     if(prtangle ==40 ) shift = shift -0.279369;
                     // if(prtangle ==45 ) shift = shift ;
                     if(prtangle ==50 ) shift = shift -0.258814;
                     //if(prtangle ==55 ) shift = shift ;
                     if(prtangle ==60 ) shift = shift -0.436126;
                     //if(prtangle ==65 ) shift = shift ;
                     if(prtangle ==70 ) shift = shift -0.570348;
                     //if(prtangle ==75 ) shift = shift ;
                     if(prtangle ==80 ) shift = shift -0.580393;
                     //if(prtangle ==85 ) shift = shift ;
                     if(prtangle ==90 ) shift = shift -0.281902;
                     //if(prtangle ==95 ) shift = shift ;
                     //if(prtangle ==100 ) shift = shift ;
                     //if(prtangle ==105 ) shift = shift ;
                     if(prtangle ==110 ) shift = shift -0.28674;
                     //if(prtangle ==115 ) shift = shift ;
                     if(prtangle == 120) shift = shift -0.173538;
                     //if(prtangle == 125) shift = shift ;
                     if(prtangle == 130) shift = shift +0.261796;
                     //if(prtangle == 135) shift = shift ;
                     if(prtangle == 140) shift = shift +0.227919;
                     //if(prtangle == 145) shift = shift ;
                     if(prtangle == 150) shift = shift +0.00865741;
                     */
                    hitTime = fHit.GetLeadTime()+ shift - t_correction_proton ;
                    //data				//======================================== dynamic cuts
                    if(fEvent->GetType()==0) {
                        Double_t cut1(0); // test 10.4
                        {   //time cuts
                            if(prtangle==20) {
                                if(hitTime<9.0 || hitTime>25 ) continue; // warning1
                                reflected = kTRUE;
                            }
                            
                            else if(prtangle==25 ) {
                                if(hitTime<9.5 || hitTime>30 ) continue;
                                reflected = kTRUE;
                            }
                            
                            else if(prtangle==30 ) {
                                if(hitTime<9.5 || hitTime>35 ) continue;
                                reflected = kTRUE;
                            }
                            else if(prtangle==35 ) {
                                if(hitTime<9.5 || hitTime>35 ) continue;
                                reflected = kTRUE;
                            }
                            else if(prtangle==40 ) {
                                if(hitTime<10 || hitTime>40 ) continue;
                                reflected = kTRUE;
                            }
                            else if(prtangle==45 ) {
                                if(hitTime<10 || hitTime>43 ) continue;
                                reflected = kTRUE;
                            }
                            else if(prtangle==50 ) {
                                if(hitTime<10 || hitTime>50 ) continue;
                                reflected = kTRUE;
                            }
                            else if(prtangle==55 ) {
                                if(hitTime<10 || hitTime>50 ) continue;
                                reflected = kTRUE;
                            }
                            
                            else if(prtangle==60 ) {
                                if(hitTime<10.5 || hitTime>50 ) continue;
                                reflected = kTRUE;
                            }
                            else if(prtangle==65 ) {
                                if(hitTime<11.0 || hitTime>50.0 ) continue;
                                reflected = kTRUE;
                            }
                            else if(prtangle==70 ) {
                                if(hitTime<11 || hitTime>50 ) continue;
                                reflected = kTRUE;
                            }
                            else if(prtangle==75 ) {
                                if(hitTime<11 || hitTime>50 ) continue;
                                reflected = kTRUE;
                            }
                            else if(prtangle==80 ) {
                                if(hitTime<11 || hitTime>50 ) continue;
                                reflected = kTRUE;
                            }
                            else if (prtangle==85.0) {
                                if(hitTime< 3.0 || hitTime>50.0 ) continue;
                                //if(hitTime< 13.5 && hitTime>10.0 ) continue;
                                if(hitTime<=12)   reflected = kFALSE;
                                if(hitTime>12) reflected = kTRUE;
                            }
                            else if (prtangle==90) {
                                if(hitTime< 3.0 || hitTime>50.0 ) continue;
                                //if(hitTime< 14.0 && hitTime>10.0 ) continue;
                                if(hitTime<=12)   reflected = kFALSE;
                                if(hitTime>12) reflected = kTRUE;
                            }
                            
                            else if(prtangle>94) {
                                if(hitTime<3 || hitTime>25) continue;
                                reflected = kFALSE;
                            }
                        }
                    }
                    //======================================== dynamic cuts
                }// end of (if false) old calculations
            }// end of else
            
            if(fVerbose==3) {
                TVector3 cd = fHit.GetMomentum();
                fHist5->Fill(cd.Theta()*TMath::Sin(cd.Phi()),cd.Theta()*TMath::Cos(cd.Phi()));
            }
            
            Int_t pixid=fHit.GetPixelId()-1;
            Int_t mcpid=fHit.GetMcpId();
            if(reflected) lenz = 2*radiatorL - lenz;
            Int_t ch = map_mpc[mcpid][pixid];
            if(prt_isBadChannel(ch)) continue;
            
            
            //if(cluster[mcpid][pixid]>8) continue;
            // Int_t x(0),y(0), piid(pixid) , nedge(0); //new
            // for(Int_t h=0; h<nHits; h++) {
            // 	Int_t pid=fEvent->GetHit(h).GetPixelId();
            // 	Int_t mid=fEvent->GetHit(h).GetMcpId();
            // 	Double_t tdif=fabs(hitTime-fEvent->GetHit(h).GetLeadTime());
            // 	if(mid!=mcpid || pid==piid || tdif>0.3) continue;
            // 	if(pid==piid-1 && piid%8!=0) y-=1;
            // 	if(pid==piid+1 && piid%8!=7) y+=1;
            // 	if(pid==piid+8 && piid<57) x-=1;
            // 	if(pid==piid-8 && piid>8)  x+=1;
            // }
            Int_t x(0),y(0), piid(pixid+1) , nedge(0); //old
            for(Int_t h_hit_loop_2=0; h_hit_loop_2<nHits; h_hit_loop_2++) {
                Int_t pid=fEvent->GetHit(h_hit_loop_2).GetPixelId();
                Int_t mid=fEvent->GetHit(h_hit_loop_2).GetMcpId();
                if(mid!=mcpid || pid==piid) continue;
                if(pid==piid-1 && piid%8!=1) x-=1;
                if(pid==piid+1 && piid%8!=0) x+=1;
                
                if(pid==piid+8 && piid<57) y+=1;
                if(pid==piid-8 && piid>8)  y-=1;
            }
            if(x== 0 && y== 0) nedge=0;
            if(x==-1 && y== 0) nedge=1;
            if(x==-1 && y== 1) nedge=2;
            if(x== 0 && y== 1) nedge=3;
            if(x== 1 && y== 1) nedge=4;
            if(x== 1 && y== 0) nedge=5;
            if(x== 1 && y==-1) nedge=6;
            if(x== 0 && y==-1) nedge=7;
            if(x==-1 && y==-1) nedge=8;
            //std::cout<< pixid << " nedge "<<nedge <<" x " <<x << "  y  "<<y<<std::endl;
            // bad channels
            Int_t sensorId = 100*mcpid+fHit.GetPixelId();
            if(sensorId==1) continue;
            //if(sensorId==231) continue;
            //if(sensorId==230) continue;
            //if(sensorId==215) continue;
            //if(sensorId==2) continue;
            Bool_t isGoodHit(false);
            Bool_t isGoodHit_true_sim(false);
            Bool_t isGoodHit_TimeCutOnly(false);
            Bool_t isCandidat(false);
            //  if(radiator==2) isGoodHit=true;
            Int_t size =fLutNode[sensorId]->Entries();
            Double_t min_time = 1000;
            Double_t min_diff_time_tangle = 0;
            
            Int_t photon_ambiguity_counter_wo=0;
            Int_t photon_ambiguity_counter_wt=0;
            Int_t photon_ambiguity_counter_wtc=0;
            for(Int_t i=0; i<size; i++) {
                
                weight = 1;//fLutNode[sensorId]->GetWeight(i);
                dird   = fLutNode[sensorId]->GetEntryCs(i,nedge); // nedge=0
                //dird   = fLutNode[sensorId]->GetEntry(i);
                evtime = fLutNode[sensorId]->GetTime(i);
                Int_t pathid = fLutNode[sensorId]->GetPathId(i);
                Bool_t samepath(false);
                if(pathid==fHit.GetPathInPrizm()) samepath=true;
                //if(fLutNode[sensorId]->GetNRefl(i)!=1 ) continue;
                //if(pathid != 130000 && pathid != 199000) continue;
                //std::cout<<"pathid "<< pathid <<std::endl;
                for(int u=0; u<4; u++) { // u<2 ,u<4
                    // if((pathid==190000 || pathid==210000) && u == 0) continue; //one from left-right
                    // if((pathid==290000 || pathid==310000) && u == 0) continue; //two from left-right
                    // if((pathid==130000 || pathid==199000) && u == 0) continue; //from up-bottom
                    if(u == 0) dir = dird;
                    if(u == 1) dir.SetXYZ( -dird.X(), dird.Y(), dird.Z());
                    if(u == 2) dir.SetXYZ( dird.X(),-dird.Y(),  dird.Z()); //no need when no divergence in vertical plane  // commint un commint 
                    if(u == 3) dir.SetXYZ( -dird.X(),-dird.Y(), dird.Z()); //no need when no divergence in vertical plane  // commint unn commint 
                    if(reflected) dir.SetXYZ( dir.X(), dir.Y(), -dir.Z());
                    if(dir.Angle(fnX1) < criticalAngle || dir.Angle(fnY1) < criticalAngle) continue; // warning
                    luttheta = dir.Theta();
                    if(luttheta > TMath::PiOver2()) luttheta = TMath::Pi()-luttheta;
                    alphAngle = luttheta/rad ;
                    //	if (alphAngle < 40.0  ) continue; // run bigger than 40
                    //	if (alphAngle > 40.0  ) continue; // run samller than 40
                    bartime = fabs(lenz/cos(luttheta)/198.);
                    //std::cout<<"@@@@@@@@"<<" lenz= "<<lenz<<std::endl;
                    //std::cout<<"@@@@@@@@"<<" cos(luttheta)= "<<cos(luttheta)<<std::endl;
                    Double_t totaltime = bartime+evtime; // photon time from the bar to the prism extracted from the LUT
                    //std::cout<<"@@@@@@@@"<<" bartime= "<<bartime<<std::endl;
                    //std::cout<<"@@@@@@@@"<<" evtime= "<<evtime<<std::endl;
                    Double_t diff_variable= totaltime-hitTime; // Ideal
                    //std::cout<<"@@@@@@@@"<<" totaltime= "<<totaltime<<std::endl;
                    //std::cout<<"@@@@@@@@"<<" hitTime= "<<hitTime<<std::endl;
                    tangle = momInBar.Angle(dir);
                    //recoAngle=fAngleP;
                    if(fabs(tangle-recoAngle)<chAngleCut) {
                        fHist1->Fill(hitTime);
                        fHist2->Fill(totaltime);
                        fHist3->Fill(fabs(totaltime),hitTime);
                        fHist0->Fill(diff_variable);
                        if(samepath)  fHist0i->Fill(diff_variable);
                        if(!samepath)  fHist0i_bg->Fill(diff_variable);
                    }
                    //falpha->Fill(alphAngle);
                    //if(samepath)  falphai->Fill(alphAngle);
                    //fHistPhotonEnergy->Fill(photonEnergy);
                    isCandidat=true;
                    if(samepath) isGoodHit_true_sim=true;
                    if(samepath) fHist_same_path_wotc->Fill(tangle ,weight);
                    ++photon_ambiguity_counter_wo;
                    // time cut
                    if(fabs(diff_variable)>timeRes) continue;
                    isGoodHit_TimeCutOnly=true;
                    ++photon_ambiguity_counter_wt;
                    
                    if(tangle >minChangle && tangle < maxChangle && tangle < 1.85) {
                        fHist->Fill(tangle ,weight);
                        if(openChCorr== 1){
                            if( prtangle==20 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.000249856;
                                if(mcpid==2) tangle += -0.000427014;
                                if(mcpid==3) tangle += 0.00231963;
                                if(mcpid==4) tangle += -0.00131774;
                                if(mcpid==5) tangle += -0.000198306;
                                if(mcpid==6) tangle += -0.0017458;
                                if(mcpid==7) tangle += -0.00787094;
                                if(mcpid==8) tangle += -0.00189573;
                            }
                            if( prtangle==30 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += 0.000745;
                                if(mcpid==1) tangle += 0.0022453;
                                if(mcpid==2) tangle += -0.00145353;
                                if(mcpid==3) tangle += 0.00150195;
                                if(mcpid==4) tangle += -0.00599335;
                                if(mcpid==5) tangle += 0.000905412;
                                if(mcpid==10) tangle += 0.00175495;
                            }
                            if( prtangle==40 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -7.17326e-05;
                                if(mcpid==1) tangle += -0.000733894;
                                if(mcpid==2) tangle += -0.00212229;
                                if(mcpid==3) tangle += -0.00126337;
                                if(mcpid==5) tangle += -0.000418482;
                                if(mcpid==6) tangle += 0.00144153;
                                if(mcpid==7) tangle += -0.00166889;
                                if(mcpid==8) tangle += 0.00505691;
                                if(mcpid==10) tangle += 0.00264607;
                            }
                            if( prtangle==50 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.00399497;
                                if(mcpid==1) tangle += -0.00172565;
                                if(mcpid==2) tangle += 0.00298899;
                                if(mcpid==3) tangle += -0.000214435;
                                if(mcpid==4) tangle += -0.000187849;
                                if(mcpid==5) tangle += -0.00322514;
                                if(mcpid==6) tangle += 0.0013623;
                                if(mcpid==7) tangle += 0.00246261;
                                if(mcpid==8) tangle += 0.00251785;
                                if(mcpid==10) tangle += 0.00521762;
                            }
                            if( prtangle==60 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.00264267;
                                if(mcpid==1) tangle += 0.00377518;
                                if(mcpid==2) tangle += 0.00622962;
                                if(mcpid==3) tangle += -0.0018508;
                                if(mcpid==4) tangle += 0.00356558;
                                if(mcpid==5) tangle += 0.00540599;
                                if(mcpid==6) tangle += 0.00246798;
                                if(mcpid==7) tangle += 0.0029743;
                                if(mcpid==8) tangle += -0.000561685;
                                if(mcpid==9) tangle += 0.00475777;
                                if(mcpid==10) tangle += 0.000285904;
                                if(mcpid==11) tangle += 0.00406557;
                            }
                            if( prtangle==70 && fEvent->GetType()==0) {
                                if(mcpid==3) tangle += -0.000519575;
                                if(mcpid==4) tangle += 0.00412223;
                                if(mcpid==5) tangle += 0.00262853;
                                if(mcpid==6) tangle += 0.00266235;
                                if(mcpid==7) tangle += 0.00306975;
                                if(mcpid==8) tangle += 0.00674443;
                                if(mcpid==9) tangle += 0.00436985;
                                if(mcpid==10) tangle += 0.00365755;
                                if(mcpid==11) tangle += 0.00172502;
                            }
                            if( prtangle==80 && fEvent->GetType()==0) {
                                if(mcpid==6) tangle += -0.00084683;
                                if(mcpid==8) tangle += 0.00231949;
                                if(mcpid==9) tangle += 0.00341753;
                                if(mcpid==10) tangle += 0.00737129;
                                if(mcpid==11) tangle += 0.0044878;
                            }
                            if( prtangle==90 && fEvent->GetType()==0) {
                                if(mcpid==6) tangle += 0.008461;
                                if(mcpid==7) tangle += 0.00414965;
                                if(mcpid==8) tangle += -0.00446982;
                                if(mcpid==9) tangle += 0.00578173;
                                if(mcpid==10) tangle += 0.00196816;
                                if(mcpid==11) tangle += 0.00868413;
                            }
                            if( prtangle==100 && fEvent->GetType()==0) {
                                if(mcpid==6) tangle += 0.00112433;
                                if(mcpid==7) tangle += 0.00856458;
                                if(mcpid==8) tangle += 0.00391686;
                                if(mcpid==9) tangle += 0.00482033;
                                if(mcpid==10) tangle += 0.00977139;
                                if(mcpid==11) tangle += 0.00458897;
                            }
                            if( prtangle==110 && fEvent->GetType()==0) {
                                if(mcpid==6) tangle += 0.00377521;
                                if(mcpid==7) tangle += 0.0082695;
                                if(mcpid==8) tangle += 0.00691154;
                                if(mcpid==9) tangle += 0.00356301;
                                if(mcpid==10) tangle += 0.0010886;
                                if(mcpid==11) tangle += 0.00247947;
                            }
                            if( prtangle==120 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.00116267;
                                if(mcpid==1) tangle += -0.000566597;
                                if(mcpid==2) tangle += 0.0059015;
                                if(mcpid==3) tangle += -0.00272398;
                                if(mcpid==4) tangle += 0.00182406;
                                if(mcpid==5) tangle += 0.00577268;
                                if(mcpid==6) tangle += 0.00124715;
                                if(mcpid==7) tangle += 0.00343227;
                                if(mcpid==8) tangle += -0.00142522;
                                if(mcpid==9) tangle += 0.0043916;
                                if(mcpid==10) tangle += 0.00343397;
                                if(mcpid==11) tangle += 0.00232637;
                            }
                            if( prtangle==130 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.002726;
                                if(mcpid==1) tangle += -0.00405196;
                                if(mcpid==2) tangle += 0.00640533;
                                if(mcpid==3) tangle += 0.000361327;
                                if(mcpid==4) tangle += 0.000339102;
                                if(mcpid==5) tangle += -0.00330543;
                                if(mcpid==6) tangle += 0.00164837;
                                if(mcpid==7) tangle += 0.00279715;
                                if(mcpid==8) tangle += 0.00304043;
                            }
                            if( prtangle==140 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.000192233;
                                if(mcpid==1) tangle += -0.00166336;
                                if(mcpid==2) tangle += -0.00135781;
                                if(mcpid==3) tangle += -0.00260317;
                                if(mcpid==5) tangle += -0.00096379;
                                if(mcpid==6) tangle += 0.00137976;
                                if(mcpid==7) tangle += -0.00136435;
                                if(mcpid==8) tangle += 0.00580495;
                                if(mcpid==10) tangle += 0.00612748;
                            }
                            if( prtangle==150 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += 0.00188749;
                                if(mcpid==1) tangle += 0.0017865;
                                if(mcpid==2) tangle += -0.000635465;
                                if(mcpid==3) tangle += 0.00218068;
                                if(mcpid==4) tangle += -0.00536816;
                                if(mcpid==5) tangle += 0.000662086;
                                if(mcpid==6) tangle += 0.00165406;
                                if(mcpid==7) tangle += -0.000691689;
                                if(mcpid==8) tangle += 0.00734245;
                                if(mcpid==10) tangle += 0.00150191;
                            }
                        }
                        
                        // correction for sim
                        if(openChCorr== 2){
                            if( prtangle==20 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.00128565;
                                if(mcpid==2) tangle += 0.000321068;
                                if(mcpid==3) tangle += 0.000101477;
                                if(mcpid==4) tangle += -0.00143456;
                                if(mcpid==5) tangle += -0.00111927;
                                if(mcpid==6) tangle += -0.000459502;
                                if(mcpid==7) tangle += -0.00677446;
                                if(mcpid==8) tangle += -0.000473108;
                            }
                            if( prtangle==30 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.000585209;
                                if(mcpid==1) tangle += -0.000298516;
                                if(mcpid==2) tangle += -0.000228988;
                                if(mcpid==3) tangle += -0.000596109;
                                if(mcpid==4) tangle += -0.00623366;
                                if(mcpid==5) tangle += 0.000409442;
                                if(mcpid==10) tangle += 0.00243757;
                            }
                            if( prtangle==40 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.00197596;
                                if(mcpid==1) tangle += -0.00213216;
                                if(mcpid==2) tangle += -0.000695968;
                                if(mcpid==3) tangle += -0.000562163;
                                if(mcpid==5) tangle += 0.000463463;
                                if(mcpid==6) tangle += 0.00206193;
                                if(mcpid==7) tangle += -0.00129395;
                                if(mcpid==8) tangle += 0.00269016;
                                if(mcpid==10) tangle += 0.00351791;
                            }
                            if( prtangle==50 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.000926266;
                                if(mcpid==1) tangle += -0.000389023;
                                if(mcpid==2) tangle += 0.000293581;
                                if(mcpid==3) tangle += 0.00167317;
                                if(mcpid==4) tangle += -8.08177e-05;
                                if(mcpid==5) tangle += 0.0012074;
                                if(mcpid==6) tangle += 0.00158284;
                                if(mcpid==7) tangle += 0.00125742;
                                if(mcpid==8) tangle += 0.0021388;
                                if(mcpid==10) tangle += 0.00484101;
                            }
                            if( prtangle==60 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += 0.00214193;
                                if(mcpid==1) tangle += -0.000518659;
                                if(mcpid==2) tangle += 0.00209713;
                                if(mcpid==3) tangle += 0.00261453;
                                if(mcpid==4) tangle += 0.00350099;
                                if(mcpid==5) tangle += 0.00291799;
                                if(mcpid==6) tangle += 0.000803695;
                                if(mcpid==7) tangle += 0.00252611;
                                if(mcpid==8) tangle += 0.00110353;
                                if(mcpid==9) tangle += 0.00253812;
                                if(mcpid==10) tangle += 0.00231193;
                                if(mcpid==11) tangle += 0.00303059;
                            }
                            if( prtangle==70 && fEvent->GetType()==0) {
                                if(mcpid==3) tangle += 0.00166194;
                                if(mcpid==4) tangle += 0.00319498;
                                if(mcpid==5) tangle += 0.00289373;
                                if(mcpid==6) tangle += 0.00303594;
                                if(mcpid==7) tangle += 0.0028858;
                                if(mcpid==8) tangle += 0.00382261;
                                if(mcpid==9) tangle += 0.00183365;
                                if(mcpid==10) tangle += 0.000560788;
                                if(mcpid==11) tangle += 0.00114054;
                            }
                            if( prtangle==80 && fEvent->GetType()==0) {
                                if(mcpid==6) tangle += 0.00037491;
                                if(mcpid==8) tangle += 1.93707e-05;
                                if(mcpid==9) tangle += 0.00238014;
                                if(mcpid==10) tangle += 0.00398277;
                                if(mcpid==11) tangle += 0.00318657;
                            }
                            if( prtangle==90 && fEvent->GetType()==0) {
                                if(mcpid==6) tangle += -0.00269974;
                                if(mcpid==7) tangle += 0.000737201;
                                if(mcpid==8) tangle += -0.00398974;
                                if(mcpid==9) tangle += 0.0037016;
                                if(mcpid==10) tangle += 0.000418994;
                                if(mcpid==11) tangle += 0.00619974;
                            }
                            if( prtangle==100 && fEvent->GetType()==0) {
                                if(mcpid==6) tangle += -0.00101829;
                                if(mcpid==7) tangle += 0.00814329;
                                if(mcpid==8) tangle += 0.00010834;
                                if(mcpid==9) tangle += 0.00342401;
                                if(mcpid==10) tangle += 0.00106718;
                                if(mcpid==11) tangle += 0.00521473;
                            }
                            if( prtangle==110 && fEvent->GetType()==0) {
                                if(mcpid==6) tangle += 0.00138531;
                                if(mcpid==7) tangle += 0.00507243;
                                if(mcpid==8) tangle += 0.00487492;
                                if(mcpid==9) tangle += -0.000456183;
                                if(mcpid==10) tangle += 0.000136379;
                                if(mcpid==11) tangle += 0.00362543;
                            }
                            if( prtangle==120 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += 0.00375808;
                                if(mcpid==1) tangle += -0.000560249;
                                if(mcpid==2) tangle += 0.00196725;
                                if(mcpid==3) tangle += 0.00381057;
                                if(mcpid==4) tangle += 0.00182874;
                                if(mcpid==5) tangle += 0.00161296;
                                if(mcpid==6) tangle += 0.00139796;
                                if(mcpid==7) tangle += 0.00337907;
                                if(mcpid==8) tangle += -0.000640502;
                                if(mcpid==9) tangle += 0.0025486;
                                if(mcpid==10) tangle += 0.00262243;
                                if(mcpid==11) tangle += 0.00178728;
                            }
                            if( prtangle==130 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.00130451;
                                if(mcpid==1) tangle += -0.00443381;
                                if(mcpid==2) tangle += 0.00506885;
                                if(mcpid==3) tangle += 0.0018132;
                                if(mcpid==4) tangle += 0.000573663;
                                if(mcpid==5) tangle += -0.000856538;
                                if(mcpid==6) tangle += 0.00241948;
                                if(mcpid==7) tangle += 0.00228173;
                                if(mcpid==8) tangle += 0.00219046;
                            }
                            if( prtangle==140 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.00301737;
                                if(mcpid==1) tangle += -0.00254585;
                                if(mcpid==2) tangle += -0.00063854;
                                if(mcpid==3) tangle += -0.000462334;
                                if(mcpid==5) tangle += 0.000916293;
                                if(mcpid==6) tangle += 0.00258798;
                                if(mcpid==7) tangle += -0.00161678;
                                if(mcpid==8) tangle += 0.00297606;
                                if(mcpid==10) tangle += 0.00269502;
                            }
                            if( prtangle==150 && fEvent->GetType()==0) {
                                if(mcpid==0) tangle += -0.000659255;
                                if(mcpid==1) tangle += 0.000568309;
                                if(mcpid==2) tangle += 0.000781704;
                                if(mcpid==3) tangle += -0.000212239;
                                if(mcpid==4) tangle += -0.0058777;
                                if(mcpid==5) tangle += 0.000338735;
                                if(mcpid==6) tangle += 0.00409483;
                                if(mcpid==7) tangle += -0.00126457;
                                if(mcpid==8) tangle += 0.0105896;
                                if(mcpid==10) tangle += 0.00132187;
                            }
                        }
                        fHist_correction->Fill(tangle ,weight);
                        fHist_copy->Fill(tangle ,weight);
                        if(samepath) fHist_same_path->Fill(tangle ,weight);
                        if(!samepath) fHist_bg->Fill(tangle ,weight);
                        if(tofPid==2212 && ((recoP==1 && recoPi==0) || fEvent->GetType()==1)) fHistMcp[mcpid]->Fill(tangle ,weight); // proton canditate
                        if(tofPid==211 && ((recoP==0 && recoPi==1)  || fEvent->GetType()==1)) fHistMcp[mcpid]->Fill(tangle ,weight); // pi candidate changed
                        if(tofPid==2212 && fEvent->GetType()==1 && samepath) fHistMcp_same_path[mcpid]->Fill(tangle ,weight); // proton canditate
                        if(tofPid==211 &&  fEvent->GetType()==1 && samepath) fHistMcp_same_path[mcpid]->Fill(tangle ,weight); // pi candidate changed
                        
                        solution_number++;
                        if(gPDF ==1) fHistCh[ch]->Fill(tangle ,weight); // not used
                        
                        if(0.7<tangle && tangle<0.9) {
                            if(fabs(tangle-recoAngle)<chAngleCut) {
                                isGoodHit=true; //default 0.04 rad = 40 mrad
                                ++photon_ambiguity_counter_wtc;
                            }
                            if(radiator==2) { // plate
                                if(fabs(tangle-recoAngle)<chAngleCut) isGoodHit=true;
                            }
                        }
                        
                        
                        
                        
                        //if(tofPid==211 &&fabs(tangle-fAnglePi)> chAngleCut) continue;
                        //if(tofPid==2212 &&fabs(tangle-fAngleP)> chAngleCut) continue;
                /*        
                         if(true && tangle>0.65 && tangle<0.95 ) {
                         sum1 += TMath::Log(gF1->Eval(tangle)+noise);
                         sum2 += TMath::Log(gF2->Eval(tangle)+noise);
                         }
                  */       
                        //~ if(method_type =! 3 && tangle>0.65 && tangle<0.95 ) {
                        //~ sum1 += TMath::Log(gr_p->Eval(tangle));
                        //~ sum2 += TMath::Log(gr_pi->Eval(tangle));
                        //~ //std::cout<<"No Problem  separation  " <<sum1<<""<<sum2<<std::endl;
                        //~ }
                        // old separation power calculations
                        
                        if(method_type == 3 && tangle>0.65 && tangle<0.95 ) {
                        //sum1 += TMath::Log(fHistCh_graph_p[ch]->Eval(tangle)); // use graphs 
                        //sum2 += TMath::Log(fHistCh_graph_pi[ch]->Eval(tangle)); // use graphs
			
			// use histograms with normalization 
			Int_t kp = fHistCh_read_p[ch]->GetXaxis()->FindBin(tangle);
			Int_t kpi = fHistCh_read_pi[ch]->GetXaxis()->FindBin(tangle);
                        //sum1 += TMath::Log(fHistCh_read_p[ch]->GetBinContent(kp));
                        //sum2 += TMath::Log(fHistCh_read_pi[ch]->GetBinContent(kpi));
                        sum1 += TMath::Log(fHistCh_read_p[ch]->GetBinContent(kp)/pdf_nph_p);
                        sum2 += TMath::Log(fHistCh_read_pi[ch]->GetBinContent(kpi)/pdf_nph_pi);
                        //std::cout<<"No Problem  separation  " <<kp<<" "<<kp<<std::endl;
                        }

                        if(fVerbose==3) {
                            TVector3 rdir = TVector3(-dir.X(),dir.Y(),dir.Z());
                            rdir.RotateUz(cz);
                            Double_t phi = rdir.Phi();
                            Double_t tt =  rdir.Theta();
                            fHist4->Fill(tt*TMath::Sin(phi),tt*TMath::Cos(phi));
                            //for cherenckov circle fit
                            gg_gr.SetPoint(gg_i,tt*TMath::Sin(phi),tt*TMath::Cos(phi));
                            gg_i++;
                        }
                    }// end of tangle if condition
                    if(fabs(diff_variable)<min_time) {
                        min_time= diff_variable;
                        min_diff_time_tangle= tangle;
                        
                    }
                }// end of bar ambiguity loop
                //std::cout<<"No Problem  ambiguity  " <<std::endl;
            }// end of lut loop
            
            histo_photon_ambiguity_wo->Fill(photon_ambiguity_counter_wo,weight);
            histo_photon_ambiguity_wt->Fill(photon_ambiguity_counter_wt,weight);
            histo_photon_ambiguity_wtc->Fill(photon_ambiguity_counter_wtc,weight);
           /*
            //fHist_copy->Fill(min_diff_time_tangle ,weight);
            if(gPDF ==1) fHistCh[ch]->Fill(min_diff_time_tangle ,weight);
            // old separation power calculations
            if(method_type == 3 && min_diff_time_tangle>0.65 && min_diff_time_tangle<0.95 && isGoodHit_TimeCutOnly && isGoodHit) {
                sum1 += TMath::Log(fHistCh_graph_p[ch]->Eval(min_diff_time_tangle));
                sum2 += TMath::Log(fHistCh_graph_pi[ch]->Eval(min_diff_time_tangle));
                solution_number_approach_selection++;
                //sum1 += TMath::Log(gF1->Eval(tangle)+noise);
                //sum2 += TMath::Log(gF2->Eval(tangle)+noise);
         //       std::cout<<"No Problem  separation  " <<sum1<<" "<<sum2<<std::endl;
            }
           */             
            
            
            //std::cout<<"No Problem  lut loop  " <<std::endl;
            if(isGoodHit) {
                nsHits++;
                prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
                nGoodPhotons++;
            }
            if(isGoodHit_true_sim)nGoodPhotons_true_sim++;
            if(isGoodHit_TimeCutOnly) nGoodPhotons_TimeCutOnly++;
            if(isCandidat) nCandidat++;
            
        }// end of hit loop
        fnHits_p_good->Fill(nGoodPhotons);
        fnHits_p->Fill(nGoodPhotons_TimeCutOnly);
        fnHits->Fill(nCandidat);
        fnHits_true_sim->Fill(nGoodPhotons_true_sim);
        //std::cout<<"@@@@@@@@@@@"<<" nGoodPhotons= "<<nGoodPhotons<<std::endl;
        
        // for(Int_t j=0; j<prt_nmcp; j++){
        //   for(Int_t i=0; i<65; i++){
        // 	mcpdata[j][i]=0;
        // 	cluster[j][i]=0;
        //   }
        // }
        
        Double_t sum = sum1-sum2;
        if(sum!=0) {
            if(tofPid==2212) hLnDiffP->Fill(sum);
            if(tofPid==211) hLnDiffPi->Fill(sum);
            //std::cout<<"@@@@@@@@@@@@@ tofPid "<<tofPid <<std::endl;
            likelihood=sum;
        }
        //std::cout<<"No Problem  likelihood  " <<std::endl;
        // if(fVerbose==1){
        //   prt_canvasAdd("ff",800,400);
        //   gF1->Draw();
        //   gF2->SetLineColor(4);
        //   gF2->Draw("same");
        
        //   prt_waitPrimitive("ff");
        //   prt_canvasDel("ff");
        //   //prt_canvasSave(1,0);
        //   //prt_canvasDel(Form("lh_%d",gg_ind));
        // }
        if(fVerbose>0 &&  fMethod==3 && nsEvents%ninfit==0) {
            if(nsHits>10) {
                // if(tofPid==2212 && sum > 0){
                //   std::cout<<"p  "<<sum1 << "   pi "<<sum2 << "  s "<< sum<<std::endl;
                //   if(fVerbose>0)  if(!FindPeak(cangle,spr, prtangle, tofPid)) continue;
                // }
                FindPeak(cangle,spr, prtangle, tofPid);
                test1 = fTest;
                distPid = FindPdg(momentum,cangle);
                nph = nsHits/(Double_t)ninfit;
                //std::cout<<"@@@@@@@@@@@@@ nph "<<nph <<std::endl;
                spr = spr*1000;
                //std::cout<<"@@@@@@@@@@@@@ spr "<<spr <<std::endl;
                trr = spr/sqrt(nph);
                //std::cout<<"@@@@@@@@@@@@@ trr "<<trr <<std::endl;
                theta = fEvent->GetAngle(); // here
                // here
                //std::cout<<"@@@@@@@@@@@@@ theta "<<theta <<std::endl;
                par3 = fEvent->GetTest1();
                //std::cout<<"@@@@@@@@@@@@@ par3 "<<par3 <<std::endl;
                tree.Fill();
                //std::cout<<"no problem 3 1"<<std::endl;
            }
            //std::cout<<"no problem 3 2"<<std::endl;
            ResetHists();
            nsHits=0;
            //std::cout<<"no problem 3 3"<<std::endl;
        }
        //std::cout<<"no problem 3 4"<<std::endl;
        if(++nsEvents>=end) break;
        //std::cout<<"no problem 3 5"<<std::endl;
    }// end of event loop
    //std::cout<<"no problem 3 6"<<std::endl;
    if(fMethod==2) {
        FindPeak(cangle,spr, prtangle);
        nph = nsHits/(Double_t)nsEvents;
        std::cout<<"@@@@@@@@@@@@@  nsEvents ="<<nsEvents <<std::endl;
        std::cout<<"@@@@@@@@@@@@@  nsHits ="<<nsHits <<std::endl;
        std::cout<<"@@@@@@@@@@@@@  nph ="<<nph <<std::endl;
        spr = spr*1000;
        trr = spr/sqrt(nph);
        theta = fEvent->GetAngle(); // here
        // here
        //theta = 60.0; // here
        par3 = fEvent->GetTest1();
        if(fVerbose) std::cout<<Form("SPR=%2.2F N=%2.2f",spr,nph)<<std::endl;
        tree.Fill();
    } else {
        TString nid = Form("_%2.0f", prtangle);
        if(fVerbose<2) gROOT->SetBatch(1);
        prt_canvasAdd("r_lhood"+nid,800,400);
        prt_normalize(hLnDiffP,hLnDiffPi);
        hLnDiffP->SetLineColor(2);
        TF1 *ff;
        Double_t m1(1),m2(1),s1(1),s2(1);
        std::cout<<"no problem 3 7"<<std::endl;
        if(hLnDiffP->GetEntries()>10) {
            hLnDiffP->Fit("gaus","S");
            ff = hLnDiffP->GetFunction("gaus");
            m1=ff->GetParameter(1);
            s1=ff->GetParameter(2);
        }
        std::cout<<"no problem 3 8"<<std::endl;
        if(hLnDiffPi->GetEntries()>10) {
            hLnDiffPi->Fit("gaus","S");
            ff = hLnDiffPi->GetFunction("gaus");
            m2=ff->GetParameter(1);
            s2=ff->GetParameter(2);
        }
        std::cout<<"no problem 3 9"<<std::endl;
        separation = (fabs(m2-m1))/(0.5*(s1+s2));
        std::cout<<"separation "<< separation <<std::endl;
        //gStyle->SetOptFit(0);
        //gStyle->SetOptStat(0);
        hLnDiffP->SetName(Form("s_%2.2f",separation));
        hLnDiffP->Draw();
        hLnDiffPi->SetLineColor(4);
        hLnDiffPi->Draw("same");
        ///////////////////////////////////////////////////////////////////////
        prt_canvasAdd("r_delta"+nid,800,400);
        hdelta_tof2tof1->SetTitle(Form("theta %3.1f", prtangle));
        hdelta_tof2tof1->SetStats(0);
        hdelta_tof2tof1->SetLineColor(kRed);
        hdelta_tof2tof1_isproton->SetLineColor(kBlue);
        hdelta_tof2tof1->GetXaxis()->SetTitle("LE TOF2 - TOF1 [ns]");
        hdelta_tof2tof1->GetYaxis()->SetTitle("count[#]");
        hdelta_tof2tof1->Draw();
        hdelta_tof2tof1_isproton->Draw("same");
        ///////////////////////////////////////////////////////////////////////
        prt_canvasAdd("r_test_hist"+nid,800,400);
        
        
        gr_pi->SetMarkerStyle(6);
        gr_pi->SetMarkerSize(21);
        gr_pi->SetLineColor(kRed);
        //gr_pi->Sort();
        
        gr_p->SetMarkerStyle(6);
        gr_p->SetMarkerSize(21);
        gr_p->SetLineColor(kBlue);
        //gr_p->Sort();
        
        
        mg->Add(gr_p);
        mg->Add(gr_pi);
        mg->SetTitle(" Modle used in Liklehood calculations ;#theta_{C} [rad];  [#]");
        mg->Draw("APL");
        
        //  gr_p->Draw("PL");
        //  gr_pi->Draw("SAMEPL");
        
        
        //prt_waitPrimitive("r_test_hist"+nid);// wait here
        prt_canvasSave(2,0);
        //prt_waitPrimitive("r_lhood","w");
        if(fVerbose) gROOT->SetBatch(0);
        tree.Fill();
    }
    std::cout<<"solution_number_approach_selection= "<<solution_number_approach_selection<<std::endl;
    std::cout<<"solution_number= "<<solution_number<<std::endl;
    tree.Write();
    fHist->Write();
    fHist_copy->Write();
    fHist_correction->Write();
    
    fHist_same_path->Write();
    fHist_bg->Write();
    fnHits->Write();
    fnHits_true_sim->Write();
    fnHits_p->Write();
    fnHits_p_good->Write();
    fHist0->Write();
    fHist0i->Write();
    fHist0i_bg->Write();
    fHist1->Write();
    fHist2->Write();
    nHits_dac_syscut_p->Write();
    nHits_dac->Write();
    hdelta_tof2tof1->Write();
    
    fHist_same_path_wotc->Write();
    
    histo_photon_ambiguity_wo->Write();
    histo_photon_ambiguity_wt->Write();
    histo_photon_ambiguity_wtc->Write();
    
    if(gPDF ==1) {
        for(Int_t i=0; i<960; i++) {
            fHistCh[i]->Write();
        }
    }
    for(Int_t i=0; i<prt_nmcp; i++) {
        fHistMcp[i]->Write();
    }
    for(Int_t i=0; i<prt_nmcp; i++) {
        fHistMcp_same_path[i]->Write();
    }
    file.Write();
    std::cout<<"no problem 3 13"<<std::endl;
    ///if(fVerbose) ResetHists();
}// end of the function

Int_t g_num =0;
Bool_t PrtLutReco::FindPeak(Double_t& cangle, Double_t& spr, Double_t a, Int_t tofpdg) {
    cangle=0;
    spr=0;
    //  gStyle->SetCanvasPreferGL(kTRUE);
    if(fHist->GetEntries()>20 || fHistPi->GetEntries()>20) {
        gROOT->SetBatch(1);
        Int_t nfound = fSpect->Search(fHist,1,"",0.9); //0.6
        if(nfound>0) cangle = fSpect->GetPositionX()[0];
        else cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
        cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
        if(cangle>0.85) cangle=0.82;
        fFit->SetParameters(100,cangle,0.010);
        fFit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
        fFit->SetParLimits(0,0.1,1E6);
        fFit->SetParLimits(1,cangle-0.04,cangle+0.04); // changed 0.04
        fFit->SetParLimits(2,0.005,0.018); // width 7-10 // changed 0.014
        //fFit->FixParameter(2,0.01);
        //fFit->FixParameter(3,0);
        //fFit->FixParameter(4,0);
        Int_t status(0);
        if(fMethod==3) status = fHist->Fit("fgaus","lq","",0.6,1);
        else status =fHist->Fit("fgaus","M","",cangle-0.06,cangle+0.06);
        Double_t chi = fFit->GetChisquare()/fFit->GetNDF();
        // if(fFit->GetParError(1)>0.0035){
        // //   // if(fFit->GetParameter(2)>0.011){
        // //   // if(fabs(chi-1<0.3 ){
        //   spr=0;
        //   cangle=0;
        //   fTest=0;
        //   return false;
        // }else{
        //   fTest=chi;
        // }
        cangle = fFit->GetParameter(1);
        spr = fFit->GetParameter(2);
        if(fVerbose>1) gROOT->SetBatch(0);
        if(fMethod==2 && fVerbose>0) {
            // Cherenkov correction
            if (false) {
                fFit->SetParLimits(2,0.004,0.008); // width 7-10
                for(Int_t i=0; i<prt_nmcp; i++) {
                    prt_canvasAdd(Form("r_tangle_%d",i),800,400);
                    fHistMcp[i]->Fit("fgaus","lq","",fAngleP-0.03,fAngleP+0.03);
                    std::cout<<"if(mcpid=="<< i<<") tangle += "<<fAngleP-fFit->GetParameter(1)<<";" <<std::endl;
                    fHistMcp[i]->Draw();
                    prt_waitPrimitive("r_tangle_"+i);// wait here
                }
            }
            /*
             fFit->SetParLimits(2,0.004,0.008); // width 7-10
             for(Int_t i=0; i<prt_nmcp; i++){
             prt_canvasAdd(Form("r_tangle_%d",i),800,400);
             fHistMcp[i]->Fit("fgaus","lq","",fAnglePi-0.03,fAnglePi+0.03);
             std::cout<<"if(mcpid=="<< i<<") tangle += "<<fAnglePi-fFit->GetParameter(1)<<";" <<std::endl;
             fHistMcp[i]->Draw();
             }
             */
            // for(Int_t i=0; i<960; i++){
            // 	prt_canvasAdd(Form("r_tangle_ch_%d",i),800,400);
            // 	fHistCh[i]->Fit("fgaus","lq","",fAngleP-0.03,fAngleP+0.03);
            // 	std::cout<<"if(ch=="<< i<<") tangle += "<<fAngleP-fFit->GetParameter(1)<<";" <<std::endl;
            // 	fHistCh[i]->Draw();
            // }
            //      TString name = Form("r_tangle_%3.1f",PrtManager::Instance()->GetTest3());
            TString nid = Form("_%2.0f",a);
            
            //~
            //~ for(Int_t pix=495; pix<500; pix++) {
            //~ prt_canvasAdd(Form("r_pix_%d",pix),800,400);
            //~ fHistCh_graph_p[pix]->Draw();
            //~ prt_canvasGet(Form("r_pix_%d",pix))->Update();
            //~ TLine *lin_ch_p_v = new TLine(0,0,0,1000);
            //~ lin_ch_p_v->SetX1(fAngleP);
            //~ lin_ch_p_v->SetX2(fAngleP);
            //~ lin_ch_p_v->SetY1(gPad->GetUymin());
            //~ lin_ch_p_v->SetY2(gPad->GetUymax());
            //~ lin_ch_p_v->SetLineColor(kRed);
            //~ lin_ch_p_v->Draw();
            //~
            //~ }
            
            
            /*
             ///////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_tangle"+nid,800,400);
             //fHist->SetTitle(Form("theta %3.1f , TOF PID = %d", a, tofpdg));
             fHist->SetTitle(Form("theta %3.1f , proton", a));
             fHist->SetMinimum(0);
             //fHist->Scale(1/fHist->GetMaximum());
             //prt_normalize(fHist,fHistPi);
             prt_normalize(fHist,fHist_same_path);
             
             fHistPi->SetLineColor(2);
             fHist_same_path->SetLineColor(2);
             fHist->Draw();
             fHist_same_path->Draw("same");
             fHist_bg->Draw("same");
             fHistPi->Draw("same");
             // gF1->Draw("same");
             // gF2->Draw("same");
             fHisti->SetLineColor(kRed+2);
             if(fHisti->GetEntries()>5) fHisti->Draw("same");
             prt_canvasGet("r_tangle"+nid)->Update();
             TLine *line = new TLine(0,0,0,1000);
             line->SetX1(fAngleP);
             line->SetX2(fAngleP);
             line->SetY1(gPad->GetUymin());
             line->SetY2(gPad->GetUymax());
             line->SetLineColor(kRed);
             line->Draw();
             TLine *line2 = new TLine(0,0,0,1000);
             line2->SetX1(fAnglePi);
             line2->SetX2(fAnglePi);
             line2->SetY1(gPad->GetUymin());
             line2->SetY2(gPad->GetUymax());
             line2->SetLineColor(kBlue);
             line2->Draw();
             std::cout<<"fAnglePi "<< fAnglePi<<std::endl;
             std::cout<<"fAngleP "<< fAngleP<<std::endl;
             //prt_waitPrimitive("r_tangle"+nid);// wait here
             
             //////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_fnHits"+nid,800,400);
             fnHits->SetTitle(Form("Theta %3.1f", a));
             //fnHits->SetStats(0);
             fnHits_p->SetLineColor(kGreen);
             fnHits_p_good->SetLineColor(kRed);
             fnHits->Draw();
             fnHits_p->Draw("same");
             fnHits_p_good->Draw("same");
             //std::cout<<"@@@@@@@@@@@@@  fnHits_p_good mean  ="<< fnHits_p_good->GetMean() <<std::endl;
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_delta"+nid,800,400);
             hdelta_tof2tof1->SetTitle(Form("theta %3.1f", a));
             hdelta_tof2tof1->SetStats(0);
             hdelta_tof2tof1->SetLineColor(kRed);
             hdelta_tof2tof1_isproton->SetLineColor(kBlue);
             hdelta_tof2tof1->GetXaxis()->SetTitle("LE TOF2 - TOF1 [ns]");
             hdelta_tof2tof1->GetYaxis()->SetTitle("count[#]");
             hdelta_tof2tof1->Draw();
             hdelta_tof2tof1_isproton->Draw("same");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htof1_le_tot"+nid,800,400);
             htof1_le_tot->SetTitle(Form("theta %3.1f", a));
             htof1_le_tot->SetStats(0);
             htof1_le_tot->Draw("colz");
             htof1_le_tot->SetXTitle("LE [ns]");
             htof1_le_tot->SetYTitle("TOT [ns]");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htof2_le_tot"+nid,800,400);
             htof2_le_tot->SetTitle(Form("theta %3.1f", a));
             htof2_le_tot->SetStats(0);
             htof2_le_tot->Draw("colz");
             htof2_le_tot->SetXTitle("LE [ns]");
             htof2_le_tot->SetYTitle("TOT [ns]");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htrg1_le_tot"+nid,800,400);
             htrg1_le_tot->SetTitle(Form("theta %3.1f", a));
             htrg1_le_tot->SetStats(0);
             htrg1_le_tot->Draw("colz");
             htrg1_le_tot->SetXTitle("LE [ns]");
             htrg1_le_tot->SetYTitle("TOT [ns]");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htrg2_le_tot"+nid,800,400);
             htrg2_le_tot->SetTitle(Form("theta %3.1f", a));
             htrg2_le_tot->SetStats(0);
             htrg2_le_tot->Draw("colz");
             htrg2_le_tot->SetXTitle("LE [ns]");
             htrg2_le_tot->SetYTitle("TOT [ns]");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htrgmzH_le_tot"+nid,800,400);
             htrgmzH_le_tot->SetTitle(Form("theta %3.1f", a));
             htrgmzH_le_tot->SetStats(0);
             htrgmzH_le_tot->Draw("colz");
             htrgmzH_le_tot->SetXTitle("LE [ns]");
             htrgmzH_le_tot->SetYTitle("TOT [ns]");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htrgmzV_le_tot"+nid,800,400);
             htrgmzV_le_tot->SetTitle(Form("theta %3.1f", a));
             htrgmzV_le_tot->SetStats(0);
             htrgmzV_le_tot->Draw("colz");
             htrgmzV_le_tot->SetXTitle("LE [ns]");
             htrgmzV_le_tot->SetYTitle("TOT [ns]");
             /////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_time"+nid,800,400);
             fHist1->SetTitle(Form("theta %3.1f", a));
             fHist1->SetLineColor(2);
             fHist1->Draw();
             fHist2->Draw("same");
             ////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_diff"+nid,800,400);
             fHist0->SetTitle(Form("theta %3.1f", a));
             fHist0->Draw();
             fHist0i->SetLineColor(kRed+2);
             if(fHist0i->GetEntries()>5)  fHist0i->Draw("same");
             //////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_measured_time"+nid,800,400);
             fHist1->SetTitle(Form("theta %3.1f", a));
             fHist1->SetLineColor(2);
             fHist1->Draw();
             //fHist2->Draw("same");
             /////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_cm"+nid,800,400);
             fHist3->SetTitle(Form("theta %3.1f", a));
             fHist3->Draw("colz");
             ////////////////////////////////////////////////////////////////////
             prt_drawDigi("m,p,v\n", 2017);//2
             //TPaveText* tit;
             //tit = new TPaveText(17.0,5,25.,5,"NB");
             //tit->SetFillColor(0);
             //tit->AddText(Form("theta = %f", a));
             prt_canvasAdd(prt_cdigi);
             prt_cdigi->SetName("r_hp"+nid);
             prt_cdigi->SetTitle(Form("theta %3.1f", a));
             //prt_cdigi->Draw();
             //tit->Draw();
             std::cout<<"no problem 4"<<std::endl;
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo"+nid,800,400);
             hodoF->SetTitle(Form("theta %3.1f", a));
             hodoF->SetStats(0);
             hodoF->SetTitle("hodo");
             hodoF->GetXaxis()->SetTitle("[mm]");
             hodoF->GetYaxis()->SetTitle("[mm]");
             hodoF->Draw("colz");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo_afterCut"+nid,800,400);
             hodo_afterCut->SetTitle(Form("theta %3.1f", a));
             hodo_afterCut->SetStats(0);
             hodo_afterCut->SetTitle("hodo");
             hodo_afterCut->GetXaxis()->SetTitle("[mm]");
             hodo_afterCut->GetYaxis()->SetTitle("[mm]");
             hodo_afterCut->Draw("colz");
             //////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo_V_pos"+nid,800,500);
             hodo_multi_withmedVfiber->SetTitle(Form("theta %3.1f", a));
             hodo_multi_withmedVfiber->SetStats(0);
             hodo_multi_withmedVfiber->SetTitle("hodo V");
             hodo_multi_withmedVfiber->GetXaxis()->SetTitle("[mm]");
             hodo_multi_withmedVfiber->GetYaxis()->SetTitle("[mm]");
             hodo_multi_withmedVfiber->Draw("colz");
             //////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo_H_pos"+nid,800,500);
             hodo_multi_withmedHfiber->SetTitle(Form("theta %3.1f", a));
             hodo_multi_withmedHfiber->SetStats(0);
             hodo_multi_withmedHfiber->SetTitle("hodo H");
             hodo_multi_withmedHfiber->GetXaxis()->SetTitle("[mm]");
             hodo_multi_withmedHfiber->GetYaxis()->SetTitle("[mm]");
             hodo_multi_withmedHfiber->Draw("colz");
             //////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_multi"+nid,800,400);
             countmulti_trg1->SetTitle(Form("theta %3.1f", a));
             countmulti_trg1->SetStats(0);
             countmulti_tof1->SetLineColor(kViolet);
             countmulti_tof2->SetLineColor(kGreen);
             countmulti_trg1->SetLineColor(kBlue);
             countmulti_trg2->SetLineColor(kCyan);
             countmulti_trgmzV->SetLineColor(kRed);
             countmulti_trgmzH->SetLineColor(kMagenta);
             countmulti_trg1->GetXaxis()->SetTitle("Multiplicity");
             countmulti_trg1->GetYaxis()->SetTitle("count [#]");
             countmulti_trg1->Draw();
             countmulti_tof2->Draw("same");
             countmulti_trg2->Draw("same");
             countmulti_tof1->Draw("same");
             countmulti_trgmzV->Draw("same");
             countmulti_trgmzH->Draw("same");
             TLegend *leg_multi = new TLegend(0.5,0.6,0.8,0.85);
             leg_multi->SetFillColor(0);
             leg_multi->SetFillStyle(0);
             leg_multi->SetBorderSize(0);
             leg_multi->SetFillStyle(0);
             leg_multi->AddEntry(countmulti_trg1,"Trigger 1","lp");
             leg_multi->AddEntry(countmulti_trg2,"Trigger 2","lp");
             leg_multi->AddEntry(countmulti_trgmzV,"Trigger MZ V","lp");
             leg_multi->AddEntry(countmulti_trgmzH,"Trigger MZ H","lp");
             leg_multi->AddEntry(countmulti_tof1,"TOF 1","lp");
             leg_multi->AddEntry(countmulti_tof2,"TOF 2","lp");
             leg_multi->Draw();
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo_H_number"+nid,800,500);
             countmulti_hodoH->SetTitle(Form("theta %3.1f", a));
             countmulti_hodoH->SetStats(0);
             countmulti_hodoH->SetLineColor(kBlue);
             countmulti_hodoH->GetXaxis()->SetTitle("number of spikes[#]");
             countmulti_hodoH->GetYaxis()->SetTitle("count [#]");
             countmulti_hodoH->Draw();
             TLegend *leg_countmulti_hodoH = new TLegend(0.261905, 0.386076 , 0.907268,  0.995781 );
             leg_countmulti_hodoH->SetFillColor(0);
             leg_countmulti_hodoH->SetFillStyle(0);
             leg_countmulti_hodoH->SetBorderSize(0);
             leg_countmulti_hodoH->SetFillStyle(0);
             leg_countmulti_hodoH->AddEntry(countmulti_hodoH,"number of 'spikes' associated with selected H fiber","lp");
             leg_countmulti_hodoH->Draw();
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo_V_number"+nid,800,500);
             countmulti_hodoV->SetTitle(Form("theta %3.1f", a));
             countmulti_hodoV->SetStats(0);
             countmulti_hodoV->SetLineColor(kViolet);
             countmulti_hodoV->GetXaxis()->SetTitle("number of spikes[#]");
             countmulti_hodoV->GetYaxis()->SetTitle("count [#]");
             countmulti_hodoV->Draw();
             TLegend *leg_countmulti_hodoV = new TLegend(0.261905, 0.386076 , 0.907268,  0.995781 );
             leg_countmulti_hodoV->SetFillColor(0);
             leg_countmulti_hodoV->SetFillStyle(0);
             leg_countmulti_hodoV->SetBorderSize(0);
             leg_countmulti_hodoV->SetFillStyle(0);
             leg_countmulti_hodoV->AddEntry(countmulti_hodoV,"number of 'spikes' associated with selected V fiber","lp");
             leg_countmulti_hodoV->Draw();
             */
            ///////////////////////////////////////////////////////////////////////
            if(false) {
                Int_t tmax, max=0;
                for(Int_t m=0; m<prt_nmcp; m++) {
                    prt_hdigi[m]->Rebin2D(8,8);
                    prt_hdigi[m]->GetXaxis()->SetNdivisions(0);
                    prt_hdigi[m]->GetYaxis()->SetNdivisions(0);
                    prt_hdigi[m]->GetXaxis()->SetTickLength(0);
                    prt_hdigi[m]->GetYaxis()->SetTickLength(0);
                    prt_hdigi[m]->GetXaxis()->SetAxisColor(1);
                    prt_hdigi[m]->GetYaxis()->SetAxisColor(1);
                    prt_hdigi[m]->SetMarkerSize(10);
                    tmax = prt_hdigi[m]->GetMaximum();
                    if(max<tmax) max = tmax;
                }
                for(Int_t m=0; m<prt_nmcp; m++) {
                    prt_hdigi[m]->Scale(1/(Double_t)max);
                }
            }
            //////////////////////////////////////////////////////////////////////
            //std::cout<<"no problem 4.1"<<std::endl;
            /*
             TCanvas * can1 = new TCanvas("can1","can1",0,0,800,400);
             fHist1-> Draw();
             fHist2-> Draw("Same");
             can1->Print("time.png");
             TCanvas * can2 = new TCanvas("can2","can2",0,0,800,400);
             fHist->Draw();
             line->Draw();
             line2->Draw();
             can2->Print("chere.png");
             */
            /*
             TCanvas * can3 = new TCanvas("can3","can3",0,0,800,400);
             hdelta_tof2tof1->Draw();
             hdelta_tof2tof1->SetLineColor(kRed+2);
             hdelta_tof2tof1_isproton->Draw("same");
             can3->Print("delta.png");
             */
            /*
             TCanvas * can3 = new TCanvas("can3","can3",0,0,800,400);
             hdelta_tof2tof1->Draw();
             hdelta_tof2tof1->SetLineColor(kRed+2);
             hdelta_tof2tof1_isproton->Draw("same");
             */
            /*
             prt_canvasAdd("r_alpha"+nid,800,400);
             falpha->SetTitle(Form("alpha - theta %3.1f", a));
             falpha->Draw();
             falphai->SetLineColor(kRed+2);
             if(falphai->GetEntries()>5)  falphai->Draw("same");
             prt_canvasAdd("r_photonEnergy"+nid,800,400);
             fHistPhotonEnergy->SetTitle(Form("PhotonEnergy - Theta %3.1f", a));
             fHistPhotonEnergy->Draw();
             */
            //std::cout<<"no problem 5"<<std::endl;
            
            /*
             TCanvas * can6 = new TCanvas("can6","can6",0,0,800,400);
             //hHodo->Draw();
             //hHodo->SetLineColor(kRed+2);
             //hHodo->SetMaximum(3100);
             //hHodo->SetMinimum(1000);
             hodoF->SetStats(0);
             hodoF->SetTitle("hodo");
             hodoF->GetXaxis()->SetTitle("[mm]");
             hodoF->GetYaxis()->SetTitle("[mm]");
             hodoF->Draw("colz");
             //can6->Print("Hodo.png");
             */
            
            /*
             TCanvas * can4 = new TCanvas("can4","can4",0,0,800,400);
             fHist0->Draw();
             can4->Print("diff.png");
             TCanvas * can5 = new TCanvas("can5","can5",0,0,800,400);
             fHist->Draw();
             can5->Print("tangle.png");
             */
            /*
             TCanvas * can7 = new TCanvas("LE","LE",0,0,800,400);
             htrg1_le->SetLineColor(kGreen);
             htof1_le->SetLineColor(kOrange);
             htof2_le->SetLineColor(kViolet);
             htrg2_le->SetLineColor(kRed);
             htrgmzH_le->SetLineColor(kCyan);
             htrgmzV_le->SetLineColor(kMagenta);
             htrg1_le->Draw();
             htof1_le->Draw("same");
             htof2_le->Draw("same");
             htrg2_le->Draw("same");
             htrgmzH_le->Draw("same");
             htrgmzV_le->Draw("same");
             can7->Print("le.png");
             */
            
            /*
             TCanvas * can8 = new TCanvas("TOT","TOT",0,0,800,400);
             htrg1_tot->SetLineColor(kGreen);
             htof1_tot->SetLineColor(kOrange);
             htof2_tot->SetLineColor(kViolet);
             htrg2_tot->SetLineColor(kRed);
             htrgmzH_tot->SetLineColor(kCyan);
             htrgmzV_tot->SetLineColor(kMagenta);
             htrg1_tot->Draw();
             htof1_tot->Draw("same");
             htof2_tot->Draw("same");
             htrg2_tot->Draw("same");
             htrgmzH_tot->Draw("same");
             htrgmzV_tot->Draw("same");
             //can8->Print("tot.png");
             */
            /*
             htof1_tot->Fill(TOT_tof1);
             htof2_le->Fill(LE_tof2);
             htof2_tot->Fill(TOT_tof2);
             htrg1_le->Fill(LE_trg1);
             htrg1_tot->Fill(TOT_trg1);
             htrg2_le->Fill(LE_trg2);
             htrg2_tot->Fill(TOT_trg2);
             htrgmzV_le->Fill(LE_trgmzV);
             htrgmzV_tot->Fill(TOT_trgmzV);
             htrgmzH_le->Fill(LE_trgmzH);
             htrgmzH_tot->Fill(TOT_trgmzH);
             */
            std::cout<<"no problem 5.1"<<std::endl;
            //prt_waitPrimitive("r_cm"+nid);// wait here
            prt_canvasSave(2,0);
            prt_canvasDel("*");
            
            
            if(fVerbose==3) {
                TCanvas* c2 = new TCanvas("c2","c2",0,0,800,400);
                c2->Divide(2,1);
                c2->cd(1);
                fHist4->SetStats(0);
                fHist4->SetTitle(Form("Calculated from LUT, #theta = %3.1f#circ", a));
                fHist4->Draw("colz");
                Double_t x0(0), y0(0), theta(cangle);
                FitRing(x0,y0,theta);
                TVector3 corr(x0,y0,1-TMath::Sqrt(x0*x0+y0*y0));
                std::cout<<"Tcorr "<< corr.Theta()*1000<< "  Pcorr "<< corr.Phi() <<std::endl;
                TLegend *leg = new TLegend(0.5,0.7,0.85,0.87);
                //leg->SetFillColor(0);
                //leg->SetFillColorAlpha(0,0.8);
                leg->SetFillStyle(0);
                //leg->SetFillStyle(4000);
                leg->SetBorderSize(0);
                leg->AddEntry((TObject*)0,Form("Entries %0.0f",fHist4->GetEntries()),"");
                leg->AddEntry((TObject*)0,Form("#Delta#theta_{c} %f [mrad]",corr.Theta()*1000),"");
                leg->AddEntry((TObject*)0,Form("#Delta#varphi_{c} %f [mrad]",corr.Phi()),"");
                leg->Draw();
                TArc *arc = new TArc(x0,y0,theta);
                arc->SetLineColor(kRed);
                arc->SetLineWidth(1);
                arc->SetFillStyle(0);
                arc->Draw();
                gg_i=0;
                gg_gr.Set(0);
                c2->cd(2);
                gStyle->SetOptStat(1110);
                fHist5->SetTitle(Form("True from MC, #theta = %d#circ", a));
                fHist5->Draw("colz");
                //c2->Print("example.pdf");
                //c2->Print(Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/prtdirc/build/spr/tcorr_%3.1f.png", a));
                //c2->Print(Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/prtdirc/build/spr/tcorr_%3.1f.root", a));
                //c2->Print(Form("spr/tcorr_%3.1f.png", a));
                std::cout<<"no problem 5.2"<<std::endl;
                c2->Modified();
                c2->Update();
                //c2->WaitPrimitive("");
                std::cout<<"no problem 6"<<std::endl;
            }
        }
    }
    if(fVerbose<2) gROOT->SetBatch(0);
    return (cangle>0 && cangle<1);
}

void PrtLutReco::ResetHists() {
    fHist->Reset();
    fHisti->Reset();
    fHist0->Reset();
    fHist0i->Reset();
    //falpha->Reset();
    //falphai->Reset();
    fHist1->Reset();
    fHist2->Reset();
    fHist3->Reset();
    fHist4->Reset();
    
    fHist->Reset();
    fnHits->Reset();
    fnHits_p->Reset();
    fnHits_p_good->Reset();
    
    
    for(Int_t m=0; m<prt_nmcp; m++) prt_hdigi[m]->Reset();
}

TF1 *lFit = new TF1("lgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.6,0.9);
TF1 *lFitPi = new TF1("lgausPi","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.6,0.9);
Double_t PrtLutReco::fillLnDiffPPi(Double_t cangle, Int_t tofPid, Double_t mom) {
    if(fHist->GetEntries()>20 ) {
        Int_t pdg[]= {11,13,211,321,2212};
        Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
        Double_t angle1(0), angle2(0), sigma(0.006),range(0.015);
        // //fHist->Scale(1/fHist->GetMaximum());
        // Double_t d1,d2, sum1(0),sum2(0);
        // Int_t sbin = fHist->FindBin(fAngleP-range);
        // Int_t ebin = fHist->FindBin(fAngleP+range);
        // // fHist->GetXaxis()->GetNbins()
        // for(Int_t i=sbin; i< ebin; i++){
        //   if(fHist->GetBinContent(i) < 0.01 ) continue;
        //   d1 = gF1->Eval(fHist->GetBinCenter(i))- fHist->GetBinContent(i);
        //   d2 = gF1->Eval(fHist->GetBinCenter(i))- fHist->GetBinContent(i);
        //   std::cout<<"f1 "<< gF1->Eval(fHist->GetBinCenter(i)) << "   f2 "<<gF2->Eval(fHist->GetBinCenter(i)) << "    v "<< fHist->GetBinContent(i) <<std::endl;
        //   // if(d1>0) sum1+=TMath::Log(d1);
        //   // if(d2>0) sum2+=TMath::Log(d2);
        //   sum1+=TMath::Log(fabs(d1));
        //   sum2+=TMath::Log(fabs(d2));
        // }
        // Double_t amin(sum1),amin2(sum2);
        // lFit->SetRange(fAngleP-range,fAngleP+range);
        // lFit->FixParameter(0,fFit->GetParameter(0));
        // lFit->FixParameter(1,fAngleP);
        // if(fFit->GetParameter(2)>sigma) sigma=fFit->GetParameter(2);
        // lFit->FixParameter(2,sigma);
        // lFit->FixParameter(3,fFit->GetParameter(3));
        // lFit->FixParameter(4,fFit->GetParameter(4));
        lFit->SetRange(fAngleP-range,fAnglePi+range);
        lFit->FixParameter(0,fHist->GetMaximum()-0.5);
        lFit->FixParameter(1,fAngleP);
        lFit->FixParameter(2,0.01);
        lFit->FixParameter(3,0);
        lFit->FixParameter(4,0.5);
        fHist->Fit("lgaus","lq","",fAngleP-range,fAnglePi+range);
        Double_t amin,amin2,edm,errdef;
        Int_t nvpar,nparx;
        TVirtualFitter *fitter = TVirtualFitter::Fitter(fHist);
        fitter->GetStats(amin,edm,errdef,nvpar,nparx);
        // lFitPi->SetRange(fAnglePi-range,fAnglePi+range);
        // lFitPi->SetLineColor(4);
        // lFitPi->FixParameter(0,fFit->GetParameter(0));
        // lFitPi->FixParameter(1,fAnglePi);
        // lFitPi->FixParameter(2,sigma);
        // lFitPi->FixParameter(3,fFit->GetParameter(3));
        // lFitPi->FixParameter(4,fFit->GetParameter(4));
        lFitPi->SetRange(fAngleP-range,fAnglePi+range);
        lFitPi->SetLineColor(4);
        lFitPi->FixParameter(0,fHist->GetMaximum()-0.5);
        lFitPi->FixParameter(1,fAnglePi);
        lFitPi->FixParameter(2,0.01);
        lFitPi->FixParameter(3,0);
        lFitPi->FixParameter(4,0.5);
        fHist->Fit("lgausPi","lq","",fAngleP-range,fAnglePi+range);
        fitter = TVirtualFitter::Fitter(fHist);
        fitter->GetStats(amin2,edm,errdef,nvpar,nparx);
        if(fVerbose) printf("tofPid %04d | %1.4f (%1.4f/%1.4f) likelihood is %1.2f/%1.2f \n",tofPid,cangle,fAngleP,fAnglePi, amin, amin2);
        gg_ind++;
        if(fVerbose==1) {
            prt_canvasAdd("ff",800,400);
            //prt_canvasAdd(Form("lh_%d",gg_ind),800,400);
            fHist->SetTitle(Form("%d",tofPid));
            fHist->Draw();
            lFit->SetLineColor(2);
            lFit->Draw("same");
            // gF1->Draw("same");
            // gF2->SetLineColor(4);
            // gF2->Draw("same");
            //if(fabs(amin-amin2)<5)
            //prt_waitPrimitive("ff"); // wait
            prt_canvasDel("ff");
            //prt_canvasSave(1,0);
            //prt_canvasDel(Form("lh_%d",gg_ind));
        }
        return amin-amin2;
    }
    return 1000;
}

Double_t PrtLutReco::fillLnDiffPPi2(Double_t cangle, Int_t tofPid, Double_t mom) {
    if(fHist->GetEntries()>20 ) {
        Int_t pdg[]= {11,13,211,321,2212};
        Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
        Double_t angle1(0), angle2(0), sigma(0.006),range(0.03);
        Double_t d1,d2, sum1(0),sum2(0);
        Int_t sbin = fHist->FindBin(fAnglePi-range);
        Int_t ebin = fHist->FindBin(fAngleP+range);
        for(Int_t i=sbin; i< ebin; i++) {
            if(fHist->GetBinContent(i)<1 ) continue;
            d1 = 10*fabs(fHist->GetBinContent(i) *(fAngleP  - fHist->GetBinCenter(i)));
            d2 = 10*fabs(fHist->GetBinContent(i) *(fAnglePi - fHist->GetBinCenter(i)));
            if(d1>0 && d2>0) {
                std::cout<<"d1  "<<d1 << "   d2    "<< d2 <<std::endl;
                sum1+=TMath::Log(d1);
                sum2+=TMath::Log(d2);
            }
        }
        if(fVerbose) printf("tofPid %04d | %1.4f (%1.4f/%1.4f) likelihood is %1.2f/%1.2f \n",tofPid,cangle,fAngleP,fAnglePi, sum1, sum2);
        gg_ind++;
        if(fVerbose==1) {
            prt_canvasAdd("ff",800,400);
            //prt_canvasAdd(Form("lh_%d",gg_ind),800,400);
            fHist->SetTitle(Form("%d",tofPid));
            fHist->Draw();
            lFit->SetLineColor(2);
            lFit->Draw("same");
            // gFp->Draw("same");
            // gFpi->SetLineColor(4);
            // gFpi->Draw("same");
            //if(fabs(amin-amin2)<5)
            //prt_waitPrimitive("ff"); // wait
            prt_canvasDel("ff");
            //prt_canvasSave(1,0);
            //prt_canvasDel(Form("lh_%d",gg_ind));
        }
        return sum1-sum2;
    }
    return 1000;
}

void circleFcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    Int_t np = gg_gr.GetN();
    f = 0;
    Double_t *x = gg_gr.GetX();
    Double_t *y = gg_gr.GetY();
    for (Int_t i=0; i<np; i++) {
        Double_t u = x[i] + par[0];
        Double_t v = y[i] + par[1];
        Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
        f += dr*dr;
    }
    std::cout<<"fcn  "<< f<<std::endl;
}

void circleFcn2(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    Int_t np = gg_gr.GetN();
    f = 0;
    Double_t *x = gg_gr.GetX();
    Double_t *y = gg_gr.GetY();
    for (Int_t i=0; i<np; i++) {
        Double_t u = x[i] + par[0];
        Double_t v = y[i] + par[1];
        Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
        if(dr>0.07) f += dr*dr;
        else f += fabs(dr);
    }
}

void PrtLutReco::FitRing(Double_t& x0, Double_t& y0, Double_t& theta) {
    TGraph ff_gr;
    Int_t ff_i(0);
    Int_t np = gg_gr.GetN();
    Double_t *x = gg_gr.GetX();
    Double_t *y = gg_gr.GetY();
    for (Int_t i=0; i<np; i++) {
        if( fabs(theta - TMath::Sqrt(x[i]*x[i]+y[i]*y[i]))<0.05) {
            ff_gr.SetPoint(ff_i,x[i],y[i]);
            ff_i++;
        }
    }
    gg_gr = ff_gr;
    //Fit a circle to the graph points
    TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
    fitter->SetPrecision(0.00000001);
    fitter->SetMaxIterations(1000);
    fitter->SetFCN(circleFcn);
    fitter->SetParameter(0, "x0",   0.03, 0.01, -0.05,0.05);
    fitter->SetParameter(1, "y0",   0, 0.01, -0.05,0.05);
    fitter->SetParameter(2, "R",    theta, 0.01, theta-0.05,theta+0.05);
    //fitter->FixParameter(0);
    //fitter->FixParameter(1);
    fitter->FixParameter(2);
    Double_t arglist[1] = {0};
    fitter->ExecuteCommand("MINIMIZE", arglist, 0);
    // fitter->SetFCN(circleFcn2);
    // fitter->ExecuteCommand("MINIMIZE", arglist, 0);
    x0 = fitter->GetParameter(0);
    y0 = fitter->GetParameter(1);
    theta = fitter->GetParameter(2);
}

Int_t PrtLutReco::FindPdg(Double_t mom, Double_t cangle) {
    // Int_t pdg[]={11,13,211,321,2212};
    // Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
    // Int_t pdg[]={211,321,2212};
    // Double_t mass[] = {0.139570,0.49368,0.9382723};
    cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
    Int_t pdg[]= {211,2212};
    Double_t mass[] = {0.139570,0.9382723};
    Double_t tdiff, diff=100;
    Int_t minid=0;
    for(Int_t i=0; i<2; i++) {
        tdiff = fabs(cangle - acos(sqrt(mom*mom + mass[i]*mass[i])/mom/1.46907)); //1.46907 - fused silica
        if(tdiff<diff) {
            diff = tdiff;
            minid = i;
        }
    }
    return pdg[minid];
}






























