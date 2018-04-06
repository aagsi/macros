#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRotation.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TLegend.h>
//#include "/u/aali/dirc/prttools/prttools.C"
#include "/Users/ahmed/dirc/prttools/prttools.C"
#define prt__sim
#include "/Users/ahmed/dirc/prtdirc/src/PrtHit.h"
#include "/Users/ahmed/dirc/prtdirc/src/PrtEvent.h"
#include "/Users/ahmed/dirc/prttools/datainfo.C"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TList.h"
#include "TGraph.h"
#include "THStack.h"
#include "TLatex.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TFrame.h"
#include "TChain.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TAxis.h"
#define PI 3.14159265

// Constants
Double_t beam_momentum= 7000; // MeV
Double_t m_proton = 938.28;
Double_t m_pi = 139.570;
Double_t nano_value = 0.000000001; // 10e-9
Double_t c = 299792458; // speed of light
Double_t measured_d_tof2tof1_plot1_co= 28.507;
// calculation of protons and pions velocity based on the momentum
Double_t v_proton = c* sqrt(1- (m_proton*m_proton) /(beam_momentum * beam_momentum+ m_proton*m_proton));
Double_t v_pi = c* sqrt(1-( m_pi*m_pi) /(beam_momentum * beam_momentum+ m_pi*m_pi));

//TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
Double_t momentum=7.0;
Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
Double_t fAngleP = acos(sqrt(momentum*momentum+ mass[4]*mass[4])/momentum/1.4738)-0.00;
Double_t fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00;

////////////////////
// proto types//////
////////////////////

void HistoStyle(TH1F *x=new TH1F(), TH1F *y=new TH1F(), TH1F *z=new TH1F());
void GraphStyle(TGraph *x=new TGraph(), TGraph *y=new TGraph(), TGraph *z=new TGraph());
// file existance
bool exists_test (const std::string& name);

Bool_t Bool_cherenkov_correction(true), Bool_separation_graph(false), Bool_photonyield_histo(false), Bool_cherenkov_PDF_histo(false);

Bool_t bool_sim(true);
Bool_t bool_data(false);

Double_t max_min_array[prt_nmcp];
Double_t max_digi(5);
Double_t min_digi(-5);

Double_t cherenkov_true_mean_p=0;
Double_t cherenkov_true_sigma_p=0;
Double_t cherenkov_true_mean_pi=0;
Double_t cherenkov_true_sigma_pi=0;

Bool_t bool_plot(false);


////////////////////
// function   //////
////////////////////
//root -b -q ../cherenkov_angle_correction.C'(0,20)' for data
//root -b -q ../cherenkov_angle_correction.C'(1,20)' for simulation
// you can run it using a script
void cherenkov_angle_correction(Int_t flag = 1, Int_t angle = 20) {
    prt_initDigi();
    gStyle->SetPalette(kGreenPink);
    //gStyle->SetPalette(kBird);
    TString nid = Form("_%2.0d", 150);
    prt_savepath="pdf";
    //std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    TFile *ffile_data_p, *ffile_data_pi;
    TGraph *fHistCh_graph_p[960], *fHistCh_graph_pi[960];
    TH1F *fHistCh_read_p[960], *fHistCh_read_pi[960], *p_cherenkov_data ;
    TH1F *hist_nph_wo_p,  *hist_nph_wt_p, *hist_nph_wtc_p, *hist_nph_wo_pi,  *hist_nph_wt_pi, *hist_nph_wtc_pi;
    
    THStack *hs_nph_p, *hs_nph_pi;
    THStack * hs_p = new THStack("hs_p","Stacked 1D histograms");
    THStack * hs_pi = new THStack("hs_pi","Stacked 1D histograms");
    
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);
    
    ///////////////////////
    // Separation power  //
    ///////////////////////
    int counter =0 ;
    TGraph *power_org = new TGraph();
    TH1F *HistMcp_pi[12], *HistMcp_p[12], *HistMcp_same_path_p[12], *HistMcp_same_path_pi[12];
    TH1F *HistMcp_pi_diff[12], *HistMcp_p_diff[12];
    TFile *ffile_data_pi_spr, *ffile_data_p_spr;
    TFile *ffile_sim_pi_spr, *ffile_sim_p_spr;
    
    TH2F * twoD_mcp =  new TH2F("twoD_mcp",";MCP number [#] ;Polar angle [degree]", 12, 0, 12, 14, 20, 160);
    TH2F * twoD_mcp_true =  new TH2F("twoD_mcp_true",";MCP number [#] ;Polar angle [degree]", 12, 0, 12, 14, 20, 160);
    
    TH2F * twoD_mcp_pi =  new TH2F("twoD_mcp_pi",";MCP number [#] ;Polar angle [degree]", 12, 0, 12, 14, 20, 160);
    TH2F * twoD_mcp_p =  new TH2F("twoD_mcp_p",";MCP number [#] ;Polar angle [degree]", 12, 0, 12, 14, 20, 160);
    
    TH2F * twoD_angle_one_bin_p =  new TH2F("twoD_angle_one_bin_p",";MCP number [#] ;All polar angle [degree]", 12, 0, 12, 1, 0, 1);
    TH2F * twoD_angle_one_bin_pi =  new TH2F("twoD_angle_one_bin_pi",";MCP number [#] ;All polar angle [degree]", 12, 0, 12, 1, 0, 1);
    TH2F * twoD_angle_one_bin_all =  new TH2F("twoD_angle_one_bin_all",";MCP number [#] ;All polar angleangle [degree]", 12, 0, 12, 1, 0, 1);
    
    TH1F *oneD_angle_one_bin_p =new TH1F("oneD_angle_one_bin_p",";MCP number [#] ;#sum #Delta#theta_{c} [mrad]", 12, 0, 12);
    TH1F *oneD_angle_one_bin_pi =new TH1F("oneD_angle_one_bin_pi",";MCP number [#] ;#sum #Delta#theta_{c} [mrad]", 12, 0, 12);
    
    TH1F *oneD_mcp_one_bin_p_p =new TH1F("oneD_mcp_one_bin_p_p",";MCP number [#] ;#sum #Delta#theta_{c} [mrad]", 14, 20, 160);
    TH1F *oneD_mcp_one_bin_p_pi =new TH1F("oneD_mcp_one_bin_p_pi",";Polar angle [degree] ;#sum #Delta#theta_{c} [mrad]", 14, 20, 160);
    
    
    TH1F *mcp_all_oneD=new TH1F("mcp_all_oneD",";MCP number [#] ;#Delta#theta_{c} p - #Delta#theta_{c} #pi [mrad]", 40, -10, 10);
    TH1F *mcp_all_oneD_true=new TH1F("mcp_all_oneD_true",";#Delta#theta_{c}p - #Delta#theta_{c}#pi [mrad]",40, -10, 10);
    
    Int_t count_plot=0;
    for (int i=20; i<=150; i+=10) { // commint
    //{ // uncommint
      //  int i = angle; // uncommint
        TString spr_data_p_path,spr_data_pi_path;
        TString spr_sim_p_path,spr_sim_pi_path;
        if (flag==1) spr_data_p_path = Form("/Users/ahmed/dirc/cherenkov_correction/%d_sph_proton_sim_spr.root", i);
        if (flag==1) spr_data_pi_path = Form("/Users/ahmed/dirc/cherenkov_correction/%d_sph_pi_sim_spr.root", i);
        if (flag==0)spr_data_p_path = Form("/Users/ahmed/dirc/cherenkov_correction/%d_test_p_data_spr.root", i);
        if (flag==0)spr_data_pi_path = Form("/Users/ahmed/dirc/cherenkov_correction/%d_test_pi_data_spr.root", i);
        spr_sim_p_path = Form("/Users/ahmed/dirc/cherenkov_correction/%d_sph_proton_sim_spr.root", i);
        spr_sim_pi_path = Form("/Users/ahmed/dirc/cherenkov_correction/%d_sph_pi_sim_spr.root", i);
        TString separation_data_path = Form("/Users/ahmed/dirc/cherenkov_correction/%d_sph_data_separation.root", i);
//        cout<<"separation data path= " <<separation_data_path<<endl;
//        cout<<"spr data p path= " <<spr_data_p_path<<endl;
//        cout<<"spr data pi path= " <<spr_data_pi_path<<endl;
//        cout<<"spr sim p path= " <<spr_sim_p_path<<endl;
//        cout<<"spr sim pi path= " <<spr_sim_pi_path<<endl;
        string path_data_separation = (string)separation_data_path;
        string path_data_p_spr = (string)spr_data_p_path;
        string path_data_pi_spr = (string)spr_data_pi_path;
        string path_sim_p_spr = (string)spr_sim_p_path;
        string path_sim_pi_spr = (string)spr_sim_pi_path;
//        cout<<"exists test(separation path data)" <<exists_test(path_data_separation)<<endl;
//        cout<<"exists test(spr path data p)" <<exists_test(path_data_p_spr)<<endl;
//        cout<<"exists test(spr path data pi)" <<exists_test(path_data_pi_spr)<<endl;
//        cout<<"exists test(spr path sim p)" <<exists_test(path_sim_p_spr)<<endl;
//        cout<<"exists test(spr path sim pi)" <<exists_test(path_sim_pi_spr)<<endl;
        
        ////////////////
        // READ Tree ///
        ////////////////
        Double_t separation(-9), separation_org(-1);
        TChain ch("dirc");
        ch.Add(separation_data_path);
        ch.SetBranchAddress("separation",&separation);
        Int_t nent = ch.GetEntries();
        ch.GetEvent(nent-1);
        separation_org=separation;
        //std::cout<<"############  separation = "<< separation  <<std::endl;
        power_org->SetPoint(counter,i,separation_org);
        counter++;
        if(Bool_cherenkov_correction) {
            ///////////////////////////
            // CH MCP by MCP  Histo ///
            ///////////////////////////
            ffile_data_p_spr  = new TFile(spr_data_p_path, "READ");
            ffile_data_pi_spr  = new TFile(spr_data_pi_path, "READ");
            ffile_sim_p_spr  = new TFile(spr_sim_p_path, "READ");
            ffile_sim_pi_spr  = new TFile(spr_sim_pi_path, "READ");
            for(Int_t mcp=0; mcp<prt_nmcp; mcp++) {
                HistMcp_p[mcp] =(TH1F*)ffile_data_p_spr->Get(Form("fHistMcp_%d",mcp));
                HistMcp_pi[mcp] =(TH1F*)ffile_data_pi_spr->Get(Form("fHistMcp_%d",mcp));
                
                HistMcp_same_path_p[mcp] =(TH1F*)ffile_sim_p_spr->Get(Form("fHistMcp_same_path_%d",mcp));
                HistMcp_same_path_pi[mcp] =(TH1F*)ffile_sim_pi_spr->Get(Form("fHistMcp_same_path_%d",mcp));
            }
            for(Int_t mcp=0; mcp<prt_nmcp; mcp++) {
                if(mcp==1 && i ==20) continue;
                if(mcp==9 && i ==20) continue;
                if(mcp==10 && i ==20) continue;
                if(mcp==11 && i ==20) continue;
                //////////////////////////////
                if(mcp==6 && i ==30) continue;
                if(mcp==7 && i ==30) continue;
                if(mcp==8 && i ==30) continue;
                if(mcp==9 && i ==30) continue;
                if(mcp==11 && i ==30) continue;
                //////////////////////////////
                if(mcp==4 && i ==40) continue;
                if(mcp==9 && i ==40) continue;
                if(mcp==11 && i ==40) continue;
                /////////////////////////////
                if(mcp==9 && i ==50) continue;
                if(mcp==11 && i ==50) continue;
                ///////////////////////////////
                if(mcp==0 && i ==70) continue;
                if(mcp==1 && i ==70) continue;
                if(mcp==2 && i ==70) continue;
                //////////////////////////////
                if(mcp==0 && i ==80) continue;
                if(mcp==1 && i ==80) continue;
                if(mcp==2 && i ==80) continue;
                if(mcp==3 && i ==80) continue;
                if(mcp==4 && i ==80) continue;
                if(mcp==5 && i ==80) continue;
                if(mcp==7 && i ==80) continue;
                //////////////////////////////
                if(mcp==0 && i ==90) continue;
                if(mcp==1 && i ==90) continue;
                if(mcp==2 && i ==90) continue;
                if(mcp==3 && i ==90) continue;
                if(mcp==4 && i ==90) continue;
                if(mcp==5 && i ==90) continue;
                //////////////////////////////
                if(mcp==0 && i ==100) continue;
                if(mcp==1 && i ==100) continue;
                if(mcp==2 && i ==100) continue;
                if(mcp==3 && i ==100) continue;
                if(mcp==4 && i ==100) continue;
                if(mcp==5 && i ==100) continue;
                ///////////////////////////////
                if(mcp==0 && i ==110) continue;
                if(mcp==1 && i ==110) continue;
                if(mcp==2 && i ==110) continue;
                if(mcp==3 && i ==110) continue;
                if(mcp==4 && i ==110) continue;
                if(mcp==5 && i ==110) continue;
                //////////////////////////////
                if(mcp==9 && i ==130) continue;
                if(mcp==10 && i ==130) continue;
                if(mcp==11 && i ==130) continue;
                ///////////////////////////////
                if(mcp==4 && i ==140) continue;
                if(mcp==9 && i ==140) continue;
                if(mcp==11 && i ==140) continue;
                ///////////////////////////////
                if(mcp==9 && i ==150) continue;
                if(mcp==11 && i ==150) continue;
                ///////////////////////////////
                
                // Truth correction
                // initialize fitting functions for proton truth
                //cout << "MCP by MCP Truth"<< endl;
                TF1 *Fit_MCP_truth_p = new TF1("Fit_MCP_truth_p","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
                Fit_MCP_truth_p->SetParameters(100,fAngleP,0.010);
                Fit_MCP_truth_p->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
                Fit_MCP_truth_p->SetParLimits(0,0.1,1E6);
                Fit_MCP_truth_p->SetParLimits(1,fAngleP-0.02,fAngleP+0.02);
                Fit_MCP_truth_p->SetParLimits(2,0.004,0.008);
                Fit_MCP_truth_p->SetLineColor(kGreen);
                
                // initialize fitting functions for pi truth
                TF1 *Fit_MCP_truth_pi = new TF1("Fit_MCP_truth_pi","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
                Fit_MCP_truth_pi->SetParameters(100,fAnglePi,0.010);
                Fit_MCP_truth_pi->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
                Fit_MCP_truth_pi->SetParLimits(0,0.1,1E6);
                Fit_MCP_truth_pi->SetParLimits(1,fAnglePi-0.02,fAnglePi+0.02);
                Fit_MCP_truth_pi->SetParLimits(2,0.004,0.008);
                Fit_MCP_truth_pi->SetLineColor(kGreen);
                
                // plot truth p
                count_plot++;
                //if(bool_plot)prt_canvasAdd(Form("r_mcp_p_true_sim_%d_prtangle_%d",mcp,i),800,400);
                if(bool_plot)prt_canvasAdd(Form("r_mcp_%d",count_plot),800,400);
                HistMcp_same_path_p[mcp]->SetTitle(Form("p truth sim mcp %d prtangle %d",mcp,i));
                //HistMcp_same_path_p[mcp]->SetLineColor(kGreen);
                HistMcp_same_path_p[mcp]->Fit("Fit_MCP_truth_p","lq","",fAngleP-0.025,fAngleP+0.025);
                if(bool_plot){
                    HistMcp_same_path_p[mcp]->Draw();
                    //prt_canvasGet(Form("r_mcp_p_true_sim_%d_prtangle_%d",mcp,i))->Update();
                    prt_canvasGet(Form("r_mcp_%d",count_plot))->Update();
                    TLine *lin_ch_p_v_truth = new TLine(0,0,0,1000);
                    lin_ch_p_v_truth->SetX1(fAngleP);
                    lin_ch_p_v_truth->SetX2(fAngleP);
                    lin_ch_p_v_truth->SetY1(gPad->GetUymin());
                    lin_ch_p_v_truth->SetY2(gPad->GetUymax());
                    lin_ch_p_v_truth->SetLineColor(kRed);
                    lin_ch_p_v_truth->Draw();
                }
                // plot truth pi
                count_plot++;
                //if(bool_plot)prt_canvasAdd(Form("r_mcp_pi_true_sim_%d_prtangle_%d",mcp,i),800,400);
                if(bool_plot)prt_canvasAdd(Form("r_mcp_%d",count_plot),800,400);
                HistMcp_same_path_pi[mcp]-> SetTitle(Form("#pi truth sim mcp %d prtangle %d",mcp,i));
                //HistMcp_same_path_pi[mcp]->SetLineColor(kGreen);
                HistMcp_same_path_pi[mcp]->Fit("Fit_MCP_truth_pi","lq","",fAnglePi-0.025,fAnglePi+0.025);
                if(bool_plot){
                    HistMcp_same_path_pi[mcp]->Draw();
                    //prt_canvasGet(Form("r_mcp_pi_true_sim_%d_prtangle_%d",mcp,i))->Update();
                    prt_canvasGet(Form("r_mcp_%d",count_plot))->Update();
                    TLine *lin_ch_pi_v_truth = new TLine(0,0,0,1000);
                    lin_ch_pi_v_truth->SetX1(fAnglePi);
                    lin_ch_pi_v_truth->SetX2(fAnglePi);
                    lin_ch_pi_v_truth->SetY1(gPad->GetUymin());
                    lin_ch_pi_v_truth->SetY2(gPad->GetUymax());
                    lin_ch_pi_v_truth->SetLineColor(kBlue);
                    lin_ch_pi_v_truth->Draw();
                }
                // calculate ranges
                cherenkov_true_mean_p = Fit_MCP_truth_p->GetParameter(1);
                cherenkov_true_sigma_p = Fit_MCP_truth_p->GetParameter(2);
                Double_t cherenkov_p_minus_5_sgma = cherenkov_true_mean_p-5*cherenkov_true_sigma_p;
                Double_t cherenkov_p_plus_5_sgma = cherenkov_true_mean_p+5*cherenkov_true_sigma_p;
                Double_t cherenkov_p_minus_3_sgma = cherenkov_true_mean_p-3*cherenkov_true_sigma_p;
                Double_t cherenkov_p_plus_3_sgma = cherenkov_true_mean_p+3*cherenkov_true_sigma_p;
                Double_t cherenkov_p_minus_2_sgma = cherenkov_true_mean_p-2*cherenkov_true_sigma_p;
                Double_t cherenkov_p_plus_2_sgma = cherenkov_true_mean_p+2*cherenkov_true_sigma_p;
                
                cherenkov_true_mean_pi = Fit_MCP_truth_pi->GetParameter(1);
                cherenkov_true_sigma_pi = Fit_MCP_truth_pi->GetParameter(2);
                Double_t cherenkov_pi_minus_5_sgma = cherenkov_true_mean_pi-5*cherenkov_true_sigma_pi;
                Double_t cherenkov_pi_plus_5_sgma = cherenkov_true_mean_pi+5*cherenkov_true_sigma_pi;
                Double_t cherenkov_pi_minus_3_sgma = cherenkov_true_mean_pi-3*cherenkov_true_sigma_pi;
                Double_t cherenkov_pi_plus_3_sgma = cherenkov_true_mean_pi+3*cherenkov_true_sigma_pi;
                Double_t cherenkov_pi_minus_2_sgma = cherenkov_true_mean_pi-2*cherenkov_true_sigma_pi;
                Double_t cherenkov_pi_plus_2_sgma = cherenkov_true_mean_pi+2*cherenkov_true_sigma_pi;
                
                // initialize fitting functions for Proton
                TF1 *Fit_MCP_p = new TF1("Fit_MCP_p","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
                //Double_t max_bin_p =  HistMcp_p[mcp]->GetXaxis()->GetBinCenter(HistMcp_p[mcp]->GetMaximumBin());
                //Fit_MCP_p->SetParameters(100,max_bin_p,0.010);
                Fit_MCP_p->SetParameters(100,fAngleP,0.010);
                Fit_MCP_p->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
                Fit_MCP_p->SetParLimits(0,0.1,1E6);
                Fit_MCP_p->SetParLimits(1,cherenkov_true_mean_p-0.02,cherenkov_true_mean_p+0.02);
                Fit_MCP_p->SetParLimits(2,0.004,0.008);
                
                // initialize fitting functions for pi
                TF1 *Fit_MCP_pi = new TF1("Fit_MCP_pi","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
                //Double_t max_bin_pi =  HistMcp_pi[mcp]->GetXaxis()->GetBinCenter(HistMcp_pi[mcp]->GetMaximumBin());
                //Fit_MCP_p->SetParameters(100,max_bin_pi,0.010);
                Fit_MCP_pi->SetParameters(100,fAnglePi,0.010);
                Fit_MCP_pi->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
                Fit_MCP_pi->SetParLimits(0,0.1,1E6);
                Fit_MCP_pi->SetParLimits(1,cherenkov_true_mean_pi-0.02,cherenkov_true_mean_pi+0.02);
                Fit_MCP_pi->SetParLimits(2,0.004,0.008);
                Fit_MCP_pi->SetLineColor(kBlue);
                
                // plot p
                count_plot++;
                //if(bool_plot)prt_canvasAdd(Form("r_mcp_p_%d_prtangle_%d",mcp,i),800,400);
                if(bool_plot)prt_canvasAdd(Form("r_mcp_%d",count_plot),800,400);
                HistMcp_p[mcp]-> SetTitle(Form("p data mcp %d prtangle %d",mcp,i));
                
                if(i==70 && mcp == 5){
                    HistMcp_p[mcp]->Fit("Fit_MCP_p","lq","",cherenkov_p_minus_3_sgma, cherenkov_p_plus_3_sgma);
                }else if(i== 150 && mcp == 8 && flag==1){
                    HistMcp_p[mcp]->Fit("Fit_MCP_p","lq","",cherenkov_p_minus_3_sgma, cherenkov_p_plus_3_sgma);
                    
                }else{
                    HistMcp_p[mcp]->Fit("Fit_MCP_p","lq","",cherenkov_p_minus_5_sgma, cherenkov_p_plus_5_sgma);
                }
                //HistMcp_p[mcp]->Fit("Fit_MCP_p","lq","",cherenkov_pi_minus_3_sgma, cherenkov_pi_plus_3_sgma);
                if(bool_plot){
                    HistMcp_p[mcp]->Draw();
                    //prt_canvasGet(Form("r_mcp_p_%d_prtangle_%d",mcp,i))->Update();
                    prt_canvasGet(Form("r_mcp_%d",count_plot))->Update();
                    TLine *lin_ch_p_v = new TLine(0,0,0,1000);
                    lin_ch_p_v->SetX1(fAngleP);
                    lin_ch_p_v->SetX2(fAngleP);
                    lin_ch_p_v->SetY1(gPad->GetUymin());
                    lin_ch_p_v->SetY2(gPad->GetUymax());
                    lin_ch_p_v->SetLineColor(kRed);
                    lin_ch_p_v->Draw();
                }
                //cout << "MCP by MCP cherenvove angle correction for Proton"<< endl;
                std::cout<<"if(mcpid=="<< mcp<<") tangle += "<<fAngleP-Fit_MCP_p->GetParameter(1)<<";" <<std::endl;
                
                
                // plot pi
                count_plot++;
                //if(bool_plot)prt_canvasAdd(Form("r_mcp_pi_%d_prtangle_%d",mcp,i),800,400);
                if(bool_plot)prt_canvasAdd(Form("r_mcp_%d",count_plot),800,400);
                HistMcp_pi[mcp]-> SetTitle(Form("#pi data mcp %d prtangle %d",mcp,i));
                if(i==100 && mcp == 11){
                    HistMcp_pi[mcp]->Fit("Fit_MCP_pi","lq","",cherenkov_pi_minus_3_sgma, cherenkov_pi_plus_3_sgma);
                }else if((i==140 && mcp == 10 && flag ==1) || (i== 150 && mcp == 8 && flag ==1)){
                    HistMcp_pi[mcp]->Fit("Fit_MCP_pi","lq","",cherenkov_pi_minus_3_sgma, cherenkov_pi_plus_3_sgma);
                    
                }else{
                    HistMcp_pi[mcp]->Fit("Fit_MCP_pi","lq","",cherenkov_pi_minus_5_sgma, cherenkov_pi_plus_5_sgma);
                }
                
                //HistMcp_pi[mcp]->Fit("Fit_MCP_pi","lq","",cherenkov_pi_minus_3_sgma, cherenkov_pi_plus_3_sgma);
                
                if(bool_plot){
                    HistMcp_pi[mcp]->Draw();
                    //prt_canvasGet(Form("r_mcp_pi_%d_prtangle_%d",mcp,i))->Update();
                    prt_canvasGet(Form("r_mcp_%d",count_plot))->Update();
                    TLine *lin_ch_pi_v = new TLine(0,0,0,1000);
                    lin_ch_pi_v->SetX1(fAnglePi);
                    lin_ch_pi_v->SetX2(fAnglePi);
                    lin_ch_pi_v->SetY1(gPad->GetUymin());
                    lin_ch_pi_v->SetY2(gPad->GetUymax());
                    lin_ch_pi_v->SetLineColor(kBlue);
                    lin_ch_pi_v->Draw();
                }
                //cout << "MCP by MCP cherenvove angle correction for PION"<< endl;
                //std::cout<<"if(mcpid=="<< mcp<<") tangle += "<<fAnglePi-Fit_MCP_pi->GetParameter(1)<<";" <<std::endl;
                
                Double_t delta_theta=fAngleP-Fit_MCP_p->GetParameter(1) - (fAnglePi-Fit_MCP_pi->GetParameter(1));
                Double_t val_1 = (fAngleP-Fit_MCP_p->GetParameter(1))*1000.0 ;
                Double_t val_2 = (fAnglePi-Fit_MCP_pi->GetParameter(1))*1000.0;
                Double_t val_3 = delta_theta*1000.0;
                Double_t delta_theta_true=fAngleP-Fit_MCP_truth_p->GetParameter(1) - (fAnglePi-Fit_MCP_truth_pi->GetParameter(1));
                Double_t val_3_true = delta_theta_true*1000.0  ;
                
                twoD_mcp_p->Fill(mcp, i,val_1);
                twoD_mcp_pi->Fill(mcp, i,val_2);
                twoD_mcp->Fill(mcp, i,val_3);
                twoD_mcp_true->Fill(mcp, i,val_3_true);
                mcp_all_oneD->Fill(val_3);
                mcp_all_oneD_true->Fill(val_3_true);
                oneD_angle_one_bin_p->Fill(mcp, val_1);
                oneD_angle_one_bin_pi->Fill(mcp, val_2);
                oneD_mcp_one_bin_p_p->Fill(i, val_1);
                oneD_mcp_one_bin_p_pi->Fill(i, val_2);
                //for digi max and min
                max_min_array[mcp]=val_1;
                
                for(Int_t m=0; m<64; m++)
                    for(Int_t n=0; n<64; n++){
                        prt_hdigi[mcp]->Fill(m,n, val_1);
                    }
            }
        }
    }
    // open plots
    if(Bool_cherenkov_correction) {
        prt_canvasAdd("r_shift_oneD",800,400);
        prt_canvasGet("r_shift_oneD")->SetGridx();
        TLegend * legend_mcp_all_oneD= new TLegend(0.576441,0.578667,0.966165,0.858667);
        legend_mcp_all_oneD->AddEntry(mcp_all_oneD_true,"#Delta#theta_{c}p - #Delta#theta_{c}#pi (true path inside prism)","l");
        legend_mcp_all_oneD->AddEntry(mcp_all_oneD,"#Delta#theta_{c}p - #Delta#theta_{c}#pi (Sim) ","l");
        mcp_all_oneD_true->SetStats(0);
        mcp_all_oneD_true->SetLineColor(kGreen);
        mcp_all_oneD_true->SetLineWidth(3);
        mcp_all_oneD->SetLineWidth(3);
        mcp_all_oneD_true->Draw();
        mcp_all_oneD->Draw("same");
        legend_mcp_all_oneD->Draw();
        
        if(false){
            prt_canvasAdd("r_shift_all",800,400);
            prt_canvasGet("r_shift_all")->SetGridx();
            prt_canvasGet("r_shift_all")->SetGridy();
            twoD_mcp->GetXaxis()->SetNdivisions(15);
            twoD_mcp->GetYaxis()->SetNdivisions(15);
            twoD_mcp->SetMaximum(11.5);
            twoD_mcp->SetMinimum(-8);
            prt_canvasGet("r_shift_all")->Update();
            twoD_mcp->SetStats(0);
            twoD_mcp-> SetTitle("MCP by MCP (#Delta#theta_{C}p - #Delta#theta_{C} #pi [mrad])" );
            twoD_mcp->SetMarkerSize(1.8);
            twoD_mcp->Draw("colztext");
            
            prt_canvasAdd("r_shift_all_true",800,400);
            prt_canvasGet("r_shift_all_true")->SetGridx();
            prt_canvasGet("r_shift_all_true")->SetGridy();
            twoD_mcp_true->GetXaxis()->SetNdivisions(15);
            twoD_mcp_true->GetYaxis()->SetNdivisions(15);
            twoD_mcp_true->SetMaximum(11.5);
            twoD_mcp_true->SetMinimum(-8);
            prt_canvasGet("r_shift_all_true")->Update();
            twoD_mcp_true->SetStats(0);
            twoD_mcp_true-> SetTitle("MCP by MCP true (#Delta#theta_{C}p - #Delta#theta_{C} #pi [mrad])" );
            twoD_mcp_true->SetMarkerSize(1.8);
            twoD_mcp_true->Draw("colztext");
            
            prt_canvasAdd("r_shift_p",800,400);
            prt_canvasGet("r_shift_p")->SetGridx();
            prt_canvasGet("r_shift_p")->SetGridy();
            twoD_mcp_p->GetXaxis()->SetNdivisions(15);
            twoD_mcp_p->GetYaxis()->SetNdivisions(15);
            twoD_mcp_p->SetMaximum(13.5);
            twoD_mcp_p->SetMinimum(-8);
            prt_canvasGet("r_shift_p")->Update();
            twoD_mcp_p->SetStats(0);
            twoD_mcp_p-> SetTitle("MCP by MCP #Delta#theta_{C}p [mrad]" );
            twoD_mcp_p->SetMarkerSize(1.8);
            twoD_mcp_p->Draw("colztext");
            
            prt_canvasAdd("r_shift_pi",800,400);
            prt_canvasGet("r_shift_pi")->SetGridx();
            prt_canvasGet("r_shift_pi")->SetGridy();
            twoD_mcp_pi->GetXaxis()->SetNdivisions(15);
            twoD_mcp_pi->GetYaxis()->SetNdivisions(15);
            twoD_mcp_pi->SetMaximum(13.5);
            twoD_mcp_pi->SetMinimum(-8);
            prt_canvasGet("r_shift_pi")->Update();
            twoD_mcp_pi->SetStats(0);
            twoD_mcp_pi-> SetTitle("MCP by MCP #Delta#theta_{C} #pi [mrad]" );
            //twoD_mcp_pi->SetMarkerColor(kRed);
            twoD_mcp_pi->SetMarkerSize(1.8);
            gStyle->SetPaintTextFormat("4.1f");
            twoD_mcp_pi->Draw("colztext");
            
            prt_canvasAdd("r_oneD_angle_one_bin",800,400);
            prt_canvasGet("r_oneD_angle_one_bin")->SetGridx();
            prt_canvasGet("r_oneD_angle_one_bin")->SetGridy();
            TLegend * legend_thetac_correction= new TLegend(0.12, 0.58,  0.33,   0.87  );
            //legend_thetac_correction->SetHeader("#pi and P #theta_{c} correction","C");
            legend_thetac_correction->AddEntry(oneD_mcp_one_bin_p_p,"P #sum #Delta#theta_{c} correction","l");
            legend_thetac_correction->AddEntry(oneD_mcp_one_bin_p_pi,"#pi #sum #Delta#theta_{c} correction ","l");
            oneD_angle_one_bin_p-> SetTitle("#sum #Delta#theta_{c} correction all angles" );
            oneD_angle_one_bin_pi-> SetTitle("#sum #Delta#theta_{c} correction all angles" );
            oneD_angle_one_bin_p->SetLineColor(kRed);
            oneD_angle_one_bin_p->SetLineStyle(1);
            oneD_angle_one_bin_p->SetLineWidth(3);
            oneD_angle_one_bin_pi->SetLineWidth(3);
            oneD_angle_one_bin_p->SetMarkerStyle(20);
            oneD_angle_one_bin_pi->SetLineColor(kBlue);
            oneD_angle_one_bin_pi->SetLineStyle(1);
            oneD_angle_one_bin_pi->SetMarkerStyle(20);
            oneD_angle_one_bin_pi->SetStats(0);
            oneD_angle_one_bin_p->SetStats(0);
            oneD_angle_one_bin_p->SetTitle("MCP by MCP #Delta#theta_{C}" );
            oneD_angle_one_bin_pi->Draw("histo");
            oneD_angle_one_bin_p->Draw("histosame");
            legend_thetac_correction->Draw();
            
            prt_canvasAdd("r_oneD_mcp_one_bin",800,400);
            prt_canvasGet("r_oneD_mcp_one_bin")->SetGridx();
            prt_canvasGet("r_oneD_mcp_one_bin")->SetGridy();
            oneD_mcp_one_bin_p_p-> SetTitle("#sum #Delta#theta_{c} correction all MCP's" );
            oneD_mcp_one_bin_p_pi-> SetTitle("#sum #Delta#theta_{c} correction all MCP's" );
            oneD_mcp_one_bin_p_p->SetLineColor(kRed);
            oneD_mcp_one_bin_p_p->SetLineStyle(1);
            oneD_mcp_one_bin_p_p->SetMarkerStyle(20);
            oneD_mcp_one_bin_p_pi->SetLineColor(kBlue);
            oneD_mcp_one_bin_p_pi->SetLineWidth(3);
            oneD_mcp_one_bin_p_p->SetLineWidth(3);
            oneD_mcp_one_bin_p_pi->SetLineStyle(1);
            oneD_mcp_one_bin_p_pi->SetMarkerStyle(20);
            
            oneD_mcp_one_bin_p_pi->SetStats(0);
            oneD_mcp_one_bin_p_p->SetStats(0);
            oneD_mcp_one_bin_p_p->SetTitle("MCP by MCP #Delta#theta_{C}" );
            oneD_mcp_one_bin_p_pi->Draw("histo");
            oneD_mcp_one_bin_p_p->Draw("histosame");
            legend_thetac_correction->Draw();
            
            
            
            
            max_digi =max_min_array[0];
            min_digi =max_min_array[0];
            
            for (Int_t i = 0; i < prt_nmcp; i++)
            {
                if (max_min_array[i] > max_digi)
                {
                    max_digi = max_min_array[i];
                }
                else if (max_min_array[i] < min_digi)
                {
                    min_digi = max_min_array[i];
                }
            }
            //cout <<"#################"<< max_digi << endl;
            //cout <<"#################"<< min_digi << endl;
            max_digi=13.5;
            min_digi=-8.0;
            /*
             if(angle==20){
             max_digi=0.54;
             min_digi=-7.41;
             }
             if(angle==30){
             max_digi=2.06;
             min_digi=-5.69;
             }
             if(angle==40){
             max_digi=4.92;
             min_digi=-2.45;
             }
             if(angle==50){
             max_digi=3.37;
             min_digi=-4.21;
             }
             if(angle==60){
             max_digi=5.95;
             min_digi=-4.68;
             }
             if(angle==70){
             max_digi=7.05;
             min_digi=-4.89;
             }
             if(angle==80){
             max_digi=10;
             min_digi=-10;
             }
             if(angle==90){
             max_digi=9.54;
             min_digi=-4.19;
             }
             if(angle==100){
             max_digi=9;
             min_digi=-5;
             }
             if(angle==110){
             max_digi=7.38;
             min_digi=-4.29;
             }
             if(angle==120){
             max_digi=6.68;
             min_digi=-6.24;
             }
             if(angle==130){
             max_digi=6.31;
             min_digi=-4.03;
             }
             if(angle==140){
             max_digi=3.5;
             min_digi=-4.52;
             }
             if(angle==150){
             max_digi=11.5;
             min_digi=-6.5;
             }
             */
            prt_drawDigi("m,p,v\n",2017, max_digi,min_digi);
            prt_cdigi->SetName(Form("hp_dataProtonS332_%d_%d",flag, angle));
            prt_canvasAdd(prt_cdigi);
            prt_cdigi_palette->Draw();
            prt_canvasSave(2,0);
        }
    }
    
    
    GraphStyle(power_org);
    /*
     //~ power_org->Sort();
     //~ power_org->SetLineColor(kBlack);
     //~ power_org->SetMarkerColor(kBlack);
     //~ power_org->SetMarkerStyle(21);
     //~ power_org->SetMarkerSize(0.7);
     //~ power_org->SetName("separation");
     //~ power_org->GetYaxis()->SetRangeUser(0,20);
     //~ power_org->GetYaxis()->SetTitle("Separation [s.d]");
     //~ power_org->GetXaxis()->SetLabelSize(0.05);
     //~ power_org->GetXaxis()->SetTitleSize(0.06);
     //~ power_org->GetXaxis()->SetTitleOffset(0.84);
     
     if(Bool_separation_graph) {
     prt_canvasAdd("r_separation",800,400);
     TLegend *leg_separation = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
     leg_separation->SetHeader("separation power ","C");
     leg_separation->SetFillColor(0);
     leg_separation->AddEntry(power_org, "Cherenkove ambiguity PDF", "lp");
     //~ leg_separation->AddEntry(power_sim, "Sim without cherenkov angle correction", "lp");
     //~ leg_separation->AddEntry(power_p_shift, "p correction", "lp");
     //~ leg_separation->AddEntry(power_custom_shift, "custom correction", "lp");
     //~ leg_separation->AddEntry(power_modle_shift, "shift model", "lp");
     TMultiGraph *mg_separation = new TMultiGraph();
     mg_separation->Add(power_org);
     //~ mg_separation->Add(power_sim);
     //~ mg_separation->Add(power_modle_shift);
     //~ mg_separation->Add(power_p_shift);
     //~ mg_separation->Add(power_custom_shift);
     mg_separation->SetTitle(" separatin power geometrical reconstruction ;#theta [degree]; separation [s.d.]");
     mg_separation->Draw("APL");
     leg_separation->Draw();
     }
     
     ///////////////////////////////
     // READ files and histograms //
     ///////////////////////////////
     TString cherenkov_data_p_path = Form("/u/aali/work/test/histo_%d_sph_p_data_spr.root", 20);
     TString cherenkov_data_pi_path = Form("/u/aali/work/test/histo_%d_sph_pi_data_spr.root", 20);
     cout<<"cherenkov_data_p_path= " <<cherenkov_data_p_path<<endl;
     cout<<"cherenkov_data_pi_path= " <<cherenkov_data_pi_path<<endl;
     ffile_data_p  = new TFile(cherenkov_data_p_path, "READ");
     ffile_data_pi  = new TFile(cherenkov_data_pi_path, "READ");
     hs_nph_p = new THStack("hs_nph_p","Stacked 1D histograms");
     hs_nph_pi = new THStack("hs_nph_pi","Stacked 1D histograms");
     
     ////////////////////////////
     // photon yield histogram //
     ////////////////////////////
     hist_nph_wo_p=(TH1F*)ffile_data_p->Get("fnHits");
     hist_nph_wt_p=(TH1F*)ffile_data_p->Get("fnHits_p");
     hist_nph_wtc_p=(TH1F*)ffile_data_p->Get("fnHits_p_good");
     hist_nph_wo_pi=(TH1F*)ffile_data_pi->Get("fnHits");
     hist_nph_wt_pi=(TH1F*)ffile_data_pi->Get("fnHits_p");
     hist_nph_wtc_pi=(TH1F*)ffile_data_pi->Get("fnHits_p_good");
     hs_nph_p->Add(hist_nph_wo_p);
     hs_nph_p->Add(hist_nph_wt_p);
     hs_nph_p->Add(hist_nph_wtc_p);
     hs_nph_pi->Add(hist_nph_wo_pi);
     hs_nph_pi->Add(hist_nph_wt_pi);
     hs_nph_pi->Add(hist_nph_wtc_pi);
     HistoStyle(hist_nph_wo_p, hist_nph_wt_p, hist_nph_wtc_p);
     HistoStyle(hist_nph_wo_pi, hist_nph_wt_pi, hist_nph_wtc_pi);
     Double_t nph_pi=hist_nph_wtc_pi->GetEntries();
     Double_t nph_p=hist_nph_wtc_p->GetEntries();
     std::cout<<"nph_p=  "<<nph_p<<std::endl;
     std::cout<<"nph_pi=  "<<nph_pi<<std::endl;
     
     //////////////////////////
     // Pi & P photon Yield  //
     //////////////////////////
     if (Bool_photonyield_histo) {
     TLegend * legend_nph_data_p= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
     legend_nph_data_p->SetHeader("photon yield (proton data)","C");
     prt_canvasAdd("r_fnHits_data_p"+nid,800,400);
     hs_nph_p->SetTitle(Form("Polar angle %d (data)", 20));
     legend_nph_data_p->AddEntry(hist_nph_wo_p,"DIRC hits without cuts ","l");
     legend_nph_data_p->AddEntry(hist_nph_wt_p,"DIRC hits with #Deltat cut ","l");
     legend_nph_data_p->AddEntry(hist_nph_wtc_p,"DIRC hits with #theta_{C} and #Deltat cut ","l");
     
     hs_nph_p->SetTitle(Form("Polar angle %d , nph=%f", 20, nph_p));
     hs_nph_p->Draw("nostack");
     hs_nph_p->GetYaxis()->SetTitle("entries [#]");
     hs_nph_p->GetXaxis()->SetTitle("number of photon per track [#]");
     legend_nph_data_p->Draw();
     
     TLegend * legend_nph_data_pi= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
     legend_nph_data_pi->SetHeader("photon yield (pion data)","C");
     prt_canvasAdd("r_fnHits_data_pi"+nid,800,400);
     hs_nph_pi->SetTitle(Form("Polar angle %d (data)", 20));
     legend_nph_data_pi->AddEntry(hist_nph_wo_pi,"DIRC hits without cuts ","l");
     legend_nph_data_pi->AddEntry(hist_nph_wt_pi,"DIRC hits with #Deltat cut ","l");
     legend_nph_data_pi->AddEntry(hist_nph_wtc_pi,"DIRC hits with #theta_{C} and #Deltat cut ","l");
     
     hs_nph_pi->SetTitle(Form("Polar angle %d , nph=%f", 20, nph_pi));
     hs_nph_pi->Draw("nostack");
     hs_nph_pi->GetYaxis()->SetTitle("entries [#]");
     hs_nph_pi->GetXaxis()->SetTitle("number of photon per track [#]");
     legend_nph_data_pi->Draw();
     }
     
     //////////////////////////
     // cherenkov PDF Histo  //
     //////////////////////////
     if (Bool_cherenkov_PDF_histo) {
     for(Int_t pix=0; pix<960; pix++) {
     //fHistCh_graph_p[pix] =new TGraph(   (TH1F*)ffile_data_p->Get(Form("fHistCh_%d",pix))   );
     //fHistCh_graph_pi[pix] =new TGraph(   (TH1F*)ffile_data_pi->Get(Form("fHistCh_%d",pix))   );
     fHistCh_read_p[pix] = (TH1F*)ffile_data_p->Get(Form("fHistCh_%d",pix));
     fHistCh_read_pi[pix] = (TH1F*)ffile_data_pi->Get(Form("fHistCh_%d",pix));
     }
     TAxis *axis_data[960], *axis_data_pi[960];
     Int_t bmin_data[960],  bmax_data[960], bmin_data_pi[960], bmax_data_pi[960];
     Double_t integral_data[960], integral_data_pi[960];
     Double_t xmin_data = 0.6;
     Double_t xmax_data = 1.0;
     for(Int_t pix=700; pix<800; pix++) {
     TMultiGraph *mg = new TMultiGraph();
     TLegend *legend_ch_match= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
     legend_ch_match->SetHeader("Ambiguity distribution per pixle","C");
     prt_canvasAdd(Form("r_20_pix_%d",pix),800,400);
     axis_data[pix] = fHistCh_read_p[pix]->GetXaxis();
     //if (i==90)xmin_data  = 0.9;
     //if (i==90)xmax_data  = 1.0;
     bmin_data[pix] = axis_data[pix]->FindBin(xmin_data);
     bmax_data[pix] = axis_data[pix]->FindBin(xmax_data);
     integral_data[pix] = fHistCh_read_p[pix]->Integral(bmin_data[pix],bmax_data[pix]);
     //fHistCh_read_p[pix]->Scale(1/integral_data[pix]);
     fHistCh_read_p[pix]->Scale(1/nph_p);
     axis_data_pi[pix] = fHistCh_read_pi[pix]->GetXaxis();
     //std::cout<<"No problem 1  "<<std::endl;
     //if (i==90)xmin_data = 0.9;
     //if (i==90)xmax_data = 1.0;
     bmin_data_pi[pix] = axis_data_pi[pix]->FindBin(xmin_data);
     bmax_data_pi[pix] = axis_data_pi[pix]->FindBin(xmax_data);
     integral_data_pi[pix] = fHistCh_read_pi[pix]->Integral(bmin_data_pi[pix],bmax_data_pi[pix]);
     //fHistCh_read_pi[pix]->Scale(1/integral_data_pi[pix]);
     fHistCh_read_pi[pix]->Scale(1/nph_pi);
     fHistCh_read_pi[pix]->SetLineColor(kRed);
     fHistCh_graph_p[pix] =new TGraph(fHistCh_read_p[pix]);
     fHistCh_graph_pi[pix] =new TGraph(fHistCh_read_pi[pix]);
     fHistCh_graph_p[pix]->SetMarkerStyle(6);
     fHistCh_graph_p[pix]->SetMarkerSize(21);
     fHistCh_graph_p[pix]->SetLineColor(kRed);
     fHistCh_graph_p[pix]->SetLineStyle(3);
     fHistCh_graph_p[pix]->Sort();
     fHistCh_graph_pi[pix]->SetMarkerStyle(6);
     fHistCh_graph_pi[pix]->SetMarkerSize(21);
     fHistCh_graph_pi[pix]->SetLineColor(kBlue);
     fHistCh_graph_pi[pix]->SetLineStyle(3);
     fHistCh_graph_pi[pix]->Sort();
     mg->Add(fHistCh_graph_p[pix]);
     mg->Add(fHistCh_graph_pi[pix]);
     //mg->SetTitle("Ambiguity distribution at 20 polar angle  ;#theta_{C} [rad];  [#]");
     mg->SetTitle(Form("Ambiguity distribution at 20 polar angle pix %d  ;#theta_{C} [rad];  [#]", pix));
     mg->Draw("APL");
     fHistCh_read_pi[pix]->SetTitle(Form("Polar angle %3.1d", pix));
     fHistCh_read_pi[pix]->Draw("samehist");
     fHistCh_read_p[pix]->Draw("samehist");
     legend_ch_match->Draw();
     prt_canvasGet(Form("r_20_pix_%d",pix))->Update();
     TLine *lin_ch_p_v = new TLine(0,0,0,1000);
     lin_ch_p_v->SetX1(fAngleP);
     lin_ch_p_v->SetX2(fAngleP);
     lin_ch_p_v->SetY1(gPad->GetUymin());
     lin_ch_p_v->SetY2(gPad->GetUymax());
     lin_ch_p_v->SetLineColor(kRed);
     lin_ch_p_v->SetLineStyle(2);
     TLine *lin_ch_pi_v = new TLine(0,0,0,1000);
     lin_ch_pi_v->SetX1(fAnglePi);
     lin_ch_pi_v->SetX2(fAnglePi);
     lin_ch_pi_v->SetY1(gPad->GetUymin());
     lin_ch_pi_v->SetY2(gPad->GetUymax());
     lin_ch_pi_v->SetLineColor(kBlue);
     lin_ch_pi_v->SetLineStyle(2);
     lin_ch_p_v->Draw();
     lin_ch_pi_v->Draw();
     legend_ch_match->AddEntry(fHistCh_graph_p[pix],"Ambiguity distribution proton","l");
     legend_ch_match->AddEntry(fHistCh_graph_pi[pix]," Ambiguity distribution pion ","l");
     legend_ch_match->AddEntry(fHistCh_read_p[pix],"Histo Ambiguity distribution proton","l");
     legend_ch_match->AddEntry(fHistCh_read_pi[pix],"Histo Ambiguity distribution pion ","l");
     legend_ch_match->AddEntry(lin_ch_p_v,"p ","l");
     legend_ch_match->AddEntry(lin_ch_pi_v,"pi ","l");
     legend_ch_match->Draw();
     }
     }
     */
    
    
    
    
    prt_canvasSave(2,0);
    prt_canvasDel("*");
    
    
    
}


/////////////////
// Histo Style //
/////////////////
void HistoStyle(TH1F *x_histo, TH1F *y_histo, TH1F *z_histo) {
    x_histo->SetLineColor(kRed);
    x_histo->SetLineStyle(1);
    x_histo->GetXaxis()->SetTitle("number of photon per track [#]");
    x_histo->GetYaxis()->SetTitle("entries [#]");
    x_histo->GetXaxis()->SetTitleSize(0.05);
    x_histo->GetYaxis()->SetTitleSize(0.05);
    x_histo->GetXaxis()->SetTitleOffset(0.9);
    x_histo->GetYaxis()->SetTitleOffset(1.0);
    y_histo->SetLineColor(kBlue);
    y_histo->SetLineStyle(1);
    y_histo->GetXaxis()->SetTitle("number of photon per track [#]");
    y_histo->GetYaxis()->SetTitle("entries [#]");
    y_histo->GetXaxis()->SetTitleSize(0.05);
    y_histo->GetYaxis()->SetTitleSize(0.05);
    y_histo->GetXaxis()->SetTitleOffset(0.9);
    y_histo->GetYaxis()->SetTitleOffset(1.0);
    z_histo->SetLineColor(kBlack);
    z_histo->SetLineStyle(1);
    z_histo->GetXaxis()->SetTitle("number of photon per track [#]");
    z_histo->GetYaxis()->SetTitle("entries [#]");
    z_histo->GetXaxis()->SetTitleSize(0.05);
    z_histo->GetYaxis()->SetTitleSize(0.05);
    z_histo->GetXaxis()->SetTitleOffset(0.9);
    z_histo->GetYaxis()->SetTitleOffset(1.0);
}
/////////////////
// Graph Style //
/////////////////
void GraphStyle(TGraph *x_graph, TGraph *y_graph, TGraph *z_graph) {
    x_graph->Sort();
    x_graph->SetLineColor(kBlack);
    x_graph->SetMarkerColor(kBlack);
    x_graph->SetMarkerStyle(21);
    x_graph->SetMarkerSize(0.7);
    x_graph->SetName("separation");
    x_graph->GetYaxis()->SetRangeUser(0,20);
    x_graph->GetYaxis()->SetTitle("Separation [s.d]");
    x_graph->GetXaxis()->SetLabelSize(0.05);
    x_graph->GetXaxis()->SetTitleSize(0.06);
    x_graph->GetXaxis()->SetTitleOffset(0.84);
    
    y_graph->Sort();
    y_graph->SetLineColor(kBlack);
    y_graph->SetMarkerColor(kBlack);
    y_graph->SetMarkerStyle(21);
    y_graph->SetMarkerSize(0.7);
    y_graph->SetName("separation");
    y_graph->GetYaxis()->SetRangeUser(0,20);
    y_graph->GetYaxis()->SetTitle("Separation [s.d]");
    y_graph->GetXaxis()->SetLabelSize(0.05);
    y_graph->GetXaxis()->SetTitleSize(0.06);
    y_graph->GetXaxis()->SetTitleOffset(0.84);
    
    z_graph->Sort();
    z_graph->SetLineColor(kBlack);
    z_graph->SetMarkerColor(kBlack);
    z_graph->SetMarkerStyle(21);
    z_graph->SetMarkerSize(0.7);
    z_graph->SetName("separation");
    z_graph->GetYaxis()->SetRangeUser(0,20);
    z_graph->GetYaxis()->SetTitle("Separation [s.d]");
    z_graph->GetXaxis()->SetLabelSize(0.05);
    z_graph->GetXaxis()->SetTitleSize(0.06);
    z_graph->GetXaxis()->SetTitleOffset(0.84);
}

//////////////////////////
// check file existance //
//////////////////////////
bool exists_test (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}
























