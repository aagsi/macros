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
#define PI 3.14159265

TCanvas *c1 = new TCanvas("c1","c1",300,200);

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
Double_t  prtangle;

////////////////////
// proto types//////
////////////////////
// file existance
bool exists_test (const std::string& name);
// histo style
void HistoStyle(TH1F *p_cherenkov_sim, TH1F *p_cherenkov_data, TH1F *p_cherenkov_data_sub, TH1F *p_cherenkov_mc_same_path, TH1F *p_cherenkov_bg_sim, TH1F *nph_sim, TH1F *p_nph_sim,TH1F *p_nph_good_sim, TH1F *p_nph_true_sim, TH1F *p_diff_time_sim, TH1F*p_diff_time_data, TH1F*p_diff_time_data_sub, TH1F*p_diff_time_mctruth, TH1F* p_diff_time_bg_sim, TH1F *p_cherenkov_data_copy, TH1F *p_photon_time_sim, TH1F *p_photon_time_data, TH1F *p_photon_time_sim_calc, TH1F *p_photon_time_data_calc, TH1F *nph_data, TH1F *p_nph_data, TH1F *p_nph_good_data, TH1F *dac_hits_data, TH1F *dac_hits_sys_cus_data);
void HistoStyleMatch(TH1F *p_cherenkov_sim, TH1F *p_cherenkov_data);
void HistoStyle_3colors(TH1F *histo1, TH1F *histo2, TH1F *histo3);

Double_t*  YieldGausFit(TH1F *photonYieldHistogram);
Double_t  YieldGausFit_double(TH1F *photonYieldHistogram) ;

// diff norm
void DiffNorm(TH1F *p_diff_time_sim, TH1F *p_diff_time_data);
// fit style
void FitStyle();
Double_t momentum=7.0;
Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
Double_t fAngleP = acos(sqrt(momentum*momentum+ mass[4]*mass[4])/momentum/1.4738)-0.00;
Double_t fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00; //-0.0014 for 160 25deg


////////////////////
// function   //////
////////////////////
void proton_plots(
                  
                  Bool_t bool_part1=false,
                  Bool_t bool_part1_1=false,
                  Bool_t bool_part1_2=false,
                  Bool_t bool_part1_3=false,
                  Bool_t bool_part1_4=false,
                  
                  Bool_t bool_part2=false,
                  Bool_t bool_part2_1=false,
                  Bool_t bool_part2_2=false,
                  Bool_t bool_part2_3=false,
                  Bool_t bool_part2_4=false,
                  Bool_t bool_part2_5=false,
                  Bool_t bool_part2_6=false,
                  Bool_t bool_part2_7=false,
                  
                  Bool_t bool_part3=false,
                  Bool_t bool_part3_1=false,
                  Bool_t bool_part3_2=false,
                  Bool_t bool_part3_3=false)
{

    vector<Double_t> chAngleCut(14);
    vector<Double_t> recoAngle(14);
    vector<Double_t> timeCut(14);
    vector<Int_t> prtangle_vector(14);
    const Int_t n = 14;
    Double_t x[n], y_spr_data_sub[n], y_spr_sim_true[n], y_spr_sim_org[n], y_diff_true[n], y_diff_data[n], y_mean_diff_data[n], y_mean_diff_sim[n], y_mean_diff_true[n], y_diff_sim[n];
    Double_t y_cangle_data_sub[n], y_cangle_sim_org[n];
    Double_t y_cangle_sim_org_true[n];
    Double_t y_spr_data_cuts[n], y_spr_sim_cuts[n], y_spr_data_sub_cuts[n], y_spr_sim_true_cuts[n], y_yield_nph_sim[n], y_yield_p_nph_sim[n], y_yield_p_nph_good_sim[n],  y_yield_p_nph_true_sim[n], y_yield_nph_data[n], y_yield_p_nph_data[n], y_yield_p_nph_good_data[n], y_yield_dac_hits_sys_cus_data[n], y_cangle_data_correction[n], y_spr_data_org[n], y_cangle_data_org[n], y_spr_data_correction[n], y_spr_data_correction_error[n], y_spr_sim_correction[n], y_cangle_sim_org_correction[n];
    
    Double_t y_yield_nph_sim_error[n], y_yield_p_nph_sim_error[n], y_yield_p_nph_good_sim_error[n], y_yield_p_nph_true_sim_error[n], y_yield_nph_data_error[n], y_yield_p_nph_data_error[n], y_yield_p_nph_good_data_error[n];
    
    int counter =0 ;
    TLegend* legend = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
    legend->SetHeader("cherenkov angle","C"); // option "C" allows to center the header
    prt_savepath="proton";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    //Int_t nf = 30;
    //TH2F * photon_yield_map =  new TH2F("photon_yield_map",";Time Cut [ns];cherenkov angle cut[mrad]", nf, 0.1, 6, nf, 0.001,0.06);
    TFile *ffile_sim, *ffile_data;
    TH1F *tof_pid, *p_cherenkov_sim, *p_cherenkov_data, *p_cherenkov_data_copy, *p_cherenkov_sim_copy, *p_cherenkov_data_sub, *p_cherenkov_mc_same_path, *p_cherenkov_bg_sim, *p_cherenkov_data_correction, *p_cherenkov_sim_correction, *nph_sim, *p_nph_sim, *p_nph_good_sim, *p_nph_true_sim;
    TH1F *p_diff_time_sim, *p_diff_time_data, *p_diff_time_data_sub, *p_diff_time_mctruth, * p_diff_time_bg_sim, *p_photon_time_sim, *p_photon_time_data,*p_photon_time_sim_calc, *p_photon_time_data_calc;
    TH1F *nph_data, *p_nph_data, *p_nph_good_data, *dac_hits_data, *dac_hits_sys_cus_data;
    
    TH1F *hist_ambiguity_sim, *histo_photon_ambiguity_wo_sim, *histo_photon_ambiguity_wt_sim, *histo_photon_ambiguity_wtc_sim;
    TH1F *hist_ambiguity_data, *histo_photon_ambiguity_wo_data, *histo_photon_ambiguity_wt_data, *histo_photon_ambiguity_wtc_data;
    

    THStack *hs, *hs2, *hs3, *hs4, *hs5, *hs6, *hs7, *hs8, *hs9, *hs10;
    std::cout<<"############"<< " no problem 0 " <<std::endl;
    TGraph *calc_mom = new TGraph();
    TGraph *calc_e_mom = new TGraph();
    TGraph *calc_tof1tof2_distance = new TGraph();
    TGraph *calc_e_tof1tof2_distance = new TGraph();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    
    TGraph * graph_spr_sim,* graph_spr_data_sub,* graph_spr_sim_true,* graph_spr_data_correction,* graph_spr_sim_correction,* graph_spr_data_org,* graph_spr_data_cuts,* graph_spr_sim_cuts,* graph_spr_data_sub_cuts,* graph_spr_sim_true_cuts,* graph_diff_true,* graph_diff_data,* graph_diff_sim,* graph_diff_data_mean,* graph_diff_sim_mean,* graph_diff_true_mean,* graph_diff_true_sim,* graph_cangle_data_sub,* graph_cangle_sim,* graph_cangle_sim_true,* graph_cangle_data_correction,* graph_cangle_sim_correction,* graph_cangle_data_org,* graph_yield_DIRC_wo_cuts_sim,* graph_yield_DIRC_wt_cuts_sim,* graph_yield_DIRC_wtc_cuts_sim,* graph_yield_DIRC_true_sim,* graph_yield_DIRC_wo_cuts_data,* graph_yield_DIRC_wt_cuts_data,* graph_yield_DIRC_wtc_cuts_data;
    
    
    for (int i=20; i<=150; i+=10) {
        prtangle= i;
        TString nid = Form("_%2.0d", i);
        cout<< "enter the if condition"<<endl;
        // proton
        TString cherenkov_data_path = Form("/Users/ahmed/perforamnce/spr_data_sim/spr_wtb_%d_sph_p_data_spr.root", i);
        TString cherenkov_sim_path = Form("/Users/ahmed/perforamnce/spr_data_sim/spr_wt_%d_sph_p_sim_spr.root", i);
        // pi
        //TString cherenkov_data_path = Form("/u/aali/work/%d_sph_pi_data_spr.root", i);
        //TString cherenkov_sim_path = Form("/u/aali/work/reco_pi_bar_3lsph_grease_theta_%d_sim_spr.root", i);
        
        cout<<"cherenkov_sim_path= " <<cherenkov_sim_path<<endl;
        cout<<"cherenkov_data_path= " <<cherenkov_data_path<<endl;
        string path_sim = (string)cherenkov_sim_path;
        string path_data = (string)cherenkov_data_path;
        cout<<"exists_test(path_sim)" <<exists_test(path_sim)<<endl;
        cout<<"exists_test(path_data)" <<exists_test(path_data)<<endl;
        if (!exists_test(path_sim)) continue;
        if (!exists_test(path_data)) continue;
        cout<<"path_sim= " <<path_sim<<endl;
        cout<<"path_data= " <<path_data<<endl;
        ///////////////////////////////
        // READ files and histograms ///
        ///////////////////////////////
        gROOT->ForceStyle(kTRUE);
        ffile_sim  = new TFile(cherenkov_sim_path, "READ");
        ffile_data  = new TFile(cherenkov_data_path, "READ");
        
        //////////////////////////
        // cherenkov histogram //
        //////////////////////////
        p_cherenkov_sim=(TH1F*)ffile_sim->Get("fHist"); // changed from fHist
        //p_cherenkov_sim->SetStats(0);
        p_cherenkov_data=(TH1F*)ffile_data->Get("fHist"); // changed from fHist
        p_cherenkov_data->SetStats(0);//kFALSE
        p_cherenkov_data_copy=(TH1F*)ffile_data->Get("fHist_copy");
        p_cherenkov_sim_copy=(TH1F*)ffile_sim->Get("fHist_copy");
        p_cherenkov_mc_same_path=(TH1F*)ffile_sim->Get("fHist_same_path");
        p_cherenkov_bg_sim=(TH1F*)ffile_sim->Get("fHist_bg");
        p_cherenkov_data_correction=(TH1F*)ffile_data->Get("fHist_correction");
        p_cherenkov_sim_correction=(TH1F*)ffile_sim->Get("fHist_correction");
        
        ////////////////////////////
        // photon yield histogram //
        ////////////////////////////
        nph_sim=(TH1F*)ffile_sim->Get("fnHits");
        p_nph_sim=(TH1F*)ffile_sim->Get("fnHits_p");
        p_nph_good_sim=(TH1F*)ffile_sim->Get("fnHits_p_good");
        p_nph_true_sim=(TH1F*)ffile_sim->Get("fnHits_true_sim");
        
        nph_data=(TH1F*)ffile_data->Get("fnHits");
        p_nph_data=(TH1F*)ffile_data->Get("fnHits_p");
        p_nph_good_data=(TH1F*)ffile_data->Get("fnHits_p_good");
        //dac_hits_data=(TH1F*)ffile_data->Get("nHist_dac");
        //dac_hits_sys_cus_data=(TH1F*)ffile_data->Get("nHist_dac_syscut_p");
        
        /////////////////////////
        // time diff histogram //
        /////////////////////////
        p_diff_time_sim=(TH1F*)ffile_sim->Get("timediff");
        p_diff_time_mctruth=(TH1F*)ffile_sim->Get("timediffi");
        p_diff_time_bg_sim=(TH1F*)ffile_sim->Get("timediffi_bg");
        p_diff_time_data=(TH1F*)ffile_data->Get("timediff");
        //p_diff_time_data_sub = (TH1F *)p_diff_time_data->Clone();
        
        ///////////////////////////
        // photon time histogram //
        ///////////////////////////
        p_photon_time_sim=(TH1F*)ffile_sim->Get("time1");
        p_photon_time_data=(TH1F*)ffile_data->Get("time1");
        p_photon_time_sim_calc=(TH1F*)ffile_sim->Get("time2");
        p_photon_time_data_calc=(TH1F*)ffile_data->Get("time2");
        p_photon_time_data->Scale(p_photon_time_sim->GetMaximum() /p_photon_time_data->GetMaximum());
        p_photon_time_data_calc->Scale(p_photon_time_sim_calc->GetMaximum() /p_photon_time_data_calc->GetMaximum());
        
        /////////////////////////
        // ambiguity histogram //
        /////////////////////////
        hist_ambiguity_sim =(TH1F*)ffile_sim->Get("hist_ambiguity");
        histo_photon_ambiguity_wo_sim=(TH1F*)ffile_sim->Get("histo_photon_ambiguity_wo");
        histo_photon_ambiguity_wt_sim=(TH1F*)ffile_sim->Get("histo_photon_ambiguity_wt");
        histo_photon_ambiguity_wtc_sim=(TH1F*)ffile_sim->Get("histo_photon_ambiguity_wtc");
        
        hist_ambiguity_data =(TH1F*)ffile_data->Get("hist_ambiguity");
        histo_photon_ambiguity_wo_data=(TH1F*)ffile_data->Get("histo_photon_ambiguity_wo");
        histo_photon_ambiguity_wt_data=(TH1F*)ffile_data->Get("histo_photon_ambiguity_wt");
        histo_photon_ambiguity_wtc_data=(TH1F*)ffile_data->Get("histo_photon_ambiguity_wtc");
        
        HistoStyle_3colors(histo_photon_ambiguity_wo_sim, histo_photon_ambiguity_wt_sim, histo_photon_ambiguity_wtc_sim);
        HistoStyle_3colors(histo_photon_ambiguity_wo_data, histo_photon_ambiguity_wt_data, histo_photon_ambiguity_wtc_data);
        
        /////////////
        // TOF PID //
        /////////////
        tof_pid=(TH1F*)ffile_data->Get("hdelta_tof2tof1");
        
        ///////////////////
        ///// part I //////
        ///////////////////
        if(bool_part1) {
            if (bool_part1_1) {
                gStyle->SetOptFit(0);
                gStyle->SetOptStat(0);
                DiffNorm(p_cherenkov_sim, p_cherenkov_data);
                HistoStyleMatch(p_cherenkov_sim, p_cherenkov_data);
                TLegend *legend_ch_match= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                legend_ch_match->SetHeader("Ambiguity distribution (proton data)","C");
                prt_canvasAdd("r_ch_match"+nid,800,400);
                p_cherenkov_sim->SetTitle(Form("Polar angle %3.1f", prtangle));
                legend_ch_match->AddEntry(p_cherenkov_sim,"Ambiguity distribution sim","f");
                legend_ch_match->AddEntry(p_cherenkov_data," Ambiguity distribution data ","f");
                p_cherenkov_sim->Draw("hist");
                p_cherenkov_data->Draw("samehist");
                legend_ch_match->Draw();
                prt_canvasGet("r_ch_match"+nid)->Update();
                TLine *lin_ch_p_v = new TLine(0,0,0,1000);
                lin_ch_p_v->SetX1(fAngleP);
                lin_ch_p_v->SetX2(fAngleP);
                lin_ch_p_v->SetY1(gPad->GetUymin());
                lin_ch_p_v->SetY2(gPad->GetUymax());
                lin_ch_p_v->SetLineColor(kRed);
                TLine *lin_ch_pi_v = new TLine(0,0,0,1000);
                lin_ch_pi_v->SetX1(fAnglePi);
                lin_ch_pi_v->SetX2(fAnglePi);
                lin_ch_pi_v->SetY1(gPad->GetUymin());
                lin_ch_pi_v->SetY2(gPad->GetUymax());
                lin_ch_pi_v->SetLineColor(kBlue);
                lin_ch_p_v->Draw();
                lin_ch_pi_v->Draw();
                prt_canvasGet("r_ch_match"+nid)->Update();
            }
            
            if(bool_part1_2) {
                prt_canvasAdd("r_correction"+nid,800,400);
                p_cherenkov_data_correction->SetTitle(Form("Polar angle %3.1f (proton data)", prtangle));
                HistoStyleMatch(p_cherenkov_data,p_cherenkov_data_correction);
                p_cherenkov_data_correction->Draw();
                p_cherenkov_data->Draw("same");
                TLegend *legend_correction= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                
                TF1 *fFit_p_cherenkov_data_correction = new TF1("fFit_p_cherenkov_data_correction","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
                fFit_p_cherenkov_data_correction->SetLineColor(1);
                Double_t cangle_cor =  p_cherenkov_data_correction->GetXaxis()->GetBinCenter(p_cherenkov_sim->GetMaximumBin());
                if(cangle_cor>0.85) cangle_cor=0.82;
                fFit_p_cherenkov_data_correction->SetParameters(100,cangle_cor,0.010);
                fFit_p_cherenkov_data_correction->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
                fFit_p_cherenkov_data_correction->SetParLimits(0,0.1,1E6);
                fFit_p_cherenkov_data_correction->SetParLimits(1,cangle_cor-0.04,cangle_cor+0.04);
                fFit_p_cherenkov_data_correction->SetParLimits(2,0.005,0.014000);
                p_cherenkov_data_correction->Fit("fFit_p_cherenkov_data_correction","M","",cangle_cor-0.06,cangle_cor+0.06);
                
                legend_correction->SetHeader("Ambiguity distribution (proton data)","C");
                legend_correction->AddEntry(p_cherenkov_data,"Ambiguity distribution","f");
                legend_correction->AddEntry(p_cherenkov_data_correction," Ambiguity distribution corrected ","f");
                legend_correction->Draw();
                prt_canvasGet("r_correction"+nid)->Update();
                TLine *line_ch_p_v = new TLine(0,0,0,1000);
                line_ch_p_v->SetX1(fAngleP);
                line_ch_p_v->SetX2(fAngleP);
                line_ch_p_v->SetY1(gPad->GetUymin());
                line_ch_p_v->SetY2(gPad->GetUymax());
                line_ch_p_v->SetLineColor(kRed);
                
                TLine *line_ch_pi_v = new TLine(0,0,0,1000);
                line_ch_pi_v->SetX1(fAnglePi);
                line_ch_pi_v->SetX2(fAnglePi);
                line_ch_pi_v->SetY1(gPad->GetUymin());
                line_ch_pi_v->SetY2(gPad->GetUymax());
                line_ch_pi_v->SetLineColor(kBlue);
                line_ch_p_v->Draw();
                line_ch_pi_v->Draw();
                prt_canvasGet("r_correction"+nid)->Update();
            }
            
            if(bool_part1_3) {
                DiffNorm(p_cherenkov_sim_correction, p_cherenkov_data_correction);
                prt_canvasAdd("r_corrected_match"+nid,800,400);
                p_cherenkov_data_correction->SetTitle(Form("Polar angle %3.1f (proton data)", prtangle));
                HistoStyleMatch(p_cherenkov_sim_correction,p_cherenkov_data_correction);
                p_cherenkov_data_correction->Draw();
                p_cherenkov_sim_correction->Draw("same");
                TLegend *legend_correction= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                
                legend_correction->SetHeader("Ambiguity distribution corrected match )","C");
                legend_correction->AddEntry(p_cherenkov_data_correction,"Ambiguity distribution (proton data corrected )","p");
                legend_correction->AddEntry(p_cherenkov_sim_correction,"Ambiguity distribution (proton sim corrected )","f");
                legend_correction->Draw();
                prt_canvasGet("r_corrected_match"+nid)->Update();
                TLine *line_ch_p_v = new TLine(0,0,0,1000);
                line_ch_p_v->SetX1(fAngleP);
                line_ch_p_v->SetX2(fAngleP);
                line_ch_p_v->SetY1(gPad->GetUymin());
                line_ch_p_v->SetY2(gPad->GetUymax());
                line_ch_p_v->SetLineColor(kRed);
                
                TLine *line_ch_pi_v = new TLine(0,0,0,1000);
                line_ch_pi_v->SetX1(fAnglePi);
                line_ch_pi_v->SetX2(fAnglePi);
                line_ch_pi_v->SetY1(gPad->GetUymin());
                line_ch_pi_v->SetY2(gPad->GetUymax());
                line_ch_pi_v->SetLineColor(kBlue);
                line_ch_p_v->Draw();
                line_ch_pi_v->Draw();
                prt_canvasGet("r_corrected_match"+nid)->Update();
            }
        }
        
        ///////////////////
        ///// part II /////
        ///////////////////
        if(true) {
            gROOT->SetBatch(1);
            TSpectrum *fSpect= new TSpectrum(10);
            //Int_t nfound = fSpect->Search(p_cherenkov_mc_same_path,1,"",0.9); //0.6
            //Double_t cangle_MC_true = fSpect->GetPositionX()[0];
            Double_t cangle_MC_true =  p_cherenkov_mc_same_path->GetXaxis()->GetBinCenter(p_cherenkov_mc_same_path->GetMaximumBin());
            /////////////////////////////
            // cherenkov normalization /
            /////////////////////////////
            TAxis *axis = p_cherenkov_sim_correction->GetXaxis();
            double xmin = 0.6;
            double xmax = 0.74;
            if (i==90)xmin = 0.9;
            if (i==90)xmax = 1.0;
            int bmin = axis->FindBin(xmin);
            int bmax = axis->FindBin(xmax);
            double integral = p_cherenkov_sim_correction->Integral(bmin,bmax);
            TAxis *axis_data = p_cherenkov_data_correction->GetXaxis();
            double xmin_data = 0.6;
            double xmax_data = 0.74;
            if (i==90)xmin_data = 0.9;
            if (i==90)xmax_data = 1.0;
            int bmin_data = axis_data->FindBin(xmin_data);
            int bmax_data = axis_data->FindBin(xmax_data);
            double integral_data = p_cherenkov_data_correction->Integral(bmin_data,bmax_data);
            Double_t norm= integral/integral_data ;
            //std::cout<<"############  norm= "<< norm <<std::endl;
            p_cherenkov_data_copy->Scale(norm);
            p_cherenkov_data_sub = (TH1F *)p_cherenkov_data_copy->Clone();
            p_cherenkov_data_sub->Add(p_cherenkov_bg_sim,-1);
            std::cout<<"############"<< " no problem 3 " <<std::endl;
            ////////////////////
            // fit ch data sub//
            ////////////////////
            TF1 *fFit_p_cherenkov_data_sub = new TF1("fFit_p_cherenkov_data_sub","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
            fFit_p_cherenkov_data_sub->SetParameters(100,cangle_MC_true,0.010);
            fFit_p_cherenkov_data_sub->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit_p_cherenkov_data_sub->SetParLimits(0,0.1,1E6);
            fFit_p_cherenkov_data_sub->SetParLimits(1,cangle_MC_true-0.04,cangle_MC_true+0.04);
            fFit_p_cherenkov_data_sub->SetParLimits(2,0.005,0.014000);
            p_cherenkov_data_sub->Fit("fFit_p_cherenkov_data_sub","M","",cangle_MC_true-0.06,cangle_MC_true+0.06);
            Double_t cangle_data_sub = fFit_p_cherenkov_data_sub->GetParameter(1);
            Double_t spr_data_sub = fFit_p_cherenkov_data_sub->GetParameter(2);
            Double_t cangle_data_sub_minus_5_sgma = cangle_data_sub-5*spr_data_sub;
            Double_t cangle_data_sub_plus_5_sgma = cangle_data_sub+5*spr_data_sub;
            Double_t cangle_data_sub_minus_3_sgma = cangle_data_sub-3*spr_data_sub;
            Double_t cangle_data_sub_plus_3_sgma = cangle_data_sub+3*spr_data_sub;
            Double_t r_min_data_sub = cangle_data_sub-8*spr_data_sub;
            Double_t r_max_data_sub = cangle_data_sub+8*spr_data_sub;
            //////////////////////////////////////////
            // fit ch data copy after normalization //
            //////////////////////////////////////////
            TF1 *fFit_p_cherenkov_data_copy = new TF1("fFit_p_cherenkov_data_copy","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
            fFit_p_cherenkov_data_copy->SetParameters(100,cangle_MC_true,0.010);
            fFit_p_cherenkov_data_copy->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit_p_cherenkov_data_copy->SetParLimits(0,0.1,1E6);
            fFit_p_cherenkov_data_copy->SetParLimits(1,cangle_MC_true-0.04,cangle_MC_true+0.04);
            fFit_p_cherenkov_data_copy->SetParLimits(2,0.005,0.014000);
            p_cherenkov_data_copy->Fit("fFit_p_cherenkov_data_copy","M","",cangle_MC_true-0.06,cangle_MC_true+0.06);
            p_cherenkov_data_copy->SetTitle(Form("Polar angle %3.1f", prtangle));
            Double_t cangle_data_copy = fFit_p_cherenkov_data_copy->GetParameter(1);
            Double_t spr_data_copy = fFit_p_cherenkov_data_copy->GetParameter(2);
            ////////////////////////////////////
            // fit Cherenkove data correction //
            ////////////////////////////////////
            prt_canvasAdd("r_ch_fit_data"+nid,800,400);
            TF1 *fFit_p_cherenkov_data_correction = new TF1("fFit_p_cherenkov_data_correction","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
            fFit_p_cherenkov_data_correction->SetParameters(100,cangle_MC_true,0.010);
            fFit_p_cherenkov_data_correction->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit_p_cherenkov_data_correction->SetParLimits(0,0.1,1E6);
            fFit_p_cherenkov_data_correction->SetParLimits(1,cangle_MC_true-0.04,cangle_MC_true+0.04);
            fFit_p_cherenkov_data_correction->SetParLimits(2,0.005,0.014000);
            p_cherenkov_data_correction->Fit("fFit_p_cherenkov_data_correction","M","",cangle_MC_true-0.06,cangle_MC_true+0.06);
            p_cherenkov_data_correction->SetTitle(Form("Polar angle %3.1f", prtangle));
            p_cherenkov_data_correction->Draw();
            Double_t cangle_data_correction = fFit_p_cherenkov_data_correction->GetParameter(1);
            Double_t spr_data_correction = fFit_p_cherenkov_data_correction->GetParameter(2);
            ///////////////////////////////////
            // fit Cherenkove sim correction //
            ///////////////////////////////////
            prt_canvasAdd("r_ch_fit_sim"+nid,800,400);
            TF1 *fFit_p_cherenkov_sim_correction = new TF1("fFit_p_cherenkov_sim_correction","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
            fFit_p_cherenkov_sim_correction->SetParameters(100,cangle_MC_true,0.010);
            fFit_p_cherenkov_sim_correction->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit_p_cherenkov_sim_correction->SetParLimits(0,0.1,1E6);
            fFit_p_cherenkov_sim_correction->SetParLimits(1,cangle_MC_true-0.04,cangle_MC_true+0.04);
            fFit_p_cherenkov_sim_correction->SetParLimits(2,0.005,0.014000);
            p_cherenkov_sim_correction->Fit("fFit_p_cherenkov_sim_correction","M","",cangle_MC_true-0.06,cangle_MC_true+0.06);
            p_cherenkov_sim_correction->SetTitle(Form("Polar angle %3.1f", prtangle));
            p_cherenkov_sim_correction->Draw();
            Double_t cangle_sim_correction = fFit_p_cherenkov_sim_correction->GetParameter(1);
            Double_t spr_sim_correction = fFit_p_cherenkov_sim_correction->GetParameter(2);
            //////////////////////////////////////
            // fit Cherenkove sim WO correction //
            //////////////////////////////////////
            TF1 *fFit_p_cherenkov_sim = new TF1("fFit_p_cherenkov_sim","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
            fFit_p_cherenkov_sim->SetParameters(100,cangle_MC_true,0.010);
            fFit_p_cherenkov_sim->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit_p_cherenkov_sim->SetParLimits(0,0.1,1E6);
            fFit_p_cherenkov_sim->SetParLimits(1,cangle_MC_true-0.04,cangle_MC_true+0.04);
            fFit_p_cherenkov_sim->SetParLimits(2,0.005,0.014000);
            p_cherenkov_sim->Fit("fFit_p_cherenkov_sim","0","",cangle_MC_true-0.06,cangle_MC_true+0.06);
            Double_t chi = fFit_p_cherenkov_sim->GetChisquare()/fFit_p_cherenkov_sim->GetNDF();
            Double_t cangle_sim_org = fFit_p_cherenkov_sim->GetParameter(1);
            Double_t spr_sim_org = fFit_p_cherenkov_sim->GetParameter(2);
            //////////////////////////////////////
            // fit Cherenkove data WO correction //
            //////////////////////////////////////
            TF1 *fFit_p_cherenkov_data = new TF1("fFit_p_cherenkov_data","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
            fFit_p_cherenkov_data->SetParameters(100,cangle_MC_true,0.010);
            fFit_p_cherenkov_data->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit_p_cherenkov_data->SetParLimits(0,0.1,1E6);
            fFit_p_cherenkov_data->SetParLimits(1,cangle_MC_true-0.04,cangle_MC_true+0.04);
            fFit_p_cherenkov_data->SetParLimits(2,0.005,0.014000);
            p_cherenkov_data->Fit("fFit_p_cherenkov_data","0","",cangle_MC_true-0.06,cangle_MC_true+0.06);
            Double_t cangle_data_org = fFit_p_cherenkov_data->GetParameter(1);
            Double_t spr_data_org = fFit_p_cherenkov_data->GetParameter(2);
            ///////////////////
            // fit true MC  ///
            ///////////////////
            TF1 *fFit_p_cherenkov_mc_same_path = new TF1("fFit_p_cherenkov_mc_same_path","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
            fFit_p_cherenkov_mc_same_path->SetParameters(100,cangle_MC_true,0.010);
            fFit_p_cherenkov_mc_same_path->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit_p_cherenkov_mc_same_path->SetParLimits(0,0.1,1E6);
            fFit_p_cherenkov_mc_same_path->SetParLimits(1,cangle_MC_true-0.04,cangle_MC_true+0.04);
            fFit_p_cherenkov_mc_same_path->SetParLimits(2,0.005,0.014000);
            p_cherenkov_mc_same_path->Fit("fFit_p_cherenkov_mc_same_path","M","",cangle_MC_true-0.06,cangle_MC_true+0.06);
            Double_t spr_sim_true = fFit_p_cherenkov_mc_same_path->GetParameter(2);
            Double_t cangle_sim_true = fFit_p_cherenkov_mc_same_path->GetParameter(1);
            Double_t cangle_minus_5_sgma = cangle_sim_true-5*spr_sim_true;
            Double_t cangle_plus_5_sgma = cangle_sim_true+5*spr_sim_true;
            Double_t cangle_minus_3_sgma = cangle_sim_true-3*spr_sim_true;
            Double_t cangle_plus_3_sgma = cangle_sim_true+3*spr_sim_true;
            
            //////////////////////
            // function call   ///
            //////////////////////

            // time diff normalization
            DiffNorm(p_diff_time_sim, p_diff_time_data);
            p_diff_time_data_sub = (TH1F *)p_diff_time_data->Clone();
            p_diff_time_data_sub->Add(p_diff_time_bg_sim,-1);
            
            // histogram style
            HistoStyle(p_cherenkov_sim, p_cherenkov_data, p_cherenkov_data_sub, p_cherenkov_mc_same_path, p_cherenkov_bg_sim, nph_sim, p_nph_sim, p_nph_good_sim, p_nph_true_sim, p_diff_time_sim, p_diff_time_data, p_diff_time_data_sub, p_diff_time_mctruth, p_diff_time_bg_sim, p_cherenkov_data_copy, p_photon_time_sim, p_photon_time_data, p_photon_time_sim_calc, p_photon_time_data_calc, nph_data, p_nph_data, p_nph_good_data, dac_hits_data, dac_hits_sys_cus_data);

            ////////////////////////////
            // fit true time diff path //
            /////////////////////////////
            Double_t diff_peak_mctruth =  p_diff_time_mctruth->GetXaxis()->GetBinCenter(p_diff_time_mctruth->GetMaximumBin());
            TF1 *fFit_p_diff_time_mctruth = new TF1("fFit_p_diff_time_mctruth","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.35,0.9);
            fFit_p_diff_time_mctruth->SetParameters(100,diff_peak_mctruth,0.010);
            fFit_p_diff_time_mctruth->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit_p_diff_time_mctruth->SetParLimits(0,0.1,1E6);
            fFit_p_diff_time_mctruth->SetParLimits(1,-5,5);
            fFit_p_diff_time_mctruth->SetParLimits(2,-5,5);
            p_diff_time_mctruth->Fit("fFit_p_diff_time_mctruth","M","",-6,6);
            Double_t mean_p_diff_time_mctruth= fFit_p_diff_time_mctruth->GetParameter(1);
            Double_t sigma_p_diff_time_mctruth= fFit_p_diff_time_mctruth->GetParameter(2);

            
            /////////////////////////////
            // fit data time diff path //
            /////////////////////////////
            Double_t diff_peak_data =  p_diff_time_data->GetXaxis()->GetBinCenter(p_diff_time_data->GetMaximumBin());
            TF1 *fFit_p_diff_time_data = new TF1("fFit_p_diff_time_data","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.35,0.9);
            fFit_p_diff_time_data->SetParameters(100,diff_peak_data,0.010);
            fFit_p_diff_time_data->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit_p_diff_time_data->SetParLimits(0,0.1,1E6);
            fFit_p_diff_time_data->SetParLimits(1,diff_peak_data-3*sigma_p_diff_time_mctruth,diff_peak_data+3*sigma_p_diff_time_mctruth);
            fFit_p_diff_time_data->SetParLimits(2,0,6);

            if (i == 120) {
                p_diff_time_data->Fit("fFit_p_diff_time_data","0","",diff_peak_data-1.5*sigma_p_diff_time_mctruth,diff_peak_data+1.0*sigma_p_diff_time_mctruth);
            }
            else if (i == 110) {
                p_diff_time_data->Fit("fFit_p_diff_time_data","0","",diff_peak_data-1.5*sigma_p_diff_time_mctruth,diff_peak_data+0.7*sigma_p_diff_time_mctruth);
            }
            else if (i == 100) {
                p_diff_time_data->Fit("fFit_p_diff_time_data","0","",diff_peak_data-2.0*sigma_p_diff_time_mctruth,diff_peak_data+2.0*sigma_p_diff_time_mctruth);
            }
            else if (i == 90) {
                p_diff_time_data->Fit("fFit_p_diff_time_data","0","",diff_peak_data-1.0*sigma_p_diff_time_mctruth,diff_peak_data+1.5*sigma_p_diff_time_mctruth);
            }
            else if (i == 50) {
                p_diff_time_data->Fit("fFit_p_diff_time_data","0","",diff_peak_data-1.5*sigma_p_diff_time_mctruth,diff_peak_data+1.5*sigma_p_diff_time_mctruth);
            }
            else if (i == 150) {
                p_diff_time_data->Fit("fFit_p_diff_time_data","0","",diff_peak_data-1.5*sigma_p_diff_time_mctruth,diff_peak_data+1.0*sigma_p_diff_time_mctruth);
            }
            else if (i == 130) {
                p_diff_time_data->Fit("fFit_p_diff_time_data","0","",diff_peak_data-1.0*sigma_p_diff_time_mctruth,diff_peak_data+1.0*sigma_p_diff_time_mctruth);
            }
            else if (i == 60) {
                p_diff_time_data->Fit("fFit_p_diff_time_data","0","",diff_peak_data-1.5*sigma_p_diff_time_mctruth,diff_peak_data+0.8*sigma_p_diff_time_mctruth);
            }
            else {
                p_diff_time_data->Fit("fFit_p_diff_time_data","0","",diff_peak_data-0.8*sigma_p_diff_time_mctruth,diff_peak_data+0.8*sigma_p_diff_time_mctruth);
                
            }
            Double_t mean_p_diff_time_data= fFit_p_diff_time_data->GetParameter(1);
            Double_t sigma_p_diff_time_data = fFit_p_diff_time_data->GetParameter(2);
            
            /////////////////////////////
            // fit sim time diff path //
            /////////////////////////////
            Double_t diff_peak_sim =  p_diff_time_sim->GetXaxis()->GetBinCenter(p_diff_time_sim->GetMaximumBin());
            TF1 *fFit_p_diff_time_sim = new TF1("fFit_p_diff_time_sim","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.35,0.9);
            fFit_p_diff_time_sim->SetParameters(100,diff_peak_sim,0.010);
            fFit_p_diff_time_sim->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            fFit_p_diff_time_sim->SetParLimits(0,0.1,1E6);
            fFit_p_diff_time_sim->SetParLimits(1,diff_peak_data-3*sigma_p_diff_time_mctruth,diff_peak_data+3*sigma_p_diff_time_mctruth);
            fFit_p_diff_time_sim->SetParLimits(2,0,6);
            p_diff_time_sim->Fit("fFit_p_diff_time_sim","0","",diff_peak_sim-1*sigma_p_diff_time_mctruth,diff_peak_data+1*sigma_p_diff_time_mctruth);
            Double_t mean_p_diff_time_sim= fFit_p_diff_time_sim->GetParameter(1);
            Double_t sigma_p_diff_time_sim = fFit_p_diff_time_sim->GetParameter(2);
            
            //////////////////////////////////
            //   mom distance calculations ///
            //////////////////////////////////
            
            Double_t calc_p_tof2tof1_plot1_co(0), distance_tof2tof1_plot1_co(0);
            TSpectrum *s_tof2tof1_plot1 = new TSpectrum(5);
            Int_t nfound_tof2tof1_plot1 = s_tof2tof1_plot1->Search(tof_pid,2,"",0.10);
            printf("Found %d candidate peaks to fit\n",nfound_tof2tof1_plot1);
            //Estimate background using TSpectrum::Background
            TH1 *hb_tof2tof1_plot1 = s_tof2tof1_plot1->Background(tof_pid,200,"same");
            //Loop on all found peaks. Eliminate peaks at the background level
            Double_t *xpeaks_tof2tof1_plot1 = s_tof2tof1_plot1->GetPositionX();
            for (Int_t p_tof2tof1_plot1=0; p_tof2tof1_plot1<nfound_tof2tof1_plot1; p_tof2tof1_plot1++) {
                Float_t xp_tof2tof1_plot1 = xpeaks_tof2tof1_plot1[p_tof2tof1_plot1];
                Int_t bin_tof2tof1_plot1 = tof_pid->GetXaxis()->FindBin(xp_tof2tof1_plot1);
                Float_t yp_tof2tof1_plot1 = tof_pid->GetBinContent(bin_tof2tof1_plot1);
            }
            Double_t par_tof2tof1_plot1[6];
            
            TF1 *g1_tof2tof1_plot1    = new TF1("g1_tof2tof1_plot1","gaus",32.4 , 32.9);
            TF1 *g2_tof2tof1_plot1    = new TF1("g2_tof2tof1_plot1","gaus",31.5 , 32.0);
            tof_pid->Fit(g1_tof2tof1_plot1,"R");
            tof_pid->Fit(g2_tof2tof1_plot1,"R+");
            g1_tof2tof1_plot1->GetParameters(&par_tof2tof1_plot1[0]);
            g2_tof2tof1_plot1->GetParameters(&par_tof2tof1_plot1[3]);
            Double_t time_proton_tof2tof1_plot1_co = g1_tof2tof1_plot1->GetParameter(1);
            Double_t time_pi_tof2tof1_plot1_co = g2_tof2tof1_plot1->GetParameter(1);
            Double_t DeltaT_tof2tof1_plot1_co= time_proton_tof2tof1_plot1_co - time_pi_tof2tof1_plot1_co;
            Double_t dmean_proton_tof2tof1_plot1_co= g1_tof2tof1_plot1->GetParError(1);
            Double_t dmean_pi_tof2tof1_plot1_co= g2_tof2tof1_plot1->GetParError(1);
            distance_tof2tof1_plot1_co= abs(((DeltaT_tof2tof1_plot1_co* nano_value)/((v_proton-v_pi)/(v_proton*v_pi))));
            
            calc_p_tof2tof1_plot1_co = sqrt((pow(DeltaT_tof2tof1_plot1_co,2) * pow(measured_d_tof2tof1_plot1_co,2) * pow(m_proton,2) * pow(nano_value,2))/(-4 * pow(DeltaT_tof2tof1_plot1_co,2) * pow(measured_d_tof2tof1_plot1_co,2) * pow(nano_value,2) + pow(DeltaT_tof2tof1_plot1_co,4) * pow(nano_value,4) * pow(c,2)) + (pow(DeltaT_tof2tof1_plot1_co,2) * pow(measured_d_tof2tof1_plot1_co,2) * pow(m_pi,2) * pow(nano_value,2))/(-4 * pow(DeltaT_tof2tof1_plot1_co,2) * pow(measured_d_tof2tof1_plot1_co,2) * pow(nano_value,2) + pow(DeltaT_tof2tof1_plot1_co,4) * pow(nano_value,4) * pow(c,2)) - (2*sqrt(pow(DeltaT_tof2tof1_plot1_co,2) * pow(measured_d_tof2tof1_plot1_co,4) * pow(nano_value,2)*(pow(measured_d_tof2tof1_plot1_co,2) * pow(m_proton,4) - 2 * pow(measured_d_tof2tof1_plot1_co,2) * pow(m_proton,2) * pow(m_pi,2) + pow(measured_d_tof2tof1_plot1_co,2) * pow(m_pi,4) + pow(DeltaT_tof2tof1_plot1_co,2) * pow(m_proton,2) * pow(m_pi,2) * pow(nano_value,2) * pow(c,2))))/(c*(-4 * pow(DeltaT_tof2tof1_plot1_co,2) * pow(measured_d_tof2tof1_plot1_co,2) * pow(nano_value,2) + pow(DeltaT_tof2tof1_plot1_co,4) * pow(nano_value,4) * pow(c,2))));

            c1->cd();
            //////////////////
            //  Fill graph ///
            /////////////////
            x[counter]=i;
            
            y_spr_sim_org[counter]=spr_sim_org*1000;
            y_spr_sim_true[counter]=spr_sim_true*1000;
            y_spr_sim_correction[counter]=spr_sim_correction*1000;
            y_spr_data_org[counter]=spr_data_org*1000;
            y_spr_data_sub[counter]=spr_data_sub*1000;
            y_spr_data_correction[counter]=spr_data_correction*1000;
            y_cangle_sim_org[counter]=cangle_sim_org;
            y_cangle_sim_org_true[counter]=cangle_sim_true;
            y_cangle_sim_org_correction[counter]=cangle_sim_correction;
            y_cangle_data_org[counter]=cangle_data_org;
            y_cangle_data_sub[counter]=cangle_data_sub;
            y_cangle_data_correction[counter]=cangle_data_correction;
            y_diff_true[counter]=sigma_p_diff_time_mctruth*5;
            y_diff_data[counter]=sigma_p_diff_time_data*5;
            y_diff_sim[counter]=sigma_p_diff_time_sim*5;
            y_mean_diff_data[counter]=mean_p_diff_time_data;
            y_mean_diff_sim[counter]=mean_p_diff_time_sim;
            y_mean_diff_true[counter]=mean_p_diff_time_mctruth;
            
            calc_mom->SetPoint(counter,i,calc_p_tof2tof1_plot1_co);
            calc_tof1tof2_distance->SetPoint(counter,i,distance_tof2tof1_plot1_co);
            
            // remove
            y_spr_data_cuts[counter]=fAngleP-spr_data_org;// warning switch between p and pi
            y_spr_sim_cuts[counter]=fAngleP-spr_sim_org;// warning
            y_spr_data_sub_cuts[counter]=fAngleP- spr_data_sub; // warning
            y_spr_sim_true_cuts[counter]= fAngleP- spr_sim_true;// warning
            

            
            std::cout<<"############  counter = "<< counter <<std::endl;
            std::cout<<"############  x[counter] = "<< x[counter] <<std::endl;
            std::cout<<"############  y_spr_data_org[counter] = "<< y_spr_data_org[counter] <<std::endl;
            gROOT->SetBatch(0);
            //method 1 return mean of the fit
            if(false) {
                y_yield_nph_sim[counter]=YieldGausFit_double(nph_sim);
                y_yield_p_nph_sim[counter]=YieldGausFit_double(p_nph_sim);
                y_yield_p_nph_good_sim[counter]=YieldGausFit_double(p_nph_good_sim);
                y_yield_p_nph_true_sim[counter]=YieldGausFit_double(p_nph_true_sim);
                y_yield_nph_data[counter]=YieldGausFit_double(nph_data);
                y_yield_p_nph_data[counter]=YieldGausFit_double(p_nph_data);
                y_yield_p_nph_good_data[counter]=YieldGausFit_double(p_nph_good_data);
                //~ y_yield_dac_hits_sys_cus_data[counter]=YieldGausFit_double(dac_hits_sys_cus_data);
            }
            //method 2 return mean of the histogram
            if(true) {
                y_yield_nph_sim[counter]=nph_sim->GetMean();
                y_yield_p_nph_sim[counter]=p_nph_sim->GetMean();
                y_yield_p_nph_good_sim[counter]=p_nph_good_sim->GetMean();
                y_yield_p_nph_true_sim[counter]=p_nph_true_sim->GetMean();
                y_yield_nph_data[counter]=nph_data->GetMean();
                y_yield_p_nph_data[counter]=p_nph_data->GetMean();
                y_yield_p_nph_good_data[counter]=p_nph_good_data->GetMean();
                
                y_yield_nph_sim_error[counter]=nph_sim->GetRMS();
                y_yield_p_nph_sim_error[counter]=p_nph_sim->GetRMS();
                y_yield_p_nph_good_sim_error[counter]=p_nph_good_sim->GetRMS();
                y_yield_p_nph_true_sim_error[counter]=p_nph_true_sim->GetRMS();
                y_yield_nph_data_error[counter]=nph_data->GetRMS();
                y_yield_p_nph_data_error[counter]=p_nph_data->GetRMS();
                y_yield_p_nph_good_data_error[counter]=p_nph_good_data->GetRMS();
                
                y_spr_data_correction_error[counter]=p_cherenkov_data_correction->GetRMSError();
            }
            counter++;
            
            recoAngle[counter]=cangle_sim_true;
            timeCut[counter]=sigma_p_diff_time_mctruth*5;
            chAngleCut[counter]=spr_sim_true*5;
            prtangle_vector[counter]= i;
            
            ////////////
            // hstack //
            ////////////
            hs = new THStack("hs","Stacked 1D histograms");
            hs2 = new THStack("hs2","Stacked 1D histograms");
            hs3 = new THStack("hs3","Stacked 1D histograms");
            hs4 = new THStack("hs4","Stacked 1D histograms");
            hs5 = new THStack("hs5","Stacked 1D histograms");
            hs6 = new THStack("hs6","Stacked 1D histograms");
            hs7 = new THStack("hs7","Stacked 1D histograms");
            hs6 = new THStack("hs6","Stacked 1D histograms");
            hs8 = new THStack("hs8","Stacked 1D histograms");
            hs9 = new THStack("hs9","Stacked 1D histograms");
            hs10 = new THStack("hs10","Stacked 1D histograms");
            
            hs->Add(p_cherenkov_bg_sim);
            hs->Add(p_cherenkov_data_sub);
            hs->Add(p_cherenkov_data_copy);
            hs2->Add(nph_sim);
            hs2->Add(p_nph_sim);
            hs2->Add(p_nph_good_sim);
            //hs2->Add(p_nph_true_sim);
            hs3->Add(p_diff_time_sim);
            hs3->Add(p_diff_time_bg_sim);
            hs3->Add(p_diff_time_data);
            hs3->Add(p_diff_time_mctruth);
            //hs3->Add(p_diff_time_data_sub);
            hs4->Add(p_photon_time_sim_calc);
            hs4->Add(p_photon_time_data_calc);
            hs5->Add(p_photon_time_sim);
            hs5->Add(p_photon_time_data);
            hs6->Add(nph_data);
            hs6->Add(p_nph_data);
            hs6->Add(p_nph_good_data);
            //hs6->Add(dac_hits_data);
            //hs6->Add(dac_hits_sys_cus_data);
            hs7->Add(p_cherenkov_bg_sim);
            hs7->Add(p_cherenkov_sim);
            hs7->Add(p_cherenkov_mc_same_path);
            hs8->Add(histo_photon_ambiguity_wo_data);
            hs8->Add(histo_photon_ambiguity_wt_data);
            hs8->Add(histo_photon_ambiguity_wtc_data);
            hs9->Add(histo_photon_ambiguity_wo_sim);
            hs9->Add(histo_photon_ambiguity_wt_sim);
            hs9->Add(histo_photon_ambiguity_wtc_sim);
            hs10->Add(hist_ambiguity_sim);
            hs10->Add(hist_ambiguity_data);
            
            ///////////////
            // TOF PID  ///
            ///////////////
            if (true) { // done
                gStyle->SetOptFit(0);
                gStyle->SetOptStat(0);
                prt_canvasAdd("r_tof_pid"+nid,800,400);
                TLegend *leg_tofpid= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                //leg_tofpid->SetHeader("cuts","C");
                tof_pid->SetTitle(Form("Polar angle %3.1f (data)", prtangle));
                leg_tofpid->AddEntry(tof_pid,"Ambiguity distribution sim","f");
                tof_pid->SetLineColor(kBlue);
                tof_pid->SetLineStyle(1);
                tof_pid->GetXaxis()->SetTitle("LE TOF2 - TOF1 [ns]");
                tof_pid->GetYaxis()->SetTitle("entries [#]");
                tof_pid->GetXaxis()->SetTitleSize(0.05);
                tof_pid->GetYaxis()->SetTitleSize(0.05);
                tof_pid->GetXaxis()->SetTitleOffset(0.9);
                tof_pid->GetYaxis()->SetTitleOffset(1.0);
                tof_pid->SetFillColor(kBlue);
                tof_pid->SetFillStyle(3003);
                tof_pid->Draw();
                //leg_tofpid->Draw();
                prt_canvasGet("r_tof_pid"+nid)->Update();
                TLine *line_proton_r = new TLine(0,0,0,1000);
                line_proton_r->SetX1(32.9);
                line_proton_r->SetX2(32.9);
                line_proton_r->SetY1(gPad->GetUymin());
                line_proton_r->SetY2(gPad->GetUymax());
                line_proton_r->SetLineColor(kRed);
                line_proton_r->Draw();
                TLine *line_proton_l = new TLine(0,0,0,1000);
                line_proton_l->SetX1(32.4);
                line_proton_l->SetX2(32.4);
                line_proton_l->SetY1(gPad->GetUymin());
                line_proton_l->SetY2(gPad->GetUymax());
                line_proton_l->SetLineColor(kRed);
                line_proton_l->Draw();
                TLine *line_pi_r = new TLine(0,0,0,1000);
                line_pi_r->SetX1(32);
                line_pi_r->SetX2(32);
                line_pi_r->SetY1(gPad->GetUymin());
                line_pi_r->SetY2(gPad->GetUymax());
                line_pi_r->SetLineColor(kBlue);
                line_pi_r->Draw();
                TLine *line_pi_l = new TLine(0,0,0,1000);
                line_pi_l->SetX1(31.5);
                line_pi_l->SetX2(31.5);
                line_pi_l->SetY1(gPad->GetUymin());
                line_pi_l->SetY2(gPad->GetUymax());
                line_pi_l->SetLineColor(kBlue);
                line_pi_l->Draw();
            }

            /////////////////////////////////
            //  ch MC truth B and Signal  ///
            /////////////////////////////////
            if(true) { // done
                gStyle->SetOptFit(0);
                gStyle->SetOptStat(0);
                gPad->UseCurrentStyle();
                TLegend *legend_signal= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                legend_signal->SetHeader("cherenkov angle (proton sim)","C");
                prt_canvasAdd("r_mctruth"+nid,800,400);
                p_cherenkov_sim->SetTitle(Form("Polar angle %3.1f", prtangle));
                legend_signal->AddEntry(p_cherenkov_sim,"Ambiguity distribution","f");
                legend_signal->AddEntry(p_cherenkov_mc_same_path," true path inside the prism ","f");
                legend_signal->AddEntry(p_cherenkov_bg_sim," combinatorial background ","f");
                hs7->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs7->Draw("nostack");
                hs7->GetXaxis()->SetTitle("#theta_{C} [rad]");
                hs7->GetYaxis()->SetTitle("entries [#]");
                legend_signal->Draw();
                TLine *line_chcut_r_copy = new TLine(0,0,0,1000);
                line_chcut_r_copy->SetX1(cangle_sim_true+5*spr_sim_true);
                line_chcut_r_copy->SetX2(cangle_sim_true+5*spr_sim_true);
                line_chcut_r_copy->SetY1(gPad->GetUymin());
                line_chcut_r_copy->SetY2(gPad->GetUymax());
                line_chcut_r_copy->SetLineColor(kMagenta);
                line_chcut_r_copy->Draw();
                TLine *line_chcut_l_copy = new TLine(0,0,0,1000);
                line_chcut_l_copy->SetX1(cangle_sim_true-5*spr_sim_true);
                line_chcut_l_copy->SetX2(cangle_sim_true-5*spr_sim_true);
                line_chcut_l_copy->SetY1(gPad->GetUymin());
                line_chcut_l_copy->SetY2(gPad->GetUymax());
                line_chcut_l_copy->SetLineColor(kMagenta);
                line_chcut_l_copy->Draw();
                
                //prt_canvasGet("r_mctruth"+nid)->Update();
                
                // TLine *line_true = new TLine(0,0,0,1000);
                // line_true->SetX1(fAngleP);
                // line_true->SetX2(fAngleP);
                // line_true->SetY1(gPad->GetUymin());
                // line_true->SetY2(gPad->GetUymax());
                // line_true->SetLineColor(kRed);
                // line_true->Draw();
                
                // TLine *line_sigm_r_true = new TLine(0,0,0,1000);
                // line_sigm_r_true->SetX1(out_array[6]);
                // line_sigm_r_true->SetX2(out_array[6]);
                // line_sigm_r_true->SetY1(gPad->GetUymin());
                // line_sigm_r_true->SetY2(gPad->GetUymax());
                // line_sigm_r_true->SetLineColor(kRed);
                // line_sigm_r_true->Draw();
                
                // TLine *line_sigm_l2_true = new TLine(0,0,0,1000);
                // line_sigm_l2_true->SetX1(out_array[5]);
                // line_sigm_l2_true->SetX2(out_array[5]);
                // line_sigm_l2_true->SetY1(gPad->GetUymin());
                // line_sigm_l2_true->SetY2(gPad->GetUymax());
                // line_sigm_l2_true->SetLineColor(kRed);
                // line_sigm_l2_true->Draw();
            }
            //////////////////////////////
            //  MC BG, Data, Data -BG  ///
            //////////////////////////////
            if(true) {
                gStyle->SetOptFit(0);
                gStyle->SetOptStat(0);
                gPad->UseCurrentStyle();
                TLegend * legend_bg_sub= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                legend_bg_sub->SetHeader("ambiguity distribution (proton)","C");
                legend_bg_sub->AddEntry(p_cherenkov_data_copy,"data","p");
                legend_bg_sub->AddEntry(p_cherenkov_data_sub,"data - BG","p");
                legend_bg_sub->AddEntry(p_cherenkov_bg_sim," BG ","f");
                prt_canvasAdd("r_bg_sub"+nid,800,400);
                hs->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs->Draw("nostack");
                hs->GetXaxis()->SetTitle("#theta_{C} [rad]");
                hs->GetYaxis()->SetTitle("entries [#]");
                legend_bg_sub->Draw();
                prt_canvasGet("r_bg_sub"+nid)->Update();
                TLine *line_chcut_r = new TLine(0,0,0,1000);
                line_chcut_r->SetX1(cangle_sim_true+5*spr_sim_true);
                line_chcut_r->SetX2(cangle_sim_true+5*spr_sim_true);
                line_chcut_r->SetY1(gPad->GetUymin());
                line_chcut_r->SetY2(gPad->GetUymax());
                line_chcut_r->SetLineColor(kMagenta);
                line_chcut_r->Draw();
                TLine *line_chcut_l = new TLine(0,0,0,1000);
                line_chcut_l->SetX1(cangle_sim_true-5*spr_sim_true);
                line_chcut_l->SetX2(cangle_sim_true-5*spr_sim_true);
                line_chcut_l->SetY1(gPad->GetUymin());
                line_chcut_l->SetY2(gPad->GetUymax());
                line_chcut_l->SetLineColor(kMagenta);
                line_chcut_l->Draw();
            }
            ///////////////////////////////////////////
            //  diff MC BG, diff Data, diffData -BG  //
            ///////////////////////////////////////////
            if (true) { // done
                gStyle->SetOptFit(0);
                gStyle->SetOptStat(0);
                TLegend *legend_diff_bg_sub= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                legend_diff_bg_sub->SetHeader("time difference distribution (proton)","C");
                legend_diff_bg_sub->AddEntry(p_diff_time_data,"data","p");
                //legend_diff_bg_sub->AddEntry(p_diff_time_data_sub,"data - BG","p");
                legend_diff_bg_sub->AddEntry(p_diff_time_bg_sim," sim BG ","f");
                legend_diff_bg_sub->AddEntry(p_diff_time_sim," sim ","f");
                legend_diff_bg_sub->AddEntry(p_diff_time_mctruth," sim true path inside prism ","f");
                prt_canvasAdd("r_time_diff"+nid,800,400);
                hs3->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs3->Draw("nostack");
                hs3->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
                hs3->GetYaxis()->SetTitle("entries [#]");
                legend_diff_bg_sub->Draw();
                prt_canvasGet("r_time_diff"+nid)->Update();
                //prt_canvasGet("r_time_diff"+nid)->Modified();
                TLine *line_sigm_diff_r = new TLine(0,0,0,1000);
                line_sigm_diff_r->SetX1(/*out_array[13]+*/sigma_p_diff_time_mctruth*5);
                line_sigm_diff_r->SetX2(/*out_array[13]+*/sigma_p_diff_time_mctruth*5);
                line_sigm_diff_r->SetY1(gPad->GetUymin());
                line_sigm_diff_r->SetY2(gPad->GetUymax());
                line_sigm_diff_r->SetLineColor(kMagenta);
                line_sigm_diff_r->Draw();
                TLine *line_sigm_diff_l = new TLine(0,0,0,1000);
                line_sigm_diff_l->SetX1(/*out_array[13]*/-sigma_p_diff_time_mctruth*5);
                line_sigm_diff_l->SetX2(/*out_array[13]*/-sigma_p_diff_time_mctruth*5);
                line_sigm_diff_l->SetY1(gPad->GetUymin());
                line_sigm_diff_l->SetY2(gPad->GetUymax());
                line_sigm_diff_l->SetLineColor(kMagenta);
                line_sigm_diff_l->Draw();
                prt_canvasGet("r_time_diff"+nid)->Update();
            }
            if (false) {
                TLegend *legend_photon_time_calc= new TLegend(0.556391, 0.712, 0.890977, 0.874667);
                legend_photon_time_calc->SetHeader("calculated time match (proton)", "C");
                legend_photon_time_calc->AddEntry(p_photon_time_sim_calc,"calculated time (sim)","f");
                legend_photon_time_calc->AddEntry(p_photon_time_data_calc,"calculated time (data) ","p");
                prt_canvasAdd("r_photon_calc_time"+nid,800,400);
                hs4->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs4->Draw("nostack");
                hs4->GetXaxis()->SetTitle("t_{calc} [ns]");
                hs4->GetYaxis()->SetTitle("entries [#]");
                legend_photon_time_calc->Draw();
            }
            if (false) {
                TLegend *legend_photon_time=  new TLegend(0.556391, 0.712, 0.890977, 0.874667);
                legend_photon_time->SetHeader("measured time match (proton)","C");
                legend_photon_time->AddEntry(p_photon_time_sim,"measured time (sim)","f");
                legend_photon_time->AddEntry(p_photon_time_data,"measured time (data) ","p");
                prt_canvasAdd("r_photon_measured_time"+nid,800,400);
                hs5->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs5->Draw("nostack");
                hs5->GetXaxis()->SetTitle("t_{measured} [ns]");
                hs5->GetYaxis()->SetTitle("entries [#]");
                legend_photon_time->Draw();
            }
            /////////////////////////
            // Draw photon yield  ///
            /////////////////////////
            if (false) {
                TLegend * legend_nph_sim= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
                legend_nph_sim->SetHeader("photon yield (proton sim)","C");
                prt_canvasAdd("r_fnHits_sim"+nid,800,400);
                nph_sim->SetTitle(Form("Polar angle %3.1f (sim)", prtangle));
                legend_nph_sim->AddEntry(nph_sim," DIRC hits without cuts  ","l");
                legend_nph_sim->AddEntry(p_nph_sim,"DIRC hits with time cut ","l");
                legend_nph_sim->AddEntry(p_nph_good_sim,"DIRC hits with #theta_{C} and time cuts","l");
                //legend_nph_sim->AddEntry(p_nph_true_sim," DIRC hits true path inside prism","l");
                hs2->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs2->Draw("nostack");
                hs2->GetYaxis()->SetTitle("entries [#]");
                hs2->GetXaxis()->SetTitle("number of photon per track [#]");
                legend_nph_sim->Draw();
            }
            if (false) {
                TLegend * legend_nph_data= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
                legend_nph_data->SetHeader("photon yield (proton data)","C");
                prt_canvasAdd("r_fnHits_data"+nid,800,400);
                hs6->SetTitle(Form("Polar angle %3.1f (data)", prtangle));
                legend_nph_data->AddEntry(nph_data,"DIRC hits without cuts ","l");
                legend_nph_data->AddEntry(p_nph_data,"DIRC hits with time cut ","l");
                legend_nph_data->AddEntry(p_nph_good_data,"DIRC hits with #theta_{C} and time cuts ","l");
                //legend_nph_data->AddEntry(dac_hits_data,"DAC hits without cuts","f");
                //legend_nph_data->AddEntry(dac_hits_sys_cus_data,"DAC hits with auxiliary detectors cuts and PID cut(proton data)","f");
                hs6->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs6->Draw("nostack");
                hs6->GetYaxis()->SetTitle("entries [#]");
                hs6->GetXaxis()->SetTitle("number of photon per track [#]");
                legend_nph_data->Draw();
            }
            if (false) {
                TLegend * legend_ambiguity_data= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
                legend_ambiguity_data->SetHeader("number of photon solutions (proton data)","C");
                prt_canvasAdd("r_ambiguity_data"+nid,800,400);
                hs8->SetTitle(Form("Polar angle %3.1f (data)", prtangle));
                legend_ambiguity_data->AddEntry(histo_photon_ambiguity_wo_sim,"number of photon solutions without cuts ","l");
                legend_ambiguity_data->AddEntry(histo_photon_ambiguity_wt_sim,"number of photon solutions with time cut ","l");
                legend_ambiguity_data->AddEntry(histo_photon_ambiguity_wtc_sim,"number of photon solutions with #theta_{C} and time cuts ","l");
                hs8->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs8->Draw("nostack");
                hs8->GetYaxis()->SetTitle("entries [#]");
                hs8->GetXaxis()->SetTitle("number of photon solutions  [#]");
                legend_ambiguity_data->Draw();
            }
            if (false) {
                TLegend * legend_ambiguity_sim= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
                legend_ambiguity_sim->SetHeader("number of photon solutions (proton data)","C");
                prt_canvasAdd("r_ambiguity_sim"+nid,800,400);
                hs9->SetTitle(Form("Polar angle %3.1f (data)", prtangle));
                legend_ambiguity_sim->AddEntry(histo_photon_ambiguity_wo_sim,"number of photon solutions without cuts ","l");
                legend_ambiguity_sim->AddEntry(histo_photon_ambiguity_wt_sim,"number of photon solutions with time cut ","l");
                legend_ambiguity_sim->AddEntry(histo_photon_ambiguity_wtc_sim,"number of photon solutions with #theta_{C} and time cuts ","l");
                hs9->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs9->Draw("nostack");
                hs9->GetYaxis()->SetTitle("entries [#]");
                hs9->GetXaxis()->SetTitle("number of photon solutions  [#]");
                legend_ambiguity_sim->Draw();
            }
            if (false) {
                TLegend * legend_lut_ambiguity= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
                legend_lut_ambiguity->SetHeader("number of LUT photon solutions for each pixel","C");
                prt_canvasAdd("r_lut_ambiguity"+nid,800,400);
                hs10->SetTitle(Form("Polar angle %3.1f (data)", prtangle));
                hs10->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs10->Draw("nostack");
                hs10->GetYaxis()->SetTitle("entries [#]");
                hs10->GetXaxis()->SetTitle("number of pixel [#]");
                legend_lut_ambiguity->Draw();
            }
             c1->cd();
        }
    }
    ///////////////////
    ///// part III ////
    ///////////////////
    if(false) {
        calc_mom->SetLineColor(kBlack);
        calc_mom->SetMarkerColor(kBlack);
        calc_mom->SetMarkerStyle(21);
        calc_mom->SetMarkerSize(0.7);
        calc_mom->GetYaxis()->SetRangeUser(0,20);
        calc_mom->GetXaxis()->SetLabelSize(0.05);
        calc_mom->GetXaxis()->SetTitleSize(0.06);
        calc_mom->GetXaxis()->SetTitleOffset(0.84);
        
        calc_tof1tof2_distance->SetLineColor(kBlack);
        calc_tof1tof2_distance->SetMarkerColor(kBlack);
        calc_tof1tof2_distance->SetMarkerStyle(21);
        calc_tof1tof2_distance->SetMarkerSize(0.7);
        calc_tof1tof2_distance->GetYaxis()->SetRangeUser(0,20);
        calc_tof1tof2_distance->GetXaxis()->SetLabelSize(0.05);
        calc_tof1tof2_distance->GetXaxis()->SetTitleSize(0.06);
        calc_tof1tof2_distance->GetXaxis()->SetTitleOffset(0.84);
        
        ////////////////////////////////////////////
        // SPR after/befor Background subtraction///
        ////////////////////////////////////////////
        
        graph_spr_sim = new TGraph(n,x,y_spr_sim_org);
        graph_spr_data_sub = new TGraph(n,x,y_spr_data_sub);
        graph_spr_sim_true = new TGraph(n,x,y_spr_sim_true);
        graph_spr_data_correction = new TGraph(n,x,y_spr_data_correction);
        graph_spr_sim_correction = new TGraph(n,x,y_spr_sim_correction);
        graph_spr_data_org = new TGraph(n,x,y_spr_data_org);
        
        graph_spr_data_cuts = new TGraph(n,x,y_spr_data_cuts);
        graph_spr_sim_cuts = new TGraph(n,x,y_spr_sim_cuts);
        graph_spr_data_sub_cuts = new TGraph(n,x,y_spr_data_sub_cuts);
        graph_spr_sim_true_cuts = new TGraph(n,x,y_spr_sim_true_cuts);
        
        graph_diff_true = new TGraph(n,x,y_diff_true);
        graph_diff_data = new TGraph(n,x,y_diff_data);
        graph_diff_sim = new TGraph(n,x,y_diff_sim);
        
        graph_diff_data_mean = new TGraph(n,x,y_mean_diff_data);
        graph_diff_sim_mean = new TGraph(n,x,y_mean_diff_sim);
        
        graph_diff_true_mean = new TGraph(n,x,y_mean_diff_true);
        graph_diff_true_sim = new TGraph(n,x,y_mean_diff_sim);
        
        graph_cangle_data_sub = new TGraph(n,x,y_cangle_data_sub);
        graph_cangle_sim = new TGraph(n,x,y_cangle_sim_org);
        graph_cangle_sim_true = new TGraph(n,x,y_cangle_sim_org_true);
        graph_cangle_data_correction = new TGraph(n,x,y_cangle_data_correction);
        graph_cangle_sim_correction = new TGraph(n,x,y_cangle_sim_org_correction);
        graph_cangle_data_org = new TGraph(n,x,y_cangle_data_org);
        
        //        TGraphErrors *graph_yield_DIRC_wo_cuts_sim = new TGraphErrors(n,x,y_yield_nph_sim,0, y_yield_nph_sim_error);
        //        TGraphErrors *graph_yield_DIRC_wt_cuts_sim = new TGraphErrors(n,x,y_yield_p_nph_sim,0, y_yield_p_nph_sim_error);
        //        TGraphErrors *graph_yield_DIRC_wtc_cuts_sim = new TGraphErrors(n,x,y_yield_p_nph_good_sim,0, y_yield_p_nph_good_sim_error);
        //        TGraphErrors *graph_yield_DIRC_true_sim = new TGraphErrors(n,x,y_yield_p_nph_true_sim,0, y_yield_p_nph_true_sim_error);
        //        TGraphErrors *graph_yield_DIRC_wo_cuts_data = new TGraphErrors(n,x,y_yield_nph_data,0, y_yield_nph_data_error);
        //        TGraphErrors *graph_yield_DIRC_wt_cuts_data = new TGraphErrors(n,x,y_yield_p_nph_data,0, y_yield_p_nph_data_error);
        //        TGraphErrors *graph_yield_DIRC_wtc_cuts_data = new TGraphErrors(n,x,y_yield_p_nph_good_data,0, y_yield_p_nph_good_data_error);

        graph_yield_DIRC_wo_cuts_sim = new TGraph(n,x,y_yield_nph_sim);
        graph_yield_DIRC_wt_cuts_sim = new TGraph(n,x,y_yield_p_nph_sim);
        graph_yield_DIRC_wtc_cuts_sim = new TGraph(n,x,y_yield_p_nph_good_sim);
        graph_yield_DIRC_true_sim     = new TGraph(n,x,y_yield_p_nph_true_sim);
        graph_yield_DIRC_wo_cuts_data = new TGraph(n,x,y_yield_nph_data);
        graph_yield_DIRC_wt_cuts_data = new TGraph(n,x,y_yield_p_nph_data);
        graph_yield_DIRC_wtc_cuts_data = new TGraph(n,x,y_yield_p_nph_good_data);
        
        //TGraphErrors *graph_spr_data_correction = new TGraphErrors(n,x,y_spr_data_correction,0,y_spr_data_correction_error);
        
        graph_spr_sim->SetLineColor(1);
        graph_spr_sim->SetLineWidth(1);
        graph_spr_sim->SetMarkerColor(1);
        graph_spr_sim->SetMarkerStyle(21);
        graph_spr_sim->SetLineStyle(1);
        
        graph_spr_data_sub->SetLineColor(2);
        graph_spr_data_sub->SetLineWidth(1);
        graph_spr_data_sub->SetMarkerColor(2);
        graph_spr_data_sub->SetMarkerStyle(21);
        
        graph_spr_sim_true->SetLineColor(2);
        graph_spr_sim_true->SetLineWidth(1);
        graph_spr_sim_true->SetMarkerColor(2);
        graph_spr_sim_true->SetMarkerStyle(21);
        graph_spr_sim_true->SetLineStyle(1);
        
        //////////////////////////////////
        graph_spr_data_cuts->SetLineColor(1);
        graph_spr_data_cuts->SetLineWidth(1);
        graph_spr_data_cuts->SetMarkerColor(1);
        graph_spr_data_cuts->SetMarkerStyle(21);
        
        graph_spr_sim_cuts->SetLineColor(1);
        graph_spr_sim_cuts->SetLineWidth(1);
        graph_spr_sim_cuts->SetMarkerColor(1);
        graph_spr_sim_cuts->SetMarkerStyle(21);
        graph_spr_sim_cuts->SetLineStyle(1);
        
        graph_spr_data_sub_cuts->SetLineColor(2);
        graph_spr_data_sub_cuts->SetLineWidth(1);
        graph_spr_data_sub_cuts->SetMarkerColor(2);
        graph_spr_data_sub_cuts->SetMarkerStyle(21);
        
        graph_spr_sim_true_cuts->SetLineColor(2);
        graph_spr_sim_true_cuts->SetLineWidth(1);
        graph_spr_sim_true_cuts->SetMarkerColor(2);
        graph_spr_sim_true_cuts->SetMarkerStyle(21);
        graph_spr_sim_true_cuts->SetLineStyle(1);
        
        graph_spr_data_correction->SetLineColor(kBlack);
        graph_spr_data_correction->SetLineWidth(1);
        graph_spr_data_correction->SetMarkerColor(kBlack);
        graph_spr_data_correction->SetMarkerStyle(21);
        graph_spr_data_correction->SetLineStyle(1);
        
        graph_spr_sim_correction->SetLineColor(kGreen);
        graph_spr_sim_correction->SetLineWidth(1);
        graph_spr_sim_correction->SetMarkerColor(kGreen);
        graph_spr_sim_correction->SetMarkerStyle(21);
        graph_spr_sim_correction->SetLineStyle(1);
        
        graph_spr_data_org->SetLineColor(kRed);
        graph_spr_data_org->SetLineWidth(1);
        graph_spr_data_org->SetMarkerColor(kRed);
        graph_spr_data_org->SetMarkerStyle(21);
        graph_spr_data_org->SetLineStyle(1);
        
        ///////////////////////////////
        graph_diff_true->SetLineColor(2);
        graph_diff_true->SetLineWidth(1);
        graph_diff_true->SetMarkerColor(2);
        graph_diff_true->SetMarkerStyle(21);
        
        graph_diff_data->SetLineColor(1);
        graph_diff_data->SetLineWidth(1);
        graph_diff_data->SetMarkerColor(1);
        graph_diff_data->SetMarkerStyle(21);
        
        graph_diff_sim->SetLineColor(kMagenta);
        graph_diff_sim->SetLineWidth(1);
        graph_diff_sim->SetMarkerColor(kMagenta);
        graph_diff_sim->SetMarkerStyle(21);
        ////////////////////////////////////
        graph_diff_data_mean->SetLineColor(1);
        graph_diff_data_mean->SetLineWidth(1);
        graph_diff_data_mean->SetMarkerColor(1);
        graph_diff_data_mean->SetMarkerStyle(21);
        
        graph_diff_sim_mean->SetLineColor(6);
        graph_diff_sim_mean->SetLineWidth(1);
        graph_diff_sim_mean->SetMarkerColor(6);
        graph_diff_sim_mean->SetMarkerStyle(21);
        
        graph_diff_true_mean->SetLineColor(2);
        graph_diff_true_mean->SetLineWidth(1);
        graph_diff_true_mean->SetMarkerColor(2);
        graph_diff_true_mean->SetMarkerStyle(21);
        
        graph_cangle_data_sub->SetLineColor(1);
        graph_cangle_data_sub->SetLineWidth(1);
        graph_cangle_data_sub->SetMarkerColor(1);
        graph_cangle_data_sub->SetMarkerStyle(21);
        
        graph_cangle_data_correction->SetLineColor(kBlack);
        graph_cangle_data_correction->SetLineWidth(1);
        graph_cangle_data_correction->SetMarkerColor(kBlack);
        graph_cangle_data_correction->SetMarkerStyle(21);
        graph_cangle_data_correction->SetLineStyle(1);
        
        graph_cangle_sim_correction->SetLineColor(kGreen);
        graph_cangle_sim_correction->SetLineWidth(1);
        graph_cangle_sim_correction->SetMarkerColor(kGreen);
        graph_cangle_sim_correction->SetMarkerStyle(21);
        graph_cangle_sim_correction->SetLineStyle(1);
        
        graph_cangle_data_org->SetLineColor(kRed);
        graph_cangle_data_org->SetLineWidth(1);
        graph_cangle_data_org->SetMarkerColor(kRed);
        graph_cangle_data_org->SetMarkerStyle(21);
        graph_cangle_data_org->SetLineStyle(1);
        
        graph_cangle_sim_true->SetLineColor(1);
        graph_cangle_sim_true->SetLineWidth(1);
        graph_cangle_sim_true->SetMarkerColor(1);
        graph_cangle_sim_true->SetMarkerStyle(21);
        graph_cangle_sim_true->SetLineStyle(1);
        
        graph_cangle_sim->SetLineColor(2);
        graph_cangle_sim->SetLineWidth(1);
        graph_cangle_sim->SetMarkerColor(2);
        graph_cangle_sim->SetMarkerStyle(21);
        graph_cangle_sim->SetLineStyle(1);
        
        graph_yield_DIRC_wo_cuts_sim->SetLineColor(kBlack);
        graph_yield_DIRC_wt_cuts_sim->SetLineColor(kBlue);
        graph_yield_DIRC_wtc_cuts_sim->SetLineColor(kMagenta);
        graph_yield_DIRC_true_sim->SetLineColor(kRed);
        graph_yield_DIRC_wo_cuts_data->SetLineColor(kRed);
        graph_yield_DIRC_wt_cuts_data->SetLineColor(kGreen);
        graph_yield_DIRC_wtc_cuts_data->SetLineColor(kCyan);
        //graph_yield_sys_wo_cuts_data->SetLineColor(8);
        
        graph_yield_DIRC_wo_cuts_sim->SetMarkerColor(kBlack);
        graph_yield_DIRC_wt_cuts_sim->SetMarkerColor(kBlue);
        graph_yield_DIRC_wtc_cuts_sim->SetMarkerColor(kMagenta);
        graph_yield_DIRC_true_sim->SetMarkerColor(kRed);
        graph_yield_DIRC_wo_cuts_data->SetMarkerColor(kRed);
        graph_yield_DIRC_wt_cuts_data->SetMarkerColor(kGreen);
        graph_yield_DIRC_wtc_cuts_data->SetMarkerColor(kCyan);
        //graph_yield_sys_wo_cuts_data->SetMarkerColor(8);
        
        graph_yield_DIRC_wo_cuts_sim->SetMarkerStyle(21);
        graph_yield_DIRC_wt_cuts_sim->SetMarkerStyle(21);
        graph_yield_DIRC_wtc_cuts_sim->SetMarkerStyle(21);
        graph_yield_DIRC_wo_cuts_data->SetMarkerStyle(21);
        graph_yield_DIRC_wt_cuts_data->SetMarkerStyle(21);
        graph_yield_DIRC_wtc_cuts_data->SetMarkerStyle(21);
        graph_yield_DIRC_true_sim->SetMarkerStyle(kRed);
        //graph_yield_sys_wo_cuts_data->SetMarkerStyle(21);
        
        graph_yield_DIRC_wo_cuts_sim->SetLineWidth(1);
        graph_yield_DIRC_wt_cuts_sim->SetLineWidth(1);
        graph_yield_DIRC_wtc_cuts_sim->SetLineWidth(1);
        graph_yield_DIRC_wo_cuts_data->SetLineWidth(1);
        graph_yield_DIRC_wt_cuts_data->SetLineWidth(1);
        graph_yield_DIRC_wtc_cuts_data->SetLineWidth(1);
        graph_yield_DIRC_true_sim->SetLineWidth(1);
        //graph_yield_sys_wo_cuts_data->SetLineWidth(1);
        
        graph_yield_DIRC_true_sim->SetLineStyle(1);
        
        if(false) { // not for website
            prt_canvasAdd("r_diff_time_shift",800,400);
            TLegend *leg_diff_mean = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_diff_mean->SetFillColor(0);
            //leg->SetHeader("test legend");
            leg_diff_mean->AddEntry(graph_diff_true_mean, "simulation true path iside prism", "lp");
            leg_diff_mean->AddEntry(graph_diff_sim_mean, "simulation", "lp");
            leg_diff_mean->AddEntry(graph_diff_data_mean, "test beam data", "lp");
            TMultiGraph *mg_diff_mean = new TMultiGraph();
            mg_diff_mean->Add(graph_diff_data_mean);
            mg_diff_mean->Add(graph_diff_sim_mean);
            mg_diff_mean->Add(graph_diff_true_mean);
            mg_diff_mean->SetTitle("time difference shift (should be before applaying any correction on hit time); #theta [degree]; t_{calc}-t_{measured} shift [ns]");
            mg_diff_mean->Draw("APL");
            leg_diff_mean->Draw();
        }
        if(false) { // not for website
            prt_canvasAdd("r_diff_time",800,400);
            TLegend *leg_diff = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_diff->SetFillColor(0);
            //leg->SetHeader("test legend");
            leg_diff->AddEntry(graph_diff_true, "true path inside prism ", "lp");
            leg_diff->AddEntry(graph_diff_data, "data ", "lp");
            leg_diff->AddEntry(graph_diff_sim, "simulation ", "lp");
            TMultiGraph *mg_diff = new TMultiGraph();
            mg_diff->Add(graph_diff_true);
            mg_diff->Add(graph_diff_data);
            mg_diff->Add(graph_diff_sim);
            mg_diff->SetTitle("3 #sigma time cut (should be befor time cut and after c cut);#theta [degree]; |time cut| [ns]");
            mg_diff->Draw("APL");
            leg_diff->Draw();
        }
        if(false) { // not for website
            prt_canvasAdd("r_cherenkov_angle_cuts",800,400);
            TLegend *leg_c_cuts = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_c_cuts->SetFillColor(0);
            leg_c_cuts->AddEntry(graph_spr_data_cuts, "data ", "lp");
            leg_c_cuts->AddEntry(graph_spr_data_sub_cuts, "Data - BG", "lp");
            leg_c_cuts->AddEntry(graph_spr_sim_cuts, "simulation ", "lp");
            leg_c_cuts->AddEntry(graph_spr_sim_true_cuts, "true path inside prism ", "lp");
            TMultiGraph *mg_c_cuts = new TMultiGraph();
            mg_c_cuts->Add(graph_spr_data_cuts);
            mg_c_cuts->Add(graph_spr_data_sub_cuts);
            mg_c_cuts->Add(graph_spr_sim_cuts);
            mg_c_cuts->Add(graph_spr_sim_true_cuts);
            mg_c_cuts->SetTitle("3 #sigma cherenkov cut (should be applayed at the end with the time cuts);#theta [degree]; |Cherenckove angle cut| [mrad]");
            mg_c_cuts->Draw("APL");
            leg_c_cuts->Draw();
        }
        //////////////
        //// cangle //
        //////////////
        if(true) {
            prt_canvasAdd("r_graph_cangle",800,400);
            TLegend *leg_cangle_data_sub = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_cangle_data_sub->SetHeader("(proton #theta_{c} data)","C");
            leg_cangle_data_sub->AddEntry(graph_cangle_data_org, "Data befor correction", "lp");
            //leg_cangle_data_sub->AddEntry(graph_cangle_data_sub, "Data - BG", "lp");
            //leg_cangle_data_sub->AddEntry(graph_cangle_sim, "Sim", "lp");
            //leg_cangle_data_sub->AddEntry(graph_cangle_sim_true, "Sim true inside prism", "lp");
            leg_cangle_data_sub->AddEntry(graph_cangle_data_correction, "Data after correction", "lp");
            
            TMultiGraph *mg_cangle_data_sub = new TMultiGraph();
            //mg_cangle_data_sub->Add(graph_cangle_data_sub);
            mg_cangle_data_sub->Add(graph_cangle_data_org);
            //mg_cangle_data_sub->Add(graph_cangle_sim);
            //mg_cangle_data_sub->Add(graph_cangle_sim_true);
            mg_cangle_data_sub->Add(graph_cangle_data_correction);
            
            mg_cangle_data_sub->SetTitle("reco cherenkov angle ; #theta [degree]; #theta_{C} [rad]");
            mg_cangle_data_sub->Draw("APL");
            mg_cangle_data_sub->GetHistogram()->GetYaxis()->SetRangeUser(0.8,0.85);
            leg_cangle_data_sub->Draw();
            prt_canvasGet("r_graph_cangle")->Update();
            TLine *line_ch_p = new TLine(0,0,0,1000);
            line_ch_p->SetY1(fAngleP);
            line_ch_p->SetY2(fAngleP);
            line_ch_p->SetX1(gPad->GetUxmin());
            line_ch_p->SetX2(gPad->GetUxmax());
            line_ch_p->SetLineColor(kRed);
            line_ch_p->Draw();
            
            TLine *line_ch_pi = new TLine(0,0,0,1000);
            line_ch_pi->SetY1(fAnglePi);
            line_ch_pi->SetY2(fAnglePi);
            line_ch_pi->SetX1(gPad->GetUxmin());
            line_ch_pi->SetX2(gPad->GetUxmax());
            line_ch_pi->SetLineColor(kBlue);
            line_ch_pi->Draw();
        }
        
        //////////////
        //// spr /////
        //////////////
        if(true){
            graph_spr_data_org->SetLineColor(kBlack);
            graph_spr_data_org->SetMarkerColor(kBlack);
            graph_spr_sim->SetLineColor(kRed);
            graph_spr_sim->SetMarkerColor(kRed);
            
            graph_spr_data_sub->SetLineColor(kBlue);
            graph_spr_data_sub->SetMarkerColor(kBlue);
            graph_spr_sim_true->SetLineColor(kGreen);
            graph_spr_sim_true->SetMarkerColor(kGreen);
            
            graph_spr_data_correction->SetLineColor(kMagenta);
            graph_spr_data_correction->SetMarkerColor(kMagenta);
            graph_spr_sim_correction->SetLineColor(kCyan);
            graph_spr_sim_correction->SetMarkerColor(kCyan);
            
            TMultiGraph *mg_spr_1 = new TMultiGraph();
            TMultiGraph *mg_spr_2 = new TMultiGraph();
            TMultiGraph *mg_spr_3 = new TMultiGraph();
            
            mg_spr_1->Add(graph_spr_data_org);
            mg_spr_1->Add(graph_spr_sim);
            TLegend *leg_spr1 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_spr1->SetHeader("photon spr (proton)","C");
            leg_spr1->AddEntry(graph_spr_data_org, "Data without corrections ", "lp");
            leg_spr1->AddEntry(graph_spr_sim, "Sim without corrections", "lp");
            
            mg_spr_2->Add(graph_spr_data_org);
            mg_spr_2->Add(graph_spr_data_correction);
            TLegend *leg_spr2 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_spr2->SetHeader("photon spr (proton)","C");
            leg_spr2->AddEntry(graph_spr_data_org, "Data without corrections ", "lp");
            leg_spr2->AddEntry(graph_spr_data_correction, "Data with corrections", "lp");
            
            mg_spr_3->Add(graph_spr_data_correction);
            mg_spr_3->Add(graph_spr_sim_correction);
            TLegend *leg_spr3 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_spr3->SetHeader("photon spr (proton)","C");
            leg_spr3->AddEntry(graph_spr_data_correction, "Data with corrections ", "lp");
            leg_spr3->AddEntry(graph_spr_sim_correction, "Sim with corrections", "lp");
            
            prt_canvasAdd("r_spr_1",800,400);
            mg_spr_1->SetTitle(" spr ;#theta [degree]; SPR [rad]");
            mg_spr_1->Draw("APL");
            leg_spr1->Draw();
            mg_spr_1->GetHistogram()->GetYaxis()->SetRangeUser(0.,20);
            gPad->Modified();
            gPad->Update();
            prt_canvasAdd("r_spr_2",800,400);
            mg_spr_2->SetTitle(" spr ;#theta [degree]; SPR [rad]");
            mg_spr_2->Draw("APL");
            leg_spr2->Draw();
            mg_spr_2->GetHistogram()->GetYaxis()->SetRangeUser(0.,20);
            gPad->Modified();
            gPad->Update();
            prt_canvasAdd("r_spr_3",800,400);
            mg_spr_3->SetTitle(" spr ;#theta [degree]; SPR [rad]");
            mg_spr_3->Draw("APL");
            leg_spr3->Draw();
            mg_spr_3->GetHistogram()->GetYaxis()->SetRangeUser(0.,20);
            gPad->Modified();
            gPad->Update();
        }
        
        //////////////
        //// yield ///
        //////////////
        if(true) {
            prt_canvasAdd("r_photon_yield_all",800,400);
            TLegend *leg_yield = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_yield->SetHeader("photon yield (proton)","C");
            leg_yield->AddEntry(graph_yield_DIRC_wo_cuts_sim, "graph_yield_DIRC_wo_cuts_sim ", "lp");
            //leg_yield->AddEntry(graph_yield_DIRC_wt_cuts_sim, "graph_yield_DIRC_wt_cuts_sim", "lp");
            leg_yield->AddEntry(graph_yield_DIRC_wtc_cuts_sim, "graph_yield_DIRC_wtc_cuts_sim ", "lp");
            leg_yield->AddEntry(graph_yield_DIRC_wo_cuts_data, "graph_yield_DIRC_wo_cuts_data", "lp");
            //leg_yield->AddEntry(graph_yield_DIRC_wt_cuts_data, "graph_yield_DIRC_wt_cuts_data", "lp");
            leg_yield->AddEntry(graph_yield_DIRC_wtc_cuts_data, "graph_yield_DIRC_wtc_cuts_data", "lp");
            //leg_yield->AddEntry(graph_yield_DIRC_true_sim, "graph_yield_DIRC_true_sim", "lp");
            //leg_yield->AddEntry(graph_yield_sys_wo_cuts_data, "graph_yield_sys_wo_cuts_data", "lp");
            
            TMultiGraph *mg_yield = new TMultiGraph();
            mg_yield->Add(graph_yield_DIRC_wo_cuts_sim);
            //mg_yield->Add(graph_yield_DIRC_wt_cuts_sim);
            mg_yield->Add(graph_yield_DIRC_wtc_cuts_sim);
            mg_yield->Add(graph_yield_DIRC_wo_cuts_data);
            //mg_yield->Add(graph_yield_DIRC_wt_cuts_data);
            mg_yield->Add(graph_yield_DIRC_wtc_cuts_data);
            //mg_yield->Add(graph_yield_DIRC_true_sim);
            //mg_yield->Add(graph_yield_sys_wo_cuts_data);
            mg_yield->SetTitle("photon yield ;#theta [degree]; count [#]");
            mg_yield->Draw("APL");
            leg_yield->Draw();
            
            TMultiGraph *mg_yield_1 = new TMultiGraph();
            TMultiGraph *mg_yield_2 = new TMultiGraph();
            TMultiGraph *mg_yield_3 = new TMultiGraph();
            
            mg_yield_1->Add(graph_yield_DIRC_wo_cuts_sim);
            mg_yield_1->Add(graph_yield_DIRC_wo_cuts_data);
            TLegend *leg_yield1 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_yield1->SetHeader("photon yield (proton)","C");
            leg_yield1->AddEntry(graph_yield_DIRC_wo_cuts_sim, "graph_yield_DIRC_wo_cuts_sim ", "lp");
            leg_yield1->AddEntry(graph_yield_DIRC_wo_cuts_data, "graph_yield_DIRC_wo_cuts_data", "lp");
            
            mg_yield_2->Add(graph_yield_DIRC_wt_cuts_sim);
            mg_yield_2->Add(graph_yield_DIRC_wt_cuts_data);
            TLegend *leg_yield2 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_yield2->SetHeader("photon yield (proton)","C");
            leg_yield2->AddEntry(graph_yield_DIRC_wt_cuts_sim, "graph_yield_DIRC_wt_cuts_sim ", "lp");
            leg_yield2->AddEntry(graph_yield_DIRC_wt_cuts_data, "graph_yield_DIRC_wt_cuts_data", "lp");
            
            mg_yield_3->Add(graph_yield_DIRC_wtc_cuts_sim);
            mg_yield_3->Add(graph_yield_DIRC_wtc_cuts_data);
            TLegend *leg_yield3 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_yield3->SetHeader("photon yield (proton)","C");
            leg_yield3->AddEntry(graph_yield_DIRC_wtc_cuts_sim, "graph_yield_DIRC_wtc_cuts_sim ", "lp");
            leg_yield3->AddEntry(graph_yield_DIRC_wtc_cuts_data, "graph_yield_DIRC_wtc_cuts_data", "lp");
            
            prt_canvasAdd("r_yield_1",800,400);
            mg_yield_1->SetTitle("photon yield ;#theta [degree]; count [#]");
            mg_yield_1->Draw("APL");
            mg_yield_1->GetHistogram()->GetYaxis()->SetRangeUser(0,120);
            leg_yield1->Draw();
            
            prt_canvasAdd("r_yield_2",800,400);
            mg_yield_2->SetTitle("photon yield ;#theta [degree]; count [#]");
            mg_yield_2->Draw("APL");
            mg_yield_2->GetHistogram()->GetYaxis()->SetRangeUser(0,120);
            leg_yield2->Draw();
            
            prt_canvasAdd("r_yield_3",800,400);
            mg_yield_3->SetTitle("photon yield ;#theta [degree]; count [#]");
            mg_yield_3->Draw("APL");
            mg_yield_3->GetHistogram()->GetYaxis()->SetRangeUser(0,120);
            leg_yield3->Draw();
            
            gPad->Modified();
            gPad->Update();
            
        }
        /////////////////////////////
        //// calc mom & distance/////
        /////////////////////////////
        if(true) {
            prt_canvasAdd("r_calc_mom",800,400);
            TMultiGraph *mg_calc_mom = new TMultiGraph();
            mg_calc_mom->Add(calc_mom);
            mg_calc_mom->SetTitle(" calculated momentum ;#theta [degree]; momentum [GeV/c]");
            mg_calc_mom->Draw("APL");
            mg_calc_mom->GetHistogram()->GetYaxis()->SetRangeUser(6800,7050);
            prt_canvasGet("r_calc_mom")->Update();
            TLine *line_mom = new TLine(0,0,0,1000);
            line_mom->SetY1(7000);
            line_mom->SetY2(7000);
            line_mom->SetX1(gPad->GetUxmin());
            line_mom->SetX2(gPad->GetUxmax());
            line_mom->SetLineColor(kRed);
            line_mom->Draw();
            
            prt_canvasAdd("r_calc_tof1tof2_distance",800,400);
            TMultiGraph *mg_calc_tof1tof2_distance = new TMultiGraph();
            mg_calc_tof1tof2_distance->Add(calc_tof1tof2_distance);
            mg_calc_tof1tof2_distance->SetTitle(" calculated distance ;#theta [degree]; tof2 tof1 distance [m]");
            mg_calc_tof1tof2_distance->Draw("APL");
            mg_calc_tof1tof2_distance->GetHistogram()->GetYaxis()->SetRangeUser(28,30);
            prt_canvasGet("r_calc_tof1tof2_distance")->Update();
            TLine *line_dis = new TLine(0,0,0,1000);
            line_dis->SetY1(28.507);
            line_dis->SetY2(28.507);
            line_dis->SetX1(gPad->GetUxmin());
            line_dis->SetX2(gPad->GetUxmax());
            line_dis->SetLineColor(kRed);
            line_dis->Draw();
        }
        // warning switch between p and pi
        gStyle->SetOptStat(0);
    }
    //////////////////////////
    // save/delete histogram//
    //////////////////////////
    
    std::cout<<"############"<< " no problem *** " <<std::endl;
    prt_canvasSave(2,0);
    prt_canvasDel("*");
    
    //for (int d=1; d<=14; d++) {
    //Printf("recoAngle%d_proton=%f", prtangle_vector[d], recoAngle[d]);
    //Printf("timeCut%d_proton=%f", prtangle_vector[d], timeCut[d]);
    //Printf("chAngleCut%d_proton=%f", prtangle_vector[d], chAngleCut[d]);
    //~ //Printf("if(prtangle == %d)gF1->SetParameter(1,%f);", prtangle_vector[d], recoAngle[d]);
    //}
    
    std::cout<<"############"<< " no problem **** " <<std::endl;
    // histo
    delete p_cherenkov_sim;
    delete p_cherenkov_data;
    delete p_cherenkov_sim_correction;
    delete p_cherenkov_data_correction;
    delete p_cherenkov_data_copy;
    delete p_cherenkov_sim_copy;
    delete p_cherenkov_mc_same_path;
    delete p_cherenkov_bg_sim;
    delete nph_sim;
    delete p_nph_sim;
    delete p_nph_good_sim;
    delete p_nph_true_sim;
    delete nph_data;
    delete p_nph_data;
    delete p_nph_good_data;
    delete p_diff_time_sim;
    delete p_diff_time_mctruth;
    delete p_diff_time_bg_sim;
    delete p_diff_time_data;
    delete p_photon_time_sim;
    delete p_photon_time_data;
    delete p_photon_time_sim_calc;
    delete p_photon_time_data_calc;
    /*
     //graphs
     delete graph_spr_sim;
     delete graph_spr_data_sub;
     delete graph_spr_sim_true;
     delete graph_spr_data_correction;
     delete graph_spr_sim_correction;
     delete graph_spr_data_org;
     delete graph_spr_data_cuts;
     delete graph_spr_sim_cuts;
     delete graph_spr_data_sub_cuts;
     delete graph_spr_sim_true_cuts;
     delete graph_diff_true;
     delete graph_diff_data;
     delete graph_diff_sim;
     delete graph_diff_data_mean;
     delete graph_diff_sim_mean;
     delete graph_diff_true_mean;
     delete graph_diff_true_sim;
     delete graph_cangle_sim;
     delete graph_cangle_sim_true;
     delete graph_cangle_data_correction;
     delete graph_cangle_sim_correction;
     delete graph_cangle_data_org;
     delete graph_yield_DIRC_wo_cuts_sim;
     delete graph_yield_DIRC_wt_cuts_sim;
     delete graph_yield_DIRC_wtc_cuts_sim;
     delete graph_yield_DIRC_true_sim;
     delete graph_yield_DIRC_wo_cuts_data;
     delete graph_yield_DIRC_wt_cuts_data;
     delete graph_yield_DIRC_wtc_cuts_data;
     
     delete calc_mom;
     delete calc_e_mom ;
     delete calc_e_tof1tof2_distance;
     */
    ffile_sim->Close();
    ffile_data->Close();
    
    std::cout<<"#################################################################"<<std::endl;
    std::cout<<"############ P cherenkov angle= "<< fAngleP <<std::endl;
    std::cout<<"############ Pi cherenkov angle= "<< fAnglePi <<std::endl;
    std::cout<<"############"<< " Macro End Succsessfully " <<std::endl;
    
   
}



//////////////////////////
// check file existance //
//////////////////////////

bool exists_test (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}


///////////////////////////////
// histo style calculation ////
///////////////////////////////
void HistoStyle(TH1F *p_cherenkov_sim, TH1F *p_cherenkov_data, TH1F *p_cherenkov_data_sub, TH1F *p_cherenkov_mc_same_path, TH1F *p_cherenkov_bg_sim, TH1F *nph_sim, TH1F *p_nph_sim,TH1F *p_nph_good_sim, TH1F *p_nph_true_sim, TH1F *p_diff_time_sim, TH1F*p_diff_time_data, TH1F*p_diff_time_data_sub, TH1F*p_diff_time_mctruth, TH1F* p_diff_time_bg_sim, TH1F *p_cherenkov_data_copy, TH1F *p_photon_time_sim, TH1F *p_photon_time_data, TH1F *p_photon_time_sim_calc, TH1F *p_photon_time_data_calc, TH1F *nph_data, TH1F *p_nph_data, TH1F *p_nph_good_data, TH1F *dac_hits_data, TH1F *dac_hits_sys_cus_data) {
    p_cherenkov_sim->SetName("MC signal");
    p_cherenkov_sim->SetLineColor(12);
    p_cherenkov_sim->SetLineStyle(1);
    p_cherenkov_sim->GetXaxis()->SetTitle("#theta_{C} [rad]");
    p_cherenkov_sim->GetYaxis()->SetTitle("entries [#]");
    p_cherenkov_sim->GetXaxis()->SetTitleSize(0.05);
    p_cherenkov_sim->GetYaxis()->SetTitleSize(0.05);
    p_cherenkov_sim->GetXaxis()->SetTitleOffset(0.9);
    p_cherenkov_sim->GetYaxis()->SetTitleOffset(1.0);
    p_cherenkov_sim->SetFillColor(12);
    p_cherenkov_sim->SetFillStyle(3003);
    
    p_cherenkov_data->SetLineColor(kBlue);
    p_cherenkov_data->SetLineStyle(1);
    p_cherenkov_data->GetXaxis()->SetTitle("#theta_{C} [rad]");
    p_cherenkov_data->GetYaxis()->SetTitle("entries [#]");
    p_cherenkov_data->GetXaxis()->SetTitleSize(0.05);
    p_cherenkov_data->GetYaxis()->SetTitleSize(0.05);
    p_cherenkov_data->GetXaxis()->SetTitleOffset(0.9);
    p_cherenkov_data->GetYaxis()->SetTitleOffset(1.0);
    p_cherenkov_data->SetFillColor(kBlue);
    p_cherenkov_data->SetFillStyle(3003);
    p_cherenkov_data->SetMarkerStyle(8);
    p_cherenkov_data->SetMarkerColor(kBlue);
    
    p_cherenkov_data_copy->SetLineColor(kBlue);
    p_cherenkov_data_copy->SetLineStyle(1);
    p_cherenkov_data_copy->GetXaxis()->SetTitle("#theta_{C} [rad]");
    p_cherenkov_data_copy->GetYaxis()->SetTitle("entries [#]");
    p_cherenkov_data_copy->GetXaxis()->SetTitleSize(0.05);
    p_cherenkov_data_copy->GetYaxis()->SetTitleSize(0.05);
    p_cherenkov_data_copy->GetXaxis()->SetTitleOffset(0.9);
    p_cherenkov_data_copy->GetYaxis()->SetTitleOffset(1.0);
    p_cherenkov_data_copy->SetFillColor(kBlue);
    p_cherenkov_data_copy->SetFillStyle(3003);
    p_cherenkov_data_copy->SetMarkerStyle(8);
    p_cherenkov_data_copy->SetMarkerColor(kBlue);
    
    p_cherenkov_data_sub->SetName("estimated data signal");
    p_cherenkov_data_sub->SetLineColor(kBlack);
    p_cherenkov_data_sub->SetLineStyle(1);
    p_cherenkov_data_sub->GetXaxis()->SetTitle("ç");
    p_cherenkov_data_sub->GetYaxis()->SetTitle("entries [#]");
    p_cherenkov_data_sub->GetXaxis()->SetTitleSize(0.05);
    p_cherenkov_data_sub->GetYaxis()->SetTitleSize(0.05);
    p_cherenkov_data_sub->GetXaxis()->SetTitleOffset(0.9);
    p_cherenkov_data_sub->GetYaxis()->SetTitleOffset(1.0);
    p_cherenkov_data_sub->SetFillColor(kMagenta);
    p_cherenkov_data_sub->SetFillStyle(3004);
    p_cherenkov_data_sub->SetMarkerStyle(8);
    p_cherenkov_data_sub->SetMarkerColor(kBlack);
    
    p_nph_good_sim->SetLineColor(kRed);
    p_nph_good_sim->SetLineStyle(1);
    p_nph_good_sim->GetXaxis()->SetTitle("number of photon per track [#]");
    p_nph_good_sim->GetYaxis()->SetTitle("entries [#]");
    p_nph_good_sim->GetXaxis()->SetTitleSize(0.05);
    p_nph_good_sim->GetYaxis()->SetTitleSize(0.05);
    p_nph_good_sim->GetXaxis()->SetTitleOffset(0.9);
    p_nph_good_sim->GetYaxis()->SetTitleOffset(1.0);
    //p_nph_good_sim->SetFillColor(kRed);
    //p_nph_good_sim->SetFillStyle(3001);
    
    
    p_nph_true_sim->SetLineColor(kRed);
    p_nph_true_sim->SetLineStyle(2);
    p_nph_true_sim->GetXaxis()->SetTitle("number of photon per track [#]");
    p_nph_true_sim->GetYaxis()->SetTitle("entries [#]");
    p_nph_true_sim->GetXaxis()->SetTitleSize(0.05);
    p_nph_true_sim->GetYaxis()->SetTitleSize(0.05);
    p_nph_true_sim->GetXaxis()->SetTitleOffset(0.9);
    p_nph_true_sim->GetYaxis()->SetTitleOffset(1.0);
    //p_nph_true_sim->SetFillColor(kRed);
    //p_nph_true_sim->SetFillStyle(3001);
    
    p_nph_sim->SetLineColor(kBlue);
    p_nph_sim->SetLineStyle(1);
    p_nph_sim->GetXaxis()->SetTitle("number of photon per track [#]");
    p_nph_sim->GetYaxis()->SetTitle("entries [#]");
    p_nph_sim->GetXaxis()->SetTitleSize(0.05);
    p_nph_sim->GetYaxis()->SetTitleSize(0.05);
    p_nph_sim->GetXaxis()->SetTitleOffset(0.9);
    p_nph_sim->GetYaxis()->SetTitleOffset(1.0);
    //p_nph_sim->SetFillColor(kBlue);
    //p_nph_sim->SetFillStyle(3003);
    
    nph_sim->SetLineColor(kBlack);
    nph_sim->SetLineStyle(1);
    nph_sim->GetXaxis()->SetTitle("number of photon per track [#]");
    nph_sim->GetYaxis()->SetTitle("entries [#]");
    nph_sim->GetXaxis()->SetTitleSize(0.05);
    nph_sim->GetYaxis()->SetTitleSize(0.05);
    nph_sim->GetXaxis()->SetTitleOffset(0.9);
    nph_sim->GetYaxis()->SetTitleOffset(1.0);
    //nph_sim->SetFillColor(kBlack);
    //nph_sim->SetFillStyle(3002);
    
    p_cherenkov_mc_same_path->SetLineColor(kMagenta);
    p_cherenkov_mc_same_path->SetLineStyle(1);
    p_cherenkov_mc_same_path->GetXaxis()->SetTitle("#theta_{C} [rad]");
    p_cherenkov_mc_same_path->GetYaxis()->SetTitle("entries [#]");
    p_cherenkov_mc_same_path->GetXaxis()->SetTitleSize(0.05);
    p_cherenkov_mc_same_path->GetYaxis()->SetTitleSize(0.05);
    p_cherenkov_mc_same_path->GetXaxis()->SetTitleOffset(0.9);
    p_cherenkov_mc_same_path->GetYaxis()->SetTitleOffset(1.0);
    p_cherenkov_mc_same_path->SetFillColor(kMagenta);
    p_cherenkov_mc_same_path->SetFillStyle(3001);
    
    p_cherenkov_bg_sim->SetLineColor(kBlack);
    p_cherenkov_bg_sim->SetLineStyle(1);
    p_cherenkov_bg_sim->GetXaxis()->SetTitle("#theta_{C} [rad]");
    p_cherenkov_bg_sim->GetYaxis()->SetTitle("entries [#]");
    p_cherenkov_bg_sim->GetXaxis()->SetTitleSize(0.05);
    p_cherenkov_bg_sim->GetYaxis()->SetTitleSize(0.05);
    p_cherenkov_bg_sim->GetXaxis()->SetTitleOffset(0.9);
    p_cherenkov_bg_sim->GetYaxis()->SetTitleOffset(1.0);
    p_cherenkov_bg_sim->SetFillColor(kBlack);
    p_cherenkov_bg_sim->SetFillStyle(4050);
    p_cherenkov_bg_sim->SetFillStyle(3001);
    p_cherenkov_bg_sim->SetMarkerStyle(4);
    
    p_diff_time_sim->SetName("MC signal");
    p_diff_time_sim->SetLineColor(12);
    p_diff_time_sim->SetLineStyle(1);
    p_diff_time_sim->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
    p_diff_time_sim->GetYaxis()->SetTitle("entries [#]");
    p_diff_time_sim->GetXaxis()->SetTitleSize(0.05);
    p_diff_time_sim->GetYaxis()->SetTitleSize(0.05);
    p_diff_time_sim->GetXaxis()->SetTitleOffset(0.9);
    p_diff_time_sim->GetYaxis()->SetTitleOffset(1.0);
    p_diff_time_sim->SetFillColor(12);
    p_diff_time_sim->SetFillStyle(3003);
    
    p_diff_time_data->SetLineColor(kBlue);
    p_diff_time_data->SetLineStyle(1);
    p_diff_time_data->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
    p_diff_time_data->GetYaxis()->SetTitle("entries [#]");
    p_diff_time_data->GetXaxis()->SetTitleSize(0.05);
    p_diff_time_data->GetYaxis()->SetTitleSize(0.05);
    p_diff_time_data->GetXaxis()->SetTitleOffset(0.9);
    p_diff_time_data->GetYaxis()->SetTitleOffset(1.0);
    p_diff_time_data->SetFillColor(kBlue);
    p_diff_time_data->SetFillStyle(3003);
    p_diff_time_data->SetMarkerStyle(8);
    p_diff_time_data->SetMarkerSize(0.5);
    p_diff_time_data->SetMarkerColor(kBlue);
    
    p_diff_time_data_sub->SetName("estimated data signal");
    p_diff_time_data_sub->SetLineColor(kRed);
    p_diff_time_data_sub->SetLineStyle(1);
    p_diff_time_data_sub->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
    p_diff_time_data_sub->GetYaxis()->SetTitle("entries [#]");
    p_diff_time_data_sub->GetXaxis()->SetTitleSize(0.05);
    p_diff_time_data_sub->GetYaxis()->SetTitleSize(0.05);
    p_diff_time_data_sub->GetXaxis()->SetTitleOffset(0.9);
    p_diff_time_data_sub->GetYaxis()->SetTitleOffset(1.0);
    p_diff_time_data_sub->SetFillColor(kRed);
    p_diff_time_data_sub->SetFillStyle(3004);
    p_diff_time_data_sub->SetMarkerStyle(8);
    p_diff_time_data_sub->SetMarkerColor(kRed);
    
    p_diff_time_mctruth->SetLineColor(kMagenta);
    p_diff_time_mctruth->SetLineStyle(1);
    p_diff_time_mctruth->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
    p_diff_time_mctruth->GetYaxis()->SetTitle("entries [#]");
    p_diff_time_mctruth->GetXaxis()->SetTitleSize(0.05);
    p_diff_time_mctruth->GetYaxis()->SetTitleSize(0.05);
    p_diff_time_mctruth->GetXaxis()->SetTitleOffset(0.9);
    p_diff_time_mctruth->GetYaxis()->SetTitleOffset(1.0);
    p_diff_time_mctruth->SetFillColor(kMagenta);
    p_diff_time_mctruth->SetFillStyle(3001);
    
    p_diff_time_bg_sim->SetLineColor(kBlack);
    p_diff_time_bg_sim->SetLineStyle(1);
    p_diff_time_bg_sim->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
    p_diff_time_bg_sim->GetYaxis()->SetTitle("entries [#]");
    p_diff_time_bg_sim->GetXaxis()->SetTitleSize(0.05);
    p_diff_time_bg_sim->GetYaxis()->SetTitleSize(0.05);
    p_diff_time_bg_sim->GetXaxis()->SetTitleOffset(0.9);
    p_diff_time_bg_sim->GetYaxis()->SetTitleOffset(1.0);
    p_diff_time_bg_sim->SetFillColor(kBlack);
    p_diff_time_bg_sim->SetFillStyle(3001);
    p_diff_time_bg_sim->SetMarkerStyle(4);
    
    p_photon_time_sim->SetLineColor(kBlue);
    p_photon_time_sim->SetLineStyle(1);
    p_photon_time_sim->GetXaxis()->SetTitle("t_{measured} [ns]");
    p_photon_time_sim->GetYaxis()->SetTitle("entries [#]");
    p_photon_time_sim->GetXaxis()->SetTitleSize(0.05);
    p_photon_time_sim->GetYaxis()->SetTitleSize(0.05);
    p_photon_time_sim->GetXaxis()->SetTitleOffset(0.9);
    p_photon_time_sim->GetYaxis()->SetTitleOffset(1.0);
    p_photon_time_sim->SetFillColor(kBlue);
    p_photon_time_sim->SetFillStyle(3003);
    
    p_photon_time_data->SetLineColor(kRed);
    p_photon_time_data->SetLineStyle(1);
    p_photon_time_data->GetXaxis()->SetTitle("t_{measured} [ns]");
    p_photon_time_data->GetYaxis()->SetTitle("entries [#]");
    p_photon_time_data->GetXaxis()->SetTitleSize(0.05);
    p_photon_time_data->GetYaxis()->SetTitleSize(0.05);
    p_photon_time_data->GetXaxis()->SetTitleOffset(0.9);
    p_photon_time_data->GetYaxis()->SetTitleOffset(1.0);
    p_photon_time_data->SetFillColor(kRed);
    p_photon_time_data->SetFillStyle(3001);
    p_photon_time_data->SetMarkerStyle(21);
    p_photon_time_data->SetMarkerColor(kRed);
    p_photon_time_data->SetMarkerSize(0.5);
    
    p_photon_time_sim_calc->SetLineColor(kBlue);
    p_photon_time_sim_calc->SetLineStyle(1);
    p_photon_time_sim_calc->GetXaxis()->SetTitle("t_{calc} [ns]");
    p_photon_time_sim_calc->GetYaxis()->SetTitle("entries [#]");
    p_photon_time_sim_calc->GetXaxis()->SetTitleSize(0.05);
    p_photon_time_sim_calc->GetYaxis()->SetTitleSize(0.05);
    p_photon_time_sim_calc->GetXaxis()->SetTitleOffset(0.9);
    p_photon_time_sim_calc->GetYaxis()->SetTitleOffset(1.0);
    p_photon_time_sim_calc->SetFillColor(kBlue);
    p_photon_time_sim_calc->SetFillStyle(3003);
    
    p_photon_time_data_calc->SetLineColor(kRed);
    p_photon_time_data_calc->SetLineStyle(1);
    p_photon_time_data_calc->GetXaxis()->SetTitle("t_{calc} [ns]");
    p_photon_time_data_calc->GetYaxis()->SetTitle("entries [#]");
    p_photon_time_data_calc->GetXaxis()->SetTitleSize(0.05);
    p_photon_time_data_calc->GetYaxis()->SetTitleSize(0.05);
    p_photon_time_data_calc->GetXaxis()->SetTitleOffset(0.9);
    p_photon_time_data_calc->GetYaxis()->SetTitleOffset(1.0);
    p_photon_time_data_calc->SetFillColor(kRed);
    p_photon_time_data_calc->SetFillStyle(3001);
    p_photon_time_data_calc->SetMarkerStyle(21);
    p_photon_time_data_calc->SetMarkerColor(kRed);
    p_photon_time_data_calc->SetMarkerSize(0.5);
    
    p_nph_good_data->SetLineColor(kRed);
    p_nph_good_data->SetLineStyle(1);
    p_nph_good_data->GetXaxis()->SetTitle("number of photon per track [#]");
    p_nph_good_data->GetYaxis()->SetTitle("entries [#]");
    p_nph_good_data->GetXaxis()->SetTitleSize(0.05);
    p_nph_good_data->GetYaxis()->SetTitleSize(0.05);
    p_nph_good_data->GetXaxis()->SetTitleOffset(0.9);
    p_nph_good_data->GetYaxis()->SetTitleOffset(1.0);
    //p_nph_good_data->SetFillColor(kRed);
    //p_nph_good_data->SetFillStyle(3001);
    
    p_nph_data->SetLineColor(kBlue);
    p_nph_data->SetLineStyle(1);
    p_nph_data->GetXaxis()->SetTitle("number of photon per track [#]");
    p_nph_data->GetYaxis()->SetTitle("entries [#]");
    p_nph_data->GetXaxis()->SetTitleSize(0.05);
    p_nph_data->GetYaxis()->SetTitleSize(0.05);
    p_nph_data->GetXaxis()->SetTitleOffset(0.9);
    p_nph_data->GetYaxis()->SetTitleOffset(1.0);
    //p_nph_data->SetFillColor(kBlue);
    //p_nph_data->SetFillStyle(3003);
    
    nph_data->SetLineColor(kBlack);
    nph_data->SetLineStyle(1);
    nph_data->GetXaxis()->SetTitle("number of photon per track [#]");
    nph_data->GetYaxis()->SetTitle("entries [#]");
    nph_data->GetXaxis()->SetTitleSize(0.05);
    nph_data->GetYaxis()->SetTitleSize(0.05);
    nph_data->GetXaxis()->SetTitleOffset(0.9);
    nph_data->GetYaxis()->SetTitleOffset(1.0);
    //nph_data->SetFillColor(kBlack);
    //nph_data->SetFillStyle(3002);
    
    /*
     dac_hits_data->SetLineColor(kGreen);
     dac_hits_data->SetLineStyle(1);
     dac_hits_data->GetXaxis()->SetTitle(" per track [#]");
     dac_hits_data->GetYaxis()->SetTitle("entries [#]");
     dac_hits_data->GetXaxis()->SetTitleSize(0.05);
     dac_hits_data->GetYaxis()->SetTitleSize(0.05);
     dac_hits_data->GetXaxis()->SetTitleOffset(0.9);
     dac_hits_data->GetYaxis()->SetTitleOffset(1.0);
     dac_hits_data->SetFillColor(kGreen);
     dac_hits_data->SetFillStyle(3002);
     
     
     dac_hits_sys_cus_data->SetLineColor(kGray);
     dac_hits_sys_cus_data->SetLineStyle(1);
     dac_hits_sys_cus_data->GetXaxis()->SetTitle(" per track [#]");
     dac_hits_sys_cus_data->GetYaxis()->SetTitle("entries [#]");
     dac_hits_sys_cus_data->GetXaxis()->SetTitleSize(0.05);
     dac_hits_sys_cus_data->GetYaxis()->SetTitleSize(0.05);
     dac_hits_sys_cus_data->GetXaxis()->SetTitleOffset(0.9);
     dac_hits_sys_cus_data->GetYaxis()->SetTitleOffset(1.0);
     dac_hits_sys_cus_data->SetFillColor(kGray);
     dac_hits_sys_cus_data->SetFillStyle(3002);
     */
    
}

void DiffNorm(TH1F *p_diff_time_sim, TH1F *p_diff_time_data) {
    
    /////////////////////////////
    // time diff normalization //
    /////////////////////////////
    // Method 1
    //     TAxis *axis_diff_right = p_diff_time_sim->GetXaxis();
    //     double xmin_diff_right = 2;
    //     double xmax_diff_right = 10;
    //     int bmin_diff_right = axis_diff_right->FindBin(xmin_diff_right);
    //     int bmax_diff_right = axis_diff_right->FindBin(xmax_diff_right);
    //     double integral_diff_right = p_diff_time_sim->Integral(bmin_diff_right,bmax_diff_right);
    //    //
    //     TAxis *axis_diff_data_right = p_diff_time_data->GetXaxis();
    //     double xmin_diff_data_right = 2;
    //     double xmax_diff_data_right = 10;
    //     int bmin_diff_data_right = axis_diff_data_right->FindBin(xmin_diff_data_right);
    //     int bmax_diff_data_right = axis_diff_data_right->FindBin(xmax_diff_data_right);
    //     double integral_data_diff_right = p_diff_time_data->Integral(bmin_diff_data_right,bmax_diff_data_right);
    //    //
    //     TAxis *axis_diff_lift = p_diff_time_sim->GetXaxis();
    //     double xmin_diff_lift = -2;
    //     double xmax_diff_lift = -10;
    //     int bmin_diff_lift = axis_diff_lift->FindBin(xmin_diff_lift);
    //     int bmax_diff_lift = axis_diff_lift->FindBin(xmax_diff_lift);
    //     double integral_diff_lift = p_diff_time_sim->Integral(bmin_diff_lift,bmax_diff_lift);
    //    //
    //     TAxis *axis_diff_data_lift = p_diff_time_data->GetXaxis();
    //     double xmin_diff_data_lift = -2;
    //     double xmax_diff_data_lift = -10;
    //     int bmin_diff_data_lift = axis_diff_data_lift->FindBin(xmin_diff_data_lift);
    //     int bmax_diff_data_lift = axis_diff_data_lift->FindBin(xmax_diff_data_lift);
    //     double integral_data_diff_lift = p_diff_time_data->Integral(bmin_diff_data_lift,bmax_diff_data_lift);
    //    //
    //     Double_t norm_diff= (integral_diff_right+integral_diff_lift)/(integral_data_diff_right+integral_data_diff_lift);
    //    //
    //     std::cout<<"############  norm_diff= "<< norm_diff <<std::endl;
    //     p_diff_time_data->Scale(norm_diff);
    
    
    // Method 2
    p_diff_time_data->Scale(p_diff_time_sim->GetMaximum() /p_diff_time_data->GetMaximum());
    
    
    // // Method 3
    //
    // Double_t norm1 = 1;
    // Double_t norm2 = 1;
    // Double_t scale1 = norm1/( p_diff_time_sim->Integral());
    // p_diff_time_sim->Scale(scale1);
    // Double_t scale2 = norm2/(p_diff_time_data->Integral());
    // p_diff_time_data->Scale(scale2);
    
    // // Method 4
    // Double_t norm1 = p_diff_time_sim->GetEntries();
    // p_diff_time_sim->Scale(1/norm1);
    // Double_t norm2 = p_diff_time_data->GetEntries();
    // p_diff_time_data->Scale(1/norm2);
    //
}

void HistoStyleMatch(TH1F *p_cherenkov_sim, TH1F *p_cherenkov_data) {
    p_cherenkov_sim->SetName("MC signal");
    p_cherenkov_sim->SetLineColor(kMagenta);
    p_cherenkov_sim->SetLineStyle(1);
    p_cherenkov_sim->GetXaxis()->SetTitle("#theta_{C} [rad]");
    p_cherenkov_sim->GetYaxis()->SetTitle("entries [#]");
    p_cherenkov_sim->GetXaxis()->SetTitleSize(0.05);
    p_cherenkov_sim->GetYaxis()->SetTitleSize(0.05);
    p_cherenkov_sim->GetXaxis()->SetTitleOffset(0.9);
    p_cherenkov_sim->GetYaxis()->SetTitleOffset(1.0);
    p_cherenkov_sim->SetFillColor(kMagenta);
    p_cherenkov_sim->SetFillStyle(3001);
    
    p_cherenkov_data->SetLineColor(kBlack);
    p_cherenkov_data->SetLineStyle(1);
    p_cherenkov_data->GetXaxis()->SetTitle("#theta_{C} [rad]");
    p_cherenkov_data->GetYaxis()->SetTitle("entries [#]");
    p_cherenkov_data->GetXaxis()->SetTitleSize(0.05);
    p_cherenkov_data->GetYaxis()->SetTitleSize(0.05);
    p_cherenkov_data->GetXaxis()->SetTitleOffset(0.9);
    p_cherenkov_data->GetYaxis()->SetTitleOffset(1.0);
    p_cherenkov_data->SetFillColor(kBlack);
    p_cherenkov_data->SetFillStyle(3003);
    p_cherenkov_data->SetMarkerStyle(8);
    p_cherenkov_data->SetMarkerColor(kBlack);
}

Double_t*  YieldGausFit(TH1F *photonYieldHistogram) {
    const int var=2;
    static Double_t return_yield_array[var];
    
    TF1 *PhotonYieldFit = new TF1("PhotonYieldFit","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,200);
    PhotonYieldFit->SetParameters(100,100,20);
    PhotonYieldFit->SetParNames("p0","mean","#sigma");
    PhotonYieldFit->SetParLimits(0,0.1,1E6);
    PhotonYieldFit->SetParLimits(1,10,300);
    PhotonYieldFit->SetParLimits(2,5,100);
    photonYieldHistogram->Fit("PhotonYieldFit","M","",0,200);
    Double_t mean_yield = PhotonYieldFit->GetParameter(1);
    Double_t sigma_yield = PhotonYieldFit->GetParameter(2);
    return_yield_array[0]=mean_yield;
    return_yield_array[1]=sigma_yield;
    
    return return_yield_array;
}


Double_t  YieldGausFit_double(TH1F *photonYieldHistogram) {
    
    Double_t return_yield=0;
    
    TF1 *PhotonYieldFit = new TF1("PhotonYieldFit","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,200);
    PhotonYieldFit->SetParameters(100,100,20);
    PhotonYieldFit->SetParNames("p0","mean","#sigma");
    PhotonYieldFit->SetParLimits(0,0.1,1E6);
    PhotonYieldFit->SetParLimits(1,10,300);
    PhotonYieldFit->SetParLimits(2,5,100);
    photonYieldHistogram->Fit("PhotonYieldFit","M","",0,200);
    Double_t mean_yield = PhotonYieldFit->GetParameter(1);
    Double_t sigma_yield = PhotonYieldFit->GetParameter(2);
    
    return_yield=mean_yield;
    
    return return_yield;
}

void HistoStyle_3colors(TH1F *histo1, TH1F *histo2, TH1F *histo3){
    
    histo1->SetLineColor(kRed);
    histo1->SetLineStyle(1);
    histo1->GetXaxis()->SetTitle("number of photon per track [#]");
    histo1->GetYaxis()->SetTitle("entries [#]");
    histo1->GetXaxis()->SetTitleSize(0.05);
    histo1->GetYaxis()->SetTitleSize(0.05);
    histo1->GetXaxis()->SetTitleOffset(0.9);
    histo1->GetYaxis()->SetTitleOffset(1.0);
    //histo1->SetFillColor(kRed);
    //histo1->SetFillStyle(3001);
    
    histo2->SetLineColor(kBlue);
    histo2->SetLineStyle(1);
    histo2->GetXaxis()->SetTitle("number of photon per track [#]");
    histo2->GetYaxis()->SetTitle("entries [#]");
    histo2->GetXaxis()->SetTitleSize(0.05);
    histo2->GetYaxis()->SetTitleSize(0.05);
    histo2->GetXaxis()->SetTitleOffset(0.9);
    histo2->GetYaxis()->SetTitleOffset(1.0);
    //histo2->SetFillColor(kBlue);
    //histo2->SetFillStyle(3003);
    
    histo3->SetLineColor(kBlack);
    histo3->SetLineStyle(1);
    histo3->GetXaxis()->SetTitle("number of photon per track [#]");
    histo3->GetYaxis()->SetTitle("entries [#]");
    histo3->GetXaxis()->SetTitleSize(0.05);
    histo3->GetYaxis()->SetTitleSize(0.05);
    histo3->GetXaxis()->SetTitleOffset(0.9);
    histo3->GetYaxis()->SetTitleOffset(1.0);
    //histo3->SetFillColor(kBlack);
    //histo3->SetFillStyle(3002);
}



