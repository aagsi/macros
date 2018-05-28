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
#include "THStack.h"
#include "TLatex.h"
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
#include <tuple>

TCanvas *c1 = new TCanvas("c1","c1",300,200);
////////////////
// constants  //
////////////////
Double_t beam_momentum= 7000; // MeV
Double_t m_proton = 938.28;
Double_t m_pi = 139.570;
Double_t nano_value = 0.000000001; // 10e-9
Double_t c = 299792458; // speed of light
//Double_t measured_d_tof2tof1_plot1_co= 28.507; //old
Double_t measured_d_tof2tof1_plot1_co= 28.428;

// calculation of protons and pions velocity based on the momentum
Double_t v_proton = c* sqrt(1- (m_proton*m_proton) /(beam_momentum * beam_momentum+ m_proton*m_proton));
Double_t v_pi = c* sqrt(1-( m_pi*m_pi) /(beam_momentum * beam_momentum+ m_pi*m_pi));
Double_t momentum=7.0;
Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
Double_t fAngleP = acos(sqrt(momentum*momentum+ mass[4]*mass[4])/momentum/1.4738)-0.00;
Double_t fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00;
Double_t  prtangle;

/////////////////
// proto types //
/////////////////
std::pair<Double_t, Double_t> FitHisto_m(TH1F *hiso_m_fit_option, Double_t cangle_true, Double_t fit_range);
std::pair<Double_t, Double_t> FitHisto_0(TH1F *hiso_0_fit_option, Double_t cangle_true, Double_t fit_range);
// file existance
bool exists_test (const std::string& name);
// histo style
void histo_style_match(TH1F *p_cherenkov_sim_org, TH1F *p_cherenkov_data_org);
void histo_style_photon_yield(TH1F *histo1, TH1F *histo2, TH1F *histo3);
void histo_style_cherenkov(TH1F *histo1, TH1F *histo2, TH1F *histo3, TH1F *histo4, TH1F *histo5);
void histo_style_time_diff(TH1F *histo1, TH1F *histo2, TH1F *histo3, TH1F *histo4, TH1F *histo5);
void histo_style_photon_time(TH1F *histo1, TH1F *histo2, TH1F *histo3, TH1F *histo4);

// diff norm
void DiffNorm(TH1F *p_diff_time_sim, TH1F *p_diff_time_data);

////////////////
// error Stuff//
////////////////
TFile *ffile_error_polar_sim_p, *ffile_error_stat_sim_p, *ffile_error_fit_range_sim_p;
TH1F *error_polar_p_cherenkov_sim_org,*error_polar_p_cherenkov_sim_corrected ,*error_polar_p_cherenkov_mc_same_path,*error_polar_yield_wo_sim,*error_polar_yield_wt_sim,*error_polar_yield_wtc_sim;
TH1F *error_stat_p_cherenkov_sim_org,*error_stat_p_cherenkov_sim_corrected ,*error_stat_p_cherenkov_mc_same_path,*error_stat_yield_wo_sim,*error_stat_yield_wt_sim,*error_stat_yield_wtc_sim;

TH1F *error_fit_range_p_cherenkov_sim_org,*error_fit_range_p_cherenkov_sim_corrected ,*error_fit_range_p_cherenkov_mc_same_path,*error_fit_range_yield_wo_sim,*error_fit_range_yield_wt_sim,*error_fit_range_yield_wtc_sim;

TH1F *histo_distribution_error_polar_spr_p_cherenkov_sim_org[14], *histo_distribution_error_stat_spr_p_cherenkov_sim_org[14], *histo_distribution_error_fit_range_spr_p_cherenkov_sim_org[14];
TH1F *histo_distribution_error_polar_spr_p_cherenkov_sim_corrected[14], *histo_distribution_error_stat_spr_p_cherenkov_sim_corrected[14], *histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected[14];
//cangle error
TH1F *histo_distribution_error_polar_cangle_p_cherenkov_sim_org[14], *histo_distribution_error_stat_cangle_p_cherenkov_sim_org[14], *histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org[14];
TH1F *histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected[14], *histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected[14], *histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected[14];

TH1F *histo_distribution_error_polar_p_yield_sim_wo[14], *histo_distribution_error_stat_p_yield_sim_wo[14], *histo_distribution_error_fit_range_p_yield_sim_wo[14];
TH1F *histo_distribution_error_polar_p_yield_sim_wt[14], *histo_distribution_error_stat_p_yield_sim_wt[14], *histo_distribution_error_fit_range_p_yield_sim_wt[14];
TH1F *histo_distribution_error_polar_p_yield_sim_wtc[14], *histo_distribution_error_stat_p_yield_sim_wtc[14], *histo_distribution_error_fit_range_p_yield_sim_wtc[14];

////////////////////
// function   //////
////////////////////
Bool_t bool_error=true;
Bool_t bool_error_histo=true;
Bool_t bool_partII=true;
Bool_t bool_partII_histo=true;
Bool_t bool_partIII=false;
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
    Double_t error_all_y_spr_sim_org[n]={0};
    Double_t error_all_y_spr_sim_corrected[n]={0};
    Double_t error_all_y_yield_wo_sim_corrected[n]={0};
    Double_t error_all_y_yield_wt_sim_corrected[n]={0};
    Double_t error_all_y_yield_wtc_sim_corrected[n]={0};
    Double_t error_all_y_cangle_sim_corrected[n]={0};
    Double_t error_all_y_cangle_sim_org[n]={0};
    
    Double_t x[n], y_spr_data_sub[n], y_spr_sim_true[n], y_spr_sim_org[n],y_spr_sim_corrected[n],  y_diff_true[n], y_diff_data[n], y_mean_diff_data[n], y_mean_diff_sim[n], y_mean_diff_true[n], y_diff_sim[n];
    Double_t y_cangle_data_sub[n], y_cangle_sim_org[n];
    Double_t y_cangle_sim_true[n];
    Double_t y_yield_wo_sim[n], y_yield_wt_sim[n], y_yield_wtc_sim[n],  y_yield_true_sim[n], y_yield_wo_data[n], y_yield_wt_data[n], y_yield_wtc_data[n], y_yield_dac_hits_sys_cus_data[n], y_cangle_data_corrected[n], y_spr_data_org[n], y_cangle_data_org[n], y_spr_data_corrected[n], y_spr_data_corrected_error[n], y_cangle_sim_corrected[n];
    
    int counter =0 ;
    prt_savepath="proton";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    TFile *ffile_sim, *ffile_data;
    TH1F *tof_pid, *p_cherenkov_sim_org, *p_cherenkov_data_org, *p_cherenkov_data_copy, *p_cherenkov_sim_copy, *p_cherenkov_data_sub, *p_cherenkov_mc_same_path, *p_cherenkov_bg_sim, *p_cherenkov_data_corrected, *p_cherenkov_sim_corrected, *p_yield_wo_sim, *p_yield_wt_sim, *p_yield_wtc_sim, *p_yield_true_sim;
    TH1F *p_diff_time_sim, *p_diff_time_data, *p_diff_time_data_sub, *p_diff_time_mctruth, * p_diff_time_bg_sim, *p_photon_time_sim, *p_photon_time_data,*p_photon_time_sim_calc, *p_photon_time_data_calc;
    TH1F *p_yield_wo_data, *p_yield_wt_data, *p_yield_wtc_data, *dac_hits_data, *dac_hits_sys_cus_data;
    TH1F *hist_ambiguity_lut_sim, *histo_photon_ambiguity_wo_sim, *histo_photon_ambiguity_wt_sim, *histo_photon_ambiguity_wtc_sim;
    TH1F *hist_ambiguity_lut_data, *histo_photon_ambiguity_wo_data, *histo_photon_ambiguity_wt_data, *histo_photon_ambiguity_wtc_data;
    THStack *hs, *hs2, *hs3, *hs4, *hs5, *hs6, *hs7, *hs8, *hs9, *hs10;
    std::cout<<"############"<< " no problem 0 " <<std::endl;
    TGraph *calc_mom = new TGraph();
    TGraph *calc_e_mom = new TGraph();
    TGraph *calc_tof1tof2_distance = new TGraph();
    TGraph *calc_e_tof1tof2_distance = new TGraph();
    //gStyle->SetOptFit(1);
    //gStyle->SetOptStat(0);
    TGraph * graph_spr_data_sub,* graph_spr_sim_true,* graph_diff_true,* graph_diff_data,* graph_diff_sim,* graph_diff_data_mean,* graph_diff_sim_mean,* graph_diff_true_mean,* graph_diff_true_sim,* graph_cangle_data_sub,* graph_cangle_sim_org,* graph_cangle_sim_true,* graph_cangle_sim_corrected;
    ///////////////////////////
    // Error Calculations   ///
    ///////////////////////////
    if(bool_error){
        ///////////////////////////
        // Error fit range      ///
        ///////////////////////////
        for (Int_t e=0; e<=14; e++) {
            // SPR
            histo_distribution_error_fit_range_spr_p_cherenkov_sim_org[e] = new TH1F(Form("histo_distribution_error_fit_range_spr_p_cherenkov_sim_org_%d",e),Form("histo_distribution_error_fit_range_spr_p_cherenkov_sim_org_%d; SPR [mrad];entries [#]",e), 50,6,14);
            histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected[e] = new TH1F(Form("histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected_%d",e),Form("histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected_%d; SPR [mrad];entries [#]",e), 50,6,14);
            //Cangle
            histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org[e] = new TH1F(Form("histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org_%d",e),Form("histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org_%d; cangle [mrad];entries [#]",e), 100,0.8,0.85);
            histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected[e] = new TH1F(Form("histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected_%d",e),Form("histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected_%d; cangle [mrad];entries [#]",e),100,0.8,0.85);
        }
        Int_t counter_error_fit_range=0;
        for (int i=20; i<=150; i+=10) {
            prtangle= i;
            TString nid = Form("_%2.0d", i);
            for (int j=1; j<=14; j+=1) { // changable
                TString calc_error_fit_range_sim_p_path = Form("/Users/ahmed/perforamnce/spr_data_sim/spr_wt_%d_sph_p_sim_spr.root", i);
                cout<<"sim path for error calc p = " <<calc_error_fit_range_sim_p_path<<endl;
                string path_calc_error_fit_range_sim_p = (string)calc_error_fit_range_sim_p_path;
                cout<<"exists_test(path_calc_error_fit_range_sim_p)" <<exists_test(path_calc_error_fit_range_sim_p)<<endl;
                if (!exists_test(path_calc_error_fit_range_sim_p)) continue;
                ffile_error_fit_range_sim_p = TFile::Open(calc_error_fit_range_sim_p_path,"read");
                error_fit_range_p_cherenkov_sim_org=(TH1F*)ffile_error_fit_range_sim_p->Get("fHist");
                error_fit_range_p_cherenkov_sim_corrected=(TH1F*)ffile_error_fit_range_sim_p->Get("fHist_correction");
                error_fit_range_p_cherenkov_mc_same_path=(TH1F*)ffile_error_fit_range_sim_p->Get("fHist_same_path");
                Double_t error_fit_range_cangle_MC_true =  error_fit_range_p_cherenkov_mc_same_path->GetXaxis()->GetBinCenter(error_fit_range_p_cherenkov_mc_same_path->GetMaximumBin());
                Double_t fit_range_value= j /100.0; // changable
                auto error_fit_range_resultPair_org= FitHisto_0( error_fit_range_p_cherenkov_sim_org , error_fit_range_cangle_MC_true, fit_range_value);
                Double_t error_fit_range_cangle_error_fit_range_p_cherenkov_sim_org = error_fit_range_resultPair_org.first;
                Double_t error_fit_range_spr_error_fit_range_p_cherenkov_sim_org = error_fit_range_resultPair_org.second;
                auto error_fit_range_resultPair_corrected= FitHisto_0( error_fit_range_p_cherenkov_sim_corrected , error_fit_range_cangle_MC_true, fit_range_value);
                Double_t error_fit_range_cangle_error_fit_range_p_cherenkov_sim_corrected = error_fit_range_resultPair_corrected.first;
                Double_t error_fit_range_spr_error_fit_range_p_cherenkov_sim_corrected = error_fit_range_resultPair_corrected.second;
                //cout<<"counter_error_fit_range = " <<counter_error_fit_range<<endl;
                histo_distribution_error_fit_range_spr_p_cherenkov_sim_org[counter_error_fit_range]->Fill(error_fit_range_spr_error_fit_range_p_cherenkov_sim_org*1000.0);
                histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected[counter_error_fit_range]->Fill(error_fit_range_spr_error_fit_range_p_cherenkov_sim_corrected*1000.0);
                
                histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org[counter_error_fit_range]->Fill(error_fit_range_cangle_error_fit_range_p_cherenkov_sim_org);
                histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected[counter_error_fit_range]->Fill(error_fit_range_cangle_error_fit_range_p_cherenkov_sim_corrected);
                
                delete error_fit_range_p_cherenkov_sim_org;
                delete error_fit_range_p_cherenkov_mc_same_path;
                delete error_fit_range_p_cherenkov_sim_corrected;
                ///////////////////
                // Close files  ///
                ///////////////////
                delete ffile_error_fit_range_sim_p;
            }
            ++counter_error_fit_range;
        }
        //    for (Int_t e=0; e<=14; e++) {
        //
        //        TString ko = Form("_%d", e);
        //        prt_canvasAdd("r_test_fit_range_org"+ko,800,400);
        //        histo_distribution_error_fit_range_spr_p_cherenkov_sim_org[e]->Draw();
        //        prt_canvasAdd("r_test_fit_range_corrected"+ko,800,400);
        //        histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected[e]->Draw();
        //    }
        
        ///////////////////////////
        // Error statistics     ///
        ///////////////////////////
        for (Int_t e=0; e<=14; e++) {
            histo_distribution_error_stat_spr_p_cherenkov_sim_org[e] = new TH1F(Form("histo_distribution_error_stat_spr_p_cherenkov_sim_org_%d",e),Form("histo_distribution_error_stat_spr_p_cherenkov_sim_org_%d; SPR [mrad];entries [#]",e), 100,6,14);
            histo_distribution_error_stat_spr_p_cherenkov_sim_corrected[e] = new TH1F(Form("histo_distribution_error_stat_spr_p_cherenkov_sim_corrected_%d",e),Form("histo_distribution_error_stat_spr_p_cherenkov_sim_corrected_%d; SPR [mrad];entries [#]",e), 100,6,14);
            //Cangle
            histo_distribution_error_stat_cangle_p_cherenkov_sim_org[e] = new TH1F(Form("histo_distribution_error_stat_cangle_p_cherenkov_sim_org_%d",e),Form("histo_distribution_error_stat_cangle_p_cherenkov_sim_org_%d; cangle [mrad];entries [#]",e), 100,0.8,0.85);
            histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected[e] = new TH1F(Form("histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected_%d",e),Form("histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected_%d; cangle [mrad];entries [#]",e), 100,0.8,0.85);
            // photon yield
            histo_distribution_error_stat_p_yield_sim_wtc[e] = new TH1F(Form("histo_distribution_error_stat_p_yield_sim_wtc_%d",e),Form("histo_distribution_error_stat_p_yield_sim_wtc_%d; SPR [mrad];entries [#]",e), 500,0,100);
            histo_distribution_error_stat_p_yield_sim_wt[e] = new TH1F(Form("histo_distribution_error_stat_p_yield_sim_wt_%d",e),Form("histo_distribution_error_stat_p_yield_sim_wt_%d; SPR [mrad];entries [#]",e), 500,0,100);
            histo_distribution_error_stat_p_yield_sim_wo[e] = new TH1F(Form("histo_distribution_error_stat_p_yield_sim_wo_%d",e),Form("histo_distribution_error_stat_p_yield_sim_wo_%d; SPR [mrad];entries [#]",e), 500,0,100);
            
        }
        Int_t counter_error_stat=0;
        for (int i=20; i<=150; i+=10) {
            prtangle= i;
            TString nid = Form("_%2.0d", i);
            for (int j=0; j<=100; j+=1) {
                TString calc_error_stat_sim_p_path = Form("/Users/ahmed/perforamnce/spr_data_sim/error_stat/statistics_error_theta_%i_seed_%d_3lsph_proton_sim_spr.root", i, j);
                cout<<"sim path for error calc p = " <<calc_error_stat_sim_p_path<<endl;
                string path_calc_error_stat_sim_p = (string)calc_error_stat_sim_p_path;
                cout<<"exists_test(path_calc_error_stat_sim_p)" <<exists_test(path_calc_error_stat_sim_p)<<endl;
                if (!exists_test(path_calc_error_stat_sim_p)) continue;
                ffile_error_stat_sim_p = TFile::Open(calc_error_stat_sim_p_path,"read");
                error_stat_p_cherenkov_sim_org=(TH1F*)ffile_error_stat_sim_p->Get("fHist");
                error_stat_p_cherenkov_sim_corrected=(TH1F*)ffile_error_stat_sim_p->Get("fHist_correction");
                error_stat_p_cherenkov_mc_same_path=(TH1F*)ffile_error_stat_sim_p->Get("fHist_same_path");
                
                error_stat_yield_wo_sim=(TH1F*)ffile_error_stat_sim_p->Get("fnHits");
                error_stat_yield_wt_sim=(TH1F*)ffile_error_stat_sim_p->Get("fnHits_p");
                error_stat_yield_wtc_sim=(TH1F*)ffile_error_stat_sim_p->Get("fnHits_p_good");
                Double_t mean_error_stat_yield_wo_sim = error_stat_yield_wo_sim->GetMean();
                Double_t mean_error_stat_yield_wt_sim = error_stat_yield_wt_sim->GetMean();
                Double_t mean_error_stat_yield_wtc_sim = error_stat_yield_wtc_sim->GetMean();
                histo_distribution_error_stat_p_yield_sim_wo[counter_error_stat]->Fill(mean_error_stat_yield_wo_sim);
                histo_distribution_error_stat_p_yield_sim_wt[counter_error_stat]->Fill(mean_error_stat_yield_wt_sim);
                histo_distribution_error_stat_p_yield_sim_wtc[counter_error_stat]->Fill(mean_error_stat_yield_wtc_sim);
                
                
                Double_t error_stat_cangle_MC_true =  error_stat_p_cherenkov_mc_same_path->GetXaxis()->GetBinCenter(error_stat_p_cherenkov_mc_same_path->GetMaximumBin());
                auto error_stat_resultPair= FitHisto_0( error_stat_p_cherenkov_sim_org , error_stat_cangle_MC_true, 0.06);
                Double_t error_stat_cangle_error_stat_p_cherenkov_sim_org = error_stat_resultPair.first;
                Double_t error_stat_spr_error_stat_p_cherenkov_sim_org = error_stat_resultPair.second;
                
                auto error_stat_resultPair_corrected= FitHisto_0( error_stat_p_cherenkov_sim_corrected , error_stat_cangle_MC_true, 0.06);
                Double_t error_stat_cangle_error_stat_p_cherenkov_sim_corrected = error_stat_resultPair_corrected.first;
                Double_t error_stat_spr_error_stat_p_cherenkov_sim_corrected = error_stat_resultPair_corrected.second;
                
                
                histo_distribution_error_stat_spr_p_cherenkov_sim_org[counter_error_stat]->Fill(error_stat_spr_error_stat_p_cherenkov_sim_org*1000.0);
                histo_distribution_error_stat_spr_p_cherenkov_sim_corrected[counter_error_stat]->Fill(error_stat_spr_error_stat_p_cherenkov_sim_corrected*1000.0);
                
                histo_distribution_error_stat_cangle_p_cherenkov_sim_org[counter_error_stat]->Fill(error_stat_cangle_error_stat_p_cherenkov_sim_org);
                histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected[counter_error_stat]->Fill(error_stat_cangle_error_stat_p_cherenkov_sim_corrected);
                
                delete error_stat_p_cherenkov_sim_org;
                delete error_stat_p_cherenkov_mc_same_path;
                delete error_stat_p_cherenkov_sim_corrected;
                
                delete error_stat_yield_wo_sim;
                delete error_stat_yield_wt_sim;
                delete error_stat_yield_wtc_sim;
                
                
                
                ///////////////////
                // Close files  ///
                ///////////////////
                delete ffile_error_stat_sim_p;
            }
            ++counter_error_stat;
        }
        //    for (Int_t e=0; e<=14; e++) {
        //        TString ko = Form("_%d", e);
        //        prt_canvasAdd("r_test_stat"+ko,800,400);
        //        histo_distribution_error_stat_spr_p_cherenkov_sim_org[e]->Draw();
        //    }
        
        ///////////////////////////
        // Error Polar angle    ///
        ///////////////////////////
        for (Int_t e=0; e<=14; e++) {
            histo_distribution_error_polar_spr_p_cherenkov_sim_org[e] = new TH1F(Form("histo_distribution_error_polar_spr_p_cherenkov_sim_org_%d",e),Form("histo_distribution_error_polar_spr_p_cherenkov_sim_org_%d; SPR [mrad];entries [#]",e), 100,6,14);
            histo_distribution_error_polar_spr_p_cherenkov_sim_corrected[e] = new TH1F(Form("histo_distribution_error_polar_spr_p_cherenkov_sim_corrected_%d",e),Form("histo_distribution_error_polar_spr_p_cherenkov_sim_corrected_%d; SPR [mrad];entries [#]",e), 100,6,14);
            // Cangle
            histo_distribution_error_polar_cangle_p_cherenkov_sim_org[e] = new TH1F(Form("histo_distribution_error_polar_cangle_p_cherenkov_sim_org_%d",e),Form("histo_distribution_error_polar_cangle_p_cherenkov_sim_org_%d; cangle [mrad];entries [#]",e), 100,0.8,0.85);
            histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected[e] = new TH1F(Form("histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected_%d",e),Form("histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected_%d; cangle [mrad];entries [#]",e),100,0.8,0.85);
            // photon yield
            histo_distribution_error_polar_p_yield_sim_wtc[e] = new TH1F(Form("histo_distribution_error_polar_p_yield_sim_wtc_%d",e),Form("histo_distribution_error_polar_p_yield_sim_wtc_%d; SPR [mrad];entries [#]",e), 500,0,100);
            histo_distribution_error_polar_p_yield_sim_wt[e] = new TH1F(Form("histo_distribution_error_polar_p_yield_sim_wt_%d",e),Form("histo_distribution_error_polar_p_yield_sim_wt_%d; SPR [mrad];entries [#]",e), 500,0,100);
            histo_distribution_error_polar_p_yield_sim_wo[e] = new TH1F(Form("histo_distribution_error_polar_p_yield_sim_wo_%d",e),Form("histo_distribution_error_polar_p_yield_sim_wo_%d; SPR [mrad];entries [#]",e), 500,0,100);
        }
        Int_t counter_error_polar=0;
        for (int i=20; i<=150; i+=10) {
            prtangle= i;
            TString nid = Form("_%2.0d", i);
            
            for (int j=-9; j<=9; j+=1) {
                //for (int j=-90; j<=90; j+=1) {
                Double_t jj = (Float_t) i+j/10.0 ;
                TString jj_string = Form("_theta_%.1f", jj);
                TString calc_error_polar_sim_p_path = "/Users/ahmed/perforamnce/spr_data_sim/error_polar/spr_polar_error"+jj_string+"_3lsph_proton_sim_spr.root";
                // for big ~ 180 file per point
                //Double_t kk = (Float_t) i+j/100.0 ;
                //TString kk_string = Form("_theta_%.2f", kk);
                //TString calc_error_polar_sim_p_path = "/Users/ahmed/perforamnce/spr_data_sim/error_polar_more/spr_polar_error"+kk_string+"_3lsph_proton_sim_spr.root";
                
                //cout<<"sim path for error calc p = " <<calc_error_polar_sim_p_path<<endl;
                string path_calc_error_polar_sim_p = (string)calc_error_polar_sim_p_path;
                cout<<"exists_test(path_calc_error_polar_sim_p)" <<exists_test(path_calc_error_polar_sim_p)<<endl;
                if (!exists_test(path_calc_error_polar_sim_p)) continue;
                ffile_error_polar_sim_p = TFile::Open(calc_error_polar_sim_p_path,"read");
                error_polar_p_cherenkov_sim_org=(TH1F*)ffile_error_polar_sim_p->Get("fHist");
                error_polar_p_cherenkov_sim_corrected=(TH1F*)ffile_error_polar_sim_p->Get("fHist_correction");
                error_polar_p_cherenkov_mc_same_path=(TH1F*)ffile_error_polar_sim_p->Get("fHist_same_path");
                
                error_polar_yield_wo_sim=(TH1F*)ffile_error_polar_sim_p->Get("fnHits");
                error_polar_yield_wt_sim=(TH1F*)ffile_error_polar_sim_p->Get("fnHits_p");
                error_polar_yield_wtc_sim=(TH1F*)ffile_error_polar_sim_p->Get("fnHits_p_good");
                Double_t mean_error_polar_yield_wo_sim = error_polar_yield_wo_sim->GetMean();
                Double_t mean_error_polar_yield_wt_sim = error_polar_yield_wt_sim->GetMean();
                Double_t mean_error_polar_yield_wtc_sim = error_polar_yield_wtc_sim->GetMean();
                histo_distribution_error_polar_p_yield_sim_wo[counter_error_polar]->Fill(mean_error_polar_yield_wo_sim);
                histo_distribution_error_polar_p_yield_sim_wt[counter_error_polar]->Fill(mean_error_polar_yield_wt_sim);
                histo_distribution_error_polar_p_yield_sim_wtc[counter_error_polar]->Fill(mean_error_polar_yield_wtc_sim);

                Double_t error_polar_cangle_MC_true =  error_polar_p_cherenkov_mc_same_path->GetXaxis()->GetBinCenter(error_polar_p_cherenkov_mc_same_path->GetMaximumBin());
                auto error_polar_resultPair= FitHisto_0( error_polar_p_cherenkov_sim_org , error_polar_cangle_MC_true, 0.06);
                Double_t error_polar_cangle_error_polar_p_cherenkov_sim_org = error_polar_resultPair.first;
                Double_t error_polar_spr_error_polar_p_cherenkov_sim_org = error_polar_resultPair.second;
                
                auto error_polar_resultPair_corrected= FitHisto_0( error_polar_p_cherenkov_sim_corrected , error_polar_cangle_MC_true, 0.06);
                Double_t error_polar_cangle_error_polar_p_cherenkov_sim_corrected = error_polar_resultPair_corrected.first;
                Double_t error_polar_spr_error_polar_p_cherenkov_sim_corrected = error_polar_resultPair_corrected.second;
                
                histo_distribution_error_polar_spr_p_cherenkov_sim_org[counter_error_polar]->Fill(error_polar_spr_error_polar_p_cherenkov_sim_org*1000.0);
                histo_distribution_error_polar_spr_p_cherenkov_sim_corrected[counter_error_polar]->Fill(error_polar_spr_error_polar_p_cherenkov_sim_corrected*1000.0);
                
                histo_distribution_error_polar_cangle_p_cherenkov_sim_org[counter_error_polar]->Fill(error_polar_cangle_error_polar_p_cherenkov_sim_org);
                histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected[counter_error_polar]->Fill(error_polar_cangle_error_polar_p_cherenkov_sim_corrected);
                
                delete error_polar_p_cherenkov_sim_org;
                delete error_polar_p_cherenkov_mc_same_path;
                delete error_polar_p_cherenkov_sim_corrected;
                
                delete error_polar_yield_wo_sim;
                delete error_polar_yield_wt_sim;
                delete error_polar_yield_wtc_sim;
                
                ///////////////////
                // Close files  ///
                ///////////////////
                delete ffile_error_polar_sim_p;
            }
            ++counter_error_polar;
        }
        
        
        //        for (Int_t e=1; e<=14; e++) {
        //            TString ko = Form("_%d", e);
        //            prt_canvasAdd("r_test_polar"+ko,800,400);
        //            histo_distribution_error_polar_spr_p_cherenkov_sim_org[e]->Draw();
        //        }
    }
    ///////////////////////////
    // Analysis             ///
    ///////////////////////////
    for (int i=20; i<=150; i+=10) {
        prtangle= i;
        TString nid = Form("_%2.0d", i);
        // proton
        TString cherenkov_data_path = Form("/Users/ahmed/perforamnce/spr_data_sim/spr_wtb_%d_sph_p_data_spr.root", i);
        TString cherenkov_sim_path = Form("/Users/ahmed/perforamnce/spr_data_sim/spr_wt_%d_sph_p_sim_spr.root", i);
        //TString cherenkov_sim_path = Form("/Users/ahmed/perforamnce/spr_data_sim/error_stat/statistics_error_theta_%i_seed_99_3lsph_proton_sim_spr.root", i);
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
        p_cherenkov_sim_org=(TH1F*)ffile_sim->Get("fHist");
        //p_cherenkov_sim_org->SetStats(0);
        p_cherenkov_data_org=(TH1F*)ffile_data->Get("fHist");
        p_cherenkov_data_org->SetStats(0);//kFALSE
        p_cherenkov_data_copy=(TH1F*)ffile_data->Get("fHist_copy");
        p_cherenkov_sim_copy=(TH1F*)ffile_sim->Get("fHist_copy");
        p_cherenkov_mc_same_path=(TH1F*)ffile_sim->Get("fHist_same_path");
        p_cherenkov_bg_sim=(TH1F*)ffile_sim->Get("fHist_bg");
        p_cherenkov_data_corrected=(TH1F*)ffile_data->Get("fHist_correction");
        p_cherenkov_sim_corrected=(TH1F*)ffile_sim->Get("fHist_correction");
        
        ////////////////////////////
        // photon yield histogram //
        ////////////////////////////
        p_yield_wo_sim=(TH1F*)ffile_sim->Get("fnHits");
        p_yield_wt_sim=(TH1F*)ffile_sim->Get("fnHits_p");
        p_yield_wtc_sim=(TH1F*)ffile_sim->Get("fnHits_p_good");
        p_yield_true_sim=(TH1F*)ffile_sim->Get("fnHits_true_sim");
        p_yield_wo_data=(TH1F*)ffile_data->Get("fnHits");
        p_yield_wt_data=(TH1F*)ffile_data->Get("fnHits_p");
        p_yield_wtc_data=(TH1F*)ffile_data->Get("fnHits_p_good");
        //dac_hits_data=(TH1F*)ffile_data->Get("nHist_dac");
        //dac_hits_sys_cus_data=(TH1F*)ffile_data->Get("nHist_dac_syscut_p");
        
        histo_style_photon_yield(p_yield_wo_sim, p_yield_wt_sim, p_yield_wtc_sim);
        histo_style_photon_yield(p_yield_wo_data, p_yield_wt_data, p_yield_wtc_data);
        
        //        hs10 = new THStack("hs10","Stacked 1D histograms");
        //        hs10->Add(p_yield_wo_sim);
        //        hs10->Add(p_yield_wt_sim);
        //        hs10->Add(p_yield_wtc_sim);
        //        prt_canvasAdd("r_test_asdsda"+nid,800,400);
        //        hs10->SetTitle(Form("Polar angle %3.1f", prtangle));
        //        hs10->Draw("nostack");
        //        hs10->GetXaxis()->SetTitle("#theta_{C} [rad]");
        //        hs10->GetYaxis()->SetTitle("entries [#]");
        //        prt_canvasGet("r_test_asdsda"+nid)->Update();
        
        
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
        hist_ambiguity_lut_sim =(TH1F*)ffile_sim->Get("hist_ambiguity");
        histo_photon_ambiguity_wo_sim=(TH1F*)ffile_sim->Get("histo_photon_ambiguity_wo");
        histo_photon_ambiguity_wt_sim=(TH1F*)ffile_sim->Get("histo_photon_ambiguity_wt");
        histo_photon_ambiguity_wtc_sim=(TH1F*)ffile_sim->Get("histo_photon_ambiguity_wtc");
        hist_ambiguity_lut_data =(TH1F*)ffile_data->Get("hist_ambiguity");
        histo_photon_ambiguity_wo_data=(TH1F*)ffile_data->Get("histo_photon_ambiguity_wo");
        histo_photon_ambiguity_wt_data=(TH1F*)ffile_data->Get("histo_photon_ambiguity_wt");
        histo_photon_ambiguity_wtc_data=(TH1F*)ffile_data->Get("histo_photon_ambiguity_wtc");
        histo_style_photon_yield(histo_photon_ambiguity_wo_sim, histo_photon_ambiguity_wt_sim, histo_photon_ambiguity_wtc_sim);
        histo_style_photon_yield(histo_photon_ambiguity_wo_data, histo_photon_ambiguity_wt_data, histo_photon_ambiguity_wtc_data);
        /////////////
        // TOF PID //
        /////////////
        tof_pid=(TH1F*)ffile_data->Get("hdelta_tof2tof1");
        
        //TSpectrum *fSpect= new TSpectrum(10);
        //Int_t nfound = fSpect->Search(p_cherenkov_mc_same_path,1,"",0.9); //0.6
        //Double_t cangle_MC_true = fSpect->GetPositionX()[0];
        Double_t cangle_MC_true =  p_cherenkov_mc_same_path->GetXaxis()->GetBinCenter(p_cherenkov_mc_same_path->GetMaximumBin());
        
        ///////////////////
        ///// part I //////
        ///////////////////
        if(bool_part1) {
            if (bool_part1_1) {
                //gStyle->SetOptFit(0);
                //gStyle->SetOptStat(0);
                DiffNorm(p_cherenkov_sim_org, p_cherenkov_data_org);
                histo_style_match(p_cherenkov_sim_org, p_cherenkov_data_org);
                TLegend *legend_ch_match= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                legend_ch_match->SetHeader("Cherenkov angle ambiguity distribution without #theta_{c} (proton)","C");
                prt_canvasAdd("r_ch_match"+nid,800,400);
                p_cherenkov_sim_org->SetTitle(Form("Polar angle %3.1f", prtangle));
                legend_ch_match->AddEntry(p_cherenkov_sim_org,"Sim","f");
                legend_ch_match->AddEntry(p_cherenkov_data_org,"Data ","f");
                p_cherenkov_sim_org->Draw("hist");
                p_cherenkov_data_org->Draw("samehist");
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
                prt_canvasAdd("r_corrected"+nid,800,400);
                p_cherenkov_data_corrected->SetTitle(Form("Polar angle %3.1f (proton data)", prtangle));
                histo_style_match(p_cherenkov_data_org,p_cherenkov_data_corrected);
                p_cherenkov_data_corrected->Draw();
                p_cherenkov_data_org->Draw("same");
                TLegend *legend_corrected= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                
                TF1 *fFit_p_cherenkov_data_corrected = new TF1("fFit_p_cherenkov_data_corrected","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
                fFit_p_cherenkov_data_corrected->SetLineColor(1);
                Double_t cangle_cor =  p_cherenkov_data_corrected->GetXaxis()->GetBinCenter(p_cherenkov_sim_org->GetMaximumBin());
                if(cangle_cor>0.85) cangle_cor=0.82;
                fFit_p_cherenkov_data_corrected->SetParameters(100,cangle_cor,0.010);
                fFit_p_cherenkov_data_corrected->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
                fFit_p_cherenkov_data_corrected->SetParLimits(0,0.1,1E6);
                fFit_p_cherenkov_data_corrected->SetParLimits(1,cangle_cor-0.04,cangle_cor+0.04);
                fFit_p_cherenkov_data_corrected->SetParLimits(2,0.005,0.014000);
                p_cherenkov_data_corrected->Fit("fFit_p_cherenkov_data_corrected","M","",cangle_cor-0.06,cangle_cor+0.06);
                
                legend_corrected->SetHeader(" Cherenkov angle ambiguity distribution (proton)","C");
                legend_corrected->AddEntry(p_cherenkov_data_org,"Data without #theta_{c} correction","f");
                legend_corrected->AddEntry(p_cherenkov_data_corrected,"Data after #theta_{c} correction","f");
                legend_corrected->Draw();
                prt_canvasGet("r_corrected"+nid)->Update();
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
                prt_canvasGet("r_corrected"+nid)->Update();
            }
            if(bool_part1_3) {
                DiffNorm(p_cherenkov_sim_corrected, p_cherenkov_data_corrected);
                prt_canvasAdd("r_corrected_match"+nid,800,400);
                p_cherenkov_data_corrected->SetTitle(Form("Polar angle %3.1f (proton data)", prtangle));
                histo_style_match(p_cherenkov_sim_corrected,p_cherenkov_data_corrected);
                p_cherenkov_data_corrected->Draw();
                p_cherenkov_sim_corrected->Draw("same");
                TLegend *legend_corrected= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                
                legend_corrected->SetHeader("Cherenkov angle ambiguity distribution after #theta_{c} correction (proton)","C");
                legend_corrected->AddEntry(p_cherenkov_data_corrected,"Data ","p");
                legend_corrected->AddEntry(p_cherenkov_sim_corrected,"Sim)","f");
                legend_corrected->Draw();
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
            /////////////////////////////
            // cherenkov normalization /
            /////////////////////////////
            TAxis *axis = p_cherenkov_sim_corrected->GetXaxis();
            double xmin = 0.6;
            double xmax = 0.74;
            if (i==90)xmin = 0.9;
            if (i==90)xmax = 1.0;
            int bmin = axis->FindBin(xmin);
            int bmax = axis->FindBin(xmax);
            double integral = p_cherenkov_sim_corrected->Integral(bmin,bmax);
            TAxis *axis_data = p_cherenkov_data_corrected->GetXaxis();
            double xmin_data = 0.6;
            double xmax_data = 0.74;
            if (i==90)xmin_data = 0.9;
            if (i==90)xmax_data = 1.0;
            int bmin_data = axis_data->FindBin(xmin_data);
            int bmax_data = axis_data->FindBin(xmax_data);
            double integral_data = p_cherenkov_data_corrected->Integral(bmin_data,bmax_data);
            Double_t norm= integral/integral_data ;
            //std::cout<<"############  norm= "<< norm <<std::endl;
            p_cherenkov_data_copy->Scale(norm);
            p_cherenkov_data_sub = (TH1F *)p_cherenkov_data_copy->Clone();
            p_cherenkov_data_sub->Add(p_cherenkov_bg_sim,-1);
            std::cout<<"############"<< " no problem 3 " <<std::endl;
            
            ////////////////////
            // fit Cherenkove //
            ///////////////////
            auto resultPair= FitHisto_m( p_cherenkov_data_corrected , cangle_MC_true, 0.06); // 0.06
            Double_t cangle_data_corrected = resultPair.first;
            Double_t spr_data_corrected = resultPair.second;
            
            resultPair= FitHisto_m( p_cherenkov_sim_corrected , cangle_MC_true, 0.06);
            Double_t cangle_sim_corrected = resultPair.first;
            Double_t spr_sim_corrected = resultPair.second;
            
            resultPair= FitHisto_m( p_cherenkov_mc_same_path , cangle_MC_true, 0.06);
            Double_t cangle_sim_true = resultPair.first;
            Double_t spr_sim_true = resultPair.second;
            
            resultPair= FitHisto_m( p_cherenkov_data_sub , cangle_MC_true, 0.06);
            Double_t cangle_data_sub = resultPair.first;
            Double_t spr_data_sub = resultPair.second;
            
            resultPair= FitHisto_m( p_cherenkov_data_copy , cangle_MC_true, 0.06);
            Double_t cangle_data_copy = resultPair.first;
            Double_t spr_data_copy = resultPair.second;
            
            resultPair= FitHisto_0( p_cherenkov_sim_org , cangle_MC_true, 0.06);
            Double_t cangle_sim_org = resultPair.first;
            Double_t spr_sim_org = resultPair.second;
            
            resultPair= FitHisto_0( p_cherenkov_data_org , cangle_MC_true, 0.06);
            Double_t cangle_data_org = resultPair.first;
            Double_t spr_data_org = resultPair.second;
            
            
            histo_style_cherenkov(p_cherenkov_data_copy, p_cherenkov_data_sub, p_cherenkov_bg_sim, p_cherenkov_mc_same_path, p_cherenkov_sim_corrected);
            
            //////////////////////////////
            // time diff normalization ///
            //////////////////////////////
            DiffNorm(p_diff_time_sim, p_diff_time_data);
            p_diff_time_data_sub = (TH1F *)p_diff_time_data->Clone();
            p_diff_time_data_sub->Add(p_diff_time_bg_sim,-1);
            //////////////////////////////
            // function histogram style///
            //////////////////////////////
            
            
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
            
            histo_style_time_diff(p_diff_time_sim, p_diff_time_data, p_diff_time_data_sub, p_diff_time_mctruth, p_diff_time_bg_sim);
            
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
            gROOT->SetBatch(0);
            //////////////////
            //  Fill graph //
            /////////////////
            
            x[counter]=i;
            if(bool_error){
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_spr_polar_error_org"+nid,800,400);
                    histo_distribution_error_polar_spr_p_cherenkov_sim_org[counter]->SetTitle(Form("SPR distibution as a result of polar angle variation between +/- 1˚ without correction @ angle %3.1f ", prtangle));
                    histo_distribution_error_polar_spr_p_cherenkov_sim_org[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_polar_spr_p_cherenkov_sim_org[counter]->GetEntries();
                    Double_t mean_error_polar_distribution = histo_distribution_error_polar_spr_p_cherenkov_sim_org[counter]->GetMean();
                    Double_t RMS_error_polar_distribution = histo_distribution_error_polar_spr_p_cherenkov_sim_org[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_polar_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_polar_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_spr_polar_error_org"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_spr_stat_error_org"+nid,800,400);
                    histo_distribution_error_stat_spr_p_cherenkov_sim_org[counter]->SetTitle(Form("SPR distibution as a result of using several statistical samples without correction @ angle %3.1f", prtangle));
                    histo_distribution_error_stat_spr_p_cherenkov_sim_org[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_stat_spr_p_cherenkov_sim_org[counter]->GetEntries();
                    Double_t mean_error_stat_distribution = histo_distribution_error_stat_spr_p_cherenkov_sim_org[counter]->GetMean();
                    Double_t RMS_error_stat_distribution = histo_distribution_error_stat_spr_p_cherenkov_sim_org[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_stat_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_stat_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_spr_stat_error_org"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_spr_fit_range_error_org"+nid,800,400);
                    histo_distribution_error_fit_range_spr_p_cherenkov_sim_org[counter]->SetTitle(Form("SPR distibution as a result of changing fitting range without correction @ angle %3.1f", prtangle));
                    histo_distribution_error_fit_range_spr_p_cherenkov_sim_org[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_fit_range_spr_p_cherenkov_sim_org[counter]->GetEntries();
                    Double_t mean_error_fit_range_distribution = histo_distribution_error_fit_range_spr_p_cherenkov_sim_org[counter]->GetMean();
                    Double_t RMS_error_fit_range_distribution = histo_distribution_error_fit_range_spr_p_cherenkov_sim_org[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_fit_range_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_fit_range_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_spr_fit_range_error_org"+nid)->Update();
                }
                ////////////////////////////
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_spr_polar_error_corrected"+nid,800,400);
                    histo_distribution_error_polar_spr_p_cherenkov_sim_corrected[counter]->SetTitle(Form("SPR distibution as a result of polar angle variation between +/- 1˚ with correction @ angle %3.1f", prtangle));
                    histo_distribution_error_polar_spr_p_cherenkov_sim_corrected[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_polar_spr_p_cherenkov_sim_corrected[counter]->GetEntries();
                    Double_t mean_error_polar_distribution = histo_distribution_error_polar_spr_p_cherenkov_sim_corrected[counter]->GetMean();
                    Double_t RMS_error_polar_distribution = histo_distribution_error_polar_spr_p_cherenkov_sim_corrected[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_polar_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_polar_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_spr_polar_error_corrected"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_spr_stat_error_corrected"+nid,800,400);
                    histo_distribution_error_stat_spr_p_cherenkov_sim_corrected[counter]->SetTitle(Form("SPR distibution as a result of using several statistical samples with correction @ angle %3.1f", prtangle));
                    histo_distribution_error_stat_spr_p_cherenkov_sim_corrected[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_stat_spr_p_cherenkov_sim_corrected[counter]->GetEntries();
                    Double_t mean_error_stat_distribution = histo_distribution_error_stat_spr_p_cherenkov_sim_corrected[counter]->GetMean();
                    Double_t RMS_error_stat_distribution = histo_distribution_error_stat_spr_p_cherenkov_sim_corrected[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_stat_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_stat_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_spr_stat_error_corrected"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_spr_fit_range_error_corrected"+nid,800,400);
                    histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected[counter]->SetTitle(Form("SPR distibution as a result of changing fitting range with correction @ angle %3.1f", prtangle));
                    histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected[counter]->GetEntries();
                    Double_t mean_error_fit_range_distribution = histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected[counter]->GetMean();
                    Double_t RMS_error_fit_range_distribution = histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_fit_range_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_fit_range_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_spr_fit_range_error_corrected"+nid)->Update();
                }
                // cangle
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_cangle_polar_error_org"+nid,800,400);
                    histo_distribution_error_polar_cangle_p_cherenkov_sim_org[counter]->SetTitle(Form("#theta_{c} distibution as a result of polar angle variation between +/- 1˚ without correction @ angle %3.1f", prtangle));
                    histo_distribution_error_polar_cangle_p_cherenkov_sim_org[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_polar_cangle_p_cherenkov_sim_org[counter]->GetEntries();
                    Double_t mean_error_polar_distribution = histo_distribution_error_polar_cangle_p_cherenkov_sim_org[counter]->GetMean();
                    Double_t RMS_error_polar_distribution = histo_distribution_error_polar_cangle_p_cherenkov_sim_org[counter]->GetRMS()*1000;
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_polar_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_polar_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_cangle_polar_error_org"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_cangle_stat_error_org"+nid,800,400);
                    histo_distribution_error_stat_cangle_p_cherenkov_sim_org[counter]->SetTitle(Form("#theta_{c} distibution as a result of using several statistical samples without correction @ angle %3.1f", prtangle));
                    histo_distribution_error_stat_cangle_p_cherenkov_sim_org[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_stat_cangle_p_cherenkov_sim_org[counter]->GetEntries();
                    Double_t mean_error_stat_distribution = histo_distribution_error_stat_cangle_p_cherenkov_sim_org[counter]->GetMean();
                    Double_t RMS_error_stat_distribution = histo_distribution_error_stat_cangle_p_cherenkov_sim_org[counter]->GetRMS()*1000;
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_stat_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_stat_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_cangle_stat_error_org"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_cangle_fit_range_error_org"+nid,800,400);
                    histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org[counter]->SetTitle(Form("#theta_{c} distibution as a result of changing fitting range without correction @ angle %3.1f", prtangle));
                    histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org[counter]->GetEntries();
                    Double_t mean_error_fit_range_distribution = histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org[counter]->GetMean();
                    Double_t RMS_error_fit_range_distribution = histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org[counter]->GetRMS()*1000;
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_fit_range_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_fit_range_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_cangle_fit_range_error_org"+nid)->Update();
                }
                ////////////////////////////
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_cangle_polar_error_corrected"+nid,800,400);
                    histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected[counter]->SetTitle(Form("#theta_{c} distibution as a result of polar angle variation between +/- 1˚ with correction @ angle %3.1f", prtangle));
                    histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected[counter]->GetEntries();
                    Double_t mean_error_polar_distribution = histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected[counter]->GetMean();
                    Double_t RMS_error_polar_distribution = histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected[counter]->GetRMS()*1000;
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_polar_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_polar_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_cangle_polar_error_corrected"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_cangle_stat_error_corrected"+nid,800,400);
                    histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected[counter]->SetTitle(Form("#theta_{c} distibution as a result of using several statistical samples with correction @ angle %3.1f", prtangle));
                    histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected[counter]->GetEntries();
                    Double_t mean_error_stat_distribution = histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected[counter]->GetMean();
                    Double_t RMS_error_stat_distribution = histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected[counter]->GetRMS()*1000;
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_stat_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_stat_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_cangle_stat_error_corrected"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_cangle_fit_range_error_corrected"+nid,800,400);
                    histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected[counter]->SetTitle(Form("#theta_{c} distibution as a result of changing fitting range with correction @ angle %3.1f", prtangle));
                    histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected[counter]->GetEntries();
                    Double_t mean_error_fit_range_distribution = histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected[counter]->GetMean();
                    Double_t RMS_error_fit_range_distribution = histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected[counter]->GetRMS()*1000;
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_fit_range_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_fit_range_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_cangle_fit_range_error_corrected"+nid)->Update();
                }
                // yield
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_yield_polar_wo"+nid,800,400);
                    histo_distribution_error_polar_p_yield_sim_wo[counter]->SetTitle(Form("photon yield distibution as a result of polar angle variation between +/- 1˚ without cuts @ angle %3.1f", prtangle));
                    histo_distribution_error_polar_p_yield_sim_wo[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_polar_p_yield_sim_wo[counter]->GetEntries();
                    Double_t mean_error_polar_distribution = histo_distribution_error_polar_p_yield_sim_wo[counter]->GetMean();
                    Double_t RMS_error_polar_distribution = histo_distribution_error_polar_p_yield_sim_wo[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_polar_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_polar_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_yield_polar_wo"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_yield_polar_wt"+nid,800,400);
                    histo_distribution_error_polar_p_yield_sim_wt[counter]->SetTitle(Form("photon yield distibution as a result of polar angle variation between +/- 1˚ after time cut @ angle %3.1f", prtangle));
                    histo_distribution_error_polar_p_yield_sim_wt[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_polar_p_yield_sim_wt[counter]->GetEntries();
                    Double_t mean_error_polar_distribution = histo_distribution_error_polar_p_yield_sim_wt[counter]->GetMean();
                    Double_t RMS_error_polar_distribution = histo_distribution_error_polar_p_yield_sim_wt[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_polar_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_polar_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_yield_polar_wt"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_yield_polar_wtc"+nid,800,400);
                    histo_distribution_error_polar_p_yield_sim_wtc[counter]->SetTitle(Form("photon yield distibution as a result of polar angle variation between +/- 1˚ after #theta_{c} & time cuts @ angle %3.1f", prtangle));
                    histo_distribution_error_polar_p_yield_sim_wtc[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_polar_p_yield_sim_wtc[counter]->GetEntries();
                    Double_t mean_error_polar_distribution = histo_distribution_error_polar_p_yield_sim_wtc[counter]->GetMean();
                    Double_t RMS_error_polar_distribution = histo_distribution_error_polar_p_yield_sim_wtc[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_polar_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_polar_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_yield_polar_wtc"+nid)->Update();
                }
                //// stat
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_yield_stat_wo"+nid,800,400);
                    histo_distribution_error_stat_p_yield_sim_wo[counter]->SetTitle(Form("photon yield distibution as a result of using several statistical samples without cuts @ angle %3.1f", prtangle));
                    histo_distribution_error_stat_p_yield_sim_wo[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_stat_p_yield_sim_wo[counter]->GetEntries();
                    Double_t mean_error_stat_distribution = histo_distribution_error_stat_p_yield_sim_wo[counter]->GetMean();
                    Double_t RMS_error_stat_distribution = histo_distribution_error_stat_p_yield_sim_wo[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_stat_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_stat_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_yield_stat_wo"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_yield_stat_wt"+nid,800,400);
                    histo_distribution_error_stat_p_yield_sim_wt[counter]->SetTitle(Form("photon yield distibution as a result of using several statistical samples after time cut @ angle %3.1f", prtangle));
                    histo_distribution_error_stat_p_yield_sim_wt[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_stat_p_yield_sim_wt[counter]->GetEntries();
                    Double_t mean_error_stat_distribution = histo_distribution_error_stat_p_yield_sim_wt[counter]->GetMean();
                    Double_t RMS_error_stat_distribution = histo_distribution_error_stat_p_yield_sim_wt[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_stat_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_stat_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_yield_stat_wt"+nid)->Update();
                }
                if (bool_error_histo){
                    gStyle->SetOptStat(0);
                    prt_canvasAdd("r_yield_stat_wtc"+nid,800,400);
                    histo_distribution_error_stat_p_yield_sim_wtc[counter]->SetTitle(Form("photon yield distibution as a result of using several statistical samples after #theta_{c} & time cuts @ angle %3.1f", prtangle));
                    histo_distribution_error_stat_p_yield_sim_wtc[counter]->Draw();
                    Int_t entries_distribution = histo_distribution_error_stat_p_yield_sim_wtc[counter]->GetEntries();
                    Double_t mean_error_stat_distribution = histo_distribution_error_stat_p_yield_sim_wtc[counter]->GetMean();
                    Double_t RMS_error_stat_distribution = histo_distribution_error_stat_p_yield_sim_wtc[counter]->GetRMS();
                    TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                    pt_data_corrected->AddText(Form("Entries =  %d [#]", entries_distribution));
                    pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_error_stat_distribution));
                    pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", RMS_error_stat_distribution));
                    pt_data_corrected->Draw();
                    prt_canvasGet("r_yield_stat_wtc"+nid)->Update();
                }
                ////////////////////////////////////////////
                //  extract sigma from histo distribution //
                ////////////////////////////////////////////
                // SPR
                Double_t error_polar_spr_sim_org= histo_distribution_error_polar_spr_p_cherenkov_sim_org[counter]->GetRMS();
                Double_t error_stat_spr_sim_org= histo_distribution_error_stat_spr_p_cherenkov_sim_org[counter]->GetRMS();
                Double_t error_fit_range_spr_sim_org= histo_distribution_error_fit_range_spr_p_cherenkov_sim_org[counter]->GetRMS();
                /////////
                Double_t error_polar_spr_sim_corrected= histo_distribution_error_polar_spr_p_cherenkov_sim_corrected[counter]->GetRMS();
                Double_t error_stat_spr_sim_corrected= histo_distribution_error_stat_spr_p_cherenkov_sim_corrected[counter]->GetRMS();
                Double_t error_fit_range_spr_sim_corrected= histo_distribution_error_fit_range_spr_p_cherenkov_sim_corrected[counter]->GetRMS();
                // cangle
                Double_t error_polar_cangle_sim_org= histo_distribution_error_polar_cangle_p_cherenkov_sim_org[counter]->GetRMS();
                Double_t error_stat_cangle_sim_org= histo_distribution_error_stat_cangle_p_cherenkov_sim_org[counter]->GetRMS();
                Double_t error_fit_range_cangle_sim_org= histo_distribution_error_fit_range_cangle_p_cherenkov_sim_org[counter]->GetRMS();
                /////////
                Double_t error_polar_cangle_sim_corrected= histo_distribution_error_polar_cangle_p_cherenkov_sim_corrected[counter]->GetRMS();
                Double_t error_stat_cangle_sim_corrected= histo_distribution_error_stat_cangle_p_cherenkov_sim_corrected[counter]->GetRMS();
                Double_t error_fit_range_cangle_sim_corrected= histo_distribution_error_fit_range_cangle_p_cherenkov_sim_corrected[counter]->GetRMS();
                // Yield
                Double_t error_polar_p_yield_sim_wo= histo_distribution_error_polar_p_yield_sim_wo[counter]->GetRMS();
                Double_t error_stat_p_yield_sim_wo= histo_distribution_error_stat_p_yield_sim_wo[counter]->GetRMS();
                /////////
                Double_t error_polar_p_yield_sim_wt= histo_distribution_error_polar_p_yield_sim_wt[counter]->GetRMS();
                Double_t error_stat_p_yield_sim_wt= histo_distribution_error_stat_p_yield_sim_wt[counter]->GetRMS();
                /////////
                Double_t error_polar_p_yield_sim_wtc= histo_distribution_error_polar_p_yield_sim_wtc[counter]->GetRMS();
                Double_t error_stat_p_yield_sim_wtc= histo_distribution_error_stat_p_yield_sim_wtc[counter]->GetRMS();
                //////////////////////
                //  Error equations //
                //////////////////////
                // spr
                error_all_y_spr_sim_org[counter]= sqrt (error_polar_spr_sim_org*error_polar_spr_sim_org + error_stat_spr_sim_org * error_stat_spr_sim_org + error_fit_range_spr_sim_org * error_fit_range_spr_sim_org ) ;
                error_all_y_spr_sim_corrected[counter]= sqrt (error_polar_spr_sim_corrected*error_polar_spr_sim_corrected + error_stat_spr_sim_corrected * error_stat_spr_sim_corrected + error_fit_range_spr_sim_corrected * error_fit_range_spr_sim_corrected ) ;
                // cangle
                error_all_y_cangle_sim_org[counter]= sqrt (error_polar_cangle_sim_org*error_polar_cangle_sim_org + error_stat_cangle_sim_org * error_stat_cangle_sim_org + error_fit_range_cangle_sim_org * error_fit_range_cangle_sim_org ) ;
                error_all_y_cangle_sim_corrected[counter]= sqrt (error_polar_cangle_sim_corrected*error_polar_cangle_sim_corrected + error_stat_cangle_sim_corrected * error_stat_cangle_sim_corrected + error_fit_range_cangle_sim_corrected * error_fit_range_cangle_sim_corrected ) ;
                // yield
                error_all_y_yield_wo_sim_corrected[counter]=sqrt(error_polar_p_yield_sim_wo*error_polar_p_yield_sim_wo + error_stat_p_yield_sim_wo *error_stat_p_yield_sim_wo);
                error_all_y_yield_wt_sim_corrected[counter]=sqrt(error_polar_p_yield_sim_wt*error_polar_p_yield_sim_wt + error_stat_p_yield_sim_wo *error_stat_p_yield_sim_wt);
                error_all_y_yield_wtc_sim_corrected[counter]=sqrt(error_polar_p_yield_sim_wtc*error_polar_p_yield_sim_wtc + error_stat_p_yield_sim_wo *error_stat_p_yield_sim_wtc);
            }
            /////////////////////////
            //  Fill graphs arrays //
            /////////////////////////
            //SPR
            // error calculated
            y_spr_sim_org[counter]=spr_sim_org*1000.0;
            y_spr_sim_corrected[counter]=spr_sim_corrected*1000.0;
            // error based on sim
            y_spr_data_org[counter]=spr_data_org*1000.0;
            y_spr_data_corrected[counter]=spr_data_corrected*1000.0;
            // no error calculations yiet
            y_spr_sim_true[counter]=spr_sim_true*1000.0;
            y_spr_data_sub[counter]=spr_data_sub*1000.0;
            // yield
            //One method: return mean of the histogram
            // error calculated
            y_yield_wo_sim[counter]=p_yield_wo_sim->GetMean();
            y_yield_wt_sim[counter]=p_yield_wt_sim->GetMean();
            y_yield_wtc_sim[counter]=p_yield_wtc_sim->GetMean();
            y_yield_true_sim[counter]=p_yield_true_sim->GetMean();
            // error based on sim
            y_yield_wo_data[counter]=p_yield_wo_data->GetMean();
            y_yield_wt_data[counter]=p_yield_wt_data->GetMean();
            y_yield_wtc_data[counter]=p_yield_wtc_data->GetMean();
            // cangle
            // error calculated
            y_cangle_sim_org[counter]=cangle_sim_org;
            y_cangle_sim_corrected[counter]=cangle_sim_corrected;
            // error based on sim
            y_cangle_data_org[counter]=cangle_data_org;
            y_cangle_data_corrected[counter]=cangle_data_corrected;
            // no error calculations yiet
            y_cangle_sim_true[counter]=cangle_sim_true;
            y_cangle_data_sub[counter]=cangle_data_sub;

            y_diff_true[counter]=sigma_p_diff_time_mctruth*5;
            y_diff_data[counter]=sigma_p_diff_time_data*5;
            y_diff_sim[counter]=sigma_p_diff_time_sim*5;
            y_mean_diff_data[counter]=mean_p_diff_time_data;
            y_mean_diff_sim[counter]=mean_p_diff_time_sim;
            y_mean_diff_true[counter]=mean_p_diff_time_mctruth;
            
            calc_mom->SetPoint(counter,i,calc_p_tof2tof1_plot1_co);
            calc_tof1tof2_distance->SetPoint(counter,i,distance_tof2tof1_plot1_co);

            std::cout<<"############  counter = "<< counter <<std::endl;
            std::cout<<"############  x[counter] = "<< x[counter] <<std::endl;
            std::cout<<"############  y_spr_data_org[counter] = "<< y_spr_data_org[counter] <<std::endl;
            
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
            //hs10 = new THStack("hs10","Stacked 1D histograms");
            
            hs->Add(p_cherenkov_bg_sim);
            hs->Add(p_cherenkov_data_sub);
            hs->Add(p_cherenkov_data_copy);
            
            hs2->Add(p_yield_wo_sim);//sim
            hs2->Add(p_yield_wt_sim);
            hs2->Add(p_yield_wtc_sim);
            
            
            hs6->Add(p_yield_wo_data);
            hs6->Add(p_yield_wt_data);
            hs6->Add(p_yield_wtc_data);
            
            //hs2->Add(p_yield_true_sim);
            hs3->Add(p_diff_time_sim);
            hs3->Add(p_diff_time_bg_sim);
            hs3->Add(p_diff_time_data);
            hs3->Add(p_diff_time_mctruth);
            //hs3->Add(p_diff_time_data_sub);
            hs4->Add(p_photon_time_sim_calc);
            hs4->Add(p_photon_time_data_calc);
            hs5->Add(p_photon_time_sim);
            hs5->Add(p_photon_time_data);
            histo_style_photon_time(p_photon_time_sim,p_photon_time_data, p_photon_time_sim_calc, p_photon_time_data_calc);
            
            //hs6->Add(dac_hits_data);
            //hs6->Add(dac_hits_sys_cus_data);
            hs7->Add(p_cherenkov_bg_sim);
            hs7->Add(p_cherenkov_sim_corrected);
            hs7->Add(p_cherenkov_mc_same_path);
            hs8->Add(histo_photon_ambiguity_wo_data);
            hs8->Add(histo_photon_ambiguity_wt_data);
            hs8->Add(histo_photon_ambiguity_wtc_data);
            hs9->Add(histo_photon_ambiguity_wo_sim);
            hs9->Add(histo_photon_ambiguity_wt_sim);
            hs9->Add(histo_photon_ambiguity_wtc_sim);
            //hs10->Add(hist_ambiguity_lut_sim);
            //hs10->Add(hist_ambiguity_lut_data);
            
            //////////////////////
            // Fit Cherenkove  ///
            //////////////////////
            if(bool_partII){
                gStyle->SetOptFit(0);
                gStyle->SetOptStat(0);
                
                prt_canvasAdd("r_ch_fit_data"+nid,800,400);
                p_cherenkov_data_corrected->SetTitle(Form("Polar angle %3.1f", prtangle));
                p_cherenkov_data_corrected->Draw();
                
                Double_t mean_p_cherenkov_data_corrected = cangle_data_corrected;
                Double_t sigma_p_cherenkov_data_corrected = spr_data_corrected* 1000;
                TPaveText *pt_data_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                pt_data_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_p_cherenkov_data_corrected));
                pt_data_corrected->AddText(Form("#sigma =  %1.3f [mrad]", sigma_p_cherenkov_data_corrected));
                pt_data_corrected->Draw();
                prt_canvasGet("r_ch_fit_data"+nid)->Update();
                
                prt_canvasAdd("r_ch_fit_sim"+nid,800,400);
                p_cherenkov_sim_corrected->SetTitle(Form("Polar angle %3.1f", prtangle));
                p_cherenkov_sim_corrected->Draw();
                Double_t mean_p_cherenkov_sim_corrected = cangle_sim_corrected;
                Double_t sigma_p_cherenkov_sim_corrected = spr_sim_corrected* 1000;
                TPaveText *pt_sim_corrected = new TPaveText(0.703008,0.685333,0.958647,0.946667, "NDC");
                pt_sim_corrected->AddText(Form("#mu =  %1.3f [rad]", mean_p_cherenkov_sim_corrected));
                pt_sim_corrected->AddText(Form("#sigma =  %1.3f [mrad]", sigma_p_cherenkov_sim_corrected));
                pt_sim_corrected->Draw();
                prt_canvasGet("r_ch_fit_sim"+nid)->Update();
            }
            
            ///////////////
            // TOF PID  ///
            ///////////////
            if (bool_partII_histo) {
                gStyle->SetOptFit(0);
                gStyle->SetOptStat(0);
                prt_canvasAdd("r_tof_pid"+nid,800,400);
                tof_pid->SetTitle(Form("Polar angle %3.1f (data)", prtangle));
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
            if(bool_partII_histo) {
                gStyle->SetOptFit(0);
                gStyle->SetOptStat(0);
                gPad->UseCurrentStyle();
                TLegend *legend_signal= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                legend_signal->SetHeader("Cherenkov angle ambiguity distribution (proton)","C");
                prt_canvasAdd("r_mctruth"+nid,800,400);
                p_cherenkov_sim_org->SetTitle(Form("Polar angle %3.1f", prtangle));
                legend_signal->AddEntry(p_cherenkov_sim_corrected,"Sim","f");
                legend_signal->AddEntry(p_cherenkov_mc_same_path,"Sim true path inside the prism ","f");
                legend_signal->AddEntry(p_cherenkov_bg_sim,"Sim combinatorial background ","f");
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
            if(bool_partII_histo) {
                gStyle->SetOptFit(0);
                gStyle->SetOptStat(0);
                gPad->UseCurrentStyle();
                TLegend * legend_bg_sub= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                legend_bg_sub->SetHeader("Cherenkov angle ambiguity distribution (proton)","C");
                legend_bg_sub->AddEntry(p_cherenkov_data_copy,"Data","p");
                legend_bg_sub->AddEntry(p_cherenkov_data_sub,"Data - BG","p");
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
            if (bool_partII_histo) { // done
                gStyle->SetOptFit(0);
                gStyle->SetOptStat(0);
                TLegend *legend_diff_bg_sub= new TLegend( 0.121554, 0.716578, 0.457393, 0.879679);
                legend_diff_bg_sub->SetHeader("Time difference distribution (proton)","C");
                legend_diff_bg_sub->AddEntry(p_diff_time_data,"Data","p");
                //legend_diff_bg_sub->AddEntry(p_diff_time_data_sub,"data - BG","p");
                legend_diff_bg_sub->AddEntry(p_diff_time_bg_sim," Sim BG ","f");
                legend_diff_bg_sub->AddEntry(p_diff_time_sim," Sim ","f");
                legend_diff_bg_sub->AddEntry(p_diff_time_mctruth," Sim true path inside prism ","f");
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
            if (bool_partII_histo) {
                TLegend *legend_photon_time_calc= new TLegend(0.556391, 0.712, 0.890977, 0.874667);
                legend_photon_time_calc->SetHeader("Calculated time match (proton)", "C");
                legend_photon_time_calc->AddEntry(p_photon_time_sim_calc,"Sim","f");
                legend_photon_time_calc->AddEntry(p_photon_time_data_calc,"Data","p");
                prt_canvasAdd("r_photon_calc_time"+nid,800,400);
                hs4->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs4->Draw("nostack");
                hs4->GetXaxis()->SetTitle("t_{calc} [ns]");
                hs4->GetYaxis()->SetTitle("entries [#]");
                legend_photon_time_calc->Draw();
            }
            if (bool_partII_histo) {
                TLegend *legend_photon_time=  new TLegend(0.556391, 0.712, 0.890977, 0.874667);
                legend_photon_time->SetHeader("Measured time match (proton)","C");
                legend_photon_time->AddEntry(p_photon_time_sim,"Sim","f");
                legend_photon_time->AddEntry(p_photon_time_data,"Data","p");
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
            if (bool_partII_histo) {
                TLegend * legend_yield_sim= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
                legend_yield_sim->SetHeader("Sim | DIRC photon solutions (proton)","C");
                prt_canvasAdd("r_fnHits_sim"+nid,800,400);
                p_yield_wo_sim->SetTitle(Form("Polar angle %3.1f (sim)", prtangle));
                legend_yield_sim->AddEntry(p_yield_wo_sim,"Without cuts  ","l");
                legend_yield_sim->AddEntry(p_yield_wt_sim,"With time cut ","l");
                legend_yield_sim->AddEntry(p_yield_wtc_sim,"With #theta_{C} and time cuts","l");
                //legend_yield_sim->AddEntry(p_yield_true_sim," DIRC hits true path inside prism","l");
                hs2->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs2->Draw("nostack");
                hs2->GetYaxis()->SetTitle("entries [#]");
                hs2->GetXaxis()->SetTitle("number of photon solutions per track [#]");
                legend_yield_sim->Draw();
            }
            if (bool_partII_histo) {
                TLegend * legend_yield_data= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
                legend_yield_data->SetHeader("Data | DIRC photon solutions (proton)","C");
                prt_canvasAdd("r_fnHits_data"+nid,800,400);
                hs6->SetTitle(Form("Polar angle %3.1f (data)", prtangle));
                legend_yield_data->AddEntry(p_yield_wo_data,"Without cuts ","l");
                legend_yield_data->AddEntry(p_yield_wt_data,"With time cut ","l");
                legend_yield_data->AddEntry(p_yield_wtc_data,"With #theta_{C} and time cuts ","l");
                //legend_yield_data->AddEntry(dac_hits_data,"DAC hits without cuts","f");
                //legend_yield_data->AddEntry(dac_hits_sys_cus_data,"DAC hits with auxiliary detectors cuts and PID cut(proton data)","f");
                hs6->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs6->Draw("nostack");
                hs6->GetYaxis()->SetTitle("entries [#]");
                hs6->GetXaxis()->SetTitle("number of photon solutions per track [#]");
                legend_yield_data->Draw();
            }
            if (bool_partII_histo) {
                TLegend * legend_ambiguity_data= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
                legend_ambiguity_data->SetHeader("Data | Number of photon solutions (proton)","C");
                prt_canvasAdd("r_ambiguity_data"+nid,800,400);
                hs8->SetTitle(Form("Polar angle %3.1f (data)", prtangle));
                legend_ambiguity_data->AddEntry(histo_photon_ambiguity_wo_data,"Without cuts ","l");
                legend_ambiguity_data->AddEntry(histo_photon_ambiguity_wt_data,"With time cut ","l");
                legend_ambiguity_data->AddEntry(histo_photon_ambiguity_wtc_data,"With #theta_{C} and time cuts ","l");
                hs8->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs8->Draw("nostack");
                hs8->GetYaxis()->SetTitle("entries [#]");
                hs8->GetXaxis()->SetTitle("number of photon solutions  [#]");
                legend_ambiguity_data->Draw();
            }
            if (bool_partII_histo) {
                TLegend * legend_ambiguity_sim= new TLegend(0.552632, 0.606952,  0.992481,   0.903743  );
                legend_ambiguity_sim->SetHeader("Sim | Number of photon solutions (proton)","C");
                prt_canvasAdd("r_ambiguity_sim"+nid,800,400);
                hs9->SetTitle(Form("Polar angle %3.1f (data)", prtangle));
                legend_ambiguity_sim->AddEntry(histo_photon_ambiguity_wo_sim,"Without cuts ","l");
                legend_ambiguity_sim->AddEntry(histo_photon_ambiguity_wt_sim,"With time cut ","l");
                legend_ambiguity_sim->AddEntry(histo_photon_ambiguity_wtc_sim,"With #theta_{C} and time cuts ","l");
                hs9->SetTitle(Form("Polar angle %3.1f", prtangle));
                hs9->Draw("nostack");
                hs9->GetYaxis()->SetTitle("entries [#]");
                hs9->GetXaxis()->SetTitle("number of photon solutions  [#]");
                legend_ambiguity_sim->Draw();
            }
            
            
            if (bool_partII_histo) {
                TLegend * legend_lut_ambiguity= new TLegend(0.135338, 0.653333,  0.454887,   0.872  );
                //legend_lut_ambiguity->SetHeader("Number of LUT photon solutions for each pixel","C");
                prt_canvasAdd("r_lut_ambiguity"+nid,800,400);
                hist_ambiguity_lut_sim->SetTitle("Number of LUT's photon solutions for each pixel");
                hist_ambiguity_lut_sim->Draw();
                //legend_lut_ambiguity->Draw();
            }
            c1->cd();
        }
    }
    
    ///////////////////
    ///// part III ////
    ///////////////////
    if(bool_partIII) {
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
        
        
        graph_spr_data_sub = new TGraph(n,x,y_spr_data_sub);
        graph_spr_sim_true = new TGraph(n,x,y_spr_sim_true);
        //graph_spr_data_corrected = new TGraph(n,x,y_spr_data_corrected);
        //graph_spr_data_org = new TGraph(n,x,y_spr_data_org);
        
        graph_diff_true = new TGraph(n,x,y_diff_true);
        graph_diff_data = new TGraph(n,x,y_diff_data);
        graph_diff_sim = new TGraph(n,x,y_diff_sim);
        
        graph_diff_data_mean = new TGraph(n,x,y_mean_diff_data);
        graph_diff_sim_mean = new TGraph(n,x,y_mean_diff_sim);
        
        graph_diff_true_mean = new TGraph(n,x,y_mean_diff_true);
        graph_diff_true_sim = new TGraph(n,x,y_mean_diff_sim);
        
        graph_cangle_data_sub = new TGraph(n,x,y_cangle_data_sub);
        graph_cangle_sim_org = new TGraph(n,x,y_cangle_sim_org);
        graph_cangle_sim_true = new TGraph(n,x,y_cangle_sim_true);
        graph_cangle_sim_corrected = new TGraph(n,x,y_cangle_sim_corrected);
        
        
        const Int_t dd = 10;
        Double_t xx[dd]  = {-0.22, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
        Double_t yy[dd]  = {1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
        
        Double_t ex[dd] = {1,1,1,1,1,1,1,1,1,1};
        
        Double_t ey[dd] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
        
        TGraphErrors *graph_spr_sim_org = new TGraphErrors(n,x,y_spr_sim_org,ex,error_all_y_spr_sim_org);
        TGraphErrors *graph_spr_data_org = new TGraphErrors(n,x,y_spr_data_org,ex,error_all_y_spr_sim_org);
        TGraphErrors *graph_spr_sim_corrected = new TGraphErrors(n,x,y_spr_sim_corrected,ex,error_all_y_spr_sim_corrected);
        TGraphErrors *graph_spr_data_corrected = new TGraphErrors(n,x,y_spr_data_corrected,ex,error_all_y_spr_sim_corrected);
        
        TGraphErrors *graph_yield_DIRC_wo_sim = new TGraphErrors(n,x,y_yield_wo_sim,ex,error_all_y_yield_wo_sim_corrected);
        TGraphErrors *graph_yield_DIRC_wt_sim = new TGraphErrors(n,x,y_yield_wt_sim,ex,error_all_y_yield_wt_sim_corrected);
        TGraphErrors *graph_yield_DIRC_wtc_sim = new TGraphErrors(n,x,y_yield_wtc_sim,ex,error_all_y_yield_wtc_sim_corrected);
        
        TGraphErrors *graph_yield_DIRC_true_sim     = new TGraphErrors(n,x,y_yield_true_sim,ex,error_all_y_yield_wtc_sim_corrected);
        TGraphErrors *graph_yield_DIRC_wo_data = new TGraphErrors(n,x,y_yield_wo_data,ex,error_all_y_yield_wo_sim_corrected);
        TGraphErrors *graph_yield_DIRC_wt_data = new TGraphErrors(n,x,y_yield_wt_data,ex,error_all_y_yield_wt_sim_corrected);
        TGraphErrors *graph_yield_DIRC_wtc_data = new TGraphErrors(n,x,y_yield_wtc_data,ex,error_all_y_yield_wtc_sim_corrected);
        
        TGraphErrors *graph_cangle_data_org = new TGraphErrors(n,x,y_cangle_data_org, ex, error_all_y_cangle_sim_org);
        TGraphErrors *graph_cangle_data_corrected = new TGraphErrors(n,x,y_cangle_data_corrected, ex, error_all_y_cangle_sim_corrected);
        
        
        
        graph_spr_sim_org->SetLineColor(1);
        graph_spr_sim_org->SetLineWidth(1);
        graph_spr_sim_org->SetMarkerColor(1);
        graph_spr_sim_org->SetMarkerStyle(4);
        graph_spr_sim_org->SetLineStyle(1);
        
        graph_spr_data_sub->SetLineColor(2);
        graph_spr_data_sub->SetLineWidth(1);
        graph_spr_data_sub->SetMarkerColor(2);
        graph_spr_data_sub->SetMarkerStyle(4);
        
        graph_spr_sim_true->SetLineColor(2);
        graph_spr_sim_true->SetLineWidth(1);
        graph_spr_sim_true->SetMarkerColor(2);
        graph_spr_sim_true->SetMarkerStyle(4);
        graph_spr_sim_true->SetLineStyle(1);
        
        //////////////////////////////////
        
        graph_spr_data_corrected->SetLineColor(kBlack);
        graph_spr_data_corrected->SetLineWidth(1);
        graph_spr_data_corrected->SetMarkerColor(kBlack);
        graph_spr_data_corrected->SetMarkerStyle(4);
        graph_spr_data_corrected->SetLineStyle(1);
        
        graph_spr_sim_corrected->SetLineColor(kMagenta);
        graph_spr_sim_corrected->SetLineWidth(1);
        graph_spr_sim_corrected->SetMarkerColor(kMagenta);
        graph_spr_sim_corrected->SetMarkerStyle(4);
        graph_spr_sim_corrected->SetLineStyle(2);
        
        graph_spr_data_org->SetLineColor(kRed);
        graph_spr_data_org->SetLineWidth(1);
        graph_spr_data_org->SetMarkerColor(kRed);
        graph_spr_data_org->SetMarkerStyle(4);
        graph_spr_data_org->SetLineStyle(1);
        
        ///////////////////////////////
        graph_diff_true->SetLineColor(2);
        graph_diff_true->SetLineWidth(1);
        graph_diff_true->SetMarkerColor(2);
        graph_diff_true->SetMarkerStyle(4);
        
        graph_diff_data->SetLineColor(1);
        graph_diff_data->SetLineWidth(1);
        graph_diff_data->SetMarkerColor(1);
        graph_diff_data->SetMarkerStyle(4);
        
        graph_diff_sim->SetLineColor(kMagenta);
        graph_diff_sim->SetLineWidth(1);
        graph_diff_sim->SetMarkerColor(kMagenta);
        graph_diff_sim->SetMarkerStyle(4);
        ////////////////////////////////////
        graph_diff_data_mean->SetLineColor(1);
        graph_diff_data_mean->SetLineWidth(1);
        graph_diff_data_mean->SetMarkerColor(1);
        graph_diff_data_mean->SetMarkerStyle(4);
        
        graph_diff_sim_mean->SetLineColor(6);
        graph_diff_sim_mean->SetLineWidth(1);
        graph_diff_sim_mean->SetMarkerColor(6);
        graph_diff_sim_mean->SetMarkerStyle(4);
        
        graph_diff_true_mean->SetLineColor(2);
        graph_diff_true_mean->SetLineWidth(1);
        graph_diff_true_mean->SetMarkerColor(2);
        graph_diff_true_mean->SetMarkerStyle(4);
        
        graph_cangle_data_sub->SetLineColor(1);
        graph_cangle_data_sub->SetLineWidth(1);
        graph_cangle_data_sub->SetMarkerColor(1);
        graph_cangle_data_sub->SetMarkerStyle(4);
        
        graph_cangle_data_corrected->SetLineColor(kBlack);
        graph_cangle_data_corrected->SetLineWidth(1);
        graph_cangle_data_corrected->SetMarkerColor(kBlack);
        graph_cangle_data_corrected->SetMarkerStyle(4);
        graph_cangle_data_corrected->SetLineStyle(1);
        
        graph_cangle_sim_corrected->SetLineColor(kBlue);
        graph_cangle_sim_corrected->SetLineWidth(1);
        graph_cangle_sim_corrected->SetMarkerColor(kBlue);
        graph_cangle_sim_corrected->SetMarkerStyle(4);
        graph_cangle_sim_corrected->SetLineStyle(2);
        
        graph_cangle_data_org->SetLineColor(kRed);
        graph_cangle_data_org->SetLineWidth(1);
        graph_cangle_data_org->SetMarkerColor(kRed);
        graph_cangle_data_org->SetMarkerStyle(4);
        graph_cangle_data_org->SetLineStyle(1);
        
        graph_cangle_sim_true->SetLineColor(1);
        graph_cangle_sim_true->SetLineWidth(1);
        graph_cangle_sim_true->SetMarkerColor(1);
        graph_cangle_sim_true->SetMarkerStyle(4);
        graph_cangle_sim_true->SetLineStyle(1);
        
        graph_cangle_sim_org->SetLineColor(2);
        graph_cangle_sim_org->SetLineWidth(1);
        graph_cangle_sim_org->SetMarkerColor(2);
        graph_cangle_sim_org->SetMarkerStyle(4);
        graph_cangle_sim_org->SetLineStyle(1);
        
        graph_yield_DIRC_wo_sim->SetLineColor(kBlack);
        graph_yield_DIRC_wt_sim->SetLineColor(kBlue);
        graph_yield_DIRC_wtc_sim->SetLineColor(kMagenta);
        graph_yield_DIRC_true_sim->SetLineColor(kRed);
        graph_yield_DIRC_wo_data->SetLineColor(kRed);
        graph_yield_DIRC_wt_data->SetLineColor(kBlue);
        graph_yield_DIRC_wtc_data->SetLineColor(kBlack);
        graph_yield_DIRC_wtc_data->SetLineStyle(2);
        
        
        graph_yield_DIRC_wo_sim->SetMarkerColor(kBlack);
        graph_yield_DIRC_wt_sim->SetMarkerColor(kBlue);
        graph_yield_DIRC_wtc_sim->SetMarkerColor(kMagenta);
        graph_yield_DIRC_true_sim->SetMarkerColor(kRed);
        graph_yield_DIRC_wo_data->SetMarkerColor(kRed);
        graph_yield_DIRC_wt_data->SetMarkerColor(kGreen);
        graph_yield_DIRC_wtc_data->SetMarkerColor(kBlack);
        
        
        graph_yield_DIRC_wo_sim->SetMarkerStyle(4);
        graph_yield_DIRC_wt_sim->SetMarkerStyle(4);
        graph_yield_DIRC_wtc_sim->SetMarkerStyle(4);
        graph_yield_DIRC_wo_data->SetMarkerStyle(4);
        graph_yield_DIRC_wt_data->SetMarkerStyle(4);
        graph_yield_DIRC_wtc_data->SetMarkerStyle(4);
        graph_yield_DIRC_true_sim->SetMarkerStyle(kRed);
        
        
        graph_yield_DIRC_wo_sim->SetLineWidth(1);
        graph_yield_DIRC_wt_sim->SetLineWidth(1);
        graph_yield_DIRC_wtc_sim->SetLineWidth(1);
        graph_yield_DIRC_wo_data->SetLineWidth(1);
        graph_yield_DIRC_wt_data->SetLineWidth(1);
        graph_yield_DIRC_wtc_data->SetLineWidth(1);
        graph_yield_DIRC_true_sim->SetLineWidth(1);
        
        graph_yield_DIRC_wt_data->SetLineStyle(2);
        
        
        graph_yield_DIRC_true_sim->SetLineStyle(1);
        
        
        //////////////
        //// cangle //
        //////////////
        
        if(bool_partIII) {
            
            graph_cangle_data_org->SetLineColor(kMagenta);
            graph_cangle_data_corrected->SetLineColor(kBlack);
            graph_cangle_sim_org->SetLineColor(kMagenta);
            graph_cangle_sim_corrected->SetLineColor(kBlack);
            graph_cangle_sim_org->SetLineStyle(2);
            graph_cangle_sim_corrected->SetLineStyle(2);
            
            
            prt_canvasAdd("r_graph_cangle_data",800,400);
            TLegend *leg_cangle_data_sub = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_cangle_data_sub->SetHeader("Data | Reconstructed #theta_{c} angle (proton)","C");
            leg_cangle_data_sub->AddEntry(graph_cangle_data_org, "Before correction", "lp");
            leg_cangle_data_sub->AddEntry(graph_cangle_data_corrected, "After correction", "lp");
            
            
            TMultiGraph *mg_cangle_data = new TMultiGraph();
            mg_cangle_data->Add(graph_cangle_data_org);
            mg_cangle_data->Add(graph_cangle_data_corrected);
            
            mg_cangle_data->SetTitle("Data | Reconstructed #theta_{c} angle (proton) ; #theta [degree]; #theta_{C} [rad]");
            mg_cangle_data->Draw("APL");
            mg_cangle_data->GetHistogram()->GetYaxis()->SetRangeUser(0.8,0.85);
            leg_cangle_data_sub->Draw();
            prt_canvasGet("r_graph_cangle_data")->Update();
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
            
            
            /////
            prt_canvasAdd("r_graph_cangle_sim",800,400);
            TLegend *leg_cangle_sim_sub = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_cangle_sim_sub->SetHeader("Sim | Reconstructed #theta_{c} angle (proton)","C");
            leg_cangle_sim_sub->AddEntry(graph_cangle_sim_org, "Before correction", "lp");
            leg_cangle_sim_sub->AddEntry(graph_cangle_sim_corrected, "After correction", "lp");
            
            TMultiGraph *mg_cangle_sim = new TMultiGraph();
            mg_cangle_sim->Add(graph_cangle_sim_org);
            mg_cangle_sim->Add(graph_cangle_sim_corrected);
            
            mg_cangle_sim->SetTitle("Sim | Reconstructed #theta_{c} angle (proton) ; #theta [degree]; #theta_{C} [rad]");
            mg_cangle_sim->Draw("APL");
            mg_cangle_sim->GetHistogram()->GetYaxis()->SetRangeUser(0.8,0.85);
            leg_cangle_sim_sub->Draw();
            prt_canvasGet("r_graph_cangle_sim")->Update();
            TLine *line_ch_p2 = new TLine(0,0,0,1000);
            line_ch_p2->SetY1(fAngleP);
            line_ch_p2->SetY2(fAngleP);
            line_ch_p2->SetX1(gPad->GetUxmin());
            line_ch_p2->SetX2(gPad->GetUxmax());
            line_ch_p2->SetLineColor(kRed);
            line_ch_p2->Draw();
            
            TLine *line_ch_pi2 = new TLine(0,0,0,1000);
            line_ch_pi2->SetY1(fAnglePi);
            line_ch_pi2->SetY2(fAnglePi);
            line_ch_pi2->SetX1(gPad->GetUxmin());
            line_ch_pi2->SetX2(gPad->GetUxmax());
            line_ch_pi2->SetLineColor(kBlue);
            line_ch_pi2->Draw();
            
        }
        
        //////////////
        //// spr /////
        //////////////
        if(bool_partIII){
            graph_spr_data_org->SetLineColor(kBlack);
            graph_spr_data_org->SetMarkerColor(kBlack);
            graph_spr_sim_org->SetLineColor(kRed);
            graph_spr_sim_org->SetMarkerColor(kRed);
            
            graph_spr_data_sub->SetLineColor(kBlue);
            graph_spr_data_sub->SetMarkerColor(kBlue);
            graph_spr_sim_true->SetLineColor(kGreen);
            graph_spr_sim_true->SetMarkerColor(kGreen);
            
            graph_spr_data_corrected->SetLineColor(kMagenta);
            graph_spr_data_corrected->SetMarkerColor(kMagenta);
            graph_spr_sim_corrected->SetLineColor(kBlack);
            graph_spr_sim_corrected->SetMarkerColor(kBlack);
            graph_spr_sim_corrected->SetLineStyle(2);
            
            TMultiGraph *mg_spr_1 = new TMultiGraph();
            TMultiGraph *mg_spr_2 = new TMultiGraph();
            TMultiGraph *mg_spr_3 = new TMultiGraph();
            
            mg_spr_1->Add(graph_spr_data_org);
            mg_spr_1->Add(graph_spr_sim_org);
            TLegend *leg_spr1 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_spr1->SetHeader("SPR (proton)","C");
            leg_spr1->AddEntry(graph_spr_data_org, "Data without corrections ", "lp");
            leg_spr1->AddEntry(graph_spr_sim_org, "Sim without corrections", "lp");
            
            mg_spr_2->Add(graph_spr_data_org);
            mg_spr_2->Add(graph_spr_data_corrected);
            TLegend *leg_spr2 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_spr2->SetHeader("SPR (proton)","C");
            leg_spr2->AddEntry(graph_spr_data_org, "Data without corrections ", "lp");
            leg_spr2->AddEntry(graph_spr_data_corrected, "Data with corrections", "lp");
            
            mg_spr_3->Add(graph_spr_data_corrected);
            mg_spr_3->Add(graph_spr_sim_corrected);
            TLegend *leg_spr3 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_spr3->SetHeader("SPR (proton)","C");
            leg_spr3->AddEntry(graph_spr_data_corrected, "Data with corrections ", "lp");
            leg_spr3->AddEntry(graph_spr_sim_corrected, "Sim with corrections", "lp");
            
            prt_canvasAdd("r_spr_1",800,400);
            graph_spr_data_org->Draw();
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
        
        if(bool_partIII) {
            
            
            graph_yield_DIRC_wt_data->SetLineColor(kBlack);
            graph_yield_DIRC_wt_sim->SetLineColor(kMagenta);
            graph_yield_DIRC_wt_sim->SetLineStyle(2);
            
            
            prt_canvasAdd("r_photon_yield_all",800,400);
            TLegend *leg_yield = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_yield->SetHeader("DIRC photon solutions (proton)","C");
            leg_yield->AddEntry(graph_yield_DIRC_wo_sim, "Sim without cuts", "lp");
            leg_yield->AddEntry(graph_yield_DIRC_wtc_sim, "Sim with #theta_{C} and time cuts", "lp");
            leg_yield->AddEntry(graph_yield_DIRC_wo_data, "Data without cuts", "lp");
            leg_yield->AddEntry(graph_yield_DIRC_wtc_data, "Data with #theta_{C} and time cuts", "lp");
            
            TMultiGraph *mg_yield = new TMultiGraph();
            mg_yield->Add(graph_yield_DIRC_wo_sim);
            mg_yield->Add(graph_yield_DIRC_wtc_sim);
            mg_yield->Add(graph_yield_DIRC_wo_data);
            mg_yield->Add(graph_yield_DIRC_wtc_data);
            
            mg_yield->SetTitle("photon yield ;#theta [degree]; count [#]");
            mg_yield->Draw("APL");
            leg_yield->Draw();
            
            TMultiGraph *mg_yield_1 = new TMultiGraph();
            TMultiGraph *mg_yield_2 = new TMultiGraph();
            TMultiGraph *mg_yield_3 = new TMultiGraph();
            
            mg_yield_1->Add(graph_yield_DIRC_wo_sim);
            mg_yield_1->Add(graph_yield_DIRC_wo_data);
            TLegend *leg_yield1 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_yield1->SetHeader("DIRC photon solutions (proton)","C");
            leg_yield1->AddEntry(graph_yield_DIRC_wo_sim, "Sim without cuts ", "lp");
            leg_yield1->AddEntry(graph_yield_DIRC_wo_data, "Data without cuts", "lp");
            
            mg_yield_2->Add(graph_yield_DIRC_wt_sim);
            mg_yield_2->Add(graph_yield_DIRC_wt_data);
            TLegend *leg_yield2 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_yield2->SetHeader("DIRC photon solutions (proton)","C");
            leg_yield2->AddEntry(graph_yield_DIRC_wt_sim, "Sim with time cut ", "lp");
            leg_yield2->AddEntry(graph_yield_DIRC_wt_data, "Data with time cut", "lp");
            
            mg_yield_3->Add(graph_yield_DIRC_wtc_sim);
            mg_yield_3->Add(graph_yield_DIRC_wtc_data);
            TLegend *leg_yield3 = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            leg_yield3->SetHeader("DIRC photon solutions (proton)","C");
            leg_yield3->AddEntry(graph_yield_DIRC_wtc_sim, "Sim with #theta_{C} and time cuts ", "lp");
            leg_yield3->AddEntry(graph_yield_DIRC_wtc_data, "Data with #theta_{C} and time cuts", "lp");
            
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
        if(bool_partIII) {
            prt_canvasAdd("r_calc_mom",800,400);
            TMultiGraph *mg_calc_mom = new TMultiGraph();
            mg_calc_mom->Add(calc_mom);
            mg_calc_mom->SetTitle(" Calculated momentum ;#theta [degree]; momentum [GeV/c]");
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
            mg_calc_tof1tof2_distance->SetTitle(" Calculated distance ;#theta [degree]; tof2 tof1 distance [m]");
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
    
    
    //    // histo
    //    delete p_cherenkov_sim_org;
    //    delete p_cherenkov_data_org;
    //    delete p_cherenkov_sim_corrected;
    //    delete p_cherenkov_data_corrected;
    //    delete p_cherenkov_data_copy;
    //    delete p_cherenkov_sim_copy;
    //    delete p_cherenkov_mc_same_path;
    //    delete p_cherenkov_bg_sim;
    //    delete p_yield_wo_sim;
    //    delete p_yield_wt_sim;
    //    delete p_yield_wtc_sim;
    //    delete p_yield_true_sim;
    //    delete p_yield_wo_data;
    //    delete p_yield_wt_data;
    //    delete p_yield_wtc_data;
    //    delete p_diff_time_sim;
    //    delete p_diff_time_mctruth;
    //    delete p_diff_time_bg_sim;
    //    delete p_diff_time_data;
    //    delete p_photon_time_sim;
    //    delete p_photon_time_data;
    //    delete p_photon_time_sim_calc;
    //    delete p_photon_time_data_calc;
    
    
    /*
     //graphs
     delete graph_spr_sim_org;
     delete graph_spr_data_sub;
     delete graph_spr_sim_true;
     delete graph_spr_data_corrected;
     delete graph_spr_sim_corrected;
     delete graph_spr_data_org;
     
     delete graph_diff_true;
     delete graph_diff_data;
     delete graph_diff_sim;
     delete graph_diff_data_mean;
     delete graph_diff_sim_mean;
     delete graph_diff_true_mean;
     delete graph_diff_true_sim;
     delete graph_cangle_sim_org;
     delete graph_cangle_sim_true;
     delete graph_cangle_data_corrected;
     delete graph_cangle_sim_corrected;
     delete graph_cangle_data_org;
     delete graph_yield_DIRC_wo_sim;
     delete graph_yield_DIRC_wt_sim;
     delete graph_yield_DIRC_wtc_sim;
     delete graph_yield_DIRC_true_sim;
     delete graph_yield_DIRC_wo_data;
     delete graph_yield_DIRC_wt_data;
     delete graph_yield_DIRC_wtc_data;
     
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
void histo_style_cherenkov(TH1F *histo1, TH1F *histo2, TH1F *histo3, TH1F *histo4, TH1F *histo5){
    //data copy
    histo1->SetLineColor(kBlue);
    histo1->SetLineStyle(1);
    histo1->GetXaxis()->SetTitle("#theta_{C} [rad]");
    histo1->GetYaxis()->SetTitle("entries [#]");
    histo1->GetXaxis()->SetTitleSize(0.05);
    histo1->GetYaxis()->SetTitleSize(0.05);
    histo1->GetXaxis()->SetTitleOffset(0.9);
    histo1->GetYaxis()->SetTitleOffset(1.0);
    histo1->SetFillColor(kBlue);
    histo1->SetFillStyle(3003);
    histo1->SetMarkerStyle(8);
    histo1->SetMarkerColor(kBlue);
    //sub
    histo2->SetName("estimated data signal");
    histo2->SetLineColor(kRed);
    histo2->SetLineStyle(1);
    histo2->GetXaxis()->SetTitle("ç");
    histo2->GetYaxis()->SetTitle("entries [#]");
    histo2->GetXaxis()->SetTitleSize(0.05);
    histo2->GetYaxis()->SetTitleSize(0.05);
    histo2->GetXaxis()->SetTitleOffset(0.9);
    histo2->GetYaxis()->SetTitleOffset(1.0);
    histo2->SetFillColor(kRed);
    histo2->SetFillStyle(3004);
    histo2->SetMarkerStyle(8);
    histo2->SetMarkerColor(kRed);
    //BG
    histo3->SetLineColor(kBlack);
    histo3->SetLineStyle(1);
    histo3->GetXaxis()->SetTitle("#theta_{C} [rad]");
    histo3->GetYaxis()->SetTitle("entries [#]");
    histo3->GetXaxis()->SetTitleSize(0.05);
    histo3->GetYaxis()->SetTitleSize(0.05);
    histo3->GetXaxis()->SetTitleOffset(0.9);
    histo3->GetYaxis()->SetTitleOffset(1.0);
    histo3->SetFillColor(kBlack);
    histo3->SetFillStyle(4050);
    histo3->SetFillStyle(3001);
    histo3->SetMarkerStyle(4);
    // true
    histo4->SetLineColor(kMagenta);
    histo4->SetLineStyle(1);
    histo4->GetXaxis()->SetTitle("#theta_{C} [rad]");
    histo4->GetYaxis()->SetTitle("entries [#]");
    histo4->GetXaxis()->SetTitleSize(0.05);
    histo4->GetYaxis()->SetTitleSize(0.05);
    histo4->GetXaxis()->SetTitleOffset(0.9);
    histo4->GetYaxis()->SetTitleOffset(1.0);
    histo4->SetFillColor(kMagenta);
    histo4->SetFillStyle(3001);
    // sim
    histo5->SetName("MC signal");
    histo5->SetLineColor(12);
    histo5->SetLineStyle(1);
    histo5->GetXaxis()->SetTitle("#theta_{C} [rad]");
    histo5->GetYaxis()->SetTitle("entries [#]");
    histo5->GetXaxis()->SetTitleSize(0.05);
    histo5->GetYaxis()->SetTitleSize(0.05);
    histo5->GetXaxis()->SetTitleOffset(0.9);
    histo5->GetYaxis()->SetTitleOffset(1.0);
    histo5->SetFillColor(12);
    histo5->SetFillStyle(3003);
}
void histo_style_time_diff(TH1F *histo1, TH1F *histo2, TH1F *histo3, TH1F *histo4, TH1F *histo5){
    // sim
    histo1->SetName("MC signal");
    histo1->SetLineColor(12);
    histo1->SetLineStyle(1);
    histo1->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
    histo1->GetYaxis()->SetTitle("entries [#]");
    histo1->GetXaxis()->SetTitleSize(0.05);
    histo1->GetYaxis()->SetTitleSize(0.05);
    histo1->GetXaxis()->SetTitleOffset(0.9);
    histo1->GetYaxis()->SetTitleOffset(1.0);
    histo1->SetFillColor(12);
    histo1->SetFillStyle(3003);
    // data
    histo2->SetLineColor(kBlue);
    histo2->SetLineStyle(1);
    histo2->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
    histo2->GetYaxis()->SetTitle("entries [#]");
    histo2->GetXaxis()->SetTitleSize(0.05);
    histo2->GetYaxis()->SetTitleSize(0.05);
    histo2->GetXaxis()->SetTitleOffset(0.9);
    histo2->GetYaxis()->SetTitleOffset(1.0);
    histo2->SetFillColor(kBlue);
    histo2->SetFillStyle(3003);
    histo2->SetMarkerStyle(8);
    histo2->SetMarkerSize(0.5);
    histo2->SetMarkerColor(kBlue);
    // sub
    histo3->SetName("estimated data signal");
    histo3->SetLineColor(kRed);
    histo3->SetLineStyle(1);
    histo3->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
    histo3->GetYaxis()->SetTitle("entries [#]");
    histo3->GetXaxis()->SetTitleSize(0.05);
    histo3->GetYaxis()->SetTitleSize(0.05);
    histo3->GetXaxis()->SetTitleOffset(0.9);
    histo3->GetYaxis()->SetTitleOffset(1.0);
    histo3->SetFillColor(kRed);
    histo3->SetFillStyle(3004);
    histo3->SetMarkerStyle(8);
    histo3->SetMarkerColor(kRed);
    // true
    histo4->SetLineColor(kMagenta);
    histo4->SetLineStyle(1);
    histo4->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
    histo4->GetYaxis()->SetTitle("entries [#]");
    histo4->GetXaxis()->SetTitleSize(0.05);
    histo4->GetYaxis()->SetTitleSize(0.05);
    histo4->GetXaxis()->SetTitleOffset(0.9);
    histo4->GetYaxis()->SetTitleOffset(1.0);
    histo4->SetFillColor(kMagenta);
    histo4->SetFillStyle(3001);
    // bg
    histo5->SetLineColor(kBlack);
    histo5->SetLineStyle(1);
    histo5->GetXaxis()->SetTitle("t_{calc}-t_{measured} [ns]");
    histo5->GetYaxis()->SetTitle("entries [#]");
    histo5->GetXaxis()->SetTitleSize(0.05);
    histo5->GetYaxis()->SetTitleSize(0.05);
    histo5->GetXaxis()->SetTitleOffset(0.9);
    histo5->GetYaxis()->SetTitleOffset(1.0);
    histo5->SetFillColor(kBlack);
    histo5->SetFillStyle(3001);
    histo5->SetMarkerStyle(4);
}



void histo_style_photon_time(TH1F *histo1, TH1F *histo2, TH1F *histo3, TH1F *histo4){
    //sim
    histo1->SetLineColor(kBlue);
    histo1->SetLineStyle(1);
    histo1->GetXaxis()->SetTitle("t_{measured} [ns]");
    histo1->GetYaxis()->SetTitle("entries [#]");
    histo1->GetXaxis()->SetTitleSize(0.05);
    histo1->GetYaxis()->SetTitleSize(0.05);
    histo1->GetXaxis()->SetTitleOffset(0.9);
    histo1->GetYaxis()->SetTitleOffset(1.0);
    histo1->SetFillColor(kBlue);
    histo1->SetFillStyle(3003);
    //data
    histo2->SetLineColor(kRed);
    histo2->SetLineStyle(1);
    histo2->GetXaxis()->SetTitle("t_{measured} [ns]");
    histo2->GetYaxis()->SetTitle("entries [#]");
    histo2->GetXaxis()->SetTitleSize(0.05);
    histo2->GetYaxis()->SetTitleSize(0.05);
    histo2->GetXaxis()->SetTitleOffset(0.9);
    histo2->GetYaxis()->SetTitleOffset(1.0);
    histo2->SetFillColor(kRed);
    histo2->SetFillStyle(3001);
    histo2->SetMarkerStyle(21);
    histo2->SetMarkerColor(kRed);
    histo2->SetMarkerSize(0.5);
    // calc sim
    histo3->SetLineColor(kBlue);
    histo3->SetLineStyle(1);
    histo3->GetXaxis()->SetTitle("t_{calc} [ns]");
    histo3->GetYaxis()->SetTitle("entries [#]");
    histo3->GetXaxis()->SetTitleSize(0.05);
    histo3->GetYaxis()->SetTitleSize(0.05);
    histo3->GetXaxis()->SetTitleOffset(0.9);
    histo3->GetYaxis()->SetTitleOffset(1.0);
    histo3->SetFillColor(kBlue);
    histo3->SetFillStyle(3003);
    // calc data
    histo4->SetLineColor(kRed);
    histo4->SetLineStyle(1);
    histo4->GetXaxis()->SetTitle("t_{calc} [ns]");
    histo4->GetYaxis()->SetTitle("entries [#]");
    histo4->GetXaxis()->SetTitleSize(0.05);
    histo4->GetYaxis()->SetTitleSize(0.05);
    histo4->GetXaxis()->SetTitleOffset(0.9);
    histo4->GetYaxis()->SetTitleOffset(1.0);
    histo4->SetFillColor(kRed);
    histo4->SetFillStyle(3001);
    histo4->SetMarkerStyle(21);
    histo4->SetMarkerColor(kRed);
    histo4->SetMarkerSize(0.5);
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

void histo_style_match(TH1F *p_cherenkov_sim_org, TH1F *p_cherenkov_data_org) {
    p_cherenkov_sim_org->SetName("MC signal");
    p_cherenkov_sim_org->SetLineColor(kMagenta);
    p_cherenkov_sim_org->SetLineStyle(1);
    p_cherenkov_sim_org->GetXaxis()->SetTitle("#theta_{C} [rad]");
    p_cherenkov_sim_org->GetYaxis()->SetTitle("entries [#]");
    p_cherenkov_sim_org->GetXaxis()->SetTitleSize(0.05);
    p_cherenkov_sim_org->GetYaxis()->SetTitleSize(0.05);
    p_cherenkov_sim_org->GetXaxis()->SetTitleOffset(0.9);
    p_cherenkov_sim_org->GetYaxis()->SetTitleOffset(1.0);
    p_cherenkov_sim_org->SetFillColor(kMagenta);
    p_cherenkov_sim_org->SetFillStyle(3001);
    
    p_cherenkov_data_org->SetLineColor(kBlack);
    p_cherenkov_data_org->SetLineStyle(1);
    p_cherenkov_data_org->GetXaxis()->SetTitle("#theta_{C} [rad]");
    p_cherenkov_data_org->GetYaxis()->SetTitle("entries [#]");
    p_cherenkov_data_org->GetXaxis()->SetTitleSize(0.05);
    p_cherenkov_data_org->GetYaxis()->SetTitleSize(0.05);
    p_cherenkov_data_org->GetXaxis()->SetTitleOffset(0.9);
    p_cherenkov_data_org->GetYaxis()->SetTitleOffset(1.0);
    p_cherenkov_data_org->SetFillColor(kBlack);
    p_cherenkov_data_org->SetFillStyle(3003);
    p_cherenkov_data_org->SetMarkerStyle(8);
    p_cherenkov_data_org->SetMarkerColor(kBlack);
}

void histo_style_photon_yield(TH1F *histo1, TH1F *histo2, TH1F *histo3){
    
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

std::pair<Double_t, Double_t> FitHisto_m(TH1F *hiso_m_fit_option, Double_t cangle_true, Double_t fit_range)
{
    TF1 *fFit_hiso_m_fit_option = new TF1("fFit_hiso_m_fit_option","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
    fFit_hiso_m_fit_option->SetParameters(100,cangle_true,0.010);
    fFit_hiso_m_fit_option->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
    fFit_hiso_m_fit_option->SetParLimits(0,0.1,1E6);
    fFit_hiso_m_fit_option->SetParLimits(1,cangle_true-fit_range,cangle_true+fit_range);
    fFit_hiso_m_fit_option->SetParLimits(2,0.005,0.014000);
    hiso_m_fit_option->Fit("fFit_hiso_m_fit_option","M","",cangle_true-fit_range,cangle_true+fit_range);
    Double_t cangle_FitHisto_m = fFit_hiso_m_fit_option->GetParameter(1);
    Double_t spr_FitHisto_m = fFit_hiso_m_fit_option->GetParameter(2);
    return std::make_pair(cangle_FitHisto_m, spr_FitHisto_m);
}

std::pair<Double_t, Double_t> FitHisto_0(TH1F *hiso_0_fit_option, Double_t cangle_true, Double_t fit_range)
{
    TF1 *fFit_hiso_0_fit_option = new TF1("fFit_hiso_0_fit_option","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
    fFit_hiso_0_fit_option->SetParameters(100,cangle_true,0.010);
    fFit_hiso_0_fit_option->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
    fFit_hiso_0_fit_option->SetParLimits(0,0.1,1E6);
    fFit_hiso_0_fit_option->SetParLimits(1,cangle_true-fit_range,cangle_true+fit_range);
    fFit_hiso_0_fit_option->SetParLimits(2,0.005,0.014000);
    hiso_0_fit_option->Fit("fFit_hiso_0_fit_option","M","",cangle_true-fit_range,cangle_true+fit_range);
    Double_t cangle_FitHisto_0 = fFit_hiso_0_fit_option->GetParameter(1);
    Double_t spr_FitHisto_0 = fFit_hiso_0_fit_option->GetParameter(2);
    return std::make_pair(cangle_FitHisto_0, spr_FitHisto_0);
}



