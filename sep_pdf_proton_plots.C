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
#include "/u/aali/dirc/prttools/prttools.C"
//#include "/Users/ahmed/dirc/prttools/prttools.C"
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


Double_t  prtangle;
// variables for tree
TTree *dirc;
Double_t kk;
Bool_t reject, reject_sd;

////////////////////
// proto types//////
////////////////////
// fitting functions
// file existance
bool exists_test (const std::string& name);

Double_t momentum=7.0;
Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
Double_t fAngleP = acos(sqrt(momentum*momentum+ mass[4]*mass[4])/momentum/1.4738)-0.00;
Double_t fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00; //-0.0014 for 160 25deg


////////////////////
// function   //////
////////////////////
void sep_pdf_proton_plots() {

    const Int_t n = 14;
    Double_t x[n];
    int counter =0 ;
    TLegend* legend = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
    legend->SetHeader("cherenkov angle","C"); // option "C" allows to center the header
    prt_savepath="proton_pdf";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    //Int_t nf = 30;
    //TH2F * photon_yield_map =  new TH2F("photon_yield_map",";Time Cut [ns];cherenkov angle cut[mrad]", nf, 0.1, 6, nf, 0.001,0.06);
    TFile *ffile_sim, *ffile_data;

    TGraph *power_org = new TGraph();
    TGraph *power_sim = new TGraph();
    TGraph *power_modle = new TGraph();
    TGraph *power_pdf = new TGraph();
    TGraph *power_2BarReflpdf = new TGraph();
    TGraph *calc_mom = new TGraph();
    TGraph *calc_e_mom = new TGraph();
    TGraph *calc_tof1tof2_distance = new TGraph();
    TGraph *calc_e_tof1tof2_distance = new TGraph();
    gStyle->SetOptFit(1);
    gStyle->SetOptStat(0);
    for (int i=20; i<=150; i+=10) {
        for (int j=1; j<=60; j+=2) {
            for (int k=1; k<=60; k+=2) {
                //cout<< "i, j, k ="<<"   "<< i<<"  "<<j<<" "<<k<<endl;
                Double_t jj = (Float_t) j/10.0 ;
                kk = (Float_t) k/1000.0;
                TString jj_string = Form("_timecut_%.1f", jj);
                TString kk_string = Form("_chcut_%.3f", kk);
                /////////////////////////////////
                /////////////////////////////////
                if (/*i==150 &&*/ j == 51 && k== 51) {
                    prtangle= i;
                    TString nid = Form("_%2.0d", i);
                    cout<< "enter the if condition"<<endl;

                    // proton
                    TString cherenkov_data_path = Form("/u/aali/work/%d_sph_proton_data_spr.root", i);
                    //TString cherenkov_sim_path = Form("/u/aali/work/reco_proton_bar_3lsph_grease_theta_%d_sim_spr.root", i);
                    TString cherenkov_sim_path = Form("/u/aali/work/%d_sph_proton_sim_spr.root", i);
                    // pi
                    //TString cherenkov_data_path = Form("/u/aali/work/%d_sph_pi_data_spr.root", i);
                    //TString cherenkov_sim_path = Form("/u/aali/work/reco_pi_bar_3lsph_grease_theta_%d_sim_spr.root", i);

                    // separation
                    TString separation_data_path = Form("/u/aali/work/new/sep_old_%d_sph_data_separation.root", i);
                    TString separation_pdf_data_path = Form("/u/aali/work/new/sep_4BarReflpdf_%d_sph_data_separation.root", i);

                    cout<<"cherenkov_sim_path= " <<cherenkov_sim_path<<endl;
                    cout<<"cherenkov_data_path= " <<cherenkov_data_path<<endl;
                    cout<<"separation_data_path= " <<separation_data_path<<endl;
                    cout<<"separation_pdf_data_path= " <<separation_pdf_data_path<<endl;


                    string path_sim = (string)cherenkov_sim_path;
                    string path_data = (string)cherenkov_data_path;
                    string path_data_separation = (string)separation_data_path;
                    string path_data_separation_pdf = (string)separation_pdf_data_path;

                    cout<<"exists_test(path_sim)" <<exists_test(path_sim)<<endl;
                    cout<<"exists_test(path_data)" <<exists_test(path_data)<<endl;
                    cout<<"exists_test(separation path_data)" <<exists_test(path_data_separation)<<endl;
                    cout<<"exists_test(separation path_data shift pdf)" <<exists_test(path_data_separation_pdf)<<endl;


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

                    ////////////////
                    // READ Tree ///
                    ////////////////
                    Double_t separation(-9), separation_org(-1),separation_pdf(-1);
                    Double_t separation1(-9);

                    TChain ch_pdf("dirc");
                    TChain ch("dirc");


                    ch_pdf.Add(separation_pdf_data_path);
                    ch.Add(separation_data_path);
                    

                    ch.SetBranchAddress("separation",&separation);
                    ch_pdf.SetBranchAddress("separation",&separation1);

                
                    Int_t nent = ch.GetEntries();
                    Int_t nent_pdf = ch_pdf.GetEntries();
 
                    ch.GetEvent(nent-1);
                    separation_org=separation;
                    std::cout<<"############  separation = "<< separation  <<std::endl;
                    

                    
                    ch_pdf.GetEvent(nent_pdf-1);
                    separation_pdf=separation1;
                    std::cout<<"############  separation_pdf = "<< separation_pdf  <<std::endl;
                    
                    


                    ///////////////////
                    ///// part II /////
                    ///////////////////
                    if(true) {
                        //////////////////
                        //  Fill graph ///
                        /////////////////
                        x[counter]=i;
                        power_org->SetPoint(counter,i,separation_org);
                        power_pdf->SetPoint(counter,i,separation_pdf);

                        counter++;
                    }
                }
            }
        }
    }
    ///////////////////
    ///// part III ////
    ///////////////////
    if(true) {
        power_org->Sort();
        power_org->SetLineColor(kBlack);
        power_org->SetMarkerColor(kBlack);
        power_org->SetMarkerStyle(21);
        power_org->SetMarkerSize(0.7);
        power_org->SetName("separation");
        power_org->GetYaxis()->SetRangeUser(0,20);
        power_org->GetYaxis()->SetTitle("Separation [s.d]");
        power_org->GetXaxis()->SetLabelSize(0.05);
        power_org->GetXaxis()->SetTitleSize(0.06);
        power_org->GetXaxis()->SetTitleOffset(0.84);


        power_pdf->Sort();
        power_pdf->SetLineColor(kBlue);
        power_pdf->SetMarkerColor(kBlue);
        power_pdf->SetMarkerStyle(21);
        power_pdf->SetMarkerSize(0.7);
        power_pdf->SetName("separation");
        power_pdf->GetYaxis()->SetRangeUser(0,20);
        power_pdf->GetYaxis()->SetTitle("Separation [s.d]");
        power_pdf->GetXaxis()->SetLabelSize(0.05);
        power_pdf->GetXaxis()->SetTitleSize(0.06);
        power_pdf->GetXaxis()->SetTitleOffset(0.84);


        /////////////////////////////
        //// separation power   /////
        /////////////////////////////
        if(true) {
            prt_canvasAdd("r_separation",800,400);
            TLegend *leg_separation = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
            //leg_separation->SetHeader("separation power ","C");
            leg_separation->SetFillColor(0);
            leg_separation->AddEntry(power_org, "Data (standared method w/o #theta_{c} correction)", "lp");
            leg_separation->AddEntry(power_pdf, " Cherenkov PDF Method", "lp");
            TMultiGraph *mg_separation = new TMultiGraph();
            mg_separation->Add(power_org);
            mg_separation->Add(power_pdf);
            mg_separation->SetTitle(" separation power geometrical reconstruction ;#theta [degree]; separation [s.d.]");
            mg_separation->Draw("APL");
            leg_separation->Draw();

        }
        gStyle->SetOptStat(0);
    }
    //////////////////////////
    // save/delete histogram//
    //////////////////////////
    prt_canvasSave(2,0);
    prt_canvasDel("*");

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

