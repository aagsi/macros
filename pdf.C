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

Bool_t Bool_cherenkov_correction(false), Bool_separation_graph(false), Bool_photonyield_histo(true), Bool_cherenkov_PDF_histo(false);
////////////////////
// function   //////
////////////////////
void pdf() {
    TString nid = Form("_%2.0d", 20);
    prt_savepath="pdf";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    TFile *ffile_data_p, *ffile_data_pi;
    TGraph *fHistCh_graph_p[960], *fHistCh_graph_pi[960];

    TH1F *fHistCh_read_p[960], *fHistCh_read_pi[960], *p_cherenkov_data ;
    TH1F *hist_nph_wo_p,  *hist_nph_wt_p, *hist_nph_wtc_p, *hist_nph_wo_pi,  *hist_nph_wt_pi, *hist_nph_wtc_pi;

    THStack *hs_nph_p, *hs_nph_pi;

    gStyle->SetOptFit(1);
    gStyle->SetOptStat(1);

    ///////////////////////
    // Separation power  //
    ///////////////////////
    int counter =0 ;
    TGraph *power_org = new TGraph();
    TH1F *HistMcp_pi[12], *HistMcp_p[12];
    TFile *ffile_data_pi_spr, *ffile_data_p_spr;
    for (int i=20; i<=150; i+=10) {
        if (Bool_cherenkov_correction== true && i >20) break;
        if (Bool_separation_graph== true && i ==70) continue;
        TString separation_data_path = Form("/u/aali/work/test/%d_sph_data_separation.root", i);
        TString spr_data_p_path = Form("/u/aali/work/test/%d_sph_p_data_spr.root", i);
        TString spr_data_pi_path = Form("/u/aali/work/test/%d_sph_pi_data_spr.root", i);
        cout<<"separation data path= " <<separation_data_path<<endl;
        cout<<"spr data p path= " <<spr_data_p_path<<endl;
        cout<<"spr data pi path= " <<spr_data_pi_path<<endl;
        string path_data_separation = (string)separation_data_path;
        string path_data_p_spr = (string)spr_data_p_path;
        string path_data_pi_spr = (string)spr_data_pi_path;
        cout<<"exists test(separation path data)" <<exists_test(path_data_separation)<<endl;
        cout<<"exists test(spr path data p)" <<exists_test(path_data_p_spr)<<endl;
        cout<<"exists test(spr path data pi)" <<exists_test(path_data_pi_spr)<<endl;
        if (!exists_test(path_data_separation)) continue;
        //if (!exists_test(path_data_p_spr)) continue;
        //if (!exists_test(path_data_pi_spr)) continue;
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
        std::cout<<"############  separation = "<< separation  <<std::endl;
        power_org->SetPoint(counter,i,separation_org);
        counter++;
        if(Bool_cherenkov_correction) {
            ///////////////////////////
            // CH MCP by MCP  Histo ///
            ///////////////////////////
            ffile_data_p_spr  = new TFile(spr_data_p_path, "READ");
            ffile_data_pi_spr  = new TFile(spr_data_pi_path, "READ");
            cout << "MCP by MCP cherenvove angle correction for PION"<< endl;
            for(Int_t mcp=0; mcp<prt_nmcp; mcp++) {
                HistMcp_p[mcp] =(TH1F*)ffile_data_p_spr->Get(Form("fHistMcp_%d",mcp));
            }
            // proton correction
            TF1 *Fit_MCP_p = new TF1("Fit_MCP_p","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
            Fit_MCP_p->SetParameters(100,fAngleP,0.010);
            Fit_MCP_p->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            Fit_MCP_p->SetParLimits(0,0.1,1E6);
            Fit_MCP_p->SetParLimits(1,fAngleP-0.02,fAngleP+0.02);
            Fit_MCP_p->SetParLimits(2,0.004,0.008);
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
                prt_canvasAdd(Form("r_mcp_p_%d_prtangle_%d",mcp,i),800,400);
                if(mcp ==3 && i ==600) {
                    HistMcp_p[mcp]->Fit("Fit_MCP_p","lq","",fAngleP-0.025,fAngleP+0.03);
                }
                else if(mcp ==9 && i ==1200) {
                    HistMcp_p[mcp]->Fit("Fit_MCP_p","lq","",fAngleP-0.025,fAngleP+0.010);
                }
                else if(mcp ==6 && i ==1400) {
                    HistMcp_p[mcp]->Fit("Fit_MCP_p","lq","",fAngleP-0.025,fAngleP+0.03);
                }
                else {
                    HistMcp_p[mcp]->Fit("Fit_MCP_p","lq","",fAngleP-0.025,fAngleP+0.025);
                }
                std::cout<<"if(mcpid=="<< mcp<<") tangle += "<<fAngleP-Fit_MCP_p->GetParameter(1)<<";" <<std::endl;
                HistMcp_p[mcp]->Draw();
                prt_canvasGet(Form("r_mcp_p_%d_prtangle_%d",mcp,i))->Update();
                TLine *lin_ch_p_v = new TLine(0,0,0,1000);
                lin_ch_p_v->SetX1(fAngleP);
                lin_ch_p_v->SetX2(fAngleP);
                lin_ch_p_v->SetY1(gPad->GetUymin());
                lin_ch_p_v->SetY2(gPad->GetUymax());
                lin_ch_p_v->SetLineColor(kRed);
                lin_ch_p_v->Draw();
            }
            cout << "MCP by MCP cherenvove angle correction for PION"<< endl;
            for(Int_t mcp=0; mcp<prt_nmcp; mcp++) {
                HistMcp_pi[mcp] =(TH1F*)ffile_data_pi_spr->Get(Form("fHistMcp_%d",mcp));
            }
            // pi correction
            TF1 *Fit_MCP_pi = new TF1("Fit_MCP_pi","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
            Fit_MCP_pi->SetParameters(100,fAnglePi,0.010);
            Fit_MCP_pi->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
            Fit_MCP_pi->SetParLimits(0,0.1,1E6);
            Fit_MCP_pi->SetParLimits(1,fAnglePi-0.02,fAnglePi+0.02);
            Fit_MCP_pi->SetParLimits(2,0.004,0.008);
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
                prt_canvasAdd(Form("r_mcp_pi_%d_prtangle_%d",mcp,i),800,400);

                if(mcp ==11 && i ==100) {
                    HistMcp_pi[mcp]->Fit("Fit_MCP_pi","lq","",fAnglePi-0.025,fAnglePi+0.03);
                }
                else if(mcp ==9 && i ==1200) {
                    HistMcp_pi[mcp]->Fit("Fit_MCP_pi","lq","",fAnglePi-0.025,fAnglePi+0.010);
                }
                else if(mcp ==6 && i ==1400) {
                    HistMcp_pi[mcp]->Fit("Fit_MCP_pi","lq","",fAnglePi-0.025,fAnglePi+0.03);
                }
                else {
                    HistMcp_pi[mcp]->Fit("Fit_MCP_pi","lq","",fAnglePi-0.025,fAnglePi+0.025);
                }
                std::cout<<"if(mcpid=="<< mcp<<") tangle += "<<fAnglePi-Fit_MCP_pi->GetParameter(1)<<";" <<std::endl;
                HistMcp_pi[mcp]->Draw();
                prt_canvasGet(Form("r_mcp_pi_%d_prtangle_%d",mcp,i))->Update();
                TLine *lin_ch_pi_v = new TLine(0,0,0,1000);
                lin_ch_pi_v->SetX1(fAnglePi);
                lin_ch_pi_v->SetX2(fAnglePi);
                lin_ch_pi_v->SetY1(gPad->GetUymin());
                lin_ch_pi_v->SetY2(gPad->GetUymax());
                lin_ch_pi_v->SetLineColor(kBlue);
                lin_ch_pi_v->Draw();
            }
        }

    }
    GraphStyle(power_org);
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






















