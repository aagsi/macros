
#include <TMath.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include "/Users/ahmed/dirc/prttools/prttools.C"


//root plot_scan_lut.C'("/Users/ahmed/dirc/prtdirc/prism_lens/opt_xstep74.0_ystep14.0_20_sph_p_data_spr.root", 20)'
void plot_scan_lut(TString inFile = "r_spr.root", Int_t angle= 20) {
    prt_savepath="opt";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    TChain ch("dirc");
    ch.Add(inFile);
    Double_t cangle,spr,trr,nph,par1,par2,par3,par4,par5,par6,test1,test2,theta,phi;
    Int_t nf = 40;
    TH2F * twoD_nph =  new TH2F("twoD_nph",";prism-lens X step [mm];prism-lens Y step[mm]", 28, 59.5, 74, 40, 12.9,17);
    
    
    
    TH1F * timeCut =  new TH1F("timecut",";Time Cut [ns];Count[#]",nf, 0.1, 6);
    TH1F * chCut =  new TH1F("chCut",";Cherenkov Angle Cut[rad];Count[#]",nf, 0.001,0.06);
    
    
    TH1F * hito_xstep =  new TH1F("hito_xstep",";X step [mm];Spr[mrad]",28, 59.5, 74);
    TH1F * hito_ystep =  new TH1F("hito_ystep",";Y step [mm];Spr[mrad]",40, 12.9, 17);
    
    
    ch.SetBranchAddress("spr",&spr);
    ch.SetBranchAddress("trr",&trr);
    ch.SetBranchAddress("nph",&nph);
    ch.SetBranchAddress("cangle",&cangle);
    //  ch.SetBranchAddress("par4",&par4);
    ch.SetBranchAddress("par5",&par5);
    ch.SetBranchAddress("par6",&par6);
    ch.SetBranchAddress("test1",&test1);
    ch.SetBranchAddress("test2",&test2);
    ch.SetBranchAddress("theta",&theta);
    ch.SetBranchAddress("phi",&phi);
    
    
    
    Int_t nent = ch.GetEntries();
    std::cout<<"# entries  "<< nent <<std::endl;
    std::cout<<"# inFile  "<< inFile <<std::endl;
    //std::cout<<"infor  "<< ch.GetTree()->GetTitle()<<std::endl;
    Int_t np = nent;
    TGraph2D *dt = new TGraph2D();
    
    dt->SetTitle("LUT optimization; prism-lens X step [mm];prism-lens Y step[mm]; spr");
    
    for (Int_t i = 0; i < nent; i++) {
        ch.GetEvent(i);
        
        twoD_nph->Fill(test1, test2, spr);
        hito_xstep->Fill(test1, spr);
        hito_ystep->Fill(test2, spr);
        //std::cout<<"# test1  "<< test1 <<std::endl;
        //std::cout<<"# test2  "<< test2 <<" "<<i<<std::endl;
        dt->SetPoint(i,test1, test2, spr);
        
        
        
    }
    
    //prt_canvasAdd("xstep",800,400);
    
    
    TCanvas* c66 = new TCanvas("c66","c66",800,500);
    c66->Divide(2,2);
    c66->cd(1);
    
    gStyle->SetPalette(1);
    dt->Draw("surf1");
    c66->cd(2);
    hito_xstep->Draw("histo");
    
    c66->cd(3);
    hito_ystep->Draw("histo");

    c66->cd(4);
    //dt->Project("xy")->Draw("colz");
    twoD_nph->SetStats(0);
    twoD_nph-> SetTitle(Form("Beam correction (SPR pion tagged data at polar angle %d)",angle) );
    twoD_nph-> Draw("colz");
    
    //~ //    c66->WriteImage("r_time_ch_opt_20.png");
    //~ TCanvas* c77 = new TCanvas("c77","c77",800,500);
    //~ twoD_sumunderCurve->SetStats(0);
    //~ twoD_sumunderCurve-> SetTitle(Form("S+B at polar angle %d",angle) );
    //~ twoD_sumunderCurve-> Draw("colz");
    //~
    //~ TCanvas* c88 = new TCanvas("c88","c88",800,500);
    //~ twoD_sumunderLine->SetStats(0);
    //~ twoD_sumunderLine-> SetTitle(Form("B at polar angle %d",angle) );
    //~ twoD_sumunderLine-> Draw("colz");
    //~
    //~ TCanvas* c99 = new TCanvas("c99","c99",800,500);
    //~ twoD_ratio->SetStats(0);
    //~ twoD_ratio-> SetTitle(Form("S/(S+B) at polar angle %d",angle) );
    //~ twoD_ratio-> Draw("colz");
    //~
    //~ TCanvas* c100 = new TCanvas("c100","c100",800,500);
    //~ twoD_ratioXnph->SetStats(0);
    //~ twoD_ratioXnph-> SetTitle(Form("number of photon solutions * S/(S+B) at polar angle %d",angle) );
    //~ twoD_ratioXnph-> Draw("colz");
    
    //////////////////////////////
    gSystem->ProcessEvents();
    /////////////////////////
    TImage *img66 = TImage::Create();
    //~ TImage *img77 = TImage::Create();
    //~ TImage *img88 = TImage::Create();
    //~ TImage *img99 = TImage::Create();
    //~ TImage *img100 = TImage::Create();
    ////////////////////
    img66->FromPad(c66);
    //~ img77->FromPad(c77);
    //~ img88->FromPad(c88);
    //~ img99->FromPad(c99);
    //~ img100->FromPad(c100);
    //////////////////////
    img66->WriteImage(Form("r_time_ch_opt_%d.png", angle) );
    //~ img77->WriteImage(Form("r_time_ch_opt_SplusB_%d.png", angle) );
    //~ img88->WriteImage(Form("r_time_ch_opt_B_%d.png", angle) );
    //~ img99->WriteImage(Form("r_time_ch_opt_ratio_%d.png", angle) );
    //~ img100->WriteImage(Form("r_time_ch_opt_ratioXnph_%d.png", angle) );
    
    
    //~ TString nid = Form("_%2.0d", angle);
    //~ prt_canvasAdd("r_timeCut"+nid,800,400);
    //~ timeCut-> Draw();
    //~
    //~ prt_canvasAdd("r_chCut"+nid,800,400);
    //~ chCut-> Draw();
    //~
    //~ //prt_waitPrimitive("test"+nid);
    //~ prt_canvasSave(2,0);
    //~ prt_canvasDel("*");
    
    //~ delete twoD_nph;
    //~ delete twoD_sumunderCurve;
    //~ delete twoD_sumunderLine;
    //~ delete twoD_ratio;
    //~ delete twoD_ratioXnph;
    //~
    //~ delete c66;
    //~ delete c77;
    //~ delete c88;
    //~ delete c99;
    //~ delete c100;
    //~
    //~ delete img66;
    //~ delete img77;
    //~ delete img88;
    //~ delete img99;
    //~ delete img100;
    
}







