#include "/u/aali/dirc/prttools/prttools.C"
//root plot_beam_correction_power.C'("/data.local/beam_correction/beam_correction/opt_20_3sph_t1_*_t2_*_data_wo_correction_separation.root", 20)'

void plot_beam_correction_power(TString inFile = "r_spr.root", Int_t angle= 20) {


//gStyle->SetPalette(kGreenPink);
gStyle->SetPalette(kRainBow);
    prt_savepath="opt";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    TChain ch("dirc");
    ch.Add(inFile);
    Double_t cangle,spr,trr,nph,par1,par2,par3,par4,par5,par6,test1,test2,theta,phi, separation;
    Int_t nf = 40;
    TH2F * twoD_nph =  new TH2F("twoD_nph",";#Delta#Theta[mrad]; #Delta#Phi[mrad]", nf, -20, 20, nf, -10,10);



    TH1F * timeCut =  new TH1F("timecut",";Time Cut [ns];Count[#]",nf, 0.1, 6);
    TH1F * chCut =  new TH1F("chCut",";Cherenkov Angle Cut[rad];Count[#]",nf, 0.001,0.06);

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
    ch.SetBranchAddress("separation",&separation);

  

    Int_t nent = ch.GetEntries();
    std::cout<<"# entries  "<< nent <<std::endl;
    std::cout<<"# inFile  "<< inFile <<std::endl;
    //std::cout<<"infor  "<< ch.GetTree()->GetTitle()<<std::endl;
    for (Int_t i = 0; i < nent; i++) {
        ch.GetEvent(i);

        twoD_nph->Fill(test1*1000, test2*1000, separation);
         //std::cout<<"# par3  "<< par3 <<std::endl;
         //std::cout<<"# test2  "<< test2 <<std::endl;


    }
            prt_canvasAdd(Form("beam_correction_power_%d", angle),800,400);
            twoD_nph->SetStats(0);
            twoD_nph-> SetTitle(Form("p/#pi Separation power %d",angle) );
            twoD_nph->SetMarkerSize(1.8);
            gStyle->SetPaintTextFormat("4.1f");
            //twoD_nph->Draw("colztext");
            twoD_nph->Draw("colz");
    
    
     prt_canvasSave(2,0);
     prt_canvasDel("*");
    
    /*
    TCanvas* c66 = new TCanvas("c66","c66",800,500);
    twoD_nph->SetStats(0);
    twoD_nph-> SetTitle(Form("p/#pi Separation power %d",angle) );
    twoD_nph-> Draw("colz");
    gSystem->ProcessEvents();
    TImage *img66 = TImage::Create();
    img66->FromPad(c66);
    img66->WriteImage(Form("beam_correction_power_%d.png", angle) );
*/

}







