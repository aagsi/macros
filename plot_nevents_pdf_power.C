#include "/u/aali/dirc/prttools/prttools.C"
//root plot_nevents_pdf_power.C'("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/sim/332/pdf_neventTest_*_20_sph_separation.root", 20)'

void plot_nevents_pdf_power(TString inFile = "r_spr.root", Int_t angle= 20) {
    int counter =0 ;
    const Int_t n = 21;
    Double_t x[n];

    TGraph *power_org = new TGraph();
    //gStyle->SetPalette(kGreenPink);
    //gStyle->SetPalette(kRainBow);
    prt_savepath="nevents_pdf";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    TChain ch("dirc");
    ch.Add(inFile);
    Double_t cangle,spr,trr,nph,par1,par2,par3,par4,par5,par6,test1,test2,theta,phi, separation;
    Int_t openChCorr;
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
    ch.SetBranchAddress("openChCorr",&openChCorr);

    Int_t nent = ch.GetEntries();
    std::cout<<"# entries  "<< nent <<std::endl;
    std::cout<<"# inFile  "<< inFile <<std::endl;
    //std::cout<<"infor  "<< ch.GetTree()->GetTitle()<<std::endl;
    for (Int_t i = 0; i < nent; i++) {
        ch.GetEvent(i);

        //twoD_nph->Fill(test1*1000, test2*1000, separation);
	//std::cout<<"# i  "<< i <<std::endl;
        //std::cout<<"# counter  "<< counter <<std::endl;
        //std::cout<<"# openChCorr  "<< openChCorr <<std::endl;
	//std::cout<<"# separation  "<< separation <<std::endl;
        if (separation>0){
	power_org->SetPoint(counter,openChCorr,separation);
	counter++;
	}


    }

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



    prt_canvasAdd("r_separation_nevent_pdf",800,400);
    TLegend *leg_separation = new TLegend( 0.607769, 0.614973, 0.887218 ,0.868984);
    //leg_separation->SetHeader("separation power ","C");
    leg_separation->SetFillColor(0);
    leg_separation->AddEntry(power_org, "PDF Events number", "lp");

    TMultiGraph *mg_separation = new TMultiGraph();
    //mg_separation->Add(power_org);
    mg_separation->Add(power_org);
    mg_separation->SetTitle(" separation power geometrical reconstruction PDF method; number of events [#]; separation [s.d.]");
    mg_separation->Draw("APL");
    leg_separation->Draw();


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









