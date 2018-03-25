#include "/Users/ahmed/dirc/prttools/prttools.C"

//root plot_beam_correction.C'("/data.local/beam_correction/db_20_3sph_t1_*_t2_*_proton_data_spr.root", 20)'
//root plot_beam_correction.C'("/Users/ahmed/dirc/beam_correction/beam_correction_20/opt_20_3sph_t1_*_t2_*_pi_data_pi_correction_spr.root", 20)'

// file existance
bool exists_test (const std::string& name);
void plot_beam_correction(/*TString inFile = "r_spr.root", Int_t angle= 20*/) {
    
    prt_savepath="opt";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    prt_setRootPalette(1);
    //TChain ch("dirc");
    //ch.Add(inFile);
    //Double_t cangle,spr,trr,nph,par1,par2,par3,par4,par5,par6,test1,test2,theta,phi;
    Int_t nf = 40;
    //TH2F * twoD_nph =  new TH2F("twoD_nph",";#Delta#Theta[mrad]; #Delta#Phi[mrad]", nf, -20, 20, nf, -10,10);
    //TH1F * timeCut =  new TH1F("timecut",";Time Cut [ns];Count[#]",nf, 0.1, 6);
    //TH1F * chCut =  new TH1F("chCut",";Cherenkov Angle Cut[rad];Count[#]",nf, 0.001,0.06);
    
    
    TH2F * twoD_nph =  new TH2F("twoD_nph",";#Delta#Theta[mrad]; #Delta#Phi[mrad]", nf, -20, 20, nf, -10,10);
    Int_t angle;
    TFile *ffile;
    TH1F *chere;
    for (int i=20; i<=150; i+=10) {
        for (int j= -200; j<=200; j+=10) {
            for (int k=-100; k<=100; k+=5) {
               // if (i == 40){
                     angle=i;
                    Double_t jj = (Float_t) j/10000.0;
                    Double_t kk = (Float_t) k/10000.0;
                    TString jj_string = Form("t1_%.3f", jj);
                    TString kk_string = Form("_t2_%.4f", kk);
                    cout<< "enter the if condition"<<endl;
                    TString cherenkov_data_path = Form("/Users/ahmed/dirc/beam_correction/beam_correction_%d/opt_%d_3sph_"+jj_string+kk_string+"_p_data_wo_correction_spr.root",i, i);
                    //cout<<"cherenkov_data_path= " <<cherenkov_data_path<<endl;
                    
                    string path = (string)cherenkov_data_path;
                    cout<<"exists_test(path)" <<exists_test(path)<<endl;
                    if (!exists_test(path)) continue;
                    cout<<"path= " <<path<<endl;
                    std::cout<<"############"<< " no problem 1 " <<std::endl;
                    ffile  = new TFile(cherenkov_data_path, "READ");
                    chere=(TH1F*)ffile->Get("fHist_correction");
                    
                    
                    
                    
                    
                    TF1 *fFit = new TF1("fFit","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
                    TSpectrum *fSpect= new TSpectrum(10);
                    Double_t cangle=0;
                    Double_t spr=0;
                    //gROOT->SetBatch(1);
                    Int_t nfound = fSpect->Search(chere,1,"",0.9); //0.6
                    if(nfound>0) cangle = fSpect->GetPositionX()[0];
                    cangle =  chere->GetXaxis()->GetBinCenter(chere->GetMaximumBin());
                    if(cangle>0.85) cangle=0.82;
                    fFit->SetParameters(100,cangle,0.010);
                    fFit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
                    fFit->SetParLimits(0,0.1,1E6);
                    fFit->SetParLimits(1,cangle-0.04,cangle+0.04);
                    fFit->SetParLimits(2,0.005,0.018); // changed 0.014
                    chere->Fit("fFit","M+","",cangle-0.1,cangle+0.1);
                    //p_cherenkov_data_copy->Fit("fFit","0","",cangle-0.06,cangle+0.06);
                    //chere->Fit("fFit","R");
                    Double_t chi = fFit->GetChisquare()/fFit->GetNDF();
                    cangle = fFit->GetParameter(1);
                    spr = fFit->GetParameter(2);
                    Double_t cangle_minus_5_sgma = cangle-5*spr;
                    Double_t cangle_plus_5_sgma = cangle+5*spr;
                    Double_t cangle_minus_3_sgma = cangle-3*spr;
                    Double_t cangle_plus_3_sgma = cangle+3*spr;
                    Double_t r_min = cangle-8*spr;
                    Double_t r_max = cangle+8*spr;
                    Double_t sumundercurve = fFit->Integral(cangle_minus_3_sgma,cangle_plus_3_sgma);
                    
                    twoD_nph->Fill(jj*1000, kk*1000, spr*1000/14.0);
                    //std::cout<<"############"<< " no problem     "<<jj*1000<<"        "<<kk*1000 <<std::endl;

                    
                    /*
                    TString nid = Form("_%2.0d", angle);
                    prt_canvasAdd("r_beam_correction_pi_pi"+nid+jj_string+kk_string,800,400);
                    chere->SetStats(0);
                    chere-> SetTitle(Form("SPR #pi data with mcp by mcp #theta_{c} correction for #pi %d",angle) );
                    //prt_canvasSave(2,0);
                    //prt_canvasDel("*");
                    */
                    
                    ffile->Close();
                    delete ffile;
               // }
            }
        }
    }
    
    
    
    TString nid = Form("_%2.0d", angle);
    prt_canvasAdd("r_beam_correction_p_wo"+nid,800,400);
    twoD_nph->SetStats(0);
    //twoD_nph-> SetTitle(Form("SPR #pi data with mcp by mcp #theta_{c} correction for P %d",angle) );
    //twoD_nph-> SetTitle("SPR #pi data with mcp by mcp #theta_{c} correction for P" );
    twoD_nph-> SetTitle("SPR P data without mcp by mcp #theta_{c} correction" );
    twoD_nph-> Draw("colz");
    
    
   prt_canvasSave(2,0);
   prt_canvasDel("*");
    

    
    
    /*
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
     for (Int_t i = 0; i < nent; i++) {
     ch.GetEvent(i);
     //if (spr< 6.0||spr> 11.3) spr=9.8 ;
     twoD_nph->Fill(test1*1000, test2*1000, spr);
     }
     
     TString nid = Form("_%2.0d", angle);
     prt_canvasAdd("r_beam_correction_pi_p"+nid,800,400);
     twoD_nph->SetStats(0);
     twoD_nph-> SetTitle(Form("SPR #pi data with mcp by mcp #theta_{c} correction for #p %d",angle) );
     twoD_nph-> Draw("colz");
     
     
     prt_canvasAdd("r_chere"+nid,800,400);
     chere->SetStats(0);
     chere-> SetTitle(Form("SPR #pi data with mcp by mcp #theta_{c} correction for #p %d",angle) );
     chere-> Draw();
     */
    
    
    
}
//////////////////////////
// check file existance //
//////////////////////////

bool exists_test (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}
