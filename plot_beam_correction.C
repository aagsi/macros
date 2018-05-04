//#include "/Users/ahmed/dirc/prttools/prttools.C"
#include "/u/aali/dirc/prttools/prttools.C"
// file existance
bool exists_test (const std::string& name);

//root plot_beam_correction.C'(20)'
void plot_beam_correction(/*TString inFile = "r_spr.root",*/ Int_t angle= 20) {




Double_t momentum=7.0;
Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
Double_t fAngleP = acos(sqrt(momentum*momentum+ mass[4]*mass[4])/momentum/1.4738)-0.00;
Double_t fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00; 
    
    prt_savepath="opt";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    prt_setRootPalette(1);
    //gStyle->SetPalette(62);
    Int_t nf = 40;
    TH2F * twoD_spr =  new TH2F("twoD_spr",";#Delta#Theta[mrad]; #Delta#Phi[mrad]", nf, -20, 20, nf, -10,10);
    TH2F * twoD_mean =  new TH2F("twoD_mean",";#Delta#Theta[mrad]; #Delta#Phi[mrad]", nf, -20, 20, nf, -10,10);
    TH2F * twoD_mean_expected =  new TH2F("twoD_mean_expected",";#Delta#Theta[mrad]; #Delta#Phi[mrad]", nf, -20, 20, nf, -10,10);
    
    //Int_t angle;
    TFile *ffile;
    TH1F *chere;
    for (int i=20; i<=150; i+=10) {
        for (int j= -200; j<=200; j+=10) {
            for (int k=-100; k<=100; k+=5) {
                if (i == angle){
                    //angle=i;
                    Double_t jj = (Float_t) j/10000.0;
                    Double_t kk = (Float_t) k/10000.0;
                    TString jj_string = Form("t1_%.3f", jj);
                    TString kk_string = Form("_t2_%.4f", kk);
                    cout<< "enter the if condition"<<endl;
                    //TString cherenkov_data_path = Form("/Users/ahmed/dirc/beam_correction/beam_correction_%d/opt_%d_3sph_"+jj_string+kk_string+"_proton_data_wo_correction_spr.root",i, i);
                    //TString cherenkov_data_path = Form("/Users/ahmed/dirc/beam_correction/beam_correction_%d/opt_%d_3sph_"+jj_string+kk_string+"_proton_data_p_correction_spr.root",i, i);
                    //TString cherenkov_data_path = Form("/Users/ahmed/dirc/beam_correction/beam_correction_%d/opt_%d_3sph_"+jj_string+kk_string+"_pi_data_pi_correction_spr.root",i, i);
                    //TString cherenkov_data_path = Form("/Users/ahmed/dirc/beam_correction/beam_correction_%d/opt_%d_3sph_"+jj_string+kk_string+"_pi_data_p_correction_spr.root",i, i);
                    TString cherenkov_data_path = Form("/data.local/beam_correction/beam_correction/beam_correction_%d/opt_%d_3sph_"+jj_string+kk_string+"_pi_data_wo_correction_spr.root",i, i);
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
                    chere->Fit("fFit","M+","",cangle-0.06,cangle+0.06);
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
                    
                    twoD_spr->Fill(jj*1000, kk*1000, spr*1000); // /14.0
                    twoD_mean->Fill(jj*1000, kk*1000, cangle-fAnglePi);
                    
                    //if( (cangle < (fAnglePi-0.0003)) || (cangle > (fAnglePi +0.0003)) ) twoD_mean_expected->Fill(jj*1000, kk*1000, cangle);
                    
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
                }
            }
        }
    }
    
    TString nid = Form("_%2.0d", angle);
    //prt_canvasAdd("r_beam_correction_p_wo",800,400);
    //prt_canvasAdd("r_beam_correction_p_wp",800,400);
    //prt_canvasAdd("r_beam_correction_pi_wpi",800,400);
    //prt_canvasAdd("r_beam_correction_pi_wp",800,400);
    //prt_canvasAdd("r_beam_correction_pi_wo",800,400);
    
    
    prt_canvasAdd("r_spr_pi_wo"+ nid,800,400);
    twoD_spr->SetStats(0);
    twoD_spr-> SetTitle(Form("SPR #pi data without mcp by mcp #theta_{c} correction %d",angle) );
    //twoD_spr-> SetTitle(Form("SPR P data without mcp by mcp #theta_{c} correction %d",angle) );
    //twoD_spr-> SetTitle("SPR P data without mcp by mcp #theta_{c} correction" );
    //twoD_spr-> SetTitle("SPR P data with mcp by mcp #theta_{c} correction" );
    //twoD_spr-> SetTitle("SPR #pi data with mcp by mcp #theta_{c} correction for #pi" );
    //twoD_spr-> SetTitle("SPR #pi data with mcp by mcp #theta_{c} correction for P" );
    //twoD_spr-> SetTitle("SPR #pi data without mcp by mcp #theta_{c} correction" );
    
    twoD_spr-> Draw("colz");
    
    prt_canvasAdd("r_mean_pi_wo"+nid,800,400);
    gStyle->SetPaintTextFormat(".3f");
    twoD_mean->SetStats(0);
    twoD_mean-> SetTitle(Form("Mean #pi data without mcp by mcp #theta_{c} correction %d",angle) );
    //twoD_mean-> SetTitle(Form("Mean P data without mcp by mcp #theta_{c} correction %d",angle) );
    //twoD_mean->SetMaximum(0.84);
    //twoD_mean->SetMinimum(0.81);
    twoD_mean-> Draw("colz");
    
   
    
    prt_canvasAdd("r_mean_pi_expected_wo"+nid,800,400);
    
   TExec *ex1 = new TExec("ex1","grey();");
   TExec *ex2 = new TExec("ex2","gStyle->SetPalette(53,0,0.4);");
   twoD_mean->Draw("col");
   ex1->Draw();
   twoD_mean->Draw("col same");
   ex2->Draw();
   twoD_mean_expected->Draw("col same");
    
    

   prt_canvasSave(2,0);
   prt_canvasDel("*");
}



void grey()
{
   Double_t Red[2]   = { 0.10, 0.10};
   Double_t Green[2] = { 0.20, 0.50};
   Double_t Blue[2]  = { 0.60, 0.90};
   Double_t Stops[2] = { 0.00, 1.00};

   Int_t nb=50;
   TColor::CreateGradientColorTable(2,Stops,Red,Green,Blue,nb);
}


//////////////////////////
// check file existance //
//////////////////////////
bool exists_test (const std::string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}
