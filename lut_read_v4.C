#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../src/PrtLutNode.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

#include "TGraph2D.h"
#include "TRotation.h"
#include "TMath.h"
#include "TROOT.h"
#include "TVector3.h"

//10m_l3_20_f_final.root
//10m_l3_180_f_0bar_final.root
//10m_l3_180_f_inbar_final.root
void lut_read_v4(TString infile="/home/a/Desktop/work/phs_reco/std/prtdirc/build_old/1m_l3_20_rot_final.root"){ //10m_l3_20_f_rot_final.root") { 10m_l3_20_f_final.root
    gStyle->SetOptStat(0);

    TH1F*  hist_phs_time = new TH1F("hist_phs_time",";time [ns];entries [#]", 250,0,50);
    TH1F*  hist_phs_time_direct = new TH1F("hist_phs_time_direct",";measured time [ns];entries [#]",   250,0,50);
    TH1F*  hist_phs_time_not_direct = new TH1F("hist_phs_time_not_direct",";measured time [ns];entries [#]",   250,0,50);


    TH1F*  hist_phs_time_m_div_x = new TH1F("hist_phs_time_m_div_x",";measured time [ns];entries [#]",   250,0,50);
    TH1F*  hist_phs_time_p_div_x = new TH1F("hist_phs_time_p_div_x",";measured time [ns];entries [#]",   250,0,50);
    TH1F*  hist_phs_time_m_div_y = new TH1F("hist_phs_time_m_div_y",";measured time [ns];entries [#]",   250,0,50);
    TH1F*  hist_phs_time_p_div_y = new TH1F("hist_phs_time_p_div_y",";measured time [ns];entries [#]",   250,0,50);
    TH1F*  hist_phs_time_m_div_z = new TH1F("hist_phs_time_m_div_z",";measured time [ns];entries [#]",   250,0,50);
    TH1F*  hist_phs_time_p_div_z = new TH1F("hist_phs_time_p_div_z",";measured time [ns];entries [#]",   250,0,50);


    TH1F*  hist_time = new TH1F("hist_time",";measured time [ns];entries [#]",   250,0,50);
    TH2F*  hist_xy = new TH2F("hist_xy",";pos y [mm];pos x [mm]", 800,-30,30, 800,-30,30);
    TH1F*  hist_z = new TH1F("hist_z",";pos z [mm];entries [#]",   1000,-1200,1200);

    // TH1F*  hist_dir_x = new TH1F("hist_dir_x",";dir x component ;entries [#]", 100,-1.0,1.0);
    // TH1F*  hist_dir_y = new TH1F("hist_dir_y",";dir y component ;entries [#]", 100,-1.0,1.0);
    // TH1F*  hist_dir_z = new TH1F("hist_dir_z",";dir z component;entries [#]", 100,-1.0,1.0);

    TH1F*  hist_pix = new TH1F("hist_pix",";oix num;entries [#]",100 ,0,100);
    TH1F*  hist_mcp = new TH1F("hist_mcp",";oix num;entries [#]",100 ,0,100);


    TH1F*  hist_angle = new TH1F("hist_angle",  "phs angle;#theta_{C} [rad];entries [#]",300,0,5 );//  ,80,0.6,1);//  ,300,0,5 );//
    TH1F*  hist_dir_z_cut = new TH1F("hist_dir_z_cut",";dir z component;entries [#]", 200,-1.5,1.5);


    TH2F*  hist_time_angle = new TH2F("hist_time_angle",";time [ns];solution angle [rad]", 250,0,50,400,0,4 );
    TH2F*  hist_time_angle_cut = new TH2F("hist_time_angle_cut",";time [ns];solution angle [rad]", 250,0,50,400,0,4 );
    TH2F*  hist_time_angle_all = new TH2F("hist_time_angle_all",";time [ns];solution angle [rad]", 250,0,50,400,0,4 );
    TH2F*  hist_time_angle_all_cut = new TH2F("hist_time_angle_all_cut",";time [ns];solution angle [rad]", 250,0,50,400,0,4 );
    //////
    TH2F*  hist_dir_x = new TH2F("hist_dir_x",";dir x component ;entries [#]", 400,-1.0,1.0, 800, 0, 4);
    TH2F*  hist_dir_y = new TH2F("hist_dir_y",";dir y component ;entries [#]", 400,-1.0,1.0, 800, 0, 4);
    TH2F*  hist_dir_z = new TH2F("hist_dir_z",";dir z component;entries [#]",  400,-1.0,1.0, 800, 0, 4);
    TH2F*  hist_dir_xy = new TH2F("hist_dir_xy",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);
    TH2F*  hist_pos_phs_x_angle = new TH2F("hist_pos_phs_x_angle",";pos x component [mm];reco angel [rad]",  800,-30,30,400, 0, 4);
    TH2F*  hist_pos_phs_y_angle = new TH2F("hist_pos_phs_y_angle",";pos y component [mm];reco angel [rad]",  800,-30,30,400, 0, 4);
    //////

    TH2F* hist_time_angle_ch[12][64], *hist_time_angle_ch2[12][64], *hist_time_angle_ch3[12][64];
    for(Int_t m=0; m<12; m++) {
        for(Int_t p=0; p<64; p++) {

            hist_time_angle_ch[m][p] = new TH2F(Form("hist_time_angle_ch_%d_%d",m,p),Form("hist_time_angle_ch_%d_%d; time [ns]; solution angle [rad]",m,p), 250,0,50,400,0,4 );
            hist_time_angle_ch2[m][p] = new TH2F(Form("hist_time_angle_ch2_%d_%d",m,p),Form("hist_time_angle_ch2_%d_%d; time [ns]; solution angle [rad]",m,p), 250,0,50,400,0,4 );
            hist_time_angle_ch3[m][p] = new TH2F(Form("hist_time_angle_ch3_%d_%d",m,p),Form("hist_time_angle_ch3_%d_%d; time [ns]; solution angle [rad]",m,p), 250,0,50,400,0,4 );
        }
    }

TH2F*  hist2_dir_xy = new TH2F("hist2_dir_xy",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);
TH2F*  hist_dir_xy_test = new TH2F("hist_dir_xy_test",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);
TH2F* hist_dir_xy_angle = new TH2F("hist_dir_xy_angle",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);

TH2F* hist_dir_xy_time = new TH2F("hist_dir_xy_time",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);
TH2F* hist_dir_xy_time_time = new TH2F("hist_dir_xy_time_time",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);



  TFile *fFileNew = TFile::Open( "test.root", "RECREATE");
  TClonesArray *fLutNew;
  TTree *fTreeNew = new TTree("prtlut","Look-up table for DIRC. Averaged");
  fLutNew = new TClonesArray("PrtLutNode");
  fTreeNew->Branch("LUT",&fLutNew,256000,2); 

  Int_t Nnodes = 5000;
  TClonesArray &fLutaNew = *fLutNew;
  for (Long64_t n=0; n<Nnodes; n++) {
    new((fLutaNew)[n]) PrtLutNode(-1);
   }
   

    Int_t kt2(-1), ka2(-1);
    Double_t  average_bin(0),average_bin_time(0), content_hist_dir_xy(0), content_hist_dir_xy_test(0), content_hist_dir_xy_time(0);


    if(!prt_init(infile,1,"data/drawHPt")) return;
    PrtHit hit;
    if(infile=="") return;
    const int nmcp = 12, npix = 64;
    TFile* f = new TFile(infile);
    TTree* t = (TTree*)f->Get("prtlut");
    TClonesArray *fLut = new TClonesArray("PrtLutNode");
    t->SetBranchAddress("LUT",&fLut);
    t->GetEntry(0);
    Double_t evtime;
    TVector3 dird;
    Int_t mcpid, pixid;
    PrtLutNode *fLutNode[5000];
    PrtLutNode *node;



    Double_t pos_x, pos_y, pos_z ;
    Double_t dir_x, dir_y, dir_z;

    Double_t rad = TMath::Pi()/180.;
    Double_t prtangle = 145;
    TVector3 momInBar(0,0,-1); // -1

    Int_t kt(-1), ka(-1), nEntries(-1);

    // in the intergration should convert the mcp and pix id to i derectly to get the the right Node

    for(int i=0; i<5000; i++) {
        node = (PrtLutNode*) fLut->At(i);
        Int_t size = node->Entries();
        if(size > 0) {
            //cout<<"node "<<i<<" has "<<size<<endl;
            mcpid = i/100;   // old
            pixid = i%100; // old


            nEntries=0;
            for(Int_t m=0; m<12; m++) {
                for(Int_t p=0; p<64; p++) {
                    hist_time_angle_ch[m][p]->Reset();
                }
            }

            for(int u=0; u<size; u++) { // size
                evtime = node->GetTime(u);
                //hPTime[mcpid][pixid]->Fill(evtime);
                prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
                //cout<<"mcpid "<<mcpid<<" pixid "<<pixid<<endl;
                //cout<<"evtime "<<evtime<<endl;
                //Int_t sensorId = 100*mcpid + pixid;
                //cout<<"##### sensorId=   "<<sensorId<<endl;


                //if (evtime > 10)continue;



                pos_x = node->GetHitPos(u).X();
                pos_y = node->GetHitPos(u).Y();
                pos_z = node->GetHitPos(u).Z();

                dird   = node->GetEntry(u).Unit();
                dir_x=dird.X();
                dir_y=dird.Y();
                dir_z=dird.Z();
//if (dir_z > 0) continue;
                //cout<<"pos_x "<<pos_x<<" pos_y  "<<pos_y<<endl;
                hist_xy->Fill(pos_y, pos_x);
                hist_z->Fill(pos_z);
                //hist_dir_x->Fill(dir_x);
                //hist_dir_y->Fill(dir_y);
                //hist_dir_z->Fill(dir_z);


                hist_pix->Fill(pixid);
                hist_mcp->Fill(mcpid);


                if (dir_x > 0) hist_phs_time_direct->Fill(evtime );
                if (dir_x < 0) hist_phs_time_not_direct->Fill(evtime );


                if (dir_x < 0)hist_phs_time_m_div_x->Fill(evtime );
                if (dir_x > 0)hist_phs_time_p_div_x->Fill(evtime );
                if (dir_y < 0)hist_phs_time_m_div_y->Fill(evtime );
                if (dir_y > 0)hist_phs_time_p_div_y->Fill(evtime );
                if (dir_z < 0)hist_phs_time_m_div_z->Fill(evtime );
                if (dir_z > 0)hist_phs_time_p_div_z->Fill(evtime );






                hist_phs_time->Fill(evtime );


                Double_t angle = momInBar.Angle(dird);
                //cout<<"##### angle=   "<<angle/rad<<endl;
                hist_angle->Fill(angle);


                if (angle < 0.85 && angle < 0.75  ) {
                    hist_dir_z_cut->Fill(dir_z);
                    // hist_dir_x->Fill(dir_x);
                    // hist_dir_y->Fill(dir_y);
                }

                hist_time_angle_all->Fill(evtime,angle);
                hist_time_angle_ch3[mcpid][pixid]->Fill(evtime,angle);

                hist_time_angle_ch[mcpid][pixid]->Fill(evtime,angle);
                kt = hist_time_angle_ch[mcpid][pixid]->GetXaxis()->FindBin(evtime);
                ka = hist_time_angle_ch[mcpid][pixid]->GetYaxis()->FindBin(angle);
                nEntries =hist_time_angle_ch[mcpid][pixid]->GetBinContent(kt,ka);
                hist_dir_xy->Fill(dir_y,dir_x);



                    if (dir_x > 0 ) continue; // for 20 deg
                    //if (time_phs<ref_point) continue;
                    if (dir_z > 0 ) continue;


// if (!(mcpid==5 && pixid==5)) continue;

                    hist2_dir_xy->Fill(dir_y, dir_x,angle);
                    hist_dir_xy_test->Fill(dir_y, dir_x);
                    hist_dir_xy_time->Fill(dir_y, dir_x,evtime);// time_phs

                    kt2 = hist2_dir_xy->GetXaxis()->FindBin(dir_x);
                    ka2 = hist2_dir_xy->GetYaxis()->FindBin(dir_y);
                    
                    content_hist_dir_xy=hist2_dir_xy->GetBinContent(ka2,kt2);
                    
                    if (content_hist_dir_xy< 15) continue;

                    content_hist_dir_xy_test=hist_dir_xy_test->GetBinContent(ka2,kt2);
                    content_hist_dir_xy_time=hist_dir_xy_time->GetBinContent(ka2,kt2);
                     
                    average_bin=content_hist_dir_xy/content_hist_dir_xy_test;
                    average_bin_time=content_hist_dir_xy_time/content_hist_dir_xy_test;

                    hist_dir_xy_angle->SetBinContent(ka2,kt2,average_bin);
                    
                    hist_dir_xy_time_time->SetBinContent(ka2,kt2,average_bin_time);



   ((PrtLutNode*)( fLutNew->At(i)))->AddEntry(node->GetDetectorId(), node->GetEntry(u).Unit(), 0,0,node->GetTime(u),node->GetHitPos(u), node->GetHitPosGlobal(u), node->GetDigiPos());
   
    //fTreeNew->Fill();


                if (nEntries>1) continue;
                //if (mcpid==5 && pixid==5) cout<<"##### mcpid=   "<<mcpid<<"  pixid= "<< pixid <<"  nEntries= "<<nEntries<<endl;
                //if (mcpid==6 && pixid==6) cout<<"##### mcpid=   "<<mcpid<<"  pixid= "<< pixid <<"  nEntries= "<<nEntries<<endl;
                hist_time_angle_ch2[mcpid][pixid]->Fill(evtime,angle);
                hist_time_angle_all_cut->Fill(evtime,angle);
                //if (pixid == 20) prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
                //if (mcpid==8 &&( pixid==27||pixid==19 )) prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
                //prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);

                //hist_dir_xy->Fill(dir_y,dir_x);

                hist_dir_x->Fill(dir_x,angle);
                hist_dir_y->Fill(dir_y,angle);
                hist_dir_z->Fill(dir_z,angle);
                hist_pos_phs_x_angle->Fill(pos_x,angle);
                hist_pos_phs_y_angle->Fill(pos_y,angle);


            }
        }
    }
    
    
    
        prt_canvasAdd("r_dirxy_angle",800,400);
    hist_dir_xy_angle->Draw("colz");
    
    
           prt_canvasAdd("r_dirxy",800,400);
    hist_dir_xy_test->Draw("colz");
    
    
    prt_canvasAdd("r_dirxy_time",800,400);
    hist_dir_xy_time_time->Draw("colz");
    
    
    
    
    /*
    prt_canvasAdd("r_time_angle_all",800,400);
    hist_time_angle_all->Draw("colz");


    prt_canvasAdd("r_time_angle_all_cut",800,400);
    hist_time_angle_all_cut->Draw("colz");

*/

/*
    {
        int mm=8;
        int pp=27;
        prt_canvasAdd("r_time_angle_8_27",800,400);
        hist_time_angle_ch3[mm][pp]->Draw("colz");
    }

    {
        int mm=8;
        int pp=19;
        prt_canvasAdd("r_time_angle_8_19_cut",800,400);
        hist_time_angle_ch3[mm][pp]->Draw("colz");
    }
*/

/*
    prt_canvasAdd("r_dir_xy",800,400);
    hist_dir_xy->Draw("colz");
    prt_canvasAdd("r_dir_x",800,400);
    hist_dir_x->Draw("colz");
    prt_canvasAdd("r_dir_y",800,400);
    hist_dir_y->Draw("colz");
    prt_canvasAdd("r_dir_z",800,400);
    hist_dir_z->Draw("colz");
    prt_canvasAdd("r_pos_x",800,400);
    hist_pos_phs_x_angle->Draw("colz");
    prt_canvasAdd("r_pos_y",800,400);
    hist_pos_phs_y_angle->Draw("colz");


*/
/*
    prt_canvasAdd("r_time",800,400);
    hist_phs_time->Draw();
    prt_canvasAdd("r_time_direct",800,400);
    hist_phs_time_direct->Draw();
    prt_canvasAdd("r_time_not_Direct",800,400);
    hist_phs_time_not_direct->Draw();
*/


/*
    prt_canvasAdd("r_hist_phs_time_m_div_x",800,400);
    hist_phs_time_m_div_x->Draw();
    prt_canvasAdd("r_hist_phs_time_p_div_x",800,400);
    hist_phs_time_p_div_x->Draw();
    prt_canvasAdd("r_hist_phs_time_m_div_y",800,400);
    hist_phs_time_m_div_y->Draw();
    prt_canvasAdd("r_hist_phs_time_p_div_y",800,400);
    hist_phs_time_p_div_y->Draw();
    prt_canvasAdd("r_hist_phs_time_m_div_z",800,400);
    hist_phs_time_m_div_z->Draw();
    prt_canvasAdd("r_hist_phs_time_p_div_z",800,400);
    hist_phs_time_p_div_z->Draw();

 */


    /*

        prt_canvasAdd("r_hist_pos_xy",800,400);
        hist_xy->Draw("colz");
        prt_canvasAdd("r_hist_pos_z",800,400);
        hist_z->Draw();
        //    prt_canvasAdd("r_graph_pos_xy",800,400);
        //    graph_pos->Draw();
        //
        //    prt_canvasAdd("r_graph_dir_xy",800,400);
        //    graph_dir->Draw();

        prt_canvasAdd("r_dir_x",800,400);
        hist_dir_x->Draw();
        prt_canvasAdd("r_dir_y",800,400);
        hist_dir_y->Draw();
        prt_canvasAdd("r_dir_z",800,400);
        hist_dir_z->Draw();
        prt_canvasAdd("r_dir_z_cut",800,400);
        hist_dir_z_cut->Draw();


        prt_canvasAdd("r_phs_angle",800,400);
        hist_angle->Draw();



        prt_canvasAdd("r_pix_num",800,400);
        hist_pix->Draw();


        prt_canvasAdd("r_mcp_num",800,400);
        hist_mcp->Draw();
    */


    prt_canvasSave(2,0);
    prt_canvasDel("*");

    //
    //    for(Int_t i=0; i<770; i++) {
    //        Int_t direction_lut =fLutNode[i]->Entries();
    //        hist_ambiguity->Fill(i, direction_lut);
    //    }
    //    Double_t pos_x,pos_y, pos_z;
    //    for (Int_t mcpid_int=0; mcpid_int<12; mcpid_int++){
    //        for (Int_t pixid_int=1; pixid_int<65; pixid_int++){
    //            Int_t sensorId_int = 100*mcpid_int + pixid_int ;
    //            pos_x   = fLutNode[sensorId_int]->GetDigiPos().X();
    //            pos_y   = fLutNode[sensorId_int]->GetDigiPos().Y();
    //            pos_z   = fLutNode[sensorId_int]->GetDigiPos().Z();
    //            lut_pix_pos_xy->Fill(pos_x ,pos_y);
    //        }
    //    }



    prt_drawDigi("m,p,v\n",2017,0,0);
    prt_cdigi->SetName(Form("hp_dataProtonS332_%d_%2.1f",(Int_t)prt_theta,prt_phi));
    prt_canvasAdd(prt_cdigi);
    prt_cdigi_palette->Draw();
    prt_canvasSave(1,0);



  fTreeNew->Fill();
  //fTreeNew->Write();
  
  fTreeNew->Print();
  fTreeNew->AutoSave();
  
  

}








