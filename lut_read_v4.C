#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../src/PrtLutNode.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

//#include "TGraph2D.h"
//#include "TRotation.h"
#include "TMath.h"
#include "TROOT.h"
#include "TVector3.h"
#include "TEllipse.h"


void lut_read_v4(TString infile="/Users/ahmed/Desktop/std/prtdirc/build/10m_l3_20_f_rot_final.root"){ //1m_l3_20_rot_final.root // test_lut_pathid.root
    gStyle->SetOptStat(0);
    
    
    double R =0.025; // circle redues 0.05, 0.025
    Int_t vsize= 500; // number of cicles per pix 500, 500
    Int_t threshold= 2;// solution threshold per circle 30  , 200
    
    
    // technical
    
    int ncircle(50);
    int counter(0);
    TVector3 dummy(1,1,1);
    
    ///////////////////////
    //// Debuging Histo ///
    ///////////////////////
    
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
    TH2F*  hist_dir_x_angle = new TH2F("hist_dir_x_angle",";dir x component ;solution angle [rad]", 400,-1.0,1.0, 800, 0, 4);
    TH2F*  hist_dir_y_angle = new TH2F("hist_dir_y_angle",";dir y component ;solution angle [rad]", 400,-1.0,1.0, 800, 0, 4);
    TH2F*  hist_dir_z_angle = new TH2F("hist_dir_z_angle",";dir z component;solution angle [rad]",  400,-1.0,1.0, 800, 0, 4);
    //////
    TH2F*  hist_pos_phs_x_angle = new TH2F("hist_pos_phs_x_angle",";pos x component [mm];reco angel [rad]",  800,-30,30,400, 0, 4);
    TH2F*  hist_pos_phs_y_angle = new TH2F("hist_pos_phs_y_angle",";pos y component [mm];reco angel [rad]",  800,-30,30,400, 0, 4);
    
    TH2F*  hist_dir_xy = new TH2F("hist_dir_xy",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);
    //////
    
    
    
    //////////////////////
    //// 2D histograms ///
    //////////////////////
    
    TH2F* hist_time_angle_ch[12][64], *hist_time_angle_ch2[12][64], *hist_time_angle_ch3[12][64];
    TH2F*  hist2_dir_xy[12][64], *hist_dir_xy_test[12][64], *hist_dir_xy_occu[12][64],*hist_dir_xy_occu_final[12][64],*hist_dir_xy_occu_final2[12][64], *hist_dir_xy_angle[12][64], *hist_dir_xy_time[12][64], *hist_dir_xy_time_time[12][64];
    
    int bin_histo= 400; //400
    for(Int_t m=0; m<12; m++) {
        for(Int_t p=0; p<64; p++) {
            
            hist_time_angle_ch[m][p] = new TH2F(Form("hist_time_angle_ch_%d_%d",m,p),Form("hist_time_angle_ch_%d_%d; time [ns]; solution angle [rad]",m,p), 250,0,50,400,0,4 );
            hist_time_angle_ch2[m][p] = new TH2F(Form("hist_time_angle_ch2_%d_%d",m,p),Form("hist_time_angle_ch2_%d_%d; time [ns]; solution angle [rad]",m,p), 250,0,50,400,0,4 );
            hist_time_angle_ch3[m][p] = new TH2F(Form("hist_time_angle_ch3_%d_%d",m,p),Form("hist_time_angle_ch3_%d_%d; time [ns]; solution angle [rad]",m,p), 250,0,50,400,0,4 );
            
            
            hist2_dir_xy[m][p] = new TH2F(Form("hist2_dir_xy_ch_%d_%d",m,p),Form("hist2_dir_xy_ch_%d_%d;dir y component;dir x component",m,p),  bin_histo,-1.0,1.0, bin_histo,-1.0,1.0);
            hist_dir_xy_test[m][p] = new TH2F(Form("hist_dir_xy_test_ch_%d_%d",m,p),Form("hist_dir_xy_test_ch_%d_%d;dir y component;dir x component",m,p),  bin_histo,-1.0,1.0, bin_histo,-1.0,1.0);
            hist_dir_xy_occu[m][p] = new TH2F(Form("hist_dir_xy_occu_ch_%d_%d",m,p),Form("hist_dir_xy_occu_ch_%d_%d;dir y component;dir x component",m,p),  bin_histo,-1.0,1.0, bin_histo,-1.0,1.0);
            hist_dir_xy_occu_final[m][p] = new TH2F(Form("hist_dir_xy_occu_final_ch_%d_%d",m,p),Form("hist_dir_xy_occu_final_ch_%d_%d;dir y component;dir x component",m,p),  bin_histo,-1.0,1.0, bin_histo/2,-1.0,0);
            hist_dir_xy_occu_final2[m][p] = new TH2F(Form("hist_dir_xy_occu_final2_ch_%d_%d",m,p),Form("hist_dir_xy_occu_final2_ch_%d_%d;dir y component;dir x component",m,p),  bin_histo,-1.0,1.0, bin_histo/2,-1.0,0);
            
            hist_dir_xy_angle[m][p] = new TH2F(Form("hist_dir_xy_angle_ch_%d_%d",m,p),Form("hist_dir_xy_angle_ch_%d_%d;dir y component;dir x component",m,p),  bin_histo,-1.0,1.0, bin_histo,-1.0,1.0);
            
            hist_dir_xy_time[m][p] = new TH2F(Form("hist_dir_xy_time_ch_%d_%d",m,p),Form("hist_dir_xy_time_ch_%d_%d;dir y component;dir x component",m,p),  bin_histo,-1.0,1.0, bin_histo,-1.0,1.0);
            hist_dir_xy_time_time[m][p] = new TH2F(Form("hist_dir_xy_time_time_ch_%d_%d",m,p),Form("hist_dir_xy_time_time_ch_%d_%d;dir y component;dir x component",m,p),  bin_histo,-1.0,1.0, bin_histo,-1.0,1.0);
        }
    }
    
    TH2F*  hist_dir_xy_occu_circle[ncircle];
    int bin_histo2= 50;
    for(Int_t c=0; c<ncircle; c++) {
        hist_dir_xy_occu_circle[c] = new TH2F(Form("hist_dir_xy_occu_circle_ch_%d",c),Form("hist_dir_xy_occu_circle_ch_%d;dir y component;dir x component",c),  bin_histo2,-1.0,1.0, bin_histo2,-1.0,1.0);
    }
    
    
    //////////////////////
    ////create New LUT ///
    //////////////////////
    
    TFile *fFileNew = TFile::Open( "test_lut.root", "RECREATE");
    TClonesArray *fLutNew;
    TTree *fTreeNew = new TTree("prtlut","Look-up table for DIRC. Averaged");
    fLutNew = new TClonesArray("PrtLutNode");
    fTreeNew->Branch("LUT",&fLutNew,256000,2);
    
    Int_t Nnodes = 5000;
    TClonesArray &fLutaNew = *fLutNew;
    for (Long64_t n=0; n<Nnodes; n++) {
        new((fLutaNew)[n]) PrtLutNode(-1);
    }
    
    /////////////////
    ////Variables ///
    /////////////////
    
    Int_t kt2(-1), ka2(-1);
    Double_t  average_bin(0),average_bin_time(0), content_hist_dir_xy(0), content_hist_dir_xy_test(0), content_hist_dir_xy_time(0);
    Double_t pathid(-1);
    
    /////////////////////////
    ////clustrung vectors ///
    /////////////////////////
    //Int_t vsize= 500;
    std::vector<TVector3> dirArrayVector[vsize];
    std::vector<TVector3> posGlobalArrayVector[vsize];
    std::vector<TVector3> posLocalArrayVector[vsize];
    
    std::vector<Double_t> tArrayVector[vsize];
    std::vector<Double_t> pathIDVector;
    TVector3  sum;
    Double_t sumt;
    
    std::vector<TVector3> circleIDVector[12][64];
    
    ////////////////////////////////
    //// circiles inisialization ///
    ////////////////////////////////
    
    //double R =0.07;
    TEllipse *el[12][64][vsize];
    for(Int_t m=0; m<12; m++) {
        for(Int_t p=0; p<64; p++) {
            for(Int_t e=0; e<vsize; e++) {
                
                el[m][p][e] = new TEllipse(-2,-2,R,R);
                
            }
        }
    }
    
    
    
    if(!prt_init(infile,1,"data/drawHPt")) return;
    PrtHit hit;
    if(infile=="") return;
    const int nmcp = 12, npix = 64;
    TFile* f = new TFile(infile);
    TTree* t = (TTree*)f->Get("prtlut");
    TClonesArray *fLut = new TClonesArray("PrtLutNode");
    t->SetBranchAddress("LUT",&fLut);
    t->GetEntry(0);
    Double_t time_solution(-1),  DetectorId(-1);
    TVector3 dird;
    
    Int_t mcpid, pixid;
    PrtLutNode *fLutNode[5000];
    PrtLutNode *node;
    
    
    
    Double_t pos_x, pos_y, pos_z ;
    Double_t dir_x, dir_y, dir_z;
    
    Double_t rad = TMath::Pi()/180.;
    Double_t prtangle = 145;
    TVector3 momInBar(0,0,-1);
    
    Int_t kt(-1), ka(-1), nEntries(-1);
    
    Int_t Xbin_pattern(-1), Ybin_pattern(-1);
    Double_t content_prt_hdigi(0);
    

    ///////////////////////
    /// Loop over nodes ///
    ///////////////////////
    
    for(int i=0; i<2000; i++) {
        node = (PrtLutNode*) fLut->At(i);
        Int_t size = node->Entries();
        if(size<1) continue;
        mcpid = i/100;
        pixid = i%100;

        
        nEntries=0;
        for(Int_t m=0; m<12; m++) {
            for(Int_t p=0; p<64; p++) {
                hist_time_angle_ch[m][p]->Reset();
            }
        }
        
        
        ///////////////////////////////////
        /// Loop over soluions per node ///
        ///////////////////////////////////
        
        for(int u=0; u<size; u++) { // size
            time_solution = node->GetTime(u);
            
            Xbin_pattern = prt_hdigi[mcpid]->GetXaxis()->FindBin(pixid%8);
            Ybin_pattern = prt_hdigi[mcpid]->GetYaxis()->FindBin(pixid/8);
            content_prt_hdigi=prt_hdigi[mcpid]->GetBinContent(Xbin_pattern,Ybin_pattern);
            if (content_prt_hdigi> 500) continue;

            
            prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
            
            ///////////////////////////
            /// Extract information ///
            ///////////////////////////
            
            //pathid = node->GetPathId(u);
            //hPTime[mcpid][pixid]->Fill(time_solution);
            //Int_t sensorId = 100*mcpid + pixid;
            DetectorId = node->GetDetectorId();
            pos_x = node->GetHitPos(u).X();
            pos_y = node->GetHitPos(u).Y();
            pos_z = node->GetHitPos(u).Z();
            
            dird   = node->GetEntry(u).Unit();
            dir_x=dird.X();
            dir_y=dird.Y();
            dir_z=dird.Z();

            
            hist_xy->Fill(pos_y, pos_x);
            hist_z->Fill(pos_z);
            //hist_dir_x->Fill(dir_x);
            //hist_dir_y->Fill(dir_y);
            //hist_dir_z->Fill(dir_z);
            
            
            hist_pix->Fill(pixid);
            hist_mcp->Fill(mcpid);
            
            
            if (dir_x > 0) hist_phs_time_direct->Fill(time_solution );
            if (dir_x < 0) hist_phs_time_not_direct->Fill(time_solution );
            
            
            if (dir_x < 0)hist_phs_time_m_div_x->Fill(time_solution );
            if (dir_x > 0)hist_phs_time_p_div_x->Fill(time_solution );
            if (dir_y < 0)hist_phs_time_m_div_y->Fill(time_solution );
            if (dir_y > 0)hist_phs_time_p_div_y->Fill(time_solution );
            if (dir_z < 0)hist_phs_time_m_div_z->Fill(time_solution );
            if (dir_z > 0)hist_phs_time_p_div_z->Fill(time_solution );
            
            
            hist_phs_time->Fill(time_solution );
            
            
            Double_t angle = momInBar.Angle(dird);
            //cout<<"##### angle=   "<<angle/rad<<endl;
            hist_angle->Fill(angle);
            
            
            if (angle < 0.85 && angle < 0.75  ) {
                hist_dir_z_cut->Fill(dir_z);
                // hist_dir_x->Fill(dir_x);
                // hist_dir_y->Fill(dir_y);
            }
            
            hist_time_angle_all->Fill(time_solution,angle);
            hist_time_angle_ch3[mcpid][pixid]->Fill(time_solution,angle);
            
            hist_time_angle_ch[mcpid][pixid]->Fill(time_solution,angle);
            kt = hist_time_angle_ch[mcpid][pixid]->GetXaxis()->FindBin(time_solution);
            ka = hist_time_angle_ch[mcpid][pixid]->GetYaxis()->FindBin(angle);
            nEntries =hist_time_angle_ch[mcpid][pixid]->GetBinContent(kt,ka);
            hist_dir_xy->Fill(dir_y,dir_x);
            
            
            
            if (dir_x > 0 ) continue; // for 20 deg
            if (dir_z > 0 ) continue;
            
            
            //if (!(mcpid==5 && pixid==5)) continue;
            
            hist2_dir_xy[mcpid][pixid]->Fill(dir_y, dir_x,angle);
            hist_dir_xy_test[mcpid][pixid]->Fill(dir_y, dir_x);
            hist_dir_xy_time[mcpid][pixid]->Fill(dir_y, dir_x,time_solution);// time_phs
            
            kt2 = hist2_dir_xy[mcpid][pixid]->GetXaxis()->FindBin(dir_x);
            ka2 = hist2_dir_xy[mcpid][pixid]->GetYaxis()->FindBin(dir_y);
            
            content_hist_dir_xy=hist2_dir_xy[mcpid][pixid]->GetBinContent(ka2,kt2);
            
            // if (content_hist_dir_xy< 1) continue;
            
            content_hist_dir_xy_test=hist_dir_xy_test[mcpid][pixid]->GetBinContent(ka2,kt2);
            content_hist_dir_xy_time=hist_dir_xy_time[mcpid][pixid]->GetBinContent(ka2,kt2);
            
            average_bin=content_hist_dir_xy/content_hist_dir_xy_test;
            average_bin_time=content_hist_dir_xy_time/content_hist_dir_xy_test;
            
            hist_dir_xy_angle[mcpid][pixid]->SetBinContent(ka2,kt2,average_bin);
            hist_dir_xy_time_time[mcpid][pixid]->SetBinContent(ka2,kt2,average_bin_time);
            hist_dir_xy_occu[mcpid][pixid]->Fill(dir_y, dir_x);
            
            TAxis *xaxis = hist_dir_xy_occu[mcpid][pixid]->GetXaxis();
            TAxis *yaxis = hist_dir_xy_occu[mcpid][pixid]->GetYaxis();
            
            Double_t binCenter_x = xaxis->GetBinCenter(ka2);
            Double_t binCenter_y = xaxis->GetBinCenter(kt2);
            
            //dird.SetX(binCenter_y);
            //dird.SetY(binCenter_x);
            
            //hist_dir_x_angle->Fill(dird.X(),angle);
            //hist_dir_y_angle->Fill(dird.Y(),angle);
            
            //((PrtLutNode*)( fLutNew->At(i)))->AddEntry(node->GetDetectorId(), node->GetEntry(u).Unit(), pathid,0,node->GetTime(u),node->GetHitPos(u), node->GetHitPosGlobal(u), node->GetDigiPos());
            // not ((PrtLutNode*)( fLutNew->At(i)))->AddEntry(node->GetDetectorId(), dird, 0,0,average_bin_time,node->GetHitPos(u), node->GetHitPosGlobal(u), node->GetDigiPos());
            
            hist_dir_xy_occu_final[mcpid][pixid]->Fill(dird.Y(), dird.X());
            
            
            /*
             if (nEntries>1) continue;
             //if (mcpid==5 && pixid==5) cout<<"##### mcpid=   "<<mcpid<<"  pixid= "<< pixid <<"  nEntries= "<<nEntries<<endl;
             //if (mcpid==6 && pixid==6) cout<<"##### mcpid=   "<<mcpid<<"  pixid= "<< pixid <<"  nEntries= "<<nEntries<<endl;
             hist_time_angle_ch2[mcpid][pixid]->Fill(time_solution,angle);
             hist_time_angle_all_cut->Fill(time_solution,angle);
             //if (pixid == 20) prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
             //if (mcpid==8 &&( pixid==27||pixid==19 )) prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
             //prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
             
             //hist_dir_xy->Fill(dir_y,dir_x);
             
             //hist_dir_x_angle->Fill(dir_x,angle);
             //hist_dir_y_angle->Fill(dir_y,angle);
             hist_dir_z_angle->Fill(dir_z,angle);
             hist_pos_phs_x_angle->Fill(pos_x,angle);
             hist_pos_phs_y_angle->Fill(pos_y,angle);
             
             */
            
            
            
            
            /////////////////
            ////clustring ///
            /////////////////
            
            bool newid = true;
            
            // loop over existed circle
            // calcluate distance between new (dir_x, dir_y) and circles center
            // if the destance smaller than reduis of one circle, flag it as false newid, register the point with taken pathid and exit check loop
            

            pathid= -1;
            
            map<int, Double_t> circle_map ;
            for(int e=0; e<circleIDVector[mcpid][pixid].size(); e++) {
                Double_t p1 = circleIDVector[mcpid][pixid][e].X()-dird.X() ;
                Double_t p2 = circleIDVector[mcpid][pixid][e].Y()-dird.Y() ;
                p1 *=p1;
                p2 *=p2;
                Double_t result=p1+p2;
                Double_t distance = sqrt(result);
                
                circle_map[e] = distance;
            }


            Double_t min_distance =100;
            Double_t order = -1;
            for(auto it = circle_map.cbegin(); it != circle_map.cend(); ++it ) {
                if (it ->second < min_distance) {
                    order = it->first;
                    min_distance = it->second;
                }
            }

            
            if( min_distance < R ){
                newid= false;
                pathid=mcpid*10000000+pixid*10000+order+1;
                
                for(int j=0; j<pathIDVector.size(); j++){
                    if(pathid == pathIDVector[j]){
                        dirArrayVector[j].push_back(dird);
                        tArrayVector[j].push_back(time_solution);
                        break;
                    }
                }
            }
            
            
            // if the point new, creat a new circle for it
            // add the circle to a list
            // register the point with new pathid
            if(newid) {
                
                el[mcpid][pixid][circleIDVector[mcpid][pixid].size()] = new TEllipse(dird.Y(),dird.X(),R,R);
                circleIDVector[mcpid][pixid].push_back(dird); //add the circle to a list
                pathid=mcpid*10000000+pixid*10000+circleIDVector[mcpid][pixid].size();
                if (pathid == 50050001)cout<<"################# 2nd The distance= "<<"j= "<<pathIDVector.size()<<"  circle number= "<< circleIDVector[mcpid][pixid].size()<<"  id= "<<fixed <<pathid<<" "<<dird.X() <<" "<<dird.Y()<<endl;
                //cout<<"@@@@@@@@@@@  2nd pathid = "<<pathid<<endl;
                dirArrayVector[pathIDVector.size()].push_back(dird);
                tArrayVector[pathIDVector.size()].push_back(time_solution);
                pathIDVector.push_back(pathid);
                //cout<<"#########################   New Solution  "<<pathid<<"  "<<pathid - pathIDVector[pathIDVector.size()-1]<<endl;
                //cout<<"@@@@@@@@@@@@  new pathIDVector[j] = "<<pathIDVector[pathIDVector.size()-1]<< " stor dir and t indix J = "<<pathIDVector.size()<<endl;
            }
            
            // test
            //((PrtLutNode*)( fLutNew->At(i)))->AddEntry(node->GetDetectorId(), node->GetEntry(u).Unit(), pathid,0,node->GetTime(u),dummy, dummy, dummy);
            //((PrtLutNode*)( fLutNew->At(i)))->AddEntry(node->GetDetectorId(), node->GetEntry(u).Unit(), pathid,0,node->GetTime(u),node->GetHitPos(u), node->GetHitPosGlobal(u), node->GetDigiPos());
        }
        
        /////////////////
        /////Fill LUT ///
        /////////////////
        
        
        for(int j=0; j<pathIDVector.size(); j++){
            sum = TVector3(0,0,0);
            sumt=0;
            
            //cout << "number of solutions   "<< pathIDVector.size() <<endl;
            //cout << "number of solutions have the same id  "<<dirArrayVector[j].size() <<endl;
            
            if (dirArrayVector[j].size() < threshold)continue;
            for(int v=0; v<(dirArrayVector[j].size()); v++) {
                sum += dirArrayVector[j][v];
                sumt += tArrayVector[j][v];
                ((PrtLutNode*)(fLutNew->At(i)))->AddEntry(node->GetDetectorId(), dirArrayVector[j][v],pathIDVector[j],0,tArrayVector[j][v],  dummy,dummy,dummy,dirArrayVector[j].size()/(Double_t)size);
                //hist_dir_xy_occu_final2[mcpid][pixid]->Fill(dirArrayVector[j][v].Y(), dirArrayVector[j][v].X());
                //if (counter< ncircle)hist_dir_xy_occu_circle[counter]->Fill(dirArrayVector[j][v].Y(), dirArrayVector[j][v].X());
                
                
            }
            
            if(dirArrayVector[j].size()<1) continue;
            Double_t weight = 1/(Double_t)dirArrayVector[j].size();
            sum *= weight;
            sumt *= weight;
            
            //            cout << "sum  "<<sum.X()<<endl;
            //            cout << "sumt  "<<sumt<<endl;
            //            cout << "d id  "<<node->GetDetectorId()<<endl;
            //            cout << "id  "<<pathIDVector[j]<<endl;
            //            cout << "w  "<<dirArrayVector[j].size()/(Double_t)size<<endl;
            
            
            //cout.precision(17);
            //cout << "Path ID " << fixed << pathIDVector[j] << endl;
            //cout << "dir y, x =" << fixed << sum.Y() << "   "<<sum.X()<<"  weight = "<<weight <<endl;
            
            //((PrtLutNode*)(fLutNew->At(i)))->AddEntry(node->GetDetectorId(), sum,pathIDVector[j],0,sumt,  node->GetDigiPos(),node->GetDigiPos(),node->GetDigiPos(),dirArrayVector[j].size()/(Double_t)size);
            //hist_dir_xy_occu_final2[mcpid][pixid]->Fill(sum.Y(), sum.X());
            ++counter;
        }
        for(int i=0; i<vsize; i++) {dirArrayVector[i].clear();  tArrayVector[i].clear();}
        pathIDVector.clear();
        
        for(int m=0; m<12; m++){
            for(int p=0; p<64; p++){
                
                circleIDVector[m][p].clear();
            }
        }
    }
    
    ///////////
    ///Draw ///
    ///////////
    
    //    {
    //        int mcpid=5;
    //        int pixid=5;
    //
    //        prt_canvasAdd("r_dirxy_angle",800,400);
    //        hist_dir_xy_angle[mcpid][pixid]->Draw("colz");
    //
    //
    //        prt_canvasAdd("r_dirxy",800,400);
    //        hist_dir_xy_test[mcpid][pixid]->Draw("colz");
    //
    //
    //        prt_canvasAdd("r_dirxy_time",800,400);
    //        hist_dir_xy_time_time[mcpid][pixid]->Draw("colz");
    //
    //        prt_canvasAdd("r_dirxy_occ",800,400);
    //        hist_dir_xy_occu_final[mcpid][pixid]->Draw("colz");
    //
    //
    //    }
    
    if(false){
        prt_canvasAdd("r_pix_phs1",800,400);
        for(Int_t m=0; m<12; m++) {
            for(Int_t p=0; p<64; p++) {
                hist_dir_xy_occu_final[m][p]->Draw("colz");
                for(Int_t e=0; e<vsize; e++) {

                    el[m][p][e]->SetFillColor(0);
                    el[m][p][e]->SetFillStyle(0);
                    el[m][p][e]->Draw("same");
                }
                prt_waitPrimitive("r_pix_phs1");
            }
        }
    }
    
    if(false){
        prt_canvasAdd("r_pix_phs2",800,400);
        for(Int_t m=0; m<12; m++) {
            for(Int_t p=0; p<64; p++) {
                hist_dir_xy_occu_final2[m][p]->Draw("colz");
                for(Int_t e=0; e<vsize; e++) {
                    el[m][p][e]->SetFillColor(0);
                    el[m][p][e]->SetFillStyle(0);
                    el[m][p][e]->Draw("same");
                }
                prt_waitPrimitive("r_pix_phs2");
            }
        }
    }
    
    
    if(false){
        int mcpid=5;
        int pixid=5;
        
        prt_canvasAdd("r_f1",800,400);
        hist_dir_xy_occu_final[mcpid][pixid]->Draw("colz");
        
        for(Int_t e=0; e<vsize; e++) {
            el[mcpid][pixid][e]->SetFillColor(0);
            el[mcpid][pixid][e]->SetFillStyle(0);
            el[mcpid][pixid][e]->Draw("same");
        }
        
        prt_canvasAdd("r_f2",800,400);
        hist_dir_xy_occu_final2[mcpid][pixid]->Draw("colz");
        for(Int_t e=0; e<vsize; e++) {
            el[mcpid][pixid][e]->SetFillColor(0);
            el[mcpid][pixid][e]->SetFillStyle(0);
            el[mcpid][pixid][e]->Draw("same");
        }
        
        
        prt_canvasAdd("r_circle",800,400);
        for(Int_t c=0; c<ncircle; c++) {
            hist_dir_xy_occu_circle[c]->Draw("colz");
            
            el[mcpid][pixid][c]->SetFillColor(0);
            el[mcpid][pixid][c]->SetFillStyle(0);
            el[mcpid][pixid][c]->Draw("same");
            
            prt_waitPrimitive("r_circle");
        }
        
        
        
        
        
        
    }
    
    
    
    
    
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
     hist_dir_x_angle->Draw("colz");
     prt_canvasAdd("r_dir_y",800,400);
     hist_dir_y_angle->Draw("colz");
     prt_canvasAdd("r_dir_z",800,400);
     hist_dir_z_angle->Draw("colz");
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
    fTreeNew->Print();
    fTreeNew->AutoSave();
    
    
    
}








