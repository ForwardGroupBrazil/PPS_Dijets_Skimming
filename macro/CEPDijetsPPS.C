#include <stdio.h>      /* printf */
#include <math.h>       /* sqrt */
//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMatrixT.h>
#include <TVectorT.h>
//
////STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
//#include <vector>
#include <map>
#include <cmath>
#define PI 3.141592653589793
using namespace std;
void CEPDijetsPPS()
{
    //------------ FWLite libraries ------
    gSystem->Load("libFWCoreFWLite.so");
    AutoLibraryLoader::enable();
    gROOT->ProcessLine("#include<vector>");
    //------------ files -----------------
    TFile *inf  = TFile::Open("ttreeCEPdijets.root");
    TFile *outf = new TFile("CEPPPS_Histos.root","RECREATE");
    TTree *tr = (TTree*)inf->Get("demo/Event");
    //----------- define the tree branch --------
    std::vector<double> *JetsPt =0;
    std::vector<double>  *JetsEta =0;
    std::vector<double>  *JetsPhi =0;
    Int_t           nVertex;
    Int_t           nTracks;
    Double_t        MinDistance;
    Double_t        MaxDistance;
    Double_t        GoldenVertexZ;
    Double_t        xiPPSArmB;
    Double_t        xiPPSArmF;
    Double_t        tPPSArmB;
    Double_t        tPPSArmF;
    Double_t        xPPSArmBDet1;
    Double_t        xPPSArmFDet1;
    Double_t        yPPSArmBDet1;
    Double_t        yPPSArmFDet1;
    Double_t        xPPSArmBDet2;
    Double_t        xPPSArmFDet2;
    Double_t        yPPSArmBDet2;
    Double_t        yPPSArmFDet2;
    Int_t           stopPPSArmFTrkDet1;
    Int_t           stopPPSArmFTrkDet2;
    Int_t           stopPPSArmBTrkDet1;
    Int_t           stopPPSArmBTrkDet2;
    Int_t           stopPPSArmFToF;
    Int_t           stopPPSArmBToF; 
    Double_t        Mjj;
    Double_t        Mpf;
    Double_t        Rjj;

    TBranch        *bJetsPt = 0;
    TBranch        *bJetsEta = 0;   
    TBranch        *bJetsPhi = 0;
    TBranch        *bnVertex;   //!
    TBranch        *bnTracks;   //!
    TBranch        *bMinDistance;   //!
    TBranch        *bMaxDistance;   //!
    TBranch        *bGoldenVertexZ;   //!
    TBranch        *bxiPPSArmB;   //!
    TBranch        *bxiPPSArmF;   //!
    TBranch        *btPPSArmB;   //!
    TBranch        *btPPSArmF;   //!
    TBranch        *bxPPSArmBDet1;   //!
    TBranch        *bxPPSArmFDet1;   //!
    TBranch        *byPPSArmBDet1;   //!
    TBranch        *byPPSArmFDet1;   //!
    TBranch        *bxPPSArmBDet2;   //!
    TBranch        *bxPPSArmFDet2;   //!
    TBranch        *byPPSArmBDet2;   //!
    TBranch        *byPPSArmFDet2;   //!
    TBranch        *bstopPPSArmFTrkDet1;   //!
    TBranch        *bstopPPSArmFTrkDet2;   //!
    TBranch        *bstopPPSArmBTrkDet1;   //!
    TBranch        *bstopPPSArmBTrkDet2;   //!
    TBranch        *bstopPPSArmFToF;   //!
    TBranch        *bstopPPSArmBToF;   //!
    TBranch        *bMjj;   //!
    TBranch        *bMpf;   //!
    TBranch        *bRjj;   //!
    tr->SetBranchAddress("JetsPt",&JetsPt,&bJetsPt);
    tr->SetBranchAddress("JetsEta",&JetsEta,&bJetsEta);
    tr->SetBranchAddress("JetsPhi",&JetsPhi,&bJetsPhi);
    tr->SetBranchAddress("nVertex", &nVertex, &bnVertex);
    tr->SetBranchAddress("nTracks", &nTracks, &bnTracks);
    tr->SetBranchAddress("MinDistance", &MinDistance, &bMinDistance);
    tr->SetBranchAddress("MaxDistance", &MaxDistance, &bMaxDistance);
    tr->SetBranchAddress("GoldenVertexZ", &GoldenVertexZ, &bGoldenVertexZ);
    tr->SetBranchAddress("xiPPSArmB", &xiPPSArmB, &bxiPPSArmB);
    tr->SetBranchAddress("xiPPSArmF", &xiPPSArmF, &bxiPPSArmF);
    tr->SetBranchAddress("tPPSArmB", &tPPSArmB, &btPPSArmB);
    tr->SetBranchAddress("tPPSArmF", &tPPSArmF, &btPPSArmF);
    tr->SetBranchAddress("xPPSArmBDet1", &xPPSArmBDet1, &bxPPSArmBDet1);
    tr->SetBranchAddress("xPPSArmFDet1", &xPPSArmFDet1, &bxPPSArmFDet1);
    tr->SetBranchAddress("yPPSArmBDet1", &yPPSArmBDet1, &byPPSArmBDet1);
    tr->SetBranchAddress("yPPSArmFDet1", &yPPSArmFDet1, &byPPSArmFDet1);
    tr->SetBranchAddress("xPPSArmBDet2", &xPPSArmBDet2, &bxPPSArmBDet2);
    tr->SetBranchAddress("xPPSArmFDet2", &xPPSArmFDet2, &bxPPSArmFDet2);
    tr->SetBranchAddress("yPPSArmBDet2", &yPPSArmBDet2, &byPPSArmBDet2);
    tr->SetBranchAddress("yPPSArmFDet2", &yPPSArmFDet2, &byPPSArmFDet2);
    tr->SetBranchAddress("stopPPSArmFTrkDet1", &stopPPSArmFTrkDet1, &bstopPPSArmFTrkDet1);
    tr->SetBranchAddress("stopPPSArmFTrkDet2", &stopPPSArmFTrkDet2, &bstopPPSArmFTrkDet2);
    tr->SetBranchAddress("stopPPSArmBTrkDet1", &stopPPSArmBTrkDet1, &bstopPPSArmBTrkDet1);
    tr->SetBranchAddress("stopPPSArmBTrkDet2", &stopPPSArmBTrkDet2, &bstopPPSArmBTrkDet2);
    tr->SetBranchAddress("stopPPSArmFToF", &stopPPSArmFToF, &bstopPPSArmFToF);
    tr->SetBranchAddress("stopPPSArmBToF", &stopPPSArmBToF, &bstopPPSArmBToF);
    tr->SetBranchAddress("Mjj", &Mjj, &bMjj);
    tr->SetBranchAddress("Mpf", &Mpf, &bMpf);
    tr->SetBranchAddress("Rjj", &Rjj, &bRjj);
    //----------- settings ---------------
    int NEVENTS(1000);
    double pTmin(100);
    double etamin(-2.5);
    double etamax(2.5);
    double MxPPS(0.0);
    double MJJ(0.0);
    double deltaeta(0.0);
    double deltaphi(0.0);
    double S(13000.0);
    double xmax(-3.1);
    double xmin(-23.1);
    double ymax( 9.0);
    double ymin(-9.0);

    bool verbose = false;

    //--------- book histos ---------------------------------
    TH1F *hLeadingJetPt = new TH1F("LeadingJetPt","Leading Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
    TH1F *hSecondJetPt = new TH1F("SecondJetPt","Second Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
    TH1F *hLeadingJetEta = new TH1F("LeadingJetEta","Leading Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
    TH1F *hSecondJetEta = new TH1F("SecondJetEta","Second Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
    TH1F *hLeadingJetPhi = new TH1F("LeadingJetPhi","Leading Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);
    TH1F *hSecondJetPhi = new TH1F("SecondJetPhi","Second Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);  
    TH1F *hDeltaEtaJets = new TH1F("DeltaEtaJets","#Delta#eta_{jj} Distribution; #Delta#eta_{jj}; N events",20,0.0,5.2);
    TH1F *hDeltaPhiJets = new TH1F("DeltaPhiJets","#Delta#phi_{jj} Distribution; #Delta#phi_{jj}; N events",20,0.0,1.2*PI);
    TH1F *hMjj = new TH1F( "Mjj" , "Mass_{JJ} Distribution; M_{jj}  [GeV]; N events" , 100, 0., 1000. );
    TH1F *hMx = new TH1F("Mx" , "Mass_{X} Distribution; M_{x}  [GeV]; N events" , 100, 0., 1000. );
    TH1F *hPPS_xiARMPlus =  new TH1F( "PPS_xiARMPlus" , "#xi_{plus} PPS; #xi_{plus}; N Event" , 100, 0., 0.4 );
    TH1F *hPPS_xiARMMinus =  new TH1F( "PPS_xiARMMinus" , "#xi_{minus} PPS; #xi_{minus}; N Event" , 100, 0., 0.4 ); 
    TH2F *hPPS_xVsy_ARMPlusDt1 =  new TH2F( "PPS_xVsy_ARMPlusDt1" , "PPS x_vs_y_{ARMPlusDt1}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 ); 
    TH2F *hPPS_xVsy_ARMPlusDt2 =  new TH2F( "PPS_xVsy_ARMPlusDt2" , "PPS x_vs_y_{ARMPlusDt2}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 );  
    TH2F *hPPS_xVsy_ARMMinusDt1 =  new TH2F( "PPS_xVsy_ARMMinusDt1" , "PPS x_vs_y_{ARMMinusDt1}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 ); 
    TH2F *hPPS_xVsy_ARMMinusDt2 =  new TH2F( "PPS_xVsy_ARMMinusDt2" , "PPS x_vs_y_{ARMMinusDt2}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 ); 
    TH2F *hMxPPSvsMjj =  new TH2F("MxvsMjj" , "Mass_{X} vs M_{JJ}  Distribution; M_{x}  [GeV];  M_{jj}  [GeV]" , 100, 0., 1000.,100, 0., 1000. );
    TH1F *hLeadingJetPtAfterPPS = new TH1F("LeadingJetPtAfterPPS","Leading Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
    TH1F *hSecondJetPtAfterPPS = new TH1F("SecondJetPtAfterPPS","Second Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
    TH1F *hLeadingJetEtaAfterPPS = new TH1F("LeadingJetEtaAfterPPS","Leading Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
    TH1F *hSecondJetEtaAfterPPS = new TH1F("SecondJetEtaAfterPPS","Second Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
    TH1F *hLeadingJetPhiAfterPPS = new TH1F("LeadingJetPhiAfterPPS","Leading Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);  
    TH1F *hSecondJetPhiAfterPPS = new TH1F("SecondJetPhiAfterPPS","Second Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);    
    TH1F *hDeltaEtaJetsAfterPPS = new TH1F("DeltaEtaJetsAfterPPS","#Delta#eta_{jj} Distribution; #Delta#eta_{jj}; N events",20,0.0,5.2);
    TH1F *hDeltaPhiJetsAfterPPS = new TH1F("DeltaPhiJetsAfterPPS","#Delta#phi_{jj} Distribution; #Delta#phi_{jj}; N events",20,0.0,1.2*PI);
    //----------- counters -----------------------
    int counter_jet(0);
    //----------- tree reading -------------------
    unsigned NEntries = tr->GetEntries();
    cout<<"Reading TREE: "<<NEntries<<" events"<<endl;
    int decade = 0;
    for(unsigned i=0;i<NEVENTS;i++) {
        //----------- progress report -------------
        double progress = 10.0*i/(1.0*NEVENTS);
        //int k = TMath::FloorNint(progress);
        Int_t k = TMath::FloorNint(progress);
        if (k > decade) 
            cout<<10*k<<" %"<<endl;
        decade = k;          
        //----------- read the event --------------
        tr->GetEntry(i);

        if(JetsPt->at(0)>pTmin && JetsPt->at(1)>pTmin){
            if((JetsEta->at(0)>etamin && JetsEta->at(0)<etamax)&&(JetsEta->at(1)>etamin && JetsEta->at(1)<etamax)){
                //------- fill some Jets  histograms ------------------
                hLeadingJetPt->Fill(JetsPt->at(0));
                hSecondJetPt->Fill(JetsPt->at(1));
                hLeadingJetEta->Fill(JetsEta->at(0));
                hSecondJetEta->Fill(JetsEta->at(1));
                hLeadingJetPhi->Fill(JetsPhi->at(0));
                hSecondJetPhi->Fill(JetsPhi->at(1));
                deltaeta = fabs(JetsEta->at(0)-JetsEta->at(1));
                deltaphi = fabs(JetsPhi->at(0)-JetsPhi->at(1));
                if (deltaphi> PI){deltaphi = 2.0*PI - deltaphi;}
                hDeltaEtaJets->Fill(deltaeta);
                hDeltaPhiJets->Fill(deltaphi);
                //--------- PPS x vs y -----------------------
                //x,det1 e det2 nos dois arms < - 3.1 e  > -23.1 mm
                //y,det1 e det2 nos dois arms > -9 e < 9 mm
                bool cutXdet1 = ((xmin<xPPSArmFDet1 && xPPSArmFDet1<xmax) && (xmin<xPPSArmBDet1 && xPPSArmBDet1<xmax));
                bool cutXdet2 = ((xmin<xPPSArmFDet2 && xPPSArmFDet2<xmax) && (xmin<xPPSArmBDet2  && xPPSArmBDet2<xmax));
                bool cutYdet1 = ((ymin<yPPSArmFDet1 && yPPSArmFDet1<ymax) && (ymin<yPPSArmBDet1 && yPPSArmBDet1<ymax));
                bool cutYdet2 = ((ymin<yPPSArmFDet2 && yPPSArmFDet2<ymax) && (ymin<yPPSArmBDet2 &&  yPPSArmBDet2<ymax));
                if(verbose)std::cout<<"Xdet1:" <<cutXdet1<<"Xdet2:"<<cutXdet2<<"Ydet1: "<<cutYdet1<<"Ydet2: "<<cutYdet2<<std::endl;
                //--------- PPS HasStoped -----------------------
                bool stopTrkDet1 = (stopPPSArmFTrkDet1==0 && stopPPSArmBTrkDet1 ==0);
                bool stopTrkDet2 = (stopPPSArmFTrkDet2==0 && stopPPSArmBTrkDet2 ==0);
                if(verbose)cout<<"stopTrkDet1:"<<stopTrkDet1<<" , "<<"stopTrkDet2:"<<stopTrkDet2<<endl;
                if(cutXdet1 && cutXdet2 && cutYdet1 && cutYdet2){
                    if(stopTrkDet1 && stopTrkDet2){ 
                        hPPS_xVsy_ARMPlusDt1->Fill(xPPSArmFDet1,yPPSArmFDet1);
                        hPPS_xVsy_ARMPlusDt2->Fill(xPPSArmFDet2,yPPSArmFDet2); 
                        hPPS_xVsy_ARMMinusDt1->Fill(xPPSArmBDet1,yPPSArmBDet1);
                        hPPS_xVsy_ARMMinusDt2->Fill(xPPSArmBDet2,yPPSArmBDet2);
                        //--------- PPS xi ----------------------------
                        if(xiPPSArmF>0.0 && xiPPSArmB>0.0){

                            hPPS_xiARMPlus->Fill(xiPPSArmF);
                            hPPS_xiARMMinus->Fill(xiPPSArmB);
                            MxPPS = S*sqrt(xiPPSArmF*xiPPSArmB);
                            if(verbose)std::cout<<"xiPPSArmF="<< xiPPSArmF<<","<<"xiPPSArmB="<<xiPPSArmB<<", "<<"S= "<<S<<std::endl;
                            if(verbose)std::cout<<"Mx = "<<sqrt(xiPPSArmF*xiPPSArmB)*S<<std::endl;
                            MJJ = Mjj; 
                            hLeadingJetPtAfterPPS->Fill(JetsPt->at(0));
                            hSecondJetPtAfterPPS->Fill(JetsPt->at(1));
                            hLeadingJetEtaAfterPPS->Fill(JetsEta->at(0));
                            hSecondJetEtaAfterPPS->Fill(JetsEta->at(1));
                            hLeadingJetPhiAfterPPS->Fill(JetsPhi->at(0));
                            hSecondJetPhiAfterPPS->Fill(JetsPhi->at(1));
                            hDeltaEtaJetsAfterPPS->Fill(deltaeta);
                            hDeltaPhiJetsAfterPPS->Fill(deltaphi);
                            //--------------------------------------------
                        } // xi cut
                        hMx->Fill(MxPPS);
                        hMjj->Fill(MJJ);
                        hMxPPSvsMjj->Fill(MxPPS,MJJ);
                    } //stopTrk
                } //x,y cut
              }//eta jet cuts
            }  //jet pt cut

        }// tree loop
        //----------------- print out some information ---------------
        cout<<"Events read:                      "<<NEVENTS<<endl;
        //----------------- save the histos to the output file ------
        outf->Write();
    }

