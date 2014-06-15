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
Bool_t switchlumiweight = true;
void CEPDijetsPPS_JetsSameVertex()
{
	//------------ FWLite libraries ------
	gSystem->Load("libFWCoreFWLite.so");
	AutoLibraryLoader::enable();
	gROOT->ProcessLine("#include<vector>");
	gROOT->ProcessLine(".exception");
	gROOT->ProcessLine(".L styleTDR.C");
	gROOT->ProcessLine("setTDRStyle()");

	//------------ files -----------------
	//TFile *inf  = TFile::Open("ttreeCEPdijetsNoOOT_NoPU.root");
	TFile *inf  = TFile::Open("ttreeCEPdijetsNoOOT_PU.root");
	//TFile *inf  = TFile::Open("ttreeCEPdijetsOOT_PU.root");
	TFile *outf = new TFile("CEPDijets_PPS_NoOOT_PU_HistosTEST.root","RECREATE");

	TTree *tr = (TTree*)inf->Get("demo/Event");
	//----------- define the tree branch --------
	std::vector<double>  *JetsPt=0;
	std::vector<double>  *JetsEta=0;
	std::vector<double>  *JetsPhi=0;
	std::vector<double>  *JetsSameVertex_pt=0;
	std::vector<double>  *JetsSameVertex_eta=0;
	std::vector<double>  *JetsSameVertex_phi=0;
	std::vector<double>  *CandidatesJets_pt=0;
	std::vector<double>  *CandidatesJets_eta=0;
	std::vector<double>  *CandidatesJets_phi=0; 
	// xi PPS
	std::vector<double>  *xiPPSArmFInfo=0;
	std::vector<double>  *xiPPSArmBInfo=0;   
	//t PPS
	std::vector<double>  *tPPSArmFInfo=0;     
	std::vector<double>  *tPPSArmBInfo=0;
	// Arm Forward det1 and det2
	std::vector<double> *xPPSArmFDet1Info=0;
	std::vector<double> *yPPSArmFDet1Info=0;

	std::vector<double> *xPPSArmFDet2Info=0;
	std::vector<double> *yPPSArmFDet2Info=0;

	// Arm Backward det1 and det2
	std::vector<double> *xPPSArmBDet1Info=0;
	std::vector<double> *yPPSArmBDet1Info=0;

	std::vector<double> *xPPSArmBDet2Info=0;
	std::vector<double> *yPPSArmBDet2Info=0;

	// ToF Informations for Forward arm
	std::vector<double> *xPPSArmFToFInfo=0;
	std::vector<double> *yPPSArmFToFInfo=0;

	// Backward arm
	std::vector<double> *xPPSArmBToFInfo=0;
	std::vector<double> *yPPSArmBToFInfo=0;
       //Stopped ToF      
        std::vector<double>  *stopPPSArmFToFInfo=0;
       std::vector<double>  *stopPPSArmBToFInfo=0;
   //Track 
       std::vector<double>  *stopPPSArmFTrkDet1Info=0;
       std::vector<double>  *stopPPSArmFTrkDet2Info=0;
       std::vector<double>  *stopPPSArmBTrkDet1Info=0;
       std::vector<double>  *stopPPSArmBTrkDet2Info=0;


        


	Double_t        VertexZPPS;
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
	Double_t        CandidatesMjj;
	Double_t        Mpf;
	Double_t        Rjj;

	// ToF detector
	Double_t deltaToF_00;
	Double_t deltaToF_0i;
	Double_t deltaToF_ii;


          //Stopped ToF      
        TBranch  *bstopPPSArmFToFInfo=0;
       TBranch  *bstopPPSArmBToFInfo=0;
   //Track 
       TBranch  *bstopPPSArmFTrkDet1Info=0;
       TBranch  *bstopPPSArmFTrkDet2Info=0;
       TBranch  *bstopPPSArmBTrkDet1Info=0;
       TBranch  *bstopPPSArmBTrkDet2Info=0;





	TBranch        *bJetsPt=0;   //!
	TBranch        *bJetsEta=0;   //!
	TBranch        *bJetsPhi=0;   //!
	TBranch        *bJetsSameVertex_pt=0 ;
	TBranch        *bJetsSameVertex_eta=0;   
	TBranch        *bJetsSameVertex_phi=0;
	TBranch        *bCandidatesJets_pt=0;
	TBranch        *bCandidatesJets_eta=0;
	TBranch        *bCandidatesJets_phi=0;

	//PPS vector branches
	TBranch *bxiPPSArmFInfo=0;
	TBranch *bxiPPSArmBInfo=0;

	TBranch *btPPSArmFInfo=0;
	TBranch *btPPSArmBInfo=0;
	// Arm Forward det1 and det2
	TBranch *bxPPSArmFDet1Info=0;
	TBranch *byPPSArmFDet1Info=0;

	TBranch *bxPPSArmFDet2Info=0;
	TBranch *byPPSArmFDet2Info=0;

	// Arm Backward det1 and det2
	TBranch *bxPPSArmBDet1Info=0;
	TBranch *byPPSArmBDet1Info=0;

	TBranch *bxPPSArmBDet2Info=0;
	TBranch *byPPSArmBDet2Info=0;

	// ToF Informations for Forward arm
	TBranch *bxPPSArmFToFInfo=0;
	TBranch *byPPSArmFToFInfo=0;

	// Backward arm
	TBranch *bxPPSArmBToFInfo=0;
	TBranch *byPPSArmBToFInfo=0;

	//Tof
	TBranch *bdeltaToF_00;
	TBranch *bdeltaToF_0i;
	TBranch *bdeltaToF_ii;

	TBranch        *bVertexZPPS; 
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
	TBranch        *bCandidatesMjj=0;
	TBranch        *bMpf;   //!
	TBranch        *bRjj;   //!


	tr->SetBranchAddress("JetsPt", &JetsPt, &bJetsPt);
	tr->SetBranchAddress("JetsEta", &JetsEta, &bJetsEta);
	tr->SetBranchAddress("JetsPhi", &JetsPhi, &bJetsPhi);
	tr->SetBranchAddress("JetsSameVertex_pt",&JetsSameVertex_pt,&bJetsSameVertex_pt);
	tr->SetBranchAddress("JetsSameVertex_eta",&JetsSameVertex_eta,&bJetsSameVertex_eta);
	tr->SetBranchAddress("JetsSameVertex_phi",&JetsSameVertex_phi,&bJetsSameVertex_phi);
	tr->SetBranchAddress("CandidatesJets_pt",&CandidatesJets_pt,&bCandidatesJets_pt);
	tr->SetBranchAddress("CandidatesJets_eta",&CandidatesJets_eta,&bCandidatesJets_eta);  
	tr->SetBranchAddress("CandidatesJets_phi",&CandidatesJets_phi,&bCandidatesJets_phi);  
	tr->SetBranchAddress("VertexZPPS", &VertexZPPS, &bVertexZPPS);
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
	tr->SetBranchAddress("CandidatesMjj",&CandidatesMjj,&bCandidatesMjj);  
	tr->SetBranchAddress("Mjj", &Mjj, &bMjj);
	tr->SetBranchAddress("Mpf", &Mpf, &bMpf);
	tr->SetBranchAddress("Rjj", &Rjj, &bRjj);
	// xi PPS
	tr->SetBranchAddress("xiPPSArmFInfo", &xiPPSArmFInfo, &bxiPPSArmFInfo);	
	tr->SetBranchAddress("xiPPSArmBInfo", &xiPPSArmBInfo, &bxiPPSArmBInfo);
	// t PPS
	tr->SetBranchAddress("tPPSArmFInfo", &tPPSArmFInfo, &btPPSArmFInfo);
	tr->SetBranchAddress("tPPSArmBInfo", &tPPSArmBInfo, &btPPSArmBInfo); 
	// Arm Fwd det1
	tr->SetBranchAddress("xPPSArmFDet1Info", &xPPSArmFDet1Info, &bxPPSArmFDet1Info);
	tr->SetBranchAddress("yPPSArmFDet1Info", &yPPSArmFDet1Info, &byPPSArmFDet1Info);
	// Arm Fwd Det2
	tr->SetBranchAddress("xPPSArmFDet2Info", &xPPSArmFDet2Info, &bxPPSArmFDet2Info);
	tr->SetBranchAddress("yPPSArmFDet2Info", &yPPSArmFDet2Info, &byPPSArmFDet2Info);
	// Arm Bk Det1
	tr->SetBranchAddress("xPPSArmBDet1Info", &xPPSArmBDet1Info, &bxPPSArmBDet1Info);
	tr->SetBranchAddress("yPPSArmBDet1Info", &yPPSArmBDet1Info, &byPPSArmBDet1Info);
	// Arm Bk Det2
	tr->SetBranchAddress("xPPSArmBDet2Info", &xPPSArmBDet2Info, &bxPPSArmBDet2Info);
	tr->SetBranchAddress("yPPSArmBDet2Info", &yPPSArmBDet2Info, &byPPSArmBDet2Info);

	//Delta ToF 
	tr->SetBranchAddress("deltaToF_00", &deltaToF_00, &bdeltaToF_00);
	tr->SetBranchAddress("deltaToF_0i", &deltaToF_0i, &bdeltaToF_0i);
	tr->SetBranchAddress("deltaToF_ii", &deltaToF_ii, &bdeltaToF_ii);  

        // Stopped 
         tr->SetBranchAddress("stopPPSArmFToFInfo", &stopPPSArmFToFInfo, &bstopPPSArmFToFInfo);
         tr->SetBranchAddress("stopPPSArmBToFInfo", &stopPPSArmBToFInfo, &bstopPPSArmBToFInfo);


       // Track
         tr->SetBranchAddress("stopPPSArmFTrkDet1Info", &stopPPSArmFTrkDet1Info, &bstopPPSArmFTrkDet1Info);
         tr->SetBranchAddress("stopPPSArmBTrkDet1Info", &stopPPSArmBTrkDet1Info, &bstopPPSArmBTrkDet1Info);      

         tr->SetBranchAddress("stopPPSArmFTrkDet2Info", &stopPPSArmFTrkDet2Info, &bstopPPSArmFTrkDet2Info);
         tr->SetBranchAddress("stopPPSArmBTrkDet2Info", &stopPPSArmBTrkDet2Info, &bstopPPSArmBTrkDet2Info);  
        


	//---------aa settings ---------------
	bool verbose = false;
	int NEVENTS = tr->GetEntries();
	cout<<"MC EVENTS= "<<NEVENTS<<endl;
	float cross_section =1700.0 ; //[fb]
	float luminosity=100.0;//[fb]^-1
	double lumiweight = (luminosity*cross_section)/NEVENTS;
	cout<<"lumiweight= "<<lumiweight<<endl;
	double mclumiweight = 1.0;
	if (switchlumiweight ){
		mclumiweight = lumiweight;
	}

	int Njets=0;
	double pTmin  = 50.0;
	double etamin = -2.0;
	double etamax = 2.0;

	double MxPPS = -999.0;
	double MxPPS_PU = -999.0;

	double MJJ = -999.0;
	double MJJ_PU = -999.0;        

	double MxPPSBeforeCuts = 0.0;
	double MxPPSBeforeCuts_PU = 0.0;

	double MJJBC = -999.0;
	double MJJBC_PU = -999.0;

	double deltaeta = 0.0;
	double deltaphi = 0.0;
	double S = 13000.0;
	double xmax = -3.15;
	double xmin =-23.15;
	double xmax2 = -2.03;
	double xmin2 = -22.03; 
	double ymax = 9.0;
	double ymin =-9.0;


	//--------- book histos ---------------------------------
	TH1F *hNJets = new TH1F("NJets","N Jets;  N Jets; N events",100,0,100);
	TH1F *hVertexZCMS = new TH1F("VertexZCMS"," Vertex Z CMS[cm];  Vertex Z [cm]; N events",25,-25.0,25.0);
	TH1F *hVertexZPPS = new TH1F("VertexZPPS","Vertex Z PPS [cm]; Vertex  Z [cm]; N events",25,-25.0,25.0);
	TH2F *hVertexZCMSPPS = new TH2F("VertexZCMSPPS","Vertex Z CMS vs Vertex Z  PPS; Vertex Z CMS [cm]; Vertex Z PPS [cm]",25,-25.0,25.0,25, -25.0,25.0);
	TH1F *hLeadingJetPt = new TH1F("LeadingJetPt","Leading Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
	TH1F *hSecondJetPt = new TH1F("SecondJetPt","Second Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
	TH1F *hLeadingJetEta = new TH1F("LeadingJetEta","Leading Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
	TH1F *hSecondJetEta = new TH1F("SecondJetEta","Second Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
	TH1F *hLeadingJetPhi = new TH1F("LeadingJetPhi","Leading Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);
	TH1F *hSecondJetPhi = new TH1F("SecondJetPhi","Second Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);  
	TH1F *hDeltaEtaJets = new TH1F("DeltaEtaJets","#Delta#eta_{jj} Distribution; #Delta#eta_{jj}; N events",20,0.0,5.2);
	TH1F *hDeltaPhiJets = new TH1F("DeltaPhiJets","#Delta#phi_{jj} Distribution; #Delta#phi_{jj}; N events",20,0.0,1.2*PI);
	TH1F *hMjj = new TH1F( "Mjj" , "Mass_{JJ} Distribution; M_{jj}  [GeV]; N events" , 100, 0., 2000. );
	TH1F *hMx = new TH1F("Mx" , "Mass_{X} Distribution; M_{x}  [GeV]; N events" , 100, 0., 2000. );
	// Signal
	TH1F *hMjjBC = new TH1F( "MjjBC" , "Mass_{JJ} Distribution; M_{jj}  [GeV]; N events" , 100, 0., 2000. );
	TH1F *hMxBC = new TH1F("MxBC" , "Mass_{X} Distribution; M_{x}  [GeV]; N events" , 100, 0., 2000. );
	TH1F *hPPS_xiARMPlus =  new TH1F( "PPS_xiARMPlus" , "#xi_{plus} PPS; #xi_{plus}; N Event" , 100, 0., 0.4 );
	TH1F *hPPS_xiARMMinus =  new TH1F( "PPS_xiARMMinus" , "#xi_{minus} PPS; #xi_{minus}; N Event" , 100, 0., 0.4 ); 
	TH2F *hPPS_xVsy_ARMPlusDt1 =  new TH2F( "PPS_xVsy_ARMPlusDt1" , "PPS x_vs_y_{ARMPlusDt1}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 ); 
	TH2F *hPPS_xVsy_ARMPlusDt2 =  new TH2F( "PPS_xVsy_ARMPlusDt2" , "PPS x_vs_y_{ARMPlusDt2}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 );  
	TH2F *hPPS_xVsy_ARMMinusDt1 =  new TH2F( "PPS_xVsy_ARMMinusDt1" , "PPS x_vs_y_{ARMMinusDt1}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 ); 
	TH2F *hPPS_xVsy_ARMMinusDt2 =  new TH2F( "PPS_xVsy_ARMMinusDt2" , "PPS x_vs_y_{ARMMinusDt2}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 ); 
	TH2F *hMxPPSvsMjj =  new TH2F("MxvsMjj" , "Mass_{X} vs M_{JJ}  Distribution; M_{x}  [GeV];  M_{jj}  [GeV]" , 100, 0., 2000.,100, 0., 2000. );
	TH2F *hMxPPSvsMjjBC =  new TH2F("MxvsMjjBC" , "Mass_{X} vs M_{JJ}  Distribution; M_{x}  [GeV];  M_{jj}  [GeV]" , 100, 0., 2000.,100, 0., 2000. );
	TH1F *hLeadingJetPtAfterPPS = new TH1F("LeadingJetPtAfterPPS","Leading Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
	TH1F *hSecondJetPtAfterPPS = new TH1F("SecondJetPtAfterPPS","Second Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
	TH1F *hLeadingJetEtaAfterPPS = new TH1F("LeadingJetEtaAfterPPS","Leading Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
	TH1F *hSecondJetEtaAfterPPS = new TH1F("SecondJetEtaAfterPPS","Second Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
	TH1F *hLeadingJetPhiAfterPPS = new TH1F("LeadingJetPhiAfterPPS","Leading Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);  
	TH1F *hSecondJetPhiAfterPPS = new TH1F("SecondJetPhiAfterPPS","Second Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);    
	TH1F *hDeltaEtaJetsAfterPPS = new TH1F("DeltaEtaJetsAfterPPS","#Delta#eta_{jj} Distribution; #Delta#eta_{jj}; N events",20,0.0,5.2);
	TH1F *hDeltaPhiJetsAfterPPS = new TH1F("DeltaPhiJetsAfterPPS","#Delta#phi_{jj} Distribution; #Delta#phi_{jj}; N events",20,0.0,1.2*PI);
        TH1F *hPPS_tARMPlus =  new TH1F( "PPS_tARMPlus" , "t_{plus} PPS; t_{plus}; N Event" , 200, 0., 20. );
        TH1F *hPPS_tARMMinus =  new TH1F( "PPS_tARMMinus" , "t_{minus} PPS; t_{minus}; N Event" , 200, 0., 20. );
	// PU Contribution
        TH1F *hMjj_PU = new TH1F( "Mjj_PU" , "Mass_{JJ} Distribution; M_{jj}  [GeV]; N events" , 100, 0., 2000. );
	TH1F *hMx_PU = new TH1F("Mx_PU" , "Mass_{X} Distribution; M_{x}  [GeV]; N events" , 100, 0., 2000. );
	TH1F *hMjjBC_PU = new TH1F( "MjjBC_PU" , "Mass_{JJ} Distribution; M_{jj}  [GeV]; N events " , 100, 0., 2000. );
	TH1F *hMxBC_PU = new TH1F("MxBC_PU" , "Mass_{X} Distribution; M_{x}  [GeV]; N events " , 100, 0., 2000. );
	TH1F *hPPS_xiARMPlus_PU =  new TH1F( "PPS_xiARMPlus_PU" , "#xi_{plus} PPS; #xi_{plus}; N Event" , 100, 0., 0.4 );
	TH1F *hPPS_xiARMMinus_PU =  new TH1F( "PPS_xiARMMinus_PU" , "#xi_{minus} PPS; #xi_{minus}; N Event" , 100, 0., 0.4 );
	TH2F *hPPS_xVsy_ARMPlusDt1_PU =  new TH2F( "PPS_xVsy_ARMPlusDt1_PU" , "PPS x_vs_y_{ARMPlusDt1}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 );
	TH2F *hPPS_xVsy_ARMPlusDt2_PU =  new TH2F( "PPS_xVsy_ARMPlusDt2_PU" , "PPS x_vs_y_{ARMPlusDt2}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 );
	TH2F *hPPS_xVsy_ARMMinusDt1_PU =  new TH2F( "PPS_xVsy_ARMMinusDt1_PU" , "PPS x_vs_y_{ARMMinusDt1}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 );
	TH2F *hPPS_xVsy_ARMMinusDt2_PU =  new TH2F( "PPS_xVsy_ARMMinusDt2_PU" , "PPS x_vs_y_{ARMMinusDt2}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 );
	TH2F *hMxPPSvsMjj_PU =  new TH2F("MxvsMjj_PU" , "Mass_{X} vs M_{JJ}  Distribution; M_{x}  [GeV];  M_{jj}  [GeV]" , 100, 0., 2000.,100, 0., 2000. );
	TH2F *hMxPPSvsMjjBC_PU =  new TH2F("MxvsMjjBC_PU" , "Mass_{X} vs M_{JJ}  Distribution; M_{x}  [GeV];  M_{jj}  [GeV]" , 100, 0., 2000.,100, 0., 2000. );
	TH1F *hLeadingJetPtAfterPPS_PU = new TH1F("LeadingJetPtAfterPPS_PU","Leading Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
	TH1F *hSecondJetPtAfterPPS_PU = new TH1F("SecondJetPtAfterPPS_PU","Second Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; N events",100,0,500);
	TH1F *hLeadingJetEtaAfterPPS_PU = new TH1F("LeadingJetEtaAfterPPS_PU","Leading Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
	TH1F *hSecondJetEtaAfterPPS_PU = new TH1F("SecondJetEtaAfterPPS_PU","Second Jet - #eta Distribution; #eta; N events",100,-5.2,5.2);
	TH1F *hLeadingJetPhiAfterPPS_PU = new TH1F("LeadingJetPhiAfterPPS_PU","Leading Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);
	TH1F *hSecondJetPhiAfterPPS_PU = new TH1F("SecondJetPhiAfterPPS_PU","Second Jet - #phi Distribution; #phi; N events",100,-1.2*PI,1.2*PI);
	TH1F *hDeltaEtaJetsAfterPPS_PU = new TH1F("DeltaEtaJetsAfterPPS_PU","#Delta#eta_{jj} Distribution; #Delta#eta_{jj}; N events",20,0.0,5.2);
	TH1F *hDeltaPhiJetsAfterPPS_PU = new TH1F("DeltaPhiJetsAfterPPS_PU","#Delta#phi_{jj} Distribution; #Delta#phi_{jj}; N events",20,0.0,1.2*PI);
     
        TH1F *hPPS_tARMPlus_PU =  new TH1F( "PPS_tARMPlus_PU" , "t_{plus} PPS; t_{plus}; N Event" , 200, 0., 20. );
        TH1F *hPPS_tARMMinus_PU =  new TH1F( "PPS_tARMMinus_PU" , "t_{minus} PPS; t_{minus}; N Event" , 200, 0., 20. );       



	//----------- counters -----------------------
	int counter_jetpT(0),counter_jetEta(0), counter_fiducial(0), counter_hasStoped(0),counter_xiArms(0);
	int counter_jetpT_PU(0),counter_jetEta_PU(0),counter_fiducial_PU(0), counter_hasStoped_PU(0),counter_xiArms_PU(0);
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
		Njets = JetsSameVertex_pt->size();
		//hNJets->Fill(Njets);

		if (JetsSameVertex_pt->size()<2)continue;
		if(verbose)cout<<"JetsSameVertex_pt->size ="<<JetsSameVertex_pt->size()<<endl;


		//------- fill some Jets  histograms ------------------
		hLeadingJetPt->Fill(JetsSameVertex_pt->at(0),mclumiweight);
		hSecondJetPt->Fill(JetsSameVertex_pt->at(1),mclumiweight);
		hLeadingJetEta->Fill(JetsSameVertex_eta->at(0),mclumiweight);
		hSecondJetEta->Fill(JetsSameVertex_eta->at(1),mclumiweight);
		hLeadingJetPhi->Fill(JetsSameVertex_phi->at(0),mclumiweight);
		hSecondJetPhi->Fill(JetsSameVertex_phi->at(1),mclumiweight);
		deltaeta = fabs(JetsSameVertex_eta->at(0)-JetsSameVertex_eta->at(1));
		deltaphi = fabs(JetsSameVertex_phi->at(0)-JetsSameVertex_phi->at(1));
		if (deltaphi> PI){deltaphi = 2.0*PI - deltaphi;}
		hDeltaEtaJets->Fill(deltaeta,mclumiweight);
		hDeltaPhiJets->Fill(deltaphi),mclumiweight;

		hVertexZCMSPPS->Fill(GoldenVertexZ,VertexZPPS,mclumiweight);        
		hVertexZCMS->Fill(GoldenVertexZ,mclumiweight);
		hVertexZPPS->Fill(VertexZPPS,mclumiweight);

		if(verbose)std::cout <<"xPPSArmFDet1Info->size(): " <<xPPSArmFDet1Info->size() << std::endl;
		if(verbose)std::cout <<"xPPSArmFDet2Info->size(): " <<xPPSArmFDet2Info->size() << std::endl;

		if(verbose)std::cout <<"xPPSArmBDet1Info->size(): " <<xPPSArmBDet1Info->size() << std::endl;
		if(verbose)std::cout <<"xPPSArmBDet2Info->size(): " <<xPPSArmBDet2Info->size() << std::endl;

		if(verbose)std::cout <<"yPPSArmFDet1Info->size(): " <<yPPSArmFDet1Info->size() << std::endl;
		if(verbose)std::cout <<"yPPSArmFDet2Info->size(): " <<yPPSArmFDet2Info->size() << std::endl;

		if(verbose)std::cout <<"yPPSArmBDet1Info->size(): " <<yPPSArmBDet1Info->size() << std::endl;
		if(verbose)std::cout <<"yPPSArmBDet2Info->size(): " <<yPPSArmBDet2Info->size() << std::endl;

		if(verbose)std::cout <<  "xiPPSArmFInfo->size(): " << xiPPSArmFInfo->size() << std::endl;
		if(verbose) std::cout <<  "xiPPSArmBInfo->size(): " << xiPPSArmBInfo->size() << std::endl;


		// Applying Fiducial Cuts for Signal index = 0
		//--------- PPS x vs y -----------------------
		//x,det1 e det2 nos dois arms < - 3.1 e  > -23.1 mm
		//y,det1 e det2 nos dois arms > -9 e < 9 mm
		bool cutXArmF = ((xmin<xPPSArmFDet1Info->at(0) && xPPSArmFDet1Info->at(0)<xmax) && (xmin2< xPPSArmFDet2Info->at(0) && xPPSArmFDet2Info->at(0)<xmax2));
		bool cutXArmB = ((xmin<xPPSArmBDet1Info->at(0) && xPPSArmBDet1Info->at(0)<xmax) && (xmin2<xPPSArmBDet2Info->at(0) && xPPSArmBDet2Info->at(0)<xmax2));
		bool cutYArmF = ((ymin<yPPSArmFDet1Info->at(0) && yPPSArmFDet1Info->at(0)<ymax) && (ymin<yPPSArmFDet2Info->at(0) && yPPSArmFDet2Info->at(0)<ymax));
		bool cutYArmB = ((ymin<yPPSArmBDet1Info->at(0) && yPPSArmBDet1Info->at(0)<ymax) && (ymin<yPPSArmBDet2Info->at(0) &&  yPPSArmBDet2Info->at(0)<ymax));
		if(verbose)std::cout<<"XArmF:" <<cutXArmF<<"XArmB:"<<cutXArmB<<"YArmF: "<<cutXArmF<<"YArmB: "<<cutYArmB<<std::endl;
		//std::cout<<"XArmF:" <<cutXArmF<<"XArmB:"<<cutXArmB<<"YArmF: "<<cutXArmF<<"YArmB: "<<cutYArmB<<std::endl; 
		//--------- PPS noHasStoped ----------------------- //
		bool stopTrkArmF = (stopPPSArmFTrkDet1Info->at(0)==0 && stopPPSArmFTrkDet2Info->at(0)==0);//Falta
		bool stopTrkArmB = (stopPPSArmBTrkDet1Info->at(0) ==0 && stopPPSArmBTrkDet2Info->at(0) ==0); //Falta
		if(verbose)cout<<"stopTrkArmF:"<<stopTrkArmF<<" , "<<"stopTrkArmB"<<stopTrkArmB<<endl;

		// Mjj and Mx before PPS sel
		if((xiPPSArmFInfo->at(0)>0.0) && (xiPPSArmBInfo->at(0)>0.0)){
			MxPPSBeforeCuts = S*sqrt((xiPPSArmFInfo->at(0))*(xiPPSArmBInfo->at(0)));
			if(verbose)cout<<"Mx BF = "<<S*sqrt((xiPPSArmFInfo->at(0))*(xiPPSArmBInfo->at(0)))<<endl;
			MJJBC = Mjj;   
		}
		hMxBC->Fill(MxPPSBeforeCuts,mclumiweight);
		hMjjBC->Fill(MJJBC,mclumiweight);                                                                                        
		hMxPPSvsMjjBC->Fill(MxPPS,MJJ,mclumiweight);


		//PPS sel
		if(cutXArmF && cutYArmF && cutXArmB && cutYArmB){
			++counter_fiducial;
			if(verbose)cout<<"x: "<<xPPSArmFDet1Info->at(0)<<", "<<xPPSArmBDet1Info->at(0)<<", "<<xPPSArmBDet1Info->at(0)<<", "<<xPPSArmBDet2Info->at(0)<<endl;
			if(verbose)cout<<"y: "<<yPPSArmFDet1Info->at(0)<<", "<<yPPSArmBDet1Info->at(0)<<", "<<yPPSArmBDet1Info->at(0)<<", "<<yPPSArmBDet2Info->at(0)<<endl;  
			if(stopTrkArmF && stopTrkArmB){
				++counter_hasStoped;
				hPPS_xVsy_ARMPlusDt1->Fill(xPPSArmFDet1Info->at(0),yPPSArmFDet1Info->at(0),mclumiweight);
				hPPS_xVsy_ARMPlusDt2->Fill(xPPSArmFDet2Info->at(0),yPPSArmFDet2Info->at(0),mclumiweight); 
				hPPS_xVsy_ARMMinusDt1->Fill(xPPSArmBDet1Info->at(0),yPPSArmBDet1Info->at(0),mclumiweight);
				hPPS_xVsy_ARMMinusDt2->Fill(xPPSArmBDet2Info->at(0),yPPSArmBDet2Info->at(0),mclumiweight);

				//--------- PPS xi ----------------------------

				if(xiPPSArmFInfo->at(0)>0.0 && xiPPSArmBInfo->at(0)>0.0){
					++counter_xiArms;
					if(JetsSameVertex_pt->at(0)>pTmin && JetsSameVertex_pt->at(1)>pTmin){
						++counter_jetpT; 

						if((JetsSameVertex_eta->at(0)>etamin && JetsSameVertex_eta->at(0)<etamax)&&(JetsSameVertex_eta->at(1)>etamin && JetsSameVertex_eta->at(1)<etamax)){
							++counter_jetEta;
							//		hPPS_xiARMPlus->Fill(xiPPSArmF,mclumiweight);
							hPPS_xiARMPlus->Fill(xiPPSArmFInfo->at(0),mclumiweight);  
							hPPS_tARMPlus->Fill(tPPSArmFInfo->at(0),mclumiweight);
							//hPPS_xiARMMinus->Fill(xiPPSArmB,mclumiweight);
							hPPS_xiARMMinus->Fill(xiPPSArmBInfo->at(0),mclumiweight);
							hPPS_tARMMinus->Fill(tPPSArmBInfo->at(0),mclumiweight);
							MxPPS = S*sqrt((xiPPSArmFInfo->at(0))*(xiPPSArmBInfo->at(0)));
							if(verbose)std::cout<<"xiPPSArmFInfo (signal)="<< xiPPSArmFInfo->at(0)<<","<<"xiPPSArmBInfo (Signal)="<<xiPPSArmBInfo->at(0)<<", "<<"S= "<<S<<std::endl;
							if(verbose)std::cout<<"Mx = "<<S*sqrt((xiPPSArmFInfo->at(0))*(xiPPSArmBInfo->at(0)))<<std::endl;
							MJJ = CandidatesMjj; 
							hLeadingJetPtAfterPPS->Fill(JetsSameVertex_pt->at(0),mclumiweight);
							hSecondJetPtAfterPPS->Fill(JetsSameVertex_pt->at(1),mclumiweight);
							hLeadingJetEtaAfterPPS->Fill(JetsSameVertex_eta->at(0),mclumiweight);
							hSecondJetEtaAfterPPS->Fill(JetsSameVertex_eta->at(1),mclumiweight);
							hLeadingJetPhiAfterPPS->Fill(JetsSameVertex_phi->at(0),mclumiweight);
							hSecondJetPhiAfterPPS->Fill(JetsSameVertex_phi->at(1),mclumiweight);
							hDeltaEtaJetsAfterPPS->Fill(deltaeta,mclumiweight);
							hDeltaPhiJetsAfterPPS->Fill(deltaphi,mclumiweight);
							//--------------------------------------------
						} // xi cut
						hMx->Fill(MxPPS,mclumiweight);
						if(verbose)std::cout<<"Mx_b = "<<S*sqrt((xiPPSArmFInfo->at(0))*(xiPPSArmBInfo->at(0)))<<std::endl;
						hMjj->Fill(MJJ,mclumiweight); 
						if(verbose)std::cout<<"Mjj = "<<Mjj<<" "<<"MJJ = "<<Mjj<<std::endl; //REver Dilson
						hMxPPSvsMjj->Fill(MxPPS,MJJ,mclumiweight);
					} //stopTrk
				} //x,y cut

			   }//eta jet cuts
			} //jet pt cut

	


		

		// PU background contributions 

		for( int k=1; k<xPPSArmFDet1Info->size(); k++ ){
		         //   if (xPPSArmFDet1Info->at(k)<!=0) continue;
                         //   if (yPPSArmFDet1Info->at(k)<!=0) continue;
                    

		//--------- PPS x vs y -----------------------
		//x,det1 e det2 nos dois arms < - 3.1 e  > -23.1 mm
		//y,det1 e det2 nos dois arms > -9 e < 9 mm
		bool cutXArmF_PU = ((xmin<xPPSArmFDet1Info->at(k) && xPPSArmFDet1Info->at(k)<xmax) && (xmin2< xPPSArmFDet2Info->at(k) && xPPSArmFDet2Info->at(k)<xmax2));
		bool cutXArmB_PU = ((xmin<xPPSArmBDet1Info->at(k) && xPPSArmBDet1Info->at(k)<xmax) && (xmin2<xPPSArmBDet2Info->at(k) && xPPSArmBDet2Info->at(k)<xmax2));
		bool cutYArmF_PU = ((ymin<yPPSArmFDet1Info->at(k) && yPPSArmFDet1Info->at(k)<ymax) && (ymin<yPPSArmFDet2Info->at(k) && yPPSArmFDet2Info->at(k)<ymax));
		bool cutYArmB_PU = ((ymin<yPPSArmBDet1Info->at(k) && yPPSArmBDet1Info->at(k)<ymax) && (ymin<yPPSArmBDet2Info->at(k) &&  yPPSArmBDet2Info->at(k)<ymax));
		if(verbose)std::cout<<"XArmF_PU:" <<cutXArmF_PU<<"XArmB_PU:"<<cutXArmB_PU<<"YArmF_PU: "<<cutXArmF_PU<<"YArmB_PU: "<<cutYArmB_PU<<std::endl;
		//std::cout<<"XArmF:" <<cutXArmF<<"XArmB:"<<cutXArmB<<"YArmF: "<<cutXArmF<<"YArmB: "<<cutYArmB<<std::endl; 
		//--------- PPS noHasStoped ----------------------- // rever com o Dilson: Precisas ser um vetor
		bool stopTrkArmF_PU =  (stopPPSArmFTrkDet1Info->at(k)==0 && stopPPSArmFTrkDet2Info->at(k)==0);
		bool stopTrkArmB_PU =  (stopPPSArmFTrkDet1Info->at(k)==0 && stopPPSArmFTrkDet2Info->at(k)==0);
		if(verbose)cout<<"stopTrkArmF_PU:"<<stopTrkArmF_PU<<" , "<<"stopTrkArmB_PU"<<stopTrkArmB_PU<<endl;

		// Mjj and Mx before PPS sel
		if((xiPPSArmFInfo->at(k)>0.0) && (xiPPSArmBInfo->at(k)>0.0)){
		MxPPSBeforeCuts_PU = S*sqrt((xiPPSArmFInfo->at(k))*(xiPPSArmBInfo->at(k)));
		if(verbose)cout<<"Mx BF PU = "<<S*sqrt((xiPPSArmFInfo->at(k))*(xiPPSArmBInfo->at(k)))<<endl;
		MJJBC_PU = Mjj; // Jet all
		}
		hMxBC_PU->Fill(MxPPSBeforeCuts_PU,mclumiweight);
		hMjjBC_PU->Fill(MJJBC_PU,mclumiweight);
		hMxPPSvsMjjBC_PU->Fill(MxPPS_PU,MJJ,mclumiweight);   
		//

		//PPS sel
		if(cutXArmF_PU && cutYArmF_PU && cutXArmB_PU && cutYArmB_PU){
		++counter_fiducial_PU;
		if(verbose)cout<<"x: "<<xPPSArmFDet1Info->at(k)<<", "<<xPPSArmBDet1Info->at(k)<<", "<<xPPSArmBDet1Info->at(k)<<", "<<xPPSArmBDet2Info->at(k)<<endl;
		if(verbose)cout<<"y: "<<yPPSArmFDet1Info->at(k)<<", "<<yPPSArmBDet1Info->at(k)<<", "<<yPPSArmBDet1Info->at(k)<<", "<<yPPSArmBDet2Info->at(k)<<endl;
		if(stopTrkArmF_PU && stopTrkArmB_PU){
		++counter_hasStoped_PU;
		hPPS_xVsy_ARMPlusDt1->Fill(xPPSArmFDet1Info->at(k),yPPSArmFDet1Info->at(k),mclumiweight);
		hPPS_xVsy_ARMPlusDt2->Fill(xPPSArmFDet2Info->at(k),yPPSArmFDet2Info->at(k),mclumiweight);
		hPPS_xVsy_ARMMinusDt1->Fill(xPPSArmBDet1Info->at(k),yPPSArmBDet1Info->at(k),mclumiweight);
		hPPS_xVsy_ARMMinusDt2->Fill(xPPSArmBDet2Info->at(k),yPPSArmBDet2Info->at(k),mclumiweight);

		//--------- PPS xi PU ----------------------------
		if(xiPPSArmFInfo->at(k)>0.0 && xiPPSArmBInfo->at(k)>0.0){
		++counter_xiArms_PU;
		if(JetsSameVertex_pt->at(0)>pTmin && JetsSameVertex_pt->at(1)>pTmin){
		++counter_jetpT_PU;

		if((JetsSameVertex_eta->at(0)>etamin && JetsSameVertex_eta->at(0)<etamax)&&(JetsSameVertex_eta->at(1)>etamin && JetsSameVertex_eta->at(1)<etamax)){
		++counter_jetEta_PU;
		//              hPPS_xiARMPlus->Fill(xiPPSArmF,mclumiweight);
		hPPS_xiARMPlus_PU->Fill(xiPPSArmFInfo->at(k),mclumiweight);
                hPPS_tARMPlus_PU->Fill(tPPSArmFInfo->at(k),mclumiweight);
		//hPPS_xiARMMinus->Fill(xiPPSArmB,mclumiweight);
		hPPS_xiARMMinus_PU->Fill(xiPPSArmBInfo->at(k),mclumiweight);
                hPPS_tARMMinus_PU->Fill(tPPSArmBInfo->at(k),mclumiweight);


		MxPPS_PU = S*sqrt((xiPPSArmFInfo->at(k))*(xiPPSArmBInfo->at(k)));
		if(verbose)std::cout<<"xiPPSArmFInfo (PU)="<< xiPPSArmFInfo->at(k)<<","<<"xiPPSArmBInfo (PU)="<<xiPPSArmBInfo->at(k)<<", "<<"S= "<<S<<std::endl;
		if(verbose)std::cout<<"Mx = "<<S*sqrt((xiPPSArmFInfo->at(k))*(xiPPSArmBInfo->at(k)))<<std::endl;
		//MJJ = CandidatesMjj; //REVER Dilson
		MJJ = Mjj;  // jetall
		hLeadingJetPtAfterPPS_PU->Fill(JetsSameVertex_pt->at(0),mclumiweight);
		hSecondJetPtAfterPPS_PU->Fill(JetsSameVertex_pt->at(1),mclumiweight);
		hLeadingJetEtaAfterPPS_PU->Fill(JetsSameVertex_eta->at(0),mclumiweight);
		hSecondJetEtaAfterPPS_PU->Fill(JetsSameVertex_eta->at(1),mclumiweight);
		hLeadingJetPhiAfterPPS_PU->Fill(JetsSameVertex_phi->at(0),mclumiweight);
		hSecondJetPhiAfterPPS_PU->Fill(JetsSameVertex_phi->at(1),mclumiweight);
		hDeltaEtaJetsAfterPPS_PU->Fill(deltaeta,mclumiweight);
		hDeltaPhiJetsAfterPPS->Fill(deltaphi,mclumiweight);
		//--------------------------------------------
		} // xi cut PU
		hMx_PU->Fill(MxPPS_PU,mclumiweight);
		if(verbose)std::cout<<"Mx_b = "<<S*sqrt((xiPPSArmFInfo->at(k))*(xiPPSArmBInfo->at(k)))<<std::endl;
		hMjj_PU->Fill(MJJ,mclumiweight);
		if(verbose)std::cout<<"Mjj = "<<Mjj<<" "<<"MJJ = "<<Mjj<<std::endl;//REVER Dilson
		hMxPPSvsMjj_PU->Fill(MxPPS_PU,MJJ_PU,mclumiweight);
        	} //stopTrk PU
           } //x,y cut PU
      }//eta jet cuts PU
   } //jet pt cut PU
 } //Loop PU

} //tree loop
			//----------------- print out some information ---------------
			cout<<"Events read:                      "<<NEVENTS<<endl;
		cout<<"MC cross-section:                 "<<cross_section<<" [fb]"<<endl;
		cout<<"Luminosity:                       "<<luminosity<<" [fb]^-1"<<endl;
		cout<<"Weight:                           "<<mclumiweight<<endl;
		cout<<"Events read:                      "<<NEVENTS<<endl;
		cout<<"-----------------Signal------------------------"<<endl; 
		cout<<"Events after fiducial:            "<<counter_fiducial<<endl; 
		cout<<"Events after hasStoped:           "<<counter_hasStoped<<endl; 
		cout<<"Events after XiArms>0(signal):    "<<counter_xiArms<<endl; 
		cout<<"Events after jetpT<50GeV:         "<<counter_jetpT<<endl;                                                                                          
		cout<<"Events after jetEta<2.0:          "<<counter_jetEta<<endl; 
		cout << "\n----------------------------------------------------" << endl;
		cout << "Numbers Normalized - Signal" << endl;
		cout << "------------------------------------------------------" << endl;
		cout<<"Events read:                      "<<NEVENTS*mclumiweight<<endl; 
		cout<<"Events after hasStoped:           "<<counter_hasStoped*mclumiweight<<endl; 
		cout<<"Events after XiArms>0 (signal):    "<<counter_xiArms*mclumiweight<<endl; 
		cout<<"Events after jetpT<50GeV:         "<<counter_jetpT*mclumiweight<<endl;                                                                                          
		cout<<"Events after jetEta<2.0:          "<<counter_jetEta*mclumiweight<<endl; 
		cout << "\n----------------------------------------------------" << endl;
		cout<<"-----------------PU background------------------------"<<endl;
		cout<<"Events after fiducial:            "<<counter_fiducial_PU<<endl;
		cout<<"Events after hasStoped:           "<<counter_hasStoped_PU<<endl;
		cout<<"Events after XiArms>0:            "<<counter_xiArms_PU<<endl;
		cout<<"Events after jetpT<50GeV:         "<<counter_jetpT_PU<<endl;
		cout<<"Events after jetEta<2.0:          "<<counter_jetEta_PU<<endl;
		cout << "\n----------------------------------------------------" << endl;
		cout << "Numbers Normalized -PU bck" << endl;
		cout << "------------------------------------------------------" << endl;
		cout<<"Events read:                      "<<NEVENTS*mclumiweight<<endl;
		cout<<"Events after hasStoped:           "<<counter_hasStoped_PU*mclumiweight<<endl;
		cout<<"Events after XiArms>0:            "<<counter_xiArms_PU*mclumiweight<<endl;
		cout<<"Events after jetpT<50GeV:         "<<counter_jetpT_PU*mclumiweight<<endl;                                                                                         
		cout<<"Events after jetEta<2.0:          "<<counter_jetEta_PU*mclumiweight<<endl;


		//----------------- save the histos to the output file ------
		outf->Write();
		}

