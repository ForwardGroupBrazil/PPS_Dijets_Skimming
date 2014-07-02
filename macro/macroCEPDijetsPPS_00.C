#include <stdio.h> /* printf */
#include <math.h> /* sqrt */
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
#include <Math/GenVector/Cartesian3D.h>
#include <utility>
#include <Math/GenVector/PxPyPzE4D.h>
////STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
#define PI 3.141592653589793
using namespace std;
Bool_t switchlumiweight = true;
void macroCEPDijetsPPS_00()
{
  //------------ FWLite libraries ------
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  gROOT->ProcessLine("#include<vector>");
  gROOT->ProcessLine(".exception");
  gStyle->SetOptStat(0); 
  //------------ files -----------------

  TFile *inf = TFile::Open("test_ttree_pomwig_NoOOT_PU.root");	TFile *outf = new TFile("test_histos_PPS_pomwig.root","RECREATE");

  TTree *tr = (TTree*)inf->Get("demo/Event");

  // Fixed size dimensions of array or collections stored in the TTree if any.
  //const Int_t kMaxJetVertex = 122;
  //const Int_t kMaxProtonsP4 = 2;
  //----------- define the tree branch --------

  std::vector<double>  *GenJetsPt=0;
  std::vector<double>  *GenJetsEta=0;
  std::vector<double>  *GenJetsPhi=0;
  std::vector<double>  *JetsPt=0;
  std::vector<double>  *JetsEta=0;
  std::vector<double>  *JetsPhi=0;
  std::vector<double>  *JetsSameVertex_pt=0;
  std::vector<double>  *JetsSameVertex_eta=0;
  std::vector<double>  *JetsSameVertex_phi=0;
  Double_t        JetsSameVertex_x;
  Double_t        JetsSameVertex_y;
  Double_t        JetsSameVertex_z;
  std::vector<double>  *CandidatesJets_pt=0;
  std::vector<double>  *CandidatesJets_eta=0;
  std::vector<double>  *CandidatesJets_phi=0;
  std::vector<int>     *TracksPerJet=0;
  std::vector<double>  *VertexCMSVector_x=0;
  std::vector<double>  *VertexCMSVector_y=0;
  std::vector<double>  *VertexCMSVector_z=0;
  std::vector<double>  *VertexGENVector_x=0;
  std::vector<double>  *VertexGENVector_y=0;
  std::vector<double>  *VertexGENVector_z=0;
  std::vector<double>*AllDiffVertexZVector=0;                                                       
  std::vector<double>  *DistanceBetweenJets=0;
  Double_t        MinDistanceZVertex;
  Double_t        MaxDistanceZVertex;
  Int_t           nVertex;
  Int_t           nTracks;
  Double_t        MinDistance;
  Double_t        MaxDistance;
  Double_t        GoldenVertexZ;
  Double_t        VertexZPPS;
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
  Double_t        xPPSArmFToF;
  Double_t        yPPSArmFToF;
  Double_t        xPPSArmBToF;
  Double_t        yPPSArmBToF;
  Int_t           stopPPSArmFTrkDet1;
  Int_t           stopPPSArmFTrkDet2;
  Int_t           stopPPSArmBTrkDet1;
  Int_t           stopPPSArmBTrkDet2;
  Int_t           stopPPSArmFToF;
  Int_t           stopPPSArmBToF;
  Double_t        Mjj;
  Double_t        GenMjj;
  Double_t        xigen_plus;
  Double_t        xigen_minus;
  Double_t        tgen_plus;
  Double_t        tgen_minus;
  Double_t        Mpf;
  Double_t        Rjj;                                                                                                                     
  Double_t        GenRjj;
  Int_t           eventNumber;
  Int_t           runNumber;
  Double_t        CandidatesMjj;
  Double_t        Mx;
  Double_t        GenMxx;
  Double_t        GeneratorWeight;
  Double_t        Pthat;
  Double_t        deltaToF_00;
  Double_t        deltaToF_0i;
  Double_t        deltaToF_ii;
  Double_t        tof30ps_deltaToF_00;
  Double_t        tof30ps_deltaToF_0i;
  Double_t        tof30ps_deltaToF_ii;
  Double_t        ToF_Mx_00;
  Double_t        ToF_Mx_0i;
  Double_t        ToF_Mx_ii;
  std::vector<double>  *xiPPSArmFInfo=0;
  std::vector<double>  *xiPPSArmBInfo=0;
  std::vector<double>  *tPPSArmFInfo=0;
  std::vector<double>  *tPPSArmBInfo=0;
  ////////////////////////////////////////
  std::vector<double> *xPPSArmBToFInfo_00=0;
  std::vector<double> *yPPSArmBToFInfo_00=0;

  std::vector<double> *xPPSArmFToFInfo_00=0;
  std::vector<double> *yPPSArmFToFInfo_00=0;

  std::vector<double> *xPPSArmBDet2Info_00=0;
  std::vector<double> *yPPSArmBDet2Info_00=0;

  std::vector<double> *xPPSArmFDet2Info_00=0;
  std::vector<double> *yPPSArmFDet2Info_00=0;

  std::vector<double> *xPPSArmBDet1Info_00=0;
  std::vector<double> *yPPSArmBDet1Info_00=0;

  std::vector<double> *xPPSArmFDet1Info_00=0;
  std::vector<double> *yPPSArmFDet1Info_00=0;

  std::vector<double> *tPPSArmFInfo_00=0;
  std::vector<double> *tPPSArmBInfo_00=0;

  std::vector<double> *xiPPSArmFInfo_00=0;
  std::vector<double> *xiPPSArmBInfo_00=0;

  ////////////////////////////////////////
  std::vector<double>  *xPPSArmFDet1Info=0;
  std::vector<double>  *yPPSArmFDet1Info=0;
  std::vector<double>  *xPPSArmFDet2Info=0;
  std::vector<double>  *yPPSArmFDet2Info=0;
  std::vector<double>  *xPPSArmBDet1Info=0;
  std::vector<double>  *yPPSArmBDet1Info=0;
  std::vector<double>  *xPPSArmBDet2Info=0;
  std::vector<double>  *yPPSArmBDet2Info=0;
  std::vector<double>  *xPPSArmFToFInfo=0;
  std::vector<double>  *yPPSArmFToFInfo=0;
  std::vector<double>  *xPPSArmBToFInfo=0;
  std::vector<double>  *yPPSArmBToFInfo=0;
  std::vector<double>  *stopPPSArmFToFInfo=0;
  std::vector<double>  *stopPPSArmBToFInfo=0;
  std::vector<double>  *stopPPSArmFTrkDet1Info=0;
  std::vector<double>  *stopPPSArmFTrkDet2Info=0;
  std::vector<double>  *stopPPSArmBTrkDet1Info=0;
  std::vector<double>  *stopPPSArmBTrkDet2Info=0;
  std::vector<double>  *VertexZPPSToF_00=0;
  std::vector<double>  *VertexZPPSToF_0i=0;
  std::vector<double>  *VertexZPPSToF_ii=0;
  std::vector<double>  *VertexZPPSToF_00_0i=0;
  std::vector<double>  *DijetsVertexZPPSToF=0;
  TBranch        *b_GenJetsPt=0;   
  TBranch        *b_GenJetsEta=0;   
  TBranch        *b_GenJetsPhi=0;   	
  TBranch        *b_JetsPt=0;   
  TBranch        *b_JetsEta=0;   
  TBranch        *b_JetsPhi=0;   
  TBranch        *b_JetVertex_;   
  TBranch        *b_JetVertex_fCoordinates_fX;   
  TBranch        *b_JetVertex_fCoordinates_fY;   
  TBranch        *b_JetVertex_fCoordinates_fZ;   
  TBranch        *b_JetsSameVertex_pt=0;   
  TBranch        *b_JetsSameVertex_eta=0;   
  TBranch        *b_JetsSameVertex_phi=0;   
  TBranch        *b_JetsSameVertex_x;   
  TBranch        *b_JetsSameVertex_y;   
  TBranch        *b_JetsSameVertex_z;   
  TBranch        *b_CandidatesJets_pt=0;   
  TBranch        *b_CandidatesJets_eta=0;   
  TBranch        *b_CandidatesJets_phi=0;   
  TBranch        *b_TracksPerJet=0;   
  TBranch        *b_VertexCMSVector_x=0;   
  TBranch        *b_VertexCMSVector_y=0;   
  TBranch        *b_VertexCMSVector_z=0;   
  TBranch        *b_VertexGENVector_x=0;   
  TBranch        *b_VertexGENVector_y=0;   
  TBranch        *b_VertexGENVector_z=0;   
  TBranch        *b_AllDiffVertexZVector=0;   
  TBranch        *b_DistanceBetweenJets=0;   
  TBranch        *b_MinDistanceZVertex;   
  TBranch        *b_MaxDistanceZVertex;   
  TBranch        *b_nVertex;   
  TBranch        *b_nTracks;   
  TBranch        *b_MinDistance;   
  TBranch        *b_MaxDistance;   
  TBranch        *b_GoldenVertexZ;   
  TBranch        *b_VertexZPPS;   
  TBranch        *b_xiPPSArmB;   
  TBranch        *b_xiPPSArmF;   
  TBranch        *b_tPPSArmB;   
  TBranch        *b_tPPSArmF;   
  TBranch        *b_xPPSArmBDet1;   
  TBranch        *b_xPPSArmFDet1;   
  TBranch        *b_yPPSArmBDet1;   
  TBranch        *b_yPPSArmFDet1;   
  TBranch        *b_xPPSArmBDet2;   
  TBranch        *b_xPPSArmFDet2;   
  TBranch        *b_yPPSArmBDet2;   
  TBranch        *b_yPPSArmFDet2;   
  TBranch        *b_xPPSArmFToF;   
  TBranch        *b_yPPSArmFToF;   
  TBranch        *b_xPPSArmBToF;   
  TBranch        *b_yPPSArmBToF;   
  TBranch        *b_stopPPSArmFTrkDet1;   
  TBranch        *b_stopPPSArmFTrkDet2;   
  TBranch        *b_stopPPSArmBTrkDet1;   
  TBranch        *b_stopPPSArmBTrkDet2;   
  TBranch        *b_stopPPSArmFToF;   
  TBranch        *b_stopPPSArmBToF;   
  TBranch        *b_Mjj;   
  TBranch        *b_GenMjj;   
  TBranch        *b_xigen_plus;   
  TBranch        *b_xigen_minus;   
  TBranch        *b_tgen_plus;   
  TBranch        *b_tgen_minus;   
  TBranch        *b_Mpf;   
  TBranch        *b_Rjj;   
  TBranch        *b_GenRjj;   
  TBranch        *b_CandidatesMjj;   
  TBranch        *b_Mx;   
  TBranch        *b_GenMxx;   
  TBranch        *b_GeneratorWeight;   
  TBranch        *b_Pthat;   
  TBranch        *b_eventNumber;   
  TBranch        *b_runNumber;   
  //TBranch        *b_FiducialCut;   
  TBranch        *b_deltaToF_00;   
  TBranch        *b_deltaToF_0i;   
  TBranch        *b_deltaToF_ii;   
  TBranch        *b_tof30ps_deltaToF_00;   
  TBranch        *b_tof30ps_deltaToF_0i;   
  TBranch        *b_tof30ps_deltaToF_ii;   
  TBranch        *b_ToF_Mx_00;   
  TBranch        *b_ToF_Mx_0i;   
  TBranch        *b_ToF_Mx_ii;   
  TBranch        *b_xiPPSArmFInfo=0;   
  TBranch        *b_xiPPSArmBInfo=0;   
  TBranch        *b_tPPSArmFInfo=0;   
  TBranch        *b_tPPSArmBInfo=0;   
  TBranch        *b_xPPSArmFDet1Info=0;   
  TBranch        *b_yPPSArmFDet1Info=0;   
  TBranch        *b_xPPSArmFDet2Info=0;   
  TBranch        *b_yPPSArmFDet2Info=0;   
  TBranch        *b_xPPSArmBDet1Info=0;   
  TBranch        *b_yPPSArmBDet1Info=0;   
  TBranch        *b_xPPSArmBDet2Info=0;   
  TBranch        *b_yPPSArmBDet2Info=0;   
  TBranch        *b_xPPSArmFToFInfo=0;   
  TBranch        *b_yPPSArmFToFInfo=0;   
  TBranch        *b_xPPSArmBToFInfo=0;   
  TBranch        *b_yPPSArmBToFInfo=0;   
  TBranch        *b_stopPPSArmFToFInfo=0;   
  TBranch        *b_stopPPSArmBToFInfo=0;   
  TBranch        *b_stopPPSArmFTrkDet1Info=0;   
  TBranch        *b_stopPPSArmFTrkDet2Info=0;   
  TBranch        *b_stopPPSArmBTrkDet1Info=0;   
  TBranch        *b_stopPPSArmBTrkDet2Info=0; 
  ///////////////////////////////////////////
  TBranch        *b_xiPPSArmFInfo_00=0;   
  TBranch        *b_xiPPSArmBInfo_00=0;   
  TBranch        *b_tPPSArmFInfo_00=0;   
  TBranch        *b_tPPSArmBInfo_00=0;   
  TBranch        *b_xPPSArmFDet1Info_00=0;   
  TBranch        *b_yPPSArmFDet1Info_00=0;   
  TBranch        *b_xPPSArmFDet2Info_00=0;   
  TBranch        *b_yPPSArmFDet2Info_00=0;   
  TBranch        *b_xPPSArmBDet1Info_00=0;   
  TBranch        *b_yPPSArmBDet1Info_00=0;   
  TBranch        *b_xPPSArmBDet2Info_00=0;   
  TBranch        *b_yPPSArmBDet2Info_00=0;   
  TBranch        *b_xPPSArmFToFInfo_00=0;   
  TBranch        *b_yPPSArmFToFInfo_00=0;   
  TBranch        *b_xPPSArmBToFInfo_00=0;   
  TBranch        *b_yPPSArmBToFInfo_00=0;   
  TBranch        *b_stopPPSArmFToFInfo_00=0;   
  TBranch        *b_stopPPSArmBToFInfo_00=0;   
  TBranch        *b_stopPPSArmFTrkDet1Info_00=0;   
  TBranch        *b_stopPPSArmFTrkDet2Info_00=0;   
  TBranch        *b_stopPPSArmBTrkDet1Info_00=0;   
  TBranch        *b_stopPPSArmBTrkDet2Info_00=0; 
  ///////////////////////////////////////////
  TBranch        *b_VertexZPPSToF=0;   
  TBranch        *b_VertexZPPSToF_00=0;   
  TBranch        *b_VertexZPPSToF_0i=0;   
  TBranch        *b_VertexZPPSToF_ii=0;   
  TBranch        *b_VertexZPPSToF_00_0i=0;   
  TBranch        *b_DijetsVertexZPPSToF=0;   

  //---------------------
  tr->SetBranchAddress("GenJetsPt", &GenJetsPt, &b_GenJetsPt);
  tr->SetBranchAddress("GenJetsEta", &GenJetsEta, &b_GenJetsEta);
  tr->SetBranchAddress("GenJetsPhi", &GenJetsPhi, &b_GenJetsPhi);
  tr->SetBranchAddress("JetsPt", &JetsPt, &b_JetsPt);
  tr->SetBranchAddress("JetsEta", &JetsEta, &b_JetsEta);
  tr->SetBranchAddress("JetsPhi", &JetsPhi, &b_JetsPhi);
  tr->SetBranchAddress("JetsSameVertex_pt", &JetsSameVertex_pt, &b_JetsSameVertex_pt);
  tr->SetBranchAddress("JetsSameVertex_eta", &JetsSameVertex_eta, &b_JetsSameVertex_eta);
  tr->SetBranchAddress("JetsSameVertex_phi", &JetsSameVertex_phi, &b_JetsSameVertex_phi);
  tr->SetBranchAddress("JetsSameVertex_x", &JetsSameVertex_x, &b_JetsSameVertex_x);
  tr->SetBranchAddress("JetsSameVertex_y", &JetsSameVertex_y, &b_JetsSameVertex_y);
  tr->SetBranchAddress("JetsSameVertex_z", &JetsSameVertex_z, &b_JetsSameVertex_z);
  tr->SetBranchAddress("CandidatesJets_pt", &CandidatesJets_pt, &b_CandidatesJets_pt);
  tr->SetBranchAddress("CandidatesJets_eta", &CandidatesJets_eta, &b_CandidatesJets_eta);
  tr->SetBranchAddress("CandidatesJets_phi", &CandidatesJets_phi, &b_CandidatesJets_phi);
  tr->SetBranchAddress("TracksPerJet", &TracksPerJet, &b_TracksPerJet);
  tr->SetBranchAddress("VertexCMSVector_x", &VertexCMSVector_x, &b_VertexCMSVector_x);
  tr->SetBranchAddress("VertexCMSVector_y", &VertexCMSVector_y, &b_VertexCMSVector_y);
  tr->SetBranchAddress("VertexCMSVector_z", &VertexCMSVector_z, &b_VertexCMSVector_z);
  tr->SetBranchAddress("VertexGENVector_x", &VertexGENVector_x, &b_VertexGENVector_x);
  tr->SetBranchAddress("VertexGENVector_y", &VertexGENVector_y, &b_VertexGENVector_y);
  tr->SetBranchAddress("VertexGENVector_z", &VertexGENVector_z, &b_VertexGENVector_z);
  tr->SetBranchAddress("AllDiffVertexZVector", &AllDiffVertexZVector, &b_AllDiffVertexZVector);
  tr->SetBranchAddress("DistanceBetweenJets", &DistanceBetweenJets, &b_DistanceBetweenJets);
  tr->SetBranchAddress("MinDistanceZVertex", &MinDistanceZVertex, &b_MinDistanceZVertex);
  tr->SetBranchAddress("MaxDistanceZVertex", &MaxDistanceZVertex, &b_MaxDistanceZVertex);
  tr->SetBranchAddress("nVertex", &nVertex, &b_nVertex);
  tr->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
  tr->SetBranchAddress("MinDistance", &MinDistance, &b_MinDistance);
  tr->SetBranchAddress("MaxDistance", &MaxDistance, &b_MaxDistance);
  tr->SetBranchAddress("GoldenVertexZ", &GoldenVertexZ, &b_GoldenVertexZ);
  tr->SetBranchAddress("VertexZPPS", &VertexZPPS, &b_VertexZPPS);
  tr->SetBranchAddress("xiPPSArmB", &xiPPSArmB, &b_xiPPSArmB);
  tr->SetBranchAddress("xiPPSArmF", &xiPPSArmF, &b_xiPPSArmF);
  tr->SetBranchAddress("tPPSArmB", &tPPSArmB, &b_tPPSArmB);
  tr->SetBranchAddress("tPPSArmF", &tPPSArmF, &b_tPPSArmF);
  tr->SetBranchAddress("xPPSArmBDet1", &xPPSArmBDet1, &b_xPPSArmBDet1);
  tr->SetBranchAddress("xPPSArmFDet1", &xPPSArmFDet1, &b_xPPSArmFDet1);
  tr->SetBranchAddress("yPPSArmBDet1", &yPPSArmBDet1, &b_yPPSArmBDet1);
  tr->SetBranchAddress("yPPSArmFDet1", &yPPSArmFDet1, &b_yPPSArmFDet1);
  tr->SetBranchAddress("xPPSArmBDet2", &xPPSArmBDet2, &b_xPPSArmBDet2);
  tr->SetBranchAddress("xPPSArmFDet2", &xPPSArmFDet2, &b_xPPSArmFDet2);
  tr->SetBranchAddress("yPPSArmBDet2", &yPPSArmBDet2, &b_yPPSArmBDet2);
  tr->SetBranchAddress("yPPSArmFDet2", &yPPSArmFDet2, &b_yPPSArmFDet2);
  tr->SetBranchAddress("xPPSArmFToF", &xPPSArmFToF, &b_xPPSArmFToF);
  tr->SetBranchAddress("yPPSArmFToF", &yPPSArmFToF, &b_yPPSArmFToF);
  tr->SetBranchAddress("xPPSArmBToF", &xPPSArmBToF, &b_xPPSArmBToF);
  tr->SetBranchAddress("yPPSArmBToF", &yPPSArmBToF, &b_yPPSArmBToF);
  tr->SetBranchAddress("stopPPSArmFTrkDet1", &stopPPSArmFTrkDet1, &b_stopPPSArmFTrkDet1);
  tr->SetBranchAddress("stopPPSArmFTrkDet2", &stopPPSArmFTrkDet2, &b_stopPPSArmFTrkDet2);
  tr->SetBranchAddress("stopPPSArmBTrkDet1", &stopPPSArmBTrkDet1, &b_stopPPSArmBTrkDet1);
  tr->SetBranchAddress("stopPPSArmBTrkDet2", &stopPPSArmBTrkDet2, &b_stopPPSArmBTrkDet2);
  tr->SetBranchAddress("stopPPSArmFToF", &stopPPSArmFToF, &b_stopPPSArmFToF);
  tr->SetBranchAddress("stopPPSArmBToF", &stopPPSArmBToF, &b_stopPPSArmBToF);
  tr->SetBranchAddress("Mjj", &Mjj, &b_Mjj);
  tr->SetBranchAddress("GenMjj", &GenMjj, &b_GenMjj);
  tr->SetBranchAddress("xigen_plus", &xigen_plus, &b_xigen_plus);
  tr->SetBranchAddress("xigen_minus", &xigen_minus, &b_xigen_minus);
  tr->SetBranchAddress("tgen_plus", &tgen_plus, &b_tgen_plus);
  tr->SetBranchAddress("tgen_minus", &tgen_minus, &b_tgen_minus);
  tr->SetBranchAddress("Mpf", &Mpf, &b_Mpf);
  tr->SetBranchAddress("Rjj", &Rjj, &b_Rjj);
  tr->SetBranchAddress("GenRjj", &GenRjj, &b_GenRjj);
  tr->SetBranchAddress("CandidatesMjj", &CandidatesMjj, &b_CandidatesMjj);
  tr->SetBranchAddress("Mx", &Mx, &b_Mx);
  tr->SetBranchAddress("GenMxx", &GenMxx, &b_GenMxx);
  tr->SetBranchAddress("GeneratorWeight", &GeneratorWeight, &b_GeneratorWeight);
  tr->SetBranchAddress("Pthat", &Pthat, &b_Pthat);
  tr->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
  tr->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
  //tr->SetBranchAddress("FiducialCut", &FiducialCut, &b_FiducialCut);
  tr->SetBranchAddress("deltaToF_00", &deltaToF_00, &b_deltaToF_00);
  tr->SetBranchAddress("deltaToF_0i", &deltaToF_0i, &b_deltaToF_0i);
  tr->SetBranchAddress("deltaToF_ii", &deltaToF_ii, &b_deltaToF_ii);
  tr->SetBranchAddress("tof30ps_deltaToF_00", &tof30ps_deltaToF_00, &b_tof30ps_deltaToF_00);
  tr->SetBranchAddress("tof30ps_deltaToF_0i", &tof30ps_deltaToF_0i, &b_tof30ps_deltaToF_0i);
  tr->SetBranchAddress("tof30ps_deltaToF_ii", &tof30ps_deltaToF_ii, &b_tof30ps_deltaToF_ii);
  tr->SetBranchAddress("ToF_Mx_00", &ToF_Mx_00, &b_ToF_Mx_00);
  tr->SetBranchAddress("ToF_Mx_0i", &ToF_Mx_0i, &b_ToF_Mx_0i);
  tr->SetBranchAddress("ToF_Mx_ii", &ToF_Mx_ii, &b_ToF_Mx_ii);
  tr->SetBranchAddress("xiPPSArmFInfo", &xiPPSArmFInfo, &b_xiPPSArmFInfo);
  tr->SetBranchAddress("xiPPSArmBInfo", &xiPPSArmBInfo, &b_xiPPSArmBInfo);
  tr->SetBranchAddress("tPPSArmFInfo", &tPPSArmFInfo, &b_tPPSArmFInfo);
  tr->SetBranchAddress("tPPSArmBInfo", &tPPSArmBInfo, &b_tPPSArmBInfo);
  tr->SetBranchAddress("xPPSArmFDet1Info", &xPPSArmFDet1Info, &b_xPPSArmFDet1Info);
  tr->SetBranchAddress("yPPSArmFDet1Info", &yPPSArmFDet1Info, &b_yPPSArmFDet1Info);
  tr->SetBranchAddress("xPPSArmFDet2Info", &xPPSArmFDet2Info, &b_xPPSArmFDet2Info);
  tr->SetBranchAddress("yPPSArmFDet2Info", &yPPSArmFDet2Info, &b_yPPSArmFDet2Info);
  tr->SetBranchAddress("xPPSArmBDet1Info", &xPPSArmBDet1Info, &b_xPPSArmBDet1Info);
  tr->SetBranchAddress("yPPSArmBDet1Info", &yPPSArmBDet1Info, &b_yPPSArmBDet1Info);
  tr->SetBranchAddress("xPPSArmBDet2Info", &xPPSArmBDet2Info, &b_xPPSArmBDet2Info);
  tr->SetBranchAddress("yPPSArmBDet2Info", &yPPSArmBDet2Info, &b_yPPSArmBDet2Info);
  tr->SetBranchAddress("xPPSArmFToFInfo", &xPPSArmFToFInfo, &b_xPPSArmFToFInfo);
  tr->SetBranchAddress("yPPSArmFToFInfo", &yPPSArmFToFInfo, &b_yPPSArmFToFInfo);
  tr->SetBranchAddress("xPPSArmBToFInfo", &xPPSArmBToFInfo, &b_xPPSArmBToFInfo);
  tr->SetBranchAddress("yPPSArmBToFInfo", &yPPSArmBToFInfo, &b_yPPSArmBToFInfo);
  tr->SetBranchAddress("stopPPSArmFToFInfo", &stopPPSArmFToFInfo, &b_stopPPSArmFToFInfo);
  tr->SetBranchAddress("stopPPSArmBToFInfo", &stopPPSArmBToFInfo, &b_stopPPSArmBToFInfo);
  tr->SetBranchAddress("stopPPSArmFTrkDet1Info", &stopPPSArmFTrkDet1Info, &b_stopPPSArmFTrkDet1Info);
  tr->SetBranchAddress("stopPPSArmFTrkDet2Info", &stopPPSArmFTrkDet2Info, &b_stopPPSArmFTrkDet2Info);
  tr->SetBranchAddress("stopPPSArmBTrkDet1Info", &stopPPSArmBTrkDet1Info, &b_stopPPSArmBTrkDet1Info);
  tr->SetBranchAddress("stopPPSArmBTrkDet2Info", &stopPPSArmBTrkDet2Info, &b_stopPPSArmBTrkDet2Info);
  ////////////////////
  tr->SetBranchAddress("xiPPSArmFInfo_00", &xiPPSArmFInfo_00, &b_xiPPSArmFInfo_00);
  tr->SetBranchAddress("xiPPSArmBInfo_00", &xiPPSArmBInfo_00, &b_xiPPSArmBInfo_00);
  tr->SetBranchAddress("tPPSArmFInfo_00", &tPPSArmFInfo_00, &b_tPPSArmFInfo_00);
  tr->SetBranchAddress("tPPSArmBInfo_00", &tPPSArmBInfo_00, &b_tPPSArmBInfo_00);
  tr->SetBranchAddress("xPPSArmFDet1Info_00", &xPPSArmFDet1Info_00, &b_xPPSArmFDet1Info_00);
  tr->SetBranchAddress("yPPSArmFDet1Info_00", &yPPSArmFDet1Info_00, &b_yPPSArmFDet1Info_00);
  tr->SetBranchAddress("xPPSArmFDet2Info_00", &xPPSArmFDet2Info_00, &b_xPPSArmFDet2Info_00);
  tr->SetBranchAddress("yPPSArmFDet2Info_00", &yPPSArmFDet2Info_00, &b_yPPSArmFDet2Info_00);
  tr->SetBranchAddress("xPPSArmBDet1Info_00", &xPPSArmBDet1Info_00, &b_xPPSArmBDet1Info_00);
  tr->SetBranchAddress("yPPSArmBDet1Info_00", &yPPSArmBDet1Info_00, &b_yPPSArmBDet1Info_00);
  tr->SetBranchAddress("xPPSArmBDet2Info_00", &xPPSArmBDet2Info_00, &b_xPPSArmBDet2Info_00);
  tr->SetBranchAddress("yPPSArmBDet2Info_00", &yPPSArmBDet2Info_00, &b_yPPSArmBDet2Info_00);
  tr->SetBranchAddress("xPPSArmFToFInfo_00", &xPPSArmFToFInfo_00, &b_xPPSArmFToFInfo_00);
  tr->SetBranchAddress("yPPSArmFToFInfo_00", &yPPSArmFToFInfo_00, &b_yPPSArmFToFInfo_00);
  tr->SetBranchAddress("xPPSArmBToFInfo_00", &xPPSArmBToFInfo_00, &b_xPPSArmBToFInfo_00);
  tr->SetBranchAddress("yPPSArmBToFInfo_00", &yPPSArmBToFInfo_00, &b_yPPSArmBToFInfo_00);
  //tr->SetBranchAddress("stopPPSArmFToFInfo_00", &stopPPSArmFToFInfo_00, &b_stopPPSArmFToFInfo_00);
  //tr->SetBranchAddress("stopPPSArmBToFInfo_00", &stopPPSArmBToFInfo_00, &b_stopPPSArmBToFInfo_00);
  //tr->SetBranchAddress("stopPPSArmFTrkDet1Info_00", &stopPPSArmFTrkDet1Info_00, &b_stopPPSArmFTrkDet1Info_00);
  //tr->SetBranchAddress("stopPPSArmFTrkDet2Info_00", &stopPPSArmFTrkDet2Info_00, &b_stopPPSArmFTrkDet2Info_00);
  //tr->SetBranchAddress("stopPPSArmBTrkDet1Info_00", &stopPPSArmBTrkDet1Info_00, &b_stopPPSArmBTrkDet1Info_00);
  //tr->SetBranchAddress("stopPPSArmBTrkDet2Info_00", &stopPPSArmBTrkDet2Info_00, &b_stopPPSArmBTrkDet2Info_00);
  ////////////////////
  tr->SetBranchAddress("VertexZPPSToF_00", &VertexZPPSToF_00, &b_VertexZPPSToF_00);
  tr->SetBranchAddress("VertexZPPSToF_0i", &VertexZPPSToF_0i, &b_VertexZPPSToF_0i);
  tr->SetBranchAddress("VertexZPPSToF_ii", &VertexZPPSToF_ii, &b_VertexZPPSToF_ii);
  tr->SetBranchAddress("VertexZPPSToF_00_0i", &VertexZPPSToF_00_0i, &b_VertexZPPSToF_00_0i);
  tr->SetBranchAddress("DijetsVertexZPPSToF", &DijetsVertexZPPSToF, &b_DijetsVertexZPPSToF);

  //----------- settings ---------------
  bool verbose = false;
  int NEVENTS = tr->GetEntries();
  cout<<"MC EVENTS= "<<NEVENTS<<endl;
  float cross_section =391.0 ; //[fb] exhume
  //float cross_section =291650.0 ; //[fb] pomwig
  float luminosity=1.0;//[fb]^-1
  double lumiweight = (luminosity*cross_section)/NEVENTS;
  cout<<"lumiweight= "<<lumiweight<<endl;
  double mclumiweight = 1.0;
  if (switchlumiweight ){
    mclumiweight = lumiweight;
  }

  int Njets=0;
  int NjetsGen=0;
  double jet1Genpt=-999.0;
  double jet2Genpt=-999.0;
  double jet1Geneta=-999.0;
  double jet2Geneta=-999.0;  
  double jet1Genphi=-999.0;
  double jet2Genphi=-999.0; 
  double pTmin = 100.0;
  double etamin = -2.0;
  double etamax = 2.0;
  double MxPPS = -999.0;
  double Mx00PPS = -999.0;
  double Mx0iPPS = -999.0;
  double MxiiPPS = -999.0;  
  double MJJ = -999.0;
  double diffMJJMx = -99999.0;
  double MxPPSBeforeCuts = 0.0;
  double MJJBC = -999.0;
  double deltaeta = 0.0;
  double deltaphi = 0.0;
  double S = 13000.0;
  double xmax = -3.15;
  double xmin =-23.15;
  double xmax2 = -2.03;
  double xmin2 = -22.03;
  double ymax = 9.0;
  double ymin =-9.0;
  double xipgen = -999.0;
  double ximgen = -999.0;
  double tpgen = -999.0;
  double tmgen = -999.0; 
  //--------- book histos ---------------------------------

  TH1F *hNJets = new TH1F("NJets","N Jets; N Jets; Events / 100 fb^{-1}",100,0,100);
  TH1F *hVertexZCMS = new TH1F("VertexZCMS"," Vertex Z CMS[cm]; Vertex Z [cm]; Events / 100 fb^{-1}",25,-25.0,25.0);
  TH1F *hVertexZPPS = new TH1F("VertexZPPS","Vertex Z PPS [cm]; Vertex Z [cm]; Events / 100 fb^{-1}",25,-25.0,25.0);
  TH2F *hVertexZCMSPPS = new TH2F("VertexZCMSPPS","Vertex Z CMS vs Vertex Z PPS; Vertex Z CMS [cm]; Vertex Z PPS [cm]",25,-25.0,25.0,25, -25.0,25.0);
  TH1F *hLeadingJetPt = new TH1F("LeadingJetPt","Leading Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; Events / 100 fb^{-1}",100,0,500);
  TH1F *hSecondJetPt = new TH1F("SecondJetPt","Second Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; Events / 100 fb^{-1}",100,0,500);
  TH1F *hLeadingJetEta = new TH1F("LeadingJetEta","Leading Jet - #eta Distribution; #eta; Events / 100 fb^{-1}",100,-5.2,5.2);
  TH1F *hSecondJetEta = new TH1F("SecondJetEta","Second Jet - #eta Distribution; #eta; Events / 100 fb^{-1}",100,-5.2,5.2);
  TH1F *hLeadingJetPhi = new TH1F("LeadingJetPhi","Leading Jet - #phi Distribution; #phi; Events / 100 fb^{-1}",100,-1.2*PI,1.2*PI);
  TH1F *hSecondJetPhi = new TH1F("SecondJetPhi","Second Jet - #phi Distribution; #phi; Events / 100 fb^{-1}",100,-1.2*PI,1.2*PI);
  TH1F *hDeltaEtaJets = new TH1F("DeltaEtaJets","#Delta#eta_{jj} Distribution; #Delta#eta_{jj}; Events / 100 fb^{-1}",20,0.0,5.2);
  TH1F *hDeltaPhiJets = new TH1F("DeltaPhiJets","#Delta#phi_{jj} Distribution; #Delta#phi_{jj}; Events / 100 fb^{-1}",20,0.0,1.2*PI);
  TH1F *hMjj = new TH1F( "Mjj" , "Mass_{JJ} Distribution; M_{jj} [GeV]; Events / 100 fb^{-1}" , 100, 0., 2000. );
  TH1F *hMx = new TH1F("Mx" , "Mass_{X} Distribution; M_{x} [GeV]; Events / 100 fb^{-1}" , 100, 0., 2000. );
  TH1F *hMjjBC = new TH1F( "MjjBC" , "Mass_{JJ} Distribution; M_{jj} [GeV]; Events / 100 fb^{-1}" , 100, 0., 2000. );
  TH1F *hMxBC = new TH1F("MxBC" , "Mass_{X} Distribution; M_{x} [GeV]; Events / 100 fb^{-1}" , 100, 0., 2000. );
  TH1F *hPPS_xiARMPlus = new TH1F( "PPS_xiARMPlus" , "#xi_{plus} PPS; #xi_{plus}; Events / 100 fb^{-1}" , 100, 0., 1.0 );
  TH1F *hPPS_xiARMMinus = new TH1F( "PPS_xiARMMinus" , "#xi_{minus} PPS; #xi_{minus}; Events / 100 fb^{-1}" , 100, 0., 1.0 );
  TH1F *hPPS_tARMPlus = new TH1F( "PPS_tARMPlus" , "|t|_{plus} PPS; |t|_{plus} [GeV^{2}]; Events / 100 fb^{-1}" , 100, 0., 4.0 );
  TH1F *hPPS_tARMMinus = new TH1F( "PPS_tARMMinus" , "|t|_{minus} PPS; |t|_{minus} [GeV^{2}]; Events / 100 fb^{-1}" , 100, 0., 4.0 );
  TH2F *hPPS_xVsy_ARMPlusDt1 = new TH2F( "PPS_xVsy_ARMPlusDt1" , "PPS x_vs_y_{ARMPlusDt1}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 );
  TH2F *hPPS_xVsy_ARMPlusDt2 = new TH2F( "PPS_xVsy_ARMPlusDt2" , "PPS x_vs_y_{ARMPlusDt2}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 );
  TH2F *hPPS_xVsy_ARMMinusDt1 = new TH2F( "PPS_xVsy_ARMMinusDt1" , "PPS x_vs_y_{ARMMinusDt1}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 );
  TH2F *hPPS_xVsy_ARMMinusDt2 = new TH2F( "PPS_xVsy_ARMMinusDt2" , "PPS x_vs_y_{ARMMinusDt2}; x [mm]; y [mm]" , 100, -50.0, 10.0,100,-10.0,10.0 );
  TH2F *hMxPPSvsMjj = new TH2F("MxvsMjj" , "Mass_{X} vs M_{JJ} Distribution; M_{x} [GeV]; M_{jj} [GeV]" , 100, 0., 2000.,100, 0., 2000. );
  TH2F *hMx00PPSvsMjj = new TH2F("Mx00vsMjj" , "Mass00_{X} vs M_{JJ} Distribution; M00_{x} [GeV]; M_{jj} [GeV]" , 100, 0., 2000.,100, 0., 2000. ); 
  //TH2F *hMx0iPPSvsMjj = new TH2F("Mx0ivsMjj" , "Mass0i_{X} vs M_{JJ} Distribution; M0i_{x} [GeV]; M_{jj} [GeV]" , 100, 0., 2000.,100, 0.,  2000. ); 
  //TH2F *hMxiiPPSvsMjj = new TH2F("MxiivsMjj" , "Massii_{X} vs M_{JJ} Distribution; Mii_{x} [GeV]; M_{jj} [GeV]" , 100, 0., 2000.,100, 0.,  2000. ); 
  TH2F *hMxPPSvsMjjBC = new TH2F("MxvsMjjBC" , "Mass_{X} vs M_{JJ} Distribution; M_{x} [GeV]; M_{jj} [GeV]" , 100, 0., 2000.,100, 0., 2000. );
  TH1F *hLeadingJetPtAfterPPS = new TH1F("LeadingJetPtAfterPPS","Leading Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; Events / 100 fb^{-1}",100,0,500);
  TH1F *hSecondJetPtAfterPPS = new TH1F("SecondJetPtAfterPPS","Second Jet - P_{T} Distribution; P_{T} [GeV.c^{-1}]; Events / 100 fb^{-1}",100,0,500);
  TH1F *hLeadingJetEtaAfterPPS = new TH1F("LeadingJetEtaAfterPPS","Leading Jet - #eta Distribution; #eta; Events / 100 fb^{-1}",100,-5.2,5.2);
  TH1F *hSecondJetEtaAfterPPS = new TH1F("SecondJetEtaAfterPPS","Second Jet - #eta Distribution; #eta; Events / 100 fb^{-1}",100,-5.2,5.2);
  TH1F *hLeadingJetPhiAfterPPS = new TH1F("LeadingJetPhiAfterPPS","Leading Jet - #phi Distribution; #phi; Events / 100 fb^{-1}",100,-1.2*PI,1.2*PI);
  TH1F *hSecondJetPhiAfterPPS = new TH1F("SecondJetPhiAfterPPS","Second Jet - #phi Distribution; #phi; Events / 100 fb^{-1}",100,-1.2*PI,1.2*PI);
  TH1F *hDeltaEtaJetsAfterPPS = new TH1F("DeltaEtaJetsAfterPPS","#Delta#eta_{jj} Distribution; #Delta#eta_{jj}; Events / 100 fb^{-1}",20,0.0,5.2);
  TH1F *hDeltaPhiJetsAfterPPS = new TH1F("DeltaPhiJetsAfterPPS","#Delta#phi_{jj} Distribution; #Delta#phi_{jj}; Events / 100 fb^{-1}",20,0.0,1.2*PI);
  TH1F *hdifMxMjj  = new TH1F("diffMjjMx" , "(Mass_{JJ}-M_{X}) Distribution; M_{jj}-M_{X} [GeV]; Events / 100 fb^{-1}" , 100, -2000.0, 2000. );
  //----------- counters -----------------------
  int counter_cand(0),counter_jetpT(0),counter_jetEta(0), counter_fiducial(0), counter_hasStoped(0),counter_xiArms(0),counter_jetpT(0), counter_vert(0), counter_gen(0),counter_jetpTGen(0),counter_jetEtaGen(0);
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
    //cout<< stopPPSArmFTrkDet1Info->size() << endl;
    Njets = JetsSameVertex_pt->size();
    hNJets->Fill(Njets);
    //cout<<"Njets"<<Njets<<endl; 
    if (JetsSameVertex_pt->size()<2)continue;
    if(verbose)cout<<"JetsSameVertex_pt->size ="<<JetsSameVertex_pt->size()<<endl;
    ++counter_cand;
    if((JetsSameVertex_eta->at(0)>etamin && JetsSameVertex_eta->at(0)<etamax)&&(JetsSameVertex_eta->at(1)>etamin && JetsSameVertex_eta->at(1)<etamax)){
      ++counter_jetEta;
      if(JetsSameVertex_pt->at(0)>pTmin && JetsSameVertex_pt->at(1)>pTmin){
	++counter_jetpT;
	cout << GoldenVertexZ << " VertexZPPSToF_00->at(0) =  " << VertexZPPSToF_00->at(0) << endl;
	if(abs(GoldenVertexZ-VertexZPPSToF_00->at(0))<0.2){
	  //cout << stopPPSArmFTrkDet1Info<< "  " << xPPSArmFDet1Info_00->size() << endl;
	  ++counter_vert; cout << counter_vert << endl;
	  //------- fill some Jets histograms ------------------
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
	  VertexZCMSPPS->Fill(GoldenVertexZ,VertexZPPSToF_00->at(0),mclumiweight);
	  hVertexZCMS->Fill(GoldenVertexZ,mclumiweight);
	  hVertexZPPS->Fill(VertexZPPSToF_00->at(0),mclumiweight);
	  //cout << VertexZPPSToF_00->size() << endl;
	  for (int k=0;k<xPPSArmFDet1Info_00->size();k++){  
	    //--------- PPS x vs y -----------------------
	    //x,det1 e det2 nos dois arms < - 3.1 e > -23.1 mm
	    //y,det1 e det2 nos dois arms > -9 e < 9 mm
	    bool cutXArmF = ((xmin<(xPPSArmFDet1Info_00->at(k)) && (xPPSArmFDet1Info_00->at(k))<xmax) && (xmin2<(xPPSArmFDet2Info_00->at(k)) && (xPPSArmFDet2Info_00->at(k))<xmax2));
	    bool cutXArmB = ((xmin<(xPPSArmBDet1Info_00->at(k)) && (xPPSArmBDet1Info_00->at(k))<xmax) && (xmin2<(xPPSArmBDet2Info_00->at(k)) && (xPPSArmBDet2Info_00->at(k))<xmax2));
	    bool cutYArmF = ((ymin<(yPPSArmFDet1Info_00->at(k)) && (yPPSArmFDet1Info_00->at(k))<ymax) && (ymin<(yPPSArmFDet2Info_00->at(k)) && (yPPSArmFDet2Info_00->at(k))<ymax));
	    bool cutYArmB = ((ymin<(yPPSArmBDet1Info_00->at(k)) && (yPPSArmBDet1Info_00->at(k))<ymax) && (ymin<(yPPSArmBDet2Info_00->at(k)) && (yPPSArmBDet2Info_00->at(k))<ymax));
	    if(verbose) std:: cout<<  "XArmF "<< cutXArmF << endl;
	    if(verbose)   std:: cout<<  "XArmB "<< cutXArmB << endl;
	    if(verbose)     std:: cout<<  "YArmF "<< cutYArmF << endl;
	    if(verbose)       std:: cout<<  "YArmB "<< cutYArmB << endl;
	      cout << xPPSArmFDet1Info_00->at(k) << " " <<  xPPSArmFDet1Info_00->at(k) << " " <<  xPPSArmFDet2Info_00->at(k) << " " <<  xPPSArmFDet2Info_00->at(k) << endl;
              cout <<  "XArmF "<< cutXArmF << " XArmB "<< cutXArmB <<  " YArmF "<< cutYArmF <<  " YArmB "<< cutYArmB << endl;

	    //      std::cout<<"XArmF:" << cutXArmF<"XArmB:"<< cutXArmB <<"YArmF: "<< cutXArmF <<"YArmB: "<< cutYArmB <std::endl;
	    //--------- PPS noHasStoped -----------------------
	    bool stopTrkArmF = (stopPPSArmFTrkDet1Info->at(0)==0 && stopPPSArmFTrkDet2Info->at(0)==0);
	    bool stopTrkArmB = (stopPPSArmBTrkDet1Info->at(0)==0 && stopPPSArmBTrkDet2Info->at(0)==0);
	    // if(verbose)cout<<"stopTrkArmF:"<<stopTrkArmFInfo->at(0)<<" , "<<"stopTrkArmB"<<stopTrkArmBInfo->at(0)<<endl;
	    // Mjj and Mx before PPS sel
	    if(xiPPSArmFInfo_00->at(k)>0.0 && xiPPSArmBInfo_00->at(k)>0.0){
	      MxPPSBeforeCuts = S*sqrt(xiPPSArmFInfo_00->at(k)*xiPPSArmBInfo_00->at(k));
	      if(verbose)cout<<"Mx BF = "<<S*sqrt(xiPPSArmFInfo_00->at(k)*xiPPSArmBInfo_00->at(k))<<endl;
	      MJJBC = Mjj; /// Eh PRA USAR CANDIDATES ?????????????
	    }
	    hMxBC->Fill(MxPPSBeforeCuts,mclumiweight);
	    hMjjBC->Fill(MJJBC,mclumiweight);
	    hMxPPSvsMjjBC->Fill(MxPPS,MJJ,mclumiweight);
	    //PPS sel
	    if(cutXArmF && cutYArmF && cutXArmB && cutYArmB){
	      ++counter_fiducial;
	      if(verbose)cout<<"x: "<<xPPSArmFDet1Info_00->at(k)<<", "<<xPPSArmBDet1Info_00->at(k)<<", "<<xPPSArmBDet1Info_00->at(k)<<", "<<xPPSArmBDet2Info_00->at(k)<<endl;
	      if(verbose)cout<<"y: "<<yPPSArmFDet1Info_00->at(k)<<", "<<yPPSArmBDet1Info_00->at(k)<<", "<<yPPSArmBDet1Info_00->at(k)<<", "<<yPPSArmBDet2Info_00->at(k)<<endl;
	      if(stopTrkArmF && stopTrkArmB){
		++counter_hasStoped;
		hPPS_xVsy_ARMPlusDt1->Fill(xPPSArmFDet1Info_00->at(k),yPPSArmFDet1Info_00->at(k),mclumiweight);
		hPPS_xVsy_ARMPlusDt2->Fill(xPPSArmFDet2Info_00->at(k),yPPSArmFDet2Info_00->at(k),mclumiweight);
		hPPS_xVsy_ARMMinusDt1->Fill(xPPSArmBDet1Info_00->at(k),yPPSArmBDet1Info_00->at(k),mclumiweight);
		hPPS_xVsy_ARMMinusDt2->Fill(xPPSArmBDet2Info_00->at(k),yPPSArmBDet2Info_00->at(k),mclumiweight);
		//--------- PPS xi ----------------------------
		if(xiPPSArmFInfo_00->at(k)>0.0 && xiPPSArmBInfo_00->at(k)>0.0){
		  ++counter_xiArms;
		  hPPS_tARMPlus->Fill(tPPSArmFInfo_00->at(k),mclumiweight);
		  hPPS_tARMMinus->Fill(tPPSArmBInfo_00->at(k),mclumiweight);
		  hPPS_xiARMPlus->Fill(xiPPSArmFInfo_00->at(k),mclumiweight);
		  hPPS_xiARMMinus->Fill(xiPPSArmBInfo_00->at(k),mclumiweight);

		  MxPPS = S*sqrt(xiPPSArmFInfo_00->at(k)*xiPPSArmBInfo_00->at(k));
		  if(verbose)std::cout<<"xiPPSArmFInfo_00->at(k)="<< xiPPSArmFInfo_00->at(k)<<","<<"xiPPSArmBInfo_00="<<xiPPSArmBInfo_00->at(k)<<", "<<"S= "<<S<<std::endl;
		  if(verbose)std::cout<<"Mx = "<<S*sqrt(xiPPSArmFInfo_00->at(k)*xiPPSArmBInfo_00->at(k))<<std::endl;
		  MJJ = Mjj;
		  diffMJJMx=(MJJ-MxPPS);
		  Mx00PPS=ToF_Mx_00;
		  if(verbose)std::cout<<"Mx00 = "<<ToF_Mx_00<<std::endl;


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
		hMjj->Fill(MJJ,mclumiweight);
		hdifMxMjj->Fill(diffMJJMx,mclumiweight);
		hMxPPSvsMjj->Fill(MxPPS,MJJ,mclumiweight);
		hMx00PPSvsMjj->Fill(Mx00PPS,MJJ,mclumiweight);

	      } //stopTrk
	    }//cut x,y 
	  } //loop
	}
      } //jet pt cut
    }//eta jet cuts
  }// tree loop
  //----------------- print out some information ---------------
  cout<<"Events read: "<<NEVENTS<<endl;
  cout<<"MC cross-section: "<<cross_section<<" [fb]"<<endl;
  cout<<"Luminosity: "<<luminosity<<" [fb]^-1"<<endl;
  cout<<"Weight: "<<mclumiweight<<endl;
  cout<<"Events read: "<<NEVENTS<<endl;
  cout<<"Events after jetEta<2.0: "<<counter_jetEta<<endl;
  cout<<"Events after jetpT>100GeV: "<<counter_jetpT<<endl;
  cout<<"Events after (cms-pps)<0.2 vertex: "<<counter_vert<<endl;
  cout<<"Events after fiducial: "<<counter_fiducial<<endl;
  cout<<"Events after hasStoped: "<<counter_hasStoped<<endl;
  cout<<"Events after XiArms>0: "<<counter_xiArms<<endl;
  cout << "\n----------------------------------------------------" << endl;
  cout << "Numbers Normalized" << endl;
  cout << "------------------------------------------------------" << endl;
  cout<<"Events read: "<<NEVENTS*mclumiweight<<endl;
  cout<<"Events after jetEta<2.0: "<<counter_jetEta*mclumiweight<<endl;
  cout<<"Events after jetpT>100GeV: "<<counter_jetpT*mclumiweight<<endl;
  cout<<"Events after (cms-pps)<0.2 vertex: "<<counter_vert*mclumiweight<<endl;
  cout<<"Events after hasStoped: "<<counter_hasStoped*mclumiweight<<endl;
  cout<<"Events after XiArms>0: "<<counter_xiArms*mclumiweight<<endl;

  //----------------- save the histos to the output file ------
  outf->Write();
}


