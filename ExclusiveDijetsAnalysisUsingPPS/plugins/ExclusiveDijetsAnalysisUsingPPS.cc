// -*- C++ -*-
//
// Package:    ExclusiveDijetsAnalysisUsingPPS
// Class:      ExclusiveDijetsAnalysisUsingPPS
//
// Authors: PPS CEP Brazilian Group
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// dataformats
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PPSObjects/interface/PPSSpectrometer.h"
#include "DataFormats/PPSObjects/interface/PPSData.h"
#include "DataFormats/PPSObjects/interface/PPSDetector.h"
#include "DataFormats/PPSObjects/interface/PPSToF.h"
#include "DataFormats/TrackReco/interface/Track.h"

// Tracks Associated with Jets
#include "DataFormats/JetReco/interface/JetTrackMatch.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"

// root
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"

// c++
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <map>

using namespace edm;
using namespace std;

//
// class declaration
//

class ExclusiveDijetsAnalysisUsingPPS : public edm::EDAnalyzer {
  public:
    explicit ExclusiveDijetsAnalysisUsingPPS(const edm::ParameterSet&);
    ~ExclusiveDijetsAnalysisUsingPPS();


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    void Init();

    // ----------member data ---------------------------

    void FillCollections(const edm::Event&, const edm::EventSetup&, bool debug);
    void SortingObjects(const edm::Event&, const edm::EventSetup&, bool debug);
    void AssociateJetsWithVertex(const edm::Event&, const edm::EventSetup&, bool debug);

    edm::InputTag jetTag_;
    edm::InputTag particleFlowTag_;
    std::string ppsTag_;
    double pTPFThresholdCharged_;
    double energyPFThresholdBar_;
    double energyPFThresholdEnd_;
    double energyPFThresholdHF_;

    int indexGold;
    std::vector<const reco::PFJet*> JetsVector;
    std::vector<const reco::Vertex*> VertexVector;
    std::vector<const reco::Track*> TracksVector;
    std::vector<const reco::PFCandidate*> PFVector;
    std::vector<const PPSSpectrometer*> PPSSpecVector;
    std::vector< std::pair<double,double> > PPSCMSVertex;
    std::vector<double> MinimumDistance;   

    std::vector<double> JetsVector_pt;
    std::vector<double> JetsVector_eta;
    std::vector<double> JetsVector_phi;
    std::vector<double> PFVector_pt;
    std::vector<double> PFVector_eta;

    TTree* eventTree_;

    int nTracks;
    int nVertex;
    double GoldenVertexZ;
    double MinDistance;
    double MaxDistance;
    double xiPPSArmF;
    double xiPPSArmB;
    double tPPSArmF;
    double tPPSArmB;
    double xPPSArmFDet1, yPPSArmFDet1;
    double xPPSArmBDet1, yPPSArmBDet1;
    double xPPSArmFDet2, yPPSArmFDet2;
    double xPPSArmBDet2, yPPSArmBDet2;
    double Mjj;
    double Mpf;
    double Rjj;

};


//
// constructors and destructor
//

ExclusiveDijetsAnalysisUsingPPS::ExclusiveDijetsAnalysisUsingPPS(const edm::ParameterSet& iConfig):
  jetTag_(iConfig.getParameter<edm::InputTag>("JetTag")),
  particleFlowTag_(iConfig.getParameter<edm::InputTag>("ParticleFlowTag")),
  ppsTag_(iConfig.getUntrackedParameter<std::string>("PPSTag","PPSReco")),
  pTPFThresholdCharged_(iConfig.getParameter<double>("pTPFThresholdCharged")),
  energyPFThresholdBar_(iConfig.getParameter<double>("energyPFThresholdBar")),
  energyPFThresholdEnd_(iConfig.getParameter<double>("energyPFThresholdEnd")),
  energyPFThresholdHF_(iConfig.getParameter<double>("energyPFThresholdHF"))
{

  edm::Service<TFileService> fs;
  eventTree_ = fs->make<TTree>("Event","Event");
  eventTree_->Branch("JetsPt",&JetsVector_pt);
  eventTree_->Branch("JetsEta",&JetsVector_eta);
  eventTree_->Branch("JetsPhi",&JetsVector_phi);
  eventTree_->Branch("PFCandidatePt",&PFVector_pt);
  eventTree_->Branch("PFCandidateEta",&PFVector_eta);
  eventTree_->Branch("nVertex",&nVertex,"nVertex/I");
  eventTree_->Branch("nTracks",&nTracks,"nTracks/I");
  eventTree_->Branch("MinDistance",&MinDistance,"MinDistance/D");
  eventTree_->Branch("MaxDistance",&MinDistance,"MaxDistance/D");
  eventTree_->Branch("GoldenVertexZ",&GoldenVertexZ,"GoldenVertexZ/D");
  eventTree_->Branch("xiPPSArmB",&xiPPSArmB,"xiPPSArmB/D");
  eventTree_->Branch("xiPPSArmF",&xiPPSArmF,"xiPPSArmF/D");
  eventTree_->Branch("tPPSArmB",&tPPSArmB,"tPPSArmB/D");
  eventTree_->Branch("tPPSArmF",&tPPSArmF,"tPPSArmF/D");
  eventTree_->Branch("xPPSArmBDet1",&xPPSArmBDet1,"xPPSArmBDet1/D");
  eventTree_->Branch("xPPSArmFDet1",&xPPSArmFDet1,"xPPSArmFDet1/D");
  eventTree_->Branch("yPPSArmBDet1",&yPPSArmBDet1,"yPPSArmBDet1/D");
  eventTree_->Branch("yPPSArmFDet1",&yPPSArmFDet1,"yPPSArmFDet1/D");
  eventTree_->Branch("xPPSArmBDet2",&xPPSArmBDet2,"xPPSArmBDet2/D");
  eventTree_->Branch("xPPSArmFDet2",&xPPSArmFDet2,"xPPSArmFDet2/D");
  eventTree_->Branch("yPPSArmBDet2",&yPPSArmBDet2,"yPPSArmBDet2/D");
  eventTree_->Branch("yPPSArmFDet2",&yPPSArmFDet2,"yPPSArmFDet2/D");
  eventTree_->Branch("Mjj",&Mjj,"Mjj/D");
  eventTree_->Branch("Mpf",&Mpf,"Mpf/D");
  eventTree_->Branch("Rjj",&Rjj,"Rjj/D");

}

ExclusiveDijetsAnalysisUsingPPS::~ExclusiveDijetsAnalysisUsingPPS()
{
}


// ------------ method called once each job just before starting event loop  ------------
void ExclusiveDijetsAnalysisUsingPPS::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void ExclusiveDijetsAnalysisUsingPPS::endJob()
{
}

// ------------ method called for each event  ------------
void ExclusiveDijetsAnalysisUsingPPS::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  Init(); // Clean Variables
  FillCollections(iEvent, iSetup, false); //true-> Print Outputs and Golden Vertex (PPS and CMS). False-> No print screen.
  SortingObjects(iEvent, iSetup, false); //true-> Print Ordered Outputs, False-> No print screen.
  AssociateJetsWithVertex(iEvent, iSetup, false); //true-> Print Ordered Outputs, False-> No print screen.
  eventTree_->Fill();

}

// ------------ Clean Variables --------------

void ExclusiveDijetsAnalysisUsingPPS::Init(){

  JetsVector.clear();
  VertexVector.clear();
  TracksVector.clear();
  PFVector.clear();
  PPSSpecVector.clear();
  PPSCMSVertex.clear();
  MinimumDistance.clear();

  JetsVector_pt.clear();
  JetsVector_eta.clear();
  JetsVector_phi.clear();
  PFVector_pt.clear();
  PFVector_eta.clear();

  nTracks = 0;
  nVertex = 0;
  GoldenVertexZ = -999.;
  MinDistance = -1.;
  MaxDistance = -1.;
  xiPPSArmF = -999.;
  xiPPSArmB = -999.;
  tPPSArmF = -999.;
  tPPSArmB = -999.;
  xPPSArmFDet1 = -999.; yPPSArmFDet1 = -999.;
  xPPSArmBDet1= -999.; yPPSArmBDet1= -999.;
  xPPSArmFDet2= -999.; yPPSArmFDet2= -999.;
  xPPSArmBDet2= -999.; yPPSArmBDet2= -999.;
  Mjj= -999.;
  Mpf= -999.;
  Rjj= -999.;

}


// ------------ Fill Vectors, All Handles  ------------
void ExclusiveDijetsAnalysisUsingPPS::FillCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug)
{

  // Debug Detailed Information, each loop
  bool debugdetails = false;

  // Fill Vertex
  Handle<edm::View<reco::Vertex> > vertex;
  iEvent.getByLabel("offlinePrimaryVertices", vertex);

  int vertexsize = vertex->size();
  int itVertex;

  if(vertex->size()>0){
    for(itVertex=0; itVertex < vertexsize; ++itVertex){
      const reco::Vertex* vertexAll = &((*vertex)[itVertex]);
      if (vertexAll->isValid()==0) continue; 
      VertexVector.push_back(vertexAll);
    }
  }

  // Fill Tracks
  Handle<edm::View<reco::Track> > tracks;
  iEvent.getByLabel("generalTracks", tracks);

  int trackssize = tracks->size();
  int itTracks;

  if(tracks->size()>0){
    for(itTracks=0; itTracks < trackssize; ++itTracks){
      const reco::Track* tracksAll = &((*tracks)[itTracks]);
      TracksVector.push_back(tracksAll);
    }
  }

  nTracks = TracksVector.size();
  nVertex = VertexVector.size();

  // Fill Jets
  Handle<edm::View<reco::PFJet> > jets;
  iEvent.getByLabel(jetTag_,jets);

  int jetsize = jets->size();
  int itJets;

  if(jets->size()>0){
    for(itJets=0; itJets < jetsize; ++itJets){
      const reco::PFJet* jetAll = &((*jets)[itJets]);
      JetsVector.push_back(jetAll);
    }
  }

  // Fill Particle Flow
  Handle <reco::PFCandidateCollection> PFCandidates;
  iEvent.getByLabel(particleFlowTag_,PFCandidates);
  reco::PFCandidateCollection::const_iterator iter;

  int pfsize = PFCandidates->size();
  int itPF;

  if(PFCandidates->size()>0){

    math::XYZTLorentzVector allCands(0.,0.,0.,0.);
    for(itPF=0; itPF < pfsize; ++itPF){
      const reco::PFCandidate* pfAll = &((*PFCandidates)[itPF]);
      double energy=pfAll->energy();
      double pt=pfAll->pt();
      double eta=pfAll->eta();
      double charge=pfAll->charge();
      if (fabs(eta)>4.7) continue;
      if ( (fabs(charge) >0 && pt > pTPFThresholdCharged_ ) ||
	  (fabs(charge) == 0 && ( (fabs(eta) <= 1.5 && energy > energyPFThresholdBar_) ||
				  (fabs(eta) > 1.5 && fabs(eta) <= 3 && energy > energyPFThresholdEnd_) ||
				  (fabs(eta) > 3 && energy >energyPFThresholdHF_) ) ) )
      { 
	allCands+=pfAll->p4();
	PFVector.push_back(pfAll);
      }
    }
    Mpf = allCands.M();
  }

  // Fill PPS Spectrometer
  Handle<PPSSpectrometer> ppsSpectrum;
  iEvent.getByLabel("ppssim",ppsTag_,ppsSpectrum);

  // Xi and t, ArmF and ArmB
  if(ppsSpectrum->ArmB.xi.size() > 0){
    xiPPSArmB = ppsSpectrum->ArmB.xi[0];
  }

  if(ppsSpectrum->ArmF.xi.size() > 0){
    xiPPSArmF = ppsSpectrum->ArmF.xi[0];
  }

  if(ppsSpectrum->ArmB.t.size() > 0){
    tPPSArmB = ppsSpectrum->ArmB.t[0];
  }

  if(ppsSpectrum->ArmF.t.size() > 0){
    tPPSArmF = ppsSpectrum->ArmF.t[0];
  }

  // ArmF and ArmB, Det1 info (x,y)
  if(ppsSpectrum->ArmF.TrkDet1.X.size() > 0){
    xPPSArmFDet1 = ppsSpectrum->ArmF.TrkDet1.X[0];
    yPPSArmFDet1 = ppsSpectrum->ArmF.TrkDet1.Y[0];
  }

  if(ppsSpectrum->ArmB.TrkDet1.X.size() > 0){
    xPPSArmBDet1 = ppsSpectrum->ArmB.TrkDet1.X[0];
    yPPSArmBDet1 = ppsSpectrum->ArmB.TrkDet1.Y[0];
  }

  // ArmF and ArmB, Det2 info (x,y)
  if(ppsSpectrum->ArmF.TrkDet2.X.size() > 0){
    xPPSArmFDet2 = ppsSpectrum->ArmF.TrkDet2.X[0];
    yPPSArmFDet2 = ppsSpectrum->ArmF.TrkDet2.Y[0];
  }

  if(ppsSpectrum->ArmB.TrkDet2.X.size() > 0){
    xPPSArmBDet2 = ppsSpectrum->ArmB.TrkDet2.X[0];
    yPPSArmBDet2 = ppsSpectrum->ArmB.TrkDet2.Y[0];
  }

  if (debug){
    cout << "\n--PPS INFO--" << endl;
    if (ppsSpectrum->vtxZ.size() > 0) cout << "vtxZ[0]: " << ppsSpectrum->vtxZ[0] << " | Vector Size: " << ppsSpectrum->vtxZ.size() << endl;
    if (ppsSpectrum->ArmF.t.size() > 0) cout << "ArmF.t[0]: " << ppsSpectrum->ArmF.t[0] << " | Vector Size: " << ppsSpectrum->ArmF.t.size() << endl;
    if (ppsSpectrum->ArmB.t.size() > 0) cout << "ArmB.t[0]: " << ppsSpectrum->ArmB.t[0] << " | Vector Size: " << ppsSpectrum->ArmB.t.size() << endl;  
    if (ppsSpectrum->ArmF.xi.size() > 0) cout << "ArmF.xi[0]: " << ppsSpectrum->ArmF.xi[0] << " | Vector Size: " << ppsSpectrum->ArmF.xi.size() << endl;
    if (ppsSpectrum->ArmB.xi.size() > 0) cout << "ArmB.xi[0]: " << ppsSpectrum->ArmB.xi[0] << " | Vector Size: " << ppsSpectrum->ArmB.xi.size() << endl;
    if (ppsSpectrum->ArmF.TrkDet1.X.size() > 0) cout << "ArmF.TrkDet1[0](x,y): (" << ppsSpectrum->ArmF.TrkDet1.X[0] << "," << ppsSpectrum->ArmF.TrkDet1.Y[0] << ") mm" << " | Vector Size: " << ppsSpectrum->ArmF.TrkDet1.X.size() <<  endl;
    if (ppsSpectrum->ArmB.TrkDet1.X.size() > 0) cout << "ArmB.TrkDet1[0](x,y): (" << ppsSpectrum->ArmB.TrkDet1.X[0] << "," << ppsSpectrum->ArmB.TrkDet1.Y[0] << ") mm" << " | Vector Size: " << ppsSpectrum->ArmB.TrkDet1.X.size() << endl;
    if (ppsSpectrum->ArmF.TrkDet2.X.size() > 0) cout << "ArmF.TrkDet2[0](x,y): (" << ppsSpectrum->ArmF.TrkDet2.X[0] << "," << ppsSpectrum->ArmF.TrkDet2.Y[0] << ") mm" << " | Vector Size: " << ppsSpectrum->ArmF.TrkDet2.X.size() << endl;
    if (ppsSpectrum->ArmB.TrkDet2.X.size() > 0) cout << "ArmB.TrkDet2[0](x,y): (" << ppsSpectrum->ArmB.TrkDet2.X[0] << "," << ppsSpectrum->ArmB.TrkDet2.Y[0] << ") mm" << " | Vector Size: " << ppsSpectrum->ArmB.TrkDet2.X.size() << endl;
    if (ppsSpectrum->ArmF.ToFDet.X.size() > 0) cout << "ArmF.ToFDet[0](x,y): (" << ppsSpectrum->ArmF.ToFDet.X[0] << "," << ppsSpectrum->ArmF.ToFDet.Y[0] << ") mm" << " | Vector Size: " << ppsSpectrum->ArmF.ToFDet.X.size() << endl;
    if (ppsSpectrum->ArmB.ToFDet.X.size() > 0) cout << "ArmB.ToFDet[0](x,y): (" << ppsSpectrum->ArmB.ToFDet.X[0] << "," << ppsSpectrum->ArmB.ToFDet.Y[0] << ") mm" << " | Vector Size: " << ppsSpectrum->ArmB.ToFDet.X.size() << endl;

    if(debugdetails){ 

      cout << "ArmF.TrkDet1: " << ppsSpectrum->ArmF.TrkDet1.X.size() << endl;
      cout << "ArmB.TrkDet1: " << ppsSpectrum->ArmB.TrkDet1.X.size() << endl;
      cout << "ArmF.TrkDet2: " << ppsSpectrum->ArmF.TrkDet2.X.size() << endl;
      cout << "ArmB.TrkDet2: " << ppsSpectrum->ArmB.TrkDet2.X.size() << endl;
      cout << "ArmF.ToFDet: " << ppsSpectrum->ArmF.ToFDet.X.size() << endl;
      cout << "ArmB.ToFDet: " << ppsSpectrum->ArmB.ToFDet.X.size() << endl;

      // Vertex PPS Info
      for (unsigned int i=0;i<ppsSpectrum->vtxZ.size();i++){
	cout << "\nPPS: vtxZ[" << i << "]: " << ppsSpectrum->vtxZ[i] << " mm" << endl;
      }

      // ARM-F
      for (unsigned int i=0;i<ppsSpectrum->ArmF.t.size();i++){
	cout << "ArmF: t[" << i << "]: " << ppsSpectrum->ArmF.t[i] << " mm" << endl;
      }

      for (unsigned int i=0;i<ppsSpectrum->ArmF.TrkDet1.X.size();i++){
	cout << "ArmF, TrkDet1(x,y): (" << ppsSpectrum->ArmF.TrkDet1.X[i] << "," << ppsSpectrum->ArmF.TrkDet1.Y[i] << ") mm" << endl;
      }

      for (unsigned int i=0;i<ppsSpectrum->ArmF.TrkDet2.X.size();i++){
	cout << "ArmF, TrkDet2(x,y): (" << ppsSpectrum->ArmF.TrkDet2.X[i] << "," << ppsSpectrum->ArmF.TrkDet2.Y[i] << ") mm" << endl;
      }

      for (unsigned int i=0;i<ppsSpectrum->ArmF.ToFDet.X.size();i++){
	cout << "ArmF, ToFDet(x,y): (" << ppsSpectrum->ArmF.ToFDet.X[i] << "," << ppsSpectrum->ArmF.ToFDet.Y[i] << ") mm" << endl;
      }


      // ARM-B
      for (unsigned int i=0;i<ppsSpectrum->ArmB.t.size();i++){
	cout << "ArmB: t[" << i << "]: " << ppsSpectrum->ArmB.t[i] << " mm" << endl;
      }

      for (unsigned int i=0;i<ppsSpectrum->ArmB.TrkDet1.X.size();i++){
	cout << "ArmB, TrkDet1(x,y): (" << ppsSpectrum->ArmB.TrkDet1.X[i] << "," << ppsSpectrum->ArmB.TrkDet1.Y[i] << ") mm" << endl;
      }

      for (unsigned int i=0;i<ppsSpectrum->ArmB.TrkDet2.X.size();i++){
	cout << "ArmB, TrkDet2(x,y): (" << ppsSpectrum->ArmB.TrkDet2.X[i] << "," << ppsSpectrum->ArmB.TrkDet2.Y[i] << ") mm" << endl;
      }

      for (unsigned int i=0;i<ppsSpectrum->ArmB.ToFDet.X.size();i++){
	cout << "ArmB, ToFDet(x,y): (" << ppsSpectrum->ArmB.ToFDet.X[i] << "," << ppsSpectrum->ArmB.ToFDet.Y[i] << ") mm\n" << endl;
      }

    }

  }

  Handle<PPSDetector> ppsDetector;
  iEvent.getByLabel("ppssim",ppsTag_,ppsDetector);

  Handle<PPSData> ppsData;
  iEvent.getByLabel("ppssim",ppsTag_,ppsData);

  //
  // Association of CMS Vertex and PPS Vertex Reconstructed by proton TOF (10 ps).
  //
  // Idea: The choosen Vertex_CMS_z will be find from the minimum | Vertex_PPS_z - Vertex_CMS_z |. 
  //

  PPSCMSVertex.clear();

  for (unsigned int i=0;i<VertexVector.size();i++){
    if (ppsSpectrum->vtxZ.size() > 0){ 
      PPSCMSVertex.push_back(std::pair<double,double>(fabs(ppsSpectrum->vtxZ[0] - VertexVector[i]->z()), VertexVector[i]->z()));
    }else{
      PPSCMSVertex.clear();
    }
  }

  // Sorting | Vertex_PPS_z - Vertex_CMS_z |
  stable_sort(PPSCMSVertex.begin(), PPSCMSVertex.end());

  if(PPSCMSVertex.size() > 0){
    for (unsigned int i=0;i<PPSCMSVertex.size();i++){
      if ((PPSCMSVertex[0].second == VertexVector[i]->z())){
	indexGold = i;
      }
    }
  }
  else{
    indexGold = -999;
  }

  if (debug){
    // Selected Vertex CMS/PPS
    if (indexGold != -999){
      cout << "\n--GOLDEN VERTEX ASSOCIATION CMS/PPS--" << endl;
      cout << "Position (x,y,z): " << VertexVector[indexGold]->position() << " mm" << endl;
    }
    cout << "--END--\n\n" << endl;
  }

}

// ------------ Sorting Vectors, Filled in FillCollections ------------
void ExclusiveDijetsAnalysisUsingPPS::SortingObjects(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug){

  // Debug Jets Details
  bool debugdetails = false;


  // Fill PF Vectors
  if(PFVector.size()>0){
    for (unsigned int i=0;i<PFVector.size();i++){
      PFVector_pt.push_back(PFVector[i]->pt());
      PFVector_eta.push_back(PFVector[i]->eta());
    }
  }

  // Ordering Jets by pT and Fill Jet Vectors

  if(JetsVector.size()>0){

    const int JetsVectorSize = (int) JetsVector.size();
    int *sortJetsVector= new int[JetsVectorSize];
    double *vjets = new double[JetsVectorSize];

    for (int i=0; i<JetsVectorSize; i++) {
      vjets[i] = JetsVector[i]->pt();
    }

    TMath::Sort(JetsVectorSize, vjets, sortJetsVector, true);

    for (unsigned int i=0;i<JetsVector.size();i++){
      JetsVector_pt.push_back(JetsVector[sortJetsVector[i]]->pt());
      JetsVector_eta.push_back(JetsVector[sortJetsVector[i]]->eta());
      JetsVector_phi.push_back(JetsVector[sortJetsVector[i]]->phi());
    }

    if (debug){
      cout << "\n--BEGIN--" << endl;

      for (unsigned int i=0;i<JetsVector.size();i++){
	if (debugdetails) {
	  cout << "ORDERED reco::PFJets[" << sortJetsVector[i] << "]\t---> " << JetsVector[sortJetsVector[i]]->print() << endl;}
	else{
	  cout << "ORDERED reco::PFJets[" << sortJetsVector[i] << "]\t---> pT [GeV]: " << JetsVector[sortJetsVector[i]]->pt() << " | eT [GeV]: " << JetsVector[sortJetsVector[i]]->et() << " | eta: " << JetsVector[sortJetsVector[i]]->eta() << " | phi: " << JetsVector[sortJetsVector[i]]->phi() << " | Vertex: " << JetsVector[sortJetsVector[i]]->vertex() << " mm" << " | Tracks Ref: " << JetsVector[sortJetsVector[i]]->getTrackRefs().size() << endl;
	}

      }

      // Vertex
      for (unsigned int i=0;i<VertexVector.size();i++){
	cout << "reco::Vertex[" << i << "]\t---> Position: " << VertexVector[i]->position() << " mm" << endl;
      }
      cout << "--END--\n\n" << endl;
    }

    // Code defense
    if(JetsVector.size()>1){
      math::XYZTLorentzVector dijetSystem(0.,0.,0.,0.);
      dijetSystem += JetsVector[sortJetsVector[0]]->p4();
      dijetSystem += JetsVector[sortJetsVector[1]]->p4();
      Mjj = dijetSystem.M();
    }
  }

  if ( (Mjj > 0.) && (Mpf > 0.)){
    Rjj = Mjj/Mpf;
  }

}

void ExclusiveDijetsAnalysisUsingPPS::AssociateJetsWithVertex(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug){

  if (indexGold != -999){
    if (debug) {
      cout << "\n--GOLDEN VERTEX ASSOCIATION CMS/PPS--" << endl;
      cout << "Position (x,y,z): " << VertexVector[indexGold]->position() << " mm" << endl;
    }
    GoldenVertexZ = VertexVector[indexGold]->z();
  }

  for (unsigned int i=0;i<TracksVector.size();i++){
    if(indexGold != -999){
      MinimumDistance.push_back(fabs(VertexVector[indexGold]->z() - TracksVector[i]->innerPosition().Z()));
    }
  }

  const int minVectorSize = (int) MinimumDistance.size();
  int *sortMinVector= new int[minVectorSize];
  double *vmin = new double[minVectorSize];

  for (int i=0; i<minVectorSize; i++) {
    vmin[i]=MinimumDistance[i];
  }

  TMath::Sort(minVectorSize, vmin, sortMinVector, false);

  if(MinimumDistance.size()>0) {
    if (debug) {
      cout << "Minimum Distance | Tracks(z)-Vertex(z) |: " << MinimumDistance[sortMinVector[0]] << " mm" << endl;
      cout << "Maximum Distance | Tracks(z)-Vertex(z) |: " << MinimumDistance[sortMinVector[MinimumDistance.size()-1]] << " mm\n" << endl;
    }
    MinDistance = MinimumDistance[sortMinVector[0]];
    MaxDistance = MinimumDistance[sortMinVector[MinimumDistance.size()-1]];
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(ExclusiveDijetsAnalysisUsingPPS);
