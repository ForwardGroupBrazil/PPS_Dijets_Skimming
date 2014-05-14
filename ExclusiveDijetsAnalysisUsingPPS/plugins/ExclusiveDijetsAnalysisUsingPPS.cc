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
using namespace HepMC;

#define NELEMS(x)  (sizeof(x) / sizeof(x[0]));

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
    void ResolutionStudies(bool debug);
    void CheckSelection();

    bool MakePlots_;
    edm::InputTag jetTag_;
    edm::InputTag particleFlowTag_;
    edm::InputTag VertexTag_;
    std::string ppsTag_;
    double pTPFThresholdCharged_;
    double energyPFThresholdBar_;
    double energyPFThresholdEnd_;
    double energyPFThresholdHF_;
    double cmsVertexResolution_;
    double PPSVertexResolution_;

    int indexGold;
    int nAssociated=0;
    int checkCounter=0;
    int CheckCounterAssociator=0;

    std::vector<const reco::PFJet*> JetsVector;
    std::vector<const reco::Vertex*> VertexVector;
    std::vector<const reco::Track*> TracksVector;
    std::vector<const reco::PFCandidate*> PFVector;
    std::vector<const PPSSpectrometer*> PPSSpecVector;
    std::vector< std::pair<double,double> > PPSCMSVertex;
    std::vector<const reco::PFJet*> PFJets;
    std::vector<double> MinimumDistance;   
    std::vector< math::XYZVector > JetVertex;
    std::vector< math::XYZVector > JetsSamePosition;

    std::vector<double> JetsVector_pt;
    std::vector<double> JetsVector_eta;
    std::vector<double> JetsVector_phi;
    std::vector<double> VertexCMSVectorX;
    std::vector<double> VertexCMSVectorY;
    std::vector<double> VertexCMSVectorZ;
    std::vector<double> VertexGENVectorX;
    std::vector<double> VertexGENVectorY;
    std::vector<double> VertexGENVectorZ;
    std::vector<double> AllDiffVertexVector;
    std::vector<int> TracksPerJetVector;
    std::vector<double> JetsSameVector_pt;
    std::vector<double> JetsSameVector_eta;
    std::vector<double> JetsSameVector_phi;
    std::vector<math::XYZTLorentzVector> JetsSameVector_p4;
    std::vector<double> CandidatesJets_pt;
    std::vector<double> CandidatesJets_eta;
    std::vector<double> CandidatesJets_phi;
    std::vector<math::XYZTLorentzVector> CandidatesJets_p4;
    double VertexGEN_x;
    double VertexGEN_y;
    double VertexGEN_z;

    TTree* eventTree_;

    double JetsSameVertex_x;
    double JetsSameVertex_y;
    double JetsSameVertex_z;
    int nTracks;
    int nVertex;
    double GoldenVertexZ;
    double VertexZPPS;
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
    double xPPSArmFToF, yPPSArmFToF;
    double xPPSArmBToF, yPPSArmBToF;
    int stopPPSArmFTrkDet1, stopPPSArmFTrkDet2;
    int stopPPSArmBTrkDet1, stopPPSArmBTrkDet2;
    int stopPPSArmBToF, stopPPSArmFToF;
    double Mjj;
    double Mpf;
    double Rjj;
    double CandidatesMjj;
    double MinDistanceZVertex;
    double MaxDistanceZVertex;
    bool acceptPrint = false;

    double resol[18]={0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5};
    int counter[18]={0};
    int size_resol = NELEMS(resol);
    TH2D *h_vertex;

};


//
// constructors and destructor
//

ExclusiveDijetsAnalysisUsingPPS::ExclusiveDijetsAnalysisUsingPPS(const edm::ParameterSet& iConfig):
  MakePlots_(iConfig.getParameter<bool>("MakePlots")),
  jetTag_(iConfig.getParameter<edm::InputTag>("JetTag")),
  particleFlowTag_(iConfig.getParameter<edm::InputTag>("ParticleFlowTag")),
  VertexTag_(iConfig.getParameter<edm::InputTag>("VertexTag")),
  ppsTag_(iConfig.getUntrackedParameter<std::string>("PPSTag","PPSReco")),
  pTPFThresholdCharged_(iConfig.getParameter<double>("pTPFThresholdCharged")),
  energyPFThresholdBar_(iConfig.getParameter<double>("energyPFThresholdBar")),
  energyPFThresholdEnd_(iConfig.getParameter<double>("energyPFThresholdEnd")),
  energyPFThresholdHF_(iConfig.getParameter<double>("energyPFThresholdHF")),
  cmsVertexResolution_(iConfig.getParameter<double>("cmsVertexResolution")),
  PPSVertexResolution_(iConfig.getParameter<double>("PPSVertexResolution"))
{

  edm::Service<TFileService> fs;
  eventTree_ = fs->make<TTree>("Event","Event");
  eventTree_->Branch("JetsPt",&JetsVector_pt);
  eventTree_->Branch("JetsEta",&JetsVector_eta);
  eventTree_->Branch("JetsPhi",&JetsVector_phi);
  eventTree_->Branch("JetsSameVertex_pt",&JetsSameVector_pt);
  eventTree_->Branch("JetsSameVertex_eta",&JetsSameVector_eta);
  eventTree_->Branch("JetsSameVertex_phi",&JetsSameVector_phi);
  eventTree_->Branch("JetsSameVertex_x",&JetsSameVertex_x,"JetsSameVertex_x/D");
  eventTree_->Branch("JetsSameVertex_y",&JetsSameVertex_y,"JetsSameVertex_y/D");
  eventTree_->Branch("JetsSameVertex_z",&JetsSameVertex_z,"JetsSameVertex_z/D");
  eventTree_->Branch("CandidatesJets_pt",&CandidatesJets_pt);
  eventTree_->Branch("CandidatesJets_eta",&CandidatesJets_eta);
  eventTree_->Branch("CandidatesJets_phi",&CandidatesJets_phi);
  eventTree_->Branch("TracksPerJet",&TracksPerJetVector);
  eventTree_->Branch("VertexCMSVector_x",&VertexCMSVectorX);
  eventTree_->Branch("VertexCMSVector_y",&VertexCMSVectorY);
  eventTree_->Branch("VertexCMSVector_z",&VertexCMSVectorZ);
  eventTree_->Branch("VertexGENVector_x",&VertexGENVectorX);
  eventTree_->Branch("VertexGENVector_y",&VertexGENVectorY);
  eventTree_->Branch("VertexGENVector_z",&VertexGENVectorZ);
  eventTree_->Branch("AllDiffVertexZVector",&AllDiffVertexVector);
  eventTree_->Branch("MinDistanceZVertex",&MinDistanceZVertex,"MinDistanceZVertex/D");
  eventTree_->Branch("MaxDistanceZVertex",&MaxDistanceZVertex,"MaxDistanceZVertex/D");
  eventTree_->Branch("nVertex",&nVertex,"nVertex/I");
  eventTree_->Branch("nTracks",&nTracks,"nTracks/I");
  eventTree_->Branch("MinDistance",&MinDistance,"MinDistance/D");
  eventTree_->Branch("MaxDistance",&MaxDistance,"MaxDistance/D");
  eventTree_->Branch("GoldenVertexZ",&GoldenVertexZ,"GoldenVertexZ/D");
  eventTree_->Branch("VertexZPPS",&VertexZPPS,"VertexZPPS/D");
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
  eventTree_->Branch("xPPSArmFToF",&xPPSArmFToF,"xPPSArmFToF/D");
  eventTree_->Branch("yPPSArmFToF",&yPPSArmFToF,"yPPSArmFToF/D");
  eventTree_->Branch("xPPSArmBToF",&xPPSArmBToF,"xPPSArmBToF/D");
  eventTree_->Branch("yPPSArmBToF",&yPPSArmBToF,"yPPSArmBToF/D");
  eventTree_->Branch("stopPPSArmFTrkDet1",&stopPPSArmFTrkDet1,"stopPPSArmFTrkDet1/I");
  eventTree_->Branch("stopPPSArmFTrkDet2",&stopPPSArmFTrkDet2,"stopPPSArmFTrkDet2/I");
  eventTree_->Branch("stopPPSArmBTrkDet1",&stopPPSArmBTrkDet1,"stopPPSArmBTrkDet1/I");
  eventTree_->Branch("stopPPSArmBTrkDet2",&stopPPSArmBTrkDet2,"stopPPSArmBTrkDet2/I");
  eventTree_->Branch("stopPPSArmFToF",&stopPPSArmFToF,"stopPPSArmFToF/I");
  eventTree_->Branch("stopPPSArmBToF",&stopPPSArmBToF,"stopPPSArmBToF/I");
  eventTree_->Branch("Mjj",&Mjj,"Mjj/D");
  eventTree_->Branch("Mpf",&Mpf,"Mpf/D");
  eventTree_->Branch("Rjj",&Rjj,"Rjj/D");
  eventTree_->Branch("CandidatesMjj",&CandidatesMjj,"CandidatesMjj/D");

  if(MakePlots_){
    TFileDirectory Dir = fs->mkdir("Info");
    h_vertex = Dir.make<TH2D>("vertex", "Vertex; #sigma_{resolution, PV}; Number of Events", 17, resol, 1000, 0., 1000.);
  }

}

ExclusiveDijetsAnalysisUsingPPS::~ExclusiveDijetsAnalysisUsingPPS()
{
}


// ------------ method called once each job just before starting event loop  ------------
void ExclusiveDijetsAnalysisUsingPPS::beginJob()
{
  cout << "\n--- P P S    C E P   A N A L Y Z E R---\n" << endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void ExclusiveDijetsAnalysisUsingPPS::endJob()
{

  bool debug = false;

  cout << "\n--- S U M M A R Y---" << endl;
  cout << "# evt (Associated Vertex, full algorithm): " << nAssociated << endl;
  if(acceptPrint){
    cout << "# evt (LeadingJets_pT >50 GeV & Jets at Tracker & Fiducial PPS): " << checkCounter << endl;
    cout << "# evt (LeadingJets_pT >50 GeV & Jets at Tracker & Fiducial PPS & PPS/CMS Vertex): " << CheckCounterAssociator << endl;
  }

  if(MakePlots_) {
    for (int i=0; i<size_resol; i++){
      h_vertex->Fill(resol[i],counter[i]);
      if (debug) cout << "Resolution PPS: " << resol[i] << " | # Events: " << counter[i] << endl;
    }
  }

  cout << "\n--- E N D---\n" << endl;

}

// ------------ method called for each event  ------------
void ExclusiveDijetsAnalysisUsingPPS::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  Init(); // Clean Variables
  FillCollections(iEvent, iSetup, false); //true-> Print Outputs and Golden Vertex (PPS and CMS). False-> No print screen.
  SortingObjects(iEvent, iSetup, false); //true-> Print Ordered Outputs, False-> No print screen.
  AssociateJetsWithVertex(iEvent, iSetup, false); //true-> Print Ordered Outputs, False-> No print screen.
  ResolutionStudies(false); //true-> Print Ordered Vertex and Check Vertices.
  CheckSelection(); // Check Algorithm.
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
  JetsSameVector_pt.clear();
  JetsSameVector_eta.clear();
  JetsSameVector_phi.clear();
  JetsSameVector_p4.clear();
  JetsSamePosition.clear();
  MinimumDistance.clear();
  PFJets.clear();
  AllDiffVertexVector.clear();

  JetsVector_pt.clear();
  JetsVector_eta.clear();
  JetsVector_phi.clear();
  VertexCMSVectorX.clear();
  VertexCMSVectorY.clear();
  VertexCMSVectorZ.clear();
  VertexGENVectorX.clear();
  VertexGENVectorY.clear();
  VertexGENVectorZ.clear();
  JetVertex.clear();
  TracksPerJetVector.clear();
  CandidatesJets_pt.clear();
  CandidatesJets_eta.clear();
  CandidatesJets_phi.clear();
  CandidatesJets_p4.clear();

  JetsSameVertex_x = -999.;
  JetsSameVertex_y = -999.;
  JetsSameVertex_z = -999.;
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
  xPPSArmFToF= -999.; yPPSArmFToF = -999.;
  xPPSArmBToF=-999.; yPPSArmBToF = -999.;
  stopPPSArmFTrkDet1 = -999; stopPPSArmFTrkDet2 = -999;
  stopPPSArmBTrkDet1 = -999; stopPPSArmBTrkDet2 = -999;
  stopPPSArmFToF = -999; stopPPSArmBToF = -999;
  Mjj= -999.;
  Mpf= -999.;
  Rjj= -999.;
  VertexZPPS = -999.;
  CandidatesMjj = -999.;
  MaxDistanceZVertex = -999.;
  MaxDistanceZVertex = -999.;

}


// ------------ Fill Vectors, All Handles  ------------
void ExclusiveDijetsAnalysisUsingPPS::FillCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug)
{

  // Debug Detailed Information, each loop
  bool debugdetails = false;

  // Fill Vertex
  Handle<edm::View<reco::Vertex> > vertex;
  iEvent.getByLabel(VertexTag_, vertex);

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

  // ArmF and ArmB, ToF info (x,y)
  if(ppsSpectrum->ArmF.ToFDet.X.size() > 0){
    xPPSArmFToF = ppsSpectrum->ArmF.ToFDet.X[0];
    yPPSArmFToF = ppsSpectrum->ArmF.ToFDet.Y[0];
  }

  if(ppsSpectrum->ArmB.ToFDet.X.size() > 0){
    xPPSArmBToF = ppsSpectrum->ArmB.ToFDet.X[0];
    yPPSArmBToF = ppsSpectrum->ArmB.ToFDet.Y[0];
  }

  // ArmF and ArmB, HasStopped Info
  if(ppsSpectrum->ArmF.ToFDet.HasStopped.size() > 0){
    stopPPSArmFToF = ppsSpectrum->ArmF.ToFDet.HasStopped[0];
  }

  if(ppsSpectrum->ArmB.ToFDet.HasStopped.size() > 0){
    stopPPSArmBToF = ppsSpectrum->ArmB.ToFDet.HasStopped[0];
  }

  if(ppsSpectrum->ArmF.TrkDet1.HasStopped.size() > 0){
    stopPPSArmFTrkDet1 = ppsSpectrum->ArmF.TrkDet1.HasStopped[0];
  }

  if(ppsSpectrum->ArmB.TrkDet1.HasStopped.size() > 0){
    stopPPSArmBTrkDet1 = ppsSpectrum->ArmB.TrkDet1.HasStopped[0];
  }

  if(ppsSpectrum->ArmF.TrkDet2.HasStopped.size() > 0){
    stopPPSArmFTrkDet2 = ppsSpectrum->ArmF.TrkDet2.HasStopped[0];
  }

  if(ppsSpectrum->ArmB.TrkDet2.HasStopped.size() > 0){
    stopPPSArmBTrkDet2 = ppsSpectrum->ArmB.TrkDet2.HasStopped[0];
  }

  // PPS Vertex
  if(ppsSpectrum->vtxZ.size() > 0){
    VertexZPPS = ppsSpectrum->vtxZ[0];
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

  for (unsigned int i=0;i<VertexVector.size();i++){

    VertexCMSVectorX.push_back(VertexVector[i]->x());
    VertexCMSVectorY.push_back(VertexVector[i]->y());
    VertexCMSVectorZ.push_back(VertexVector[i]->z());

    if (ppsSpectrum->vtxZ.size() > 0){ 
      PPSCMSVertex.push_back(std::pair<double,double>(fabs(ppsSpectrum->vtxZ[0] - VertexVector[i]->z()), VertexVector[i]->z()));
    }else{
      PPSCMSVertex.clear();
    }
  }

  // Fill GEN Vertex Info
  Handle<PPSSpectrometer> ppsGEN;
  iEvent.getByLabel("ppssim","PPSGen",ppsGEN);

  if (ppsGEN->vtxZ.size() > 0){
    for (unsigned int i=0;i<ppsGEN->vtxZ.size();i++){
      VertexGENVectorX.push_back(ppsGEN->vtxX[i]);
      VertexGENVectorY.push_back(ppsGEN->vtxY[i]);
      VertexGENVectorZ.push_back(ppsGEN->vtxZ[i]);
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
      cout << "Position (x,y,z): " << VertexVector[indexGold]->position() << " cm" << endl;
    }
    cout << "--END--\n\n" << endl;
  }

}

// ------------ Sorting Vectors, Filled in FillCollections ------------
void ExclusiveDijetsAnalysisUsingPPS::SortingObjects(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug){

  // Debug Jets Details
  bool debugdetails = false;


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
	  cout << "ORDERED reco::PFJets[" << sortJetsVector[i] << "]\t---> pT [GeV]: " << JetsVector[sortJetsVector[i]]->pt() << " | eT [GeV]: " << JetsVector[sortJetsVector[i]]->et() << " | eta: " << JetsVector[sortJetsVector[i]]->eta() << " | phi: " << JetsVector[sortJetsVector[i]]->phi() << " | Vertex: " << JetsVector[sortJetsVector[i]]->vertex() << " cm" << " | Tracks Ref: " << JetsVector[sortJetsVector[i]]->getTrackRefs().size() << endl;
	}

      }

      // Vertex
      for (unsigned int i=0;i<VertexVector.size();i++){
	cout << "reco::Vertex[" << i << "]\t---> Position: " << VertexVector[i]->position() << " cm" << endl;
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

  bool debugdeep = false;

  if (indexGold != -999){
    if (debug) {
      cout << "\n--GOLDEN VERTEX ASSOCIATION CMS/PPS--" << endl;
      cout << "Position (x,y,z): " << VertexVector[indexGold]->position() << " cm" << endl;
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
      cout << "Minimum Distance | Tracks(z)-Vertex(z) |: " << MinimumDistance[sortMinVector[0]] << " cm" << endl;
      cout << "Maximum Distance | Tracks(z)-Vertex(z) |: " << MinimumDistance[sortMinVector[MinimumDistance.size()-1]] << " cm\n" << endl;
    }
    MinDistance = MinimumDistance[sortMinVector[0]];
    MaxDistance = MinimumDistance[sortMinVector[MinimumDistance.size()-1]];
  }

  // A S S O C I A T I O N   O F   J E T S / V E R T E X`
  double vx_mean = -999.;
  double vy_mean = -999.;
  double vz_mean = -999.;
  double pt2_x = 0.;
  double pt2_y = 0.;
  double pt2_z = 0.;
  double pt2 = 0.;

  edm::Handle<std::vector<reco::PFJet> > pfJetCollection;
  iEvent.getByLabel(jetTag_, pfJetCollection);

  for( unsigned k=0; k<pfJetCollection->size(); k++ ) {
    const reco::PFJet* pfjet = &(pfJetCollection->at(k));
    reco::TrackRefVector tracks = pfjet->getTrackRefs();

    if (debugdeep){
      cout << "-- Info Jet[" << k << "] | Jet pT[GeV]: " << pfjet->pt() << " | Tracks per jet: " << tracks.size() <<  endl; 
    }

    int checktrack = 0;
    for (reco::TrackRefVector::const_iterator itTrack = tracks.begin(); itTrack != tracks.end(); ++itTrack) {

      // Only three tracks
      if(checktrack > 3) break;
      ++checktrack;

      pt2 += (*itTrack)->pt()*(*itTrack)->pt();
      pt2_x += (*itTrack)->pt()*(*itTrack)->pt()*(*itTrack)->vx();
      pt2_y += (*itTrack)->pt()*(*itTrack)->pt()*(*itTrack)->vy();
      pt2_z += (*itTrack)->pt()*(*itTrack)->pt()*(*itTrack)->vz();
      if (debugdeep){
	cout << "Track-> pT [GeV]: " << (*itTrack)->pt() << " | Eta: " << (*itTrack)->eta() << " | phi: " << (*itTrack)->phi() << " | Vertex Origin (x,y,z) [cm]: " << (*itTrack)->vx() << ", " << (*itTrack)->vy() << ", " << (*itTrack)->vz() << endl;  
      }
    }

    if (pt2 > 0.) {
      vx_mean = pt2_x/pt2;
      vy_mean = pt2_y/pt2;
      vz_mean = pt2_z/pt2;
    }

    math::XYZVector coord(vx_mean,vy_mean,vz_mean);
    JetVertex.push_back(coord);
    TracksPerJetVector.push_back(tracks.size());

  }

  if(JetsVector.size()>0 && VertexVector.size()>0){

    math::XYZVector coordleadingjet(JetVertex[0].X(),JetVertex[0].Y(),JetVertex[0].Z());
    JetsSamePosition.push_back(coordleadingjet);
    JetsSameVector_pt.push_back(JetsVector[0]->pt());
    JetsSameVector_eta.push_back(JetsVector[0]->eta());
    JetsSameVector_phi.push_back(JetsVector[0]->phi());
    JetsSameVector_p4.push_back(JetsVector[0]->p4());

    JetsSameVertex_x = JetVertex[0].X();
    JetsSameVertex_y = JetVertex[0].Y();
    JetsSameVertex_z = JetVertex[0].Z();

    for(unsigned int i=1;i<JetsVector.size();i++){

      if(debug){
	cout << "\nJet pT: " << JetsVector[i]->pt() <<endl;
	cout << "<vx,vy,vz> [cm]: (" << JetVertex[i].X() << ", " << JetVertex[i].Y() << ", " << JetVertex[i].Z() << ")\n" << endl;
      }

      if(fabs(JetVertex[0].Z()-JetVertex[i].Z()) < cmsVertexResolution_){
	math::XYZVector coordjets(JetVertex[i].X(),JetVertex[i].Y(),JetVertex[i].Z());
	JetsSamePosition.push_back(coordjets);
	JetsSameVector_pt.push_back(JetsVector[i]->pt());
	JetsSameVector_eta.push_back(JetsVector[i]->eta());
	JetsSameVector_phi.push_back(JetsVector[i]->phi());
	JetsSameVector_p4.push_back(JetsVector[i]->p4());
      }

    }
  }

  // Getting at least Dijets Events
  if(JetsSameVector_pt.size()>1){

    for(unsigned int i=0;i<JetsSameVector_pt.size();i++){
      if (debug) cout << "Jet["<< i << "], pT [GeV]: " << JetsSameVector_pt[i] << " | Position (x,y,z) cm: " << JetsSamePosition[i].X() << ", " << JetsSamePosition[i].Y() << ", " << JetsSamePosition[i].Z() << " | PPS Vertex z [cm]: " << VertexZPPS << endl;

      // Fill at least one jet from the same PPS/CMS associated vertex
      if(fabs(JetsSamePosition[0].Z() - VertexZPPS) < PPSVertexResolution_){
	CandidatesJets_pt.push_back(JetsSameVector_pt[i]);
	CandidatesJets_eta.push_back(JetsSameVector_eta[i]);
	CandidatesJets_phi.push_back(JetsSameVector_phi[i]);
	CandidatesJets_p4.push_back(JetsSameVector_p4[i]);
      }
    }

    // Counter Number of Associated Events CMS/PPS Vertex
    if(fabs(JetsSamePosition[0].Z() - VertexZPPS) < PPSVertexResolution_){
      ++nAssociated;
    }

    for(int i=0; i < size_resol; i++){
      // Resolution Studies
      if(fabs(JetsSamePosition[0].Z() - VertexZPPS) < resol[i]){
	counter[i]++;
      }
    }

  }

  // Fill Mjj leading dijets from candidate CMS/PPS associated vertex
  if(CandidatesJets_p4.size()>1){
    math::XYZTLorentzVector dijetSystemCand(0.,0.,0.,0.);
    dijetSystemCand += CandidatesJets_p4[0];
    dijetSystemCand += CandidatesJets_p4[1];
    CandidatesMjj = dijetSystemCand.M();
  }

}

void ExclusiveDijetsAnalysisUsingPPS::ResolutionStudies(bool debug){

  // Put Vertex in order of Z.
  const int VertexVectorSize = (int) VertexVector.size();
  int *sortVertexVector= new int[VertexVectorSize];
  double *vvertex = new double[VertexVectorSize];

  for (int i=0; i<VertexVectorSize; i++) {
    vvertex[i] = VertexVector[i]->z();
  }

  TMath::Sort(VertexVectorSize, vvertex, sortVertexVector, true);

  const int  size = (int) VertexVector.size();

  if(VertexVector.size()>1){
    for (int i=0; i<(size-1); i++) {
      AllDiffVertexVector.push_back(fabs(VertexVector[sortVertexVector[i+1]]->z()-VertexVector[sortVertexVector[i]]->z()));
    }
  }

  // see neighboring.
  if(VertexVector.size()>0){
    for(unsigned int i=0;i<VertexVector.size();i++){
      if(debug) cout << "Vertex Ordered Z [cm]: " <<  VertexVector[sortVertexVector[i]]->z() << endl;
    }
  }

  // Distances between neighboring vertexes
  if(AllDiffVertexVector.size()>1){
    for(unsigned int i=0;i<AllDiffVertexVector.size();i++){
      if(debug) cout << "Difference Vertex Z [cm]: " <<  AllDiffVertexVector[i] << endl;
    }
  }

  // Sorting | Vertex_PPS_z - Vertex_CMS_z |
  stable_sort(AllDiffVertexVector.begin(), AllDiffVertexVector.end());

  if(AllDiffVertexVector.size()>1){
    for(unsigned int i=0;i<AllDiffVertexVector.size();i++){
      if(debug) cout << "Ordered Difference Vertex Z [cm]: " <<  AllDiffVertexVector[i] << endl;
    }

    MinDistanceZVertex = AllDiffVertexVector[0];
    MaxDistanceZVertex = AllDiffVertexVector[AllDiffVertexVector.size()-1];

  }

}

void ExclusiveDijetsAnalysisUsingPPS::CheckSelection(){

  acceptPrint = true;

  double xmax = -3.15;
  double xmin = -23.15;
  double xmax2 = -2.03;
  double xmin2 = -22.03;
  double ymax = 9.0;
  double ymin = -9.0;

  bool cutXdet1 = ((xmin<xPPSArmFDet1 && xPPSArmFDet1<xmax) && (xmin<xPPSArmBDet1 && xPPSArmBDet1<xmax));
  bool cutXdet2 = ((xmin2<xPPSArmFDet2 && xPPSArmFDet2<xmax2) && (xmin2<xPPSArmBDet2 && xPPSArmBDet2<xmax2));
  bool cutYdet1 = ((ymin<yPPSArmFDet1 && yPPSArmFDet1<ymax) && (ymin<yPPSArmBDet1 && yPPSArmBDet1<ymax));
  bool cutYdet2 = ((ymin<yPPSArmFDet2 && yPPSArmFDet2<ymax) && (ymin<yPPSArmBDet2 && yPPSArmBDet2<ymax));

  bool stopTrkDet1 = (stopPPSArmFTrkDet1==0 && stopPPSArmBTrkDet1 ==0);
  bool stopTrkDet2 = (stopPPSArmFTrkDet2==0 && stopPPSArmBTrkDet2 ==0);

  if(JetsVector.size()>1){
    if(JetsVector_pt[0]>50. && JetsVector_pt[1]>50.){
      if(JetsVector_eta[0]>-2. && JetsVector_eta[1]<2.){
	if(cutXdet1 && cutXdet2 && cutYdet1 && cutYdet2){
	  if(stopTrkDet1 && stopTrkDet2){
	    ++checkCounter;
	    if(fabs(VertexVector[0]->z()-VertexZPPS)<PPSVertexResolution_){
	      ++CheckCounterAssociator;
	    }
	  }
	}
      }
    }
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(ExclusiveDijetsAnalysisUsingPPS);
