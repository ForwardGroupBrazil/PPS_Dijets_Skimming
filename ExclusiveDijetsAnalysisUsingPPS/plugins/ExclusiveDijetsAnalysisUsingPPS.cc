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

// root
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

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


    // ----------member data ---------------------------

    void FillCollections(const edm::Event&, const edm::EventSetup&, bool debug);
    void SortingObjects(const edm::Event&, const edm::EventSetup&, bool debug);

    edm::InputTag jetTag_;
    edm::InputTag particleFlowTag_;
    std::string ppsTag_;

    std::vector<const reco::Jet*> JetsVector;
    std::vector<const reco::Vertex*> VertexVector;
    std::vector<const PPSSpectrometer*> PPSSpecVector;

};

//
// constructors and destructor
//

ExclusiveDijetsAnalysisUsingPPS::ExclusiveDijetsAnalysisUsingPPS(const edm::ParameterSet& iConfig):
  jetTag_(iConfig.getParameter<edm::InputTag>("JetTag")),
  particleFlowTag_(iConfig.getParameter<edm::InputTag>("ParticleFlowTag")),
  ppsTag_(iConfig.getUntrackedParameter<std::string>("PPSTag","PPSReco"))
{

  edm::Service<TFileService> fs;

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

  FillCollections(iEvent, iSetup, true); //true-> Print Ordered Outputs, False-> No print screen. Debugger. There is a deep debugger inside the function. 
  SortingObjects(iEvent, iSetup, false); //true-> Print Ordered Outputs, False-> No print screen. Debugger.

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
  VertexVector.clear();

  if(vertex->size()>0){
    for(itVertex=0; itVertex < vertexsize; ++itVertex){
      const reco::Vertex* vertexAll = &((*vertex)[itVertex]);
      VertexVector.push_back(vertexAll);
    }
  }


  // Fill Jets
  Handle<edm::View<reco::Jet> > jets;
  iEvent.getByLabel(jetTag_,jets);

  int jetsize = jets->size();
  int itJets;
  JetsVector.clear();

  if(jets->size()>0){
    for(itJets=0; itJets < jetsize; ++itJets){
      const reco::Jet* jetAll = &((*jets)[itJets]);
      JetsVector.push_back(jetAll);
    }
  }

  // Fill PPS Spectrometer
  Handle<PPSSpectrometer> ppsSpectrum;
  iEvent.getByLabel("ppssim",ppsTag_,ppsSpectrum);

  if (debug){
    cout << "\n--PPS INFO--" << endl;
    if (ppsSpectrum->vtxZ.size() > 0) cout << "vtxZ[0]: " << ppsSpectrum->vtxZ[0] << " | Vector Size: " << ppsSpectrum->vtxZ.size() << endl;
    if (ppsSpectrum->ArmF.t.size() > 0) cout << "ArmF.t[0]: " << ppsSpectrum->ArmF.t[0] << " | Vector Size: " << ppsSpectrum->ArmF.t.size() << endl;
    if (ppsSpectrum->ArmB.t.size() > 0) cout << "ArmB.t[0]: " << ppsSpectrum->ArmB.t[0] << " | Vector Size: " << ppsSpectrum->ArmB.t.size() << endl;  
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

    cout << "--END--\n\n" << endl;

  }


  Handle<PPSDetector> ppsDetector;
  iEvent.getByLabel("ppssim",ppsTag_,ppsDetector);

  Handle<PPSData> ppsData;
  iEvent.getByLabel("ppssim",ppsTag_,ppsData);


}

// ------------ Sorting Vectors, Filled in FillCollections ------------
void ExclusiveDijetsAnalysisUsingPPS::SortingObjects(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug){

  // Ordering Jets by pT
  if(JetsVector.size()>0){

    const int JetsVectorSize = (int) JetsVector.size();
    int *sortJetsVector= new int[JetsVectorSize];
    double *vjets = new double[JetsVectorSize];

    for (int i=0; i<JetsVectorSize; i++) {
      vjets[i] = JetsVector[i]->pt();
    }

    TMath::Sort(JetsVectorSize, vjets, sortJetsVector, true);

    if (debug){
      cout << "\n--BEGIN--" << endl;

      for (unsigned int i=0;i<JetsVector.size();i++){
	cout << "ORDERED reco::Jets[" << sortJetsVector[i] << "]\t---> pT [GeV]: " << JetsVector[sortJetsVector[i]]->pt() << " | eT [GeV]: " << JetsVector[sortJetsVector[i]]->et() << " | eta: " << JetsVector[sortJetsVector[i]]->eta() << " | phi: " << JetsVector[sortJetsVector[i]]->phi() << " | Vertex: " << JetsVector[sortJetsVector[i]]->vertex() << " mm" << endl;
      }

      // Vertex
      for (unsigned int i=0;i<VertexVector.size();i++){
	cout << "reco::Vertex[" << i << "]\t---> Position: " << VertexVector[i]->position() << " mm" << endl;
      }

      cout << "--END--\n\n" << endl;
    }
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(ExclusiveDijetsAnalysisUsingPPS);
