// -*- C++ -*-
//
// Package:    ExclusiveDijetsAnalysisUsingPPS
// Class:      ExclusiveDijetsAnalysisUsingPPS
// 
/**\class ExclusiveDijetsAnalysisUsingPPS ExclusiveDijetsAnalysisUsingPPS.cc PPS_Dijets_Skimming/ExclusiveDijetsAnalysisUsingPPS/plugins/ExclusiveDijetsAnalysisUsingPPS.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
 */
//
// Original Author:  Sandro Fonseca De Souza
//         Created:  Sun, 13 Apr 2014 01:53:33 GMT
// $Id$
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

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "ForwardAnalysis/Utilities/interface/LargestGenRapidityGap.h"
//#include "ForwardAnalysis/Utilities/interface/CastorEnergy.h"
//#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/ExclusiveDijetsEvent.h"
//#include "ForwardAnalysis/ForwardTTreeAnalysis/interface/FWLiteTools.h"




//PPS
#include "DataFormats/PPSObjects/interface/PPSSpectrometer.h"


#include "TH1F.h"
#include "TH2F.h"

#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <map>


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


		edm::InputTag jetTag_;
		edm::InputTag particleFlowTag_;
		 std::string ppsTag_;
		//edm::InputTag ppsTag_;
		TH1F * h_pt_jet1, * h_eta_jet1, * h_phi_jet1,*h_dijetMass, *h_RJJ,*h_MassJets,* h_MxFromJets,* h_MxFromPFCands ;
		TH1F * h_pt_jet2, * h_eta_jet2, * h_phi_jet2,*PPS_xiARMPlus;	

};

////////////////////////////////////////////////////////////////////////
// constructors and destructor
//
ExclusiveDijetsAnalysisUsingPPS::ExclusiveDijetsAnalysisUsingPPS(const edm::ParameterSet& iConfig):
	jetTag_(iConfig.getParameter<edm::InputTag>("JetTag")),
	particleFlowTag_(iConfig.getParameter<edm::InputTag>("ParticleFlowTag")),
//	ppsTag_(iConfig.getParameter<edm::InputTag>("PPSTag"))
 	ppsTag_(iConfig.getUntrackedParameter<std::string>("PPSTag","PPSReco"))
{
	edm::Service<TFileService> fs;
	h_pt_jet1    = fs->make<TH1F>( "jet1_pt"  , "p_{t}", 100,  0., 300. );
	h_eta_jet1   = fs->make<TH1F>( "jet1_eta" , "#eta" , 50, -3, 3 );
	h_phi_jet1   = fs->make<TH1F>( "jet1_phi" , "#phi" , 50, -M_PI, M_PI );
	h_pt_jet2    = fs->make<TH1F>( "jet2_pt"  , "p_{t}", 100,  0., 300. );
	h_eta_jet2   = fs->make<TH1F>( "jet2_eta" , "#eta" , 50, -3, 3 );
	h_phi_jet2   = fs->make<TH1F>( "jet2_phi" , "#phi" , 50, -M_PI, M_PI );
	h_MassJets   = fs->make<TH1F>( "JetsMass" , "jets Mass" , 100, 0., 1000. );
	h_dijetMass   = fs->make<TH1F>( "diJetsMass" , "Dijets Mass" , 100, 0., 1000. );
	h_RJJ   = fs->make<TH1F>( "RJJ" , "Dijets Mass" , 100, 0., 2.0 );
	h_MxFromJets =  fs->make<TH1F>( "MxFromJets" , "MxFromJets" , 100, 0., 20. );
	h_MxFromPFCands =  fs->make<TH1F>( "MxFromPFCands" , "MxFromPFCands" , 100, 0., 20. );
	PPS_xiARMPlus =  fs->make<TH1F>( "PPS_xiARMPlus" , "PPS_xiARMPlus" , 100, -2.0, 2.0 ); 		
}

////////////////////////////////////////////////////////////////////////////////////////////////
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
void ExclusiveDijetsAnalysisUsingPPS::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace std;


	Handle<edm::View<reco::Jet> > jetCollectionH;
	iEvent.getByLabel(jetTag_,jetCollectionH);

	Handle<edm::View<reco::PFCandidate> > particleFlowCollectionH;
	iEvent.getByLabel(particleFlowTag_,particleFlowCollectionH);

//	Handle<edm::View<PPSSpectrometer> > ppsCollectionH;
//	iEvent.getByLabel(ppsTag_,ppsCollectionH);		
       	Handle<PPSSpectrometer> ppsCollectionH;
// 	iEvent.getByLabel("ppssim","PPSReco",ppsCollectionH);
	iEvent.getByLabel("ppssim",ppsTag_,ppsCollectionH);
         PPS_xiARMPlus->Fill(ppsCollectionH->ArmF.t[0]);
	//cout << " \t Nvtx = " << ppsCollectionH->Nvtx << endl;


	if(jetCollectionH->size() > 1){
		const reco::Jet& jet1 = (*jetCollectionH)[0];// they come out ordered right?
		const reco::Jet& jet2 = (*jetCollectionH)[1];

		h_pt_jet1->Fill(jet1.pt()) ;
		h_eta_jet1->Fill(jet1.eta() );
		h_phi_jet1->Fill(jet1.phi() );	

		h_pt_jet2->Fill(jet2.pt()) ;
		h_eta_jet2->Fill(jet2.eta() );
		h_phi_jet2->Fill(jet2.phi() );

		//	cout<<" jet1.pt(): "<<jet1.pt()<<" jet1.eta(): "<< jet1.eta() <<" jet1.phi(): "<<jet1.phi()<<endl;
		//	cout<<" jet2.pt(): "<<jet1.pt()<<" jet2.eta(): "<< jet2.eta() <<" jet2.phi(): "<<jet2.phi()<<endl;  

		//Using Lorentz Vector	
		math::XYZTLorentzVector jetSystem(0.,0.,0.,0.);
		jetSystem += jet1.p4();
		h_MassJets->Fill(jetSystem.M());
		// Di-jet mass
		math::XYZTLorentzVector dijetSystem(0.,0.,0.,0.);
		dijetSystem += jet1.p4();
		dijetSystem += jet2.p4();
		h_dijetMass->Fill(dijetSystem.M());

		// M_{X}
		math::XYZTLorentzVector allJets(0.,0.,0.,0.);
		for(edm::View<reco::Jet>::const_iterator jet = jetCollectionH->begin();
				jet != jetCollectionH->end(); ++jet) allJets += jet->p4();


		h_MxFromJets->Fill(allJets.M());

		math::XYZTLorentzVector allPFCands(0.,0.,0.,0.);
		for(edm::View<reco::PFCandidate>::const_iterator pfCand = particleFlowCollectionH->begin();
				pfCand != particleFlowCollectionH->end();
				++pfCand) allPFCands += pfCand->p4();


		h_MxFromPFCands->Fill(allPFCands.M());


		if(jetCollectionH->size() > 1){
			//      double RjjFromJets = Rjj(*jetCollectionH,*jetCollectionH);
			////h_RJJ->Fill(RjjFromJets);
		}
		///////////////////////////////////////////



	}//Jets selection


}//end analysis
/////////////////////////////////////////////////////

//define this as a plug-in
DEFINE_FWK_MODULE(ExclusiveDijetsAnalysisUsingPPS);
