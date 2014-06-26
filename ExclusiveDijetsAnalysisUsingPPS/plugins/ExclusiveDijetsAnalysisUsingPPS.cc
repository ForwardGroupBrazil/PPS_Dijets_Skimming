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


//GenJets
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
// Gen Particles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

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
#include "DataFormats/Math/interface/LorentzVector.h"

//PPS
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
//Event Info
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"



// root
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"
//#include "TLorentzVector.h"

// c++
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <map>

//random gauss to ToF smearing
#include "TRandom3.h"

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
                void MCGenInfo(const edm::Event&, const edm::EventSetup&, bool debug);
		void FillCollections(const edm::Event&, const edm::EventSetup&, bool debug);
		void SortingObjects(const edm::Event&, const edm::EventSetup&, bool debug);
		void AssociateJetsWithVertex(const edm::Event&, const edm::EventSetup&, bool debug);
		void ResolutionStudies(bool debug);
		void PPSSelection();
		void GenCollections(const edm::Event&, const edm::EventSetup&, bool debug);

		bool MakePlots_;
		bool runWithWeightGen_;
		edm::InputTag  genjets_ ;
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
		double EBeam_;

		int indexGold,indexGoldToF_00,indexGoldToF_0i,indexGoldToF_ii;

		int nAssociated=0;
		int checkCounter=0;
		int CheckCounterAssociator=0;

		std::vector<const reco::GenJet*> GenJetsVector;
		std::vector<const reco::PFJet*> JetsVector;
		std::vector<const reco::Vertex*> VertexVector;
		std::vector<const reco::Track*> TracksVector;
		std::vector<const reco::PFCandidate*> PFVector;
		std::vector<const PPSSpectrometer*> PPSSpecVector;
		std::vector< std::pair<double,double> > PPSCMSVertex;

		std::vector< std::pair<double,double> >PPSCMSVertexToF_00;
		std::vector< std::pair<double,double> >PPSCMSVertexToF_0i;
		std::vector< std::pair<double,double> >PPSCMSVertexToF_ii;

		std::vector< std::pair<double,double> >PPSCMSVertexToF_00_0i;
		std::vector< std::pair<double,double> >PPSCMSVertexToF_00_0i_ii;

		std::vector<const reco::PFJet*> PFJets;
		std::vector<double> MinimumDistance;
		std::vector<double> DistanceBetweenJets;   
		std::vector< math::XYZVector > JetVertex;
		std::vector< math::XYZVector > JetsSamePosition;

		std::vector<double> GenJetsVector_pt;
		std::vector<double> GenJetsVector_eta;
		std::vector<double> GenJetsVector_phi;

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
		// Vector 
		std::vector<double> xiPPSArmFInfo;
		std::vector<double> xiPPSArmBInfo;
		std::vector<double> tPPSArmFInfo;
		std::vector<double> tPPSArmBInfo;
		// Arm Forward det1 and det2
		std::vector<double> xPPSArmFDet1Info;
		std::vector<double> yPPSArmFDet1Info;

		std::vector<double> xPPSArmFDet2Info;
		std::vector<double> yPPSArmFDet2Info;

		// Arm  Backward det1 and det2
		std::vector<double> xPPSArmBDet1Info;
		std::vector<double> yPPSArmBDet1Info;

		std::vector<double> xPPSArmBDet2Info;
		std::vector<double> yPPSArmBDet2Info;

		// ToF Informations for Forward arm
		std::vector<double> xPPSArmFToFInfo;
		std::vector<double> yPPSArmFToFInfo;

		// Backward arm
		std::vector<double> xPPSArmBToFInfo;
		std::vector<double> yPPSArmBToFInfo;


		std::vector<double>  stopPPSArmFToFInfo;
		std::vector<double>  stopPPSArmBToFInfo;

		std::vector<double>  stopPPSArmFTrkDet1Info;
		std::vector<double>  stopPPSArmFTrkDet2Info;
		std::vector<double>  stopPPSArmBTrkDet1Info;
		std::vector<double>  stopPPSArmBTrkDet2Info;

		// Gen Proton Info
		std::vector<const reco::GenParticle*> GenProtonVectorInfo;
		std::vector< math::XYZTLorentzVector > protonLorentzVector;

		// Vertex position on Z using ToF
		std::vector<double> VertexZPPSToF_00;
		std::vector<double> VertexZPPSToF_0i;
		std::vector<double> VertexZPPSToF_ii;
		std::vector<double> VertexZPPSToF;
		std::vector<double> VertexZPPSToF_00_0i;

		std::vector<double> DijetsVertexZPPSToF;

		// Gen Vertex 
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
		double GenRjj,Rjj;
		double GenMjj;
		double CandidatesMjj;
		double MinDistanceZVertex;
		double MaxDistanceZVertex;
		bool acceptPrint = false;
		double GenMxx,Mx;
		bool FiducialCut = false;

		double resol[18]={0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5};
		int counter[18]={0};
		int size_resol = NELEMS(resol);
		TH2D *h_vertex;

		// xi and t gen 
		double xigen_plus,xigen_minus,tgen_plus,tgen_minus;                

		// ToF detector
		double deltaToF_00;
		double deltaToF_0i;
		double deltaToF_ii;

		double tof30ps_deltaToF_00;
		double tof30ps_deltaToF_0i;
		double tof30ps_deltaToF_ii;

		double ToF_Mx_00;
		double ToF_Mx_0i;
		double ToF_Mx_ii;

		double VertexCombineToF_00;
		double VertexCombineToF_0i;
		double VertexCombineToF_ii;

		map<int,double> BeamX;
		map<int,double> BeamY;
		map<int,double> BeamX_RMS;
		map<int,double> BeamY_RMS;
		float ToFXPosSigma;

		int ToFCellId(double x, double y);
		void setToFGeometry(std::string);
		void setQuartic();
		//void setDiamond();
		void setToFParameters(std::string geo,double xpos);
		void setBeamParameters(int pos, double x,double y,double rms_x,double rms_y);

		std::string ToFGeometryType="";
		double      ToFXPos;
		vector<pair<double,double> > ToFCellColumn;
		vector<pair<double,double> > ToFCellRow;
		const double RPWindowThickness=0.3;
		const int   NYCell = 4; // number of cell in Y for Quartic
		const int   NXCell = 5; // number of cell in X for Quartic
		const int   NXCellLowerRow = 1; // number of cell in X for lower row for Diamond
		const int   NXCellUpperRow = 1; // number of cell in X for upper row for Diamond
		const int   NXCellCentralRow = 14; // number of cell in X for central row for Diamond
		const int   DiamondNYCell = 3; // number of cell in Y for Diamond
		const float celws[14]= {0.3, 0.4, 0.45, 0.55, 0.65, 0.75, 0.9, 1.0, 1.3, 1.7, 2, 2.4, 3, 4.6};
		const float DiamondCellHeight = 5.;
		const float QuarticCellHeight = 3.;
		const float QuarticCellWidth = 3.;
		const float DiamondUpperCellW = 10;
		//const vector<float> DiamondCentralCellW(celws,0.3);
		//const vector<float> DiamondCentralCellW(celws,celws + sizeof(celws) / sizeof(float));
		const float DiamondLowerCellW = 10;
		int n1pu =0; int n2pu = 0;

                double GeneratorWeight, Pthat;
                int eventNumber, runNumber;

};


//
// constructors and destructor
//

ExclusiveDijetsAnalysisUsingPPS::ExclusiveDijetsAnalysisUsingPPS(const edm::ParameterSet& iConfig):
	MakePlots_(iConfig.getParameter<bool>("MakePlots")),
	runWithWeightGen_(iConfig.getParameter<bool>("RunWithWeightGen")),
	genjets_(iConfig.getParameter<edm::InputTag>("GenJets")), 
	jetTag_(iConfig.getParameter<edm::InputTag>("JetTag")),
	particleFlowTag_(iConfig.getParameter<edm::InputTag>("ParticleFlowTag")),
	VertexTag_(iConfig.getParameter<edm::InputTag>("VertexTag")),
	ppsTag_(iConfig.getUntrackedParameter<std::string>("PPSTag","PPSReco")),
	pTPFThresholdCharged_(iConfig.getParameter<double>("pTPFThresholdCharged")),
	energyPFThresholdBar_(iConfig.getParameter<double>("energyPFThresholdBar")),
	energyPFThresholdEnd_(iConfig.getParameter<double>("energyPFThresholdEnd")),
	energyPFThresholdHF_(iConfig.getParameter<double>("energyPFThresholdHF")),
	cmsVertexResolution_(iConfig.getParameter<double>("cmsVertexResolution")),
	PPSVertexResolution_(iConfig.getParameter<double>("PPSVertexResolution")),
	EBeam_(iConfig.getParameter<double>("EBeam"))
{

	edm::Service<TFileService> fs;
	eventTree_ = fs->make<TTree>("Event","Event");
	//eventTree_->Branch("GenJetsVector",&GenJetsVector);

	eventTree_->Branch("GenJetsPt",&GenJetsVector_pt);
	eventTree_->Branch("GenJetsEta",&GenJetsVector_eta);
	eventTree_->Branch("GenJetsPhi",&GenJetsVector_phi);

	eventTree_->Branch("JetsPt",&JetsVector_pt);
	eventTree_->Branch("JetsEta",&JetsVector_eta);
	eventTree_->Branch("JetsPhi",&JetsVector_phi);

	eventTree_->Branch("JetVertex",&JetVertex);
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
	eventTree_->Branch("DistanceBetweenJets",&DistanceBetweenJets);
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

	eventTree_->Branch("GenMjj",&GenMjj,"GenMjj/D");
	eventTree_->Branch("xigen_plus",&xigen_plus,"xigen_plus/D");
	eventTree_->Branch("xigen_minus",&xigen_minus,"xigen_minus/D");
	eventTree_->Branch("tgen_plus",&tgen_plus,"tgen_plus/D");
	eventTree_->Branch("tgen_minus",&tgen_minus,"tgen_minus/D");

	eventTree_->Branch("Mpf",&Mpf,"Mpf/D");
	eventTree_->Branch("Rjj",&Rjj,"Rjj/D");

	eventTree_->Branch("GenRjj",&Rjj,"GenRjj/D");
	eventTree_->Branch("CandidatesMjj",&CandidatesMjj,"CandidatesMjj/D");
	eventTree_->Branch("Mx",&Mx,"Mx/D");

	eventTree_->Branch("GenMxx",&GenMxx,"GenMxx/D");

        eventTree_->Branch("GeneratorWeight",&GeneratorWeight,"GeneratorWeight/D");        
        eventTree_->Branch("Pthat",&Pthat,"Pthat/D");
        
        eventTree_->Branch("eventNumber",&eventNumber,"eventNumber/I");  
        eventTree_->Branch("runNumber",&runNumber,"runNumber/I");             

	eventTree_->Branch("FiducialCut",&FiducialCut,"FiducialCut/B");
	eventTree_->Branch("deltaToF_00",&deltaToF_00,"deltaToF_00/D");
	eventTree_->Branch("deltaToF_0i",&deltaToF_0i,"deltaToF_0i/D");
	eventTree_->Branch("deltaToF_ii",&deltaToF_ii,"deltaToF_ii/D");
	eventTree_->Branch("tof30ps_deltaToF_00",&tof30ps_deltaToF_00,"tof30ps_deltaToF_00/D");
	eventTree_->Branch("tof30ps_deltaToF_0i",&tof30ps_deltaToF_0i,"tof30ps_deltaToF_0i/D");
	eventTree_->Branch("tof30ps_deltaToF_ii",&tof30ps_deltaToF_ii,"tof30ps_deltaToF_ii/D");
	eventTree_->Branch("ToF_Mx_00",&ToF_Mx_00,"ToF_Mx_00/D");
	eventTree_->Branch("ToF_Mx_0i",&ToF_Mx_0i,"ToF_Mx_0i/D");
	eventTree_->Branch("ToF_Mx_ii",&ToF_Mx_ii,"ToF_Mx_ii/D");

	///vector to storage of information of Arm Forward for signal and PU protons
	// xi 
	eventTree_->Branch("xiPPSArmFInfo",&xiPPSArmFInfo);
	eventTree_->Branch("xiPPSArmBInfo",&xiPPSArmBInfo);
	// t 
	eventTree_->Branch("tPPSArmFInfo",&tPPSArmFInfo);
	eventTree_->Branch("tPPSArmBInfo",&tPPSArmBInfo);
	// x  and y position in detector Det1 Arm F
	eventTree_->Branch("xPPSArmFDet1Info",&xPPSArmFDet1Info);
	eventTree_->Branch("yPPSArmFDet1Info",&yPPSArmFDet1Info);
	// x  and y position in detector Det2 Arm F
	eventTree_->Branch("xPPSArmFDet2Info",&xPPSArmFDet2Info);
	eventTree_->Branch("yPPSArmFDet2Info",&yPPSArmFDet2Info);
	// x  and y position in detector Det1 Arm B
	eventTree_->Branch("xPPSArmBDet1Info",&xPPSArmBDet1Info);
	eventTree_->Branch("yPPSArmBDet1Info",&yPPSArmBDet1Info);
	// x  and y position in detector Det2 Arm B
	eventTree_->Branch("xPPSArmBDet2Info",&xPPSArmBDet2Info);
	eventTree_->Branch("yPPSArmBDet2Info",&yPPSArmBDet2Info);
	// Tof info Arm Forward
	eventTree_->Branch("xPPSArmFToFInfo",&xPPSArmFToFInfo);
	eventTree_->Branch("yPPSArmFToFInfo",&yPPSArmFToFInfo);
	// Arm Backward
	eventTree_->Branch("xPPSArmBToFInfo",&xPPSArmBToFInfo);
	eventTree_->Branch("yPPSArmBToFInfo",&yPPSArmBToFInfo);

	//Stop ToF    
	eventTree_->Branch("stopPPSArmFToFInfo",&stopPPSArmFToFInfo);
	eventTree_->Branch("stopPPSArmBToFInfo",&stopPPSArmBToFInfo);
	// Arm TRK FWD
	eventTree_->Branch("stopPPSArmFTrkDet1Info",&stopPPSArmFTrkDet1Info);
	eventTree_->Branch("stopPPSArmFTrkDet2Info",&stopPPSArmFTrkDet2Info);
	//Arm TRK BCW
	eventTree_->Branch("stopPPSArmBTrkDet1Info",&stopPPSArmBTrkDet1Info);
	eventTree_->Branch("stopPPSArmBTrkDet2Info",&stopPPSArmBTrkDet2Info);

	// ToF Position combination on Z axis
	//PPSCMSVertexToF_00 sinal + sinal  PPSCMSVertexToF_0i signal + PU PPSCMSVertexToF_ii PU + PU
	eventTree_->Branch("PPSCMSVertexToF_00",&PPSCMSVertexToF_00);
	eventTree_->Branch("PPSCMSVertexToF_0i",&PPSCMSVertexToF_0i);
	eventTree_->Branch("PPSCMSVertexToF_ii",&PPSCMSVertexToF_ii);
	eventTree_->Branch("PPSCMSVertexToF_00_0i",&PPSCMSVertexToF_00_0i);
	eventTree_->Branch("PPSCMSVertexToF_00_0i_ii",&PPSCMSVertexToF_00_0i_ii);

	// VertexZPPSToF
	eventTree_->Branch("VertexZPPSToF",&VertexZPPSToF);
	eventTree_->Branch("VertexZPPSToF_00",&VertexZPPSToF_00);
	eventTree_->Branch("VertexZPPSToF_0i",&VertexZPPSToF_0i);
	eventTree_->Branch("VertexZPPSToF_ii",&VertexZPPSToF_ii);
	eventTree_->Branch("VertexZPPSToF_00_0i",&VertexZPPSToF_00_0i);

	//
	//DijetsVertexZPPSToF (Distance between jet vertex and ToF Z vertex
	eventTree_->Branch("DijetsVertexZPPSToF",&DijetsVertexZPPSToF);

	// Proton Info      
	//	eventTree_->Branch("GenProtonVectorInfo",&GenProtonVectorInfo);
	eventTree_->Branch("ProtonsP4",&protonLorentzVector);

        

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
		cout << "# evt (LeadingJets_pT >100 GeV & Jets at Tracker & Fiducial PPS): " << checkCounter << endl;
		cout << "# evt (LeadingJets_pT >100 GeV & Jets at Tracker & Fiducial PPS & PPS/CMS Vertex): " << CheckCounterAssociator << endl;
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
        MCGenInfo(iEvent, iSetup, false);
	GenCollections(iEvent, iSetup, false);
	FillCollections(iEvent, iSetup, false); //true-> Print Outputs and Golden Vertex (PPS and CMS). False-> No print screen.
	SortingObjects(iEvent, iSetup, false); //true-> Print Ordered Outputs, False-> No print screen.
	AssociateJetsWithVertex(iEvent, iSetup, false); //true-> Print Ordered Outputs, False-> No print screen.
	ResolutionStudies(false); //true-> Print Ordered Vertex and Check Vertices.
	PPSSelection(); // Check Algorithm.
	eventTree_->Fill();

}

// ------------ Clean Variables --------------

void ExclusiveDijetsAnalysisUsingPPS::Init(){

	GenJetsVector.clear();
	JetsVector.clear();
	VertexVector.clear();
	TracksVector.clear();
	PFVector.clear();
	PPSSpecVector.clear();
	PPSCMSVertex.clear();

	//ToF Vertex position in Z axis combinations : sig + sign (00), sign + PU(0i), PU+PU (ii) 
	PPSCMSVertexToF_00.clear();
	PPSCMSVertexToF_0i.clear();
	PPSCMSVertexToF_ii.clear();

	PPSCMSVertexToF_00_0i.clear();
	PPSCMSVertexToF_00_0i_ii.clear();


	JetsSameVector_pt.clear();
	JetsSameVector_eta.clear();
	JetsSameVector_phi.clear();
	JetsSameVector_p4.clear();
	JetsSamePosition.clear();
	MinimumDistance.clear();
	PFJets.clear();
	AllDiffVertexVector.clear();
	DistanceBetweenJets.clear();

	GenJetsVector_pt.clear();
	GenJetsVector_eta.clear();
	GenJetsVector_phi.clear();

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
	deltaToF_00 = -999.;
	deltaToF_0i = -999.;
	deltaToF_ii = -999.;
	ToF_Mx_00 =-999.;
	ToF_Mx_0i =-999.;
	ToF_Mx_ii =-999.;

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
	// Mjj gen jet level
	GenMjj = -999.;
	GenRjj = -999.;

	VertexZPPS = -999.;
	CandidatesMjj = -999.;
	MaxDistanceZVertex = -999.;
	MaxDistanceZVertex = -999.;
	Mx = -999.;
	GenMxx = -999.;

	FiducialCut = false;
	//t and xi gen level
	xigen_plus = -999.;
	xigen_minus = -999.;
	tgen_plus = -999.;
	tgen_minus = -999.;

	VertexCombineToF_00 = -999.;
	VertexCombineToF_0i = -999.;
	VertexCombineToF_ii = -999.;

	// Clearing vectors
	xiPPSArmFInfo.clear();
	xiPPSArmBInfo.clear();
	tPPSArmFInfo.clear();
	tPPSArmBInfo.clear();

	xPPSArmFDet1Info.clear();
	yPPSArmFDet1Info.clear();

	xPPSArmFDet2Info.clear();
	yPPSArmFDet2Info.clear();

	xPPSArmBDet1Info.clear();
	yPPSArmBDet1Info.clear();

	xPPSArmBDet2Info.clear();
	yPPSArmBDet2Info.clear();

	//
	xPPSArmFToFInfo.clear();
	yPPSArmFToFInfo.clear();
	xPPSArmBToFInfo.clear();
	yPPSArmBToFInfo.clear();

	stopPPSArmFToFInfo.clear(); 
	stopPPSArmBToFInfo.clear();

	stopPPSArmFTrkDet1Info.clear();
	stopPPSArmFTrkDet2Info.clear();
	stopPPSArmBTrkDet1Info.clear();
	stopPPSArmBTrkDet2Info.clear();
	// Gen particle information
	GenProtonVectorInfo.clear();
	protonLorentzVector.clear(); 
	// allGenParticles.clear();

	VertexZPPSToF_00.clear();
	VertexZPPSToF_0i.clear();
	VertexZPPSToF_ii.clear();
	VertexZPPSToF.clear();
	VertexZPPSToF_00_0i.clear();
	DijetsVertexZPPSToF.clear();
      
        Pthat = -1.;
        GeneratorWeight = -1.;
        runNumber = -1;
        eventNumber = -1;

}
//----------------------------MC Gen Ifnfo---------------------
// MCGenInfo
void ExclusiveDijetsAnalysisUsingPPS::MCGenInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug)
{
unsigned int eventNumber_ = iEvent.id().event();
unsigned int runNumber_ = iEvent.id().run();

eventNumber = eventNumber_;
runNumber = runNumber_;


if( runWithWeightGen_ ){
      edm::Handle<GenEventInfoProduct> genEventInfoH;
      iEvent.getByLabel("generator", genEventInfoH);
      Pthat = genEventInfoH->binningValues()[0] ;
      GeneratorWeight= genEventInfoH->weight() ;
   } else {
      Pthat= -1. ;
      GeneratorWeight= -1. ;
   }

}


//------------------Fill Gen Jet information-------------------------------------------
void ExclusiveDijetsAnalysisUsingPPS::GenCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug)
{

	math::XYZTLorentzVector allGenParticles(0.,0.,0.,0.);

	// Accessing Gen Particles Information

	Handle<reco::GenParticleCollection> genParticle;
	iEvent.getByLabel("genParticles",genParticle);
	// const reco::GenParticleCollection& genParticles = *genParticlesCollectionH; 

	int gensize = genParticle->size();
	int itGen;

	if(genParticle->size()>0){
		for(itGen=0; itGen < gensize; ++itGen){
			const reco::GenParticle* genAll = &((*genParticle)[itGen]);
			if (genAll->status() != 1) continue;
			allGenParticles += genAll->p4();
			if (genAll->pdgId() != 2212) continue;
			GenProtonVectorInfo.push_back(genAll);
		}
		GenMxx = allGenParticles.mass();
	}
	// Saving One or Two Leading Protons
	if (GenProtonVectorInfo.size()==1){
		protonLorentzVector.push_back(GenProtonVectorInfo[0]->p4());
	}
	if (GenProtonVectorInfo.size()>1){
		protonLorentzVector.push_back(GenProtonVectorInfo[0]->p4());
		protonLorentzVector.push_back(GenProtonVectorInfo[1]->p4());
	}
	if (protonLorentzVector.size() > 0 && protonLorentzVector.at(0).pz()>0.){

		xigen_plus = ( 1. - (protonLorentzVector.at(0).pz()/EBeam_) );
		math::XYZTLorentzVector vec_pi(0.,0.,EBeam_,EBeam_);
		math::XYZTLorentzVector vec_pf(protonLorentzVector.at(0).px(),protonLorentzVector.at(0).py(),protonLorentzVector.at(0).pz(),protonLorentzVector.at(0).energy());
		math::XYZTLorentzVector vec_t = (vec_pf - vec_pi);

		tgen_plus = vec_t.mag2(); 
		tgen_plus = abs(tgen_plus);
		if(debug){
			cout << "--> xi, plus: " << xigen_plus << endl;
			cout << "--> t, plus (abs) in Z positive : " << tgen_plus << endl;
			cout << " proton pZ ( proton 1):-> "  << protonLorentzVector.at(0).pz() << endl;
			cout << " proton py (proton 1):-> "  << protonLorentzVector.at(0).py() << endl;
			cout << " proton px (proton 1):-> "  << protonLorentzVector.at(0).px() << endl;

		}
	} // Gen proton plus side

	///////////////////
	////////Gen protons minus side

	if (protonLorentzVector.size() > 0 && protonLorentzVector.at(1).pz()<0.){

		xigen_minus = ( 1. + (protonLorentzVector.at(1).pz()/EBeam_) );
		math::XYZTLorentzVector vec_pi(0.,0.,-EBeam_,EBeam_);
		math::XYZTLorentzVector vec_pf(protonLorentzVector.at(1).px(),protonLorentzVector.at(1).py(),protonLorentzVector.at(1).pz(),protonLorentzVector.at(1).energy());
		math::XYZTLorentzVector vec_t = (vec_pf - vec_pi);
		tgen_minus = vec_t.mag2();
		tgen_minus = abs(tgen_minus);
		if(debug){
			cout << "--> xi, minus: " << xigen_minus << endl;
			cout << "--> t, minus (abs) in Z negative: " << tgen_minus << endl;
			cout << " proton pZ ( proton 2):-> "  << protonLorentzVector.at(1).pz() << endl;
			cout << " proton py (proton 2):-> "  << protonLorentzVector.at(1).py() << endl;
			cout << " proton px (proton 2):-> "  << protonLorentzVector.at(1).px() << endl;
		}
	} // Gen proton minus side
	////////////////////////////////////////////////////////////

	//Gen Jets Information
	Handle<reco::GenJetCollection> genjets;
	iEvent.getByLabel(genjets_,genjets);
	int Genjetsize = genjets->size();
	int itGenJets;

	if(debug)  cout<< "Gen Jets size: "<< Genjetsize << endl;
	if(genjets->size()>0){
		for(itGenJets=0; itGenJets < Genjetsize; ++itGenJets){
			const reco::GenJet* GenjetAll = &((*genjets)[itGenJets]);
			GenJetsVector.push_back(GenjetAll);
		}
	}

	// Ordering Jets by pT and Fill Jet Vectors
	if(GenJetsVector.size()>0){

		const int GenJetsVectorSize = (int) GenJetsVector.size();
		int *sortGenJetsVector= new int[GenJetsVectorSize];
		double *genvjets = new double[GenJetsVectorSize];

		for (int i=0; i< GenJetsVectorSize; i++) {
			genvjets[i] = GenJetsVector[i]->pt();
		}

		TMath::Sort(GenJetsVectorSize, genvjets, sortGenJetsVector, true);

		for (unsigned int i=0;i<GenJetsVector.size();i++){
			GenJetsVector_pt.push_back(GenJetsVector[sortGenJetsVector[i]]->pt());
			GenJetsVector_eta.push_back(GenJetsVector[sortGenJetsVector[i]]->eta());
			GenJetsVector_phi.push_back(GenJetsVector[sortGenJetsVector[i]]->phi());
		}

		if (debug) {   
			cout << "\n--BEGIN--" << endl;

			for (unsigned int i=0;i<GenJetsVector.size();i++){
				if (debug) {
					cout << "ORDERED reco::GenJets[" << sortGenJetsVector[i] << "]\t---> " << GenJetsVector[sortGenJetsVector[i]]->print() << endl;}
				else{
					cout << "ORDERED reco::GenJets[" << sortGenJetsVector[i] << "]\t---> pT [GeV]: " << GenJetsVector[sortGenJetsVector[i]]->pt() << " | eT [GeV]: " << GenJetsVector[sortGenJetsVector[i]]->et() << " | eta: " << GenJetsVector[sortGenJetsVector[i]]->eta() << " | phi: " << GenJetsVector[sortGenJetsVector[i]]->phi() << " | Vertex: " << GenJetsVector[sortGenJetsVector[i]]->vertex() << " cm" << endl;
				}
			} // Loop debugger
		}
	}

	//////////Gen Mjj ////////////////////////////////////////////////////
	if(Genjetsize < 2) return;
	const reco::GenJet* genJet1 = &(*genjets)[0];
	const reco::GenJet* genJet2 = &(*genjets)[1];

	if(genJet1&&genJet2){
		math::XYZTLorentzVector dijetGenSystem(0.,0.,0.,0.);
		dijetGenSystem += genJet1->p4();
		dijetGenSystem += genJet2->p4();
		double massGen = dijetGenSystem.M();
		GenMjj=massGen;
		if(debug) cout << ">>> Leading Jet pt,eta: " << genJet1->pt() << " , " << genJet1->eta() << endl;

		if(debug) cout << ">>> Second leading Jet pt,eta: " << genJet2->pt() << " , " << genJet2->eta() << endl;

		if(debug) cout << ">>> Dijets Gen Mass " << GenMjj << endl; 
	}


	////////////////////////////////////////////////////////////

	if ( (GenMjj > 0.) && (GenMxx > 0.))
	{
		GenRjj = GenMjj/GenMxx;

		if(debug) cout << ">>> Gen RJJ " << GenRjj << endl;
	}
}//end function


// ------------ Fill Vectors, All Handles  ------------
void ExclusiveDijetsAnalysisUsingPPS::FillCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug)
{
	// Debug Detailed Information, each loop
	bool debugdetails = true;

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
		for (unsigned int k=0;k<ppsSpectrum->ArmB.xi.size();k++){
			xiPPSArmBInfo.push_back(ppsSpectrum->ArmB.xi[k]);
			// cout <<"ppsSpectrum->ArmB.xi[k]: " << ppsSpectrum->ArmB.xi[k] << " k: "<< k << endl;
		}
	}
	if(ppsSpectrum->ArmF.xi.size() > 0){
		xiPPSArmF = ppsSpectrum->ArmF.xi[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmF.xi.size();k++){
			xiPPSArmFInfo.push_back(ppsSpectrum->ArmF.xi[k]);
			// cout <<"ppsSpectrum->ArmF.xi[k]: " << ppsSpectrum->ArmF.xi[k] << " k: "<< k << endl;
			//cout <<  "xiPPSArmFInfo: " << " k " << k  << " " <<  xiPPSArmFInfo.at(k) << endl;
		}
	}
	if(ppsSpectrum->ArmB.t.size() > 0){
		tPPSArmB = ppsSpectrum->ArmB.t[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmB.t.size();k++){
			tPPSArmBInfo.push_back(ppsSpectrum->ArmB.t[k]);
			// cout <<"ppsSpectrum->ArmF.xi[k]: " << ppsSpectrum->ArmF.xi[k] << " k: "<< k << endl;
			//cout <<  "xiPPSArmFInfo: " << " k " << k  << " " <<  xiPPSArmFInfo.at(k) << endl;
		}
	}
	if(ppsSpectrum->ArmF.t.size() > 0){
		tPPSArmF = ppsSpectrum->ArmF.t[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmF.t.size();k++){
			tPPSArmFInfo.push_back(ppsSpectrum->ArmF.t[k]);
			// cout <<"ppsSpectrum->ArmF.xi[k]: " << ppsSpectrum->ArmF.xi[k] << " k: "<< k << endl;
			//cout <<  "xiPPSArmFInfo: " << " k " << k  << " " <<  xiPPSArmFInfo.at(k) << endl;
		}
	}
	// ArmF and ArmB, Det1 info (x,y)
	if(ppsSpectrum->ArmF.TrkDet1.X.size() > 0){
		xPPSArmFDet1 = ppsSpectrum->ArmF.TrkDet1.X[0];
		yPPSArmFDet1 = ppsSpectrum->ArmF.TrkDet1.Y[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmF.TrkDet1.X.size();k++){
			xPPSArmFDet1Info.push_back(ppsSpectrum->ArmF.TrkDet1.X[k]);
			yPPSArmFDet1Info.push_back(ppsSpectrum->ArmF.TrkDet1.Y[k]);
		}
	}
	if(ppsSpectrum->ArmB.TrkDet1.X.size() > 0){
		xPPSArmBDet1 = ppsSpectrum->ArmB.TrkDet1.X[0];
		yPPSArmBDet1 = ppsSpectrum->ArmB.TrkDet1.Y[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmB.TrkDet1.X.size();k++){
			xPPSArmBDet1Info.push_back(ppsSpectrum->ArmB.TrkDet1.X[k]);
			yPPSArmBDet1Info.push_back(ppsSpectrum->ArmB.TrkDet1.Y[k]);
		}
	}
	// ArmF and ArmB, Det2 info (x,y)
	if(ppsSpectrum->ArmF.TrkDet2.X.size() > 0){
		xPPSArmFDet2 = ppsSpectrum->ArmF.TrkDet2.X[0];
		yPPSArmFDet2 = ppsSpectrum->ArmF.TrkDet2.Y[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmF.TrkDet2.X.size();k++){
			xPPSArmFDet2Info.push_back(ppsSpectrum->ArmF.TrkDet2.X[k]);
			yPPSArmFDet2Info.push_back(ppsSpectrum->ArmF.TrkDet2.Y[k]);
		}
	}
	if(ppsSpectrum->ArmB.TrkDet2.X.size() > 0){
		xPPSArmBDet2 = ppsSpectrum->ArmB.TrkDet2.X[0];
		yPPSArmBDet2 = ppsSpectrum->ArmB.TrkDet2.Y[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmB.TrkDet2.X.size();k++){
			xPPSArmBDet2Info.push_back(ppsSpectrum->ArmB.TrkDet2.X[k]);
			yPPSArmBDet2Info.push_back(ppsSpectrum->ArmB.TrkDet2.Y[k]);
		}
	}
	// ArmF and ArmB, ToF info (x,y)
	if(ppsSpectrum->ArmF.ToFDet.X.size() > 0){
		xPPSArmFToF = ppsSpectrum->ArmF.ToFDet.X[0];
		yPPSArmFToF = ppsSpectrum->ArmF.ToFDet.Y[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmF.ToFDet.X.size();k++){
			xPPSArmFToFInfo.push_back(ppsSpectrum->ArmF.ToFDet.X[k]);
			yPPSArmFToFInfo.push_back(ppsSpectrum->ArmF.ToFDet.Y[k]);    
		}
	}
	if(ppsSpectrum->ArmB.ToFDet.X.size() > 0){
		xPPSArmBToF = ppsSpectrum->ArmB.ToFDet.X[0];
		yPPSArmBToF = ppsSpectrum->ArmB.ToFDet.Y[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmB.ToFDet.X.size();k++){
			xPPSArmBToFInfo.push_back(ppsSpectrum->ArmB.ToFDet.X[k]);
			yPPSArmBToFInfo.push_back(ppsSpectrum->ArmB.ToFDet.Y[k]);
		}
	}
	// ArmF and ArmB, ToF vs ptLeading
	//*****************************ToF tof30ps ******************************//
	// Gaussian Random numbers for smearing
	double smearToF1 = gRandom->Gaus(0,30e-3);
	double smearToF2 = gRandom->Gaus(0,30e-3);
	Handle<PPSSpectrometer> ppsSim;
	iEvent.getByLabel("ppssim","PPSSim",ppsSim);
	//cout << smearToF <<  "  " <<  ppsSim->ArmB.ToF[0] << "  " <<  ppsSim->ArmB.ToF[0]+smearToF <<    " " <<  ppsSpectrum->ArmB.ToF[0] << endl;
	//*****************************ToF tof30ps ******************************//

	//cout << "============== SETTING PPS BEAM PARAMETERS ==============" << endl;
	setBeamParameters(1,0.00013,-0.00062,0.186,0.495);
	setBeamParameters(2,-0.00020,0.00011,0.115,0.436);
	setBeamParameters(3,-0.00011,-0.00258,0.113,0.419);
	setToFParameters("Quartic",15);  //cout << " ---- Quartic --- " << endl;
	if(ppsSpectrum->ArmB.ToFDet.HasStopped.size() > 0 && ppsSpectrum->ArmF.ToFDet.HasStopped.size() > 0){
		unsigned int nhits =  ppsSpectrum->ArmB.ToFDet.NHits; 
		//cout << " nhitsB = " <<  ppsSpectrum->ArmB.ToFDet.NHits<< " nhitsF = " <<  ppsSpectrum->ArmF.ToFDet.NHits << endl;
		//cout << " XB = " <<  ppsSpectrum->ArmB.ToFDet.X.size() << " XF = " <<  ppsSpectrum->ArmF.ToFDet.X.size() << endl;
		//cout << " YB = " <<  ppsSpectrum->ArmB.ToFDet.Y.size() << " YF = " <<  ppsSpectrum->ArmF.ToFDet.Y.size() << endl;
		int n1pu2 =0; int n2pu2 = 0;
		for (unsigned int iB=0;iB<nhits;iB++){
			int cellIdB = ToFCellId(ppsSpectrum->ArmB.ToFDet.X.at(iB),ppsSpectrum->ArmB.ToFDet.Y.at(iB)); 
			for (unsigned int iF=0;iF<nhits;iF++){
				int cellIdF = ToFCellId(ppsSpectrum->ArmF.ToFDet.X.at(iF),ppsSpectrum->ArmF.ToFDet.Y.at(iF));
				if(cellIdB==0 || cellIdF==0) continue;  
				if(iB == 0 && iF == 0 && ppsSpectrum->ArmB.ToF.at(0) != 0 && ppsSpectrum->ArmF.ToF.at(0) != 0 ){ 
					deltaToF_00 = (ppsSpectrum->ArmF.ToF.at(iF)) - (ppsSpectrum->ArmB.ToF.at(iB));  
					tof30ps_deltaToF_00 = (ppsSim->ArmF.ToF.at(iF)+smearToF1) - (ppsSim->ArmB.ToF.at(iB)+smearToF2);  
					if(ppsSpectrum->ArmB.xi[iB] > 0. && ppsSpectrum->ArmF.xi[iF] > 0.){ 
						ToF_Mx_00 = EBeam_*TMath::Sqrt(ppsSpectrum->ArmB.xi[iB]*ppsSpectrum->ArmF.xi[iF]);
					}
					continue;
				}
				if((iB==0&&iF > 0 && ppsSpectrum->ArmB.ToF.at(0) != 0 && ppsSpectrum->ArmF.ToF.at(iF) != 0 )|| (iB > 0 && iF == 0 && ppsSpectrum->ArmB.ToF.at(iB) != 0 && ppsSpectrum->ArmF.ToF.at(0) != 0 )) {	  
					deltaToF_0i=(ppsSpectrum->ArmF.ToF.at(iF)-ppsSpectrum->ArmB.ToF.at(iB)) ; 
					tof30ps_deltaToF_0i = (ppsSim->ArmF.ToF.at(iF)+smearToF1) - (ppsSim->ArmB.ToF.at(iB)+smearToF2);  
					if(ppsSpectrum->ArmB.xi[iB] > 0. && ppsSpectrum->ArmF.xi[iF] > 0.){
						ToF_Mx_0i = EBeam_*TMath::Sqrt(ppsSpectrum->ArmB.xi[iB]*ppsSpectrum->ArmF.xi[iF]);
					} 
					//cout << "pelo menos 1 PU pronton in one arm " << ++n1pu2 << endl;
					continue; 
				}   
				if (iB > 0 && iF > 0 && ppsSpectrum->ArmB.ToF.at(iB) != 0 && ppsSpectrum->ArmF.ToF.at(iF) != 0 ) {	  
					deltaToF_ii=(ppsSpectrum->ArmF.ToF.at(iF)-ppsSpectrum->ArmB.ToF.at(iB));
					tof30ps_deltaToF_ii = (ppsSim->ArmF.ToF.at(iF)+smearToF1) - (ppsSim->ArmB.ToF.at(iB)+smearToF2);  
					if(ppsSpectrum->ArmB.xi[iB] > 0. && ppsSpectrum->ArmF.xi[iF] > 0.){
						ToF_Mx_ii = EBeam_*TMath::Sqrt(ppsSpectrum->ArmB.xi[iB]*ppsSpectrum->ArmF.xi[iF]);}
					//cout << "pelo menos 1 PU pronton in both arm " << ++n2pu2 << endl;
					continue;
				}   
			}
		}
		// cout << FiducialCut << endl;
		if (FiducialCut && n1pu2>0) ++n1pu;
		if (FiducialCut && n2pu2>0) ++n2pu;
		//cout << n1pu << " " << n2pu << endl;
	}
	// ArmF and ArmB, HasStopped Info
	if(ppsSpectrum->ArmF.ToFDet.HasStopped.size() > 0){
		stopPPSArmFToF = ppsSpectrum->ArmF.ToFDet.HasStopped[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmF.ToFDet.HasStopped.size();k++){
			stopPPSArmFToFInfo.push_back(ppsSpectrum->ArmF.ToFDet.HasStopped[k]);
		} 
	}
	if(ppsSpectrum->ArmB.ToFDet.HasStopped.size() > 0){
		stopPPSArmBToF = ppsSpectrum->ArmB.ToFDet.HasStopped[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmB.ToFDet.HasStopped.size();k++){
			stopPPSArmBToFInfo.push_back(ppsSpectrum->ArmB.ToFDet.HasStopped[k]);
		}
	}
	if(ppsSpectrum->ArmF.TrkDet1.HasStopped.size() > 0){
		stopPPSArmFTrkDet1 = ppsSpectrum->ArmF.TrkDet1.HasStopped[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmF.TrkDet1.HasStopped.size();k++){
			stopPPSArmFTrkDet1Info.push_back(ppsSpectrum->ArmF.TrkDet1.HasStopped[k]);
		}
	}
	if(ppsSpectrum->ArmB.TrkDet1.HasStopped.size() > 0){
		stopPPSArmBTrkDet1 = ppsSpectrum->ArmB.TrkDet1.HasStopped[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmB.TrkDet1.HasStopped.size();k++){
			stopPPSArmBTrkDet1Info.push_back(ppsSpectrum->ArmB.TrkDet1.HasStopped[k]);
		}
	}
	if(ppsSpectrum->ArmF.TrkDet2.HasStopped.size() > 0){
		stopPPSArmFTrkDet2 = ppsSpectrum->ArmF.TrkDet2.HasStopped[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmF.TrkDet2.HasStopped.size();k++){
			stopPPSArmFTrkDet2Info.push_back(ppsSpectrum->ArmF.TrkDet2.HasStopped[k]);
		}
	}
	if(ppsSpectrum->ArmB.TrkDet2.HasStopped.size() > 0){
		stopPPSArmBTrkDet2 = ppsSpectrum->ArmB.TrkDet2.HasStopped[0];
		for (unsigned int k=0;k<ppsSpectrum->ArmB.TrkDet2.HasStopped.size();k++){
			stopPPSArmBTrkDet2Info.push_back(ppsSpectrum->ArmB.TrkDet2.HasStopped[k]);
		}
	}
	// PPS Vertex
	if(ppsSpectrum->vtxZ.size() > 0){
		VertexZPPS = ppsSpectrum->vtxZ[0];//mudar isso para a combinação do ToF position para todos os casos
	}
	//  cout << "VertexZPPS: " << VertexZPPS << endl; 
	// PPS Mx
	if(ppsSpectrum->ArmF.xi[0] > 0. && ppsSpectrum->ArmB.xi[0] > 0.){
		Mx = EBeam_*TMath::Sqrt(ppsSpectrum->ArmF.xi[0]*ppsSpectrum->ArmB.xi[0]);
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
	double c = 299792458.0 ;// m / s
	//double c = 1.0;
	VertexCombineToF_00 = deltaToF_00 * c; // Units: (nano meter) 
	VertexCombineToF_0i = deltaToF_0i * c; //Units: (nano meter) 
	VertexCombineToF_ii = deltaToF_ii * c; //Units: (nano meter) 

	VertexCombineToF_00 = (VertexCombineToF_00 * pow(10,-7))/2.0; // In cm
	VertexCombineToF_0i = (VertexCombineToF_0i * pow(10,-7))/2.0; // in cm
	VertexCombineToF_ii = (VertexCombineToF_ii * pow(10,-7))/2.0; // in cm

	if(debug)  cout << "VertexZPPS: " << VertexZPPS << endl;
	if(debug)  cout << " VertexCombineToF_00: " << VertexCombineToF_00 << endl;
	if(debug) cout << " VertexCombineToF_0i: " << VertexCombineToF_0i << endl;
	if (debug) cout << " VertexCombineToF_ii: " << VertexCombineToF_ii << endl;

	//   VertexZPPSToF_all_Combination 
	if ( VertexCombineToF_00 != -999){
		VertexZPPSToF_00.push_back(VertexCombineToF_00);
	}
	if ( VertexCombineToF_0i != -999){
		VertexZPPSToF_0i.push_back(VertexCombineToF_0i);
	}
	if ( VertexCombineToF_ii != -999){
		VertexZPPSToF_ii.push_back(VertexCombineToF_ii);
	}
	if(debug) cout << "VertexVector.size(): " << VertexVector.size() << endl;
	if(debug) cout << " VertexZPPSToF_00.size: " <<  VertexZPPSToF_00.size() << endl;
	if(debug) cout << " VertexZPPSToF_0i.size: " <<  VertexZPPSToF_0i.size() << endl;
	if(debug) cout << " VertexZPPSToF_ii.size: " <<  VertexZPPSToF_ii.size() << endl;

	for (unsigned int i=0;i<VertexVector.size();i++){
		VertexCMSVectorX.push_back(VertexVector[i]->x());
		VertexCMSVectorY.push_back(VertexVector[i]->y());
		VertexCMSVectorZ.push_back(VertexVector[i]->z());

		if (ppsSpectrum->vtxZ.size() > 0){
			PPSCMSVertex.push_back(std::pair<double,double>(fabs(ppsSpectrum->vtxZ[0] - VertexVector[i]->z()), VertexVector[i]->z()));
		}else{
			PPSCMSVertex.clear();
		}


		///00 (signal+ signal)        
		if (VertexZPPSToF_00.size() > 0){ 
			for (unsigned int j=0;i<PPSCMSVertexToF_00.size();j++){
				//			PPSCMSVertex.push_back(std::pair<double,double>(fabs(ppsSpectrum->vtxZ[0] - VertexVector[i]->z()), VertexVector[i]->z()));
				PPSCMSVertexToF_00.push_back(std::pair<double,double>(fabs(VertexZPPSToF_00[j] - VertexVector[i]->z()), VertexVector[i]->z()));
				//cout << "VertexCombineToF_00-> " << VertexCombineToF_00 << " Vertex CMS in z axis: " <<  VertexVector[i]->z() << endl;
			}
		}else{
			//			PPSCMSVertex.clear();
			PPSCMSVertexToF_00.clear();
		}
		/// 0i (signal + PU)
		if (VertexZPPSToF_0i.size() > 0){
			for (unsigned int j=0;i<PPSCMSVertexToF_0i.size();j++){ 
				PPSCMSVertexToF_0i.push_back(std::pair<double,double>(fabs(VertexZPPSToF_0i[j] - VertexVector[i]->z()), VertexVector[i]->z()));
				//cout << "VertexCombineToF_00-> " << VertexCombineToF_00 << " Vertex CMS in z axis: " <<  VertexVector[i]->z() << endl;
			}
		}else{
			PPSCMSVertexToF_0i.clear();
		}
		// ii (PU + PU)
		if (VertexZPPSToF_ii.size() > 0){
			for (unsigned int j=0;i<PPSCMSVertexToF_ii.size();j++){	
				PPSCMSVertexToF_ii.push_back(std::pair<double,double>(fabs(VertexZPPSToF_ii[j] - VertexVector[i]->z()), VertexVector[i]->z()));     
				//cout << "VertexCombineToF_00-> " << VertexCombineToF_00 << " Vertex CMS in z axis: " <<  VertexVector[i]->z() << endl;
			}
		}else{
			PPSCMSVertexToF_ii.clear();
		}
		///////////////////////////////////////////////////////
	} //Vertex Loop

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
	//////////////////
	// Sorting | Vertex_PPSToF_z - Vertex_CMS_z | for sinal+ sinal
	stable_sort(PPSCMSVertexToF_00.begin(), PPSCMSVertexToF_00.end());
	if(PPSCMSVertexToF_00.size() > 0){
		for (unsigned int i=0;i<PPSCMSVertexToF_00.size();i++){
			if ((PPSCMSVertexToF_00[0].second == VertexVector[i]->z())){
				indexGoldToF_00 = i;
			}
		}
	}
	else{
		indexGoldToF_00 = -999;
	}
	if (debug){
		// Selected Vertex CMS/PPS
		if (indexGoldToF_00 != -999){
			cout << "\n--GOLDEN VERTEX ASSOCIATION CMS/PPS using Tof combination for signal + signal--" << endl;
			cout << "Position (x,y,z): " << VertexVector[indexGoldToF_00]->position() << " cm" << " indexGoldToF_00: " << indexGoldToF_00 <<  endl;
		}
		cout << "--END--\n\n" << endl;
	}
	/////////////////////////////
	//////////////////
	// Sorting | Vertex_PPSToF_z - Vertex_CMS_z | for signal + PU
	stable_sort(PPSCMSVertexToF_0i.begin(), PPSCMSVertexToF_0i.end());

	if(PPSCMSVertexToF_0i.size() > 0){
		for (unsigned int i=0;i<PPSCMSVertexToF_0i.size();i++){
			if ((PPSCMSVertexToF_0i[0].second == VertexVector[i]->z())){
				indexGoldToF_0i = i;
			}
		}
	}
	else{
		indexGoldToF_0i = -999;
	}

	if (debug){
		// Selected Vertex CMS/PPS
		if (indexGoldToF_0i != -999){
			cout << "\n--GOLDEN VERTEX ASSOCIATION CMS/PPS using Tof combination for signal + PU bkg --" << endl;
			cout << "Position (x,y,z): " << VertexVector[indexGoldToF_0i]->position() << " cm" << " indexGoldToF_0i: " << indexGoldToF_0i << endl;
		}
		cout << "--END--\n\n" << endl;
	}
	/////////////////////////////
	//////////////////
	// Sorting | Vertex_PPSToF_z - Vertex_CMS_z | for PU + PU
	stable_sort(PPSCMSVertexToF_ii.begin(), PPSCMSVertexToF_ii.end());

	if(PPSCMSVertexToF_ii.size() > 0){
		for (unsigned int i=0;i<PPSCMSVertexToF_ii.size();i++){
			if ((PPSCMSVertexToF_ii[0].second == VertexVector[i]->z())){
				indexGoldToF_ii = i;
			}
		}
	}
	else{
		indexGoldToF_ii = -999;
	}
	if (debug){
		// Selected Vertex CMS/PPS
		if (indexGoldToF_ii != -999){
			cout << "\n--GOLDEN VERTEX ASSOCIATION CMS/PPS using Tof combination for PU bkg + PU bkg --" << endl;
			cout << "Position (x,y,z): " << VertexVector[indexGoldToF_ii]->position() << " cm"  << " indexGoldToF_ii: " << indexGoldToF_ii << endl;
		}
		cout << "--END--\n\n" << endl;
	}
	/////////////////////////////
	if (debug)cout <<" PPSCMSVertex.size(): " << PPSCMSVertex.size() << endl;
	if (debug)cout <<" PPSCMSVertexToF_00.size(): " << PPSCMSVertexToF_00.size() << endl;
	if (debug)cout <<" PPSCMSVertexToF_0i.size(): " << PPSCMSVertexToF_0i.size() << endl;
	if (debug)cout <<" PPSCMSVertexToF_ii.size(): " << PPSCMSVertexToF_ii.size() << endl;
	//VertexZPPSToF_all_Combination
	// The main idea should merge of information for 3 cases signal + signal, signal+ PU and PU + PU
	// Using function merge : http://www.cplusplus.com/reference/algorithm/merge/
	//http://stackoverflow.com/questions/20694791/how-to-union-two-sorted-vectors-and-combine-overlapping-elements

	//  First step Merge signal + signal and signal + PU =  PPSCMSVertexToF_00_0i
	// Merge last new vector ( PPSCMSVertexToF_00_0i) + PPSCMSVertexToF_ii = PPSCMSVertexToF_00_0i_ii
	// Step1:
	merge(begin(PPSCMSVertexToF_00),end(PPSCMSVertexToF_00),begin(PPSCMSVertexToF_0i),end(PPSCMSVertexToF_0i),inserter(PPSCMSVertexToF_00_0i,PPSCMSVertexToF_00_0i.begin()));
	//merge(begin(PPSCMSVertexToF_00),end(PPSCMSVertexToF_00),begin(PPSCMSVertexToF_0i),end(PPSCMSVertexToF_0i),PPSCMSVertexToF_00_0i.begin());    
	stable_sort(PPSCMSVertexToF_00_0i.begin(), PPSCMSVertexToF_00_0i.end());
	//Step2
	merge(begin(PPSCMSVertexToF_00_0i),end(PPSCMSVertexToF_00_0i),begin(PPSCMSVertexToF_ii),end(PPSCMSVertexToF_ii),inserter(PPSCMSVertexToF_00_0i_ii,PPSCMSVertexToF_00_0i_ii.begin()));
	stable_sort(PPSCMSVertexToF_00_0i_ii.begin(), PPSCMSVertexToF_00_0i_ii.end());

	//////////////
	// Merging
	//VertexZPPSToF_ii
	//VertexZPPSToF_00
	//VertexZPPSToF_0i
	// merge 1
	merge(begin(VertexZPPSToF_00),end(VertexZPPSToF_00),begin(VertexZPPSToF_0i),end(VertexZPPSToF_0i),inserter(VertexZPPSToF_00_0i,VertexZPPSToF_00_0i.begin()));
	stable_sort(VertexZPPSToF_00_0i.begin(), VertexZPPSToF_00_0i.end());
	//merge 2
	merge(begin(VertexZPPSToF_00_0i),end(VertexZPPSToF_00_0i),begin(VertexZPPSToF_ii),end(VertexZPPSToF_ii),inserter(VertexZPPSToF,VertexZPPSToF.begin()));
	stable_sort(VertexZPPSToF.begin(), VertexZPPSToF.end());

}

// ------------ Sorting Vectors, Filled in FillCollections ------------
void ExclusiveDijetsAnalysisUsingPPS::SortingObjects(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool debug){

	// Debug Jets Details
	bool debugdetails = false;
	// REsoluaçnao 30 ps
	// resoluçao 2 mm

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

	Handle<edm::View<reco::Vertex> > vertexCollection;
	iEvent.getByLabel(VertexTag_, vertexCollection);

	for( unsigned k=0; k<pfJetCollection->size(); k++ ) {
		const reco::PFJet* pfjet = &(pfJetCollection->at(k));
		reco::TrackRefVector tracks = pfjet->getTrackRefs();

		if (debugdeep){
			cout << "-- Info Jet[" << k << "] | Jet pT[GeV]: " << pfjet->pt() << " | Tracks per jet: " << tracks.size() <<  endl; 
		}
		// A.P.
		if (debugdeep)
			cout << "-- Primary vertex (x,y,z) [cm]: " << VertexVector[0]->x() << ", " << VertexVector[0]->y() << ", " << VertexVector[0]->z() <<  endl; 

		vx_mean = -999.;
		vy_mean = -999.;
		vz_mean = -999.;
		pt2_x = 0.;
		pt2_y = 0.;
		pt2_z = 0.;
		pt2 = 0.;
		int checktrack = 0;
		for (reco::TrackRefVector::const_iterator itTrack = tracks.begin(); itTrack != tracks.end(); ++itTrack) {

			// A.P.
			// Check from which vertex track comes from
			reco::TrackBaseRef trackBaseRef( *itTrack );
			float bestweight = 0.;
			int i_vtx_track = -1;

			unsigned idx_vertex = 0;
			edm::View<reco::Vertex>::const_iterator vertex = vertexCollection->begin();
			edm::View<reco::Vertex>::const_iterator vertex_end = vertexCollection->end();
			for(; vertex != vertex_end; ++vertex, ++idx_vertex){
				reco::Vertex::trackRef_iterator it_trk = vertex->tracks_begin();
				for(; it_trk != vertex->tracks_end(); ++it_trk){
					const reco::TrackBaseRef& baseRef = *it_trk;
					if( baseRef == trackBaseRef ) {
						float w = vertex->trackWeight(baseRef);
						if (w > bestweight){
							bestweight = w;
							i_vtx_track = idx_vertex;
						}
					}
				} 
			}

			if( i_vtx_track == -1) continue;

			// Only three tracks
			++checktrack;
			if(checktrack > 3) break;

			pt2 += (*itTrack)->pt()*(*itTrack)->pt();
			pt2_x += (*itTrack)->pt()*(*itTrack)->pt()*(*itTrack)->vx();
			pt2_y += (*itTrack)->pt()*(*itTrack)->pt()*(*itTrack)->vy();
			pt2_z += (*itTrack)->pt()*(*itTrack)->pt()*(*itTrack)->vz();

			if (debugdeep){
				cout << "Track-> pT [GeV]: " << (*itTrack)->pt() << " | Eta: " << (*itTrack)->eta() << " | phi: " << (*itTrack)->phi() << " | Vertex Origin (x,y,z) [cm]: " << (*itTrack)->vx() << ", " << (*itTrack)->vy() << ", " << (*itTrack)->vz() << " | offline vertex index: " << i_vtx_track << endl;  
			}
		}

		// A.P.
		if (debugdeep){
			cout << "-- Jet pt2: " << pt2 << endl
				<< "   pt2_x, pt2_y, pt2_z: " << pt2_x << ", " << pt2_y << ", " << pt2_z << endl;
		}

		if (pt2 > 0.) {
			vx_mean = pt2_x/pt2;
			vy_mean = pt2_y/pt2;
			vz_mean = pt2_z/pt2;
		}

		// A.P.
		if (debugdeep)
			cout << "-- Jet vertex (x,y,z) [cm]: " << vx_mean << ", " << vy_mean << ", " << vz_mean << endl; 

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

			DistanceBetweenJets.push_back(fabs(JetVertex[0].Z()-JetVertex[i].Z()));      

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
	///////////////////
	// Getting at least Dijets Events
	if(JetsSameVector_pt.size()>1){

		for(unsigned int i=0;i<JetsSameVector_pt.size();i++){
			for( unsigned int j=0; j < VertexZPPSToF.size(); j++){
				if (debug) cout << "Jet["<< i << "], pT [GeV]: " << JetsSameVector_pt[i] << " | Position (x,y,z) cm: " << JetsSamePosition[i].X() << ", " << JetsSamePosition[i].Y() << ", " << JetsSamePosition[i].Z() << " | PPS Vertex z [cm]: " << VertexZPPSToF[j] << endl;
				// Fill at least one jet from the same PPS/CMS associated vertex
				DijetsVertexZPPSToF.push_back(fabs(JetsSamePosition[0].Z() - VertexZPPSToF[j]));
			} //vretex loop
		}//jetloop

	}//dijets
	//http://www.cplusplus.com/reference/algorithm/stable_sort/
	stable_sort(DijetsVertexZPPSToF.begin(), DijetsVertexZPPSToF.end());
	// Getting Dijets candidates Events
	if(JetsSameVector_pt.size()>1){
		for(unsigned int i=0;i<JetsSameVector_pt.size();i++){
			// Fill at least one jet from the same PPS/CMS associated vertex
			if( DijetsVertexZPPSToF[0] < PPSVertexResolution_){
				CandidatesJets_pt.push_back(JetsSameVector_pt[i]);
				CandidatesJets_eta.push_back(JetsSameVector_eta[i]);
				CandidatesJets_phi.push_back(JetsSameVector_phi[i]);
				CandidatesJets_p4.push_back(JetsSameVector_p4[i]);
			}
		} //Jetloop
		// Counter Number of Associated Events CMS/PPS Vertex
		if(DijetsVertexZPPSToF[0] < PPSVertexResolution_){
			++nAssociated;
		}
		for(int i=0; i < size_resol; i++){
			// Resolution Studies
			if(DijetsVertexZPPSToF[0] < resol[i]){    	
				counter[i]++;
			}
		}
	}//dijJet selection
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

void ExclusiveDijetsAnalysisUsingPPS::PPSSelection(){

	acceptPrint = true;
	double pt_cut = 100.; //GeV 

	// Fiducial Cuts
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

	// Check Some Selectiona
	if(JetsVector.size()>1){
		if(JetsVector_pt[0]>pt_cut && JetsVector_pt[1]>pt_cut){
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

	// Applying PPS Fiducial if Requested
	if(cutXdet1 && cutXdet2 && cutYdet1 && cutYdet2){
		if(stopTrkDet1 && stopTrkDet2){
			FiducialCut = true;
		}
	}
}

void ExclusiveDijetsAnalysisUsingPPS::setToFGeometry(std::string geotype)
{
	ToFCellRow.clear();
	ToFCellColumn.clear();
	if (geotype=="Quartic") {
		ToFGeometryType = "Quartic";
		setQuartic();
	} else if (geotype=="Diamond") {
		ToFGeometryType = "Diamond";
		//setDiamond();
	} else {
		std::cout << "Unknown ToF geometry." << std::endl;
	}
}
void ExclusiveDijetsAnalysisUsingPPS::setQuartic()
{
	// the vertical positions starts from the negative(bottom) to the positive(top) corner
	for(int i=0;i<NYCell;i++) {
		// vector index points to the row number from below
		double y1=QuarticCellHeight*(i-NYCell/2.);
		double y2=QuarticCellHeight*(i-NYCell/2.+1);
		ToFCellRow.push_back(pair<double,double>(y1,y2));
	}
	// vector index points to the column number
	for(int i=0;i<NXCell;i++) {
		double x1 = ToFXPos-(QuarticCellWidth*i);
		double x2 = ToFXPos-(QuarticCellWidth*(i+1));
		ToFCellColumn.push_back(pair<double,double>(x1,x2));
	}
}
int ExclusiveDijetsAnalysisUsingPPS::ToFCellId(double x, double y)
{
	int y_idx,x_idx;
	// first, get the row number
	unsigned int i;
	unsigned int start_idx=0;
	unsigned int end_idx=ToFCellRow.size();
	for(i=0;i<ToFCellRow.size();i++){
		if (y>=ToFCellRow.at(i).first&&y<=ToFCellRow.at(i).second) break;
	}
	if (i>=ToFCellRow.size()) return 0;
	y_idx = i+1;
	if (ToFGeometryType=="Quartic") {
		start_idx=0;end_idx=ToFCellColumn.size();
		// cout <<  start_idx + end_idx << endl;
	}
	/*   if (ToFGeometryType=="Diamond") {
	     switch (y_idx) {
	     case 1: start_idx=0; end_idx=1;break;
	     case 3: start_idx=DiamondCentralCellW.size()+1;end_idx=ToFCellColumn.size();break;
	     case 2: start_idx=1;end_idx=1+DiamondCentralCellW.size();break;
	     default: cout  << "ERROR: unknown ToF row index" << endl;return 0;
	     }
	     }*/
	for(i=start_idx;i<end_idx;i++) {
		if (x<=ToFCellColumn.at(i).first&&x>ToFCellColumn.at(i).second) break;
	}
	if (i>=end_idx) return 0;
	x_idx=i+1-start_idx;
	return 100*y_idx+x_idx;
}
void ExclusiveDijetsAnalysisUsingPPS::setToFParameters(std::string geo,double xpos)
{
	//     if (!fInitialized) setDataFormat();
	if (BeamX.count(1)==0||BeamX.count(2)==0||BeamX.count(3)==0||
			BeamY.count(1)==0||BeamY.count(2)==0||BeamY.count(3)==0||
			BeamX_RMS.count(1)==0||BeamX_RMS.count(2)==0||BeamX_RMS.count(3)==0||
			BeamY_RMS.count(1)==0||BeamY_RMS.count(2)==0||BeamY_RMS.count(3)==0) {
		std::cout << "ERROR: Beam parameters not given..." << std::endl;
		return;
	}
	ToFXPos=-(xpos*BeamX_RMS[3]+RPWindowThickness);
	ToFXPosSigma=xpos;
	setToFGeometry(geo);
}
void ExclusiveDijetsAnalysisUsingPPS::setBeamParameters(int pos, double x,double y,double rms_x,double rms_y)
{
	if (pos<1||pos>3) {
		cout << "ERROR: unknown detector position index: " << pos << endl;
		cout << "       It should be 1 for Tracker 1, 2 for Tracker 2 and 3 for ToF."<<endl;
		return;
	}
	BeamX[pos]=x;
	BeamY[pos]=y;
	BeamX_RMS[pos]=rms_x;
	BeamY_RMS[pos]=rms_y;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExclusiveDijetsAnalysisUsingPPS);
