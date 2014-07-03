'''
>>--------------------<<
CEP JJ NTuple Producer
>>--------------------<<

Goal:
Produce CEP JJ ntuple.

Usage:
cmsRun DijetsAnalyzer_cfg.py

Example:
cmsRun DijetsAnalyzer_cfg.py Run=MC_NO_PU

Optional arguments:
Run = MC_OOT_PU, MC_NO_OOT_PU or MC_NO_OOT_NO_PU

Authors: Brazilian/PPS Group
'''

# Loading Modules and Python Libraries
import FWCore.ParameterSet.Config as cms
import os, sys
import atexit
from FWCore.ParameterSet.VarParsing import VarParsing

# Setting Input Parameters from Line Command
options = VarParsing ('analysis')
options.register('Run','MC_NO_OOT_PU_QCD',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: MC with or not PU")
options.parseArguments()

# Some Variables
MC_NO_OOT_PU = False
MC_NO_OOT_PU_POMWIG = False
MC_NO_OOT_PU_QCD = False
jettag = "ak5PFJetsCHS" #ak5PFJets or ak5PFJetsCHS

# Setting Code Options
print("")
print("########################################")
print("Running CEP NTuple ProducerMC no Pile-up")
print("########################################")
print("")

if options.Run == "MC_NO_OOT_PU_QCD":
 print("")
 print(">> Running: MC No OOT and Pile Up - QCD")
 print("")
 MC_NO_OOT_PU_QCD = True
 fileout = 'test_ttree_QCD_w_t_PPSInfo.root'

elif options.Run == "MC_NO_OOT_PU_POMWIG":
 print("")
 print(">> Running: MC No OOT and Pile Up - POMWIG")
 print("")
 MC_NO_OOT_PU_POMWIG = True
 fileout = 'test_ttree_pomwig_NoOOT_PU.root'

elif options.Run == "MC_NO_OOT_PU":
 print("")
 print(">> Running: MC No OOT and Pile Up - ExHume")
 print("")
 MC_NO_OOT_PU = True
 fileout = 'ttree_exhume_NoOOT_PU.root'

else:
  print("")
  print("")
  raise RuntimeError, "Unknown option. EXIT! YOU NEED TO SETUP WITH ONE OF THE CORRECT OPTIONS."
  print("")


# CMSSW Default Configuration
process = cms.Process("Demo")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Openning files with or without Pile-Up

if MC_NO_OOT_PU_QCD:
 process.source = cms.Source("PoolSource",
  duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
  fileNames = cms.untracked.vstring(
                 'file:/storage2/sfonseca/CMSSW/PPS_Analysis/git3/CMSSW_6_2_5/src/PPS_Dijets_Skimming/ExclusiveDijetsAnalysisUsingPPS/QCD_PU40_BX25.root'

	 )
 )




if MC_NO_OOT_PU_POMWIG:
 process.source = cms.Source("PoolSource",
  duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
  fileNames = cms.untracked.vstring(
		'file:/storage2/polme/PW/Pomwig_PU_1.root',
		'file:/storage2/polme/PW/Pomwig_PU_2.root',
		'file:/storage2/polme/PW/Pomwig_PU_3.root',
		'file:/storage2/polme/PW/Pomwig_PU_4.root',
		'file:/storage2/polme/PW/Pomwig_PU_5.root',
		'file:/storage2/polme/PW/Pomwig_PU_6.root',
		'file:/storage2/polme/PW/Pomwig_PU_9.root',
		'file:/storage2/polme/PW/Pomwig_PU_10.root',
		'file:/storage2/polme/PW/Pomwig_PU_11.root',
		'file:/storage2/polme/PW/Pomwig_PU_12.root',
		'file:/storage2/polme/PW/Pomwig_PU_13.root',
		'file:/storage2/polme/PW/Pomwig_PU_14.root',
		'file:/storage2/polme/PW/Pomwig_PU_15.root',
		'file:/storage2/polme/PW/Pomwig_PU_16.root',
		'file:/storage2/polme/PW/Pomwig_PU_17.root',
		'file:/storage2/polme/PW/Pomwig_PU_18.root',
		'file:/storage2/polme/PW/Pomwig_PU_20.root'
  )
 )

if MC_NO_OOT_PU:
 process.source = cms.Source("PoolSource",
  duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
  fileNames = cms.untracked.vstring(
		'file:/storage2/polme/GG/GG_PU_1.root',
		'file:/storage2/polme/GG/GG_PU_2.root',
		'file:/storage2/polme/GG/GG_PU_3.root',
		'file:/storage2/polme/GG/GG_PU_4.root',
 	 	'file:/storage2/polme/GG/GG_PU_5.root',
		'file:/storage2/polme/GG/GG_PU_6.root',
		'file:/storage2/polme/GG/GG_PU_7.root',
		'file:/storage2/polme/GG/GG_PU_8.root',
		'file:/storage2/polme/GG/GG_PU_9.root',
		'file:/storage2/polme/GG/GG_PU_10.root',
		'file:/storage2/polme/GG/GG_PU_11.root',
		'file:/storage2/polme/GG/GG_PU_12.root',
		'file:/storage2/polme/GG/GG_PU_13.root',
		'file:/storage2/polme/GG/GG_PU_14.root',
		'file:/storage2/polme/GG/GG_PU_15.root',
		'file:/storage2/polme/GG/GG_PU_16.root',
		'file:/storage2/polme/GG/GG_PU_17.root',
		'file:/storage2/polme/GG/GG_PU_18.root',
		'file:/storage2/polme/GG/GG_PU_19.root',
		'file:/storage2/polme/GG/GG_PU_20.root',
		'file:/storage2/polme/GG/GG_PU_21.root',
		'file:/storage2/polme/GG/GG_PU_22.root',
		'file:/storage2/polme/GG/GG_PU_23.root',
		'file:/storage2/polme/GG/GG_PU_24.root',
		'file:/storage2/polme/GG/GG_PU_25.root',
		'file:/storage2/polme/GG/GG_PU_26.root',
		'file:/storage2/polme/GG/GG_PU_27.root',
		'file:/storage2/polme/GG/GG_PU_28.root',
		'file:/storage2/polme/GG/GG_PU_29.root',
		'file:/storage2/polme/GG/GG_PU_30.root',
		'file:/storage2/polme/GG/GG_PU_31.root',
		'file:/storage2/polme/GG/GG_PU_32.root',
		'file:/storage2/polme/GG/GG_PU_33.root',
		'file:/storage2/polme/GG/GG_PU_34.root',
		'file:/storage2/polme/GG/GG_PU_35.root',
		'file:/storage2/polme/GG/GG_PU_36.root',
		'file:/storage2/polme/GG/GG_PU_37.root',
		'file:/storage2/polme/GG/GG_PU_38.root',
		'file:/storage2/polme/GG/GG_PU_39.root',
		'file:/storage2/polme/GG/GG_PU_40.root',
		'file:/storage2/polme/GG/GG_PU_41.root',
		'file:/storage2/polme/GG/GG_PU_42.root',
		'file:/storage2/polme/GG/GG_PU_43.root',
		'file:/storage2/polme/GG/GG_PU_44.root',
		'file:/storage2/polme/GG/GG_PU_45.root',
		'file:/storage2/polme/GG/GG_PU_46.root',
		'file:/storage2/polme/GG/GG_PU_47.root',
		'file:/storage2/polme/GG/GG_PU_48.root',
		'file:/storage2/polme/GG/GG_PU_49.root',
		'file:/storage2/polme/GG/GG_PU_50.root'


	   )
   )
#  process.source.dropDescendantsOfDroppedBranches = cms.untracked.bool(True)

process.ak5PFJetsCHS = cms.EDFilter(
    "PFJetSelector",
    src = cms.InputTag("ak5PFJetsCHS"),
    cut = cms.string("pt > 100"),
    filter = cms.bool(True)
)

#process.filter = cms.EDFilter("NJetsMC",
#      GenTag = cms.untracked.InputTag('ak5GenJets'),
#      MinPt  = cms.double(100.),
#      Njets  = cms.int32(1)
#   )


# Adding Good Primary Vertex
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

# Adding Tracks Associator with Vertex Collection
process.ak5JetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
       tracks = cms.InputTag("generalTracks"),
                                                         jets = cms.InputTag("ak5PFJets"),
                                                         coneSize = cms.double(0.5)
                                                       )

# EDAnalyzer Parameters
process.demo = cms.EDAnalyzer('ExclusiveDijetsAnalysisUsingPPS_QCD',
                               MakePlots = cms.bool(True),
			       PPS_Flag = cms.bool(False),	
                               RunWithWeightGen = cms.bool(True),
                               GenJets = cms.InputTag('ak5GenJets'),#Gen information
                               JetTag = cms.InputTag(jettag),
ParticleFlowTag = cms.InputTag("particleFlow"),
                               VertexTag = cms.InputTag("goodOfflinePrimaryVertices"),
PPSTag = cms.untracked.string("PPSReco"),
                               pTPFThresholdCharged = cms.double(0.1),
                               energyPFThresholdBar = cms.double(1.5),
                               energyPFThresholdEnd = cms.double(3.5),
                               energyPFThresholdHF = cms.double(4.0),
                               cmsVertexResolution = cms.double(0.005), #cm
                               PPSVertexResolution = cms.double(0.2), #cm
                               EBeam = cms.double(13000.)
                               )

# TFileService
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(fileout)
                                   )

####################################################################################################


#
#process.f = cms.Path(process.filter)
#
# Path, Run modules in order.
#process.p = cms.Path(process.ak5PFJetsCHS*process.goodOfflinePrimaryVertices*process.ak5JetTracksAssociatorAtVertex*process.demo)
process.p = cms.Path(process.ak5PFJetsCHS+process.goodOfflinePrimaryVertices*process.ak5JetTracksAssociatorAtVertex*process.demo)



