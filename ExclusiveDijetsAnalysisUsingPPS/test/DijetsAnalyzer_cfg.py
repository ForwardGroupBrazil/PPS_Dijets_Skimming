#flake8: noqa

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
Run = MC_PU or MC_NO_PU

Authors: Brazilian/PPS Group
'''

# Loading Modules and Python Libraries
import FWCore.ParameterSet.Config as cms
import os, sys
import atexit
from FWCore.ParameterSet.VarParsing import VarParsing

# Setting Input Parameters from Line Command
options = VarParsing ('analysis')
options.register('Run','MC_PU',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: MC with or not PU")
options.parseArguments()

# Some Variables
MCNOPU = False
MCPU = False

# Setting Code Options
print("")
print("########################################")
print("Running CEP NTuple ProducerMC no Pile-up")
print("########################################")
print("")

if options.Run == "MC_NO_PU":
  print("")
  print(">> Running: MC no Pile-up")
  print("")
  MCNOPU = True
  fileout = 'ttreeCEPdijetsNoPU.root'

elif options.Run == "MC_PU":
  print("")
  print(">> Running: MC with Pile Up")
  print("")
  MCPU = True
  fileout = 'ttreeCEPdijetsPU.root'

else:
  print("")
  print("")
  raise RuntimeError, "Unknown option. EXIT! YOU NEED TO SETUP WITH ONE OF THE CORRECT OPTIONS."
  print("")


# CMSSW Default Configuration
process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Openning files with or without Pile-Up
if MCPU:
   process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:///storage1/dmf/Samples/PPS/GG_new2_1.root',
          'file:///storage1/dmf/Samples/PPS/GG_new2_2.root',
          'file:///storage1/dmf/Samples/PPS/GG_new2_3.root',
          'file:///storage1/dmf/Samples/PPS/GG_new2_4.root',
          'file:///storage1/dmf/Samples/PPS/GG_new2_5.root'
     )
   )

if MCNOPU:
   process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:/afs/cern.ch/work/p/polme/public/PPS/GG/GG_noOOT_noPU_1.root',
          'file:/afs/cern.ch/work/p/polme/public/PPS/GG/GG_noOOT_noPU_2.root',
          'file:/afs/cern.ch/work/p/polme/public/PPS/GG/GG_noOOT_noPU_3.root',
          'file:/afs/cern.ch/work/p/polme/public/PPS/GG/GG_noOOT_noPU_4.root'
      )
   )

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
                                                         jets   = cms.InputTag("ak5PFJets"),
                                                         coneSize = cms.double(0.5)
                                                       )

# EDAnalyzer Parameters
process.demo = cms.EDAnalyzer('ExclusiveDijetsAnalysisUsingPPS',
                               MakePlots = cms.bool(True),
			       JetTag = cms.InputTag("ak5PFJets"),
                               #JetTag = cms.InputTag("AK5PFchs"), #Testing
			       ParticleFlowTag = cms.InputTag("particleFlow"),
                               VertexTag = cms.InputTag("goodOfflinePrimaryVertices"),
			       PPSTag = cms.untracked.string("PPSReco"),
                               pTPFThresholdCharged = cms.double(0.1),
                               energyPFThresholdBar = cms.double(1.5),
                               energyPFThresholdEnd = cms.double(3.5),
                               energyPFThresholdHF = cms.double(4.0),
                               cmsVertexResolution = cms.double(0.2), #cm
                               PPSVertexResolution = cms.double(0.2) #cm
                               )

# TFileService
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(fileout)
                                   )

# Path, Run modules in order.
process.p = cms.Path(process.goodOfflinePrimaryVertices*process.ak5JetTracksAssociatorAtVertex*process.demo)
