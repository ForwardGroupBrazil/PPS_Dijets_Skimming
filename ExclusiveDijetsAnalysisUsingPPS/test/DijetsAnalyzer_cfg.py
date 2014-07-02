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
options.register('Run','MC_NO_OOT_PU',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: MC with or not PU")
options.parseArguments()


# Some Variables
MC_OOT_PU = False
MC_NO_OOT_NO_PU = False
MC_NO_OOT_PU = False
MC_NO_OOT_PU_DPE = False
MC_NO_OOT_PU_SDPPU = False
MC_NO_OOT_PU_SDMPU = False
jettag = "ak5PFJetsCHS" #ak5PFJets or ak5PFJetsCHS


# Setting Code Options
print("")
print("########################################")
print("Running CEP NTuple ProducerMC no Pile-up")
print("########################################")
print("")

if options.Run == "MC_OOT_PU":
  print("")
  print(">> Running: MC OOT and Pile-up")
  print("")
  MC_OOT_PU = True
  fileout = 'ttreeCEPdijetsOOT_PU.root'

elif options.Run == "MC_NO_OOT_NO_PU":
  print("")
  print(">> Running: MC No OOT and No Pile Up")
  print("")
  MC_NO_OOT_NO_PU = True
  fileout = 'ttreeCEPdijetsNoOOT_NoPU.root'

elif options.Run == "MC_NO_OOT_PU":
  print("")
  print(">> Running: MC No OOT and Pile Up")
  print("")
  MC_NO_OOT_PU = True
  fileout = 'ttreeCEPdijetsNoOOT_PU.root'

elif options.Run == "MC_NO_OOT_PU_DPE":
  print("")
  print(">> Running: MC No OOT and Pile Up")
  print("")
  MC_NO_OOT_PU_DPE = True
  fileout = 'ttreeCEPdijetsDPENoOOT_PU.root'

elif options.Run == "MC_NO_OOT_PU_SDPPU":
  print("")
  print(">> Running: MC No OOT and Pile Up")
  print("")
  MC_NO_OOT_PU_SDPPU = True
  fileout = 'ttreeCEPdijetsSDPPUNoOOT_PU_2.root'

elif options.Run == "MC_NO_OOT_PU_SDMPU":
  print("")
  print(">> Running: MC No OOT and Pile Up")
  print("")
  MC_NO_OOT_PU_SDMPU = True
  fileout = 'ttreeCEPdijetsSDMPUNoOOT_PU.root'
  
elif options.Run == "MC_NO_OOT_PU_SDPOMWIG":
  print("")
  print(">> Running: MC No OOT and Pile Up")
  print("")
  MC_NO_OOT_PU_SDMPU = True
  fileout = 'ttreeCEPdijetsSDPOMWIGNoOOT_PU.root'

else:
  print("")
  print("")
  raise RuntimeError, "Unknown option. EXIT! YOU NEED TO SETUP WITH ONE OF THE CORRECT OPTIONS."
  print("")


# CMSSW Default Configuration
process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(
        SkipEvent = cms.untracked.vstring('ProductNotFound'),
        wantSummary = cms.untracked.bool(True)
        )
# Openning files with or without Pile-Up
if MC_OOT_PU:
   process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:/storage2/polme/GG/GG_new2_1.root',
          'file:/storage2/polme/GG/GG_new2_2.root',
          'file:/storage2/polme/GG/GG_new2_3.root',
          'file:/storage2/polme/GG/GG_new2_4.root',
          'file:/storage2/polme/GG/GG_new2_5.root'
     )
   )

if MC_NO_OOT_NO_PU:
   process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:/storage2/polme/GG/GG_noOOT_noPU_1.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_10a.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_10b.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_2a.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_2b.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_3a.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_3b.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_4.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_5.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_6.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_7.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_8.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_9a.root',
          'file:/storage2/polme/GG/GG_noOOT_noPU_9b.root'
      )
   )
if MC_NO_OOT_PU_SDMPU:
   process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:/storage2/polme/PWSD/PomwigSDM_PU_1.root',
          'file:/storage2/polme/PWSD/PomwigSDM_PU_2.root',
          'file:/storage2/polme/PWSD/PomwigSDM_PU_3.root',
          'file:/storage2/polme/PWSD/PomwigSDM_PU_4.root',
          'file:/storage2/polme/PWSD/PomwigSDM_PU_5.root',
          #'file:/storage2/polme/PWSD/PomwigSDM_PU_6.root',
          'file:/storage2/polme/PWSD/PomwigSDM_PU_7.root',
          'file:/storage2/polme/PWSD/PomwigSDM_PU_8.root',
          'file:/storage2/polme/PWSD/PomwigSDM_PU_10.root'
      )
   )

if MC_NO_OOT_PU_SDPPU:
   process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:/storage2/polme/PWSD/PomwigSDP_PU_1.root',
          'file:/storage2/polme/PWSD/PomwigSDP_PU_2.root',
          'file:/storage2/polme/PWSD/PomwigSDP_PU_3.root',
          'file:/storage2/polme/PWSD/PomwigSDP_PU_4.root',
          'file:/storage2/polme/PWSD/PomwigSDP_PU_5.root',
          'file:/storage2/polme/PWSD/PomwigSDP_PU_6.root',
          'file:/storage2/polme/PWSD/PomwigSDP_PU_7.root',
          'file:/storage2/polme/PWSD/PomwigSDP_PU_8.root',
          'file:/storage2/polme/PWSD/PomwigSDP_PU_10.root'
    )
 )


if MC_NO_OOT_PU_POMWIGSD:
 process.source = cms.Source("PoolSource",
  duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
  fileNames = cms.untracked.vstring(
                'file:/storage2/polme/PWSD/PomwigSDM_PU_1.root',
                'file:/storage2/polme/PWSD/PomwigSDM_PU_2.root',
                'file:/storage2/polme/PWSD/PomwigSDM_PU_3.root',
                'file:/storage2/polme/PWSD/PomwigSDM_PU_4.root',
                'file:/storage2/polme/PWSD/PomwigSDM_PU_5.root',
                'file:/storage2/polme/PWSD/PomwigSDM_PU_7.root',
                'file:/storage2/polme/PWSD/PomwigSDM_PU_8.root',
                #'file:/storage2/polme/PWSD/PomwigSDM_PU_9.root',
                'file:/storage2/polme/PWSD/PomwigSDM_PU_10.root',
                'file:/storage2/polme/PWSD/PomwigSDP_PU_1.root',
                'file:/storage2/polme/PWSD/PomwigSDP_PU_2.root',
                'file:/storage2/polme/PWSD/PomwigSDP_PU_3.root',
                'file:/storage2/polme/PWSD/PomwigSDP_PU_4.root',
                'file:/storage2/polme/PWSD/PomwigSDP_PU_5.root',
                'file:/storage2/polme/PWSD/PomwigSDP_PU_6.root',
                'file:/storage2/polme/PWSD/PomwigSDP_PU_7.root',
                'file:/storage2/polme/PWSD/PomwigSDP_PU_8.root',
                'file:/storage2/polme/PWSD/PomwigSDP_PU_9.root',
                'file:/storage2/polme/PWSD/PomwigSDP_PU_10.root'
  )
 )



if MC_NO_OOT_PU_DPE:
   process.source = cms.Source("PoolSource",
      duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
      fileNames = cms.untracked.vstring(
          'file:/storage2/polme/PW/Pomwig_PU_1.root',
          'file:/storage2/polme/PW/Pomwig_PU_2.root',
          'file:/storage2/polme/PW/Pomwig_PU_3.root',
          'file:/storage2/polme/PW/Pomwig_PU_4.root',
          'file:/storage2/polme/PW/Pomwig_PU_5.root',
          'file:/storage2/polme/PW/Pomwig_PU_6.root',
          'file:/storage2/polme/PW/Pomwig_PU_7.root',
          'file:/storage2/polme/PW/Pomwig_PU_8.root',
          'file:/storage2/polme/PW/Pomwig_PU_10.root',
          'file:/storage2/polme/PW/Pomwig_PU_11.root',
          'file:/storage2/polme/PW/Pomwig_PU_12.root',
          'file:/storage2/polme/PW/Pomwig_PU_13.root',
          'file:/storage2/polme/PW/Pomwig_PU_14.root',
          'file:/storage2/polme/PW/Pomwig_PU_15.root',
          'file:/storage2/polme/PW/Pomwig_PU_16.root',
          'file:/storage2/polme/PW/Pomwig_PU_17.root',
          'file:/storage2/polme/PW/Pomwig_PU_18.root',
          'file:/storage2/polme/PW/Pomwig_PU_19.root',
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
process.demo = cms.EDAnalyzer('ExclusiveDijetsAnalysisUsingPPSv14',
                               MakePlots = cms.bool(True),
                               PPS_Flag = cms.bool(False),
                               RunWithWeightGen = cms.bool(False),
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

# Path, Run modules in order.
process.p = cms.Path(process.goodOfflinePrimaryVertices*process.ak5JetTracksAssociatorAtVertex*process.demo)

