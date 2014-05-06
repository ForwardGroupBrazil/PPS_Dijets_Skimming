import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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

process.ak5JetTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
      	                                                 tracks = cms.InputTag("generalTracks"),
                                                         jets   = cms.InputTag("ak5PFJets"),
                                                         coneSize = cms.double(0.5)
                                                       )


process.demo = cms.EDAnalyzer('ExclusiveDijetsAnalysisUsingPPS',
			       JetTag = cms.InputTag("ak5PFJets"),
			       ParticleFlowTag = cms.InputTag("particleFlow"),
			       PPSTag = cms.untracked.string("PPSReco"),
                               pTPFThresholdCharged = cms.double(0.1),
                               energyPFThresholdBar = cms.double(1.5),
                               energyPFThresholdEnd = cms.double(3.5),
                               energyPFThresholdHF = cms.double(4.0)
                               )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('ttreeCEPdijets.root')
                                   )

process.p = cms.Path(process.ak5JetTracksAssociatorAtVertex*process.demo)
