import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
          'file:/afs/cern.ch/work/p/polme/public/PPS/CMSSW_6_2_0/src/GG_newtest_1.root'
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
                                       fileName = cms.string('histo_dijets.root')
                                   )

process.p = cms.Path(process.ak5JetTracksAssociatorAtVertex*process.demo)
