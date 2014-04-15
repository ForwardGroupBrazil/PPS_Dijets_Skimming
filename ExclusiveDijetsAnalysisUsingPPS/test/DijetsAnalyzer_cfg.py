import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/storage2/sfonseca/CMSSW/PPS_Analysis/CMSSW_6_2_0/src/PPS_Dijets_Skimming/patTuple_PF2PAT.root'
#    'file:/afs/cern.ch/work/p/polme/public/PPS/CMSSW_6_2_0/src/test_GG_exhume.root'
    )
)

process.demo = cms.EDAnalyzer('ExclusiveDijetsAnalysisUsingPPS',
			       ## jet collections ###########################
			       JetTag = cms.InputTag("ak5PFJets"),
			       ParticleFlowTag = cms.InputTag("particleFlow"),
			       PPSTag = cms.untracked.string("PPSReco")			
                               )




process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('histo_dijets.root')
                                   )



process.p = cms.Path(process.demo)
