import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
#'file:/user/amkalsi/MetFilterStudy/ECALDead/CMSSW_10_2_13/src/MetScanning/skim/python/crab_projects_Data_v1/16963797-0937-E811-ABE2-008CFAE45134.root'
'file:pickevents_DM_17Spet_BadMuon.root',
'file:pickevents_DM_PromptRecoD_BadMuon.root'
        )
                            )

process.TFileService = cms.Service("TFileService", fileName = cms.string("output.root") )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.GlobalTag.globaltag="102X_dataRun2_v8"
process.GlobalTag.globaltag="102X_upgrade2018_realistic_v15"

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")


process.load('Configuration.StandardSequences.Reconstruction_cff')


process.ntuplemakerminiaod = cms.EDAnalyzer('METScanningNtupleMakerMINIAOD',
                                            METFiltersPAT = cms.InputTag("TriggerResults::PAT"),
                                            METFiltersRECO = cms.InputTag("TriggerResults::RECO"),
                                            ECALBadCalibFilterUpdate=cms.InputTag("ecalBadCalibReducedMINIAOD2019Filter"),
                                            ECALLaserCorrFilterUpdate=cms.InputTag("ecalLaserCorrFilter"),
                                            ECALDeadCellBoundaryEnergyFilterUpdate=cms.InputTag("ecalDeadCellBoundaryEnergyFilterUpdate"),
                                            BadChargedCandidateFilterUpdate=cms.InputTag("BadChargedCandidateFilterUpdate"),
                                            BadPFMuonFilterUpdate =cms.InputTag("BadPFMuonFilterUpdate"),
                                            Vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                            Jets=cms.InputTag("slimmedJets"),
                                            PFCandCollection=cms.InputTag("packedPFCandidates"),
                                            PFMet=cms.InputTag("slimmedMETs"),
                                            PuppiMet=cms.InputTag("slimmedMETsPuppi"),
					    electronCollection=cms.InputTag("slimmedElectrons"),
                                            muonCollection=cms.InputTag("slimmedMuons"),
                                            file_INPUT=cms.string("list_2tev_4tev.txt")
                              )

#Rerunning the ecalbadcalibration filter
from RecoMET.METFilters.ecalBadCalibFilter_cfi import ecalBadCalibFilter

baddetEcallistnew2019 = cms.vuint32(
    [872439604,872422825,872420274,872423218,872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,872437052,872420649,872421950,872437185,
     872422564,872421566,872421695,872421955,872421567,872437184,872421951,872421694,
     872437056,872437057,872437313,872438182,872438951,872439990,872439864,872439609,
     872437181,872437182,872437053,872436794,872436667,872436536,872421541,872421413,
     872421414,872421031,872423083,872421439,872423224,872421438,872420397,872421566,
     872422589,872423096,872422717,872423214,872421415,872422311,872421926,872439469,
     872438567,872436659,872439731,872438311,872438078,872438438,872439601,872437951,
     872437950,872439729,872436792,872438183,872439468,872436663,872439728,872439727,
     872437694,872437823,872438845,872438973,872439354,872438566,872439733,872436530,
     872436655,872439600,872439730]
    )

process.ecalBadCalibReducedMINIAOD2019Filter = ecalBadCalibFilter.clone(
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallistnew2019,
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )

#Rerunning the laser correction filter
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
process.ecalLaserCorrFilter = cms.EDFilter(
    "EcalLaserCorrFilter",
    EBRecHitSource = cms.InputTag("reducedEgamma:reducedEBRecHits"),
    EERecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    EBLaserMIN     = cms.double(0.3),
    EELaserMIN     = cms.double(0.3),
    EBLaserMAX     = cms.double(5.0), #this was updated wrt default
    EELaserMAX     = cms.double(100.0), #this was updated wrt default
    EBEnegyMIN     = cms.double(10.0),
    EEEnegyMIN     = cms.double(10.0),
    taggingMode    = cms.bool(True), #updated wrt default
    Debug          = cms.bool(False)
    )

#Rerunning EcalDeadCellBoundaryEnergyFilter
from RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi import EcalDeadCellBoundaryEnergyFilter
process.ecalDeadCellBoundaryEnergyFilterUpdate=EcalDeadCellBoundaryEnergyFilter.clone(
    recHitsEB = cms.InputTag("reducedEgamma:reducedEBRecHits"),
    recHitsEE = cms.InputTag("reducedEgamma:reducedEERecHits"),
    cutBoundEnergyDeadCellsEE=cms.untracked.double(10),
    taggingMode    = cms.bool(True)
    )

#Rerunning BadChargedCandidateFilter
from RecoMET.METFilters.BadChargedCandidateFilter_cfi import BadChargedCandidateFilter 
process.BadChargedCandidateFilterUpdate=BadChargedCandidateFilter.clone(
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    taggingMode    = cms.bool(True)
)

from RecoMET.METFilters.BadPFMuonFilter_cfi import BadPFMuonFilter
process.BadPFMuonFilterUpdate=BadPFMuonFilter.clone(
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    minDxyBestTrack = cms.double(0.2),
    minDzBestTrack = cms.double(0.5), 
    taggingMode    = cms.bool(True)
)


import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
JSONfile = 'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
process.source.lumisToProcess.extend(myLumis)

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,era='2018-Prompt')


process.p = cms.Path(
    process.ecalBadCalibReducedMINIAOD2019Filter *
    process.ecalLaserCorrFilter *
    process.ecalDeadCellBoundaryEnergyFilterUpdate *
    process.BadChargedCandidateFilterUpdate *
    process.BadPFMuonFilterUpdate *
#    process.egammaPostRecoSeq *    
    process.ntuplemakerminiaod
    )

