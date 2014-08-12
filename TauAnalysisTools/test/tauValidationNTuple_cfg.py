#import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.coreTools import *

#runOnMC = False
runOnMC = True
if runOnMC:
    print "Running on MC"
else:
    print "Running on Data"


process = cms.Process("tauNTuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.extend(['TauValidationNTupleProd'])
process.MessageLogger.cerr.default.limit = -1
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

#-- Calibration tag -----------------------------------------------------------

process.load("Configuration.Geometry.GeometryExtended2019Reco_cff")
process.load("Configuration.Geometry.GeometryExtended2019_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'DES19_62_V8::All'
process.load("Configuration.StandardSequences.MagneticField_38T_PostLS1_cff")

#-- PAT standard config -------------------------------------------------------

process.options.allowUnscheduled = cms.untracked.bool( True )

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#process.load("RecoVertex.Configuration.RecoVertex_cff")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

## Dummy output for PAT. Not used in the analysis ##
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName       = cms.untracked.string('dummy.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    dropMetaData   = cms.untracked.string('DROPPED'),
    outputCommands = cms.untracked.vstring('keep *')
    )

process.source = cms.Source(
    'PoolSource',
    skipEvents = cms.untracked.uint32( 0 ),
    fileNames = cms.untracked.vstring(
        #'file:////disk1/knutzen//TauPOG/CrossCheck_Upgrade/testRun/ACABF955-DAF7-E311-A024-0025905A608E.root'
        #'file:////disk1/knutzen//TauPOG/CrossCheck_Upgrade/testRun/CE034383-D0F7-E311-A247-002618943983.root'
        #'/store/mc/Muon2023Upg14DR/DYToMuMu_M-20_TuneZ2star_14TeV-pythia6-tauola/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/001DC0AD-54E5-E311-87E5-0025905A60E4.root'
        #'/store/mc/GEM2019Upg14DR/DYToTauTau_M-20_TuneZ2star_14TeV-pythia6-tauola/AODSIM/age1k_PU140bx25_PH1_1K_FB_V1-v1/00000/001601DE-58E5-E311-A498-002618FDA28E.root'
        #'/store/mc/GEM2019Upg14DR/DYToTauTau_M-20_TuneZ2star_14TeV-pythia6-tauola/AODSIM/PU50bx25_DES19_62_V8-v1/00000/00D0FBCE-4BE3-E311-A444-00261894386B.root'
        #'/store/mc/GEM2019Upg14DR/QCD_Pt-15to3000_Tune4C_Flat_14TeV_pythia8/AODSIM/NoPileUp_DES19_62_V8-v1/00000/CEB07121-01F8-E311-BF5F-0025905B85EE.root'
        #'/store/mc/GEM2019Upg14DR/QCD_Pt-15to3000_Tune4C_Flat_14TeV_pythia8/AODSIM/NoPileUp_DES19_62_V8-v1/00000/CEB07121-01F8-E311-BF5F-0025905B85EE.root'
        '/store/mc/GEM2019Upg14DR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/AODSIM/final_phase1_PU50bx25_DES19_62_V8-v1/00000/06C9D927-741B-E411-8A09-0025905A609A.root'
        )
    )

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )

#process.load("PhysicsTools/PatAlgos/patSequences_cff")
#
#from PhysicsTools.PatAlgos.tools.tauTools import *
#switchToPFTauHPS(process) #create HPS Taus from the pat default sequence


  # --------------------Modifications for taus--------------------
process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)
  
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

# --------------------Modifications for MET--------------------
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1')        
    )

# switch on PAT trigger
#from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
#switchOnTrigger(process) #create pat trigger objects

#process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

#removeMCMatching(process, ['All'])

from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools import pfTools

jetCorrFactors = cms.vstring( 'L1FastJet', 'L2Relative', 'L3Absolute' )

pfTools.usePF2PAT( process,
                   runPF2PAT = True,
                   jetAlgo = 'AK5',
                   jetCorrections = ( 'AK5PFchs', jetCorrFactors, "" ),
                   #jetCorrections = ['L1FastJet', 'L2Relative', 'L3Absolute'],
                   runOnMC = runOnMC,
                   postfix = "PFlow",
                   pvCollection = cms.InputTag( 'goodOfflinePrimaryVertices' ),
                   outputModules = []
                   )

#switchJetCollection(
#        process,
#        cms.InputTag('ak5PFJets')
#        #doJTA = True
#        #doBTagging = True
#        #jetCorrLabel = _jetCorrections,
#        #doType1MET = False,
#        #doJetID = True,
#        #jetIdLabel = "ak5"
#)

process.TFileService = cms.Service("TFileService", fileName = cms.string('TauValidationNTuple.root') #output file
                                                                                  )
###############################################
#############    User input  ##################
###############################################

# Enter a list with all trigger filter names you want to investigate.
# A bool with the same name will be created for each filter which denotes if a filter object is matched to the tag tau
filterName = [
               ]

#Enter a list of HPS discriminators you want to store in the output tree for the tag tau
IDName = [  
               "decayModeFindingNewDMs", 
               "decayModeFindingOldDMs", 
               "decayModeFinding", 
               "byLooseIsolation", 
               "byVLooseCombinedIsolationDeltaBetaCorr", 
               "byLooseCombinedIsolationDeltaBetaCorr", 
               "byMediumCombinedIsolationDeltaBetaCorr", 
               "byTightCombinedIsolationDeltaBetaCorr", 
               "byCombinedIsolationDeltaBetaCorrRaw", 
               "byLooseCombinedIsolationDeltaBetaCorr3Hits", 
               "byMediumCombinedIsolationDeltaBetaCorr3Hits", 
               "byTightCombinedIsolationDeltaBetaCorr3Hits", 
               "byCombinedIsolationDeltaBetaCorrRaw3Hits", 
               "chargedIsoPtSum", 
               "neutralIsoPtSum", 
               "puCorrPtSum", 
               "rhoCorrPtSum", 
               "neutralHadronIsoPtSum", 
               "weightedNeutralIsoPtSum1", 
               "weightedNeutralIsoPtSum2", 
               "weightedNeutralIsoPtSum1NQ", 
               "weightedNeutralIsoPtSum2NQ", 
               "weightedNeutralHadronIsoPtSum1", 
               "weightedNeutralHadronIsoPtSum2", 
               "weightedNeutralHadronIsoPtSum1NQ", 
               "weightedNeutralHadronIsoPtSum2NQ", 
               "byIsolationMVA3oldDMwoLTraw", 
               "byVLooseIsolationMVA3oldDMwoLT", 
               "byLooseIsolationMVA3oldDMwoLT", 
               "byMediumIsolationMVA3oldDMwoLT", 
               "byTightIsolationMVA3oldDMwoLT", 
               "byVTightIsolationMVA3oldDMwoLT", 
               "byVVTightIsolationMVA3oldDMwoLT", 
               "byIsolationMVA3oldDMwLTraw", 
               "byVLooseIsolationMVA3oldDMwLT", 
               "byLooseIsolationMVA3oldDMwLT", 
               "byMediumIsolationMVA3oldDMwLT", 
               "byTightIsolationMVA3oldDMwLT", 
               "byVTightIsolationMVA3oldDMwLT", 
               "byVVTightIsolationMVA3oldDMwLT", 
               "byIsolationMVA3newDMwoLTraw", 
               "byVLooseIsolationMVA3newDMwoLT", 
               "byLooseIsolationMVA3newDMwoLT", 
               "byMediumIsolationMVA3newDMwoLT", 
               "byTightIsolationMVA3newDMwoLT", 
               "byVTightIsolationMVA3newDMwoLT", 
               "byVVTightIsolationMVA3newDMwoLT", 
               "byIsolationMVA3newDMwLTraw", 
               "byVLooseIsolationMVA3newDMwLT", 
               "byLooseIsolationMVA3newDMwLT", 
               "byMediumIsolationMVA3newDMwLT", 
               "byTightIsolationMVA3newDMwLT", 
               "byVTightIsolationMVA3newDMwLT", 
               "byVVTightIsolationMVA3newDMwLT", 
               "againstElectronLoose", 
               "againstElectronMedium", 
               "againstElectronTight", 
               "againstElectronMedium2", 
               "againstElectronMVA5raw", 
               "againstElectronMVA5category", 
               "againstElectronVLooseMVA5", 
               "againstElectronLooseMVA5", 
               "againstElectronMediumMVA5", 
               "againstElectronTightMVA5", 
               "againstElectronVTightMVA5", 
               "againstElectronDeadECAL", 
               "againstMuonLoose", 
               "againstMuonMedium", 
               "againstMuonTight", 
               "againstMuonLoose2", 
               "againstMuonMedium2", 
               "againstMuonTight2", 
               "againstMuonLoose3", 
               "againstMuonTight3", 
               "againstMuonMVAraw", 
               "againstMuonLooseMVA", 
               "againstMuonMediumMVA", 
               "againstMuonTightMVA"
               ]

common_ntuple_branches = cms.PSet(
    index = cms.string("index"), # Index of reco object in the event
    nRecoObjects = cms.string("nTotalObjects"), # Number of reco objects in the event

    nvtx = cms.string("Nvtx"),

    PV_x = cms.string("getPV.x"),
    PV_y = cms.string("getPV.y"),
    PV_z = cms.string("getPV.z"),
    PV_xError = cms.string("getPV.xError"),
    PV_yError = cms.string("getPV.yError"),
    PV_zError = cms.string("getPV.zError"),
    PV_ndof = cms.string("getPV.ndof"),
    PV_nTracks = cms.string("getPV.nTracks"),
    PV_isValid = cms.string("getPV.isValid"),



    RunNr = cms.string("RunNr"),
    EvtNr = cms.string("EvtNr"),
    LumiSec = cms.string("LumiSec"),

    hasValidEleMVA = cms.string("hasValidEleMVA"),

    EleMVA_Tau_GSFTracklnPt = cms.string("Tau_GSFTracklnPt"),
    EleMVA_Tau_GSFTrackEta = cms.string("Tau_GSFTrackEta"),
    EleMVA_Tau_VisMass = cms.string("Tau_VisMass"),
    EleMVA_Tau_EmFraction = cms.string("Tau_EmFraction"),
    EleMVA_Tau_GSFTrackResol = cms.string("Tau_GSFTrackResol"),
    EleMVA_Tau_Pt = cms.string("Tau_Pt"),
    EleMVA_Tau_GammaEnFrac = cms.string("Tau_GammaEnFrac"),
    EleMVA_Tau_LeadChargedPFCandEtaAtEcalEntrance = cms.string("Tau_LeadChargedPFCandEtaAtEcalEntrance"),
    EleMVA_Tau_LeadChargedPFCandPt = cms.string("Tau_LeadChargedPFCandPt"),
    EleMVA_Tau_NumGammaCands = cms.string("Tau_NumGammaCands"),

    EleMVA_Tau_HadrHoP = cms.string("Tau_HadrHoP"),
    EleMVA_Tau_GammaEtaMom = cms.string("Tau_GammaEtaMom"),
    EleMVA_Tau_GammaPhiMom = cms.string("Tau_GammaPhiMom"),
    EleMVA_Tau_HadrEoP = cms.string("Tau_HadrEoP"),
    EleMVA_Tau_EtaAtEcalEntrance = cms.string("Tau_EtaAtEcalEntrance"),
    EleMVA_Tau_GSFChi2 = cms.string("Tau_GSFChi2"),
    EleMVA_Tau_GSFNumHits = cms.string("Tau_GSFNumHits"),
    EleMVA_Tau_KFNumHits = cms.string("Tau_KFNumHits"),
    EleMVA_Tau_HadrMva = cms.string("Tau_HadrMva"),
    EleMVA_Tau_HasGsf = cms.string("Tau_HasGsf"),

    EleMVA_Elec_GSFTrackEta = cms.string("Elec_GSFTrackEta"),
    EleMVA_Elec_EtotOverPin = cms.string("Elec_EtotOverPin"),
    EleMVA_Elec_Fbrem = cms.string("Elec_Fbrem"),
    EleMVA_Elec_GSFNumHits = cms.string("Elec_GSFNumHits"),
    EleMVA_Elec_Chi2GSF = cms.string("Elec_Chi2GSF"),
    EleMVA_Elec_GSFTracklnPt = cms.string("Elec_GSFTracklnPt"),
    EleMVA_Elec_EgammaOverPdif = cms.string("Elec_EgammaOverPdif"),
    EleMVA_Elec_GSFTrackResol = cms.string("Elec_GSFTrackResol"),

    recoTauCand_Pt = cms.string("recoTauCand.pt"),
    recoTauCand_Eta = cms.string("recoTauCand.eta"),
    recoTauCand_Phi = cms.string("recoTauCand.phi"),
    recoTauCand_DecayMode = cms.string("recoTauCand.decayMode"),
    recoTauCand_Vz = cms.string("recoTauCand.vz"),
    recoTauCand_dxy_PCA_X = cms.string("recoTauCand.dxy_PCA.X"), 
    recoTauCand_dxy_PCA_Y = cms.string("recoTauCand.dxy_PCA.Y"), 
    recoTauCand_dxy_PCA_Z = cms.string("recoTauCand.dxy_PCA.Z"), 
    recoTauCand_dxy = cms.string("recoTauCand.dxy"),
    recoTauCand_dxy_error = cms.string("recoTauCand.dxy_error"),
    recoTauCand_dxy_Sig = cms.string("recoTauCand.dxy_Sig"),
    recoTauCand_hasSecVtx = cms.string("recoTauCand.hasSecondaryVertex"),
    recoTauCand_flightLength_X = cms.string("recoTauCand.flightLength.X"),
    recoTauCand_flightLength_Y = cms.string("recoTauCand.flightLength.Y"),
    recoTauCand_flightLength_Z = cms.string("recoTauCand.flightLength.Z"),
    recoTauCand_flightLength_Cov_00 = cms.string("recoTauCand.flightLengthCov.At(0,0)"),
    recoTauCand_flightLength_Cov_10 = cms.string("recoTauCand.flightLengthCov.At(1,0)"),
    recoTauCand_flightLength_Cov_20 = cms.string("recoTauCand.flightLengthCov.At(2,0)"),
    recoTauCand_flightLength_Cov_01 = cms.string("recoTauCand.flightLengthCov.At(0,1)"),
    recoTauCand_flightLength_Cov_11 = cms.string("recoTauCand.flightLengthCov.At(1,1)"),
    recoTauCand_flightLength_Cov_21 = cms.string("recoTauCand.flightLengthCov.At(2,1)"),
    recoTauCand_flightLength_Cov_02 = cms.string("recoTauCand.flightLengthCov.At(0,2)"),
    recoTauCand_flightLength_Cov_12 = cms.string("recoTauCand.flightLengthCov.At(1,2)"),
    recoTauCand_flightLength_Cov_22 = cms.string("recoTauCand.flightLengthCov.At(2,2)"),
    recoTauCand_flightLengthSig = cms.string("recoTauCand.flightLengthSig"),
    recoTauCand_primaryVertex_X = cms.string("recoTauCand.primaryVertexPos().X"),
    recoTauCand_primaryVertex_Y = cms.string("recoTauCand.primaryVertexPos().Y"),
    recoTauCand_primaryVertex_Z = cms.string("recoTauCand.primaryVertexPos().Z"),
    recoTauCand_primaryVertexCov_00 = cms.string("recoTauCand.primaryVertexCov.At(0,0)"),
    recoTauCand_primaryVertexCov_10 = cms.string("recoTauCand.primaryVertexCov.At(1,0)"),
    recoTauCand_primaryVertexCov_20 = cms.string("recoTauCand.primaryVertexCov.At(2,0)"),
    recoTauCand_primaryVertexCov_01 = cms.string("recoTauCand.primaryVertexCov.At(0,1)"),
    recoTauCand_primaryVertexCov_11 = cms.string("recoTauCand.primaryVertexCov.At(1,1)"),
    recoTauCand_primaryVertexCov_21 = cms.string("recoTauCand.primaryVertexCov.At(2,1)"),
    recoTauCand_primaryVertexCov_02 = cms.string("recoTauCand.primaryVertexCov.At(0,2)"),
    recoTauCand_primaryVertexCov_12 = cms.string("recoTauCand.primaryVertexCov.At(1,2)"),
    recoTauCand_primaryVertexCov_22 = cms.string("recoTauCand.primaryVertexCov.At(2,2)"),
    recoTauCand_secVertex_X = cms.string("recoTauCand.secondaryVertexPos().X"),
    recoTauCand_secVertex_Y = cms.string("recoTauCand.secondaryVertexPos().Y"),
    recoTauCand_secVertex_Z = cms.string("recoTauCand.secondaryVertexPos().Z"),
    recoTauCand_secVertexCov_00 = cms.string("recoTauCand.secondaryVertexCov.At(0,0)"),
    recoTauCand_secVertexCov_10 = cms.string("recoTauCand.secondaryVertexCov.At(1,0)"),
    recoTauCand_secVertexCov_20 = cms.string("recoTauCand.secondaryVertexCov.At(2,0)"),
    recoTauCand_secVertexCov_01 = cms.string("recoTauCand.secondaryVertexCov.At(0,1)"),
    recoTauCand_secVertexCov_11 = cms.string("recoTauCand.secondaryVertexCov.At(1,1)"),
    recoTauCand_secVertexCov_21 = cms.string("recoTauCand.secondaryVertexCov.At(2,1)"),
    recoTauCand_secVertexCov_02 = cms.string("recoTauCand.secondaryVertexCov.At(0,2)"),
    recoTauCand_secVertexCov_12 = cms.string("recoTauCand.secondaryVertexCov.At(1,2)"),
    recoTauCand_secVertexCov_22 = cms.string("recoTauCand.secondaryVertexCov.At(2,2)"),

    pfCand_dz = cms.string("pfCand_dz"),

    hasValidTrack = cms.string("hasValidTrack"),
    TransImpPara = cms.string("TransImpPara"),
    TransImpParaError = cms.string("TransImpParaError"),

    isAltTauObjMatched = cms.string("isAltTauObjMatched"),
    altTauObj_Pt = cms.string("altTauObj.pt"),
    altTauObj_Eta = cms.string("altTauObj.eta"),
    altTauObj_Phi = cms.string("altTauObj.phi"),

    tauSeedPfJet_Pt = cms.string("recoTauCand.pfJetRef().pt"),
    tauSeedPfJet_Eta = cms.string("recoTauCand.pfJetRef().eta"),
    tauSeedPfJet_Phi = cms.string("recoTauCand.pfJetRef().phi"),
    tauSeedPfJet_Vz = cms.string("recoTauCand.pfJetRef().vz"),

    pfJet_Pt = cms.string("PfJet.pt"),
    pfJet_Eta = cms.string("PfJet.eta"),
    pfJet_Phi = cms.string("PfJet.phi"),
    pfJetVz = cms.string("PfJet.vz"),
    isPfJetMatched = cms.string("isPfJetMatched"),

    genJet_Pt = cms.string("GenJet.pt"),
    genJet_Eta = cms.string("GenJet.eta"),
    genJet_Phi = cms.string("GenJet.phi"),
    genJetVz = cms.string("GenJet.vz"),
    isGenJetMatched = cms.string("isGenJetMatched"),
)
if runOnMC:
    common_ntuple_branches.genWeight = cms.string("genInfo.weight")

    common_ntuple_branches.isGenParticleMatched = cms.string("isGenParticelMatched")
    # Careful! Only use GenTauMatch (returns the values of the generator particle matched to the tagTau) values if "bool TauTrigMatch::GenTauMatchTest()" returns "true". Otherwise it contains (unrealsitic) default values
    common_ntuple_branches.GenParticle_Pt = cms.string("GenParticle.pt")
    common_ntuple_branches.GenParticle_Eta = cms.string("GenParticle.eta")
    common_ntuple_branches.GenParticle_Phi = cms.string("GenParticle.phi")
    common_ntuple_branches.GenParticel_pdgId = cms.string("GenParticle.pdgId")
    common_ntuple_branches.isTauGenJetMatched = cms.string("isTauGenJetMatched")
    common_ntuple_branches.isTauGenJetMatched = cms.string("isTauGenJetMatched")
    # Careful! Only use GenTauJet (returns the values of the generated tau Jet) values if "bool TauTrigMatch::GenHadTauMatch()" returns "true". Otherwise it contains (unrealsitic) default values
    common_ntuple_branches.GenTauJet_Pt = cms.string("GenTauJet.pt")
    common_ntuple_branches.GenTauJet_Eta = cms.string("GenTauJet.eta")
    common_ntuple_branches.GenTauJet_Phi = cms.string("GenTauJet.phi")
    common_ntuple_branches.GenTauJet_DecayMode = cms.string("genDecayMode") # Decay Modes encoded as: 5 * (Prong - 1) + numberOfPi0 (example: oneProng2Pi0 = 2, threeProng1Pi0 = 11) cf. recoDecayModes. Additional decays: oneProngOther = -1, threeProngOther = -10, rare = -100 and NoGenJet = -999 


process.tauNTuple = cms.EDAnalyzer('TauValidationNTupleProd',
      tauTag         = cms.InputTag("patTaus"),
      altTauTag      = cms.InputTag("patTaus"),
      trigTag        = cms.InputTag("patTriggerEvent"),
      runOnMC        = cms.bool(runOnMC),
      ntuple         = common_ntuple_branches,
      maxDR          = cms.double(0.5), #The DeltaR parameter used for the trigger matching
      filterNames    = cms.vstring(),
)

###############################################

for j in range(len(filterName)):
    setattr(common_ntuple_branches, filterName[j], cms.string( "isTrigObjMatched(%i)"%j) )

for j in range(len(IDName)):
    setattr(common_ntuple_branches, IDName[j], cms.string( "recoTauCandID(\"%s\")"%IDName[j]) )

process.tauNTuple.filterNames = cms.vstring( filterName )


process.p = cms.Path(
        process.tauNTuple
        )
