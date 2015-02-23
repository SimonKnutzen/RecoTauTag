
//
// Original Author:  Simon Knutzen
//         Created:  Mon Apr 15 18:03:26 CEST 2013
// $Id: TauValidationNTupleProd.cc,v 1.5 2013/05/14 13:25:46 knutzen Exp $
//
//


#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include <math.h>
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h" 
#include "RecoTauTag/TauAnalysisTools/interface/ExpressionNtuple.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"

#include "RecoTauTag/TauAnalysisTools/interface/TauInfoContainer.h"


class TauValidationNTupleProd : public edm::EDAnalyzer {
   public:
      explicit TauValidationNTupleProd(const edm::ParameterSet&);
      ~TauValidationNTupleProd();


   private:

      uct::ExpressionNtuple<TauInfoContainer> ntuple_;

      edm::InputTag tauSrc_;
      edm::InputTag altTauSrc_;
      edm::InputTag triggerSrc_;
      double maxDR_;
      unsigned int Nvtx_;
      std::vector< std::string > filtNames ;
      TH1F * h_NumEvents;
      bool isMC_;
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      const reco::Vertex* getVertexCollection(const edm::Event& evt);
      const GenEventInfoProduct* getGenEvtInfo(const edm::Event& evt);
      std::vector<const reco::GenParticle*> getGenParticleCollection(const edm::Event& evt);
      std::vector<const pat::Tau*> getRecoCandCollections(const edm::Event& evt, const edm::InputTag& collection);
      std::vector<const reco::Candidate*> getTrigObjCandCollections(const edm::Event& evt, const edm::InputTag& collection, const std::string& filtername);
      std::vector<const reco::GsfElectron*> getGsfElectronCollections(const edm::Event& evt);
      const reco::GsfElectron* findBestgsfEleMatch(const pat::Tau* TagTauObj, std::vector<const reco::GsfElectron*>& gsfEle, double maxDR);
      const reco::Candidate* findBestMatch(const pat::Tau* TagTauObj,std::vector<const reco::Candidate*>& FilterSelection, double maxDR);
      const reco::GenParticle* findBestGenMatch(const pat::Tau* TagTauObj,std::vector<const reco::GenParticle*>& GenPart, double maxDR);
      std::vector<const reco::Candidate*> getRecoJetCollections(const edm::Event& evt);
      std::vector<const reco::Candidate*> getGenJetCollections(const edm::Event& evt);
      const reco::Candidate* findBestRecoObjMatch(const pat::Tau* TagTauObj,std::vector<const reco::Candidate* >& RecoObjSel, double maxDR);
      const pat::Tau* findBestTauObjMatch(const pat::Tau* TagTauObj,std::vector<const pat::Tau* >& TauObjSel, double maxDR);
      //std::vector<const GEMRecHitCollection*> getGemRecHitsCollections(const edm::Event& evt);
};



TauValidationNTupleProd::TauValidationNTupleProd(const edm::ParameterSet& iConfig):
    ntuple_(iConfig.getParameterSet("ntuple"))
{
      edm::Service<TFileService> fs;

      Nvtx_=0;

      ntuple_.initialize(*fs);
      h_NumEvents    = fs->make<TH1F>( "counter"  , "counter", 1,  0., 1. );

      tauSrc_           =   iConfig.getParameter<edm::InputTag>("tauTag");
      altTauSrc_        =   iConfig.getParameter<edm::InputTag>("altTauTag");
      triggerSrc_       =   iConfig.getParameter<edm::InputTag>("trigTag");
      maxDR_            =   iConfig.getParameter<double>("maxDR");
      filtNames         =   iConfig.getParameter< std::vector<std::string> >("filterNames");
      isMC_             =   iConfig.getParameter<bool>("runOnMC");

}

TauValidationNTupleProd::~TauValidationNTupleProd()
{ 

}

// Get vertex collection 
 
edm::Handle<const reco::Vertex* > vertices;

const reco::Vertex* TauValidationNTupleProd::getVertexCollection(const edm::Event& evt) {
    const reco::Vertex* output;
    edm::Handle< std::vector<reco::Vertex > > handle;
    evt.getByLabel("offlinePrimaryVertices", handle);    
    Nvtx_ = handle->size();
    //std::cout << handle->size() << std::endl;
    // Loop over objects in current collection

    const reco::Vertex& object = handle->at(0); // The primary vertex collection is sorted according to the sum of the Pt squared of the tracks associated to each vertex, such that the vertex with largest sum, likely to be the "signal" vertex, appears first. Some justification of this can be found in the Higgs note CMS-AN-11-129.
    output = &object;

  return output;
}

const GenEventInfoProduct* TauValidationNTupleProd::getGenEvtInfo(const edm::Event& evt) {
    //GenEventInfoProduct* const output;
    //edm::Handle<GenEventInfoProduct> evt_info;
    //evt.getByType(evt_info);
    //
    edm::Handle<GenEventInfoProduct> evt_info;
    evt.getByLabel("generator",evt_info);
    GenEventInfoProduct const *output = &(*evt_info);
    //output = 0; 
    return output;
}

// Get collection of generator particles with status 2

std::vector<const reco::GenParticle*> TauValidationNTupleProd::getGenParticleCollection(const edm::Event& evt) {
    std::vector<const reco::GenParticle*> output;
    edm::Handle< std::vector<reco::GenParticle> > handle;
    evt.getByLabel("genParticles", handle);
    // Loop over objects in current collection
    for (size_t j = 0; j < handle->size(); ++j) {
      const reco::GenParticle& object = handle->at(j);
      //if(fabs(object.pdgId())==15 && object.status() == 2) output.push_back(&object);
      if(object.pt()>15. && fabs(object.eta())< 2.5){
         if(object.status() == 2) output.push_back(&object);
      }
    }
  return output;
}


// Get collection of pat::taus
//
std::vector<const pat::Tau*> TauValidationNTupleProd::getRecoCandCollections(const edm::Event& evt, const edm::InputTag& collection) {
    std::vector<const pat::Tau*> output;
    edm::Handle< std::vector<pat::Tau> > handle;
    evt.getByLabel(collection, handle);
    // Loop over objects in current collection


    for (size_t j = 0; j < handle->size(); ++j) {
      const pat::Tau& object = handle->at(j);
 //     if(object.genJet() != 0) std::cout << "Gen decay Mode: " << (JetMCTagUtils::genTauDecayMode(*(object.genJet()))).c_str() << std::endl;
     if(object.pt()>15. && fabs(object.eta())< 2.5){
        output.push_back(&object);
      }
    }
  return output;
}
// Get gemRecHits colelction:
//
//std::vector<const GEMRecHitCollection*> TauValidationNTupleProd::getGemRecHitsCollections(const edm::Event& evt) {
//    std::vector<const GEMRecHitCollection*> output;
//    edm::Handle<GEMRecHitCollection> gemRecHits;
//    evt.getByLabel("gemRecHits","",gemRecHits);
//    for (size_t j = 0; j < gemRecHits->size(); ++j) {
//        output.push_back(gemRecHits->at(j));
//    }
//  return output;
//}
//
//std::vector<const reco::GenJet*> TauValidationNTupleProd::getGenJetCollections(const edm::Event& evt) {
//    std::vector<const reco::GenJet*> output;
//    edm::Handle< std::vector<reco::GenJet> > handle;
//    evt.getByLabel("ak5GenJets", handle);
//    for (size_t j = 0; j < handle->size(); ++j) {
//      const ak5GenJets& object =   (handle->at(j));
//        output.push_back(&object);
//    }
//  return output;
//}

std::vector<const reco::Candidate*> TauValidationNTupleProd::getGenJetCollections(const edm::Event& evt) {
    std::vector<const reco::Candidate*> output;
    edm::Handle< std::vector<reco::GenJet> > handle;
    evt.getByLabel("ak5GenJets", handle);
    for (size_t j = 0; j < handle->size(); ++j) {
      const reco::Candidate& object =  dynamic_cast< const reco::Candidate& > (handle->at(j));
      if(object.pt()>15. && fabs(object.eta())< 2.5){
        output.push_back(&object);
      }
    }
  return output;
}
std::vector<const reco::Candidate*> TauValidationNTupleProd::getRecoJetCollections(const edm::Event& evt) {
    std::vector<const reco::Candidate*> output;
    edm::Handle< std::vector<pat::Jet> > handle;
    evt.getByLabel("patJetsPFlow", handle);
    for (size_t j = 0; j < handle->size(); ++j) {
      const reco::Candidate& object =  dynamic_cast< const reco::Candidate& > (handle->at(j));
      if(object.pt()>15. && fabs(object.eta())< 2.5){
        output.push_back(&object);
      }
    }
  return output;
}
//Get gsfElectronCollection
//
std::vector<const reco::GsfElectron*> TauValidationNTupleProd::getGsfElectronCollections(const edm::Event& evt) {
    std::vector<const reco::GsfElectron*> output;
    edm::Handle<reco::GsfElectronCollection> gsfElectrons;
    evt.getByLabel("gsfElectrons",gsfElectrons);

    for (size_t j = 0; j < gsfElectrons->size(); ++j) {
      const reco::GsfElectron& object = (gsfElectrons->at(j));
        output.push_back(&object);
    }
  return output;
}

//Match gsfEle to tagTau
const reco::GsfElectron* TauValidationNTupleProd::findBestgsfEleMatch(const pat::Tau* TagTauObj,
    std::vector<const reco::GsfElectron*>& gsfEle, double maxDR) {
  const reco::GsfElectron*  output = NULL;
  double maxPt = 0;
  for (size_t i = 0; i < gsfEle.size(); ++i) {
    double deltaR = reco::deltaR(*TagTauObj, *gsfEle[i]);
    double PT = gsfEle[i]->pt();
    if (deltaR < maxDR && gsfEle[i]->pt() > 10. ) {
      if (!output || ( PT > maxPt )) {
        output = gsfEle[i];
        maxPt = PT;
      }
    }
  }
  return output;
}

//Match pfJet to tagTau
const reco::Candidate* TauValidationNTupleProd::findBestRecoObjMatch(const pat::Tau* TagTauObj,
    std::vector<const reco::Candidate* >& RecoObjSel, double maxDR) {
  const reco::Candidate*  output = NULL;
  double bestDeltaR = -1;
  for (size_t i = 0; i < RecoObjSel.size(); ++i) {
    double deltaR = reco::deltaR(*TagTauObj, *RecoObjSel[i]);
    if (deltaR < maxDR) {
      if (!output || deltaR < bestDeltaR) {
        output = RecoObjSel[i];
        bestDeltaR = deltaR;
      }
    }
  }
  return output;
}

const pat::Tau* TauValidationNTupleProd::findBestTauObjMatch(const pat::Tau* TagTauObj,
    std::vector<const pat::Tau* >& TauObjSel, double maxDR) {
  const pat::Tau*  output = NULL;
  double bestDeltaR = -1;
  for (size_t i = 0; i < TauObjSel.size(); ++i) {
    double deltaR = reco::deltaR(*TagTauObj, *TauObjSel[i]);
    if (deltaR < maxDR) {
      if (!output || deltaR < bestDeltaR) {
        output = TauObjSel[i];
        bestDeltaR = deltaR;
      }
    }
  }
  return output;
}

// Get vector of collections of TriggerFilterObjects
//
std::vector<const reco::Candidate*> TauValidationNTupleProd::getTrigObjCandCollections(const edm::Event& evt, const edm::InputTag& collection, const std::string& filtername) {
  std::vector<const reco::Candidate*> output;
    edm::Handle<pat::TriggerEvent> triggerEv;
    evt.getByLabel(collection, triggerEv);
    pat::TriggerObjectRefVector FilterObjects = triggerEv->filterObjects(filtername);

    for (unsigned int i = 0; i < FilterObjects.size(); i++) {
      const reco::Candidate &object = dynamic_cast< const reco::Candidate& >( *(FilterObjects.at(i)) );
      output.push_back(&object);
    }
  return output;
}

// Method to find the best match between tag tau and trigger filter object. The best matched filter object will be returned. If there is no match within a DR < 0.5, a null pointer is returned
const reco::Candidate*  TauValidationNTupleProd::findBestMatch(const pat::Tau* TagTauObj,
    std::vector<const reco::Candidate*>& FilterSelection, double maxDR) {
  const reco::Candidate* output = NULL;
  double bestDeltaR = -1;
  for (size_t i = 0; i < FilterSelection.size(); ++i) {
    double deltaR = reco::deltaR(*TagTauObj, *FilterSelection[i]);
    if (deltaR < maxDR) {
      if (!output || deltaR < bestDeltaR) {
        output = FilterSelection[i];
        bestDeltaR = deltaR;
      }
    }
  }
  return output;
}

// Method to find the best match between tag tau and gen object. The best matched gen tau object will be returned. If there is no match within a DR < 0.5, a null pointer is returned
const reco::GenParticle* TauValidationNTupleProd::findBestGenMatch(const pat::Tau* TagTauObj,
    std::vector<const reco::GenParticle*>& GenPart, double maxDR) {
  const reco::GenParticle* output = NULL;
  double bestDeltaR = -1;
  for (size_t i = 0; i < GenPart.size(); ++i) {
    double deltaR = reco::deltaR(*TagTauObj, *GenPart[i]);
    if (deltaR < maxDR) {
      if (!output || deltaR < bestDeltaR) {
        output = GenPart[i];
        bestDeltaR = deltaR;
      }
    }
  }
  return output;
}

void
TauValidationNTupleProd::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    unsigned int evtInfo[3] = {iEvent.run(),iEvent.luminosityBlock(),(iEvent.eventAuxiliary()).event()};

    std::vector<const pat::Tau*> tauObjects = getRecoCandCollections(iEvent, tauSrc_);
    std::vector<const pat::Tau*> altTauObjects = getRecoCandCollections(iEvent, altTauSrc_);
    std::vector<const reco::Candidate*> pfJets = getRecoJetCollections(iEvent);
    std::vector<const reco::Candidate*> genJets = getGenJetCollections(iEvent);
    std::vector<const reco::GsfElectron*> gsfElectrons = getGsfElectronCollections(iEvent);
    std::vector<const reco::GenParticle*> GenObjects;
    if(isMC_) GenObjects = getGenParticleCollection(iEvent);
    const reco::Vertex* Vertex = getVertexCollection(iEvent);
    std::vector<std::vector<const reco::Candidate*>> allTrigObjects;


    const GenEventInfoProduct* genInfo;
    if(isMC_){
      genInfo = getGenEvtInfo(iEvent);
    }
    else{
      genInfo = NULL;
    }
    //for(unsigned int i = 0; i < filtNames.size(); i++){

    //     std::vector<const reco::Candidate*> trigObjects = getTrigObjCandCollections(iEvent, triggerSrc_, filtNames[i]);
    //     allTrigObjects.push_back(trigObjects); // enter collection of trigger objects for each filter in a vector
    //}

    std::vector<TauInfoContainer*> matches;
    TauInfoContainer *theMatch = NULL;

    for(unsigned int i = 0; i<tauObjects.size(); i++){

            const pat::Tau* TagTau = tauObjects[i];
            std::vector<const reco::Candidate* > allBestFilterMatches;
            const reco::Candidate* bestGenMatch = findBestGenMatch(TagTau,GenObjects, maxDR_) ;
            const reco::Candidate* bestJetMatch = findBestRecoObjMatch(TagTau, pfJets, maxDR_);
            const reco::Candidate* bestGenJetMatch = findBestRecoObjMatch(TagTau, genJets, maxDR_);
            const pat::Tau* bestAltTauObjMatch = findBestTauObjMatch(TagTau, altTauObjects, maxDR_);
            const reco::GsfElectron* bestEle = findBestgsfEleMatch(TagTau, gsfElectrons, 0.3);
            
            EleMVAContainer* EleMVA = NULL;
            if( bestEle != 0 &&  TagTau->tauID("decayModeFindingNewDMs") > 0.5 ) EleMVA = new EleMVAContainer( bestEle, TagTau );
            else EleMVA = new EleMVAContainer();

            //if( bestEle != 0 ) std::cout << bestEle->pt() << std::endl;
            //const reco::Candidate*  jetTest = NULL; // For test

            allBestFilterMatches.push_back( NULL ); // For test

            //for(unsigned int j = 0; j < allTrigObjects.size(); j++){

            //    const reco::Candidate* bestFilterMatch = findBestMatch(TagTau, allTrigObjects[j], maxDR_);
            //    allBestFilterMatches.push_back(bestFilterMatch); // enter the best matched trigger object for each filter into a vector 
            //       
            //}    

            //theMatch = new TauInfoContainer(TagTau,bestAltTauObjMatch,&allBestFilterMatches,bestGenMatch,matches.size(),tauObjects.size(), Nvtx_, jetTest, Vertex); // create a TauInfoContainer object for each tag tau
            //theMatch = new TauInfoContainer(TagTau,bestAltTauObjMatch,&allBestFilterMatches,bestGenMatch,matches.size(),tauObjects.size(), Nvtx_, &iEvent, bestJetMatch, Vertex); // create a TauInfoContainer object for each tag tau
            theMatch = new TauInfoContainer(TagTau,bestAltTauObjMatch,&allBestFilterMatches,bestGenMatch,matches.size(),tauObjects.size(), genInfo, Nvtx_ ,evtInfo , bestJetMatch, Vertex, EleMVA, bestGenJetMatch); // create a TauInfoContainer object for each tag tau

            theMatch->genDecayMode();

            matches.push_back(theMatch); 

    }


    for (size_t i = 0; i < matches.size(); ++i) {
       ntuple_.fill(*matches.at(i));  // create TTree
    }

    h_NumEvents->Fill(0);
}


void 
TauValidationNTupleProd::beginJob()
{
}

void 
TauValidationNTupleProd::endJob() 
{
}

DEFINE_FWK_MODULE(TauValidationNTupleProd);
