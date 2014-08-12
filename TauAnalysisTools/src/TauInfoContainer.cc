
#include "RecoTauTag/TauAnalysisTools/interface/TauInfoContainer.h"
#include <TLorentzVector.h>
#include "Math/GenVector/LorentzVector.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

using namespace edm;



EleMVAContainer::EleMVAContainer(){
   hasValidEleMVA_ = 0;
   Tau_GammaEnFrac_ = -1;    
   Tau_GSFTracklnPt_ = -1;
   Tau_GSFTrackEta_ = -1;
   Tau_VisMass_ = -1;
   Tau_EmFraction_ = -1;
   Tau_GSFTrackResol_ = -1;
   Tau_Pt_ = -1;
   Tau_LeadChargedPFCandEtaAtEcalEntrance_ = -1;
   Tau_LeadChargedPFCandPt_ = -1;
   Tau_NumGammaCands_ = -1;
   Tau_HadrHoP_ = -1;
   Tau_GammaEtaMom_ = -1;
   Tau_GammaPhiMom_ = -1;
   Tau_GammaPhiMom_ = -1;
   Tau_HadrEoP_ = -1;
   Tau_EtaAtEcalEntrance_ = -1;
   Tau_GSFChi2_ = -1;
   Tau_GSFNumHits_ = -1;
   Tau_KFNumHits_ = -1;
   Tau_HadrMva_ = -1;

    Elec_GSFTrackEta_ = -1;
    Elec_EtotOverPin_ = -1;
    Elec_Fbrem_ = -1;
    Elec_GSFNumHits_ = -1;
    Elec_Chi2GSF_ = -1;
    Elec_GSFTracklnPt_ = -1;
    Elec_EgammaOverPdif_ = -1;
    Elec_GSFTrackResol_ = -1;
}

EleMVAContainer::~EleMVAContainer(){}

EleMVAContainer::EleMVAContainer( const reco::GsfElectron* bestEle, const pat::Tau* recoTauCand ):
bestEle_(bestEle), recoTauCand_(recoTauCand)
{
   hasValidEleMVA_ = 1;
   calcMVAvar();
}

const reco::GsfElectron* EleMVAContainer::bestElectron() const {
   return bestEle_;
}


void EleMVAContainer::calcMVAvar(){

   float sumEtaTimesEnergy = 0.;
   float sumPhiTimesEnergy = 0.;
   float sumEnergy = 0.;
   float sumEtaTimesEnergyEcalEnWeighted = 0.;
   float sumPhiTimesEnergyEcalEnWeighted = 0.;
   float sumEnergyEcalEnWeighted = 0.;
   const std::vector<reco::PFCandidatePtr>& signalPFCands = recoTauCand_->signalPFCands();
   for ( std::vector<reco::PFCandidatePtr>::const_iterator pfCandidate = signalPFCands.begin();
   pfCandidate != signalPFCands.end(); ++pfCandidate ) {
     sumEtaTimesEnergy += ((*pfCandidate)->positionAtECALEntrance().eta()*(*pfCandidate)->energy());
     sumPhiTimesEnergy += ((*pfCandidate)->positionAtECALEntrance().phi()*(*pfCandidate)->energy());
     sumEnergy += (*pfCandidate)->energy();
     sumEtaTimesEnergyEcalEnWeighted += ((*pfCandidate)->positionAtECALEntrance().eta()*(*pfCandidate)->ecalEnergy());
     sumPhiTimesEnergyEcalEnWeighted += ((*pfCandidate)->positionAtECALEntrance().phi()*(*pfCandidate)->ecalEnergy());
     sumEnergyEcalEnWeighted += (*pfCandidate)->ecalEnergy();
   }
   if ( sumEnergy > 0. ) {
     Tau_EtaAtEcalEntrance_ = sumEtaTimesEnergy/sumEnergy;
     Tau_PhiAtEcalEntrance_ = sumPhiTimesEnergy/sumEnergy;
   }
   if ( sumEnergyEcalEnWeighted > 0. ) {
     Tau_EtaAtEcalEntranceEcalEnWeighted_ = sumEtaTimesEnergyEcalEnWeighted/sumEnergyEcalEnWeighted;
     Tau_PhiAtEcalEntranceEcalEnWeighted_ = sumPhiTimesEnergyEcalEnWeighted/sumEnergyEcalEnWeighted;
   }
   for ( std::vector<reco::PFCandidatePtr>::const_iterator pfCandidate = signalPFCands.begin();
   pfCandidate != signalPFCands.end(); ++pfCandidate ) {
     const reco::Track* track = 0;
     if ( (*pfCandidate)->trackRef().isNonnull() ) track = (*pfCandidate)->trackRef().get();
     else if ( (*pfCandidate)->muonRef().isNonnull() && (*pfCandidate)->muonRef()->innerTrack().isNonnull() ) track = (*pfCandidate)->muonRef()->innerTrack().get();
     else if ( (*pfCandidate)->muonRef().isNonnull() && (*pfCandidate)->muonRef()->globalTrack().isNonnull() ) track = (*pfCandidate)->muonRef()->globalTrack().get();
     else if ( (*pfCandidate)->muonRef().isNonnull() && (*pfCandidate)->muonRef()->outerTrack().isNonnull() ) track = (*pfCandidate)->muonRef()->outerTrack().get();
     else if ( (*pfCandidate)->gsfTrackRef().isNonnull() ) track = (*pfCandidate)->gsfTrackRef().get();
     if ( track ) {
        if ( track->pt() > Tau_LeadChargedPFCandPt_ ) {
           Tau_LeadChargedPFCandEtaAtEcalEntrance_ = (*pfCandidate)->positionAtECALEntrance().eta();
           Tau_LeadChargedPFCandPhiAtEcalEntrance_ = (*pfCandidate)->positionAtECALEntrance().phi();
           Tau_LeadChargedPFCandPt_ = track->pt();
        }
     } else {
         if ( (*pfCandidate)->pt() > Tau_LeadNeutralPFCandPt_ ) {
           Tau_LeadNeutralPFCandEtaAtEcalEntrance_ = (*pfCandidate)->positionAtECALEntrance().eta();
           Tau_LeadNeutralPFCandPhiAtEcalEntrance_ = (*pfCandidate)->positionAtECALEntrance().phi();
           Tau_LeadNeutralPFCandPt_ = (*pfCandidate)->pt();
         }
     }
   }
   Tau_Eta_ = recoTauCand_->eta();
   Tau_Pt_ = recoTauCand_->pt();
   Tau_Phi_ = recoTauCand_->phi();
   Tau_EmFraction_ = TMath::Max(recoTauCand_->emFraction(), float(0.));
   Tau_NumChargedCands_ = recoTauCand_->signalPFChargedHadrCands().size();
   Tau_NumGammaCands_ = recoTauCand_->signalPFGammaCands().size();
   
   if ( recoTauCand_->leadPFChargedHadrCand().isNonnull() ) {
     Tau_LeadHadronPt_ = recoTauCand_->leadPFChargedHadrCand()->pt();
     Tau_HasGsf_ = (recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull();
   }
   if ( (recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull() ) {
     Tau_GSFChi2_ = recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()->normalizedChi2();
     Tau_GSFNumHits_ = recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()->numberOfValidHits();
     Tau_GSFNumPixelHits_ = recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()->hitPattern().numberOfValidPixelHits();
     Tau_GSFNumStripHits_ = recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()->hitPattern().numberOfValidStripHits();
     Tau_GSFTrackResol_ = recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()->ptError()/recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()->pt();
     Tau_GSFTracklnPt_ = log(recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()->pt())*TMath::Ln10();
     Tau_GSFTrackEta_ = recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()->eta();
   }
   
   if ( recoTauCand_->leadPFChargedHadrCand().isNonnull() ) {
     Tau_HasKF_ = (recoTauCand_->leadPFChargedHadrCand()->trackRef()).isNonnull();
   }
   if ( (recoTauCand_->leadPFChargedHadrCand()->trackRef()).isNonnull() ) {
     Tau_KFChi2_ = recoTauCand_->leadPFChargedHadrCand()->trackRef()->normalizedChi2();
     Tau_KFNumHits_ = recoTauCand_->leadPFChargedHadrCand()->trackRef()->numberOfValidHits();
     Tau_KFNumPixelHits_ = recoTauCand_->leadPFChargedHadrCand()->trackRef()->hitPattern().numberOfValidPixelHits();
     Tau_KFNumStripHits_ = recoTauCand_->leadPFChargedHadrCand()->trackRef()->hitPattern().numberOfValidStripHits();
     Tau_KFTrackResol_ = recoTauCand_->leadPFChargedHadrCand()->trackRef()->ptError()/recoTauCand_->leadPFChargedHadrCand()->trackRef()->pt();
     Tau_KFTracklnPt_ = log(recoTauCand_->leadPFChargedHadrCand()->trackRef()->pt())*TMath::Ln10();
     Tau_KFTrackEta_ = recoTauCand_->leadPFChargedHadrCand()->trackRef()->eta();
   }
   
   if ( recoTauCand_->leadPFChargedHadrCand().isNonnull() ) {
     Tau_HadrHoP_ = recoTauCand_->leadPFChargedHadrCand()->hcalEnergy()/recoTauCand_->leadPFChargedHadrCand()->p();
     Tau_HadrEoP_ = recoTauCand_->leadPFChargedHadrCand()->ecalEnergy()/recoTauCand_->leadPFChargedHadrCand()->p();
   }
   
   GammasdEta_.clear();
   GammasdPhi_.clear();
   GammasPt_.clear();
   const std::vector<reco::PFCandidatePtr>& signalPFGammaCands = recoTauCand_->signalPFGammaCands();
   for ( std::vector<reco::PFCandidatePtr>::const_iterator pfGamma = signalPFGammaCands.begin();
     pfGamma != signalPFGammaCands.end(); ++pfGamma ) {
     if ( recoTauCand_->leadPFChargedHadrCand().isNonnull() ) {
        GammasdEta_.push_back((*pfGamma)->eta() - recoTauCand_->leadPFChargedHadrCand()->eta());
        GammasdPhi_.push_back((*pfGamma)->phi() - recoTauCand_->leadPFChargedHadrCand()->phi());
     } else {
        GammasdEta_.push_back((*pfGamma)->eta() - recoTauCand_->eta());
        GammasdPhi_.push_back((*pfGamma)->phi() - recoTauCand_->phi());
     }
     GammasPt_.push_back((*pfGamma)->pt());
   }
   
   float sumPt = 0.;
   float dEta = 0.;
   float dEta2 = 0.;
   float dPhi = 0.;
   float dPhi2 = 0.;
   float sumPt2 = 0.;
   size_t numPFGammas = GammasPt_.size();
   assert(GammasdEta_.size() == numPFGammas);
   assert(GammasdPhi_.size() == numPFGammas);
   for ( size_t idxPFGamma = 0; idxPFGamma < numPFGammas; ++idxPFGamma ) {
     float gamma_pt = GammasPt_[idxPFGamma];
     float gamma_dPhi = GammasdPhi_[idxPFGamma];
     if ( gamma_dPhi > TMath::Pi() ) gamma_dPhi -= 2.*TMath::Pi();
     else if ( gamma_dPhi < -TMath::Pi() ) gamma_dPhi += 2.*TMath::Pi();
     float gamma_dEta = GammasdEta_[idxPFGamma];
     sumPt += gamma_pt;
     sumPt2 += (gamma_pt*gamma_pt);
     dEta += (gamma_pt*gamma_dEta);
     dEta2 += (gamma_pt*gamma_dEta*gamma_dEta);
     dPhi += (gamma_pt*gamma_dPhi);
     dPhi2 += (gamma_pt*gamma_dPhi*gamma_dPhi);
   }
   
   float gammadPt = sumPt/recoTauCand_->pt();
   
   if ( sumPt > 0. ) {
     dEta /= sumPt;
     dPhi /= sumPt;
     dEta2 /= sumPt;
     dPhi2 /= sumPt;
   }
   
   Tau_GammaEtaMom_ = TMath::Sqrt(dEta2)*TMath::Sqrt(gammadPt)*recoTauCand_->pt();
   Tau_GammaPhiMom_ = TMath::Sqrt(dPhi2)*TMath::Sqrt(gammadPt)*recoTauCand_->pt();
   Tau_GammaEnFrac_ = gammadPt;
   Tau_VisMass_ = recoTauCand_->mass();
   Tau_HadrMva_ = TMath::Max(recoTauCand_->electronPreIDOutput(), float(-1.));

   Elec_HasSC_ = 0;
   reco::SuperClusterRef pfSuperCluster = bestEle_->pflowSuperCluster();
   if ( pfSuperCluster.isNonnull() && pfSuperCluster.isAvailable() ) {
      Elec_HasSC_ = 1;
      Elec_Ee_ = 0.;
      Elec_Egamma_ = 0.;
      for ( reco::CaloCluster_iterator pfCluster = pfSuperCluster->clustersBegin();
      pfCluster != pfSuperCluster->clustersEnd(); ++pfCluster ) {
      float pfClusterEn = (*pfCluster)->energy();
      if ( pfCluster == pfSuperCluster->clustersBegin() ) Elec_Ee_ += pfClusterEn;
      else Elec_Egamma_ += pfClusterEn;
      }
      Elec_Pin_ = TMath::Sqrt(bestEle_->trackMomentumAtVtx().Mag2());
      Elec_Pout_ = TMath::Sqrt(bestEle_->trackMomentumOut().Mag2());
      Elec_EtotOverPin_ = (Elec_Ee_ + Elec_Egamma_)/Elec_Pin_;
      Elec_EeOverPout_ = Elec_Ee_/Elec_Pout_;
      Elec_EgammaOverPdif_ = Elec_Egamma_/(Elec_Pin_ - Elec_Pout_);
      Elec_HoHplusE_ = bestEle_->mvaInput().hadEnergy/(bestEle_->mvaInput().hadEnergy+Elec_Ee_) ;
   }
   Elec_HasKF_ = 0;
   if ( bestEle_->closestCtfTrackRef().isNonnull() ) {
      Elec_HasKF_ = 1;
      Elec_Chi2KF_ = bestEle_->closestCtfTrackRef()->normalizedChi2();
      Elec_KFNumHits_ = bestEle_->closestCtfTrackRef()->numberOfValidHits();
      Elec_KFNumPixelHits_ = bestEle_->closestCtfTrackRef()->hitPattern().numberOfValidPixelHits();
      Elec_KFNumStripHits_ = bestEle_->closestCtfTrackRef()->hitPattern().numberOfValidStripHits();
           Elec_KFTrackResol_ = bestEle_->closestCtfTrackRef()->ptError()/bestEle_->closestCtfTrackRef()->pt();
      Elec_KFTracklnPt_ = log(bestEle_->closestCtfTrackRef()->pt())*TMath::Ln10();
      Elec_KFTrackEta_ = bestEle_->closestCtfTrackRef()->eta();
   }
      
   // Variables related to the GsfTrack
   Elec_HasGSF_ = 0;
   if ( bestEle_->gsfTrack().isNonnull() ) {
      Elec_HasGSF_ = 1;
      Elec_Chi2GSF_ = bestEle_->gsfTrack()->normalizedChi2();
      Elec_GSFNumHits_ = bestEle_->gsfTrack()->numberOfValidHits();
      Elec_GSFNumPixelHits_ = bestEle_->gsfTrack()->hitPattern().numberOfValidPixelHits();
      Elec_GSFNumStripHits_ = bestEle_->gsfTrack()->hitPattern().numberOfValidStripHits();
      Elec_GSFTrackResol_ = bestEle_->gsfTrack()->ptError()/bestEle_->gsfTrack()->pt();
      Elec_GSFTracklnPt_ = log(bestEle_->gsfTrack()->pt())*TMath::Ln10();
      Elec_GSFTrackEta_ = bestEle_->gsfTrack()->eta();
   }


}
////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////



//TauInfoContainer::TauInfoContainer(const pat::Tau* recoTauCand, const pat::Tau* altTauObj, std::vector<const reco::Candidate*>* trigObj, const reco::Candidate* GenParticle ,unsigned int index, unsigned int nTotalObjects, const GenEventInfoProduct* GenInfo, unsigned int NVTX, const edm::Event* evt, const reco::Candidate* pfJet, const reco::Vertex* Vertex ):
//  recoTauCand_(recoTauCand),altTauObj_(altTauObj), trigObj_(trigObj),GenParticle_(GenParticle),index_(index), nTotalObjects_(nTotalObjects), genInfo_(GenInfo),Nvtx_(NVTX),Evt_(evt), pfJet_(pfJet), Vertex_(Vertex){
//  
//        // Create a dummy reco::Candidate Object with unrealistic LorentzVector values as a default output to return in case of a failed matching.  
//        dummyCandidate_ = dynamic_cast<reco::Candidate* >( recoTauCand->clone());
//        math::XYZTLorentzVector *v = new math::XYZTLorentzVector();
//        v->SetPxPyPzE(-999.,-999.,-9999.,-999.);
//        dummyCandidate_->setP4(((const math::XYZTLorentzVector)*v)); 
//  
//        // Create a dummy reco::Candidate Object with unrealistic LorentzVector values as a default output to return in case of a failed matching.  
//        dummyCandidateTau_ = dynamic_cast<pat::Tau* >( recoTauCand->clone());
//        dummyCandidateTau_->setP4(((const math::XYZTLorentzVector)*v)); 
//  }
TauInfoContainer::TauInfoContainer(const pat::Tau* recoTauCand, const pat::Tau* altTauObj, std::vector<const reco::Candidate*>* trigObj, const reco::Candidate* GenParticle ,unsigned int index, unsigned int nTotalObjects, const GenEventInfoProduct* GenInfo, unsigned int NVTX, unsigned int* evtInfo, const reco::Candidate* pfJet, const reco::Vertex* Vertex, EleMVAContainer* eleMVA, const reco::Candidate* genJet ):
  recoTauCand_(recoTauCand),altTauObj_(altTauObj), trigObj_(trigObj),GenParticle_(GenParticle),index_(index), nTotalObjects_(nTotalObjects), genInfo_(GenInfo), Nvtx_(NVTX), evtInfo_(evtInfo), pfJet_(pfJet), Vertex_(Vertex), eleMVA_(eleMVA), genJet_(genJet){
         // Create a dummy reco::Candidate Object with unrealistic LorentzVector values as a default output to return in case of a failed matching.  
        dummyCandidate_ = dynamic_cast<reco::Candidate* >( recoTauCand->clone());
        math::XYZTLorentzVector *v = new math::XYZTLorentzVector();
        v->SetPxPyPzE(-999.,-999.,-9999.,-999.);
        dummyCandidate_->setP4(((const math::XYZTLorentzVector)*v)); 
         //Create a dummy reco::Candidate Object with unrealistic LorentzVector values as a default output to return in case of a failed matching.  
        dummyCandidateTau_ = dynamic_cast<pat::Tau* >( recoTauCand->clone());
        dummyCandidateTau_->setP4(((const math::XYZTLorentzVector)*v)); 
  }
TauInfoContainer::TauInfoContainer(){}
TauInfoContainer::~TauInfoContainer(){}

const reco::Vertex* TauInfoContainer::getPV() const{
   return Vertex_;
}

//float TauInfoContainer::Tau_GammaEnFrac() const{
//   return eleMVA_->Tau_GammaEnFrac_;
//}

const EleMVAContainer* TauInfoContainer::getEleMVA() const{
   return eleMVA_;
}

unsigned int TauInfoContainer::index() const {
   return index_;
}

unsigned int TauInfoContainer::nTotalObjects() const {
   return nTotalObjects_;
}

const pat::Tau* TauInfoContainer::recoTauCand() const {
   return recoTauCand_;
}

const pat::Tau* TauInfoContainer::altTauObj() const {
   if( altTauObj_!= NULL) return altTauObj_;
   else return dummyCandidateTau_; // Careful! Method return dummy object to ensure successfull termination of program. Only use GenParticle values if "bool TauInfoContainer::isGenParticelMatched()" returns "true"
}

const reco::Candidate* TauInfoContainer::PfJet() const {
   if(pfJet_ != NULL) return pfJet_;
   else return dummyCandidate_; // Careful! Method return dummy object to ensure successfull termination of program. Only use GenParticle values if "bool TauInfoContainer::isPfJetMatched()" returns "true"
}
const reco::Candidate* TauInfoContainer::GenJet() const {
   if(genJet_ != NULL) return genJet_;
   else return dummyCandidate_; // Careful! Method return dummy object to ensure successfull termination of program. Only use GenParticle values if "bool TauInfoContainer::isPfJetMatched()" returns "true"
}

bool TauInfoContainer::isPfJetMatched() const{
   return pfJet_ != NULL;
}
bool TauInfoContainer::isGenJetMatched() const{
   return genJet_ != NULL;
}
bool TauInfoContainer::hasValidTrack() const{
   if( recoTauCand_->leadPFChargedHadrCand().isNonnull() ){
      if(recoTauCand_->leadPFChargedHadrCand().get() != 0){
          if( recoTauCand_->leadPFChargedHadrCand().get()->trackRef().isNonnull() ){

            return true;
         }
        }
     }
   return false;
   
}

bool TauInfoContainer::Tau_HasGsf() const{
   bool output = 0;
   if( &(*(recoTauCand_->leadPFChargedHadrCand())) != 0 ){
      output = (recoTauCand_->leadPFChargedHadrCand()->gsfTrackRef()).isNonnull();
   }else{
      output = 0;
   }
   return output;
}
//const reco::Track* TauInfoContainer::Track() const {
//   if( recoTauCand_->leadPFChargedHadrCand().isNonnull() ) return &(*recoTauCand_->leadPFChargedHadrCand().get()->trackRef());
//   else return &(*dummyCandidatePF_->trackRef());
//}

const GenEventInfoProduct* TauInfoContainer::genInfo() const {
   return genInfo_;
}

double TauInfoContainer::RunNr() const {
   return evtInfo_[0];
}

double TauInfoContainer::EvtNr() const {
   return evtInfo_[2];
}

double TauInfoContainer::LumiSec() const {
   return evtInfo_[1];
}

double TauInfoContainer::pfCand_dz() const {
       std::vector<reco::PFCandidatePtr> pfCands = recoTauCand_->pfJetRef()->getPFConstituents();
         double temppt = -999;
         double out_dz = -999;
        for ( std::vector<reco::PFCandidatePtr>::const_iterator jetConstituent = pfCands.begin(); jetConstituent != pfCands.end(); ++jetConstituent ) {
         if((*jetConstituent)->trackRef().isNonnull()) {
 
           if(((*jetConstituent)->trackRef()->pt()) > temppt ) {  
           temppt = ((*jetConstituent)->trackRef()->pt());
           out_dz = ((*jetConstituent)->trackRef()->dz(Vertex_->position()) );  
      } 
   }
 
  }

   return out_dz;
}

double TauInfoContainer::TransImpPara() const {
     if( recoTauCand_->leadPFChargedHadrCand().isNonnull() ){
        if(recoTauCand_->leadPFChargedHadrCand().get() != 0){
           if( recoTauCand_->leadPFChargedHadrCand().get()->trackRef().isNonnull() ){
               return (*recoTauCand_->leadPFChargedHadrCand().get()->trackRef()).d0();
           }
        }
     }
     return -9;
}
double TauInfoContainer::TransImpParaError() const {
     if( recoTauCand_->leadPFChargedHadrCand().isNonnull() ){
        if(recoTauCand_->leadPFChargedHadrCand().get() != 0){
           if( recoTauCand_->leadPFChargedHadrCand().get()->trackRef().isNonnull() ){
               return (*recoTauCand_->leadPFChargedHadrCand().get()->trackRef()).d0Error();
           }
        }
     }
     return -9;
}

const reco::Candidate* TauInfoContainer::GenParticle() const {
   if(GenParticle_ != NULL) return GenParticle_;
   else return dummyCandidate_; // Careful! Method return dummy object to ensure successfull termination of program. Only use GenParticle values if "bool TauInfoContainer::isGenParticelMatched()" returns "true"
}

const reco::Candidate* TauInfoContainer::GenTauJet() const {
   if(recoTauCand_->genJet() != NULL) return recoTauCand_->genJet();
   else return dummyCandidate_; // Careful!  Method return dummy object to ensure successfull termination of program. Only use GenTauJet values if "bool TauInfoContainer::isTauGenJetMatched()" returns "true"

}

bool TauInfoContainer::isGenParticelMatched() const {
   return GenParticle_ != NULL;
}

bool TauInfoContainer::isAltTauObjMatched() const {
   return altTauObj_ != NULL;
}

bool TauInfoContainer::isTauGenJetMatched() const {
   return recoTauCand_->genJet() != NULL;
}

bool TauInfoContainer::isTrigObjMatched(int a) const {
   return trigObj_->at(a) != NULL;
}

double TauInfoContainer::recoTauCandID(std::string DiscriminatorName) const{

   return recoTauCand_->tauID(DiscriminatorName);  

} 

int TauInfoContainer::genDecayMode() const{

   std::string genDecayMode = "NoGenJet";

   if(recoTauCand_->genJet() != 0) genDecayMode = JetMCTagUtils::genTauDecayMode(*(recoTauCand_->genJet()));

   int prong = -999;
   int pi0 = -999;
   int posPi0 = -1;
   int result = -999;

   posPi0 = genDecayMode.find("Prong") + 5;


   if(genDecayMode.at(0)== 'o') prong = 0;
   else if(genDecayMode.at(0)=='t') prong = 10;

   if(posPi0!=4 && genDecayMode.at(posPi0)!='O' ) pi0 = genDecayMode.at(posPi0) - '0';

   result = prong + pi0;
   
   if(genDecayMode.at(0)== 'o' && genDecayMode.at(posPi0)=='O') result = -1; 
   if(genDecayMode.at(0)== 't' && genDecayMode.at(posPi0)=='O') result = -10;
   if(genDecayMode.at(0)== 'r') result = -100;
   if(genDecayMode.at(0)== 'N') result = -999;

//   if(result!=-999)std::cout << result <<" : " << genDecayMode << std::endl;

   return result;  

}    

int TauInfoContainer::Nvtx() const{
   return Nvtx_;
}
