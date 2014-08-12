#ifndef TauInfoContainer_h
#define TauInfoContainer_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

class EleMVAContainer {
   public: 
      EleMVAContainer();
      ~EleMVAContainer();
      EleMVAContainer( const reco::GsfElectron* bestEle, const pat::Tau* recoTauCand );

      const reco::GsfElectron* bestElectron() const;

      void calcMVAvar();

      bool hasValidEleMVA_;
      std::vector<float> GammasdEta_;
      std::vector<float> GammasdPhi_;
      std::vector<float> GammasPt_;
      int Tau_GsfEleMatch_;
      int Tau_GenEleMatch_;
      int Tau_GenHadMatch_;
      float Tau_Eta_;
      float Tau_EtaAtEcalEntrance_;
      float Tau_PhiAtEcalEntrance_;
      float Tau_EtaAtEcalEntranceEcalEnWeighted_;
      float Tau_PhiAtEcalEntranceEcalEnWeighted_;
      float Tau_LeadNeutralPFCandEtaAtEcalEntrance_;
      float Tau_LeadNeutralPFCandPhiAtEcalEntrance_;
      float Tau_LeadNeutralPFCandPt_;
      float Tau_LeadChargedPFCandEtaAtEcalEntrance_;
      float Tau_LeadChargedPFCandPhiAtEcalEntrance_;
      float Tau_LeadChargedPFCandPt_;
      float Tau_Pt_;
      float Tau_LeadHadronPt_;
      float Tau_Phi_;
      int Tau_HasGsf_;
      float Tau_GSFChi2_;
      int Tau_GSFNumHits_;
      int Tau_GSFNumPixelHits_;
      int Tau_GSFNumStripHits_;
      float Tau_GSFTrackResol_;
      float Tau_GSFTracklnPt_;
      float Tau_GSFTrackEta_;
      int Tau_HasKF_;
      float Tau_KFChi2_;
      int Tau_KFNumHits_;
      int Tau_KFNumPixelHits_;
      int Tau_KFNumStripHits_;
      float Tau_KFTrackResol_;
      float Tau_KFTracklnPt_;
      float Tau_KFTrackEta_;
      float Tau_EmFraction_;
      int Tau_NumChargedCands_;
      int Tau_NumGammaCands_;
      float Tau_HadrHoP_;
      float Tau_HadrEoP_;
      float Tau_VisMass_;
      float Tau_GammaEtaMom_;
      float Tau_GammaPhiMom_;
      float Tau_GammaEnFrac_;
      float Tau_HadrMva_;
      int Tau_DecayMode_;
      int Tau_MatchElePassVeto_;
      float Tau_VtxZ_;
      float Tau_zImpact_;
      
      int Elec_GenEleMatch_;
      int Elec_GenEleFromZMatch_;
      int Elec_GenEleFromZTauTauMatch_;
      int Elec_GenHadMatch_;
      int Elec_GenJetMatch_;
      float Elec_AbsEta_;
      float Elec_Pt_;
      int Elec_HasSC_;
      float Elec_PFMvaOutput_;
      float Elec_Ee_;
      float Elec_Egamma_;
      float Elec_Pin_;
      float Elec_Pout_;
      float Elec_EtotOverPin_;
      float Elec_EeOverPout_;
      float Elec_EgammaOverPdif_;
      int Elec_EarlyBrem_;
      int Elec_LateBrem_;
      float Elec_Logsihih_;
      float Elec_DeltaEta_;
      float Elec_HoHplusE_;
      float Elec_Fbrem_;
      int Elec_HasKF_;
      float Elec_Chi2KF_;
      int Elec_KFNumHits_;
      int Elec_KFNumPixelHits_;
      int Elec_KFNumStripHits_;
      float Elec_KFTrackResol_;
      float Elec_KFTracklnPt_;
      float Elec_KFTrackEta_;
      int Elec_HasGSF_;
      float Elec_Chi2GSF_;
      int Elec_GSFNumHits_;
      int Elec_GSFNumPixelHits_;
      int Elec_GSFNumStripHits_;
      float Elec_GSFTrackResol_;
      float Elec_GSFTracklnPt_;
      float Elec_GSFTrackEta_;
      
      int ElecVeto_N_;
      float ElecVeto_Pt_;
      float ElecVeto_Eta_;
      float ElecVeto_Phi_;

   private:
      const reco::GsfElectron* bestEle_;
      const pat::Tau* recoTauCand_;
      
};

class TauInfoContainer {
  public:
    // Default needed for persistency
    TauInfoContainer(); 
    ~TauInfoContainer();

    //TauInfoContainer( const pat::Tau* recoTauCand, const pat::Tau* altTauObj, std::vector<const reco::Candidate*>* trigObj, const reco::Candidate* GenParticle, unsigned int index, unsigned int nTotalObjects, const GenEventInfoProduct* GenInfo, unsigned int NVTX, const edm::Event* evt, const reco::Candidate* pfJet, const reco::Vertex* Vertex );

    //TauInfoContainer( const pat::Tau* recoTauCand, const pat::Tau* altTauObj, std::vector<const reco::Candidate*>* trigObj, const reco::Candidate* GenParticle, unsigned int index, unsigned int nTotalObjects, unsigned int NVTX, const edm::Event* evt, const reco::Candidate* pfJet, const reco::Vertex* Vertex );
    //
    TauInfoContainer( const pat::Tau* recoTauCand, const pat::Tau* altTauObj, std::vector<const reco::Candidate*>* trigObj, const reco::Candidate* GenParticle, unsigned int index, unsigned int nTotalObjects, const GenEventInfoProduct* GenInfo, unsigned int NVTX, unsigned int* evtInfo, const reco::Candidate* pfJet, const reco::Vertex* Vertex, EleMVAContainer * eleMVA, const reco::Candidate* genJet  );
    
    // Get tag tau object
    const pat::Tau* recoTauCand() const;

    const EleMVAContainer* getEleMVA() const;

    const pat::Tau* altTauObj() const;

    //const reco::Track* Track() const;

    const reco::Candidate* PfJet() const;

    const reco::Candidate* GenJet() const;

    const reco::Candidate* GenParticle() const; 

    bool isGenParticelMatched() const; 

    const reco::Vertex* getPV() const;

    // return true if pat::tau is matched to a hadronically decaying Gen Tau
    bool isTauGenJetMatched() const;

    bool isAltTauObjMatched() const;

    bool isPfJetMatched() const;
    bool isGenJetMatched() const;

    bool hasValidTrack() const;

    double pfCand_dz() const;

    // Get match status of trigger filter object
    bool isTrigObjMatched(int a) const;

    // Get status of Discriminator 
    double recoTauCandID(std::string DiscriminatorName) const;

    // Get the index of this match in the event.
    unsigned int index() const;
    // Get the total number of reco objects in this event.
    unsigned int nTotalObjects() const;

    const reco::Candidate* GenTauJet() const; 

    const GenEventInfoProduct* genInfo() const;

    int genDecayMode() const;

    int Nvtx() const;

    //const edm::Event* Evt() const;

    double RunNr() const;

    double EvtNr() const;

    double LumiSec() const;

    double TransImpPara() const;

    double TransImpParaError() const;

    //float Tau_GammaEnFrac() const;
    bool hasValidEleMVA() const{ return eleMVA_->hasValidEleMVA_; }
    float Tau_GSFTracklnPt() const{ return eleMVA_->Tau_GSFTracklnPt_; }
    float Tau_GSFTrackEta() const{ return eleMVA_->Tau_GSFTrackEta_; }
    float Tau_VisMass() const{ return eleMVA_->Tau_VisMass_; }
    float Tau_EmFraction() const{ return eleMVA_->Tau_EmFraction_; }
    float Tau_GSFTrackResol() const{ return eleMVA_->Tau_GSFTrackResol_; }
    float Tau_Pt() const{ return eleMVA_->Tau_Pt_; }
    float Tau_GammaEnFrac() const{ return eleMVA_->Tau_GammaEnFrac_; }
    float Tau_LeadChargedPFCandEtaAtEcalEntrance() const{ return eleMVA_->Tau_LeadChargedPFCandEtaAtEcalEntrance_; }
    float Tau_LeadChargedPFCandPt() const{ return eleMVA_->Tau_LeadChargedPFCandPt_; }
    int Tau_NumGammaCands() const{ return eleMVA_->Tau_NumGammaCands_; }
    float Tau_HadrHoP() const{ return eleMVA_->Tau_HadrHoP_; }
    float Tau_GammaEtaMom() const{ return eleMVA_->Tau_GammaEtaMom_; }
    float Tau_GammaPhiMom() const{ return eleMVA_->Tau_GammaPhiMom_; }
    float Tau_HadrEoP() const{ return eleMVA_->Tau_HadrEoP_; }
    float Tau_EtaAtEcalEntrance() const{ return eleMVA_->Tau_EtaAtEcalEntrance_; }
    float Tau_GSFChi2() const{ return eleMVA_->Tau_GSFChi2_; }
    int Tau_GSFNumHits() const{ return eleMVA_->Tau_GSFNumHits_; }
    int Tau_KFNumHits() const{ return eleMVA_->Tau_KFNumHits_; }
    float Tau_HadrMva() const{ return eleMVA_->Tau_HadrMva_; }

    bool Tau_HasGsf() const;


    float Elec_GSFTrackEta() const{ return eleMVA_->Elec_GSFTrackEta_; }
    float Elec_EtotOverPin() const{ return eleMVA_->Elec_EtotOverPin_; }
    float Elec_Fbrem() const{ return eleMVA_->Elec_Fbrem_; }
    float Elec_GSFNumHits() const{ return eleMVA_->Elec_GSFNumHits_; }
    float Elec_Chi2GSF() const{ return eleMVA_->Elec_Chi2GSF_; }
    float Elec_GSFTracklnPt() const{ return eleMVA_->Elec_GSFTracklnPt_; }
    float Elec_EgammaOverPdif() const{ return eleMVA_->Elec_EgammaOverPdif_; }
    float Elec_GSFTrackResol() const{ return eleMVA_->Elec_GSFTrackResol_; }

  private:

    const pat::Tau* recoTauCand_;
    const pat::Tau* altTauObj_;
    std::vector<const reco::Candidate*>* trigObj_;
    const reco::Candidate* GenParticle_;
    reco::Candidate* dummyCandidate_;
    pat::Tau* dummyCandidateTau_;
    unsigned int index_;
    unsigned int nTotalObjects_;
    const GenEventInfoProduct* genInfo_; 
    unsigned int Nvtx_;
    //const edm::Event* Evt_;
    unsigned int* evtInfo_;
    const reco::Candidate* pfJet_;
    const reco::Vertex* Vertex_;
    const EleMVAContainer* eleMVA_;
    const reco::Candidate* genJet_;
};

#endif /* end of include guard: TauInfoContainer _h */
