#ifndef TreeProducer_GUARD
#define TreeProducer_GUARD

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include <string>
#include <vector>

class TreeProducer : public edm::EDAnalyzer {
    
public:
    
    // constructor
	TreeProducer(const edm::ParameterSet &config);
    
    // destructor
	~TreeProducer();
    
private:
    
    // high-level trigger matching
    bool FindMatch(const trigger::TriggerEvent &triggerEvent, edm::InputTag &filterTag,
                   TLorentzVector &p4);
    
    // analyze event
	void analyze(const edm::Event &event, const edm::EventSetup &setup);
	
	// input variables
    bool isMC;
	std::string outFileName;
    edm::InputTag fsrWeightTag;
    edm::InputTag pdfWeightTag_ct10, pdfWeightTag_mstw, pdfWeightTag_nnpdf;
    edm::InputTag genParticleTag;
    edm::InputTag vertexTag;
	edm::InputTag triggerResultsTag, triggerSummaryTag;
    edm::InputTag electronTag;
    edm::InputTag isoRhoTag;
    edm::InputTag conversionTag;
    edm::InputTag beamSpotTag;
    edm::InputTag hfElectronTag;
    edm::InputTag hfClusterMapTag;
    edm::InputTag pfJetTag;
    edm::InputTag pfMetTag;
    
    // output variables
	TFile *outFile;
	TTree *tree;
    // event information
    unsigned int ev_run, ev_lumi, ev_nr, ev_bunch;
    float ev_weightPU, ev_weightFSR;
    float ev_weightPDF_ct10[53], ev_weightPDF_mstw[41], ev_weightPDF_nnpdf[101];
    // generated particles (e, Z, Z')
    unsigned int gen_size, gen_numberOfDaughters[99];
    int gen_charge[99], gen_status[99], gen_pdgId[99], gen_pdgIdMother[99], gen_pdgIdGrandMother[99];
    float gen_pt[99], gen_phi[99], gen_eta[99], gen_mass[99];
    float gen_vx[99], gen_vy[99], gen_vz[99];
    // primary vertex information
    unsigned int pv_size, pv_nTracks[99];
    float pv_chi2[99], pv_normalizedChi2[99], pv_ndof[99], pv_x[99], pv_y[99], pv_z[99], pv_rho[99];
    // trigger information
    bool hlt_ele8[7], hlt_ele17[7], hlt_ele17ele8[9];
    bool hlt_ele23hft30[9], hlt_ele27hft15[9], hlt_ele32sc17m50[7];
    bool hlt_ele27[8], hlt_ele27wp80[9];
    // electron information
    unsigned int ele_size;
    bool ele_ele8[99], ele_ele17[99], ele_ele17ele8_f1[99], ele_ele17ele8_f2[99];
    bool ele_ele23hft30[99], ele_ele27hft15[99];
    bool ele_ele32sc17m50_f1[99], ele_ele32sc17m50_f2[99], ele_ele32sc17m50_f3[99];
    bool ele_ele27[99], ele_ele27wp80[99];
    int ele_charge[99];
    float ele_pt[99], ele_phi[99], ele_eta[99], ele_energy[99], ele_calibError[99], ele_etaSC[99];
    float ele_dEtaIn[99], ele_dPhiIn[99], ele_sigmaIEtaIEta[99], ele_HoE[99];
    float ele_d0[99], ele_dZ[99], ele_IoEmIoP[99];
    float ele_e1x5[99], ele_e2x5[99], ele_e5x5[99];
    float ele_isoCH[99], ele_isoNH[99], ele_isoEM[99], ele_isoNeutralNoPU[99], ele_iso[99];
    float ele_isoDetEM[99], ele_isoDetHad1[99], ele_isoDetTrk[99];
    bool ele_hasConversion[99], ele_isEcalDriven[99];
    unsigned int ele_missingHits[99];
    bool ele_goodVetoId[99], ele_goodLooseId[99], ele_goodMediumId[99], ele_goodTightId[99];
    bool ele_goodHeepId[99];
    float ele_mvaTrig[99], ele_mvaNonTrig[99];
    float ele_vx[99], ele_vy[99], ele_vz[99];
    // HF electron information
    unsigned int hfEle_size;
    bool hfEle_ele23hft30[99], hfEle_ele27hft15[99];
    float hfEle_pt[99], hfEle_phi[99], hfEle_eta[99], hfEle_energy[99];
    float hfEle_e1x1[99], hfEle_e3x3[99], hfEle_e5x5[99];
    float hfEle_eLong1x1[99], hfEle_eLong3x3[99], hfEle_eLong5x5[99];
    float hfEle_eShort1x1[99], hfEle_eShort3x3[99], hfEle_eShort5x5[99];
    float hfEle_e9e25[99], hfEle_eCOREe9[99], hfEle_eCore[99], hfEle_eSeL[99];
    float hfEle_vx[99], hfEle_vy[99], hfEle_vz[99];
    // PF jet informaton
    unsigned int pfJet_size;
    float pfJet_pt[99], pfJet_phi[99], pfJet_eta[99], pfJet_energy[99];
    unsigned int pfJet_nConstituents[99], pfJet_chargedMultiplicity[99];
    float pfJet_chargedEmEnergyFraction[99], pfJet_chargedHadronEnergyFraction[99];
    float pfJet_neutralEmEnergyFraction[99], pfJet_neutralHadronEnergyFraction[99];
    bool pfJet_goodLooseId[99], pfJet_goodMediumId[99], pfJet_goodTightId[99];
    float pfJet_jetArea[99];
    float pfJet_vx[99], pfJet_vy[99], pfJet_vz[99];
    // PF MET information
    float pfMet_et, pfMet_sumEt, pfMet_phi;
    
};

#endif
