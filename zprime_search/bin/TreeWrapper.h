#ifndef TreeWrapper_GUARD
#define TreeWrapper_GUARD

#include "TChain.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <string>
#include <vector>

// event information -------------------------------------------------------------------------------
struct EventInfo {
    
    // auxiliary information
    struct AuxInfo {
        unsigned int run, lumi, bunch, event;
        double weightMC, weightPU, weightFSR;
        std::vector<double> weightPDF_ct10, weightPDF_mstw, weightPDF_nnpdf;
    } auxInfo;
    
    // MC particles
    struct Particle {
        unsigned int numberOfDaughters;
        int charge, status, pdgId, pdgIdMother, pdgIdGrandMother;
        TLorentzVector p4;
        TVector3 vtx;
    };
    std::vector<Particle> particles;
    
    // primary vertices
    struct PrimVertex {
        unsigned int nTracks;
        double chi2, normalizedChi2, ndof, rho;
        TVector3 v3;
    };
    std::vector<PrimVertex> primVertices;
    
    // triggers
    std::vector<bool> ele8, ele17, ele17ele8, ele23hft30, ele27hft15, ele32sc17m50, ele27, ele27wp80;
    
    // electrons
    struct Electron {
        bool ele8, ele17, ele17ele8_f1, ele17ele8_f2, ele23hft30, ele27hft15;
        bool ele32sc17m50_f1, ele32sc17m50_f2, ele32sc17m50_f3, ele27, ele27wp80;
        int charge;
        TLorentzVector p4;
        double calibError, etaSC, dEtaIn, dPhiIn, sigmaIEtaIEta, HoE, d0, dZ, IoEmIoP, e1x5, e2x5;
        double e5x5, isoCH, isoNH, isoEM, isoNeutralNoPU, iso, isoDetEM, isoDetHad1, isoDetTrk;
        bool hasConversion, isEcalDriven;
        unsigned int missingHits;
        bool goodVetoId, goodLooseId, goodMediumId, goodTightId, goodHeepId;
        double mvaTrig, mvaNonTrig;
        TVector3 vtx;
    };
    std::vector<Electron> electrons;
    
    // HF electrons
    struct HfElectron {
        bool ele23hft30, ele27hft15;
        TLorentzVector p4;
        double smrf;
        double e1x1, e3x3, e5x5, eLong1x1, eLong3x3, eLong5x5, eShort1x1, eShort3x3, eShort5x5;
        double e9e25, eCOREe9, eCore, eSeL;
        TVector3 vtx;
    };
    std::vector<HfElectron> hfElectrons;
    
    // PF jets
    struct PfJet {
        TLorentzVector p4;
        unsigned int nConstituents, chargedMultiplicity;
        double chargedEmEnergyFraction, chargedHadronEnergyFraction;
        double neutralEmEnergyFraction, neutralHadronEnergyFraction, jetArea;
        bool goodLooseId, goodMediumId, goodTightId;
        TVector3 vtx;
    };
    std::vector<PfJet> pfJets;
    
    // PF MET
    struct PfMet { double et, sumEt, phi; } pfMet;
    
};

// tree wrapper class ------------------------------------------------------------------------------
class TreeWrapper {
    
public:
    
    // constructor
    TreeWrapper(std::vector<std::string> &fileNames);
    
    // destructor
    ~TreeWrapper();
    
    // get number of entries
    unsigned int Size() const;
    
    // read the nth entry
    void Read(unsigned int n);
    
    // get event
    EventInfo& Event();
    
private:
    
    // private objects
    TChain *chain;
    EventInfo event;
    bool isMC, isSG;
    
    // event information
    unsigned int ev_run, ev_lumi, ev_nr, ev_bunch;
    float ev_weightPU, ev_weightFSR;
    float ev_weightPDF_ct10[53], ev_weightPDF_mstw[41], ev_weightPDF_nnpdf[101];
    
    // generated particles
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
    float hfEle_pt[99], hfEle_phi[99], hfEle_eta[99], hfEle_energy[99], hfEle_smrf[99];
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
