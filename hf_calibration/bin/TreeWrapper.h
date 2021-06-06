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
    } auxInfo;
    
    // primary vertices
    struct PrimVertex {
        unsigned int nTracks;
        double chi2, normalizedChi2, ndof, rho;
        TVector3 v3;
    };
    std::vector<PrimVertex> primVertices;
    
    // electrons
    struct Electron {
        int charge;
        TLorentzVector p4;
        double etaSC, dEtaIn, dPhiIn, sigmaIEtaIEta, full5x5SigmaIEtaIEta;
        double HoE, d0, dZ, IoEmIoP;
        double iso;
        bool hasConversion;
        unsigned int missingHits;
    };
    std::vector<Electron> electrons;
    
    // HF electrons
    struct HfElectron {
        TLorentzVector p4;
        double e1x1, e3x3, e5x5, eLong1x1, eLong3x3, eLong5x5, eShort1x1, eShort3x3, eShort5x5;
        double e9e25, eCOREe9, eCore, eSeL;
    };
    std::vector<HfElectron> hfElectrons;
    
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
    
    // event information
    unsigned int ev_run, ev_lumi, ev_nr, ev_bunch;
    
    // primary vertex information
    unsigned int pv_size, pv_nTracks[99];
    float pv_chi2[99], pv_normalizedChi2[99], pv_ndof[99], pv_x[99], pv_y[99], pv_z[99], pv_rho[99];
    
    // electron information
    unsigned int ele_size;
    int ele_charge[99];
    float ele_pt[99], ele_phi[99], ele_eta[99], ele_energy[99], ele_etaSC[99];
    float ele_dEtaIn[99], ele_dPhiIn[99], ele_sigmaIEtaIEta[99], ele_full5x5SigmaIEtaIEta[99];
    float ele_HoE[99], ele_d0[99], ele_dZ[99], ele_IoEmIoP[99];
    float ele_iso[99];
    bool ele_hasConversion[99];
    unsigned int ele_missingHits[99];
    
    // HF electron information
    unsigned int hfEle_size;
    float hfEle_pt[99], hfEle_phi[99], hfEle_eta[99], hfEle_energy[99];
    float hfEle_e1x1[99], hfEle_e3x3[99], hfEle_e5x5[99];
    float hfEle_eLong1x1[99], hfEle_eLong3x3[99], hfEle_eLong5x5[99];
    float hfEle_eShort1x1[99], hfEle_eShort3x3[99], hfEle_eShort5x5[99];
    float hfEle_e9e25[99], hfEle_eCOREe9[99], hfEle_eCore[99], hfEle_eSeL[99];
    
};

#endif
