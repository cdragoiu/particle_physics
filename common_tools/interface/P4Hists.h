#ifndef P4Hists_GUARD
#define P4Hists_GUARD

#include "TH1.h"
#include "TLorentzVector.h"
#include <string>

class P4Hists {
    
public:
    
    // constructors
    P4Hists() {}
    P4Hists(std::string name);
    
    // destructor
    ~P4Hists();
    
    // fill histograms
    void Fill(const TLorentzVector &p4, double weight = 1.0);
    
    // write histograms
    void Write();
    
private:
    
    TH1D *mass, *energy, *et, *pt, *phi, *eta, *rapidity;
    
};

#endif
