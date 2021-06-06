#include "MyAnalyses/CommonTools/interface/P4Hists.h"
#include "MyAnalyses/CommonTools/interface/Utils.h"
#include "TMath.h"

// constructor -------------------------------------------------------------------------------------
P4Hists::P4Hists(std::string name) {
    mass = Book1D(name + "_mass", 4000, 0.0, 4000.0);
    energy = Book1D(name + "_energy", 4000, 0.0, 4000.0);
    et = Book1D(name + "_et", 4000, 0.0, 4000.0);
    pt = Book1D(name + "_pt", 4000, 0.0, 4000.0);
    phi = Book1D(name + "_phi", 200, -TMath::Pi(), TMath::Pi());
    eta = Book1D(name + "_eta", 200, -5.0, 5.0);
    rapidity = Book1D(name + "_rapidity", 200, -5.0, 5.0);
}

// destructor --------------------------------------------------------------------------------------
P4Hists::~P4Hists() {
    delete mass;
    delete energy;
    delete et;
    delete pt;
    delete phi;
    delete eta;
    delete rapidity;
}

// fill histograms ---------------------------------------------------------------------------------
void P4Hists::Fill(const TLorentzVector &p4, double weight) {
    mass->Fill(p4.M(), weight);
    energy->Fill(p4.Energy(), weight);
    et->Fill(p4.Et(), weight);
    pt->Fill(p4.Pt(), weight);
    phi->Fill(p4.Phi(), weight);
    eta->Fill(p4.Eta(), weight);
    rapidity->Fill(p4.Rapidity(), weight);
    return;
}

// write histograms --------------------------------------------------------------------------------
void P4Hists::Write() {
    mass->Write();
    energy->Write();
    et->Write();
    pt->Write();
    phi->Write();
    eta->Write();
    rapidity->Write();
    return;
}
