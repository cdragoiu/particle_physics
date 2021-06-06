#include "HCALCalibration/HFCalibration/bin/Calibration.h"
#include "TMath.h"
#include <cmath>

// constructor -------------------------------------------------------------------------------------
Calibration::Calibration() : electron(0), hfElectron(0), goodElectrons(false) {
    BookHistograms();
}

// destructor --------------------------------------------------------------------------------------
Calibration::~Calibration() {
    delete ele_pt;
    delete ele_eta;
    delete ele_phi;
    delete hfEle_et;
    delete hfEle_eta;
    delete hfEle_phi;
    delete hfEle_eLong3x3;
    delete hfEle_eLong5x5;
    delete hfEle_eLCeL9;
    delete hfEle_eS9eL9;
    delete hfEle_eL9eL25;
    delete hfEle_cut2D;
    delete ll_pt;
    delete ll_mass;
    delete ll_rapidity;
    delete ll_phi;
    delete ele_Mcut_pt;
    delete ele_Mcut_eta;
    delete ele_Mcut_phi;
    delete hfEle_Mcut_et;
    delete hfEle_Mcut_eta;
    delete hfEle_Mcut_phi;
    delete hfEle_etaX_corrY;
}

// process event -----------------------------------------------------------------------------------
void Calibration::ProcessEvent(const EventInfo &event) {
    SelectElectrons(event.electrons, event.hfElectrons);
    if (!goodElectrons) return;
    FillHistograms();
    return;
}

// write histograms --------------------------------------------------------------------------------
void Calibration::Write() {
    ele_pt->Write();
    ele_eta->Write();
    ele_phi->Write();
    hfEle_et->Write();
    hfEle_eta->Write();
    hfEle_phi->Write();
    hfEle_eLong3x3->Write();
    hfEle_eLong5x5->Write();
    hfEle_eLCeL9->Write();
    hfEle_eS9eL9->Write();
    hfEle_eL9eL25->Write();
    hfEle_cut2D->Write();
    ll_pt->Write();
    ll_mass->Write();
    ll_rapidity->Write();
    ll_phi->Write();
    ele_Mcut_pt->Write();
    ele_Mcut_eta->Write();
    ele_Mcut_phi->Write();
    hfEle_Mcut_et->Write();
    hfEle_Mcut_eta->Write();
    hfEle_Mcut_phi->Write();
    hfEle_etaX_corrY->Write();
    return;
}

// select electrons --------------------------------------------------------------------------------
void Calibration::SelectElectrons(const std::vector<EventInfo::Electron> &electrons,
                                  const std::vector<EventInfo::HfElectron> &hfElectrons) {
    goodElectrons = false;
    std::vector<const EventInfo::Electron*> passElectrons;
    std::vector<const EventInfo::HfElectron*> passHfElectrons;
    for (unsigned int i = 0; i != electrons.size(); ++i) {
        if (!IsGoodEleP4(electrons.at(i).p4)) continue;
        passElectrons.push_back(&electrons.at(i));
    }
    for (unsigned int i = 0; i != hfElectrons.size(); ++i) {
        if (!IsGoodHfEleP4(hfElectrons.at(i).p4)) continue;
        passHfElectrons.push_back(&hfElectrons.at(i));
    }
    if (passElectrons.size() != 1 || passHfElectrons.size() != 1) return;
    if (!IsGoodEleID(*passElectrons.at(0)) || !IsGoodHfEleID(*passHfElectrons.at(0))) return;
    electron = passElectrons.at(0);
    hfElectron = passHfElectrons.at(0);
    goodElectrons = true;
    return;
}

// fill histograms ---------------------------------------------------------------------------------
void Calibration::FillHistograms() {
    ele_pt->Fill(electron->p4.Pt());
    ele_eta->Fill(electron->p4.Eta());
    ele_phi->Fill(electron->p4.Phi());
    hfEle_et->Fill(hfElectron->p4.Et());
    hfEle_eta->Fill(hfElectron->p4.Eta());
    hfEle_phi->Fill(hfElectron->p4.Phi());
    hfEle_eLong3x3->Fill(hfElectron->eLong3x3);
    hfEle_eLong5x5->Fill(hfElectron->eLong5x5);
    hfEle_eLCeL9->Fill(hfElectron->eCOREe9);
    hfEle_eS9eL9->Fill(hfElectron->eSeL);
    hfEle_eL9eL25->Fill(hfElectron->eLong3x3 / hfElectron->eLong5x5);
    hfEle_cut2D->Fill(hfElectron->eCOREe9 - 1.125 * hfElectron->eSeL);
    TLorentzVector ll = electron->p4 + hfElectron->p4;
    ll_pt->Fill(ll.Pt());
    ll_mass->Fill(ll.M());
    ll_rapidity->Fill(ll.Rapidity());
    ll_phi->Fill(ll.Phi());
    if (ll.M() < 50.0 || ll.M() > 90.0) return;
    double M = 91.1876;
    double emEta = electron->p4.Eta();
    double hfEta = hfElectron->p4.Eta();
    double emPhi = electron->p4.Phi();
    double hfPhi = hfElectron->p4.Phi();
    double emE = electron->p4.Energy();
    double predE = (pow(M,2)*cosh(emEta)*cosh(hfEta))/(2*emE*(cosh(emEta-hfEta)-cos(emPhi-hfPhi)));
    ele_Mcut_pt->Fill(electron->p4.Pt());
    ele_Mcut_eta->Fill(electron->p4.Eta());
    ele_Mcut_phi->Fill(electron->p4.Phi());
    hfEle_Mcut_et->Fill(hfElectron->p4.Et());
    hfEle_Mcut_eta->Fill(hfElectron->p4.Eta());
    hfEle_Mcut_phi->Fill(hfElectron->p4.Phi());
    hfEle_etaX_corrY->Fill(hfElectron->p4.Eta(), predE / hfElectron->p4.Energy());
    return;
}

// book histograms ---------------------------------------------------------------------------------
void Calibration::BookHistograms() {
    double pi = TMath::Pi();
    ele_pt = new TH1D("ele_pt", "", 1000, 0.0, 1000.0);
    ele_pt->Sumw2();
    ele_eta = new TH1D("ele_eta", "", 208, -5.2, 5.2);
    ele_eta->Sumw2();
    ele_phi = new TH1D("ele_phi", "", 200, -pi, pi);
    ele_phi->Sumw2();
    hfEle_et = new TH1D("hfEle_et", "", 4000, 0.0, 4000.0);
    hfEle_et->Sumw2();
    hfEle_eta = new TH1D("hfEle_eta", "", 208, -5.2, 5.2);
    hfEle_eta->Sumw2();
    hfEle_phi = new TH1D("hfEle_phi", "", 200, -pi, pi);
    hfEle_phi->Sumw2();
    hfEle_eLong3x3 = new TH1D("hfEle_eLong3x3", "", 2000, 0.0, 2000.0);
    hfEle_eLong3x3->Sumw2();
    hfEle_eLong5x5 = new TH1D("hfEle_eLong5x5", "", 2000, 0.0, 2000.0);
    hfEle_eLong5x5->Sumw2();
    hfEle_eLCeL9 = new TH1D("hfEle_eLCeL9", "", 1100, 0.0, 1.1);
    hfEle_eLCeL9->Sumw2();
    hfEle_eS9eL9 = new TH1D("hfEle_eS9eL9", "", 1100, 0.0, 1.1);
    hfEle_eS9eL9->Sumw2();
    hfEle_eL9eL25 = new TH1D("hfEle_eL9eL25", "", 1100, 0.0, 1.1);
    hfEle_eL9eL25->Sumw2();
    hfEle_cut2D = new TH1D("hfEle_cut2D", "", 1100, 0.0, 1.1);
    hfEle_cut2D->Sumw2();
    ll_pt = new TH1D("ll_pt", "", 1000, 0.0, 1000.0);
    ll_pt->Sumw2();
    ll_mass = new TH1D("ll_mass", "", 4000, 0.0, 4000.0);
    ll_mass->Sumw2();
    ll_rapidity = new TH1D("ll_rapidity", "", 208, -5.2, 5.2);
    ll_rapidity->Sumw2();
    ll_phi = new TH1D("ll_phi", "", 200, -pi, pi);
    ll_phi->Sumw2();
    ele_Mcut_pt = new TH1D("ele_Mcut_pt", "", 1000, 0.0, 1000.0);
    ele_Mcut_pt->Sumw2();
    ele_Mcut_eta = new TH1D("ele_Mcut_eta", "", 208, -5.2, 5.2);
    ele_Mcut_eta->Sumw2();
    ele_Mcut_phi = new TH1D("ele_Mcut_phi", "", 200, -pi, pi);
    ele_Mcut_phi->Sumw2();
    hfEle_Mcut_et = new TH1D("hfEle_Mcut_et", "", 4000, 0.0, 4000.0);
    hfEle_Mcut_et->Sumw2();
    hfEle_Mcut_eta = new TH1D("hfEle_Mcut_eta", "", 208, -5.2, 5.2);
    hfEle_Mcut_eta->Sumw2();
    hfEle_Mcut_phi = new TH1D("hfEle_Mcut_phi", "", 200, -pi, pi);
    hfEle_Mcut_phi->Sumw2();
    double bins[28] = {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013,
                       -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853,
                        2.853,  2.964,  3.139,  3.314,  3.489,  3.664,  3.839,
                        4.013,  4.191,  4.363,  4.538,  4.716,  4.889,  5.191};
    hfEle_etaX_corrY = new TH2D("hfEle_etaX_corrY", "", 27, bins, 200, 0.0, 2.0);
    hfEle_etaX_corrY->Sumw2();
    return;
}

// test electron kinematics ------------------------------------------------------------------------
bool Calibration::IsGoodEleP4(const TLorentzVector &p4) const {
    bool testP4 = (fabs(p4.Eta()) < 1.442 || fabs(p4.Eta()) > 1.556) && (p4.Pt() > 25.0);
    return testP4;
}

// test electron ID --------------------------------------------------------------------------------
bool Calibration::IsGoodEleID(const EventInfo::Electron &ele) const {
    bool testID = (
                   (fabs(ele.etaSC) < 1.442)     &&
                   (ele.sigmaIEtaIEta < 0.01)    &&
                   (fabs(ele.dEtaIn) < 0.008)    &&
                   (fabs(ele.dPhiIn) < 0.03)     &&
                   (ele.HoE < 0.06)              &&
                   (ele.iso < 0.10)              &&
                   (ele.IoEmIoP < 0.15)          &&
                   (fabs(ele.d0) < 0.01)         &&
                   (fabs(ele.dZ) < 0.07)         &&
                   (ele.missingHits <= 1)        &&
                   (!ele.hasConversion)
                  )                              ||
                  (
                   (fabs(ele.etaSC) > 1.566)     &&
                   (ele.sigmaIEtaIEta < 0.03)    &&
                   (fabs(ele.dEtaIn) < 0.009)    &&
                   (fabs(ele.dPhiIn) < 0.04)     &&
                   (ele.HoE < 0.10)              &&
                   (ele.iso < 0.12)              &&
                   (ele.IoEmIoP < 0.14)          &&
                   (fabs(ele.d0) < 0.05)         &&
                   (fabs(ele.dZ) < 0.18)         &&
                   (ele.missingHits <= 1)        &&
                   (!ele.hasConversion)
                  );
    return testID;
}

// test HF electron kinematics ---------------------------------------------------------------------
bool Calibration::IsGoodHfEleP4(const TLorentzVector &p4) const {
    bool testP4 = p4.Et() > 10.0;
    return testP4;
}

// test HF electron ID -----------------------------------------------------------------------------
bool Calibration::IsGoodHfEleID(const EventInfo::HfElectron &hfEle) const {
    bool testID = ((hfEle.eLong3x3 / hfEle.eLong5x5) > 0.96) &&
                  ((hfEle.eCOREe9 - 1.125 * hfEle.eSeL) > 0.4);
    return testID;
}
