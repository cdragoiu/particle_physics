#include "MyAnalyses/CommonTools/interface/Utils.h"
#include "MyAnalyses/ZprimeSearch/bin/GenInfo.h"
#include "MyAnalyses/ZprimeSearch/bin/Zprime.h"
#include "TLorentzVector.h"
#include "TMath.h"

// constructor -------------------------------------------------------------------------------------
GenInfo::GenInfo(bool _doHist) : doHist(_doHist), ele(0), pos(0), Z(0), goodParticles(false) {
    if (doHist) BookHistograms();
}

// destructor --------------------------------------------------------------------------------------
GenInfo::~GenInfo() {
    if (doHist) {
        delete gen_ele_p4;
        delete gen_pos_p4;
        delete gen_Z_p4;
        delete gen_dPhi;
        delete gen_dEta;
        delete gen_dPt;
        delete gen_dE;
        delete gen_ele_etaX_phiY;
        delete gen_pos_etaX_phiY;
        delete gen_ll_massX_cosThetaY_0y1;
        delete gen_ll_massX_cosThetaY_1y1p25;
        delete gen_ll_massX_cosThetaY_1p25y1p5;
        delete gen_ll_massX_cosThetaY_1p5y2p4;
        delete gen_ll_massX_cosThetaY_2p4y5;
    }
}

// get electron ------------------------------------------------------------------------------------
const EventInfo::Particle& GenInfo::GetEle() const {
    return *ele;
}

// get positron ------------------------------------------------------------------------------------
const EventInfo::Particle& GenInfo::GetPos() const {
    return *pos;
}

// get Z boson -------------------------------------------------------------------------------------
const EventInfo::Particle& GenInfo::GetZ() const {
    return *Z;
}

// get status --------------------------------------------------------------------------------------
bool GenInfo::IsGood() const {
    return goodParticles;
}

// reset information -------------------------------------------------------------------------------
void GenInfo::Reset() {
    ele = 0;
    pos = 0;
    Z = 0;
    goodParticles = false;
    return;
}

// select particles --------------------------------------------------------------------------------
void GenInfo::SelectParticles(const std::vector<EventInfo::Particle> &particles) {
    Reset();
    bool goodEle(false), goodPos(false), goodZ(false);
    for (unsigned int i = 0; i != particles.size(); ++i) {
        if (!goodEle && particles.at(i).pdgId == 11 && particles.at(i).pdgIdMother == 23) {
            ele = &particles.at(i);
            goodEle = true;
        }
        if (!goodPos && particles.at(i).pdgId == -11 && particles.at(i).pdgIdMother == 23) {
            pos = &particles.at(i);
            goodPos = true;
        }
        if (!goodZ && particles.at(i).pdgId == 23) {
            Z = &particles.at(i);
            goodZ = true;
        }
        if (goodEle && goodPos && goodZ) {
            goodParticles = true;
            break;
        }
    }
    return;
}

// fill raw histograms -----------------------------------------------------------------------------
void GenInfo::FillRawHistograms(double weight) {
    TLorentzVector ll = ele->p4 + pos->p4;
    double cosTheta = Zprime::GetCosTheta(ele->p4, pos->p4);
    if (fabs(ll.Rapidity()) < 1.0) {
        gen_ll_massX_cosThetaY_0y1->Fill(ll.M(), cosTheta, weight);
    }
    else if (fabs(ll.Rapidity()) < 1.25) {
        gen_ll_massX_cosThetaY_1y1p25->Fill(ll.M(), cosTheta, weight);
    }
    else if (fabs(ll.Rapidity()) < 1.5) {
        gen_ll_massX_cosThetaY_1p25y1p5->Fill(ll.M(), cosTheta, weight);
    }
    else if (fabs(ll.Rapidity()) < 2.4) {
        gen_ll_massX_cosThetaY_1p5y2p4->Fill(ll.M(), cosTheta, weight);
    }
    else if (fabs(ll.Rapidity()) < 5.0) {
        gen_ll_massX_cosThetaY_2p4y5->Fill(ll.M(), cosTheta, weight);
    }
    return;
}

// fill histograms ---------------------------------------------------------------------------------
void GenInfo::FillHistograms(double weight) {
    gen_ele_p4->Fill(ele->p4, weight);
    gen_pos_p4->Fill(pos->p4, weight);
    gen_Z_p4->Fill(Z->p4, weight);
    gen_dPhi->Fill(ele->p4.DeltaPhi(pos->p4), weight);
    gen_dEta->Fill(ele->p4.Eta() - pos->p4.Eta(), weight);
    gen_dPt->Fill(ele->p4.Pt() / pos->p4.Pt(), weight);
    gen_dE->Fill(ele->p4.Energy() / pos->p4.Energy(), weight);
    gen_ele_etaX_phiY->Fill(ele->p4.Eta(), ele->p4.Phi(), weight);
    gen_pos_etaX_phiY->Fill(pos->p4.Eta(), pos->p4.Phi(), weight);
    return;
}

// write histograms --------------------------------------------------------------------------------
void GenInfo::Write() {
    gen_ele_p4->Write();
    gen_pos_p4->Write();
    gen_Z_p4->Write();
    gen_dPhi->Write();
    gen_dEta->Write();
    gen_dPt->Write();
    gen_dE->Write();
    gen_ele_etaX_phiY->Write();
    gen_pos_etaX_phiY->Write();
    gen_ll_massX_cosThetaY_0y1->Write();
    gen_ll_massX_cosThetaY_1y1p25->Write();
    gen_ll_massX_cosThetaY_1p25y1p5->Write();
    gen_ll_massX_cosThetaY_1p5y2p4->Write();
    gen_ll_massX_cosThetaY_2p4y5->Write();
    return;
}

// book histograms ---------------------------------------------------------------------------------
void GenInfo::BookHistograms() {
    gen_ele_p4 = new P4Hists("gen_ele");
    gen_pos_p4 = new P4Hists("gen_pos");
    gen_Z_p4 = new P4Hists("gen_Z");
    gen_dPhi = Book1D("gen_dPhi", 400, -TMath::Pi(), TMath::Pi());
    gen_dEta = Book1D("gen_dEta", 400, -10.0, 10.0);
    gen_dPt = Book1D("gen_dPt", 400, 0.0, 4.0);
    gen_dE = Book1D("gen_dE", 400, 0.0, 4.0);
    gen_ele_etaX_phiY = Book2D("gen_ele_etaX_phiY", 200, -5.0, 5.0, 200, -TMath::Pi(), TMath::Pi());
    gen_pos_etaX_phiY = Book2D("gen_pos_etaX_phiY", 200, -5.0, 5.0, 200, -TMath::Pi(), TMath::Pi());
    gen_ll_massX_cosThetaY_0y1 = new TH2D("gen_ll_massX_cosThetaY_0y1", "",
                                          16, Zprime::massBins_ele, 200, -1.0, 1.0);
    gen_ll_massX_cosThetaY_0y1->Sumw2();
    gen_ll_massX_cosThetaY_1y1p25 = new TH2D("gen_ll_massX_cosThetaY_1y1p25", "",
                                             16, Zprime::massBins_ele, 200, -1.0, 1.0);
    gen_ll_massX_cosThetaY_1y1p25->Sumw2();
    gen_ll_massX_cosThetaY_1p25y1p5 = new TH2D("gen_ll_massX_cosThetaY_1p25y1p5", "",
                                               16, Zprime::massBins_ele, 200, -1.0, 1.0);
    gen_ll_massX_cosThetaY_1p25y1p5->Sumw2();
    gen_ll_massX_cosThetaY_1p5y2p4 = new TH2D("gen_ll_massX_cosThetaY_1p5y2p4", "",
                                              16, Zprime::massBins_ele, 200, -1.0, 1.0);
    gen_ll_massX_cosThetaY_1p5y2p4->Sumw2();
    gen_ll_massX_cosThetaY_2p4y5 = new TH2D("gen_ll_massX_cosThetaY_2p4y5", "",
                                              9, Zprime::massBins_hf, 200, -1.0, 1.0);
    gen_ll_massX_cosThetaY_2p4y5->Sumw2();
    return;
}
