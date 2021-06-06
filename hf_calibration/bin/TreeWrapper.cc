#include "HCALCalibration/HFCalibration/bin/TreeWrapper.h"

// constructor -------------------------------------------------------------------------------------
TreeWrapper::TreeWrapper(std::vector<std::string> &fileNames) : event() {
    chain = new TChain("Events");
    for (unsigned int i = 0; i != fileNames.size(); ++i) chain->Add(fileNames.at(i).c_str());
    chain->SetBranchAddress("ev_run", &ev_run);
    chain->SetBranchAddress("ev_lumi", &ev_lumi);
    chain->SetBranchAddress("ev_nr", &ev_nr);
    chain->SetBranchAddress("ev_bunch", &ev_bunch);
    chain->SetBranchAddress("pv_size", &pv_size);
    chain->SetBranchAddress("pv_nTracks", &pv_nTracks);
    chain->SetBranchAddress("pv_chi2", &pv_chi2);
    chain->SetBranchAddress("pv_normalizedChi2", &pv_normalizedChi2);
    chain->SetBranchAddress("pv_ndof", &pv_ndof);
    chain->SetBranchAddress("pv_x", &pv_x);
    chain->SetBranchAddress("pv_y", &pv_y);
    chain->SetBranchAddress("pv_z", &pv_z);
    chain->SetBranchAddress("pv_rho", &pv_rho);
    chain->SetBranchAddress("ele_size", &ele_size);
    chain->SetBranchAddress("ele_charge", &ele_charge);
    chain->SetBranchAddress("ele_pt", &ele_pt);
    chain->SetBranchAddress("ele_phi", &ele_phi);
    chain->SetBranchAddress("ele_eta", &ele_eta);
    chain->SetBranchAddress("ele_energy", &ele_energy);
    chain->SetBranchAddress("ele_etaSC", &ele_etaSC);
    chain->SetBranchAddress("ele_dEtaIn", &ele_dEtaIn);
    chain->SetBranchAddress("ele_dPhiIn", &ele_dPhiIn);
    chain->SetBranchAddress("ele_sigmaIEtaIEta", &ele_sigmaIEtaIEta);
    chain->SetBranchAddress("ele_full5x5SigmaIEtaIEta", &ele_full5x5SigmaIEtaIEta);
    chain->SetBranchAddress("ele_HoE", &ele_HoE);
    chain->SetBranchAddress("ele_d0", &ele_d0);
    chain->SetBranchAddress("ele_dZ", &ele_dZ);
    chain->SetBranchAddress("ele_IoEmIoP", &ele_IoEmIoP);
    chain->SetBranchAddress("ele_iso", &ele_iso);
    chain->SetBranchAddress("ele_hasConversion", &ele_hasConversion);
    chain->SetBranchAddress("ele_missingHits", &ele_missingHits);
    chain->SetBranchAddress("hfEle_size", &hfEle_size);
    chain->SetBranchAddress("hfEle_pt", &hfEle_pt);
    chain->SetBranchAddress("hfEle_phi", &hfEle_phi);
    chain->SetBranchAddress("hfEle_eta", &hfEle_eta);
    chain->SetBranchAddress("hfEle_energy", &hfEle_energy);
    chain->SetBranchAddress("hfEle_e1x1", &hfEle_e1x1);
    chain->SetBranchAddress("hfEle_e3x3", &hfEle_e3x3);
    chain->SetBranchAddress("hfEle_e5x5", &hfEle_e5x5);
    chain->SetBranchAddress("hfEle_eLong1x1", &hfEle_eLong1x1);
    chain->SetBranchAddress("hfEle_eLong3x3", &hfEle_eLong3x3);
    chain->SetBranchAddress("hfEle_eLong5x5", &hfEle_eLong5x5);
    chain->SetBranchAddress("hfEle_eShort1x1", &hfEle_eShort1x1);
    chain->SetBranchAddress("hfEle_eShort3x3", &hfEle_eShort3x3);
    chain->SetBranchAddress("hfEle_eShort5x5", &hfEle_eShort5x5);
    chain->SetBranchAddress("hfEle_e9e25", &hfEle_e9e25);
    chain->SetBranchAddress("hfEle_eCOREe9", &hfEle_eCOREe9);
    chain->SetBranchAddress("hfEle_eCore", &hfEle_eCore);
    chain->SetBranchAddress("hfEle_eSeL", &hfEle_eSeL);
}

// destructor --------------------------------------------------------------------------------------
TreeWrapper::~TreeWrapper() {
    delete chain;
}

// get number of entries ---------------------------------------------------------------------------
unsigned int TreeWrapper::Size() const {
    return chain->GetEntries();
}

// read the nth entry ------------------------------------------------------------------------------
void TreeWrapper::Read(unsigned int n) {
    event.primVertices.clear();
    event.electrons.clear();
    event.hfElectrons.clear();
    chain->GetEntry(n);
    event.auxInfo.run = ev_run;
    event.auxInfo.lumi = ev_lumi;
    event.auxInfo.bunch = ev_bunch;
    event.auxInfo.event = ev_nr;
    for (unsigned int i = 0; i != pv_size && i != 99; ++i) {
        EventInfo::PrimVertex pv;
        pv.nTracks = pv_nTracks[i];
        pv.chi2 = pv_chi2[i];
        pv.normalizedChi2 = pv_normalizedChi2[i];
        pv.ndof = pv_ndof[i];
        pv.rho = pv_rho[i];
        pv.v3.SetXYZ(pv_x[i], pv_y[i], pv_z[i]);
        event.primVertices.push_back(pv);
    }
    for (unsigned int i = 0; i != ele_size && i != 99; ++i) {
        EventInfo::Electron ele;
        ele.charge = ele_charge[i];
        ele.p4.SetPtEtaPhiE(ele_pt[i], ele_eta[i], ele_phi[i], ele_energy[i]);
        ele.etaSC = ele_etaSC[i];
        ele.dEtaIn = ele_dEtaIn[i];
        ele.dPhiIn = ele_dPhiIn[i];
        ele.sigmaIEtaIEta = ele_sigmaIEtaIEta[i];
        ele.full5x5SigmaIEtaIEta = ele_full5x5SigmaIEtaIEta[i];
        ele.HoE = ele_HoE[i];
        ele.d0 = ele_d0[i];
        ele.dZ = ele_dZ[i];
        ele.IoEmIoP = ele_IoEmIoP[i];
        ele.iso = ele_iso[i];
        ele.hasConversion = ele_hasConversion[i];
        ele.missingHits = ele_missingHits[i];
        event.electrons.push_back(ele);
    }
    for (unsigned int i = 0; i != hfEle_size && i != 99; ++i) {
        EventInfo::HfElectron hfEle;
        hfEle.p4.SetPtEtaPhiE(hfEle_pt[i], hfEle_eta[i], hfEle_phi[i], hfEle_energy[i]);
        hfEle.e1x1 = hfEle_e1x1[i];
        hfEle.e3x3 = hfEle_e3x3[i];
        hfEle.e5x5 = hfEle_e5x5[i];
        hfEle.eLong1x1 = hfEle_eLong1x1[i];
        hfEle.eLong3x3 = hfEle_eLong3x3[i];
        hfEle.eLong5x5 = hfEle_eLong5x5[i];
        hfEle.eShort1x1 = hfEle_eShort1x1[i];
        hfEle.eShort3x3 = hfEle_eShort3x3[i];
        hfEle.eShort5x5 = hfEle_eShort5x5[i];
        hfEle.e9e25 = hfEle_e9e25[i];
        hfEle.eCOREe9 = hfEle_eCOREe9[i];
        hfEle.eCore = hfEle_eCore[i];
        hfEle.eSeL = hfEle_eSeL[i];
        event.hfElectrons.push_back(hfEle);
    }
    return;
}

// get event ---------------------------------------------------------------------------------------
EventInfo& TreeWrapper::Event() {
    return event;
}
