#include "MyAnalyses/ZprimeSearch/bin/TreeWrapper.h"

// constructor -------------------------------------------------------------------------------------
TreeWrapper::TreeWrapper(std::vector<std::string> &fileNames) : event(), isMC(false), isSG(false) {
    chain = new TChain("Events");
    for (unsigned int i = 0; i != fileNames.size(); ++i) chain->Add(fileNames.at(i).c_str());
    chain->SetBranchAddress("ev_run", &ev_run);
    chain->SetBranchAddress("ev_lumi", &ev_lumi);
    chain->SetBranchAddress("ev_nr", &ev_nr);
    chain->SetBranchAddress("ev_bunch", &ev_bunch);
    chain->SetBranchAddress("ev_weightPU", &ev_weightPU);
    if (chain->GetBranch("ev_weightFSR")) {
        isSG = true;
        chain->SetBranchAddress("ev_weightFSR", &ev_weightFSR);
        chain->SetBranchAddress("ev_weightPDF_ct10", &ev_weightPDF_ct10);
        chain->SetBranchAddress("ev_weightPDF_mstw", &ev_weightPDF_mstw);
        chain->SetBranchAddress("ev_weightPDF_nnpdf", &ev_weightPDF_nnpdf);
    }
    if (chain->GetBranch("gen_size")) {
        isMC = true;
        chain->SetBranchAddress("gen_size", &gen_size);
        chain->SetBranchAddress("gen_numberOfDaughters", &gen_numberOfDaughters);
        chain->SetBranchAddress("gen_charge", &gen_charge);
        chain->SetBranchAddress("gen_status", &gen_status);
        chain->SetBranchAddress("gen_pdgId", &gen_pdgId);
        chain->SetBranchAddress("gen_pdgIdMother", &gen_pdgIdMother);
        chain->SetBranchAddress("gen_pdgIdGrandMother", &gen_pdgIdGrandMother);
        chain->SetBranchAddress("gen_pt", &gen_pt);
        chain->SetBranchAddress("gen_phi", &gen_phi);
        chain->SetBranchAddress("gen_eta", &gen_eta);
        chain->SetBranchAddress("gen_mass", &gen_mass);
        chain->SetBranchAddress("gen_vx", &gen_vx);
        chain->SetBranchAddress("gen_vy", &gen_vy);
        chain->SetBranchAddress("gen_vz", &gen_vz);
    }
    chain->SetBranchAddress("pv_size", &pv_size);
    chain->SetBranchAddress("pv_nTracks", &pv_nTracks);
    chain->SetBranchAddress("pv_chi2", &pv_chi2);
    chain->SetBranchAddress("pv_normalizedChi2", &pv_normalizedChi2);
    chain->SetBranchAddress("pv_ndof", &pv_ndof);
    chain->SetBranchAddress("pv_x", &pv_x);
    chain->SetBranchAddress("pv_y", &pv_y);
    chain->SetBranchAddress("pv_z", &pv_z);
    chain->SetBranchAddress("pv_rho", &pv_rho);
    chain->SetBranchAddress("hlt_ele8", &hlt_ele8);
    chain->SetBranchAddress("hlt_ele17", &hlt_ele17);
    chain->SetBranchAddress("hlt_ele17ele8", &hlt_ele17ele8);
    chain->SetBranchAddress("hlt_ele23hft30", &hlt_ele23hft30);
    chain->SetBranchAddress("hlt_ele27hft15", &hlt_ele27hft15);
    chain->SetBranchAddress("hlt_ele32sc17m50", &hlt_ele32sc17m50);
    chain->SetBranchAddress("hlt_ele27", &hlt_ele27);
    chain->SetBranchAddress("hlt_ele27wp80", &hlt_ele27wp80);
    chain->SetBranchAddress("ele_size", &ele_size);
    chain->SetBranchAddress("ele_ele8", &ele_ele8);
    chain->SetBranchAddress("ele_ele17", &ele_ele17);
    chain->SetBranchAddress("ele_ele17ele8_f1", &ele_ele17ele8_f1);
    chain->SetBranchAddress("ele_ele17ele8_f2", &ele_ele17ele8_f2);
    chain->SetBranchAddress("ele_ele23hft30", &ele_ele23hft30);
    chain->SetBranchAddress("ele_ele27hft15", &ele_ele27hft15);
    chain->SetBranchAddress("ele_ele32sc17m50_f1", &ele_ele32sc17m50_f1);
    chain->SetBranchAddress("ele_ele32sc17m50_f2", &ele_ele32sc17m50_f2);
    chain->SetBranchAddress("ele_ele32sc17m50_f3", &ele_ele32sc17m50_f3);
    chain->SetBranchAddress("ele_ele27", &ele_ele27);
    chain->SetBranchAddress("ele_ele27wp80", &ele_ele27wp80);
    chain->SetBranchAddress("ele_charge", &ele_charge);
    chain->SetBranchAddress("ele_pt", &ele_pt);
    chain->SetBranchAddress("ele_phi", &ele_phi);
    chain->SetBranchAddress("ele_eta", &ele_eta);
    chain->SetBranchAddress("ele_energy", &ele_energy);
    chain->SetBranchAddress("ele_calibError", &ele_calibError);
    chain->SetBranchAddress("ele_etaSC", &ele_etaSC);
    chain->SetBranchAddress("ele_dEtaIn", &ele_dEtaIn);
    chain->SetBranchAddress("ele_dPhiIn", &ele_dPhiIn);
    chain->SetBranchAddress("ele_sigmaIEtaIEta", &ele_sigmaIEtaIEta);
    chain->SetBranchAddress("ele_HoE", &ele_HoE);
    chain->SetBranchAddress("ele_d0", &ele_d0);
    chain->SetBranchAddress("ele_dZ", &ele_dZ);
    chain->SetBranchAddress("ele_IoEmIoP", &ele_IoEmIoP);
    chain->SetBranchAddress("ele_e1x5", &ele_e1x5);
    chain->SetBranchAddress("ele_e2x5", &ele_e2x5);
    chain->SetBranchAddress("ele_e5x5", &ele_e5x5);
    chain->SetBranchAddress("ele_isoCH", &ele_isoCH);
    chain->SetBranchAddress("ele_isoNH", &ele_isoNH);
    chain->SetBranchAddress("ele_isoEM", &ele_isoEM);
    chain->SetBranchAddress("ele_isoNeutralNoPU", &ele_isoNeutralNoPU);
    chain->SetBranchAddress("ele_iso", &ele_iso);
    chain->SetBranchAddress("ele_isoDetEM", &ele_isoDetEM);
    chain->SetBranchAddress("ele_isoDetHad1", &ele_isoDetHad1);
    chain->SetBranchAddress("ele_isoDetTrk", &ele_isoDetTrk);
    chain->SetBranchAddress("ele_hasConversion", &ele_hasConversion);
    chain->SetBranchAddress("ele_isEcalDriven", &ele_isEcalDriven);
    chain->SetBranchAddress("ele_missingHits", &ele_missingHits);
    chain->SetBranchAddress("ele_goodVetoId", &ele_goodVetoId);
    chain->SetBranchAddress("ele_goodLooseId", &ele_goodLooseId);
    chain->SetBranchAddress("ele_goodMediumId", &ele_goodMediumId);
    chain->SetBranchAddress("ele_goodTightId", &ele_goodTightId);
    chain->SetBranchAddress("ele_goodHeepId", &ele_goodHeepId);
    chain->SetBranchAddress("ele_mvaTrig", &ele_mvaTrig);
    chain->SetBranchAddress("ele_mvaNonTrig", &ele_mvaNonTrig);
    chain->SetBranchAddress("ele_vx", &ele_vx);
    chain->SetBranchAddress("ele_vy", &ele_vy);
    chain->SetBranchAddress("ele_vz", &ele_vz);
    chain->SetBranchAddress("hfEle_size", &hfEle_size);
    chain->SetBranchAddress("hfEle_ele23hft30", &hfEle_ele23hft30);
    chain->SetBranchAddress("hfEle_ele27hft15", &hfEle_ele27hft15);
    chain->SetBranchAddress("hfEle_pt", &hfEle_pt);
    chain->SetBranchAddress("hfEle_phi", &hfEle_phi);
    chain->SetBranchAddress("hfEle_eta", &hfEle_eta);
    chain->SetBranchAddress("hfEle_energy", &hfEle_energy);
    if (chain->GetBranch("hfEle_smrf")) chain->SetBranchAddress("hfEle_smrf", &hfEle_smrf);
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
    chain->SetBranchAddress("hfEle_vx", &hfEle_vx);
    chain->SetBranchAddress("hfEle_vy", &hfEle_vy);
    chain->SetBranchAddress("hfEle_vz", &hfEle_vz);
    chain->SetBranchAddress("pfJet_size", &pfJet_size);
    chain->SetBranchAddress("pfJet_pt", &pfJet_pt);
    chain->SetBranchAddress("pfJet_phi", &pfJet_phi);
    chain->SetBranchAddress("pfJet_eta", &pfJet_eta);
    chain->SetBranchAddress("pfJet_energy", &pfJet_energy);
    chain->SetBranchAddress("pfJet_nConstituents", &pfJet_nConstituents);
    chain->SetBranchAddress("pfJet_chargedMultiplicity", &pfJet_chargedMultiplicity);
    chain->SetBranchAddress("pfJet_chargedEmEnergyFraction", &pfJet_chargedEmEnergyFraction);
    chain->SetBranchAddress("pfJet_chargedHadronEnergyFraction", &pfJet_chargedHadronEnergyFraction);
    chain->SetBranchAddress("pfJet_neutralEmEnergyFraction", &pfJet_neutralEmEnergyFraction);
    chain->SetBranchAddress("pfJet_neutralHadronEnergyFraction", &pfJet_neutralHadronEnergyFraction);
    chain->SetBranchAddress("pfJet_goodLooseId", &pfJet_goodLooseId);
    chain->SetBranchAddress("pfJet_goodMediumId", &pfJet_goodMediumId);
    chain->SetBranchAddress("pfJet_goodTightId", &pfJet_goodTightId);
    chain->SetBranchAddress("pfJet_jetArea", &pfJet_jetArea);
    chain->SetBranchAddress("pfJet_vx", &pfJet_vx);
    chain->SetBranchAddress("pfJet_vy", &pfJet_vy);
    chain->SetBranchAddress("pfJet_vz", &pfJet_vz);
    chain->SetBranchAddress("pfMet_et", &pfMet_et);
    chain->SetBranchAddress("pfMet_sumEt", &pfMet_sumEt);
    chain->SetBranchAddress("pfMet_phi", &pfMet_phi);
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
    event.auxInfo.weightPDF_ct10.clear();
    event.auxInfo.weightPDF_mstw.clear();
    event.auxInfo.weightPDF_nnpdf.clear();
    event.particles.clear();
    event.primVertices.clear();
    event.ele8.clear();
    event.ele17.clear();
    event.ele17ele8.clear();
    event.ele23hft30.clear();
    event.ele27hft15.clear();
    event.ele32sc17m50.clear();
    event.ele27.clear();
    event.ele27wp80.clear();
    event.electrons.clear();
    event.hfElectrons.clear();
    event.pfJets.clear();
    chain->GetEntry(n);
    event.auxInfo.run = ev_run;
    event.auxInfo.lumi = ev_lumi;
    event.auxInfo.bunch = ev_bunch;
    event.auxInfo.event = ev_nr;
    event.auxInfo.weightMC = chain->GetWeight();
    event.auxInfo.weightPU = ev_weightPU;
    if (isSG) {
        event.auxInfo.weightFSR = ev_weightFSR;
        for (unsigned int i = 0; i != 53; ++i) {
            event.auxInfo.weightPDF_ct10.push_back(ev_weightPDF_ct10[i]);
        }
        for (unsigned int i = 0; i != 41; ++i) {
            event.auxInfo.weightPDF_mstw.push_back(ev_weightPDF_mstw[i]);
        }
        for (unsigned int i = 0; i != 101; ++i) {
            event.auxInfo.weightPDF_nnpdf.push_back(ev_weightPDF_nnpdf[i]);
        }
    }
    if (isMC) {
        for (unsigned int i = 0; i != gen_size && i != 99; ++i) {
            EventInfo::Particle gen;
            gen.numberOfDaughters = gen_numberOfDaughters[i];
            gen.charge = gen_charge[i];
            gen.status = gen_status[i];
            gen.pdgId = gen_pdgId[i];
            gen.pdgIdMother = gen_pdgIdMother[i];
            gen.pdgIdGrandMother = gen_pdgIdGrandMother[i];
            gen.p4.SetPtEtaPhiM(gen_pt[i], gen_eta[i], gen_phi[i], gen_mass[i]);
            gen.vtx.SetXYZ(gen_vx[i], gen_vy[i], gen_vz[i]);
            event.particles.push_back(gen);
        }
    }
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
    for (unsigned int i = 0; i != 7; ++i) event.ele8.push_back(hlt_ele8[i]);
    for (unsigned int i = 0; i != 7; ++i) event.ele17.push_back(hlt_ele17[i]);
    for (unsigned int i = 0; i != 9; ++i) event.ele17ele8.push_back(hlt_ele17ele8[i]);
    for (unsigned int i = 0; i != 9; ++i) event.ele23hft30.push_back(hlt_ele23hft30[i]);
    for (unsigned int i = 0; i != 9; ++i) event.ele27hft15.push_back(hlt_ele27hft15[i]);
    for (unsigned int i = 0; i != 7; ++i) event.ele32sc17m50.push_back(hlt_ele32sc17m50[i]);
    for (unsigned int i = 0; i != 8; ++i) event.ele27.push_back(hlt_ele27[i]);
    for (unsigned int i = 0; i != 9; ++i) event.ele27wp80.push_back(hlt_ele27wp80[i]);
    for (unsigned int i = 0; i != ele_size && i != 99; ++i) {
        EventInfo::Electron ele;
        ele.ele8 = ele_ele8[i];
        ele.ele17 = ele_ele17[i];
        ele.ele17ele8_f1 = ele_ele17ele8_f1[i];
        ele.ele17ele8_f2 = ele_ele17ele8_f2[i];
        ele.ele23hft30 = ele_ele23hft30[i];
        ele.ele27hft15 = ele_ele27hft15[i];
        ele.ele32sc17m50_f1 = ele_ele32sc17m50_f1[i];
        ele.ele32sc17m50_f2 = ele_ele32sc17m50_f2[i];
        ele.ele32sc17m50_f3 = ele_ele32sc17m50_f3[i];
        ele.ele27 = ele_ele27[i];
        ele.ele27wp80 = ele_ele27wp80[i];
        ele.charge = ele_charge[i];
        ele.p4.SetPtEtaPhiE(ele_pt[i], ele_eta[i], ele_phi[i], ele_energy[i]);
        ele.calibError = ele_calibError[i];
        ele.etaSC = ele_etaSC[i];
        ele.dEtaIn = ele_dEtaIn[i];
        ele.dPhiIn = ele_dPhiIn[i];
        ele.sigmaIEtaIEta = ele_sigmaIEtaIEta[i];
        ele.HoE = ele_HoE[i];
        ele.d0 = ele_d0[i];
        ele.dZ = ele_dZ[i];
        ele.IoEmIoP = ele_IoEmIoP[i];
        ele.e1x5 = ele_e1x5[i];
        ele.e2x5 = ele_e2x5[i];
        ele.e5x5 = ele_e5x5[i];
        ele.isoCH = ele_isoCH[i];
        ele.isoNH = ele_isoNH[i];
        ele.isoEM = ele_isoEM[i];
        ele.isoNeutralNoPU = ele_isoNeutralNoPU[i];
        ele.iso = ele_iso[i];
        ele.isoDetEM = ele_isoDetEM[i];
        ele.isoDetHad1 = ele_isoDetHad1[i];
        ele.isoDetTrk = ele_isoDetTrk[i];
        ele.hasConversion = ele_hasConversion[i];
        ele.isEcalDriven = ele_isEcalDriven[i];
        ele.missingHits = ele_missingHits[i];
        ele.goodVetoId = ele_goodVetoId[i];
        ele.goodLooseId = ele_goodLooseId[i];
        ele.goodMediumId = ele_goodMediumId[i];
        ele.goodTightId = ele_goodTightId[i];
        ele.goodHeepId = ele_goodHeepId[i];
        ele.mvaTrig = ele_mvaTrig[i];
        ele.mvaNonTrig = ele_mvaNonTrig[i];
        ele.vtx.SetXYZ(ele_vx[i], ele_vy[i], ele_vz[i]);
        event.electrons.push_back(ele);
    }
    for (unsigned int i = 0; i != hfEle_size && i != 99; ++i) {
        EventInfo::HfElectron hfEle;
        hfEle.ele23hft30 = hfEle_ele23hft30[i];
        hfEle.ele27hft15 = hfEle_ele27hft15[i];
        hfEle.p4.SetPtEtaPhiE(hfEle_pt[i], hfEle_eta[i], hfEle_phi[i], hfEle_energy[i]);
        if (isMC) hfEle.smrf = hfEle_smrf[i];
        else hfEle.smrf = 1.0;
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
        hfEle.vtx.SetXYZ(hfEle_vx[i], hfEle_vy[i], hfEle_vz[i]);
        event.hfElectrons.push_back(hfEle);
    }
    for (unsigned int i = 0; i != pfJet_size && i != 99; ++i) {
        EventInfo::PfJet pfJet;
        pfJet.p4.SetPtEtaPhiE(pfJet_pt[i], pfJet_eta[i], pfJet_phi[i], pfJet_energy[i]);
        pfJet.nConstituents = pfJet_nConstituents[i];
        pfJet.chargedMultiplicity = pfJet_chargedMultiplicity[i];
        pfJet.chargedEmEnergyFraction = pfJet_chargedEmEnergyFraction[i];
        pfJet.chargedHadronEnergyFraction = pfJet_chargedHadronEnergyFraction[i];
        pfJet.neutralEmEnergyFraction = pfJet_neutralEmEnergyFraction[i];
        pfJet.neutralHadronEnergyFraction = pfJet_neutralHadronEnergyFraction[i];
        pfJet.jetArea = pfJet_jetArea[i];
        pfJet.goodLooseId = pfJet_goodLooseId[i];
        pfJet.goodMediumId = pfJet_goodMediumId[i];
        pfJet.goodTightId = pfJet_goodTightId[i];
        pfJet.vtx.SetXYZ(pfJet_vx[i], pfJet_vy[i], pfJet_vz[i]);
        event.pfJets.push_back(pfJet);
    }
    event.pfMet.et = pfMet_et;
    event.pfMet.sumEt = pfMet_sumEt;
    event.pfMet.phi = pfMet_phi;
    return;
}

// get event ---------------------------------------------------------------------------------------
EventInfo& TreeWrapper::Event() {
    return event;
}
