import math, numpy, os, random, ROOT
from PlotStyle import *
from SharedData import *

# check forward/backward statistics ----------------------------------------------------------------
def CheckStat():
    file = ROOT.TFile.Open('root/results_ele/histograms_ele_DoubleElectron.root')
    # (0y1, 1y1p25, 1p25y1p5, 1p5y2p4) (_noBkg, noBkg_unf)
    hist = file.Get('ll_massX_cosThetaY_1p5y2p4_noBkg_unf')
    hist.RebinY(100)
    histB = hist.ProjectionX('B', 1, 1, 'e')
    histF = hist.ProjectionX('F', 2, 2, 'e')
    print
    print '-------------------------------'
    print '    Mass          N_B       N_F'
    print '-------------------------------'
    for b in range(2, histB.GetNbinsX()):
        mMin = histB.GetBinLowEdge(b)
        mMax = histB.GetBinLowEdge(b+1)
        nB = histB.GetBinContent(b)
        nF = histF.GetBinContent(b)
        print '{:>4n} - {:<4n}   {:>7.0f}   {:>7.0f}'.format(mMin, mMax, nB, nF)
    print '-------------------------------'
    print

# check unfolding iteration ------------------------------------------------------------------------
def CheckUnfIter():
    file4 = ROOT.TFile.Open('root/checks/unfIter/I4/histograms_ele_DoubleElectron.root')
    hist4 = file4.Get('afb_0y1_noBkg_unf')
    hist4.SetName('I4')
    file2 = ROOT.TFile.Open('root/checks/unfIter/I2/histograms_ele_DoubleElectron.root')
    hist2 = file2.Get('afb_0y1_noBkg_unf')
    hist2.SetName('I2')
    file3 = ROOT.TFile.Open('root/checks/unfIter/I3/histograms_ele_DoubleElectron.root')
    hist3 = file3.Get('afb_0y1_noBkg_unf')
    hist3.SetName('I3')
    file5 = ROOT.TFile.Open('root/checks/unfIter/I5/histograms_ele_DoubleElectron.root')
    hist5 = file5.Get('afb_0y1_noBkg_unf')
    hist5.SetName('I5')
    hist = hist5 # <---
    canvas = ROOT.TCanvas('canvas')
    canvas.SetLogx()
    hist4.SetStats(0)
    hist4.GetXaxis().SetTitle('M')
    hist4.GetYaxis().SetTitle('A_{FB}')
    hist4.GetXaxis().SetRangeUser(41.0, 1999.0)
    hist4.SetLineColor(ROOT.kBlack)
    hist4.SetMarkerStyle(24)
    hist4.SetMarkerSize(0.6)
    hist4.SetMarkerColor(ROOT.kBlack)
    hist4.Draw('e1 x0 p')
    hist.SetLineColor(ROOT.kRed)
    hist.SetMarkerStyle(7)
    hist.SetMarkerColor(ROOT.kRed)
    hist.Draw('same e1 x0 p')
    legend = ROOT.TLegend(0.15, 0.75, 0.25, 0.85)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.AddEntry(hist4, ' N = 4', 'l')
    legend.AddEntry(hist, ' N = ' + hist.GetName(), 'l')
    legend.Draw()
    canvas.Update()
    raw_input('...')

# get acceptance correction ------------------------------------------------------------------------
def GetAccCorr():
    # 0y1, 1y1p25, 1p25y1p5, 1p5y2p4, 2p4y5
    yb = '0y1'
    tag = 'ele'
    fileA = ROOT.TFile.Open('root/checks/accCorr/histograms_' + tag + '_DYtoEE.root')
    fileB = ROOT.TFile.Open('root/results_' + tag + '/histograms_' + tag + '_DYtoEE.root')
    histA = fileA.Get('gen_afb_' + yb)
    histA.SetName('A')
    histB = fileB.Get('gen_afb_' + yb)
    histB.SetName('B')
    ratio = histA.Clone('ratio')
    ratio.Divide(histA, histB, 1.0, 1.0)
    print
    print '-------------------------------'
    print '    Mass               Corr'
    print '-------------------------------'
    for b in range(2, ratio.GetNbinsX()):
        mMin = ratio.GetBinLowEdge(b)
        mMax = ratio.GetBinLowEdge(b+1)
        cor = ratio.GetBinContent(b)
        err = ratio.GetBinError(b)
        print '{:>4n} - {:<4n}       {:>1.2f} +/- {:<1.2f}'.format(mMin, mMax, cor, err)
    print '-------------------------------'
    print
    canvas = ROOT.TCanvas('canvas')
    canvas.SetLogx()
    histB.SetStats(0)
    histB.GetXaxis().SetTitle('M')
    histB.GetYaxis().SetTitle('A_{FB}')
    xMin = histB.GetBinLowEdge(2)
    xMax = histB.GetBinLowEdge(histB.GetNbinsX())
    histB.GetXaxis().SetRangeUser(xMin+1, xMax-1)
    yMin = min(histA.GetMinimum(), histB.GetMinimum())
    yMax = max(histA.GetMaximum(), histB.GetMaximum())
    histB.GetYaxis().SetRangeUser(1.2*yMin, 1.2*yMax)
    histB.SetLineColor(ROOT.kBlack)
    histB.SetMarkerStyle(20)
    histB.SetMarkerColor(ROOT.kBlack)
    histB.Draw('e x0 p')
    histA.SetLineColor(ROOT.kRed)
    histA.SetMarkerStyle(20)
    histA.SetMarkerColor(ROOT.kRed)
    histA.Draw('same e x0 p')
    legend = ROOT.TLegend(0.15, 0.75, 0.25, 0.85)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)
    legend.SetTextSize(0.04)
    legend.AddEntry(histB, ' gen (no cuts)', 'p')
    legend.AddEntry(histA, ' gen (with cuts)', 'p')
    legend.Draw()
    canvas.Update()
    raw_input('...')

# check unfold errors (pull) -----------------------------------------------------------------------
def CheckUnfErr():
    ROOT.gSystem.Load(os.environ['CMSSW_BASE'] + '/src/RooUnfold/libRooUnfold.so')
    histP = ROOT.TH1D('hP', '', 100, -5.0, 5.0)
    histP.Sumw2()
    for i in range(200):
        histG = ROOT.TH2D('hG_'+str(i), '', 12, 0.0, 12.0, 12, 0.0, 12.0)
        histG.Sumw2()
        histR = ROOT.TH2D('hR_'+str(i), '', 12, 0.0, 12.0, 12, 0.0, 12.0)
        histR.Sumw2()
        resp = ROOT.RooUnfoldResponse(histG, histG, 'resp')
        for j in range(2000):
            xg = 6.0 + math.pow(-1, j) * 0.5 * random.gauss(1, 0.4)
            yg = 6.0 + math.pow(-1, j) * 0.5 * random.gauss(1, 0.4)
            histG.Fill(xg, yg)
        for j in range(2000):
            xg = 6.0 + math.pow(-1, j) * 0.5 * random.gauss(1, 0.4)
            yg = 6.0 + math.pow(-1, j) * 0.5 * random.gauss(1, 0.4)
            xr = xg * random.gauss(1, 0.2) + 0.2
            yr = yg * random.gauss(1, 0.2) - 0.2
            histR.Fill(xr, yr)
        for j in range(2000):
            xg = 6.0 + math.pow(-1, j) * 0.5 * random.gauss(1, 0.4)
            yg = 6.0 + math.pow(-1, j) * 0.5 * random.gauss(1, 0.4)
            xr = xg * random.gauss(1, 0.2) + 0.2
            yr = yg * random.gauss(1, 0.2) - 0.2
            resp.Fill(xr, yr, xg, yg)
        bayes = ROOT.RooUnfoldBayes(resp, histR, 4)
        histU = bayes.Hreco(1)
        b = histU.GetBin(6,6)
        pull = (histU.GetBinContent(b) - histG.GetBinContent(b)) / histU.GetBinError(b)
        histP.Fill(pull)
        hist = histU
        hist.SetStats(0)
        hist.SetMarkerSize(1.5)
        ROOT.gStyle.SetPaintTextFormat('4.0f')
#        hist.Draw('col text')
#        raw_input('...')
    histP.Draw('')
    raw_input('...')

# get electron charge fake rate --------------------------------------------------------------------
def GetChargeFakeRate():
    file = ROOT.TFile.Open('root/checks/chargeFakeRate/histograms_ele_DYtoEE.root')
    histA = file.Get('ele_mass_cfr_miss')
    histB = file.Get('ele_mass_cfr_all')
    ratio = histA.Clone('ratio')
    ratio.Divide(histA, histB, 1.0, 1.0, 'B')
    ratio.Scale(100)
    canvas = ROOT.TCanvas('canvas')
    canvas.SetLogx()
    canvas.SetGridy()
    ratio.SetStats(0)
    ratio.GetXaxis().SetTitle('M')
    ratio.GetYaxis().SetTitle('Charge Fake Rate (%)')
    ratio.GetXaxis().SetRangeUser(41.0, 1999.0)
    ratio.GetYaxis().SetRangeUser(0.0, 1.7)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.Draw('e x0 p')
    raw_input('...')

# get the HF MC overall scale factor uncertainty ---------------------------------------------------
def GetHfMcSfErr():
    fileD = ROOT.TFile.Open('root/checks/hfMcSfErr/histograms_hf_SingleElectron.root')
    fileM = ROOT.TFile.Open('root/checks/hfMcSfErr/histograms_hf_DYtoEE.root')
    histD_2D = fileD.Get('hfEleEtX_llMassY')
    histM_2D = fileM.Get('hfEleEtX_llMassY')
    bins = [20.0, 30.0, 40.0, 50.0, 70.0]
    hist = ROOT.TH1D('hist', '', len(bins)-1, numpy.array(bins))
    sumEV = 0.0
    sumE = 0.0
    for b in range(len(bins)-1):
        histD = histD_2D.ProjectionY('data', int(bins[b])+1, int(bins[b+1]), 'e')
        histM = histM_2D.ProjectionY('mc', int(bins[b])+1, int(bins[b+1]), 'e')
        histM.Scale(Lumi('hf'))
        ratio = histD.Clone('ratio')
        ratio.Divide(histD, histM, 1.0, 1.0)
        func = ROOT.TF1('func', '[0]', 80.0, 100.0)
        ratio.Fit('func', 'RQN')
        ratio.Fit('func', 'RQN')
        v = func.GetParameter(0)
        e = func.GetParError(0)
        hist.SetBinContent(b+1, v)
        hist.SetBinError(b+1, e)
        sumEV = sumEV + e * abs(McSF('hf') - v)
        sumE = sumE + e
    print 'err = ' + str(sumEV / sumE)
    canvas = ROOT.TCanvas('canvas', '', 440, 100, GetW(), GetH())
    SetCanvas(canvas)
    SetAxes(hist)
    hist.SetStats(0)
    hist.SetLineColor(ROOT.kBlack)
    hist.SetLineWidth(2)
    hist.GetXaxis().SetTitle('E_{T}(HF)')
    hist.GetYaxis().SetTitle('HF MC SF')
    hist.GetYaxis().SetRangeUser(0.21, 0.79)
    hist.Draw()
    line = ROOT.TLine(bins[0], McSF('hf'), bins[-1], McSF('hf'))
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(2)
    line.Draw('same')
    canvas.Update()
    raw_input('...')

# get inclusive cos(theta) distributions -----------------------------------------------------------
def GetCosTheta(mb1, mb2):
    path = 'root/results_ele/'
    names = ['histograms_ele_DYtoEE.root',
             'histograms_ele_DYtoTauTau.root',
             'histograms_ele_DoubleElectron.root',
             'histograms_ele_QCD.root',
             'histograms_ele_TT.root',
             'histograms_ele_TW.root',
             'histograms_ele_VV.root',
             'histograms_ele_WJets.root']
    for name in names:
        file = ROOT.TFile.Open(path + name)
        h1_2d = file.Get('ll_massX_cosThetaY_0y1')
        h2_2d = file.Get('ll_massX_cosThetaY_1y1p25')
        h3_2d = file.Get('ll_massX_cosThetaY_1p25y1p5')
        h4_2d = file.Get('ll_massX_cosThetaY_1p5y2p4')
        h1 = h1_2d.ProjectionY('h1', mb1, mb2, 'e')
        h2 = h2_2d.ProjectionY('h2', mb1, mb2, 'e')
        h3 = h3_2d.ProjectionY('h3', mb1, mb2, 'e')
        h4 = h4_2d.ProjectionY('h4', mb1, mb2, 'e')
        hist = h1.Clone('cosTheta')
        hist.Reset()
        hist.Add(h1, 1.0)
        hist.Add(h2, 1.0)
        hist.Add(h3, 1.0)
        hist.Add(h4, 1.0)
        ofile = ROOT.TFile.Open(name, 'recreate')
        ofile.cd()
        hist.Write()
        ofile.Close()
        file.Close()

#CheckStat()

#CheckUnfIter()

#GetAccCorr()

#CheckUnfErr()

#GetChargeFakeRate()

#GetHfMcSfErr()

#GetCosTheta(10, 10)
