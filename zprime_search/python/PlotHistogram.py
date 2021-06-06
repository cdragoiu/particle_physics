import ROOT, sys, uuid
from PlotStyle import *
from SharedData import *

# plot histograms ----------------------------------------------------------------------------------
def Plot(data, mcs, drawP, tag):
    canvas = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 100, GetW(), GetH())
    SetCanvas(canvas)
    if 'x' in drawP.log:
        canvas.SetLogx()
    if 'y' in drawP.log:
        canvas.SetLogy()
    legend = ROOT.TLegend(0.76, 0.94-0.04*(1+len(mcs)), 0.88, 0.94)
    SetLegend(legend)
    fileDT = ROOT.TFile.Open(data.fileName)
    histDT = fileDT.Get(data.histName)
    histDT.SetName(str(uuid.uuid4()))
    if histDT.GetSumw2N() == 0:
        histDT.Sumw2()
    histDT.Rebin(drawP.bins)
    histDT.SetLineColor(data.color)
    histDT.SetMarkerColor(data.color)
    histDT.SetMarkerStyle(20)
    histDT.SetMarkerSize(1.0)
    histDT.GetXaxis().SetRangeUser(drawP.minX, drawP.maxX)
    histDT.GetYaxis().SetRangeUser(drawP.minY, drawP.maxY)
    legend.AddEntry(histDT, data.legendName, 'p')
    filesMC = []
    histsMC = []
    hstack = ROOT.THStack(str(uuid.uuid4()), '')
    for mc in mcs:
        filesMC.append(ROOT.TFile.Open(mc.fileName))
        histsMC.append(filesMC[-1].Get(mc.histName))
        histsMC[-1].SetName(str(uuid.uuid4()))
        if histsMC[-1].GetSumw2N() == 0:
            histsMC[-1].Sumw2()
        histsMC[-1].Rebin(drawP.bins)
        if 'QCD.root' in mc.fileName:
            histsMC[-1].Scale(QcdSF(tag))
        else:
            histsMC[-1].Scale(McSF(tag)*Lumi(tag))
        histsMC[-1].SetFillStyle(1001)
        histsMC[-1].SetFillColor(mc.color)
        histsMC[-1].SetLineColor(ROOT.kGray+3)
        histsMC[-1].GetXaxis().SetRangeUser(drawP.minX, drawP.maxX)
        histsMC[-1].GetYaxis().SetRangeUser(drawP.minY, drawP.maxY)
        legend.AddEntry(histsMC[-1], mc.legendName, 'f')
        hstack.Add(histsMC[-1])
    histTMP = histDT.Clone(str(uuid.uuid4()))
    histTMP.Reset()
    histTMP.SetStats(0)
    SetAxes(histTMP)
    histTMP.GetXaxis().SetTitle(drawP.titleX)
    histTMP.GetXaxis().SetRangeUser(drawP.minX, drawP.maxX)
    histTMP.GetYaxis().SetTitle(drawP.titleY)
    histTMP.GetYaxis().SetRangeUser(drawP.minY, drawP.maxY)
    histTMP.Draw()
    hstack.Draw('same hist')
    histDT.Draw('same e p')
    histTMP.Draw('same axis')
    legend.Draw()
    canvas.Update()
    histSUM = histDT.Clone(str(uuid.uuid4()))
    histSUM.Reset()
    for hist in histsMC:
        histSUM.Add(hist)
    ratio = histDT.Clone(str(uuid.uuid4()))
    ratio.Divide(histDT, histSUM)
    canvasR = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 130+GetH(), GetW(), GetH('R'))
    SetCanvas(canvasR, 'R')
    canvasR.SetGridy()
    if 'x' in drawP.log:
        canvasR.SetLogx()
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(1.0)
    ratio.SetStats(0)
    SetAxes(ratio, 'R')
    ratio.GetXaxis().SetTitle(drawP.titleX)
    ratio.GetXaxis().SetRangeUser(drawP.minX, drawP.maxX)
    ratio.GetYaxis().SetTitle('Data / MC')
    ratio.GetYaxis().SetRangeUser(0.5, 1.5)
    ratio.Draw('e p')
    canvasR.Update()
    raw_input('...')

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'usage: ' + sys.argv[0] + ' path histName runType [command]'
        print '       path     = path to histograms with prefix'
        print '       histName = name of histogram to plot'
        print '       runType  = ele, hf'
        print '       command  = noQCD'
        sys.exit()
    path = sys.argv[1]
    name = sys.argv[2]
    tag = sys.argv[3]
    if len(sys.argv) == 5 and sys.argv[4] == 'noQCD':
        addQCD = False
    else:
        addQCD = True
    if tag == 'ele':
        data = Hist(name, path + 'DoubleElectron.root', 'Data', ROOT.kBlack)
    elif tag == 'hf':
        data = Hist(name, path + 'SingleElectron.root', 'Data', ROOT.kBlack)
    mcs = []
    if tag == 'ele':
        mcs.append(Hist(name, path + 'WJets.root', ' W+jets', ROOT.kCyan-3))
        mcs.append(Hist(name, path + 'TW.root', ' tW,#bar{t}W', ROOT.kMagenta-3))
        mcs.append(Hist(name, path + 'DYtoTauTau.root', ' Z#rightarrow#tau#tau', ROOT.kYellow-3))
        if addQCD:
            mcs.append(Hist(name, path + 'QCD.root', ' jets', ROOT.kGray+1))
        mcs.append(Hist(name, path + 'TT.root', ' t#bar{t}', ROOT.kRed-3))
        mcs.append(Hist(name, path + 'VV.root', ' WW,WZ,ZZ', ROOT.kGreen-3))
        mcs.append(Hist(name, path + 'DYtoEE.root', ' Z#rightarrowee', ROOT.kAzure-3))
    elif tag == 'hf':
        mcs.append(Hist(name, path + 'TW.root', ' tW,#bar{t}W', ROOT.kMagenta-3))
        mcs.append(Hist(name, path + 'DYtoTauTau.root', ' Z#rightarrow#tau#tau', ROOT.kYellow-3))
        mcs.append(Hist(name, path + 'TT.root', ' t#bar{t}', ROOT.kRed-3))
        mcs.append(Hist(name, path + 'VV.root', ' WW,WZ,ZZ', ROOT.kGreen-3))
        mcs.append(Hist(name, path + 'WJets.root', ' W+jets', ROOT.kCyan-3))
        if addQCD:
            mcs.append(Hist(name, path + 'QCD.root', ' jets', ROOT.kGray+1))
        mcs.append(Hist(name, path + 'DYtoEE.root', ' Z#rightarrowee', ROOT.kAzure-3))
    Plot(data, mcs, GetDrawP(tag+'_'+name), tag)
