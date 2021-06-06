import numpy, ROOT, sys, uuid
from PlotStyle import *

# extract trigger efficiency -----------------------------------------------------------------------
def GetEff(file, nameT, nameTnP, etaMin, etaMax, ptBins):
    histT = file.Get(nameT)
    histTnP = file.Get(nameTnP)
    histT_PX = histT.ProjectionX(str(uuid.uuid4()), 1+int(1000*etaMin), int(1000*etaMax), 'e')
    histTnP_PX = histTnP.ProjectionX(str(uuid.uuid4()), 1+int(1000*etaMin), int(1000*etaMax), 'e')
    histT_RB = histT_PX.Rebin(len(ptBins)-1, str(uuid.uuid4()), numpy.array(ptBins))
    histTnP_RB = histTnP_PX.Rebin(len(ptBins)-1, str(uuid.uuid4()), numpy.array(ptBins))
    ratio = histT_RB.Clone(str(uuid.uuid4()))
    ratio.Reset()
    ratio.Divide(histTnP_RB, histT_RB, 1.0, 1.0, 'B')
    ratio.Scale(100.0)
    return ratio

# extract trigger scale factors --------------------------------------------------------------------
def GetSF(path, runType, trgType):
    ptBins = [20.0, 30.0, 40.0, 50.0, 70.0, 250.0]
    etaBins = [0.0, 0.8, 1.442, 1.556, 2.0, 2.5]
    nameT = 'hltEle27_ptX_etaY_T'
    nameTnP = 'hltEle' + trgType + '_ptX_etaY_TnP'
    for e in range(len(etaBins)-1):
        fileDT = ROOT.TFile.Open(path + 'histograms_' + runType + '_SingleElectron.root')
        fileMC = ROOT.TFile.Open(path + 'histograms_' + runType + '_DYtoEE.root')
        effDT = GetEff(fileDT, nameT, nameTnP, etaBins[e], etaBins[e+1], ptBins)
        effMC = GetEff(fileMC, nameT, nameTnP, etaBins[e], etaBins[e+1], ptBins)
        ratio = effDT.Clone(str(uuid.uuid4()))
        ratio.Reset()
        ratio.Divide(effDT, effMC, 1.0, 1.0)
        for r in range(1, ratio.GetNbinsX()+1):
            print 'trgEle' + trgType + 'SF.push_back(ScaleFactor(' + \
                  str(ratio.GetBinLowEdge(r)) + ', ' + str(ratio.GetBinLowEdge(r+1)) + ', ' + \
                  str(etaBins[e]) + ', ' + str(etaBins[e+1]) + ', ' + \
                  str(ratio.GetBinContent(r)) + ', ' + str(ratio.GetBinError(r)) + '));'

# plot trigger efficiencies ------------------------------------------------------------------------
def PlotEff(path, runType, trgType):
    ptBins = [4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0,
              35.0, 40.0, 50.0, 60.0, 80.0, 100.0, 150.0, 200.0]
    etaBins = [0.0, 0.8, 1.442, 1.556, 2.0, 2.5]
    nameT = 'hltEle27_ptX_etaY_T'
    nameTnP = 'hltEle' + trgType + '_ptX_etaY_TnP'
    if trgType == '8': ptMin = 6.0
    elif trgType == '17': ptMin = 10.0
    elif trgType == '27': ptMin = 20.0
    for e in range(len(etaBins)-1):
        fileDT = ROOT.TFile.Open(path + 'histograms_' + runType + '_SingleElectron.root')
        fileMC = ROOT.TFile.Open(path + 'histograms_' + runType + '_DYtoEE.root')
        effDT = GetEff(fileDT, nameT, nameTnP, etaBins[e], etaBins[e+1], ptBins)
        effMC = GetEff(fileMC, nameT, nameTnP, etaBins[e], etaBins[e+1], ptBins)
        canvasName = 'canvas_' + str(etaBins[e]) + '_' + str(etaBins[e+1])
        canvas = ROOT.TCanvas(canvasName, '', 440, 100, GetW(), GetH())
        SetCanvas(canvas)
        canvas.SetLogx()
        canvas.SetGridy()
        SetAxes(effMC)
        effMC.SetStats(0)
        effMC.GetXaxis().SetTitle('p_{T}')
        effMC.GetXaxis().SetRangeUser(ptMin, -1.0)
        effMC.GetYaxis().SetTitle('#varepsilon (%)')
        effMC.GetYaxis().SetRangeUser(0.0, 105.0)
        effMC.SetLineColor(ROOT.kRed-3)
        effMC.SetMarkerColor(ROOT.kRed-3)
        effMC.SetMarkerStyle(21)
        effMC.Draw('e x0 p')
        effDT.SetLineColor(ROOT.kBlack)
        effDT.SetMarkerColor(ROOT.kBlack)
        effDT.SetMarkerStyle(20)
        effDT.Draw('same e x0 p')
        legend = ROOT.TLegend(0.80, 0.20, 0.92, 0.28)
        SetLegend(legend)
        legend.AddEntry(effDT, 'Data', 'p')
        legend.AddEntry(effMC, 'Z#rightarrowee', 'p')
        legend.Draw()
        canvas.Update()
        raw_input('...')

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 5:
        print 'usage: ' + sys.argv[0] + ' path command runType trigger'
        print '       path    = path to histograms'
        print '       command = print, plot'
        print '       runType = ele, hf'
        print '       trigger = 8, 17, 27'
        sys.exit()
    if 'print' in sys.argv[2]:
        GetSF(sys.argv[1], sys.argv[3], sys.argv[4])
    if 'plot' in sys.argv[2]:
        PlotEff(sys.argv[1], sys.argv[3], sys.argv[4])
