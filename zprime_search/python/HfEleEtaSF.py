import numpy, ROOT, sys, uuid
from PlotStyle import *

# calculate scale factors --------------------------------------------------------------------------
def ExtractSF(path, flag, plot):
    if flag == 'pos':
        bins = [3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 5.0]
    elif flag == 'neg':
        bins = [-5.0, -4.0, -3.8, -3.6, -3.4, -3.2, -3.0]
    histName = 'hfEle_Mcut_eta'
    fileD = ROOT.TFile.Open(path + 'histograms_hf_SingleElectron.root', 'READ')
    histRawD = fileD.Get(histName)
    histRawD.SetName(str(uuid.uuid4()))
    fileM = ROOT.TFile.Open(path + 'histograms_hf_DYtoEE.root', 'READ')
    histRawM = fileM.Get(histName)
    histRawM.SetName(str(uuid.uuid4()))
    histRawM.Scale(histRawD.Integral() / histRawM.Integral())
    histD = histRawD.Rebin(len(bins)-1, str(uuid.uuid4()), numpy.array(bins))
    histM = histRawM.Rebin(len(bins)-1, str(uuid.uuid4()), numpy.array(bins))
    ratio = histD.Clone(str(uuid.uuid4()))
    ratio.Reset()
    ratio.Divide(histD, histM, 1.0, 1.0)
    for r in range(1, 1+ratio.GetNbinsX()):
        print 'hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, ' + \
              str(ratio.GetBinLowEdge(r)) + ', ' + str(ratio.GetBinLowEdge(r+1)) + ', ' + \
              str(ratio.GetBinContent(r)) + ', ' + str(ratio.GetBinError(r)) + '));'
    if plot:
        canvas = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 100, GetW(), GetH())
        SetCanvas(canvas)
        if flag == 'pos':
            legend = ROOT.TLegend(0.76, 0.80, 0.88, 0.88)
        elif flag == 'neg':
            legend = ROOT.TLegend(0.16, 0.80, 0.28, 0.88)
        SetLegend(legend)
        histD.SetLineColor(ROOT.kBlack)
        histD.SetMarkerColor(ROOT.kBlack)
        histD.SetMarkerStyle(20)
        histD.SetMarkerSize(1.0)
        legend.AddEntry(histD, ' Data', 'p')
        histM.SetFillStyle(1001)
        histM.SetFillColor(ROOT.kAzure-3)
        histM.SetLineColor(ROOT.kGray+3)
        legend.AddEntry(histM, ' Z#rightarrowee', 'f')
        histM.SetStats(0)
        SetAxes(histM)
        histM.GetXaxis().SetTitle('#eta(HF)')
        histM.GetYaxis().SetRangeUser(1.0, 9900)
        histM.Draw('hist')
        histD.Draw('same e p')
        legend.Draw()
        canvas.Update()
        canvasR = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 130+GetH(), GetW(), GetH('R'))
        SetCanvas(canvasR, 'R')
        canvasR.SetGridy()
        ratio.SetLineColor(ROOT.kBlack)
        ratio.SetMarkerColor(ROOT.kBlack)
        ratio.SetMarkerStyle(20)
        ratio.SetMarkerSize(1.0)
        ratio.SetStats(0)
        SetAxes(ratio, 'R')
        ratio.GetXaxis().SetTitle('#eta(HF)')
        ratio.GetYaxis().SetTitle('Data / MC')
        ratio.GetYaxis().SetRangeUser(0.5, 1.5)
        ratio.Draw('e p')
        canvasR.Update()
        raw_input('...')

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'usage: ' + sys.argv[0] + ' path detSide [command]'
        print '       path    = path to histograms'
        print '       detSide = pos, neg'
        print '       command = plot'
        sys.exit()
    if len(sys.argv) is 4 and sys.argv[3] == 'plot':
        ExtractSF(sys.argv[1], sys.argv[2], True)
    else:
        ExtractSF(sys.argv[1], sys.argv[2], False)
