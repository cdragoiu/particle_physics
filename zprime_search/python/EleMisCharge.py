import numpy, ROOT, sys, uuid
from PlotStyle import *

# plot electron charge fake rate -------------------------------------------------------------------
def PlotFakeRate(path):
    bins = [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 90.0, 120.0, 200.0]
    histName = 'ele1_Mcut_pt'
    fileDT_OS = ROOT.TFile.Open(path + 'histograms_OS_DoubleElectron.root')
    histDT_OS = fileDT_OS.Get(histName)
    histDT_OS.SetName(str(uuid.uuid4()))
    fileDT_SS = ROOT.TFile.Open(path + 'histograms_SS_DoubleElectron.root')
    histDT_SS = fileDT_SS.Get(histName)
    histDT_SS.SetName(str(uuid.uuid4()))
    fileMC_OS = ROOT.TFile.Open(path + 'histograms_OS_DYtoEE.root')
    histMC_OS = fileMC_OS.Get(histName)
    histMC_OS.SetName(str(uuid.uuid4()))
    fileMC_SS = ROOT.TFile.Open(path + 'histograms_SS_DYtoEE.root')
    histMC_SS = fileMC_SS.Get(histName)
    histMC_SS.SetName(str(uuid.uuid4()))
    histDTr_OS = histDT_OS.Rebin(len(bins)-1, str(uuid.uuid4()), numpy.array(bins))
    histDTr_SS = histDT_SS.Rebin(len(bins)-1, str(uuid.uuid4()), numpy.array(bins))
    histMCr_OS = histMC_OS.Rebin(len(bins)-1, str(uuid.uuid4()), numpy.array(bins))
    histMCr_SS = histMC_SS.Rebin(len(bins)-1, str(uuid.uuid4()), numpy.array(bins))
    ratioDT = histDTr_OS.Clone(str(uuid.uuid4()))
    ratioDT.Divide(histDTr_SS, histDTr_OS, 1.0, 1.0)
    ratioDT.Scale(100)
    ratioMC = histMCr_OS.Clone(str(uuid.uuid4()))
    ratioMC.Divide(histMCr_SS, histMCr_OS, 1.0, 1.0)
    ratioMC.Scale(100)
    canvas = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 100, GetW(), GetH())
    SetCanvas(canvas)
    canvas.SetGridy()
    canvas.SetLogx()
    ratioMC.SetStats(0)
    SetAxes(ratioMC)
    ratioMC.GetXaxis().SetTitle('p_{T}')
    ratioMC.GetYaxis().SetTitle('Fakes (%)')
    ratioMC.GetYaxis().SetRangeUser(0.0, 1.3)
    ratioMC.SetLineColor(ROOT.kRed)
    ratioMC.SetMarkerColor(ROOT.kRed)
    ratioMC.SetMarkerStyle(22)
    ratioMC.SetMarkerSize(1.2)
    ratioMC.Draw('e p')
    ratioDT.SetLineColor(ROOT.kBlack)
    ratioDT.SetMarkerColor(ROOT.kBlack)
    ratioDT.SetMarkerStyle(23)
    ratioDT.SetMarkerSize(1.2)
    ratioDT.Draw('same e p')
    legend = ROOT.TLegend(0.76, 0.26, 0.88, 0.37)
    SetLegend(legend)
    legend.AddEntry(ratioDT, ' Data', 'p')
    legend.AddEntry(ratioMC, ' MC', 'p')
    legend.Draw()
    canvas.Update()
    raw_input('...')

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'usage: ' + sys.argv[0] + ' path'
        print '       path = path to histograms'
        sys.exit()
    PlotFakeRate(sys.argv[1])
