import numpy, os, ROOT, sys, uuid
from PlotStyle import *

# get covariance matrix ----------------------------------------------------------------------------
def CovMatrix(path, runType, ybin):
    if 'ele' in runType:
        data = 'DoubleElectron'
    elif 'hf' in runType:
        data = 'SingleElectron'
    ROOT.gSystem.Load(os.environ['CMSSW_BASE'] + '/src/RooUnfold/libRooUnfold.so')
    fileRESP = ROOT.TFile.Open(path + 'responses_unf_' + runType + '.root')
    fileHIST = ROOT.TFile.Open(path + 'histograms_' + runType + '_' + data + '.root')
    resp = fileRESP.Get('response_' + ybin)
    hist = fileHIST.Get('ll_massX_cosThetaY_' + ybin + '_noBkg')
    hist.RebinY(100)
    bayes = ROOT.RooUnfoldBayes(resp, hist, 4)
    covm = bayes.Ereco()
    canvas = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 100, GetW(), GetH())
    SetCanvas(canvas)
    canvas.SetMargin(0.12, 0.11, 0.15, 0.05)
    covm.Draw()
    covh = ROOT.gROOT.FindObject('TMatrixDBase')
    if ybin == '2p4y5':
        vlim = 9e7
    else:
        vlim = 9e5
    for x in range(covh.GetNbinsX()):
        for y in range(covh.GetNbinsY()):
            b = covh.GetBin(1+x, 1+y)
            val = covh.GetBinContent(b)
            if val < -vlim:
                covh.SetBinContent(b, -vlim)
            if val > vlim:
                covh.SetBinContent(b, vlim)
    covh.SetStats(0)
    SetAxes(covh)
    pn = 999
    ps = [0.0, 0.49, 0.5, 0.51, 1.0]
    pr = [0.3, 0.0 , 1.0, 1.0 , 0.9]
    pg = [0.8, 0.0 , 1.0, 0.0 , 0.8]
    pb = [0.9, 1.0 , 1.0, 0.0 , 0.3]
    ROOT.TColor.CreateGradientColorTable(len(ps), numpy.array(ps), numpy.array(pr),
                                         numpy.array(pg), numpy.array(pb), pn)
    ROOT.gStyle.SetNumberContours(pn)
    covh.GetZaxis().SetRangeUser(-vlim, vlim)
    ROOT.TGaxis.SetMaxDigits(3)
    covh.Draw('colz')
    if ybin == '2p4y5':
        bins = [0, 40, 76, 86, 96, 106, 120, 150, 320, 8000]
    else:
        bins = [0, 40, 50, 60, 76, 86, 96, 106, 120, 133, 150, 171, 200, 320, 500, 2000, 8000]
    nb = len(bins)-1
    for b in range(nb):
        covh.GetXaxis().SetBinLabel(b+1, str(bins[b])+'-'+str(bins[b+1]))
        covh.GetXaxis().SetBinLabel(b+1+nb, str(bins[b])+'-'+str(bins[b+1]))
        covh.GetYaxis().SetBinLabel(b+1, str(bins[b])+'-'+str(bins[b+1]))
        covh.GetYaxis().SetBinLabel(b+1+nb, str(bins[b])+'-'+str(bins[b+1]))
    covh.GetXaxis().LabelsOption('v')
    covh.GetYaxis().LabelsOption('h')
    lineH = ROOT.TLine(0, nb, 2*nb, nb)
    lineH.SetLineStyle(3)
    lineH.Draw('same')
    lineV = ROOT.TLine(nb, 0, nb, 2*nb)
    lineV.SetLineStyle(3)
    lineV.Draw('same')
    canvas.Update()
    raw_input('...')
    fileHIST.Close()
    fileRESP.Close()

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print
        print 'usage: ' + sys.argv[0] + ' path runType ybin'
        print '       path    = path to histograms and responses'
        print '       runType = ele, hf'
        print '       ybin    = 0y1, 1y1p25, 1p25y1p5, 1p5y2p4, 2p4y5'
        print
        sys.exit()
    CovMatrix(sys.argv[1], sys.argv[2], sys.argv[3])
