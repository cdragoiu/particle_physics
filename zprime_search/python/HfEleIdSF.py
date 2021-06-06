import numpy, ROOT, sys, uuid
from PlotStyle import *

# plot efficiency vs et/eta ------------------------------------------------------------------------
def PlotEff(fileName, runType):
    file = ROOT.TFile.Open(fileName, 'READ')
    if 'reco' in runType:
        histName = 'raw_hfEle_etX_etaY'
    elif 'gen' in runType:
        histName = 'hfEleGD_etX_etaY'
    hist_2D = file.Get(histName)
    histID_2D = file.Get(histName + '_ID')
    if 'eta' in runType:
        hist = hist_2D.ProjectionY('hist')
        histID = histID_2D.ProjectionY('histID')
        nb = 10
    elif 'et' in runType:
        hist = hist_2D.ProjectionX('hist')
        histID = histID_2D.ProjectionX('histID')
        nb = 5
    hist.Rebin(nb)
    histID.Rebin(nb)
    ratio = hist.Clone('ratio')
    ratio.Reset()
    ratio.Divide(histID, hist, 1.0, 1.0, 'B')
    ratio.Scale(100.0)
    canvas = ROOT.TCanvas('canvas', '', 440, 100, GetW(), GetH())
    SetCanvas(canvas)
    canvas.SetGridy()
    SetAxes(ratio)
    ratio.SetStats(0)
    ratio.SetLineColor(ROOT.kRed)
    ratio.SetMarkerColor(ROOT.kRed)
    ratio.SetMarkerStyle(20)
    if 'eta' in runType:
        ratio.GetXaxis().SetRangeUser(2.8, 5.0)
        ratio.GetXaxis().SetTitle('#eta')
    elif 'et' in runType:
        ratio.GetXaxis().SetRangeUser(25.0, 155.0)
        ratio.GetXaxis().SetTitle('E_{T}')
    ratio.GetYaxis().SetRangeUser(5.0, 105.0)
    ratio.GetYaxis().SetTitle('#varepsilon (%)')
    ratio.Draw('e p')
    canvas.Update()
    raw_input('...')

# calculate scale factors --------------------------------------------------------------------------
def ExtractSF(path):
    etBins = [20.0, 30.0, 40.0, 50.0, 70.0, 150.0]
    etaBins = [3.0, 3.5, 4.0, 5.0]
    histName = 'raw_hfEle_etX_etaY'
    fileD = ROOT.TFile.Open(path + 'histograms_hf_SingleElectron.root', 'READ')
    histD_2D = fileD.Get(histName)
    histD_2D.SetName(str(uuid.uuid4()))
    histD_ID_2D = fileD.Get(histName + '_ID')
    histD_ID_2D.SetName(str(uuid.uuid4()))
    fileM = ROOT.TFile.Open(path + 'histograms_hf_DYtoEE.root', 'READ')
    histM_2D = fileM.Get(histName)
    histM_2D.SetName(str(uuid.uuid4()))
    histM_ID_2D = fileM.Get(histName + '_ID')
    histM_ID_2D.SetName(str(uuid.uuid4()))
    for b in range(len(etBins)-1):
        histD_PY = histD_2D.ProjectionY(str(uuid.uuid4()), int(etBins[b])+1, int(etBins[b+1]))
        histD_ID_PY = histD_ID_2D.ProjectionY(str(uuid.uuid4()), int(etBins[b])+1, int(etBins[b+1]))
        histM_PY = histM_2D.ProjectionY(str(uuid.uuid4()), int(etBins[b])+1, int(etBins[b+1]))
        histM_ID_PY = histM_ID_2D.ProjectionY(str(uuid.uuid4()), int(etBins[b])+1, int(etBins[b+1]))
        histD = histD_PY.Rebin(len(etaBins)-1, str(uuid.uuid4()), numpy.array(etaBins))
        histD_ID = histD_ID_PY.Rebin(len(etaBins)-1, str(uuid.uuid4()), numpy.array(etaBins))
        histM = histM_PY.Rebin(len(etaBins)-1, str(uuid.uuid4()), numpy.array(etaBins))
        histM_ID = histM_ID_PY.Rebin(len(etaBins)-1, str(uuid.uuid4()), numpy.array(etaBins))
        ratioD = histD.Clone(str(uuid.uuid4()))
        ratioD.Reset()
        ratioD.Divide(histD_ID, histD, 1.0, 1.0, 'B')
        ratioM = histM.Clone(str(uuid.uuid4()))
        ratioM.Reset()
        ratioM.Divide(histM_ID, histM, 1.0, 1.0, 'B')
        ratio = ratioD.Clone(str(uuid.uuid4()))
        ratio.Reset()
        ratio.Divide(ratioD, ratioM, 1.0, 1.0)
        for r in range(1, ratio.GetNbinsX()+1):
            print 'hfEleIdSF.push_back(ScaleFactor(' + str(etBins[b]) + ', ' + str(etBins[b+1]) + \
                  ', ' + str(ratio.GetBinLowEdge(r)) + ', ' + str(ratio.GetBinLowEdge(r+1)) + \
                  ', ' + str(ratio.GetBinContent(r)) + ', ' + str(ratio.GetBinError(r)) + '));'

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) not in [2, 4]:
        print 'usage:'
        print '    get SF : ' + sys.argv[0] + ' path'
        print '             path = path to histograms'
        print '    plot SF: ' + sys.argv[0] + ' rootFile runType effType'
        print '             rootFile = input root file'
        print '             runType  = reco, gen'
        print '             effType  = et, eta'
        sys.exit()
    if len(sys.argv) is 2:
        ExtractSF(sys.argv[1])
    else:
        PlotEff(sys.argv[1], sys.argv[2] + '_' + sys.argv[3])
