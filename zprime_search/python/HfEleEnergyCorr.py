import numpy, ROOT, sys, uuid
from PlotStyle import *

# get energy correction ----------------------------------------------------------------------------
def GetCorr(fileName, runType):
    file = ROOT.TFile.Open(fileName)
    if 'gen' in runType:
        hist2D = file.Get('hfEleGD_etX_corrY')
    else:
        hist2D = file.Get('hfEle_etX_corrY')
    bins = [20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0,
            40.0, 42.0, 44.0, 46.0, 48.0, 50.0, 55.0, 60.0, 70.0, 90.0, 120.0]
    hist = ROOT.TH1D(str(uuid.uuid4()), '', len(bins)-1, numpy.array(bins))
    for i in range(len(bins)-1):
        histPY = hist2D.ProjectionY(str(uuid.uuid4()), int(bins[i])+1, int(bins[i+1]))
        histPY.Rebin(4)
        fit = ROOT.TF1(str(uuid.uuid4()), 'gaus', 0.0, 2.0)
        histPY.Fit(fit, 'RQN')
        histPY.Fit(fit, 'RQN')
        hist.SetBinContent(i+1, fit.GetParameter(1))
        hist.SetBinError(i+1, fit.GetParError(1))
    if 'mc' in runType or 'gen' in runType:
        expr = '[0] - [1]*exp(-[2]*(x-[3])^2) + [4]*TMath::Erf(-[5]*(x-[6])^2)'
    else:
        expr = '[0]*exp(-[1]*(x-[2])^2) + [3]*TMath::Erf([4]*x-[5])'
    func = ROOT.TF1('func', expr, 20.0, 120.0)
    if 'mc' in runType or 'gen' in runType:
        func.SetParameters(1.3, 0.3, 0.001, 3.0, 1.0, 0.001, 10.0)
    else:
        func.SetParameters(1.0, 0.001, 1.0, 1.0, 0.01, 0.1)
    hist.Fit('func', 'RQN')
    hist.Fit('func', 'RQN')
    print
    print 'f(x): ' + expr
    print '       range: 20 <= Et <= 120 GeV'
    print '       Chi2/NDF = ' + str(func.GetChisquare()) + ' / ' + str(func.GetNDF())
    print '       [0] = ' + str(func.GetParameter(0)) + ' +/- ' + str(func.GetParError(0))
    print '       [1] = ' + str(func.GetParameter(1)) + ' +/- ' + str(func.GetParError(1))
    print '       [2] = ' + str(func.GetParameter(2)) + ' +/- ' + str(func.GetParError(2))
    print '       [3] = ' + str(func.GetParameter(3)) + ' +/- ' + str(func.GetParError(3))
    print '       [4] = ' + str(func.GetParameter(4)) + ' +/- ' + str(func.GetParError(4))
    print '       [5] = ' + str(func.GetParameter(5)) + ' +/- ' + str(func.GetParError(5))
    if 'mc' in runType or 'gen' in runType:
        print '       [6] = ' + str(func.GetParameter(6)) + ' +/- ' + str(func.GetParError(6))
    print
    canvas = ROOT.TCanvas('canvas', '', 440, 100, GetW(), GetH())
    SetCanvas(canvas)
    SetAxes(hist)
    hist.SetStats(0)
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(1.0)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetLineColor(ROOT.kBlack)
    hist.GetXaxis().SetTitle('E_{T}')
    hist.GetYaxis().SetTitle('E_{pred} / E_{reco}')
    hist.GetYaxis().SetRangeUser(0.76, 1.19)
    hist.Draw('e p')
    func.Draw('same')
    canvas.Update()
    raw_input('...')

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'usage: ' + sys.argv[0] + ' rootFile corrType'
        print '       rootFile = input root file'
        print '       corrType = data, mc, gen'
        sys.exit()
    GetCorr(sys.argv[1], sys.argv[2])
