import ROOT, sys, uuid
from PlotStyle import *
from SharedData import *
from QcdShape import *

# find scale factor --------------------------------------------------------------------------------
def FindSF(path, runType, nb, command='', sfPlot=0.0):
    histName = 'llRapidityX_llMassY'
    if 'ele' in runType:
        tag = 'ele'
        data = 'DoubleElectron'
        massMin = 41
        massMax = 60
    elif 'hf' in runType:
        tag = 'hf'
        data = 'SingleElectron'
        massMin = 121
        massMax = 300
    fileDT = ROOT.TFile.Open(path + 'histograms_' + runType + '_' + data + '.root')
    histDT_2D = fileDT.Get(histName)
    histDT = histDT_2D.ProjectionX(str(uuid.uuid4()), massMin, massMax, 'e')
    histDT.Rebin(nb)
    histMC = histDT.Clone(str(uuid.uuid4()))
    histMC.Reset()
    for mc in ['DYtoEE', 'DYtoTauTau', 'TT', 'TW', 'VV', 'WJets']:
        file = ROOT.TFile.Open(path + 'histograms_' + runType + '_' + mc + '.root')
        hist_2D = file.Get(histName)
        hist = hist_2D.ProjectionX(str(uuid.uuid4()), massMin, massMax, 'e')
        hist.Rebin(nb)
        hist.Scale(McSF(tag)*Lumi(tag))
        histMC.Add(hist)
    QcdShape(path, histName, runType)
    fileQCD = ROOT.TFile.Open(path + 'histograms_' + runType + '_QCD.root')
    histQCD_2D = fileQCD.Get(histName)
    histQCD = histQCD_2D.ProjectionX(str(uuid.uuid4()), massMin, massMax, 'e')
    histQCD.Rebin(nb)
    histTST = histDT.Clone(str(uuid.uuid4()))
    val = []
    err = []
    for i in range(1, 300):
        sf = 0.01 * i
        histTST.Reset()
        histTST.Add(histMC, histQCD, 1.0, sf)
        ratio = histDT.Clone(str(uuid.uuid4()))
        ratio.Divide(histDT, histTST, 1.0, 1.0)
        func = ROOT.TF1(str(uuid.uuid4()), '[0]', -5.0, 5.0)
        ratio.Fit(func, 'RQN')
        ratio.Fit(func, 'RQN')
        val.append(abs(1.0-func.GetParameter(0)))
        err.append(func.GetParError(0))
        if command == 'print':
            print 'SF = ' + str(sf) + '\t |1-val| = ' + str(val[-1]) + '\t err = ' + str(err[-1])
        elif command == 'plot' and abs(sf - sfPlot) < 0.001:
            canvas = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 100, GetW(), GetH())
            SetCanvas(canvas)
            canvas.SetLogy()
            legend = ROOT.TLegend(0.14, 0.78, 0.26, 0.94)
            SetLegend(legend)
            histDT.SetLineColor(ROOT.kBlack)
            histDT.SetMarkerColor(ROOT.kBlack)
            histDT.SetMarkerStyle(20)
            histDT.SetMarkerSize(1.0)
            legend.AddEntry(histDT, ' Data', 'p')
            histMC.SetFillStyle(1001)
            histMC.SetFillColor(ROOT.kAzure-3)
            histMC.SetLineColor(ROOT.kGray+3)
            legend.AddEntry(histMC, ' MC Sim', 'f')
            histQCD.Scale(sf)
            histQCD.SetFillStyle(1001)
            histQCD.SetFillColor(ROOT.kRed-3)
            histQCD.SetLineColor(ROOT.kGray+3)
            legend.AddEntry(histQCD, ' QCD', 'f')
            hstack = ROOT.THStack(str(uuid.uuid4()), '')
            hstack.Add(histQCD)
            hstack.Add(histMC)
            histTST.Reset()
            histTST.SetStats(0)
            SetAxes(histTST)
            histTST.GetXaxis().SetTitle('Y(ee)')
            histTST.GetYaxis().SetRangeUser(1, 1e4)
            histTST.Draw()
            hstack.Draw('same hist')
            histDT.Draw('same e p')
            histTST.Draw('same axis')
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
            ratio.GetXaxis().SetTitle('Y(ee)')
            ratio.GetYaxis().SetTitle('Data / MC')
            ratio.GetYaxis().SetRangeUser(0.5, 1.5)
            ratio.Draw('e p')
            canvasR.Update()
            raw_input('...')
            return
    k = 0
    for i in range(1, len(val)):
        if val[i] < val[k]:
            k = i
    kmin = k
    for i in range(k-1, 0, -1):
        if val[i] < val[k] + err[k]:
            continue
        kmin = i
        break
    kmax = k
    for i in range(k+1, len(val)):
        if val[i] < val[k] + err[k]:
            continue
        kmax = i
        break
    print
    sf = 0.01 * (kmin+1)
    print 'SF = ' + str(sf) + '\t |1-val| = ' + str(val[kmin]) + '\t err = ' + str(err[kmin])
    sf = 0.01 * (k+1)
    print 'SF = ' + str(sf) + '\t |1-val| = ' + str(val[k]) + '\t err = ' + str(err[k])
    sf = 0.01 * (kmax+1)
    print 'SF = ' + str(sf) + '\t |1-val| = ' + str(val[kmax]) + '\t err = ' + str(err[kmax])
    print

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'usage: ' + sys.argv[0] + ' path runType rebinFactor [command [SFtoPlot]]'
        print '       path        = path to histograms'
        print '       runType     = runZprime runType list (NO qcd, unf)'
        print '       rebinFactor = histogram rebining factor'
        print '       command     = print, plot'
        print '       SFtoPlot    = scale factor to plot'
        sys.exit()
    if len(sys.argv) == 6:
        FindSF(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], float(sys.argv[5]))
    elif len(sys.argv) == 5:
        FindSF(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
    else:
        FindSF(sys.argv[1], sys.argv[2], int(sys.argv[3]))
