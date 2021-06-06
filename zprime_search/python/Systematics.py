import math, ROOT, sys
from PlotStyle import *

# estimate systematics -----------------------------------------------------------------------------
def GetSys(basePath, baseRunType, sysPath, sysRunType, N, ybin):
    if 'ele' in baseRunType:
        data = 'DoubleElectron'
    elif 'hf' in baseRunType:
        data = 'SingleElectron'
    dirs = []
    for i in range(N):
        if 'pdf' in sysRunType:
            dirs.append(sysRunType + '_' + str(i+1) + '/')
        else:
            dirs.append(sysRunType + '_max_' + str(i+1) + '/')
            dirs.append(sysRunType + '_min_' + str(i+1) + '/')
    if N == 0:
        dirs.append('')
    histName = 'afb_' + ybin + '_noBkg_unf'
    fileB = ROOT.TFile.Open(basePath + 'histograms_' + baseRunType + '_' + data + '.root')
    histB = fileB.Get(histName)
    histP = histB.Clone('histP')
    histP.Reset()
    histM = histB.Clone('histM')
    histM.Reset()
    file = ROOT.TFile.Open(sysPath + 'histograms_' + sysRunType + '_' + data + '.root')
    if 'sysNnpdf' in sysRunType:
        hists = []
        for b in range(2, histB.GetNbinsX()):
            hists.append(ROOT.TH1D('sys_'+str(b), '', 1000, -0.5, 0.5))
        for dir in dirs:
            hist = file.Get(dir + histName)
            hist.Add(histB, -1.0)
            for b in range(2, histB.GetNbinsX()):
                hists[b-2].Fill(hist.GetBinContent(b))
        for b in range(2, histB.GetNbinsX()):
            sys = hists[b-2].GetRMS()
            histP.SetBinContent(b, sys)
            histM.SetBinContent(b, -sys)
    else:
        for dir in dirs:
            hist = file.Get(dir + histName)
            hist.Add(histB, -1.0)
            for b in range(2, histB.GetNbinsX()):
                sys = hist.GetBinContent(b)
                if sys > 0.0:
                    sysP = histP.GetBinContent(b)
                    histP.SetBinContent(b, math.sqrt(math.pow(sysP,2)+math.pow(sys,2)))
                else:
                    sysM = histM.GetBinContent(b)
                    histM.SetBinContent(b, -math.sqrt(math.pow(sysM,2)+math.pow(sys,2)))
    file.Close()
    if 'sysCt10pdf' in sysRunType:
        for b in range(2, histB.GetNbinsX()):
            sys = (histP.GetBinContent(b) - histM.GetBinContent(b)) / 2
            histP.SetBinContent(b, sys)
            histM.SetBinContent(b, -sys)
    for b in range(2, histB.GetNbinsX()):
        sysP = histP.GetBinContent(b)
        sysM = histM.GetBinContent(b)
        if -sysM < 0.5*sysP:
            histM.SetBinContent(b, -sysP)
        elif sysP < -0.5*sysM:
            histP.SetBinContent(b, -sysM)
    fileS = ROOT.TFile.Open(sysPath + 'systematics_' + sysRunType + '.root', 'update')
    fileS.cd()
    histP.Write('sys_max_' + ybin, 6)
    histM.Write('sys_min_' + ybin, 6)
    fileS.Close()
    fileB.Close()

# plot systematics ---------------------------------------------------------------------------------
def PlotSys(basePath, baseRunType, sysPath, sysRunType, ybin):
    if 'ele' in baseRunType:
        data = 'DoubleElectron'
        xmax = 2000.0
        lx1 = 0.4
    elif 'hf' in baseRunType:
        data = 'SingleElectron'
        xmax = 320.0
        lx1 = 0.1
    fileB = ROOT.TFile.Open(basePath + 'histograms_' + baseRunType + '_' + data + '.root')
    histB = fileB.Get('afb_' + ybin + '_noBkg_unf')
    histErrP = histB.Clone('histErrP')
    histErrM = histB.Clone('histErrM')
    for b in range(2, histB.GetNbinsX()):
        err = histB.GetBinError(b)
        histErrP.SetBinContent(b, err)
        histErrM.SetBinContent(b, -err)
    fileS = ROOT.TFile.Open(sysPath + 'systematics_' + sysRunType + '.root')
    histP = fileS.Get('sys_max_' + ybin)
    histM = fileS.Get('sys_min_' + ybin)
    canvas = ROOT.TCanvas('canvas', '', 440, 130, GetW(), GetH('S'))
    SetCanvas(canvas, 'S')
    canvas.SetLogx()
    ROOT.gStyle.SetHistMinimumZero()
    for hist in [histErrP, histErrM, histP, histM]:
        hist.SetBarWidth(1.0)
        hist.SetBarOffset(0.0)
        hist.GetXaxis().SetRangeUser(40.0, xmax)
    for hist in [histErrP, histErrM]:
        hist.SetFillColor(ROOT.kGray)
        hist.SetFillStyle(1001)
    for hist in [histP, histM]:
        hist.SetFillColor(ROOT.kRed)
        hist.SetFillStyle(3001)
    histErrP.SetStats(0)
    histErrP.Draw('hist bar')
    for hist in [histErrM, histP, histM]:
        hist.Draw('same hist bar')
    histErrP.Draw('same axis')
    SetAxes(histErrP, 'S')
    histErrP.GetXaxis().SetTitle('M(ee)')
    ymax = 1.5 * max(histErrP.GetMaximum(), -histErrM.GetMinimum())
    histErrP.GetYaxis().SetRangeUser(-ymax, ymax)
    legendErr = ROOT.TLegend(lx1, 0.02, lx1+0.08, 0.08)
    SetLegend(legendErr, 'S')
    legendErr.AddEntry(histErrP, ' STAT', 'f')
    legendErr.Draw('same')
    legend = ROOT.TLegend(lx1+0.12, 0.02, lx1+0.20, 0.08)
    SetLegend(legend, 'S')
    legend.AddEntry(histP, ' SYS', 'f')
    legend.Draw('same')
    canvas.Update()
    raw_input('...')

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) not in [6,7]:
        cmd = sys.argv[0]
        print 'usage:'
        print '    get SYS : ' + cmd + ' basePath baseRunType sysPath sysRunType iterations ybin'
        print '    plot SYS: ' + cmd + ' basePath baseRunType sysPath sysRunType ybin'
        print '              basePath    = path to base histograms'
        print '              baseRunType = runZprime runType list (NO qcd, unf)'
        print '              sysPath     = path to systematic histograms'
        print '              sysRunType  = runZprime runType list (NO qcd, unf)'
        print '              iterations  = number of systematic variations'
        print '              ybin        = 0y1, 1y1p25, 1p25y1p5, 1p5y2p4, 2p4y5'
        sys.exit()
    if len(sys.argv) == 7:
        GetSys(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6])
    else:
        PlotSys(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
