import math, numpy, ROOT, sys, uuid
from PlotStyle import *
from SharedData import *

# calculate the AFB --------------------------------------------------------------------------------
def Afb(hist2D):
    histB = hist2D.ProjectionX(str(uuid.uuid4()), 1, 1, 'e')
    histF = hist2D.ProjectionX(str(uuid.uuid4()), 2, 2, 'e')
    histDIF = histF.Clone(str(uuid.uuid4()))
    histDIF.Add(histB, -1.0)
    histSUM = histF.Clone(str(uuid.uuid4()))
    histSUM.Add(histB, 1.0)
    histAFB = histDIF.Clone(str(uuid.uuid4()))
    histAFB.Divide(histDIF, histSUM, 1.0, 1.0)
    return histAFB

# extract the AFB distribution ---------------------------------------------------------------------
def ExtractAfb(fileName, histName, runType, N=0):
    dirs = []
    for i in range(N):
        if 'pdf' in runType:
            dirs.append(runType + '_' + str(i+1) + '/')
        else:
            dirs.append(runType + '_max_' + str(i+1) + '/')
            dirs.append(runType + '_min_' + str(i+1) + '/')
    if N == 0:
        dirs.append('')
    if 'ele' in runType:
        tag = 'ele'
    elif 'hf' in runType:
        tag = 'hf'
    file = ROOT.TFile.Open(fileName, 'update')
    for dir in dirs:
        hist = file.Get(dir + histName)
        if 'unf' not in histName:
            hist.RebinY(100)
        if 'DoubleElectron' not in file.GetName() and 'SingleElectron' not in file.GetName():
            hist.Scale(McSF(tag)*Lumi(tag))
        afb = Afb(hist)
        file.cd(dir)
        afb.Write(histName.replace('ll_massX_cosThetaY', 'afb'), 6)
    file.Close()

# plot the ratio -----------------------------------------------------------------------------------
def PlotRatio(ratio):
    canvas = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 130+GetH(), GetW(), GetH('R'))
    SetCanvas(canvas, 'R')
    canvas.SetLogx()
    ratio.SetStats(0)
    SetAxes(ratio, 'R')
    ratio.GetXaxis().SetTitle('M(ee)')
    xMin = ratio.GetBinLowEdge(2)
    xMax = ratio.GetBinLowEdge(ratio.GetNbinsX())
    ratio.GetXaxis().SetRangeUser(xMin, xMax)
    ratio.GetYaxis().SetTitle('(Data - MC) / #sigma')
    ratio.GetYaxis().SetRangeUser(-2.9, 2.9)
    ratio.GetYaxis().SetNdivisions(207)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(2)
    ratio.Draw('axis')
    x = [xMin, xMax, xMax, xMin]
    y = [2.0, 2.0, -2.0, -2.0]
    graph2S = ROOT.TGraph(len(x), numpy.array(x), numpy.array(y))
    graph2S.SetFillColor(ROOT.kYellow)
    graph2S.Draw('same f')
    y = [1.0, 1.0, -1.0, -1.0]
    graph1S = ROOT.TGraph(len(x), numpy.array(x), numpy.array(y))
    graph1S.SetFillColor(ROOT.kGreen)
    graph1S.Draw('same f')
    line = ROOT.TLine(xMin, 0.0, xMax, 0.0)
    line.SetLineColor(ROOT.kGray+1)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw('same')
    ratio.Draw('same hist l')
    ratio.Draw('same axis')
    canvas.Update()
    raw_input('...')

# plot the AFB distribution ------------------------------------------------------------------------
def PlotAfb(fileNameDT, histNameDT, fileNameMC, histNameMC):
    fileDT = ROOT.TFile.Open(fileNameDT)
    histDT = fileDT.Get(histNameDT)
    histDT.SetName(str(uuid.uuid4()))
    fileMC = ROOT.TFile.Open(fileNameMC)
    histMC = fileMC.Get(histNameMC)
    histMC.SetName(str(uuid.uuid4()))
    canvas = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 100, GetW(), GetH())
    SetCanvas(canvas)
    canvas.SetLogx()
    histMC.SetStats(0)
    SetAxes(histMC)
    histMC.GetXaxis().SetTitle('M(ee)')
    xMin = histMC.GetBinLowEdge(2)
    xMax = histMC.GetBinLowEdge(histMC.GetNbinsX())
    histMC.GetXaxis().SetRangeUser(xMin, xMax)
    histMC.GetYaxis().SetTitle('A_{FB}')
    histMC.GetYaxis().SetRangeUser(-0.39, 0.89)
    histMC.SetLineColor(ROOT.kRed)
    histMC.SetFillColor(ROOT.kRed)
    histMC.SetFillStyle(1001)
    histMC.Draw('e5')
    histDT.SetLineColor(ROOT.kBlack)
    histDT.SetMarkerColor(ROOT.kBlack)
    histDT.SetMarkerStyle(20)
    histDT.SetMarkerSize(1.0)
    histDT.Draw('same e p x0')
    legend = ROOT.TLegend(0.76, 0.26, 0.88, 0.38)
    SetLegend(legend)
    legend.AddEntry(histDT, ' Data', 'p')
    legend.AddEntry(histMC, ' Z#rightarrowee', 'f')
    legend.Draw()
    canvas.Update()
    ratio = histDT.Clone(str(uuid.uuid4()))
    ratio.Reset()
    for b in range(2, ratio.GetNbinsX()):
        err = math.sqrt(math.pow(histDT.GetBinError(b),2) + math.pow(histMC.GetBinError(b),2))
        dif = histDT.GetBinContent(b) - histMC.GetBinContent(b)
        try:
            ratio.SetBinContent(b, dif/err)
        except ZeroDivisionError:
            print '--> ZeroDivisionError (' + str(b) + ')'
            ratio.SetBinContent(b, 0.0)
    PlotRatio(ratio)

# plot the AFB distribution with systematics -------------------------------------------------------
def PlotAfbSys(path, runType, sysPath, sysRunType, ybin, command):
    if 'ele' in runType:
        data = 'DoubleElectron'
    elif 'hf' in runType:
        data = 'SingleElectron'
    fileDT = ROOT.TFile.Open(path + 'histograms_' + runType + '_' + data + '.root')
    histDT = fileDT.Get('afb_' + ybin + '_noBkg_unf')
    fileMC = ROOT.TFile.Open(path + 'histograms_' + runType + '_DYtoEE.root')
    histMC = fileMC.Get('gen_afb_' + ybin)
    fileS = ROOT.TFile.Open(sysPath + 'systematics_' + sysRunType + '.root')
    histP = fileS.Get('sys_max_' + ybin)
    histM = fileS.Get('sys_min_' + ybin)
    if command == 'print':
        info = '${:n}$ & ${:n}$ & ${:1.4f} \pm {:1.4f} ^{}{:+1.4f}{} _{}{:+1.4f}{}$ & ' + \
               '${:1.4f} \pm {:1.4f}$ \\\ [0.5ex]'
        for b in range(2, histDT.GetNbinsX()):
            mMin = histDT.GetBinLowEdge(b)
            mMax = histDT.GetBinLowEdge(b+1)
            afbDT = histDT.GetBinContent(b)
            errDT = histDT.GetBinError(b)
            sysP = histP.GetBinContent(b)
            sysM = histM.GetBinContent(b)
            afbMC = histMC.GetBinContent(b)
            errMC = histMC.GetBinError(b)
            print info.format(mMin,mMax,afbDT,errDT,'{',sysP,'}','{',sysM,'}',afbMC,errMC)
    x = []
    y = []
    exl = []
    exh = []
    eyl = []
    eyh = []
    for b in range(2, histDT.GetNbinsX()):
        x.append(histDT.GetBinLowEdge(b) + 0.5 * histDT.GetBinWidth(b))
        y.append(histDT.GetBinContent(b))
        exl.append(0.0)
        exh.append(0.0)
        err = histDT.GetBinError(b)
        sysP = histP.GetBinContent(b)
        sysM = histM.GetBinContent(b)
        eyl.append(math.sqrt(math.pow(err,2)+math.pow(sysM,2)))
        eyh.append(math.sqrt(math.pow(err,2)+math.pow(sysP,2)))
    graph = ROOT.TGraphAsymmErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(exl),
                                   numpy.array(exh), numpy.array(eyl), numpy.array(eyh))
    canvas = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 100, GetW(), GetH())
    SetCanvas(canvas)
    canvas.SetLogx()
    histMC.SetStats(0)
    SetAxes(histMC)
    histMC.GetXaxis().SetTitle('M(ee)')
    xMin = histMC.GetBinLowEdge(2)
    xMax = histMC.GetBinLowEdge(histMC.GetNbinsX())
    histMC.GetXaxis().SetRangeUser(xMin, xMax)
    histMC.GetYaxis().SetTitle('A_{FB}')
    histMC.GetYaxis().SetRangeUser(-0.39, 0.89)
    histMC.SetLineColor(ROOT.kRed)
    histMC.SetFillColor(ROOT.kRed)
    histMC.SetFillStyle(1001)
    histMC.Draw('e5')
    graph.SetLineColor(ROOT.kBlack)
    graph.SetMarkerColor(ROOT.kBlack)
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.0)
    graph.Draw('same p z')
    legend = ROOT.TLegend(0.76, 0.26, 0.88, 0.38)
    SetLegend(legend)
    legend.AddEntry(graph, ' Data', 'p')
    legend.AddEntry(histMC, ' Z#rightarrowee', 'f')
    legend.Draw()
    canvas.Update()
    ratio = histDT.Clone(str(uuid.uuid4()))
    ratio.Reset()
    for b in range(2, ratio.GetNbinsX()):
        err = math.sqrt(math.pow(histDT.GetBinError(b),2) + math.pow(histMC.GetBinError(b),2))
        dif = histDT.GetBinContent(b) - histMC.GetBinContent(b)
        if dif > 0.0:
            err = math.sqrt(math.pow(err,2) + math.pow(histM.GetBinContent(b),2))
        else:
            err = math.sqrt(math.pow(err,2) + math.pow(histP.GetBinContent(b),2))
        try:
            ratio.SetBinContent(b, dif/err)
        except ZeroDivisionError:
            print '--> ZeroDivisionError (' + str(b) + ')'
            ratio.SetBinContent(b, 0.0)
    PlotRatio(ratio)

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) not in [4,5,6,7]:
        print
        print 'usage:'
        print '    get AFB     : ' + sys.argv[0] + ' file histogram runType [iterations]'
        print '                  file       = name of histogram file including path'
        print '                  histogram  = name of histogram to process'
        print '                  runType    = runZprime runType list'
        print '                  iterations = number of systematic variations'
        print
        print '    plot AFB    : ' + sys.argv[0] + ' dataFile dataHistogram mcFile mcHistogram'
        print '                  dataFile      = name of DATA file including path'
        print '                  dataHistogram = name of DATA histogram to plot'
        print '                  mcFile        = name of MC file including path'
        print '                  mcHistogram   = name of MC histogram to plot'
        print
        print '    plot AFB+SYS: ' + sys.argv[0] + ' path runType sysPath sysRunType ybin [command]'
        print '                  path       = path to data and signal MC files'
        print '                  runType    = runZprime runType list'
        print '                  sysPath    = path to systematics'
        print '                  sysRunType = runZprime runType list'
        print '                  ybin       = 0y1, 1y1p25, 1p25y1p5, 1p5y2p4, 2p4y5'
        print '                  command    = print'
        print
        sys.exit()
    if len(sys.argv) is 7:
        PlotAfbSys(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    elif len(sys.argv) is 6:
        PlotAfbSys(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], '')
    elif len(sys.argv) is 5:
        try:
            val = int(sys.argv[4])
            plot = False
        except ValueError:
            plot = True
        if plot:
            PlotAfb(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
        else:
            ExtractAfb(sys.argv[1], sys.argv[2], sys.argv[3], val)
    else:
        ExtractAfb(sys.argv[1], sys.argv[2], sys.argv[3])
