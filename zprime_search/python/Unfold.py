import os, ROOT, sys, uuid
from PlotStyle import *

# unfold histogram ---------------------------------------------------------------------------------
def Unfold(path, fileTag, respName, histName, runType, N=0):
    ROOT.gSystem.Load(os.environ['CMSSW_BASE'] + '/src/RooUnfold/libRooUnfold.so')
    dirs = []
    for i in range(N):
        if 'pdf' in runType:
            dirs.append('unf_' + runType + '_' + str(i+1) + '/')
        else:
            dirs.append('unf_' + runType + '_max_' + str(i+1) + '/')
            dirs.append('unf_' + runType + '_min_' + str(i+1) + '/')
    if N == 0:
        dirs.append('')
    fileRESP = ROOT.TFile.Open(path + 'responses_unf_' + runType + '.root')
    fileHIST = ROOT.TFile.Open(path + 'histograms_' + runType + '_' + fileTag + '.root', 'update')
    for dir in dirs:
        resp = fileRESP.Get(dir + respName)
        hist = fileHIST.Get(dir.replace('unf_', '') + histName)
        hist.RebinY(100)
        bayes = ROOT.RooUnfoldBayes(resp, hist, 4)
        histUNF = bayes.Hreco(1)
        fileHIST.cd(dir.replace('unf_', ''))
        histUNF.SetTitle('')
        histUNF.Write(histName + '_unf', 6)
    fileHIST.Close()
    fileRESP.Close()

# plot response ------------------------------------------------------------------------------------
def Plot(fileName, respName):
    ROOT.gROOT.LoadMacro(os.environ['CMSSW_BASE'] + '/src/RooUnfold/libRooUnfold.so')
    file = ROOT.TFile.Open(fileName)
    resp = file.Get(respName)
    hist = resp.Hresponse()
    canvas = ROOT.TCanvas(str(uuid.uuid4()), '', 440, 100, GetW(), GetH())
    SetCanvas(canvas)
    canvas.SetMargin(0.12, 0.10, 0.15, 0.04)
    canvas.SetLogz()
    hist.SetStats(0)
    hist.SetTitle('')
    SetAxes(hist)
    hist.GetXaxis().SetTitle('RECO')
    hist.GetXaxis().SetTitleOffset(1.6)
    hist.GetYaxis().SetTitle('GEN')
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.SetContour(99)
    hist.Draw('colz')
    if '2p4y5' in respName:
        bins = [0, 40, 76, 86, 96, 106, 120, 150, 320, 8000]
    else:
        bins = [0, 40, 50, 60, 76, 86, 96, 106, 120, 133, 150, 171, 200, 320, 500, 2000, 8000]
    N = len(bins)-1
    for b in range(N):
        hist.GetXaxis().SetBinLabel(b+1, str(bins[b])+'-'+str(bins[b+1]))
        hist.GetXaxis().SetBinLabel(b+1+N, str(bins[b])+'-'+str(bins[b+1]))
        hist.GetYaxis().SetBinLabel(b+1, str(bins[b])+'-'+str(bins[b+1]))
        hist.GetYaxis().SetBinLabel(b+1+N, str(bins[b])+'-'+str(bins[b+1]))
    hist.GetXaxis().LabelsOption('v')
    hist.GetYaxis().LabelsOption('h')
    lineH = ROOT.TLine(0, N, 2*N, N)
    lineH.SetLineStyle(3)
    lineH.Draw('same')
    lineV = ROOT.TLine(N, 0, N, 2*N)
    lineV.SetLineStyle(3)
    lineV.Draw('same')
    canvas.Update()
    raw_input('...')

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) not in [3,6,7]:
        cmd = sys.argv[0]
        print 'usage:'
        print '    unfold data  : ' + cmd + ' path fileTag respName histName runType [iterations]'
        print '                   path       = path to histograms and responses'
        print '                   fileTag    = DYtoEE, DoubleElectron, SingleElectron'
        print '                   respName   = name of response matrix to use'
        print '                   histName   = name of histogram to unfold'
        print '                   runType    = runZprime runType list (NO qcd, unf)'
        print '                   iterations = number of systematic variations'
        print '    plot response: ' + cmd + ' file response'
        print '                   file     = name of response file including path'
        print '                   response = name of response matrix to plot'
        sys.exit()
    if len(sys.argv) is 3:
        Plot(sys.argv[1], sys.argv[2])
    elif len(sys.argv) is 7:
        Unfold(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6]))
    else:
        Unfold(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
