import ROOT, sys
from SharedData import *

# get QCD shape ------------------------------------------------------------------------------------
def QcdShape(path, histName, runType, N=0):
    dirs = []
    for i in range(N):
        if 'pdf' in runType:
            dirs.append('qcd_' + runType + '_' + str(i+1) + '/')
        else:
            dirs.append('qcd_' + runType + '_max_' + str(i+1) + '/')
            dirs.append('qcd_' + runType + '_min_' + str(i+1) + '/')
    if N == 0:
        dirs.append('')
    if 'ele' in runType:
        tag = 'ele'
        data = 'DoubleElectron'
    elif 'hf' in runType:
        tag = 'hf'
        data = 'SingleElectron'
    fileDT = ROOT.TFile.Open(path + 'histograms_qcd_' + runType + '_' + data + '.root')
    filesMC = []
    for mc in ['DYtoEE', 'DYtoTauTau', 'TT', 'TW', 'VV', 'WJets']:
        filesMC.append(ROOT.TFile.Open(path + 'histograms_qcd_' + runType + '_' + mc + '.root'))
    fileQCD = ROOT.TFile.Open(path + 'histograms_' + runType + '_QCD.root', 'update')
    for dir in dirs:
        histQCD = fileDT.Get(dir + histName)
        for fileMC in filesMC:
            histMC = fileMC.Get(dir + histName)
            histMC.Scale(McSF(tag)*Lumi(tag))
            histQCD.Add(histMC, -1.0)
        if not fileQCD.GetDirectory(dir.replace('qcd_','')):
            fileQCD.mkdir(dir.replace('qcd_',''))
        fileQCD.cd(dir.replace('qcd_',''))
        histQCD.Write(histName, 6)
    fileQCD.Close()
    for fileMC in filesMC: fileMC.Close()
    fileDT.Close()

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'usage: ' + sys.argv[0] + ' path histName runType [iterations]'
        print '       path       = path to histograms'
        print '       histName   = name of histogram to process'
        print '       runType    = runZprime runType list (NO qcd, unf)'
        print '       iterations = number of systematic variations'
        sys.exit()
    if len(sys.argv) == 5:
        QcdShape(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        QcdShape(sys.argv[1], sys.argv[2], sys.argv[3])
