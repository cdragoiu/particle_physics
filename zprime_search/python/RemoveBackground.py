import ROOT, sys, uuid
from SharedData import *

# remove backgrounds from data ---------------------------------------------------------------------
def RemoveBkgs(path, histName, runType, N=0):
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
        data = 'DoubleElectron'
    elif 'hf' in runType:
        tag = 'hf'
        data = 'SingleElectron'
    fileDT = ROOT.TFile.Open(path + 'histograms_' + runType + '_' + data + '.root', 'update')
    filesMC = []
    for mc in ['DYtoTauTau', 'QCD', 'TT', 'TW', 'VV', 'WJets']:
        filesMC.append(ROOT.TFile.Open(path + 'histograms_' + runType + '_' + mc + '.root'))
    for dir in dirs:
        histDT = fileDT.Get(dir + histName)
        histDT.SetName(str(uuid.uuid4()))
        for file in filesMC:
            hist = file.Get(dir + histName)
            if 'QCD' in file.GetName():
                hist.Scale(QcdSF(tag))
            else:
                hist.Scale(McSF(tag)*Lumi(tag))
            histDT.Add(hist, -1.0)
        fileDT.cd(dir)
        histDT.Write(histName + '_noBkg', 6)
    for file in filesMC: file.Close()
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
        RemoveBkgs(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        RemoveBkgs(sys.argv[1], sys.argv[2], sys.argv[3])
