import math, ROOT, string, sys

# combine systematics ------------------------------------------------------------------------------
def Combine(path, outFileName, ybin):
    if 'ele' in path:
        if 'sysPdf' in path:
            sysTags = ['mc_ele_sysCt10pdf', 'mc_ele_sysNnpdf']
        else:
            sysTags = ['mc_ele_sysEleErg', 'mc_ele_sysEleId', 'mc_ele_sysFsr',
                       'mc_ele_sysLumi', 'mc_ele_sysPdf', 'mc_ele_sysPu', 'mc_ele_sysQcd',
                       'mc_ele_sysTrg8', 'mc_ele_sysTrg17', 'mc_ele_sysUnf']
    elif 'hf' in path:
        if 'sysPdf' in path:
            sysTags = ['mc_hf_sysCt10pdf', 'mc_hf_sysNnpdf']
        else:
            sysTags = ['mc_hf_sysEleErg', 'data_hf_sysHfErg', 'mc_hf_sysHfErg',
                       'mc_hf_sysEleId', 'mc_hf_sysFsr', 'mc_hf_sysHfEta', 'mc_hf_sysHfId',
                       'mc_hf_sysLumi', 'mc_hf_sysMcNorm', 'mc_hf_sysPdf', 'mc_hf_sysPu',
                       'mc_hf_sysQcd', 'mc_hf_sysTrg27', 'mc_hf_sysUnf']
    files = []
    histsP = []
    histsM = []
    for i in range(len(sysTags)):
        fileTag = string.replace(string.replace(sysTags[i], 'mc_', ''), 'data_', '')
        fileName = path + 'results_' + sysTags[i] + '/systematics_' + fileTag + '.root'
        files.append(ROOT.TFile.Open(fileName))
        histsP.append(files[-1].Get('sys_max_' + ybin))
        histsM.append(files[-1].Get('sys_min_' + ybin))
    histP = histsP[0].Clone('histP')
    histP.Reset()
    histM = histsM[0].Clone('histM')
    histM.Reset()
    for i in range(len(files)):
        for b in range(2, histP.GetNbinsX()):
            sysP = histsP[i].GetBinContent(b)
            sysM = histsM[i].GetBinContent(b)
            if 'sysPdf' in path:
                if sysP > histP.GetBinContent(b):
                    histP.SetBinContent(b, sysP)
                if sysM < histM.GetBinContent(b):
                    histM.SetBinContent(b, sysM)
            else:
                sys = histP.GetBinContent(b)
                histP.SetBinContent(b, math.sqrt(math.pow(sys,2)+math.pow(sysP,2)))
                sys = histM.GetBinContent(b)
                histM.SetBinContent(b, -math.sqrt(math.pow(sys,2)+math.pow(sysM,2)))
    fileOUT = ROOT.TFile.Open(outFileName, 'update')
    fileOUT.cd()
    histP.Write('sys_max_' + ybin, 6)
    histM.Write('sys_min_' + ybin, 6)
    fileOUT.Close()
    for i in range(len(files)):
        files[i].Close()

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'usage: ' + sys.argv[0] + ' path outFile ybin'
        print '       path    = common path to files'
        print '       outFile = output root file'
        print '       ybin    = 0y1, 1y1p25, 1p25y1p5, 1p5y2p4, 2p4y5'
        sys.exit()
    Combine(sys.argv[1], sys.argv[2], sys.argv[3])
