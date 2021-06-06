import os, ROOT, sys

# make files ---------------------------------------------------------------------------------------
def MakeFiles(basePath, baseRunType, sysPath, sysRunType, N, sysClass):
    ROOT.gSystem.Load(os.environ['CMSSW_BASE'] + '/src/RooUnfold/libRooUnfold.so')
    if 'ele' in baseRunType:
        data = '_DoubleElectron'
    elif 'hf' in baseRunType:
        data = '_SingleElectron'
    if sysClass == 'data':
        names = ['', '_DYtoEE', '_DYtoTauTau', '_TT', '_TW', '_VV', '_WJets']
    elif sysClass == 'mc':
        names = [data]
    dirs = []
    for i in range(N):
        if 'pdf' in sysRunType:
            dirs.append(sysRunType + '_' + str(i+1))
        else:
            dirs.append(sysRunType + '_max_' + str(i+1))
            dirs.append(sysRunType + '_min_' + str(i+1))
    for name in names:
        for tag in ['', 'qcd_', 'unf_']:
            if name in ['', '_DYtoEE'] and tag == '': continue
            if name == '' and tag == 'qcd_': continue
            if name != '' and tag == 'unf_': continue
            if name == '' and tag == 'unf_':
                fileTag = 'responses_'
                objTag = 'response_'
            else:
                fileTag = 'histograms_'
                objTag = 'll_massX_cosThetaY_'
            baseFile = ROOT.TFile.Open(basePath+fileTag+tag+baseRunType+name+'.root')
            sysFile = ROOT.TFile.Open(sysPath+fileTag+tag+sysRunType+name+'.root','recreate')
            for dir in dirs:
                sysFile.mkdir(tag + dir)
                for ybin in ['0y1', '1y1p25', '1p25y1p5', '1p5y2p4', '2p4y5']:
                    obj = baseFile.Get(objTag + ybin)
                    sysFile.cd(tag + dir)
                    obj.Write(objTag + ybin, 6)
            sysFile.Close()
            baseFile.Close()

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 7:
        print 'usage: '+sys.argv[0]+' basePath baseRunType sysPath sysRunType iterations sysClass'
        print '       basePath    = path to base histograms'
        print '       baseRunType = base runType (NO qcd, unf)'
        print '       sysPath     = path to systematic histograms'
        print '       sysRunType  = systematic runType (NO qcd, unf)'
        print '       iterations  = number of systematic variations'
        print '       sysClass    = data, mc'
        sys.exit()
    MakeFiles(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6])
