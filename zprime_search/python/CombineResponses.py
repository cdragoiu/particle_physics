import os, ROOT, sys, uuid

# combine responses --------------------------------------------------------------------------------
def Combine(inFilePrefix, outFileName, runType='', N=0):
    ROOT.gSystem.Load(os.environ['CMSSW_BASE'] + '/src/RooUnfold/libRooUnfold.so')
    ybins = ['0y1', '1y1p25', '1p25y1p5', '1p5y2p4', '2p4y5']
    dirs = []
    for i in range(N):
        if 'pdf' in runType:
            dirs.append(runType + '_' + str(i+1) + '/')
        else:
            dirs.append(runType + '_max_' + str(i+1) + '/')
            dirs.append(runType + '_min_' + str(i+1) + '/')
    if N == 0:
        dirs.append('')
    responses = []
    for i in range(len(dirs)):
        rlist = []
        for j in range(len(ybins)):
            rlist.append(ROOT.RooUnfoldResponse(str(uuid.uuid4()), ''))
        responses.append(rlist)
    for name in os.listdir('.'):
        if inFilePrefix in name:
            file = ROOT.TFile.Open(name)
            for i in range(len(dirs)):
                for j in range(len(ybins)):
                    resp = file.Get(dirs[i] + 'response_' + ybins[j])
                    responses[i][j].Add(resp)
            file.Close()
    fileOUT = ROOT.TFile.Open(outFileName, 'update')
    for i in range(len(dirs)):
        fileOUT.mkdir(dirs[i])
        fileOUT.cd(dirs[i])
        for j in range(len(ybins)):
            responses[i][j].Write('response_' + ybins[j], 6)
    fileOUT.Close()

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'usage: ' + sys.argv[0] + ' inFilePrefix outFile [runType iterations]'
        print '       inFilePrefix = common name of input files'
        print '       outFile      = output root file'
        print '       runType      = runZprime runType list'
        print '       iterations   = number of systematic variations'
        sys.exit()
    if len(sys.argv) == 5:
        Combine(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        Combine(sys.argv[1], sys.argv[2])
