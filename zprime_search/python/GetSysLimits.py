import ROOT, sys

# get systematic limits ----------------------------------------------------------------------------
def GetLimits(fileName, ybin):
    file = ROOT.TFile.Open(fileName)
    histP = file.Get('sys_max_' + ybin)
    histM = file.Get('sys_min_' + ybin)
    sysMin = 9.9
    sysMax = 0.0
    for b in range(2, histP.GetNbinsX()):
        sysP = abs(histP.GetBinContent(b))
        sysM = abs(histM.GetBinContent(b))
        sys = min(sysP, sysM)
        if sys < sysMin:
            sysMin = sys
        sys = max(sysP, sysM)
        if sys > sysMax:
            sysMax = sys
    print
    print 'min: ' + str(sysMin)
    print 'max: ' + str(sysMax)
    print

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print
        print 'usage: ' + sys.argv[0] + ' file ybin'
        print '       file = input file'
        print '       ybin = 0y1, 1y1p25, 1p25y1p5, 1p5y2p4, 2p4y5'
        print
        sys.exit()
    GetLimits(sys.argv[1], sys.argv[2])
