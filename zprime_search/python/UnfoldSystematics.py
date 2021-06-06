import ROOT, sys

# estimate unfold systematics ----------------------------------------------------------------------
def GetSys(mcFileName, outFileName, ybin):
    file = ROOT.TFile.Open(mcFileName)
    histU = file.Get('afb_' + ybin + '_unf')
    histG = file.Get('gen_afb_' + ybin)
    histP = histG.Clone('histP')
    histP.Reset()
    histM = histG.Clone('histM')
    histM.Reset()
    for b in range(2, histG.GetNbinsX()):
        sys = abs(histG.GetBinContent(b) - histU.GetBinContent(b))
        histP.SetBinContent(b, sys)
        histM.SetBinContent(b, -sys)
    fileOUT = ROOT.TFile.Open(outFileName, 'update')
    fileOUT.cd()
    histP.Write('sys_max_' + ybin, 6)
    histM.Write('sys_min_' + ybin, 6)
    fileOUT.Close()
    file.Close()

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'usage: ' + sys.argv[0] + ' mcFile outFile ybin'
        print '       mcFile  = signal MC file'
        print '       outFile = output file for systematics'
        print '       ybin    = 0y1, 1y1p25, 1p25y1p5, 1p5y2p4, 2p4y5'
        sys.exit()
    GetSys(sys.argv[1], sys.argv[2], sys.argv[3])
