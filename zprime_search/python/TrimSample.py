import ROOT, sys

# trim sample --------------------------------------------------------------------------------------
def Trim(inFileName, outFileName, massCut, weight):
    inFile = ROOT.TFile.Open(inFileName)
    inTree = inFile.Get('Events')
    outFile = ROOT.TFile(outFileName, 'RECREATE')
    outTree = inTree.CloneTree(0)
    print '--> processing ' + str(inTree.GetEntries()) + ' events'
    for evt in range(inTree.GetEntries()):
        inTree.GetEntry(evt)
        if inTree.gen_mass[0] < massCut:
            outTree.Fill()
    print '--> ' + str(outTree.GetEntries()) + ' events selected'
    outTree.SetWeight(weight)
    outFile.cd()
    outTree.Write()
    outFile.Close()
    inFile.Close()

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 5:
        print 'usage: ' + sys.argv[0] + ' inFile outFile massCut weight'
        print '       inFile  = input root file'
        print '       outFile = output root file'
        print '       massCut = upper bound on generated Z mass'
        print '       weight  = new tree weight'
        sys.exit()
    Trim(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))
