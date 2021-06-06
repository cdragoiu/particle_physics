import ROOT, sys

# save MC pileup information -----------------------------------------------------------------------
def GetPileup(fileNames):
    ROOT.gSystem.Load('libFWCoreFWLite')
    ROOT.gSystem.Load('libDataFormatsFWLite')
    ROOT.AutoLibraryLoader.enable()
    chain = ROOT.TChain('Events')
    for name in fileNames:
        chain.Add(name)
    hist = ROOT.TH1D('pileup', '', 60, 0.0, 60.0)
    print '--> using ' + str(chain.GetEntries()) + ' events'
    chain.Draw('PileupSummaryInfos_addPileupInfo__HLT.obj.getTrueNumInteractions() >> +pileup',
               'PileupSummaryInfos_addPileupInfo__HLT.obj.getBunchCrossing() == 0')
    hist.SaveAs('pileup_MC.root')
    raw_input('...')

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'usage: ' + sys.argv[0] + ' txtFile'
        print '       txtFile = file with list of input files'
        sys.exit()
    fileNames = []
    txtFile = open(sys.argv[1], 'r')
    for fileName in txtFile:
        fileNames.append(fileName.rstrip('\n'))
    GetPileup(fileNames)
