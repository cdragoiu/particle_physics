import ROOT, sys

# set tree weight ----------------------------------------------------------------------------------
def SetWeight(fileName, treeName, weight):
    file = ROOT.TFile.Open(fileName, 'UPDATE')
    tree = file.Get(treeName)
    tree.SetWeight(weight)
    tree.Write()
    file.Close()

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'usage: ' + sys.argv[0] + ' fileName treeName weight'
        print '       fileName = name of input file'
        print '       treeName = name of root tree'
        print '       weight   = new tree weight'
        sys.exit()
    SetWeight(sys.argv[1], sys.argv[2], float(sys.argv[3]))
