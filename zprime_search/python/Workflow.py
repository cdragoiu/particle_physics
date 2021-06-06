from AFB import *
from QcdShape import *
from RemoveBackground import *
from Unfold import *

# workflow -----------------------------------------------------------------------------------------
def WorkFlow(path, runType, N=0):
    if 'ele' in runType:
        data = 'DoubleElectron'
    elif 'hf' in runType:
        data = 'SingleElectron'
    for ybin in ['0y1', '1y1p25', '1p25y1p5', '1p5y2p4', '2p4y5']:
        histName = 'll_massX_cosThetaY_' + ybin
        QcdShape(path, histName, runType, N)
        RemoveBkgs(path, histName, runType, N)
        respName = 'response_' + ybin
        for name in ['DYtoEE', data]:
            if 'sys' in runType and name == 'DYtoEE': continue
            for tag in ['', '_noBkg']:
                if name == 'DYtoEE' and tag == '_noBkg': continue
                Unfold(path, name, respName, histName+tag, runType, N)
            fileName = path + 'histograms_' + runType + '_' + name + '.root'
            for tag in ['', '_unf', '_noBkg', '_noBkg_unf', 'gen_']:
                if name == 'DYtoEE' and tag in ['_noBkg', '_noBkg_unf']: continue
                if name == data and tag == 'gen_': continue
                if tag == 'gen_':
                    ExtractAfb(fileName, tag+histName, runType, N)
                else:
                    ExtractAfb(fileName, histName+tag, runType, N)

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'usage: ' + sys.argv[0] + ' path runType [iterations]'
        print '       path       = path to histograms and responses'
        print '       runType    = runZprime runType list (NO qcd, unf)'
        print '       iterations = number of systematic variations'
        sys.exit()
    if len(sys.argv) is 4:
        WorkFlow(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    else:
        WorkFlow(sys.argv[1], sys.argv[2])
