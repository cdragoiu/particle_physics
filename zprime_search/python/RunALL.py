import string, subprocess, sys
from RunOnCondor import *

# run parameters -----------------------------------------------------------------------------------
class Run:
    def __init__(self, _name, _size, _tag, _type):
        self.name = _name
        self.size = _size
        self.tag  = _tag
        self.type = _type

# execute jobs -------------------------------------------------------------------------------------
def Execute(command, runClass, runType, N=0):
    if 'ele' in runType:
        data = Run('DoubleElectron', 18, 'data', '')
    elif 'hf' in runType:
        data = Run('SingleElectron', 52, 'data', '')
    runs = []
    if runClass == 'base':
        for qcdTag in ['', 'qcd_']:
            runs.append(Run('DYtoEE', 19, 'mc', qcdTag + runType))
            if qcdTag == '':
                runs.append(Run('DYtoEE', 19, 'mc', qcdTag + 'unf_' + runType))
            runs.append(Run('DYtoTauTau', 12, 'mc', qcdTag + runType))
            runs.append(Run(data.name, data.size, data.tag, qcdTag + runType))
            runs.append(Run('TT', 7, 'mc', qcdTag + runType))
            runs.append(Run('TW', 2, 'mc', qcdTag + runType))
            runs.append(Run('VV', 9, 'mc', qcdTag + runType))
            runs.append(Run('WJets', 35, 'mc', qcdTag + runType))
    elif runClass == 'sysDATA':
        for qcdTag in ['', 'qcd_']:
            runs.append(Run(data.name, data.size, data.tag, qcdTag + runType))
    elif runClass == 'sysMC':
        for qcdTag in ['', 'qcd_']:
            if qcdTag == '':
                runs.append(Run('DYtoEE', 19, 'mc', qcdTag + 'unf_' + runType))
            else:
                runs.append(Run('DYtoEE', 19, 'mc', qcdTag + runType))
            runs.append(Run('DYtoTauTau', 12, 'mc', qcdTag + runType))
            runs.append(Run('TT', 7, 'mc', qcdTag + runType))
            runs.append(Run('TW', 2, 'mc', qcdTag + runType))
            runs.append(Run('VV', 9, 'mc', qcdTag + runType))
            runs.append(Run('WJets', 35, 'mc', qcdTag + runType))
    for run in runs:
        if command == 'submit':
            RunJobs('data/files_' + run.name + '.txt', run.tag, run.type, N)
        elif command == 'check':
            for i in range(run.size):
                fileName = 'histograms_' + run.type + '_' + run.name + '_' + str(i+1) + '.root'
                try:
                    file = open(fileName, 'r')
                    file.close()
                except IOError:
                    print 'missing: ' + fileName
        elif command == 'hadd':
            fileTag = 'histograms_' + run.type + '_' + run.name
            if 'unf' in run.type:
                callArgs = 'python python/CombineResponses.py ' + fileTag + ' responses_' + \
                           run.type + '.root ' + run.type + ' ' + str(N)
            else:
                callArgs = 'hadd ' + fileTag + '.root ' + fileTag + '_*'
            subprocess.call(callArgs, shell = True)
            subprocess.call('rm ' + fileTag + '_*', shell = True)

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'usage: ' + sys.argv[0] + ' command runClass runType [iterations]'
        print '       command    = submit, check, hadd'
        print '       runClass   = base, sysDATA, sysMC'
        print '       runType    = runZprime runType list (NO qcd, unf)'
        print '       iterations = number of systematic variations'
        sys.exit()
    if len(sys.argv) == 5:
        Execute(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        Execute(sys.argv[1], sys.argv[2], sys.argv[3])
