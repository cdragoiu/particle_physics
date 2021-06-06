import os, string, sys
from SubmitToCondor import *

# run jobs -----------------------------------------------------------------------------------------
def RunJobs(txtFileName, classType, runType, N=0):
    tag = txtFileName[txtFileName.find('files_')+6:txtFileName.find('.txt')]
    txtFile = open(txtFileName, 'r')
    counter = 0
    for name in txtFile:
        counter = counter + 1
        newTxtFileName = 'files_' + tag + '_' + str(counter) + '.txt'
        newTxtFile = open(newTxtFileName, 'w')
        newTxtFile.write(name.rstrip('\n'))
        newTxtFile.close()
        newTag = runType + '_' + tag + '_' + str(counter)
        outFileName = 'histograms_' + newTag + '.root'
        arguments = newTxtFileName+' '+outFileName+' '+classType+' '+runType+' '+str(N)
        scriptFileName = 'script_' + newTag + '.sh'
        scriptFile = open(scriptFileName, 'w')
        scriptFile.write('LD_LIBRARY_PATH=$LD_LIBRARY_PATH":$PWD"' + '\n')
        scriptFile.write('./runZprime ' + arguments + '\n')
        scriptFile.close()
        cmssw_base = os.environ['CMSSW_BASE']
        scram_arch = os.environ['SCRAM_ARCH']
        transferFiles = cmssw_base + '/lib/' + scram_arch + '/libRooUnfold.so,' + \
                        cmssw_base + '/lib/' + scram_arch + '/libMyAnalysesCommonTools.so,' + \
                        cmssw_base + '/bin/' + scram_arch + '/runZprime,' + \
                        newTxtFileName
        SubmitToCondor(scriptFileName, transferFiles, '', newTag)
    txtFile.close()

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'usage: ' + sys.argv[0] + ' txtFile classType runType [iterations]'
        print '       txtFile    = file with list of input files'
        print '       classType  = data, mc'
        print '       runType    = runZprime runType list'
        print '       iterations = number of systematic variations'
        sys.exit()
    if len(sys.argv) == 5:
        RunJobs(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
    else:
        RunJobs(sys.argv[1], sys.argv[2], sys.argv[3])
