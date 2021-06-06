import os, subprocess, sys

# submit to condor ---------------------------------------------------------------------------------
def SubmitToCondor(executable, transferFiles, arguments, tag):
    configFile = open('condor.cfg', 'w')
    configFile.write('universe = vanilla' + '\n')
    configFile.write('requirements = OpSys == "LINUX" && Arch != "DUMMY"' + '\n')
    configFile.write('request_disk = 3000000' + '\n')
    configFile.write('request_memory = 300' + '\n')
    configFile.write('getenv = True' + '\n')
    configFile.write('executable = ' + executable + '\n')
    configFile.write('transfer_input_files = ' + transferFiles + '\n')
    configFile.write('arguments = ' + arguments + '\n')
    configFile.write('should_transfer_files = YES' + '\n')
    configFile.write('when_to_transfer_output = ON_EXIT' + '\n')
    configFile.write('notification = Error' + '\n')
    configFile.write('notify_user = cdragoiu@fnal.gov' + '\n')
    configFile.write('log = condor_log_' + tag + '.txt' + '\n')
    configFile.write('output = condor_out_' + tag + '.txt' + '\n')
    configFile.write('error = condor_err_' + tag + '.txt' + '\n')
    configFile.write('queue' + '\n')
    configFile.close()
    subprocess.call('condor_submit condor.cfg', shell = True)
    os.remove('condor.cfg')

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'usage: ' + sys.argv[0] + ' txtFile runType'
        print '       txtFile = file with list of input files'
        print '       runType = data, mc'
        sys.exit()
    txtFile = open(sys.argv[1], 'r')
    runType = sys.argv[2]
    transferFiles = 'TreeProducer.py'
    if 'mc' in runType:
        transferFiles = transferFiles + ', pileup_MC.root, pileup_DATA.root'
    counter = 0
    for name in txtFile:
        counter = counter + 1
        arguments = 'TreeProducer.py ' + runType + ' ' + name.rstrip('\n') + \
                    ' ntuple_' + str(counter) + '.root'
        SubmitToCondor('cmsRun', transferFiles, arguments, str(counter))
    txtFile.close()
