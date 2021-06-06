from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'DYtoEE_13TeV'
config.General.workArea = '/uscms_data/d2/cdragoiu/'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'TreeProducer.py'
config.JobType.outputFiles = ['ntuple.root']

config.section_('Data')
ds = '/DYToEE_M-50_Tune4C_13TeV-pythia8/Phys14DR-PU20bx25_tsg_castor_PHYS14_25_V1-v1/AODSIM'
config.Data.inputDataset = ds
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.outLFN = '/store/user/dragoiu/DYtoEE_13TeV'

config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'
