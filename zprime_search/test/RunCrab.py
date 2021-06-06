import os, subprocess

# crab job characteristics -------------------------------------------------------------------------
class Job:
    def __init__(self, scheduler_, site_, lumiMask_, dataset_, outDir_, runType_, nrJobs_):
        self.scheduler = scheduler_
        self.site = site_
        self.lumiMask = lumiMask_
        self.dataset = dataset_
        self.outDir = outDir_
        self.runType = runType_
        self.nrJobs = nrJobs_

# run crab -----------------------------------------------------------------------------------------
def RunCrab(job):
    cfgFile = open('crab.cfg', 'w')
    cfgFile.write('[CRAB]' + '\n\n')
    cfgFile.write('jobtype = cmssw' + '\n')
    cfgFile.write('use_server = 0' + '\n')
    cfgFile.write('scheduler = ' + job.scheduler + '\n\n')
    cfgFile.write('[CMSSW]' + '\n\n')
    cfgFile.write('datasetpath = ' + job.dataset + '\n')
    cfgFile.write('pset = TreeProducer.py' + '\n')
    cfgFile.write('pycfg_params = TreeProducer.py ' + job.runType + '\n')
    if job.runType == 'data':
        cfgFile.write('lumi_mask = ' + job.lumiMask + '\n')
        cfgFile.write('total_number_of_lumis = -1' + '\n')
    else:
        cfgFile.write('total_number_of_events = -1' + '\n')
    cfgFile.write('number_of_jobs = ' + job.nrJobs + '\n')
    cfgFile.write('output_file = ntuple.root' + '\n\n')
    cfgFile.write('[USER]' + '\n\n')
    cfgFile.write('ui_working_dir = /uscms_data/d2/cdragoiu/' + job.outDir + '\n')
    if job.runType == 'mc':
        cfgFile.write('additional_input_files = pileup_MC.root, pileup_DATA.root' + '\n')
    cfgFile.write('copy_data = 1' + '\n')
    cfgFile.write('storage_element = cmseos.fnal.gov' + '\n')
    cfgFile.write('storage_path = /srm/v2/server?SFN=/eos/uscms/store/user/cdragoiu/' + '\n')
    cfgFile.write('user_remote_dir = ' + job.outDir + '\n\n')
    cfgFile.write('[GRID]' + '\n\n')
    cfgFile.write('ce_white_list = ' + job.site + '\n')
    cfgFile.write('se_white_list = ' + job.site + '\n')
    cfgFile.close()
    subprocess.call('crab -create', shell = True)
    os.remove('crab.cfg')

# run as the main program only ---------------------------------------------------------------------
if __name__ == '__main__':
    
    jobs = [
            Job('condor', 'fnal', '',
                '/DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'DYToEE/M20', 'mc', '2100'),
            Job('condor', 'fnal', '',
                '/DYToEE_M-200_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'DYToEE/M200', 'mc', '5'),
            Job('condor', 'fnal', '',
                '/DYToEE_M-400_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'DYToEE/M400', 'mc', '5'),
            Job('condor', 'fnal', '',
                '/DYToEE_M-500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'DYToEE/M500', 'mc', '5'),
            Job('condor', 'fnal', '',
                '/DYToEE_M-700_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'DYToEE/M700', 'mc', '5'),
            Job('condor', 'fnal', '',
                '/DYToEE_M-800_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'DYToEE/M800', 'mc', '5'),
            Job('condor', 'fnal', '',
                '/DYToEE_M-1000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'DYToEE/M1000', 'mc', '5'),
            Job('condor', 'fnal', '',
                '/DYToEE_M-1500_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'DYToEE/M1500', 'mc', '5'),
            Job('condor', 'fnal', '',
                '/DYToEE_M-2000_CT10_TuneZ2star_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'DYToEE/M2000', 'mc', '5'),
            Job('condor', 'fnal', '',
                '/DYToTauTau_M-20_CT10_TuneZ2star_v2_8TeV-powheg-tauola-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
                'DYToTauTau/M20', 'mc', '2400'),
            Job('condor', 'fnal', 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12.json',
                '/DoubleElectron/Run2012A-22Jan2013-v1/AOD',
                'DoubleElectron/Run2012A', 'data', '600'),
            Job('condor', 'fnal', 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12.json',
                '/DoubleElectron/Run2012B-22Jan2013-v1/AOD',
                'DoubleElectron/Run2012B', 'data', '1100'),
            Job('condor', 'fnal', 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12.json',
                '/DoubleElectron/Run2012C-22Jan2013-v1/AOD',
                'DoubleElectron/Run2012C', 'data', '1600'),
            Job('condor', 'fnal', 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12.json',
                '/DoubleElectron/Run2012D-22Jan2013-v1/AOD',
                'DoubleElectron/Run2012D', 'data', '1700'),
            Job('condor', 'fnal', '',
                '/G_Pt-15to30_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/15pt30', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/G_Pt-30to50_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/30pt50', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/G_Pt-50to80_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/50pt80', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/G_Pt-80to120_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/80pt120', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/G_Pt-120to170_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/120pt170', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/G_Pt-170to300_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/170pt300', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/G_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/300pt470', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/G_Pt-470to800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/470pt800', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/G_Pt-800to1400_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/800pt1400', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/G_Pt-1400to1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/1400pt1800', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/G_Pt-1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Photon/1800pt8000', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/QCD_Pt_20_30_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_BC/20pt30', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/QCD_Pt_30_80_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_BC/30pt80', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/QCD_Pt_80_170_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_BC/80pt170', 'mc', '100'),
            Job('glite', 'infn', '',
                '/QCD_Pt_170_250_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_BC/170pt250', 'mc', '100'),
            Job('glite', 'unl', '',
                '/QCD_Pt_250_350_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_BC/250pt350', 'mc', '100'),
            Job('glite', 'infn', '',
                '/QCD_Pt_350_BCtoE_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
                'QCD_BC/350pt8000', 'mc', '100'),
            Job('condor', 'fnal', '',
                '/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_EM/20pt30', 'mc', '1700'),
            Job('condor', 'fnal', '',
                '/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_EM/30pt80', 'mc', '1700'),
            Job('condor', 'fnal', '',
                '/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_EM/80pt170', 'mc', '1700'),
            Job('condor', 'fnal', '',
                '/QCD_Pt_170_250_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_EM/170pt250', 'mc', '1700'),
            Job('condor', 'fnal', '',
                '/QCD_Pt_250_350_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_EM/250pt350', 'mc', '1700'),
            Job('condor', 'fnal', '',
                '/QCD_Pt_350_EMEnriched_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'QCD_EM/350pt8000', 'mc', '1700'),
            Job('condor', 'fnal', 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12.json',
                '/SingleElectron/Run2012A-22Jan2013-v1/AOD',
                'SingleElectron/Run2012A', 'data', '1000'),
            Job('condor', 'fnal', 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12.json',
                '/SingleElectron/Run2012B-22Jan2013-v1/AOD',
                'SingleElectron/Run2012B', 'data', '3300'),
            Job('condor', 'fnal', 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12.json',
                '/SingleElectron/Run2012C-22Jan2013-v1/AOD',
                'SingleElectron/Run2012C', 'data', '5000'),
            Job('condor', 'fnal', 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12.json',
                '/SingleElectron/Run2012D-22Jan2013-v1/AOD',
                'SingleElectron/Run2012D', 'data', '5300'),
            Job('glite', 'desy', '',
                '/TT_CT10_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'TT/v1', 'mc', '300'),
            Job('condor', 'fnal', '',
                '/TT_CT10_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
                'TT/v2', 'mc', '1000'),
            Job('condor', 'fnal', '',
                '/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'T_tW', 'mc', '20'),
            Job('condor', 'fnal', '',
                '/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'Tbar_tW', 'mc', '20'),
            Job('condor', 'fnal', '',
                '/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
                'WJetsToLNu', 'mc', '2800'),
            Job('glite', 'wisc', '',
                '/WJetsToLNu_PtW-50To70_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'WJetsToLNu/50pt70', 'mc', '2400'),
            Job('glite', 'infn', '',
                '/WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'WJetsToLNu/70pt100', 'mc', '1100'),
            Job('condor', 'fnal', '',
                '/WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'WJetsToLNu/100pt8000', 'mc', '600'),
            Job('condor', 'fnal', '',
                '/WW_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'WW', 'mc', '500'),
            Job('condor', 'fnal', '',
                '/WZ_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'WZ', 'mc', '500'),
            Job('condor', 'fnal', '',
                '/ZZ_TuneZ2star_8TeV_pythia6_tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
                'ZZ', 'mc', '500')
           ]

for job in jobs:
    RunCrab(job)
