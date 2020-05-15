import sys
import os
import subprocess
import readline
import string
import Files_2016
import Files_2017
import Files_2018
SAMPLES = {}
#SAMPLES ['DYM10to50'] = ['address', 'data/mc','dataset','year', 'run', 'cross section','lumi','Neventsraw']
mc_2016 = True
data_2016 = True
mc_2017 = True
data_2017 = True
mc_2018 = True
data_2018 = True

if mc_2016:
    SAMPLES.update(Files_2016.mc2016_samples)
if data_2016:
    SAMPLES.update(Files_2016.data2016_samples)
if mc_2017:
    SAMPLES.update(Files_2017.mc2017_samples)
if data_2017:
    SAMPLES.update(Files_2017.data2017_samples)
if mc_2018:
    SAMPLES.update(Files_2018.mc2018_samples)
if data_2018:
    SAMPLES.update(Files_2018.data2018_samples)

#print SAMPLES

print "Now printing gfal-copy commands"
for key, value in SAMPLES.items():
    #Use gfal-copy to move files from --loc to a newly created EOS directory --name  
    #Example use :
    #python makeGFAL-COPYtxtfiles.py --name 2016_DYM10to50_ext0 --loc srm://maite.iihe.ac.be:8443/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3/191130_081246/0000/
 
    # If there is more than one sample, label them ext0, ext1, ext2 etc....
    
    
    tkey = key
    #print value
    #print len(value[0])
    
    ext = '_ext'
  
    ### value[0] is a list of file locations
    for iv, v in enumerate(value[0]) :
        print v
        extt = ext + str(iv)
        if len(value[0]) > 1 :
            tkey = key + extt 

        osc = "python makeGFAL-COPYtxtfiles.py --name  "+ tkey + " --loc  srm://maite.iihe.ac.be:8443"+v
        print osc
        os.system( osc)   
#[['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_vector_emutu/'], 'mc', '', '2017', '', '0.515', '35.92', '494000']
#1




