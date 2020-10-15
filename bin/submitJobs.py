import sys
import os
import subprocess
import readline
import string
import csv, subprocess
from GFAL_GetROOTfiles import *


import argparse
# set up an argument parser                                                                                                                                                                                        
parser = argparse.ArgumentParser()

parser.add_argument('--v', dest='VERBOSE', default=True)
parser.add_argument('--l', dest = 'LOCATION', default= '/afs/cern.ch/user/a/asparker/public/LFVTopCode_MyFork/old_for_comparison/TopLFV/')
parser.add_argument('--n', dest = 'NAMETAG', default= 'SMEFTfr' )

ARGS = parser.parse_args()

verbose = ARGS.VERBOSE
loc = ARGS.LOCATION
name = ARGS.NAMETAG

dire = loc + 'hists/'

import Files_2016
import Files_2017
import Files_2018
import Files_2017_A

SAMPLES = {}
mc_2016 = False
data_2016 = False
mc_2017 = False
data_2017 = False
mc_2018 = False
data_2018 = False

SAMPLES.update(Files_2017.mc2017_samples)
SAMPLES.update(Files_2017.data2017_samples)

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

submit = 'universe = vanilla\n' ##writing .sub file
submit += 'arguments = "$(argument)"\n'
submit += 'output = submit01.out\n'
submit += 'error = submit01.err\n'
submit += 'log = submit01.log\n'
submit += '+JobFlavour = "tomorrow"\n' ##finish writing .sh file
submit += 'queue\n'
submitName = 'submit01.sub'
sub1 = open('Jobs/'+submitName,'wt')
sub1.write(submit+'\n')
sub1.close()

for key, value in SAMPLES.items():
    #if name  not in key:
    #   continue
    year = value[3]
    nf = 40
    for idx, S in enumerate(value[0]):
        if value[1]=='data':
            nf = 255
        #for subdir, dirs, files in os.walk(S):
        filelist = GFAL_GetROOTfiles( S ,"srm://maite.iihe.ac.be:8443/srm/managerv2?SFN=" )#"srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2?SFN=/storage/data/cms")           
        #for files in filelist:
        sequance = [filelist[i:i+nf] for i in range(0,len(filelist),nf)]
        for num,  seq in enumerate(sequance):
                f = key +'_' + str(idx) +'_' + str(num)
                subprocess.call('rm '+ dire + year + '/' + f + '.root', shell=True)
                qsub = "condor_submit Jobs/"+ submitName +" executable=Jobs/"+ f + '.sh'
                subprocess.call(qsub, shell=True)
         #   break

