import sys
import os
import subprocess
import readline
import string
import Files_2016
import Files_2017
import Files_2018
import Files_2017_A
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


SAMPLES = {}
#SAMPLES ['DYM10to50'] = ['address', 'data/mc','dataset','year', 'run', 'cross section','lumi','Neventsraw']
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
    SAMPLES.update(Files_2017_A.mc2017_samples)
if data_2017:
    SAMPLES.update(Files_2017_A.data2017_samples)
if mc_2018:
    SAMPLES.update(Files_2018.mc2018_samples)
if data_2018:
    SAMPLES.update(Files_2018.data2018_samples)

rootlib1 = subprocess.check_output("root-config --cflags", shell=True)
rootlib11="".join([s for s in rootlib1.strip().splitlines(True) if s.strip()])
rootlib2 = subprocess.check_output("root-config --glibs", shell=True)
rootlib22="".join([s for s in rootlib2.strip().splitlines(True) if s.strip()])

dire = loc+'bin'
dire_h = loc+'hists/'
print dire
nf =40

for key, value in SAMPLES.items():
#########################################
    #if name  not in key:
    #   continue
    nf = 40
    for idx, S in enumerate(value[0]):
        #print S
        filelist = GFAL_GetROOTfiles( S ,"srm://maite.iihe.ac.be:8443/srm/managerv2?SFN=" )#"srm://ingrid-se02.cism.ucl.ac.be:8444/srm/managerv2?SFN=/storage/data/cms")
        #for files in filelist:
        #for subdir, dirs, files in os.walk(S):
        if value[1]=='data': 
           nf = 255
        sequance = [filelist[i:i+nf] for i in range(0,len(filelist),nf)]
        for num,  seq in enumerate(sequance):
###############################
#                if num<18:
#                    continue
#############################
                text = ''
                text += '    TChain* ch    = new TChain("IIHEAnalysis") ;\n'
                for filename in seq:
                    text += '    ch ->Add("' + filename + '");\n'
                text += '    MyAnalysis t1(ch);\n'
                text += '    t1.Loop("'+dire_h+ value[3] + '/' + key +'_' + str(idx) +'_' +str(num)  + '.root", "' + value[1] + '" , "'+ value[2] + '" , "'+ value[3] + '" , "'+ value[4] + '" , ' + value[5] + ' , '+ value[6] + ' , '+ value[7] + ');\n'
                SHNAME1 = key +'_' + str(idx) +'_' +str(num) + '.C'
                SHFILE1='#include "MyAnalysis.h"\n' +\
                'main(){\n' +\
                text +\
                '}'
                open('Jobs/'+SHNAME1, 'wt').write(SHFILE1)
#                os.system('g++ -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2  -I./../include   '+ rootlib11 +' -ldl  -o ' + SHNAME1.split('.')[0] + ' ' + SHNAME1+ ' ../lib/main.so ' + rootlib22 + '  -lMinuit -lMinuit2 -lTreePlayer -lGenVector')

                SHNAME = key +'_' + str(idx) +'_' + str(num) +'.sh'
                SHFILE="#!/bin/bash\n" +\
                "cd "+ dire + "\n"+\
                'export X509_USER_PROXY=/afs/cern.ch/user/j/jingyan/x509up_u122075' +"\n"+\
                'g++ -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2  -I./../include   '+ rootlib11 +' -ldl  -o ' + SHNAME1.split('.')[0] + ' Jobs/' + SHNAME1+ ' ../lib/main.so ' + rootlib22 + '  -lMinuit -lTreePlayer' + "\n"+\
                "./" + SHNAME1.split('.')[0]+ "\n"+\
                'FILE='+ dire_h + value[3] + '/' + key +'_' + str(idx) +'_' +str(num)  + '.root'+ "\n"+\
                'if [ -f "$FILE" ]; then'+ "\n"+\
                '    rm  ' + SHNAME1.split('.')[0] + "\n"+\
                'fi'
                open('Jobs/'+SHNAME, 'wt').write(SHFILE)
                os.system("chmod +x "+'Jobs/'+SHNAME)
#                os.system("qsub -q localgrid  -o STDOUT/" + SHNAME1.split('.')[0] + ".stdout -e STDERR/" + SHNAME1.split('.')[0] + ".stderr " + SHNAME)
           #break
    if verbose : 
        print key + ' jobs are made'
   
 
