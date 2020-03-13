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
data_2016 = False
mc_2017 = True
data_2017 = False
mc_2018 = True
data_2018 = False

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

rootlib1 = subprocess.check_output("root-config --cflags", shell=True)
rootlib11="".join([s for s in rootlib1.strip().splitlines(True) if s.strip()])
rootlib2 = subprocess.check_output("root-config --glibs", shell=True)
rootlib22="".join([s for s in rootlib2.strip().splitlines(True) if s.strip()])

dire = '/user/rgoldouz/NewAnalysis2020/Analysis/bin'
cms = '/user/rgoldouz/CMSSW_9_3_4/src/'
nf =40

for key, value in SAMPLES.items():
    nf = 40
    for idx, S in enumerate(value[0]):
        for subdir, dirs, files in os.walk(S):
            if value[1]=='data': 
                nf = 250
            sequance = [files[i:i+nf] for i in range(0,len(files),nf)]
            for num,  seq in enumerate(sequance):
                if key!='2017_DYM10to50' or num<11:
                    continue
                text = ''
                text += '    TChain* ch    = new TChain("IIHEAnalysis") ;\n'
                for filename in seq:
                    text += '    ch ->Add("' + S+ filename + '");\n'
                text += '    MyAnalysis t1(ch);\n'
                text += '    t1.Loop("/user/rgoldouz/NewAnalysis2020/Analysis/hists/' + value[3] + '/' + key +'_' + str(idx) +'_' +str(num)  + '.root", "' + value[1] + '" , "'+ value[2] + '" , "'+ value[3] + '" , "'+ value[4] + '" , ' + value[5] + ' , '+ value[6] + ' , '+ value[7] + ');\n'
                SHNAME1 = key +'_' + str(idx) +'_' +str(num) + '.C'
                SHFILE1='#include "MyAnalysis.h"\n' +\
                'main(){\n' +\
                text +\
                '}'
                open(SHNAME1, 'wt').write(SHFILE1)
                os.system('g++ -fPIC -fno-var-tracking -Wno-deprecated -D_GNU_SOURCE -O2  -I./../include   '+ rootlib11 +' -ldl  -o ' + SHNAME1.split('.')[0] + ' ' + SHNAME1+ ' ../lib/main.so ' + rootlib22 + '  -lMinuit -lMinuit2 -lTreePlayer -lGenVector')

                SHNAME = key +'_' + str(idx) +'_' + str(num) +'.sh'
                SHFILE="#!/bin/bash\n" +\
                "cd "+ cms + "\n"+\
                "eval `scramv1 runtime -sh`\n"+\
                "cd "+ dire + "\n"+\
                "./" + SHNAME1.split('.')[0]+ "\n"+\
                'FILE='+'/user/rgoldouz/NewAnalysis2020/Analysis/hists/' + value[3] + '/' + key +'_' + str(idx) +'_' +str(num)  + '.root'+ "\n"+\
                'if [ -f "$FILE" ]; then'+ "\n"+\
                '    rm  ' + SHNAME1.split('.')[0] + '*'+ "\n"+\
                'fi'
                open(SHNAME, 'wt').write(SHFILE)
                os.system("chmod +x "+SHNAME)
                os.system("qsub -q localgrid  -o STDOUT/" + SHNAME1.split('.')[0] + ".stdout -e STDERR/" + SHNAME1.split('.')[0] + ".stderr " + SHNAME)
    print key + ' jobs are submitted'
   
 
