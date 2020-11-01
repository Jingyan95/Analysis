import sys
import os
import subprocess
import readline
import string


sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'bin'))

import Files_2016
import Files_2017
import Files_2018
from GFAL_GetROOTfiles import *

SAMPLES = {}
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


addedFilesData = {"2016": [], "2017": [], "2018": []}
addedFilesMc = {"2016": [], "2017": [], "2018": []}
addedFilesTTV = {"2016": [], "2017": [], "2018": []}
addedFilesWZ = {"2016": [], "2017": [], "2018": []}
addedFilesZZ = {"2016": [], "2017": [], "2018": []}
addedFilesTTbar = {"2016": [], "2017": [], "2018": []}

for key, value in SAMPLES.items():
    year = value[3]
    os.system('rm '+ key + '.root')
    os.system('rm '+ year + '/' +key + '.root')
    nf = 40
    hadd='hadd ' + year + '/' + key + '.root '
    if value[1]=='data':
        addedFilesData[year].append( year + '/' + key + '.root ')
    elif ('TTWJetsToLNu' in key) or ('TTZToLLNuNu' in  key):
        addedFilesTTV[year].append(  year + '/' + key + '.root ')
    elif ('WZTo3LNu' in key):
        addedFilesWZ[year].append(  year + '/' + key + '.root ')
    elif ('ZZTo4L' in key):
        addedFilesZZ[year].append(  year + '/' + key + '.root ')
    elif ('TTTo2L2Nu' in key):
        addedFilesTTbar[year].append(  year + '/' + key + '.root ')
    elif ('LFV' not in key):
        addedFilesMc[year].append(  year + '/' + key + '.root ')
    else:
        hadd='hadd ' + key + '.root '
    for idx, S in enumerate(value[0]):
        filelist = GFAL_GetROOTfiles( S ,"srm://maite.iihe.ac.be:8443/srm/managerv2?SFN="
        )
        if value[1]=='data':
            nf = 255
        #for subdir, dirs, files in os.walk(S):
        sequance = [filelist[i:i+nf] for i in range(0,len(filelist),nf)]
        for num,  seq in enumerate(sequance):
            hadd +=  year + '/' + key +'_' + str(idx) +'_' + str(num) + '.root ' 
    os.system(hadd)

os.system('rm *_data.root')
os.system('rm *_others.root')
os.system('rm *_TTV.root')
os.system('rm *_WZ.root')
os.system('rm *_ZZ.root')
os.system('rm *_TTbar.root')
hadddata_2016 ='hadd 2016_data' + '.root ' + ' '.join(addedFilesData['2016'])
hadddata_2017 ='hadd 2017_data' + '.root ' + ' '.join(addedFilesData['2017'])
hadddata_2018 ='hadd 2018_data' + '.root ' + ' '.join(addedFilesData['2018'])

haddmc_2016 ='hadd 2016_others' + '.root ' + ' '.join(addedFilesMc['2016'])
haddmc_2017 ='hadd 2017_others' + '.root ' + ' '.join(addedFilesMc['2017'])
haddmc_2018 ='hadd 2018_others' + '.root ' + ' '.join(addedFilesMc['2018'])

haddTTV_2016 ='hadd 2016_TTV' + '.root ' + ' '.join(addedFilesTTV['2016'])
haddTTV_2017 ='hadd 2017_TTV' + '.root ' + ' '.join(addedFilesTTV['2017'])
haddTTV_2018 ='hadd 2018_TTV' + '.root ' + ' '.join(addedFilesTTV['2018'])

haddWZ_2016 ='hadd 2016_WZ' + '.root ' + ' '.join(addedFilesWZ['2016'])
haddWZ_2017 ='hadd 2017_WZ' + '.root ' + ' '.join(addedFilesWZ['2017'])
haddWZ_2018 ='hadd 2018_WZ' + '.root ' + ' '.join(addedFilesWZ['2018'])

haddZZ_2016 ='hadd 2016_ZZ' + '.root ' + ' '.join(addedFilesZZ['2016'])
haddZZ_2017 ='hadd 2017_ZZ' + '.root ' + ' '.join(addedFilesZZ['2017'])
haddZZ_2018 ='hadd 2018_ZZ' + '.root ' + ' '.join(addedFilesZZ['2018'])

haddTTbar_2016 ='hadd 2016_TTbar' + '.root ' + ' '.join(addedFilesTTbar['2016'])
haddTTbar_2017 ='hadd 2017_TTbar' + '.root ' + ' '.join(addedFilesTTbar['2017'])
haddTTbar_2018 ='hadd 2018_TTbar' + '.root ' + ' '.join(addedFilesTTbar['2018'])

os.system(haddmc_2016)
os.system(haddmc_2017)
os.system(haddmc_2018)

os.system(hadddata_2016)
os.system(hadddata_2017)
os.system(hadddata_2018)

os.system(haddTTV_2016)
os.system(haddTTV_2017)
os.system(haddTTV_2018)

os.system(haddWZ_2016)
os.system(haddWZ_2017)
os.system(haddWZ_2018)

os.system(haddZZ_2016)
os.system(haddZZ_2017)
os.system(haddZZ_2018)

os.system(haddTTbar_2016)
os.system(haddTTbar_2017)
os.system(haddTTbar_2018)

Y = ['2016','2017','2018']
allSamples = ['data','TTV','WZ', 'ZZ', 'TTbar', 'others', 'LFVStVecU', 'LFVTtVecU']

for s in allSamples:
    os.system('rm All_' + s + '.root ')
    haddall='hadd All_' + s + '.root '
    for y in Y:
        haddall +=  y + '_' + s + '.root '
    os.system(haddall)

