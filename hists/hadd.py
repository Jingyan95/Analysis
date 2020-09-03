import sys
import os
import subprocess
import readline
import string


sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'bin'))

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

SAMPLES.update(Files_2017_A.mc2017_samples)
SAMPLES.update(Files_2017_A.data2017_samples)

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


addedFilesData = {"2017": []}
addedFilesMc = {"2017": []}
addedFilesTTV = {"2017": []}
addedFilesWZ = {"2017": []}
addedFilesZZ = {"2017": []}
addedFilesTTbar = {"2017": []}

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
    elif ('SMEFTfr' not in key):
        addedFilesMc[year].append(  year + '/' + key + '.root ')
    else:
        hadd='hadd ' + key + '.root '
    for idx, S in enumerate(value[0]):
        if value[1]=='data':
            nf = 255
        for subdir, dirs, files in os.walk(S):
            sequance = [files[i:i+nf] for i in range(0,len(files),nf)]
            for num,  seq in enumerate(sequance):
                hadd +=  year + '/' + key +'_' + str(idx) +'_' + str(num) + '.root '
            break
    os.system(hadd)

os.system('rm *_data.root')
os.system('rm *_others.root')
os.system('rm *_TTV.root')
os.system('rm *_WZ.root')
os.system('rm *_ZZ.root')
os.system('rm *_TTbar.root')
#hadddata_2016 ='hadd 2016_data' + '.root ' + ' '.join(addedFilesData['2016'])
hadddata_2017 ='hadd 2017_data' + '.root ' + ' '.join(addedFilesData['2017'])
#hadddata_2018 ='hadd 2018_data' + '.root ' + ' '.join(addedFilesData['2018'])

#haddmc_2016 ='hadd 2016_others' + '.root ' + ' '.join(addedFilesMc['2016'])
haddmc_2017 ='hadd 2017_others' + '.root ' + ' '.join(addedFilesMc['2017'])
#haddmc_2018 ='hadd 2018_others' + '.root ' + ' '.join(addedFilesMc['2018'])
haddTTV_2017 ='hadd 2017_TTV' + '.root ' + ' '.join(addedFilesTTV['2017'])
haddWZ_2017 ='hadd 2017_WZ' + '.root ' + ' '.join(addedFilesWZ['2017'])
haddZZ_2017 ='hadd 2017_ZZ' + '.root ' + ' '.join(addedFilesZZ['2017'])
haddTTbar_2017 ='hadd 2017_TTbar' + '.root ' + ' '.join(addedFilesTTbar['2017'])

#os.system(haddmc_2016)
os.system(haddmc_2017)
#os.system(haddmc_2018)

#os.system(hadddata_2016)
os.system(hadddata_2017)
#os.system(hadddata_2018)
os.system(haddTTV_2017)
os.system(haddWZ_2017)
os.system(haddZZ_2017)
os.system(haddTTbar_2017)

#os.system('rm *_DY.root')
#os.system('hadd 2016_DY.root 2016/2016_DYM50.root 2016/2016_DYM10to50.root')
#os.system('hadd 2017_DY.root 2017/2017_DYM50.root 2017/2017_DYM10to50.root')
#os.system('hadd 2018_DY.root 2018/2018_DYM50.root 2018/2018_DYM10to50.root')

#os.system('rm *tW.root')
#os.system('hadd 2016_ST_tW.root 2016/2016_ST_tW.root 2016/2016_ST_atW.root')
#os.system('hadd 2017_ST_tW.root 2017/2017_ST_tW.root 2017/2017_ST_atW.root')
#os.system('hadd 2018_ST_tW.root 2018/2018_ST_tW.root 2018/2018_ST_atW.root')

#os.system('rm *_LFV*')
#os.system('hadd 2016_LFVVecC.root 2017/2016_LFVTtVecC.root 2017/2016_LFVStVecC.root')
#os.system('hadd 2017_LFVVecC.root 2017/2017_LFVTtVecC.root 2017/2017_LFVStVecC.root')
#os.system('hadd 2018_LFVVecC.root 2017/2018_LFVTtVecC.root 2017/2018_LFVStVecC.root')
#
#os.system('hadd 2016_LFVVecU.root 2017/2016_LFVTtVecU.root 2017/2016_LFVStVecU.root')
#os.system('hadd 2017_LFVVecU.root 2017/2017_LFVTtVecU.root 2017/2017_LFVStVecU.root')
#os.system('hadd 2018_LFVVecU.root 2017/2018_LFVTtVecU.root 2017/2018_LFVStVecU.root')


#Y = ['2017']
#Sam = ['TTTo2L2Nu','WJetsToLNu']
#
#for y in Y:
#    for s in Sam:
#        txt = 'cp ' + y + '/' + y + '_' + s + '.root .'
#        os.system(txt)
#
#allSamples = ['TTTo2L2Nu','ST_tW','WJetsToLNu','data', 'others', 'DY', 'SMEFTfr_ST_vector_emutc', 'SMEFTfr_ST_vector_emutu', 'SMEFTfr_TT_vector_emutc', 'SMEFTfr_TT_vector_emutu', 'SMEFTfr_ST_clequ1_emutc', 'SMEFTfr_ST_clequ3_emutu', 'SMEFTfr_ST_clequ3_emutc', 'SMEFTfr_TT_clequ1_emutc', 'SMEFTfr_TT_clequ1_emutu', 'SMEFTfr_TT_clequ3_emutu']
#allMCSamples = ['TTTo2L2Nu','ST_tW','WJetsToLNu', 'others', 'DY']
#
#for s in allSamples:
#    os.system('rm All_' + s + '.root ')
#    haddall='hadd All_' + s + '.root '
#    for y in Y:
#        haddall +=  y + '_' + s + '.root '
#    os.system(haddall)
#
#for y in Y:
#    os.system('rm mc_' + y + '.root ')
#    haddall='hadd mc_' + y + '.root '
#    for s in allMCSamples:
#        haddall +=  y + '_' + s + '.root '
#    os.system(haddall)
