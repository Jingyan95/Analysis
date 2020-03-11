import sys
import os
import subprocess
import readline
import string


SAMPLES = {}
#SAMPLES ['DYM10to50'] = ['address', 'data/mc','dataset','year', 'run', 'cross section','lumi','Neventsraw']
mc_2016 = True
data_2016 = True

#2016 MC
if mc_2016:
    SAMPLES['2016_DYM10to50'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3/191130_081246/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3_ext1/191130_081323/0000/'], 'mc','','2016', '','18610','35.92','78843820']
    SAMPLES['2016_DYM50'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M50_v3_ext1/191130_075804/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_v3_ext2/191130_081130/0000/'], 'mc','','2016', '','5765.4','35.92','146280395']
    SAMPLES['2016_TTTo2L2Nu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/crab_TTTo2L2Nu/191130_081559/0000/'], 'mc','','2016', '','87.31','35.92','79140880']
    SAMPLES['2016_ST_tW'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_antitop_v3_ext1/191130_092526/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/crab_ST_tW_top_v3_ext1/191130_093520/0000/'], 'mc','','2016', '','38.94','35.92','6368069']
    SAMPLES['2016_WWTo2L2Nu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WWTo2L2Nu_13TeV-powheg/crab_WWTo2L2Nu/191130_090952/0000/'], 'mc','','2016', '','12.178','35.92','1999000']
    SAMPLES['2016_ZZTo2L2Nu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ZZTo2L2Nu_13TeV_powheg_pythia8/crab_ZZTo2L2Nu/191130_091725/0000/'], 'mc','','2016', '','0.564','35.92','8931750']
    SAMPLES['2016_ZZTo4L'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/ZZTo4L_13TeV_powheg_pythia8/crab_ZZTo4L/191130_092020/0000/'], 'mc','','2016', '','1.212','35.92','6669988']
    SAMPLES['2016_WZTo2L2Q'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo2L2Q/191130_091455/0000/'], 'mc','','2016', '','5.595','35.92','15879472']
    SAMPLES['2016_WZTo3LNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/crab_WZTo3LNu_v3_ext1_v1/191201_062001/0000/'], 'mc','','2016', '','4.43','35.92','18000000']
    SAMPLES['2016_WJetsToLNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_v3_ext2_v2/191130_113603/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu_v3_v2/191130_115257/0000/'], 'mc','','2016', '','61526.7','35.92','86916455']
    SAMPLES['2016_TTWJetsToQQ'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToQQ/191130_110154/0000/'], 'mc','','2016', '','0.4062','35.92','430310']
    SAMPLES['2016_TTWJetsToLNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_v3_ext1_v2/191201_070704/0000/','/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu/191130_110349/0000/'], 'mc','','2016', '','0.2043','35.92','2716249']
    SAMPLES['2016_TTZToLLNuNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuN/191130_110706/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_v3_ext1_v2/191201_071154/0000/', '/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_v3_ext3_v1/191201_071416/0000/'], 'mc','','2016', '','0.2529','35.92','6420825']
    SAMPLES['2016_TTZToQQ'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTZToQQ/191130_110848/0000/'], 'mc','','2016', '','0.5297','35.92','351164']

#data 2016
if data_2016:
    SAMPLES['2016_B_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runB/191203_070326/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runB/191203_070326/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runB/191203_070326/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runB/191203_070326/0003/'], 'data','MuonEG','2016', 'B','1','1','1']
    SAMPLES['2016_C_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runC/191203_070345/0001/'], 'data','MuonEG','2016', 'C','1','1','1']
    SAMPLES['2016_D_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runD/191203_070404/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runD/191203_070404/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runD/191203_070404/0002/'], 'data','MuonEG','2016', 'D','1','1','1']
    SAMPLES['2016_E_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runE/191203_071509/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runE/191203_071509/0001/'], 'data','MuonEG','2016', 'E','1','1','1']
    SAMPLES['2016_F_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runF/191203_070443/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runF/191203_070443/0001/'], 'data','MuonEG','2016', 'F','1','1','1']
    SAMPLES['2016_G_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runG/191203_070502/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runG/191203_070502/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runG/191203_070502/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runG/191203_070502/0003/'], 'data','MuonEG','2016', 'G','1','1','1']
    SAMPLES['2016_H_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runH/191203_070522/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runH/191203_070522/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runH/191203_070522/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/MuonEG/crab_MuonEG_runH/191203_070522/0003/'], 'data','MuonEG','2016', 'H','1','1','1']

    SAMPLES['2016_B_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runB/191203_065806/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runB/191203_065806/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runB/191203_065806/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runB/191203_065806/0003/'], 'data','SingleMuon','2016', 'B','1','1','1']
    SAMPLES['2016_C_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runC/191203_065920/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runC/191203_065920/0001/'], 'data','SingleMuon','2016', 'C','1','1','1']
    SAMPLES['2016_D_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runD/191203_065939/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runD/191203_065939/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runD/191203_065939/0002/'], 'data','SingleMuon','2016', 'D','1','1','1']
    SAMPLES['2016_E_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runD/191203_065939/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runD/191203_065939/0001/'], 'data','SingleMuon','2016', 'E','1','1','1']
    SAMPLES['2016_F_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runF/191203_070017/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runF/191203_070017/0001/'], 'data','SingleMuon','2016', 'F','1','1','1']
    SAMPLES['2016_G_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runG/191203_070036/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runG/191203_070036/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runG/191203_070036/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runG/191203_070036/0003/'], 'data','SingleMuon','2016', 'G','1','1','1']
    SAMPLES['2016_H_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runH/191203_070054/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runH/191203_070054/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runH/191203_070054/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleMuon/crab_SingleMuon_runH/191203_070054/0003/'], 'data','SingleMuon','2016', 'H','1','1','1']

    SAMPLES['2016_B_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runB/191203_070113/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runB/191203_070113/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runB/191203_070113/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runB/191203_070113/0003/'], 'data','SingleElectron','2016', 'B','1','1','1']
    SAMPLES['2016_C_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runC/191203_070132/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runC/191203_070132/0001/'], 'data','SingleElectron','2016', 'C','1','1','1']
    SAMPLES['2016_D_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runD/191203_070151/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runD/191203_070151/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runD/191203_070151/0002/'], 'data','SingleElectron','2016', 'D','1','1','1']
    SAMPLES['2016_E_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runE/191203_070210/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runE/191203_070210/0001/'], 'data','SingleElectron','2016', 'E','1','1','1']
    SAMPLES['2016_F_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runF/191203_070229/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runF/191203_070229/0001/'], 'data','SingleElectron','2016', 'F','1','1','1']
    SAMPLES['2016_G_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runG/191203_070247/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runG/191203_070247/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runG/191203_070247/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runG/191203_070247/0003/'], 'data','SingleElectron','2016', 'G','1','1','1']
    SAMPLES['2016_H_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runH/191203_070307/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runH/191203_070307/0001/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runH/191203_070307/0002/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/SingleElectron/crab_SingleElectron_runH/191203_070307/0003/'], 'data','SingleElectron','2016', 'H','1','1','1']

    SAMPLES['2016_B_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runB/200306_134908/0000/', '/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runB/200306_134908/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runB/200306_134908/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runB/200306_134908/0003/'], 'data','DoubleEG','2016', 'B','1','1','1']
    SAMPLES['2016_C_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runC/200306_134951/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runC/200306_134951/0001/'], 'data','DoubleEG','2016', 'C','1','1','1']
    SAMPLES['2016_D_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runD/200306_135009/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runD/200306_135009/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runD/200306_135009/0002/'], 'data','DoubleEG','2016', 'D','1','1','1']
    SAMPLES['2016_E_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runE/200306_135026/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runE/200306_135026/0001/'], 'data','DoubleEG','2016', 'E','1','1','1']
    SAMPLES['2016_F_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runF/200306_135047/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runF/200306_135047/0001/'], 'data','DoubleEG','2016', 'F','1','1','1']
    SAMPLES['2016_G_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runG/200306_135106/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runG/200306_135106/0001','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runG/200306_135106/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runG/200306_135106/0003/'], 'data','DoubleEG','2016', 'G','1','1','1']
    SAMPLES['2016_H_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runH/200306_135125/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runH/200306_135125/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runH/200306_135125/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleEG/crab_DoubleEG_runH/200306_135125/0003/'], 'data','DoubleEG','2016', 'H','1','1','1']

    SAMPLES['2016_B_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runB/200306_135144/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runB/200306_135144/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runB/200306_135144/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runB/200306_135144/0003/'], 'data','DoubleMu','2016', 'B','1','1','1']
    SAMPLES['2016_C_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runC/200306_135205/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runC/200306_135205/0001/'], 'data','DoubleMu','2016', 'C','1','1','1']
    SAMPLES['2016_D_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runD/200306_135224/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runD/200306_135224/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runD/200306_135224/0002/'], 'data','DoubleMu','2016', 'D','1','1','1']
    SAMPLES['2016_E_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runE/200306_135245/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runE/200306_135245/0001/'], 'data','DoubleMu','2016', 'E','1','1','1']
    SAMPLES['2016_F_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runF/200306_135305/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runF/200306_135305/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runF/200306_135305/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runF/200306_135305/0003/'], 'data','DoubleMu','2016', 'F','1','1','1']
    SAMPLES['2016_G_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runG/200306_135327/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runG/200306_135327/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runG/200306_135327/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runG/200306_135327/0003/'], 'data','DoubleMu','2016', 'G','1','1','1']
    SAMPLES['2016_H_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runH/200306_135349/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runH/200306_135349/0001/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runH/200306_135349/0002/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2016/DoubleMuon/crab_DoubleMuon_runH/200306_135349/0003/'], 'data','DoubleMu','2016', 'H','1','1','1']


addedFilesData = []

addedFilesMc = []

for key, value in SAMPLES.items():
    os.system('rm '+key + '.root')
    nf = 40
#    os.system('rm ' + key + '.root')
    if value[1]=='data':
        addedFilesData.append( key + '.root ')
    else:
        addedFilesMc.append( key + '.root ')
    hadd='hadd ' + key + '.root '
    for idx, S in enumerate(value[0]):
        if value[1]=='data':
            nf = 250
        for subdir, dirs, files in os.walk(S):
            sequance = [files[i:i+nf] for i in range(0,len(files),nf)]
            for num,  seq in enumerate(sequance):
                hadd += key +'_' + str(idx) +'_' + str(num) + '.root '

    os.system(hadd)

STEP1=True
if STEP1:
    os.system('rm *_data.root')
    os.system('rm *_others.root')
    hadddata_2016 ='hadd 2016_data' + '.root '
    hadddata_2017 ='hadd 2017_data' + '.root '
    hadddata_2018 ='hadd 2018_data' + '.root '
    
    haddmc_2016 ='hadd 2016_others' + '.root '
    haddmc_2017 ='hadd 2017_others' + '.root '
    haddmc_2018 ='hadd 2018_others' + '.root '
    
    for num in addedFilesData:
        if num.split("_")[0] == '2016':
            hadddata_2016 += num + ' '
        if num.split("_")[0] == '2017':
            hadddata_2017 += num + ' '
        if num.split("_")[0] == '2018':
            hadddata_2018 += num + ' '
    
    for num in addedFilesMc:
        if ('TTTo2L2Nu' in num) or ('DY' in  num) or ('ST_tW' in  num) or ('WJetsToLNu' in  num):
            continue   
        if num.split("_")[0] == '2016':
            haddmc_2016 += num + ' '
        if num.split("_")[0] == '2017':
            haddmc_2017 += num + ' '
        if num.split("_")[0] == '2018':
            haddmc_2018 += num + ' '
    
    os.system(haddmc_2016)
    os.system(hadddata_2016)
    os.system('rm *_DY.root')
    os.system('hadd 2016_DY.root 2016_DYM50.root 2016_DYM10to50.root')
















