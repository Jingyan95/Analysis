import sys
import os
import subprocess
import readline
import string

data2017_samples = {}
mc2017_samples = {}
#data2017_samples ['DYM10to50'] = ['address', 'data/mc','dataset','year', 'run', 'cross section','lumi','Neventsraw']

mc2017_samples['2017_DYM10to50'] = [['/eos/cms/store/user/asparker/TopLFV/2017_DYM10to50_ext0/', '/eos/cms/store/user/asparker/TopLFV/2017_DYM10to50_ext1/'], 'mc','','2017', '','18610','41.53','78994955']
mc2017_samples['2017_DYM50'] = [['/eos/user/a/asparker/TopLFV/2017_DYM50_ext0/', '/eos/user/a/asparker/TopLFV/2017_DYM50_ext1/', '/eos/user/a/asparker/TopLFV/2017_DYM50_ext2/'], 'mc','','2017', '','5765.4','41.53','84597528']
mc2017_samples['2017_TTTo2L2Nu'] = [['/eos/user/a/asparker/TopLFV/2017_TTTo2L2Nu/'], 'mc','','2017', '','87.31','41.53','8926992']
mc2017_samples['2017_ST_tW'] = [['/eos/user/a/asparker/TopLFV/2017_ST_tW/'], 'mc','','2017', '','19.47','41.53','4917240']
mc2017_samples['2017_ST_atW'] = [['/eos/user/a/asparker/TopLFV/2017_ST_atW/'], 'mc','','2017', '','19.47','41.53','5592819']
mc2017_samples['2017_WWTo2L2Nu'] = [['/eos/user/a/asparker/TopLFV/2017_WWTo2L2Nu/'], 'mc','','2017', '','12.178','41.53','1992522']
mc2017_samples['2017_ZZTo2L2Nu'] = [['/eos/user/a/asparker/TopLFV/2017_ZZTo2L2Nu/'], 'mc','','2017', '','0.564','41.53','8733658']
mc2017_samples['2017_ZZTo4L'] = [['/eos/user/a/asparker/TopLFV/2017_ZZTo4L/'], 'mc','','2017', '','1.212','41.53','6893887']
mc2017_samples['2017_WZTo2L2Q'] = [['/eos/user/a/asparker/TopLFV/2017_WZTo2L2Q/'], 'mc','','2017', '','5.595','41.53','16620982']
mc2017_samples['2017_WZTo3LNu'] = [['/eos/user/a/asparker/TopLFV/2017_WZTo3LNu/'], 'mc','','2017', '','4.43','41.53','6887413']
mc2017_samples['2017_WJetsToLNu'] = [['/eos/user/a/asparker/TopLFV/2017_WJetsToLNu/'], 'mc','','2017', '','61526.7','41.53','33043732']
mc2017_samples['2017_TTWJetsToQQ'] = [['/eos/user/a/asparker/TopLFV/2017_TTWJetsToQQ/'], 'mc','','2017', '','0.4062','41.53','441560']
mc2017_samples['2017_TTWJetsToLNu'] = [['/eos/user/a/asparker/TopLFV/2017_TTWJetsToLNu/'], 'mc','','2017', '','0.2043','41.53','2678775']
mc2017_samples['2017_TTZToLLNuNu'] = [['/eos/user/a/asparker/TopLFV/2017_TTZToLLNuNu/'], 'mc','','2017', '','0.2529','41.53','3570720']
mc2017_samples['2017_TTZToQQ'] = [['/eos/user/a/asparker/TopLFV/2017_TTZToQQ/'], 'mc','','2017', '','0.5297','41.53','356286']

#mc2017_samples['2017_LFVStVecC'] = [['/eos/cms/store/user/asparker/TopLFV/2017_LFVStVecC/'], 'mc','','2017', '','0.0512' ,'41.53','500000']
#vector interaction
mc2017_samples['2017_SMEFTfr_ST_vector_emutc'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_ST_vector_emutc/'], 'mc','','2017', '','0.04' ,'41.53','500000']
mc2017_samples['2017_SMEFTfr_ST_vector_emutu'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_ST_vector_emutu/'], 'mc','','2017', '','0.452' ,'41.53','494000']
mc2017_samples['2017_SMEFTfr_TT_vector_emutc'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_TT_vector_emutc/'], 'mc','','2017', '','0.032'  ,'41.53','500000']
mc2017_samples['2017_SMEFTfr_TT_vector_emutu'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_TT_vector_emutu/'], 'mc','','2017', '','0.032','41.53','498000']

#tensor interaction
mc2017_samples['2017_SMEFTfr_ST_clequ3_emutu'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_ST_clequ3_emutu/'], 'mc','','2017', '','1.9' ,'41.53','500000']
mc2017_samples['2017_SMEFTfr_ST_clequ3_emutc'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_ST_clequ3_emutc/'], 'mc','','2017', '','0.187' ,'41.53','500000']
mc2017_samples['2017_SMEFTfr_TT_clequ3_emutu'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_TT_clequ3_emutu/'], 'mc','','2017', '','0.1876','41.53','500000']
mc2017_samples['2017_SMEFTfr_TT_clequ3_emutc'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_TT_clequ3_emutc/'], 'mc','','2017', '','0.1876','41.53','500000']

#scalar interaction
mc2017_samples['2017_SMEFTfr_ST_clequ1_emutc'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_ST_clequ1_emutc/'], 'mc','','2017', '','0.008' ,'41.53','500000']
mc2017_samples['2017_SMEFTfr_ST_clequ1_emutu'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_ST_clequ1_emutu/'], 'mc','','2017', '','0.102' ,'41.53','500000']
mc2017_samples['2017_SMEFTfr_TT_clequ1_emutc'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_TT_clequ1_emutc/'], 'mc','','2017', '','0.004'  ,'41.53','500000']
mc2017_samples['2017_SMEFTfr_TT_clequ1_emutu'] = [['/eos/cms/store/user/asparker/TopLFV/2017_SMEFTfr_TT_clequ1_emutu/'], 'mc','','2017', '','0.004','41.53','500000']



data2017_samples['2017_B_MuonEG'] = [['/eos/user/a/asparker/TopLFV/2017_B_MuonEG/'], 'data','MuonEG','2017', 'B','1','1','1']
data2017_samples['2017_C_MuonEG'] = [['/eos/user/a/asparker/TopLFV/2017_C_MuonEG_ext0/', '/eos/user/a/asparker/TopLFV/2017_C_MuonEG_ext1/'], 'data','MuonEG','2017', 'C','1','1','1']
data2017_samples['2017_D_MuonEG'] = [['/eos/user/a/asparker/TopLFV/2017_D_MuonEG/'], 'data','MuonEG','2017', 'D','1','1','1']
data2017_samples['2017_E_MuonEG'] = [['/eos/user/a/asparker/TopLFV/2017_E_MuonEG_ext0/', '/eos/user/a/asparker/TopLFV/2017_E_MuonEG_ext1/'], 'data','MuonEG','2017', 'E','1','1','1']
data2017_samples['2017_F_MuonEG'] = [['/eos/user/a/asparker/TopLFV/2017_F_MuonEG_ext0/', '/eos/user/a/asparker/TopLFV/2017_F_MuonEG_ext1/'], 'data','MuonEG','2017', 'F','1','1','1']

data2017_samples['2017_B_SingleMuon'] = [['/eos/user/a/asparker/TopLFV/2017_B_SingleMuon/'], 'data','SingleMuon','2017', 'B','1','1','1']
data2017_samples['2017_C_SingleMuon'] = [['/eos/user/a/asparker/TopLFV/2017_C_SingleMuon_ext0/', '/eos/cms/store/user/asparker/TopLFV/2017_C_SingleMuon_ext1/'], 'data','SingleMuon','2017', 'C','1','1','1']
data2017_samples['2017_D_SingleMuon'] = [['/eos/user/a/asparker/TopLFV/2017_D_SingleMuon/'], 'data','SingleMuon','2017', 'D','1','1','1']
data2017_samples['2017_E_SingleMuon'] = [['/eos/user/a/asparker/TopLFV/2017_E_SingleMuon_ext0/', '/eos/user/a/asparker/TopLFV/2017_E_SingleMuon_ext1/'], 'data','SingleMuon','2017', 'E','1','1','1']
data2017_samples['2017_F_SingleMuon'] = [['/eos/user/a/asparker/TopLFV/2017_F_SingleMuon_ext0/', '/eos/user/a/asparker/TopLFV/2017_F_SingleMuon_ext1/'], 'data','SingleMuon','2017', 'F','1','1','1']

data2017_samples['2017_B_SingleElectron'] = [['/eos/user/a/asparker/TopLFV/2017_B_SingleElectron/'], 'data','SingleElectron','2017', 'B','1','1','1']
data2017_samples['2017_C_SingleElectron'] = [['/eos/cms/store/user/asparker/TopLFV/2017_C_SingleElectron_ext0/', '/eos/cms/store/user/asparker/TopLFV/2017_C_SingleElectron_ext1/'], 'data','SingleElectron','2017', 'C','1','1','1']
data2017_samples['2017_D_SingleElectron'] = [['/eos/user/a/asparker/TopLFV/2017_D_SingleElectron/'], 'data','SingleElectron','2017', 'D','1','1','1']
data2017_samples['2017_E_SingleElectron'] = [['/eos/user/a/asparker/TopLFV/2017_E_SingleElectron_ext0/', '/eos/user/a/asparker/TopLFV/2017_E_SingleElectron_ext1/'], 'data','SingleElectron','2017', 'E','1','1','1']
data2017_samples['2017_F_SingleElectron'] = [['/eos/user/a/asparker/TopLFV/2017_F_SingleElectron_ext0/', '/eos/user/a/asparker/TopLFV/2017_F_SingleElectron_ext1/'], 'data','SingleElectron','2017', 'F','1','1','1']

data2017_samples['2017_B_DoubleEG'] = [['/eos/user/a/asparker/TopLFV/2017_B_DoubleEG/'], 'data','DoubleEG','2017', 'B','1','1','1']
data2017_samples['2017_C_DoubleEG'] = [['/eos/user/a/asparker/TopLFV/2017_C_DoubleEG_ext0/', '/eos/user/a/asparker/TopLFV/2017_C_DoubleEG_ext1/'], 'data','DoubleEG','2017', 'C','1','1','1']
data2017_samples['2017_D_DoubleEG'] = [['/eos/cms/store/user/asparker/TopLFV/2017_D_DoubleEG/'], 'data','DoubleEG','2017', 'D','1','1','1']
data2017_samples['2017_E_DoubleEG'] = [['/eos/user/a/asparker/TopLFV/2017_E_DoubleEG_ext0/', '/eos/user/a/asparker/TopLFV/2017_E_DoubleEG_ext1/'], 'data','DoubleEG','2017', 'E','1','1','1']
data2017_samples['2017_F_DoubleEG'] = [['/eos/user/a/asparker/TopLFV/2017_F_DoubleEG_ext0/', '/eos/user/a/asparker/TopLFV/2017_F_DoubleEG_ext1/'], 'data','DoubleEG','2017', 'F','1','1','1']

data2017_samples['2017_B_DoubleMu'] = [['/eos/cms/store/user/asparker/TopLFV/2017_B_DoubleMu/'], 'data','DoubleMu','2017', 'B','1','1','1']
data2017_samples['2017_C_DoubleMu'] = [['/eos/user/a/asparker/TopLFV/2017_C_DoubleMu_ext0/', '/eos/user/a/asparker/TopLFV/2017_C_DoubleMu_ext1/'], 'data','DoubleMu','2017', 'C','1','1','1']
data2017_samples['2017_D_DoubleMu'] = [['/eos/user/a/asparker/TopLFV/2017_D_DoubleMu/'], 'data','DoubleMu','2017', 'D','1','1','1']
data2017_samples['2017_E_DoubleMu'] = [['/eos/user/a/asparker/TopLFV/2017_E_DoubleMu_ext0/', '/eos/user/a/asparker/TopLFV/2017_E_DoubleMu_ext1/'], 'data','DoubleMu','2017', 'E','1','1','1']
data2017_samples['2017_F_DoubleMu'] = [['/eos/user/a/asparker/TopLFV/2017_F_DoubleMu_ext0/', '/eos/user/a/asparker/TopLFV/2017_F_DoubleMu_ext1/'], 'data','DoubleMu','2017', 'F','1','1','1']


