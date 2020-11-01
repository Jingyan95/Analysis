import sys
import os
import subprocess
import readline
import string

data2017_samples = {}
mc2017_samples = {}
#data2017_samples ['DYM10to50'] = ['address', 'data/mc','dataset','year', 'run', 'cross section','lumi','Neventsraw']

mc2017_samples['2017_DYM10to50'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-10to50_v14-v1/191201_123859/0000/','/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-10to50_v14_ext1_v2/191201_123656/0000/'], 'mc','','2017', '','18610','41.53','78994955']
mc2017_samples['2017_DYM50'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50_v14_ext1-v1/191201_122119/0000/','/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50_v14_ext1-v1/191201_122119/0001/','/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50_v14_ext1-v1/191201_122119/0002/'], 'mc','','2017', '','5765.4','41.53','84597528']
mc2017_samples['2017_TTTo2L2Nu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/crab_TTTo2L2Nu_v14_v1/191201_124029/0000/'], 'mc','','2017', '','87.31','41.53','8926992']
mc2017_samples['2017_ST_tW'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/crab_ST_tW_top/191201_124931/0000/'], 'mc','','2017', '','19.47','41.53','4917240']
mc2017_samples['2017_ST_atW'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/crab_ST_tW_antitop/191201_125119/0000/'], 'mc','','2017', '','19.47','41.53','5592819']
mc2017_samples['2017_WWTo2L2Nu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/crab_WWTo2L2Nu/191201_124202/0000/'], 'mc','','2017', '','12.178','41.53','1992522']
mc2017_samples['2017_ZZTo2L2Nu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/ZZTo2L2Nu_13TeV_powheg_pythia8/crab_ZZTo2L2Nu/191201_124641/0000/'], 'mc','','2017', '','0.564','41.53','8733658']
mc2017_samples['2017_ZZTo4L'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/ZZTo4L_13TeV_powheg_pythia8/crab_ZZTo4L/191201_124808/0000/'], 'mc','','2017', '','1.212','41.53','6893887']
mc2017_samples['2017_WZTo2L2Q'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/crab_WZTo2L2Q/191201_124509/0000/'], 'mc','','2017', '','5.595','41.53','16620982']
mc2017_samples['2017_WZTo3LNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/crab_WZTo3LNu/191201_124341/0000/'], 'mc','','2017', '','4.43','41.53','6887413']
mc2017_samples['2017_WJetsToLNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/crab_WJetsToLNu/191201_132144/0000/'], 'mc','','2017', '','61526.7','41.53','33043732']
mc2017_samples['2017_TTWJetsToQQ'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToQQ_v14_v1/191201_131707/0000/'], 'mc','','2017', '','0.4062','41.53','441560']
mc2017_samples['2017_TTWJetsToLNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTWJetsToLNu_v14_v1/191201_125315/0000/'], 'mc','','2017', '','0.2043','41.53','2678775']
mc2017_samples['2017_TTZToLLNuNu'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/crab_TTZToLLNuNu_v14_v1/191201_131835/0000/'], 'mc','','2017', '','0.2529','41.53','3570720']
mc2017_samples['2017_TTZToQQ'] = [['/pnfs/iihe/cms/store/user/schenara/MC_RunII_2017/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/crab_TTZToQQ_v14_v1/191201_132958/0000/'], 'mc','','2017', '','0.5297','41.53','356286']

mc2017_samples['2017_LFVStVecC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_vector_emutc/'], 'mc','','2017', '','0.04' ,'41.53','500000']
mc2017_samples['2017_LFVStVecU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_vector_emutu/'], 'mc','','2017', '','0.452' ,'41.53','494000']
mc2017_samples['2017_LFVTtVecC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_TT_vector_emutc/'], 'mc','','2017', '','0.032'  ,'41.53','500000']
mc2017_samples['2017_LFVTtVecU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_TT_vector_emutu/'], 'mc','','2017', '','0.032','41.53','498000']

mc2017_samples['2017_LFVStScalarC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ1_emutc/'], 'mc','','2017', '','0.008' ,'41.53','500000']
mc2017_samples['2017_LFVStScalarU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ1_emutu/'], 'mc','','2017', '','0.102' ,'41.53','500000']
mc2017_samples['2017_LFVTtScalarC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_TT_clequ1_emutc/'], 'mc','','2017', '','0.004'  ,'41.53','500000']
mc2017_samples['2017_LFVTtScalarU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_TT_clequ1_emutu/'], 'mc','','2017', '','0.004','41.53','500000']

mc2017_samples['2017_LFVStTensorC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ3_emutc/'], 'mc','','2017', '','0.187' ,'41.53','500000']
mc2017_samples['2017_LFVStTensorU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_ST_clequ3_emutu/'], 'mc','','2017', '','1.9' ,'41.53','500000']
mc2017_samples['2017_LFVTtTensorC'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_TT_clequ3_emutc/'], 'mc','','2017', '','0.1876'  ,'41.53','500000']
mc2017_samples['2017_LFVTtTensorU'] = [['/pnfs/iihe/cms/store/user/rgoldouz/TopLfvFullSim/2017/IIHE_Ntuple/ntuple_SMEFTfr_TT_clequ3_emutu/'], 'mc','','2017', '','0.1876','41.53','500000']

data2017_samples['2017_B_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/MuonEG/crab_MuonEG_runB/191204_020507/0000/'], 'data','MuonEG','2017', 'B','1','1','1']
data2017_samples['2017_C_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/MuonEG/crab_MuonEG_runC/191204_020526/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/MuonEG/crab_MuonEG_runC/191204_020526/0001/'], 'data','MuonEG','2017', 'C','1','1','1']
data2017_samples['2017_D_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/MuonEG/crab_MuonEG_runD/191204_020545/0000/'], 'data','MuonEG','2017', 'D','1','1','1']
data2017_samples['2017_E_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/MuonEG/crab_MuonEG_runE/191204_020606/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/MuonEG/crab_MuonEG_runE/191204_020606/0001/'], 'data','MuonEG','2017', 'E','1','1','1']
data2017_samples['2017_F_MuonEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/MuonEG/crab_MuonEG_runF/191204_020625/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/MuonEG/crab_MuonEG_runF/191204_020625/0001/'], 'data','MuonEG','2017', 'F','1','1','1']

data2017_samples['2017_B_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleMuon/crab_SingleMuon_runB/191204_020201/0000/'], 'data','SingleMuon','2017', 'B','1','1','1']
data2017_samples['2017_C_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleMuon/crab_SingleMuon_runC/191204_020219/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleMuon/crab_SingleMuon_runC/191204_020219/0001/'], 'data','SingleMuon','2017', 'C','1','1','1']
data2017_samples['2017_D_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleMuon/crab_SingleMuon_runD/191204_020238/0000/'], 'data','SingleMuon','2017', 'D','1','1','1']
data2017_samples['2017_E_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleMuon/crab_SingleMuon_runE/191204_020257/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleMuon/crab_SingleMuon_runE/191204_020257/0001/'], 'data','SingleMuon','2017', 'E','1','1','1']
data2017_samples['2017_F_SingleMuon'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleMuon/crab_SingleMuon_runF/191204_020316/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleMuon/crab_SingleMuon_runF/191204_020316/0001/'], 'data','SingleMuon','2017', 'F','1','1','1']

data2017_samples['2017_B_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleElectron/crab_SingleElectron_runB/191204_020335/0000/'], 'data','SingleElectron','2017', 'B','1','1','1']
data2017_samples['2017_C_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleElectron/crab_SingleElectron_runC/191204_020353/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleElectron/crab_SingleElectron_runC/191204_020353/0001/'], 'data','SingleElectron','2017', 'C','1','1','1']
data2017_samples['2017_D_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleElectron/crab_SingleElectron_runD/191204_020411/0000/'], 'data','SingleElectron','2017', 'D','1','1','1']
data2017_samples['2017_E_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleElectron/crab_SingleElectron_runE/191204_020430/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleElectron/crab_SingleElectron_runE/191204_020430/0001/'], 'data','SingleElectron','2017', 'E','1','1','1']
data2017_samples['2017_F_SingleElectron'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleElectron/crab_SingleElectron_runF/191204_020448/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/SingleElectron/crab_SingleElectron_runF/191204_020448/0001/'], 'data','SingleElectron','2017', 'F','1','1','1']

data2017_samples['2017_B_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleEG/crab_DoubleEG_runB/200306_135425/0000/'], 'data','DoubleEG','2017', 'B','1','1','1']
data2017_samples['2017_C_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleEG/crab_DoubleEG_runC/200306_135450/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleEG/crab_DoubleEG_runC/200306_135450/0001/'], 'data','DoubleEG','2017', 'C','1','1','1']
data2017_samples['2017_D_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleEG/crab_DoubleEG_runD/200306_135510/0000/'], 'data','DoubleEG','2017', 'D','1','1','1']
data2017_samples['2017_E_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleEG/crab_DoubleEG_runE/200306_135529/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleEG/crab_DoubleEG_runE/200306_135529/0001/'], 'data','DoubleEG','2017', 'E','1','1','1']
data2017_samples['2017_F_DoubleEG'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleEG/crab_DoubleEG_runF/200306_135549/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleEG/crab_DoubleEG_runF/200306_135549/0001/'], 'data','DoubleEG','2017', 'F','1','1','1']

data2017_samples['2017_B_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleMuon/crab_DoubleMuon_runB/200306_135611/0000/'], 'data','DoubleMu','2017', 'B','1','1','1']
data2017_samples['2017_C_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleMuon/crab_DoubleMuon_runC/200306_135630/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleMuon/crab_DoubleMuon_runC/200306_135630/0001/'], 'data','DoubleMu','2017', 'C','1','1','1']
data2017_samples['2017_D_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleMuon/crab_DoubleMuon_runD/200306_135648/0000/'], 'data','DoubleMu','2017', 'D','1','1','1']
data2017_samples['2017_E_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleMuon/crab_DoubleMuon_runE/200306_135706/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleMuon/crab_DoubleMuon_runE/200306_135706/0001/'], 'data','DoubleMu','2017', 'E','1','1','1']
data2017_samples['2017_F_DoubleMu'] = [['/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleMuon/crab_DoubleMuon_runF/200306_135724/0000/','/pnfs/iihe/cms/store/user/xgao/samples-20191203/data/2017/DoubleMuon/crab_DoubleMuon_runF/200306_135724/0001/'], 'data','DoubleMu','2017', 'F','1','1','1']


