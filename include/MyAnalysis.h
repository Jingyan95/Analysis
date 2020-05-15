//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb 22 14:33:31 2020 by ROOT version 6.10/05
// from TTree IIHEAnalysis/IIHEAnalysis
// found on file: /pnfs/iihe/cms/store/user/schenara/SYS_2018/TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8/crab_TTTo2L2Nu_hdampUP_2018/200211_122927/0000/outfile_5_2016_71.root
//////////////////////////////////////////////////////////

#ifndef MyAnalysis_h
#define MyAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

using namespace std;
class MyAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       ev_event;
   ULong64_t       ev_run;
   ULong64_t       ev_luminosityBlock;
   UInt_t          ev_time;
   UInt_t          ev_time_unixTime;
   UInt_t          ev_time_microsecondOffset;
   Float_t         ev_fixedGridRhoAll;
   Float_t         ev_fixedGridRhoFastjetAll;
   Float_t         ev_fixedGridRhoFastjetAllCalo;
   Float_t         ev_fixedGridRhoFastjetCentralCalo;
   Float_t         ev_fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         ev_fixedGridRhoFastjetCentralNeutral;
   Float_t         ev_prefiringweight;
   Float_t         ev_prefiringweightup;
   Float_t         ev_prefiringweightdown;
   vector<float>   *LHE_Pt;
   vector<float>   *LHE_Eta;
   vector<float>   *LHE_Phi;
   vector<float>   *LHE_E;
   vector<int>     *LHE_pdgid;
   vector<int>     *LHE_status;
   Float_t         LHE_weight_nominal;
   vector<float>   *LHE_weight_sys;
   UInt_t          mc_n;
   Int_t           mc_nMEPartons;
   Int_t           mc_nMEPartonsFiltered;
   vector<float>   *mc_DJRValues;
   Float_t         mc_weight;
   Float_t         mc_w_sign;
   Int_t           mc_id_first;
   Int_t           mc_id_second;
   Float_t         mc_x_first;
   Float_t         mc_x_second;
   Float_t         mc_xPDF_first;
   Float_t         mc_xPDF_second;
   Float_t         mc_scalePDF;
   vector<int>     *mc_index;
   vector<int>     *mc_pdgId;
   vector<int>     *mc_charge;
   vector<int>     *mc_status;
   vector<int>     *mc_status_flags;
   vector<int>     *mc_status_tau_flags;
   vector<int>     *mc_tau_charge;
   vector<int>     *mc_tau_pdgId;
   vector<int>     *mc_tau_decay;
   vector<int>     *mc_tau_had_status;
   vector<int>     *mc_tau_had_charge;
   vector<int>     *mc_tau_had_pdgId;
   vector<float>   *mc_mass;
   vector<float>   *mc_px;
   vector<float>   *mc_py;
   vector<float>   *mc_pz;
   vector<float>   *mc_pt;
   vector<float>   *mc_eta;
   vector<float>   *mc_phi;
   vector<float>   *mc_energy;
   vector<float>   *mc_tau_pt;
   vector<float>   *mc_tau_eta;
   vector<float>   *mc_tau_phi;
   vector<float>   *mc_tau_energy;
   vector<float>   *mc_tau_had_pt;
   vector<float>   *mc_tau_had_eta;
   vector<float>   *mc_tau_had_phi;
   vector<float>   *mc_tau_had_energy;
   vector<unsigned int> *mc_numberOfDaughters;
   vector<unsigned int> *mc_numberOfMothers;
   vector<vector<int> > *mc_mother_index;
   vector<vector<int> > *mc_mother_pdgId;
   vector<vector<float> > *mc_mother_px;
   vector<vector<float> > *mc_mother_py;
   vector<vector<float> > *mc_mother_pz;
   vector<vector<float> > *mc_mother_pt;
   vector<vector<float> > *mc_mother_eta;
   vector<vector<float> > *mc_mother_phi;
   vector<vector<float> > *mc_mother_energy;
   vector<vector<float> > *mc_mother_mass;
   Int_t           mc_trueNumInteractions;
   Int_t           mc_PU_NumInteractions;
   vector<float>   *genjet_pt;
   vector<float>   *genjet_eta;
   vector<float>   *genjet_phi;
   vector<float>   *genjet_energy;
   vector<float>   *genjetAK8_pt;
   vector<float>   *genjetAK8_eta;
   vector<float>   *genjetAK8_phi;
   vector<float>   *genjetAK8_energy;
   vector<float>   *gen_weight_sys;
   UInt_t          pv_n;
   vector<float>   *pv_x;
   vector<float>   *pv_y;
   vector<float>   *pv_z;
   vector<float>   *pv_ndof;
   vector<float>   *pv_normalizedChi2;
   vector<int>     *pv_isValid;
   vector<int>     *pv_isFake;
   UInt_t          gsf_n;
   vector<int>     *gsf_classification;
   vector<float>   *gsf_ecalTrkEnergyPostCorr;
   vector<float>   *gsf_ecalEnergyPostCorr;
   vector<float>   *gsf_energy;
   vector<float>   *gsf_p;
   vector<float>   *gsf_pt;
   vector<float>   *gsf_et;
   vector<float>   *gsf_scE1x5;
   vector<float>   *gsf_scE5x5;
   vector<float>   *gsf_scE2x5Max;
   vector<float>   *gsf_full5x5_e5x5;
   vector<float>   *gsf_full5x5_e1x5;
   vector<float>   *gsf_full5x5_e2x5Max;
   vector<float>   *gsf_full5x5_sigmaIetaIeta;
   vector<float>   *gsf_full5x5_hcalOverEcal;
   vector<float>   *gsf_eta;
   vector<float>   *gsf_phi;
   vector<float>   *gsf_theta;
   vector<float>   *gsf_px;
   vector<float>   *gsf_py;
   vector<float>   *gsf_pz;
   vector<float>   *gsf_caloEnergy;
   vector<float>   *gsf_deltaEtaSuperClusterTrackAtVtx;
   vector<float>   *gsf_deltaPhiSuperClusterTrackAtVtx;
   vector<float>   *gsf_hadronicOverEm;
   vector<float>   *gsf_hcalDepth1OverEcal;
   vector<float>   *gsf_hcalDepth2OverEcal;
   vector<float>   *gsf_dr03TkSumPt;
   vector<float>   *gsf_heepTrkPtIso;
   vector<float>   *gsf_dr03EcalRecHitSumEt;
   vector<float>   *gsf_dr03HcalDepth1TowerSumEt;
   vector<float>   *gsf_dr03HcalDepth2TowerSumEt;
   vector<float>   *gsf_sigmaIetaIeta;
   vector<int>     *gsf_charge;
   vector<int>     *gsf_ecaldrivenSeed;
   vector<int>     *gsf_trackerdrivenSeed;
   vector<int>     *gsf_isEB;
   vector<int>     *gsf_isEE;
   vector<int>     *gsf_isPF;
   vector<int>     *gsf_passConversionVeto;
   vector<int>     *gsf_VID_cutBasedElectronID_Fall17_94X_V2_loose;
   vector<int>     *gsf_VID_cutBasedElectronID_Fall17_94X_V2_medium;
   vector<int>     *gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight;
   vector<int>     *gsf_VID_cutBasedElectronID_Fall17_94X_V2_veto;
   vector<int>     *gsf_VID_heepElectronID_HEEPV70;
   vector<int>     *gsf_VID_mvaEleID_Fall17_iso_V2_wp80;
   vector<int>     *gsf_VID_mvaEleID_Fall17_iso_V2_wp90;
   vector<int>     *gsf_VID_mvaEleID_Fall17_iso_V2_wpLoose;
   vector<int>     *gsf_VID_mvaEleID_Fall17_noIso_V2_wp80;
   vector<int>     *gsf_VID_mvaEleID_Fall17_noIso_V2_wp90;
   vector<int>     *gsf_VID_mvaEleID_Fall17_noIso_V2_wpLoose;
   vector<float>   *gsf_deltaEtaSeedClusterTrackAtCalo;
   vector<float>   *gsf_deltaPhiSeedClusterTrackAtCalo;
   vector<float>   *gsf_ecalEnergy;
   vector<float>   *gsf_eSuperClusterOverP;
   vector<float>   *gsf_dxy;
   vector<float>   *gsf_dxy_beamSpot;
   vector<float>   *gsf_dxy_firstPVtx;
   vector<float>   *gsf_dxyError;
   vector<float>   *gsf_dz;
   vector<float>   *gsf_dz_beamSpot;
   vector<float>   *gsf_dz_firstPVtx;
   vector<float>   *gsf_dzError;
   vector<float>   *gsf_vz;
   vector<int>     *gsf_numberOfValidHits;
   vector<int>     *gsf_nLostInnerHits;
   vector<int>     *gsf_nLostOuterHits;
   vector<int>     *gsf_convFlags;
   vector<float>   *gsf_convDist;
   vector<float>   *gsf_convDcot;
   vector<float>   *gsf_convRadius;
   vector<float>   *gsf_fBrem;
   vector<float>   *gsf_e1x5;
   vector<float>   *gsf_e2x5Max;
   vector<float>   *gsf_e5x5;
   vector<float>   *gsf_r9;
   vector<float>   *gsf_deltaEtaSeedClusterTrackAtVtx;
   vector<float>   *gsf_effArea;
   vector<float>   *gsf_sumChargedHadronPt;
   vector<float>   *gsf_sumNeutralHadronEt;
   vector<float>   *gsf_sumPhotonEt;
   vector<float>   *gsf_sumPUPt;
   vector<float>   *gsf_ooEmooP;
   vector<vector<int> > *gsf_hitsinfo;
   vector<float>   *gsf_pixelMatch_dPhi1;
   vector<float>   *gsf_pixelMatch_dPhi2;
   vector<float>   *gsf_pixelMatch_dRz1;
   vector<float>   *gsf_pixelMatch_dRz2;
   vector<int>     *gsf_pixelMatch_subDetector1;
   vector<int>     *gsf_pixelMatch_subDetector2;
   vector<float>   *gsf_mc_bestDR;
   vector<int>     *gsf_mc_index;
   vector<float>   *gsf_mc_ERatio;
   vector<float>   *gsf_sc_energy;
   vector<float>   *gsf_sc_seed_eta;
   vector<float>   *gsf_sc_eta;
   vector<float>   *gsf_sc_etacorr;
   vector<float>   *gsf_sc_theta;
   vector<float>   *gsf_sc_thetacorr;
   vector<float>   *gsf_sc_et;
   vector<float>   *gsf_sc_phi;
   vector<float>   *gsf_sc_px;
   vector<float>   *gsf_sc_py;
   vector<float>   *gsf_sc_pz;
   vector<float>   *gsf_sc_x;
   vector<float>   *gsf_sc_y;
   vector<float>   *gsf_sc_z;
   vector<float>   *gsf_sc_phiWidth;
   vector<float>   *gsf_sc_etaWidth;
   vector<int>     *gsf_sc_seed_rawId;
   vector<int>     *gsf_sc_seed_ieta;
   vector<int>     *gsf_sc_seed_iphi;
   vector<int>     *gsf_sc_seed_kHasSwitchToGain6;
   vector<int>     *gsf_sc_seed_kHasSwitchToGain1;
   Int_t           EHits_isSaturated;
   vector<int>     *EBHits_rawId;
   vector<int>     *EBHits_iRechit;
   vector<int>     *EBHits_ieta;
   vector<int>     *EBHits_iphi;
   vector<int>     *EBHits_RecoFlag;
   vector<int>     *EBHits_kSaturated;
   vector<int>     *EBHits_kLeadingEdgeRecovered;
   vector<int>     *EBHits_kNeighboursRecovered;
   vector<int>     *EBHits_kWeird;
   vector<float>   *EBHits_energy;
   vector<int>     *EEHits_rawId;
   vector<int>     *EEHits_iRechit;
   vector<int>     *EEHits_ieta;
   vector<int>     *EEHits_iphi;
   vector<int>     *EEHits_RecoFlag;
   vector<int>     *EEHits_kSaturated;
   vector<int>     *EEHits_kLeadingEdgeRecovered;
   vector<int>     *EEHits_kNeighboursRecovered;
   vector<int>     *EEHits_kWeird;
   vector<float>   *EEHits_energy;
   UInt_t          mu_n;
   vector<float>   *mu_gt_qoverp;
   vector<int>     *mu_gt_charge;
   vector<float>   *mu_gt_pt;
   vector<float>   *mu_gt_eta;
   vector<float>   *mu_gt_phi;
   vector<float>   *mu_gt_p;
   vector<float>   *mu_gt_px;
   vector<float>   *mu_gt_py;
   vector<float>   *mu_gt_pz;
   vector<float>   *mu_gt_theta;
   vector<float>   *mu_gt_lambda;
   vector<float>   *mu_gt_d0;
   vector<float>   *mu_gt_dz;
   vector<float>   *mu_gt_dz_beamspot;
   vector<float>   *mu_gt_dz_firstPVtx;
   vector<float>   *mu_gt_dxy;
   vector<float>   *mu_gt_dxy_beamspot;
   vector<float>   *mu_gt_dxy_firstPVtx;
   vector<float>   *mu_gt_dsz;
   vector<float>   *mu_gt_vx;
   vector<float>   *mu_gt_vy;
   vector<float>   *mu_gt_vz;
   vector<float>   *mu_gt_qoverpError;
   vector<float>   *mu_gt_ptError;
   vector<float>   *mu_gt_thetaError;
   vector<float>   *mu_gt_lambdaError;
   vector<float>   *mu_gt_phiError;
   vector<float>   *mu_gt_dxyError;
   vector<float>   *mu_gt_d0Error;
   vector<float>   *mu_gt_dszError;
   vector<float>   *mu_gt_dzError;
   vector<float>   *mu_gt_etaError;
   vector<float>   *mu_gt_chi2;
   vector<float>   *mu_gt_ndof;
   vector<float>   *mu_gt_normalizedChi2;
   vector<float>   *mu_ot_qoverp;
   vector<int>     *mu_ot_charge;
   vector<float>   *mu_ot_pt;
   vector<float>   *mu_ot_eta;
   vector<float>   *mu_ot_phi;
   vector<float>   *mu_ot_p;
   vector<float>   *mu_ot_px;
   vector<float>   *mu_ot_py;
   vector<float>   *mu_ot_pz;
   vector<float>   *mu_ot_theta;
   vector<float>   *mu_ot_lambda;
   vector<float>   *mu_ot_d0;
   vector<float>   *mu_ot_dz;
   vector<float>   *mu_ot_dz_beamspot;
   vector<float>   *mu_ot_dz_firstPVtx;
   vector<float>   *mu_ot_dxy;
   vector<float>   *mu_ot_dxy_beamspot;
   vector<float>   *mu_ot_dxy_firstPVtx;
   vector<float>   *mu_ot_dsz;
   vector<float>   *mu_ot_vx;
   vector<float>   *mu_ot_vy;
   vector<float>   *mu_ot_vz;
   vector<float>   *mu_ot_qoverpError;
   vector<float>   *mu_ot_ptError;
   vector<float>   *mu_ot_thetaError;
   vector<float>   *mu_ot_lambdaError;
   vector<float>   *mu_ot_phiError;
   vector<float>   *mu_ot_dxyError;
   vector<float>   *mu_ot_d0Error;
   vector<float>   *mu_ot_dszError;
   vector<float>   *mu_ot_dzError;
   vector<float>   *mu_ot_etaError;
   vector<float>   *mu_ot_chi2;
   vector<float>   *mu_ot_ndof;
   vector<float>   *mu_ot_normalizedChi2;
   vector<float>   *mu_it_qoverp;
   vector<int>     *mu_it_charge;
   vector<float>   *mu_it_pt;
   vector<float>   *mu_it_eta;
   vector<float>   *mu_it_phi;
   vector<float>   *mu_it_p;
   vector<float>   *mu_it_px;
   vector<float>   *mu_it_py;
   vector<float>   *mu_it_pz;
   vector<float>   *mu_it_theta;
   vector<float>   *mu_it_lambda;
   vector<float>   *mu_it_d0;
   vector<float>   *mu_it_dz;
   vector<float>   *mu_it_dz_beamspot;
   vector<float>   *mu_it_dz_firstPVtx;
   vector<float>   *mu_it_dxy;
   vector<float>   *mu_it_dxy_beamspot;
   vector<float>   *mu_it_dxy_firstPVtx;
   vector<float>   *mu_it_dsz;
   vector<float>   *mu_it_vx;
   vector<float>   *mu_it_vy;
   vector<float>   *mu_it_vz;
   vector<float>   *mu_it_qoverpError;
   vector<float>   *mu_it_ptError;
   vector<float>   *mu_it_thetaError;
   vector<float>   *mu_it_lambdaError;
   vector<float>   *mu_it_phiError;
   vector<float>   *mu_it_dxyError;
   vector<float>   *mu_it_d0Error;
   vector<float>   *mu_it_dszError;
   vector<float>   *mu_it_dzError;
   vector<float>   *mu_it_etaError;
   vector<float>   *mu_it_chi2;
   vector<float>   *mu_it_ndof;
   vector<float>   *mu_it_normalizedChi2;
   vector<float>   *mu_ibt_qoverp;
   vector<int>     *mu_ibt_charge;
   vector<float>   *mu_ibt_pt;
   vector<float>   *mu_ibt_eta;
   vector<float>   *mu_ibt_phi;
   vector<float>   *mu_ibt_p;
   vector<float>   *mu_ibt_px;
   vector<float>   *mu_ibt_py;
   vector<float>   *mu_ibt_pz;
   vector<float>   *mu_ibt_theta;
   vector<float>   *mu_ibt_lambda;
   vector<float>   *mu_ibt_d0;
   vector<float>   *mu_ibt_dz;
   vector<float>   *mu_ibt_dz_beamspot;
   vector<float>   *mu_ibt_dz_firstPVtx;
   vector<float>   *mu_ibt_dxy;
   vector<float>   *mu_ibt_dxy_beamspot;
   vector<float>   *mu_ibt_dxy_firstPVtx;
   vector<float>   *mu_ibt_dsz;
   vector<float>   *mu_ibt_vx;
   vector<float>   *mu_ibt_vy;
   vector<float>   *mu_ibt_vz;
   vector<float>   *mu_ibt_qoverpError;
   vector<float>   *mu_ibt_ptError;
   vector<float>   *mu_ibt_thetaError;
   vector<float>   *mu_ibt_lambdaError;
   vector<float>   *mu_ibt_phiError;
   vector<float>   *mu_ibt_dxyError;
   vector<float>   *mu_ibt_d0Error;
   vector<float>   *mu_ibt_dszError;
   vector<float>   *mu_ibt_dzError;
   vector<float>   *mu_ibt_etaError;
   vector<float>   *mu_ibt_chi2;
   vector<float>   *mu_ibt_ndof;
   vector<float>   *mu_ibt_normalizedChi2;
   vector<int>     *mu_isGlobalMuon;
   vector<int>     *mu_isStandAloneMuon;
   vector<int>     *mu_isTrackerMuon;
   vector<int>     *mu_isPFMuon;
   vector<int>     *mu_isPFIsolationValid;
   vector<int>     *mu_isGoodMuonTMLastStationLoose;
   vector<int>     *mu_isGoodMuonTMLastStationTight;
   vector<int>     *mu_isGoodMuonTM2DCompatibilityLoose;
   vector<int>     *mu_isGoodMuonTM2DCompatibilityTight;
   vector<int>     *mu_isGoodMuonTMOneStationLoose;
   vector<int>     *mu_isGoodMuonTMOneStationTight;
   vector<int>     *mu_isGoodMuonTMLastStationOptimizedLowPtLoose;
   vector<int>     *mu_isGoodMuonTMLastStationOptimizedLowPtTight;
   vector<int>     *mu_isTightMuon;
   vector<int>     *mu_isMediumMuon;
   vector<int>     *mu_isLooseMuon;
   vector<int>     *mu_isSoftMuon;
   vector<int>     *mu_isHighPtMuonOld;
   vector<int>     *mu_isHighPtMuon;
   vector<int>     *mu_isTrackerHighPtMuon;
   vector<int>     *mu_CutBasedIdLoose;
   vector<int>     *mu_CutBasedIdMedium;
   vector<int>     *mu_CutBasedIdMediumPrompt;
   vector<int>     *mu_CutBasedIdTight;
   vector<int>     *mu_CutBasedIdGlobalHighPt;
   vector<int>     *mu_CutBasedIdTrkHighPt;
   vector<int>     *mu_PFIsoVeryLoose;
   vector<int>     *mu_PFIsoLoose;
   vector<int>     *mu_PFIsoMedium;
   vector<int>     *mu_PFIsoTight;
   vector<int>     *mu_PFIsoVeryTight;
   vector<int>     *mu_PFIsoVeryVeryTight;
   vector<int>     *mu_TkIsoLoose;
   vector<int>     *mu_TkIsoTight;
   vector<int>     *mu_SoftCutBasedId;
   vector<int>     *mu_SoftMvaId;
   vector<int>     *mu_MvaLoose;
   vector<int>     *mu_MvaMedium;
   vector<int>     *mu_MvaTight;
   vector<int>     *mu_MiniIsoLoose;
   vector<int>     *mu_MiniIsoMedium;
   vector<int>     *mu_MiniIsoTight;
   vector<int>     *mu_MiniIsoVeryTight;
   vector<int>     *mu_TriggerIdLoose;
   vector<int>     *mu_InTimeMuon;
   vector<int>     *mu_MultiIsoLoose;
   vector<int>     *mu_MultiIsoMedium;
   vector<int>     *mu_numberOfMatchedStations;
   vector<int>     *mu_numberOfValidPixelHits;
   vector<int>     *mu_trackerLayersWithMeasurement;
   vector<int>     *mu_numberOfValidMuonHits;
   vector<int>     *mu_pixelLayersWithMeasurement;
   vector<float>   *mu_innerTrack_validFraction;
   vector<float>   *mu_combinedQuality_trkKink;
   vector<float>   *mu_combinedQuality_chi2LocalPosition;
   vector<float>   *mu_segmentCompatibility;
   vector<float>   *mu_dB;
   vector<float>   *mu_pt_default;
   vector<float>   *mu_isolationR03_sumPt;
   vector<float>   *mu_isolationR03_trackerVetoPt;
   vector<float>   *mu_isolationR03_emEt;
   vector<float>   *mu_isolationR03_emVetoEt;
   vector<float>   *mu_isolationR03_hadEt;
   vector<float>   *mu_isolationR03_hadVetoEt;
   vector<float>   *mu_isolationR05_sumPt;
   vector<float>   *mu_isolationR05_trackerVetoPt;
   vector<float>   *mu_isolationR05_emEt;
   vector<float>   *mu_isolationR05_emVetoEt;
   vector<float>   *mu_isolationR05_hadEt;
   vector<float>   *mu_isolationR05_hadVetoEt;
   vector<float>   *mu_pfIsolationR03_sumChargedHadronPt;
   vector<float>   *mu_pfIsolationR03_sumNeutralHadronEt;
   vector<float>   *mu_pfIsolationR03_sumChargedParticlePt;
   vector<float>   *mu_pfIsolationR03_sumPhotonEt;
   vector<float>   *mu_pfIsolationR03_sumNeutralHadronEtHighThreshold;
   vector<float>   *mu_pfIsolationR03_sumPhotonEtHighThreshold;
   vector<float>   *mu_pfIsolationR03_sumPUPt;
   vector<float>   *mu_pfIsolationR04_sumChargedHadronPt;
   vector<float>   *mu_pfIsolationR04_sumNeutralHadronEt;
   vector<float>   *mu_pfIsolationR04_sumChargedParticlePt;
   vector<float>   *mu_pfIsolationR04_sumPhotonEt;
   vector<float>   *mu_pfIsolationR04_sumNeutralHadronEtHighThreshold;
   vector<float>   *mu_pfIsolationR04_sumPhotonEtHighThreshold;
   vector<float>   *mu_pfIsolationR04_sumPUPt;
   vector<float>   *mu_pfIsoDbCorrected03;
   vector<float>   *mu_pfIsoDbCorrected04;
   vector<float>   *mu_isoTrackerBased03;
   vector<float>   *mu_mc_bestDR;
   vector<int>     *mu_mc_index;
   vector<float>   *mu_mc_ERatio;
   UInt_t          jet_n;
   vector<float>   *jet_px;
   vector<float>   *jet_py;
   vector<float>   *jet_pz;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_theta;
   vector<float>   *jet_phi;
   vector<float>   *jet_energy;
   vector<float>   *jet_mass;
   vector<float>   *jet_chargedEmEnergyFraction;
   vector<float>   *jet_neutralHadronEnergyFraction;
   vector<float>   *jet_neutralEmEnergyFraction;
   vector<float>   *jet_chargedHadronEnergyFraction;
   vector<float>   *jet_muonEnergyFraction;
   vector<int>     *jet_chargedMultiplicity;
   vector<int>     *jet_neutralMultiplicity;
   vector<int>     *jet_partonFlavour;
   vector<int>     *jet_hadronFlavour;
   vector<float>   *jet_CSVv2;
   vector<float>   *jet_CvsL;
   vector<float>   *jet_CvsB;
   vector<float>   *jet_MVA2BJets;
   vector<float>   *jet_CvsB_DeepJet_charm_tagger;
   vector<float>   *jet_CvsL_DeepJet_charm_tagger;
   vector<float>   *jet_CvsB_DeepCSV_charm_tagger;
   vector<float>   *jet_CvsL_DeepCSV_charm_tagger;
   vector<float>   *jet_DeepJet;
   vector<float>   *jet_DeepCSV;
   vector<int>     *jet_isJetIDLoose_2016;
   vector<int>     *jet_isJetIDTight_2016;
   vector<int>     *jet_isJetIDTightLepVeto_2016;
   vector<int>     *jet_isJetID_2017;
   vector<int>     *jet_isJetIDLepVeto_2017;
   vector<int>     *jet_isJetID_2018;
   vector<int>     *jet_isJetIDLepVeto_2018;
   vector<float>   *jet_Smeared_pt;
   vector<float>   *jet_SmearedJetResUp_pt;
   vector<float>   *jet_SmearedJetResDown_pt;
   vector<float>   *jet_SmearedJetEnUp_pt;
   vector<float>   *jet_SmearedJetEnDown_pt;
   Float_t         BtagSF_Deepcsv_loose;
   Float_t         BtagSF_Deepcsv_Up_loose;
   Float_t         BtagSF_Deepcsv_Down_loose;
   Float_t         BtagSF_Deepcsv_medium;
   Float_t         BtagSF_Deepcsv_Up_medium;
   Float_t         BtagSF_Deepcsv_Down_medium;
   Float_t         BtagSF_Deepcsv_tight;
   Float_t         BtagSF_Deepcsv_Up_tight;
   Float_t         BtagSF_Deepcsv_Down_tight;
   Float_t         MET_nominal_Pt;
   Float_t         MET_nominal_phi;
   Float_t         MET_nominal_significance;
   Float_t         MET_Pt;
   Float_t         MET_phi;
   Float_t         MET_significance;
   Float_t         MET_T1_Pt;
   Float_t         MET_T1_phi;
   Float_t         MET_T1_significance;
   Float_t         MET_gen_pt;
   Float_t         MET_gen_phi;
   vector<float>   *MET_Type1Unc_Px;
   vector<float>   *MET_Type1Unc_Py;
   vector<float>   *MET_Type1Unc_Pt;
   Float_t         MET_T1JetEnDown_Pt;
   Float_t         MET_T1JetEnDown_phi;
   Float_t         MET_T1JetEnDown_significance;
   Float_t         MET_T1JetEnUp_Pt;
   Float_t         MET_T1JetEnUp_phi;
   Float_t         MET_T1JetEnUp_significance;
   Float_t         MET_T1Smear_Pt;
   Float_t         MET_T1Smear_phi;
   Float_t         MET_T1Smear_significance;
   Float_t         MET_T1SmearJetEnDown_Pt;
   Float_t         MET_T1SmearJetEnDown_phi;
   Float_t         MET_T1SmearJetEnDown_significance;
   Float_t         MET_T1SmearJetEnUp_Pt;
   Float_t         MET_T1SmearJetEnUp_phi;
   Float_t         MET_T1SmearJetEnUp_significance;
   Float_t         MET_T1SmearJetResDown_Pt;
   Float_t         MET_T1SmearJetResDown_phi;
   Float_t         MET_T1SmearJetResDown_significance;
   Float_t         MET_T1SmearJetResUp_Pt;
   Float_t         MET_T1SmearJetResUp_phi;
   Float_t         MET_T1SmearJetResUp_significance;
   Float_t         MET_T1Txy_Pt;
   Float_t         MET_T1Txy_phi;
   Float_t         MET_T1Txy_significance;
   Float_t         MET_FinalCollection_Pt;
   Float_t         MET_FinalCollection_phi;
   Float_t         MET_FinalCollection_significance;
   Int_t           trig_Flag_ecalBadCalibReduced;
   Int_t           trig_Flag_HBHENoiseFilter_accept;
   Int_t           trig_Flag_HBHENoiseIsoFilter_accept;
   Int_t           trig_Flag_CSCTightHaloFilter_accept;
   Int_t           trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept;
   Int_t           trig_Flag_CSCTightHalo2015Filter_accept;
   Int_t           trig_Flag_globalTightHalo2016Filter_accept;
   Int_t           trig_Flag_globalSuperTightHalo2016Filter_accept;
   Int_t           trig_Flag_HcalStripHaloFilter_accept;
   Int_t           trig_Flag_hcalLaserEventFilter_accept;
   Int_t           trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept;
   Int_t           trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept;
   Int_t           trig_Flag_ecalBadCalibFilter_accept;
   Int_t           trig_Flag_goodVertices_accept;
   Int_t           trig_Flag_eeBadScFilter_accept;
   Int_t           trig_Flag_ecalLaserCorrFilter_accept;
   Int_t           trig_Flag_trkPOGFilters_accept;
   Int_t           trig_Flag_chargedHadronTrackResolutionFilter_accept;
   Int_t           trig_Flag_muonBadTrackFilter_accept;
   Int_t           trig_Flag_BadChargedCandidateFilter_accept;
   Int_t           trig_Flag_BadPFMuonFilter_accept;
   Int_t           trig_Flag_BadChargedCandidateSummer16Filter_accept;
   Int_t           trig_Flag_BadPFMuonSummer16Filter_accept;
   Int_t           trig_Flag_trkPOG_manystripclus53X_accept;
   Int_t           trig_Flag_trkPOG_toomanystripclus53X_accept;
   Int_t           trig_Flag_trkPOG_logErrorTooManyClusters_accept;
   Int_t           trig_Flag_METFilters_accept;
   Int_t           trig_raw2digi_step_accept;
   Int_t           trig_L1Reco_step_accept;
   Int_t           trig_reconstruction_step_accept;
   Int_t           trig_recosim_step_accept;
   Int_t           trig_eventinterpretaion_step_accept;
   Int_t           trig_HLT_DoubleEle25_CaloIdL_MW_accept;
   Int_t           trig_HLT_DoubleEle27_CaloIdL_MW_accept;
   Int_t           trig_HLT_DoubleEle33_CaloIdL_MW_accept;
   Int_t           trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele27_Ele37_CaloIdL_MW_accept;
   Int_t           trig_HLT_Mu27_Ele37_CaloIdL_MW_accept;
   Int_t           trig_HLT_Mu37_Ele27_CaloIdL_MW_accept;
   Int_t           trig_HLT_Mu37_TkMu27_accept;
   Int_t           trig_HLT_DoubleMu4_3_Bs_accept;
   Int_t           trig_HLT_Ele20_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele15_WPLoose_Gsf_accept;
   Int_t           trig_HLT_Ele17_WPLoose_Gsf_accept;
   Int_t           trig_HLT_Ele20_WPLoose_Gsf_accept;
   Int_t           trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept;
   Int_t           trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept;
   Int_t           trig_HLT_Ele27_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele28_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele30_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele32_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele35_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept;
   Int_t           trig_HLT_Ele38_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele40_WPTight_Gsf_accept;
   Int_t           trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept;
   Int_t           trig_HLT_IsoMu20_accept;
   Int_t           trig_HLT_IsoMu24_accept;
   Int_t           trig_HLT_IsoTkMu24_accept;
   Int_t           trig_HLT_IsoMu24_eta2p1_accept;
   Int_t           trig_HLT_IsoMu27_accept;
   Int_t           trig_HLT_IsoMu30_accept;
   Int_t           trig_HLT_L1SingleMu18_accept;
   Int_t           trig_HLT_L1SingleMu25_accept;
   Int_t           trig_HLT_L2Mu10_accept;
   Int_t           trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept;
   Int_t           trig_HLT_L2Mu10_NoVertex_NoBPTX_accept;
   Int_t           trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept;
   Int_t           trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept;
   Int_t           trig_HLT_L2Mu50_accept;
   Int_t           trig_HLT_L2Mu23NoVtx_2Cha_accept;
   Int_t           trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept;
   Int_t           trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept;
   Int_t           trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept;
   Int_t           trig_HLT_DoubleL2Mu50_accept;
   Int_t           trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept;
   Int_t           trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept;
   Int_t           trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept;
   Int_t           trig_HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_accept;
   Int_t           trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept;
   Int_t           trig_HLT_Mu25_TkMu0_Onia_accept;
   Int_t           trig_HLT_Mu30_TkMu0_Psi_accept;
   Int_t           trig_HLT_Mu20_TkMu0_Phi_accept;
   Int_t           trig_HLT_Mu25_TkMu0_Phi_accept;
   Int_t           trig_HLT_Mu12_accept;
   Int_t           trig_HLT_Mu15_accept;
   Int_t           trig_HLT_Mu20_accept;
   Int_t           trig_HLT_Mu27_accept;
   Int_t           trig_HLT_Mu50_accept;
   Int_t           trig_HLT_Mu55_accept;
   Int_t           trig_HLT_OldMu100_accept;
   Int_t           trig_HLT_TkMu100_accept;
   Int_t           trig_HLT_PFHT500_PFMET100_PFMHT100_IDTight_accept;
   Int_t           trig_HLT_PFHT500_PFMET110_PFMHT110_IDTight_accept;
   Int_t           trig_HLT_PFHT700_PFMET85_PFMHT85_IDTight_accept;
   Int_t           trig_HLT_PFHT700_PFMET95_PFMHT95_IDTight_accept;
   Int_t           trig_HLT_PFHT800_PFMET75_PFMHT75_IDTight_accept;
   Int_t           trig_HLT_PFHT800_PFMET85_PFMHT85_IDTight_accept;
   Int_t           trig_HLT_PFMET110_PFMHT110_IDTight_accept;
   Int_t           trig_HLT_PFMET120_PFMHT120_IDTight_accept;
   Int_t           trig_HLT_PFMET130_PFMHT130_IDTight_accept;
   Int_t           trig_HLT_PFMET140_PFMHT140_IDTight_accept;
   Int_t           trig_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_accept;
   Int_t           trig_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_accept;
   Int_t           trig_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_accept;
   Int_t           trig_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_accept;
   Int_t           trig_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_accept;
   Int_t           trig_HLT_PFMET120_PFMHT120_IDTight_PFHT60_accept;
   Int_t           trig_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_accept;
   Int_t           trig_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_accept;
   Int_t           trig_HLT_PFMETTypeOne110_PFMHT110_IDTight_accept;
   Int_t           trig_HLT_PFMETTypeOne120_PFMHT120_IDTight_accept;
   Int_t           trig_HLT_PFMETTypeOne130_PFMHT130_IDTight_accept;
   Int_t           trig_HLT_PFMETTypeOne140_PFMHT140_IDTight_accept;
   Int_t           trig_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_accept;
   Int_t           trig_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_accept;
   Int_t           trig_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_accept;
   Int_t           trig_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_accept;
   Int_t           trig_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_accept;
   Int_t           trig_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_accept;
   Int_t           trig_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_accept;
   Int_t           trig_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_accept;
   Int_t           trig_HLT_CaloMET80_NotCleaned_accept;
   Int_t           trig_HLT_CaloMET90_NotCleaned_accept;
   Int_t           trig_HLT_CaloMET100_NotCleaned_accept;
   Int_t           trig_HLT_CaloMET110_NotCleaned_accept;
   Int_t           trig_HLT_CaloMET250_NotCleaned_accept;
   Int_t           trig_HLT_CaloMET70_HBHECleaned_accept;
   Int_t           trig_HLT_CaloMET80_HBHECleaned_accept;
   Int_t           trig_HLT_CaloMET90_HBHECleaned_accept;
   Int_t           trig_HLT_CaloMET100_HBHECleaned_accept;
   Int_t           trig_HLT_CaloMET250_HBHECleaned_accept;
   Int_t           trig_HLT_CaloMET300_HBHECleaned_accept;
   Int_t           trig_HLT_CaloMET350_HBHECleaned_accept;
   Int_t           trig_HLT_PFMET200_NotCleaned_accept;
   Int_t           trig_HLT_PFMET200_HBHECleaned_accept;
   Int_t           trig_HLT_PFMET250_HBHECleaned_accept;
   Int_t           trig_HLT_PFMET300_HBHECleaned_accept;
   Int_t           trig_HLT_PFMET200_HBHE_BeamHaloCleaned_accept;
   Int_t           trig_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_accept;
   Int_t           trig_HLT_MET105_IsoTrk50_accept;
   Int_t           trig_HLT_MET120_IsoTrk50_accept;
   Int_t           trig_HLT_Photon300_NoHE_accept;
   Int_t           trig_HLT_Mu8_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu17_TrkIsoVVL_accept;
   Int_t           trig_HLT_Mu19_TrkIsoVVL_accept;
   Int_t           trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;
   Int_t           trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_HLT_Photon20_accept;
   Int_t           trig_HLT_Photon33_accept;
   Int_t           trig_HLT_Photon50_accept;
   Int_t           trig_HLT_Photon75_accept;
   Int_t           trig_HLT_Photon90_accept;
   Int_t           trig_HLT_Photon120_accept;
   Int_t           trig_HLT_Photon150_accept;
   Int_t           trig_HLT_Photon175_accept;
   Int_t           trig_HLT_Photon200_accept;
   Int_t           trig_HLT_Photon100EB_TightID_TightIso_accept;
   Int_t           trig_HLT_Photon110EB_TightID_TightIso_accept;
   Int_t           trig_HLT_Photon120EB_TightID_TightIso_accept;
   Int_t           trig_HLT_Photon100EBHE10_accept;
   Int_t           trig_HLT_Photon100EEHE10_accept;
   Int_t           trig_HLT_Photon100EE_TightID_TightIso_accept;
   Int_t           trig_HLT_Photon50_R9Id90_HE10_IsoM_accept;
   Int_t           trig_HLT_Photon75_R9Id90_HE10_IsoM_accept;
   Int_t           trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_accept;
   Int_t           trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_accept;
   Int_t           trig_HLT_Photon90_R9Id90_HE10_IsoM_accept;
   Int_t           trig_HLT_Photon120_R9Id90_HE10_IsoM_accept;
   Int_t           trig_HLT_Photon165_R9Id90_HE10_IsoM_accept;
   Int_t           trig_HLT_Photon90_CaloIdL_PFHT700_accept;
   Int_t           trig_HLT_Photon35_TwoProngs35_accept;
   Int_t           trig_HLT_IsoMu24_TwoProngs35_accept;
   Int_t           trig_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_accept;
   Int_t           trig_HLT_DoubleMu3_DZ_PFMET70_PFMHT70_accept;
   Int_t           trig_HLT_DoubleMu3_DZ_PFMET90_PFMHT90_accept;
   Int_t           trig_HLT_DoubleMu43NoFiltersNoVtx_accept;
   Int_t           trig_HLT_DoubleMu48NoFiltersNoVtx_accept;
   Int_t           trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept;
   Int_t           trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept;
   Int_t           trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_accept;
   Int_t           trig_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_accept;
   Int_t           trig_HLT_DoubleMu33NoFiltersNoVtxDisplaced_accept;
   Int_t           trig_HLT_DoubleMu40NoFiltersNoVtxDisplaced_accept;
   Int_t           trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept;
   Int_t           trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept;
   Int_t           trig_HLT_DiJet110_35_Mjj650_PFMET110_accept;
   Int_t           trig_HLT_DiJet110_35_Mjj650_PFMET120_accept;
   Int_t           trig_HLT_DiJet110_35_Mjj650_PFMET130_accept;
   Int_t           trig_HLT_TripleJet110_35_35_Mjj650_PFMET110_accept;
   Int_t           trig_HLT_TripleJet110_35_35_Mjj650_PFMET120_accept;
   Int_t           trig_HLT_TripleJet110_35_35_Mjj650_PFMET130_accept;
   Int_t           trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept;
   Int_t           trig_HLT_Ele15_IsoVVVL_PFHT450_PFMET50_accept;
   Int_t           trig_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_accept;
   Int_t           trig_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_accept;
   Int_t           trig_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_accept;
   Int_t           trig_HLT_Mu15_IsoVVVL_PFHT450_PFMET50_accept;
   Int_t           trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_accept;
   Int_t           trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_accept;
   Int_t           trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_accept;
   Int_t           trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_accept;
   Int_t           trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_accept;
   Int_t           trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_accept;
   Int_t           trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_accept;
   Int_t           trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_accept;
   Int_t           trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept;
   Int_t           trig_HLT_Dimuon24_Phi_noCorrL1_accept;
   Int_t           trig_HLT_DoubleIsoMu20_eta2p1_accept;
   Int_t           trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept;
   Int_t           trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept;
   Int_t           trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept;
   Int_t           trig_HLT_Mu8_accept;
   Int_t           trig_HLT_Mu17_accept;
   Int_t           trig_HLT_Mu19_accept;
   Int_t           trig_HLT_Mu17_Photon30_IsoCaloId_accept;
   Int_t           trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept;
   Int_t           trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept;
   Int_t           trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept;
   Int_t           trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_accept;
   Int_t           trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept;
   Int_t           trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept;
   Int_t           trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept;
   Int_t           trig_HLT_Photon20_HoverELoose_accept;
   Int_t           trig_HLT_Photon30_HoverELoose_accept;
   Int_t           trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept;
   Int_t           trig_AlCa_RPCMuonNormalisation_accept;
   Int_t           trig_MC_PFMET_accept;
   Int_t           trig_MC_CaloMET_accept;
   Int_t           trig_MC_CaloMET_JetIdCleaned_accept;
   Int_t           trig_MC_DoubleEle5_CaloIdL_MW_accept;
   Int_t           trig_MC_Ele5_WPTight_Gsf_accept;
   Int_t           trig_MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_accept;
   Int_t           trig_MC_IsoMu_accept;
   Int_t           trig_MC_DoubleMu_TrkIsoVVL_DZ_accept;
   Int_t           trig_MC_DoubleMuNoFiltersNoVtx_accept;
   Int_t           trig_HLT_IsoMu27_MET90_accept;
   Int_t           trig_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_accept;
   Int_t           trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_accept;
   Int_t           trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_accept;
   Int_t           trig_HLT_PFMET100_PFMHT100_IDTight_PFHT60_accept;
   Int_t           trig_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_accept;
   Int_t           trig_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_accept;
   Int_t           trig_HLT_Mu18_Mu9_SameSign_accept;
   Int_t           trig_HLT_Mu18_Mu9_SameSign_DZ_accept;
   Int_t           trig_HLT_Mu18_Mu9_accept;
   Int_t           trig_HLT_Mu18_Mu9_DZ_accept;
   Int_t           trig_HLT_Mu20_Mu10_SameSign_accept;
   Int_t           trig_HLT_Mu20_Mu10_SameSign_DZ_accept;
   Int_t           trig_HLT_Mu20_Mu10_accept;
   Int_t           trig_HLT_Mu20_Mu10_DZ_accept;
   Int_t           trig_HLT_Mu23_Mu12_SameSign_accept;
   Int_t           trig_HLT_Mu23_Mu12_SameSign_DZ_accept;
   Int_t           trig_HLT_Mu23_Mu12_accept;
   Int_t           trig_HLT_Mu23_Mu12_DZ_accept;
   Int_t           trig_HLT_DoubleMu3_DCA_PFMET50_PFMHT60_accept;
   Int_t           trig_HLT_TrkMu6NoFiltersNoVtx_accept;
   Int_t           trig_HLT_TrkMu16NoFiltersNoVtx_accept;
   Int_t           trig_HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_accept;

   // List of branches
   TBranch        *b_ev_event;   //!
   TBranch        *b_ev_run;   //!
   TBranch        *b_ev_luminosityBlock;   //!
   TBranch        *b_ev_time;   //!
   TBranch        *b_ev_time_unixTime;   //!
   TBranch        *b_ev_time_microsecondOffset;   //!
   TBranch        *b_ev_fixedGridRhoAll;   //!
   TBranch        *b_ev_fixedGridRhoFastjetAll;   //!
   TBranch        *b_ev_fixedGridRhoFastjetAllCalo;   //!
   TBranch        *b_ev_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_ev_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_ev_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_ev_prefiringweight;   //!
   TBranch        *b_ev_prefiringweightup;   //!
   TBranch        *b_ev_prefiringweightdown;   //!
   TBranch        *b_LHE_Pt;   //!
   TBranch        *b_LHE_Eta;   //!
   TBranch        *b_LHE_Phi;   //!
   TBranch        *b_LHE_E;   //!
   TBranch        *b_LHE_pdgid;   //!
   TBranch        *b_LHE_status;   //!
   TBranch        *b_LHE_weight_nominal;   //!
   TBranch        *b_LHE_weight_sys;   //!
   TBranch        *b_mc_n;   //!
   TBranch        *b_mc_nMEPartons;   //!
   TBranch        *b_mc_nMEPartonsFiltered;   //!
   TBranch        *b_mc_DJRValues;   //!
   TBranch        *b_mc_weight;   //!
   TBranch        *b_mc_w_sign;   //!
   TBranch        *b_mc_id_first;   //!
   TBranch        *b_mc_id_second;   //!
   TBranch        *b_mc_x_first;   //!
   TBranch        *b_mc_x_second;   //!
   TBranch        *b_mc_xPDF_first;   //!
   TBranch        *b_mc_xPDF_second;   //!
   TBranch        *b_mc_scalePDF;   //!
   TBranch        *b_mc_index;   //!
   TBranch        *b_mc_pdgId;   //!
   TBranch        *b_mc_charge;   //!
   TBranch        *b_mc_status;   //!
   TBranch        *b_mc_status_flags;   //!
   TBranch        *b_mc_status_tau_flags;   //!
   TBranch        *b_mc_tau_charge;   //!
   TBranch        *b_mc_tau_pdgId;   //!
   TBranch        *b_mc_tau_decay;   //!
   TBranch        *b_mc_tau_had_status;   //!
   TBranch        *b_mc_tau_had_charge;   //!
   TBranch        *b_mc_tau_had_pdgId;   //!
   TBranch        *b_mc_mass;   //!
   TBranch        *b_mc_px;   //!
   TBranch        *b_mc_py;   //!
   TBranch        *b_mc_pz;   //!
   TBranch        *b_mc_pt;   //!
   TBranch        *b_mc_eta;   //!
   TBranch        *b_mc_phi;   //!
   TBranch        *b_mc_energy;   //!
   TBranch        *b_mc_tau_pt;   //!
   TBranch        *b_mc_tau_eta;   //!
   TBranch        *b_mc_tau_phi;   //!
   TBranch        *b_mc_tau_energy;   //!
   TBranch        *b_mc_tau_had_pt;   //!
   TBranch        *b_mc_tau_had_eta;   //!
   TBranch        *b_mc_tau_had_phi;   //!
   TBranch        *b_mc_tau_had_energy;   //!
   TBranch        *b_mc_numberOfDaughters;   //!
   TBranch        *b_mc_numberOfMothers;   //!
   TBranch        *b_mc_mother_index;   //!
   TBranch        *b_mc_mother_pdgId;   //!
   TBranch        *b_mc_mother_px;   //!
   TBranch        *b_mc_mother_py;   //!
   TBranch        *b_mc_mother_pz;   //!
   TBranch        *b_mc_mother_pt;   //!
   TBranch        *b_mc_mother_eta;   //!
   TBranch        *b_mc_mother_phi;   //!
   TBranch        *b_mc_mother_energy;   //!
   TBranch        *b_mc_mother_mass;   //!
   TBranch        *b_mc_trueNumInteractions;   //!
   TBranch        *b_mc_PU_NumInteractions;   //!
   TBranch        *b_genjet_pt;   //!
   TBranch        *b_genjet_eta;   //!
   TBranch        *b_genjet_phi;   //!
   TBranch        *b_genjet_energy;   //!
   TBranch        *b_genjetAK8_pt;   //!
   TBranch        *b_genjetAK8_eta;   //!
   TBranch        *b_genjetAK8_phi;   //!
   TBranch        *b_genjetAK8_energy;   //!
   TBranch        *b_gen_weight_sys;   //!
   TBranch        *b_pv_n;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_pv_normalizedChi2;   //!
   TBranch        *b_pv_isValid;   //!
   TBranch        *b_pv_isFake;   //!
   TBranch        *b_gsf_n;   //!
   TBranch        *b_gsf_classification;   //!
   TBranch        *b_gsf_ecalTrkEnergyPostCorr;   //!
   TBranch        *b_gsf_ecalEnergyPostCorr;   //!
   TBranch        *b_gsf_energy;   //!
   TBranch        *b_gsf_p;   //!
   TBranch        *b_gsf_pt;   //!
   TBranch        *b_gsf_et;   //!
   TBranch        *b_gsf_scE1x5;   //!
   TBranch        *b_gsf_scE5x5;   //!
   TBranch        *b_gsf_scE2x5Max;   //!
   TBranch        *b_gsf_full5x5_e5x5;   //!
   TBranch        *b_gsf_full5x5_e1x5;   //!
   TBranch        *b_gsf_full5x5_e2x5Max;   //!
   TBranch        *b_gsf_full5x5_sigmaIetaIeta;   //!
   TBranch        *b_gsf_full5x5_hcalOverEcal;   //!
   TBranch        *b_gsf_eta;   //!
   TBranch        *b_gsf_phi;   //!
   TBranch        *b_gsf_theta;   //!
   TBranch        *b_gsf_px;   //!
   TBranch        *b_gsf_py;   //!
   TBranch        *b_gsf_pz;   //!
   TBranch        *b_gsf_caloEnergy;   //!
   TBranch        *b_gsf_deltaEtaSuperClusterTrackAtVtx;   //!
   TBranch        *b_gsf_deltaPhiSuperClusterTrackAtVtx;   //!
   TBranch        *b_gsf_hadronicOverEm;   //!
   TBranch        *b_gsf_hcalDepth1OverEcal;   //!
   TBranch        *b_gsf_hcalDepth2OverEcal;   //!
   TBranch        *b_gsf_dr03TkSumPt;   //!
   TBranch        *b_gsf_heepTrkPtIso;   //!
   TBranch        *b_gsf_dr03EcalRecHitSumEt;   //!
   TBranch        *b_gsf_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_gsf_dr03HcalDepth2TowerSumEt;   //!
   TBranch        *b_gsf_sigmaIetaIeta;   //!
   TBranch        *b_gsf_charge;   //!
   TBranch        *b_gsf_ecaldrivenSeed;   //!
   TBranch        *b_gsf_trackerdrivenSeed;   //!
   TBranch        *b_gsf_isEB;   //!
   TBranch        *b_gsf_isEE;   //!
   TBranch        *b_gsf_isPF;   //!
   TBranch        *b_gsf_passConversionVeto;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Fall17_94X_V2_loose;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Fall17_94X_V2_medium;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight;   //!
   TBranch        *b_gsf_VID_cutBasedElectronID_Fall17_94X_V2_veto;   //!
   TBranch        *b_gsf_VID_heepElectronID_HEEPV70;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_iso_V2_wp80;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_iso_V2_wp90;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_iso_V2_wpLoose;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_noIso_V2_wp80;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_noIso_V2_wp90;   //!
   TBranch        *b_gsf_VID_mvaEleID_Fall17_noIso_V2_wpLoose;   //!
   TBranch        *b_gsf_deltaEtaSeedClusterTrackAtCalo;   //!
   TBranch        *b_gsf_deltaPhiSeedClusterTrackAtCalo;   //!
   TBranch        *b_gsf_ecalEnergy;   //!
   TBranch        *b_gsf_eSuperClusterOverP;   //!
   TBranch        *b_gsf_dxy;   //!
   TBranch        *b_gsf_dxy_beamSpot;   //!
   TBranch        *b_gsf_dxy_firstPVtx;   //!
   TBranch        *b_gsf_dxyError;   //!
   TBranch        *b_gsf_dz;   //!
   TBranch        *b_gsf_dz_beamSpot;   //!
   TBranch        *b_gsf_dz_firstPVtx;   //!
   TBranch        *b_gsf_dzError;   //!
   TBranch        *b_gsf_vz;   //!
   TBranch        *b_gsf_numberOfValidHits;   //!
   TBranch        *b_gsf_nLostInnerHits;   //!
   TBranch        *b_gsf_nLostOuterHits;   //!
   TBranch        *b_gsf_convFlags;   //!
   TBranch        *b_gsf_convDist;   //!
   TBranch        *b_gsf_convDcot;   //!
   TBranch        *b_gsf_convRadius;   //!
   TBranch        *b_gsf_fBrem;   //!
   TBranch        *b_gsf_e1x5;   //!
   TBranch        *b_gsf_e2x5Max;   //!
   TBranch        *b_gsf_e5x5;   //!
   TBranch        *b_gsf_r9;   //!
   TBranch        *b_gsf_deltaEtaSeedClusterTrackAtVtx;   //!
   TBranch        *b_gsf_effArea;   //!
   TBranch        *b_gsf_sumChargedHadronPt;   //!
   TBranch        *b_gsf_sumNeutralHadronEt;   //!
   TBranch        *b_gsf_sumPhotonEt;   //!
   TBranch        *b_gsf_sumPUPt;   //!
   TBranch        *b_gsf_ooEmooP;   //!
   TBranch        *b_gsf_hitsinfo;   //!
   TBranch        *b_gsf_pixelMatch_dPhi1;   //!
   TBranch        *b_gsf_pixelMatch_dPhi2;   //!
   TBranch        *b_gsf_pixelMatch_dRz1;   //!
   TBranch        *b_gsf_pixelMatch_dRz2;   //!
   TBranch        *b_gsf_pixelMatch_subDetector1;   //!
   TBranch        *b_gsf_pixelMatch_subDetector2;   //!
   TBranch        *b_gsf_mc_bestDR;   //!
   TBranch        *b_gsf_mc_index;   //!
   TBranch        *b_gsf_mc_ERatio;   //!
   TBranch        *b_gsf_sc_energy;   //!
   TBranch        *b_gsf_sc_seed_eta;   //!
   TBranch        *b_gsf_sc_eta;   //!
   TBranch        *b_gsf_sc_etacorr;   //!
   TBranch        *b_gsf_sc_theta;   //!
   TBranch        *b_gsf_sc_thetacorr;   //!
   TBranch        *b_gsf_sc_et;   //!
   TBranch        *b_gsf_sc_phi;   //!
   TBranch        *b_gsf_sc_px;   //!
   TBranch        *b_gsf_sc_py;   //!
   TBranch        *b_gsf_sc_pz;   //!
   TBranch        *b_gsf_sc_x;   //!
   TBranch        *b_gsf_sc_y;   //!
   TBranch        *b_gsf_sc_z;   //!
   TBranch        *b_gsf_sc_phiWidth;   //!
   TBranch        *b_gsf_sc_etaWidth;   //!
   TBranch        *b_gsf_sc_seed_rawId;   //!
   TBranch        *b_gsf_sc_seed_ieta;   //!
   TBranch        *b_gsf_sc_seed_iphi;   //!
   TBranch        *b_gsf_sc_seed_kHasSwitchToGain6;   //!
   TBranch        *b_gsf_sc_seed_kHasSwitchToGain1;   //!
   TBranch        *b_EHits_isSaturated;   //!
   TBranch        *b_EBHits_rawId;   //!
   TBranch        *b_EBHits_iRechit;   //!
   TBranch        *b_EBHits_ieta;   //!
   TBranch        *b_EBHits_iphi;   //!
   TBranch        *b_EBHits_RecoFlag;   //!
   TBranch        *b_EBHits_kSaturated;   //!
   TBranch        *b_EBHits_kLeadingEdgeRecovered;   //!
   TBranch        *b_EBHits_kNeighboursRecovered;   //!
   TBranch        *b_EBHits_kWeird;   //!
   TBranch        *b_EBHits_energy;   //!
   TBranch        *b_EEHits_rawId;   //!
   TBranch        *b_EEHits_iRechit;   //!
   TBranch        *b_EEHits_ieta;   //!
   TBranch        *b_EEHits_iphi;   //!
   TBranch        *b_EEHits_RecoFlag;   //!
   TBranch        *b_EEHits_kSaturated;   //!
   TBranch        *b_EEHits_kLeadingEdgeRecovered;   //!
   TBranch        *b_EEHits_kNeighboursRecovered;   //!
   TBranch        *b_EEHits_kWeird;   //!
   TBranch        *b_EEHits_energy;   //!
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_gt_qoverp;   //!
   TBranch        *b_mu_gt_charge;   //!
   TBranch        *b_mu_gt_pt;   //!
   TBranch        *b_mu_gt_eta;   //!
   TBranch        *b_mu_gt_phi;   //!
   TBranch        *b_mu_gt_p;   //!
   TBranch        *b_mu_gt_px;   //!
   TBranch        *b_mu_gt_py;   //!
   TBranch        *b_mu_gt_pz;   //!
   TBranch        *b_mu_gt_theta;   //!
   TBranch        *b_mu_gt_lambda;   //!
   TBranch        *b_mu_gt_d0;   //!
   TBranch        *b_mu_gt_dz;   //!
   TBranch        *b_mu_gt_dz_beamspot;   //!
   TBranch        *b_mu_gt_dz_firstPVtx;   //!
   TBranch        *b_mu_gt_dxy;   //!
   TBranch        *b_mu_gt_dxy_beamspot;   //!
   TBranch        *b_mu_gt_dxy_firstPVtx;   //!
   TBranch        *b_mu_gt_dsz;   //!
   TBranch        *b_mu_gt_vx;   //!
   TBranch        *b_mu_gt_vy;   //!
   TBranch        *b_mu_gt_vz;   //!
   TBranch        *b_mu_gt_qoverpError;   //!
   TBranch        *b_mu_gt_ptError;   //!
   TBranch        *b_mu_gt_thetaError;   //!
   TBranch        *b_mu_gt_lambdaError;   //!
   TBranch        *b_mu_gt_phiError;   //!
   TBranch        *b_mu_gt_dxyError;   //!
   TBranch        *b_mu_gt_d0Error;   //!
   TBranch        *b_mu_gt_dszError;   //!
   TBranch        *b_mu_gt_dzError;   //!
   TBranch        *b_mu_gt_etaError;   //!
   TBranch        *b_mu_gt_chi2;   //!
   TBranch        *b_mu_gt_ndof;   //!
   TBranch        *b_mu_gt_normalizedChi2;   //!
   TBranch        *b_mu_ot_qoverp;   //!
   TBranch        *b_mu_ot_charge;   //!
   TBranch        *b_mu_ot_pt;   //!
   TBranch        *b_mu_ot_eta;   //!
   TBranch        *b_mu_ot_phi;   //!
   TBranch        *b_mu_ot_p;   //!
   TBranch        *b_mu_ot_px;   //!
   TBranch        *b_mu_ot_py;   //!
   TBranch        *b_mu_ot_pz;   //!
   TBranch        *b_mu_ot_theta;   //!
   TBranch        *b_mu_ot_lambda;   //!
   TBranch        *b_mu_ot_d0;   //!
   TBranch        *b_mu_ot_dz;   //!
   TBranch        *b_mu_ot_dz_beamspot;   //!
   TBranch        *b_mu_ot_dz_firstPVtx;   //!
   TBranch        *b_mu_ot_dxy;   //!
   TBranch        *b_mu_ot_dxy_beamspot;   //!
   TBranch        *b_mu_ot_dxy_firstPVtx;   //!
   TBranch        *b_mu_ot_dsz;   //!
   TBranch        *b_mu_ot_vx;   //!
   TBranch        *b_mu_ot_vy;   //!
   TBranch        *b_mu_ot_vz;   //!
   TBranch        *b_mu_ot_qoverpError;   //!
   TBranch        *b_mu_ot_ptError;   //!
   TBranch        *b_mu_ot_thetaError;   //!
   TBranch        *b_mu_ot_lambdaError;   //!
   TBranch        *b_mu_ot_phiError;   //!
   TBranch        *b_mu_ot_dxyError;   //!
   TBranch        *b_mu_ot_d0Error;   //!
   TBranch        *b_mu_ot_dszError;   //!
   TBranch        *b_mu_ot_dzError;   //!
   TBranch        *b_mu_ot_etaError;   //!
   TBranch        *b_mu_ot_chi2;   //!
   TBranch        *b_mu_ot_ndof;   //!
   TBranch        *b_mu_ot_normalizedChi2;   //!
   TBranch        *b_mu_it_qoverp;   //!
   TBranch        *b_mu_it_charge;   //!
   TBranch        *b_mu_it_pt;   //!
   TBranch        *b_mu_it_eta;   //!
   TBranch        *b_mu_it_phi;   //!
   TBranch        *b_mu_it_p;   //!
   TBranch        *b_mu_it_px;   //!
   TBranch        *b_mu_it_py;   //!
   TBranch        *b_mu_it_pz;   //!
   TBranch        *b_mu_it_theta;   //!
   TBranch        *b_mu_it_lambda;   //!
   TBranch        *b_mu_it_d0;   //!
   TBranch        *b_mu_it_dz;   //!
   TBranch        *b_mu_it_dz_beamspot;   //!
   TBranch        *b_mu_it_dz_firstPVtx;   //!
   TBranch        *b_mu_it_dxy;   //!
   TBranch        *b_mu_it_dxy_beamspot;   //!
   TBranch        *b_mu_it_dxy_firstPVtx;   //!
   TBranch        *b_mu_it_dsz;   //!
   TBranch        *b_mu_it_vx;   //!
   TBranch        *b_mu_it_vy;   //!
   TBranch        *b_mu_it_vz;   //!
   TBranch        *b_mu_it_qoverpError;   //!
   TBranch        *b_mu_it_ptError;   //!
   TBranch        *b_mu_it_thetaError;   //!
   TBranch        *b_mu_it_lambdaError;   //!
   TBranch        *b_mu_it_phiError;   //!
   TBranch        *b_mu_it_dxyError;   //!
   TBranch        *b_mu_it_d0Error;   //!
   TBranch        *b_mu_it_dszError;   //!
   TBranch        *b_mu_it_dzError;   //!
   TBranch        *b_mu_it_etaError;   //!
   TBranch        *b_mu_it_chi2;   //!
   TBranch        *b_mu_it_ndof;   //!
   TBranch        *b_mu_it_normalizedChi2;   //!
   TBranch        *b_mu_ibt_qoverp;   //!
   TBranch        *b_mu_ibt_charge;   //!
   TBranch        *b_mu_ibt_pt;   //!
   TBranch        *b_mu_ibt_eta;   //!
   TBranch        *b_mu_ibt_phi;   //!
   TBranch        *b_mu_ibt_p;   //!
   TBranch        *b_mu_ibt_px;   //!
   TBranch        *b_mu_ibt_py;   //!
   TBranch        *b_mu_ibt_pz;   //!
   TBranch        *b_mu_ibt_theta;   //!
   TBranch        *b_mu_ibt_lambda;   //!
   TBranch        *b_mu_ibt_d0;   //!
   TBranch        *b_mu_ibt_dz;   //!
   TBranch        *b_mu_ibt_dz_beamspot;   //!
   TBranch        *b_mu_ibt_dz_firstPVtx;   //!
   TBranch        *b_mu_ibt_dxy;   //!
   TBranch        *b_mu_ibt_dxy_beamspot;   //!
   TBranch        *b_mu_ibt_dxy_firstPVtx;   //!
   TBranch        *b_mu_ibt_dsz;   //!
   TBranch        *b_mu_ibt_vx;   //!
   TBranch        *b_mu_ibt_vy;   //!
   TBranch        *b_mu_ibt_vz;   //!
   TBranch        *b_mu_ibt_qoverpError;   //!
   TBranch        *b_mu_ibt_ptError;   //!
   TBranch        *b_mu_ibt_thetaError;   //!
   TBranch        *b_mu_ibt_lambdaError;   //!
   TBranch        *b_mu_ibt_phiError;   //!
   TBranch        *b_mu_ibt_dxyError;   //!
   TBranch        *b_mu_ibt_d0Error;   //!
   TBranch        *b_mu_ibt_dszError;   //!
   TBranch        *b_mu_ibt_dzError;   //!
   TBranch        *b_mu_ibt_etaError;   //!
   TBranch        *b_mu_ibt_chi2;   //!
   TBranch        *b_mu_ibt_ndof;   //!
   TBranch        *b_mu_ibt_normalizedChi2;   //!
   TBranch        *b_mu_isGlobalMuon;   //!
   TBranch        *b_mu_isStandAloneMuon;   //!
   TBranch        *b_mu_isTrackerMuon;   //!
   TBranch        *b_mu_isPFMuon;   //!
   TBranch        *b_mu_isPFIsolationValid;   //!
   TBranch        *b_mu_isGoodMuonTMLastStationLoose;   //!
   TBranch        *b_mu_isGoodMuonTMLastStationTight;   //!
   TBranch        *b_mu_isGoodMuonTM2DCompatibilityLoose;   //!
   TBranch        *b_mu_isGoodMuonTM2DCompatibilityTight;   //!
   TBranch        *b_mu_isGoodMuonTMOneStationLoose;   //!
   TBranch        *b_mu_isGoodMuonTMOneStationTight;   //!
   TBranch        *b_mu_isGoodMuonTMLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_mu_isGoodMuonTMLastStationOptimizedLowPtTight;   //!
   TBranch        *b_mu_isTightMuon;   //!
   TBranch        *b_mu_isMediumMuon;   //!
   TBranch        *b_mu_isLooseMuon;   //!
   TBranch        *b_mu_isSoftMuon;   //!
   TBranch        *b_mu_isHighPtMuonOld;   //!
   TBranch        *b_mu_isHighPtMuon;   //!
   TBranch        *b_mu_isTrackerHighPtMuon;   //!
   TBranch        *b_mu_CutBasedIdLoose;   //!
   TBranch        *b_mu_CutBasedIdMedium;   //!
   TBranch        *b_mu_CutBasedIdMediumPrompt;   //!
   TBranch        *b_mu_CutBasedIdTight;   //!
   TBranch        *b_mu_CutBasedIdGlobalHighPt;   //!
   TBranch        *b_mu_CutBasedIdTrkHighPt;   //!
   TBranch        *b_mu_PFIsoVeryLoose;   //!
   TBranch        *b_mu_PFIsoLoose;   //!
   TBranch        *b_mu_PFIsoMedium;   //!
   TBranch        *b_mu_PFIsoTight;   //!
   TBranch        *b_mu_PFIsoVeryTight;   //!
   TBranch        *b_mu_PFIsoVeryVeryTight;   //!
   TBranch        *b_mu_TkIsoLoose;   //!
   TBranch        *b_mu_TkIsoTight;   //!
   TBranch        *b_mu_SoftCutBasedId;   //!
   TBranch        *b_mu_SoftMvaId;   //!
   TBranch        *b_mu_MvaLoose;   //!
   TBranch        *b_mu_MvaMedium;   //!
   TBranch        *b_mu_MvaTight;   //!
   TBranch        *b_mu_MiniIsoLoose;   //!
   TBranch        *b_mu_MiniIsoMedium;   //!
   TBranch        *b_mu_MiniIsoTight;   //!
   TBranch        *b_mu_MiniIsoVeryTight;   //!
   TBranch        *b_mu_TriggerIdLoose;   //!
   TBranch        *b_mu_InTimeMuon;   //!
   TBranch        *b_mu_MultiIsoLoose;   //!
   TBranch        *b_mu_MultiIsoMedium;   //!
   TBranch        *b_mu_numberOfMatchedStations;   //!
   TBranch        *b_mu_numberOfValidPixelHits;   //!
   TBranch        *b_mu_trackerLayersWithMeasurement;   //!
   TBranch        *b_mu_numberOfValidMuonHits;   //!
   TBranch        *b_mu_pixelLayersWithMeasurement;   //!
   TBranch        *b_mu_innerTrack_validFraction;   //!
   TBranch        *b_mu_combinedQuality_trkKink;   //!
   TBranch        *b_mu_combinedQuality_chi2LocalPosition;   //!
   TBranch        *b_mu_segmentCompatibility;   //!
   TBranch        *b_mu_dB;   //!
   TBranch        *b_mu_pt_default;   //!
   TBranch        *b_mu_isolationR03_sumPt;   //!
   TBranch        *b_mu_isolationR03_trackerVetoPt;   //!
   TBranch        *b_mu_isolationR03_emEt;   //!
   TBranch        *b_mu_isolationR03_emVetoEt;   //!
   TBranch        *b_mu_isolationR03_hadEt;   //!
   TBranch        *b_mu_isolationR03_hadVetoEt;   //!
   TBranch        *b_mu_isolationR05_sumPt;   //!
   TBranch        *b_mu_isolationR05_trackerVetoPt;   //!
   TBranch        *b_mu_isolationR05_emEt;   //!
   TBranch        *b_mu_isolationR05_emVetoEt;   //!
   TBranch        *b_mu_isolationR05_hadEt;   //!
   TBranch        *b_mu_isolationR05_hadVetoEt;   //!
   TBranch        *b_mu_pfIsolationR03_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfIsolationR03_sumNeutralHadronEt;   //!
   TBranch        *b_mu_pfIsolationR03_sumChargedParticlePt;   //!
   TBranch        *b_mu_pfIsolationR03_sumPhotonEt;   //!
   TBranch        *b_mu_pfIsolationR03_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR03_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR03_sumPUPt;   //!
   TBranch        *b_mu_pfIsolationR04_sumChargedHadronPt;   //!
   TBranch        *b_mu_pfIsolationR04_sumNeutralHadronEt;   //!
   TBranch        *b_mu_pfIsolationR04_sumChargedParticlePt;   //!
   TBranch        *b_mu_pfIsolationR04_sumPhotonEt;   //!
   TBranch        *b_mu_pfIsolationR04_sumNeutralHadronEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR04_sumPhotonEtHighThreshold;   //!
   TBranch        *b_mu_pfIsolationR04_sumPUPt;   //!
   TBranch        *b_mu_pfIsoDbCorrected03;   //!
   TBranch        *b_mu_pfIsoDbCorrected04;   //!
   TBranch        *b_mu_isoTrackerBased03;   //!
   TBranch        *b_mu_mc_bestDR;   //!
   TBranch        *b_mu_mc_index;   //!
   TBranch        *b_mu_mc_ERatio;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_jet_px;   //!
   TBranch        *b_jet_py;   //!
   TBranch        *b_jet_pz;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_theta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_energy;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_jet_chargedEmEnergyFraction;   //!
   TBranch        *b_jet_neutralHadronEnergyFraction;   //!
   TBranch        *b_jet_neutralEmEnergyFraction;   //!
   TBranch        *b_jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_jet_muonEnergyFraction;   //!
   TBranch        *b_jet_chargedMultiplicity;   //!
   TBranch        *b_jet_neutralMultiplicity;   //!
   TBranch        *b_jet_partonFlavour;   //!
   TBranch        *b_jet_hadronFlavour;   //!
   TBranch        *b_jet_CSVv2;   //!
   TBranch        *b_jet_CvsL;   //!
   TBranch        *b_jet_CvsB;   //!
   TBranch        *b_jet_MVA2BJets;   //!
   TBranch        *b_jet_CvsB_DeepJet_charm_tagger;   //!
   TBranch        *b_jet_CvsL_DeepJet_charm_tagger;   //!
   TBranch        *b_jet_CvsB_DeepCSV_charm_tagger;   //!
   TBranch        *b_jet_CvsL_DeepCSV_charm_tagger;   //!
   TBranch        *b_jet_DeepJet;   //!
   TBranch        *b_jet_DeepCSV;   //!
   TBranch        *b_jet_isJetIDLoose_2016;   //!
   TBranch        *b_jet_isJetIDTight_2016;   //!
   TBranch        *b_jet_isJetIDTightLepVeto_2016;   //!
   TBranch        *b_jet_isJetID_2017;   //!
   TBranch        *b_jet_isJetIDLepVeto_2017;   //!
   TBranch        *b_jet_isJetID_2018;   //!
   TBranch        *b_jet_isJetIDLepVeto_2018;   //!
   TBranch        *b_jet_Smeared_pt;   //!
   TBranch        *b_jet_SmearedJetResUp_pt;   //!
   TBranch        *b_jet_SmearedJetResDown_pt;   //!
   TBranch        *b_jet_SmearedJetEnUp_pt;   //!
   TBranch        *b_jet_SmearedJetEnDown_pt;   //!
   TBranch        *b_BtagSF_Deepcsv_loose;   //!
   TBranch        *b_BtagSF_Deepcsv_Up_loose;   //!
   TBranch        *b_BtagSF_Deepcsv_Down_loose;   //!
   TBranch        *b_BtagSF_Deepcsv_medium;   //!
   TBranch        *b_BtagSF_Deepcsv_Up_medium;   //!
   TBranch        *b_BtagSF_Deepcsv_Down_medium;   //!
   TBranch        *b_BtagSF_Deepcsv_tight;   //!
   TBranch        *b_BtagSF_Deepcsv_Up_tight;   //!
   TBranch        *b_BtagSF_Deepcsv_Down_tight;   //!
   TBranch        *b_MET_nominal_Pt;   //!
   TBranch        *b_MET_nominal_phi;   //!
   TBranch        *b_MET_nominal_significance;   //!
   TBranch        *b_MET_Pt;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_significance;   //!
   TBranch        *b_MET_T1_Pt;   //!
   TBranch        *b_MET_T1_phi;   //!
   TBranch        *b_MET_T1_significance;   //!
   TBranch        *b_MET_gen_pt;   //!
   TBranch        *b_MET_gen_phi;   //!
   TBranch        *b_MET_Type1Unc_Px;   //!
   TBranch        *b_MET_Type1Unc_Py;   //!
   TBranch        *b_MET_Type1Unc_Pt;   //!
   TBranch        *b_MET_T1JetEnDown_Pt;   //!
   TBranch        *b_MET_T1JetEnDown_phi;   //!
   TBranch        *b_MET_T1JetEnDown_significance;   //!
   TBranch        *b_MET_T1JetEnUp_Pt;   //!
   TBranch        *b_MET_T1JetEnUp_phi;   //!
   TBranch        *b_MET_T1JetEnUp_significance;   //!
   TBranch        *b_MET_T1Smear_Pt;   //!
   TBranch        *b_MET_T1Smear_phi;   //!
   TBranch        *b_MET_T1Smear_significance;   //!
   TBranch        *b_MET_T1SmearJetEnDown_Pt;   //!
   TBranch        *b_MET_T1SmearJetEnDown_phi;   //!
   TBranch        *b_MET_T1SmearJetEnDown_significance;   //!
   TBranch        *b_MET_T1SmearJetEnUp_Pt;   //!
   TBranch        *b_MET_T1SmearJetEnUp_phi;   //!
   TBranch        *b_MET_T1SmearJetEnUp_significance;   //!
   TBranch        *b_MET_T1SmearJetResDown_Pt;   //!
   TBranch        *b_MET_T1SmearJetResDown_phi;   //!
   TBranch        *b_MET_T1SmearJetResDown_significance;   //!
   TBranch        *b_MET_T1SmearJetResUp_Pt;   //!
   TBranch        *b_MET_T1SmearJetResUp_phi;   //!
   TBranch        *b_MET_T1SmearJetResUp_significance;   //!
   TBranch        *b_MET_T1Txy_Pt;   //!
   TBranch        *b_MET_T1Txy_phi;   //!
   TBranch        *b_MET_T1Txy_significance;   //!
   TBranch        *b_MET_FinalCollection_Pt;   //!
   TBranch        *b_MET_FinalCollection_phi;   //!
   TBranch        *b_MET_FinalCollection_significance;   //!
   TBranch        *b_trig_Flag_ecalBadCalibReduced;   //!
   TBranch        *b_trig_Flag_HBHENoiseFilter_accept;   //!
   TBranch        *b_trig_Flag_HBHENoiseIsoFilter_accept;   //!
   TBranch        *b_trig_Flag_CSCTightHaloFilter_accept;   //!
   TBranch        *b_trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept;   //!
   TBranch        *b_trig_Flag_CSCTightHalo2015Filter_accept;   //!
   TBranch        *b_trig_Flag_globalTightHalo2016Filter_accept;   //!
   TBranch        *b_trig_Flag_globalSuperTightHalo2016Filter_accept;   //!
   TBranch        *b_trig_Flag_HcalStripHaloFilter_accept;   //!
   TBranch        *b_trig_Flag_hcalLaserEventFilter_accept;   //!
   TBranch        *b_trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept;   //!
   TBranch        *b_trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept;   //!
   TBranch        *b_trig_Flag_ecalBadCalibFilter_accept;   //!
   TBranch        *b_trig_Flag_goodVertices_accept;   //!
   TBranch        *b_trig_Flag_eeBadScFilter_accept;   //!
   TBranch        *b_trig_Flag_ecalLaserCorrFilter_accept;   //!
   TBranch        *b_trig_Flag_trkPOGFilters_accept;   //!
   TBranch        *b_trig_Flag_chargedHadronTrackResolutionFilter_accept;   //!
   TBranch        *b_trig_Flag_muonBadTrackFilter_accept;   //!
   TBranch        *b_trig_Flag_BadChargedCandidateFilter_accept;   //!
   TBranch        *b_trig_Flag_BadPFMuonFilter_accept;   //!
   TBranch        *b_trig_Flag_BadChargedCandidateSummer16Filter_accept;   //!
   TBranch        *b_trig_Flag_BadPFMuonSummer16Filter_accept;   //!
   TBranch        *b_trig_Flag_trkPOG_manystripclus53X_accept;   //!
   TBranch        *b_trig_Flag_trkPOG_toomanystripclus53X_accept;   //!
   TBranch        *b_trig_Flag_trkPOG_logErrorTooManyClusters_accept;   //!
   TBranch        *b_trig_Flag_METFilters_accept;   //!
   TBranch        *b_trig_raw2digi_step_accept;   //!
   TBranch        *b_trig_L1Reco_step_accept;   //!
   TBranch        *b_trig_reconstruction_step_accept;   //!
   TBranch        *b_trig_recosim_step_accept;   //!
   TBranch        *b_trig_eventinterpretaion_step_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle25_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle27_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle33_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele27_Ele37_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_Mu27_Ele37_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_Mu37_Ele27_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_HLT_Mu37_TkMu27_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu4_3_Bs_accept;   //!
   TBranch        *b_trig_HLT_Ele20_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele15_WPLoose_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele17_WPLoose_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele20_WPLoose_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept;   //!
   TBranch        *b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept;   //!
   TBranch        *b_trig_HLT_Ele27_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele28_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele30_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept;   //!
   TBranch        *b_trig_HLT_Ele38_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele40_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept;   //!
   TBranch        *b_trig_HLT_IsoMu20_accept;   //!
   TBranch        *b_trig_HLT_IsoMu24_accept;   //!
   TBranch        *b_trig_HLT_IsoTkMu24_accept;   //!
   TBranch        *b_trig_HLT_IsoMu24_eta2p1_accept;   //!
   TBranch        *b_trig_HLT_IsoMu27_accept;   //!
   TBranch        *b_trig_HLT_IsoMu30_accept;   //!
   TBranch        *b_trig_HLT_L1SingleMu18_accept;   //!
   TBranch        *b_trig_HLT_L1SingleMu25_accept;   //!
   TBranch        *b_trig_HLT_L2Mu10_accept;   //!
   TBranch        *b_trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept;   //!
   TBranch        *b_trig_HLT_L2Mu10_NoVertex_NoBPTX_accept;   //!
   TBranch        *b_trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept;   //!
   TBranch        *b_trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept;   //!
   TBranch        *b_trig_HLT_L2Mu50_accept;   //!
   TBranch        *b_trig_HLT_L2Mu23NoVtx_2Cha_accept;   //!
   TBranch        *b_trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu50_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_accept;   //!
   TBranch        *b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept;   //!
   TBranch        *b_trig_HLT_Mu25_TkMu0_Onia_accept;   //!
   TBranch        *b_trig_HLT_Mu30_TkMu0_Psi_accept;   //!
   TBranch        *b_trig_HLT_Mu20_TkMu0_Phi_accept;   //!
   TBranch        *b_trig_HLT_Mu25_TkMu0_Phi_accept;   //!
   TBranch        *b_trig_HLT_Mu12_accept;   //!
   TBranch        *b_trig_HLT_Mu15_accept;   //!
   TBranch        *b_trig_HLT_Mu20_accept;   //!
   TBranch        *b_trig_HLT_Mu27_accept;   //!
   TBranch        *b_trig_HLT_Mu50_accept;   //!
   TBranch        *b_trig_HLT_Mu55_accept;   //!
   TBranch        *b_trig_HLT_OldMu100_accept;   //!
   TBranch        *b_trig_HLT_TkMu100_accept;   //!
   TBranch        *b_trig_HLT_PFHT500_PFMET100_PFMHT100_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFHT500_PFMET110_PFMHT110_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFHT700_PFMET85_PFMHT85_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFHT700_PFMET95_PFMHT95_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFHT800_PFMET75_PFMHT75_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFHT800_PFMET85_PFMHT85_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMET110_PFMHT110_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMET120_PFMHT120_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMET130_PFMHT130_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMET140_PFMHT140_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_accept;   //!
   TBranch        *b_trig_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_accept;   //!
   TBranch        *b_trig_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_accept;   //!
   TBranch        *b_trig_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_accept;   //!
   TBranch        *b_trig_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_accept;   //!
   TBranch        *b_trig_HLT_PFMET120_PFMHT120_IDTight_PFHT60_accept;   //!
   TBranch        *b_trig_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_accept;   //!
   TBranch        *b_trig_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_accept;   //!
   TBranch        *b_trig_HLT_PFMETTypeOne110_PFMHT110_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMETTypeOne120_PFMHT120_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMETTypeOne130_PFMHT130_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMETTypeOne140_PFMHT140_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_accept;   //!
   TBranch        *b_trig_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_accept;   //!
   TBranch        *b_trig_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_accept;   //!
   TBranch        *b_trig_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_accept;   //!
   TBranch        *b_trig_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_accept;   //!
   TBranch        *b_trig_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_accept;   //!
   TBranch        *b_trig_HLT_CaloMET80_NotCleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET90_NotCleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET100_NotCleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET110_NotCleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET250_NotCleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET70_HBHECleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET80_HBHECleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET90_HBHECleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET100_HBHECleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET250_HBHECleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET300_HBHECleaned_accept;   //!
   TBranch        *b_trig_HLT_CaloMET350_HBHECleaned_accept;   //!
   TBranch        *b_trig_HLT_PFMET200_NotCleaned_accept;   //!
   TBranch        *b_trig_HLT_PFMET200_HBHECleaned_accept;   //!
   TBranch        *b_trig_HLT_PFMET250_HBHECleaned_accept;   //!
   TBranch        *b_trig_HLT_PFMET300_HBHECleaned_accept;   //!
   TBranch        *b_trig_HLT_PFMET200_HBHE_BeamHaloCleaned_accept;   //!
   TBranch        *b_trig_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_accept;   //!
   TBranch        *b_trig_HLT_MET105_IsoTrk50_accept;   //!
   TBranch        *b_trig_HLT_MET120_IsoTrk50_accept;   //!
   TBranch        *b_trig_HLT_Photon300_NoHE_accept;   //!
   TBranch        *b_trig_HLT_Mu8_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu17_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Mu19_TrkIsoVVL_accept;   //!
   TBranch        *b_trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept;   //!
   TBranch        *b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_HLT_Photon20_accept;   //!
   TBranch        *b_trig_HLT_Photon33_accept;   //!
   TBranch        *b_trig_HLT_Photon50_accept;   //!
   TBranch        *b_trig_HLT_Photon75_accept;   //!
   TBranch        *b_trig_HLT_Photon90_accept;   //!
   TBranch        *b_trig_HLT_Photon120_accept;   //!
   TBranch        *b_trig_HLT_Photon150_accept;   //!
   TBranch        *b_trig_HLT_Photon175_accept;   //!
   TBranch        *b_trig_HLT_Photon200_accept;   //!
   TBranch        *b_trig_HLT_Photon100EB_TightID_TightIso_accept;   //!
   TBranch        *b_trig_HLT_Photon110EB_TightID_TightIso_accept;   //!
   TBranch        *b_trig_HLT_Photon120EB_TightID_TightIso_accept;   //!
   TBranch        *b_trig_HLT_Photon100EBHE10_accept;   //!
   TBranch        *b_trig_HLT_Photon100EEHE10_accept;   //!
   TBranch        *b_trig_HLT_Photon100EE_TightID_TightIso_accept;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_accept;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_accept;   //!
   TBranch        *b_trig_HLT_Photon90_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon120_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon165_R9Id90_HE10_IsoM_accept;   //!
   TBranch        *b_trig_HLT_Photon90_CaloIdL_PFHT700_accept;   //!
   TBranch        *b_trig_HLT_Photon35_TwoProngs35_accept;   //!
   TBranch        *b_trig_HLT_IsoMu24_TwoProngs35_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu3_DZ_PFMET70_PFMHT70_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu3_DZ_PFMET90_PFMHT90_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu43NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu48NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept;   //!
   TBranch        *b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept;   //!
   TBranch        *b_trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_accept;   //!
   TBranch        *b_trig_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu33NoFiltersNoVtxDisplaced_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu40NoFiltersNoVtxDisplaced_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept;   //!
   TBranch        *b_trig_HLT_DiJet110_35_Mjj650_PFMET110_accept;   //!
   TBranch        *b_trig_HLT_DiJet110_35_Mjj650_PFMET120_accept;   //!
   TBranch        *b_trig_HLT_DiJet110_35_Mjj650_PFMET130_accept;   //!
   TBranch        *b_trig_HLT_TripleJet110_35_35_Mjj650_PFMET110_accept;   //!
   TBranch        *b_trig_HLT_TripleJet110_35_35_Mjj650_PFMET120_accept;   //!
   TBranch        *b_trig_HLT_TripleJet110_35_35_Mjj650_PFMET130_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept;   //!
   TBranch        *b_trig_HLT_Ele15_IsoVVVL_PFHT450_PFMET50_accept;   //!
   TBranch        *b_trig_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_accept;   //!
   TBranch        *b_trig_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_accept;   //!
   TBranch        *b_trig_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_accept;   //!
   TBranch        *b_trig_HLT_Mu15_IsoVVVL_PFHT450_PFMET50_accept;   //!
   TBranch        *b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_accept;   //!
   TBranch        *b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_accept;   //!
   TBranch        *b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_accept;   //!
   TBranch        *b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_accept;   //!
   TBranch        *b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_accept;   //!
   TBranch        *b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_accept;   //!
   TBranch        *b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_accept;   //!
   TBranch        *b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_accept;   //!
   TBranch        *b_trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept;   //!
   TBranch        *b_trig_HLT_Dimuon24_Phi_noCorrL1_accept;   //!
   TBranch        *b_trig_HLT_DoubleIsoMu20_eta2p1_accept;   //!
   TBranch        *b_trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_Mu8_accept;   //!
   TBranch        *b_trig_HLT_Mu17_accept;   //!
   TBranch        *b_trig_HLT_Mu19_accept;   //!
   TBranch        *b_trig_HLT_Mu17_Photon30_IsoCaloId_accept;   //!
   TBranch        *b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept;   //!
   TBranch        *b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_accept;   //!
   TBranch        *b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept;   //!
   TBranch        *b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept;   //!
   TBranch        *b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept;   //!
   TBranch        *b_trig_HLT_Photon20_HoverELoose_accept;   //!
   TBranch        *b_trig_HLT_Photon30_HoverELoose_accept;   //!
   TBranch        *b_trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept;   //!
   TBranch        *b_trig_AlCa_RPCMuonNormalisation_accept;   //!
   TBranch        *b_trig_MC_PFMET_accept;   //!
   TBranch        *b_trig_MC_CaloMET_accept;   //!
   TBranch        *b_trig_MC_CaloMET_JetIdCleaned_accept;   //!
   TBranch        *b_trig_MC_DoubleEle5_CaloIdL_MW_accept;   //!
   TBranch        *b_trig_MC_Ele5_WPTight_Gsf_accept;   //!
   TBranch        *b_trig_MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_accept;   //!
   TBranch        *b_trig_MC_IsoMu_accept;   //!
   TBranch        *b_trig_MC_DoubleMu_TrkIsoVVL_DZ_accept;   //!
   TBranch        *b_trig_MC_DoubleMuNoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_IsoMu27_MET90_accept;   //!
   TBranch        *b_trig_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_accept;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_accept;   //!
   TBranch        *b_trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_accept;   //!
   TBranch        *b_trig_HLT_PFMET100_PFMHT100_IDTight_PFHT60_accept;   //!
   TBranch        *b_trig_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_accept;   //!
   TBranch        *b_trig_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_accept;   //!
   TBranch        *b_trig_HLT_Mu18_Mu9_SameSign_accept;   //!
   TBranch        *b_trig_HLT_Mu18_Mu9_SameSign_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu18_Mu9_accept;   //!
   TBranch        *b_trig_HLT_Mu18_Mu9_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu20_Mu10_SameSign_accept;   //!
   TBranch        *b_trig_HLT_Mu20_Mu10_SameSign_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu20_Mu10_accept;   //!
   TBranch        *b_trig_HLT_Mu20_Mu10_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu23_Mu12_SameSign_accept;   //!
   TBranch        *b_trig_HLT_Mu23_Mu12_SameSign_DZ_accept;   //!
   TBranch        *b_trig_HLT_Mu23_Mu12_accept;   //!
   TBranch        *b_trig_HLT_Mu23_Mu12_DZ_accept;   //!
   TBranch        *b_trig_HLT_DoubleMu3_DCA_PFMET50_PFMHT60_accept;   //!
   TBranch        *b_trig_HLT_TrkMu6NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_TrkMu16NoFiltersNoVtx_accept;   //!
   TBranch        *b_trig_HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_accept;   //!

   MyAnalysis(TTree *tree=0);
   virtual ~MyAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString, TString, TString, TString, TString, float,float,float);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyAnalysis_cxx
MyAnalysis::MyAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/pnfs/iihe/cms/store/user/schenara/SYS_2018/TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8/crab_TTTo2L2Nu_hdampUP_2018/200211_122927/0000/outfile_5_2016_71.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/pnfs/iihe/cms/store/user/schenara/SYS_2018/TTTo2L2Nu_hdampUP_TuneCP5_13TeV-powheg-pythia8/crab_TTTo2L2Nu_hdampUP_2018/200211_122927/0000/outfile_5_2016_71.root");
      }
      f->GetObject("IIHEAnalysis",tree);

   }
   Init(tree);
}

MyAnalysis::~MyAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   LHE_Pt = 0;
   LHE_Eta = 0;
   LHE_Phi = 0;
   LHE_E = 0;
   LHE_pdgid = 0;
   LHE_status = 0;
   LHE_weight_sys = 0;
   mc_DJRValues = 0;
   mc_index = 0;
   mc_pdgId = 0;
   mc_charge = 0;
   mc_status = 0;
   mc_status_flags = 0;
   mc_status_tau_flags = 0;
   mc_tau_charge = 0;
   mc_tau_pdgId = 0;
   mc_tau_decay = 0;
   mc_tau_had_status = 0;
   mc_tau_had_charge = 0;
   mc_tau_had_pdgId = 0;
   mc_mass = 0;
   mc_px = 0;
   mc_py = 0;
   mc_pz = 0;
   mc_pt = 0;
   mc_eta = 0;
   mc_phi = 0;
   mc_energy = 0;
   mc_tau_pt = 0;
   mc_tau_eta = 0;
   mc_tau_phi = 0;
   mc_tau_energy = 0;
   mc_tau_had_pt = 0;
   mc_tau_had_eta = 0;
   mc_tau_had_phi = 0;
   mc_tau_had_energy = 0;
   mc_numberOfDaughters = 0;
   mc_numberOfMothers = 0;
   mc_mother_index = 0;
   mc_mother_pdgId = 0;
   mc_mother_px = 0;
   mc_mother_py = 0;
   mc_mother_pz = 0;
   mc_mother_pt = 0;
   mc_mother_eta = 0;
   mc_mother_phi = 0;
   mc_mother_energy = 0;
   mc_mother_mass = 0;
   genjet_pt = 0;
   genjet_eta = 0;
   genjet_phi = 0;
   genjet_energy = 0;
   genjetAK8_pt = 0;
   genjetAK8_eta = 0;
   genjetAK8_phi = 0;
   genjetAK8_energy = 0;
   gen_weight_sys = 0;
   pv_x = 0;
   pv_y = 0;
   pv_z = 0;
   pv_ndof = 0;
   pv_normalizedChi2 = 0;
   pv_isValid = 0;
   pv_isFake = 0;
   gsf_classification = 0;
   gsf_ecalTrkEnergyPostCorr = 0;
   gsf_ecalEnergyPostCorr = 0;
   gsf_energy = 0;
   gsf_p = 0;
   gsf_pt = 0;
   gsf_et = 0;
   gsf_scE1x5 = 0;
   gsf_scE5x5 = 0;
   gsf_scE2x5Max = 0;
   gsf_full5x5_e5x5 = 0;
   gsf_full5x5_e1x5 = 0;
   gsf_full5x5_e2x5Max = 0;
   gsf_full5x5_sigmaIetaIeta = 0;
   gsf_full5x5_hcalOverEcal = 0;
   gsf_eta = 0;
   gsf_phi = 0;
   gsf_theta = 0;
   gsf_px = 0;
   gsf_py = 0;
   gsf_pz = 0;
   gsf_caloEnergy = 0;
   gsf_deltaEtaSuperClusterTrackAtVtx = 0;
   gsf_deltaPhiSuperClusterTrackAtVtx = 0;
   gsf_hadronicOverEm = 0;
   gsf_hcalDepth1OverEcal = 0;
   gsf_hcalDepth2OverEcal = 0;
   gsf_dr03TkSumPt = 0;
   gsf_heepTrkPtIso = 0;
   gsf_dr03EcalRecHitSumEt = 0;
   gsf_dr03HcalDepth1TowerSumEt = 0;
   gsf_dr03HcalDepth2TowerSumEt = 0;
   gsf_sigmaIetaIeta = 0;
   gsf_charge = 0;
   gsf_ecaldrivenSeed = 0;
   gsf_trackerdrivenSeed = 0;
   gsf_isEB = 0;
   gsf_isEE = 0;
   gsf_isPF = 0;
   gsf_passConversionVeto = 0;
   gsf_VID_cutBasedElectronID_Fall17_94X_V2_loose = 0;
   gsf_VID_cutBasedElectronID_Fall17_94X_V2_medium = 0;
   gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight = 0;
   gsf_VID_cutBasedElectronID_Fall17_94X_V2_veto = 0;
   gsf_VID_heepElectronID_HEEPV70 = 0;
   gsf_VID_mvaEleID_Fall17_iso_V2_wp80 = 0;
   gsf_VID_mvaEleID_Fall17_iso_V2_wp90 = 0;
   gsf_VID_mvaEleID_Fall17_iso_V2_wpLoose = 0;
   gsf_VID_mvaEleID_Fall17_noIso_V2_wp80 = 0;
   gsf_VID_mvaEleID_Fall17_noIso_V2_wp90 = 0;
   gsf_VID_mvaEleID_Fall17_noIso_V2_wpLoose = 0;
   gsf_deltaEtaSeedClusterTrackAtCalo = 0;
   gsf_deltaPhiSeedClusterTrackAtCalo = 0;
   gsf_ecalEnergy = 0;
   gsf_eSuperClusterOverP = 0;
   gsf_dxy = 0;
   gsf_dxy_beamSpot = 0;
   gsf_dxy_firstPVtx = 0;
   gsf_dxyError = 0;
   gsf_dz = 0;
   gsf_dz_beamSpot = 0;
   gsf_dz_firstPVtx = 0;
   gsf_dzError = 0;
   gsf_vz = 0;
   gsf_numberOfValidHits = 0;
   gsf_nLostInnerHits = 0;
   gsf_nLostOuterHits = 0;
   gsf_convFlags = 0;
   gsf_convDist = 0;
   gsf_convDcot = 0;
   gsf_convRadius = 0;
   gsf_fBrem = 0;
   gsf_e1x5 = 0;
   gsf_e2x5Max = 0;
   gsf_e5x5 = 0;
   gsf_r9 = 0;
   gsf_deltaEtaSeedClusterTrackAtVtx = 0;
   gsf_effArea = 0;
   gsf_sumChargedHadronPt = 0;
   gsf_sumNeutralHadronEt = 0;
   gsf_sumPhotonEt = 0;
   gsf_sumPUPt = 0;
   gsf_ooEmooP = 0;
   gsf_hitsinfo = 0;
   gsf_pixelMatch_dPhi1 = 0;
   gsf_pixelMatch_dPhi2 = 0;
   gsf_pixelMatch_dRz1 = 0;
   gsf_pixelMatch_dRz2 = 0;
   gsf_pixelMatch_subDetector1 = 0;
   gsf_pixelMatch_subDetector2 = 0;
   gsf_mc_bestDR = 0;
   gsf_mc_index = 0;
   gsf_mc_ERatio = 0;
   gsf_sc_energy = 0;
   gsf_sc_seed_eta = 0;
   gsf_sc_eta = 0;
   gsf_sc_etacorr = 0;
   gsf_sc_theta = 0;
   gsf_sc_thetacorr = 0;
   gsf_sc_et = 0;
   gsf_sc_phi = 0;
   gsf_sc_px = 0;
   gsf_sc_py = 0;
   gsf_sc_pz = 0;
   gsf_sc_x = 0;
   gsf_sc_y = 0;
   gsf_sc_z = 0;
   gsf_sc_phiWidth = 0;
   gsf_sc_etaWidth = 0;
   gsf_sc_seed_rawId = 0;
   gsf_sc_seed_ieta = 0;
   gsf_sc_seed_iphi = 0;
   gsf_sc_seed_kHasSwitchToGain6 = 0;
   gsf_sc_seed_kHasSwitchToGain1 = 0;
   EBHits_rawId = 0;
   EBHits_iRechit = 0;
   EBHits_ieta = 0;
   EBHits_iphi = 0;
   EBHits_RecoFlag = 0;
   EBHits_kSaturated = 0;
   EBHits_kLeadingEdgeRecovered = 0;
   EBHits_kNeighboursRecovered = 0;
   EBHits_kWeird = 0;
   EBHits_energy = 0;
   EEHits_rawId = 0;
   EEHits_iRechit = 0;
   EEHits_ieta = 0;
   EEHits_iphi = 0;
   EEHits_RecoFlag = 0;
   EEHits_kSaturated = 0;
   EEHits_kLeadingEdgeRecovered = 0;
   EEHits_kNeighboursRecovered = 0;
   EEHits_kWeird = 0;
   EEHits_energy = 0;
   mu_gt_qoverp = 0;
   mu_gt_charge = 0;
   mu_gt_pt = 0;
   mu_gt_eta = 0;
   mu_gt_phi = 0;
   mu_gt_p = 0;
   mu_gt_px = 0;
   mu_gt_py = 0;
   mu_gt_pz = 0;
   mu_gt_theta = 0;
   mu_gt_lambda = 0;
   mu_gt_d0 = 0;
   mu_gt_dz = 0;
   mu_gt_dz_beamspot = 0;
   mu_gt_dz_firstPVtx = 0;
   mu_gt_dxy = 0;
   mu_gt_dxy_beamspot = 0;
   mu_gt_dxy_firstPVtx = 0;
   mu_gt_dsz = 0;
   mu_gt_vx = 0;
   mu_gt_vy = 0;
   mu_gt_vz = 0;
   mu_gt_qoverpError = 0;
   mu_gt_ptError = 0;
   mu_gt_thetaError = 0;
   mu_gt_lambdaError = 0;
   mu_gt_phiError = 0;
   mu_gt_dxyError = 0;
   mu_gt_d0Error = 0;
   mu_gt_dszError = 0;
   mu_gt_dzError = 0;
   mu_gt_etaError = 0;
   mu_gt_chi2 = 0;
   mu_gt_ndof = 0;
   mu_gt_normalizedChi2 = 0;
   mu_ot_qoverp = 0;
   mu_ot_charge = 0;
   mu_ot_pt = 0;
   mu_ot_eta = 0;
   mu_ot_phi = 0;
   mu_ot_p = 0;
   mu_ot_px = 0;
   mu_ot_py = 0;
   mu_ot_pz = 0;
   mu_ot_theta = 0;
   mu_ot_lambda = 0;
   mu_ot_d0 = 0;
   mu_ot_dz = 0;
   mu_ot_dz_beamspot = 0;
   mu_ot_dz_firstPVtx = 0;
   mu_ot_dxy = 0;
   mu_ot_dxy_beamspot = 0;
   mu_ot_dxy_firstPVtx = 0;
   mu_ot_dsz = 0;
   mu_ot_vx = 0;
   mu_ot_vy = 0;
   mu_ot_vz = 0;
   mu_ot_qoverpError = 0;
   mu_ot_ptError = 0;
   mu_ot_thetaError = 0;
   mu_ot_lambdaError = 0;
   mu_ot_phiError = 0;
   mu_ot_dxyError = 0;
   mu_ot_d0Error = 0;
   mu_ot_dszError = 0;
   mu_ot_dzError = 0;
   mu_ot_etaError = 0;
   mu_ot_chi2 = 0;
   mu_ot_ndof = 0;
   mu_ot_normalizedChi2 = 0;
   mu_it_qoverp = 0;
   mu_it_charge = 0;
   mu_it_pt = 0;
   mu_it_eta = 0;
   mu_it_phi = 0;
   mu_it_p = 0;
   mu_it_px = 0;
   mu_it_py = 0;
   mu_it_pz = 0;
   mu_it_theta = 0;
   mu_it_lambda = 0;
   mu_it_d0 = 0;
   mu_it_dz = 0;
   mu_it_dz_beamspot = 0;
   mu_it_dz_firstPVtx = 0;
   mu_it_dxy = 0;
   mu_it_dxy_beamspot = 0;
   mu_it_dxy_firstPVtx = 0;
   mu_it_dsz = 0;
   mu_it_vx = 0;
   mu_it_vy = 0;
   mu_it_vz = 0;
   mu_it_qoverpError = 0;
   mu_it_ptError = 0;
   mu_it_thetaError = 0;
   mu_it_lambdaError = 0;
   mu_it_phiError = 0;
   mu_it_dxyError = 0;
   mu_it_d0Error = 0;
   mu_it_dszError = 0;
   mu_it_dzError = 0;
   mu_it_etaError = 0;
   mu_it_chi2 = 0;
   mu_it_ndof = 0;
   mu_it_normalizedChi2 = 0;
   mu_ibt_qoverp = 0;
   mu_ibt_charge = 0;
   mu_ibt_pt = 0;
   mu_ibt_eta = 0;
   mu_ibt_phi = 0;
   mu_ibt_p = 0;
   mu_ibt_px = 0;
   mu_ibt_py = 0;
   mu_ibt_pz = 0;
   mu_ibt_theta = 0;
   mu_ibt_lambda = 0;
   mu_ibt_d0 = 0;
   mu_ibt_dz = 0;
   mu_ibt_dz_beamspot = 0;
   mu_ibt_dz_firstPVtx = 0;
   mu_ibt_dxy = 0;
   mu_ibt_dxy_beamspot = 0;
   mu_ibt_dxy_firstPVtx = 0;
   mu_ibt_dsz = 0;
   mu_ibt_vx = 0;
   mu_ibt_vy = 0;
   mu_ibt_vz = 0;
   mu_ibt_qoverpError = 0;
   mu_ibt_ptError = 0;
   mu_ibt_thetaError = 0;
   mu_ibt_lambdaError = 0;
   mu_ibt_phiError = 0;
   mu_ibt_dxyError = 0;
   mu_ibt_d0Error = 0;
   mu_ibt_dszError = 0;
   mu_ibt_dzError = 0;
   mu_ibt_etaError = 0;
   mu_ibt_chi2 = 0;
   mu_ibt_ndof = 0;
   mu_ibt_normalizedChi2 = 0;
   mu_isGlobalMuon = 0;
   mu_isStandAloneMuon = 0;
   mu_isTrackerMuon = 0;
   mu_isPFMuon = 0;
   mu_isPFIsolationValid = 0;
   mu_isGoodMuonTMLastStationLoose = 0;
   mu_isGoodMuonTMLastStationTight = 0;
   mu_isGoodMuonTM2DCompatibilityLoose = 0;
   mu_isGoodMuonTM2DCompatibilityTight = 0;
   mu_isGoodMuonTMOneStationLoose = 0;
   mu_isGoodMuonTMOneStationTight = 0;
   mu_isGoodMuonTMLastStationOptimizedLowPtLoose = 0;
   mu_isGoodMuonTMLastStationOptimizedLowPtTight = 0;
   mu_isTightMuon = 0;
   mu_isMediumMuon = 0;
   mu_isLooseMuon = 0;
   mu_isSoftMuon = 0;
   mu_isHighPtMuonOld = 0;
   mu_isHighPtMuon = 0;
   mu_isTrackerHighPtMuon = 0;
   mu_CutBasedIdLoose = 0;
   mu_CutBasedIdMedium = 0;
   mu_CutBasedIdMediumPrompt = 0;
   mu_CutBasedIdTight = 0;
   mu_CutBasedIdGlobalHighPt = 0;
   mu_CutBasedIdTrkHighPt = 0;
   mu_PFIsoVeryLoose = 0;
   mu_PFIsoLoose = 0;
   mu_PFIsoMedium = 0;
   mu_PFIsoTight = 0;
   mu_PFIsoVeryTight = 0;
   mu_PFIsoVeryVeryTight = 0;
   mu_TkIsoLoose = 0;
   mu_TkIsoTight = 0;
   mu_SoftCutBasedId = 0;
   mu_SoftMvaId = 0;
   mu_MvaLoose = 0;
   mu_MvaMedium = 0;
   mu_MvaTight = 0;
   mu_MiniIsoLoose = 0;
   mu_MiniIsoMedium = 0;
   mu_MiniIsoTight = 0;
   mu_MiniIsoVeryTight = 0;
   mu_TriggerIdLoose = 0;
   mu_InTimeMuon = 0;
   mu_MultiIsoLoose = 0;
   mu_MultiIsoMedium = 0;
   mu_numberOfMatchedStations = 0;
   mu_numberOfValidPixelHits = 0;
   mu_trackerLayersWithMeasurement = 0;
   mu_numberOfValidMuonHits = 0;
   mu_pixelLayersWithMeasurement = 0;
   mu_innerTrack_validFraction = 0;
   mu_combinedQuality_trkKink = 0;
   mu_combinedQuality_chi2LocalPosition = 0;
   mu_segmentCompatibility = 0;
   mu_dB = 0;
   mu_pt_default = 0;
   mu_isolationR03_sumPt = 0;
   mu_isolationR03_trackerVetoPt = 0;
   mu_isolationR03_emEt = 0;
   mu_isolationR03_emVetoEt = 0;
   mu_isolationR03_hadEt = 0;
   mu_isolationR03_hadVetoEt = 0;
   mu_isolationR05_sumPt = 0;
   mu_isolationR05_trackerVetoPt = 0;
   mu_isolationR05_emEt = 0;
   mu_isolationR05_emVetoEt = 0;
   mu_isolationR05_hadEt = 0;
   mu_isolationR05_hadVetoEt = 0;
   mu_pfIsolationR03_sumChargedHadronPt = 0;
   mu_pfIsolationR03_sumNeutralHadronEt = 0;
   mu_pfIsolationR03_sumChargedParticlePt = 0;
   mu_pfIsolationR03_sumPhotonEt = 0;
   mu_pfIsolationR03_sumNeutralHadronEtHighThreshold = 0;
   mu_pfIsolationR03_sumPhotonEtHighThreshold = 0;
   mu_pfIsolationR03_sumPUPt = 0;
   mu_pfIsolationR04_sumChargedHadronPt = 0;
   mu_pfIsolationR04_sumNeutralHadronEt = 0;
   mu_pfIsolationR04_sumChargedParticlePt = 0;
   mu_pfIsolationR04_sumPhotonEt = 0;
   mu_pfIsolationR04_sumNeutralHadronEtHighThreshold = 0;
   mu_pfIsolationR04_sumPhotonEtHighThreshold = 0;
   mu_pfIsolationR04_sumPUPt = 0;
   mu_pfIsoDbCorrected03 = 0;
   mu_pfIsoDbCorrected04 = 0;
   mu_isoTrackerBased03 = 0;
   mu_mc_bestDR = 0;
   mu_mc_index = 0;
   mu_mc_ERatio = 0;
   jet_px = 0;
   jet_py = 0;
   jet_pz = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_theta = 0;
   jet_phi = 0;
   jet_energy = 0;
   jet_mass = 0;
   jet_chargedEmEnergyFraction = 0;
   jet_neutralHadronEnergyFraction = 0;
   jet_neutralEmEnergyFraction = 0;
   jet_chargedHadronEnergyFraction = 0;
   jet_muonEnergyFraction = 0;
   jet_chargedMultiplicity = 0;
   jet_neutralMultiplicity = 0;
   jet_partonFlavour = 0;
   jet_hadronFlavour = 0;
   jet_CSVv2 = 0;
   jet_CvsL = 0;
   jet_CvsB = 0;
   jet_MVA2BJets = 0;
   jet_CvsB_DeepJet_charm_tagger = 0;
   jet_CvsL_DeepJet_charm_tagger = 0;
   jet_CvsB_DeepCSV_charm_tagger = 0;
   jet_CvsL_DeepCSV_charm_tagger = 0;
   jet_DeepJet = 0;
   jet_DeepCSV = 0;
   jet_isJetIDLoose_2016 = 0;
   jet_isJetIDTight_2016 = 0;
   jet_isJetIDTightLepVeto_2016 = 0;
   jet_isJetID_2017 = 0;
   jet_isJetIDLepVeto_2017 = 0;
   jet_isJetID_2018 = 0;
   jet_isJetIDLepVeto_2018 = 0;
   jet_Smeared_pt = 0;
   jet_SmearedJetResUp_pt = 0;
   jet_SmearedJetResDown_pt = 0;
   jet_SmearedJetEnUp_pt = 0;
   jet_SmearedJetEnDown_pt = 0;
   MET_Type1Unc_Px = 0;
   MET_Type1Unc_Py = 0;
   MET_Type1Unc_Pt = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev_event", &ev_event, &b_ev_event);
   fChain->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
   fChain->SetBranchAddress("ev_luminosityBlock", &ev_luminosityBlock, &b_ev_luminosityBlock);
   fChain->SetBranchAddress("ev_time", &ev_time, &b_ev_time);
   fChain->SetBranchAddress("ev_time_unixTime", &ev_time_unixTime, &b_ev_time_unixTime);
   fChain->SetBranchAddress("ev_time_microsecondOffset", &ev_time_microsecondOffset, &b_ev_time_microsecondOffset);
   fChain->SetBranchAddress("ev_fixedGridRhoAll", &ev_fixedGridRhoAll, &b_ev_fixedGridRhoAll);
   fChain->SetBranchAddress("ev_fixedGridRhoFastjetAll", &ev_fixedGridRhoFastjetAll, &b_ev_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("ev_fixedGridRhoFastjetAllCalo", &ev_fixedGridRhoFastjetAllCalo, &b_ev_fixedGridRhoFastjetAllCalo);
   fChain->SetBranchAddress("ev_fixedGridRhoFastjetCentralCalo", &ev_fixedGridRhoFastjetCentralCalo, &b_ev_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("ev_fixedGridRhoFastjetCentralChargedPileUp", &ev_fixedGridRhoFastjetCentralChargedPileUp, &b_ev_fixedGridRhoFastjetCentralChargedPileUp);
   fChain->SetBranchAddress("ev_fixedGridRhoFastjetCentralNeutral", &ev_fixedGridRhoFastjetCentralNeutral, &b_ev_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("ev_prefiringweight", &ev_prefiringweight, &b_ev_prefiringweight);
   fChain->SetBranchAddress("ev_prefiringweightup", &ev_prefiringweightup, &b_ev_prefiringweightup);
   fChain->SetBranchAddress("ev_prefiringweightdown", &ev_prefiringweightdown, &b_ev_prefiringweightdown);
   fChain->SetBranchAddress("LHE_Pt", &LHE_Pt, &b_LHE_Pt);
   fChain->SetBranchAddress("LHE_Eta", &LHE_Eta, &b_LHE_Eta);
   fChain->SetBranchAddress("LHE_Phi", &LHE_Phi, &b_LHE_Phi);
   fChain->SetBranchAddress("LHE_E", &LHE_E, &b_LHE_E);
   fChain->SetBranchAddress("LHE_pdgid", &LHE_pdgid, &b_LHE_pdgid);
   fChain->SetBranchAddress("LHE_status", &LHE_status, &b_LHE_status);
   fChain->SetBranchAddress("LHE_weight_nominal", &LHE_weight_nominal, &b_LHE_weight_nominal);
   fChain->SetBranchAddress("LHE_weight_sys", &LHE_weight_sys, &b_LHE_weight_sys);
   fChain->SetBranchAddress("mc_n", &mc_n, &b_mc_n);
   fChain->SetBranchAddress("mc_nMEPartons", &mc_nMEPartons, &b_mc_nMEPartons);
   fChain->SetBranchAddress("mc_nMEPartonsFiltered", &mc_nMEPartonsFiltered, &b_mc_nMEPartonsFiltered);
   fChain->SetBranchAddress("mc_DJRValues", &mc_DJRValues, &b_mc_DJRValues);
   fChain->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight);
   fChain->SetBranchAddress("mc_w_sign", &mc_w_sign, &b_mc_w_sign);
   fChain->SetBranchAddress("mc_id_first", &mc_id_first, &b_mc_id_first);
   fChain->SetBranchAddress("mc_id_second", &mc_id_second, &b_mc_id_second);
   fChain->SetBranchAddress("mc_x_first", &mc_x_first, &b_mc_x_first);
   fChain->SetBranchAddress("mc_x_second", &mc_x_second, &b_mc_x_second);
   fChain->SetBranchAddress("mc_xPDF_first", &mc_xPDF_first, &b_mc_xPDF_first);
   fChain->SetBranchAddress("mc_xPDF_second", &mc_xPDF_second, &b_mc_xPDF_second);
   fChain->SetBranchAddress("mc_scalePDF", &mc_scalePDF, &b_mc_scalePDF);
   fChain->SetBranchAddress("mc_index", &mc_index, &b_mc_index);
   fChain->SetBranchAddress("mc_pdgId", &mc_pdgId, &b_mc_pdgId);
   fChain->SetBranchAddress("mc_charge", &mc_charge, &b_mc_charge);
   fChain->SetBranchAddress("mc_status", &mc_status, &b_mc_status);
   fChain->SetBranchAddress("mc_status_flags", &mc_status_flags, &b_mc_status_flags);
   fChain->SetBranchAddress("mc_status_tau_flags", &mc_status_tau_flags, &b_mc_status_tau_flags);
   fChain->SetBranchAddress("mc_tau_charge", &mc_tau_charge, &b_mc_tau_charge);
   fChain->SetBranchAddress("mc_tau_pdgId", &mc_tau_pdgId, &b_mc_tau_pdgId);
   fChain->SetBranchAddress("mc_tau_decay", &mc_tau_decay, &b_mc_tau_decay);
   fChain->SetBranchAddress("mc_tau_had_status", &mc_tau_had_status, &b_mc_tau_had_status);
   fChain->SetBranchAddress("mc_tau_had_charge", &mc_tau_had_charge, &b_mc_tau_had_charge);
   fChain->SetBranchAddress("mc_tau_had_pdgId", &mc_tau_had_pdgId, &b_mc_tau_had_pdgId);
   fChain->SetBranchAddress("mc_mass", &mc_mass, &b_mc_mass);
   fChain->SetBranchAddress("mc_px", &mc_px, &b_mc_px);
   fChain->SetBranchAddress("mc_py", &mc_py, &b_mc_py);
   fChain->SetBranchAddress("mc_pz", &mc_pz, &b_mc_pz);
   fChain->SetBranchAddress("mc_pt", &mc_pt, &b_mc_pt);
   fChain->SetBranchAddress("mc_eta", &mc_eta, &b_mc_eta);
   fChain->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi);
   fChain->SetBranchAddress("mc_energy", &mc_energy, &b_mc_energy);
   fChain->SetBranchAddress("mc_tau_pt", &mc_tau_pt, &b_mc_tau_pt);
   fChain->SetBranchAddress("mc_tau_eta", &mc_tau_eta, &b_mc_tau_eta);
   fChain->SetBranchAddress("mc_tau_phi", &mc_tau_phi, &b_mc_tau_phi);
   fChain->SetBranchAddress("mc_tau_energy", &mc_tau_energy, &b_mc_tau_energy);
   fChain->SetBranchAddress("mc_tau_had_pt", &mc_tau_had_pt, &b_mc_tau_had_pt);
   fChain->SetBranchAddress("mc_tau_had_eta", &mc_tau_had_eta, &b_mc_tau_had_eta);
   fChain->SetBranchAddress("mc_tau_had_phi", &mc_tau_had_phi, &b_mc_tau_had_phi);
   fChain->SetBranchAddress("mc_tau_had_energy", &mc_tau_had_energy, &b_mc_tau_had_energy);
   fChain->SetBranchAddress("mc_numberOfDaughters", &mc_numberOfDaughters, &b_mc_numberOfDaughters);
   fChain->SetBranchAddress("mc_numberOfMothers", &mc_numberOfMothers, &b_mc_numberOfMothers);
   fChain->SetBranchAddress("mc_mother_index", &mc_mother_index, &b_mc_mother_index);
   fChain->SetBranchAddress("mc_mother_pdgId", &mc_mother_pdgId, &b_mc_mother_pdgId);
   fChain->SetBranchAddress("mc_mother_px", &mc_mother_px, &b_mc_mother_px);
   fChain->SetBranchAddress("mc_mother_py", &mc_mother_py, &b_mc_mother_py);
   fChain->SetBranchAddress("mc_mother_pz", &mc_mother_pz, &b_mc_mother_pz);
   fChain->SetBranchAddress("mc_mother_pt", &mc_mother_pt, &b_mc_mother_pt);
   fChain->SetBranchAddress("mc_mother_eta", &mc_mother_eta, &b_mc_mother_eta);
   fChain->SetBranchAddress("mc_mother_phi", &mc_mother_phi, &b_mc_mother_phi);
   fChain->SetBranchAddress("mc_mother_energy", &mc_mother_energy, &b_mc_mother_energy);
   fChain->SetBranchAddress("mc_mother_mass", &mc_mother_mass, &b_mc_mother_mass);
   fChain->SetBranchAddress("mc_trueNumInteractions", &mc_trueNumInteractions, &b_mc_trueNumInteractions);
   fChain->SetBranchAddress("mc_PU_NumInteractions", &mc_PU_NumInteractions, &b_mc_PU_NumInteractions);
   fChain->SetBranchAddress("genjet_pt", &genjet_pt, &b_genjet_pt);
   fChain->SetBranchAddress("genjet_eta", &genjet_eta, &b_genjet_eta);
   fChain->SetBranchAddress("genjet_phi", &genjet_phi, &b_genjet_phi);
   fChain->SetBranchAddress("genjet_energy", &genjet_energy, &b_genjet_energy);
   fChain->SetBranchAddress("genjetAK8_pt", &genjetAK8_pt, &b_genjetAK8_pt);
   fChain->SetBranchAddress("genjetAK8_eta", &genjetAK8_eta, &b_genjetAK8_eta);
   fChain->SetBranchAddress("genjetAK8_phi", &genjetAK8_phi, &b_genjetAK8_phi);
   fChain->SetBranchAddress("genjetAK8_energy", &genjetAK8_energy, &b_genjetAK8_energy);
   fChain->SetBranchAddress("gen_weight_sys", &gen_weight_sys, &b_gen_weight_sys);
   fChain->SetBranchAddress("pv_n", &pv_n, &b_pv_n);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("pv_normalizedChi2", &pv_normalizedChi2, &b_pv_normalizedChi2);
   fChain->SetBranchAddress("pv_isValid", &pv_isValid, &b_pv_isValid);
   fChain->SetBranchAddress("pv_isFake", &pv_isFake, &b_pv_isFake);
   fChain->SetBranchAddress("gsf_n", &gsf_n, &b_gsf_n);
   fChain->SetBranchAddress("gsf_classification", &gsf_classification, &b_gsf_classification);
   fChain->SetBranchAddress("gsf_ecalTrkEnergyPostCorr", &gsf_ecalTrkEnergyPostCorr, &b_gsf_ecalTrkEnergyPostCorr);
   fChain->SetBranchAddress("gsf_ecalEnergyPostCorr", &gsf_ecalEnergyPostCorr, &b_gsf_ecalEnergyPostCorr);
   fChain->SetBranchAddress("gsf_energy", &gsf_energy, &b_gsf_energy);
   fChain->SetBranchAddress("gsf_p", &gsf_p, &b_gsf_p);
   fChain->SetBranchAddress("gsf_pt", &gsf_pt, &b_gsf_pt);
   fChain->SetBranchAddress("gsf_et", &gsf_et, &b_gsf_et);
   fChain->SetBranchAddress("gsf_scE1x5", &gsf_scE1x5, &b_gsf_scE1x5);
   fChain->SetBranchAddress("gsf_scE5x5", &gsf_scE5x5, &b_gsf_scE5x5);
   fChain->SetBranchAddress("gsf_scE2x5Max", &gsf_scE2x5Max, &b_gsf_scE2x5Max);
   fChain->SetBranchAddress("gsf_full5x5_e5x5", &gsf_full5x5_e5x5, &b_gsf_full5x5_e5x5);
   fChain->SetBranchAddress("gsf_full5x5_e1x5", &gsf_full5x5_e1x5, &b_gsf_full5x5_e1x5);
   fChain->SetBranchAddress("gsf_full5x5_e2x5Max", &gsf_full5x5_e2x5Max, &b_gsf_full5x5_e2x5Max);
   fChain->SetBranchAddress("gsf_full5x5_sigmaIetaIeta", &gsf_full5x5_sigmaIetaIeta, &b_gsf_full5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("gsf_full5x5_hcalOverEcal", &gsf_full5x5_hcalOverEcal, &b_gsf_full5x5_hcalOverEcal);
   fChain->SetBranchAddress("gsf_eta", &gsf_eta, &b_gsf_eta);
   fChain->SetBranchAddress("gsf_phi", &gsf_phi, &b_gsf_phi);
   fChain->SetBranchAddress("gsf_theta", &gsf_theta, &b_gsf_theta);
   fChain->SetBranchAddress("gsf_px", &gsf_px, &b_gsf_px);
   fChain->SetBranchAddress("gsf_py", &gsf_py, &b_gsf_py);
   fChain->SetBranchAddress("gsf_pz", &gsf_pz, &b_gsf_pz);
   fChain->SetBranchAddress("gsf_caloEnergy", &gsf_caloEnergy, &b_gsf_caloEnergy);
   fChain->SetBranchAddress("gsf_deltaEtaSuperClusterTrackAtVtx", &gsf_deltaEtaSuperClusterTrackAtVtx, &b_gsf_deltaEtaSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("gsf_deltaPhiSuperClusterTrackAtVtx", &gsf_deltaPhiSuperClusterTrackAtVtx, &b_gsf_deltaPhiSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("gsf_hadronicOverEm", &gsf_hadronicOverEm, &b_gsf_hadronicOverEm);
   fChain->SetBranchAddress("gsf_hcalDepth1OverEcal", &gsf_hcalDepth1OverEcal, &b_gsf_hcalDepth1OverEcal);
   fChain->SetBranchAddress("gsf_hcalDepth2OverEcal", &gsf_hcalDepth2OverEcal, &b_gsf_hcalDepth2OverEcal);
   fChain->SetBranchAddress("gsf_dr03TkSumPt", &gsf_dr03TkSumPt, &b_gsf_dr03TkSumPt);
   fChain->SetBranchAddress("gsf_heepTrkPtIso", &gsf_heepTrkPtIso, &b_gsf_heepTrkPtIso);
   fChain->SetBranchAddress("gsf_dr03EcalRecHitSumEt", &gsf_dr03EcalRecHitSumEt, &b_gsf_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("gsf_dr03HcalDepth1TowerSumEt", &gsf_dr03HcalDepth1TowerSumEt, &b_gsf_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("gsf_dr03HcalDepth2TowerSumEt", &gsf_dr03HcalDepth2TowerSumEt, &b_gsf_dr03HcalDepth2TowerSumEt);
   fChain->SetBranchAddress("gsf_sigmaIetaIeta", &gsf_sigmaIetaIeta, &b_gsf_sigmaIetaIeta);
   fChain->SetBranchAddress("gsf_charge", &gsf_charge, &b_gsf_charge);
   fChain->SetBranchAddress("gsf_ecaldrivenSeed", &gsf_ecaldrivenSeed, &b_gsf_ecaldrivenSeed);
   fChain->SetBranchAddress("gsf_trackerdrivenSeed", &gsf_trackerdrivenSeed, &b_gsf_trackerdrivenSeed);
   fChain->SetBranchAddress("gsf_isEB", &gsf_isEB, &b_gsf_isEB);
   fChain->SetBranchAddress("gsf_isEE", &gsf_isEE, &b_gsf_isEE);
   fChain->SetBranchAddress("gsf_isPF", &gsf_isPF, &b_gsf_isPF);
   fChain->SetBranchAddress("gsf_passConversionVeto", &gsf_passConversionVeto, &b_gsf_passConversionVeto);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V2_loose", &gsf_VID_cutBasedElectronID_Fall17_94X_V2_loose, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V2_loose);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V2_medium", &gsf_VID_cutBasedElectronID_Fall17_94X_V2_medium, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V2_medium);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight", &gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight);
   fChain->SetBranchAddress("gsf_VID_cutBasedElectronID_Fall17_94X_V2_veto", &gsf_VID_cutBasedElectronID_Fall17_94X_V2_veto, &b_gsf_VID_cutBasedElectronID_Fall17_94X_V2_veto);
   fChain->SetBranchAddress("gsf_VID_heepElectronID_HEEPV70", &gsf_VID_heepElectronID_HEEPV70, &b_gsf_VID_heepElectronID_HEEPV70);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_iso_V2_wp80", &gsf_VID_mvaEleID_Fall17_iso_V2_wp80, &b_gsf_VID_mvaEleID_Fall17_iso_V2_wp80);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_iso_V2_wp90", &gsf_VID_mvaEleID_Fall17_iso_V2_wp90, &b_gsf_VID_mvaEleID_Fall17_iso_V2_wp90);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_iso_V2_wpLoose", &gsf_VID_mvaEleID_Fall17_iso_V2_wpLoose, &b_gsf_VID_mvaEleID_Fall17_iso_V2_wpLoose);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_noIso_V2_wp80", &gsf_VID_mvaEleID_Fall17_noIso_V2_wp80, &b_gsf_VID_mvaEleID_Fall17_noIso_V2_wp80);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_noIso_V2_wp90", &gsf_VID_mvaEleID_Fall17_noIso_V2_wp90, &b_gsf_VID_mvaEleID_Fall17_noIso_V2_wp90);
   fChain->SetBranchAddress("gsf_VID_mvaEleID_Fall17_noIso_V2_wpLoose", &gsf_VID_mvaEleID_Fall17_noIso_V2_wpLoose, &b_gsf_VID_mvaEleID_Fall17_noIso_V2_wpLoose);
   fChain->SetBranchAddress("gsf_deltaEtaSeedClusterTrackAtCalo", &gsf_deltaEtaSeedClusterTrackAtCalo, &b_gsf_deltaEtaSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("gsf_deltaPhiSeedClusterTrackAtCalo", &gsf_deltaPhiSeedClusterTrackAtCalo, &b_gsf_deltaPhiSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("gsf_ecalEnergy", &gsf_ecalEnergy, &b_gsf_ecalEnergy);
   fChain->SetBranchAddress("gsf_eSuperClusterOverP", &gsf_eSuperClusterOverP, &b_gsf_eSuperClusterOverP);
   fChain->SetBranchAddress("gsf_dxy", &gsf_dxy, &b_gsf_dxy);
   fChain->SetBranchAddress("gsf_dxy_beamSpot", &gsf_dxy_beamSpot, &b_gsf_dxy_beamSpot);
   fChain->SetBranchAddress("gsf_dxy_firstPVtx", &gsf_dxy_firstPVtx, &b_gsf_dxy_firstPVtx);
   fChain->SetBranchAddress("gsf_dxyError", &gsf_dxyError, &b_gsf_dxyError);
   fChain->SetBranchAddress("gsf_dz", &gsf_dz, &b_gsf_dz);
   fChain->SetBranchAddress("gsf_dz_beamSpot", &gsf_dz_beamSpot, &b_gsf_dz_beamSpot);
   fChain->SetBranchAddress("gsf_dz_firstPVtx", &gsf_dz_firstPVtx, &b_gsf_dz_firstPVtx);
   fChain->SetBranchAddress("gsf_dzError", &gsf_dzError, &b_gsf_dzError);
   fChain->SetBranchAddress("gsf_vz", &gsf_vz, &b_gsf_vz);
   fChain->SetBranchAddress("gsf_numberOfValidHits", &gsf_numberOfValidHits, &b_gsf_numberOfValidHits);
   fChain->SetBranchAddress("gsf_nLostInnerHits", &gsf_nLostInnerHits, &b_gsf_nLostInnerHits);
   fChain->SetBranchAddress("gsf_nLostOuterHits", &gsf_nLostOuterHits, &b_gsf_nLostOuterHits);
   fChain->SetBranchAddress("gsf_convFlags", &gsf_convFlags, &b_gsf_convFlags);
   fChain->SetBranchAddress("gsf_convDist", &gsf_convDist, &b_gsf_convDist);
   fChain->SetBranchAddress("gsf_convDcot", &gsf_convDcot, &b_gsf_convDcot);
   fChain->SetBranchAddress("gsf_convRadius", &gsf_convRadius, &b_gsf_convRadius);
   fChain->SetBranchAddress("gsf_fBrem", &gsf_fBrem, &b_gsf_fBrem);
   fChain->SetBranchAddress("gsf_e1x5", &gsf_e1x5, &b_gsf_e1x5);
   fChain->SetBranchAddress("gsf_e2x5Max", &gsf_e2x5Max, &b_gsf_e2x5Max);
   fChain->SetBranchAddress("gsf_e5x5", &gsf_e5x5, &b_gsf_e5x5);
   fChain->SetBranchAddress("gsf_r9", &gsf_r9, &b_gsf_r9);
   fChain->SetBranchAddress("gsf_deltaEtaSeedClusterTrackAtVtx", &gsf_deltaEtaSeedClusterTrackAtVtx, &b_gsf_deltaEtaSeedClusterTrackAtVtx);
   fChain->SetBranchAddress("gsf_effArea", &gsf_effArea, &b_gsf_effArea);
   fChain->SetBranchAddress("gsf_sumChargedHadronPt", &gsf_sumChargedHadronPt, &b_gsf_sumChargedHadronPt);
   fChain->SetBranchAddress("gsf_sumNeutralHadronEt", &gsf_sumNeutralHadronEt, &b_gsf_sumNeutralHadronEt);
   fChain->SetBranchAddress("gsf_sumPhotonEt", &gsf_sumPhotonEt, &b_gsf_sumPhotonEt);
   fChain->SetBranchAddress("gsf_sumPUPt", &gsf_sumPUPt, &b_gsf_sumPUPt);
   fChain->SetBranchAddress("gsf_ooEmooP", &gsf_ooEmooP, &b_gsf_ooEmooP);
   fChain->SetBranchAddress("gsf_hitsinfo", &gsf_hitsinfo, &b_gsf_hitsinfo);
   fChain->SetBranchAddress("gsf_pixelMatch_dPhi1", &gsf_pixelMatch_dPhi1, &b_gsf_pixelMatch_dPhi1);
   fChain->SetBranchAddress("gsf_pixelMatch_dPhi2", &gsf_pixelMatch_dPhi2, &b_gsf_pixelMatch_dPhi2);
   fChain->SetBranchAddress("gsf_pixelMatch_dRz1", &gsf_pixelMatch_dRz1, &b_gsf_pixelMatch_dRz1);
   fChain->SetBranchAddress("gsf_pixelMatch_dRz2", &gsf_pixelMatch_dRz2, &b_gsf_pixelMatch_dRz2);
   fChain->SetBranchAddress("gsf_pixelMatch_subDetector1", &gsf_pixelMatch_subDetector1, &b_gsf_pixelMatch_subDetector1);
   fChain->SetBranchAddress("gsf_pixelMatch_subDetector2", &gsf_pixelMatch_subDetector2, &b_gsf_pixelMatch_subDetector2);
   fChain->SetBranchAddress("gsf_mc_bestDR", &gsf_mc_bestDR, &b_gsf_mc_bestDR);
   fChain->SetBranchAddress("gsf_mc_index", &gsf_mc_index, &b_gsf_mc_index);
   fChain->SetBranchAddress("gsf_mc_ERatio", &gsf_mc_ERatio, &b_gsf_mc_ERatio);
   fChain->SetBranchAddress("gsf_sc_energy", &gsf_sc_energy, &b_gsf_sc_energy);
   fChain->SetBranchAddress("gsf_sc_seed_eta", &gsf_sc_seed_eta, &b_gsf_sc_seed_eta);
   fChain->SetBranchAddress("gsf_sc_eta", &gsf_sc_eta, &b_gsf_sc_eta);
   fChain->SetBranchAddress("gsf_sc_etacorr", &gsf_sc_etacorr, &b_gsf_sc_etacorr);
   fChain->SetBranchAddress("gsf_sc_theta", &gsf_sc_theta, &b_gsf_sc_theta);
   fChain->SetBranchAddress("gsf_sc_thetacorr", &gsf_sc_thetacorr, &b_gsf_sc_thetacorr);
   fChain->SetBranchAddress("gsf_sc_et", &gsf_sc_et, &b_gsf_sc_et);
   fChain->SetBranchAddress("gsf_sc_phi", &gsf_sc_phi, &b_gsf_sc_phi);
   fChain->SetBranchAddress("gsf_sc_px", &gsf_sc_px, &b_gsf_sc_px);
   fChain->SetBranchAddress("gsf_sc_py", &gsf_sc_py, &b_gsf_sc_py);
   fChain->SetBranchAddress("gsf_sc_pz", &gsf_sc_pz, &b_gsf_sc_pz);
   fChain->SetBranchAddress("gsf_sc_x", &gsf_sc_x, &b_gsf_sc_x);
   fChain->SetBranchAddress("gsf_sc_y", &gsf_sc_y, &b_gsf_sc_y);
   fChain->SetBranchAddress("gsf_sc_z", &gsf_sc_z, &b_gsf_sc_z);
   fChain->SetBranchAddress("gsf_sc_phiWidth", &gsf_sc_phiWidth, &b_gsf_sc_phiWidth);
   fChain->SetBranchAddress("gsf_sc_etaWidth", &gsf_sc_etaWidth, &b_gsf_sc_etaWidth);
   fChain->SetBranchAddress("gsf_sc_seed_rawId", &gsf_sc_seed_rawId, &b_gsf_sc_seed_rawId);
   fChain->SetBranchAddress("gsf_sc_seed_ieta", &gsf_sc_seed_ieta, &b_gsf_sc_seed_ieta);
   fChain->SetBranchAddress("gsf_sc_seed_iphi", &gsf_sc_seed_iphi, &b_gsf_sc_seed_iphi);
   fChain->SetBranchAddress("gsf_sc_seed_kHasSwitchToGain6", &gsf_sc_seed_kHasSwitchToGain6, &b_gsf_sc_seed_kHasSwitchToGain6);
   fChain->SetBranchAddress("gsf_sc_seed_kHasSwitchToGain1", &gsf_sc_seed_kHasSwitchToGain1, &b_gsf_sc_seed_kHasSwitchToGain1);
   fChain->SetBranchAddress("EHits_isSaturated", &EHits_isSaturated, &b_EHits_isSaturated);
   fChain->SetBranchAddress("EBHits_rawId", &EBHits_rawId, &b_EBHits_rawId);
   fChain->SetBranchAddress("EBHits_iRechit", &EBHits_iRechit, &b_EBHits_iRechit);
   fChain->SetBranchAddress("EBHits_ieta", &EBHits_ieta, &b_EBHits_ieta);
   fChain->SetBranchAddress("EBHits_iphi", &EBHits_iphi, &b_EBHits_iphi);
   fChain->SetBranchAddress("EBHits_RecoFlag", &EBHits_RecoFlag, &b_EBHits_RecoFlag);
   fChain->SetBranchAddress("EBHits_kSaturated", &EBHits_kSaturated, &b_EBHits_kSaturated);
   fChain->SetBranchAddress("EBHits_kLeadingEdgeRecovered", &EBHits_kLeadingEdgeRecovered, &b_EBHits_kLeadingEdgeRecovered);
   fChain->SetBranchAddress("EBHits_kNeighboursRecovered", &EBHits_kNeighboursRecovered, &b_EBHits_kNeighboursRecovered);
   fChain->SetBranchAddress("EBHits_kWeird", &EBHits_kWeird, &b_EBHits_kWeird);
   fChain->SetBranchAddress("EBHits_energy", &EBHits_energy, &b_EBHits_energy);
   fChain->SetBranchAddress("EEHits_rawId", &EEHits_rawId, &b_EEHits_rawId);
   fChain->SetBranchAddress("EEHits_iRechit", &EEHits_iRechit, &b_EEHits_iRechit);
   fChain->SetBranchAddress("EEHits_ieta", &EEHits_ieta, &b_EEHits_ieta);
   fChain->SetBranchAddress("EEHits_iphi", &EEHits_iphi, &b_EEHits_iphi);
   fChain->SetBranchAddress("EEHits_RecoFlag", &EEHits_RecoFlag, &b_EEHits_RecoFlag);
   fChain->SetBranchAddress("EEHits_kSaturated", &EEHits_kSaturated, &b_EEHits_kSaturated);
   fChain->SetBranchAddress("EEHits_kLeadingEdgeRecovered", &EEHits_kLeadingEdgeRecovered, &b_EEHits_kLeadingEdgeRecovered);
   fChain->SetBranchAddress("EEHits_kNeighboursRecovered", &EEHits_kNeighboursRecovered, &b_EEHits_kNeighboursRecovered);
   fChain->SetBranchAddress("EEHits_kWeird", &EEHits_kWeird, &b_EEHits_kWeird);
   fChain->SetBranchAddress("EEHits_energy", &EEHits_energy, &b_EEHits_energy);
   fChain->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
   fChain->SetBranchAddress("mu_gt_qoverp", &mu_gt_qoverp, &b_mu_gt_qoverp);
   fChain->SetBranchAddress("mu_gt_charge", &mu_gt_charge, &b_mu_gt_charge);
   fChain->SetBranchAddress("mu_gt_pt", &mu_gt_pt, &b_mu_gt_pt);
   fChain->SetBranchAddress("mu_gt_eta", &mu_gt_eta, &b_mu_gt_eta);
   fChain->SetBranchAddress("mu_gt_phi", &mu_gt_phi, &b_mu_gt_phi);
   fChain->SetBranchAddress("mu_gt_p", &mu_gt_p, &b_mu_gt_p);
   fChain->SetBranchAddress("mu_gt_px", &mu_gt_px, &b_mu_gt_px);
   fChain->SetBranchAddress("mu_gt_py", &mu_gt_py, &b_mu_gt_py);
   fChain->SetBranchAddress("mu_gt_pz", &mu_gt_pz, &b_mu_gt_pz);
   fChain->SetBranchAddress("mu_gt_theta", &mu_gt_theta, &b_mu_gt_theta);
   fChain->SetBranchAddress("mu_gt_lambda", &mu_gt_lambda, &b_mu_gt_lambda);
   fChain->SetBranchAddress("mu_gt_d0", &mu_gt_d0, &b_mu_gt_d0);
   fChain->SetBranchAddress("mu_gt_dz", &mu_gt_dz, &b_mu_gt_dz);
   fChain->SetBranchAddress("mu_gt_dz_beamspot", &mu_gt_dz_beamspot, &b_mu_gt_dz_beamspot);
   fChain->SetBranchAddress("mu_gt_dz_firstPVtx", &mu_gt_dz_firstPVtx, &b_mu_gt_dz_firstPVtx);
   fChain->SetBranchAddress("mu_gt_dxy", &mu_gt_dxy, &b_mu_gt_dxy);
   fChain->SetBranchAddress("mu_gt_dxy_beamspot", &mu_gt_dxy_beamspot, &b_mu_gt_dxy_beamspot);
   fChain->SetBranchAddress("mu_gt_dxy_firstPVtx", &mu_gt_dxy_firstPVtx, &b_mu_gt_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_gt_dsz", &mu_gt_dsz, &b_mu_gt_dsz);
   fChain->SetBranchAddress("mu_gt_vx", &mu_gt_vx, &b_mu_gt_vx);
   fChain->SetBranchAddress("mu_gt_vy", &mu_gt_vy, &b_mu_gt_vy);
   fChain->SetBranchAddress("mu_gt_vz", &mu_gt_vz, &b_mu_gt_vz);
   fChain->SetBranchAddress("mu_gt_qoverpError", &mu_gt_qoverpError, &b_mu_gt_qoverpError);
   fChain->SetBranchAddress("mu_gt_ptError", &mu_gt_ptError, &b_mu_gt_ptError);
   fChain->SetBranchAddress("mu_gt_thetaError", &mu_gt_thetaError, &b_mu_gt_thetaError);
   fChain->SetBranchAddress("mu_gt_lambdaError", &mu_gt_lambdaError, &b_mu_gt_lambdaError);
   fChain->SetBranchAddress("mu_gt_phiError", &mu_gt_phiError, &b_mu_gt_phiError);
   fChain->SetBranchAddress("mu_gt_dxyError", &mu_gt_dxyError, &b_mu_gt_dxyError);
   fChain->SetBranchAddress("mu_gt_d0Error", &mu_gt_d0Error, &b_mu_gt_d0Error);
   fChain->SetBranchAddress("mu_gt_dszError", &mu_gt_dszError, &b_mu_gt_dszError);
   fChain->SetBranchAddress("mu_gt_dzError", &mu_gt_dzError, &b_mu_gt_dzError);
   fChain->SetBranchAddress("mu_gt_etaError", &mu_gt_etaError, &b_mu_gt_etaError);
   fChain->SetBranchAddress("mu_gt_chi2", &mu_gt_chi2, &b_mu_gt_chi2);
   fChain->SetBranchAddress("mu_gt_ndof", &mu_gt_ndof, &b_mu_gt_ndof);
   fChain->SetBranchAddress("mu_gt_normalizedChi2", &mu_gt_normalizedChi2, &b_mu_gt_normalizedChi2);
   fChain->SetBranchAddress("mu_ot_qoverp", &mu_ot_qoverp, &b_mu_ot_qoverp);
   fChain->SetBranchAddress("mu_ot_charge", &mu_ot_charge, &b_mu_ot_charge);
   fChain->SetBranchAddress("mu_ot_pt", &mu_ot_pt, &b_mu_ot_pt);
   fChain->SetBranchAddress("mu_ot_eta", &mu_ot_eta, &b_mu_ot_eta);
   fChain->SetBranchAddress("mu_ot_phi", &mu_ot_phi, &b_mu_ot_phi);
   fChain->SetBranchAddress("mu_ot_p", &mu_ot_p, &b_mu_ot_p);
   fChain->SetBranchAddress("mu_ot_px", &mu_ot_px, &b_mu_ot_px);
   fChain->SetBranchAddress("mu_ot_py", &mu_ot_py, &b_mu_ot_py);
   fChain->SetBranchAddress("mu_ot_pz", &mu_ot_pz, &b_mu_ot_pz);
   fChain->SetBranchAddress("mu_ot_theta", &mu_ot_theta, &b_mu_ot_theta);
   fChain->SetBranchAddress("mu_ot_lambda", &mu_ot_lambda, &b_mu_ot_lambda);
   fChain->SetBranchAddress("mu_ot_d0", &mu_ot_d0, &b_mu_ot_d0);
   fChain->SetBranchAddress("mu_ot_dz", &mu_ot_dz, &b_mu_ot_dz);
   fChain->SetBranchAddress("mu_ot_dz_beamspot", &mu_ot_dz_beamspot, &b_mu_ot_dz_beamspot);
   fChain->SetBranchAddress("mu_ot_dz_firstPVtx", &mu_ot_dz_firstPVtx, &b_mu_ot_dz_firstPVtx);
   fChain->SetBranchAddress("mu_ot_dxy", &mu_ot_dxy, &b_mu_ot_dxy);
   fChain->SetBranchAddress("mu_ot_dxy_beamspot", &mu_ot_dxy_beamspot, &b_mu_ot_dxy_beamspot);
   fChain->SetBranchAddress("mu_ot_dxy_firstPVtx", &mu_ot_dxy_firstPVtx, &b_mu_ot_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_ot_dsz", &mu_ot_dsz, &b_mu_ot_dsz);
   fChain->SetBranchAddress("mu_ot_vx", &mu_ot_vx, &b_mu_ot_vx);
   fChain->SetBranchAddress("mu_ot_vy", &mu_ot_vy, &b_mu_ot_vy);
   fChain->SetBranchAddress("mu_ot_vz", &mu_ot_vz, &b_mu_ot_vz);
   fChain->SetBranchAddress("mu_ot_qoverpError", &mu_ot_qoverpError, &b_mu_ot_qoverpError);
   fChain->SetBranchAddress("mu_ot_ptError", &mu_ot_ptError, &b_mu_ot_ptError);
   fChain->SetBranchAddress("mu_ot_thetaError", &mu_ot_thetaError, &b_mu_ot_thetaError);
   fChain->SetBranchAddress("mu_ot_lambdaError", &mu_ot_lambdaError, &b_mu_ot_lambdaError);
   fChain->SetBranchAddress("mu_ot_phiError", &mu_ot_phiError, &b_mu_ot_phiError);
   fChain->SetBranchAddress("mu_ot_dxyError", &mu_ot_dxyError, &b_mu_ot_dxyError);
   fChain->SetBranchAddress("mu_ot_d0Error", &mu_ot_d0Error, &b_mu_ot_d0Error);
   fChain->SetBranchAddress("mu_ot_dszError", &mu_ot_dszError, &b_mu_ot_dszError);
   fChain->SetBranchAddress("mu_ot_dzError", &mu_ot_dzError, &b_mu_ot_dzError);
   fChain->SetBranchAddress("mu_ot_etaError", &mu_ot_etaError, &b_mu_ot_etaError);
   fChain->SetBranchAddress("mu_ot_chi2", &mu_ot_chi2, &b_mu_ot_chi2);
   fChain->SetBranchAddress("mu_ot_ndof", &mu_ot_ndof, &b_mu_ot_ndof);
   fChain->SetBranchAddress("mu_ot_normalizedChi2", &mu_ot_normalizedChi2, &b_mu_ot_normalizedChi2);
   fChain->SetBranchAddress("mu_it_qoverp", &mu_it_qoverp, &b_mu_it_qoverp);
   fChain->SetBranchAddress("mu_it_charge", &mu_it_charge, &b_mu_it_charge);
   fChain->SetBranchAddress("mu_it_pt", &mu_it_pt, &b_mu_it_pt);
   fChain->SetBranchAddress("mu_it_eta", &mu_it_eta, &b_mu_it_eta);
   fChain->SetBranchAddress("mu_it_phi", &mu_it_phi, &b_mu_it_phi);
   fChain->SetBranchAddress("mu_it_p", &mu_it_p, &b_mu_it_p);
   fChain->SetBranchAddress("mu_it_px", &mu_it_px, &b_mu_it_px);
   fChain->SetBranchAddress("mu_it_py", &mu_it_py, &b_mu_it_py);
   fChain->SetBranchAddress("mu_it_pz", &mu_it_pz, &b_mu_it_pz);
   fChain->SetBranchAddress("mu_it_theta", &mu_it_theta, &b_mu_it_theta);
   fChain->SetBranchAddress("mu_it_lambda", &mu_it_lambda, &b_mu_it_lambda);
   fChain->SetBranchAddress("mu_it_d0", &mu_it_d0, &b_mu_it_d0);
   fChain->SetBranchAddress("mu_it_dz", &mu_it_dz, &b_mu_it_dz);
   fChain->SetBranchAddress("mu_it_dz_beamspot", &mu_it_dz_beamspot, &b_mu_it_dz_beamspot);
   fChain->SetBranchAddress("mu_it_dz_firstPVtx", &mu_it_dz_firstPVtx, &b_mu_it_dz_firstPVtx);
   fChain->SetBranchAddress("mu_it_dxy", &mu_it_dxy, &b_mu_it_dxy);
   fChain->SetBranchAddress("mu_it_dxy_beamspot", &mu_it_dxy_beamspot, &b_mu_it_dxy_beamspot);
   fChain->SetBranchAddress("mu_it_dxy_firstPVtx", &mu_it_dxy_firstPVtx, &b_mu_it_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_it_dsz", &mu_it_dsz, &b_mu_it_dsz);
   fChain->SetBranchAddress("mu_it_vx", &mu_it_vx, &b_mu_it_vx);
   fChain->SetBranchAddress("mu_it_vy", &mu_it_vy, &b_mu_it_vy);
   fChain->SetBranchAddress("mu_it_vz", &mu_it_vz, &b_mu_it_vz);
   fChain->SetBranchAddress("mu_it_qoverpError", &mu_it_qoverpError, &b_mu_it_qoverpError);
   fChain->SetBranchAddress("mu_it_ptError", &mu_it_ptError, &b_mu_it_ptError);
   fChain->SetBranchAddress("mu_it_thetaError", &mu_it_thetaError, &b_mu_it_thetaError);
   fChain->SetBranchAddress("mu_it_lambdaError", &mu_it_lambdaError, &b_mu_it_lambdaError);
   fChain->SetBranchAddress("mu_it_phiError", &mu_it_phiError, &b_mu_it_phiError);
   fChain->SetBranchAddress("mu_it_dxyError", &mu_it_dxyError, &b_mu_it_dxyError);
   fChain->SetBranchAddress("mu_it_d0Error", &mu_it_d0Error, &b_mu_it_d0Error);
   fChain->SetBranchAddress("mu_it_dszError", &mu_it_dszError, &b_mu_it_dszError);
   fChain->SetBranchAddress("mu_it_dzError", &mu_it_dzError, &b_mu_it_dzError);
   fChain->SetBranchAddress("mu_it_etaError", &mu_it_etaError, &b_mu_it_etaError);
   fChain->SetBranchAddress("mu_it_chi2", &mu_it_chi2, &b_mu_it_chi2);
   fChain->SetBranchAddress("mu_it_ndof", &mu_it_ndof, &b_mu_it_ndof);
   fChain->SetBranchAddress("mu_it_normalizedChi2", &mu_it_normalizedChi2, &b_mu_it_normalizedChi2);
   fChain->SetBranchAddress("mu_ibt_qoverp", &mu_ibt_qoverp, &b_mu_ibt_qoverp);
   fChain->SetBranchAddress("mu_ibt_charge", &mu_ibt_charge, &b_mu_ibt_charge);
   fChain->SetBranchAddress("mu_ibt_pt", &mu_ibt_pt, &b_mu_ibt_pt);
   fChain->SetBranchAddress("mu_ibt_eta", &mu_ibt_eta, &b_mu_ibt_eta);
   fChain->SetBranchAddress("mu_ibt_phi", &mu_ibt_phi, &b_mu_ibt_phi);
   fChain->SetBranchAddress("mu_ibt_p", &mu_ibt_p, &b_mu_ibt_p);
   fChain->SetBranchAddress("mu_ibt_px", &mu_ibt_px, &b_mu_ibt_px);
   fChain->SetBranchAddress("mu_ibt_py", &mu_ibt_py, &b_mu_ibt_py);
   fChain->SetBranchAddress("mu_ibt_pz", &mu_ibt_pz, &b_mu_ibt_pz);
   fChain->SetBranchAddress("mu_ibt_theta", &mu_ibt_theta, &b_mu_ibt_theta);
   fChain->SetBranchAddress("mu_ibt_lambda", &mu_ibt_lambda, &b_mu_ibt_lambda);
   fChain->SetBranchAddress("mu_ibt_d0", &mu_ibt_d0, &b_mu_ibt_d0);
   fChain->SetBranchAddress("mu_ibt_dz", &mu_ibt_dz, &b_mu_ibt_dz);
   fChain->SetBranchAddress("mu_ibt_dz_beamspot", &mu_ibt_dz_beamspot, &b_mu_ibt_dz_beamspot);
   fChain->SetBranchAddress("mu_ibt_dz_firstPVtx", &mu_ibt_dz_firstPVtx, &b_mu_ibt_dz_firstPVtx);
   fChain->SetBranchAddress("mu_ibt_dxy", &mu_ibt_dxy, &b_mu_ibt_dxy);
   fChain->SetBranchAddress("mu_ibt_dxy_beamspot", &mu_ibt_dxy_beamspot, &b_mu_ibt_dxy_beamspot);
   fChain->SetBranchAddress("mu_ibt_dxy_firstPVtx", &mu_ibt_dxy_firstPVtx, &b_mu_ibt_dxy_firstPVtx);
   fChain->SetBranchAddress("mu_ibt_dsz", &mu_ibt_dsz, &b_mu_ibt_dsz);
   fChain->SetBranchAddress("mu_ibt_vx", &mu_ibt_vx, &b_mu_ibt_vx);
   fChain->SetBranchAddress("mu_ibt_vy", &mu_ibt_vy, &b_mu_ibt_vy);
   fChain->SetBranchAddress("mu_ibt_vz", &mu_ibt_vz, &b_mu_ibt_vz);
   fChain->SetBranchAddress("mu_ibt_qoverpError", &mu_ibt_qoverpError, &b_mu_ibt_qoverpError);
   fChain->SetBranchAddress("mu_ibt_ptError", &mu_ibt_ptError, &b_mu_ibt_ptError);
   fChain->SetBranchAddress("mu_ibt_thetaError", &mu_ibt_thetaError, &b_mu_ibt_thetaError);
   fChain->SetBranchAddress("mu_ibt_lambdaError", &mu_ibt_lambdaError, &b_mu_ibt_lambdaError);
   fChain->SetBranchAddress("mu_ibt_phiError", &mu_ibt_phiError, &b_mu_ibt_phiError);
   fChain->SetBranchAddress("mu_ibt_dxyError", &mu_ibt_dxyError, &b_mu_ibt_dxyError);
   fChain->SetBranchAddress("mu_ibt_d0Error", &mu_ibt_d0Error, &b_mu_ibt_d0Error);
   fChain->SetBranchAddress("mu_ibt_dszError", &mu_ibt_dszError, &b_mu_ibt_dszError);
   fChain->SetBranchAddress("mu_ibt_dzError", &mu_ibt_dzError, &b_mu_ibt_dzError);
   fChain->SetBranchAddress("mu_ibt_etaError", &mu_ibt_etaError, &b_mu_ibt_etaError);
   fChain->SetBranchAddress("mu_ibt_chi2", &mu_ibt_chi2, &b_mu_ibt_chi2);
   fChain->SetBranchAddress("mu_ibt_ndof", &mu_ibt_ndof, &b_mu_ibt_ndof);
   fChain->SetBranchAddress("mu_ibt_normalizedChi2", &mu_ibt_normalizedChi2, &b_mu_ibt_normalizedChi2);
   fChain->SetBranchAddress("mu_isGlobalMuon", &mu_isGlobalMuon, &b_mu_isGlobalMuon);
   fChain->SetBranchAddress("mu_isStandAloneMuon", &mu_isStandAloneMuon, &b_mu_isStandAloneMuon);
   fChain->SetBranchAddress("mu_isTrackerMuon", &mu_isTrackerMuon, &b_mu_isTrackerMuon);
   fChain->SetBranchAddress("mu_isPFMuon", &mu_isPFMuon, &b_mu_isPFMuon);
   fChain->SetBranchAddress("mu_isPFIsolationValid", &mu_isPFIsolationValid, &b_mu_isPFIsolationValid);
   fChain->SetBranchAddress("mu_isGoodMuonTMLastStationLoose", &mu_isGoodMuonTMLastStationLoose, &b_mu_isGoodMuonTMLastStationLoose);
   fChain->SetBranchAddress("mu_isGoodMuonTMLastStationTight", &mu_isGoodMuonTMLastStationTight, &b_mu_isGoodMuonTMLastStationTight);
   fChain->SetBranchAddress("mu_isGoodMuonTM2DCompatibilityLoose", &mu_isGoodMuonTM2DCompatibilityLoose, &b_mu_isGoodMuonTM2DCompatibilityLoose);
   fChain->SetBranchAddress("mu_isGoodMuonTM2DCompatibilityTight", &mu_isGoodMuonTM2DCompatibilityTight, &b_mu_isGoodMuonTM2DCompatibilityTight);
   fChain->SetBranchAddress("mu_isGoodMuonTMOneStationLoose", &mu_isGoodMuonTMOneStationLoose, &b_mu_isGoodMuonTMOneStationLoose);
   fChain->SetBranchAddress("mu_isGoodMuonTMOneStationTight", &mu_isGoodMuonTMOneStationTight, &b_mu_isGoodMuonTMOneStationTight);
   fChain->SetBranchAddress("mu_isGoodMuonTMLastStationOptimizedLowPtLoose", &mu_isGoodMuonTMLastStationOptimizedLowPtLoose, &b_mu_isGoodMuonTMLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("mu_isGoodMuonTMLastStationOptimizedLowPtTight", &mu_isGoodMuonTMLastStationOptimizedLowPtTight, &b_mu_isGoodMuonTMLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("mu_isTightMuon", &mu_isTightMuon, &b_mu_isTightMuon);
   fChain->SetBranchAddress("mu_isMediumMuon", &mu_isMediumMuon, &b_mu_isMediumMuon);
   fChain->SetBranchAddress("mu_isLooseMuon", &mu_isLooseMuon, &b_mu_isLooseMuon);
   fChain->SetBranchAddress("mu_isSoftMuon", &mu_isSoftMuon, &b_mu_isSoftMuon);
   fChain->SetBranchAddress("mu_isHighPtMuonOld", &mu_isHighPtMuonOld, &b_mu_isHighPtMuonOld);
   fChain->SetBranchAddress("mu_isHighPtMuon", &mu_isHighPtMuon, &b_mu_isHighPtMuon);
   fChain->SetBranchAddress("mu_isTrackerHighPtMuon", &mu_isTrackerHighPtMuon, &b_mu_isTrackerHighPtMuon);
   fChain->SetBranchAddress("mu_CutBasedIdLoose", &mu_CutBasedIdLoose, &b_mu_CutBasedIdLoose);
   fChain->SetBranchAddress("mu_CutBasedIdMedium", &mu_CutBasedIdMedium, &b_mu_CutBasedIdMedium);
   fChain->SetBranchAddress("mu_CutBasedIdMediumPrompt", &mu_CutBasedIdMediumPrompt, &b_mu_CutBasedIdMediumPrompt);
   fChain->SetBranchAddress("mu_CutBasedIdTight", &mu_CutBasedIdTight, &b_mu_CutBasedIdTight);
   fChain->SetBranchAddress("mu_CutBasedIdGlobalHighPt", &mu_CutBasedIdGlobalHighPt, &b_mu_CutBasedIdGlobalHighPt);
   fChain->SetBranchAddress("mu_CutBasedIdTrkHighPt", &mu_CutBasedIdTrkHighPt, &b_mu_CutBasedIdTrkHighPt);
   fChain->SetBranchAddress("mu_PFIsoVeryLoose", &mu_PFIsoVeryLoose, &b_mu_PFIsoVeryLoose);
   fChain->SetBranchAddress("mu_PFIsoLoose", &mu_PFIsoLoose, &b_mu_PFIsoLoose);
   fChain->SetBranchAddress("mu_PFIsoMedium", &mu_PFIsoMedium, &b_mu_PFIsoMedium);
   fChain->SetBranchAddress("mu_PFIsoTight", &mu_PFIsoTight, &b_mu_PFIsoTight);
   fChain->SetBranchAddress("mu_PFIsoVeryTight", &mu_PFIsoVeryTight, &b_mu_PFIsoVeryTight);
   fChain->SetBranchAddress("mu_PFIsoVeryVeryTight", &mu_PFIsoVeryVeryTight, &b_mu_PFIsoVeryVeryTight);
   fChain->SetBranchAddress("mu_TkIsoLoose", &mu_TkIsoLoose, &b_mu_TkIsoLoose);
   fChain->SetBranchAddress("mu_TkIsoTight", &mu_TkIsoTight, &b_mu_TkIsoTight);
   fChain->SetBranchAddress("mu_SoftCutBasedId", &mu_SoftCutBasedId, &b_mu_SoftCutBasedId);
   fChain->SetBranchAddress("mu_SoftMvaId", &mu_SoftMvaId, &b_mu_SoftMvaId);
   fChain->SetBranchAddress("mu_MvaLoose", &mu_MvaLoose, &b_mu_MvaLoose);
   fChain->SetBranchAddress("mu_MvaMedium", &mu_MvaMedium, &b_mu_MvaMedium);
   fChain->SetBranchAddress("mu_MvaTight", &mu_MvaTight, &b_mu_MvaTight);
   fChain->SetBranchAddress("mu_MiniIsoLoose", &mu_MiniIsoLoose, &b_mu_MiniIsoLoose);
   fChain->SetBranchAddress("mu_MiniIsoMedium", &mu_MiniIsoMedium, &b_mu_MiniIsoMedium);
   fChain->SetBranchAddress("mu_MiniIsoTight", &mu_MiniIsoTight, &b_mu_MiniIsoTight);
   fChain->SetBranchAddress("mu_MiniIsoVeryTight", &mu_MiniIsoVeryTight, &b_mu_MiniIsoVeryTight);
   fChain->SetBranchAddress("mu_TriggerIdLoose", &mu_TriggerIdLoose, &b_mu_TriggerIdLoose);
   fChain->SetBranchAddress("mu_InTimeMuon", &mu_InTimeMuon, &b_mu_InTimeMuon);
   fChain->SetBranchAddress("mu_MultiIsoLoose", &mu_MultiIsoLoose, &b_mu_MultiIsoLoose);
   fChain->SetBranchAddress("mu_MultiIsoMedium", &mu_MultiIsoMedium, &b_mu_MultiIsoMedium);
   fChain->SetBranchAddress("mu_numberOfMatchedStations", &mu_numberOfMatchedStations, &b_mu_numberOfMatchedStations);
   fChain->SetBranchAddress("mu_numberOfValidPixelHits", &mu_numberOfValidPixelHits, &b_mu_numberOfValidPixelHits);
   fChain->SetBranchAddress("mu_trackerLayersWithMeasurement", &mu_trackerLayersWithMeasurement, &b_mu_trackerLayersWithMeasurement);
   fChain->SetBranchAddress("mu_numberOfValidMuonHits", &mu_numberOfValidMuonHits, &b_mu_numberOfValidMuonHits);
   fChain->SetBranchAddress("mu_pixelLayersWithMeasurement", &mu_pixelLayersWithMeasurement, &b_mu_pixelLayersWithMeasurement);
   fChain->SetBranchAddress("mu_innerTrack_validFraction", &mu_innerTrack_validFraction, &b_mu_innerTrack_validFraction);
   fChain->SetBranchAddress("mu_combinedQuality_trkKink", &mu_combinedQuality_trkKink, &b_mu_combinedQuality_trkKink);
   fChain->SetBranchAddress("mu_combinedQuality_chi2LocalPosition", &mu_combinedQuality_chi2LocalPosition, &b_mu_combinedQuality_chi2LocalPosition);
   fChain->SetBranchAddress("mu_segmentCompatibility", &mu_segmentCompatibility, &b_mu_segmentCompatibility);
   fChain->SetBranchAddress("mu_dB", &mu_dB, &b_mu_dB);
   fChain->SetBranchAddress("mu_pt_default", &mu_pt_default, &b_mu_pt_default);
   fChain->SetBranchAddress("mu_isolationR03_sumPt", &mu_isolationR03_sumPt, &b_mu_isolationR03_sumPt);
   fChain->SetBranchAddress("mu_isolationR03_trackerVetoPt", &mu_isolationR03_trackerVetoPt, &b_mu_isolationR03_trackerVetoPt);
   fChain->SetBranchAddress("mu_isolationR03_emEt", &mu_isolationR03_emEt, &b_mu_isolationR03_emEt);
   fChain->SetBranchAddress("mu_isolationR03_emVetoEt", &mu_isolationR03_emVetoEt, &b_mu_isolationR03_emVetoEt);
   fChain->SetBranchAddress("mu_isolationR03_hadEt", &mu_isolationR03_hadEt, &b_mu_isolationR03_hadEt);
   fChain->SetBranchAddress("mu_isolationR03_hadVetoEt", &mu_isolationR03_hadVetoEt, &b_mu_isolationR03_hadVetoEt);
   fChain->SetBranchAddress("mu_isolationR05_sumPt", &mu_isolationR05_sumPt, &b_mu_isolationR05_sumPt);
   fChain->SetBranchAddress("mu_isolationR05_trackerVetoPt", &mu_isolationR05_trackerVetoPt, &b_mu_isolationR05_trackerVetoPt);
   fChain->SetBranchAddress("mu_isolationR05_emEt", &mu_isolationR05_emEt, &b_mu_isolationR05_emEt);
   fChain->SetBranchAddress("mu_isolationR05_emVetoEt", &mu_isolationR05_emVetoEt, &b_mu_isolationR05_emVetoEt);
   fChain->SetBranchAddress("mu_isolationR05_hadEt", &mu_isolationR05_hadEt, &b_mu_isolationR05_hadEt);
   fChain->SetBranchAddress("mu_isolationR05_hadVetoEt", &mu_isolationR05_hadVetoEt, &b_mu_isolationR05_hadVetoEt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumChargedHadronPt", &mu_pfIsolationR03_sumChargedHadronPt, &b_mu_pfIsolationR03_sumChargedHadronPt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumNeutralHadronEt", &mu_pfIsolationR03_sumNeutralHadronEt, &b_mu_pfIsolationR03_sumNeutralHadronEt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumChargedParticlePt", &mu_pfIsolationR03_sumChargedParticlePt, &b_mu_pfIsolationR03_sumChargedParticlePt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumPhotonEt", &mu_pfIsolationR03_sumPhotonEt, &b_mu_pfIsolationR03_sumPhotonEt);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumNeutralHadronEtHighThreshold", &mu_pfIsolationR03_sumNeutralHadronEtHighThreshold, &b_mu_pfIsolationR03_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumPhotonEtHighThreshold", &mu_pfIsolationR03_sumPhotonEtHighThreshold, &b_mu_pfIsolationR03_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR03_sumPUPt", &mu_pfIsolationR03_sumPUPt, &b_mu_pfIsolationR03_sumPUPt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumChargedHadronPt", &mu_pfIsolationR04_sumChargedHadronPt, &b_mu_pfIsolationR04_sumChargedHadronPt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumNeutralHadronEt", &mu_pfIsolationR04_sumNeutralHadronEt, &b_mu_pfIsolationR04_sumNeutralHadronEt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumChargedParticlePt", &mu_pfIsolationR04_sumChargedParticlePt, &b_mu_pfIsolationR04_sumChargedParticlePt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumPhotonEt", &mu_pfIsolationR04_sumPhotonEt, &b_mu_pfIsolationR04_sumPhotonEt);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumNeutralHadronEtHighThreshold", &mu_pfIsolationR04_sumNeutralHadronEtHighThreshold, &b_mu_pfIsolationR04_sumNeutralHadronEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumPhotonEtHighThreshold", &mu_pfIsolationR04_sumPhotonEtHighThreshold, &b_mu_pfIsolationR04_sumPhotonEtHighThreshold);
   fChain->SetBranchAddress("mu_pfIsolationR04_sumPUPt", &mu_pfIsolationR04_sumPUPt, &b_mu_pfIsolationR04_sumPUPt);
   fChain->SetBranchAddress("mu_pfIsoDbCorrected03", &mu_pfIsoDbCorrected03, &b_mu_pfIsoDbCorrected03);
   fChain->SetBranchAddress("mu_pfIsoDbCorrected04", &mu_pfIsoDbCorrected04, &b_mu_pfIsoDbCorrected04);
   fChain->SetBranchAddress("mu_isoTrackerBased03", &mu_isoTrackerBased03, &b_mu_isoTrackerBased03);
   fChain->SetBranchAddress("mu_mc_bestDR", &mu_mc_bestDR, &b_mu_mc_bestDR);
   fChain->SetBranchAddress("mu_mc_index", &mu_mc_index, &b_mu_mc_index);
   fChain->SetBranchAddress("mu_mc_ERatio", &mu_mc_ERatio, &b_mu_mc_ERatio);
   fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
   fChain->SetBranchAddress("jet_px", &jet_px, &b_jet_px);
   fChain->SetBranchAddress("jet_py", &jet_py, &b_jet_py);
   fChain->SetBranchAddress("jet_pz", &jet_pz, &b_jet_pz);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_theta", &jet_theta, &b_jet_theta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_energy", &jet_energy, &b_jet_energy);
   fChain->SetBranchAddress("jet_mass", &jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("jet_chargedEmEnergyFraction", &jet_chargedEmEnergyFraction, &b_jet_chargedEmEnergyFraction);
   fChain->SetBranchAddress("jet_neutralHadronEnergyFraction", &jet_neutralHadronEnergyFraction, &b_jet_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("jet_neutralEmEnergyFraction", &jet_neutralEmEnergyFraction, &b_jet_neutralEmEnergyFraction);
   fChain->SetBranchAddress("jet_chargedHadronEnergyFraction", &jet_chargedHadronEnergyFraction, &b_jet_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("jet_muonEnergyFraction", &jet_muonEnergyFraction, &b_jet_muonEnergyFraction);
   fChain->SetBranchAddress("jet_chargedMultiplicity", &jet_chargedMultiplicity, &b_jet_chargedMultiplicity);
   fChain->SetBranchAddress("jet_neutralMultiplicity", &jet_neutralMultiplicity, &b_jet_neutralMultiplicity);
   fChain->SetBranchAddress("jet_partonFlavour", &jet_partonFlavour, &b_jet_partonFlavour);
   fChain->SetBranchAddress("jet_hadronFlavour", &jet_hadronFlavour, &b_jet_hadronFlavour);
   fChain->SetBranchAddress("jet_CSVv2", &jet_CSVv2, &b_jet_CSVv2);
   fChain->SetBranchAddress("jet_CvsL", &jet_CvsL, &b_jet_CvsL);
   fChain->SetBranchAddress("jet_CvsB", &jet_CvsB, &b_jet_CvsB);
   fChain->SetBranchAddress("jet_MVA2BJets", &jet_MVA2BJets, &b_jet_MVA2BJets);
   fChain->SetBranchAddress("jet_CvsB_DeepJet_charm_tagger", &jet_CvsB_DeepJet_charm_tagger, &b_jet_CvsB_DeepJet_charm_tagger);
   fChain->SetBranchAddress("jet_CvsL_DeepJet_charm_tagger", &jet_CvsL_DeepJet_charm_tagger, &b_jet_CvsL_DeepJet_charm_tagger);
   fChain->SetBranchAddress("jet_CvsB_DeepCSV_charm_tagger", &jet_CvsB_DeepCSV_charm_tagger, &b_jet_CvsB_DeepCSV_charm_tagger);
   fChain->SetBranchAddress("jet_CvsL_DeepCSV_charm_tagger", &jet_CvsL_DeepCSV_charm_tagger, &b_jet_CvsL_DeepCSV_charm_tagger);
   fChain->SetBranchAddress("jet_DeepJet", &jet_DeepJet, &b_jet_DeepJet);
   fChain->SetBranchAddress("jet_DeepCSV", &jet_DeepCSV, &b_jet_DeepCSV);
   fChain->SetBranchAddress("jet_isJetIDLoose_2016", &jet_isJetIDLoose_2016, &b_jet_isJetIDLoose_2016);
   fChain->SetBranchAddress("jet_isJetIDTight_2016", &jet_isJetIDTight_2016, &b_jet_isJetIDTight_2016);
   fChain->SetBranchAddress("jet_isJetIDTightLepVeto_2016", &jet_isJetIDTightLepVeto_2016, &b_jet_isJetIDTightLepVeto_2016);
   fChain->SetBranchAddress("jet_isJetID_2017", &jet_isJetID_2017, &b_jet_isJetID_2017);
   fChain->SetBranchAddress("jet_isJetIDLepVeto_2017", &jet_isJetIDLepVeto_2017, &b_jet_isJetIDLepVeto_2017);
   fChain->SetBranchAddress("jet_isJetID_2018", &jet_isJetID_2018, &b_jet_isJetID_2018);
   fChain->SetBranchAddress("jet_isJetIDLepVeto_2018", &jet_isJetIDLepVeto_2018, &b_jet_isJetIDLepVeto_2018);
   fChain->SetBranchAddress("jet_Smeared_pt", &jet_Smeared_pt, &b_jet_Smeared_pt);
   fChain->SetBranchAddress("jet_SmearedJetResUp_pt", &jet_SmearedJetResUp_pt, &b_jet_SmearedJetResUp_pt);
   fChain->SetBranchAddress("jet_SmearedJetResDown_pt", &jet_SmearedJetResDown_pt, &b_jet_SmearedJetResDown_pt);
   fChain->SetBranchAddress("jet_SmearedJetEnUp_pt", &jet_SmearedJetEnUp_pt, &b_jet_SmearedJetEnUp_pt);
   fChain->SetBranchAddress("jet_SmearedJetEnDown_pt", &jet_SmearedJetEnDown_pt, &b_jet_SmearedJetEnDown_pt);
//   fChain->SetBranchAddress("BtagSF_Deepcsv_loose", &BtagSF_Deepcsv_loose, &b_BtagSF_Deepcsv_loose);
//   fChain->SetBranchAddress("BtagSF_Deepcsv_Up_loose", &BtagSF_Deepcsv_Up_loose, &b_BtagSF_Deepcsv_Up_loose);
//   fChain->SetBranchAddress("BtagSF_Deepcsv_Down_loose", &BtagSF_Deepcsv_Down_loose, &b_BtagSF_Deepcsv_Down_loose);
//   fChain->SetBranchAddress("BtagSF_Deepcsv_medium", &BtagSF_Deepcsv_medium, &b_BtagSF_Deepcsv_medium);
//   fChain->SetBranchAddress("BtagSF_Deepcsv_Up_medium", &BtagSF_Deepcsv_Up_medium, &b_BtagSF_Deepcsv_Up_medium);
//   fChain->SetBranchAddress("BtagSF_Deepcsv_Down_medium", &BtagSF_Deepcsv_Down_medium, &b_BtagSF_Deepcsv_Down_medium);
//   fChain->SetBranchAddress("BtagSF_Deepcsv_tight", &BtagSF_Deepcsv_tight, &b_BtagSF_Deepcsv_tight);
//   fChain->SetBranchAddress("BtagSF_Deepcsv_Up_tight", &BtagSF_Deepcsv_Up_tight, &b_BtagSF_Deepcsv_Up_tight);
//   fChain->SetBranchAddress("BtagSF_Deepcsv_Down_tight", &BtagSF_Deepcsv_Down_tight, &b_BtagSF_Deepcsv_Down_tight);
   fChain->SetBranchAddress("MET_nominal_Pt", &MET_nominal_Pt, &b_MET_nominal_Pt);
   fChain->SetBranchAddress("MET_nominal_phi", &MET_nominal_phi, &b_MET_nominal_phi);
   fChain->SetBranchAddress("MET_nominal_significance", &MET_nominal_significance, &b_MET_nominal_significance);
   fChain->SetBranchAddress("MET_Pt", &MET_Pt, &b_MET_Pt);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   fChain->SetBranchAddress("MET_T1_Pt", &MET_T1_Pt, &b_MET_T1_Pt);
   fChain->SetBranchAddress("MET_T1_phi", &MET_T1_phi, &b_MET_T1_phi);
   fChain->SetBranchAddress("MET_T1_significance", &MET_T1_significance, &b_MET_T1_significance);
   fChain->SetBranchAddress("MET_gen_pt", &MET_gen_pt, &b_MET_gen_pt);
   fChain->SetBranchAddress("MET_gen_phi", &MET_gen_phi, &b_MET_gen_phi);
   fChain->SetBranchAddress("MET_Type1Unc_Px", &MET_Type1Unc_Px, &b_MET_Type1Unc_Px);
   fChain->SetBranchAddress("MET_Type1Unc_Py", &MET_Type1Unc_Py, &b_MET_Type1Unc_Py);
   fChain->SetBranchAddress("MET_Type1Unc_Pt", &MET_Type1Unc_Pt, &b_MET_Type1Unc_Pt);
   fChain->SetBranchAddress("MET_T1JetEnDown_Pt", &MET_T1JetEnDown_Pt, &b_MET_T1JetEnDown_Pt);
   fChain->SetBranchAddress("MET_T1JetEnDown_phi", &MET_T1JetEnDown_phi, &b_MET_T1JetEnDown_phi);
   fChain->SetBranchAddress("MET_T1JetEnDown_significance", &MET_T1JetEnDown_significance, &b_MET_T1JetEnDown_significance);
   fChain->SetBranchAddress("MET_T1JetEnUp_Pt", &MET_T1JetEnUp_Pt, &b_MET_T1JetEnUp_Pt);
   fChain->SetBranchAddress("MET_T1JetEnUp_phi", &MET_T1JetEnUp_phi, &b_MET_T1JetEnUp_phi);
   fChain->SetBranchAddress("MET_T1JetEnUp_significance", &MET_T1JetEnUp_significance, &b_MET_T1JetEnUp_significance);
   fChain->SetBranchAddress("MET_T1Smear_Pt", &MET_T1Smear_Pt, &b_MET_T1Smear_Pt);
   fChain->SetBranchAddress("MET_T1Smear_phi", &MET_T1Smear_phi, &b_MET_T1Smear_phi);
   fChain->SetBranchAddress("MET_T1Smear_significance", &MET_T1Smear_significance, &b_MET_T1Smear_significance);
   fChain->SetBranchAddress("MET_T1SmearJetEnDown_Pt", &MET_T1SmearJetEnDown_Pt, &b_MET_T1SmearJetEnDown_Pt);
   fChain->SetBranchAddress("MET_T1SmearJetEnDown_phi", &MET_T1SmearJetEnDown_phi, &b_MET_T1SmearJetEnDown_phi);
   fChain->SetBranchAddress("MET_T1SmearJetEnDown_significance", &MET_T1SmearJetEnDown_significance, &b_MET_T1SmearJetEnDown_significance);
   fChain->SetBranchAddress("MET_T1SmearJetEnUp_Pt", &MET_T1SmearJetEnUp_Pt, &b_MET_T1SmearJetEnUp_Pt);
   fChain->SetBranchAddress("MET_T1SmearJetEnUp_phi", &MET_T1SmearJetEnUp_phi, &b_MET_T1SmearJetEnUp_phi);
   fChain->SetBranchAddress("MET_T1SmearJetEnUp_significance", &MET_T1SmearJetEnUp_significance, &b_MET_T1SmearJetEnUp_significance);
   fChain->SetBranchAddress("MET_T1SmearJetResDown_Pt", &MET_T1SmearJetResDown_Pt, &b_MET_T1SmearJetResDown_Pt);
   fChain->SetBranchAddress("MET_T1SmearJetResDown_phi", &MET_T1SmearJetResDown_phi, &b_MET_T1SmearJetResDown_phi);
   fChain->SetBranchAddress("MET_T1SmearJetResDown_significance", &MET_T1SmearJetResDown_significance, &b_MET_T1SmearJetResDown_significance);
   fChain->SetBranchAddress("MET_T1SmearJetResUp_Pt", &MET_T1SmearJetResUp_Pt, &b_MET_T1SmearJetResUp_Pt);
   fChain->SetBranchAddress("MET_T1SmearJetResUp_phi", &MET_T1SmearJetResUp_phi, &b_MET_T1SmearJetResUp_phi);
   fChain->SetBranchAddress("MET_T1SmearJetResUp_significance", &MET_T1SmearJetResUp_significance, &b_MET_T1SmearJetResUp_significance);
   fChain->SetBranchAddress("MET_T1Txy_Pt", &MET_T1Txy_Pt, &b_MET_T1Txy_Pt);
   fChain->SetBranchAddress("MET_T1Txy_phi", &MET_T1Txy_phi, &b_MET_T1Txy_phi);
   fChain->SetBranchAddress("MET_T1Txy_significance", &MET_T1Txy_significance, &b_MET_T1Txy_significance);
   fChain->SetBranchAddress("MET_FinalCollection_Pt", &MET_FinalCollection_Pt, &b_MET_FinalCollection_Pt);
   fChain->SetBranchAddress("MET_FinalCollection_phi", &MET_FinalCollection_phi, &b_MET_FinalCollection_phi);
   fChain->SetBranchAddress("MET_FinalCollection_significance", &MET_FinalCollection_significance, &b_MET_FinalCollection_significance);
   fChain->SetBranchAddress("trig_Flag_ecalBadCalibReduced", &trig_Flag_ecalBadCalibReduced, &b_trig_Flag_ecalBadCalibReduced);
   fChain->SetBranchAddress("trig_Flag_HBHENoiseFilter_accept", &trig_Flag_HBHENoiseFilter_accept, &b_trig_Flag_HBHENoiseFilter_accept);
   fChain->SetBranchAddress("trig_Flag_HBHENoiseIsoFilter_accept", &trig_Flag_HBHENoiseIsoFilter_accept, &b_trig_Flag_HBHENoiseIsoFilter_accept);
   fChain->SetBranchAddress("trig_Flag_CSCTightHaloFilter_accept", &trig_Flag_CSCTightHaloFilter_accept, &b_trig_Flag_CSCTightHaloFilter_accept);
   fChain->SetBranchAddress("trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept", &trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept, &b_trig_Flag_CSCTightHaloTrkMuUnvetoFilter_accept);
   fChain->SetBranchAddress("trig_Flag_CSCTightHalo2015Filter_accept", &trig_Flag_CSCTightHalo2015Filter_accept, &b_trig_Flag_CSCTightHalo2015Filter_accept);
   fChain->SetBranchAddress("trig_Flag_globalTightHalo2016Filter_accept", &trig_Flag_globalTightHalo2016Filter_accept, &b_trig_Flag_globalTightHalo2016Filter_accept);
   fChain->SetBranchAddress("trig_Flag_globalSuperTightHalo2016Filter_accept", &trig_Flag_globalSuperTightHalo2016Filter_accept, &b_trig_Flag_globalSuperTightHalo2016Filter_accept);
   fChain->SetBranchAddress("trig_Flag_HcalStripHaloFilter_accept", &trig_Flag_HcalStripHaloFilter_accept, &b_trig_Flag_HcalStripHaloFilter_accept);
   fChain->SetBranchAddress("trig_Flag_hcalLaserEventFilter_accept", &trig_Flag_hcalLaserEventFilter_accept, &b_trig_Flag_hcalLaserEventFilter_accept);
   fChain->SetBranchAddress("trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept", &trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept, &b_trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept);
   fChain->SetBranchAddress("trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept", &trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept, &b_trig_Flag_EcalDeadCellBoundaryEnergyFilter_accept);
   fChain->SetBranchAddress("trig_Flag_ecalBadCalibFilter_accept", &trig_Flag_ecalBadCalibFilter_accept, &b_trig_Flag_ecalBadCalibFilter_accept);
   fChain->SetBranchAddress("trig_Flag_goodVertices_accept", &trig_Flag_goodVertices_accept, &b_trig_Flag_goodVertices_accept);
   fChain->SetBranchAddress("trig_Flag_eeBadScFilter_accept", &trig_Flag_eeBadScFilter_accept, &b_trig_Flag_eeBadScFilter_accept);
   fChain->SetBranchAddress("trig_Flag_ecalLaserCorrFilter_accept", &trig_Flag_ecalLaserCorrFilter_accept, &b_trig_Flag_ecalLaserCorrFilter_accept);
   fChain->SetBranchAddress("trig_Flag_trkPOGFilters_accept", &trig_Flag_trkPOGFilters_accept, &b_trig_Flag_trkPOGFilters_accept);
   fChain->SetBranchAddress("trig_Flag_chargedHadronTrackResolutionFilter_accept", &trig_Flag_chargedHadronTrackResolutionFilter_accept, &b_trig_Flag_chargedHadronTrackResolutionFilter_accept);
   fChain->SetBranchAddress("trig_Flag_muonBadTrackFilter_accept", &trig_Flag_muonBadTrackFilter_accept, &b_trig_Flag_muonBadTrackFilter_accept);
   fChain->SetBranchAddress("trig_Flag_BadChargedCandidateFilter_accept", &trig_Flag_BadChargedCandidateFilter_accept, &b_trig_Flag_BadChargedCandidateFilter_accept);
   fChain->SetBranchAddress("trig_Flag_BadPFMuonFilter_accept", &trig_Flag_BadPFMuonFilter_accept, &b_trig_Flag_BadPFMuonFilter_accept);
   fChain->SetBranchAddress("trig_Flag_BadChargedCandidateSummer16Filter_accept", &trig_Flag_BadChargedCandidateSummer16Filter_accept, &b_trig_Flag_BadChargedCandidateSummer16Filter_accept);
   fChain->SetBranchAddress("trig_Flag_BadPFMuonSummer16Filter_accept", &trig_Flag_BadPFMuonSummer16Filter_accept, &b_trig_Flag_BadPFMuonSummer16Filter_accept);
   fChain->SetBranchAddress("trig_Flag_trkPOG_manystripclus53X_accept", &trig_Flag_trkPOG_manystripclus53X_accept, &b_trig_Flag_trkPOG_manystripclus53X_accept);
   fChain->SetBranchAddress("trig_Flag_trkPOG_toomanystripclus53X_accept", &trig_Flag_trkPOG_toomanystripclus53X_accept, &b_trig_Flag_trkPOG_toomanystripclus53X_accept);
   fChain->SetBranchAddress("trig_Flag_trkPOG_logErrorTooManyClusters_accept", &trig_Flag_trkPOG_logErrorTooManyClusters_accept, &b_trig_Flag_trkPOG_logErrorTooManyClusters_accept);
   fChain->SetBranchAddress("trig_Flag_METFilters_accept", &trig_Flag_METFilters_accept, &b_trig_Flag_METFilters_accept);
   fChain->SetBranchAddress("trig_raw2digi_step_accept", &trig_raw2digi_step_accept, &b_trig_raw2digi_step_accept);
//   fChain->SetBranchAddress("trig_L1Reco_step_accept", &trig_L1Reco_step_accept, &b_trig_L1Reco_step_accept);
   fChain->SetBranchAddress("trig_reconstruction_step_accept", &trig_reconstruction_step_accept, &b_trig_reconstruction_step_accept);
   fChain->SetBranchAddress("trig_recosim_step_accept", &trig_recosim_step_accept, &b_trig_recosim_step_accept);
   fChain->SetBranchAddress("trig_eventinterpretaion_step_accept", &trig_eventinterpretaion_step_accept, &b_trig_eventinterpretaion_step_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle25_CaloIdL_MW_accept", &trig_HLT_DoubleEle25_CaloIdL_MW_accept, &b_trig_HLT_DoubleEle25_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle27_CaloIdL_MW_accept", &trig_HLT_DoubleEle27_CaloIdL_MW_accept, &b_trig_HLT_DoubleEle27_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle33_CaloIdL_MW_accept", &trig_HLT_DoubleEle33_CaloIdL_MW_accept, &b_trig_HLT_DoubleEle33_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept", &trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept, &b_trig_HLT_DoubleEle24_eta2p1_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele27_Ele37_CaloIdL_MW_accept", &trig_HLT_Ele27_Ele37_CaloIdL_MW_accept, &b_trig_HLT_Ele27_Ele37_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_Mu27_Ele37_CaloIdL_MW_accept", &trig_HLT_Mu27_Ele37_CaloIdL_MW_accept, &b_trig_HLT_Mu27_Ele37_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_Mu37_Ele27_CaloIdL_MW_accept", &trig_HLT_Mu37_Ele27_CaloIdL_MW_accept, &b_trig_HLT_Mu37_Ele27_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_HLT_Mu37_TkMu27_accept", &trig_HLT_Mu37_TkMu27_accept, &b_trig_HLT_Mu37_TkMu27_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu4_3_Bs_accept", &trig_HLT_DoubleMu4_3_Bs_accept, &b_trig_HLT_DoubleMu4_3_Bs_accept);
   fChain->SetBranchAddress("trig_HLT_Ele20_WPTight_Gsf_accept", &trig_HLT_Ele20_WPTight_Gsf_accept, &b_trig_HLT_Ele20_WPTight_Gsf_accept);
//   fChain->SetBranchAddress("trig_HLT_Ele15_WPLoose_Gsf_accept", &trig_HLT_Ele15_WPLoose_Gsf_accept, &b_trig_HLT_Ele15_WPLoose_Gsf_accept);
//   fChain->SetBranchAddress("trig_HLT_Ele17_WPLoose_Gsf_accept", &trig_HLT_Ele17_WPLoose_Gsf_accept, &b_trig_HLT_Ele17_WPLoose_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele20_WPLoose_Gsf_accept", &trig_HLT_Ele20_WPLoose_Gsf_accept, &b_trig_HLT_Ele20_WPLoose_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept", &trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept, &b_trig_HLT_Ele20_eta2p1_WPLoose_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept", &trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept, &b_trig_HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_accept);
   fChain->SetBranchAddress("trig_HLT_Ele27_WPTight_Gsf_accept", &trig_HLT_Ele27_WPTight_Gsf_accept, &b_trig_HLT_Ele27_WPTight_Gsf_accept);
//   fChain->SetBranchAddress("trig_HLT_Ele28_WPTight_Gsf_accept", &trig_HLT_Ele28_WPTight_Gsf_accept, &b_trig_HLT_Ele28_WPTight_Gsf_accept);
//   fChain->SetBranchAddress("trig_HLT_Ele30_WPTight_Gsf_accept", &trig_HLT_Ele30_WPTight_Gsf_accept, &b_trig_HLT_Ele30_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_accept", &trig_HLT_Ele32_WPTight_Gsf_accept, &b_trig_HLT_Ele32_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_accept", &trig_HLT_Ele35_WPTight_Gsf_accept, &b_trig_HLT_Ele35_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept", &trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept, &b_trig_HLT_Ele35_WPTight_Gsf_L1EGMT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele38_WPTight_Gsf_accept", &trig_HLT_Ele38_WPTight_Gsf_accept, &b_trig_HLT_Ele38_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele40_WPTight_Gsf_accept", &trig_HLT_Ele40_WPTight_Gsf_accept, &b_trig_HLT_Ele40_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept", &trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept, &b_trig_HLT_Ele32_WPTight_Gsf_L1DoubleEG_accept);
   fChain->SetBranchAddress("trig_HLT_IsoMu20_accept", &trig_HLT_IsoMu20_accept, &b_trig_HLT_IsoMu20_accept);
   fChain->SetBranchAddress("trig_HLT_IsoMu24_accept", &trig_HLT_IsoMu24_accept, &b_trig_HLT_IsoMu24_accept);
//   fChain->SetBranchAddress("trig_HLT_IsoTkMu24_accept", &trig_HLT_IsoTkMu24_accept, &b_trig_HLT_IsoTkMu24_accept);
   fChain->SetBranchAddress("trig_HLT_IsoMu24_eta2p1_accept", &trig_HLT_IsoMu24_eta2p1_accept, &b_trig_HLT_IsoMu24_eta2p1_accept);
   fChain->SetBranchAddress("trig_HLT_IsoMu27_accept", &trig_HLT_IsoMu27_accept, &b_trig_HLT_IsoMu27_accept);
   fChain->SetBranchAddress("trig_HLT_IsoMu30_accept", &trig_HLT_IsoMu30_accept, &b_trig_HLT_IsoMu30_accept);
   fChain->SetBranchAddress("trig_HLT_L1SingleMu18_accept", &trig_HLT_L1SingleMu18_accept, &b_trig_HLT_L1SingleMu18_accept);
   fChain->SetBranchAddress("trig_HLT_L1SingleMu25_accept", &trig_HLT_L1SingleMu25_accept, &b_trig_HLT_L1SingleMu25_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu10_accept", &trig_HLT_L2Mu10_accept, &b_trig_HLT_L2Mu10_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept", &trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept, &b_trig_HLT_L2Mu10_NoVertex_NoBPTX3BX_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu10_NoVertex_NoBPTX_accept", &trig_HLT_L2Mu10_NoVertex_NoBPTX_accept, &b_trig_HLT_L2Mu10_NoVertex_NoBPTX_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept", &trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept, &b_trig_HLT_L2Mu45_NoVertex_3Sta_NoBPTX3BX_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept", &trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept, &b_trig_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_accept);
   fChain->SetBranchAddress("trig_HLT_L2Mu50_accept", &trig_HLT_L2Mu50_accept, &b_trig_HLT_L2Mu50_accept);
//   fChain->SetBranchAddress("trig_HLT_L2Mu23NoVtx_2Cha_accept", &trig_HLT_L2Mu23NoVtx_2Cha_accept, &b_trig_HLT_L2Mu23NoVtx_2Cha_accept);
//   fChain->SetBranchAddress("trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept", &trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept, &b_trig_HLT_L2Mu23NoVtx_2Cha_CosmicSeed_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept", &trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu30NoVtx_2Cha_CosmicSeed_Eta2p4_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept", &trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu30NoVtx_2Cha_Eta2p4_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu50_accept", &trig_HLT_DoubleL2Mu50_accept, &b_trig_HLT_DoubleL2Mu50_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept", &trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept, &b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_accept", &trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_accept, &b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_NoL2Matched_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_NoL2Matched_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_CosmicSeed_Eta2p4_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept", &trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept, &b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_accept", &trig_HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_accept, &b_trig_HLT_DoubleL2Mu23NoVtx_2Cha_NoL2Matched_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_NoL2Matched_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept", &trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept, &b_trig_HLT_DoubleL2Mu25NoVtx_2Cha_Eta2p4_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept", &trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept, &b_trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept", &trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept, &b_trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept", &trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept, &b_trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept", &trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept, &b_trig_HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_accept);
   fChain->SetBranchAddress("trig_HLT_Mu25_TkMu0_Onia_accept", &trig_HLT_Mu25_TkMu0_Onia_accept, &b_trig_HLT_Mu25_TkMu0_Onia_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu30_TkMu0_Psi_accept", &trig_HLT_Mu30_TkMu0_Psi_accept, &b_trig_HLT_Mu30_TkMu0_Psi_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_TkMu0_Phi_accept", &trig_HLT_Mu20_TkMu0_Phi_accept, &b_trig_HLT_Mu20_TkMu0_Phi_accept);
   fChain->SetBranchAddress("trig_HLT_Mu25_TkMu0_Phi_accept", &trig_HLT_Mu25_TkMu0_Phi_accept, &b_trig_HLT_Mu25_TkMu0_Phi_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu12_accept", &trig_HLT_Mu12_accept, &b_trig_HLT_Mu12_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu15_accept", &trig_HLT_Mu15_accept, &b_trig_HLT_Mu15_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_accept", &trig_HLT_Mu20_accept, &b_trig_HLT_Mu20_accept);
   fChain->SetBranchAddress("trig_HLT_Mu27_accept", &trig_HLT_Mu27_accept, &b_trig_HLT_Mu27_accept);
   fChain->SetBranchAddress("trig_HLT_Mu50_accept", &trig_HLT_Mu50_accept, &b_trig_HLT_Mu50_accept);
   fChain->SetBranchAddress("trig_HLT_Mu55_accept", &trig_HLT_Mu55_accept, &b_trig_HLT_Mu55_accept);
   fChain->SetBranchAddress("trig_HLT_OldMu100_accept", &trig_HLT_OldMu100_accept, &b_trig_HLT_OldMu100_accept);
   fChain->SetBranchAddress("trig_HLT_TkMu100_accept", &trig_HLT_TkMu100_accept, &b_trig_HLT_TkMu100_accept);
   fChain->SetBranchAddress("trig_HLT_PFHT500_PFMET100_PFMHT100_IDTight_accept", &trig_HLT_PFHT500_PFMET100_PFMHT100_IDTight_accept, &b_trig_HLT_PFHT500_PFMET100_PFMHT100_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFHT500_PFMET110_PFMHT110_IDTight_accept", &trig_HLT_PFHT500_PFMET110_PFMHT110_IDTight_accept, &b_trig_HLT_PFHT500_PFMET110_PFMHT110_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFHT700_PFMET85_PFMHT85_IDTight_accept", &trig_HLT_PFHT700_PFMET85_PFMHT85_IDTight_accept, &b_trig_HLT_PFHT700_PFMET85_PFMHT85_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFHT700_PFMET95_PFMHT95_IDTight_accept", &trig_HLT_PFHT700_PFMET95_PFMHT95_IDTight_accept, &b_trig_HLT_PFHT700_PFMET95_PFMHT95_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFHT800_PFMET75_PFMHT75_IDTight_accept", &trig_HLT_PFHT800_PFMET75_PFMHT75_IDTight_accept, &b_trig_HLT_PFHT800_PFMET75_PFMHT75_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFHT800_PFMET85_PFMHT85_IDTight_accept", &trig_HLT_PFHT800_PFMET85_PFMHT85_IDTight_accept, &b_trig_HLT_PFHT800_PFMET85_PFMHT85_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET110_PFMHT110_IDTight_accept", &trig_HLT_PFMET110_PFMHT110_IDTight_accept, &b_trig_HLT_PFMET110_PFMHT110_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET120_PFMHT120_IDTight_accept", &trig_HLT_PFMET120_PFMHT120_IDTight_accept, &b_trig_HLT_PFMET120_PFMHT120_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET130_PFMHT130_IDTight_accept", &trig_HLT_PFMET130_PFMHT130_IDTight_accept, &b_trig_HLT_PFMET130_PFMHT130_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET140_PFMHT140_IDTight_accept", &trig_HLT_PFMET140_PFMHT140_IDTight_accept, &b_trig_HLT_PFMET140_PFMHT140_IDTight_accept);
//   fChain->SetBranchAddress("trig_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_accept", &trig_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_accept, &b_trig_HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_accept);
//   fChain->SetBranchAddress("trig_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_accept", &trig_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_accept, &b_trig_HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_accept);
//   fChain->SetBranchAddress("trig_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_accept", &trig_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_accept, &b_trig_HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_accept);
//   fChain->SetBranchAddress("trig_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_accept", &trig_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_accept, &b_trig_HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_accept);
//   fChain->SetBranchAddress("trig_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_accept", &trig_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_accept, &b_trig_HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET120_PFMHT120_IDTight_PFHT60_accept", &trig_HLT_PFMET120_PFMHT120_IDTight_PFHT60_accept, &b_trig_HLT_PFMET120_PFMHT120_IDTight_PFHT60_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_accept", &trig_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_accept, &b_trig_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_accept", &trig_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_accept, &b_trig_HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETTypeOne110_PFMHT110_IDTight_accept", &trig_HLT_PFMETTypeOne110_PFMHT110_IDTight_accept, &b_trig_HLT_PFMETTypeOne110_PFMHT110_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETTypeOne120_PFMHT120_IDTight_accept", &trig_HLT_PFMETTypeOne120_PFMHT120_IDTight_accept, &b_trig_HLT_PFMETTypeOne120_PFMHT120_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETTypeOne130_PFMHT130_IDTight_accept", &trig_HLT_PFMETTypeOne130_PFMHT130_IDTight_accept, &b_trig_HLT_PFMETTypeOne130_PFMHT130_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETTypeOne140_PFMHT140_IDTight_accept", &trig_HLT_PFMETTypeOne140_PFMHT140_IDTight_accept, &b_trig_HLT_PFMETTypeOne140_PFMHT140_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_accept", &trig_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_accept, &b_trig_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_accept", &trig_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_accept, &b_trig_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_accept", &trig_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_accept, &b_trig_HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_accept", &trig_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_accept, &b_trig_HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_accept", &trig_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_accept, &b_trig_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_accept", &trig_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_accept, &b_trig_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_accept", &trig_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_accept, &b_trig_HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_accept", &trig_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_accept, &b_trig_HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET80_NotCleaned_accept", &trig_HLT_CaloMET80_NotCleaned_accept, &b_trig_HLT_CaloMET80_NotCleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET90_NotCleaned_accept", &trig_HLT_CaloMET90_NotCleaned_accept, &b_trig_HLT_CaloMET90_NotCleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET100_NotCleaned_accept", &trig_HLT_CaloMET100_NotCleaned_accept, &b_trig_HLT_CaloMET100_NotCleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET110_NotCleaned_accept", &trig_HLT_CaloMET110_NotCleaned_accept, &b_trig_HLT_CaloMET110_NotCleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET250_NotCleaned_accept", &trig_HLT_CaloMET250_NotCleaned_accept, &b_trig_HLT_CaloMET250_NotCleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET70_HBHECleaned_accept", &trig_HLT_CaloMET70_HBHECleaned_accept, &b_trig_HLT_CaloMET70_HBHECleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET80_HBHECleaned_accept", &trig_HLT_CaloMET80_HBHECleaned_accept, &b_trig_HLT_CaloMET80_HBHECleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET90_HBHECleaned_accept", &trig_HLT_CaloMET90_HBHECleaned_accept, &b_trig_HLT_CaloMET90_HBHECleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET100_HBHECleaned_accept", &trig_HLT_CaloMET100_HBHECleaned_accept, &b_trig_HLT_CaloMET100_HBHECleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET250_HBHECleaned_accept", &trig_HLT_CaloMET250_HBHECleaned_accept, &b_trig_HLT_CaloMET250_HBHECleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET300_HBHECleaned_accept", &trig_HLT_CaloMET300_HBHECleaned_accept, &b_trig_HLT_CaloMET300_HBHECleaned_accept);
   fChain->SetBranchAddress("trig_HLT_CaloMET350_HBHECleaned_accept", &trig_HLT_CaloMET350_HBHECleaned_accept, &b_trig_HLT_CaloMET350_HBHECleaned_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET200_NotCleaned_accept", &trig_HLT_PFMET200_NotCleaned_accept, &b_trig_HLT_PFMET200_NotCleaned_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET200_HBHECleaned_accept", &trig_HLT_PFMET200_HBHECleaned_accept, &b_trig_HLT_PFMET200_HBHECleaned_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET250_HBHECleaned_accept", &trig_HLT_PFMET250_HBHECleaned_accept, &b_trig_HLT_PFMET250_HBHECleaned_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET300_HBHECleaned_accept", &trig_HLT_PFMET300_HBHECleaned_accept, &b_trig_HLT_PFMET300_HBHECleaned_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET200_HBHE_BeamHaloCleaned_accept", &trig_HLT_PFMET200_HBHE_BeamHaloCleaned_accept, &b_trig_HLT_PFMET200_HBHE_BeamHaloCleaned_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_accept", &trig_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_accept, &b_trig_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_accept);
   fChain->SetBranchAddress("trig_HLT_MET105_IsoTrk50_accept", &trig_HLT_MET105_IsoTrk50_accept, &b_trig_HLT_MET105_IsoTrk50_accept);
   fChain->SetBranchAddress("trig_HLT_MET120_IsoTrk50_accept", &trig_HLT_MET120_IsoTrk50_accept, &b_trig_HLT_MET120_IsoTrk50_accept);
   fChain->SetBranchAddress("trig_HLT_Photon300_NoHE_accept", &trig_HLT_Photon300_NoHE_accept, &b_trig_HLT_Photon300_NoHE_accept);
   fChain->SetBranchAddress("trig_HLT_Mu8_TrkIsoVVL_accept", &trig_HLT_Mu8_TrkIsoVVL_accept, &b_trig_HLT_Mu8_TrkIsoVVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_TrkIsoVVL_accept", &trig_HLT_Mu17_TrkIsoVVL_accept, &b_trig_HLT_Mu17_TrkIsoVVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_TrkIsoVVL_accept", &trig_HLT_Mu19_TrkIsoVVL_accept, &b_trig_HLT_Mu19_TrkIsoVVL_accept);
//   fChain->SetBranchAddress("trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_accept);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept", &trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept, &b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept);
//   fChain->SetBranchAddress("trig_HLT_Photon20_accept", &trig_HLT_Photon20_accept, &b_trig_HLT_Photon20_accept);
   fChain->SetBranchAddress("trig_HLT_Photon33_accept", &trig_HLT_Photon33_accept, &b_trig_HLT_Photon33_accept);
   fChain->SetBranchAddress("trig_HLT_Photon50_accept", &trig_HLT_Photon50_accept, &b_trig_HLT_Photon50_accept);
   fChain->SetBranchAddress("trig_HLT_Photon75_accept", &trig_HLT_Photon75_accept, &b_trig_HLT_Photon75_accept);
   fChain->SetBranchAddress("trig_HLT_Photon90_accept", &trig_HLT_Photon90_accept, &b_trig_HLT_Photon90_accept);
   fChain->SetBranchAddress("trig_HLT_Photon120_accept", &trig_HLT_Photon120_accept, &b_trig_HLT_Photon120_accept);
   fChain->SetBranchAddress("trig_HLT_Photon150_accept", &trig_HLT_Photon150_accept, &b_trig_HLT_Photon150_accept);
   fChain->SetBranchAddress("trig_HLT_Photon175_accept", &trig_HLT_Photon175_accept, &b_trig_HLT_Photon175_accept);
   fChain->SetBranchAddress("trig_HLT_Photon200_accept", &trig_HLT_Photon200_accept, &b_trig_HLT_Photon200_accept);
//   fChain->SetBranchAddress("trig_HLT_Photon100EB_TightID_TightIso_accept", &trig_HLT_Photon100EB_TightID_TightIso_accept, &b_trig_HLT_Photon100EB_TightID_TightIso_accept);
//   fChain->SetBranchAddress("trig_HLT_Photon110EB_TightID_TightIso_accept", &trig_HLT_Photon110EB_TightID_TightIso_accept, &b_trig_HLT_Photon110EB_TightID_TightIso_accept);
//   fChain->SetBranchAddress("trig_HLT_Photon120EB_TightID_TightIso_accept", &trig_HLT_Photon120EB_TightID_TightIso_accept, &b_trig_HLT_Photon120EB_TightID_TightIso_accept);
//   fChain->SetBranchAddress("trig_HLT_Photon100EBHE10_accept", &trig_HLT_Photon100EBHE10_accept, &b_trig_HLT_Photon100EBHE10_accept);
//   fChain->SetBranchAddress("trig_HLT_Photon100EEHE10_accept", &trig_HLT_Photon100EEHE10_accept, &b_trig_HLT_Photon100EEHE10_accept);
//   fChain->SetBranchAddress("trig_HLT_Photon100EE_TightID_TightIso_accept", &trig_HLT_Photon100EE_TightID_TightIso_accept, &b_trig_HLT_Photon100EE_TightID_TightIso_accept);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon50_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_accept);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon75_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_accept);
//   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_accept", &trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_accept, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3_accept);
//   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_accept", &trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_accept, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3_accept);
   fChain->SetBranchAddress("trig_HLT_Photon90_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon90_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon90_R9Id90_HE10_IsoM_accept);
   fChain->SetBranchAddress("trig_HLT_Photon120_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon120_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon120_R9Id90_HE10_IsoM_accept);
   fChain->SetBranchAddress("trig_HLT_Photon165_R9Id90_HE10_IsoM_accept", &trig_HLT_Photon165_R9Id90_HE10_IsoM_accept, &b_trig_HLT_Photon165_R9Id90_HE10_IsoM_accept);
   fChain->SetBranchAddress("trig_HLT_Photon90_CaloIdL_PFHT700_accept", &trig_HLT_Photon90_CaloIdL_PFHT700_accept, &b_trig_HLT_Photon90_CaloIdL_PFHT700_accept);
//   fChain->SetBranchAddress("trig_HLT_Photon35_TwoProngs35_accept", &trig_HLT_Photon35_TwoProngs35_accept, &b_trig_HLT_Photon35_TwoProngs35_accept);
//   fChain->SetBranchAddress("trig_HLT_IsoMu24_TwoProngs35_accept", &trig_HLT_IsoMu24_TwoProngs35_accept, &b_trig_HLT_IsoMu24_TwoProngs35_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_accept", &trig_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_accept, &b_trig_HLT_DoubleMu3_DZ_PFMET50_PFMHT60_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu3_DZ_PFMET70_PFMHT70_accept", &trig_HLT_DoubleMu3_DZ_PFMET70_PFMHT70_accept, &b_trig_HLT_DoubleMu3_DZ_PFMET70_PFMHT70_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu3_DZ_PFMET90_PFMHT90_accept", &trig_HLT_DoubleMu3_DZ_PFMET90_PFMHT90_accept, &b_trig_HLT_DoubleMu3_DZ_PFMET90_PFMHT90_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu43NoFiltersNoVtx_accept", &trig_HLT_DoubleMu43NoFiltersNoVtx_accept, &b_trig_HLT_DoubleMu43NoFiltersNoVtx_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu48NoFiltersNoVtx_accept", &trig_HLT_DoubleMu48NoFiltersNoVtx_accept, &b_trig_HLT_DoubleMu48NoFiltersNoVtx_accept);
   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept", &trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept, &b_trig_HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_accept);
   fChain->SetBranchAddress("trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept", &trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept, &b_trig_HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_accept", &trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_accept, &b_trig_HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_accept", &trig_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_accept, &b_trig_HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleMu33NoFiltersNoVtxDisplaced_accept", &trig_HLT_DoubleMu33NoFiltersNoVtxDisplaced_accept, &b_trig_HLT_DoubleMu33NoFiltersNoVtxDisplaced_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleMu40NoFiltersNoVtxDisplaced_accept", &trig_HLT_DoubleMu40NoFiltersNoVtxDisplaced_accept, &b_trig_HLT_DoubleMu40NoFiltersNoVtxDisplaced_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept", &trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept, &b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept", &trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept, &b_trig_HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_accept);
   fChain->SetBranchAddress("trig_HLT_DiJet110_35_Mjj650_PFMET110_accept", &trig_HLT_DiJet110_35_Mjj650_PFMET110_accept, &b_trig_HLT_DiJet110_35_Mjj650_PFMET110_accept);
   fChain->SetBranchAddress("trig_HLT_DiJet110_35_Mjj650_PFMET120_accept", &trig_HLT_DiJet110_35_Mjj650_PFMET120_accept, &b_trig_HLT_DiJet110_35_Mjj650_PFMET120_accept);
   fChain->SetBranchAddress("trig_HLT_DiJet110_35_Mjj650_PFMET130_accept", &trig_HLT_DiJet110_35_Mjj650_PFMET130_accept, &b_trig_HLT_DiJet110_35_Mjj650_PFMET130_accept);
   fChain->SetBranchAddress("trig_HLT_TripleJet110_35_35_Mjj650_PFMET110_accept", &trig_HLT_TripleJet110_35_35_Mjj650_PFMET110_accept, &b_trig_HLT_TripleJet110_35_35_Mjj650_PFMET110_accept);
   fChain->SetBranchAddress("trig_HLT_TripleJet110_35_35_Mjj650_PFMET120_accept", &trig_HLT_TripleJet110_35_35_Mjj650_PFMET120_accept, &b_trig_HLT_TripleJet110_35_35_Mjj650_PFMET120_accept);
   fChain->SetBranchAddress("trig_HLT_TripleJet110_35_35_Mjj650_PFMET130_accept", &trig_HLT_TripleJet110_35_35_Mjj650_PFMET130_accept, &b_trig_HLT_TripleJet110_35_35_Mjj650_PFMET130_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept", &trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept, &b_trig_HLT_DoubleMu20_7_Mass0to30_Photon23_accept);
   fChain->SetBranchAddress("trig_HLT_Ele15_IsoVVVL_PFHT450_PFMET50_accept", &trig_HLT_Ele15_IsoVVVL_PFHT450_PFMET50_accept, &b_trig_HLT_Ele15_IsoVVVL_PFHT450_PFMET50_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_accept", &trig_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_accept, &b_trig_HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_accept);
   fChain->SetBranchAddress("trig_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_accept", &trig_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_accept, &b_trig_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_accept);
   fChain->SetBranchAddress("trig_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_accept", &trig_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_accept, &b_trig_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_accept);
   fChain->SetBranchAddress("trig_HLT_Mu15_IsoVVVL_PFHT450_PFMET50_accept", &trig_HLT_Mu15_IsoVVVL_PFHT450_PFMET50_accept, &b_trig_HLT_Mu15_IsoVVVL_PFHT450_PFMET50_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_accept", &trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_accept, &b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_accept", &trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_accept, &b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_accept", &trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_accept, &b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_accept", &trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_accept, &b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_accept", &trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_accept, &b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_accept", &trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_accept, &b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_accept", &trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_accept, &b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_accept);
//   fChain->SetBranchAddress("trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_accept", &trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_accept, &b_trig_HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_accept);
   fChain->SetBranchAddress("trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept", &trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept, &b_trig_HLT_Dimuon14_Phi_Barrel_Seagulls_accept);
   fChain->SetBranchAddress("trig_HLT_Dimuon24_Phi_noCorrL1_accept", &trig_HLT_Dimuon24_Phi_noCorrL1_accept, &b_trig_HLT_Dimuon24_Phi_noCorrL1_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleIsoMu20_eta2p1_accept", &trig_HLT_DoubleIsoMu20_eta2p1_accept, &b_trig_HLT_DoubleIsoMu20_eta2p1_accept);
   fChain->SetBranchAddress("trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept", &trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_accept);
   fChain->SetBranchAddress("trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept", &trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu16_DoubleTrkMu6NoFiltersNoVtx_accept);
   fChain->SetBranchAddress("trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept", &trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx_accept);
   fChain->SetBranchAddress("trig_HLT_Mu8_accept", &trig_HLT_Mu8_accept, &b_trig_HLT_Mu8_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_accept", &trig_HLT_Mu17_accept, &b_trig_HLT_Mu17_accept);
   fChain->SetBranchAddress("trig_HLT_Mu19_accept", &trig_HLT_Mu19_accept, &b_trig_HLT_Mu19_accept);
   fChain->SetBranchAddress("trig_HLT_Mu17_Photon30_IsoCaloId_accept", &trig_HLT_Mu17_Photon30_IsoCaloId_accept, &b_trig_HLT_Mu17_Photon30_IsoCaloId_accept);
   fChain->SetBranchAddress("trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele135_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele145_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele200_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele250_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept", &trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept, &b_trig_HLT_Ele300_CaloIdVT_GsfTrkIdT_accept);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_accept);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_accept);
   fChain->SetBranchAddress("trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_accept", &trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_accept, &b_trig_HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350MinPFJet15_accept);
//   fChain->SetBranchAddress("trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept", &trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept, &b_trig_HLT_SinglePhoton10_Eta3p1ForPPRef_accept);
//   fChain->SetBranchAddress("trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept", &trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept, &b_trig_HLT_SinglePhoton20_Eta3p1ForPPRef_accept);
//   fChain->SetBranchAddress("trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept", &trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept, &b_trig_HLT_SinglePhoton30_Eta3p1ForPPRef_accept);
   fChain->SetBranchAddress("trig_HLT_Photon20_HoverELoose_accept", &trig_HLT_Photon20_HoverELoose_accept, &b_trig_HLT_Photon20_HoverELoose_accept);
   fChain->SetBranchAddress("trig_HLT_Photon30_HoverELoose_accept", &trig_HLT_Photon30_HoverELoose_accept, &b_trig_HLT_Photon30_HoverELoose_accept);
   fChain->SetBranchAddress("trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept", &trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept, &b_trig_HLT_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142_accept);
   fChain->SetBranchAddress("trig_AlCa_RPCMuonNormalisation_accept", &trig_AlCa_RPCMuonNormalisation_accept, &b_trig_AlCa_RPCMuonNormalisation_accept);
   fChain->SetBranchAddress("trig_MC_PFMET_accept", &trig_MC_PFMET_accept, &b_trig_MC_PFMET_accept);
   fChain->SetBranchAddress("trig_MC_CaloMET_accept", &trig_MC_CaloMET_accept, &b_trig_MC_CaloMET_accept);
   fChain->SetBranchAddress("trig_MC_CaloMET_JetIdCleaned_accept", &trig_MC_CaloMET_JetIdCleaned_accept, &b_trig_MC_CaloMET_JetIdCleaned_accept);
   fChain->SetBranchAddress("trig_MC_DoubleEle5_CaloIdL_MW_accept", &trig_MC_DoubleEle5_CaloIdL_MW_accept, &b_trig_MC_DoubleEle5_CaloIdL_MW_accept);
   fChain->SetBranchAddress("trig_MC_Ele5_WPTight_Gsf_accept", &trig_MC_Ele5_WPTight_Gsf_accept, &b_trig_MC_Ele5_WPTight_Gsf_accept);
   fChain->SetBranchAddress("trig_MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_accept", &trig_MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_accept, &b_trig_MC_Ele15_Ele10_CaloIdL_TrackIdL_IsoVL_DZ_accept);
   fChain->SetBranchAddress("trig_MC_IsoMu_accept", &trig_MC_IsoMu_accept, &b_trig_MC_IsoMu_accept);
   fChain->SetBranchAddress("trig_MC_DoubleMu_TrkIsoVVL_DZ_accept", &trig_MC_DoubleMu_TrkIsoVVL_DZ_accept, &b_trig_MC_DoubleMu_TrkIsoVVL_DZ_accept);
   fChain->SetBranchAddress("trig_MC_DoubleMuNoFiltersNoVtx_accept", &trig_MC_DoubleMuNoFiltersNoVtx_accept, &b_trig_MC_DoubleMuNoFiltersNoVtx_accept);
//   fChain->SetBranchAddress("trig_HLT_IsoMu27_MET90_accept", &trig_HLT_IsoMu27_MET90_accept, &b_trig_HLT_IsoMu27_MET90_accept);
   fChain->SetBranchAddress("trig_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_accept", &trig_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_accept, &b_trig_HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50_accept);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_accept", &trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_accept, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_accept);
   fChain->SetBranchAddress("trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_accept", &trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_accept, &b_trig_HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3_accept);
   fChain->SetBranchAddress("trig_HLT_PFMET100_PFMHT100_IDTight_PFHT60_accept", &trig_HLT_PFMET100_PFMHT100_IDTight_PFHT60_accept, &b_trig_HLT_PFMET100_PFMHT100_IDTight_PFHT60_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_accept", &trig_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_accept, &b_trig_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_accept);
   fChain->SetBranchAddress("trig_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_accept", &trig_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_accept, &b_trig_HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_accept);
   fChain->SetBranchAddress("trig_HLT_Mu18_Mu9_SameSign_accept", &trig_HLT_Mu18_Mu9_SameSign_accept, &b_trig_HLT_Mu18_Mu9_SameSign_accept);
   fChain->SetBranchAddress("trig_HLT_Mu18_Mu9_SameSign_DZ_accept", &trig_HLT_Mu18_Mu9_SameSign_DZ_accept, &b_trig_HLT_Mu18_Mu9_SameSign_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu18_Mu9_accept", &trig_HLT_Mu18_Mu9_accept, &b_trig_HLT_Mu18_Mu9_accept);
   fChain->SetBranchAddress("trig_HLT_Mu18_Mu9_DZ_accept", &trig_HLT_Mu18_Mu9_DZ_accept, &b_trig_HLT_Mu18_Mu9_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_Mu10_SameSign_accept", &trig_HLT_Mu20_Mu10_SameSign_accept, &b_trig_HLT_Mu20_Mu10_SameSign_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_Mu10_SameSign_DZ_accept", &trig_HLT_Mu20_Mu10_SameSign_DZ_accept, &b_trig_HLT_Mu20_Mu10_SameSign_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_Mu10_accept", &trig_HLT_Mu20_Mu10_accept, &b_trig_HLT_Mu20_Mu10_accept);
   fChain->SetBranchAddress("trig_HLT_Mu20_Mu10_DZ_accept", &trig_HLT_Mu20_Mu10_DZ_accept, &b_trig_HLT_Mu20_Mu10_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_Mu12_SameSign_accept", &trig_HLT_Mu23_Mu12_SameSign_accept, &b_trig_HLT_Mu23_Mu12_SameSign_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_Mu12_SameSign_DZ_accept", &trig_HLT_Mu23_Mu12_SameSign_DZ_accept, &b_trig_HLT_Mu23_Mu12_SameSign_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_Mu12_accept", &trig_HLT_Mu23_Mu12_accept, &b_trig_HLT_Mu23_Mu12_accept);
   fChain->SetBranchAddress("trig_HLT_Mu23_Mu12_DZ_accept", &trig_HLT_Mu23_Mu12_DZ_accept, &b_trig_HLT_Mu23_Mu12_DZ_accept);
   fChain->SetBranchAddress("trig_HLT_DoubleMu3_DCA_PFMET50_PFMHT60_accept", &trig_HLT_DoubleMu3_DCA_PFMET50_PFMHT60_accept, &b_trig_HLT_DoubleMu3_DCA_PFMET50_PFMHT60_accept);
//   fChain->SetBranchAddress("trig_HLT_TrkMu6NoFiltersNoVtx_accept", &trig_HLT_TrkMu6NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu6NoFiltersNoVtx_accept);
//   fChain->SetBranchAddress("trig_HLT_TrkMu16NoFiltersNoVtx_accept", &trig_HLT_TrkMu16NoFiltersNoVtx_accept, &b_trig_HLT_TrkMu16NoFiltersNoVtx_accept);
//   fChain->SetBranchAddress("trig_HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_accept", &trig_HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_accept, &b_trig_HLT_DoubleTrkMu_16_6_NoFiltersNoVtx_accept);
   Notify();
}

Bool_t MyAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyAnalysis_cxx
