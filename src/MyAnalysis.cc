#define MyAnalysis_cxx
#include "../include/MyAnalysis.h"
#include "../include/PU_reWeighting.h"
#include "../include/lepton_candidate.h"
#include "../include/jet_candidate.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TEntryList.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <vector>
#include "../include/RoccoR.h"
#include "../include/BTagCalibrationStandalone.h"

void displayProgress(long current, long max){
  using std::cerr;
  if (max<2500) return;
  if (current%(max/2500)!=0 && current<max-1) return;

  int width = 52; // Hope the terminal is at least that wide.
  int barWidth = width - 2;
  cerr << "\x1B[2K"; // Clear line
  cerr << "\x1B[2000D"; // Cursor left
  cerr << '[';
  for(int i=0 ; i<barWidth ; ++i){ if(i<barWidth*current/max){ cerr << '=' ; }else{ cerr << ' ' ; } }
  cerr << ']';
  cerr << " " << Form("%8d/%8d (%5.2f%%)", (int)current, (int)max, 100.0*current/max) ;
  cerr.flush();
}

Double_t deltaPhi(Double_t phi1, Double_t phi2) {
  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
  return dPhi;
}


Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);
  return sqrt(dEta*dEta+dPhi*dPhi);
}

bool ComparePtLep(lepton_candidate *a, lepton_candidate *b) { return a->pt_ > b->pt_; }
bool CompareChargeLep(lepton_candidate *a, lepton_candidate *b) { return a->charge_ < b->charge_; }
bool ComparePtJet(jet_candidate *a, jet_candidate *b) { return a->pt_ > b->pt_; }

float scale_factor( TH2F* h, float X, float Y , TString uncert){
  int NbinsX=h->GetXaxis()->GetNbins();
  int NbinsY=h->GetYaxis()->GetNbins();
  float x_min=h->GetXaxis()->GetBinLowEdge(1);
  float x_max=h->GetXaxis()->GetBinLowEdge(NbinsX)+h->GetXaxis()->GetBinWidth(NbinsX);
  float y_min=h->GetYaxis()->GetBinLowEdge(1);
  float y_max=h->GetYaxis()->GetBinLowEdge(NbinsY)+h->GetYaxis()->GetBinWidth(NbinsY);
  TAxis *Xaxis = h->GetXaxis();
  TAxis *Yaxis = h->GetYaxis();
  Int_t binx=1;
  Int_t biny=1;
  if(x_min < X && X < x_max) binx = Xaxis->FindBin(X);
  else binx= (X<=x_min) ? 1 : NbinsX ;
  if(y_min < Y && Y < y_max) biny = Yaxis->FindBin(Y);
  else biny= (Y<=y_min) ? 1 : NbinsY ;
  if(uncert=="up") return (h->GetBinContent(binx, biny)+h->GetBinError(binx, biny));
  else if(uncert=="down") return (h->GetBinContent(binx, biny)-h->GetBinError(binx, biny));
  else return  h->GetBinContent(binx, biny);
}

float topPt(float pt){
  return (0.973 - (0.000134 * pt) + (0.103 * exp(pt * (-0.0118))));  
}

void MyAnalysis::Loop(TString fname, TString data, TString dataset ,TString year, TString run, float xs, float lumi, float Nevent)
{

  typedef vector<TH2F*> Dim1;
  typedef vector<Dim1> Dim2;

  std::vector<TString> regions{"ll","llMetg20","llMetg20Jetgeq1Bleq1","llJetgeq2Bleq1","llMetg20Bgeq1","llJetgeq1B0","llMetl20Jetgeq1B0","llMetg20B2"};
  std::vector<TString> channels{"D_e","D_mu","N_e","N_mu"};
    
  Double_t ptBinsFF[9] = {20., 30., 40., 50., 60., 70., 80., 90., 100.};
  Double_t etaBinsFF[9] = {0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4};

  Dim2 HistsMR(channels.size(),Dim1(regions.size()));
  std::stringstream name;
  TH2F *h_test;
  for (int i=0;i<(int)channels.size();++i){
    for (int k=0;k<(int)regions.size();++k){
        name<<channels[i]<<"_"<<regions[k];
        h_test = new TH2F((name.str()).c_str(),(name.str()).c_str(),8,etaBinsFF,8,ptBinsFF);
        h_test->StatOverflows(kTRUE);
        h_test->Sumw2(kTRUE);
        HistsMR[i][k] = h_test;
        name.str("");
    }
  }

//Get scale factor and weight histograms
  TH2F  sf_Ele_Reco_H;
  TH2F  sf_Ele_ID_H;
  TH2F  sf_Mu_ID_H;
  TH2F  sf_Mu_ISO_H;
  TH2F  sf_triggeree_H;
  TH2F  sf_triggeremu_H;
  TH2F  sf_triggermumu_H;
  TH2F  btagEff_b_H;
  TH2F  btagEff_c_H;
  TH2F  btagEff_udsg_H;
  PU wPU;
  std::string rochesterFile;
  std::string btagFile;
  BTagCalibrationReader reader(BTagEntry::OP_MEDIUM, "central", {"up", "down"});

  std::string CMSSW_BASE(getenv("CMSSW_BASE"));
  RoccoR  rc;
  if(year == "2016")    rochesterFile = CMSSW_BASE + std::string("/src/data/TopLFV/input/RoccoR2016.txt");
  if(year == "2017")    rochesterFile = CMSSW_BASE + std::string("/src/data/TopLFV/input/RoccoR2017.txt");
  if(year == "2018")    rochesterFile = CMSSW_BASE + std::string("/src/data/TopLFV/input/RoccoR2016.txt");
  rc.init(rochesterFile);

  if(data == "mc"){
    TFile *f_btagEff_Map = new TFile(  ( CMSSW_BASE  + std::string("/src/data/TopLFV/input/btagEff.root") ).c_str()  );
    if(year == "2016"){
      btagEff_b_H = *(TH2F*)f_btagEff_Map->Get("2016_h2_BTaggingEff_b");
      btagEff_c_H = *(TH2F*)f_btagEff_Map->Get("2016_h2_BTaggingEff_c");
      btagEff_udsg_H = *(TH2F*)f_btagEff_Map->Get("2016_h2_BTaggingEff_udsg");
    }
    if(year == "2017"){
      btagEff_b_H = *(TH2F*)f_btagEff_Map->Get("2017_h2_BTaggingEff_b");
      btagEff_c_H = *(TH2F*)f_btagEff_Map->Get("2017_h2_BTaggingEff_c");
      btagEff_udsg_H = *(TH2F*)f_btagEff_Map->Get("2017_h2_BTaggingEff_udsg");
    }
    if(year == "2018"){
      btagEff_b_H = *(TH2F*)f_btagEff_Map->Get("2018_h2_BTaggingEff_b");
      btagEff_c_H = *(TH2F*)f_btagEff_Map->Get("2018_h2_BTaggingEff_c");
      btagEff_udsg_H = *(TH2F*)f_btagEff_Map->Get("2018_h2_BTaggingEff_udsg");
    }

    if(year == "2016")    btagFile = CMSSW_BASE + std::string("/src/data/TopLFV/input/DeepCSV_2016LegacySF_V1.csv");
    if(year == "2017")    btagFile = CMSSW_BASE + std::string("/src/data/TopLFV/input/DeepCSV_94XSF_V4_B_F.csv");
    if(year == "2018")    btagFile = CMSSW_BASE + std::string("/src/data/TopLFV/input/DeepCSV_102XSF_V1.csv");

    BTagCalibration calib("DeepCSV",btagFile);
    reader.load(calib,BTagEntry::FLAV_B,"comb"); 
    reader.load(calib,BTagEntry::FLAV_C,"comb");
    reader.load(calib,BTagEntry::FLAV_UDSG,"comb");

    if(year == "2016"){
      TFile *f_Ele_Reco_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root")).c_str() );
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2016LegacyReReco_ElectronTight_Fall17V2.root")).c_str() );
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map_1 = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2016_RunBCDEF_SF_ID.root")).c_str() ) ;
      TH2F *sf_Mu_ID_H_1 = (TH2F*)f_Mu_ID_Map_1->Get("NUM_TightID_DEN_genTracks_eta_pt");
      TFile *f_Mu_ID_Map_2 = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2016_RunGH_SF_ID.root")).c_str() ) ;
      TH2F *sf_Mu_ID_H_2 = (TH2F*)f_Mu_ID_Map_2->Get("NUM_TightID_DEN_genTracks_eta_pt");
      sf_Mu_ID_H_1->Scale(0.55);
      sf_Mu_ID_H_2->Scale(0.45);
      sf_Mu_ID_H_1->Add(sf_Mu_ID_H_2);
      sf_Mu_ID_H = *sf_Mu_ID_H_1;

      TFile *f_Mu_ISO_Map_1 = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2016_RunBCDEF_SF_ISO.root")).c_str() ) ;
      TH2F *sf_Mu_ISO_H_1 = (TH2F*)f_Mu_ISO_Map_1->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
      TFile *f_Mu_ISO_Map_2 = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2016_RunGH_SF_ISO.root")).c_str() ) ;
      TH2F *sf_Mu_ISO_H_2 = (TH2F*)f_Mu_ISO_Map_2->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
      sf_Mu_ISO_H_1->Scale(0.55);
      sf_Mu_ISO_H_2->Scale(0.45);
      sf_Mu_ISO_H_1->Add(sf_Mu_ISO_H_2);
      sf_Mu_ISO_H = *sf_Mu_ISO_H_1;

      TFile *f_triggeree = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/TriggerSF_ee2016_pt.root")).c_str() ) ;
      sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
      TFile *f_triggeremu = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/TriggerSF_emu2016_pt.root")).c_str() ) ;
      sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("h_lep1Pt_lep2Pt_Step3");
      TFile *f_triggermumu = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/TriggerSF_mumu2016_pt.root")).c_str() ) ;
      sf_triggermumu_H = *(TH2F*)f_triggermumu->Get("h_lep1Pt_lep2Pt_Step9");

      f_Ele_Reco_Map->Close();
      f_Ele_ID_Map->Close();
      f_Mu_ID_Map_1->Close();
      f_Mu_ID_Map_2->Close();
      f_Mu_ISO_Map_1->Close();
      f_Mu_ISO_Map_2->Close();
      f_triggeree->Close();
      f_triggeremu->Close();
      f_triggermumu->Close();
    }
    if(year == "2017"){
        TFile *f_Ele_Reco_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root")).c_str() ) ;
        sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

        TFile *f_Ele_ID_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2017_ElectronTight.root")).c_str() ) ;
        sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

        TFile *f_Mu_ID_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2017_RunBCDEF_SF_ID_syst.root")).c_str() ) ;
        sf_Mu_ID_H = *(TH2F*)f_Mu_ID_Map->Get("NUM_TightID_DEN_genTracks_pt_abseta");

        TFile *f_Mu_ISO_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2017_RunBCDEF_SF_ISO_syst.root")).c_str() ) ;
        sf_Mu_ISO_H = *(TH2F*)f_Mu_ISO_Map->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

        TFile *f_triggeree = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/TriggerSF_ee2017_pt.root")).c_str() ) ;
        sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
        TFile *f_triggeremu = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/TriggerSF_emu2017_pt.root")).c_str() ) ;
        sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("h_lep1Pt_lep2Pt_Step3");
        TFile *f_triggermumu = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/TriggerSF_mumu2017_pt.root")).c_str() ) ;
        sf_triggermumu_H = *(TH2F*)f_triggermumu->Get("h_lep1Pt_lep2Pt_Step9");

      f_Ele_Reco_Map->Close();
      f_Ele_ID_Map->Close();
      f_Mu_ID_Map->Close();
      f_Mu_ISO_Map->Close();
      f_triggeree->Close();
      f_triggeremu->Close();
      f_triggermumu->Close();
    }
    if(year == "2018"){
        TFile *f_Ele_Reco_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/egammaEffi.txt_EGM2D_updatedAll.root")).c_str() ) ;
        sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

        TFile *f_Ele_ID_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2018_ElectronTight.root")).c_str() ) ;
        sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

        TFile *f_Mu_ID_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2018_RunABCD_SF_ID.root")).c_str() ) ;
        sf_Mu_ID_H = *(TH2F*)f_Mu_ID_Map->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta");

        TFile *f_Mu_ISO_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2018_RunABCD_SF_ISO.root")).c_str() ) ;
        sf_Mu_ISO_H = *(TH2F*)f_Mu_ISO_Map->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

        TFile *f_triggeree = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/TriggerSF_ee2018_pt.root")).c_str() ) ;
        sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
        TFile *f_triggeremu = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/TriggerSF_emu2018_pt.root")).c_str() ) ;
        sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("h_lep1Pt_lep2Pt_Step3");
        TFile *f_triggermumu = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/TriggerSF_mumu2018_pt.root")).c_str() ) ;
        sf_triggermumu_H = *(TH2F*)f_triggermumu->Get("h_lep1Pt_lep2Pt_Step9");

      f_Ele_Reco_Map->Close();
      f_Ele_ID_Map->Close();
      f_Mu_ID_Map->Close();
      f_Mu_ISO_Map->Close();
      f_triggeree->Close();
      f_triggeremu->Close();
      f_triggermumu->Close();
    }
  }

  TFile file_out (fname,"RECREATE");
  TTree tree_out("analysis","main analysis") ;

//    cout<<"ev_event"<<"   "<<"sf_Ele_Reco"<<"   "<<"sf_Ele_ID"<<"      "<<"sf_Mu_ID"<<"   "<<"sf_Mu_ISO"<<"   "<<"sf_trigger"<<"   "<<"PU weight"<<endl;
  std::vector<lepton_candidate*> *selectedLeptons;
  std::vector<jet_candidate*> *selectedJets;
  TLorentzVector wp, wm, b, ab;
  bool triggerPassEE;
  bool triggerPassEMu;
  bool triggerPassMuMu;
  bool triggerPassSingleE;
  bool triggerPassSingleMu;
  bool metFilterPass;
  bool ifTopPt=false;
  int ch;
  float sf_Ele_Reco;
  float sf_Ele_ID;
  float sf_Mu_ID;
  float sf_Mu_ISO;
  float sf_Trigger;
  float weight_PU;
  float weight_Lumi;
  float weight_lep;
  float weight_lepB;
  float weight_prefiring;
  float weight_topPt;
  float elePt;
  float eleEta;
  double muPtSFRochester;
  double P_bjet_data;
  double P_bjet_mc;
  int nAccept=0;
  int nbjet;
  int nTight; //# of tight leptons
  int nLoose; //# of leptons we found in events. These leptons must pass loose requirement.
  bool isTight;//is the fakable lepton tight ?
  float Zmass;//Z candidate mass
  bool OnZ;
  bool VRveto;//When f=0, we require same-sign lepton pair. When f>0, we require three/four same-flavor leptons.
  int anti_index;//This is useful only when there are two tight leptons and one anti-selected lepton.
  

  if (fname.Contains("TTTo2L2Nu")) ifTopPt=true;

  if (fChain == 0) return;
  Long64_t nbytes = 0, nb = 0;
  Long64_t ntr = eList->GetN();
  for ( Long64_t  kentry=0; kentry<ntr;kentry++) {
    Long64_t   jentry = fChain->GetEntryNumber(kentry);//global entry number
    if (jentry < 0) break;
    Long64_t   ientry = LoadTree(jentry);//local entry number
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    displayProgress(kentry, ntr) ;

    triggerPassEE = false;
    triggerPassEMu = false;
    triggerPassMuMu = false;
    triggerPassSingleE = false;
    triggerPassSingleMu = false;
    metFilterPass = false;
    ch =10;
    muPtSFRochester = 1.;
    sf_Ele_Reco =1;
    sf_Ele_ID =1;
    sf_Mu_ID =1;
    sf_Mu_ISO =1;
    sf_Trigger =1;
    weight_PU =1;
    weight_Lumi =1;
    weight_lep =1;
    weight_lepB =1;
    weight_prefiring =1;
    weight_topPt =1;
    P_bjet_data =1;
    P_bjet_mc =1;
//MET filters

      if(data == "mc"){
        if(year == "2016" || year == "2018" ){
          if ( Flag_goodVertices==1  &&  Flag_globalSuperTightHalo2016Filter==1 && Flag_HBHENoiseFilter==1 &&  Flag_HBHENoiseIsoFilter==1 && Flag_EcalDeadCellTriggerPrimitiveFilter==1 && Flag_BadPFMuonFilter==1)
        metFilterPass = true;
        }
        if(year == "2017"){
          if ( Flag_goodVertices==1  &&  Flag_globalSuperTightHalo2016Filter==1 && Flag_HBHENoiseFilter==1 &&  Flag_HBHENoiseIsoFilter==1 && Flag_EcalDeadCellTriggerPrimitiveFilter==1 && Flag_BadPFMuonFilter==1 && Flag_ecalBadCalibFilter ==1)
        metFilterPass = true;
        }
      }

      if(data == "data"){
        if(year == "2016" || year == "2018"){
          if ( Flag_goodVertices==1  &&  Flag_globalSuperTightHalo2016Filter==1 && Flag_HBHENoiseFilter==1 &&  Flag_HBHENoiseIsoFilter==1 && Flag_EcalDeadCellTriggerPrimitiveFilter==1 && Flag_BadPFMuonFilter==1 &&  Flag_eeBadScFilter==1)
        metFilterPass = true;
        }
        if(year == "2017"){
          if ( Flag_goodVertices==1  &&  Flag_globalSuperTightHalo2016Filter==1 && Flag_HBHENoiseFilter==1 &&  Flag_HBHENoiseIsoFilter==1 && Flag_EcalDeadCellTriggerPrimitiveFilter==1 && Flag_BadPFMuonFilter==1 && Flag_ecalBadCalibFilter ==1)
        metFilterPass = true;
        }
      }

//trigger
////MC
      if(data == "mc" && year == "2016"){
        if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele35_WPTight_Gsf ) triggerPassEE =true;
        if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele35_WPTight_Gsf || HLT_IsoMu27 ) triggerPassEMu =true;
        if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL || HLT_IsoMu27 ) triggerPassMuMu =true;
        if(HLT_Ele35_WPTight_Gsf) triggerPassSingleE =true;
        if(HLT_IsoMu27) triggerPassSingleMu =true;
      }
      
      if(data == "mc" && year == "2017"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele35_WPTight_Gsf) triggerPassEE =true;
          if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele35_WPTight_Gsf || HLT_IsoMu27 ) triggerPassEMu =true;
          if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_IsoMu27) triggerPassMuMu =true;
          if(HLT_Ele35_WPTight_Gsf) triggerPassSingleE =true;
          if(HLT_IsoMu27) triggerPassSingleMu =true;
      }
      
      if(data == "mc" && year == "2018"){
        if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele35_WPTight_Gsf) triggerPassEE =true;
        if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele35_WPTight_Gsf || HLT_IsoMu27) triggerPassEMu =true;
        if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 || HLT_IsoMu27) triggerPassMuMu =true;
        if(HLT_Ele35_WPTight_Gsf) triggerPassSingleE =true;
        if(HLT_IsoMu27) triggerPassSingleMu =true;
      }
  
////DATA
      if(data == "data"){
          
          if(year == "2016"){
            if(run == "H"){
              if(dataset=="MuonEG"){
                if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEMu =true;
//                  if(!(HLT_Ele38_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&HLT_Ele27_WPTight_Gsf) {
//                      triggerPassSingleE=true;
//                  }
//                  if(!(HLT_IsoMu27||HLT_IsoMu24)&&HLT_IsoMu20) {
//                      triggerPassSingleMu =true;
//                  }
              }
              if(dataset=="SingleElectron"){
                if(!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_Ele35_WPTight_Gsf) triggerPassEE =true;
                if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_Ele35_WPTight_Gsf) triggerPassEMu =true;
                  if(HLT_Ele35_WPTight_Gsf){
                      triggerPassSingleE=true;
                  }
              }
              if(dataset=="SingleMuon"){
                if(!HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ && HLT_IsoMu27) triggerPassMuMu =true;
                if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele35_WPTight_Gsf) && HLT_IsoMu27 ) triggerPassEMu =true;
                  if(HLT_IsoMu27){
                      triggerPassSingleMu=true;
                  }
              }
              if(dataset=="DoubleEG"){
                if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEE =true;
//                  if(!HLT_Ele38_WPTight_Gsf&&HLT_Ele35_WPTight_Gsf){
//                      triggerPassSingleE=true;
//                  }
              }
              if(dataset=="DoubleMu"){
                if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ) triggerPassMuMu =true;
//                  if(!HLT_IsoMu27&&HLT_IsoMu24){
//                      triggerPassSingleMu=true;
//                  }
              }
            }
            if(run != "H"){
              if(dataset=="MuonEG"){
                if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEMu =true;
//                  if(!(HLT_Ele38_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&HLT_Ele27_WPTight_Gsf) {
//                      triggerPassSingleE=true;
//                  }
//                  if(!(HLT_IsoMu27||HLT_IsoMu24)&&HLT_IsoMu20) {
//                      triggerPassSingleMu =true;
//                  }
              }
              if(dataset=="SingleElectron"){
                if(!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_Ele35_WPTight_Gsf) triggerPassEE =true;
                if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_Ele35_WPTight_Gsf) triggerPassEMu =true;
                  if(HLT_Ele35_WPTight_Gsf){
                      triggerPassSingleE=true;
                  }
              }
              if(dataset=="SingleMuon"){
                if(!HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL && HLT_IsoMu24 ) triggerPassMuMu =true;
                if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele35_WPTight_Gsf) && HLT_IsoMu27 ) triggerPassEMu =true;
                  if(HLT_IsoMu27){
                      triggerPassSingleMu=true;
                  }
              }
              if(dataset=="DoubleEG"){
                if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEE =true;
//                  if(!HLT_Ele38_WPTight_Gsf&&HLT_Ele35_WPTight_Gsf){
//                      triggerPassSingleE=true;
//                  }
              }
              if(dataset=="DoubleMu"){
                if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ) triggerPassMuMu =true;
//                  if(!HLT_IsoMu27&&HLT_IsoMu24){
//                      triggerPassSingleMu=true;
//                  }
              }
            }
          }
          
          if(year == "2017"){
              if(dataset=="MuonEG"){
                  if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ) triggerPassEMu =true;
//                  if(!(HLT_Ele38_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&HLT_Ele27_WPTight_Gsf) {
//                      triggerPassSingleE=true;
//                  }
//                  if(!(HLT_IsoMu27||HLT_IsoMu24)&&HLT_IsoMu20) {
//                      triggerPassSingleMu =true;
//                  }
              }
              if(dataset=="SingleElectron"){
                  if(!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL && HLT_Ele35_WPTight_Gsf)) triggerPassEE =true;
                  if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_Ele35_WPTight_Gsf) triggerPassEMu =true;
                  if(HLT_Ele35_WPTight_Gsf){
                      triggerPassSingleE=true;
                  }
              }
              if(dataset=="SingleMuon"){
                  if(!HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 && HLT_IsoMu27) triggerPassMuMu =true;
                  if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ|| HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele35_WPTight_Gsf) && HLT_IsoMu27) triggerPassEMu =true;
                  if(HLT_IsoMu27){
                      triggerPassSingleMu=true;
                  }
              }
              if(dataset=="DoubleEG"){
                  if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL) triggerPassEE =true;
//                  if(!HLT_Ele38_WPTight_Gsf&&HLT_Ele35_WPTight_Gsf){
//                      triggerPassSingleE=true;
//                  }
              }
              if(dataset=="DoubleMu"){
                  if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)triggerPassMuMu =true;
//                  if(!HLT_IsoMu27&&HLT_IsoMu24){
//                      triggerPassSingleMu=true;
//                  }
              }
          }
          
          if(year == "2018"){
            if(dataset=="MuonEG"){
              if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) triggerPassEMu =true;
//                if(!(HLT_Ele38_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&HLT_Ele27_WPTight_Gsf) {
//                    triggerPassSingleE=true;
//                }
//                if(!(HLT_IsoMu27||HLT_IsoMu24)&&HLT_IsoMu20) {
//                    triggerPassSingleMu =true;
//                }
            }
            if(dataset=="EGamma"){
              if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele35_WPTight_Gsf) triggerPassEE =true;
              if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_Ele35_WPTight_Gsf) triggerPassEMu =true;
                if(HLT_Ele35_WPTight_Gsf){
                    triggerPassSingleE=true;
                }
            }
            if(dataset=="SingleMuon"){
              if(!HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 && HLT_IsoMu24) triggerPassMuMu =true;
              if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele35_WPTight_Gsf) && HLT_IsoMu27) triggerPassEMu =true;
                if(HLT_IsoMu27){
                    triggerPassSingleMu=true;
                }
            }
            if(dataset=="DoubleMu"){
              if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8) triggerPassMuMu =true;
//            if(!HLT_IsoMu27&&HLT_IsoMu24){
//                triggerPassSingleMu=true;
//            }
            }
          }
          
      }
 
    if(!(triggerPassEE || triggerPassMuMu)) continue;
    if(triggerPassEMu||triggerPassSingleE||triggerPassSingleMu) triggerPassEMu=false;
    if(!metFilterPass) continue;

// lepton selection
  selectedLeptons = new std::vector<lepton_candidate*>();
  isTight=false;
  nTight=0;
  nLoose=0;
  Zmass=0;
  OnZ=false;
  VRveto=false;
  anti_index=0;
// Muon
  int genMuIdx =0;
  for (UInt_t l=0;l< nMuon ;l++){
      if(data == "data"){
          muPtSFRochester = rc.kScaleDT(Muon_charge[l], Muon_pt[l],Muon_eta[l],Muon_phi[l], 0, 0);
      }
      if (data == "mc"){
          genMuIdx = Muon_genPartIdx[l];
          if (genMuIdx!=-1 && abs(GenPart_pdgId[genMuIdx]) == 13) muPtSFRochester = rc.kSpreadMC(Muon_charge[l], Muon_pt[l],Muon_eta[l],Muon_phi[l], GenPart_pt[genMuIdx],0, 0);
          if (genMuIdx<0 &&  Muon_nTrackerLayers[l] < 30 ) {
           muPtSFRochester = rc.kSmearMC(Muon_charge[l], Muon_pt[l] , Muon_eta[l] , Muon_phi[l], Muon_nTrackerLayers[l] , gRandom->Rndm(),0, 0);
          }
      }
      if( (muPtSFRochester * Muon_pt[l] <20) || (abs(Muon_eta[l]) > 2.4) ) continue;
      if(Muon_pfRelIso04_all[l] > 0.15) continue;//Loose Muon ID, this is used to enhance the presence of fake muons
      nLoose++;
      if (nTight==1){//If there is already a tight lepton, then we probe the second lepton found in event.
          selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          if (Muon_mediumId[l] && Muon_mvaTOP[l]>0.65) isTight = true;//Depending on the fraction of leptons passing this cut, we calculate fake rate (FR)
          continue;
      }
      if (Muon_mediumId[l] && Muon_mvaTOP[l]>0.65){//Looking for the first tight lepton
          selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
          (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
          nTight++;
      }
  }
      cout<<"////////////////////////////////"<<endl;
      cout<<"nLoose = "<<nLoose<<"; nTight = "<<nTight<<endl;
  if (!(nLoose==2&&nTight==1)){
      for (int l=0;l<(int)selectedLeptons->size();l++){
          delete (*selectedLeptons)[l];
      }
      selectedLeptons->clear();
      selectedLeptons->shrink_to_fit();
      delete selectedLeptons;
      selectedLeptons = new std::vector<lepton_candidate*>();
      nTight=0;
      nLoose=0;
      isTight=false;
// electron
    for (UInt_t l=0;l< nElectron ;l++){
      elePt = Electron_pt[l]  ;
      eleEta = Electron_eta[l] + Electron_deltaEtaSC[l];
      cout<<"Looking at "<<l<<"/"<<nElectron<<" th Electron"<<endl;
      if (elePt <20 || abs(Electron_eta[l]) > 2.4 || (abs(eleEta)> 1.4442 && (abs(eleEta)< 1.566))) continue;
      nLoose++;
      cout<<l<<" th Electron passes kinematic cuts; nLoose = "<<nLoose<<endl;
      if (nTight==1){
          cout<<"Probing "<<l<<" th Electron; nTight = "<<nTight<<endl;
      selectedLeptons->push_back(new lepton_candidate(elePt,Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
      if ((int)Electron_cutBased[l] >= 1 && Electron_mvaTOP[l] > 0.9) isTight = true;//Tight electron ID
      continue;
      }
      if ((int)Electron_cutBased[l] >= 1 && Electron_mvaTOP[l] > 0.9){
         selectedLeptons->push_back(new lepton_candidate(elePt,Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
         (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
         nTight++;
         cout<<l<<" th Electron passes Tight ID; nTight = "<<nTight<<endl;
      }
    }
  }
      if (selectedLeptons->size()>=2){
         Zmass=((*selectedLeptons)[0]->p4_+(*selectedLeptons)[1]->p4_).M();
      }
      if (Zmass>76&&Zmass<106) OnZ=true;
      
    sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
    if (selectedLeptons->size()>=2){
        if (abs((*selectedLeptons)[0]->charge_+(*selectedLeptons)[1]->charge_)>0||nLoose!=2||nTight!=1) VRveto=true;
    }
// dilepton/trilepton selection
//cout<<ev_event<<"  "<<triggerPass<<"  "<<metFilterPass<<"  "<<(int)selectedLeptons->size()<<endl;
    if((int)selectedLeptons->size()!=2 ||
      ((*selectedLeptons)[0]->pt_ <30) ||
       VRveto) {
      for (int l=0;l<(int)selectedLeptons->size();l++){
        delete (*selectedLeptons)[l];  
      }
      selectedLeptons->clear();
      selectedLeptons->shrink_to_fit();
      delete selectedLeptons;
      continue;
    }
      
    if ((*selectedLeptons)[1]->lep_ == 1) ch = 0;
    if ((*selectedLeptons)[1]->lep_ ==10) ch = 1;
    
    anti_index=(*selectedLeptons)[1]->isTag?0:1;
      
    if(ch == 0 && !triggerPassEE) continue;
    if(ch == 1 && !triggerPassMuMu) continue;
//jets
    selectedJets = new std::vector<jet_candidate*>();
    bool jetlepfail;
    for (UInt_t l=0;l< nJet;l++){
      if(data == "mc" && ((Jet_pt)[l] <30 || abs((Jet_eta)[l]) > 2.4)) continue;
      if(data == "data" && ((Jet_pt)[l] <30 || abs((Jet_eta)[l]) > 2.4)) continue;
//      if(year == "2016" && !(*jet_isJetIDTightLepVeto_2016)[l]) continue;
//      if(year == "2017" && !(*jet_isJetIDLepVeto_2017)[l]) continue;
//      if(year == "2018" && !(*jet_isJetIDLepVeto_2018)[l]) continue;
      if( (int) (Jet_jetId)[l] < 6 ) continue;
      jetlepfail = false;
      for (int i=0;i<(int)selectedLeptons->size();i++){
        if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,(Jet_eta)[l],(Jet_phi)[l]) < 0.4 ) jetlepfail=true;
      }
      if(jetlepfail) continue;
        
        float JetEnergy;
        TLorentzVector* jet_temp = new TLorentzVector() ;
        jet_temp->SetPtEtaPhiM( (Jet_pt)[l],(Jet_eta)[l],(Jet_phi)[l], (Jet_mass)[l] );
        JetEnergy = jet_temp->Energy() ;
        
      if(data == "mc"){
        selectedJets->push_back(new jet_candidate((Jet_pt)[l],(Jet_eta)[l],(Jet_phi)[l],JetEnergy,(Jet_btagDeepB)[l], year,(Jet_partonFlavour)[l]));
      }
      if(data == "data"){
        selectedJets->push_back(new jet_candidate((Jet_pt)[l],(Jet_eta)[l],(Jet_phi)[l],JetEnergy,(Jet_btagDeepB)[l],year,0));
      }
    }

    sort(selectedJets->begin(), selectedJets->end(), ComparePtJet);
    nbjet=0;
    for (int l=0;l<(int)selectedJets->size();l++){
      if((*selectedJets)[l]->btag_) nbjet++;
      if(data == "data") continue;
      if( abs((*selectedJets)[l]->flavor_) == 5){
        if( (*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"");
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),""));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }  
      }
      if( abs((*selectedJets)[l]->flavor_) == 4){
        if( (*selectedJets)[l]->btag_) {
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"");
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),""));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }
      if( abs((*selectedJets)[l]->flavor_) != 4 && abs((*selectedJets)[l]->flavor_) != 5){
        if( (*selectedJets)[l]->btag_) {
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"");
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),""));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }
    }

    if (ch==0) {
        sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");}
    if (ch==1) {
        sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");}
    if (data == "mc" && year == "2016") weight_PU = wPU.PU_2016(Pileup_nTrueInt,"nominal");
    if (data == "mc" && year == "2017") weight_PU = wPU.PU_2017(Pileup_nTrueInt,"nominal");
    if (data == "mc" && year == "2018") weight_PU = wPU.PU_2018(Pileup_nTrueInt,"nominal");
    if (data == "mc") weight_Lumi = (1000*xs*lumi)/Nevent;
    if (data == "mc" && (year == "2016" || year == "2017")){
        float weightPre;
        weightPre = 1;
        if ( L1PreFiringWeight_Nom >=0 && L1PreFiringWeight_Nom <= 1.0 ){
            weightPre =  L1PreFiringWeight_Nom;
        }
        weight_prefiring = weightPre;
    }
      if (data == "mc" && ifTopPt) {
        
        for (int l=0;l< nGenPart ;l++){
          //float GenPart_energy;
          //GenPart_energy = GenPart_pt[l];
          float GenEnergy;
          TLorentzVector* gen_temp = new TLorentzVector() ;
          gen_temp->SetPtEtaPhiM( GenPart_pt[l],(GenPart_eta)[l],(GenPart_phi)[l], (GenPart_mass)[l] );
          GenEnergy = gen_temp->Energy() ;

          if(GenPart_status[l]<30 && GenPart_status[l]>20 && GenPart_pdgId[l]==24) wp.SetPtEtaPhiE(GenPart_pt[l], GenPart_eta[l], (GenPart_phi)[l], GenEnergy) ;
          if(GenPart_status[l]<30 && GenPart_status[l]>20 && GenPart_pdgId[l]==-24) wm.SetPtEtaPhiE(GenPart_pt[l], (GenPart_eta)[l], (GenPart_phi)[l], GenEnergy) ;
          if(GenPart_status[l]<30 && GenPart_status[l]>20 && GenPart_pdgId[l]==5) b.SetPtEtaPhiE(GenPart_pt[l], (GenPart_eta)[l], (GenPart_phi)[l], GenEnergy) ;
          if(GenPart_status[l]<30 && GenPart_status[l]>20 && GenPart_pdgId[l]==-5) ab.SetPtEtaPhiE(GenPart_pt[l], (GenPart_eta)[l], (GenPart_phi)[l], GenEnergy ) ;
        }
        weight_topPt = sqrt(topPt((wp + b).Pt()) * topPt((wm + ab).Pt()));
      }
      
    float mc_w_sign;
    mc_w_sign = 1.;

    if (data == "mc") weight_lep = sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_Lumi  * mc_w_sign * weight_prefiring * weight_topPt;
    if (data == "mc") weight_lepB = sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_Lumi * mc_w_sign *  weight_prefiring * weight_topPt * (P_bjet_data/P_bjet_mc);

 if(OnZ){//Z
     cout <<"Lep Eta = "<<abs((*selectedLeptons)[anti_index]->eta_)<<endl;
     cout <<"Lep Pt = "<<(*selectedLeptons)[anti_index]->pt_<<endl;
    HistsMR[ch][0]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lep);
    if(isTight) HistsMR[ch+2][0]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lep);
      
    if ((MET_pt)>20){
      HistsMR[ch][1]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lep);
      if(isTight) HistsMR[ch+2][1]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lep);
    }
      
    if ((MET_pt)>20&&selectedJets->size()>=1&&nbjet<=1){
      HistsMR[ch][2]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lep);
      if(isTight) HistsMR[ch+2][2]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lep);
    }
      
    if (selectedJets->size()>=2&&nbjet<=1){
      HistsMR[ch][3]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
      if(isTight) HistsMR[ch+2][3]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    }
     
    if ((MET_pt)>20 && nbjet >=1){
      HistsMR[ch][4]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
      if(isTight) HistsMR[ch+2][4]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    }
     
    if (selectedJets->size()>=1 && nbjet ==0){
      HistsMR[ch][5]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
      if(isTight) HistsMR[ch+2][5]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    }
     
    if ((MET_pt)<20&&selectedJets->size()>=1 && nbjet ==0){
      HistsMR[ch][6]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
      if(isTight) HistsMR[ch+2][6]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    }
     
    if ((MET_pt)>20&&nbjet ==2){
      HistsMR[ch][7]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
      if(isTight) HistsMR[ch+2][7]->Fill(abs((*selectedLeptons)[anti_index]->eta_),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    }

     nAccept++;
  }

    for (int l=0;l<(int)selectedLeptons->size();l++){
      delete (*selectedLeptons)[l];
    }
    for (int l=0;l<(int)selectedJets->size();l++){
      delete (*selectedJets)[l];
    }
    selectedLeptons->clear();
    selectedLeptons->shrink_to_fit();
    delete selectedLeptons;
    selectedJets->clear();
    selectedJets->shrink_to_fit();
    delete selectedJets;
  } //end of event loop

  cout<<"from "<<ntr<<" evnets, "<<nAccept<<" events are accepted"<<endl;

  for (int i=0;i<(int)channels.size();++i){
    for (int k=0;k<(int)regions.size();++k){
            HistsMR[i][k] ->Write("",TObject::kOverwrite);
    }
  }

    
  for (int i=0;i<(int)channels.size();++i){
    for (int k=0;k<(int)regions.size();++k){
              delete HistsMR[i][k];
    }
  }

  file_out.Close() ;
}
