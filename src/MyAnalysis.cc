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

  typedef vector<TH1F*> Dim1;
  typedef vector<Dim1> Dim2;
  typedef vector<Dim2> Dim3;
  typedef vector<Dim3> Dim4;

  std::vector<TString> regions{"ll","llMetg20","llMetg20Jetgeq1Bleq1","llJetgeq2Bleq1","llMetg20Bgeq1","llJetgeq1B0","llMetl20Jetgeq1B0","llMetg20B2"};
  std::vector<TString> channels{"MR_e","MR_mu"};
  std::vector<TString> etaregs{"barrel","transition","endcap"};
  std::vector<TString> vars   {"FlepPt","FlepEta","FlepPhi","TlepPt","TlepEta","TlepPhi"};
  std::vector<int>    nbins   {30      ,20       ,25       ,30      ,20       ,25       };
  std::vector<float> lowEdge  {0       ,-3       ,-4       ,0       ,-3       ,-4       };
  std::vector<float> highEdge {300     ,3        ,4        ,300     ,3        ,4        };
    
  Double_t ptBinsFF[8] = {20., 23., 27., 32., 38., 45., 65., 100.};

  Dim4 HistsMR(channels.size(),Dim3(regions.size(),Dim2(etaregs.size(),Dim1(vars.size()))));
  std::stringstream name;
  TH1F *h_test;
  for (int i=0;i<(int)channels.size();++i){
    for (int k=0;k<(int)regions.size();++k){
       for (int j=0;j<(int)etaregs.size();++j){
          for (int l=0;l<(int)vars.size();++l){
            name<<channels[i]<<"_"<<regions[k]<<"_"<<etaregs[j]<<"_"<<vars[l];
            if (l==0||l==3){
            h_test = new TH1F((name.str()).c_str(),(name.str()).c_str(),7,ptBinsFF);
            }
            else{
            h_test = new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]);
            }
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsMR[i][k][j][l] = h_test;
            name.str("");
          }
      }
    }
  }

 std::vector<TString> FF{"VR","AR1","AR2"};
 std::vector<TString> regionsFF{"lll","lllMetl20","lllMetg20","lllOnZ","lllOffZ","lllOffZMetg20Jetgeq1Bleq1"};
 std::vector<TString> channelsFF{"eee","mumumu"};
 std::vector<TString> varsFF   {"lep1Pt","lep1Eta","jet1Pt","jet1Eta","njet","nbjet","Met","nVtx"};
 std::vector<int>    nbinsFF   {12      ,15       ,15      ,15       ,10    ,6      ,15   ,70};
 std::vector<float> lowEdgeFF  {30      ,-3       ,30      ,-3       ,0     ,0      ,0    ,0};
 std::vector<float> highEdgeFF {300     ,3        ,300     ,3        ,10    ,6      ,210  ,70};
    
 Dim4 HistsFF(FF.size(),Dim3((int)channelsFF.size(),Dim2(regionsFF.size(),Dim1(varsFF.size()))));
 for (int j=0;j<(int)FF.size();++j){
    for (int i=0;i<(int)channelsFF.size();++i){
        for (int k=0;k<(int)regionsFF.size();++k){
            for (int l=0;l<(int)varsFF.size();++l){
                name<<FF[j]<<"_"<<channelsFF[i]<<"_"<<regionsFF[k]<<"_"<<varsFF[l];
                h_test = new TH1F((name.str()).c_str(),(name.str()).c_str(),nbinsFF[l],lowEdgeFF[l],highEdgeFF[l]);
                h_test->StatOverflows(kTRUE);
                h_test->Sumw2(kTRUE);
                HistsFF[j][i][k][l] = h_test;
                name.str("");
            }
        }
    }
 }
    
 std::vector<TString> channelsZZ{"eeee","mumumumu"};
 Dim3 HistsZZ((int)channelsZZ.size(),Dim2(regionsFF.size(),Dim1(varsFF.size())));
 for (int i=0;i<(int)channelsZZ.size();++i){
    for (int j=0;j<(int)regionsFF.size();++j){
        for (int k=0;k<(int)varsFF.size();++k){
            name<<channelsZZ[i]<<"_"<<regionsFF[j]<<"_"<<varsFF[k];
            h_test = new TH1F((name.str()).c_str(),(name.str()).c_str(),nbinsFF[k],lowEdgeFF[k],highEdgeFF[k]);
            h_test->StatOverflows(kTRUE);
            h_test->Sumw2(kTRUE);
            HistsZZ[i][j][k] = h_test;
            name.str("");
        }
    }
 }
    
    
 typedef vector<TH2F*> Dim21;
 typedef vector<Dim21> Dim22;
 typedef vector<Dim22> Dim23;
 typedef vector<Dim23> Dim24;
 typedef vector<TH3F*> Dim31;
 typedef vector<Dim31> Dim32;
 typedef vector<Dim32> Dim33;
 typedef vector<Dim33> Dim34;
 typedef vector<Dim34> Dim35;
    
 Dim24 HistsAR1((int)channelsFF.size(),Dim23(regionsFF.size(),Dim22(etaregs.size(),Dim21(varsFF.size()))));
 Dim35 HistsAR2((int)channelsFF.size(),Dim34(regionsFF.size(),Dim33(etaregs.size(),Dim32(etaregs.size(),Dim31(varsFF.size())))));
    
 TH2F *h2_test;
 TH3F *h3_test;
 for (int i=0;i<(int)channelsFF.size();++i){
     for (int k=0;k<(int)regionsFF.size();++k){
         for (int j=0;j<(int)etaregs.size();++j){
             for (int l=0;l<(int)varsFF.size();++l){
                 //Fill Bin Array
                 Double_t BinArray[nbinsFF[l]+1];
                 for (int n=0;n<=nbinsFF[l];n++){
                     BinArray[n]=lowEdgeFF[l]+n*(highEdgeFF[l]-lowEdgeFF[l])/nbinsFF[l];
                 }
                 //Fill 2D/3D Histograms
                 name<<"AR1_2D_"<<channelsFF[i]<<"_"<<regionsFF[k]<<"_"<<etaregs[j]<<"_"<<varsFF[l];
                 h2_test = new TH2F((name.str()).c_str(),(name.str()).c_str(),nbinsFF[l],BinArray,7,ptBinsFF);
                 h2_test->StatOverflows(kTRUE);
                 h2_test->Sumw2(kTRUE);
                 HistsAR1[i][k][j][l] = h2_test;
                 name.str("");
                 for (int m=0;m<(int)etaregs.size();++m){
                     name<<"AR2_3D_"<<channelsFF[i]<<"_"<<regionsFF[k]<<"_"<<etaregs[j]<<"_"<<etaregs[m]<<"_"<<varsFF[l];
                     h3_test = new TH3F((name.str()).c_str(),(name.str()).c_str(),nbinsFF[l],BinArray,7,ptBinsFF,7,ptBinsFF);
                     h3_test->StatOverflows(kTRUE);
                     h3_test->Sumw2(kTRUE);
                     HistsAR2[i][k][j][m][l] = h3_test;
                     name.str("");
                 }
             }
         }
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
  float mZ=91.2;
  int nLeptonCut;//# of leptons we require
  bool VRveto;//When f=0, we require same-sign lepton pair. When f>0, we require three/four same-flavor leptons.
  int anti_index;//This is useful only when there are two tight leptons and one anti-selected lepton.
  int etabin1;//eta region of the first fakeable lepton
  int etabin2;//eta region of the second fakeable lepton
  

  if (fname.Contains("TTTo2L2Nu")) ifTopPt=true;

  if (fChain == 0) return;
  Long64_t nbytes = 0, nb = 0;
  Long64_t ntr = eList->GetN();
for (int f=0;f<2;f++){//f=0:MR;f=1:VR+AR1+AR2;f=2:ZZ->4l CR
  for ( Long64_t  kentry=0; kentry<ntr;kentry++) {
    Long64_t   jentry = fChain->GetEntryNumber(kentry);//global entry number
    if (jentry < 0) break;
    Long64_t   ientry = LoadTree(jentry);//local entry number
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    displayProgress(kentry+f*ntr, 2*ntr) ;

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
//      if(data == "mc" && year == "2016"){
//          if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || HLT_Ele27_WPTight_Gsf ) triggerPassEE =true;
//          if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || HLT_Ele27_WPTight_Gsf || HLT_IsoMu24 || trig_HLT_IsoTkMu24_accept) triggerPassEMu =true;
//          if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept || HLT_IsoMu24 || trig_HLT_IsoTkMu24_accept) triggerPassMuMu =true;
//          if(HLT_Ele27_WPTight_Gsf) triggerPassSingleE =true;
//          if(HLT_IsoMu24) triggerPassSingleMu =true;
//      }
      
      if(data == "mc" && year == "2017"){
          if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele35_WPTight_Gsf) triggerPassEE =true;
          if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele35_WPTight_Gsf || HLT_IsoMu27 ) triggerPassEMu =true;
          if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_IsoMu27) triggerPassMuMu =true;
          if(HLT_Ele35_WPTight_Gsf) triggerPassSingleE =true;
          if(HLT_IsoMu24) triggerPassSingleMu =true;
      }
      
//      if(data == "mc" && year == "2018"){
//          if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || HLT_Ele38_WPTight_Gsf) triggerPassEE =true;
//          if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || HLT_Ele38_WPTight_Gsf || HLT_IsoMu24) triggerPassEMu =true;
//          if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept || HLT_IsoMu24) triggerPassMuMu =true;
//          if(HLT_Ele27_WPTight_Gsf) triggerPassSingleE =true;
//          if(HLT_IsoMu24) triggerPassSingleMu =true;
//      }
  
////DATA
      if(data == "data"){
//          if(year == "2016"){
//              if(run == "H"){
//                  if(dataset=="MuonEG"){
//                      if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEMu =true;
//                      if(!(HLT_Ele38_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&HLT_Ele27_WPTight_Gsf) {
//                          triggerPassSingleE=true;
//                      }
//                      if(!(HLT_IsoMu27||HLT_IsoMu24)&&HLT_IsoMu20) {
//                          triggerPassSingleMu =true;
//                      }
//                  }
//                  if(dataset=="SingleElectron"){
//                      if(!(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && HLT_Ele27_WPTight_Gsf) triggerPassEE =true;
//                      if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && HLT_Ele27_WPTight_Gsf) triggerPassEMu =true;
//                      if(HLT_Ele38_WPTight_Gsf){
//                          triggerPassSingleE=true;
//                      }
//                  }
//                  if(dataset=="SingleMuon"){
//                      if(!(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) && (HLT_IsoMu24 || trig_HLT_IsoTkMu24_accept)) triggerPassMuMu =true;
//                      if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || HLT_Ele27_WPTight_Gsf) && (HLT_IsoMu24 || trig_HLT_IsoTkMu24_accept)) triggerPassEMu =true;
//                      if(HLT_IsoMu27){
//                          triggerPassSingleMu=true;
//                      }
//                  }
//                  if(dataset=="DoubleEG"){
//                      if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEE =true;
//                      if(!HLT_Ele38_WPTight_Gsf&&HLT_Ele27_WPTight_Gsf){
//                          triggerPassSingleE=true;
//                      }
//
//                  }
//                  if(dataset=="DoubleMu"){
//                      if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) triggerPassMuMu =true;
//                      if(!HLT_IsoMu27&&HLT_IsoMu24){
//                          triggerPassSingleMu=true;
//                      }
//
//                  }
//              }
//              if(run != "H"){
//                  if(dataset=="MuonEG"){
//                      if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) triggerPassEMu =true;
//                      if(!(HLT_Ele38_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&HLT_Ele27_WPTight_Gsf) {
//                          triggerPassSingleE=true;
//                      }
//                      if(!(HLT_IsoMu27||HLT_IsoMu24)&&HLT_IsoMu20) {
//                          triggerPassSingleMu =true;
//                      }
//                  }
//                  if(dataset=="SingleElectron"){
//                      if(!(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && HLT_Ele27_WPTight_Gsf) triggerPassEE =true;
//                      if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && HLT_Ele27_WPTight_Gsf) triggerPassEMu =true;
//                      if(HLT_Ele38_WPTight_Gsf){
//                          triggerPassSingleE=true;
//                      }
//                  }
//                  if(dataset=="SingleMuon"){
//                      if(!(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept) && (HLT_IsoMu24 || trig_HLT_IsoTkMu24_accept)) triggerPassMuMu =true;
//                      if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || HLT_Ele27_WPTight_Gsf) && (HLT_IsoMu24 || trig_HLT_IsoTkMu24_accept)) triggerPassEMu =true;
//                      if(HLT_IsoMu27){
//                          triggerPassSingleMu=true;
//                      }
//                  }
//                  if(dataset=="DoubleEG"){
//                      if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEE =true;
//                      if(!HLT_Ele38_WPTight_Gsf&&HLT_Ele27_WPTight_Gsf){
//                          triggerPassSingleE=true;
//                      }
//                  }
//                  if(dataset=="DoubleMu"){
//                      if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) triggerPassMuMu =true;
//                      if(!HLT_IsoMu27&&HLT_IsoMu24){
//                          triggerPassSingleMu=true;
//                      }
//                  }
//              }
//          }
          if(year == "2017"){
              if(dataset=="MuonEG"){
                  if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ) triggerPassEMu =true;
                  if(!(HLT_Ele38_WPTight_Gsf||HLT_Ele35_WPTight_Gsf)&&HLT_Ele27_WPTight_Gsf) {
                      triggerPassSingleE=true;
                  }
                  if(!(HLT_IsoMu27||HLT_IsoMu24)&&HLT_IsoMu20) {
                      triggerPassSingleMu =true;
                  }
              }
              if(dataset=="SingleElectron"){
                  if(!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL && HLT_Ele35_WPTight_Gsf)) triggerPassEE =true;
                  if(!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_Ele35_WPTight_Gsf) triggerPassEMu =true;
                  if(HLT_Ele38_WPTight_Gsf){
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
                  if(!HLT_Ele38_WPTight_Gsf&&HLT_Ele35_WPTight_Gsf){
                      triggerPassSingleE=true;
                  }
              }
              if(dataset=="DoubleMu"){
                  if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)triggerPassMuMu =true;
                  if(!HLT_IsoMu27&&HLT_IsoMu24){
                      triggerPassSingleMu=true;
                  }
              }
          }
//          if(year == "2018"){
//              if(dataset=="MuonEG"){
//                  if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) triggerPassEMu =true;
//                  if(!(HLT_Ele38_WPTight_Gsf||HLT_Ele27_WPTight_Gsf)&&HLT_Ele27_WPTight_Gsf) {
//                      triggerPassSingleE=true;
//                  }
//                  if(!(HLT_IsoMu27||HLT_IsoMu24)&&HLT_IsoMu20) {
//                      triggerPassSingleMu =true;
//                  }
//              }
//              if(dataset=="EGamma"){
//                  if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || HLT_Ele38_WPTight_Gsf) triggerPassEE =true;
//                  if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && HLT_Ele38_WPTight_Gsf) triggerPassEMu =true;
//                  if(HLT_Ele38_WPTight_Gsf){
//                      triggerPassSingleE=true;
//                  }
//              }
//              if(dataset=="SingleMuon"){
//                  if(!trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept && HLT_IsoMu24) triggerPassMuMu =true;
//                  if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || HLT_Ele38_WPTight_Gsf) && HLT_IsoMu24) triggerPassEMu =true;
//                  if(HLT_IsoMu27){
//                      triggerPassSingleMu=true;
//                  }
//              }
//              if(dataset=="DoubleMu"){
//                  if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept) triggerPassMuMu =true;
//                  if(!HLT_IsoMu27&&HLT_IsoMu24){
//                      triggerPassSingleMu=true;
//                  }
//              }
//          }
      }
 


//cout<<ev_event<<"  "<<trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept <<"  "<< trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept <<"  "<< HLT_Ele27_WPTight_Gsf  <<"  "<<trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept <<"  "<< trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept <<"  "<< trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept<<"  "<<HLT_IsoMu24<<"  "<<trig_HLT_IsoTkMu24_accept<<endl;
//cout<<"trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept "<< "trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept "<< "HLT_Ele27_WPTight_Gsf " <<"trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept "<< "trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept "<< "trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept "  <<"HLT_IsoMu24 "<<"trig_HLT_IsoTkMu24_accept"<<endl;
    if(!(triggerPassSingleE || triggerPassSingleMu)&&f==0) continue;
    //if(!(triggerPassEE || triggerPassEMu || triggerPassMuMu)&&f>0) continue;
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
  etabin1=0;
  etabin2=0;
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
      if(!Muon_mediumId[l]||Muon_pfRelIso04_all[l] > 0.15) continue;//Loose Muon ID, this is used to enhance the presence of fake muons
      nLoose++;
      if (f==0){//f=0 -> MR
          if (nTight==1){//If there is already a tight lepton, then we probe the second lepton found in event.
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_pt[l], abs(Muon_eta[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_pt[l], abs(Muon_eta[l]),"");
              if (Muon_mvaTOP[l]>0.65) isTight = true;//Depending on the fraction of leptons passing this cut, we calculate fake rate (FR)
              continue;
          }
          if (Muon_mvaTOP[l]>0.65){//Looking for the first tight lepton
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
              (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_pt[l], abs(Muon_eta[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_pt[l], abs(Muon_eta[l]),"");
              nTight++;
          }
      }
      else{//f>0 -> VR+AR
          if (Muon_mvaTOP[l]>0.65){
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
              (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_pt[l], abs(Muon_eta[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_pt[l], abs(Muon_eta[l]),"");
              nTight++;
          }
          else{
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_pt[l], abs(Muon_eta[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_pt[l], abs(Muon_eta[l]),"");
          }
      }
  }
// electron
    for (UInt_t l=0;l< nElectron ;l++){
      elePt = Electron_pt[l]  ;
      eleEta = Electron_eta[l] + Electron_deltaEtaSC[l];
      if (elePt <20 || abs(Electron_eta[l]) > 2.4 || (abs(eleEta)> 1.4442 && (abs(eleEta)< 1.566))) continue;
      if ((int)Electron_cutBased[l] < 1) continue;//Loose electron ID
      nLoose++;
      if (f==0){
         if (nTight==1){
         selectedLeptons->push_back(new lepton_candidate(elePt,Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
         if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,eleEta,elePt,"");
         if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,eleEta,elePt,"");
         if (Electron_mvaTOP[l] > 0.9) isTight = true;//Tight electron ID
         continue;
         }
         if (Electron_mvaTOP[l] > 0.9){
            selectedLeptons->push_back(new lepton_candidate(elePt,Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
            (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
            if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,eleEta,elePt,"");
            if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,eleEta,elePt,"");
            nTight++;
         }
      }
      else{
      if (Electron_mvaTOP[l] > 0.9){
         selectedLeptons->push_back(new lepton_candidate(elePt,Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
         (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
         if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,eleEta,elePt,"");
         if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,eleEta,elePt,"");
         nTight++;
      }
      else{
      selectedLeptons->push_back(new lepton_candidate(elePt,Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
      if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,eleEta,elePt,"");
      if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,eleEta,elePt,"");
      }
      }
    }
    if (f==0&&nLoose==2){//This solves the problem concerning the order of lepton loops
      sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
      if (selectedLeptons->size()!=2 ||
         !((*selectedLeptons)[0]->isTag)){
          for (int l=0;l<(int)selectedLeptons->size();l++){
              delete (*selectedLeptons)[l];
          }
          selectedLeptons->clear();
          selectedLeptons->shrink_to_fit();
          delete selectedLeptons;
          selectedLeptons = new std::vector<lepton_candidate*>();
          nTight=0;
          isTight=false;
          sf_Ele_Reco=1;
          sf_Ele_ID=1;
          sf_Mu_ID=1;
          sf_Mu_ISO=1;
      // electron
      for (int l=nElectron-1;l>-1;l--){
          elePt = Electron_pt[l]  ;
          eleEta = Electron_eta[l] + Electron_deltaEtaSC[l];
          if (elePt <20 || abs(Electron_eta[l]) > 2.4 || (abs(eleEta)> 1.4442 && (abs(eleEta)< 1.566))) continue;
          if ((int)Electron_cutBased[l] < 1) continue;
          if (nTight==1){
             selectedLeptons->push_back(new lepton_candidate(elePt,Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
             if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,eleEta,elePt,"");
             if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,eleEta,elePt,"");
             if (Electron_mvaTOP[l] > 0.9) isTight = true;
             continue;
          }
          if (Electron_mvaTOP[l] > 0.9){
             selectedLeptons->push_back(new lepton_candidate(elePt,Electron_eta[l],Electron_phi[l],Electron_charge[l],l,1));
             (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
             if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,eleEta,elePt,"");
             if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,eleEta,elePt,"");
             nTight++;
          }
      }
      // Muon
      for (int l=nMuon-1;l>-1;l--){
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
          if(muPtSFRochester * Muon_pt[l] <20 || abs(Muon_eta[l]) > 2.4) continue;
          if(!Muon_mediumId[l]||Muon_pfRelIso04_all[l] > 0.15) continue;
          if (nTight==1){
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_pt[l], abs(Muon_eta[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_pt[l], abs(Muon_eta[l]),"");
              if (Muon_mvaTOP[l]>0.65) isTight = true;
              continue;
          }
          if (Muon_mvaTOP[l]>0.65){
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
              (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_eta[l], Muon_pt[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_pt[l], abs(Muon_eta[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_pt[l], abs(Muon_eta[l]),"");
              nTight++;
          }
      }
    }
    }
      if (selectedLeptons->size()==2&&(*selectedLeptons)[0]->lep_+(*selectedLeptons)[1]->lep_==2){
         Zmass=((*selectedLeptons)[0]->p4_+(*selectedLeptons)[1]->p4_).M();
      }
      if (selectedLeptons->size()==3 && abs((*selectedLeptons)[0]->charge_+(*selectedLeptons)[1]->charge_+(*selectedLeptons)[2]->charge_)==1){
          sort(selectedLeptons->begin(), selectedLeptons->end(), CompareChargeLep);
          if ((*selectedLeptons)[0]->isTag&&(*selectedLeptons)[2]->isTag){
          Zmass=((*selectedLeptons)[0]->p4_+(*selectedLeptons)[2]->p4_).M();
          }
          if ((*selectedLeptons)[0]->charge_+(*selectedLeptons)[1]->charge_==0){
              if ((*selectedLeptons)[0]->isTag&&(*selectedLeptons)[1]->isTag&&abs(Zmass-mZ)>abs(((*selectedLeptons)[0]->p4_+(*selectedLeptons)[1]->p4_).M()-mZ)){
              Zmass=((*selectedLeptons)[0]->p4_+(*selectedLeptons)[1]->p4_).M();
              }
          }
          else{
              if ((*selectedLeptons)[1]->isTag&&(*selectedLeptons)[2]->isTag&&abs(Zmass-mZ)>abs(((*selectedLeptons)[1]->p4_+(*selectedLeptons)[2]->p4_).M()-mZ)){
              Zmass=((*selectedLeptons)[1]->p4_+(*selectedLeptons)[2]->p4_).M();
              }
          }
      }
      if (selectedLeptons->size()==4){//4-lepton ZZ control region
          sort(selectedLeptons->begin(), selectedLeptons->end(), CompareChargeLep);
          Zmass=((*selectedLeptons)[0]->p4_+(*selectedLeptons)[2]->p4_).M();
          if (abs(Zmass-mZ)>abs(((*selectedLeptons)[0]->p4_+(*selectedLeptons)[3]->p4_).M()-mZ)){
              Zmass=((*selectedLeptons)[0]->p4_+(*selectedLeptons)[3]->p4_).M();
          }
          if (abs(Zmass-mZ)>abs(((*selectedLeptons)[1]->p4_+(*selectedLeptons)[2]->p4_).M()-mZ)){
              Zmass=((*selectedLeptons)[1]->p4_+(*selectedLeptons)[2]->p4_).M();
          }
          if (abs(Zmass-mZ)>abs(((*selectedLeptons)[1]->p4_+(*selectedLeptons)[3]->p4_).M()-mZ)){
              Zmass=((*selectedLeptons)[1]->p4_+(*selectedLeptons)[3]->p4_).M();
          }
      }
      if (Zmass>76&&Zmass<106) OnZ=true;
      
    sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
    nLeptonCut=f+2;
    if (selectedLeptons->size()==2){
        if ((*selectedLeptons)[0]->charge_+(*selectedLeptons)[1]->charge_==0||nLoose!=2) VRveto=true;
    }
    if (selectedLeptons->size()==3){
        if ((*selectedLeptons)[0]->lep_+(*selectedLeptons)[1]->lep_+(*selectedLeptons)[2]->lep_>3&&
            (*selectedLeptons)[0]->lep_+(*selectedLeptons)[1]->lep_+(*selectedLeptons)[2]->lep_<30) VRveto=true;
    }
    if (selectedLeptons->size()==4){
        if ((*selectedLeptons)[0]->charge_+(*selectedLeptons)[1]->charge_+(*selectedLeptons)[2]->charge_+(*selectedLeptons)[3]->charge_==0){
            if ((*selectedLeptons)[0]->lep_+(*selectedLeptons)[1]->lep_+(*selectedLeptons)[2]->lep_+(*selectedLeptons)[3]->lep_>4&&
                (*selectedLeptons)[0]->lep_+(*selectedLeptons)[1]->lep_+(*selectedLeptons)[2]->lep_+(*selectedLeptons)[3]->lep_<40){
                VRveto=true;
            }
        }
        else{
        VRveto=true;
        }
    }
// dilepton/trilepton selection
//cout<<ev_event<<"  "<<triggerPass<<"  "<<metFilterPass<<"  "<<(int)selectedLeptons->size()<<endl;
    if((int)selectedLeptons->size()!=nLeptonCut ||
      ((*selectedLeptons)[0]->pt_ <30) ||
       VRveto || !((*selectedLeptons)[0]->isTag)) {
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
    
    anti_index=(*selectedLeptons)[1]->isTag?2:1;
    //Eta bins
    if(f==1&&nTight==1){
       if(abs((*selectedLeptons)[1]->eta_)>1.6){
          etabin1=2;
        }
       else if(abs((*selectedLeptons)[1]->eta_)>0.8){
          etabin1=1;
       }
        
       if(abs((*selectedLeptons)[2]->eta_)>1.6){
          etabin2=2;
       }
       else if(abs((*selectedLeptons)[2]->eta_)>0.8){
            etabin2=1;
       }
    }
    else{
        if(abs((*selectedLeptons)[anti_index]->eta_)>1.6){
            etabin1=2;
        }
        else if(abs((*selectedLeptons)[anti_index]->eta_)>0.8){
            etabin1=1;
        }
    }
      
    if (nTight==1){
       if ((*selectedLeptons)[0]->lep_ == 1 && !triggerPassSingleE) continue;
       if ((*selectedLeptons)[0]->lep_ == 10 && !triggerPassSingleMu) continue;
      }
    else{
    if((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 2 && !triggerPassEE) continue;
    if((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 11 && !triggerPassEMu) continue;
    if((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 20 && !triggerPassMuMu) continue;
    }
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

    if (data == "mc" && (*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 2) {
        sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");}
    if (data == "mc" && (*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 11) {
        sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");}
    if (data == "mc" && (*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 20) {
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

 if(f==0&&!OnZ){//Z veto by default
    HistsMR[ch][0][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    HistsMR[ch][0][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    HistsMR[ch][0][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    if(isTight) HistsMR[ch][0][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    if(isTight) HistsMR[ch][0][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    if(isTight) HistsMR[ch][0][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
      
    if ((MET_pt)>20){
      HistsMR[ch][1][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
      HistsMR[ch][1][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
      HistsMR[ch][1][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
      if(isTight) HistsMR[ch][1][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
      if(isTight) HistsMR[ch][1][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
      if(isTight) HistsMR[ch][1][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    }
      
    if ((MET_pt)>20&&selectedJets->size()>=1&&nbjet<=1){
      HistsMR[ch][2][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
      HistsMR[ch][2][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
      HistsMR[ch][2][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
      if(isTight) HistsMR[ch][2][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
      if(isTight) HistsMR[ch][2][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
      if(isTight) HistsMR[ch][2][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    }
      
    if (selectedJets->size()>=2&&nbjet<=1){
      HistsMR[ch][3][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsMR[ch][3][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsMR[ch][3][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      if(isTight) HistsMR[ch][3][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      if(isTight) HistsMR[ch][3][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      if(isTight) HistsMR[ch][3][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    }
     
    if ((MET_pt)>20 && nbjet >=1){
      HistsMR[ch][4][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsMR[ch][4][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsMR[ch][4][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      if(isTight) HistsMR[ch][4][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      if(isTight) HistsMR[ch][4][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      if(isTight) HistsMR[ch][4][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    }
     
    if (selectedJets->size()>=1 && nbjet ==0){
      HistsMR[ch][5][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsMR[ch][5][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsMR[ch][5][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      if(isTight) HistsMR[ch][5][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      if(isTight) HistsMR[ch][5][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      if(isTight) HistsMR[ch][5][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    }
     
    if ((MET_pt)<20&&selectedJets->size()>=1 && nbjet ==0){
      HistsMR[ch][6][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsMR[ch][6][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsMR[ch][6][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      if(isTight) HistsMR[ch][6][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      if(isTight) HistsMR[ch][6][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      if(isTight) HistsMR[ch][6][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    }
     
    if ((MET_pt)>20&&nbjet ==2){
      HistsMR[ch][7][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsMR[ch][7][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsMR[ch][7][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      if(isTight) HistsMR[ch][7][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      if(isTight) HistsMR[ch][7][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      if(isTight) HistsMR[ch][7][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    }

     nAccept++;
  }
  else if (f==1){
      
  HistsFF[3-nTight][ch][0][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][0][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][0][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][0][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][0][4]->Fill(selectedJets->size(),weight_lep);
  HistsFF[3-nTight][ch][0][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][0][6]->Fill((MET_pt),weight_lep);
  HistsFF[3-nTight][ch][0][7]->Fill(Pileup_nTrueInt,weight_lep);
    
  if ((MET_pt)<20){
  HistsFF[3-nTight][ch][1][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][1][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][1][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][1][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][1][4]->Fill(selectedJets->size(),weight_lep);
  HistsFF[3-nTight][ch][1][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][1][6]->Fill((MET_pt),weight_lep);
  HistsFF[3-nTight][ch][1][7]->Fill(Pileup_nTrueInt,weight_lep);
  }
  else{
  HistsFF[3-nTight][ch][2][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][2][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][2][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][2][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][2][4]->Fill(selectedJets->size(),weight_lep);
  HistsFF[3-nTight][ch][2][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][2][6]->Fill((MET_pt),weight_lep);
  HistsFF[3-nTight][ch][2][7]->Fill(Pileup_nTrueInt,weight_lep);
  }
      
  if (OnZ){
  HistsFF[3-nTight][ch][3][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][3][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][3][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][3][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][3][4]->Fill(selectedJets->size(),weight_lep);
  HistsFF[3-nTight][ch][3][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][3][6]->Fill((MET_pt),weight_lep);
  HistsFF[3-nTight][ch][3][7]->Fill(Pileup_nTrueInt,weight_lep);
  }
  else{
  HistsFF[3-nTight][ch][4][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][4][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][4][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][4][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][4][4]->Fill(selectedJets->size(),weight_lep);
  HistsFF[3-nTight][ch][4][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][4][6]->Fill((MET_pt),weight_lep);
  HistsFF[3-nTight][ch][4][7]->Fill(Pileup_nTrueInt,weight_lep);
  }
      
  if ((MET_pt)>20&&!OnZ&&nbjet<=1&&selectedJets->size()>=1){
  HistsFF[3-nTight][ch][5][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
  HistsFF[3-nTight][ch][5][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][5][2]->Fill((*selectedJets)[0]->pt_,weight_lepB);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][5][3]->Fill((*selectedJets)[0]->eta_,weight_lepB);
  HistsFF[3-nTight][ch][5][4]->Fill(selectedJets->size(),weight_lepB);
  HistsFF[3-nTight][ch][5][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][5][6]->Fill((MET_pt),weight_lepB);
  HistsFF[3-nTight][ch][5][7]->Fill(Pileup_nTrueInt,weight_lepB);
  }
  //AR1:ll\overline{l}
  if (nTight==2){
      
  HistsAR1[ch][0][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][0][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][0][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][0][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][0][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][0][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][0][etabin1][6]->Fill((MET_pt),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][0][etabin1][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
    
  if ((MET_pt)<20){
  HistsAR1[ch][1][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][1][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][1][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][1][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][1][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][1][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][1][etabin1][6]->Fill((MET_pt),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][1][etabin1][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  }
  else{
  HistsAR1[ch][2][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][2][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][2][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][2][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][2][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][2][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][2][etabin1][6]->Fill((MET_pt),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][2][etabin1][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  }
      
  if (OnZ){
  HistsAR1[ch][3][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][3][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][3][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][3][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][3][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][3][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][3][etabin1][6]->Fill((MET_pt),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][3][etabin1][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  }
  else{
  HistsAR1[ch][4][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][4][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][4][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][4][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][4][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][4][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][4][etabin1][6]->Fill((MET_pt),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][4][etabin1][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  }
      
  if ((MET_pt)>20&&!OnZ&&nbjet<=1&&selectedJets->size()>=1){
  HistsAR1[ch][5][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][5][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  if(selectedJets->size()>0) HistsAR1[ch][5][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  if(selectedJets->size()>0) HistsAR1[ch][5][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][5][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][5][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][5][etabin1][6]->Fill((MET_pt),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][5][etabin1][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  }
}
  //AR2:l\overline{l}\overline{l}
  if (nTight==1){
  
  HistsAR2[ch][0][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][0][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][0][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][0][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][0][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][0][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][0][etabin1][etabin2][6]->Fill((MET_pt),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][0][etabin1][etabin2][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  
  if ((MET_pt)<20){
  HistsAR2[ch][1][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][1][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][1][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][1][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][1][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][1][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][1][etabin1][etabin2][6]->Fill((MET_pt),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][1][etabin1][etabin2][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  }
  else{
  HistsAR2[ch][2][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][2][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][2][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][2][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][2][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][2][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][2][etabin1][etabin2][6]->Fill((MET_pt),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][2][etabin1][etabin2][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  }
      
  if (OnZ){
  HistsAR2[ch][3][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][3][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][3][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][3][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][3][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][3][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][3][etabin1][etabin2][6]->Fill((MET_pt),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][3][etabin1][etabin2][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  }
  else{
  HistsAR2[ch][4][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][4][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][4][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][4][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][4][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][4][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][4][etabin1][etabin2][6]->Fill((MET_pt),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][4][etabin1][etabin2][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  }
      
  if ((MET_pt)>20&&!OnZ&&nbjet<=1&&selectedJets->size()>=1){
  HistsAR2[ch][5][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][5][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  if(selectedJets->size()>0) HistsAR2[ch][5][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  if(selectedJets->size()>0) HistsAR2[ch][5][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][5][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][5][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][5][etabin1][etabin2][6]->Fill((MET_pt),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][5][etabin1][etabin2][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  }
      
  }
      
 }
  else if (nTight==4){
      
    HistsZZ[ch][0][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    HistsZZ[ch][0][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) HistsZZ[ch][0][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) HistsZZ[ch][0][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
    HistsZZ[ch][0][4]->Fill(selectedJets->size(),weight_lep);
    HistsZZ[ch][0][5]->Fill(nbjet,weight_lepB);
    HistsZZ[ch][0][6]->Fill((MET_pt),weight_lep);
    HistsZZ[ch][0][7]->Fill(Pileup_nTrueInt,weight_lep);

    if ((MET_pt)<20){
    HistsZZ[ch][1][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    HistsZZ[ch][1][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) HistsZZ[ch][1][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) HistsZZ[ch][1][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
    HistsZZ[ch][1][4]->Fill(selectedJets->size(),weight_lep);
    HistsZZ[ch][1][5]->Fill(nbjet,weight_lepB);
    HistsZZ[ch][1][6]->Fill((MET_pt),weight_lep);
    HistsZZ[ch][1][7]->Fill(Pileup_nTrueInt,weight_lep);
    }
    else{
    HistsZZ[ch][2][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    HistsZZ[ch][2][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) HistsZZ[ch][2][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) HistsZZ[ch][2][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
    HistsZZ[ch][2][4]->Fill(selectedJets->size(),weight_lep);
    HistsZZ[ch][2][5]->Fill(nbjet,weight_lepB);
    HistsZZ[ch][2][6]->Fill((MET_pt),weight_lep);
    HistsZZ[ch][2][7]->Fill(Pileup_nTrueInt,weight_lep);
    }

    if (OnZ){
    HistsZZ[ch][3][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    HistsZZ[ch][3][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) HistsZZ[ch][3][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) HistsZZ[ch][3][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
    HistsZZ[ch][3][4]->Fill(selectedJets->size(),weight_lep);
    HistsZZ[ch][3][5]->Fill(nbjet,weight_lepB);
    HistsZZ[ch][3][6]->Fill((MET_pt),weight_lep);
    HistsZZ[ch][3][7]->Fill(Pileup_nTrueInt,weight_lep);
    }
    else{
    HistsZZ[ch][4][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    HistsZZ[ch][4][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) HistsZZ[ch][4][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) HistsZZ[ch][4][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
    HistsZZ[ch][4][4]->Fill(selectedJets->size(),weight_lep);
    HistsZZ[ch][4][5]->Fill(nbjet,weight_lepB);
    HistsZZ[ch][4][6]->Fill((MET_pt),weight_lep);
    HistsZZ[ch][4][7]->Fill(Pileup_nTrueInt,weight_lep);
    }

    if ((MET_pt)>20&&!OnZ&&nbjet<=1&&selectedJets->size()>=1){
    HistsZZ[ch][5][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    HistsZZ[ch][5][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) HistsZZ[ch][5][2]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) HistsZZ[ch][5][3]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    HistsZZ[ch][5][4]->Fill(selectedJets->size(),weight_lepB);
    HistsZZ[ch][5][5]->Fill(nbjet,weight_lepB);
    HistsZZ[ch][5][6]->Fill((MET_pt),weight_lepB);
    HistsZZ[ch][5][7]->Fill(Pileup_nTrueInt,weight_lepB);
    }
      
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
    
}
  cout<<"from "<<ntr<<" evnets, "<<nAccept<<" events are accepted"<<endl;

  for (int i=0;i<(int)channels.size();++i){
    for (int k=0;k<(int)regions.size();++k){
      for (int j=0;j<(int)etaregs.size();++j){
        for (int l=0;l<(int)vars.size();++l){
            HistsMR[i][k][j][l]  ->Write("",TObject::kOverwrite);
        }
      }
    }
  }
    
  for (int j=0;j<(int)FF.size();++j){
    for (int i=0;i<(int)channelsFF.size();++i){
        for (int k=0;k<(int)regionsFF.size();++k){
            for (int l=0;l<(int)varsFF.size();++l){
                HistsFF[j][i][k][l]  ->Write("",TObject::kOverwrite);
            }
        }
    }
  }
    
  for (int i=0;i<(int)channelsZZ.size();++i){
        for (int j=0;j<(int)regionsFF.size();++j){
            for (int k=0;k<(int)varsFF.size();++k){
                HistsZZ[i][j][k]  ->Write("",TObject::kOverwrite);
        }
    }
  }

  for (int i=0;i<(int)channelsFF.size();++i){
      for (int k=0;k<(int)regionsFF.size();++k){
          for (int j=0;j<(int)etaregs.size();++j){
              for (int l=0;l<(int)varsFF.size();++l){
                  HistsAR1[i][k][j][l]  ->Write("",TObject::kOverwrite);
                  for (int m=0;m<(int)etaregs.size();++m){
                      HistsAR2[i][k][j][m][l]  ->Write("",TObject::kOverwrite);
                 }
              }
          }
      }
  }
    
  for (int i=0;i<(int)channels.size();++i){
    for (int k=0;k<(int)regions.size();++k){
      for (int j=0;j<(int)etaregs.size();++j){
          for (int l=0;l<(int)vars.size();++l){
              delete HistsMR[i][k][j][l];
          }
      }
    }
  }
    
  for (int j=0;j<(int)FF.size();++j){
      for (int i=0;i<(int)channelsFF.size();++i){
          for (int k=0;k<(int)regionsFF.size();++k){
              for (int l=0;l<(int)varsFF.size();++l){
              delete HistsFF[j][i][k][l];
              }
          }
      }
  }
    
  for (int i=0;i<(int)channelsZZ.size();++i){
        for (int j=0;j<(int)regionsFF.size();++j){
            for (int k=0;k<(int)varsFF.size();++k){
            delete HistsZZ[i][j][k];
            }
        }
    }
    
  for (int i=0;i<(int)channelsFF.size();++i){
      for (int k=0;k<(int)regionsFF.size();++k){
         for (int j=0;j<(int)etaregs.size();++j){
             for (int l=0;l<(int)varsFF.size();++l){
                 delete HistsAR1[i][k][j][l];
                 for (int m=0;m<(int)etaregs.size();++m){
                     delete HistsAR2[i][k][j][m][l];
                 }
             }
          }
      }
  }
    
  file_out.Close() ;
}
