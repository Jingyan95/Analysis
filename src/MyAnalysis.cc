#define MyAnalysis_cxx
#include "MyAnalysis.h"
#include "PU_reWeighting.h"
#include "lepton_candidate.h"
#include "jet_candidate.h"
#include "TRandom.h"
#include "TRandom3.h"
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
#include "RoccoR.h"
#include "BTagCalibrationStandalone.h"

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

  std::vector<TString> regions{"ll","llMetg20","llMetg20Jetgeq1Bleq1","llMetg20Jetgeq2Bleq1","llMetg20Bgeq1"};
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
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
       for (int j=0;j<etaregs.size();++j){
          for (int l=0;l<vars.size();++l){
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
 std::vector<TString> regionsFF{"lll","lllMetl20","lllMetg20","lllOnZ","lllOffZ","lllMetg20Jetgeq1Bleq1"};
 std::vector<TString> channelsFF{"eee","mumumu"};
 std::vector<TString> varsFF   {"lep1Pt","lep1Eta","jet1Pt","jet1Eta","njet","nbjet","Met","nVtx"};
 std::vector<int>    nbinsFF   {12      ,15       ,15      ,15       ,10    ,6      ,15   ,70};
 std::vector<float> lowEdgeFF  {30      ,-3       ,30      ,-3       ,0     ,0      ,0    ,0};
 std::vector<float> highEdgeFF {300     ,3        ,300     ,3        ,10    ,6      ,210  ,70};
    
 Dim4 HistsFF(FF.size(),Dim3(channelsFF.size(),Dim2(regionsFF.size(),Dim1(varsFF.size()))));
 for (int j=0;j<FF.size();++j){
    for (int i=0;i<channelsFF.size();++i){
        for (int k=0;k<regionsFF.size();++k){
            for (int l=0;l<varsFF.size();++l){
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
    
 typedef vector<TH2F*> Dim21;
 typedef vector<Dim21> Dim22;
 typedef vector<Dim22> Dim23;
 typedef vector<Dim23> Dim24;
 typedef vector<TH3F*> Dim31;
 typedef vector<Dim31> Dim32;
 typedef vector<Dim32> Dim33;
 typedef vector<Dim33> Dim34;
 typedef vector<Dim34> Dim35;
    
 Dim24 HistsAR1(channelsFF.size(),Dim23(regionsFF.size(),Dim22(etaregs.size(),Dim21(varsFF.size()))));
 Dim35 HistsAR2(channelsFF.size(),Dim34(regionsFF.size(),Dim33(etaregs.size(),Dim32(etaregs.size(),Dim31(varsFF.size())))));
    
 TH2F *h2_test;
 TH3F *h3_test;
 for (int i=0;i<channelsFF.size();++i){
     for (int k=0;k<regionsFF.size();++k){
         for (int j=0;j<etaregs.size();++j){
             for (int l=0;l<varsFF.size();++l){
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
                 for (int m=0;m<etaregs.size();++m){
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

  RoccoR  rc;
  if(year == "2016")    rochesterFile = "input/RoccoR2016.txt";
  if(year == "2017")    rochesterFile = "input/RoccoR2017.txt";
  if(year == "2018")    rochesterFile = "input/RoccoR2018.txt";
  rc.init(rochesterFile);

  if(data == "mc"){
    TFile *f_btagEff_Map = new TFile("input/btagEff.root");
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

    if(year == "2016")    btagFile = "input/DeepCSV_2016LegacySF_V1.csv";
    if(year == "2017")    btagFile = "input/DeepCSV_94XSF_V4_B_F.csv";
    if(year == "2018")    btagFile = "input/DeepCSV_102XSF_V1.csv";

    BTagCalibration calib("DeepCSV",btagFile);
    reader.load(calib,BTagEntry::FLAV_B,"comb"); 
    reader.load(calib,BTagEntry::FLAV_C,"comb");
    reader.load(calib,BTagEntry::FLAV_UDSG,"comb");

    if(year == "2016"){
      TFile *f_Ele_Reco_Map = new TFile("input/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("input/2016LegacyReReco_ElectronTight_Fall17V2.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map_1 = new TFile("input/2016_RunBCDEF_SF_ID.root");
      TH2F *sf_Mu_ID_H_1 = (TH2F*)f_Mu_ID_Map_1->Get("NUM_TightID_DEN_genTracks_eta_pt");
      TFile *f_Mu_ID_Map_2 = new TFile("input/2016_RunGH_SF_ID.root");
      TH2F *sf_Mu_ID_H_2 = (TH2F*)f_Mu_ID_Map_2->Get("NUM_TightID_DEN_genTracks_eta_pt");
      sf_Mu_ID_H_1->Scale(0.55);
      sf_Mu_ID_H_2->Scale(0.45);
      sf_Mu_ID_H_1->Add(sf_Mu_ID_H_2);
      sf_Mu_ID_H = *sf_Mu_ID_H_1;

      TFile *f_Mu_ISO_Map_1 = new TFile("input/2016_RunBCDEF_SF_ISO.root");
      TH2F *sf_Mu_ISO_H_1 = (TH2F*)f_Mu_ISO_Map_1->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
      TFile *f_Mu_ISO_Map_2 = new TFile("input/2016_RunGH_SF_ISO.root");
      TH2F *sf_Mu_ISO_H_2 = (TH2F*)f_Mu_ISO_Map_2->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
      sf_Mu_ISO_H_1->Scale(0.55);
      sf_Mu_ISO_H_2->Scale(0.45);
      sf_Mu_ISO_H_1->Add(sf_Mu_ISO_H_2);
      sf_Mu_ISO_H = *sf_Mu_ISO_H_1;

      TFile *f_triggeree = new TFile("input/TriggerSF_ee2016_pt.root");
      sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
      TFile *f_triggeremu = new TFile("input/TriggerSF_emu2016_pt.root");
      sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("h_lep1Pt_lep2Pt_Step3");
      TFile *f_triggermumu = new TFile("input/TriggerSF_mumu2016_pt.root");
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
      TFile *f_Ele_Reco_Map = new TFile("input/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("input/2017_ElectronTight.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map = new TFile("input/2017_RunBCDEF_SF_ID_syst.root");
      sf_Mu_ID_H = *(TH2F*)f_Mu_ID_Map->Get("NUM_TightID_DEN_genTracks_pt_abseta");

      TFile *f_Mu_ISO_Map = new TFile("input/2017_RunBCDEF_SF_ISO_syst.root");
      sf_Mu_ISO_H = *(TH2F*)f_Mu_ISO_Map->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

      TFile *f_triggeree = new TFile("input/TriggerSF_ee2017_pt.root");
      sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
      TFile *f_triggeremu = new TFile("input/TriggerSF_emu2017_pt.root");
      sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("h_lep1Pt_lep2Pt_Step3");
      TFile *f_triggermumu = new TFile("input/TriggerSF_mumu2017_pt.root");
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
      TFile *f_Ele_Reco_Map = new TFile("input/egammaEffi.txt_EGM2D_updatedAll.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("input/2018_ElectronTight.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map = new TFile("input/2018_RunABCD_SF_ID.root");
      sf_Mu_ID_H = *(TH2F*)f_Mu_ID_Map->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta");

      TFile *f_Mu_ISO_Map = new TFile("input/2018_RunABCD_SF_ISO.root");
      sf_Mu_ISO_H = *(TH2F*)f_Mu_ISO_Map->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

      TFile *f_triggeree = new TFile("input/TriggerSF_ee2018_pt.root");
      sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
      TFile *f_triggeremu = new TFile("input/TriggerSF_emu2018_pt.root");
      sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("h_lep1Pt_lep2Pt_Step3");
      TFile *f_triggermumu = new TFile("input/TriggerSF_mumu2018_pt.root");
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
  int nLeptonCut;//# of leptons we require
  bool VRveto;//When f=0, we require same-sign lepton pair. When f>0, we require three same-flavor leptons.
  int anti_index;//This is useful only when there are two tight leptons and one anti-selected lepton.
  int etabin1;//eta region of the first fakeable lepton
  int etabin2;//eta region of the second fakeable lepton
  

  if (fname.Contains("TTTo2L2Nu")) ifTopPt=true;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  Long64_t ntr = fChain->GetEntries ();
for (int f=0;f<2;f++){//f=0:MR;f=1:VR+AR1+AR2
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//  for (Long64_t jentry=0; jentry<100;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    displayProgress(jentry+f*ntr, 2*ntr) ;

    triggerPassEE = false;
    triggerPassEMu = false;
    triggerPassMuMu = false;
    triggerPassSingleE = false;
    triggerPassSingleMu = false;
    metFilterPass = false;
    ch =10;
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
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1)
        metFilterPass = true;
        }
        if(year == "2017"){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1 && trig_Flag_ecalBadCalibReduced ==1)
        metFilterPass = true;
        }
    }

    if(data == "data"){
      if(year == "2016" || year == "2018"){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1 &&  trig_Flag_eeBadScFilter_accept==1)
        metFilterPass = true;
        }
        if(year == "2017"){
        if ( trig_Flag_goodVertices_accept==1  &&  trig_Flag_globalSuperTightHalo2016Filter_accept==1 && trig_Flag_HBHENoiseFilter_accept==1 &&  trig_Flag_HBHENoiseIsoFilter_accept==1 && trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept==1 && trig_Flag_BadPFMuonFilter_accept==1 && trig_Flag_ecalBadCalibReduced ==1)
        metFilterPass = true;
        }
    }

//trigger
////MC
      if(data == "mc" && year == "2016"){
          if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele27_WPTight_Gsf_accept ) triggerPassEE =true;
          if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele27_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept) triggerPassEMu =true;
          if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept || trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept) triggerPassMuMu =true;
          if(trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassSingleE =true;
          if(trig_HLT_IsoMu24_accept) triggerPassSingleMu =true;
      }
      
      if(data == "mc" && year == "2017"){
          if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept) triggerPassEE =true;
          if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept || trig_HLT_IsoMu27_accept) triggerPassEMu =true;
          if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept || trig_HLT_IsoMu27_accept) triggerPassMuMu =true;
          if(trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassSingleE =true;
          if(trig_HLT_IsoMu24_accept) triggerPassSingleMu =true;
      }
      
      if(data == "mc" && year == "2018"){
          if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEE =true;
          if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept) triggerPassEMu =true;
          if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept || trig_HLT_IsoMu24_accept) triggerPassMuMu =true;
          if(trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassSingleE =true;
          if(trig_HLT_IsoMu24_accept) triggerPassSingleMu =true;
      }
  
////DATA
      if(data == "data"){
          if(year == "2016"){
              if(run == "H"){
                  if(dataset=="MuonEG"){
                      if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEMu =true;
                      if(!(trig_HLT_Ele32_WPTight_Gsf_accept||trig_HLT_Ele27_WPTight_Gsf_accept)&&trig_HLT_Ele20_WPTight_Gsf_accept) {
                          triggerPassSingleE=true;
                      }
                      if(!(trig_HLT_IsoMu30_accept||trig_HLT_IsoMu24_accept)&&trig_HLT_IsoMu20_accept) {
                          triggerPassSingleMu =true;
                      }
                  }
                  if(dataset=="SingleElectron"){
                      if(!(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEE =true;
                      if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEMu =true;
                      if(trig_HLT_Ele32_WPTight_Gsf_accept){
                          triggerPassSingleE=true;
                      }
                  }
                  if(dataset=="SingleMuon"){
                      if(!(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassMuMu =true;
                      if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele27_WPTight_Gsf_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassEMu =true;
                      if(trig_HLT_IsoMu30_accept){
                          triggerPassSingleMu=true;
                      }
                  }
                  if(dataset=="DoubleEG"){
                      if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEE =true;

                  }
                  if(dataset=="DoubleMu"){
                      if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) triggerPassMuMu =true;

                  }
              }
              if(run != "H"){
                  if(dataset=="MuonEG"){
                      if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) triggerPassEMu =true;
                      if(!(trig_HLT_Ele32_WPTight_Gsf_accept||trig_HLT_Ele27_WPTight_Gsf_accept)&&trig_HLT_Ele20_WPTight_Gsf_accept) {
                          triggerPassSingleE=true;
                      }
                      if(!(trig_HLT_IsoMu30_accept||trig_HLT_IsoMu24_accept)&&trig_HLT_IsoMu20_accept) {
                          triggerPassSingleMu =true;
                      }
                  }
                  if(dataset=="SingleElectron"){
                      if(!(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEE =true;
                      if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEMu =true;
                      if(trig_HLT_Ele32_WPTight_Gsf_accept){
                          triggerPassSingleE=true;
                      }
                  }
                  if(dataset=="SingleMuon"){
                      if(!(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassMuMu =true;
                      if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele27_WPTight_Gsf_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassEMu =true;
                      if(trig_HLT_IsoMu30_accept){
                          triggerPassSingleMu=true;
                      }
                  }
                  if(dataset=="DoubleEG"){
                      if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEE =true;
                      if(!trig_HLT_Ele32_WPTight_Gsf_accept&&trig_HLT_Ele27_WPTight_Gsf_accept){
                          triggerPassSingleE=true;
                      }
                  }
                  if(dataset=="DoubleMu"){
                      if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) triggerPassMuMu =true;
                      if(!trig_HLT_IsoMu30_accept&&trig_HLT_IsoMu24_accept){
                          triggerPassSingleMu=true;
                      }
                  }
              }
          }
          if(year == "2017"){
              if(dataset=="MuonEG"){
                  if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept ) triggerPassEMu =true;
                  if(!(trig_HLT_Ele32_WPTight_Gsf_accept||trig_HLT_Ele27_WPTight_Gsf_accept)&&trig_HLT_Ele20_WPTight_Gsf_accept) {
                      triggerPassSingleE=true;
                  }
                  if(!(trig_HLT_IsoMu30_accept||trig_HLT_IsoMu24_accept)&&trig_HLT_IsoMu20_accept) {
                      triggerPassSingleMu =true;
                  }
              }
              if(dataset=="SingleElectron"){
                  if(!trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept && trig_HLT_Ele35_WPTight_Gsf_accept) triggerPassEE =true;
                  if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && trig_HLT_Ele35_WPTight_Gsf_accept) triggerPassEMu =true;
                  if(trig_HLT_Ele32_WPTight_Gsf_accept){
                      triggerPassSingleE=true;
                  }
              }
              if(dataset=="SingleMuon"){
                  if(!trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept && trig_HLT_IsoMu27_accept) triggerPassMuMu =true;
                  if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept) && trig_HLT_IsoMu27_accept) triggerPassEMu =true;
                  if(trig_HLT_IsoMu30_accept){
                      triggerPassSingleMu=true;
                  }
              }
              if(dataset=="DoubleEG"){
                  if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept) triggerPassEE =true;
                  if(!trig_HLT_Ele32_WPTight_Gsf_accept&&trig_HLT_Ele27_WPTight_Gsf_accept){
                      triggerPassSingleE=true;
                  }
              }
              if(dataset=="DoubleMu"){
                  if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept)triggerPassMuMu =true;
                  if(!trig_HLT_IsoMu30_accept&&trig_HLT_IsoMu24_accept){
                      triggerPassSingleMu=true;
                  }
              }
          }
          if(year == "2018"){
              if(dataset=="MuonEG"){
                  if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) triggerPassEMu =true;
                  if(!(trig_HLT_Ele32_WPTight_Gsf_accept||trig_HLT_Ele27_WPTight_Gsf_accept)&&trig_HLT_Ele20_WPTight_Gsf_accept) {
                      triggerPassSingleE=true;
                  }
                  if(!(trig_HLT_IsoMu30_accept||trig_HLT_IsoMu24_accept)&&trig_HLT_IsoMu20_accept) {
                      triggerPassSingleMu =true;
                  }
              }
              if(dataset=="EGamma"){
                  if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEE =true;
                  if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEMu =true;
                  if(trig_HLT_Ele32_WPTight_Gsf_accept){
                      triggerPassSingleE=true;
                  }
              }
              if(dataset=="SingleMuon"){
                  if(!trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept && trig_HLT_IsoMu24_accept) triggerPassMuMu =true;
                  if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) && trig_HLT_IsoMu24_accept) triggerPassEMu =true;
                  if(trig_HLT_IsoMu30_accept){
                      triggerPassSingleMu=true;
                  }
              }
              if(dataset=="DoubleMu"){
                  if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept) triggerPassMuMu =true;
                  if(!trig_HLT_IsoMu30_accept&&trig_HLT_IsoMu24_accept){
                      triggerPassSingleMu=true;
                  }
              }
          }
      }
 


//cout<<ev_event<<"  "<<trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept <<"  "<< trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept <<"  "<< trig_HLT_Ele27_WPTight_Gsf_accept  <<"  "<<trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept <<"  "<< trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept <<"  "<< trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept<<"  "<<trig_HLT_IsoMu24_accept<<"  "<<trig_HLT_IsoTkMu24_accept<<endl;
//cout<<"trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept "<< "trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept "<< "trig_HLT_Ele27_WPTight_Gsf_accept " <<"trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept "<< "trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept "<< "trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept "  <<"trig_HLT_IsoMu24_accept "<<"trig_HLT_IsoTkMu24_accept"<<endl;
    if(!(triggerPassSingleE || triggerPassSingleMu)&&f==0) continue;
    if(!(triggerPassEE || triggerPassEMu || triggerPassMuMu)&&f>0) continue;
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
  for (int l=0;l<mu_gt_pt->size();l++){
      if(data == "data"){
          muPtSFRochester = rc.kScaleDT((*mu_gt_charge)[l], (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l], 0, 0);
      }
      if (data == "mc"){
          if ((*mu_mc_index)[l]!=-1 && abs((*mc_pdgId)[(*mu_mc_index)[l]]) == 13) muPtSFRochester = rc.kSpreadMC((*mu_gt_charge)[l], (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l], (*mc_pt)[(*mu_mc_index)[l]],0, 0);
          if ((*mu_mc_index)[l]<0) muPtSFRochester = rc.kSmearMC((*mu_gt_charge)[l], (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l], (*mu_trackerLayersWithMeasurement)[l] , gRandom->Rndm(),0, 0);
      }
      if(muPtSFRochester * (*mu_gt_pt)[l] <20 || abs((*mu_gt_eta)[l]) > 2.4) continue;
      if(!(*mu_CutBasedIdMedium)[l]||(*mu_pfIsoDbCorrected04)[l] > 0.15) continue;
      //if((*mu_pfIsoDbCorrected04)[l] > 0.15 && (*mu_pfIsoDbCorrected04)[l] < 0.25) continue;
      nLoose++;
      if (f==0){
          if (nTight==1){
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l],(*mu_gt_charge)[l],l,10));
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              if ((*mu_MvaMedium)[l]) isTight = true;
              continue;
          }
          if ((*mu_MvaMedium)[l]){
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l],(*mu_gt_charge)[l],l,10));
              (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              nTight++;
          }
      }
      else{
          if ((*mu_MvaMedium)[l]){
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l],(*mu_gt_charge)[l],l,10));
              (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              nTight++;
          }
          else{
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l],(*mu_gt_charge)[l],l,10));
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
          }
      }
  }
// electron
    for (int l=0;l<gsf_pt->size();l++){
      elePt = (*gsf_ecalTrkEnergyPostCorr)[l]*sin(2.*atan(exp(-1.*(*gsf_eta)[l]))) ;
      if(elePt <20 || abs((*gsf_eta)[l]) > 2.4 || (abs((*gsf_sc_eta)[l])> 1.4442 && (abs((*gsf_sc_eta)[l])< 1.566))) continue;
      if(!(*gsf_VID_mvaEleID_Fall17_noIso_V2_wp80)[l]) continue;
      nLoose++;
      if (f==0){
         if (nTight==1){
         selectedLeptons->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
         if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
         if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
         if ((*gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight)[l]) isTight = true;
         continue;
         }
         if ((*gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight)[l]){
            selectedLeptons->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
            (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
            if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
            if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
            nTight++;
         }
      }
      else{
      if ((*gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight)[l]){
         selectedLeptons->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
         (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
         if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
         if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
         nTight++;
      }
      else{
      selectedLeptons->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
      if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
      if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
      }
      }
    }
    if (f==0&&nLoose==2){
      sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
      if (selectedLeptons->size()!=2 ||
         (*selectedLeptons)[0]->charge_+(*selectedLeptons)[1]->charge_==0 ||
         (*selectedLeptons)[0]->pt_ <30||
         !((*selectedLeptons)[0]->isTag)){
          for (int l=0;l<selectedLeptons->size();l++){
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
      for (int l=gsf_pt->size()-1;l>-1;l--){
          elePt = (*gsf_ecalTrkEnergyPostCorr)[l]*sin(2.*atan(exp(-1.*(*gsf_eta)[l]))) ;
          if(elePt <20 || abs((*gsf_eta)[l]) > 2.4 || (abs((*gsf_sc_eta)[l])> 1.4442 && (abs((*gsf_sc_eta)[l])< 1.566))) continue;
          if(!(*gsf_VID_mvaEleID_Fall17_noIso_V2_wp80)[l]) continue;
          if (nTight==1){
             selectedLeptons->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
             if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
             if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
             if ((*gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight)[l]) isTight = true;
             continue;
          }
          if ((*gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight)[l]){
             selectedLeptons->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
             (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
             if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
             if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
             nTight++;
          }
      }
      // Muon
      for (int l=mu_gt_pt->size()-1;l>-1;l--){
          if(data == "data"){
              muPtSFRochester = rc.kScaleDT((*mu_gt_charge)[l], (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l], 0, 0);
          }
          if (data == "mc"){
              if ((*mu_mc_index)[l]!=-1 && abs((*mc_pdgId)[(*mu_mc_index)[l]]) == 13) muPtSFRochester = rc.kSpreadMC((*mu_gt_charge)[l], (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l], (*mc_pt)[(*mu_mc_index)[l]],0, 0);
              if ((*mu_mc_index)[l]<0) muPtSFRochester = rc.kSmearMC((*mu_gt_charge)[l], (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l], (*mu_trackerLayersWithMeasurement)[l] , gRandom->Rndm(),0, 0);
          }
          if(muPtSFRochester * (*mu_gt_pt)[l] <20 || abs((*mu_gt_eta)[l]) > 2.4) continue;
          if(!(*mu_CutBasedIdMedium)[l]||(*mu_pfIsoDbCorrected04)[l] > 0.15) continue;
          //if((*mu_pfIsoDbCorrected04)[l] > 0.15 && (*mu_pfIsoDbCorrected04)[l] < 0.25) continue;
          if (nTight==1){
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l],(*mu_gt_charge)[l],l,10));
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              if ((*mu_MvaMedium)[l]) isTight = true;
              continue;
          }
          if ((*mu_MvaMedium)[l]){
              selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l],(*mu_gt_charge)[l],l,10));
              (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
              if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
              if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
              nTight++;
          }
      }
    }
    }
      if (selectedLeptons->size()==3 && abs((*selectedLeptons)[0]->charge_+(*selectedLeptons)[1]->charge_+(*selectedLeptons)[2]->charge_)==1){
          if (((*selectedLeptons)[0]->charge_+(*selectedLeptons)[1]->charge_)==0){
              Zmass=((*selectedLeptons)[0]->p4_+(*selectedLeptons)[1]->p4_).M();
          }
          else{
              Zmass=((*selectedLeptons)[0]->p4_+(*selectedLeptons)[2]->p4_).M();
          }
      }
      if (Zmass>76&&Zmass<106) OnZ=true;
    sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
    nLeptonCut=f<1?2:3;
    if (selectedLeptons->size()==2){
        if ((*selectedLeptons)[0]->charge_+(*selectedLeptons)[1]->charge_==0||nLoose!=2) VRveto=true;
    }
    if (selectedLeptons->size()==3){
        if ((*selectedLeptons)[0]->lep_+(*selectedLeptons)[1]->lep_+(*selectedLeptons)[2]->lep_>3&&
            (*selectedLeptons)[0]->lep_+(*selectedLeptons)[1]->lep_+(*selectedLeptons)[2]->lep_<30) VRveto=true;
    }
// dilepton/trilepton selection
//cout<<ev_event<<"  "<<triggerPass<<"  "<<metFilterPass<<"  "<<selectedLeptons->size()<<endl;
    if(selectedLeptons->size()!=nLeptonCut ||
      ((*selectedLeptons)[0]->pt_ <30) ||
       VRveto || !((*selectedLeptons)[0]->isTag)) {
      for (int l=0;l<selectedLeptons->size();l++){
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
//    if(ch==1) etabin1=0;
//    if(ch==1) etabin2=0;
      
    if (f==0){
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
    for (int l=0;l<jet_pt->size();l++){
      if(data == "mc" && ((*jet_Smeared_pt)[l] <30 || abs((*jet_eta)[l]) > 2.4)) continue;
      if(data == "data" && ((*jet_pt)[l] <30 || abs((*jet_eta)[l]) > 2.4)) continue;
      if(year == "2016" && !(*jet_isJetIDTightLepVeto_2016)[l]) continue;
      if(year == "2017" && !(*jet_isJetIDLepVeto_2017)[l]) continue;
      if(year == "2018" && !(*jet_isJetIDLepVeto_2018)[l]) continue;
      jetlepfail = false;
      for (int i=0;i<selectedLeptons->size();i++){
        if(deltaR((*selectedLeptons)[i]->eta_,(*selectedLeptons)[i]->phi_,(*jet_eta)[l],(*jet_phi)[l]) < 0.4 ) jetlepfail=true;
      }
      if(jetlepfail) continue; 
      if(data == "mc"){
        selectedJets->push_back(new jet_candidate((*jet_Smeared_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,(*jet_partonFlavour)[l]));
      }
      if(data == "data"){
        selectedJets->push_back(new jet_candidate((*jet_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l],year,0));
      }
    }

    sort(selectedJets->begin(), selectedJets->end(), ComparePtJet);
    nbjet=0;
    for (int l=0;l<selectedJets->size();l++){
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
    if (data == "mc" && year == "2016") weight_PU = wPU.PU_2016(mc_trueNumInteractions,"nominal");
    if (data == "mc" && year == "2017") weight_PU = wPU.PU_2017(mc_trueNumInteractions,"nominal");
    if (data == "mc" && year == "2018") weight_PU = wPU.PU_2018(mc_trueNumInteractions,"nominal");
    if (data == "mc") weight_Lumi = (1000*xs*lumi)/Nevent;
    if (data == "mc" && (year == "2016" || year == "2017")) weight_prefiring = ev_prefiringweight;
    if (data == "mc" && ifTopPt) {
      for (int l=0;l<mc_status->size();l++){
        if((*mc_status)[l]<30 && (*mc_status)[l]>20 && (*mc_pdgId)[l]==24) wp.SetPtEtaPhiE((*mc_pt)[l], (*mc_eta)[l], (*mc_phi)[l], (*mc_energy)[l]) ;
        if((*mc_status)[l]<30 && (*mc_status)[l]>20 && (*mc_pdgId)[l]==-24) wm.SetPtEtaPhiE((*mc_pt)[l], (*mc_eta)[l], (*mc_phi)[l], (*mc_energy)[l]) ;
        if((*mc_status)[l]<30 && (*mc_status)[l]>20 && (*mc_pdgId)[l]==5) b.SetPtEtaPhiE((*mc_pt)[l], (*mc_eta)[l], (*mc_phi)[l], (*mc_energy)[l]) ;
        if((*mc_status)[l]<30 && (*mc_status)[l]>20 && (*mc_pdgId)[l]==-5) ab.SetPtEtaPhiE((*mc_pt)[l], (*mc_eta)[l], (*mc_phi)[l], (*mc_energy)[l]) ;
      }
    weight_topPt = sqrt(topPt((wp + b).Pt()) * topPt((wm + ab).Pt()));
    }

    if (data == "mc") weight_lep = sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_Lumi  * mc_w_sign * weight_prefiring * weight_topPt;
    if (data == "mc") weight_lepB = sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_Lumi * mc_w_sign *  weight_prefiring * weight_topPt * (P_bjet_data/P_bjet_mc);
//     cout<<ev_event<<"   "<<sf_Ele_Reco<<"   "<<sf_Ele_ID<<"      "<<sf_Mu_ID<<"   "<<sf_Mu_ISO<<"   "<<sf_Trigger<<"   "<<weight_PU<<endl;
//    if(selectedJets->size()<3 || MET_FinalCollection_Pt>30 || nbjet !=1) continue;
      
    //          Debug
//    cout<<endl<<" ch= "<<ch<<endl;
//    cout<<" f = "<<f<<" nTight = "<<nTight<<" isTight = "<<isTight<<" anti-index = "<<anti_index<<endl;
//    cout<<" etabin1 = "<<etabin1<<" etabin2 = "<<etabin2<<endl;
//    cout<<" 1st lep eta = "<<(*selectedLeptons)[0]->eta_<<" 2nd lep = "<<(*selectedLeptons)[1]->eta_<<" 3rd lep = "<<(*selectedLeptons)[f<1?1:2]->eta_<<endl;
//    cout<<" 1st lep flavor = "<<(*selectedLeptons)[0]->lep_<<" 2nd lep = "<<(*selectedLeptons)[1]->lep_<<" 3rd lep = "<<(*selectedLeptons)[f<1?1:2]->lep_<<endl;
//    cout<<" 1st lep charge = "<<(*selectedLeptons)[0]->charge_<<" 2nd lep = "<<(*selectedLeptons)[1]->charge_<<" 3rd lep = "<<(*selectedLeptons)[f<1?1:2]->charge_<<endl;
//    cout<<" 1st lep isTag = "<<(*selectedLeptons)[0]->isTag<<" 2nd lep = "<<(*selectedLeptons)[1]->isTag<<" 3rd lep = "<<(*selectedLeptons)[f<1?1:2]->isTag<<endl;

 if(f==0){
    HistsMR[ch][0][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    HistsMR[ch][0][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    HistsMR[ch][0][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    if(isTight) HistsMR[ch][0][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    if(isTight) HistsMR[ch][0][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    if(isTight) HistsMR[ch][0][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
      
    if (MET_FinalCollection_Pt>20){
      HistsMR[ch][1][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
      HistsMR[ch][1][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
      HistsMR[ch][1][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
      if(isTight) HistsMR[ch][1][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
      if(isTight) HistsMR[ch][1][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
      if(isTight) HistsMR[ch][1][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    }
      
    if (MET_FinalCollection_Pt>20&&selectedJets->size()>=1&&nbjet<=1){
      HistsMR[ch][2][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
      HistsMR[ch][2][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
      HistsMR[ch][2][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
      if(isTight) HistsMR[ch][2][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
      if(isTight) HistsMR[ch][2][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
      if(isTight) HistsMR[ch][2][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    }
      
    if (MET_FinalCollection_Pt>20&&selectedJets->size()>=2&&nbjet<=1){
      HistsMR[ch][3][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsMR[ch][3][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsMR[ch][3][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      if(isTight) HistsMR[ch][3][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      if(isTight) HistsMR[ch][3][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      if(isTight) HistsMR[ch][3][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    }
     
    if (MET_FinalCollection_Pt>20 && nbjet >=1){
      HistsMR[ch][4][etabin1][0]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      HistsMR[ch][4][etabin1][1]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      HistsMR[ch][4][etabin1][2]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
      if(isTight) HistsMR[ch][4][etabin1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
      if(isTight) HistsMR[ch][4][etabin1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
      if(isTight) HistsMR[ch][4][etabin1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    }

     nAccept++;
  }
  else{
      
  HistsFF[3-nTight][ch][0][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][0][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][0][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][0][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][0][4]->Fill(selectedJets->size(),weight_lepB);
  HistsFF[3-nTight][ch][0][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][0][6]->Fill(MET_FinalCollection_Pt,weight_lep);
  HistsFF[3-nTight][ch][0][7]->Fill(pv_n,weight_lep);
    
  if (MET_FinalCollection_Pt<20){
  HistsFF[3-nTight][ch][1][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][1][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][1][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][1][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][1][4]->Fill(selectedJets->size(),weight_lep);
  HistsFF[3-nTight][ch][1][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][1][6]->Fill(MET_FinalCollection_Pt,weight_lep);
  HistsFF[3-nTight][ch][1][7]->Fill(pv_n,weight_lep);
  }
  else{
  HistsFF[3-nTight][ch][2][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][2][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][2][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][2][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][2][4]->Fill(selectedJets->size(),weight_lep);
  HistsFF[3-nTight][ch][2][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][2][6]->Fill(MET_FinalCollection_Pt,weight_lep);
  HistsFF[3-nTight][ch][2][7]->Fill(pv_n,weight_lep);
  }
      
  if (OnZ){
  HistsFF[3-nTight][ch][3][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][3][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][3][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][3][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][3][4]->Fill(selectedJets->size(),weight_lep);
  HistsFF[3-nTight][ch][3][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][3][6]->Fill(MET_FinalCollection_Pt,weight_lep);
  HistsFF[3-nTight][ch][3][7]->Fill(pv_n,weight_lep);
  }
  else{
  HistsFF[3-nTight][ch][4][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][4][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][4][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][4][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][4][4]->Fill(selectedJets->size(),weight_lep);
  HistsFF[3-nTight][ch][4][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][4][6]->Fill(MET_FinalCollection_Pt,weight_lep);
  HistsFF[3-nTight][ch][4][7]->Fill(pv_n,weight_lep);
  }
      
  if (MET_FinalCollection_Pt>20&&!OnZ&&nbjet<=1&&selectedJets->size()>=1){
  HistsFF[3-nTight][ch][5][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
  HistsFF[3-nTight][ch][5][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][5][2]->Fill((*selectedJets)[0]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsFF[3-nTight][ch][5][3]->Fill((*selectedJets)[0]->eta_,weight_lep);
  HistsFF[3-nTight][ch][5][4]->Fill(selectedJets->size(),weight_lep);
  HistsFF[3-nTight][ch][5][5]->Fill(nbjet,weight_lepB);
  HistsFF[3-nTight][ch][5][6]->Fill(MET_FinalCollection_Pt,weight_lep);
  HistsFF[3-nTight][ch][5][7]->Fill(pv_n,weight_lep);
  }
  //AR1:ll\overline{l}
  if (nTight==2){
      
  HistsAR1[ch][0][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][0][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][0][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][0][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][0][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][0][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][0][etabin1][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][0][etabin1][7]->Fill(pv_n,(*selectedLeptons)[anti_index]->pt_,weight_lep);
    
  if (MET_FinalCollection_Pt<20){
  HistsAR1[ch][1][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][1][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][1][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][1][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][1][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][1][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][1][etabin1][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][1][etabin1][7]->Fill(pv_n,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  }
  else{
  HistsAR1[ch][2][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][2][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][2][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][2][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][2][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][2][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][2][etabin1][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][2][etabin1][7]->Fill(pv_n,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  }
      
  if (OnZ){
  HistsAR1[ch][3][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][3][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][3][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][3][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][3][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][3][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][3][etabin1][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][3][etabin1][7]->Fill(pv_n,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  }
  else{
  HistsAR1[ch][4][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][4][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][4][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][4][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][4][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][4][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][4][etabin1][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][4][etabin1][7]->Fill(pv_n,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  }
      
  if (MET_FinalCollection_Pt>20&&!OnZ&&nbjet<=1&&selectedJets->size()>=1){
  HistsAR1[ch][5][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][5][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][5][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR1[ch][5][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][5][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][5][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
  HistsAR1[ch][5][etabin1][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[anti_index]->pt_,weight_lep);
  HistsAR1[ch][5][etabin1][7]->Fill(pv_n,(*selectedLeptons)[anti_index]->pt_,weight_lep);
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
  HistsAR2[ch][0][etabin1][etabin2][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][0][etabin1][etabin2][7]->Fill(pv_n,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  
  if (MET_FinalCollection_Pt<20){
  HistsAR2[ch][1][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][1][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][1][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][1][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][1][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][1][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][1][etabin1][etabin2][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][1][etabin1][etabin2][7]->Fill(pv_n,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  }
  else{
  HistsAR2[ch][2][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][2][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][2][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][2][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][2][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][2][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][2][etabin1][etabin2][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][2][etabin1][etabin2][7]->Fill(pv_n,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  }
      
  if (OnZ){
  HistsAR2[ch][3][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][3][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][3][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][3][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][3][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][3][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][3][etabin1][etabin2][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][3][etabin1][etabin2][7]->Fill(pv_n,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  }
  else{
  HistsAR2[ch][4][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][4][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][4][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][4][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][4][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][4][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][4][etabin1][etabin2][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][4][etabin1][etabin2][7]->Fill(pv_n,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  }
      
  if (MET_FinalCollection_Pt>20&&!OnZ&&nbjet<=1&&selectedJets->size()>=1){
  HistsAR2[ch][5][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][5][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][5][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  if(selectedJets->size()>0) HistsAR2[ch][5][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][5][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][5][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
  HistsAR2[ch][5][etabin1][etabin2][6]->Fill(MET_FinalCollection_Pt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  HistsAR2[ch][5][etabin1][etabin2][7]->Fill(pv_n,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lep);
  }
      
  }
      
 }

    for (int l=0;l<selectedLeptons->size();l++){
      delete (*selectedLeptons)[l];
    }
    for (int l=0;l<selectedJets->size();l++){
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

  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int j=0;j<etaregs.size();++j){
        for (int l=0;l<vars.size();++l){
            HistsMR[i][k][j][l]  ->Write("",TObject::kOverwrite);
        }
      }
    }
  }
    
  for (int j=0;j<FF.size();++j){
    for (int i=0;i<channelsFF.size();++i){
        for (int k=0;k<regionsFF.size();++k){
            for (int l=0;l<varsFF.size();++l){
                HistsFF[j][i][k][l]  ->Write("",TObject::kOverwrite);
            }
        }
    }
  }

  for (int i=0;i<channelsFF.size();++i){
      for (int k=0;k<regionsFF.size();++k){
          for (int j=0;j<etaregs.size();++j){
              for (int l=0;l<varsFF.size();++l){
                  HistsAR1[i][k][j][l]  ->Write("",TObject::kOverwrite);
                  for (int m=0;m<etaregs.size();++m){
                      HistsAR2[i][k][j][m][l]  ->Write("",TObject::kOverwrite);
                 }
              }
          }
      }
  }
    
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int j=0;j<etaregs.size();++j){
          for (int l=0;l<vars.size();++l){
              delete HistsMR[i][k][j][l];
          }
      }
    }
  }
    
  for (int j=0;j<FF.size();++j){
      for (int i=0;i<channelsFF.size();++i){
          for (int k=0;k<regionsFF.size();++k){
              for (int l=0;l<varsFF.size();++l){
              delete HistsFF[j][i][k][l];
              }
          }
      }
  }
    
  for (int i=0;i<channelsFF.size();++i){
      for (int k=0;k<regionsFF.size();++k){
         for (int j=0;j<etaregs.size();++j){
             for (int l=0;l<varsFF.size();++l){
                 delete HistsAR1[i][k][j][l];
                 for (int m=0;m<etaregs.size();++m){
                     delete HistsAR2[i][k][j][m][l];
                 }
             }
          }
      }
  }
    
  file_out.Close() ;
}
