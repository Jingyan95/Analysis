#define MyAnalysis_cxx
#include "MyAnalysis.h"
#include "PU_reWeighting.h"
#include "lepton_candidate.h"
#include "jet_candidate.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <TH2.h>
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
bool CompareFlavourLep(lepton_candidate *a, lepton_candidate *b) { return a->lep_ < b->lep_; }
bool CompareBaLep(lepton_candidate *a, lepton_candidate *b) { return a->isbalep < b->isbalep; }
bool ComparePtJet(jet_candidate *a, jet_candidate *b) { return a->pt_ > b->pt_; }
bool CompareBtagJet(jet_candidate *a, jet_candidate *b) { return a->bt_ > b->bt_; }

float getTopmass(lepton_candidate *a, jet_candidate *b, float MET, float phi){
    float mW=80.2;
    float Topmass=-99999;
    TLorentzVector n;
    TLorentzVector l=a->p4_;
    TLorentzVector bjet=b->p4_;
    float pz;
    float plx=l.Px();
    float ply=l.Py();
    float plz=l.Pz();
    float El=l.E();
    float x=mW*mW-El*El+plx*plx+ply*ply+plz*plz+2*plx*MET*sin(phi)+2*ply*MET*cos(phi);
    float A=4*(El*El-plz*plz);
    float B=-4*x*plz;
    float C=4*El*El*MET*MET-x*x;
    float delta=B*B-4*A*C; //quadratic formula
    if(delta<0){
        pz=-B/(2*A);
        n.SetPxPyPzE(MET*sin(phi),MET*cos(phi),pz,sqrt(MET*MET+pz*pz));
        Topmass=(n+l+bjet).M();
    }
    else{
       float sol1=(-B-sqrt(delta))/(2*A);
       float sol2=(-B+sqrt(delta))/(2*A);
       if(abs(sol1-plz)<abs(sol2-plz)){
       pz=sol1; //pick the one closest to lepton pz
        }
       else{
       pz=sol2;
       }
       n.SetPxPyPzE(MET*sin(phi),MET*cos(phi),pz,sqrt(MET*MET+pz*pz));
       Topmass=(n+l+bjet).M();
    }
    return Topmass;
}

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

   Double_t ptBins[11] = {30., 40., 60., 80., 100., 150., 200., 300., 400., 500., 1000.};
   Double_t etaBins [4]= {0., 0.6, 1.2, 2.4};
   TH2D *h2_BTaggingEff_Denom_b    = new TH2D("h2_BTaggingEff_Denom_b"   , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Denom_c    = new TH2D("h2_BTaggingEff_Denom_c"   , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Denom_udsg = new TH2D("h2_BTaggingEff_Denom_udsg", ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Num_b      = new TH2D("h2_BTaggingEff_Num_b"     , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Num_c      = new TH2D("h2_BTaggingEff_Num_c"     , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins);
   TH2D *h2_BTaggingEff_Num_udsg   = new TH2D("h2_BTaggingEff_Num_udsg"  , ";p_{T} [GeV];#eta", 10 , ptBins, 3 , etaBins); 


  typedef vector<TH1F*> Dim1;
  typedef vector<Dim1> Dim2;
  typedef vector<Dim2> Dim3;

  std::vector<TString> regions{"lll","lllOffZ","lllB1", "lllBg1", "lllMetl20", "lllMetg20", "lllMetl20Jet1B1", "lllMetl20Jetg1B1", "lllMetl20Jetg2Bg1", "lllMetg20Jet1B1", "lllMetg20Jetg1B1", "lllMetg20Jetg2Bg1"};
  std::vector<TString> channels{"eee", "emul", "mumumu"};
  std::vector<TString> vars   {"lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","lep3Pt","lep3Eta","lep3Phi",
        "LFVePt","LFVeEta","LFVePhi","LFVmuPt","LFVmuEta","LFVmuPhi","balPt","balEta","balPhi","Topmass",
        "llM","llPt","llDr","llDphi","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx", "llMZw"};

  std::vector<int>    nbins   {30      ,20       ,25       ,20      ,20       ,25       ,15      ,20       ,25       ,
        20      ,20       ,25       ,20       ,20        ,25        ,20     ,20      ,25      ,12       ,
        30   ,20    ,25    ,15      ,20      ,20       ,25       ,10    ,6      ,30   ,20      ,70    ,80};

  std::vector<float> lowEdge  {0       ,-3       ,-4       ,0       ,-3       ,-4       ,0       ,-3       ,-4       ,
        0       ,-3       ,-4       ,0        ,-3        ,-4        ,0      ,-3      ,-4      ,100      ,
        0    ,0     ,0     ,0       ,0       ,-3       ,-4       ,0     ,0      ,0    ,-4      ,0     ,70};

  std::vector<float> highEdge {300     ,3        ,4        ,200     ,3        ,4        ,150     ,3        ,4        ,
        200     ,3        ,4        ,200      ,3         ,4         ,200    ,3       ,4       ,340      ,
        500  ,200   ,7     ,4       ,300     ,3        ,4        ,10    ,6      ,210  ,4       ,70    ,110};

  Dim3 Hists(channels.size(),Dim2(regions.size(),Dim1(vars.size())));  
  std::stringstream name;
  TH1F *h_test;
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l];
        h_test = new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]);
        h_test->StatOverflows(kTRUE);
        h_test->Sumw2(kTRUE);
        Hists[i][k][l] = h_test;
        name.str("");
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
  if(year == "2016")    rochesterFile = "../TopLFV/input/RoccoR2016.txt";
  if(year == "2017")    rochesterFile = "../TopLFV/input/RoccoR2017.txt";
  if(year == "2018")    rochesterFile = "../TopLFV/input/RoccoR2018.txt";
  rc.init(rochesterFile);

  if(data == "mc"){
    TFile *f_btagEff_Map = new TFile("../TopLFV/input/btagEff.root");
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

    if(year == "2016")    btagFile = "../TopLFV/input/DeepCSV_2016LegacySF_V1.csv";
    if(year == "2017")    btagFile = "../TopLFV/input/DeepCSV_94XSF_V4_B_F.csv";
    if(year == "2018")    btagFile = "../TopLFV/input/DeepCSV_102XSF_V1.csv";

    BTagCalibration calib("DeepCSV",btagFile);
    reader.load(calib,BTagEntry::FLAV_B,"comb"); 
    reader.load(calib,BTagEntry::FLAV_C,"comb");
    reader.load(calib,BTagEntry::FLAV_UDSG,"comb");

    if(year == "2016"){
      TFile *f_Ele_Reco_Map = new TFile("../TopLFV/input/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("../TopLFV/input/2016LegacyReReco_ElectronTight_Fall17V2.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map_1 = new TFile("../TopLFV/input/2016_RunBCDEF_SF_ID.root");
      TH2F *sf_Mu_ID_H_1 = (TH2F*)f_Mu_ID_Map_1->Get("NUM_TightID_DEN_genTracks_eta_pt");
      TFile *f_Mu_ID_Map_2 = new TFile("../TopLFV/input/2016_RunGH_SF_ID.root");
      TH2F *sf_Mu_ID_H_2 = (TH2F*)f_Mu_ID_Map_2->Get("NUM_TightID_DEN_genTracks_eta_pt");
      sf_Mu_ID_H_1->Scale(0.55);
      sf_Mu_ID_H_2->Scale(0.45);
      sf_Mu_ID_H_1->Add(sf_Mu_ID_H_2);
      sf_Mu_ID_H = *sf_Mu_ID_H_1;

      TFile *f_Mu_ISO_Map_1 = new TFile("../TopLFV/input/2016_RunBCDEF_SF_ISO.root");
      TH2F *sf_Mu_ISO_H_1 = (TH2F*)f_Mu_ISO_Map_1->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
      TFile *f_Mu_ISO_Map_2 = new TFile("../TopLFV/input/2016_RunGH_SF_ISO.root");
      TH2F *sf_Mu_ISO_H_2 = (TH2F*)f_Mu_ISO_Map_2->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
      sf_Mu_ISO_H_1->Scale(0.55);
      sf_Mu_ISO_H_2->Scale(0.45);
      sf_Mu_ISO_H_1->Add(sf_Mu_ISO_H_2);
      sf_Mu_ISO_H = *sf_Mu_ISO_H_1;

      TFile *f_triggeree = new TFile("../TopLFV/input/TriggerSF_ee2016_pt.root");
      sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
      TFile *f_triggeremu = new TFile("../TopLFV/input/TriggerSF_emu2016_pt.root");
      sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("h_lep1Pt_lep2Pt_Step3");
      TFile *f_triggermumu = new TFile("../TopLFV/input/TriggerSF_mumu2016_pt.root");
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
      TFile *f_Ele_Reco_Map = new TFile("../TopLFV/input/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("../TopLFV/input/2017_ElectronTight.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map = new TFile("../TopLFV/input/2017_RunBCDEF_SF_ID_syst.root");
      sf_Mu_ID_H = *(TH2F*)f_Mu_ID_Map->Get("NUM_TightID_DEN_genTracks_pt_abseta");

      TFile *f_Mu_ISO_Map = new TFile("../TopLFV/input/2017_RunBCDEF_SF_ISO_syst.root");
      sf_Mu_ISO_H = *(TH2F*)f_Mu_ISO_Map->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

      TFile *f_triggeree = new TFile("../TopLFV/input/TriggerSF_ee2017_pt.root");
      sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
      TFile *f_triggeremu = new TFile("../TopLFV/input/TriggerSF_emu2017_pt.root");
      sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("h_lep1Pt_lep2Pt_Step3");
      TFile *f_triggermumu = new TFile("../TopLFV/input/TriggerSF_mumu2017_pt.root");
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
      TFile *f_Ele_Reco_Map = new TFile("../TopLFV/input/egammaEffi.txt_EGM2D_updatedAll.root");
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile("../TopLFV/input/2018_ElectronTight.root");
      sf_Ele_ID_H = *(TH2F*)f_Ele_ID_Map->Get("EGamma_SF2D");

      TFile *f_Mu_ID_Map = new TFile("../TopLFV/input/2018_RunABCD_SF_ID.root");
      sf_Mu_ID_H = *(TH2F*)f_Mu_ID_Map->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta");

      TFile *f_Mu_ISO_Map = new TFile("../TopLFV/input/2018_RunABCD_SF_ISO.root");
      sf_Mu_ISO_H = *(TH2F*)f_Mu_ISO_Map->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

      TFile *f_triggeree = new TFile("../TopLFV/input/TriggerSF_ee2018_pt.root");
      sf_triggeree_H = *(TH2F*)f_triggeree->Get("h_lep1Pt_lep2Pt_Step6");
      TFile *f_triggeremu = new TFile("../TopLFV/input/TriggerSF_emu2018_pt.root");
      sf_triggeremu_H = *(TH2F*)f_triggeremu->Get("h_lep1Pt_lep2Pt_Step3");
      TFile *f_triggermumu = new TFile("../TopLFV/input/TriggerSF_mumu2018_pt.root");
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
  std::vector<lepton_candidate*> *selectedLeptons_copy;
  std::vector<jet_candidate*> *selectedJets;
  std::vector<jet_candidate*> *selectedJets_copy;
  TLorentzVector wp, wm, b, ab;
  bool triggerPassEE;
  bool triggerPassEMu;
  bool triggerPassMuMu;
  bool metFilterPass;
  bool ifTopPt=false;
  int ch;
  int ch1;
  bool compete=false;
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
  float t1,t2;
  float Topmass=0;
  float mT=172;

  if (fname.Contains("TTTo2L2Nu")) ifTopPt=true;

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  Long64_t ntr = fChain->GetEntries ();
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//  for (Long64_t jentry=0; jentry<100;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    displayProgress(jentry, ntr) ;

    triggerPassEE = false;
    triggerPassEMu = false;
    triggerPassMuMu = false;
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
      }

      if(data == "mc" && year == "2017"){
        if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept) triggerPassEE =true;
        if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept || trig_HLT_IsoMu27_accept) triggerPassEMu =true;
        if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept || trig_HLT_IsoMu27_accept) triggerPassMuMu =true;
      }

      if(data == "mc" && year == "2018"){
        if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEE =true;
        if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept) triggerPassEMu =true;
        if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept || trig_HLT_IsoMu24_accept) triggerPassMuMu =true;
      } 

////DATA
    if(data == "data"){
      if(year == "2016"){
        if(run == "H"){
          if(dataset=="MuonEG"){
            if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleElectron"){
            if(!(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEE =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleMuon"){
            if(!(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassMuMu =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele27_WPTight_Gsf_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassEMu =true; 
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
          }
          if(dataset=="SingleElectron"){
            if(!(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEE =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleMuon"){
            if(!(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassMuMu =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele27_WPTight_Gsf_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassEMu =true;
          }
          if(dataset=="DoubleEG"){
            if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEE =true;
          }
          if(dataset=="DoubleMu"){
            if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) triggerPassMuMu =true;
          }
        }
      }
      if(year == "2017"){
        if(dataset=="MuonEG"){
          if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept ) triggerPassEMu =true;
        }
        if(dataset=="SingleElectron"){
          if(!trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept && trig_HLT_Ele35_WPTight_Gsf_accept) triggerPassEE =true;
          if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && trig_HLT_Ele35_WPTight_Gsf_accept) triggerPassEMu =true;
        }
        if(dataset=="SingleMuon"){
           if(!trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept && trig_HLT_IsoMu27_accept) triggerPassMuMu =true;
          if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele35_WPTight_Gsf_accept) && trig_HLT_IsoMu27_accept) triggerPassEMu =true;
        }
        if(dataset=="DoubleEG"){
          if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept) triggerPassEE =true;
        }
        if(dataset=="DoubleMu"){
          if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_accept)triggerPassMuMu =true;
        }
      }
      if(year == "2018"){
        if(dataset=="MuonEG"){
          if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) triggerPassEMu =true;
        }
        if(dataset=="EGamma"){
          if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEE =true;
          if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept) && trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEMu =true;
      }
      if(dataset=="SingleMuon"){
        if(!trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept && trig_HLT_IsoMu24_accept) triggerPassMuMu =true;
        if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) && trig_HLT_IsoMu24_accept) triggerPassEMu =true;
      }
      if(dataset=="DoubleMu"){
        if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept) triggerPassMuMu =true;
      }
    }
  }
 


//cout<<ev_event<<"  "<<trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept <<"  "<< trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept <<"  "<< trig_HLT_Ele27_WPTight_Gsf_accept  <<"  "<<trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept <<"  "<< trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept <<"  "<< trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept<<"  "<<trig_HLT_IsoMu24_accept<<"  "<<trig_HLT_IsoTkMu24_accept<<endl;
//cout<<"trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept "<< "trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_accept "<< "trig_HLT_Ele27_WPTight_Gsf_accept " <<"trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept "<< "trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept "<< "trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept "  <<"trig_HLT_IsoMu24_accept "<<"trig_HLT_IsoTkMu24_accept"<<endl;
    if(!(triggerPassEE || triggerPassEMu || triggerPassMuMu)) continue;
    if(!metFilterPass) continue;

// lepton selection
  selectedLeptons = new std::vector<lepton_candidate*>();//typlical ordered by pT
  selectedLeptons_copy = new std::vector<lepton_candidate*>();// ordered by [e, mu , bachelor lepton ]
// electron
    for (int l=0;l<gsf_pt->size();l++){
      elePt = (*gsf_ecalTrkEnergyPostCorr)[l]*sin(2.*atan(exp(-1.*(*gsf_eta)[l]))) ;
      if(elePt <20 || abs((*gsf_eta)[l]) > 2.4 || (abs((*gsf_sc_eta)[l])> 1.4442 && (abs((*gsf_sc_eta)[l])< 1.566))) continue;
      if(!(*gsf_VID_cutBasedElectronID_Fall17_94X_V2_tight)[l]) continue;
      selectedLeptons->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
      selectedLeptons_copy->push_back(new lepton_candidate(elePt,(*gsf_eta)[l],(*gsf_phi)[l],(*gsf_charge)[l],l,1));
      if (data == "mc") sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
      if (data == "mc") sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,(*gsf_sc_eta)[l],(*gsf_pt)[l],"");
    }
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
      if(!(*mu_isTightMuon)[l]) continue;
      if((*mu_pfIsoDbCorrected04)[l] > 0.15) continue;
      selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l],(*mu_gt_charge)[l],l,10));
      selectedLeptons_copy->push_back(new lepton_candidate(muPtSFRochester * (*mu_gt_pt)[l],(*mu_gt_eta)[l],(*mu_gt_phi)[l],(*mu_gt_charge)[l],l,10));
      if (data == "mc" && year == "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
      if (data == "mc" && year == "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_eta)[l], (*mu_gt_pt)[l],"");
      if (data == "mc" && year != "2016") sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");
      if (data == "mc" && year != "2016") sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, (*mu_gt_pt)[l], abs((*mu_gt_eta)[l]),"");

    }
    sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
// trilepton selection
//cout<<ev_event<<"  "<<triggerPass<<"  "<<metFilterPass<<"  "<<selectedLeptons->size()<<endl;
    if(selectedLeptons->size()!=3 ||
      ((*selectedLeptons)[0]->pt_ <27) ||
      (abs((*selectedLeptons)[0]->charge_ + (*selectedLeptons)[1]->charge_ + (*selectedLeptons)[2]->charge_) != 1)) {
//      ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ != 11 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()<106 && ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()>76) ||
//      ((*selectedLeptons)[0]->p4_ + (*selectedLeptons)[1]->p4_).M()<20) {
      for (int l=0;l<selectedLeptons->size();l++){
        delete (*selectedLeptons)[l];
        delete (*selectedLeptons_copy)[l];
      }
      selectedLeptons->clear();
      selectedLeptons->shrink_to_fit();
      delete selectedLeptons;
      selectedLeptons_copy->clear();
      selectedLeptons_copy->shrink_to_fit();
      delete selectedLeptons_copy;
      continue;
    }

    if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ + (*selectedLeptons)[2]->lep_ == 3) ch = 0;//eee channel
    if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ + (*selectedLeptons)[2]->lep_ == 12) {
        ch = 1;//emul channel
        if ((*selectedLeptons)[0]->charge_ + (*selectedLeptons)[1]->charge_ + (*selectedLeptons)[2]->charge_ ==1){
            ch1 = 0; //eemu subchannel with ++- charge
        }
        else{
            ch1 = 1; //eemu subchannel with +-- charge
        }
    }
    if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ + (*selectedLeptons)[2]->lep_ == 21) {
        ch = 1;
        if ((*selectedLeptons)[0]->charge_ + (*selectedLeptons)[1]->charge_ + (*selectedLeptons)[2]->charge_ ==1){
            ch1 = 2; //emumu subchannel with ++- charge
        }
        else{
            ch1 = 3; //emumu subchannel with +-- charge
        }
      }
    if ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ + (*selectedLeptons)[2]->lep_ == 30) ch = 2; //mumumu channel

    if(ch ==0 && !triggerPassEE) continue;
    //if(ch ==1 && !triggerPassEMu) continue;
    // For now we use ee+emu+mumu trigger for emul channel
    if(ch ==1){
      sort(selectedLeptons_copy->begin(), selectedLeptons_copy->end(), CompareFlavourLep);
      if(ch1==0){
          if ((*selectedLeptons_copy)[2]->charge_<0){
              compete=true;// compete=true means that we have two pairs of leptons of Opposite sign&&Opposite flavour (OSOF), we need to make a decision based on kinematic reconstruction
          }
          else if ((*selectedLeptons_copy)[1]->charge_<0){
              compete=false;
              (*selectedLeptons_copy)[0]->setbalep();// This lepton is set as "Bachelor Lepton" meaning it comes from standard top
          }
          else {
              compete=false;
              (*selectedLeptons_copy)[1]->setbalep();
          }
      }
      if(ch1==1){
          if ((*selectedLeptons_copy)[2]->charge_>0){
              compete=true;
          }
          else if ((*selectedLeptons_copy)[1]->charge_>0){
              compete=false;
              (*selectedLeptons_copy)[0]->setbalep();
          }
          else {
              compete=false;
              (*selectedLeptons_copy)[1]->setbalep();
          }
      }
      if(ch1==2){
          if ((*selectedLeptons_copy)[0]->charge_<0){
              compete=true;
          }
          else if ((*selectedLeptons_copy)[1]->charge_<0){
              compete=false;
              (*selectedLeptons_copy)[2]->setbalep();
          }
          else {
              compete=false;
              (*selectedLeptons_copy)[1]->setbalep();
          }
      }
      if(ch1==3){
          if ((*selectedLeptons_copy)[0]->charge_>0){
              compete=true;
          }
          else if ((*selectedLeptons_copy)[1]->charge_>0){
              compete=false;
              (*selectedLeptons_copy)[2]->setbalep();
          }
          else {
              compete=false;
              (*selectedLeptons_copy)[1]->setbalep();
          }
      }
    }
    if(ch ==2 && !triggerPassMuMu) continue;
//jets
    selectedJets = new std::vector<jet_candidate*>();
    selectedJets_copy = new std::vector<jet_candidate*>();
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
        selectedJets_copy->push_back(new jet_candidate((*jet_Smeared_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l], year,(*jet_partonFlavour)[l]));
      }
      if(data == "data"){
        selectedJets->push_back(new jet_candidate((*jet_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l],year,0));
        selectedJets_copy->push_back(new jet_candidate((*jet_pt)[l],(*jet_eta)[l],(*jet_phi)[l],(*jet_energy)[l],(*jet_DeepCSV)[l],year,0));
      }
    }

    sort(selectedJets->begin(), selectedJets->end(), ComparePtJet);
    sort(selectedJets_copy->begin(), selectedJets_copy->end(), CompareBtagJet);// Orderd by b_tagging score
    nbjet=0;
    for (int l=0;l<selectedJets->size();l++){
      if((*selectedJets)[l]->btag_) nbjet++;
      if(data == "data") continue;
      if( abs((*selectedJets)[l]->flavor_) == 5){
        h2_BTaggingEff_Denom_b->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_ ) {
          h2_BTaggingEff_Num_b->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"");
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),""));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }  
      }
      if( abs((*selectedJets)[l]->flavor_) == 4){
        h2_BTaggingEff_Denom_c->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_) {
          h2_BTaggingEff_Num_c->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"");
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),""));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }
      if( abs((*selectedJets)[l]->flavor_) != 4 && abs((*selectedJets)[l]->flavor_) != 5){
        h2_BTaggingEff_Denom_udsg->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
        if( (*selectedJets)[l]->btag_) {
          h2_BTaggingEff_Num_udsg->Fill((*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_));
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"");
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),""));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }
    }
      
    //Two OSOFs ? let's compete
    if(compete&&selectedJets_copy->size()>0){
      (*selectedJets_copy)[0]->setbajet();
      if (ch1==0||ch1==1){
          t1=getTopmass((*selectedLeptons_copy)[0],(*selectedJets_copy)[0],MET_FinalCollection_Pt,MET_FinalCollection_phi);
          t2=getTopmass((*selectedLeptons_copy)[1],(*selectedJets_copy)[0],MET_FinalCollection_Pt,MET_FinalCollection_phi);
          if (t1<0&&t2<0) continue;
          if (abs(t1-mT)>abs(t2-mT)){// the one gives the better standard top mass wins
              Topmass=t2;
              (*selectedLeptons_copy)[1]->setbalep();
          }
          else{
              Topmass=t1;
              (*selectedLeptons_copy)[0]->setbalep();
          }
      }
      else{
          t1=getTopmass((*selectedLeptons_copy)[1],(*selectedJets_copy)[0],MET_FinalCollection_Pt,MET_FinalCollection_phi);
          t2=getTopmass((*selectedLeptons_copy)[2],(*selectedJets_copy)[0],MET_FinalCollection_Pt,MET_FinalCollection_phi);
          if (t1<0&&t2<0) continue;
          if (abs(t1-mT)>abs(t2-mT)){
              Topmass=t2;
              (*selectedLeptons_copy)[2]->setbalep();
          }
          else{
              Topmass=t1;
              (*selectedLeptons_copy)[1]->setbalep();
          }
      }
    }
      
    sort(selectedLeptons_copy->begin(), selectedLeptons_copy->end(), CompareBaLep);
    sort(selectedLeptons_copy->begin(), selectedLeptons_copy->begin()+2, CompareFlavourLep);// [e,mu,ba-lep]
    
//    if(ch==1){//debug
//    cout<<endl<<"ch1 = "<<ch1<<endl;
//    cout<<"compete = "<<compete<<endl;
//    cout<<"lep1 flavour = "<<(*selectedLeptons_copy)[0]->lep_<<" lep1 charge = "<<(*selectedLeptons_copy)[0]->charge_<<" ba = "<<(*selectedLeptons_copy)[0]->isbalep<<endl;
//    cout<<"lep2 flavour = "<<(*selectedLeptons_copy)[1]->lep_<<" lep1 charge = "<<(*selectedLeptons_copy)[1]->charge_<<" ba = "<<(*selectedLeptons_copy)[1]->isbalep<<endl;
//    cout<<"lep3 flavour = "<<(*selectedLeptons_copy)[2]->lep_<<" lep1 charge = "<<(*selectedLeptons_copy)[2]->charge_<<" ba = "<<(*selectedLeptons_copy)[2]->isbalep<<endl;
//      }

    if (data == "mc" && ch==0) sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");
    if (data == "mc" && ch==1) sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");
    if (data == "mc" && ch==2) sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");
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

    Hists[ch][0][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    Hists[ch][0][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    Hists[ch][0][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
    Hists[ch][0][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    Hists[ch][0][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    Hists[ch][0][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    Hists[ch][0][6]->Fill((*selectedLeptons)[2]->pt_,weight_lep);
    Hists[ch][0][7]->Fill((*selectedLeptons)[2]->eta_,weight_lep);
    Hists[ch][0][8]->Fill((*selectedLeptons)[2]->phi_,weight_lep);
    Hists[ch][0][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lep);
    Hists[ch][0][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lep);
    Hists[ch][0][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lep);
    Hists[ch][0][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lep);
    Hists[ch][0][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lep);
    Hists[ch][0][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lep);
    Hists[ch][0][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lep);
    Hists[ch][0][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lep);
    Hists[ch][0][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lep);
    Hists[ch][0][18]->Fill(Topmass,weight_lep);
    Hists[ch][0][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lep);
    Hists[ch][0][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lep);
    Hists[ch][0][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lep);
    Hists[ch][0][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lep);
    if(selectedJets->size()>0) Hists[ch][0][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][0][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][0][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lep);
    Hists[ch][0][26]->Fill(selectedJets->size(),weight_lep);
    Hists[ch][0][27]->Fill(nbjet,weight_lepB);
    Hists[ch][0][28]->Fill(MET_FinalCollection_Pt,weight_lep);
    Hists[ch][0][29]->Fill(MET_FinalCollection_phi,weight_lep);
    Hists[ch][0][30]->Fill(pv_n,weight_lep);
    Hists[ch][0][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lep);

    if (ch==1&&!compete){
    TLorentzVector zcand;
    if(ch1==0||ch1==1){
      zcand=(*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[2]->p4_;
    }
    else{
      zcand=(*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_;
    }
    if (zcand.M()>76&&zcand.M()<106) continue;//Off Z

    }

    Hists[ch][1][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    Hists[ch][1][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    Hists[ch][1][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
    Hists[ch][1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    Hists[ch][1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    Hists[ch][1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    Hists[ch][1][6]->Fill((*selectedLeptons)[2]->pt_,weight_lep);
    Hists[ch][1][7]->Fill((*selectedLeptons)[2]->eta_,weight_lep);
    Hists[ch][1][8]->Fill((*selectedLeptons)[2]->phi_,weight_lep);
    Hists[ch][1][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lep);
    Hists[ch][1][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lep);
    Hists[ch][1][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lep);
    Hists[ch][1][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lep);
    Hists[ch][1][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lep);
    Hists[ch][1][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lep);
    Hists[ch][1][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lep);
    Hists[ch][1][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lep);
    Hists[ch][1][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lep);
    Hists[ch][1][18]->Fill(Topmass,weight_lep);
    Hists[ch][1][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lep);
    Hists[ch][1][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lep);
    Hists[ch][1][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lep);
    Hists[ch][1][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lep);
    if(selectedJets->size()>0) Hists[ch][1][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][1][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][1][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lep);
    Hists[ch][1][26]->Fill(selectedJets->size(),weight_lep);
    Hists[ch][1][27]->Fill(nbjet,weight_lepB);
    Hists[ch][1][28]->Fill(MET_FinalCollection_Pt,weight_lep);
    Hists[ch][1][29]->Fill(MET_FinalCollection_phi,weight_lep);
    Hists[ch][1][30]->Fill(pv_n,weight_lep);
    Hists[ch][1][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lep);

    if(nbjet==1){
    Hists[ch][2][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][2][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][2][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][2][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][2][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][2][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][2][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][2][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][2][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][2][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][2][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][2][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][2][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][2][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][2][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][2][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][2][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][2][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][2][18]->Fill(Topmass,weight_lepB);
    Hists[ch][2][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][2][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][2][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][2][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][2][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][2][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][2][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lepB);
    Hists[ch][2][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][2][27]->Fill(nbjet,weight_lepB);
    Hists[ch][2][28]->Fill(MET_FinalCollection_Pt,weight_lepB);
    Hists[ch][2][29]->Fill(MET_FinalCollection_phi,weight_lepB);
    Hists[ch][2][30]->Fill(pv_n,weight_lepB);
    Hists[ch][2][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    }
    if(nbjet>1){
    Hists[ch][3][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][3][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][3][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][3][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][3][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][3][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][3][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][3][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][3][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][3][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][3][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][3][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][3][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][3][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][3][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][3][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][3][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][3][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][3][18]->Fill(Topmass,weight_lepB);
    Hists[ch][3][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][3][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][3][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][3][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][3][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][3][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][3][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lepB);
    Hists[ch][3][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][3][27]->Fill(nbjet,weight_lepB);
    Hists[ch][3][28]->Fill(MET_FinalCollection_Pt,weight_lepB);
    Hists[ch][3][29]->Fill(MET_FinalCollection_phi,weight_lepB);
    Hists[ch][3][30]->Fill(pv_n,weight_lepB);
    Hists[ch][3][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    }
    if(MET_FinalCollection_Pt<20){
    Hists[ch][4][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    Hists[ch][4][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    Hists[ch][4][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
    Hists[ch][4][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    Hists[ch][4][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    Hists[ch][4][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    Hists[ch][4][6]->Fill((*selectedLeptons)[2]->pt_,weight_lep);
    Hists[ch][4][7]->Fill((*selectedLeptons)[2]->eta_,weight_lep);
    Hists[ch][4][8]->Fill((*selectedLeptons)[2]->phi_,weight_lep);
    Hists[ch][4][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lep);
    Hists[ch][4][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lep);
    Hists[ch][4][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lep);
    Hists[ch][4][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lep);
    Hists[ch][4][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lep);
    Hists[ch][4][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lep);
    Hists[ch][4][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lep);
    Hists[ch][4][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lep);
    Hists[ch][4][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lep);
    Hists[ch][4][18]->Fill(Topmass,weight_lep);
    Hists[ch][4][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lep);
    Hists[ch][4][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lep);
    Hists[ch][4][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lep);
    Hists[ch][4][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lep);
    if(selectedJets->size()>0) Hists[ch][4][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][4][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][4][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lep);
    Hists[ch][4][26]->Fill(selectedJets->size(),weight_lep);
    Hists[ch][4][27]->Fill(nbjet,weight_lepB);
    Hists[ch][4][28]->Fill(MET_FinalCollection_Pt,weight_lep);
    Hists[ch][4][29]->Fill(MET_FinalCollection_phi,weight_lep);
    Hists[ch][4][30]->Fill(pv_n,weight_lep);
    Hists[ch][4][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lep);
    }

    if(MET_FinalCollection_Pt>20){
    Hists[ch][5][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    Hists[ch][5][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    Hists[ch][5][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
    Hists[ch][5][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    Hists[ch][5][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    Hists[ch][5][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    Hists[ch][5][6]->Fill((*selectedLeptons)[2]->pt_,weight_lep);
    Hists[ch][5][7]->Fill((*selectedLeptons)[2]->eta_,weight_lep);
    Hists[ch][5][8]->Fill((*selectedLeptons)[2]->phi_,weight_lep);
    Hists[ch][5][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lep);
    Hists[ch][5][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lep);
    Hists[ch][5][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lep);
    Hists[ch][5][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lep);
    Hists[ch][5][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lep);
    Hists[ch][5][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lep);
    Hists[ch][5][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lep);
    Hists[ch][5][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lep);
    Hists[ch][5][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lep);
    Hists[ch][5][18]->Fill(Topmass,weight_lep);
    Hists[ch][5][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lep);
    Hists[ch][5][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lep);
    Hists[ch][5][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lep);
    Hists[ch][5][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lep);
    if(selectedJets->size()>0) Hists[ch][5][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][5][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][5][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lep);
    Hists[ch][5][26]->Fill(selectedJets->size(),weight_lep);
    Hists[ch][5][27]->Fill(nbjet,weight_lepB);
    Hists[ch][5][28]->Fill(MET_FinalCollection_Pt,weight_lep);
    Hists[ch][5][29]->Fill(MET_FinalCollection_phi,weight_lep);
    Hists[ch][5][30]->Fill(pv_n,weight_lep);
    Hists[ch][5][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lep);
    }

    if(nbjet==1 && MET_FinalCollection_Pt<20 && selectedJets->size()>0){
    Hists[ch][6][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][6][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][6][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][6][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][6][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][6][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][6][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][6][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][6][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][6][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][6][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][6][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][6][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][6][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][6][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][6][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][6][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][6][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][6][18]->Fill(Topmass,weight_lepB);
    Hists[ch][6][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][6][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][6][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][6][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][6][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][6][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][6][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lepB);
    Hists[ch][6][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][6][27]->Fill(nbjet,weight_lepB);
    Hists[ch][6][28]->Fill(MET_FinalCollection_Pt,weight_lepB);
    Hists[ch][6][29]->Fill(MET_FinalCollection_phi,weight_lepB);
    Hists[ch][6][30]->Fill(pv_n,weight_lepB);
    Hists[ch][6][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    }

    if(nbjet==1 && MET_FinalCollection_Pt<20 && selectedJets->size()>1){
    Hists[ch][7][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][7][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][7][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][7][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][7][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][7][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][7][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][7][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][7][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][7][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][7][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][7][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][7][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][7][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][7][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][7][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][7][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][7][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][7][18]->Fill(Topmass,weight_lepB);
    Hists[ch][7][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][7][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][7][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][7][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][7][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][7][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][7][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lepB);
    Hists[ch][7][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][7][27]->Fill(nbjet,weight_lepB);
    Hists[ch][7][28]->Fill(MET_FinalCollection_Pt,weight_lepB);
    Hists[ch][7][29]->Fill(MET_FinalCollection_phi,weight_lepB);
    Hists[ch][7][30]->Fill(pv_n,weight_lepB);
    Hists[ch][7][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    }

    if(nbjet>1 && MET_FinalCollection_Pt<20 && selectedJets->size()>2){
    Hists[ch][8][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][8][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][8][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][8][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][8][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][8][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][8][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][8][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][8][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][8][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][8][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][8][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][8][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][8][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][8][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][8][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][8][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][8][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][8][18]->Fill(Topmass,weight_lepB);
    Hists[ch][8][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][8][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][8][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][8][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][8][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][8][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][8][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lepB);
    Hists[ch][8][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][8][27]->Fill(nbjet,weight_lepB);
    Hists[ch][8][28]->Fill(MET_FinalCollection_Pt,weight_lepB);
    Hists[ch][8][29]->Fill(MET_FinalCollection_phi,weight_lepB);
    Hists[ch][8][30]->Fill(pv_n,weight_lepB);
    Hists[ch][8][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    }

    if(nbjet==1 && MET_FinalCollection_Pt>20 && selectedJets->size()>0){
    Hists[ch][9][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][9][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][9][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][9][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][9][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][9][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][9][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][9][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][9][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][9][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][9][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][9][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][9][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][9][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][9][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][9][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][9][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][9][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][9][18]->Fill(Topmass,weight_lepB);
    Hists[ch][9][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][9][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][9][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][9][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][9][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][9][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][9][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lepB);
    Hists[ch][9][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][9][27]->Fill(nbjet,weight_lepB);
    Hists[ch][9][28]->Fill(MET_FinalCollection_Pt,weight_lepB);
    Hists[ch][9][29]->Fill(MET_FinalCollection_phi,weight_lepB);
    Hists[ch][9][30]->Fill(pv_n,weight_lepB);
    Hists[ch][9][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    }
    if(nbjet==1 && MET_FinalCollection_Pt>20 && selectedJets->size()>1){
    Hists[ch][10][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][10][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][10][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][10][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][10][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][10][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][10][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][10][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][10][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][10][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][10][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][10][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][10][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][10][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][10][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][10][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][10][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][10][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][10][18]->Fill(Topmass,weight_lepB);
    Hists[ch][10][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][10][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][10][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][10][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][10][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][10][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][10][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lepB);
    Hists[ch][10][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][10][27]->Fill(nbjet,weight_lepB);
    Hists[ch][10][28]->Fill(MET_FinalCollection_Pt,weight_lepB);
    Hists[ch][10][29]->Fill(MET_FinalCollection_phi,weight_lepB);
    Hists[ch][10][30]->Fill(pv_n,weight_lepB);
    Hists[ch][10][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    }
    if(nbjet>1 && MET_FinalCollection_Pt>20 && selectedJets->size()>2){
    Hists[ch][11][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][11][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][11][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][11][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][11][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][11][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][11][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][11][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][11][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][11][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][11][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][11][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][11][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][11][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][11][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][11][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][11][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][11][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][11][18]->Fill(Topmass,weight_lepB);
    Hists[ch][11][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][11][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][11][21]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][11][22]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][11][23]->Fill((*selectedJets_copy)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][11][24]->Fill((*selectedJets_copy)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][11][25]->Fill((*selectedJets_copy)[0]->phi_,weight_lepB);
    Hists[ch][11][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][11][27]->Fill(nbjet,weight_lepB);
    Hists[ch][11][28]->Fill(MET_FinalCollection_Pt,weight_lepB);
    Hists[ch][11][29]->Fill(MET_FinalCollection_phi,weight_lepB);
    Hists[ch][11][30]->Fill(pv_n,weight_lepB);
    Hists[ch][11][31]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    }

    for (int l=0;l<selectedLeptons->size();l++){
      delete (*selectedLeptons)[l];
      delete (*selectedLeptons_copy)[l];
    }
    for (int l=0;l<selectedJets->size();l++){
      delete (*selectedJets)[l];
      delete (*selectedJets_copy)[l];
    }
    selectedLeptons->clear();
    selectedLeptons->shrink_to_fit();
    delete selectedLeptons;
    selectedLeptons_copy->clear();
    selectedLeptons_copy->shrink_to_fit();
    delete selectedLeptons_copy;
    selectedJets->clear();
    selectedJets->shrink_to_fit();
    delete selectedJets;
    selectedJets_copy->clear();
    selectedJets_copy->shrink_to_fit();
    delete selectedJets_copy;

    nAccept++;
  } //end of event loop
  cout<<endl<<"from "<<ntr<<" events, "<<nAccept<<" events are accepted"<<endl;

  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        Hists[i][k][l]  ->Write("",TObject::kOverwrite);
      }
    }
  }
  for (int i=0;i<channels.size();++i){
    for (int k=0;k<regions.size();++k){
      for (int l=0;l<vars.size();++l){
        delete Hists[i][k][l];
      }
    }
  }

   h2_BTaggingEff_Denom_b   ->Write("",TObject::kOverwrite);
   h2_BTaggingEff_Denom_c   ->Write("",TObject::kOverwrite);
   h2_BTaggingEff_Denom_udsg->Write("",TObject::kOverwrite);
   h2_BTaggingEff_Num_b     ->Write("",TObject::kOverwrite);
   h2_BTaggingEff_Num_c     ->Write("",TObject::kOverwrite);
   h2_BTaggingEff_Num_udsg  ->Write("",TObject::kOverwrite);

   delete h2_BTaggingEff_Denom_b;
   delete h2_BTaggingEff_Denom_c;
   delete h2_BTaggingEff_Denom_udsg;
   delete h2_BTaggingEff_Num_b;
   delete h2_BTaggingEff_Num_c;
   delete h2_BTaggingEff_Num_udsg;
  file_out.Close() ;
}
