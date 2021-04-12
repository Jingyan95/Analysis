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
// test
// test2
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

bool verbose=false;

bool ComparePtLep(lepton_candidate *a, lepton_candidate *b) { return a->pt_ > b->pt_; }
bool CompareFlavourLep(lepton_candidate *a, lepton_candidate *b) { return a->lep_ < b->lep_; }
bool CompareBaLep(lepton_candidate *a, lepton_candidate *b) { return a->isbalep < b->isbalep; }
bool CompareChargeLep(lepton_candidate *a, lepton_candidate *b) { return a->charge_ < b->charge_; }
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

float getLFVTopmass(lepton_candidate *a,lepton_candidate *b,std::vector<jet_candidate*> *selectedJets){

  struct LFV_top_quark{
    int ij;
    float mass;
    float masserr;
  }temp,min_tq;
  min_tq.masserr=std::numeric_limits<float>::infinity();

  TLorentzVector lv_elmujet;

  for(int ij = 1;ij < (int)selectedJets->size();ij++){
    lv_elmujet = a->p4_ + b->p4_ + (*selectedJets)[ij]->p4_;
    temp.ij=ij;
    temp.mass=lv_elmujet.M();
    temp.masserr=fabs(lv_elmujet.M()-172);
    if(temp.masserr<min_tq.masserr){
      min_tq=temp;
    }
  }
  if(min_tq.masserr==std::numeric_limits<float>::infinity()){
    return -1;
  }
  else{
    return min_tq.mass;
  }
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

  typedef vector<TH1F*> Dim1;
  typedef vector<Dim1> Dim2;
  typedef vector<Dim2> Dim3;
  typedef vector<Dim3> Dim4;
		
  std::vector<TString> regions{"lll","lllOnZ","lllOffZ","lllOffZB0","lllOffZB1", "lllOffZBgeq2", "lllOffZMetl20", "lllOffZMetg20", "lllOffZMetg20Jetgeq1", "lllOffZMetg20Jet1B1", "lllOffZMetg20Jet2B1","lllOffZMetg20Jetgeq1Bleq1","lllOnZMetg20Jetgeq1Bleq1"};
  std::vector<TString> channels{"eee", "emul", "mumumu"};	
  std::vector<TString> vars   {"lep1Pt","lep1Eta","lep1Phi","lep2Pt","lep2Eta","lep2Phi","lep3Pt","lep3Eta","lep3Phi",	
"LFVePt","LFVeEta","LFVePhi","LFVmuPt","LFVmuEta","LFVmuPhi","balPt","balEta","balPhi","Topmass",	
"lllM","lllPt","lllHt","lllMt","jet1Pt","jet1Eta","jet1Phi","njet","nbjet","Met","MetPhi","nVtx",	
"llZM","llZPt","llZDr","llZDphi","LFVTopmass","llM","llPt","llDr","llDphi",	
"Ht","Ms","ZlDr","ZlDphi","JeDr","JmuDr","tM"};	

  std::vector<int>    nbins   {12      ,15       ,15       ,12      ,15       ,15       ,12      ,15       ,15       ,	
12      ,15       ,15       ,12       ,15        ,15        ,12     ,15      ,15      ,12       ,	
20    ,20     ,20     ,20     ,15      ,15       ,15       ,10    ,6      ,15   ,15      ,70    ,	
20    ,12     ,25     ,15       ,12          ,20   ,12    ,25    ,15      ,	
20  ,10  ,25    ,15      ,25    ,25     ,15};	

  std::vector<float> lowEdge  {30      ,-3       ,-4       ,20      ,-3       ,-4       ,20      ,-3       ,-4       ,	
20      ,-3       ,-4       ,20       ,-3        ,-4        ,20     ,-3      ,-4      ,100      ,	
0     ,0      ,0      ,0      ,30      ,-3       ,-4       ,0     ,0      ,0    ,-4      ,0     ,	
50    ,0      ,0      ,0        ,100         ,50   ,0     ,0     ,0       ,	
0   ,0   ,0     ,0       ,0     ,0      ,0};	

  std::vector<float> highEdge {300     ,3        ,4        ,200     ,3        ,4        ,150     ,3        ,4        ,	
240     ,3        ,4        ,240      ,3         ,4         ,120    ,3       ,4       ,340      ,	
500   ,500    ,500    ,500    ,300     ,3        ,4        ,10    ,6      ,210  ,4       ,70    ,	
130   ,240    ,7      ,4        ,340         ,130  ,240   ,7     ,4       ,	
600 ,100 ,7     ,4       ,7     ,7      ,150};
	
  Dim3 Hists(channels.size(),Dim2(regions.size(),Dim1(vars.size())));  
  std::stringstream name;
  TH1F *h_test;
  for (int i=0;i<(int)channels.size();++i){
    for (int k=0;k<(int)regions.size();++k){
      for (int l=0;l<(int)vars.size();++l){
        name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l];
        h_test = new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]);
        h_test->StatOverflows(kTRUE);
        h_test->Sumw2(kTRUE);
        Hists[i][k][l] = h_test;
        name.str("");
      }
    }
  }
    
  std::vector<TString> sys{"eleRecoSf", "eleIDSf", "muIdSf", "muIsoSf", "bcTagSF", "udsgTagSF","pu", "prefiring"};
  Dim4 HistsSysUp(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sys.size()))));
  Dim4 HistsSysDown(channels.size(),Dim3(regions.size(),Dim2(vars.size(),Dim1(sys.size()))));

  for (int i=0;i<(int)channels.size();++i){
    for (int k=0;k<(int)regions.size();++k){
      for (int l=0;l<(int)vars.size();++l){
        for (int n=0;n<(int)sys.size();++n){
              name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_"<<sys[n]<<"_Up";
              h_test = new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]);
              h_test->StatOverflows(kTRUE);
              h_test->Sumw2(kTRUE);
              HistsSysUp[i][k][l][n] = h_test;
              name.str("");
              name<<channels[i]<<"_"<<regions[k]<<"_"<<vars[l]<<"_"<<sys[n]<<"_Down";
              h_test = new TH1F((name.str()).c_str(),(name.str()).c_str(),nbins[l],lowEdge[l],highEdge[l]);
              h_test->StatOverflows(kTRUE);
              h_test->Sumw2(kTRUE);
              HistsSysDown[i][k][l][n] = h_test;
              name.str("");
        }
      }
    }
  }
    
    //Fake Factor Implementation
    Double_t ptBinsFF[8] = {20., 23., 27., 32., 38., 45., 65., 100.};

    std::vector<TString> channelsFF{"e","mu","emu"};
    std::vector<TString> etaregs{"barrel","transition","endcap"};
    std::vector<TString> varsFF   {"lep1Pt","lep1Eta","jet1Pt","jet1Eta","njet","nbjet","Met","nVtx"};
    std::vector<int>    nbinsFF   {12      ,15       ,15      ,15       ,10    ,6      ,15   ,70};
    std::vector<float> lowEdgeFF  {30      ,-3       ,30      ,-3       ,0     ,0      ,0    ,0};
    std::vector<float> highEdgeFF {300     ,3        ,300     ,3        ,10    ,6      ,210  ,70};

    typedef vector<TH2F*> Dim21;
    typedef vector<Dim21> Dim22;
    typedef vector<Dim22> Dim23;
    typedef vector<Dim23> Dim24;
    typedef vector<TH3F*> Dim31;
    typedef vector<Dim31> Dim32;
    typedef vector<Dim32> Dim33;
    typedef vector<Dim33> Dim34;
    typedef vector<Dim34> Dim35;

    Dim24 HistsAR1(channelsFF.size(),Dim23(2,Dim22(etaregs.size(),Dim21(varsFF.size()))));
    Dim35 HistsAR2(channelsFF.size(),Dim34(2,Dim33(etaregs.size(),Dim32(etaregs.size(),Dim31(varsFF.size())))));

    TH2F *h2_test;
    TH3F *h3_test;
    for (int i=0;i(int)<channelsFF.size();++i){
        for (int k=0;k<2;++k){
            for (int j=0;j<(int)etaregs.size();++j){
                for (int l=0;l<(int)varsFF.size();++l){
                    //Fill Bin Array
                    Double_t BinArray[nbinsFF[l]+1];
                    for (int n=0;n<=(int)nbinsFF[l];n++){
                        BinArray[n]=lowEdgeFF[l]+n*(highEdgeFF[l]-lowEdgeFF[l])/nbinsFF[l];
                    }
                    //Fill 2D/3D Histograms
                    name<<"AR1_2D_"<<channelsFF[i]<<"_"<<regions[regions.size()+k-2]<<"_"<<etaregs[j]<<"_"<<varsFF[l];
                    h2_test = new TH2F((name.str()).c_str(),(name.str()).c_str(),nbinsFF[l],BinArray,7,ptBinsFF);
                    h2_test->StatOverflows(kTRUE);
                    h2_test->Sumw2(kTRUE);
                    HistsAR1[i][k][j][l] = h2_test;
                    name.str("");
                    for (int m=0;m<(int)etaregs.size();++m){
                        name<<"AR2_3D_"<<channelsFF[i]<<"_"<<regions[regions.size()+k-2]<<"_"<<etaregs[j]<<"_"<<etaregs[m]<<"_"<<varsFF[l];
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
  //char* cmschar(getenv("CMSSW_BASE"));
  RoccoR  rc;
  if(year == "2016")    rochesterFile = CMSSW_BASE + std::string("/src/data/TopLFV/input/RoccoR2016.txt");
  if(year == "2017")    rochesterFile = CMSSW_BASE + std::string("/src/data/TopLFV/input/RoccoR2017.txt");
  if(year == "2018")    rochesterFile = CMSSW_BASE + std::string("/src/data/TopLFV/input/RoccoR2018.txt");
  rc.init(rochesterFile);

  if(data == "mc"){
    //std::string btf = CMSSW_BASE  + std::string("/src/data/TopLFV/input/btagEff.root");
    //const char* btagrootFile;
    //btagrootFile =  (const char*)btf ;
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

    if(year == "2016")    btagFile =  CMSSW_BASE + std::string("/src/data/TopLFV/input/DeepCSV_2016LegacySF_V1.csv");
    if(year == "2017")    btagFile =  CMSSW_BASE + std::string("/src/data/TopLFV/input/DeepCSV_94XSF_V4_B_F.csv");
    if(year == "2018")    btagFile =  CMSSW_BASE + std::string("/src/data/TopLFV/input/DeepCSV_102XSF_V1.csv");

    BTagCalibration calib("DeepCSV",btagFile);
    reader.load(calib,BTagEntry::FLAV_B,"comb"); 
    reader.load(calib,BTagEntry::FLAV_C,"comb");
    reader.load(calib,BTagEntry::FLAV_UDSG,"comb");

    if(year == "2016"){
      TFile *f_Ele_Reco_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root")).c_str() );
      sf_Ele_Reco_H = *(TH2F*)f_Ele_Reco_Map->Get("EGamma_SF2D");

      TFile *f_Ele_ID_Map = new TFile(( CMSSW_BASE  + std::string("/src/data/TopLFV/input/2016LegacyReReco_ElectronTight_Fall17V2.root")).c_str() ) ;
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

  std::vector<lepton_candidate*> *selectedLeptons;
  std::vector<lepton_candidate*> *selectedLeptons_copy;
  std::vector<jet_candidate*> *selectedJets;
  std::vector<jet_candidate*> *selectedJets_copy;
  TLorentzVector wp, wm, b, ab;
  std::vector<float> nominalWeights;
  nominalWeights.assign(sys.size(), 1);
  std::vector<float> sysUpWeights;
  sysUpWeights.assign(sys.size(), 1);
  std::vector<float> sysDownWeights;
  sysDownWeights.assign(sys.size(), 1);
  bool triggerPassEE;
  bool triggerPassEMu;
  bool triggerPassMuMu;
  bool triggerPassSingleE;
  bool triggerPassSingleMu;
  bool metFilterPass;
  bool ifTopPt=false;
  bool compete=false;
  bool OnZ=false;//Opposite Sign&&Same Flavor (OSSF) pair present in event
  bool ARveto;//For application, we only study emul channel
  int ch;
  int ch1;
  int nAccept=0;
  int nbjet;
  int anti_index;//This is useful only when there are two tight leptons and one anti-selected lepton.
  int anti_flavor;//
  int etabin1;//eta region of the first fakeable lepton
  int etabin2;//eta region of the second fakeable lepton
  int nTight; //# of tight leptons
  float sf_Ele_Reco;
  float sf_Ele_ID;
  float sf_Mu_ID;
  float sf_Mu_ISO;
  float sf_Trigger;
  float weight_PU;
  float weight_Lumi;
  float weight_lep;
  float weight_lepB;
  float weight_lepC;
  float weight_prefiring;
  float weight_topPt;
  float elePt;
  float eleEta;
  float t1,t2,t3,t4;
  float Topmass=0;
  float LFVTopmass=0;
  float Zmass;
  float Zpt;//Z candidate pT
  float ZDr;//dr between two final state leptons of Z candidate
  float ZDphi;//dphi between two final state leptons of Z candidate
  float ZlDr;//dr between Z cancidate and the 3rd lepton
  float ZlDphi;//dphi between Z cancidate and the 3rd lepton
  float JeDr;//dr between LFV electron and light jet
  float JmuDr;//dr between LFV muon and light jet
  float mT=173.07;
  float mZ=91.2;
  float Ht;//Scalar sum of pT of all objects
  float Ms;//Scalar sum of mass of all objects
  float tM;//Transverse mass of the W candidate
  double muPtSFRochester;
  double P_bjet_data;
  double P_bjet_mc;

  if (fname.Contains("TTTo2L2Nu")) ifTopPt=true;

  if (fChain == 0) return;
  Long64_t nbytes = 0, nb = 0;
  Long64_t ntr = eList->GetN();
  // loop over number of entries
  // Long64_t
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
    ch1 =10;
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
    weight_lepC =1;
    weight_prefiring =1;
    weight_topPt =1;
    P_bjet_data =1;
    P_bjet_mc =1;
    for (int n=0;n<(int)sys.size();++n){
      nominalWeights[n] =1;
      sysUpWeights[n] =1;
      sysDownWeights[n] =1;
    }
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
/* UNCOMMENT and switch to new branch names

    //trigger
    ////MC
    if(data == "mc" && year == "2016"){
      if(trig_HLT_Ele23_Ele12_CaloId
      L_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele27_WPTight_Gsf_accept ) triggerPassEE =true;
      if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele27_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept) triggerPassEMu =true;
      if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept || trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept) triggerPassMuMu =true;
    }
*/
    if(data == "mc" && year == "2017"){
      if(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele35_WPTight_Gsf) triggerPassEE =true;
      if(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Ele35_WPTight_Gsf || HLT_IsoMu27 ) triggerPassEMu =true;
      if(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_IsoMu27) triggerPassMuMu =true;
      if(HLT_Ele35_WPTight_Gsf) triggerPassSingleE =true;
      if(HLT_IsoMu24) triggerPassSingleMu =true;
    }
/*
    if(data == "mc" && year == "2018"){
      if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEE =true;
      if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele32_WPTight_Gsf_accept || trig_HLT_IsoMu24_accept) triggerPassEMu =true;
      if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept || trig_HLT_IsoMu24_accept) triggerPassMuMu =true;
    } 
*/
    ////DATA
    if(data == "data"){

      /*
      if(year == "2016"){
        if(run == "H"){
          if(dataset=="MuonEG"){
            if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_DZ_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleElectron"){
            if(!(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEE =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleMuon"){
            if(!(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassMuMu =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_DZ_accept || trig_HLT_Ele27_WPTight_Gsf_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassEMu =true;
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
            if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleElectron"){
            if(!(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEE =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele27_WPTight_Gsf_accept) triggerPassEMu =true;
          }
          if(dataset=="SingleMuon"){
            if(!(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassMuMu =true;
            if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele27_WPTight_Gsf_accept) && (trig_HLT_IsoMu24_accept || trig_HLT_IsoTkMu24_accept)) triggerPassEMu =true;
          }
          if(dataset=="DoubleEG"){
            if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEE =true;
          }
          if(dataset=="DoubleMu"){
            if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_accept || trig_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_accept) triggerPassMuMu =true;
          }
        }
      }
      */
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
/*
      if(year == "2018"){
        if(dataset=="MuonEG"){
          if(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) triggerPassEMu =true;
        }
        if(dataset=="EGamma"){
          if(trig_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_accept || trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEE =true;
          if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept) && trig_HLT_Ele32_WPTight_Gsf_accept) triggerPassEMu =true;
	}
	if(dataset=="SingleMuon"){
	  if(!trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept && trig_HLT_IsoMu24_accept) triggerPassMuMu =true;
	  if(!(trig_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_accept || trig_HLT_Ele32_WPTight_Gsf_accept) && trig_HLT_IsoMu24_accept) triggerPassEMu =true;
	}
	if(dataset=="DoubleMu"){
	  if(trig_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_accept) triggerPassMuMu =true;
	}
      }
*/
    }

    if (verbose ) {
    cout << ".............................................................................................." << endl;
    cout << "event " << kentry << endl;
    }
    if(!(triggerPassEE || triggerPassEMu || triggerPassMuMu)) continue;
    if(!metFilterPass) continue;
    if (verbose ) {
        cout<< "event passed Trigger & MET filters! "<<endl;
    };
      
// lepton selection
  selectedLeptons = new std::vector<lepton_candidate*>();//typlical ordered by pT
  selectedLeptons_copy = new std::vector<lepton_candidate*>();// ordered by [e, mu , bachelor lepton ]
  Ht=0;
  Ms=0;
  anti_index=0;
  anti_flavor=0;
  etabin1=0;
  etabin2=0;
  nTight=0;
  ARveto=false;

  if (verbose ){
      cout << "There are  " << nElectron << " Electrons and " << nMuon  << " Muons as well as  " <<  nJet << " Jets "<< endl ;
      if (nElectron>0){
      cout <<endl<< "Electron loop begins"  << endl;
      }
  }

// electron
    for (UInt_t l=0;l< nElectron ;l++){
      // nano electron pt is already corrected
      //elePt = (*gsf_ecalTrkEnergyPostCorr)[l]*sin(2.*atan(exp(-1.*(*gsf_eta)[l]))) ;

      elePt = Electron_pt[l]  ;   
      eleEta = Electron_eta[l] + Electron_deltaEtaSC[l]; 

      if (verbose ){
      cout << "loop over Electron number  " << l << " has pt  " << elePt << " and SC eta " <<  eleEta<< endl;
      cout << "mass  " << Electron_mass[l] << " phi   " << Electron_phi[l] << endl  ;
       }
        
      if(elePt <20 || abs(eleEta) > 2.4 || (abs(eleEta)> 1.4442 && (abs(eleEta)< 1.566))) continue;
      if((int) Electron_cutBased[l] < 1) continue; //  0 = veto for cut based ID
      selectedLeptons->push_back(new lepton_candidate(elePt,eleEta,Electron_phi[l],Electron_charge[l],l,1));
      selectedLeptons_copy->push_back(new lepton_candidate(elePt,eleEta,Electron_phi[l],Electron_charge[l],l,1));
      Ht+=((*selectedLeptons)[selectedLeptons->size()-1]->p4_).Pt();
      Ms+=((*selectedLeptons)[selectedLeptons->size()-1]->p4_).M();
      if (Electron_mvaTOP[l]>0.9){
        (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
        nTight++;
      }
        
      if (verbose ){
       cout << "Electron number  " << l << " is selected !" << endl ;
      }

      if (data == "mc"){
          sf_Ele_Reco = sf_Ele_Reco * scale_factor(&sf_Ele_Reco_H ,eleEta,Electron_pt[l],"");
          nominalWeights[0] = nominalWeights[0] * scale_factor(&sf_Ele_Reco_H ,eleEta,Electron_pt[l],"");
          sysUpWeights[0] = sysUpWeights[0] * scale_factor(&sf_Ele_Reco_H ,eleEta,Electron_pt[l],"up");
          sysDownWeights[0] = sysDownWeights[0] * scale_factor(&sf_Ele_Reco_H ,eleEta,Electron_pt[l],"down");
            
          sf_Ele_ID = sf_Ele_ID * scale_factor(&sf_Ele_ID_H ,eleEta,Electron_pt[l],"");
          nominalWeights[1] = nominalWeights[1] * scale_factor(&sf_Ele_ID_H ,eleEta,Electron_pt[l],"");
          sysUpWeights[1] = sysUpWeights[1] * scale_factor(&sf_Ele_ID_H ,eleEta,Electron_pt[l],"up");
          sysDownWeights[1] = sysDownWeights[1] * scale_factor(&sf_Ele_ID_H ,eleEta,Electron_pt[l],"down");
      }
    }

    if (verbose && nMuon>0){
          cout <<endl<< "Muon loop begins"  << endl;
    }
      
// Muon
    int genMuIdx =0;
    for (UInt_t l=0;l< nMuon ;l++){
        
      if(data == "data"){
        muPtSFRochester = rc.kScaleDT(Muon_charge[l], Muon_pt[l],Muon_eta[l],Muon_phi[l], 0, 0);
      }
      if (data == "mc"){
	     genMuIdx = Muon_genPartIdx[l];
         /// ??? I think mu_mc_index = Muon_genPartIdx
        
         if (genMuIdx!=-1 && abs(GenPart_pdgId[genMuIdx]) == 13) muPtSFRochester = rc.kSpreadMC(Muon_charge[l], Muon_pt[l],Muon_eta[l],Muon_phi[l], GenPart_pt[genMuIdx],0, 0);

         //if (Muon_nTrackerLayers[l] > 30)  cout << " !!!!!!!!!!!!!!!!!! set muPtSFRochester = 1 ,  Muon has charge " << Muon_charge[l] << " , phi=  " << Muon_phi[l] << " , Muon_nTrackerLayers[l] " <<  Muon_nTrackerLayers[l] << endl;   

         if (genMuIdx<0 &&  Muon_nTrackerLayers[l] < 30 ) {
	   //cout << "no gen level muons =(   " << endl;
	   //cout << " Muon has charge " << Muon_charge[l] << " , phi=  " << Muon_phi[l] << " , Muon_nTrackerLayers[l] " <<  Muon_nTrackerLayers[l] << endl;   

          muPtSFRochester = rc.kSmearMC(Muon_charge[l], Muon_pt[l] , Muon_eta[l] , Muon_phi[l], Muon_nTrackerLayers[l] , gRandom->Rndm(),0, 0);
         }
      }
      
      if (verbose ){
         cout << "loop over Muon  number  " << l << " has pt  " << muPtSFRochester * Muon_pt[l] << " and eta " <<  Muon_eta[l] << endl;
         //cout << "Muon_genPartIdx[l] = genMuIdx = " << genMuIdx  << "  abs(GenPart_pdgId[Muon_genPartIdx[l]])   " <<  abs(GenPart_pdgId[Muon_genPartIdx[l]]) << endl;
         cout << "Muon_pfRelIso04_all[l] = " << Muon_pfRelIso04_all[l] <<endl;
      }
        
      if( (muPtSFRochester * Muon_pt[l] <20) || (abs(Muon_eta[l]) > 2.4) ) continue;

      ///  Muon_mvaId == 2 is mu_MvaMedium
      /// https://github.com/Fedespring/cmssw/blob/3f7b3c37caeaaf058bb1c7461b9c3c91a0672f68/PhysicsTools/NanoAOD/python/muons_cff.py#L138
      //if ((!(*mu_MvaMedium)[l]) || (!(*mu_CutBasedIdMedium)[l])) continue;
      if (!(Muon_mediumId[l]) || Muon_pfRelIso04_all[l] > 0.15) continue;
    
      selectedLeptons->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
      selectedLeptons_copy->push_back(new lepton_candidate(muPtSFRochester * Muon_pt[l],Muon_eta[l],Muon_phi[l],Muon_charge[l],l,10));
      Ht+=((*selectedLeptons)[selectedLeptons->size()-1]->p4_).Pt();
      Ms+=((*selectedLeptons)[selectedLeptons->size()-1]->p4_).M();
      if (Muon_mvaTOP[l]>0.7){
         (*selectedLeptons)[selectedLeptons->size()-1]->setTag();
         nTight++;
      }
      
      if (verbose ){
	     cout << "Muon number  " << l << " is selected !" << endl ;
      }
        
      if (data == "mc" && year == "2016") {
	sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_eta[l], Muon_pt[l],"");
	nominalWeights[2] = nominalWeights[2] * scale_factor(&sf_Mu_ID_H, Muon_eta[l], Muon_pt[l],"");
	sysUpWeights[2] = sysUpWeights[2] * scale_factor(&sf_Mu_ID_H, Muon_eta[l], Muon_pt[l],"up");
	sysDownWeights[2] = sysDownWeights[2] * scale_factor(&sf_Mu_ID_H, Muon_eta[l], Muon_pt[l],"down");
            
	sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_eta[l], Muon_pt[l],"");
	nominalWeights[3] = nominalWeights[3] * scale_factor(&sf_Mu_ISO_H, Muon_eta[l], Muon_pt[l],"");
	sysUpWeights[3] = sysUpWeights[3] * scale_factor(&sf_Mu_ISO_H, Muon_eta[l], Muon_pt[l],"up");
	sysDownWeights[3] = sysDownWeights[3] * scale_factor(&sf_Mu_ISO_H, Muon_eta[l], Muon_pt[l],"down");
      }
      if (data == "mc" && year != "2016") {
	sf_Mu_ID = sf_Mu_ID * scale_factor(&sf_Mu_ID_H, Muon_pt[l], abs(Muon_eta[l]),"");
	nominalWeights[2] = nominalWeights[2] * scale_factor(&sf_Mu_ID_H, Muon_pt[l], abs(Muon_eta[l]),"");
	sysUpWeights[2] = sysUpWeights[2] * scale_factor(&sf_Mu_ID_H, Muon_pt[l], abs(Muon_eta[l]),"up");
	sysDownWeights[2] = sysDownWeights[2] * scale_factor(&sf_Mu_ID_H, Muon_pt[l], abs(Muon_eta[l]),"down");
            
	sf_Mu_ISO = sf_Mu_ISO * scale_factor(&sf_Mu_ISO_H, Muon_pt[l], abs(Muon_eta[l]),"");
	nominalWeights[3] = nominalWeights[3] * scale_factor(&sf_Mu_ISO_H, Muon_pt[l], abs(Muon_eta[l]),"");
	sysUpWeights[3] = sysUpWeights[3] * scale_factor(&sf_Mu_ISO_H, Muon_pt[l], abs(Muon_eta[l]),"up");
	sysDownWeights[3] = sysDownWeights[3] * scale_factor(&sf_Mu_ISO_H, Muon_pt[l], abs(Muon_eta[l]),"down");
      }
    }
    sort(selectedLeptons->begin(), selectedLeptons->end(), ComparePtLep);
    if (selectedLeptons->size()==3 && ((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ + (*selectedLeptons)[2]->lep_ == 3 ||
        (*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ + (*selectedLeptons)[2]->lep_ == 30) && nTight<3){
        ARveto=true;
    }
      
    if (verbose ){
          cout <<endl<< "event has X leptons passing full leption selection, X = " <<  selectedLeptons->size()  << endl;
    }
     
    // trilepton selection
    if(selectedLeptons->size()!=3 || ARveto ||
       ((*selectedLeptons)[0]->pt_ <30) || !((*selectedLeptons)[0]->isTag) ||
       (abs((*selectedLeptons)[0]->charge_ + (*selectedLeptons)[1]->charge_ + (*selectedLeptons)[2]->charge_) != 1)) {
      for (int l=0;l<(int)selectedLeptons->size();l++){
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

    compete=false;
    if (nTight==1){
       if ((*selectedLeptons)[0]->lep_ == 1 && !triggerPassSingleE) continue;
       if ((*selectedLeptons)[0]->lep_ == 10 && !triggerPassSingleMu) continue;
    }
    else{
       if((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 2 && !triggerPassEE) continue;
       if((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 11 && !triggerPassEMu) continue;
       if((*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 20 && !triggerPassMuMu) continue;
    }
    if(ch == 1){
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
    if(ch == 0||ch == 2) sort(selectedLeptons_copy->begin(), selectedLeptons_copy->end(), CompareChargeLep);
    //jets
    selectedJets = new std::vector<jet_candidate*>();
    selectedJets_copy = new std::vector<jet_candidate*>();
    bool jetlepfail;
      
    if (verbose && nJet>0){
       cout <<endl<< "Jet loop begins"  << endl;
    }
      
    for (UInt_t l=0;l< nJet;l++){
      if (verbose ){
          cout << "loop over Jet number  " << l << " has pt  " << Jet_pt[l] << " and eta " <<  Jet_eta[l] << endl;
          cout << "mass  " << Jet_mass[l] << " phi   " << Jet_phi[l] << endl  ;   
      }
      if(data == "mc" && ((Jet_pt)[l] <30 || abs((Jet_eta)[l]) > 2.4)) continue;
      /// ??? xwhat is jet_smeared_pt in nano? 
      //if(data == "mc" && ((Jet_pt)[l] <30 || abs((Jet_eta)[l]) > 2.4)) continue;


      if(data == "data" && ((Jet_pt)[l] <30 || abs((Jet_eta)[l]) > 2.4)) continue;

      //if(year == "2016" && !(Jet_isJetIDTightLepVeto_2016)[l]) continue;
      //if(year == "2017" && !(Jet_isJetIDLepVeto_2017)[l]) continue;
      //if(year == "2018" && !(Jet_isJetIDLepVeto_2018)[l]) continue;
      if( (int) (Jet_jetId)[l] < 6 ) continue;
      // jet ID of 6 is tight with lep veto
      // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD


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
        /// ??? to-do : smear the jet pts in MC
        ///selectedJets->push_back(new jet_candidate((Jet_Smeared_pt)[l],(Jet_eta)[l],(Jet_phi)[l],(Jet_energy)[l],(Jet_DeepCSV)[l], year,(Jet_partonFlavour)[l]));
        //selectedJets_copy->push_back(new jet_candidate((Jet_Smeared_pt)[l],(Jet_eta)[l],(Jet_phi)[l],(Jet_energy)[l],(Jet_DeepCSV)[l], year,(Jet_partonFlavour)[l]));
        // ??? jet_Smeared_pt = Jet_pt in nano ? ... is jet_DeepCSV = Jet_btagDeepB?
        selectedJets->push_back(new jet_candidate((Jet_pt)[l],(Jet_eta)[l],(Jet_phi)[l],JetEnergy,(Jet_btagDeepB)[l], year,(Jet_partonFlavour)[l]));
        selectedJets_copy->push_back(new jet_candidate((Jet_pt)[l],(Jet_eta)[l],(Jet_phi)[l],JetEnergy,(Jet_btagDeepB)[l], year,(Jet_partonFlavour)[l]));


      }
      if(data == "data"){
        selectedJets->push_back(new jet_candidate((Jet_pt)[l],(Jet_eta)[l],(Jet_phi)[l], JetEnergy ,(Jet_btagDeepB)[l],year,0));
        selectedJets_copy->push_back(new jet_candidate((Jet_pt)[l],(Jet_eta)[l],(Jet_phi)[l], JetEnergy ,(Jet_btagDeepB)[l],year,0));
      }
      Ht+=((*selectedJets)[selectedJets->size()-1]->p4_).Pt();
      Ms+=((*selectedJets)[selectedJets->size()-1]->p4_).M();
        
      if (verbose ){
           cout << "Jet number  " << l << " is selected !" << endl ;
      }
        
    }

    sort(selectedJets->begin(), selectedJets->end(), ComparePtJet);
    sort(selectedJets_copy->begin(), selectedJets_copy->end(), CompareBtagJet);// Orderd by b_tagging score

    nbjet=0;
    for (int l=0;l<(int)selectedJets->size();l++){
      if((*selectedJets)[l]->btag_) nbjet++;
      if(data == "data") continue;
      // B-tagging systematics
      if( abs((*selectedJets)[l]->flavor_) == 5){
        if( (*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"");
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          nominalWeights[4] = nominalWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[4] = sysUpWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("up", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[4] = sysDownWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("down", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
            
          nominalWeights[5] = nominalWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[5] = sysUpWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("up", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[5] = sysDownWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("down", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),""));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          nominalWeights[4] = nominalWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[4] = sysUpWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("up", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[4] = sysDownWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("down", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
            
          nominalWeights[5] = nominalWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[5] = sysUpWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[5] = sysDownWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_B,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }  
      }
      if( abs((*selectedJets)[l]->flavor_) == 4){
        if( (*selectedJets)[l]->btag_) {
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"");
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          nominalWeights[4] = nominalWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[4] = sysUpWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("up", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[4] = sysDownWeights[4] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("down", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
            
          nominalWeights[5] = nominalWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[5] = sysUpWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[5] = sysDownWeights[5] * scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),""));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_c_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          nominalWeights[4] = nominalWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[4] = sysUpWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("up", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[4] = sysDownWeights[4]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("down", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
            
          nominalWeights[5] = nominalWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[5] = sysUpWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[5] = sysDownWeights[5]* (1- (scale_factor(&btagEff_b_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_C,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }
      if( abs((*selectedJets)[l]->flavor_) != 4 && abs((*selectedJets)[l]->flavor_) != 5){
        if( (*selectedJets)[l]->btag_) {
          P_bjet_mc = P_bjet_mc * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"");
          P_bjet_data = P_bjet_data * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          nominalWeights[4] = nominalWeights[4]* scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[4] = sysUpWeights[4]* scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[4] = sysDownWeights[4]* scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
            
          nominalWeights[5] = nominalWeights[5] * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysUpWeights[5] = sysUpWeights[5] * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
          sysDownWeights[5] = sysDownWeights[5] * scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_);
        }
        if( !(*selectedJets)[l]->btag_ ) {
          P_bjet_mc = P_bjet_mc * (1 - scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),""));
          P_bjet_data = P_bjet_data * (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          nominalWeights[4] = nominalWeights[4]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[4] = sysUpWeights[4]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[4] = sysDownWeights[4]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
            
          nominalWeights[5] = nominalWeights[5]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysUpWeights[5] = sysUpWeights[5]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
          sysDownWeights[5] = sysDownWeights[5]* (1- (scale_factor(&btagEff_udsg_H, (*selectedJets)[l]->pt_, abs((*selectedJets)[l]->eta_),"") * reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG,  abs((*selectedJets)[l]->eta_), (*selectedJets)[l]->pt_)));
        }
      }
    }
      
    //Two OSOFs ? let's compete
    if(ch==1&&compete&&(int)selectedJets_copy->size()>0){
      (*selectedJets_copy)[0]->setbajet();
      if (ch1==0||ch1==1){
	t1=getTopmass((*selectedLeptons_copy)[0],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	t2=getTopmass((*selectedLeptons_copy)[1],(*selectedJets_copy)[0],(MET_pt),(MET_phi)); /// ??? which MET do we want?
	t3=getLFVTopmass((*selectedLeptons_copy)[0],(*selectedLeptons_copy)[2],selectedJets_copy);
	t4=getLFVTopmass((*selectedLeptons_copy)[1],(*selectedLeptons_copy)[2],selectedJets_copy);
	if (t1<0&&t2<0) continue;
	if (abs(t1-mT)>abs(t2-mT)){// the one gives the better standard top mass wins
	  Topmass=t2;
	  LFVTopmass=t3;
	  (*selectedLeptons_copy)[1]->setbalep();
	}
	else{
	  Topmass=t1;
	  LFVTopmass=t4;
	  (*selectedLeptons_copy)[0]->setbalep();
	}
      }
      else{
	t1=getTopmass((*selectedLeptons_copy)[1],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	t2=getTopmass((*selectedLeptons_copy)[2],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	t3=getLFVTopmass((*selectedLeptons_copy)[1],(*selectedLeptons_copy)[0],selectedJets_copy);
	t4=getLFVTopmass((*selectedLeptons_copy)[2],(*selectedLeptons_copy)[0],selectedJets_copy);
	if (t1<0&&t2<0) continue;
	if (abs(t1-mT)>abs(t2-mT)){
	  Topmass=t2;
	  LFVTopmass=t3;
	  (*selectedLeptons_copy)[2]->setbalep();
	}
	else{
	  Topmass=t1;
	  LFVTopmass=t4;
	  (*selectedLeptons_copy)[1]->setbalep();
	}
      }
    }
    
    if(ch==1&&!compete&&selectedJets_copy->size()>0){
      if (ch1==0){
	if ((*selectedLeptons_copy)[0]->charge_>0){
	  Topmass=getTopmass((*selectedLeptons_copy)[0],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	  LFVTopmass=getLFVTopmass((*selectedLeptons_copy)[1],(*selectedLeptons_copy)[2],selectedJets_copy);
	}
	else{
	  Topmass=getTopmass((*selectedLeptons_copy)[1],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	  LFVTopmass=getLFVTopmass((*selectedLeptons_copy)[0],(*selectedLeptons_copy)[2],selectedJets_copy);
	}
      }
      if (ch1==1){
	if ((*selectedLeptons_copy)[0]->charge_>0){
	  Topmass=getTopmass((*selectedLeptons_copy)[1],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	  LFVTopmass=getLFVTopmass((*selectedLeptons_copy)[0],(*selectedLeptons_copy)[2],selectedJets_copy);
	}
	else{
	  Topmass=getTopmass((*selectedLeptons_copy)[0],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	  LFVTopmass=getLFVTopmass((*selectedLeptons_copy)[1],(*selectedLeptons_copy)[2],selectedJets_copy);
	}
      }
      if (ch1==2){
	if ((*selectedLeptons_copy)[1]->charge_>0){
	  Topmass=getTopmass((*selectedLeptons_copy)[1],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	  LFVTopmass=getLFVTopmass((*selectedLeptons_copy)[0],(*selectedLeptons_copy)[2],selectedJets_copy);
	}
	else{
	  Topmass=getTopmass((*selectedLeptons_copy)[2],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	  LFVTopmass=getLFVTopmass((*selectedLeptons_copy)[0],(*selectedLeptons_copy)[1],selectedJets_copy);
	}
      }
      if (ch1==3){
	if ((*selectedLeptons_copy)[1]->charge_>0){
	  Topmass=getTopmass((*selectedLeptons_copy)[2],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	  LFVTopmass=getLFVTopmass((*selectedLeptons_copy)[0],(*selectedLeptons_copy)[1],selectedJets_copy);
	}
	else{
	  Topmass=getTopmass((*selectedLeptons_copy)[1],(*selectedJets_copy)[0],(MET_pt),(MET_phi));
	  LFVTopmass=getLFVTopmass((*selectedLeptons_copy)[0],(*selectedLeptons_copy)[2],selectedJets_copy);
	}
      }
    }
      
    OnZ=false;
    Zmass=0;
    Zpt=0;
    ZDr=9;//Good for ML?
    ZDphi=6;
    ZlDr=9;
    ZlDphi=6;
    JeDr=9;
    JmuDr=9;
    tM=0;
    //Reconstruct Z candidate
    if (ch==1&&!compete){
      if(ch1==0||ch1==1){
	    Zmass=((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[2]->p4_).M();
        Zpt=((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[1]->p4_).Pt();
        ZDphi=deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_);
        ZDr=deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_);
        ZlDr=deltaR(((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[1]->p4_).Eta(),((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[1]->p4_).Phi(),(*selectedLeptons_copy)[2]->eta_,(*selectedLeptons_copy)[2]->phi_);
        ZlDphi=deltaPhi(((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[1]->p4_).Phi(),(*selectedLeptons_copy)[2]->phi_);
        tM=sqrt(2*(MET_pt)*((*selectedLeptons_copy)[2]->p4_).Et()*(1-cos(deltaPhi((*selectedLeptons_copy)[2]->phi_,(MET_phi)))));
      }
      else{
	    Zmass=((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).M();
        Zpt=((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).Pt();
        ZDphi=deltaPhi((*selectedLeptons_copy)[1]->phi_,(*selectedLeptons_copy)[2]->phi_);
        ZDr=deltaR((*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_,(*selectedLeptons_copy)[2]->eta_,(*selectedLeptons_copy)[2]->phi_);
        ZlDr=deltaR(((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).Eta(),((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).Phi(),(*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_);
        ZlDphi=deltaPhi(((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).Phi(),(*selectedLeptons_copy)[0]->phi_);
        tM=sqrt(2*(MET_pt)*((*selectedLeptons_copy)[0]->p4_).Et()*(1-cos(deltaPhi((*selectedLeptons_copy)[0]->phi_,(MET_phi)))));
      }
    }
    if (ch==0||ch==2){
        Zmass=((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[2]->p4_).M();
        Zpt=((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[2]->p4_).Pt();
        ZDphi=deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[2]->phi_);
        ZDr=deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[2]->eta_,(*selectedLeptons_copy)[2]->phi_);
        ZlDr=deltaR(((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[2]->p4_).Eta(),((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[2]->p4_).Phi(),(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_);
        ZlDphi=deltaPhi(((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[2]->p4_).Phi(),(*selectedLeptons_copy)[1]->phi_);
        tM=sqrt(2*(MET_pt)*((*selectedLeptons_copy)[1]->p4_).Et()*(1-cos(deltaPhi((*selectedLeptons_copy)[1]->phi_,(MET_phi)))));
        if(((*selectedLeptons)[0]->charge_ + (*selectedLeptons)[1]->charge_ + (*selectedLeptons)[2]->charge_)>0){
          if (abs(((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[1]->p4_).M()-mZ)<abs(Zmass-mZ)){
              Zmass=((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[1]->p4_).M();
              Zpt=((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[1]->p4_).Pt();
              ZDphi=deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_);
              ZDr=deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_);
              ZlDr=deltaR(((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[1]->p4_).Eta(),((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[1]->p4_).Phi(),(*selectedLeptons_copy)[2]->eta_,(*selectedLeptons_copy)[2]->phi_);
              ZlDphi=deltaPhi(((*selectedLeptons_copy)[0]->p4_+(*selectedLeptons_copy)[1]->p4_).Phi(),(*selectedLeptons_copy)[2]->phi_);
              tM=sqrt(2*(MET_pt)*((*selectedLeptons_copy)[2]->p4_).Et()*(1-cos(deltaPhi((*selectedLeptons_copy)[2]->phi_,(MET_phi)))));
              }
        }
        else if (abs(((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).M()-mZ)<abs(Zmass-mZ)){
          Zmass=((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).M();
          Zpt=((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).Pt();
          ZDphi=deltaPhi((*selectedLeptons_copy)[1]->phi_,(*selectedLeptons_copy)[2]->phi_);
          ZDr=deltaR((*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_,(*selectedLeptons_copy)[2]->eta_,(*selectedLeptons_copy)[2]->phi_);
          ZlDr=deltaR(((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).Eta(),((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).Phi(),(*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_);
          ZlDphi=deltaPhi(((*selectedLeptons_copy)[1]->p4_+(*selectedLeptons_copy)[2]->p4_).Phi(),(*selectedLeptons_copy)[0]->phi_);
          tM=sqrt(2*(MET_pt)*((*selectedLeptons_copy)[0]->p4_).Et()*(1-cos(deltaPhi((*selectedLeptons_copy)[0]->phi_,(MET_phi)))));
        }
    }
      
    if (Zmass>76&&Zmass<106) {
      OnZ=true;
    }
      
    if(ch == 1){
      sort(selectedLeptons_copy->begin(), selectedLeptons_copy->end(), CompareBaLep);
      sort(selectedLeptons_copy->begin(), selectedLeptons_copy->begin()+2, CompareFlavourLep);// [e,mu,ba-lep]
    }
    else{
    sort(selectedLeptons_copy->begin(), selectedLeptons_copy->end(), ComparePtLep);
    }
          
    if (selectedJets->size()>=2){
        JeDr=deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedJets_copy)[1]->eta_,(*selectedJets_copy)[1]->phi_);
        JmuDr=deltaR((*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_,(*selectedJets_copy)[1]->eta_,(*selectedJets_copy)[1]->phi_);
    }

    if (data == "mc" && (*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 2) {
      sf_Trigger = scale_factor(&sf_triggeree_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");}
    if (data == "mc" && (*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 11) {
      sf_Trigger = scale_factor(&sf_triggeremu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");}
    if (data == "mc" && (*selectedLeptons)[0]->lep_ + (*selectedLeptons)[1]->lep_ == 20) {
      sf_Trigger = scale_factor(&sf_triggermumu_H, (*selectedLeptons)[0]->pt_, (*selectedLeptons)[1]->pt_,"");}

    //PU reweighting
    if (data == "mc" && year == "2016") { /// is  mc_trueNumInteractions = Pileup_nTrueInt
      weight_PU = wPU.PU_2016(Pileup_nTrueInt ,"nominal");
      nominalWeights[6] = wPU.PU_2016(Pileup_nTrueInt,"nominal");
      sysUpWeights[6] = wPU.PU_2016(Pileup_nTrueInt,"up");
      sysDownWeights[6] = wPU.PU_2016(Pileup_nTrueInt,"down");
    }
    if (data == "mc" && year == "2017") {
      weight_PU = wPU.PU_2017(Pileup_nTrueInt,"nominal");
      nominalWeights[6] = wPU.PU_2017(Pileup_nTrueInt,"nominal");
      sysUpWeights[6] = wPU.PU_2017(Pileup_nTrueInt,"up");
      sysDownWeights[6] = wPU.PU_2017(Pileup_nTrueInt,"down");
    }
    if (data == "mc" && year == "2018") {
      weight_PU = wPU.PU_2018(Pileup_nTrueInt,"nominal");
      nominalWeights[6] = wPU.PU_2018(Pileup_nTrueInt,"nominal");
      sysUpWeights[6] = wPU.PU_2018(Pileup_nTrueInt,"up");
      sysDownWeights[6] = wPU.PU_2018(Pileup_nTrueInt,"down");
    }
    //lumi xs weights
    if (data == "mc") weight_Lumi = (1000*xs*lumi)/Nevent;
    
    if (data == "mc" && (year == "2016" || year == "2017")) { 
      float weightPre;
      weightPre = 1;
      if ( L1PreFiringWeight_Nom >=0 && L1PreFiringWeight_Nom <= 1.0 ){
      	weightPre =  L1PreFiringWeight_Nom;
      } 
      weight_prefiring = weightPre;
      nominalWeights[7] = weightPre;
      sysUpWeights[7] = L1PreFiringWeight_Up;
      sysDownWeights[7] = L1PreFiringWeight_Dn;
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

    //  \\ fix this ??? get sign of the the gen W ?
    float mc_w_sign;
    mc_w_sign = 1.;

    if (data == "mc") weight_lep = sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_Lumi  * mc_w_sign * weight_prefiring * weight_topPt;
    if (data == "mc") weight_lepB = sf_Ele_Reco * sf_Ele_ID * sf_Mu_ID * sf_Mu_ISO * sf_Trigger * weight_PU * weight_Lumi * mc_w_sign *  weight_prefiring * weight_topPt * (P_bjet_data/P_bjet_mc);
    
    weight_lepC = weight_lep;
    if (ch==1&&compete) weight_lepC = weight_lepB;

    // ??? fix this? which MET do we want?
    float MET_pt0;
    MET_pt0 = (MET_pt);
    float MET_phi0;
    MET_phi0 = (MET_phi);
        
    if(nTight==3){

    Hists[ch][0][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    Hists[ch][0][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    Hists[ch][0][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
    Hists[ch][0][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    Hists[ch][0][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    Hists[ch][0][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    Hists[ch][0][6]->Fill((*selectedLeptons)[2]->pt_,weight_lep);
    Hists[ch][0][7]->Fill((*selectedLeptons)[2]->eta_,weight_lep);
    Hists[ch][0][8]->Fill((*selectedLeptons)[2]->phi_,weight_lep);
    Hists[ch][0][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepC);
    Hists[ch][0][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepC);
    Hists[ch][0][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepC);
    Hists[ch][0][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepC);
    Hists[ch][0][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepC);
    Hists[ch][0][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepC);
    Hists[ch][0][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepC);
    Hists[ch][0][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepC);
    Hists[ch][0][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepC);
    Hists[ch][0][18]->Fill(Topmass,weight_lepB);
    Hists[ch][0][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepC);
    Hists[ch][0][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepC);
    Hists[ch][0][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepC);
    Hists[ch][0][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepC);
    if(selectedJets->size()>0) Hists[ch][0][23]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][0][24]->Fill((*selectedJets)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][0][25]->Fill((*selectedJets)[0]->phi_,weight_lep);
    Hists[ch][0][26]->Fill(selectedJets->size(),weight_lep);
    Hists[ch][0][27]->Fill(nbjet,weight_lepB);
    Hists[ch][0][28]->Fill(MET_pt0,weight_lep);
    Hists[ch][0][29]->Fill(MET_phi0,weight_lep);
    Hists[ch][0][30]->Fill(Pileup_nTrueInt,weight_lep);
    Hists[ch][0][31]->Fill(Zmass,weight_lepC);
    Hists[ch][0][32]->Fill(Zpt,weight_lepC);
    Hists[ch][0][33]->Fill(ZDr,weight_lepC);
    Hists[ch][0][34]->Fill(ZDphi,weight_lepC);
    Hists[ch][0][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][0][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepC);
    Hists[ch][0][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepC);
    Hists[ch][0][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][0][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][0][40]->Fill(Ht,weight_lepC);
    Hists[ch][0][41]->Fill(Ms,weight_lepC);
    Hists[ch][0][42]->Fill(ZlDr,weight_lepC);
    Hists[ch][0][43]->Fill(ZlDphi,weight_lepC);
    Hists[ch][0][44]->Fill(JeDr,weight_lepB);
    Hists[ch][0][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][0][46]->Fill(tM,weight_lepB);

    if(OnZ){
    Hists[ch][1][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    Hists[ch][1][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    Hists[ch][1][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
    Hists[ch][1][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    Hists[ch][1][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    Hists[ch][1][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    Hists[ch][1][6]->Fill((*selectedLeptons)[2]->pt_,weight_lep);
    Hists[ch][1][7]->Fill((*selectedLeptons)[2]->eta_,weight_lep);
    Hists[ch][1][8]->Fill((*selectedLeptons)[2]->phi_,weight_lep);
    Hists[ch][1][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepC);
    Hists[ch][1][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepC);
    Hists[ch][1][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepC);
    Hists[ch][1][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepC);
    Hists[ch][1][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepC);
    Hists[ch][1][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepC);
    Hists[ch][1][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepC);
    Hists[ch][1][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepC);
    Hists[ch][1][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepC);
    Hists[ch][1][18]->Fill(Topmass,weight_lepB);
    Hists[ch][1][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepC);
    Hists[ch][1][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepC);
    Hists[ch][1][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepC);
    Hists[ch][1][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepC);
    if(selectedJets->size()>0) Hists[ch][1][23]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][1][24]->Fill((*selectedJets)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][1][25]->Fill((*selectedJets)[0]->phi_,weight_lep);
    Hists[ch][1][26]->Fill(selectedJets->size(),weight_lep);
    Hists[ch][1][27]->Fill(nbjet,weight_lepB);
    Hists[ch][1][28]->Fill(MET_pt0,weight_lep);
    Hists[ch][1][29]->Fill(MET_phi0,weight_lep);
    Hists[ch][1][30]->Fill(Pileup_nTrueInt,weight_lep);
    Hists[ch][1][31]->Fill(Zmass,weight_lepC);
    Hists[ch][1][32]->Fill(Zpt,weight_lepC);
    Hists[ch][1][33]->Fill(ZDr,weight_lepC);
    Hists[ch][1][34]->Fill(ZDphi,weight_lepC);
    Hists[ch][1][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][1][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepC);
    Hists[ch][1][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepC);
    Hists[ch][1][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][1][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][1][40]->Fill(Ht,weight_lepC);
    Hists[ch][1][41]->Fill(Ms,weight_lepC);
    Hists[ch][1][42]->Fill(ZlDr,weight_lepC);
    Hists[ch][1][43]->Fill(ZlDphi,weight_lepC);
    Hists[ch][1][44]->Fill(JeDr,weight_lepB);
    Hists[ch][1][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][1][46]->Fill(tM,weight_lepB);
    }
    //Off Z

    if(!OnZ){
    Hists[ch][2][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    Hists[ch][2][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    Hists[ch][2][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
    Hists[ch][2][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    Hists[ch][2][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    Hists[ch][2][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    Hists[ch][2][6]->Fill((*selectedLeptons)[2]->pt_,weight_lep);
    Hists[ch][2][7]->Fill((*selectedLeptons)[2]->eta_,weight_lep);
    Hists[ch][2][8]->Fill((*selectedLeptons)[2]->phi_,weight_lep);
    Hists[ch][2][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepC);
    Hists[ch][2][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepC);
    Hists[ch][2][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepC);
    Hists[ch][2][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepC);
    Hists[ch][2][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepC);
    Hists[ch][2][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepC);
    Hists[ch][2][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepC);
    Hists[ch][2][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepC);
    Hists[ch][2][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepC);
    Hists[ch][2][18]->Fill(Topmass,weight_lepB);
    Hists[ch][2][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepC);
    Hists[ch][2][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepC);
    Hists[ch][2][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepC);
    Hists[ch][2][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepC);
    if(selectedJets->size()>0) Hists[ch][2][23]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][2][24]->Fill((*selectedJets)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][2][25]->Fill((*selectedJets)[0]->phi_,weight_lep);
    Hists[ch][2][26]->Fill(selectedJets->size(),weight_lep);
    Hists[ch][2][27]->Fill(nbjet,weight_lepB);
    Hists[ch][2][28]->Fill(MET_pt0,weight_lep);
    Hists[ch][2][29]->Fill(MET_phi0,weight_lep);
    Hists[ch][2][30]->Fill(Pileup_nTrueInt,weight_lep);
    Hists[ch][2][31]->Fill(Zmass,weight_lepC);
    Hists[ch][2][32]->Fill(Zpt,weight_lepC);
    Hists[ch][2][33]->Fill(ZDr,weight_lepC);
    Hists[ch][2][34]->Fill(ZDphi,weight_lepC);
    Hists[ch][2][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][2][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepC);
    Hists[ch][2][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepC);
    Hists[ch][2][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][2][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][2][40]->Fill(Ht,weight_lepC);
    Hists[ch][2][41]->Fill(Ms,weight_lepC);
    Hists[ch][2][42]->Fill(ZlDr,weight_lepC);
    Hists[ch][2][43]->Fill(ZlDphi,weight_lepC);
    Hists[ch][2][44]->Fill(JeDr,weight_lepB);
    Hists[ch][2][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][2][46]->Fill(tM,weight_lepB);
    }

    if(nbjet==0 && !OnZ){
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
    Hists[ch][3][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB);
    Hists[ch][3][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][3][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][3][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][3][23]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][3][24]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][3][25]->Fill((*selectedJets)[0]->phi_,weight_lepB);
    Hists[ch][3][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][3][27]->Fill(nbjet,weight_lepB);
    Hists[ch][3][28]->Fill(MET_pt0,weight_lepB);
    Hists[ch][3][29]->Fill(MET_phi0,weight_lepB);
    Hists[ch][3][30]->Fill(Pileup_nTrueInt,weight_lepB);
    Hists[ch][3][31]->Fill(Zmass,weight_lepB);
    Hists[ch][3][32]->Fill(Zpt,weight_lepB);
    Hists[ch][3][33]->Fill(ZDr,weight_lepB);
    Hists[ch][3][34]->Fill(ZDphi,weight_lepB);
    Hists[ch][3][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][3][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][3][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][3][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][3][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][3][40]->Fill(Ht,weight_lepB);
    Hists[ch][3][41]->Fill(Ms,weight_lepB);
    Hists[ch][3][42]->Fill(ZlDr,weight_lepB);
    Hists[ch][3][43]->Fill(ZlDphi,weight_lepB);
    Hists[ch][3][44]->Fill(JeDr,weight_lepB);
    Hists[ch][3][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][3][46]->Fill(tM,weight_lepB);
    }
    if(nbjet==1 && !OnZ){
    Hists[ch][4][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][4][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][4][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][4][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][4][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][4][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][4][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][4][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][4][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][4][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][4][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][4][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][4][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][4][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][4][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][4][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][4][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][4][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][4][18]->Fill(Topmass,weight_lepB);
    Hists[ch][4][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB);
    Hists[ch][4][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][4][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][4][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][4][23]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][4][24]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][4][25]->Fill((*selectedJets)[0]->phi_,weight_lepB);
    Hists[ch][4][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][4][27]->Fill(nbjet,weight_lepB);
    Hists[ch][4][28]->Fill(MET_pt0,weight_lepB);
    Hists[ch][4][29]->Fill(MET_phi0,weight_lepB);
    Hists[ch][4][30]->Fill(Pileup_nTrueInt,weight_lepB);
    Hists[ch][4][31]->Fill(Zmass,weight_lepB);
    Hists[ch][4][32]->Fill(Zpt,weight_lepB);
    Hists[ch][4][33]->Fill(ZDr,weight_lepB);
    Hists[ch][4][34]->Fill(ZDphi,weight_lepB);
    Hists[ch][4][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][4][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][4][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][4][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][4][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][4][40]->Fill(Ht,weight_lepB);
    Hists[ch][4][41]->Fill(Ms,weight_lepB);
    Hists[ch][4][42]->Fill(ZlDr,weight_lepB);
    Hists[ch][4][43]->Fill(ZlDphi,weight_lepB);
    Hists[ch][4][44]->Fill(JeDr,weight_lepB);
    Hists[ch][4][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][4][46]->Fill(tM,weight_lepB);
    }
    if(nbjet>=2 && !OnZ){
    Hists[ch][5][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][5][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][5][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][5][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][5][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][5][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][5][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][5][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][5][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][5][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][5][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][5][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][5][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][5][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][5][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][5][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][5][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][5][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][5][18]->Fill(Topmass,weight_lepB);
    Hists[ch][5][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB);
    Hists[ch][5][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][5][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][5][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][5][23]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][5][24]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][5][25]->Fill((*selectedJets)[0]->phi_,weight_lepB);
    Hists[ch][5][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][5][27]->Fill(nbjet,weight_lepB);
    Hists[ch][5][28]->Fill(MET_pt0,weight_lepB);
    Hists[ch][5][29]->Fill(MET_phi0,weight_lepB);
    Hists[ch][5][30]->Fill(Pileup_nTrueInt,weight_lepB);
    Hists[ch][5][31]->Fill(Zmass,weight_lepB);
    Hists[ch][5][32]->Fill(Zpt,weight_lepB);
    Hists[ch][5][33]->Fill(ZDr,weight_lepB);
    Hists[ch][5][34]->Fill(ZDphi,weight_lepB);
    Hists[ch][5][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][5][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepC);
    Hists[ch][5][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepC);
    Hists[ch][5][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][5][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][5][40]->Fill(Ht,weight_lepC);
    Hists[ch][5][41]->Fill(Ms,weight_lepC);
    Hists[ch][5][42]->Fill(ZlDr,weight_lepC);
    Hists[ch][5][43]->Fill(ZlDphi,weight_lepC);
    Hists[ch][5][44]->Fill(JeDr,weight_lepB);
    Hists[ch][5][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][5][46]->Fill(tM,weight_lepB);
    }
    if(MET_pt0<20 && !OnZ){
    Hists[ch][6][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    Hists[ch][6][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    Hists[ch][6][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
    Hists[ch][6][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    Hists[ch][6][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    Hists[ch][6][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    Hists[ch][6][6]->Fill((*selectedLeptons)[2]->pt_,weight_lep);
    Hists[ch][6][7]->Fill((*selectedLeptons)[2]->eta_,weight_lep);
    Hists[ch][6][8]->Fill((*selectedLeptons)[2]->phi_,weight_lep);
    Hists[ch][6][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepC);
    Hists[ch][6][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepC);
    Hists[ch][6][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepC);
    Hists[ch][6][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepC);
    Hists[ch][6][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepC);
    Hists[ch][6][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepC);
    Hists[ch][6][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepC);
    Hists[ch][6][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepC);
    Hists[ch][6][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepC);
    Hists[ch][6][18]->Fill(Topmass,weight_lepB);
    Hists[ch][6][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepC);
    Hists[ch][6][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepC);
    Hists[ch][6][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepC);
    Hists[ch][6][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepC);
    if(selectedJets->size()>0) Hists[ch][6][23]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][6][24]->Fill((*selectedJets)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][6][25]->Fill((*selectedJets)[0]->phi_,weight_lep);
    Hists[ch][6][26]->Fill(selectedJets->size(),weight_lep);
    Hists[ch][6][27]->Fill(nbjet,weight_lepB);
    Hists[ch][6][28]->Fill(MET_pt0,weight_lep);
    Hists[ch][6][29]->Fill(MET_phi0,weight_lep);
    Hists[ch][6][30]->Fill(Pileup_nTrueInt,weight_lep);
    Hists[ch][6][31]->Fill(Zmass,weight_lepC);
    Hists[ch][6][32]->Fill(Zpt,weight_lepC);
    Hists[ch][6][33]->Fill(ZDr,weight_lepC);
    Hists[ch][6][34]->Fill(ZDphi,weight_lepC);
    Hists[ch][6][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][6][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepC);
    Hists[ch][6][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepC);
    Hists[ch][6][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][6][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][6][40]->Fill(Ht,weight_lepC);
    Hists[ch][6][41]->Fill(Ms,weight_lepC);
    Hists[ch][6][42]->Fill(ZlDr,weight_lepC);
    Hists[ch][6][43]->Fill(ZlDphi,weight_lepC);
    Hists[ch][6][44]->Fill(JeDr,weight_lepB);
    Hists[ch][6][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][6][46]->Fill(tM,weight_lepB);
    }

    if(MET_pt0>20 && !OnZ){
    Hists[ch][7][0]->Fill((*selectedLeptons)[0]->pt_,weight_lep);
    Hists[ch][7][1]->Fill((*selectedLeptons)[0]->eta_,weight_lep);
    Hists[ch][7][2]->Fill((*selectedLeptons)[0]->phi_,weight_lep);
    Hists[ch][7][3]->Fill((*selectedLeptons)[1]->pt_,weight_lep);
    Hists[ch][7][4]->Fill((*selectedLeptons)[1]->eta_,weight_lep);
    Hists[ch][7][5]->Fill((*selectedLeptons)[1]->phi_,weight_lep);
    Hists[ch][7][6]->Fill((*selectedLeptons)[2]->pt_,weight_lep);
    Hists[ch][7][7]->Fill((*selectedLeptons)[2]->eta_,weight_lep);
    Hists[ch][7][8]->Fill((*selectedLeptons)[2]->phi_,weight_lep);
    Hists[ch][7][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepC);
    Hists[ch][7][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepC);
    Hists[ch][7][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepC);
    Hists[ch][7][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepC);
    Hists[ch][7][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepC);
    Hists[ch][7][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepC);
    Hists[ch][7][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepC);
    Hists[ch][7][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepC);
    Hists[ch][7][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepC);
    Hists[ch][7][18]->Fill(Topmass,weight_lepB);
    Hists[ch][7][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepC);
    Hists[ch][7][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepC);
    Hists[ch][7][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepC);
    Hists[ch][7][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepC);
    if(selectedJets->size()>0) Hists[ch][7][23]->Fill((*selectedJets)[0]->pt_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][7][24]->Fill((*selectedJets)[0]->eta_,weight_lep);
    if(selectedJets->size()>0) Hists[ch][7][25]->Fill((*selectedJets)[0]->phi_,weight_lep);
    Hists[ch][7][26]->Fill(selectedJets->size(),weight_lep);
    Hists[ch][7][27]->Fill(nbjet,weight_lepB);
    Hists[ch][7][28]->Fill(MET_pt0,weight_lep);
    Hists[ch][7][29]->Fill(MET_phi0,weight_lep);
    Hists[ch][7][30]->Fill(Pileup_nTrueInt,weight_lep);
    Hists[ch][7][31]->Fill(Zmass,weight_lepC);
    Hists[ch][7][32]->Fill(Zpt,weight_lepC);
    Hists[ch][7][33]->Fill(ZDr,weight_lepC);
    Hists[ch][7][34]->Fill(ZDphi,weight_lepC);
    Hists[ch][7][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][7][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepC);
    Hists[ch][7][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepC);
    Hists[ch][7][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][7][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepC);
    Hists[ch][7][40]->Fill(Ht,weight_lepC);
    Hists[ch][7][41]->Fill(Ms,weight_lepC);
    Hists[ch][7][42]->Fill(ZlDr,weight_lepC);
    Hists[ch][7][43]->Fill(ZlDphi,weight_lepC);
    Hists[ch][7][44]->Fill(JeDr,weight_lepB);
    Hists[ch][7][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][7][46]->Fill(tM,weight_lepB);
    }

    if(MET_pt0>20 && !OnZ && selectedJets->size()>0){
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
    Hists[ch][8][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB);
    Hists[ch][8][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][8][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][8][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][8][23]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][8][24]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][8][25]->Fill((*selectedJets)[0]->phi_,weight_lepB);
    Hists[ch][8][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][8][27]->Fill(nbjet,weight_lepB);
    Hists[ch][8][28]->Fill(MET_pt0,weight_lepB);
    Hists[ch][8][29]->Fill(MET_phi0,weight_lepB);
    Hists[ch][8][30]->Fill(Pileup_nTrueInt,weight_lepB);
    Hists[ch][8][31]->Fill(Zmass,weight_lepB);
    Hists[ch][8][32]->Fill(Zpt,weight_lepB);
    Hists[ch][8][33]->Fill(ZDr,weight_lepB);
    Hists[ch][8][34]->Fill(ZDphi,weight_lepB);
    Hists[ch][8][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][8][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][8][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][8][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][8][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][8][40]->Fill(Ht,weight_lepB);
    Hists[ch][8][41]->Fill(Ms,weight_lepB);
    Hists[ch][8][42]->Fill(ZlDr,weight_lepB);
    Hists[ch][8][43]->Fill(ZlDphi,weight_lepB);
    Hists[ch][8][44]->Fill(JeDr,weight_lepB);
    Hists[ch][8][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][8][46]->Fill(tM,weight_lepB);
    }

    if(nbjet==1 && MET_pt0>20 && selectedJets->size()==1 && !OnZ){
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
    Hists[ch][9][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB);
    Hists[ch][9][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][9][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][9][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][9][23]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][9][24]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][9][25]->Fill((*selectedJets)[0]->phi_,weight_lepB);
    Hists[ch][9][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][9][27]->Fill(nbjet,weight_lepB);
    Hists[ch][9][28]->Fill(MET_pt0,weight_lepB);
    Hists[ch][9][29]->Fill(MET_phi0,weight_lepB);
    Hists[ch][9][30]->Fill(Pileup_nTrueInt,weight_lepB);
    Hists[ch][9][31]->Fill(Zmass,weight_lepB);
    Hists[ch][9][32]->Fill(Zpt,weight_lepB);
    Hists[ch][9][33]->Fill(ZDr,weight_lepB);
    Hists[ch][9][34]->Fill(ZDphi,weight_lepB);
    Hists[ch][9][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][9][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][9][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][9][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][9][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][9][40]->Fill(Ht,weight_lepB);
    Hists[ch][9][41]->Fill(Ms,weight_lepB);
    Hists[ch][9][42]->Fill(ZlDr,weight_lepB);
    Hists[ch][9][43]->Fill(ZlDphi,weight_lepB);
    Hists[ch][9][44]->Fill(JeDr,weight_lepB);
    Hists[ch][9][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][9][46]->Fill(tM,weight_lepB);
    }

    if(nbjet==1 && MET_pt0>20 && selectedJets->size()==2 && !OnZ){
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
    Hists[ch][10][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB);
    Hists[ch][10][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][10][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][10][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][10][23]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][10][24]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][10][25]->Fill((*selectedJets)[0]->phi_,weight_lepB);
    Hists[ch][10][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][10][27]->Fill(nbjet,weight_lepB);
    Hists[ch][10][28]->Fill(MET_pt0,weight_lepB);
    Hists[ch][10][29]->Fill(MET_phi0,weight_lepB);
    Hists[ch][10][30]->Fill(Pileup_nTrueInt,weight_lepB);
    Hists[ch][10][31]->Fill(Zmass,weight_lepB);
    Hists[ch][10][32]->Fill(Zpt,weight_lepB);
    Hists[ch][10][33]->Fill(ZDr,weight_lepB);
    Hists[ch][10][34]->Fill(ZDphi,weight_lepB);
    Hists[ch][10][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][10][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][10][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][10][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][10][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][10][40]->Fill(Ht,weight_lepB);
    Hists[ch][10][41]->Fill(Ms,weight_lepB);
    Hists[ch][10][42]->Fill(ZlDr,weight_lepB);
    Hists[ch][10][43]->Fill(ZlDphi,weight_lepB);
    Hists[ch][10][44]->Fill(JeDr,weight_lepB);
    Hists[ch][10][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][10][46]->Fill(tM,weight_lepB);
    }

    if(MET_pt0>20&&!OnZ&&nbjet<=1&&selectedJets->size()>=1){
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
    Hists[ch][11][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB);
    Hists[ch][11][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][11][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][11][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][11][23]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][11][24]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][11][25]->Fill((*selectedJets)[0]->phi_,weight_lepB);
    Hists[ch][11][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][11][27]->Fill(nbjet,weight_lepB);
    Hists[ch][11][28]->Fill(MET_pt0,weight_lepB);
    Hists[ch][11][29]->Fill(MET_phi0,weight_lepB);
    Hists[ch][11][30]->Fill(Pileup_nTrueInt,weight_lepB);
    Hists[ch][11][31]->Fill(Zmass,weight_lepB);
    Hists[ch][11][32]->Fill(Zpt,weight_lepB);
    Hists[ch][11][33]->Fill(ZDr,weight_lepB);
    Hists[ch][11][34]->Fill(ZDphi,weight_lepB);
    Hists[ch][11][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][11][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][11][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][11][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][11][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][11][40]->Fill(Ht,weight_lepB);
    Hists[ch][11][41]->Fill(Ms,weight_lepB);
    Hists[ch][11][42]->Fill(ZlDr,weight_lepB);
    Hists[ch][11][43]->Fill(ZlDphi,weight_lepB);
    Hists[ch][11][44]->Fill(JeDr,weight_lepB);
    Hists[ch][11][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][11][46]->Fill(tM,weight_lepB);
    }

    if(MET_pt0>20&&OnZ&&nbjet<=1&&selectedJets->size()>=1){
    Hists[ch][12][0]->Fill((*selectedLeptons)[0]->pt_,weight_lepB);
    Hists[ch][12][1]->Fill((*selectedLeptons)[0]->eta_,weight_lepB);
    Hists[ch][12][2]->Fill((*selectedLeptons)[0]->phi_,weight_lepB);
    Hists[ch][12][3]->Fill((*selectedLeptons)[1]->pt_,weight_lepB);
    Hists[ch][12][4]->Fill((*selectedLeptons)[1]->eta_,weight_lepB);
    Hists[ch][12][5]->Fill((*selectedLeptons)[1]->phi_,weight_lepB);
    Hists[ch][12][6]->Fill((*selectedLeptons)[2]->pt_,weight_lepB);
    Hists[ch][12][7]->Fill((*selectedLeptons)[2]->eta_,weight_lepB);
    Hists[ch][12][8]->Fill((*selectedLeptons)[2]->phi_,weight_lepB);
    Hists[ch][12][9]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB);
    Hists[ch][12][10]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB);
    Hists[ch][12][11]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB);
    Hists[ch][12][12]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB);
    Hists[ch][12][13]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB);
    Hists[ch][12][14]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB);
    Hists[ch][12][15]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB);
    Hists[ch][12][16]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB);
    Hists[ch][12][17]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB);
    Hists[ch][12][18]->Fill(Topmass,weight_lepB);
    Hists[ch][12][19]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB);
    Hists[ch][12][20]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][12][21]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB);
    Hists[ch][12][22]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB);
    if(selectedJets->size()>0) Hists[ch][12][23]->Fill((*selectedJets)[0]->pt_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][12][24]->Fill((*selectedJets)[0]->eta_,weight_lepB);
    if(selectedJets->size()>0) Hists[ch][12][25]->Fill((*selectedJets)[0]->phi_,weight_lepB);
    Hists[ch][12][26]->Fill(selectedJets->size(),weight_lepB);
    Hists[ch][12][27]->Fill(nbjet,weight_lepB);
    Hists[ch][12][28]->Fill(MET_pt0,weight_lepB);
    Hists[ch][12][29]->Fill(MET_phi0,weight_lepB);
    Hists[ch][12][30]->Fill(Pileup_nTrueInt,weight_lepB);
    Hists[ch][12][31]->Fill(Zmass,weight_lepB);
    Hists[ch][12][32]->Fill(Zpt,weight_lepB);
    Hists[ch][12][33]->Fill(ZDr,weight_lepB);
    Hists[ch][12][34]->Fill(ZDphi,weight_lepB);
    Hists[ch][12][35]->Fill(LFVTopmass,weight_lepB);
    Hists[ch][12][36]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB);
    Hists[ch][12][37]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB);
    Hists[ch][12][38]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][12][39]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB);
    Hists[ch][12][40]->Fill(Ht,weight_lepB);
    Hists[ch][12][41]->Fill(Ms,weight_lepB);
    Hists[ch][12][42]->Fill(ZlDr,weight_lepB);
    Hists[ch][12][43]->Fill(ZlDphi,weight_lepB);
    Hists[ch][12][44]->Fill(JeDr,weight_lepB);
    Hists[ch][12][45]->Fill(JmuDr,weight_lepB);
    Hists[ch][12][46]->Fill(tM,weight_lepB);
    }

    if (nbjet<=1&&MET_pt0>20&&selectedJets->size()>0){
    for (int n=0;n<8;++n){
    if(!OnZ){
    HistsSysUp[ch][0][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][6][n]->Fill((*selectedLeptons)[2]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][7][n]->Fill((*selectedLeptons)[2]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][8][n]->Fill((*selectedLeptons)[2]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][9][n]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][10][n]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][11][n]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][12][n]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][13][n]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][14][n]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][15][n]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][16][n]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][17][n]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][18][n]->Fill(Topmass,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][19][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][20][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][21][n]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][22][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysUp[ch][0][23][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysUp[ch][0][24][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysUp[ch][0][25][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][26][n]->Fill(selectedJets->size(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][27][n]->Fill(nbjet,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][28][n]->Fill(MET_pt0,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][29][n]->Fill(MET_phi0,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][30][n]->Fill(Pileup_nTrueInt,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][31][n]->Fill(Zmass,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][32][n]->Fill(Zpt,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][33][n]->Fill(ZDr,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][34][n]->Fill(ZDphi,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][35][n]->Fill(LFVTopmass,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][36][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][37][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][38][n]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][39][n]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][40][n]->Fill(Ht,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][41][n]->Fill(Ms,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][42][n]->Fill(ZlDr,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][43][n]->Fill(ZlDphi,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][44][n]->Fill(JeDr,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][45][n]->Fill(JmuDr,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][0][46][n]->Fill(tM,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));

    HistsSysDown[ch][0][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][6][n]->Fill((*selectedLeptons)[2]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][7][n]->Fill((*selectedLeptons)[2]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][8][n]->Fill((*selectedLeptons)[2]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][9][n]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][10][n]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][11][n]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][12][n]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][13][n]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][14][n]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][15][n]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][16][n]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][17][n]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][18][n]->Fill(Topmass,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][19][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][20][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][21][n]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][22][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysDown[ch][0][23][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysDown[ch][0][24][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysDown[ch][0][25][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][26][n]->Fill(selectedJets->size(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][27][n]->Fill(nbjet,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][28][n]->Fill(MET_pt0,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][29][n]->Fill(MET_phi0,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][30][n]->Fill(Pileup_nTrueInt,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][31][n]->Fill(Zmass,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][32][n]->Fill(Zpt,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][33][n]->Fill(ZDr,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][34][n]->Fill(ZDphi,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][35][n]->Fill(LFVTopmass,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][36][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][37][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][38][n]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][39][n]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][40][n]->Fill(Ht,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][41][n]->Fill(Ms,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][42][n]->Fill(ZlDr,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][43][n]->Fill(ZlDphi,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][44][n]->Fill(JeDr,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][45][n]->Fill(JmuDr,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][0][46][n]->Fill(tM,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    }
    else{
    HistsSysUp[ch][1][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][6][n]->Fill((*selectedLeptons)[2]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][7][n]->Fill((*selectedLeptons)[2]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][8][n]->Fill((*selectedLeptons)[2]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][9][n]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][10][n]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][11][n]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][12][n]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][13][n]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][14][n]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][15][n]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][16][n]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][17][n]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][18][n]->Fill(Topmass,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][19][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][20][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][21][n]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][22][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysUp[ch][1][23][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysUp[ch][1][24][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysUp[ch][1][25][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][26][n]->Fill(selectedJets->size(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][27][n]->Fill(nbjet,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][28][n]->Fill(MET_pt0,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][29][n]->Fill(MET_phi0,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][30][n]->Fill(Pileup_nTrueInt,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][31][n]->Fill(Zmass,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][32][n]->Fill(Zpt,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][33][n]->Fill(ZDr,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][34][n]->Fill(ZDphi,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][35][n]->Fill(LFVTopmass,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][36][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][37][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][38][n]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][39][n]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][40][n]->Fill(Ht,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][41][n]->Fill(Ms,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][42][n]->Fill(ZlDr,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][43][n]->Fill(ZlDphi,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][44][n]->Fill(JeDr,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][45][n]->Fill(JmuDr,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));
    HistsSysUp[ch][1][46][n]->Fill(tM,weight_lepB * (sysUpWeights[n]/nominalWeights[n]));

    HistsSysDown[ch][1][0][n]->Fill((*selectedLeptons)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][1][n]->Fill((*selectedLeptons)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][2][n]->Fill((*selectedLeptons)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][3][n]->Fill((*selectedLeptons)[1]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][4][n]->Fill((*selectedLeptons)[1]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][5][n]->Fill((*selectedLeptons)[1]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][6][n]->Fill((*selectedLeptons)[2]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][7][n]->Fill((*selectedLeptons)[2]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][8][n]->Fill((*selectedLeptons)[2]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][9][n]->Fill((*selectedLeptons_copy)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][10][n]->Fill((*selectedLeptons_copy)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][11][n]->Fill((*selectedLeptons_copy)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][12][n]->Fill((*selectedLeptons_copy)[1]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][13][n]->Fill((*selectedLeptons_copy)[1]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][14][n]->Fill((*selectedLeptons_copy)[1]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][15][n]->Fill((*selectedLeptons_copy)[2]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][16][n]->Fill((*selectedLeptons_copy)[2]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][17][n]->Fill((*selectedLeptons_copy)[2]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][18][n]->Fill(Topmass,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][19][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][20][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][21][n]->Fill(((*selectedLeptons_copy)[0]->p4_).Pt()+((*selectedLeptons_copy)[1]->p4_).Pt()+((*selectedLeptons_copy)[2]->p4_).Pt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][22][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_ + (*selectedLeptons_copy)[2]->p4_).Mt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysDown[ch][1][23][n]->Fill((*selectedJets)[0]->pt_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysDown[ch][1][24][n]->Fill((*selectedJets)[0]->eta_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    if(selectedJets->size()>0) HistsSysDown[ch][1][25][n]->Fill((*selectedJets)[0]->phi_,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][26][n]->Fill(selectedJets->size(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][27][n]->Fill(nbjet,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][28][n]->Fill(MET_pt0,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][29][n]->Fill(MET_phi0,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][30][n]->Fill(Pileup_nTrueInt,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][31][n]->Fill(Zmass,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][32][n]->Fill(Zpt,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][33][n]->Fill(ZDr,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][34][n]->Fill(ZDphi,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][35][n]->Fill(LFVTopmass,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][36][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).M(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][37][n]->Fill(((*selectedLeptons_copy)[0]->p4_ + (*selectedLeptons_copy)[1]->p4_).Pt(),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][38][n]->Fill(deltaR((*selectedLeptons_copy)[0]->eta_,(*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->eta_,(*selectedLeptons_copy)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][39][n]->Fill(deltaPhi((*selectedLeptons_copy)[0]->phi_,(*selectedLeptons_copy)[1]->phi_),weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][40][n]->Fill(Ht,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][41][n]->Fill(Ms,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][42][n]->Fill(ZlDr,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][43][n]->Fill(ZlDphi,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][44][n]->Fill(JeDr,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][45][n]->Fill(JmuDr,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    HistsSysDown[ch][1][46][n]->Fill(tM,weight_lepB * (sysDownWeights[n]/nominalWeights[n]));
    }
    }
    }
      
    }

    //AR1:ll\overline{l}
    if (nTight==2){

    anti_index=(*selectedLeptons)[1]->isTag?2:1;
    anti_flavor=(*selectedLeptons)[anti_index]->lep_>1?1:0;
    if(abs((*selectedLeptons)[anti_index]->eta_)>1.6){
      etabin1=2;
    }
    else if(abs((*selectedLeptons)[anti_index]->eta_)>0.8){
      etabin1=1;
    }

    if (MET_pt0>20&&!OnZ&&nbjet<=1&&selectedJets->size()>=1){
    HistsAR1[anti_flavor][0][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    HistsAR1[anti_flavor][0][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    if(selectedJets->size()>0) HistsAR1[anti_flavor][0][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    if(selectedJets->size()>0) HistsAR1[anti_flavor][0][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    HistsAR1[anti_flavor][0][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    HistsAR1[anti_flavor][0][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    HistsAR1[anti_flavor][0][etabin1][6]->Fill(MET_pt0,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    HistsAR1[anti_flavor][0][etabin1][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    }

    if (MET_pt0>20&&OnZ&&nbjet<=1&&selectedJets->size()>=1){
    HistsAR1[anti_flavor][1][etabin1][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    HistsAR1[anti_flavor][1][etabin1][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    if(selectedJets->size()>0) HistsAR1[anti_flavor][1][etabin1][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    if(selectedJets->size()>0) HistsAR1[anti_flavor][1][etabin1][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    HistsAR1[anti_flavor][1][etabin1][4]->Fill(selectedJets->size(),(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    HistsAR1[anti_flavor][1][etabin1][5]->Fill(nbjet,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    HistsAR1[anti_flavor][1][etabin1][6]->Fill(MET_pt0,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    HistsAR1[anti_flavor][1][etabin1][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[anti_index]->pt_,weight_lepB);
    }
    }
    //AR2:l\overline{l}\overline{l}
    if (nTight==1){

    if ((*selectedLeptons)[1]->lep_+(*selectedLeptons)[2]->lep_==2||(*selectedLeptons)[1]->lep_+(*selectedLeptons)[2]->lep_==20){
      anti_flavor=(*selectedLeptons)[1]->lep_>1?1:0;
    }
    else{
      anti_flavor=2;
      sort(selectedLeptons->begin()+1, selectedLeptons->end(), CompareFlavourLep);
    }

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

    if (MET_pt0>20&&!OnZ&&nbjet<=1&&selectedJets->size()>=1){
    HistsAR2[anti_flavor][0][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    HistsAR2[anti_flavor][0][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    if(selectedJets->size()>0) HistsAR2[anti_flavor][0][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    if(selectedJets->size()>0) HistsAR2[anti_flavor][0][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    HistsAR2[anti_flavor][0][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    HistsAR2[anti_flavor][0][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    HistsAR2[anti_flavor][0][etabin1][etabin2][6]->Fill(MET_pt0,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    HistsAR2[anti_flavor][0][etabin1][etabin2][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    }

    if (MET_pt0>20&&OnZ&&nbjet<=1&&selectedJets->size()>=1){
    HistsAR2[anti_flavor][1][etabin1][etabin2][0]->Fill((*selectedLeptons)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    HistsAR2[anti_flavor][1][etabin1][etabin2][1]->Fill((*selectedLeptons)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    if(selectedJets->size()>0) HistsAR2[anti_flavor][1][etabin1][etabin2][2]->Fill((*selectedJets)[0]->pt_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    if(selectedJets->size()>0) HistsAR2[anti_flavor][1][etabin1][etabin2][3]->Fill((*selectedJets)[0]->eta_,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    HistsAR2[anti_flavor][1][etabin1][etabin2][4]->Fill(selectedJets->size(),(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    HistsAR2[anti_flavor][1][etabin1][etabin2][5]->Fill(nbjet,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    HistsAR2[anti_flavor][1][etabin1][etabin2][6]->Fill(MET_pt0,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    HistsAR2[anti_flavor][1][etabin1][etabin2][7]->Fill(Pileup_nTrueInt,(*selectedLeptons)[1]->pt_,(*selectedLeptons)[2]->pt_,weight_lepB);
    }
    }
      
    for (int l=0;l<(int)selectedLeptons->size();l++){
      delete (*selectedLeptons)[l];
      delete (*selectedLeptons_copy)[l];
    }
    for (int l=0;l<(int)selectedJets->size();l++){
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

  for (int i=0;i<(int)channels.size();++i){
    for (int k=0;k<(int)regions.size();++k){
      for (int l=0;l<(int)vars.size();++l){
        Hists[i][k][l]  ->Write("",TObject::kOverwrite);
        for (int n=0;n<(int)sys.size();++n){
	    if (k==0){//Only write the first set of histograms to save space
            HistsSysUp[i][k][l][n]->Write("",TObject::kOverwrite);
            HistsSysDown[i][k][l][n]->Write("",TObject::kOverwrite);
	     }
        }
      }
    }
  }
    
  for (int i=0;i<(int)channelsFF.size();++i){
      for (int k=0;k<2;++k){
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
      for (int l=0;l<(int)vars.size();++l){
          delete Hists[i][k][l];
        for (int n=0;n<(int)sys.size();++n){
	        delete HistsSysUp[i][k][l][n];
	        delete HistsSysDown[i][k][l][n];
        }
      }
    }
  }
    
  for (int i=0;i<(int)channelsFF.size();++i){
      for (int k=0;k<2;++k){
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
