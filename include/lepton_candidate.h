#ifndef MY_lepton_candidate
#define MY_lepton_candidate

#include<cmath>
#include<string>
#include<iostream>
#include<vector>
#include<complex>
#include <TLorentzVector.h>

using namespace std;
//using namespace math;
class lepton_candidate {
  
public:
  lepton_candidate(float, float, float, int, int, int );
  ~lepton_candidate();
  float pt_;
  float eta_;
  float phi_;
  int charge_;
  int indice_;
  int lep_;
  int isbalep;//bachelor lepton: lepton comeing from standard top
  int isantilep;//anti-selected lepton
  void setbalep(){
       isbalep=1;
  }
  void setantilep(){
       isantilep=1;
  }
  TLorentzVector p4_;


private:
  
};

#endif

