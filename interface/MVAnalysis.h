#ifndef __MVAnalysis__h
#define __MVAnalysis__h

#include <fstream>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

typedef struct  
{
  float phopt;
  //  float phoeta;
  //  float phophi;
  //float photon_r9;
  float mu1pt;
  //float mu1eta;
  //float mu1phi;
  float mu2pt;
  // float mu2eta;
  //float mu2phi;
  float mu1iso;
  //float mu2iso;
  //float mu1_sip3d;
  //float mu2_sip3d;
  //float dphijpsipho;
  float mod_dphijpsipho;
  //float dphimu1pho;
  ////  float mod_dphimu1pho;
  //float dphimu2pho ;
  ////float mod_dphimu2pho;
  // float dphimu1mu2;
  ////  float mod_dphimu1mu2;
  float drmu1pho;
  float drmu2pho ;
  float drjpsipho;
  ////  float drmu1mu2;
  // float jpsipt;
  //float Hm;
  //float Hmpt ;
  //  float ptratio;
  //  float etratio;
  //float jpsim;
  //float tlumi_wt;
  //float tlumi_wt_hasGen;
  //float jpsiSSpt;
  //  float HmSSpt;
}  InputVariables;


class MVAnalysis {
    
public:

  MVAnalysis(const std::string& mva_algo, const std::string& xmlfile);
  virtual ~MVAnalysis() {}

  double evaluate(const std::string& tag, const InputVariables& varList);

  InputVariables varList_;
  std::unique_ptr<TMVA::Reader> reader_;
};
#endif
