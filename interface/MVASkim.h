#ifndef __MVASkim__h
#define __MVASkim__h

#include <fstream>
#include <string>

class TTree;
class TFile;

typedef struct  
{
  float phopt;
  float phoeta;
  float phophi;
  float photon_r9;
  float mu1pt;
  float mu1eta;
  float mu1phi;
  float mu2pt;
  float mu2eta;
  float mu2phi;
  float mu1iso;
  float mu2iso;
  float mu1_sip3d;
  float mu2_sip3d;
  float dphijpsipho;
  float mod_dphijpsipho; 
  float dphimu1pho;
  float mod_dphimu1pho;
  float dphimu2pho ;
  float mod_dphimu2pho;
  float dphimu1mu2;
  float mod_dphimu1mu2;
  float drmu1pho;
  float drmu2pho ;
  float drjpsipho;
  float drmu1mu2;
  float jpsim;
  float jpsipt;
  float Hm;
  float Hmpt ;
  float ptratio;
  float etratio;
  float jpsiSSpt;
  float HmSSpt;
  float tlumi_wt;
  float tlumi_wt_hasGen;
} TreeVariables;

class MVASkim {
    
public:

  MVASkim(const std::string& filename);
  virtual ~MVASkim();

  void fill(const TreeVariables& varList);
  void close();

  std::unique_ptr<TFile> mvaFile_; // The output file with a flat tree for MVA study
  TTree* tree_;

  TreeVariables varList_;
};
#endif
