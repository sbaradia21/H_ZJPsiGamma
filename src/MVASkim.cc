#include <iostream>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "MVASkim.h"

using std::string;
using std::cout;
using std::endl;

MVASkim::MVASkim(const string& filename) {
  TFile* f = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  mvaFile_.reset(std::move(f));
  tree_ = new TTree("FlatTree", "Flat Tree for pre-selected events");
  tree_->Branch("phopt",&varList_.phopt,"phopt/F");
  tree_->Branch("phoeta",&varList_.phoeta,"phoeta/F");
  tree_->Branch("phophi",&varList_.phophi,"phophi/F");
  tree_->Branch("photon_r9",&varList_.photon_r9,"photon_r9/F");
  tree_->Branch("mu1pt",&varList_.mu1pt,"mu1pt/F");
  tree_->Branch("mu1eta",&varList_.mu1eta,"mu1eta/F");
  tree_->Branch("mu1phi",&varList_.mu1phi,"mu1phi/F");
  tree_->Branch("mu2pt",&varList_.mu2pt,"mu2pt/F");
  tree_->Branch("mu2eta",&varList_.mu2eta,"mu2eta/F");
  tree_->Branch("mu1iso",&varList_.mu1iso,"mu1iso/F");
  tree_->Branch("mu2iso",&varList_.mu2iso,"mu2iso/F");
  tree_->Branch("mu1_sip3d",&varList_.mu1_sip3d,"mu1_sip3d/F");
  tree_->Branch("mu2_sip3d",&varList_.mu2_sip3d,"mu2_sip3d/F");
  tree_->Branch("dphijpsipho",&varList_.dphijpsipho,"dphijpsipho/F");
  tree_->Branch("mod_dphijpsipho",&varList_.mod_dphijpsipho,"mod_dphijpsipho/F");
  tree_->Branch("dphimu1pho",&varList_.dphimu1pho,"dphimu1pho/F");
  tree_->Branch("mod_dphimu1pho",&varList_.mod_dphimu1pho,"mod_dphimu1pho/F");
  tree_->Branch("dphimu2pho",&varList_.dphimu2pho,"dphimu2pho/F");
  tree_->Branch("mod_dphimu2pho",&varList_.mod_dphimu2pho,"mod_dphimu2pho/F");
  tree_->Branch("dphimu1mu2",&varList_.dphimu1mu2,"dphimu1mu2/F");
  tree_->Branch("mod_dphimu1mu2",&varList_.mod_dphimu1mu2,"mod_dphimu1mu2/F");
  tree_->Branch("drmu1pho",&varList_.drmu1pho,"drmu1pho/F");
  tree_->Branch("drmu2pho",&varList_.drmu2pho,"drmu2pho/F");
  tree_->Branch("drjpsipho",&varList_.drjpsipho,"drjpsipho/F");
  tree_->Branch("drmu1mu2",&varList_.drmu1mu2,"drmu1mu2/F");
  tree_->Branch("jpsim",&varList_.jpsim,"jpsim/F");
  tree_->Branch("jpsipt",&varList_.jpsipt,"jpsipt/F");
  tree_->Branch("Hm",&varList_.Hm,"Hm/F");
  tree_->Branch("Hmpt",&varList_.Hmpt,"Hmpt/F");
  tree_->Branch("ptratio",&varList_.ptratio,"ptratio/F");
  tree_->Branch("etratio",&varList_.etratio,"etratio/F");
  tree_->Branch("jpsiSSpt",&varList_.jpsiSSpt,"jpsiSSpt/F");
  tree_->Branch("HmSSpt",&varList_.HmSSpt,"HmSSpt/F");
  tree_->Branch("tlumi_wt",&varList_.tlumi_wt,"tlumi_wt/F");
  tree_->Branch("tlumi_wt_hasGen",&varList_.tlumi_wt_hasGen,"tlumi_wt_hasGen/F");
  mvaFile_->ls();
}
MVASkim::~MVASkim() {
}
void MVASkim::fill(const TreeVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList));
  mvaFile_->cd();
  tree_->Fill();  
}
void MVASkim::close() {
  mvaFile_->cd();
  tree_->Print();
  tree_->Write();
  mvaFile_->Save();
  mvaFile_->Write();
  mvaFile_->Close();
}
