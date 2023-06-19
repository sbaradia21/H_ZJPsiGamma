#include <iostream>
#include <memory>
#include "MVAnalysis.h"
using std::string;
using std::cout;
using std::endl;

MVAnalysis::MVAnalysis(const string& mva_algo, const string& xmlfile) 
{
  reader_ = std::make_unique<TMVA::Reader>("!Color:!Silent");
  reader_->AddVariable("phopt",&varList_.phopt);
  //  reader_->AddVariable("phoeta",&varList_.phoeta);
  //  reader_->AddVariable("phophi",&varList_.phophi);
  //reader_->AddVariable("photon_r9",&varList_.photon_r9);
  reader_->AddVariable("mu1pt",&varList_.mu1pt);
  //reader_->AddVariable("mu1eta",&varList_.mu1eta);
  // reader_->AddVariable("mu1phi",&varList_.mu1phi);
  reader_->AddVariable("mu2pt",&varList_.mu2pt);
  //reader_->AddVariable("mu2eta",&varList_.mu2eta);
  //reader_->AddVariable("mu2phi",&varList_.mu2phi);
  reader_->AddVariable("mu1iso",&varList_.mu1iso);
  //reader_->AddVariable("mu2iso",&varList_.mu2iso);
  // reader_->AddVariable("mu1_sip3d",&varList_.mu1_sip3d);
  //reader_->AddVariable("mu2_sip3d",&varList_.mu2_sip3d);
  //reader_->AddVariable("dphijpsipho",&varList_.dphijpsipho);
  reader_->AddVariable("mod_dphijpsipho",&varList_.mod_dphijpsipho);
  //reader_->AddVariable("dphimu1pho",&varList_.dphimu1pho);
  ////  reader_->AddVariable("mod_dphimu1pho",&varList_.mod_dphimu1pho);
  //reader_->AddVariable("dphimu2pho",&varList_.dphimu2pho);
  ////reader_->AddVariable("mod_dphimu2pho",&varList_.mod_dphimu2pho);
  //reader_->AddVariable("dphimu1mu2",&varList_.dphimu1mu2);                                                                                                                                          
  ////reader_->AddVariable("mod_dphimu1mu2",&varList_.mod_dphimu1mu2); 
  reader_->AddVariable("drmu1pho",&varList_.drmu1pho);
  reader_->AddVariable("drmu2pho",&varList_.drmu2pho);
  reader_->AddVariable("drjpsipho",&varList_.drjpsipho);
  //  reader_->AddVariable("drmu1mu2",&varList_.drmu1mu2);
  // reader_->AddVariable("jpsipt",&varList_.jpsipt);
  //reader_->AddVariable("Hm",&varList_.Hm);
  //reader_->AddVariable("Hmpt",&varList_.Hmpt);
  //  reader_->AddVariable("ptratio",&varList_.ptratio);
  //  reader_->AddVariable("etratio",&varList_.etratio);
  //reader_->AddVariable("jpsim",&varList_.jpsim);
  //reader_->AddVariable("tlumi_wt",&varList_.tlumi_wt);
  //reader_->AddVariable("tlumi_wt_hasGen",&varList_.tlumi_wt_hasGen);
  //  reader_->AddVariable("jpsiSSpt",&varList_.jpsiSSpt);
  //reader_->AddVariable("HmSSpt",&varList_.HmSSpt);
  reader_->BookMVA(mva_algo.c_str(), xmlfile);
}
double MVAnalysis::evaluate(const string& mva_algo, const InputVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList)); // use move syntax here
  return reader_->EvaluateMVA(mva_algo.c_str());
}
