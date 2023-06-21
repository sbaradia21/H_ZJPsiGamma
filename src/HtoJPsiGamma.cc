#define HtoJPsiGamma_cxx
#include <chrono>
#include "HtoJPsiGamma.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TDirectory.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TString.h>
#include <assert.h>
#include <TFile.h>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1K.h"
#include "TObject.h"
#include "TH1D.h"
#include "TSystem.h"
#include "configana.h"
#include "correction.h"
#include "correctionlib_version.h"
#include "RoccoR.h"
#include "json.hpp"
using namespace std;
//using json = nlohmann::json;

void HtoJPsiGamma::tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
  // Skip delimiters at beginning.                                                                                                         
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // Find first "non-delimiter".                                                                                                            
  string::size_type pos = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)  {
    // Found a token, add it to the vector.                                                                                                 
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    
    // Skip delimiters.  Note the "not_of"                                                                                                  
    lastPos = str.find_first_not_of(delimiters, pos);
    
    // Find next "non-delimiter"                                                                                                            
    pos = str.find_first_of(delimiters, lastPos);
  }
}
void HtoJPsiGamma::setJson(const std::string gjsonF) {
  std::ifstream fin(gjsonF.c_str());
  fin >> datajson_; 
  fin.close();
}


bool HtoJPsiGamma::isGoodLumi(unsigned int run, unsigned int lumi) {
  std::string rkey = std::to_string(run);
  // find a run
  if (datajson_.find(rkey) != datajson_.end()) {//first check if run in json
    for( auto& lsL : datajson_.at(rkey)) {//these are the lumi blocks
      unsigned int low = lsL.at(0);
      unsigned int up = lsL.at(1);
      if(lumi >= low && lumi <=up)//Lumi is in range
        return true;//there is no else part since one has loop over all LS blocks in a Run
    }//lsb loop
  }//if block 
  return false;
}

void HtoJPsiGamma::buildList(const vector<string>& tokens, vector<string>& list) {
  for (vector<string>::const_iterator it  = tokens.begin()+1; it != tokens.end(); ++it){
    list.push_back(*it);
  }
}

int HtoJPsiGamma::setInputFile(const string& fname){
  auto found = fname.find("root:");
  if (found == string::npos && gSystem->AccessPathName(fname.c_str())){
    cerr << ">>> Warning: File <<" << fname << ">> was not found!!" << endl;
    return static_cast<int>(chain->GetEntries());
  }
  chain->Add(fname.c_str(), -1);
  chain1->Add(fname.c_str(), -1);
  return static_cast<int>(chain->GetEntries());

}


TH1* HtoJPsiGamma::getHist1D(const char* hname) {
  //  outfile->cd("HtoJPsiGamma");
  TObject *obj = gDirectory->GetList()->FindObject(hname);
  if (obj == nullptr) {
    cerr << "**** getHist1D: Histogram for <" << hname
	 << "> not found! ("
	 << __FILE__ << ":" << __LINE__ << ")"
	 << endl;
    return nullptr;
  }
  TH1* h = nullptr;
  if (obj->InheritsFrom("TH1D"))
    h = dynamic_cast<TH1D*>(obj);
  else if (obj->InheritsFrom("TH1C"))
    h = dynamic_cast<TH1C*>(obj);
  else if (obj->InheritsFrom("TH1K"))
    h = dynamic_cast<TH1K*>(obj);
  else if (obj->InheritsFrom("TH1S"))
    h = dynamic_cast<TH1S*>(obj);
  else if (obj->InheritsFrom("TH1I"))
    h = dynamic_cast<TH1I*>(obj);
  else
    h = dynamic_cast<TH1F*>(obj);
  if (h == nullptr) {
    cerr << "**** getHist1D: <" << hname
	 << "> may not be a 1D Histogram! ("
	 << __FILE__ << ":" << __LINE__ << ")"
	 << endl;
  }
  return h;
}

TH1* HtoJPsiGamma::getHist1D(const string& hname) {
  return getHist1D(hname.c_str());
}



TH2* HtoJPsiGamma::getHist2D(const char* hname) {
  TObject *obj = gDirectory->GetList()->FindObject(hname);
  if (obj == nullptr) {
    cerr << "**** getHist2D: Histogram for <" << hname
	 << "> not found! ("
	 << __FILE__ << ":" << __LINE__ << ")"
	 << endl;
    return nullptr;
  }
  TH2* h = nullptr;
  if (obj->InheritsFrom("TH2D"))
    h = dynamic_cast<TH2D*>(obj);
  else if (obj->InheritsFrom("TH2C"))
    h = dynamic_cast<TH2C*>(obj);
  else if (obj->InheritsFrom("TH2S"))
    h = dynamic_cast<TH2S*>(obj);
  else if (obj->InheritsFrom("TH2I"))
    h = dynamic_cast<TH2I*>(obj);
  else
    h = dynamic_cast<TH2F*>(obj);
  if (h == nullptr) {
    cerr << "**** getHist2D: <" << hname
	 << "> may not be a 2D Histogram! ("
	 << __FILE__ << ":" << __LINE__ << ")"
	 << endl;
  }
  return h;
}

TH2* getHist2D(const string& hname) {
  return getHist2D(hname.c_str());
}
double HtoJPsiGamma::lumiWt(double evtWeightSum, int nent) const
{
  double nevt = (evtWeightSum > -1) ? evtWeightSum : nent;
  //  if (!verbose)
  //  xsec = xsec_;
  
  std::cout << "-- intLumi: " << intLumi_
	    << " xsec: " << xsec_
	    << " nevt: " << nevt  <<std::endl;
  return (intLumi_ * xsec_ / nevt);
}

// taken from HZZ4lUtil to show cutflow efficiency
void HtoJPsiGamma::showEfficiency(const string& hname,const std::vector<std::string>& slist,const string& header,const string& tag,std::ostream& os)
{
  os << ">>> " << header << " Efficiency" << endl;
  TH1 *h = HtoJPsiGamma::getHist1D(hname);
  if (h != nullptr) {
    os << setw(64) << "CutFlow"
       << setw(20) << tag
       << setw(20) << "AbsEff"
       << setw(20) << "RelEff"
       << endl;
    os.precision(3);
    int nbins = h->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
      double cont  = h->GetBinContent(1);
      double conti = h->GetBinContent(i);
      double contj = h->GetBinContent(i-1);
      os << setw(64) << slist[i-1]
	 << setprecision(8)
	 << setw(20) << conti
	 << setprecision(5)
	 << setw(20) << ((conti > 0) ? conti/cont : 0.0)
	 << setw(20) << ( i == 1 ? 1.0 :(contj > 0) ? conti/contj : 0.0)
	 << endl;
    }
  }
}

void HtoJPsiGamma::objectEfficiency() {
  
  // Muon  
  vector<string> muLabels {
    "All",
      "looseid",
      "tightid"
      };
  HtoJPsiGamma::showEfficiency("muCutFlow", muLabels, "Muon Selection","Muons");

  // Photon
  vector<string> phoLabels {
    "All",
      "pt",
      "ScEta",
      "veto",
      "mvaid"
      };
  HtoJPsiGamma::showEfficiency("photonCutFlow", phoLabels, "Photon Selection","Photons");
}

void HtoJPsiGamma::eventEfficiency(){
  // outfile->cd();
  //outfile->cd("HtoJPsiGamma");
  std::vector<std::string> evlabels {
    "Events processed",
      //      "Gen level filter : Minimum 2 gen Muons & Minimum 1 gen photon",
      //"Gen level filter after detector cuts, eta =< 2.5",
      "Events with > 0 good vertex",
      "Events passing trigger",
      "Atleast 2 muons after muon selection",
      //"Atleast 1 photon after photon selection",     
      "Relative Isolation of leading mu",
      "leading mu pt",
      "subleading mu pt",
      "Atleast 1 photon after photon selection",
      "deltaR between mu1 & pho",
      "deltaR between mu2 & pho",
      "deltaR between dimu & pho",
      "deltaphi between dimu & pho",
      "JPsi mass cut",
      "H mass cut 100-150",
      "ratio cut for pt of dimu",
      "ratio cut for Et of pho",
      "scale factor"
      };
  
  HtoJPsiGamma::showEfficiency("evtCutFlow",evlabels, "Event Selection","Eventflow");

  double lumiFac = 1.0;
  if (useLumiWt_)
    lumiFac = lumiWt(evtWeightSum_, 200000);
  cout << endl
       << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
       << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
       << "      useLumiWt_: " << setw(10) << setprecision(5) <<lumiFac
       << endl;
  
  //  HtoJPsiGamma::scaleHistogram("evtCutFlowWt", lumiFac);

  TH1 *h = HtoJPsiGamma::getHist1D("evtCutFlowWt");
  if (h != nullptr) {
    int nbins = h->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
      double cont = h->GetBinContent(i) * lumiFac;
      double err  = h->GetBinError(i) * lumiFac;
      h->SetBinContent(i, cont);
      h->SetBinError(i, err);
    }
  }

  HtoJPsiGamma::showEfficiency("evtCutFlowWt", evlabels, "Event Selection (Weighted)", "EventFlow");

      std::vector<std::string> catlabels {
      "total final yield",
	"cat1",
	"cat2",
	"cat3"
	//"cat4"
      };

        HtoJPsiGamma::showEfficiency("nPhoCat",catlabels, "category of events", "category");

  double lumiFac_cat = 1.0;
  if (useLumiWt_)
    lumiFac_cat = lumiWt(evtWeightSum_, 200000);
  cout << endl
       << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
       << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac_cat
       << endl;
  outfile->cd("HtoJPsiGamma");
  std::cout<<"the Lumi fac is:"<<lumiFac_cat<<std::endl;
  TH1D* lumi_fac=new TH1D("lumi_fac", "total lumi wt", 1, -0.5, 0.5);
  lumi_fac->SetBinContent(1,lumiFac_cat);
  lumi_fac->Write();

  TH1 *hcat = HtoJPsiGamma::getHist1D("nPhoCatWt");
  if (hcat != nullptr) {
    int nbins_cat = hcat->GetNbinsX();
    for (int j = 1; j <= nbins_cat; ++j) {
      double cont_cat = hcat->GetBinContent(j) * lumiFac_cat;
      double err_cat  = hcat->GetBinError(j) * lumiFac_cat;
      hcat->SetBinContent(j, cont_cat);
      hcat->SetBinError(j, err_cat);
    }
    hcat->Write();
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt",catlabels, "category of weighted events", "Category");

  TH1 *hcat1 = HtoJPsiGamma::getHist1D("nPhoCatWt_triggP");
  if (hcat1 != nullptr) {
    int nbins_cat1 = hcat1->GetNbinsX();
    for (int j = 1; j <= nbins_cat1; ++j) {
      double cont_cat1 = hcat1->GetBinContent(j) * lumiFac_cat;
      double err_cat1  = hcat1->GetBinError(j) * lumiFac_cat;
      hcat1->SetBinContent(j, cont_cat1);
      hcat1->SetBinError(j, err_cat1);
    }
    hcat1->Write();                                                                                        
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt_triggP",catlabels, "category of weighted events triggP", "Category");

  TH1 *hcat2 = HtoJPsiGamma::getHist1D("nPhoCatWt_triggN");
  if (hcat2 != nullptr) {
    int nbins_cat2 = hcat2->GetNbinsX();
    for (int j = 1; j <= nbins_cat2; ++j) {
      double cont_cat2 = hcat2->GetBinContent(j) * lumiFac_cat;
      double err_cat2  = hcat2->GetBinError(j) * lumiFac_cat;
      hcat2->SetBinContent(j, cont_cat2);
      hcat2->SetBinError(j, err_cat2);
    }
    hcat2->Write();                                                                                         
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt_triggN",catlabels, "category of weighted events triggN", "Category");

  TH1 *hcat3 = HtoJPsiGamma::getHist1D("nPhoCatWt_muP");
  if (hcat3 != nullptr) {
    int nbins_cat3 = hcat3->GetNbinsX();
    for (int j = 1; j <= nbins_cat3; ++j) {
      double cont_cat3 = hcat3->GetBinContent(j) * lumiFac_cat;
      double err_cat3  = hcat3->GetBinError(j) * lumiFac_cat;
      hcat3->SetBinContent(j, cont_cat3);
      hcat3->SetBinError(j, err_cat3);
    }
    hcat3->Write();
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt_muP",catlabels, "category of weighted events muP", "Category");

  TH1 *hcat4 = HtoJPsiGamma::getHist1D("nPhoCatWt_muN");
  if (hcat4 != nullptr) {
    int nbins_cat4 = hcat4->GetNbinsX();
    for (int j = 1; j <= nbins_cat4; ++j) {
      double cont_cat4 = hcat4->GetBinContent(j) * lumiFac_cat;
      double err_cat4  = hcat4->GetBinError(j) * lumiFac_cat;
      hcat4->SetBinContent(j, cont_cat4);
      hcat4->SetBinError(j, err_cat4);
    }
    hcat4->Write();
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt_muN",catlabels, "category of weighted events muN", "Category");



  TH1 *hcat5 = HtoJPsiGamma::getHist1D("nPhoCatWt_phoP");
  if (hcat5 != nullptr) {
    int nbins_cat5 = hcat5->GetNbinsX();
    for (int j = 1; j <= nbins_cat5; ++j) {
      double cont_cat5 = hcat5->GetBinContent(j) * lumiFac_cat;
      double err_cat5 = hcat5->GetBinError(j) * lumiFac_cat;
      hcat5->SetBinContent(j, cont_cat5);
      hcat5->SetBinError(j, err_cat5);
    }
    hcat5->Write();
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt_phoP",catlabels, "category of weighted events phoP", "Category");


  TH1 *hcat6 = HtoJPsiGamma::getHist1D("nPhoCatWt_phoN");
  if (hcat6 != nullptr) {
    int nbins_cat6 = hcat6->GetNbinsX();
    for (int j = 1; j <= nbins_cat6; ++j) {
      double cont_cat6 = hcat6->GetBinContent(j) * lumiFac_cat;
      double err_cat6 = hcat6->GetBinError(j) * lumiFac_cat;
      hcat6->SetBinContent(j, cont_cat6);
      hcat6->SetBinError(j, err_cat6);
    }
    hcat6->Write();
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt_phoN",catlabels, "category of weighted events phoN", "Category");

  TH1 *hcat7 = HtoJPsiGamma::getHist1D("nPhoCatWt_PuUp");
  if (hcat7 != nullptr) {
    int nbins_cat7 = hcat7->GetNbinsX();
    for (int j = 1; j <= nbins_cat7; ++j) {
      double cont_cat7 = hcat7->GetBinContent(j) * lumiFac_cat;
      double err_cat7 = hcat7->GetBinError(j) * lumiFac_cat;
      hcat7->SetBinContent(j, cont_cat7);
      hcat7->SetBinError(j, err_cat7);
    }
    hcat7->Write();
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt_PuUp",catlabels, "category of weighted events PuUp", "Category");

  TH1 *hcat8 = HtoJPsiGamma::getHist1D("nPhoCatWt_PuDown");
  if (hcat8 != nullptr) {
    int nbins_cat8 = hcat8->GetNbinsX();
    for (int j = 1; j <= nbins_cat8; ++j) {
      double cont_cat8 = hcat8->GetBinContent(j) * lumiFac_cat;
      double err_cat8 = hcat8->GetBinError(j) * lumiFac_cat;
      hcat8->SetBinContent(j, cont_cat8);
      hcat8->SetBinError(j, err_cat8);
    }
    hcat8->Write();
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt_PuDown",catlabels, "category of weighted events PuDown", "Category");


  TH1 *hcat9 = HtoJPsiGamma::getHist1D("nPhoCatWt_eleVetoP");
  if (hcat9 != nullptr) {
    int nbins_cat9 = hcat9->GetNbinsX();
    for (int j = 1; j <= nbins_cat9; ++j) {
      double cont_cat9 = hcat9->GetBinContent(j) * lumiFac_cat;
      double err_cat9 = hcat9->GetBinError(j) * lumiFac_cat;
      hcat9->SetBinContent(j, cont_cat9);
      hcat9->SetBinError(j, err_cat9);
    }
    hcat9->Write();
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt_eleVetoP",catlabels, "category of weighted events eleVetoP", "Category");


  TH1 *hcat10 = HtoJPsiGamma::getHist1D("nPhoCatWt_eleVetoN");
  if (hcat10 != nullptr) {
    int nbins_cat10 = hcat10->GetNbinsX();
    for (int j = 1; j <= nbins_cat10; ++j) {
      double cont_cat10 = hcat10->GetBinContent(j) * lumiFac_cat;
      double err_cat10 = hcat10->GetBinError(j) * lumiFac_cat;
      hcat10->SetBinContent(j, cont_cat10);
      hcat10->SetBinError(j, err_cat10);
    }
    hcat10->Write();
  }
  HtoJPsiGamma::showEfficiency("nPhoCatWt_eleVetoN",catlabels, "category of weighted events eleVetoN", "Category");

}
      
void HtoJPsiGamma::showCategoryYield(){
  // just for test, please ignore
}

// Helicity angle between J/Psi and mu+ 
float HtoJPsiGamma::costheta1(TLorentzVector particle_orig, TLorentzVector parent_orig){
  TLorentzVector particle = particle_orig; // make a dummy
  TLorentzVector parent = parent_orig;
  TVector3 boosttoparent = -(parent.BoostVector());
  particle.Boost(boosttoparent);
  parent.Boost(boosttoparent);
 
  TVector3 particle3 = particle.Vect();
  TVector3 parent3 = parent.Vect();
  Float_t numerator = particle3.Dot(parent3);
  Float_t denominator = (particle3.Mag())*(parent3.Mag());
  Float_t temp = numerator/denominator;

  return temp;
}

float HtoJPsiGamma::costhetastar(TLorentzVector particle_orig, TLorentzVector parent_orig){
  TLorentzVector particle = particle_orig; // make a dummy
  TLorentzVector parent = parent_orig;
  TVector3 boosttoparent = -(parent.BoostVector());
  particle.Boost(boosttoparent);
  parent.Boost(boosttoparent);
  //  TVector3 particle3 = particle.Vect();
  //  TVector3 parent3 = parent.Vect();
  Float_t temp = particle.CosTheta();
  return temp;
}

void HtoJPsiGamma::bookHistograms() {
  outfile->cd();
  outfile->mkdir("ObjectSelection");
  outfile->cd("ObjectSelection");
  new TH1D("muCutFlow", "Muon Cut Flow", 3, -0.5, 2.5);
  new TH1D("photonCutFlow", "Photon Cut Flow", 5, -0.5, 4.5);
  new TH1D("photonCutFlowWt", "Photon Cut Flow Wt", 5, -0.5, 4.5);

  
  new TH1F("muPt", "Muon p_{T}", 200, 0., 100.);//mu loose id
  new TH1F("muEta", "Muon #eta", 200, -4., 4.);
  new TH1F("muDxy", "Muon d_{xy}", 50, -10., 10.);
  new TH1F("muDz", "Muon d_z", 50, -10., 10.);
  new TH1F("muSIP3d", "Muon sip3d", 50, -10., 10.);
  new TH1F("muSIP3dincm", "Muon sip3d in cm", 50, -10., 10.);
  new TH1F("muLooseid", "Muon looseid", 2, 0., 1.);
  new TH1F("muIspfcand", "Muon isPFcand", 2, 0., 1.);//mutightid
  new TH1F("muHighptid", "Muon isHighptid", 2, 0., 1.);
  new TH1F("muTightid", "Muon tightId", 2, 0., 1.);

  new TH1F("phoMvaidwp90", "Photon is Mvaidwp90", 2, 0., 1.);//photon
  //  new TH1F("phoPt", "Photon p_{T}", 400, 0., 200.);
  new TH1F("phoPt", "Photon p_{T} for leading photons passing MVA", 400, 0., 200.);
  new TH1F("phoPt_MVA", "Photon p_{T} for Highest Pt photon passing MVA", 400, 0., 200.);
  new TH1F("Photon_isScEtaEB", "Photon isScEtaEB", 2, 0., 1.);
  new TH1F("Photon_isScEtaEE", "Photon isScEtaEE", 2, 0., 1.);
  new TH1D("nPhoScEta", "#Photons in ScEta regions", 2, -0.5, 1.5);
  new TH1D("nPhoScEta_", "#Photons in ScEta regions after veto cut", 2, -0.5, 1.5);
  new TH1D("nPhoScEta_W", "#Photons in ScEta regions after veto cut W", 2, -0.5, 1.5);

  new TH1D("phoR9", "Photon R9", 40, 0.4, 1.0);
  new TH1D("phoR9_ScEB", "Photon R9 for ScEta in EB", 40, 0.4, 1.0);
  new TH1D("phoR9_ScEE", "Photon R9 for ScEta in EE", 40, 0.4, 1.0);

  outfile->cd();
  outfile->mkdir("HtoJPsiGamma");
  outfile->cd("HtoJPsiGamma");
  new TH1F("generatorweight", "H MC generator weight", 100, 8.5e16, 9e+16);
  new TH1F("generatorweight_S", "H MC generator weight scaled by gen event sum",100, 0.95e-6, 1.5e-6);
  new TH1F("puweight", "PU weight", 200, 0., 2.);
  new TH2F("pudist","Pu weight distribution", 80, 0., 80., 100., 0., 2. );
  /*
  new TH1F("genweight", "generator weight", 3, -1.5, 1.5);
  new TH1F("genMu1Pt","leading gen mu p_T ", 50, 0.,100.);
  new TH1F("genMu2Pt","subleading gen mu p_T ", 50,0., 100.);
  new TH1F("genPhoPt"," gen photon p_T ", 40, 0., 200.);
  new TH1F("genMu1Eta","leading gen mu eta ", 50, -2.5, 2.5);
  new TH1F("genMu2Eta","subleading gen mu eta ", 50, -2.5 , 2.5);
  new TH1F("genPhoEta"," gen photon eta ", 50, -2.5, 2.5);
  new TH1F("gendRMu1Mu2","dR between gen mu1 & mu2", 50, 0., 1.);
  new TH1F("gendRdiMuPho","dR between gen diMu & Pho", 50, 0., 5.);
  new TH1F("gendRMu1Pho","dR between gen diMu1 & Pho", 50, 0., 5.);
  new TH1F("gendRMu2Pho","dR between gen diMu2 & Pho", 50, 0., 5.);
  */
  //Preselection
  new TH1F("nPV", "No. of reconstructed primary vertices", 15, -0.5, 14.5);
  new TH1F("nPVgood", "No. of good reconstructed primary vertices, passing all the criteria", 15, -0.5, 14.5);

  //Helicity
  new TH1F("costheta1","Helicity angle between JPsi and mu",80, -1.0, 1.0);
  new TH1F("costheta1_rewt","Helicity angle between JPsi and mu (reweighted)",80, -1.0, 1.0);
  new TH1F("costhetas","Helicity angle between JPsi and beam",80, -1.0, 1.0);
  new TH1F("costhetas_rewt","Helicity angle between JPsi and beam (reweighted)",80, -1.0, 1.0);

  new TH1F("nLoosemuons", "No. of loose muons", 10, -0.5, 9.5);
  new TH1F("nTightmuons", "No. of tight muons", 10, -0.5, 9.5);
  new TH1F("nPhotons", "No. of photons", 5, -0.5, 4.5);

  new TH1F("muPfRelIso03all", "muPfRelIso03all", 100, 0., 1.);
  new TH1F("muPfRelIso03chg", "muPfRelIso03chg", 100, 0., 1.);
  new TH1F("muPfRelIso04all", "muPfRelIso04all", 100, 0., 1.);
  new TH1F("muPfRelIso03allByMuSubleadpt", "muPfRelIso03all", 80, 0., 2.);
  new TH1F("muLeadpt", "Leading Mu pt", 20, 0., 100.);
  new TH1F("muSubleadpt", "Subleading Mu pt", 20, 0., 100.);
  
  new TH1F("dRmu1pho", "dR between mu1 & photon", 100, 0., 5.);
  new TH1F("dRmu2pho", "dR between mu2 & photon", 100, 0., 5.);
  new TH1F("dRdimupho", "dR between diMuon & photon", 120, 0., 8.);
  new TH1F("dPhidimupho", "dPhi between diMuon & photon", 120, -4., 4.);
  new TH1F("ratiocutdimu", "ratio cut for diMuon", 100, 0., 1.);
  new TH1F("ratiocutpho", "ratio cut for photon", 100, 0., 1.);

  new TH1F("mvaoutput_bdtg_tlumi_wt_hasGen","mvaoutput_bdtg_tlumi_wt_hasGen,",40, -1., 1.);
  new TH1F("mvaoutput_bdtg_pu_gen_wt","mvaoutput_bdtg_pu_gen_wt",40, -1., 1.);
  new TH1F("mvaoutput_bdtg","mvaoutput_bdtg",40, -1., 1.);

  new TH1F("jpsimass_bmasscut","jpsimass before mass cut", 60, 2.8, 3.4);
  new TH1F("jpsimass_amasscut","jpsimass after mass cut 120-130", 60, 2.8, 3.4);
  new TH1F("jpsimass_afinalcut","jpsimass after final cut", 60, 2.8, 3.4);
  new TH1F("zmass_bmasscut","zmass before Z mass cut 120_130", 64, 40., 200.);
  new TH1F("zmass_bratiocut1","zmass after all cut", 64, 40., 200.);
  new TH1F("zmass_bratiocut2","zmass after all cut", 64, 40., 200.);
  new TH1F("zmass_bscorecut","zmass after all cut", 64, 40., 200.);


  new TH1F("zmass_amasscut","zmass after Z mass cut", 64, 40., 200.);
  new TH1F("zmass_afinalcut","zmass after final cut", 64, 40., 200.);

  new TH1F("vtxdist_x","DiMuon vertex position (X)", 40, -0.2, 0.2);
  new TH1F("vtxdist_y","DiMuon vertex position (Y)", 40, -0.2, 0.2);
  new TH1F("vtxdist_z","DiMuon vertex position (Z)", 40, -0.2, 0.2);
  
  new TH1D("evtCutFlow", "Event CutFlow", 17, -0.5, 16.5); // event cutflow
  new TH1D("evtCutFlowWt", "Event CutFlow(Weighted)", 17, -0.5, 16.5); // event cutflow weighted
  
  // Event category for Z signal
 
  new TH1D("nPhoCat", "Event Category", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt", "Event Category weighted", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt_triggP", "Event Category weighted triggP", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt_triggN", "Event Category weighted triggN", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt_muP", "Event Category weighted muP", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt_muN", "Event Category weighted muN", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt_phoP", "Event Category weighted phoP", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt_phoN", "Event Category weighted phoN", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt_eleVetoP", "Event Category weighted eleVetoP", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt_eleVetoN", "Event Category weighted eleVetoN", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt_PuUp", "Event Category weighted PuUp", 4, -0.5, 3.5);
  new TH1D("nPhoCatWt_PuDown", "Event Category weighted PuDown", 4, -0.5, 3.5);
  new TH2F("scoreVsInv","bdt score vs INV mass",64,40.,200.,40,-1,1);
  /*
  new TH1F("patdRmu1mu2_bHLT","dR between pat mu1 & mu2 before HLT", 50, 0., 1.);
  new TH1F("patMmu1mu2_bHLT","Invariant Mass of pat mu1 & mu2 before HLT", 30, 2.8, 3.4);
  new TH1F("patPTmu1mu2_bHLT","Pt of pat mu1 & mu2 before HLT", 45, 10., 100.);
  new TH1F("patdRmu1mu2_aHLT","dR between pat mu1 & mu2 after HLT", 50, 0., 1.);
  new TH1F("patMmu1mu2_aHLT","Invariant Mass of pat mu1 & mu2 after HLT", 30, 2.8, 3.4);
  new TH1F("patPTmu1mu2_aHLT","Pt of pat mu1 & mu2 after HLT", 45, 10., 100.);
  */
  // categorywise final plot 
  new TH1F("muLeadpt_cat1", "Leading Mu pt for cat1 ", 32, 0., 160.);
  new TH1F("muSubleadpt_cat1", "Subleading Mu pt for cat1", 32, 0., 160.);
  new TH1F("phoPt_cat1", "Leading Mu pt for cat1 ", 32, 0., 160.);
  new TH1F("muLeadeta_cat1","leading Mu eta for cat1 ", 25, -2.5, 2.5);
  new TH1F("muSubleadeta_cat1","Subleading Mu eta for cat1 ", 25, -2.5, 2.5);
  new TH1F("phoEta_cat1","leading photon eta for cat1 ", 25, -2.5, 2.5);
  new TH1F("dRmu1pho_cat1", "dR between mu1 & photon for cat1", 50, 0., 5.);
  new TH1F("dRmu2pho_cat1", "dR between mu2 & photon for cat1", 50, 0., 5.);
  new TH1F("dRdimupho_cat1", "dR between diMuon & photon for cat1", 80, 0., 8.);
  new TH1F("dPhidimupho_cat1", "dPhi between diMuon & photon for cat1", 80, -4., 4.);
  new TH1F("dRmu1mu2_cat1", "dR between mu1 & mu2 for cat1", 100, 0., 1.);
  new TH1F("jpsimass_cat1","jpsimass for cat1", 40, 2.9, 3.3);  
  new TH1F("diMuPt_cat1","DiMuon Pt for cat1", 32, 0. , 160.);  
  new TH1F("finalZPt_cat1","Pt for final dimuon+gamma for cat1", 25, 0., 100.);  
  new TH1F("zmass_cat1","Zmass for cat1", 64, 40., 200.);  

  new TH1F("muLeadpt_cat2", "Leading Mu pt for cat2 ", 32, 0., 160.);
  new TH1F("muSubleadpt_cat2", "Subleading Mu pt for cat2", 32, 0., 160.);
  new TH1F("phoPt_cat2", "Leading Mu pt for cat2 ", 32, 0., 160.);
  new TH1F("muLeadeta_cat2","leading Mu eta for cat2 ", 25, -2.5, 2.5);
  new TH1F("muSubleadeta_cat2","Subleading Mu eta for cat2 ", 25, -2.5, 2.5);
  new TH1F("phoEta_cat2","leading photon eta for cat2 ", 25, -2.5, 2.5);
  new TH1F("dRmu1pho_cat2", "dR between mu1 & photon for cat2", 50, 0., 5.);
  new TH1F("dRmu2pho_cat2", "dR between mu2 & photon for cat2", 50, 0., 5.);
  new TH1F("dRdimupho_cat2", "dR between diMuon & photon for cat2", 80, 0., 8.);
  new TH1F("dPhidimupho_cat2", "dPhi between diMuon & photon for cat2", 80, -4., 4.);
  new TH1F("dRmu1mu2_cat2", "dR between mu1 & mu2 for cat2", 100, 0., 1.);
  new TH1F("jpsimass_cat2","jpsimass for cat2", 40, 2.9, 3.3);
  new TH1F("diMuPt_cat2","DiMuon Pt for cat2", 32, 0. , 160.);
  new TH1F("finalZPt_cat2","Pt for final dimuon+gamma for cat2", 25, 0., 100.);
  new TH1F("zmass_cat2","Zmass for cat2", 64, 40., 200.);

  new TH1F("muLeadpt_cat3", "Leading Mu pt for cat3 ", 32, 0., 160.);
  new TH1F("muSubleadpt_cat3", "Subleading Mu pt for cat3", 32, 0., 160.);
  new TH1F("phoPt_cat3", "Leading Mu pt for cat3 ", 32, 0., 160.);
  new TH1F("muLeadeta_cat3","leading Mu eta for cat3 ", 25, -2.5, 2.5);
  new TH1F("muSubleadeta_cat3","Subleading Mu eta for cat3 ", 25, -2.5, 2.5);
  new TH1F("phoEta_cat3","leading photon eta for cat3 ", 25, -2.5, 2.5);
  new TH1F("dRmu1pho_cat3", "dR between mu1 & photon for cat3", 50, 0., 5.);
  new TH1F("dRmu2pho_cat3", "dR between mu2 & photon for cat3", 50, 0., 5.);
  new TH1F("dRdimupho_cat3", "dR between diMuon & photon for cat3", 80, 0., 8.);
  new TH1F("dPhidimupho_cat3", "dPhi between diMuon & photon for cat3", 80, -4., 4.);
  new TH1F("dRmu1mu2_cat3", "dR between mu1 & mu2 for cat3", 100, 0., 1.);
  new TH1F("jpsimass_cat3","jpsimass for cat3", 40, 2.9, 3.3);
  new TH1F("diMuPt_cat3","DiMuon Pt for cat3", 32, 0. , 160.);
  new TH1F("finalZPt_cat3","Pt for final dimuon+gamma for cat3", 25, 0., 100.);
  new TH1F("zmass_cat3","Zmass for cat3", 64, 40., 200.);

  //New Histograms.....                                                                                                                                  
  new TH1F("Mu1_jetIDX","muon1 has jet index", 2, -0.5, 1.5);
  new TH1F("Mu2_jetIDX","muon2 has jet index", 2, -0.5, 1.5);
  new TH1F("Pho_jetIDX","photon has jet index ", 2, -0.5, 1.5);
  new TH1F("Pho_eleIDX","photon has electron index ", 2, -0.5, 1.5);
  new TH1F("Mu1_jetIDX_Wt","muon1 has jet index Wt", 2, -0.5, 1.5);
  new TH1F("Mu2_jetIDX_Wt","muon2 has jet index Wt", 2, -0.5, 1.5);
  new TH1F("Pho_jetIDX_Wt","photon has jet index Wt", 2, -0.5, 1.5);
  new TH1F("Pho_eleIDX_Wt","photon has electron index Wt", 2, -0.5, 1.5);


  // To select which cut to choose
  //  new TH1F("dRMu1Cut", "# of Events for different dR(mu1, #gamma)", )
    
  bookedHistograms = true;
  
}

void HtoJPsiGamma::saveHistograms() {
  cout << "About to close the histogram file" << endl;
  if (bookedHistograms) {
    //    outfile->cd();
    //outfile->cd("HtoJPsiGamma");
    outfile->Write();  
    if (skimObj_ != nullptr) skimObj_->close();
    //    if (skimObj_ != nullptr) skimObj_->Write();
    //    outfile->Close();
  }
  cout << "histogram file not closed" << endl;

}
// ***********************************************


bool HtoJPsiGamma::readJob(const std::string& jobFile, int& nFiles){
  
  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), std::ios::in);
  if (!fin) {
    cerr << "==> Input File: <<" << jobFile << ">> could not be opened " << endl;
    return false;
  }
  // eventFileList_.clear();

  static constexpr int BUF_SIZE = 256;
  char buf[BUF_SIZE];
  while (fin.getline(buf, BUF_SIZE, '\n')) {
    // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;

    // enable '#' and '//' style comments                           
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    // Split the line into words

    vector<string> tokens;
    HtoJPsiGamma::tokenize(line, tokens);
    int vsize = tokens.size();
    assert(vsize > 1);

    const string& key   = tokens.at(0);
    const string& value = tokens.at(1);

    if (key == "dataType") {
      string vtmp(value);
      std::transform(vtmp.begin(), vtmp.end(), vtmp.begin(), ::toupper);
      vector<string> dt;
      tokenize(vtmp, dt, "#");
      if (dt.size()) {
        isMC_ = (dt.at(0) == "MC") ? true : false;
        if (isMC_ && dt.size() > 1) {
          isSignal_ = (dt.at(1) == "SIGNAL") ? true : false;
        }
      }
    }
    else if(key== "jsonFile") {
      jsonFile_ = value;
      setJson(jsonFile_);
    }
    else if (key == "useLumiWt")
      useLumiWt_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "intLumi")
      intLumi_ = std::stof(value.c_str());
    else if (key == "xsec")
      xsec_ = std::stof(value.c_str());
    else if (key == "maxEvent")
      maxEvent_ = std::stoi(value.c_str());
    else if (key == "histFile")
      histFile_ = value;
    else if (key == "hlt")
      hlt_ = std::stoi(value.c_str());
    else if (key == "year")
      year_ = std::stoi(value.c_str());
    else if (key == "isPre") {
      string vtmp(value);
      isPre_ = (vtmp == "True") ? true:false;
    }

    else if(key == "find_SU") {
      find_SU = std::stoi(value.c_str()) > 0 ? true:false;
    }
    else if (key == "usePu_wt") {
      ispu_wt = std::stoi(value.c_str()) > 0 ? true:false;
    }
    else if (key == "useGen_wt") {
      isgen_wt= std::stoi(value.c_str()) > 0 ? true:false;
    }
    else if (key == "usePU_gen_wt") {
      ispu_gen_wt= std::stoi(value.c_str()) > 0  ? true:false;
    }
    else if (key == "useTlumi_wt") {
      istlumi_wt= std::stoi(value.c_str()) > 0 ? true:false;
    }
    else if (key == "useTlumi_wt_hasGen") {
      istlumi_wt_hasGen= std::stoi(value.c_str()) > 0 ? true:false;
    }
    else if (key == "createMVATree")
      createMVATree_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "mvaSkimFile")
      mvaSkimFile_ = value;
    else if (key == "readMVA")
      readMVA_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "mvaXMLFile")
      mvaXMLFile_ = value;
    else if (key == "mvAlgo")
      mvAlgo_ = value;
    else if (key == "inputFile")
      HtoJPsiGamma::buildList(tokens, fileList_);
    }
  
  
  // Close the file                                                  
  fin.close();
  
  //if (!isMC_) usePUWt_ = false;
  //if (!isSignal_) readGenInfo_ = false;

  // Build the chain of root files
  for (const auto& fname: fileList_) {
    cout << ">>> INFO. Adding input file " << fname << " to TChain " << endl;
    ++nFiles;
    int nevt = setInputFile(fname);
    if (maxEvent_ > 0 && nevt >= maxEvent_) break;
  }

  if (!nFiles) {
    cerr << ">>> WARN. Input Root file list is empty! exiting ..." << endl;
    return false;
  }

  return true;



}

// ************************************

void HtoJPsiGamma::Loop()
{
  Double_t geneventsumw=0.;
  if (chain1 == 0) return;
  int nentries1 = static_cast<int>(chain1->GetEntries());
  std::cout << "Total number of entries in Run tree: " << nentries1 << std::endl;
  Long64_t nbytes1 = 0, nb1 = 0;
  for (Long64_t jentry=0; jentry<nentries1;jentry++) {
    genEventSumw=0;
    genEventSumw_=0;
    Long64_t ientry = chain1->LoadTree(jentry);
    if (ientry < 0) break;
    nb1 = chain1->GetEntry(jentry);   nbytes1 += nb1;
    if(isMC()) {
      //std::cout<<jentry<<","<<genEventSumw<<","<<genEventSumw_<<std::endl;                                           
      if(genEventSumw!=0.) {
        geneventsumw += genEventSumw;
      }
      else   {
        geneventsumw += genEventSumw_;
      }
    }
  }

  Float_t zm1 =0,zm2=0,zm3=0, mu1pt=0, mu2pt=0, phopt=0, muiso=0, drmu1pho=0, drmu2pho=0, drjpsipho=0, dphijpsipho=0, jpsim=0, Hm=0, ptratio=0, etratio=0, finalm=0, pu_wt=0, gen_wt=0, pu_gen_wt=0, tlumi_wt=0, tlumi_wt_hasGen=0, pu_wt1=0, gen_wt1=0, pu_gen_wt1=0, tlumi_wt1=0, tlumi_wt_hasGen1=0, pu_wt2=0, gen_wt2=0, pu_gen_wt2=0, tlumi_wt2=0, tlumi_wt_hasGen2=0, pu_wt3=0, gen_wt3=0, pu_gen_wt3=0, tlumi_wt3=0, tlumi_wt_hasGen3=0;
  
  Float_t phoeta=0,phophi=0, mu1_sip3d=0,mu2_sip3d=0,mu1eta=0,mu1phi=0,mu2eta=0,mu2phi=0,mu1iso=0,mu2iso=0,mod_dphijpsipho=0,dphimu1pho=0,mod_dphimu1pho=0,dphimu2pho=0,mod_dphimu2pho=0,dphimu1mu2=0, mod_dphimu1mu2=0, drmu1mu2=0, jpsipt=0,Hmpt=0,photon_r9=0, jpsiSSpt=0, HmSSpt=0;
  
  
  Float_t pu_gen_wt_triggP = 0, pu_gen_wt_triggN =0, pu_gen_wt_muP =0, pu_gen_wt_muN =0, pu_gen_wt_phoP = 0, pu_gen_wt_phoN = 0, pu_gen_wt_eleVetoP = 0, pu_gen_wt_eleVetoN = 0;/*, pu_gen_wt_triggP1 = 0, pu_gen_wt_triggN1 =0, pu_gen_wt_muP1 =0, pu_gen_wt_muN1 =0, pu_gen_wt_phoP1 = 0, pu_gen_wt_phoN1 = 0, pu_gen_wt_eleVetoP1 = 0, pu_gen_wt_eleVetoN1 = 0, pu_gen_wt_triggP2 = 0, pu_gen_wt_triggN2 = 0, pu_gen_wt_muP2 =0, pu_gen_wt_muN2 =0, pu_gen_wt_phoP2 = 0, pu_gen_wt_phoN2 = 0, pu_gen_wt_eleVetoP2 = 0, pu_gen_wt_eleVetoN2 = 0, pu_gen_wt_triggP3 = 0, pu_gen_wt_triggN3 =0, pu_gen_wt_muP3 =0, pu_gen_wt_muN3 =0, pu_gen_wt_phoP3 = 0, pu_gen_wt_phoN3 = 0, pu_gen_wt_eleVetoP3 = 0, pu_gen_wt_eleVetoN3 = 0 ;*/
  
    /*
  Float_t pu_wt_triggP = 0, pu_wt_triggN =0, pu_wt_muP =0, pu_wt_muN =0, pu_wt_phoP = 0, pu_wt_phoN = 0, pu_wt_eleVetoP = 0, pu_wt_eleVetoN = 0, pu_wt_triggP1 = 0, pu_wt_triggN1 =0, pu_wt_muP1 =0, pu_wt_muN1 =0, pu_wt_phoP1 = 0, pu_wt_phoN1 = 0, pu_wt_eleVetoP1 = 0, pu_wt_eleVetoN1 = 0, pu_wt_triggP2 = 0, pu_wt_triggN2 = 0, pu_wt_muP2 =0, pu_wt_muN2 =0, pu_wt_phoP2 = 0, pu_wt_phoN2 = 0, pu_wt_eleVetoP2 = 0, pu_wt_eleVetoN2 = 0, pu_wt_triggP3 = 0, pu_wt_triggN3 =0, pu_wt_muP3 =0, pu_wt_muN3 =0, pu_wt_phoP3 = 0, pu_wt_phoN3= 0, pu_wt_eleVetoP3 = 0, pu_wt_eleVetoN3 = 0 ;
  */
  Float_t pu_wt_=0, pu_gen_wt_=0, pu_wt_1=0, pu_gen_wt_1=0, pu_wt_2=0, pu_gen_wt_2=0,pu_wt_3=0, pu_gen_wt_3=0, LSigma =0 , LSigma1=0, LSigma2=0, LSigma3=0;
 
  Float_t pu_wtUp=0, pu_wtDown=0, pu_gen_wtUp=0, pu_gen_wtDown=0;/* pu_wtUp1=0, pu_wtDown1=0, pu_gen_wtUp1=0, pu_gen_wtDown1=0, pu_wtUp2=0, pu_wtDown2=0, pu_gen_wtUp2=0, pu_gen_wtDown2=0, pu_wtUp3=0, pu_wtDown3=0, pu_gen_wtUp3=0, pu_gen_wtDown3=0;*/

  
  struct ComparePt {
    bool operator()(const TLorentzVector& vec1, const TLorentzVector& vec2) const {
      return vec1.Pt() > vec2.Pt();
    }
  };

  if(chain == 0) return;
  //Long64_t nentries = chain->GetEntriesFast();                                                                                                                                                     
  int nentries = static_cast<int>(chain->GetEntries());
  std::cout << "Total number of entries: " << nentries<< std::endl;
  if (maxEvent_ > 0) nentries = std::min(nentries, maxEvent_);
  //std::cout << "max event: " << maxEvt_ << std::endl;                                                                                                                                         
  std::cout << "Total number of entries running for: " << nentries<< std::endl;
  outfile = TFile::Open(histFile_.c_str(), "RECREATE");
  
  TTree *T = new TTree("T","test z mass tree");
  TTree *T1 = new TTree("T1","cat1 mass tree");
  TTree *T2 = new TTree("T2","cat2 mass tree");
  TTree *T3 = new TTree("T3","cat3 mass tree");
  T1->Branch("zm1",&zm1,"zm1/F");
  T1->Branch("pu_wt1",& pu_wt1," pu_wt1/F");
  T1->Branch("pu_wt_1",& pu_wt_1," pu_wt_1/F");
  T1->Branch("pu_gen_wt1",&pu_gen_wt1,"pu_gen_wt1/F");
  T1->Branch("pu_gen_wt_1",&pu_gen_wt_1,"pu_gen_wt_1/F");
  T1->Branch("gen_wt1",&gen_wt1,"gen_wt1/F");
  T1->Branch("tlumi_wt1",&tlumi_wt1,"tlumi_wt1/F");
  T1->Branch("tlumi_wt_hasGen1",&tlumi_wt_hasGen1,"tlumi_wt_hasGen1/F");
  T1->Branch("LSigma1",&LSigma1,"LSigma1/F");


  T2->Branch("zm2",&zm2,"zm2/F");
  T2->Branch("pu_wt2",& pu_wt2," pu_wt2/F");
  T2->Branch("pu_wt_2",& pu_wt_2," pu_wt_2/F");
  T2->Branch("pu_gen_wt2",&pu_gen_wt2,"pu_gen_wt2/F");
  T2->Branch("pu_gen_wt_2",&pu_gen_wt_2,"pu_gen_wt_2/F");
  T2->Branch("gen_wt2",&gen_wt2,"gen_wt2/F");
  T2->Branch("tlumi_wt2",&tlumi_wt2,"tlumi_wt2/F");
  T2->Branch("tlumi_wt_hasGen2",&tlumi_wt_hasGen2,"tlumi_wt_hasGen2/F");
  T2->Branch("LSigma2",&LSigma2,"LSigma2/F");

  T3->Branch("zm3",&zm3,"zm3/F");
  T3->Branch("pu_wt3",& pu_wt3," pu_wt3/F");
  T3->Branch("pu_wt_3",& pu_wt_3," pu_wt_3/F");
  T3->Branch("pu_gen_wt3",&pu_gen_wt3,"pu_gen_wt3/F");
  T3->Branch("pu_gen_wt_3",&pu_gen_wt_3,"pu_gen_wt_3/F");
  T3->Branch("gen_wt3",&gen_wt3,"gen_wt3/F");
  T3->Branch("tlumi_wt3",&tlumi_wt3,"tlumi_wt3/F");
  T3->Branch("tlumi_wt_hasGen3",&tlumi_wt_hasGen3,"tlumi_wt_hasGen3/F");
  T3->Branch("LSigma3",&LSigma3,"LSigma3/F");

  T->Branch("mu1pt",&mu1pt,"mu1pt/F");
  T->Branch("mu2pt",&mu2pt,"mu2pt/F");
  T->Branch("phopt",&phopt,"phopt/F");    
  T->Branch("muiso",&muiso,"muiso/F");
  T->Branch("drmu1pho",&drmu1pho,"drmu1pho/F");
  T->Branch("drmu2pho",&drmu2pho,"drmu2pho/F");
  T->Branch("drjpsipho",&drjpsipho,"drjpsipho/F");
  T->Branch("dphijpsipho",&dphijpsipho,"dphijpsipho/F");
  T->Branch("jpsim",&jpsim,"jpsim/F");
  T->Branch("Hm",&Hm,"Hm/F");
  T->Branch("finalm",&finalm,"finalm/F");
  T->Branch("ptratio",&ptratio,"ptratio/F");
  T->Branch("etratio",&etratio,"etratio/F");
  T->Branch("pu_wt",& pu_wt," pu_wt/F");
  T->Branch("pu_wt_",& pu_wt_," pu_wt_/F");
  T->Branch("pu_gen_wt",&pu_gen_wt,"pu_gen_wt/F");
  T->Branch("pu_gen_wt_",&pu_gen_wt_,"pu_gen_wt_/F");
  T->Branch("gen_wt",&gen_wt,"gen_wt/F");
  T->Branch("tlumi_wt",&tlumi_wt,"tlumi_wt/F");
  T->Branch("tlumi_wt_hasGen",&tlumi_wt_hasGen,"tlumi_wt_hasGen/F");
  T->Branch("LSigma",&LSigma,"LSigma/F");

  //xsec = xsec_;
  bookHistograms();
  if (createMVATree_)
    skimObj_ = std::make_unique<MVASkim>(mvaSkimFile_);
  else if (readMVA_){
    mvAna_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
  }
  TFile *mu_scale_factors;
  TFile *pho_scale_factors;
  TFile *trigg_scale_factors;
  string photon_ptCorr_file;
  string photon_vetoSF_file;
  string muon_ptCorr_file;
  string era;
  if (year() == 16) {
    mu_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/final_HZZ_SF_2016UL_mupogsysts_newLoose.root");
    if(isPre())  {
      pho_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/egammaEffi.txt_EGM2D_Pho_wp90_UL16.root");
      trigg_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/trigg_SF/Mu17Pho30SF_2016preVFP.root");
      muon_ptCorr_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/RoccoR2016aUL.txt";
      era = "2016preVFP";
      photon_ptCorr_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/EGM_ScaleUnc_2016Pre.json";
      photon_vetoSF_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/photon_2016Pre.json";
    }
    else  {
      pho_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/egammaEffi.txt_EGM2D_Pho_MVA90_UL16_postVFP.root");
      trigg_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/trigg_SF/Mu17Pho30SF_2016postVFP.root");
      muon_ptCorr_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/RoccoR2016bUL.txt";
      era = "2016postVFP";
      photon_ptCorr_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/EGM_ScaleUnc_2016Post.json";
      photon_vetoSF_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/photon_2016Post.json";
    }
  }
  
  if(year() == 17) {
    mu_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/final_HZZ_SF_2017UL_mupogsysts_newLoose.root");
    pho_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/egammaEffi.txt_EGM2D_PHO_MVA90_UL17.root");
    trigg_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/trigg_SF/Mu17Pho30SF_2017.root");
    muon_ptCorr_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/RoccoR2017UL.txt";
    era = "2017";
    photon_ptCorr_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/EGM_ScaleUnc_2017.json";
    photon_vetoSF_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/photon_2017.json";
  }

  if(year() == 18) {
    mu_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/final_HZZ_SF_2018UL_mupogsysts_newLoose.root");
    pho_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/egammaEffi.txt_EGM2D_Pho_wp90.root_UL18.root");
    trigg_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/trigg_SF/Mu17Pho30SF_2018.root");
    muon_ptCorr_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/RoccoR2018UL.txt";
    era = "2018";
    photon_ptCorr_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/EGM_ScaleUnc_2018.json";
    photon_vetoSF_file = "/afs/cern.ch/work/s/sbaradia/public/Analysis/HZtoJPsiGamma_final/object_SF/photon_2018.json";

  }

  RoccoR rc;
  rc.init(muon_ptCorr_file);
  
  TH2F *mu_scale;
  TH2F *pho_scale;
  TH2F *trigg_scale_EB;
  TH2F *trigg_scale_EE;

  if(isMC()) mu_scale = (TH2F*)mu_scale_factors->Get("FINAL");
  if(isMC()) pho_scale = (TH2F*)pho_scale_factors->Get("EGamma_SF2D");
  if(isMC()) trigg_scale_EB = (TH2F*)trigg_scale_factors->Get("Mu17Pho30_SFs_EB");
  if(isMC()) trigg_scale_EE = (TH2F*)trigg_scale_factors->Get("Mu17Pho30_SFs_EE");
  
  int num_pu_vtx = 0;
  
  Long64_t nbytes = 0, nb = 0;
  std::cout<<"No. of entries is"<<nentries<<std::endl;
  auto t_start = std::chrono::high_resolution_clock::now();
  auto t_end = std::chrono::high_resolution_clock::now();
  
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    //for (Long64_t jentry=0; jentry<100000; jentry++) {
    HLT_Mu17_Photon30_CaloIdL_L1ISO =0;
    HLT_Mu17_Photon30_IsoCaloId = 0;
    bool isGood = 1;
    int count_=0;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = chain->GetEntry(jentry);   nbytes += nb;
    if(jentry%500000 == 0 )
      {
	t_end = std::chrono::high_resolution_clock::now();
	std::cout<<"Processing Entry in event loop (Rate) : "<<jentry<<" / "<<nentries<<"  [ "<<100.0*jentry/nentries<<"  % ]  "
		 << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
		 <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries - jentry)/(1e-9 + jentry)* 0.001
		 <<std::endl;
      }
    if(!isMC()) {
      isGood = isGoodLumi(run,luminosityBlock);
    }
    if(!isGood) continue;
    //    std::cout<<isGoodLumi(run,luminosityBlock)<<std::endl;
    //if(!isMC() && !isGoodLumi(run,luminosityBlock)) continue;
    float weight_scalefactor_mu1 = 1.;
    float weight_scalefactor_mu2 = 1.;
    float weight_scalefactor_mu = 1.;
    float weight_scalefactor_pho = 1.;
    float weight_scalefactor_eleVeto = 1.;
    float weight_scalefactor_trigg =1.;

    float weight_scalefactor_SU_mu1 = 1.;
    float weight_scalefactor_SU_mu2 = 1.;
    //    float weight_scalefactor_SU_mu = 1.;
    float weight_scalefactor_SU_pho = 1.;
    float weight_scalefactor_SU_trigg = 1.;

    float weight_scalefactor_mu1P =0.;
    float weight_scalefactor_mu1N =0.;
    float weight_scalefactor_mu2P =0.;
    float weight_scalefactor_mu2N =0.;
    float weight_scalefactor_muP =0.;
    float weight_scalefactor_muN =0.;
    float weight_scalefactor_phoP =0.;
    float weight_scalefactor_phoN =0.;
    float weight_scalefactor_eleVetoP =0.;
    float weight_scalefactor_eleVetoN =0.;
    float weight_scalefactor_triggP =0.;
    float weight_scalefactor_triggN =0.;


    std::vector<float> vec_genmupt;
    std::vector<float> vec_genphopt;
    std::vector<float> vec_genmueta;
    std::vector<float> vec_genphoeta;
    TLorentzVector genmuP4,  genphoP4;
    std::vector<TLorentzVector> vec_genmuP4;
    std::vector<TLorentzVector> vec_genphoP4;
    //std::cout << "debug3 : " << std::endl;
    
    double weight=1.0;
    num_pu_vtx = Pileup_nPU;
    
    if (isMC()) {
      pu_wt = puWeight;
      gen_wt = Generator_weight;
      pu_gen_wt = puWeight*Generator_weight;

      tlumi_wt = (pu_wt * intLumi_ * xsec_)/nentries;
      tlumi_wt_hasGen = (pu_gen_wt * intLumi_ * xsec_)/geneventsumw;
      LSigma = intLumi_ * xsec_;
      
      if(find_SU) {
	pu_wtUp = puWeightUp;
	pu_wtDown = puWeightDown;
	
	pu_gen_wtUp = puWeightUp*Generator_weight;
	pu_gen_wtDown = puWeightDown*Generator_weight;
      }
    }
    else {
      pu_wt=1.0;
      gen_wt = 1.0;
      pu_gen_wt = 1.0;

      tlumi_wt = 1.0;
      tlumi_wt_hasGen = 1.0;
      if(find_SU) {
	pu_wtUp = 1.0;
	pu_wtDown = 1.0;
	pu_gen_wtUp = 1.0;
	pu_gen_wtDown = 1.0;
      }
    }
    
    double evtweight = -99.;
    
    if (isMC() && Generator_weight != 1.0) evtweight = Generator_weight;
    else evtweight = 1.0;
    evtWeightSum_ += evtweight;
    
    if(isMC()) {
      if (ispu_wt) weight = pu_wt;
      else if (isgen_wt) weight = gen_wt;
      else if (ispu_gen_wt) weight = pu_gen_wt;
       else if (istlumi_wt) weight = tlumi_wt;
       else if (istlumi_wt_hasGen) weight = tlumi_wt_hasGen;
    }
    
    outfile->cd();
    outfile->cd("HtoJPsiGamma");
     
     fillHist1D("generatorweight", Generator_weight);
     double scale =Generator_weight/geneventsumw;
     fillHist1D("generatorweight_S", scale);
     
     fillHist1D("puweight", weight);
     fillHist2D("pudist", num_pu_vtx, weight);
     
     fillHist1D("evtCutFlow", 0);
     fillHist1D("evtCutFlowWt", 0, pu_gen_wt);
     fillHist1D("nPV", PV_npvs, weight);
     fillHist1D("nPVgood", PV_npvsGood, weight);
     if (PV_npvsGood < 1) continue;
     fillHist1D("evtCutFlow", 1);
     fillHist1D("evtCutFlowWt", 1, pu_gen_wt);
     
     bool HLT;
     if (hlt() == 0) HLT = HLT_Mu17_Photon30_IsoCaloId;
     else if (hlt() == 1) HLT = HLT_Dimuon25_Jpsi;
     else if (hlt() == 2) HLT = HLT_DoubleMu20_7_Mass0to30_Photon23;
     else if (hlt() == 4) HLT = HLT_Mu17_Photon30_CaloIdL_L1ISO;
     else std::cout << "HLT not found" << std::endl;
     
     if (!HLT) continue;
     fillHist1D("evtCutFlow", 2);
     fillHist1D("evtCutFlowWt", 2, pu_gen_wt);

     TLorentzVector muLP4, muTP4, muTIsoP4, phoP4;
     std::vector<TLorentzVector> vec_muLP4, vec_muTP4, vec_muTIsoP4;

     std::vector<int> vec_muid;
     std::vector<TLorentzVector> vec_phoP4;
     std::map<TLorentzVector,int, ComparePt>map_muTP4;

     outfile->cd();
     outfile->cd("ObjectSelection");
     float mupt_max = -999.;
     int mu1id = -999;
     for (int i=0; i < nMuon; i++){
       fillHist1D("muCutFlow", 0);
       
       fillHist1D("muPt", Muon_pt[i], weight);
       fillHist1D("muEta", Muon_eta[i], weight);
       fillHist1D("muDxy", Muon_dxy[i], weight);
       fillHist1D("muDz", Muon_dz[i], weight);
       fillHist1D("muSIP3d", Muon_sip3d[i], weight);
       fillHist1D("muSIP3dincm", Muon_ip3d[i], weight);

       if (!Muon_looseId[i]) continue; // HZZ4L Loose ID
       fillHist1D("muCutFlow", 1);
       muLP4.SetPtEtaPhiM( Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
       vec_muLP4.push_back(muLP4);
       
       if (!Muon_tightId[i]) continue;
       fillHist1D("muCutFlow", 2);
       
       //Implementing Rochester Muon energy scale and reso correction                                                                            
       double muon_corr_SF =1;
       if(isMC()) {
	 
	 //	 if (Muon_genPartIdx[i]>= 0 && Muon_genPartIdx[i] < (sizeof(GenPart_pt)/sizeof(GenPart_pt[0]))) muon_corr_SF = rc.kSpreadMC(Muon_charge[i], Muon_pt[i] , Muon_eta[i], Muon_phi[i],GenPart_pt[Muon_genPartIdx[i]]);
	 if (static_cast<long unsigned int>(Muon_genPartIdx[i])>= 0 && static_cast<long unsigned int>(Muon_genPartIdx[i]) < (sizeof(GenPart_pt)/sizeof(GenPart_pt[0]))) muon_corr_SF = rc.kSpreadMC(Muon_charge[i], Muon_pt[i] , Muon_eta[i], Muon_phi[i],GenPart_pt[Muon_genPartIdx[i]]);
         else  muon_corr_SF = rc.kSmearMC(Muon_charge[i], Muon_pt[i] , Muon_eta[i], Muon_phi[i],Muon_nTrackerLayers[i],gRandom->Rndm());
       }
       else
         muon_corr_SF = rc.kScaleDT(Muon_charge[i], Muon_pt[i] , Muon_eta[i], Muon_phi[i]);
       
       Muon_pt[i]=Muon_pt[i]*muon_corr_SF;
       muTP4.SetPtEtaPhiM( Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
       vec_muTP4.push_back(muTP4);
       map_muTP4.insert(std::make_pair(muTP4, i));
       if (Muon_pt[i] > mupt_max){
	 mupt_max = Muon_pt[i];
	 mu1id = i;
       }
       
     }
     
     /*     int mu2id = -999.;
     int mu2ndmaxpt = -999.;
     if (nMuon > 0) {
       for (int j=0; j < nMuon; j++){
	 if (j==mu1id) continue;
	 if (Muon_pt[j] > mu2ndmaxpt){
	   mu2ndmaxpt = Muon_pt[j];
	   mu2id = j;
	 }
       }
       }*/

     outfile->cd();
     outfile->cd("HtoJPsiGamma");
     
     fillHist1D("nLoosemuons", vec_muLP4.size(), weight);
     fillHist1D("nTightmuons", vec_muTP4.size(), weight);
     if ( vec_muTP4.size() < 2 ) continue;
     fillHist1D("evtCutFlow", 3);
     fillHist1D("evtCutFlowWt", 3, pu_gen_wt);
     std::sort(vec_muTP4.begin(), vec_muTP4.end(), PtComparatorTL<TLorentzVector>());
     
     auto it1 = map_muTP4.begin();
     auto it2 = std::next(it1, 1); // Move the iterator to the second element
     //     int mu1id_ = it1->second;
     int mu2id = it2->second;
     //     std::cout<<"MU1ID:"<<mu1id<<", MU1ID_:"<<mu1id_<<", MU2ID:"<<mu2id<<std::endl;

     fillHist1D("muPfRelIso03all", Muon_pfRelIso03_all[mu1id], weight);
     fillHist1D("muPfRelIso03chg", Muon_pfRelIso03_chg[mu1id], weight);
     fillHist1D("muPfRelIso04all", Muon_pfRelIso04_all[mu1id], weight);
     //     fillHist1D("muPfRelIso03all", Muon_pfRelIso03_all[mu1id], wt);
     fillHist1D("muPfRelIso03allByMuSubleadpt", Muon_pfRelIso03_all[mu1id] / vec_muTP4[1].Pt(), weight);
     if (Muon_pfRelIso03_all[mu1id] >= 0.35) continue;
     muiso = Muon_pfRelIso03_all[mu1id];
     mu1_sip3d = Muon_sip3d[mu1id]; //VAR
     mu1iso = Muon_pfRelIso03_all[mu1id];  //VAR
     fillHist1D("evtCutFlow", 4);
     fillHist1D("evtCutFlowWt", 4, pu_gen_wt);
     
     fillHist1D("muLeadpt", vec_muTP4[0].Pt(), weight);
     if (vec_muTP4[0].Pt() <= 20) continue;
     mu1pt = vec_muTP4[0].Pt();                   //VAR
     mu1eta = vec_muTP4[0].Eta();  //VAR
     mu1phi = vec_muTP4[0].Phi(); //VAR
     fillHist1D("evtCutFlow", 5);
     fillHist1D("evtCutFlowWt", 5, pu_gen_wt);

     fillHist1D("muSubleadpt", vec_muTP4[1].Pt(), weight);     
     if (vec_muTP4[1].Pt() <= 4) continue;
     mu2pt = vec_muTP4[1].Pt(); //VAR
     mu2eta = vec_muTP4[1].Eta();  //VAR                                              
     mu2phi = vec_muTP4[1].Phi();  //VAR
     fillHist1D("evtCutFlow", 6);
     fillHist1D("evtCutFlowWt", 6, pu_gen_wt);

     outfile->cd();
     outfile->cd("ObjectSelection");
     
     float phopt_max = -999.;
     int pho1id = -999;
     
     for (int i = 0 ; i < nPhoton ; i++){
       fillHist1D("photonCutFlow", 0);
       fillHist1D("photonCutFlowWt", 0, pu_gen_wt);

       if (Photon_pt[i] <= 33) continue;
       fillHist1D("photonCutFlow", 1);
       fillHist1D("photonCutFlowWt", 1, pu_gen_wt);
       
       if (Photon_isScEtaEB[i]) fillHist1D("nPhoScEta", 0, weight); // 0 = ScEtaEB, 1 = ScEtaEE                                                              
       else if (Photon_isScEtaEE[i]) fillHist1D("nPhoScEta", 1, weight);
       if (! (Photon_isScEtaEB[i] || Photon_isScEtaEE[i])) continue;
       fillHist1D("photonCutFlow", 2);
       fillHist1D("photonCutFlowWt",2, pu_gen_wt);       

       //       if(!Photon_electronVeto[i]) continue;
       if(Photon_pixelSeed[i]) continue;
       fillHist1D("photonCutFlow", 3);
       fillHist1D("photonCutFlowWt", 3, pu_gen_wt);

       if (Photon_isScEtaEB[i]) { 
	 fillHist1D("nPhoScEta_", 0); // 0 = ScEtaEB, 1 = ScEtaEE                                                            
	 fillHist1D("nPhoScEta_W", 0, pu_gen_wt);
       }
       else if (Photon_isScEtaEE[i]) {
	 fillHist1D("nPhoScEta_", 1);
	 fillHist1D("nPhoScEta_W", 1, pu_gen_wt);
       }
       
       if (Photon_mvaID_WP90[i] == false) continue;
       count_++;
       fillHist1D("photonCutFlow", 4);
       fillHist1D("photonCutFlowWt", 4, pu_gen_wt);

       if(i==0)
         fillHist1D("phoPt", Photon_pt[i], weight);
       if(count_ ==1)
         fillHist1D("phoPt_MVA", Photon_pt[i], weight);
       

       fillHist1D("phoR9", Photon_r9[i], weight);
       if (Photon_isScEtaEB[i]) fillHist1D("phoR9_ScEB", Photon_r9[i], weight);
       if (Photon_isScEtaEE[i]) fillHist1D("phoR9_ScEE", Photon_r9[i], weight);
       
       phoP4.SetPtEtaPhiM(Photon_pt[i], Photon_eta[i], Photon_phi[i], Photon_mass[i]);
       vec_phoP4.push_back(phoP4);
       
       if (Photon_pt[i] > phopt_max){
         phopt_max = Photon_pt[i];
         pho1id = i;
       }
     }
     
     phopt = Photon_pt[pho1id];   //VAR
     phoeta = Photon_eta[pho1id];  //VAR
     phophi = Photon_phi[pho1id];  //VAR
     photon_r9=Photon_r9[pho1id]; //VAR
     
     outfile->cd();
     outfile->cd("HtoJPsiGamma");
     fillHist1D("nPhotons", vec_phoP4.size(), weight);
     if (vec_phoP4.size() < 1) continue;
     fillHist1D("evtCutFlow", 7);
     fillHist1D("evtCutFlowWt", 7, pu_gen_wt);
     std::sort(vec_phoP4.begin(), vec_phoP4.end(), PtComparatorTL<TLorentzVector>());
     
     drmu1pho = vec_muTP4[0].DeltaR(vec_phoP4[0]);  //VAR     
     fillHist1D("dRmu1pho", drmu1pho, weight);
     if (drmu1pho <= 1 ) continue;
     fillHist1D("evtCutFlow", 8);
     fillHist1D("evtCutFlowWt", 8, pu_gen_wt);
     
     
     drmu2pho = vec_muTP4[1].DeltaR(vec_phoP4[0]); //VAR
     fillHist1D("dRmu2pho",drmu2pho, weight);
     if (drmu2pho <= 1) continue;
     fillHist1D("evtCutFlow", 9);     
     fillHist1D("evtCutFlowWt", 9, pu_gen_wt);
     
     
     drjpsipho = (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]); //VAR
     fillHist1D("dRdimupho", drjpsipho, weight);
     if (drjpsipho <= 2 ) continue;
     fillHist1D("evtCutFlow", 10);  
     fillHist1D("evtCutFlowWt", 10, pu_gen_wt);
     
     

     drmu1mu2 = vec_muTP4[0].DeltaR(vec_muTP4[1]); //VAR
     dphimu1mu2 = vec_muTP4[0].DeltaPhi(vec_muTP4[1]); //VAR
     mod_dphimu1mu2 = std::abs(dphimu1mu2); //VAR
     dphijpsipho = (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]); //VAR
     mod_dphijpsipho = std::abs(dphijpsipho); //VAR     

     dphimu1pho = (vec_muTP4[0]).DeltaPhi(vec_phoP4[0]); //VAR                                                                                                                                         
     mod_dphimu1pho = std::abs(dphimu1pho); //VAR                                                                                                                                               
     dphimu2pho = (vec_muTP4[1]).DeltaPhi(vec_phoP4[0]); //VAR                                                                                                                                    
     mod_dphimu2pho = std::abs(dphimu2pho); //VAR   


     fillHist1D("dPhidimupho",dphijpsipho, weight);
     if (mod_dphijpsipho <= 1.5 ) continue;
     fillHist1D("evtCutFlow", 11);
     fillHist1D("evtCutFlowWt", 11, pu_gen_wt);





     jpsim = (vec_muTP4[0] + vec_muTP4[1]).M(); //VAR                                                                                                                                                    
     jpsipt = (vec_muTP4[0] + vec_muTP4[1]).Pt(); //VAR   
     fillHist1D("jpsimass_bmasscut",jpsim, weight);
     if (jpsim < 3.0 || jpsim > 3.2) continue;  //For Signal MVA        
     fillHist1D("evtCutFlow", 12);
     fillHist1D("evtCutFlowWt", 12, pu_gen_wt);
     fillHist1D("jpsimass_amasscut",jpsim, weight);


     Hm = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(); //VAR                                                                
     Hmpt = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).Pt();
     ptratio = (vec_muTP4[0] + vec_muTP4[1]).Pt()/ (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(); //VAR
     etratio = vec_phoP4[0].Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();  //VAR
     jpsiSSpt = mu1pt + mu2pt;
     HmSSpt = mu1pt+mu2pt+phopt;

     /*
     float mvaOutput = -999.9;

     // Fill the flat tree here                                                                                                                       
     if (createMVATree_ && skimObj_ != nullptr) {
       //       std::cout<<"I am filling"<<std::endl;
       TreeVariables varList;
       varList.phopt = phopt;
       varList.phoeta = phoeta;
       varList.phophi = phophi; //NEW
       varList.photon_r9 = photon_r9;
       varList.mu1pt = mu1pt;
       varList.mu1eta = mu1eta;
       varList.mu1phi = mu1phi;
       varList.mu2pt = mu2pt;
       varList.mu2eta = mu2eta;
       varList.mu2phi = mu2phi;
       varList.mu1iso = mu1iso;
       varList.mu2iso = mu2iso;
       varList.mu1_sip3d = mu1_sip3d;
       varList.mu2_sip3d = mu2_sip3d;
       varList.dphijpsipho = dphijpsipho;
       varList.mod_dphijpsipho = mod_dphijpsipho;
       varList.dphimu1pho = dphimu1pho;
       varList.mod_dphimu1pho = mod_dphimu1pho;
       varList.dphimu2pho = dphimu2pho;
       varList.mod_dphimu2pho = mod_dphimu2pho;
       varList.dphimu1mu2 = dphimu1mu2; //NEW                                                                                                                                          
       varList.mod_dphimu1mu2 = mod_dphimu1mu2; //NEW   
       varList.drmu1pho = drmu1pho;
       varList.drmu2pho = drmu2pho;
       varList.drjpsipho = drjpsipho;
       varList.drmu1mu2 = drmu1mu2; //NEW
       varList.jpsim = jpsim;
       varList.jpsipt = jpsipt;
       varList.Hm = Hm;
       varList.Hmpt = Hmpt;
       varList.ptratio = ptratio;
       varList.etratio = etratio;
       varList.jpsiSSpt = jpsiSSpt; //NEW
       varList.HmSSpt = HmSSpt; //NEW
       varList.tlumi_wt = tlumi_wt;
       varList.tlumi_wt_hasGen = tlumi_wt_hasGen;
       skimObj_->fill(varList);
     }
     else if (readMVA_ && (mvAna_ != nullptr)) {
       InputVariables varList;
       varList.phopt = phopt;
       //       varList.phoeta = phoeta;
       //       varList.phophi = phophi;
       //varList.photon_r9 = photon_r9;
       varList.mu1pt = mu1pt;
       //varList.mu1eta = mu1eta;
       //       varList.mu1phi = mu1phi;
       varList.mu2pt = mu2pt;
       //varList.mu2eta = mu2eta;
       //       varList.mu2phi = mu2phi;
       varList.mu1iso = mu1iso;
       //       varList.mu2iso = mu2iso;
       //       varList.mu1_sip3d = mu1_sip3d;
       //       varList.mu2_sip3d = mu2_sip3d;
       //       varList.dphijpsipho = dphijpsipho;
       varList.mod_dphijpsipho = mod_dphijpsipho;
       //       varList.dphimu1pho = dphimu1pho;
       // varList.mod_dphimu1pho = mod_dphimu1pho;
       //       varList.dphimu2pho = dphimu2pho;
       // varList.mod_dphimu2pho = mod_dphimu2pho;
       //       varList.dphimu1mu2 = dphimu1mu2;                                                                                                                                             
       //varList.mod_dphimu1mu2 = mod_dphimu1mu2; 
       varList.drmu1pho = drmu1pho;
       varList.drmu2pho = drmu2pho;
       varList.drjpsipho = drjpsipho;
       //       varList.drmu1mu2 = drmu1mu2;
       //       varList.jpsim = jpsim;
       //       varList.jpsipt = jpsipt;
       //       varList.Hm = Hm;
       //       varList.Hmpt = Hmpt;
       // varList.ptratio = ptratio;
       //varList.etratio = etratio;
       //       varList.jpsiSSpt = mu1pt+mu2pt;
       //varList.HmSSpt = mu1pt+mu2pt+phopt;
       mvaOutput = mvAna_->evaluate(mvAlgo_, varList);
       //       varList.tlumi_wt = tlumi_wt;
       //       varList.tlumi_wt_hasGen = tlumi_wt_hasGen;
     }
     outfile->cd("HtoJPsiGamma");
     //     std::cout<<"the value of weight is "<<weight<<std::endl;
     if(readMVA_) {
       fillHist1D("mvaoutput_bdtg_tlumi_wt_hasGen", mvaOutput, tlumi_wt_hasGen);
       fillHist1D("mvaoutput_bdtg_pu_gen_wt", mvaOutput, pu_gen_wt);
       fillHist1D("mvaoutput_bdtg", mvaOutput);
       fillHist2D("scoreVsInv",(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(),mvaOutput);
       fillHist1D("zmass_bscorecut", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), weight);
       if(mvaOutput < -0.30) continue; 
       }*/
     //  fillHist1D("evtCutFlow", 14);
     //fillHist1D("evtCutFlowWt", 14, pu_gen_wt);
     
     outfile->cd("HtoJPsiGamma");
     fillHist1D("zmass_bmasscut",Hm, weight);
     if (Hm <= 100.0 || Hm >= 150.0) continue;
     fillHist1D("evtCutFlow", 13);
     fillHist1D("evtCutFlowWt", 13, pu_gen_wt);
     
     fillHist1D("zmass_amasscut",Hm, weight);
     fillHist1D("zmass_bratiocut1",Hm, weight);
     fillHist1D("ratiocutdimu",ptratio, weight);
     if (ptratio <= 0.28) continue;
     fillHist1D("evtCutFlow", 14);
     fillHist1D("evtCutFlowWt", 14, pu_gen_wt);

     fillHist1D("zmass_bratiocut2",Hm, weight);
     fillHist1D("ratiocutpho",etratio, weight);
     if (etratio<= 0.28) continue;
     fillHist1D("evtCutFlow", 15);
     fillHist1D("evtCutFlowWt", 15, pu_gen_wt);

     if(isMC()) {
       string tag_eleVeto_SF;
       auto cset1 = correction::CorrectionSet::from_file(photon_vetoSF_file);
       if(Photon_isScEtaEB[pho1id] == 1 ) tag_eleVeto_SF = "EBInc";
       else tag_eleVeto_SF = "EEInc";
       
       auto Map1 = cset1->at("UL-Photon-ID-SF");     //Using from json
       //       auto Map2 = cset1->at("UL-Photon-CSEV-SF");
       auto Map2 = cset1->at("UL-Photon-PixVeto-SF");
       
       weight_scalefactor_pho *= Map1->evaluate({era,"sf","wp90",vec_phoP4[0].Eta(),vec_phoP4[0].Pt()});
       weight_scalefactor_eleVeto *= Map2->evaluate({era,"sf","MVA",tag_eleVeto_SF});
       
       weight_scalefactor_mu1 *= HtoJPsiGamma::SFCalc(mu_scale, vec_muTP4[0].Pt(), vec_muTP4[0].Eta(), 1);
       weight_scalefactor_mu2 *= HtoJPsiGamma::SFCalc(mu_scale, vec_muTP4[1].Pt(), vec_muTP4[1].Eta(), 1);
       weight_scalefactor_mu *= weight_scalefactor_mu1 * weight_scalefactor_mu2;
       
       if(find_SU) {
         weight_scalefactor_phoP = Map1->evaluate({era, "sfup","wp90",vec_phoP4[0].Eta(),vec_phoP4[0].Pt()});
	 weight_scalefactor_phoN = Map1->evaluate({era, "sfdown","wp90",vec_phoP4[0].Eta(),vec_phoP4[0].Pt()});

	 weight_scalefactor_eleVetoP = Map2->evaluate({era,"sfup","MVA",tag_eleVeto_SF});
         weight_scalefactor_eleVetoN = Map2->evaluate({era,"sfdown","MVA",tag_eleVeto_SF});
	 

	 weight_scalefactor_SU_mu1 *= HtoJPsiGamma::SFCalcSU(mu_scale, vec_muTP4[0].Pt(), vec_muTP4[0].Eta(), 1);
         weight_scalefactor_SU_mu2 *= HtoJPsiGamma::SFCalcSU(mu_scale, vec_muTP4[1].Pt(), vec_muTP4[1].Eta(), 1);
         weight_scalefactor_mu1P = weight_scalefactor_mu1 + weight_scalefactor_SU_mu1;
         weight_scalefactor_mu1N = weight_scalefactor_mu1 - weight_scalefactor_SU_mu1;
         weight_scalefactor_mu2P = weight_scalefactor_mu2 + weight_scalefactor_SU_mu2;
         weight_scalefactor_mu2N = weight_scalefactor_mu2 - weight_scalefactor_SU_mu2;

         weight_scalefactor_muP = weight_scalefactor_mu1P * weight_scalefactor_mu2P;
         weight_scalefactor_muN= weight_scalefactor_mu1N * weight_scalefactor_mu2N;
       }
       if(Photon_isScEtaEB[pho1id]) {
         weight_scalefactor_trigg *= HtoJPsiGamma::triggSFCalc(trigg_scale_EB,vec_muTP4[0].Pt(),vec_phoP4[0].Pt());
         if(find_SU) {
           weight_scalefactor_SU_trigg *= HtoJPsiGamma::triggSFCalcSU(trigg_scale_EB,vec_muTP4[0].Pt(),vec_phoP4[0].Pt());
           weight_scalefactor_triggP =weight_scalefactor_trigg +weight_scalefactor_SU_trigg;
           weight_scalefactor_triggN =weight_scalefactor_trigg -weight_scalefactor_SU_trigg;
         }
       }
       else if(Photon_isScEtaEE[pho1id]) {
         weight_scalefactor_trigg *= HtoJPsiGamma::triggSFCalc(trigg_scale_EE,vec_muTP4[0].Pt(),vec_phoP4[0].Pt());
         if(find_SU) {
           weight_scalefactor_SU_trigg *= HtoJPsiGamma::triggSFCalcSU(trigg_scale_EE,vec_muTP4[0].Pt(),vec_phoP4[0].Pt());
           weight_scalefactor_triggP =weight_scalefactor_trigg +weight_scalefactor_SU_trigg;
           weight_scalefactor_triggN =weight_scalefactor_trigg -weight_scalefactor_SU_trigg;
         }
       }
     
       if (weight_scalefactor_mu  > 0. && weight_scalefactor_pho > 0. && weight_scalefactor_trigg > 0. && weight_scalefactor_triggP > 0. && weight_scalefactor_triggN > 0. && weight_scalefactor_muP > 0. && weight_scalefactor_muN > 0. && weight_scalefactor_phoP > 0. && weight_scalefactor_phoN > 0. && weight_scalefactor_eleVeto > 0.) {

         pu_gen_wt_triggP= pu_gen_wt * weight_scalefactor_mu *  weight_scalefactor_pho * weight_scalefactor_eleVeto * weight_scalefactor_triggP;
         pu_gen_wt_triggN= pu_gen_wt * weight_scalefactor_mu *  weight_scalefactor_pho * weight_scalefactor_eleVeto * weight_scalefactor_triggN;

         pu_gen_wt_muP= pu_gen_wt * weight_scalefactor_muP * weight_scalefactor_pho * weight_scalefactor_eleVeto * weight_scalefactor_trigg;
         pu_gen_wt_muN= pu_gen_wt * weight_scalefactor_muN * weight_scalefactor_pho * weight_scalefactor_eleVeto * weight_scalefactor_trigg;

         pu_gen_wt_phoP= pu_gen_wt * weight_scalefactor_mu * weight_scalefactor_phoP * weight_scalefactor_eleVeto * weight_scalefactor_trigg;
         pu_gen_wt_phoN= pu_gen_wt * weight_scalefactor_mu * weight_scalefactor_phoN * weight_scalefactor_eleVeto * weight_scalefactor_trigg;

	 pu_gen_wt_eleVetoP= pu_gen_wt * weight_scalefactor_mu * weight_scalefactor_pho * weight_scalefactor_eleVetoP * weight_scalefactor_trigg;
	 pu_gen_wt_eleVetoN= pu_gen_wt * weight_scalefactor_mu * weight_scalefactor_pho * weight_scalefactor_eleVetoN * weight_scalefactor_trigg;


         gen_wt*= weight_scalefactor_mu* weight_scalefactor_pho * weight_scalefactor_eleVeto *  weight_scalefactor_trigg;
         pu_wt*= weight_scalefactor_mu* weight_scalefactor_pho *  weight_scalefactor_eleVeto * weight_scalefactor_trigg;
         pu_gen_wt*= weight_scalefactor_mu* weight_scalefactor_pho *  weight_scalefactor_eleVeto * weight_scalefactor_trigg;

         tlumi_wt*= weight_scalefactor_mu* weight_scalefactor_pho * weight_scalefactor_eleVeto *  weight_scalefactor_trigg;
         tlumi_wt_hasGen*= weight_scalefactor_mu * weight_scalefactor_pho * weight_scalefactor_eleVeto *  weight_scalefactor_trigg;

         pu_wtUp*=weight_scalefactor_mu * weight_scalefactor_pho * weight_scalefactor_eleVeto *  weight_scalefactor_trigg;
         pu_wtDown*=weight_scalefactor_mu * weight_scalefactor_pho * weight_scalefactor_eleVeto *  weight_scalefactor_trigg;

         pu_gen_wtUp*=weight_scalefactor_mu* weight_scalefactor_pho * weight_scalefactor_eleVeto *  weight_scalefactor_trigg;
         pu_gen_wtDown*=weight_scalefactor_mu* weight_scalefactor_pho * weight_scalefactor_eleVeto *  weight_scalefactor_trigg;
       }
     }
     //     std::cout<<"weight_scalefactor_eleVeto" << weight_scalefactor_eleVeto<<" , "<<"weight_scalefactor_pho"<<weight_scalefactor_pho<<" , "<<"weight_scalefactor_eleVetoP"<<weight_scalefactor_eleVetoP<<" , "<<"weight_scalefactor_eleVetoN"<<weight_scalefactor_eleVetoN<<" , "<<"weight_scalefactor_phoP"<<weight_scalefactor_phoP<<std::endl;


     fillHist1D("evtCutFlow", 16);
     fillHist1D("evtCutFlowWt", 16, pu_gen_wt);
     
     fillHist1D("nPhoCat", 0); 
     fillHist1D("nPhoCatWt", 0, pu_gen_wt);                                        

     if(Muon_jetIdx[mu1id] < 0) {
       fillHist1D("Mu1_jetIDX", 0);
       fillHist1D("Mu1_jetIDX_Wt", 0,pu_gen_wt);
     }
     else {
       fillHist1D("Mu1_jetIDX", 1);
       fillHist1D("Mu1_jetIDX_Wt",1, pu_gen_wt);
     }

     if(Muon_jetIdx[mu2id] < 0) {
       fillHist1D("Mu2_jetIDX", 0);
       fillHist1D("Mu2_jetIDX_Wt", 0,pu_gen_wt);
     }
     else {
       fillHist1D("Mu2_jetIDX", 1);
       fillHist1D("Mu2_jetIDX_Wt",1, pu_gen_wt);
     }

     if(Photon_jetIdx[pho1id] < 0) {
       fillHist1D("Pho_jetIDX", 0);
       fillHist1D("Pho_jetIDX_Wt", 0,pu_gen_wt);
     }
     else {
       fillHist1D("Pho_jetIDX", 1);
       fillHist1D("Pho_jetIDX_Wt",1, pu_gen_wt);
     }

     if(Photon_electronIdx[pho1id] < 0) {
       fillHist1D("Pho_eleIDX", 0);
       fillHist1D("Pho_eleIDX_Wt",0, pu_gen_wt);
     }
     else {
       fillHist1D("Pho_eleIDX", 1);
       fillHist1D("Pho_eleIDX_Wt",1, pu_gen_wt);
     }



     fillHist1D("nPhoCatWt_triggP", 0, pu_gen_wt_triggP);
     fillHist1D("nPhoCatWt_triggN", 0, pu_gen_wt_triggN);
     fillHist1D("nPhoCatWt_muP", 0, pu_gen_wt_muP);
     fillHist1D("nPhoCatWt_muN", 0, pu_gen_wt_muN);
     fillHist1D("nPhoCatWt_phoP", 0, pu_gen_wt_phoP);
     fillHist1D("nPhoCatWt_phoN", 0, pu_gen_wt_phoN);
     fillHist1D("nPhoCatWt_eleVetoP", 0, pu_gen_wt_eleVetoP);
     fillHist1D("nPhoCatWt_eleVetoN", 0, pu_gen_wt_eleVetoN);
     fillHist1D("nPhoCatWt_PuUp", 0, pu_gen_wtUp);
     fillHist1D("nPhoCatWt_PuDown", 0, pu_gen_wtDown);
     //     photon_r9=Photon_r9[pho1id]; //VAR
     if (Photon_isScEtaEB[pho1id] && Photon_r9[pho1id] >= 0.96 )   { 
       fillHist1D("nPhoCat", 1); 
       fillHist1D("nPhoCatWt", 1, pu_gen_wt);
       fillHist1D("nPhoCatWt_triggP", 1, pu_gen_wt_triggP);
       fillHist1D("nPhoCatWt_triggN", 1, pu_gen_wt_triggN);
       fillHist1D("nPhoCatWt_muP", 1, pu_gen_wt_muP);
       fillHist1D("nPhoCatWt_muN", 1, pu_gen_wt_muN);
       fillHist1D("nPhoCatWt_phoP", 1, pu_gen_wt_phoP);
       fillHist1D("nPhoCatWt_phoN", 1, pu_gen_wt_phoN);
       fillHist1D("nPhoCatWt_eleVetoP", 1, pu_gen_wt_eleVetoP);
       fillHist1D("nPhoCatWt_eleVetoN", 1, pu_gen_wt_eleVetoN);
       fillHist1D("nPhoCatWt_PuUp", 1, pu_gen_wtUp);
       fillHist1D("nPhoCatWt_PuDown", 1, pu_gen_wtDown);
       fillHist1D("muLeadpt_cat1", vec_muTP4[0].Pt(), weight);
       fillHist1D("muLeadeta_cat1", vec_muTP4[0].Eta(), weight);
       fillHist1D("muSubleadpt_cat1", vec_muTP4[1].Pt(), weight);
       fillHist1D("muSubleadeta_cat1", vec_muTP4[1].Eta(), weight);
       fillHist1D("phoPt_cat1", vec_phoP4[0].Pt(), weight);
       fillHist1D("phoEta_cat1", vec_phoP4[0].Eta(), weight);
       fillHist1D("dRmu1pho_cat1", vec_muTP4[0].DeltaR(vec_phoP4[0]), weight);
       fillHist1D("dRmu2pho_cat1", vec_muTP4[1].DeltaR(vec_phoP4[0]), weight);
       fillHist1D("dRdimupho_cat1", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), weight);
       fillHist1D("dPhidimupho_cat1", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), weight);
       fillHist1D("jpsimass_cat1", (vec_muTP4[0] + vec_muTP4[1]).M(), weight);
       fillHist1D("dRmu1mu2_cat1", vec_muTP4[0].DeltaR(vec_muTP4[1]), weight); 
       fillHist1D("diMuPt_cat1", (vec_muTP4[0] + vec_muTP4[1]).Pt(), weight);
       fillHist1D("finalZPt_cat1", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).Pt(), weight);
       fillHist1D("zmass_cat1", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), weight);
       fillHist1D("vtxdist_x" , PV_x - (vec_muTP4[0] + vec_muTP4[1]).X(), weight);
       fillHist1D("vtxdist_y" , PV_y - (vec_muTP4[0] + vec_muTP4[1]).Y(), weight);
       fillHist1D("vtxdist_z" , PV_z - (vec_muTP4[0] + vec_muTP4[1]).Z(), weight);
       zm1 = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
       pu_wt1 = pu_wt;
       gen_wt1 = gen_wt;
       pu_gen_wt1 = pu_gen_wt;
       tlumi_wt1 = tlumi_wt;
       tlumi_wt_hasGen1 = tlumi_wt_hasGen;
       LSigma1 = LSigma;
       T1->Fill();
     }
     else if (Photon_isScEtaEB[pho1id] && Photon_r9[pho1id] < 0.96 ){
       fillHist1D("nPhoCat", 2); 
       fillHist1D("nPhoCatWt", 2, pu_gen_wt); 
       fillHist1D("nPhoCatWt_triggP", 2, pu_gen_wt_triggP);
       fillHist1D("nPhoCatWt_triggN", 2, pu_gen_wt_triggN);
       fillHist1D("nPhoCatWt_muP", 2, pu_gen_wt_muP);
       fillHist1D("nPhoCatWt_muN", 2, pu_gen_wt_muN);
       fillHist1D("nPhoCatWt_phoP", 2, pu_gen_wt_phoP);
       fillHist1D("nPhoCatWt_phoN", 2, pu_gen_wt_phoN);
       fillHist1D("nPhoCatWt_eleVetoP", 2, pu_gen_wt_eleVetoP);
       fillHist1D("nPhoCatWt_eleVetoN", 2, pu_gen_wt_eleVetoN);
       fillHist1D("nPhoCatWt_PuUp", 2, pu_gen_wtUp);
       fillHist1D("nPhoCatWt_PuDown", 2, pu_gen_wtDown);
       fillHist1D("muLeadpt_cat2", vec_muTP4[0].Pt(), weight);
       fillHist1D("muLeadeta_cat2", vec_muTP4[0].Eta(), weight);
       fillHist1D("muSubleadpt_cat2", vec_muTP4[1].Pt(), weight);
       fillHist1D("muSubleadeta_cat2", vec_muTP4[1].Eta(), weight);
       fillHist1D("phoPt_cat2", vec_phoP4[0].Pt(), weight);
       fillHist1D("phoEta_cat2", vec_phoP4[0].Eta(), weight);
       fillHist1D("dRmu1pho_cat2", vec_muTP4[0].DeltaR(vec_phoP4[0]), weight);
       fillHist1D("dRmu2pho_cat2", vec_muTP4[1].DeltaR(vec_phoP4[0]), weight);
       fillHist1D("dRdimupho_cat2", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), weight);
       fillHist1D("dPhidimupho_cat2", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), weight);
       fillHist1D("jpsimass_cat2", (vec_muTP4[0] + vec_muTP4[1]).M(), weight);
       fillHist1D("dRmu1mu2_cat2", vec_muTP4[0].DeltaR(vec_muTP4[1]), weight);
       fillHist1D("diMuPt_cat2", (vec_muTP4[0] + vec_muTP4[1]).Pt(), weight);
       fillHist1D("finalZPt_cat2", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).Pt(), weight);
       fillHist1D("zmass_cat2", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), weight);
       zm2 = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
       pu_wt2= pu_wt;
       gen_wt2= gen_wt;
       pu_gen_wt2 = pu_gen_wt;
       tlumi_wt2 = tlumi_wt;
       tlumi_wt_hasGen2 = tlumi_wt_hasGen;
       LSigma2 = LSigma;
       T2->Fill();
     } 
     else if (Photon_isScEtaEE[pho1id] ) {
       fillHist1D("nPhoCat", 3); 
       fillHist1D("nPhoCatWt", 3, pu_gen_wt);
       fillHist1D("nPhoCatWt_triggP", 3, pu_gen_wt_triggP);
       fillHist1D("nPhoCatWt_triggN", 3, pu_gen_wt_triggN);
       fillHist1D("nPhoCatWt_muP", 3, pu_gen_wt_muP);
       fillHist1D("nPhoCatWt_muN", 3, pu_gen_wt_muN);
       fillHist1D("nPhoCatWt_phoP", 3, pu_gen_wt_phoP);
       fillHist1D("nPhoCatWt_phoN", 3, pu_gen_wt_phoN);
       fillHist1D("nPhoCatWt_eleVetoP", 3, pu_gen_wt_eleVetoP);
       fillHist1D("nPhoCatWt_eleVetoN", 3, pu_gen_wt_eleVetoN);
       fillHist1D("nPhoCatWt_PuUp", 3, pu_gen_wtUp);
       fillHist1D("nPhoCatWt_PuDown", 3, pu_gen_wtDown);
       fillHist1D("muLeadpt_cat3", vec_muTP4[0].Pt(), weight);
       fillHist1D("muLeadeta_cat3", vec_muTP4[0].Eta(), weight);
       fillHist1D("muSubleadpt_cat3", vec_muTP4[1].Pt(), weight);
       fillHist1D("muSubleadeta_cat3", vec_muTP4[1].Eta(), weight);
       fillHist1D("phoPt_cat3", vec_phoP4[0].Pt(), weight);
       fillHist1D("phoEta_cat3", vec_phoP4[0].Eta(), weight);
       fillHist1D("dRmu1pho_cat3", vec_muTP4[0].DeltaR(vec_phoP4[0]), weight);
       fillHist1D("dRmu2pho_cat3", vec_muTP4[1].DeltaR(vec_phoP4[0]), weight);
       fillHist1D("dRdimupho_cat3", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), weight);
       fillHist1D("dPhidimupho_cat3", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), weight);
       fillHist1D("jpsimass_cat3", (vec_muTP4[0] + vec_muTP4[1]).M(), weight);
       fillHist1D("dRmu1mu2_cat3", vec_muTP4[0].DeltaR(vec_muTP4[1]), weight);
       fillHist1D("diMuPt_cat3", (vec_muTP4[0] + vec_muTP4[1]).Pt(), weight);
       fillHist1D("finalZPt_cat3", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).Pt(), weight);
       fillHist1D("zmass_cat3", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), weight);
       zm3 = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
       pu_wt3= pu_wt;
       gen_wt3= gen_wt;
       pu_gen_wt3 = pu_gen_wt;
       tlumi_wt3 = tlumi_wt;
       tlumi_wt_hasGen3 = tlumi_wt_hasGen;
       LSigma3 = LSigma;
       T3->Fill();
     }

     T->Fill();
     
     fillHist1D("jpsimass_afinalcut", (vec_muTP4[0] + vec_muTP4[1]).M(), weight);

     fillHist1D("zmass_afinalcut", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), weight);
   }
   saveHistograms();
   outfile->cd();
   outfile->cd("ObjectSelection");
   objectEfficiency();
   showCategoryYield();
   outfile->cd();
   outfile->cd("HtoJPsiGamma");
   T->Write();
   T1->Write();
   T2->Write();
   T3->Write();
   eventEfficiency();
   outfile->Close();
   cout << "Now the ouput file is closed" << endl;
   cout<<"total geneventsumwt from Runs:"<<geneventsumw<<endl;
   cout<<"total geneventsumw from sum mc genWeight:"<<evtWeightSum_<< endl;
}



double HtoJPsiGamma::SFCalc(TH2F* h, double pt, double eta, int id){
  double tmp_pt = pt;
  if(id == 1) {
    if (pt > 200)
      tmp_pt = 199.; // Overflow Bin // Mu SF                                                                                                    
  }
  if(id == 2) {
    if (pt > 500)
      tmp_pt = 499.; // Overflow Bin // Pho SF                                                                                            
  }

  int biny = h->GetYaxis()->FindBin(tmp_pt);
  int binx = h->GetXaxis()->FindBin(eta);
  if ( h->GetBinContent(binx,biny) && h->GetBinContent(binx,biny)>0.) return h->GetBinContent(binx,biny);
  else return 1.0;
}

double HtoJPsiGamma::SFCalcSU(TH2F* h, double pt, double eta, int id){
  double tmp_pt = pt;
  if(id == 1) {
    if (pt > 200)
      tmp_pt = 199.; // Overflow Bin // Mu SF                                                                                                 
  }
  if(id == 2) {
    if (pt > 500)
      tmp_pt = 499.; // Overflow Bin // Pho SF                                                                                                    
  }
  int biny = h->GetYaxis()->FindBin(tmp_pt);
  int binx = h->GetXaxis()->FindBin(eta);
  if (h->GetBinError(binx,biny) && h->GetBinError(binx,biny)>0.) return h->GetBinError(binx,biny);
  else return 0.0;
}

double HtoJPsiGamma::triggSFCalc(TH2F* h, double pt_mu, double pt_pho) {
  double tmp_pt_mu = pt_mu;
  double tmp_pt_pho = pt_pho;
  if (pt_mu > 100) tmp_pt_mu = 99;
  if (pt_mu < 20) tmp_pt_mu = 21;
  if (pt_pho > 100) tmp_pt_pho = 99;
  if (pt_pho < 30) tmp_pt_pho = 31;

  int biny = h->GetYaxis()->FindBin(tmp_pt_pho);
  int binx = h->GetXaxis()->FindBin(tmp_pt_mu);
  if ( h->GetBinContent(binx,biny))
    return h->GetBinContent(binx,biny);
  else {
    std::cout<<"Returning 1"<<std::endl;
    return 1.0;
  }
}


double HtoJPsiGamma::triggSFCalcSU(TH2F* h, double pt_mu, double pt_pho) {
  double tmp_pt_mu = pt_mu;
  double tmp_pt_pho = pt_pho;
  if (pt_mu > 100) tmp_pt_mu = 99;
  if (pt_mu < 20) tmp_pt_mu = 21;
  if (pt_pho > 100) tmp_pt_pho = 99;
  if (pt_pho < 30) tmp_pt_pho = 31;
  
  int biny = h->GetYaxis()->FindBin(tmp_pt_pho);
  int binx = h->GetXaxis()->FindBin(tmp_pt_mu);

  if (h->GetBinError(binx,biny))
    {
      return (h->GetBinError(binx,biny));
    }
  else {
    std::cout<<"Returning 0"<<std::endl;
    return 0.;
  }
}






