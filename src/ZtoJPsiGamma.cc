#define ZtoJPsiGamma_cxx
#include "ZtoJPsiGamma.h"
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
#include <chrono>
using namespace std;

void ZtoJPsiGamma::tokenize(const string& str, vector<string>& tokens, const string& delimiters) {
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

void ZtoJPsiGamma::buildList(const vector<string>& tokens, vector<string>& list) {
  for (vector<string>::const_iterator it  = tokens.begin()+1; it != tokens.end(); ++it){
    list.push_back(*it);
  }
}

int ZtoJPsiGamma::setInputFile(const string& fname){
  auto found = fname.find("root:");
  if (found == string::npos && gSystem->AccessPathName(fname.c_str())){
    //    std::cout<<fname<<"setting input inside warning"<<std::endl;
    cerr << ">>> Warning: File <<" << fname << ">> was not found!!" << endl;
    return static_cast<int>(chain->GetEntries());
  }
  // std::cout<<fname<<"setting input outside warning"<<std::endl;

  chain->Add(fname.c_str(), -1);
  std::cout<<fname<<"adding to chain"<<std::endl;
  return static_cast<int>(chain->GetEntries());

}

TH1* ZtoJPsiGamma::getHist1D(const char* hname) {
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

TH1* ZtoJPsiGamma::getHist1D(const string& hname) {
  return getHist1D(hname.c_str());
}



TH2* ZtoJPsiGamma::getHist2D(const char* hname) {
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

double ZtoJPsiGamma::lumiWt(double evtWeightSum, int nent) const
//double ZtoJPsiGamma::lumiWt(double evtWeightSum)
{
  double nevt = (evtWeightSum > -1) ? evtWeightSum : nent;
  std::cout << "-- intLumi: " << intLumi_
            << " xsec: " << xsec_
            << " nevt: " << nevt  <<std::endl;
  
  return (intLumi_ * xsec_ / nevt);

}

// taken from HZZ4lUtil to show cutflow efficiency
void ZtoJPsiGamma::showEfficiency(const string& hname,const std::vector<std::string>& slist,const string& header,const string& tag,std::ostream& os)
{
  os << ">>> " << header << " Efficiency" << endl;
  TH1 *h = ZtoJPsiGamma::getHist1D(hname);
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

void ZtoJPsiGamma::objectEfficiency() {
  // Object selection Efficiency
  // outfile->cd();
  // outfile->cd("ObjectSelection");
  
  // Muon  
  vector<string> muLabels {
    "All",
      "looseid",
      "tightid"
      };
  ZtoJPsiGamma::showEfficiency("muCutFlow", muLabels, "Muon Selection","Muons");

  // Photon
  vector<string> phoLabels {
    "All",
      "mvaid",
      "pt",
      "ScEta"
      };
  ZtoJPsiGamma::showEfficiency("photonCutFlow", phoLabels, "Photon Selection","Photons");
}

void ZtoJPsiGamma::eventEfficiency(){
  //outfile->cd();
  //outfile->cd("ZtoJPsiGamma");
  std::vector<std::string> evlabels {
    "Events processed",
      //      "Gen level filter : Minimum 2 gen Muons & Minimum 1 gen photon",
      //"Gen level filter after detector cuts, eta =< 2.5",
      "Events with > 0 good vertex",
      "Events passing trigger",
      "Atleast 2 muons after muon selection",
      "Relative Isolation of leading mu",
      "leading mu pt",
      "subleading mu pt",
      "Atleast 1 photon after photon selection",
      "deltaR between mu1 & pho",
      "deltaR between mu2 & pho",
      "deltaR between dimu & pho",
      "deltaphi between dimu & pho",
      "JPsi mass cut",
      "Z mass cut",
      "ratio cut for pt of dimu",
      "ratio cut for Et of pho",
      "scale factor"
      };
  
  ZtoJPsiGamma::showEfficiency("evtCutFlow",evlabels, "Event Selection","Eventflow");

  double lumiFac = 1.0;
  if (useLumiWt_)
    lumiFac = lumiWt(evtWeightSum_,200000);
  cout << endl
       << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
       << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
       << endl;

  //  ZtoJPsiGamma::scaleHistogram("evtCutFlowWt", lumiFac);

  TH1 *h = ZtoJPsiGamma::getHist1D("evtCutFlowWt");
  if (h != nullptr) {
    int nbins = h->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
      double cont = h->GetBinContent(i) * lumiFac;
      double err  = h->GetBinError(i) * lumiFac;
      h->SetBinContent(i, cont);
      h->SetBinError(i, err);
    }
  }

  ZtoJPsiGamma::showEfficiency("evtCutFlowWt", evlabels, "Event Selection (Weighted)", "EventFlow");

      std::vector<std::string> catlabels {
      "total final yield",
	"cat1",
	"cat2",
	"cat3"
	//"cat4"
      };

        ZtoJPsiGamma::showEfficiency("nPhoCat",catlabels, "category of events", "category");

  double lumiFac_cat = 1.0;
  if (useLumiWt_) 
    lumiFac_cat = lumiWt(evtWeightSum_, 200000);
  cout << endl
       << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
       << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac_cat
       << endl;
  
  TH1 *hcat = ZtoJPsiGamma::getHist1D("nPhoCatWt");
  if (hcat != nullptr) {
    int nbins_cat = hcat->GetNbinsX();
    for (int j = 1; j <= nbins_cat; ++j) {
      double cont_cat = hcat->GetBinContent(j) * lumiFac_cat;
      double err_cat  = hcat->GetBinError(j) * lumiFac_cat;
      hcat->SetBinContent(j, cont_cat);
      hcat->SetBinError(j, err_cat);
    }
  }

  ZtoJPsiGamma::showEfficiency("nPhoCatWt",catlabels, "category of weighted events", "Category");

}
      
void ZtoJPsiGamma::showCategoryYield(){
  // just for test, please ignore

}

// Helicity angle between J/Psi and mu+ 
float ZtoJPsiGamma::costheta1(TLorentzVector particle_orig, TLorentzVector parent_orig){
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

float ZtoJPsiGamma::costhetastar(TLorentzVector particle_orig, TLorentzVector parent_orig){
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

void ZtoJPsiGamma::bookHistograms() {
  //TString of =  "ZtoJPsiGamma_hist.root";
  //outfile = new TFile(of,"RECREATE");
  outfile->cd();
  outfile->mkdir("ObjectSelection");
  outfile->cd("ObjectSelection");
  new TH1D("muCutFlow", "Muon Cut Flow", 3, -0.5, 2.5);
  new TH1D("photonCutFlow", "Photon Cut Flow", 4, -0.5, 3.5);
  
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
  new TH1F("phoPt", "Photon p_{T}", 400, 0., 200.);
  new TH1F("Photon_isScEtaEB", "Photon isScEtaEB", 2, 0., 1.);
  new TH1F("Photon_isScEtaEE", "Photon isScEtaEE", 2, 0., 1.);
  new TH1D("nPhoScEta", "#Photons in ScEta regions", 2, -0.5, 1.5);
  new TH1D("phoR9", "Photon R9", 40, 0.4, 1.0);
  new TH1D("phoR9_ScEB", "Photon R9 for ScEta in EB", 40, 0.4, 1.0);
  new TH1D("phoR9_ScEE", "Photon R9 for ScEta in EE", 40, 0.4, 1.0);

  outfile->cd();
  outfile->mkdir("ZtoJPsiGamma");
  outfile->cd("ZtoJPsiGamma");
  new TH1F("generatorweight", "MC generator weight", 3, -1.5, 1.5);
  new TH1F("puweight", "PU weight", 200, 0., 2.);
  new TH2F("pudist","Pu weight distribution", 80, 0., 80., 100., 0., 2. );
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


  new TH1F("jpsimass_bmasscut","jpsimass before mass cut", 60, 2.8, 3.4);
  new TH1F("jpsimass_amasscut","jpsimass after mass cut", 60, 2.8, 3.4);
  new TH1F("jpsimass_afinalcut","jpsimass after final cut", 60, 2.8, 3.4);
  new TH1F("zmass_bmasscut","zmass before Z mass cut", 64, 40., 200.);
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
  

  // To select which cut to choose
  //  new TH1F("dRMu1Cut", "# of Events for different dR(mu1, #gamma)", )
    
  bookedHistograms = true;
  
}

void ZtoJPsiGamma::saveHistograms() {
  cout << "About to close the histogram file" << endl;
  if (bookedHistograms) {
    //    outfile->cd();
    //outfile->cd("ZtoJPsiGamma");
    outfile->Write();  
    //    outfile->Close();
  }
  cout << "histogram file not closed" << endl;

}
// ***********************************************


bool ZtoJPsiGamma::readJob(const std::string& jobFile, int& nFiles){
  
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
    ZtoJPsiGamma::tokenize(line, tokens);
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
    else if (key == "useLumiWt")
      useLumiWt_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "intLumi")
      intLumi_ = std::stof(value.c_str());
    else if (key == "xsec")
      xsec_ = std::stof(value.c_str());
    else if (key == "nevents")
      nevents_ = std::stoi(value.c_str()); 
    else if (key == "year") 
      year_ = std::stoi(value.c_str());
    else if (key == "isPre") {
      string vtmp(value);
      isPre_ = (vtmp == "True") ? true:false;
    }
    else if (key == "histFile")
      histFile_ = value;
    else if (key == "hlt")
      hlt_ = std::stoi(value.c_str());
    else if (key == "inputFile") {
      ZtoJPsiGamma::buildList(tokens, fileList_);
    }
  }
  
  // Close the file                                                  
  fin.close();
  // Build the chain of root files
  for (const auto& fname: fileList_) {
    cout << ">>> INFO. Adding input file " << fname << " to TChain " << endl;
    ++nFiles;
    std::cout<<"Inside"<<std::endl;
    int nevt = setInputFile(fname);
    //    if (maxEvt_ > 0 && nevt >= maxEvt_) break;
  }

  if (!nFiles) {
    cerr << ">>> WARN. Input Root file list is empty! exiting ..." << endl;
    return false;
  }

  return true;



}

// ************************************

void ZtoJPsiGamma::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ZtoJPsiGamma.C
//      root> ZtoJPsiGamma t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    chain->SetBranchStatus("*",0);  // disable all branches
//    chain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    chain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  Float_t zm1,zm2,zm3, mu1pt, mu2pt, phopt, muiso, drmu1pho, drmu2pho, drjpsipho, dphijpsipho, jpsim, Hm, ptratio, etratio, finalm;


   if (chain == 0) return;
   //Long64_t nentries = chain->GetEntriesFast();
   int nentries = static_cast<int>(chain->GetEntries());

   //   if (maxEvt_ > 0) nentries = std::min(nentries, maxEvt_);
   //std::cout << "max event: " << maxEvt_ << std::endl;
   std::cout << "Total number of entries: " << nentries << std::endl;
   
   outfile = TFile::Open(histFile_.c_str(), "RECREATE");
   
   //Float_t finaljpsim,finalm, finalm1,finalm2,finalm3;
   TTree *T = new TTree("T","test final mass tree");
   TTree *T1 = new TTree("T1","cat1 mass tree");
   TTree *T2 = new TTree("T2","cat2 mass tree");
   TTree *T3 = new TTree("T3","cat3 mass tree");

   T1->Branch("zm1",&zm1,"zm1/F");
   T2->Branch("zm2",&zm2,"zm2/F");
   T3->Branch("zm3",&zm3,"zm3/F");

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


   bookHistograms();
   //TFile *mu_scale_factors;
   TFile *pho_scale_factors;

   if (year() == 16) {
     //mu_scale_factors = new TFile("../SF/final_HZZ_muon_SF_2016RunB2H_legacy_newLoose_newIso_paper.root");// Twiki      
     if(isPre()) {
       //       std::cout<<"gotPre"<<std::endl;
       pho_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/H_Z_toJPsiGamma/SF/egammaEffi.txt_EGM2D_Pho_wp90_UL16.root");
     }
     else   pho_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/H_Z_toJPsiGamma/SF/egammaEffi.txt_EGM2D_Pho_MVA90_UL16_postVFP.root");
   }
   
   if(year() == 17) {
     // mu_scale_factors = new TFile("../SF/final_HZZ_muon_SF_2017_newLooseIso_mupogSysts_paper.root");
     pho_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/H_Z_toJPsiGamma/SF/egammaEffi.txt_EGM2D_PHO_MVA90_UL17.root");
   }
   
   if(year() == 18) {
     //  mu_scale_factors = new TFile("../SF/final_HZZ_muon_SF_2018RunA2D_ER_newLoose_newIso_paper.root");
     pho_scale_factors = new TFile("/afs/cern.ch/work/s/sbaradia/public/Analysis/H_Z_toJPsiGamma/SF/egammaEffi.txt_EGM2D_Pho_wp90.root_UL18.root");
   }
   TH2F *pho_scale;
   // TH2F *mu_scale = (TH2F*)mu_scale_factors->Get("FINAL");
   if(isMC()) pho_scale = (TH2F*)pho_scale_factors->Get("EGamma_SF2D");

     int num_pu_vtx = 0;

   Long64_t nbytes = 0, nb = 0;
   std::cout<<"No. of entries is"<<nentries<<std::endl;   
   auto t_start = std::chrono::high_resolution_clock::now();
   auto t_end = std::chrono::high_resolution_clock::now();

   //   for (Long64_t jentry=0; jentry<5000; jentry++) {
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = chain->GetEntry(jentry);   nbytes += nb;
     
     if(jentry%100000 == 0 )
       {
         t_end = std::chrono::high_resolution_clock::now();
	 std::cout<<"Processing Entry in event loop (Rate) : "<<jentry<<" / "<<nentries<<"  [ "<<100.0*jentry/nentries<<"  % ]  "
                  << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                  <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries - jentry)/(1e-9 + jentry)* 0.001
                  <<std::endl;
       }

     float weight_scalefactor_mu = 1.;
     float weight_scalefactor_pho = 1.;


     std::vector<float> vec_genmupt;
     std::vector<float> vec_genphopt;
     std::vector<float> vec_genmueta;
     std::vector<float> vec_genphoeta;
     TLorentzVector genmuP4,  genphoP4;
     std::vector<TLorentzVector> vec_genmuP4;
     std::vector<TLorentzVector> vec_genphoP4;
          
     double wt = 1.0;
     num_pu_vtx = Pileup_nPU;
     if (isMC()) wt = wt * Generator_weight * puWeight;
     else wt = 1.0;
     double evtweight = -99.;

     if (isMC() && Generator_weight != 1.0) evtweight = Generator_weight;
     else evtweight = 1.0;
     
     evtWeightSum_ += evtweight;
     

     outfile->cd();
     outfile->cd("ZtoJPsiGamma");
     
     fillHist1D("generatorweight", Generator_weight);
     fillHist1D("puweight", puWeight);
     fillHist2D("pudist", num_pu_vtx, puWeight);
     fillHist1D("evtCutFlow", 0);
     fillHist1D("evtCutFlowWt", 0, wt);
     fillHist1D("nPV", PV_npvs, wt);
     fillHist1D("nPVgood", PV_npvsGood, wt);
     if (PV_npvsGood < 1) continue;
     fillHist1D("evtCutFlow", 1);
     fillHist1D("evtCutFlowWt", 1, wt);

     bool HLT;
     //if (hlt() == 0) HLT = HLT_Mu17_Photon30_CaloIdL_L1ISO;
     if (hlt() == 0) HLT = HLT_Mu17_Photon30_IsoCaloId;
     else if (hlt() == 1) HLT = HLT_Dimuon25_Jpsi;
     else if (hlt() == 2) HLT = HLT_DoubleMu20_7_Mass0to30_Photon23;
     else if (hlt() == 3) HLT = HLT_Dimuon20_Jpsi;
     else if (hlt() == 4) HLT = HLT_Mu17_Photon30_CaloIdL_L1ISO;
     else std::cout << "HLT not found" << std::endl;


     if (!HLT) continue;
     fillHist1D("evtCutFlow", 2);
     fillHist1D("evtCutFlowWt", 2, wt);
     

     //     if (nMuon < 2 || nPhoton < 1 ) continue; // Just for random grid search

     //std::vector<float> vec_mupt;
       //std::vector<float> vec_phopt;
     TLorentzVector muLP4, muTP4, muTIsoP4, phoP4;
     std::vector<TLorentzVector> vec_muLP4, vec_muTP4, vec_muTIsoP4;
     std::vector<int> vec_muid;
     std::vector<TLorentzVector> vec_phoP4;
     
     outfile->cd();
     outfile->cd("ObjectSelection");
     float mupt_max = -999.;
     int mu1id = -999;
     for (int i=0; i < nMuon; i++){
       fillHist1D("muCutFlow", 0);
       
       fillHist1D("muPt", Muon_pt[i], wt);
       fillHist1D("muEta", Muon_eta[i], wt);
       fillHist1D("muDxy", Muon_dxy[i], wt);
       fillHist1D("muDz", Muon_dz[i], wt);
       fillHist1D("muSIP3d", Muon_sip3d[i], wt);
       fillHist1D("muSIP3dincm", Muon_ip3d[i], wt);
       if (!Muon_looseId[i]) continue; // Muon POG Loose ID
       fillHist1D("muCutFlow", 1);
       muLP4.SetPtEtaPhiM( Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
       vec_muLP4.push_back(muLP4);

       //std::cout << "Muon_highPtId : " << Muon_highPtId[i] << std::endl; 
       
       //if (!(Muon_isPFcand[i] || (Muon_highPtId[i] == 1 && Muon_pt[i] > 200))) continue; // HZZ4L Tight ID
       if (!Muon_tightId[i]) continue;
       fillHist1D("muCutFlow", 2);
       
       muTP4.SetPtEtaPhiM( Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
       vec_muTP4.push_back(muTP4);

       if (Muon_pt[i] > mupt_max){
	 mupt_max = Muon_pt[i];
	 mu1id = i;
       }
       
     }

     int mu2id = -999.;
     int mu2ndmaxpt = -999.;
     if (nMuon > 0) {
       for (int j=0; j < nMuon; j++){
	 if (j==mu1id) continue;
	 if (Muon_pt[j] > mu2ndmaxpt){
	   mu2ndmaxpt = Muon_pt[j];
	   mu2id = j;
	 }
       }
     }
     outfile->cd();
     outfile->cd("ObjectSelection");     
     
     float phopt_max = -999.;
     int pho1id = -999;     

     for (int i = 0 ; i < nPhoton ; i++){
       fillHist1D("photonCutFlow", 0);

       //fillHist1D("phoMvaidwp90", Photon_mvaID_WP90[i], wt);     
       if (Photon_mvaID_WP90[i] == false) continue;
       fillHist1D("photonCutFlow", 1);
       
       fillHist1D("phoPt", Photon_pt[i], wt);

       if (Photon_pt[i] <= 33) continue;
       phopt = Photon_pt[i];
       fillHist1D("photonCutFlow", 2);
       
       if (Photon_isScEtaEB[i]) fillHist1D("nPhoScEta", 0, wt); // 0 = ScEtaEB, 1 = ScEtaEE
       else if (Photon_isScEtaEE[i]) fillHist1D("nPhoScEta", 1, wt);

       //       fillHist1D("Photon_isScEtaEB", Photon_isScEtaEB[i], wt);
       //fillHist1D("Photon_isScEtaEE", Photon_isScEtaEE[i], wt);
       if (! (Photon_isScEtaEB[i] || Photon_isScEtaEE[i])) continue;
       fillHist1D("photonCutFlow", 3);

       fillHist1D("phoR9", Photon_r9[i], wt);
       if (Photon_isScEtaEB[i]) fillHist1D("phoR9_ScEB", Photon_r9[i], wt);
       if (Photon_isScEtaEE[i]) fillHist1D("phoR9_ScEE", Photon_r9[i], wt);

       phoP4.SetPtEtaPhiM(Photon_pt[i], Photon_eta[i], Photon_phi[i], Photon_mass[i]);
       vec_phoP4.push_back(phoP4);
       
       if (Photon_pt[i] > phopt_max){
	 phopt_max = Photon_pt[i];
	 pho1id = i;
       }


     }
     //     std::cout << "# photons  : " << vec_phoP4.size() << std::endl;

     outfile->cd();
     outfile->cd("ZtoJPsiGamma");
     
     //if (vec_muTP4.size()>0) std::cout << "lead mu pt from P4 : " << vec_muTP4[0].Pt() << "lead mu pt from index : " << Muon_pt[mu1id] << std::endl;
     fillHist1D("nLoosemuons", vec_muLP4.size(), wt);
     fillHist1D("nTightmuons", vec_muTP4.size(), wt);
     if ( vec_muTP4.size() < 2 ) continue;
     fillHist1D("evtCutFlow", 3);
     fillHist1D("evtCutFlowWt", 3, wt);
     
     std::sort(vec_muTP4.begin(), vec_muTP4.end(), PtComparatorTL<TLorentzVector>());
     fillHist1D("muPfRelIso03all", Muon_pfRelIso03_all[mu1id], wt);
     fillHist1D("muPfRelIso03chg", Muon_pfRelIso03_chg[mu1id], wt);
     fillHist1D("muPfRelIso04all", Muon_pfRelIso04_all[mu1id], wt);

     fillHist1D("muPfRelIso03allByMuSubleadpt", Muon_pfRelIso03_all[mu1id] / vec_muTP4[1].Pt(), wt);
     if (Muon_pfRelIso03_all[mu1id] >= 0.35) continue;
     muiso = Muon_pfRelIso03_all[mu1id];
     fillHist1D("evtCutFlow", 4);
     fillHist1D("evtCutFlowWt", 4, wt);
     
     fillHist1D("muLeadpt", vec_muTP4[0].Pt(), wt);
     //     fillHist1D("muLeadpt", vec_muTP4[0].Pt(), wt); //test
 
     if (vec_muTP4[0].Pt() <= 20) continue;
     mu1pt = vec_muTP4[0].Pt();
     fillHist1D("evtCutFlow", 5);
     fillHist1D("evtCutFlowWt", 5, wt);

     fillHist1D("muSubleadpt", vec_muTP4[1].Pt(), wt);     
     if (vec_muTP4[1].Pt() <= 4) continue;
     mu2pt = vec_muTP4[1].Pt();
     fillHist1D("evtCutFlow", 6);
     fillHist1D("evtCutFlowWt", 6, wt);


     fillHist1D("nPhotons", vec_phoP4.size(), wt);
     if ( vec_phoP4.size() < 1) continue;
     fillHist1D("evtCutFlow", 7);
     fillHist1D("evtCutFlowWt", 7, wt);
     std::sort(vec_phoP4.begin(), vec_phoP4.end(), PtComparatorTL<TLorentzVector>());

     outfile->cd();
     outfile->cd("ZtoJPsiGamma");

     fillHist1D("dRmu1pho", vec_muTP4[0].DeltaR(vec_phoP4[0]), wt);
     if (vec_muTP4[0].DeltaR(vec_phoP4[0]) <= 1 ) continue;
     drmu1pho = vec_muTP4[0].DeltaR(vec_phoP4[0]);
     fillHist1D("evtCutFlow", 8);
     fillHist1D("evtCutFlowWt", 8, wt);

     fillHist1D("dRmu2pho", vec_muTP4[1].DeltaR(vec_phoP4[0]), wt);
     if (vec_muTP4[1].DeltaR(vec_phoP4[0]) <= 1 ) continue;
     drmu2pho = vec_muTP4[1].DeltaR(vec_phoP4[0]);
     fillHist1D("evtCutFlow", 9);     
     fillHist1D("evtCutFlowWt", 9, wt);

     fillHist1D("dRdimupho", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), wt);
     if  ((vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]) <= 2 ) continue;
     drjpsipho = (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]);
     fillHist1D("evtCutFlow", 10);  
     fillHist1D("evtCutFlowWt", 10, wt);

     fillHist1D("dPhidimupho", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), wt);
     if (std::fabs((vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0])) <= 1.5 ) continue;
     dphijpsipho = (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]);
     fillHist1D("evtCutFlow", 11);
     fillHist1D("evtCutFlowWt", 11, wt);

     fillHist1D("jpsimass_bmasscut", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);
     
     if ( (vec_muTP4[0] + vec_muTP4[1]).M() < 3.0 || (vec_muTP4[0] + vec_muTP4[1]).M() > 3.2) continue;
     jpsim = (vec_muTP4[0] + vec_muTP4[1]).M(); 
     fillHist1D("evtCutFlow", 12);
     fillHist1D("evtCutFlowWt", 12, wt);

     fillHist1D("jpsimass_amasscut", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);     

     fillHist1D("zmass_bmasscut", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt );
     
     if ((vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() <= 70.0 || (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() >= 120.0) continue;
     //if((vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() <= ZtoJPsiGamma::cutValue(SRCutMap(), "low") || (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() >= ZtoJPsiGamma::cutValue(SRCutMap(), "high")) continue;
     Hm = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
     fillHist1D("evtCutFlow", 13);
     fillHist1D("evtCutFlowWt", 13, wt);

     fillHist1D("zmass_amasscut", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt );

     fillHist1D("ratiocutdimu", (vec_muTP4[0] + vec_muTP4[1]).Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
     if ((vec_muTP4[0] + vec_muTP4[1]).Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() <= 0.384) continue;//
     //     if ((vec_muTP4[0] + vec_muTP4[1]).Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() <= ZtoJPsiGamma::cutValue(ratioCutMap(), "ptratio")) continue;
     ptratio = (vec_muTP4[0] + vec_muTP4[1]).Pt()/ (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
     fillHist1D("evtCutFlow", 14);
     fillHist1D("evtCutFlowWt", 14, wt);

     fillHist1D("ratiocutpho", vec_phoP4[0].Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
     if (vec_phoP4[0].Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() <= 0.384) continue; 
     //     if (vec_phoP4[0].Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() <= ZtoJPsiGamma::cutValue(ratioCutMap(), "etratio")) continue;
     etratio = vec_phoP4[0].Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
     fillHist1D("evtCutFlow", 15);          
     fillHist1D("evtCutFlowWt", 15, wt);

     //     if (vec_phoP4[0].Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() > ZtoJPsiGamma::cutValue(ratioCutMap(), "etratio")) { fillHist1D("nPhoCat", 0); fillHist1D("nPhoCatWt", 0, wt);}

     //     if (vec_phoP4[0].Pt()/(vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M() > 0.384) {
     // fillHist1D("nPhoCat", 0); fillHist1D("nPhoCatWt", 0, wt);
     // }

     if(isMC()) {
       //     weight_scalefactor_mu *= SFCalc(mu_scale, vec_muTP4[0].Pt(), vec_muTP4[0].Eta() ) * SFCalc(mu_scale, vec_muTP4[1].Pt(), vec_muTP4[1].Eta() );                                
       weight_scalefactor_pho *= ZtoJPsiGamma::SFCalc(pho_scale, vec_phoP4[0].Pt(), vec_phoP4[0].Eta());                          
       
       //     std::cout << "weight_scalefactor mu = " << weight_scalefactor_mu <<"weight_scalefactor pho"<<weight_scalefactor_pho<< "\t" << "Wt before scalefactor = " << wt << std::endl;
       
       if (weight_scalefactor_mu > 0. &&  weight_scalefactor_pho > 0.) wt *= weight_scalefactor_mu * weight_scalefactor_pho;
       //       std::cout<< "Wt after scalefactor = " << wt << std::endl;
     }
     fillHist1D("nPhoCat", 0); fillHist1D("nPhoCatWt", 0, wt);
     
     fillHist1D("evtCutFlow", 16);
     fillHist1D("evtCutFlowWt", 16, wt);
     
     if (Photon_isScEtaEB[pho1id] && Photon_r9[pho1id] >= 0.94 )   { 
       fillHist1D("nPhoCat", 1); 
       fillHist1D("nPhoCatWt", 1, wt);
       fillHist1D("muLeadpt_cat1", vec_muTP4[0].Pt(), wt);
       fillHist1D("muLeadeta_cat1", vec_muTP4[0].Eta(), wt);
       fillHist1D("muSubleadpt_cat1", vec_muTP4[1].Pt(), wt);
       fillHist1D("muSubleadeta_cat1", vec_muTP4[1].Eta(), wt);
       fillHist1D("phoPt_cat1", vec_phoP4[0].Pt(), wt);
       fillHist1D("phoEta_cat1", vec_phoP4[0].Eta(), wt);
       fillHist1D("dRmu1pho_cat1", vec_muTP4[0].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRmu2pho_cat1", vec_muTP4[1].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRdimupho_cat1", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dPhidimupho_cat1", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), wt);
       fillHist1D("jpsimass_cat1", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);
       fillHist1D("dRmu1mu2_cat1", vec_muTP4[0].DeltaR(vec_muTP4[1]), wt); 
       fillHist1D("diMuPt_cat1", (vec_muTP4[0] + vec_muTP4[1]).Pt(), wt);
       fillHist1D("finalZPt_cat1", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).Pt(), wt);
       fillHist1D("zmass_cat1", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
       fillHist1D("vtxdist_x" , PV_x - (vec_muTP4[0] + vec_muTP4[1]).X(), wt);
       fillHist1D("vtxdist_y" , PV_y - (vec_muTP4[0] + vec_muTP4[1]).Y(), wt);
       fillHist1D("vtxdist_z" , PV_z - (vec_muTP4[0] + vec_muTP4[1]).Z(), wt);
       zm1 = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
       T1->Fill();
     }
     else if (Photon_isScEtaEB[pho1id] && Photon_r9[pho1id] < 0.94 ){
       fillHist1D("nPhoCat", 2); 
       fillHist1D("nPhoCatWt", 2, wt); 
       fillHist1D("muLeadpt_cat2", vec_muTP4[0].Pt(), wt);
       fillHist1D("muLeadeta_cat2", vec_muTP4[0].Eta(), wt);
       fillHist1D("muSubleadpt_cat2", vec_muTP4[1].Pt(), wt);
       fillHist1D("muSubleadeta_cat2", vec_muTP4[1].Eta(), wt);
       fillHist1D("phoPt_cat2", vec_phoP4[0].Pt(), wt);
       fillHist1D("phoEta_cat2", vec_phoP4[0].Eta(), wt);
       fillHist1D("dRmu1pho_cat2", vec_muTP4[0].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRmu2pho_cat2", vec_muTP4[1].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRdimupho_cat2", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dPhidimupho_cat2", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), wt);
       fillHist1D("jpsimass_cat2", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);
       fillHist1D("dRmu1mu2_cat2", vec_muTP4[0].DeltaR(vec_muTP4[1]), wt);
       fillHist1D("diMuPt_cat2", (vec_muTP4[0] + vec_muTP4[1]).Pt(), wt);
       fillHist1D("finalZPt_cat2", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).Pt(), wt);
       fillHist1D("zmass_cat2", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
       zm2 = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
       T2->Fill();
     } 
     else if (Photon_isScEtaEE[pho1id] ) {
       fillHist1D("nPhoCat", 3); 
       fillHist1D("nPhoCatWt", 3, wt);
       fillHist1D("muLeadpt_cat3", vec_muTP4[0].Pt(), wt);
       fillHist1D("muLeadeta_cat3", vec_muTP4[0].Eta(), wt);
       fillHist1D("muSubleadpt_cat3", vec_muTP4[1].Pt(), wt);
       fillHist1D("muSubleadeta_cat3", vec_muTP4[1].Eta(), wt);
       fillHist1D("phoPt_cat3", vec_phoP4[0].Pt(), wt);
       fillHist1D("phoEta_cat3", vec_phoP4[0].Eta(), wt);
       fillHist1D("dRmu1pho_cat3", vec_muTP4[0].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRmu2pho_cat3", vec_muTP4[1].DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dRdimupho_cat3", (vec_muTP4[0] + vec_muTP4[1]).DeltaR(vec_phoP4[0]), wt);
       fillHist1D("dPhidimupho_cat3", (vec_muTP4[0] + vec_muTP4[1]).DeltaPhi(vec_phoP4[0]), wt);
       fillHist1D("jpsimass_cat3", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);
       fillHist1D("dRmu1mu2_cat3", vec_muTP4[0].DeltaR(vec_muTP4[1]), wt);
       fillHist1D("diMuPt_cat3", (vec_muTP4[0] + vec_muTP4[1]).Pt(), wt);
       fillHist1D("finalZPt_cat3", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).Pt(), wt);
       fillHist1D("zmass_cat3", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
       zm3 = (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M();
       T3->Fill();
     }
     //final plots
     
     T->Fill();
     
     fillHist1D("jpsimass_afinalcut", (vec_muTP4[0] + vec_muTP4[1]).M(), wt);
     fillHist1D("zmass_afinalcut", (vec_muTP4[0] + vec_muTP4[1] + vec_phoP4[0]).M(), wt);
   
   }
   saveHistograms();
   outfile->cd();
   outfile->cd("ObjectSelection");
   objectEfficiency();
   showCategoryYield();
   outfile->cd();
   outfile->cd("ZtoJPsiGamma");
   T->Write();
   T1->Write();
   T2->Write();
   T3->Write();
   eventEfficiency();
   outfile->Close();
   cout << "Now the ouput file is closed" << endl;
}


double ZtoJPsiGamma::SFCalc(TH2F* h, double pt, double eta){
  double tmp_pt = pt;
  if (pt > 200) tmp_pt = 199.; // Overflow Bin // Mu SF                                                                                                                                              
  int biny = h->GetYaxis()->FindBin(tmp_pt);
  int binx = h->GetXaxis()->FindBin(eta);
  if ( h->GetBinContent(binx,biny) && h->GetBinContent(binx,biny)>0.) return h->GetBinContent(binx,biny);
  else return 1.0;
}




