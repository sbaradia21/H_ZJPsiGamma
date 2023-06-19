#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <utility>
#include <typeinfo>
#include <memory>
#include <cmath>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include "TSpline.h"
#include "TString.h"

#include "configana.h"
#include "HZZ4lUtil.h"
#include "GenAnalysis.h"
#include "DiHiggs4l.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;

using namespace vhtm;
using ZSpace::ZCandidate;
using ZSpace::ZType;

// -----------
// Constructor
// -----------
DiHiggs4l::DiHiggs4l()
  : PhysicsObjSelector()
{
}
// ----------
// Destructor
// ----------
DiHiggs4l::~DiHiggs4l()
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool DiHiggs4l::beginJob()
{
  // Open the output ROOT file (in AnaBase)
  if (!PhysicsObjSelector::beginJob()) return false;

  histf()->cd();
  histf()->mkdir("DiHiggs4l");
  
  bookHistograms();
  if (readGenInfo() && genAna_ != nullptr) {
    genAna_ = std::make_unique<GenAnalysis>();
    genAna_->bookHistograms(histf());
  }
  
  // Optionally write selected events to a flat tree for MVA studies
  // or prepare for MVA Analysis
  
  
  
  if (createMVATree_)
    skimObj_ = std::make_unique<MVASkim>(mvaSkimFile_);
  else if (readMVA_ && (isSig_ == 1) ){
    mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
    
    /*    if (sigtype_ == "grav") {
      mvaXMLFile_ = "dataset" + std::to_string(masspoint_) + "/weights/TMVAClassification_BDTG.weights.xml";
      if (masspoint_ == 260) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 300) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 400) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 500) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 600) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 700) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 750) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 800) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 900) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 1000) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
    }
    else if (sigtype_ == "rad") {
      mvaXMLFile_ = "datasetR" + std::to_string(masspoint_) + "/weights/TMVAClassification_BDTG.weights.xml";
      if (masspoint_ == 260) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 300) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 400) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 500) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 600) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 700) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 750) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 800) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 900) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      else if (masspoint_ == 1000) mvAna_   = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFile_);
      }*/
  }
  else if (readMVA_ && (isSig_ ==0) ){
    if (sigtype_ == "grav"){
      mvaXMLFileBkg260_ ="dataset260/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg300_ ="dataset300/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg400_ ="dataset400/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg500_ ="dataset500/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg600_ ="dataset600/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg700_ ="dataset700/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg750_ ="dataset750/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg800_ ="dataset800/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg900_ ="dataset900/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg1000_ ="dataset1000/weights/TMVAClassification_BDTG.weights.xml";
    }
    else if (sigtype_ == "rad"){
      mvaXMLFileBkg260_ ="dataset260R/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg300_ ="dataset300R/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg400_ ="dataset400R/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg500_ ="dataset500R/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg600_ ="dataset600R/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg700_ ="dataset700R/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg750_ ="dataset750R/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg800_ ="dataset800R/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg900_ ="dataset900R/weights/TMVAClassification_BDTG.weights.xml";
      mvaXMLFileBkg1000_ ="dataset1000R/weights/TMVAClassification_BDTG.weights.xml";
    }

    std::cout << "mvaXMLFileBkg1000_ : " << mvaXMLFileBkg1000_ << std::endl;

      mvAnaBkg260_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFileBkg260_);
      mvAnaBkg300_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFileBkg300_);
      mvAnaBkg400_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFileBkg400_);    
      mvAnaBkg500_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFileBkg500_);    
      mvAnaBkg600_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFileBkg600_); 
      mvAnaBkg700_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFileBkg700_);    
      mvAnaBkg750_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFileBkg750_);    
      mvAnaBkg800_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFileBkg800_);    
      mvAnaBkg900_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFileBkg900_);
      mvAnaBkg1000_ = std::make_unique<MVAnalysis>(mvAlgo_, mvaXMLFileBkg1000_);
      
      if (mvAnaBkg1000_ == nullptr)  std::cout << "mvAnaBkg1000 is null "<< std::endl;
      else std::cout << "mvAnaBkg1000 is not null "<< std::endl;
    
  }

#ifdef SKIP_DUPLICATE_ALL
  eventIdStore_.clear();
#endif

  return true;
  //fin} // don't know why this line is here !!!
}
// ---------------
// Book histograms
// ---------------
void DiHiggs4l::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();

  histf()->cd();
  histf()->cd("DiHiggs4l");

  // book histograms to be filled at different stages
  new TH1D("nvtx_stage0", "Number of Good vertices (stage 0)", 60, 0, 60);
  new TH1D("nvtx_stage1", "Number of Good vertices (stage 1)", 60, 0, 60);
  new TH1D("nvtx_stage2", "Number of Good vertices (stage 2)", 60, 0, 60);
  new TH1D("nvtx_stage3", "Number of Good vertices (stage 3)", 60, 0, 60);
  new TH1D("nvtx_stage4", "Number of Good vertices (stage 4)", 60, 0, 60);
  new TH1D("nvtx_stage5", "Number of Good vertices (stage 5)", 60, 0, 60);
  new TH1D("nvtx_stage6", "Number of Good vertices (stage 6)", 60, 0, 60);
  new TH1D("nvtx_stage7", "Number of Good vertices (stage 7)", 60, 0, 60);
  if (isMC()) {
    new TH1D("puweight", "PU reweight factor", 100, 0, 2);
    new TH1D("evtweight", "Event weight factor (MC events)", 100, -5, 5);
  }

  // ------- Object PLots -----------------------------------------------
  new TH1D("nGoodMuon", "Number of Good muons", 6, -0.5, 5.5);
  new TH1D("nGoodEle",  "Number of Good electrons", 6, -0.5, 5.5);
  new TH1D("nJets",     "Number of Tight Jets", 10, -0.5, 9.5);
  new TH1D("nbJets",    "Number of b-Jets", 10, -0.5, 9.5);

  new TH1D("zzCutFlow", "ZZ selection CutFlow", 8, -0.5, 7.5);
  new TH1D("evtCutFlow", "Event CutFlow", 17, -0.5, 16.5);
  if (isMC())
    new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 17, -0.5, 16.5);

  new TH1D("isTriggered", "Event triggered", 2, -0.5, 1.5);

  new TH1D("nLooseLeptons",    "Total tight isolated leptons", 10, -0.5, 9.5);
  new TH1D("nTightIsoLeptons", "Total tight isolated leptons", 8, -0.5, 7.5);
  new TH1D("hLepPt_stage0",    "Highest pT lepton pT (stage0)", 300, 0, 300);
  new TH1D("hLepPt",           "Highest pT lepton pT", 300, 0, 300);
  new TH1D("met",              "Missing Transverse Energy", 150, 0, 150);
  new TH1D("corrmet",          "Missing Transverse Energy (Type-I Corrected)", 150, 0, 150);
  new TH1D("puppimet",         "Missing Transverse Energy (PUPPI)", 150, 0, 150);

  // jet properties
  new TH1F("bjetDiscriminator_alljets", "bjet discriminator for all jets", 100, 0., 1.);  
  new TH1F("bjetDiscriminator_preselectedjets", "bjet discriminator for preselcted jets(#p_{T} > 20, |#eta|< 2.4)", 100, 0., 1.);
  new TH1F("bjetDiscriminator_cleanjets", "bjet discriminator for cleaned selected jets(#p_{T} > 20, |#eta|< 2.4, #Delta R(lep, j) > 0.3)", 100, 0., 1.);

  // Z Candidates
  new TH1F("nZcand", "Number of selected Z candidates", 10, 0, 10);
  new TH1F("massZcand", "Mass of selected Z candidates", 100, 0., 200.);

  new TH1F("Z1mass", "ll invariant mass of the Z1 candidate", 100, 20. , 120.);
  new TH1F("Z2mass", "ll invariant mass of the Z2 candidate", 30, 0., 200.);
  new TH1F("fourLmass", "Invariant mass of the 4-lepton candidate", 200, 0., 200.);
  new TH1F("mass2l", "Invariant mass of di-lepton candidate (CR) ", 200, 0., 200.);
   

  // 4-lepton properties---------------------------------------------------------------------------------------------------------------------------------------------
  new TH1F("mass4l",  "4l invariant mass(all)", 50, 0., 200.);
  new TH1F("mass4l_2j",  "4l invariant mass(all) (4l plus atleast 2jets)", 120, 80., 140.);
  new TH1F("mass4l_2j_mHc",  "4l invariant mass(al-l) (4l plus atleast 2jets plus |m4l-mH| <= 10)", 60, 110., 140.);
  new TH1F("mass4l_2j_mHc_1b",  "4l invariant mas-s(all) (4l plus atleast 2jets plus |m4l-mH| <= 10 plus 1bjet)", 60, 110., 140.);

  // Resonant mass 
  new TH1F("resomass",  "Invariant mass of Resonant particle after 4l plus atleast 2jets plus |m4l-mH| <= 10", 140, 100., 1500.);
  new TH1F("resopt",  "Pt of Resonant particle after 4l plus atleast 2jets plus |m4l-mH| <= 10", 140, 100., 1500.);
  new TH1F("resomass_csv",  "Invariant mass of Resonant particle after 4l plus atleast 2jets plus |m4l-mH| <= 10 for highest b-tagged jets", 140, 100., 1500.);
  new TH1F("resopt_csv",  "Pt of Resonant particle after 4l plus atleast 2jets plus |m4l-mH| <= 10 for highest b-tagged jets", 140, 100., 1500.);

  // Lepton properties (pT, eta etc.)
  new TH1F("Z1l1Pt", "pT of the first lepton from Z1", 300, 0, 300);
  new TH1F("Z1l2Pt", "pT of the second lepton from Z1", 300, 0, 300);
  new TH1F("Z2l1Pt", "pT of the first lepton from Z2", 300, 0, 300);
  new TH1F("Z2l2Pt", "pT of the second lepton from Z2", 300, 0, 300);

  new TH1F("Z1l1Eta", "Eta of the first lepton from Z1", 200, -4., 4);
  new TH1F("Z1l2Eta", "Eta of the second lepton from Z1", 200, -4, 4);
  new TH1F("Z2l1Eta", "Eta of the first lepton from Z2", 200, -4., 4);
  new TH1F("Z2l2Eta", "Eta of the second lepton from Z2", 200, -4, 4);
 
  //new TH1F("Zl1Phi", "Phi of the first lepton from Z", 200, -M_PI, M_PI);
  //new TH1F("Zl2Phi", "Phi of the second lepton from Z", 200, -M_PI, M_PI);

  new TH1F("nEventsPerType", "Selected events per category", 6, -1.5, 4.5);
  new TH1F("nEventsPerType_standard4l", "Selected events per category after standard 4l cut", 6, -1.5, 4.5);
  new TH1F("nEventsPerType_4l+2j", "Selected events per category after starndard 4l+2j cut", 6, -1.5, 4.5);
  new TH1F("nEventsPerType_mass4l", "Selected events per category after mass4l cut", 6, -1.5, 4.5);
  new TH1F("nEventsPerType_tmva", "Selected events per category after tmva cut", 6, -1.5, 4.5);
  //new TH1F("nEventsPerType_bjet", "Selected events per category after bjet cut", 5, -1.5, 3.5);
  if (isMC()){
    new TH1F("nEventsPerTypeWt", "Selected events per category (Weighted)", 6, -1.5, 4.5);
    new TH1F("nEventsPerTypeWt_standard4l", "Selected events per category after starndard 4l cut(Weighted)", 6, -1.5, 4.5);
new TH1F("nEventsPerTypeWt_4l+2j", "Selected events per category after starndard 4l+2j cut(Weighted)", 6, -1.5, 4.5);
    new TH1F("nEventsPerTypeWt_mass4l", "Selected events per category after mass4l cut(Weighted)", 6, -1.5, 4.5);
    new TH1F("nEventsPerTypeWt_tmva", "Selected events per category after tmva cut(Weighted)", 6, -1.5, 4.5);
    //new TH1F("nEventsPerTypeWt_bjet", "Selected events per category after tmva cut(Weighted)", 5, -1.5, 3.5);
  }


  //MVA histograms
  new TH1F("mvaoutput_bdt",  "MVA output bdt value is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_bkg260",  "MVA output bdt value for 260 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_bkg300",  "MVA output bdt value for 300 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_bkg400",  "MVA output bdt value for 400 is : ", 20, -1., 1.); 
  new TH1F("mvaoutput_bdt_bkg500",  "MVA output bdt value for 500 is : ", 20, -1., 1.); 
  new TH1F("mvaoutput_bdt_bkg600",  "MVA output bdt value for 600 is : ", 20, -1., 1.); 
  new TH1F("mvaoutput_bdt_bkg700",  "MVA output bdt value for 700 is : ", 20, -1., 1.); 
  new TH1F("mvaoutput_bdt_bkg750",  "MVA output bdt value for 750 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_bkg800",  "MVA output bdt value for 800 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_bkg900",  "MVA output bdt value for 900 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_bkg1000",  "MVA output bdt value for 1000 is : ", 20, -1., 1.);



  new TH1F("mvaoutput_bdt_4mu",  "MVA output bdt value for 4mu channel is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4mu_bkg260",  "MVA output bdt value for 4mu for 260 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4mu_bkg300",  "MVA output bdt value for 4mu for 300 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4mu_bkg400",  "MVA output bdt value for 4mu for 400 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4mu_bkg500",  "MVA output bdt value for 4mu for 500 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4mu_bkg600",  "MVA output bdt value for 4mu for 600 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4mu_bkg700",  "MVA output bdt value for 4mu for 700 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4mu_bkg750",  "MVA output bdt value for 4mu for 750 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4mu_bkg800",  "MVA output bdt value for 4mu for 800 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4mu_bkg900",  "MVA output bdt value for 4mu for 900 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4mu_bkg1000",  "MVA output bdt value for 4mu for 1000 is : ", 20, -1., 1.);

  new TH1F("mvaoutput_bdt_4e",  "MVA output bdt value for 4e channel is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4e_bkg260",  "MVA output bdt value for 4e for 260 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4e_bkg300",  "MVA output bdt value for 4e for 300 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4e_bkg400",  "MVA output bdt value for 4e for 400 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4e_bkg500",  "MVA output bdt value for 4e for 500 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4e_bkg600",  "MVA output bdt value for 4e for 600 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4e_bkg700",  "MVA output bdt value for 4e for 700 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4e_bkg750",  "MVA output bdt value for 4e for 750 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4e_bkg800",  "MVA output bdt value for 4e for 800 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4e_bkg900",  "MVA output bdt value for 4e for 900 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_4e_bkg1000",  "MVA output bdt value for 4e for 1000 is : ", 20, -1., 1.);

  new TH1F("mvaoutput_bdt_2e2mu",  "MVA output bdt value for 2e2mu channel is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_2e2mu_bkg260",  "MVA output bdt value for 2e2mu for 260 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_2e2mu_bkg300",  "MVA output bdt value for 2e2mu for 300 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_2e2mu_bkg400",  "MVA output bdt value for 2e2mu for 400 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_2e2mu_bkg500",  "MVA output bdt value for 2e2mu for 500 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_2e2mu_bkg600",  "MVA output bdt value for 2e2mu for 600 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_2e2mu_bkg700",  "MVA output bdt value for 2e2mu for 700 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_2e2mu_bkg750",  "MVA output bdt value for 2e2mu for 750 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_2e2mu_bkg800",  "MVA output bdt value for 2e2mu for 800 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_2e2mu_bkg900",  "MVA output bdt value for 2e2mu for 900 is : ", 20, -1., 1.);
  new TH1F("mvaoutput_bdt_2e2mu_bkg1000",  "MVA output bdt value for 2e2mu for 1000 is : ", 20, -1., 1.);


  //  new TH1F("mvaoutput_bdtg",  "MVA output bdtg value is : ", 100, -1.0, 1.0);
  //  new TH1F("mvaoutput_mlp",  "MVA output mlp value is : ", 50, 0., 1.0);

  new TH1F("mass4l_4mu",  "4l invariant mass(4mu)", 200, 0., 200.);
  new TH1F("mass4l_4e",  "4l invariant mass(4e)", 200, 0., 200.);
  new TH1F("mass4l_2e2mu",  "4l invariant mass(2e2mu)", 200, 0., 200.);

  
  //new TH1F("jet1Mass",  "Jet1 mass (4l + #jets >=2) plus |m4l-mH| <= 10", 20----------------------------------------------------------------------------------------------------------------0, -1000., 1000.);
  new TH1F("jet1Pt",  "Jet1 pt (4l + #jets >=2) plus |m4l-mH| <= 10", 100, 0., 1000.);
  new TH1F("jet1Eta",  "Jet1 Eta (4l + #jets >=2) plus |m4l-mH| <= 10", 60, -3., 3.);
  new TH1F("jet1Phi",  "Jet1 phi (4l + #jets >=2) plus |m4l-mH| <= 10", 64, -3.2, 3.2);
  new TH1F("jet1Energy",  "Jet1 energy (4l + #jets >=2) plus |m4l-mH| <= 10", 100, 0., 1000.);
  //new TH2F("jet1EnergyvsPt",  "Jet1 energyvspt (4l + #jets >=2) plus |m4l-mH| <= 10", 200, -1000., 1000., 200, -1000., 1000.);

  //new TH1F("jet2Mass",  "Jet2 mass (4l + #jets >=2) plus |m4l-mH| <= 10", 200, -1000., 1000.);  
  new TH1F("jet2Pt",  "Jet2 pt (4l + #jets >=2) plus |m4l-mH| <= 10", 100, 0., 1000.);
  new TH1F("jet2Eta",  "Jet2 Eta (4l + #jets >=2) plus |m4l-mH| <= 10", 60, -3., 3.);
  new TH1F("jet2Phi",  "Jet2 phi (4l + #jets >=2) plus |m4l-mH| <= 10", 64, -3.2, 3.2);
  new TH1F("jet2Energy",  "Jet2 energy (4l + #jets >=2) plus |m4l-mH| <= 10", 100, 0., 1000.);
  //new TH2F("jet2EnergyvsPt",  "Jet2 energyvspt (4l + #jets >=2) plus |m4l-mH| <= 10", 200, -1000., 1000., 200, -1000., 1000.);

  // COntrol region
  new TH1F("mass2l_CR", "Invariant mass of di-lepton candidate in CR", 200, 0., 200.);
  new TH1F("dijetMass_CR",  "DiJet mass in CR", 40., 0., 200.);
  new TH1F("jet1pt_CR",  "Jet1 pt in CR", 80., 0., 400.);
  new TH1F("jet2pt_CR",  "Jet2 pt in CR", 80., 0., 400.);
  new TH1F("btag1_CR",  "b-tag score of first bJet in CR", 25, 0., 1.);
  new TH1F("btag2_CR",  "b-tag score of second bJet in CR", 25, 0., 1.);
  new TH1F("Met_CR",  "Met in CR", 80, 0., 400.);
  new TH1F("fourLPt_CR",  "four L pt in CR", 50, 0., 200.);

 
  new TH1F("dijetMass",  "DiJet mass (4l + #jets >=2) plus |m4l-mH| <= 10", 100, 0., 200.);  
  new TH1F("dijetPt",  "DiJet pt (4l + #jets >=2) plus |m4l-mH| <= 10", 100, 0., 1000.);
  new TH1F("dijetMass_csv",  "DiJet mass (4l + #jets >=2) plus |m4l-mH| <= 10 for highest b-tagged jets", 100, 0., 200.);  
  new TH1F("dijetPt_csv",  "DiJet pt (4l + #jets >=2) plus |m4l-mH| <= 10 for highest b-tagged jets", 100, 0., 1000.);
  new TH1F("dijetEta",  "DiJet Eta (4l + #jets >=2) plus |m4l-mH| <= 10", 60, -3., 3.);
  //new TH1F("dijetPhi",  "DiJet phi (4l + #jets >=2) plus |m4l-mH| <= 10", 64, -3.2, 3.2);
  new TH1F("dijetEnergy",  "DiJet energy (4l + #jets >=2) plus |m4l-mH| <= 10", 200, -1000., 1000.);  
  //new TH2F("dijetEnergyvsPt",  "dijet energyvspt (4l + #jets >=2) plus |m4l-mH| <= 10", 200, -1000., 1000., 200, -1000., 1000.);

  new TH1F("dR_JJ",  "dR_JJ", 60, 0., 6.);
  new TH1F("dPhi_JJ",  "dPhi_JJ", 60, -3.14, 3.14);
  new TH1F("dEta_JJ",  "dEta_JJ", 100, 0., 10.);
  new TH1F("Product_Eta_JJ",  "Product_Eta_JJ", 800, -16., 16.);
  new TH1F("dR_JJ_csv",  "dR_JJ for highest b-tagged jets", 60, 0., 6.);
  new TH1F("dPhi_JJ_csv",  "dPhi_JJ for highest b-tagged jets", 60, -3.14, 3.14);
  new TH1F("dEta_JJ_csv",  "dEta_JJ for highest b-tagged jets", 100, 0., 10.);
  new TH1F("Product_Eta_JJ_csv",  "Product_Eta_JJ for highest b-tagged jets", 800, -16., 16.);
  new TH1F("dR_ZZ",  "dR_ZZ", 60, 0., 6.);
  new TH1F("dPhi_ZZ",  "dPhi_ZZ", 60, -3.14, 3.14);
  
  new TH1F("Met",  "Met", 150, 0., 150.);
  new TH1F("MetPhi",  "MetPhi", 60, -3.14, 3.14);

  new TH1F("b1hbb_hel",  "b1hbb_helicity", 20, -1., 1.);
  new TH1F("b2hbb_hel",  "b2hbb_helicity", 20, -1., 1.);
  new TH1F("z1hzz_hel",  "z1hzz_helicity", 20, -1., 1.);
  new TH1F("z2hzz_hel",  "z2hzz_helicity", 20, -1., 1.);
  new TH1F("h_hbb_hel",  "h_hbb_helicity", 20, -1., 1.);
  new TH1F("h_hzz_hel",  "h_hzz_helicity", 20, -1., 1.);

  new TH1F("jTotalscalarpt",  "Scalar Sum of Pt of JJ", 150, 0., 1500.);
  new TH1F("jTotalscalarpt_csv",  "Scalar Sum of Pt of JJ for highest b-tagged jets", 150, 0., 1500.);
  new TH1F("m4l_PtAll",  "scalar sum of 4l pt", 50, 0., 500.);
  new TH1F("hh_PtDiff",  "HH_Ptdifference", 80, 0., 200.);
  new TH1F("hh_MDiff",  "HH_Mdifference", 40, 0., 100.);
  new TH1F("scalarPtAll",  "scalarPtAll", 150, 0., 1500.);
  new TH1F("ZZ_massdiff",  "ZZ_massdifference", 40, 0., 100.);
  new TH1F("dR_h4l_hjj",  "dR_h4l_hjj", 50, 0., 10.);
  new TH1F("dPhi_h4l_hjj",  "dPhi_h4l_hjj", 60, -3.14, 3.14);
  new TH1F("hh_PtDiff_csv",  "HH_Ptdifference for highest b-tagged jets", 80, 0., 200.);
  new TH1F("hh_MDiff_csv",  "HH_Mdifference for highest b-tagged jets", 40, 0., 100.);
  new TH1F("scalarPtAll_csv",  "scalarPtAll for highest b-tagged jets", 150, 0., 1500.);
  new TH1F("dR_h4l_hjj_csv",  "dR_h4l_hjj for highest b-tagged jets", 50, 0., 10.);
  new TH1F("dPhi_h4l_hjj_csv",  "dPhi_h4l_hjj for highest b-tagged jets", 60, -3.14, 3.14);
  //  new TH1I("nLooseLeptons",  "#LooseLeptons", 40, 0, 4);
  //new TH1I("nbJets",  "#bJets", 40, 0, 4);
  new TH1F("btag1",  "b-tag score of first bJet", 50, 0., 1.);
  new TH1F("btag2",  "b-tag score of second bJet", 50, 0., 1.);
  new TH1F("btag1_csv",  "b-tag score of first bJet for highest b-tagged jets", 50, 0., 1.);
  new TH1F("btag2_csv",  "b-tag score of second bJet for highest b-tagged jets", 50, 0., 1.);

  new TH1F("Z1massC0", "ll invariant mass of first Z candidate (mmmm)", 200, 0., 200.);
  new TH1F("Z1massC1", "ll invariant mass of first Z candidate (mmee)", 200, 0., 200.);
  new TH1F("Z1massC2", "ll invariant mass of first Z candidate (eemm)", 200, 0., 200.);
  new TH1F("Z1massC3", "ll invariant mass of first Z candidate (eeee)", 200, 0., 200.);

  new TH1F("Z2massC0", "ll invariant mass of second Z candidate (mmmm)", 200, 0., 200.);
  new TH1F("Z2massC1", "ll invariant mass of second Z candidate (mmee)", 200, 0., 200.);
  new TH1F("Z2massC2", "ll invariant mass of second Z candidate (eemm)", 200, 0., 200.);
  new TH1F("Z2massC3", "ll invariant mass of second Z candidate (eeee)", 200, 0., 200.);

  new TH1F("dRZ1Z2C0", "dR(Z1, Z2) (mmmm)", 100, 0, 2 * M_PI);
  new TH1F("dRZ1Z2C1", "dR(Z1, Z2) (mmee)", 100, 0, 2 * M_PI);
  new TH1F("dRZ1Z2C2", "dR(Z1, Z2) (eemm)", 100, 0, 2 * M_PI);
  new TH1F("dRZ1Z2C3", "dR(Z1, Z2) (eeee)", 100, 0, 2 * M_PI);

  histf()->cd();
  histf()->ls();
}
// -------------------------------
// Clear vectors before event loop
// -------------------------------
void DiHiggs4l::clearLists() {
  PhysicsObjSelector::clear();
  vtxList_.clear();
}
// -------------------
// The main event loop
// -------------------
void DiHiggs4l::eventLoop()
{
  int nPrint = std::max(10000, nEvents()/1000);

  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;
  /*  int n_mm_ini = 0;
  int n_ee_ini = 0;
  int n_me_ini = 0;*/
  /*int n_mu_chnl = 0;
  int n_det_mu_chnl = 0;
  int n_el_chnl = 0;
  int n_det_el_chnl = 0;
  int n_muel_chnl = 0;
  int n_det_muel_chnl = 0;*/
  
  //  int n_mm_fin = 0;
  // int n_ee_fin = 0;
  // int n_me_fin = 0;

  /*  float b1hbb_hel = 99.;
  float b2hbb_hel = 99.;
  float z1hzz_hel = 99.;
  float z2hzz_hel = 99.;
  float h_hbb_hel = 99.;
  float h_hzz_hel = 99.;*/

  // Muon scale factor for 2017
  //  TFile *mu_scale_factors = new TFile("ScaleFactors_mu_Moriond2018_final.root"); // Nicola
  TFile *mu_scale_factors = new TFile("/lustrehome/ppalit/CMSSW_9_4_9/src/offlineAnalysis_v2/final_HZZ_muon_SF_2017_newLooseIso_mupogSysts_paper.root");// Twiki
  TH2F *mu_scale_2017 = (TH2F*)gDirectory->Get("FINAL");

  //  TFile *ele_scale_factors2017 = new TFile("egammaEffi_txt_EGM2D_Moriond2018v1.root"); // Nicola
  TFile *ele_scale_factors2017 = new TFile("/lustrehome/ppalit/CMSSW_9_4_9/src/offlineAnalysis_v2/ElectronSF_Legacy_2017_NoGap.root"); //Twiki 
  TH2F *ele_scale_2017 = (TH2F*)gDirectory->Get("EGamma_SF2D");
  //  TFile *ele_scale_factors2017_gap = new TFile("egammaEffi_txt_EGM2D_Moriond2018v1_gap.root"); // Nicola
  TFile *ele_scale_factors2017_gap = new TFile("/lustrehome/ppalit/CMSSW_9_4_9/src/offlineAnalysis_v2/ElectronSF_Legacy_2017_Gap.root"); // Twiki
  TH2F *ele_scale_2017_gap = (TH2F*)gDirectory->Get("EGamma_SF2D");

  //  TFile *ele_RecoEff_2017_lowEt = new TFile("egammaEffi_txt_EGM2D_runBCDEF_passingRECO_lowEt.root"); //Low Et<20 // Nicola               
  TFile *ele_RecoEff_2017_lowEt = new TFile("/lustrehome/ppalit/CMSSW_9_4_9/src/offlineAnalysis_v2/Ele_Reco_LowEt_2017.root"); // Twiki                                                         
  TH2F *ele_Reco_eff_2017_lowEt = (TH2F*)gDirectory->Get("EGamma_SF2D");

  //  TFile *ele_RecoEff_2017_highEt = new TFile("egammaEffi_txt_EGM2D_runBCDEF_passingRECO.root"); //high Et>20 // Nicola           
  TFile *ele_RecoEff_2017_highEt = new TFile("/lustrehome/ppalit/CMSSW_9_4_9/src/offlineAnalysis_v2/Ele_Reco_2017.root"); // Twiki
  TH2F *ele_Reco_eff_2017_highEt = (TH2F*)gDirectory->Get("EGamma_SF2D");

  //***************************************************************
  TString strSystTitle[9] ={
    "Nominal",
    "PDFScaleDn",
    "PDFScaleUp",
    "QCDScaleDn",
    "QCDScaleUp",
    "AsDn",
    "AsUp",
    "PDFReplicaDn",
    "PDFReplicaUp"
  };
  
  // Open the file
  TFile* fin = TFile::Open("/lustrehome/ppalit/CMSSW_9_4_9/src/offlineAnalysis_v2/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  TSpline3* ggZZ_kf[9]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
  for(int f=0;f<9;f++){
    ggZZ_kf[f] = (TSpline3*)fin->Get(Form("sp_kfactorratio_%s", strSystTitle[f].Data()));
  }   
  //  fin->Close();
  //*****************************************************************************

  // --------------------
  // Start the event loop
  // -------------------- 
  string lastFile;
  int fevt = firstEvent() > -1 ? firstEvent() : 0;
  int levt = lastEvent() > -1 ? lastEvent() : nEvents();
  cout << ">>> Event range: [" << fevt << ", " << levt -1 << "]" << endl;
  int nEventSel = 0;
  evtWeightSum_ = 0;
  for (int ev = fevt; ev < levt; ++ev) {
    clearEvent(); // reset tree variables
    clearLists(); // reset analysis related lists  for each event
    int lflag = chain()->LoadTree(ev);
    int nbytes = getEntry(lflag);    // returns total bytes read
    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName()));

    histf()->cd();

    // For data or for MC without pileup
    double puevWt {1.0};
    //******************************
    
    Double_t eff_weight = 1.;

    // Show status of the run
    const Event& evt = eventColl()->at(0);
    unsigned long run   = evt.run;
    unsigned long event = evt.event;
    unsigned long lumis = evt.lumis;
    int nv = evt.npuvtx;
    //float evweight = evt.evwt;

    if (isMC() && usePUWt()) {
      //int npu {0};
      //puevWt = wtPileUp(npu);
      puevWt = dopuweight(nv,0);
    }

    //***********************************************************8
    // ggZZ kfactor

    //std::cout << " Debug_v1  : " << std::endl;
    // qqZZ kfactor
    /*
    int finalState=-999;
    float weight_factor = 1.;*/
    //if (isMC() && isQQtoZZ()){
    //  for 
    //}
    

    //*********************************************************************


    // Show status of the run

    /*
    if (currentFile != lastFile)
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << chain()->GetCurrentFile()->GetName()
	   << " <<< Run# " << setw(8) << run
	   << " Lumis# " << setw(6) << lumis
	   << " Event# " << setw(12) << event << " >>> "
	   << " Events proc. " << setw(9) << ev << "(of " << setw(9) << levt-1 << ")"
	   << endl;
    lastFile = currentFile;
    

    // Show the status
    if (ev%nPrint == 0 || firstEvent() > -1)
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << currentFile
	   << " <<< Run# " << setw(8) << run
	   << " Lumis# " << setw(6) << lumis
	   << " Event# " << setw(12) << event << " >>> "
	   << " Events proc. " << setw(8) << ((firstEvent() > -1) ? ev - firstEvent() : ev)
	   << endl;
    */
 
    //   std::cout << "use event list : " << useEventList_ << std::endl;
    // std::cout << "eventIdMap().size() : " << eventIdMap().size() << std::endl;
      
    // Select a set of events by [run, event]
    if (useEventList_ && eventIdMap().size() > 0) {
      std::ostringstream mkey;
      mkey << run << "-" << lumis << "-" << event;
      //if (eventIdMap().find(mkey.str()) != eventIdMap().end()) continue; 
      if (eventIdMap().find(mkey.str()) == eventIdMap().end()) continue;
    }

    histf()->cd();
    histf()->cd("DiHiggs4l");

    if (isMC()) {
      AnaUtil::fillHist1D("puweight", puevWt);
      //**********************
      
      // only when exclusive jets datasets are used (for np = 0)
      //      if (0) std::cout << "== nMEPartons: " << genEventColl()->at(0).nMEPartons
      //		       << ", lheNOutPartons: " << genEventColl()->at(0).lheNOutPartons
      //	       << std::endl;
      if (selectPM_ && genEventColl()->at(0).nMEPartons != nMEPartons_) continue;

      double wt = genEventColl()->size() ? genEventColl()->at(0).evtWeight : 1;
      //      double wt = evt.evwt>0 ? evt.evwt : 1.; 
      //      std::cout << "wt for the evnt is:    " << wt << std::endl;
      AnaUtil::fillHist1D("evtweight", wt);
      evtWeightSum_ += wt;   // this is used for final normalization
    }
    AnaUtil::fillHist1D("evtCutFlow", 0);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 0, puevWt);

    // Check this later
    if (useGenFilter_ && genAna_ != nullptr) {
      genAna_->setEvent(genParticleColl());
      if (!genAna_->filter()) continue;
    }

    AnaUtil::fillHist1D("evtCutFlow", 1);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 1, puevWt);

    int ngen_mm = 0,
      ngen_tt = 0,
      ngen_ee = 0;
    int npat_mm = 0,
      npat_tt = 0,
      npat_mm_det = 0;
    int ngen_det_mm = 0,
      ngen_det_ee = 0; //ngen_det_tt = 0;
    bool gen_mumu_pass = false;
    bool gen_elel_pass = false;                                     
    bool gen_muel_pass = false;
    bool pat_mumu_pass = false;
    bool pat_elel_pass = false;                                     
    bool pat_muel_pass = false;
    bool gen_det_mumu_pass = false;
    bool gen_det_elel_pass = false;
    bool gen_det_muel_pass = false;
    bool pat_mumu_pass_det = false;
    bool pat_elel_pass_det = false;                                  
    bool pat_muel_pass_det = false;    

    std::vector<int> vec_diZmassgtZero;
    std::vector<float> vec_diZmass;
    bool diZmassgtZero = true;
    //    std::vector<vhtm::GenParticle> genMu_vec;
    //std::vector<vhtm::GenParticle> genbJet_vec;
    float genzzmass = 1.;
    float  weight_kfactor=1.;
    float weight_qqzzkfactor=1.;
    double weight_scalefactor = 1.;
    double weight_elerecofactor = 1.;
    for (const auto& gp: *genParticleColl()) {
      /*      int pdg = std::abs(gp.pdgId);
      int st  = gp.status;
      float pt  = gp.pt;
      float eta = gp.eta;
      float vx  = gp.vx;
      float vy  = gp.vy;
      float vz  = gp.vz;
      if ((st == 1 && (pdg == 11 || pdg == 13 )) || (st==2 && pdg == 15)) {
        int mid = -1;
        int indx = getMotherId(gp, mid);
        if (indx < 0) continue;
        const auto& mgp = genParticleColl()->at(indx);
        if (std::abs(mgp.pdgId) == 23) { // Z
          int mmid = -1;
          int indx = getMotherId(mgp, mmid);
          if (indx < 0) continue;
          const auto& mmgp = genParticleColl()->at(indx);
          if (std::abs(mmgp.pdgId) == 25) { // H
            if (pdg == 13) ngen_mm++;  
            else if (pdg == 11) ngen_ee++;
            else if (pdg == 15) ngen_tt++;
            if (pt > 5 && std::abs(eta) < 2.4 && pdg == 13) ngen_det_mm++;
            else if (pt > 7 && std::abs(eta) < 2.5 && pdg == 11) ngen_det_ee++;
          }    
        }
	}*/

      //      std::cout << " gp info ====> pdgid = " << gp.lepPdgId << " : pt = " << gp.lepPt << " : eta = " << gp.lepEta << std::endl;   

      /*      
      if (gp.ZZpdgid == 13 || gp.ZZpdgid == -13) ngen_mm++;
      if ((gp.ZZpdgid == 13 || gp.ZZpdgid == -13) && std::fabs(gp.ZZeta) < 2.5) ngen_det_mm++;
      //      if ((gp.ZZpdgid == 13 || gp.ZZpdgid == -13) && gp.ZZpt > 5 && std::fabs(gp.ZZeta) < 2.4) ngen_det_mm++;
      if (gp.ZZpdgid == 11 || gp.ZZpdgid == -11) ngen_ee++;
      if ((gp.ZZpdgid == 11 || gp.ZZpdgid == -11) && std::fabs(gp.ZZeta) < 2.5) ngen_det_ee++;*/
      //      if ((gp.ZZpdgid == 11 || gp.ZZpdgid == -11) && gp.ZZpt > 7 && std::fabs(gp.ZZeta) < 2.5) ngen_det_ee++;

      if (gp.lepPdgId == 13 || gp.lepPdgId == -13) ngen_mm++;
      if ((gp.lepPdgId == 13 || gp.lepPdgId == -13) && std::fabs(gp.lepEta) < 2.4) {
	ngen_det_mm++;
	//	genMu_vec.push_back(gp);
      }
      //      if ((gp.ZZpdgid == 13 || gp.ZZpdgid == -13) && gp.ZZpt > 5 && std::fabs(gp.ZZeta) < 2.4) ngen_det_mm++;
      if (gp.lepPdgId == 11 || gp.lepPdgId == -11) ngen_ee++;
      if ((gp.lepPdgId == 11 || gp.lepPdgId == -11) && std::fabs(gp.lepEta) < 2.4) ngen_det_ee++;
      //std::cout << "debug version 3" << std::endl;  
    //      if ((gp.ZZpdgid == 11 || gp.ZZpdgid == -11) && gp.ZZpt > 7 && std::fabs(gp.ZZeta) < 2.5) ngen_det_ee++;
      //std::cout << "gp_pdgid : " << gp.lepPdgId << std::endl;
      /*
      if (isMC() && isQQtoZZ()){
	if (gp.ZZmass > 0. && (gp.ZZpdgid == 13 || gp.ZZpdgid == -13) && (gp.ZZpdgid == 11 || gp.ZZpdgid == -11) && (gp.ZZpdgid == 15 || gp.ZZpdgid == -15)) finalState=1;
	else finalState=2;
	weight_factor = kfactor_qqZZ_qcd_M(gp.ZZmass,finalState);
	}*/
      /*std::cout << "The kfactor for qqZZ is= " << weight_factor<< std::endl;
      puevWt = puevWt*weight_factor;
      */
      
      double ggzz_kf_wgt[9];

      if (isMC() && isGluGluHToZZ()){
	//if (isMC()){ 
	for(int f=0;f<9;f++) {ggzz_kf_wgt[f] = ggZZ_kf[f]->Eval(gp.diZmass); // Evaluate at the true m4l
	  //	   std::cout << "diZMass : " << gp.diZmass << "f :  " << f << "\t" << "ggzz_kf_wgt = " << ggzz_kf_wgt[f] << std::endl;
	}
	

	weight_kfactor=ggzz_kf_wgt[0]; // Using the nominal one
	//std::cout << "The kfactor for ggZZ is= " << weight_kfactor<< std::endl;
		//
      }
      
      genzzmass = gp.diZmass;
      //      std::cout << "genzzmass within gen loop = " << genzzmass << std::endl;
    }

      //    std::cout << "genzzmass = " << genzzmass << std::endl;
    
      // std::cout << "pu weight = " << puevWt << "\t" << "k factor = " << weight_kfactor << std::endl;
    puevWt = puevWt*weight_kfactor;
    //std::cout << "pu weight*weight_kfactor = " << puevWt << std::endl;
      /*
      if (gp.diZmass > 0.) {
	vec_diZmassgtZero.push_back(1);
	vec_diZmass.push_back(gp.diZmass);
      }
      else vec_diZmassgtZero.push_back(0);

    }
      */
    //   std::cout << "eventId " << run << " " << lumis << " " << event << std::endl;
    //    std::cout << "n_genmuons : " << ngen_mm << "n_geneles : " << ngen_ee << std::endl;
    if (ngen_mm >= 4){
      gen_mumu_pass = true;
      //n_mu_chnl++;
    }
    else if (ngen_ee >= 4){
      gen_elel_pass = true;
      //n_el_chnl++;
    }
    else if (ngen_mm >= 2 && ngen_ee >= 2){
      gen_muel_pass = true;
      //n_muel_chnl++;
    }

    double qqzz_kf_wgt;
    int finalState=-999;
    if (isMC() && isqqZZ()) {
      if (gen_mumu_pass == true || gen_elel_pass == true) finalState=1;
      else if (gen_muel_pass ==  true ) finalState=2;
      weight_qqzzkfactor = kfactor_qqZZ_qcd_M(genzzmass, finalState);
      //     std::cout << "The kfactor for qqZZ is= " << weight_qqzzkfactor<< std::endl;
      //std::cout << "finalState = " << finalState<< std::endl;
    }
 
    //std::cout << "pu weight = " << puevWt << "\t" << "k factor for qqzz= " << weight_qqzzkfactor << std::endl;
    puevWt = puevWt*weight_qqzzkfactor;
    // std::cout << "pu weight*weight_qqzzkfactor = " << puevWt << std::endl;

    /*
    for (int i = 0; i < vec_diZmassgtZero.size(); i++){
      if (vec_diZmassgtZero.at(i) == 0) diZmassgtZero = false; 
    }

    std::cout << "crash after this : " << std::endl;
    int finalState=-999;
    if (isMC() && isqqZZ()){
      if (gp.diZmass > 0.) {
        vec_diZmass.push_back(gp.diZmass);
      }
      if ((gen_mumu_pass == true || gen_elel_pass == true) && diZmassgtZero == true ) finalState=1;
      else if (gen_muel_pass ==  true && diZmassgtZero == true ) finalState=2;
      weight_kfactor = kfactor_qqZZ_qcd_M(vec_diZmass[0],finalState);
      std::cout << "crash after this v2: " << std::endl;
      std::cout << " final state  : " << finalState << std::endl;
      std::cout << " vec_diZmass  : " << vec_diZmass[0] << std::endl;
      std::cout << " weight_kfactor  : " << weight_kfactor << std::endl;
      std::cout << " puevwt before weight_kfactor  : " << puevWt << std::endl;
      puevWt = puevWt*weight_kfactor;
      std::cout << " puevwt after weight_kfactor  : " << puevWt << std::endl;
    }
    */
    //    if (!gen_mumu_pass) continue;
    //if (!gen_elel_pass) continue;
    //if (!gen_muel_pass) continue;

    AnaUtil::fillHist1D("evtCutFlow", 2);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 2, puevWt);

    if (ngen_det_mm >= 4){
      gen_det_mumu_pass = true;
      //n_det_mu_chnl++;
    }
    else if (ngen_det_ee >= 4){
      gen_det_elel_pass = true;
      //n_det_el_chnl++;
    }
    else if (ngen_det_mm >= 2 && ngen_det_ee >= 2){
      gen_det_muel_pass = true;
      //n_det_muel_chnl++;
    }
 
    //if (!gen_det_mumu_pass) continue;
    //    if (!gen_det_elel_pass) continue;
    //if (!gen_det_muel_pass) continue;

    
 
    histf()->cd();
    histf()->cd("DiHiggs4l");

    /*
    std::vector<vhtm::Muon> muPat;
    std::vector<vhtm::Muon> muPat_det;
    std::vector<vhtm::Electron> elePat;
    std::vector<vhtm::Electron> elePat_det;    

    int npat_bef_det= 0;
    int npat_aft_det= 0;
    for (const auto& mu: *muonColl()){
      muPat.push_back(mu);
      npat_bef_det++;
      if (std::abs(mu.eta) < 2.4) {muPat_det.push_back(mu);npat_aft_det++;}
    }

    std::cout << "npat_bef_det = " << npat_bef_det << std::endl;
    std::cout << "npat_aft_det = " << npat_aft_det << std::endl;
    for (const auto& ele: *electronColl()){
      elePat.push_back(ele);
      if (std::abs(ele.eta) < 2.5) elePat_det.push_back(ele);
    }
 
    if (muPat.size() >= 4){
      pat_mumu_pass = true;
      //n_mu_chnl++;
    }
    else if (elePat.size() >= 4){
      pat_elel_pass = true;
      //n_el_chnl++;
    }
    else if (muPat.size() >= 2 && elePat.size() >= 2){
      pat_muel_pass = true;
      //n_muel_chnl++;
    }

    if (!pat_mumu_pass) continue;
    //if (!pat_elel_pass) continue;
    //if (!pat_muel_pass) continue;

    */
    AnaUtil::fillHist1D("evtCutFlow", 3);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 3, puevWt);

    //AnaUtil::fillHist1D("evtCutFlow", 3);
    //if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 3, puevWt);

    // at least 1 good PV
    // good vertex finding
    op.verbose = (logOption() >> 1 & 0x1);
    findVtxInfo(vtxList_, op, fLog());
    size_t ngoodVtx = vtxList_.size();
    // std::cout << "ngoodVtx ===" << ngoodVtx << std::endl;
    //    if (ngoodVtx < 1) continue;
    //    if (nv < 1) continue;

    /*
    if (muPat_det.size() >= 4){
      pat_mumu_pass_det = true;
      //n_mu_chnl++;
    }
    else if (elePat_det.size() >= 4){
      pat_elel_pass_det = true;
      //n_el_chnl++;
    }
    else if (muPat_det.size() >= 2 && elePat_det.size() >= 2){
      pat_muel_pass_det = true;
      //n_muel_chnl++;
    }

    if (!pat_mumu_pass_det) continue;
    //if (!pat_elel_pass) continue;
    //if (!pat_muel_pass) continue;

    */
    AnaUtil::fillHist1D("evtCutFlow", 4);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 4, puevWt);

#ifdef SKIP_DUPLICATE_ALL
    // Duplicate event removal
    if (!isMC()) {
      std::ostringstream mkey;
      mkey << run << "-" << lumis << "-" << event;
      std::string evs {mkey.str()};
      if (eventIdStore_.find(evs) != eventIdStore_.end()) {
	if (1) cout << "DuplicateAll: " << evs << endl;
	continue;
      }
      else {
	eventIdStore_.insert({evs, 1});
      }
    }
#endif
    //    AnaUtil::fillHist1D("isTriggered", (isTriggered(true, false)?1:0), puevWt);

    // is event triggered?
    //    if (0) dumpTriggerPaths(std::cout, true);
    // if (useTrigger() && !isTriggered(true, false)) continue;

    AnaUtil::fillHist1D("evtCutFlow", 5);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 5, puevWt);

    // main analysis object selection
    double vz = (vtxList_.size() > 0) ? vtxList_.at(0).z : -999;
    //    dumpEvent(vz, false, true);

    //std::cout << "Vertex Z is ++++++++" << vz << std::endl;
    findObjects(vz, 1); // do not use event weight

    //    if (genAna_ != nullptr && dumpGenInfo_) genAna_->dumpEvent();
 
    // access selected objects
    // default H4l analysis
    const auto& elePhotonPairList = getTightIsoElePhotonPairList();
    int nEle = elePhotonPairList.size();

    //    std::cout << "elec size is =====> " << nEle << std::endl;

    const auto& muPhotonPairList  = getTightIsoMuPhotonPairList();
    int nMuon = muPhotonPairList.size();//number of isolated muons
 
    //   std::cout << "Muon size is ======> " << nMuon << std::endl;

    //    const auto& tauList = getIsoTauList();
    // size_t nTau = tauList.size();

    bool mm_ch = false,
      ee_ch = false,
      me_ch = false;

    if (nMuon >= 4) {
      mm_ch = true;
      //nmm_ini++;                                                   
    }
    else if (nEle >= 4) {
      ee_ch = true;
      //nee_ini++;                                                   
    }
    else if (nMuon >= 2 && nEle >= 2 ) {
      me_ch = true;
      //nme_ini++;
    }

    histf()->cd();
    histf()->cd("DiHiggs4l");
    //    std::cout << " Debug_v7_1  : " << std::endl;
    
    AnaUtil::fillHist1D("nGoodMuon", nMuon, puevWt);
    AnaUtil::fillHist1D("nGoodEle", nEle, puevWt);
    // std::cout << " Debug_v7_2  : " << std::endl;
    //    if (!mm_ch || ee_ch || me_ch) continue;
    //    if (!ee_ch) continue;
    //    if (!me_ch) continue;

    //    std::cout << " Debug_v8  : " << std::endl;    
    AnaUtil::fillHist1D("evtCutFlow", 6);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 6, puevWt);
    
    // require at least 2 leptons of the same flavor (electron/muon)
    if ((nMuon + nEle) < 2) continue;
    //    if (nMuon < 2 && nEle < 2) continue;
    
    AnaUtil::fillHist1D("evtCutFlow", 7);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 7, puevWt);

    AnaUtil::fillHist1D("nTightIsoLeptons", nMuon + nEle , puevWt);

    // Highest pT lepton pT > 20
    double hLepPt = (nMuon > 0) ? muPhotonPairList[0].first.pt : 0.0;
    if (nEle > 0 && elePhotonPairList[0].first.pt > hLepPt) 
      hLepPt = elePhotonPairList[0].first.pt;
    AnaUtil::fillHist1D("hLepPt_stage0", hLepPt, puevWt);

    if (hLepPt < AnaUtil::cutValue(evselCutMap(), "hLepPtMin")) continue;
    AnaUtil::fillHist1D("evtCutFlow", 8);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 8, puevWt);

    // find Z candidates
    std::vector<ZCandidate> ZCandTmp;
    //   std::cout << "nmuon size : " << nMuon << "\t" << "nEle size : " << nEle << std::endl;

    if (nMuon >= 2) ZSelector<vhtm::Muon>(muPhotonPairList, ZCandTmp);
    if (nEle >= 2)  ZSelector<vhtm::Electron>(elePhotonPairList, ZCandTmp);
    AnaUtil::fillHist1D("nZcand", ZCandTmp.size(), puevWt);


    //    std::cout << "ZCandTmp size : " << ZCandTmp.size() << std::endl;
    if (ZCandTmp.empty()) continue;
    //if (ZCandTmp.size() < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 9);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 9, puevWt);
    //std::cout << "ZCandTmp size_v2 : " << ZCandTmp.size() << std::endl;

    if (ZCandTmp.size() > 1)
      std::sort(ZCandTmp.begin(), ZCandTmp.end(), dmComparator);

    // For all possible ZZ pairs, require both Z candidate masses(computed including the FSR photons if present) 
    // to be 12 < m(ll(g)) < 120 GeV 
    std::vector<ZCandidate> ZCandList;
    for (const auto& ZCand: ZCandTmp) {
      //std:;cout << "Z mass : " << ZCand.mass << std::endl;
      if (!((ZCand.mass > AnaUtil::cutValue(evselCutMap(), "ZMassLow") &&
	     ZCand.mass < AnaUtil::cutValue(evselCutMap(), "ZMassHigh")))) continue;
      AnaUtil::fillHist1D("massZcand", ZCand.mass, puevWt);
      ZCandList.push_back(ZCand);
    }
    AnaUtil::fillHist1D("nZcand", ZCandList.size(), puevWt);    
    if (verbosity_) std::cout << "#Z Candidates passing mass cut: " << ZCandList.size() << std::endl;

    //std::cout << "#Z Candidates passing mass cut: " << ZCandList.size() << std::endl;
    //    cout << "4. fsr photons: " << getFSRPhotonList().size() << endl;
 
    //dumpEvent(vz, false, true);
    //for (size_t nZ = 0; nZ < ZCandList.size(); ++nZ) 
      //      HZZ4lUtil::printZCandidate(ZCandList[nZ], "Z Index: " + std::to_string(nZ));
    
    if (ZCandList.empty()) continue;
    //if (ZCandList.size()<2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 10);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 10, puevWt);   

    // Now, add lepton isolation to the Z candidates found above
    addLeptonIsolation(elePhotonPairList, muPhotonPairList, ZCandList);
    /*
    if (verbosity_ > 1) {
     for (size_t nZ = 0; nZ < ZCandList.size(); ++nZ) 
       HZZ4lUtil::printZCandidate(ZCandList[nZ], "Z Index: " + std::to_string(nZ));
       }*/

#ifdef SKIP_DUPLICATE_ZMASS
    // Duplicate event removal
    if (!isMC()) {
      if (skipDuplicate_) { 
	std::ostringstream mkey;
	mkey << run << "-" << lumis << "-" << event;
	if (eventIdStore_.find(mkey.str()) != eventIdStore_.end()) {
	  if (0) cout << "Duplicate Event: " << mkey.str() << endl;
	  continue;
	}
      }
      evLog() << run << " " << lumis << " " << event << std::endl;
    }
#endif

    // Form Z Pairs
    zzvector ZZCandList;
    //std::cout << run << " " << lumis << " " << event << std::endl; 
    ZZCandidateSelector(ZCandList, ZZCandList);
    
    if (ZZCandList.empty()) continue;
    AnaUtil::fillHist1D("evtCutFlow", 11);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 11, puevWt);
    
    if (verbosity_ > 0) {
      std::cout << "=> Z candidate pair info: " << endl;  
      for (size_t i = 0; i < ZZCandList.size(); ++i) {
	double sumpt = ZZCandList[i].first.l1P4.Pt()  
	  + ZZCandList[i].first.l2P4.Pt() 
	  + ZZCandList[i].second.l1P4.Pt() 
	  + ZZCandList[i].second.l2P4.Pt();  
	std::cout << "Z Pair index: " << i << std::endl; 
	HZZ4lUtil::printZCandidate(ZZCandList[i].first, "Z1");
	HZZ4lUtil::printZCandidate(ZZCandList[i].second, "Z2");
	//std::cout << "=> SumPt of Z pair:" << sumpt << std::endl;
      }    
    }

    // Selection of best ZZ candidate
    // if more than two Z pairs are presented in the event, the ZZ pair with higher pt are considered 
    // to the computation of the Higgs invariant mass
    std::pair<ZCandidate, ZCandidate> BestZZCand = ZZCandList[0];
    int bestZZidx = 0;
    double BestZZCand_sumpt = 
        BestZZCand.first.l1RP4.Pt()  
      + BestZZCand.first.l2RP4.Pt() 
      + BestZZCand.second.l1RP4.Pt() 
      + BestZZCand.second.l2RP4.Pt();
    double BestZZCand_Z1mass_diff = fabs ((BestZZCand.first.l1RP4 + BestZZCand.first.l2RP4).M() - ZSpace::MZnominal);
    double BestZZCand_Z2lep_scalarsum = BestZZCand.second.l1RP4.Pt() + BestZZCand.second.l2RP4.Pt();
    //std::cout << "sumpt of ZZcand0 : " << BestZZCand_sumpt << "  mass of 4l of ZZcand0 : " << (BestZZCand.first.l1RP4+BestZZCand.first.l2RP4+BestZZCand.second.l1RP4+BestZZCand.second.l2RP4).M() << "  Z1 mass : " << (BestZZCand.first.l1RP4+BestZZCand.first.l2RP4).M() << "  Z2 mass : " << (BestZZCand.second.l1RP4+BestZZCand.second.l2RP4).M() << std::endl;
    if (ZZCandList.size() > 1) {
      //std::cout << "ZZcandlist.size final : " << ZZCandList.size() << std::endl;
      for (size_t i = 1; i < ZZCandList.size(); ++i) {
        double sumpt = ZZCandList[i].first.l1RP4.Pt()  
	             + ZZCandList[i].first.l2RP4.Pt() 
                     + ZZCandList[i].second.l1RP4.Pt() 
	             + ZZCandList[i].second.l2RP4.Pt();
	double Z1mass_diff = fabs ((ZZCandList[i].first.l1RP4 + ZZCandList[i].first.l2RP4).M() - ZSpace::MZnominal);
	double Z2lep_scalarsum = ZZCandList[i].second.l1RP4.Pt() + ZZCandList[i].second.l2RP4.Pt();
	//std::cout << "sumpt : " << sumpt << "  mass of 4l : " << (ZZCandList[i].first.l1RP4+ZZCandList[i].first.l2RP4+ZZCandList[i].second.l1RP4+ZZCandList[i].second.l2RP4).M() << "  Z1 mass : " << (ZZCandList[i].first.l1RP4+ZZCandList[i].first.l2RP4).M() << "  Z2 mass : " << (ZZCandList[i].second.l1RP4+ZZCandList[i].second.l2RP4).M() << std::endl;
	
	if (Z1mass_diff < BestZZCand_Z1mass_diff) {
          if (verbosity_ > 0) std::cout << "=> ZZ at idx = " << i << " has more sumpt than Z at idx = " << bestZZidx << std::endl;
          BestZZCand = ZZCandList[i];
          BestZZCand_Z1mass_diff = Z1mass_diff;
          bestZZidx = i;
        }
	else if (Z1mass_diff == BestZZCand_Z1mass_diff) {
	  if (Z2lep_scalarsum > BestZZCand_Z2lep_scalarsum){
          if (verbosity_ > 0) std::cout << "=> ZZ at idx = " << i << " has more sumpt than Z at idx = " << bestZZidx << std::endl;
          BestZZCand = ZZCandList[i];
          BestZZCand_Z2lep_scalarsum = Z2lep_scalarsum;
          bestZZidx = i;
	  }
        }
	/*
        if (sumpt > BestZZCand_sumpt) {
          if (verbosity_ > 0) std::cout << "=> ZZ at idx = " << i << " has more sumpt than Z at idx = " << bestZZidx << std::endl;
          BestZZCand = ZZCandList[i];
          BestZZCand_sumpt = sumpt;
          bestZZidx = i;
	  }*/
      }
    } 

    //std::cout << "best ZZ candidate :   " << BestZZCand.first.l1RP4.Pt() << "\t" << BestZZCand.first.l1RP4.Eta() << "\t" <<  BestZZCand.first.l2RP4.Pt() << "\t" << BestZZCand.first.l2RP4.Eta() << "\t" << BestZZCand.second.l1RP4.Pt() << "\t" << BestZZCand.second.l1RP4.Eta() << "\t" << BestZZCand.second.l2RP4.Pt() << "\t" << BestZZCand.second.l2RP4.Eta() << "\t" << std::endl;

    int type0 = static_cast<int>(setEventType(BestZZCand.first, BestZZCand.second));
    AnaUtil::fillHist1D("nEventsPerType_standard4l", type0);
    TH2F *hz1l1 , *hz1l2, *hz2l1, *hz2l2;
    //std::cout <<  " | Z1 candiate isgap : " << BestZZCand.first.l1isGap << "\t" << BestZZCand.first.l2isGap << std::endl;
    //std::cout <<  " | Z2 candiate isgap : " << BestZZCand.second.l1isGap << "\t" << BestZZCand.second.l2isGap << std::endl;

    //std::cout <<  " | Z1 candiate sceta : " << BestZZCand.first.l1sceta << "\t" << BestZZCand.first.l2sceta << std::endl;
    //std::cout <<  " | Z2 candiate sceta : " << BestZZCand.second.l1sceta << "\t" << BestZZCand.second.l2sceta << std::endl;

    if (isMC()){
      if (BestZZCand.first.l1isGap == 0 ) hz1l1 = ele_scale_2017;
      else if (BestZZCand.first.l1isGap == 1 ) hz1l1 = ele_scale_2017_gap;
      else hz1l1 = mu_scale_2017;
      
      if (BestZZCand.first.l2isGap == 0) hz1l2 = ele_scale_2017;
      else if (BestZZCand.first.l2isGap == 1) hz1l2 = ele_scale_2017_gap;
      else hz1l2 = mu_scale_2017;
      
      if (BestZZCand.second.l1isGap == 0) hz2l1 = ele_scale_2017;
      else if (BestZZCand.second.l1isGap == 1) hz2l1 = ele_scale_2017_gap;
      else hz2l1 = mu_scale_2017;
      
      if (BestZZCand.second.l2isGap == 0) hz2l2 = ele_scale_2017;
      else if (BestZZCand.second.l2isGap == 1) hz2l2 = ele_scale_2017_gap;
      else hz2l2 = mu_scale_2017;
      
      weight_scalefactor *= ScaleFactor(BestZZCand.first , BestZZCand.second, hz1l1, hz1l2, hz2l1, hz2l2);
      
      //  std::cout << "weight_scalefactor = " << weight_scalefactor << "\t" << "puevWt before scalefactor = " << puevWt << std::endl; 
    }
    if (weight_scalefactor > 0.) puevWt = puevWt * weight_scalefactor;

    //std::cout<< "puevWt after scalefactor = " << puevWt << std::endl;

    //Electron Reconstruction efficiency factor
    TH2F *hz1l1_ele , *hz1l2_ele, *hz2l1_ele, *hz2l2_ele;
     
    if (isMC()){
      if (BestZZCand.first.l1RP4.Pt() < 20) hz1l1_ele = ele_Reco_eff_2017_lowEt;
      else if (BestZZCand.first.l1RP4.Pt() >= 20) hz1l1_ele = ele_Reco_eff_2017_highEt;
      
      if (BestZZCand.first.l2RP4.Pt() < 20) hz1l2_ele = ele_Reco_eff_2017_lowEt;
      else if (BestZZCand.first.l2RP4.Pt() >= 20) hz1l2_ele = ele_Reco_eff_2017_highEt;
      
      if (BestZZCand.second.l1RP4.Pt() < 20) hz2l1_ele = ele_Reco_eff_2017_lowEt;
      else if (BestZZCand.second.l1RP4.Pt() >= 20) hz2l1_ele = ele_Reco_eff_2017_highEt;
    
      if (BestZZCand.second.l2RP4.Pt() < 20) hz2l2_ele = ele_Reco_eff_2017_lowEt;
      else if (BestZZCand.second.l2RP4.Pt() >= 20) hz2l2_ele = ele_Reco_eff_2017_highEt;
      
      weight_elerecofactor *=   EleRecoFactor(BestZZCand.first , BestZZCand.second, hz1l1_ele, hz1l2_ele, hz2l1_ele, hz2l2_ele);
    }
    if (weight_elerecofactor > 0.) puevWt = puevWt * weight_elerecofactor;
    
    if (isMC()) AnaUtil::fillHist1D("nEventsPerTypeWt_standard4l", type0, puevWt);    

    AnaUtil::fillHist1D("evtCutFlow", 12);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 12, puevWt);

    
    double mass_4l_SR = (BestZZCand.first.l1RP4 + BestZZCand.first.l2RP4 + BestZZCand.second.l1RP4 + BestZZCand.second.l2RP4).M();
    //std::cout << "mass_4l in SR : " << mass_4l_SR << std::endl;
    //if (mass_4l_SR < 115 || mass_4l_SR > 135) continue; // SR Cut
    AnaUtil::fillHist1D("evtCutFlow", 13);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 13, puevWt);
    
    /*if (BestZZCand.first.flavour == ZType::mumu && BestZZCand.second.flavour == ZType::mumu)
      std::cout << "4mu : After Standard 4l selection : eventId=" << run << ":" << lumis << ":" << event << std::endl;
    else if (BestZZCand.first.flavour == ZType::ee && BestZZCand.second.flavour == ZType::ee)
      std::cout << "4e : After Standard 4l selection : eventId=" << run << ":" << lumis << ":" << event << std::endl;
    else if ((BestZZCand.first.flavour == ZType::mumu && BestZZCand.second.flavour == ZType::ee) || (BestZZCand.first.flavour == ZType::ee && BestZZCand.second.flavour == ZType::mumu)) 
      std::cout << "2e2mu : After Standard 4l selection : eventId=" << run << ":" << lumis << ":" << event << std::endl;
    */

    // this piece is is important to understand Z1 flavours and isGap & sceta
    /*
    string flavz1 = "e-mu";
    if (BestZZCand.first.flavour == ZType::mumu) flavz1 = "muon";
    else if (BestZZCand.first.flavour == ZType::ee) flavz1 = "electron";

    if (flavz1 == "electron") std::cout << "Z1 candiate flavours : " << flavz1 << " | Z1 candiate sc etas : " << BestZZCand.first.l1sceta << "\t" << BestZZCand.first.l2sceta << std::endl;
    if (flavz1 == "electron") std::cout << "Z1 candiate flavours : " << flavz1 << " | Z1 candiate isgap : " << BestZZCand.first.l1isGap << "\t" << BestZZCand.first.l2isGap << std::endl;
    */

    TLorentzVector fourlP4 = BestZZCand.first.l1RP4  
                           + BestZZCand.first.l2RP4 
                           + BestZZCand.second.l1RP4 
                           + BestZZCand.second.l2RP4;
    double m4l_bestZZ = fourlP4.M(); 
    AnaUtil::fillHist1D("Z1mass", BestZZCand.first.mass); 
    AnaUtil::fillHist1D("Z2mass", BestZZCand.second.mass); 
    AnaUtil::fillHist1D("fourLmass", m4l_bestZZ);
    AnaUtil::fillHist1D("mass4l", m4l_bestZZ);

    // Do separate jet selection since the cuts are different compared to HZZ
    std::vector<vhtm::Jet> selectedJets;
    int nbJets = 0;
    int nbnegJets = 0;
    //std::cout << "#jets before selection: " << jetColl()->size() << std::endl;
    for (const auto& jt: *jetColl()) {
      AnaUtil::fillHist1D("bjetDiscriminator_alljets", jt.pfCombinedInclusiveSecondaryVertexV2BJetTags);

      //      std::cout << "jet loop started" << std::endl;
      // std::cout << "jet btag score : " << jt.pfCombinedInclusiveSecondaryVertexV2BJetTags << std::endl;
      //      if (jt.pt <= 20.) continue;
      //if (std::fabs(jt.eta) >= 2.4) continue;      
      //  if (std::fabs(jt.eta) > 4.7) continue;
      //if (!HZZ4lUtil::isLooseJet(jt)) continue;

      //std::cout << "jet pt : " << jt.pt << "  jt eta : " << jt.eta << "  pujetid : " << jt.pujetid << std::endl;
      bool jet_cut = jt.pt > 20 && std::fabs(jt.eta) < 2.4 && (jt.pujetid == 1 || jt.pt > 50);
      if (!jet_cut) continue;

      AnaUtil::fillHist1D("bjetDiscriminator_preselectedjets", jt.pfCombinedInclusiveSecondaryVertexV2BJetTags);
      //if (!jetLeptonCleaning(jt, 0.3)) std::cout << "jet is not lepton cleaned" << std::endl; 
      //std::cout << "now ignore the jet-lepton dR dump" << std::endl;
      if (!jetLeptonCleaning(jt, 0.3)) continue;

      AnaUtil::fillHist1D("bjetDiscriminator_cleanjets", jt.pfCombinedInclusiveSecondaryVertexV2BJetTags);
      //      std::cout<< "clean jet btag score : " << jt.pfCombinedInclusiveSecondaryVertexV2BJetTags << std::endl;
      // Medium working point
      if (jt.pfCombinedInclusiveSecondaryVertexV2BJetTags > AnaUtil::cutValue(jetCutMap(), "btagFactor"))
        nbJets++;

      if (jt.pfCombinedInclusiveSecondaryVertexV2BJetTags < 0) nbnegJets++;
      selectedJets.push_back(jt);
    }
   
    double jet1Pt = 0., jet2Pt = 0.;
    size_t nJets = selectedJets.size();
    if (nJets > 0) jet1Pt = selectedJets[0].pt; 
    if (nJets > 1) jet2Pt = selectedJets[1].pt; 

    //std::cout << "#jets : " << nJets << std::endl;
    //if (nJets > 0) std::cout << selectedJets[0].pfCombinedInclusiveSecondaryVertexV2BJetTags << std::endl;
    //if (nJets > 1) std::cout << selectedJets[1].pfCombinedInclusiveSecondaryVertexV2BJetTags << std::endl;

    //    std::cout << "Pikabooh" << std::endl;

    //    HZZ4lUtil::syncDumper(run, lumis, event, BestZZCand.first, BestZZCand.second, nJets,
    //                      jet1Pt, jet2Pt, syncDumpf_);

    if (nJets < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 14);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 14, puevWt);    
    AnaUtil::fillHist1D("mass4l_2j", m4l_bestZZ);

    int type = static_cast<int>(setEventType(BestZZCand.first, BestZZCand.second));
    AnaUtil::fillHist1D("nEventsPerType_4l+2j", type);
    if (isMC()) AnaUtil::fillHist1D("nEventsPerTypeWt_4l+2j", type, puevWt);

    //std::cout << "mass_4l in SR : " << mass_4l_SR << "  " << m4l_bestZZ << std::endl;
    // |mass4l - massH| < 10 GeV cut
    if (std::fabs(m4l_bestZZ - 125) > 10) continue; 
    AnaUtil::fillHist1D("evtCutFlow", 15);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 15, puevWt);    
    AnaUtil::fillHist1D("mass4l_2j_mHc", m4l_bestZZ);

    /*    if (BestZZCand.first.flavour == ZType::mumu && BestZZCand.second.flavour == ZType::mumu)
      std::cout << "4mu : After Standard 4l + 2j selection : eventId=" << run << ":" << lumis << ":" << event << std::endl;
    else if (BestZZCand.first.flavour == ZType::ee && BestZZCand.second.flavour == ZType::ee)
      std::cout << "4e : After Standard 4l + 2j selection : eventId=" << run << ":" << lumis << ":" << event << std::endl;
    else if ((BestZZCand.first.flavour == ZType::mumu && BestZZCand.second.flavour == ZType::ee) || (BestZZCand.first.flavour == ZType::ee && BestZZCand.second.flavour == ZType::mumu)) 
    std::cout << "2e2mu : After Standard 4l + 2j selection : eventId=" << run << ":" << lumis << ":" << event << std::endl;*/
    
    int type1 = static_cast<int>(setEventType(BestZZCand.first, BestZZCand.second));
    AnaUtil::fillHist1D("nEventsPerType_mass4l", type1);
    if (isMC()) AnaUtil::fillHist1D("nEventsPerTypeWt_mass4l", type1, puevWt);

    /*     if (type > -1) {
      // fill histograms for different event categories (e.g mmmm, mmee etc.)
      std::ostringstream zcat;
      zcat << "Z1massC" << setw(1) << type;
      AnaUtil::fillHist1D(zcat.str(), BestZZCand.first.mass, puevWt);
      
      std::ostringstream hcat;
      hcat << "Z2massC" << setw(1) << type;
      AnaUtil::fillHist1D(hcat.str(), BestZZCand.second.mass, puevWt);
      
      std::ostringstream drcat;
      drcat << "dRZ1Z2C" << setw(1) << type;
      AnaUtil::fillHist1D(drcat.str(), BestZZCand.first.p4.DeltaR(BestZZCand.second.p4), puevWt);
      }
    */

    // At this point dump event information for synchronization
     //       HZZ4lUtil::syncDumper(run, lumis, event, BestZZCand.first, BestZZCand.second, nJets,
     //                       jet1Pt, jet2Pt, syncDumpf_);

    // di-jet mass
    //AnaUtil::fillHist1D("jet1Mass", (HZZ4lUtil::getP4(selectedJets[0])).M());
    AnaUtil::fillHist1D("jet1Pt", (HZZ4lUtil::getP4(selectedJets[0])).Pt());
    //    AnaUtil::fillHist1D("jet1Eta", (HZZ4lUtil::getP4(selectedJets[0])).Eta());
    //AnaUtil::fillHist1D("jet1Phi", (HZZ4lUtil::getP4(selectedJets[0])).Phi());
    //AnaUtil::fillHist1D("jet1Energy", (HZZ4lUtil::getP4(selectedJets[0])).E());
    //AnaUtil::fillHist2D("jet1EnergyvsPt", (HZZ4lUtil::getP4(selectedJets[0])).E(), HZZ4lUtil::getP4(selectedJets[0]).Pt());

    //AnaUtil::fillHist1D("jet2Mass", (HZZ4lUtil::getP4(selectedJets[1])).M());
    AnaUtil::fillHist1D("jet2Pt", (HZZ4lUtil::getP4(selectedJets[1])).Pt());
    //AnaUtil::fillHist1D("jet2Eta", (HZZ4lUtil::getP4(selectedJets[1])).Eta());
    //AnaUtil::fillHist1D("jet2Phi", (HZZ4lUtil::getP4(selectedJets[1])).Phi());
    //AnaUtil::fillHist1D("jet2Energy", (HZZ4lUtil::getP4(selectedJets[1])).E());
    // AnaUtil::fillHist2D("jet2EnergyvsPt", (HZZ4lUtil::getP4(selectedJets[1])).E(), HZZ4lUtil::getP4(selectedJets[1]).Pt());

    TLorentzVector diJetP4 = HZZ4lUtil::getP4(selectedJets[0]) + HZZ4lUtil::getP4(selectedJets[1]);
    AnaUtil::fillHist1D("dijetMass", diJetP4.M());
    AnaUtil::fillHist1D("dijetPt", diJetP4.Pt());
    //AnaUtil::fillHist1D("dijetEta", diJetP4.Eta());
    //AnaUtil::fillHist1D("dijetPhi", diJetP4.Phi());
    //AnaUtil::fillHist1D("dijetEnergy", diJetP4.E());
    //AnaUtil::fillHist2D("dijetEnergyvsPt", diJetP4.E(), diJetP4.Pt());

    AnaUtil::fillHist1D("resomass", (fourlP4 + diJetP4).M());
    AnaUtil::fillHist1D("resopt", (fourlP4 + diJetP4).Pt());

    float dR_jj = (HZZ4lUtil::getP4(selectedJets[0])).DeltaR(HZZ4lUtil::getP4(selectedJets[1]));
    float dPhi_jj = (HZZ4lUtil::getP4(selectedJets[0])).DeltaPhi(HZZ4lUtil::getP4(selectedJets[1]));
    float dEta_jj = std::fabs((HZZ4lUtil::getP4(selectedJets[0])).Eta() - (HZZ4lUtil::getP4(selectedJets[1])).Eta());
    float prod_Eta_jj = ((HZZ4lUtil::getP4(selectedJets[0])).Eta())*((HZZ4lUtil::getP4(selectedJets[1])).Eta());

    //    TVector3 hzz_Boost = fourlP4.BoostVector();
  //TVector3 b1_boost = (HZZ4l::getP4(selectedJets[0])).BoostVector();-
    //TVector3 hbb_Boost = diJetP4.BoostVector();
    TLorentzVector j1 = HZZ4lUtil::getP4(selectedJets[0]);
    TLorentzVector j2 = HZZ4lUtil::getP4(selectedJets[1]);
 
    //float b1hbb_hel = -10.;
    //float b2hbb_hel = -10.;
    //float z1hzz_hel = -10.;
    //float z2hzz_hel = -10.;
    //float h_hbb_hel = -10.;
    //float h_hzz_hel = -10.;

    //b1hbb_hel = getHelicity( j1 , hbb_Boost );
    //b2hbb_hel = getHelicity( j2 , hbb_Boost );

    TLorentzVector z1 = BestZZCand.first.p4;
    TLorentzVector z2 = BestZZCand.second.p4;
    //z1hzz_hel = getHelicity( z1 , hzz_Boost );
    //z2hzz_hel = getHelicity( z2 , hzz_Boost );

    TLorentzVector diHiggs =  diJetP4 + fourlP4;
    //TVector3 h_hbb_hzz_Boost = diHiggs.BoostVector();
    
    //h_hbb_hel = getHelicity( diJetP4 , h_hbb_hzz_Boost );
    // h_hzz_hel = getHelicity( fourlP4 , h_hbb_hzz_Boost );
    /*
    AnaUtil::fillHist1D("b1hbb_hel", b1hbb_hel);
    AnaUtil::fillHist1D("b2hbb_hel", b2hbb_hel);
    AnaUtil::fillHist1D("z1hzz_hel", z1hzz_hel);
    AnaUtil::fillHist1D("z2hzz_hel", z2hzz_hel);
    AnaUtil::fillHist1D("h_hbb_hel", h_hbb_hel);
    AnaUtil::fillHist1D("h_hzz_hel", h_hzz_hel);
    */
   
    AnaUtil::fillHist1D("dR_JJ", dR_jj);
    AnaUtil::fillHist1D("dPhi_JJ", dPhi_jj);
    AnaUtil::fillHist1D("dEta_JJ", dEta_jj);
    AnaUtil::fillHist1D("Product_Eta_JJ", prod_Eta_jj);

    AnaUtil::fillHist1D("Met", metColl()->at(0).met);
    AnaUtil::fillHist1D("MetPhi", metColl()->at(0).metphi);
    
    AnaUtil::fillHist1D("jTotalscalarpt", selectedJets[0].pt + selectedJets[1].pt);

    AnaUtil::fillHist1D("nbJets", nbJets);

    AnaUtil::fillHist1D("btag1", selectedJets[0].pfCombinedInclusiveSecondaryVertexV2BJetTags);

    AnaUtil::fillHist1D("btag2", selectedJets[1].pfCombinedInclusiveSecondaryVertexV2BJetTags);

    //    float dR_zz = (HZZ4lUtil::getP4(BestZZCand.first)).DeltaR(HZZ4lUtil::getP4(BestZZCand.second));
    // float dPhi_zz = (HZZ4lUtil::getP4(BestZZCand.first)).DeltaPhi(HZZ4lUtil::getP4(BestZZCand.second));

    float dR_zz = (BestZZCand.first.p4).DeltaR((BestZZCand.second.p4));
    AnaUtil::fillHist1D("dR_ZZ", dR_zz);
    float dPhi_zz = (BestZZCand.first.p4).DeltaPhi((BestZZCand.second.p4));
    AnaUtil::fillHist1D("dPhi_ZZ", dPhi_zz);

    float m4l_PtAll = BestZZCand.first.l1P4.Pt() + BestZZCand.first.l2P4.Pt() + BestZZCand.second.l1P4.Pt() + BestZZCand.second.l2P4.Pt();
    AnaUtil::fillHist1D("m4l_PtAll", m4l_PtAll);
    //float jj_ptAll  = selectedJets[0].pt + selectedJets[1].pt;                                                                              
    float hh_PtDiff = std::fabs(m4l_PtAll - diJetP4.M());
    AnaUtil::fillHist1D("hh_PtDiff", hh_PtDiff);
    float hh_MDiff  = std::fabs(m4l_bestZZ - diJetP4.M());
    AnaUtil::fillHist1D("hh_MDiff", hh_MDiff);
    float scalarPtAll = m4l_PtAll + selectedJets[0].pt + selectedJets[1].pt;
    AnaUtil::fillHist1D("scalarPtAll", scalarPtAll);
    float ZZ_massdiff = std::fabs(BestZZCand.first.mass - BestZZCand.second.mass);
    
    AnaUtil::fillHist1D("ZZ_massdiff", ZZ_massdiff);
    AnaUtil::fillHist1D("dR_h4l_hjj", fourlP4.DeltaR(diJetP4));
    AnaUtil::fillHist1D("dPhi_h4l_hjj", fourlP4.DeltaPhi(diJetP4));

    int nLooseLeptons = getTrueLooseLeptons();
    //    AnaUtil::fillHist1D("nLooseLeptons", nLooseLeptons);

    //    float bTag1 = -999.9;
    // float bTag2 = -999.9;

    //  float bTag1 = selectedJets[0].pfCombinedInclusiveSecondaryVertexV2BJetTags;
    //  float bTag2 = selectedJets[1].pfCombinedInclusiveSecondaryVertexV2BJetTags;

    if (selectedJets[0].pfCombinedInclusiveSecondaryVertexV2BJetTags < 0 || selectedJets[1].pfCombinedInclusiveSecondaryVertexV2BJetTags < 0) continue;

    //     if (selectedJets[0].pfCombinedInclusiveSecondaryVertexV2BJetTags >= 0) bTag1 = selectedJets[0].pfCombinedInclusiveSecondaryVertexV2BJetTags;
    // if (selectedJets[1].pfCombinedInclusiveSecondaryVertexV2BJetTags >= 0) bTag2 = selectedJets[1].pfCombinedInclusiveSecondaryVertexV2BJetTags;

    float bTag1 = selectedJets[0].pfCombinedInclusiveSecondaryVertexV2BJetTags;
    float bTag2 = selectedJets[1].pfCombinedInclusiveSecondaryVertexV2BJetTags;


    // select the H->bb as the two highest b-tag scored jets
    std::vector<vhtm::Jet> selectedJets_csv;

    for (int i=0; i<selectedJets.size(); i++)  selectedJets_csv.push_back(selectedJets[i]);
    
    std::sort(selectedJets_csv.begin(), selectedJets_csv.end(), btagComparator);

    //    std::cout << "btag vector size : " << selectedJets_csv.size() << std::endl;
    //std::cout << "btag scores : " << selectedJets_csv[0].pfCombinedInclusiveSecondaryVertexV2BJetTags << "\t" << selectedJets_csv[1].pfCombinedInclusiveSecondaryVertexV2BJetTags << "\t" << selectedJets_csv[2].pfCombinedInclusiveSecondaryVertexV2BJetTags << "\t" << selectedJets_csv[3].pfCombinedInclusiveSecondaryVertexV2BJetTags << "\t" << std::endl;
    TLorentzVector diJetP4_csv = HZZ4lUtil::getP4(selectedJets_csv[0]) + HZZ4lUtil::getP4(selectedJets_csv[1]);
    AnaUtil::fillHist1D("dijetMass_csv", diJetP4_csv.M());
    AnaUtil::fillHist1D("dijetPt_csv", diJetP4_csv.Pt());
    AnaUtil::fillHist1D("resomass_csv", (fourlP4 + diJetP4_csv).M());
    AnaUtil::fillHist1D("resopt_csv", (fourlP4 + diJetP4_csv).Pt());
    
    float dR_jj_csv = (HZZ4lUtil::getP4(selectedJets_csv[0])).DeltaR(HZZ4lUtil::getP4(selectedJets_csv[1]));
    float dPhi_jj_csv = (HZZ4lUtil::getP4(selectedJets_csv[0])).DeltaPhi(HZZ4lUtil::getP4(selectedJets_csv[1]));
    float dEta_jj_csv = std::fabs((HZZ4lUtil::getP4(selectedJets_csv[0])).Eta() - (HZZ4lUtil::getP4(selectedJets_csv[1])).Eta());
    float prod_Eta_jj_csv = ((HZZ4lUtil::getP4(selectedJets_csv[0])).Eta())*((HZZ4lUtil::getP4(selectedJets_csv[1])).Eta());

    TLorentzVector j1_csv = HZZ4lUtil::getP4(selectedJets_csv[0]);
    TLorentzVector j2_csv = HZZ4lUtil::getP4(selectedJets_csv[1]);

    TLorentzVector diHiggs_csv =  diJetP4_csv + fourlP4;

    AnaUtil::fillHist1D("dR_JJ_csv", dR_jj_csv);
    AnaUtil::fillHist1D("dPhi_JJ_csv", dPhi_jj_csv);
    AnaUtil::fillHist1D("dEta_JJ_csv", dEta_jj_csv);
    AnaUtil::fillHist1D("Product_Eta_JJ_csv", prod_Eta_jj_csv);
    AnaUtil::fillHist1D("jTotalscalarpt_csv", selectedJets_csv[0].pt + selectedJets_csv[1].pt);
    AnaUtil::fillHist1D("btag1_csv", selectedJets_csv[0].pfCombinedInclusiveSecondaryVertexV2BJetTags);
    AnaUtil::fillHist1D("btag2_csv", selectedJets_csv[1].pfCombinedInclusiveSecondaryVertexV2BJetTags);
    float hh_PtDiff_csv = std::fabs(m4l_PtAll - diJetP4_csv.M());
    AnaUtil::fillHist1D("hh_PtDiff_csv", hh_PtDiff_csv);
    float hh_MDiff_csv  = std::fabs(m4l_bestZZ - diJetP4_csv.M());
    AnaUtil::fillHist1D("hh_MDiff_csv", hh_MDiff_csv);
    float scalarPtAll_csv = m4l_PtAll + selectedJets_csv[0].pt + selectedJets_csv[1].pt;
    AnaUtil::fillHist1D("scalarPtAll_csv", scalarPtAll_csv);
    AnaUtil::fillHist1D("dR_h4l_hjj_csv", fourlP4.DeltaR(diJetP4_csv));
    AnaUtil::fillHist1D("dPhi_h4l_hjj_csv", fourlP4.DeltaPhi(diJetP4_csv));

    if (selectedJets_csv[0].pfCombinedInclusiveSecondaryVertexV2BJetTags < 0 || selectedJets_csv[1].pfCombinedInclusiveSecondaryVertexV2BJetTags < 0) continue;
    float bTag1_csv = selectedJets_csv[0].pfCombinedInclusiveSecondaryVertexV2BJetTags;
    float bTag2_csv = selectedJets_csv[1].pfCombinedInclusiveSecondaryVertexV2BJetTags;

    float mvaOutput = -999.9;
    float mvaOutputBkg260 = -999.0;
    float mvaOutputBkg300 = -999.0;
    float mvaOutputBkg400 = -999.0;
    float mvaOutputBkg500 = -999.0;
    float mvaOutputBkg600 = -999.0;
    float mvaOutputBkg700 = -999.0;
    float mvaOutputBkg750 = -999.0;
    float mvaOutputBkg800 = -999.0;
    float mvaOutputBkg900 = -999.0;
    float mvaOutputBkg1000 = -999.0;
    // Fill the flat tree here
    if (createMVATree_ && skimObj_ != nullptr) {
      TreeVariables varList;
      if (mm_ch) {
      varList.z1l1Pt    = BestZZCand.first.l1P4.Pt();
      varList.z1l2Pt    = BestZZCand.first.l2P4.Pt();
      varList.z2l1Pt    = BestZZCand.second.l1P4.Pt();
      varList.z2l2Pt    = BestZZCand.second.l2P4.Pt();
      varList.m4lptAll  = m4l_PtAll;
      varList.scalarptAll  = scalarPtAll;
      varList.zzmassdiff = ZZ_massdiff;
      varList.m4l       = m4l_bestZZ;
      varList.dr_zz     = dR_zz;
      varList.dphi_zz   = dPhi_zz;
      varList.dRh4l_hjj = fourlP4.DeltaR(diJetP4);
      varList.dPhih4l_hjj = fourlP4.DeltaPhi(diJetP4);
      varList.dr_jj     = dR_jj;
      varList.dphi_jj   = dPhi_jj;
      varList.deta_jj   = dEta_jj;
      varList.prod_eta_jj   = prod_Eta_jj;
      varList.dijetpt   = diJetP4.Pt();
      varList.dijetmass   = diJetP4.M();
      varList.hhptdiff  = hh_PtDiff; 
      varList.hhmassdiff  = hh_MDiff;
      varList.Met       = metColl()->at(0).met;
      varList.MetPhi    = metColl()->at(0).metphi;
      varList.resomass  = (fourlP4 + diJetP4).M();
      varList.resopt    = (fourlP4 + diJetP4).Pt();
      //varList.b1hbb_hel  = b1hbb_hel;
      //varList.b2hbb_hel  = b2hbb_hel;
      //varList.z1hzz_hel  = z1hzz_hel;
      //varList.z2hzz_hel  = z2hzz_hel;
      //varList.h_hbb_hel  = h_hbb_hel;
      //varList.h_hzz_hel  = h_hzz_hel;
      // varList.nlooseleptons    = nLooseLeptons;
      varList.nbjets    = nbJets;
      //      varList.jTotalscalarpt = selectedJets[0].pt + selectedJets[1].pt;
      varList.j1Pt      = selectedJets[0].pt;
      varList.j2Pt      = selectedJets[1].pt;
      varList.btag1     = bTag1;
      varList.btag2     = bTag2;
      varList.wtall     = puevWt;
      varList.categ       = 1;

      // jet with highest b-tags
      varList.scalarptAll_csv = scalarPtAll_csv;
      varList.dRh4l_hjj_csv   = fourlP4.DeltaR(diJetP4_csv);
      varList.dPhih4l_hjj_csv = fourlP4.DeltaPhi(diJetP4_csv);
      varList.dr_jj_csv       = dR_jj_csv;
      varList.dphi_jj_csv     = dPhi_jj_csv;
      varList.deta_jj_csv     = dEta_jj_csv;
      varList.prod_eta_jj_csv = prod_Eta_jj_csv;
      varList.dijetpt_csv     = diJetP4_csv.Pt();
      varList.dijetmass_csv   = diJetP4_csv.M();
      varList.hhptdiff_csv  = hh_PtDiff_csv;
      varList.hhmassdiff_csv  = hh_MDiff_csv;
      varList.resomass_csv  = (fourlP4 + diJetP4_csv).M();
      varList.resopt_csv    = (fourlP4 + diJetP4_csv).Pt();
      varList.j1Pt_csv      = selectedJets_csv[0].pt;
      varList.j2Pt_csv      = selectedJets_csv[1].pt;
      varList.btag1_csv     = bTag1_csv;
      varList.btag2_csv     = bTag2_csv;
      }

      else if (ee_ch) {
      varList.z1l1Pt    = BestZZCand.first.l1P4.Pt();
      varList.z1l2Pt    = BestZZCand.first.l2P4.Pt();
      varList.z2l1Pt    = BestZZCand.second.l1P4.Pt();
      varList.z2l2Pt    = BestZZCand.second.l2P4.Pt();
      varList.m4lptAll  = m4l_PtAll;
      varList.scalarptAll  = scalarPtAll;
      varList.zzmassdiff = ZZ_massdiff;
      varList.m4l       = m4l_bestZZ;
      varList.dr_zz     = dR_zz;
      varList.dphi_zz   = dPhi_zz;
      varList.dRh4l_hjj = fourlP4.DeltaR(diJetP4);
      varList.dPhih4l_hjj = fourlP4.DeltaPhi(diJetP4);
      varList.dr_jj     = dR_jj;
      varList.dphi_jj   = dPhi_jj;
      varList.deta_jj   = dEta_jj;
      varList.prod_eta_jj   = prod_Eta_jj;
      varList.dijetpt   = diJetP4.Pt();
      varList.dijetmass   = diJetP4.M();
      varList.hhptdiff  = hh_PtDiff; 
      varList.hhmassdiff  = hh_MDiff;
      varList.Met       = metColl()->at(0).met;
      varList.MetPhi    = metColl()->at(0).metphi;
      varList.resomass  = (fourlP4 + diJetP4).M();
      varList.resopt    = (fourlP4 + diJetP4).Pt();
      //varList.b1hbb_hel  = b1hbb_hel;
      //varList.b2hbb_hel  = b2hbb_hel;
      //varList.z1hzz_hel  = z1hzz_hel;
      //varList.z2hzz_hel  = z2hzz_hel;
      //varList.h_hbb_hel  = h_hbb_hel;
      //varList.h_hzz_hel  = h_hzz_hel;
      // varList.nlooseleptons    = nLooseLeptons;
      varList.nbjets    = nbJets;
      //      varList.jTotalscalarpt = selectedJets[0].pt + selectedJets[1].pt;
      varList.j1Pt      = selectedJets[0].pt;
      varList.j2Pt      = selectedJets[1].pt;
      varList.btag1     = bTag1;
      varList.btag2     = bTag2;
      varList.wtall     = puevWt;
      varList.categ       = 2;

      // jet with highest b-tags
      varList.scalarptAll_csv = scalarPtAll_csv;
      varList.dRh4l_hjj_csv   = fourlP4.DeltaR(diJetP4_csv);
      varList.dr_jj_csv       = dR_jj_csv;
      varList.dphi_jj_csv     = dPhi_jj_csv;
      varList.deta_jj_csv     = dEta_jj_csv;
      varList.prod_eta_jj_csv = prod_Eta_jj_csv;
      varList.dijetpt_csv     = diJetP4_csv.Pt();
      varList.dijetmass_csv   = diJetP4_csv.M();
      varList.hhptdiff_csv  = hh_PtDiff_csv;
      varList.hhmassdiff_csv  = hh_MDiff_csv;
      varList.resomass_csv  = (fourlP4 + diJetP4_csv).M();
      varList.resopt_csv    = (fourlP4 + diJetP4_csv).Pt();
      varList.j1Pt_csv      = selectedJets_csv[0].pt;
      varList.j2Pt_csv      = selectedJets_csv[1].pt;
      varList.btag1_csv     = bTag1_csv;
      varList.btag2_csv     = bTag2_csv;
      }

      else if (me_ch) {
      varList.z1l1Pt    = BestZZCand.first.l1P4.Pt();
      varList.z1l2Pt    = BestZZCand.first.l2P4.Pt();
      varList.z2l1Pt    = BestZZCand.second.l1P4.Pt();
      varList.z2l2Pt    = BestZZCand.second.l2P4.Pt();
      varList.m4lptAll  = m4l_PtAll;
      varList.scalarptAll  = scalarPtAll;
      varList.zzmassdiff = ZZ_massdiff;
      varList.m4l       = m4l_bestZZ;
      varList.dr_zz     = dR_zz;
      varList.dphi_zz   = dPhi_zz;
      varList.dRh4l_hjj = fourlP4.DeltaR(diJetP4);
      varList.dPhih4l_hjj = fourlP4.DeltaPhi(diJetP4);
      varList.dr_jj     = dR_jj;
      varList.dphi_jj   = dPhi_jj;
      varList.deta_jj   = dEta_jj;
      varList.prod_eta_jj   = prod_Eta_jj;
      varList.dijetpt   = diJetP4.Pt();
      varList.dijetmass   = diJetP4.M();
      varList.hhptdiff  = hh_PtDiff; 
      varList.hhmassdiff  = hh_MDiff;
      varList.Met       = metColl()->at(0).met;
      varList.MetPhi    = metColl()->at(0).metphi;
      varList.resomass  = (fourlP4 + diJetP4).M();
      varList.resopt    = (fourlP4 + diJetP4).Pt();
      //varList.b1hbb_hel  = b1hbb_hel;
      //varList.b2hbb_hel  = b2hbb_hel;
      //varList.z1hzz_hel  = z1hzz_hel;
      //varList.z2hzz_hel  = z2hzz_hel;
      //varList.h_hbb_hel  = h_hbb_hel;
      //varList.h_hzz_hel  = h_hzz_hel;
      // varList.nlooseleptons    = nLooseLeptons;
      varList.nbjets    = nbJets;
      //      varList.jTotalscalarpt = selectedJets[0].pt + selectedJets[1].pt;
      varList.j1Pt      = selectedJets[0].pt;
      varList.j2Pt      = selectedJets[1].pt;
      varList.btag1     = bTag1;
      varList.btag2     = bTag2;
      varList.wtall     = puevWt;
      varList.categ       = 3;

      // jet with highest b-tags
      varList.scalarptAll_csv = scalarPtAll_csv;
      varList.dRh4l_hjj_csv   = fourlP4.DeltaR(diJetP4_csv);
      varList.dr_jj_csv       = dR_jj_csv;
      varList.dphi_jj_csv     = dPhi_jj_csv;
      varList.deta_jj_csv     = dEta_jj_csv;
      varList.prod_eta_jj_csv = prod_Eta_jj_csv;
      varList.dijetpt_csv     = diJetP4_csv.Pt();
      varList.dijetmass_csv   = diJetP4_csv.M();
      varList.hhptdiff_csv  = hh_PtDiff_csv;
      varList.hhmassdiff_csv  = hh_MDiff_csv;
      varList.resomass_csv  = (fourlP4 + diJetP4_csv).M();
      varList.resopt_csv    = (fourlP4 + diJetP4_csv).Pt();
      varList.j1Pt_csv      = selectedJets_csv[0].pt;
      varList.j2Pt_csv      = selectedJets_csv[1].pt;
      varList.btag1_csv     = bTag1_csv;
      varList.btag2_csv     = bTag2_csv;
      }
    
      skimObj_->fill(varList);
    }
    else if (readMVA_ && (mvAna_ != nullptr || mvAnaBkg260_ != nullptr || mvAnaBkg300_ != nullptr || mvAnaBkg400_ != nullptr || mvAnaBkg500_ != nullptr || mvAnaBkg600_ != nullptr || mvAnaBkg700_ != nullptr || mvAnaBkg750_ != nullptr || mvAnaBkg800_ != nullptr || mvAnaBkg900_ != nullptr || mvAnaBkg1000_ != nullptr)) {
      InputVariables varList;

      //float mz1 = BestZZCand.first.mass, 
      //mz2 = BestZZCand.second.mass;

      /*varList.z1l1Pt    = BestZZCand.first.l1P4.Pt()/mz1;
      varList.z1l2Pt    = BestZZCand.first.l2P4.Pt()/mz1;
      varList.z2l1Pt    = BestZZCand.second.l1P4.Pt()/mz2;
      varList.z2l2Pt    = BestZZCand.second.l1P4.Pt()/mz2;*/
      varList.m4lptAll  = m4l_PtAll;
      varList.scalarptAll = scalarPtAll;
      varList.zzmassdiff = ZZ_massdiff;
      varList.m4l       = m4l_bestZZ;
      varList.dr_zz     = dR_zz;
      varList.dphi_zz   = dPhi_zz;
      varList.dRh4l_hjj = fourlP4.DeltaR(diJetP4);
      varList.dPhih4l_hjj = fourlP4.DeltaPhi(diJetP4);
      varList.dr_jj     = dR_jj;
      varList.dphi_jj   = dPhi_jj;
      varList.deta_jj   = dEta_jj;
      varList.prod_eta_jj   = prod_Eta_jj;
      varList.j1Pt      = selectedJets[0].pt;
      varList.j2Pt      = selectedJets[1].pt;
      varList.Met    = metColl()->at(0).met;
      varList.MetPhi    = metColl()->at(0).metphi;
      varList.resomass  = (fourlP4 + diJetP4).M();
      varList.resopt  = (fourlP4 + diJetP4).Pt();
      varList.dijetpt   = diJetP4.Pt();
      varList.hhptdiff  = hh_PtDiff;
      //varList.b1hbb_hel    = b1hbb_hel;
      //varList.z2hzz_hel  = z2hzz_hel;
      varList.dijetmass = diJetP4.M();
      varList.hhmassdiff  = hh_MDiff;
      //      varList.nbjets    = nbJets;
      //varList.jTotalscalarpt = selectedJets[0].pt + selectedJets[1].pt;
      //varList.j1Pt      = selectedJets[0].pt/m4l_bestZZ;
      //varList.j2Pt      = selectedJets[1].pt/m4l_bestZZ;
      varList.btag1     = bTag1;
      varList.btag2     = bTag2;
//   varList.btag1     = selectedJets[0].pfCombinedInclusiveSecondaryVertexV2BJetTags;
      //  varList.btag2     = selectedJets[1].pfCombinedInclusiveSecondaryVertexV2BJetTags;

      if (isSig_==1) mvaOutput = mvAna_->evaluate(mvAlgo_, varList);
      //  	std::cout << "***************MVA output value wthnloop is: " << mvaOutput << std::endl;
      else if (isSig_ == 0){
	mvaOutputBkg260 = mvAnaBkg260_->evaluate(mvAlgo_, varList);
	mvaOutputBkg300 = mvAnaBkg300_->evaluate(mvAlgo_, varList);
	mvaOutputBkg400 = mvAnaBkg400_->evaluate(mvAlgo_, varList);
	mvaOutputBkg500 = mvAnaBkg500_->evaluate(mvAlgo_, varList);
	mvaOutputBkg600 = mvAnaBkg600_->evaluate(mvAlgo_, varList);
	mvaOutputBkg700 = mvAnaBkg700_->evaluate(mvAlgo_, varList);
	mvaOutputBkg750 = mvAnaBkg750_->evaluate(mvAlgo_, varList);
	mvaOutputBkg800 = mvAnaBkg800_->evaluate(mvAlgo_, varList);
	mvaOutputBkg900 = mvAnaBkg900_->evaluate(mvAlgo_, varList);
	mvaOutputBkg1000 = mvAnaBkg1000_->evaluate(mvAlgo_, varList);
	
      }

    }

    //    std::cout << "***************MVA output value is: " << mvaOutput << std::endl;
    //std::cout << "***************MVA output value for 1000 is: " << mvaOutputBkg1000 << std::endl;
    if (isSig_==1){
      AnaUtil::fillHist1D("mvaoutput_bdt", mvaOutput, puevWt);  
      if (mm_ch) AnaUtil::fillHist1D("mvaoutput_bdt_4mu", mvaOutput, puevWt);
      else if (ee_ch) AnaUtil::fillHist1D("mvaoutput_bdt_4e", mvaOutput, puevWt);
      else if (me_ch) AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu", mvaOutput, puevWt);
    }
    else if (isSig_==0){
      AnaUtil::fillHist1D("mvaoutput_bdt_bkg260", mvaOutputBkg260, puevWt); 
      AnaUtil::fillHist1D("mvaoutput_bdt_bkg300", mvaOutputBkg300, puevWt); 
      AnaUtil::fillHist1D("mvaoutput_bdt_bkg400", mvaOutputBkg400, puevWt); 
      AnaUtil::fillHist1D("mvaoutput_bdt_bkg500", mvaOutputBkg500, puevWt); 
      AnaUtil::fillHist1D("mvaoutput_bdt_bkg600", mvaOutputBkg600, puevWt); 
      AnaUtil::fillHist1D("mvaoutput_bdt_bkg700", mvaOutputBkg700, puevWt); 
      AnaUtil::fillHist1D("mvaoutput_bdt_bkg750", mvaOutputBkg750, puevWt); 
      AnaUtil::fillHist1D("mvaoutput_bdt_bkg800", mvaOutputBkg800, puevWt); 
      AnaUtil::fillHist1D("mvaoutput_bdt_bkg900", mvaOutputBkg900, puevWt); 
      AnaUtil::fillHist1D("mvaoutput_bdt_bkg1000", mvaOutputBkg1000, puevWt); 

      if (mm_ch) {
	AnaUtil::fillHist1D("mvaoutput_bdt_4mu_bkg260", mvaOutputBkg260, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_4mu_bkg300", mvaOutputBkg300, puevWt);	
	AnaUtil::fillHist1D("mvaoutput_bdt_4mu_bkg400", mvaOutputBkg400, puevWt);	
	AnaUtil::fillHist1D("mvaoutput_bdt_4mu_bkg500", mvaOutputBkg500, puevWt);	
	AnaUtil::fillHist1D("mvaoutput_bdt_4mu_bkg600", mvaOutputBkg600, puevWt);	
	AnaUtil::fillHist1D("mvaoutput_bdt_4mu_bkg700", mvaOutputBkg700, puevWt);	
	AnaUtil::fillHist1D("mvaoutput_bdt_4mu_bkg750", mvaOutputBkg750, puevWt);	
	AnaUtil::fillHist1D("mvaoutput_bdt_4mu_bkg800", mvaOutputBkg800, puevWt);	
	AnaUtil::fillHist1D("mvaoutput_bdt_4mu_bkg900", mvaOutputBkg900, puevWt);	
	AnaUtil::fillHist1D("mvaoutput_bdt_4mu_bkg1000", mvaOutputBkg1000, puevWt);	
      }      
      else if (ee_ch) {
	AnaUtil::fillHist1D("mvaoutput_bdt_4e_bkg260", mvaOutputBkg260, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_4e_bkg300", mvaOutputBkg300, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_4e_bkg400", mvaOutputBkg400, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_4e_bkg500", mvaOutputBkg500, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_4e_bkg600", mvaOutputBkg600, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_4e_bkg700", mvaOutputBkg700, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_4e_bkg750", mvaOutputBkg750, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_4e_bkg800", mvaOutputBkg800, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_4e_bkg900", mvaOutputBkg900, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_4e_bkg1000", mvaOutputBkg1000, puevWt);
      }
      else if (me_ch) {
	AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu_bkg260", mvaOutputBkg260, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu_bkg300", mvaOutputBkg300, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu_bkg400", mvaOutputBkg400, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu_bkg500", mvaOutputBkg500, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu_bkg600", mvaOutputBkg600, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu_bkg700", mvaOutputBkg700, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu_bkg750", mvaOutputBkg750, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu_bkg800", mvaOutputBkg800, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu_bkg900", mvaOutputBkg900, puevWt);
	AnaUtil::fillHist1D("mvaoutput_bdt_2e2mu_bkg1000", mvaOutputBkg1000, puevWt);
       } 
    }

    //    AnaUtil::fillHist1D("mvaoutput_bdtg", mvaOutput);
    //AnaUtil::fillHist1D("mvaoutput_mlp", mvaOutput);

    //    if (mvaOutput < -0.02) continue;//bdt_gg+zz
    //if (mvaOutput < 0.0) continue;//bdt_all
    //    if (mvaOutput < -0.07) continue;//bdt_all
    //    if (mvaOutput < -0.05) continue;//bdtg_gg+zz
    //    if (mvaOutput < -0.30) continue;//bdtg_all
    //    if (mvaOutput < 0.67) continue;//mlp_gg+zz
    //    if (mvaOutput < 0.41) continue;//mlp_all
    AnaUtil::fillHist1D("evtCutFlow", 16);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 16, puevWt);

    //    if (nEle >=4) n_ee_fin++;
    // if (nMuon>=4) n_mm_fin++;
    // if (nMuon>=2 && nEle >=2 ) n_me_fin++;

        // Event type
    int type2 = static_cast<int>(setEventType(BestZZCand.first, BestZZCand.second));
    AnaUtil::fillHist1D("nEventsPerType_tmva", type2);
    if (isMC()) AnaUtil::fillHist1D("nEventsPerTypeWt_tmva", type2, puevWt);

    if (type2 > -1) {
      // fill histograms for different event categories (e.g mmmm, mmee etc.)
      std::ostringstream zcat;
      zcat << "Z1massC" << setw(1) << type2;
      AnaUtil::fillHist1D(zcat.str(), BestZZCand.first.mass, puevWt);
      
      std::ostringstream hcat;
      hcat << "Z2massC" << setw(1) << type2;
      AnaUtil::fillHist1D(hcat.str(), BestZZCand.second.mass, puevWt);
      
      std::ostringstream drcat;
      drcat << "dRZ1Z2C" << setw(1) << type2;
      AnaUtil::fillHist1D(drcat.str(), BestZZCand.first.p4.DeltaR(BestZZCand.second.p4), puevWt);
      }

    //    histf()->cd();
    // histf()->cd("DiHiggs4l");

    // selection of bJets
    //if (nbJets == 0) continue;
    //AnaUtil::fillHist1D("evtCutFlow", 12);
    //if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 12, puevWt);    
    AnaUtil::fillHist1D("mass4l_2j_mHc_1b", m4l_bestZZ);

    // Print only the first n events; n configurable
    if (isMC() && dumpEventCount_ > -1 && ++nEventSel >= dumpEventCount_) continue;

    // Dump information
    //dumpEvent(vz, false, true);
    
    /*
    cout << ">>> "
         << "<nLooseMuon>: " << getLooseMuList().size()
         << ", <nMuon>: " << muPhotonPairList.size()
         << ", <nLooseElectron>: " << getLooseEleList().size()
         << ", <nElectron>: " << elePhotonPairList.size()
      //<< ", <nTau>: " << getTauList().size()
      // << ", <nIsoTau>: " << nTau
         << ", <nLooseJet>: " << looseJetList().size()
         << ", <nTightJet>: " << tightJetList().size()
         << ", <nZCandidate>: " << ZCandList.size()
	 << endl;
    */
    //   std::cout << "=> Best ZZ Candidate Info, index: " << bestZZidx << endl;
    // HZZ4lUtil::printZCandidate(BestZZCand.first, "BestZZCandidate:Z1");
    // HZZ4lUtil::printZCandidate(BestZZCand.second, "BestZZCandidate:Z2");
  }

  fin->Close();
 
     //  std::cout << "gen mu chnnl : " << n_mu_chnl << " | gen el chnl : " << n_el_chnl << " | gen el-mu chnl : " << n_muel_chnl << std::endl;
    //std::cout << "gen mu chnnl after det cuts: " << n_det_mu_chnl << " | gen el chnl after det cuts:" << n_det_el_chnl << " | gen el-mu chnl after det cuts:" << n_det_muel_chnl << std::endl;

    //  std::cout << "muon-muon channel yield: " << n_mm_ini << std::endl;
    //std::cout << "e-e channel yield: " << n_ee_ini <<  std::endl;
    // std::cout << "m-e channel yield: " << n_me_ini <<  std::endl;

}
// -------------------------------------------------------------------------------------------------------------------------
// Among all ZZ pairs, require that:
//  - dRCleaning: dR(eta,phi) > 0.02 between each of the four leptons (to remove ghosts)
//  - lepPtcleaning: the two highest-pT leptons of the four pass pT > 20 and 10 GeV
//  - QCD suppression: m(ll)>4 GeV on all four OS pairs that can be made with the four leptons (regardless of flavour). 
//    FSR photons are not used in computing m(ll); since a qcd-induced low mass dilepton (eg jpsi) may have photons nearby 
//    (e.g. from pi0) and so it's safer to keep the cut on m(ll) alone. define the Z1 as the one with mass closest to the 
//    nominal mZ; require mZ1 > 40 GeV. The other Z is the Z2.
//  - additional "smart cut" to reject 4mu/4e pairs where the alternative pairing looks like a on-shell Z+low-mass ll. 
//    (This is required to avoid a background increase when selecting the 
//    "best candidate" after all cuts; cf. Simon's slides. The requirement is: !( |mZa-mZ| < |mZ1-mZ| && mZb<12), 
//    where Za and Zb are the mass-sorted alternative pairing Z candidates (Za being the one closest to the nominal Z mass). 
//    FSR photons associated to leptons are included in the mZa/mZb computation.
//  - m4lmasscut m(4l) > 70 GeV 
// -------------------------------------------------------------------------------------------------------------------------
void DiHiggs4l::ZZCandidateSelector(const std::vector<ZCandidate>& ZCandList, zzvector& ZZCandList, int verbosity) {
  histf()->cd();
  histf()->cd("DiHiggs4l");

  if (verbosity > 0) std::cout << "ZZ selector combination: " << endl
			    << std::setw(6) << "idx" 
			    << std::setw(6) << "jdx" 
			    << std::setw(6) << "dR"
			    << std::setw(6) << "lPt"
			    << std::setw(6) << "qcd"
			    << std::setw(6) << "m4l"
			    << std::setw(6) << "m4l"
		            << std::setw(6) << "z1z2"
			    << std::endl;
  //  std::cout << "# zcandidates present " << ZCandList.size() << std::endl;

  for (size_t i = 0; i < ZCandList.size(); ++i) {
    const auto& zi = ZCandList[i];
    
    for (size_t j = i+1; j < ZCandList.size(); ++j) {
      AnaUtil::fillHist1D("zzCutFlow", 0);

      const auto& zj = ZCandList[j];

      //  std::cout << "leptons pt are : before dRcleaning  " << zi.l1P4.Pt() << "\t" << zi.l1P4.Eta() << "\t" <<  zi.l2P4.Pt() << "\t" << zi.l2P4.Eta() << "\t" << zj.l1P4.Pt() << "\t" << zj.l1P4.Eta() << "\t" << zj.l2P4.Pt() << "\t" << zj.l2P4.Eta() << "\t" << std::endl;
 
      //      std::cout << "dR_cleaning : " << zi.l1P4.DeltaR(zi.l2P4) << "\t" << zj.l1P4.DeltaR(zj.l2P4) << "\t" << zi.l1P4.DeltaR(zj.l1P4) << "\t" << zi.l1P4.DeltaR(zj.l2P4) << "\t" << zi.l2P4.DeltaR(zj.l1P4) << "\t" << zi.l2P4.DeltaR(zj.l2P4) << std::endl;

      if ( zi.l1RP4.Pt() == zj.l1RP4.Pt() || zi.l1RP4.Pt() == zj.l2RP4.Pt()) continue;
      if ( zi.l2RP4.Pt() == zj.l1RP4.Pt() || zi.l2RP4.Pt() == zj.l2RP4.Pt()) continue; // l1p4

      AnaUtil::fillHist1D("zzCutFlow", 1);


      bool p_dRcleaning = zi.l1RP4.DeltaR(zi.l2RP4) > 0.02
	&& zj.l1RP4.DeltaR(zj.l2RP4) > 0.02
	&& zi.l1RP4.DeltaR(zj.l1RP4) > 0.02
	&& zi.l1RP4.DeltaR(zj.l2RP4) > 0.02
	&& zi.l2RP4.DeltaR(zj.l1RP4) > 0.02 
	&& zi.l2RP4.DeltaR(zj.l2RP4) > 0.02;  //l1p4
      

      // the two highest-pT leptons of the four pass pT > 20 and 10 GeV 
      std::vector<TLorentzVector> lepP4List {zi.l1RP4, zi.l2RP4, zj.l1RP4, zj.l2RP4}; // l1p4

      //std::cout << "leptons pt are : after dRcleaning  " << lepP4List[0].Pt() << "\t" << lepP4List[0].Eta() << "\t" <<  lepP4List[1].Pt() << "\t" << lepP4List[1].Eta() << "\t" << lepP4List[2].Pt() << "\t" << lepP4List[2].Eta() << "\t" << lepP4List[3].Pt() << "\t" << lepP4List[3].Eta() << "\t" << std::endl;
      std::sort(lepP4List.begin(), lepP4List.end(), PtComparatorTL<TLorentzVector>());

      if (0) std::cout << "Lepton Pt:" 
                           << lepP4List[0].Pt() 
                   << "\t" << lepP4List[1].Pt() 
                   << "\t" << lepP4List[2].Pt() 
                   << "\t" << lepP4List[3].Pt() 
                   << std::endl;

      // not chechking zi.l1P4 + zi.l2P4 mass since z mass > 12 GeV already required
      // m4l cut 
      // smart cut --> do we need smart cut for diHiggs
      // Full cut flow
      // Z1, Z2 selection
      // QCD suppression
       
      //std::cout << "zi charges : " << zi.l1Charge << "\t" << zi.l2Charge << "\t" << "  zj charges : " << zj.l1Charge << "\t" << zj.l2Charge << std::endl;
      //std::cout << "qcd mass pairs before cut: " << (zi.l1P4 + zj.l1P4).M() << "\t" << (zi.l2P4 + zj.l2P4).M() << "\t" << (zi.l1P4 + zj.l2P4).M() << "\t" << (zi.l2P4 + zj.l1P4).M() << std::endl;
      //std::cout << "qcd mass pairs within same Z before cut : " << (zi.l1P4 + zi.l2P4).M() << "\t" << (zj.l1P4 + zj.l2P4).M() << std::endl;

      
      // OSOZ : Opposite sign opposite Z, OSSZ : Opposite sign same Z
      bool OSOZ11 = (zi.l1Charge * zj.l1Charge < 0) ? (zi.l1P4 + zj.l1P4).M() > 4. && (zi.l2P4 + zj.l2P4).M() > 4. : false; 
      bool OSOZ12 = (zi.l1Charge * zj.l2Charge < 0) ? (zi.l1P4 + zj.l2P4).M() > 4. && (zi.l2P4 + zj.l1P4).M() > 4. : false;                           
      bool OSSZ = (zi.l1Charge * zi.l2Charge < 0 && zj.l1Charge * zj.l2Charge < 0) ? (zi.l1P4 + zi.l2P4).M() > 4. && (zj.l1P4 + zj.l2P4).M() > 4 : false;
      bool p_qcdSuppress = (OSOZ11 || OSOZ12) && OSSZ;




      //bool p_qcdSuppress = (!(zi.l1Charge == zj.l1Charge)) ? (zi.l1P4 + zj.l1P4).M() > 4. && (zi.l2P4 + zj.l2P4).M() > 4.
      //                                                   : (zi.l1P4 + zj.l2P4).M() > 4. && (zi.l2P4 + zj.l1P4).M() > 4.;
      
      bool p_lepPt = lepP4List[0].Pt() > 20. && lepP4List[1].Pt() > 10.; // Pt criteria for dilepton
      bool p_mzpair = (zi.p4 + zj.p4).M() > 70.;
      bool p_m4l = (zi.l1P4 + zi.l2P4 + zj.l1P4 + zj.l2P4).M() > 70.;
      bool p_m4l_fsr = (zi.l1RP4 + zi.l2RP4 + zj.l1RP4 + zj.l2RP4).M() > 70.;
      bool selZ1Z2 = false, 
           selZ2Z1 = false;

      double mZ1, mZ2, ptZ1, ptZ2;
      if (zi.massDiff < zj.massDiff) {
	mZ1 = zi.mass;
	ptZ1 = (zi.l1RP4 + zi.l2RP4).Pt();
	mZ2 = zj.mass;
	ptZ2 = (zj.l1RP4 + zj.l2RP4).Pt();
	if (zi.mass > 40.) selZ1Z2 = true;
      }
      else if (zj.massDiff < zi.massDiff){
	mZ1 = zj.mass;
	ptZ1 = (zj.l1RP4 + zj.l2RP4).Pt();
	mZ2 = zi.mass;
	ptZ2 = (zi.l1RP4 + zi.l2RP4).Pt();
	if (zj.mass > 40.)  selZ2Z1 = true;
      }

      bool p_zpairing = (selZ1Z2 || selZ2Z1);

      // smart cut

      TLorentzVector finalZP4 = zi.l1RP4 + zi.l2RP4 + zj.l1RP4 + zj.l2RP4; // this is wrong (this is H, not Z)
      TLorentzVector Za, Zb;
      double massZa, massZb;

      // std::cout << "zi charges : " << zi.l1Charge << "\t" << zi.l2Charge << "\t" << "  zj charges : " << zj.l1Charge << "\t" << zj.l2Charge << std::endl;
      
      bool p_2e2mu = ((zi.flavour == ZType::mumu && zj.flavour == ZType::ee) || (zi.flavour == ZType::ee && zj.flavour == ZType::mumu)) ? true : false;
      bool p_4mu = (zi.flavour == ZType::mumu && zj.flavour == ZType::mumu) ? true : false;
      bool p_4e = (zi.flavour == ZType::ee && zj.flavour == ZType::ee) ? true : false;


      if (zi.l1Charge * zj.l1Charge < 0){
	//Za = zi.l1RP4 + zj.l1RP4;
	//Zb = zi.l2RP4 + zj.l2RP4;
	Za = zi.l1RP4 + zj.l1RP4;
	Zb = zi.l2RP4 + zj.l2RP4; 
	if (fabs (Za.M() - ZSpace::MZnominal) < fabs (Zb.M() - ZSpace::MZnominal)){
	    massZa = Za.M();                                                                                                                
	    massZb = Zb.M();
	}
	else  { 
	  massZa = Zb.M();
	  massZb = Za.M();
	}
      }
      else if (zi.l1Charge * zj.l2Charge < 0){
	//	Za = zi.l1RP4 + zj.l2RP4;
	//Zb = zi.l2RP4 + zj.l1RP4;
	Za = zi.l1RP4 + zj.l2RP4;
	Zb = zi.l2RP4 + zj.l1RP4;
	if (fabs (Za.M() - ZSpace::MZnominal) < fabs (Zb.M() - ZSpace::MZnominal)){
	  massZa = Za.M();
	  massZb = Zb.M();
        }
	else {
          massZa = Zb.M();
          massZb = Za.M();
        }
      }

      //std::cout << "massZa : " << massZa << "  massZb : " << massZb  << "\t" << "massZ1 : " << mZ1 << " massZ2 : " << mZ2 << std::endl;
      //std::cout << "m_4l : " << (zi.l1P4 + zi.l2P4 + zj.l1P4 + zj.l2P4).M() << std::endl;
      //std::cout << "m_4l_fsr : " << (zi.l1RP4 + zi.l2RP4 + zj.l1RP4 + zj.l2RP4).M() << std::endl;
      
      //std::cout << "Z1 Pt => " <<  ptZ1 << "\t" << "Z2 Pt => " << ptZ2 << std::endl;

      bool p_smartcut = false;
      if (p_4mu || p_4e) p_smartcut = fabs( massZa - ZSpace::MZnominal ) < fabs ( mZ1 - ZSpace::MZnominal ) && massZb < 12 ; 
      else if (p_2e2mu) p_smartcut = false; 

      if (verbosity > 0)
      std::cout << std::setw(6) << i  
		<< std::setw(6) << j 
		<< std::setw(6) << p_dRcleaning
		<< std::setw(6) << p_lepPt
		<< std::setw(6) << p_qcdSuppress
		<< std::setw(6) << p_mzpair
                << std::setw(6) << p_m4l
                << std::setw(6) << p_zpairing
		<< std::setw(6) << p_smartcut
		<< std::endl;

      

      if (!p_dRcleaning) continue;
      AnaUtil::fillHist1D("zzCutFlow", 2);

      if (!p_lepPt) continue;
      AnaUtil::fillHist1D("zzCutFlow", 3);

      if (!p_qcdSuppress) continue;
      AnaUtil::fillHist1D("zzCutFlow", 4);

      //std::cout << "qcd mass pairs after cut: " << (zi.l1P4 + zj.l1P4).M() << "\t" << (zi.l2P4 + zj.l2P4).M() << "\t" << (zi.l1P4 + zj.l2P4).M() << "\t" << (zi.l2P4 + zj.l1P4).M() << std::endl;

      if (!p_zpairing) continue; 
      AnaUtil::fillHist1D("zzCutFlow", 5);
      

      if (p_smartcut) continue;
      AnaUtil::fillHist1D("zzCutFlow", 6);

      //if (!p_mzpair)  continue;
      if (!p_m4l_fsr)  continue; 
      AnaUtil::fillHist1D("zzCutFlow", 7);
      ((selZ1Z2) ? ZZCandList.push_back(std::make_pair(zi, zj))
                 : ZZCandList.push_back(std::make_pair(zj, zi)));

    }
  }
}

double DiHiggs4l::SFCalc(TH2F* h, double pt, double sc_eta, int type){
  double tmp_pt = pt;
  if (type == 0 && pt > 200) tmp_pt = 199.; // Overflow Bin // Mu SF
  if (type == 1 && pt > 500) tmp_pt = 499.; // Ele no Gap & Gap SF
  //  if (!h) std::cout << "histo not found, do again" << std::endl;
  int biny = h->GetYaxis()->FindBin(tmp_pt);
  //std::cout << "biny : " << biny << std::endl;
  int binx = h->GetXaxis()->FindBin(sc_eta);
  //std::cout << "binx : " << binx << std::endl;
  // std::cout << "SFCalc eta : " << sc_eta << std::endl;
  if ( h->GetBinContent(binx,biny) && h->GetBinContent(binx,biny)>0.) return h->GetBinContent(binx,biny);
  else return 1.0;
}


double DiHiggs4l::ScaleFactor(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, TH2F *hZ1l1, TH2F *hZ1l2, TH2F *hZ2l1, TH2F *hZ2l2){
  if (Z1Cand.flavour == ZType::mumu && Z2Cand.flavour == ZType::mumu) return SFCalc(hZ1l1, Z1Cand.l1P4.Pt(), Z1Cand.l1P4.Eta(), 0)* SFCalc(hZ1l2, Z1Cand.l2P4.Pt(), Z1Cand.l2P4.Eta(), 0) * SFCalc(hZ2l1, Z2Cand.l1P4.Pt(), Z2Cand.l1P4.Eta(), 0) * SFCalc(hZ2l2, Z2Cand.l2P4.Pt(), Z2Cand.l2P4.Eta(), 0);
  else if (Z1Cand.flavour == ZType::ee && Z2Cand.flavour == ZType::ee) return SFCalc(hZ1l1, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta, 1)* SFCalc(hZ1l2, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta, 1) * SFCalc(hZ2l1, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta, 1) * SFCalc(hZ2l2, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta, 1);
  else if (Z1Cand.flavour == ZType::mumu && Z2Cand.flavour == ZType::ee) return SFCalc(hZ1l1, Z1Cand.l1P4.Pt(), Z1Cand.l1P4.Eta(), 0)* SFCalc(hZ1l2, Z1Cand.l2P4.Pt(), Z1Cand.l2P4.Eta(), 0) * SFCalc(hZ2l1, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta, 1) * SFCalc(hZ2l2, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta, 1);
  else if (Z1Cand.flavour == ZType::ee && Z2Cand.flavour == ZType::mumu) return SFCalc(hZ1l1, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta, 1)* SFCalc(hZ1l2, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta, 1) * SFCalc(hZ2l1, Z2Cand.l1P4.Pt(), Z2Cand.l1P4.Eta(), 0) * SFCalc(hZ2l2, Z2Cand.l2P4.Pt(), Z2Cand.l2P4.Eta(), 0);
}

double DiHiggs4l::EleRecoFactor(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, TH2F *hZ1l1, TH2F *hZ1l2, TH2F *hZ2l1, TH2F *hZ2l2){
    if (Z1Cand.flavour == ZType::mumu && Z2Cand.flavour == ZType::mumu) return 1.;
    else if (Z1Cand.flavour == ZType::ee && Z2Cand.flavour == ZType::ee) return SFCalc(hZ1l1, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta, 1)* SFCalc(hZ1l2, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta, 1) * SFCalc(hZ2l1, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta, 1) * SFCalc(hZ2l2, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta, 1);
    else if (Z1Cand.flavour == ZType::mumu && Z2Cand.flavour == ZType::ee) return 1.* SFCalc(hZ2l1, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta, 1) * SFCalc(hZ2l2, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta, 1);
    else if (Z1Cand.flavour == ZType::ee && Z2Cand.flavour == ZType::mumu) return SFCalc(hZ1l1, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta, 1)* SFCalc(hZ1l2, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta, 1) * 1.;
}

/*
double DiHiggs4l::ScaleFactor(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, TH2F *h1, TH2F* h2, TH2F* hGap){
  if (Z1Cand.flavour == ZType::mumu && Z2Cand.flavour == ZType::mumu) return SFCalc(h1, Z1Cand.l1P4.Pt(), Z1Cand.l1P4.Eta())* SFCalc(h1, Z1Cand.l2P4.Pt(), Z1Cand.l2P4.Eta()) * SFCalc(h1, Z2Cand.l1P4.Pt(), Z2Cand.l1P4.Eta()) * SFCalc(h1, Z2Cand.l2P4.Pt(), Z2Cand.l2P4.Eta());
  else if (Z1Cand.flavour == ZType::ee && Z2Cand.flavour == ZType::ee) {
    if (Z1Cand.l1isGap == 0 && Z1Cand.l2isGap == 0 && Z2Cand.l1isGap == 0 && Z2Cand.l2isGap == 0) return SFCalc(h2, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta)* SFCalc(h2, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta) * SFCalc(h2, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta) * SFCalc(h2, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta);
    else if (Z1Cand.l1isGap == 1 && Z1Cand.l2isGap == 1 && Z2Cand.l1isGap == 1 && Z2Cand.l2isGap == 1) return SFCalc(hGap, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta)* SFCalc(hGap, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta) * SFCalc(hGap, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta) * SFCalc(hGap, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta);
    //else if (Z1Cand.l2isGap == 0 && Z2Cand.l1isGap == 0)
    else if (Z1Cand.l1isGap == 1 && Z1Cand.l2isGap == 0 && Z2Cand.l1isGap == 0 && Z2Cand.l2isGap == 0) return SFCalc(hGap, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta)* SFCalc(h2, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta) * SFCalc(h2, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta) * SFCalc(h2, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta);
    else if (Z1Cand.l1isGap == 0 && Z1Cand.l2isGap == 1 && Z2Cand.l1isGap == 0 && Z2Cand.l2isGap == 0) return SFCalc(hGap, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta)* SFCalc(hGap, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta) * SFCalc(h2, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta) * SFCalc(h2, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta);
    else if (Z1Cand.l1isGap == 0 && Z1Cand.l2isGap == 0 && Z2Cand.l1isGap == 1 && Z2Cand.l2isGap == 0) return SFCalc(hGap, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta)* SFCalc(h2, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta) * SFCalc(hGap, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta) * SFCalc(h2, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta);
    else if (Z1Cand.l1isGap == 0 && Z1Cand.l2isGap == 0 && Z2Cand.l1isGap == 0 && Z2Cand.l2isGap == 1) return SFCalc(hGap, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta)* SFCalc(h2, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta) * SFCalc(h2, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta) * SFCalc(hGap, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta);
    else if (Z1Cand.l1isGap == 1 && Z1Cand.l2isGap == 1 && Z2Cand.l1isGap == 0 && Z2Cand.l2isGap == 0) return SFCalc(hGap, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta)* SFCalc(hGap, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta) * SFCalc(h2, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta) * SFCalc(h2, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta);
    else if (Z1Cand.l1isGap == 1 && Z1Cand.l2isGap == 0 && Z2Cand.l1isGap == 0 && Z2Cand.l2isGap == 0) return SFCalc(hGap, Z1Cand.l1P4.Pt(), Z1Cand.l1sceta)* SFCalc(h2, Z1Cand.l2P4.Pt(), Z1Cand.l2sceta) * SFCalc(h2, Z2Cand.l1P4.Pt(), Z2Cand.l1sceta) * SFCalc(h2, Z2Cand.l2P4.Pt(), Z2Cand.l2sceta);    
    




  }
  else if (Z1Cand.flavour == ZType::mumu && Z2Cand.flavour == ZType::ee) return SFCalc(h1, Z1Cand.l1P4.Pt(), Z1Cand.l1P4.Eta())* SFCalc(h1, Z1Cand.l2P4.Pt(), Z1Cand.l2P4.Eta()) * SFCalc(h2, Z2Cand.l1P4.Pt(), Z2Cand.l1P4.Eta()) * SFCalc(h2, Z2Cand.l2P4.Pt(), Z2Cand.l2P4.Eta());
  else if (Z1Cand.flavour == ZType::ee && Z2Cand.flavour == ZType::mumu) return SFCalc(h2, Z1Cand.l1P4.Pt(), Z1Cand.l1P4.Eta())* SFCalc(h2, Z1Cand.l2P4.Pt(), Z1Cand.l2P4.Eta()) * SFCalc(h1, Z2Cand.l1P4.Pt(), Z2Cand.l1P4.Eta()) * SFCalc(h1, Z2Cand.l2P4.Pt(), Z2Cand.l2P4.Eta());
  
} 

*/
DiHiggs4l::EventType DiHiggs4l::setEventType(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand) {
  EventType type = EventType::unkwn;
  if (0) cout << "Z1Cand.flavour: " << static_cast<int>(Z1Cand.flavour)
	      << "Z2Cand.flavour: " << static_cast<int>(Z2Cand.flavour)
	      << endl;
  if (Z1Cand.flavour == ZType::mumu) {
    if      (Z2Cand.flavour == ZType::mumu) type = EventType::mmmm;
    else if (Z2Cand.flavour == ZType::ee)   type = EventType::mmee;
  }
  else if (Z1Cand.flavour == ZType::ee) {
    if      (Z2Cand.flavour == ZType::mumu) type = EventType::eemm;
    else if (Z2Cand.flavour == ZType::ee)   type = EventType::eeee;
  }
  if ((Z1Cand.flavour == ZType::mumu && Z2Cand.flavour == ZType::ee) || (Z1Cand.flavour == ZType::ee && Z2Cand.flavour == ZType::mumu)) 
    type = EventType::em;

  return type;
}
// Function to calculate number of extra tight leptons passing isolation apart from leptons in
// Z-h1 candidate
int DiHiggs4l::findExtraLeptons(const ZCandidate& Z1Cand, const ZCandidate& Z2Cand) {
  int nExtLep = 0;
  for (const auto& mu: getTightIsoMuList()) {
    TLorentzVector muP4(HZZ4lUtil::getP4(mu));
    if (AnaUtil::sameObject(muP4, Z1Cand.l1P4) ||
	AnaUtil::sameObject(muP4, Z1Cand.l2P4) ||
	AnaUtil::sameObject(muP4, Z2Cand.l1P4) ||
	AnaUtil::sameObject(muP4, Z2Cand.l2P4)) continue;
    nExtLep++;
  }
  for (const auto& ele: getTightIsoEleList()) {
    TLorentzVector eleP4(HZZ4lUtil::getP4(ele));
    if (AnaUtil::sameObject(eleP4, Z1Cand.l1P4) ||
	AnaUtil::sameObject(eleP4, Z1Cand.l2P4) ||
	AnaUtil::sameObject(eleP4, Z2Cand.l1P4) ||
	AnaUtil::sameObject(eleP4, Z2Cand.l2P4)) continue;
    nExtLep++;
  }
  return nExtLep;
}
float DiHiggs4l::getHelicity( TLorentzVector& p4 , TVector3& boost ){
  float hel = 1e10;
  p4.Boost( -boost );
  hel = TMath::Cos( p4.Vect().Angle( boost ) );
  return hel;
}

float DiHiggs4l::kfactor_qqZZ_qcd_M(double GENmassZZ, int finalState){
  // finalState=1 : 4e/4mu/4tau
  // finalState=2 : 2e2mu/2mutau/2e2tau

  float k=0.0;

  if (finalState==1) {
    k+=1.23613311013*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
    k+=1.17550314639*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
    k+=1.17044565911*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
    k+=1.03141209689*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
    k+=1.05285574912*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
    k+=1.11287217794*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
    k+=1.13361441158*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
    k+=1.10355603327*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
    k+=1.10053981637*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
    k+=1.10972676811*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
    k+=1.12069120525*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
    k+=1.11589101635*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
    k+=1.13906170314*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
    k+=1.14854594271*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
    k+=1.14616229031*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
    k+=1.14573157789*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
    k+=1.13829430515*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
    k+=1.15521193686*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
    k+=1.13679822698*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
    k+=1.13223956942*(abs(GENmassZZ)>475.0);
  }

  if (finalState==2) {
    k+=1.25094466582*(abs(GENmassZZ)>0.0&&abs(GENmassZZ)<=25.0);
    k+=1.22459455362*(abs(GENmassZZ)>25.0&&abs(GENmassZZ)<=50.0);
    k+=1.19287368979*(abs(GENmassZZ)>50.0&&abs(GENmassZZ)<=75.0);
    k+=1.04597506451*(abs(GENmassZZ)>75.0&&abs(GENmassZZ)<=100.0);
    k+=1.08323413771*(abs(GENmassZZ)>100.0&&abs(GENmassZZ)<=125.0);
    k+=1.09994968030*(abs(GENmassZZ)>125.0&&abs(GENmassZZ)<=150.0);
    k+=1.16698455800*(abs(GENmassZZ)>150.0&&abs(GENmassZZ)<=175.0);
    k+=1.10399053155*(abs(GENmassZZ)>175.0&&abs(GENmassZZ)<=200.0);
    k+=1.10592664340*(abs(GENmassZZ)>200.0&&abs(GENmassZZ)<=225.0);
    k+=1.10690381480*(abs(GENmassZZ)>225.0&&abs(GENmassZZ)<=250.0);
    k+=1.11194928918*(abs(GENmassZZ)>250.0&&abs(GENmassZZ)<=275.0);
    k+=1.13522586553*(abs(GENmassZZ)>275.0&&abs(GENmassZZ)<=300.0);
    k+=1.11895090244*(abs(GENmassZZ)>300.0&&abs(GENmassZZ)<=325.0);
    k+=1.13898508615*(abs(GENmassZZ)>325.0&&abs(GENmassZZ)<=350.0);
    k+=1.15463977506*(abs(GENmassZZ)>350.0&&abs(GENmassZZ)<=375.0);
    k+=1.17341664594*(abs(GENmassZZ)>375.0&&abs(GENmassZZ)<=400.0);
    k+=1.20093349763*(abs(GENmassZZ)>400.0&&abs(GENmassZZ)<=425.0);
    k+=1.18915554919*(abs(GENmassZZ)>425.0&&abs(GENmassZZ)<=450.0);
    k+=1.18546007375*(abs(GENmassZZ)>450.0&&abs(GENmassZZ)<=475.0);
    k+=1.12864505708*(abs(GENmassZZ)>475.0);
  }

  if (k==0.0) return 1.1;
  else return k; // if something goes wrong return inclusive k-factor
}

int DiHiggs4l::findEventCategory(int nleptons, const std::vector<vhtm::Jet>& jetList, int nbjets,
				 const ZCandidate& Z1Cand, const ZCandidate& Z2Cand, bool verbose) {
  int njets = jetList.size();

  TLorentzVector j1P4, j2P4;
  j1P4.SetPtEtaPhiE(0.,0.,0.,0.);
  j2P4.SetPtEtaPhiE(0.,0.,0.,0.);
  if (njets) j1P4 = HZZ4lUtil::getP4(jetList[0]);
  if (njets > 1) j2P4 = HZZ4lUtil::getP4(jetList[1]);

  double djet = -1., mjj = 0.;
  TLorentzVector final4lP4 = Z1Cand.l1P4 + Z1Cand.l2P4 + Z1Cand.fsrPhoP4 + Z2Cand.l1P4 + Z2Cand.l2P4 + Z2Cand.fsrPhoP4;
  if (njets >= 2) {
    mjj = (j1P4 + j2P4).M();
    djet = 0.18 * std::fabs(j1P4.Eta() - j2P4.Eta()) + 1.92E-04 * mjj;
  }
  histf()->cd();
  histf()->cd("DiHiggs4l");
  AnaUtil::fillHist1D("djet", djet, 1);

  int cat = 0;
  if (nleptons == 4 && njets >= 2 && nbjets <= 1 && djet > 0.5)
    cat = 2;
  else if ( (nleptons == 4 &&
        njets >= 2 && DiHiggs4l::hasJetPair(jetList) &&
        final4lP4.Pt() > final4lP4.M() ) ||
      (nleptons == 4 && njets == 2 && nbjets == 2) )
    cat = 4;
  else if (njets <= 2 && nbjets == 0 && nleptons >= 5)
    cat = 3;
  else if ( (njets >= 3 && nbjets >= 1) || nleptons >= 5 )
    cat = 5;
  else if (njets >= 1)
    cat = 1;

  if (verbose) {
    cout << "---- Event Category" << endl;
    cout << "  nlep  njet nbjet   jet1Pt  jet1Eta   jet2Pt  jet2Eta      mjj     4lPt      4lM     djet category"
      << endl;
    cout << setw(6) << nleptons
      << setw(6) << njets
      << setw(6) << nbjets
      << setw(9) << j1P4.Pt()
      << setw(9) << (j1P4.Pt() > 0 ? j1P4.Eta() : 99)
      << setw(9) << j2P4.Pt()
      << setw(9) << (j2P4.Pt() > 0 ? j2P4.Eta() : 99)
      << setw(9) << mjj
      << setw(9) << final4lP4.Pt()
      << setw(9) << final4lP4.M()
      << setw(9) << djet
      << setw(9) << cat
      << endl << endl;
  }
  return cat;
}
bool DiHiggs4l::hasJetPair(const std::vector<vhtm::Jet>& jetList) {
  for (size_t i = 0; i < jetList.size(); ++i) {
    const auto& j1 = jetList[i];
    TLorentzVector j1P4(HZZ4lUtil::getP4(j1));
    for (size_t j = i+1; j < jetList.size(); ++j) {
      const auto& j2 = jetList[j];
      TLorentzVector j2P4(HZZ4lUtil::getP4(j2));
      double mjj = (j1P4 + j2P4).M();
      if (0) cout << "mjj[" << i << ", " << j << "] = " << mjj << endl;
      if ( (std::fabs(j1.eta) < 2.4 && j1.pt > 40.) &&
	   (std::fabs(j2.eta) < 2.4 && j2.pt > 40.) &&
	   (mjj > 60. && mjj < 120.) ) return true;
    }
  }
  return false;
}
void DiHiggs4l::endJob() {
  PhysicsObjSelector::endJob();

  histf()->cd();
  histf()->cd("DiHiggs4l");

  std::vector<std::string> zzlabels {
    "All",
    "Removal of  duplicate entries",
    "#DeltaR cleaning",
    "LeptonPtcut",
    "QCD Suppression",
      "Z1Z2 selected",
      "smart cut",
      "#m_{4l} > 70"
  };
  HZZ4lUtil::showEfficiency("zzCutFlow", zzlabels, "ZZ Selection", "ZZCandidates");

  std::vector<std::string> evlabels {
    "Events processed",
      "Gen Filter",
      "Gen level channel filter",
      "Reco level channel filter", //after detector cuts",
      "Events with > 0 good vertex",
      "Events passing trigger",
      "Detector level channel filter",
      "# SF leptons >= 2",
      "hLepPt >= 20 GeV",
      "# of Z Candidates > 0",
      "mlow < Z mass < mhigh",
      "# ZZ candidate",
      "Standard 4l selection (scale factor)",
      "Signal Region cut : 115 < mass(4l) < 135",
      "# jets >= 2",
      "|mass4l - massH| < 10 GeV",
      //      ""
      "tmva cut"
      //"at least 1 b-jet at medium WP",
      // "",
  };
  HZZ4lUtil::showEfficiency("evtCutFlow", evlabels, "Event Selection");

  double lumiFac = 1.0;
  if (isMC()) {
    lumiFac = lumiWt(evtWeightSum_);
    cout << endl
         << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
         << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
         << endl;
  }

  if (isMC()) {
    HZZ4lUtil::scaleHistogram("evtCutFlowWt", lumiFac);
    HZZ4lUtil::showEfficiency("evtCutFlowWt", evlabels, "Event Selection (Weighted)", "Events");
  }

  // Event category
  std::vector<std::string> evCatList {
    "unknwn",
      "mmmm",
      "mmee",
      "eemm",
      "eeee",
      "2e2m"
  };
  HZZ4lUtil::showCount("nEventsPerType", evCatList, "Events Per type");
  HZZ4lUtil::showCount("nEventsPerType_standard4l", evCatList, "Events Per type after standard 4l cut");
  HZZ4lUtil::showCount("nEventsPerType_4l+2j", evCatList, "Events Per type after 4l+2j cut");
  HZZ4lUtil::showCount("nEventsPerType_mass4l", evCatList, "Events Per type after mass4l cut");
  HZZ4lUtil::showCount("nEventsPerType_tmva", evCatList, "Events Per type after tmva cut");
  //  HZZ4lUtil::showCount("nEventsPerType_bjet", evCatList, "Events Per type after bjet cut");

  if (isMC()) {
    HZZ4lUtil::scaleHistogram("nEventsPerTypeWt", lumiFac);
    HZZ4lUtil::scaleHistogram("nEventsPerTypeWt_standard4l", lumiFac);
    HZZ4lUtil::scaleHistogram("nEventsPerTypeWt_4l+2j", lumiFac);
    HZZ4lUtil::scaleHistogram("nEventsPerTypeWt_mass4l", lumiFac);
    HZZ4lUtil::scaleHistogram("nEventsPerTypeWt_tmva", lumiFac);
    //  HZZ4lUtil::scaleHistogram("nEventsPerTypeWt_bjet", lumiFac);                                                                           

    HZZ4lUtil::showCount("nEventsPerTypeWt", evCatList, "Events Per type (Weighted)", 1);
    HZZ4lUtil::showCount("nEventsPerTypeWt_standard4l", evCatList, "Events Per type (Weighted) after standard 4l cut" , 1);
    HZZ4lUtil::showCount("nEventsPerTypeWt_4l+2j", evCatList, "Events Per type (Weighted) after 4l+2j cut" , 1);
    HZZ4lUtil::showCount("nEventsPerTypeWt_mass4l", evCatList, "Events Per type (Weighted) after mass4l cut" , 1);
    HZZ4lUtil::showCount("nEventsPerTypeWt_tmva", evCatList, "Events Per type (Weighted) after tmva cut" , 1);
    //  HZZ4lUtil::showCount("nEventsPerTypeWt_bjet", evCatList, "Events Per type (Weighted) after bjet cut" , 1);

  }
}
void DiHiggs4l::closeFiles() {
  AnaBase::closeFiles();

  // Take care of local stuff first
  if (skimObj_ != nullptr) skimObj_->close();
  if (syncDumpf_.is_open()) syncDumpf_.close();
}
// -------------------------------------------------------------------------------
// Poor man's way of a datacard. Each line between the 'START' and 'END' tags
// is read in turn, split into words, where the first element is the 'key' and
// the rest the value(s). If more than one values are present they are
// stored in a vector. No safety mechanism is in place. Any line with an unknown
// key is skipped. Comments lines should start with either '#' or '//', preferably
// in the first column. Empty lines are skipped. The file containing the datacards
// is passed as the only argument of the program, there is no default
// -------------------------------------------------------------------------------
bool DiHiggs4l::readJob(const string& jobFile, int& nFiles)
{
  if (!PhysicsObjSelector::readJob(jobFile, nFiles)) return false;

  // Open the same file containing the datacards again to read analysis specific cards
  ifstream fin(jobFile.c_str(), std::ios::in);
  if (!fin) {
    cerr << "Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }

  eventFileList_.clear();

  static constexpr int BUF_SIZE = 256;
  char buf[BUF_SIZE];
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;

    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    // Split the line into words
    vector<string> tokens;
    AnaUtil::tokenize(line, tokens);
    assert(tokens.size() > 1);
    const string& key = tokens[0];
    const string& value = tokens[1];
    if (key == "useEventList")
      useEventList_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "verbosity")
      verbosity_ = std::stoi(value.c_str());
    else if (key == "skipDuplicate")
      skipDuplicate_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "dumpGenInfo")
      dumpGenInfo_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "useGenFilter")
      useGenFilter_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "inputEventFile")
      eventFileList_.push_back(value);
    else if (key == "syncDumpFile")
      syncFilename_ = value;
    else if (key == "dumpEventMax")
      dumpEventCount_ = std::stoi(value.c_str());
    else if (key == "selectPartons")
      selectPM_ = (std::stoi(value.c_str()) > 0) ? true : false;
    else if (key == "nMEPartons")
      nMEPartons_ = std::stoi(value.c_str());
    else if (key == "createMVATree")
      createMVATree_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "mvaSkimFile")
      mvaSkimFile_ = value;
    else if (key == "readMVA")
      readMVA_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "sigtype")
      sigtype_ = value;
    else if (key == "isSig")
      isSig_   = std::stoi(value.c_str());
    else if (key == "masspoint")
      masspoint_ = std::stoi(value.c_str());
    //else if (key == "")
    else if (key == "mvaXMLFile")
    mvaXMLFile_ = value;
    else if (key == "mvAlgo")
      mvAlgo_ = value;
  }
  // Close the file
  fin.close();

  if (!syncFilename_.empty()) {
    syncDumpf_.open(syncFilename_, std::ios::out);
    if (!syncDumpf_) {
      cerr << "Output File: " << syncFilename_ << " could not be opened!" << endl;
      return false;
    }
  }
#ifdef SKIP_DUPLICATE_ZMASS
  if (skipDuplicate_ && !eventFilelist_.empty()) {
    eventIdStore_.clear();
    for (const auto& f: eventFilelist_) {
      cout << ">>> Reading file: " << f << endl;
      ifstream fin(f, std::ios::in);
      if (!fin) {
	cerr << "Input file: " << f << " could not be opened!" << endl;
	continue;
      }
      char buf[BUF_SIZE];
      while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
	string line(buf);
	if (line.empty()) continue;   
	
	// enable '#' and '//' style comments
	if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    
	// Split the line into words
	vector<string> tokens;
	AnaUtil::tokenize(line, tokens);
	assert(tokens.size() > 2);
	string key = tokens.at(0) + "-" + tokens.at(1) + "-" + tokens.at(2);
	eventIdStore_.insert({key, 1});
      }
      // Close the file
      fin.close();
    }
    cout << ">>> Total events present: " << eventMap_.size() << endl;
  }  
#endif
  printJob();

  return true;
}
void DiHiggs4l::printJob(ostream& os) const
{
  AnaBase::printJob(os);
  os << endl;
  os << "   useEventList: " << std::boolalpha << useEventList_ << endl
     << "  skipDuplicate: " << std::boolalpha << skipDuplicate_ << endl
     << " dumpEventCount: " << dumpEventCount_ << endl
     << "   syncDumpFile: " << syncFilename_ << endl
     << "   dumpEventMax: " << dumpEventCount_ << endl
     << "  selectPartons: " << std::boolalpha << selectPM_ << endl
     << "     nMEPartons: " << nMEPartons_ << endl
     << "  createMVATree: " << createMVATree_ << endl
     << "        readMVA: " << readMVA_ << endl;
  if (createMVATree_) os << "    mvaSkimFile: " << mvaSkimFile_ << endl;
  if (readMVA_) {
    os << "     mvaXMLFile: " << mvaXMLFile_ << endl;
    os << "         mvAlgo: " << mvAlgo_ << endl;
  }
  if (isMC()) {
    os << "   useGenFilter: " << std::boolalpha << useGenFilter_ << endl;
    os << "    dumpGenInfo: " << std::boolalpha << dumpGenInfo_ << endl;
  }
  AnaUtil::showList(eventFileList_, ">>> INFO. Input event files:", os);
}
