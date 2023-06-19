#ifndef __DiHiggs4l__hh
#define __DiHiggs4l__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <bitset>

#include "TLorentzVector.h"
#include "TVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

#include "PhysicsObjects.h"
#include "AnaUtil.h"
#include "PhysicsObjSelector.h"
#include "ZCandidate.h"

#include "MVASkim.h"
#include "MVAnalysis.h"
#include "GenAnaBase.h"

using zzvector = std::vector<std::pair<ZSpace::ZCandidate, ZSpace::ZCandidate>>;

class DiHiggs4l: public PhysicsObjSelector {
    
public:
  enum class EventType: int {
    unkwn = -1, mmmm = 0, mmee = 1, eemm = 2, eeee = 3, em = 4
  };
  DiHiggs4l();
  virtual ~DiHiggs4l();
  
  virtual void eventLoop() final;  // the main analysis 
  virtual bool beginJob() override;
  virtual void endJob() override;
  virtual bool readJob(const std::string& jobFile, int& nFiles) override;
  virtual void printJob(std::ostream& os=std::cout) const override;
  virtual void bookHistograms() override;
  virtual void closeFiles() override;
    
  void clearLists();
  int findExtraLeptons(const ZSpace::ZCandidate& Z1, const ZSpace::ZCandidate& Z2);
  float getHelicity( TLorentzVector& p4 , TVector3& boost );
  float kfactor_qqZZ_qcd_M(double GENmassZZ, int finalState);
  int findEventCategory(int nleptons, const std::vector<vhtm::Jet>& jetList, int nbjets,
			const ZSpace::ZCandidate& Z1Cand, const ZSpace::ZCandidate& Z2Cand, bool verbose=true);

  double SFCalc(TH2F* h, double pt, double eta, int type);

  double ScaleFactor(const ZSpace::ZCandidate& Z1Cand, const ZSpace::ZCandidate& Z2Cand, TH2F *hZ1l1, TH2F* hZ1l2, TH2F *hZ2l1, TH2F *hZ2l2);

  double EleRecoFactor(const ZSpace::ZCandidate& Z1Cand, const ZSpace::ZCandidate& Z2Cand, TH2F *hZ1l1, TH2F* hZ1l2, TH2F *hZ2l1, TH2F *hZ2l2);

  EventType setEventType(const ZSpace::ZCandidate& Z1Cand, const ZSpace::ZCandidate& Z2Cand);
  static bool hasJetPair(const std::vector<vhtm::Jet>& jetList);
  static bool dmComparator(const ZSpace::ZCandidate& a, const ZSpace::ZCandidate& b) {
    return (a.massDiff < b.massDiff);
  }

  static bool btagComparator(const vhtm::Jet& a, const vhtm::Jet& b) {
    return (a.pfCombinedInclusiveSecondaryVertexV2BJetTags > b.pfCombinedInclusiveSecondaryVertexV2BJetTags);
  }


  static bool massComparator(const ZSpace::ZCandidate& a, const ZSpace::ZCandidate& b) {
    return (a.mass > b.mass);
  }
  void ZZCandidateSelector(const std::vector<ZSpace::ZCandidate>& ZCandList, zzvector& ZZCandList, int verbosity=0);

private:
  std::vector<vhtm::Vertex> vtxList_;

   int verbosity_ {1};
  bool useGenFilter_ {false};
  bool dumpGenInfo_ {false};
  bool useEventList_ {false};
  bool isGluGluToZZ_{false};
  bool skipDuplicate_ {false};
  bool selectPM_ {false};
  int nMEPartons_ {-1};
  ofstream syncDumpf_;
  std::string syncFilename_ {"syncDumpFile.txt"};
  int dumpEventCount_ {0};
  std::vector<std::string> eventFileList_;
  /*  float b1hbb_hel{99.};
  float b2hbb_hel{99.};
  float z1hzz_hel{99.};
  float z2hzz_hel{99.};
  float h_hbb_hel{99.};
  float h_hzz_hel{99.};*/

  double evtWeightSum_ {0};

  bool createMVATree_ {false};
  bool readMVA_ {false};
  std::string mvaSkimFile_ {""};
  std::string sigtype_ {""};
  int isSig_ {-99};
  int masspoint_ {-999};
  std::string mvaXMLFile_ {""};
  std::string mvaXMLFileBkg260_ {""}; 
  std::string mvaXMLFileBkg300_ {""}; 
  std::string mvaXMLFileBkg400_ {""}; 
  std::string mvaXMLFileBkg500_ {""}; 
  std::string mvaXMLFileBkg600_ {""}; 
  std::string mvaXMLFileBkg700_ {""}; 
  std::string mvaXMLFileBkg750_ {""}; 
  std::string mvaXMLFileBkg800_ {""}; 
  std::string mvaXMLFileBkg900_ {""}; 
  std::string mvaXMLFileBkg1000_ {""}; 

  std::string mvAlgoo_ {"BDTG"};

  std::unique_ptr<GenAnaBase> genAna_ {nullptr};
  std::unique_ptr<MVASkim> skimObj_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAna_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAnaBkg260_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAnaBkg300_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAnaBkg400_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAnaBkg500_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAnaBkg600_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAnaBkg700_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAnaBkg750_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAnaBkg800_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAnaBkg900_ {nullptr};
  std::unique_ptr<MVAnalysis> mvAnaBkg1000_ {nullptr};

  std::map<std::string, int> eventIdStore_;
};
#endif
