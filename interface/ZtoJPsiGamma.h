//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 20 21:42:57 2019 by ROOT version 6.10/09
// from TTree Events/Events
// found on file: root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv5/ZToJPsiGamma_TuneCP5_13TeV-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/250000/8EEB3BA0-9E04-AC44-AA78-ED82C143F4E3.root
//////////////////////////////////////////////////////////

#ifndef ZtoJPsiGamma_h
#define ZtoJPsiGamma_h
#include <iostream>
#include<iomanip>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <TString.h>
#include "TLorentzVector.h"
#include "MVASkim.h"
#include "MVAnalysis.h"


using std::ostream;


// Header file for the classes stored in the TTree if any.

class ZtoJPsiGamma {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   Float_t         puWeight;
   Float_t         puWeightUp;
   Float_t         puWeightDown;
   Float_t         HTXS_Higgs_pt;
   Float_t         HTXS_Higgs_y;
   Int_t           HTXS_stage1_1_cat_pTjet25GeV;
   Int_t           HTXS_stage1_1_cat_pTjet30GeV;
   Int_t           HTXS_stage1_1_fine_cat_pTjet25GeV;
   Int_t           HTXS_stage1_1_fine_cat_pTjet30GeV;
   Int_t           HTXS_stage_0;
   Int_t           HTXS_stage_1_pTjet25;
   Int_t           HTXS_stage_1_pTjet30;
   UChar_t         HTXS_njets25;
   UChar_t         HTXS_njets30;
   Float_t         btagWeight_CSVV2;
   Float_t         btagWeight_DeepCSVB;
   Float_t         CaloMET_phi;
   Float_t         CaloMET_pt;
   Float_t         CaloMET_sumEt;
   Float_t         ChsMET_phi;
   Float_t         ChsMET_pt;
   Float_t         ChsMET_sumEt;
   UInt_t          nCorrT1METJet;
   Float_t         CorrT1METJet_area[27];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_eta[27];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_muonSubtrFactor[27];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_phi[27];   //[nCorrT1METJet]
   Float_t         CorrT1METJet_rawPt[27];   //[nCorrT1METJet]
   UInt_t          nElectron;
   Float_t         Electron_deltaEtaSC[5];   //[nElectron]
   Float_t         Electron_dr03EcalRecHitSumEt[5];   //[nElectron]
   Float_t         Electron_dr03HcalDepth1TowerSumEt[5];   //[nElectron]
   Float_t         Electron_dr03TkSumPt[5];   //[nElectron]
   Float_t         Electron_dr03TkSumPtHEEP[5];   //[nElectron]
   Float_t         Electron_dxy[5];   //[nElectron]
   Float_t         Electron_dxyErr[5];   //[nElectron]
   Float_t         Electron_dz[5];   //[nElectron]
   Float_t         Electron_dzErr[5];   //[nElectron]
   Float_t         Electron_eCorr[5];   //[nElectron]
   Float_t         Electron_eInvMinusPInv[5];   //[nElectron]
   Float_t         Electron_energyErr[5];   //[nElectron]
   Float_t         Electron_eta[5];   //[nElectron]
   Float_t         Electron_hoe[5];   //[nElectron]
   Float_t         Electron_ip3d[5];   //[nElectron]
   Float_t         Electron_jetPtRelv2[5];   //[nElectron]
   Float_t         Electron_jetRelIso[5];   //[nElectron]
   Float_t         Electron_mass[5];   //[nElectron]
   Float_t         Electron_miniPFRelIso_all[5];   //[nElectron]
   Float_t         Electron_miniPFRelIso_chg[5];   //[nElectron]
   Float_t         Electron_mvaFall17V1Iso[5];   //[nElectron]
   Float_t         Electron_mvaFall17V1noIso[5];   //[nElectron]
   Float_t         Electron_mvaFall17V2Iso[5];   //[nElectron]
   Float_t         Electron_mvaFall17V2noIso[5];   //[nElectron]
   Float_t         Electron_pfRelIso03_all[5];   //[nElectron]
   Float_t         Electron_pfRelIso03_chg[5];   //[nElectron]
   Float_t         Electron_phi[5];   //[nElectron]
   Float_t         Electron_pt[5];   //[nElectron]
   Float_t         Electron_r9[5];   //[nElectron]
   Float_t         Electron_sieie[5];   //[nElectron]
   Float_t         Electron_sip3d[5];   //[nElectron]
   Float_t         Electron_mvaTTH[5];   //[nElectron]
   Int_t           Electron_charge[5];   //[nElectron]
   Int_t           Electron_cutBased[5];   //[nElectron]
   Int_t           Electron_cutBased_Fall17_V1[5];   //[nElectron]
   Int_t           Electron_jetIdx[5];   //[nElectron]
   Int_t           Electron_pdgId[5];   //[nElectron]
   Int_t           Electron_photonIdx[5];   //[nElectron]
   Int_t           Electron_tightCharge[5];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmap[5];   //[nElectron]
   Bool_t          Electron_convVeto[5];   //[nElectron]
   Bool_t          Electron_cutBased_HEEP[5];   //[nElectron]
   Bool_t          Electron_isPFcand[5];   //[nElectron]
   UChar_t         Electron_lostHits[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WP80[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WP90[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WPL[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WP80[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WP90[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WPL[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP80[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP90[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WPL[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP80[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP90[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WPL[5];   //[nElectron]
   UChar_t         Electron_seedGain[5];   //[nElectron]
   Bool_t          Flag_ecalBadCalibFilterV2;
   UInt_t          nFatJet;
   Float_t         FatJet_area[3];   //[nFatJet]
   Float_t         FatJet_btagCMVA[3];   //[nFatJet]
   Float_t         FatJet_btagCSVV2[3];   //[nFatJet]
   Float_t         FatJet_btagDDBvL[3];   //[nFatJet]
   Float_t         FatJet_btagDDCvB[3];   //[nFatJet]
   Float_t         FatJet_btagDDCvL[3];   //[nFatJet]
   Float_t         FatJet_btagDeepB[3];   //[nFatJet]
   Float_t         FatJet_btagHbb[3];   //[nFatJet]
   Float_t         FatJet_deepTagMD_H4qvsQCD[3];   //[nFatJet]
   Float_t         FatJet_deepTagMD_HbbvsQCD[3];   //[nFatJet]
   Float_t         FatJet_deepTagMD_TvsQCD[3];   //[nFatJet]
   Float_t         FatJet_deepTagMD_WvsQCD[3];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZHbbvsQCD[3];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZHccvsQCD[3];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZbbvsQCD[3];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ZvsQCD[3];   //[nFatJet]
   Float_t         FatJet_deepTagMD_bbvsLight[3];   //[nFatJet]
   Float_t         FatJet_deepTagMD_ccvsLight[3];   //[nFatJet]
   Float_t         FatJet_deepTag_H[3];   //[nFatJet]
   Float_t         FatJet_deepTag_QCD[3];   //[nFatJet]
   Float_t         FatJet_deepTag_QCDothers[3];   //[nFatJet]
   Float_t         FatJet_deepTag_TvsQCD[3];   //[nFatJet]
   Float_t         FatJet_deepTag_WvsQCD[3];   //[nFatJet]
   Float_t         FatJet_deepTag_ZvsQCD[3];   //[nFatJet]
   Float_t         FatJet_eta[3];   //[nFatJet]
   Float_t         FatJet_mass[3];   //[nFatJet]
   Float_t         FatJet_msoftdrop[3];   //[nFatJet]
   Float_t         FatJet_n2b1[3];   //[nFatJet]
   Float_t         FatJet_n3b1[3];   //[nFatJet]
   Float_t         FatJet_phi[3];   //[nFatJet]
   Float_t         FatJet_pt[3];   //[nFatJet]
   Float_t         FatJet_rawFactor[3];   //[nFatJet]
   Float_t         FatJet_tau1[3];   //[nFatJet]
   Float_t         FatJet_tau2[3];   //[nFatJet]
   Float_t         FatJet_tau3[3];   //[nFatJet]
   Float_t         FatJet_tau4[3];   //[nFatJet]
   Int_t           FatJet_jetId[3];   //[nFatJet]
   Int_t           FatJet_subJetIdx1[3];   //[nFatJet]
   Int_t           FatJet_subJetIdx2[3];   //[nFatJet]
   UInt_t          nGenJetAK8;
   Float_t         GenJetAK8_eta[5];   //[nGenJetAK8]
   Float_t         GenJetAK8_mass[5];   //[nGenJetAK8]
   Float_t         GenJetAK8_phi[5];   //[nGenJetAK8]
   Float_t         GenJetAK8_pt[5];   //[nGenJetAK8]
   UInt_t          nGenJet;
   Float_t         GenJet_eta[16];   //[nGenJet]
   Float_t         GenJet_mass[16];   //[nGenJet]
   Float_t         GenJet_phi[16];   //[nGenJet]
   Float_t         GenJet_pt[16];   //[nGenJet]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[82];   //[nGenPart]
   Float_t         GenPart_mass[82];   //[nGenPart]
   Float_t         GenPart_phi[82];   //[nGenPart]
   Float_t         GenPart_pt[82];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[82];   //[nGenPart]
   Int_t           GenPart_pdgId[82];   //[nGenPart]
   Int_t           GenPart_status[82];   //[nGenPart]
   Int_t           GenPart_statusFlags[82];   //[nGenPart]
   UInt_t          nSubGenJetAK8;
   Float_t         SubGenJetAK8_eta[10];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_mass[10];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_phi[10];   //[nSubGenJetAK8]
   Float_t         SubGenJetAK8_pt[10];   //[nSubGenJetAK8]
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   UInt_t          nGenVisTau;
   Float_t         GenVisTau_eta[2];   //[nGenVisTau]
   Float_t         GenVisTau_mass[2];   //[nGenVisTau]
   Float_t         GenVisTau_phi[2];   //[nGenVisTau]
   Float_t         GenVisTau_pt[2];   //[nGenVisTau]
   Int_t           GenVisTau_charge[2];   //[nGenVisTau]
   Int_t           GenVisTau_genPartIdxMother[2];   //[nGenVisTau]
   Int_t           GenVisTau_status[2];   //[nGenVisTau]
   Float_t         genWeight;
   UInt_t          nPSWeight;
   Float_t         PSWeight[1];   //[nPSWeight]
   UInt_t          nIsoTrack;
   Float_t         IsoTrack_dxy[4];   //[nIsoTrack]
   Float_t         IsoTrack_dz[4];   //[nIsoTrack]
   Float_t         IsoTrack_eta[4];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_all[4];   //[nIsoTrack]
   Float_t         IsoTrack_pfRelIso03_chg[4];   //[nIsoTrack]
   Float_t         IsoTrack_phi[4];   //[nIsoTrack]
   Float_t         IsoTrack_pt[4];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_all[4];   //[nIsoTrack]
   Float_t         IsoTrack_miniPFRelIso_chg[4];   //[nIsoTrack]
   Int_t           IsoTrack_fromPV[4];   //[nIsoTrack]
   Int_t           IsoTrack_pdgId[4];   //[nIsoTrack]
   Bool_t          IsoTrack_isHighPurityTrack[4];   //[nIsoTrack]
   Bool_t          IsoTrack_isPFcand[4];   //[nIsoTrack]
   Bool_t          IsoTrack_isFromLostTrack[4];   //[nIsoTrack]
   UInt_t          nJet;
   Float_t         Jet_area[28];   //[nJet]
   Float_t         Jet_btagCMVA[28];   //[nJet]
   Float_t         Jet_btagCSVV2[28];   //[nJet]
   Float_t         Jet_btagDeepB[28];   //[nJet]
   Float_t         Jet_btagDeepC[28];   //[nJet]
   Float_t         Jet_btagDeepFlavB[28];   //[nJet]
   Float_t         Jet_btagDeepFlavC[28];   //[nJet]
   Float_t         Jet_chEmEF[28];   //[nJet]
   Float_t         Jet_chHEF[28];   //[nJet]
   Float_t         Jet_eta[28];   //[nJet]
   Float_t         Jet_jercCHF[28];   //[nJet]
   Float_t         Jet_jercCHPUF[28];   //[nJet]
   Float_t         Jet_mass[28];   //[nJet]
   Float_t         Jet_muEF[28];   //[nJet]
   Float_t         Jet_muonSubtrFactor[28];   //[nJet]
   Float_t         Jet_neEmEF[28];   //[nJet]
   Float_t         Jet_neHEF[28];   //[nJet]
   Float_t         Jet_phi[28];   //[nJet]
   Float_t         Jet_pt[28];   //[nJet]
   Float_t         Jet_qgl[28];   //[nJet]
   Float_t         Jet_rawFactor[28];   //[nJet]
   Float_t         Jet_bRegCorr[28];   //[nJet]
   Float_t         Jet_bRegRes[28];   //[nJet]
   Int_t           Jet_electronIdx1[28];   //[nJet]
   Int_t           Jet_electronIdx2[28];   //[nJet]
   Int_t           Jet_jetId[28];   //[nJet]
   Int_t           Jet_muonIdx1[28];   //[nJet]
   Int_t           Jet_muonIdx2[28];   //[nJet]
   Int_t           Jet_nConstituents[28];   //[nJet]
   Int_t           Jet_nElectrons[28];   //[nJet]
   Int_t           Jet_nMuons[28];   //[nJet]
   Int_t           Jet_puId[28];   //[nJet]
   Float_t         L1PreFiringWeight_Dn;
   Float_t         L1PreFiringWeight_Nom;
   Float_t         L1PreFiringWeight_Up;
   Float_t         METFixEE2017_MetUnclustEnUpDeltaX;
   Float_t         METFixEE2017_MetUnclustEnUpDeltaY;
   Float_t         METFixEE2017_covXX;
   Float_t         METFixEE2017_covXY;
   Float_t         METFixEE2017_covYY;
   Float_t         METFixEE2017_phi;
   Float_t         METFixEE2017_pt;
   Float_t         METFixEE2017_significance;
   Float_t         METFixEE2017_sumEt;
   Float_t         GenMET_phi;
   Float_t         GenMET_pt;
   Float_t         MET_MetUnclustEnUpDeltaX;
   Float_t         MET_MetUnclustEnUpDeltaY;
   Float_t         MET_covXX;
   Float_t         MET_covXY;
   Float_t         MET_covYY;
   Float_t         MET_phi;
   Float_t         MET_pt;
   Float_t         MET_significance;
   Float_t         MET_sumEt;
   UInt_t          nMuon;
   Float_t         Muon_dxy[7];   //[nMuon]
   Float_t         Muon_dxyErr[7];   //[nMuon]
   Float_t         Muon_dz[7];   //[nMuon]
   Float_t         Muon_dzErr[7];   //[nMuon]
   Float_t         Muon_eta[7];   //[nMuon]
   Float_t         Muon_ip3d[7];   //[nMuon]
   Float_t         Muon_jetPtRelv2[7];   //[nMuon]
   Float_t         Muon_jetRelIso[7];   //[nMuon]
   Float_t         Muon_mass[7];   //[nMuon]
   Float_t         Muon_miniPFRelIso_all[7];   //[nMuon]
   Float_t         Muon_miniPFRelIso_chg[7];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[7];   //[nMuon]
   Float_t         Muon_pfRelIso03_chg[7];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[7];   //[nMuon]
   Float_t         Muon_phi[7];   //[nMuon]
   Float_t         Muon_pt[7];   //[nMuon]
   Float_t         Muon_ptErr[7];   //[nMuon]
   Float_t         Muon_segmentComp[7];   //[nMuon]
   Float_t         Muon_sip3d[7];   //[nMuon]
   Float_t         Muon_tkRelIso[7];   //[nMuon]
   Float_t         Muon_tunepRelPt[7];   //[nMuon]
   Float_t         Muon_mvaLowPt[7];   //[nMuon]
   Float_t         Muon_mvaTTH[7];   //[nMuon]
   Int_t           Muon_charge[7];   //[nMuon]
   Int_t           Muon_jetIdx[7];   //[nMuon]
   Int_t           Muon_nStations[7];   //[nMuon]
   Int_t           Muon_nTrackerLayers[7];   //[nMuon]
   Int_t           Muon_pdgId[7];   //[nMuon]
   Int_t           Muon_tightCharge[7];   //[nMuon]
   UChar_t         Muon_highPtId[7];   //[nMuon]
   Bool_t          Muon_inTimeMuon[7];   //[nMuon]
   Bool_t          Muon_isGlobal[7];   //[nMuon]
   Bool_t          Muon_isPFcand[7];   //[nMuon]
   Bool_t          Muon_isTracker[7];   //[nMuon]
   Bool_t          Muon_looseId[7];   //[nMuon]
   Bool_t          Muon_mediumId[7];   //[nMuon]
   Bool_t          Muon_mediumPromptId[7];   //[nMuon]
   UChar_t         Muon_miniIsoId[7];   //[nMuon]
   UChar_t         Muon_multiIsoId[7];   //[nMuon]
   UChar_t         Muon_mvaId[7];   //[nMuon]
   UChar_t         Muon_pfIsoId[7];   //[nMuon]
   Bool_t          Muon_softId[7];   //[nMuon]
   Bool_t          Muon_softMvaId[7];   //[nMuon]
   Bool_t          Muon_tightId[7];   //[nMuon]
   UChar_t         Muon_tkIsoId[7];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[7];   //[nMuon]
   UInt_t          nPhoton;
   Float_t         Photon_eCorr[6];   //[nPhoton]
   Float_t         Photon_energyErr[6];   //[nPhoton]
   Float_t         Photon_eta[6];   //[nPhoton]
   Float_t         Photon_hoe[6];   //[nPhoton]
   Float_t         Photon_mass[6];   //[nPhoton]
   Float_t         Photon_mvaID[6];   //[nPhoton]
   Float_t         Photon_mvaIDV1[6];   //[nPhoton]
   Float_t         Photon_pfRelIso03_all[6];   //[nPhoton]
   Float_t         Photon_pfRelIso03_chg[6];   //[nPhoton]
   Float_t         Photon_phi[6];   //[nPhoton]
   Float_t         Photon_pt[6];   //[nPhoton]
   Float_t         Photon_r9[6];   //[nPhoton]
   Float_t         Photon_sieie[6];   //[nPhoton]
   Int_t           Photon_charge[6];   //[nPhoton]
   Int_t           Photon_cutBasedBitmap[6];   //[nPhoton]
   Int_t           Photon_cutBasedV1Bitmap[6];   //[nPhoton]
   Int_t           Photon_electronIdx[6];   //[nPhoton]
   Int_t           Photon_jetIdx[6];   //[nPhoton]
   Int_t           Photon_pdgId[6];   //[nPhoton]
   Int_t           Photon_vidNestedWPBitmap[6];   //[nPhoton]
   Bool_t          Photon_electronVeto[6];   //[nPhoton]
   Bool_t          Photon_isScEtaEB[6];   //[nPhoton]
   Bool_t          Photon_isScEtaEE[6];   //[nPhoton]
   Bool_t          Photon_mvaID_WP80[6];   //[nPhoton]
   Bool_t          Photon_mvaID_WP90[6];   //[nPhoton]
   Bool_t          Photon_pixelSeed[6];   //[nPhoton]
   UChar_t         Photon_seedGain[6];   //[nPhoton]
   Float_t         Photon_dEsigmaUp[6];
   Float_t         Photon_dEsigmaDown[6];
   Float_t         Pileup_nTrueInt;
   Float_t         Pileup_pudensity;
   Float_t         Pileup_gpudensity;
   Int_t           Pileup_nPU;
   Int_t           Pileup_sumEOOT;
   Int_t           Pileup_sumLOOT;
   Float_t         PuppiMET_phi;
   Float_t         PuppiMET_pt;
   Float_t         PuppiMET_sumEt;
   Float_t         RawMET_phi;
   Float_t         RawMET_pt;
   Float_t         RawMET_sumEt;
   Float_t         fixedGridRhoFastjetAll;
   Float_t         fixedGridRhoFastjetCentral;
   Float_t         fixedGridRhoFastjetCentralCalo;
   Float_t         fixedGridRhoFastjetCentralChargedPileUp;
   Float_t         fixedGridRhoFastjetCentralNeutral;
   UInt_t          nGenDressedLepton;
   Float_t         GenDressedLepton_eta[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_mass[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_phi[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_pt[4];   //[nGenDressedLepton]
   Int_t           GenDressedLepton_pdgId[4];   //[nGenDressedLepton]
   Bool_t          GenDressedLepton_hasTauAnc[4];   //[nGenDressedLepton]
   UInt_t          nSoftActivityJet;
   Float_t         SoftActivityJet_eta[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_phi[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJet_pt[6];   //[nSoftActivityJet]
   Float_t         SoftActivityJetHT;
   Float_t         SoftActivityJetHT10;
   Float_t         SoftActivityJetHT2;
   Float_t         SoftActivityJetHT5;
   Int_t           SoftActivityJetNjets10;
   Int_t           SoftActivityJetNjets2;
   Int_t           SoftActivityJetNjets5;
   UInt_t          nSubJet;
   Float_t         SubJet_btagCMVA[4];   //[nSubJet]
   Float_t         SubJet_btagCSVV2[4];   //[nSubJet]
   Float_t         SubJet_btagDeepB[4];   //[nSubJet]
   Float_t         SubJet_eta[4];   //[nSubJet]
   Float_t         SubJet_mass[4];   //[nSubJet]
   Float_t         SubJet_n2b1[4];   //[nSubJet]
   Float_t         SubJet_n3b1[4];   //[nSubJet]
   Float_t         SubJet_phi[4];   //[nSubJet]
   Float_t         SubJet_pt[4];   //[nSubJet]
   Float_t         SubJet_rawFactor[4];   //[nSubJet]
   Float_t         SubJet_tau1[4];   //[nSubJet]
   Float_t         SubJet_tau2[4];   //[nSubJet]
   Float_t         SubJet_tau3[4];   //[nSubJet]
   Float_t         SubJet_tau4[4];   //[nSubJet]
   UInt_t          nTau;
   Float_t         Tau_chargedIso[4];   //[nTau]
   Float_t         Tau_dxy[4];   //[nTau]
   Float_t         Tau_dz[4];   //[nTau]
   Float_t         Tau_eta[4];   //[nTau]
   Float_t         Tau_leadTkDeltaEta[4];   //[nTau]
   Float_t         Tau_leadTkDeltaPhi[4];   //[nTau]
   Float_t         Tau_leadTkPtOverTauPt[4];   //[nTau]
   Float_t         Tau_mass[4];   //[nTau]
   Float_t         Tau_neutralIso[4];   //[nTau]
   Float_t         Tau_phi[4];   //[nTau]
   Float_t         Tau_photonsOutsideSignalCone[4];   //[nTau]
   Float_t         Tau_pt[4];   //[nTau]
   Float_t         Tau_puCorr[4];   //[nTau]
   Float_t         Tau_rawAntiEle[4];   //[nTau]
   Float_t         Tau_rawAntiEle2018[4];   //[nTau]
   Float_t         Tau_rawIso[4];   //[nTau]
   Float_t         Tau_rawIsodR03[4];   //[nTau]
   Float_t         Tau_rawMVAnewDM2017v2[4];   //[nTau]
   Float_t         Tau_rawMVAoldDM[4];   //[nTau]
   Float_t         Tau_rawMVAoldDM2017v1[4];   //[nTau]
   Float_t         Tau_rawMVAoldDM2017v2[4];   //[nTau]
   Float_t         Tau_rawMVAoldDMdR032017v2[4];   //[nTau]
   Int_t           Tau_charge[4];   //[nTau]
   Int_t           Tau_decayMode[4];   //[nTau]
   Int_t           Tau_jetIdx[4];   //[nTau]
   Int_t           Tau_rawAntiEleCat[4];   //[nTau]
   Int_t           Tau_rawAntiEleCat2018[4];   //[nTau]
   UChar_t         Tau_idAntiEle[4];   //[nTau]
   UChar_t         Tau_idAntiEle2018[4];   //[nTau]
   UChar_t         Tau_idAntiMu[4];   //[nTau]
   Bool_t          Tau_idDecayMode[4];   //[nTau]
   Bool_t          Tau_idDecayModeNewDMs[4];   //[nTau]
   UChar_t         Tau_idMVAnewDM2017v2[4];   //[nTau]
   UChar_t         Tau_idMVAoldDM[4];   //[nTau]
   UChar_t         Tau_idMVAoldDM2017v1[4];   //[nTau]
   UChar_t         Tau_idMVAoldDM2017v2[4];   //[nTau]
   UChar_t         Tau_idMVAoldDMdR032017v2[4];   //[nTau]
   Float_t         TkMET_phi;
   Float_t         TkMET_pt;
   Float_t         TkMET_sumEt;
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[26];   //[nTrigObj]
   Float_t         TrigObj_eta[26];   //[nTrigObj]
   Float_t         TrigObj_phi[26];   //[nTrigObj]
   Float_t         TrigObj_l1pt[26];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[26];   //[nTrigObj]
   Float_t         TrigObj_l2pt[26];   //[nTrigObj]
   Int_t           TrigObj_id[26];   //[nTrigObj]
   Int_t           TrigObj_l1iso[26];   //[nTrigObj]
   Int_t           TrigObj_l1charge[26];   //[nTrigObj]
   Int_t           TrigObj_filterBits[26];   //[nTrigObj]
   Int_t           genTtbarId;
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   Float_t         SV_dlen[11];   //[nSV]
   Float_t         SV_dlenSig[11];   //[nSV]
   Float_t         SV_pAngle[11];   //[nSV]
   Int_t           Electron_genPartIdx[5];   //[nElectron]
   UChar_t         Electron_genPartFlav[5];   //[nElectron]
   Int_t           GenJetAK8_partonFlavour[5];   //[nGenJetAK8]
   UChar_t         GenJetAK8_hadronFlavour[5];   //[nGenJetAK8]
   Int_t           GenJet_partonFlavour[16];   //[nGenJet]
   UChar_t         GenJet_hadronFlavour[16];   //[nGenJet]
   Int_t           Jet_genJetIdx[28];   //[nJet]
   Int_t           Jet_hadronFlavour[28];   //[nJet]
   Int_t           Jet_partonFlavour[28];   //[nJet]
   Int_t           Muon_genPartIdx[7];   //[nMuon]
   UChar_t         Muon_genPartFlav[7];   //[nMuon]
   Int_t           Photon_genPartIdx[6];   //[nPhoton]
   UChar_t         Photon_genPartFlav[6];   //[nPhoton]
   Float_t         MET_fiducialGenPhi;
   Float_t         MET_fiducialGenPt;
   UChar_t         Electron_cleanmask[5];   //[nElectron]
   UChar_t         Jet_cleanmask[28];   //[nJet]
   UChar_t         Muon_cleanmask[7];   //[nMuon]
   UChar_t         Photon_cleanmask[6];   //[nPhoton]
   UChar_t         Tau_cleanmask[4];   //[nTau]
   Float_t         SV_chi2[11];   //[nSV]
   Float_t         SV_eta[11];   //[nSV]
   Float_t         SV_mass[11];   //[nSV]
   Float_t         SV_ndof[11];   //[nSV]
   Float_t         SV_phi[11];   //[nSV]
   Float_t         SV_pt[11];   //[nSV]
   Float_t         SV_x[11];   //[nSV]
   Float_t         SV_y[11];   //[nSV]
   Float_t         SV_z[11];   //[nSV]
   Int_t           Tau_genPartIdx[4];   //[nTau]
   UChar_t         Tau_genPartFlav[4];   //[nTau]
   Bool_t          L1simulation_step; 

   //HLT
   //TLorentzVector  TrigObj;
   Bool_t          HLTriggerFirstPath;
   Bool_t          HLT_DoubleMu20_7_Mass0to30_Photon23;
   Bool_t          HLT_Dimuon20_Jpsi;
   Bool_t          HLT_Dimuon25_Jpsi;
   Bool_t          HLT_Mu17_Photon30_IsoCaloId;
   Bool_t          HLT_Mu17_Photon30_CaloIdL_L1ISO;
   Bool_t          HLT_SingleMu24;
   Bool_t          HLT_SingleMu27;
   Bool_t          HLT_DoubleMu4_3_Jpsi;
   Bool_t          HLTriggerFinalPath;
   
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;
   Bool_t          L1_AlwaysTrue;
   Bool_t          L1_BPTX_AND_Ref1_VME;
   Bool_t          L1_BPTX_AND_Ref3_VME;
   Bool_t          L1_BPTX_AND_Ref4_VME;
   Bool_t          L1_BPTX_BeamGas_B1_VME;
   Bool_t          L1_BPTX_BeamGas_B2_VME;
   Bool_t          L1_BPTX_BeamGas_Ref1_VME;
   Bool_t          L1_BPTX_BeamGas_Ref2_VME;
   Bool_t          L1_BPTX_NotOR_VME;
   Bool_t          L1_BPTX_OR_Ref3_VME;
   Bool_t          L1_BPTX_OR_Ref4_VME;
   Bool_t          L1_BPTX_RefAND_VME;
   Bool_t          L1_BptxMinus;
   Bool_t          L1_BptxOR;
   Bool_t          L1_BptxPlus;
   Bool_t          L1_BptxXOR;
   Bool_t          L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;
   Bool_t          L1_DoubleEG6_HTT240er;
   Bool_t          L1_DoubleEG6_HTT250er;
   Bool_t          L1_DoubleEG6_HTT255er;
   Bool_t          L1_DoubleEG6_HTT270er;
   Bool_t          L1_DoubleEG6_HTT300er;
   Bool_t          L1_DoubleEG8er2p6_HTT255er;
   Bool_t          L1_DoubleEG8er2p6_HTT270er;
   Bool_t          L1_DoubleEG8er2p6_HTT300er;
   Bool_t          L1_DoubleEG_15_10;
   Bool_t          L1_DoubleEG_18_17;
   Bool_t          L1_DoubleEG_20_18;
   Bool_t          L1_DoubleEG_22_10;
   Bool_t          L1_DoubleEG_22_12;
   Bool_t          L1_DoubleEG_22_15;
   Bool_t          L1_DoubleEG_23_10;
   Bool_t          L1_DoubleEG_24_17;
   Bool_t          L1_DoubleEG_25_12;
   Bool_t          L1_DoubleEG_25_13;
   Bool_t          L1_DoubleEG_25_14;
   Bool_t          L1_DoubleEG_LooseIso23_10;
   Bool_t          L1_DoubleEG_LooseIso24_10;
   Bool_t          L1_DoubleIsoTau28er2p1;
   Bool_t          L1_DoubleIsoTau30er2p1;
   Bool_t          L1_DoubleIsoTau32er2p1;
   Bool_t          L1_DoubleIsoTau33er2p1;
   Bool_t          L1_DoubleIsoTau34er2p1;
   Bool_t          L1_DoubleIsoTau35er2p1;
   Bool_t          L1_DoubleIsoTau36er2p1;
   Bool_t          L1_DoubleIsoTau38er2p1;
   Bool_t          L1_DoubleJet100er2p3_dEta_Max1p6;
   Bool_t          L1_DoubleJet100er2p7;
   Bool_t          L1_DoubleJet112er2p3_dEta_Max1p6;
   Bool_t          L1_DoubleJet112er2p7;
   Bool_t          L1_DoubleJet120er2p7;
   Bool_t          L1_DoubleJet150er2p7;
   Bool_t          L1_DoubleJet30_Mass_Min300_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min320_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min340_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min360_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min380_dEta_Max1p5;
   Bool_t          L1_DoubleJet30_Mass_Min400_Mu10;
   Bool_t          L1_DoubleJet30_Mass_Min400_Mu6;
   Bool_t          L1_DoubleJet30_Mass_Min400_dEta_Max1p5;
   Bool_t          L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450;
   Bool_t          L1_DoubleJet40er2p7;
   Bool_t          L1_DoubleJet50er2p7;
   Bool_t          L1_DoubleJet60er2p7;
   Bool_t          L1_DoubleJet60er2p7_ETM100;
   Bool_t          L1_DoubleJet60er2p7_ETM60;
   Bool_t          L1_DoubleJet60er2p7_ETM70;
   Bool_t          L1_DoubleJet60er2p7_ETM80;
   Bool_t          L1_DoubleJet60er2p7_ETM90;
   Bool_t          L1_DoubleJet80er2p7;
   Bool_t          L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;
   Bool_t          L1_DoubleJet_100_35_DoubleJet35_Mass_Min620;
   Bool_t          L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;
   Bool_t          L1_DoubleJet_110_40_DoubleJet40_Mass_Min620;
   Bool_t          L1_DoubleJet_115_35_DoubleJet35_Mass_Min620;
   Bool_t          L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;
   Bool_t          L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;
   Bool_t          L1_DoubleLooseIsoEG22er2p1;
   Bool_t          L1_DoubleLooseIsoEG24er2p1;
   Bool_t          L1_DoubleMu0;
   Bool_t          L1_DoubleMu0_ETM40;
   Bool_t          L1_DoubleMu0_ETM55;
   Bool_t          L1_DoubleMu0_ETM60;
   Bool_t          L1_DoubleMu0_ETM65;
   Bool_t          L1_DoubleMu0_ETM70;
   Bool_t          L1_DoubleMu0_SQ;
   Bool_t          L1_DoubleMu0_SQ_OS;
   Bool_t          L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er1p4_dEta_Max1p8_OS;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS;
   Bool_t          L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;
   Bool_t          L1_DoubleMu0er1p5_SQ_dR_Max1p4;
   Bool_t          L1_DoubleMu0er2_SQ_dR_Max1p4;
   Bool_t          L1_DoubleMu18er2p1;
   Bool_t          L1_DoubleMu22er2p1;
   Bool_t          L1_DoubleMu3_OS_DoubleEG7p5Upsilon;
   Bool_t          L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30;
   Bool_t          L1_DoubleMu3_SQ_HTT100er;
   Bool_t          L1_DoubleMu3_SQ_HTT200er;
   Bool_t          L1_DoubleMu3_SQ_HTT220er;
   Bool_t          L1_DoubleMu3_SQ_HTT240er;
   Bool_t          L1_DoubleMu4_OS_EG12;
   Bool_t          L1_DoubleMu4_SQ_OS;
   Bool_t          L1_DoubleMu4_SQ_OS_dR_Max1p2;
   Bool_t          L1_DoubleMu4p5_SQ;
   Bool_t          L1_DoubleMu4p5_SQ_OS;
   Bool_t          L1_DoubleMu4p5_SQ_OS_dR_Max1p2;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS;
   Bool_t          L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18;
   Bool_t          L1_DoubleMu5Upsilon_OS_DoubleEG3;
   Bool_t          L1_DoubleMu5_OS_EG12;
   Bool_t          L1_DoubleMu5_SQ_OS;
   Bool_t          L1_DoubleMu5_SQ_OS_Mass7to18;
   Bool_t          L1_DoubleMu6_SQ_OS;
   Bool_t          L1_DoubleMu7_EG7;
   Bool_t          L1_DoubleMu7_SQ_EG7;
   Bool_t          L1_DoubleMu8_SQ;
   Bool_t          L1_DoubleMu_10_0_dEta_Max1p8;
   Bool_t          L1_DoubleMu_11_4;
   Bool_t          L1_DoubleMu_12_5;
   Bool_t          L1_DoubleMu_12_8;
   Bool_t          L1_DoubleMu_13_6;
   Bool_t          L1_DoubleMu_15_5;
   Bool_t          L1_DoubleMu_15_5_SQ;
   Bool_t          L1_DoubleMu_15_7;
   Bool_t          L1_DoubleMu_15_7_SQ;
   Bool_t          L1_DoubleMu_15_7_SQ_Mass_Min4;
   Bool_t          L1_DoubleMu_20_2_SQ_Mass_Max20;
   Bool_t          L1_DoubleTau50er2p1;
   Bool_t          L1_DoubleTau70er2p1;
   Bool_t          L1_EG25er2p1_HTT125er;
   Bool_t          L1_EG27er2p1_HTT200er;
   Bool_t          L1_ETM100;
   Bool_t          L1_ETM100_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM105;
   Bool_t          L1_ETM110;
   Bool_t          L1_ETM110_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM115;
   Bool_t          L1_ETM120;
   Bool_t          L1_ETM150;
   Bool_t          L1_ETM30;
   Bool_t          L1_ETM40;
   Bool_t          L1_ETM50;
   Bool_t          L1_ETM60;
   Bool_t          L1_ETM70;
   Bool_t          L1_ETM75;
   Bool_t          L1_ETM75_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM80;
   Bool_t          L1_ETM80_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM85;
   Bool_t          L1_ETM90;
   Bool_t          L1_ETM90_Jet60_dPhi_Min0p4;
   Bool_t          L1_ETM95;
   Bool_t          L1_ETMHF100;
   Bool_t          L1_ETMHF100_HTT60er;
   Bool_t          L1_ETMHF100_Jet60_OR_DiJet30woTT28;
   Bool_t          L1_ETMHF100_Jet60_OR_DoubleJet30;
   Bool_t          L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF110;
   Bool_t          L1_ETMHF110_HTT60er;
   Bool_t          L1_ETMHF110_Jet60_OR_DiJet30woTT28;
   Bool_t          L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF120;
   Bool_t          L1_ETMHF120_HTT60er;
   Bool_t          L1_ETMHF120_Jet60_OR_DiJet30woTT28;
   Bool_t          L1_ETMHF150;
   Bool_t          L1_ETMHF70;
   Bool_t          L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF80;
   Bool_t          L1_ETMHF80_HTT60er;
   Bool_t          L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETMHF90;
   Bool_t          L1_ETMHF90_HTT60er;
   Bool_t          L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30;
   Bool_t          L1_ETT100_BptxAND;
   Bool_t          L1_ETT110_BptxAND;
   Bool_t          L1_ETT40_BptxAND;
   Bool_t          L1_ETT50_BptxAND;
   Bool_t          L1_ETT60_BptxAND;
   Bool_t          L1_ETT70_BptxAND;
   Bool_t          L1_ETT75_BptxAND;
   Bool_t          L1_ETT80_BptxAND;
   Bool_t          L1_ETT85_BptxAND;
   Bool_t          L1_ETT90_BptxAND;
   Bool_t          L1_ETT95_BptxAND;
   Bool_t          L1_FirstBunchAfterTrain;
   Bool_t          L1_FirstBunchInTrain;
   Bool_t          L1_FirstCollisionInOrbit;
   Bool_t          L1_FirstCollisionInTrain;
   Bool_t          L1_HTT120er;
   Bool_t          L1_HTT160er;
   Bool_t          L1_HTT200er;
   Bool_t          L1_HTT220er;
   Bool_t          L1_HTT240er;
   Bool_t          L1_HTT250er_QuadJet_70_55_40_35_er2p5;
   Bool_t          L1_HTT255er;
   Bool_t          L1_HTT270er;
   Bool_t          L1_HTT280er;
   Bool_t          L1_HTT280er_QuadJet_70_55_40_35_er2p5;
   Bool_t          L1_HTT300er;
   Bool_t          L1_HTT300er_QuadJet_70_55_40_35_er2p5;
   Bool_t          L1_HTT320er;
   Bool_t          L1_HTT320er_QuadJet_70_55_40_40_er2p4;
   Bool_t          L1_HTT320er_QuadJet_70_55_40_40_er2p5;
   Bool_t          L1_HTT320er_QuadJet_70_55_45_45_er2p5;
   Bool_t          L1_HTT340er;
   Bool_t          L1_HTT340er_QuadJet_70_55_40_40_er2p5;
   Bool_t          L1_HTT340er_QuadJet_70_55_45_45_er2p5;
   Bool_t          L1_HTT380er;
   Bool_t          L1_HTT400er;
   Bool_t          L1_HTT450er;
   Bool_t          L1_HTT500er;
   Bool_t          L1_IsoEG33_Mt40;
   Bool_t          L1_IsoEG33_Mt44;
   Bool_t          L1_IsoEG33_Mt48;
   Bool_t          L1_IsoTau40er_ETM100;
   Bool_t          L1_IsoTau40er_ETM105;
   Bool_t          L1_IsoTau40er_ETM110;
   Bool_t          L1_IsoTau40er_ETM115;
   Bool_t          L1_IsoTau40er_ETM120;
   Bool_t          L1_IsoTau40er_ETM80;
   Bool_t          L1_IsoTau40er_ETM85;
   Bool_t          L1_IsoTau40er_ETM90;
   Bool_t          L1_IsoTau40er_ETM95;
   Bool_t          L1_IsoTau40er_ETMHF100;
   Bool_t          L1_IsoTau40er_ETMHF110;
   Bool_t          L1_IsoTau40er_ETMHF120;
   Bool_t          L1_IsoTau40er_ETMHF80;
   Bool_t          L1_IsoTau40er_ETMHF90;
   Bool_t          L1_IsolatedBunch;
   Bool_t          L1_LastCollisionInTrain;
   Bool_t          L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3;
   Bool_t          L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7;
   Bool_t          L1_LooseIsoEG26er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3;
   Bool_t          L1_LooseIsoEG28er2p1_HTT100er;
   Bool_t          L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3;
   Bool_t          L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3;
   Bool_t          L1_MU20_EG15;
   Bool_t          L1_MinimumBiasHF0_AND_BptxAND;
   Bool_t          L1_MinimumBiasHF0_OR_BptxAND;
   Bool_t          L1_Mu10er2p1_ETM30;
   Bool_t          L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;
   Bool_t          L1_Mu12_EG10;
   Bool_t          L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;
   Bool_t          L1_Mu14er2p1_ETM30;
   Bool_t          L1_Mu15_HTT100er;
   Bool_t          L1_Mu18_HTT100er;
   Bool_t          L1_Mu18_Jet24er2p7;
   Bool_t          L1_Mu18er2p1_IsoTau26er2p1;
   Bool_t          L1_Mu18er2p1_Tau24er2p1;
   Bool_t          L1_Mu20_EG10;
   Bool_t          L1_Mu20_EG17;
   Bool_t          L1_Mu20_LooseIsoEG6;
   Bool_t          L1_Mu20er2p1_IsoTau26er2p1;
   Bool_t          L1_Mu20er2p1_IsoTau27er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau28er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau30er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau32er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau33er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau34er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau35er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau36er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau38er2p1;
   Bool_t          L1_Mu22er2p1_IsoTau40er2p1;
   Bool_t          L1_Mu22er2p1_Tau50er2p1;
   Bool_t          L1_Mu22er2p1_Tau70er2p1;
   Bool_t          L1_Mu23_EG10;
   Bool_t          L1_Mu23_LooseIsoEG10;
   Bool_t          L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4;
   Bool_t          L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4;
   Bool_t          L1_Mu3_Jet30er2p5;
   Bool_t          L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4;
   Bool_t          L1_Mu5_EG15;
   Bool_t          L1_Mu5_EG20;
   Bool_t          L1_Mu5_EG23;
   Bool_t          L1_Mu5_LooseIsoEG18;
   Bool_t          L1_Mu5_LooseIsoEG20;
   Bool_t          L1_Mu6_DoubleEG10;
   Bool_t          L1_Mu6_DoubleEG17;
   Bool_t          L1_Mu6_HTT200er;
   Bool_t          L1_Mu6_HTT240er;
   Bool_t          L1_Mu6_HTT250er;
   Bool_t          L1_Mu7_EG23;
   Bool_t          L1_Mu7_LooseIsoEG20;
   Bool_t          L1_Mu7_LooseIsoEG23;
   Bool_t          L1_Mu8_HTT150er;
   Bool_t          L1_NotBptxOR;
   Bool_t          L1_QuadJet36er2p7_IsoTau52er2p1;
   Bool_t          L1_QuadJet36er2p7_Tau52;
   Bool_t          L1_QuadJet40er2p7;
   Bool_t          L1_QuadJet50er2p7;
   Bool_t          L1_QuadJet60er2p7;
   Bool_t          L1_QuadMu0;
   Bool_t          L1_SingleEG10;
   Bool_t          L1_SingleEG15;
   Bool_t          L1_SingleEG18;
   Bool_t          L1_SingleEG24;
   Bool_t          L1_SingleEG26;
   Bool_t          L1_SingleEG28;
   Bool_t          L1_SingleEG2_BptxAND;
   Bool_t          L1_SingleEG30;
   Bool_t          L1_SingleEG32;
   Bool_t          L1_SingleEG34;
   Bool_t          L1_SingleEG34er2p1;
   Bool_t          L1_SingleEG36;
   Bool_t          L1_SingleEG36er2p1;
   Bool_t          L1_SingleEG38;
   Bool_t          L1_SingleEG38er2p1;
   Bool_t          L1_SingleEG40;
   Bool_t          L1_SingleEG42;
   Bool_t          L1_SingleEG45;
   Bool_t          L1_SingleEG5;
   Bool_t          L1_SingleEG50;
   Bool_t          L1_SingleIsoEG18;
   Bool_t          L1_SingleIsoEG18er2p1;
   Bool_t          L1_SingleIsoEG20;
   Bool_t          L1_SingleIsoEG20er2p1;
   Bool_t          L1_SingleIsoEG22;
   Bool_t          L1_SingleIsoEG22er2p1;
   Bool_t          L1_SingleIsoEG24;
   Bool_t          L1_SingleIsoEG24er2p1;
   Bool_t          L1_SingleIsoEG26;
   Bool_t          L1_SingleIsoEG26er2p1;
   Bool_t          L1_SingleIsoEG28;
   Bool_t          L1_SingleIsoEG28er2p1;
   Bool_t          L1_SingleIsoEG30;
   Bool_t          L1_SingleIsoEG30er2p1;
   Bool_t          L1_SingleIsoEG32;
   Bool_t          L1_SingleIsoEG32er2p1;
   Bool_t          L1_SingleIsoEG33er2p1;
   Bool_t          L1_SingleIsoEG34;
   Bool_t          L1_SingleIsoEG34er2p1;
   Bool_t          L1_SingleIsoEG35;
   Bool_t          L1_SingleIsoEG35er2p1;
   Bool_t          L1_SingleIsoEG36;
   Bool_t          L1_SingleIsoEG36er2p1;
   Bool_t          L1_SingleIsoEG37;
   Bool_t          L1_SingleIsoEG38;
   Bool_t          L1_SingleIsoEG38er2p1;
   Bool_t          L1_SingleIsoEG40;
   Bool_t          L1_SingleIsoEG40er2p1;
   Bool_t          L1_SingleJet120;
   Bool_t          L1_SingleJet120_FWD;
   Bool_t          L1_SingleJet12_BptxAND;
   Bool_t          L1_SingleJet140;
   Bool_t          L1_SingleJet150;
   Bool_t          L1_SingleJet16;
   Bool_t          L1_SingleJet160;
   Bool_t          L1_SingleJet170;
   Bool_t          L1_SingleJet180;
   Bool_t          L1_SingleJet20;
   Bool_t          L1_SingleJet200;
   Bool_t          L1_SingleJet20er2p7_NotBptxOR;
   Bool_t          L1_SingleJet20er2p7_NotBptxOR_3BX;
   Bool_t          L1_SingleJet35;
   Bool_t          L1_SingleJet35_FWD;
   Bool_t          L1_SingleJet35_HFm;
   Bool_t          L1_SingleJet35_HFp;
   Bool_t          L1_SingleJet43er2p7_NotBptxOR_3BX;
   Bool_t          L1_SingleJet46er2p7_NotBptxOR_3BX;
   Bool_t          L1_SingleJet60;
   Bool_t          L1_SingleJet60_FWD;
   Bool_t          L1_SingleJet60_HFm;
   Bool_t          L1_SingleJet60_HFp;
   Bool_t          L1_SingleJet90;
   Bool_t          L1_SingleJet90_FWD;
   Bool_t          L1_SingleMu0_BMTF;
   Bool_t          L1_SingleMu0_EMTF;
   Bool_t          L1_SingleMu0_OMTF;
   Bool_t          L1_SingleMu10_LowQ;
   Bool_t          L1_SingleMu11_LowQ;
   Bool_t          L1_SingleMu12_LowQ_BMTF;
   Bool_t          L1_SingleMu12_LowQ_EMTF;
   Bool_t          L1_SingleMu12_LowQ_OMTF;
   Bool_t          L1_SingleMu14er2p1;
   Bool_t          L1_SingleMu16;
   Bool_t          L1_SingleMu16er2p1;
   Bool_t          L1_SingleMu18;
   Bool_t          L1_SingleMu18er2p1;
   Bool_t          L1_SingleMu20;
   Bool_t          L1_SingleMu20er2p1;
   Bool_t          L1_SingleMu22;
   Bool_t          L1_SingleMu22_BMTF;
   Bool_t          L1_SingleMu22_EMTF;
   Bool_t          L1_SingleMu22_OMTF;
   Bool_t          L1_SingleMu22er2p1;
   Bool_t          L1_SingleMu25;
   Bool_t          L1_SingleMu3;
   Bool_t          L1_SingleMu30;
   Bool_t          L1_SingleMu5;
   Bool_t          L1_SingleMu7;
   Bool_t          L1_SingleMuCosmics;
   Bool_t          L1_SingleMuCosmics_BMTF;
   Bool_t          L1_SingleMuCosmics_EMTF;
   Bool_t          L1_SingleMuCosmics_OMTF;
   Bool_t          L1_SingleMuOpen;
   Bool_t          L1_SingleMuOpen_NotBptxOR;
   Bool_t          L1_SingleMuOpen_NotBptxOR_3BX;
   Bool_t          L1_SingleTau100er2p1;
   Bool_t          L1_SingleTau120er2p1;
   Bool_t          L1_SingleTau130er2p1;
   Bool_t          L1_SingleTau140er2p1;
   Bool_t          L1_SingleTau20;
   Bool_t          L1_SingleTau80er2p1;
   Bool_t          L1_TripleEG_14_10_8;
   Bool_t          L1_TripleEG_18_17_8;
   Bool_t          L1_TripleEG_LooseIso20_10_5;
   Bool_t          L1_TripleJet_100_85_72_VBF;
   Bool_t          L1_TripleJet_105_85_76_VBF;
   Bool_t          L1_TripleJet_84_68_48_VBF;
   Bool_t          L1_TripleJet_88_72_56_VBF;
   Bool_t          L1_TripleJet_92_76_64_VBF;
   Bool_t          L1_TripleJet_98_83_71_VBF;
   Bool_t          L1_TripleMu0;
   Bool_t          L1_TripleMu0_OQ;
   Bool_t          L1_TripleMu3;
   Bool_t          L1_TripleMu3_SQ;
   Bool_t          L1_TripleMu_4_4_4;
   Bool_t          L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ;
   Bool_t          L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;
   Bool_t          L1_TripleMu_5_0_0;
   Bool_t          L1_TripleMu_5_3_3;
   Bool_t          L1_TripleMu_5_3p5_2p5;
   Bool_t          L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;
   Bool_t          L1_TripleMu_5_5_3;
   Bool_t          L1_UnpairedBunchBptxMinus;
   Bool_t          L1_UnpairedBunchBptxPlus;
   Bool_t          L1_ZeroBias;
   Bool_t          L1_ZeroBias_copy;
   Double_t        genEventSumw=0.;
   Double_t        genEventSumw_=0.;

   // List of branches
   TBranch        *b_genEventSumw;
   TBranch        *b_genEventSumw_;
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_puWeight;
   TBranch        *b_puWeightUp;
   TBranch        *b_puWeightDown;
   TBranch        *b_HTXS_Higgs_pt;   //!
   TBranch        *b_HTXS_Higgs_y;   //!
   TBranch        *b_HTXS_stage1_1_cat_pTjet25GeV;   //!
   TBranch        *b_HTXS_stage1_1_cat_pTjet30GeV;   //!
   TBranch        *b_HTXS_stage1_1_fine_cat_pTjet25GeV;   //!
   TBranch        *b_HTXS_stage1_1_fine_cat_pTjet30GeV;   //!
   TBranch        *b_HTXS_stage_0;   //!
   TBranch        *b_HTXS_stage_1_pTjet25;   //!
   TBranch        *b_HTXS_stage_1_pTjet30;   //!
   TBranch        *b_HTXS_njets25;   //!
   TBranch        *b_HTXS_njets30;   //!
   TBranch        *b_btagWeight_CSVV2;   //!
   TBranch        *b_btagWeight_DeepCSVB;   //!
   TBranch        *b_CaloMET_phi;   //!
   TBranch        *b_CaloMET_pt;   //!
   TBranch        *b_CaloMET_sumEt;   //!
   TBranch        *b_ChsMET_phi;   //!
   TBranch        *b_ChsMET_pt;   //!
   TBranch        *b_ChsMET_sumEt;   //!
   TBranch        *b_nCorrT1METJet;   //!
   TBranch        *b_CorrT1METJet_area;   //!
   TBranch        *b_CorrT1METJet_eta;   //!
   TBranch        *b_CorrT1METJet_muonSubtrFactor;   //!
   TBranch        *b_CorrT1METJet_phi;   //!
   TBranch        *b_CorrT1METJet_rawPt;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_Electron_deltaEtaSC;   //!
   TBranch        *b_Electron_dr03EcalRecHitSumEt;   //!
   TBranch        *b_Electron_dr03HcalDepth1TowerSumEt;   //!
   TBranch        *b_Electron_dr03TkSumPt;   //!
   TBranch        *b_Electron_dr03TkSumPtHEEP;   //!
   TBranch        *b_Electron_dxy;   //!
   TBranch        *b_Electron_dxyErr;   //!
   TBranch        *b_Electron_dz;   //!
   TBranch        *b_Electron_dzErr;   //!
   TBranch        *b_Electron_eCorr;   //!
   TBranch        *b_Electron_eInvMinusPInv;   //!
   TBranch        *b_Electron_energyErr;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_hoe;   //!
   TBranch        *b_Electron_ip3d;   //!
   TBranch        *b_Electron_jetPtRelv2;   //!
   TBranch        *b_Electron_jetRelIso;   //!
   TBranch        *b_Electron_mass;   //!
   TBranch        *b_Electron_miniPFRelIso_all;   //!
   TBranch        *b_Electron_miniPFRelIso_chg;   //!
   TBranch        *b_Electron_mvaFall17V1Iso;   //!
   TBranch        *b_Electron_mvaFall17V1noIso;   //!
   TBranch        *b_Electron_mvaFall17V2Iso;   //!
   TBranch        *b_Electron_mvaFall17V2noIso;   //!
   TBranch        *b_Electron_pfRelIso03_all;   //!
   TBranch        *b_Electron_pfRelIso03_chg;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_r9;   //!
   TBranch        *b_Electron_sieie;   //!
   TBranch        *b_Electron_sip3d;   //!
   TBranch        *b_Electron_mvaTTH;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_cutBased;   //!
   TBranch        *b_Electron_cutBased_Fall17_V1;   //!
   TBranch        *b_Electron_jetIdx;   //!
   TBranch        *b_Electron_pdgId;   //!
   TBranch        *b_Electron_photonIdx;   //!
   TBranch        *b_Electron_tightCharge;   //!
   TBranch        *b_Electron_vidNestedWPBitmap;   //!
   TBranch        *b_Electron_convVeto;   //!
   TBranch        *b_Electron_cutBased_HEEP;   //!
   TBranch        *b_Electron_isPFcand;   //!
   TBranch        *b_Electron_lostHits;   //!
   TBranch        *b_Electron_mvaFall17V1Iso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V1Iso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V1Iso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V1noIso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V2Iso_WPL;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WP80;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WP90;   //!
   TBranch        *b_Electron_mvaFall17V2noIso_WPL;   //!
   TBranch        *b_Electron_seedGain;   //!
   TBranch        *b_Photon_dEscaleUp;
   TBranch        *b_Photon_dEscaleDown;
   TBranch        *b_Flag_ecalBadCalibFilterV2;   //!
   TBranch        *b_nFatJet;   //!
   TBranch        *b_FatJet_area;   //!
   TBranch        *b_FatJet_btagCMVA;   //!
   TBranch        *b_FatJet_btagCSVV2;   //!
   TBranch        *b_FatJet_btagDDBvL;   //!
   TBranch        *b_FatJet_btagDDCvB;   //!
   TBranch        *b_FatJet_btagDDCvL;   //!
   TBranch        *b_FatJet_btagDeepB;   //!
   TBranch        *b_FatJet_btagHbb;   //!
   TBranch        *b_FatJet_deepTagMD_H4qvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_HbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_TvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_WvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZHbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZHccvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZbbvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_ZvsQCD;   //!
   TBranch        *b_FatJet_deepTagMD_bbvsLight;   //!
   TBranch        *b_FatJet_deepTagMD_ccvsLight;   //!
   TBranch        *b_FatJet_deepTag_H;   //!
   TBranch        *b_FatJet_deepTag_QCD;   //!
   TBranch        *b_FatJet_deepTag_QCDothers;   //!
   TBranch        *b_FatJet_deepTag_TvsQCD;   //!
   TBranch        *b_FatJet_deepTag_WvsQCD;   //!
   TBranch        *b_FatJet_deepTag_ZvsQCD;   //!
   TBranch        *b_FatJet_eta;   //!
   TBranch        *b_FatJet_mass;   //!
   TBranch        *b_FatJet_msoftdrop;   //!
   TBranch        *b_FatJet_n2b1;   //!
   TBranch        *b_FatJet_n3b1;   //!
   TBranch        *b_FatJet_phi;   //!
   TBranch        *b_FatJet_pt;   //!
   TBranch        *b_FatJet_rawFactor;   //!
   TBranch        *b_FatJet_tau1;   //!
   TBranch        *b_FatJet_tau2;   //!
   TBranch        *b_FatJet_tau3;   //!
   TBranch        *b_FatJet_tau4;   //!
   TBranch        *b_FatJet_jetId;   //!
   TBranch        *b_FatJet_subJetIdx1;   //!
   TBranch        *b_FatJet_subJetIdx2;   //!
   TBranch        *b_nGenJetAK8;   //!
   TBranch        *b_GenJetAK8_eta;   //!
   TBranch        *b_GenJetAK8_mass;   //!
   TBranch        *b_GenJetAK8_phi;   //!
   TBranch        *b_GenJetAK8_pt;   //!
   TBranch        *b_nGenJet;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_mass;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_mass;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_genPartIdxMother;   //!
   TBranch        *b_GenPart_pdgId;   //!
   TBranch        *b_GenPart_status;   //!
   TBranch        *b_GenPart_statusFlags;   //!
   TBranch        *b_nSubGenJetAK8;   //!
   TBranch        *b_SubGenJetAK8_eta;   //!
   TBranch        *b_SubGenJetAK8_mass;   //!
   TBranch        *b_SubGenJetAK8_phi;   //!
   TBranch        *b_SubGenJetAK8_pt;   //!
   TBranch        *b_Generator_binvar;   //!
   TBranch        *b_Generator_scalePDF;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Generator_x1;   //!
   TBranch        *b_Generator_x2;   //!
   TBranch        *b_Generator_xpdf1;   //!
   TBranch        *b_Generator_xpdf2;   //!
   TBranch        *b_Generator_id1;   //!
   TBranch        *b_Generator_id2;   //!
   TBranch        *b_nGenVisTau;   //!
   TBranch        *b_GenVisTau_eta;   //!
   TBranch        *b_GenVisTau_mass;   //!
   TBranch        *b_GenVisTau_phi;   //!
   TBranch        *b_GenVisTau_pt;   //!
   TBranch        *b_GenVisTau_charge;   //!
   TBranch        *b_GenVisTau_genPartIdxMother;   //!
   TBranch        *b_GenVisTau_status;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
   TBranch        *b_nIsoTrack;   //!
   TBranch        *b_IsoTrack_dxy;   //!
   TBranch        *b_IsoTrack_dz;   //!
   TBranch        *b_IsoTrack_eta;   //!
   TBranch        *b_IsoTrack_pfRelIso03_all;   //!
   TBranch        *b_IsoTrack_pfRelIso03_chg;   //!
   TBranch        *b_IsoTrack_phi;   //!
   TBranch        *b_IsoTrack_pt;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_all;   //!
   TBranch        *b_IsoTrack_miniPFRelIso_chg;   //!
   TBranch        *b_IsoTrack_fromPV;   //!
   TBranch        *b_IsoTrack_pdgId;   //!
   TBranch        *b_IsoTrack_isHighPurityTrack;   //!
   TBranch        *b_IsoTrack_isPFcand;   //!
   TBranch        *b_IsoTrack_isFromLostTrack;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_Jet_area;   //!
   TBranch        *b_Jet_btagCMVA;   //!
   TBranch        *b_Jet_btagCSVV2;   //!
   TBranch        *b_Jet_btagDeepB;   //!
   TBranch        *b_Jet_btagDeepC;   //!
   TBranch        *b_Jet_btagDeepFlavB;   //!
   TBranch        *b_Jet_btagDeepFlavC;   //!
   TBranch        *b_Jet_chEmEF;   //!
   TBranch        *b_Jet_chHEF;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_jercCHF;   //!
   TBranch        *b_Jet_jercCHPUF;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_muEF;   //!
   TBranch        *b_Jet_muonSubtrFactor;   //!
   TBranch        *b_Jet_neEmEF;   //!
   TBranch        *b_Jet_neHEF;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_qgl;   //!
   TBranch        *b_Jet_rawFactor;   //!
   TBranch        *b_Jet_bRegCorr;   //!
   TBranch        *b_Jet_bRegRes;   //!
   TBranch        *b_Jet_electronIdx1;   //!
   TBranch        *b_Jet_electronIdx2;   //!
   TBranch        *b_Jet_jetId;   //!
   TBranch        *b_Jet_muonIdx1;   //!
   TBranch        *b_Jet_muonIdx2;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_nElectrons;   //!
   TBranch        *b_Jet_nMuons;   //!
   TBranch        *b_Jet_puId;   //!
   TBranch        *b_L1PreFiringWeight_Dn;   //!
   TBranch        *b_L1PreFiringWeight_Nom;   //!
   TBranch        *b_L1PreFiringWeight_Up;   //!
   TBranch        *b_METFixEE2017_MetUnclustEnUpDeltaX;   //!
   TBranch        *b_METFixEE2017_MetUnclustEnUpDeltaY;   //!
   TBranch        *b_METFixEE2017_covXX;   //!
   TBranch        *b_METFixEE2017_covXY;   //!
   TBranch        *b_METFixEE2017_covYY;   //!
   TBranch        *b_METFixEE2017_phi;   //!
   TBranch        *b_METFixEE2017_pt;   //!
   TBranch        *b_METFixEE2017_significance;   //!
   TBranch        *b_METFixEE2017_sumEt;   //!
   TBranch        *b_GenMET_phi;   //!
   TBranch        *b_GenMET_pt;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaX;   //!
   TBranch        *b_MET_MetUnclustEnUpDeltaY;   //!
   TBranch        *b_MET_covXX;   //!
   TBranch        *b_MET_covXY;   //!
   TBranch        *b_MET_covYY;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_pt;   //!
   TBranch        *b_MET_significance;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_nMuon;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_dxyErr;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_dzErr;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_ip3d;   //!
   TBranch        *b_Muon_jetPtRelv2;   //!
   TBranch        *b_Muon_jetRelIso;   //!
   TBranch        *b_Muon_mass;   //!
   TBranch        *b_Muon_miniPFRelIso_all;   //!
   TBranch        *b_Muon_miniPFRelIso_chg;   //!
   TBranch        *b_Muon_pfRelIso03_all;   //!
   TBranch        *b_Muon_pfRelIso03_chg;   //!
   TBranch        *b_Muon_pfRelIso04_all;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptErr;   //!
   TBranch        *b_Muon_segmentComp;   //!
   TBranch        *b_Muon_sip3d;   //!
   TBranch        *b_Muon_tkRelIso;   //!
   TBranch        *b_Muon_tunepRelPt;   //!
   TBranch        *b_Muon_mvaLowPt;   //!
   TBranch        *b_Muon_mvaTTH;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_jetIdx;   //!
   TBranch        *b_Muon_nStations;   //!
   TBranch        *b_Muon_nTrackerLayers;   //!
   TBranch        *b_Muon_pdgId;   //!
   TBranch        *b_Muon_tightCharge;   //!
   TBranch        *b_Muon_highPtId;   //!
   TBranch        *b_Muon_inTimeMuon;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isPFcand;   //!
   TBranch        *b_Muon_isTracker;   //!
   TBranch        *b_Muon_looseId;   //!
   TBranch        *b_Muon_mediumId;   //!
   TBranch        *b_Muon_mediumPromptId;   //!
   TBranch        *b_Muon_miniIsoId;   //!
   TBranch        *b_Muon_multiIsoId;   //!
   TBranch        *b_Muon_mvaId;   //!
   TBranch        *b_Muon_pfIsoId;   //!
   TBranch        *b_Muon_softId;   //!
   TBranch        *b_Muon_softMvaId;   //!
   TBranch        *b_Muon_tightId;   //!
   TBranch        *b_Muon_tkIsoId;   //!
   TBranch        *b_Muon_triggerIdLoose;   //!
   TBranch        *b_nPhoton;   //!
   TBranch        *b_Photon_eCorr;   //!
   TBranch        *b_Photon_energyErr;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_hoe;   //!
   TBranch        *b_Photon_mass;   //!
   TBranch        *b_Photon_mvaID;   //!
   TBranch        *b_Photon_mvaIDV1;   //!
   TBranch        *b_Photon_pfRelIso03_all;   //!
   TBranch        *b_Photon_pfRelIso03_chg;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_sieie;   //!
   TBranch        *b_Photon_charge;   //!
   TBranch        *b_Photon_cutBasedBitmap;   //!
   TBranch        *b_Photon_cutBasedV1Bitmap;   //!
   TBranch        *b_Photon_electronIdx;   //!
   TBranch        *b_Photon_jetIdx;   //!
   TBranch        *b_Photon_pdgId;   //!
   TBranch        *b_Photon_vidNestedWPBitmap;   //!
   TBranch        *b_Photon_electronVeto;   //!
   TBranch        *b_Photon_isScEtaEB;   //!
   TBranch        *b_Photon_isScEtaEE;   //!
   TBranch        *b_Photon_mvaID_WP80;   //!
   TBranch        *b_Photon_mvaID_WP90;   //!
   TBranch        *b_Photon_pixelSeed;   //!
   TBranch        *b_Photon_seedGain;   //!
   TBranch        *b_Photon_dEsigmaDown;
   TBranch        *b_Photon_dEsigmaUp;
   TBranch        *b_Pileup_nTrueInt;   //!
   TBranch        *b_Pileup_pudensity;   //!
   TBranch        *b_Pileup_gpudensity;   //!
   TBranch        *b_Pileup_nPU;   //!
   TBranch        *b_Pileup_sumEOOT;   //!
   TBranch        *b_Pileup_sumLOOT;   //!
   TBranch        *b_PuppiMET_phi;   //!
   TBranch        *b_PuppiMET_pt;   //!
   TBranch        *b_PuppiMET_sumEt;   //!
   TBranch        *b_RawMET_phi;   //!
   TBranch        *b_RawMET_pt;   //!
   TBranch        *b_RawMET_sumEt;   //!
   TBranch        *b_fixedGridRhoFastjetAll;   //!
   TBranch        *b_fixedGridRhoFastjetCentral;   //!
   TBranch        *b_fixedGridRhoFastjetCentralCalo;   //!
   TBranch        *b_fixedGridRhoFastjetCentralChargedPileUp;   //!
   TBranch        *b_fixedGridRhoFastjetCentralNeutral;   //!
   TBranch        *b_nGenDressedLepton;   //!
   TBranch        *b_GenDressedLepton_eta;   //!
   TBranch        *b_GenDressedLepton_mass;   //!
   TBranch        *b_GenDressedLepton_phi;   //!
   TBranch        *b_GenDressedLepton_pt;   //!
   TBranch        *b_GenDressedLepton_pdgId;   //!
   TBranch        *b_GenDressedLepton_hasTauAnc;   //!
   TBranch        *b_nSoftActivityJet;   //!
   TBranch        *b_SoftActivityJet_eta;   //!
   TBranch        *b_SoftActivityJet_phi;   //!
   TBranch        *b_SoftActivityJet_pt;   //!
   TBranch        *b_SoftActivityJetHT;   //!
   TBranch        *b_SoftActivityJetHT10;   //!
   TBranch        *b_SoftActivityJetHT2;   //!
   TBranch        *b_SoftActivityJetHT5;   //!
   TBranch        *b_SoftActivityJetNjets10;   //!
   TBranch        *b_SoftActivityJetNjets2;   //!
   TBranch        *b_SoftActivityJetNjets5;   //!
   TBranch        *b_nSubJet;   //!
   TBranch        *b_SubJet_btagCMVA;   //!
   TBranch        *b_SubJet_btagCSVV2;   //!
   TBranch        *b_SubJet_btagDeepB;   //!
   TBranch        *b_SubJet_eta;   //!
   TBranch        *b_SubJet_mass;   //!
   TBranch        *b_SubJet_n2b1;   //!
   TBranch        *b_SubJet_n3b1;   //!
   TBranch        *b_SubJet_phi;   //!
   TBranch        *b_SubJet_pt;   //!
   TBranch        *b_SubJet_rawFactor;   //!
   TBranch        *b_SubJet_tau1;   //!
   TBranch        *b_SubJet_tau2;   //!
   TBranch        *b_SubJet_tau3;   //!
   TBranch        *b_SubJet_tau4;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_Tau_chargedIso;   //!
   TBranch        *b_Tau_dxy;   //!
   TBranch        *b_Tau_dz;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_leadTkDeltaEta;   //!
   TBranch        *b_Tau_leadTkDeltaPhi;   //!
   TBranch        *b_Tau_leadTkPtOverTauPt;   //!
   TBranch        *b_Tau_mass;   //!
   TBranch        *b_Tau_neutralIso;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_photonsOutsideSignalCone;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_puCorr;   //!
   TBranch        *b_Tau_rawAntiEle;   //!
   TBranch        *b_Tau_rawAntiEle2018;   //!
   TBranch        *b_Tau_rawIso;   //!
   TBranch        *b_Tau_rawIsodR03;   //!
   TBranch        *b_Tau_rawMVAnewDM2017v2;   //!
   TBranch        *b_Tau_rawMVAoldDM;   //!
   TBranch        *b_Tau_rawMVAoldDM2017v1;   //!
   TBranch        *b_Tau_rawMVAoldDM2017v2;   //!
   TBranch        *b_Tau_rawMVAoldDMdR032017v2;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_decayMode;   //!
   TBranch        *b_Tau_jetIdx;   //!
   TBranch        *b_Tau_rawAntiEleCat;   //!
   TBranch        *b_Tau_rawAntiEleCat2018;   //!
   TBranch        *b_Tau_idAntiEle;   //!
   TBranch        *b_Tau_idAntiEle2018;   //!
   TBranch        *b_Tau_idAntiMu;   //!
   TBranch        *b_Tau_idDecayMode;   //!
   TBranch        *b_Tau_idDecayModeNewDMs;   //!
   TBranch        *b_Tau_idMVAnewDM2017v2;   //!
   TBranch        *b_Tau_idMVAoldDM;   //!
   TBranch        *b_Tau_idMVAoldDM2017v1;   //!
   TBranch        *b_Tau_idMVAoldDM2017v2;   //!
   TBranch        *b_Tau_idMVAoldDMdR032017v2;   //!
   TBranch        *b_TkMET_phi;   //!
   TBranch        *b_TkMET_pt;   //!
   TBranch        *b_TkMET_sumEt;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_TrigObj_pt;   //!
   TBranch        *b_TrigObj_eta;   //!
   TBranch        *b_TrigObj_phi;   //!
   TBranch        *b_TrigObj_l1pt;   //!
   TBranch        *b_TrigObj_l1pt_2;   //!
   TBranch        *b_TrigObj_l2pt;   //!
   TBranch        *b_TrigObj_id;   //!
   TBranch        *b_TrigObj_l1iso;   //!
   TBranch        *b_TrigObj_l1charge;   //!
   TBranch        *b_TrigObj_filterBits;   //!
   TBranch        *b_genTtbarId;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_OtherPV_z;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_SV_dlen;   //!
   TBranch        *b_SV_dlenSig;   //!
   TBranch        *b_SV_pAngle;   //!
   TBranch        *b_Electron_genPartIdx;   //!
   TBranch        *b_Electron_genPartFlav;   //!
   TBranch        *b_GenJetAK8_partonFlavour;   //!
   TBranch        *b_GenJetAK8_hadronFlavour;   //!
   TBranch        *b_GenJet_partonFlavour;   //!
   TBranch        *b_GenJet_hadronFlavour;   //!
   TBranch        *b_Jet_genJetIdx;   //!
   TBranch        *b_Jet_hadronFlavour;   //!
   TBranch        *b_Jet_partonFlavour;   //!
   TBranch        *b_Muon_genPartIdx;   //!
   TBranch        *b_Muon_genPartFlav;   //!
   TBranch        *b_Photon_genPartIdx;   //!
   TBranch        *b_Photon_genPartFlav;   //!
   TBranch        *b_MET_fiducialGenPhi;   //!
   TBranch        *b_MET_fiducialGenPt;   //!
   TBranch        *b_Electron_cleanmask;   //!
   TBranch        *b_Jet_cleanmask;   //!
   TBranch        *b_Muon_cleanmask;   //!
   TBranch        *b_Photon_cleanmask;   //!
   TBranch        *b_Tau_cleanmask;   //!
   TBranch        *b_SV_chi2;   //!
   TBranch        *b_SV_eta;   //!
   TBranch        *b_SV_mass;   //!
   TBranch        *b_SV_ndof;   //!
   TBranch        *b_SV_phi;   //!
   TBranch        *b_SV_pt;   //!
   TBranch        *b_SV_x;   //!
   TBranch        *b_SV_y;   //!
   TBranch        *b_SV_z;   //!
   TBranch        *b_Tau_genPartIdx;   //!
   TBranch        *b_Tau_genPartFlav;   //!
   TBranch        *b_L1simulation_step;   //!

   //TBranch        *b_TrigObj;
   TBranch        *b_HLTriggerFirstPath;   //!
   TBranch        *b_HLT_DoubleMu20_7_Mass0to30_Photon23;  //!
   TBranch        *b_HLT_SingleMu24;
   TBranch        *b_HLT_SingleMu27;
   TBranch        *b_HLT_DoubleMu4_3_Jpsi;
   TBranch        *b_HLT_Dimuon20_Jpsi;
   TBranch        *b_HLT_Dimuon25_Jpsi; //!
   TBranch        *b_HLT_Mu17_Photon30_IsoCaloId;   //!
   TBranch        *b_HLT_Mu17_Photon30_CaloIdL_L1ISO;
   TBranch        *b_HLTriggerFinalPath;   //!

   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_globalTightHalo2016Filter;   //!
   TBranch        *b_Flag_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_ecalBadCalibFilter;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_Flag_muonBadTrackFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateFilter;   //!
   TBranch        *b_Flag_BadPFMuonFilter;   //!
   TBranch        *b_Flag_BadChargedCandidateSummer16Filter;   //!
   TBranch        *b_Flag_BadPFMuonSummer16Filter;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_METFilters;   //!
   TBranch        *b_L1_AlwaysTrue;   //!
   TBranch        *b_L1_BPTX_AND_Ref1_VME;   //!
   TBranch        *b_L1_BPTX_AND_Ref3_VME;   //!
   TBranch        *b_L1_BPTX_AND_Ref4_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_B1_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_B2_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_Ref1_VME;   //!
   TBranch        *b_L1_BPTX_BeamGas_Ref2_VME;   //!
   TBranch        *b_L1_BPTX_NotOR_VME;   //!
   TBranch        *b_L1_BPTX_OR_Ref3_VME;   //!
   TBranch        *b_L1_BPTX_OR_Ref4_VME;   //!
   TBranch        *b_L1_BPTX_RefAND_VME;   //!
   TBranch        *b_L1_BptxMinus;   //!
   TBranch        *b_L1_BptxOR;   //!
   TBranch        *b_L1_BptxPlus;   //!
   TBranch        *b_L1_BptxXOR;   //!
   TBranch        *b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142;   //!
   TBranch        *b_L1_DoubleEG6_HTT240er;   //!
   TBranch        *b_L1_DoubleEG6_HTT250er;   //!
   TBranch        *b_L1_DoubleEG6_HTT255er;   //!
   TBranch        *b_L1_DoubleEG6_HTT270er;   //!
   TBranch        *b_L1_DoubleEG6_HTT300er;   //!
   TBranch        *b_L1_DoubleEG8er2p6_HTT255er;   //!
   TBranch        *b_L1_DoubleEG8er2p6_HTT270er;   //!
   TBranch        *b_L1_DoubleEG8er2p6_HTT300er;   //!
   TBranch        *b_L1_DoubleEG_15_10;   //!
   TBranch        *b_L1_DoubleEG_18_17;   //!
   TBranch        *b_L1_DoubleEG_20_18;   //!
   TBranch        *b_L1_DoubleEG_22_10;   //!
   TBranch        *b_L1_DoubleEG_22_12;   //!
   TBranch        *b_L1_DoubleEG_22_15;   //!
   TBranch        *b_L1_DoubleEG_23_10;   //!
   TBranch        *b_L1_DoubleEG_24_17;   //!
   TBranch        *b_L1_DoubleEG_25_12;   //!
   TBranch        *b_L1_DoubleEG_25_13;   //!
   TBranch        *b_L1_DoubleEG_25_14;   //!
   TBranch        *b_L1_DoubleEG_LooseIso23_10;   //!
   TBranch        *b_L1_DoubleEG_LooseIso24_10;   //!
   TBranch        *b_L1_DoubleIsoTau28er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau30er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau32er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau33er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau34er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau35er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau36er2p1;   //!
   TBranch        *b_L1_DoubleIsoTau38er2p1;   //!
   TBranch        *b_L1_DoubleJet100er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleJet100er2p7;   //!
   TBranch        *b_L1_DoubleJet112er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_DoubleJet112er2p7;   //!
   TBranch        *b_L1_DoubleJet120er2p7;   //!
   TBranch        *b_L1_DoubleJet150er2p7;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min300_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min320_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min340_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min360_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min380_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min400_Mu10;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min400_Mu6;   //!
   TBranch        *b_L1_DoubleJet30_Mass_Min400_dEta_Max1p5;   //!
   TBranch        *b_L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450;   //!
   TBranch        *b_L1_DoubleJet40er2p7;   //!
   TBranch        *b_L1_DoubleJet50er2p7;   //!
   TBranch        *b_L1_DoubleJet60er2p7;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM100;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM60;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM70;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM80;   //!
   TBranch        *b_L1_DoubleJet60er2p7_ETM90;   //!
   TBranch        *b_L1_DoubleJet80er2p7;   //!
   TBranch        *b_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_100_35_DoubleJet35_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_110_40_DoubleJet40_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_115_35_DoubleJet35_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620;   //!
   TBranch        *b_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620;   //!
   TBranch        *b_L1_DoubleLooseIsoEG22er2p1;   //!
   TBranch        *b_L1_DoubleLooseIsoEG24er2p1;   //!
   TBranch        *b_L1_DoubleMu0;   //!
   TBranch        *b_L1_DoubleMu0_ETM40;   //!
   TBranch        *b_L1_DoubleMu0_ETM55;   //!
   TBranch        *b_L1_DoubleMu0_ETM60;   //!
   TBranch        *b_L1_DoubleMu0_ETM65;   //!
   TBranch        *b_L1_DoubleMu0_ETM70;   //!
   TBranch        *b_L1_DoubleMu0_SQ;   //!
   TBranch        *b_L1_DoubleMu0_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er1p4_dEta_Max1p8_OS;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er1p5_SQ_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu0er2_SQ_dR_Max1p4;   //!
   TBranch        *b_L1_DoubleMu18er2p1;   //!
   TBranch        *b_L1_DoubleMu22er2p1;   //!
   TBranch        *b_L1_DoubleMu3_OS_DoubleEG7p5Upsilon;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT100er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT200er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT220er;   //!
   TBranch        *b_L1_DoubleMu3_SQ_HTT240er;   //!
   TBranch        *b_L1_DoubleMu4_OS_EG12;   //!
   TBranch        *b_L1_DoubleMu4_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4_SQ_OS_dR_Max1p2;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18;   //!
   TBranch        *b_L1_DoubleMu5Upsilon_OS_DoubleEG3;   //!
   TBranch        *b_L1_DoubleMu5_OS_EG12;   //!
   TBranch        *b_L1_DoubleMu5_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu5_SQ_OS_Mass7to18;   //!
   TBranch        *b_L1_DoubleMu6_SQ_OS;   //!
   TBranch        *b_L1_DoubleMu7_EG7;   //!
   TBranch        *b_L1_DoubleMu7_SQ_EG7;   //!
   TBranch        *b_L1_DoubleMu8_SQ;   //!
   TBranch        *b_L1_DoubleMu_10_0_dEta_Max1p8;   //!
   TBranch        *b_L1_DoubleMu_11_4;   //!
   TBranch        *b_L1_DoubleMu_12_5;   //!
   TBranch        *b_L1_DoubleMu_12_8;   //!
   TBranch        *b_L1_DoubleMu_13_6;   //!
   TBranch        *b_L1_DoubleMu_15_5;   //!
   TBranch        *b_L1_DoubleMu_15_5_SQ;   //!
   TBranch        *b_L1_DoubleMu_15_7;   //!
   TBranch        *b_L1_DoubleMu_15_7_SQ;   //!
   TBranch        *b_L1_DoubleMu_15_7_SQ_Mass_Min4;   //!
   TBranch        *b_L1_DoubleMu_20_2_SQ_Mass_Max20;   //!
   TBranch        *b_L1_DoubleTau50er2p1;   //!
   TBranch        *b_L1_DoubleTau70er2p1;   //!
   TBranch        *b_L1_EG25er2p1_HTT125er;   //!
   TBranch        *b_L1_EG27er2p1_HTT200er;   //!
   TBranch        *b_L1_ETM100;   //!
   TBranch        *b_L1_ETM100_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM105;   //!
   TBranch        *b_L1_ETM110;   //!
   TBranch        *b_L1_ETM110_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM115;   //!
   TBranch        *b_L1_ETM120;   //!
   TBranch        *b_L1_ETM150;   //!
   TBranch        *b_L1_ETM30;   //!
   TBranch        *b_L1_ETM40;   //!
   TBranch        *b_L1_ETM50;   //!
   TBranch        *b_L1_ETM60;   //!
   TBranch        *b_L1_ETM70;   //!
   TBranch        *b_L1_ETM75;   //!
   TBranch        *b_L1_ETM75_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM80;   //!
   TBranch        *b_L1_ETM80_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM85;   //!
   TBranch        *b_L1_ETM90;   //!
   TBranch        *b_L1_ETM90_Jet60_dPhi_Min0p4;   //!
   TBranch        *b_L1_ETM95;   //!
   TBranch        *b_L1_ETMHF100;   //!
   TBranch        *b_L1_ETMHF100_HTT60er;   //!
   TBranch        *b_L1_ETMHF100_Jet60_OR_DiJet30woTT28;   //!
   TBranch        *b_L1_ETMHF100_Jet60_OR_DoubleJet30;   //!
   TBranch        *b_L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF110;   //!
   TBranch        *b_L1_ETMHF110_HTT60er;   //!
   TBranch        *b_L1_ETMHF110_Jet60_OR_DiJet30woTT28;   //!
   TBranch        *b_L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF120;   //!
   TBranch        *b_L1_ETMHF120_HTT60er;   //!
   TBranch        *b_L1_ETMHF120_Jet60_OR_DiJet30woTT28;   //!
   TBranch        *b_L1_ETMHF150;   //!
   TBranch        *b_L1_ETMHF70;   //!
   TBranch        *b_L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF80;   //!
   TBranch        *b_L1_ETMHF80_HTT60er;   //!
   TBranch        *b_L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETMHF90;   //!
   TBranch        *b_L1_ETMHF90_HTT60er;   //!
   TBranch        *b_L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30;   //!
   TBranch        *b_L1_ETT100_BptxAND;   //!
   TBranch        *b_L1_ETT110_BptxAND;   //!
   TBranch        *b_L1_ETT40_BptxAND;   //!
   TBranch        *b_L1_ETT50_BptxAND;   //!
   TBranch        *b_L1_ETT60_BptxAND;   //!
   TBranch        *b_L1_ETT70_BptxAND;   //!
   TBranch        *b_L1_ETT75_BptxAND;   //!
   TBranch        *b_L1_ETT80_BptxAND;   //!
   TBranch        *b_L1_ETT85_BptxAND;   //!
   TBranch        *b_L1_ETT90_BptxAND;   //!
   TBranch        *b_L1_ETT95_BptxAND;   //!
   TBranch        *b_L1_FirstBunchAfterTrain;   //!
   TBranch        *b_L1_FirstBunchInTrain;   //!
   TBranch        *b_L1_FirstCollisionInOrbit;   //!
   TBranch        *b_L1_FirstCollisionInTrain;   //!
   TBranch        *b_L1_HTT120er;   //!
   TBranch        *b_L1_HTT160er;   //!
   TBranch        *b_L1_HTT200er;   //!
   TBranch        *b_L1_HTT220er;   //!
   TBranch        *b_L1_HTT240er;   //!
   TBranch        *b_L1_HTT250er_QuadJet_70_55_40_35_er2p5;   //!
   TBranch        *b_L1_HTT255er;   //!
   TBranch        *b_L1_HTT270er;   //!
   TBranch        *b_L1_HTT280er;   //!
   TBranch        *b_L1_HTT280er_QuadJet_70_55_40_35_er2p5;   //!
   TBranch        *b_L1_HTT300er;   //!
   TBranch        *b_L1_HTT300er_QuadJet_70_55_40_35_er2p5;   //!
   TBranch        *b_L1_HTT320er;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_40_40_er2p4;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_40_40_er2p5;   //!
   TBranch        *b_L1_HTT320er_QuadJet_70_55_45_45_er2p5;   //!
   TBranch        *b_L1_HTT340er;   //!
   TBranch        *b_L1_HTT340er_QuadJet_70_55_40_40_er2p5;   //!
   TBranch        *b_L1_HTT340er_QuadJet_70_55_45_45_er2p5;   //!
   TBranch        *b_L1_HTT380er;   //!
   TBranch        *b_L1_HTT400er;   //!
   TBranch        *b_L1_HTT450er;   //!
   TBranch        *b_L1_HTT500er;   //!
   TBranch        *b_L1_IsoEG33_Mt40;   //!
   TBranch        *b_L1_IsoEG33_Mt44;   //!
   TBranch        *b_L1_IsoEG33_Mt48;   //!
   TBranch        *b_L1_IsoTau40er_ETM100;   //!
   TBranch        *b_L1_IsoTau40er_ETM105;   //!
   TBranch        *b_L1_IsoTau40er_ETM110;   //!
   TBranch        *b_L1_IsoTau40er_ETM115;   //!
   TBranch        *b_L1_IsoTau40er_ETM120;   //!
   TBranch        *b_L1_IsoTau40er_ETM80;   //!
   TBranch        *b_L1_IsoTau40er_ETM85;   //!
   TBranch        *b_L1_IsoTau40er_ETM90;   //!
   TBranch        *b_L1_IsoTau40er_ETM95;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF100;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF110;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF120;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF80;   //!
   TBranch        *b_L1_IsoTau40er_ETMHF90;   //!
   TBranch        *b_L1_IsolatedBunch;   //!
   TBranch        *b_L1_LastCollisionInTrain;   //!
   TBranch        *b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7;   //!
   TBranch        *b_L1_LooseIsoEG26er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG28er2p1_HTT100er;   //!
   TBranch        *b_L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3;   //!
   TBranch        *b_L1_MU20_EG15;   //!
   TBranch        *b_L1_MinimumBiasHF0_AND_BptxAND;   //!
   TBranch        *b_L1_MinimumBiasHF0_OR_BptxAND;   //!
   TBranch        *b_L1_Mu10er2p1_ETM30;   //!
   TBranch        *b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_Mu12_EG10;   //!
   TBranch        *b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6;   //!
   TBranch        *b_L1_Mu14er2p1_ETM30;   //!
   TBranch        *b_L1_Mu15_HTT100er;   //!
   TBranch        *b_L1_Mu18_HTT100er;   //!
   TBranch        *b_L1_Mu18_Jet24er2p7;   //!
   TBranch        *b_L1_Mu18er2p1_IsoTau26er2p1;   //!
   TBranch        *b_L1_Mu18er2p1_Tau24er2p1;   //!
   TBranch        *b_L1_Mu20_EG10;   //!
   TBranch        *b_L1_Mu20_EG17;   //!
   TBranch        *b_L1_Mu20_LooseIsoEG6;   //!
   TBranch        *b_L1_Mu20er2p1_IsoTau26er2p1;   //!
   TBranch        *b_L1_Mu20er2p1_IsoTau27er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau28er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau30er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau32er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau33er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau34er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau35er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau36er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau38er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_IsoTau40er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_Tau50er2p1;   //!
   TBranch        *b_L1_Mu22er2p1_Tau70er2p1;   //!
   TBranch        *b_L1_Mu23_EG10;   //!
   TBranch        *b_L1_Mu23_LooseIsoEG10;   //!
   TBranch        *b_L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4;   //!
   TBranch        *b_L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4;   //!
   TBranch        *b_L1_Mu3_Jet30er2p5;   //!
   TBranch        *b_L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4;   //!
   TBranch        *b_L1_Mu5_EG15;   //!
   TBranch        *b_L1_Mu5_EG20;   //!
   TBranch        *b_L1_Mu5_EG23;   //!
   TBranch        *b_L1_Mu5_LooseIsoEG18;   //!
   TBranch        *b_L1_Mu5_LooseIsoEG20;   //!
   TBranch        *b_L1_Mu6_DoubleEG10;   //!
   TBranch        *b_L1_Mu6_DoubleEG17;   //!
   TBranch        *b_L1_Mu6_HTT200er;   //!
   TBranch        *b_L1_Mu6_HTT240er;   //!
   TBranch        *b_L1_Mu6_HTT250er;   //!
   TBranch        *b_L1_Mu7_EG23;   //!
   TBranch        *b_L1_Mu7_LooseIsoEG20;   //!
   TBranch        *b_L1_Mu7_LooseIsoEG23;   //!
   TBranch        *b_L1_Mu8_HTT150er;   //!
   TBranch        *b_L1_NotBptxOR;   //!
   TBranch        *b_L1_QuadJet36er2p7_IsoTau52er2p1;   //!
   TBranch        *b_L1_QuadJet36er2p7_Tau52;   //!
   TBranch        *b_L1_QuadJet40er2p7;   //!
   TBranch        *b_L1_QuadJet50er2p7;   //!
   TBranch        *b_L1_QuadJet60er2p7;   //!
   TBranch        *b_L1_QuadMu0;   //!
   TBranch        *b_L1_SingleEG10;   //!
   TBranch        *b_L1_SingleEG15;   //!
   TBranch        *b_L1_SingleEG18;   //!
   TBranch        *b_L1_SingleEG24;   //!
   TBranch        *b_L1_SingleEG26;   //!
   TBranch        *b_L1_SingleEG28;   //!
   TBranch        *b_L1_SingleEG2_BptxAND;   //!
   TBranch        *b_L1_SingleEG30;   //!
   TBranch        *b_L1_SingleEG32;   //!
   TBranch        *b_L1_SingleEG34;   //!
   TBranch        *b_L1_SingleEG34er2p1;   //!
   TBranch        *b_L1_SingleEG36;   //!
   TBranch        *b_L1_SingleEG36er2p1;   //!
   TBranch        *b_L1_SingleEG38;   //!
   TBranch        *b_L1_SingleEG38er2p1;   //!
   TBranch        *b_L1_SingleEG40;   //!
   TBranch        *b_L1_SingleEG42;   //!
   TBranch        *b_L1_SingleEG45;   //!
   TBranch        *b_L1_SingleEG5;   //!
   TBranch        *b_L1_SingleEG50;   //!
   TBranch        *b_L1_SingleIsoEG18;   //!
   TBranch        *b_L1_SingleIsoEG18er2p1;   //!
   TBranch        *b_L1_SingleIsoEG20;   //!
   TBranch        *b_L1_SingleIsoEG20er2p1;   //!
   TBranch        *b_L1_SingleIsoEG22;   //!
   TBranch        *b_L1_SingleIsoEG22er2p1;   //!
   TBranch        *b_L1_SingleIsoEG24;   //!
   TBranch        *b_L1_SingleIsoEG24er2p1;   //!
   TBranch        *b_L1_SingleIsoEG26;   //!
   TBranch        *b_L1_SingleIsoEG26er2p1;   //!
   TBranch        *b_L1_SingleIsoEG28;   //!
   TBranch        *b_L1_SingleIsoEG28er2p1;   //!
   TBranch        *b_L1_SingleIsoEG30;   //!
   TBranch        *b_L1_SingleIsoEG30er2p1;   //!
   TBranch        *b_L1_SingleIsoEG32;   //!
   TBranch        *b_L1_SingleIsoEG32er2p1;   //!
   TBranch        *b_L1_SingleIsoEG33er2p1;   //!
   TBranch        *b_L1_SingleIsoEG34;   //!
   TBranch        *b_L1_SingleIsoEG34er2p1;   //!
   TBranch        *b_L1_SingleIsoEG35;   //!
   TBranch        *b_L1_SingleIsoEG35er2p1;   //!
   TBranch        *b_L1_SingleIsoEG36;   //!
   TBranch        *b_L1_SingleIsoEG36er2p1;   //!
   TBranch        *b_L1_SingleIsoEG37;   //!
   TBranch        *b_L1_SingleIsoEG38;   //!
   TBranch        *b_L1_SingleIsoEG38er2p1;   //!
   TBranch        *b_L1_SingleIsoEG40;   //!
   TBranch        *b_L1_SingleIsoEG40er2p1;   //!
   TBranch        *b_L1_SingleJet120;   //!
   TBranch        *b_L1_SingleJet120_FWD;   //!
   TBranch        *b_L1_SingleJet12_BptxAND;   //!
   TBranch        *b_L1_SingleJet140;   //!
   TBranch        *b_L1_SingleJet150;   //!
   TBranch        *b_L1_SingleJet16;   //!
   TBranch        *b_L1_SingleJet160;   //!
   TBranch        *b_L1_SingleJet170;   //!
   TBranch        *b_L1_SingleJet180;   //!
   TBranch        *b_L1_SingleJet20;   //!
   TBranch        *b_L1_SingleJet200;   //!
   TBranch        *b_L1_SingleJet20er2p7_NotBptxOR;   //!
   TBranch        *b_L1_SingleJet20er2p7_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet35;   //!
   TBranch        *b_L1_SingleJet35_FWD;   //!
   TBranch        *b_L1_SingleJet35_HFm;   //!
   TBranch        *b_L1_SingleJet35_HFp;   //!
   TBranch        *b_L1_SingleJet43er2p7_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet46er2p7_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleJet60;   //!
   TBranch        *b_L1_SingleJet60_FWD;   //!
   TBranch        *b_L1_SingleJet60_HFm;   //!
   TBranch        *b_L1_SingleJet60_HFp;   //!
   TBranch        *b_L1_SingleJet90;   //!
   TBranch        *b_L1_SingleJet90_FWD;   //!
   TBranch        *b_L1_SingleMu0_BMTF;   //!
   TBranch        *b_L1_SingleMu0_EMTF;   //!
   TBranch        *b_L1_SingleMu0_OMTF;   //!
   TBranch        *b_L1_SingleMu10_LowQ;   //!
   TBranch        *b_L1_SingleMu11_LowQ;   //!
   TBranch        *b_L1_SingleMu12_LowQ_BMTF;   //!
   TBranch        *b_L1_SingleMu12_LowQ_EMTF;   //!
   TBranch        *b_L1_SingleMu12_LowQ_OMTF;   //!
   TBranch        *b_L1_SingleMu14er2p1;   //!
   TBranch        *b_L1_SingleMu16;   //!
   TBranch        *b_L1_SingleMu16er2p1;   //!
   TBranch        *b_L1_SingleMu18;   //!
   TBranch        *b_L1_SingleMu18er2p1;   //!
   TBranch        *b_L1_SingleMu20;   //!
   TBranch        *b_L1_SingleMu20er2p1;   //!
   TBranch        *b_L1_SingleMu22;   //!
   TBranch        *b_L1_SingleMu22_BMTF;   //!
   TBranch        *b_L1_SingleMu22_EMTF;   //!
   TBranch        *b_L1_SingleMu22_OMTF;   //!
   TBranch        *b_L1_SingleMu22er2p1;   //!
   TBranch        *b_L1_SingleMu25;   //!
   TBranch        *b_L1_SingleMu3;   //!
   TBranch        *b_L1_SingleMu30;   //!
   TBranch        *b_L1_SingleMu5;   //!
   TBranch        *b_L1_SingleMu7;   //!
   TBranch        *b_L1_SingleMuCosmics;   //!
   TBranch        *b_L1_SingleMuCosmics_BMTF;   //!
   TBranch        *b_L1_SingleMuCosmics_EMTF;   //!
   TBranch        *b_L1_SingleMuCosmics_OMTF;   //!
   TBranch        *b_L1_SingleMuOpen;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR;   //!
   TBranch        *b_L1_SingleMuOpen_NotBptxOR_3BX;   //!
   TBranch        *b_L1_SingleTau100er2p1;   //!
   TBranch        *b_L1_SingleTau120er2p1;   //!
   TBranch        *b_L1_SingleTau130er2p1;   //!
   TBranch        *b_L1_SingleTau140er2p1;   //!
   TBranch        *b_L1_SingleTau20;   //!
   TBranch        *b_L1_SingleTau80er2p1;   //!
   TBranch        *b_L1_TripleEG_14_10_8;   //!
   TBranch        *b_L1_TripleEG_18_17_8;   //!
   TBranch        *b_L1_TripleEG_LooseIso20_10_5;   //!
   TBranch        *b_L1_TripleJet_100_85_72_VBF;   //!
   TBranch        *b_L1_TripleJet_105_85_76_VBF;   //!
   TBranch        *b_L1_TripleJet_84_68_48_VBF;   //!
   TBranch        *b_L1_TripleJet_88_72_56_VBF;   //!
   TBranch        *b_L1_TripleJet_92_76_64_VBF;   //!
   TBranch        *b_L1_TripleJet_98_83_71_VBF;   //!
   TBranch        *b_L1_TripleMu0;   //!
   TBranch        *b_L1_TripleMu0_OQ;   //!
   TBranch        *b_L1_TripleMu3;   //!
   TBranch        *b_L1_TripleMu3_SQ;   //!
   TBranch        *b_L1_TripleMu_4_4_4;   //!
   TBranch        *b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0OQ;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9;   //!
   TBranch        *b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9;   //!
   TBranch        *b_L1_TripleMu_5_0_0;   //!
   TBranch        *b_L1_TripleMu_5_3_3;   //!
   TBranch        *b_L1_TripleMu_5_3p5_2p5;   //!
   TBranch        *b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17;   //!
   TBranch        *b_L1_TripleMu_5_5_3;   //!
   TBranch        *b_L1_UnpairedBunchBptxMinus;   //!
   TBranch        *b_L1_UnpairedBunchBptxPlus;   //!
   TBranch        *b_L1_ZeroBias;   //!
   TBranch        *b_L1_ZeroBias_copy;   //!


   ZtoJPsiGamma(TTree *tree = 0);
   virtual ~ZtoJPsiGamma();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   TChain* chain = new TChain();
   TChain* chain1 = new TChain();
   // *******************************************************
   /*
   namespace AnaUtil {
     // Templated functioned must be defined in the header itself                                                                            
     void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters=" ");
     double cutValue(const std::map<std::string, double>& m, const std::string& cname);
     const std::map<std::string, double>& cutMap(const std::map<std::string, std::map<std::string, double>>& hmap, const std::string& mkey);
     void buildList(const std::vector<std::string>& tokens, std::vector<std::string>& list);
     void buildMap(const std::vector<std::string>& tokens, std::map<std::string, int>& hmap);
     void buildMap(const std::vector<std::string>& tokens, std::unordered_map<std::string, int>& hmap);
     void storeCuts(const std::vector<std::string>& tokens, std::map<std::string, std::map<std::string, double>>& hmap);
     void showCuts(const std::map<std::string, std::map<std::string, double> >& hmap, ostream& os=std::cout);
     
     // print_list_elements()                                                                                                           
      // - prints optional C-string optcstr followed by                    * - all elements of the collection coll                           
 
      
     template <class T>
       void showList(const T& coll, const char* optcstr="", std::ostream& os=std::cout) {
       os << optcstr << ", Total # = " << coll.size() << ":" << std::endl;
       for (auto const& v: coll)
	 os << v << std::endl;
     }
     template <class T1, class T2>
       void showMap(const std::map<T1,T2>& m, const char* optcstr="", std::ostream& os=std::cout) {
       os << optcstr << std::endl;
       for (auto const& k: m)
	 os << k.first << std::endl;
     }
     template <class T1, class T2>
       void showMap(const std::unordered_map<T1,T2>& m, const char* optcstr="", std::ostream& os=std::cout) {
       os << optcstr << std::endl;
       for (auto const& k: m)
	 os << k.first << std::endl;
     }
     //copyList 
     template <class T>
       void copyList (const T& sourceColl, T& destColl) {
       destColl.clear();
       for (auto const& v: sourceColl)
	 destColl.push_back(v);
     }
   }


   virtual bool readJob(const string& jobFile);
*/

   //**************************************

   int setInputFile(const std::string& fname);
   virtual bool readJob(const std::string& jobFile, int& nFiles);
   //const std::unique_ptr<TFile>& histf() const {return histf_;}
   //std::unique_ptr<TFile>& histf() {return histf_;}
   bool isMC() const {return isMC_;}
   bool isSignal() const {return isSignal_;}
   int hlt() const {return hlt_;}
   int year() const{return year_;}
   bool isPre() const {return isPre_;}
   
   //double xsection() const {return xsec_;}

   //double cutValue(const std::map<std::string, double>& m, const std::string& cname);
   //std::map<std::string, double>& cutMap(const std::map<std::string, std::map<std::string, double>>& hmap, const std::string& mkey);
   //std::map<std::string, double>& lumiWtMap()  {return cutMap(hmap_, "lumiWtList");}
   //std::unique_ptr<TFile> histf_;       // The output file with histograms;
   bool isMC_ {false};
   bool isSignal_ {false};
   std::vector<std::string> fileList_;
   std::string histFile_ ;//{"default.root"};
   float den;
   int hlt_{-1};
   int year_{-1};
   float intLumi_{-1.};
   float xsec_{-1.};
   int nevents_{-1};
   bool useLumiWt_{false};
   bool isPre_{false};
   int maxEvt_;
   int maxEvent_;
   bool ispu_wt{false};
   bool isgen_wt{false};
   bool ispu_gen_wt{false};
   bool istlumi_wt{false};
   bool istlumi_wt_hasGen{false};
   bool find_SU{false};
   bool createMVATree_ {false};
   bool readMVA_ {false};
   std::string mvaSkimFile_ {""};
   std::string muon_SUfile;
   std::string mvAlgo_ {"BDTG"};
   std::string mvaXMLFile_ {""};
   std::unique_ptr<MVASkim> skimObj_ {nullptr};
   std::unique_ptr<MVAnalysis> mvAna_ {nullptr};

   std::map<std::string, std::map<std::string, double>> hmap_;
   double cutValue(const std::map<std::string, double>& m, const std::string& cname);                                           
   void storeCuts(const std::vector<std::string>& tokens, std::map<std::string, std::map<std::string, double>>& hmap);     

   const std::map<std::string, double>& cutMap(const std::map<std::string, std::map<std::string, double>>& hmap, const std::string& mkey);
   const std::map<std::string, double>& lumiWtMap();

   //double xsec_{-99.};
   
   //std::unordered_map<std::string, int> eventIdMap_;
   //std::map<std::string, std::map<std::string, double>> hmap_;



   float costheta1(TLorentzVector particle_orig, TLorentzVector parent_orig);
   float costhetastar(TLorentzVector particle_orig, TLorentzVector parent_orig);
   void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters=" ");
   void buildList(const std::vector<std::string>& tokens, std::vector<std::string>& list);
   //   void storeCuts(const std::vector<std::string>& tokens, std::map<std::string, std::map<std::string, double>>& hmap);
   //void buildMap(const vector<string>& tokens, map<string, int>& hmap);
   
   // for filling histograms

   // double cutValue(const std::map<std::string, double>& m, const std::string& cname);
   //const std::map<std::string, double>& cutMap(const std::map<std::string, std::map<std::string, double>>& hmap, const std::string& mkey);

   TH1* getHist1D(const char* hname);
   TH1* getHist1D(const std::string& hname);
   template <class T>
     bool fillHist1D(const char* hname, T value, double w=1.0) {
     TH1* h = getHist1D(hname);
     if (h == nullptr) return false;
     h->Fill(value, w);
     return true;
   }
   
   TH2* getHist2D(const char* hname);
   TH2* getHist2D(const std::string& hname);
   template <class T1, class T2>
     bool fillHist2D(const char* hname, T1 xvalue, T2 yvalue, double w=1.0){
     TH2* h = getHist2D(hname);
     if (h == nullptr) return false;
     h->Fill(xvalue, yvalue, w);
     return true;
   }

   //const std::map<std::string, double>& lumiWtMap() {return cutMap(hmap_, "lumiWtList");}
   //const std::map<std::string, double>& lumiWtMap();
   //  const std::map<std::string, double>& ratioCutMap();
   // const std::map<std::string, double>& SRCutMap();
   void showEfficiency(const std::string& hname,const std::vector<std::string>& slist,const std::string& header,const std::string& tag,std::ostream& os = std::cout);
   //void showEfficiency(const string& hname,const std::vector<std::string>& slist,const string& header,const string& tag);
   void objectEfficiency();
   void eventEfficiency();
   void showCategoryYield();
   
   void bookHistograms();
   void saveHistograms();

   double evtWeightSum_ {0.};
   //double intLumi = 27.070; // fb-1 (for D,E,F) for total : 41557
   //double intLumi = 35.9;
   //double xsec;
   //double xsec = 55130.0; // fb (?, from xsdb) for h signal sample
   //double xsec = 78.23; // fb drellyan
   //double xsec = 79.260 ; //pb // low mass Zg, as said by Ming
   double lumiWt(double evtWeightSum, int nent) const;
   //double lumiWt(double evtWeightSum);
   void scaleHistogram(const std::string& hname, double fac);


   //TFile *outfile = new TFile("ZtoJPsiGamma_hist_t3.root","RECREATE");
   //TFile *outfile = new TFile("DrellYan_hist.root","RECREATE");;   
   TFile *outfile ;//= TFile::Open(histFile_.c_str(), "RECREATE");
   // histf_.reset(std::move(f));

     //PileUp file
   //TString pufilename = "/afs/cern.ch/work/s/sarkar/public/Exotic/offlineAnalysis/MC/puWeightsMoriond17_v2.root";// 2016
   //TString pufilename = "mcPileup2017.root";
   //TString pufilename = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/PileupHistogram-goldenJSON-13tev-2017-69200ub.root";
   //TString pufilename =   "/afs/cern.ch/work/p/ppalit/public/Nicola_Analysis/CMSSW_9_4_9/src/HiggsAnalysis/HiggsToZZ4Leptons/test/macros/puWeightsMoriond17_v2.root";
   //TFile* fPU = new TFile(pufilename);
   // TString puhistname = "pu_mc";
   //TString puhistname = "pileup";
   //TString puhistname = "weights";// 2016 & 2017 Nicola
   //double dopuweight(int& num_pu_vt) const;
   bool bookedHistograms;
   double SFCalc(TH2F* h, double pt, double eta, int id);
   double triggSFCalc(TH2F* h, double pt_mu, double pt_pho);
   double triggSFCalcSU(TH2F* h, double pt_mu, double pt_pho);
   double SFCalcSU(TH2F* h, double pt, double eta, int id);
   template <class T>
     class PtComparatorTL {
   public:
     bool operator()(const T& a, const T& b) const {
       return a.Pt() > b.Pt();
     }
   };

};

#endif

#ifdef ZtoJPsiGamma_cxx
ZtoJPsiGamma::ZtoJPsiGamma(TTree *tree) : chain(0),chain1(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  
  if (tree == 0) {
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/user/p/ppalit/public/hztoJpsiGamma/CMSSW_9_4_8/src/offlineAnalysis_17/8EEB3BA0-9E04-AC44-AA78-ED82C143F4E3.root"); // Z signal before skimming     
    TFile *f = TFile::Open("/eos/user/s/sroychow/jpsi/anaInput/ZToJPsiGamma/8EEB3BA0-9E04-AC44-AA78-ED82C143F4E3_Skim.root"); // Z signal after skimming     
    //TFile *f = TFile::Open("/eos/user/g/gsaha4/Exotic/jpsippalit/Jan2021/y2016/crab_HToJPsiG_ToMuMuG_allProdMode_M125_TuneCUETP8M1_13TeV_Pythia8_2016_2016V2_task_2/tree_1.root"); // Z signal after skimming
 
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/user/p/ppalit/public/hztoJpsiGamma/CMSSW_9_4_8/src/offlineAnalysis_17/A897EC1C-6A36-4149-BC29-6EE7BE5DC1E7.root");   // DY background   
    if (!f || !f->IsOpen()) {
      //f = new TFile("/afs/cern.ch/user/p/ppalit/public/hztoJpsiGamma/CMSSW_9_4_8/src/offlineAnalysis_17/8EEB3BA0-9E04-AC44-AA78-ED82C143F4E3.root"); // Z signal before skimming
      f = new TFile("/eos/user/s/sroychow/jpsi/anaInput/ZToJPsiGamma/8EEB3BA0-9E04-AC44-AA78-ED82C143F4E3_Skim.root"); // Z signal after skimming
	   // f = new TFile("/eos/user/g/gsaha4/Exotic/jpsippalit/Jan2021/y2016/crab_HToJPsiG_ToMuMuG_allProdMode_M125_TuneCUETP8M1_13TeV_Pythia8_2016_2016V2_task_2/tree_1.root"); // Z signal after skimming
      //f = new TFile("/afs/cern.ch/user/p/ppalit/public/hztoJpsiGamma/CMSSW_9_4_8/src/offlineAnalysis_17/A897EC1C-6A36-4149-BC29-6EE7BE5DC1E7.root");  // DY Background
    }
  

    f->GetObject("Events",tree);
    f->GetObject("Runs",tree);
  }

   Init(tree);
   bookedHistograms = false;
}

ZtoJPsiGamma::~ZtoJPsiGamma()
{
   if (!chain) return;
   delete chain->GetCurrentFile();
   if (!chain1) return;
   delete chain1->GetCurrentFile();

}

Int_t ZtoJPsiGamma::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!chain) return 0;
   return      chain->GetEntry(entry);
}

Long64_t ZtoJPsiGamma::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!chain) return -5;
   Long64_t centry = chain->LoadTree(entry);
   if (centry < 0) return centry;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ZtoJPsiGamma::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   chain  = new TChain("Events");
   //chain = tree;
   fCurrent = -1;
   chain->SetMakeClass(1);

   chain->SetBranchAddress("run", &run, &b_run);
   chain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   chain->SetBranchAddress("event", &event, &b_event);
   chain->SetBranchAddress("puWeight", &puWeight, &b_puWeight);
   chain->SetBranchAddress("puWeightUp", &puWeightUp, &b_puWeightUp);
   chain->SetBranchAddress("puWeightDown", &puWeightDown, &b_puWeightDown);
   chain->SetBranchAddress("HTXS_Higgs_pt", &HTXS_Higgs_pt, &b_HTXS_Higgs_pt);
   chain->SetBranchAddress("HTXS_Higgs_y", &HTXS_Higgs_y, &b_HTXS_Higgs_y);
   chain->SetBranchAddress("HTXS_stage1_1_cat_pTjet25GeV", &HTXS_stage1_1_cat_pTjet25GeV, &b_HTXS_stage1_1_cat_pTjet25GeV);
   chain->SetBranchAddress("HTXS_stage1_1_cat_pTjet30GeV", &HTXS_stage1_1_cat_pTjet30GeV, &b_HTXS_stage1_1_cat_pTjet30GeV);
   chain->SetBranchAddress("HTXS_stage1_1_fine_cat_pTjet25GeV", &HTXS_stage1_1_fine_cat_pTjet25GeV, &b_HTXS_stage1_1_fine_cat_pTjet25GeV);
   chain->SetBranchAddress("HTXS_stage1_1_fine_cat_pTjet30GeV", &HTXS_stage1_1_fine_cat_pTjet30GeV, &b_HTXS_stage1_1_fine_cat_pTjet30GeV);
   chain->SetBranchAddress("HTXS_stage_0", &HTXS_stage_0, &b_HTXS_stage_0);
   chain->SetBranchAddress("HTXS_stage_1_pTjet25", &HTXS_stage_1_pTjet25, &b_HTXS_stage_1_pTjet25);
   chain->SetBranchAddress("HTXS_stage_1_pTjet30", &HTXS_stage_1_pTjet30, &b_HTXS_stage_1_pTjet30);
   chain->SetBranchAddress("HTXS_njets25", &HTXS_njets25, &b_HTXS_njets25);
   chain->SetBranchAddress("HTXS_njets30", &HTXS_njets30, &b_HTXS_njets30);
   chain->SetBranchAddress("btagWeight_CSVV2", &btagWeight_CSVV2, &b_btagWeight_CSVV2);
   chain->SetBranchAddress("btagWeight_DeepCSVB", &btagWeight_DeepCSVB, &b_btagWeight_DeepCSVB);
   chain->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi);
   chain->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt);
   chain->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt);
   chain->SetBranchAddress("ChsMET_phi", &ChsMET_phi, &b_ChsMET_phi);
   chain->SetBranchAddress("ChsMET_pt", &ChsMET_pt, &b_ChsMET_pt);
   chain->SetBranchAddress("ChsMET_sumEt", &ChsMET_sumEt, &b_ChsMET_sumEt);
   chain->SetBranchAddress("nCorrT1METJet", &nCorrT1METJet, &b_nCorrT1METJet);
   chain->SetBranchAddress("CorrT1METJet_area", CorrT1METJet_area, &b_CorrT1METJet_area);
   chain->SetBranchAddress("CorrT1METJet_eta", CorrT1METJet_eta, &b_CorrT1METJet_eta);
   chain->SetBranchAddress("CorrT1METJet_muonSubtrFactor", CorrT1METJet_muonSubtrFactor, &b_CorrT1METJet_muonSubtrFactor);
   chain->SetBranchAddress("CorrT1METJet_phi", CorrT1METJet_phi, &b_CorrT1METJet_phi);
   chain->SetBranchAddress("CorrT1METJet_rawPt", CorrT1METJet_rawPt, &b_CorrT1METJet_rawPt);
   chain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   chain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   chain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   chain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   chain->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   chain->SetBranchAddress("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP, &b_Electron_dr03TkSumPtHEEP);
   chain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   chain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   chain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   chain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   chain->SetBranchAddress("Electron_eCorr", Electron_eCorr, &b_Electron_eCorr);
   chain->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   chain->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
   chain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   chain->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
   chain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   chain->SetBranchAddress("Electron_jetPtRelv2", Electron_jetPtRelv2, &b_Electron_jetPtRelv2);
   chain->SetBranchAddress("Electron_jetRelIso", Electron_jetRelIso, &b_Electron_jetRelIso);
   chain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   chain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   chain->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   chain->SetBranchAddress("Electron_mvaFall17V1Iso", Electron_mvaFall17V1Iso, &b_Electron_mvaFall17V1Iso);
   chain->SetBranchAddress("Electron_mvaFall17V1noIso", Electron_mvaFall17V1noIso, &b_Electron_mvaFall17V1noIso);
   chain->SetBranchAddress("Electron_mvaFall17V2Iso", Electron_mvaFall17V2Iso, &b_Electron_mvaFall17V2Iso);
   chain->SetBranchAddress("Electron_mvaFall17V2noIso", Electron_mvaFall17V2noIso, &b_Electron_mvaFall17V2noIso);
   chain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   chain->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   chain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   chain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   chain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   chain->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
   chain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   chain->SetBranchAddress("Electron_mvaTTH", Electron_mvaTTH, &b_Electron_mvaTTH);
   chain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   chain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   chain->SetBranchAddress("Electron_cutBased_Fall17_V1", Electron_cutBased_Fall17_V1, &b_Electron_cutBased_Fall17_V1);
   chain->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
   chain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   chain->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
   chain->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
   chain->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
   chain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   chain->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   chain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
   chain->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
   chain->SetBranchAddress("Electron_mvaFall17V1Iso_WP80", Electron_mvaFall17V1Iso_WP80, &b_Electron_mvaFall17V1Iso_WP80);
   chain->SetBranchAddress("Electron_mvaFall17V1Iso_WP90", Electron_mvaFall17V1Iso_WP90, &b_Electron_mvaFall17V1Iso_WP90);
   chain->SetBranchAddress("Electron_mvaFall17V1Iso_WPL", Electron_mvaFall17V1Iso_WPL, &b_Electron_mvaFall17V1Iso_WPL);
   chain->SetBranchAddress("Electron_mvaFall17V1noIso_WP80", Electron_mvaFall17V1noIso_WP80, &b_Electron_mvaFall17V1noIso_WP80);
   chain->SetBranchAddress("Electron_mvaFall17V1noIso_WP90", Electron_mvaFall17V1noIso_WP90, &b_Electron_mvaFall17V1noIso_WP90);
   chain->SetBranchAddress("Electron_mvaFall17V1noIso_WPL", Electron_mvaFall17V1noIso_WPL, &b_Electron_mvaFall17V1noIso_WPL);
   chain->SetBranchAddress("Electron_mvaFall17V2Iso_WP80", Electron_mvaFall17V2Iso_WP80, &b_Electron_mvaFall17V2Iso_WP80);
   chain->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", Electron_mvaFall17V2Iso_WP90, &b_Electron_mvaFall17V2Iso_WP90);
   chain->SetBranchAddress("Electron_mvaFall17V2Iso_WPL", Electron_mvaFall17V2Iso_WPL, &b_Electron_mvaFall17V2Iso_WPL);
   chain->SetBranchAddress("Electron_mvaFall17V2noIso_WP80", Electron_mvaFall17V2noIso_WP80, &b_Electron_mvaFall17V2noIso_WP80);
   chain->SetBranchAddress("Electron_mvaFall17V2noIso_WP90", Electron_mvaFall17V2noIso_WP90, &b_Electron_mvaFall17V2noIso_WP90);
   chain->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", Electron_mvaFall17V2noIso_WPL, &b_Electron_mvaFall17V2noIso_WPL);
   chain->SetBranchAddress("Electron_seedGain", Electron_seedGain, &b_Electron_seedGain);
   chain->SetBranchAddress("Flag_ecalBadCalibFilterV2", &Flag_ecalBadCalibFilterV2, &b_Flag_ecalBadCalibFilterV2);
   chain->SetBranchAddress("nFatJet", &nFatJet, &b_nFatJet);
   chain->SetBranchAddress("FatJet_area", FatJet_area, &b_FatJet_area);
   chain->SetBranchAddress("FatJet_btagCMVA", FatJet_btagCMVA, &b_FatJet_btagCMVA);
   chain->SetBranchAddress("FatJet_btagCSVV2", FatJet_btagCSVV2, &b_FatJet_btagCSVV2);
   chain->SetBranchAddress("FatJet_btagDDBvL", FatJet_btagDDBvL, &b_FatJet_btagDDBvL);
   chain->SetBranchAddress("FatJet_btagDDCvB", FatJet_btagDDCvB, &b_FatJet_btagDDCvB);
   chain->SetBranchAddress("FatJet_btagDDCvL", FatJet_btagDDCvL, &b_FatJet_btagDDCvL);
   chain->SetBranchAddress("FatJet_btagDeepB", FatJet_btagDeepB, &b_FatJet_btagDeepB);
   chain->SetBranchAddress("FatJet_btagHbb", FatJet_btagHbb, &b_FatJet_btagHbb);
   chain->SetBranchAddress("FatJet_deepTagMD_H4qvsQCD", FatJet_deepTagMD_H4qvsQCD, &b_FatJet_deepTagMD_H4qvsQCD);
   chain->SetBranchAddress("FatJet_deepTagMD_HbbvsQCD", FatJet_deepTagMD_HbbvsQCD, &b_FatJet_deepTagMD_HbbvsQCD);
   chain->SetBranchAddress("FatJet_deepTagMD_TvsQCD", FatJet_deepTagMD_TvsQCD, &b_FatJet_deepTagMD_TvsQCD);
   chain->SetBranchAddress("FatJet_deepTagMD_WvsQCD", FatJet_deepTagMD_WvsQCD, &b_FatJet_deepTagMD_WvsQCD);
   chain->SetBranchAddress("FatJet_deepTagMD_ZHbbvsQCD", FatJet_deepTagMD_ZHbbvsQCD, &b_FatJet_deepTagMD_ZHbbvsQCD);
   chain->SetBranchAddress("FatJet_deepTagMD_ZHccvsQCD", FatJet_deepTagMD_ZHccvsQCD, &b_FatJet_deepTagMD_ZHccvsQCD);
   chain->SetBranchAddress("FatJet_deepTagMD_ZbbvsQCD", FatJet_deepTagMD_ZbbvsQCD, &b_FatJet_deepTagMD_ZbbvsQCD);
   chain->SetBranchAddress("FatJet_deepTagMD_ZvsQCD", FatJet_deepTagMD_ZvsQCD, &b_FatJet_deepTagMD_ZvsQCD);
   chain->SetBranchAddress("FatJet_deepTagMD_bbvsLight", FatJet_deepTagMD_bbvsLight, &b_FatJet_deepTagMD_bbvsLight);
   chain->SetBranchAddress("FatJet_deepTagMD_ccvsLight", FatJet_deepTagMD_ccvsLight, &b_FatJet_deepTagMD_ccvsLight);
   chain->SetBranchAddress("FatJet_deepTag_H", FatJet_deepTag_H, &b_FatJet_deepTag_H);
   chain->SetBranchAddress("FatJet_deepTag_QCD", FatJet_deepTag_QCD, &b_FatJet_deepTag_QCD);
   chain->SetBranchAddress("FatJet_deepTag_QCDothers", FatJet_deepTag_QCDothers, &b_FatJet_deepTag_QCDothers);
   chain->SetBranchAddress("FatJet_deepTag_TvsQCD", FatJet_deepTag_TvsQCD, &b_FatJet_deepTag_TvsQCD);
   chain->SetBranchAddress("FatJet_deepTag_WvsQCD", FatJet_deepTag_WvsQCD, &b_FatJet_deepTag_WvsQCD);
   chain->SetBranchAddress("FatJet_deepTag_ZvsQCD", FatJet_deepTag_ZvsQCD, &b_FatJet_deepTag_ZvsQCD);
   chain->SetBranchAddress("FatJet_eta", FatJet_eta, &b_FatJet_eta);
   chain->SetBranchAddress("FatJet_mass", FatJet_mass, &b_FatJet_mass);
   chain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop, &b_FatJet_msoftdrop);
   chain->SetBranchAddress("FatJet_n2b1", FatJet_n2b1, &b_FatJet_n2b1);
   chain->SetBranchAddress("FatJet_n3b1", FatJet_n3b1, &b_FatJet_n3b1);
   chain->SetBranchAddress("FatJet_phi", FatJet_phi, &b_FatJet_phi);
   chain->SetBranchAddress("FatJet_pt", FatJet_pt, &b_FatJet_pt);
   chain->SetBranchAddress("FatJet_rawFactor", FatJet_rawFactor, &b_FatJet_rawFactor);
   chain->SetBranchAddress("FatJet_tau1", FatJet_tau1, &b_FatJet_tau1);
   chain->SetBranchAddress("FatJet_tau2", FatJet_tau2, &b_FatJet_tau2);
   chain->SetBranchAddress("FatJet_tau3", FatJet_tau3, &b_FatJet_tau3);
   chain->SetBranchAddress("FatJet_tau4", FatJet_tau4, &b_FatJet_tau4);
   chain->SetBranchAddress("FatJet_jetId", FatJet_jetId, &b_FatJet_jetId);
   chain->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1, &b_FatJet_subJetIdx1);
   chain->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2, &b_FatJet_subJetIdx2);
   chain->SetBranchAddress("nGenJetAK8", &nGenJetAK8, &b_nGenJetAK8);
   chain->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta, &b_GenJetAK8_eta);
   chain->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass, &b_GenJetAK8_mass);
   chain->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi, &b_GenJetAK8_phi);
   chain->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt, &b_GenJetAK8_pt);
   chain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   chain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
   chain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
   chain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
   chain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
   chain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   chain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   chain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   chain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   chain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   chain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   chain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   chain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   chain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   chain->SetBranchAddress("nSubGenJetAK8", &nSubGenJetAK8, &b_nSubGenJetAK8);
   chain->SetBranchAddress("SubGenJetAK8_eta", SubGenJetAK8_eta, &b_SubGenJetAK8_eta);
   chain->SetBranchAddress("SubGenJetAK8_mass", SubGenJetAK8_mass, &b_SubGenJetAK8_mass);
   chain->SetBranchAddress("SubGenJetAK8_phi", SubGenJetAK8_phi, &b_SubGenJetAK8_phi);
   chain->SetBranchAddress("SubGenJetAK8_pt", SubGenJetAK8_pt, &b_SubGenJetAK8_pt);
   chain->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar);
   chain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   chain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   chain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   chain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   chain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   chain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   chain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   chain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   chain->SetBranchAddress("nGenVisTau", &nGenVisTau, &b_nGenVisTau);
   chain->SetBranchAddress("GenVisTau_eta", GenVisTau_eta, &b_GenVisTau_eta);
   chain->SetBranchAddress("GenVisTau_mass", GenVisTau_mass, &b_GenVisTau_mass);
   chain->SetBranchAddress("GenVisTau_phi", GenVisTau_phi, &b_GenVisTau_phi);
   chain->SetBranchAddress("GenVisTau_pt", GenVisTau_pt, &b_GenVisTau_pt);
   chain->SetBranchAddress("GenVisTau_charge", GenVisTau_charge, &b_GenVisTau_charge);
   chain->SetBranchAddress("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother, &b_GenVisTau_genPartIdxMother);
   chain->SetBranchAddress("GenVisTau_status", GenVisTau_status, &b_GenVisTau_status);
   chain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   chain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   chain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
   chain->SetBranchAddress("nIsoTrack", &nIsoTrack, &b_nIsoTrack);
   chain->SetBranchAddress("IsoTrack_dxy", IsoTrack_dxy, &b_IsoTrack_dxy);
   chain->SetBranchAddress("IsoTrack_dz", IsoTrack_dz, &b_IsoTrack_dz);
   chain->SetBranchAddress("IsoTrack_eta", IsoTrack_eta, &b_IsoTrack_eta);
   chain->SetBranchAddress("IsoTrack_pfRelIso03_all", IsoTrack_pfRelIso03_all, &b_IsoTrack_pfRelIso03_all);
   chain->SetBranchAddress("IsoTrack_pfRelIso03_chg", IsoTrack_pfRelIso03_chg, &b_IsoTrack_pfRelIso03_chg);
   chain->SetBranchAddress("IsoTrack_phi", IsoTrack_phi, &b_IsoTrack_phi);
   chain->SetBranchAddress("IsoTrack_pt", IsoTrack_pt, &b_IsoTrack_pt);
   chain->SetBranchAddress("IsoTrack_miniPFRelIso_all", IsoTrack_miniPFRelIso_all, &b_IsoTrack_miniPFRelIso_all);
   chain->SetBranchAddress("IsoTrack_miniPFRelIso_chg", IsoTrack_miniPFRelIso_chg, &b_IsoTrack_miniPFRelIso_chg);
   chain->SetBranchAddress("IsoTrack_fromPV", IsoTrack_fromPV, &b_IsoTrack_fromPV);
   chain->SetBranchAddress("IsoTrack_pdgId", IsoTrack_pdgId, &b_IsoTrack_pdgId);
   chain->SetBranchAddress("IsoTrack_isHighPurityTrack", IsoTrack_isHighPurityTrack, &b_IsoTrack_isHighPurityTrack);
   chain->SetBranchAddress("IsoTrack_isPFcand", IsoTrack_isPFcand, &b_IsoTrack_isPFcand);
   chain->SetBranchAddress("IsoTrack_isFromLostTrack", IsoTrack_isFromLostTrack, &b_IsoTrack_isFromLostTrack);
   chain->SetBranchAddress("nJet", &nJet, &b_nJet);
   chain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   chain->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA, &b_Jet_btagCMVA);
   chain->SetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2, &b_Jet_btagCSVV2);
   chain->SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB, &b_Jet_btagDeepB);
   chain->SetBranchAddress("Jet_btagDeepC", Jet_btagDeepC, &b_Jet_btagDeepC);
   chain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   chain->SetBranchAddress("Jet_btagDeepFlavC", Jet_btagDeepFlavC, &b_Jet_btagDeepFlavC);
   chain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   chain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   chain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   chain->SetBranchAddress("Jet_jercCHF", Jet_jercCHF, &b_Jet_jercCHF);
   chain->SetBranchAddress("Jet_jercCHPUF", Jet_jercCHPUF, &b_Jet_jercCHPUF);
   chain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   chain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   chain->SetBranchAddress("Jet_muonSubtrFactor", Jet_muonSubtrFactor, &b_Jet_muonSubtrFactor);
   chain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   chain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   chain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   chain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   chain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   chain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   chain->SetBranchAddress("Jet_bRegCorr", Jet_bRegCorr, &b_Jet_bRegCorr);
   chain->SetBranchAddress("Jet_bRegRes", Jet_bRegRes, &b_Jet_bRegRes);
   chain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   chain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   chain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   chain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   chain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   chain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   chain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   chain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   chain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   chain->SetBranchAddress("L1PreFiringWeight_Dn", &L1PreFiringWeight_Dn, &b_L1PreFiringWeight_Dn);
   chain->SetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom, &b_L1PreFiringWeight_Nom);
   chain->SetBranchAddress("L1PreFiringWeight_Up", &L1PreFiringWeight_Up, &b_L1PreFiringWeight_Up);
   chain->SetBranchAddress("METFixEE2017_MetUnclustEnUpDeltaX", &METFixEE2017_MetUnclustEnUpDeltaX, &b_METFixEE2017_MetUnclustEnUpDeltaX);
   chain->SetBranchAddress("METFixEE2017_MetUnclustEnUpDeltaY", &METFixEE2017_MetUnclustEnUpDeltaY, &b_METFixEE2017_MetUnclustEnUpDeltaY);
   chain->SetBranchAddress("METFixEE2017_covXX", &METFixEE2017_covXX, &b_METFixEE2017_covXX);
   chain->SetBranchAddress("METFixEE2017_covXY", &METFixEE2017_covXY, &b_METFixEE2017_covXY);
   chain->SetBranchAddress("METFixEE2017_covYY", &METFixEE2017_covYY, &b_METFixEE2017_covYY);
   chain->SetBranchAddress("METFixEE2017_phi", &METFixEE2017_phi, &b_METFixEE2017_phi);
   chain->SetBranchAddress("METFixEE2017_pt", &METFixEE2017_pt, &b_METFixEE2017_pt);
   chain->SetBranchAddress("METFixEE2017_significance", &METFixEE2017_significance, &b_METFixEE2017_significance);
   chain->SetBranchAddress("METFixEE2017_sumEt", &METFixEE2017_sumEt, &b_METFixEE2017_sumEt);
   chain->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi);
   chain->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt);
   chain->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
   chain->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
   chain->SetBranchAddress("MET_covXX", &MET_covXX, &b_MET_covXX);
   chain->SetBranchAddress("MET_covXY", &MET_covXY, &b_MET_covXY);
   chain->SetBranchAddress("MET_covYY", &MET_covYY, &b_MET_covYY);
   chain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   chain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   chain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   chain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   chain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   chain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   chain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   chain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   chain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   chain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   chain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   chain->SetBranchAddress("Muon_jetPtRelv2", Muon_jetPtRelv2, &b_Muon_jetPtRelv2);
   chain->SetBranchAddress("Muon_jetRelIso", Muon_jetRelIso, &b_Muon_jetRelIso);
   chain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   chain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   chain->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   chain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   chain->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   chain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   chain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   chain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   chain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   chain->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
   chain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   chain->SetBranchAddress("Muon_tkRelIso", Muon_tkRelIso, &b_Muon_tkRelIso);
   chain->SetBranchAddress("Muon_tunepRelPt", Muon_tunepRelPt, &b_Muon_tunepRelPt);
   chain->SetBranchAddress("Muon_mvaLowPt", Muon_mvaLowPt, &b_Muon_mvaLowPt);
   chain->SetBranchAddress("Muon_mvaTTH", Muon_mvaTTH, &b_Muon_mvaTTH);
   chain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   chain->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
   chain->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
   chain->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   chain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   chain->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
   chain->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
   chain->SetBranchAddress("Muon_inTimeMuon", Muon_inTimeMuon, &b_Muon_inTimeMuon);
   chain->SetBranchAddress("Muon_isGlobal", Muon_isGlobal, &b_Muon_isGlobal);
   chain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   chain->SetBranchAddress("Muon_isTracker", Muon_isTracker, &b_Muon_isTracker);
   chain->SetBranchAddress("Muon_looseId", Muon_looseId, &b_Muon_looseId);
   chain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   chain->SetBranchAddress("Muon_mediumPromptId", Muon_mediumPromptId, &b_Muon_mediumPromptId);
   chain->SetBranchAddress("Muon_miniIsoId", Muon_miniIsoId, &b_Muon_miniIsoId);
   chain->SetBranchAddress("Muon_multiIsoId", Muon_multiIsoId, &b_Muon_multiIsoId);
   chain->SetBranchAddress("Muon_mvaId", Muon_mvaId, &b_Muon_mvaId);
   chain->SetBranchAddress("Muon_pfIsoId", Muon_pfIsoId, &b_Muon_pfIsoId);
   chain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   chain->SetBranchAddress("Muon_softMvaId", Muon_softMvaId, &b_Muon_softMvaId);
   chain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   chain->SetBranchAddress("Muon_tkIsoId", Muon_tkIsoId, &b_Muon_tkIsoId);
   chain->SetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose, &b_Muon_triggerIdLoose);
   chain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   chain->SetBranchAddress("Photon_eCorr", Photon_eCorr, &b_Photon_eCorr);
   chain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   chain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   chain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   chain->SetBranchAddress("Photon_mass", Photon_mass, &b_Photon_mass);
   chain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   chain->SetBranchAddress("Photon_mvaIDV1", Photon_mvaIDV1, &b_Photon_mvaIDV1);
   chain->SetBranchAddress("Photon_pfRelIso03_all", Photon_pfRelIso03_all, &b_Photon_pfRelIso03_all);
   chain->SetBranchAddress("Photon_pfRelIso03_chg", Photon_pfRelIso03_chg, &b_Photon_pfRelIso03_chg);
   chain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   chain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   chain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   chain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   chain->SetBranchAddress("Photon_charge", Photon_charge, &b_Photon_charge);
   chain->SetBranchAddress("Photon_cutBasedBitmap", Photon_cutBasedBitmap, &b_Photon_cutBasedBitmap);
   chain->SetBranchAddress("Photon_cutBasedV1Bitmap", Photon_cutBasedV1Bitmap, &b_Photon_cutBasedV1Bitmap);
   chain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   chain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   chain->SetBranchAddress("Photon_pdgId", Photon_pdgId, &b_Photon_pdgId);
   chain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   chain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   chain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   chain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   chain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   chain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   chain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   chain->SetBranchAddress("Photon_seedGain", Photon_seedGain, &b_Photon_seedGain);
   chain->SetBranchAddress("Photon_dEsigmaDown",Photon_dEsigmaDown, &b_Photon_dEsigmaDown);
   chain->SetBranchAddress("Photon_dEsigmaUp", Photon_dEsigmaUp, &b_Photon_dEsigmaUp);
   chain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   chain->SetBranchAddress("Pileup_pudensity", &Pileup_pudensity, &b_Pileup_pudensity);
   chain->SetBranchAddress("Pileup_gpudensity", &Pileup_gpudensity, &b_Pileup_gpudensity);
   chain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   chain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   chain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   chain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   chain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   chain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   chain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   chain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   chain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
   chain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   chain->SetBranchAddress("fixedGridRhoFastjetCentral", &fixedGridRhoFastjetCentral, &b_fixedGridRhoFastjetCentral);
   chain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   chain->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, &b_fixedGridRhoFastjetCentralChargedPileUp);
   chain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   chain->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton, &b_nGenDressedLepton);
   chain->SetBranchAddress("GenDressedLepton_eta", GenDressedLepton_eta, &b_GenDressedLepton_eta);
   chain->SetBranchAddress("GenDressedLepton_mass", GenDressedLepton_mass, &b_GenDressedLepton_mass);
   chain->SetBranchAddress("GenDressedLepton_phi", GenDressedLepton_phi, &b_GenDressedLepton_phi);
   chain->SetBranchAddress("GenDressedLepton_pt", GenDressedLepton_pt, &b_GenDressedLepton_pt);
   chain->SetBranchAddress("GenDressedLepton_pdgId", GenDressedLepton_pdgId, &b_GenDressedLepton_pdgId);
   chain->SetBranchAddress("GenDressedLepton_hasTauAnc", GenDressedLepton_hasTauAnc, &b_GenDressedLepton_hasTauAnc);
   chain->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet, &b_nSoftActivityJet);
   chain->SetBranchAddress("SoftActivityJet_eta", SoftActivityJet_eta, &b_SoftActivityJet_eta);
   chain->SetBranchAddress("SoftActivityJet_phi", SoftActivityJet_phi, &b_SoftActivityJet_phi);
   chain->SetBranchAddress("SoftActivityJet_pt", SoftActivityJet_pt, &b_SoftActivityJet_pt);
   chain->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT, &b_SoftActivityJetHT);
   chain->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10, &b_SoftActivityJetHT10);
   chain->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2, &b_SoftActivityJetHT2);
   chain->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5, &b_SoftActivityJetHT5);
   chain->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10, &b_SoftActivityJetNjets10);
   chain->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2, &b_SoftActivityJetNjets2);
   chain->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5, &b_SoftActivityJetNjets5);
   chain->SetBranchAddress("nSubJet", &nSubJet, &b_nSubJet);
   chain->SetBranchAddress("SubJet_btagCMVA", SubJet_btagCMVA, &b_SubJet_btagCMVA);
   chain->SetBranchAddress("SubJet_btagCSVV2", SubJet_btagCSVV2, &b_SubJet_btagCSVV2);
   chain->SetBranchAddress("SubJet_btagDeepB", SubJet_btagDeepB, &b_SubJet_btagDeepB);
   chain->SetBranchAddress("SubJet_eta", SubJet_eta, &b_SubJet_eta);
   chain->SetBranchAddress("SubJet_mass", SubJet_mass, &b_SubJet_mass);
   chain->SetBranchAddress("SubJet_n2b1", SubJet_n2b1, &b_SubJet_n2b1);
   chain->SetBranchAddress("SubJet_n3b1", SubJet_n3b1, &b_SubJet_n3b1);
   chain->SetBranchAddress("SubJet_phi", SubJet_phi, &b_SubJet_phi);
   chain->SetBranchAddress("SubJet_pt", SubJet_pt, &b_SubJet_pt);
   chain->SetBranchAddress("SubJet_rawFactor", SubJet_rawFactor, &b_SubJet_rawFactor);
   chain->SetBranchAddress("SubJet_tau1", SubJet_tau1, &b_SubJet_tau1);
   chain->SetBranchAddress("SubJet_tau2", SubJet_tau2, &b_SubJet_tau2);
   chain->SetBranchAddress("SubJet_tau3", SubJet_tau3, &b_SubJet_tau3);
   chain->SetBranchAddress("SubJet_tau4", SubJet_tau4, &b_SubJet_tau4);
   chain->SetBranchAddress("nTau", &nTau, &b_nTau);
   chain->SetBranchAddress("Tau_chargedIso", Tau_chargedIso, &b_Tau_chargedIso);
   chain->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
   chain->SetBranchAddress("Tau_dz", Tau_dz, &b_Tau_dz);
   chain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   chain->SetBranchAddress("Tau_leadTkDeltaEta", Tau_leadTkDeltaEta, &b_Tau_leadTkDeltaEta);
   chain->SetBranchAddress("Tau_leadTkDeltaPhi", Tau_leadTkDeltaPhi, &b_Tau_leadTkDeltaPhi);
   chain->SetBranchAddress("Tau_leadTkPtOverTauPt", Tau_leadTkPtOverTauPt, &b_Tau_leadTkPtOverTauPt);
   chain->SetBranchAddress("Tau_mass", Tau_mass, &b_Tau_mass);
   chain->SetBranchAddress("Tau_neutralIso", Tau_neutralIso, &b_Tau_neutralIso);
   chain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   chain->SetBranchAddress("Tau_photonsOutsideSignalCone", Tau_photonsOutsideSignalCone, &b_Tau_photonsOutsideSignalCone);
   chain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   chain->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
   chain->SetBranchAddress("Tau_rawAntiEle", Tau_rawAntiEle, &b_Tau_rawAntiEle);
   chain->SetBranchAddress("Tau_rawAntiEle2018", Tau_rawAntiEle2018, &b_Tau_rawAntiEle2018);
   chain->SetBranchAddress("Tau_rawIso", Tau_rawIso, &b_Tau_rawIso);
   chain->SetBranchAddress("Tau_rawIsodR03", Tau_rawIsodR03, &b_Tau_rawIsodR03);
   chain->SetBranchAddress("Tau_rawMVAnewDM2017v2", Tau_rawMVAnewDM2017v2, &b_Tau_rawMVAnewDM2017v2);
   chain->SetBranchAddress("Tau_rawMVAoldDM", Tau_rawMVAoldDM, &b_Tau_rawMVAoldDM);
   chain->SetBranchAddress("Tau_rawMVAoldDM2017v1", Tau_rawMVAoldDM2017v1, &b_Tau_rawMVAoldDM2017v1);
   chain->SetBranchAddress("Tau_rawMVAoldDM2017v2", Tau_rawMVAoldDM2017v2, &b_Tau_rawMVAoldDM2017v2);
   chain->SetBranchAddress("Tau_rawMVAoldDMdR032017v2", Tau_rawMVAoldDMdR032017v2, &b_Tau_rawMVAoldDMdR032017v2);
   chain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   chain->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
   chain->SetBranchAddress("Tau_jetIdx", Tau_jetIdx, &b_Tau_jetIdx);
   chain->SetBranchAddress("Tau_rawAntiEleCat", Tau_rawAntiEleCat, &b_Tau_rawAntiEleCat);
   chain->SetBranchAddress("Tau_rawAntiEleCat2018", Tau_rawAntiEleCat2018, &b_Tau_rawAntiEleCat2018);
   chain->SetBranchAddress("Tau_idAntiEle", Tau_idAntiEle, &b_Tau_idAntiEle);
   chain->SetBranchAddress("Tau_idAntiEle2018", Tau_idAntiEle2018, &b_Tau_idAntiEle2018);
   chain->SetBranchAddress("Tau_idAntiMu", Tau_idAntiMu, &b_Tau_idAntiMu);
   chain->SetBranchAddress("Tau_idDecayMode", Tau_idDecayMode, &b_Tau_idDecayMode);
   chain->SetBranchAddress("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs, &b_Tau_idDecayModeNewDMs);
   chain->SetBranchAddress("Tau_idMVAnewDM2017v2", Tau_idMVAnewDM2017v2, &b_Tau_idMVAnewDM2017v2);
   chain->SetBranchAddress("Tau_idMVAoldDM", Tau_idMVAoldDM, &b_Tau_idMVAoldDM);
   chain->SetBranchAddress("Tau_idMVAoldDM2017v1", Tau_idMVAoldDM2017v1, &b_Tau_idMVAoldDM2017v1);
   chain->SetBranchAddress("Tau_idMVAoldDM2017v2", Tau_idMVAoldDM2017v2, &b_Tau_idMVAoldDM2017v2);
   chain->SetBranchAddress("Tau_idMVAoldDMdR032017v2", Tau_idMVAoldDMdR032017v2, &b_Tau_idMVAoldDMdR032017v2);
   chain->SetBranchAddress("TkMET_phi", &TkMET_phi, &b_TkMET_phi);
   chain->SetBranchAddress("TkMET_pt", &TkMET_pt, &b_TkMET_pt);
   chain->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt, &b_TkMET_sumEt);
   chain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   chain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   chain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   chain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   chain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   chain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   chain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   chain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   chain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   chain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   chain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   chain->SetBranchAddress("genTtbarId", &genTtbarId, &b_genTtbarId);
   chain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   chain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   chain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   chain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   chain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   chain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   chain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   chain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   chain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   chain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   chain->SetBranchAddress("nSV", &nSV, &b_nSV);
   chain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   chain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   chain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   chain->SetBranchAddress("Electron_genPartIdx", Electron_genPartIdx, &b_Electron_genPartIdx);
   chain->SetBranchAddress("Electron_genPartFlav", Electron_genPartFlav, &b_Electron_genPartFlav);
   chain->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour, &b_GenJetAK8_partonFlavour);
   chain->SetBranchAddress("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour, &b_GenJetAK8_hadronFlavour);
   chain->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);
   chain->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour, &b_GenJet_hadronFlavour);
   chain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
   chain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   chain->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
   chain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
   chain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
   chain->SetBranchAddress("Photon_genPartIdx", Photon_genPartIdx, &b_Photon_genPartIdx);
   chain->SetBranchAddress("Photon_genPartFlav", Photon_genPartFlav, &b_Photon_genPartFlav);
   chain->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi, &b_MET_fiducialGenPhi);
   chain->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt, &b_MET_fiducialGenPt);
   chain->SetBranchAddress("Electron_cleanmask", Electron_cleanmask, &b_Electron_cleanmask);
   chain->SetBranchAddress("Jet_cleanmask", Jet_cleanmask, &b_Jet_cleanmask);
   chain->SetBranchAddress("Muon_cleanmask", Muon_cleanmask, &b_Muon_cleanmask);
   chain->SetBranchAddress("Photon_cleanmask", Photon_cleanmask, &b_Photon_cleanmask);
   chain->SetBranchAddress("Tau_cleanmask", Tau_cleanmask, &b_Tau_cleanmask);
   chain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   chain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   chain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   chain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   chain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   chain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   chain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   chain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   chain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   chain->SetBranchAddress("Tau_genPartIdx", Tau_genPartIdx, &b_Tau_genPartIdx);
   chain->SetBranchAddress("Tau_genPartFlav", Tau_genPartFlav, &b_Tau_genPartFlav);
   chain->SetBranchAddress("L1simulation_step", &L1simulation_step, &b_L1simulation_step);

   //chain->SetBranchAddress("TrigObj", &TrigObj, &b_TrigObj);
   chain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   chain->SetBranchAddress("HLT_DoubleMu20_7_Mass0to30_Photon23", &HLT_DoubleMu20_7_Mass0to30_Photon23, &b_HLT_DoubleMu20_7_Mass0to30_Photon23);
   chain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi", &HLT_DoubleMu4_3_Jpsi, &b_HLT_DoubleMu4_3_Jpsi);
   chain->SetBranchAddress("HLT_Mu17_Photon30_IsoCaloId", &HLT_Mu17_Photon30_IsoCaloId, &b_HLT_Mu17_Photon30_IsoCaloId);
   chain->SetBranchAddress("HLT_Mu17_Photon30_CaloIdL_L1ISO", &HLT_Mu17_Photon30_CaloIdL_L1ISO, &b_HLT_Mu17_Photon30_CaloIdL_L1ISO);
   chain->SetBranchAddress("HLT_SingleMu24", &HLT_SingleMu24, &b_HLT_SingleMu24);
   chain->SetBranchAddress("HLT_SingleMu27", &HLT_SingleMu27, &b_HLT_SingleMu27);
   //   if(hlt()==3) chain->SetBranchAddress("HLT_Dimuon20_Jpsi", &HLT_Dimuon20_Jpsi, &b_HLT_Dimuon20_Jpsi);
   chain->SetBranchAddress("HLT_Dimuon25_Jpsi", &HLT_Dimuon25_Jpsi, &b_HLT_Dimuon25_Jpsi);
      chain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   
      chain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   chain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   chain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   chain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   chain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   chain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   chain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   chain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   chain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   chain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   chain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   chain->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter, &b_Flag_ecalBadCalibFilter);
   chain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   chain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   chain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   chain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   chain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   chain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   chain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   chain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   chain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   chain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   chain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   chain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   chain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   chain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   chain->SetBranchAddress("L1_AlwaysTrue", &L1_AlwaysTrue, &b_L1_AlwaysTrue);
   chain->SetBranchAddress("L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME, &b_L1_BPTX_AND_Ref1_VME);
   chain->SetBranchAddress("L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME, &b_L1_BPTX_AND_Ref3_VME);
   chain->SetBranchAddress("L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME, &b_L1_BPTX_AND_Ref4_VME);
   chain->SetBranchAddress("L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME, &b_L1_BPTX_BeamGas_B1_VME);
   chain->SetBranchAddress("L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME, &b_L1_BPTX_BeamGas_B2_VME);
   chain->SetBranchAddress("L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME, &b_L1_BPTX_BeamGas_Ref1_VME);
   chain->SetBranchAddress("L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME, &b_L1_BPTX_BeamGas_Ref2_VME);
   chain->SetBranchAddress("L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME, &b_L1_BPTX_NotOR_VME);
   chain->SetBranchAddress("L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME, &b_L1_BPTX_OR_Ref3_VME);
   chain->SetBranchAddress("L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME, &b_L1_BPTX_OR_Ref4_VME);
   chain->SetBranchAddress("L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME, &b_L1_BPTX_RefAND_VME);
   chain->SetBranchAddress("L1_BptxMinus", &L1_BptxMinus, &b_L1_BptxMinus);
   chain->SetBranchAddress("L1_BptxOR", &L1_BptxOR, &b_L1_BptxOR);
   chain->SetBranchAddress("L1_BptxPlus", &L1_BptxPlus, &b_L1_BptxPlus);
   chain->SetBranchAddress("L1_BptxXOR", &L1_BptxXOR, &b_L1_BptxXOR);
   chain->SetBranchAddress("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142, &b_L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142);
   chain->SetBranchAddress("L1_DoubleEG6_HTT240er", &L1_DoubleEG6_HTT240er, &b_L1_DoubleEG6_HTT240er);
   chain->SetBranchAddress("L1_DoubleEG6_HTT250er", &L1_DoubleEG6_HTT250er, &b_L1_DoubleEG6_HTT250er);
   chain->SetBranchAddress("L1_DoubleEG6_HTT255er", &L1_DoubleEG6_HTT255er, &b_L1_DoubleEG6_HTT255er);
   chain->SetBranchAddress("L1_DoubleEG6_HTT270er", &L1_DoubleEG6_HTT270er, &b_L1_DoubleEG6_HTT270er);
   chain->SetBranchAddress("L1_DoubleEG6_HTT300er", &L1_DoubleEG6_HTT300er, &b_L1_DoubleEG6_HTT300er);
   chain->SetBranchAddress("L1_DoubleEG8er2p6_HTT255er", &L1_DoubleEG8er2p6_HTT255er, &b_L1_DoubleEG8er2p6_HTT255er);
   chain->SetBranchAddress("L1_DoubleEG8er2p6_HTT270er", &L1_DoubleEG8er2p6_HTT270er, &b_L1_DoubleEG8er2p6_HTT270er);
   chain->SetBranchAddress("L1_DoubleEG8er2p6_HTT300er", &L1_DoubleEG8er2p6_HTT300er, &b_L1_DoubleEG8er2p6_HTT300er);
   chain->SetBranchAddress("L1_DoubleEG_15_10", &L1_DoubleEG_15_10, &b_L1_DoubleEG_15_10);
   chain->SetBranchAddress("L1_DoubleEG_18_17", &L1_DoubleEG_18_17, &b_L1_DoubleEG_18_17);
   chain->SetBranchAddress("L1_DoubleEG_20_18", &L1_DoubleEG_20_18, &b_L1_DoubleEG_20_18);
   chain->SetBranchAddress("L1_DoubleEG_22_10", &L1_DoubleEG_22_10, &b_L1_DoubleEG_22_10);
   chain->SetBranchAddress("L1_DoubleEG_22_12", &L1_DoubleEG_22_12, &b_L1_DoubleEG_22_12);
   chain->SetBranchAddress("L1_DoubleEG_22_15", &L1_DoubleEG_22_15, &b_L1_DoubleEG_22_15);
   chain->SetBranchAddress("L1_DoubleEG_23_10", &L1_DoubleEG_23_10, &b_L1_DoubleEG_23_10);
   chain->SetBranchAddress("L1_DoubleEG_24_17", &L1_DoubleEG_24_17, &b_L1_DoubleEG_24_17);
   chain->SetBranchAddress("L1_DoubleEG_25_12", &L1_DoubleEG_25_12, &b_L1_DoubleEG_25_12);
   chain->SetBranchAddress("L1_DoubleEG_25_13", &L1_DoubleEG_25_13, &b_L1_DoubleEG_25_13);
   chain->SetBranchAddress("L1_DoubleEG_25_14", &L1_DoubleEG_25_14, &b_L1_DoubleEG_25_14);
   chain->SetBranchAddress("L1_DoubleEG_LooseIso23_10", &L1_DoubleEG_LooseIso23_10, &b_L1_DoubleEG_LooseIso23_10);
   chain->SetBranchAddress("L1_DoubleEG_LooseIso24_10", &L1_DoubleEG_LooseIso24_10, &b_L1_DoubleEG_LooseIso24_10);
   chain->SetBranchAddress("L1_DoubleIsoTau28er2p1", &L1_DoubleIsoTau28er2p1, &b_L1_DoubleIsoTau28er2p1);
   chain->SetBranchAddress("L1_DoubleIsoTau30er2p1", &L1_DoubleIsoTau30er2p1, &b_L1_DoubleIsoTau30er2p1);
   chain->SetBranchAddress("L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1, &b_L1_DoubleIsoTau32er2p1);
   chain->SetBranchAddress("L1_DoubleIsoTau33er2p1", &L1_DoubleIsoTau33er2p1, &b_L1_DoubleIsoTau33er2p1);
   chain->SetBranchAddress("L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1, &b_L1_DoubleIsoTau34er2p1);
   chain->SetBranchAddress("L1_DoubleIsoTau35er2p1", &L1_DoubleIsoTau35er2p1, &b_L1_DoubleIsoTau35er2p1);
   chain->SetBranchAddress("L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1, &b_L1_DoubleIsoTau36er2p1);
   chain->SetBranchAddress("L1_DoubleIsoTau38er2p1", &L1_DoubleIsoTau38er2p1, &b_L1_DoubleIsoTau38er2p1);
   chain->SetBranchAddress("L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6, &b_L1_DoubleJet100er2p3_dEta_Max1p6);
   chain->SetBranchAddress("L1_DoubleJet100er2p7", &L1_DoubleJet100er2p7, &b_L1_DoubleJet100er2p7);
   chain->SetBranchAddress("L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6, &b_L1_DoubleJet112er2p3_dEta_Max1p6);
   chain->SetBranchAddress("L1_DoubleJet112er2p7", &L1_DoubleJet112er2p7, &b_L1_DoubleJet112er2p7);
   chain->SetBranchAddress("L1_DoubleJet120er2p7", &L1_DoubleJet120er2p7, &b_L1_DoubleJet120er2p7);
   chain->SetBranchAddress("L1_DoubleJet150er2p7", &L1_DoubleJet150er2p7, &b_L1_DoubleJet150er2p7);
   chain->SetBranchAddress("L1_DoubleJet30_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30_Mass_Min300_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min300_dEta_Max1p5);
   chain->SetBranchAddress("L1_DoubleJet30_Mass_Min320_dEta_Max1p5", &L1_DoubleJet30_Mass_Min320_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min320_dEta_Max1p5);
   chain->SetBranchAddress("L1_DoubleJet30_Mass_Min340_dEta_Max1p5", &L1_DoubleJet30_Mass_Min340_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min340_dEta_Max1p5);
   chain->SetBranchAddress("L1_DoubleJet30_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30_Mass_Min360_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min360_dEta_Max1p5);
   chain->SetBranchAddress("L1_DoubleJet30_Mass_Min380_dEta_Max1p5", &L1_DoubleJet30_Mass_Min380_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min380_dEta_Max1p5);
   chain->SetBranchAddress("L1_DoubleJet30_Mass_Min400_Mu10", &L1_DoubleJet30_Mass_Min400_Mu10, &b_L1_DoubleJet30_Mass_Min400_Mu10);
   chain->SetBranchAddress("L1_DoubleJet30_Mass_Min400_Mu6", &L1_DoubleJet30_Mass_Min400_Mu6, &b_L1_DoubleJet30_Mass_Min400_Mu6);
   chain->SetBranchAddress("L1_DoubleJet30_Mass_Min400_dEta_Max1p5", &L1_DoubleJet30_Mass_Min400_dEta_Max1p5, &b_L1_DoubleJet30_Mass_Min400_dEta_Max1p5);
   chain->SetBranchAddress("L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450", &L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450, &b_L1_DoubleJet35_rmovlp_IsoTau45_Mass_Min450);
   chain->SetBranchAddress("L1_DoubleJet40er2p7", &L1_DoubleJet40er2p7, &b_L1_DoubleJet40er2p7);
   chain->SetBranchAddress("L1_DoubleJet50er2p7", &L1_DoubleJet50er2p7, &b_L1_DoubleJet50er2p7);
   chain->SetBranchAddress("L1_DoubleJet60er2p7", &L1_DoubleJet60er2p7, &b_L1_DoubleJet60er2p7);
   chain->SetBranchAddress("L1_DoubleJet60er2p7_ETM100", &L1_DoubleJet60er2p7_ETM100, &b_L1_DoubleJet60er2p7_ETM100);
   chain->SetBranchAddress("L1_DoubleJet60er2p7_ETM60", &L1_DoubleJet60er2p7_ETM60, &b_L1_DoubleJet60er2p7_ETM60);
   chain->SetBranchAddress("L1_DoubleJet60er2p7_ETM70", &L1_DoubleJet60er2p7_ETM70, &b_L1_DoubleJet60er2p7_ETM70);
   chain->SetBranchAddress("L1_DoubleJet60er2p7_ETM80", &L1_DoubleJet60er2p7_ETM80, &b_L1_DoubleJet60er2p7_ETM80);
   chain->SetBranchAddress("L1_DoubleJet60er2p7_ETM90", &L1_DoubleJet60er2p7_ETM90, &b_L1_DoubleJet60er2p7_ETM90);
   chain->SetBranchAddress("L1_DoubleJet80er2p7", &L1_DoubleJet80er2p7, &b_L1_DoubleJet80er2p7);
   chain->SetBranchAddress("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620, &b_L1_DoubleJet_100_30_DoubleJet30_Mass_Min620);
   chain->SetBranchAddress("L1_DoubleJet_100_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_100_35_DoubleJet35_Mass_Min620, &b_L1_DoubleJet_100_35_DoubleJet35_Mass_Min620);
   chain->SetBranchAddress("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620, &b_L1_DoubleJet_110_35_DoubleJet35_Mass_Min620);
   chain->SetBranchAddress("L1_DoubleJet_110_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_110_40_DoubleJet40_Mass_Min620, &b_L1_DoubleJet_110_40_DoubleJet40_Mass_Min620);
   chain->SetBranchAddress("L1_DoubleJet_115_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_115_35_DoubleJet35_Mass_Min620, &b_L1_DoubleJet_115_35_DoubleJet35_Mass_Min620);
   chain->SetBranchAddress("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620, &b_L1_DoubleJet_115_40_DoubleJet40_Mass_Min620);
   chain->SetBranchAddress("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620, &b_L1_DoubleJet_90_30_DoubleJet30_Mass_Min620);
   chain->SetBranchAddress("L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1, &b_L1_DoubleLooseIsoEG22er2p1);
   chain->SetBranchAddress("L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1, &b_L1_DoubleLooseIsoEG24er2p1);
   chain->SetBranchAddress("L1_DoubleMu0", &L1_DoubleMu0, &b_L1_DoubleMu0);
   chain->SetBranchAddress("L1_DoubleMu0_ETM40", &L1_DoubleMu0_ETM40, &b_L1_DoubleMu0_ETM40);
   chain->SetBranchAddress("L1_DoubleMu0_ETM55", &L1_DoubleMu0_ETM55, &b_L1_DoubleMu0_ETM55);
   chain->SetBranchAddress("L1_DoubleMu0_ETM60", &L1_DoubleMu0_ETM60, &b_L1_DoubleMu0_ETM60);
   chain->SetBranchAddress("L1_DoubleMu0_ETM65", &L1_DoubleMu0_ETM65, &b_L1_DoubleMu0_ETM65);
   chain->SetBranchAddress("L1_DoubleMu0_ETM70", &L1_DoubleMu0_ETM70, &b_L1_DoubleMu0_ETM70);
   chain->SetBranchAddress("L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ, &b_L1_DoubleMu0_SQ);
   chain->SetBranchAddress("L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS, &b_L1_DoubleMu0_SQ_OS);
   chain->SetBranchAddress("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4);
   chain->SetBranchAddress("L1_DoubleMu0er1p4_dEta_Max1p8_OS", &L1_DoubleMu0er1p4_dEta_Max1p8_OS, &b_L1_DoubleMu0er1p4_dEta_Max1p8_OS);
   chain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS, &b_L1_DoubleMu0er1p5_SQ_OS);
   chain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4, &b_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4);
   chain->SetBranchAddress("L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4, &b_L1_DoubleMu0er1p5_SQ_dR_Max1p4);
   chain->SetBranchAddress("L1_DoubleMu0er2_SQ_dR_Max1p4", &L1_DoubleMu0er2_SQ_dR_Max1p4, &b_L1_DoubleMu0er2_SQ_dR_Max1p4);
   chain->SetBranchAddress("L1_DoubleMu18er2p1", &L1_DoubleMu18er2p1, &b_L1_DoubleMu18er2p1);
   chain->SetBranchAddress("L1_DoubleMu22er2p1", &L1_DoubleMu22er2p1, &b_L1_DoubleMu22er2p1);
   chain->SetBranchAddress("L1_DoubleMu3_OS_DoubleEG7p5Upsilon", &L1_DoubleMu3_OS_DoubleEG7p5Upsilon, &b_L1_DoubleMu3_OS_DoubleEG7p5Upsilon);
   chain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF40_Jet60_OR_DoubleJet30);
   chain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF50_Jet60_OR_DoubleJet30);
   chain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF60_Jet60_OR_DoubleJet30);
   chain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF70_Jet60_OR_DoubleJet30);
   chain->SetBranchAddress("L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30", &L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30, &b_L1_DoubleMu3_SQ_ETMHF80_Jet60_OR_DoubleJet30);
   chain->SetBranchAddress("L1_DoubleMu3_SQ_HTT100er", &L1_DoubleMu3_SQ_HTT100er, &b_L1_DoubleMu3_SQ_HTT100er);
   chain->SetBranchAddress("L1_DoubleMu3_SQ_HTT200er", &L1_DoubleMu3_SQ_HTT200er, &b_L1_DoubleMu3_SQ_HTT200er);
   chain->SetBranchAddress("L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er, &b_L1_DoubleMu3_SQ_HTT220er);
   chain->SetBranchAddress("L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er, &b_L1_DoubleMu3_SQ_HTT240er);
   chain->SetBranchAddress("L1_DoubleMu4_OS_EG12", &L1_DoubleMu4_OS_EG12, &b_L1_DoubleMu4_OS_EG12);
   chain->SetBranchAddress("L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS, &b_L1_DoubleMu4_SQ_OS);
   chain->SetBranchAddress("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2, &b_L1_DoubleMu4_SQ_OS_dR_Max1p2);
   chain->SetBranchAddress("L1_DoubleMu4p5_SQ", &L1_DoubleMu4p5_SQ, &b_L1_DoubleMu4p5_SQ);
   chain->SetBranchAddress("L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS, &b_L1_DoubleMu4p5_SQ_OS);
   chain->SetBranchAddress("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2, &b_L1_DoubleMu4p5_SQ_OS_dR_Max1p2);
   chain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS, &b_L1_DoubleMu4p5er2p0_SQ_OS);
   chain->SetBranchAddress("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18, &b_L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18);
   chain->SetBranchAddress("L1_DoubleMu5Upsilon_OS_DoubleEG3", &L1_DoubleMu5Upsilon_OS_DoubleEG3, &b_L1_DoubleMu5Upsilon_OS_DoubleEG3);
   chain->SetBranchAddress("L1_DoubleMu5_OS_EG12", &L1_DoubleMu5_OS_EG12, &b_L1_DoubleMu5_OS_EG12);
   chain->SetBranchAddress("L1_DoubleMu5_SQ_OS", &L1_DoubleMu5_SQ_OS, &b_L1_DoubleMu5_SQ_OS);
   chain->SetBranchAddress("L1_DoubleMu5_SQ_OS_Mass7to18", &L1_DoubleMu5_SQ_OS_Mass7to18, &b_L1_DoubleMu5_SQ_OS_Mass7to18);
   chain->SetBranchAddress("L1_DoubleMu6_SQ_OS", &L1_DoubleMu6_SQ_OS, &b_L1_DoubleMu6_SQ_OS);
   chain->SetBranchAddress("L1_DoubleMu7_EG7", &L1_DoubleMu7_EG7, &b_L1_DoubleMu7_EG7);
   chain->SetBranchAddress("L1_DoubleMu7_SQ_EG7", &L1_DoubleMu7_SQ_EG7, &b_L1_DoubleMu7_SQ_EG7);
   chain->SetBranchAddress("L1_DoubleMu8_SQ", &L1_DoubleMu8_SQ, &b_L1_DoubleMu8_SQ);
   chain->SetBranchAddress("L1_DoubleMu_10_0_dEta_Max1p8", &L1_DoubleMu_10_0_dEta_Max1p8, &b_L1_DoubleMu_10_0_dEta_Max1p8);
   chain->SetBranchAddress("L1_DoubleMu_11_4", &L1_DoubleMu_11_4, &b_L1_DoubleMu_11_4);
   chain->SetBranchAddress("L1_DoubleMu_12_5", &L1_DoubleMu_12_5, &b_L1_DoubleMu_12_5);
   chain->SetBranchAddress("L1_DoubleMu_12_8", &L1_DoubleMu_12_8, &b_L1_DoubleMu_12_8);
   chain->SetBranchAddress("L1_DoubleMu_13_6", &L1_DoubleMu_13_6, &b_L1_DoubleMu_13_6);
   chain->SetBranchAddress("L1_DoubleMu_15_5", &L1_DoubleMu_15_5, &b_L1_DoubleMu_15_5);
   chain->SetBranchAddress("L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ, &b_L1_DoubleMu_15_5_SQ);
   chain->SetBranchAddress("L1_DoubleMu_15_7", &L1_DoubleMu_15_7, &b_L1_DoubleMu_15_7);
   chain->SetBranchAddress("L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ, &b_L1_DoubleMu_15_7_SQ);
   chain->SetBranchAddress("L1_DoubleMu_15_7_SQ_Mass_Min4", &L1_DoubleMu_15_7_SQ_Mass_Min4, &b_L1_DoubleMu_15_7_SQ_Mass_Min4);
   chain->SetBranchAddress("L1_DoubleMu_20_2_SQ_Mass_Max20", &L1_DoubleMu_20_2_SQ_Mass_Max20, &b_L1_DoubleMu_20_2_SQ_Mass_Max20);
   chain->SetBranchAddress("L1_DoubleTau50er2p1", &L1_DoubleTau50er2p1, &b_L1_DoubleTau50er2p1);
   chain->SetBranchAddress("L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1, &b_L1_DoubleTau70er2p1);
   chain->SetBranchAddress("L1_EG25er2p1_HTT125er", &L1_EG25er2p1_HTT125er, &b_L1_EG25er2p1_HTT125er);
   chain->SetBranchAddress("L1_EG27er2p1_HTT200er", &L1_EG27er2p1_HTT200er, &b_L1_EG27er2p1_HTT200er);
   chain->SetBranchAddress("L1_ETM100", &L1_ETM100, &b_L1_ETM100);
   chain->SetBranchAddress("L1_ETM100_Jet60_dPhi_Min0p4", &L1_ETM100_Jet60_dPhi_Min0p4, &b_L1_ETM100_Jet60_dPhi_Min0p4);
   chain->SetBranchAddress("L1_ETM105", &L1_ETM105, &b_L1_ETM105);
   chain->SetBranchAddress("L1_ETM110", &L1_ETM110, &b_L1_ETM110);
   chain->SetBranchAddress("L1_ETM110_Jet60_dPhi_Min0p4", &L1_ETM110_Jet60_dPhi_Min0p4, &b_L1_ETM110_Jet60_dPhi_Min0p4);
   chain->SetBranchAddress("L1_ETM115", &L1_ETM115, &b_L1_ETM115);
   chain->SetBranchAddress("L1_ETM120", &L1_ETM120, &b_L1_ETM120);
   chain->SetBranchAddress("L1_ETM150", &L1_ETM150, &b_L1_ETM150);
   chain->SetBranchAddress("L1_ETM30", &L1_ETM30, &b_L1_ETM30);
   chain->SetBranchAddress("L1_ETM40", &L1_ETM40, &b_L1_ETM40);
   chain->SetBranchAddress("L1_ETM50", &L1_ETM50, &b_L1_ETM50);
   chain->SetBranchAddress("L1_ETM60", &L1_ETM60, &b_L1_ETM60);
   chain->SetBranchAddress("L1_ETM70", &L1_ETM70, &b_L1_ETM70);
   chain->SetBranchAddress("L1_ETM75", &L1_ETM75, &b_L1_ETM75);
   chain->SetBranchAddress("L1_ETM75_Jet60_dPhi_Min0p4", &L1_ETM75_Jet60_dPhi_Min0p4, &b_L1_ETM75_Jet60_dPhi_Min0p4);
   chain->SetBranchAddress("L1_ETM80", &L1_ETM80, &b_L1_ETM80);
   chain->SetBranchAddress("L1_ETM80_Jet60_dPhi_Min0p4", &L1_ETM80_Jet60_dPhi_Min0p4, &b_L1_ETM80_Jet60_dPhi_Min0p4);
   chain->SetBranchAddress("L1_ETM85", &L1_ETM85, &b_L1_ETM85);
   chain->SetBranchAddress("L1_ETM90", &L1_ETM90, &b_L1_ETM90);
   chain->SetBranchAddress("L1_ETM90_Jet60_dPhi_Min0p4", &L1_ETM90_Jet60_dPhi_Min0p4, &b_L1_ETM90_Jet60_dPhi_Min0p4);
   chain->SetBranchAddress("L1_ETM95", &L1_ETM95, &b_L1_ETM95);
   chain->SetBranchAddress("L1_ETMHF100", &L1_ETMHF100, &b_L1_ETMHF100);
   chain->SetBranchAddress("L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er, &b_L1_ETMHF100_HTT60er);
   chain->SetBranchAddress("L1_ETMHF100_Jet60_OR_DiJet30woTT28", &L1_ETMHF100_Jet60_OR_DiJet30woTT28, &b_L1_ETMHF100_Jet60_OR_DiJet30woTT28);
   chain->SetBranchAddress("L1_ETMHF100_Jet60_OR_DoubleJet30", &L1_ETMHF100_Jet60_OR_DoubleJet30, &b_L1_ETMHF100_Jet60_OR_DoubleJet30);
   chain->SetBranchAddress("L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF100_Jet90_OR_DoubleJet45_OR_TripleJet30);
   chain->SetBranchAddress("L1_ETMHF110", &L1_ETMHF110, &b_L1_ETMHF110);
   chain->SetBranchAddress("L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er, &b_L1_ETMHF110_HTT60er);
   chain->SetBranchAddress("L1_ETMHF110_Jet60_OR_DiJet30woTT28", &L1_ETMHF110_Jet60_OR_DiJet30woTT28, &b_L1_ETMHF110_Jet60_OR_DiJet30woTT28);
   chain->SetBranchAddress("L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF110_Jet90_OR_DoubleJet45_OR_TripleJet30);
   chain->SetBranchAddress("L1_ETMHF120", &L1_ETMHF120, &b_L1_ETMHF120);
   chain->SetBranchAddress("L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er, &b_L1_ETMHF120_HTT60er);
   chain->SetBranchAddress("L1_ETMHF120_Jet60_OR_DiJet30woTT28", &L1_ETMHF120_Jet60_OR_DiJet30woTT28, &b_L1_ETMHF120_Jet60_OR_DiJet30woTT28);
   chain->SetBranchAddress("L1_ETMHF150", &L1_ETMHF150, &b_L1_ETMHF150);
   chain->SetBranchAddress("L1_ETMHF70", &L1_ETMHF70, &b_L1_ETMHF70);
   chain->SetBranchAddress("L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF70_Jet90_OR_DoubleJet45_OR_TripleJet30);
   chain->SetBranchAddress("L1_ETMHF80", &L1_ETMHF80, &b_L1_ETMHF80);
   chain->SetBranchAddress("L1_ETMHF80_HTT60er", &L1_ETMHF80_HTT60er, &b_L1_ETMHF80_HTT60er);
   chain->SetBranchAddress("L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF80_Jet90_OR_DoubleJet45_OR_TripleJet30);
   chain->SetBranchAddress("L1_ETMHF90", &L1_ETMHF90, &b_L1_ETMHF90);
   chain->SetBranchAddress("L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er, &b_L1_ETMHF90_HTT60er);
   chain->SetBranchAddress("L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30", &L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30, &b_L1_ETMHF90_Jet90_OR_DoubleJet45_OR_TripleJet30);
   chain->SetBranchAddress("L1_ETT100_BptxAND", &L1_ETT100_BptxAND, &b_L1_ETT100_BptxAND);
   chain->SetBranchAddress("L1_ETT110_BptxAND", &L1_ETT110_BptxAND, &b_L1_ETT110_BptxAND);
   chain->SetBranchAddress("L1_ETT40_BptxAND", &L1_ETT40_BptxAND, &b_L1_ETT40_BptxAND);
   chain->SetBranchAddress("L1_ETT50_BptxAND", &L1_ETT50_BptxAND, &b_L1_ETT50_BptxAND);
   chain->SetBranchAddress("L1_ETT60_BptxAND", &L1_ETT60_BptxAND, &b_L1_ETT60_BptxAND);
   chain->SetBranchAddress("L1_ETT70_BptxAND", &L1_ETT70_BptxAND, &b_L1_ETT70_BptxAND);
   chain->SetBranchAddress("L1_ETT75_BptxAND", &L1_ETT75_BptxAND, &b_L1_ETT75_BptxAND);
   chain->SetBranchAddress("L1_ETT80_BptxAND", &L1_ETT80_BptxAND, &b_L1_ETT80_BptxAND);
   chain->SetBranchAddress("L1_ETT85_BptxAND", &L1_ETT85_BptxAND, &b_L1_ETT85_BptxAND);
   chain->SetBranchAddress("L1_ETT90_BptxAND", &L1_ETT90_BptxAND, &b_L1_ETT90_BptxAND);
   chain->SetBranchAddress("L1_ETT95_BptxAND", &L1_ETT95_BptxAND, &b_L1_ETT95_BptxAND);
   chain->SetBranchAddress("L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain, &b_L1_FirstBunchAfterTrain);
   chain->SetBranchAddress("L1_FirstBunchInTrain", &L1_FirstBunchInTrain, &b_L1_FirstBunchInTrain);
   chain->SetBranchAddress("L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit, &b_L1_FirstCollisionInOrbit);
   chain->SetBranchAddress("L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain, &b_L1_FirstCollisionInTrain);
   chain->SetBranchAddress("L1_HTT120er", &L1_HTT120er, &b_L1_HTT120er);
   chain->SetBranchAddress("L1_HTT160er", &L1_HTT160er, &b_L1_HTT160er);
   chain->SetBranchAddress("L1_HTT200er", &L1_HTT200er, &b_L1_HTT200er);
   chain->SetBranchAddress("L1_HTT220er", &L1_HTT220er, &b_L1_HTT220er);
   chain->SetBranchAddress("L1_HTT240er", &L1_HTT240er, &b_L1_HTT240er);
   chain->SetBranchAddress("L1_HTT250er_QuadJet_70_55_40_35_er2p5", &L1_HTT250er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT250er_QuadJet_70_55_40_35_er2p5);
   chain->SetBranchAddress("L1_HTT255er", &L1_HTT255er, &b_L1_HTT255er);
   chain->SetBranchAddress("L1_HTT270er", &L1_HTT270er, &b_L1_HTT270er);
   chain->SetBranchAddress("L1_HTT280er", &L1_HTT280er, &b_L1_HTT280er);
   chain->SetBranchAddress("L1_HTT280er_QuadJet_70_55_40_35_er2p5", &L1_HTT280er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT280er_QuadJet_70_55_40_35_er2p5);
   chain->SetBranchAddress("L1_HTT300er", &L1_HTT300er, &b_L1_HTT300er);
   chain->SetBranchAddress("L1_HTT300er_QuadJet_70_55_40_35_er2p5", &L1_HTT300er_QuadJet_70_55_40_35_er2p5, &b_L1_HTT300er_QuadJet_70_55_40_35_er2p5);
   chain->SetBranchAddress("L1_HTT320er", &L1_HTT320er, &b_L1_HTT320er);
   chain->SetBranchAddress("L1_HTT320er_QuadJet_70_55_40_40_er2p4", &L1_HTT320er_QuadJet_70_55_40_40_er2p4, &b_L1_HTT320er_QuadJet_70_55_40_40_er2p4);
   chain->SetBranchAddress("L1_HTT320er_QuadJet_70_55_40_40_er2p5", &L1_HTT320er_QuadJet_70_55_40_40_er2p5, &b_L1_HTT320er_QuadJet_70_55_40_40_er2p5);
   chain->SetBranchAddress("L1_HTT320er_QuadJet_70_55_45_45_er2p5", &L1_HTT320er_QuadJet_70_55_45_45_er2p5, &b_L1_HTT320er_QuadJet_70_55_45_45_er2p5);
   chain->SetBranchAddress("L1_HTT340er", &L1_HTT340er, &b_L1_HTT340er);
   chain->SetBranchAddress("L1_HTT340er_QuadJet_70_55_40_40_er2p5", &L1_HTT340er_QuadJet_70_55_40_40_er2p5, &b_L1_HTT340er_QuadJet_70_55_40_40_er2p5);
   chain->SetBranchAddress("L1_HTT340er_QuadJet_70_55_45_45_er2p5", &L1_HTT340er_QuadJet_70_55_45_45_er2p5, &b_L1_HTT340er_QuadJet_70_55_45_45_er2p5);
   chain->SetBranchAddress("L1_HTT380er", &L1_HTT380er, &b_L1_HTT380er);
   chain->SetBranchAddress("L1_HTT400er", &L1_HTT400er, &b_L1_HTT400er);
   chain->SetBranchAddress("L1_HTT450er", &L1_HTT450er, &b_L1_HTT450er);
   chain->SetBranchAddress("L1_HTT500er", &L1_HTT500er, &b_L1_HTT500er);
   chain->SetBranchAddress("L1_IsoEG33_Mt40", &L1_IsoEG33_Mt40, &b_L1_IsoEG33_Mt40);
   chain->SetBranchAddress("L1_IsoEG33_Mt44", &L1_IsoEG33_Mt44, &b_L1_IsoEG33_Mt44);
   chain->SetBranchAddress("L1_IsoEG33_Mt48", &L1_IsoEG33_Mt48, &b_L1_IsoEG33_Mt48);
   chain->SetBranchAddress("L1_IsoTau40er_ETM100", &L1_IsoTau40er_ETM100, &b_L1_IsoTau40er_ETM100);
   chain->SetBranchAddress("L1_IsoTau40er_ETM105", &L1_IsoTau40er_ETM105, &b_L1_IsoTau40er_ETM105);
   chain->SetBranchAddress("L1_IsoTau40er_ETM110", &L1_IsoTau40er_ETM110, &b_L1_IsoTau40er_ETM110);
   chain->SetBranchAddress("L1_IsoTau40er_ETM115", &L1_IsoTau40er_ETM115, &b_L1_IsoTau40er_ETM115);
   chain->SetBranchAddress("L1_IsoTau40er_ETM120", &L1_IsoTau40er_ETM120, &b_L1_IsoTau40er_ETM120);
   chain->SetBranchAddress("L1_IsoTau40er_ETM80", &L1_IsoTau40er_ETM80, &b_L1_IsoTau40er_ETM80);
   chain->SetBranchAddress("L1_IsoTau40er_ETM85", &L1_IsoTau40er_ETM85, &b_L1_IsoTau40er_ETM85);
   chain->SetBranchAddress("L1_IsoTau40er_ETM90", &L1_IsoTau40er_ETM90, &b_L1_IsoTau40er_ETM90);
   chain->SetBranchAddress("L1_IsoTau40er_ETM95", &L1_IsoTau40er_ETM95, &b_L1_IsoTau40er_ETM95);
   chain->SetBranchAddress("L1_IsoTau40er_ETMHF100", &L1_IsoTau40er_ETMHF100, &b_L1_IsoTau40er_ETMHF100);
   chain->SetBranchAddress("L1_IsoTau40er_ETMHF110", &L1_IsoTau40er_ETMHF110, &b_L1_IsoTau40er_ETMHF110);
   chain->SetBranchAddress("L1_IsoTau40er_ETMHF120", &L1_IsoTau40er_ETMHF120, &b_L1_IsoTau40er_ETMHF120);
   chain->SetBranchAddress("L1_IsoTau40er_ETMHF80", &L1_IsoTau40er_ETMHF80, &b_L1_IsoTau40er_ETMHF80);
   chain->SetBranchAddress("L1_IsoTau40er_ETMHF90", &L1_IsoTau40er_ETMHF90, &b_L1_IsoTau40er_ETMHF90);
   chain->SetBranchAddress("L1_IsolatedBunch", &L1_IsolatedBunch, &b_L1_IsolatedBunch);
   chain->SetBranchAddress("L1_LastCollisionInTrain", &L1_LastCollisionInTrain, &b_L1_LastCollisionInTrain);
   chain->SetBranchAddress("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3, &b_L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3);
   chain->SetBranchAddress("L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er, &b_L1_LooseIsoEG24er2p1_HTT100er);
   chain->SetBranchAddress("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3, &b_L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3);
   chain->SetBranchAddress("L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3", &L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3, &b_L1_LooseIsoEG24er2p1_Jet26er2p7_dR_Min0p3);
   chain->SetBranchAddress("L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7", &L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7, &b_L1_LooseIsoEG24er2p1_TripleJet_26er2p7_26_26er2p7);
   chain->SetBranchAddress("L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er, &b_L1_LooseIsoEG26er2p1_HTT100er);
   chain->SetBranchAddress("L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3, &b_L1_LooseIsoEG26er2p1_Jet34er2p7_dR_Min0p3);
   chain->SetBranchAddress("L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er, &b_L1_LooseIsoEG28er2p1_HTT100er);
   chain->SetBranchAddress("L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3, &b_L1_LooseIsoEG28er2p1_Jet34er2p7_dR_Min0p3);
   chain->SetBranchAddress("L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3, &b_L1_LooseIsoEG30er2p1_Jet34er2p7_dR_Min0p3);
   chain->SetBranchAddress("L1_MU20_EG15", &L1_MU20_EG15, &b_L1_MU20_EG15);
   chain->SetBranchAddress("L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND, &b_L1_MinimumBiasHF0_AND_BptxAND);
   chain->SetBranchAddress("L1_MinimumBiasHF0_OR_BptxAND", &L1_MinimumBiasHF0_OR_BptxAND, &b_L1_MinimumBiasHF0_OR_BptxAND);
   chain->SetBranchAddress("L1_Mu10er2p1_ETM30", &L1_Mu10er2p1_ETM30, &b_L1_Mu10er2p1_ETM30);
   chain->SetBranchAddress("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6, &b_L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6);
   chain->SetBranchAddress("L1_Mu12_EG10", &L1_Mu12_EG10, &b_L1_Mu12_EG10);
   chain->SetBranchAddress("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6, &b_L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6);
   chain->SetBranchAddress("L1_Mu14er2p1_ETM30", &L1_Mu14er2p1_ETM30, &b_L1_Mu14er2p1_ETM30);
   chain->SetBranchAddress("L1_Mu15_HTT100er", &L1_Mu15_HTT100er, &b_L1_Mu15_HTT100er);
   chain->SetBranchAddress("L1_Mu18_HTT100er", &L1_Mu18_HTT100er, &b_L1_Mu18_HTT100er);
   chain->SetBranchAddress("L1_Mu18_Jet24er2p7", &L1_Mu18_Jet24er2p7, &b_L1_Mu18_Jet24er2p7);
   chain->SetBranchAddress("L1_Mu18er2p1_IsoTau26er2p1", &L1_Mu18er2p1_IsoTau26er2p1, &b_L1_Mu18er2p1_IsoTau26er2p1);
   chain->SetBranchAddress("L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1, &b_L1_Mu18er2p1_Tau24er2p1);
   chain->SetBranchAddress("L1_Mu20_EG10", &L1_Mu20_EG10, &b_L1_Mu20_EG10);
   chain->SetBranchAddress("L1_Mu20_EG17", &L1_Mu20_EG17, &b_L1_Mu20_EG17);
   chain->SetBranchAddress("L1_Mu20_LooseIsoEG6", &L1_Mu20_LooseIsoEG6, &b_L1_Mu20_LooseIsoEG6);
   chain->SetBranchAddress("L1_Mu20er2p1_IsoTau26er2p1", &L1_Mu20er2p1_IsoTau26er2p1, &b_L1_Mu20er2p1_IsoTau26er2p1);
   chain->SetBranchAddress("L1_Mu20er2p1_IsoTau27er2p1", &L1_Mu20er2p1_IsoTau27er2p1, &b_L1_Mu20er2p1_IsoTau27er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_IsoTau28er2p1", &L1_Mu22er2p1_IsoTau28er2p1, &b_L1_Mu22er2p1_IsoTau28er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_IsoTau30er2p1", &L1_Mu22er2p1_IsoTau30er2p1, &b_L1_Mu22er2p1_IsoTau30er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1, &b_L1_Mu22er2p1_IsoTau32er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_IsoTau33er2p1", &L1_Mu22er2p1_IsoTau33er2p1, &b_L1_Mu22er2p1_IsoTau33er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1, &b_L1_Mu22er2p1_IsoTau34er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_IsoTau35er2p1", &L1_Mu22er2p1_IsoTau35er2p1, &b_L1_Mu22er2p1_IsoTau35er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1, &b_L1_Mu22er2p1_IsoTau36er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_IsoTau38er2p1", &L1_Mu22er2p1_IsoTau38er2p1, &b_L1_Mu22er2p1_IsoTau38er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1, &b_L1_Mu22er2p1_IsoTau40er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_Tau50er2p1", &L1_Mu22er2p1_Tau50er2p1, &b_L1_Mu22er2p1_Tau50er2p1);
   chain->SetBranchAddress("L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1, &b_L1_Mu22er2p1_Tau70er2p1);
   chain->SetBranchAddress("L1_Mu23_EG10", &L1_Mu23_EG10, &b_L1_Mu23_EG10);
   chain->SetBranchAddress("L1_Mu23_LooseIsoEG10", &L1_Mu23_LooseIsoEG10, &b_L1_Mu23_LooseIsoEG10);
   chain->SetBranchAddress("L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4", &L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4, &b_L1_Mu3_Jet120er2p7_dEta_Max0p4_dPhi_Max0p4);
   chain->SetBranchAddress("L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4", &L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4, &b_L1_Mu3_Jet16er2p7_dEta_Max0p4_dPhi_Max0p4);
   chain->SetBranchAddress("L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5, &b_L1_Mu3_Jet30er2p5);
   chain->SetBranchAddress("L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4", &L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4, &b_L1_Mu3_Jet60er2p7_dEta_Max0p4_dPhi_Max0p4);
   chain->SetBranchAddress("L1_Mu5_EG15", &L1_Mu5_EG15, &b_L1_Mu5_EG15);
   chain->SetBranchAddress("L1_Mu5_EG20", &L1_Mu5_EG20, &b_L1_Mu5_EG20);
   chain->SetBranchAddress("L1_Mu5_EG23", &L1_Mu5_EG23, &b_L1_Mu5_EG23);
   chain->SetBranchAddress("L1_Mu5_LooseIsoEG18", &L1_Mu5_LooseIsoEG18, &b_L1_Mu5_LooseIsoEG18);
   chain->SetBranchAddress("L1_Mu5_LooseIsoEG20", &L1_Mu5_LooseIsoEG20, &b_L1_Mu5_LooseIsoEG20);
   chain->SetBranchAddress("L1_Mu6_DoubleEG10", &L1_Mu6_DoubleEG10, &b_L1_Mu6_DoubleEG10);
   chain->SetBranchAddress("L1_Mu6_DoubleEG17", &L1_Mu6_DoubleEG17, &b_L1_Mu6_DoubleEG17);
   chain->SetBranchAddress("L1_Mu6_HTT200er", &L1_Mu6_HTT200er, &b_L1_Mu6_HTT200er);
   chain->SetBranchAddress("L1_Mu6_HTT240er", &L1_Mu6_HTT240er, &b_L1_Mu6_HTT240er);
   chain->SetBranchAddress("L1_Mu6_HTT250er", &L1_Mu6_HTT250er, &b_L1_Mu6_HTT250er);
   chain->SetBranchAddress("L1_Mu7_EG23", &L1_Mu7_EG23, &b_L1_Mu7_EG23);
   chain->SetBranchAddress("L1_Mu7_LooseIsoEG20", &L1_Mu7_LooseIsoEG20, &b_L1_Mu7_LooseIsoEG20);
   chain->SetBranchAddress("L1_Mu7_LooseIsoEG23", &L1_Mu7_LooseIsoEG23, &b_L1_Mu7_LooseIsoEG23);
   chain->SetBranchAddress("L1_Mu8_HTT150er", &L1_Mu8_HTT150er, &b_L1_Mu8_HTT150er);
   chain->SetBranchAddress("L1_NotBptxOR", &L1_NotBptxOR, &b_L1_NotBptxOR);
   chain->SetBranchAddress("L1_QuadJet36er2p7_IsoTau52er2p1", &L1_QuadJet36er2p7_IsoTau52er2p1, &b_L1_QuadJet36er2p7_IsoTau52er2p1);
   chain->SetBranchAddress("L1_QuadJet36er2p7_Tau52", &L1_QuadJet36er2p7_Tau52, &b_L1_QuadJet36er2p7_Tau52);
   chain->SetBranchAddress("L1_QuadJet40er2p7", &L1_QuadJet40er2p7, &b_L1_QuadJet40er2p7);
   chain->SetBranchAddress("L1_QuadJet50er2p7", &L1_QuadJet50er2p7, &b_L1_QuadJet50er2p7);
   chain->SetBranchAddress("L1_QuadJet60er2p7", &L1_QuadJet60er2p7, &b_L1_QuadJet60er2p7);
   chain->SetBranchAddress("L1_QuadMu0", &L1_QuadMu0, &b_L1_QuadMu0);
   chain->SetBranchAddress("L1_SingleEG10", &L1_SingleEG10, &b_L1_SingleEG10);
   chain->SetBranchAddress("L1_SingleEG15", &L1_SingleEG15, &b_L1_SingleEG15);
   chain->SetBranchAddress("L1_SingleEG18", &L1_SingleEG18, &b_L1_SingleEG18);
   chain->SetBranchAddress("L1_SingleEG24", &L1_SingleEG24, &b_L1_SingleEG24);
   chain->SetBranchAddress("L1_SingleEG26", &L1_SingleEG26, &b_L1_SingleEG26);
   chain->SetBranchAddress("L1_SingleEG28", &L1_SingleEG28, &b_L1_SingleEG28);
   chain->SetBranchAddress("L1_SingleEG2_BptxAND", &L1_SingleEG2_BptxAND, &b_L1_SingleEG2_BptxAND);
   chain->SetBranchAddress("L1_SingleEG30", &L1_SingleEG30, &b_L1_SingleEG30);
   chain->SetBranchAddress("L1_SingleEG32", &L1_SingleEG32, &b_L1_SingleEG32);
   chain->SetBranchAddress("L1_SingleEG34", &L1_SingleEG34, &b_L1_SingleEG34);
   chain->SetBranchAddress("L1_SingleEG34er2p1", &L1_SingleEG34er2p1, &b_L1_SingleEG34er2p1);
   chain->SetBranchAddress("L1_SingleEG36", &L1_SingleEG36, &b_L1_SingleEG36);
   chain->SetBranchAddress("L1_SingleEG36er2p1", &L1_SingleEG36er2p1, &b_L1_SingleEG36er2p1);
   chain->SetBranchAddress("L1_SingleEG38", &L1_SingleEG38, &b_L1_SingleEG38);
   chain->SetBranchAddress("L1_SingleEG38er2p1", &L1_SingleEG38er2p1, &b_L1_SingleEG38er2p1);
   chain->SetBranchAddress("L1_SingleEG40", &L1_SingleEG40, &b_L1_SingleEG40);
   chain->SetBranchAddress("L1_SingleEG42", &L1_SingleEG42, &b_L1_SingleEG42);
   chain->SetBranchAddress("L1_SingleEG45", &L1_SingleEG45, &b_L1_SingleEG45);
   chain->SetBranchAddress("L1_SingleEG5", &L1_SingleEG5, &b_L1_SingleEG5);
   chain->SetBranchAddress("L1_SingleEG50", &L1_SingleEG50, &b_L1_SingleEG50);
   chain->SetBranchAddress("L1_SingleIsoEG18", &L1_SingleIsoEG18, &b_L1_SingleIsoEG18);
   chain->SetBranchAddress("L1_SingleIsoEG18er2p1", &L1_SingleIsoEG18er2p1, &b_L1_SingleIsoEG18er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG20", &L1_SingleIsoEG20, &b_L1_SingleIsoEG20);
   chain->SetBranchAddress("L1_SingleIsoEG20er2p1", &L1_SingleIsoEG20er2p1, &b_L1_SingleIsoEG20er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG22", &L1_SingleIsoEG22, &b_L1_SingleIsoEG22);
   chain->SetBranchAddress("L1_SingleIsoEG22er2p1", &L1_SingleIsoEG22er2p1, &b_L1_SingleIsoEG22er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG24", &L1_SingleIsoEG24, &b_L1_SingleIsoEG24);
   chain->SetBranchAddress("L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1, &b_L1_SingleIsoEG24er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG26", &L1_SingleIsoEG26, &b_L1_SingleIsoEG26);
   chain->SetBranchAddress("L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1, &b_L1_SingleIsoEG26er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG28", &L1_SingleIsoEG28, &b_L1_SingleIsoEG28);
   chain->SetBranchAddress("L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1, &b_L1_SingleIsoEG28er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG30", &L1_SingleIsoEG30, &b_L1_SingleIsoEG30);
   chain->SetBranchAddress("L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1, &b_L1_SingleIsoEG30er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG32", &L1_SingleIsoEG32, &b_L1_SingleIsoEG32);
   chain->SetBranchAddress("L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1, &b_L1_SingleIsoEG32er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG33er2p1", &L1_SingleIsoEG33er2p1, &b_L1_SingleIsoEG33er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG34", &L1_SingleIsoEG34, &b_L1_SingleIsoEG34);
   chain->SetBranchAddress("L1_SingleIsoEG34er2p1", &L1_SingleIsoEG34er2p1, &b_L1_SingleIsoEG34er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG35", &L1_SingleIsoEG35, &b_L1_SingleIsoEG35);
   chain->SetBranchAddress("L1_SingleIsoEG35er2p1", &L1_SingleIsoEG35er2p1, &b_L1_SingleIsoEG35er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG36", &L1_SingleIsoEG36, &b_L1_SingleIsoEG36);
   chain->SetBranchAddress("L1_SingleIsoEG36er2p1", &L1_SingleIsoEG36er2p1, &b_L1_SingleIsoEG36er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG37", &L1_SingleIsoEG37, &b_L1_SingleIsoEG37);
   chain->SetBranchAddress("L1_SingleIsoEG38", &L1_SingleIsoEG38, &b_L1_SingleIsoEG38);
   chain->SetBranchAddress("L1_SingleIsoEG38er2p1", &L1_SingleIsoEG38er2p1, &b_L1_SingleIsoEG38er2p1);
   chain->SetBranchAddress("L1_SingleIsoEG40", &L1_SingleIsoEG40, &b_L1_SingleIsoEG40);
   chain->SetBranchAddress("L1_SingleIsoEG40er2p1", &L1_SingleIsoEG40er2p1, &b_L1_SingleIsoEG40er2p1);
   chain->SetBranchAddress("L1_SingleJet120", &L1_SingleJet120, &b_L1_SingleJet120);
   chain->SetBranchAddress("L1_SingleJet120_FWD", &L1_SingleJet120_FWD, &b_L1_SingleJet120_FWD);
   chain->SetBranchAddress("L1_SingleJet12_BptxAND", &L1_SingleJet12_BptxAND, &b_L1_SingleJet12_BptxAND);
   chain->SetBranchAddress("L1_SingleJet140", &L1_SingleJet140, &b_L1_SingleJet140);
   chain->SetBranchAddress("L1_SingleJet150", &L1_SingleJet150, &b_L1_SingleJet150);
   chain->SetBranchAddress("L1_SingleJet16", &L1_SingleJet16, &b_L1_SingleJet16);
   chain->SetBranchAddress("L1_SingleJet160", &L1_SingleJet160, &b_L1_SingleJet160);
   chain->SetBranchAddress("L1_SingleJet170", &L1_SingleJet170, &b_L1_SingleJet170);
   chain->SetBranchAddress("L1_SingleJet180", &L1_SingleJet180, &b_L1_SingleJet180);
   chain->SetBranchAddress("L1_SingleJet20", &L1_SingleJet20, &b_L1_SingleJet20);
   chain->SetBranchAddress("L1_SingleJet200", &L1_SingleJet200, &b_L1_SingleJet200);
   chain->SetBranchAddress("L1_SingleJet20er2p7_NotBptxOR", &L1_SingleJet20er2p7_NotBptxOR, &b_L1_SingleJet20er2p7_NotBptxOR);
   chain->SetBranchAddress("L1_SingleJet20er2p7_NotBptxOR_3BX", &L1_SingleJet20er2p7_NotBptxOR_3BX, &b_L1_SingleJet20er2p7_NotBptxOR_3BX);
   chain->SetBranchAddress("L1_SingleJet35", &L1_SingleJet35, &b_L1_SingleJet35);
   chain->SetBranchAddress("L1_SingleJet35_FWD", &L1_SingleJet35_FWD, &b_L1_SingleJet35_FWD);
   chain->SetBranchAddress("L1_SingleJet35_HFm", &L1_SingleJet35_HFm, &b_L1_SingleJet35_HFm);
   chain->SetBranchAddress("L1_SingleJet35_HFp", &L1_SingleJet35_HFp, &b_L1_SingleJet35_HFp);
   chain->SetBranchAddress("L1_SingleJet43er2p7_NotBptxOR_3BX", &L1_SingleJet43er2p7_NotBptxOR_3BX, &b_L1_SingleJet43er2p7_NotBptxOR_3BX);
   chain->SetBranchAddress("L1_SingleJet46er2p7_NotBptxOR_3BX", &L1_SingleJet46er2p7_NotBptxOR_3BX, &b_L1_SingleJet46er2p7_NotBptxOR_3BX);
   chain->SetBranchAddress("L1_SingleJet60", &L1_SingleJet60, &b_L1_SingleJet60);
   chain->SetBranchAddress("L1_SingleJet60_FWD", &L1_SingleJet60_FWD, &b_L1_SingleJet60_FWD);
   chain->SetBranchAddress("L1_SingleJet60_HFm", &L1_SingleJet60_HFm, &b_L1_SingleJet60_HFm);
   chain->SetBranchAddress("L1_SingleJet60_HFp", &L1_SingleJet60_HFp, &b_L1_SingleJet60_HFp);
   chain->SetBranchAddress("L1_SingleJet90", &L1_SingleJet90, &b_L1_SingleJet90);
   chain->SetBranchAddress("L1_SingleJet90_FWD", &L1_SingleJet90_FWD, &b_L1_SingleJet90_FWD);
   chain->SetBranchAddress("L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF, &b_L1_SingleMu0_BMTF);
   chain->SetBranchAddress("L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF, &b_L1_SingleMu0_EMTF);
   chain->SetBranchAddress("L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF, &b_L1_SingleMu0_OMTF);
   chain->SetBranchAddress("L1_SingleMu10_LowQ", &L1_SingleMu10_LowQ, &b_L1_SingleMu10_LowQ);
   chain->SetBranchAddress("L1_SingleMu11_LowQ", &L1_SingleMu11_LowQ, &b_L1_SingleMu11_LowQ);
   chain->SetBranchAddress("L1_SingleMu12_LowQ_BMTF", &L1_SingleMu12_LowQ_BMTF, &b_L1_SingleMu12_LowQ_BMTF);
   chain->SetBranchAddress("L1_SingleMu12_LowQ_EMTF", &L1_SingleMu12_LowQ_EMTF, &b_L1_SingleMu12_LowQ_EMTF);
   chain->SetBranchAddress("L1_SingleMu12_LowQ_OMTF", &L1_SingleMu12_LowQ_OMTF, &b_L1_SingleMu12_LowQ_OMTF);
   chain->SetBranchAddress("L1_SingleMu14er2p1", &L1_SingleMu14er2p1, &b_L1_SingleMu14er2p1);
   chain->SetBranchAddress("L1_SingleMu16", &L1_SingleMu16, &b_L1_SingleMu16);
   chain->SetBranchAddress("L1_SingleMu16er2p1", &L1_SingleMu16er2p1, &b_L1_SingleMu16er2p1);
   chain->SetBranchAddress("L1_SingleMu18", &L1_SingleMu18, &b_L1_SingleMu18);
   chain->SetBranchAddress("L1_SingleMu18er2p1", &L1_SingleMu18er2p1, &b_L1_SingleMu18er2p1);
   chain->SetBranchAddress("L1_SingleMu20", &L1_SingleMu20, &b_L1_SingleMu20);
   chain->SetBranchAddress("L1_SingleMu20er2p1", &L1_SingleMu20er2p1, &b_L1_SingleMu20er2p1);
   chain->SetBranchAddress("L1_SingleMu22", &L1_SingleMu22, &b_L1_SingleMu22);
   chain->SetBranchAddress("L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF, &b_L1_SingleMu22_BMTF);
   chain->SetBranchAddress("L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF, &b_L1_SingleMu22_EMTF);
   chain->SetBranchAddress("L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF, &b_L1_SingleMu22_OMTF);
   chain->SetBranchAddress("L1_SingleMu22er2p1", &L1_SingleMu22er2p1, &b_L1_SingleMu22er2p1);
   chain->SetBranchAddress("L1_SingleMu25", &L1_SingleMu25, &b_L1_SingleMu25);
   chain->SetBranchAddress("L1_SingleMu3", &L1_SingleMu3, &b_L1_SingleMu3);
   chain->SetBranchAddress("L1_SingleMu30", &L1_SingleMu30, &b_L1_SingleMu30);
   chain->SetBranchAddress("L1_SingleMu5", &L1_SingleMu5, &b_L1_SingleMu5);
   chain->SetBranchAddress("L1_SingleMu7", &L1_SingleMu7, &b_L1_SingleMu7);
   chain->SetBranchAddress("L1_SingleMuCosmics", &L1_SingleMuCosmics, &b_L1_SingleMuCosmics);
   chain->SetBranchAddress("L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF, &b_L1_SingleMuCosmics_BMTF);
   chain->SetBranchAddress("L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF, &b_L1_SingleMuCosmics_EMTF);
   chain->SetBranchAddress("L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF, &b_L1_SingleMuCosmics_OMTF);
   chain->SetBranchAddress("L1_SingleMuOpen", &L1_SingleMuOpen, &b_L1_SingleMuOpen);
   chain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR, &b_L1_SingleMuOpen_NotBptxOR);
   chain->SetBranchAddress("L1_SingleMuOpen_NotBptxOR_3BX", &L1_SingleMuOpen_NotBptxOR_3BX, &b_L1_SingleMuOpen_NotBptxOR_3BX);
   chain->SetBranchAddress("L1_SingleTau100er2p1", &L1_SingleTau100er2p1, &b_L1_SingleTau100er2p1);
   chain->SetBranchAddress("L1_SingleTau120er2p1", &L1_SingleTau120er2p1, &b_L1_SingleTau120er2p1);
   chain->SetBranchAddress("L1_SingleTau130er2p1", &L1_SingleTau130er2p1, &b_L1_SingleTau130er2p1);
   chain->SetBranchAddress("L1_SingleTau140er2p1", &L1_SingleTau140er2p1, &b_L1_SingleTau140er2p1);
   chain->SetBranchAddress("L1_SingleTau20", &L1_SingleTau20, &b_L1_SingleTau20);
   chain->SetBranchAddress("L1_SingleTau80er2p1", &L1_SingleTau80er2p1, &b_L1_SingleTau80er2p1);
   chain->SetBranchAddress("L1_TripleEG_14_10_8", &L1_TripleEG_14_10_8, &b_L1_TripleEG_14_10_8);
   chain->SetBranchAddress("L1_TripleEG_18_17_8", &L1_TripleEG_18_17_8, &b_L1_TripleEG_18_17_8);
   chain->SetBranchAddress("L1_TripleEG_LooseIso20_10_5", &L1_TripleEG_LooseIso20_10_5, &b_L1_TripleEG_LooseIso20_10_5);
   chain->SetBranchAddress("L1_TripleJet_100_85_72_VBF", &L1_TripleJet_100_85_72_VBF, &b_L1_TripleJet_100_85_72_VBF);
   chain->SetBranchAddress("L1_TripleJet_105_85_76_VBF", &L1_TripleJet_105_85_76_VBF, &b_L1_TripleJet_105_85_76_VBF);
   chain->SetBranchAddress("L1_TripleJet_84_68_48_VBF", &L1_TripleJet_84_68_48_VBF, &b_L1_TripleJet_84_68_48_VBF);
   chain->SetBranchAddress("L1_TripleJet_88_72_56_VBF", &L1_TripleJet_88_72_56_VBF, &b_L1_TripleJet_88_72_56_VBF);
   chain->SetBranchAddress("L1_TripleJet_92_76_64_VBF", &L1_TripleJet_92_76_64_VBF, &b_L1_TripleJet_92_76_64_VBF);
   chain->SetBranchAddress("L1_TripleJet_98_83_71_VBF", &L1_TripleJet_98_83_71_VBF, &b_L1_TripleJet_98_83_71_VBF);
   chain->SetBranchAddress("L1_TripleMu0", &L1_TripleMu0, &b_L1_TripleMu0);
   chain->SetBranchAddress("L1_TripleMu0_OQ", &L1_TripleMu0_OQ, &b_L1_TripleMu0_OQ);
   chain->SetBranchAddress("L1_TripleMu3", &L1_TripleMu3, &b_L1_TripleMu3);
   chain->SetBranchAddress("L1_TripleMu3_SQ", &L1_TripleMu3_SQ, &b_L1_TripleMu3_SQ);
   chain->SetBranchAddress("L1_TripleMu_4_4_4", &L1_TripleMu_4_4_4, &b_L1_TripleMu_4_4_4);
   chain->SetBranchAddress("L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17, &b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17);
   chain->SetBranchAddress("L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14", &L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14, &b_L1_TripleMu_5OQ_3p5OQ_2p5OQ_DoubleMu_5_2p5_OQ_OS_Mass_8to14);
   chain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ, &b_L1_TripleMu_5SQ_3SQ_0OQ);
   chain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9);
   chain->SetBranchAddress("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9, &b_L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9);
   chain->SetBranchAddress("L1_TripleMu_5_0_0", &L1_TripleMu_5_0_0, &b_L1_TripleMu_5_0_0);
   chain->SetBranchAddress("L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3, &b_L1_TripleMu_5_3_3);
   chain->SetBranchAddress("L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5, &b_L1_TripleMu_5_3p5_2p5);
   chain->SetBranchAddress("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   chain->SetBranchAddress("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17, &b_L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17);
   chain->SetBranchAddress("L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3, &b_L1_TripleMu_5_5_3);
   chain->SetBranchAddress("L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus, &b_L1_UnpairedBunchBptxMinus);
   chain->SetBranchAddress("L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus, &b_L1_UnpairedBunchBptxPlus);
   chain->SetBranchAddress("L1_ZeroBias", &L1_ZeroBias, &b_L1_ZeroBias);
   chain->SetBranchAddress("L1_ZeroBias_copy", &L1_ZeroBias_copy, &b_L1_ZeroBias_copy);
   chain1  = new TChain("Runs");
   chain1->SetMakeClass(1);
   chain1->SetBranchAddress("genEventSumw", &genEventSumw, &b_genEventSumw);
   chain1->SetBranchAddress("genEventSumw_", &genEventSumw_, &b_genEventSumw_);

   Notify();
}

Bool_t ZtoJPsiGamma::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZtoJPsiGamma::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!chain) return;
   chain->Show(entry);
}
Int_t ZtoJPsiGamma::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ZtoJPsiGamma_cxx
