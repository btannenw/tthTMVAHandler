//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 27 06:58:18 2018 by ROOT version 6.10/09
// from TTree tree/tree
// found on file: out.root
//////////////////////////////////////////////////////////

#ifndef base_h
#define base_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class base {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           EventNumber;
   Int_t           EventCategory;
   Float_t         TotalEventWeight;
   Int_t           id_additionalJetEventId;
   Int_t           TruePV;
   Int_t           GenPV;
   Int_t           N_Vtx;
   Float_t         wgt_TOTAL;
   Float_t         wgt_PU;
   Float_t         wgt_btag;
   Float_t         wgt_trigMu;
   Float_t         wgt_trigEl;
   Float_t         wgt_EGZtvx;
   Float_t         wgt_muID;
   Float_t         wgt_elID;
   Float_t         wgt_elReco;
   Float_t         wgt_muIso;
   std::vector<float>   *periodwgt_PU;
   std::vector<float>   *periodwgt_btag;
   std::vector<float>   *periodwgt_trigMu;
   std::vector<float>   *periodwgt_trigEl;
   std::vector<float>   *periodwgt_EGZtvx;
   std::vector<float>   *periodwgt_muID;
   std::vector<float>   *periodwgt_elID;
   std::vector<float>   *periodwgt_elReco;
   std::vector<float>   *periodwgt_muIso;
   Float_t         wgt_trigMuUP;
   Float_t         wgt_trigElUP;
   Float_t         wgt_muIDUP;
   Float_t         wgt_elIDUP;
   Float_t         wgt_elRecoUP;
   Float_t         wgt_muIsoUP;
   Float_t         wgt_trigMuDOWN;
   Float_t         wgt_trigElDOWN;
   Float_t         wgt_muIDDOWN;
   Float_t         wgt_elIDDOWN;
   Float_t         wgt_elRecoDOWN;
   Float_t         wgt_muIsoDOWN;
   Float_t         wgt_MCEventWeight;
   Float_t         wgt_topPT;
   std::vector<double>  *mcWeight_value;
   Int_t           pass_goodVtx;
   Int_t           pass_El;
   Int_t           pass_Mu;
   Int_t           pass_ElEl;
   Int_t           pass_ElMu;
   Int_t           pass_MuMu;
   Int_t           pass_TrigEl;
   Int_t           pass_TrigMu;
   Int_t           pass_TrigElEl;
   Int_t           pass_TrigElMu;
   Int_t           pass_TrigMuMu;
   Int_t           nNonIsoEl;
   Int_t           nNonIsoMu;
   Int_t           nTightLep;
   Int_t           nLooseLep;
   Float_t         lepton_pt;
   Float_t         lepton_eta;
   Float_t         lepton_phi;
   Float_t         lepton_m;
   Int_t           lepton_isMuon;
   Int_t           lepton_charge;
   Float_t         lepton2_pt;
   Float_t         lepton2_eta;
   Float_t         lepton2_phi;
   Float_t         lepton2_m;
   Int_t           lepton2_isMuon;
   Int_t           lepton2_charge;
   Int_t           nJet;
   Int_t           nBJet;
   Float_t         met_abs;
   Float_t         met_phi;
   Float_t         met_type1xy_abs;
   Float_t         met_type1xy_phi;
   std::vector<float>   *jet_pt;
   std::vector<float>   *jet_eta;
   std::vector<float>   *jet_phi;
   std::vector<float>   *jet_m;
   std::vector<float>   *jet_btag;
   std::vector<int>     *jet_flav;
   Float_t         Mll;
   Float_t         lepjet_minDr;
   std::vector<float>   *fatjet_p;
   std::vector<float>   *fatjet_pt;
   std::vector<float>   *fatjet_eta;
   std::vector<float>   *fatjet_phi;
   std::vector<float>   *fatjet_sdmass;
   std::vector<float>   *fatjet_tau21;
   std::vector<float>   *fatjet_tau32;
   std::vector<float>   *fatjet_doubleb;
   std::vector<int>     *fatjet_nsubjets;
   std::vector<int>     *fatjet_subjet_nloosebtag;
   std::vector<int>     *fatjet_subjet_nmediumbtag;
   std::vector<float>   *fatjet_subjet_meanbtagger;
   std::vector<float>   *fatjet_subjet_minbtagger;
   Float_t         truth_higgs_pt;
   Float_t         truth_higgs_eta;
   Float_t         truth_higgs_phi;

   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_EventCategory;   //!
   TBranch        *b_TotalEventWeight;   //!
   TBranch        *b_id_additionalJetEventId;   //!
   TBranch        *b_TruePV;   //!
   TBranch        *b_GenPV;   //!
   TBranch        *b_N_Vtx;   //!
   TBranch        *b_wgt_TOTAL;   //!
   TBranch        *b_wgt_PU;   //!
   TBranch        *b_wgt_btag;   //!
   TBranch        *b_wgt_trigMu;   //!
   TBranch        *b_wgt_trigEl;   //!
   TBranch        *b_wgt_EGZtvx;   //!
   TBranch        *b_wgt_muID;   //!
   TBranch        *b_wgt_elID;   //!
   TBranch        *b_wgt_elReco;   //!
   TBranch        *b_wgt_muIso;   //!
   TBranch        *b_periodwgt_PU;   //!
   TBranch        *b_periodwgt_btag;   //!
   TBranch        *b_periodwgt_trigMu;   //!
   TBranch        *b_periodwgt_trigEl;   //!
   TBranch        *b_periodwgt_EGZtvx;   //!
   TBranch        *b_periodwgt_muID;   //!
   TBranch        *b_periodwgt_elID;   //!
   TBranch        *b_periodwgt_elReco;   //!
   TBranch        *b_periodwgt_muIso;   //!
   TBranch        *b_wgt_trigMuUP;   //!
   TBranch        *b_wgt_trigElUP;   //!
   TBranch        *b_wgt_muIDUP;   //!
   TBranch        *b_wgt_elIDUP;   //!
   TBranch        *b_wgt_elRecoUP;   //!
   TBranch        *b_wgt_muIsoUP;   //!
   TBranch        *b_wgt_trigMuDOWN;   //!
   TBranch        *b_wgt_trigElDOWN;   //!
   TBranch        *b_wgt_muIDDOWN;   //!
   TBranch        *b_wgt_elIDDOWN;   //!
   TBranch        *b_wgt_elRecoDOWN;   //!
   TBranch        *b_wgt_muIsoDOWN;   //!
   TBranch        *b_wgt_MCEventWeight;   //!
   TBranch        *b_wgt_topPT;   //!
   TBranch        *b_mcWeight_key;   //!
   TBranch        *b_mcWeight_value;   //!
   TBranch        *b_pass_goodVtx;   //!
   TBranch        *b_pass_El;   //!
   TBranch        *b_pass_Mu;   //!
   TBranch        *b_pass_ElEl;   //!
   TBranch        *b_pass_ElMu;   //!
   TBranch        *b_pass_MuMu;   //!
   TBranch        *b_pass_TrigEl;   //!
   TBranch        *b_pass_TrigMu;   //!
   TBranch        *b_pass_TrigElEl;   //!
   TBranch        *b_pass_TrigElMu;   //!
   TBranch        *b_pass_TrigMuMu;   //!
   TBranch        *b_nNonIsoEl;   //!
   TBranch        *b_nNonIsoMu;   //!
   TBranch        *b_nTightLep;   //!
   TBranch        *b_nLooseLep;   //!
   TBranch        *b_lepton_pt;   //!
   TBranch        *b_lepton_eta;   //!
   TBranch        *b_lepton_phi;   //!
   TBranch        *b_lepton_m;   //!
   TBranch        *b_lepton_isMuon;   //!
   TBranch        *b_lepton_charge;   //!
   TBranch        *b_lepton2_pt;   //!
   TBranch        *b_lepton2_eta;   //!
   TBranch        *b_lepton2_phi;   //!
   TBranch        *b_lepton2_m;   //!
   TBranch        *b_lepton2_isMuon;   //!
   TBranch        *b_lepton2_charge;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_nBJet;   //!
   TBranch        *b_met_abs;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_type1xy_abs;   //!
   TBranch        *b_met_type1xy_phi;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_m;   //!
   TBranch        *b_jet_btag;   //!
   TBranch        *b_jet_flav;   //!
   TBranch        *b_Mll;   //!
   TBranch        *b_lepjet_minDr;   //!
   TBranch        *b_fatjet_p;   //!
   TBranch        *b_fatjet_pt;   //!
   TBranch        *b_fatjet_eta;   //!
   TBranch        *b_fatjet_phi;   //!
   TBranch        *b_fatjet_sdmass;   //!
   TBranch        *b_fatjet_tau21;   //!
   TBranch        *b_fatjet_tau32;   //!
   TBranch        *b_fatjet_doubleb;   //!
   TBranch        *b_fatjet_nsubjets;   //!
   TBranch        *b_fatjet_subjet_nloosebtag;   //!
   TBranch        *b_fatjet_subjet_nmediumbtag;   //!
   TBranch        *b_fatjet_subjet_meanbtagger;   //!
   TBranch        *b_fatjet_subjet_minbtagger;   //!
   TBranch        *b_truth_higgs_pt;   //!
   TBranch        *b_truth_higgs_eta;   //!
   TBranch        *b_truth_higgs_phi;   //!

   base(TTree *tree=0);
   virtual ~base();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef base_cxx
base::base(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("out.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("out.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

base::~base()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t base::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t base::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void base::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   periodwgt_PU = 0;
   periodwgt_btag = 0;
   periodwgt_trigMu = 0;
   periodwgt_trigEl = 0;
   periodwgt_EGZtvx = 0;
   periodwgt_muID = 0;
   periodwgt_elID = 0;
   periodwgt_elReco = 0;
   periodwgt_muIso = 0;
   mcWeight_value = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_m = 0;
   jet_btag = 0;
   jet_flav = 0;
   fatjet_p = 0;
   fatjet_pt = 0;
   fatjet_eta = 0;
   fatjet_phi = 0;
   fatjet_sdmass = 0;
   fatjet_tau21 = 0;
   fatjet_tau32 = 0;
   fatjet_doubleb = 0;
   fatjet_nsubjets = 0;
   fatjet_subjet_nloosebtag = 0;
   fatjet_subjet_nmediumbtag = 0;
   fatjet_subjet_meanbtagger = 0;
   fatjet_subjet_minbtagger = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("EventCategory", &EventCategory, &b_EventCategory);
   fChain->SetBranchAddress("TotalEventWeight", &TotalEventWeight, &b_TotalEventWeight);
   fChain->SetBranchAddress("id_additionalJetEventId", &id_additionalJetEventId, &b_id_additionalJetEventId);
   fChain->SetBranchAddress("TruePV", &TruePV, &b_TruePV);
   fChain->SetBranchAddress("GenPV", &GenPV, &b_GenPV);
   fChain->SetBranchAddress("N_Vtx", &N_Vtx, &b_N_Vtx);
   fChain->SetBranchAddress("wgt_TOTAL", &wgt_TOTAL, &b_wgt_TOTAL);
   fChain->SetBranchAddress("wgt_PU", &wgt_PU, &b_wgt_PU);
   fChain->SetBranchAddress("wgt_btag", &wgt_btag, &b_wgt_btag);
   fChain->SetBranchAddress("wgt_trigMu", &wgt_trigMu, &b_wgt_trigMu);
   fChain->SetBranchAddress("wgt_trigEl", &wgt_trigEl, &b_wgt_trigEl);
   fChain->SetBranchAddress("wgt_EGZtvx", &wgt_EGZtvx, &b_wgt_EGZtvx);
   fChain->SetBranchAddress("wgt_muID", &wgt_muID, &b_wgt_muID);
   fChain->SetBranchAddress("wgt_elID", &wgt_elID, &b_wgt_elID);
   fChain->SetBranchAddress("wgt_elReco", &wgt_elReco, &b_wgt_elReco);
   fChain->SetBranchAddress("wgt_muIso", &wgt_muIso, &b_wgt_muIso);
   fChain->SetBranchAddress("periodwgt_PU", &periodwgt_PU, &b_periodwgt_PU);
   fChain->SetBranchAddress("periodwgt_btag", &periodwgt_btag, &b_periodwgt_btag);
   fChain->SetBranchAddress("periodwgt_trigMu", &periodwgt_trigMu, &b_periodwgt_trigMu);
   fChain->SetBranchAddress("periodwgt_trigEl", &periodwgt_trigEl, &b_periodwgt_trigEl);
   fChain->SetBranchAddress("periodwgt_EGZtvx", &periodwgt_EGZtvx, &b_periodwgt_EGZtvx);
   fChain->SetBranchAddress("periodwgt_muID", &periodwgt_muID, &b_periodwgt_muID);
   fChain->SetBranchAddress("periodwgt_elID", &periodwgt_elID, &b_periodwgt_elID);
   fChain->SetBranchAddress("periodwgt_elReco", &periodwgt_elReco, &b_periodwgt_elReco);
   fChain->SetBranchAddress("periodwgt_muIso", &periodwgt_muIso, &b_periodwgt_muIso);
   fChain->SetBranchAddress("wgt_trigMuUP", &wgt_trigMuUP, &b_wgt_trigMuUP);
   fChain->SetBranchAddress("wgt_trigElUP", &wgt_trigElUP, &b_wgt_trigElUP);
   fChain->SetBranchAddress("wgt_muIDUP", &wgt_muIDUP, &b_wgt_muIDUP);
   fChain->SetBranchAddress("wgt_elIDUP", &wgt_elIDUP, &b_wgt_elIDUP);
   fChain->SetBranchAddress("wgt_elRecoUP", &wgt_elRecoUP, &b_wgt_elRecoUP);
   fChain->SetBranchAddress("wgt_muIsoUP", &wgt_muIsoUP, &b_wgt_muIsoUP);
   fChain->SetBranchAddress("wgt_trigMuDOWN", &wgt_trigMuDOWN, &b_wgt_trigMuDOWN);
   fChain->SetBranchAddress("wgt_trigElDOWN", &wgt_trigElDOWN, &b_wgt_trigElDOWN);
   fChain->SetBranchAddress("wgt_muIDDOWN", &wgt_muIDDOWN, &b_wgt_muIDDOWN);
   fChain->SetBranchAddress("wgt_elIDDOWN", &wgt_elIDDOWN, &b_wgt_elIDDOWN);
   fChain->SetBranchAddress("wgt_elRecoDOWN", &wgt_elRecoDOWN, &b_wgt_elRecoDOWN);
   fChain->SetBranchAddress("wgt_muIsoDOWN", &wgt_muIsoDOWN, &b_wgt_muIsoDOWN);
   fChain->SetBranchAddress("wgt_MCEventWeight", &wgt_MCEventWeight, &b_wgt_MCEventWeight);
   fChain->SetBranchAddress("wgt_topPT", &wgt_topPT, &b_wgt_topPT);
   fChain->SetBranchAddress("mcWeight_value", &mcWeight_value, &b_mcWeight_value);
   fChain->SetBranchAddress("pass_goodVtx", &pass_goodVtx, &b_pass_goodVtx);
   fChain->SetBranchAddress("pass_El", &pass_El, &b_pass_El);
   fChain->SetBranchAddress("pass_Mu", &pass_Mu, &b_pass_Mu);
   fChain->SetBranchAddress("pass_ElEl", &pass_ElEl, &b_pass_ElEl);
   fChain->SetBranchAddress("pass_ElMu", &pass_ElMu, &b_pass_ElMu);
   fChain->SetBranchAddress("pass_MuMu", &pass_MuMu, &b_pass_MuMu);
   fChain->SetBranchAddress("pass_TrigEl", &pass_TrigEl, &b_pass_TrigEl);
   fChain->SetBranchAddress("pass_TrigMu", &pass_TrigMu, &b_pass_TrigMu);
   fChain->SetBranchAddress("pass_TrigElEl", &pass_TrigElEl, &b_pass_TrigElEl);
   fChain->SetBranchAddress("pass_TrigElMu", &pass_TrigElMu, &b_pass_TrigElMu);
   fChain->SetBranchAddress("pass_TrigMuMu", &pass_TrigMuMu, &b_pass_TrigMuMu);
   fChain->SetBranchAddress("nNonIsoEl", &nNonIsoEl, &b_nNonIsoEl);
   fChain->SetBranchAddress("nNonIsoMu", &nNonIsoMu, &b_nNonIsoMu);
   fChain->SetBranchAddress("nTightLep", &nTightLep, &b_nTightLep);
   fChain->SetBranchAddress("nLooseLep", &nLooseLep, &b_nLooseLep);
   fChain->SetBranchAddress("lepton_pt", &lepton_pt, &b_lepton_pt);
   fChain->SetBranchAddress("lepton_eta", &lepton_eta, &b_lepton_eta);
   fChain->SetBranchAddress("lepton_phi", &lepton_phi, &b_lepton_phi);
   fChain->SetBranchAddress("lepton_m", &lepton_m, &b_lepton_m);
   fChain->SetBranchAddress("lepton_isMuon", &lepton_isMuon, &b_lepton_isMuon);
   fChain->SetBranchAddress("lepton_charge", &lepton_charge, &b_lepton_charge);
   fChain->SetBranchAddress("lepton2_pt", &lepton2_pt, &b_lepton2_pt);
   fChain->SetBranchAddress("lepton2_eta", &lepton2_eta, &b_lepton2_eta);
   fChain->SetBranchAddress("lepton2_phi", &lepton2_phi, &b_lepton2_phi);
   fChain->SetBranchAddress("lepton2_m", &lepton2_m, &b_lepton2_m);
   fChain->SetBranchAddress("lepton2_isMuon", &lepton2_isMuon, &b_lepton2_isMuon);
   fChain->SetBranchAddress("lepton2_charge", &lepton2_charge, &b_lepton2_charge);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("nBJet", &nBJet, &b_nBJet);
   fChain->SetBranchAddress("met_abs", &met_abs, &b_met_abs);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_type1xy_abs", &met_type1xy_abs, &b_met_type1xy_abs);
   fChain->SetBranchAddress("met_type1xy_phi", &met_type1xy_phi, &b_met_type1xy_phi);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_m", &jet_m, &b_jet_m);
   fChain->SetBranchAddress("jet_btag", &jet_btag, &b_jet_btag);
   fChain->SetBranchAddress("jet_flav", &jet_flav, &b_jet_flav);
   fChain->SetBranchAddress("Mll", &Mll, &b_Mll);
   fChain->SetBranchAddress("lepjet_minDr", &lepjet_minDr, &b_lepjet_minDr);
   fChain->SetBranchAddress("fatjet_p", &fatjet_p, &b_fatjet_p);
   fChain->SetBranchAddress("fatjet_pt", &fatjet_pt, &b_fatjet_pt);
   fChain->SetBranchAddress("fatjet_eta", &fatjet_eta, &b_fatjet_eta);
   fChain->SetBranchAddress("fatjet_phi", &fatjet_phi, &b_fatjet_phi);
   fChain->SetBranchAddress("fatjet_sdmass", &fatjet_sdmass, &b_fatjet_sdmass);
   fChain->SetBranchAddress("fatjet_tau21", &fatjet_tau21, &b_fatjet_tau21);
   fChain->SetBranchAddress("fatjet_tau32", &fatjet_tau32, &b_fatjet_tau32);
   fChain->SetBranchAddress("fatjet_doubleb", &fatjet_doubleb, &b_fatjet_doubleb);
   fChain->SetBranchAddress("fatjet_nsubjets", &fatjet_nsubjets, &b_fatjet_nsubjets);
   fChain->SetBranchAddress("fatjet_subjet_nloosebtag", &fatjet_subjet_nloosebtag, &b_fatjet_subjet_nloosebtag);
   fChain->SetBranchAddress("fatjet_subjet_nmediumbtag", &fatjet_subjet_nmediumbtag, &b_fatjet_subjet_nmediumbtag);
   fChain->SetBranchAddress("fatjet_subjet_meanbtagger", &fatjet_subjet_meanbtagger, &b_fatjet_subjet_meanbtagger);
   fChain->SetBranchAddress("fatjet_subjet_minbtagger", &fatjet_subjet_minbtagger, &b_fatjet_subjet_minbtagger);
   fChain->SetBranchAddress("truth_higgs_pt", &truth_higgs_pt, &b_truth_higgs_pt);
   fChain->SetBranchAddress("truth_higgs_eta", &truth_higgs_eta, &b_truth_higgs_eta);
   fChain->SetBranchAddress("truth_higgs_phi", &truth_higgs_phi, &b_truth_higgs_phi);
   Notify();
}

Bool_t base::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void base::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t base::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef base_cxx
