#ifndef MY_ANALYZR
#define MY_ANALYZR


#include <vector>
#include <string>
#include <map>

#include <TH1D.h>
#include <TH2D.h>

#include "reader.h"

#include "ttHYggdrasilScaleFactors.h"


#define ENABLE_COMMON_CLASSFIER

#ifdef ENABLE_COMMON_CLASSFIER
#include "TTH/CommonClassifier/interface/BDTClassifier.h"
#endif



//#define  TRUTH_INFO_STUDY

#include "TMVA/Reader.h"


class analyzer {

 private : 

  std::vector<std::string> filelist  ;
  std::string output ; 
  reader * r ; 

  double totalEvent ;

  bool isMC ;
  bool isMuonStream ;
  
  long n_th_event;

#ifdef ENABLE_COMMON_CLASSFIER
  BDTClassifier bdt;
#endif

  TMVA::Reader * tmva_reader ;
  float tmva_ak8jet_pt ;	 
  float tmva_ak8jet_eta ; 	 
  float tmva_ak8jet_sdmass ;	 
  float tmva_ak8jet_deepcsv_min ; 
  float tmva_ak8jet_doubleb ;    



 public : 

  analyzer ();
  ~analyzer ();

  void AnalyzeEvent();
  void SetPUReweightingFile( std::string data_path , std::string mc_path );

  void SetReader( reader * _r );

  void init();
  void postProcess();

  inline void SetOutputfileName( std::string outputpath ){ output.assign(outputpath); }
  inline void AddNumberOfProcessedEvent( double N ){ totalEvent += N; }
  void SetDataModeOn();

  void SetIsMuonStream( bool _isMuonStream  );

  void SetFakeEstimationModeOn();


  enum TtbarAdditionalJetID{
    noRequirement , 
    ttbarPlusB,
    ttbarPlus2B,
    ttbarPlusBBbar,
    ttbarPlusCCbar,
    ttbarOther
  };

  void setTtbarAdditionalJetIDCut( TtbarAdditionalJetID id );

  void SetPUReweightSyst( int PU_syst );

  void EnableSkippingOddEventNumber();

  void setSFSystematic( int syst );

  void SetDileptonFakeLeptonAnalysisMode();

  void SetMCPileupChannel( std::string name );

  void SetPeriod( int p ); // starting from 0 = Period-B, negative value = inclysive analysis.

 private :

  void FillHistogram( TH1D * h , double val , double weight );
  void FillHistogram2D( TH2D * h , double val_x , double val_y,  double weight );

  bool _passTtbarAdditionalJetID();

  std::string _NameofCategory( int cate );

  double _getFakeSF( bool mu_channel = true ); // In case of data, this flag is not used.

  float _calcDR2( float eta1, float eta2, float phi1, float phi2 );
  float miniDR( float basis_eta ,  std::vector<float> * etas, float basis_phi , std::vector<float> * phis );

  bool b_FakeEstimation ;
  TH2D * h_FakeSF_Mu ;
  TH2D * h_FakeSF_El ;

  TFile * f_out ;
  TDirectory * main_directory ;

  void work_for_truthInfo();

  void work_for_singleLepton();
  void work_for_fakeevaluation();
  void work_for_fakeeval_Dilepton();

  int npT ;

  TTree * outtree_odd ; 
  TTree * outtree_even ; 

  TH1D * h_totalevent; 
  TH1D * h_totalevent_noWeight; 

  TH1D * h_pileup_noWeight;
  TH1D * h_pileup_PUWeight;


  TH1D * h_SL_cutflow[2]; // el, mu channels.

  TH1D * h_SL_EventCategorizationBasedonNjetNBtagJet ; 
  TH1D * h_SL_EventCategorizationBasedonNjetNBtagJet_El ; 
  TH1D * h_SL_EventCategorizationBasedonNjetNBtagJet_Mu ; 
  TH1D * h_SL_EventCategorizationBasedonNjetNBtagJet_FakeCR_Mu ; 
  TH1D * h_SL_EventCategorizationBasedonNjetNBtagJet_FakeCR_El ; 

  TH1D * h_SL_EventCategorizationBasedonNjetNBtagJet_withMetCut ; 
  TH1D * h_SL_EventCategorizationBasedonNjetNBtagJet_withMetCut_el ; 
  TH1D * h_SL_EventCategorizationBasedonNjetNBtagJet_withMetCut_mu ; 


  TH1D * h_SL_Much_JetPtCategorization4to6 ;
  TH1D * h_SL_Much_JetPtCategorization3to6 ; 
  TH1D * h_SL_Much_JetPtCategorization2to4and3to2;
  TH1D * h_SL_Elch_JetPtCategorization4to6 ;
  TH1D * h_SL_Elch_JetPtCategorization3to6 ; 
  TH1D * h_SL_Elch_JetPtCategorization2to4and3to2;

  TH1D * h_SL_MET_AllButJetRequirement_jet0_mu ;
  TH1D * h_SL_MET_AllButJetRequirement_jet1_mu ;
  TH1D * h_SL_MET_AllButJetRequirement_jet2_mu ;
  TH1D * h_SL_MET_AllButJetRequirement_jet3_mu ;
  TH1D * h_SL_MET_AllButJetRequirement_jet4_mu ;
  TH1D * h_SL_MET_AllButJetRequirement_jet5_mu ;

  TH1D * h_SL_N_Vtx_AfterSingleLepTrigAndTightLep ;
  TH1D * h_SL_N_VtxWITHOUTPU_AfterSingleLepTrigAndTightLep ;

  TH1D * h_SL_N_Vtx_AfterAllSelection ;
  TH1D * h_SL_N_VtxWITHOUTPU_AfterAllSelection ;

  TH1D * h_SL_N_ManuallyReducedVtx_AfterAllSelection;

  TH1D * h_SL_isMuonAfterAllSelection ;
  TH1D * h_SL_MuonPtAfterAllSelection ;
  TH1D * h_SL_ElPtAfterAllSelection   ;
  TH1D * h_SL_MuonEtaAfterAllSelection;
  TH1D * h_SL_MuonPhiAfterAllSelection;
  TH1D * h_SL_ElEtaAfterAllSelection  ;
  TH1D * h_SL_ElPhiAfterAllSelection  ;
  
  TH1D * h_SL_MuonPtAfterAllSelection_ge3btag ;
  TH1D * h_SL_ElPtAfterAllSelection_ge3btag  ;

  TH1D * h_SL_ElPtAfterAllSelection_noElSF ; 
  TH1D * h_SL_ElEtaAfterAllSelection_noElSF  ;

  TH1D * h_SL_MuonPt_NoTopPTWgt_AfterAllSelection ;
  TH1D * h_SL_ElPt_NoTopPTWgt_AfterAllSelection   ;
  TH1D * h_SL_MuonEta_NoTopPTWgt_AfterAllSelection;
  TH1D * h_SL_ElEta_NoTopPTWgt_AfterAllSelection  ;

  TH1D * h_SL_MuonPt_NoMuTrigWgt_AfterAllSelection ;
  TH1D * h_SL_MuonEta_NoMuTrigWgt_AfterAllSelection;
  TH1D * h_SL_ElEta_NoElTrigWgt_AfterAllSelection  ;
  TH1D * h_SL_ElPt_NoElTrigWgt_AfterAllSelection   ;


  TH1D * h_SL_nJet_AfterAllSelection  ;
  TH1D * h_SL_nJet_AfterAllSelection_el  ;
  TH1D * h_SL_nJet_AfterAllSelection_mu  ;
  TH1D * h_SL_nJet_NoTopPTWgt_AfterAllSelection  ;
  TH1D * h_SL_nJetWITHOUTPUWGT_AfterAllSelection  ;

  TH1D *  h_SL_nJet_AllButJetRequirement  ;
  TH1D *  h_SL_nJet_mu_AllButJetRequirement;
  TH1D *  h_SL_nJet_el_AllButJetRequirement;


  TH1D * h_SL_Entries_nJnBCategoryAfterLeptonCut_MetCut_Mu ; 
  TH1D * h_SL_MuonPT_nJnBCategoryAfterLeptonCut_MetCut [7][5];// within 6j4b category
  TH1D * h_SL_Met_nJnBCategoryAfterLeptonCut_Mu[7][5];// within 6j4b category

  TH1D *  h_SL_MuonPT_AllButJetRequirement;
  TH1D *  h_SL_ElPT_AllButJetRequirement  ;
  TH1D *  h_SL_Met_AllButJetRequirement_el;
  TH1D *  h_SL_Met_AllButJetRequirement_mu;
  
  TH1D * h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet0 ; 
  TH1D * h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet1 ; 
  TH1D * h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet2 ; 
  TH1D * h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet3 ; 
  TH1D * h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet4 ; 
  TH1D * h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet5more ; 

  TH1D *  h_SL_MuonPT_AllButNBJetRequirement;
  TH1D *  h_SL_ElPT_AllButNBJetRequirement  ;
  TH1D *  h_SL_Met_AllButNBJetRequirement_el;
  TH1D *  h_SL_Met_AllButNBJetRequirement_mu;

  TH1D * h_SL_MetType1xy_AllButJetRequirement_el;
  TH1D * h_SL_MetType1xy_AllButJetRequirement_mu;

  TH1D * h_SL_nBtagJet_AfterAllSelection ;
  TH1D * h_SL_nBtagJet_AfterAllSelection_el ;
  TH1D * h_SL_nBtagJet_AfterAllSelection_mu ;
  TH1D * h_SL_nBtagJetWITHOUTBTAGSF_AfterAllSelection ;

  TH1D * h_SL_MetAbs_AfterAllSelection   ;
  TH1D * h_SL_MetAbs_AfterAllSelection_el;
  TH1D * h_SL_MetAbs_AfterAllSelection_mu;
  TH1D * h_SL_MetPhi_AfterAllSelection   ;
  TH1D * h_SL_MetPhi_AfterAllSelection_el;
  TH1D * h_SL_MetPhi_AfterAllSelection_mu;

  TH1D * h_SL_MetType1xyAbs_AfterAllSelection   ;
  TH1D * h_SL_MetType1xyAbs_AfterAllSelection_el;
  TH1D * h_SL_MetType1xyAbs_AfterAllSelection_mu;
  TH1D * h_SL_MetType1xyPhi_AfterAllSelection   ;
  TH1D * h_SL_MetType1xyPhi_AfterAllSelection_el;
  TH1D * h_SL_MetType1xyPhi_AfterAllSelection_mu;

  TH1D * h_SL_MetPhi_abs0to30GeV_AfterAllSelection;
  TH1D * h_SL_MetPhi_abs30to100GeV_AfterAllSelection;
  TH1D * h_SL_MetPhi_abs100GeVto200_AfterAllSelection;
  TH1D * h_SL_MetPhi_abs200GeVtoinf_AfterAllSelection;

  TH1D * h_SL_MetEx_AfterAllSelection   ;
  TH1D * h_SL_MetEx_AfterAllSelection_el;
  TH1D * h_SL_MetEx_AfterAllSelection_mu;

  TH1D * h_SL_MetEy_AfterAllSelection   ;
  TH1D * h_SL_MetEy_AfterAllSelection_el;
  TH1D * h_SL_MetEy_AfterAllSelection_mu;


  TH1D * h_SL_jet_pt_AfterAllSelection [8] ;
  TH1D * h_SL_jet_eta_AfterAllSelection [8];
  TH1D * h_SL_jetBtag_AfterAllSelection[8];

  TH1D * h_SL_jet_pt_AfterAllSelection_ge3btag [8] ;

  TH1D * h_SL_jet_pt_AfterAllSelection_mu [8] ;
  TH1D * h_SL_jet_eta_AfterAllSelection_mu [8];
  TH1D * h_SL_jetBtag_AfterAllSelection_mu [8];

  TH1D * h_SL_jet_pt_AfterAllSelection_el [8] ;
  TH1D * h_SL_jet_eta_AfterAllSelection_el [8];
  TH1D * h_SL_jetBtag_AfterAllSelection_el[8];


  TH1D * h_SL_jet_pt_WOtopPTwgt_AfterAllSelection[8];
  TH1D * h_SL_jet_eta_WOtopPTwgt_AfterAllSelection[8];
  TH1D * h_SL_jetBtagWOBtagSF_AfterAllSelection[8];
  TH1D * h_SL_jetBtagWOBtagSF_AfterAllSelection_mu[8];
  TH1D * h_SL_jetBtagWOBtagSF_AfterAllSelection_el[8];

  TH1D * h_SL_jet_pt_AfterAllButBtagJetRequirement[8];
  TH1D * h_SL_jetBtag_AfterAllButBtagJetRequirement[8];


  std::vector< TH1D * > histogram_list ;

  long _EventCateBasedOnNjetNBtagJet( long nJ , long nB );

  TH1D * h_wgt_TOTAL ;
  TH1D * h_wgt_PU    ;
  TH1D * h_wgt_btag  ;
  TH1D * h_wgt_trigMu;
  TH1D * h_wgt_trigEl;
  TH1D * h_wgt_muID  ;
  TH1D * h_wgt_elID  ;
  TH1D * h_wgt_elReco;
  TH1D * h_wgt_muIso ;
  TH1D * h_wgt_MCEven;
  TH1D * h_wgt_topPT ;
  TH1D * h_wgt_EGZvtx;

  TH1D *  h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons ;
  TH1D *  h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_noElSF ;
  TH1D *  h_SL_Mmm_AfterSingleTriggerAndExactTwoTightleptons ;
  TH1D *  h_SL_Mmm_AfterSingleTriggerAndExactTwoTightleptons_noMuSF ;

  TH1D *  h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElRecoSF ;
  TH1D *  h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElTrigSF ;
  TH1D *  h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElIDSF   ;

  TH1D * h_SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut ;
  TH1D * h_SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut ;

  TH1D * h_SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut_noElSF ;
  TH1D * h_SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut_noElSF ;

  TH1D * h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut ;
  TH1D * h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut ;

  TH1D * h_SL_nJet_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut;
  TH1D * h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut ;
  TH1D * h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut ;
  TH1D * h_SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut ;
  TH1D * h_SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut ;
  TH1D * h_SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut ;

  TH2D * h2_SL_NPV_NJet_AfterLeptonEventSelection_mu;

  TH1D * h_SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;
  TH1D * h_SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;
  TH1D * h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;
  TH1D * h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;

  TH1D * h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuSF ;
  TH1D * h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuSF ;

  TH1D * h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuTrigSF ;
  TH1D * h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuTrigSF ;

  TH1D *  h_SL_InvMass_AfterMuMuSameSignSelection             ;
  TH1D *  h_SL_MuonPT__AfterMuMuSameSignSelection_ZmassCut    ;
  TH1D *  h_SL_MuonEta_AfterMuMuSameSignSelection_ZmassCut    ;

  TH1D *  h_SL_InvMass_AfterElElSameSignSelection 	    ;	
  TH1D *  h_SL_ElectronPT__AfterElElSameSignSelection_ZmassCut;
  TH1D *  h_SL_ElectronEta_AfterElElSameSignSelection_ZmassCut;
  
  TH1D * h_SL_El1_PT_WWenrichRegion;
  TH1D * h_SL_El2_PT_WWenrichRegion;
  TH1D * h_SL_Mu1_PT_WWenrichRegion;
  TH1D * h_SL_Mu2_PT_WWenrichRegion;

  TH1D * h_SL_LeptopnPT_WWemuEnrichRegionLoose ;
  TH1D * h_MET_WWemuEnrichRegionLoose ;
  TH1D * h_SL_LeptopnPT_WWemuEnrichRegionTight ;
  TH1D * h_MET_WWemuEnrichRegionTight ;

  TH1D *  h_SL_MuonPT_WWemuSameSignEnrichRegionLoose  ;
  TH1D *  h_SL_ElectronPT_WWemuSameSignEnrichRegionLoose  ;
  TH1D *  h_MET_WWemuSameSignEnrichRegionLoose ;
  TH1D *  h_nJet_WWemuSameSignEnrichRegionLoose;
  TH1D *  h_nBJet_WWemuSameSignEnrichRegionLoose;

  TH1D * h_SL_nJet_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut;
  TH1D * h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;
  TH1D * h_SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;
  TH1D * h_SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;
  TH1D * h_SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;
  TH1D * h_SL_JetPt_5_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;
  
  TH1D * h_SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;
  TH1D * h_SL_MetPhi_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut ;
  
  TH1D * h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet0     ;
  TH1D * h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet1     ;
  TH1D * h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet2     ;
  TH1D * h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet3     ;
  TH1D * h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet4     ;
  TH1D * h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet5more ;

  TH1D * h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut ;
  TH1D * h_SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut ;
  TH1D * h_SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut ;
  TH1D * h_SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut ;

  TH1D * h_SL_jetPt_upto8th_BtagLFControlRegion ;
  TH1D * h_SL_jetPt_upto8th_BtagHFControlRegion ;
  TH1D * h_SL_jetBtag_upto8th_BtagLFControlRegion ;
  TH1D * h_SL_jetBtag_upto8th_BtagHFControlRegion ;
  TH1D * h_SL_jetBtag_upto8th_BtagLFControlRegionWOBtagSF ;
  TH1D * h_SL_jetBtag_upto8th_BtagHFControlRegionWOBtagSF ;

  #ifdef ENABLE_COMMON_CLASSFIER
  TH1D * h_SL_BDT[10];
  TH1D * h_SL_BDT_metcut[10];
  TH1D * h_SL_BDT_metcut_el[10];
  TH1D * h_SL_BDT_metcut_mu[10];
  TH1D * h_SL_BDT_nometcut_el[10];
  TH1D * h_SL_BDT_nometcut_mu[10];
  TH1D * h_SL_ElPtAfterCate_MetCut[10];
  TH1D * h_SL_MuonPtAfterCate_MetCut[10];
  #endif 

  TH1D * h_SL_fatjet_BDT_metcut_withInfoSDMassinHiggsWindow[10];
  TH1D * h_SL_fatjet_BDT_metcut_SDMass_in_HiggsWindow[10];
  TH1D * h_SL_fatjet_BDT_metcut_SDMass_out_HiggsWindow[10];

  TH1D * h_SL_fatjet_BDT_ExactOneBBFatJetPass_SDMassPass [ 10 ] ;
  TH1D * h_SL_fatjet_BDT_ExactOneBBFatJetPass_SDMassFail [ 10 ] ;
  TH1D * h_SL_fatjet_BDT_ExactOneBBFatJetPass [ 10 ] ;
  TH1D * h_SL_fatjet_BDT_ExactOneBBFatjetFail [ 10 ] ;


  TH1D* h_SL_MuPt_AfterNjNbCategory[10]     ;
  TH1D* h_SL_MuEta_AfterNjNbCategory[10]     ;
  TH1D* h_SL_MinMuJetDr_AfterNjNbCategory[10]     ;

  TH1D* h_SL_ElPt_AfterNjNbCategory[10]     ;
  TH1D* h_SL_ElEta_AfterNjNbCategory[10]     ;
  TH1D* h_SL_MinElJetDr_AfterNjNbCategory[10]     ;


  TH1D*  h_SL_El_LeadingJetPT_AfterNjNbCategory[10]    ;
  TH1D*  h_SL_El_LeastJetPT_AfterNjNbCategory[10]      ;
  TH1D*  h_SL_El_LeadingJetBtag_AfterNjNbCategory[10]  ;
  TH1D*  h_SL_El_LeastJetBtag_AfterNjNbCategory[10]    ;
  TH1D*  h_SL_El_LargestBtagValue_AfterNjNbCategory[10];
  TH1D*  h_SL_El_LeastBtagValue_AfterNjNbCategory[10]  ;

  TH1D*  h_SL_Mu_LeadingJetPT_AfterNjNbCategory[10]    ;
  TH1D*  h_SL_Mu_LeastJetPT_AfterNjNbCategory[10]      ;
  TH1D*  h_SL_Mu_LeadingJetBtag_AfterNjNbCategory[10]  ;
  TH1D*  h_SL_Mu_LeastJetBtag_AfterNjNbCategory[10]    ;
  TH1D*  h_SL_Mu_LargestBtagValue_AfterNjNbCategory[10];
  TH1D*  h_SL_Mu_LeastBtagValue_AfterNjNbCategory[10]  ;


  TH1D* h_SL_JetPt_AfterNjNbCategory[10][6]      ;
  TH1D*  h_SL_JetBtag_AfterNjNbCategory[10][6]   ;
  TH1D*  h_SL_Btagorder_AfterNjNbCategory[10][6] ; 


  TH1D* h_SL_JetPt_AfterNjNbCategory_el[10][6]      ;
  TH1D*  h_SL_JetBtag_AfterNjNbCategory_el[10][6]   ;
  TH1D*  h_SL_Btagorder_AfterNjNbCategory_el[10][6] ; 

  TH1D* h_SL_JetPt_AfterNjNbCategory_mu[10][6]      ;
  TH1D*  h_SL_JetBtag_AfterNjNbCategory_mu[10][6]   ;
  TH1D*  h_SL_Btagorder_AfterNjNbCategory_mu[10][6] ; 

  TH1D * h_FakeEval_Met_NonIsoEl_njnb[ 7 ][ 7 ] ;
  TH1D * h_FakeEval_Met_IsoEl_njnb   [ 7 ][ 7 ] ;
  TH1D * h_FakeEval_Met_NonIsoMu_njnb[ 7 ][ 7 ]	;
  TH1D * h_FakeEval_Met_IsoMu_njnb   [ 7 ][ 7 ] ;

  TH1D * h_FakeEval_dR_Jet_and_Mu_onejetcut_noFakeWgt_1j0b ;
  TH1D * h_FakeEval_dR_Jet_and_Mu_onejetcut_noFakeWgt_1j1b ; 
  TH1D * h_FakeEval_dR_Jet_and_El_onejetcut_noFakeWgt_1j0b ;
  TH1D * h_FakeEval_dR_Jet_and_El_onejetcut_noFakeWgt_1j1b ; 


  TH1D * h_FakeEval_dR_Jet_and_Mu_onejetcut_noFakeWgt [7][7];
  TH1D * h_FakeEval_dR_Btag_and_Mu_onejetcut_noFakeWgt[7][7];

  TH1D * h_FakeEval_dR_Jet_and_El_onejetcut_noFakeWgt [7][7];
  TH1D * h_FakeEval_dR_Btag_and_El_onejetcut_noFakeWgt[7][7];

  TH1D * h_FakeEval_nJet_El_onejetcut_noFakeWgt ;
  TH1D * h_FakeEval_nBJet_El_onejetcut_noFakeWgt;
  TH1D * h_FakeEval_Met_El_njnb_noFakeWgt    [7][7];
  TH1D * h_FakeEval_NEntries_El_LowMet_noFakeWgt[7][7];
  TH1D * h_FakeEval_nJet_El_onejetcut_withFakeWgt ;
  TH1D * h_FakeEval_nBJet_El_onejetcut_withFakeWgt;
  TH1D * h_FakeEval_Met_El_njnb_withFakeWgt    [ 7 ][ 7 ];
  TH1D * h_FakeEval_NEntries_El_LowMet_withFakeWgt[ 7 ][ 7 ] ; 

  TH1D * h_FakeEval_nJet_Mu_onejetcut_noFakeWgt ;
  TH1D * h_FakeEval_nBJet_Mu_onejetcut_noFakeWgt;
  TH1D * h_FakeEval_Met_Mu_njnb_noFakeWgt    [7][7];
  TH1D * h_FakeEval_NEntries_Mu_LowMet_noFakeWgt[7][7];
  TH1D * h_FakeEval_NEntriesVsNPV_Mu_LowMet_noFakeWgt[7][7];
  TH1D * h_FakeEval_nJet_Mu_onejetcut_withFakeWgt ;
  TH1D * h_FakeEval_nBJet_Mu_onejetcut_withFakeWgt;
  TH1D * h_FakeEval_Met_Mu_njnb_withFakeWgt    [ 7 ][ 7 ];
  TH1D * h_FakeEval_NEntries_Mu_LowMet_withFakeWgt[ 7 ][ 7 ] ; 

  TH2D * h2_FakeEval_NEntries_El_nJnB_withFakeWgt;
  TH2D * h2_FakeEval_NEntries_El_nJnB_noFakeWgt  ;
  TH2D * h2_FakeEval_NEntries_Mu_nJnB_withFakeWgt;
  TH2D * h2_FakeEval_NEntries_Mu_nJnB_noFakeWgt  ;

  TH2D * h2_FakeEval_NEntries_LowMetSide_El_nJnB_withFakeWgt ;
  TH2D * h2_FakeEval_NEntries_LowMetSide_El_nJnB_noFakeWgt   ;
  TH2D * h2_FakeEval_NEntries_LowMetSide_Mu_nJnB_withFakeWgt ;
  TH2D * h2_FakeEval_NEntries_LowMetSide_Mu_nJnB_noFakeWgt   ;

  TH2D * h2_FakeEval_NEntries_HighMetSide_El_nJnB_withFakeWgt;
  TH2D * h2_FakeEval_NEntries_HighMetSide_El_nJnB_noFakeWgt  ;
  TH2D * h2_FakeEval_NEntries_HighMetSide_Mu_nJnB_withFakeWgt;
  TH2D * h2_FakeEval_NEntries_HighMetSide_Mu_nJnB_noFakeWgt  ;

  TH2D * h2_FakeEval_NonisoEl_njetNbin ;
  TH2D * h2_FakeEval_NonisoEl_njetNbinCutOff ;
  TH2D * h2_FakeEval_NonisoEl_njetNbin_LowMet ;
  TH2D * h2_FakeEval_NonisoEl_njetNbinCutOff_LowMet ;
  TH2D * h2_FakeEval_NonisoEl_njetNbin_notLowMet ;
  TH2D * h2_FakeEval_NonisoEl_njetNbinCutOff_notLowMet ;

  TH2D * h2_FakeEval_IsoEl_njetNbin ;
  TH2D * h2_FakeEval_IsoEl_njetNbinCutOff ;
  TH2D * h2_FakeEval_IsoEl_njetNbin_LowMet ;
  TH2D * h2_FakeEval_IsoEl_njetNbinCutOff_LowMet ;
  TH2D * h2_FakeEval_IsoEl_njetNbin_notLowMet ;
  TH2D * h2_FakeEval_IsoEl_njetNbinCutOff_notLowMet ;

  TH2D * h2_FakeEval_NonisoMu_njetNbin ;
  TH2D * h2_FakeEval_NonisoMu_njetNbinCutOff ;
  TH2D * h2_FakeEval_NonisoMu_njetNbin_LowMet ;
  TH2D * h2_FakeEval_NonisoMu_njetNbinCutOff_LowMet ;
  TH2D * h2_FakeEval_NonisoMu_njetNbin_notLowMet ;
  TH2D * h2_FakeEval_NonisoMu_njetNbinCutOff_notLowMet ;

  TH2D * h2_FakeEval_IsoMu_njetNbin ;
  TH2D * h2_FakeEval_IsoMu_njetNbinCutOff ;
  TH2D * h2_FakeEval_IsoMu_njetNbin_LowMet ;
  TH2D * h2_FakeEval_IsoMu_njetNbinCutOff_LowMet ;
  TH2D * h2_FakeEval_IsoMu_njetNbin_notLowMet ;
  TH2D * h2_FakeEval_IsoMu_njetNbinCutOff_notLowMet ;

  TtbarAdditionalJetID ttbarAdditionalJetID; 
  

  ttHYggdrasilScaleFactors ttHSF ;

  bool flag_skip_oddEventNumber ; 

  int syst_muTrigSF; 
  int syst_elTrigSF; 
  int syst_muIDSF  ; 
  int syst_elIDSF  ;
  int syst_muIsoSF;
  int syst_elRecoSF;
  int syst_FakeMuStat ; 
  int syst_FakeElStat ; 


  TH1D * h_ee_SameSign_SanityCheck_NlooseLepMinusTight;
  TH1D * h_ee_SameSign_SanityCheck_NlooseLepMinusNonIsoTight;
  TH1D * h_ee_SameSign_BtagjetMultiplicity;
  TH1D * h_ee_SameSign_BtagjetMultiplicity_MllCut;
  TH1D * h_ee_SameSign_mll ; 
  TH1D * h_ee_SameSign_ZmassWindow  [2][5];
  TH1D * h_ee_SameSign_ElectronPT   [2][5];
  TH1D * h_ee_SameSign_LeadingJetPT [2][5];


  bool b_DileptonFakeLeptonAnalysisMode;

  TH1D *  h_SL_eventCate_Mu_NoWeight;
  TH1D *  h_SL_eventCate_El_NoWeight;
  TH1D *  h_SL_eventCate_Mu_NoWeight_noMetCut;
  TH1D *  h_SL_eventCate_El_NoWeight_noMetCut;

  int period ;

  ttHYggdrasilScaleFactors::MCWeightKey MCGenWeightKey ; 

  TH1D *  h_SL_fatjet_pt_AfterAllSelection_el                  ; 
  TH1D *  h_SL_fatjet_eta_AfterAllSelection_el		     ; 
  TH1D *  h_SL_fatjet_phi_AfterAllSelection_el		     ; 
  TH1D *  h_SL_fatjet_tau21_AfterAllSelection_el		     ; 
  TH1D *  h_SL_fatjet_tau32_AfterAllSelection_el		     ; 
  TH1D *  h_SL_fatjet_sdmass_AfterAllSelection_el		     ; 
  TH1D *  h_SL_fatjet_subjet_loosetagged_AfterAllSelection_el  ; 
  TH1D *  h_SL_fatjet_subjet_mediumtagged_AfterAllSelection_el ; 
  TH1D *  h_SL_fatjet_subjet_mean_btagger_AfterAllSelection_el ; 
  TH1D *  h_SL_N_fatjet_AfterAllSelection_el		     ; 

  TH1D *  h_SL_fatjet_pt_AfterAllSelection_mu		     ; 
  TH1D *  h_SL_fatjet_eta_AfterAllSelection_mu		     ; 
  TH1D *  h_SL_fatjet_phi_AfterAllSelection_mu		     ; 
  TH1D *  h_SL_fatjet_tau21_AfterAllSelection_mu		     ; 
  TH1D *  h_SL_fatjet_tau32_AfterAllSelection_mu		     ; 
  TH1D *  h_SL_fatjet_sdmass_AfterAllSelection_mu		     ; 
  TH1D *  h_SL_fatjet_subjet_loosetagged_AfterAllSelection_mu  ; 
  TH1D *  h_SL_fatjet_subjet_mediumtagged_AfterAllSelection_mu ; 
  TH1D *  h_SL_fatjet_subjet_mean_btagger_AfterAllSelection_mu ; 
  TH1D *  h_SL_N_fatjet_AfterAllSelection_mu                   ; 


  TH1D * h_SL_fatjet_AmbitiousSelection_Cate6j4b_TwoFatJets   ;

  TH1D * h_SL_fatjet_IsSDMassInHiggsWindow[10];

  TH1D * h_SL_fatjet_SDMass_ofLargeSDmassJet_Cate6j4b_TwoFatJets   ;
  TH1D * h_SL_fatjet_SDMass_ofSmallSDmassJet_Cate6j4b_TwoFatJets   ;
  TH1D * h_SL_fatjet_pt_ofLargeSDmassJet_Cate6j4b_TwoFatJets   ;
  TH1D * h_SL_fatjet_pt_ofSmallSDmassJet_Cate6j4b_TwoFatJets   ;
  TH1D * h_SL_fatjet_eta_ofLargeSDmassJet_Cate6j4b_TwoFatJets   ;
  TH1D * h_SL_fatjet_eta_ofSmallSDmassJet_Cate6j4b_TwoFatJets   ;
  TH1D * h_SL_fatjet_tau32_ofLargeSDmassJet_Cate6j4b_TwoFatJets    ;
  TH1D * h_SL_fatjet_tau32_ofSmallSDmassJet_Cate6j4b_TwoFatJets    ;
  TH1D * h_SL_fatjet_tau21_ofLargeSDmassJet_Cate6j4b_TwoFatJets    ;
  TH1D * h_SL_fatjet_tau21_ofSmallSDmassJet_Cate6j4b_TwoFatJets    ;
  TH1D * h_SL_fatjet_NlooseB_ofLargeSDmassJet_Cate6j4b_TwoFatJets  ;
  TH1D * h_SL_fatjet_NlooseB_ofSmallSDmassJet_Cate6j4b_TwoFatJets  ;
  TH1D * h_SL_fatjet_NMediumB_ofLargeSDmassJet_Cate6j4b_TwoFatJets ;
  TH1D * h_SL_fatjet_NMediumB_ofSmallSDmassJet_Cate6j4b_TwoFatJets ;
  TH1D * h_SL_fatjet_dRlep_ofLargeSDmassJet_Cate6j4b_TwoFatJets ;
  TH1D * h_SL_fatjet_dRlep_ofSmallSDmassJet_Cate6j4b_TwoFatJets ;




  TH1D *  h_SL_fatjet_pt_NbNjCate                 [10] ; 
  TH1D *  h_SL_fatjet_eta_NbNjCate		     [10]; 
  TH1D *  h_SL_fatjet_phi_NbNjCate		     [10]; 
  TH1D *  h_SL_fatjet_tau21_NbNjCate		     [10]; 
  TH1D *  h_SL_fatjet_tau32_NbNjCate		     [10]; 
  TH1D *  h_SL_fatjet_sdmass_NbNjCate	     [10]  ; 
  TH1D *  h_SL_fatjet_subjet_loosetagged_NbNjCate [10] ; 
  TH1D *  h_SL_fatjet_subjet_mediumtagged_NbNjCate[10] ; 
  TH1D *  h_SL_fatjet_subjet_mean_btagger_NbNjCate[10] ; 
  TH1D *  h_SL_N_fatjet_NbNjCate		     [10]; 
  TH1D *  h_SL_N_LooseBBTaggedfatjet_NbNjCate		     [10]; 
  TH1D *  h_SL_N_MediumBBTaggedfatjet_NbNjCate		     [10]; 


  TH1D *  h_SL_fatjet_pt_NbNjCate_SDmassCut                 [10] ; 
  TH1D *  h_SL_fatjet_eta_NbNjCate_SDmassCut		     [10]; 
  TH1D *  h_SL_fatjet_phi_NbNjCate_SDmassCut		     [10]; 
  TH1D *  h_SL_fatjet_tau21_NbNjCate_SDmassCut		     [10]; 
  TH1D *  h_SL_fatjet_tau32_NbNjCate_SDmassCut		     [10]; 
  TH1D *  h_SL_fatjet_subjet_loosetagged_NbNjCate_SDmassCut [10] ; 
  TH1D *  h_SL_fatjet_subjet_mediumtagged_NbNjCate_SDmassCut[10] ; 
  TH1D *  h_SL_fatjet_subjet_mean_btagger_NbNjCate_SDmassCut[10] ; 
  TH1D *  h_SL_N_fatjet_NbNjCate_SDmassCut		     [10]; 


  TH1D * h_SL_fatjet_Pt_TheFatJet_Cate6j4b_OneLooseBBtagFatJet      ;
  TH1D * h_SL_fatjet_Eta_TheFatJet_Cate6j4b_OneLooseBBtagFatJet     ;
  TH1D * h_SL_fatjet_AbsEta_TheFatJet_Cate6j4b_OneLooseBBtagFatJet     ;
  TH1D * h_SL_fatjet_SDMass_TheFatJet_Cate6j4b_OneLooseBBtagFatJet  ;
  TH1D * h_SL_fatjet_dRLep_TheFatJet_Cate6j4b_OneLooseBBtagFatJet   ;
  TH1D * h_SL_fatjet_dEtaLep_TheFatJet_Cate6j4b_OneLooseBBtagFatJet ;
  TH1D * h_SL_fatjet_dPhiLep_TheFatJet_Cate6j4b_OneLooseBBtagFatJet ;

  TH1D * h_SL_fatjet_Pt_TheFatJet_Cate6j4b_OneMediumBBtagFatJet      ;
  TH1D * h_SL_fatjet_Eta_TheFatJet_Cate6j4b_OneMediumBBtagFatJet     ;
  TH1D * h_SL_fatjet_AbsEta_TheFatJet_Cate6j4b_OneMediumBBtagFatJet     ;
  TH1D * h_SL_fatjet_SDMass_TheFatJet_Cate6j4b_OneMediumBBtagFatJet  ;
  TH1D * h_SL_fatjet_dRLep_TheFatJet_Cate6j4b_OneMediumBBtagFatJet   ;
  TH1D * h_SL_fatjet_dEtaLep_TheFatJet_Cate6j4b_OneMediumBBtagFatJet ;
  TH1D * h_SL_fatjet_dPhiLep_TheFatJet_Cate6j4b_OneMediumBBtagFatJet ;


  TH1D * h_SL_N_MediumBBTaggedfatjet_Cate5moreJ_2moreB;
  TH1D * h_SL_fatjet_Pt_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet      ;
  TH1D * h_SL_fatjet_Eta_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet     ;
  TH1D * h_SL_fatjet_AbsEta_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet     ;
  TH1D * h_SL_fatjet_SDMass_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet  ;
  TH1D * h_SL_fatjet_dRLep_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet   ;
  TH1D * h_SL_fatjet_dEtaLep_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet ;
  TH1D * h_SL_fatjet_dPhiLep_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet ;


  TH1D *  h_SL_fatjet_pt_NbNjCate_el                 [10] ; 
  TH1D *  h_SL_fatjet_eta_NbNjCate_el		     [10]; 
  TH1D *  h_SL_fatjet_phi_NbNjCate_el		     [10]; 
  TH1D *  h_SL_fatjet_tau21_NbNjCate_el		     [10]; 
  TH1D *  h_SL_fatjet_tau32_NbNjCate_el		     [10]; 
  TH1D *  h_SL_fatjet_sdmass_NbNjCate_el	     [10]  ; 
  TH1D *  h_SL_fatjet_subjet_loosetagged_NbNjCate_el [10] ; 
  TH1D *  h_SL_fatjet_subjet_mediumtagged_NbNjCate_el[10] ; 
  TH1D *  h_SL_fatjet_subjet_mean_btagger_NbNjCate_el[10] ; 
  TH1D *  h_SL_N_fatjet_NbNjCate_el		     [10]; 

  TH1D *  h_SL_fatjet_pt_NbNjCate_mu		     [10]; 
  TH1D *  h_SL_fatjet_eta_NbNjCate_mu		     [10]; 
  TH1D *  h_SL_fatjet_phi_NbNjCate_mu		     [10]; 
  TH1D *  h_SL_fatjet_tau21_NbNjCate_mu		     [10]; 
  TH1D *  h_SL_fatjet_tau32_NbNjCate_mu		     [10]; 
  TH1D *  h_SL_fatjet_sdmass_NbNjCate_mu	     [10]	     ; 
  TH1D *  h_SL_fatjet_subjet_loosetagged_NbNjCate_mu [10] ; 
  TH1D *  h_SL_fatjet_subjet_mediumtagged_NbNjCate_mu[10] ; 
  TH1D *  h_SL_fatjet_subjet_mean_btagger_NbNjCate_mu[10] ; 
  TH1D *  h_SL_N_fatjet_NbNjCate_mu                  [10] ; 

  int  syst_SDMassScale ; 
  int  syst_SDMassResolution ; 

  TH1D *  h_MinBTagDiscri_NoCut [2] ;
  TH1D *  h_DoubleB_NoCut [2] ;

  TH1D *  h_MinBTagDiscri_6j4b [2] ;
  TH1D *  h_DoubleB_6j4b [2] ;
  TH1D *  h_HandMergedBtag_6j4b [2] ;

  TH1D *  h_ak8BDT_6j4b [2] ;

  TH1D *  h_MinBTagDiscri_5j4b [2] ;
  TH1D *  h_DoubleB_5j4b [2] ;

  TH1D *  h_MinBTagDiscri_4j4b [2] ;
  TH1D *  h_DoubleB_4j4b [2] ;

  TH1D *  h_MinBTagDiscri_6j2b [2] ;
  TH1D *  h_DoubleB_6j2b [2] ;

  TH1D *  h_MinBTagDiscri_6j3b [2] ;
  TH1D *  h_DoubleB_6j3b [2] ;



  TH1D *  h_SL_MinBTagDiscri_NoCut [2] ;
  TH1D *  h_SL_DoubleB_NoCut [2] ;
  TH1D *  h_SL_MinBTagDiscri_6j4b [2] ;
  TH1D *  h_SL_DoubleB_6j4b [2] ;
  TH1D *  h_SL_MinBTagDiscri_5j4b [2] ;
  TH1D *  h_SL_DoubleB_5j4b [2] ;
  TH1D *  h_SL_MinBTagDiscri_4j4b [2] ;
  TH1D *  h_SL_DoubleB_4j4b [2] ;
  TH1D *  h_SL_MinBTagDiscri_6j2b [2] ;
  TH1D *  h_SL_DoubleB_6j2b [2] ;
  TH1D *  h_SL_MinBTagDiscri_6j3b [2] ;
  TH1D *  h_SL_DoubleB_6j3b [2] ;

  TH1D *  h_SL_HandMergedBtag_6j4b [2] ;


  TH2D * h2_MinB_DoubleB_noCut [2];
  TH2D * h2_MinB_DoubleB_6j4b [2];
  TH2D * h2_MinB_DoubleB_5j4b [2];
  TH2D * h2_MinB_DoubleB_4j4b [2];
  TH2D * h2_MinB_DoubleB_6j3b [2];
  TH2D * h2_MinB_DoubleB_6j2b [2];

  TH2D * h2_MinB_DoubleB_6j4b_twosubjet [2][2];

  TH1D * h_NFatJetPassingKinCriteria [10];


  long   eventnumber ; 
  long   eventcategory ; 
  float  eventweight ;
  long  ak8jet_genhiggsmatch ;
  float ak8jet_pt     ;
  float ak8jet_eta    ;
  float ak8jet_sdmass ;
  float ak8jet_deepcsv_min ; ;
  float ak8jet_doubleb ; 
  float ak8jet_tau21 ;
  float ak8jet_tau32 ;

  TH1D *   h_SLFatjetStudy_BDTScore_6j4b   ;
  TH1D *   h_SLFatjetStudy_BDTScore_6j4b_NoHiggsAk8Candidate    ;
  TH1D *   h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate     ;
  TH1D *   h_SLFatjetStudy_ak8BDTScore_6j4b    ;
  TH1D *   h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_highak8bdt    ;
  TH1D *   h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_lowak8bdt     ;

  TH1D *   h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_high0ak8bdt    ;
  TH1D *   h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_low0ak8bdt     ;

  TH1D *   h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_high2ak8bdt    ;
  TH1D *   h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_low2ak8bdt     ;

};

#endif
