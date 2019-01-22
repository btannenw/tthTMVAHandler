
#include "analyzer.h"

#include <iostream>
#include <math.h>
#include <algorithm>

#include <TLorentzVector.h>
#include <TRandom3.h>
#include <assert.h>

#include "reader.h"

#define MAX_N_NPTPLOT 100


analyzer::analyzer ()
 : output("out.root")
 , totalEvent ( 0 )
 , isMC( true )
 , isMuonStream( true )
 , n_th_event(0)
 , b_FakeEstimation(false)
 , h_FakeSF_Mu (0)
 , h_FakeSF_El (0)
 , ttbarAdditionalJetID ( noRequirement )
 , flag_skip_oddEventNumber( false )
 , syst_muTrigSF( 0 ) 
 , syst_elTrigSF( 0 ) 
 , syst_muIDSF  ( 0 ) 
 , syst_elIDSF  ( 0 )
 , syst_muIsoSF  ( 0 ) 
 , syst_elRecoSF  ( 0 )
 , syst_FakeMuStat( 0 )
 , syst_FakeElStat( 0 )
 , b_DileptonFakeLeptonAnalysisMode( false )
 , period ( -1 )
 , MCGenWeightKey ( ttHYggdrasilScaleFactors::DEFAULT )  
 , syst_SDMassScale (0)
 , syst_SDMassResolution (0)
{
}


void analyzer::setSFSystematic( int syst ){

  if( syst == 1 ){ syst_muTrigSF = + 1 ; };
  if( syst == 2 ){ syst_muTrigSF = - 1 ; };
  if( syst == 3 ){ syst_elTrigSF = + 1 ; };
  if( syst == 4 ){ syst_elTrigSF = - 1 ; };

  if( syst == 5 ){ syst_muIDSF = + 1 ; };
  if( syst == 6 ){ syst_muIDSF = - 1 ; };
  if( syst == 7 ){ syst_elIDSF = + 1 ; };
  if( syst == 8 ){ syst_elIDSF = - 1 ; };

  if( syst == 9  ){ syst_muIsoSF  =  +1 ; };
  if( syst == 10 ){ syst_muIsoSF  =  -1 ; };
  if( syst == 11 ){ syst_elRecoSF =  +1 ; };
  if( syst == 12 ){ syst_elRecoSF =  -1 ; };

  if( syst == 20 ){ syst_SDMassScale = +1  ;}
  if( syst == 21 ){ syst_SDMassScale = -1  ;}
  if( syst == 22 ){ syst_SDMassResolution = -1 ;}

  if( syst == 101 ){ syst_FakeMuStat = +1; }
  if( syst == 102 ){ syst_FakeMuStat = -1; }
  if( syst == 103 ){ syst_FakeElStat = +1; }
  if( syst == 104 ){ syst_FakeElStat = -1; }

  if( syst == 202 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::ME_muR_DOUBLE_muF_NOM     ; } 
  if( syst == 203 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::ME_muR_HALF_muF_NOM	     ; } 
  if( syst == 204 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::ME_muR_NOM_muF_DOUBLE     ; } 
  if( syst == 205 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::ME_muR_DOUBLE_muF_DOUBLE  ; } 
  if( syst == 206 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::ME_muR_HALF_muF_DOUBLE    ; } 
  if( syst == 207 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::ME_muR_NOM_muF_HALF	     ; } 
  if( syst == 208 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::ME_muR_DOUBLE_muF_HALF    ; } 
  if( syst == 209 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::ME_muR_HALF_muF_HALF      ; } 

  // PS uncertainty 
  if( syst == 301 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_REDUCED_ISRUP        ;}
  if( syst == 302 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_REDUCED_FSRUP        ;}
  if( syst == 303 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_REDUCED_ISRDOWN      ;}
  if( syst == 304 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_REDUCED_FSRDOWN      ;}
  if( syst == 305 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_DEFAULT_ISRUP        ;}
  if( syst == 306 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_DEFAULT_FSRUP        ;}
  if( syst == 307 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_DEFAULT_ISRDOWN      ;}
  if( syst == 308 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_DEFAULT_FSRDOWN      ;}
  if( syst == 309 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_CONSERVATIVE_ISRUP   ;} 
  if( syst == 310 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_CONSERVATIVE_FSRUP   ;} 
  if( syst == 311 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_CONSERVATIVE_ISRDOWN ;} 
  if( syst == 312 ){ MCGenWeightKey = ttHYggdrasilScaleFactors::PS_CONSERVATIVE_FSRDOWN ;} 

}


analyzer::~analyzer (){

}

void analyzer::AnalyzeEvent(){
  
  if( r -> EventNumber == 78430374 ){
    std::cout <<"[SPECIAL] the event 78430374 is ignored."  << std::endl ; 
    return ; 
  }


  n_th_event++ ;
  if( n_th_event % 10000 == 0 ){
    std::cout << "Event-"<<n_th_event << " is analyzed."<< std::endl ; 
  }
  

  // And we only use the events with odd event IDs for the training
  // (where this is the ID you can get from asking

  if( flag_skip_oddEventNumber && ( (r -> EventNumber ) % 2 == 1 ) ){
    return ; 
  }


  AddNumberOfProcessedEvent( 1.0 );

  npT = r -> TruePV ;
  npT = 
    ( npT  <  0 ) ?  0 :
    ( npT  > 49 ) ? 49 : npT ; 

  FillHistogram ( h_totalevent , 0 ,  isMC ? r->wgt_MCEventWeight : 1.0  );
  FillHistogram ( h_totalevent_noWeight , 0  , 1.0 );


//  // -- Overwrite scale factors.
//  double recalc_PUwgt = ttHSF . get_pu_wgt( r -> TruePV ) ;  
//  r -> wgt_PU = recalc_PUwgt ; 
//

  if ( syst_muTrigSF > 0 ){ r -> wgt_trigMu = r -> wgt_trigMuUP   ; }
  if ( syst_muTrigSF < 0 ){ r -> wgt_trigMu = r -> wgt_trigMuDOWN ; }
  if ( syst_elTrigSF > 0 ){ r -> wgt_trigEl = r -> wgt_trigElUP   ; }
  if ( syst_elTrigSF < 0 ){ r -> wgt_trigEl = r -> wgt_trigElDOWN ; }
  if ( syst_muIDSF   > 0 ){ r -> wgt_muID   = r -> wgt_muIDUP ; }
  if ( syst_muIDSF   < 0 ){ r -> wgt_muID   = r -> wgt_muIDDOWN ; }
  if ( syst_elIDSF   > 0 ){ r -> wgt_elID   = r -> wgt_elIDUP ; }
  if ( syst_elIDSF   < 0 ){ r -> wgt_elID   = r -> wgt_elIDDOWN ; }
  if ( syst_muIsoSF   > 0 ){ r -> wgt_muIso   = r -> wgt_muIsoUP ; }
  if ( syst_muIsoSF   < 0 ){ r -> wgt_muIso   = r -> wgt_muIsoDOWN ; }
  if ( syst_elRecoSF   > 0 ){ r -> wgt_elReco   = r -> wgt_elRecoUP ; }
  if ( syst_elRecoSF   < 0 ){ r -> wgt_elReco   = r -> wgt_elRecoDOWN ; }


  if( MCGenWeightKey != ttHYggdrasilScaleFactors::DEFAULT ){
    // Use the key to find the index in the weight vector.
    double sf = ( r -> mcWeight_value ->at(  MCGenWeightKey ) ) / ( r -> mcWeight_value ->at( ttHYggdrasilScaleFactors::ME_LHC_ORIGINAL_XWGTUP ) );
    r-> wgt_MCEventWeight *= sf ; 

  }


  // Overwrite event weight with period-dependent weight 

  if( period >= 0 ){
    r -> wgt_PU     = ( r -> periodwgt_PU    -> at( period ) ) ;
    r -> wgt_btag   = ( r -> periodwgt_btag  -> at( period ) ) ;
    r -> wgt_trigMu = ( r -> periodwgt_trigMu-> at( period ) ) ;
    r -> wgt_trigEl = ( r -> periodwgt_trigEl-> at( period ) ) ;
    r -> wgt_EGZtvx = ( r -> periodwgt_EGZtvx-> at( period ) ) ;
    r -> wgt_muID   = ( r -> periodwgt_muID  -> at( period ) ) ;
    r -> wgt_elID   = ( r -> periodwgt_elID  -> at( period ) ) ;
    r -> wgt_elReco = ( r -> periodwgt_elReco-> at( period ) ) ;
    r -> wgt_muIso  = ( r -> periodwgt_muIso -> at( period ) ) ;
  }

  // overwright event weights which are not ready for 2017-data analysis.
  r -> wgt_trigEl = 1 ; 
  r -> wgt_muID   = 1 ; 
  r -> wgt_muIso  = 1 ; 
  r -> wgt_topPT  = 1 ; 
  
  //- -- - - - - - - - - - - - - - - 

  r -> wgt_TOTAL = 
      r -> wgt_PU 
    * r -> wgt_btag  
    * r -> wgt_trigMu
    * r -> wgt_trigEl
    * r -> wgt_muID  
    * r -> wgt_elID  
    * r -> wgt_elReco
    * r -> wgt_muIso 
    * r -> wgt_MCEventWeight 
    * r -> wgt_topPT 
    ;


  FillHistogram ( h_wgt_TOTAL ,  r-> wgt_TOTAL  , 1.0 );
  FillHistogram ( h_wgt_PU    ,  r-> wgt_PU     , 1.0 );
  FillHistogram ( h_wgt_btag  ,  r-> wgt_btag   , 1.0 );
  FillHistogram ( h_wgt_trigMu,  r-> wgt_trigMu , 1.0 );
  FillHistogram ( h_wgt_trigEl,  r-> wgt_trigEl , 1.0 );
  FillHistogram ( h_wgt_muID  ,  r-> wgt_muID   , 1.0 );
  FillHistogram ( h_wgt_elID  ,  r-> wgt_elID   , 1.0 );
  FillHistogram ( h_wgt_elReco,  r-> wgt_elReco , 1.0 );
  FillHistogram ( h_wgt_muIso ,  r-> wgt_muIso  , 1.0 );
  FillHistogram ( h_wgt_MCEven,  r-> wgt_MCEventWeight , 1.0 );
  FillHistogram ( h_wgt_topPT ,  r-> wgt_topPT  , 1.0 );

  if(  ttbarAdditionalJetID  != noRequirement ){
    
    if( ! _passTtbarAdditionalJetID() ){ return ; }

  }


  FillHistogram ( h_pileup_noWeight ,  r -> TruePV    ,  1.0  );
  FillHistogram ( h_pileup_PUWeight ,  r -> TruePV    ,  r -> wgt_PU  );


#ifdef TRUTH_INFO_STUDY

  work_for_truthInfo();

#else
  
  work_for_singleLepton();
  
  work_for_fakeevaluation();


  //(*) : wgt_Total includes Z_vtx weight assuming the analysis is single lepton channel only.
  //      You need to modify the code to apply Z_vtx appropriately for di-lepton channel
  //       before you re-enable this part.
  //(*) work_for_fakeeval_Dilepton();

  //outtree   ->Fill();  

#endif

}


void analyzer::work_for_fakeeval_Dilepton(){

//  [1] di-lepton trigger,
//  [2] two tight leptons
//   (for fake, confirm that one isolated, one non-isolated.)
//  [3] lepton same sign,
//  [4] 4jets (to forcus on ttbar SL events with a fake lepton)
//

  const double wgt_total = ( isMC ? r-> wgt_TOTAL / r->wgt_trigMu / r->wgt_trigEl : 1.0 );

  FillHistogram( h_ee_SameSign_SanityCheck_NlooseLepMinusTight ,  r -> nLooseLep - r -> nTightLep , wgt_total ) ;
  //  FillHistogram( h_ee_SameSign_SanityCheck_NlooseLepMinusNonIsoTight ,  r -> nLooseLep - r -> nTightButNonIsoLep , wgt_total ) ;

  if( r -> pass_goodVtx != 1 ) return ; 
  if( r -> nLooseLep != 2 ) return ; 
  if( r -> nTightLep != 2 ) return ;

  // Same sign requirement.
  if( r -> lepton_charge != r -> lepton2_charge ) return ; 

  // In case of "b_DileptonFakeLeptonAnalysisMode", 
  //  - One lepton should be iso
  //  - The other should be non-iso
  if( b_DileptonFakeLeptonAnalysisMode ){
    //    if ( r -> nTightButNonIsoLep != 1 || r -> nTightLep != 2 ) return ; 
  }

  // -- Only di-electoon channel :
  if( r -> lepton_isMuon != 0 ) return ; 
  if( r -> lepton2_isMuon != 0 ) return ; 

  // matching btw trigger and reconstructed lepton
  if( r -> lepton_isMuon == 1 && r -> lepton2_isMuon == 1 && r -> pass_TrigMuMu != 1 ) return ; 
  if( r -> lepton_isMuon      != r -> lepton2_isMuon      && r -> pass_TrigElMu != 1 ) return ; 
  if( r -> lepton_isMuon == 0 && r -> lepton2_isMuon == 0 && r -> pass_TrigElEl != 1 ) return ;

  if( r -> nJet != 4 ) return ; 

  //  categorize :
//  [A] b-tagged jet multiplicity
//  [B] in/out Zmass window (check charge flip contamination)


// Histogram : 
// - Amount in each category.
// - Electron PT
// - Electron Eta

  

  const int nB = r->nBJet > 4 ? 4 : r->nBJet ; 

  FillHistogram( h_ee_SameSign_BtagjetMultiplicity , nB , wgt_total ) ;

  TLorentzVector l1, l2 ;
  l1 . SetPtEtaPhiM( r-> lepton_pt , r->lepton_eta, r->lepton_phi, r->lepton_m);
  l2 . SetPtEtaPhiM( r-> lepton2_pt , r->lepton2_eta, r->lepton2_phi, r->lepton2_m);
  
  const double mll = (l1+l2) .M();

  FillHistogram( h_ee_SameSign_mll , mll , wgt_total ) ;

  const int isInMllWindow = fabs ( mll - 91 ) < 10 ? 0 : 1 ;

  if(  isInMllWindow ==1 ){
    FillHistogram( h_ee_SameSign_BtagjetMultiplicity_MllCut , nB , wgt_total ) ;
  }

  FillHistogram( h_ee_SameSign_ZmassWindow  [isInMllWindow][ nB ] , 0.5 , wgt_total) ;
  FillHistogram( h_ee_SameSign_ElectronPT   [isInMllWindow][ nB ] , r -> lepton_pt  , wgt_total) ;
  FillHistogram( h_ee_SameSign_ElectronPT   [isInMllWindow][ nB ] , r -> lepton2_pt , wgt_total) ;
  FillHistogram( h_ee_SameSign_LeadingJetPT [isInMllWindow][ nB ] , r -> jet_pt->at(0)  , wgt_total) ;

  return ;
}


void analyzer::work_for_fakeevaluation(){


  bool minimal_requirement = 
    r-> pass_goodVtx == 1 
    &&
    ( r-> pass_TrigEl == 1 ||  r-> pass_TrigMu == 1 ) 
    &&
    // Matching of flavor of the DataStream and the trigger 
    ( isMC || (isMuonStream && (r-> pass_TrigMu == 1 ) ) || ( ! isMuonStream && r-> pass_TrigEl == 1 ) )
    && 
    ( r-> lepton_pt > 0 ) // at least one lepton
    ;

  if( ! minimal_requirement ) return ; 
  
  bool NeedToApply_EGamma_ZVtxWgt =
    isMC 
    &&
    r-> pass_TrigEl == 1 
    && 
    r-> pass_TrigMu == 0 ;



  ///bool single_electron_event =  r -> nLooseLep == 1 && r -> nTightLep == 1 &&  r -> lepton_isMuon == 0 ;
  const bool single_muon_event =  r -> nLooseLep == 1 && r -> nTightLep == 1 &&  r -> lepton_isMuon == 1 ;

  const double wgt_FakeSF = ( ! b_FakeEstimation ) ? 1.0 : 
    ( isMC ? -1.0 : 1.0 ) * _getFakeSF( single_muon_event ); // negative because NonIsolatedMCEvents are  to be subtracted from data.

  const double wgt_total_noFakeSF = ! isMC ? 1.0 : (
    (  r-> wgt_TOTAL  )
    * 
    ( NeedToApply_EGamma_ZVtxWgt ? r -> wgt_EGZtvx : 1.0 ) ) ;

  const double wgt_total = 
    wgt_total_noFakeSF
    *  wgt_FakeSF ;

  // GoodVtx and trigger were already required.

  const double MET_CUT_THRESHOLD = 20 ; 

  const int bin_nJ =  r->nJet > 6 ? 6 : r->nJet ; 
  const int bin_nB =  r->nBJet > 4 ? 4 : r->nBJet ; 

  if( 
     ( isMC || (! isMuonStream ) )
     &&
     r -> nLooseLep == 1 && r -> nTightLep == 1 &&  r -> lepton_isMuon == 0  && r-> pass_TrigEl == 1 && r->nJet != 0 ){
    
    FillHistogram( h_FakeEval_nJet_El_onejetcut_noFakeWgt   ,  r->nJet   , wgt_total_noFakeSF );
    FillHistogram( h_FakeEval_nBJet_El_onejetcut_noFakeWgt   ,  r->nBJet , wgt_total_noFakeSF );
    FillHistogram( h_FakeEval_Met_El_njnb_noFakeWgt    [ bin_nJ ][ bin_nB ] ,  r->met_abs , wgt_total_noFakeSF );
    if( r->met_abs < MET_CUT_THRESHOLD ){
      FillHistogram( h_FakeEval_NEntries_El_LowMet_noFakeWgt[ bin_nJ ][ bin_nB ] , 1 , wgt_total_noFakeSF );
      FillHistogram2D( h2_FakeEval_NEntries_El_nJnB_noFakeWgt , bin_nJ, bin_nB , wgt_total_noFakeSF );
    }

    FillHistogram( h_FakeEval_nJet_El_onejetcut_withFakeWgt   ,  r->nJet   , wgt_total );
    FillHistogram( h_FakeEval_nBJet_El_onejetcut_withFakeWgt   ,  r->nBJet , wgt_total );
    FillHistogram( h_FakeEval_Met_El_njnb_withFakeWgt    [ bin_nJ ][ bin_nB ] ,  r->met_abs , wgt_total );
    if( r->met_abs < MET_CUT_THRESHOLD ){
      FillHistogram( h_FakeEval_NEntries_El_LowMet_withFakeWgt[ bin_nJ ][ bin_nB ] , 1 , wgt_total );
      FillHistogram2D( h2_FakeEval_NEntries_El_nJnB_withFakeWgt , bin_nJ, bin_nB , wgt_total );
    }

    //For Fake systematic uncertainty : 
    if( 0 < r->met_abs && r->met_abs < MET_CUT_THRESHOLD / 2.0 ){ // low side 
      FillHistogram2D( h2_FakeEval_NEntries_LowMetSide_El_nJnB_noFakeWgt , bin_nJ, bin_nB , wgt_total_noFakeSF );
      FillHistogram2D( h2_FakeEval_NEntries_LowMetSide_El_nJnB_withFakeWgt , bin_nJ, bin_nB , wgt_total );
    }else  if( MET_CUT_THRESHOLD / 2.0 < r->met_abs && r->met_abs < MET_CUT_THRESHOLD ){ // highMet side 
      FillHistogram2D( h2_FakeEval_NEntries_HighMetSide_El_nJnB_noFakeWgt , bin_nJ, bin_nB , wgt_total_noFakeSF );
      FillHistogram2D( h2_FakeEval_NEntries_HighMetSide_El_nJnB_withFakeWgt , bin_nJ, bin_nB , wgt_total );
    }

    // For a study of fake-lepton 
    const double min_dr_jet_lep = miniDR( r -> lepton_eta , r -> jet_eta  , r -> lepton_phi , r -> jet_phi ) ;
    FillHistogram( h_FakeEval_dR_Jet_and_El_onejetcut_noFakeWgt  [ bin_nJ ][ bin_nB ] , min_dr_jet_lep  , wgt_total_noFakeSF );

    std::vector<float> b_eta;
    std::vector<float> b_phi;
    for( unsigned int i = 0 ; i <  r -> jet_phi -> size() ; i ++ ){
      if(   r -> jet_btag -> at( i ) < 0.8484 ) continue ; // btag cut 
      b_eta . push_back( r -> jet_eta ->at(i) );
      b_phi . push_back( r -> jet_phi ->at(i) );
    }
    const double min_dr_btag_lep = miniDR( r -> lepton_eta , & b_eta , r -> lepton_phi , & b_phi ) ;
    FillHistogram( h_FakeEval_dR_Btag_and_El_onejetcut_noFakeWgt [ bin_nJ ][ bin_nB ]  , min_dr_btag_lep , wgt_total_noFakeSF );
    
    if( bin_nJ == 1 && bin_nB == 0 ){
      FillHistogram( h_FakeEval_dR_Jet_and_El_onejetcut_noFakeWgt_1j0b , min_dr_jet_lep , wgt_total_noFakeSF );
    }
    if( bin_nJ == 1 && bin_nB == 1 ){
      FillHistogram( h_FakeEval_dR_Jet_and_El_onejetcut_noFakeWgt_1j1b , min_dr_jet_lep , wgt_total_noFakeSF );
    }
      

  }


  if( 
     ( isMC || (isMuonStream ) )
     &&
     r -> nLooseLep == 1 && r -> nTightLep == 1 &&  r -> lepton_isMuon == 1  && r-> pass_TrigMu == 1 && r->nJet != 0 ){
    
    FillHistogram( h_FakeEval_nJet_Mu_onejetcut_noFakeWgt   ,  r->nJet   , wgt_total_noFakeSF );
    FillHistogram( h_FakeEval_nBJet_Mu_onejetcut_noFakeWgt   ,  r->nBJet , wgt_total_noFakeSF );
    FillHistogram( h_FakeEval_Met_Mu_njnb_noFakeWgt    [ bin_nJ ][ bin_nB ] ,  r->met_abs , wgt_total_noFakeSF );
    if( r->met_abs < MET_CUT_THRESHOLD ){
      FillHistogram( h_FakeEval_NEntries_Mu_LowMet_noFakeWgt[ bin_nJ ][ bin_nB ] , 1 , wgt_total_noFakeSF );
      FillHistogram2D( h2_FakeEval_NEntries_Mu_nJnB_noFakeWgt , bin_nJ, bin_nB , wgt_total_noFakeSF );

      FillHistogram( h_FakeEval_NEntriesVsNPV_Mu_LowMet_noFakeWgt[ bin_nJ ][ bin_nB ] ,  r -> N_Vtx  , wgt_total_noFakeSF );
    }

    FillHistogram( h_FakeEval_nJet_Mu_onejetcut_withFakeWgt   ,  r->nJet   , wgt_total );
    FillHistogram( h_FakeEval_nBJet_Mu_onejetcut_withFakeWgt   ,  r->nBJet , wgt_total );
    FillHistogram( h_FakeEval_Met_Mu_njnb_withFakeWgt    [ bin_nJ ][ bin_nB ] ,  r->met_abs , wgt_total );
    if( r->met_abs < MET_CUT_THRESHOLD ){
      FillHistogram( h_FakeEval_NEntries_Mu_LowMet_withFakeWgt[ bin_nJ ][ bin_nB ] , 1 , wgt_total );
      FillHistogram2D( h2_FakeEval_NEntries_Mu_nJnB_withFakeWgt , bin_nJ, bin_nB , wgt_total );
    }

    //For Fake systematic uncertainty : 
    if( 0 < r->met_abs && r->met_abs < MET_CUT_THRESHOLD / 2.0 ){ // low side 
      FillHistogram2D( h2_FakeEval_NEntries_LowMetSide_Mu_nJnB_noFakeWgt , bin_nJ, bin_nB , wgt_total_noFakeSF );
      FillHistogram2D( h2_FakeEval_NEntries_LowMetSide_Mu_nJnB_withFakeWgt , bin_nJ, bin_nB , wgt_total );
    }else  if( MET_CUT_THRESHOLD / 2.0 < r->met_abs && r->met_abs < MET_CUT_THRESHOLD ){ // highMet side 
      FillHistogram2D( h2_FakeEval_NEntries_HighMetSide_Mu_nJnB_noFakeWgt , bin_nJ, bin_nB , wgt_total_noFakeSF );
      FillHistogram2D( h2_FakeEval_NEntries_HighMetSide_Mu_nJnB_withFakeWgt , bin_nJ, bin_nB , wgt_total );
    }

    // For a study of fake-lepton 
    const double min_dr_jet_lep = miniDR( r -> lepton_eta , r -> jet_eta  , r -> lepton_phi , r -> jet_phi ) ;
    FillHistogram( h_FakeEval_dR_Jet_and_Mu_onejetcut_noFakeWgt  [ bin_nJ ][ bin_nB ] , min_dr_jet_lep  , wgt_total_noFakeSF );

    std::vector<float> b_eta;
    std::vector<float> b_phi;
    for( unsigned int i = 0 ; i <  r -> jet_phi -> size() ; i ++ ){
      if(   r -> jet_btag -> at( i ) < 0.8484 ) continue ; // btag cut 
      b_eta . push_back( r -> jet_eta ->at(i) );
      b_phi . push_back( r -> jet_phi ->at(i) );
    }
    const double min_dr_btag_lep = miniDR( r -> lepton_eta , & b_eta , r -> lepton_phi , & b_phi ) ;
    FillHistogram( h_FakeEval_dR_Btag_and_Mu_onejetcut_noFakeWgt [ bin_nJ ][ bin_nB ]  , min_dr_btag_lep , wgt_total_noFakeSF );

    if( bin_nJ == 1 && bin_nB == 0 ){
      FillHistogram( h_FakeEval_dR_Jet_and_Mu_onejetcut_noFakeWgt_1j0b , min_dr_jet_lep , wgt_total_noFakeSF );
    }
    if( bin_nJ == 1 && bin_nB == 1 ){
      FillHistogram( h_FakeEval_dR_Jet_and_Mu_onejetcut_noFakeWgt_1j1b , min_dr_jet_lep , wgt_total_noFakeSF );
    }

  }



  const double wgt_nonIso = 
    b_FakeEstimation ? 0.0 // <- does not fill.
    :
    (
    isMC ? 
    r-> wgt_PU *
    r-> wgt_MCEventWeight *
    r-> wgt_trigMu * 
    r-> wgt_trigEl *
    r-> wgt_muID *
    //    r-> wgt_elID *
    r-> wgt_elReco *
    //r-> wgt_muIso 
    r->wgt_btag
    :
    1.0 
     ) ;

  const double wgt_iso = 
    wgt_FakeSF
    * 
    (
    isMC ? 
    r-> wgt_PU *
    r-> wgt_MCEventWeight *
    r-> wgt_trigMu * 
    r-> wgt_trigEl *
    r-> wgt_muID *
    r-> wgt_elID *
    r-> wgt_elReco *
    r-> wgt_muIso *
    r->wgt_btag
    :
    1.0 
     )
    ;



  if( 
     ( isMC || (! isMuonStream ) )
     &&
     r -> nLooseLep == 1 && r -> nNonIsoEl == 1 && r -> nNonIsoMu == 0 && r-> pass_TrigEl == 1 && r->nJet != 0 ){
    // memo : nNonIsoEl is inclusive : contains both isolated- and non-isolated-lepton.

    if( r-> nTightLep == 0 ){ // El fake control region
      
      FillHistogram( h_FakeEval_Met_NonIsoEl_njnb[ r->nJet > 6 ? 6 : r->nJet ][ r->nBJet > 4 ? 4 : r->nBJet ] ,  r->met_abs , wgt_nonIso);
      
      h2_FakeEval_NonisoEl_njetNbin       -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
      h2_FakeEval_NonisoEl_njetNbinCutOff -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
      if( r->met_abs < 20 ){ // to select QCD, require met to be small
	h2_FakeEval_NonisoEl_njetNbin_LowMet -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
	h2_FakeEval_NonisoEl_njetNbinCutOff_LowMet -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
      }else{
	h2_FakeEval_NonisoEl_njetNbin_notLowMet -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
	h2_FakeEval_NonisoEl_njetNbinCutOff_notLowMet -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
      }

    }else if ( r-> nTightLep == 1 && r->lepton_isMuon == 0 ){

      FillHistogram( h_FakeEval_Met_IsoEl_njnb[ r->nJet > 6 ? 6 : r->nJet ][ r->nBJet > 4 ? 4 : r->nBJet ] ,  r->met_abs , wgt_iso );

      h2_FakeEval_IsoEl_njetNbin       -> Fill(  r->nJet  , r->nBJet , wgt_iso);
      h2_FakeEval_IsoEl_njetNbinCutOff -> Fill(  r->nJet  , r->nBJet , wgt_iso);
      if( r->met_abs < 20 ){ // to select QCD, require met to be small  
	h2_FakeEval_IsoEl_njetNbin_LowMet -> Fill(  r->nJet  , r->nBJet , wgt_iso);
	h2_FakeEval_IsoEl_njetNbinCutOff_LowMet -> Fill(  r->nJet  , r->nBJet , wgt_iso);
      }else{
	h2_FakeEval_IsoEl_njetNbin_notLowMet -> Fill(  r->nJet  , r->nBJet , wgt_iso);
	h2_FakeEval_IsoEl_njetNbinCutOff_notLowMet -> Fill(  r->nJet  , r->nBJet , wgt_iso);
      }


    }

  }
  
  if(
     ( isMC || ( isMuonStream ) )
     &&
     r -> nLooseLep == 1 && r -> nNonIsoEl == 0 && r -> nNonIsoMu == 1 && r-> pass_TrigMu == 1 && r->nJet != 0 ){
    // memo : nNonIsoMu is inclusive : contains both isolated- and non-isolated-lepton.

    if( r-> nTightLep == 0 ){ // Mu fake control region
      
      FillHistogram( h_FakeEval_Met_NonIsoMu_njnb[ r->nJet > 6 ? 6 : r->nJet ][ r->nBJet > 4 ? 4 : r->nBJet ] ,  r->met_abs , wgt_nonIso );

      h2_FakeEval_NonisoMu_njetNbin       -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
      h2_FakeEval_NonisoMu_njetNbinCutOff -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
      if( r->met_abs < 20 ){
	h2_FakeEval_NonisoMu_njetNbin_LowMet -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
	h2_FakeEval_NonisoMu_njetNbinCutOff_LowMet -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
      }else{
	h2_FakeEval_NonisoMu_njetNbin_notLowMet -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
	h2_FakeEval_NonisoMu_njetNbinCutOff_notLowMet -> Fill(  r->nJet  , r->nBJet , wgt_nonIso);
      }


    }else if ( r-> nTightLep == 1 && r->lepton_isMuon == 1 ){

      FillHistogram( h_FakeEval_Met_IsoMu_njnb[ r->nJet > 6 ? 6 : r->nJet ][ r->nBJet > 4 ? 4 : r->nBJet ] ,  r->met_abs , wgt_iso );

      h2_FakeEval_IsoMu_njetNbin       -> Fill(  r->nJet  , r->nBJet , wgt_iso);
      h2_FakeEval_IsoMu_njetNbinCutOff -> Fill(  r->nJet  , r->nBJet , wgt_iso);
      if( r->met_abs < 20 ){
	h2_FakeEval_IsoMu_njetNbin_LowMet -> Fill(  r->nJet  , r->nBJet , wgt_iso);
	h2_FakeEval_IsoMu_njetNbinCutOff_LowMet -> Fill(  r->nJet  , r->nBJet , wgt_iso);
      }else{
	h2_FakeEval_IsoMu_njetNbin_notLowMet -> Fill(  r->nJet  , r->nBJet , wgt_iso);
	h2_FakeEval_IsoMu_njetNbinCutOff_notLowMet -> Fill(  r->nJet  , r->nBJet , wgt_iso);
      }

    }

  }
  

}

void analyzer::work_for_singleLepton(){


  bool minimal_requirement = 
    r-> pass_goodVtx == 1 
    &&
    ( r-> pass_TrigEl == 1 ||  r-> pass_TrigMu == 1 ) 
    &&
    // Matching of flavor of the DataStream and the trigger 
    ( isMC || (isMuonStream && (r-> pass_TrigMu == 1 ) ) || ( ! isMuonStream && r-> pass_TrigEl == 1 ) )
    && 
    ( r-> lepton_pt > 0 ) // at least one lepton
    ;

  if( ! minimal_requirement ) return ; 
  
  bool NeedToApply_EGamma_ZVtxWgt =
    isMC 
    &&
    r-> pass_TrigEl == 1 
    && 
    r-> pass_TrigMu == 0 ;

  const bool single_muon_event =  r -> nLooseLep == 1 && r -> nTightLep == 1 &&  r -> lepton_isMuon == 1 ;

  const double wgt_FakeSF = ( ! b_FakeEstimation ) ? 1.0 :
    ( isMC ? -1.0 : 1.0 ) * _getFakeSF( single_muon_event ); // negative because NonIsolatedMCEvents are  to be subtracted from data.                                                  


  const double wgt_total = ( isMC ? r-> wgt_TOTAL : 1.0 ) 
    * wgt_FakeSF 
    * ( NeedToApply_EGamma_ZVtxWgt ? r -> wgt_EGZtvx : 1.0 ) ;

  const double wgt_afterLepTrig = 
    wgt_FakeSF
    * 
    (
    isMC ? 
    r-> wgt_PU *
    r-> wgt_MCEventWeight *
    r-> wgt_trigMu * 
    r-> wgt_trigEl *
    ( NeedToApply_EGamma_ZVtxWgt ? r -> wgt_EGZtvx : 1.0 ) * 
    r-> wgt_muID *
    r-> wgt_elID *
    r-> wgt_elReco *
    r-> wgt_muIso :
    1.0 
     )
    ;
  // Point : no b-tag SF.

  FillHistogram( h_wgt_EGZvtx , NeedToApply_EGamma_ZVtxWgt ? r -> wgt_EGZtvx : 1.0 , 1.0 ) ;

  const double wgt_total_noPU = ( 
				 ! isMC ? 
				 1.0 :
				  r -> wgt_btag
				  * r -> wgt_trigMu
				  * r -> wgt_trigEl
				  * r -> wgt_muID
				  * r -> wgt_elID
				  * r -> wgt_elReco
				  * r -> wgt_muIso
				  * r -> wgt_MCEventWeight
				  * r -> wgt_topPT 
				  )
    * ( NeedToApply_EGamma_ZVtxWgt ? r -> wgt_EGZtvx : 1.0 ) 
    *  wgt_FakeSF ;

  const double wgt_el_related = 
    isMC ?
    r-> wgt_trigEl * r-> wgt_elID * r-> wgt_elReco 
    * ( NeedToApply_EGamma_ZVtxWgt ? r -> wgt_EGZtvx : 1.0 )
    : 1.0;

  const double wgt_mu_related = 
    isMC ?
    r-> wgt_trigMu * r-> wgt_muIso * r-> wgt_muID 
    : 1.0;

  FillHistogram( h_SL_N_Vtx_AfterSingleLepTrigAndTightLep  , r -> N_Vtx , wgt_afterLepTrig );
  FillHistogram( h_SL_N_VtxWITHOUTPU_AfterSingleLepTrigAndTightLep  , r -> N_Vtx , wgt_afterLepTrig / r->wgt_PU);


  // Jet kinematics check after lepton reuirement 
  {
    // [1] electron
    bool elchannel =  r-> pass_goodVtx == 1 && r -> pass_TrigEl == 1 && r-> nTightLep == 1 && r->nLooseLep == 1 && r->lepton_isMuon == 0  && ( isMC || ! isMuonStream );
    if( elchannel ){
      FillHistogram( h_SL_nJet_el_AllButJetRequirement , r->nJet , wgt_afterLepTrig );
      FillHistogram( h_SL_ElPT_AllButJetRequirement    , r-> lepton_pt , wgt_afterLepTrig ) ;
      FillHistogram( h_SL_Met_AllButJetRequirement_el  , r -> met_abs , wgt_afterLepTrig ) ;
      FillHistogram( h_SL_MetType1xy_AllButJetRequirement_el  , r -> met_type1xy_abs , wgt_afterLepTrig ) ;
    }
    // [2] muon
    bool muchannel =  r-> pass_goodVtx == 1 && r -> pass_TrigMu == 1 && r-> nTightLep == 1 && r->nLooseLep == 1 && r->lepton_isMuon == 1  && ( isMC || isMuonStream );
    if( muchannel ){
      FillHistogram( h_SL_nJet_mu_AllButJetRequirement , r->nJet , wgt_afterLepTrig );
      FillHistogram( h_SL_MuonPT_AllButJetRequirement  , r-> lepton_pt , wgt_afterLepTrig ) ;
      FillHistogram( h_SL_Met_AllButJetRequirement_mu  , r -> met_abs  , wgt_afterLepTrig ) ;
      FillHistogram( h_SL_MetType1xy_AllButJetRequirement_mu  , r -> met_type1xy_abs  , wgt_afterLepTrig ) ;


      {
	long nJet  = r->nJet  < 6 ? r-> nJet : 6 ;
	long nBJet = r->nBJet < 4 ? r->nBJet : 4 ;

	FillHistogram( h_SL_Met_nJnBCategoryAfterLeptonCut_Mu [nJet] [nBJet], 
		       r-> met_abs , wgt_total ) ; // here the evet weight includes b-tagging weight as I use nBJet.

	if( r-> met_abs > 20 ) {
	  FillHistogram( h_SL_MuonPT_nJnBCategoryAfterLeptonCut_MetCut [nJet][nBJet] , 
			 r-> lepton_pt , wgt_total ) ; // here the evet weight includes b-tagging weight as I use nBJet.

	  const long bin = ( nJet + 1 )* nJet / 2 + ( nBJet + 1 ) ;
	  FillHistogram( h_SL_Entries_nJnBCategoryAfterLeptonCut_MetCut_Mu , bin , wgt_total );

	}

      }



      if(  r->nJet == 0 ) FillHistogram( h_SL_MET_AllButJetRequirement_jet0_mu  , r -> met_abs  , wgt_afterLepTrig ) ;
      if(  r->nJet == 1 ) FillHistogram( h_SL_MET_AllButJetRequirement_jet1_mu  , r -> met_abs  , wgt_afterLepTrig ) ;
      if(  r->nJet == 2 ) FillHistogram( h_SL_MET_AllButJetRequirement_jet2_mu  , r -> met_abs  , wgt_afterLepTrig ) ;
      if(  r->nJet == 3 ) FillHistogram( h_SL_MET_AllButJetRequirement_jet3_mu  , r -> met_abs  , wgt_afterLepTrig ) ;
      if(  r->nJet == 4 ) FillHistogram( h_SL_MET_AllButJetRequirement_jet4_mu  , r -> met_abs  , wgt_afterLepTrig ) ;
      if(  r->nJet  > 4 ) FillHistogram( h_SL_MET_AllButJetRequirement_jet5_mu  , r -> met_abs  , wgt_afterLepTrig ) ;

      FillHistogram2D( h2_SL_NPV_NJet_AfterLeptonEventSelection_mu , r -> N_Vtx , r->nJet , wgt_afterLepTrig ) ;  

    }
    if( elchannel || muchannel ){
      FillHistogram( h_SL_nJet_AllButJetRequirement , r->nJet , wgt_afterLepTrig );

      if( r->nJet >= 4 ){ // scope for jot 

	if( elchannel ){
	  FillHistogram( h_SL_ElPT_AllButNBJetRequirement   , r-> lepton_pt , wgt_afterLepTrig ) ;
	  FillHistogram( h_SL_Met_AllButNBJetRequirement_el , r -> met_abs  , wgt_afterLepTrig ) ;
	}else{
	  FillHistogram( h_SL_MuonPT_AllButNBJetRequirement , r-> lepton_pt , wgt_afterLepTrig ) ;
	  FillHistogram( h_SL_Met_AllButNBJetRequirement_mu , r -> met_abs  , wgt_afterLepTrig ) ;
	}

	for( int i = 0 ; i < 8 && i < r->nJet ; i ++ ){
	  FillHistogram ( h_SL_jet_pt_AfterAllButBtagJetRequirement  [i], r-> jet_pt ->at(i)   , wgt_afterLepTrig ); 
	  FillHistogram ( h_SL_jetBtag_AfterAllButBtagJetRequirement [i], r-> jet_btag->at(i)  , wgt_afterLepTrig );
	}
	
      }// scope of jet : closed      
      
    }

  }//end of scope 

  /// e-mu channel
  if( ( ( r -> pass_TrigEl == 1 && (isMC || !isMuonStream) )
	||
	( r -> pass_TrigMu == 1 && (isMC ||  isMuonStream) ) ) // trigger-stream matching in case of data.
      && r-> nTightLep == 2 && r->nLooseLep == 2  // exact two leptons
      && ( r->lepton_isMuon != r->lepton2_isMuon  ) // Opposite flav
      && ( r->lepton_charge != r->lepton2_charge  ) // Opposite charge
      ){
    
    // To avoid event overlap in MuonDataStream and ElectronDataStream,
    //  decide that only events MuTrig=fail && ElTrig=pass is allow in ElectronStrema.
    if( isMC || 
	isMuonStream ||
	( r -> pass_TrigMu == 0 && r -> pass_TrigEl == 1 ) )
      {
     
    TLorentzVector v1,v2;
    v1 . SetPtEtaPhiM( r -> lepton_pt,  r -> lepton_eta , r -> lepton_phi , 0 );
    v2 . SetPtEtaPhiM( r -> lepton2_pt, r -> lepton2_eta, r -> lepton2_phi, 0 );
    const double mass = (v1 + v2 ).M();
    if(  r-> nJet == 0  &&  mass > 60 && r -> met_abs > 30  ){
      
      FillHistogram( h_SL_LeptopnPT_WWemuEnrichRegionLoose , r-> lepton_pt  , wgt_afterLepTrig );
      FillHistogram( h_SL_LeptopnPT_WWemuEnrichRegionLoose , r-> lepton2_pt , wgt_afterLepTrig );
      FillHistogram( h_MET_WWemuEnrichRegionLoose , r->  met_abs , wgt_afterLepTrig );
      if( mass > 90 && r -> met_abs > 40  ){
	FillHistogram( h_SL_LeptopnPT_WWemuEnrichRegionTight , r-> lepton_pt  , wgt_afterLepTrig );
	FillHistogram( h_SL_LeptopnPT_WWemuEnrichRegionTight , r-> lepton2_pt , wgt_afterLepTrig );
	FillHistogram( h_MET_WWemuEnrichRegionTight , r->  met_abs , wgt_afterLepTrig );
      }
      
    }
    }// end : avoiding overlap in MuonDataStream and ElectronDataStream

  }



  /// e-mu same sign
  if( ( ( r -> pass_TrigEl == 1 && (isMC || !isMuonStream ) )
	||
	( r -> pass_TrigMu == 1 && (isMC ||  isMuonStream ) ) )
      && r-> nTightLep == 2 && r->nLooseLep == 2  // exact two leptons
      && ( r->lepton_isMuon != r->lepton2_isMuon  ) // Opposite flav
      && ( r->lepton_charge == r->lepton2_charge  ) // same charge
      ){

    // To avoid event overlap in MuonDataStream and ElectronDataStream,
    //  decide that only events MuTrig=fail && ElTrig=pass is allow in ElectronStrema.
    if( isMC || 
	isMuonStream ||
	( r -> pass_TrigMu == 0 && r -> pass_TrigEl == 1 ) )
      {
    
//    TLorentzVector v1,v2;
//    v1 . SetPtEtaPhiM( r -> lepton_pt,  r -> lepton_eta , r -> lepton_phi , 0 );
//    v2 . SetPtEtaPhiM( r -> lepton2_pt, r -> lepton2_eta, r -> lepton2_phi, 0 );
//    const double mass = (v1 + v2 ).M();

    if(  r->lepton_isMuon == 1 ){ 
      FillHistogram( h_SL_MuonPT_WWemuSameSignEnrichRegionLoose,     r-> lepton_pt  , wgt_afterLepTrig );
      FillHistogram( h_SL_ElectronPT_WWemuSameSignEnrichRegionLoose, r-> lepton2_pt  , wgt_afterLepTrig );
    }else{
      FillHistogram( h_SL_MuonPT_WWemuSameSignEnrichRegionLoose,     r-> lepton2_pt  , wgt_afterLepTrig );
      FillHistogram( h_SL_ElectronPT_WWemuSameSignEnrichRegionLoose, r-> lepton_pt  , wgt_afterLepTrig );
    }
    FillHistogram( h_MET_WWemuSameSignEnrichRegionLoose  , r-> met_abs , wgt_afterLepTrig );
    FillHistogram( h_nJet_WWemuSameSignEnrichRegionLoose , r-> nJet    , wgt_afterLepTrig );
    FillHistogram( h_nBJet_WWemuSameSignEnrichRegionLoose , r-> nBJet   , wgt_total );

    }// end : avoiding overlap in MuonDataStream and ElectronDataStream

  }


  // Two electron channels.
  if( r -> pass_TrigEl == 1 && r-> nTightLep == 2 && r->nLooseLep == 2 && r->lepton_isMuon == 0 && r-> Mll > 0 
      && ( isMC || ! isMuonStream ) // datastream-trigger macthing.
      && ( r->lepton_isMuon == r->lepton2_isMuon  ) // same flav
      && ( r->lepton_charge != r->lepton2_charge  ) // Opposite charge	
      ){
    FillHistogram ( h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons ,  r-> Mll  , wgt_afterLepTrig );
    FillHistogram ( h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_noElSF ,  r-> Mll  , wgt_afterLepTrig / wgt_el_related );

    FillHistogram ( h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElTrigSF ,  r-> Mll  , wgt_afterLepTrig / wgt_el_related * r-> wgt_trigEl );
    FillHistogram ( h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElRecoSF ,  r-> Mll  , wgt_afterLepTrig / wgt_el_related * r-> wgt_elReco );
    FillHistogram ( h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElIDSF   ,  r-> Mll  , wgt_afterLepTrig / wgt_el_related * r-> wgt_elID   );

    if( r-> Mll > 130 && r-> nJet == 0 && r -> met_abs > 50 ){
      // ww control region 
      FillHistogram( h_SL_El1_PT_WWenrichRegion , r-> lepton_pt  , wgt_afterLepTrig );
      FillHistogram( h_SL_El2_PT_WWenrichRegion , r-> lepton2_pt , wgt_afterLepTrig );
    }

    if( fabs( r-> Mll - 91.0 ) < 10 ){
      FillHistogram( h_SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut , r-> lepton_pt  , wgt_afterLepTrig );
      FillHistogram( h_SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut , r-> lepton2_pt , wgt_afterLepTrig );

      FillHistogram( h_SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut_noElSF , r-> lepton_pt  , wgt_afterLepTrig / wgt_el_related );
      FillHistogram( h_SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut_noElSF , r-> lepton2_pt , wgt_afterLepTrig / wgt_el_related );

      FillHistogram( h_SL_nJet_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut , r-> nJet , wgt_afterLepTrig );
      
      if( r-> nJet > 0 ){
	FillHistogram( h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut , r-> jet_pt -> at(0) , wgt_afterLepTrig );
      }
      
      if( r-> nJet > 3 ){
	FillHistogram( h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut , r-> jet_pt -> at(0) , wgt_afterLepTrig );
	FillHistogram( h_SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut , r-> jet_pt -> at(1) , wgt_afterLepTrig );
	FillHistogram( h_SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut , r-> jet_pt -> at(2) , wgt_afterLepTrig );
	FillHistogram( h_SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut , r-> jet_pt -> at(3) , wgt_afterLepTrig );
      }

    }

  }

  // Two muon  channels.
  if( r -> pass_TrigMu == 1 && r-> nTightLep == 2 && r->nLooseLep == 2 && r->lepton_isMuon == 1 && r-> Mll > 0 
      && ( isMC || isMuonStream ) // datastream-trigger macthing.
      && ( r->lepton_isMuon == r->lepton2_isMuon  ) // Same flav 
      && ( r->lepton_charge != r->lepton2_charge  ) // Opposite charge	
      ){
    FillHistogram ( h_SL_Mmm_AfterSingleTriggerAndExactTwoTightleptons ,  r-> Mll  , wgt_afterLepTrig );

    FillHistogram ( h_SL_Mmm_AfterSingleTriggerAndExactTwoTightleptons_noMuSF ,  r-> Mll  , wgt_afterLepTrig / wgt_mu_related );

    if( r-> Mll > 130 && r-> nJet == 0 && r -> met_abs > 50 ){
      // ww control region 
      FillHistogram( h_SL_Mu1_PT_WWenrichRegion , r-> lepton_pt  , wgt_afterLepTrig );
      FillHistogram( h_SL_Mu2_PT_WWenrichRegion , r-> lepton2_pt , wgt_afterLepTrig );
    }

    if( fabs( r-> Mll - 91.0 ) < 10 ){
      FillHistogram( h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut , r-> lepton_pt  , wgt_afterLepTrig );
      FillHistogram( h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut , r-> lepton2_pt , wgt_afterLepTrig );

      FillHistogram( h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuSF , r-> lepton_pt  , wgt_afterLepTrig / wgt_mu_related);
      FillHistogram( h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuSF , r-> lepton2_pt , wgt_afterLepTrig / wgt_mu_related);

      FillHistogram( h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuTrigSF , r-> lepton_pt  , wgt_afterLepTrig / r->wgt_trigMu );
      FillHistogram( h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuTrigSF , r-> lepton2_pt , wgt_afterLepTrig / r->wgt_trigMu );


      FillHistogram( h_SL_nJet_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut , r-> nJet , wgt_afterLepTrig );

      FillHistogram( h_SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut , r-> met_abs , wgt_afterLepTrig );
      FillHistogram( h_SL_MetPhi_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut , r-> met_phi , wgt_afterLepTrig );
      
      if( r-> nJet > 0 ){
	FillHistogram( h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut , r-> jet_pt ->at(0) , wgt_afterLepTrig );
	if( r-> nJet > 1 ){
	  FillHistogram( h_SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut , r-> jet_pt ->at(1) , wgt_afterLepTrig );
	  if( r-> nJet > 2 ){
	    FillHistogram( h_SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut , r-> jet_pt ->at(2) , wgt_afterLepTrig );
	    if( r-> nJet > 3 ){
	      FillHistogram( h_SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut , r-> jet_pt ->at(3) , wgt_afterLepTrig );
	      if( r-> nJet > 4 ){
		FillHistogram( h_SL_JetPt_5_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut , r-> jet_pt ->at(4) , wgt_afterLepTrig );
	      }
	    }
	  }
	}
      }


      if( r -> nJet == 0 )FillHistogram( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet0     , r -> met_abs , wgt_afterLepTrig );
      if( r -> nJet == 1 )FillHistogram( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet1     , r -> met_abs , wgt_afterLepTrig );
      if( r -> nJet == 2 )FillHistogram( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet2     , r -> met_abs , wgt_afterLepTrig );
      if( r -> nJet == 3 )FillHistogram( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet3     , r -> met_abs , wgt_afterLepTrig );
      if( r -> nJet == 4 )FillHistogram( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet4     , r -> met_abs , wgt_afterLepTrig );
      if( r -> nJet  > 4 )FillHistogram( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet5more , r -> met_abs , wgt_afterLepTrig );


      if( r -> nJet == 0 )FillHistogram( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet0     , r -> Mll , wgt_afterLepTrig );
      if( r -> nJet == 1 )FillHistogram( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet1     , r -> Mll , wgt_afterLepTrig );
      if( r -> nJet == 2 )FillHistogram( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet2     , r -> Mll , wgt_afterLepTrig );
      if( r -> nJet == 3 )FillHistogram( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet3     , r -> Mll , wgt_afterLepTrig );
      if( r -> nJet == 4 )FillHistogram( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet4     , r -> Mll , wgt_afterLepTrig );
      if( r -> nJet  > 4 )FillHistogram( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet5more , r -> Mll , wgt_afterLepTrig );
      
      if( r-> nJet > 3 ){
	FillHistogram( h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut , r-> jet_pt -> at(0) , wgt_afterLepTrig );
	FillHistogram( h_SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut , r-> jet_pt -> at(1) , wgt_afterLepTrig );
	FillHistogram( h_SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut , r-> jet_pt -> at(2) , wgt_afterLepTrig );
	FillHistogram( h_SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut , r-> jet_pt -> at(3) , wgt_afterLepTrig );
      }
      
    }
  }



  /// Charge flip study
  if( r -> pass_TrigMu == 1 && r-> nTightLep == 2 && r->nLooseLep == 2 && r->lepton_isMuon == 1 
      && ( isMC || isMuonStream ) // datastream-trigger macthing.
      && ( r->lepton_isMuon == r->lepton2_isMuon  ) // Same flav 
      && ( r->lepton_charge == r->lepton2_charge  ) // Same flav
      ){

    TLorentzVector v1,v2;
    v1 . SetPtEtaPhiM( r -> lepton_pt,  r -> lepton_eta , r -> lepton_phi , 0 );
    v2 . SetPtEtaPhiM( r -> lepton2_pt, r -> lepton2_eta, r -> lepton2_phi, 0 );
    const double mass = (v1 + v2 ).M();

    FillHistogram ( h_SL_InvMass_AfterMuMuSameSignSelection ,  mass  , wgt_afterLepTrig );

    if( fabs( mass - 91.0  ) < 10 ){
      FillHistogram ( h_SL_MuonPT__AfterMuMuSameSignSelection_ZmassCut ,   r -> lepton_pt  , wgt_afterLepTrig );
      FillHistogram ( h_SL_MuonEta_AfterMuMuSameSignSelection_ZmassCut ,   r -> lepton_eta , wgt_afterLepTrig );
      FillHistogram ( h_SL_MuonPT__AfterMuMuSameSignSelection_ZmassCut ,   r -> lepton2_pt  , wgt_afterLepTrig );
      FillHistogram ( h_SL_MuonEta_AfterMuMuSameSignSelection_ZmassCut ,   r -> lepton2_eta , wgt_afterLepTrig );
    }

  }
  /// Charge flip study
  if( r -> pass_TrigEl == 1 && r-> nTightLep == 2 && r->nLooseLep == 2 && r->lepton_isMuon == 0 
      && ( isMC || ! isMuonStream ) // datastream-trigger macthing.
      && ( r->lepton_isMuon == r->lepton2_isMuon  ) // Same flav 
      && ( r->lepton_charge == r->lepton2_charge  ) // Same flav
      ){

    TLorentzVector v1,v2;
    v1 . SetPtEtaPhiM( r -> lepton_pt,  r -> lepton_eta , r -> lepton_phi , 0 );
    v2 . SetPtEtaPhiM( r -> lepton2_pt, r -> lepton2_eta, r -> lepton2_phi, 0 );
    const double mass = (v1 + v2 ).M();

    FillHistogram ( h_SL_InvMass_AfterElElSameSignSelection ,  mass  , wgt_afterLepTrig );

    if( fabs( mass - 91.0  ) < 10 ){
      FillHistogram ( h_SL_ElectronPT__AfterElElSameSignSelection_ZmassCut ,   r -> lepton_pt  , wgt_afterLepTrig );
      FillHistogram ( h_SL_ElectronEta_AfterElElSameSignSelection_ZmassCut ,   r -> lepton_eta , wgt_afterLepTrig );
      FillHistogram ( h_SL_ElectronPT__AfterElElSameSignSelection_ZmassCut ,   r -> lepton2_pt  , wgt_afterLepTrig );
      FillHistogram ( h_SL_ElectronEta_AfterElElSameSignSelection_ZmassCut ,   r -> lepton2_eta , wgt_afterLepTrig );
    }

  }




  { 

  // Check btag scale factor 
  // [1] LF entich, di-lepton region
  if( r -> pass_TrigMu == 1 && r-> nTightLep == 2 && r->nLooseLep == 2 && r->lepton_isMuon == 1 && r->lepton2_isMuon == 1 
      && ( isMC || isMuonStream ) // datastream-trigger macthing.
      && ( r->lepton_charge != r->lepton2_charge  ) // Opposite charge	
      && std::fabs( r-> Mll -91 ) < 10 
      && r-> nJet >= 2 ){
    
    // Compared to HeavyFlav enrich region below, we do not have to care TrigELSF because I require two of the leptons to be muons. No ElTriggerSF is in the total weight.

    for( int i = 0 ; i < 8 && i < r->nJet ; i ++ ){
      FillHistogram ( h_SL_jetPt_upto8th_BtagLFControlRegion   , r->jet_pt -> at(i)  , wgt_afterLepTrig ); 
      FillHistogram ( h_SL_jetBtag_upto8th_BtagLFControlRegion , r->jet_btag -> at(i) , wgt_total        ); // <- In order to see effect of btagSF, use weight with btagSF.
      FillHistogram ( h_SL_jetBtag_upto8th_BtagLFControlRegionWOBtagSF , r->jet_btag -> at(i) , wgt_total / r->wgt_btag );
    }
    
  }
  if( r -> pass_TrigMu == 1 && r-> nTightLep == 2 && r->nLooseLep == 2 && r->lepton_isMuon == 1 && r->lepton2_isMuon == 1 
      && ( isMC || isMuonStream ) // datastream-trigger macthing.
      && ( r->lepton_charge != r->lepton2_charge  ) // Opposite charge	
      && ( r-> Mll - 91 ) > 10  //(only take high side)
      && r-> nJet >= 2 ){
    for( int i = 0 ; i < 8  && i < r->nJet ; i ++ ){
      FillHistogram ( h_SL_jetPt_upto8th_BtagHFControlRegion   ,  r->jet_pt -> at(i)   , wgt_afterLepTrig ); 
      FillHistogram ( h_SL_jetBtag_upto8th_BtagHFControlRegion ,  r->jet_btag ->at(i) , wgt_total       ); // <- In order to see effect of btagSF, use weight with btagSF.
      FillHistogram ( h_SL_jetBtag_upto8th_BtagHFControlRegionWOBtagSF , r->jet_btag ->at(i) , wgt_total / r->wgt_btag );
    }     
 
  }
  }// end of scope.

  
  
  // - - Fill Cut flow histogram 
  if( r-> pass_goodVtx == 1 && ( r-> pass_TrigEl == 1 || r-> pass_TrigMu ) && r-> nTightLep == 1 ) {
    // = Good Vtx, lepton trigger, nTight=1.
    
    if( isMC 
	|| 
	( r->lepton_isMuon == 0  && ! isMuonStream )  // or or offline electron in electron data stream 
	||
	( r->lepton_isMuon == 1  &&   isMuonStream )  // or offline electron in muon data stream 
	){

      const int idx_lepflav = r->lepton_isMuon ; 
	
      int CutflowStep = 0 ; 
      FillHistogram( h_SL_cutflow[idx_lepflav] , CutflowStep , wgt_afterLepTrig );

      if(  r->nLooseLep == 1 ){
	FillHistogram( h_SL_cutflow[idx_lepflav] , ++ CutflowStep , wgt_afterLepTrig );
	if( r->nJet > 0 ){
	  FillHistogram( h_SL_cutflow[idx_lepflav] , ++ CutflowStep , wgt_afterLepTrig );
	  if( r->nJet > 1 ){
	    FillHistogram( h_SL_cutflow[idx_lepflav] , ++ CutflowStep , wgt_afterLepTrig );
	    if( r->nJet > 2 ){
	      FillHistogram( h_SL_cutflow[idx_lepflav] , ++ CutflowStep , wgt_afterLepTrig );
	      if( r->nJet > 3 ){
		FillHistogram( h_SL_cutflow[idx_lepflav] , ++ CutflowStep , wgt_afterLepTrig );
		if( r -> nBJet > 0 ){
		  FillHistogram( h_SL_cutflow[idx_lepflav] , ++ CutflowStep , wgt_total );
		  if( r -> nBJet > 1 ){
		    FillHistogram( h_SL_cutflow[idx_lepflav] , ++ CutflowStep , wgt_total );
		  }
		}
	      }
	    }
	  }
	}

      }

    }
  } // End of filling cut-flow histogram 
  


  // ** require full event selection
  if( r-> pass_El != 1 &&  r-> pass_Mu != 1 ) return ;
  if( ! isMC && ! isMuonStream && r -> lepton_isMuon != 0 ) return ; // datastream-flavour macthing.
  if( ! isMC &&   isMuonStream && r -> lepton_isMuon != 1 ) return ; // datastream-flavour macthing.

  double bdt_score = -10 ; 
  #ifdef ENABLE_COMMON_CLASSFIER
  {
    std::vector<TLorentzVector> selectedLeptonP4 ; 
    std::vector<TLorentzVector> selectedJetP4 ; 
    std::vector<double> selectedJetCSV ;
    //    std::vector<TLorentzVector> looseSelectedJetP4 ; 
    //    std::vector<double> looseSelectedJetCSV ; 
    TLorentzVector metP4;

    {
      TLorentzVector v;
      v. SetPtEtaPhiM( r->lepton_pt, r->lepton_eta, r->lepton_phi, r->lepton_m );
      selectedLeptonP4.push_back( v );
    }
    for( int j = 0 ; j < r->nJet ; j++ ){
      TLorentzVector v;
      v. SetPtEtaPhiM( r->jet_pt ->at(j), r->jet_eta->at(j), r->jet_phi->at(j), r->jet_m->at(j) );
      selectedJetP4 .push_back( v );
      selectedJetCSV. push_back( r -> jet_btag ->at(j) );
    }
    
    metP4 . SetPtEtaPhiM( r->met_abs , 0,  r->met_phi , 0 );
    bdt_score = bdt.GetBDTOutput(
				 selectedLeptonP4, 
				 selectedJetP4 ,
				 selectedJetCSV , 
				 metP4 );

    //    std::cout <<"BDT_Score " << bdt_score << std::endl ; 
    const long cate = _EventCateBasedOnNjetNBtagJet(  r -> nJet , r -> nBJet );

    if(                       r-> pass_Mu == 1 ) FillHistogram( h_SL_eventCate_Mu_NoWeight_noMetCut , cate , 1.0 );
    if(                       r-> pass_El == 1 ) FillHistogram( h_SL_eventCate_El_NoWeight_noMetCut , cate , 1.0 );
    if(  r -> met_abs > 20 && r-> pass_Mu == 1 ) FillHistogram( h_SL_eventCate_Mu_NoWeight          , cate , 1.0 );
    if(  r -> met_abs > 20 && r-> pass_El == 1 ) FillHistogram( h_SL_eventCate_El_NoWeight	       , cate , 1.0 );
    

    if(  r -> met_abs > 20 && r-> pass_El == 1 ){ FillHistogram(  h_SL_ElPtAfterCate_MetCut [cate]   ,  r->lepton_pt , wgt_total ) ; }
    if(  r -> met_abs > 20 && r-> pass_Mu == 1 ){ FillHistogram(  h_SL_MuonPtAfterCate_MetCut [cate] ,  r->lepton_pt , wgt_total ) ; }

    if( isMC || b_FakeEstimation ){// This is blinding data. FakeData is processed.
      FillHistogram(  h_SL_BDT[cate] , bdt_score , wgt_total ) ; 

      // before met cut 
      if(   r-> pass_El == 1 ){ FillHistogram(  h_SL_BDT_nometcut_el [cate] , bdt_score , wgt_total ) ; }
      if(   r-> pass_Mu == 1 ){ FillHistogram(  h_SL_BDT_nometcut_mu [cate] , bdt_score , wgt_total ) ; }


      // after met cut 
      if(  r -> met_abs > 20 ){ FillHistogram(  h_SL_BDT_metcut [cate] , bdt_score , wgt_total ) ; }

      if(  r -> met_abs > 20 && r-> pass_El == 1 ){ FillHistogram(  h_SL_BDT_metcut_el [cate] , bdt_score , wgt_total ) ; }
      if(  r -> met_abs > 20 && r-> pass_Mu == 1 ){ FillHistogram(  h_SL_BDT_metcut_mu [cate] , bdt_score , wgt_total ) ; }
    }
  }
  #endif 


  { 
    // 
    const int N_rank = 4 ;
    double thre_jet[6][N_rank]= { { 80, 150, 250, 1000000 }
				  , { 60, 110, 200, 1000000 }
				  , { 40,  60, 100, 1000000 }
				  , { 40,  60, 100, 1000000 }
				  , { 0, 35 , 40, 1000000 } // Take into account "nJet<5."
				  , { 0, 35 , 40, 1000000 }};
    
    int rank_of_jet_pt[6]={0,0,0,0,0,0};
    for( unsigned int jet = 0 ; jet < 6 && jet < r -> jet_pt -> size() ; jet ++ ){
      for( int iRank = 1 ; iRank < 4 ; iRank ++ ){
	if(    r -> jet_pt ->at( jet ) > thre_jet[jet][ iRank - 1] 
	    && r -> jet_pt ->at( jet ) < thre_jet[jet][ iRank ] ){
	  rank_of_jet_pt [jet] = iRank ; 
	}
      }
    }

    long total_rank_4to6 = 0 ;
    for( unsigned int jet = 0 ; jet < 6 && jet < r -> jet_pt -> size() ; jet ++ ){
      total_rank_4to6 += rank_of_jet_pt[jet] * pow( 4 , jet );
    }

    long total_rank_3to6 = 0 ;
    for( unsigned int jet = 0 ; jet < 6 && jet < r -> jet_pt -> size() ; jet ++ ){
      total_rank_3to6 += ( rank_of_jet_pt[jet] == 3 ? 2 : rank_of_jet_pt[jet] ) * pow( 3 , jet );
    }

    long total_rank_2to4and3to2 = 0 ;
    for( unsigned int jet = 0 ; jet < 4 && jet < r -> jet_pt -> size() ; jet ++ ){ /// 1,2,3,4 jet
      total_rank_2to4and3to2 += ( rank_of_jet_pt[jet] > 0 ? 1 : 0 ) * pow( 2 , jet );
    }
    for( unsigned int jet = 4 ; jet < 5 && jet < r -> jet_pt -> size() ; jet ++ ){// 5 and 6 jet
      total_rank_2to4and3to2 +=  ( rank_of_jet_pt[jet] == 3 ? 2 : rank_of_jet_pt[jet] ) * pow( 3 , jet );
    }

    if( r->lepton_pt > 100 ){
      total_rank_4to6        +=  pow( 4, 6 );
      total_rank_3to6        +=  pow( 3, 6 );
      total_rank_2to4and3to2 +=  pow( 2, 4 ) * pow( 3, 2 );
    }

    if( false ){
      std::cout <<"Rank " 
		<< rank_of_jet_pt [0] << " " 
		<< rank_of_jet_pt [1] << " " 
		<< rank_of_jet_pt [2] << " " 
		<< rank_of_jet_pt [3] << " " 
		<< rank_of_jet_pt [4] << " " 
		<< rank_of_jet_pt [5] << std::endl;
      std::cout << "kepton = " <<  r->lepton_pt << std::endl ; 
      std::cout << "total rank " << total_rank_4to6  << std::endl ; 
      std::cout << "           " << total_rank_3to6  << std::endl ; 
      std::cout << "           " << total_rank_2to4and3to2  << std::endl ; 
    }

    if( r -> lepton_isMuon == 0 ){
     FillHistogram( h_SL_Elch_JetPtCategorization4to6        , total_rank_4to6        , wgt_total );
     FillHistogram( h_SL_Elch_JetPtCategorization3to6        , total_rank_3to6        , wgt_total );
     FillHistogram( h_SL_Elch_JetPtCategorization2to4and3to2 , total_rank_2to4and3to2 , wgt_total );
    }else{
     FillHistogram( h_SL_Much_JetPtCategorization4to6        , total_rank_4to6        , wgt_total );
     FillHistogram( h_SL_Much_JetPtCategorization3to6        , total_rank_3to6        , wgt_total );
     FillHistogram( h_SL_Much_JetPtCategorization2to4and3to2 , total_rank_2to4and3to2 , wgt_total );
    }

      
  }



  const long cate = _EventCateBasedOnNjetNBtagJet(  r -> nJet , r -> nBJet );
  {
    if( cate == 0 ){ std::cout <<"[problem] Event category based on nJet/nBjet was zero, which is not supposed to be after the event selection. Need to check..." << std::endl; }
    FillHistogram( h_SL_EventCategorizationBasedonNjetNBtagJet , cate , wgt_total ) ; 

    if(  r -> met_abs > 20 ){ FillHistogram( h_SL_EventCategorizationBasedonNjetNBtagJet_withMetCut , cate , wgt_total ) ; }//todo

    if(  r-> pass_Mu  ){
      FillHistogram( h_SL_EventCategorizationBasedonNjetNBtagJet_Mu , cate , wgt_total ) ; 
    }

    if(  r-> pass_El  ){
      FillHistogram( h_SL_EventCategorizationBasedonNjetNBtagJet_El , cate , wgt_total ) ; 
    }

    if(  r-> pass_Mu && r -> met_abs < 20 ){
      FillHistogram( h_SL_EventCategorizationBasedonNjetNBtagJet_FakeCR_Mu , cate , wgt_total ) ; 
    }

    if(  r-> pass_El && r -> met_abs < 20 ){
      FillHistogram( h_SL_EventCategorizationBasedonNjetNBtagJet_FakeCR_El , cate , wgt_total ) ; 
    }

    if( r -> met_abs > 20 ){

      int njet_category = r -> nJet < 6 ?  r -> nJet : 6 ;

      std::vector<long> idxFatJet ; 
      std::vector<long> idxLooseBBtaggedFatJet ; 
      std::vector<long> idxMediumBBtaggedFatJet ; 
      int nFatJet_SDmasscut = 0 ; 
      
      for(unsigned int iFat = 0 ;  iFat < r->fatjet_pt->size() ; iFat ++ ){

	// eta cut 
	if( fabs(  r->fatjet_eta ->at( iFat ) ) > 2.0 ) continue;  
	
	// dR(lepton) cut
	const double dR2 = 
	  _calcDR2( r->fatjet_eta ->at( iFat ),
		    r->lepton_eta , 
		    r->fatjet_phi ->at( iFat ),
		    r->lepton_phi ) ; 
	if( dR2 < 0.8 * 0.8 ) continue ; // (overlapping with leptons)

	idxFatJet . push_back( iFat ) ;  

	if( r -> fatjet_subjet_nloosebtag  ->at(iFat ) == 2 ){ idxLooseBBtaggedFatJet . push_back( iFat ) ; } ;
	if( r -> fatjet_subjet_nmediumbtag ->at(iFat ) == 2 ){ idxMediumBBtaggedFatJet. push_back( iFat ) ; } ;

	FillHistogram( h_SL_fatjet_pt_NbNjCate	[cate], r -> fatjet_pt     ->at(iFat)     , wgt_total ) ;  
	FillHistogram( h_SL_fatjet_eta_NbNjCate	[cate], r -> fatjet_eta    ->at(iFat)     , wgt_total ) ;  
	FillHistogram( h_SL_fatjet_phi_NbNjCate	[cate], r -> fatjet_phi    ->at(iFat)     , wgt_total ) ;  
	FillHistogram( h_SL_fatjet_sdmass_NbNjCate	[cate], r -> fatjet_sdmass ->at(iFat)     , wgt_total ) ;  
	FillHistogram( h_SL_fatjet_tau21_NbNjCate		    [cate], r -> fatjet_tau21  ->at(iFat)     , wgt_total ) ;  
	FillHistogram( h_SL_fatjet_tau32_NbNjCate		    [cate], r -> fatjet_tau32  ->at(iFat)     , wgt_total ) ;  
	FillHistogram( h_SL_fatjet_subjet_loosetagged_NbNjCate   [cate], r -> fatjet_subjet_nloosebtag  ->at(iFat), wgt_total ) ;  
	FillHistogram( h_SL_fatjet_subjet_mediumtagged_NbNjCate  [cate], r -> fatjet_subjet_nmediumbtag ->at(iFat), wgt_total ) ;  
	FillHistogram( h_SL_fatjet_subjet_mean_btagger_NbNjCate  [cate], r -> fatjet_subjet_meanbtagger ->at(iFat), wgt_total ) ;  
	
	if( r -> fatjet_sdmass ->at(iFat) > 40 ){ 
	  nFatJet_SDmasscut ++; 
	  FillHistogram( h_SL_fatjet_pt_NbNjCate_SDmassCut	[cate], r -> fatjet_pt     ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_eta_NbNjCate_SDmassCut	[cate], r -> fatjet_eta    ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_phi_NbNjCate_SDmassCut	[cate], r -> fatjet_phi    ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_tau21_NbNjCate_SDmassCut      [cate], r -> fatjet_tau21  ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_tau32_NbNjCate_SDmassCut      [cate], r -> fatjet_tau32  ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_subjet_loosetagged_NbNjCate_SDmassCut   [cate], r -> fatjet_subjet_nloosebtag  ->at(iFat), wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_subjet_mediumtagged_NbNjCate_SDmassCut  [cate], r -> fatjet_subjet_nmediumbtag ->at(iFat), wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_subjet_mean_btagger_NbNjCate_SDmassCut  [cate], r -> fatjet_subjet_meanbtagger ->at(iFat), wgt_total ) ;  
	} // end if SD mass cut.

      } // Fat jet loop ends for el/mu combined

      FillHistogram( h_SL_N_fatjet_NbNjCate           [cate], idxFatJet . size()  , wgt_total ) ;  
      FillHistogram( h_SL_N_fatjet_NbNjCate_SDmassCut [cate], nFatJet_SDmasscut   , wgt_total ) ;  

      FillHistogram( h_SL_N_LooseBBTaggedfatjet_NbNjCate      [cate], idxLooseBBtaggedFatJet . size()  , wgt_total ) ;  
      FillHistogram( h_SL_N_MediumBBTaggedfatjet_NbNjCate      [cate], idxMediumBBtaggedFatJet . size()  , wgt_total ) ;  


      // cate == 9 : 6j4b.
      if( cate == 9 && idxFatJet.size()  == 2  ){

	int idx_largeSDmassJet = idxFatJet[0] ; 
	int idx_smallSDmassJet = idxFatJet[1] ; 
	if( r-> fatjet_sdmass ->at( idx_largeSDmassJet ) < r-> fatjet_sdmass ->at( idx_smallSDmassJet ) ){
	  idx_largeSDmassJet = idxFatJet[1] ; 
	  idx_smallSDmassJet = idxFatJet[0] ; 
	}

	FillHistogram( h_SL_fatjet_SDMass_ofLargeSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_sdmass ->at( idx_largeSDmassJet)     , wgt_total ) ;
	FillHistogram( h_SL_fatjet_SDMass_ofSmallSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_sdmass ->at( idx_smallSDmassJet)     , wgt_total ) ;

	FillHistogram( h_SL_fatjet_pt_ofLargeSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_pt ->at( idx_largeSDmassJet)     , wgt_total ) ;
	FillHistogram( h_SL_fatjet_pt_ofSmallSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_pt ->at( idx_smallSDmassJet)     , wgt_total ) ;

	FillHistogram( h_SL_fatjet_eta_ofLargeSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_eta ->at( idx_largeSDmassJet)     , wgt_total ) ;
	FillHistogram( h_SL_fatjet_eta_ofSmallSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_eta ->at( idx_smallSDmassJet)     , wgt_total ) ;


	FillHistogram( h_SL_fatjet_tau32_ofLargeSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_tau32 ->at( idx_largeSDmassJet)     , wgt_total ) ;
	FillHistogram( h_SL_fatjet_tau32_ofSmallSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_tau32 ->at( idx_smallSDmassJet)     , wgt_total ) ;

	FillHistogram( h_SL_fatjet_tau21_ofLargeSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_tau21 ->at( idx_largeSDmassJet)     , wgt_total ) ;
	FillHistogram( h_SL_fatjet_tau21_ofSmallSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_tau21 ->at( idx_smallSDmassJet)     , wgt_total ) ;

	FillHistogram( h_SL_fatjet_NlooseB_ofLargeSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_subjet_nloosebtag ->at( idx_largeSDmassJet)     , wgt_total ) ;
	FillHistogram( h_SL_fatjet_NlooseB_ofSmallSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_subjet_nloosebtag ->at( idx_smallSDmassJet)     , wgt_total ) ;

	FillHistogram( h_SL_fatjet_NMediumB_ofLargeSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_subjet_nmediumbtag ->at( idx_largeSDmassJet)     , wgt_total ) ;
	FillHistogram( h_SL_fatjet_NMediumB_ofSmallSDmassJet_Cate6j4b_TwoFatJets , r -> fatjet_subjet_nmediumbtag ->at( idx_smallSDmassJet)     , wgt_total ) ;


	{
	  long int fatjet_cutflow_step = 0 ; 
	  FillHistogram( h_SL_fatjet_AmbitiousSelection_Cate6j4b_TwoFatJets , fatjet_cutflow_step ++ , wgt_total );
	  if( r -> fatjet_sdmass ->at( idx_largeSDmassJet) > 20 ){
	    FillHistogram( h_SL_fatjet_AmbitiousSelection_Cate6j4b_TwoFatJets , fatjet_cutflow_step ++ , wgt_total );
	    if( r -> fatjet_sdmass ->at( idx_smallSDmassJet) > 20 ){
	      FillHistogram( h_SL_fatjet_AmbitiousSelection_Cate6j4b_TwoFatJets , fatjet_cutflow_step ++ , wgt_total );
	      if( r -> fatjet_subjet_nloosebtag ->at( idx_smallSDmassJet) >= 2 ){
		FillHistogram( h_SL_fatjet_AmbitiousSelection_Cate6j4b_TwoFatJets , fatjet_cutflow_step ++ , wgt_total );
		if( r -> fatjet_subjet_nloosebtag ->at( idx_largeSDmassJet) < 2 ){ // 0 or 1
		  FillHistogram( h_SL_fatjet_AmbitiousSelection_Cate6j4b_TwoFatJets , fatjet_cutflow_step ++ , wgt_total );
		}
	      }
	    }
	  }
	}// end of scope.

	//todo
	{
	  const double dR2 = 
	    _calcDR2( r->fatjet_eta ->at( idx_largeSDmassJet),
		      r->lepton_eta , 
		      r->fatjet_phi ->at( idx_largeSDmassJet),
		      r->lepton_phi ) ; 
	  const double dR_with_lep = sqrt( dR2 ) ; 
	  FillHistogram( h_SL_fatjet_dRlep_ofLargeSDmassJet_Cate6j4b_TwoFatJets ,dR_with_lep  , wgt_total ) ;
	}
	{
	  const double dR2 = 
	    _calcDR2( r->fatjet_eta ->at( idx_smallSDmassJet),
		      r->lepton_eta , 
		      r->fatjet_phi ->at( idx_smallSDmassJet),
		      r->lepton_phi ) ; 
	  const double dR_with_lep = sqrt( dR2 ) ; 
	  FillHistogram( h_SL_fatjet_dRlep_ofSmallSDmassJet_Cate6j4b_TwoFatJets ,dR_with_lep  , wgt_total ) ;
	}

      }// end of 6j4b && FatJet==2.

      const bool b_ExactOnelooseBBFatJet = (idxLooseBBtaggedFatJet.size() == 1);
      bool SDmass_in_HiggsWindow_looseBBFatJet = false; 

      if(  b_ExactOnelooseBBFatJet ){
	unsigned long idx = idxLooseBBtaggedFatJet[0];

	{
	  SDmass_in_HiggsWindow_looseBBFatJet = fabs( r -> fatjet_sdmass ->at( idx ) - 120 ) < 20 ;
	}

	if( cate == 9 ){
	
	FillHistogram( h_SL_fatjet_Pt_TheFatJet_Cate6j4b_OneLooseBBtagFatJet    ,  r -> fatjet_pt     ->at( idx ) , wgt_total ) ;
	FillHistogram( h_SL_fatjet_Eta_TheFatJet_Cate6j4b_OneLooseBBtagFatJet   ,  r -> fatjet_eta    ->at( idx ) , wgt_total ) ;
	FillHistogram( h_SL_fatjet_AbsEta_TheFatJet_Cate6j4b_OneLooseBBtagFatJet,fabs(r -> fatjet_eta ->at( idx )), wgt_total ) ;
	FillHistogram( h_SL_fatjet_SDMass_TheFatJet_Cate6j4b_OneLooseBBtagFatJet,  r -> fatjet_sdmass ->at( idx ) , wgt_total ) ;

	{
	  const double dR2 = 
	    _calcDR2( r->fatjet_eta ->at( idx ),
		      r->lepton_eta , 
		      r->fatjet_phi ->at( idx ),
		      r->lepton_phi ) ; 
	  const double dR_with_lep = sqrt( dR2 ) ; 
	  FillHistogram( h_SL_fatjet_dRLep_TheFatJet_Cate6j4b_OneLooseBBtagFatJet ,dR_with_lep  , wgt_total ) ;
	}
	{
	  const double abs_dEta = fabs( r->fatjet_eta ->at( idx ) -  r->lepton_eta  ) ; 
	  FillHistogram( h_SL_fatjet_dEtaLep_TheFatJet_Cate6j4b_OneLooseBBtagFatJet , abs_dEta  , wgt_total ) ;
	}
	{
	  double d_phi = fabs( r->fatjet_phi ->at( idx ) -    r->lepton_phi ) ;
	  d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ;
	  FillHistogram( h_SL_fatjet_dPhiLep_TheFatJet_Cate6j4b_OneLooseBBtagFatJet , fabs( d_phi ) , wgt_total ) ;
	  
	}


      }// end 6j4b && LooseBBFatJet==1
      } // end if looseBB fat jet == 1 

      {  // Fill BDT histogram using information of SD mass of loose-BB fatjet 
	  if(  r -> met_abs > 20 ){


	    // Divide events based on existence of bb-tagged fat jet 
	    // 
	    //  [A]+[B]=total.
	    //  [A] = [a21] + [a-2].
	    // so, those histogram can be used to evaluate analysis performance in combination of 
	    //  * [A] + [B]
	    //  or
	    //  * [a-1] + [a-2] + [B]
	    //
	    if( b_ExactOnelooseBBFatJet ){
	      FillHistogram(  h_SL_fatjet_BDT_ExactOneBBFatJetPass [cate] , bdt_score   , wgt_total ) ;  // --[A]
	      if( SDmass_in_HiggsWindow_looseBBFatJet ){
		FillHistogram(  h_SL_fatjet_BDT_ExactOneBBFatJetPass_SDMassPass [cate] , bdt_score   , wgt_total ) ;  //[a-1]
	      }else{
		FillHistogram(  h_SL_fatjet_BDT_ExactOneBBFatJetPass_SDMassFail [cate] , bdt_score   , wgt_total ) ;  //[a-2]
	      }

	    }else{
	      FillHistogram(  h_SL_fatjet_BDT_ExactOneBBFatjetFail [cate] , bdt_score   , wgt_total ) ; // --[B]
	    }


	    const double bdt_score_truncatedOverflow = (   bdt_score < -0.8 
							   ?  - 0.799 
							   : ( bdt_score > +0.8 ? 0.799 : bdt_score )
							   );

	    const double bdt_score_with_offset = 
	      ! SDmass_in_HiggsWindow_looseBBFatJet  
	      ? bdt_score_truncatedOverflow
	      : ( bdt_score_truncatedOverflow < 0 ?  0.9 : 1.1 ) ; 

	    
	    if( SDmass_in_HiggsWindow_looseBBFatJet){
	      FillHistogram(  h_SL_fatjet_BDT_metcut_SDMass_in_HiggsWindow [cate] , bdt_score                  , wgt_total ) ; 
	    }else{
	      FillHistogram(  h_SL_fatjet_BDT_metcut_SDMass_out_HiggsWindow [cate] , bdt_score                 , wgt_total ) ; 
	    }
	    FillHistogram(  h_SL_fatjet_IsSDMassInHiggsWindow [cate]                  , SDmass_in_HiggsWindow_looseBBFatJet ? 1 : 0 , wgt_total ) ; 
	    FillHistogram(  h_SL_fatjet_BDT_metcut_withInfoSDMassinHiggsWindow [cate] , bdt_score_with_offset                       , wgt_total ) ; 

	  }
      }


      if(  cate == 9 && idxMediumBBtaggedFatJet.size() == 1 ){
	unsigned long idx = idxMediumBBtaggedFatJet[0];
	
	FillHistogram( h_SL_fatjet_Pt_TheFatJet_Cate6j4b_OneMediumBBtagFatJet    ,  r -> fatjet_pt     ->at( idx ) , wgt_total ) ;
	FillHistogram( h_SL_fatjet_Eta_TheFatJet_Cate6j4b_OneMediumBBtagFatJet   ,  r -> fatjet_eta    ->at( idx ) , wgt_total ) ;
	FillHistogram( h_SL_fatjet_AbsEta_TheFatJet_Cate6j4b_OneMediumBBtagFatJet,fabs(r -> fatjet_eta ->at( idx )), wgt_total ) ;
	FillHistogram( h_SL_fatjet_SDMass_TheFatJet_Cate6j4b_OneMediumBBtagFatJet,  r -> fatjet_sdmass ->at( idx ) , wgt_total ) ;

	{
	  const double dR2 = 
	    _calcDR2( r->fatjet_eta ->at( idx ),
		      r->lepton_eta , 
		      r->fatjet_phi ->at( idx ),
		      r->lepton_phi ) ; 
	  const double dR_with_lep = sqrt( dR2 ) ; 
	  FillHistogram( h_SL_fatjet_dRLep_TheFatJet_Cate6j4b_OneMediumBBtagFatJet ,dR_with_lep  , wgt_total ) ;
	}
	{
	  const double abs_dEta = fabs( r->fatjet_eta ->at( idx ) -  r->lepton_eta  ) ; 
	  FillHistogram( h_SL_fatjet_dEtaLep_TheFatJet_Cate6j4b_OneMediumBBtagFatJet , abs_dEta  , wgt_total ) ;
	}
	{
	  double d_phi = fabs( r->fatjet_phi ->at( idx ) -    r->lepton_phi ) ;
	  d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ;
	  FillHistogram( h_SL_fatjet_dPhiLep_TheFatJet_Cate6j4b_OneMediumBBtagFatJet , fabs( d_phi ) , wgt_total ) ;
	  
	}
      }// end 6j4b && MediumBBFatJet==1



      if( cate == 9 &&  r -> met_abs > 20 ){ // Study 6j4b-ak8-bdt performance

	int   N_candidate_ak8_jets = 0 ; 
	float maximum_ak8bdtscore = -10 ;
	
	for( unsigned int iFat = 0 ; iFat < ( r -> fatjet_pt ->size() ) ; iFat ++ ){
	  
	  if( ( r ->  fatjet_pt -> at(iFat ) )  < 250 ){ continue ; }
	  if( fabs ( r ->  fatjet_eta -> at(iFat ) )  > 2.0 ){ continue ; }
	  const double dR2 =
	    _calcDR2( r->fatjet_eta ->at( iFat ),
		      r->lepton_eta ,
		      r->fatjet_phi ->at( iFat ),
		      r->lepton_phi ) ;
	  if( dR2 < 0.8 * 0.8 ) continue ; // (overlapping with leptons)          

	  // -> Pass all conditions on ak8 jets.
	  N_candidate_ak8_jets ++ ; 

	  // Variables for MVA(BDT) INPUT 
	  tmva_ak8jet_pt          = r -> fatjet_pt              -> at(iFat);
	  tmva_ak8jet_eta         = r -> fatjet_eta             -> at(iFat) ;
	  tmva_ak8jet_sdmass      = r -> fatjet_sdmass          -> at(iFat) ;
	  {
	    double minb = ( r -> fatjet_subjet_minbtagger  -> at(iFat) )  ;
	    minb = minb < 0 ? 0 : minb ;
	    tmva_ak8jet_deepcsv_min =  minb ;
	  }
	  tmva_ak8jet_doubleb     = r -> fatjet_doubleb  -> at(iFat) ;
	  const double my_ak8bdt_value = tmva_reader -> EvaluateMVA("SatoshiAk8BDT");
	  
	  maximum_ak8bdtscore = maximum_ak8bdtscore > my_ak8bdt_value ? maximum_ak8bdtscore : my_ak8bdt_value ; 

	}// looof of fat jet ends.


	FillHistogram(  h_SLFatjetStudy_BDTScore_6j4b  , bdt_score   , wgt_total ) ;

	if( N_candidate_ak8_jets == 0 ){
	  // no candidate aj8 jet
	  
	  FillHistogram(  h_SLFatjetStudy_BDTScore_6j4b_NoHiggsAk8Candidate  , bdt_score   , wgt_total ) ;

	}else{
	  // at least one candidate ak8 jet.

	  FillHistogram(  h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate   , bdt_score   , wgt_total ) ; // BDT of events with good ak8.
	  FillHistogram(  h_SLFatjetStudy_ak8BDTScore_6j4b  , maximum_ak8bdtscore   , wgt_total ) ;// Score itlself

	  if( maximum_ak8bdtscore > 0.05 ){
	    // boosted higgs like event 
	    FillHistogram(  h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_highak8bdt   , bdt_score   , wgt_total ) ; 
	  }else{
	    // nonboosted higgs like event
	    FillHistogram(  h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_lowak8bdt   , bdt_score   , wgt_total ) ; 
	  }


	  if( maximum_ak8bdtscore > 0.1 ){
	    // boosted higgs like event 
	    FillHistogram(  h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_high2ak8bdt   , bdt_score   , wgt_total ) ; 
	  }else{
	    // nonboosted higgs like event
	    FillHistogram(  h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_low2ak8bdt   , bdt_score   , wgt_total ) ; 
	  }

	  if( maximum_ak8bdtscore > 0.0 ){
	    // boosted higgs like event 
	    FillHistogram(  h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_high0ak8bdt   , bdt_score   , wgt_total ) ; 
	  }else{
	    // nonboosted higgs like event
	    FillHistogram(  h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_low0ak8bdt   , bdt_score   , wgt_total ) ; 
	  }

	}

      }


      bool b_Event_5orMoreJet_2orMoreB = ( r -> nJet >= 5 ) && ( r -> nBJet >= 2 ) && (cate != 9 ) ; // exclude 6j4b cate.
      if( b_Event_5orMoreJet_2orMoreB ){

	FillHistogram( h_SL_N_MediumBBTaggedfatjet_Cate5moreJ_2moreB , idxMediumBBtaggedFatJet . size()  , wgt_total ) ;

	if( idxMediumBBtaggedFatJet.size() == 1 ){

	  unsigned long idx = idxMediumBBtaggedFatJet[0];

	  FillHistogram( h_SL_fatjet_Pt_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet    ,  r -> fatjet_pt     ->at( idx ) , wgt_total ) ;
	  FillHistogram( h_SL_fatjet_Eta_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet   ,  r -> fatjet_eta    ->at( idx ) , wgt_total ) ;
	  FillHistogram( h_SL_fatjet_AbsEta_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet,fabs(r -> fatjet_eta ->at( idx )), wgt_total ) ;
	  FillHistogram( h_SL_fatjet_SDMass_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet,  r -> fatjet_sdmass ->at( idx ) , wgt_total ) ;

	  {
	    const double dR2 = 
	      _calcDR2( r->fatjet_eta ->at( idx ),
			r->lepton_eta , 
			r->fatjet_phi ->at( idx ),
			r->lepton_phi ) ; 
	    const double dR_with_lep = sqrt( dR2 ) ; 
	    FillHistogram( h_SL_fatjet_dRLep_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet ,dR_with_lep  , wgt_total ) ;
	  }
	  {
	    const double abs_dEta = fabs( r->fatjet_eta ->at( idx ) -  r->lepton_eta  ) ; 
	    FillHistogram( h_SL_fatjet_dEtaLep_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet , abs_dEta  , wgt_total ) ;
	  }
	  {
	    double d_phi = fabs( r->fatjet_phi ->at( idx ) -    r->lepton_phi ) ;
	    d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ;
	    FillHistogram( h_SL_fatjet_dPhiLep_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet , fabs( d_phi ) , wgt_total ) ;
	  
	  }

	} // end if one fat BB-tagged jet.

      }// end 5 more jet, 2 more b-jet
	


      if( r->pass_El == 1 ){
	FillHistogram( h_SL_ElPt_AfterNjNbCategory[cate]    ,  r-> lepton_pt  , wgt_total );
	FillHistogram( h_SL_ElEta_AfterNjNbCategory[cate]   ,  r-> lepton_eta , wgt_total );
	FillHistogram( h_SL_MinElJetDr_AfterNjNbCategory[cate]  , r-> lepjet_minDr , wgt_total );

	FillHistogram( h_SL_El_LeadingJetPT_AfterNjNbCategory[cate]       , ( r -> jet_pt -> at( 0                 )  ), wgt_total ); 
	FillHistogram( h_SL_El_LeastJetPT_AfterNjNbCategory[cate]         , ( r -> jet_pt -> at( njet_category - 1 )  ), wgt_total );
	FillHistogram( h_SL_El_LeadingJetBtag_AfterNjNbCategory[cate]     , ( r -> jet_btag -> at(0                ) ), wgt_total ); 
	FillHistogram( h_SL_El_LeastJetBtag_AfterNjNbCategory[cate]       , ( r -> jet_btag -> at(njet_category - 1) ), wgt_total ); 
	FillHistogram( h_SL_El_LargestBtagValue_AfterNjNbCategory[cate]   ,  0 , wgt_total ); 
	FillHistogram( h_SL_El_LeastBtagValue_AfterNjNbCategory[cate]     ,  0 , wgt_total ); 

	FillHistogram( h_SL_N_fatjet_NbNjCate_el                    [cate], r->fatjet_pt->size()  , wgt_total ) ;  
	for(unsigned int iFat = 0 ;  iFat < r->fatjet_pt->size() ; iFat ++ ){
	  FillHistogram( h_SL_fatjet_pt_NbNjCate_el		    [cate], r -> fatjet_pt     ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_eta_NbNjCate_el		    [cate], r -> fatjet_eta    ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_phi_NbNjCate_el		    [cate], r -> fatjet_phi    ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_sdmass_NbNjCate_el	            [cate], r -> fatjet_sdmass ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_tau21_NbNjCate_el		    [cate], r -> fatjet_tau21  ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_tau32_NbNjCate_el		    [cate], r -> fatjet_tau32  ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_subjet_loosetagged_NbNjCate_el   [cate], r -> fatjet_subjet_nloosebtag  ->at(iFat), wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_subjet_mediumtagged_NbNjCate_el  [cate], r -> fatjet_subjet_nmediumbtag ->at(iFat), wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_subjet_mean_btagger_NbNjCate_el  [cate], r -> fatjet_subjet_meanbtagger ->at(iFat), wgt_total ) ;  
	}


      }else{
	FillHistogram( h_SL_MuPt_AfterNjNbCategory[cate]    ,  r-> lepton_pt  , wgt_total );
	FillHistogram( h_SL_MuEta_AfterNjNbCategory[cate]   ,  r-> lepton_eta , wgt_total );
	FillHistogram( h_SL_MinMuJetDr_AfterNjNbCategory[cate]  , r-> lepjet_minDr , wgt_total );

	FillHistogram( h_SL_Mu_LeadingJetPT_AfterNjNbCategory[cate]       , (  r -> jet_pt -> at( 0                 )), wgt_total ); 
	FillHistogram( h_SL_Mu_LeastJetPT_AfterNjNbCategory[cate]         , (  r -> jet_pt -> at( njet_category - 1 )), wgt_total );
	FillHistogram( h_SL_Mu_LeadingJetBtag_AfterNjNbCategory[cate]     , (  r -> jet_btag -> at( 0                ) ), wgt_total ); 
	FillHistogram( h_SL_Mu_LeastJetBtag_AfterNjNbCategory[cate]       , (  r -> jet_btag -> at( njet_category - 1) ), wgt_total ); 
	FillHistogram( h_SL_Mu_LargestBtagValue_AfterNjNbCategory[cate]   ,  0 , wgt_total ); 
	FillHistogram( h_SL_Mu_LeastBtagValue_AfterNjNbCategory[cate]     ,  0 , wgt_total ); 

	FillHistogram( h_SL_N_fatjet_NbNjCate_mu                    [cate], r->fatjet_pt->size()  , wgt_total ) ;  
	for(unsigned int iFat = 0 ;  iFat < r->fatjet_pt->size() ; iFat ++ ){
	  FillHistogram( h_SL_fatjet_pt_NbNjCate_mu		    [cate], r -> fatjet_pt     ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_eta_NbNjCate_mu		    [cate], r -> fatjet_eta    ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_phi_NbNjCate_mu		    [cate], r -> fatjet_phi    ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_sdmass_NbNjCate_mu	            [cate], r -> fatjet_sdmass ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_tau21_NbNjCate_mu		    [cate], r -> fatjet_tau21  ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_tau32_NbNjCate_mu		    [cate], r -> fatjet_tau32  ->at(iFat)     , wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_subjet_loosetagged_NbNjCate_mu   [cate], r -> fatjet_subjet_nloosebtag  ->at(iFat), wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_subjet_mediumtagged_NbNjCate_mu  [cate], r -> fatjet_subjet_nmediumbtag ->at(iFat), wgt_total ) ;  
	  FillHistogram( h_SL_fatjet_subjet_mean_btagger_NbNjCate_mu  [cate], r -> fatjet_subjet_meanbtagger ->at(iFat), wgt_total ) ;  
	}

      }
    }


  }

  FillHistogram( h_SL_N_Vtx_AfterAllSelection  , r -> N_Vtx , wgt_total );
  FillHistogram( h_SL_N_VtxWITHOUTPU_AfterAllSelection  , r -> N_Vtx , wgt_total_noPU );
  {
    TRandom3 rand(0); // 0 : If seed is 0, the seed is automatically computed via a TUUID object. In this case the seed is guaranteed to be unique in space and time.
    long reduced_NPV = 0 ; 
    for( int i = 0 ; i < r -> N_Vtx ; i ++ ){
      if( rand.Rndm() < 0.85 ) reduced_NPV++ ;
    }
    if( isMC ){
      FillHistogram( h_SL_N_ManuallyReducedVtx_AfterAllSelection  ,  reduced_NPV , wgt_total );
    }else{
      FillHistogram( h_SL_N_ManuallyReducedVtx_AfterAllSelection  ,  r -> N_Vtx , wgt_total );
    }
  }

  
  FillHistogram( h_SL_isMuonAfterAllSelection  , r -> lepton_isMuon , wgt_total );


  FillHistogram ( h_SL_MetAbs_AfterAllSelection, r -> met_abs , wgt_total ) ;
  FillHistogram ( h_SL_MetPhi_AfterAllSelection, r -> met_phi , wgt_total ) ;

  FillHistogram ( h_SL_MetType1xyAbs_AfterAllSelection, r -> met_type1xy_abs , wgt_total ) ;
  FillHistogram ( h_SL_MetType1xyPhi_AfterAllSelection, r -> met_type1xy_phi , wgt_total ) ;

  if(  r -> met_abs < 30 ){
    FillHistogram ( h_SL_MetPhi_abs0to30GeV_AfterAllSelection, r -> met_phi , wgt_total ) ;
  }else if(  r -> met_abs < 100 ){
    FillHistogram ( h_SL_MetPhi_abs30to100GeV_AfterAllSelection, r -> met_phi , wgt_total ) ;
  }else if(  r -> met_abs < 200 ){
    FillHistogram ( h_SL_MetPhi_abs100GeVto200_AfterAllSelection, r -> met_phi , wgt_total ) ;
  } else {
    FillHistogram ( h_SL_MetPhi_abs200GeVtoinf_AfterAllSelection, r -> met_phi , wgt_total ) ;
  }
  
  FillHistogram ( h_SL_MetEx_AfterAllSelection, r -> met_abs * cos (r -> met_phi) , wgt_total ) ;
  FillHistogram ( h_SL_MetEy_AfterAllSelection, r -> met_abs * sin (r -> met_phi) , wgt_total ) ;

  if( r -> lepton_isMuon == 1 ){
    FillHistogram ( h_SL_MuonPtAfterAllSelection, r -> lepton_pt , wgt_total ) ;
    FillHistogram ( h_SL_MuonEtaAfterAllSelection, r -> lepton_eta , wgt_total ) ;
    FillHistogram ( h_SL_MuonPhiAfterAllSelection, r -> lepton_phi , wgt_total ) ;
    FillHistogram ( h_SL_MuonPt_NoTopPTWgt_AfterAllSelection, r -> lepton_pt , wgt_total / r->wgt_topPT ) ;
    FillHistogram ( h_SL_MuonEta_NoTopPTWgt_AfterAllSelection, r -> lepton_eta , wgt_total/ r->wgt_topPT ) ;
    FillHistogram ( h_SL_MuonPt_NoMuTrigWgt_AfterAllSelection, r -> lepton_pt , wgt_total / r->wgt_trigMu ) ;
    FillHistogram ( h_SL_MuonEta_NoMuTrigWgt_AfterAllSelection, r -> lepton_eta , wgt_total/ r->wgt_trigMu ) ;

    FillHistogram ( h_SL_MetAbs_AfterAllSelection_mu, r -> met_abs , wgt_total ) ;
    FillHistogram ( h_SL_MetPhi_AfterAllSelection_mu, r -> met_phi , wgt_total ) ;

    FillHistogram ( h_SL_MetType1xyAbs_AfterAllSelection_mu, r -> met_type1xy_abs, wgt_total ) ;
    FillHistogram ( h_SL_MetType1xyPhi_AfterAllSelection_mu, r -> met_type1xy_phi, wgt_total ) ;

    FillHistogram ( h_SL_MetEx_AfterAllSelection_mu, r -> met_abs * cos (r -> met_phi) , wgt_total ) ;
    FillHistogram ( h_SL_MetEy_AfterAllSelection_mu, r -> met_abs * sin (r -> met_phi) , wgt_total ) ;

    if( r -> nBJet >=3 ){
      FillHistogram ( h_SL_MuonPtAfterAllSelection_ge3btag, r -> lepton_pt , wgt_total ) ;
    }


    if(  r -> met_abs > 50  ){
      // For the moment, require MET cut to remove QCD.
      
      FillHistogram ( h_SL_N_fatjet_AfterAllSelection_mu                  , r->fatjet_pt->size()  , wgt_total ); ; 
      for(unsigned int i = 0 ;  i < r->fatjet_pt->size() ; i ++ ){                                                  
	FillHistogram ( h_SL_fatjet_pt_AfterAllSelection_mu		  , r -> fatjet_pt     ->at(i)  , wgt_total );  ; 
	FillHistogram ( h_SL_fatjet_eta_AfterAllSelection_mu	          , r -> fatjet_eta    ->at(i)  , wgt_total );     ; 
	FillHistogram ( h_SL_fatjet_phi_AfterAllSelection_mu	          , r -> fatjet_phi    ->at(i)  , wgt_total );     ; 
	FillHistogram ( h_SL_fatjet_sdmass_AfterAllSelection_mu		  , r -> fatjet_sdmass ->at(i)  , wgt_total );     ; 
	FillHistogram ( h_SL_fatjet_tau21_AfterAllSelection_mu		  , r -> fatjet_tau21  ->at(i)  , wgt_total );     ; 
	FillHistogram ( h_SL_fatjet_tau32_AfterAllSelection_mu		  , r -> fatjet_tau32  ->at(i)  , wgt_total );     ; 
	FillHistogram ( h_SL_fatjet_subjet_loosetagged_AfterAllSelection_mu , r -> fatjet_subjet_nloosebtag  ->at(i)  , wgt_total ); ; 
	FillHistogram ( h_SL_fatjet_subjet_mediumtagged_AfterAllSelection_mu, r -> fatjet_subjet_nmediumbtag ->at(i)  , wgt_total ); ; 
	FillHistogram ( h_SL_fatjet_subjet_mean_btagger_AfterAllSelection_mu, r -> fatjet_subjet_meanbtagger ->at(i)  , wgt_total ); ; 
      } // fat jet loop ends.
    }//  Met-cut IF end.
    
  }else{
    FillHistogram ( h_SL_ElPtAfterAllSelection, r -> lepton_pt , wgt_total ) ;
    FillHistogram ( h_SL_ElEtaAfterAllSelection, r -> lepton_eta , wgt_total ) ;
    FillHistogram ( h_SL_ElPhiAfterAllSelection, r -> lepton_phi , wgt_total ) ;
    FillHistogram ( h_SL_ElPt_NoTopPTWgt_AfterAllSelection, r -> lepton_pt , wgt_total / r->wgt_topPT ) ;
    FillHistogram ( h_SL_ElEta_NoTopPTWgt_AfterAllSelection, r -> lepton_eta , wgt_total/ r->wgt_topPT ) ;
    FillHistogram ( h_SL_ElPt_NoElTrigWgt_AfterAllSelection, r -> lepton_pt , wgt_total / r->wgt_trigEl ) ;
    FillHistogram ( h_SL_ElEta_NoElTrigWgt_AfterAllSelection, r -> lepton_eta , wgt_total/ r->wgt_trigEl ) ;

    FillHistogram ( h_SL_ElPtAfterAllSelection_noElSF, r -> lepton_pt , wgt_total / wgt_el_related) ;
    FillHistogram ( h_SL_ElEtaAfterAllSelection_noElSF, r -> lepton_eta , wgt_total / wgt_el_related  ) ;

    FillHistogram ( h_SL_MetAbs_AfterAllSelection_el, r -> met_abs , wgt_total ) ;
    FillHistogram ( h_SL_MetPhi_AfterAllSelection_el, r -> met_phi , wgt_total ) ;

    FillHistogram ( h_SL_MetType1xyAbs_AfterAllSelection_el, r -> met_type1xy_abs , wgt_total ) ;
    FillHistogram ( h_SL_MetType1xyPhi_AfterAllSelection_el, r -> met_type1xy_phi , wgt_total ) ;

    FillHistogram ( h_SL_MetEx_AfterAllSelection_el, r -> met_abs * cos (r -> met_phi) , wgt_total ) ;
    FillHistogram ( h_SL_MetEy_AfterAllSelection_el, r -> met_abs * sin (r -> met_phi) , wgt_total ) ;

    if( r -> nBJet >=3 ){
      FillHistogram ( h_SL_ElPtAfterAllSelection_ge3btag, r -> lepton_pt , wgt_total ) ;
    }

    if(  r -> met_abs > 50  ){
      // For the moment, require MET cut to remove QCD.
      
      FillHistogram ( h_SL_N_fatjet_AfterAllSelection_el                  , r->fatjet_pt->size()  , wgt_total ); ; 
      for(unsigned int i = 0 ;  i < r->fatjet_pt->size() ; i ++ ){
	FillHistogram ( h_SL_fatjet_pt_AfterAllSelection_el		  , r -> fatjet_pt     ->at(i)  , wgt_total );  ; 
	FillHistogram ( h_SL_fatjet_eta_AfterAllSelection_el	          , r -> fatjet_eta    ->at(i)  , wgt_total );     ; 
	FillHistogram ( h_SL_fatjet_phi_AfterAllSelection_el	          , r -> fatjet_phi    ->at(i)  , wgt_total );     ; 
	FillHistogram ( h_SL_fatjet_sdmass_AfterAllSelection_el		  , r -> fatjet_sdmass ->at(i)  , wgt_total );     ; 
	FillHistogram ( h_SL_fatjet_tau21_AfterAllSelection_el		  , r -> fatjet_tau21  ->at(i)  , wgt_total );     ; 
	FillHistogram ( h_SL_fatjet_tau32_AfterAllSelection_el		  , r -> fatjet_tau32  ->at(i)  , wgt_total );     ; 
	FillHistogram ( h_SL_fatjet_subjet_loosetagged_AfterAllSelection_el , r -> fatjet_subjet_nloosebtag  ->at(i)  , wgt_total ); ; 
	FillHistogram ( h_SL_fatjet_subjet_mediumtagged_AfterAllSelection_el, r -> fatjet_subjet_nmediumbtag ->at(i)  , wgt_total ); ; 
	FillHistogram ( h_SL_fatjet_subjet_mean_btagger_AfterAllSelection_el, r -> fatjet_subjet_meanbtagger ->at(i)  , wgt_total ); ; 
      } // fat jet loop ends.
    }//  Met-cut IF end.

  }


  FillHistogram ( h_SL_nJet_AfterAllSelection, r -> nJet , wgt_total ) ;
  FillHistogram ( h_SL_nJet_NoTopPTWgt_AfterAllSelection, r -> nJet , wgt_total / r->wgt_topPT ) ;
  FillHistogram ( h_SL_nJetWITHOUTPUWGT_AfterAllSelection, r -> nJet , wgt_total / r->wgt_PU ) ;

  FillHistogram ( h_SL_nBtagJet_AfterAllSelection, r -> nBJet , wgt_total ) ;
  FillHistogram ( h_SL_nBtagJetWITHOUTBTAGSF_AfterAllSelection, r -> nBJet , wgt_total / r->wgt_btag ) ;

  if( r -> lepton_isMuon == 1 ){
    FillHistogram ( h_SL_nJet_AfterAllSelection_mu , r -> nJet , wgt_total ) ;
    FillHistogram ( h_SL_nBtagJet_AfterAllSelection_mu, r -> nBJet , wgt_total ) ;
  }else{
    FillHistogram ( h_SL_nJet_AfterAllSelection_el , r -> nJet , wgt_total ) ;
    FillHistogram ( h_SL_nBtagJet_AfterAllSelection_el, r -> nBJet , wgt_total ) ;
  }

  {
    for( int i = 0 ; i < 8 && i < r->nJet ; i ++ ){
      FillHistogram ( h_SL_jet_pt_AfterAllSelection  [i], ( r -> jet_pt->at(i))   , wgt_total ); 
      FillHistogram ( h_SL_jet_eta_AfterAllSelection [i], ( r -> jet_eta->at(i))  , wgt_total ); 
      FillHistogram ( h_SL_jetBtag_AfterAllSelection [i], ( r -> jet_btag->at(i)) ,wgt_total );

      if( r -> nBJet >= 3 ){
	FillHistogram ( h_SL_jet_pt_AfterAllSelection_ge3btag  [i], ( r -> jet_pt->at(i))   , wgt_total ); 
      }

      if(  r -> lepton_isMuon == 1  ){
	FillHistogram ( h_SL_jet_pt_AfterAllSelection_mu  [i], ( r -> jet_pt->at(i))   , wgt_total ); 
	FillHistogram ( h_SL_jet_eta_AfterAllSelection_mu [i], ( r -> jet_eta->at(i))  , wgt_total ); 
	FillHistogram ( h_SL_jetBtag_AfterAllSelection_mu [i], ( r -> jet_btag->at(i)) ,wgt_total );
	FillHistogram ( h_SL_jetBtagWOBtagSF_AfterAllSelection_mu [i], ( r -> jet_btag->at(i)) ,wgt_total  / r->wgt_btag );
      }else{
	FillHistogram ( h_SL_jet_pt_AfterAllSelection_el  [i], ( r -> jet_pt->at(i))   , wgt_total ); 
	FillHistogram ( h_SL_jet_eta_AfterAllSelection_el [i], ( r -> jet_eta->at(i))  , wgt_total ); 
	FillHistogram ( h_SL_jetBtag_AfterAllSelection_el [i], ( r -> jet_btag->at(i)) ,wgt_total );
      FillHistogram ( h_SL_jetBtagWOBtagSF_AfterAllSelection_el [i], ( r -> jet_btag->at(i)) ,wgt_total  / r->wgt_btag );
      }


      FillHistogram ( h_SL_jet_pt_WOtopPTwgt_AfterAllSelection  [i], ( r -> jet_pt->at(i))   , wgt_total / r->wgt_topPT ); 
      FillHistogram ( h_SL_jet_eta_WOtopPTwgt_AfterAllSelection [i], ( r -> jet_eta->at(i))  , wgt_total / r->wgt_topPT ); 
      FillHistogram ( h_SL_jetBtagWOBtagSF_AfterAllSelection [i], ( r -> jet_btag->at(i)) ,wgt_total  / r->wgt_btag );
    }

  }
  
  





}


void analyzer::init(){


  ttHSF . init_all();

  f_out = TFile::Open( output.c_str() ,"recreate");
  main_directory = gDirectory;



  tmva_reader = new TMVA::Reader( "!Color:!Silent" );

  tmva_reader ->AddVariable( "ak8jet_pt",                                        & tmva_ak8jet_pt ) ;	 
  tmva_reader ->AddVariable( "ak8jet_eta",					 & tmva_ak8jet_eta  ) ; 	 
  tmva_reader ->AddVariable( "ak8jet_sdmass",					 & tmva_ak8jet_sdmass ) ;	 
  tmva_reader ->AddVariable( "ak8jet_deepcsv_min < 0 ? 0 : ak8jet_deepcsv_min ", & tmva_ak8jet_deepcsv_min ) ;
  tmva_reader ->AddVariable( "ak8jet_doubleb", 					 & tmva_ak8jet_doubleb ) ;    

  tmva_reader -> BookMVA("SatoshiAk8BDT" , "satoshi_ak8_bdt/TMVAClassification_BDT.weights.xml" ) ;


  histogram_list.push_back( h_totalevent = new TH1D("h_totalevent", "h_totalevent", 1, -10 , 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("total event");

  histogram_list.push_back( h_totalevent_noWeight = new TH1D("h_totalevent_NoWeight", "h_totalevent_NoWeight", 1, -10 , 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("total event (mc w/o weight)");


  histogram_list.push_back( h_pileup_noWeight = new TH1D("h_pileup_noWeight",
							 "h_pileup_noWeight", 
							 102, -2 , 100 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Pileup. noWgt");

  histogram_list.push_back( h_pileup_PUWeight = new TH1D("h_pileup_PUWeight",
							 "h_pileup_PUWeight", 
							 102, -2 , 100 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Pileup. with PUWgt");



  histogram_list.push_back( h_wgt_TOTAL  = new TH1D("wgt_TOTAL", "wgt_TOTAL", 100, -10 , 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight total");

  histogram_list.push_back( h_wgt_PU     = new TH1D("wgt_PU", "wgt_PU", 120, -1 , 5 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight PU");

  histogram_list.push_back( h_wgt_btag   = new TH1D("wgt_btag", "wgt_btag", 200, -1 , 3 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight btag");

  histogram_list.push_back( h_wgt_trigMu = new TH1D("wgt_trigMu", "wgt_trig Mu", 100, 0 , 2 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight trigMu");

  histogram_list.push_back( h_wgt_trigEl = new TH1D("wgt_trigEl", "wgt_trigEl", 60, 0.5 , 1.1 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight trig El");

  histogram_list.push_back( h_wgt_muID   = new TH1D("wgt_muID", "wgt_muID", 100, 0 , 2 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight Mu ID");

  histogram_list.push_back( h_wgt_elID   = new TH1D("wgt_elID", "wgt_elID", 100, 0 , 2 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight El ID");

  histogram_list.push_back( h_wgt_elReco = new TH1D("wgt_elReco", "wgt_elReco", 100, 0 , 2 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight El Reco");

  histogram_list.push_back( h_wgt_muIso  = new TH1D("wgt_muIso", "wgt_muIso", 100, 0 , 2 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight Mu Iso");

  histogram_list.push_back( h_wgt_MCEven = new TH1D("wgt_MCEven", "wgt_MCEven", 20, -2. , 2. ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight MCEvent");

  histogram_list.push_back( h_wgt_topPT  = new TH1D("wgt_topPT", "wgt_topPT", 10, -10 , 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight topPT");

  histogram_list.push_back( h_wgt_EGZvtx  = new TH1D("wgt_EGZvtx", "wgt_EGZvtx", 200, 0.9 , 1.1 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Weight EG Zvtx");


#ifdef TRUTH_INFO_STUDY


  outtree_odd = new TTree( "treeOdd", "treeOdd");
  outtree_odd -> SetAutoFlush( - 3 * 100 * 1000 ) ; // 0.3 MB auto flush instead of 30MB. (hearder threshold because of more(~70) systematics.)

  outtree_odd -> Branch( "EventNumber" , & eventnumber          , "EventNumber/I" );
  outtree_odd -> Branch( "eventcategory" , & eventcategory          , "eventcategory/I" );
  outtree_odd -> Branch( "EventWeight" , & eventweight , "EventWeight/F" );
  outtree_odd -> Branch( "ak8jet_genhiggsmatch" , & ak8jet_genhiggsmatch , "ak8jet_genhiggsmatch/I" );
  outtree_odd -> Branch( "ak8jet_pt" , & ak8jet_pt , "ak8jet_pt/F" );
  outtree_odd -> Branch( "ak8jet_eta" , & ak8jet_eta , "ak8jet_eta/F" );
  outtree_odd -> Branch( "ak8jet_sdmass" , & ak8jet_sdmass , "ak8jet_sdmass/F" );
  outtree_odd -> Branch( "ak8jet_deepcsv_min" , & ak8jet_deepcsv_min , "ak8jet_deepcsv_min/F" );
  outtree_odd -> Branch( "ak8jet_doubleb" , & ak8jet_doubleb  , "ak8jet_doubleb/F" );
  outtree_odd -> Branch( "ak8jet_tau21" , & ak8jet_tau21  , "ak8jet_tau21/F" );
  outtree_odd -> Branch( "ak8jet_tau32" , & ak8jet_tau32  , "ak8jet_tau32/F" );


  outtree_even = new TTree( "treeEven", "treeEven");
  outtree_even -> SetAutoFlush( - 3 * 100 * 1000 ) ; // 0.3 MB auto flush instead of 30MB. (hearder threshold because of more(~70) systematics.)

  outtree_even -> Branch( "EventNumber" , & eventnumber          , "EventNumber/I" );
  outtree_even -> Branch( "eventcategory" , & eventcategory          , "eventcategory/I" );
  outtree_even -> Branch( "EventWeight" , & eventweight , "EventWeight/F" );
  outtree_even -> Branch( "ak8jet_genhiggsmatch" , & ak8jet_genhiggsmatch , "ak8jet_genhiggsmatch/I" );
  outtree_even -> Branch( "ak8jet_pt" , & ak8jet_pt , "ak8jet_pt/F" );
  outtree_even -> Branch( "ak8jet_eta" , & ak8jet_eta , "ak8jet_eta/F" );
  outtree_even -> Branch( "ak8jet_sdmass" , & ak8jet_sdmass , "ak8jet_sdmass/F" );
  outtree_even -> Branch( "ak8jet_deepcsv_min" , & ak8jet_deepcsv_min , "ak8jet_deepcsv_min/F" );
  outtree_even -> Branch( "ak8jet_doubleb" , & ak8jet_doubleb  , "ak8jet_doubleb/F" );
  outtree_even -> Branch( "ak8jet_tau21" , & ak8jet_tau21  , "ak8jet_tau21/F" );
  outtree_even -> Branch( "ak8jet_tau32" , & ak8jet_tau32  , "ak8jet_tau32/F" );

  for( unsigned int i = 0 ; i < 2 ; i ++ ){

    char name[100];

    sprintf( name , "MinB_NoCut_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_MinBTagDiscri_NoCut  [i]    = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "Min(DeepCSV), NoCut, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "DoubleB_NoCut_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_DoubleB_NoCut  [i]  = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "DoubleB, NoCut, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name , "MinB_6j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_MinBTagDiscri_6j4b  [i]    = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "Min(DeepCSV), 6j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "DoubleB_6j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_DoubleB_6j4b  [i]  = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "DoubleB, 6j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name , "AK8BDT_6j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_ak8BDT_6j4b  [i]  = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "AK8 BDT, 6j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);



    sprintf( name , "HandMergeBtag_6j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_HandMergedBtag_6j4b [i]  = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "DoubleB, 6j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);




    sprintf( name , "MinB_6j3b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_MinBTagDiscri_6j3b  [i]    = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "Min(DeepCSV), 6j3b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "DoubleB_6j3b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_DoubleB_6j3b  [i]  = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "DoubleB, 6j3b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    sprintf( name , "MinB_6j2b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_MinBTagDiscri_6j2b  [i]    = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "Min(DeepCSV), 6j2b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "DoubleB_6j2b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_DoubleB_6j2b  [i]  = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "DoubleB, 6j2b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name , "MinB_5j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_MinBTagDiscri_5j4b  [i]    = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "Min(DeepCSV), 5j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "DoubleB_5j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_DoubleB_5j4b  [i]  = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "DoubleB, 5j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name , "MinB_4j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_MinBTagDiscri_4j4b  [i]    = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "Min(DeepCSV), 4j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "DoubleB_4j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_DoubleB_4j4b  [i]  = new TH1D (name,name, 1000, -2, 2 ) ) ;
    sprintf( name , "DoubleB, 4j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);





    sprintf( name , "SL_MinB_NoCut_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_MinBTagDiscri_NoCut  [i]    = new TH1D (name,name, 20, 0, 1 ) ) ;
    sprintf( name , "Min(DeepCSV), NoCut, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "SL_DoubleB_NoCut_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_DoubleB_NoCut  [i]  = new TH1D (name,name, 20, -1, 1 ) ) ;
    sprintf( name , "DoubleB, NoCut, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name , "SL_MinB_6j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_MinBTagDiscri_6j4b  [i]    = new TH1D (name,name, 20, 0, 1 ) ) ;
    sprintf( name , "Min(DeepCSV), 6j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "SL_DoubleB_6j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_DoubleB_6j4b  [i]  = new TH1D (name,name, 20, -1, 1 ) ) ;
    sprintf( name , "DoubleB, 6j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name , "SL_MinB_6j3b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_MinBTagDiscri_6j3b  [i]    = new TH1D (name,name, 20, 0, 1 ) ) ;
    sprintf( name , "Min(DeepCSV), 6j3b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "SL_DoubleB_6j3b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_DoubleB_6j3b  [i]  = new TH1D (name,name, 20, -1, 1 ) ) ;
    sprintf( name , "DoubleB, 6j3b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    sprintf( name , "SL_MinB_6j2b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_MinBTagDiscri_6j2b  [i]    = new TH1D (name,name, 20, 0, 1 ) ) ;
    sprintf( name , "Min(DeepCSV), 6j2b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "SL_DoubleB_6j2b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_DoubleB_6j2b  [i]  = new TH1D (name,name, 20, -1, 1 ) ) ;
    sprintf( name , "DoubleB, 6j2b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name , "SL_MinB_5j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_MinBTagDiscri_5j4b  [i]    = new TH1D (name,name, 20, 0, 1 ) ) ;
    sprintf( name , "Min(DeepCSV), 5j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "SL_DoubleB_5j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_DoubleB_5j4b  [i]  = new TH1D (name,name, 20, -1, 1 ) ) ;
    sprintf( name , "DoubleB, 5j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name , "SL_MinB_4j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_MinBTagDiscri_4j4b  [i]    = new TH1D (name,name, 20, 0, 1 ) ) ;
    sprintf( name , "Min(DeepCSV), 4j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
    
    sprintf( name , "SL_DoubleB_4j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_DoubleB_4j4b  [i]  = new TH1D (name,name, 20, -1, 1 ) ) ;
    sprintf( name , "DoubleB, 4j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name , "SL_HandMergedBtag_6j4b_Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list . push_back(  h_SL_HandMergedBtag_6j4b [i]  = new TH1D (name,name, 23, -1, 1.3 ) ) ;
    sprintf( name , "MergedBtag, 6j4b, Higgs%s" , i == 0 ? "Match" : "NotMatch" ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);











#else 


  if( b_FakeEstimation ){

    TFile * tf_FakeSF = TFile::Open("FakeScaleFactor_nominal.root");
    tf_FakeSF -> GetObject( "Mu_FakeSF" , h_FakeSF_Mu );
    tf_FakeSF -> GetObject( "El_FakeSF" , h_FakeSF_El );

    if( h_FakeSF_Mu == 0 ){
      std::cout <<"[critical] analyzer.cc : Tried to get the Muon Fake SF from the SF file, but failed."<< std::endl ;
    }
    if( h_FakeSF_El == 0 ){
      std::cout <<"[critical] analyzer.cc : Tried to get the Electron Fake SF from the SF file, but failed."<< std::endl ;
    }

    assert(  h_FakeSF_Mu != 0 && h_FakeSF_El != 0 ) ;

  }

  //  outtree   = new TTree("tree","tree");
  //  outtree -> SetAutoFlush( - 3 * 1000 * 1000 ) ; // 3 MB auto flush instead of 30MB.
  //  outtree -> SetAutoFlush( - 3 * 1000 * 100 ) ; // 0.3 MB auto flush instead of 30MB.
  


  histogram_list.push_back( h_SL_cutflow[0] = new TH1D("h_cutflow_el", "h_cutflow_el", 8, -0.5 , 7.5 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Cutflow, el");
  histogram_list.push_back( h_SL_cutflow[1] = new TH1D("h_cutflow_mu", "h_cutflow_mu", 8, -0.5 , 7.5 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Cutflow, mu");


  histogram_list.push_back( h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons = new TH1D("SL_Mee_AfterSingleTriggerAndExactTwoTightleptons","SL_Mee_AfterSingleTriggerAndExactTwoTightleptons", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mee [GeV]");

  histogram_list.push_back( h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_noElSF = new TH1D("SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_noElSF","SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_noElSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mee [GeV], w/o EL SSs");


  histogram_list.push_back( h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElRecoSF = new TH1D("SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElRecoSF","SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElRecoSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mee [GeV], w/o EL SFs, with El Reco SF");

  histogram_list.push_back( h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElTrigSF = new TH1D("SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElTrigSF","SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElTrigSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mee [GeV], w/o EL SFs, with El Trig SF");

  histogram_list.push_back( h_SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElIDSF = new TH1D("SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElIDSF","SL_Mee_AfterSingleTriggerAndExactTwoTightleptons_onlyElIDSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mee [GeV], w/o EL SFs, with El  ID SF");



  histogram_list.push_back( h_SL_Mmm_AfterSingleTriggerAndExactTwoTightleptons = new TH1D("SL_Mmm_AfterSingleTriggerAndExactTwoTightleptons","SL_Mmm_AfterSingleTriggerAndExactTwoTightleptons", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("M mumu [GeV]");

  histogram_list.push_back( h_SL_Mmm_AfterSingleTriggerAndExactTwoTightleptons_noMuSF = new TH1D("SL_Mmm_AfterSingleTriggerAndExactTwoTightleptons_noMuSF","SL_Mmm_AfterSingleTriggerAndExactTwoTightleptons_noMuSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("M mumu [GeV], no Muon SFs");

  histogram_list.push_back( h_SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut 
			    = new TH1D("SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut",
				       "SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Electron Pt[GeV]");

  histogram_list.push_back( h_SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut 
			    = new TH1D("SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut",
				       "SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Leading Electron Pt[GeV]");

  histogram_list.push_back( h_SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut_noElSF
			    = new TH1D("SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut_noElSF",
				       "SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut_noElSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Electron Pt[GeV], no El SFs");

  histogram_list.push_back( h_SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut_noElSF 
			    = new TH1D("SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut_noElSF",
				       "SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut_noElSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Leading Electron Pt[GeV], no El SFs");

  histogram_list.push_back( h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut 
			    = new TH1D("SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut",
				       "SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Muon Pt[GeV]");

  histogram_list.push_back( h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut 
			    = new TH1D("SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut",
				       "SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Leading Muon Pt[GeV]");

  histogram_list.push_back( h_SL_nJet_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut
			    = new TH1D("SL_nJet_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut",
				       "SL_nJet_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut", 10, 0, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("nJet");  

  histogram_list.push_back( h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut
			    = new TH1D("SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut",
				       "SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCut", 80, 0, 400 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Jet Pt[GeV]");    

  histogram_list.push_back( h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut
			    = new TH1D("SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut",
				       "SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut", 80, 0, 400 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Jet Pt[GeV]");      

  histogram_list.push_back( h_SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut
			    = new TH1D("SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut",
				       "SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut", 80, 0, 400 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Leading Jet Pt[GeV]");      

  histogram_list.push_back( h_SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut
			    = new TH1D("SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut",
				       "SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut", 80, 0, 400 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Third Leading Jet Pt[GeV]");    
  
  histogram_list.push_back( h_SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut
			    = new TH1D("SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut",
				       "SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMeeCutAND4JetCut", 80, 0, 400 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Fourth Jet Pt[GeV]");    


  histogram_list.push_back( h_SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut 
			    = new TH1D("SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_El1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Electron PT[GeV]");    
  histogram_list.push_back( h_SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut 
			    = new TH1D("SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_El2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Electron PT[GeV]");    

  histogram_list.push_back( h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut 
			    = new TH1D("SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Muon PT[GeV]");    
  histogram_list.push_back( h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut 
			    = new TH1D("SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Leading Muon PT[GeV]");    


  histogram_list.push_back( h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuSF 
			    = new TH1D("SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuSF",
				       "SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Muon PT[GeV], no Mu SFs");    
  histogram_list.push_back( h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuSF 
			    = new TH1D("SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuSF",
				       "SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Leading Muon PT[GeV], no Mu SFs");    

  histogram_list.push_back( h_SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuTrigSF
			    = new TH1D("SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuTrigSF",
				       "SL_Mu1_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuTrigSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Muon PT[GeV], no Mu Trig SF");    
  histogram_list.push_back( h_SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuTrigSF
			    = new TH1D("SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuTrigSF",
				       "SL_Mu2_PT_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_noMuTrigSF", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Leading Muon PT[GeV], no Mu Trig SF");    

  histogram_list.push_back( h_SL_nJet_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut
			    = new TH1D("SL_nJet_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_nJet_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 10, 0, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("nJet");    

  histogram_list.push_back( h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut
			    = new TH1D("SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 80, 0, 400 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Jet PT[GeV]");      


  histogram_list.push_back( h_SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut
			    = new TH1D("SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 80, 0, 300 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("2nd leading Jet PT[GeV]");      


  histogram_list.push_back( h_SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut
			    = new TH1D("SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 80, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("3rd leading Jet PT[GeV]");      


  histogram_list.push_back( h_SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut
			    = new TH1D("SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 40, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("4th leading Jet PT[GeV]");      


  histogram_list.push_back( h_SL_JetPt_5_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut
			    = new TH1D("SL_JetPt_5_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_JetPt_5_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 20, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("5th leading Jet PT[GeV]");      


  histogram_list.push_back( h_SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut 
			    = new TH1D("SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 20, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met abs[GeV]");      



  histogram_list.push_back( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet0
			    = new TH1D("SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet0",
				       "SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet0", 20, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met abs[GeV], 0Jet");      

  histogram_list.push_back( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet1
			    = new TH1D("SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet1",
				       "SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet1", 20, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met abs[GeV], 1Jet");      


  histogram_list.push_back( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet2
			    = new TH1D("SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet2",
				       "SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet2", 20, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met abs[GeV], 2Jet");      


  histogram_list.push_back( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet3
			    = new TH1D("SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet3",
				       "SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet3", 20, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met abs[GeV], 3Jet");      


  histogram_list.push_back( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet4
			    = new TH1D("SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet4",
				       "SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet4", 20, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met abs[GeV], 4Jet");      


  histogram_list.push_back( h_SL_Met_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet5more
			    = new TH1D("SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet5",
				       "SL_MetAbs_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_jet5", 20, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met abs[GeV], 5/moreJet");      






  histogram_list.push_back( h_SL_MetPhi_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut 
			    = new TH1D("SL_MetPhi_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut",
				       "SL_MetPhi_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut", 20, -M_PI , M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi");      



  histogram_list.push_back( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet0
			    = new TH1D("SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet0",
				       "SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet0", 24, 91 - 12, 91+12) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mmm[GeV], 0-jet event");      


  histogram_list.push_back( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet1
			    = new TH1D("SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet1",
				       "SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet1", 24, 91 - 12, 91+12) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mmm[GeV], 1-jet event");      

  histogram_list.push_back( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet2
			    = new TH1D("SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet2",
				       "SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet2", 24, 91 - 12, 91+12) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mmm[GeV], 2-jet event");      


  histogram_list.push_back( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet3
			    = new TH1D("SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet3",
				       "SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet3", 24, 91 - 12, 91+12) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mmm[GeV], 3-jet event");      


  histogram_list.push_back( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet4
			    = new TH1D("SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet4",
				       "SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet4", 24, 91 - 12, 91+12) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mmm[GeV], 4-jet event");      

  histogram_list.push_back( h_SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet5more
			    = new TH1D("SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet5",
				       "SL_Mmumu_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCut_NJet5", 24, 91 - 12, 91+12) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mmm[GeV], 5/more-jet event");      


  histogram_list.push_back( h_SL_InvMass_AfterMuMuSameSignSelection              = new TH1D("SL_InvMass_AfterMuMuSameSignSelection"              ,"SL_InvMass_AfterMuMuSameSignSelection"              , 20 , 41.0 , 141.0 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("InvMass Mmumu[GeV]");      
  histogram_list.push_back( h_SL_MuonPT__AfterMuMuSameSignSelection_ZmassCut 	 = new TH1D("SL_MuonPT__AfterMuMuSameSignSelection_ZmassCut" 	 ,"SL_MuonPT__AfterMuMuSameSignSelection_ZmassCut"     , 16 , 20, 100 ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon PT[GeV]");      
  histogram_list.push_back( h_SL_MuonEta_AfterMuMuSameSignSelection_ZmassCut 	 = new TH1D("SL_MuonEta_AfterMuMuSameSignSelection_ZmassCut" 	 ,"SL_MuonEta_AfterMuMuSameSignSelection_ZmassCut"     , 12 , -2.4, 2.4)); 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon Eta");      

  histogram_list.push_back( h_SL_InvMass_AfterElElSameSignSelection 		 = new TH1D("SL_InvMass_AfterElElSameSignSelection" 		 ,"SL_InvMass_AfterElElSameSignSelection" 	       , 20 , 41.0 , 141.0 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("InvMass Mee[GeV]");      
  histogram_list.push_back( h_SL_ElectronPT__AfterElElSameSignSelection_ZmassCut = new TH1D("SL_ElectronPT__AfterElElSameSignSelection_ZmassCut" ,"SL_ElectronPT__AfterElElSameSignSelection_ZmassCut" , 16 , 20, 100 ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron PT[GeV]");      
  histogram_list.push_back( h_SL_ElectronEta_AfterElElSameSignSelection_ZmassCut = new TH1D("SL_ElectronEta_AfterElElSameSignSelection_ZmassCut" ,"SL_ElectronEta_AfterElElSameSignSelection_ZmassCut" , 12 , -2.4, 2.4)); 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron Eta");      

  histogram_list.push_back( h_SL_Mu1_PT_WWenrichRegion
			    = new TH1D("SL_Mu1_PT_WWenrichRegion",
				       "SL_Mu1_PT_WWenrichRegion", 20, 0, 200) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Muon PT[GeV]");    
  histogram_list.push_back( h_SL_Mu2_PT_WWenrichRegion
			    = new TH1D("SL_Mu2_PT_WWenrichRegion",
				       "SL_Mu2_PT_WWenrichRegion", 20, 0, 200) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Leading Muon PT[GeV]");    


  histogram_list.push_back( h_SL_El1_PT_WWenrichRegion
			    = new TH1D("SL_El1_PT_WWenrichRegion",
				       "SL_El1_PT_WWenrichRegion", 20, 0, 200) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Electron PT[GeV]");    
  histogram_list.push_back( h_SL_El2_PT_WWenrichRegion
			    = new TH1D("SL_El2_PT_WWenrichRegion",
				       "SL_El2_PT_WWenrichRegion", 20, 0, 200) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Leading Electron PT[GeV]");    


  histogram_list.push_back( h_SL_LeptopnPT_WWemuEnrichRegionLoose  = new TH1D("SL_LeptonPT_WWemuEnrichRegionLoose","SL_LeptonPT_WWemuEnrichRegionLoose", 10,0,200) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Lepton PT[GeV]");      

  histogram_list.push_back( h_MET_WWemuEnrichRegionLoose = new TH1D("SL_MET_WWemuEnrichRegionLoose","SL_MET_WWemuEnrichRegionLoose", 10,0,200) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met[GeV]");      

  histogram_list.push_back( h_SL_LeptopnPT_WWemuEnrichRegionTight  = new TH1D("SL_LeptonPT_WWemuEnrichRegionTight","SL_LeptonPT_WWemuEnrichRegionTight", 10,0,200) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Lepton PT[GeV]");      

  histogram_list.push_back( h_MET_WWemuEnrichRegionTight = new TH1D("SL_MET_WWemuEnrichRegionTight","SL_MET_WWemuEnrichRegionTight", 10,0,200) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met[GeV]");      




  histogram_list.push_back( h_SL_MuonPT_WWemuSameSignEnrichRegionLoose  = new TH1D("SL_MuonPT_WWemuSameSignEnrichRegionLoose","SL_MuonPT_WWemuSameSignEnrichRegionLoose", 10,0,200) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon PT[GeV]");      

  histogram_list.push_back( h_SL_ElectronPT_WWemuSameSignEnrichRegionLoose  = new TH1D("SL_ElectronPT_WWemuSameSignEnrichRegionLoose","SL_ElectronPT_WWemuSameSignEnrichRegionLoose", 10,0,200) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron PT[GeV]");      

  histogram_list.push_back( h_MET_WWemuSameSignEnrichRegionLoose = new TH1D("SL_MET_WWemuSameSignEnrichRegionLoose","SL_MET_WWemuSameSignEnrichRegionLoose", 10,0,200) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met[GeV]");      

  histogram_list.push_back( h_nJet_WWemuSameSignEnrichRegionLoose = new TH1D("SL_nJet_WWemuSameSignEnrichRegionLoose","SL_nJet_WWemuSameSignEnrichRegionLoose", 10,0,10) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("nJet");      

  histogram_list.push_back( h_nBJet_WWemuSameSignEnrichRegionLoose = new TH1D("SL_nBJet_WWemuSameSignEnrichRegionLoose","SL_nBJet_WWemuSameSignEnrichRegionLoose", 5,0,5) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("n-Btagged-Jet");      



  histogram_list.push_back( h_SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut
			    = new TH1D("SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut",
				       "SL_JetPt_1_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut", 80, 0, 400 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Jet PT[GeV]");      
  
  histogram_list.push_back( h_SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut
			    = new TH1D("SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut",
				       "SL_JetPt_2_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut", 80, 0, 400 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Second Leading Jet PT[GeV]");      
  
  histogram_list.push_back( h_SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut
			    = new TH1D("SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut",
				       "SL_JetPt_3_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut", 80, 0, 400 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Third Leading Jet PT[GeV]");        

  histogram_list.push_back( h_SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut
			    = new TH1D("SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut",
				       "SL_JetPt_4_AfterSingleTriggerAndExactTwoTightleptonsAndMmmCutAND4JetCut", 80, 0, 400 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Fourth Leading Jet PT[GeV]");        


  
  histogram_list.push_back( h_SL_N_Vtx_AfterSingleLepTrigAndTightLep = new TH1D("SL_N_Vtx_AfterSingleLepTrigAndTightLep", "SL_N_Vtx_AfterSingleLepTrigAndTightLep", 50, 0 , 50 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Number of vertexes");        
  histogram_list.push_back( h_SL_N_VtxWITHOUTPU_AfterSingleLepTrigAndTightLep
			    = new TH1D("SL_N_VtxWITHOUTPU_AfterSingleLepTrigAndTightLep", "SL_N_VtxWITHOUTPU_AfterSingleLepTrigAndTightLep", 50, 0 , 50 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Number of vertexes, no PU reweighting");        

  histogram_list.push_back( h_SL_EventCategorizationBasedonNjetNBtagJet
			    = new TH1D("SL_EventCategorizationBasedonNjetNBtagJet", "SL_EventCategorizationBasedonNjetNBtagJet"
				       , 9 , 1, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category");        

  histogram_list.push_back( h_SL_EventCategorizationBasedonNjetNBtagJet_withMetCut
			    = new TH1D("SL_EventCategorizationBasedonNjetNBtagJet_withMetCut", "SL_EventCategorizationBasedonNjetNBtagJet_withmetCut"
				       , 9 , 1, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category(with met cut)");        



  histogram_list.push_back( h_SL_EventCategorizationBasedonNjetNBtagJet_withMetCut_el
			    = new TH1D("SL_EventCategorizationBasedonNjetNBtagJet_withMetCut_el", "SL_EventCategorizationBasedonNjetNBtagJet_withmetCut_el"
				       , 9 , 1, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category(el, with met cut)");        



  histogram_list.push_back( h_SL_EventCategorizationBasedonNjetNBtagJet_withMetCut_mu
			    = new TH1D("SL_EventCategorizationBasedonNjetNBtagJet_withMetCut_mu", "SL_EventCategorizationBasedonNjetNBtagJet_withmetCut_mu"
				       , 9 , 1, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category(mu, with met cut)");        



  histogram_list.push_back( h_SL_Much_JetPtCategorization4to6
			    = new TH1D("SL_Much_JetPtCategorization4to6","SL_Much_JetPtCategorization4to6",
				       2*pow(4,6) , 0, 2* pow(4,6 ) ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mu-ch JetPtCategorization 4to6");        

  histogram_list.push_back( h_SL_Elch_JetPtCategorization4to6
			    = new TH1D("SL_Elch_JetPtCategorization4to6","SL_Elch_JetPtCategorization4to6",
				       2* pow(4,6) , 0, 2* pow(4,6 ) ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("El-ch JetPtCategorization 4to6");        



  histogram_list.push_back( h_SL_Much_JetPtCategorization3to6
			    = new TH1D("SL_Much_JetPtCategorization3to6","SL_Much_JetPtCategorization3to6",
				       2*pow(3,6) , 0, 2*pow(3,6 ) )  );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mu-ch JetPtCategorization 3to6");        

  histogram_list.push_back( h_SL_Elch_JetPtCategorization3to6
			    = new TH1D("SL_Elch_JetPtCategorization3to6","SL_Elch_JetPtCategorization3to6",
				       2* pow(3,6) , 0, 2*pow(3,6 ) ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("El-ch JetPtCategorization 3to6");        




  histogram_list.push_back( h_SL_Much_JetPtCategorization2to4and3to2
			    = new TH1D("SL_Much_JetPtCategorization2to4and3to2","SL_Much_JetPtCategorization2to4and3to2",
				       2* pow(2,4)*pow(3,2) , 0, 2* pow(2,4)*pow(3,2) ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mu-ch JetPtCategorization 2to4and3to2");        

  histogram_list.push_back( h_SL_Elch_JetPtCategorization2to4and3to2
			    = new TH1D("SL_Elch_JetPtCategorization2to4and3to2","SL_Elch_JetPtCategorization2to4and3to2",
				       2* pow(2,4)*pow(3,2) , 0, 2*pow(2,4)*pow(3,2) ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("El-ch JetPtCategorization 2to4and3to2");        



  histogram_list.push_back( h_SL_eventCate_Mu_NoWeight
                            = new TH1D("SL_eventCate_Mu_NoWeight", "SL_eventCate_Mu_NoWeight"  , 9 , 1, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category");

  histogram_list.push_back( h_SL_eventCate_El_NoWeight
                            = new TH1D("SL_eventCate_El_NoWeight", "SL_eventCate_El_NoWeight"  , 9 , 1, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category");

  histogram_list.push_back( h_SL_eventCate_Mu_NoWeight_noMetCut
                            = new TH1D("SL_eventCate_Mu_NoWeight_noMetCut", "SL_eventCate_Mu_NoWeight_noMetCut"  , 9 , 1, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category");

  histogram_list.push_back( h_SL_eventCate_El_NoWeight_noMetCut
                            = new TH1D("SL_eventCate_El_NoWeight_noMetCut", "SL_eventCate_El_NoWeight_noMetCut"  , 9 , 1, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category");
  



  histogram_list.push_back( h_SL_EventCategorizationBasedonNjetNBtagJet_Mu
			    = new TH1D("SL_EventCategorizationBasedonNjetNBtagJet_Mu", "SL_EventCategorizationBasedonNjetNBtagJet_Mu"
				       , 10 , 0, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category");        


  histogram_list.push_back( h_SL_EventCategorizationBasedonNjetNBtagJet_El
			    = new TH1D("SL_EventCategorizationBasedonNjetNBtagJet_El", "SL_EventCategorizationBasedonNjetNBtagJet_El"
				       , 10 , 0, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category");        


  histogram_list.push_back( h_SL_EventCategorizationBasedonNjetNBtagJet_FakeCR_El
			    = new TH1D("SL_EventCategorizationBasedonNjetNBtagJet_FakeCR_El", "SL_EventCategorizationBasedonNjetNBtagJet_FakeCR_El"
				       , 10 , 0, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category");        


  histogram_list.push_back( h_SL_EventCategorizationBasedonNjetNBtagJet_FakeCR_Mu
			    = new TH1D("SL_EventCategorizationBasedonNjetNBtagJet_FakeCR_Mu", "SL_EventCategorizationBasedonNjetNBtagJet_FakeCR_Mu"
				       , 10 , 0, 10 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NjNb event category");        

  histogram_list.push_back( h_SL_N_Vtx_AfterAllSelection = new TH1D("SL_N_Vtx_AfterAllSelection", "SL_N_Vtx_AfterAllSelection", 80, 0 , 80 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Number of vertexes");        

  histogram_list.push_back( h_SL_N_ManuallyReducedVtx_AfterAllSelection = new TH1D("h_SL_N_ManuallyReducedVtx_AfterAllSelection", "h_SL_N_ManuallyReducedVtx_AfterAllSelection", 80, 0 , 80 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("RandomlyRedycedNVertex");        

  histogram_list.push_back( h_SL_N_VtxWITHOUTPU_AfterAllSelection
			    = new TH1D("SL_N_VtxWITHOUTPU_AfterAllSelection", "SL_N_VtxWITHOUTPU_AfterAllSelection", 80, 0 , 80 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Number of vertexes, no PU reweighting");

  histogram_list.push_back( h_SL_isMuonAfterAllSelection = new TH1D("SL_isMuon_AfterAllSelection","h_SL_isMuonAfterAllSelection",  2, 0, 2 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Lepton flavour : 0=Electron, 1=Muon");

  histogram_list.push_back( h_SL_MetAbs_AfterAllSelection    = new TH1D("h_SL_MetAbs_AfterAllSelection"   ,"h_SL_MetAbs_AfterAllSelection",  40, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Absolute MET [GeV]");
  histogram_list.push_back( h_SL_MetAbs_AfterAllSelection_el = new TH1D("h_SL_MetAbs_AfterAllSelection_el","h_SL_MetAbs_AfterAllSelection_el",  40, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Absolute MET [GeV], el channel");
  histogram_list.push_back( h_SL_MetAbs_AfterAllSelection_mu = new TH1D("h_SL_MetAbs_AfterAllSelection_mu","h_SL_MetAbs_AfterAllSelection_mu",  40, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Absolute MET [GeV], mu channel");


  histogram_list.push_back( h_SL_MetType1xyAbs_AfterAllSelection    = new TH1D("h_SL_MetType1xyAbs_AfterAllSelection"   ,"h_SL_MetType1xyAbs_AfterAllSelection",  40, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Absolute MET (type1xy)[GeV]");
  histogram_list.push_back( h_SL_MetType1xyAbs_AfterAllSelection_el = new TH1D("h_SL_MetType1xyAbs_AfterAllSelection_el","h_SL_MetType1xyAbs_AfterAllSelection_el",  40, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Absolute MET (type1xy)[GeV], el channel");
  histogram_list.push_back( h_SL_MetType1xyAbs_AfterAllSelection_mu = new TH1D("h_SL_MetType1xyAbs_AfterAllSelection_mu","h_SL_MetType1xyAbs_AfterAllSelection_mu",  40, 0, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Absolute MET (type1xy)[GeV], mu channel");
 
  histogram_list.push_back( h_SL_MetPhi_abs0to30GeV_AfterAllSelection    = new TH1D("h_SL_MetPhi_abs0to30GeV_AfterAllSelection"    ,"h_SL_MetPhi_abs0to30GeV_AfterAllSelection"   ,  40, -M_PI ,    M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi (AbsMet0-30GeV)");
  histogram_list.push_back( h_SL_MetPhi_abs30to100GeV_AfterAllSelection  = new TH1D("h_SL_MetPhi_abs30to100GeV_AfterAllSelection"  ,"h_SL_MetPhi_abs30to100GeV_AfterAllSelection" ,  40, -M_PI ,    M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi (AbsMet30-100GeV)");
  histogram_list.push_back( h_SL_MetPhi_abs100GeVto200_AfterAllSelection = new TH1D("h_SL_MetPhi_abs100GeVto200_AfterAllSelection" ,"h_SL_MetPhi_abs100GeVto200_AfterAllSelection",  40, -M_PI ,    M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi (AbsMet100-200GeV)");
  histogram_list.push_back( h_SL_MetPhi_abs200GeVtoinf_AfterAllSelection = new TH1D("h_SL_MetPhi_abs200GeVtoInf_AfterAllSelection" ,"h_SL_MetPhi_abs200GeVtoInf_AfterAllSelection",  40, -M_PI ,    M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi (AbsMet200GeV-)");

  histogram_list.push_back( h_SL_MetPhi_AfterAllSelection    = new TH1D("h_SL_MetPhi_AfterAllSelection"   ,"h_SL_MetPhi_AfterAllSelection",     40, -M_PI ,    M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi");
  histogram_list.push_back( h_SL_MetPhi_AfterAllSelection_el = new TH1D("h_SL_MetPhi_AfterAllSelection_el","h_SL_MetPhi_AfterAllSelection_el",  40, -M_PI , M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi, el channel");
  histogram_list.push_back( h_SL_MetPhi_AfterAllSelection_mu = new TH1D("h_SL_MetPhi_AfterAllSelection_mu","h_SL_MetPhi_AfterAllSelection_mu",  40, -M_PI , M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi, mu channel");


  histogram_list.push_back( h_SL_MetType1xyPhi_AfterAllSelection    = new TH1D("h_SL_MetType1xyPhi_AfterAllSelection"   ,"h_SL_MetType1xyPhi_AfterAllSelection",     40, -M_PI ,    M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi(type1xy)");
  histogram_list.push_back( h_SL_MetType1xyPhi_AfterAllSelection_el = new TH1D("h_SL_MetType1xyPhi_AfterAllSelection_el","h_SL_MetType1xyPhi_AfterAllSelection_el",  40, -M_PI , M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi(type1xy), el channel");
  histogram_list.push_back( h_SL_MetType1xyPhi_AfterAllSelection_mu = new TH1D("h_SL_MetType1xyPhi_AfterAllSelection_mu","h_SL_MetType1xyPhi_AfterAllSelection_mu",  40, -M_PI , M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met phi(type1xy), mu channel");


  histogram_list.push_back( h_SL_MetEx_AfterAllSelection    = new TH1D("h_SL_MetEx_AfterAllSelection"   ,"h_SL_MetEx_AfterAllSelection",     40, -200, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met Ex [GeV]");
  histogram_list.push_back( h_SL_MetEx_AfterAllSelection_el = new TH1D("h_SL_MetEx_AfterAllSelection_el","h_SL_MetEx_AfterAllSelection_el",  40, -200, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met Ex [GeV], el channel");
  histogram_list.push_back( h_SL_MetEx_AfterAllSelection_mu = new TH1D("h_SL_MetEx_AfterAllSelection_mu","h_SL_MetEx_AfterAllSelection_mu",  40, -200, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met Ex [GeV], mu channel");

  histogram_list.push_back( h_SL_MetEy_AfterAllSelection    = new TH1D("h_SL_MetEy_AfterAllSelection"   ,"h_SL_MetEy_AfterAllSelection",     40, -200, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met Ey [GeV]");
  histogram_list.push_back( h_SL_MetEy_AfterAllSelection_el = new TH1D("h_SL_MetEy_AfterAllSelection_el","h_SL_MetEy_AfterAllSelection_el",  40, -200, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met Ey [GeV], el channel");
  histogram_list.push_back( h_SL_MetEy_AfterAllSelection_mu = new TH1D("h_SL_MetEy_AfterAllSelection_mu","h_SL_MetEy_AfterAllSelection_mu",  40, -200, 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Met Ey [GeV], mu channel");




  histogram_list.push_back( h_SL_MuonPtAfterAllSelection = new TH1D("SL_MuonPT_AfterAllSelection","SL_MuonPT_AfterAllSelection", 28, 20, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon PT[GeV]");
  histogram_list.push_back( h_SL_ElPtAfterAllSelection   = new TH1D("SL_ElPT_AfterAllSelection","SL_ElPT_AfterAllSelection",  29, 30, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron PT[GeV]");


  histogram_list.push_back( h_SL_ElPtAfterAllSelection_ge3btag   = new TH1D("SL_ElPT_AfterAllSelection_ge3btag","SL_ElPT_AfterAllSelection_ge3btag",  29, 30, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron PT[GeV] (ge 3btag jet)");
  histogram_list.push_back( h_SL_MuonPtAfterAllSelection_ge3btag = new TH1D("SL_MuonPT_AfterAllSelection_ge3btag","SL_MuonPT_AfterAllSelection_ge3btag", 28, 20, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon PT[GeV] (ge 3 btag jet");



  histogram_list.push_back( h_SL_ElPtAfterAllSelection_noElSF   = new TH1D("SL_ElPT_AfterAllSelection_noElSF","SL_ElPT_AfterAllSelection_noElSF",  100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron PT[GeV], no EL SFs");

  histogram_list.push_back( h_SL_MuonPt_NoTopPTWgt_AfterAllSelection = new TH1D("SL_MuonPT_NoTopPTWgt_AfterAllSelection","SL_MuonPT_NoTopPTWgt_AfterAllSelection", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon PT[GeV], no top reweight SF");

  histogram_list.push_back( h_SL_ElPt_NoTopPTWgt_AfterAllSelection   = new TH1D("SL_ElPT_NoTopPTWgt_AfterAllSelection","SL_ElPT_NoTopPTWgt_AfterAllSelection",  100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron PT[GeV], no top reweight SF");

  histogram_list.push_back( h_SL_MuonEtaAfterAllSelection = new TH1D("SL_MuonEtaAfterAllSelection","SL_MuonEta_AfterAllSelection",  50 , -2.5, 2.5 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon Eta");

  histogram_list.push_back( h_SL_MuonPhiAfterAllSelection = new TH1D("SL_MuonPhiAfterAllSelection","SL_MuonPhi_AfterAllSelection",  50 , -M_PI, M_PI ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon Phi");

  histogram_list.push_back( h_SL_ElEtaAfterAllSelection   = new TH1D("SL_ElEta_AfterAllSelection","SL_ElEta_AfterAllSelection",   50 , -2.5, 2.5 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron Eta");

  histogram_list.push_back( h_SL_ElPhiAfterAllSelection   = new TH1D("SL_ElPhi_AfterAllSelection","SL_ElPhi_AfterAllSelection",   50 , -M_PI, +M_PI ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron Phi");

  histogram_list.push_back( h_SL_ElEtaAfterAllSelection_noElSF   = new TH1D("SL_ElEta_AfterAllSelection_noElSF","SL_ElEta_AfterAllSelection_noElSF",   50 , -2.5, 2.5 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron Eta, no EL SFs");

  histogram_list.push_back( h_SL_MuonEta_NoTopPTWgt_AfterAllSelection = new TH1D("SL_MuonEta_NoTopPTWgt_AfterAllSelection","SL_MuonEta_NoTopPTWgt_AfterAllSelection",  50 , -2.5, 2.5 ));
  histogram_list.push_back( h_SL_ElEta_NoTopPTWgt_AfterAllSelection   = new TH1D("SL_ElEta_NoTopPTWgt_AfterAllSelection","SL_ElEta_NoTopPTWgt_AfterAllSelection",   50 , -2.5, 2.5 ));


  histogram_list.push_back( h_SL_MuonPt_NoMuTrigWgt_AfterAllSelection = new TH1D("SL_MuonPT_NoMuTrigWgt_AfterAllSelection","SL_MuonPT_NoMuTrigWgt_AfterAllSelection", 100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon PT[GeV], no mu trig SF");
  histogram_list.push_back( h_SL_MuonEta_NoMuTrigWgt_AfterAllSelection = new TH1D("SL_MuonEta_NoMuTrigWgt_AfterAllSelection","SL_MuonEta_NoMuTrigWgt_AfterAllSelection",  50 , -2.5, 2.5 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon Eta, no mu trig SF");
  histogram_list.push_back( h_SL_ElEta_NoElTrigWgt_AfterAllSelection   = new TH1D("SL_ElEta_NoElTrigWgt_AfterAllSelection","SL_ElEta_NoElTrigWgt_AfterAllSelection",   50 , -2.5, 2.5 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron Eta, no el trig SF");
  histogram_list.push_back( h_SL_ElPt_NoElTrigWgt_AfterAllSelection   = new TH1D("SL_ElPT_NoElTrigWgt_AfterAllSelection","SL_ElPT_NoElTrigWgt_AfterAllSelection",  100, 0, 300) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron PT[GeV], no el trig SF");

  histogram_list.push_back( h_SL_nJet_AfterAllSelection   = new TH1D("SL_nJet_AfterAllSelection","SL_nJet_AfterAllSelection",  10 , 4, 14 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("nJet");

  histogram_list.push_back( h_SL_nJet_AfterAllSelection_el   = new TH1D("SL_nJet_AfterAllSelection_el","SL_nJet_AfterAllSelection_el",  10 , 4, 14 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("nJet, el channel");
  histogram_list.push_back( h_SL_nJet_AfterAllSelection_mu  = new TH1D("SL_nJet_AfterAllSelection_mu","SL_nJet_AfterAllSelection_mu",  10 , 4, 14 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("nJet, mu channel");

  histogram_list.push_back( h_SL_nJet_NoTopPTWgt_AfterAllSelection   = new TH1D("SL_nJet_NoTopPTWgt_AfterAllSelection","SL_nJet_NoTopPTWgt_AfterAllSelection",  10 , 4, 14 ));
  histogram_list.push_back( h_SL_nJetWITHOUTPUWGT_AfterAllSelection   = new TH1D("SL_nJetWITHOUTPUWGT_AfterAllSelection","SL_nJetWITHOUTPUWGT_AfterAllSelection",  10 , 4, 14 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("nJet, no PU weight");


  histogram_list.push_back( h_SL_MuonPT_AllButJetRequirement  = new TH1D("SL_MuonPT_AllButJetRequirement","SL_MuonPT_AllButJetRequirement",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon PT[GeV]");
  histogram_list.push_back( h_SL_ElPT_AllButJetRequirement  = new TH1D("SL_ElPT_AllButJetRequirement","SL_ElPT_AllButJetRequirement",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron PT[GeV]");
  histogram_list.push_back( h_SL_Met_AllButJetRequirement_el  = new TH1D("SL_Met_AllButJetRequirement_el","SL_Met_AllButJetRequirement_el",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV], el-ch");

  histogram_list.push_back( h_SL_Met_AllButJetRequirement_mu  = new TH1D("SL_Met_AllButJetRequirement_mu","SL_Met_AllButJetRequirement_mu",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV], mu-ch");


  histogram_list.push_back( h_SL_MET_AllButJetRequirement_jet0_mu =  new TH1D("SL_Met_AllButJetRequirement_jet0_mu","SL_Met_AllButJetRequirement_jet0_mu",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV], 0jet, mu-ch");

  histogram_list.push_back( h_SL_MET_AllButJetRequirement_jet1_mu =  new TH1D("SL_Met_AllButJetRequirement_jet1_mu","SL_Met_AllButJetRequirement_jet1_mu",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV], 1jet, mu-ch");

  histogram_list.push_back( h_SL_MET_AllButJetRequirement_jet2_mu =  new TH1D("SL_Met_AllButJetRequirement_jet2_mu","SL_Met_AllButJetRequirement_jet2_mu",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV], 2jet, mu-ch");

  histogram_list.push_back( h_SL_MET_AllButJetRequirement_jet3_mu =  new TH1D("SL_Met_AllButJetRequirement_jet3_mu","SL_Met_AllButJetRequirement_jet3_mu",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV], 3jet, mu-ch");

  histogram_list.push_back( h_SL_MET_AllButJetRequirement_jet4_mu =  new TH1D("SL_Met_AllButJetRequirement_jet4_mu","SL_Met_AllButJetRequirement_jet4_mu",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV], 4jet, mu-ch");

  histogram_list.push_back( h_SL_MET_AllButJetRequirement_jet5_mu =  new TH1D("SL_Met_AllButJetRequirement_jet5_mu","SL_Met_AllButJetRequirement_jet5_mu",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV], 5/more jet, mu-ch");


  histogram_list.push_back( h_SL_MetType1xy_AllButJetRequirement_el  = new TH1D("SL_MetType1xy_AllButJetRequirement_el","SL_MetType1xy_AllButJetRequirement_el",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET Type1xy[GeV], el-ch");
  histogram_list.push_back( h_SL_MetType1xy_AllButJetRequirement_mu  = new TH1D("SL_MetType1xy_AllButJetRequirement_mu","SL_MetType1xy_AllButJetRequirement_mu",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET Type1xy[GeV], mu-ch");

  histogram_list.push_back( h_SL_MuonPT_AllButNBJetRequirement  = new TH1D("SL_MuonPT_AllButNBJetRequirement","SL_MuonPT_AllButNBJetRequirement",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon PT[GeV]");
  histogram_list.push_back( h_SL_ElPT_AllButNBJetRequirement  = new TH1D("SL_ElPT_AllButNBJetRequirement","SL_ElPT_AllButNBJetRequirement",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron PT[GeV]");
  histogram_list.push_back( h_SL_Met_AllButNBJetRequirement_el  = new TH1D("SL_Met_AllButNBJetRequirement_el","SL_Met_AllButNBJetRequirement_el",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV], el-ch");
  histogram_list.push_back( h_SL_Met_AllButNBJetRequirement_mu  = new TH1D("SL_Met_AllButNBJetRequirement_mu","SL_Met_AllButNBJetRequirement_mu",  30 , 0, 300 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV], mu-ch");


  histogram_list.push_back( h_SL_nJet_AllButJetRequirement  = new TH1D("SL_nJet_AllButJetRequirement","SL_nJet_AllButJetRequirement",  10 , 0, 10 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("nJet");
  histogram_list.push_back( h_SL_nJet_mu_AllButJetRequirement  = new TH1D("SL_nJet_mu_AllButJetRequirement","SL_nJet_mu_AllButJetRequirement",  10 , 0, 10 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("nJet, mu channel");
  histogram_list.push_back( h_SL_nJet_el_AllButJetRequirement  = new TH1D("SL_nJet_el_AllButJetRequirement","SL_nJet_el_AllButJetRequirement",  10 , 0, 10 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("nJet, el channel");

  histogram_list.push_back( h_SL_nBtagJet_AfterAllSelection   = new TH1D("SL_nBtagJet_AfterAllSelection","SL_nBtagJet_AfterAllSelection",  7 , 0, 7 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("n-bjet");

  histogram_list.push_back( h_SL_nBtagJet_AfterAllSelection_el   = new TH1D("SL_nBtagJet_AfterAllSelection_el","SL_nBtagJet_AfterAllSelection_el",  7 , 0, 7 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("n-bjet, el channel");
  histogram_list.push_back( h_SL_nBtagJet_AfterAllSelection_mu   = new TH1D("SL_nBtagJet_AfterAllSelection_mu","SL_nBtagJet_AfterAllSelection_mu",  7 , 0, 7 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("n-bjet, mu channel");
  histogram_list.push_back( h_SL_nBtagJetWITHOUTBTAGSF_AfterAllSelection   = new TH1D("SL_nBtagJetWITHOUTBTAGSF_AfterAllSelection","SL_nBtagJetWITHOUTBTAGSF_AfterAllSelection",  7, 0, 7 ));
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("n-bjet, without b-tag SF");




  // histogram within 6j 4b category :
  {

    histogram_list.push_back(  h_SL_Entries_nJnBCategoryAfterLeptonCut_MetCut_Mu
			       = new TH1D ( "SL_Entries_nJnBCategoryAfterLeptonCut_MetCut_Mu", "SL_Entries_nJnBCategoryAfterLeptonCut_MetCut_Mu",
					    26 , 0.5 , 26.5 ) ) ; 
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( "nJmB category, Met>20GeV" );

    char name[100];
    
    for( int j = 0 ; j <= 6 ; j ++ ) {
      for( int b = 0 ; b <= j && b <= 4 ; b++ ){
	sprintf( name , "SL_MuonPT_category%dj%db_AfterLeptonAndMetCut" , j , b ); 
	histogram_list.push_back(  h_SL_MuonPT_nJnBCategoryAfterLeptonCut_MetCut [j][b] 
				   = new TH1D ( name , name , 20, 0 , 200 ) ) ;  
	sprintf( name , "Muon PT[GeV](%dj%db, Met>20GeV)", j , b ) ; 
	histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( name );
	
	// ----- 

	sprintf( name , "SL_MET_category%dj%db_AfterLepton" , j , b ); 
	histogram_list.push_back(  h_SL_Met_nJnBCategoryAfterLeptonCut_Mu [j][b] 
				   = new TH1D ( name , name , 20, 0 , 200 ) ) ;  
	sprintf( name , "Met[GeV](%dj%db)", j , b ) ; 
	histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( name );
      }
    }
  }



  const double pt_range[8] = {
    400, 400, 300, 200, 120, 120, 80, 80
  };
  for( int i = 0 ; i < 8 ; i++){
    char name [100];
    sprintf(name , "SL_jet%d_pt_AfterAllSelection", i +1 );
    histogram_list.push_back( h_SL_jet_pt_AfterAllSelection [i]  = new TH1D( name , name , 20, 30, pt_range[i] ));
    sprintf(name , "%d-th Leading jet PT[GeV]", i +1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf(name , "SL_jet%d_pt_AfterAllSelection_ge3btag", i +1 );
    histogram_list.push_back( h_SL_jet_pt_AfterAllSelection_ge3btag [i]  = new TH1D( name , name , 20, 30, pt_range[i] ));
    sprintf(name , "%d-th Leading jet PT[GeV] (ge3btag)", i +1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf(name , "SL_jet%d_pt_WoTopPTegt_AfterAllSelection", i +1 );
    histogram_list.push_back( h_SL_jet_pt_WOtopPTwgt_AfterAllSelection [i]  = new TH1D( name , name , 20, 30 , pt_range[i] ));

    sprintf(name , "SL_jet%d_pt_AfterAllButBtagJetRequirement", i +1 );
    histogram_list.push_back( h_SL_jet_pt_AfterAllButBtagJetRequirement [i]  = new TH1D( name , name , 20, 30 , pt_range[i] ));
    sprintf(name , "%d-th Leading jet PT[GeV]", i +1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    sprintf(name , "SL_jet%d_pt_AfterAllSelection_el", i +1 );
    histogram_list.push_back( h_SL_jet_pt_AfterAllSelection_el [i]  = new TH1D( name , name , 20, 30, pt_range[i] ));
    sprintf(name , "%d-th Leading jet PT[GeV], el-ch", i +1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    sprintf(name , "SL_jet%d_pt_AfterAllSelection_mu", i +1 );
    histogram_list.push_back( h_SL_jet_pt_AfterAllSelection_mu [i]  = new TH1D( name , name , 20, 30, pt_range[i] ));
    sprintf(name , "%d-th Leading jet PT[GeV], mu-ch", i +1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);
  }
  for( int i = 0 ; i < 8 ; i++){
    char name [100];
    sprintf(name , "SL_jet%d_eta_AfterAllSelection", i +1 );
    histogram_list.push_back( h_SL_jet_eta_AfterAllSelection [i]  = new TH1D( name , name , 25, -2.5, 2.5 ) );
    sprintf(name , "%d-th Leading jet Eta", i +1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    sprintf(name , "SL_jet%d_eta_WoTopPTegt_AfterAllSelection", i +1 );
    histogram_list.push_back( h_SL_jet_eta_WOtopPTwgt_AfterAllSelection [i]  = new TH1D( name , name , 25, -2.5, 2.5 ) );

    sprintf(name , "SL_jet%d_eta_AfterAllSelection_el", i +1 );
    histogram_list.push_back( h_SL_jet_eta_AfterAllSelection_el [i]  = new TH1D( name , name , 25, -2.5, 2.5 ) );
    sprintf(name , "%d-th Leading jet Eta, el-ch", i +1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    sprintf(name , "SL_jet%d_eta_AfterAllSelection_mu", i +1 );
    histogram_list.push_back( h_SL_jet_eta_AfterAllSelection_mu [i]  = new TH1D( name , name , 25, -2.5, 2.5 ) );
    sprintf(name , "%d-th Leading jet Eta, mu-ch", i +1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

  }
  for( int iJet = 0 ; iJet < 8 ; iJet++ ){
    char name[100];
    sprintf( name, "SL_btagDiscriminant_%dthLeading_AfterAllSelection", iJet + 1 );
    histogram_list.push_back( h_SL_jetBtag_AfterAllSelection[ iJet ] = new TH1D( name , name , 20, 0 , 1) ) ;
    sprintf(name , "b-discri of %d-th leading jet", iJet + 1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    sprintf( name, "SL_btagDiscriminantWOBtagSF_%dthLeading_AfterAllSelection", iJet + 1 );
    histogram_list.push_back( h_SL_jetBtagWOBtagSF_AfterAllSelection[ iJet ] = new TH1D( name , name , 20, 0 , 1) ) ;
    sprintf(name , "b-discri of %d-th leading jet. nobtagSF", iJet + 1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    sprintf( name, "SL_btagDiscriminant_%dthLeading_AfterAllButBtagJetRequirement", iJet + 1 );
    histogram_list.push_back( h_SL_jetBtag_AfterAllButBtagJetRequirement[ iJet ] = new TH1D( name , name , 20, 0 , 1) ) ;
    sprintf(name , "b-discri of %d-th leading jet", iJet + 1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name, "SL_btagDiscriminant_%dthLeading_AfterAllSelection_el", iJet + 1 );
    histogram_list.push_back( h_SL_jetBtag_AfterAllSelection_el[ iJet ] = new TH1D( name , name , 20, 0 , 1) ) ;
    sprintf(name , "b-discri of %d-th leading jet, el-ch", iJet + 1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    sprintf( name, "SL_btagDiscriminant_%dthLeading_AfterAllSelection_mu", iJet + 1 );
    histogram_list.push_back( h_SL_jetBtag_AfterAllSelection_mu[ iJet ] = new TH1D( name , name , 20, 0 , 1) ) ;
    sprintf(name , "b-discri of %d-th leading jet, mu-ch", iJet + 1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


    sprintf( name, "SL_btagDiscriminantWOBtagSF_%dthLeading_AfterAllSelection_el", iJet + 1 );
    histogram_list.push_back( h_SL_jetBtagWOBtagSF_AfterAllSelection_el[ iJet ] = new TH1D( name , name , 20, 0 , 1) ) ;
    sprintf(name , "b-discri of %d-th leading jet, without btag SF, el-ch", iJet + 1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    sprintf( name, "SL_btagDiscriminantWOBtagSF_%dthLeading_AfterAllSelection_mu", iJet + 1 );
    histogram_list.push_back( h_SL_jetBtagWOBtagSF_AfterAllSelection_mu[ iJet ] = new TH1D( name , name , 20, 0 , 1) ) ;
    sprintf(name , "b-discri of %d-th leading jet, without btag SF, mu-ch", iJet + 1 );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

  }



  histogram_list.push_back(  h_SL_jetPt_upto8th_BtagLFControlRegion  =  new TH1D ( "SL_jetPt_upto8th_BtagLFControlRegion","SL_jetPt_upto8th_BtagLFControlRegion", 80, 0, 400  ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Jet PT[GeV], 1-8th leading jets. LF");
  histogram_list.push_back(  h_SL_jetPt_upto8th_BtagHFControlRegion  =  new TH1D ( "SL_jetPt_upto8th_BtagHFControlRegion","SL_jetPt_upto8th_BtagHFControlRegion", 80, 0, 400  ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Jet PT[GeV], 1-8th leading jets. HF");

  histogram_list.push_back(  h_SL_jetBtag_upto8th_BtagLFControlRegion = new TH1D( "SL_jetBtag_upto8th_BtagLFControlRegion","SL_jetBtag_upto8th_BtagLFControlRegion",20, 0 , 1 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("b-discri, 1-8th leading jets. LF");
  histogram_list.push_back(  h_SL_jetBtag_upto8th_BtagHFControlRegion = new TH1D( "SL_jetBtag_upto8th_BtagHFControlRegion","SL_jetBtag_upto8th_BtagHFControlRegion",20, 0 , 1 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("b-discri, 1-8th leading jets. HF");

  histogram_list.push_back(  h_SL_jetBtag_upto8th_BtagLFControlRegionWOBtagSF =new TH1D( "SL_jetBtag_upto8th_BtagLFControlRegionWOBtagSF","SL_jetBtag_upto8th_BtagLFControlRegionWOBtagSF",20, 0 , 1 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("b-discri, 1-8th leading jets. HF, noBtagSF");
  histogram_list.push_back(  h_SL_jetBtag_upto8th_BtagHFControlRegionWOBtagSF =new TH1D( "SL_jetBtag_upto8th_BtagHFControlRegionWOBtagSF","SL_jetBtag_upto8th_BtagHFControlRegionWOBtagSF",20, 0 , 1 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("b-discri, 1-8th leading jets. HF, noBtagSF");


  for( int cate = 1 ; cate < 10 ; cate++  ){

    std::string catename = _NameofCategory( cate );
    char histname [100];

    sprintf( histname , "SL_ElPT_Category%s" , catename.c_str() );
    histogram_list.push_back(h_SL_ElPt_AfterNjNbCategory[cate]     = new TH1D(histname,histname, 20, 0 , 300) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron PT[GeV] (Met>20GeV)");    

    sprintf( histname , "SL_ElEta_Category%s" ,  catename.c_str() );
    histogram_list.push_back(h_SL_ElEta_AfterNjNbCategory[cate]     = new TH1D(histname,histname, 25, -2.5 , 2.5 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Electron Eta (Met>20GeV)");    

    sprintf( histname , "SL_MuPT_Category%s" , catename.c_str() );
    histogram_list.push_back(h_SL_MuPt_AfterNjNbCategory[cate]     = new TH1D(histname,histname, 20, 0 , 300) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon PT[GeV] (Met>20GeV)");        

    sprintf( histname , "SL_MuEta_Category%s" ,  catename.c_str() );
    histogram_list.push_back(h_SL_MuEta_AfterNjNbCategory[cate]     = new TH1D(histname,histname, 25, -2.5 , 2.5 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Muon Eta (Met>20GeV)");    

    sprintf( histname , "SL_MinElJetDr_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back(h_SL_MinElJetDr_AfterNjNbCategory[cate]     = new TH1D(histname,histname, 20, 0 , 2.0) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MinDR(EL,Jet) (Met>20GeV)");    

    sprintf( histname , "SL_MinMuJetDr_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back(h_SL_MinMuJetDr_AfterNjNbCategory[cate]     = new TH1D(histname,histname, 20, 0 , 2.0) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MinDR(Mu,Jet) (Met>20GeV)");    




    sprintf( histname , "SL_El_LeadingJetPT_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back(h_SL_El_LeadingJetPT_AfterNjNbCategory  [cate]     = new TH1D(histname,histname, 20, 0 , 300) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(El ch) Leading jet PT[GeV]");    


    sprintf( histname , "SL_El_LeastJetPT_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_El_LeastJetPT_AfterNjNbCategory[cate]     = new TH1D(histname,histname, 15, 0 , 150) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(El ch) Least Leading jet PT[GeV]");    


    sprintf( histname , "SL_El_LeadingJetBtag_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_El_LeadingJetBtag_AfterNjNbCategory[cate]     = new TH1D(histname,histname,40, 0 , 1 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(El ch) Leading Jet B-tagger");    

    sprintf( histname , "SL_El_LeastJetBtag_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_El_LeastJetBtag_AfterNjNbCategory[cate]     = new TH1D(histname,histname,40, 0 , 1 ) ) ;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(El ch) Least Leading jet B-tagger");    


    sprintf( histname , "SL_El_LargestBtagValue_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_El_LargestBtagValue_AfterNjNbCategory[cate]     = new TH1D(histname,histname,40, 0 , 1) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(El ch) Biggest B-tagger");    

    sprintf( histname , "SL_El_LeastBtagValue_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_El_LeastBtagValue_AfterNjNbCategory[cate]     = new TH1D(histname,histname,40, 0 , 1) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(El ch) Smallest B-tagger");    


    long nBin[10];
    nBin [1] = 20 ; // "4j2b" ;
    nBin [2] = 15 ; // "4j3b" ;
    nBin [3] = 10 ; // "4j4b" ;
    nBin [4] = 20 ; // "5j2b" ;
    nBin [5] = 15 ; // "5j3b" ;
    nBin [6] = 10 ; // "5j4b" ;
    nBin [7] = 20 ; // "6j2b" ;
    nBin [8] = 15 ; // "6j3b" ;
    nBin [9] = 10 ; // "6j4b" ;


    sprintf( histname , "SL_Mu_LeadingJetPT_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back(h_SL_Mu_LeadingJetPT_AfterNjNbCategory  [cate]     = new TH1D(histname,histname, nBin[cate], 30 , 200) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(Mu ch) Leading jet PT[GeV]");    

    sprintf( histname , "SL_Mu_LeastJetPT_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_Mu_LeastJetPT_AfterNjNbCategory[cate]     = new TH1D(histname,histname, nBin[cate], 30 , 150) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(Mu ch) Least Leading jet PT[GeV]");    


    sprintf( histname , "SL_Mu_LeadingJetBtag_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_Mu_LeadingJetBtag_AfterNjNbCategory[cate]     = new TH1D(histname,histname,40, 0 , 1 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(Mu ch) Leading Jet B-tagger");    

    sprintf( histname , "SL_Mu_LeastJetBtag_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_Mu_LeastJetBtag_AfterNjNbCategory[cate]     = new TH1D(histname,histname,40, 0 , 1 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(Mu ch) Least Leading jet B-tagger");    


    sprintf( histname , "SL_Mu_LargestBtagValue_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_Mu_LargestBtagValue_AfterNjNbCategory[cate]     = new TH1D(histname,histname,40, 0 , 1) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(Mu ch) Biggest B-tagger");    

    sprintf( histname , "SL_Mu_LeastBtagValue_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_Mu_LeastBtagValue_AfterNjNbCategory[cate]     = new TH1D(histname,histname,40, 0 , 1) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("(Mu ch) Smallest B-tagger");    



//disableCategoryPlot    for( int i = 0 ; i < 6 ; i ++ ){
//disableCategoryPlot
//disableCategoryPlot      sprintf( histname , "SL_JetPT_%d_Category%s" , i+1 ,  catename.c_str() );
//disableCategoryPlot      histogram_list.push_back(h_SL_JetPt_AfterNjNbCategory[cate][i]     = new TH1D(histname,histname, 20, 0 , 300) );
//disableCategoryPlot
//disableCategoryPlot      sprintf( histname , "SL_JetBtag_%d_Category%s" , i+1 ,  catename.c_str() );
//disableCategoryPlot      histogram_list.push_back( h_SL_JetBtag_AfterNjNbCategory[cate][i] = new TH1D(histname,histname, 20, 0 , 1 ) );
//disableCategoryPlot
//disableCategoryPlot      sprintf( histname , "SL_BtagInBtagOrder_%d_Category%s" , i+1 ,  catename.c_str() );
//disableCategoryPlot      histogram_list.push_back( h_SL_Btagorder_AfterNjNbCategory[cate][i] = new TH1D(histname,histname, 20, 0 , 1 ) ); 
//disableCategoryPlot
//disableCategoryPlot      sprintf( histname , "SL_JetPT_%d_Category%s_el" , i+1 ,  catename.c_str() );
//disableCategoryPlot      histogram_list.push_back(h_SL_JetPt_AfterNjNbCategory_el[cate][i]     = new TH1D(histname,histname, 20, 0 , 300) );
//disableCategoryPlot
//disableCategoryPlot      sprintf( histname , "SL_JetBtag_%d_Category%s_el" , i+1 ,  catename.c_str() );
//disableCategoryPlot      histogram_list.push_back( h_SL_JetBtag_AfterNjNbCategory_el[cate][i] = new TH1D(histname,histname, 20, 0 , 1 ) );
//disableCategoryPlot
//disableCategoryPlot      sprintf( histname , "SL_BtagInBtagOrder_%d_Category%s_el" , i+1 ,  catename.c_str() );
//disableCategoryPlot      histogram_list.push_back( h_SL_Btagorder_AfterNjNbCategory_el[cate][i] = new TH1D(histname,histname, 20, 0 , 1 ) ); 
//disableCategoryPlot
//disableCategoryPlot      sprintf( histname , "SL_JetPT_%d_Category%s_mu" , i+1 ,  catename.c_str() );
//disableCategoryPlot      histogram_list.push_back(h_SL_JetPt_AfterNjNbCategory_mu[cate][i]     = new TH1D(histname,histname, 20, 0 , 300) );
//disableCategoryPlot
//disableCategoryPlot      sprintf( histname , "SL_JetBtag_%d_Category%s_mu" , i+1 ,  catename.c_str() );
//disableCategoryPlot      histogram_list.push_back( h_SL_JetBtag_AfterNjNbCategory_mu[cate][i] = new TH1D(histname,histname, 20, 0 , 1 ) );
//disableCategoryPlot
//disableCategoryPlot      sprintf( histname , "SL_BtagInBtagOrder_%d_Category%s_mu" , i+1 ,  catename.c_str() );
//disableCategoryPlot      histogram_list.push_back( h_SL_Btagorder_AfterNjNbCategory_mu[cate][i] = new TH1D(histname,histname, 20, 0 , 1 ) ); 
//disableCategoryPlot
//disableCategoryPlot    } // jet loop



    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_pt"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_pt_NbNjCate [cate] = new TH1D( histname  , histname ,  10, 200, 1000 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet PT[GeV]"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_eta"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_eta_NbNjCate [cate] = new TH1D( histname  , histname ,  10, -4, 4 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet #eta"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_phi"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_phi_NbNjCate [cate] = new TH1D( histname  , histname ,  10, -M_PI, M_PI ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet #phi"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_sdmass"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_sdmass_NbNjCate [cate] = new TH1D( histname  , histname ,  10, 0 , 200 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet SDmass[GeV]"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_loosetagged"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_loosetagged_NbNjCate [cate] = new TH1D( histname  , histname ,  5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet n-loose-b"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_mediumtagged"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_mediumtagged_NbNjCate [cate] = new TH1D( histname  , histname ,  5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet n-medim-b"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_meanbtagger"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_mean_btagger_NbNjCate [cate] = new TH1D( histname  , histname ,   10, 0 , 1 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet mean b-tagger"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_Nfat"  , catename.c_str() );
    histogram_list.push_back( h_SL_N_fatjet_NbNjCate [cate] = new TH1D( histname  , histname ,   5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s" , "N AK8 jets"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_tau21"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_tau21_NbNjCate [cate] = new TH1D( histname  , histname ,  5, 0 , 1 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet tau21"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_tau32"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_tau32_NbNjCate [cate] = new TH1D( histname  , histname ,  5, 0 , 1 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet tau32"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);



    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_NLooseBBfat"  , catename.c_str() );
    histogram_list.push_back( h_SL_N_LooseBBTaggedfatjet_NbNjCate [cate] = new TH1D( histname  , histname ,   5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s" , "N AK8 LooseBB jets"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_AfterNjNbCategory%s" , "SL_fatjet_NMediumBBfat"  , catename.c_str() );
    histogram_list.push_back( h_SL_N_MediumBBTaggedfatjet_NbNjCate [cate] = new TH1D( histname  , histname ,   5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s" , "N AK8 MedBB jets"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);






    sprintf( histname , "%s_SDmassCut_AfterNjNbCategory%s" , "SL_fatjet_pt"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_pt_NbNjCate_SDmassCut [cate] = new TH1D( histname  , histname ,  10, 200, 1000 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet PT[GeV]"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_SDmassCut_AfterNjNbCategory%s" , "SL_fatjet_eta"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_eta_NbNjCate_SDmassCut [cate] = new TH1D( histname  , histname ,  10, -4, 4 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet #eta"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_SDmassCut_AfterNjNbCategory%s" , "SL_fatjet_phi"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_phi_NbNjCate_SDmassCut [cate] = new TH1D( histname  , histname ,  10, -M_PI, M_PI ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet #phi"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_SDmassCut_AfterNjNbCategory%s" , "SL_fatjet_loosetagged"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_loosetagged_NbNjCate_SDmassCut [cate] = new TH1D( histname  , histname ,  5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet n-loose-b"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_SDmassCut_AfterNjNbCategory%s" , "SL_fatjet_mediumtagged"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_mediumtagged_NbNjCate_SDmassCut [cate] = new TH1D( histname  , histname ,  5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet n-medim-b"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_SDmassCut_AfterNjNbCategory%s" , "SL_fatjet_meanbtagger"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_mean_btagger_NbNjCate_SDmassCut [cate] = new TH1D( histname  , histname ,   10, 0 , 1 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet mean b-tagger"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_SDmassCut_AfterNjNbCategory%s" , "SL_fatjet_Nfat"  , catename.c_str() );
    histogram_list.push_back( h_SL_N_fatjet_NbNjCate_SDmassCut [cate] = new TH1D( histname  , histname ,   5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s" , "N AK8 jets"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_SDmassCut_AfterNjNbCategory%s" , "SL_fatjet_tau21"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_tau21_NbNjCate_SDmassCut [cate] = new TH1D( histname  , histname ,  5, 0 , 1 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet tau21"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_SDmassCut_AfterNjNbCategory%s" , "SL_fatjet_tau32"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_tau32_NbNjCate_SDmassCut [cate] = new TH1D( histname  , histname ,  5, 0 , 1 ) );
    sprintf( histname , "%s, Cate%s" , "AK8 Jet tau32"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);





    sprintf( histname , "%s_el_AfterNjNbCategory%s" , "SL_fatjet_pt"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_pt_NbNjCate_el [cate] = new TH1D( histname  , histname ,  10, 200, 1000 ) );
    sprintf( histname , "%s, Cate%s, el" , "AK8 Jet PT[GeV]"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_el_AfterNjNbCategory%s" , "SL_fatjet_eta"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_eta_NbNjCate_el [cate] = new TH1D( histname  , histname ,  10, -4, 4 ) );
    sprintf( histname , "%s, Cate%s, el" , "AK8 Jet #eta"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_el_AfterNjNbCategory%s" , "SL_fatjet_phi"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_phi_NbNjCate_el [cate] = new TH1D( histname  , histname ,  10, -M_PI, M_PI ) );
    sprintf( histname , "%s, Cate%s, el" , "AK8 Jet #phi"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_el_AfterNjNbCategory%s" , "SL_fatjet_sdmass"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_sdmass_NbNjCate_el [cate] = new TH1D( histname  , histname ,  10, 0 , 200 ) );
    sprintf( histname , "%s, Cate%s, el" , "AK8 Jet SDmass[GeV]"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_el_AfterNjNbCategory%s" , "SL_fatjet_loosetagged"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_loosetagged_NbNjCate_el [cate] = new TH1D( histname  , histname ,  5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s, el" , "AK8 Jet n-loose-b"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_el_AfterNjNbCategory%s" , "SL_fatjet_mediumtagged"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_mediumtagged_NbNjCate_el [cate] = new TH1D( histname  , histname ,  5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s, el" , "AK8 Jet n-medim-b"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_el_AfterNjNbCategory%s" , "SL_fatjet_meanbtagger"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_mean_btagger_NbNjCate_el [cate] = new TH1D( histname  , histname ,   10, 0 , 1 ) );
    sprintf( histname , "%s, Cate%s, el" , "AK8 Jet mean b-tagger"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_el_AfterNjNbCategory%s" , "SL_fatjet_Nfat"  , catename.c_str() );
    histogram_list.push_back( h_SL_N_fatjet_NbNjCate_el [cate] = new TH1D( histname  , histname ,   5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, Cate%s, el" , "N AK8 jets"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_el_AfterNjNbCategory%s" , "SL_fatjet_tau21"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_tau21_NbNjCate_el [cate] = new TH1D( histname  , histname ,  5, 0 , 1 ) );
    sprintf( histname , "%s, Cate%s, el" , "AK8 Jet tau21"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_el_AfterNjNbCategory%s" , "SL_fatjet_tau32"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_tau32_NbNjCate_el [cate] = new TH1D( histname  , histname ,  5, 0 , 1 ) );
    sprintf( histname , "%s, Cate%s, el" , "AK8 Jet tau32"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);




    sprintf( histname , "%s_mu_AfterNjNbCategory%s" , "SL_fatjet_pt"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_pt_NbNjCate_mu [cate] = new TH1D( histname  , histname ,  10, 200, 1000 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet PT after all selection (mu)[GeV]");
    sprintf( histname , "%s, Cate%s, mu" , "AK8 Jet Pt"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_mu_AfterNjNbCategory%s" , "SL_fatjet_eta"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_eta_NbNjCate_mu [cate] = new TH1D( histname  , histname ,   10, -4, 4 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet #eta after all selection (mu)");
    sprintf( histname , "%s, Cate%s, mu" , "AK8 Jet #eta"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_mu_AfterNjNbCategory%s" , "SL_fatjet_phi"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_phi_NbNjCate_mu [cate] = new TH1D( histname  , histname ,  10, -M_PI, M_PI ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet #phi after all selection (mu)");
    sprintf( histname , "%s, Cate%s, mu" , "AK8 Jet #phi"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_mu_AfterNjNbCategory%s" , "SL_fatjet_sdmass"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_sdmass_NbNjCate_mu [cate] = new TH1D( histname  , histname ,  10, 0 , 200 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet SD mass after mu selection[GeV]");
    sprintf( histname , "%s, Cate%s, mu" , "AK8 Jet SDmass[GeV]"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_mu_AfterNjNbCategory%s" , "SL_fatjet_loosetagged"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_loosetagged_NbNjCate_mu [cate] = new TH1D( histname  , histname ,   5, -0.5 , 4.5 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("N of loose b-tagged subjets (after mu selection)");
    sprintf( histname , "%s, Cate%s, mu" , "AK8 Jet n-loose-b"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_mu_AfterNjNbCategory%s" , "SL_fatjet_mediumtagged"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_mediumtagged_NbNjCate_mu [cate] = new TH1D( histname  , histname ,  5, -0.5 , 4.5 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("N of medium b-tagged subjets (after mu selection)");
    sprintf( histname , "%s, Cate%s, mu" , "AK8 Jet n-medium-b"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_mu_AfterNjNbCategory%s" , "SL_fatjet_meanbtagger"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_subjet_mean_btagger_NbNjCate_mu [cate] = new TH1D( histname  , histname ,  10, 0 , 1 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mean of b-tagger of subjets (after mu selection)");
    sprintf( histname , "%s, Cate%s, mu" , "AK8 Jet mean b-tagger"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_mu_AfterNjNbCategory%s" , "SL_fatjet_Nfat"  , catename.c_str() );
    histogram_list.push_back( h_SL_N_fatjet_NbNjCate_mu [cate] = new TH1D( histname  , histname ,   5, -0.5 , 4.5 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("N of AK8 jets after all selection (mu)");
    sprintf( histname , "%s, Cate%s, mu" , "N AK8 jets"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_mu_AfterNjNbCategory%s" , "SL_fatjet_tau21"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_tau21_NbNjCate_mu [cate] = new TH1D( histname  , histname ,  5, 0 , 1 )) ;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("tau21 of AK8 jet(after mu selection)");
    sprintf( histname , "%s, Cate%s, mu" , "AK8 Jet tau21"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "%s_mu_AfterNjNbCategory%s" , "SL_fatjet_tau32"  , catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_tau32_NbNjCate_mu [cate] = new TH1D( histname  , histname ,   5, 0 , 1 )) ;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("tau32 of AK8 jet(after mu selection)");
    sprintf( histname , "%s, Cate%s, mu" , "AK8 Jet tau32"  , catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);



  }// cate loop



    

  {

    char histname [200];

    sprintf( histname , "SL_fatjet_%s_ofLargeSDmassJet_Cate6j4b_TwoFatJets" , "SDMass"  );;
    histogram_list.push_back( h_SL_fatjet_SDMass_ofLargeSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  10, 0 , 200 ) );
    sprintf( histname , "largeSDmassAk8Jet %s" , "SDMass[GeV]"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofSmallSDmassJet_Cate6j4b_TwoFatJets" , "SDMass"  );;
    histogram_list.push_back( h_SL_fatjet_SDMass_ofSmallSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  10, 0 , 200 ) );
    sprintf( histname , "SmallSDmassAk8Jet %s" , "SDMass[GeV]"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "SL_fatjet_%s_ofLargeSDmassJet_Cate6j4b_TwoFatJets" , "PT"  );;
    histogram_list.push_back( h_SL_fatjet_pt_ofLargeSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  10, 200 , 500 ) );
    sprintf( histname , "largeSDmassAk8Jet %s" , "PT[GeV]"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofSmallSDmassJet_Cate6j4b_TwoFatJets" , "PT"  );;
    histogram_list.push_back( h_SL_fatjet_pt_ofSmallSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  10, 200 , 500 ) );
    sprintf( histname , "SmallSDmassAk8Jet %s" , "PT[GeV]"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofLargeSDmassJet_Cate6j4b_TwoFatJets" , "eta"  );;
    histogram_list.push_back( h_SL_fatjet_eta_ofLargeSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  12 , -2.4 , 2.4 ) );
    sprintf( histname , "largeSDmassAk8Jet %s" , "#eta"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofSmallSDmassJet_Cate6j4b_TwoFatJets" , "eta"  );;
    histogram_list.push_back( h_SL_fatjet_eta_ofSmallSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  12 , -2.4 , 2.4 ) );
    sprintf( histname , "SmallSDmassAk8Jet %s" , "#eta"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);



    sprintf( histname , "SL_fatjet_%s_ofLargeSDmassJet_Cate6j4b_TwoFatJets" , "tau32"  );;
    histogram_list.push_back( h_SL_fatjet_tau32_ofLargeSDmassJet_Cate6j4b_TwoFatJets  = new TH1D( histname  , histname ,  10, 0 , 1 ) );
    sprintf( histname , "largeSDmassAk8Jet %s" , "tau32"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofSmallSDmassJet_Cate6j4b_TwoFatJets" , "tau32"  );;
    histogram_list.push_back( h_SL_fatjet_tau32_ofSmallSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  10, 0 , 1 ) );
    sprintf( histname , "SmallSDmassAk8Jet %s" , "tau32"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofLargeSDmassJet_Cate6j4b_TwoFatJets" , "tau21"  );;
    histogram_list.push_back( h_SL_fatjet_tau21_ofLargeSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  10, 0 , 1 ) );
    sprintf( histname , "largeSDmassAk8Jet %s" , "tau21"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofSmallSDmassJet_Cate6j4b_TwoFatJets" , "tau21"  );;
    histogram_list.push_back( h_SL_fatjet_tau21_ofSmallSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  10, 0 , 1 ) );
    sprintf( histname , "SmallSDmassAk8Jet %s" , "tau21"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "SL_fatjet_%s_ofLargeSDmassJet_Cate6j4b_TwoFatJets" , "NLooseB"  );;
    histogram_list.push_back( h_SL_fatjet_NlooseB_ofLargeSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  3, -0.5 , 2.5 ) );
    sprintf( histname , "largeSDmassAk8Jet %s" , "NLooseB"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofSmallSDmassJet_Cate6j4b_TwoFatJets" , "NLooseB"  );;
    histogram_list.push_back( h_SL_fatjet_NlooseB_ofSmallSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  3, -0.5 , 2.5 ) );
    sprintf( histname , "SmallSDmassAk8Jet %s" , "NLooseB"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofLargeSDmassJet_Cate6j4b_TwoFatJets" , "NMediumB"  );;
    histogram_list.push_back( h_SL_fatjet_NMediumB_ofLargeSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname , 3, -0.5 , 2.5 ) );
    sprintf( histname , "largeSDmassAk8Jet %s" , "NMediumB"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofSmallSDmassJet_Cate6j4b_TwoFatJets" , "NMediumB"  );;
    histogram_list.push_back( h_SL_fatjet_NMediumB_ofSmallSDmassJet_Cate6j4b_TwoFatJets = new TH1D( histname  , histname , 3, -0.5 , 2.5 ) );
    sprintf( histname , "SmallSDmassAk8Jet %s" , "NMediumB"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_ofLargeSDmassJet_Cate6j4b_TwoFatJets" , "dRLep"  );;
    histogram_list.push_back( h_SL_fatjet_dRlep_ofLargeSDmassJet_Cate6j4b_TwoFatJets  = new TH1D( histname  , histname ,  30 , 0 , 3 ) );
    sprintf( histname , "largeSDmassAk8Jet %s" , "dRLep"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "SL_fatjet_%s_ofSmallSDmassJet_Cate6j4b_TwoFatJets" , "dRLep"  );;
    histogram_list.push_back( h_SL_fatjet_dRlep_ofSmallSDmassJet_Cate6j4b_TwoFatJets  = new TH1D( histname  , histname ,  30 , 0 , 3 ) );
    sprintf( histname , "SmallSDmassAk8Jet %s" , "dRLep"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_AmbitiousSelection_Cate6j4b_TwoFatJet"  );;
    histogram_list.push_back( h_SL_fatjet_AmbitiousSelection_Cate6j4b_TwoFatJets = new TH1D( histname  , histname ,  10 , 0 , 10 ) );
    sprintf( histname , "CutFlow With FatJetInfo"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);



    
    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneLooseBBtagFatJet" , "Pt" );;
    histogram_list.push_back(   h_SL_fatjet_Pt_TheFatJet_Cate6j4b_OneLooseBBtagFatJet   = new TH1D( histname  , histname ,  10 , 200 , 500 ) );
    sprintf( histname , "LooseBB-AK8Jet %s" , "PT[GeV]"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneLooseBBtagFatJet" , "Eta" );;
    histogram_list.push_back(  h_SL_fatjet_Eta_TheFatJet_Cate6j4b_OneLooseBBtagFatJet      = new TH1D( histname  , histname ,  12 , -2.4 , 2.4) );
    sprintf( histname , "LooseBB-AK8Jet %s" , "#eta"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneLooseBBtagFatJet" , "AbsEta" );;
    histogram_list.push_back(  h_SL_fatjet_AbsEta_TheFatJet_Cate6j4b_OneLooseBBtagFatJet      = new TH1D( histname  , histname ,  6 , 0 , 2.4) );
    sprintf( histname , "LooseBB-AK8Jet %s" , "|#eta|"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneLooseBBtagFatJet" , "SDmass" );;
    histogram_list.push_back(  h_SL_fatjet_SDMass_TheFatJet_Cate6j4b_OneLooseBBtagFatJet    = new TH1D( histname  , histname , 
													10,
													0 , 200 
													) );
    sprintf( histname , "LooseBB-AK8Jet %s" , "SDmass[GeV]"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    
    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneLooseBBtagFatJet" , "dRLep" );;
    histogram_list.push_back(  h_SL_fatjet_dRLep_TheFatJet_Cate6j4b_OneLooseBBtagFatJet = new TH1D( histname  , histname ,  11 , 0.8 , 5.2 ) );
    sprintf( histname , "LooseBB-AK8Jet %s" , "#DeltaR(lep)"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);
    
    
    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneLooseBBtagFatJet" , "dEtaLep" );;
    histogram_list.push_back(      h_SL_fatjet_dEtaLep_TheFatJet_Cate6j4b_OneLooseBBtagFatJet   = new TH1D( histname  , histname ,  10 , 0 , 4 ) );
    sprintf( histname , "LooseBB-AK8Jet %s" , "#Delta#eta(lep)"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneLooseBBtagFatJet" , "dPhiLep" );;
    histogram_list.push_back(  h_SL_fatjet_dPhiLep_TheFatJet_Cate6j4b_OneLooseBBtagFatJet  = new TH1D( histname  , histname ,  10 , 0 , M_PI ) );
    sprintf( histname , "LooseBB-AK8Jet %s" , "#Delta#phi(lep)"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);



    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneMediumBBtagFatJet" , "Pt" );;
    histogram_list.push_back(   h_SL_fatjet_Pt_TheFatJet_Cate6j4b_OneMediumBBtagFatJet   = new TH1D( histname  , histname ,  10 , 200 , 500 ) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "PT[GeV]"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneMediumBBtagFatJet" , "Eta" );;
    histogram_list.push_back(  h_SL_fatjet_Eta_TheFatJet_Cate6j4b_OneMediumBBtagFatJet      = new TH1D( histname  , histname ,  12 , -2.4 , 2.4) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "#eta"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneMediumBBtagFatJet" , "AbsEta" );;
    histogram_list.push_back(  h_SL_fatjet_AbsEta_TheFatJet_Cate6j4b_OneMediumBBtagFatJet      = new TH1D( histname  , histname ,  6 , 0 , 2.4) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "|#eta|"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);



    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneMediumBBtagFatJet" , "SDmass" );;
    histogram_list.push_back(  h_SL_fatjet_SDMass_TheFatJet_Cate6j4b_OneMediumBBtagFatJet    = new TH1D( histname  , histname , 
													 10 , 0 , 200 
													) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "SDmass[GeV]"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    
    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneMediumBBtagFatJet" , "dRLep" );;
    histogram_list.push_back(  h_SL_fatjet_dRLep_TheFatJet_Cate6j4b_OneMediumBBtagFatJet = new TH1D( histname  , histname ,  11 , 0.8 , 5.2 ) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "#DeltaR(lep)"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);
    
    
    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneMediumBBtagFatJet" , "dEtaLep" );;
    histogram_list.push_back(      h_SL_fatjet_dEtaLep_TheFatJet_Cate6j4b_OneMediumBBtagFatJet   = new TH1D( histname  , histname ,  10 , 0 , 4 ) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "#Delta#eta(lep)"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate6j4b_OneMediumBBtagFatJet" , "dPhiLep" );;
    histogram_list.push_back(  h_SL_fatjet_dPhiLep_TheFatJet_Cate6j4b_OneMediumBBtagFatJet  = new TH1D( histname  , histname ,  10 , 0 , M_PI ) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "#Delta#phi(lep)"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);



    sprintf( histname , "%s_%s" , "SL_fatjet_NMediumBBfat"  ,  "Cate5moreJ_2moreB"  );
    histogram_list.push_back( h_SL_N_MediumBBTaggedfatjet_Cate5moreJ_2moreB = new TH1D( histname  , histname ,   5, -0.5 , 4.5 ) );
    sprintf( histname , "%s, cate(%s)" , "N AK8 MedBB jets"  , ">=5j,>=2b" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet" , "Pt" );;
    histogram_list.push_back(   h_SL_fatjet_Pt_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet   = new TH1D( histname  , histname ,  10 , 200 , 500 ) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "PT[GeV]"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet" , "Eta" );;
    histogram_list.push_back(  h_SL_fatjet_Eta_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet      = new TH1D( histname  , histname ,  12 , -2.4 , 2.4) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "#eta"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet" , "AbsEta" );;
    histogram_list.push_back(  h_SL_fatjet_AbsEta_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet      = new TH1D( histname  , histname ,  6 , 0 , 2.4) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "|#eta|"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);



    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet" , "SDmass" );;
    histogram_list.push_back(  h_SL_fatjet_SDMass_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet    = new TH1D( histname  , histname , 
														  20, 0, 200 
														  ) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "SDmass[GeV]"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);

    
    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet" , "dRLep" );;
    histogram_list.push_back(  h_SL_fatjet_dRLep_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet = new TH1D( histname  , histname ,  11 , 0.8 , 5.2 ) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "#DeltaR(lep)"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);
    
    
    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet" , "dEtaLep" );;
    histogram_list.push_back(      h_SL_fatjet_dEtaLep_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet   = new TH1D( histname  , histname ,  10 , 0 , 4 ) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "#Delta#eta(lep)"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);


    sprintf( histname , "SL_fatjet_%s_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet" , "dPhiLep" );;
    histogram_list.push_back(  h_SL_fatjet_dPhiLep_TheFatJet_Cate5moreJ_2moreB_OneMediumBBtagFatJet  = new TH1D( histname  , histname ,  10 , 0 , M_PI ) );
    sprintf( histname , "MediumBB-AK8Jet %s" , "#Delta#phi(lep)"  );;
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(histname);








  }



  // Note : for the moment, I require met cut. Add comment to the x-axis label.
 
  histogram_list.push_back( h_SL_fatjet_pt_AfterAllSelection_el
			    = new TH1D("SL_fatjet_pt_AfterAllSelection_el","SL_fatjet_pt_AfterAllSelection_el",  40, 200, 1000 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet PT after all selection (el)[GeV] w/ met cut");

  histogram_list.push_back( h_SL_fatjet_eta_AfterAllSelection_el
			    = new TH1D("SL_fatjet_eta_AfterAllSelection_el","SL_fatjet_eta_AfterAllSelection_el",  40, -4, 4 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet #eta after all selection (el) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_phi_AfterAllSelection_el
			    = new TH1D("SL_fatjet_phi_AfterAllSelection_el","SL_fatjet_phi_AfterAllSelection_el",  30, -M_PI, M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet #phi after all selection (el) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_sdmass_AfterAllSelection_el
			    = new TH1D("SL_fatjet_sdmass_AfterAllSelection_el","SL_fatjet_sdmass_AfterAllSelection_el", 20, 0 , 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet SD mass after el selection[GeV] w/ met cut");

  histogram_list.push_back( h_SL_fatjet_subjet_loosetagged_AfterAllSelection_el
			    = new TH1D("SL_fatjet_subjet_loosetagged_AfterAllSelection_el","SL_fatjet_subjet_loosetagged_AfterAllSelection_el", 5, -0.5 , 4.5 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("N of loose b-tagged subjets (after el selection) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_subjet_mediumtagged_AfterAllSelection_el
			    = new TH1D("SL_fatjet_subjet_mediumtagged_AfterAllSelection_el","SL_fatjet_subjet_mediumtagged_AfterAllSelection_el", 5, -0.5 , 4.5 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("N of medium b-tagged subjets (after el selection) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_subjet_mean_btagger_AfterAllSelection_el
			    = new TH1D("SL_fatjet_subjet_mean_btagger_AfterAllSelection_el","SL_fatjet_subjet_mean_btagger_AfterAllSelection_el", 20, 0 , 1 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mean of b-tagger of subjets (after el selection) w/ met cut");

  histogram_list.push_back( h_SL_N_fatjet_AfterAllSelection_el
			    = new TH1D("SL_N_fatjet_AfterAllSelection_el","SL_N_fatjet_AfterAllSelection_el", 5, -0.5 , 4.5 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("N of AK8 jets after all selection (el) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_tau21_AfterAllSelection_el
			    = new TH1D("SL_fatjet_tau21_AfterAllSelection_el","SL_fatjet_tau21_btagger_AfterAllSelection_el", 20, 0 , 1 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("tau21 of AK8 jet(after el selection) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_tau32_AfterAllSelection_el
			    = new TH1D("SL_fatjet_tau32_AfterAllSelection_el","SL_fatjet_tau32_btagger_AfterAllSelection_el", 20, 0 , 1 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("tau32 of AK8 jet(after el selection) w/ met cut");



  histogram_list.push_back( h_SL_fatjet_pt_AfterAllSelection_mu
			    = new TH1D("SL_fatjet_pt_AfterAllSelection_mu","SL_fatjet_pt_AfterAllSelection_mu",  40, 200, 1000 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet PT after all selection (mu)[GeV] w/ met cut");

  histogram_list.push_back( h_SL_fatjet_eta_AfterAllSelection_mu
			    = new TH1D("SL_fatjet_eta_AfterAllSelection_mu","SL_fatjet_eta_AfterAllSelection_mu",  40, -4, 4 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet #eta after all selection (mu) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_phi_AfterAllSelection_mu
			    = new TH1D("SL_fatjet_phi_AfterAllSelection_mu","SL_fatjet_phi_AfterAllSelection_mu",  30, -M_PI, M_PI ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet #phi after all selection (mu) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_sdmass_AfterAllSelection_mu
			    = new TH1D("SL_fatjet_sdmass_AfterAllSelection_mu","SL_fatjet_sdmass_AfterAllSelection_mu", 20, 0 , 200 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("AK8 Jet SD mass after mu selection[GeV] w/ met cut");

  histogram_list.push_back( h_SL_fatjet_subjet_loosetagged_AfterAllSelection_mu
			    = new TH1D("SL_fatjet_subjet_loosetagged_AfterAllSelection_mu","SL_fatjet_subjet_loosetagged_AfterAllSelection_mu", 5, -0.5 , 4.5 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("N of loose b-tagged subjets (after mu selection) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_subjet_mediumtagged_AfterAllSelection_mu
			    = new TH1D("SL_fatjet_subjet_mediumtagged_AfterAllSelection_mu","SL_fatjet_subjet_mediumtagged_AfterAllSelection_mu", 5, -0.5 , 4.5 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("N of medium b-tagged subjets (after mu selection) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_subjet_mean_btagger_AfterAllSelection_mu
			    = new TH1D("SL_fatjet_subjet_mean_btagger_AfterAllSelection_mu","SL_fatjet_subjet_mean_btagger_AfterAllSelection_mu", 20, 0 , 1 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Mean of b-tagger of subjets (after mu selection) w/ met cut");

  histogram_list.push_back( h_SL_N_fatjet_AfterAllSelection_mu
			    = new TH1D("SL_N_fatjet_AfterAllSelection_mu","SL_N_fatjet_AfterAllSelection_mu", 5, -0.5 , 4.5 ) );
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("N of AK8 jets after all selection (mu) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_tau21_AfterAllSelection_mu
			    = new TH1D("SL_fatjet_tau21_AfterAllSelection_mu","SL_fatjet_tau21_btagger_AfterAllSelection_mu", 20, 0 , 1 )) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("tau21 of AK8 jet(after mu selection) w/ met cut");

  histogram_list.push_back( h_SL_fatjet_tau32_AfterAllSelection_mu
			    = new TH1D("SL_fatjet_tau32_AfterAllSelection_mu","SL_fatjet_tau32_btagger_AfterAllSelection_mu", 20, 0 , 1 )) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("tau32 of AK8 jet(after mu selection) w/ met cut");



  for( int cate = 1 ; cate < 10 ; cate++  ){
    std::string catename = _NameofCategory( cate );
    char histname [100];
    sprintf( histname , "SL_fatjet_BDT_AfterNjNbCategory%s_withMetCut_HiggsWindown" ,  catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_BDT_metcut_withInfoSDMassinHiggsWindow [cate]   = new TH1D(histname, histname , 10, -0.8 , 1.2 ) );
    sprintf( histname , "BDT, %scate (8+2 bins)" ,  catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );   


    sprintf( histname , "SL_fatjet_BDT_AfterNjNbCategory%s_withMetCut_InHiggsWindown" ,  catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_BDT_metcut_SDMass_in_HiggsWindow [cate]   = new TH1D(histname, histname , 8, -0.8  ,  0.8 ) );
    sprintf( histname , "BDT, %scate, SDmass in Higgs window" ,  catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );   

    sprintf( histname , "SL_fatjet_BDT_AfterNjNbCategory%s_withMetCut_OutHiggsWindown" ,  catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_BDT_metcut_SDMass_out_HiggsWindow [cate]   = new TH1D(histname, histname , 8, -0.8 ,  0.8 ) );
    sprintf( histname , "BDT, %scate, SDmass out Higgs window" ,  catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );   


    sprintf( histname , "SL_fatjet_isSDMassInHiggsWidnwo_AfterNjNbCategory%s_withMetCut_HiggsWindown" ,  catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_IsSDMassInHiggsWindow [cate]   = new TH1D(histname, histname , 2, -0.5 , 1.5 ) );
    sprintf( histname , "Is SD mass in Higgs window, cate%s" ,  catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    


    sprintf( histname , "SL_fatjet_BDT_category%s__ExactOneBBFatJetPass" ,  catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_BDT_ExactOneBBFatJetPass  [cate]   = new TH1D(histname, histname , 8, -0.8 , 0.8 ) );
    sprintf( histname , "BDT, cate%s with one BB-ak8" ,  catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    

    sprintf( histname , "SL_fatjet_BDT_category%s__ExactOneBBFatJetFail" ,  catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_BDT_ExactOneBBFatjetFail  [cate]   = new TH1D(histname, histname , 8, -0.8 , 0.8 ) );
    sprintf( histname , "BDT, cate%s not exact one BB-ak8" ,  catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    


    sprintf( histname , "SL_fatjet_BDT_category%s__ExactOneBBFatJetPass_SDMassPass" ,  catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_BDT_ExactOneBBFatJetPass_SDMassPass  [cate]   = new TH1D(histname, histname , 2, -0.8 , 0.8 ) ); // Only 2-bins on purpose, because of statistics.
    sprintf( histname , "BDT, cate%s with one BB-ak8" ,  catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    

    sprintf( histname , "SL_fatjet_BDT_category%s__ExactOneBBFatJetPass_SDMassFail" ,  catename.c_str() );
    histogram_list.push_back( h_SL_fatjet_BDT_ExactOneBBFatJetPass_SDMassFail  [cate]   = new TH1D(histname, histname , 8, -0.8 , 0.8 ) );
    sprintf( histname , "BDT, cate%s with one BB-ak8" ,  catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    


  }


  
#ifdef ENABLE_COMMON_CLASSFIER

  for( int cate = 1 ; cate < 10 ; cate++  ){

    std::string catename = _NameofCategory( cate );
    char histname [100];
    sprintf( histname , "SL_BDT_AfterNjNbCategory%s" ,  catename.c_str() );
    histogram_list.push_back( h_SL_BDT[cate]   = new TH1D(histname, histname , 8, -0.8 , 0.8 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    

    sprintf( histname , "SL_BDT_AfterNjNbCategory%s_withMetCut" ,  catename.c_str() );
    histogram_list.push_back( h_SL_BDT_metcut[cate]   = new TH1D(histname, histname , 8, -0.8 , 0.8 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    


    sprintf( histname , "SL_BDT_AfterNjNbCategory%s_withMetCut_el" ,  catename.c_str() );
    histogram_list.push_back( h_SL_BDT_metcut_el[cate]   = new TH1D(histname, histname , 8, -0.8 , 0.8 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    

    sprintf( histname , "SL_BDT_AfterNjNbCategory%s_withMetCut_mu" ,  catename.c_str() );
    histogram_list.push_back( h_SL_BDT_metcut_mu[cate]   = new TH1D(histname, histname , 8, -0.8 , 0.8 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    


    sprintf( histname , "SL_BDT_AfterNjNbCategory%s_noMetCut_el" ,  catename.c_str() );
    histogram_list.push_back( h_SL_BDT_nometcut_el[cate]   = new TH1D(histname, histname , 8, -0.8 , 0.8 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    

    sprintf( histname , "SL_BDT_AfterNjNbCategory%s_noMetCut_mu" ,  catename.c_str() );
    histogram_list.push_back( h_SL_BDT_nometcut_mu[cate]   = new TH1D(histname, histname , 8, -0.8 , 0.8 ) );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    


    sprintf( histname , "SL_MuonPT_AfterNjNbCategory%s_withMetCut" ,  catename.c_str() );
    histogram_list.push_back( h_SL_MuonPtAfterCate_MetCut[ cate ] = new TH1D( histname, histname, 9, 20, 200) ) ;
    sprintf( histname , "Cate%s MuonPT(GeV) withMetCut" ,  catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    

    sprintf( histname , "SL_ElPT_AfterNjNbCategory%s_withMetCut" ,  catename.c_str() );
    if( catename == "4j4b" || catename == "5j4b" || catename == "6j4b" ){ 
      histogram_list.push_back( h_SL_ElPtAfterCate_MetCut[ cate ] = new TH1D( histname, histname, 6, 30, 210) ) ;
    }else if( catename == "4j3b" || catename == "5j3b" || catename == "6j3b" ){
      histogram_list.push_back( h_SL_ElPtAfterCate_MetCut[ cate ] = new TH1D( histname, histname, 10, 30, 330) ) ;
    }else{
      histogram_list.push_back( h_SL_ElPtAfterCate_MetCut[ cate ] = new TH1D( histname, histname, 18, 30, 300) ) ;
    }
    sprintf( histname , "Cate%s ElonPT(GeV) withMetCut" ,  catename.c_str() );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );    
    // todo

  }

#endif

  h2_SL_NPV_NJet_AfterLeptonEventSelection_mu = new TH2D("h2_SL_NPV_NJet_AfterLeptonEventSelection_mu", "h2_SL_NPV_NJet_AfterLeptonEventSelection_mu", 
							 10 , 0, 50, 11, -0.5, 10.5 );
  h2_SL_NPV_NJet_AfterLeptonEventSelection_mu -> GetXaxis() ->SetTitle("NPV");
  h2_SL_NPV_NJet_AfterLeptonEventSelection_mu -> GetYaxis() ->SetTitle("NJet");

  h2_FakeEval_NEntries_El_nJnB_withFakeWgt = new TH2D( "FakeEval_NEntries_El_nJnB_withFakeWgt","FakeEval_NEntries_El_nJnB_withFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 
  h2_FakeEval_NEntries_El_nJnB_noFakeWgt   = new TH2D( "FakeEval_NEntries_El_nJnB_noFakeWgt","FakeEval_NEntries_El_nJnB_noFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 
  h2_FakeEval_NEntries_Mu_nJnB_withFakeWgt = new TH2D( "FakeEval_NEntries_Mu_nJnB_withFakeWgt","FakeEval_NEntries_Mu_nJnB_withFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 
  h2_FakeEval_NEntries_Mu_nJnB_noFakeWgt   = new TH2D( "FakeEval_NEntries_Mu_nJnB_noFakeWgt","FakeEval_NEntries_Mu_nJnB_noFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 

  h2_FakeEval_NEntries_LowMetSide_El_nJnB_withFakeWgt = new TH2D( "FakeEval_NEntries_LowMetSide_El_nJnB_withFakeWgt","FakeEval_NEntries_LowMetSide_El_nJnB_withFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 
  h2_FakeEval_NEntries_LowMetSide_El_nJnB_noFakeWgt   = new TH2D( "FakeEval_NEntries_LowMetSide_El_nJnB_noFakeWgt"  ,"FakeEval_NEntries_LowMetSide_El_nJnB_noFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 
  h2_FakeEval_NEntries_LowMetSide_Mu_nJnB_withFakeWgt = new TH2D( "FakeEval_NEntries_LowMetSide_Mu_nJnB_withFakeWgt","FakeEval_NEntries_LowMetSide_Mu_nJnB_withFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 
  h2_FakeEval_NEntries_LowMetSide_Mu_nJnB_noFakeWgt   = new TH2D( "FakeEval_NEntries_LowMetSide_Mu_nJnB_noFakeWgt"  ,"FakeEval_NEntries_LowMetSide_Mu_nJnB_noFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 
						      
  h2_FakeEval_NEntries_HighMetSide_El_nJnB_withFakeWgt = new TH2D( "FakeEval_NEntries_HighMetSide_El_nJnB_withFakeWgt","FakeEval_NEntries_HighMetSideEl_nJnB_withFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 
  h2_FakeEval_NEntries_HighMetSide_El_nJnB_noFakeWgt   = new TH2D( "FakeEval_NEntries_HighMetSide_El_nJnB_noFakeWgt"  ,"FakeEval_NEntries_HighMetSideEl_nJnB_noFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 
  h2_FakeEval_NEntries_HighMetSide_Mu_nJnB_withFakeWgt = new TH2D( "FakeEval_NEntries_HighMetSide_Mu_nJnB_withFakeWgt","FakeEval_NEntries_HighMetSideMu_nJnB_withFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 
  h2_FakeEval_NEntries_HighMetSide_Mu_nJnB_noFakeWgt   = new TH2D( "FakeEval_NEntries_HighMetSide_Mu_nJnB_noFakeWgt"  ,"FakeEval_NEntries_HighMetSideMu_nJnB_noFakeWgt", 7, 0 , 7 , 5 , 0, 5 ) ; 

  h2_FakeEval_NonisoEl_njetNbin = new TH2D( "FakeEval_NonisoEl_njetNbin","FakeEval_NonisoEl_njetNbin", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_NonisoEl_njetNbinCutOff = new TH2D( "FakeEval_NonisoEl_njetNbin_cutoff","FakeEval_NonisoEl_njetNbin_cutoff", 7, 0 ,7 , 5, 0, 5 );
  h2_FakeEval_NonisoEl_njetNbin_LowMet = new TH2D( "FakeEval_NonisoEl_njetNbin_met30", "FakeEval_NonisoEl_njetNbin_met30", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_NonisoEl_njetNbinCutOff_LowMet = new TH2D( "FakeEval_NonisoEl_njetNbin_met30_cutoff","FakeEval_NonisoEl_njetNbin_met30_cutoff" , 7, 0 ,7 , 5, 0, 5 );
  h2_FakeEval_NonisoEl_njetNbin_notLowMet = new TH2D( "FakeEval_NonisoEl_njetNbin_notmet30", "FakeEval_NonisoEl_njetNbin_notmet30", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_NonisoEl_njetNbinCutOff_notLowMet = new TH2D( "FakeEval_NonisoEl_njetNbin_notmet30_cutoff","FakeEval_NonisoEl_njetNbin_notmet30_cutoff" , 7, 0 ,7 , 5, 0, 5 );

  h2_FakeEval_IsoEl_njetNbin = new TH2D( "FakeEval_IsoEl_njetNbin","FakeEval_IsoEl_njetNbin", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_IsoEl_njetNbinCutOff = new TH2D( "FakeEval_IsoEl_njetNbin_cutoff","FakeEval_IsoEl_njetNbin_cutoff", 7, 0 ,7 , 5, 0, 5 );
  h2_FakeEval_IsoEl_njetNbin_LowMet = new TH2D( "FakeEval_IsoEl_njetNbin_met30", "FakeEval_IsoEl_njetNbin_met30", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_IsoEl_njetNbinCutOff_LowMet = new TH2D( "FakeEval_IsoEl_njetNbin_met30_cutoff","FakeEval_IsoEl_njetNbin_met30_cutoff" , 7, 0 ,7 , 5, 0, 5 );
  h2_FakeEval_IsoEl_njetNbin_notLowMet = new TH2D( "FakeEval_IsoEl_njetNbin_notmet30", "FakeEval_IsoEl_njetNbin_notmet30", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_IsoEl_njetNbinCutOff_notLowMet = new TH2D( "FakeEval_IsoEl_njetNbin_notmet30_cutoff","FakeEval_IsoEl_njetNbin_notmet30_cutoff" , 7, 0 ,7 , 5, 0, 5 );

  h2_FakeEval_NonisoMu_njetNbin = new TH2D( "FakeEval_NonisoMu_njetNbin","FakeEval_NonisoMu_njetNbin", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_NonisoMu_njetNbinCutOff = new TH2D( "FakeEval_NonisoMu_njetNbin_cutoff","FakeEval_NonisoMu_njetNbin_cutoff", 7, 0 ,7 , 5, 0, 5 );
  h2_FakeEval_NonisoMu_njetNbin_LowMet = new TH2D( "FakeEval_NonisoMu_njetNbin_met30", "FakeEval_NonisoMu_njetNbin_met30", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_NonisoMu_njetNbinCutOff_LowMet = new TH2D( "FakeEval_NonisoMu_njetNbin_met30_cutoff","FakeEval_NonisoMu_njetNbin_met30_cutoff" , 7, 0 ,7 , 5, 0, 5 );
  h2_FakeEval_NonisoMu_njetNbin_notLowMet = new TH2D( "FakeEval_NonisoMu_njetNbin_notmet30", "FakeEval_NonisoMu_njetNbin_notmet30", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_NonisoMu_njetNbinCutOff_notLowMet = new TH2D( "FakeEval_NonisoMu_njetNbin_notmet30_cutoff","FakeEval_NonisoMu_njetNbin_notmet30_cutoff" , 7, 0 ,7 , 5, 0, 5 );

  h2_FakeEval_IsoMu_njetNbin = new TH2D( "FakeEval_IsoMu_njetNbin","FakeEval_IsoMu_njetNbin", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_IsoMu_njetNbinCutOff = new TH2D( "FakeEval_IsoMu_njetNbin_cutoff","FakeEval_IsoMu_njetNbin_cutoff", 7, 0 ,7 , 5, 0, 5 );
  h2_FakeEval_IsoMu_njetNbin_LowMet = new TH2D( "FakeEval_IsoMu_njetNbin_met30", "FakeEval_IsoMu_njetNbin_met30", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_IsoMu_njetNbinCutOff_LowMet = new TH2D( "FakeEval_IsoMu_njetNbin_met30_cutoff","FakeEval_IsoMu_njetNbin_met30_cutoff" , 7, 0 ,7 , 5, 0, 5 );
  h2_FakeEval_IsoMu_njetNbin_notLowMet = new TH2D( "FakeEval_IsoMu_njetNbin_notmet30", "FakeEval_IsoMu_njetNbin_notmet30", 10, 0 , 10 , 6, 0, 6 );
  h2_FakeEval_IsoMu_njetNbinCutOff_notLowMet = new TH2D( "FakeEval_IsoMu_njetNbin_notmet30_cutoff","FakeEval_IsoMu_njetNbin_notmet30_cutoff" , 7, 0 ,7 , 5, 0, 5 );

  histogram_list.push_back( h_FakeEval_nJet_El_onejetcut_noFakeWgt     = new TH1D( "FakeEval_nJet_el_oneJetCut_noFakeWgt","FakeEval_nJet_el_oneJetCut_noFakeWgt", 10 , 0 , 10 ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Jet multiplicity");

  histogram_list.push_back( h_FakeEval_nBJet_El_onejetcut_noFakeWgt    = new TH1D( "FakeEval_nBJet_El_onejetcut_noFakeWgt", "FakeEval_nBJet_El_onejetcut_noFakeWgt", 6 , 0 , 6 ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Btagged jet multiplicity");

  histogram_list.push_back( h_FakeEval_nJet_El_onejetcut_withFakeWgt   = new TH1D( "FakeEval_nJet_El_onejetcut_withFakeWgt", "FakeEval_nJet_El_onejetcut_withFakeWgt", 10 , 0 , 10 ) ) ; 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Jet multiplicity");

  histogram_list.push_back( h_FakeEval_nBJet_El_onejetcut_withFakeWgt  = new TH1D( "FakeEval_nBJet_El_onejetcut_withFakeWgt", "FakeEval_nBJet_El_onejetcut_withFakeWgt", 6 , 0 , 6 ) ) ; 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Btagged jet multiplicity");

  histogram_list.push_back( h_FakeEval_nJet_Mu_onejetcut_noFakeWgt     = new TH1D( "FakeEval_nJet_Mu_onejetcut_noFakeWgt", "FakeEval_nJet_Mu_onejetcut_noFakeWgt", 10 , 0 , 10 ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Jet multiplicity");

  histogram_list.push_back( h_FakeEval_nBJet_Mu_onejetcut_noFakeWgt    = new TH1D( "FakeEval_nBJet_Mu_onejetcut_noFakeWgt", "FakeEval_nBJet_Mu_onejetcut_noFakeWgt", 6 , 0 , 6 ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Btagged jet multiplicity");

  histogram_list.push_back( h_FakeEval_nJet_Mu_onejetcut_withFakeWgt   = new TH1D( "FakeEval_nJet_Mu_onejetcut_withFakeWgt", "FakeEval_nJet_Mu_onejetcut_withFakeWgt" , 10 , 0 , 10 ) ) ; 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Jet multiplicity");

  histogram_list.push_back( h_FakeEval_nBJet_Mu_onejetcut_withFakeWgt  = new TH1D( "FakeEval_nBJet_Mu_onejetcut_withFakeWgt","FakeEval_nBJet_Mu_onejetcut_withFakeWgt", 6 , 0 , 6 ) ) ; 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Btagged jet multiplicity");

  for (int j = 1 ; j <= 6 ; j ++){
    for (int b = 0 ; b <= j ; b ++){
      char name [100];
      sprintf( name , "FakeEval_Met_NonIsoEl_%dj%db" , j , b );
      histogram_list.push_back( h_FakeEval_Met_NonIsoEl_njnb[ j ][ b ] = new TH1D( name,name, 20 , 0 , 100 ) );
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV]");

      sprintf( name , "FakeEval_Met_IsoEl_%dj%db" , j , b );
      histogram_list.push_back( h_FakeEval_Met_IsoEl_njnb[ j ][ b ] = new TH1D( name,name, 20 , 0 , 100 ) );
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV]");

      sprintf( name , "FakeEval_Met_NonIsoMu_%dj%db" , j , b );
      histogram_list.push_back( h_FakeEval_Met_NonIsoMu_njnb[ j ][ b ] = new TH1D( name,name, 20 , 0 , 100 ) );
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV]");

      sprintf( name , "FakeEval_Met_IsoMu_%dj%db" , j , b );
      histogram_list.push_back( h_FakeEval_Met_IsoMu_njnb[ j ][ b ] = new TH1D( name,name, 20 , 0 , 100 ) );
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV]");


      sprintf( name , "FakeEval_Met_El_Cate%dj%db_NoFakeWgt" , j , b );
      histogram_list.push_back( h_FakeEval_Met_El_njnb_noFakeWgt   [ j ][ b ] = new TH1D( name,name, 20 , 0 , 100 ) ) ;
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV]");

      sprintf( name , "FakeEval_NentriesInLowMet_El_Cate%dj%db_NoFakeWgt" , j , b );
      histogram_list.push_back( h_FakeEval_NEntries_El_LowMet_noFakeWgt  [ j ][ b ] = new TH1D( name,name, 1, -10 , 10 ) ) ;
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Entries");

      sprintf( name , "FakeEval_Met_El_Cate%dj%db_WithFakeWgt" , j , b );
      if( b==4 ){ 
	histogram_list.push_back( h_FakeEval_Met_El_njnb_withFakeWgt       [ j ][ b ] = new TH1D( name,name, 5 , 0 , 200 ) ) ;
      }else if( b==3 ){ 
	histogram_list.push_back( h_FakeEval_Met_El_njnb_withFakeWgt       [ j ][ b ] = new TH1D( name,name, 10 , 0 , 200 ) ) ;
      }else{
	histogram_list.push_back( h_FakeEval_Met_El_njnb_withFakeWgt       [ j ][ b ] = new TH1D( name,name, 20 , 0 , 200 ) ) ;
      }
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV]");

      sprintf( name , "FakeEval_NentriesInLowMet_El_Cate%dj%db_WithFakeWgt" , j , b );
      histogram_list.push_back( h_FakeEval_NEntries_El_LowMet_withFakeWgt[ j ][ b ] = new TH1D( name,name, 1, -10 , 10 ) );
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Entries");      

      sprintf( name , "FakeEval_Met_Mu_Cate%dj%db_NoFakeWgt" , j , b );
      histogram_list.push_back( h_FakeEval_Met_Mu_njnb_noFakeWgt         [ j ][ b ] = new TH1D( name,name, 20 , 0 , 100 ) ) ;
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV]");

      sprintf( name , "FakeEval_NentriesInLowMet_Mu_Cate%dj%db_NoFakeWgt" , j , b );
      histogram_list.push_back( h_FakeEval_NEntries_Mu_LowMet_noFakeWgt  [ j ][ b ] = new TH1D( name,name, 1, -10 , 10 ) ) ;
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Entries");

      {
	sprintf( name , "FakeEval_NentriesvsNPV_InLowMet_Mu_Cate%dj%db_NoFakeWgt" , j , b );
	const int n_bin = 6  ;
	const double binning [ n_bin + 1 ] = { 0, 10, 20, 25, 30, 40, 60 };
	histogram_list . push_back( h_FakeEval_NEntriesVsNPV_Mu_LowMet_noFakeWgt [j][b] = new TH1D( name, name , n_bin, binning ) ) ;
	histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NPV");
      }

      sprintf( name , "FakeEval_Met_Mu_Cate%dj%db_WithFakeWgt" , j , b );
      if( b == 3 || b == 4 ){
	histogram_list.push_back( h_FakeEval_Met_Mu_njnb_withFakeWgt       [ j ][ b ] = new TH1D( name,name, 10 , 0 , 200 ) ) ;
      }else{
	histogram_list.push_back( h_FakeEval_Met_Mu_njnb_withFakeWgt       [ j ][ b ] = new TH1D( name,name, 20 , 0 , 200 ) ) ;
      }
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("MET[GeV]");

      sprintf( name , "FakeEval_NentriesInLowMet_Mu_Cate%dj%db_WithFakeWgt" , j , b );
      histogram_list.push_back( h_FakeEval_NEntries_Mu_LowMet_withFakeWgt[ j ][ b ] = new TH1D( name,name, 1, -10 , 10 ) );
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Entries");
    

      sprintf( name , "FakeEval_dR_Jet_and_Mu_onejetcut_Cate%dj%db_noFakeWgt" , j , b );
      histogram_list.push_back( h_FakeEval_dR_Jet_and_Mu_onejetcut_noFakeWgt [j][b] = new TH1D( name,name, 10 , 0.4 , 2.4 ) ) ;
      sprintf( name , "minimum dR(Jet,Mu), %dj%db" , j , b );
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

      sprintf( name , "FakeEval_dR_BJet_and_Mu_onejetcut_Cate%dj%db_noFakeWgt" , j , b );
      histogram_list.push_back( h_FakeEval_dR_Btag_and_Mu_onejetcut_noFakeWgt[j][b] = new TH1D( name,name, 10 , 0.4 , 2.4 ) ) ;
      sprintf( name , "minimum dR(BJet,Mu), %dj%db" , j , b );
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


      sprintf( name , "FakeEval_dR_Jet_and_El_onejetcut_Cate%dj%db_noFakeWgt" , j , b );
      histogram_list.push_back( h_FakeEval_dR_Jet_and_El_onejetcut_noFakeWgt [j][b] = new TH1D( name,name, 10 , 0.4 , 2.4 ) ) ;
      sprintf( name , "minimum dR(Jet,El), %dj%db" , j , b );
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);


      sprintf( name , "FakeEval_dR_BJet_and_El_onejetcut_Cate%dj%db_noFakeWgt" , j , b );
      histogram_list.push_back( h_FakeEval_dR_Btag_and_El_onejetcut_noFakeWgt[j][b] = new TH1D( name,name, 10 , 0.4 , 2.4 ) ) ;
      sprintf( name , "minimum dR(BJet,El), %dj%db" , j , b );
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle(name);

    }
  }

  histogram_list.push_back( h_FakeEval_dR_Jet_and_Mu_onejetcut_noFakeWgt_1j0b  = new TH1D( "FakeEval_dR_BJet_and_Mu_1j0b_noFakeWgt", "FakeEval_dR_BJet_and_Mu_1j0b_noFakeWgt", 30 , 0 , 6 ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("dR(jet,mu), 1j0b");

  histogram_list.push_back( h_FakeEval_dR_Jet_and_Mu_onejetcut_noFakeWgt_1j1b  = new TH1D( "FakeEval_dR_BJet_and_Mu_1j1b_noFakeWgt", "FakeEval_dR_BJet_and_Mu_1j1b_noFakeWgt", 30 , 0 , 6 ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("dR(jet,mu), 1j1b");

  histogram_list.push_back( h_FakeEval_dR_Jet_and_El_onejetcut_noFakeWgt_1j0b  = new TH1D( "FakeEval_dR_BJet_and_El_1j0b_noFakeWgt", "FakeEval_dR_BJet_and_El_1j0b_noFakeWgt", 30 , 0 , 6 ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("dR(jet,el), 1j0b");
  
  histogram_list.push_back( h_FakeEval_dR_Jet_and_El_onejetcut_noFakeWgt_1j1b  = new TH1D( "FakeEval_dR_BJet_and_El_1j1b_noFakeWgt", "FakeEval_dR_BJet_and_El_1j1b_noFakeWgt", 30 , 0 , 6 ) ) ;
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("dR(jet,el), 1j1b");



  histogram_list.push_back( h_ee_SameSign_SanityCheck_NlooseLepMinusTight = new TH1D( "ee_samesign_SanityCheck_NlooseLepMinusTight" , "ee_samesign_SanityCheck_NlooseLepMinusTight ", 10, -5 , 5 ) ) ; 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NLooseLep-NTightLep");

  histogram_list.push_back( h_ee_SameSign_SanityCheck_NlooseLepMinusNonIsoTight = new TH1D( "ee_samesign_SanityCheck_NlooseLepMinusNonIsoTight" , "ee_samesign_SanityCheck_NlooseLepMinusNonIsoTight ", 10, -5 , 5 ) ) ; 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("NLooseLep-NNonIsoTightLep");




  histogram_list.push_back( h_ee_SameSign_BtagjetMultiplicity  = new TH1D( "ee_samesign_btagjetmultipilicity", "ee_samesign_btagjetmultipilicity", 5, 0 , 5 ) ) ; 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("ee_samesign btag multiplicity");
  histogram_list[histogram_list.size()-1]->GetYaxis()->SetTitle("Entries(no Met cut)");


  histogram_list.push_back( h_ee_SameSign_BtagjetMultiplicity_MllCut  = new TH1D( "ee_samesign_btagjetmultipilicity_MllCut", "ee_samesign_btagjetmultipilicity_MllCut", 5, 0 , 5 ) ) ; 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("ee_samesign btag multiplicity");
  histogram_list[histogram_list.size()-1]->GetYaxis()->SetTitle("Entries(after Met cut)");


  histogram_list.push_back( h_ee_SameSign_mll                  = new TH1D( "ee_samesign_mll", "ee_samesign_mll", 20, 0 , 200 ) ) ; 
  histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("ee_samesign Mee[GeV]");
  
  char name [100];
  for( int i = 0 ; i < 2 ; i ++  ){
    for( int jet = 0 ; jet <= 4 ; jet ++ ){
      sprintf( name, "ee_samesign__NEvent__%sMllCut_%dBJetCate" , i == 0 ? "In" : "Out", jet );
      histogram_list.push_back(h_ee_SameSign_ZmassWindow  [i][jet] =  new TH1D(name, name, 1, 0 , 1 ) ) ;
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Entry");
      
      sprintf( name, "ee_samesign__ElPt__%sMllCut_%dBJetCate" , i == 0 ? "In" : "Out", jet );
      histogram_list.push_back(h_ee_SameSign_ElectronPT   [i][jet] =  new TH1D(name, name, 20, 0 , 200 ) ) ;
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("El Pt[GeV]");
      
      sprintf( name, "ee_samesign__LeadingLetPt__%sMllCut_%dBJetCate" , i == 0 ? "In" : "Out", jet );
      histogram_list.push_back(h_ee_SameSign_LeadingJetPT [i][jet] =  new TH1D(name, name, 20, 0 , 200 ) ) ;
      histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle("Leading Jet Pt[GeV]");
    }
  }


#endif
  // -> this is end-if of the "else" of "TRUTH_INFO_STUDY"


  {

    char histname[100];

    sprintf( histname , "SLFatjetStudy_BDTScore_6j4b");
    histogram_list.push_back( h_SLFatjetStudy_BDTScore_6j4b = new TH1D( histname, histname, 8, -0.8, 0.8) ) ;
    sprintf( histname , "BDT, 6j4b, w/ 20GeV Metcut" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );


    sprintf( histname , "SLFatjetStudy_BDTScore_6j4b_nohiggsak8candidate");
    histogram_list.push_back( h_SLFatjetStudy_BDTScore_6j4b_NoHiggsAk8Candidate = new TH1D( histname, histname, 8, -0.8, 0.8) ) ;
    sprintf( histname , "BDT, 6j4b/Metcut/NoBoostedAk8" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );

    sprintf( histname , "SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate");
    histogram_list.push_back( h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate = new TH1D( histname, histname, 8, -0.8, 0.8) ) ;
    sprintf( histname , "BDT, 6j4b/Metcut/WithBoostedAk8" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );

    sprintf( histname , "SLFatjetStudy_ak8BDTScore_6j4b");
    histogram_list.push_back( h_SLFatjetStudy_ak8BDTScore_6j4b = new TH1D( histname, histname, 16, -0.4, 0.4) ) ;
    sprintf( histname , "AK8BDT, 6j4b/Metcut" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );

    sprintf( histname , "SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_highak8bdt");
    histogram_list.push_back( h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_highak8bdt  = new TH1D( histname, histname, 4, -0.8, 0.8) ) ;
    sprintf( histname , "BDT, 6j4b/Metcut/HighAK8BDT" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );

    sprintf( histname , "SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_lowak8bdt");
    histogram_list.push_back( h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_lowak8bdt  = new TH1D( histname, histname, 8, -0.8, 0.8) ) ;
    sprintf( histname , "BDT, 6j4b/Metcut/LowAK8BDT" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );


    sprintf( histname , "SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_highak8bdt_0p1");
    histogram_list.push_back( h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_high2ak8bdt  = new TH1D( histname, histname, 4, -0.8, 0.8) ) ;
    sprintf( histname , "BDT, 6j4b/Metcut/HighAK8BDT" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );

    sprintf( histname , "SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_lowak8bdt_0p1");
    histogram_list.push_back( h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_low2ak8bdt  = new TH1D( histname, histname, 8, -0.8, 0.8) ) ;
    sprintf( histname , "BDT, 6j4b/Metcut/LowAK8BDT" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );


    sprintf( histname , "SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_highak8bdt_0p0");
    histogram_list.push_back( h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_high0ak8bdt  = new TH1D( histname, histname, 4, -0.8, 0.8) ) ;
    sprintf( histname , "BDT, 6j4b/Metcut/HighAK8BDT" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );

    sprintf( histname , "SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_lowak8bdt_0p0");
    histogram_list.push_back( h_SLFatjetStudy_BDTScore_6j4b_WithHiggsCandidate_low0ak8bdt  = new TH1D( histname, histname, 8, -0.8, 0.8) ) ;
    sprintf( histname , "BDT, 6j4b/Metcut/LowAK8BDT" );
    histogram_list[histogram_list.size()-1]->GetXaxis()->SetTitle( histname );


  }




}

void analyzer::FillHistogram( TH1D * h , double val , double weight ){
  
  TAxis * axis = h->GetXaxis() ;

  const double max = axis->GetBinUpEdge ( axis->GetNbins() );
  val = val < max ? val : axis->GetBinCenter( axis -> GetNbins() )  ;

  const double min = axis->GetBinLowEdge ( 1 );
  val = val > min ? val : axis->GetBinCenter( 1 );

  h -> Fill( val , ( isMC || b_FakeEstimation ? weight : 1.0 ) ) ;

}

void analyzer::FillHistogram2D( TH2D * h , double val_x , double val_y,  double weight ){
  
  {
    TAxis * axis = h->GetXaxis() ;
    const double max = axis->GetBinUpEdge ( axis->GetNbins() );
    val_x = val_x < max ? val_x : axis->GetBinCenter( axis -> GetNbins() )  ;
    const double min = axis->GetBinLowEdge ( 1 );
    val_x = val_x > min ? val_x : axis->GetBinCenter( 1 );
  }
  {
    TAxis * axis = h->GetYaxis() ;
    const double max = axis->GetBinUpEdge ( axis->GetNbins() );
    val_y = val_y < max ? val_y : axis->GetBinCenter( axis -> GetNbins() )  ;
    const double min = axis->GetBinLowEdge ( 1 );
    val_y = val_y > min ? val_y : axis->GetBinCenter( 1 );
  }
  h -> Fill( val_x ,  val_y, ( isMC || b_FakeEstimation ? weight : 1.0 ) ) ;

}




void analyzer::postProcess(){

  main_directory -> cd();

  // writing histograms, etc...

  //  outtree ->Write();

  for( std::vector<TH1D *>::iterator h = histogram_list.begin();
       h != histogram_list.end();
       h++ ){
    (*h)->Write();
  }



#ifdef TRUTH_INFO_STUDY

#else

  h2_SL_NPV_NJet_AfterLeptonEventSelection_mu -> Write();

  h2_FakeEval_NEntries_El_nJnB_withFakeWgt -> Write() ;
  h2_FakeEval_NEntries_El_nJnB_noFakeWgt   -> Write() ;
  h2_FakeEval_NEntries_Mu_nJnB_withFakeWgt -> Write() ;
  h2_FakeEval_NEntries_Mu_nJnB_noFakeWgt   -> Write() ;

  h2_FakeEval_NEntries_LowMetSide_El_nJnB_withFakeWgt -> Write() ; 
  h2_FakeEval_NEntries_LowMetSide_El_nJnB_noFakeWgt   -> Write() ; 
  h2_FakeEval_NEntries_LowMetSide_Mu_nJnB_withFakeWgt -> Write() ; 
  h2_FakeEval_NEntries_LowMetSide_Mu_nJnB_noFakeWgt   -> Write() ; 

  h2_FakeEval_NEntries_HighMetSide_El_nJnB_withFakeWgt-> Write() ; 
  h2_FakeEval_NEntries_HighMetSide_El_nJnB_noFakeWgt  -> Write() ; 
  h2_FakeEval_NEntries_HighMetSide_Mu_nJnB_withFakeWgt-> Write() ; 
  h2_FakeEval_NEntries_HighMetSide_Mu_nJnB_noFakeWgt  -> Write() ; 

  h2_FakeEval_NonisoEl_njetNbin ->Write();
  h2_FakeEval_NonisoEl_njetNbinCutOff ->Write();
  h2_FakeEval_NonisoEl_njetNbin_LowMet ->Write();
  h2_FakeEval_NonisoEl_njetNbinCutOff_LowMet ->Write();
  h2_FakeEval_NonisoEl_njetNbin_notLowMet ->Write();
  h2_FakeEval_NonisoEl_njetNbinCutOff_notLowMet ->Write();

  h2_FakeEval_IsoEl_njetNbin ->Write();
  h2_FakeEval_IsoEl_njetNbinCutOff ->Write();
  h2_FakeEval_IsoEl_njetNbin_LowMet ->Write();
  h2_FakeEval_IsoEl_njetNbinCutOff_LowMet ->Write();
  h2_FakeEval_IsoEl_njetNbin_notLowMet ->Write();
  h2_FakeEval_IsoEl_njetNbinCutOff_notLowMet ->Write();

  h2_FakeEval_NonisoMu_njetNbin ->Write();
  h2_FakeEval_NonisoMu_njetNbinCutOff ->Write();
  h2_FakeEval_NonisoMu_njetNbin_LowMet ->Write();
  h2_FakeEval_NonisoMu_njetNbinCutOff_LowMet ->Write();
  h2_FakeEval_NonisoMu_njetNbin_notLowMet ->Write();
  h2_FakeEval_NonisoMu_njetNbinCutOff_notLowMet ->Write();

  h2_FakeEval_IsoMu_njetNbin ->Write();
  h2_FakeEval_IsoMu_njetNbinCutOff ->Write();
  h2_FakeEval_IsoMu_njetNbin_LowMet ->Write();
  h2_FakeEval_IsoMu_njetNbinCutOff_LowMet ->Write();
  h2_FakeEval_IsoMu_njetNbin_notLowMet ->Write();
  h2_FakeEval_IsoMu_njetNbinCutOff_notLowMet ->Write();

#endif

  f_out -> Close();

}


void analyzer::SetDataModeOn(){
  isMC = false;
};

void analyzer::SetIsMuonStream( bool _isMuonStream  ){

  isMuonStream = _isMuonStream ;
  SetDataModeOn();

}


long analyzer::_EventCateBasedOnNjetNBtagJet( long nJ , long nB ){


  if ( nJ <  4 ) return 0 ; 
  if ( nB <  2 ) return 0 ; 
  if ( nJ == 4 && nB == 2 ) return 1 ; 
  if ( nJ == 4 && nB == 3 ) return 2 ; 
  if ( nJ == 4 && nB == 4 ) return 3 ; 
  if ( nJ == 5 && nB == 2 ) return 4 ; 
  if ( nJ == 5 && nB == 3 ) return 5 ; 
  if ( nJ == 5 && nB >= 4 ) return 6 ; 
  if ( nJ >= 6 && nB == 2 ) return 7 ; 
  if ( nJ >= 6 && nB == 3 ) return 8 ; 
  if ( nJ >= 6 && nB >= 4 ) return 9 ; 

  return -1 ; 

}

std::string analyzer::_NameofCategory( int cate ){
  if( cate == 1  ) return "4j2b" ; 
  if( cate == 2  ) return "4j3b" ; 
  if( cate == 3  ) return "4j4b" ; 
  if( cate == 4  ) return "5j2b" ; 
  if( cate == 5  ) return "5j3b" ; 
  if( cate == 6  ) return "5j4b" ; 
  if( cate == 7  ) return "6j2b" ; 
  if( cate == 8  ) return "6j3b" ; 
  if( cate == 9  ) return "6j4b" ; 
  return "NoCategoryDefined" ;
}

void analyzer::SetFakeEstimationModeOn(){
  std::cout <<"analyzer.cc : Fake estimation is set on."<<std::endl;
  b_FakeEstimation = true; 
}

double analyzer::_getFakeSF( bool mu_channel  ){


  if( r->nJet == 0 ) return 0 ; 

  double nJ = r-> nJet > 6 ? 6 : r->nJet ;
  double nB = r->nBJet > 4 ? 4 : r->nBJet ;

  double sf = 0 ; 
  double sf_err = 0 ; 
  double sf_err_sign = 0 ; 

  if( ( ! isMC && isMuonStream ) || ( isMC && mu_channel ) ){ 
    // = data muon stream 
    sf     =  h_FakeSF_Mu ->GetBinContent( nJ+1 , nB+1 );
    sf_err =  h_FakeSF_Mu ->GetBinError  ( nJ+1 , nB+1 );
    sf_err_sign = syst_FakeMuStat == 0 ? 0 : ( syst_FakeMuStat > 0 ? +1 : -1 );
  }else{
    // = data electron stream, or MC_elChannel. 
    sf     =  h_FakeSF_El ->GetBinContent( nJ+1 , nB+1 );
    sf_err =  h_FakeSF_El ->GetBinError  ( nJ+1 , nB+1 );
    sf_err_sign = syst_FakeElStat == 0 ? 0 : ( syst_FakeElStat > 0 ? +1 : -1 );
  }

  return sf + sf_err_sign * sf_err ; 
}


bool analyzer::_passTtbarAdditionalJetID(){

  if( ttbarAdditionalJetID == noRequirement ) return true ; 

//  ttbb (xxx53 - xxx55): at least two additional b-jets
//  ttb (xxx51): exactly one b-jet, containing exactly one B-hadron
//  tt2b (xxx52): exactly one b-jet, containing at least two B-hadrons, mainly from collinear gluon splitting, very different treatment in simulation
//  ttcc (xxx41 - xxx45): at least one additional c-jet
//  ttlf (xxx00): additional jets are light-flavored jets

  const long id = r->id_additionalJetEventId % 100 ; 
  if( ttbarAdditionalJetID == ttbarPlusBBbar && (  53 <= id && id <= 55 ) ) return true ;
  if( ttbarAdditionalJetID == ttbarPlusB     && (  id  == 51 ) ) return true ;
  if( ttbarAdditionalJetID == ttbarPlus2B    && (  id  == 52 ) ) return true ;
  if( ttbarAdditionalJetID == ttbarPlusCCbar && (  41 <= id && id <= 45 ) ) return true ;
  if( ttbarAdditionalJetID == ttbarOther     && (  id  == 0 ) ) return true ;

  return false ; 

}

void analyzer::setTtbarAdditionalJetIDCut( TtbarAdditionalJetID id ){
  ttbarAdditionalJetID = id ; 
}

void analyzer::SetPUReweightSyst( int PU_syst ){

  if( PU_syst > 0 ){
    ttHSF . SetupDataPileupFile( "PileupHistogram_Moriond17_MinBiasPlus4p6_72383.root" );
  }else if ( PU_syst < 0 ){
    ttHSF . SetupDataPileupFile( "PileupHistogram_Moriond17_MinBiasMinus4p6_66017.root" );
  }


}


void analyzer::SetMCPileupChannel( std::string name ){

  ttHSF . SetMCPileupChannel( name );

}



void analyzer::EnableSkippingOddEventNumber(){

  flag_skip_oddEventNumber = true ;
  return ; 
}


float analyzer::_calcDR2( float eta1, float eta2, float phi1, float phi2 ){
  
  float d_eta = eta1 - eta2 ;
  float d_phi = fabs( phi1 - phi2 ) ;
  d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ;

  return  d_eta*d_eta + d_phi*d_phi ;

}


float analyzer::miniDR( float basis_eta ,  std::vector<float> * etas, float basis_phi , std::vector<float> * phis ){


  float mindr2 = 10000 ; 

  for( unsigned int i = 0 ; i < etas -> size() ; i++ ){
    
    float dr2 = _calcDR2( basis_eta , etas->at(i) , basis_phi , phis->at(i) );

    mindr2 = mindr2 < dr2 ? mindr2 : dr2 ; 

  }

  return sqrt( mindr2 ) ; 
}


void analyzer::SetDileptonFakeLeptonAnalysisMode(){

  b_DileptonFakeLeptonAnalysisMode = true ;
}

void analyzer::SetPeriod( int p ){

  period = p ; 

}


void analyzer::SetReader( reader * _r ){
  r = _r ;

  if( syst_SDMassResolution == -1 ){
    r -> SetSDMassResolition( -1 );
  }

  if( syst_SDMassScale != 0 ){
    r -> SetSDMassScale( syst_SDMassScale > 0 ? +1 : -1 );
  }

} 



void analyzer::work_for_truthInfo(){


  long N_Fatjet_PassingKinematicCriteria = 0 ;

  for( unsigned int iFat = 0 ; iFat < ( r -> fatjet_pt ->size() ) ; iFat ++ ){

    if( ( r ->  fatjet_pt -> at(iFat ) )  < 250 ){ continue ; }
    if( fabs ( r ->  fatjet_eta -> at(iFat ) )  > 2.0 ){ continue ; }

    int associateWithHiggs = 
      (
       r -> truth_higgs_pt > 200 
       &&
       _calcDR2( r -> fatjet_eta -> at( iFat ),  r-> truth_higgs_eta,  
		 r -> fatjet_phi -> at( iFat ),  r-> truth_higgs_phi  ) < 0.4 * 0.4 )
      ? 0 : 1 ; // 0 = associated

    FillHistogram( h_MinBTagDiscri_NoCut [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
    FillHistogram( h_DoubleB_NoCut       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );

    FillHistogram( h_SL_MinBTagDiscri_NoCut [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
    FillHistogram( h_SL_DoubleB_NoCut       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );
    FillHistogram2D ( h2_MinB_DoubleB_noCut [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r->fatjet_doubleb ->at(iFat) ,  r-> wgt_TOTAL ) ;

    // dR(lepton) cut                                                                                                                                                          
    const double dR2 =
      _calcDR2( r->fatjet_eta ->at( iFat ),
		r->lepton_eta ,
		r->fatjet_phi ->at( iFat ),
		r->lepton_phi ) ;
    if( dR2 < 0.8 * 0.8 ) continue ; // (overlapping with leptons)                                                                                                             
    // (Memo : LJ selection is required in the followings, where only one tight lepton exist.)


    N_Fatjet_PassingKinematicCriteria ++ ;

    const double handMergedBtagger = 
      ( r->fatjet_doubleb  ->at(iFat) < 0.3 ) 
      ? r->fatjet_doubleb  ->at(iFat) 
      : 0.3 + ( r->fatjet_subjet_minbtagger->at(iFat) < 0 ? 0 : r->fatjet_subjet_minbtagger->at(iFat) ) ; 


     eventnumber    = r -> EventNumber ; 
     eventcategory  = _EventCateBasedOnNjetNBtagJet(  r -> nJet , r -> nBJet ); 
     eventweight    = r-> wgt_TOTAL ;
     ak8jet_genhiggsmatch = ( associateWithHiggs == 0 ? 1 : 0 ) ; // In tree, flag==1 means matched.

     ak8jet_pt            = r -> fatjet_pt              -> at(iFat)  ;
     ak8jet_eta           = r -> fatjet_eta             -> at(iFat)  ;
     ak8jet_sdmass        = r -> fatjet_sdmass          -> at(iFat)  ;
     ak8jet_deepcsv_min   = r -> fatjet_subjet_minbtagger  -> at(iFat)  ;
     ak8jet_doubleb       = r -> fatjet_doubleb  -> at(iFat)  ; 
     ak8jet_tau21         = r -> fatjet_tau21  -> at(iFat)  ;
     ak8jet_tau32         = r -> fatjet_tau32  -> at(iFat)  ;

     if( (r -> EventNumber ) % 2 == 1 ) {
       outtree_odd -> Fill();
     }else{
       outtree_even -> Fill();
     }

     // Variables for MVA(BDT) INPUT
     tmva_ak8jet_pt          = r -> fatjet_pt              -> at(iFat);	 
     tmva_ak8jet_eta         = r -> fatjet_eta             -> at(iFat) ; 	 
     tmva_ak8jet_sdmass      = r -> fatjet_sdmass          -> at(iFat) ;	 
     {
       double minb = ( r -> fatjet_subjet_minbtagger  -> at(iFat) )  ; 
       minb = minb < 0 ? 0 : minb ; 
       tmva_ak8jet_deepcsv_min =  minb ; 
     }
     tmva_ak8jet_doubleb     = r -> fatjet_doubleb  -> at(iFat) ;    
     const double my_ak8bdt_value = tmva_reader -> EvaluateMVA("SatoshiAk8BDT");


    if( ( r -> pass_El || r -> pass_Mu ) && ( r -> nJet == 4 ) && ( r -> nBJet >= 4 ) ){
      FillHistogram( h_MinBTagDiscri_4j4b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_DoubleB_4j4b       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );

      FillHistogram( h_SL_MinBTagDiscri_4j4b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_SL_DoubleB_4j4b       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );

      FillHistogram2D ( h2_MinB_DoubleB_4j4b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r->fatjet_doubleb ->at(iFat) ,  r-> wgt_TOTAL ) ;
    }

    if( ( r -> pass_El || r -> pass_Mu ) && ( r -> nJet == 5 ) && ( r -> nBJet >= 4 ) ){
      FillHistogram( h_MinBTagDiscri_5j4b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_DoubleB_5j4b       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );

      FillHistogram( h_SL_MinBTagDiscri_5j4b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_SL_DoubleB_5j4b       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );

      FillHistogram2D ( h2_MinB_DoubleB_5j4b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r->fatjet_doubleb ->at(iFat) ,  r-> wgt_TOTAL ) ;
    }

    if( ( r -> pass_El || r -> pass_Mu ) && ( r -> nJet >= 6 ) && ( r -> nBJet == 2 ) ){
      FillHistogram( h_MinBTagDiscri_6j2b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_DoubleB_6j2b       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );

      FillHistogram( h_SL_MinBTagDiscri_6j2b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_SL_DoubleB_6j2b       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );

      FillHistogram2D ( h2_MinB_DoubleB_6j2b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r->fatjet_doubleb ->at(iFat) ,  r-> wgt_TOTAL ) ;
    }

    if( ( r -> pass_El || r -> pass_Mu ) && ( r -> nJet >= 6 ) && ( r -> nBJet == 3 ) ){
      FillHistogram( h_MinBTagDiscri_6j3b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_DoubleB_6j3b       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );

      FillHistogram( h_SL_MinBTagDiscri_6j3b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_SL_DoubleB_6j3b       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram2D ( h2_MinB_DoubleB_6j3b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r->fatjet_doubleb ->at(iFat) ,  r-> wgt_TOTAL ) ;
    }

    if( ( r -> pass_El || r -> pass_Mu ) && ( r -> nJet >= 6 ) && ( r -> nBJet >= 4 ) ){
      FillHistogram( h_MinBTagDiscri_6j4b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_DoubleB_6j4b       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_HandMergedBtag_6j4b    [associateWithHiggs] , handMergedBtagger ,  r-> wgt_TOTAL );

      if( (r -> EventNumber ) % 2 != 1 ){
	FillHistogram( h_ak8BDT_6j4b [associateWithHiggs] , my_ak8bdt_value , r-> wgt_TOTAL );
      }

      FillHistogram( h_SL_MinBTagDiscri_6j4b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_SL_DoubleB_6j4b       [associateWithHiggs] ,  r->fatjet_doubleb          ->at(iFat) ,  r-> wgt_TOTAL );
      FillHistogram( h_SL_HandMergedBtag_6j4b  [associateWithHiggs] , handMergedBtagger ,  r-> wgt_TOTAL );

      FillHistogram2D ( h2_MinB_DoubleB_6j4b [associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r->fatjet_doubleb ->at(iFat) ,  r-> wgt_TOTAL ) ;


      if( r -> fatjet_nsubjets -> at (iFat) == 2 ){
	FillHistogram2D ( h2_MinB_DoubleB_6j4b_twosubjet [0][associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r->fatjet_doubleb ->at(iFat) ,  r-> wgt_TOTAL ) ;
      }else {
	FillHistogram2D ( h2_MinB_DoubleB_6j4b_twosubjet [1][associateWithHiggs] ,  r->fatjet_subjet_minbtagger->at(iFat) ,  r->fatjet_doubleb ->at(iFat) ,  r-> wgt_TOTAL ) ;
      }

    }
    
  }


  if( ( r -> pass_El || r -> pass_Mu ) ) {

    const long cate = _EventCateBasedOnNjetNBtagJet(  r -> nJet , r -> nBJet );
    FillHistogram ( h_NFatJetPassingKinCriteria[cate] , N_Fatjet_PassingKinematicCriteria ,  r-> wgt_TOTAL ) ;
    
  }

  return ; 
}
