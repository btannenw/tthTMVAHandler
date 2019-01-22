
#include "controller.h"

#include <iostream>
#include <math.h>
#include <algorithm>

#include <TLorentzVector.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "systematics.h"

controller::controller ()
  : output("out.root")
  , outputdir("./")
  , JetPU("CHS")
  , isMC ( true )
  , isMuonStream( false )
  , b_FakeEstimationMode( false )
  , ttbarAdditionalJetID( analyzer::TtbarAdditionalJetID::noRequirement )
  , PU_syst(0)
  , flag_skip_oddEventNumber( false  )
  , idx_SFsystematic(0)
  , b_DileptonFakeLeptonAnalysisMode (false)
  , MCPileupChannelName("")
  , period( -1 ) 
  , N_maximum_processed_events( -1 )
{

}


void controller::setSFSystematic( int syst ){

  idx_SFsystematic = syst ; 

}


void controller::addSystematics( int syst ){

  Systematics.push_back( GetSystematicsName( syst ) );

}

controller::~controller (){
  
}

void controller::do_analyses(){
  
  if( Systematics.size() == 0 ){
    Systematics.push_back("");// nominal.
  }

  init();
  long total_processed_events = 0 ; 
  
  for( std::vector<std::string>::iterator f = filelist.begin();
       f != filelist.end();
       f++ ){
    
    std::cout << "Input file = \""<<(*f)<<"\" is opend." << std::endl ; 
    
    TFile * tf = TFile::Open( f->c_str() );

    std::vector< TTree *  > mytree ; 
    std::vector< reader * > readers ;
    for( unsigned int iSyst = 0 ; iSyst < Systematics.size() ; iSyst++ ){

      TTree * t = 0 ;
      char treename[100];
      sprintf(treename , "tree%s" ,( Systematics[iSyst] == "" ? std::string("") : ( "_" + Systematics[iSyst] )  ).c_str() );

      tf->GetObject(treename , t ) ;
      std::cout <<"Get TTree named = " << treename << " from file = " << (*f)
		<< "\n -> result : " << (t==0? "null pointer after GetObject. Something wrong." : "looks OK")
		<< std::endl;
      reader * r = new reader( t ) ;
      r -> setJetPUName( JetPU );
      r -> init();
      ana[ iSyst ] -> SetReader( r );

      mytree  . push_back( t );
      readers . push_back( r );
    }

    const long NENTRIES = mytree[0]->GetEntries();
    for (Long64_t jEntry=0; jEntry < NENTRIES 
	   &&
	   ( N_maximum_processed_events < 0 
	     ||
	     total_processed_events < N_maximum_processed_events 
	     )
	   ; jEntry ++ ){

      total_processed_events ++ ;

      for( unsigned int iSyst = 0 ; iSyst < readers.size() ; iSyst++ ){
	Long64_t ientry = readers[iSyst] -> LoadTree(jEntry);
	if (ientry < 0) break;
	readers[iSyst] -> InitEvent (jEntry);
      }

      for( unsigned int iSyst = 0 ; iSyst < readers.size() ; iSyst++ ){
	ana[iSyst]->AnalyzeEvent();
      }

    } // Events-loop ends

    for( unsigned int iSyst = 0 ; iSyst < Systematics.size() ; iSyst++ ){
      delete readers[iSyst];
      readers[iSyst] = 0 ; 
    }

    // (??? not needed because 'delete reader' kills TBranch? see 'base::~base()',, )
    // tf ->Close();

  } // file loop ends.

  postProcess();
}



void controller::init(){

  for( unsigned int iSyst = 0 ; iSyst < Systematics.size() ; iSyst++ ){

    std::cout <<"[controller] : Preparing an analyzer for systematics = "
	      << ( Systematics[iSyst] == ""  ? "nominal" : Systematics[iSyst] ) << std::endl ; 

    ana.push_back(  new analyzer() ) ;

    std::string syst_output_prefix = outputdir +"/" + "output" + Systematics[iSyst] + getFileNameExtention();
    ana[iSyst] -> SetOutputfileName( ( syst_output_prefix + "_" + output  +".root").c_str() );

    if( ! isMC ){
      ana[iSyst]->SetIsMuonStream( isMuonStream ) ;
    }
    if( b_FakeEstimationMode ){
      std::cout << "[controller.cc] Fake estimation mode will be on"<<std::endl;
      ana[iSyst]->SetFakeEstimationModeOn();
    }

    if( ttbarAdditionalJetID != analyzer::TtbarAdditionalJetID::noRequirement ){
      ana[iSyst]->setTtbarAdditionalJetIDCut( ttbarAdditionalJetID );
    }
    
    if( PU_syst != 0 ){
      ana[iSyst] -> SetPUReweightSyst ( PU_syst );
    }

    if( flag_skip_oddEventNumber ){
      ana[iSyst] -> EnableSkippingOddEventNumber();
    }

    if( idx_SFsystematic != 0 ){
      ana[iSyst] -> setSFSystematic( idx_SFsystematic );
    }

    if( b_DileptonFakeLeptonAnalysisMode ){
      ana[iSyst] -> SetDileptonFakeLeptonAnalysisMode () ; 
    }

    if( MCPileupChannelName != "" ){
      ana[iSyst] -> SetMCPileupChannel( MCPileupChannelName ) ; 
    }

    if( period >= 0 ){
      ana[iSyst] -> SetPeriod( period ) ;
    }


    ana[iSyst] -> init();

  }

}

void controller::postProcess(){

  for( unsigned int iSyst = 0 ; iSyst < Systematics.size() ; iSyst++ ){
    ana[iSyst] -> postProcess();
  }

}

void controller::useDataAnalysisMode(){
  isMC = false ; 
}

void controller::SetIsMuonStream( bool isMu ){
  isMuonStream = isMu ; 
  useDataAnalysisMode();
}

void controller::SetFakeEstimationModeOn(){
  b_FakeEstimationMode = true ; 
}


void controller::SetTtbarAdditionalJetIDCut( analyzer::TtbarAdditionalJetID id ){
  ttbarAdditionalJetID = id ; 
}

void controller::SetPUsyst( int _pu_syst ){
  PU_syst = _pu_syst;
}


std::string controller::getFileNameExtention(){

  if( PU_syst > 0 ) return std::string("PUUp");
  if( PU_syst < 0 ) return std::string("PUDown");

  if( idx_SFsystematic == 1  ){ return std::string("muTrigSFUP"  ); }
  if( idx_SFsystematic == 2  ){ return std::string("muTrigSFDown"); }
  if( idx_SFsystematic == 3  ){ return std::string("elTrigSFUP"  ); } 
  if( idx_SFsystematic == 4  ){ return std::string("elTrigSFDown"); } 
  if( idx_SFsystematic == 5  ){ return std::string("muIDSFUP"  ); }
  if( idx_SFsystematic == 6  ){ return std::string("muIDSFDown"); }
  if( idx_SFsystematic == 7  ){ return std::string("elIDSFUP"  ); }
  if( idx_SFsystematic == 8  ){ return std::string("elIDSFDown"); }
  if( idx_SFsystematic == 9  ){ return std::string("muIsoSFUP"  ); } 
  if( idx_SFsystematic == 10 ){ return std::string("muIsoSFDown"); } 
  if( idx_SFsystematic == 11 ){ return std::string("elRecoSFUP"  ); }
  if( idx_SFsystematic == 12 ){ return std::string("elRecoSFDown"); }

  if( idx_SFsystematic == 20 ){ return std::string("SDMassScaleUp"); }
  if( idx_SFsystematic == 21 ){ return std::string("SDMassScaleDown"); }
  if( idx_SFsystematic == 22 ){ return std::string("SDMassResolutionWorse"); }

  if( idx_SFsystematic == 101 ){ return std::string("FakeMuStatUp"); }	
  if( idx_SFsystematic == 102 ){ return std::string("FakeMuStatDown"); }
  if( idx_SFsystematic == 103 ){ return std::string("FakeElStatUp"); }	
  if( idx_SFsystematic == 104 ){ return std::string("FakeElStatDown"); }


  // ME uncertainty                
  if( idx_SFsystematic == 202 ){ return std::string( "ME_muR_DOUBLE_muF_NOM" ) ; } 
  if( idx_SFsystematic == 203 ){ return std::string( "ME_muR_HALF_muF_NOM" ) ; } 	   
  if( idx_SFsystematic == 204 ){ return std::string( "ME_muR_NOM_muF_DOUBLE" ) ; }    
  if( idx_SFsystematic == 205 ){ return std::string( "ME_muR_DOUBLE_muF_DOUBLE" ) ; } 
  if( idx_SFsystematic == 206 ){ return std::string( "ME_muR_HALF_muF_DOUBLE" ) ; }   
  if( idx_SFsystematic == 207 ){ return std::string( "ME_muR_NOM_muF_HALF" ) ; } 	   
  if( idx_SFsystematic == 208 ){ return std::string( "ME_muR_DOUBLE_muF_HALF" ) ; }   
  if( idx_SFsystematic == 209 ){ return std::string( "ME_muR_HALF_muF_HALF" ) ; }     

  // PS uncertainty 
  if( idx_SFsystematic == 301 ){ return std::string("PS_REDUCED_ISRUP" ) ; }  	  
  if( idx_SFsystematic == 302 ){ return std::string("PS_REDUCED_FSRUP" ) ; }  	    	  
  if( idx_SFsystematic == 303 ){ return std::string("PS_REDUCED_ISRDOWN" ) ; }  	  	  
  if( idx_SFsystematic == 304 ){ return std::string("PS_REDUCED_FSRDOWN" ) ; }  	  	  
  if( idx_SFsystematic == 305 ){ return std::string("PS_DEFAULT_ISRUP" ) ; }  	    	  
  if( idx_SFsystematic == 306 ){ return std::string("PS_DEFAULT_FSRUP" ) ; }  	    	  
  if( idx_SFsystematic == 307 ){ return std::string("PS_DEFAULT_ISRDOWN" ) ; }  	  	  
  if( idx_SFsystematic == 308 ){ return std::string("PS_DEFAULT_FSRDOWN" ) ; }  	  	  
  if( idx_SFsystematic == 309 ){ return std::string("PS_CONSERVATIVE_ISRUP" ) ; }  	    
  if( idx_SFsystematic == 310 ){ return std::string("PS_CONSERVATIVE_FSRUP" ) ; }  	    
  if( idx_SFsystematic == 311 ){ return std::string("PS_CONSERVATIVE_ISRDOWN" ) ; }  	  
  if( idx_SFsystematic == 312 ){ return std::string("PS_CONSERVATIVE_FSRDOWN" ) ; }  	  





  return std::string("");

}


void controller::EnableSkippingOddEventNumber(){

  flag_skip_oddEventNumber = true ;
  return ; 
}

void controller::SetDileptonFakeLeptonAnalysisMode(){
  b_DileptonFakeLeptonAnalysisMode = true ; 
  return ; 
}

void controller::SetMCPileupChannel( std::string _name ){

  MCPileupChannelName = _name ;

}

void controller::SetPeriod( int p ) {

  period = p ; 

}

void controller::SetMaxNEvents( long _N_maximum_processed_events ){

  N_maximum_processed_events = _N_maximum_processed_events ;
  std::cout <<"Maximum number of events to be proceesed is set to : " << N_maximum_processed_events << std::endl ; 
  
}
