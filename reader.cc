

#include "reader.h"

#include <iostream>
#include <algorithm>

Int_t reader::InitEvent( Long64_t event){


  const Int_t val = GetEntry( event ); // <= function of base.h

  // Prepare variable which are calculated from input files, if you have any.

  if( EventCategory >= 0 ){
    //    MVAValue = TMVAreader[ EventCategory ]->EvaluateMVA( "myBDT" ) ; 
  }else{
    //    MVAValue = -2;
  }


  if( syst_SDMassScale != 0 || syst_SDMassResolution != 0 ){


    // Init SD-mass scale factors
    double scalefactors[ fatjet_sdmass -> size() ];
    for( unsigned int i = 0 ; i < fatjet_sdmass -> size() ; i ++ ){
      scalefactors[i] = 1 ; 
    }   

    if( syst_SDMassScale != 0 ){ // Take care SD mass scale factor 

      for( unsigned int i = 0 ; i < fatjet_sdmass -> size() ; i ++ ){
	scalefactors[i] = 1 + 0.0094 * syst_SDMassScale ;  // The value taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging?rev=57
      }

    }else{  // Take care SD mass resolution
      // reaching here emans syst_SDMassResolution is not zero.

      for( unsigned int i = 0 ; i < fatjet_sdmass -> size() ; i ++ ){

	if( fatjet_sdmass ->at(i) != 0 ){

	  myrandom . SetSeed( i + EventNumber + ( (long) ( met_phi * 1000 ) ) + ( (long) ( fatjet_eta -> at(i) * 1000 ) ) ) ;
	  const double sigma = 10.1 / fatjet_sdmass ->at(i);  // The value taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging?rev=57
	  const double gauss = myrandom.Gaus( 0 , sigma ) ; 
	  
	  const double scalefactor = 0.66 ;  
	  // which comes from sqrt ( max ( SF^2 - 1  , 0 )  where SF = 1(nominal) + 0.2(uncertainty) = 1.2.
	  // The value taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging?rev=57
	  
	  const double scale = 1 + scalefactor * gauss ; 
	  scalefactors[i] = scale > 0 ? scale : 0 ;
	}

      }

    } 
      
    // Apply scale factors to each jets
    for( unsigned int i = 0 ; i < fatjet_sdmass -> size() ; i ++ ){
      fatjet_sdmass ->at(i) *= scalefactors[i] ;
    }   


  }


  return val ; 
}

double reader::getMVAValue(){
  return MVAValue ; 
}


reader::reader(TTree *tree) : base( tree )
			    , JetPU( "CHS" )
			    , syst_SDMassScale (0)
			    , syst_SDMassResolution (0)

{

}


void reader::SetSDMassResolition( int sys ){
  syst_SDMassResolution = sys ; 
}

void reader::SetSDMassScale( int sys ){
  syst_SDMassScale = sys ; 
}


void reader::init(){

//  for( int i = 0 ; i < 8 ; i ++ ){
//    TMVAreader[i] = new TMVA::Reader();
//
//    TMVAreader[i]->AddVariable( "lepton_eta" , & lepton_eta );
//    TMVAreader[i]->AddVariable( "lepton_pt"  , & lepton_pt );
//    TMVAreader[i]->AddVariable( "lepton_phi" , & lepton_phi );
//
//
//    TMVAreader[i]->AddVariable( "jet_1_eta" , & jet_1_eta );
//    TMVAreader[i]->AddVariable( "jet_2_eta" , & jet_2_eta );
//    TMVAreader[i]->AddVariable( "jet_3_eta" , & jet_3_eta );
//    TMVAreader[i]->AddVariable( "jet_4_eta" , & jet_4_eta );
//    if( i >= 3 ){
//      TMVAreader[i]->AddVariable( "jet_5_eta" , & jet_5_eta );
//      if( i >=  5 ){
//	TMVAreader[i]->AddVariable( "jet_6_eta" , & jet_6_eta );
//      }
//    }
//    TMVAreader[i]->AddVariable( "jet_1_pt"  , & jet_1_pt  ); 
//    TMVAreader[i]->AddVariable( "jet_2_pt"  , & jet_2_pt  ); 
//    TMVAreader[i]->AddVariable( "jet_3_pt"  , & jet_3_pt  ); 
//    TMVAreader[i]->AddVariable( "jet_4_pt"  , & jet_4_pt  ); 
//    if( i >= 3 ){
//      TMVAreader[i]->AddVariable( "jet_5_pt"  , & jet_5_pt  ); 
//      if( i >=  5 ){
//	TMVAreader[i]->AddVariable( "jet_6_pt"  , & jet_6_pt  ); 
//      }
//    }
//
//    TMVAreader[i]->AddVariable( "jet_1_phi"  , & jet_1_phi  ); 
//    TMVAreader[i]->AddVariable( "jet_2_phi"  , & jet_2_phi  ); 
//    TMVAreader[i]->AddVariable( "jet_3_phi"  , & jet_3_phi  ); 
//    TMVAreader[i]->AddVariable( "jet_4_phi"  , & jet_4_phi  ); 
//    if( i >= 3 ){
//      TMVAreader[i]->AddVariable( "jet_5_phi"  , & jet_5_phi  ); 
//      if( i >=  5 ){
//	TMVAreader[i]->AddVariable( "jet_6_phi"  , & jet_6_phi  ); 
//      }
//    }
//
//    TMVAreader[i]->AddVariable( "jet_1_m"  , & jet_1_m  ); 
//    TMVAreader[i]->AddVariable( "jet_2_m"  , & jet_2_m  ); 
//    TMVAreader[i]->AddVariable( "jet_3_m"  , & jet_3_m  ); 
//    TMVAreader[i]->AddVariable( "jet_4_m"  , & jet_4_m  ); 
//    if( i >= 3 ){
//      TMVAreader[i]->AddVariable( "jet_5_m"  , & jet_5_m  ); 
//      if( i >=  5 ){
//	TMVAreader[i]->AddVariable( "jet_6_m"  , & jet_6_m  ); 
//      }
//    }
//
//    if( false ) // skip this part because the DTB is not used. But I keep this for future reference.
//    {//
//    char filepath[100];
//    sprintf( filepath , 
//	     "weightFiles/2016_04_28_MCStatBugFixed/TMVAClassification_%s_ttHcategory%d_BDTG.weights.xml",
//	     JetPU.c_str(),
//	     i );
//    TMVAreader[i]->BookMVA( "myBDT", filepath );
//    }//
//
//  }

}
