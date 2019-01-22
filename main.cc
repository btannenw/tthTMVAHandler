
#include "controller.h"

#include <unistd.h>
#include <iostream>
#include <stdlib.h>

#include <fstream>

int main(int argc,char *argv[]){

  controller my_controller;

  TH1::SetDefaultSumw2();

  {
    int c ;
    while(( c = getopt(argc,argv,"i:o:Ds:S:d:EMFt:pPOIQ:x:N:")) != -1 ){
      
      switch( c){
      case 'x' :
	{
	  my_controller . SetPeriod( std::stoi(optarg) );
	  break;
	}
      case 'N' :
	{
	  my_controller . SetMaxNEvents( std::stoi(optarg) );
	  break;
	}
      case 'Q' :
	{
	  my_controller . SetMCPileupChannel( std::string(optarg) );
	  break;
	}
      case 'i' :
	{
	  my_controller . addInputFileName( std::string(optarg) );
	  break;
	}
      case 'o' :
	{
	  my_controller . setOutputFileName( std::string(optarg) );
	  break ; 
	}
      case 'D' :
	{
	  my_controller . useDataAnalysisMode();
	  break ; 
	}
      case 'd' :
	{
	  my_controller . setOutputDir(  std::string(optarg) );
	  break ; 
	}
      case 's' :
	{
	  // systematic variation which use YggdrasilSyst Ttree.
	  my_controller . addSystematics( std::stoi(optarg) );
	  break ; 
	}

      case 'S' :
	{
	  // Systematic variation which does not use YggdrasilSyst Ttree.
	  my_controller . setSFSystematic( std::stoi(optarg) );
	  break ; 
	}

      case 'E' :
	{
	  my_controller . SetIsMuonStream(false);
	  break ; 
	}
      case 'M' :
	{
	  my_controller . SetIsMuonStream(true);
	  break ; 
	}
      case 'F' :
	{
	  my_controller . SetFakeEstimationModeOn( );
	  break ; 
	}
      case 't' :
	{
	  const int ttbarId = std::stoi(optarg) ;
	  if( ttbarId == 1 ){ my_controller . SetTtbarAdditionalJetIDCut( analyzer::TtbarAdditionalJetID::ttbarPlusB     ) ; }
	  if( ttbarId == 2 ){ my_controller . SetTtbarAdditionalJetIDCut( analyzer::TtbarAdditionalJetID::ttbarPlus2B    ) ; }
	  if( ttbarId == 3 ){ my_controller . SetTtbarAdditionalJetIDCut( analyzer::TtbarAdditionalJetID::ttbarPlusBBbar ) ; }
	  if( ttbarId == 4 ){ my_controller . SetTtbarAdditionalJetIDCut( analyzer::TtbarAdditionalJetID::ttbarPlusCCbar ) ; }
	  if( ttbarId == 5 ){ my_controller . SetTtbarAdditionalJetIDCut( analyzer::TtbarAdditionalJetID::ttbarOther     ) ; }
	  break ; 
	}
      case 'p' :
	{
	  my_controller . SetPUsyst( + 1 );
	  break ; 
	}
      case 'P' :
	{
	  my_controller . SetPUsyst( - 1 );
	  break ; 
	}
      case 'O' :
	{
	  my_controller . EnableSkippingOddEventNumber();
	  break ; 
	}
      case 'I' :
	{
	  my_controller . SetDileptonFakeLeptonAnalysisMode();
	  break ; 
	}
      }//end switch
    }//end while-Loop
  } //end scope.

  my_controller.do_analyses ();

  return 0 ; 

}
