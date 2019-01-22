
// Macro to draw TH1F histograms.
//
// usage :
//    root 'DrawTH1FPlot.C+("output_LPC/output")'

#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>


int DrawTH1FPlot( std::string dirpath ){
  
  bool detail_stat_info_inLegend = false;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  std::vector<std::string> inputs;
  std::vector<std::string> names ; 

  names.push_back( "CHS ttH(bb)" );inputs.push_back( dirpath + "../output_CHS/tth.root");
  names.push_back( "CHS ttbar  " );inputs.push_back( dirpath + "../output_CHS/ttbar.root");
  names.push_back( "PUPPI ttH(bb)" );inputs.push_back( dirpath + "../output_PUPPI/tth.root");
  names.push_back( "PUPPI ttbar  " );inputs.push_back( dirpath + "../output_PUPPI/ttbar.root");

  long myColor[] = { 
    kBlue,
    kRed,
    kBlue,
    kRed
  };

  long myLineStyle[] = { 
    1, // solid
    1,
    2, // dash
    2
  };
  
  struct hist{
    std::string histname;
    std::string epsname;
    std::string xAxisName;
    
    hist( const char * v1, const char * v2 , const char *v3 ){
      histname .assign( v1 );
      epsname  .assign( v2 );
      xAxisName.assign( v3 );
    };
  };
  std::vector< hist > histograms ; 

  histograms . push_back( hist( "BDT" , "BDT", "BDT output" ) );

  histograms . push_back( hist( "CSVv2_1thLeadingJet" , "CSVv2_1thLeadingJet", "1st Leading Jet CSVv2" ) );
  histograms . push_back( hist( "CSVv2_2thLeadingJet" , "CSVv2_2thLeadingJet", "2nd Leading Jet CSVv2" ) );
  histograms . push_back( hist( "CSVv2_3thLeadingJet" , "CSVv2_3thLeadingJet", "3rd Leading Jet CSVv2" ) );
  histograms . push_back( hist( "CSVv2_4thLeadingJet" , "CSVv2_4thLeadingJet", "4th Leading Jet CSVv2" ) );
  histograms . push_back( hist( "CSVv2_5thLeadingJet" , "CSVv2_5thLeadingJet", "5th Leading Jet CSVv2" ) );
  histograms . push_back( hist( "CSVv2_6thLeadingJet" , "CSVv2_6thLeadingJet", "6th Leading Jet CSVv2" ) );
  histograms . push_back( hist( "CSVv2_7thLeadingJet" , "CSVv2_7thLeadingJet", "7th Leading Jet CSVv2" ) );
  histograms . push_back( hist( "CSVv2_8thLeadingJet" , "CSVv2_8thLeadingJet", "8th Leading Jet CSVv2" ) );

  histograms . push_back( hist( "cMVA_1thLeadingJet" , "cMVA_1thLeadingJet", "1st Leading Jet cMVA" ) );
  histograms . push_back( hist( "cMVA_2thLeadingJet" , "cMVA_2thLeadingJet", "2nd Leading Jet cMVA" ) );
  histograms . push_back( hist( "cMVA_3thLeadingJet" , "cMVA_3thLeadingJet", "3rd Leading Jet cMVA" ) );
  histograms . push_back( hist( "cMVA_4thLeadingJet" , "cMVA_4thLeadingJet", "4th Leading Jet cMVA" ) );
  histograms . push_back( hist( "cMVA_5thLeadingJet" , "cMVA_5thLeadingJet", "5th Leading Jet cMVA" ) );
  histograms . push_back( hist( "cMVA_6thLeadingJet" , "cMVA_6thLeadingJet", "6th Leading Jet cMVA" ) );
  histograms . push_back( hist( "cMVA_7thLeadingJet" , "cMVA_7thLeadingJet", "7th Leading Jet cMVA" ) );
  histograms . push_back( hist( "cMVA_8thLeadingJet" , "cMVA_8thLeadingJet", "8th Leading Jet cMVA" ) );

  std::vector< TFile * > files ; 
  for( std::vector<std::string>::iterator filepath =inputs.begin() ;
       filepath != inputs.end() ;
       filepath ++ ){
    std::cout <<"Opening file : " << (*filepath) << std::endl ; 
    files .push_back( TFile::Open( filepath->c_str() ) );
    std::cout <<"Opening file : " << (*filepath) <<" done." << std::endl ; 
  }
  std::cout <<"Opening all files done. \n"<< std::endl ; 


  TCanvas c;

  for( std::vector<hist>::iterator h = histograms.begin()  ;
       h != histograms.end();
       h ++ ){

    std::cout << h->histname << " " << h->xAxisName << std::endl ; 
    TH1F * _h_[ inputs.size() ] ;

    int idx_MaxHist = 0 ; 
    for( int i =0 ; i < inputs.size() ; i++ ){
      files[i] -> GetObject(  h->histname.c_str()  ,  _h_[ i ] ) ;
      _h_[i]->SetLineColor( myColor[i] );
      _h_[i]->SetLineStyle( myLineStyle[i] );

      if( _h_[ i ]           -> DrawNormalized () ->GetMaximum() 
	  >
	  _h_[ idx_MaxHist ] -> DrawNormalized () ->GetMaximum() 
	  ){  idx_MaxHist = i ; }
    }


    c.Clear();

    TLegend leg(0.4,0.6,0.75,0.9);
    _h_[ idx_MaxHist ] -> GetXaxis() -> SetTitle( h->xAxisName.c_str() );
    _h_[ idx_MaxHist ] -> DrawNormalized("H");
    for( int i = 0 ; i < inputs.size() ; i++  ){

      _h_[ i ] -> DrawNormalized("sameH");

      leg.AddEntry( _h_[i], names[i].c_str() ,"le");

      if( detail_stat_info_inLegend ){
	TObject * nullobj = 0 ; 
	
	char info[100];
	float mean   =  _h_[i] -> GetMean();
	float mean_e =  _h_[i] -> GetMeanError();
	sprintf( info, "(Mean = %3.2f +- %3.2f)", mean , mean_e );
	leg.AddEntry(  nullobj , info ,"");
	float rms   =  _h_[i] -> GetRMS();
	float rms_e =  _h_[i] -> GetRMSError();
	sprintf( info, "(RMS  = %3.2f +- %3.2f)",  rms , rms_e) ;
	leg.AddEntry(  nullobj , info ,"");
      }
    }
    leg.Draw();

    c.Print(  (h->epsname +".eps") .c_str() );

  }


  return 0 ; 
}
