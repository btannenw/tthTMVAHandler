


// Marco to draw a graph of the mean/RMS of a set of histograms.
//
//  For example, you make histograms of something for each N_PV :
//    you have TH1F Mll_NPVbin[x], where x is NPV(e.g. 0-20),
//  and would like to draw a graph of RMS as a function of N_PV.
//
// usage :
//    root 'Draw_ProfilePlot.C+("output_LPC/output")'


#include "analyzer.cc"
#include "base.cc"
#include "base_event.cc"
#include "base_met.cc"
#include "base_muon.cc"
#include "reader.cc"
#include "reader_muon.cc"

#include <iostream>

#include <TFile.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>

int Draw_ProfilePlot( std::string dirpath ){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  std::vector<std::string> filepaths, names ;

  names.push_back( "PFMET" );filepaths.push_back( dirpath + "/PFMET_AK4PFPUPPI.root");
  names.push_back( "PFCHSMET" );filepaths.push_back( dirpath + "/PFCHSMET_AK4PFPUPPI.root");

  names.push_back( "PFCHSMET uncorrect" );filepaths.push_back( "output_UncorrectedMet/output/PFCHSMET_AK4PFPUPPI.root");

  names.push_back( "76x PUPPI" );filepaths.push_back( dirpath + "/76xTuning_AK4PFPUPPI.root");
  names.push_back( "BarrelOnlyTuned" );filepaths.push_back( dirpath + "/BarrelOnlyTuned_AK4PFPUPPI.root");
  names.push_back( "BarrelTuningAndForwardExpSF" );filepaths.push_back( dirpath + "/BarrelTuningAndForwardExpSF_AK4PFPUPPI.root");
  names.push_back( "quarterTuning" );filepaths.push_back( dirpath + "/quorter80xTuning_AK4PFPUPPI.root");
  names.push_back( "halfTuning" );filepaths.push_back( dirpath + "/half80xTuning_AK4PFPUPPI.root");
  names.push_back( "80x PUPPI" );filepaths.push_back( dirpath + "/80xTuning_AK4PFPUPPI.root");

  long myColor[] = { 
    kBlue-10 , 
    kBlack,
    kGray,
    kGreen,
    kYellow,
    kOrange,
    kRed,
    kMagenta,
    kBlue
  };


  enum UsedInfo {
    mean ,
    rms 
  };
  
  struct plotstyle{
    std::string plotTH1Fname;
    std::string plotEPSname;
    std::string x_axis ; 
    std::string y_axis ; 

    UsedInfo info ;

    double frame_x_min, frame_x_max ;
    double frame_y_min, frame_y_max ;

    double legend_x_min, legend_x_max ;
    double legend_y_min, legend_y_max ;
    
    long NbinsX;
    double x_bin[100];
    double x_bin_err[100];
  };


  const long NPLOT = 5 ; 
  plotstyle plots[ NPLOT ];

  const int idx_paraResolution = 0;
  plots[0].plotTH1Fname.assign("recoil_para");
  plots[0].plotEPSname.assign("MetResponse_Parallel");
  plots[0].x_axis.assign("Z_{P_{T}}[GeV]");
  plots[0].y_axis.assign("RMS of - u_{#parallel} - P_{T}^{Z} [GeV]");
  plots[0].frame_x_min = analyzer::Zpt_Binning[0] ;
  plots[0].frame_x_max = analyzer::Zpt_Binning[ analyzer::N_Zpt_Binning ];
  plots[0].frame_y_min = 8 ;
  plots[0].frame_y_max = 40 ; 
  plots[0].legend_x_min = 0.1 ; 
  plots[0].legend_x_max = 0.8 ;
  plots[0].legend_y_min = 0.7 ;
  plots[0].legend_y_max = 0.9 ;
  plots[0].info = rms ; 
  plots[0].NbinsX = analyzer::N_Zpt_Binning ;
  for( int i = 0 ; i <  plots[0].NbinsX ; i ++ ){ 
    plots[0].x_bin[i]     =  ( analyzer::Zpt_Binning[ i + 1 ] + analyzer::Zpt_Binning[ i ] ) / 2;
    plots[0].x_bin_err[i] =  ( analyzer::Zpt_Binning[ i + 1 ] - analyzer::Zpt_Binning[ i ] ) / 2 ;
  }




  plots[1].plotTH1Fname.assign("recoil_perp");
  plots[1].plotEPSname.assign("MetResponse_Perpendicular");
  plots[1].x_axis.assign("Z_{P_{T}}[GeV]");
  plots[1].y_axis.assign("RMS of u_{#perp} [GeV]");
  plots[1].frame_x_min = analyzer::Zpt_Binning[0] ;
  plots[1].frame_x_max = analyzer::Zpt_Binning[ analyzer::N_Zpt_Binning ];
  plots[1].frame_y_min = 8 ;
  plots[1].frame_y_max = 30 ; 
  plots[1].legend_x_min = 0.1 ; 
  plots[1].legend_x_max = 0.8 ;
  plots[1].legend_y_min = 0.7 ;
  plots[1].legend_y_max = 0.9 ;
  plots[1].info = rms ; 
  plots[1].NbinsX = analyzer::N_Zpt_Binning ;
  for( int i = 0 ; i <  plots[1].NbinsX ; i ++ ){ 
    plots[1].x_bin[i]     =  ( analyzer::Zpt_Binning[ i + 1 ] + analyzer::Zpt_Binning[ i ] ) / 2;
    plots[1].x_bin_err[i] =  ( analyzer::Zpt_Binning[ i + 1 ] - analyzer::Zpt_Binning[ i ] ) / 2 ;
  }

  const int idx_response = 2;
  plots[2].plotTH1Fname.assign("response");
  plots[2].plotEPSname.assign("MetResponse_ZptDependence");
  plots[2].x_axis.assign("Z_{P_{T}}[GeV]");
  plots[2].y_axis.assign("Mean of (-u_{#parallel})/Z_{P_{T}} [GeV]");
  plots[2].frame_x_min = analyzer::Zpt_Binning[0] ;
  plots[2].frame_x_max = analyzer::Zpt_Binning[ analyzer::N_Zpt_Binning ];
  plots[2].frame_y_min = 0 ;
  plots[2].frame_y_max = 1.1 ; 
  plots[2].legend_x_min = 0.5 ; 
  plots[2].legend_x_max = 0.9 ;
  plots[2].legend_y_min = 0.1 ;
  plots[2].legend_y_max = 0.4 ;
  plots[2].info = mean ; 
  plots[2].NbinsX = analyzer::N_Zpt_Binning ;
  for( int i = 0 ; i <  plots[2].NbinsX ; i ++ ){ 
    plots[2].x_bin[i]     =  ( analyzer::Zpt_Binning[ i + 1 ] + analyzer::Zpt_Binning[ i ] ) / 2;
    plots[2].x_bin_err[i] =  ( analyzer::Zpt_Binning[ i + 1 ] - analyzer::Zpt_Binning[ i ] ) / 2 ;
  }


  plots[3].plotTH1Fname.assign("recoil_para_Npv");
  plots[3].plotEPSname.assign("MetResponse_Parallel_NpvDependence");
  plots[3].x_axis.assign("NPV");
  plots[3].y_axis.assign("RMS of - u_{#parallel} - P_{T}^{Z} [GeV]");
  plots[3].frame_x_min = analyzer::Npv_Binning[0] ;
  plots[3].frame_x_max = analyzer::Npv_Binning[ analyzer::N_Npv_Binning ];
  plots[3].frame_y_min = 0 ;
  plots[3].frame_y_max = 20 ; 
  plots[3].legend_x_min = 0.5 ; 
  plots[3].legend_x_max = 0.9 ;
  plots[3].legend_y_min = 0.1 ;
  plots[3].legend_y_max = 0.4 ;
  plots[3].info = rms ; 
  plots[3].NbinsX = analyzer::N_Npv_Binning ;
  for( int i = 0 ; i <  plots[3].NbinsX ; i ++ ){ 
    plots[3].x_bin[i]     =  ( analyzer::Npv_Binning[ i + 1 ] + analyzer::Npv_Binning[ i ] ) / 2;
    plots[3].x_bin_err[i] =  ( analyzer::Npv_Binning[ i + 1 ] - analyzer::Npv_Binning[ i ] ) / 2 ;
  }


  plots[4].plotTH1Fname.assign("response_Npv");
  plots[4].plotEPSname.assign("MetResponse_NpvDependence");
  plots[4].x_axis.assign("NPV");
  plots[4].y_axis.assign("Mean of (-u_{#perp})/Z_{P_{T}} [GeV]");
  plots[4].frame_x_min = analyzer::Npv_Binning[0] ;
  plots[4].frame_x_max = analyzer::Npv_Binning[ analyzer::N_Npv_Binning ];
  plots[4].frame_y_min = 0 ;
  plots[4].frame_y_max = 1.1 ; 
  plots[4].legend_x_min = 0.5 ; 
  plots[4].legend_x_max = 0.9 ;
  plots[4].legend_y_min = 0.1 ;
  plots[4].legend_y_max = 0.4 ;
  plots[4].info = mean ; 
  plots[4].NbinsX = analyzer::N_Npv_Binning ;
  for( int i = 0 ; i <  plots[4].NbinsX ; i ++ ){ 
    plots[4].x_bin[i]     =  ( analyzer::Npv_Binning[ i + 1 ] + analyzer::Npv_Binning[ i ] ) / 2;
    plots[4].x_bin_err[i] =  ( analyzer::Npv_Binning[ i + 1 ] - analyzer::Npv_Binning[ i ] ) / 2 ;
  }

  TCanvas c ;
  c.Print( "AllPlots_response.pdf[" ) ;

  std::vector< std::vector< TGraphErrors * > > all_graphs ;

  for( int i_plot = 0 ; i_plot < NPLOT ; i_plot ++ ){

  std::vector< TGraphErrors * > graphs ;

  for( std::vector<std::string>::iterator file = filepaths.begin();
       file != filepaths.end();
       file ++ ){

    double y_bins[ analyzer::N_Zpt_Binning ], y_bins_err[ analyzer::N_Zpt_Binning ];

    TFile * openfile = TFile::Open( file->c_str() );
    for( int i = 0 ; i < analyzer::N_Zpt_Binning ; i++ ){
      TH1F * h ; 
      char name[100];
      sprintf( name , "%s_bin%d", plots[i_plot].plotTH1Fname.c_str() , i );
      std::cout <<"GetObject " << name << " : from file = "  << *file << std::endl ; 
      openfile -> GetObject( name , h ) ; 

      if( plots[ i_plot ].info == mean ){
	y_bins[ i ]     = h -> GetMean();
	y_bins_err[ i ] = h -> GetMeanError();
      }else if( plots[ i_plot ].info == rms ){
	y_bins[ i ]     = h -> GetRMS();
	y_bins_err[ i ] = h -> GetRMSError();
      }else{
	std::cout << "ERROR : plot i="<<i_plot << " does not have appropriate \"UsedInfo\". Fill zero into the graph..." << std::endl ;
	y_bins[ i ]   = 0;
	y_bins_err[ i ] = 0;
      }
      std::cout << y_bins[ i ] <<" +- " << y_bins_err[ i ] << std::endl ; 
    }
    
    graphs . push_back( new TGraphErrors( plots[ i_plot ].NbinsX , plots[i_plot].x_bin, y_bins, plots[i_plot].x_bin_err, y_bins_err ) );

  }// File Loop ends.


  TH2F frame ( "f","f" , 1 , plots[i_plot].frame_x_min , plots[i_plot].frame_x_max , 1 , plots[i_plot].frame_y_min , plots[i_plot].frame_y_max );
  frame.GetXaxis()->SetTitle( plots[i_plot].x_axis.c_str() );
  frame.GetYaxis()->SetTitle( plots[i_plot].y_axis.c_str() );
  frame.Draw();
  TLegend leg( plots[i_plot].legend_x_min , plots[i_plot].legend_y_min ,  plots[i_plot].legend_x_max , plots[i_plot].legend_y_max );
  leg . SetNColumns(2);
  for( int iFile = 0 ; iFile < filepaths.size() ; iFile++ ){
    graphs[iFile] -> SetLineColor( myColor[iFile]  );
    graphs[iFile] -> Draw("P");
    leg.AddEntry( graphs[iFile] , names[iFile].c_str() ,"le");
  }
  leg.Draw();


  c.Print( ( plots[i_plot].plotEPSname +".eps" ).c_str());
  c.Print( "AllPlots_response.pdf" ) ;

  all_graphs . push_back( graphs ) ;

  }




  // - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - Draw Responce-corrected resolution plot :
  // - - - - - - - - - - - - - - - - - - - - - - - -
  {

    plotstyle ResponseCorrectedResolutionP;

    ResponseCorrectedResolutionP.plotTH1Fname.assign("");// not used
    ResponseCorrectedResolutionP.plotEPSname.assign("ResponceCorrectedMetResponse_Parallel");
    ResponseCorrectedResolutionP.x_axis.assign("Z_{P_{T}}[GeV]");
    ResponseCorrectedResolutionP.y_axis.assign(" RMS(-u_{#parallel} - P_{T}^{Z}) / Mean( -u_{#parallel}/Z_{P_{T}} ) [GeV]");
    ResponseCorrectedResolutionP.frame_x_min = analyzer::Zpt_Binning[0] ;
    ResponseCorrectedResolutionP.frame_x_max = analyzer::Zpt_Binning[ analyzer::N_Zpt_Binning ];
    ResponseCorrectedResolutionP.frame_y_min = 10 ;
    ResponseCorrectedResolutionP.frame_y_max = 45 ; 
    ResponseCorrectedResolutionP.legend_x_min = 0.1 ; 
    ResponseCorrectedResolutionP.legend_x_max = 0.8 ;
    ResponseCorrectedResolutionP.legend_y_min = 0.7 ;
    ResponseCorrectedResolutionP.legend_y_max = 0.9 ;
    ResponseCorrectedResolutionP.info = rms ;  // not used.


    TH2F frame ( "f","f" ,
		 1 , ResponseCorrectedResolutionP.frame_x_min , ResponseCorrectedResolutionP.frame_x_max ,
		 1 , ResponseCorrectedResolutionP.frame_y_min , ResponseCorrectedResolutionP.frame_y_max );
    frame.GetXaxis()->SetTitle( ResponseCorrectedResolutionP.x_axis.c_str() );
    frame.GetYaxis()->SetTitle( ResponseCorrectedResolutionP.y_axis.c_str() );
    frame.Draw();
    
    TLegend leg( ResponseCorrectedResolutionP.legend_x_min , ResponseCorrectedResolutionP.legend_y_min ,
		 ResponseCorrectedResolutionP.legend_x_max , ResponseCorrectedResolutionP.legend_y_max );
    
    leg . SetNColumns(2);
    
    double * X  = all_graphs[ idx_paraResolution ] [0] -> GetX ();
    double * Xe = all_graphs[ idx_paraResolution ] [0] -> GetEX ();

    for( unsigned int iFile = 0 ; iFile < filepaths.size() ; iFile++  ){

      double * Resolution   = all_graphs[ idx_paraResolution ] [ iFile ] -> GetY ();
      double * Resolution_e = all_graphs[ idx_paraResolution ] [ iFile ] -> GetEY ();
      
      double * Responce   = all_graphs[ idx_response ] [iFile ] -> GetY ();
      double * Responce_e = all_graphs[ idx_response ] [iFile ] -> GetEY ();
      
      double ResponseCorrected_Resolution   [ analyzer::N_Zpt_Binning ] ;
      double ResponseCorrected_Resolution_e [ analyzer::N_Zpt_Binning ] ;

      for( int iBin = 0 ; iBin < analyzer::N_Zpt_Binning ; iBin ++ ){

	ResponseCorrected_Resolution[ iBin ] = 
	  Resolution[ iBin ] / Responce [ iBin ]; 

	std::cout <<"ResponseCorrected_Resolution " << Resolution[ iBin ]  << " and " << Responce [ iBin ]  << std::endl ; 
	std::cout <<"ResponseCorrected_Resolution " << iBin << " = " << ResponseCorrected_Resolution[ iBin ] << std::endl ; 

	// Assuming the errors of Resolution and Responce are independent,
	//   e ( x/y ) = sqrt ( ( e(x)/y )^2 + ( x/(y^2) e(y))^2 ).
	ResponseCorrected_Resolution_e[ iBin ] = 
	  sqrt( 
	       pow ( Resolution_e[ iBin ] / Responce [ iBin ] , 2 )
	       +
	       pow ( Resolution[ iBin ] / pow( Responce [ iBin ] , 2 )  * Responce_e [ iBin ]  ,2 )
		);

      }

      TGraphErrors *g = new TGraphErrors( analyzer::N_Zpt_Binning , X , ResponseCorrected_Resolution , Xe, ResponseCorrected_Resolution_e );

      g -> SetLineColor( myColor[iFile]  );
      g -> Draw("P");
      leg.AddEntry( g , names[iFile].c_str() ,"le");
    } // end iFile Loop

    leg.Draw();
    c.Print( ( ResponseCorrectedResolutionP.plotEPSname +".eps" ).c_str() );
    c.Print( "AllPlots_response.pdf" ) ;

  } //end of scope for drawing of responce corrected resolution.

  c.Print( "AllPlots_response.pdf]" ) ;
  
  return 0 ; 
}
