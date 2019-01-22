
// Macro to draw ratio plot, such as efficiency.
// The input ROOT files is expected to have two histograms -- nominator and denominator ;
//  e.g. h_AllMuons, h_SelectedMuons
// 
//
// usage :
//    root 'macro_DrawRatioPlot.C+("hoge/")'


#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>

int DrawRatioPlot( std::string dirpath ){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  std::vector<std::string> inputs;
  std::vector<std::string> names ; 

  names.push_back( "CHS" );inputs.push_back( dirpath + "output/CHS.root");
  names.push_back( "PUPPI" );inputs.push_back( dirpath + "output/PUPPI.root");

  long myColor[] = { 
    kGreen,
    kOrange,
    kYellow,
    kRed,
    kMagenta,
    kBlue,
    kBlack,
    kGray
  };
  
  struct graphInput {
    std::string eps ;
    std::string total ;
    std::string pass ;
    std::string xAxisName;
    std::string yAxisName;
    std::pair<double, double> range_x , range_y ; 
    std::pair<double, double> legend_x , legend_y ; 

    graphInput( const char *epsname , const char * v_total, const char * v_pass  ,
		const char *v_Xaxis,
		const char *v_Yaxis,
		double x_min , double x_max , 
		double y_min , double y_max ,
		double leg_x_min , double leg_x_max , 
		double leg_y_min , double leg_y_max 
		){
      eps .assign( epsname );
      total .assign( v_total );
      pass  .assign( v_pass );
      xAxisName.assign( v_Xaxis );
      yAxisName.assign( v_Yaxis );
      
      range_x  = std::pair<double, double> ( x_min , x_max );
      range_y  = std::pair<double, double> ( y_min , y_max );
      legend_x = std::pair<double, double> ( leg_x_min , leg_x_max );
      legend_y = std::pair<double, double> ( leg_y_min , leg_y_max );

    };
  };

  std::vector< graphInput > graphInputs ; 
  graphInputs. push_back( graphInput( 
				     "JetRecon_Purity_eta" , 
				     "AllRecoJetWithCut_eta",
				     "GenAssociatedRecoJetWithCut_eta",
				     "Reconstructed jet #eta",
				     "Purity",
				     -5, 5,  // x range
				     0.8, 1, // y range 
				     0.4,  0.6, 0.1, 0.2// legend
				      ) );

  graphInputs. push_back( graphInput( 
				     "JetRecon_Efficiency_eta" , 
				     "AllGenJetWithCut_eta",
				     "RecoAssociatedGenJetWithCut_eta",
				     "Gen jet #eta",
				     "Efficiency",
				     -5, 5,  // x range
				     0.7, 1 , // y range 
				     0.4,  0.6, 0.1, 0.2// legend
				      ) );



  graphInputs. push_back( graphInput( 
				     "JetRecon_ForwardJet_Efficiency_pt" , 
				     "ForwardGenJetWithCut_pt",
				     "RecoAssociatedForwardGenJetWithCut_pt",
				     "Forward gen jet Pt",
				     "Efficiency",
				     0, 200,  // x range
				     0.7, 1, // y range 
				     0.4,  0.6, 0.1, 0.2// legend
				      ) );


  graphInputs. push_back( graphInput( 
				     "JetRecon_BarrelJet_Efficiency_pt" , 
				     "BarrelGenJetWithCut_pt",
				     "RecoAssociatedBarrelGenJetWithCut_pt",
				     "Barrel gen jet Pt",
				     "Efficiency",
				     0, 200,  // x range
				     0.7, 1, // y range 
				     0.4,  0.6, 0.1, 0.2// legend
				      ) );



  graphInputs. push_back( graphInput( 
				     "JetRecon_ForwardJet_Efficiency_NPV" , 
				     "ForwardGenJetWithCut_npv",
				     "RecoAssociatedForwardGenJetWithCut_npv",
				     "NPV",
				     "Forward Jet Efficiency",
				     0, 40,  // x range
				     0.7, 1 , // y range 
				     0.7,  0.9, 0.1, 0.2// legend
				      ) );


  graphInputs. push_back( graphInput( 
				     "JetRecon_BarrelJet_Efficiency_NPV" , 
				     "BarrelGenJetWithCut_npv",
				     "RecoAssociatedBarrelGenJetWithCut_npv",
				     "NPV",
				     "Barrel Jet Efficiency",
				     0, 40,  // x range
				     0.7, 1 , // y range 
				     0.7,  0.9, 0.1, 0.2// legend
				      ) );


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

  c.Print("allRatioPliot.pdf[");

  for( std::vector< graphInput >::iterator setting = graphInputs.begin()  ;
       setting != graphInputs.end();
       setting ++ ){

    c.Clear();
    
    TH2F h_frame ("","", 
		  1,  setting->range_x.first, setting->range_x.second,
		  1,  setting->range_y.first, setting->range_y.second
		  );
    h_frame . GetXaxis() -> SetTitle(setting -> xAxisName.c_str() ) ;
    h_frame . GetYaxis() -> SetTitle(setting -> yAxisName.c_str() ) ;
    h_frame . Draw();
    
    TLegend leg( setting->legend_x.first, setting->legend_y.first , setting->legend_x.second, setting->legend_y.second );
    // for each input file, draw one TGraph.
    
    long iFile = -1 ; 
    for( std::vector< TFile * >::iterator f = files.begin(); f != files.end(); f++ ){
      iFile ++ ; 
	TH1F * pass ;
	TH1F * total ;

	(*f)->GetObject( setting->total.c_str() , total);
	(*f)->GetObject( setting->pass .c_str() , pass );

	TGraphAsymmErrors * g = new TGraphAsymmErrors( pass , total );
	
	g->SetLineColor( myColor[iFile ] );
	g->Draw("P");
	leg.AddEntry( g, names[iFile].c_str() ,"lep");

    }

    leg.Draw();
    c.Print( ( setting ->eps +".eps") .c_str() );
    c.Print("allRatioPliot.pdf");

  }

  c.Print("allRatioPliot.pdf]");
  return 0 ; 
}
