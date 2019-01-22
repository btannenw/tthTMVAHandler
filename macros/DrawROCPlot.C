
// Macro to draw ROC curbe.
// - Two files input. Each has a TH1F histogram with the same name.
//
// usage :
//    root 'DrawROCPlot.C+("../output_CHS/")'


#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>

int DrawROCPlot( std::string dirpath ){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  std::vector<std::string> inputs;
  std::vector<std::string> names ; 

  names.push_back( "ttH(bb)" );inputs.push_back( dirpath + "tth.root");
  names.push_back( "ttbar  " );inputs.push_back( dirpath + "ttbar.root");

  long myColor[] = { 
    kBlue,
    kMagenta,
  };
  

  // Definition of the direction of the two axies.
  enum AxisDefinition{
    efficiency  ,
    rejection_power // = 1 - efficiency.
  };

  struct graphInput {
    std::string eps ;
    std::string histogramname ;
    std::string xAxisName;
    AxisDefinition xAxisDef ;
    std::string yAxisName;
    AxisDefinition yAxisDef ;
    std::pair<double, double> range_x , range_y ; 

    TGraph * roc ; // Graph to be created by this program.

    graphInput( const char *epsname , const char * hist ,
		const char *v_Xaxis,
		const AxisDefinition xAxisDef_ ,
		const char *v_Yaxis,
		const AxisDefinition yAxisDef_ ,
		double x_min , double x_max , 
		double y_min , double y_max 
		){
      eps .assign( epsname );
      histogramname .assign( hist );
      xAxisName.assign( v_Xaxis );
      xAxisDef = xAxisDef_ ;
      yAxisName.assign( v_Yaxis );
      yAxisDef = yAxisDef_ ;
      range_x  = std::pair<double, double> ( x_min , x_max );
      range_y  = std::pair<double, double> ( y_min , y_max );

    };
  };

  std::vector< graphInput > graphInputs ; 
  graphInputs. push_back( graphInput( 
				     "BDT_ROC" , 
				     "BDT_for_ROC",
				     "ttH signal",
				     efficiency , // X axis
				     "ttbar rejection",
				     rejection_power , // Y axix
				     0.5, 1,  // x range
				     0.5, 1  // y range 
				      ) );


  for( int i = 0 ; i < 8 ; i ++  ){ 
    char histname [100];
    char epsname [100];
    sprintf( histname, "BDT_Category%d_forROC", i );
    sprintf( epsname, "BDT_ROC_category%d", i );
    graphInputs. push_back( graphInput( epsname , histname,
					"ttH signal", // X axis
					efficiency , // X axis
					"ttbar rejection",
					rejection_power , // Y axix
					0.5, 1,  // x range
					0.5, 1  // y range 
					) );
  }


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
    
    std::cout <<"TFile : Get Object " << setting->histogramname << std::endl;
    TH1F * h_1 ; 
    files.at(0)->GetObject( setting->histogramname.c_str() , h_1);

    TH1F * h_2 ; 
    files.at(1)->GetObject( setting->histogramname.c_str() , h_2);
    
    const long NBINS = h_1 -> GetNbinsX() ; 
    double roc_x [ NBINS +2 ] , roc_y[ NBINS +2 ]; // "+2" for overflow and underflow bins.
    double bin_values [ NBINS ];

    for( long iBin = 0 ; iBin <= NBINS +1 ; iBin ++ ){ 
      // bin = 0 : underflow, 
      // bin = NBINS +1  : overflow

      if(  iBin != 0 && iBin != NBINS +1 ){
	// Skip underflow-bin and overflow-bin;
	bin_values[iBin-1] = h_1 -> GetBinLowEdge( iBin );
      }

      { // X-axis
	double total    = h_1 -> Integral( );
	double integral = h_1 -> Integral( iBin , NBINS+1 ); // sum the right side of the histogram.
	roc_x[iBin] = ( setting->xAxisDef == efficiency ? ( integral ) : ( total - integral) ) / total;
      }
      { // Y-axis
	double total    = h_2 -> Integral( );
	double integral = h_2 -> Integral( iBin , NBINS+1 ); // sum the right side of the histogram.N_
	roc_y[iBin] = ( setting->yAxisDef == efficiency ? ( integral ) : ( total - integral) ) / total;
      }

    }
    TGraph * g = new TGraph( NBINS + 2 , roc_x , roc_y ) ;
    setting -> roc = g ; 

    c.Clear();
    TH2F h_frame_allRange ("","",  1, 0, 1, 1, 0 , 1 );
    h_frame_allRange . GetXaxis() -> SetTitle(setting -> xAxisName.c_str() ) ;
    h_frame_allRange . GetYaxis() -> SetTitle(setting -> yAxisName.c_str() ) ;
    h_frame_allRange . Draw();
    g->Draw("l");
    c.Print( ( setting ->eps +"_FullRange.eps") .c_str() );
    c.Print("allRatioPliot.pdf");

    c.Clear();
    TH2F h_frame ("","", 
		  1,  setting->range_x.first, setting->range_x.second,
		  1,  setting->range_y.first, setting->range_y.second
		  );
    h_frame . GetXaxis() -> SetTitle(setting -> xAxisName.c_str() ) ;
    h_frame . GetYaxis() -> SetTitle(setting -> yAxisName.c_str() ) ;
    h_frame . Draw();
    g->Draw("l");
    c.Print( ( setting ->eps +".eps") .c_str() );
    c.Print("allRatioPliot.pdf");


    TH2F h_frame_2 ("","", 
		    1, bin_values[0] , bin_values[ NBINS-1 ] ,
		    1, 0, 1
		    );
    h_frame_2 . Draw();
    TGraph * g_1 = new TGraph( NBINS  , bin_values , &( roc_x[1] ) ) ;
    TGraph * g_2 = new TGraph( NBINS  , bin_values , &( roc_y[1] ) ) ;
    g_1->SetLineColor( myColor[0]);
    g_2->SetLineColor( myColor[1]);
    g_1->Draw("l");
    g_2->Draw("l");

    TLegend leg(0.75,0.5,0.9,0.7);
    leg.AddEntry(g_1, names[0].c_str() ,"l");
    leg.AddEntry(g_2, names[1].c_str() ,"l");
    leg.Draw();

    c.Print( ( setting ->eps +".eps") .c_str() );

  }

  // - - - - - - - - - - - - - - - - - - - 
  // - - - - - Draw extra plot - - - - - -
  // - - - - - - - - - - - - - - - - - - -
  {

    // Draw the graphs 1-8 ( = category_0-7) into one canvas.
    long myColor[9] = {
      0, //  dummy
      kRed + 3 , 
      kRed  -7 ,
      kMagenta -7 ,
      kViolet +2 , 
      kBlue -4 ,
      kAzure + 6,
      kTeal  - 3 , 
      kGreen +2  
    } ;
    
//    TH2F h_frame ("","", 
//		  1,  graphInputs[1].range_x.first, graphInputs[1].range_x.second,
//		  1,  graphInputs[1].range_y.first, graphInputs[1].range_y.second
//		  );
    TH2F h_frame ("","", 
		  1, 0 , 1 , 
		  1, 0 , 1 
		  );
    h_frame . GetXaxis() -> SetTitle( graphInputs[1] . xAxisName.c_str() ) ;
    h_frame . GetYaxis() -> SetTitle( graphInputs[1] . yAxisName.c_str() ) ;
    h_frame . Draw();
    TLegend leg(0.1,0.1,0.3,0.5);
    for( int i = 1 ; i <= 8 ; i++ ){
      graphInputs[i] . roc -> SetLineColor( myColor[i] );
      graphInputs[i] . roc -> Draw("l");

      char legendname[100];
      sprintf( legendname , "Category %d", i -1  );
      leg.AddEntry( graphInputs[i] . roc , legendname  ,"l");
    }
    leg.Draw();
    c.Print("ROC_allCategories.eps");
    c.Print("allRatioPliot.pdf");
  }
  // - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - -

  c.Print("allRatioPliot.pdf]");
  return 0 ; 
}
