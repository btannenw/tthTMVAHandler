#ifndef MY_CONTROLLER
#define MY_CONTROLLER


#include <vector>
#include <string>
#include <map>

#include <TH1F.h>
#include <TH2F.h>

#include "analyzer.h"

class controller {

 private : 

  std::vector<std::string> filelist  ;
  std::string output ; 
  std::string outputdir;
  const char * JetPU ;
  std::string PUReweighting_DataFile, PUReweighting_MCFile;

 public : 

  controller ();
  ~controller ();

  inline void addInputFileName ( std::string name ){ filelist.push_back(name); };
  inline void setOutputFileName( std::string name ){ output.assign(name); };
  inline void setOutputDir( std::string name ){ outputdir.assign(name); };
  inline void setJetPU( const char * name ){ JetPU = name ;};
  void addSystematics( int syst );

  void do_analyses();

  void useDataAnalysisMode();
  void SetIsMuonStream( bool isMu );

  void SetFakeEstimationModeOn();

  void SetTtbarAdditionalJetIDCut( analyzer::TtbarAdditionalJetID id );

  void SetPUsyst( int _pu_syst );

  void EnableSkippingOddEventNumber();

  void setSFSystematic( int syst );

  void SetDileptonFakeLeptonAnalysisMode();

  void SetMCPileupChannel( std::string _name );

  void SetPeriod( int p ) ;

  void SetMaxNEvents( long _N_maximum_processed_events );

 private :

  void init();
  void postProcess();

  std::vector< analyzer * > ana ; 

  bool isMC;
  bool isMuonStream;

  std::vector<std::string> Systematics;

  bool b_FakeEstimationMode ;

  analyzer::TtbarAdditionalJetID ttbarAdditionalJetID;

  int PU_syst;

  std::string getFileNameExtention();

  bool flag_skip_oddEventNumber ; 

  int idx_SFsystematic ; 
  

  bool b_DileptonFakeLeptonAnalysisMode ; 

  std::string MCPileupChannelName;

  int period ;

  long N_maximum_processed_events ;

};

#endif
