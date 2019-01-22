
#ifndef MY_READER_H
#define MY_READER_H

#include "base.h"

#include <TMVA/Reader.h>

#include <TRandom3.h>

class reader : public base {

 private :

  std::string JetPU ;

 public :

 reader(TTree *tree=0);
  ~reader (){} ;

  void init();

  Int_t InitEvent( Long64_t event ) ;

  double getMVAValue();

  inline void setJetPUName( const char * name ){ JetPU.assign(name); };

  void SetSDMassResolition( int sys );
  void SetSDMassScale( int sys );

 private :

  double MVAValue;

  TMVA::Reader * TMVAreader[8] ; 

  int syst_SDMassScale ;
  int syst_SDMassResolution ;

  TRandom3 myrandom ; 

};

#endif
