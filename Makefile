
# .SUFFIX = .so

CPP             = g++
CPPFLAGS        = -g -O1 -Wall -fPIC -D_REENTRANT -Wno-deprecated -I. 

ROOTCFLAGS      := $(shell root-config --cflags)
ROOTLIBS        := $(shell root-config --libs) -lMinuit -lEG -lTMVA #-lg2c
# ROOTLIBS        += $(shell echo ${MY_LDLib})
ROOTGLIBS       := $(shell root-config --glibs)

CPPFLAGS        += $(ROOTCFLAGS)
LIBS            = $(ROOTLIBS) -lm -L.
GLIBS           = $(ROOTGLIBS)

TARGET = main

OBJS4CLASSFIRE= MEMClassifier.o BDTClassifier.o CommonBDTvars.o AngularVariables.o Parameters.o Integrand.o JetLikelihood.o
LIB4CLASSFIER= \
	-L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gsl/2.2.1/lib  -lgsl -lgslcblas \
	-L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/openloops/1.2.4-ikhhed/lib  -lopenloops  \
        -L/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/lib   -lLHAPDF \
	-L../Cuba-4.2 -lcuba \
	-lMathMore 
# versions were taken form MEIntegratorStandalone/deps/gsl.xml :  <environment name="GSL_BASE" default="/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gsl/2.2.1"/>
# libMathMore : from ROOT package, for "VegasParameters"

# LHAPDF's path was determined from the environmental variable LHAPDF_DATA_PATH (echo ${LHAPDF_DATA_PATH})

INC4CLASSFIER= -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gsl/2.2.1/include -I../ 
CPPFLAGS += $(INC4CLASSFIER)

OBJS =   base.o reader.o main.o analyzer.o controller.o  $(OBJS4CLASSFIRE) ttHYggdrasilScaleFactors.o ttHYggdrasilEventSelection.o



all: ${TARGET} ttHYggdrasilScaleFactors_data
	echo "- - - - - - - - - - - - - - \n- - - - Successfully - - - - \n- - - - - - - - - - - - - -"

${TARGET}: $(OBJS)
	$(CPP) -o $@ $(CPPFLAGS) $(OBJS) $(LIBS)  $(LIB4CLASSFIER) $(EXTERNAL_OBJS)

main.o :  main.cc controller.h analyzer.h  base.h reader.h ttHYggdrasilScaleFactors.h ttHYggdrasilEventSelection.h
	$(CPP) -c $(CPPFLAGS) main.cc

controller.o : controller.cc controller.h analyzer.h base.h reader.h systematics.h ttHYggdrasilScaleFactors.h ttHYggdrasilEventSelection.h
	$(CPP) -c $(CPPFLAGS) controller.cc

analyzer.o : analyzer.cc  analyzer.h  base.h reader.h ttHYggdrasilScaleFactors.h ttHYggdrasilEventSelection.h
	$(CPP) -c $(CPPFLAGS)  analyzer.cc

reader.o : reader.cc base.h reader.h
	$(CPP) -c $(CPPFLAGS) reader.cc

base.o : base.cc base.h 
	$(CPP) -c $(CPPFLAGS) base.cc

systematics.h : ../ttHAnalyzer/systematics.h
	ln -s ../ttHAnalyzer/systematics.h systematics.h

clean: 
	-rm ${TARGET}
	-rm ${OBJS}
	-rm *~
	-rm systematics.h
	-rm ttHYggdrasilScaleFactors_data
	-rm ttHYggdrasilScaleFactors.cc 
	-rm ttHYggdrasilScaleFactors.h
	-rm ttHYggdrasilEventSelection.cc
	-rm ttHYggdrasilEventSelection.h

ttHYggdrasilScaleFactors_data : 
	ln -s  ../ttH-LeptonPlusJets/YggdrasilTreeMaker/data ttHYggdrasilScaleFactors_data

ttHYggdrasilScaleFactors.cc :
	ln -s  ../ttH-LeptonPlusJets/YggdrasilTreeMaker/plugins/ttHYggdrasilScaleFactors.cc ttHYggdrasilScaleFactors.cc

ttHYggdrasilScaleFactors.h :
	ln -s  ../ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilScaleFactors.h  ttHYggdrasilScaleFactors.h

ttHYggdrasilEventSelection.cc : 
	ln -s ../ttH-LeptonPlusJets/YggdrasilTreeMaker/plugins/ttHYggdrasilEventSelection.cc ttHYggdrasilEventSelection.cc

ttHYggdrasilEventSelection.h :
	ln -s ../ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilEventSelection.h ttHYggdrasilEventSelection.h


#BDTClassifier.o : ../TTH/CommonClassifier/src/BDTClassifier.cpp
#	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) ../TTH/CommonClassifier/src/BDTClassifier.cpp	
#
#MEMClassifier.o : ../TTH/CommonClassifier/src/MEMClassifier.cc
#	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) ../TTH/CommonClassifier/src/MEMClassifier.cc
#

AngularVariables.o   : ../TTH/CommonClassifier/src/AngularVariables.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

BDTClassifier.o	      : ../TTH/CommonClassifier/src/BDTClassifier.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

BDTClassifierHCRetrained.o   : ../TTH/CommonClassifier/src/BDTClassifierHCRetrained.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

BlrBDTClassifier.o	      : ../TTH/CommonClassifier/src/BlrBDTClassifier.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

CommonBDTvars.o	      : ../TTH/CommonClassifier/src/CommonBDTvars.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

DLBDTClassifier.o	      : ../TTH/CommonClassifier/src/DLBDTClassifier.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

DLBDTVars.o		      : ../TTH/CommonClassifier/src/DLBDTVars.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

DNNClassifier.o	      : ../TTH/CommonClassifier/src/DNNClassifier.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

MEMClassifier.o		      : ../TTH/CommonClassifier/src/MEMClassifier.cc
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

MemBDTClassifier.o	      : ../TTH/CommonClassifier/src/MemBDTClassifier.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

MemBDTClassifierV2.o	      : ../TTH/CommonClassifier/src/MemBDTClassifierV2.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

RecoLikelihoodVariables.o    : ../TTH/CommonClassifier/src/RecoLikelihoodVariables.cpp   
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

Parameters.o : ../TTH/MEIntegratorStandalone/src/Parameters.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

Integrand.o : ../TTH/MEIntegratorStandalone/src/Integrand.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<

JetLikelihood.o : ../TTH/MEIntegratorStandalone/src/JetLikelihood.cpp
	$(CPP) -c $(CPPFLAGS) $(INC4CLASSFIER) $<