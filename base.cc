#define base_cxx
#include "base.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <math.h>

#include <TH2F.h>

#include <TStyle.h>

#include <iostream>

void base::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


   }
}

