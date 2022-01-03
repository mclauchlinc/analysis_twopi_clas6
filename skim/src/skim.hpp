/**************************************/
/*										*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#ifndef SKIM_H_GUARD
#define SKIM_H_GUARD
#include <iostream>
#include "TChain.h"
#include "TFile.h"
#include "constants.hpp"

class Skim {
 private:
  TChain *chain;
  std::string fout;
  std::string fin;
  TFile *RootOutputFile;
  int GPART = 1;

 public:
  Skim(std::string input, std::string output);
  void Basic();
  ~Skim();
};
#endif