#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "histogram.hpp"
#include "constants.hpp"
#include "functions.hpp"
//#include <vector>
#include "THnSparse.h"

TFile *exp_pim_file;
TFile *sim_pim_file;
TFile *empty_pim_file;
TFile *sim_no_rad_pim_file;
TFile *exp_pro_file;
TFile *sim_pro_file;
TFile *empty_pro_file;
TFile *sim_no_rad_pro_file;
TFile *exp_pip_file;
TFile *sim_pip_file;
TFile *empty_pip_file;
TFile *sim_no_rad_pip_file;
TFile *hole_file;
TFile *hole_pos_file;
TFile *hole_neg_file;



#endif