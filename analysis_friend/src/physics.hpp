#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD

#include "constants.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "histogram.hpp"
#include "THnSparse.h"
#include "TH1.h"
#include "TMath.h"


namespace physics {
double Virtual_Photon_Flux(double W_, double Q2_, double E_);
double Luminosity(double Qtot_, double corr_factor_=1.0);
//double Radiative_Corr();
double Error(double N_gen_, double N_rec_, double weight_sum_);//Weighted
double Error(double N_gen_, double N_rec_);//Unweighted
};


#endif