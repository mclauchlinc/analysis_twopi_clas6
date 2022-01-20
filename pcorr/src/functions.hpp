#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "constants.hpp"
#include "branches.hpp"
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <TLorentzVector.h>
#include "TVector3.h"

namespace fun {
	float theta_calc(float theta_p_, float beam_energy_);
	float delta_theta(float theta_e_, float theta_p_, float beam_energy_);
	float poly_4(float x_, float a_, float b_, float c_, float d_, float e_);
	float poly_3(float x_, float a_, float b_, float c_, float d_);
	float poly_2(float x_, float a_, float b_, float c_);
	float theta(std::shared_ptr<Branches> data_, int idx_, bool radians_=false);
	int get_sector(float phi_);
	int get_sector(std::shared_ptr<Branches> data_, int idx_);
	float phi_center( float phi_);
	float phi(std::shared_ptr<Branches> data_, int idx_, bool center_=false, bool radians_=false);
	float p_calc(float theta_e_, float beam_energy_);
	float delta_p_e(float p_e_, float theta_e_, float beam_energy_);
	float W(TLorentzVector k_mu_prime_, int run_);
	TLorentzVector Make_4Vector(float p, float theta, float phi, float m);
}

#endif