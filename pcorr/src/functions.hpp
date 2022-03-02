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
#include <math.h>
#include "flags.hpp"
#include "TFile.h"
#include "TChain.h"

namespace fun {
	std::shared_ptr<TFile> Name_File(std::shared_ptr<Flags> flag_);
	float theta_calc(float theta_p_, float beam_energy_, bool radians_=false);
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
	TLorentzVector Set_k_mu(int run_);
	float W(TLorentzVector k_mu_prime_, int run_);
	TLorentzVector Make_4Vector(float p, float theta, float phi, float m);
	void loadChain(std::shared_ptr<TChain> chain_, std::string file_, int thread_id_, int max_, std::shared_ptr<Flags> flags_);
	std::vector<std::string> read_file_list(std::string path, int thread_num, std::shared_ptr<Flags> flags_);
	void print_vector_idx(std::vector<int> vec_);
	float p_calc_e(float theta_e_, float beam_);
	float delta_p(float theta_e_, float beam_, float p_e_);
}

#endif