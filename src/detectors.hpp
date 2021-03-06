#ifndef DETECTORS_HPP
#define DETECTORS_HPP

#include "physics.hpp"
#include "branches.hpp"
#include "TVector.h"
#include "constants.hpp"


const static int _num_cc_segments_ = 18;


namespace detect{
	//float Acc = -0.000785;
	//float Bcc = 0.0;
	//float Ccc = -0.00168;
	//float Dcc = 1.0;
//CC functions
	int cc_segment(int cc_segm);
	int cc_lrc(int cc_segm);
	TVector3 cc_sector_center(TVector3 v_);
	TVector3 cc_sector_center(float x_, float y_, float z_);
	float cc_theta(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_);
	float cc_theta(std::shared_ptr<Branches> data_, int idx_);
	float cc_phi(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_);
	float cc_phi(std::shared_ptr<Branches> data_, int idx_);
	float cc_x(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_, int sec_);
	float cc_y(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_, int sec_);
	float cc_x(std::shared_ptr<Branches> data_, int idx_);
	float cc_y(std::shared_ptr<Branches> data_, int idx_);

	//SC functions
	float sc_theta(std::shared_ptr<Branches> data_, int idx_);
	float sc_phi(std::shared_ptr<Branches> data_, int idx_);

	//EC functions
	float ec_theta(std::shared_ptr<Branches> data_, int idx_);
	float ec_phi(std::shared_ptr<Branches> data_, int idx_);

}



#endif