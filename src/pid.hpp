#ifndef PID_HPP
#define PID_HPP

#include "constants.hpp"
#include "branches.hpp"
#include "cuts.hpp"


namespace pid {
	//Particle ID
	std::vector<bool> pid(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool pid_ele(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool pid_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool pid_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool pid_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	//ID Cut
	bool id_bank(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_,const char* _species_);
	//Fiducial Cuts
	std::vector<bool> fid(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool fid_ele(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool fid_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool fid_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool fid_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	//Delta T Cuts
	std::vector<bool> delta_t(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool delta_t_ele(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool delta_t_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool delta_t_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool delta_t_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	//Sampling Fraction
	bool sf(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	//Min CC
	bool min_cc(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	//Min EC
	bool min_ec(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	//Sanity Cuts
	bool sanity(int particle_, int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool sanity_ele(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool sanity_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool sanity_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	bool sanity_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	//Vertex Cuts
	bool vertex_e(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	//Efficiency Cuts
	bool sc_eff(int par_, int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
}

#endif