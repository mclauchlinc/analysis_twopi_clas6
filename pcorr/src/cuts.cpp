#include "cuts.hpp"

bool cuts::W_cut(float W_, int run_){
	if(W_ > _W_cut_param_[0] && W_ < _W_cut_param_[1]){
		return true;
	}
	return false;
}

float cuts::fid_e_theta_cut(int run_, int sector_, float p_){ 
	return _c1e_[run_][sector_-1] + _c2e_[run_][sector_-1] / ((float)p_+_p_shift_e_[run_][sector_-1]);
}
float cuts::fid_e_expon_pos(int run_, int sector_, float p_){ 
	return _c3e_[run_][sector_-1][0] * TMath::Power((float)p_,_factor_e_[run_][sector_-1][0]);
}
float cuts::fid_e_expon_neg(int run_, int sector_, float p_){ 
	return _c3e_[run_][sector_-1][1] * TMath::Power((float)p_,_factor_e_[run_][sector_-1][1]);
}
float cuts::fid_e_del_phi_pos(int run_, int sector_, float p_, float theta_){ 
	return _c4e_[run_][sector_-1][0] * TMath::Power((TMath::Sin((theta_-cuts::fid_e_theta_cut(run_,sector_,p_))/_degree_)),cuts::fid_e_expon_pos(run_,sector_,p_));
}
float cuts::fid_e_del_phi_neg(int run_, int sector_, float p_, float theta_){ 
	return -_c4e_[run_][sector_-1][1] * TMath::Power((TMath::Sin((theta_-cuts::fid_e_theta_cut(run_,sector_,p_))/_degree_)),cuts::fid_e_expon_neg(run_,sector_,p_));
}

bool cuts::fid_e (float p_, float theta_, float phi_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	if(flags_->Flags::Fid_Cut(0) && !isnan(phi_) && !isnan(theta_)){
		//Calculate angles and sector
		float phi_c = fun::phi_center(phi_);
		int sector = fun::get_sector(phi_);

		//Application of Cut
		pass = true; 
		pass &= (phi_c<=cuts::fid_e_del_phi_pos(flags_->Flags::Run(), sector, p_, theta_));
		pass &= (phi_c >=cuts::fid_e_del_phi_neg(flags_->Flags::Run(), sector, p_, theta_));
		pass &= theta_>=fid_e_theta_cut(flags_->Flags::Run(),sector,p_);
	}
	return pass;
}