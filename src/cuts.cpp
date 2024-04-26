#include "cuts.hpp"

float cuts::fid_e_theta_cut(int run_, int sector_, float p_, int sim_idx_){ 
	return _c1e_[run_][sim_idx_][sector_-1] + _c2e_[run_][sim_idx_][sector_-1] / ((float)p_+_p_shift_e_[run_][sim_idx_][sector_-1]);
}
/*
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
}*/

float cuts::fid_e_del(int run_, int sector_, float p_, float theta_, bool sim_, int side_){
	float c4 = 0.0;
	float c3 = 0.0;
	int sim_idx = 0;
	if(sim_){sim_idx = 1;}
	for(int i=0; i<3; i++){
		c4 += _c4e_[run_][sim_idx][sector_-1][side_][i]*pow(p_,i);
		c3 += _c3e_[run_][sim_idx][sector_-1][side_][i]*pow(p_,i);
	}
	return c4*TMath::Power((TMath::Sin((theta_-cuts::fid_e_theta_cut(run_,sector_,p_,sim_))/_degree_)),c3);
}

bool cuts::fid_e (float p_, float theta_, float phi_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	if(flags_->Flags::Fid_Cut(0) && !isnan(phi_) && !isnan(theta_)){
		//Calculate angles and sector
		float phi_c = physics::phi_center(phi_);
		int sector = physics::get_sector(phi_);
		int sim_idx = 0;
		if(flags_->Flags::Sim()){sim_idx = 1;}
		//Application of Cut
		pass = true; 
		//pass &= (phi_c<=cuts::fid_e_del_phi_pos(flags_->Flags::Run(), sector, p_, theta_));
		//pass &= (phi_c >=cuts::fid_e_del_phi_neg(flags_->Flags::Run(), sector, p_, theta_));
		pass &= (phi_c <= cuts::fid_e_del(flags_->Flags::Run(), sector, p_, theta_, flags_->Flags::Sim(), 0));
		pass &= (phi_c >= cuts::fid_e_del(flags_->Flags::Run(), sector, p_, theta_, flags_->Flags::Sim(), 1));
		pass &= theta_>=fid_e_theta_cut(flags_->Flags::Run(),sector,p_,sim_idx);
	}
	return pass;
}

float cuts::phi_min(int hadron_, float theta_, int run_, int sector_, int sim_){
	//return -(_a0mh_[run_][hadron_][sector_-1]*(1.0-TMath::Exp(-_a1mh_[run_][hadron_][sector_-1]*(theta_-_a2mh_[run_][hadron_][sector_-1])))-_a3mh_[run_][hadron_][sector_-1]);
	return (_a0mh_[run_][sim_][hadron_][sector_-1]*(1.0-TMath::Exp(-_a1mh_[run_][sim_][hadron_][sector_-1]*(theta_-_a2mh_[run_][sim_][hadron_][sector_-1])))+_a3mh_[run_][sim_][hadron_][sector_-1]);
}

float cuts::phi_max(int hadron_, float theta_, int run_, int sector_, int sim_){
	return (_a0xh_[run_][sim_][hadron_][sector_-1]*(1.0-TMath::Exp(-_a1xh_[run_][sim_][hadron_][sector_-1]*(theta_-_a2xh_[run_][sim_][hadron_][sector_-1])))+_a3xh_[run_][sim_][hadron_][sector_-1]);
}

bool cuts::fid_h (int hadron_, float theta_, float phi_, std::shared_ptr<Flags> flags_){ //{0,1}->{pro,pip}
	bool pass = false;
	//Calculate angles and sector
	//Phi is centered for the sector
	float phi_c = physics::phi_center(phi_);
	int sector = physics::get_sector(phi_);
	int sim_stat = -1;
	if(flags_->Flags::Sim()){
		sim_stat = 1;
	}else{
		sim_stat = 0;
	}

	//Actual application of the cut
	if(phi_c>=cuts::phi_min(hadron_, theta_, flags_->Flags::Run(), sector, sim_stat) && phi_c<=cuts::phi_max(hadron_, theta_, flags_->Flags::Run(), sector, sim_stat) && !isnan(phi_) && !isnan(theta_))
	{
		pass = (theta_ >= _hadron_min_theta_[flags_->Flags::Run()][sim_stat][hadron_][sector-1]);
	}
	return pass;
}
float cuts::fid_h_theta_cut(int run_, int sector_, float p_){ 
	return _c1h_[run_][sector_-1] + _c2h_[run_][sector_-1] / ((float)p_+_p_shift_h_[run_][sector_-1]);
}
float cuts::fid_h_expon_pos(int run_, int sector_, float p_){ 
	return _c3h_[run_][sector_-1][0] * TMath::Power((float)p_,_factor_h_[run_][sector_-1][0]);
}
float cuts::fid_h_expon_neg(int run_, int sector_, float p_){ 
	return _c3h_[run_][sector_-1][1] * TMath::Power((float)p_,_factor_h_[run_][sector_-1][1]);
}
float cuts::fid_h_del_phi_pos(int run_, int sector_, float p_, float theta_){ 
	return _c4h_[run_][sector_-1][0] * TMath::Power((TMath::Sin((theta_-cuts::fid_h_theta_cut(run_,sector_,p_))/_degree_)),cuts::fid_h_expon_pos(run_,sector_,p_));
}
float cuts::fid_h_del_phi_neg(int run_, int sector_, float p_, float theta_){ 
	return -_c4h_[run_][sector_-1][1] * TMath::Power((TMath::Sin((theta_-cuts::fid_h_theta_cut(run_,sector_,p_))/_degree_)),cuts::fid_h_expon_neg(run_,sector_,p_));
}

bool cuts::fid_pim (float p_, float theta_, float phi_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	if(flags_->Flags::Fid_Cut(3)){
		//Calculate angles and sector
		float phi_c = physics::phi_center(phi_);
		int sector = physics::get_sector(phi_);

		//Application of Cut
		pass = true; 
		pass &= (phi_c<=cuts::fid_h_del_phi_pos(flags_->Flags::Run(), sector, p_, theta_));
		pass &= (phi_c >=cuts::fid_h_del_phi_neg(flags_->Flags::Run(), sector, p_, theta_));
		pass &= theta_>=fid_h_theta_cut(flags_->Flags::Run(),sector,p_);
	}
	return pass;
}

bool cuts::fid_cut(int part_, float p_, float theta_, float phi_, std::shared_ptr<Flags> flags_){
	bool pass = false; 
	switch(part_){
		case 0:
			pass = cuts::fid_e(p_,theta_,phi_,flags_);
		break;
		case 1:
			pass = cuts::fid_h(part_-1,theta_,phi_,flags_);
		break;
		case 2:
			pass = cuts::fid_h(part_-1,theta_,phi_,flags_);
		break;
		case 3:
			pass = cuts::fid_pim(p_,theta_,phi_,flags_);
		break;
		default:
			std::cout<<"Particle Index Error in fid cut: " <<part_ <<"\n";
		break;
	}
	//if(part ==0){//Electron
	//	pass = fid_e(p,cx,cy,cz);
	//}else{
	//	pass = fid_h(p,cx,cy,cz);
	//}
	return pass;
}

//Proton Formulae from Arjun
float cuts::delta_low(int part_, int run_, bool sim_, float p_){
	float d0 = NAN;
	float d1 = NAN;
	float d2 = NAN;
	float d3 = NAN;
	int sim_idx = 0;
	if(sim_){
		sim_idx = 1;
	}	
	switch(part_){
		/*case 0:
			d0=_dt_ele_low_[run_][sim_idx][0];
			d1=_dt_ele_low_[run_][sim_idx][1];
			d2=_dt_ele_low_[run_][sim_idx][2];
			d3=_dt_ele_low_[run_][sim_idx][3];
		break;*/
		case 1:
			d0=_dt_pro_low_[run_][sim_idx][0];
			d1=_dt_pro_low_[run_][sim_idx][1];
			d2=_dt_pro_low_[run_][sim_idx][2];
			d3=_dt_pro_low_[run_][sim_idx][3];
		break;
		case 2:
			d0=_dt_pip_low_[run_][sim_idx][0];
			d1=_dt_pip_low_[run_][sim_idx][1];
			d2=_dt_pip_low_[run_][sim_idx][2];
			d3=_dt_pip_low_[run_][sim_idx][3];
		break;
		case 3:
			d0=_dt_pim_low_[run_][sim_idx][0];
			d1=_dt_pim_low_[run_][sim_idx][1];
			d2=_dt_pim_low_[run_][sim_idx][2];
			d3=_dt_pim_low_[run_][sim_idx][3];
		break;
		default:
			std::cout<<"Particle Index Error in Delta Low cut: " <<part_ <<"\n";
		break;
	}
	return d0+d1*p_+d2*p_*p_+d3*p_*p_*p_;
}

float cuts::delta_high(int part_, int run_, bool sim_, float p_){
	float d0 = NAN;
	float d1 = NAN;
	float d2 = NAN;
	float d3 = NAN;
	int sim_idx = 0;
	if(sim_){
		sim_idx = 1;
	}	
	switch(part_){
		/*case 0:
			d0=_dt_ele_top_[run_][sim_idx][0];
			d1=_dt_ele_top_[run_][sim_idx][1];
			d2=_dt_ele_top_[run_][sim_idx][2];
			d3=_dt_ele_top_[run_][sim_idx][3];
		break;*/
		case 1:
			d0=_dt_pro_top_[run_][sim_idx][0];
			d1=_dt_pro_top_[run_][sim_idx][1];
			d2=_dt_pro_top_[run_][sim_idx][2];
			d3=_dt_pro_top_[run_][sim_idx][3];
		break;
		case 2:
			d0=_dt_pip_top_[run_][sim_idx][0];
			d1=_dt_pip_top_[run_][sim_idx][1];
			d2=_dt_pip_top_[run_][sim_idx][2];
			d3=_dt_pip_top_[run_][sim_idx][3];
		break;
		case 3:
			d0=_dt_pim_top_[run_][sim_idx][0];
			d1=_dt_pim_top_[run_][sim_idx][1];
			d2=_dt_pim_top_[run_][sim_idx][2];
			d3=_dt_pim_top_[run_][sim_idx][3];
		break;
		default:
			std::cout<<"Particle Index Error in Delta High cut: " <<part_ <<"\n";
		break;
	}
	return d0+d1*p_+d2*p_*p_+d3*p_*p_*p_;
}

bool cuts::delta_t_cut(int species_, float p_, float dt_, std::shared_ptr<Flags> flags_){
	//std::cout<<"\tPerforming Delta Cut for: " <<_species_[species_] <<"\n";
	bool pass = false;
	if(flags_->Flags::Delta_Cut(species_)){
		if(dt_>cuts::delta_low(species_,flags_->Flags::Run(),flags_->Flags::Sim(),p_) && dt_<cuts::delta_high(species_,flags_->Flags::Run(),flags_->Flags::Sim(),p_)){
			pass = true;
		}
	}
	return pass;
}

bool cuts::min_cc(int cc_segm_, int cc_sect_, int nphe_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	if(flags_->Flags::CC_Cut()){
		int seg = detect::cc_segment(cc_segm_);
		int side_2 = -1;
		//if(sector == -1){
		//	if(nphe > 35){//Just a flat cut for now
		//		pass = true;
		//	}
		//}
		if(cc_sect_ >= 0 && cc_segm_ >= 0 && nphe_ > 0){
			int side = detect::cc_lrc(cc_segm_);
			switch(side){//Due to having a different layout for the layers in my cut determination program, the indexes need to be changed to reflect that 1/19/21
				case 0:
					side_2 = 1;
				break;
				case 1:
					side_2 = 2;
				break;
				case 2:
					side_2 = 0;
				break;
				default:
					side_2 = -1;
				break;
			}
			//std::cout<<std::endl <<"CC Sec:" <<cc_sect <<"|| cc seg:" <<seg <<"|| side:" <<side <<"->" <<side_2 <<"|| nphe:" <<nphe <<"|| MinCC: " <<MinCC_Cut_New[cc_sect-1][side_2][seg] <<std::endl;
			if(side_2 != -1){
				if(nphe_ > _cc_min_[flags_->Flags::Run()][cc_sect_-1][side_2][seg]){
					pass = true;
				}
			}
		}
	}
	return pass; 
}


bool cuts::min_ec(Float_t etot_, int sector_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	int sim_idx = 0;
	if(flags_->Flags::Sim()){
		sim_idx = 1; 
	}
	if(etot_ >= _ec_min_cut_[flags_->Flags::Run()][sim_idx][sector_-1]){
		pass = true;
	}
	return pass;
}

float cuts::sf_low(int run_, int sim_, Float_t p_, int sector_){
	//std::cout<<"sf_low: "<<_sf_low_e16_[run_][sim_][sector_-1][0] + _sf_low_e16_[run_][sim_][sector_-1][1]*p_ + _sf_low_e16_[run_][sim_][sector_-1][2]*p_*p_ + _sf_low_e16_[run_][sim_][sector_-1][3]*p_*p_*p_ <<"\n";
	return _sf_low_e16_[run_][sim_][sector_-1][0] + _sf_low_e16_[run_][sim_][sector_-1][1]*p_ + _sf_low_e16_[run_][sim_][sector_-1][2]*p_*p_ + _sf_low_e16_[run_][sim_][sector_-1][3]*p_*p_*p_;
}

float cuts::sf_top(int run_, int sim_, Float_t p_, int sector_){
	//std::cout<<"sf_top: "<<_sf_top_e16_[run_][sim_][sector_-1][0] + _sf_top_e16_[run_][sim_][sector_-1][1]*p_ + _sf_top_e16_[run_][sim_][sector_-1][2]*p_*p_ + _sf_top_e16_[run_][sim_][sector_-1][3]*p_*p_*p_ <<"\n";
	return _sf_top_e16_[run_][sim_][sector_-1][0] + _sf_top_e16_[run_][sim_][sector_-1][1]*p_ + _sf_top_e16_[run_][sim_][sector_-1][2]*p_*p_ + _sf_top_e16_[run_][sim_][sector_-1][3]*p_*p_*p_;
}
//Electron Sampling Fraction cut using the EC
bool cuts::sf_cut(float p_, float sf_, float phi_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	int sim_idx = 0; 
	if(flags_->Flags::Sim()){
		sim_idx = 1; 
	}
	int sector = physics::get_sector(phi_);
	//std::cout<<"sf: " <<sf_ <<"\n";
	if(sf_ >= cuts::sf_low(flags_->Flags::Run(),sim_idx,p_,sector) && sf_ <= cuts::sf_top(flags_->Flags::Run(),sim_idx,p_,sector)){
		pass = true;
	}
	return pass;
}
//Not utilized yet, but will be in time 
float cuts::mm_top(int run_, int sim_, int top_, int sector_, float W_, int cut_width_){
	float output = 0.0;
	/*
	run {0,1} -> e16,e1f
	sim {0,1} -> sim,exp
	top {0,1,2,3} -> {mpro,mpip,mpim,mzero}
	*/
	//std::cout<<"\tgetting mm_top for W: " <<W_ <<"\n";
	//std::cout<<"\t\tsim : " <<sim_ <<" run: " <<run_ <<" top: " <<top_ <<" sector: " <<sector_-1 <<" top/bot: " <<1 <<" slope/inter: " <<0 <<"\n";
	//std::cout<<"\t\t\tslope: " <<_mm2_par_[sim_][run_][top_][sector_-1][1][0] <<"\tintercept: " <<_mm2_par_[sim_][run_][top_][sector_-1][1][1] <<"\n";
	//std::cout<<"\t\tupper: " <<_mm2_par_[sim_][run_][top_][sector_-1][1][0]*W_ + 						_mm2_par_[sim_][run_][top_][sector_-1][1][1] <<"\n";
	//std::cout<<"\t\t\tCheck: " <<_mm2_par_[sim_][run_][top_][sector_-1][1][0]*W_ <<" + " <<	_mm2_par_[sim_][run_][top_][sector_-1][1][1] <<"\n";
	
	return ((_mm2_var_[run_][sim_][top_][sector_-1][1][cut_width_][0]*W_) + _mm2_var_[run_][sim_][top_][sector_-1][1][cut_width_][1]);
	//return ((_mm2_par_[sim_][run_][top_][sector_-1][1][0]*W_) + _mm2_par_[sim_][run_][top_][sector_-1][1][1]);
}
//Not utilized yet, but will be in time 
float cuts::mm_bot(int run_, int sim_, int top_, int sector_, float W_, int cut_width_){
	/*
	run {0,1} -> e16,e1f
	sim {0,1} -> sim,exp
	top {0,1,2,3} -> {mpro,mpip,mpim,mzero}
	*/
	//std::cout<<"\tgetting mm_bot: for W: " <<W_ <<"\n";
	//std::cout<<"\t\tsim : " <<sim_ <<" run: " <<run_ <<" top: " <<top_ <<" sector: " <<sector_-1 <<" top/bot: " <<0 <<" slope/inter: " <<0 <<"\n";
	//std::cout<<"\t\t\tslope" <<_mm2_par_[sim_][run_][top_][sector_-1][0][0] <<"\tintercept" <<_mm2_par_[sim_][run_][top_][sector_-1][0][1] <<"\n";
	//std::cout<<"\t\tlower: " <<_mm2_par_[sim_][run_][top_][sector_-1][0][0]*W_ + _mm2_par_[sim_][run_][top_][sector_-1][0][1] <<"\n";
	//std::cout<<"\t\t\tCheck: " <<_mm2_par_[sim_][run_][top_][sector_-1][0][0]*W_ <<" + " <<_mm2_par_[sim_][run_][top_][sector_-1][0][1] <<"\n";
	
	return ((_mm2_var_[run_][sim_][top_][sector_-1][0][cut_width_][0]*W_) + _mm2_var_[run_][sim_][top_][sector_-1][0][cut_width_][1]);
	//return ((_mm2_par_[sim_][run_][top_][sector_-1][0][0]*W_) + _mm2_par_[sim_][run_][top_][sector_-1][0][1]);
}
//Will need some modification to have W dependence, but this will work for now
bool cuts::MM_cut(int top_, float mm_, int sector_, float W_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	if(!flags_->Flags::MM_Cut(top_)){
		return false;
	}
	int sim_idx = 0; 
	if(flags_->Flags::Sim()){
		sim_idx = 1; 
	}
	//std::cout<<"Trying to cut MM for: " <<_top_[top_] <<" with W: " <<W_ <<"\n";
	//std::cout<<"\tMissing Mass: " <<mm_ <<" with top: " <<cuts::mm_top(flags_->Flags::Run(),sim_idx,top_,sector_,W_) <<" and bot: " <<cuts::mm_bot(flags_->Flags::Run(),sim_idx,top_,sector_,W_) <<"\n";

	if(mm_ < cuts::mm_top(flags_->Flags::Run(),sim_idx,top_,sector_,W_) && mm_ > cuts::mm_bot(flags_->Flags::Run(),sim_idx,top_,sector_,W_)){
		pass = true;
	}
	return pass; 


	/*//New comment out
	int sim_idx = 0; 
	if(flags_->Flags::Sim()){
		sim_idx = 1; 
	}
	//
	//Current non-W dependent method
	if(mm_ > _mm_bounds_[flags_->Flags::Run()][sim_idx][top_][0] && mm_ < _mm_bounds_[flags_->Flags::Run()][sim_idx][top_][1]){
		pass = true;
	}
	if(pass){
		//std::cout<<"Top: " <<_top_[top_] <<" MM: " <<mm_ <<" bounds: " <<_mm_bounds_[flags_->Flags::Run()][sim_idx][top_][0] <<" - " <<_mm_bounds_[flags_->Flags::Run()][sim_idx][top_][1] <<"  Pass:" <<pass <<"\n";
	}

	//Future W Dependent Method
	/*if(mm_ > cuts::mm_bot(flags_->Flags::Run(),sim_idx,top_,W_) && mm_ < cuts::mm_top(flags_->Flags::Run(),sim_idx,top_,W_)){
		pass = true;
	}*/
	//return pass; 
}

//Putting together the e_sanity cuts based on environment we set
bool cuts::in_range(float W_, float Q2_){
	bool pass = false;
	if(W_ >= _WminAna_ && W_ <= _WmaxAna_ ){//Checking to see if the particle is in the relevant W Q2 region 
		if(Q2_ >= _Q2minAna_ && Q2_ <= _Q2maxAna_){
			pass = true;
		}
	}
	return pass; 
}

bool cuts::in_range(float W_){
	bool pass = false;
	if(W_ >= _WminAna_ && W_ <= _WmaxAna_ ){//Checking to see if the particle is in the relevant W Q2 region 
		pass = true;
	}
	return pass; 
}

bool cuts::e_sanity(int dc_, int sc_, int ec_, int cc_, int stat_){
  bool pass = (dc_>0);
  pass &= (sc_>0);
  pass &= (ec_>0);
  pass &= (cc_>0);
  pass &= (stat_>0);
  return pass; 
}

bool cuts::pro_sanity(int dc_, int sc_, int stat_){
	bool pass = (dc_>0);
  pass &= (sc_>0);
  pass &= (stat_>0);
  return pass;
}
bool cuts::pip_sanity(int dc_, int sc_, int stat_){
	bool pass = (dc_>0);
  pass &= (sc_>0);
  pass &= (stat_>0);
  return pass;
}
bool cuts::pim_sanity(int dc_, int sc_, int stat_){
	bool pass = (dc_>0);
  pass &= (sc_>0);
  pass &= (stat_>0);
  return pass;
}

bool cuts::h_sanity(char* species_, int q_, int dc_, int sc_, int stat_){
  bool pass = (dc_>0);
  pass &= (sc_>0);
  pass &= (stat_>0);
  if(species_ == _pim_){
  	pass &= (q_ < 0.0);
  }else if(species_ == _pro_ || species_ == _pip_){
  	pass &= (q_ > 0.0);
  }
  return pass; 
}

bool cuts::vertex_cut(float vz_, int run_){
	bool pass = false;
	if(vz_ > _vz_bot_[run_] && vz_ < _vz_top_[run_]){
		pass = true;
	}
	return pass;
}
//These all need to be built up 1/18/23
bool cuts::sc_eff_ele_cut(float p_, float theta_, int run_){
	bool pass = false;
	return pass;
}
bool cuts::sc_eff_pro_cut(float p_, float theta_, int run_){
	bool pass = false;
	return pass;
}
bool cuts::sc_eff_pip_cut(float p_, float theta_, int run_){
	bool pass = false;
	return pass;
}
bool cuts::sc_eff_pim_cut(float p_, float theta_, int run_){
	bool pass = false;
	return pass;

}
bool cuts::sc_eff_cut(float p_, float theta_, int run_, int par_){
	bool pass = false;
	switch(par_){
		case 0:
			pass = sc_eff_ele_cut(p_,theta_,run_);
		break;
		case 1:
			pass = sc_eff_pro_cut(p_,theta_,run_);
		break;
		case 2:
			pass = sc_eff_pim_cut(p_,theta_,run_);
		break;
		case 3:
			pass = sc_eff_pim_cut(p_,theta_,run_);
		break;
		default:
			pass = sc_eff_ele_cut(p_,theta_,run_);
		break;
	}
	return pass;
}

