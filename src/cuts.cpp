#include "cuts.hpp"

float cuts::fid_e_theta_cut(int run_, int sector_, float p_, int sim_idx_, int cut_width_){ 
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

float cuts::fid_e_del(int run_, int sector_, float p_, float theta_, bool sim_, int side_, int cut_width_){
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

float cuts::phi_min(int hadron_, float theta_, int run_, int sector_, int sim_, int cut_width_){
	//return -(_a0mh_[run_][hadron_][sector_-1]*(1.0-TMath::Exp(-_a1mh_[run_][hadron_][sector_-1]*(theta_-_a2mh_[run_][hadron_][sector_-1])))-_a3mh_[run_][hadron_][sector_-1]);
	return (_a0mh_[run_][sim_][hadron_][sector_-1]*(1.0-TMath::Exp(-_a1mh_[run_][sim_][hadron_][sector_-1]*(theta_-_a2mh_[run_][sim_][hadron_][sector_-1])))+_a3mh_[run_][sim_][hadron_][sector_-1]);
}

float cuts::phi_max(int hadron_, float theta_, int run_, int sector_, int sim_, int cut_width_){
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
float cuts::fid_h_theta_cut(int run_, int sector_, float p_, int cut_width_){ 
	return _c1h_[run_][sector_-1] + _c2h_[run_][sector_-1] / ((float)p_+_p_shift_h_[run_][sector_-1]);
}
float cuts::fid_h_expon_pos(int run_, int sector_, float p_, int cut_width_){ 
	return _c3h_[run_][sector_-1][0] * TMath::Power((float)p_,_factor_h_[run_][sector_-1][0]);
}
float cuts::fid_h_expon_neg(int run_, int sector_, float p_, int cut_width_){ 
	return _c3h_[run_][sector_-1][1] * TMath::Power((float)p_,_factor_h_[run_][sector_-1][1]);
}
float cuts::fid_h_del_phi_pos(int run_, int sector_, float p_, float theta_, int cut_width_){ 
	return _c4h_[run_][sector_-1][0] * TMath::Power((TMath::Sin((theta_-cuts::fid_h_theta_cut(run_,sector_,p_))/_degree_)),cuts::fid_h_expon_pos(run_,sector_,p_));
}
float cuts::fid_h_del_phi_neg(int run_, int sector_, float p_, float theta_, int cut_width_){ 
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
float cuts::delta_low(int part_, int run_, bool sim_, float p_, int sec_, int cut_width_){
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
			d0=_dt_pro_low_[run_][sim_idx][sec_-1][cut_width_][0];
			d1=_dt_pro_low_[run_][sim_idx][sec_-1][cut_width_][1];
			d2=_dt_pro_low_[run_][sim_idx][sec_-1][cut_width_][2];
			d3=_dt_pro_low_[run_][sim_idx][sec_-1][cut_width_][3];
		break;
		case 2:
			d0=_dt_pip_low_[run_][sim_idx][sec_-1][cut_width_][0];
			d1=_dt_pip_low_[run_][sim_idx][sec_-1][cut_width_][1];
			d2=_dt_pip_low_[run_][sim_idx][sec_-1][cut_width_][2];
			d3=_dt_pip_low_[run_][sim_idx][sec_-1][cut_width_][3];
		break;
		case 3:
			d0=_dt_pim_low_[run_][sim_idx][sec_-1][cut_width_][0];
			d1=_dt_pim_low_[run_][sim_idx][sec_-1][cut_width_][1];
			d2=_dt_pim_low_[run_][sim_idx][sec_-1][cut_width_][2];
			d3=_dt_pim_low_[run_][sim_idx][sec_-1][cut_width_][3];
		break;
		default:
			std::cout<<"Particle Index Error in Delta Low cut: " <<part_ <<"\n";
		break;
	}
	return d0+d1*p_+d2*p_*p_+d3*p_*p_*p_;
}

float cuts::delta_high(int part_, int run_, bool sim_, float p_,int sec_, int cut_width_){
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
			d0=_dt_pro_top_[run_][sim_idx][sec_-1][cut_width_][0];
			d1=_dt_pro_top_[run_][sim_idx][sec_-1][cut_width_][1];
			d2=_dt_pro_top_[run_][sim_idx][sec_-1][cut_width_][2];
			d3=_dt_pro_top_[run_][sim_idx][sec_-1][cut_width_][3];
		break;
		case 2:
			d0=_dt_pip_top_[run_][sim_idx][sec_-1][cut_width_][0];
			d1=_dt_pip_top_[run_][sim_idx][sec_-1][cut_width_][1];
			d2=_dt_pip_top_[run_][sim_idx][sec_-1][cut_width_][2];
			d3=_dt_pip_top_[run_][sim_idx][sec_-1][cut_width_][3];
		break;
		case 3:
			d0=_dt_pim_top_[run_][sim_idx][sec_-1][cut_width_][0];
			d1=_dt_pim_top_[run_][sim_idx][sec_-1][cut_width_][1];
			d2=_dt_pim_top_[run_][sim_idx][sec_-1][cut_width_][2];
			d3=_dt_pim_top_[run_][sim_idx][sec_-1][cut_width_][3];
		break;
		default:
			std::cout<<"Particle Index Error in Delta High cut: " <<part_ <<"\n";
		break;
	}
	return d0+d1*p_+d2*p_*p_+d3*p_*p_*p_;
}

bool cuts::delta_t_cut(int species_, float p_, float dt_, int sec_, std::shared_ptr<Flags> flags_){
	//std::cout<<"\tPerforming Delta Cut for: " <<_species_[species_] <<"\n";
	bool pass = false;
	int cut_width = flags_->Flags::Delta_Cut_Width(species_);
	if(flags_->Flags::Delta_Cut(species_)){
		//std::cout<<"delta t " <<species_ <<"  dt:" <<dt_ <<"  low:" <<cuts::delta_low(species_,flags_->Flags::Run(),flags_->Flags::Sim(),p_,sec_,cut_width) <<"  top:" <<cuts::delta_high(species_,flags_->Flags::Run(),flags_->Flags::Sim(),p_,sec_,cut_width) <<"\n";
		if(dt_>cuts::delta_low(species_,flags_->Flags::Run(),flags_->Flags::Sim(),p_,sec_,cut_width) && dt_<cuts::delta_high(species_,flags_->Flags::Run(),flags_->Flags::Sim(),p_,sec_,cut_width)){
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
		int width = flags_->Flags::Min_CC_Cut_Width();
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
				if(nphe_ > _cc_min_[flags_->Flags::Run()][cc_sect_-1][side_2][seg][width]){
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

float cuts::sf_low(int run_, int sim_, Float_t p_, int sector_, int cut_width_){
	//std::cout<<"sf_low: "<<_sf_low_[run_][sim_][sector_-1][0] + _sf_low_e16_[run_][sim_][sector_-1][1]*p_ + _sf_low_e16_[run_][sim_][sector_-1][2]*p_*p_ + _sf_low_e16_[run_][sim_][sector_-1][3]*p_*p_*p_ <<"\n";
	//return _sf_low_e16_[run_][sim_][sector_-1][0] + _sf_low_e16_[run_][sim_][sector_-1][1]*p_ + _sf_low_e16_[run_][sim_][sector_-1][2]*p_*p_ + _sf_low_e16_[run_][sim_][sector_-1][3]*p_*p_*p_;
	return _sf_low_[run_][sim_][sector_-1][cut_width_][0] + _sf_low_[run_][sim_][sector_-1][cut_width_][1]*p_ + _sf_low_[run_][sim_][sector_-1][cut_width_][2]*p_*p_ + _sf_low_[run_][sim_][sector_-1][cut_width_][3]*p_*p_*p_;
}

float cuts::sf_top(int run_, int sim_, Float_t p_, int sector_, int cut_width_){
	//std::cout<<"sf_top: "<<_sf_top_e16_[run_][sim_][sector_-1][0] + _sf_top_e16_[run_][sim_][sector_-1][1]*p_ + _sf_top_e16_[run_][sim_][sector_-1][2]*p_*p_ + _sf_top_e16_[run_][sim_][sector_-1][3]*p_*p_*p_ <<"\n";
	//return _sf_top_e16_[run_][sim_][sector_-1][0] + _sf_top_e16_[run_][sim_][sector_-1][1]*p_ + _sf_top_e16_[run_][sim_][sector_-1][2]*p_*p_ + _sf_top_e16_[run_][sim_][sector_-1][3]*p_*p_*p_;
	return _sf_top_[run_][sim_][sector_-1][cut_width_][0] + _sf_top_[run_][sim_][sector_-1][cut_width_][1]*p_ + _sf_top_[run_][sim_][sector_-1][cut_width_][2]*p_*p_ + _sf_top_[run_][sim_][sector_-1][cut_width_][3]*p_*p_*p_;
}
//Electron Sampling Fraction cut using the EC
bool cuts::sf_cut(float p_, float sf_, float phi_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	int sim_idx = 0; 
	if(flags_->Flags::Sim()){
		sim_idx = 1; 
	}
	int cut_width = fun::cut_width(_ele_,_sf_cut_,flags_);
	//std::cout<<"cut_width:" <<cut_width <<"\n";
	int sector = physics::get_sector(phi_);
	//std::cout<<"sf: " <<sf_ <<"   low:" <<cuts::sf_low(flags_->Flags::Run(),sim_idx,p_,sector,cut_width) <<"  top:" <<cuts::sf_top(flags_->Flags::Run(),sim_idx,p_,sector,cut_width) <<"\n";
	if(sf_ >= cuts::sf_low(flags_->Flags::Run(),sim_idx,p_,sector,cut_width) && sf_ <= cuts::sf_top(flags_->Flags::Run(),sim_idx,p_,sector,cut_width)){
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
	int cut_width = flags_->Flags::MM_Cut_Width(top_);

	if(mm_ < cuts::mm_top(flags_->Flags::Run(),sim_idx,top_,sector_,W_) && mm_ >= cuts::mm_bot(flags_->Flags::Run(),sim_idx,top_,sector_,W_)){
		pass = true;
	}
	/*if(pass){
		std::cout<<"\tPassed MM cut for: " <<_top_[top_] <<" with W: " <<W_ <<" pass:" <<pass <<"\n";
		std::cout<<"\t\tMissing Mass: " <<mm_ <<" with top: " <<cuts::mm_top(flags_->Flags::Run(),sim_idx,top_,sector_,W_) <<" and bot: " <<cuts::mm_bot(flags_->Flags::Run(),sim_idx,top_,sector_,W_) <<"\n";
	}else{
		std::cout<<"\tFailed MM cut for: " <<_top_[top_] <<" with W: " <<W_ <<" pass:" <<pass <<"\n";
		std::cout<<"\t\tMissing Mass: " <<mm_ <<" with top: " <<cuts::mm_top(flags_->Flags::Run(),sim_idx,top_,sector_,W_) <<" and bot: " <<cuts::mm_bot(flags_->Flags::Run(),sim_idx,top_,sector_,W_) <<"\n";
	}*/

	
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

bool cuts::vertex_cut(float vz_, int run_, int sector_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	int sim_idx = 0;
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	int cut_width = fun::cut_width(_ele_,_vertex_cut_,flags_);
	//std::cout<<"vertex ele:" <<vz_ <<"  low:" <<_vz_cut_par_[run_][sim_idx][sector_-1][cut_width][0] <<"  top:" <<_vz_cut_par_[run_][sim_idx][sector_-1][cut_width][1] <<"\n";
	if(vz_ >= _vz_cut_par_[run_][sim_idx][sector_-1][cut_width][0] && vz_ < _vz_cut_par_[run_][sim_idx][sector_-1][cut_width][1]){
		pass = true;
	}
	return pass;
}
//These all need to be built up 1/18/23
bool cuts::sc_eff_ele_cut(float p_, float theta_, int run_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	return pass;
}
bool cuts::sc_eff_pro_cut(float p_, float theta_, int run_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	return pass;
}
bool cuts::sc_eff_pip_cut(float p_, float theta_, int run_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	return pass;
}
bool cuts::sc_eff_pim_cut(float p_, float theta_, int run_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	return pass;

}
bool cuts::sc_eff_cut(float p_, float theta_, int run_, int par_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	switch(par_){
		case 0:
			pass = sc_eff_ele_cut(p_,theta_,run_,flags_);
		break;
		case 1:
			pass = sc_eff_pro_cut(p_,theta_,run_,flags_);
		break;
		case 2:
			pass = sc_eff_pim_cut(p_,theta_,run_,flags_);
		break;
		case 3:
			pass = sc_eff_pim_cut(p_,theta_,run_,flags_);
		break;
		default:
			pass = sc_eff_ele_cut(p_,theta_,run_,flags_);
		break;
	}
	return pass;
}


//Geometric CC Cuts
//edge Cuts for Right and Left {e16,e1f}{exp,sim}{sec}{tight,mid,loose}{right,left}{intercept,slope}
//static float _geo_cc_side_edge_pars_[2][2][6][3][2][2]
//Center Cuts for Right and Left {e16,e1f}{exp,sim}{sec}{tight,mid,loose}{right,left}
//static float _geo_cc_side_center_pars_[2][2][6][3][2]
//static float _geo_cc_mid_x_pars_[2][2][6][_num_cc_segments_][3][2]
//{run group}{exp/sim}{seg}{sec}{side}{tight,mid,loose}{bot,top}
//static float _geo_seg_bounds_[2][2][_num_cc_segments_][6][3][3][2]
bool cuts::cc_geo_left_cut(int par_, float x_, float y_,int sec_,int side_, int seg_, std::shared_ptr<Flags> flags_){
	bool pass = true;
	float x_cen = detect::x_det_center(x_,y_,sec_+1);
	float y_cen = detect::y_det_center(x_,y_,sec_+1);
	int side_idx = side_;
	int cut_width = flags_->Flags::Fid_Geo_Cut_Width(par_,0);
	int run = flags_->Flags::Run();
	int sim_idx = 0;
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	if(side_ == -1){return false;}
	if(seg_ == -1){return false;}

	pass &= (y_cen > (_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][0]+_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][1]*x_cen));
	//std::cout<<"CC_Geo_Left_Cut\n";
	//std::cout<<"par:" <<par_ <<" x_cen:" <<x_cen <<" y_cen:" <<y_cen <<" sec:" <<sec_ <<" side:" <<side_ <<" seg:" <<seg_ <<"\n";
	/*
	if( side_ == 0){
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (_geo_cc_side_center_pars_[run][sim_idx][sec_][cut_width][side_idx] !=0.0 );
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (x_cen > _geo_cc_side_center_pars_[run][sim_idx][sec_][cut_width][side_idx]);
	//	std::cout<<"x_cen:" <<x_cen <<" > " <<_geo_cc_side_center_pars_[run][sim_idx][sec_][cut_width][1] <<"\n";
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
	}else if(side_ == 1){
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (_geo_cc_mid_x_pars_[run][sim_idx][sec_][seg_][cut_width][1] != 0.0);
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (x_cen > _geo_cc_mid_x_pars_[run][sim_idx][sec_][seg_][cut_width][0]);
	//	std::cout<<"x_cen:" <<x_cen <<" > " <<_geo_cc_side_center_pars_[run][sim_idx][sec_][cut_width][1] <<"\n";
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
	}else if(side_ == 2){
		side_idx = 1; 
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][0] != 0.0);
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][1] != 0.0);
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (y_cen > (_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][0]+_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][1]*x_cen));
	//	std::cout<<"y_cen:" <<y_cen <<" > " <<(_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][0]+_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][1][0]*x_cen) <<"\n";
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
	}
	*/
	return pass;
}
bool cuts::cc_geo_right_cut(int par_, float x_, float y_,int sec_,int side_, int seg_, std::shared_ptr<Flags> flags_){
	bool pass = true;
	float x_cen = detect::x_det_center(x_,y_,sec_+1);
	float y_cen = detect::y_det_center(x_,y_,sec_+1);
	int side_idx = side_;
	int cut_width = flags_->Flags::Fid_Geo_Cut_Width(par_,0);
	int run = flags_->Flags::Run();
	int sim_idx = 0;
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	if(side_ == -1){return false;}
	if(seg_ == -1){return false;}
	//std::cout<<"CC_Geo_Right_Cut\n";
	//std::cout<<"par:" <<par_ <<" x_cen:" <<x_cen <<" y_cen:" <<y_cen <<" sec:" <<sec_ <<" side:" <<side_ <<" seg:" <<seg_ <<"\n";
	pass &= (y_cen > (_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][0]+_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][1]*x_cen));
	/*
	if( side_ == 2){
		side_idx = 1; 
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (_geo_cc_side_center_pars_[run][sim_idx][sec_][cut_width][0] != 0.0);
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (x_cen < _geo_cc_side_center_pars_[run][sim_idx][sec_][cut_width][side_idx]);
	//	std::cout<<"x_cen:" <<x_cen <<" < " <<_geo_cc_side_center_pars_[run][sim_idx][sec_][cut_width][0] <<"\n";
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
	}else if(side_ == 1){
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (_geo_cc_mid_x_pars_[run][sim_idx][sec_][seg_][cut_width][0] != 0.0);
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (x_cen < _geo_cc_mid_x_pars_[run][sim_idx][sec_][seg_][cut_width][1]);
	//	std::cout<<"x_cen:" <<x_cen <<" < " <<_geo_cc_mid_x_pars_[run][sim_idx][sec_][seg_][cut_width][0] <<"\n";
	//	std::cout<<"\t\tpass:" <<pass <<"\n";
	}else if(side_ == 0){
		side_idx = 0;
		//std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][0] != 0);
		//std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][1] !=0);
		//std::cout<<"\t\tpass:" <<pass <<"\n";
		pass &= (y_cen > (_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][0]+_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][1]*x_cen));
		//std::cout<<"ycen:" <<y_cen <<" > " <<(_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][0]+_geo_cc_side_edge_pars_[run][sim_idx][sec_][cut_width][side_idx][0]*x_cen) <<"\n";
		//std::cout<<"\t\tpass:" <<pass <<"\n";
	}
	*/
	return pass;
}
bool cuts::cc_geo_seg_cut(int par_, float x_, float y_,int sec_,int side_, int seg_, std::shared_ptr<Flags> flags_){
	bool pass = true;
	if(seg_ == -1){return pass;}
	float x_cen = detect::x_det_center(x_,y_,sec_+1);
	float y_cen = detect::y_det_center(x_,y_,sec_+1);
	int side_idx = side_;
	int cut_width = flags_->Flags::Fid_Geo_Cut_Width(par_,0);
	int run = flags_->Flags::Run();
	int sim_idx = 0;
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	if(side_ == -1){return false;}
	if(seg_ == -1){return false;}
	//std::cout<<"CC_Geo_Seg_Cut\n";
	//std::cout<<"par:" <<par_ <<" x_cen:" <<x_cen <<" y_cen:" <<y_cen <<" sec:" <<sec_ <<" side:" <<side_ <<" seg:" <<seg_ <<"\n";
	//std::cout<<"\t\tpass:" <<pass <<"\n";
	//pass &= (_geo_seg_bounds_[run][sim_idx][seg_][sec_][side_][cut_width][0] != 0.0);
	pass &= (_cc_geo_seg_bounds_[run][sim_idx][seg_][sec_][cut_width][0] != 0.0);
	//std::cout<<"\t\tpass:" <<pass <<"\n";
	//pass &= (_geo_seg_bounds_[run][sim_idx][seg_][sec_][side_][cut_width][1] != 0.0);
	pass &= (_cc_geo_seg_bounds_[run][sim_idx][seg_][sec_][cut_width][1] != 0.0);
	//std::cout<<"\t\tpass:" <<pass <<"\n";
	//pass &= (y_cen >  _geo_seg_bounds_[run][sim_idx][seg_][sec_][side_][cut_width][0]);
	pass &= (y_cen >  _cc_geo_seg_bounds_[run][sim_idx][seg_][sec_][cut_width][0]);
	//std::cout<<"y_cen:" <<y_cen <<" > " <<_geo_seg_bounds_[run][sim_idx][seg_][sec_][side_][cut_width][0] <<"\n";
	//std::cout<<"\t\tpass:" <<pass <<"\n";
	//pass &= (y_cen <  _geo_seg_bounds_[run][sim_idx][seg_][sec_][side_][cut_width][1]);
	pass &= (y_cen <  _cc_geo_seg_bounds_[run][sim_idx][seg_][sec_][cut_width][1]);
	//std::cout<<"y_cen:" <<y_cen <<" < " <<_geo_seg_bounds_[run][sim_idx][seg_][sec_][side_][cut_width][1] <<"\n";
	//std::cout<<"\t\tpass:" <<pass <<"\n";
	return pass;
}
bool cuts::cc_geo_cut(int par_, float x_, float y_, int sec_,int side_, int seg_, std::shared_ptr<Flags> flags_){
	bool pass = true;
	int sim_idx = 0;
	if(flags_->Flags::Sim()){ sim_idx=1;}
	if(side_ == -1){return false;}
	if(seg_ == -1){return false;}
	pass &= cuts::cc_geo_left_cut(par_,x_,y_,sec_,side_,seg_,flags_);
	pass &= cuts::cc_geo_right_cut(par_,x_,y_,sec_,side_,seg_,flags_);
	pass &=  cuts::cc_geo_seg_cut(par_,x_,y_,sec_,side_,seg_,flags_);
	return pass; 
}
bool cuts::cc_good_segment(int par_, int sec_, int seg_, std::shared_ptr<Flags> flags_){
	return _cc_good_segments[flags_->Flags::Run()][seg_][sec_];
}

bool cuts::sc_good_paddle(int par_, int sec_, int pad_, std::shared_ptr<Flags> flags_){
	return _cc_good_segments[flags_->Flags::Run()][pad_][sec_];
}

bool cuts::sc_geo_left_cut(int par_, float x_, float y_, int sec_,int pad_, std::shared_ptr<Flags> flags_){
	bool pass = true;
	float x_cen = detect::x_det_center(x_,y_,1);
	float y_cen = detect::y_det_center(x_,y_,1);
	int cut_width = flags_->Flags::Fid_Geo_Cut_Width(par_,1);
	int run = flags_->Flags::Run();
	int sim_idx = 0;
	//std::cout<<"\tSc geo left cut\n";
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	//for(int i=0; i<1; i++){
	//	std::cout<<"\t" <<_geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][1][i];
	//}
	//std::cout<<"y_cen=" <<y_cen <<" is it greater than the line at " <<x_cen <<" : " <<(_geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][1][0] + _geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][1][1]*x_cen) <<"\n";
	if(par_ == 0 || par_==3){
		pass &= (y_cen > (_geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][1][0] + _geo_ec_side_pars_[run][sim_idx][par_][sec_][cut_width][1][1]*x_cen));
	}else{
		float side_min = 0.0;
		for(int i=0; i<6; i++){
			side_min += _geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][1][i]*corr::power(x_cen,i);
		}
		pass &= (y_cen > side_min);
	}
	//pass &= (y_cen > (_geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][1][0] + _geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][1][1]*x_cen));
	return pass;
}
bool cuts::sc_geo_right_cut(int par_, float x_, float y_, int sec_,int pad_, std::shared_ptr<Flags> flags_){
	bool pass = true;
	float x_cen = detect::x_det_center(x_,y_,1);
	float y_cen = detect::y_det_center(x_,y_,1);
	int cut_width = flags_->Flags::Fid_Geo_Cut_Width(par_,1);
	int run = flags_->Flags::Run();
	int sim_idx = 0;
	//std::cout<<"\tSc geo right cut\n";
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	//for(int i=0; i<1; i++){
		//std::cout<<"\t" <<_geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][0][i];
	//}
	//std::cout<<"y_cen=" <<y_cen <<" is it greater than the line at " <<x_cen <<" : " <<(_geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][0][0] + _geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][0][1]*x_cen) <<"\n";
	if(par_ == 0 || par_==3){
		pass &= (y_cen > (_geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][0][0] + _geo_ec_side_pars_[run][sim_idx][par_][sec_][cut_width][0][1]*x_cen));
	}else{
		float side_min = 0.0;
		for(int i=0; i<6; i++){
			side_min += _geo_sc_side_pars_[run][sim_idx][par_][sec_][cut_width][0][i]*corr::power(x_cen,i);
		}
		pass &= (y_cen > side_min);
	}
	
	
	return pass;
}
bool cuts::sc_geo_pad_cut(int par_, float x_, float y_, int sec_,int pad_, std::shared_ptr<Flags> flags_){
	bool pass = true;
	float x_cen = detect::x_det_center(x_,y_,1);
	float y_cen = detect::y_det_center(x_,y_,1);
	int cut_width = flags_->Flags::Fid_Geo_Cut_Width(par_,1);
	int run = flags_->Flags::Run();
	int sim_idx = 0;
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	
	if(pad_ == -1){return false;}
	//std::cout<<"CC_Geo_Seg_Cut\n";
	
	//std::cout<<"\t\tpass:" <<pass <<"\n";
	pass &= (_geo_sc_seg_bounds_[run][sim_idx][par_][pad_][sec_][cut_width][0] != 0.0);
	//std::cout<<"\t\tpass:" <<pass <<"\n";
	pass &= (_geo_sc_seg_bounds_[run][sim_idx][par_][pad_][sec_][cut_width][1] != 0.0);
	//std::cout<<"\t\tpass:" <<pass <<"\n";
	pass &= (y_cen >  _geo_sc_seg_bounds_[run][sim_idx][par_][pad_][sec_][cut_width][0]);
	//std::cout<<"y_cen:" <<y_cen <<" > " <<_geo_sc_seg_bounds_[run][sim_idx][par_][pad_][sec_][cut_width][0] <<"\n";
	//std::cout<<"\t\tpass:" <<pass <<"\n";
	pass &= (y_cen <  _geo_sc_seg_bounds_[run][sim_idx][par_][pad_][sec_][cut_width][1]);
	//std::cout<<"y_cen:" <<y_cen <<" < " <<_geo_sc_seg_bounds_[run][sim_idx][par_][pad_][sec_][cut_width][1] <<"\n";
	//std::cout<<"\t\tpass:" <<pass <<"\n";
	//std::cout<<"par:" <<par_ <<" x_cen:" <<x_cen <<" y_cen:" <<y_cen <<" sec:" <<sec_  <<" pad:" <<pad_ <<" cut_width:" <<cut_width <<"\n";
	//std::cout<<"\t" <<_geo_sc_seg_bounds_[run][sim_idx][par_][pad_][sec_][cut_width][0] <<" < " <<y_cen <<" < " <<_geo_sc_seg_bounds_[run][sim_idx][par_][pad_][sec_][cut_width][1] <<"  pass:" <<pass <<"\n";
	return pass;
}
bool cuts::sc_geo_cut(int par_, float x_, float y_, int sec_,int pad_, std::shared_ptr<Flags> flags_){
	if(isnan(x_) && isnan(y_)){ return false;}
	bool pass = true;
	int sim_idx = 0;
	float side_min=0.0;
	int max_idx = 0;
	if(par_==0 || par_==3){
		max_idx = 2;
	}else{
		max_idx = 6;
	}
	if(flags_->Flags::Sim()){ sim_idx=1;}
	//std::cout<<"par:" <<_species_[par_] <<" with x:" <<x_ <<" y:" <<y_ <<" sec:" <<sec_ <<" pad:" <<pad_ <<"***\n";
	//std::cout<<"\tX should become: " <<detect::x_det_center(x_,y_,1) <<"\n";
	//std::cout<<"\tY should become: " <<detect::y_det_center(x_,y_,1) <<"\n";
	pass &= cuts::sc_geo_left_cut(par_,x_,y_,sec_,pad_,flags_);
	pass &= cuts::sc_geo_right_cut(par_,x_,y_,sec_,pad_,flags_);
	pass &=  cuts::sc_geo_pad_cut(par_,x_,y_,sec_,pad_,flags_);
	
	/*if(!pass && !isnan(x_) && !isnan(y_) && pad_>=0){
	//if(_species_[par_]==_pim_){
		std::cout<<"Perf SC_Geo Cut | par:" <<_species_[par_]<<" x:" <<x_ <<" y:" <<y_ <<" sec:" <<sec_ <<" pad:" <<pad_ <<"***\n";
		std::cout<<"\tX should become: " <<detect::x_det_center(x_,y_,1) <<"\n";
		std::cout<<"\tY should become: " <<detect::y_det_center(x_,y_,1) <<"\n";
		for(int i=0; i<max_idx; i++){
			side_min += _geo_sc_side_pars_[flags_->Flags::Run()][sim_idx][par_][sec_][flags_->Flags::Fid_Geo_Cut_Width(par_,1)][0][i]*corr::power(detect::x_det_center(x_,y_,1),i);
		}
		for(int i=0; i<max_idx; i++){
			std::cout<<"\t\tright par " <<i <<":" <<_geo_sc_side_pars_[flags_->Flags::Run()][sim_idx][par_][sec_][flags_->Flags::Fid_Geo_Cut_Width(par_,1)][0][i] <<"\n";
		}
		std::cout<<"\tRight Cut:" <<cuts::sc_geo_right_cut(par_,x_,y_,sec_,pad_,flags_) <<" min_val:" <<side_min <<"\n";
		
		side_min = 0;
		for(int i=0; i<max_idx; i++){
			side_min += _geo_sc_side_pars_[flags_->Flags::Run()][sim_idx][par_][sec_][flags_->Flags::Fid_Geo_Cut_Width(par_,1)][1][i]*corr::power(detect::x_det_center(x_,y_,1),i);
		}
		for(int i=0; i<max_idx; i++){
			std::cout<<"\t\tleft par " <<i <<":" <<_geo_sc_side_pars_[flags_->Flags::Run()][sim_idx][par_][sec_][flags_->Flags::Fid_Geo_Cut_Width(par_,1)][1][i] <<"\n";
		}
		std::cout<<"\tLeft Cut:" <<cuts::sc_geo_left_cut(par_,x_,y_,sec_,pad_,flags_) <<" min_val:" <<side_min <<"\n";
		std::cout<<"\tSeG Bounds:" <<_geo_sc_seg_bounds_[flags_->Flags::Run()][sim_idx][par_][pad_][sec_][flags_->Flags::Fid_Geo_Cut_Width(par_,1)][0] <<" < " <<	detect::y_det_center(x_,y_,1) <<" < " <<_geo_sc_seg_bounds_[flags_->Flags::Run()][sim_idx][par_][pad_][sec_][flags_->Flags::Fid_Geo_Cut_Width(par_,1)][1] <<"\n";
		std::cout<<"\tSeg Cut:" <<cuts::sc_geo_pad_cut(par_,x_,y_,sec_,pad_,flags_) <<"\n"; 
	}*/
	
	return pass; 
}
bool cuts::ec_geo_left_cut(int par_, float x_, float y_, int sec_,std::shared_ptr<Flags> flags_){
	bool pass = true;
	float x_cen = detect::x_det_center(x_,y_,sec_+1);
	float y_cen = detect::y_det_center(x_,y_,sec_+1);
	int cut_width = flags_->Flags::Fid_Geo_Cut_Width(par_,2);
	int run = flags_->Flags::Run();
	int sim_idx = 0;
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	pass &= (y_cen > (_geo_ec_side_pars_[run][sim_idx][par_][sec_][cut_width][1][0] + _geo_ec_side_pars_[run][sim_idx][par_][sec_][cut_width][1][1]*x_cen));
	return pass;
}
bool cuts::ec_geo_right_cut(int par_, float x_, float y_, int sec_,std::shared_ptr<Flags> flags_){
	bool pass = true;
	float x_cen = detect::x_det_center(x_,y_,sec_+1);
	float y_cen = detect::y_det_center(x_,y_,sec_+1);
	int cut_width = flags_->Flags::Fid_Geo_Cut_Width(par_,2);
	int run = flags_->Flags::Run();
	int sim_idx = 0;
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	pass &= (y_cen > (_geo_ec_side_pars_[run][sim_idx][par_][sec_][cut_width][0][0] + _geo_ec_side_pars_[run][sim_idx][par_][sec_][cut_width][0][1]*x_cen));
	return pass;
}
bool cuts::ec_geo_top_cut(int par_, float x_, float y_, int sec_,std::shared_ptr<Flags> flags_){
	bool pass = true;
	float x_cen = detect::x_det_center(x_,y_,sec_+1);
	float y_cen = detect::y_det_center(x_,y_,sec_+1);
	int cut_width = flags_->Flags::Fid_Geo_Cut_Width(par_,2);
	int run = flags_->Flags::Run();
	int sim_idx = 0;
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	pass &= (y_cen < _geo_ec_top_cut_[run][sim_idx][par_][sec_][cut_width]);
	return pass;
}
bool cuts::ec_geo_eff_cut(int par_, float x_, float y_, int sec_,std::shared_ptr<Flags> flags_){
	bool pass = true;
	float x_cen = detect::x_det_center(x_,y_,sec_+1);
	float y_cen = detect::y_det_center(x_,y_,sec_+1);
	int cut_width = flags_->Flags::Fid_Geo_Cut_Width(par_,2);
	int run = flags_->Flags::Run();
	int sim_idx = 0;
	if(flags_->Flags::Sim()){
		sim_idx = 1;
	}
	if(sec_ != 4){return pass;}
	//The diagonal EXP cut
	bool low1 = (y_cen < (_e16_sec5_eff_cut1_[sim_idx][par_][cut_width][0][0] + _e16_sec5_eff_cut1_[sim_idx][par_][cut_width][0][1]*x_cen));
	bool high1 = (y_cen > (_e16_sec5_eff_cut1_[sim_idx][par_][cut_width][1][0] + _e16_sec5_eff_cut1_[sim_idx][par_][cut_width][1][1]*x_cen));
	//Horizontal SIM cut
	bool low2 = (y_cen < _e16_sec5_eff_cut2_[sim_idx][par_][cut_width][0]);
	bool high2 = (y_cen > _e16_sec5_eff_cut2_[sim_idx][par_][cut_width][1]);
	
	pass &= (low1 || high1 );
	pass &= (low2 || high2 );
	return pass;
}
bool cuts::ec_geo_cut(int par_, float x_, float y_, int sec_,std::shared_ptr<Flags> flags_){
	bool pass = true;
	int sim_idx = 0;
	if(flags_->Flags::Sim()){ sim_idx=1;}
	pass &= cuts::ec_geo_left_cut(par_,x_,y_,sec_,flags_);
	pass &= cuts::ec_geo_right_cut(par_,x_,y_,sec_,flags_);
	//pass &=  cuts::ec_geo_top_cut(par_,x_,y_,sec_,flags_);
	//pass &=  cuts::ec_geo_eff_cut(par_,x_,y_,sec_,flags_);
	return pass; 
}

bool cuts::kin_eff_cut(int par_, int sec_, float p_, float theta_, std::shared_ptr<Flags> flags_){
	return true;
}