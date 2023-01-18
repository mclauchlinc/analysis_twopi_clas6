#include "particle.hpp"

//Initialize the particle by filling it with everything we want
Particle::Particle(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_, bool thrown_){
	_idx = idx_; 
	_run = flags_->Flags::Run();
	_thrown = thrown_;
	_sim = flags_->Flags::Sim();
	if(_sim){
		_weight = data_->Branches::weight();
	}else{
		_weight = 1.0;
	}
	if(_thrown){
		Particle::PID_Thrown(idx_,data_);
	}else{
		Particle::PID_Recon(idx_,data_,flags_);
	}
}

void Particle::PID_Thrown(int idx_, std::shared_ptr<Branches> data_){
	//std::cout<<"Identifying Thrown Particle\n";
	//std::cout<<"\tMomentum\n";
	_p = data_->Branches::mcp(idx_);
	//std::cout<<"\tTheta\n";
	_theta = data_->Branches::mctheta(idx_);
	//std::cout<<"\tPhi\n";
	_phi = data_->Branches::mcphi(idx_);
	//std::cout<<"\tID" <<data_->Branches::mcid(idx_) <<"\n";
	switch(data_->Branches::mcid(idx_)){
		case _ELECTRON_:
			_pid = {true,false,false,false};
		break;
		case _PROTON_:
			_pid = {false,true,false,false};
		break;
		case _PION_:
			_pid = {false,false,true,false};
		break;
		case -_PION_:
			_pid = {false,false,false,true};
		break;
		default:
			std::cout<<"Unrecognized Particle ID for Thrown Particle\n";
		break;
	}
	//std::cout<<"\tParticle Identified: " <<data_->Branches::mcid(idx_) <<"\n";
}

void Particle::PID_Recon(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
	//Kinematic Quantities
	if(idx_==0 && flags_->E_PCorr()){
		//if(flags_->E_Theta_Corr()){
			_p = corr::p_corr_e(data_->Branches::p(idx_),corr::theta_e_corr(physics::get_theta(data_->Branches::cz(idx_)),physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),flags_->Run()),physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),flags_->Run());
		//}else{
			//_p = corr::p_corr_e(data_->Branches::p(idx_),physics::get_theta(data_->Branches::cz(idx_)),physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),flags_->Run());
		//}
	}
	_p = data_->Branches::p(idx_);
	
	if(idx_==0 && (flags_->E_Theta_Corr() || flags_->E_PCorr())){
		_theta = corr::theta_e_corr(physics::get_theta(data_->Branches::cz(idx_)),physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),flags_->Run());
	}else{
		_theta = physics::get_theta(data_->Branches::cz(idx_));
	}
	if(flags_->Flags::Plot_SC_Eff() && data_->Branches::sc(idx_)>0){
		_sc_pd = data_->Branches::sc_pd(idx_);
	}
	_phi = physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_));
	_q = data_->Branches::q(idx_);
	for(int i=0; i<4; i++){
		_dt[i] = physics::delta_t(i,data_, idx_);
	}
	//Particle ID
	_pid = pid::pid(idx_,data_,flags_);
	_ided = true;
	if(idx_ == 0){//_Electron_
		//Detector Quantities
		_vz = data_->Branches::vz(idx_);
		_vx = data_->Branches::vx(idx_);
		_vy = data_->Branches::vy(idx_);
		_etot = data_->Branches::etot(idx_);
		_sf = (_etot/_p);
		_nphe = data_->Branches::nphe(idx_);
		_cc_seg = detect::cc_segment(data_->Branches::cc_segm(idx_));
		_cc_lrc = detect::cc_lrc(data_->Branches::cc_segm(idx_));
		if(_pid[0]){
			_sf_pass = true;
			_min_ec_pass = true;
			_cc_pass = true;
			_sanity_pass[0] = true;
			_dt_pass[0] = true;
			_fid_pass[0] = true;
			_vertex_pass = true;
			//_beta_pass[0] = true;
			_id_pass[0] = true;
		}else{
			_sanity_pass[0] = pid::sanity_ele(idx_,data_,flags_);
			if(_sanity_pass[0]){
				_sf_pass = pid::sf(idx_,data_,flags_);
				_min_ec_pass = pid::min_ec(idx_,data_,flags_);
				_cc_pass = pid::min_cc(idx_,data_,flags_);
				_vertex_pass = pid::vertex_e(idx_,data_,flags_);
				//_dt_pass[0] = pid::delta_t_ele(idx_,data_,flags_);//Maybe implement later
				_fid_pass[0] = pid::fid_ele(idx_,data_,flags_);
				//_beta_pass[0] = pid::
				_sc_eff_pass[0] = pid::sc_eff(0, idx_, data_, flags_);
				_id_pass[0] = pid::id_bank(idx_,data_,flags_,_ele_);
			}
		}
	}else{//Hadron
		if(_q > 0){
			if(_pid[1]){//Proton
				_sanity_pass[1] = true;
				_dt_pass[1] = true;
				_fid_pass[1] = true;
				_sc_eff_pass[1] = true;
				_id_pass[1] = true;
			}else{
				_sanity_pass[1] = pid::sanity_pro(idx_,data_,flags_);
				if(_sanity_pass[1]){
					_dt_pass[1] = pid::delta_t_pro(idx_,data_,flags_);
					_fid_pass[1] = pid::fid_pro(idx_,data_,flags_);
					_sc_eff_pass[1] = pid::sc_eff(1, idx_, data_, flags_);
				}
				_id_pass[1] = pid::id_bank(idx_,data_,flags_,_pro_);
			}
			if(_pid[2]){//Pip
				_sanity_pass[2] = true;
				_dt_pass[2] = true;
				_fid_pass[2] = true;
				_id_pass[2] = true;
				_sc_eff_pass[2] = true;
			}else{
				_sanity_pass[2] = pid::sanity_pip(idx_,data_,flags_);
				if(_sanity_pass[2]){
					_dt_pass[2] = pid::delta_t_pip(idx_,data_,flags_);
					_fid_pass[2] = pid::fid_pip(idx_,data_,flags_);
					_sc_eff_pass[2] = pid::sc_eff(2, idx_, data_, flags_);
				}
				_id_pass[2] = pid::id_bank(idx_,data_,flags_,_pip_);
			}
		}else if(_q < 0){
			if(_pid[3]){//Pim
				_sanity_pass[3] = true;
				_dt_pass[3] = true;
				_fid_pass[3] = true;
				_id_pass[3] = true;
				_sc_eff_pass[0] = true;
			}else{
				_sanity_pass[3] = pid::sanity_pim(idx_,data_,flags_);
				if(_sanity_pass[3]){
					_dt_pass[3] = pid::delta_t_pim(idx_,data_,flags_);
					_fid_pass[3] = pid::fid_pim(idx_,data_,flags_);
					_sc_eff_pass[3] = pid::sc_eff(3, idx_, data_, flags_);
				}
				_id_pass[3] = pid::id_bank(idx_,data_,flags_,_pim_);
			}
		}
	}
}

bool Particle::Pass_Sanity(int i){
	return _sanity_pass[i];
}
bool Particle::Pass_ec(){
	return _min_ec_pass;
}
bool Particle::Pass_fid(int i){
	return _fid_pass[i];
}
bool Particle::Pass_sf(){
	return _sf_pass;
}
bool Particle::Pass_cc(){
	return _cc_pass;
}
bool Particle::Pass_dt(int i){
	return _dt_pass[i];
}
bool Particle::Pass_id(int i){
	return _id_pass[i];
}
bool Particle::Pass_pid(int i){
	return _pid[i];
}
bool Particle::Pass_vertex(){
	return _vertex_pass;
}
bool Particle::Pass_SC_Eff(int i){
	return _sc_eff_pass[i];
}
bool Particle::Corr_p(){
	return _p_corr; 
}
int Particle::ID_crisis(){
	return _id_crisis; 
}
bool Particle::IDed(){
	return _ided;
}
bool Particle::Is_Sim(){
	return _sim;
}

bool Particle::Is_Thrown(){
	return _thrown;
}


bool Particle::Is_Elec(){
	return _pid[0];
}
bool Particle::Is_Pro(){
	return _pid[1];
}
bool Particle::Is_Pip(){
	return _pid[2];
}
bool Particle::Is_Pim(){
	return _pid[3];
}

float Particle::Get_p(){
	return _p; 
}
float Particle::Get_theta(){
	return _theta; 
}
float Particle::Get_phi(){
	return _phi; 
}

int Particle::Get_run(){
	return _run;
}

int Particle::Get_idx(){
	return _idx;
}

float Particle::Get_Weight(){
	return _weight;
}

int Particle::Get_q(){
	return _q;
}


int Particle::Get_sc_pd(){
	return _sc_pd;
}

float Particle::W(){//Will assume whatever particle is there is an electron
	TLorentzVector k_mu_prime = physics::Make_4Vector(true, _p, _theta, _phi, _me_);
	return physics::W(k_mu_prime,_run);
}

float Particle::Q2(){
	TLorentzVector k_mu_prime = physics::Make_4Vector(true, _p, _theta, _phi, _me_);
	return physics::Q2(k_mu_prime,_run);
}

int Particle::Sector(){
	return physics::get_sector(_phi);
}

TLorentzVector Particle::Get_4Vec(int i){
	TLorentzVector output; 
	if(_pid[i]){
		switch(i){
			case 0:
				output = physics::Make_4Vector(true, _p, _theta, _phi, _me_);
			break;
			case 1:
				output = physics::Make_4Vector(true, _p, _theta, _phi, _mp_);
			break;
			case 2:
				output = physics::Make_4Vector(true, _p, _theta, _phi, _mpi_);
			break;
			case 3:
				output = physics::Make_4Vector(true, _p, _theta, _phi, _mpi_);
			break;
			default:
				std::cout<<"Called for bad particle idx\n";
			break;
		}
	}else{
		std::cout<<"Called for a bad Particle 4Vec\tidx: "<<i <<"\n";
	}
	return output;
}

float Particle::Get_sf(){
	return _sf;
}
float Particle::Get_etot(){
	return _etot;
}
int Particle::Get_cc_seg(){
	return _cc_seg;
}
int Particle::Get_cc_lrc(){
	return _cc_lrc;
}
int Particle::Get_nphe(){
	return _nphe;
}
float Particle::Get_vz(){
	return _vz;
}
float Particle::Get_vx(){
	return _vx;
}
float Particle::Get_vy(){
	return _vy;
}
float Particle::Get_delta(int par_){
	return _dt[par_];
}
