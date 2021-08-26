#include "event.hpp"

Event::Event(int top_, Particle p0_, Particle p1_, Particle p2_, std::shared_ptr<Flags> flags_, float weight_, int hel_, bool thrown_){
	if(Event::Check_Particles(top_,p0_,p1_,p2_,flags_) && top_!= 3){
		_sim = flags_->Flags::Sim();
		_run = flags_->Flags::Run();
		if(_sim){
			_thrown = thrown_;
		}
		Event::Extract_Particles(top_,p0_,p1_,p2_,flags_);
		_weight = weight_;
		_top[top_] = true;
		if(cuts::MM_cut(top_,_MM,_W,flags_)){
			std::cout<<"\tPassed " <<_top_[top_] <<"\n";
			_pass[top_] = true;
			_pass_top = top_;
			_hel = hel_; 
			Event::Calc_Error();
		}
	}
}

Event::Event(int top_, Particle p0_, Particle p1_, Particle p2_, Particle p3_, std::shared_ptr<Flags> flags_, float weight_, int hel_, bool thrown_){
	if(Event::Check_Particles(top_,p0_,p1_,p2_,p3_,flags_) && top_== 3){
		_sim = flags_->Flags::Sim();
		_run = flags_->Flags::Run();
		if(_sim){
			_thrown = thrown_;
		}
		Event::Extract_Particles(top_,p0_,p1_,p2_,p3_,flags_);
		_weight = weight_;
		_top[top_] = true;
		if(cuts::MM_cut(top_,_MM,_W,flags_) && _MM != 0.0){
			std::cout<<"\tPassed " <<_top_[top_] <<"\n";
			_pass[top_] = true;
			_pass_top = top_;
			_hel = hel_; 
			Event::Calc_Error();
		}
	}
}

bool Event::Pass_Top(int i){
	return _pass[i];
}

bool Event::Pass(){
	if(_pass_top>-1){
		return _pass[_pass_top];
	}else{
		return false;
	}
}

int Event::Top(){
	int out = -1;
	for(int i=0; i<4; i++){
		if(_top[i]){
			out = i;
		}
	}
	return out;
}


bool Event::Check_Particles(int top_, Particle p0_, Particle p1_, Particle p2_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	if(p0_.Particle::Is_Elec()){
		//std::cout<<"\tChecking particles in Event ";
		_W = physics::W(p0_.Particle::Get_4Vec(0),flags_->Flags::Run());
		_Q2 = physics::Q2(p0_.Particle::Get_4Vec(0),flags_->Flags::Run());
		switch(top_){
			case 0:
				if(p1_.Particle::Is_Pip() && p2_.Particle::Is_Pim()){
					pass = true;
				}
			break;
			case 1:
				if(p1_.Particle::Is_Pro() && p2_.Particle::Is_Pim()){
					pass = true;
				}
			break;
			case 2:
				if(p1_.Particle::Is_Pro() && p2_.Particle::Is_Pip()){
					pass = true;
				}
			break;
			default:
				std::cout<<"You've entered the wrong topology\n";
			break;
		}
		//std::cout<<pass <<"\n";
	}else{
		//std::cout<<"\tBAD ELECTRON\n";
	}
	return pass;
}
bool Event::Check_Particles(int top_, Particle p0_, Particle p1_, Particle p2_, Particle p3_, std::shared_ptr<Flags> flags_){
	bool pass = false;
	if(top_ == 3){
		if(p0_.Particle::Is_Elec() && p1_.Particle::Is_Pro() && p2_.Particle::Is_Pip() && p3_.Particle::Is_Pim()){
			_W = physics::W(p0_.Particle::Get_4Vec(0),flags_->Flags::Run());
			_Q2 = physics::Q2(p0_.Particle::Get_4Vec(0),flags_->Flags::Run());
			pass = true;
		}
	}
	return pass;
}

void Event::Extract_Particles(int top_, Particle p0_, Particle p1_, Particle p2_, std::shared_ptr<Flags> flags_){
	_k1_lab = physics::Make_4Vector(_beam_energy_[_run],0.0,0.0,1.0,_me_);
	switch(top_){
			case 0:
				_vec_lab[0] = p0_.Particle::Get_4Vec(0);
				_vec_lab[2] = p1_.Particle::Get_4Vec(2);
				_vec_lab[3] = p2_.Particle::Get_4Vec(3);
				_MM = physics::MM_event(flags_->Flags::Run(),0,_vec_lab[0],_vec_lab[2],_vec_lab[3]);
				_MM2 = _MM*_MM;
			break;
			case 1:
				_vec_lab[0] = p0_.Particle::Get_4Vec(0);
				_vec_lab[1] = p1_.Particle::Get_4Vec(1);
				_vec_lab[3] = p2_.Particle::Get_4Vec(3);
				_MM = physics::MM_event(flags_->Flags::Run(),0,_vec_lab[0],_vec_lab[1],_vec_lab[3]);
				_MM2 = _MM*_MM;
			break;
			case 2:
				_vec_lab[0] = p0_.Particle::Get_4Vec(0);
				_vec_lab[1] = p1_.Particle::Get_4Vec(1);
				_vec_lab[2] = p2_.Particle::Get_4Vec(2);
				_MM = physics::MM_event(flags_->Flags::Run(),0,_vec_lab[0],_vec_lab[1],_vec_lab[2]);
				_MM2 = _MM*_MM;
			break;
			default:
				std::cout<<"You've entered the wrong topology\n";
			break;
	}
	Event::Get_Angles();
}
void Event::Extract_Particles(int top_, Particle p0_, Particle p1_, Particle p2_, Particle p3_, std::shared_ptr<Flags> flags_){
	_k1_lab = physics::Make_4Vector(_beam_energy_[_run],0.0,0.0,1.0,_me_);
	if(top_==3){
		_vec_lab[0] = p0_.Particle::Get_4Vec(0);
		_vec_lab[1] = p1_.Particle::Get_4Vec(1);
		_vec_lab[2] = p2_.Particle::Get_4Vec(2);
		_vec_lab[3] = p3_.Particle::Get_4Vec(3);
		_MM = physics::MM_event(flags_->Flags::Run(),0,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
	}
	Event::Get_Angles();
}

void Event::Calc_Error(){
	float m = -99.9; 
	switch(_pass_top){
		case 0:
			m = _mp_;
		break;
		case 1:
			m = _mpi_;
		break;
		case 2:
			m = _mpi_;
		break;
		case 3:
			m = 0.0;
		break;
		default:
		break;
	}
	_error =  TMath::Power((m-_MM)/m,2.0);
}

float Event::Error(){
	return _error;
}

float Event::Weight(){
	return _weight;
}
float Event::MM(){
	return _MM;
}

float Event::MM2(){
	return _MM2;
}

float Event::Theta(int particle_){
	return _theta_lab[particle_];
}
float Event::Phi(int particle_){
	return _phi_lab[particle_];
}

float Event::W(){
	return _W;
}

float Event::Q2(){
	return _Q2;
}

float Event::MMb(int i){
	return _MMb[i];
}

float Event::MM2b(int i){
	return _MM2b[i];
}

float Event::Thetab(int i){
	return _thetab[i];
}

float Event::Phib(int i){
	return _phib[i];
}

float Event::Alphab(int i){
	return _alphab[i];
}

void Event::Missing_Hadron(){
	switch(_pass_top){
		case 0:
			_vec_lab[1] = _k1_lab + _p_mu_ - _vec_lab[0] - _vec_lab[2] - _vec_lab[3];
			_full_event = true;
		break;
		case 1:
			_vec_lab[2] = _k1_lab + _p_mu_ - _vec_lab[0] - _vec_lab[1] - _vec_lab[3];
			_full_event = true;
		break;
		case 2:
			_vec_lab[3] = _k1_lab + _p_mu_ - _vec_lab[0] - _vec_lab[1] - _vec_lab[2];
			_full_event = true;
		break;
		case 3:
			_full_event = true;
		break;
		default:
		break;
	}
}


void Event::COM(){
	if(_full_event){
			  //physics::COM_gp(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
		_vec[0] = physics::COM_gp(0,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
		_vec[1] = physics::COM_gp(1,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
		_vec[2] = physics::COM_gp(2,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
		_vec[3] = physics::COM_gp(3,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
		_k1 = physics::COM_gp(4,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
		_p1 = physics::COM_gp(5,_k1,_vec_lab[0],_vec_lab[1],_vec_lab[2],_vec_lab[3]);
		_COM = true;
	}
	Event::Get_Angles(_COM);
}

void Event::Vars(){
	for(int i = 0; i<3; i++){
		if(_COM){
			_alphab[i] = physics::alpha(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
			_thetab[i] = physics::Ev_Theta(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
			_MMb[i] = physics::Ev_MM(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
			_MM2b[i] = physics::Ev_MM2(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
			_phib[i] = physics::Ev_Phi(i, _k1, _vec[0], _vec[1], _vec[2], _vec[3], true);
		}else{
			std::cout<<"You need to get your 4-Momenta into center of mass" <<std::endl;
		}
	}
}

int Event::Helicity(){
	return _hel;
}

float Event::Get_Px(int particle_, bool COM_){
	if(COM_){
		return _vec[particle_][0];
	}else{
		return _vec_lab[particle_][0];
	}
}
float Event::Get_Py(int particle_, bool COM_){
	if(COM_){
		return _vec[particle_][1];
	}else{
		return _vec_lab[particle_][1];
	}
}
float Event::Get_Pz(int particle_, bool COM_){
	if(COM_){
		return _vec[particle_][2];
	}else{
		return _vec_lab[particle_][2];
	}
}
float Event::Get_P0(int particle_, bool COM_){
	if(COM_){
		return _vec[particle_][3];
	}else{
		return _vec_lab[particle_][3];
	}
}

float Event::Get_Beam_Comp(int component_, bool COM_){
	if(COM_){
		return _k1[component_];
	}else{
		return _k1_lab[component_];
	}
}

float Event::Get_Target_Comp(int component_, bool COM_){
	if(COM_){
		return _p1[component_];
	}else{
		return _p_mu_[component_];
	}
}

int Event::Get_PID(int particle_){
	switch(particle_){
		case 0:
			return _ELECTRON_;
		break;
		case 1:
			return _PROTON_;
		break;
		case 2:
			return _PION_;
		break;
		case 3:
			return -_PION_;
		break;
		default:
		break;
	}
}

void Event::Get_Angles(bool COM_){
	for(int i=0; i<4; i++){
		if(_part[i]){
			if(COM_){
				_phi[i] = physics::get_phi(_vec[i][0],_vec[i][1]);
			}else{
				_phi_lab[i] = physics::get_phi(_vec_lab[i][0],_vec_lab[i][1]);
			}
		}
	}
}

bool Event::Was_COM(){
	return _COM;
}

int Event::Run(){
	return _run;
}

int Event::Sector(int particle_){
	return physics::get_sector(_phi_lab[particle_]);
}