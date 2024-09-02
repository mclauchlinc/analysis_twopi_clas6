#include "pid.hpp"

	std::vector<bool> pid::pid(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		//std::cout<<"Begin Particle ID\n";
		std::vector<bool> pass;//[4]={false,false,false,false};
		if(idx_==0){
			//std::cout<<"Particle ID for Electron\n";
			pass.push_back(pid::pid_ele(idx_,data_,flags_));
			pass.push_back(false);//Pro
			pass.push_back(false);//Pip
			pass.push_back(false);//Pim
		}else{
			pass.push_back(false);//Ele
			//std::cout<<"Particle ID for Hadrons idx:" <<idx_ <<"\n";
			//std::cout<<"Particle Charge: " <<data_->Branches::q(idx_) <<std::endl;
			if(data_->Branches::q(idx_) > 0){
				//std::cout<<"\tParticle ID for Proton\n";
				pass.push_back(pid::pid_pro(idx_,data_,flags_));
				//std::cout<<"\tParticle ID for Pi+\n";
				pass.push_back(pid::pid_pip(idx_,data_,flags_));
			}else{
				pass.push_back(false);//Pro
				pass.push_back(false);//Pip
			}
			if(data_->Branches::q(idx_) < 0){
				//std::cout<<"\tParticle ID for Pi-\n";
				pass.push_back(pid::pid_pim(idx_,data_,flags_));
			}else{
				pass.push_back(false);//Pim
			}
		}
		return pass;
	}
	bool pid::pid_ele(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		//std::cout<<"\t\tIn Particle ID for Electron\n";
		bool pass = false;
		if(idx_==0){
			pass = pid::sanity_ele(idx_,data_,flags_);
			//pass &= pid::fid_ele(idx_,data_,flags_);
			pass &= pid::sf(idx_,data_,flags_);
			pass &= pid::min_ec(idx_,data_,flags_);
			//pass &= pid::sc_eff(0, idx_, data_, flags_);
			if(!flags_->Flags::Sim()){
				pass &= pid::min_cc(idx_,data_,flags_);
				pass &= pid::vertex_e(idx_,data_,flags_);
			}
			
			//pass &= pid::id_bank(idx_,data_,flags_,_ele_);
			pass &= pid::geo_cc_cut(0,idx_,data_,flags_);
			pass &= pid::geo_sc_cut(0,idx_,data_,flags_);
			pass &= pid::geo_ec_cut(0,idx_,data_,flags_);
			pass &= pid::good_sc_paddle(0,idx_,data_,flags_);
			pass &= pid::good_cc_segment(0,idx_,data_,flags_);
			pass &= pid::kin_eff_cut(0,idx_,data_,flags_);
			/*
			std::cout<<"ele sanity:" <<pid::sanity_ele(idx_,data_,flags_) <<"\n";
			std::cout<<"ele fid_ele:" <<pid::fid_ele(idx_,data_,flags_) <<"\n";
			std::cout<<"ele sf:" <<pid::sf(idx_,data_,flags_) <<"\n";
			std::cout<<"ele min ec:" <<pid::min_ec(idx_,data_,flags_) <<"\n";
			//std::cout<<"ele sanity:" <<pid::sc_eff(0, idx_, data_, flags_);
			std::cout<<"ele min cc:" <<pid::min_cc(idx_,data_,flags_) <<"\n";
			std::cout<<"ele vertex:" <<pid::vertex_e(idx_,data_,flags_) <<"\n";
			std::cout<<"ele id:" <<pid::id_bank(idx_,data_,flags_,_ele_) <<"\n";
			std::cout<<"ele fid geo cc:" <<pid::geo_cc_cut(0,idx_,data_,flags_) <<"\n";
			std::cout<<"ele fid geo sc:" <<pid::geo_sc_cut(0,idx_,data_,flags_) <<"\n";
			std::cout<<"ele fid geo ec:" <<pid::geo_ec_cut(0,idx_,data_,flags_) <<"\n";
			std::cout<<"ele kin_eff_cut:" <<pid::kin_eff_cut(0,idx_,data_,flags_) <<"\n";
			*/
		}
		//std::cout<<"\t\tIn Particle ID for ELE| pass:" <<pass <<"\n";
		return pass; 
	}
	bool pid::pid_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		//std::cout<<"\t\tIn Particle ID for Proton\n";
		bool pass = false;
		if(idx_ > 0){
			pass = pid::sanity_pro(idx_,data_,flags_);
			pass &= pid::delta_t_pro(idx_,data_,flags_);
			//pass &= pid::sc_eff(1, idx_, data_, flags_);
			//pass &= pid::fid_pro(idx_,data_,flags_);
			//pass &= pid::id_bank(idx_,data_,flags_,_pro_);
			//pass &= pid::geo_cc_cut(1,idx_,data_,flags_);
			pass &= pid::geo_sc_cut(1,idx_,data_,flags_);
			//pass &= pid::geo_ec_cut(1,idx_,data_,flags_);
			pass &= pid::kin_eff_cut(1,idx_,data_,flags_);
			pass &= pid::good_sc_paddle(1,idx_,data_,flags_);
			/*
			std::cout<<"pro sanity:" <<pid::sanity_pro(idx_,data_,flags_) <<"\n";
			std::cout<<"pro delta t:" <<pid::delta_t_pro(idx_,data_,flags_) <<"\n";
			//std::cout<<"pro sanity:" <<pid::sc_eff(3, idx_, data_, flags_);
			std::cout<<"pro fid:" <<pid::fid_pro(idx_,data_,flags_) <<"\n";
			std::cout<<"pro id:" <<pid::id_bank(idx_,data_,flags_,_pro_) <<"\n";
			//std::cout<<"pro sanity:" <<pid::geo_cc_cut(3,idx_,data_,flags_);
			std::cout<<"pro fid_geo_sc:" <<pid::geo_sc_cut(1,idx_,data_,flags_) <<"\n";
			//std::cout<<"pro sanity:" <<pid::geo_ec_cut(3,idx_,data_,flags_);
			std::cout<<"pro kin_eff:" <<pid::kin_eff_cut(1,idx_,data_,flags_) <<"\n";
			*/
		}
		//std::cout<<"\t\tIn Particle ID for PRO| pass:" <<pass <<"\n";
		return pass;
	}
	bool pid::pid_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		
		bool pass = false;
		if(idx_ > 0){
			pass = pid::sanity_pip(idx_,data_,flags_);
			pass &= pid::delta_t_pip(idx_,data_,flags_);
			//pass &= pid::sc_eff(2, idx_, data_, flags_);
			//pass &= pid::fid_pip(idx_,data_,flags_);
			//pass &= pid::id_bank(idx_,data_,flags_,_pip_);
			//pass &= pid::geo_cc_cut(2,idx_,data_,flags_);
			pass &= pid::geo_sc_cut(2,idx_,data_,flags_);
			//pass &= pid::geo_ec_cut(2,idx_,data_,flags_);
			pass &= pid::kin_eff_cut(2,idx_,data_,flags_);
			pass &= pid::good_sc_paddle(2,idx_,data_,flags_);
			/*
			std::cout<<"pip sanity:" <<pid::sanity_pip(idx_,data_,flags_) <<"\n";
			std::cout<<"pip delta t:" <<pid::delta_t_pip(idx_,data_,flags_) <<"\n";
			//std::cout<<"pip sanity:" <<pid::sc_eff(3, idx_, data_, flags_);
			std::cout<<"pip fid:" <<pid::fid_pip(idx_,data_,flags_) <<"\n";
			std::cout<<"pip id:" <<pid::id_bank(idx_,data_,flags_,_pip_) <<"\n";
			//std::cout<<"pip sanity:" <<pid::geo_cc_cut(3,idx_,data_,flags_);
			std::cout<<"pip fid_geo_sc:" <<pid::geo_sc_cut(2,idx_,data_,flags_) <<"\n";
			//std::cout<<"pip sanity:" <<pid::geo_ec_cut(3,idx_,data_,flags_);
			std::cout<<"pip kin_eff:" <<pid::kin_eff_cut(2,idx_,data_,flags_) <<"\n";
			*/
		}
		//std::cout<<"\t\tIn Particle ID for PIP| pass:" <<pass <<"\n";
		return pass;
	}
	bool pid::pid_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		//std::cout<<"\t\tIn Particle ID for PIM\n";
		bool pass = false;
		if(idx_ > 0){
			pass = pid::sanity_pim(idx_,data_,flags_);
			pass &= pid::delta_t_pim(idx_,data_,flags_);
			//pass &= pid::sc_eff(3, idx_, data_, flags_);
			//pass &= pid::fid_pim(idx_,data_,flags_);
			pass &= pid::id_bank(idx_,data_,flags_,_pim_);
			//pass &= pid::geo_cc_cut(3,idx_,data_,flags_);
			pass &= pid::geo_sc_cut(3,idx_,data_,flags_);
			//pass &= pid::geo_ec_cut(3,idx_,data_,flags_);
			pass &= pid::kin_eff_cut(3,idx_,data_,flags_);
			pass &= pid::good_sc_paddle(3,idx_,data_,flags_);
			/*std::cout<<"pim sanity:" <<pid::sanity_pim(idx_,data_,flags_) <<"\n";
			std::cout<<"pim delta t:" <<pid::delta_t_pim(idx_,data_,flags_) <<"\n";
			//std::cout<<"pim sanity:" <<pid::sc_eff(3, idx_, data_, flags_);
			std::cout<<"pim fid:" <<pid::fid_pim(idx_,data_,flags_) <<"\n";
			std::cout<<"pim id:" <<pid::id_bank(idx_,data_,flags_,_pim_) <<"\n";
			//std::cout<<"pim sanity:" <<pid::geo_cc_cut(3,idx_,data_,flags_);
			std::cout<<"pim fid_geo_sc:" <<pid::geo_sc_cut(3,idx_,data_,flags_) <<"\n";
			//std::cout<<"pim sanity:" <<pid::geo_ec_cut(3,idx_,data_,flags_);
			std::cout<<"pim kin_eff:" <<pid::kin_eff_cut(3,idx_,data_,flags_) <<"\n";
			*/
		}
		//std::cout<<"\t\tIn Particle ID for PIM| pass:" <<pass <<"\n";
		return pass;
	}
	bool pid::id_bank(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_,const char* _species_){
		if(!flags_->Flags::ID_Cut()){
			return true;
		}
		if(_species_ == _ele_ && data_->Branches::id(idx_) == _ELECTRON_){
			return true;
		}else if(_species_ == _pro_ && data_->Branches::id(idx_) == _PROTON_){
			return true;
		}else if(_species_ == _pip_ && data_->Branches::id(idx_) == _PION_){
			return true;
		}else if(_species_ == _pim_ && data_->Branches::id(idx_) == -_PION_){
			return true;
		}else{
			return false;
		}
	}
	std::vector<bool> pid::fid(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		std::vector<bool> pass;//[4] = {false,false,false,false};
		if(idx_ == 0){
			pass.push_back(pid::fid_ele(idx_,data_,flags_));
			pass.push_back(false);
			pass.push_back(false);
			pass.push_back(false);
		}else{
			pass.push_back(false);
			if(data_->Branches::q(idx_) > 0){
				pass.push_back(pid::fid_pro(idx_,data_,flags_));
				pass.push_back(pid::fid_pip(idx_,data_,flags_));
			}else{
				pass.push_back(false);
				pass.push_back(false);
			}
			if(data_->Branches::q(idx_) < 0){
				pass.push_back(pid::fid_pim(idx_,data_,flags_));
			}else{
				pass.push_back(false);
			}
		}
		return pass;
	}
	bool pid::fid_ele(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(!flags_->Flags::Fid_Cut(0)){return true;}
		if(idx_ == 0){
			pass = cuts::fid_e(data_->Branches::p(idx_), physics::get_theta(data_->Branches::cz(idx_)), physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)), flags_);
		}
		return pass;
	}
	bool pid::fid_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(!flags_->Flags::Fid_Cut(1)){return true;}
		if(idx_ > 0){
			pass = cuts::fid_h(0, physics::get_theta(data_->Branches::cz(idx_)), physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),  flags_);
		}
		return pass;
	}
	bool pid::fid_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(!flags_->Flags::Fid_Cut(2)){return true;}
		if(idx_ > 0){
			pass = cuts::fid_h(1, physics::get_theta(data_->Branches::cz(idx_)), physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),  flags_);
		}
		return pass;
	}
	bool pid::fid_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(!flags_->Flags::Fid_Cut(3)){return true;}
		if(idx_ > 0){
			pass = cuts::fid_pim(data_->Branches::p(idx_), physics::get_theta(data_->Branches::cz(idx_)), physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)), flags_);
		}
		return pass;
	}
	std::vector<bool> pid::delta_t(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		std::cout<<"\tPerforming Delta t cut on: " <<idx_ <<"\n";
		std::vector<bool> pass;// = {false,false,false,false};
		if(idx_ == 0){
			pass.push_back(pid::delta_t_ele(idx_,data_,flags_));
			pass.push_back(false);
			pass.push_back(false);
			pass.push_back(false);
		}else{
			pass.push_back(false);
			if(data_->Branches::q(idx_) > 0){
				pass.push_back(pid::delta_t_pro(idx_,data_,flags_));
				pass.push_back(pid::delta_t_pip(idx_,data_,flags_));
			}else{
				pass.push_back(false);
				pass.push_back(false);
			}
			if(data_->Branches::q(idx_) < 0){
				pass.push_back(pid::delta_t_pim(idx_,data_,flags_));
			}else{
				pass.push_back(false);
			}
		}
		return pass;
	}
	bool pid::delta_t_ele(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		std::cout<<"\t\tDelta t for electron\n";
		bool pass = false;
		if(!flags_->Flags::Delta_Cut(0)){return true;}
		if(idx_ == 0){
			int run = flags_->Flags::Run();
			int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
			bool sim = flags_->Flags::Sim();
			float p = data_->Branches::p(idx_);
			float dt = physics::delta_t(0,data_,idx_);
			pass = cuts::delta_t_cut(0,p, dt,sector, flags_);
		}
		return pass;
	}
	bool pid::delta_t_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(!flags_->Flags::Delta_Cut(1)){return true;}
		if(idx_ > 0){
			int run = flags_->Flags::Run();
			int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
			bool sim = flags_->Flags::Sim();
			float p = data_->Branches::p(idx_);
			float dt = physics::delta_t(1,data_,idx_);
			pass = cuts::delta_t_cut(1,p, dt,sector, flags_);
		}
		return pass;
	}
	bool pid::delta_t_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(!flags_->Flags::Delta_Cut(2)){return true;}
		if(idx_ > 0){
			int run = flags_->Flags::Run();
			int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
			bool sim = flags_->Flags::Sim();
			float p = data_->Branches::p(idx_);
			float dt = physics::delta_t(2,data_,idx_);
			pass = cuts::delta_t_cut(2,p, dt,sector, flags_);
		}
		return pass;
	}
	bool pid::delta_t_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(!flags_->Flags::Delta_Cut(3)){return true;}
		if(idx_ > 0){
			int run = flags_->Flags::Run();
			int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
			bool sim = flags_->Flags::Sim();
			float p = data_->Branches::p(idx_);
			float dt = physics::delta_t(3,data_,idx_);
			pass = cuts::delta_t_cut(3,p, dt,sector, flags_);
		}
		return pass;
	}
	bool pid::sf(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(!flags_->Flags::SF_Cut()){return true;}
		if(idx_ == 0){
			int run = flags_->Flags::Run();
			int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
			bool sim = flags_->Flags::Sim();
			float p = data_->Branches::p(idx_);
			float sf = physics::sf(data_->Branches::p(idx_),data_->Branches::etot(idx_));
			pass = cuts::sf_cut(p,sf,physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),flags_);
		}
		return pass;
	}
	bool pid::min_cc(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(!flags_->Flags::CC_Cut()){return true;}
		if(idx_ == 0){
			bool sim = flags_->Flags::Sim();
			if(!sim){
				int run = flags_->Flags::Run();
				int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
				int cc_segm = data_->Branches::cc_segm(idx_);
				int cc_sect = data_->Branches::cc_sect(idx_);
				int nphe = data_->Branches::nphe(idx_);
				pass = cuts::min_cc(cc_segm, cc_sect, nphe, flags_);
			}
		}
		return pass;
	}
	bool pid::min_ec(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::EC_Cut()){return true;}
		if(idx_ == 0){
			int run = flags_->Flags::Run();
			bool sim = flags_->Flags::Sim(); 
			int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
			float etot = data_->Branches::etot(idx_);
			pass = cuts::min_ec(etot,sector,flags_);
		}
		
		return pass;
	}
	bool pid::sanity(int particle_, int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass[4] = {false,false,false,false};
		if(particle_ == 0 && idx_ == 0){
			pass[0] = pid::sanity_ele(idx_,data_,flags_);
		}
		if(particle_ > 0 && idx_ > 0){
			if(data_->Branches::q(idx_) > 0){
				pass[1] = pid::sanity_pro(idx_,data_,flags_);
				pass[2] = pid::sanity_pip(idx_,data_,flags_);
			}else if(data_->Branches::q(idx_) < 0){
				pass[3] = pid::sanity_pim(idx_,data_,flags_);
			}
		}
		return pass;
	}
	bool pid::sanity_ele(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = (idx_==0);
		int dc = data_->Branches::dc(idx_);
		int sc = data_->Branches::sc(idx_);
		int ec = data_->Branches::ec(idx_);
		int cc = data_->Branches::cc(idx_);
		int stat = data_->Branches::stat(idx_);
		pass &= cuts::e_sanity(dc,sc,ec,cc,stat);
		return pass;
	}
	bool pid::sanity_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = (idx_ >0);
		pass &= cuts::pro_sanity(data_->Branches::dc(idx_),data_->Branches::sc(idx_),data_->Branches::stat(idx_));
		return pass;
	}
	bool pid::sanity_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = (idx_ >0);
		pass &= cuts::pip_sanity(data_->Branches::dc(idx_),data_->Branches::sc(idx_),data_->Branches::stat(idx_));
		return pass;
	}
	bool pid::sanity_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = (idx_ >0);
		pass &= cuts::pim_sanity(data_->Branches::dc(idx_),data_->Branches::sc(idx_),data_->Branches::stat(idx_));
		return pass;
	}

	bool pid::vertex_e(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		if(!flags_->Flags::Vertex_Cut()){return true;}
		int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
		return cuts::vertex_cut(data_->Branches::vz(idx_),flags_->Flags::Run(),sector,flags_);
	}

	//Efficiency Cuts
	bool pid::sc_eff(int par_, int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		if(!flags_->Flags::SC_Eff()){return true;}
		return cuts::sc_eff_cut(data_->Branches::p(idx_), physics::get_theta(data_->Branches::cz(idx_)), flags_->Flags::Run(), par_, flags_);
	}

	bool pid::geo_cc_cut(int par_, int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		if(!flags_->Flags::Geo_Cut(par_,0)){return true;}
		bool pass = true;
		int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
		pass&=cuts::cc_geo_cut(par_, detect::cc_x(data_,idx_),  detect::cc_y(data_,idx_),  sector-1,  detect::cc_lrc(data_->Branches::cc_segm(idx_)),  detect::cc_segment(data_->Branches::cc_segm(idx_)), flags_);
		
		return pass;
	}
	bool pid::geo_sc_cut(int par_, int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		if(!flags_->Flags::Geo_Cut(par_,1)){return true;}
		bool pass = true;
		int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
		
		pass &= cuts::sc_geo_cut(par_, data_->Branches::dc_xsc(idx_), data_->Branches::dc_ysc(idx_), sector-1,data_->Branches::sc_pd(idx_),flags_);
		//if(!pass && !isnan(data_->Branches::dc_xsc(idx_)) && !isnan(data_->Branches::dc_ysc(idx_))){
		//if(_species_[par_]==_pim_){
		//	std::cout<<"PID sc_geo_cut: species:" <<_species_[par_] <<" x:" <<data_->Branches::dc_xsc(idx_) <<" y:" <<data_->Branches::dc_ysc(idx_) <<" sec:" <<sector-1 <<" pad:" <<data_->Branches::sc_pd(idx_) <<" pass:" <<pass <<"\n";
		//}
		return pass;
	}
	bool pid::geo_ec_cut(int par_, int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		if(!flags_->Flags::Geo_Cut(par_,2)){return true;}
		bool pass = true;
		int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
		pass &= cuts::ec_geo_cut(par_,data_->Branches::ech_x(idx_),data_->Branches::ech_y(idx_),sector-1,flags_);
		return pass;
	}
	bool pid::kin_eff_cut(int par_, int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		if(!flags_->Flags::Kin_Eff_Cut(par_)){return true;}
		bool pass = true;
		float p = data_->Branches::p(idx_);
		float theta = physics::get_theta(data_->Branches::cz(idx_));
		if(par_ == 0){
			if(flags_->E_Theta_Corr()){
				theta = corr::theta_e_corr(physics::get_theta(data_->Branches::cz(idx_)),physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),flags_->Run());
				if(flags_->E_PCorr()){
					p = corr::p_corr_e(data_->Branches::p(idx_),corr::theta_e_corr(physics::get_theta(data_->Branches::cz(idx_)),physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),flags_->Run()),physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),flags_->Run());
				}
			}
		}
		int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
		pass &= cuts::kin_eff_cut(par_,sector,p,theta,flags_);
		return pass;
	}

	bool pid::good_sc_paddle(int par_, int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
		return cuts::sc_good_paddle(par_,sector-1,data_->Branches::sc_pd(idx_),flags_);
		

	}
	bool pid::good_cc_segment(int par_, int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
		return cuts::cc_good_segment(par_,sector-1,detect::cc_segment(data_->Branches::cc_segm(idx_)),flags_);
	}