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
			//std::cout<<"\t\tSanity Electron: ";
			pass = pid::sanity_ele(idx_,data_,flags_);
			//std::cout<<pass <<"\n\t\tFid Electron: ";
			pass &= pid::fid_ele(idx_,data_,flags_);
			//std::cout<<pass <<"\n\t\tSF Electron: ";
			pass &= pid::sf(idx_,data_,flags_);
			//std::cout<<pass <<"\n\t\tEC Electron: ";
			pass &= pid::min_ec(idx_,data_,flags_);
			//pass &= pid::min_p(idx_,data_,flags_);
			//std::cout<<pass <<"\n\t\tCC Electron: ";
			pass &= pid::min_cc(idx_,data_,flags_);
			//std::cout<<pass <<"\n\t\tVertex Electron: ";
			pass &= pid::vertex_e(idx_,data_,flags_);
			//std::cout<<pass <<"\n";
		}
		return pass; 
	}
	bool pid::pid_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		//std::cout<<"\t\tIn Particle ID for Proton\n";
		bool pass = false;
		if(idx_ > 0){
			//std::cout<<"\t\tSanity for Proton\n";
			pass = pid::sanity_pro(idx_,data_,flags_);
			//std::cout<<"\t\tDelta for Proton\n";
			pass &= pid::delta_t_pro(idx_,data_,flags_);
			//std::cout<<"\t\tFid for Proton\n";
			pass &= pid::fid_pro(idx_,data_,flags_);
		}
		return pass;
	}
	bool pid::pid_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		//std::cout<<"\t\tIn Particle ID for Proton\n";
		bool pass = false;
		if(idx_ > 0){
			//std::cout<<"\t\tSanity for Pi+\n";
			pass = pid::sanity_pip(idx_,data_,flags_);
			//std::cout<<"\t\tDelta for Pi+\n";
			pass &= pid::delta_t_pip(idx_,data_,flags_);
			//std::cout<<"\t\tFid for Pi+\n";
			pass &= pid::fid_pip(idx_,data_,flags_);
		}
		return pass;
	}
	bool pid::pid_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(idx_ > 0){
			pass = pid::sanity_pim(idx_,data_,flags_);
			pass &= pid::delta_t_pim(idx_,data_,flags_);
			pass &= pid::fid_pim(idx_,data_,flags_);
		}
		return pass;
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
		if(flags_->Flags::Fid_Cut(0)){
			if(idx_ == 0){
				pass = cuts::fid_e(data_->Branches::p(idx_), physics::get_theta(data_->Branches::cz(idx_)), physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)), flags_);
			}
		}else{
			pass = true;
		}
		return pass;
	}
	bool pid::fid_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::Fid_Cut(1)){
			if(idx_ > 0){
				pass = cuts::fid_h(0, physics::get_theta(data_->Branches::cz(idx_)), physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),  flags_);
			}
		}else{
			pass = true;
		}
		return pass;
	}
	bool pid::fid_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::Fid_Cut(2)){
			if(idx_ > 0){
				pass = cuts::fid_h(1, physics::get_theta(data_->Branches::cz(idx_)), physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),  flags_);
			}
		}else{
			pass = true;
		}
		return pass;
	}
	bool pid::fid_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::Fid_Cut(3)){
			if(idx_ > 0){
				pass = cuts::fid_pim(data_->Branches::p(idx_), physics::get_theta(data_->Branches::cz(idx_)), physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)), flags_);
			}
		}else{
			pass = true;
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
		if(flags_->Flags::Delta_Cut(0)){
			if(idx_ == 0){
				int run = flags_->Flags::Run();
				int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
				bool sim = flags_->Flags::Sim();
				float p = data_->Branches::p(idx_);
				float dt = physics::delta_t(0,data_,idx_);
				pass = cuts::delta_t_cut(0,p, dt, flags_);
			}
		}else{
			pass = true;
		}
		return pass;
	}
	bool pid::delta_t_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::Delta_Cut(1)){
			if(idx_ > 0){
				int run = flags_->Flags::Run();
				int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
				bool sim = flags_->Flags::Sim();
				float p = data_->Branches::p(idx_);
				float dt = physics::delta_t(1,data_,idx_);
				pass = cuts::delta_t_cut(1,p, dt, flags_);
			}
		}else{
			pass = true;
		}
		return pass;
	}
	bool pid::delta_t_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::Delta_Cut(2)){
			if(idx_ > 0){
				int run = flags_->Flags::Run();
				int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
				bool sim = flags_->Flags::Sim();
				float p = data_->Branches::p(idx_);
				float dt = physics::delta_t(2,data_,idx_);
				pass = cuts::delta_t_cut(2,p, dt, flags_);
			}
		}else{
			pass = true;
		}
		return pass;
	}
	bool pid::delta_t_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::Delta_Cut(3)){
			if(idx_ > 0){
				int run = flags_->Flags::Run();
				int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
				bool sim = flags_->Flags::Sim();
				float p = data_->Branches::p(idx_);
				float dt = physics::delta_t(3,data_,idx_);
				pass = cuts::delta_t_cut(3,p, dt, flags_);
			}
		}else{
			pass = true;
		}
		return pass;
	}
	bool pid::sf(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::SF_Cut()){
			if(idx_ == 0){
				int run = flags_->Flags::Run();
				int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
				bool sim = flags_->Flags::Sim();
				float p = data_->Branches::p(idx_);
				float sf = physics::sf(data_->Branches::p(idx_),data_->Branches::etot(idx_));
				pass = cuts::sf_cut(p,sf,physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)),flags_);
			}
		}else{
			pass = true;
		}
		return pass;
	}
	bool pid::min_cc(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::CC_Cut()){
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
		}else{
			pass = true;
		}
		return pass;
	}
	bool pid::min_ec(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::EC_Cut()){
			if(idx_ == 0){
				int run = flags_->Flags::Run();
				bool sim = flags_->Flags::Sim(); 
				int sector = physics::get_sector(physics::get_phi(data_->Branches::cx(idx_),data_->Branches::cy(idx_)));
				float etot = data_->Branches::etot(idx_);
				pass = cuts::min_ec(etot,sector,flags_);
			}
		}else{
			pass = true;
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
		bool pass = false;
		if(idx_ == 0){
			int dc = data_->Branches::dc(idx_);
			int sc = data_->Branches::sc(idx_);
			int ec = data_->Branches::ec(idx_);
			int cc = data_->Branches::cc(idx_);
			int stat = data_->Branches::stat(idx_);
			pass = cuts::e_sanity(dc,sc,ec,cc,stat);
		}
		return pass;
	}
	bool pid::sanity_pro(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(idx_ > 0){
			int dc = data_->Branches::dc(idx_);
			int sc = data_->Branches::sc(idx_);
			int stat = data_->Branches::stat(idx_);
			pass = cuts::pro_sanity(dc,sc,stat);
		}
		return pass;
	}
	bool pid::sanity_pip(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(idx_ > 0){
			int dc = data_->Branches::dc(idx_);
			int sc = data_->Branches::sc(idx_);
			int stat = data_->Branches::stat(idx_);
			pass = cuts::pip_sanity(dc,sc,stat);
		}
		return pass;
	}
	bool pid::sanity_pim(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(idx_ > 0){
			int dc = data_->Branches::dc(idx_);
			int sc = data_->Branches::sc(idx_);
			int stat = data_->Branches::stat(idx_);
			pass = cuts::pim_sanity(dc,sc,stat);
		}
		return pass;
	}

	bool pid::vertex_e(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_){
		bool pass = false;
		if(flags_->Flags::Vertex_Cut()){
			if(idx_ == 0){
				int run = flags_->Flags::Run();
				bool sim = flags_->Flags::Sim();
				float vz = data_->Branches::vz(idx_);
				pass = cuts::vertex_cut(vz,run);
			}
		}else{
			pass = true;
		}
		return pass;
	}