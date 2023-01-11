#include "plot.hpp"

void plot::plot_pid(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	if(particle_.Get_idx() == 0){
		plot::plot_ele(particle_, hist_, flags_);
	}else{
		//std::cout<<"q: " <<particle_.Particle::Get_q() <<"\n";
		if(particle_.Particle::Get_q() > 0){
			plot::plot_pro(particle_,hist_,flags_);
			plot::plot_pip(particle_,hist_,flags_);
		}else if(particle_.Particle::Get_q() < 0){
			plot::plot_pim(particle_,hist_,flags_);
		}
	}
}
void plot::plot_ele(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//No Cuts
	plot::plot_no_cut(particle_, hist_, flags_, 0);
	//Sanity
	plot::plot_sanity_cut(particle_, hist_, flags_, 0);
	//Vertex
	plot::plot_vertex_cut(particle_, hist_, flags_, 0);
	//Fid
	plot::plot_fid_cut(particle_, hist_, flags_, 0);
	//Delta T
		//plot::plot_delta_cut(particle_, hist_, flags_, 0);
	//SF
	plot::plot_sf_cut(particle_, hist_, flags_, 0);
	//CC
	plot::plot_cc_cut(particle_, hist_, flags_, 0);
	//EC
	plot::plot_ec_cut(particle_, hist_, flags_, 0);
	//ID Bank
	plot::plot_id_cut(particle_, hist_, flags_, 0);
	//PID
	plot::plot_pid_cut(particle_, hist_, flags_, 0);
}
void plot::plot_pro(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//No Cuts
	plot::plot_no_cut(particle_, hist_, flags_, 1);
	//Sanity
	plot::plot_sanity_cut(particle_, hist_, flags_, 1);
	//Fid
	plot::plot_fid_cut(particle_, hist_, flags_, 1);
	//Delta T
	plot::plot_delta_cut(particle_, hist_, flags_, 1);
	//ID Bank
	plot::plot_id_cut(particle_, hist_, flags_, 1);
	//PID
	plot::plot_pid_cut(particle_, hist_, flags_, 1);
}
void plot::plot_pip(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//No Cuts
	plot::plot_no_cut(particle_, hist_, flags_, 2);
	//Sanity
	plot::plot_sanity_cut(particle_, hist_, flags_, 2);
	//Fid
	plot::plot_fid_cut(particle_, hist_, flags_, 2);
	//Delta T
	plot::plot_delta_cut(particle_, hist_, flags_, 2);
	//ID Bank
	plot::plot_id_cut(particle_, hist_, flags_, 2);
	//PID
	plot::plot_pid_cut(particle_, hist_, flags_, 2);
}
void plot::plot_pim(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//std::cout<<"Plotting Pim!\n";
	//No Cuts
	plot::plot_no_cut(particle_, hist_, flags_, 3);
	//Sanity
	plot::plot_sanity_cut(particle_, hist_, flags_, 3);
	//Fid
	plot::plot_fid_cut(particle_, hist_, flags_, 3);
	//Delta T
	plot::plot_delta_cut(particle_, hist_, flags_, 3);
	//ID Bank
	plot::plot_id_cut(particle_, hist_, flags_, 3);
	//PID
	plot::plot_pid_cut(particle_, hist_, flags_, 3);
}

void plot::plot_thrown(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//std::cout<<"Trying to plot a thrown particle\n";
	if(flags_->Flags::Sim() && particle_.Particle::Is_Thrown()){
		if(particle_.Particle::Is_Elec()){
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_event_,_cut_applied_,_mzero_,_thrown_,flags_,particle_.Particle::Get_Weight());
		}
	}
}


void plot::plot_no_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(par_ == 0){
		//std::cout<<"No Cut SF\n";
		hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _none_, _no_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
		//std::cout<<"No Cut CC\n";
		hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _none_, _no_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
		//hist_->Fill_EC();
		//std::cout<<"No Cut Vertex\n";
		hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _none_, _no_cut_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
		//std::cout<<"No Cut WQ2\n";
		hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_none_,_no_cut_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
		//std::cout<<"No Cut CC Eff\n";
		hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _none_, _no_cut_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
	}
	hist_->Histogram::SC_Eff_Fill(particle_.Particle::Get_sc_pd(), particle_.Particle::Get_Weight(), _species_[par_], _none_, _sector_[particle_.Particle::Sector()-1], flags_);
	//std::cout<<"No Cut Fid\n";
	hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_none_,_sector_[particle_.Particle::Sector()-1],_no_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
	//std::cout<<"No Cut Delta\n";
	hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _none_, _no_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
}
void plot::plot_sanity_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(particle_.Pass_Sanity(par_)){
		if(par_ == 0){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _sanity_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _sanity_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _sanity_, _cut_applied_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_sanity_,_cut_applied_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _sanity_, _cut_applied_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
		hist_->Histogram::SC_Eff_Fill(particle_.Particle::Get_sc_pd(), particle_.Particle::Get_Weight(), _species_[par_], _sanity_, _sector_[particle_.Particle::Sector()-1], flags_);
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_sanity_,_sector_[particle_.Particle::Sector()-1],_cut_[0],_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
		hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _sanity_, _cut_applied_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
	}else{
		if(par_ == 0){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _sanity_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _sanity_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _sanity_, _anti_cut_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_sanity_,_anti_cut_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _sanity_, _anti_cut_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_sanity_,_sector_[particle_.Particle::Sector()-1],_anti_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
		hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _sanity_, _anti_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
	}
}
void plot::plot_vertex_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(particle_.Pass_Sanity(par_)){
		if(par_ == 0){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _vertex_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _vertex_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _vertex_cut_, _cut_applied_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_vertex_cut_,_cut_applied_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _vertex_cut_, _cut_applied_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
		hist_->Histogram::SC_Eff_Fill(particle_.Particle::Get_sc_pd(), particle_.Particle::Get_Weight(), _species_[par_], _vertex_cut_, _sector_[particle_.Particle::Sector()-1], flags_);
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_vertex_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[0],_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
		hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _vertex_cut_, _cut_applied_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
	}else{
		if(par_ == 0){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _vertex_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _vertex_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _vertex_cut_, _anti_cut_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_vertex_cut_,_anti_cut_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _vertex_cut_, _anti_cut_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_vertex_cut_,_sector_[particle_.Particle::Sector()-1],_anti_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
		hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _vertex_cut_, _anti_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
	}
}
void plot::plot_fid_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(particle_.Pass_fid(par_)){
		if(par_ == 0){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _fid_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _fid_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _fid_cut_, _cut_applied_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_fid_cut_,_cut_applied_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _fid_cut_, _cut_applied_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
		hist_->Histogram::SC_Eff_Fill(particle_.Particle::Get_sc_pd(), particle_.Particle::Get_Weight(), _species_[par_], _fid_cut_, _sector_[particle_.Particle::Sector()-1], flags_);
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_fid_cut_,_sector_[particle_.Particle::Sector()-1],_cut_applied_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
		hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _fid_cut_, _cut_applied_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
	}else{
		if(par_ == 0){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _fid_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _fid_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _fid_cut_, _anti_cut_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_fid_cut_,_anti_cut_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _fid_cut_, _anti_cut_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_fid_cut_,_sector_[particle_.Particle::Sector()-1],_anti_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
		hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _fid_cut_, _anti_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
	}
}
void plot::plot_delta_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(particle_.Pass_dt(par_)){
		if(par_ == 0){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _delta_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _delta_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _delta_cut_, _cut_applied_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_delta_cut_,_cut_applied_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _delta_cut_, _cut_applied_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
		hist_->Histogram::SC_Eff_Fill(particle_.Particle::Get_sc_pd(), particle_.Particle::Get_Weight(), _species_[par_], _delta_cut_, _sector_[particle_.Particle::Sector()-1], flags_);
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_delta_cut_,_sector_[particle_.Particle::Sector()-1],_cut_applied_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
		hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _delta_cut_, _cut_applied_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
	}else{
		if(par_ == 0){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _delta_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _delta_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _delta_cut_, _anti_cut_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_delta_cut_,_anti_cut_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _delta_cut_, _anti_cut_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_delta_cut_,_sector_[particle_.Particle::Sector()-1],_anti_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
		hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _delta_cut_, _anti_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
	}
}
void plot::plot_sf_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(par_==0){
		if(particle_.Pass_sf()){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _sf_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _sf_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _sf_cut_, _cut_applied_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_sf_cut_,_sector_[particle_.Particle::Sector()-1],_cut_applied_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _sf_cut_, _cut_applied_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_sf_cut_,_cut_applied_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _sf_cut_, _cut_applied_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
			hist_->Histogram::SC_Eff_Fill(particle_.Particle::Get_sc_pd(), particle_.Particle::Get_Weight(), _species_[par_], _sf_cut_, _sector_[particle_.Particle::Sector()-1], flags_);
		}else{
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _sf_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _sf_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _sf_cut_, _anti_cut_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_sf_cut_,_sector_[particle_.Particle::Sector()-1],_anti_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _sf_cut_, _anti_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_sf_cut_,_anti_cut_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _sf_cut_, _anti_cut_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
	}
}
void plot::plot_cc_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(par_==0){
		if(particle_.Pass_sf()){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _cc_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _cc_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _cc_cut_, _cut_applied_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_cc_cut_,_sector_[particle_.Particle::Sector()-1],_cut_applied_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _cc_cut_, _cut_applied_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_cc_cut_,_cut_applied_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _cc_cut_, _cut_applied_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
			hist_->Histogram::SC_Eff_Fill(particle_.Particle::Get_sc_pd(), particle_.Particle::Get_Weight(), _species_[par_], _cc_cut_, _sector_[particle_.Particle::Sector()-1], flags_);
		}else{
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _cc_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _cc_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _cc_cut_, _anti_cut_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_cc_cut_,_sector_[particle_.Particle::Sector()-1],_anti_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _cc_cut_, _anti_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_cc_cut_,_anti_cut_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _cc_cut_, _anti_cut_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
	}
}
void plot::plot_ec_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(par_==0){
		if(particle_.Pass_ec()){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _ec_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _ec_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _ec_cut_, _cut_applied_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_ec_cut_,_sector_[particle_.Particle::Sector()-1],_cut_applied_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _ec_cut_, _cut_applied_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_ec_cut_,_cut_applied_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _ec_cut_, _cut_applied_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
			hist_->Histogram::SC_Eff_Fill(particle_.Particle::Get_sc_pd(), particle_.Particle::Get_Weight(), _species_[par_], _ec_cut_, _sector_[particle_.Particle::Sector()-1], flags_);
		}else{
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _ec_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _ec_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _ec_cut_, _anti_cut_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_ec_cut_,_sector_[particle_.Particle::Sector()-1],_anti_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _ec_cut_, _anti_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_ec_cut_,_anti_cut_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _ec_cut_, _anti_cut_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
	}
}
/*void plot::plot_beta_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(particle_.Pass_beta(par_)){
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
		}
		hist_->Histogram::Fid_Fill(par_,particle_.Particle::Get_theta(),particle_.Particle::Get_phi(),_beta_cut_,_sector_[particle_.Particle::Sector()],_cut_applied_,fun::top_idx(_none_),flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}else{
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
		}
		hist_->Histogram::Fid_Fill(par_,particle_.Particle::Get_theta(),particle_.Particle::Get_phi(),_beta_cut_,_sector_[particle_.Particle::Sector()],_anti_cut_,fun::top_idx(_none_),flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}
}*/
void plot::plot_id_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(!flags_->Flags::ID_Cut()){
		return;
	}
	if(par_==0){
		if(particle_.Pass_id(par_)){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _id_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _id_cut_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _id_cut_, _cut_applied_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_id_cut_,_sector_[particle_.Particle::Sector()-1],_cut_applied_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _id_cut_, _cut_applied_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_id_cut_,_cut_applied_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _id_cut_, _cut_applied_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
			hist_->Histogram::SC_Eff_Fill(particle_.Particle::Get_sc_pd(), particle_.Particle::Get_Weight(), _species_[par_], _id_cut_, _sector_[particle_.Particle::Sector()-1], flags_);
		}else{
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _id_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _id_cut_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _id_cut_, _anti_cut_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_id_cut_,_sector_[particle_.Particle::Sector()-1],_anti_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _id_cut_, _anti_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_id_cut_,_anti_cut_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _id_cut_, _anti_cut_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
	}else{
		if(particle_.Pass_id(par_)){
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_id_cut_,_sector_[particle_.Particle::Sector()-1],_cut_applied_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _id_cut_, _cut_applied_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
		}else{
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_id_cut_,_sector_[particle_.Particle::Sector()-1],_anti_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _id_cut_, _anti_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
		}
	}
}
void plot::plot_pid_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(particle_.Pass_pid(par_)){
		if(par_ == 0){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _pid_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _pid_, _cut_applied_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _pid_, _cut_applied_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_pid_,_cut_applied_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _pid_, _cut_applied_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
		hist_->Histogram::SC_Eff_Fill(particle_.Particle::Get_sc_pd(), particle_.Particle::Get_Weight(), _species_[par_], _pid_, _sector_[particle_.Particle::Sector()-1], flags_);
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_pid_,_sector_[particle_.Particle::Sector()-1],_cut_applied_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
		hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _pid_, _cut_applied_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
	}else{
		if(par_ == 0){
			hist_->Histogram::SF_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_sf(), particle_.Particle::W(), _pid_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], particle_.Particle::Get_Weight(), flags_);
			hist_->Histogram::CC_Fill(particle_.Particle::Get_nphe(), particle_.Particle::Get_cc_seg(), _pid_, _anti_cut_, _mnone_, _sector_[particle_.Particle::Sector()-1], _cc_sides_[particle_.Particle::Get_cc_lrc()], flags_);
			//hist_->Fill_EC();
			hist_->Histogram::Vertex_Fill(particle_.Particle::Get_vz(), particle_.Particle::Get_Weight(), _pid_, _anti_cut_,_mnone_,_sector_[particle_.Particle::Sector()-1],flags_);
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_pid_,_anti_cut_,_mnone_,_nthrown_,flags_,particle_.Particle::Get_Weight());
			hist_->Histogram::CC_Eff_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_theta(), particle_.Particle::Get_Weight(), _pid_, _anti_cut_, _sector_[particle_.Particle::Sector()-1], _mnone_, flags_);
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_pid_,_sector_[particle_.Particle::Sector()-1],_anti_cut_,_mnone_,particle_.Particle::Get_p(),flags_,particle_.Particle::Get_Weight());
		hist_->Histogram::Delta_Fill(particle_.Particle::Get_p(), particle_.Particle::Get_delta(par_), particle_.Particle::Get_Weight(), particle_.Particle::W(),_species_[par_], _pid_, _anti_cut_, _mnone_,_sector_[particle_.Particle::Sector()-1], flags_ );
	}
}





void plot::plot_event(Event event_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, bool thrown_){
	//std::cout<<"Plotting Event \n";
	if(thrown_){
		hist_->Histogram::WQ2_Fill(event_.Event::W(),event_.Event::Q2(),_event_,_cut_applied_,_mzero_,_thrown_,flags_,event_.Event::Weight());
		for(int j=0; j<3; j++){
			if(flags_->Sim()){
				hist_->Histogram::Friend_Fill(_mzero_, event_.Event::W(), event_.Event::Q2(), event_.Event::MMb(j), event_.Event::MM2b(j), event_.Event::Thetab(j), event_.Event::Alphab(j), event_.Event::Phib(j) , j, thrown_, event_.Event::Weight(), event_.Event::Helicity(), event_.Event::Virtual_Photon_Flux(), flags_);
			}else{
				hist_->Histogram::Friend_Fill(_mzero_, event_.Event::W(), event_.Event::Q2(), event_.Event::MMb(j), event_.Event::MM2b(j), event_.Event::Thetab(j), event_.Event::Alphab(j), event_.Event::Phib(j) , j, thrown_, event_.Event::Weight(), event_.Event::Helicity(), event_.Event::CC_eff()*event_.Event::Virtual_Photon_Flux(), flags_);
			}
			
		}
	}else{
		//if(event_.Event::Top()==3){
			//std::cout<<"Plotting Event: " <<_top_[event_.Event::Top()] <<" w/ MM^2: " <<event_.Event::MM2() <<" Pass: " <<_truth_[fun::truth_idx(event_.Event::Pass())] <<"\n";// <<"\n";
		//}
		if(event_.Event::Pass()){
			//std::cout<<"\tPassed\n";
			//PID Plots with Event Selection
			for(int j=0; j<4; j++){//Particle
				if(j == 0){//Electron
					//hist_->Fill_SF();
					//hist_->Fill_CC();
					//hist_->Fill_EC();
				}
				//hist_->Histogram::Fid_Fill(_species_[j],event_.Event::Theta(j),event_.Event::Phi(j),_event_,_sector_[event_.Event::Sector(j)-1],_cut_applied_,_top_[event_.Event::Top()],flags_,event_.Event::Weight());
				//hist_->Fill_Delta();
			}
			hist_->Histogram::MM_Fill(_top_[event_.Event::Top()],_cut_applied_,_dirty_,_sector_[event_.Event::Sector(0)-1],event_.Event::MM2(),event_.Event::W(),event_.Event::Weight(),flags_);
			//hist_->Fill_MM(i);
			hist_->Histogram::WQ2_Fill(event_.Event::W(),event_.Event::Q2(),_event_,_cut_applied_,_top_[event_.Event::Top()],_nthrown_,flags_,event_.Event::Weight());
		}else{
			//std::cout<<"\tFailed\n";
			//PID Plots with Event Selection
			for(int j=0; j<4; j++){//Particle
				if(j == 0){//Electron
					//hist_->Fill_SF();
					//hist_->Fill_CC();
					//hist_->Fill_EC();
				}
				
				//hist_->Histogram::Fid_Fill(_species_[j],event_.Event::Theta(j),event_.Event::Phi(j),_event_,_sector_[event_.Event::Sector(j)-1],_anti_cut_,_top_[event_.Event::Top()],flags_,event_.Event::Weight());
				//hist_->Fill_Delta();
			}
			hist_->Histogram::MM_Fill(_top_[event_.Event::Top()],_anti_cut_,_dirty_,_sector_[event_.Event::Sector(0)-1],event_.Event::MM2(),event_.Event::W(),event_.Event::Weight(),flags_);
			//hist_->Fill_MM(i);
			hist_->Histogram::WQ2_Fill(event_.Event::W(),event_.Event::Q2(),_event_,_anti_cut_,_top_[event_.Event::Top()],_nthrown_,flags_,event_.Event::Weight());
		}
		hist_->Histogram::MM_Fill(_top_[event_.Event::Top()],_no_cut_,_dirty_,_sector_[event_.Event::Sector(0)-1],event_.Event::MM2(),event_.Event::W(),event_.Event::Weight(),flags_);
	}
}


void plot::plot_clean_event(Event event_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	int top = event_.Event::Top();
	//std::cout<<"Plotting Clean Event: " <<_top_[event_.Event::Top()] <<"\n";
	if(event_.Event::Pass()){
		//std::cout<<"Event Passed \n";
		for(int k=0; k<3; k++){
			//std::cout<<"Friend Fill \n";
			if(flags_->Flags::Sim()){
				hist_->Histogram::Friend_Fill(_top_[event_.Event::Top()], event_.Event::W(), event_.Event::Q2(), event_.Event::MMb(k), event_.Event::MM2b(k), event_.Event::Thetab(k), event_.Event::Alphab(k), event_.Event::Phib(k) , k, false, event_.Event::Weight(), event_.Event::Helicity(), event_.Event::CC_eff(), flags_);
			}else{
				hist_->Histogram::Friend_Fill(_top_[event_.Event::Top()], event_.Event::W(), event_.Event::Q2(), event_.Event::MMb(k), event_.Event::MM2b(k), event_.Event::Thetab(k), event_.Event::Alphab(k), event_.Event::Phib(k) , k, false, event_.Event::Weight(), event_.Event::Helicity(), event_.Event::CC_eff()*event_.Event::Virtual_Photon_Flux(), flags_);
			}
		}
		//std::cout<<"particle Filling\n";
		for(int j=0; j<4; j++){//Particle
			//std::cout<<"Filling particle " <<_species_[j] <<"\n";
			if(j == 0){//Electron
				hist_->Histogram::SF_Fill(event_.Event::P(j), event_.Event::SF(), event_.Event::W(), _event_, _cut_applied_, _top_[event_.Event::Top()], _sector_[event_.Event::Sector(j)-1], event_.Event::Weight(), flags_);
				hist_->Histogram::CC_Fill(event_.Event::nphe(), event_.Event::CC_seg(), _event_, _cut_applied_, _top_[top], _sector_[event_.Event::Sector(j)-1], _cc_sides_[event_.Event::CC_side()], flags_);
				//hist_->Fill_EC();
				hist_->Histogram::Vertex_Fill(event_.Event::Vz(), event_.Event::Weight(), _event_, _cut_applied_, _top_[event_.Event::Top()], _sector_[event_.Event::Sector(j)-1],flags_);
				hist_->Histogram::CC_Eff_Fill(event_.Event::P(j), event_.Event::Theta(j), event_.Event::Weight(), _event_, _cut_applied_, _sector_[event_.Event::Sector(j)-1], _top_[event_.Event::Top()], flags_);
			}
			//Do need to modify the plots here
			hist_->Histogram::Fid_Fill(_species_[j],event_.Event::Theta(j),physics::phi_center(event_.Event::Phi(j)),_event_,_sector_[event_.Event::Sector(j)-1],_cut_applied_,_top_[top],event_.Event::P(j),flags_,event_.Event::Weight());
			//hist_->Histogram::Fid_Fill(_species_[j],event_.Event::Theta(j),event_.Event::Theta(j),_event_,_sector_[event_.Event::Sector(j)-1],_cut_applied_,_top_[top],flags_,event_.Event::Weight());
			hist_->Histogram::Delta_Fill(event_.Event::P(j), event_.Event::Delta(j), event_.Event::Weight(), event_.Event::W(),_species_[j], _event_, _cut_applied_, _top_[top], _sector_[event_.Event::Sector(j)-1], flags_ );
		}
		hist_->Histogram::MM_Fill(_top_[top],_cut_applied_,_clean_,_sector_[event_.Event::Sector(0)-1],event_.Event::MM2(),event_.Event::W(),event_.Event::Weight(),flags_);
	}else{
		//std::cout<<"Event Failed \n";
		for(int j=0; j<4; j++){//Particle
			//std::cout<<"Filling particle " <<_species_[j] <<"\n";
			if(j == 0){//Electron
				//std::cout<<"Filling SF failure\n";
				hist_->Histogram::SF_Fill(event_.Event::P(j), event_.Event::SF(), event_.Event::W(), _event_, _anti_cut_, _top_[event_.Event::Top()], _sector_[event_.Event::Sector(j)-1], event_.Event::Weight(), flags_);
				//std::cout<<"Filling CC failure\n";
				hist_->Histogram::CC_Fill(event_.Event::nphe(), event_.Event::CC_seg(), _event_, _cut_applied_, _top_[top], _sector_[event_.Event::Sector(j)-1], _cc_sides_[event_.Event::CC_side()], flags_);
				//hist_->Fill_EC();
				//std::cout<<"Filling Vertex failure\n";
				hist_->Histogram::Vertex_Fill(event_.Event::Vz(), event_.Event::Weight(), _event_, _anti_cut_, _top_[event_.Event::Top()], _sector_[event_.Event::Sector(j)-1],flags_);
				//std::cout<<"Filling CC_Eff failure\n";
				hist_->Histogram::CC_Eff_Fill(event_.Event::P(j), event_.Event::Theta(j), event_.Event::Weight(), _event_, _anti_cut_, _sector_[event_.Event::Sector(j)-1], _top_[event_.Event::Top()], flags_);
			}
			//Do need to modify the plots here
			//std::cout<<"Filling Fid failure\n";
			hist_->Histogram::Fid_Fill(_species_[j],event_.Event::Theta(j),physics::phi_center(event_.Event::Phi(j)),_event_,_sector_[event_.Event::Sector(j)-1],_anti_cut_,_top_[top],event_.Event::P(j),flags_,event_.Event::Weight());
			//std::cout<<"Filling Delta failure\n";
			hist_->Histogram::Delta_Fill(event_.Event::P(j), event_.Event::Delta(j), event_.Event::Weight(), event_.Event::W(),_species_[j], _event_, _anti_cut_, _top_[top], _sector_[event_.Event::Sector(j)-1], flags_ );
		}
		//std::cout<<"Filling MM for failure\n";
		hist_->Histogram::MM_Fill(_top_[top],_anti_cut_,_clean_,_sector_[event_.Event::Sector(0)-1],event_.Event::MM2(),event_.Event::W(),event_.Event::Weight(),flags_);
	}
	//std::cout<<"Filling MM for no cut\n";
	hist_->Histogram::MM_Fill(_top_[top],_no_cut_,_clean_,_sector_[event_.Event::Sector(0)-1],event_.Event::MM2(),event_.Event::W(),event_.Event::Weight(),flags_);
	if(flags_->Plot_Check()){
		if(flags_->E_PCorr() && !flags_->Sim()){
			hist_->Histogram::PCorr_Check_Fill(event_.Event::MM2(),event_.Event::Sector(0),_top_[top],_ele_p_corr_,flags_);
		}
		if(flags_->E_Theta_Corr()&& !flags_->Sim()){
			hist_->Histogram::PCorr_Check_Fill(event_.Event::MM2(),event_.Event::Sector(0),_top_[top],_ele_angle_corr_,flags_);
		}else if(!flags_->Sim()){
			hist_->Histogram::PCorr_Check_Fill(event_.Event::MM2(),event_.Event::Sector(0),_top_[top],_no_corr_,flags_);
		}
	}
	
	//hist_->Histogram::WQ2_Fill(event_.Event::W(),event_.Event::Q2(),_event_,_cut_applied_,event_.Event::Top(),_nthrown_,flags_,event_.Event::Weight());
}

void plot::plot_isolated_event(Event event_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//std::cout<<"\tPlotting Isolated event " <<_top_[event_.Event::Top()] <<" " <<_truth_[event_.Event::Pass()] <<"\n";
	//int top_idx = -1; 
	hist_->Histogram::WQ2_Fill(event_.Event::W(),event_.Event::Q2(),_event_,_cut_applied_,_mall_,_nthrown_,flags_,event_.Event::Weight());
	if(event_.Event::Pass()){
		for(int j=0; j<3; j++){
			if(flags_->Flags::Sim()){
				hist_->Histogram::Friend_Fill(_mall_, event_.Event::W(), event_.Event::Q2(), event_.Event::MMb(j), event_.Event::MM2b(j), event_.Event::Thetab(j), event_.Event::Alphab(j), event_.Event::Phib(j) , j, false, event_.Event::Weight(), event_.Event::Helicity(), event_.Event::Virtual_Photon_Flux(), flags_);
			}else{
				hist_->Histogram::Friend_Fill(_mall_, event_.Event::W(), event_.Event::Q2(), event_.Event::MMb(j), event_.Event::MM2b(j), event_.Event::Thetab(j), event_.Event::Alphab(j), event_.Event::Phib(j) , j, false, event_.Event::Weight(), event_.Event::Helicity(), event_.Event::CC_eff()*event_.Event::Virtual_Photon_Flux(), flags_);
			}	
		}
		for(int i=0; i<4; i++){//Species Loop
			if(_species_[i]==_ele_){
				hist_->Histogram::Vertex_Fill(event_.Event::Vz(), event_.Event::Weight(), _event_, _cut_applied_, _mall_, _sector_[event_.Event::Sector(i)-1],flags_);
				hist_->Histogram::Fid_Fill(_species_[i],event_.Event::Theta(i),physics::phi_center(event_.Event::Phi(i)),_event_,_sector_[event_.Event::Sector(i)-1],_cut_applied_,_mall_,event_.Event::P(i),flags_,event_.Event::Weight());
				hist_->Histogram::SF_Fill(event_.Event::P(i), event_.Event::SF(), event_.Event::W(), _event_, _cut_applied_, _mall_, _sector_[event_.Event::Sector(i)-1], event_.Event::Weight(), flags_);
				hist_->Histogram::CC_Fill(event_.Event::nphe(), event_.Event::CC_seg(), _event_, _cut_applied_, _mall_, _sector_[event_.Event::Sector(i)-1], _cc_sides_[event_.Event::CC_side()], flags_);
				hist_->Histogram::CC_Eff_Fill(event_.Event::P(i), event_.Event::Theta(i), event_.Event::Weight(), _event_, _cut_applied_, _sector_[event_.Event::Sector(i)-1], _mall_, flags_);
				//std::cout<<"Filling SC_Eff " <<_species_[i] <<" event with " <<event_.Event::SC_pd(i) <<"\n";
				hist_->Histogram::SC_Eff_Fill(event_.Event::SC_pd(i), event_.Event::Weight(), _species_[i], _event_, _sector_[event_.Event::Sector(i)-1], flags_);
			}else{
				switch(event_.Event::Top()){
					case 0:
						if(_species_[i]!=_pro_){
							hist_->Histogram::Fid_Fill(_species_[i],event_.Event::Theta(i),physics::phi_center(event_.Event::Phi(i)),_event_,_sector_[event_.Event::Sector(i)-1],_cut_applied_,_mall_,event_.Event::P(i),flags_,event_.Event::Weight());
							hist_->Histogram::Delta_Fill(event_.Event::P(i), event_.Event::Delta(i), event_.Event::Weight(), event_.Event::W(),_species_[i], _event_, _cut_applied_, _mall_,_sector_[event_.Event::Sector(i)-1], flags_ );
							//std::cout<<"Filling SC_Eff " <<_species_[i] <<" event with " <<event_.Event::SC_pd(i) <<"\n";
							hist_->Histogram::SC_Eff_Fill(event_.Event::SC_pd(i), event_.Event::Weight(), _species_[i], _event_, _sector_[event_.Event::Sector(i)-1], flags_);
						}
					break;
					case 1:
						if(_species_[i]!=_pip_){
							hist_->Histogram::Fid_Fill(_species_[i],event_.Event::Theta(i),physics::phi_center(event_.Event::Phi(i)),_event_,_sector_[event_.Event::Sector(i)-1],_cut_applied_,_mall_,event_.Event::P(i),flags_,event_.Event::Weight());
							hist_->Histogram::Delta_Fill(event_.Event::P(i), event_.Event::Delta(i), event_.Event::Weight(), event_.Event::W(),_species_[i], _event_, _cut_applied_, _mall_,_sector_[event_.Event::Sector(i)-1], flags_ );
							//std::cout<<"Filling SC_Eff " <<_species_[i] <<" event with " <<event_.Event::SC_pd(i) <<"\n";
							hist_->Histogram::SC_Eff_Fill(event_.Event::SC_pd(i), event_.Event::Weight(), _species_[i], _event_, _sector_[event_.Event::Sector(i)-1], flags_);
						}
					break;
					case 2:
						if(_species_[i]!=_pim_){
							hist_->Histogram::Fid_Fill(_species_[i],event_.Event::Theta(i),physics::phi_center(event_.Event::Phi(i)),_event_,_sector_[event_.Event::Sector(i)-1],_cut_applied_,_mall_,event_.Event::P(i),flags_,event_.Event::Weight());
							hist_->Histogram::Delta_Fill(event_.Event::P(i), event_.Event::Delta(i), event_.Event::Weight(), event_.Event::W(),_species_[i], _event_, _cut_applied_, _mall_,_sector_[event_.Event::Sector(i)-1], flags_ );
							//std::cout<<"Filling SC_Eff " <<_species_[i] <<" event with " <<event_.Event::SC_pd(i) <<"\n";
							hist_->Histogram::SC_Eff_Fill(event_.Event::SC_pd(i), event_.Event::Weight(), _species_[i], _event_, _sector_[event_.Event::Sector(i)-1], flags_);
						}
					break;
					case 3:
						hist_->Histogram::Fid_Fill(_species_[i],event_.Event::Theta(i),physics::phi_center(event_.Event::Phi(i)),_event_,_sector_[event_.Event::Sector(i)-1],_cut_applied_,_mall_,event_.Event::P(i),flags_,event_.Event::Weight());
						hist_->Histogram::Delta_Fill(event_.Event::P(i), event_.Event::Delta(i), event_.Event::Weight(), event_.Event::W(),_species_[i], _event_, _cut_applied_, _mall_,_sector_[event_.Event::Sector(i)-1], flags_ );
						//std::cout<<"Filling SC_Eff " <<_species_[i] <<" event with " <<event_.Event::SC_pd(i) <<"\n";
						hist_->Histogram::SC_Eff_Fill(event_.Event::SC_pd(i), event_.Event::Weight(), _species_[i], _event_, _sector_[event_.Event::Sector(i)-1], flags_);
					break;
					default:
						std::cout<<"Plotting Isolated Event without valid Topology\t" <<event_.Event::Top() <<"\n";
					break;
				}
			}
		}
	}
	/*for(int i=0; i<4; i++){//Topology
		if(event_.Event::Top()==i && event_.Event::Pass()){
			//std::cout<<"found the top: ";
			//PID Plots with Event Selection
			for(int j=0; j<4; j++){//Particle
				if(j == 0){//Electron
					//hist_->Fill_SF();
					//hist_->Fill_CC();
					//hist_->Fill_EC();
				}
				//Need to modify this
				//hist_->Histogram::Fid_Fill(j,event_.Event::Get_Theta(j),event_.Event::Get_Phi(j),_event_,_sector_[event_.Event::Sector(j)-1],_cut_applied_,top,flags_,event_.Event::Get_Weight());
				//hist_->Fill_Delta();
			}
			//hist_->Fill_MM(i);
			std::cout<<"\ttop: " <<_top_[i] <<"  cut: " <<_cut_applied_ <<"   MM^2: " <<event_.Event::MM2() <<"   W: " <<event_.Event::W() <<"\n";
			//hist_->Histogram::MM_Fill(_top_[i],_cut_applied_,_isolated_,event_.Event::MM(),event_.Event::W(),event_.Event::Weight(),flags_);
			hist_->Histogram::WQ2_Fill(event_.Event::W(),event_.Event::Q2(),_event_,_cut_applied_,_mall_,_nthrown_,flags_,event_.Event::Weight());
		}
	}*/
	
}