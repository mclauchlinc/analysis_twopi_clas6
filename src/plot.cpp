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
	//PID
	plot::plot_pid_cut(particle_, hist_, flags_, 3);
}

void plot::plot_no_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(par_ == 0){
		//hist_->Fill_SF();
		//hist_->Fill_CC();
		//hist_->Fill_EC();
		hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_none_,_cut_[0],_mnone_,false,flags_,particle_.Particle::Get_Weight());
	}
	hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_none_,_sector_[particle_.Particle::Sector()-1],_cut_[0],_mnone_,flags_,particle_.Particle::Get_Weight());
	//hist_->Fill_Delta();
}
void plot::plot_sanity_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(particle_.Pass_Sanity(par_)){
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_sanity_,_cut_[0],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_sanity_,_sector_[particle_.Particle::Sector()-1],_cut_[0],_mnone_,flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}else{
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_sanity_,_cut_[1],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_sanity_,_sector_[particle_.Particle::Sector()-1],_cut_[1],_mnone_,flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}
}
void plot::plot_fid_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(particle_.Pass_fid(par_)){
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_fid_cut_,_cut_[0],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_fid_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[0],_mnone_,flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}else{
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_fid_cut_,_cut_[1],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_fid_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[1],_mnone_,flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}
}
void plot::plot_delta_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(particle_.Pass_dt(par_)){
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_delta_cut_,_cut_[0],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_delta_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[0],_mnone_,flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}else{
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_delta_cut_,_cut_[1],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_delta_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[1],_mnone_,flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}
}
void plot::plot_sf_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(par_==0){
		if(particle_.Pass_sf()){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_sf_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[0],_mnone_,flags_,particle_.Particle::Get_Weight());
			//hist_->Fill_Delta();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_sf_cut_,_cut_[0],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}else{
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_sf_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[1],_mnone_,flags_,particle_.Particle::Get_Weight());
			//hist_->Fill_Delta();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_sf_cut_,_cut_[1],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}
	}
}
void plot::plot_cc_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(par_==0){
		if(particle_.Pass_sf()){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_cc_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[0],_mnone_,flags_,particle_.Particle::Get_Weight());
			//hist_->Fill_Delta();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_cc_cut_,_cut_[0],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}else{
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_cc_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[1],_mnone_,flags_,particle_.Particle::Get_Weight());
			//hist_->Fill_Delta();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_cc_cut_,_cut_[1],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}
	}
}
void plot::plot_ec_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(par_==0){
		if(particle_.Pass_sf()){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_ec_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[0],_mnone_,flags_,particle_.Particle::Get_Weight());
			//hist_->Fill_Delta();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_ec_cut_,_cut_[0],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}else{
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),physics::phi_center(particle_.Particle::Get_phi()),_ec_cut_,_sector_[particle_.Particle::Sector()-1],_cut_[1],_mnone_,flags_,particle_.Particle::Get_Weight());
			//hist_->Fill_Delta();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_ec_cut_,_cut_[1],_mnone_,false,flags_,particle_.Particle::Get_Weight());
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
		hist_->Histogram::Fid_Fill(par_,particle_.Particle::Get_theta(),particle_.Particle::Get_phi(),_beta_cut_,_sector_[particle_.Particle::Sector()],_cut_[0],fun::top_idx(_none_),flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}else{
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
		}
		hist_->Histogram::Fid_Fill(par_,particle_.Particle::Get_theta(),particle_.Particle::Get_phi(),_beta_cut_,_sector_[particle_.Particle::Sector()],_cut_[1],fun::top_idx(_none_),flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}
}*/
void plot::plot_pid_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_){
	if(particle_.Pass_pid(par_)){
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_pid_,_cut_[0],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),particle_.Particle::Get_phi(),_pid_,_sector_[particle_.Particle::Sector()-1],_cut_[0],_mnone_,flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}else{
		if(par_ == 0){
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
			hist_->Histogram::WQ2_Fill(particle_.Particle::W(),particle_.Particle::Q2(),_pid_,_cut_[1],_mnone_,false,flags_,particle_.Particle::Get_Weight());
		}
		hist_->Histogram::Fid_Fill(_species_[par_],particle_.Particle::Get_theta(),particle_.Particle::Get_phi(),_pid_,_sector_[particle_.Particle::Sector()-1],_cut_[1],_mnone_,flags_,particle_.Particle::Get_Weight());
		//hist_->Fill_Delta();
	}
}





void plot::plot_event(Event event_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, bool thrown_){
	if(thrown_){

	}else{
		//std::cout<<"Plotting Event: " <<_top_[event_.Event::Top()] <<"\n";
		for(int i=0; i<4; i++){//Topology
			if(event_.Event::Top()==i){
				if(event_.Event::Pass()){
					//PID Plots with Event Selection
					for(int j=0; j<4; j++){//Particle
						if(j == 0){//Electron
							//hist_->Fill_SF();
							//hist_->Fill_CC();
							//hist_->Fill_EC();
						}
						hist_->Histogram::Fid_Fill(_species_[j],event_.Event::Theta(j),event_.Event::Phi(j),_event_,_sector_[event_.Event::Sector(j)-1],_cut_[0],_top_[i],flags_,event_.Event::Weight());
						//hist_->Fill_Delta();
					}
					hist_->Histogram::MM_Fill(_top_[i],_cut_applied_,_clean_,event_.Event::MM2(),event_.Event::W(),event_.Event::Weight(),flags_);
					//hist_->Fill_MM(i);
					hist_->Histogram::WQ2_Fill(event_.Event::W(),event_.Event::Q2(),_event_,_cut_[0],_top_[event_.Event::Top()],thrown_,flags_,event_.Event::Weight());
				}else{
					//PID Plots with Event Selection
					for(int j=0; j<4; j++){//Particle
						if(j == 0){//Electron
							//hist_->Fill_SF();
							//hist_->Fill_CC();
							//hist_->Fill_EC();
						}
						
						hist_->Histogram::Fid_Fill(_species_[j],event_.Event::Theta(j),event_.Event::Phi(j),_event_,_sector_[event_.Event::Sector(j)-1],_cut_[1],_top_[i],flags_,event_.Event::Weight());
						//hist_->Fill_Delta();
					}
					hist_->Histogram::MM_Fill(_top_[i],_anti_cut_,_clean_,event_.Event::MM2(),event_.Event::W(),event_.Event::Weight(),flags_);
					//hist_->Fill_MM(i);
					hist_->Histogram::WQ2_Fill(event_.Event::W(),event_.Event::Q2(),_event_,_cut_[1],_top_[event_.Event::Top()],thrown_,flags_,event_.Event::Weight());
				}
				hist_->Histogram::MM_Fill(_top_[i],_no_cut_,_clean_,event_.Event::MM2(),event_.Event::W(),event_.Event::Weight(),flags_);
			}
		}
	}
}

void plot::plot_clean_event(Event event_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	int top = event_.Event::Top();
	//std::cout<<"Plotting Clean Event: " <<_top_[event_.Event::Top()] <<"\n";
	for(int j=0; j<4; j++){//Particle
		if(j == 0){//Electron
			//hist_->Fill_SF();
			//hist_->Fill_CC();
			//hist_->Fill_EC();
		}
		//Do need to modify the plots here
		//hist_->Histogram::Fid_Fill(j,event_.Event::Theta(j),event_.Event::Theta(j),_event_,_sector_[event_.Event::Sector(j)-1],_cut_[0],top,flags_,event_.Event::Weight());
		//hist_->Fill_Delta();
	}
	if(event_.Event::Pass()){
		hist_->Histogram::MM_Fill(_top_[top],_cut_applied_,_clean_,event_.Event::MM(),event_.Event::W(),event_.Event::Weight(),flags_);
	}else{
		hist_->Histogram::MM_Fill(_top_[top],_anti_cut_,_clean_,event_.Event::MM(),event_.Event::W(),event_.Event::Weight(),flags_);
	}
	hist_->Histogram::MM_Fill(_top_[top],_no_cut_,_clean_,event_.Event::MM(),event_.Event::W(),event_.Event::Weight(),flags_);
	//hist_->Fill_MM(top);
	//hist_->Histogram::WQ2_Fill(event_.Event::W(),event_.Event::Q2(),_event_,_cut_[0],event_.Event::Top(),false,flags_,event_.Event::Weight());
}

void plot::plot_isolated_event(Event event_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//std::cout<<"\tPlotting Isolated event " <<event_.Event::Top() <<" " <<event_.Event::Pass() <<"\n";
	int top_idx = -1; 
	for(int i=0; i<4; i++){//Topology
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
				//hist_->Histogram::Fid_Fill(j,event_.Event::Get_Theta(j),event_.Event::Get_Phi(j),_event_,_sector_[event_.Event::Sector(j)-1],_cut_[0],top,flags_,event_.Event::Get_Weight());
				//hist_->Fill_Delta();
			}
			//hist_->Fill_MM(i);
			//std::cout<<"\ttop: " <<_top_[i] <<"  cut: " <<_cut_applied_ <<"   MM^2: " <<event_.Event::MM2() <<"   W: " <<event_.Event::W() <<"\n";
			//hist_->Histogram::MM_Fill(_top_[i],_cut_applied_,_isolated_,event_.Event::MM(),event_.Event::W(),event_.Event::Weight(),flags_);
			hist_->Histogram::WQ2_Fill(event_.Event::W(),event_.Event::Q2(),_event_,_cut_applied_,_mall_,false,flags_,event_.Event::Weight());
		}
	}
	
}