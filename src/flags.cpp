#include "flags.hpp"


Flags::Flags(){}

void Flags::Help(){
	std::cout<<"Welcome to the Help Menu to assist in your usage of this \n";
	std::cout<<"ele = electron\npro = proton\npip = pi+\npim=pi-\n";
	std::cout<<"--------------\n";
	std::cout<<"	*Run Information* \n";
	std::cout<<"-t=*";
	std::cout<<"	e16 or e1f 	(run group)\n";
	std::cout<<"	filled 		(target filled)\n";
	std::cout<<"	sim      	(simulated data?)\n";
	std::cout<<"	hel			(Track beam Helicity)\n";
	std::cout<<"--------------\n";
	std::cout<<"	*Cut Information*			\n";
	std::cout<<"-cut_all		(perform all cuts)			\n";
	std::cout<<"-cut=*			(perform these cuts)			\n";
	std::cout<<"-ncut=*			(do not perform these cuts [needs -cut_all to be initiated first])\n";
	std::cout<<"	fid_ele		(fiducial cut electron)\n";
	std::cout<<"	fid_pro		(fiducial cut proton)\n";
	std::cout<<"	fid_pip		(fiducial cut pi+)\n";
	std::cout<<"	fid_pim		(fiducial cut pi-)\n";
	std::cout<<"	dt_ele		(delta t cut electron)\n";
	std::cout<<"	dt_pro		(delta t cut proton)\n";
	std::cout<<"	dt_pip		(delta t cut pi+)\n";
	std::cout<<"	dt_pim		(delta t cut pi-)\n";
	std::cout<<"	beta_ele	(beta cut electron)\n";
	std::cout<<"	beta_pro	(beta cut proton)\n";
	std::cout<<"	beta_pip	(beta cut pi+)\n";
	std::cout<<"	beta_pim	(beta cut pi-)\n";
	std::cout<<"	sf			(SF cut ele)\n";
	std::cout<<"	cc			(Min CC cut ele)\n";
	std::cout<<"	ec			(Min EC cut ele)\n";
	std::cout<<"	vert		(vertex cut ele)\n";
	std::cout<<"	~Event Selection~	\n";
	std::cout<<"	mm_pro		(Missing Pro)\n";
	std::cout<<"	mm_pip		(Missing Pip)\n";
	std::cout<<"	mm_pim		(Missing Pim)\n";
	std::cout<<"	mm_zero		(Exclusive)\n";
	std::cout<<"--------------\n";
	std::cout<<"	*Plotting Information*			\n";
	std::cout<<"-plot_all		(Make/Fill all Plots)			\n";
	std::cout<<"-plot=*			(Make/Fill these Plots)			\n";
	std::cout<<"-nplot=*		(do not Make/Fill these Plots [needs -plot_all to be initiated first])\n";
	std::cout<<"	wq2 		(W vs Q2)\n";
}

void Flags::Read_Flags(int argc, char** argv){
	if(argc > 1){
		for(int i=0; i<argc-1; i++){
			std::cout<<"Reading Flag: " <<argv[i+1] <<std::endl;
			Flags::ID_Flag(std::string(argv[i+1]));
		}
	}
}
void Flags::ID_Flag(std::string argi){
	std::string str=argi.c_str();
	int str_len = str.length();
	std::cout<<"IDing Flags:\n\t"<<"flag:" <<str <<" | len:" <<str_len <<std::endl;
	if(str_len > 3){
		//Check for Target
		//std::cout<<"Checking Length 3 "<<str.substr(0,3) <<"\n";
		if(str.substr(0,3) == _target_){
			//std::cout<<"\t\tInfo on the target\n";
			Flags::Run_Flag(str.substr(3,str_len));
			//std::cout<<"\t\t"<<str.substr(3,str_len) <<std::endl;
		}
		if(str.substr(0,3) == _num_files_){
			//std::cout<<"\t\tInfo on the number of files\n";
			_num_files = stoi(str.substr(3,str_len));
			//std::cout<<"\t\t"<<str.substr(3,str_len) <<std::endl;
		}
	}if(str_len > 5){
		//std::cout<<"Checking Length 5 "<<str.substr(0,5) <<"\n";
		if(str.substr(0,5) == _cut_f_){
			//std::cout<<"\t\tInfo on the cuts\n";
			Flags::Cut_Flag(str.substr(5,str_len),true);
			//std::cout<<"\t\t"<<str.substr(5,str_len) <<std::endl;
		}
		if(str.substr(0,5) == _eff_cut_){
			//std::cout<<"Info on the efficiencies\n";
			Flags::Eff_Flag(str.substr(5,str_len),true);
			//std::cout<<"\t"<<str.substr(5,str_len) <<std::endl;
		}
		if(str.substr(0,5) == _geo_cut_){
			//std::cout<<"Info on the Geo cuts\n";
			//Flags::Geo_Flag(str.substr(5,str_len),true);
			//std::cout<<"\t"<<str.substr(5,str_len) <<std::endl;
		}
		if(str.substr(0,5) == _location_){
			//std::cout<<"Info on the file locations\n";
			_location = str.substr(5,str_len);
			//std::cout<<"\t"<<str.substr(5,str_len) <<std::endl;
		}
	}if(str_len > 6){
		//std::cout<<"Checking Length 6 "<<str.substr(0,6) <<"\n";
		if(str.substr(0,6) == _ncut_){
			//std::cout<<"Info on the what not to cut\n";
			Flags::Cut_Flag(str.substr(6,str_len),false);
		}else if(str.substr(0,6) == _neff_cut_){
			//std::cout<<"Info on the what efficiencies to NOT do\n";
			Flags::Eff_Flag(str.substr(6,str_len),false);
		}else if(str.substr(0,6) == _ngeo_cut_){
			//std::cout<<"Info on the Not Geo \n";
			//Flags::Geo_Flag(str.substr(6,str_len),false);
		}else if(str.substr(0,6) == _corr_){
			//std::cout<<"Info on the Corrections\n";
			Flags::Corr_Flag(str.substr(6,str_len),true);
		}else if(str.substr(0,6) == _plot_){
			std::cout<<"Info on the plotting\n";
			Flags::Plot_Flag(str.substr(6,str_len),true);
		}else if(str.substr(0,6) == _output_name_){
			//std::cout<<"Info on the output name\n";
			_output_name = str.substr(6,str_len);
			std::cout<<"Output file named: " <<_output_name <<"\n";
		}else if(str.substr(0,6) == _base_){
			_file_base = str.substr(6,str_len);
		}
	}if(str_len > 7){
		//std::cout<<"Checking Length 7 "<<str.substr(0,7) <<"\n";
		if(str.substr(0,7) == _ncorr_){
			//std::cout<<"Info on the NOT Corrections\n";
			Flags::Corr_Flag(str.substr(7,str_len),false);
		}else if(str.substr(0,7) == _nplot_){
			//std::cout<<"Info on the NOT plot\n";
			Flags::Plot_Flag(str.substr(7,str_len),false);
		}else if(str.substr(0,7) == _make_image_){
			//std::cout<<"Info on the image name\n";
			_make_image = true;
			_image_name = str.substr(7,str_len);
		}else if(str.substr(0,7) == _num_cores_){
			_num_cores = stoi(str.substr(7,str_len));
		}else if(str.substr(0,7) == _loose_cut_){
			_loose_cut = stoi(str.substr(7,str_len));
		}else if(str.substr(0,7) == _tight_cut_){
			_tight_cut = stoi(str.substr(7,str_len));
		}
	}if(str_len > 8){
		if(str.substr(0,8) == _make_friend_){
			_make_friend = true;
			_friend_name = str.substr(8,str_len);
			std::cout<<"\tName friend: " <<_friend_name <<"\n";
		}
	}
}
void Flags::Cut_Flag(std::string str, bool include){
	if(str==_all_f_){
		std::cout<<"Will do all cuts\n";
		for(int j=0; j<4; j++){
			_fid_cut[j] = true;
			_delta_cut[j] = true;
			_beta_cut[j] = true;
			_mm_cut[j] = true;
			for(int i=0; i<3; i++){
				_fid_geo_cut[j][i] = true;
			}
			_kin_eff_cut[j] = true;
		}
		_sf_cut = true;
		_cc_cut = true;
		_ec_cut = true;
		_vertex_cut=true;
		_wq2_cut=true;
		_cut_all = true;

	}
	if(include && !_cut_all){
		for(int i=0; i<4; i++){
			if(str==_fid_cut_f_[i]){
				std::cout<<"Will do Fid Cut: " <<i <<std::endl;
				_fid_cut[i] = true;
			}
			if(str==_delta_cut_f_[i]){
				std::cout<<"Will do Delta Cut: " <<i <<std::endl;
				_delta_cut[i] = true;
			}
			if(str==_beta_cut_f_[i]){
				std::cout<<"Will do Beta Cut: " <<i <<std::endl;
				_beta_cut[i] = true;
			}
			if(str==_mm_cut_f_[i]){
				std::cout<<"Will do MM Cut: " <<i <<std::endl;
				_mm_cut[i]=true;
			}
			if(str==_kin_eff_cut_f_[i]){
				std::cout<<"Will do Kin Eff Cut: " <<i <<std::endl;
				_kin_eff_cut[i]=true;
			}
			for(int j=0; j<3; j++){
				if(str == _fid_geo_cut_f_[i][j]){
					std::cout<<"Will cut fid geo " <<i <<" " <<j <<"\n";
					_fid_geo_cut[i][j] = true;
				}
			}
		}
		if(str==_sf_cut_f_){
			std::cout<<"Will do SF Cut\n";
			_sf_cut=true;
		}
		if(str==_id_cut_f_){
			std::cout<<"Will do ID Cut\n";
			_id_cut=true;
		}
		if(str==_cc_cut_f_){
			std::cout<<"Will do CC Cut\n";
			_cc_cut=true;
		}
		if(str==_ec_cut_f_){
			std::cout<<"Will do EC Cut\n";
			_ec_cut=true;
		}
		if(str==_vertex_cut_f_){
			std::cout<<"Will do Vertex Cut\n";
			_vertex_cut=true;
		}
		if(str==_wq2_cut_f_){
			std::cout<<"Will do W Q^2 Cut\n";
			_wq2_cut=true;
		}
		
	}
	if(!include && _cut_all){
		for(int i=0; i<4; i++){
			if(str==_fid_cut_f_[i]){
				std::cout<<"Will NOT do Fid Cut: " <<i <<"\n";
				_fid_cut[i] = false;
			}
			if(str==_delta_cut_f_[i]){
				std::cout<<"Will NOT do Delta Cut: " <<i <<"\n";
				_delta_cut[i] = false;
			}
			if(str==_beta_cut_f_[i]){
				std::cout<<"Will NOT do Beta Cut: " <<i <<"\n";
				_beta_cut[i] = false;
			}
			if(str==_mm_cut_f_[i]){
				std::cout<<"Will NOT do MM Cut: " <<i <<"\n";
				_mm_cut[i]=false;
			}
			if(str==_kin_eff_cut_f_[i]){
				std::cout<<"Will NOT do Kin Eff Cut: " <<i <<std::endl;
				_kin_eff_cut[i]=false;
			}
			for(int j=0; j<3; j++){
				if(str == _fid_geo_cut_f_[i][j]){
					std::cout<<"Will NOT Cut fid geo " <<i <<" " <<j <<"\n";
					_fid_geo_cut[i][j] = false;
				}
			}
		}
		if(str==_sf_cut_f_){
			std::cout<<"Will NOT do SF Cut\n";
			_sf_cut=false;
		}
		if(str==_cc_cut_f_){
			std::cout<<"Will NOT do CC Cut\n";
			_cc_cut=false;
		}
		if(str==_ec_cut_f_){
			std::cout<<"Will NOT do EC Cut\n";
			_ec_cut=false;
		}
		if(str==_vertex_cut_f_){
			std::cout<<"Will NOT do Vertex Cut\n";
			_vertex_cut=false;
		}
		if(str==_wq2_cut_f_){
			std::cout<<"Will NOT do W Q^2 Cut\n";
			_wq2_cut=false;
		}
	}
}
void Flags::Run_Flag(std::string str){
	for(int i = 0; i<2; i++){
		if( str == _run_group_[i]){
			std::cout<<"Named Run Group: " <<i <<std::endl;
			_run_group = i;
		}	
	}
	if(str==_filled_){
		std::cout<<"Target Filled\n";
		_filled = true;
	}
	if( str == _sim_){
		std::cout<<"Simulated Events\n";
		_sim = true;
	}
	if( str == _helicity_){
		std::cout<<"Track Helicity\n";
		_helicity = true;
	}
}
void Flags::Plot_Flag(std::string str, bool include){
	std::cout<<"\tPlot ID: "<<str <<"\n";
	if(str==_all_f_){
		std::cout<<"Will do all plots\n";
		_plot_wq2 = true;
		for(int i=0; i<4; i++){
			_plot_fid[i] = true;
			_plot_dt[i] = true;
			_plot_mm[i] = true;
			_plot_beta[i] = true;
			for(int j=0; j<3; j++){
					_plot_fid_geo[i][j] = true;
			}
			_plot_kin_eff[i] = true;
		}
		_plot_sf = true;
		_plot_cc = true;
		_plot_cc_eff = true;
		_plot_sc_eff = true;
		_plot_ec_eff = true;
		_plot_dc_eff = true;
		_plot_vertex = true;
		_plot_e_pcorr = true;
		_plot_elastic = true;
		_plot_check = true;
		_plot_all = true;
	}
	if(include && !_plot_all){
		if(str == _plot_wq2_){
			std::cout<<"Info on WQ2 plot\n";
			_plot_wq2 = true;
		}
		for(int i=0; i<4; i++){
			if(str == _plot_fid_[i]){
				std::cout<<"Will Plot Fid:" <<i <<"\n";
				_plot_fid[i] = true;
			}
			if(str == _plot_dt_[i]){
				std::cout<<"Will Plot Delta T:" <<i <<"\n";
				_plot_dt[i] = true;
			}
			if(str == _plot_mm_[i]){
				std::cout<<"Will Plot MM:" <<i <<"\n";
				_plot_mm[i] = true;
			}
			if(str == _plot_beta_[i]){
				std::cout<<"Will Plot Beta:" <<i <<"\n";
				_plot_beta[i] = true;
			}
			if(str == _plot_kin_eff_[i]){
				std::cout<<"Will Plot Kinematic Efficiency:" <<i <<"\n";
				_plot_kin_eff[i] = true;
			}
			for(int j=0; j<3; j++){
				if(str == _plot_fid_geo_[i][j]){
					std::cout<<"Will Plot fid geo " <<i <<" " <<j <<"\n";
					_plot_fid_geo[i][j] = true;
				}
			}
		}
		if(str == _plot_sf_){
			std::cout<<"Will Plot SF\n";
			_plot_sf = true;
		}
		if(str == _plot_cc_){
			std::cout<<"Will Plot cc\n";
			_plot_cc = true;
		}
		if(str == _plot_cc_eff_){
			std::cout<<"Will Plot cc eff\n";
			_plot_cc_eff = true;
		}
		if(str == _plot_sc_eff_){
			std::cout<<"Will Plot sc eff\n";
			_plot_sc_eff = true;
		}
		if(str == _plot_ec_eff_){
			std::cout<<"Will Plot ec eff\n";
			_plot_ec_eff = true;
		}
		if(str == _plot_dc_eff_){
			std::cout<<"Will Plot dc eff\n";
			_plot_dc_eff = true;
		}
		if(str == _plot_vertex_){
			std::cout<<"Will Plot vertex\n";
			_plot_vertex = true;
		}
		if(str == _plot_e_pcorr_){
			std::cout<<"Will Plot e pcorr\n";
			_plot_e_pcorr = true;
		}
		if(str == _plot_elastic_){
			std::cout<<"Will Plot Elastic\n";
			_plot_elastic = true;
		}
		if(str == _plot_check_){
			std::cout<<"Will Plot plot check\n";
			_plot_check = true;
		}
		if(str == _plot_bin_centering_){
			std::cout<<"Will Plot Bin Centering\n";
			_plot_bin_centering = true;
		}
		if(str == _plot_bin_centering2_){
			std::cout<<"Will Plot Bin Centering2\n";
			_plot_bin_centering2 = true;
		}
		if(str == _plot_ecorr_angle_){
			std::cout<<"Will Plot Electron Angle Corrections\n";
			_plot_ecorr_angle = true;
		}
		if(str == _plot_ecorr_mag_){
			std::cout<<"Will Plot Electron Mag Corrections\n";
			_plot_ecorr_mag = true;
		}
	}
	if(!include && _plot_all){
		if(str == _plot_wq2_){
			_plot_wq2 = false;
		}
		for(int i=0; i<4; i++){
			if(str == _plot_fid_[i]){
				std::cout<<"Will NOT Plot Fid " <<i <<"\n";
				_plot_fid[i] = false;
			}
			if(str == _plot_dt_[i]){
				std::cout<<"Will NOT Plot DT " <<i <<"\n";
				_plot_dt[i] = false;
			}
			if(str == _plot_mm_[i]){
				std::cout<<"Will NOT Plot MM " <<i <<"\n";
				_plot_mm[i] = false;
			}
			if(str == _plot_beta_[i]){
				std::cout<<"Will NOT Plot Beta " <<i <<"\n";
				_plot_beta[i] = false;
			}
			if(str == _plot_kin_eff_[i]){
				std::cout<<"Will NOT Plot kin_eff " <<i <<"\n";
				_plot_kin_eff[i] = false;
			}
			for(int j=0; j<3; j++){
				if(str == _plot_fid_geo_[i][j]){
					std::cout<<"Will Not Plot fid geo " <<i <<" " <<j <<"\n";
					_plot_fid_geo[i][j] = false;
				}
			}
		}
		if(str == _plot_sf_){
			_plot_sf = false;
		}
		if(str == _plot_cc_){
			_plot_cc = false;
		}
		if(str == _plot_cc_eff_){
			_plot_cc_eff = false;
		}
		if(str == _plot_sc_eff_){
			_plot_sc_eff = false;
		}
		if(str == _plot_ec_eff_){
			_plot_ec_eff = false;
		}
		if(str == _plot_dc_eff_){
			_plot_dc_eff = false;
		}
		if(str == _plot_vertex_){
			_plot_vertex = false;
		}
		if(str == _plot_e_pcorr_){
			_plot_e_pcorr = false;
		}
		if(str == _plot_elastic_){
			_plot_elastic = false;
		}
		if(str == _plot_check_){
			_plot_check = false;
		}
		if(str == _plot_bin_centering_){
			_plot_bin_centering = false;
		}
		if(str == _plot_bin_centering2_){
			_plot_bin_centering2 = false;
		}
		if(str == _plot_ecorr_angle_){
			_plot_ecorr_angle = false;
		}
		if(str == _plot_ecorr_angle_){
			_plot_ecorr_angle = false;
		}
	}
}
void Flags::Eff_Flag(std::string str, bool include){
	std::cout<<"Eff ID: " <<str <<"\n";
	if(str==_all_f_){
		_cc_eff_cut = true;
		_sc_eff_cut = true;
		_ec_eff_cut = true;
		_eff_all = true;
	}
	if(include && _eff_all){
		if(str == _eff_cc_){
			_cc_eff_cut = true;
		}
		if(str == _eff_dc_){
			_dc_eff_cut = true;
		}
		if(str == _eff_ec_){
			_ec_eff_cut = true;
		}
		if(str == _eff_sc_){
			_sc_eff_cut = true;
		}
	}
	if(!include && _eff_all){
		if(str == _eff_cc_){
			_cc_eff_cut = false;
		}
		if(str == _eff_dc_){
			_dc_eff_cut = false;
		}
		if(str == _eff_ec_){
			_ec_eff_cut = false;
		}
		if(str == _eff_sc_){
			_sc_eff_cut = false;
		}
	}
}
void Flags::Geo_Flag(std::string str, bool include){
	std::cout<<"Geo ID: " <<str <<"\n";
	if(str==_all_f_){
		_cc_geo_cut = true;
		_sc_geo_cut = true;
		_ec_geo_cut = true;
		_geo_all = true;
	}
	if(include && !_geo_all){
		if(str == _eff_cc_){
			_cc_geo_cut = true;
		}
		if(str == _eff_ec_){
			_ec_geo_cut = true;
		}
		if(str == _eff_sc_){
			_sc_geo_cut = true;
		}
	}
	if(!include && _geo_all){
		if(str == _eff_cc_){
			_cc_geo_cut = false;
		}
		if(str == _eff_ec_){
			_ec_geo_cut = false;
		}
		if(str == _eff_sc_){
			_sc_geo_cut = false;
		}
	}
}
void Flags::Corr_Flag(std::string str, bool include){
	std::cout<<"Corr ID: " <<str <<"\n";
	if(str == _all_f_){
		_p_corr = true;
		_en_corr = true;
		_v_corr = true;
		_corr_all = true;
	}
	if(include && !_corr_all){
		if(str == _p_corr_){
			_p_corr = true;
			_e_theta_corr = true;
			_e_p_corr = true;
		}
		if(str == _e_theta_corr_){
			_e_theta_corr = true;
		}
		if(str == _e_p_corr_){
			_e_p_corr = true;
		}
		if(str == _en_corr_){
			_en_corr = true;
		}
		if(str == _v_corr_){
			_v_corr = true;
		}
	}
	if(!include && _corr_all){
		if(str == _p_corr_){
			_p_corr = false;
		}
		if(str == _e_theta_corr_){
			_e_theta_corr = false;
		}
		if(str == _e_p_corr_){
			_e_p_corr = false;
		}
		if(str == _en_corr_){
			_en_corr = false;
		}
		if(str == _v_corr_){
			_v_corr = false;
		}
	}
}
//Run Info
int Flags::Run(){
	return _run_group;
}
bool Flags::Sim(){
	return _sim;
}
bool Flags::Fill(){
	return _filled;
}
bool Flags::Helicity(){
	return _helicity;
}
std::string Flags::Files(){
	return _location;
}
std::string Flags::Base(){
	return _file_base;
}
int Flags::Num_Files(){
	return _num_files;
}
//PID Cuts 
bool Flags::Fid_Cut(int particle){
	return _fid_cut[particle];
}
bool Flags::Delta_Cut(int particle){
	return _delta_cut[particle];
}
bool Flags::SF_Cut(){
	return _sf_cut;
}
bool Flags::CC_Cut(){
	return _cc_cut; 
}
bool Flags::EC_Cut(){
	return _ec_cut; 
}
bool Flags::ID_Cut(){
	return _id_cut; 
}
bool Flags::Beta_Cut(int particle){
	return _beta_cut[particle];
}
bool Flags::Geo_Cut(int particle_, int detector_){
	return _fid_geo_cut[particle_][detector_];
}
bool Flags::Kin_Eff_Cut(int particle_){
	return _kin_eff_cut[particle_];
}

//Other Cuts
bool Flags::Vertex_Cut(){
	return _vertex_cut;
}
bool Flags::WQ2_Cut(){
	return _wq2_cut;
}
//Event Selection 
bool Flags::MM_Cut(int topology){
	return _mm_cut[topology];
}
//Efficiency Cuts
bool Flags::CC_Eff(){
	return _cc_eff_cut;
}
bool Flags::EC_Eff(){
	return _ec_eff_cut;
}
bool Flags::SC_Eff(){
	return _sc_eff_cut;
}
bool Flags::DC_Eff(){
	return _dc_eff_cut;
}
bool Flags::CC_Geo(){
	return _cc_geo_cut;
}
bool Flags::EC_Geo(){
	return _ec_geo_cut;
}
bool Flags::SC_Geo(){
	return _sc_geo_cut;
}
//Corrections
bool Flags::P_Corr(){
	return _p_corr;
}
bool Flags::E_Theta_Corr(){
	return _e_theta_corr;
}
bool Flags::E_PCorr(){
	return _e_p_corr;
}
bool Flags::E_Corr(){
	return _en_corr;
}
bool Flags::Vertex_Corr(){
	return _v_corr;
}
//Histograms
bool Flags::Plot_WQ2(){
	return _plot_wq2;
}
bool Flags::Plot_Fid(int particle){
	return _plot_fid[particle];
}
bool Flags::Plot_SF(){
	return _plot_sf;
}
bool Flags::Plot_CC(){
	return _plot_cc;
}
bool Flags::Plot_CC_Eff(){
	return _plot_cc_eff;
}
bool Flags::Plot_SC_Eff(){
	return _plot_sc_eff;
}
bool Flags::Plot_Delta(int particle){
	return _plot_dt[particle];
}
bool Flags::Plot_Fid_Geo(int particle_, int detector_){
	return _plot_fid_geo[particle_][detector_];
}
bool Flags::Plot_MM(int topology){
	return _plot_mm[topology];
}
bool Flags::Plot_Beta(int particle){
	return _plot_beta[particle];
}
bool Flags::Plot_Vertex(){
	return _plot_vertex;
}
bool Flags::Plot_E_PCorr(){
	return _plot_e_pcorr;
}
bool Flags::Plot_Elastic(){
	return _plot_elastic;
}
bool Flags::Plot_Check(){
	return _plot_check;
}
bool Flags::Plot_Kin_Eff(int particle_){
	return _plot_kin_eff[particle_];
}
bool Flags::Plot_Bin_Centering(){
	return _plot_bin_centering;
}
bool Flags::Plot_Bin_Centering2(){
	return _plot_bin_centering2;
}
bool Flags::Plot_Electron_Angle_Corr(){
	return _plot_ecorr_angle;
}
bool Flags::Plot_Electron_Mag_Corr(){
	return _plot_ecorr_mag;
}

//Histogram Separation

//THnSparse
bool Flags::Make_Friend(){
	return _make_friend;
}

bool Flags::Make_Image(){
	return _make_image;
}

std::string Flags::Output_Name(){
	return _output_name;
}
std::string Flags::Friend_Name(){
	return _friend_name;
}

std::string Flags::Image_Name(){
	return _image_name;
}
int Flags::Num_Cores(){
	if(_num_cores > 0){
		return _num_cores;
	}else{
		return 1;
	}
}

void Flags::Print_Flags(){
	std::cout<<"Printing all Flag Variables\n";
	std::cout<<"Run: " <<Run()<<"\n";
	std::cout<<"Sim: " <<Sim()<<"\n";
	std::cout<<"Fill: " <<Fill()<<"\n";
	std::cout<<"Helicity: " <<Helicity()<<"\n";
	std::cout<<"Files: " <<Files()<<"\n";
	std::cout<<"Base: " <<Base()<<"\n";
	std::cout<<"Num Files: " <<Num_Files()<<"\n";
	//PID Cuts 
	std::cout<<"Fid Cut:\n" ;
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<i <<Fid_Cut(i)<<"\n";
	}
	std::cout<<"Delta Cut:\n"; 
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<i <<Delta_Cut(i)<<"\n";
	}
	std::cout<<"SF Cut: " <<SF_Cut()<<"\n";
	std::cout<<"CC Cut: " <<CC_Cut()<<"\n";
	std::cout<<"EC Cut: " <<EC_Cut()<<"\n";
	std::cout<<"ID Cut: " <<ID_Cut()<<"\n";
	std::cout<<"Beta Cut:\n";
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<i <<Beta_Cut(i)<<"\n";
	} 
	std::cout<<"Geo Cut:\n";
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++){
			std::cout<<"\t" <<i <<" " <<j <<Geo_Cut(i,j) <<"\n";
		}
	}
	std::cout<<"Kinematic Efficiency Cut:\n";
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<i <<Kin_Eff_Cut(i)<<"\n";
	} 
	//Other Cuts
	std::cout<<"Vertex Cut: " <<Vertex_Cut()<<"\n";
	std::cout<<"WQ2 Cut: " <<WQ2_Cut()<<"\n";
	//Event Selection 
	std::cout<<"MM Cut:\n" ;
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<i <<MM_Cut(i)<<"\n";
	}
	//Efficiency Cuts
	std::cout<<"CC Eff: " <<CC_Eff()<<"\n";
	std::cout<<"EC Eff: " <<EC_Eff()<<"\n";
	std::cout<<"SC Eff: " <<SC_Eff()<<"\n";
	std::cout<<"DC Eff: " <<DC_Eff()<<"\n";
	std::cout<<"CC Geo: " <<CC_Geo()<<"\n";
	std::cout<<"EC Geo: " <<EC_Geo()<<"\n";
	std::cout<<"SC Geo: " <<SC_Geo()<<"\n";
	//Corrections
	std::cout<<"P Corr: " <<P_Corr()<<"\n";
	std::cout<<"E Theta Corr: " <<E_Theta_Corr()<<"\n";
	std::cout<<"E PCorr: " <<E_PCorr()<<"\n";
	//std::cout<<"Delta Corr: " <<Delta_Corr()<<"\n";
	std::cout<<"e Corr: " <<E_Corr()<<"\n";
	std::cout<<"Vertex Corr: " <<Vertex_Corr()<<"\n";
	//Histograms ,
	std::cout<<"Plot WQ2: " <<Plot_WQ2()<<"\n";
	std::cout<<"Plot Fid:\n"; 
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<i <<Plot_Fid(i)<<"\n";
	}
	std::cout<<"Plot SF: " <<Plot_SF()<<"\n";
	std::cout<<"Plot CC: " <<Plot_CC()<<"\n";
	std::cout<<"Plot CC Eff: " <<Plot_CC_Eff()<<"\n";
	std::cout<<"Plot SC Eff: " <<Plot_SC_Eff()<<"\n";
	//std::cout<<"Plot CC Eff: " <<Plot_CC_Eff()<<"\n";
	std::cout<<"Plot Delta:\n" ;
	for(int i= 0; i<4; i++){
		std::cout<<"\t" <<i <<Plot_Delta(i) <<"\n";
	}
	std::cout<<"Plot Geo:\n";
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++){
			std::cout<<"\t" <<i <<" " <<j <<" perf:"<<Plot_Fid_Geo(i,j) <<"\n";
		}
	}
	std::cout<<"Plot Kinematic Efficiency:\n";
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<i <<Plot_Kin_Eff(i)<<"\n";
	}
	std::cout<<"Plot MM:\n";
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<i <<Plot_MM(i) <<"\n";
	}
	std::cout<<"Plot Beta:\n" ;
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<i <<Plot_Beta(i) <<"\n";
	}
	std::cout<<"Plot Vertex: " <<Plot_Vertex() <<"\n";
	std::cout<<"Plot E PCorr: " <<Plot_E_PCorr()<<"\n";
	std::cout<<"Plot Elastic: " <<Plot_Elastic()<<"\n";
	std::cout<<"Plot Check: " <<Plot_Check()<<"\n";
	//Portions of Histograms
	//Ele_Cut(const char * ecut_);

	//Histogram Separation

	//THnSparse
	std::cout<<"Make Friend: " <<Make_Friend()<<"\n";
	std::cout<<"Make Image: " <<Make_Image()<<"\n";

	//Output File
	std::cout<<"Output Name: " <<Output_Name()<<"\n";
	std::cout<<"Friend Name: " <<Friend_Name()<<"\n";
	std::cout<<"Image Name: " <<Image_Name()<<"\n";

	std::cout<<"Num Cores: " <<Num_Cores()<<"\n";
}

int Flags::Fid_Cut_Width(int i_){
	if(i_<4){
		if(Flags::Plot_Fid(i_)){
			if(_loose_cut==_fid_cut_f_[i_]){
				return 2;
			}else if(_tight_cut==_fid_cut_f_[i_]){
				return 0;
			}else{
				return 1;
			}
		}else{
			return 1;
		}
	}else{
		return 1;
	}
}

int Flags::Delta_Cut_Width(int i_){
	if(i_<4){
		if(Flags::Plot_Delta(i_)){
			if(_loose_cut==_delta_cut_f_[i_]){
				return 2;
			}else if(_tight_cut==_delta_cut_f_[i_]){
				return 0;
			}else{
				return 1;
			}
		}else{
			return 1;
		}
	}else{
		return 1;
	}
}
int Flags::Fid_Geo_Cut_Width(int i_, int detector_){
	if(i_<4 && detector_ < 3){
		if(_loose_cut==_fid_geo_cut_f_[i_][detector_]){
			return 2;
		}else if(_tight_cut==_fid_geo_cut_f_[i_][detector_]){
			return 0;
		}else{
			return 1;
		}
	}else{
		return 1;
	}
}
int Flags::SF_Cut_Width(){
	std::string cut = _sf_cut_f_;
	if(_loose_cut==cut){
		return 2;
	}else if(_tight_cut==cut){
		return 0;
	}else{
		return 1;
	}
}
int Flags::Min_CC_Cut_Width(){
	std::string cut = _cc_cut_f_;
	if(_loose_cut==cut){
		return 2;
	}else if(_tight_cut==cut){
		return 0;
	}else{
		return 1;
	}
}
int Flags::MM_Cut_Width(int top_){
	std::string cut = "";
	cut = _mm_cut_f_[top_];
	if(_loose_cut==cut || _loose_cut == "mall"){
		return 2;
	}else if(_tight_cut==cut || _tight_cut == "mall"){
		return 0;
	}else{
		return 1;
	}
}
int Flags::Vertex_Cut_Width(){
	std::string cut = _vertex_cut_f_;
	if(_loose_cut==cut){
		return 2;
	}else if(_tight_cut==cut){
		return 0;
	}else{
		return 1;
	}
}
int Flags::Kin_Eff_Cut_Width(int species_){
	std::string cut = _kin_eff_cut_f_[species_];
	if(_loose_cut==cut){
		return 2;
	}else if(_tight_cut==cut){
		return 0;
	}else{
		return 1;
	}
}