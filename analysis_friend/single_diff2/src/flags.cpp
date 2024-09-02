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

	void Flags::Read_Flags(int argc_, char** argv_){
		if(argc_ > 1){
			for(int i=0; i<argc_-1; i++){
				//std::cout<<"Reading Flag: " <<argv[i+1] <<std::endl;
				Flags::ID_Flag(std::string(argv_[i+1]));
			}
		}
	}
	void Flags::ID_Flag(std::string argi_){
		std::string str=argi_.c_str();
		int str_len = str.length();
		std::cout<<"IDing Flags:\n\t"<<"flag:" <<str <<" | len:" <<str_len <<std::endl;
		if(str_len > 3){
			if(str.substr(0,3) == _give_top_){
				_input_top = str.substr(3,str_len);
				std::cout<<"\tIDed: input top is " <<_input_top <<"\n";
			}
			if(str.substr(0,3) == _real_top_){
				_real_top = str.substr(3,str_len);
				std::cout<<"\tIDed: real top is " <<_real_top <<"\n"; 
			}
			if(str.substr(0,3) == _histograms_){
				Flags::Plot_Flag(str.substr(3,str_len));
				std::cout<<"\tIDed: Plots to Make " <<str.substr(3,str_len) <<"\n"; 
			}
			if(str.substr(0,3) == _info_){
				Flags::Info_Flag(str.substr(3,str_len));
				std::cout<<"\tIDed: Info added " <<str.substr(3,str_len) <<"\n";
			}
			if(str.substr(0,3) == _lumen_){
				std::cout<<"\tIDed: Luminosity added " <<str.substr(3,str_len) <<"\n";
				_luminosity=stof(str.substr(3,str_len));
			}
		}if(str_len >4){
			if(str.substr(0,4) == _qratio_){
				std::cout<<"\tIDed: Qratio: " <<str.substr(4,str_len) <<"\n";
				_Qr = stof(str.substr(4,str_len));
			}
			if(str.substr(0,4) == _charge1_){
				std::cout<<"\tIDed: Charge1: " <<str.substr(4,str_len) <<"\n";
				_charge1 = stof(str.substr(4,str_len));
			}
			if(str.substr(0,4) == _charge2_){
				std::cout<<"\tIDed: Charge2: " <<str.substr(4,str_len) <<"\n";
				_charge2 = stof(str.substr(4,str_len));
			}
			if(str.substr(0,4) == _lumen2_){
				std::cout<<"\tIDed: Luminosity added " <<str.substr(4,str_len) <<"\n";
				_luminosity2=stof(str.substr(4,str_len));
			}
		}if(str_len > 5){
			if(str.substr(0,5) == _exp_file_){
				_exp_loc = str.substr(5,str_len);
				_has_exp = true;
				std::cout<<"\tIDed: exp file is " <<_exp_loc <<"\n"; 
			}
			if(str.substr(0,5) == _sim_file_){
				_sim_loc = str.substr(5,str_len);
				_has_sim = true;
				std::cout<<"\tIDed: sim file is " <<_sim_loc <<"\n"; 
			}
			if(str.substr(0,5) == _qratio2_){
				std::cout<<"\tIDed: Qratio: " <<str.substr(5,str_len) <<"\n";
				_Qr2 = stof(str.substr(5,str_len));
			}
			/*if(str.substr(0,5) == _exp_pos_file_){
				_exp_pos_loc = str.substr(5,str_len);
				_helicity = true;
				std::cout<<"\tIDed: exp pos file is " <<_exp_pos_loc <<"\n"; 
			}
			if(str.substr(0,5) == _exp_neg_file_){
				_exp_neg_loc = str.substr(5,str_len);
				_helicity = true;
				std::cout<<"\tIDed: exp neg file is " <<_exp_neg_loc <<"\n"; 
			}*/
		}if(str_len > 6){
			if(str.substr(0,6) == _exp_file2_){
				_exp_loc2 = str.substr(6,str_len);
				_has_exp = true;
				std::cout<<"\tIDed: exp file2 is " <<_exp_loc2 <<"\n"; 
			}
			if(str.substr(0,6) == _sim_file2_){
				_sim_loc2 = str.substr(6,str_len);
				_has_sim = true;
				std::cout<<"\tIDed: sim file2 is " <<_sim_loc2 <<"\n"; 
			}
			if(str.substr(0,6) == _output_name_){
				_output_name = str.substr(6,str_len);
				std::cout<<"\tIDed: output file is " <<_output_name <<"\n"; 
			}
			/*if(str.substr(0,6) == _exp_pos_file2_){
				_exp_pos_loc2 = str.substr(6,str_len);
				_helicity = true;
				std::cout<<"\tIDed: exp pos file2 is " <<_exp_pos_loc2 <<"\n"; 
			}
			if(str.substr(0,6) == _exp_neg_file2_){
				_exp_neg_loc2 = str.substr(6,str_len);
				_helicity = true;
				std::cout<<"\tIDed: exp neg file2 is " <<_exp_neg_loc2 <<"\n"; 
			}*/
		}if(str_len > 7){
			if(str.substr(0,7) == _output_name1_){
				_output_name1 = str.substr(7,str_len);
				std::cout<<"\tIDed: output file is " <<_output_name1 <<"\n"; 
			}
			if(str.substr(0,7) == _output_name2_){
				_output_name2 = str.substr(7,str_len);
				std::cout<<"\tIDed: output file is " <<_output_name2 <<"\n"; 
			}
			if(str.substr(0,7) == _output_name3_){
				_output_name3 = str.substr(7,str_len);
				std::cout<<"\tIDed: output file is " <<_output_name3 <<"\n"; 
			}
			if(str.substr(0,7) == _empty_file_){
				_empty_loc = str.substr(7,str_len);
				_has_empty = true;
				std::cout<<"\tIDed: empty file is " <<_empty_loc <<"\n"; 
			}
			if(str.substr(0,7) == _image_name_){
				_image_name = str.substr(7,str_len);
				_make_image = true;
				std::cout<<"\tIDed: Image file is " <<_image_name <<"\n"; 
			}
			if(str.substr(0,7) == _loc_holes1_){
				_localized_holes_name = str.substr(7,str_len);
				_localized_holes = true;
				_has_localized_holes=true;
				std::cout<<"\tIDed: Local Holes file 1 is " <<_localized_holes_name <<"\n"; 
			}
			
		}if(str_len > 8){
			/*if(str.substr(0,8) == _weight_file_){
				_weight_loc = str.substr(8,str_len);
				_has_weight = true;
				std::cout<<"\tIDed: weight file is " <<_weight_loc <<"\n"; 
			}*/
			if(str.substr(0,8) == _empty_file_){
				_empty_loc2 = str.substr(8,str_len);
				_has_empty = true;
				std::cout<<"\tIDed: empty file2 is " <<_empty_loc <<"\n"; 
			}
			if(str.substr(0,8) == _sim_no_rad_file_){
				_sim_no_rad_loc = str.substr(8,str_len);
				_rad_corr = true;
				std::cout<<"\tIDed: Sim No Rad File is " <<_sim_no_rad_loc <<"\n"; 
			}
			if(str.substr(0,8) == _loc_holes2_){
				_localized_holes_name2 = str.substr(8,str_len);
				_localized_holes = true;
				_has_localized_holes=true;
				std::cout<<"\tIDed: Local Holes file 2 is " <<_localized_holes_name2 <<"\n"; 
			}
			
		}
		if(str_len > 9){
			/*if(str.substr(0,9) == _weight_file_){
				_weight_loc2 = str.substr(9,str_len);
				_has_weight = true;
				std::cout<<"\tIDed: weight file2 is " <<_weight_loc2 <<"\n"; 
			}*/
			if(str.substr(0,9) == _sim_no_rad_file2_){
				_sim_no_rad_loc2 = str.substr(9,str_len);
				_rad_corr = true;
				std::cout<<"\tIDed: Sim No Rad File2 is " <<_sim_no_rad_loc2 <<"\n"; 
			}
		}
        if(str_len >= 13){
            if(str.substr(0,13) == _no_loc_holes_){
				_has_localized_holes = false;
				std::cout<<"\tIDed: Sim No Rad File2 is " <<_sim_no_rad_loc2 <<"\n"; 
			}
			if(str.substr(0,13) == _acc_rel_cut_){
				if(str.substr(13,str_len) == _tight_){
					_acc_rel_cut = 0;
				}else if(str.substr(13,str_len) == _mid_){
					_acc_rel_cut = 1;
				}else if(str.substr(13,str_len) == _loose_){
					_acc_rel_cut = 2;
				}else{
					_acc_rel_cut = 3;
				}
				std::cout<<"\tIDed: Acceptance Relative Error Cut idx is " <<_acc_rel_cut <<"\n"; 
			}
        }
		if(str_len >= 14){
			if(str.substr(0,14) == _nonlocal_holes_){
				//_output_name1 = str.substr(7,str_len);
				_nonlocal_holes = true;
				std::cout<<"\tIDed: output file is " <<"nonlocal holes" <<"\n"; 
			}
		}
	}

	void Flags::Plot_Flag(std::string str_){
		if(str_ == _all_){
			_plot_all = true;
			_helicity=true;
		}else if(str_ == _polarization_){
			_plot_polarization = true;
		}else if(str_ == _single_diff_){
			_plot_single_diff = true;
		}else if(str_ == _beam_spin_){
			_plot_beam_spin = true;
			_helicity=true;
		}else if(str_ == _eff_ ){
			_plot_eff = true;
		}else if(str_ == _err_ ){
			_plot_err = true;
		}else if(str_ == _wq2_ ){
			_plot_wq2 = true;
		}else if(str_ == _acceptance_ ){
			_plot_acceptance = true;
		}else if(str_ == _pol_mm1_){
			_plot_pol_mm1 = true;
		}else if(str_ == _pol_mm2_){
			_plot_pol_mm2 = true;
		}else if(str_ == _pol_theta_){
			_plot_pol_theta = true;
		}else if(str_ == _pol_alpha_){
			_plot_pol_alpha = true;
		}else if(str_ == _error_){
			_plot_error = true;
		}else if(str_ == _local_holes_){
			_localized_holes = true;
		}else if(str_ == _plot_local_holes_){
			_plot_localized_holes = true;
			//_helicity=true;//in order to make helicity dependent localized holes!
		}
	}

	void Flags::Info_Flag(std::string str_){
		//std::cout<<"Current Run Info:\n\trun group: " <<_run_group <<"\n\tvar set: " <<_var_set <<"\n\tvar idx: "<<_var_idx <<"\n";
		if(str_ == _run_e16_){
			_run_group = 0;
			std::cout<<"Run group e16\n";
		}else if(str_ == _run_e1f_){
			_run_group = 1;
			std::cout<<"Run group e1f\n";
		}else if(str_ == _both_){
			_run_group = 2;
			std::cout<<"Run group both\n";
		}else if(str_ == _var_pim_){
			_var_set = "pim";
			_var_idx = 0;
		}
		else if(str_ == _var_pro_){
			_var_set = "pro";
			_var_idx = 1;
		}else if(str_ == _var_pip_){
			_var_set = "pip";
			_var_idx = 2;
		}else if(str_ == _flux_included_){
			_flux_inc = true;
		}else if(str_ == _helicity_ ){
			_helicity=true;
		}
		//std::cout<<"New Run Info:\n\trun group: " <<_run_group <<"\n\tvar set: " <<_var_set <<"\n\tvar idx: "<<_var_idx <<"\n";
	}

	int Flags::Run(){
		return _run_group;
	}

	std::string Flags::Run(std::string str_){
		if(_run_group == 0){
			return "e16";
		}else if(_run_group == 1){
			return "e1f";
		}else if(_run_group == 2){
			return "both";
		}else{
			return "improper";
		}
	}

	//Plotting
	bool Flags::Plot_Pol(){
		return (_plot_all || _plot_polarization);
	}
	bool Flags::Plot_Single_Diff(){
		return (_plot_all || _plot_single_diff);
	}
	bool Flags::Plot_Beam_Spin(){
		return (_plot_all || _plot_beam_spin);
	}
	bool Flags::Plot_Eff(){
		return (_plot_all || _plot_eff);
	}
	bool Flags::Plot_Err(){
		return (_plot_all || _plot_err);
	}
	bool Flags::Plot_WQ2(){
		return (_plot_all || _plot_wq2);
	}
	bool Flags::Plot_Acceptance(){
		return (_plot_all || _plot_acceptance);
	}

	//File Names
	std::string Flags::Sim_File(){
		return _sim_loc;
	}
	std::string Flags::Exp_File(){
		return _exp_loc;
	}
	std::string Flags::Empty_File(){
		return _empty_loc;
	}
	std::string Flags::Weight_File(){
		return _weight_loc;
	}
	std::string Flags::Sim_File2(){
		return _sim_loc2;
	}
	std::string Flags::Exp_File2(){
		return _exp_loc2;
	}
	std::string Flags::Empty_File2(){
		return _empty_loc2;
	}
	std::string Flags::Weight_File2(){
		return _weight_loc2;
	}
	std::string Flags::Output_File(){
		return _output_name;
		/*switch(i_){
			case 0:
				return _output_name1;
			break;
			case 1:
				return _output_name2;
			break;
			case 2:
				return _output_name3;
			break;
			default:
				std::cout<<"You should specifcy which rootfile\n";
				return _output_name;
			break;
		}*/
	}
	std::string Flags::Image_File(){
		return _image_name;
	}
	std::string Flags::Exp_Hel(){
		return _exp_loc;//The two helicity histograms are in the same root file now 1/9/23
	}
	std::string Flags::Exp_Hel2(){
		return _exp_loc2;
	}
	std::string Flags::Sim_No_Rad(){
		return _sim_no_rad_loc;
	}
	std::string Flags::Sim_No_Rad2(){
		return _sim_no_rad_loc2;
	}
	std::string Flags::Holes_File(){
		return _localized_holes_name;
	}
	std::string Flags::Holes_File2(){
		return _localized_holes_name2;
	}
	//Have Things
	bool Flags::Has_Weight(){
		return _has_weight;
	}
	bool Flags::Has_Sim(){
		return _has_sim;
	}
	bool Flags::Has_Exp(){
		return _has_exp;
	}
	bool Flags::Has_Empty(){
		return _has_empty;
	}
	bool Flags::Make_Image(){
		return _make_image;
	}

	int Flags::Var_idx(){
		return _var_idx;
	}
	std::string Flags::Var_Set(){
		return _var_set;
	}

	std::string Flags::Top(){
		return _real_top;
	}
	int Flags::Top_idx(){
		for(int i=0; i<5; i++){
			if(_top_[i] == _real_top){
				return i;
			}
		}
	}
	bool Flags::Helicity(){
		return _helicity;
	}
	bool Flags::Rad_Corr(){
		return _rad_corr;
	}
	bool Flags::Flux_Included(){
		return _flux_inc;
	}

	float Flags::Qr(){
		return _Qr;
	}
	float Flags::Qr2(){
		return _Qr2;
	}

	double Flags::L(int i_){
        if(i_==0){
			return _luminosity; 
		}else if(i_==1){
			return _luminosity2;
		}
		/*double l_base = 1.36252*pow(10.0,12);// (1/microbarns)*(1/Coulombs) 
		if(i_==0){
			return _charge1*l_base; 
		}else if(i_==1){
			return _charge2*l_base;
		}*/
	}

	bool Flags::Plot_Polarization(int i_){
		switch(i_){
			case 0:
				return _plot_pol_mm1;
			break;
			case 1:
				return _plot_pol_mm2;
			break;
			case 2:
				return _plot_pol_theta;
			break;
			case 3:
				return _plot_pol_alpha;
			break;
			default:
				std::cout<<"Not a valid polarization part\n";
				return false;
			break;
		}
	}

	bool Flags::Error(){
		return _plot_error;
	}

	bool Flags::Localized_Holes(){
		return _localized_holes;
	}

	bool Flags::Plot_Localized_Holes(){
		return _plot_localized_holes;
	}

	bool Flags::Has_Localized_Holes(){
		return _has_localized_holes;
	}

	bool Flags::Nonlocal_Holes(){
		return _nonlocal_holes;
	}

	int Flags::Acc_Rel_Error_Cut(){
		return _acc_rel_cut;
	}