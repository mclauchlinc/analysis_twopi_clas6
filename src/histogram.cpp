#include "histogram.hpp"


Histogram::Histogram(std::shared_ptr<Flags> flags_){
	/*std::cout<<"Test!\nP min for bins:\n";
	std::cout<<"\tmin: " <<_p_bin_min_ <<"\n";
	for(int i=0; i<_p_bin_bins_; i++){
		std::cout<<"\t" <<Histogram::P_Min(i) << " center-> " <<Histogram::P_center(i) <<"\n";
	}*/
	std::cout<<"Num W bins: " <<Histogram::W_bins() <<"\n";
	//_RootOutputFile = fun::Name_File(output_file);
	//def = new TCanvas("def");
	Histogram::WQ2_Make(flags_);
	Histogram::Fid_Make(flags_);
	Histogram::SF_Make(flags_);
	Histogram::Delta_Make(flags_);
	Histogram::CC_Make(flags_);
	Histogram::Vertex_Make(flags_);
	Histogram::MM_Make(flags_);
	Histogram::Kinematic_Eff_Make(flags_);
	Histogram::Kin_Eff_SC_Seg_Make(flags_);
	Histogram::Friend_Make(flags_);
	Histogram::PCorr_Check_Make(flags_);
	Histogram::SC_Eff_Make(flags_);
	Histogram::Elastic_Peak_Make(flags_);
	//Histogram::Cross_Make(flags_);
	//Histogram::XY_Make(flags_);
	//Histogram::Fid_Det_Make(flags_);
	//Histogram::Charge_Make(flags_);
	//Histogram::Thrown_Make(flags_);
	//Histogram::WQ2_sf_Make(_envi);
	Histogram::Bin_Centering_Make(flags_);
	std::cout<<"made it past here1\n";
	Histogram::Bin_Centering2_Make(flags_);
	std::cout<<"made it past here2\n";
	Histogram::Ele_Angle_Corr_Make(flags_);
	Histogram::Ele_Mag_Corr_Make(flags_);
	Histogram::Make_Top_Check(flags_);
	Histogram::Geo_Fid_Make(flags_);
	Histogram::Sim_Vertex_Make(flags_);
	Histogram::Hel_Ratio_Make(flags_);
	Histogram::Golden_Make(flags_);
}

bool Histogram::OK_Idx(std::vector<int> idx_){
	bool pass = true;//Assume it's ok. Maybe dangerous, but I'll take that gamble for now
	for(int i=0; i<idx_.size(); i++){
		if(idx_[i]<0){
			pass = false;
		}
	}
	return pass;
}

//Histogram::~Histogram() { this->Write(); }

void Histogram::Write(std::shared_ptr<Flags> flags_){
	std::cout<<"Writing Histogram\n\tFile Named: "<<flags_->Flags::Output_Name() <<"\n";
	std::cout<<"Making File:";
	_RootOutputFile = fun::Name_File(flags_);
	std::cout<<"Success\n";
	def = new TCanvas("def");
	std::cout<< "Writing Plots" <<std::endl;
	_RootOutputFile->cd();
	std::cout<<"Thinking about Writing WQ2 Histograms\n";
	Histogram::WQ2_Write(flags_);
	std::cout<<"Thinking about Writing Fid Histograms\n";
	Histogram::Fid_Write(flags_);
	std::cout<<"Thinking about Writing SF Histograms\n";
	Histogram::SF_Write(flags_);
	std::cout<<"Thinking about Writing Delta Histograms\n";
	Histogram::Delta_Write(flags_);
	std::cout<<"Thinking about Writing CC Histograms\n";
	Histogram::CC_Write(flags_);
	std::cout<<"Thinking about Writing Vertex Histograms\n";
	Histogram::Vertex_Write(flags_);
	std::cout<<"Thinking about Writing MM Histograms\n";
	Histogram::MM_Write(flags_);
	std::cout<<"Thinking about Writing Kinematic Eff Histograms\n";
	Histogram::Kinematic_Eff_Write(flags_);
	std::cout<<"Thinking about writing Other Kinematic Eff Histograms\n";
	Histogram::Kin_Eff_SC_Seg_Write(flags_);
	std::cout<<"Thinking about Writing PCorr Histograms\n";
	Histogram::PCorr_Check_Write(flags_);
	std::cout<<"Thinking about Writing SC Eff Histograms\n";
	Histogram::SC_Eff_Write(flags_);
	std::cout<<"Thinking about Writing Bin Centering Histograms\n";
	Histogram::Bin_Centering_Write(flags_);
	std::cout<<"Thinking about Writing Bin Centering2 Histograms\n";
	Histogram::Bin_Centering2_Write(flags_);
	std::cout<<"Thinking about Writing Ele angle Corr Histograms\n";
	Histogram::Ele_Angle_Corr_Write(flags_);
	std::cout<<"Thinking about Writing Ele Mag Corr Histograms\n";
	Histogram::Ele_Mag_Corr_Write(flags_);
	std::cout<<"Thinking about Writing Elastic Peak Histograms\n";
	Histogram::Elastic_Peak_Write(flags_);
	std::cout<<"Thinking about Writing Geo Fid Histograms\n";
	Histogram::Geo_Fid_Write(flags_);
	std::cout<<"Thinking about writing Top Check histograms\n";
	Histogram::Write_Top_Check(flags_);
	std::cout<<"Thinking about writing Sim Vertex Runs\n";
	Histogram::Sim_Vertex_Write(flags_);
	std::cout<<"Thinking about writing Helicity Ratio Histograms\n";
	Histogram::Hel_Ratio_Write(flags_);
	std::cout<<"Thinking about writing Golden Histograms\n";
	Histogram::Golden_Write(flags_);
	//Histogram::XY_Write(_envi);
	//Histogram::Fid_Det_Write(_envi);
	//Friend_Write(_envi);
	//Histogram::Cross_Write(_envi);
	//Histogram::Charge_Write(_envi);
	//Histogram::WQ2_sf_Write(_envi);
	_RootOutputFile->Close();
	if(flags_->Flags::Make_Friend()){
		std::cout<<"Writing THnSparse Histograms: " <<flags_->Flags::Friend_Name()<<"\n";
		_SparseFile = fun::Name_Sparse(flags_);
		Histogram::Friend_Write(flags_);
		_SparseFile->Close();
		std::cout<<" Done\n";
	}
	std::cout<<"Histograms Done!" <<std::endl;
}

/*
void Histogram::Print(const std::string& output_dir, std::shared_ptr<Environment> envi_){
	int _check_ = -1;
	
	//if(envi_->was_print()>0){
		std::string _curr_dir_ = fun::get_current_dir();
		std::cout<<"Current Directory: " <<_curr_dir_ <<std::endl;
		std::string _out_dir_ = "$curr/$name/plots";
		/*if(envi_->was_cluster()){
			_out_dir_ = "$curr/$name/Plots";
		}
		else{
			_out_dir_ = "$curr/$name/Plots";
		}*/ /*
		std::cout<<"Printing Plots \n";
		//HistImageFile = fun::Name_Image_File(output_dir);
		std::cout<<"Making Plot Directory:";
		fun::replace(_out_dir_, "$curr", _curr_dir_);
		fun::replace(_out_dir_, "$name", output_dir);
		_check_ = fun::Make_Dir(_out_dir_);
		std::cout<<_out_dir_ <<std::endl;
		std::cout<<" Done" <<std::endl;
		//chdir(_out_dir_.c_str());
		//std::cout<<" Done \n 	Printing Histograms";
	//}
	//if(_check_ == 0){
		//std::cout<<"Pass the check?";
		//Histogram::WQ2_Print(output_dir, envi_);
		//Histogram::Fid_Print(const std::string& output_dir, envi_);
		//Histogram::SF_Print(const std::string& output_dir, envi_);
		//Histogram::DT_Print(const std::string& output_dir, envi_);
		//Histogram::CC_Print(const std::string& output_dir, envi_);
		//Histogram::MM_Print(const std::string& output_dir, envi_);
		//Histogram::XY_Print(const std::string& output_dir,envi_);
		//Histogram::Fid_Det_Print(const std::string& output_dir,envi_);
		//Friend_Print(const std::string& output_dir,envi_);
		//Histogram::Cross_Print(const std::string& output_dir,envi_);
	//}	
}*/

int Histogram::P_bin(float p_){
	if(p_<=_p_bin_min_ && p_ > _p_bin_max_){
		return -1;
	}else{
		float p_res = ((float)_p_bin_max_-(float)_p_bin_min_)/((float)_p_bin_bins_);
		for(int i=0; i<_p_bin_bins_; i++){
			if(p_>(_p_bin_min_+((float)i*p_res)) && p_<=(_p_bin_min_+(((float)i+1.0)*p_res))){
				return i;
			}
		}
	}
}

float Histogram::P_Min(int p_bin_){
	if(p_bin_==-1){
		return NAN;
	}else{
		float p_res = ((float)_p_bin_max_-(float)_p_bin_min_)/((float)_p_bin_bins_);
		return _p_bin_min_ + (float)p_bin_*p_res;
	}
}

float Histogram::P_Max(int p_bin_){
	if(p_bin_==-1){
		return NAN;
	}else{
		float p_res = ((float)_p_bin_max_-(float)_p_bin_min_)/((float)_p_bin_bins_);
		return _p_bin_min_ + ((float)p_bin_+1.0)*p_res;
	}
}

float Histogram::P_center(int p_bin_){
	if(p_bin_>=0){
		float p_res = ((float)_p_bin_max_-(float)_p_bin_min_)/((float)_p_bin_bins_);
		return _p_bin_min_ + ((float)p_bin_+0.5)*p_res;
	}else{
		return NAN; 
	}
}

int Histogram::W_bins(){
	float W_top = _W_min_;
	int W_bins = 0; 
	while(W_top < _W_max_){
		W_top += _W_res_;
		if(W_top < _W_max_){
			W_bins++;
		}
	}
	return W_bins;
}

int Histogram::W_bin(float W_){
	float W_bot = _W_min_;
	int W_bin = -1; 
	//std::cout<<"W_bot: " <<W_bot <<" W: " <<W_ <<"  W_top: "
	for(int i=0; i<Histogram::W_bins(); i++){
		if(W_ >= (W_bot+i*_W_res_) && W_ < W_bot+((1+i)*_W_res_)){
			W_bin = i;
		}
	}
	return W_bin;
}

float Histogram::W_bot(int bin_){
	return _W_min_ + bin_*_W_res_;
}

float Histogram::W_top(int bin_){
	return _W_min_ + (bin_+1)*_W_res_;
}

float Histogram::W_center(int bin_){
	return  _W_min_ + ((float)bin_+0.5)*_W_res_;
}

int Histogram::Q2_bin(float Q2_){
	int Q2_bin = -1;
	for(int i=0; i<5; i++){
		if(Q2_>= _Q2_bins_[i] && Q2_ < _Q2_bins_[i+1]){
			Q2_bin = i;
		}
	}
	return Q2_bin;
}

float Histogram::Q2_bot(int bin_){
	return _Q2_bins_[bin_];
}

float Histogram::Q2_top(int bin_){
	return _Q2_bins_[bin_+1];
}	


//W Qsquared plots
int Histogram::W_binning(float W_){
  int j = 0;
  float top, bot; 
  for(int i = 1; i < (xc_bins[1]+1); i++){
    top = Wbin_start + i*Wbin_res;//constants.hpp
    bot = top - Wbin_res; 
    if(W_ < top && W_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::p_binning(float p_){
  int j = 0;
  float top, bot; 
  for(int i = 1; i < 26; i++){//Changed from 12 to 26 7/7/20
    top = pbin_start + i*pbin_res;//constants.hpp
    bot = top - pbin_res; 
    if(p_ < top && p_ >= bot){
      j = i; 
    }
  }
  return j; 
}

char Histogram::Part_cut(int species, int cut){
	char the_cut[100]; 
	switch(species){
		case 0:
		sprintf(the_cut,"%s",eid_cut[cut]);
		break;
		case 1:
		if(cut >= 6){
			sprintf(the_cut,"INVALD_CUT");
		}else{
			sprintf(the_cut,"%s",hid_cut[cut]);
		}
		break;
	}
	//return the_cut; 
}
int Histogram::ecut_idx(char * ecut_, std::shared_ptr<Flags> flags_){
	int idx = -1; 
	int idx_tmp = -1; 
	for(int i = 0; i<std::distance(std::begin(_ecuts_), std::end(_ecuts_)); i++){
		if(ecut_ == _ecuts_[i]){
			idx_tmp += i+1;
		}
	}
	switch(idx_tmp){
		case 0:
			idx = idx_tmp;
		break;
		case 1:
			if(flags_->Flags::Fid_Cut(0)){
				idx = idx_tmp; 
			}
		break;
		case 2:
			if(flags_->Flags::SF_Cut()){
				if(flags_->Flags::Fid_Cut(0)){
					idx = idx_tmp; 
				}else{
					idx = idx_tmp-1; 
				}
			}
		break;
		case 3:
			if(flags_->Flags::CC_Cut()){
				idx = idx_tmp ;
				if(!flags_->Flags::SF_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::Fid_Cut(0)){
					idx += -1; 
				}
			}
		break;
		case 4:
			if(flags_->Flags::EC_Cut()){
				idx = idx_tmp ;
				if(!flags_->Flags::CC_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::SF_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::Fid_Cut(0)){
					idx += -1; 
				}
			}
		break;
		case 5:
			if(flags_->Flags::Vertex_Cut()){
				idx = idx_tmp ;
				if(!flags_->Flags::EC_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::CC_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::SF_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::Fid_Cut(0)){
					idx += -1; 
				}
			}
		break;
		case 6:
			if(flags_->Flags::Vertex_Cut()){
				idx = idx_tmp; 
				if(!flags_->Flags::EC_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::CC_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::SF_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::Fid_Cut(0)){
					idx += -1; 
				}
			}
		break;
		case 7:
			if(flags_->Flags::Vertex_Cut()){
				idx = idx_tmp; 
				if(!flags_->Flags::EC_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::CC_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::SF_Cut()){
					idx += -1; 
				}
				if(!flags_->Flags::Fid_Cut(0)){
					idx += -1; 
				}
			}
		break;
		default:
			std::cout<<"\tInvalid EID index\n";
		break;
	}
	return idx; 
}

int Histogram::ecut_idx(const char * ecut_, std::shared_ptr<Flags> flags_){
	int idx = fun::ecut_idx(ecut_); 
	int idx_tmp = fun::ecut_idx(ecut_);
	//std::cout<<"ecut_idx try: " <<ecut_ <<" "<<idx_tmp <<"\n";
	for(int i=0; i<idx_tmp; i++){
		if(_ecuts_[i] == _sanity_){
		}else if(_ecuts_[i] == _fid_cut_){
			if(!flags_->Flags::Fid_Cut(0)){
				idx += -1; 
			}
		}else if(_ecuts_[i] == _sf_cut_){
			if(!flags_->Flags::SF_Cut()){
				idx += -1; 
			}
		}else if(_ecuts_[i] == _cc_cut_){
			if(!flags_->Flags::CC_Cut()){
				idx += -1; 
			}
		}else if(_ecuts_[i] == _ec_cut_){
			if(!flags_->Flags::EC_Cut()){
				idx += -1; 
			}
		}else if(_ecuts_[i] == _vertex_cut_){
			if(!flags_->Flags::Vertex_Cut()){
				idx += -1; 
			}
		}
	}
	//std::cout<<"\tarrived at: "<<idx <<"\n";
	return idx; 
}
//*------------------------------- Start W Q2 ---------------------------------*
void Histogram::WQ2_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Plot_WQ2()){
		std::cout<<"Making WQ2 Histograms\n";
		TH2F_ptr_4d plot_4d;
		TH2F_ptr_3d plot_3d;
		TH2F_ptr_2d plot_2d;
		TH2F_ptr_1d plot_1d;
		Bool_4d made_4d;
		Bool_3d made_3d;
		Bool_2d made_2d;
		Bool_1d made_1d;
		std::vector<long> space_dims(5);
		space_dims[4] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));//Electron Cuts
		space_dims[3] = std::distance(std::begin(_cut_), std::end(_cut_));//Cut and anti cut
		space_dims[2] = std::distance(std::begin(_top_), std::end(_top_));//{Topologies  + all combined + None}
		space_dims[1] = std::distance(std::begin(_weight_), std::end(_weight_));//Not Weight vs. Weighted
		space_dims[0] = std::distance(std::begin(_recon_), std::end(_recon_));//Reconstructed vs. Thrown
		CartesianGenerator cart(space_dims);
		CartesianGenerator cart2(space_dims);
		char hname[100];//For naming histograms
		Double_t Q2_bin[6] = {_Q2_bins_[0],_Q2_bins_[1],_Q2_bins_[2],_Q2_bins_[3],_Q2_bins_[4],_Q2_bins_[5]};

		for(int i=0; i<5; i++){
			std::cout<<_top_[i] <<" perform?: " <<fun::top_perform(_top_[i], flags_) <<"\n";
			if(fun::top_perform(_top_[i], flags_)){
				sprintf(hname,"W Q2 Yield Top:%s",_top_[i]);
				_WQ2_yield_hist.push_back(new TH2F(hname,hname,Histogram::W_bins(),_W_min_,_W_max_,5,Q2_bin));
			}
		}
		if(flags_->Flags::Sim()){
			sprintf(hname,"W Q2 Yield Thrown");
			_WQ2_yield_hist.push_back(new TH2F(hname,hname,Histogram::W_bins(),_W_min_,_W_max_,5,Q2_bin));
		}
		std::cout<<"\nYield W-Q^2 hists made";
		
		//std::cout<<"Making WQ2 made histogram\n";
		while(cart2.GetNextCombination()){
			/*std::cout<<"\t\tIndex: ";
			for(int i=0; i<5; i++){
				std::cout<<cart2[i];
			}
			std::cout<<"\n";*/
			//if(_cut_[cart2[3]] != _no_cut_){
				made_1d.push_back(false);
			//}
			if(cart2[0]==space_dims[0]-1 && made_1d.size()>0){
				made_2d.push_back(made_1d);
				//std::cout<<"\tLength of made_1d: " <<made_1d.size() <<"\n";
				made_1d.clear();
			}
			if(cart2[1]==space_dims[1]-1 && made_2d.size()>0 && cart2[0]==space_dims[0]-1){
				made_3d.push_back(made_2d);
				//std::cout<<"\tLength of made_2d: " <<made_2d.size() <<"\n";
				made_2d.clear();
			}
			if(cart2[2]==space_dims[2]-1 && made_3d.size()>0 && cart2[1]==space_dims[1]-1 && cart2[0]==space_dims[0]-1){
				made_4d.push_back(made_3d);
				//std::cout<<"\tLength of made_3d: " <<made_3d.size() <<"\n";
				made_3d.clear();
			}
			if(cart2[3]==space_dims[3]-1 && made_4d.size()>0 && cart2[1]==space_dims[1]-1 && cart2[0]==space_dims[0]-1 && cart2[2]==space_dims[2]-1){
				_WQ2_made.push_back(made_4d);
				//std::cout<<"\tLength of made_4d: " <<made_4d.size() <<"\n";
				made_4d.clear();
			}
			if(cart2[4]==space_dims[4]-1 ){
				//std::cout<<"\tLength of made: " <<_WQ2_made.size() <<"\n";
			}
		}
		while(cart.GetNextCombination()){
			/*std::cout<<"\t\tIndex: ";
			for(int i=0; i<5; i++){
				std::cout<<cart[i];
			}
			std::cout<<"\n";*/
			//std::cout<<"Index: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<"\n";
			if((flags_->Flags::Sim() || _weight_[cart[1]]==_nweighted_)){
				if(_recon_[cart[0]]==_thrown_){//Thrown Stuff
					if(flags_->Flags::Sim()){
						if(_ecuts_[cart[4]] == _event_ && _cut_[cart[3]] == _cut_applied_ && _top_[cart[2]] == _mzero_){
							std::cout<<"made WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							
							sprintf(hname,"W_Q2_%s_%s",_recon_[cart[0]],_weight_[cart[1]]);
							std::cout<<"hist name: " <<hname <<"\n";
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
						}
					}
				}else{//Reconstructed Stuff
					if(_ecuts_[cart[4]]== _event_ && _top_[cart[2]]!=_mnone_ && _cut_[cart[3]]!=_no_cut_){//Event Selection
						std::cout<<"\t\tMaking Event Selection WQ2\n";
						std::cout<<"Filled WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
						sprintf(hname,"W_Q2_%s_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]]);
						plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
						_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
						//fun::print_vector_idx(cart);
						fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
					}else if(_ecuts_[cart[4]] != _event_ && _top_[cart[2]]==_mnone_){//Electron ID
						//std::cout<<"\t\tMaking PID WQ2\n";
						if(_ecuts_[cart[4]]== _none_ && _cut_[cart[3]]==_no_cut_){//&& _cut_[cart[3]]==_no_cut_
							//std::cout<<"Filled WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
						if(_ecuts_[cart[4]]== _sanity_ && _cut_[cart[3]]!=_no_cut_){
							//std::cout<<"\t Index: ";
							//std::cout<<"Filled WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
						if(_ecuts_[cart[4]]== _fid_cut_ && flags_->Flags::Fid_Cut(0) && _cut_[cart[3]]!=_no_cut_){
							//std::cout<<"\t Index: ";
							//std::cout<<"Filled WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
						if(_ecuts_[cart[4]]== _sf_cut_ && flags_->Flags::SF_Cut() && _cut_[cart[3]]!=_no_cut_){
							//std::cout<<"\t Index: ";
							//std::cout<<"Filled WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
						if(_ecuts_[cart[4]]== _cc_cut_ && flags_->Flags::CC_Cut() && _cut_[cart[3]]!=_no_cut_){
							//std::cout<<"\t Index: ";
							//std::cout<<"Filled WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
						if(_ecuts_[cart[4]]== _ec_cut_ && flags_->Flags::EC_Cut() && _cut_[cart[3]]!=_no_cut_){
							//std::cout<<"\t Index: ";
							//std::cout<<"Filled WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
						if(_ecuts_[cart[4]]== _vertex_cut_ && flags_->Flags::Vertex_Cut() && _cut_[cart[3]]!=_no_cut_){
							//std::cout<<"\t Index: ";
							//std::cout<<"Filled WQ2 idx: " <<_WQ2_hist.size() <<" " <<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
						if(_ecuts_[cart[4]]== _geo_cc_cut_ && flags_->Flags::Geo_Cut(0,0) && _cut_[cart[3]]!=_no_cut_){
							//std::cout<<"\t Index: ";
							//std::cout<<"Filled WQ2 idx: " <<_WQ2_hist.size() <<" " <<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
						if(_ecuts_[cart[4]]== _geo_sc_cut_ && flags_->Flags::Geo_Cut(0,1) && _cut_[cart[3]]!=_no_cut_){
							//std::cout<<"\t Index: ";
							//std::cout<<"Filled WQ2 idx: " <<_WQ2_hist.size() <<" " <<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
						if(_ecuts_[cart[4]]== _geo_ec_cut_ && flags_->Flags::Geo_Cut(0,2) && _cut_[cart[3]]!=_no_cut_){
							//std::cout<<"\t Index: ";
							//std::cout<<"Filled WQ2 idx: " <<_WQ2_hist.size() <<" " <<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
						if(_ecuts_[cart[4]]== _pid_ && _cut_[cart[3]]!=_no_cut_){
							//std::cout<<"\t Index: ";
							//std::cout<<"Filled WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
							//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						}
					}
				}
			}
			if(cart[0] == space_dims[0]-1){
				//std::cout<<"\t\tPush back 1d\n";
				if(plot_1d.size()>0){
					plot_2d.push_back(plot_1d);
					//std::cout<<"Size of 1d: " <<plot_1d.size() <<" | Size of 2d: " <<plot_2d.size()<<"\n";
					plot_1d.clear();
				}
				if(cart[1] == space_dims[1]-1){
					if(plot_2d.size()>0){
						plot_3d.push_back(plot_2d);
						plot_2d.clear();
					}
					if(cart[2] == space_dims[2]-1){
						if(plot_3d.size()>0){
							plot_4d.push_back(plot_3d);
							plot_3d.clear();
						}
						if(cart[3] == space_dims[3]-1){
							if(plot_4d.size()>0){
								_WQ2_hist.push_back(plot_4d);
								plot_4d.clear();
							}
						}
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::WQ2_cut_idx(const char* ecut_, const char * cut_, const char* top_, const char * weight_, const char * recon_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	int ecut_idx = -1; 
	int cut_idx = -1; 
	//std::cout<<"input to WQ2 Cut idx: " <<ecut_ <<"\t" <<cut_ <<"\t" <<top_ <<"\t" <<weight_ <<"\t" <<recon_  <<"\n";
	//if(cut_!=_no_cut_){
		if(recon_==_thrown_ && flags_->Flags::Sim()){
			if(ecut_ == _event_ && cut_ == _cut_applied_ && top_==_mzero_){//Thrown 
				//std::cout<<"input to WQ2 Thrown Cut idx: " <<ecut_ <<"\n\t" <<cut_ <<"\n\t" <<top_ <<"\n\t" <<weight_ <<"\n\t" <<recon_  <<"\n";
				idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
				idx.push_back(0);
				idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
				idx.push_back(fun::weight_idx(weight_));
				idx.push_back(fun::recon_idx(recon_));
			}else{
				idx.push_back(-1);
				idx.push_back(-1);
				idx.push_back(-1);
				idx.push_back(-1);
				idx.push_back(-1);
			}
		}else if(recon_==_nthrown_){//Reconstructed
			if(ecut_==_event_ && top_!=_mnone_ && cut_!=_no_cut_){
				idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
				idx.push_back(fun::cut_idx(cut_));
				idx.push_back(fun::top_idx(top_)+ fun::top_offset(top_,flags_));
			}else if(ecut_!=_event_ && top_==_mnone_ && cut_!=_no_cut_){
				idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
				idx.push_back(fun::cut_idx(cut_));
				idx.push_back(0);
			}else if(ecut_==_none_ && cut_==_no_cut_ && top_==_mnone_){
				idx.push_back(fun::ecut_idx(ecut_));
				idx.push_back(0);
				idx.push_back(0);
			}else{
				idx.push_back(-1);
				idx.push_back(-1);
				idx.push_back(-1);
			}
			if(flags_->Sim() || weight_==_nweighted_){
				idx.push_back(fun::weight_idx(weight_));
			}else{
				idx.push_back(-1);
			}
			idx.push_back(fun::recon_idx(recon_));
		}else{
			idx.push_back(-1);
			idx.push_back(-1);
			idx.push_back(-1);
			idx.push_back(-1);
			idx.push_back(-1);
		}
	/*}else{
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
	}*/
	//std::cout<<"WQ2 idx: ";
	//fun::print_vector_idx(idx);
	return idx;
}

bool Histogram::Made_WQ2_idx(const char* ecut_, const char * cut_, const char* top_, const char * weight_, const char * recon_){
	if(_WQ2_made.size() > fun::ecut_idx(ecut_) && fun::ecut_idx(ecut_)>=0){
		if(_WQ2_made[fun::ecut_idx(ecut_)].size() > fun::cut_idx(cut_) && fun::cut_idx(cut_)>=0){
			if(_WQ2_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)].size() > fun::top_idx(top_) && fun::top_idx(top_)>=0){
				if(_WQ2_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)].size() > fun::weight_idx(weight_) && fun::weight_idx(weight_)>=0){
					if(_WQ2_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::weight_idx(weight_)].size() > fun::recon_idx(recon_) && fun::recon_idx(recon_)>=0){
						if(_WQ2_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::weight_idx(weight_)][fun::recon_idx(recon_)]){
							//std::cout<<"Made hist for: "<<ecut_ <<" " <<cut_ <<" "<<top_ <<" "<<weight_ <<" "<<recon_ <<"\n";
						}
						return _WQ2_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::weight_idx(weight_)][fun::recon_idx(recon_)];
					}else{
						std::cout<<"Looked for a recon not there: "<<recon_<<"\n";
					}
				}else{
					std::cout<<"Looked for a weight not there: "<<weight_<<"\n";
				}
			}else{
				std::cout<<"Looked for a top not there: " <<top_ <<"\n";
			}
		}else{
			std::cout<<"Looked for a cut not there: "<<cut_ <<"\n";
		}
	}else{
		std::cout<<"Looked for a ecut not there: "<<ecut_ <<"\n";
	}
	
}

void Histogram::WQ2_Fill(float W_, float Q2_, const char* ecut_, const char *cut_, const char * top_, const char * thrown_, std::shared_ptr<Flags> flags_, float weight_){
	if(!flags_->Plot_WQ2()){return;}
	if(isnan(W_)){return;}
	if(isnan(Q2_)){return;}
	//std::cout<<"** Trying to fill WQ2: " <<ecut_ <<" " <<cut_ <<" " <<top_ <<" " <<thrown_ <<"\n";
	std::vector<int> idx;
	if(!fun::top_perform(top_,flags_)){return;}
	if(flags_->Flags::Sim()){//Weighted
		idx = Histogram::WQ2_cut_idx(ecut_,cut_,top_,_weighted_,thrown_,flags_);
		if(Histogram::OK_Idx(idx)){
			if(ecut_==_event_){
				//if(thrown_ == _thrown_){
					//std::cout<<"Filling WQ2 Yield thrown " <<_WQ2_yield_hist.size()-1 <<"\n";
					//_WQ2_yield_hist[_WQ2_yield_hist.size()-1]->Fill(W_,Q2_,1.0);
				//}else{
				if(thrown_ != _thrown_){
					//std::cout<<"Filling WQ2 Yield " <<fun::top_idx(top_)+fun::top_offset(top_,flags_) <<"\n";
					_WQ2_yield_hist[fun::top_idx(top_)+fun::top_offset(top_,flags_)]->Fill(W_,Q2_,1.0);
					//_WQ2_yield_hist[fun::top_idx(top_)+fun::top_offset(_mall_,flags_)]->Fill(W_,Q2_,1.0);
				}
			}
			//if(Made_WQ2_idx(ecut_,cut_,top_,_nweighted_,thrown_)){	
			//if(Histogram::Good_WQ2_idx(ecut_,cut_,top_,_nweighted_,_recon_[fun::truth_idx(thrown_)],flags_)){
			//std::cout<<"Trying to fill WQ2_hist at \n";
			//fun::print_vector_idx(idx);
			_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(W_,Q2_,weight_);
			//}
		}
		idx.clear();
	}
	//if(Histogram::OK_Idx(idx)){
	//if(Made_WQ2_idx(ecut_,cut_,top_,_nweighted_,thrown_)){
	idx = Histogram::WQ2_cut_idx(ecut_,cut_,top_,_nweighted_,thrown_,flags_);
	if(Histogram::OK_Idx(idx)){
		if(!flags_->Flags::Sim() && ecut_==_event_){
			//std::cout<<"Filling WQ2 Yield " <<fun::top_idx(top_)+fun::top_offset(top_,flags_) <<"\n";
			_WQ2_yield_hist[fun::top_idx(top_)+fun::top_offset(top_,flags_)]->Fill(W_,Q2_,1.0);
			//_WQ2_yield_hist[fun::top_idx(top_)+fun::top_offset(_mall_,flags_)]->Fill(W_,Q2_,1.0);
		}
		//std::cout<<"Trying to fill WQ2_hist at \n";
		//fun::print_vector_idx(idx);
		_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(W_,Q2_,1.0);
	}
	idx.clear();
	//}
		
		//std::cout<<"input to WQ2 Filling: " <<ecut_ <<"\n\t" <<cut_ <<"\n\t" <<top_ <<"\n\t" <<weight_ <<"\n\t" <<_recon_[fun::truth_idx(thrown_)]  <<"\n";
		//std::cout<<"Trying to fill WQ2:" <<ecut_ <<" " <<cut_ <<" " <<top_ <<" " <<_recon_[fun::truth_idx(thrown_)] <<"\n";
		
	
}

void Histogram::WQ2_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Plot_WQ2()){
		std::cout<<"Writing WQ2\n";
		char dir_name[100];
		TDirectory* dir_WQ2 = _RootOutputFile->mkdir("W vs. Q^{2}");
		dir_WQ2->cd();
		TDirectory* dir_WQ2_sub[2];
		sprintf(dir_name,"W Q2 Ele Cuts");
		dir_WQ2_sub[0] = dir_WQ2->mkdir(dir_name);
		sprintf(dir_name,"W Q2 Event Selection");
		dir_WQ2_sub[1] = dir_WQ2->mkdir(dir_name);
		TDir_ptr_1d dir_1d;
		TDir_ptr_2d dir_2d;
		TDir_ptr_3d dir_3d;

		std::vector<int> idx;
		std::vector<long> space_dims(5);
		space_dims[4] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));//Electron Cuts
		space_dims[3] = std::distance(std::begin(_cut_), std::end(_cut_));//Cut and anti cut
		space_dims[2] = std::distance(std::begin(_top_), std::end(_top_));//{Topologies  + all combined + None}
		space_dims[1] = std::distance(std::begin(_weight_), std::end(_weight_));//Not Weight vs. Weighted
		space_dims[0] = std::distance(std::begin(_recon_), std::end(_recon_));//Reconstructed vs. Thrown
		CartesianGenerator cart(space_dims);

		dir_WQ2->cd();
		for(int i=0; i<_WQ2_yield_hist.size(); i++){
			_WQ2_yield_hist[i]->SetXTitle("W (GeV)");
			_WQ2_yield_hist[i]->SetYTitle("Q2 (GeV^2)");
			_WQ2_yield_hist[i]->Write();
		}
		while(cart.GetNextCombination()){
			if(_cut_[cart[3]]!=_no_cut_){	
				idx = Histogram::WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_);
				if(Made_WQ2_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]])){
			//if(Histogram::Good_WQ2_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_)){	
					//std::cout<<"OK index: ";
					//fun::print_vector_idx(idx);
					//std::cout<<"\n";
					if(_ecuts_[cart[4]]==_event_){//Event Selection
						//std::cout<<"Event Writing\n";
						//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						dir_WQ2_sub[1]->cd();
						//std::cout<<"Went into sub directory1\n";
						if(_WQ2_hist.size()>idx[0]){
							if(_WQ2_hist[idx[0]].size()>idx[1]){
								if(_WQ2_hist[idx[0]][idx[1]].size()>idx[2]){
									if(_WQ2_hist[idx[0]][idx[1]][idx[2]].size()>idx[3]){
										if(_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]].size()>idx[4]){
											_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("W (GeV)");
											_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Q^{2} (GeV^{2}");
											_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetOption("Colz");
											_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
										}else{
											std::cout<<"Exceeded recon: " <<_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]].size() <<" tried: " <<_recon_[idx[4]]<<"\n";
										}
									}else{
										std::cout<<"Exceeded weight: " <<_WQ2_hist[idx[0]][idx[1]][idx[2]].size() <<" tried: " <<_weight_[idx[3]]<<"\n";
									}
								}else{
									std::cout<<"Exceeded top "<<_WQ2_hist[idx[0]][idx[1]].size() <<" tried: " <<_top_[idx[2]]<<"\n";
								}
							}else{
								std::cout<<"Exceeded cut: "<<_WQ2_hist[idx[0]].size() <<" tried: " <<_cut_[idx[1]]<<"\n";
							}
						}else{
							std::cout<<"Exceeded ecut: "<<_WQ2_hist.size() <<" tried: " <<_ecuts_[idx[0]]<<"\n";
						}
					}else{//PID Plots
						//std::cout<<"PID Writing\n";
						//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
						dir_WQ2_sub[0]->cd();
						//std::cout<<"Went into sub directory0\n";
						if(_WQ2_hist.size()>idx[0]){
							if(_WQ2_hist[idx[0]].size()>idx[1]){
								if(_WQ2_hist[idx[0]][idx[1]].size()>idx[2]){
									if(_WQ2_hist[idx[0]][idx[1]][idx[2]].size()>idx[3]){
										if(_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]].size()>idx[4]){
											_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("W (GeV)");
											_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Q^{2} (GeV^{2}");
											_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetOption("Colz");
											_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
										}else{
											std::cout<<"Exceeded recon: " <<_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]].size() <<" tried: " <<_recon_[idx[4]]<<"\n";
										}
									}else{
										std::cout<<"Exceeded weight: " <<_WQ2_hist[idx[0]][idx[1]][idx[2]].size() <<" tried: " <<_weight_[idx[3]]<<"\n";
									}
								}else{
									std::cout<<"Exceeded top "<<_WQ2_hist[idx[0]][idx[1]].size() <<" tried: " <<_top_[idx[2]]<<"\n";
								}
							}else{
								std::cout<<"Exceeded cut: "<<_WQ2_hist[idx[0]].size() <<" tried: " <<_cut_[idx[1]]<<"\n";
							}
						}else{
							std::cout<<"Exceeded ecut: "<<_WQ2_hist.size() <<" tried: " <<_ecuts_[idx[0]]<<"\n";
						}
						
					}
				/*}else{
					std::cout<<"\tWasn't a Made index for WQ2\n";
					fun::print_vector_idx(idx);
				}*/
				}else{
					//std::cout<<"\tWasn't a made idx for WQ2\n";
					//fun::print_vector_idx(idx);
				}
				idx.clear();
			}
		}
	}
}
 //*--------------------------------End W Q2 Plot----------------------------*

 //*-------------------------------Start Fid Plot----------------------------*

void Histogram::Fid_Make(std::shared_ptr<Flags> flags_){
	std::cout<<"Fid Plot: " <<flags_->Plot_Fid(0) <<" " <<flags_->Plot_Fid(1) <<" " <<flags_->Plot_Fid(2) <<" " <<flags_->Plot_Fid(3) <<"\n";
	if(flags_->Plot_Fid(0) || flags_->Plot_Fid(1) || flags_->Plot_Fid(2) || flags_->Plot_Fid(3)){
		std::cout<<"Making Fid Histograms\n";
		TH2F_ptr_6d plot_6d;
		TH2F_ptr_5d plot_5d;
		TH2F_ptr_4d plot_4d;
		TH2F_ptr_3d plot_3d;
		TH2F_ptr_2d plot_2d;
		TH2F_ptr_1d plot_1d;
		Bool_5d made_5d;
		Bool_4d made_4d;
		Bool_3d made_3d;
		Bool_2d made_2d;
		Bool_1d made_1d;
		std::vector<long> space_dims(7);
		space_dims[6] = std::distance(std::begin(_weight_), std::end(_weight_));//Not Weight vs. Weighted;//_p_bin_bins_+1;
		space_dims[5] = std::distance(std::begin(_species_), std::end(_species_));//{ele,pro,pip,pim}
		space_dims[4] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));//Electron Cuts (electrons have more cuts than hadrons)
		space_dims[3]	= std::distance(std::begin(_sector_), std::end(_sector_));//Sectors
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));//Cut and anti cut
		space_dims[1] = std::distance(std::begin(_top_), std::end(_top_));//{Topologies  + all combined + None}
		space_dims[0] = _p_bin_bins_+1;//std::distance(std::begin(_weight_), std::end(_weight_));//Not Weight vs. Weighted
		CartesianGenerator cart(space_dims);
		//CartesianGenerator cart2(space_dims);
		char hname[100];//For naming histograms
		int idx = 0; 


		while(cart.GetNextCombination()){
			if(flags_->Plot_Fid(cart[5]) && (flags_->Flags::Sim() || _weight_[cart[6]]==_nweighted_)){
				if(_species_[cart[5]] == _ele_){
					if(_ecuts_[cart[4]] == _event_ && _top_[cart[1]]!=_mnone_ && _cut_[cart[2]]!=_no_cut_ && fun::top_perform(_top_[cart[1]],flags_)){
						if(cart[0]==_p_bin_bins_){
							sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[6]],_top_[cart[1]]);

							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
						}else{
							sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
						}
					}else if(_ecuts_[cart[4]] != _event_ && _top_[cart[1]]==_mnone_ ){//No Event Selection
						if(_ecuts_[cart[4]]== _none_ && _cut_[cart[2]]==_no_cut_){// && _cut_[cart[2]]==_no_cut_){
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[6]]);
								//std::cout<<"hname:" <<hname <<"\n";
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
						}
						if(_ecuts_[cart[4]]== _sanity_ && _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
						if(_ecuts_[cart[4]]== _fid_cut_ && flags_->Flags::Fid_Cut(0)&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
						if(_ecuts_[cart[4]]== _sf_cut_ && flags_->Flags::SF_Cut()&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
						if(_ecuts_[cart[4]]== _cc_cut_ && flags_->Flags::CC_Cut()&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
						if(_ecuts_[cart[4]]== _ec_cut_ && flags_->Flags::EC_Cut()&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
						if(_ecuts_[cart[4]]== _vertex_cut_ && flags_->Flags::Vertex_Cut()&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
						if(_ecuts_[cart[4]]== _geo_cc_cut_ && flags_->Flags::Geo_Cut(cart[5],0)&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
						if(_ecuts_[cart[4]]== _geo_sc_cut_ && flags_->Flags::Geo_Cut(cart[5],1)&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
						if(_ecuts_[cart[4]]== _geo_ec_cut_ && flags_->Flags::Geo_Cut(cart[5],2)&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
						if(_ecuts_[cart[4]]== _id_cut_ && flags_->Flags::ID_Cut()&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
						if(_ecuts_[cart[4]]== _pid_ && _cut_[cart[2]]!=_no_cut_){
						//	std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}
					}
				}else{// if((_species_[cart[5]]!=_pim_ && cart[0]==_p_bin_bins_)||(_species_[cart[5]]==_pim_))
					if(cart[4]<(std::distance(std::begin(_hcuts_), std::end(_hcuts_)))){//Make sure we aren't going over hcuts
						//std::cout<<"\tFid for Hadrons\n";
						if(_hcuts_[cart[4]] == _event_ && _top_[cart[1]]!=_mnone_ && _cut_[cart[2]]!=_no_cut_ && fun::top_perform(_top_[cart[1]],flags_)){
							//std::cout<<"\tFid Event for Hadrons\n";
							//std::cout<<"\t\tFid Fid Hadrons " <<_weight_[ cart[6]] <<" "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" pbin:" <<cart[0] <<"\n";
							//std::cout<<"\t\t\tFid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							if(cart[0]==_p_bin_bins_){
								if( _species_[cart[5]]==_pim_){
									sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]]);
									plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								}
							}else{
								sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							}
							//sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]]);
							//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
						}else if(_hcuts_[cart[4]] != _event_ && _top_[cart[1]]==_mnone_){//No Event Selection
							//std::cout<<"\tFid PID Hadrons\n";
							if(_hcuts_[cart[4]]== _none_ && _cut_[cart[2]]==_no_cut_ ){
								//std::cout<<"\t\tFid None Hadrons "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_weight_[ cart[6]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\tFid Fid Hadrons " <<_weight_[ cart[6]] <<" "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\t\tFid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								if(cart[0]==_p_bin_bins_){
									if( _species_[cart[5]]==_pim_){
										sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
										plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
									}
								}else{
									sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
									plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								}
								//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
							}
							if(_hcuts_[cart[4]]== _sanity_ && _cut_[cart[2]]!=_no_cut_){
								//std::cout<<"\t\tFid Sanity Hadrons "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_weight_[ cart[6]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\tFid Fid Hadrons " <<_weight_[ cart[6]] <<" "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\t\tFid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								if(cart[0]==_p_bin_bins_){
									if( _species_[cart[5]]==_pim_){
										sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
										plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
									}
								}else{
									sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
									plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								}
								//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
							}
							if(_hcuts_[cart[4]]== _fid_cut_ && flags_->Flags::Fid_Cut(cart[5]) && _cut_[cart[2]]!=_no_cut_){
								//std::cout<<"\t\tFid Fid Hadrons " <<_weight_[ cart[6]] <<" "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\t\tFid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								if(cart[0]==_p_bin_bins_){
									if( _species_[cart[5]]==_pim_){
										sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
										plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
									}
								}else{
									sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
									plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								}
								//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
							}
							if(_hcuts_[cart[4]]== _delta_cut_ && flags_->Flags::Delta_Cut(cart[5]) && _cut_[cart[2]]!=_no_cut_){
								//std::cout<<"\t\tFid delta Hadrons " <<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_weight_[ cart[6]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\tFid Fid Hadrons " <<_weight_[ cart[6]] <<" "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\t\tFid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								if(cart[0]==_p_bin_bins_){
									if( _species_[cart[5]]==_pim_){
										sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
										plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
									}
								}else{
									sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
									plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								}
								//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
							}
							if(_hcuts_[cart[4]]== _geo_sc_cut_ && flags_->Flags::Geo_Cut(cart[5],1) && _cut_[cart[2]]!=_no_cut_){
								//std::cout<<"\t\tFid delta Hadrons " <<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_weight_[ cart[6]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\tFid Fid Hadrons " <<_weight_[ cart[6]] <<" "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\t\tFid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								if(cart[0]==_p_bin_bins_){
									if( _species_[cart[5]]==_pim_){
										sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
										plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
									}
								}else{
									sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
									plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								}
								//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
							}
							if(_hcuts_[cart[4]]== _geo_ec_cut_ && flags_->Flags::Geo_Cut(cart[5],2) && _cut_[cart[2]]!=_no_cut_){
								//std::cout<<"\t\tFid delta Hadrons " <<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_weight_[ cart[6]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\tFid Fid Hadrons " <<_weight_[ cart[6]] <<" "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\t\tFid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								if(cart[0]==_p_bin_bins_){
									if( _species_[cart[5]]==_pim_){
										sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
										plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
									}
								}else{
									sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
									plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								}
								//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
							}
							if(_hcuts_[cart[4]]== _pid_  && _cut_[cart[2]]!=_no_cut_){
								///std::cout<<"\t\tFid PID "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_weight_[ cart[6]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\tFid Fid Hadrons " <<_weight_[ cart[6]] <<" "<<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" pbin:" <<cart[0] <<"\n";
								//std::cout<<"\t\t\tFid Hist idx: " <<_Fid_hist.size() <<" " <<plot_6d.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								if(cart[0]==_p_bin_bins_){
									if(_species_[cart[5]]==_pim_){
										sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
										plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
									}
								}else{
									sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s_p:%.3f-%.3f",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[6]],_top_[cart[1]],Histogram::P_Min(cart[0]),Histogram::P_Max(cart[0]));
									plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								}
								//sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[ cart[6]]);
								//plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								//_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][ cart[6]] = true;
							}
						}
					}
				}
			}
			if(cart[0] == space_dims[0]-1){
				if(plot_1d.size()>0){
					plot_2d.push_back(plot_1d);
					plot_1d.clear();
				}
				if(cart[1] == space_dims[1]-1){ 
					if(plot_2d.size()>0){
						plot_3d.push_back(plot_2d);
						plot_2d.clear();
					}
					if(cart[2] == space_dims[2]-1){
						if(plot_3d.size()>0){
							plot_4d.push_back(plot_3d);
							plot_3d.clear();
						}
						if(cart[3] == space_dims[3]-1){
							if(plot_4d.size()>0){
								plot_5d.push_back(plot_4d);
								plot_4d.clear();
							}
							if(cart[4] == space_dims[4]-1){
								if(plot_5d.size()>0){
									plot_6d.push_back(plot_5d);
									plot_5d.clear();
								}
								if(cart[5] == space_dims[5]-1){
									if(plot_6d.size()>0){
										_Fid_hist.push_back(plot_6d);
										plot_6d.clear();
									}
									
								}
							}
						}
					}
				}
			}
		}
	}
}

int Histogram::hcut_idx(char * hcut_, int hadron_, std::shared_ptr<Flags> flags_){
	int idx = -1; 
	int idx_tmp = -1; 
	for(int i = 0; i<std::distance(std::begin(_hcuts_), std::end(_hcuts_)); i++){
		if(hcut_ == _hcuts_[i]){
			idx_tmp += i+1;
		}
	}
	switch(idx_tmp){
		case 0://None
			idx = idx_tmp;
		break;
		case 1://Sanity
			if(flags_->Flags::Fid_Cut(hadron_)){
				idx = idx_tmp; 
			}
		break;
		case 2://Fid
			if(flags_->Flags::Delta_Cut(hadron_)){
				if(flags_->Flags::Fid_Cut(hadron_)){
					idx = idx_tmp; 
				}else{
					idx = idx_tmp-1; 
				}
			}
		break;
		case 3://Delta t
			if(flags_->Flags::Delta_Cut(hadron_)){
				if(flags_->Flags::Fid_Cut(hadron_)){
					idx = idx_tmp; 
				}else{
					idx = idx_tmp-1; 
				}
			}
		break;
		case 4://sc_eff cut
			if(flags_->Flags::SC_Eff()){
				if(flags_->Flags::Delta_Cut(hadron_)){
					if(flags_->Flags::Fid_Cut(hadron_)){
						idx = idx_tmp; 
					}else{
						idx = idx_tmp-1; 
					}
				}else{
					if(flags_->Flags::Fid_Cut(hadron_)){
						idx = idx_tmp-1; 
					}else{
						idx = idx_tmp-2; 
					}
				}
			}
		break;
		case 5://id cut
			if(flags_->Flags::ID_Cut()){
				if(flags_->Flags::SC_Eff()){
					if(flags_->Flags::Delta_Cut(hadron_)){
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp; 
						}else{
							idx = idx_tmp-1; 
						}
					}else{
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp-1; 
						}else{
							idx = idx_tmp-2; 
						}
					}
				}
			}else{
				if(flags_->Flags::SC_Eff()){
					if(flags_->Flags::Delta_Cut(hadron_)){
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp-1; 
						}else{
							idx = idx_tmp-2; 
						}
					}else{
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp-2; 
						}else{
							idx = idx_tmp-3; 
						}
					}
				}
			}
		break;
		case 6://pid
			if(flags_->Flags::ID_Cut()){
				if(flags_->Flags::SC_Eff()){
					if(flags_->Flags::Delta_Cut(hadron_)){
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp; 
						}else{
							idx = idx_tmp-1; 
						}
					}else{
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp-1; 
						}else{
							idx = idx_tmp-2; 
						}
					}
				}
			}else{
				if(flags_->Flags::SC_Eff()){
					if(flags_->Flags::Delta_Cut(hadron_)){
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp-1; 
						}else{
							idx = idx_tmp-2; 
						}
					}else{
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp-2; 
						}else{
							idx = idx_tmp-3; 
						}
					}
				}
			}
		break;
		case 7://event
			if(flags_->Flags::ID_Cut()){
				if(flags_->Flags::SC_Eff()){
					if(flags_->Flags::Delta_Cut(hadron_)){
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp; 
						}else{
							idx = idx_tmp-1; 
						}
					}else{
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp-1; 
						}else{
							idx = idx_tmp-2; 
						}
					}
				}
			}else{
				if(flags_->Flags::SC_Eff()){
					if(flags_->Flags::Delta_Cut(hadron_)){
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp-1; 
						}else{
							idx = idx_tmp-2; 
						}
					}else{
						if(flags_->Flags::Fid_Cut(hadron_)){
							idx = idx_tmp-2; 
						}else{
							idx = idx_tmp-3; 
						}
					}
				}
			}
		break;
		default:
			std::cout<<"\tInvalid HID index\n";
		break;
	}
	return idx; 
}
int Histogram::hcut_idx(const char * hcut_, int hadron_, std::shared_ptr<Flags> flags_){
	int idx = fun::hcut_idx(hcut_); 
	int idx_tmp = fun::hcut_idx(hcut_);
	//std::cout<<"hcut_idx try: " <<hcut_ <<" "<<idx_tmp <<"\n";
	for(int i=0; i<idx_tmp; i++){
		if(_hcuts_[i] == _sanity_){
		}else if(_hcuts_[i] == _fid_cut_){
			if(!flags_->Flags::Fid_Cut(hadron_)){
				idx += -1; 
			}
		}else if(_hcuts_[i] == _delta_cut_){
			if(!flags_->Flags::Delta_Cut(hadron_)){
				idx += -1; 
			}
		}
	}
	return idx; 
}
	

/*
bool Histogram::Made_Fid_idx(const char* species_, const char* pcut_, const char* sector_, const char * cut_, const char * top_, const char * weight_){
	if(species_==_ele_){
		return _Fid_made[fun::species_idx(species_)][fun::ecut_idx(pcut_)+fun::ecut_offset(pcut_,flags_)][fun::sector_idx(sector_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::weight_idx(weight_)];
	}else{
		return _Fid_made[fun::species_idx(species_)][fun::hcut_idx(pcut_)+fun::hcut_offset(species_,pcut_,flags_)][fun::sector_idx(sector_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::weight_idx(weight_)];
	}
}
*/

std::vector<int> Histogram::Fid_cut_idx(const char* species_, const char* pcut_, const char* sector_, const char * cut_, const char * top_, const char * weight_, const char * p_dep_, float p_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;// = {-1,-1,-1,-1,-1};
	int species_idx= fun::species_idx(species_)+fun::species_offset(species_,_fid_cut_,flags_); 
	int pcut_idx = -1;  
	int cut_idx = -1;
	//idx [pbin][species][pcut][sector][cut][top]
	//Weight
	if(flags_->Flags::Sim() || weight_==_nweighted_){
		idx.push_back(fun::weight_idx(weight_));
	}else{
		idx.push_back(-1);
	}
	if(species_ == _ele_ && flags_->Flags::Plot_Fid(0)){
		pcut_idx = fun::ecut_idx(pcut_)+fun::ecut_offset(pcut_,flags_);
	}else if(species_ != _ele_ && flags_->Flags::Plot_Fid(fun::species_idx(species_))){
		pcut_idx = fun::hcut_idx(pcut_)+fun::hcut_offset(species_,pcut_,flags_);
	}
	if(pcut_idx >= 0 && fun::pcut_perform(species_,pcut_,flags_)){
		//Species
		idx.push_back(species_idx);
		/*if(species_==_pim_ && p_dep_==_p_dep_){
			idx.push_back(1);
		}else{
			idx.push_back(species_idx);
		}*/
		//PID Cut
		if(pcut_ == _event_ && top_!=_mnone_ && cut_!=_no_cut_){
			idx.push_back(pcut_idx);
		}else if(pcut_ != _event_ && top_==_mnone_ && cut_!=_no_cut_){
			idx.push_back(pcut_idx);
		}else if(pcut_ == _none_ && cut_==_no_cut_){
			idx.push_back(pcut_idx);
		}else{
			idx.push_back(-1);
		}
		//Sector
		idx.push_back(fun::sector_idx(sector_));
		//Cut or Anti
		if(pcut_!=_none_ && cut_!=_no_cut_){
			idx.push_back(fun::cut_idx(cut_));
		}else if(pcut_==_none_ && cut_==_no_cut_){
			idx.push_back(0);
		}else{
			idx.push_back(-1);
		}
		
		//Topology
		if(pcut_==_event_ && top_!=_mnone_){
			idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
		}else if(pcut_!=_event_ && top_==_mnone_){
			idx.push_back(0);
		}else{
			idx.push_back(-1);
		}
		if((species_==_ele_ || species_==_pim_)&& p_dep_==_p_dep_){
			idx.push_back(Histogram::P_bin(p_));
		}else if(p_dep_ == _no_p_dep_){
			if(species_==_ele_ || species_ == _pim_){
				idx.push_back(_p_bin_bins_);
			}else{
				if(flags_->Flags::Plot_Fid(3)){
					idx.push_back(_p_bin_bins_);
				}else{
					idx.push_back(0);
				}
			}
		}else if(species_==_ele_ || species_ == _pim_){
			idx.push_back(Histogram::P_bin(p_));
		}else{
			idx.push_back(-1);
		}
	}else{
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
	}
	//std::cout<<"\t\tFid idx";
	//if(top_!=_mnone_ && cut_ == _cut_applied_){
	//	std::cout<<"\tFid Idx being called: " <<weight_ <<" " <<species_ <<" " <<pcut_ <<" " <<sector_ <<" " <<cut_ <<" " <<top_ <<" "  <<p_dep_ <<" " <<p_ <<"\n";
	//	fun::print_vector_idx(idx);
	//}
	return idx;
}

void Histogram::Fid_Fill(const char * species_, float theta_, float phi_, const char* pcut_, const char* sector_, const char *cut_, const char * top_, float p_, std::shared_ptr<Flags> flags_, float weight_){
	if(flags_->Plot_Fid(fun::species_idx(species_)) && !isnan(theta_) && !isnan(phi_)){
		//std::cout<<"Trying to fill Fid:" <<species_ <<" " <<pcut_ <<" "<<sector_ <<" " <<cut_ <<" " <<top_ <<" " <<weight_ <<" " <<p_ <<"\n";
		int w_idx = -1;
		if(weight_!=1.0){
			w_idx = 1;
		}else{
			w_idx = 0;
		}
		std::vector<int> idx;
		if(flags_->Flags::Sim()){//Weighted Fiducial Plots
			idx = Histogram::Fid_cut_idx(species_,pcut_,sector_,cut_,top_,_weighted_,_no_p_dep_,p_,flags_);;
			if(Histogram::OK_Idx(idx)){
				//if(Histogram::Made_Fid_idx(species_,pcut_,sector_,cut_,top_,_weighted_)){
					//fun::print_vector_idx(idx);
					//std::cout<<"Filling " <<_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->GetName() <<"\n";
					_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Fill(phi_,theta_,weight_);
				//}
			}
			idx.clear();
			idx = Histogram::Fid_cut_idx(species_,pcut_,_sec_all_,cut_,top_,_weighted_,_no_p_dep_,p_,flags_);;
			if(Histogram::OK_Idx(idx)){
				//if(Histogram::Made_Fid_idx(species_,pcut_,_sec_all_,cut_,top_,_weighted_)){
					//fun::print_vector_idx(idx);
					//std::cout<<"Filling " <<_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->GetName() <<"\n";
					_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Fill(phi_,theta_,weight_);
				//}
			}
			idx.clear();
			if(species_==_ele_ || species_==_pim_){
				idx = Histogram::Fid_cut_idx(species_,pcut_,_sec_all_,cut_,top_,_weighted_,_p_dep_,p_,flags_);;
				if(Histogram::OK_Idx(idx)){
					//if(Histogram::Made_Fid_idx(species_,pcut_,_sec_all_,cut_,top_,_weighted_)){
						//fun::print_vector_idx(idx);
						//std::cout<<"Filling " <<_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->GetName() <<"\n";
						_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Fill(phi_,theta_,weight_);
					//}
				}
			}
		}else{
			//Non Weighted Histograms

			//No Momentum Dependence
			//Sector Specific
			idx = Histogram::Fid_cut_idx(species_,pcut_,sector_,cut_,top_,_nweighted_,_no_p_dep_,p_,flags_);
			
			//if(pcut_==_event_){
				//std::cout<<"Trying to fill Fid: phi:" <<phi_ <<" theta:" <<theta_ <<" " <<species_ <<" " <<pcut_ <<" "<<sector_ <<" " <<cut_ <<" " <<top_ <<" " <<weight_ <<"\n";
			//}
			//fun::print_vector_idx(idx);
			if(Histogram::OK_Idx(idx)){
				//if(Histogram::Made_Fid_idx(species_,pcut_,sector_,cut_,top_,_nweighted_)){
					//if(pcut_ == _event_){
					//	fun::print_vector_idx(idx);
					//	std::cout<<species_ <<" " <<pcut_ <<" " <<sector_ <<" " <<cut_ <<" " <<top_ <<" " <<_nweighted_ <<"\n";
					//}
					//fun::print_vector_idx(idx);
					//std::cout<<"Filling " <<_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->GetName() <<"\n";
					_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Fill(phi_,theta_,1.0);
				//}
			}
			idx.clear();
			//Fill All Sectors. No Momentum Dependence
			idx = Histogram::Fid_cut_idx(species_,pcut_,_sec_all_,cut_,top_,_nweighted_,_no_p_dep_,p_,flags_);
			//fun::print_vector_idx(idx);
			if(Histogram::OK_Idx(idx)){
				//if(Histogram::Made_Fid_idx(species_,pcut_,_sec_all_,cut_,top_,_nweighted_)){
					//if(pcut_ == _event_){
						//fun::print_vector_idx(idx);
					//}
					//fun::print_vector_idx(idx);
					//std::cout<<"Filling " <<_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->GetName() <<"\n";
					_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Fill(phi_,theta_,1.0);
				//}
			}
			idx.clear();
			//std::vector<int> idx = Histogram::Fid_cut_idx(species_,pcut_,cut_,top_,w_idx,flags_);
			if(species_==_ele_ || species_==_pim_){
				idx = Histogram::Fid_cut_idx(species_,pcut_,sector_,cut_,top_,_nweighted_,_p_dep_,p_,flags_);
				//fun::print_vector_idx(idx);
				if(Histogram::OK_Idx(idx)){
					//if(Histogram::Made_Fid_idx(species_,pcut_,_sec_all_,cut_,top_,_nweighted_)){
						//if(pcut_ == _event_){
							//fun::print_vector_idx(idx);
						//}
						//fun::print_vector_idx(idx);
						//std::cout<<"Filling " <<_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->GetName() <<"\n";
						_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Fill(phi_,theta_,1.0);
					//}
				}
				idx.clear();
				idx = Histogram::Fid_cut_idx(species_,pcut_,_sec_all_,cut_,top_,_nweighted_,_p_dep_,p_,flags_);
				//fun::print_vector_idx(idx);
				if(Histogram::OK_Idx(idx)){
					//if(Histogram::Made_Fid_idx(species_,pcut_,_sec_all_,cut_,top_,_nweighted_)){
						//if(pcut_ == _event_){
							//fun::print_vector_idx(idx);
						//}
						//fun::print_vector_idx(idx);
						//std::cout<<"Filling " <<_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->GetName() <<"\n";
						_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Fill(phi_,theta_,1.0);
					//}
				}
				idx.clear();
			}
		}
	}
}

void Histogram::Fid_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Plot_Fid(0) || flags_->Plot_Fid(1) || flags_->Plot_Fid(2) || flags_->Plot_Fid(3)){
		std::cout<<"Writing Fiducial Plots\n";
		char dir_name[100];
		TDirectory* dir_fid = _RootOutputFile->mkdir("Fiducial");
		dir_fid->cd();
		std::cout<<"\tMaking Directories\n";
		TDirectory* dir_fid_sub[4][2][std::distance(std::begin(_ecuts_), std::end(_ecuts_))+1][2][std::distance(std::begin(_sector_), std::end(_sector_))+1];//{species,pid/event,pid_cuts,p_dep,sector}
		for(int i=0; i<std::distance(std::begin(_species_), std::end(_species_)); i++){
			if(flags_->Flags::Plot_Fid(i)){
				sprintf(dir_name,"Fid %s PID Cuts",_species_[i]);
				dir_fid_sub[i][0][0][0][0] = dir_fid->mkdir(dir_name);
				//std::cout<<"\t\tMade Dir: " <<i <<" 0 0 0 0 " <<dir_name <<"\n ";// <<j+1 <<" 0\n";
				if(_species_[i]==_ele_){
					for(int j=0; j<std::distance(std::begin(_ecuts_), std::end(_ecuts_));j++){
						if(fun::ecut_perform(_ecuts_[j],flags_)){
							sprintf(dir_name,"Fid %s %s Cut",_species_[i],_ecuts_[j]);
							dir_fid_sub[i][0][j+1][0][0] = dir_fid_sub[i][0][0][0][0]->mkdir(dir_name);
							//std::cout<<"\t\tMade Dir: " <<i <<" 0 " <<j+1 <<" 0 " <<0 <<" " <<dir_name <<"\n";
							sprintf(dir_name,"Fid %s %s %s Cut",_species_[i],_ecuts_[j],_p_dep_);
							dir_fid_sub[i][0][j+1][1][0] = dir_fid_sub[i][0][j+1][0][0]->mkdir(dir_name);
							//std::cout<<"\t\tMade Dir: " <<i <<" 0 " <<j+1 <<" 1 " <<0 <<" " <<dir_name <<"\n";
							for(int k=0; k<std::distance(std::begin(_sector_), std::end(_sector_)); k++){
								sprintf(dir_name,"Fid %s %s %s Cut",_species_[i],_ecuts_[j],_sector_[k]);
								dir_fid_sub[i][0][j+1][0][k+1] = dir_fid_sub[i][0][j+1][0][0]->mkdir(dir_name);
								//std::cout<<"\t\tMade Dir: " <<i <<" 0 " <<j+1 <<" 0 " <<k+1 <<" " <<dir_name <<"\n";
								sprintf(dir_name,"Fid %s %s %s %s Cut",_species_[i],_ecuts_[j],_p_dep_,_sector_[k]);
								dir_fid_sub[i][0][j+1][1][k+1] = dir_fid_sub[i][0][j+1][1][0]->mkdir(dir_name);
							//	std::cout<<"\t\tMade Dir: " <<i <<" 0 " <<j+1 <<" 1 " <<k+1 <<" " <<dir_name <<"\n";
							}
						}
					}
				}else{
					for(int j=0; j<std::distance(std::begin(_hcuts_), std::end(_hcuts_));j++){
						if(fun::hcut_perform(_species_[i],_hcuts_[j],flags_)){
							sprintf(dir_name,"Fid %s %s Cut",_species_[i],_hcuts_[j]);
							dir_fid_sub[i][0][j+1][0][0] = dir_fid_sub[i][0][0][0][0]->mkdir(dir_name);
						//	std::cout<<"\t\tMade Dir: " <<i <<" 0 " <<j+1 <<" 0 " <<0 <<" " <<dir_name <<"\n";
							for(int k=0; k<std::distance(std::begin(_sector_), std::end(_sector_)); k++){
								sprintf(dir_name,"Fid %s %s %s Cut",_species_[i],_hcuts_[j],_sector_[k]);
								dir_fid_sub[i][0][j+1][0][k+1] = dir_fid_sub[i][0][j+1][0][0]->mkdir(dir_name);
								//std::cout<<"\t\tMade Dir: " <<i <<" 0 " <<j+1 <<" 0 " <<k+1 <<" " <<dir_name <<"\n";
							}
							if(_species_[i]==_pim_){
								sprintf(dir_name,"Fid %s %s %s Cut",_species_[i],_hcuts_[j],_p_dep_);
								dir_fid_sub[i][0][j+1][1][0] = dir_fid_sub[i][0][j+1][0][0]->mkdir(dir_name);
								//std::cout<<"Made Dir: " <<i <<" 0 " <<j+1 <<" 1 " <<0 <<"\n";
								for(int k=0; k<std::distance(std::begin(_sector_), std::end(_sector_)); k++){
									sprintf(dir_name,"Fid %s %s %s %s Cut",_species_[i],_hcuts_[j],_p_dep_,_sector_[k]);
									dir_fid_sub[i][0][j+1][1][k+1] = dir_fid_sub[i][0][j+1][1][0]->mkdir(dir_name);
									//std::cout<<"\t\tMade Dir: " <<i <<" 0 " <<j+1 <<" 1 " <<k+1 <<"\n";
								}
							}
						}
					}
				}
				sprintf(dir_name,"Fid Event %s Selection",_species_[i]);
				dir_fid_sub[i][1][0][0][0] = dir_fid->mkdir(dir_name);
				for(int j=0; j<5; j++){
					if(fun::top_perform(_top_[j],flags_)){
						sprintf(dir_name,"Fid Event %s %s Selection",_species_[i],_top_[j]);
						dir_fid_sub[i][1][j+1][0][0] = dir_fid_sub[i][1][0][0][0]->mkdir(dir_name);
						//std::cout<<"\t\tMade Dir: " <<i <<" 1 " <<j+1 <<" 0 " <<0 <<" " <<dir_name <<"\n";
						for(int k=0; k<std::distance(std::begin(_sector_), std::end(_sector_)); k++){
							sprintf(dir_name,"Fid Event %s %s %s Selection",_species_[i],_top_[j],_sector_[k]);
							dir_fid_sub[i][1][j+1][0][k+1] = dir_fid_sub[i][1][j+1][0][0]->mkdir(dir_name);
							//std::cout<<"\t\tMade Dir: " <<i <<" 1 " <<j+1 <<" 0 " <<k+1 <<" " <<dir_name <<"\n";
						}
						if(_species_[i]==_ele_ || _species_[i]==_pim_){
							sprintf(dir_name,"Fid Event %s %s %s Selection",_species_[i],_top_[j],_p_dep_);
							dir_fid_sub[i][1][j+1][1][0] = dir_fid_sub[i][1][0][0][0]->mkdir(dir_name);
							std::cout<<"\t\tMade Dir: " <<i <<" 1 " <<j+1 <<" 1 " <<0 <<" " <<dir_name <<"\n";
							for(int k=0; k<std::distance(std::begin(_sector_), std::end(_sector_)); k++){
								sprintf(dir_name,"Fid Event %s %s %s %s Selection",_species_[i],_top_[j],_p_dep_,_sector_[k]);
								dir_fid_sub[i][1][j+1][1][k+1] = dir_fid_sub[i][1][j+1][1][0]->mkdir(dir_name);
								//std::cout<<"\t\tMade Dir: " <<i <<" 1 " <<j+1 <<" 1 " <<k+1 <<" " <<dir_name <<"\n";
							}
						}
					}
				}
			}
		}
		std::cout<<"\tDirectories Made\n";

		std::vector<long> space_dims(7);
		space_dims[0] = _p_bin_bins_+1;
		space_dims[5] = std::distance(std::begin(_species_), std::end(_species_));//{ele,pro,pip,pim}
		space_dims[4] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));//Electron Cuts (electrons have more cuts than hadrons)
		space_dims[3] = std::distance(std::begin(_sector_), std::end(_sector_));
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));//Cut and anti cut
		space_dims[1] = std::distance(std::begin(_top_), std::end(_top_));//{Topologies  + all combined + None}
		space_dims[6] = std::distance(std::begin(_weight_), std::end(_weight_));//Not Weight vs. Weighted
		CartesianGenerator cart(space_dims);

		std::vector<int> idx;
		std::cout<<"\tWriting Fid Histograms\n";
		while(cart.GetNextCombination()){
			if(flags_->Flags::Plot_Fid(cart[5])){
				if(_cut_[cart[2]]!=_no_cut_){
					if(_species_[cart[5]] == _ele_){
						if(cart[0]==_p_bin_bins_){
							idx = Histogram::Fid_cut_idx(_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[6]],_no_p_dep_,Histogram::P_center(cart[0]),flags_);
						}else{
							idx = Histogram::Fid_cut_idx(_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[6]],_p_dep_,Histogram::P_center(cart[0]),flags_);
						}
						//idx = Histogram::Fid_cut_idx(_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[6]],_p_look_[cart[0]],Histogram::P_center(cart[0]),flags_);
						//if(Histogram::OK_Idx(idx)){
						if(Histogram::OK_Idx(idx)){ //&& Histogram::Made_Fid_idx(_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[6]])){
							//fun::print_vector_idx(idx);
							//idx = Histogram::W2_cut_idx(_ecuts_[cart[5]],_cut_[cart[4]],_top_[cart[2]],cart[6],flags_,cart[1]);
							if(_ecuts_[cart[4]] == _event_){//Event Selection
								if(cart[0]==_p_bin_bins_){
									//std::cout<<"Looking at Fid Dir: " <<cart[5] <<" 1 " <<cart[1]+1 <<" 0 " <<cart[3]+1 <<"\n";
									dir_fid_sub[cart[5]][1][cart[1]+1][0][cart[3]+1]->cd();
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetXTitle("W (GeV)");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetYTitle("Q^{2} (GeV^{2}");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetOption("Colz");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Write();
								}else{
									//std::cout<<"Looking at Fid Dir: " <<cart[5] <<" 1 " <<cart[1]+1 <<" 1 " <<cart[3]+1 <<"\n";
									dir_fid_sub[cart[5]][1][cart[1]+1][1][cart[3]+1]->cd();
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetXTitle("W (GeV)");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetYTitle("Q^{2} (GeV^{2}");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetOption("Colz");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Write();
								}
							}else{//PID Plots
								if(cart[0]==_p_bin_bins_){
									//std::cout<<"Looking at Fid Dir: " <<cart[5] <<" 0 " <<cart[4]+1 <<" 0 " <<cart[3]+1 <<"\n";
									dir_fid_sub[cart[5]][0][cart[4]+1][0][cart[3]+1]->cd();
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetXTitle("W (GeV)");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetYTitle("Q^{2} (GeV^{2}");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetOption("Colz");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Write();
								}else{
									//std::cout<<"Looking at Fid Dir: " <<cart[5] <<" 0 " <<cart[4]+1 <<" 1 " <<cart[3]+1 <<"\n";
									dir_fid_sub[cart[5]][0][cart[4]+1][1][cart[3]+1]->cd();
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetXTitle("W (GeV)");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetYTitle("Q^{2} (GeV^{2}");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetOption("Colz");
									_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Write();
								}
							}
						}
					}else{
						if(cart[4] < sizeof(_hcuts_) && Histogram::OK_Idx(idx)){// && Histogram::Made_Fid_idx(_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[0]])){
							if((_hcuts_[cart[4]]==_event_ && _top_[cart[1]]!=_mnone_ && ((_species_[cart[5]]==_pro_ && _top_[cart[1]]!=_mpro_) || (_species_[cart[5]]==_pip_ && _top_[cart[1]]!=_mpip_) || (_species_[cart[5]]==_pim_ && _top_[cart[1]]!=_mpim_))) || (_hcuts_[cart[4]]!=_event_ && _top_[cart[1]]==_mnone_)){
								if(cart[0]==_p_bin_bins_){
									idx = Histogram::Fid_cut_idx(_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[6]],_no_p_dep_,Histogram::P_center(cart[0]),flags_);
								}else if(_species_[cart[5]]==_pim_){
									idx = Histogram::Fid_cut_idx(_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[6]],_p_dep_,Histogram::P_center(cart[0]),flags_);
								}
								if(Histogram::OK_Idx(idx)){
									if(_hcuts_[cart[4]]==_event_){//Event Selection
										if(cart[0]==_p_bin_bins_){
											//std::cout<<"writing hist: " <<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" " <<_weight_[cart[6]] <<" " <<_no_p_dep_ <<"\n";
											//std::cout<<"Looking at Fid Dir: " <<cart[5] <<" 1 " <<cart[1]+1 <<" 0 " <<cart[3]+1 <<"\n";
											//std::cout<<"\t" <<dir_fid_sub[cart[5]][1][cart[1]+1][0][cart[3]+1]->GetName() <<"\n";
											//std::cout<<"writing hist: " <<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" " <<_weight_[cart[6]] <<" " <<_no_p_dep_ <<"\n";
											//fun::print_vector_idx(idx);
											dir_fid_sub[cart[5]][1][cart[1]+1][0][cart[3]+1]->cd();
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetXTitle("W (GeV)");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetYTitle("Q^{2} (GeV^{2}");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetOption("Colz");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Write();
										}else if(_species_[cart[5]]==_pim_){
											//std::cout<<"writing hist: " <<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" " <<_weight_[cart[6]] <<" " <<_no_p_dep_ <<"\n";
											//std::cout<<"Looking at Fid Dir: " <<cart[5] <<" 1 " <<cart[1]+1 <<" 1 " <<cart[3]+1 <<"\n";
											//std::cout<<"writing hist: " <<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" " <<_weight_[cart[6]] <<" " <<_no_p_dep_ <<"\n";
											//fun::print_vector_idx(idx);
											dir_fid_sub[cart[5]][1][cart[1]+1][1][cart[3]+1]->cd();
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetXTitle("W (GeV)");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetYTitle("Q^{2} (GeV^{2}");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetOption("Colz");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Write();
										}
									}else if(fun::hcut_perform(_species_[cart[5]],_hcuts_[cart[4]],flags_)){//PID Plots
										if(cart[0]==_p_bin_bins_){
											//std::cout<<"writing hist: " <<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" " <<_weight_[cart[6]] <<" " <<_no_p_dep_ <<"\n";
											//std::cout<<"Looking at Fid Dir: " <<cart[5] <<" 0 " <<cart[4]+1 <<" 0 " <<cart[3]+1 <<"\n";
											//std::cout<<"\t" <<dir_fid_sub[cart[5]][0][cart[4]+1][0][cart[3]+1]->GetName() <<"\n";
											
											//fun::print_vector_idx(idx);
											dir_fid_sub[cart[5]][0][cart[4]+1][0][cart[3]+1]->cd();
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetXTitle("W (GeV)");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetYTitle("Q^{2} (GeV^{2}");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetOption("Colz");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Write();
										}else if(_species_[cart[5]]==_pim_){
											//std::cout<<"writing hist: " <<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" " <<_weight_[cart[6]] <<" " <<_no_p_dep_ <<"\n";
											//std::cout<<"Looking at Fid Dir: " <<cart[5] <<" 0 " <<cart[4]+1 <<" 1 " <<cart[3]+1 <<"\n";
											//std::cout<<"writing hist: " <<_species_[cart[5]] <<" " <<_hcuts_[cart[4]] <<" " <<_sector_[cart[3]] <<" " <<_cut_[cart[2]] <<" " <<_top_[cart[1]] <<" " <<_weight_[cart[6]] <<" " <<_no_p_dep_ <<"\n";
											//fun::print_vector_idx(idx);
											dir_fid_sub[cart[5]][0][cart[4]+1][1][cart[3]+1]->cd();
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetXTitle("W (GeV)");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetYTitle("Q^{2} (GeV^{2}");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->SetOption("Colz");
											_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]][idx[6]]->Write();
										}
									}
								}
							}
						}
					}
					idx.clear();
				}
			}
		}
		std::cout<<"\tDone Writing Histograms\n";
	}
}
 //*--------------------------------End Fid Plot----------------------------*

 //*-------------------------------Start SF Plot----------------------------*

void Histogram::SF_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_SF()){
		std::cout<<"Making SF Histograms\n";
		Bool_1d made_1d;
		Bool_2d made_2d;
		Bool_3d made_3d;
		Bool_4d made_4d;
		TH2F_ptr_1d plot_1d;
		TH2F_ptr_2d plot_2d;
		TH2F_ptr_3d plot_3d;
		TH2F_ptr_4d plot_4d;
		//std::cout<<"Start Combo Stuff\n";
		std::vector<long> space_dims(5);
		space_dims[4] = std::distance(std::begin(_ecuts_),std::end(_ecuts_));//Ecuts
		space_dims[3] = std::distance(std::begin(_cut_),std::end(_cut_));//Cut applied?
		space_dims[2] = std::distance(std::begin(_top_),std::end(_top_));//Topology
		space_dims[1] = std::distance(std::begin(_sector_),std::end(_sector_));//sector
		space_dims[0] = Histogram::W_bins()+2;//W Dependence + exp range + full range
		char hname[100]; 
		CartesianGenerator cart2(space_dims);//For Made Histograms
		CartesianGenerator cart(space_dims);//For Histograms
		//std::cout<<"Making 'Made' Array\n";
		//Make "Made" Vectors
		/*while(cart2.GetNextCombination()){
			made_1d.push_back(false);
			if(cart2[0] == space_dims[0]-1 && made_1d.size()>0){
				made_2d.push_back(made_1d);
				//std::cout<<"size of 1d: " <<made_1d.size() <<"\t size of 2d: " <<made_2d.size() <<"\n";
				made_1d.clear();
			}
			if(cart2[1] == space_dims[1]-1 && made_2d.size()>0 && cart2[0] == space_dims[0]-1){
				made_3d.push_back(made_2d);
				//std::cout<<"size of 2d: " <<made_2d.size() <<"\t size of 3d: " <<made_3d.size()<<"\n";
				made_2d.clear();
			}
			if(cart2[2] == space_dims[2]-1 && made_3d.size()>0 && cart2[1] == space_dims[1]-1 && cart2[0] == space_dims[0]-1){
				made_4d.push_back(made_3d);
				//std::cout<<"size of 3d: " <<made_3d.size() <<"\t size of 4d: " <<made_4d.size()<<"\n";
				made_3d.clear();
			}
			if(cart2[3] == space_dims[3]-1 && made_4d.size()>0 && cart2[1] == space_dims[1]-1 && cart2[0] == space_dims[0]-1 && cart2[2] == space_dims[2]-1){
				_SF_made.push_back(made_4d);
				//std::cout<<"size of 4d: " <<made_4d.size() <<"\t size of 5d: " <<_SF_made.size()<<"\n";
				made_4d.clear();
			}
		}*/
		//std::cout<<"Making Actual Histograms\n";
		//Make Actual Histograms
		while(cart.GetNextCombination()){
			if(cart[0]>Histogram::W_bins()-1){
				//std::cout<<"SF trying to make: " <<_ecuts_[cart[4]] <<" " <<_cut_[cart[3]] <<" " <<_top_[cart[2]] <<" " <<_sector_[cart[1]] <<" " <<cart[0] <<"\n";
				//std::cout<<"\ttop perform? " <<fun::top_perform(_top_[cart[2]],flags_) <<"\n";
			}
			if(_ecuts_[cart[4]]==_event_ && _top_[cart[2]]!=_mnone_ && fun::top_perform(_top_[cart[2]],flags_)){//Event Selection
				if(_cut_[cart[3]]!=_no_cut_){
					//if(cart[0]>Histogram::W_bins()-1){
						//std::cout<<"Trying to make an event histogram: " <<cart[0] <<"\n";
					//}
					if(cart[0] < Histogram::W_bins()){//W Dependence
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"SF_%s_%s_%s_%s_W:%f-%f",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],Histogram::W_bot(cart[0]),Histogram::W_top(cart[0]));
						//std::cout<<"\t" <<hname <<"\n";
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						//_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}else if(cart[0] == Histogram::W_bins()){//Experimental W range
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"SF_%s_%s_%s_%s_W:%s",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],"in_range");
						//std::cout<<"\t" <<hname <<"\n";
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						//_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}else{//Full W range
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"SF_%s_%s_%s_%s_W:%s",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],"all");
						//std::cout<<"\t" <<hname <<"\n";
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						//_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}
				}
			}else if(_ecuts_[cart[4]]!=_event_ && _top_[cart[2]]==_mnone_ && ((_cut_[cart[3]]!=_no_cut_ && _ecuts_[cart[4]]!=_none_) || ((_cut_[cart[3]]==_no_cut_ && _ecuts_[cart[4]]==_none_)))){//Electron ID Cuts
				if(fun::ecut_perform(_ecuts_[cart[4]],flags_)){
					if(cart[0]>Histogram::W_bins()-1){
						//std::cout<<"Trying to make a PID histogram: " <<cart[0] <<"\n";
					}
					if(cart[0] < Histogram::W_bins()){//W Dependence
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"SF_%s_%s_%s_%s_W:%f-%f",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],Histogram::W_bot(cart[0]),Histogram::W_top(cart[0]));
						//std::cout<<"\t" <<hname <<"\n";
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						//_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}else if(cart[0] == Histogram::W_bins()){//Experimental W range
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"SF_%s_%s_%s_%s_W:%s",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],"in_range");
						//std::cout<<"\t" <<hname <<"\n";
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						//_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}else{//Full W range
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n"; 
						sprintf(hname,"SF_%s_%s_%s_%s_W:%s",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],"all");
						//std::cout<<"\t" <<hname <<"\n";
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						//_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}
				}
			}
			if(cart[0]==space_dims[0]-1 && plot_1d.size()>0){
				plot_2d.push_back(plot_1d);
				//std::cout<<"size of 1d: " <<plot_1d.size() <<"\t size of 2d: " <<plot_2d.size()<<"\n";
				plot_1d.clear();
			}
			if(cart[1]==space_dims[1]-1 && plot_2d.size()>0 && cart[0]==space_dims[0]-1){
				plot_3d.push_back(plot_2d);
				//std::cout<<"size of 2d: " <<plot_2d.size() <<"\t size of 3d: " <<plot_3d.size()<<"\n";
				plot_2d.clear();
			}
			if(cart[2]==space_dims[2]-1 && plot_3d.size()>0 && cart[0]==space_dims[0]-1 && cart[1]==space_dims[1]-1){
				plot_4d.push_back(plot_3d);
				//std::cout<<"size of 3d: " <<plot_3d.size() <<"\t size of 4d: " <<plot_4d.size()<<"\n";
				plot_3d.clear();
			}
			if(cart[3]==space_dims[3]-1 && plot_4d.size()>0 && cart[0]==space_dims[0]-1 && cart[1]==space_dims[1]-1 && cart[2]==space_dims[2]-1){
				_SF_hist.push_back(plot_4d);
				//std::cout<<"size of 4d: " <<plot_4d.size() <<"\t size of 5d: " <<_SF_hist.size()<<"\n";
				plot_4d.clear();
			}
		}
	}
}

std::vector<int> Histogram::SF_idx(const char* ecut_, const char* cut_, const char* top_, const char* sector_, const char * W_dep_, std::shared_ptr<Flags> flags_, float W_){
	std::vector<int> idx; 
	if(flags_->Flags::Plot_SF()){
		if(fun::ecut_perform(ecut_,flags_)){
			if(cut_==_no_cut_ && ecut_ == _none_ && top_==_mnone_){
				idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
			}else if(cut_ != _no_cut_ && ecut_!=_none_){
				idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
			}else{
				idx.push_back(-1);
			}
		}else{
			idx.push_back(-1);
		}
		if(cut_==_no_cut_){
			idx.push_back(0);
		}else{
			idx.push_back(fun::cut_idx(cut_));
		}
		if(ecut_==_event_ && top_ != _mnone_ ){
			idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
		}else if(ecut_!=_event_ && top_==_mnone_){
			//idx.push_back(fun::top_idx(top_));
			idx.push_back(0);
		}else{
			idx.push_back(-1);
		}
		idx.push_back(fun::sector_idx(sector_));
		if(W_dep_==_W_var_){
			idx.push_back(Histogram::W_bin(W_));
		}else if(W_dep_ == _W_range_){
			idx.push_back(Histogram::W_bins());
		}else if(W_dep_ == _W_all_){
			idx.push_back(Histogram::W_bins()+1);
		}else{
			idx.push_back(-1);
		}
	}else{
		for(int i=0; i<5; i++){
			idx.push_back(-1);
		}
	}
	//std::cout<<"SF idx: ";
	//std::cout<<"SF input: " <<ecut_ <<" " <<cut_ <<" " <<top_ <<" " <<sector_ <<" " <<W_dep_ <<"\n";
	//fun::print_vector_idx(idx);
	return idx; 
}

void Histogram::SF_Fill(float p_, float sf_, float W_, const char* ecut_, const char* cut_, const char* top_, const char* sector_, float weight_, std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_SF()){
		std::vector<int> idx = Histogram::SF_idx(ecut_,cut_,top_,sector_,_W_var_,flags_,W_);
		//if(ecut_==_event_){
		//	std::cout<<"Filling SF for an event: " <<ecut_ <<" " <<cut_ <<" " <<top_ <<" " <<sector_ <<" p:" <<p_ <<" sf:" <<sf_ <<" W:" <<W_ <<" idx:";
		//	fun::print_vector_idx(idx);
		//}
		//std::cout<<"Filling SF| W:" <<W_ <<"\t" <<ecut_ <<" " <<cut_ <<" " <<top_ <<" " <<sector_ <<" " <<" idx:";
		
		std::vector<int> idx2 = Histogram::SF_idx(ecut_,cut_,top_,_sec_all_,_W_var_,flags_,W_);
		if(Histogram::W_bin(W_) >=0){
			
			//std::cout<<"Looking at _SF_made idx: " <<fun::ecut_idx(ecut_) <<" " <<fun::cut_idx(cut_)<<" "<<fun::top_idx(top_)<<" "<<fun::sector_idx(sector_)<<" "<<Histogram::W_bin(W_)<<"\n";
			//if(_SF_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::sector_idx(sector_)][Histogram::W_bin(W_)]){
				if(Histogram::OK_Idx(idx)){
					//fun::print_vector_idx(idx);
					_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,sf_,weight_);
				}
				if(Histogram::OK_Idx(idx2)){
					//fun::print_vector_idx(idx2);
					_SF_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]][idx2[4]]->Fill(p_,sf_,weight_);
				}
			//}
		}
		idx.clear();
		idx2.clear();
		//std::cout<<"Looking at _SF_made idx: " <<fun::ecut_idx(ecut_)<<" " <<fun::cut_idx(cut_)<<" "<<fun::top_idx(top_)<<" "<<fun::sector_idx(sector_)<<" "<<Histogram::W_bins()<<"\n";
		//if(_SF_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::sector_idx(sector_)][Histogram::W_bins()] && 
		if(cuts::in_range(W_)){
			idx = Histogram::SF_idx(ecut_,cut_,top_,sector_,_W_range_,flags_,W_);
			//fun::print_vector_idx(idx);
			idx2 = Histogram::SF_idx(ecut_,cut_,top_,_sec_all_,_W_range_,flags_,W_);
			if(Histogram::OK_Idx(idx)){
				//fun::print_vector_idx(idx);
				_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,sf_,weight_);
			}
			if(Histogram::OK_Idx(idx2)){
				//fun::print_vector_idx(idx2);
				_SF_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]][idx2[4]]->Fill(p_,sf_,weight_);
			}
			idx.clear();
			idx2.clear();
		}
		//std::cout<<"Looking at _SF_made idx: " <<fun::ecut_idx(ecut_)<<" " <<fun::cut_idx(cut_)<<" "<<fun::top_idx(top_)<<" "<<fun::sector_idx(sector_)<<" "<<Histogram::W_bins()+1<<"\n";
		//if(_SF_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::sector_idx(sector_)][Histogram::W_bins()+1] && 
		if(W_>0){
			idx = Histogram::SF_idx(ecut_,cut_,top_,sector_,_W_all_,flags_,W_);
			//fun::print_vector_idx(idx);
			idx2 = Histogram::SF_idx(ecut_,cut_,top_,_sec_all_,_W_all_,flags_,W_);
			if(Histogram::OK_Idx(idx)){
				//fun::print_vector_idx(idx);
				_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,sf_,weight_);
			}
			if(Histogram::OK_Idx(idx2)){
				//fun::print_vector_idx(idx2);
				_SF_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]][idx2[4]]->Fill(p_,sf_,weight_);
			}
			idx.clear();
			idx2.clear();
		}
	}
}

std::vector<int> Histogram::SF_dir_idx(const char* ecut_, const char* cut_, const char* top_, const char* sector_, const char* W_dep_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx; 
	if((ecut_ != _event_ && top_ == _mnone_) || (ecut_ == _event_ && top_ != _mnone_)){
		idx.push_back(fun::ecut_idx(ecut_) + fun::ecut_offset(ecut_,flags_));
	}else{
		idx.push_back(-1);
	}
	if((cut_==_no_cut_ && ecut_==_none_) || (cut_!=_no_cut_ && ecut_!=_none_)){
		idx.push_back(fun::cut_idx(cut_)+1);
	}else{
		idx.push_back(-1);
	}
	if(fun::top_perform(top_,flags_)){
		idx.push_back(fun::top_idx(top_)+1);
	}else{
		idx.push_back(-1);
	}
	idx.push_back(fun::sector_idx(sector_));
	if(W_dep_==_W_var_){
		idx.push_back(1);
	}else{
		idx.push_back(2);
	}
	return idx;
}

void Histogram::SF_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_SF()){
		std::cout<<"Writing SF Plots\n";
		TDirectory* dir_sf = _RootOutputFile->mkdir("SF");
		dir_sf->cd();
		TDirectory* dir_sf_sub[std::distance(std::begin(_ecuts_),std::end(_ecuts_))][std::distance(std::begin(_cut_),std::end(_cut_))+1][std::distance(std::begin(_top_),std::end(_top_))+1][std::distance(std::begin(_sector_),std::end(_sector_))+1][2+1];
		char dir_name[100];
		for(int i=0; i<std::distance(std::begin(_ecuts_),std::end(_ecuts_)); i++){//ECUT Index
			if(fun::ecut_perform(_ecuts_[i],flags_)){
				sprintf(dir_name,"sf_%s",_ecuts_[i]);
				dir_sf_sub[i][0][0][0][0] = dir_sf->mkdir(dir_name);
				for(int j=0; j<std::distance(std::begin(_cut_),std::end(_cut_)); j++){
					if(_cut_[j]==_no_cut_ && _ecuts_[i]==_none_){
						sprintf(dir_name,"sf_%s_%s",_ecuts_[i],_cut_[j]);
						dir_sf_sub[i][j+1][0][0][0] = dir_sf_sub[i][0][0][0][0]->mkdir(dir_name);
						//for(int k=0; k<std::distance(std::begin(_top_),std::end(_top_)); k++){
							//if(fun::top_perform(_top_[k],flags_)){
						sprintf(dir_name,"sf_%s_%s_%s",_ecuts_[i],_cut_[j],_mnone_);
						dir_sf_sub[i][j+1][fun::top_idx(_mnone_)+1][0][0] = dir_sf_sub[i][j+1][0][0][0]->mkdir(dir_name);
						for(int l=0; l<std::distance(std::begin(_sector_),std::end(_sector_)); l++){
							sprintf(dir_name,"sf_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_mnone_,_sector_[l]);
							dir_sf_sub[i][j+1][fun::top_idx(_mnone_)+1][l+1][0] = dir_sf_sub[i][j+1][fun::top_idx(_mnone_)+1][0][0]->mkdir(dir_name);
							sprintf(dir_name,"sf_%s_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_mnone_,_sector_[l],_W_var_);
							dir_sf_sub[i][j+1][fun::top_idx(_mnone_)+1][l+1][1] = dir_sf_sub[i][j+1][fun::top_idx(_mnone_)+1][l+1][0]->mkdir(dir_name);
							sprintf(dir_name,"sf_%s_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_mnone_,_sector_[l],"W_range");
							dir_sf_sub[i][j+1][fun::top_idx(_mnone_)+1][l+1][2] = dir_sf_sub[i][j+1][fun::top_idx(_mnone_)+1][l+1][0]->mkdir(dir_name);
						}
					}else if(_cut_[j]!=_no_cut_ && _ecuts_[i]!=_none_){
						sprintf(dir_name,"sf_%s_%s",_ecuts_[i],_cut_[j]);
						dir_sf_sub[i][j+1][0][0][0] = dir_sf_sub[i][0][0][0][0]->mkdir(dir_name);
						for(int k=0; k<std::distance(std::begin(_top_),std::end(_top_)); k++){
							if((_ecuts_[i]==_event_ && _top_[k]!=_mnone_) || (_ecuts_[i]!=_event_ && _top_[k]==_mnone_)){
								if(fun::top_perform(_top_[k],flags_)){
									sprintf(dir_name,"sf_%s_%s_%s",_ecuts_[i],_cut_[j],_top_[k]);
									dir_sf_sub[i][j+1][k+1][0][0] = dir_sf_sub[i][j+1][0][0][0]->mkdir(dir_name);
									for(int l=0; l<std::distance(std::begin(_sector_),std::end(_sector_)); l++){
										sprintf(dir_name,"sf_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_top_[k],_sector_[l]);
										dir_sf_sub[i][j+1][k+1][l+1][0] = dir_sf_sub[i][j+1][k+1][0][0]->mkdir(dir_name);
										sprintf(dir_name,"sf_%s_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_top_[k],_sector_[l],_W_var_);
										dir_sf_sub[i][j+1][k+1][l+1][1] = dir_sf_sub[i][j+1][k+1][l+1][0]->mkdir(dir_name);
										sprintf(dir_name,"sf_%s_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_top_[k],_sector_[l],"W_range");
										dir_sf_sub[i][j+1][k+1][l+1][2] = dir_sf_sub[i][j+1][k+1][l+1][0]->mkdir(dir_name);
									}
								}
							}
						}
					}
				}
			}
		}
		std::vector<long> space_dims(5);
		space_dims[4] = std::distance(std::begin(_ecuts_),std::end(_ecuts_));//Ecuts
		space_dims[3] = std::distance(std::begin(_cut_),std::end(_cut_));//Cut applied?
		space_dims[2] = std::distance(std::begin(_top_),std::end(_top_));//Topology
		space_dims[1] = std::distance(std::begin(_sector_),std::end(_sector_));//sector
		space_dims[0] = Histogram::W_bins() + 2;//W Dependence + exp range + full range
		CartesianGenerator cart(space_dims);
		std::vector<int> idx;
		while(cart.GetNextCombination()){
			//if(_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] && 
			if(fun::ecut_perform(_ecuts_[cart[4]],flags_) && fun::top_perform(_top_[cart[2]],flags_)){
				if(cart[0]<Histogram::W_bins()){
					//std::cout<<"Will be trying to write a W varying SF histogram\n";
					//std::cout<<"\t" <<_ecuts_[cart[4]] <<" " <<_cut_[cart[3]] <<" " <<_top_[cart[2]] <<" " <<_sector_[cart[1]] <<" " <<_W_var_ <<"\n";
					idx = Histogram::SF_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],_W_var_,flags_,Histogram::W_center(cart[0]));
					//fun::print_vector_idx(idx);
				}else if(cart[0]==Histogram::W_bins()){
					idx = Histogram::SF_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],_W_range_,flags_);
				}else if(cart[0]==Histogram::W_bins()+1){
					idx = Histogram::SF_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],_W_all_,flags_);
				}
				//idx[4] = (Histogram::W_bins()+1)-((Histogram::W_bins()+1)-cart[0]);
				//fun::print_vector_idx(idx);
				if(Histogram::OK_Idx(idx)){
					if(cart[0]<Histogram::W_bins()){
						//std::cout<<"\tTrying to write a W varying SF histogram\n";
						dir_sf_sub[cart[4]][cart[3]+1][cart[2]+1][cart[1]+1][1]->cd();
						_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
						_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Sampling Fraction");
						_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
					}else{
						dir_sf_sub[cart[4]][cart[3]+1][cart[2]+1][cart[1]+1][2]->cd();
						_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
						_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Sampling Fraction");
						_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
					}
				}
				idx.clear();
			}
		}
	}
}
 //*-------------------------------End SF Plot------------------------------*

 //*-------------------------------Start CC Plot----------------------------*
void Histogram::CC_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_CC()){
		std::cout<<"Making CC Histograms\n";
		Bool_5d made_5d;
		Bool_4d made_4d;
		Bool_3d made_3d;
		Bool_2d made_2d;
		Bool_1d made_1d;
		TH1F_ptr_5d plot_5d;
		TH1F_ptr_4d plot_4d;
		TH1F_ptr_3d plot_3d;
		TH1F_ptr_2d plot_2d;
		TH1F_ptr_1d plot_1d;
		std::vector<long> space_dims(6); //{ecut,top,sector,side,segment}
		space_dims[5] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		space_dims[4] = std::distance(std::begin(_cut_),std::end(_cut_));
		space_dims[3] = std::distance(std::begin(_top_), std::end(_top_));
		space_dims[2] = std::distance(std::begin(_sector_), std::end(_sector_));
		space_dims[1] = std::distance(std::begin(_cc_sides_), std::end(_cc_sides_));
		space_dims[0] = 18;

		CartesianGenerator cart2(space_dims);
		CartesianGenerator cart(space_dims);
		char hname[100];

		while(cart2.GetNextCombination()){
			made_1d.push_back(false);
			if(cart2[0]== space_dims[0]-1 && made_1d.size()>0){
				made_2d.push_back(made_1d);
				//std::cout<<"\tLength of made_1d: " <<made_1d.size()<<"\tLength of made_2d: " <<made_2d.size() <<"\n";
				made_1d.clear();
			}
			if(cart2[1] == space_dims[1]-1 && made_2d.size()>0 && cart2[0]== space_dims[0]-1){
					made_3d.push_back(made_2d);
					//std::cout<<"\tLength of made_2d: " <<made_2d.size()<<"\tLength of made_3d: " <<made_3d.size() <<"\n";
					made_2d.clear();
			}
			if(cart2[2] == space_dims[2]-1 && made_3d.size()>0 && cart2[1] == space_dims[1]-1 && cart2[0]== space_dims[0]-1){
					made_4d.push_back(made_3d);
					//std::cout<<"\tLength of made_1d: " <<made_1d.size()<<"\tLength of made_2d: " <<made_2d.size() <<"\n";
					made_3d.clear();
			}
			if(cart2[3] == space_dims[3]-1 && made_4d.size()>0 && cart2[1] == space_dims[1]-1 && cart2[0]== space_dims[0]-1 && cart2[2] == space_dims[2]-1){
					made_5d.push_back(made_4d);
					made_4d.clear();
			}
			if(cart2[4] == space_dims[4]-1 && made_5d.size()>0 && cart2[1] == space_dims[1]-1 && cart2[0]== space_dims[0]-1 && cart2[2] == space_dims[2]-1 && cart2[3] == space_dims[3]-1){
					_CC_made.push_back(made_5d);
					made_5d.clear();
			}

		}
		int ecut_idx = -1;
		int cut_idx = -1;
		int top_idx = -1;
		int sec_idx = -1;
		int side_idx = -1;
		int seg_idx = -1;
		while(cart.GetNextCombination()){
			ecut_idx = cart[5];
			cut_idx = cart[4];
			top_idx = cart[3];
			sec_idx = cart[2];
			side_idx = cart[1];
			seg_idx = cart[0];
			if(fun::ecut_perform(_ecuts_[ecut_idx],flags_) && fun::top_perform(_top_[top_idx],flags_)){
				if(_ecuts_[ecut_idx] != _event_ && _top_[top_idx] == _mnone_){
					if(_ecuts_[ecut_idx] == _none_ && _cut_[cut_idx]==_no_cut_){
						//std::cout<<"CC hist Made @: " <<_CC_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"cc_%s_%s_%s_%s_%s_seg%d",_ecuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_cc_sides_[side_idx],seg_idx+1);
						plot_1d.push_back(new TH1F(hname,hname,_cc_xbin_,_cc_xmin_,_cc_xmax_));
						//_CC_made[ecut_idx][cut_idx][top_idx][sec_idx][side_idx][seg_idx]=true;
					}else if(_ecuts_[ecut_idx] != _none_ && _cut_[cut_idx]!=_no_cut_){
						//std::cout<<"CC hist Made @: " <<_CC_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"cc_%s_%s_%s_%s_%s_seg%d",_ecuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_cc_sides_[side_idx],seg_idx+1);
						plot_1d.push_back(new TH1F(hname,hname,_cc_xbin_,_cc_xmin_,_cc_xmax_));
						//_CC_made[ecut_idx][cut_idx][top_idx][sec_idx][side_idx][seg_idx]=true;
					}
				}else if(_ecuts_[ecut_idx] == _event_ && _top_[top_idx] != _mnone_ && _cut_[cut_idx]!=_no_cut_){
					//std::cout<<"CC hist Made @: " <<_CC_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
					sprintf(hname,"cc_%s_%s_%s_%s_%s_seg%d",_ecuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_cc_sides_[side_idx],seg_idx+1);
					plot_1d.push_back(new TH1F(hname,hname,_cc_xbin_,_cc_xmin_,_cc_xmax_));
					//_CC_made[ecut_idx][cut_idx][top_idx][sec_idx][side_idx][seg_idx]=true;
				}
			}
			if(cart[0] == space_dims[0]-1 ){
				if(plot_1d.size()>0){
					plot_2d.push_back(plot_1d);
					plot_1d.clear();
				}
				if(cart[1] == space_dims[1]-1){
					if(plot_2d.size()>0){
						plot_3d.push_back(plot_2d);
						plot_2d.clear();
					}
					if(cart[2] == space_dims[2]-1){
						if(plot_3d.size()>0){
							plot_4d.push_back(plot_3d);
							plot_3d.clear();
						}
						if(cart[3] == space_dims[3]-1){
							if(plot_4d.size()>0){
								plot_5d.push_back(plot_4d);
								plot_4d.clear();
							}
							if(cart[4] == space_dims[4]-1){
								if(plot_5d.size()>0){
									_CC_hist.push_back(plot_5d);
									plot_5d.clear();
								}
							}
						}
					}
				}
			}
			ecut_idx = -1;
			cut_idx = -1;
			top_idx = -1;
			sec_idx = -1;
			side_idx = -1;
			seg_idx = -1;
		}
	}
}

std::vector<int> Histogram::CC_idx(const char * ecut_, const char * cut_, const char* top_, const char * sector_, const char* side_, int seg_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	//std::cout<<"\tLooking at CC idx: "  <<ecut_ <<" " <<cut_ <<" " <<top_ <<" " <<sector_ <<" " <<side_ <<" seg:" <<seg_;
	if(ecut_==_none_ && cut_ == _no_cut_){
		idx.push_back(0);
		idx.push_back(0);
	}else if(ecut_!=_none_ && cut_ != _no_cut_){
		if(fun::ecut_perform(ecut_,flags_)){
			idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
		}else{
			idx.push_back(-1);
		}
		idx.push_back(fun::cut_idx(cut_));
	}else{
		idx.push_back(-1);
		idx.push_back(-1);
	}
	if(ecut_==_event_ && top_!=_mnone_){
		if(fun::top_perform(top_,flags_)){
			idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
		}else{
			idx.push_back(-1);
		}
	}else if(ecut_!=_event_ && top_==_mnone_){
		idx.push_back(0);
	}else{
		idx.push_back(-1);
	}
	idx.push_back(fun::sector_idx(sector_));
	idx.push_back(fun::cc_side_idx(side_));
	idx.push_back(seg_);
	//fun::print_vector_idx(idx);
	return idx;
}

void Histogram::CC_Fill(int nphe_, int seg_, const char * ecut_, const char* cut_, const char* top_, const char * sector_, const char* side_, std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_CC() && !flags_->Flags::Sim()){
		std::vector<int> idx = Histogram::CC_idx(ecut_,cut_,top_,sector_,side_,seg_,flags_);
		//std::cout<<"Filling CC: " <<nphe_ <<"\t" <<ecut_ <<" " <<cut_ <<" " <<top_ <<" " <<sector_ <<" " <<side_ <<" seg:" <<seg_;
		//fun::print_vector_idx(idx);
		//std::cout<<"\n";
		if(Histogram::OK_Idx(idx)){
			//if(_CC_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::sector_idx(sector_)][fun::cc_side_idx(side_)][seg_-1]){
				_CC_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(nphe_);
		//	}
		}
		idx.clear();
		idx = Histogram::CC_idx(ecut_,cut_,top_,_sec_all_,side_,seg_,flags_);
		if(Histogram::OK_Idx(idx)){
			//if(_CC_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::sector_idx(_sec_all_)][fun::cc_side_idx(side_)][seg_-1]){
				_CC_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(nphe_);
			//}
		}
		idx.clear();
	}
}

void Histogram::CC_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_CC()){
		std::cout<<"Writing CC Plots\n";
		int ecut_len = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		int cut_len = std::distance(std::begin(_cut_),std::end(_cut_));
		int top_len = std::distance(std::begin(_top_), std::end(_top_));
		int sec_len = std::distance(std::begin(_sector_), std::end(_sector_));
		int side_len = std::distance(std::begin(_cc_sides_), std::end(_cc_sides_));
		TDirectory* dir_cc = _RootOutputFile->mkdir("CC");
		dir_cc->cd();
		TDirectory* dir_cc_sub[ecut_len][cut_len+1][top_len+1][sec_len+1][side_len+1];
		char dir_name[100];
		for(int i=0; i<ecut_len; i++){

			if(fun::ecut_perform(_ecuts_[i],flags_)){
				if(_ecuts_[i] == _none_){
					//std::cout<<"\t" <<i <<" 0 0 0 0" <<"\n";
					sprintf(dir_name,"min_cc_%s",_ecuts_[i]);
					dir_cc_sub[i][0][0][0][0] = dir_cc->mkdir(dir_name);
					for(int l=0; l<sec_len; l++){//Sector
							//std::cout<<"\t" <<i <<" " <<fun::cut_idx(_no_cut_)+1 <<" " <<fun::top_idx(_mnone_)+1 <<" " <<l+1 <<" " <<0 <<"\n";
							sprintf(dir_name,"min_cc_%s_%s_%s_%s",_ecuts_[i],_no_cut_,_mnone_,_sector_[l]);
							dir_cc_sub[i][fun::cut_idx(_no_cut_)+1][fun::top_idx(_mnone_)+1][l+1][0] = dir_cc_sub[i][0][0][0][0]->mkdir(dir_name);
							for(int m=0; m<side_len; m++){//Side
								//std::cout<<"\t" <<i <<" " <<fun::cut_idx(_no_cut_)+1 <<" " <<fun::top_idx(_mnone_)+1 <<" " <<l+1 <<" " <<m+1 <<"\n";
								sprintf(dir_name,"min_cc_%s_%s_%s_%s_%s",_ecuts_[i],_no_cut_,_mnone_,_sector_[l],_cc_sides_[m]);
								dir_cc_sub[i][fun::cut_idx(_no_cut_)+1][fun::top_idx(_mnone_)+1][l+1][m+1] = dir_cc_sub[i][fun::cut_idx(_no_cut_)+1][fun::top_idx(_mnone_)+1][l+1][0]->mkdir(dir_name);
							}
						}
				}else if(_ecuts_[i] == _event_){
					sprintf(dir_name,"min_cc_%s",_ecuts_[i]);
					//std::cout<<"\t" <<i <<" " <<0 <<" " <<0 <<" " <<0 <<" " <<0 <<"\n";
					dir_cc_sub[i][0][0][0][0] = dir_cc->mkdir(dir_name);
					for(int j=0; j<cut_len; j++){//Cut
						if(_cut_[j]!=_no_cut_){
							//std::cout<<"\t" <<i <<" " <<j+1 <<" " <<0 <<" " <<0 <<" " <<0 <<"\n";
							sprintf(dir_name,"min_cc_%s_%s",_ecuts_[i],_cut_[j]);
							dir_cc_sub[i][j+1][0][0][0] = dir_cc_sub[i][0][0][0][0]->mkdir(dir_name);
							for(int k=0; k<top_len; k++){//Top
								if(_top_[k]!=_mnone_ && fun::top_perform(_top_[k],flags_)){	
									//std::cout<<"\t" <<i <<" " <<j+1 <<" " <<k+1 <<" " <<0 <<" " <<0 <<"\n";
									sprintf(dir_name,"min_cc_%s_%s_%s",_ecuts_[i],_cut_[j],_top_[k]);
									dir_cc_sub[i][j+1][k+1][0][0] = dir_cc_sub[i][j+1][0][0][0]->mkdir(dir_name);
									for(int l=0; l<sec_len; l++){//Sector
										//std::cout<<"\t" <<i <<" " <<j+1 <<" " <<k+1 <<" " <<l+1 <<" " <<0 <<"\n";
										sprintf(dir_name,"min_cc_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_top_[k],_sector_[l]);
										dir_cc_sub[i][j+1][k+1][l+1][0] = dir_cc_sub[i][j+1][k+1][0][0]->mkdir(dir_name);
										for(int m=0; m<side_len; m++){//Side
											//std::cout<<"\t" <<i <<" " <<j+1 <<" " <<k+1 <<" " <<l+1 <<" " <<m+1 <<"\n";
											sprintf(dir_name,"min_cc_%s_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_top_[k],_sector_[l],_cc_sides_[m]);
											dir_cc_sub[i][j+1][k+1][l+1][m+1] = dir_cc_sub[i][j+1][k+1][l+1][0]->mkdir(dir_name);
										}
									}
								}
							}
						}
					}
				}else{
					//std::cout<<"\t" <<i <<" 0 0 0 0" <<"\n";
					sprintf(dir_name,"min_cc_%s",_ecuts_[i]);
					dir_cc_sub[i][0][0][0][0] = dir_cc->mkdir(dir_name);
					for(int j=0; j<cut_len; j++){//Cut
						if(_cut_[j]!=_no_cut_){
						//	std::cout<<"\t" <<i <<" " <<j+1 <<" " <<0 <<" " <<0 <<" " <<0 <<"\n";
							sprintf(dir_name,"min_cc_%s_%s",_ecuts_[i],_cut_[j]);
							dir_cc_sub[i][j+1][0][0][0] = dir_cc_sub[i][0][0][0][0]->mkdir(dir_name);
							//std::cout<<"\t" <<i <<" " <<j+1 <<" " <<fun::top_idx(_mnone_)+1 <<" " <<0 <<" " <<0 <<"\n";
							sprintf(dir_name,"min_cc_%s_%s_%s",_ecuts_[i],_cut_[j],_mnone_);
							dir_cc_sub[i][j+1][fun::top_idx(_mnone_)+1][0][0] = dir_cc_sub[i][j+1][0][0][0]->mkdir(dir_name);
							for(int l=0; l<sec_len; l++){//Sector
								//std::cout<<"\t" <<i <<" " <<j+1 <<" " <<fun::top_idx(_mnone_)+1 <<" " <<l+1 <<" " <<0 <<"\n";
								sprintf(dir_name,"min_cc_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_mnone_,_sector_[l]);
								dir_cc_sub[i][j+1][fun::top_idx(_mnone_)+1][l+1][0] = dir_cc_sub[i][j+1][fun::top_idx(_mnone_)+1][0][0]->mkdir(dir_name);
								for(int m=0; m<side_len; m++){//Side
									//std::cout<<"\t" <<i <<" " <<j+1 <<" " <<fun::top_idx(_mnone_)+1 <<" " <<l+1 <<" " <<m+1 <<"\n";
									sprintf(dir_name,"min_cc_%s_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_mnone_,_sector_[l],_cc_sides_[m]);
									dir_cc_sub[i][j+1][fun::top_idx(_mnone_)+1][l+1][m+1] = dir_cc_sub[i][j+1][fun::top_idx(_mnone_)+1][l+1][0]->mkdir(dir_name);
								}
							}
						}
					}
				}
			}
		}
		std::vector<long> space_dims(6); //{ecut,top,sector,side,segment}
		space_dims[5] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		space_dims[4] = std::distance(std::begin(_cut_),std::end(_cut_));
		space_dims[3] = std::distance(std::begin(_top_), std::end(_top_));
		space_dims[2] = std::distance(std::begin(_sector_), std::end(_sector_));
		space_dims[1] = std::distance(std::begin(_cc_sides_), std::end(_cc_sides_));
		space_dims[0] = 18;
		CartesianGenerator cart(space_dims);
		std::vector<int> idx; 
		while(cart.GetNextCombination()){
			//std::cout<<"\tdirectory idx: " <<cart[5] <<" " <<cart[4] <<" " <<cart[3] <<" " <<cart[2] <<" " <<cart[1] <<" " <<cart[0];
			idx = Histogram::CC_idx(_ecuts_[cart[5]],_cut_[cart[4]],_top_[cart[3]],_sector_[cart[2]],_cc_sides_[cart[1]],cart[0],flags_);
			//fun::print_vector_idx(idx);
			if(Histogram::OK_Idx(idx)){
				//std::cout<<"\tdirectory idx: " <<cart[5] <<" " <<cart[4] <<" " <<cart[3] <<" " <<cart[2] <<" " <<cart[1] <<" " <<cart[0];
				//std::cout<<"\n\tHist idx:";
				//fun::print_vector_idx(idx);
				dir_cc_sub[cart[5]][cart[4]+1][cart[3]+1][cart[2]+1][cart[1]+1]->cd();
				_CC_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("nphe");
				_CC_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Yield");
				_CC_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
			}
		}
	}
}
 //*-------------------------------End CC Plot----------------------------*

 //*-------------------------------Start EC Plot----------------------------*
 //*-------------------------------End EC Plot------------------------------*

 //*-------------------------------Start Vertex Plot----------------------------*
void Histogram::Vertex_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_Vertex()){
		std::cout<<"Making Vertex Histograms\n";
		TH1F_ptr_1d plot_1d;
		TH1F_ptr_2d plot_2d;
		TH1F_ptr_3d plot_3d;

		std::vector<long> space_dims(4); //{ecuts,cut,top}
		space_dims[3] = std::distance(std::begin(_sector_), std::end(_sector_));
		space_dims[2] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		space_dims[1] = std::distance(std::begin(_cut_), std::end(_cut_));
		space_dims[0] = std::distance(std::begin(_top_), std::end(_top_));
		char hname[100];
		CartesianGenerator cart(space_dims);
		while(cart.GetNextCombination()){
			if(_ecuts_[cart[2]] == _none_ && _cut_[cart[1]] == _no_cut_ && _top_[cart[0]]==_mnone_){
				//std::cout<<"Making Vertex IDX: " <<_Vertex_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
				sprintf(hname,"vertex_%s_%s_%s",_none_,_no_cut_,_sector_[cart[3]]);
				plot_1d.push_back(new TH1F(hname,hname,_vertex_bin_[flags_->Flags::Run()],_vertex_min_[flags_->Flags::Run()],_vertex_max_[flags_->Flags::Run()]));
			}else if(_ecuts_[cart[2]] != _none_ && _cut_[cart[1]] != _no_cut_){
				if(_ecuts_[cart[2]]==_event_ && _top_[cart[0]]!=_mnone_ && fun::top_perform(_top_[cart[0]],flags_)){
					//std::cout<<"Making Vertex IDX: " <<_Vertex_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
					sprintf(hname,"vertex_%s_%s_%s_%s",_ecuts_[cart[2]],_cut_[cart[1]],_top_[cart[0]],_sector_[cart[3]]);
					plot_1d.push_back(new TH1F(hname,hname,_vertex_bin_[flags_->Flags::Run()],_vertex_min_[flags_->Flags::Run()],_vertex_max_[flags_->Flags::Run()]));
				}else if(_ecuts_[cart[2]]!=_event_ && _top_[cart[0]]==_mnone_ && fun::ecut_perform(_ecuts_[cart[2]],flags_)){
					//std::cout<<"Making Vertex IDX: " <<_Vertex_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
					sprintf(hname,"vertex_%s_%s_%s_%s",_ecuts_[cart[2]],_cut_[cart[1]],_top_[cart[0]],_sector_[cart[3]]);
					plot_1d.push_back(new TH1F(hname,hname,_vertex_bin_[flags_->Flags::Run()],_vertex_min_[flags_->Flags::Run()],_vertex_max_[flags_->Flags::Run()]));
				}
			}
			if(cart[0] == space_dims[0]-1 ){
				if(plot_1d.size()>0){
					plot_2d.push_back(plot_1d);
					plot_1d.clear();
				}
				if(cart[1] == space_dims[1]-1){
					if(plot_2d.size()>0){
						plot_3d.push_back(plot_2d);
						plot_2d.clear();
					}
					if(cart[2] == space_dims[2]-1){
						if(plot_3d.size()>0){
							_Vertex_hist.push_back(plot_3d);
							plot_3d.clear();
						}
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::Vertex_idx(const char* ecut_, const char* cut_, const char* top_, const char* sector_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	idx.push_back(fun::sector_idx(sector_));
	if(ecut_ == _none_ && cut_ == _no_cut_){
		idx.push_back(fun::ecut_idx(ecut_) + fun::ecut_offset(ecut_,flags_));
		idx.push_back(0);
	}else if(ecut_ != _none_ && cut_ != _no_cut_ && fun::ecut_perform(ecut_,flags_)){
		idx.push_back(fun::ecut_idx(ecut_) + fun::ecut_offset(ecut_,flags_));
		idx.push_back(fun::cut_idx(cut_));
	}else{
		idx.push_back(-1);
		idx.push_back(-1);
	}
	if(ecut_ == _event_ && top_ != _mnone_ && fun::top_perform(top_,flags_)){
		idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
	}else if(ecut_ != _event_ && top_ == _mnone_){
		idx.push_back(0);
	}else{
		idx.push_back(-1);
	}
	//fun::print_vector_idx(idx);
	return idx;
}

void Histogram::Vertex_Fill(float vz_, float weight_, const char* ecut_, const char* cut_, const char* top_, const char* sector_, std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_Vertex()){
		std::vector<int> idx = Histogram::Vertex_idx(ecut_,cut_,top_,sector_,flags_);
		if(Histogram::OK_Idx(idx)){
			_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(vz_, weight_);
		}
		idx.clear();
		idx = Histogram::Vertex_idx(ecut_,cut_,top_,_sec_all_,flags_);
		if(Histogram::OK_Idx(idx)){
			_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(vz_, weight_);
		}
		idx.clear();
	}
}

void Histogram::Vertex_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_Vertex()){
		std::cout<<"Writing Vertex Histograms\n";
		char dir_name[100];
		TDirectory* dir_vert = _RootOutputFile->mkdir("Vertex");
		dir_vert->cd();
		TDirectory* dir_vert_sub[std::distance(std::begin(_ecuts_), std::end(_ecuts_))][std::distance(std::begin(_cut_), std::end(_cut_))+1][std::distance(std::begin(_top_), std::end(_top_))+1][std::distance(std::begin(_sector_),std::end(_sector_))+1];
		for(int i=0; i<std::distance(std::begin(_ecuts_), std::end(_ecuts_)); i++){
			if(fun::ecut_perform(_ecuts_[i],flags_)){
				sprintf(dir_name,"vertex_%s",_ecuts_[i]);
				dir_vert_sub[i][0][0][0] = dir_vert->mkdir(dir_name);
				if(_ecuts_[i] != _none_ ){
					for(int j=0; j<std::distance(std::begin(_cut_), std::end(_cut_)); j++){
						if(_cut_[j]!=_no_cut_){
							sprintf(dir_name,"vertex_%s_%s",_ecuts_[i],_cut_[j]);
							dir_vert_sub[i][j+1][0][0] = dir_vert_sub[i][0][0][0]->mkdir(dir_name);
							if(_ecuts_[i]==_event_){
								for(int k=0; k<std::distance(std::begin(_top_), std::end(_top_)); k++){
									if(_top_[k]!=_mnone_){
										sprintf(dir_name,"vertex_%s_%s_%s",_ecuts_[i],_cut_[j],_top_[k]);
										dir_vert_sub[i][j+1][k+1][0] = dir_vert_sub[i][j+1][0][0]->mkdir(dir_name);
										for(int l=0; l<std::distance(std::begin(_sector_),std::end(_sector_)); l++){
											sprintf(dir_name,"vertex_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_top_[k],_sector_[l]);
											dir_vert_sub[i][j+1][k+1][l+1] = dir_vert_sub[i][j+1][k+1][0]->mkdir(dir_name);
										}
									}
								}
							}else{
								for(int l=0; l<std::distance(std::begin(_sector_),std::end(_sector_)); l++){
									sprintf(dir_name,"vertex_%s_%s_%s_%s",_ecuts_[i],_cut_[j],_mnone_,_sector_[l]);
									dir_vert_sub[i][j+1][0][l+1] = dir_vert_sub[i][j+1][0][0]->mkdir(dir_name);
								}
							}
						}
					}
				}else{
					for(int l=0; l<std::distance(std::begin(_sector_),std::end(_sector_)); l++){
						sprintf(dir_name,"vertex_%s_%s_%s_%s",_ecuts_[i],_no_cut_,_mnone_,_sector_[l]);
						dir_vert_sub[i][0][0][l+1] = dir_vert_sub[i][0][0][0]->mkdir(dir_name);
					}
				}
			}
		}
		std::vector<long> space_dims(4); //{ecuts,cut,top}
		space_dims[3] = std::distance(std::begin(_sector_),std::end(_sector_));
		space_dims[2] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		space_dims[1] = std::distance(std::begin(_cut_), std::end(_cut_));
		space_dims[0] = std::distance(std::begin(_top_), std::end(_top_));
		CartesianGenerator cart(space_dims);
		std::vector<int> idx;
		while(cart.GetNextCombination()){
			if(_ecuts_[cart[2]]==_none_ && _cut_[cart[1]] == _no_cut_ && _top_[cart[0]]==_mnone_){
				idx = Histogram::Vertex_idx(_ecuts_[cart[2]],_cut_[cart[1]],_top_[cart[0]],_sector_[cart[3]],flags_);
				if(Histogram::OK_Idx(idx)){
					dir_vert_sub[cart[2]][0][0][cart[3]+1]->cd();
					_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetXTitle("Vz (nm)");
					_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetYTitle("Yield");
					_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Write();
				}
				idx.clear();
			}else if(_ecuts_[cart[2]]!=_none_ && _cut_[cart[1]] != _no_cut_){
				if(_ecuts_[cart[2]]==_event_ && _top_[cart[0]]!=_mnone_){
					if(fun::top_perform(_top_[cart[0]],flags_)){
						idx = Histogram::Vertex_idx(_ecuts_[cart[2]],_cut_[cart[1]],_top_[cart[0]],_sector_[cart[3]],flags_);
						if(Histogram::OK_Idx(idx)){
							dir_vert_sub[cart[2]][cart[1]+1][cart[0]+1][cart[3]+1]->cd();
							_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetXTitle("Vz (nm)");
							_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetYTitle("Yield");
							_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Write();
						}
						idx.clear();
					}
				}else if(_ecuts_[cart[2]]!=_event_ && _top_[cart[0]]==_mnone_){
					if(fun::ecut_perform(_ecuts_[cart[2]],flags_)){
						idx = Histogram::Vertex_idx(_ecuts_[cart[2]],_cut_[cart[1]],_top_[cart[0]],_sector_[cart[3]],flags_);
						if(Histogram::OK_Idx(idx)){
							dir_vert_sub[cart[2]][cart[1]+1][0][cart[3]+1]->cd();
							_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetXTitle("Vz (nm)");
							_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetYTitle("Yield");
							_Vertex_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Write();
						}
						idx.clear();
					}
				}
			}
		}
	}
}
 //*-------------------------------End Vertex Plot------------------------------*

 //*-------------------------------Start Delta T Plot----------------------------*
void Histogram::Delta_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_Delta(0) || flags_->Flags::Plot_Delta(1) || flags_->Flags::Plot_Delta(2) || flags_->Flags::Plot_Delta(3)){
		std::cout<<"Making Delta Histograms\n";
		TH2F_ptr_5d plot_5d;
		TH2F_ptr_4d plot_4d;
		TH2F_ptr_3d plot_3d;
		TH2F_ptr_2d plot_2d;
		TH2F_ptr_1d plot_1d;

		std::vector<long> space_dims(6);//{species,pcut,cut,top,w_dep}
		space_dims[5] = std::distance(std::begin(_sector_), std::end(_sector_));
		space_dims[4] = std::distance(std::begin(_species_), std::end(_species_));
		space_dims[3] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));
		space_dims[1] = std::distance(std::begin(_top_), std::end(_top_));
		space_dims[0] = std::distance(std::begin(_W_dep_), std::end(_W_dep_));
		char hname[100];
		int species_idx = -1;
		int pcut_idx = -1;
		int cut_idx = -1; 
		int top_idx = -1;
		int w_idx = -1; 
		int sec_idx = -1;
		CartesianGenerator cart(space_dims);
		while(cart.GetNextCombination()){
			species_idx = cart[4];
			pcut_idx = cart[3];
			cut_idx = cart[2];
			top_idx = cart[1];
			w_idx = cart[0];
			sec_idx = cart[5];
			//std::cout<<"\tsector: " <<_sector_[sec_idx] <<"\n";
			if(_species_[species_idx]==_ele_ && flags_->Plot_Delta(species_idx)){
				//std::cout<<_sector_[sec_idx] <<" " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
				if(_ecuts_[pcut_idx]==_none_ && _cut_[cut_idx]==_no_cut_ && _top_[top_idx]==_mnone_){
					if(_W_dep_[w_idx]==_W_var_){
						for(int i=0; i<Histogram::W_bins(); i++){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s_W:%.3f-%.3f",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],Histogram::W_bot(i),Histogram::W_top(i));
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}else if(_W_dep_[w_idx]==_W_range_){
						//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_range_);
						plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
					}
					else if(_W_dep_[w_idx]==_W_all_){
						//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_all_);
						plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
					}
				}else if(_ecuts_[pcut_idx]!=_none_ && _cut_[cut_idx]!=_no_cut_){
					if(_ecuts_[pcut_idx]==_event_ && _top_[top_idx]!=_mnone_ && fun::top_perform(_top_[top_idx],flags_)){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"delta_%s_%s_%s_%s_%s_W:%.3f-%.3f",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],Histogram::W_bot(i),Histogram::W_top(i));
								plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
							}
						}else if(_W_dep_[w_idx]==_W_range_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_range_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
						else if(_W_dep_[w_idx]==_W_all_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_all_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}else if(_ecuts_[pcut_idx]!=_event_ && _top_[top_idx]==_mnone_ && fun::ecut_perform(_ecuts_[pcut_idx],flags_)){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"delta_%s_%s_%s_%s_%s_W:%.3f-%.3f",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],Histogram::W_bot(i),Histogram::W_top(i));
								plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
							}
						}
						if(_W_dep_[w_idx]==_W_range_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_range_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}else if(_W_dep_[w_idx]==_W_all_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_all_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}
				}
			}else if(pcut_idx<std::distance(std::begin(_hcuts_), std::end(_hcuts_)) && flags_->Flags::Plot_Delta(species_idx)){
				//std::cout<<_sector_[sec_idx] <<" " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
				if(_hcuts_[pcut_idx]==_none_ && _cut_[cut_idx]==_no_cut_ && _top_[top_idx]==_mnone_){
					if(_W_dep_[w_idx]==_W_var_){
						for(int i=0; i<Histogram::W_bins(); i++){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s_W:%.3f-%.3f",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],Histogram::W_bot(i),Histogram::W_top(i));
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}else if(_W_dep_[w_idx]==_W_range_){
						//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_range_);
						plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
					}
					else if(_W_dep_[w_idx]==_W_all_){
						//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_all_);
						plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
					}
				}else if(_hcuts_[pcut_idx]!=_none_ && _cut_[cut_idx]!=_no_cut_){
					if(_hcuts_[pcut_idx]==_event_ && _top_[top_idx]!=_mnone_ && fun::top_perform(_top_[top_idx],flags_)){
						//std::cout<<"\tMaking Delta Event Histograms\n";
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"delta_%s_%s_%s_%s_%s_W:%.3f-%.3f",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],Histogram::W_bot(i),Histogram::W_top(i));
								plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
							}
						}else if(_W_dep_[w_idx]==_W_range_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_range_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
						else if(_W_dep_[w_idx]==_W_all_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_all_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}else if(_hcuts_[pcut_idx]!=_event_ && _top_[top_idx]==_mnone_ && fun::hcut_perform(_species_[species_idx],_hcuts_[pcut_idx],flags_)){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"delta_%s_%s_%s_%s_%s_W:%.3f-%.3f",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],Histogram::W_bot(i),Histogram::W_top(i));
								plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
							}
						}else if(_W_dep_[w_idx]==_W_range_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_range_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
						else if(_W_dep_[w_idx]==_W_all_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_sector_[sec_idx],_W_all_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}
				}
			}
			if(cart[0] == space_dims[0]-1){
				if(plot_1d.size()>0){
					plot_2d.push_back(plot_1d);
					plot_1d.clear();
				}
				if(cart[1] == space_dims[1]-1){
					if(plot_2d.size()>0){
						plot_3d.push_back(plot_2d);
						plot_2d.clear();
					}
					if(cart[2] == space_dims[2]-1){
						if(plot_3d.size()>0){
							plot_4d.push_back(plot_3d);
							plot_3d.clear();
						}
						if(cart[3] == space_dims[3]-1){
							if(plot_4d.size()>0){
								plot_5d.push_back(plot_4d);
								plot_4d.clear();
							}
							if(cart[4] == space_dims[4]-1){
								if(plot_5d.size()>0){
									_Delta_hist.push_back(plot_5d);
									plot_5d.clear();
								}
							}
						}
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::Delta_idx(float W_, const char * sector_, const char* species_, const char* pcut_, const char* cut_, const char* top_, const char* W_dep_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	idx.push_back(fun::sector_idx(sector_));
	idx.push_back(fun::species_idx(species_)+fun::species_offset(species_,_delta_cut_,flags_));
	if(species_==_ele_){
		if(pcut_==_none_ && cut_ == _no_cut_ && top_==_mnone_){
			idx.push_back(0);
			idx.push_back(0);
		}else if(pcut_!=_none_ && cut_ != _no_cut_ && fun::ecut_perform(pcut_,flags_)){
			idx.push_back(fun::ecut_idx(pcut_)+fun::ecut_offset(pcut_,flags_));
			idx.push_back(fun::cut_idx(cut_));
		}else{
			idx.push_back(-1);
			idx.push_back(-1);
		}
	}else{
		if(pcut_==_none_ && cut_ == _no_cut_ && top_==_mnone_){
			idx.push_back(0);
			idx.push_back(0);
		}else if(pcut_!=_none_ && cut_ != _no_cut_ && fun::hcut_perform(species_,pcut_,flags_)){
			idx.push_back(fun::hcut_idx(pcut_)+fun::hcut_offset(species_,pcut_,flags_));
			idx.push_back(fun::cut_idx(cut_));
		}else{
			idx.push_back(-1);
			idx.push_back(-1);
		}
	}
	if(pcut_==_event_ && top_!=_mnone_ && fun::top_perform(top_,flags_)){
		idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
	}else if(pcut_!=_event_ && top_==_mnone_){
		idx.push_back(0);
	}else{
		idx.push_back(-1);
	}
	if(W_dep_==_W_all_){
		idx.push_back(Histogram::W_bins()+1);
	}else if(cuts::in_range(W_)){
		if(W_dep_==_W_var_){
			idx.push_back(Histogram::W_bin(W_));
		}else if(W_dep_==_W_range_){
			idx.push_back(Histogram::W_bins());
		}else{
			idx.push_back(-1);
		}
	}else{
		idx.push_back(-1);
	}
	//std::cout<<"DT idx: " <<sector_ <<" " <<species_ <<" " <<pcut_ <<" " <<cut_ <<" " <<top_ <<" " <<W_dep_ <<" " <<W_ <<"\n";
	//fun::print_vector_idx(idx);
	return idx;
}


void Histogram::Delta_Fill(float p_, float dt_, float weight_, float W_, const char* species_, const char* pcut_, const char* cut_, const char* top_, const char* sector_, std::shared_ptr<Flags> flags_ ){
	if(flags_->Flags::Plot_Delta(fun::species_idx(species_))){
		std::vector<int> idx;
		//std::cout<<"Filling Delta with p:" <<p_ <<" dt:" <<dt_ <<" " <<sector_ <<" " <<species_ <<" " <<pcut_ <<" " <<cut_  <<" " <<top_  <<" W:" <<W_ <<"\n";
		if(cuts::in_range(W_)){
			//W Variance
			idx = Histogram::Delta_idx(W_,sector_,species_,pcut_,cut_,top_,_W_var_,flags_);
			//fun::print_vector_idx(idx); 
			if(Histogram::OK_Idx(idx)){
				//fun::print_vector_idx(idx);
				_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(p_,dt_,weight_);
			}
			idx.clear();
			idx = Histogram::Delta_idx(W_,_sec_all_,species_,pcut_,cut_,top_,_W_var_,flags_);
			//fun::print_vector_idx(idx); 
			if(Histogram::OK_Idx(idx)){
				//fun::print_vector_idx(idx);
				_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(p_,dt_,weight_);
			}
			idx.clear();
			//W Range
			//Sectors
			idx = Histogram::Delta_idx(W_,sector_,species_,pcut_,cut_,top_,_W_range_,flags_);
			//fun::print_vector_idx(idx); 
			if(Histogram::OK_Idx(idx)){
				//fun::print_vector_idx(idx); 
				_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(p_,dt_,weight_);
			}
			idx.clear();
			//All Sectors W Range
			idx = Histogram::Delta_idx(W_,_sec_all_,species_,pcut_,cut_,top_,_W_range_,flags_);
			//fun::print_vector_idx(idx); 
			if(Histogram::OK_Idx(idx)){
				//fun::print_vector_idx(idx); 
				_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(p_,dt_,weight_);
			}
			idx.clear();
		}
		idx = Histogram::Delta_idx(W_,sector_,species_,pcut_,cut_,top_,_W_all_,flags_);
		//fun::print_vector_idx(idx); 
		if(Histogram::OK_Idx(idx)){
			//fun::print_vector_idx(idx);
			_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(p_,dt_,weight_);
		}
		idx.clear();
		idx = Histogram::Delta_idx(W_,_sec_all_,species_,pcut_,cut_,top_,_W_all_,flags_);
		//fun::print_vector_idx(idx); 
		if(Histogram::OK_Idx(idx)){
			//fun::print_vector_idx(idx); 
			_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(p_,dt_,weight_);
		}
		idx.clear();
	}
}

void Histogram::Delta_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_Delta(0) || flags_->Flags::Plot_Delta(1) || flags_->Flags::Plot_Delta(2) || flags_->Flags::Plot_Delta(3)){
		std::cout<<"Writing Delta Histograms\n";
		char dir_name[100];
		TDirectory* dir_delta = _RootOutputFile->mkdir("Delta");
		dir_delta->cd();
		TDirectory* dir_delta_sub[std::distance(std::begin(_species_), std::end(_species_))][std::distance(std::begin(_ecuts_), std::end(_ecuts_))+1][std::distance(std::begin(_cut_), std::end(_cut_))+1][std::distance(std::begin(_top_), std::end(_top_))+1][std::distance(std::begin(_sector_),std::end(_sector_))+1][2+1];//Topology,cut,clean event, W dependence 
		//std::cout<<"\tMaking Directories\n";
		for(int i=0; i<std::distance(std::begin(_species_), std::end(_species_)); i++){
			if(flags_->Flags::Plot_Delta(i)){
				//std::cout<<"\tMaking Directory: " <<i<<" " <<0<<" " <<0<<" " <<0<<" " <<0<<" " <<0 <<"\n";
				sprintf(dir_name,"delta_%s",_species_[i]);
				dir_delta_sub[i][0][0][0][0][0] = dir_delta->mkdir(dir_name);
				if(_species_[i]==_ele_){
					for(int j=0; j<std::distance(std::begin(_ecuts_), std::end(_ecuts_)); j++){
						if(_ecuts_[j]==_none_){
							sprintf(dir_name,"delta_%s_%s_%s",_species_[i],_ecuts_[j],_no_cut_);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0<<" " <<0 <<"\n";
							dir_delta_sub[i][j+1][0][0][0][0] = dir_delta_sub[i][0][0][0][0][0]->mkdir(dir_name);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0<<" " <<1 <<"\n";
							sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_ecuts_[j],_no_cut_,_W_var_);
							dir_delta_sub[i][j+1][0][0][0][1] = dir_delta_sub[i][j+1][0][0][0][0]->mkdir(dir_name);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0<<" " <<2 <<"\n";
							sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_ecuts_[j],_no_cut_,_W_range_);
							dir_delta_sub[i][j+1][0][0][0][2] = dir_delta_sub[i][j+1][0][0][0][0]->mkdir(dir_name);
						}else if(fun::ecut_perform(_ecuts_[j],flags_)){
							sprintf(dir_name,"delta_%s_%s",_species_[i],_ecuts_[j]);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0<<" " <<0 <<"\n";
							dir_delta_sub[i][j+1][0][0][0][0] = dir_delta_sub[i][0][0][0][0][0]->mkdir(dir_name);
							for(int m=0; m<std::distance(std::begin(_sector_),std::end(_sector_)); m++){
								sprintf(dir_name,"delta_%s_%s_%s",_species_[i],_ecuts_[j],_sector_[m]);
								//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<m+1<<" " <<0 <<"\n";
								dir_delta_sub[i][j+1][0][0][m+1][0] = dir_delta_sub[i][j+1][0][0][0][0]->mkdir(dir_name);
								for(int k=0; k<std::distance(std::begin(_cut_), std::end(_cut_)); k++){
									if(_cut_[k]!=_no_cut_){
										//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<0<<" " <<m+1<<" " <<0 <<"\n";
										sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_ecuts_[j],_sector_[m],_cut_[k]);
										dir_delta_sub[i][j+1][k+1][0][m+1][0] = dir_delta_sub[i][j+1][0][0][m+1][0]->mkdir(dir_name);
										for(int l=0; l<std::distance(std::begin(_top_), std::end(_top_)); l++){
											if(fun::top_perform(_top_[l],flags_)){
												if((_ecuts_[j]==_event_ && _top_[l]!=_mnone_) || (_ecuts_[j]!=_event_ && _top_[l]==_mnone_)){
													//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<m+1<<" " <<0 <<"\n";
													sprintf(dir_name,"delta_%s_%s_%s_%s_%s",_species_[i],_ecuts_[j],_sector_[m],_cut_[k],_top_[l]);
													dir_delta_sub[i][j+1][k+1][l+1][m+1][0] = dir_delta_sub[i][j+1][k+1][0][m+1][0]->mkdir(dir_name);
													//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<m+1<<" " <<1 <<"\n";
													sprintf(dir_name,"delta_%s_%s_%s_%s_%s_%s",_species_[i],_ecuts_[j],_sector_[m],_cut_[k],_top_[l],_W_var_);
													dir_delta_sub[i][j+1][k+1][l+1][m+1][1] = dir_delta_sub[i][j+1][k+1][l+1][m+1][0]->mkdir(dir_name);
													//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<m+1<<" " <<2 <<"\n";
													sprintf(dir_name,"delta_%s_%s_%s_%s_%s_%s",_species_[i],_ecuts_[j],_sector_[m],_cut_[k],_top_[l],_W_range_);
													dir_delta_sub[i][j+1][k+1][l+1][m+1][2] = dir_delta_sub[i][j+1][k+1][l+1][m+1][0]->mkdir(dir_name);
												}
											}
										}
									}
								}
							}
						}
					}
				}else{
					for(int j=0; j<std::distance(std::begin(_hcuts_), std::end(_hcuts_)); j++){
						if(_hcuts_[j]==_none_){
							sprintf(dir_name,"delta_%s_%s_%s",_species_[i],_hcuts_[j],_no_cut_);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0<<" " <<0 <<"\n";
							dir_delta_sub[i][j+1][0][0][0][0] = dir_delta_sub[i][0][0][0][0][0]->mkdir(dir_name);
							sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_hcuts_[j],_no_cut_,_W_var_);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0<<" " <<1 <<"\n";
							dir_delta_sub[i][j+1][0][0][0][1] = dir_delta_sub[i][j+1][0][0][0][0]->mkdir(dir_name);
							sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_hcuts_[j],_no_cut_,_W_range_);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0 <<" " <<0<<" " <<0<<" " <<2 <<"\n";
							dir_delta_sub[i][j+1][0][0][0][2] = dir_delta_sub[i][j+1][0][0][0][0]->mkdir(dir_name);
							for(int e=0; e<std::distance(std::begin(_sector_), std::end(_sector_)); e++){
								//sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_hcuts_[j],_no_cut_,_sector_[e]);
								//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<e+1<<" " <<0 <<"\n";
								//dir_delta_sub[i][j+1][0][0][e+1][0] = dir_delta_sub[i][0][0][0][0][0]->mkdir(dir_name);
								sprintf(dir_name,"delta_%s_%s_%s_%s_%s",_species_[i],_hcuts_[j],_no_cut_,_sector_[e],_W_var_);
								//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<e+1<<" " <<1 <<"\n";
								dir_delta_sub[i][j+1][0][0][e+1][1] = dir_delta_sub[i][j+1][0][0][0][1]->mkdir(dir_name);
								sprintf(dir_name,"delta_%s_%s_%s_%s_%s",_species_[i],_hcuts_[j],_no_cut_,_sector_[e],_W_range_);
								//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0 <<" " <<0<<" " <<e+1<<" " <<2 <<"\n";
								dir_delta_sub[i][j+1][0][0][e+1][2] = dir_delta_sub[i][j+1][0][0][0][2]->mkdir(dir_name);
							}
						}else if(fun::hcut_perform(_species_[i],_hcuts_[j],flags_)){
							sprintf(dir_name,"delta_%s_%s",_species_[i],_hcuts_[j]);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0<<" " <<0 <<"\n";
							dir_delta_sub[i][j+1][0][0][0][0] = dir_delta_sub[i][0][0][0][0][0]->mkdir(dir_name);
							for(int m=0; m<std::distance(std::begin(_sector_),std::end(_sector_)); m++){
								sprintf(dir_name,"delta_%s_%s_%s",_species_[i],_hcuts_[j],_sector_[m]);
								//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<m+1<<" " <<0 <<"\n";
								dir_delta_sub[i][j+1][0][0][m+1][0] = dir_delta_sub[i][j+1][0][0][0][0]->mkdir(dir_name);
								for(int k=0; k<std::distance(std::begin(_cut_), std::end(_cut_)); k++){
									if(_cut_[k]!=_no_cut_){
										sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_hcuts_[j],_sector_[m],_cut_[k]);
										//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<0<<" " <<m+1<<" " <<0 <<"\n";
										dir_delta_sub[i][j+1][k+1][0][m+1][0] = dir_delta_sub[i][j+1][0][0][m+1][0]->mkdir(dir_name);
										for(int l=0; l<std::distance(std::begin(_top_), std::end(_top_)); l++){
											if(fun::top_perform(_top_[l],flags_)){
												if((_hcuts_[j]==_event_ && _top_[l]!=_mnone_) || (_hcuts_[j]!=_event_ && _top_[l]==_mnone_)){
													sprintf(dir_name,"delta_%s_%s_%s_%s_%s",_species_[i],_hcuts_[j],_sector_[m],_cut_[k],_top_[l]);
													//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<m+1<<" " <<0 <<"\n";
													dir_delta_sub[i][j+1][k+1][l+1][m+1][0] = dir_delta_sub[i][j+1][k+1][0][m+1][0]->mkdir(dir_name);
													sprintf(dir_name,"delta_%s_%s_%s_%s_%s_%s",_species_[i],_hcuts_[j],_sector_[m],_cut_[k],_top_[l],_W_var_);
													//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<m+1<<" " <<1 <<"\n";
													dir_delta_sub[i][j+1][k+1][l+1][m+1][1] = dir_delta_sub[i][j+1][k+1][l+1][m+1][0]->mkdir(dir_name);
													sprintf(dir_name,"delta_%s_%s_%s_%s_%s_%s",_species_[i],_hcuts_[j],_sector_[m],_cut_[k],_top_[l],_W_range_);
													//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<m+1<<" " <<2 <<"\n";
													dir_delta_sub[i][j+1][k+1][l+1][m+1][2] = dir_delta_sub[i][j+1][k+1][l+1][m+1][0]->mkdir(dir_name);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		std::vector<long> space_dims(6);//{species,pcut,cut,top,w_dep}
		space_dims[5] = std::distance(std::begin(_sector_),std::end(_sector_));
		space_dims[4] = std::distance(std::begin(_species_), std::end(_species_));
		space_dims[3] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));
		space_dims[1] = std::distance(std::begin(_top_), std::end(_top_));
		space_dims[0] = std::distance(std::begin(_W_dep_), std::end(_W_dep_));
		int species_idx = -1;
		int pcut_idx = -1;
		int cut_idx = -1; 
		int top_idx = -1;
		int w_idx = -1; 
		int sec_idx = -1;
		CartesianGenerator cart(space_dims);
		std::vector<int> idx;
		//std::cout<<"\tWriting Histograms\n";
		while(cart.GetNextCombination()){
			species_idx = cart[4];
			pcut_idx = cart[3];
			cut_idx = cart[2];
			top_idx = cart[1];
			w_idx = cart[0];
			sec_idx = cart[5];
			//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
			if(flags_->Plot_Delta(species_idx)){
				if(_species_[species_idx]==_ele_){
					if(_ecuts_[pcut_idx]==_none_ && _cut_[cut_idx]==_no_cut_ && _top_[top_idx]==_mnone_){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								idx = Histogram::Delta_idx(Histogram::W_center(i),_sector_[sec_idx],_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
									//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<0 <<" " <<0 <<" " <<sec_idx+1 <<" " <<1 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][0][0][sec_idx+1][1]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
								}
								idx.clear();
							}
						}else if(_W_dep_[w_idx]==_W_range_){
							idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
							if(Histogram::OK_Idx(idx)){
								//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
								//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<0 <<" " <<0 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
								dir_delta_sub[species_idx][pcut_idx+1][0][0][sec_idx+1][2]->cd();
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
							}
							idx.clear();
						}else if(_W_dep_[w_idx]==_W_all_){
							idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
							if(Histogram::OK_Idx(idx)){
								//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
								dir_delta_sub[species_idx][pcut_idx+1][0][0][sec_idx+1][2]->cd();
								//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<0 <<" " <<0 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
							}
							idx.clear();
						}
					}else if(_ecuts_[pcut_idx]!=_none_ && _cut_[cut_idx]!=_no_cut_){
						if(_ecuts_[pcut_idx]==_event_ && _top_[top_idx]!=_mnone_ && fun::top_perform(_top_[top_idx],flags_)){
							if(_W_dep_[w_idx]==_W_var_){
								for(int i=0; i<Histogram::W_bins(); i++){
									idx = Histogram::Delta_idx(Histogram::W_center(i),_sector_[sec_idx],_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
									if(Histogram::OK_Idx(idx)){
										//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
										//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<1 <<"\n";
										//fun::print_vector_idx(idx);
										//std::cout<<"\tPre directory: "<<species_idx<<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<1 <<"\n";
										dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][1]->cd();
										//std::cout<<"\tPost directory\n";
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
										//std::cout<<"\tPost Writing\n";
									}
									idx.clear();
								}
							}else if(_W_dep_[w_idx]==_W_range_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									//std::cout<<"\tPre directory: "<<species_idx<<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<2 <<"\n";
									//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][2]->cd();
									//std::cout<<"\tPost directory\n";
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
									//std::cout<<"\tPost writing\n";
								}
								idx.clear();
							}else if(_W_dep_[w_idx]==_W_all_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									//std::cout<<"\tPre directory: "<<species_idx<<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<1 <<"\n";
									//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][2]->cd();
									//std::cout<<"\tPost directory\n";
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
									//std::cout<<"\tPost write\n";
								}
								idx.clear();
							}
						}else if(_ecuts_[pcut_idx]!=_event_ && _top_[top_idx]==_mnone_ && fun::ecut_perform(_ecuts_[pcut_idx],flags_)){
							if(_W_dep_[w_idx]==_W_var_){
								for(int i=0; i<Histogram::W_bins(); i++){
									idx = Histogram::Delta_idx(Histogram::W_center(i),_sector_[sec_idx],_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
									if(Histogram::OK_Idx(idx)){
										//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
										//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<1 <<"\n";
										dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][1]->cd();
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
									}
									idx.clear();
								}
							}else if(_W_dep_[w_idx]==_W_range_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
								}
								idx.clear();
							}else if(_W_dep_[w_idx]==_W_all_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
								}
								idx.clear();
							}
						}
					}
				}else if(pcut_idx<std::distance(std::begin(_hcuts_), std::end(_hcuts_))){
					if(_hcuts_[pcut_idx]==_none_ && _cut_[cut_idx]==_no_cut_ && _top_[top_idx]==_mnone_){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								idx = Histogram::Delta_idx(Histogram::W_center(i),_sector_[sec_idx],_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
									//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<0 <<" " <<0 <<" " <<sec_idx+1 <<" " <<1 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][0][0][sec_idx+1][1]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
								}
								idx.clear();
							}
						}else if(_W_dep_[w_idx]==_W_range_){
							idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
							if(Histogram::OK_Idx(idx)){
								//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
								//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<0 <<" " <<0 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
								dir_delta_sub[species_idx][pcut_idx+1][0][0][sec_idx+1][2]->cd();
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
							}
							idx.clear();
						}else if(_W_dep_[w_idx]==_W_var_){
							idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
							if(Histogram::OK_Idx(idx)){
								//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
								//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<0 <<" " <<0 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
								dir_delta_sub[species_idx][pcut_idx+1][0][0][sec_idx+1][2]->cd();
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
							}
							idx.clear();
						}
					}else if(_hcuts_[pcut_idx]!=_none_ && _cut_[cut_idx]!=_no_cut_){
						if(_hcuts_[pcut_idx]==_event_ && _top_[top_idx]!=_mnone_ && fun::top_perform(_top_[top_idx],flags_)){
							if(_W_dep_[w_idx]==_W_var_){
								for(int i=0; i<Histogram::W_bins(); i++){
									idx = Histogram::Delta_idx(Histogram::W_center(i),_sector_[sec_idx],_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
									if(Histogram::OK_Idx(idx)){
										//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
										//fun::print_vector_idx(idx);
										//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<1 <<"\n";
										dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][1]->cd();
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
									}
									idx.clear();
								}
							}else if(_W_dep_[w_idx]==_W_range_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
								}
								idx.clear();
							}else if(_W_dep_[w_idx]==_W_all_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
								}
								idx.clear();
							}
						}else if(_hcuts_[pcut_idx]!=_event_ && _top_[top_idx]==_mnone_ && fun::hcut_perform(_species_[species_idx],_hcuts_[pcut_idx],flags_)){
							if(_W_dep_[w_idx]==_W_var_){
								for(int i=0; i<Histogram::W_bins(); i++){
									idx = Histogram::Delta_idx(Histogram::W_center(i),_sector_[sec_idx],_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
									if(Histogram::OK_Idx(idx)){
										//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
										//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<1 <<"\n";
										dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][1]->cd();
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
									}
									idx.clear();
								}
							}else if(_W_dep_[w_idx]==_W_range_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
								}
								idx.clear();
							}else if(_W_dep_[w_idx]==_W_all_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_sector_[sec_idx],_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									//std::cout<<"\t" <<species_idx <<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<sec_idx+1 <<" " <<2 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][sec_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
								}
								idx.clear();
							}
						}
					}
				}
			}
		}
	}
}

 //*-------------------------------End Delta T Plot------------------------------*

 //*-------------------------------Start MM Plot----------------------------*
void Histogram::MM_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_MM(0) || flags_->Flags::Plot_MM(1) || flags_->Flags::Plot_MM(2) || flags_->Flags::Plot_MM(3)){
		std::cout<<"Making MM histograms\n";
		Bool_3d made_3d;
		Bool_2d made_2d;
		Bool_1d made_1d;
		TH1F_ptr_4d plot_4d;
		TH1F_ptr_3d plot_3d;
		TH1F_ptr_2d plot_2d;
		TH1F_ptr_1d plot_1d;
		TH1F_ptr_4d plot_4d_lin;
		TH1F_ptr_3d plot_3d_lin;
		TH1F_ptr_2d plot_2d_lin;
		TH1F_ptr_1d plot_1d_lin;

		std::vector<long> space_dims(5);
		space_dims[4] = std::distance(std::begin(_sector_), std::end(_sector_));
		space_dims[3] = std::distance(std::begin(_top_), std::end(_top_));//{ele,pro,pip,pim}
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));//Electron Cuts (electrons have more cuts than hadrons)
		space_dims[1] = std::distance(std::begin(_clean_event_),std::end(_clean_event_));//Clean event or not
		space_dims[0]	= (Histogram::W_bins()+2);//W Binning + All in range + All period
		//std::cout<<"W bins+2: " <<space_dims[0];
 		CartesianGenerator cart(space_dims);
		CartesianGenerator cart2(space_dims);
		char hname[100];//For naming histograms
		//std::cout<<"W Bins: " <<Histogram::W_bins() <<"\n";
		/*while(cart2.GetNextCombination()){
			//if(_top_[cart2[3]]!=_mnone_ || _top_[cart2[3]]!=_mall_){
				made_1d.push_back(false);
				if(cart2[0]== space_dims[0]-1 && made_1d.size()>0){
					made_2d.push_back(made_1d);
					made_1d.clear();
				}
				if(cart2[1] == space_dims[1]-1 && made_2d.size()>0 && cart2[0]== space_dims[0]-1){
						made_3d.push_back(made_2d);
						made_2d.clear();
				}
				if(cart2[2] == space_dims[2]-1 && made_3d.size()>0 && cart2[1] == space_dims[1]-1 && cart2[0]== space_dims[0]-1){
						_MM_made.push_back(made_3d);
						made_3d.clear();
				}
			//}
		}*/
		while(cart.GetNextCombination()){
			if(fun::top_perform(_top_[cart[3]],flags_)){
				if(_top_[cart[3]]!=_mnone_ && _top_[cart[3]]!=_mall_ && _clean_event_[cart[1]]!=_isolated_){
					if(cart[0]< Histogram::W_bins()){
						//std::cout<<"MM idx: " <<_MM_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"MM_%s_%s_%s_%s_W:%.3f-%.3f",_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],_sector_[cart[4]],Histogram::W_bot(cart[0]),Histogram::W_top(cart[0]));
						plot_1d.push_back(new TH1F(hname,hname,_mm2_bin_[cart[3]],_mm2_min_[cart[3]],_mm2_max_[cart[3]]));
						sprintf(hname,"MM_linear_%s_%s_%s_%s_W:%.3f-%.3f",_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],_sector_[cart[4]],Histogram::W_bot(cart[0]),Histogram::W_top(cart[0]));
						plot_1d_lin.push_back(new TH1F(hname,hname,_mm_bin_[cart[3]],_mm_min_[cart[3]],_mm_max_[cart[3]]));
						//_MM_made[cart[3]][cart[2]][cart[1]][cart[0]]=true;
					}else if(cart[0] == Histogram::W_bins()){
						//std::cout<<"MM idx: " <<_MM_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"MM_%s_%s_%s_%s_W:%s",_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],_sector_[cart[4]],"in_range");
						plot_1d.push_back(new TH1F(hname,hname,_mm2_bin_[cart[3]],_mm2_min_[cart[3]],_mm2_max_[cart[3]]));
						sprintf(hname,"MM_linear_%s_%s_%s_%s_W:%s",_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],_sector_[cart[4]],"in_range");
						plot_1d_lin.push_back(new TH1F(hname,hname,_mm_bin_[cart[3]],_mm_min_[cart[3]],_mm_max_[cart[3]]));
						//_MM_made[cart[3]][cart[2]][cart[1]][cart[0]]=true;
					}else if(cart[0] == Histogram::W_bins()+1){
						//std::cout<<"MM idx: " <<_MM_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"MM_%s_%s_%s_%s_W:%s",_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],_sector_[cart[4]],"all");
						plot_1d.push_back(new TH1F(hname,hname,_mm2_bin_[cart[3]],_mm2_min_[cart[3]],_mm2_max_[cart[3]]));
						sprintf(hname,"MM_linear_%s_%s_%s_%s_W:%s",_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],_sector_[cart[4]],"all");
						plot_1d_lin.push_back(new TH1F(hname,hname,_mm_bin_[cart[3]],_mm_min_[cart[3]],_mm_max_[cart[3]]));
						//_MM_made[cart[3]][cart[2]][cart[1]][cart[0]]=true;
					}
				}
			}
			if(cart[0] == space_dims[0]-1){
				if(plot_1d.size()>0){
					plot_2d.push_back(plot_1d);
					plot_1d.clear();
				}
				if(plot_1d_lin.size()>0){
					plot_2d_lin.push_back(plot_1d_lin);
					plot_1d_lin.clear();
				}
				if(cart[1] == space_dims[1]-1){
					if(plot_2d.size()>0){
						plot_3d.push_back(plot_2d);
						plot_2d.clear();
					}
					if(plot_2d_lin.size()>0){
						plot_3d_lin.push_back(plot_2d_lin);
						plot_2d_lin.clear();
					}
					if(cart[2] == space_dims[2]-1){
						if(plot_3d.size()>0){
							plot_4d.push_back(plot_3d);
							plot_3d.clear();
						}
						if(plot_3d_lin.size()>0){
							plot_4d_lin.push_back(plot_3d_lin);
							plot_3d_lin.clear();
						}
						if(cart[3]== space_dims[3]-1){
							if(plot_4d.size()>0){
								_MM_hist.push_back(plot_4d);
								plot_4d.clear();
							}
							if(plot_4d_lin.size()>0){
								_MM_lin_hist.push_back(plot_4d_lin);
								plot_4d_lin.clear();
							}
						}
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::MM_idx(const char* top_, const char* cut_, const char * clean_, const char * sector_, const char * W_dep_, float W_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx; 
	idx.push_back(fun::sector_idx(sector_));
	if(top_!=_mnone_ &&  top_!=_mall_){
		idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
		idx.push_back(fun::cut_idx(cut_));
		idx.push_back(fun::clean_idx(clean_));
		if(W_dep_ == _W_var_){
			idx.push_back(Histogram::W_bin(W_));
		}else if(W_dep_ == _W_range_){
			idx.push_back(Histogram::W_bins());
		}else if(W_dep_ == _W_all_){
			idx.push_back(Histogram::W_bins()+1);
		}else{
			idx.push_back(-1);
		}
	}else{
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
	}
	//std::cout<<"MM idx: ";
	//fun::print_vector_idx(idx);
	return idx; 
}

void Histogram::MM_Fill(const char* top_, const char* cut_, const char * clean_, const char * sector_, float MM_, float W_, float weight_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_MM(fun::top_idx(top_))){return;}
	//std::cout<<"\tTrying to fill MM top:" <<top_ <<" cut:" <<cut_ <<" clean:" <<clean_ <<" MM:" <<MM_ <<" W:" <<W_ <<"\n";
	std::vector<int> idx; 
	float MM_lin = -99.9;
	if(MM_ > 0.0){
		MM_lin = TMath::Sqrt(MM_);
	}else{
		MM_lin = -TMath::Sqrt(-MM_);
	}
	if(flags_->Flags::Plot_MM(fun::top_idx(top_)) && fun::top_perform(top_,flags_)){
		if(top_ != _mnone_ && top_ != _mall_ && clean_!=_isolated_){
			//std::cout<<"W: " <<W_ <<"  => W Bin: " <<Histogram::W_bin(W_) <<"\n";
			if(Histogram::W_bin(W_)>=0){
				//if(_MM_made[fun::top_idx(top_)][fun::cut_idx(cut_)][fun::clean_idx(clean_)][Histogram::W_bin(W_)]){//W Binning
					idx = Histogram::MM_idx(top_,cut_,clean_,sector_,_W_var_,W_,flags_);
					if(Histogram::OK_Idx(idx)){
						_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_,weight_);
						_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_lin,weight_);
					}
					idx.clear();
					idx = Histogram::MM_idx(top_,cut_,clean_,_sec_all_,_W_var_,W_,flags_);
					if(Histogram::OK_Idx(idx)){
						_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_,weight_);
						_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_lin,weight_);
					}
					idx.clear();
				//}
			}
			//if(_MM_made[fun::top_idx(top_)][fun::cut_idx(cut_)][fun::clean_idx(clean_)][Histogram::W_bins()]){//W in range
				idx = Histogram::MM_idx(top_,cut_,clean_,sector_,_W_range_,W_,flags_);
				if(Histogram::OK_Idx(idx)){
					_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_,weight_);
					_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_lin,weight_);
				}
				idx.clear();
				idx = Histogram::MM_idx(top_,cut_,clean_,_sec_all_,_W_range_,W_,flags_);
				if(Histogram::OK_Idx(idx)){
					_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_,weight_);
					_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_lin,weight_);
				}
				idx.clear();
			//}
			//if(_MM_made[fun::top_idx(top_)][fun::cut_idx(cut_)][fun::clean_idx(clean_)][Histogram::W_bins()+1]){//All W
				idx = Histogram::MM_idx(top_,cut_,clean_,sector_,_W_all_,W_,flags_);
				if(Histogram::OK_Idx(idx)){
					_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_,weight_);
					_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_lin,weight_);
				}
				idx.clear();
				idx = Histogram::MM_idx(top_,cut_,clean_,_sec_all_,_W_all_,W_,flags_);
				if(Histogram::OK_Idx(idx)){
					_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_,weight_);
					_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(MM_lin,weight_);
				}
				idx.clear();
			//}
		}
	}
}

void Histogram::MM_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_MM(0) || flags_->Flags::Plot_MM(1) || flags_->Flags::Plot_MM(2) || flags_->Flags::Plot_MM(3)){
		std::cout<<"Writing MM Plots\n";
		char dir_name[100];
		TDirectory* dir_mm = _RootOutputFile->mkdir("MM");
		dir_mm->cd();
		TDirectory* dir_mm_sub[4][std::distance(std::begin(_cut_), std::end(_cut_))+1][2+1][2+1][std::distance(std::begin(_sector_), std::end(_sector_))+1];//Topology,cut,clean event, W dependence 
		for(int i=0; i<4; i++){
			if(fun::top_perform(_top_[i],flags_) && flags_->Plot_MM(i)){
				sprintf(dir_name,"MM_%s",_top_[i]);
				//std::cout<<"\tMaking Dir: " <<dir_name <<" " <<i <<" " <<0 <<" " <<0 <<" " <<0 <<" " <<0 <<"\n";
				dir_mm_sub[i][0][0][0][0] = dir_mm->mkdir(dir_name);
				for(int j=0; j<std::distance(std::begin(_cut_), std::end(_cut_)); j++){
					sprintf(dir_name,"MM_%s_%s",_top_[i],_cut_[j]);
					//std::cout<<"\t\tMaking Dir: " <<dir_name <<" " <<i <<" " <<j+1 <<" " <<0 <<" " <<0 <<" " <<0 <<"\n";
					dir_mm_sub[i][1+j][0][0][0] = dir_mm_sub[i][0][0][0][0]->mkdir(dir_name);
					for(int l=0; l<std::distance(std::begin(_sector_), std::end(_sector_)); l++){
						sprintf(dir_name,"MM_%s_%s_%s",_top_[i],_cut_[j],_sector_[l]);
						//std::cout<<"\t\t\tMaking Dir: " <<dir_name <<" " <<i <<" " <<j+1 <<" " <<0 <<" " <<0 <<" " <<l+1 <<"\n";
						dir_mm_sub[i][1+j][0][0][l+1] = dir_mm_sub[i][j+1][0][0][0]->mkdir(dir_name);
						for(int k=0; k<(std::distance(std::begin(_clean_event_),std::end(_clean_event_))-1); k++){
							sprintf(dir_name,"MM_%s_%s_%s_%s",_top_[i],_cut_[j],_sector_[l],_clean_event_[k]);
							//std::cout<<"\t\t\t\tMaking Dir: " <<dir_name <<" " <<i <<" " <<j+1 <<" " <<k+1 <<" " <<0 <<" " <<l+1 <<"\n";
							dir_mm_sub[i][j+1][k+1][0][l+1] = dir_mm_sub[i][j+1][0][0][l+1]->mkdir(dir_name); 
							sprintf(dir_name,"MM_%s_%s_%s_%s_W_dep",_top_[i],_cut_[j],_sector_[l],_clean_event_[k]);
							//std::cout<<"\t\t\t\tMaking Dir: " <<dir_name <<" " <<i <<" " <<j1 <<" " <<k+1 <<" " <<1 <<" " <<l+1 <<"\n";
							dir_mm_sub[i][j+1][k+1][1][l+1] = dir_mm_sub[i][j+1][k+1][0][l+1]->mkdir(dir_name);
							sprintf(dir_name,"MM_%s_%s_%s_%s_W_Range",_top_[i],_cut_[j],_sector_[l],_clean_event_[k]);
							//std::cout<<"\t\t\t\tMaking Dir: " <<dir_name <<" " <<i <<" " <<j+1 <<" " <<k+1 <<" " <<2 <<" " <<l+1 <<"\n";
							dir_mm_sub[i][j+1][k+1][2][l+1] = dir_mm_sub[i][j+1][k+1][0][l+1]->mkdir(dir_name);
						}
					}
				}
			}
		}

		std::vector<long> space_dims(5);
		space_dims[4] = std::distance(std::begin(_sector_), std::end(_sector_));
		space_dims[3] = std::distance(std::begin(_top_), std::end(_top_));//{pro,pip,pim,zero}
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));//Electron Cuts (electrons have more cuts than hadrons)
		space_dims[1] = std::distance(std::begin(_clean_event_),std::end(_clean_event_));//Clean event or not
		space_dims[0]	= Histogram::W_bins()+2;//W Binning + All in range + All period 
		CartesianGenerator cart(space_dims);

		std::vector<int> idx;
		while(cart.GetNextCombination()){
			if(flags_->Flags::Plot_MM(cart[3])){
				if(fun::top_perform(_top_[cart[3]],flags_)){
					if(_top_[cart[3]]!=_mnone_ && _top_[cart[3]]!=_mall_ && _clean_event_[cart[1]]!=_isolated_){
						if(cart[0] >= Histogram::W_bins()){// && _MM_made[cart[3]][cart[2]][cart[1]][cart[0]]){
							//std::cout<<"Writing in dir: " <<cart[3] <<" " <<cart[2]+1 <<" " <<cart[1]+1 <<" " <<2 <<" " <<cart[4]+1 <<"\n";
							dir_mm_sub[cart[3]][cart[2]+1][cart[1]+1][2][cart[4]+1]->cd();
							idx=Histogram::MM_idx(_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],_sector_[cart[4]],_W_dep_[cart[0]-Histogram::W_bins()],Histogram::W_center(cart[0]),flags_);
							if(Histogram::OK_Idx(idx)){
								_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("MM^2 (GeV^{2})");
								_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("MM (GeV)");
								if(flags_->Flags::Sim()){
									_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Weighted Yield");
									_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Weighted Yield");
								}else{
									_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Yield");
									_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Yield");
								}
								_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
								_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
							}
							idx.clear();
						}else {// if(_MM_made[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]){
							//std::cout<<"Writing in dir: " <<cart[3] <<" " <<cart[2]+1 <<" " <<cart[1]+1 <<" " <<1 <<" " <<cart[4]+1 <<"\n";
							dir_mm_sub[cart[3]][cart[2]+1][cart[1]+1][1][cart[4]+1]->cd();
							idx=Histogram::MM_idx(_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],_sector_[cart[4]],_W_var_,Histogram::W_center(cart[0]),flags_);
							if(Histogram::OK_Idx(idx)){
								_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("MM^2 (GeV^{2})");
								_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("MM (GeV)");
								if(flags_->Flags::Sim()){
									_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Weighted Yield");
									_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Weighted Yield");
								}else{
									_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Yield");
									_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Yield");
								}
								_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
								_MM_lin_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
							}
							idx.clear();
						}
					}
				}
			}
		}
	}
}


 //*-------------------------------End MM Plot------------------------------*

 //*-------------------------------Start Kinematic Eff Plot------------------------------*
void Histogram::Kinematic_Eff_Make(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Kin_Eff(0) && !flags_->Flags::Plot_Kin_Eff(1) && !flags_->Flags::Plot_Kin_Eff(2) && !flags_->Flags::Plot_Kin_Eff(3)){ return;}
	std::cout<<"Making Kinematic Efficiency Plots\n";
	/*
	Histogram Splits
	- Electron Cuts
	- Cut Status
	- Topology
	*/
	TH2F_ptr_1d plot_1d_1;
	TH2F_ptr_2d plot_2d_1;
	TH2F_ptr_3d plot_3d_1;
	TH2F_ptr_4d plot_4d_1;
	TH1F_ptr_1d plot_1d_2;
	TH1F_ptr_2d plot_2d_2;
	TH1F_ptr_3d plot_3d_2;

	std::vector<long> space_dims(5);
	space_dims[4] = std::distance(std::begin(_species_), std::end(_species_));
	space_dims[3] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
	space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));
	space_dims[1] = std::distance(std::begin(_sector_),std::end(_sector_));
	space_dims[0] = std::distance(std::begin(_top_), std::end(_top_));
	char hname[100];
	CartesianGenerator cart(space_dims);
	int ecut_idx;
	int cut_idx;
	int sec_idx;
	int top_idx;
	while(cart.GetNextCombination()){
		ecut_idx = cart[3];
		cut_idx = cart[2];
		sec_idx = cart[1];
		top_idx = cart[0];
		//std::cout<<"Looking idx: " <<_ecuts_[ecut_idx] <<" " <<_cut_[cut_idx] <<" " <<_sector_[sec_idx] <<" " <<_top_[top_idx] <<" " <<fun::ecut_perform(_ecuts_[ecut_idx],flags_) <<" " <<(_ecuts_[ecut_idx] == _none_)  <<(_cut_[cut_idx] == _no_cut_)  <<(_top_[top_idx]==_none_)<<"\n";
		if(_species_[cart[4]]==_ele_ && flags_->Flags::Plot_Kin_Eff(0)){	
			if(fun::ecut_perform(_ecuts_[ecut_idx],flags_)){
				if(_ecuts_[ecut_idx] == _none_ && _cut_[cut_idx] == _no_cut_ && _top_[top_idx]==_mnone_){
					//std::cout<<"K Eff1 Idx: " <<_Kinematic_Eff1_hist.size() <<" " <<plot_4d_1.size() <<" " <<plot_3d_1.size() <<" " <<plot_2d_1.size() <<" " <<plot_1d_1.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
					sprintf(hname,"Kinematic_eff_%s_%s_%s_%s_%s",_species_[cart[4]],_ecuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx]);
					plot_1d_1.push_back(new TH2F(hname,hname,_cc_eff1_xbin_,_cc_eff1_xmin_,_cc_eff1_xmax_,_cc_eff1_ybin_,_cc_eff1_ymin_,_cc_eff1_ymax_));
					if(_sector_[sec_idx]==_sec_all_){
						//std::cout<<"\tK Eff2 Idx: " <<_Kinematic_Eff2_hist.size() <<" " <<plot_3d_2.size() <<" " <<plot_2d_2.size() <<" " <<plot_1d_2.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
						sprintf(hname,"Kinematic_eff_sec_yields_%s_%s_%s_%s",_species_[cart[4]],_ecuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx]);
						plot_1d_2.push_back(new TH1F(hname,hname,_cc_eff2_xbin_,_cc_eff2_xmin_,_cc_eff2_xmax_));
					}
				}else if(_ecuts_[ecut_idx] != _none_ && _cut_[cut_idx] != _no_cut_){
					if(_ecuts_[ecut_idx] != _event_ && _top_[top_idx] == _mnone_){
						//std::cout<<"K Eff1 Idx: " <<_Kinematic_Eff1_hist.size() <<" " <<plot_4d_1.size() <<" " <<plot_3d_1.size() <<" " <<plot_2d_1.size() <<" " <<plot_1d_1.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
						sprintf(hname,"Kinematic_eff_%s_%s_%s_%s_%s",_species_[cart[4]],_ecuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx]);
						plot_1d_1.push_back(new TH2F(hname,hname,_cc_eff1_xbin_,_cc_eff1_xmin_,_cc_eff1_xmax_,_cc_eff1_ybin_,_cc_eff1_ymin_,_cc_eff1_ymax_));
						if(_sector_[sec_idx]==_sec_all_){
						//	std::cout<<"\tK Eff2 Idx: " <<_Kinematic_Eff2_hist.size() <<" " <<plot_3d_2.size() <<" " <<plot_2d_2.size() <<" " <<plot_1d_2.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
							sprintf(hname,"Kinematic_eff_sec_yields_%s_%s_%s_%s",_species_[cart[4]],_ecuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx]);
							plot_1d_2.push_back(new TH1F(hname,hname,_cc_eff2_xbin_,_cc_eff2_xmin_,_cc_eff2_xmax_));
						}
					}else if(_ecuts_[ecut_idx] == _event_ && _top_[top_idx] != _mnone_){
						//std::cout<<"K Eff1 Idx: " <<_Kinematic_Eff1_hist.size() <<" " <<plot_4d_1.size() <<" " <<plot_3d_1.size() <<" " <<plot_2d_1.size() <<" " <<plot_1d_1.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
						sprintf(hname,"Kinematic_eff_%s_%s_%s_%s_%s",_species_[cart[4]],_ecuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx]);
						plot_1d_1.push_back(new TH2F(hname,hname,_cc_eff1_xbin_,_cc_eff1_xmin_,_cc_eff1_xmax_,_cc_eff1_ybin_,_cc_eff1_ymin_,_cc_eff1_ymax_));
						if(_sector_[sec_idx]==_sec_all_){
							//std::cout<<"\tK Eff2 Idx: " <<_Kinematic_Eff2_hist.size() <<" " <<plot_3d_2.size() <<" " <<plot_2d_2.size() <<" " <<plot_1d_2.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
							sprintf(hname,"Kinematic_eff_sec_yields_%s_%s_%s_%s",_species_[cart[4]],_ecuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx]);
							plot_1d_2.push_back(new TH1F(hname,hname,_cc_eff2_xbin_,_cc_eff2_xmin_,_cc_eff2_xmax_));
						}
					}
				}
			}
		}else if(ecut_idx<std::distance(std::begin(_hcuts_), std::end(_hcuts_)) && flags_->Plot_Kin_Eff(cart[4])){
			if(fun::hcut_perform(_species_[cart[4]],_hcuts_[ecut_idx],flags_)){
				if(_hcuts_[ecut_idx] == _none_ && _cut_[cut_idx] == _no_cut_ && _top_[top_idx]==_mnone_){
					//std::cout<<"K Eff1 Idx: " <<_Kinematic_Eff1_hist.size() <<" " <<plot_4d_1.size() <<" " <<plot_3d_1.size() <<" " <<plot_2d_1.size() <<" " <<plot_1d_1.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
					sprintf(hname,"Kinematic_eff_%s_%s_%s_%s_%s",_species_[cart[4]],_hcuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx]);
					plot_1d_1.push_back(new TH2F(hname,hname,_cc_eff1_xbin_,_cc_eff1_xmin_,_cc_eff1_xmax_,_cc_eff1_ybin_,_cc_eff1_ymin_,_cc_eff1_ymax_));
					if(_sector_[sec_idx]==_sec_all_){
						//std::cout<<"\tK Eff2 Idx: " <<_Kinematic_Eff2_hist.size() <<" " <<plot_3d_2.size() <<" " <<plot_2d_2.size() <<" " <<plot_1d_2.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
						sprintf(hname,"Kinematic_eff_sec_yields_%s_%s_%s_%s",_species_[cart[4]],_hcuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx]);
						plot_1d_2.push_back(new TH1F(hname,hname,_cc_eff2_xbin_,_cc_eff2_xmin_,_cc_eff2_xmax_));
					}
				}else if(_hcuts_[ecut_idx] != _none_ && _cut_[cut_idx] != _no_cut_){
					if(_hcuts_[ecut_idx] != _event_ && _top_[top_idx] == _mnone_){
						//std::cout<<"K Eff1 Idx: " <<_Kinematic_Eff1_hist.size() <<" " <<plot_4d_1.size() <<" " <<plot_3d_1.size() <<" " <<plot_2d_1.size() <<" " <<plot_1d_1.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
						sprintf(hname,"Kinematic_eff_%s_%s_%s_%s_%s",_species_[cart[4]],_hcuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx]);
						plot_1d_1.push_back(new TH2F(hname,hname,_cc_eff1_xbin_,_cc_eff1_xmin_,_cc_eff1_xmax_,_cc_eff1_ybin_,_cc_eff1_ymin_,_cc_eff1_ymax_));
						if(_sector_[sec_idx]==_sec_all_){
							//std::cout<<"\tK Eff2 Idx: " <<_Kinematic_Eff2_hist.size() <<" " <<plot_3d_2.size() <<" " <<plot_2d_2.size() <<" " <<plot_1d_2.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
							sprintf(hname,"Kinematic_eff_sec_yields_%s_%s_%s_%s",_species_[cart[4]],_hcuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx]);
							plot_1d_2.push_back(new TH1F(hname,hname,_cc_eff2_xbin_,_cc_eff2_xmin_,_cc_eff2_xmax_));
						}
					}else if(_hcuts_[ecut_idx] == _event_ && _top_[top_idx] != _mnone_ && ((_species_[cart[4]]==_pro_ && _top_[top_idx]!=_mpro_) || (_species_[cart[4]]==_pip_ && _top_[top_idx]!=_mpip_) || (_species_[cart[4]]==_pim_ && _top_[top_idx]!=_mpim_))){
						//std::cout<<"K Eff1 Idx: " <<_Kinematic_Eff1_hist.size() <<" " <<plot_4d_1.size() <<" " <<plot_3d_1.size() <<" " <<plot_2d_1.size() <<" " <<plot_1d_1.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
						sprintf(hname,"Kinematic_eff_%s_%s_%s_%s_%s",_species_[cart[4]],_hcuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx]);
						plot_1d_1.push_back(new TH2F(hname,hname,_cc_eff1_xbin_,_cc_eff1_xmin_,_cc_eff1_xmax_,_cc_eff1_ybin_,_cc_eff1_ymin_,_cc_eff1_ymax_));
						if(_sector_[sec_idx]==_sec_all_){
							//std::cout<<"\tK Eff2 Idx: " <<_Kinematic_Eff2_hist.size() <<" " <<plot_3d_2.size() <<" " <<plot_2d_2.size() <<" " <<plot_1d_2.size() <<"  " <<_ecuts_[ecut_idx] <<"\n";
							sprintf(hname,"Kinematic_eff_sec_yields_%s_%s_%s_%s_%s",_species_[cart[4]],_hcuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx]);
							plot_1d_2.push_back(new TH1F(hname,hname,_cc_eff2_xbin_,_cc_eff2_xmin_,_cc_eff2_xmax_));
						}
					}
				}
			}
		}
		if(cart[0] == space_dims[0]-1){
			if(plot_1d_1.size()>0){
				plot_2d_1.push_back(plot_1d_1);
				plot_1d_1.clear();
			}
			if(plot_1d_2.size()>0){
				plot_2d_2.push_back(plot_1d_2);
				plot_1d_2.clear();
			}
			if(cart[1] == space_dims[1]-1){
				if(plot_2d_1.size()>0){
					plot_3d_1.push_back(plot_2d_1);
					plot_2d_1.clear();
				}
				if(cart[2] == space_dims[2]-1){
					if(plot_3d_1.size()>0){
						plot_4d_1.push_back(plot_3d_1);
						plot_3d_1.clear();
					}
					if(plot_2d_2.size()>0){
						plot_3d_2.push_back(plot_2d_2);
						plot_2d_2.clear();
					}
					if(cart[3]==space_dims[3]-1){
						if(plot_3d_2.size()>0){
							_Kinematic_Eff2_hist.push_back(plot_3d_2);
							plot_3d_2.clear();
						}
						if(plot_4d_1.size()>0){
							_Kinematic_Eff1_hist.push_back(plot_4d_1);
							plot_4d_1.clear();
						}
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::Kinematic_Eff_idx(int which_, const char* species_, const char * ecut_, const char* cut_, const char* sector_, const char* top_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	if(!flags_->Flags::Plot_Kin_Eff(fun::species_idx(species_))){
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		if(which_==1){
			idx.push_back(-1);
		}
		return idx;
	}
	if(species_ == _pro_ && ecut_==_event_ && top_ == _mpro_){
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		if(which_==1){
			idx.push_back(-1);
		}
		return idx;
	}
	if(species_ == _pip_ && ecut_==_event_ && top_ == _mpip_){
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		if(which_==1){
			idx.push_back(-1);
		}
		return idx;
	}
	if(species_ == _pim_ && ecut_==_event_ && top_ == _mpim_){
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		if(which_==1){
			idx.push_back(-1);
		}
		return idx;
	}
	

	idx.push_back(fun::species_idx(species_)+fun::species_offset(species_,_kin_eff_cut_,flags_));
	if(species_==_ele_){
		if(fun::ecut_perform(ecut_,flags_)){
			if(which_==2){//Yield by sector
				if(ecut_==_none_ && cut_==_no_cut_ && top_==_mnone_){
					idx.push_back(0);
					idx.push_back(0);
					idx.push_back(0);
				}else if(cut_!=_no_cut_){
					if(ecut_ != _event_ && top_ == _mnone_){
						idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
						idx.push_back(fun::cut_idx(cut_));
						idx.push_back(0);
					}else if(ecut_ == _event_ && top_ != _mnone_){
						idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
						idx.push_back(fun::cut_idx(cut_));
						idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
					}else{
						idx.push_back(-1);
						idx.push_back(-1);
						idx.push_back(-1);
					}
				}else{
					idx.push_back(-1);
					idx.push_back(-1);
					idx.push_back(-1);
				}
			}else if(which_==1){//Distribution per sector by p and theta
				if(ecut_==_none_ && cut_==_no_cut_ && top_==_mnone_){
					idx.push_back(0);
					idx.push_back(0);
					idx.push_back(fun::sector_idx(sector_));
					idx.push_back(0);
				}else if(cut_!=_no_cut_){
					if(ecut_ != _event_ && top_ == _mnone_){
						idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
						idx.push_back(fun::cut_idx(cut_));
						idx.push_back(fun::sector_idx(sector_));
						idx.push_back(0);
					}else if(ecut_ == _event_ && top_ != _mnone_){
						idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
						idx.push_back(fun::cut_idx(cut_));
						idx.push_back(fun::sector_idx(sector_));
						idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
					}else{
						idx.push_back(-1);
						idx.push_back(-1);
						idx.push_back(-1);
						idx.push_back(-1);
						//idx.push_back(-1);
					}
				}else{
					idx.push_back(-1);
					idx.push_back(-1);
					idx.push_back(-1);
					idx.push_back(-1);
					//idx.push_back(-1);
				}
			}
		}else{ 
			idx.push_back(-1);
			idx.push_back(-1);
			idx.push_back(-1);
			idx.push_back(-1);
			//idx.push_back(-1);
		}
	}else{
		if(fun::hcut_perform(species_,ecut_,flags_)){
			if(which_==2){//Yield by sector
				if(ecut_==_none_ && cut_==_no_cut_ && top_==_mnone_){
					idx.push_back(0);
					idx.push_back(0);
					idx.push_back(0);
				}else if(cut_!=_no_cut_){
					if(ecut_ != _event_ && top_ == _mnone_){
						idx.push_back(fun::hcut_idx(ecut_)+fun::hcut_offset(species_,ecut_,flags_));
						idx.push_back(fun::cut_idx(cut_));
						idx.push_back(0);
					}else if(ecut_ == _event_ && top_ != _mnone_){
						idx.push_back(fun::hcut_idx(ecut_)+fun::hcut_offset(species_,ecut_,flags_));
						idx.push_back(fun::cut_idx(cut_));
						idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
					}else{
						idx.push_back(-1);
						idx.push_back(-1);
						idx.push_back(-1);
					}
				}else{
					idx.push_back(-1);
					idx.push_back(-1);
					idx.push_back(-1);
				}
			}else if(which_==1){//Distribution per sector by p and theta
				if(ecut_==_none_ && cut_==_no_cut_ && top_==_mnone_){
					idx.push_back(0);
					idx.push_back(0);
					idx.push_back(fun::sector_idx(sector_));
					idx.push_back(0);
				}else if(cut_!=_no_cut_){
					if(ecut_ != _event_ && top_ == _mnone_){
						idx.push_back(fun::hcut_idx(ecut_)+fun::hcut_offset(species_,ecut_,flags_));
						idx.push_back(fun::cut_idx(cut_));
						idx.push_back(fun::sector_idx(sector_));
						idx.push_back(0);
					}else if(ecut_ == _event_ && top_ != _mnone_){
						idx.push_back(fun::hcut_idx(ecut_)+fun::hcut_offset(species_,ecut_,flags_));
						idx.push_back(fun::cut_idx(cut_));
						idx.push_back(fun::sector_idx(sector_));
						idx.push_back(fun::top_idx(top_)+fun::top_offset(top_,flags_));
					}else{
						idx.push_back(-1);
						idx.push_back(-1);
						idx.push_back(-1);
						idx.push_back(-1);
					}
				}else{
					idx.push_back(-1);
					idx.push_back(-1);
					idx.push_back(-1);
					idx.push_back(-1);
				}
			}
		}else{ 
			idx.push_back(-1);
			idx.push_back(-1);
			idx.push_back(-1);
			idx.push_back(-1);
		}
	}
	//if(ecut_=_delta_cut_){
//		std::cout<<"Doing a Delta cut Kineff\n";
//		fun::print_vector_idx(idx);
//	}
	return idx;
	//Problem at 	Filling Kinematic Eff Histogram1 for pro event anti all mpim
	//				Index: 1 5 1 6 2 
}

void Histogram::Kinematic_Eff_Fill(float p_, float theta_, float weight_, const char* species_, const char* ecut_, const char* cut_, const char* sector_, const char* top_, std::shared_ptr<Flags> flags_){
	bool cont = false;
	for(int i=0; i<4; i++){
		if(flags_->Flags::Plot_Kin_Eff(i)){cont=true;}
	}
	if(!cont){return;}
	//if(!flags_->Flags::Plot_Kin_Eff(0) && !flags_->Flags::Plot_Kin_Eff(1) && !flags_->Flags::Plot_Kin_Eff(2) && !flags_->Flags::Plot_Kin_Eff(3)){ return;}
	std::vector<int> idx = Histogram::Kinematic_Eff_idx(1,species_,ecut_,cut_,sector_,top_,flags_);
	//fun::print_vector_idx(idx1);
	if(Histogram::OK_Idx(idx)){
		//std::cout<<"Part 1 Filling Kinematic Eff Histogram1 for "<<species_ <<" " <<ecut_ <<" " <<cut_ <<" " <<sector_ <<" " <<top_ <<"\n";
		//fun::print_vector_idx(idx);
	}
	//if(species_ == _pro_ && ecut_ == _delta_cut_){
//			std::cout<<"Part 1 Filling Kinematic Eff Histogram1 for "<<species_ <<" " <<ecut_ <<" " <<cut_ <<" " <<sector_ <<" " <<top_ <<"\n";
//			fun::print_vector_idx(idx);
//		}
	if(Histogram::OK_Idx(idx)){
		//if(species_ == _pro_ && ecut_ == _event_ && cut_ == _cut_applied_){
		//	std::cout<<"Part 2 Filling Kinematic Eff Histogram1 for "<<species_ <<" " <<ecut_ <<" " <<cut_ <<" " <<sector_ <<" " <<top_ <<"\n";
		//	std::cout<<"\tp:"<< p_ <<" theta:" <<theta_ <<" weight:" <<weight_ <<"\n";
		//	fun::print_vector_idx(idx);
		//}
		_Kinematic_Eff1_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,theta_,weight_);
		idx.clear();
		idx = Histogram::Kinematic_Eff_idx(1,species_,ecut_,cut_,_sec_all_,top_,flags_);
		//std::cout<<"Filling Kinematic Eff Histogram1 for " <<species_ <<" "<<ecut_ <<" " <<cut_ <<" " <<_sec_all_ <<" " <<top_ <<"\n";
		//fun::print_vector_idx(idx);
		if(Histogram::OK_Idx(idx)){
			_Kinematic_Eff1_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,theta_,weight_);
		}
	}
	idx.clear();
	idx = Histogram::Kinematic_Eff_idx(2,species_,ecut_,cut_,sector_,top_,flags_);
	//std::cout<<"Filling Kinematic Eff Histogram2 for " <<species_ <<" "<<ecut_ <<" " <<cut_ <<" " <<_sec_all_ <<" " <<top_ <<"\n";
	//fun::print_vector_idx(idx);
	//fun::print_vector_idx(idx2);
	if(Histogram::OK_Idx(idx)){
		_Kinematic_Eff2_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(fun::sector_idx(sector_)+1,weight_);
	}
}

void Histogram::Kinematic_Eff_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_Kin_Eff(0) || flags_->Flags::Plot_Kin_Eff(1) || flags_->Flags::Plot_Kin_Eff(2) || flags_->Flags::Plot_Kin_Eff(3)){
		std::cout<<"Writing Kinematic Efficiency\n";
		char dirname[100];
		TDirectory* dir_kin_eff = _RootOutputFile->mkdir("Kinematic Efficiency");
		TDirectory* dir_kin_eff_sub[std::distance(std::begin(_species_), std::end(_species_))][std::distance(std::begin(_ecuts_), std::end(_ecuts_))+1][std::distance(std::begin(_cut_), std::end(_cut_))+1][std::distance(std::begin(_top_), std::end(_top_))+1];
		for(int h=0; h<std::distance(std::begin(_species_), std::end(_species_)); h++){
			//std::cout<<"making dir at " <<h <<" 0 0 0\n";
			if(flags_->Flags::Plot_Kin_Eff(h)){
				sprintf(dirname,"kinematic_eff_%s",_species_[h]);
				dir_kin_eff_sub[h][0][0][0] = dir_kin_eff->mkdir(dirname);
				if(_species_[h]==_ele_){
					for(int i=0; i<std::distance(std::begin(_ecuts_), std::end(_ecuts_)); i++){
						if(fun::ecut_perform(_ecuts_[i],flags_)){
							//std::cout<<"making dir at "<<h <<" " <<i <<" " <<0 <<" " <<0 <<" putting in " <<h <<" 0 0 0\n";
							sprintf(dirname,"kinematic_eff_%s_%s",_species_[h],_ecuts_[i]);
							dir_kin_eff_sub[h][i+1][0][0] = dir_kin_eff_sub[h][0][0][0]->mkdir(dirname);
							if(_ecuts_[i]==_event_){
								for(int j=0; j<std::distance(std::begin(_top_), std::end(_top_)); j++){
									if(_top_[j]!=_mnone_ && fun::top_perform(_top_[j],flags_)){
										//std::cout<<"\tmaking dir at "<<h <<" " <<i <<" " <<0 <<" " <<j+1 <<" in " <<h <<" " <<i <<" " <<0 <<" " <<0 <<"\n";
										sprintf(dirname,"kinematic_eff_%s_%s_%s",_species_[h],_ecuts_[i],_top_[j]);
										dir_kin_eff_sub[h][i+1][0][j+1] = dir_kin_eff_sub[h][i+1][0][0]->mkdir(dirname);
										for(int k=0; k<std::distance(std::begin(_cut_), std::end(_cut_)); k++){
											if(_cut_[k]!=_no_cut_){
												//std::cout<<"\t\tmaking dir at "<<h <<" " <<i <<" " <<k+1 <<" " <<j+1 <<" in " <<h <<" " <<i <<" " <<0 <<" " <<j+1 <<"\n";
												sprintf(dirname,"kinematic_eff_%s_%s_%s_%s",_species_[h],_ecuts_[i],_top_[j],_cut_[k]);
												dir_kin_eff_sub[h][i+1][k+1][j+1] = dir_kin_eff_sub[h][i+1][0][j+1]->mkdir(dirname);
											}
										}
									}
								}
							}else{
								for(int k=0; k<std::distance(std::begin(_cut_), std::end(_cut_)); k++){
									if(_cut_[k]!=_no_cut_ && _ecuts_[i]!=_none_){
										//std::cout<<"\tmaking dir at "<<h <<" " <<i <<" " <<k+1 <<" " <<0 <<" in " <<h <<" " <<i <<" " <<0 <<" " <<0 <<"\n";
										sprintf(dirname,"kinematic_eff_%s_%s_%s",_species_[h],_ecuts_[i],_cut_[k]);
										dir_kin_eff_sub[h][i+1][k+1][0] = dir_kin_eff_sub[h][i+1][0][0]->mkdir(dirname);
									}
								}
							}
						}
					}
				}else{
					for(int i=0; i<std::distance(std::begin(_hcuts_), std::end(_hcuts_)); i++){
						if(fun::hcut_perform(_species_[h],_hcuts_[i],flags_)){
							//std::cout<<"making dir at "<<h <<" " <<i <<" " <<0 <<" " <<0 <<" putting in " <<h <<" 0 0 0\n";
							sprintf(dirname,"kinematic_eff_%s_%s",_species_[h],_hcuts_[i]);
							dir_kin_eff_sub[h][i+1][0][0] = dir_kin_eff_sub[h][0][0][0]->mkdir(dirname);
							if(_hcuts_[i]==_event_){
								for(int j=0; j<std::distance(std::begin(_top_), std::end(_top_)); j++){
									if(_top_[j]!=_mnone_ && fun::top_perform(_top_[j],flags_)){
										//if((_species_[h] == _pro_ && _top_[j]!=_mpro_) || (_species_[h] == _pip_ && _top_[j]!=_mpip_) || (_species_[h] == _pim_ && _top_[j]!=_mpim_)){
										//std::cout<<"\tmaking dir at "<<h <<" " <<i <<" " <<0 <<" " <<j+1 <<" in " <<h <<" " <<i <<" " <<0 <<" " <<0 <<"\n";
										sprintf(dirname,"kinematic_eff_%s_%s_%s",_species_[h],_hcuts_[i],_top_[j]);
										dir_kin_eff_sub[h][i+1][0][j+1] = dir_kin_eff_sub[h][i+1][0][0]->mkdir(dirname);
										for(int k=0; k<std::distance(std::begin(_cut_), std::end(_cut_)); k++){
											if(_cut_[k]!=_no_cut_){
												//std::cout<<"\t\tmaking dir at "<<h <<" " <<i <<" " <<k+1 <<" " <<j+1 <<" in " <<h <<" " <<i <<" " <<0 <<" " <<j+1 <<"\n";
												sprintf(dirname,"kinematic_eff_%s_%s_%s_%s",_species_[h],_hcuts_[i],_top_[j],_cut_[k]);
												dir_kin_eff_sub[h][i+1][k+1][j+1] = dir_kin_eff_sub[h][i+1][0][j+1]->mkdir(dirname);
											}
										}
										//}
									}
								}
							}else{
								for(int k=0; k<std::distance(std::begin(_cut_), std::end(_cut_)); k++){
									if(_cut_[k]!=_no_cut_ && _hcuts_[i]!=_none_){
										sprintf(dirname,"kinematic_eff_%s_%s_%s",_species_[h],_hcuts_[i],_cut_[k]);
										dir_kin_eff_sub[h][i+1][k+1][0] = dir_kin_eff_sub[h][i+1][0][0]->mkdir(dirname);
										//std::cout<<"\tmaking dir at "<<h <<" " <<i <<" " <<k+1 <<" " <<0 <<" in " <<h <<" " <<i <<" " <<0 <<" " <<0 <<"\n";
										
									}
								}
							}
						}
					}
				}
			}
		}
		std::vector<long> space_dims(5);
		space_dims[4] = std::distance(std::begin(_species_), std::end(_species_));
		space_dims[3] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));
		space_dims[1] = std::distance(std::begin(_sector_),std::end(_sector_));
		space_dims[0] = std::distance(std::begin(_top_), std::end(_top_));
		char hname[100];
		std::vector<int> idx1;
		std::vector<int> idx2;
		int pcut_idx;
		int cut_idx;
		int sec_idx;
		int top_idx;
		CartesianGenerator cart(space_dims);
		while(cart.GetNextCombination()){
			pcut_idx = cart[3];
			cut_idx = cart[2];
			sec_idx = cart[1];
			top_idx = cart[0];
			/*if(_species_[cart[4]]==_ele_){
				std::cout<<"Looking at 1" <<_species_[cart[4]] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<"\n";
			}else if(pcut_idx < std::distance(std::begin(_hcuts_), std::end(_hcuts_))){
				std::cout<<"Looking at 1" <<_species_[cart[4]] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<"\n";
			}*/
			if(flags_->Flags::Plot_Kin_Eff(cart[4])){
				if(_species_[cart[4]]==_ele_){	
					idx1 = Histogram::Kinematic_Eff_idx(1,_species_[cart[4]],_ecuts_[pcut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx],flags_);
					idx2 = Histogram::Kinematic_Eff_idx(2,_species_[cart[4]],_ecuts_[pcut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx],flags_);
					if(fun::ecut_perform(_ecuts_[pcut_idx],flags_)){
						if(_ecuts_[pcut_idx] == _none_ && _cut_[cut_idx] == _no_cut_ && _top_[top_idx]==_none_){
						//std::cout<<"Writing In Directory 2| " <<fun::species_idx(_species_[cart[4]]) <<" " <<fun::ecut_idx(_ecuts_[pcut_idx]) <<" 0 0\n"; 
						dir_kin_eff_sub[fun::species_idx(_species_[cart[4]])][fun::ecut_idx(_ecuts_[pcut_idx])+1][0][0]->cd();
							if(Histogram::OK_Idx(idx1)){
								//fun::print_vector_idx(idx1);
								_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetXTitle("Momentum (GeV)");
								_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetYTitle("Theta (Degrees)");
								_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->Write();
							}
							if(_sector_[sec_idx]==_sec_all_){
								if(Histogram::OK_Idx(idx2)){
									//fun::print_vector_idx(idx2);
									_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetXTitle("Sector");
									_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetYTitle("Yield");
									_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->Write();
								}
							}
						}else if(_ecuts_[pcut_idx] != _none_ && _cut_[cut_idx] != _no_cut_){
							if(_ecuts_[pcut_idx] != _event_ && _top_[top_idx] == _mnone_){
								//std::cout<<"Writing In Directory 3| " <<fun::species_idx(_species_[cart[4]]) <<" " <<fun::ecut_idx(_ecuts_[pcut_idx]) <<" " <<fun::cut_idx(_cut_[cut_idx])+1 <<" 0\n";
								dir_kin_eff_sub[fun::species_idx(_species_[cart[4]])][fun::ecut_idx(_ecuts_[pcut_idx])+1][fun::cut_idx(_cut_[cut_idx])+1][0]->cd();
								if(Histogram::OK_Idx(idx1)){
									//fun::print_vector_idx(idx1);
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetXTitle("Momentum (GeV)");
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetYTitle("Theta (Degrees)");
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->Write();
								}
								if(_sector_[sec_idx]==_sec_all_){
									if(Histogram::OK_Idx(idx2)){
										//fun::print_vector_idx(idx2);
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetXTitle("Sector");
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetYTitle("Yield");
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->Write();
									}
								}
							}else if(_ecuts_[pcut_idx] == _event_ && _top_[top_idx] != _mnone_){
								std::cout<<"Writing In Directory 4| " <<fun::species_idx(_species_[cart[4]]) <<" " <<fun::ecut_idx(_ecuts_[pcut_idx]) <<" " <<fun::cut_idx(_cut_[cut_idx])+1 <<" " <<fun::top_idx(_top_[top_idx])+1 <<"\n";
								dir_kin_eff_sub[fun::species_idx(_species_[cart[4]])][fun::ecut_idx(_ecuts_[pcut_idx])+1][fun::cut_idx(_cut_[cut_idx])+1][fun::top_idx(_top_[top_idx])+1]->cd();
								if(Histogram::OK_Idx(idx1)){
									//fun::print_vector_idx(idx1);
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetXTitle("Momentum (GeV)");
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetYTitle("Theta (Degrees)");
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->Write();
								}
								if(_sector_[sec_idx]==_sec_all_){
									if(Histogram::OK_Idx(idx2)){
										//fun::print_vector_idx(idx1);
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetXTitle("Sector");
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetYTitle("Yield");
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->Write();
									}
								}
							}
						}
					}
				}else if(pcut_idx < std::distance(std::begin(_hcuts_), std::end(_hcuts_))){
					idx1 = Histogram::Kinematic_Eff_idx(1,_species_[cart[4]],_hcuts_[pcut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx],flags_);
					idx2 = Histogram::Kinematic_Eff_idx(2,_species_[cart[4]],_hcuts_[pcut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx],flags_);
					if(fun::hcut_perform(_species_[cart[4]],_hcuts_[pcut_idx],flags_)){
						if(_hcuts_[pcut_idx] == _none_ && _cut_[cut_idx] == _no_cut_ && _top_[top_idx]==_none_){
						//std::cout<<"Writing In Directory 5| " <<fun::species_idx(_species_[cart[4]]) <<" " <<fun::hcut_idx(_hcuts_[pcut_idx]) <<" 0 0\n"; 
						dir_kin_eff_sub[fun::species_idx(_species_[cart[4]])][fun::hcut_idx(_hcuts_[pcut_idx])+1][0][0]->cd();
							if(Histogram::OK_Idx(idx1)){
								//fun::print_vector_idx(idx1);
								_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetXTitle("Momentum (GeV)");
								_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetYTitle("Theta (Degrees)");
								_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->Write();
							}
							if(_sector_[sec_idx]==_sec_all_){
								if(Histogram::OK_Idx(idx2)){
									//fun::print_vector_idx(idx2);
									_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetXTitle("Sector");
									_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetYTitle("Yield");
									_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->Write();
								}
							}
						}else if(_hcuts_[pcut_idx] != _none_ && _cut_[cut_idx] != _no_cut_){
							if(_hcuts_[pcut_idx] != _event_ && _top_[top_idx] == _mnone_){
								//std::cout<<"Writing In Directory 6| " <<fun::species_idx(_species_[cart[4]]) <<" " <<fun::hcut_idx(_hcuts_[pcut_idx]) <<" " <<fun::cut_idx(_cut_[cut_idx])+1 <<" 0\n";
								dir_kin_eff_sub[fun::species_idx(_species_[cart[4]])][fun::hcut_idx(_hcuts_[pcut_idx])+1][fun::cut_idx(_cut_[cut_idx])+1][0]->cd();
								if(Histogram::OK_Idx(idx1)){
									//fun::print_vector_idx(idx1);
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetXTitle("Momentum (GeV)");
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetYTitle("Theta (Degrees)");
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->Write();
								}
								if(_sector_[sec_idx]==_sec_all_){
									if(Histogram::OK_Idx(idx2)){
										//fun::print_vector_idx(idx2);
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetXTitle("Sector");
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetYTitle("Yield");
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->Write();
									}
								}
							}else if(_hcuts_[pcut_idx] == _event_ && _top_[top_idx] != _mnone_ && fun::top_perform(_top_[top_idx],flags_)){
								//std::cout<<"Writing In Directory 7| " <<fun::species_idx(_species_[cart[4]]) <<" " <<fun::hcut_idx(_hcuts_[pcut_idx]) <<" " <<fun::cut_idx(_cut_[cut_idx])+1 <<" " <<fun::top_idx(_top_[top_idx])+1 <<"\n";
								//if((_top_[top_idx]==_mpro_ && _species_[cart[4]]!=_pro_) || (_top_[top_idx]==_mpip_ && _species_[cart[4]]!=_pip_) || (_top_[top_idx]==_mpim_ && _species_[cart[4]]!=_pim_) ){
								dir_kin_eff_sub[fun::species_idx(_species_[cart[4]])][fun::hcut_idx(_hcuts_[pcut_idx])+1][fun::cut_idx(_cut_[cut_idx])+1][fun::top_idx(_top_[top_idx])+1]->cd();
								if(Histogram::OK_Idx(idx1)){
									//fun::print_vector_idx(idx1);
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetXTitle("Momentum (GeV)");
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->SetYTitle("Theta (Degrees)");
									_Kinematic_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]][idx1[4]]->Write();
								}
								if(_sector_[sec_idx]==_sec_all_){
									if(Histogram::OK_Idx(idx2)){
										//fun::print_vector_idx(idx1);
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetXTitle("Sector");
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->SetYTitle("Yield");
										_Kinematic_Eff2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]]->Write();
									}
								}
								//}
							}
						}
					}
				}
			}
			idx1.clear();
			idx2.clear();
		}
	}
}

//*-------------------------------End Kinematic Eff Plot------------------------------*
//*-------------------------------Start SC Eff Plot----------------------------*
void Histogram::SC_Eff_Make(std::shared_ptr<Flags> flags_){//Go by 
	if(!flags_->Flags::Plot_SC_Eff()){
		return;
	}
	std::cout<<"\tMaking SC Efficiency Histograms\n";
	std::vector<long> space_dims(4);
	space_dims[3] = std::distance(std::begin(_species_),std::end(_species_));
	space_dims[2] = std::distance(std::begin(_sector_),std::end(_sector_));
	space_dims[1] = std::distance(std::begin(_ecuts_),std::end(_ecuts_));
	space_dims[0] = std::distance(std::begin(_cut_),std::end(_cut_));//Don't include "no cut"
	

	TH1F_ptr_1d plot_1d;
	TH1F_ptr_2d plot_2d;
	TH1F_ptr_3d plot_3d;
	char hname[100];
	CartesianGenerator cart(space_dims);
	int species_idx = -1;
	int sec_idx = -1;
	int pcut_idx = -1;
	int cut_idx = -1;
	while(cart.GetNextCombination()){
		species_idx = cart[3];
		sec_idx = cart[2];
		pcut_idx = cart[1];
		cut_idx = cart[0];
		if(_species_[species_idx]==_ele_){
			if(fun::ecut_perform(_ecuts_[pcut_idx],flags_) && ((_ecuts_[pcut_idx]==_none_ && _cut_[cut_idx]==_no_cut_) || ((_ecuts_[pcut_idx]!=_none_ && _cut_[cut_idx]!=_no_cut_)))){
				sprintf(hname,"SC_Paddle_Yield_%s_%s_%s_%s",_species_[species_idx],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx]);
				//std::cout<<"hist name: " <<hname <<" |for species: " <<_species_[species_idx] <<" sector: " <<_sector_[sec_idx] <<" ecut: " <<_ecuts_[pcut_idx] <<" cut: " <<_cut_[cut_idx] <<" => idx: " <<_SC_Eff_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
				plot_1d.push_back(new TH1F(hname,hname,_sc_eff_xbin_,_sc_eff_xmin_,_sc_eff_xmax_));
			}
		}else if(pcut_idx< std::distance(std::begin(_hcuts_), std::end(_hcuts_)) && ((_hcuts_[pcut_idx]==_none_ && _cut_[cut_idx]==_no_cut_) || ((_hcuts_[pcut_idx]!=_none_ && _cut_[cut_idx]!=_no_cut_)))){
			if(fun::hcut_perform(_species_[species_idx],_hcuts_[pcut_idx],flags_)){
				sprintf(hname,"SC_Paddle_Yield_%s_%s_%s_%s",_species_[species_idx],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx]);
				//std::cout<<"hist name: " <<hname <<" |for species: " <<_species_[species_idx] <<" sector: " <<_sector_[sec_idx] <<" hcut: " <<_hcuts_[pcut_idx] <<" cut: " <<_cut_[cut_idx] <<" => idx: " <<_SC_Eff_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
				plot_1d.push_back(new TH1F(hname,hname,_sc_eff_xbin_,_sc_eff_xmin_,_sc_eff_xmax_));
			}
		}
		if(cut_idx==space_dims[0]-1){
			if(plot_1d.size()>0){
				plot_2d.push_back(plot_1d);
				plot_1d.clear();
			}
			if(pcut_idx==space_dims[1]-1){
				if(plot_2d.size()>0){
					plot_3d.push_back(plot_2d);
					plot_2d.clear();
				}
				if(sec_idx==space_dims[2]-1){
					if(plot_3d.size()>0){
						_SC_Eff_hist.push_back(plot_3d);
						plot_3d.clear();
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::SC_Eff_idx(const char * species_, const char * sector_, const char * pcut_, const char* cut_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	if(!flags_->Flags::Plot_SC_Eff()){
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		return idx;
	}
	idx.push_back(fun::species_idx(species_));
	idx.push_back(fun::sector_idx(sector_));
	if(fun::pcut_perform(species_,pcut_,flags_) && ((pcut_==_none_ && cut_==_no_cut_)||(pcut_!=_none_ && cut_!=_no_cut_))){
		if(species_==_ele_){
			idx.push_back(fun::ecut_idx(pcut_)+fun::ecut_offset(pcut_,flags_));
		}else{
			idx.push_back(fun::hcut_idx(pcut_)+fun::hcut_offset(species_,pcut_,flags_));
		}
	}else{
		idx.push_back(-1);
	}
	if(pcut_==_none_){
		idx.push_back(0);
	}else{
		idx.push_back(fun::cut_idx(cut_));
	}
	
	//std::cout<<"Looking for index of " <<species_ <<" " <<sector_ <<" " <<pcut_ <<" => " <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<"\n";
	return idx;
}
void Histogram::SC_Eff_Fill(int sc_pd_, float weight_, const char* species_, const char* pcut_, const char* cut_, const char * sector_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_SC_Eff()){
		return;
	}
	//std::cout<<"Filling for " <<species_ <<" " <<sector_ <<" " <<pcut_ <<" " <<cut_ <<"\n";
	std::vector<int> idx=SC_Eff_idx(species_,sector_,pcut_,cut_,flags_);
	//fun::print_vector_idx(idx);
	if(Histogram::OK_Idx(idx)){
		_SC_Eff_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(sc_pd_,weight_);
	}
	idx.clear();
	idx=SC_Eff_idx(species_,_sec_all_,pcut_,cut_,flags_);
	//fun::print_vector_idx(idx);
	if(Histogram::OK_Idx(idx)){
		_SC_Eff_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(sc_pd_,weight_);
	}
}
void Histogram::SC_Eff_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_SC_Eff()){
		return;
	}
	std::cout<<"Writing SC Efficiency Histograms\n";
	std::vector<long> space_dims(4);
	space_dims[3] = std::distance(std::begin(_species_),std::end(_species_));
	space_dims[2] = std::distance(std::begin(_sector_),std::end(_sector_));
	space_dims[1] = std::distance(std::begin(_ecuts_),std::end(_ecuts_));
	space_dims[0] = std::distance(std::begin(_cut_),std::end(_cut_))-1;//Don't include "no cut";
	TDirectory* Dir_SC = _RootOutputFile->mkdir("SC_Eff");
	TDirectory* Dir_SC_sub[space_dims[3]][space_dims[2]+1][space_dims[1]+1];
	CartesianGenerator cart(space_dims);
	char dirname[100];
	std::vector<int> idx; 
	for(int i=0; i<space_dims[3]; i++){//species
		sprintf(dirname,"SC_Eff_%s",_species_[i]);
		Dir_SC_sub[i][0][0] = Dir_SC->mkdir(dirname);
		for(int j=0; j<space_dims[2]; j++){//Sector
			sprintf(dirname,"SC_Eff_%s_%s",_species_[i],_sector_[j]);
			Dir_SC_sub[i][j+1][0] = Dir_SC_sub[i][0][0]->mkdir(dirname);
			for(int k=0; k<space_dims[1]; k++){//Pcut
				if(_species_[i]==_ele_ && fun::ecut_perform(_ecuts_[k],flags_)){
					sprintf(dirname,"SC_Eff_%s_%s_%s",_species_[i],_sector_[j],_ecuts_[k]);
					Dir_SC_sub[i][j+1][k+1] = Dir_SC_sub[i][j+1][0]->mkdir(dirname);
				}else if((k < std::distance(std::begin(_hcuts_),std::end(_hcuts_))) && fun::hcut_perform(_species_[i],_hcuts_[k],flags_)){
					sprintf(dirname,"SC_Eff_%s_%s_%s",_species_[i],_sector_[j],_hcuts_[k]);
					Dir_SC_sub[i][j+1][k+1] = Dir_SC_sub[i][j+1][0]->mkdir(dirname);
				}
			}
		}
	}
	while(cart.GetNextCombination()){
		if(_species_[cart[3]]==_ele_){
			idx = SC_Eff_idx(_species_[cart[3]],_sector_[cart[2]],_ecuts_[cart[1]],_cut_[cart[0]],flags_);
		}else{
			idx = SC_Eff_idx(_species_[cart[3]],_sector_[cart[2]],_hcuts_[cart[1]],_cut_[cart[0]],flags_);
		}
		if(Histogram::OK_Idx(idx)){
			Dir_SC_sub[cart[3]][cart[2]+1][cart[1]+1]->cd();
			_SC_Eff_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetXTitle("SC Paddle");
			_SC_Eff_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetYTitle("Yield");
			_SC_Eff_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Write();
		}
		idx.clear();
	} 
}
//*-------------------------------End SC Eff Plot------------------------------*
//*------------------------------- Start Friend Plot --------------------------*
std::vector<int> Histogram::Friend_Bin_Sizes(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Make_Friend()){
		std::vector<int> bin_sizes(7);
		bin_sizes[0] = Histogram::W_bins();
		bin_sizes[1] = std::distance(std::begin(_Q2_bins_), std::end(_Q2_bins_));
		bin_sizes[2] = _MM_bins_+2*_MM_wider_;
		bin_sizes[3] = _MM_bins_+2*_MM_wider_;
		bin_sizes[4] = _theta_bins_;
		bin_sizes[5] = _alpha_bins_;
		bin_sizes[6] = _phi_bins_;
		return bin_sizes;
	}
}

void Histogram::Friend_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Make_Friend()){
		std::cout<<"Making Friend Histograms\n";
		char hname[100];
		char hname2[100];
		Int_t bins[7];
		Int_t bins_5d[5];
		std::vector<int> bins_vec = Histogram::Friend_Bin_Sizes(flags_);//A Bit redundant, but I believe it needs to be an Int array rather than a vector
		for(int k=0; k<7; k++){
			bins[k] = bins_vec[k];
			if(k>1){
				bins_5d[k-2] = bins_vec[k];
			}
		}
		double xmin[5];
		double xmax[5];
		xmin[2]=_theta_min_;
		xmin[3]=_alpha_min_;
		xmin[4]=_phi_min_;
		xmax[2]=_theta_max_;
		xmax[3]=_alpha_max_;
		xmax[4]=_phi_max_;
		//int _wider_ = 5;
		for(int i=0; i<3; i++){//Variable Sets
			//Double_t xmin[7] = {_W_min_,_Q2_min_,_MM_min_[i],_MM2_min_[i],_theta_min_,_alpha_min_,_phi_min_};
			//Double_t xmax[7] = {_W_max_,_Q2_max_,_MM_max_[i],_MM2_max_[i],_theta_max_,_alpha_max_,_phi_max_};
			//Double_t xmin[5] = {_MM_min_[i],_MM2_min_[i],_theta_min_,_alpha_min_,_phi_min_};
			//Double_t xmax[5] = {_MM_max_[i],_MM2_max_[i],_theta_max_,_alpha_max_,_phi_max_};
			for(int j=0; j<29; j++){//W
				//float MM_res = ((Histogram::MM_max(j,i)-_MM_min_[i])/_MM_bins_);
				//float MM2_res = ((Histogram::MM2_max(j,i)-_MM2_min_[i])/_MM_bins_);
				xmin[0] = _MM_min_[i]-_MM_wider_*((_MM_max_[i][j]-_MM_min_[i])/_MM_bins_);
				xmin[1] = _MM2_min_[i]-_MM_wider_*((_MM2_max_[i][j]-_MM2_min_[i])/_MM_bins_);
				xmax[0] = _MM_max_[i][j]+_MM_wider_*((_MM_max_[i][j]-_MM_min_[i])/_MM_bins_);
				xmax[1] = _MM2_max_[i][j]+_MM_wider_*((_MM2_max_[i][j]-_MM2_min_[i])/_MM_bins_);
				//std::cout<<"W Bin: " <<j <<" var set: " <<i <<"\n";
				for(int l=0; l<5; l++){
					//std::cout<<"xmin[" <<l <<"]: " <<xmin[l] <<"  xmax[" <<l <<"]: " <<xmax[l] <<"\n";
				}
				for(int k=0; k<5; k++){//Q2
					sprintf(hname,"2#pi_off_proton_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));
					_Friend[i][j][k] = new THnSparseD(hname,hname,5,bins_5d,xmin,xmax);
					_Friend[i][j][k]->Sumw2();//Weights as normal
					sprintf(hname,"sys_err_2#pi_off_proton_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));
					_Duo[i][j][k] = new THnSparseD(hname,hname,5,bins_5d,xmin,xmax);
					_Duo[i][j][k]->Sumw2();//Weights as normal
					if(flags_->Flags::Helicity()){
						sprintf(hname,"2#pi_off_proton_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_pos",_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));
						_Friend1[i][j][k] = new THnSparseD(hname,hname,5,bins_5d,xmin,xmax);
						_Friend1[i][j][k]->Sumw2();//Normal Weights
						sprintf(hname,"2#pi_off_proton_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_neg",_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));
						_Friend2[i][j][k] = new THnSparseD(hname,hname,5,bins_5d,xmin,xmax);
						_Friend2[i][j][k]->Sumw2();//Normal Weights
						sprintf(hname,"sys_err_2#pi_off_proton_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_pos",_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));
						_Duo1[i][j][k] = new THnSparseD(hname,hname,5,bins_5d,xmin,xmax);
						_Duo1[i][j][k]->Sumw2();//Normal Weights
						sprintf(hname,"sys_err_2#pi_off_proton_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_neg",_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));
						_Duo2[i][j][k] = new THnSparseD(hname,hname,5,bins_5d,xmin,xmax);
						_Duo2[i][j][k]->Sumw2();//Normal Weights
					}
				
					for(int l=0; l<5; l++){
						if(l>1){
							sprintf(hname,"%s_Friend_Dist_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_friend_pars_[l+2],_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));//_top_[j]);
							sprintf(hname2,"Thrown_%s_Friend_Dist_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_friend_pars_[l+2],_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));//_top_[j]);
						}else if(l==0){
							sprintf(hname,"%s_Friend_Dist_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_MM1_[i],_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));//_top_[j]);
							sprintf(hname2,"Thrown_%s_Friend_Dist_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_MM1_[i],_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));//_top_[j]);
						}else if(l==1){
							sprintf(hname,"%s_Friend_Dist_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_MM2_[i],_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));//_top_[j]);
							sprintf(hname2,"Thrown_%s_Friend_Dist_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_MM2_[i],_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));//_top_[j]);
						}
							
						//used to be [i][j] for topologies
						switch(l){
							/*case 0: 
								_W_Dist[i] = new TH1D(hname,hname,Histogram::W_bins(),_W_min_,_W_max_);
							break;
							case 1: 
								//_Q2_Dist[i] = new TH1D(hname,hname,5,_Q2_min_,_Q2_max_);
								_Q2_Dist[i] = new TH1D(hname,hname,5,_Q2_bins_);
							break;*/
							//case 2:
							case 0: 
								//_MM1_Dist[i][j][k] = new TH1D(hname,hname,_MM_bins_,_MM_min_[i],Histogram::MM_max(j,i));
								_MM1_Dist[i][j][k] = new TH1D(hname,hname,140,_MM_min_[i]-_MM_wider_*((_MM_max_[i][j]-_MM_min_[i])/_MM_bins_),_MM_max_[i][j]+_MM_wider_*((_MM_max_[i][j]-_MM_min_[i])/_MM_bins_));
								if(flags_->Flags::Sim()){
									//_MM1_Dist_thr[i][j][k] = new TH1D(hname2,hname2,_MM_bins_,_MM_min_[i],_MM_max_[i][j]);
									_MM1_Dist_thr[i][j][k] = new TH1D(hname2,hname2,140,_MM_min_[i]-_MM_wider_*((_MM_max_[i][j]-_MM_min_[i])/_MM_bins_),_MM_max_[i][j]+_MM_wider_*((_MM_max_[i][j]-_MM_min_[i])/_MM_bins_));
								}
							break;
							//case 3:
							case 1: 
								//_MM2_Dist[i][j][k] = new TH1D(hname,hname,_MM_bins_,_MM2_min_[i],Histogram::MM2_max(j,i));
								_MM2_Dist[i][j][k] = new TH1D(hname,hname,140,_MM2_min_[i]-_MM_wider_*((_MM2_max_[i][j]-_MM2_min_[i])/_MM_bins_),_MM2_max_[i][j]+_MM_wider_*((_MM2_max_[i][j]-_MM2_min_[i])/_MM_bins_));
								if(flags_->Flags::Sim()){
									//_MM2_Dist_thr[i][j][k] = new TH1D(hname2,hname2,_MM_bins_,_MM2_min_[i],_MM2_max_[i][j]);
									_MM2_Dist_thr[i][j][k] = new TH1D(hname2,hname2,140,_MM2_min_[i]-_MM_wider_*((_MM2_max_[i][j]-_MM2_min_[i])/_MM_bins_),_MM2_max_[i][j]+_MM_wider_*((_MM2_max_[i][j]-_MM2_min_[i])/_MM_bins_));
								}
							break;
							//case 4: 
							case 2:
								_Theta_Dist[i][j][k] = new TH1D(hname,hname,_theta_bins_,_theta_min_,_theta_max_);
								if(flags_->Flags::Sim()){
									_Theta_Dist_thr[i][j][k] = new TH1D(hname2,hname2,_theta_bins_,_theta_min_,_theta_max_);
								}
							break;
							//case 5:
							case 3: 
								_Alpha_Dist[i][j][k] = new TH1D(hname,hname,_alpha_bins_,_alpha_min_,_alpha_max_);
								if(flags_->Flags::Sim()){
									_Alpha_Dist_thr[i][j][k] = new TH1D(hname2,hname2,_alpha_bins_,_alpha_min_,_alpha_max_);
								}
							break;
							//case 6:
							case 4: 
								_Phi_Dist[i][j][k] = new TH1D(hname,hname,_phi_bins_,_phi_min_,_phi_max_);
								if(flags_->Flags::Sim()){
									_Phi_Dist_thr[i][j][k] = new TH1D(hname2,hname2,_phi_bins_,_phi_min_,_phi_max_);
								}
							break;
						}
					}
					if(flags_->Flags::Sim()){
						sprintf(hname,"Thrown_2#pi_off_proton_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_var_names_[i],Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));
						_Thrown[i][j][k] = new THnSparseD(hname,hname,5,bins_5d,xmin,xmax);
						//_Thrown[i][j][k]->GetAxis(1)->Set(5,_Q2_bins_);
						_Thrown[i][j][k]->Sumw2();//Allow Weights with virtual photon flux, etc. 
						//sprintf(hname,"Scaled_Thrown_2#pi_off_proton_%s",_var_names_[i]);
						//_W_Thrown[i] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
						//_W_Thrown[i]->GetAxis(1)->Set(5,_Q2_bins_);
						//_W_Thrown[i]->Sumw2();//No modification to weights with virtual photon flux
					}
				}
			}/*
			//for(int j=0; j<5; j++){//Topologies
			sprintf(hname,"2#pi_off_proton_%s_%s",_var_names_[i],_mall_);//_top_[j]);
				//_Friend[i][j] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
				//_Friend[i][j]->GetAxis(1)->Set(5,_Q2_bins_);
				//_Friend[i][j]->Sumw2();//Weights as normal
			_Friend[i] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
			_Friend[i]->GetAxis(1)->Set(5,_Q2_bins_);
			_Friend[i]->Sumw2();//Weights as normal
			if(flags_->Flags::Helicity()){
				sprintf(hname,"2#pi_off_proton_%s_%s_pos",_var_names_[i],_mall_);//_top_[j]);
					//_Friend1[i][j] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
					//_Friend1[i][j]->GetAxis(1)->Set(5,_Q2_bins_);
					//_Friend1[i][j]->Sumw2();//Normal Weights
				_Friend1[i] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
				_Friend1[i]->GetAxis(1)->Set(5,_Q2_bins_);
				_Friend1[i]->Sumw2();//Normal Weights
					//sprintf(hname,"Scaled_2#pi_off_proton_%s_%s_pos",_var_names_[i],_top_[j]);
					//_W_Friend1[i][j] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
					//_W_Friend1[i][j]->GetAxis(1)->Set(5,_Q2_bins_);
					//_W_Friend1[i][j]->Sumw2();//Normal weights scaled with virtual photon flux (and cc efficiency for exp)
				sprintf(hname,"2#pi_off_proton_%s_%s_neg",_var_names_[i],_mall_);//_top_[j]);
					//_Friend2[i][j] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
					//_Friend2[i][j]->GetAxis(1)->Set(5,_Q2_bins_);
					//_Friend2[i][j]->Sumw2();//Normal Weights
				_Friend2[i] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
				_Friend2[i]->GetAxis(1)->Set(5,_Q2_bins_);
				_Friend2[i]->Sumw2();//Normal Weights
					//sprintf(hname,"Scaled_2#pi_off_proton_%s_%s_neg",_var_names_[i],_top_[j]);
					//_W_Friend2[i][j] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
					//_W_Friend2[i][j]->GetAxis(1)->Set(5,_Q2_bins_);
					//_W_Friend2[i][j]->Sumw2();//Normal weights scaled with virtual photon flux (and cc efficiency for exp)
			}
			for(int k=0; k<5; k++){
				if(k!=0 && k!=1){
					sprintf(hname,"%s_Friend_Dist_%s_%s_W:%.2f-%.2f_Q2:%.2f-%.2f",_friend_pars_[k],_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));//_top_[j]);
				}else{
					if(k==2){
						sprintf(hname,"%s_Friend_Dist_%s_%s_W:%.2f-%.2f_Q2:%.2f-%.2f",_MM1_[i],_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));//_top_[j]);
					}else if(k==3){
						sprintf(hname,"%s_Friend_Dist_%s_%s_W:%.2f-%.2f_Q2:%.2f-%.2f",_MM2_[i],_var_names_[i],_mall_,Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));//_top_[j]);
					}
				}	
				//used to be [i][j] for topologies
				switch(k){
					case 0: 
						_W_Dist[i] = new TH1D(hname,hname,Histogram::W_bins(),_W_min_,_W_max_);
					break;
					case 1: 
						//_Q2_Dist[i] = new TH1D(hname,hname,5,_Q2_min_,_Q2_max_);
						_Q2_Dist[i] = new TH1D(hname,hname,5,_Q2_bins_);
					break;
					case 2: 
						_MM1_Dist[i] = new TH1D(hname,hname,_MM_bins_,_MM_min_[i],_MM_max_[i]);
					break;
					case 3: 
						_MM2_Dist[i] = new TH1D(hname,hname,_MM_bins_,_MM2_min_[i],_MM2_max_[i]);
					break;
					case 4: 
						_Theta_Dist[i] = new TH1D(hname,hname,_theta_bins_,_theta_min_,_theta_max_);
					break;
					case 5: 
						_Alpha_Dist[i] = new TH1D(hname,hname,_alpha_bins_,_alpha_min_,_alpha_max_);
					break;
					case 6: 
						_Phi_Dist[i] = new TH1D(hname,hname,_phi_bins_,_phi_min_,_phi_max_);
					break;
				}
			}
			if(flags_->Flags::Sim()){
				sprintf(hname,"Thrown_2#pi_off_proton_%s_W:%.2f-%.2f_Q2:%.2f-%.2f",_var_names_[i],Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(k),Histogram::Q2_top(k));
				_Thrown[i] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
				_Thrown[i]->GetAxis(1)->Set(5,_Q2_bins_);
				_Thrown[i]->Sumw2();//Allow Weights with virtual photon flux, etc. 
				//sprintf(hname,"Scaled_Thrown_2#pi_off_proton_%s",_var_names_[i]);
				//_W_Thrown[i] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
				//_W_Thrown[i]->GetAxis(1)->Set(5,_Q2_bins_);
				//_W_Thrown[i]->Sumw2();//No modification to weights with virtual photon flux
			}*/
		} //Making 5d histograms now uggg
	}
}


int Histogram::Friend_W_idx(float W_){
	return Histogram::W_bin(W_); 
}

int Histogram::Friend_Q2_idx(float Q2_){
  int j = -1;
  for(int i = 0; i < std::distance(std::begin(_Q2_bins_), std::end(_Q2_bins_)); i++){//constants.hpp
    if(Q2_ < _Q2_bins_[i+1] && Q2_ >= _Q2_bins_[i]){
      j = i; 
    }
  }
  return j; 
}
/*
int Histogram::Friend_MM_idx(float MM_, int var_, int Wbin_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _MM_bins_; i++){//constants.hpp
    top = _MM_min_[var_] + (i+1)*((_MM_max_[var_][Wbin_]-_MM_min_[var_])/_MM_bins_);//constants.hpp
    bot = top - ((_MM_max_[var_][Wbin_]-_MM_min_[var_])/_MM_bins_); 
    if(MM_ < top && MM_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_MM2_idx(float MM_, int var_, int Wbin_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _MM_bins_; i++){//constants.hpp
    top = _MM2_min_[var_] + (i+1)*((_MM2_max_[var_][Wbin_]-_MM2_min_[var_])/_MM_bins_);//constants.hpp
    bot = top - ((_MM2_max_[var_][Wbin_]-_MM2_min_[var_])/_MM_bins_); 
    if(MM_ < top && MM_ >= bot){
      j = i; 
    }
  }
  return j; 
}*/

int Histogram::Friend_theta_idx(float theta_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _theta_bins_; i++){//constants.hpp
    top = _theta_min_ + (i+1)*((_theta_max_-_theta_min_)/_theta_bins_);//constants.hpp
    bot = top - ((_theta_max_-_theta_min_)/_theta_bins_); 
    if(theta_ < top && theta_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_alpha_idx(float alpha_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _alpha_bins_; i++){//constants.hpp
    top = _alpha_min_ + (i+1)*((_alpha_max_-_alpha_min_)/_alpha_bins_);//constants.hpp
    bot = top - ((_alpha_max_-_alpha_min_)/_alpha_bins_); 
    if(alpha_ < top && alpha_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_phi_idx(float phi_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _phi_bins_; i++){//constants.hpp
    top = _phi_min_ + (i+1)*((_phi_max_-_phi_min_)/_phi_bins_);//constants.hpp
    bot = top - ((_phi_max_-_phi_min_)/_phi_bins_); 
    if(phi_ < top && phi_ >= bot){
      j = i; 
    }
  }
  return j; 
}


int Histogram::Friend_MM_idx(float MM_, int var_, int W_bin_){
  int j = -1;
  float top, bot; 
  float res = ((_MM_max_[var_][W_bin_]-_MM_min_[var_])/_MM_bins_);
  for(int i = 0; i < (_MM_bins_+2*_MM_wider_); i++){//constants.hpp
    bot = _MM_min_[var_] + ((i)- _MM_wider_)*res;
	top = bot + res;
	//top = _MM_min_[var_]-(_MM_wider_-(i+1))*res;//constants.hpp
    //bot = top - res; 
    if(MM_ < top && MM_ >= bot){
      j = i; 
    }
  }
  return j; 
}
int Histogram::Friend_MM2_idx(float MM_, int var_, int W_bin_){
  int j = -1;
  float top, bot; 
  float res = ((_MM2_max_[var_][W_bin_]-_MM2_min_[var_])/_MM_bins_);
  for(int i = 0; i < (_MM_bins_+2*_MM_wider_); i++){//constants.hpp
    bot = _MM2_min_[var_] + ((i)- _MM_wider_)*res;
	top = bot + res;
	//top = _MM2_min_[var_]-(_MM_wider_-(i+1))*res;//constants.hpp
    //bot = top - res; 
    if(MM_ < top && MM_ >= bot){
      j = i; 
    }
  }
  return j; 
}


std::vector<int>  Histogram::Friend_idx( float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_){
	std::vector<int> x(7); 
	int test = 0; 
	bool in_region = false; 
	//x[0] = top; //{pmiss,pipmiss,pimmiss,zeromiss,all}
	x[0] = Histogram::Friend_W_idx(W_);
	x[1] = Histogram::Friend_Q2_idx(Q2_);
	x[2] = Histogram::Friend_MM_idx(MM_,var_,x[0]);
	x[3] = Histogram::Friend_MM2_idx(MM2_,var_,x[0]);
	x[4] = Histogram::Friend_theta_idx(theta_);
	x[5] = Histogram::Friend_alpha_idx(alpha_);
	x[6] = Histogram::Friend_phi_idx(phi_);
	/*std::cout<<std::endl <<"Filling Friend in channel " <<channel <<std::endl; 
	std::cout<<"top = " <<top <<" bin: " <<x[0] <<std::endl;
	std::cout<<"W = " <<W_ <<" bin: " <<x[1] <<std::endl;
	std::cout<<"Q2 = " <<Q2_ <<" bin: " <<x[2] <<std::endl;
	std::cout<<"MM = " <<MM_ <<" bin: " <<x[3] <<std::endl;
	std::cout<<"theta = " <<theta_ <<" bin: " <<x[4] <<std::endl;
	std::cout<<"alpha = " <<alpha_ <<" bin: " <<x[5] <<std::endl;
	std::cout<<"phi = " <<phi_ <<" bin: " <<x[6] <<std::endl;
	*/
	return x; 
}

bool Histogram::Good_Friend_Idx(float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_){
	bool pass = true;
	pass &= (W_ >= _W_min_);
	pass &= (W_ < _W_max_);
	pass &= (Q2_ >= _Q2_min_);
	pass &= (Q2_ < _Q2_max_);
	pass &= (MM_ >= _MM_min_[var_]-_MM_wider_*((_MM_max_[var_][Histogram::W_bin(W_)]-_MM_min_[var_])/_MM_bins_));
	pass &= (MM_ < _MM_max_[var_][Histogram::W_bin(W_)]+_MM_wider_*((_MM_max_[var_][Histogram::W_bin(W_)]-_MM_min_[var_])/_MM_bins_));
	pass &= (MM2_ >= _MM2_min_[var_]-_MM_wider_*((_MM2_max_[var_][Histogram::W_bin(W_)]-_MM2_min_[var_])/_MM_bins_));
	pass &= (MM2_ < _MM2_max_[var_][Histogram::W_bin(W_)]+_MM_wider_*((_MM2_max_[var_][Histogram::W_bin(W_)]-_MM2_min_[var_])/_MM_bins_));
	pass &= (theta_ >= _theta_min_);
	pass &= (theta_ < _theta_max_);
	pass &= (alpha_ >= _alpha_min_);
	pass &= (alpha_ < _alpha_max_);
	pass &= (phi_ >= _phi_min_);
	pass &= (phi_ < _phi_max_);
	/*if(!pass){
		std::cout<<"----fail----\n";
		std::cout<<"W:" <<W_ <<" Q2:" <<Q2_ <<" MM1:" <<MM_ <<" MM2:" <<MM2_ <<" theta:" <<theta_ <<" alpha:" <<alpha_ <<" phi:" <<phi_ <<"\n";
		std::cout<<"W min:" <<(W_ >= _W_min_) <<"\n";
		std::cout<<"W max:" <<(W_ < _W_max_) <<"\n";
		std::cout<<"Q2 min:" <<(Q2_ >= _Q2_min_) <<"\n";
		std::cout<<"Q2 max:" <<(Q2_ < _Q2_max_) <<"\n";
		std::cout<<"MM min:" <<(MM_ >= _MM_min_[var_]-2*((Histogram::MM_max(Histogram::W_bin(W_),var_)-_MM_min_[var_])/_MM_bins_)) <<"\n";
		std::cout<<"MM max:" <<(MM_ < _MM_max_[var_]+2*((Histogram::MM_max(Histogram::W_bin(W_),var_)-_MM_min_[var_])/_MM_bins_)) <<"\n";
		std::cout<<"MM2 min:" <<(MM2_ >= _MM2_min_[var_]-2*((Histogram::MM2_max(Histogram::W_bin(W_),var_)-_MM2_min_[var_])/_MM_bins_)) <<"\n";
		std::cout<<"MM2 max:" <<(MM2_ < _MM2_max_[var_]+2*((Histogram::MM2_max(Histogram::W_bin(W_),var_)-_MM2_min_[var_])/_MM_bins_)) <<"\n";
		std::cout<<"theta min:" <<(theta_ >= _theta_min_) <<"\n";
		std::cout<<"theta max:" <<(theta_ < _theta_max_) <<"\n";
		std::cout<<"alpha min:" <<(alpha_ >= _alpha_min_) <<"\n";
		std::cout<<"alpha max:" <<(alpha_ < _alpha_max_) <<"\n";
		std::cout<<"phi min:" <<(phi_ >= _phi_min_) <<"\n";
		std::cout<<"phi max:" <<(phi_ < _phi_max_) <<"\n";
		std::cout<<"\tMM1: " <<MM_ <<" bounds:" <<_MM_min_[var_]-2*((Histogram::MM_max(Histogram::W_bin(W_),var_)-_MM_min_[var_])/_MM_bins_) <<"-" <<_MM_max_[var_]+2*((Histogram::MM_max(Histogram::W_bin(W_),var_)-_MM_min_[var_])/_MM_bins_) <<"\n\tMM2: " <<MM2_ <<" bounds:" <<_MM2_min_[var_]-2*((Histogram::MM2_max(Histogram::W_bin(W_),var_)-_MM2_min_[var_])/_MM_bins_) <<"-" <<_MM2_max_[var_]+2*((Histogram::MM2_max(Histogram::W_bin(W_),var_)-_MM2_min_[var_])/_MM_bins_) <<"\n\tpass:" <<pass <<"\n";}
	*/
	//if(pass){
	//	std::cout<<"++++pass++++\n";
	//	std::cout<<"\tMM1: " <<MM_ <<" bounds:" <<_MM_min_[var_]-2*((Histogram::MM_max(Histogram::W_bin(W_),var_)-_MM_min_[var_])/_MM_bins_) <<"-" <<_MM_max_[var_]+2*((Histogram::MM_max(Histogram::W_bin(W_),var_)-_MM_min_[var_])/_MM_bins_) <<"\n\tMM2: " <<MM2_ <<" bounds:" <<_MM2_min_[var_]-2*((Histogram::MM2_max(Histogram::W_bin(W_),var_)-_MM2_min_[var_])/_MM_bins_) <<"-" <<_MM2_max_[var_]+2*((Histogram::MM2_max(Histogram::W_bin(W_),var_)-_MM2_min_[var_])/_MM_bins_) <<"\n\tpass:" <<pass <<"\n";}
	return pass;
}

void Histogram::Print_Friend_Bin(float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_, int var_){
	std::cout<<std::endl <<"--Printing Friend idx--" <<std::endl;
	std::cout<<"W: " <<W_ <<" in bin: " <<Histogram::Friend_W_idx(W_) <<std::endl;
	std::cout<<"Q2: " <<Q2_ <<" in bin: " <<Histogram::Friend_Q2_idx(Q2_) <<std::endl;
	std::cout<<"MM: " <<MM_ <<" in bin: " <<Histogram::Friend_MM_idx(MM_,var_,Histogram::W_bin(W_)) <<std::endl;
	std::cout<<"MM2: " <<MM2_ <<" in bin: " <<Histogram::Friend_MM2_idx(MM2_,var_,Histogram::W_bin(W_)) <<std::endl;
	std::cout<<"Theta: " <<theta_ <<" in bin: " <<Histogram::Friend_theta_idx(theta_) <<std::endl;
	std::cout<<"Alpha: " <<alpha_ <<" in bin: " <<Histogram::Friend_alpha_idx(alpha_) <<std::endl;
	std::cout<<"Phi: " <<phi_ <<" in bin: " <<Histogram::Friend_phi_idx(phi_) <<std::endl;
	//std::cout<<"Total Bin: " <<Histogram::Friend_idx(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_) <<std::endl;
}

void Histogram::Friend_Fill(const char* top_, float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_, bool thrown_, float weight_, int helicity_, float plus_weight_, int top_passed_, std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Make_Friend() && ((fun::top_perform(top_,flags_)) || (top_==_mzero_ && thrown_))){
		_n_attempted_evnts+=1;
		bool check_vals = true;
		bool alpha_loop = false;
		bool phi_loop = false;
		float phi = 0.0;
		float alpha = 0.0;
		check_vals &= (W_ >= _W_min_);
		check_vals &= (W_ < _W_max_);
		//if(!check_vals){
		//	std::cout<<"top:" <<top_ <<" var:" <<var_ <<" W:" <<W_ <<" thrown:" <<thrown_ <<"\n";
		//}
		check_vals = true;
		check_vals &= (Q2_ >= _Q2_min_);
		check_vals &= (Q2_ < _Q2_max_);
		/*if(!check_vals){
			//_bvar[var_]+=1;
			std::cout<<"top:" <<top_ <<" var:" <<var_ <<" Q2:" <<Q2_ <<" thrown:" <<thrown_ <<"\n";
		}*/
		check_vals = true;
		check_vals &= (MM_ >= _MM_min_[var_]-_MM_wider_*((_MM_max_[var_][Histogram::W_bin(W_)]-_MM_min_[var_])/_MM_bins_));
		check_vals &= (MM_ < _MM_max_[var_][Histogram::W_bin(W_)]+_MM_wider_*((_MM_max_[var_][Histogram::W_bin(W_)]-_MM_min_[var_])/_MM_bins_));
		/*if(!check_vals){
			_bvar[var_]+=1;
			std::cout<<"top:" <<top_ <<" var:" <<var_ <<" MM:" <<MM_ <<" MM bounds:" <<_MM_min_[var_]-_MM_wider_*((_MM_max_[var_][Histogram::W_bin(W_)]-_MM_min_[var_])/_MM_bins_) <<"-" <<_MM_max_[var_][Histogram::W_bin(W_)]+_MM_wider_*((_MM_max_[var_][Histogram::W_bin(W_)]-_MM_min_[var_])/_MM_bins_) <<" thrown:" <<thrown_ <<"\n";
		}*/
		check_vals = true;
		check_vals &= (MM2_ >= _MM2_min_[var_]-_MM_wider_*((_MM2_max_[var_][Histogram::W_bin(W_)]-_MM2_min_[var_])/_MM_bins_));
		check_vals &= (MM2_ < _MM2_max_[var_][Histogram::W_bin(W_)]+_MM_wider_*((_MM2_max_[var_][Histogram::W_bin(W_)]-_MM2_min_[var_])/_MM_bins_));
		/*if(!check_vals){
			_bvar[var_]+=1;
			std::cout<<"top:" <<top_ <<" var:" <<var_ <<" MM2:" <<MM2_ <<" MM2 bounds:" <<_MM2_min_[var_]-_MM_wider_*((_MM2_max_[var_][Histogram::W_bin(W_)]-_MM2_min_[var_])/_MM_bins_) <<"-" <<_MM2_max_[var_][Histogram::W_bin(W_)]+_MM_wider_*((_MM2_max_[var_][Histogram::W_bin(W_)]-_MM2_min_[var_])/_MM_bins_) <<" thrown:" <<thrown_ <<"\n";
		}*/
		check_vals = true;
		check_vals &= (theta_ >= _theta_min_);
		check_vals &= (theta_ < _theta_max_);
		//if(!check_vals){
		//	std::cout<<"top:" <<top_ <<" var:" <<var_ <<" theta:" <<theta_ <<" theta bounds:" <<_theta_min_ <<"-" <<_theta_max_ <<" thrown:" <<thrown_ <<"\n";
		//}
		check_vals = true;
		check_vals &= (alpha_ >= _alpha_min_);
		check_vals &= (alpha_ < _alpha_max_);
		if(alpha_ == _alpha_max_){
			alpha_loop = true;
			check_vals = true;
		}
		/*if(!check_vals){
			std::cout<<"top:" <<top_ <<" var:" <<var_ <<" alpha:" <<alpha_ <<" alpha bounds:" <<_alpha_min_ <<"-" <<_alpha_max_ <<" thrown:" <<thrown_ <<"\n";
		}*/
		check_vals = true;
		check_vals &= (phi_ >= _phi_min_);
		check_vals &= (phi_ < _phi_max_);
		if(phi_ == _phi_max_){
			check_vals = true;
			phi_loop = true;
		}
		if(!check_vals){
			_bvar[var_]+=1;
			//std::cout<<"top:" <<top_ <<" var:" <<var_  <<" phi:" <<phi_ <<" phi bounds:" <<_phi_min_ <<"-" <<_phi_max_ <<" weight:" <<weight_ <<" helicity:" <<helicity_ <<" thrown:" <<thrown_ <<"\n";
		}
		if(!std::isnan(W_) && !std::isnan(Q2_) && !std::isnan(MM_) && !std::isnan(MM2_) && !std::isnan(theta_) && !std::isnan(alpha_) && !std::isnan(phi_) && !std::isnan(weight_)){
			if(phi_loop){
				phi = 0.0;
			}else{
				phi = phi_;
			}
			if(alpha_loop){
				alpha = 0.0;
			}else{
				alpha = alpha_;
			}
			if(Histogram::Good_Friend_Idx(W_,Q2_,MM_,MM2_,theta_,alpha,phi,var_)){
			//if(Histogram::OK_Idx(Histogram::Friend_idx(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_))){
			//int *y = Friend_binning(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_);
			//std::cout<<std::endl <<"Y binning: " <<y;
				//Histogram::Print_Friend_Bin(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_);
				Double_t x[5] ;//= { (double)W_, (double)Q2_, (double)MM_, (double)MM2_, (double)theta_, (double)alpha_, (double)phi_};
				//x[0] = (double)W_;
				//x[1] = (double)Q2_;
				x[0] = (double)MM_;
				x[1] = (double)MM2_;
				x[2] = (double)theta_;
				x[3] = (double)alpha;
				x[4] = (double)phi;
				if(thrown_){
					//std::lock_guard<std::mutex> lk(std::mutex);//Muting the multithreading for THnSparse filling
					TThread::Lock();
					_Thrown[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(x,weight_/plus_weight_);
					_MM1_Dist_thr[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(MM_,weight_/plus_weight_);
					_MM2_Dist_thr[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(MM2_,weight_/plus_weight_);
					_Theta_Dist_thr[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(theta_,weight_/plus_weight_);
					_Alpha_Dist_thr[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(alpha,weight_/plus_weight_);
					_Phi_Dist_thr[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(phi,weight_/plus_weight_);
					TThread::UnLock();
					//std::cout<<"\tThrown " <<var_  <<" " <<Histogram::W_bin(W_) <<" " <<Histogram::Q2_bin(Q2_) <<"\n"; 
				}else{
					_nevnts+=1;
					_nvar[var_]+=1;
					//std::cout<<"Filling Friend!\n";
					//std::lock_guard<std::mutex> lk(std::mutex);//Muting the multithreading for THnSparse filling
					if(flags_->Flags::Helicity()){//If taking Helicity into account then we'll have two output THnSparse
						if(helicity_ == 1){
							//std::cout<<"\tPos " <<var_ <<" " <<fun::top_idx(top_) <<" " <<Histogram::W_bin(W_) <<" " <<Histogram::Q2_bin(Q2_) <<"\n"; 
							TThread::Lock();
							//_Friend1[var_][fun::top_idx(top_)]->Fill(x,weight_);
							_Friend1[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(x,weight_/plus_weight_);
							TThread::UnLock();
						}else if(helicity_ == -1){
							//std::cout<<"\tNeg " <<var_ <<" " <<fun::top_idx(top_) <<" " <<Histogram::W_bin(W_) <<" " <<Histogram::Q2_bin(Q2_) <<"\n"; 
							TThread::Lock();
							//_Friend2[var_][fun::top_idx(top_)]->Fill(x,weight_);
							_Friend2[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(x,weight_/plus_weight_);
							TThread::UnLock();
							
						}
					}
					TThread::Lock();
					//used to have [var_][fun::top_idx(top_)]
					_Friend[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(x,weight_/plus_weight_);
					_n_wq2[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)][top_passed_]+=1;
					//_W_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(W_,weight_);
					//_Q2_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(Q2_,weight_);
					_MM1_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(MM_,weight_/plus_weight_);
					_MM2_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(MM2_,weight_/plus_weight_);
					_Theta_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(theta_,weight_/plus_weight_);
					_Alpha_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(alpha,weight_/plus_weight_);
					_Phi_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(phi,weight_/plus_weight_);
					TThread::UnLock();
					//std::cout<<"\tNormal " <<var_ <<" " <<fun::top_idx(top_) <<" " <<Histogram::W_bin(W_) <<" " <<Histogram::Q2_bin(Q2_) <<"\n"; 
				}
				//std::cout<<std::endl <<"Filling Friend with " <<x <<" with weight " <<weight_;
			}else{
				_event_npass[var_][top_passed_]+=1;
				if(theta_ >= _theta_max_){
					_nbad_angles[var_][top_passed_][0]+=1;
				}else if(isnan(theta_)){
					_nbad_angles[var_][top_passed_][0]+=1;
				}
				if(alpha_ >= _alpha_max_){
					_nbad_angles[var_][top_passed_][1]+=1;
				}else if(isnan(alpha_)){
					_nbad_angles[var_][top_passed_][1]+=1;
				}
				if(phi_ >= _phi_max_){
					_nbad_angles[var_][top_passed_][2]+=1;
				}else if(isnan(phi_)){
					_nbad_angles[var_][top_passed_][2]+=1;
				}
			}
		}else{
			_event_npass[var_][top_passed_]+=1;
		}
	}
}


void Histogram::Duo_Fill(const char* top_, float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_, bool thrown_, float weight_, int helicity_, float plus_weight_, int top_passed_, std::shared_ptr<Flags> flags_){
	return ;
	if(flags_->Flags::Make_Friend() && ((fun::top_perform(top_,flags_)) || (top_==_mzero_ && thrown_))){
		bool check_vals = true;
		bool alpha_loop = false;
		bool phi_loop = false;
		float phi = 0.0;
		float alpha = 0.0;
		double tmp_duo = 0.0;
		double tmp_duo1 = 0.0;
		double tmp_duo2 = 0.0;
		long tmp_idx = -1;
		long tmp_idx1 = -1;
		long tmp_idx2 = -1;
		check_vals &= (W_ >= _W_min_);
		check_vals &= (W_ < _W_max_);
		if(!check_vals){
			std::cout<<"top:" <<top_ <<" var:" <<var_ <<" W:" <<W_ <<" thrown:" <<thrown_ <<"\n";
		}
		check_vals = true;
		check_vals &= (Q2_ >= _Q2_min_);
		check_vals &= (Q2_ < _Q2_max_);
		if(!check_vals){
			std::cout<<"top:" <<top_ <<" var:" <<var_ <<" Q2:" <<Q2_ <<" thrown:" <<thrown_ <<"\n";
		}
		check_vals = true;
		check_vals &= (MM_ >= _MM_min_[var_]-_MM_wider_*((_MM_max_[var_][Histogram::W_bin(W_)]-_MM_min_[var_])/_MM_bins_));
		check_vals &= (MM_ < _MM_max_[var_][Histogram::W_bin(W_)]+_MM_wider_*((_MM_max_[var_][Histogram::W_bin(W_)]-_MM_min_[var_])/_MM_bins_));
		if(!check_vals){
			std::cout<<"top:" <<top_ <<" var:" <<var_ <<" MM:" <<MM_ <<" MM bounds:" <<_MM_min_[var_]-_MM_wider_*((_MM_max_[var_][Histogram::W_bin(W_)]-_MM_min_[var_])/_MM_bins_) <<"-" <<_MM_max_[var_][Histogram::W_bin(W_)]+_MM_wider_*((_MM_max_[var_][Histogram::W_bin(W_)]-_MM_min_[var_])/_MM_bins_) <<" thrown:" <<thrown_ <<"\n";
		}
		check_vals = true;
		check_vals &= (MM2_ >= _MM2_min_[var_]-_MM_wider_*((_MM2_max_[var_][Histogram::W_bin(W_)]-_MM2_min_[var_])/_MM_bins_));
		check_vals &= (MM2_ < _MM2_max_[var_][Histogram::W_bin(W_)]+_MM_wider_*((_MM2_max_[var_][Histogram::W_bin(W_)]-_MM2_min_[var_])/_MM_bins_));
		if(!check_vals){
			std::cout<<"top:" <<top_ <<" var:" <<var_ <<" MM2:" <<MM2_ <<" MM2 bounds:" <<_MM2_min_[var_]-_MM_wider_*((_MM2_max_[var_][Histogram::W_bin(W_)]-_MM2_min_[var_])/_MM_bins_) <<"-" <<_MM2_max_[var_][Histogram::W_bin(W_)]+_MM_wider_*((_MM2_max_[var_][Histogram::W_bin(W_)]-_MM2_min_[var_])/_MM_bins_) <<" thrown:" <<thrown_ <<"\n";
		}
		check_vals = true;
		check_vals &= (theta_ >= _theta_min_);
		check_vals &= (theta_ < _theta_max_);
		if(!check_vals){
			std::cout<<"top:" <<top_ <<" var:" <<var_ <<" theta:" <<theta_ <<" theta bounds:" <<_theta_min_ <<"-" <<_theta_max_ <<" thrown:" <<thrown_ <<"\n";
		}
		check_vals = true;
		check_vals &= (alpha_ >= _alpha_min_);
		check_vals &= (alpha_ < _alpha_max_);
		if(alpha_ == _alpha_max_){
			alpha_loop = true;
			check_vals = true;
		}
		if(!check_vals){
			std::cout<<"top:" <<top_ <<" var:" <<var_ <<" alpha:" <<alpha_ <<" alpha bounds:" <<_alpha_min_ <<"-" <<_alpha_max_ <<" thrown:" <<thrown_ <<"\n";
		}
		check_vals = true;
		check_vals &= (phi_ >= _phi_min_);
		check_vals &= (phi_ < _phi_max_);
		if(phi_ == _phi_max_){
			check_vals = true;
			phi_loop = true;
		}
		if(!check_vals){
			std::cout<<"top:" <<top_ <<" var:" <<var_  <<" phi:" <<phi_ <<" phi bounds:" <<_phi_min_ <<"-" <<_phi_max_ <<" weight:" <<weight_ <<" helicity:" <<helicity_ <<" thrown:" <<thrown_ <<"\n";
		}
		if(!std::isnan(W_) && !std::isnan(Q2_) && !std::isnan(MM_) && !std::isnan(MM2_) && !std::isnan(theta_) && !std::isnan(alpha_) && !std::isnan(phi_) && !std::isnan(weight_)){
			if(phi_loop){
				phi = 0.0;
			}else{
				phi = phi_;
			}
			if(alpha_loop){
				alpha = 0.0;
			}else{
				alpha = alpha_;
			}
			if(Histogram::Good_Friend_Idx(W_,Q2_,MM_,MM2_,theta_,alpha,phi,var_)){
			//if(Histogram::OK_Idx(Histogram::Friend_idx(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_))){
			//int *y = Friend_binning(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_);
			//std::cout<<std::endl <<"Y binning: " <<y;
				//Histogram::Print_Friend_Bin(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_);
				Double_t x[5] ;//= { (double)W_, (double)Q2_, (double)MM_, (double)MM2_, (double)theta_, (double)alpha_, (double)phi_};
				//x[0] = (double)W_;
				//x[1] = (double)Q2_;
				x[0] = (double)MM_;
				x[1] = (double)MM2_;
				x[2] = (double)theta_;
				x[3] = (double)alpha;
				x[4] = (double)phi;
				if(thrown_){
					//std::lock_guard<std::mutex> lk(std::mutex);//Muting the multithreading for THnSparse filling
					//TThread::Lock();
					//_Thrown[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(x,weight_/plus_weight_);
					//_MM1_Dist_thr[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(MM_,weight_/plus_weight_);
					//_MM2_Dist_thr[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(MM2_,weight_/plus_weight_);
					//_Theta_Dist_thr[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(theta_,weight_/plus_weight_);
					//_Alpha_Dist_thr[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(alpha,weight_/plus_weight_);
					//_Phi_Dist_thr[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(phi,weight_/plus_weight_);
					//TThread::UnLock();
					//std::cout<<"\tThrown " <<var_  <<" " <<Histogram::W_bin(W_) <<" " <<Histogram::Q2_bin(Q2_) <<"\n"; 
				}else{
					//std::cout<<"Filling Friend!\n";
					//std::lock_guard<std::mutex> lk(std::mutex);//Muting the multithreading for THnSparse filling
					if(flags_->Flags::Helicity()){//If taking Helicity into account then we'll have two output THnSparse
						if(helicity_ == 1){
							//std::cout<<"\tPos " <<var_ <<" " <<fun::top_idx(top_) <<" " <<Histogram::W_bin(W_) <<" " <<Histogram::Q2_bin(Q2_) <<"\n"; 
							TThread::Lock();
							//_Friend1[var_][fun::top_idx(top_)]->Fill(x,weight_);
							tmp_idx1 = _Duo1[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->GetBin(x);
							tmp_duo1 = _Duo1[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->GetBinContent(tmp_idx1);
							_Duo1[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->SetBinContent(tmp_idx1,TMath::Sqrt(tmp_duo1*tmp_duo1 + (weight_/(2*plus_weight_))*(weight_/(2*plus_weight_))));
							//_Duo1[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(x,weight_/plus_weight_);
							TThread::UnLock();
						}else if(helicity_ == -1){
							//std::cout<<"\tNeg " <<var_ <<" " <<fun::top_idx(top_) <<" " <<Histogram::W_bin(W_) <<" " <<Histogram::Q2_bin(Q2_) <<"\n"; 
							TThread::Lock();
							//_Friend2[var_][fun::top_idx(top_)]->Fill(x,weight_);
							tmp_idx2 = _Duo2[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->GetBin(x);
							tmp_duo2 = _Duo2[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->GetBinContent(tmp_idx2);
							_Duo2[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->SetBinContent(tmp_idx2,TMath::Sqrt(tmp_duo2*tmp_duo2 + (weight_/(2*plus_weight_))*(weight_/(2*plus_weight_))));
							//_Duo2[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(x,weight_/plus_weight_);
							TThread::UnLock();
							
						}
					}
					TThread::Lock();
					//used to have [var_][fun::top_idx(top_)]
					//_Friend[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(x,weight_/plus_weight_);
					tmp_idx = _Duo[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->GetBin(x);
					tmp_duo = _Duo[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->GetBinContent(tmp_idx);
					_Duo[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->SetBinContent(tmp_idx,TMath::Sqrt(tmp_duo*tmp_duo + (weight_/(2*plus_weight_))*(weight_/(2*plus_weight_))));
					//_n_wq2[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)][top_passed_]+=1;
					//_W_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(W_,weight_);
					//_Q2_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(Q2_,weight_);
					//_MM1_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(MM_,weight_/plus_weight_);
					//_MM2_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(MM2_,weight_/plus_weight_);
					//_Theta_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(theta_,weight_/plus_weight_);
					//_Alpha_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(alpha,weight_/plus_weight_);
					//_Phi_Dist[var_][Histogram::W_bin(W_)][Histogram::Q2_bin(Q2_)]->Fill(phi,weight_/plus_weight_);
					TThread::UnLock();
					//std::cout<<"\tNormal " <<var_ <<" " <<fun::top_idx(top_) <<" " <<Histogram::W_bin(W_) <<" " <<Histogram::Q2_bin(Q2_) <<"\n"; 
				}
				//std::cout<<std::endl <<"Filling Friend with " <<x <<" with weight " <<weight_;
			}/*else{
				_event_npass[var_][top_passed_]+=1;
				if(theta_ >= _theta_max_){
					_nbad_angles[var_][top_passed_][0]+=1;
				}else if(isnan(theta_)){
					_nbad_angles[var_][top_passed_][0]+=1;
				}
				if(alpha_ >= _alpha_max_){
					_nbad_angles[var_][top_passed_][1]+=1;
				}else if(isnan(alpha_)){
					_nbad_angles[var_][top_passed_][1]+=1;
				}
				if(phi_ >= _phi_max_){
					_nbad_angles[var_][top_passed_][2]+=1;
				}else if(isnan(phi_)){
					_nbad_angles[var_][top_passed_][2]+=1;
				}
			}*/
		}/*else{
			_event_npass[var_][top_passed_]+=1;
		}*/
	}
}



void Histogram::Friend_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Make_Friend()){
		std::cout<<"Writing Friend\n";
		_SparseFile->cd();
		TH1D* check_7d[3][29][5][5];
		TH1D* thr_check_7d[3][29][5][5];
		char hname[100];
		//TH1D* check_5d[5][Histogram::W_bins()][5];
		for(int i = 0; i <3; i++){//Variable Set
			for(int j=0; j<Histogram::W_bins(); j++){//W
				for(int l=0; l<5; l++){//Q2
					if(flags_->Flags::Sim()){
						//std::cout<<"Thrown "  <<i <<" integral: " <<_Thrown[i]->ComputeIntegral();
						
						_Thrown[i][j][l]->Write();
						_MM1_Dist_thr[i][j][l]->Write();
						_MM2_Dist_thr[i][j][l]->Write();
						_Theta_Dist_thr[i][j][l]->Write();
						_Alpha_Dist_thr[i][j][l]->Write();
						_Phi_Dist_thr[i][j][l]->Write();
						//_Thrown[i][k][l]->Write();
						for(int k=0; k<5; k++){
							thr_check_7d[i][j][l][k] = _Thrown[i][j][l]->Projection(k,"E");
							if(k==0){
								sprintf(hname,"Thrown_2#pi_off_proton_%s_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_var_names_[i],_mall_,_MM1_[i],Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(l),Histogram::Q2_top(l));
							}else if(k==1){
								sprintf(hname,"Thrown_2#pi_off_proton_%s_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_var_names_[i],_mall_,_MM2_[i],Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(l),Histogram::Q2_top(l));
							}else{
								sprintf(hname,"Thrown_2#pi_off_proton_%s_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_var_names_[i],_mall_,_friend_pars_[k+2],Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(l),Histogram::Q2_top(l));
							}
							thr_check_7d[i][j][l][k]->SetNameTitle(hname,hname);
							thr_check_7d[i][j][l][k]->Write();
						}
					}
					//for(int j = 0; j<5; j++){//Topology
					if(fun::top_perform(_mall_,flags_)){
							//_Friend[i][j][k][l]->Write();
							//std::cout<<"Friend "  <<i <<" " <<j <<" integral: " <<_Friend[i][j]->ComputeIntegral();
							
						_Friend[i][j][l]->Write();
						//_Duo[i][j][l]->Write();
						//_W_Dist[i][j][l]->Write();
						//_Q2_Dist[i][j][l]->Write();
						_MM1_Dist[i][j][l]->Write();
						_MM2_Dist[i][j][l]->Write();
						_Theta_Dist[i][j][l]->Write();
						_Alpha_Dist[i][j][l]->Write();
						_Phi_Dist[i][j][l]->Write();
						for(int k=0; k<5; k++){
							check_7d[i][j][l][k] = _Friend[i][j][l]->Projection(k,"E");
							if(k==0){
								sprintf(hname,"2#pi_off_proton_%s_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_var_names_[i],_mall_,_MM1_[i],Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(l),Histogram::Q2_top(l));
							}else if(k==1){
								sprintf(hname,"2#pi_off_proton_%s_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_var_names_[i],_mall_,_MM2_[i],Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(l),Histogram::Q2_top(l));
							}else{
								sprintf(hname,"2#pi_off_proton_%s_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_var_names_[i],_mall_,_friend_pars_[k+2],Histogram::W_bot(j),Histogram::W_top(j),Histogram::Q2_bot(l),Histogram::Q2_top(l));
							}
							check_7d[i][j][l][k]->SetNameTitle(hname,hname);
							check_7d[i][j][l][k]->Write();
						}
							
							
						if(flags_->Flags::Helicity()){
							
							_Friend1[i][j][l]->Write();
							
							_Friend2[i][j][l]->Write();
							//_Duo1[i][j][l]->Write();
							//_Duo2[i][j][l]->Write();
							
							//_Friend1[i][j][k][l]->Write();
							//_Friend2[i][j][k][l]->Write();
						}
							//for(int m=0; m<4; m++){
								//for(int n=0; n<_Friend[0][0][0][0]->GetAxis(m)->GetNbins(); n++){//Xij other than Phi
									//_Friend_Phi[i][fun::top_idx(_top_[j])+fun::top_offset(_top_[j],flags_)][k][l][m][n]->Write();
									//if(flags_->Flags::Helicity()){
									//	_Friend1_Phi[i][fun::top_idx(_top_[j])+fun::top_offset(_top_[j],flags_)][k][l][m][n]->Write();
									//	_Friend2_Phi[i][fun::top_idx(_top_[j])+fun::top_offset(_top_[j],flags_)][k][l][m][n]->Write();
									//}
									//if(flags_->Sim() && j==0){
									//	_Thrown_Phi[i][k][l][m][n]->Write();
									//}
								//}
							//}
					}
				}
			}
		}
	}
}

//*------------------------------- Start Check ---------------------------------*
void Histogram::PCorr_Check_Make(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Check()){ return ;}

	TH1F_ptr_1d plot_1d;
	TH1F_ptr_2d plot_2d;


	std::vector<long> space_dims(3);
	space_dims[0] = 6; //Sectors
	space_dims[1] = 4; //Missing Mass
	space_dims[2] = 3; //Corrections {none,e_theta,e_p}
	CartesianGenerator cart(space_dims);
	char hname[100];

	for(int i=0; i<6; i++){
		sprintf(hname,"Elastic Theta Relation Sec:%d",i+1);
		_Elastic_Theta_hist[i] = new TH2F(hname,hname,300,0.0,55.0,300,0.0,55.0);
	}
	sprintf(hname,"Elastic Theta Relation Sec:All");
	_Elastic_Theta_hist[6] = new TH2F(hname,hname,300,0.0,55.0,300,0.0,55.0);

	while(cart.GetNextCombination()){
		//std::cout<<"Trying to make: sector:" <<cart[0]+1 <<" ele_corr:" <<_ele_corr_[cart[2]] <<" pro_thres:" <<_proton_threshold_[cart[1]] <<"\n";
		//std::cout<<"Making Check Histogram: " <<_Check_hist.size() <<" " <<plot_1d.size() <<"\n";
		if((!flags_->Flags::E_Theta_Corr() && _ele_corr_[cart[2]]==_ele_angle_corr_) || (!flags_->Flags::E_PCorr() && _ele_corr_[cart[2]]==_ele_p_corr_)){
			//
		}else{
			sprintf(hname,"Ele_PCorr_Check_%s_%s_Sector:%d",_top_[cart[1]],_ele_corr_[cart[2]],cart[0]+1);
			plot_1d.push_back(new TH1F(hname,hname,_mm_bin_[cart[1]],_mm_min_[cart[1]],_mm_max_[cart[1]]));
			if(cart[0] == space_dims[0]-1){
				if(plot_1d.size()>0){
					plot_2d.push_back(plot_1d);
					plot_1d.clear();
				}
				if(cart[1]== space_dims[1]-1){
					if(plot_2d.size()>0){
						_PCorr_Check_hist.push_back(plot_2d);
						plot_2d.clear();
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::PCorr_Check_idx(int sector_, const char* top_, const char* corr_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	//std::cout<<"Trying to fill in " <<sector_ <<" " <<check_ <<"\n";
	if(!flags_->Plot_Check()){
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		return idx;
	}
	if(flags_->Flags::E_Theta_Corr() && corr_ == _ele_angle_corr_){
		idx.push_back(1);
	}else if(flags_->Flags::E_PCorr() && corr_ == _ele_p_corr_){
		idx.push_back(2);
	}else if(corr_== _no_corr_){
		idx.push_back(0);
	}else{
		idx.push_back(-1);
	}
	idx.push_back(fun::top_idx(top_));
	idx.push_back(sector_-1);
	//std::cout<<"idx: " <<idx[0] <<" " <<idx[1] <<"\n";
	return idx;
}

void Histogram::PCorr_Check_Fill(float MM2_, int sector_, const char* top_, const char* corr_, std::shared_ptr<Flags> flags_){
	//std::cout<<"\tFilling Check Hist for " <<xval_ <<" " <<yval_ <<" " <<sector_ <<" " <<check_ <<"\n";
	std::vector<int> idx = PCorr_Check_idx(sector_,top_, corr_,flags_);
	if(OK_Idx(idx) && flags_->Flags::Plot_Check()){
		//std::cout<<"\t\tGood Index! Let's fill: " <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<"\n";
		_PCorr_Check_hist[idx[0]][idx[1]][idx[2]]->Fill(MM2_);
	}
}

void Histogram::PCorr_Elast_Theta_Fill(float theta_e, float theta_p, int sector, std::shared_ptr<Flags> flags_){
	//std::cout<<"trying to Pcorr_Elast_Theta_Fill with " <<theta_e <<" " <<theta_p <<" and " <<sector <<"\n";
	if(!flags_->Flags::Plot_Check()){ return;}
	if(isnan(theta_e)){return ;}
	if(isnan(theta_p)){return ;}
	_Elastic_Theta_hist[sector-1]->Fill(theta_p, theta_e);
	_Elastic_Theta_hist[6]->Fill(theta_p, theta_e);
	//std::cout<<"Success \n";
}

void Histogram::PCorr_Check_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Check()){ return ;}
	std::cout<<"Writing PCorr Check Plots\n";
	char dir_name[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_check = _RootOutputFile->mkdir("PCorr Check Plots");
	//std::cout<<" Done\n";
	dir_check->cd();
	//std::cout<<"Making Sub Directories: ";
	TDirectory* dir_sub[6][5];//{sector,proton_thesh}
	for(int i=0; i<6; i++){
		sprintf(dir_name,"PCorr Check Plots Sector:%d",i+1);
		dir_sub[i][0] = dir_check->mkdir(dir_name);
		for(int j=0; j<3; j++){
			sprintf(dir_name,"PCorr Check Plots Sector:%d %s",i+1,_ele_corr_[j]);
			dir_sub[i][j+1] = dir_sub[i][0]->mkdir(dir_name);
		}
	}
	std::vector<long> space_dims(3);
	space_dims[0] = 6; //Sectors
	space_dims[1] = 4; //Topologies
	space_dims[2] = 3; //corrections 
	CartesianGenerator cart(space_dims);
	char hname[100];
	std::vector<int> idx;

	for(int i=0; i<6; i++){
		dir_sub[i][0]->cd();
		_Elastic_Theta_hist[i]->SetXTitle("Theta p (deg)");
		_Elastic_Theta_hist[i]->SetYTitle("Theta e (deg)");
		_Elastic_Theta_hist[i]->Write();
	}
	dir_check->cd();
	_Elastic_Theta_hist[6]->SetXTitle("Theta p (deg)");
	_Elastic_Theta_hist[6]->SetYTitle("Theta e (deg)");
	_Elastic_Theta_hist[6]->Write();

	while(cart.GetNextCombination()){
		idx = PCorr_Check_idx(cart[0]+1,_top_[cart[1]],_ele_corr_[cart[2]],flags_);
		if(OK_Idx(idx)){
			dir_sub[cart[0]][cart[2]+1]->cd();
			_PCorr_Check_hist[idx[0]][idx[1]][idx[2]]->SetXTitle("MM^2 (GeV)");
			_PCorr_Check_hist[idx[0]][idx[1]][idx[2]]->SetYTitle("Yield");
			_PCorr_Check_hist[idx[0]][idx[1]][idx[2]]->Write();
		}
	}
}

//*------------------------------- End Check2 ---------------------------------*
//*------------------------------Start Bin Centering Corrections------------------*
double Histogram::Xij_Bin_Min(int bin_, int Xij_, int W_bin_, int var_set_){
	if(Xij_>1){
		if(Xij_ == 2){
			return bin_*18.0;
		}else{
			return bin_*36.0;
		}
	}else{
		if(Xij_==0){
			return _MM_min_[var_set_]-(_MM_wider_-bin_)*((_MM_max_[var_set_][W_bin_]-_MM_min_[var_set_])/_MM_bins_);
		}else if(Xij_==1){
			return _MM2_min_[var_set_]-(_MM_wider_-bin_)*((_MM2_max_[var_set_][W_bin_]-_MM2_min_[var_set_])/_MM_bins_);
		}
	}
}

double Histogram::Xij_Bin_Max(int bin_, int Xij_, int W_bin_, int var_set_){
	if(Xij_>1){
		if(Xij_ == 2){
			return (bin_+1)*18.0;
		}else{
			return (bin_+1)*36.0;
		}
	}else{
		if(Xij_==0){
			return _MM_min_[var_set_]-(_MM_wider_-bin_-1)*((_MM_max_[var_set_][W_bin_]-_MM_min_[var_set_])/_MM_bins_);
		}else if(Xij_==1){
			return _MM2_min_[var_set_]-(_MM_wider_-bin_-1)*((_MM2_max_[var_set_][W_bin_]-_MM2_min_[var_set_])/_MM_bins_);
		}
	}
}

void Histogram::Bin_Centering_Make(std::shared_ptr<Flags> flags_){
	//[W][Q2][Var Set][Xij][Bin Xij][W,Q2,Xij projection]
	if(!flags_->Plot_Bin_Centering()){
		return;
	}
	std::cout<<"Making Bin Centering Correction Histograms\n";
	std::vector<long> space_dims(4);
	space_dims[3] = 29; //W Bins
	space_dims[2] = 5; //Q2 Bins
	space_dims[1] = 3; //Var Set
	space_dims[0] = 5; //Xij
	CartesianGenerator cart(space_dims);
	char hname[100];
	TH1D_ptr_1d hist_1d;
	TH1D_ptr_2d hist_2d;
	TH1D_ptr_3d hist_3d;
	TH1D_ptr_4d hist_4d;
	TH1D_ptr_5d hist_5d;
	double bot;
	double top;
	int num_bins_xij[5] = {_MM_bins_+2*_MM_wider_,_MM_bins_+2*_MM_wider_,_theta_bins_,_alpha_bins_,_phi_bins_};
	int num_bins_fin = 0;
	int Wbin = 0;
	int Q2bin = 0;
	int Xij=0;
	int var=0;
	char proj_name[100];
	char* xij_names[] = {"MM1","MM2","theta","alpha","phi"};
	while(cart.GetNextCombination()){
		Wbin = cart[3];
		Q2bin=cart[2];
		Xij=cart[0];
		var=cart[1];
		for(int j=0; j<num_bins_xij[Xij]; j++){
			for(int i=0; i<3; i++){
				switch(i){
					case 0: 
						num_bins_fin = 29;//W
						sprintf(proj_name,"W");
						bot = Histogram::W_bot(Wbin);
						top = Histogram::W_top(Wbin);
					break;
					case 1: 
						num_bins_fin = 5;//Q2
						sprintf(proj_name,"Q2");
						bot = Histogram::Q2_bot(Q2bin);
						top = Histogram::Q2_top(Q2bin);
					break;
					case 2: 
						num_bins_fin = num_bins_xij[Xij];//Xij
						sprintf(proj_name,"%s_%s",_var_names_[var],xij_names[Xij]);
						bot = Histogram::Xij_Bin_Min(j,Xij,Wbin,var);
						top = Histogram::Xij_Bin_Max(j,Xij,Wbin,var);
					break;
					default:
						std::cout<<"Improper variable choice for bin centering\n";
					break;
				}
				sprintf(hname,"%s_Bin_Center_W:%.3f-%.3f_Q2:%.2f-%.2f_%s_%s:%.4f-%.4f",proj_name,Histogram::W_bot(Wbin),Histogram::W_top(Wbin),Histogram::Q2_bot(Q2bin),Histogram::Q2_top(Q2bin),_var_names_[var],xij_names[Xij],Histogram::Xij_Bin_Min(j,Xij,Wbin,var),Histogram::Xij_Bin_Max(j,Xij,Wbin,var));
				hist_1d.push_back(new TH1D(hname,hname,11,bot,top));
			}
			hist_2d.push_back(hist_1d);
			hist_1d.clear();
		}
		if(hist_2d.size()>0){
			hist_3d.push_back(hist_2d);
			hist_2d.clear();
		}
		if(cart[0] == space_dims[0]-1){
			if(hist_3d.size()>0){
				hist_4d.push_back(hist_3d);
				hist_3d.clear();
			}
			if(cart[1] == space_dims[1]-1){
				if(hist_4d.size()>0){
					hist_5d.push_back(hist_4d);
					hist_4d.clear();
				}
				if(cart[2] == space_dims[2]-1){
					if(hist_5d.size()>0){
						_Bin_Center_hist.push_back(hist_5d);
						hist_5d.clear();
					}
				}
			}
		}
	}
	std::cout<<"\tFinished Making Bin Centering Histograms\n";
}

std::vector<int> Histogram::Bin_Centering_idx(float W_, float Q2_, int var_set_, int Xij_, double Xij_val_, int variable_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	idx.push_back(Histogram::W_bin(W_));
	idx.push_back(Histogram::Q2_bin(Q2_));
	idx.push_back(var_set_);
	idx.push_back(Xij_);
	switch(Xij_){
		case 0:
			idx.push_back(Histogram::Friend_MM_idx(Xij_val_,var_set_,Histogram::W_bin(W_)));
		break;
		case 1:
			idx.push_back(Histogram::Friend_MM2_idx(Xij_val_,var_set_,Histogram::W_bin(W_)));
		break;
		case 2:
			idx.push_back(Histogram::Friend_theta_idx(Xij_val_));
		break;
		case 3:
			idx.push_back(Histogram::Friend_alpha_idx(Xij_val_));
		break;
		case 4:
			idx.push_back(Histogram::Friend_phi_idx(Xij_val_));
		break;
		default:
			idx.push_back(-1);
		break;
	}
	idx.push_back(variable_);
	return idx;
}

void Histogram::Bin_Centering_Fill(float W_, float Q2_, int var_set_, int Xij_, double Xij_val_, int variable_, double weight_, std::shared_ptr<Flags> flags_){
	if(!flags_->Plot_Bin_Centering()){
		return;
	}
	std::vector<int> idx = Histogram::Bin_Centering_idx(W_,Q2_,var_set_,Xij_,Xij_val_,variable_,flags_);
	//std::cout<<"Filling Bin Centering Corrections with W:" <<W_ <<" Q2:" <<Q2_ <<" var_set:" <<var_set_ <<" Xij:" <<Xij_ <<" Xij_val:" <<Xij_val_ <<" variable:" <<variable_ <<" weight:" <<weight_ <<"\n";
	if(Histogram::OK_Idx(idx)){
		switch(variable_){
			case 0:
				_Bin_Center_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(W_,weight_);
			break;
			case 1:
				_Bin_Center_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(Q2_,weight_);
			break;
			case 2:
				_Bin_Center_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(Xij_val_,weight_);
			break;
			default:
				std::cout<<"Improper Variable designation for Bin Centering Corrections\n";
			break;
		}
	}
}

void Histogram::Bin_Centering_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Plot_Bin_Centering()){ return; }
	std::cout<<"Writing Bin Centering Plots\n";
	char dir_name[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_bc = _RootOutputFile->mkdir("Bin Centering Plots");
	//std::cout<<" Done\n";
	dir_bc->cd();
	//std::cout<<"Making Sub Directories: ";
	TDirectory* dir_bc_sub[29][5+1][3+1][5+1][3+1];//{sector,proton_thesh}
	//Histogram::W_bot(Wbin),Histogram::W_top(Wbin),Histogram::Q2_bot(Q2bin),Histogram::Q2_top(Q2bin),Histogram::Xij_Bin_Min(j,Xij,Wbin,var),Histogram::Xij_Bin_Max(j,Xij,Wbin,var)
	char* xij_names[] = {"MM1","MM2","theta","alpha","phi"};
	char proj_name[100];
	for(int i=0; i<29; i++){//W
		sprintf(dir_name,"Bin Centering W:%.3f-%.3f",Histogram::W_bot(i),Histogram::W_top(i));
		dir_bc_sub[i][0][0][0][0] = dir_bc->mkdir(dir_name);
		dir_bc_sub[i][0][0][0][0]->cd();
		for(int j=0; j<5; j++){//Q2
			sprintf(dir_name,"Bin Centering W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_bot(i),Histogram::W_top(i),Histogram::Q2_bot(j),Histogram::Q2_top(j));
			dir_bc_sub[i][j+1][0][0][0] = dir_bc_sub[i][0][0][0][0]->mkdir(dir_name);
			dir_bc_sub[i][j+1][0][0][0]->cd();
			for(int k=0; k<3; k++){//var set
				sprintf(dir_name,"Bin Centering W:%.3f-%.3f Q2:%.2f-%.2f var:%s",Histogram::W_bot(i),Histogram::W_top(i),Histogram::Q2_bot(j),Histogram::Q2_top(j),_var_names_[k]);
				dir_bc_sub[i][j+1][k+1][0][0] = dir_bc_sub[i][j+1][0][0][0]->mkdir(dir_name);
				dir_bc_sub[i][j+1][k+1][0][0]->cd();
				for(int l=0; l<5; l++){//Xij
					sprintf(dir_name,"Bin Centering W:%.3f-%.3f Q2:%.2f-%.2f var:%s Xij:%s",Histogram::W_bot(i),Histogram::W_top(i),Histogram::Q2_bot(j),Histogram::Q2_top(j),_var_names_[k],xij_names[l]);
					dir_bc_sub[i][j+1][k+1][l+1][0] = dir_bc_sub[i][j+1][k+1][0][0]->mkdir(dir_name);
					dir_bc_sub[i][j+1][k+1][l+1][0]->cd();
					for(int m=0; m<3; m++){//{W,Q2,Xij}
						switch(m){
							case 0:
								sprintf(proj_name,"W");
							break;
							case 1:
								sprintf(proj_name,"Q2");
							break;
							case 2:
								sprintf(proj_name,"%s",xij_names[l]);
							break;
							default:
								sprintf(proj_name,"wrong");
							break;
						}
						sprintf(dir_name,"Bin Centering W:%.3f-%.3f Q2:%.2f-%.2f var:%s Xij:%s proj:%s",Histogram::W_bot(i),Histogram::W_top(i),Histogram::Q2_bot(j),Histogram::Q2_top(j),_var_names_[k],xij_names[l],proj_name);
						dir_bc_sub[i][j+1][k+1][l+1][m+1] = dir_bc_sub[i][j+1][k+1][l+1][0]->mkdir(dir_name);
					}
				}
			}
		}
	}
	std::vector<long> space_dims(4);
	space_dims[3] = 29; //W Bins
	space_dims[2] = 5; //Q2 Bins
	space_dims[1] = 3; //Var Set
	space_dims[0] = 5; //Xij
	CartesianGenerator cart(space_dims);
	char hname[100];
	std::vector<int> idx;
	int Wbin, Q2bin, Xij, var;
	double val, Xij_val;
	int num_bins_xij[5] = {_MM_bins_+2*_MM_wider_,_MM_bins_+2*_MM_wider_,_theta_bins_,_alpha_bins_,_phi_bins_};
	char * xij_units[5] = {"GeV","GeV","Degrees","Degrees","Degrees"};
	double top,bot;
	char xlabel[100];
	while(cart.GetNextCombination()){
		Wbin = cart[3];
		Q2bin=cart[2];
		Xij=cart[0];
		var=cart[1];
		//std::cout<<"W bin:" <<Wbin <<" Q2 bin:" <<Q2bin <<" var:" <<var <<" Xij:" <<Xij <<" going into loop\n";
		for(int j=0; j<num_bins_xij[Xij]; j++){
			Xij_val = (Xij_Bin_Min(j,Xij,Wbin,var)+Xij_Bin_Max(j,Xij,Wbin,var))/2.0;
			
			for(int i=0; i<3; i++){//Variable projection
				switch(i){
					case 0: 
						bot = Histogram::W_bot(Wbin);
						top = Histogram::W_top(Wbin);
						sprintf(xlabel,"W (GeV)");
					break;
					case 1: 
						bot = Histogram::Q2_bot(Q2bin);
						top = Histogram::Q2_top(Q2bin);
						sprintf(xlabel,"Q2 (GeV^2)");
					break;
					case 2: 
						bot = Histogram::Xij_Bin_Min(j,Xij,Wbin,var);
						top = Histogram::Xij_Bin_Max(j,Xij,Wbin,var);
						sprintf(xlabel,"%s %s",xij_names[Xij],xij_units[Xij]);
					break;
				}
				val = (top+bot)/2.0;
				//std::cout<<"W: " <<W_center(Wbin) <<" Q2:" <<(Q2_bot(Q2bin)+Q2_top(Q2bin))/2.0 <<" var:" <<var <<" Xij:" <<Xij <<" Xij_val:" <<Xij_val <<" variable:" <<i <<"\n";
				idx = Histogram::Bin_Centering_idx(W_center(Wbin),(Q2_bot(Q2bin)+Q2_top(Q2bin))/2.0,var,Xij,Xij_val,i,flags_);
				//fun::print_vector_idx(idx);
				if(Histogram::OK_Idx(idx)){
					dir_bc_sub[Wbin][Q2bin+1][var+1][Xij+1][i+1]->cd();
					_Bin_Center_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle(xlabel);
					_Bin_Center_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Yield");
					_Bin_Center_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
				}
			}
		}
	}
	std::cout<<"\tFinished Writing Bin Centering Plots\n";
}

//*------------------------------End Bin Centering Corrections------------------*
void Histogram::Clean_Top_Increment(int i_){
	_ctop[i_]+=1;
}

void Histogram::Isolated_Top_Increment(int i_){
	_itop[i_]+=1;
}

long Histogram::Clean_Top(int i_){
	return _ctop[i_];
}

long Histogram::Isolated_Top(int i_){
	return _itop[i_];
}
void Histogram::Clean_Not_Isolated_Top_Increment(int i_){
	_cnitop[i_]+=1;
}
long Histogram::Clean_Not_Isolated_Top(int i_){
	return _cnitop[i_];
}
void Histogram::Clean_and_Isolated_Top_Increment(int i_){
	_caitop[i_]+=1;
}
long Histogram::Clean_and_Isolated_Top(int i_){
	return _caitop[i_];
}

void Histogram::Top_Increment(int i_){
	_ntop[i_]+=1;
}

void Histogram::Top_Pot_Increment(int i_){
	_nptop[i_]+=1;
}
long Histogram::NTop(int i_){
	return _ntop[i_];
}

long Histogram::NPTop(int i_){
	return _nptop[i_];
}

long Histogram::N_Yield(int var_, int Wbin_, int Q2bin_, int top_){
	return _n_wq2[var_][Wbin_][Q2bin_][top_];
}

long Histogram::Bad_Angles(int var_, int top_, int angle_){
	return _nbad_angles[var_][top_][angle_];
}

long Histogram::Ev_No_Pass(int var_, int top_){
	return _event_npass[var_][top_];
}

//*------------------------------Start Electron Angle Corrections------------------*
void Histogram::Ele_Angle_Corr_Make(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Electron_Angle_Corr()){ return;}
	std::cout<<"Making Electron Angle Correction Histograms\n";
	std::vector<long> space_dims(3);
	space_dims[2] = 6; //Sector
	space_dims[1] = (int)((max_theta_pecorr-min_theta_pecorr)/res_theta_pecorr); //Theta
	space_dims[0] = (int)((max_phi_pecorr-min_phi_pecorr)/res_phi_pecorr); //phi
	char hname[100];
	CartesianGenerator cart(space_dims);
	TH1F_ptr_1d tmp_1d;
	TH1F_ptr_2d tmp_2d;

	TH2F_ptr_1d tmp2_1d;
	for(int sec=0; sec<7; sec++){
		sprintf(hname,"Electron-Proton Elastic Angles %s",_sector_[sec]);
		_Pecorr_Angle_Dist_hist.push_back(new TH2F(hname,hname,_fid_ybin_,_fid_ymin_,_fid_ymax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
		for(int top=0; top<5; top++){
			sprintf(hname,"Electron-Proton Elastic Angles %s top:%s",_sector_[sec],_top_[top]);
			tmp2_1d.push_back(new TH2F(hname,hname,_fid_ybin_,_fid_ymin_,_fid_ymax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
		}
		_Pecorr_Angle_Dist2_hist.push_back(tmp2_1d);
		tmp2_1d.clear();
	}
	while(cart.GetNextCombination()){
		sprintf(hname,"Electron_Angle_Diff_Theta:%.1f-%.1f_Phi:%.1f-%.1f_Sec:%d",min_theta_pecorr+cart[1]*res_theta_pecorr,min_theta_pecorr+(cart[1]+1)*res_theta_pecorr,min_phi_pecorr+cart[0]*res_phi_pecorr,min_phi_pecorr+(cart[0]+1)*res_phi_pecorr,cart[2]+1);
		tmp_1d.push_back(new TH1F(hname,hname,bins_delta_theta,min_delta_theta,max_delta_theta));
		if(cart[0] == space_dims[0]-1){
			if(tmp_1d.size()>0){
				tmp_2d.push_back(tmp_1d);
				tmp_1d.clear();
			}
			if(cart[1] == space_dims[1]-1){
				if(tmp_2d.size()>0){
					_Pecorr_Angle_hist.push_back(tmp_2d);
					tmp_2d.clear();
				}
			}
		}
	}
	std::cout<<"\tFinished Making Electron Angle Correction Histograms\n";
}
int Histogram::Ele_Angle_Corr_Theta_idx(float theta_){
	int output = -1;
	for(int i=0; i<(int)((max_theta_pecorr-min_theta_pecorr)/res_theta_pecorr); i++){
		if(theta_ >= min_theta_pecorr+i*res_theta_pecorr && theta_ < min_theta_pecorr+(i+1)*res_theta_pecorr){
			output = i;
		}
	}
	return output;
}
int Histogram::Ele_Angle_Corr_Phi_idx(float phi_){
	int output = -1;
	for(int i=0; i<(int)((max_phi_pecorr-min_phi_pecorr)/res_phi_pecorr); i++){
		if(phi_ >= min_phi_pecorr+i*res_phi_pecorr && phi_ < min_phi_pecorr+(i+1)*res_phi_pecorr){
			output = i;
		}
	}
	return output;
}
std::vector<int> Histogram::Ele_Angle_Corr_idx(int sector_, float theta_, float phi_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx; 
	if(!flags_->Flags::Plot_Electron_Angle_Corr()){ 
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		return idx;
	}
	idx.push_back(sector_-1);
	idx.push_back(Histogram::Ele_Angle_Corr_Theta_idx(theta_));
	idx.push_back(Histogram::Ele_Angle_Corr_Phi_idx(phi_));
	return idx;
}
void Histogram::Ele_Angle_Corr_Fill(int sector_, float theta_, float phi_, float W_, float thetap_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Electron_Angle_Corr()){ return;}
	float min_proton_theta = 15.0;//5/17/24
	if(W_>= 0.7 && W_ < 1.05){
		_Pecorr_Angle_Dist_hist[sector_-1]->Fill(theta_,thetap_);
		_Pecorr_Angle_Dist_hist[6]->Fill(theta_,thetap_);
		if(thetap_ >= min_proton_theta){
			std::vector<int> idx = Histogram::Ele_Angle_Corr_idx(sector_,theta_,phi_,flags_);
			if(Histogram::OK_Idx(idx)){
				//std::cout<<"filling Ele_Mag_Corr W:" <<W_ <<" theta:" <<theta_ <<" phi:" <<phi_ <<" sector: " <<sector_ <<" to get index:" <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<"\n";
				_Pecorr_Angle_hist[idx[0]][idx[1]][idx[2]]->Fill(physics::delta_theta_e(theta_, thetap_, flags_->Flags::Run()));
			}
		}
	}
}

void Histogram::Ele_Pro_Angle_Dist_Fill(const char* sector_, float theta_, float thetap_, const char* top_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Electron_Angle_Corr()){ return;}
	int sec_idx = fun::sector_idx(sector_);
	_Pecorr_Angle_Dist2_hist[sec_idx][fun::top_idx(top_)]->Fill(theta_,thetap_);
	_Pecorr_Angle_Dist2_hist[6][fun::top_idx(top_)]->Fill(theta_,thetap_);
}

void Histogram::Ele_Angle_Corr_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Electron_Angle_Corr()){ return;}
	std::cout<<"Writing Electron Angle Correction Histograms\n";
	char dirname[100];
	TDirectory* dir_pecorr_angle = _RootOutputFile->mkdir("Electron Angle Corrections");
	//std::cout<<" Done\n";
	dir_pecorr_angle->cd();
	//std::cout<<"Making Sub Directories: ";
	TDirectory* dir_pecorr_angle_sub[6][(int)((max_theta_pecorr-min_theta_pecorr)/res_theta_pecorr)+1];//{sector,proton_thesh}
	//Histogram::W_bot(Wbin),Histogram::W_top(Wbin),Histogram::Q2_bot(Q2bin),Histogram::Q2_top(Q2bin),Histogram::Xij_Bin_Min(j,Xij,Wbin,var),Histogram::Xij_Bin_Max(j,Xij,Wbin,var)
	for(int i=0; i<6; i++){
		sprintf(dirname,"Ele Angle Corr Sector:%d",i+1);
		dir_pecorr_angle_sub[i][0] = dir_pecorr_angle->mkdir(dirname);
		for(int j=0; j<(int)((max_theta_pecorr-min_theta_pecorr)/res_theta_pecorr); j++){
			sprintf(dirname,"Ele Angle Corr Sector:%d Theta:%.1f-%.1f",i+1,min_theta_pecorr+j*res_theta_pecorr,min_theta_pecorr+(j+1)*res_theta_pecorr);
			dir_pecorr_angle_sub[i][j+1] = dir_pecorr_angle_sub[i][0]->mkdir(dirname);
		}
	}


	std::vector<long> space_dims(3);
	space_dims[2] = 6; //Sector
	space_dims[1] = (int)((max_theta_pecorr-min_theta_pecorr)/res_theta_pecorr); //Theta
	space_dims[0] = (int)((max_phi_pecorr-min_phi_pecorr)/res_phi_pecorr); //phi
	char hname[100];
	CartesianGenerator cart(space_dims);
	TH1F_ptr_1d tmp_1d;
	TH1F_ptr_2d tmp_2d;
	dir_pecorr_angle->cd();
	for(int sec=0; sec<7; sec++){
		_Pecorr_Angle_Dist_hist[sec]->SetXTitle("Electron Theta (deg)");
		_Pecorr_Angle_Dist_hist[sec]->SetYTitle("Proton Theta (deg)");
		_Pecorr_Angle_Dist_hist[sec]->Write();
		for(int top = 0; top<5; top++){
			_Pecorr_Angle_Dist2_hist[sec][top]->SetXTitle("Electron Theta (deg)");
			_Pecorr_Angle_Dist2_hist[sec][top]->SetYTitle("Proton Theta (deg)");
			_Pecorr_Angle_Dist2_hist[sec][top]->Write();
		}
	}

	while(cart.GetNextCombination()){
		dir_pecorr_angle_sub[cart[2]][cart[1]+1]->cd();
		_Pecorr_Angle_hist[cart[2]][cart[1]][cart[0]]->SetXTitle("Delta Theta (deg)");
		_Pecorr_Angle_hist[cart[2]][cart[1]][cart[0]]->SetYTitle("Yield");
		_Pecorr_Angle_hist[cart[2]][cart[1]][cart[0]]->Write();
	}
	std::cout<<"Finished Writing Electron Angle Correction Histograms\n";
}
//*------------------------------End Electron Angle Corrections------------------*
//*------------------------------Start Electron Momentum Magnitude Corrections------------------*
void Histogram::Ele_Mag_Corr_Make(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Electron_Mag_Corr()){ return;}
	std::cout<<"Making Electron Mag Correction Histograms\n";
	std::vector<long> space_dims(3);
	space_dims[2] = 6; //Sector
	space_dims[1] = (int)((max_theta_pecorr-min_theta_pecorr)/res_theta_pecorr); //Theta
	space_dims[0] = (int)((max_phi_pecorr-min_phi_pecorr)/res_phi_pecorr); //phi
	char hname[100];
	CartesianGenerator cart(space_dims);
	TH1F_ptr_1d tmp_1d;
	TH1F_ptr_2d tmp_2d;
	while(cart.GetNextCombination()){
		sprintf(hname,"Electron_Mag_Diff_Theta:%.1f-%.1f_Phi:%.1f-%.1f_Sec:%d",min_theta_pecorr+cart[1]*res_theta_pecorr,min_theta_pecorr+(cart[1]+1)*res_theta_pecorr,min_phi_pecorr+cart[0]*res_phi_pecorr,min_phi_pecorr+(cart[0]+1)*res_phi_pecorr,cart[2]+1);
		tmp_1d.push_back(new TH1F(hname,hname,bins_delta_p,min_delta_p,max_delta_p));
		if(cart[0] == space_dims[0]-1){
			if(tmp_1d.size()>0){
				tmp_2d.push_back(tmp_1d);
				tmp_1d.clear();
			}
			if(cart[1] == space_dims[1]-1){
				if(tmp_2d.size()>0){
					_Pecorr_Mag_hist.push_back(tmp_2d);
					tmp_2d.clear();
				}
			}
		}
	}
	std::cout<<"\tFinished Making Electron Mag Correction Histograms\n";
}
int Histogram::Ele_Mag_Corr_Theta_idx(float theta_){
	int output = -1;
	for(int i=0; i<(int)((max_theta_pecorr-min_theta_pecorr)/res_theta_pecorr); i++){
		if(theta_ >= min_theta_pecorr+i*res_theta_pecorr && theta_ < min_theta_pecorr+(i+1)*res_theta_pecorr){
			output = i;
		}
	}
	return output;
}
int Histogram::Ele_Mag_Corr_Phi_idx(float phi_){
	int output = -1;
	for(int i=0; i<(int)((max_phi_pecorr-min_phi_pecorr)/res_phi_pecorr); i++){
		if(phi_ >= min_phi_pecorr+i*res_phi_pecorr && phi_ < min_phi_pecorr+(i+1)*res_phi_pecorr){
			output = i;
		}
	}
	return output;
}
std::vector<int> Histogram::Ele_Mag_Corr_idx(int sector_, float theta_, float phi_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx; 
	if(!flags_->Flags::Plot_Electron_Mag_Corr()){ 
		idx.push_back(-1);
		idx.push_back(-1);
		idx.push_back(-1);
		return idx;
	}
	idx.push_back(sector_-1);
	idx.push_back(Histogram::Ele_Mag_Corr_Theta_idx(theta_));
	idx.push_back(Histogram::Ele_Mag_Corr_Phi_idx(phi_));
	return idx;
}
void Histogram::Ele_Mag_Corr_Fill(float pe_, int sector_, float theta_, float phi_, float W_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Electron_Mag_Corr()){ return;}
	//std::cout<<"\tEntered Mag Corr Filling with pe:" <<pe_ <<" sec:" <<sector_ <<" theta:" <<theta_ <<" W:" <<W_ <<"\n";
	if(W_>= 0.7 && W_ < 1.05){
		std::vector<int> idx = Histogram::Ele_Mag_Corr_idx(sector_,theta_,phi_,flags_);
		if(Histogram::OK_Idx(idx)){
			//std::cout<<"filling Ele_Mag_Corr with delta p: " <<physics::delta_p_e(pe_,theta_, flags_->Flags::Run()) <<"W:" <<W_ <<" theta:" <<theta_ <<" phi:" <<phi_ <<" sector: " <<sector_ <<" to get index:" <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<"\n";
			_Pecorr_Mag_hist[idx[0]][idx[1]][idx[2]]->Fill(physics::delta_p_e(pe_,theta_, flags_->Flags::Run()));
		}
	}
}

void Histogram::Ele_Mag_Corr_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Electron_Mag_Corr()){ return;}
	std::cout<<"Writing Electron Mag Correction Histograms\n";
	char dirname[100];
	TDirectory* dir_pecorr_Mag = _RootOutputFile->mkdir("Electron Mag Corrections");
	//std::cout<<" Done\n";
	dir_pecorr_Mag->cd();
	//std::cout<<"Making Sub Directories: ";
	TDirectory* dir_pecorr_Mag_sub[6][(int)((max_theta_pecorr-min_theta_pecorr)/res_theta_pecorr)+1];//{sector,proton_thesh}
	//Histogram::W_bot(Wbin),Histogram::W_top(Wbin),Histogram::Q2_bot(Q2bin),Histogram::Q2_top(Q2bin),Histogram::Xij_Bin_Min(j,Xij,Wbin,var),Histogram::Xij_Bin_Max(j,Xij,Wbin,var)
	for(int i=0; i<6; i++){
		sprintf(dirname,"Ele Mag Corr Sector:%d",i+1);
		dir_pecorr_Mag_sub[i][0] = dir_pecorr_Mag->mkdir(dirname);
		for(int j=0; j<(int)((max_theta_pecorr-min_theta_pecorr)/res_theta_pecorr); j++){
			sprintf(dirname,"Ele Mag Corr Sector:%d Theta:%.1f-%.1f",i+1,min_theta_pecorr+j*res_theta_pecorr,min_theta_pecorr+(j+1)*res_theta_pecorr);
			dir_pecorr_Mag_sub[i][j+1] = dir_pecorr_Mag_sub[i][0]->mkdir(dirname);
		}
	}


	std::vector<long> space_dims(3);
	space_dims[2] = 6; //Sector
	space_dims[1] = (int)((max_theta_pecorr-min_theta_pecorr)/res_theta_pecorr); //Theta
	space_dims[0] = (int)((max_phi_pecorr-min_phi_pecorr)/res_phi_pecorr); //phi
	char hname[100];
	CartesianGenerator cart(space_dims);
	TH1F_ptr_1d tmp_1d;
	TH1F_ptr_2d tmp_2d;
	while(cart.GetNextCombination()){
		dir_pecorr_Mag_sub[cart[2]][cart[1]+1]->cd();
		_Pecorr_Mag_hist[cart[2]][cart[1]][cart[0]]->SetXTitle("Delta Theta (deg)");
		_Pecorr_Mag_hist[cart[2]][cart[1]][cart[0]]->SetYTitle("Yield");
		_Pecorr_Mag_hist[cart[2]][cart[1]][cart[0]]->Write();
	}
	std::cout<<"Finished Writing Electron Mag Correction Histograms\n";
}
//*------------------------------End Electron Momentum Magnitude Corrections------------------*
//*------------------------------Start Proton Energy Loss Corrections------------------*

//*------------------------------End Proton Energy Loss Corrections------------------*
//*------------------------------Start Elastic Peak------------------*
void Histogram::Elastic_Peak_Make(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Elastic()){ return;}
	std::cout<<"Making Elastic Peak Histograms\n";
	char hname[100];
	char * corr_stuff[3] = {"no_corr","e_theta_corr","e_pcorr"};
	for(int i=0; i<6; i++){
		for(int j=0; j<3; j++){
			sprintf(hname,"Elastic_Peak_Sec:%d_%s",i,corr_stuff[j]);
			_Elastic_Peak_hist[i][j] = new TH1F(hname,hname,300,0.6,1.2);
		}
	}
}

void Histogram::Elastic_Peak_Fill(float W_, int corr_, int sector_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Elastic()){ return;}
	//std::cout<<"Filling Elastic Peak with W" <<W_ <<" corr:" <<corr_ <<" sector:" <<sector_ <<"\n";
	_Elastic_Peak_hist[sector_-1][corr_]->Fill(W_);
}

void Histogram::Elastic_Peak_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Elastic()){ return; }
	//char dirname[100];
	std::cout<<"Writing Elastic Peak Histogram\n";
	TDirectory* dir_elast = _RootOutputFile->mkdir("Elastic Peak");
	dir_elast->cd();
	for(int i=0; i<6; i++){
		for(int j=0; j<3; j++){
			_Elastic_Peak_hist[i][j]->SetXTitle("W (GeV)");
			_Elastic_Peak_hist[i][j]->SetYTitle("Yield");
			_Elastic_Peak_hist[i][j]->Write();
		}
	}
}
//*------------------------------End Elastic Peak------------------*

//*------------------------------Start Detector Geometric Cut Peak------------------*
int Histogram::Geo_Fid_Paddle_Idx(int det_, int seg_idx_){
	int seg_idx = -1;
	if(detectors[det_+1]=="cc"){
		if(seg_idx_ < _geo_fid_cc_segments_*3){
			for(int i=0; i<_geo_fid_cc_sides_; i++){
				if(seg_idx_ >= _geo_fid_cc_segments_*i){
					seg_idx = seg_idx_ - _geo_fid_cc_segments_*i;
				}
			}
		}else{
			seg_idx = _geo_fid_cc_segments_;
		}
	}else if(detectors[det_+1]=="sc"){
		if(seg_idx_ <= _num_sc_paddles_){
			seg_idx = seg_idx_;
		}
	}else if(detectors[det_+1]=="ec"){
		if(seg_idx_ == 0){
			seg_idx = seg_idx_;
		}
	}
	return seg_idx; 
}
int Histogram::Geo_Fid_Side_Idx(int det_, int seg_idx_){
	int side_idx = -1;
	if(detectors[det_+1]=="cc"){
		if(seg_idx_ < _geo_fid_cc_segments_*3){
			for(int i=0; i<_geo_fid_cc_sides_; i++){
				if(seg_idx_ >= _geo_fid_cc_segments_*i){
					side_idx = i;
				}
			}
		}else{
			side_idx = seg_idx_ - _geo_fid_cc_segments_*3;
		}
	}else if(detectors[det_+1]=="sc"){
		if(seg_idx_ <= _num_sc_paddles_){
			side_idx = 0;
		}
	}else if(detectors[det_+1]=="ec"){
		if(seg_idx_ == 0){
			side_idx = 0;
		}
	}
	return side_idx; 
}

int Histogram::Geo_Fid_Seg_Idx(int det_, int pad_idx_, int side_idx_){
	int seg_idx = -1;
	if(detectors[det_+1]==_cc_){
		if(side_idx_ == _geo_fid_cc_sides_){
			if(pad_idx_ == _geo_fid_cc_segments_){
				seg_idx = (_geo_fid_cc_segments_+1)*_geo_fid_cc_sides_;
			}
		}else{
			if(pad_idx_ == _geo_fid_cc_segments_){
				seg_idx = (_geo_fid_cc_segments_)*_geo_fid_cc_sides_+side_idx_;
			}else{
				seg_idx = _geo_fid_cc_segments_*side_idx_+pad_idx_;
			}
		}
	}else if(detectors[det_+1]==_sc_){
		if(side_idx_ == 0){
			if(pad_idx_ >= 0){
				seg_idx = pad_idx_;
			}
		}
	}else if(detectors[det_+1]==_ec_){
		if(side_idx_ == 0){
			if(pad_idx_ >= 0){
				seg_idx = pad_idx_;
			}
		}
	}
	return seg_idx;
}


void Histogram::Geo_Fid_Make(std::shared_ptr<Flags> flags_){
	bool plot_geo = false;
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++){
			if(flags_->Flags::Plot_Fid_Geo(i,j)){
				plot_geo=true;
			}
		}
	}
	bool plot_cc = false;
	bool plot_sc = false;
	bool plot_ec = false;
	for(int i=0; i<4; i++){
		if(flags_->Flags::Plot_Fid_Geo(i,0)){
			plot_cc=true;
		}
		if(flags_->Flags::Plot_Fid_Geo(i,1)){
			plot_sc=true;
		}
		if(flags_->Flags::Plot_Fid_Geo(i,2)){
			plot_ec=true;
		}
	}
	if(!plot_geo){ return;}
	std::cout<<"Making Geometric Fiducial Plots\n";
	std::vector<long> space_dims(6);
	space_dims[5] = std::distance(std::begin(_species_), std::end(_species_)); //Particle Species
	space_dims[4] = 3; //CC, SC, EC
	space_dims[3] = std::distance(std::begin(_sector_), std::end(_sector_)); //Sector
	space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));//Cut and anti cut //Cut Status
	space_dims[1] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));//ID Cut
	space_dims[0] = 54+4;//total number of CC segment/side combinations. + 4  SC paddles go up to 18
	char hname[100];
	char hname2[100];

	CartesianGenerator cart(space_dims);
	TH2F_ptr_1d tmp_1d;
	TH2F_ptr_2d tmp_2d;
	TH2F_ptr_3d tmp_3d;
	TH2F_ptr_4d tmp_4d;
	TH2F_ptr_5d tmp_5d;
	TH2F_ptr_1d tmp2_1d;
	TH2F_ptr_2d tmp2_2d;
	TH2F_ptr_3d tmp2_3d;
	TH2F_ptr_4d tmp2_4d;
	TH2F_ptr_5d tmp2_5d;
	int species_idx = -1;
	int det_idx = -1;
	int sec_idx = -1;
	int cut_idx = -1;
	int pcut_idx = -1;
	int seg_idx = -1;
	char* pcut_name;
	bool pass = false;
	int cc_seg_idx = -1;
	int cc_side_idx = -1;
	while(cart.GetNextCombination()){
		pass = false;
		species_idx = cart[5];
		det_idx = cart[4];
		sec_idx = cart[3];
		cut_idx = cart[2];
		pcut_idx = cart[1];
		seg_idx = cart[0];
		if((pcut_idx == 0 && cut_idx == 2) || (pcut_idx >0 && cut_idx < 2)){
			if(_species_[species_idx]==_ele_){
				if(fun::pcut_perform(_species_[species_idx],_ecuts_[pcut_idx],flags_)){
					cc_seg_idx = Histogram::Geo_Fid_Paddle_Idx(det_idx,seg_idx);
					cc_side_idx = Histogram::Geo_Fid_Side_Idx(det_idx,seg_idx);
					if(detectors[det_idx+1]=="cc"){
						if(cc_seg_idx == _geo_fid_cc_segments_){
							if(cc_side_idx == _geo_fid_cc_sides_){
								sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s_Seg:All_Side:All",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx]);
								sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s_Seg:All_Side:All",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx]);
								pass = true;
							}else{
								sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s_Seg:All_Side:%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_cc_sides_[cc_side_idx]);
								sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s_Seg:All_Side:%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_cc_sides_[cc_side_idx]);
								pass = true;
							}
						}else{
							sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s_Seg:%d_Side:%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx],cc_seg_idx,_cc_sides_[cc_side_idx]);
							sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s_Seg:%d_Side:%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx],cc_seg_idx,_cc_sides_[cc_side_idx]);
							pass = true;
						}
					}else if(detectors[det_idx+1]=="sc"){
						if(cc_side_idx == 0){
							if(cc_seg_idx == _geo_fid_sc_paddles_){
								sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s_Paddle:All",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx]);
								sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s_Paddle:All",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx]);
							pass = true;	
							}else if(cc_seg_idx < _geo_fid_sc_paddles_){
								sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s_Paddle:%d",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx],seg_idx);
								sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s_Paddle:%d",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx],seg_idx);
								pass = true;
							}else{
								pass = false;
							}
						}else{
							pass = false;
						}
					}else if(detectors[det_idx+1]=="ec"){
						if(seg_idx == 0){
							sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx]);
							sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx]);
							pass = true;
						}else{
							pass == false;
						}
					}
					//sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_ecuts_[pcut_idx],_cut_[cut_idx]);
					//pass = true;
					//std::cout<<_Geo_Fid_hist.size() <<" " <<tmp_5d.size() <<" " <<tmp_4d.size() <<" " <<tmp_3d.size() <<" "<<tmp_2d.size() <<" "<<tmp_1d.size() <<" " <<hname <<"\n";
				}
			}else if(pcut_idx < std::distance(std::begin(_hcuts_), std::end(_hcuts_))){
				if(fun::pcut_perform(_species_[species_idx],_hcuts_[pcut_idx],flags_)){
					if(detectors[det_idx+1]=="cc"){
						if(seg_idx < _geo_fid_cc_segments_*_geo_fid_cc_sides_){
							for(int i=0; i<3; i++){
								if(seg_idx >= _geo_fid_cc_segments_*i){
									cc_side_idx = i;
									cc_seg_idx = seg_idx-(_geo_fid_cc_segments_*i);
								}
							}
							sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s_Seg:%d_Side:%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx],cc_seg_idx,_cc_sides_[cc_side_idx]);
							sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s_Seg:%d_Side:%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx],cc_seg_idx,_cc_sides_[cc_side_idx]);
							pass = true;
						}else if(seg_idx < _geo_fid_cc_segments_*_geo_fid_cc_sides_+3){
							cc_side_idx = seg_idx - (_geo_fid_cc_segments_*_geo_fid_cc_sides_);
							sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s_Seg:All_Side:%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_cc_sides_[cc_side_idx]);
							sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s_Seg:All_Side:%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_cc_sides_[cc_side_idx]);
							pass = true;
						}else{
							sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s_Seg:All_Side:All",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx]);
							sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s_Seg:All_Side:All",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx]);
							pass = true;
						}
					}else if(detectors[det_idx+1]=="sc"){
						if(seg_idx < _geo_fid_sc_paddles_){
							sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s_Paddle:%d",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx],seg_idx);
							sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s_Paddle:%d",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx],seg_idx);
							pass = true;
						}else if(seg_idx == _geo_fid_sc_paddles_){
							sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s_Paddle:All",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx]);
							sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s_Paddle:All",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx]);
							pass = true;
						}else{
							pass = false;
						}
					}else if(detectors[det_idx+1]=="ec"){
						if(seg_idx == 0){
							sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx]);
							sprintf(hname2,"%s_%s_Geometric_Fid2_Sec:%s_Cut:%s_%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx]);
							pass = true;
						}else{
							pass == false;
						}
					}
					//sprintf(hname,"%s_%s_Geometric_Fid_Sec:%s_Cut:%s_%s",_species_[species_idx],detectors[det_idx+1],_sector_[sec_idx],_hcuts_[pcut_idx],_cut_[cut_idx]);
					//pass = true;
					
					//std::cout<<_Geo_Fid_hist.size() <<" " <<tmp_4d.size() <<" " <<tmp_3d.size() <<" "<<tmp_2d.size() <<" "<<tmp_1d.size() <<" " <<hname <<"\n";
				}
			}
		}
		if(!plot_cc && detectors[det_idx+1]=="cc"){
			pass = false;
		}
		if(!plot_sc && detectors[det_idx+1]=="sc"){
			pass = false;
		}
		if(!plot_ec && detectors[det_idx+1]=="ec"){
			pass = false;
		}
		if(pass){
			//std::cout<<_Geo_Fid_hist.size() <<" " <<tmp_5d.size() <<" " <<tmp_4d.size() <<" " <<tmp_3d.size() <<" "<<tmp_2d.size() <<" "<<tmp_1d.size() <<" " <<hname <<"\n";
			tmp_1d.push_back(new TH2F(hname,hname,_geo_fid_xbin_[det_idx],_geo_fid_xmin_[det_idx],_geo_fid_xmax_[det_idx],_geo_fid_ybin_[det_idx],_geo_fid_ymin_[det_idx],_geo_fid_ymax_[det_idx]));
			tmp2_1d.push_back(new TH2F(hname2,hname2,_geo_fid2_xbin_[det_idx],_geo_fid2_xmin_[det_idx],_geo_fid2_xmax_[det_idx],_geo_fid2_ybin_[det_idx],_geo_fid2_ymin_[det_idx],_geo_fid2_ymax_[det_idx]));
		}
		if(cart[0] == space_dims[0]-1){
			if(tmp_1d.size()>0){
				tmp_2d.push_back(tmp_1d);
				tmp_1d.clear();
			}
			if(tmp2_1d.size()>0){
				tmp2_2d.push_back(tmp2_1d);
				tmp2_1d.clear();
			}
			if(cart[1] == space_dims[1]-1){
				if(tmp_2d.size()>0){
					tmp_3d.push_back(tmp_2d);
					tmp_2d.clear();
				}
				if(tmp2_2d.size()>0){
					tmp2_3d.push_back(tmp2_2d);
					tmp2_2d.clear();
				}
				if(cart[2] == space_dims[2]-1){
					if(tmp_3d.size()>0){
						tmp_4d.push_back(tmp_3d);
						tmp_3d.clear();
					}
					if(tmp2_3d.size()>0){
						tmp2_4d.push_back(tmp2_3d);
						tmp2_3d.clear();
					}
					if(cart[3] == space_dims[3]-1){
						if(tmp_4d.size()>0){
							tmp_5d.push_back(tmp_4d);
							tmp_4d.clear();
						}
						if(tmp2_4d.size()>0){
							tmp2_5d.push_back(tmp2_4d);
							tmp2_4d.clear();
						}
						if(cart[4] == space_dims[4]-1){
							if(tmp_5d.size()>0){
								_Geo_Fid_hist.push_back(tmp_5d);
								tmp_5d.clear();
							}
							if(tmp2_5d.size()>0){
								_Geo_Fid2_hist.push_back(tmp2_5d);
								tmp2_5d.clear();
							}
						}
					}
				}
			}
		}
	}
	std::cout<<"\tFinished Making Geometric Fiducial Histograms\n";
}
std::vector<int> Histogram::Geo_Fid_idx(const char* species_,const char* detector_, const char* sec_, const char* cut_, const char* pcut_, int det_seg_, int det_side_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	int species_idx = fun::species_idx(species_);//+fun::species_offset(species_,pcut_,flags_);
	idx.push_back(species_idx);
	if(detector_==_cc_){
		if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),0)){
			idx.push_back(0);
		}else{
			idx.push_back(-1);
		}
	}else if(detector_==_sc_){
		if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),0)){
			if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),1)){
				idx.push_back(1);
			}else{
				idx.push_back(-1);
			}
		}else{
			if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),1)){
				idx.push_back(0);
			}else{
				idx.push_back(-1);
			}
		}
	}else if(detector_==_ec_){
		if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),0)){
			if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),1)){
				if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),2)){
					idx.push_back(2);
				}else{
					idx.push_back(-1);
				}
			}else{
				if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),2)){
					idx.push_back(1);
				}else{
					idx.push_back(-1);
				}
			}
		}else{
			if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),1)){
				if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),2)){
					idx.push_back(1);
				}else{
					idx.push_back(-1);
				}
			}else{
				if(flags_->Flags::Plot_Fid_Geo(fun::species_idx(species_),2)){
					idx.push_back(0);
				}else{
					idx.push_back(-1);
				}
			}
		}
	}else{
		idx.push_back(-1);
	}
	idx.push_back(fun::sector_idx(sec_));
	idx.push_back(fun::cut_idx(cut_));
	
	if(species_==_ele_){
		if(fun::ecut_perform(pcut_,flags_)){
			if(pcut_==_none_ && cut_==_no_cut_){
				idx.push_back(fun::ecut_idx(pcut_)+fun::ecut_offset(pcut_,flags_));
			}else if(cut_!=_no_cut_){
				idx.push_back(fun::ecut_idx(pcut_)+fun::ecut_offset(pcut_,flags_)-1);
			}else{
				idx.push_back(-1);
			}
			//idx.push_back(fun::ecut_idx(pcut_)+fun::ecut_offset(pcut_,flags_));
		}else{
			idx.push_back(-1);
		}
	}else{
		if(fun::hcut_perform(species_, pcut_,flags_)){
			if(pcut_==_none_ && cut_==_no_cut_){
				idx.push_back(fun::hcut_idx(pcut_)+fun::hcut_offset(species_,pcut_,flags_));
			}else if(cut_!=_no_cut_){
				idx.push_back(fun::hcut_idx(pcut_)+fun::hcut_offset(species_,pcut_,flags_)-1);
			}else{
				idx.push_back(-1);
			}
			//idx.push_back(fun::hcut_idx(pcut_)+fun::hcut_offset(species_,pcut_,flags_));
		}else{
			idx.push_back(-1);
		}
	}
	if(detector_=="cc"){
		idx.push_back(Histogram::Geo_Fid_Seg_Idx(0, det_seg_, det_side_));
	}else if(detector_=="sc"){
		idx.push_back(Histogram::Geo_Fid_Seg_Idx(1, det_seg_, det_side_));
	}else if(detector_=="ec"){
		idx.push_back(Histogram::Geo_Fid_Seg_Idx(2, det_seg_, det_side_));
	}else{
		idx.push_back(-1);
	}
	return idx;
}

void Histogram::Geo_Fid_Fill(float x_, float y_, float weight_, const char* species_, const char* detector_, const char* sec_, const char* cut_, const char* pcut_, int det_seg_, int det_side_, std::shared_ptr<Flags> flags_){
	//std::cout<<"Attempting to fill Geo Fid\n";
	//if(!flags_->CC_Geo() && !flags_->SC_Geo() && !flags_->EC_Geo()){ return;}
	int sp_idx = fun::species_idx(species_);
	int det_idx = fun::geo_det_idx(detector_);

	if(!flags_->Flags::Plot_Fid_Geo(sp_idx,det_idx)){ return;}
	if(detector_ == _sc_ && det_seg_==-1){return;}
	if(detector_ == _cc_ && det_seg_==-1){return;}
	if(isnan(x_) || isnan(y_)){return;}
	if(x_==0){return;}
	
	std::vector<int> idx = Histogram::Geo_Fid_idx(species_,detector_,sec_,cut_,pcut_,det_seg_,det_side_,flags_);
	
	
	//std::cout<<"Filling Fid Geo " <<species_ <<" " <<detector_ <<" with: ";
	//std::cout<<"\tx:" <<x_ <<" y:" <<y_ <<" weight:" <<weight_ <<" sec:" <<sec_ <<" cut:" <<cut_ <<" pcut_:" <<pcut_ <<" det_seg:" <<det_seg_ <<" det_side:" <<det_side_ <<"\n";
	//fun::print_vector_idx(idx);
	//if(fun::check_sec(sec_,x_,y_) || (detector_ == _sc_ && fun::check_sec(_sec1_,x_,y_))){ 
	//std::cout<<"x:" <<x_ <<" " <<isnan(x_) <<" y:" <<y_ <<" " <<isnan(y_) <<"\n";
	if(!isnan(x_) && !isnan(y_)){

		float x_cent = detect::x_det_center(x_,y_,fun::sector_idx(sec_)+1);
		float y_cent = detect::y_det_center(x_,y_,fun::sector_idx(sec_)+1);
		if(detector_ == _sc_){//All the x-y coordinates for SC are sector centered, but not in the way I want. 
			x_cent = detect::x_det_center(x_,y_,fun::sector_idx(_sec1_)+1);
			y_cent = detect::y_det_center(x_,y_,fun::sector_idx(_sec1_)+1);
			//if(species_ == _pim_){
		//		std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" paddle:" <<det_seg_ <<" side:" <<det_side_ <<" x:" <<x_ <<" y:" <<y_ <<" x_c:" <<x_cent <<" y_c:" <<y_cent <<"\n";
		//		fun::print_vector_idx(idx);
		//	}
		}
		//std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" paddle:" <<det_seg_ <<" side:" <<det_side_ <<" x:" <<x_ <<" y:" <<y_ <<" x_c:" <<x_cent <<" y_c:" <<y_cent <<"\n";
		//fun::print_vector_idx(idx);
		
		if(Histogram::OK_Idx(idx)){
			//std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" paddle:" <<det_seg_ <<" side:" <<det_side_ <<" x:" <<x_ <<" y:" <<y_ <<" x_c:" <<x_cent <<" y_c:" <<y_cent <<"\n";
			//fun::print_vector_idx(idx);
			_Geo_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(x_cent,y_cent,weight_);//Sector
			_Geo_Fid2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(x_,y_,weight_);//Sector
			idx = Histogram::Geo_Fid_idx(species_,detector_,_sec_all_,cut_,pcut_,det_seg_,det_side_,flags_);
			//fun::print_vector_idx(idx);
			_Geo_Fid_hist[idx[0]][idx[1]][6][idx[3]][idx[4]][idx[5]]->Fill(x_cent,y_cent,weight_);//All Sectors
			_Geo_Fid2_hist[idx[0]][idx[1]][6][idx[3]][idx[4]][idx[5]]->Fill(x_,y_,weight_);//All Sectors
		}
		if(detector_ == _cc_){
			//Sum over All Segments, split by side
			std::vector<int> idx2 = Histogram::Geo_Fid_idx(species_,detector_,sec_,cut_,pcut_,_num_cc_segments_,det_side_,flags_);
			//std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" x:" <<x_ <<" y:" <<y_ <<"\n";
			//fun::print_vector_idx(idx);
			if(Histogram::OK_Idx(idx2)){
				_Geo_Fid_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]][idx2[4]][idx2[5]]->Fill(x_cent,y_cent,weight_);//Sector
				_Geo_Fid_hist[idx2[0]][idx2[1]][6][idx2[3]][idx2[4]][idx2[5]]->Fill(x_cent,y_cent,weight_);//All Sectors
				_Geo_Fid2_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]][idx2[4]][idx2[5]]->Fill(x_,y_,weight_);//Sector
				_Geo_Fid2_hist[idx2[0]][idx2[1]][6][idx2[3]][idx2[4]][idx2[5]]->Fill(x_,y_,weight_);//All Sectors
			}
			//Sum over all segments, add sides together
			std::vector<int> idx3 = Histogram::Geo_Fid_idx(species_,detector_,sec_,cut_,pcut_,_num_cc_segments_,3,flags_);
			//std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" x:" <<x_ <<" y:" <<y_ <<"\n";
			//fun::print_vector_idx(idx);
			if(Histogram::OK_Idx(idx3)){
				_Geo_Fid_hist[idx3[0]][idx3[1]][idx3[2]][idx3[3]][idx3[4]][idx3[5]]->Fill(x_cent,y_cent,weight_);//Sector
				_Geo_Fid_hist[idx3[0]][idx3[1]][6][idx3[3]][idx3[4]][idx3[5]]->Fill(x_cent,y_cent,weight_);//All Sectors
				_Geo_Fid2_hist[idx3[0]][idx3[1]][idx3[2]][idx3[3]][idx3[4]][idx3[5]]->Fill(x_,y_,weight_);//Sector
				_Geo_Fid2_hist[idx3[0]][idx3[1]][6][idx3[3]][idx3[4]][idx3[5]]->Fill(x_,y_,weight_);//All Sectors
			}
		}else if(detector_ == _sc_){
			//Sum over all Paddles
			std::vector<int> idx4 = Histogram::Geo_Fid_idx(species_,detector_,sec_,cut_,pcut_,_num_sc_paddles_,det_side_,flags_);
			//std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" x:" <<x_ <<" y:" <<y_ <<"\n";
			//fun::print_vector_idx(idx);
			if(Histogram::OK_Idx(idx4)){
				_Geo_Fid_hist[idx4[0]][idx4[1]][idx4[2]][idx4[3]][idx4[4]][idx4[5]]->Fill(x_cent,y_cent,weight_);//Sector
				_Geo_Fid_hist[idx4[0]][idx4[1]][6][idx4[3]][idx4[4]][idx4[5]]->Fill(x_cent,y_cent,weight_);//All Sectors
				_Geo_Fid2_hist[idx4[0]][idx4[1]][idx4[2]][idx4[3]][idx4[4]][idx4[5]]->Fill(x_,y_,weight_);//Sector
				_Geo_Fid2_hist[idx4[0]][idx4[1]][6][idx4[3]][idx4[4]][idx4[5]]->Fill(x_,y_,weight_);//All Sectors
			}
		}/*else if(detector_ == _ec_){
			std::vector<int> idx5 = Histogram::Geo_Fid_idx(species_,detector_,sec_,cut_,pcut_,0,0,flags_);
			//std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" x:" <<x_ <<" y:" <<y_ <<"\n";
			//fun::print_vector_idx(idx);
			if(Histogram::OK_Idx(idx5)){
				_Geo_Fid_hist[idx5[0]][idx5[1]][idx5[2]][idx5[3]][idx5[4]][idx5[5]]->Fill(x_cent,y_cent,weight_);//Sector
				_Geo_Fid_hist[idx5[0]][idx5[1]][6][idx5[3]][idx5[4]][idx5[5]]->Fill(x_cent,y_cent,weight_);//All Sectors
			}
		}*/
		/*if(Histogram::OK_Idx(idx)){
			_Geo_Fid2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(x_,y_,weight_);//Sector
			_Geo_Fid2_hist[idx[0]][idx[1]][6][idx[3]][idx[4]][idx[5]]->Fill(x_,y_,weight_);//All Sectors
		}
		//Full range histograms
		if(detector_ == _cc_){
			//Sum over All Segments, split by side
			std::vector<int> idx22 = Histogram::Geo_Fid_idx(species_,detector_,sec_,cut_,pcut_,_num_cc_segments_,det_side_,flags_);
			//std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" x:" <<x_ <<" y:" <<y_ <<"\n";
			//fun::print_vector_idx(idx);
			if(Histogram::OK_Idx(idx22)){
				_Geo_Fid2_hist[idx22[0]][idx22[1]][idx22[2]][idx22[3]][idx22[4]][idx22[5]]->Fill(x_,y_,weight_);//Sector
				_Geo_Fid2_hist[idx22[0]][idx22[1]][6][idx22[3]][idx22[4]][idx22[5]]->Fill(x_,y_,weight_);//All Sectors
			}
			//Sum over all segments, add sides together
			std::vector<int> idx23 = Histogram::Geo_Fid_idx(species_,detector_,sec_,cut_,pcut_,_num_cc_segments_,3,flags_);
			//std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" x:" <<x_ <<" y:" <<y_ <<"\n";
			//fun::print_vector_idx(idx);
			if(Histogram::OK_Idx(idx23)){
				_Geo_Fid2_hist[idx23[0]][idx23[1]][idx23[2]][idx23[3]][idx23[4]][idx23[5]]->Fill(x_,y_,weight_);//Sector
				_Geo_Fid2_hist[idx23[0]][idx23[1]][6][idx23[3]][idx23[4]][idx23[5]]->Fill(x_,y_,weight_);//All Sectors
			}
		}else if(detector_ == _sc_){
			//Sum over all Paddles
			std::vector<int> idx24 = Histogram::Geo_Fid_idx(species_,detector_,sec_,cut_,pcut_,_num_sc_paddles_,det_side_,flags_);
			//std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" x:" <<x_ <<" y:" <<y_ <<"\n";
			//fun::print_vector_idx(idx);
			if(Histogram::OK_Idx(idx24)){
				_Geo_Fid2_hist[idx24[0]][idx24[1]][idx24[2]][idx24[3]][idx24[4]][idx24[5]]->Fill(x_,y_,weight_);//Sector
				_Geo_Fid2_hist[idx24[0]][idx24[1]][6][idx24[3]][idx24[4]][idx24[5]]->Fill(x_,y_,weight_);//All Sectors
			}
		}else if(detector_ == _ec_){
			//Sum over all Paddles
			std::vector<int> idx25 = Histogram::Geo_Fid_idx(species_,detector_,sec_,cut_,pcut_,_num_sc_paddles_,det_side_,flags_);
			std::cout<<"Filling " <<species_ <<" " <<detector_ <<" " <<sec_ <<" " <<cut_ <<" " <<pcut_ <<" x:" <<x_ <<" y:" <<y_ <<"\n";
			fun::print_vector_idx(idx25);
			if(Histogram::OK_Idx(idx25)){
				_Geo_Fid2_hist[idx25[0]][idx25[1]][idx25[2]][idx25[3]][idx25[4]][idx25[5]]->Fill(x_,y_,weight_);//Sector
				idx25 = Histogram::Geo_Fid_idx(species_,detector_,_sec_all_,cut_,pcut_,_num_sc_paddles_,det_side_,flags_);
				fun::print_vector_idx(idx25);
				_Geo_Fid2_hist[idx25[0]][idx25[1]][6][idx25[3]][idx25[4]][idx25[5]]->Fill(x_,y_,weight_);//All Sectors
			}
		}*/
	}
}

void Histogram::Geo_Fid_Write(std::shared_ptr<Flags> flags_){
	//std::cout<<"Trying to Write Geometric Fiducial Plots\n";
	bool plot_geo = false;
	bool plot_cc[4] = {false,false,false,false};
	bool plot_sc[4] = {false,false,false,false};
	bool plot_ec[4] = {false,false,false,false};
	bool plot_par[4] = {false,false,false,false};
	for(int i=0; i<4; i++){
		if(flags_->Flags::Plot_Fid_Geo(i,0)){
			plot_cc[i] = true;
		}
		if(flags_->Flags::Plot_Fid_Geo(i,1)){
			plot_sc[i] = true;
		}
		if(flags_->Flags::Plot_Fid_Geo(i,2)){
			plot_ec[i] = true;
		}
		for(int j=0; j<3; j++){
			if(flags_->Flags::Plot_Fid_Geo(i,j)){
				plot_geo=true;
				plot_par[i] = true;
			}
		}
	}
	if(!plot_geo){ return;}
	std::cout<<"Writing Geometric Fiducial Plots\n";
	std::vector<long> space_dims(6);
	space_dims[5] = std::distance(std::begin(_species_), std::end(_species_)); //Particle Species
	space_dims[4] = 3; //CC, SC, EC
	space_dims[3] = std::distance(std::begin(_sector_), std::end(_sector_)); //Sector
	space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));//Cut and anti cut //Cut Status
	space_dims[1] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));//ID Cut
	space_dims[0] = 54+4;//total number of CC segment/side combinations. + 4  SC paddles go up to 18

	char dirname[100];
	std::cout<<"\tMaking TDirectories\n";
	TDirectory* dir_geo_fid= _RootOutputFile->mkdir("Geometric Distributions");
	//std::cout<<" Done\n";
	dir_geo_fid->cd();
	//std::cout<<"Making Sub Directories: ";
	TDirectory* dir_geo_fid_sub[4][3+1][7+1][std::distance(std::begin(_ecuts_), std::end(_ecuts_))+1];//{sector,proton_thesh}
	//Histogram::W_bot(Wbin),Histogram::W_top(Wbin),Histogram::Q2_bot(Q2bin),Histogram::Q2_top(Q2bin),Histogram::Xij_Bin_Min(j,Xij,Wbin,var),Histogram::Xij_Bin_Max(j,Xij,Wbin,var)
	for(int i=0; i<4; i++){
		if(plot_par[i]){
			sprintf(dirname,"Geo Fid %s",species[i]);
			dir_geo_fid_sub[i][0][0][0] = dir_geo_fid->mkdir(dirname);
			//std::cout<<"\tMade " <<i <<" " <<0 <<" " <<0 <<" " <<0 <<"\n";
			for(int j=0; j<3; j++){
				if(flags_->Plot_Fid_Geo(i,j)){
					sprintf(dirname,"Geo Fid %s %s",species[i],detectors[j+1]);
					//std::cout<<"\tMade " <<i <<" " <<j+1 <<" " <<0 <<" " <<0 <<"\n";
					dir_geo_fid_sub[i][j+1][0][0] = dir_geo_fid_sub[i][0][0][0]->mkdir(dirname);
					for(int k=0; k<7; k++){
						sprintf(dirname,"Geo Fid %s %s %s",species[i],detectors[j+1],_sector_[k]);
						dir_geo_fid_sub[i][j+1][k+1][0] = dir_geo_fid_sub[i][j+1][0][0]->mkdir(dirname);
						//std::cout<<"\tMade " <<i <<" " <<j+1 <<" " <<k+1 <<" " <<0 <<"\n";
						if(species[i]==_ele_){
							for(int l=0; l<std::distance(std::begin(_ecuts_), std::end(_ecuts_)); l++){
								if(fun::ecut_perform(_ecuts_[l],flags_)){
									sprintf(dirname,"Geo Fid %s %s %s %s",species[i],detectors[j+1],_sector_[k],_ecuts_[l]);
									//std::cout<<"\tMade " <<i <<" " <<j+1 <<" " <<k+1 <<" " <<l+1 <<"\t" <<species[i] <<" " <<detectors[j+1] <<" " <<_sector_[k] <<" " <<_ecuts_[l] <<"\n";
									dir_geo_fid_sub[i][j+1][k+1][l+1] = dir_geo_fid_sub[i][j+1][k+1][0]->mkdir(dirname);
								}
							}
						}else{
							for(int l=0; l<std::distance(std::begin(_hcuts_), std::end(_hcuts_)); l++){
								if(fun::hcut_perform(species[i],_hcuts_[l],flags_)){
									sprintf(dirname,"Geo Fid %s %s %s %s",species[i],detectors[j+1],_sector_[k],_hcuts_[l]);
									//std::cout<<"\tMade " <<i <<" " <<j+1 <<" " <<k+1 <<" " <<l+1 <<"\t" <<species[i] <<" " <<detectors[j+1] <<" " <<_sector_[k] <<" " <<_hcuts_[l] <<"\n";
									dir_geo_fid_sub[i][j+1][k+1][l+1] = dir_geo_fid_sub[i][j+1][k+1][0]->mkdir(dirname);
								}
							}
						}
					}
				}
			}
		}
	}
	std::cout<<"Directories Made\n";


	CartesianGenerator cart(space_dims);
	int species_idx = -1;
	int det_idx = -1;
	int sec_idx = -1;
	int cut_idx = -1;
	int pcut_idx = -1;
	int seg_idx = -1;
	int pad_idx = -1;
	int side_idx = -1;
	char* pcut_name;
	bool pass = false;
	std::vector<int> hist_idx;
	while(cart.GetNextCombination()){
		pass = false;
		species_idx = cart[5];
		det_idx = cart[4];
		sec_idx = cart[3];
		cut_idx = cart[2];
		pcut_idx = cart[1];
		seg_idx = cart[0];
		//std::cout<<"********** " <<species_idx <<" " <<det_idx <<" " <<sec_idx <<" " <<cut_idx <<" " <<pcut_idx <<"\n";
		pad_idx = Histogram::Geo_Fid_Paddle_Idx(det_idx,seg_idx);
		side_idx = Histogram::Geo_Fid_Side_Idx(det_idx,seg_idx);
		if(flags_->Flags::Plot_Fid_Geo(species_idx,det_idx)){
			if(species[species_idx]==_ele_){
				hist_idx = Histogram::Geo_Fid_idx(species[species_idx],detectors[det_idx+1],_sector_[sec_idx],_cut_[cut_idx],_ecuts_[pcut_idx],pad_idx,side_idx,flags_);
				//std::cout<<"\t" <<species[species_idx] <<" " <<detectors[det_idx+1] <<" " <<_sector_[sec_idx] <<" " <<_cut_[cut_idx] <<" " <<_ecuts_[pcut_idx] <<"\n";
			}else if(pcut_idx < std::distance(std::begin(_hcuts_), std::end(_hcuts_))){
				hist_idx = Histogram::Geo_Fid_idx(species[species_idx],detectors[det_idx+1],_sector_[sec_idx],_cut_[cut_idx],_hcuts_[pcut_idx],pad_idx,side_idx,flags_);
				//std::cout<<"\t" <<species[species_idx] <<" " <<detectors[det_idx+1] <<" " <<_sector_[sec_idx] <<" " <<_cut_[cut_idx] <<" " <<_hcuts_[pcut_idx] <<"\n";
			}else{
				//hist_idx = Histogram::Geo_Fid_idx("no","no","no","no","no",pad_idx,side_idx,flags_);
				hist_idx = Histogram::Geo_Fid_idx("no","no","no","no","no",-1,-1,flags_);
			}
			if(Histogram::OK_Idx(hist_idx)){
				//std::cout<<"\tcd()-> " <<species_idx <<" " <<det_idx+1 <<" " <<sec_idx+1 <<" " <<pcut_idx+1 <<"\n";
				if(species_idx == 0){
					//std::cout<<"\tcd()-> " <<species[species_idx] <<" " <<detectors[det_idx+1] <<" " <<_sector_[sec_idx] <<" " <<_ecuts_[pcut_idx] <<"\n";
				}else{
					//std::cout<<"\tcd()-> " <<species[species_idx] <<" " <<detectors[det_idx+1] <<" " <<_sector_[sec_idx] <<" " <<_ecuts_[pcut_idx] <<"\n";
				}
				//fun::print_vector_idx(hist_idx);
				dir_geo_fid_sub[species_idx][det_idx+1][sec_idx+1][pcut_idx+1]->cd();
				//std::cout<<"\t" <<_Geo_Fid_hist[hist_idx[0]][hist_idx[1]][hist_idx[2]][hist_idx[3]][hist_idx[4]]->GetName() <<"\n";
				_Geo_Fid_hist[hist_idx[0]][hist_idx[1]][hist_idx[2]][hist_idx[3]][hist_idx[4]][hist_idx[5]]->Write();
				_Geo_Fid2_hist[hist_idx[0]][hist_idx[1]][hist_idx[2]][hist_idx[3]][hist_idx[4]][hist_idx[5]]->Write();
			}
		}
	}

}

//*------------------------------End Detector Geometric Cut Peak------------------*

//*------------------------------Start Bin Centering Corrections for Polarization------------------*
	//double Xij_Bin_Min(int bin_, int Xij_, int W_bin_, int var_set_);
	//double Xij_Bin_Max(int bin_, int Xij_, int W_bin_, int var_set_);
	void Histogram::Bin_Centering2_Make(std::shared_ptr<Flags> flags_){
	//[W][Q2][Var Set][Xij][Bin Xij][W,Q2,Xij projection]
		if(!flags_->Plot_Bin_Centering2()){return;}
		std::cout<<"Making Bin Centering Correction Histograms\n";
		std::vector<long> space_dims(5);
		space_dims[3] = 29; //W Bins
		space_dims[2] = 5; //Q2 Bins
		space_dims[1] = 3; //Var Set
		space_dims[0] = 4; //Xij
		space_dims[4] = 4;//variable
		CartesianGenerator cart(space_dims);
		char hname[200];
		TH1D_ptr_1d hist_1d;
		TH1D_ptr_2d hist_2d;
		TH1D_ptr_3d hist_3d;
		TH1D_ptr_4d hist_4d;
		TH1D_ptr_5d hist_5d;
		double bot;
		double top;
		int num_bins_xij[5] = {_MM_bins_+2*_MM_wider_,_MM_bins_+2*_MM_wider_,_theta_bins_,_alpha_bins_,_phi_bins_};
		int num_bins_fin = 0;
		int Wbin = 0;
		int Q2bin = 0;
		int Xij=0;
		int var=0;
		char proj_name[100];
		char* xij_names[] = {"MM1","MM2","theta","alpha","phi"};
		while(cart.GetNextCombination()){
			Wbin = cart[3];
			Q2bin=cart[2];
			Xij=cart[0];
			var=cart[1];
			for(int k=0; k<_phi_bins_; k++){
			//for(int k=0; k<3; k++){
				for(int j=0; j<num_bins_xij[Xij]; j++){
					switch(cart[4]){
						case 0: 
							num_bins_fin = 29;//W
							sprintf(proj_name,"W");
							bot = Histogram::W_bot(Wbin);
							top = Histogram::W_top(Wbin);
							//std::cout<<cart[4] <<" " <<cart[3] <<" " <<cart[2] <<" " <<cart[1] <<" " <<cart[0] <<" " <<<<"\n";
							//std::cout<<"Making BC2 Hist1: " <<_Bin_Center2_hist1.size() <<" " <<hist_5d.size() <<" " <<hist_4d.size() <<" " <<hist_3d.size() <<" " <<hist_2d.size() <<" " <<hist_1d.size() <<"\n";
						break;
						case 1: 
							num_bins_fin = 5;//Q2
							sprintf(proj_name,"Q2");
							bot = Histogram::Q2_bot(Q2bin);
							top = Histogram::Q2_top(Q2bin);
							//std::cout<<"Making BC2 Hist2: " <<_Bin_Center2_hist2.size() <<" " <<hist_5d.size() <<" " <<hist_4d.size() <<" " <<hist_3d.size() <<" " <<hist_2d.size() <<" " <<hist_1d.size() <<"\n";
						break;
						case 2: 
							num_bins_fin = num_bins_xij[Xij];//Xij
							sprintf(proj_name,"%s_%s",_var_names_[var],xij_names[Xij]);
							bot = Histogram::Xij_Bin_Min(j,Xij,Wbin,var);
							top = Histogram::Xij_Bin_Max(j,Xij,Wbin,var);
							//std::cout<<"Making BC2 Hist3: " <<_Bin_Center2_hist3.size() <<" " <<hist_5d.size() <<" " <<hist_4d.size() <<" " <<hist_3d.size() <<" " <<hist_2d.size() <<" " <<hist_1d.size() <<"\n";
						break;
						case 3: 
							num_bins_fin = num_bins_xij[4];//phi
							sprintf(proj_name,"%s_%s",_var_names_[var],xij_names[4]);
							bot = Histogram::Xij_Bin_Min(k,4,Wbin,var);
							top = Histogram::Xij_Bin_Max(k,4,Wbin,var);
							//std::cout<<"Making BC2 Hist4: " <<_Bin_Center2_hist4.size() <<" " <<hist_5d.size() <<" " <<hist_4d.size() <<" " <<hist_3d.size() <<" " <<hist_2d.size() <<" " <<hist_1d.size() <<"\n";
						break;
						default:
							std::cout<<"Improper variable choice for bin centering\n";
						break;
					}
					sprintf(hname,"%s_Bin_Center_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_%s_%s:[%.4f,%.4f)_Phi:[%.4f,%.4f)",proj_name,Histogram::W_bot(Wbin),Histogram::W_top(Wbin),Histogram::Q2_bot(Q2bin),Histogram::Q2_top(Q2bin),_var_names_[var],xij_names[Xij],Histogram::Xij_Bin_Min(j,Xij,Wbin,var),Histogram::Xij_Bin_Max(j,Xij,Wbin,var),Histogram::Xij_Bin_Min(k,4,Wbin,var),Histogram::Xij_Bin_Max(k,4,Wbin,var));
					hist_1d.push_back(new TH1D(hname,hname,11,bot,top));
				}
				if(hist_1d.size()>0){
					hist_2d.push_back(hist_1d);
					hist_1d.clear();
				}
			}
			if(hist_2d.size()>0){
				hist_3d.push_back(hist_2d);
				hist_2d.clear();
			}
			if(cart[0] == space_dims[0]-1){
				if(hist_3d.size()>0){
					hist_4d.push_back(hist_3d);
					hist_3d.clear();

					//7 1 0 0 5 6 
				}
				if(cart[1] == space_dims[1]-1){
					if(hist_4d.size()>0){
						hist_5d.push_back(hist_4d);
						hist_4d.clear();
					}
					if(cart[2] == space_dims[2]-1){
						if(hist_5d.size()>0){
							if(cart[4]==0){
								_Bin_Center2_hist1.push_back(hist_5d);
							}else if(cart[4]==1){
								_Bin_Center2_hist2.push_back(hist_5d);
							}else if(cart[4]==2){
								_Bin_Center2_hist3.push_back(hist_5d);
							}else if(cart[4]==3){
								_Bin_Center2_hist4.push_back(hist_5d);
							}
							//_Bin_Center2_hist.push_back(hist_5d);
							hist_5d.clear();
						}
					}
				}
			}
		}
		std::cout<<"\tFinished Making Bin Centering2 Histograms\nbop bop bop\n";
	}

	std::vector<int> Histogram::Bin_Centering2_idx(float W_, float Q2_, int var_set_, int Xij_, double Xij_val_, double phi_val_, int variable_, std::shared_ptr<Flags> flags_){
		std::vector<int> idx;
		idx.push_back(Histogram::W_bin(W_));
		idx.push_back(Histogram::Q2_bin(Q2_));
		idx.push_back(var_set_);
		idx.push_back(Xij_);
		idx.push_back(Histogram::Friend_phi_idx(phi_val_));
		switch(Xij_){
			case 0:
				idx.push_back(Histogram::Friend_MM_idx(Xij_val_,var_set_,Histogram::W_bin(W_)));
			break;
			case 1:
				idx.push_back(Histogram::Friend_MM2_idx(Xij_val_,var_set_,Histogram::W_bin(W_)));
			break;
			case 2:
				idx.push_back(Histogram::Friend_theta_idx(Xij_val_));
			break;
			case 3:
				idx.push_back(Histogram::Friend_alpha_idx(Xij_val_));
			break;
			case 4:
				idx.push_back(Histogram::Friend_phi_idx(Xij_val_));
			break;
			default:
				idx.push_back(-1);
			break;
		}
		//idx.push_back(variable_);
		return idx;
	}
	void Histogram::Bin_Centering2_Fill(float W_, float Q2_, int var_set_, int Xij_, double Xij_val_, double phi_val_, int variable_, double weight_, std::shared_ptr<Flags> flags_){
		if(!flags_->Plot_Bin_Centering2()){return;}
		std::vector<int> idx = Histogram::Bin_Centering2_idx(W_,Q2_,var_set_,Xij_,Xij_val_,phi_val_,variable_,flags_);
		//std::cout<<"Filling Bin Centering Corrections2 with W:" <<W_ <<" Q2:" <<Q2_ <<" var_set:" <<var_set_ <<" Xij:" <<Xij_ <<" Xij_val:" <<Xij_val_ <<" phi_val:" <<phi_val_ <<" variable:" <<variable_ <<" weight:" <<weight_ <<"\n";
		//fun::print_vector_idx(idx);
		if(Histogram::OK_Idx(idx)){
			switch(variable_){
				case 0:
					_Bin_Center2_hist1[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(W_,weight_);
				break;
				case 1:
					_Bin_Center2_hist2[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(Q2_,weight_);
				break;
				case 2:
					_Bin_Center2_hist3[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(Xij_val_,weight_);
				break;
				case 3:
					_Bin_Center2_hist4[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(phi_val_,weight_);
				break;
				default:
					std::cout<<"Improper Variable designation for Bin Centering Corrections\n";
				break;
			}
		}
	}
	void Histogram::Bin_Centering2_Write(std::shared_ptr<Flags> flags_){
		if(!flags_->Plot_Bin_Centering2()){ return; }
		std::cout<<"Writing Bin Centering Plots for Polarization\n";
		char dir_name[100];
		//std::cout<<"Making Large Directory:";
		TDirectory* dir_bc2 = _RootOutputFile->mkdir("Bin Centering2 Plots");
		//std::cout<<" Done\n";
		dir_bc2->cd();
		//std::cout<<"Making Sub Directories: ";
		TDirectory* dir_bc2_sub[29][5+1][3+1][4+1][4+1];//[W][Q2][var set][xij][W,Q2,Xij,phi]
		//Histogram::W_bot(Wbin),Histogram::W_top(Wbin),Histogram::Q2_bot(Q2bin),Histogram::Q2_top(Q2bin),Histogram::Xij_Bin_Min(j,Xij,Wbin,var),Histogram::Xij_Bin_Max(j,Xij,Wbin,var)
		char* xij_names[] = {"MM1","MM2","theta","alpha","phi"};
		char proj_name[100];
		for(int i=0; i<29; i++){//W
			sprintf(dir_name,"Bin Centering2 W:%.3f-%.3f",Histogram::W_bot(i),Histogram::W_top(i));
			dir_bc2_sub[i][0][0][0][0] = dir_bc2->mkdir(dir_name);
			dir_bc2_sub[i][0][0][0][0]->cd();
			for(int j=0; j<5; j++){//Q2
				sprintf(dir_name,"Bin Centering2 W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_bot(i),Histogram::W_top(i),Histogram::Q2_bot(j),Histogram::Q2_top(j));
				dir_bc2_sub[i][j+1][0][0][0] = dir_bc2_sub[i][0][0][0][0]->mkdir(dir_name);
				dir_bc2_sub[i][j+1][0][0][0]->cd();
				for(int k=0; k<3; k++){//var set
					sprintf(dir_name,"Bin Centering2 W:%.3f-%.3f Q2:%.2f-%.2f var:%s",Histogram::W_bot(i),Histogram::W_top(i),Histogram::Q2_bot(j),Histogram::Q2_top(j),_var_names_[k]);
					dir_bc2_sub[i][j+1][k+1][0][0] = dir_bc2_sub[i][j+1][0][0][0]->mkdir(dir_name);
					dir_bc2_sub[i][j+1][k+1][0][0]->cd();
					for(int l=0; l<4; l++){//Xij
						sprintf(dir_name,"Bin Centering2 W:%.3f-%.3f Q2:%.2f-%.2f var:%s Xij:%s",Histogram::W_bot(i),Histogram::W_top(i),Histogram::Q2_bot(j),Histogram::Q2_top(j),_var_names_[k],xij_names[l]);
						dir_bc2_sub[i][j+1][k+1][l+1][0] = dir_bc2_sub[i][j+1][k+1][0][0]->mkdir(dir_name);
						dir_bc2_sub[i][j+1][k+1][l+1][0]->cd();
						for(int m=0; m<4; m++){//{W,Q2,Xij,phi}
							switch(m){
								case 0:
									sprintf(proj_name,"W");
								break;
								case 1:
									sprintf(proj_name,"Q2");
								break;
								case 2:
									sprintf(proj_name,"%s",xij_names[l]);
								break;
								case 3:
									sprintf(proj_name,"%s",xij_names[4]);
								break;
								default:
									sprintf(proj_name,"wrong");
								break;
							}
							sprintf(dir_name,"Bin Centering W:%.3f-%.3f Q2:%.2f-%.2f var:%s Xij:%s proj:%s",Histogram::W_bot(i),Histogram::W_top(i),Histogram::Q2_bot(j),Histogram::Q2_top(j),_var_names_[k],xij_names[l],proj_name);
							dir_bc2_sub[i][j+1][k+1][l+1][m+1] = dir_bc2_sub[i][j+1][k+1][l+1][0]->mkdir(dir_name);
						}
					}
				}
			}
		}
		std::vector<long> space_dims(5);
		space_dims[4] = 29; //W Bins
		space_dims[3] = 5; //Q2 Bins
		space_dims[2] = 3; //Var Set
		space_dims[1] = 4; //Xij
		space_dims[0] = 4; //variable
		CartesianGenerator cart(space_dims);
		char hname[100];
		std::vector<int> idx;
		int Wbin, Q2bin, Xij, var;
		double val, Xij_val, phi_val;
		int num_bins_xij[5] = {_MM_bins_+2*_MM_wider_,_MM_bins_+2*_MM_wider_,_theta_bins_,_alpha_bins_,_phi_bins_};
		char * xij_units[5] = {"GeV","GeV","Degrees","Degrees","Degrees"};
		double top,bot;
		char xlabel[100];
		while(cart.GetNextCombination()){
			Wbin = cart[4];
			Q2bin=cart[3];
			Xij=cart[1];
			var=cart[2];
			//std::cout<<"W bin:" <<Wbin <<" Q2 bin:" <<Q2bin <<" var:" <<var <<" Xij:" <<Xij <<" going into loop\n";
			for(int k=0; k<num_bins_xij[4]; k++){
				phi_val = (Xij_Bin_Min(k,4,Wbin,var)+Xij_Bin_Max(k,4,Wbin,var))/2.0;
				for(int j=0; j<num_bins_xij[Xij]; j++){
					Xij_val = (Xij_Bin_Min(j,Xij,Wbin,var)+Xij_Bin_Max(j,Xij,Wbin,var))/2.0;
				
					switch(cart[0]){
						case 0: 
							bot = Histogram::W_bot(Wbin);
							top = Histogram::W_top(Wbin);
							sprintf(xlabel,"W (GeV)");
						break;
						case 1: 
							bot = Histogram::Q2_bot(Q2bin);
							top = Histogram::Q2_top(Q2bin);
							sprintf(xlabel,"Q2 (GeV^2)");
						break;
						case 2: 
							bot = Histogram::Xij_Bin_Min(j,Xij,Wbin,var);
							top = Histogram::Xij_Bin_Max(j,Xij,Wbin,var);
							sprintf(xlabel,"%s %s",xij_names[Xij],xij_units[Xij]);
						break;
						case 3: 
							bot = Histogram::Xij_Bin_Min(k,4,Wbin,var);
							top = Histogram::Xij_Bin_Max(k,4,Wbin,var);
							sprintf(xlabel,"%s %s",xij_names[4],xij_units[4]);
						break;
					}
					val = (top+bot)/2.0;
					//std::cout<<"W: " <<W_center(Wbin) <<" Q2:" <<(Q2_bot(Q2bin)+Q2_top(Q2bin))/2.0 <<" var:" <<var <<" Xij:" <<Xij <<" Xij_val:" <<Xij_val <<" variable:" <<i <<"\n";
					idx = Histogram::Bin_Centering2_idx(W_center(Wbin),(Q2_bot(Q2bin)+Q2_top(Q2bin))/2.0,var,Xij,Xij_val,phi_val,cart[0],flags_);
					//fun::print_vector_idx(idx);
					if(Histogram::OK_Idx(idx)){
						dir_bc2_sub[Wbin][Q2bin+1][var+1][Xij+1][cart[0]+1]->cd();
						if(cart[0]==0){
							_Bin_Center2_hist1[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle(xlabel);
							_Bin_Center2_hist1[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Yield");
							_Bin_Center2_hist1[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
						}else if(cart[0]==1){
							_Bin_Center2_hist2[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle(xlabel);
							_Bin_Center2_hist2[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Yield");
							_Bin_Center2_hist2[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
						}else if(cart[0]==2){
							_Bin_Center2_hist3[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle(xlabel);
							_Bin_Center2_hist3[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Yield");
							_Bin_Center2_hist3[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
						}else if(cart[0]==3){
							_Bin_Center2_hist4[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle(xlabel);
							_Bin_Center2_hist4[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Yield");
							_Bin_Center2_hist4[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
						}
					}
				}
			}
		}
		std::cout<<"\tFinished Writing Bin Centering Plots for Polarization\n";
	}
	//*------------------------------End Bin Centering Corrections for Polarization------------------*

void Histogram::Make_Top_Check(std::shared_ptr<Flags> flags_){
	bool run = false;
	for(int i=0; i<4; i++){
		if(flags_->Flags::Plot_Top_Check(i)){
			if(fun::top_perform(_top_[i],flags_)){run = true;}
		}
	}
	if(!run){return;}
	std::cout<<"Making Top Check Histograms\n";
	char hname[100];
	sprintf(hname,"Topology_Cross_Check");
	_Top_Check_hist= new TH1D(hname,hname,_top_check_bins_,_top_check_min_,_top_check_max_);
	std::cout<<"\tMade " <<hname <<"\n";
	for(int i=0; i<4; i++){
		if(fun::top_perform(_top_[i],flags_)){
			sprintf(hname,"%s_Cross_Check",_top_[i]);
			_Top_Check_hist2[i]= new TH1D(hname,hname,_top_check_bins_,_top_check_min_,_top_check_max_);
			std::cout<<"\tMade " <<hname <<"\n";
		}
	}
}
void Histogram::Fill_Top_Check(bool mpro_pass_, bool mpip_pass_, bool mpim_pass_, bool mzero_pass_, int npro_, int npip_, int npim_, int nzero_, std::shared_ptr<Flags> flags_){
	bool run = false;
	for(int i=0; i<4; i++){
		if(flags_->Flags::Plot_Top_Check(i)){
			if(fun::top_perform(_top_[i],flags_)){run = true;}
		}
	}
	if(!run){return;}
	int input = 0;
	if(mpro_pass_ && fun::top_perform(_top_[0],flags_)){
		input += 1;
	}
	if(mpip_pass_ && fun::top_perform(_top_[1],flags_)){
		input += 2;
	}
	if(mpim_pass_ && fun::top_perform(_top_[2],flags_)){
		input += 4;
	}
	if(mzero_pass_ && fun::top_perform(_top_[3],flags_)){
		input += 8;
	}
	_Top_Check_hist->Fill(input);
	if(fun::top_perform(_top_[0],flags_)){
		_Top_Check_hist2[0]->Fill(npro_);
	}
	if(fun::top_perform(_top_[1],flags_)){
		_Top_Check_hist2[1]->Fill(npip_);
	}
	if(fun::top_perform(_top_[2],flags_)){
		_Top_Check_hist2[2]->Fill(npim_);
	}
	if(fun::top_perform(_top_[3],flags_)){
		_Top_Check_hist2[3]->Fill(nzero_);
	}
	
	
	
	
}

void Histogram::Write_Top_Check(std::shared_ptr<Flags> flags_){
	bool run = false;
	for(int i=0; i<4; i++){
		if(flags_->Flags::Plot_Top_Check(i)){
			if(fun::top_perform(_top_[i],flags_)){run = true;}
		}
	}
	if(!run){return;}
	std::cout<<"Writing Top Check Histograms\n";
	char dir_name[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_tc = _RootOutputFile->mkdir("Topology Check");
	//std::cout<<" Done\n";
	dir_tc->cd();
	_Top_Check_hist->SetXTitle("Topology Combination");
	_Top_Check_hist->SetYTitle("Num Events");
	_Top_Check_hist->Write();
	for(int i=0; i<4; i++){
		if(fun::top_perform(_top_[i],flags_)){
			_Top_Check_hist2[i]->SetXTitle("Num of simultaneous tops");
			_Top_Check_hist2[i]->SetYTitle("Num of Events");
			_Top_Check_hist2[i]->Write();
		}
	}
}

//********** Event Idx *********
/*int Histogram::Friend_Event_Idx(Event event_, int i_, int var_){
	switch(i_){
		case 0:
			return Histogram::Friend_W_idx(event_.Event::W());
		break;
		case 1:
			return Histogram::Friend_Q2_idx(event_.Event::Q2());
		break;
		case 2:
			return Histogram::Friend_MM_idx(event_.Event::MMb(var_),var_,event_.Event::W());
		break;
		case 3:
			return Histogram::Friend_MM2_idx(event_.Event::MM2b(var_),var_,event_.Event::W());
		break;
		case 4:
			return Histogram::Friend_theta_idx(event_.Event::Thetab(var_));
		break;
		case 5:
			return Histogram::Friend_alpha_idx(event_.Event::Alphab(var_));
		break;
		case 6:
			return Histogram::Friend_phi_idx(event_.Event::Phib(var_));
		break;
		default:
			std::cout<<"Incorrect Index of Indexes\n";
			return -1;
		break;
	}
}*/

int Histogram::Friend_Event_Idx1(float W_){
	return Histogram::Friend_W_idx(W_);
}
int Histogram::Friend_Event_Idx2(float Q2_){
	return Histogram::Friend_Q2_idx(Q2_);
}
int Histogram::Friend_Event_Idx3(float MM_, float W_,  int var_){
	return Histogram::Friend_MM_idx(MM_,var_,Histogram::Friend_W_idx(W_));
}
int Histogram::Friend_Event_Idx4(float MM2_, float W_, int var_){
	return Histogram::Friend_MM2_idx(MM2_,var_,Histogram::Friend_W_idx(W_));
}
int Histogram::Friend_Event_Idx5(float theta_){
	return Histogram::Friend_theta_idx(theta_);
}
int Histogram::Friend_Event_Idx6(float alpha_){
	return Histogram::Friend_alpha_idx(alpha_);
}
int Histogram::Friend_Event_Idx7(float phi_){
	return Histogram::Friend_phi_idx(phi_);
}


long Histogram::Number_of_Events(){
	return _nevnts;
}

long Histogram::Number_of_Attempted_Events(){
	return _n_attempted_evnts;
}

long Histogram::Number_of_Filled_Var(int i_){
	return _nvar[i_];
}

long Histogram::Number_of_Unfilled_Var(int i_){
	return _bvar[i_];
}

//*--------------------SC Segment Separated Kinematic Efficiency Plots ---------------------*
void Histogram::Kin_Eff_SC_Seg_Make(std::shared_ptr<Flags> flags_){
	bool plots = false;
	for(int i=0; i<4; i++){
		if(flags_->Flags::Plot_Kin_Eff(i)){
			plots = true;
		}
	}
	if(!plots){ return; }
	
	std::vector<long> space_dims(3);
	space_dims[2] = 4; //species
	space_dims[1] = 6; //Sector
	space_dims[0] = 48; //Segment
	char hname[200];
	CartesianGenerator cart(space_dims);
	TH2F_ptr_1d hist_1d;
	TH2F_ptr_2d hist_2d;

	int species_idx = -1;
	int sec_idx = -1;
	int seg_idx = -1;
	while(cart.GetNextCombination()){
		species_idx = cart[2];
		sec_idx = cart[1];
		seg_idx = cart[0];
		if(flags_->Flags::Plot_Kin_Eff(species_idx)){
			sprintf(hname,"Kinematic_Efficiency_%s_%s_Seg:%d",_species_[species_idx],_sector_[sec_idx],seg_idx+1);
			hist_1d.push_back(new TH2F(hname,hname,_cc_eff1_xbin_,_cc_eff1_xmin_,_cc_eff1_xmax_,_cc_eff1_ybin_,_cc_eff1_ymin_,_cc_eff1_ymax_));
		}
		if(cart[0] == space_dims[0]-1){
			if(hist_1d.size()>0){
				hist_2d.push_back(hist_1d);
				hist_1d.clear();
			}
			if(cart[1] == space_dims[1]-1){
				if(hist_2d.size()>0){
					_Kinematic_Eff_SC_Seg_hist.push_back(hist_2d);
					hist_2d.clear();
				}
			}
		}
	}
}

std::vector<int> Histogram::Kin_Eff_SC_Seg_idx(const char* species_, const char* sector_, int seg_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	if(!flags_->Flags::Plot_Kin_Eff(fun::species_idx(species_))){ 
		for(int i=0; i<3; i++){
			idx.push_back(-1);
		}
		return idx;
	}
	idx.push_back(fun::species_idx(species_)+fun::species_offset(species_,_kin_eff_cut_,flags_));
	idx.push_back(fun::sector_idx(sector_));
	idx.push_back(seg_);
	return idx;
}

void Histogram::Kin_Eff_SC_Seg_Fill(float p_, float theta_, float weight_, const char* species_, const char* sector_, int seg_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Kin_Eff(fun::species_idx(species_))){ return; }
	std::vector<int> idx = Histogram::Kin_Eff_SC_Seg_idx(species_,sector_,seg_,flags_);
	if(!Histogram::OK_Idx(idx)){ return; }
	_Kinematic_Eff_SC_Seg_hist[idx[0]][idx[1]][idx[2]]->Fill(p_,theta_,weight_);
}

void Histogram::Kin_Eff_SC_Seg_Write(std::shared_ptr<Flags> flags_){
	bool plots = false;
	for(int i=0; i<4; i++){
		if(flags_->Flags::Plot_Kin_Eff(i)){
			plots = true;
		}
	}
	if(!plots){ return; }
	std::cout<<"\tWriting Other Kinematic Efficiency Histograms\n";
	char dir_name[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_kin_eff2 = _RootOutputFile->mkdir("Kinematic Efficiency by Segment");
	dir_kin_eff2->cd();
	TDirectory* dir_kin_eff2_sub1[4];
	TDirectory* dir_kin_eff2_sub2[4][6];
	for(int j=0; j<4; j++){
		if(flags_->Flags::Plot_Kin_Eff(j)){
			sprintf(dir_name,"Kinematic Efficiency %s",_species_[j]);
			dir_kin_eff2_sub1[j] = dir_kin_eff2->mkdir(dir_name);
			for(int i=0; i<6; i++){
				sprintf(dir_name,"Kinematic Efficiency %s Sector:%d",_species_[j],i+1);
				dir_kin_eff2_sub2[j][i] = dir_kin_eff2_sub1[j]->mkdir(dir_name);
			}
		}
	}

	std::vector<long> space_dims(3);
	space_dims[2] = 4; //species
	space_dims[1] = 6; //Sector
	space_dims[0] = 48; //Segment
	CartesianGenerator cart(space_dims);

	int species_idx = -1;
	int sec_idx = -1;
	int seg_idx = -1;
	std::vector<int> idx;
	while(cart.GetNextCombination()){
		species_idx = cart[2];
		sec_idx = cart[1];
		seg_idx = cart[0];
		if(flags_->Flags::Plot_Kin_Eff(species_idx)){
			idx = Histogram::Kin_Eff_SC_Seg_idx(_species_[species_idx],_sector_[sec_idx],seg_idx,flags_);
			if(Histogram::OK_Idx(idx)){
				dir_kin_eff2_sub2[species_idx][sec_idx]->cd();
				_Kinematic_Eff_SC_Seg_hist[idx[0]][idx[1]][idx[2]]->SetXTitle("Momentum (GeV)");
				_Kinematic_Eff_SC_Seg_hist[idx[0]][idx[1]][idx[2]]->SetXTitle("Theta (degrees)");
				_Kinematic_Eff_SC_Seg_hist[idx[0]][idx[1]][idx[2]]->Write();
			}
		}
	}
}


//*----------- Sim Vertex Determination ------ 
void Histogram::Sim_Vertex_Make(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Sim()){ return;}
	char hname[100];
	sprintf(hname,"Sim_Vertex_Runs");
	_Sim_Vertex_Runs = new TH1D(hname,hname,600,106000,36000000);
}
void Histogram::Sim_Vertex_Fill(int run_num_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Sim()){ return;}
	_Sim_Vertex_Runs->Fill(run_num_);
}
void Histogram::Sim_Vertex_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Sim()){ return;}
	std::cout<<"\tWriting Sim Vertex Runs\n";
	char dir_name[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_sim_vert= _RootOutputFile->mkdir("Sim Vertex Runs");
	dir_sim_vert->cd();
	_Sim_Vertex_Runs->SetXTitle("Sim Run Number");
	_Sim_Vertex_Runs->SetYTitle("N Events passed Vertex Cut");
	_Sim_Vertex_Runs->Write();
}


//*---------- Helicity Ratio --------
void Histogram::Hel_Ratio_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Sim()){ return;}
	if(!flags_->Flags::Helicity()){ return;}
	std::cout<<"Making Helicity Ratio Histograms\n";
	Double_t Q2_bin[6] = {_Q2_bins_[0],_Q2_bins_[1],_Q2_bins_[2],_Q2_bins_[3],_Q2_bins_[4],_Q2_bins_[5]};
	char hname[100];
	sprintf(hname,"Positive Helicity Raw Yield");
	_Pos_Hel_wq2_hist[0] = new TH2D(hname,hname,Histogram::W_bins(),_W_min_,_W_max_,5,Q2_bin);
	sprintf(hname,"Positive Helicity Event Yield");
	_Pos_Hel_wq2_hist[1] = new TH2D(hname,hname,Histogram::W_bins(),_W_min_,_W_max_,5,Q2_bin);
	sprintf(hname,"Negative Helicity Raw Yield");
	_Neg_Hel_wq2_hist[0] = new TH2D(hname,hname,Histogram::W_bins(),_W_min_,_W_max_,5,Q2_bin);
	sprintf(hname,"Negative Helicity Event Yield");
	_Neg_Hel_wq2_hist[1] = new TH2D(hname,hname,Histogram::W_bins(),_W_min_,_W_max_,5,Q2_bin);
	sprintf(hname,"Pos/Neg Helicity Ratio Raw");
	_Ratio_Hel_wq2_hist[0] = new TH2D(hname,hname,Histogram::W_bins(),_W_min_,_W_max_,5,Q2_bin);
	sprintf(hname,"Pos/Neg Helicity Ratio Event");
	_Ratio_Hel_wq2_hist[1] = new TH2D(hname,hname,Histogram::W_bins(),_W_min_,_W_max_,5,Q2_bin);

}
void Histogram::Hel_Ratio_Fill(float weight_, int hel_idx_, float W_, float Q2_, int event_stat_,std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Sim()){ return;}
	if(!flags_->Flags::Helicity()){ return;}
	//if(event_stat_==1){
	//	std::cout<<"Filling Hel Ratio " <<weight_ <<" hel_idx:" <<hel_idx_ <<" W:" <<W_ <<" Q2:" <<Q2_ <<" event_stat:" <<event_stat_ <<"\n";
	//}
	if(hel_idx_==0){
		_Pos_Hel_wq2_hist[event_stat_]->Fill(W_,Q2_,weight_);
	}else if(hel_idx_==1){
		_Neg_Hel_wq2_hist[event_stat_]->Fill(W_,Q2_,weight_);
	}
}
void Histogram::Hel_Ratio_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Sim()){ return;}
	if(!flags_->Flags::Helicity()){ return;}
	std::cout<<"\tWriting Helicity Yields\n";
	char dir_name[100];
	char hname[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_hel_ratio= _RootOutputFile->mkdir("Helicity Yields and Ratios");
	dir_hel_ratio->cd();
	char* Rname[2] = {"Raw","Event"};
	for(int i=0; i< 2; i++){
		_Pos_Hel_wq2_hist[i]->SetXTitle("W (GeV)");
		_Pos_Hel_wq2_hist[i]->SetYTitle("Q2 (GeV^2)");
		_Pos_Hel_wq2_hist[i]->Write();
		_Neg_Hel_wq2_hist[i]->SetXTitle("W (GeV)");
		_Neg_Hel_wq2_hist[i]->SetYTitle("Q2 (GeV^2)");
		_Neg_Hel_wq2_hist[i]->Write();
		_Ratio_Hel_wq2_hist[i] = (TH2D*) _Pos_Hel_wq2_hist[i]->Clone();
		_Ratio_Hel_wq2_hist[i]->Divide(_Neg_Hel_wq2_hist[i]);
		_Ratio_Hel_wq2_hist[i]->SetXTitle("W (GeV)");
		_Ratio_Hel_wq2_hist[i]->SetYTitle("Q2 (GeV^2)");
		sprintf(hname,"Pos/Neg Helicity Ratio %s",Rname[i]);
		_Ratio_Hel_wq2_hist[i]->SetNameTitle(hname,hname);
		_Ratio_Hel_wq2_hist[i]->Write();
	}
}
	
void Histogram::Golden_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Sim()){ return;}
	if(!flags_->Flags::Golden_Run()){ return ;}
	char hname[100];
	sprintf(hname,"Normalized Integrated Faraday Cup Charge by Run E1-6");
	_Golden_Run1 = new TH1D(hname,hname,945,30539.5,31484.5);

	sprintf(hname,"Normalized Integrated Faraday Cup Charge Distribution E1-6");
	_Golden_Run2 = new TH1D(hname,hname,90,0.0,10.0);

	sprintf(hname,"Integrated Faraday Cup Charge by Run E1-6");
	_Golden_Run_Charge = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"Integrated Faraday Cup Charge by Run E1-6 Pos");
	_Golden_Run_Charge_pos = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"Integrated Faraday Cup Charge by Run E1-6 Neg");
	_Golden_Run_Charge_neg = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"Integrated Faraday Cup Charge by Run E1-6 None");
	_Golden_Run_Charge_none = new TH1D(hname,hname,945,30539.5,31484.5);

	sprintf(hname,"Integrated Faraday Cup Events by Run E1-6");
	_Golden_Run_Events = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"Integrated Faraday Cup Events by Run E1-6 Pos");
	_Golden_Run_Events_pos = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"Integrated Faraday Cup Events by Run E1-6 Neg");
	_Golden_Run_Events_neg = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"Integrated Faraday Cup Events by Run E1-6 None");
	_Golden_Run_Events_none = new TH1D(hname,hname,945,30539.5,31484.5);

	

	sprintf(hname,"2 Integrated Faraday Cup Charge by Run E1-6 ");
	_Golden_Run_Charge2 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"2 Integrated Faraday Cup Charge by Run E1-6 Pos");
	_Golden_Run_Charge_pos2 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"2 Integrated Faraday Cup Charge by Run E1-6 Neg");
	_Golden_Run_Charge_neg2 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"2 Integrated Faraday Cup Charge by Run E1-6 None");
	_Golden_Run_Charge_none2 = new TH1D(hname,hname,945,30539.5,31484.5);

	sprintf(hname,"2 Integrated Faraday Cup Events by Run E1-6");
	_Golden_Run_Events2 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"2 Integrated Faraday Cup Events by Run E1-6 Pos");
	_Golden_Run_Events_pos2 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"2 Integrated Faraday Cup Events by Run E1-6 Neg");
	_Golden_Run_Events_neg2 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"2 Integrated Faraday Cup Events by Run E1-6 None");
	_Golden_Run_Events_none2 = new TH1D(hname,hname,945,30539.5,31484.5);

	sprintf(hname,"3 Integrated Faraday Cup Charge by Run E1-6 ");
	_Golden_Run_Charge3 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"3 Integrated Faraday Cup Charge by Run E1-6 Pos");
	_Golden_Run_Charge_pos3 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"3 Integrated Faraday Cup Charge by Run E1-6 Neg");
	_Golden_Run_Charge_neg3 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"3 Integrated Faraday Cup Charge by Run E1-6 None");
	_Golden_Run_Charge_none3 = new TH1D(hname,hname,945,30539.5,31484.5);

	sprintf(hname,"3 Integrated Faraday Cup Events by Run E1-6");
	_Golden_Run_Events3 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"3 Integrated Faraday Cup Events by Run E1-6 Pos");
	_Golden_Run_Events_pos3 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"3 Integrated Faraday Cup Events by Run E1-6 Neg");
	_Golden_Run_Events_neg3 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"3 Integrated Faraday Cup Events by Run E1-6 None");
	_Golden_Run_Events_none3 = new TH1D(hname,hname,945,30539.5,31484.5);


	sprintf(hname,"4 Integrated Faraday Cup Charge by Run E1-6 ");
	_Golden_Run_Charge4 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"4 Integrated Faraday Cup Charge by Run E1-6 Pos");
	_Golden_Run_Charge_pos4 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"4 Integrated Faraday Cup Charge by Run E1-6 Neg");
	_Golden_Run_Charge_neg4 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"4 Integrated Faraday Cup Charge by Run E1-6 None");
	_Golden_Run_Charge_none4 = new TH1D(hname,hname,945,30539.5,31484.5);

	sprintf(hname,"4 Integrated Faraday Cup Events by Run E1-6");
	_Golden_Run_Events4 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"4 Integrated Faraday Cup Events by Run E1-6 Pos");
	_Golden_Run_Events_pos4 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"4 Integrated Faraday Cup Events by Run E1-6 Neg");
	_Golden_Run_Events_neg4 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"4 Integrated Faraday Cup Events by Run E1-6 None");
	_Golden_Run_Events_none4 = new TH1D(hname,hname,945,30539.5,31484.5);

	sprintf(hname,"5 Integrated Faraday Cup Charge by Run E1-6 ");
	_Golden_Run_Charge5 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"5 Integrated Faraday Cup Charge by Run E1-6 Pos");
	_Golden_Run_Charge_pos5 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"5 Integrated Faraday Cup Charge by Run E1-6 Neg");
	_Golden_Run_Charge_neg5 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"5 Integrated Faraday Cup Charge by Run E1-6 None");
	_Golden_Run_Charge_none5 = new TH1D(hname,hname,945,30539.5,31484.5);

	sprintf(hname,"5 Integrated Faraday Cup Events by Run E1-6");
	_Golden_Run_Events5 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"5 Integrated Faraday Cup Events by Run E1-6 Pos");
	_Golden_Run_Events_pos5 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"5 Integrated Faraday Cup Events by Run E1-6 Neg");
	_Golden_Run_Events_neg5 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"5 Integrated Faraday Cup Events by Run E1-6 None");
	_Golden_Run_Events_none5 = new TH1D(hname,hname,945,30539.5,31484.5);


	sprintf(hname,"6 Integrated Faraday Cup Charge by Run E1-6 ");
	_Golden_Run_Charge6 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"6 Integrated Faraday Cup Charge by Run E1-6 Pos");
	_Golden_Run_Charge_pos6 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"6 Integrated Faraday Cup Charge by Run E1-6 Neg");
	_Golden_Run_Charge_neg6 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"6 Integrated Faraday Cup Charge by Run E1-6 None");
	_Golden_Run_Charge_none6 = new TH1D(hname,hname,945,30539.5,31484.5);

	sprintf(hname,"6 Integrated Faraday Cup Events by Run E1-6");
	_Golden_Run_Events6 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"6 Integrated Faraday Cup Events by Run E1-6 Pos");
	_Golden_Run_Events_pos6 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"6 Integrated Faraday Cup Events by Run E1-6 Neg");
	_Golden_Run_Events_neg6 = new TH1D(hname,hname,945,30539.5,31484.5);
	sprintf(hname,"6 Integrated Faraday Cup Events by Run E1-6 None");
	_Golden_Run_Events_none6 = new TH1D(hname,hname,945,30539.5,31484.5);



}
void Histogram::Golden_Fill(int run_num_, float q_sum_, int run_evts_){
	_Golden_Run1->Fill(run_num_,q_sum_/(float)run_evts_);
	_Golden_Run2->Fill(q_sum_/(float)run_evts_);

	if(q_sum_/(float)run_evts_ > 10.0){
		std::cout<<"Amount Exceeded for run " <<run_num_ <<" with " <<q_sum_/(float)run_evts_ <<"\n";
	}
}

void Histogram::Golden_Indiv_Fill(int run_num_, float delta_q_, float hel_){
	if(hel_ != 0.0){
		if(hel_ == 1.0){
			_Golden_Run_Charge_pos->Fill(run_num_,delta_q_);
			_Golden_Run_Events_pos->Fill(run_num_);
		}else if(hel_ == -1.0){
			_Golden_Run_Charge_neg->Fill(run_num_,delta_q_);
			_Golden_Run_Events_neg->Fill(run_num_);
		}else{
			_Golden_Run_Charge_none->Fill(run_num_,delta_q_);
			_Golden_Run_Events_none->Fill(run_num_);
		}
	}else{
		_Golden_Run_Charge_none->Fill(run_num_,delta_q_);
		_Golden_Run_Events_none->Fill(run_num_);
	}
	_Golden_Run_Charge->Fill(run_num_,delta_q_);
	_Golden_Run_Events->Fill(run_num_);
}
void Histogram::Golden_Indiv_Fill2(int run_num_, float delta_q_, float hel_){
	if(hel_ != 0.0){
		if(hel_ == 1.0){
			_Golden_Run_Charge_pos2->Fill(run_num_,delta_q_);
			_Golden_Run_Events_pos2->Fill(run_num_);
		}else if(hel_ == -1.0){
			_Golden_Run_Charge_neg2->Fill(run_num_,delta_q_);
			_Golden_Run_Events_neg2->Fill(run_num_);
		}else{
			_Golden_Run_Charge_none2->Fill(run_num_,delta_q_);
			_Golden_Run_Events_none2->Fill(run_num_);
		}
	}else{
		_Golden_Run_Charge_none2->Fill(run_num_,delta_q_);
		_Golden_Run_Events_none2->Fill(run_num_);
	}
	_Golden_Run_Charge2->Fill(run_num_,delta_q_);
	_Golden_Run_Events2->Fill(run_num_);
}

void Histogram::Golden_Indiv_Fill3(int run_num_, float delta_q_, float hel_){
	if(hel_ != 0.0){
		if(hel_ == 1.0){
			_Golden_Run_Charge_pos3->Fill(run_num_,delta_q_);
			_Golden_Run_Events_pos3->Fill(run_num_);
		}else if(hel_ == -1.0){
			_Golden_Run_Charge_neg3->Fill(run_num_,delta_q_);
			_Golden_Run_Events_neg3->Fill(run_num_);
		}else{
			_Golden_Run_Charge_none3->Fill(run_num_,delta_q_);
			_Golden_Run_Events_none3->Fill(run_num_);
		}
	}else{
		_Golden_Run_Charge_none3->Fill(run_num_,delta_q_);
		_Golden_Run_Events_none3->Fill(run_num_);
	}
	_Golden_Run_Charge3->Fill(run_num_,delta_q_);
	_Golden_Run_Events3->Fill(run_num_);
}

void Histogram::Golden_Indiv_Fill4(int run_num_, float delta_q_, float hel_){
	if(hel_ != 0.0){
		if(hel_ == 1.0){
			_Golden_Run_Charge_pos4->Fill(run_num_,delta_q_);
			_Golden_Run_Events_pos4->Fill(run_num_);
		}else if(hel_ == -1.0){
			_Golden_Run_Charge_neg4->Fill(run_num_,delta_q_);
			_Golden_Run_Events_neg4->Fill(run_num_);
		}else{
			_Golden_Run_Charge_none4->Fill(run_num_,delta_q_);
			_Golden_Run_Events_none4->Fill(run_num_);
		}
	}else{
		_Golden_Run_Charge_none4->Fill(run_num_,delta_q_);
		_Golden_Run_Events_none4->Fill(run_num_);
	}
	_Golden_Run_Charge4->Fill(run_num_,delta_q_);
	_Golden_Run_Events4->Fill(run_num_);
}
void Histogram::Golden_Indiv_Fill5(int run_num_, float delta_q_, float hel_){
	if(hel_ != 0.0){
		if(hel_ == 1.0){
			_Golden_Run_Charge_pos5->Fill(run_num_,delta_q_);
			_Golden_Run_Events_pos5->Fill(run_num_);
		}else if(hel_ == -1.0){
			_Golden_Run_Charge_neg5->Fill(run_num_,delta_q_);
			_Golden_Run_Events_neg5->Fill(run_num_);
		}else{
			_Golden_Run_Charge_none5->Fill(run_num_,delta_q_);
			_Golden_Run_Events_none5->Fill(run_num_);
		}
	}else{
		_Golden_Run_Charge_none5->Fill(run_num_,delta_q_);
		_Golden_Run_Events_none5->Fill(run_num_);
	}
	_Golden_Run_Charge5->Fill(run_num_,delta_q_);
	_Golden_Run_Events5->Fill(run_num_);
}
void Histogram::Golden_Indiv_Fill6(int run_num_, float delta_q_, float hel_){
	if(hel_ != 0.0){
		if(hel_ == 1.0){
			_Golden_Run_Charge_pos6->Fill(run_num_,delta_q_);
			_Golden_Run_Events_pos6->Fill(run_num_);
		}else if(hel_ == -1.0){
			_Golden_Run_Charge_neg6->Fill(run_num_,delta_q_);
			_Golden_Run_Events_neg6->Fill(run_num_);
		}else{
			_Golden_Run_Charge_none6->Fill(run_num_,delta_q_);
			_Golden_Run_Events_none6->Fill(run_num_);
		}
	}else{
		_Golden_Run_Charge_none6->Fill(run_num_,delta_q_);
		_Golden_Run_Events_none6->Fill(run_num_);
	}
	_Golden_Run_Charge6->Fill(run_num_,delta_q_);
	_Golden_Run_Events6->Fill(run_num_);
}

void Histogram::Golden_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Sim()){ return;}
	if(!flags_->Flags::Golden_Run()){ return ;}
	char dir_name[100];
	char hname[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_gold= _RootOutputFile->mkdir("Golden Run Determination");
	dir_gold->cd();
	_Golden_Run1->SetXTitle("Run Number");
	_Golden_Run1->SetYTitle("Normalized Faraday Cup Charge");
	_Golden_Run1->Write();

	_Golden_Run2->SetXTitle("Normalized Faraday Cup Charge");
	_Golden_Run2->SetYTitle("Number of Runs");
	_Golden_Run2->Write();

	_Golden_Run_Charge->SetXTitle("Run Number");
	_Golden_Run_Charge->SetYTitle("Faraday Cup Charge");
	_Golden_Run_Charge->Write();

	_Golden_Run_Events->SetXTitle("Run Number");
	_Golden_Run_Events->SetYTitle("Events");
	_Golden_Run_Events->Write();

	_Golden_Run_Charge2->SetXTitle("Run Number");
	_Golden_Run_Charge2->SetYTitle("Faraday Cup Charge");
	_Golden_Run_Charge2->Write();

	_Golden_Run_Events2->SetXTitle("Run Number");
	_Golden_Run_Events2->SetYTitle("Events");
	_Golden_Run_Events2->Write();

	_Golden_Run_Charge3->SetXTitle("Run Number");
	_Golden_Run_Charge3->SetYTitle("Faraday Cup Charge");
	_Golden_Run_Charge3->Write();

	_Golden_Run_Events3->SetXTitle("Run Number");
	_Golden_Run_Events3->SetYTitle("Events");
	_Golden_Run_Events3->Write();

	_Golden_Run_Charge4->SetXTitle("Run Number");
	_Golden_Run_Charge4->SetYTitle("Faraday Cup Charge");
	_Golden_Run_Charge4->Write();

	_Golden_Run_Events4->SetXTitle("Run Number");
	_Golden_Run_Events4->SetYTitle("Events");
	_Golden_Run_Events4->Write();

	_Golden_Run_Charge5->SetXTitle("Run Number");
	_Golden_Run_Charge5->SetYTitle("Faraday Cup Charge");
	_Golden_Run_Charge5->Write();

	_Golden_Run_Events5->SetXTitle("Run Number");
	_Golden_Run_Events5->SetYTitle("Events");
	_Golden_Run_Events5->Write();

	_Golden_Run_Charge6->SetXTitle("Run Number");
	_Golden_Run_Charge6->SetYTitle("Faraday Cup Charge");
	_Golden_Run_Charge6->Write();

	_Golden_Run_Events6->SetXTitle("Run Number");
	_Golden_Run_Events6->SetYTitle("Events");
	_Golden_Run_Events6->Write();

	if(flags_->Flags::Helicity()){
		_Golden_Run_Charge_pos->SetXTitle("Run Number");
		_Golden_Run_Charge_pos->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_pos->Write();

		_Golden_Run_Events_pos->SetXTitle("Run Number");
		_Golden_Run_Events_pos->SetYTitle("Events");
		_Golden_Run_Events_pos->Write();

		_Golden_Run_Charge_neg->SetXTitle("Run Number");
		_Golden_Run_Charge_neg->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_neg->Write();

		_Golden_Run_Events_neg->SetXTitle("Run Number");
		_Golden_Run_Events_neg->SetYTitle("Events");
		_Golden_Run_Events_neg->Write();

		_Golden_Run_Charge_none->SetXTitle("Run Number");
		_Golden_Run_Charge_none->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_none->Write();

		_Golden_Run_Events_none->SetXTitle("Run Number");
		_Golden_Run_Events_none->SetYTitle("Events");
		_Golden_Run_Events_none->Write();


		_Golden_Run_Charge_pos2->SetXTitle("Run Number");
		_Golden_Run_Charge_pos2->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_pos2->Write();

		_Golden_Run_Events_pos2->SetXTitle("Run Number");
		_Golden_Run_Events_pos2->SetYTitle("Events");
		_Golden_Run_Events_pos2->Write();

		_Golden_Run_Charge_neg2->SetXTitle("Run Number");
		_Golden_Run_Charge_neg2->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_neg2->Write();

		_Golden_Run_Events_neg2->SetXTitle("Run Number");
		_Golden_Run_Events_neg2->SetYTitle("Events");
		_Golden_Run_Events_neg2->Write();

		_Golden_Run_Charge_none2->SetXTitle("Run Number");
		_Golden_Run_Charge_none2->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_none2->Write();

		_Golden_Run_Events_none2->SetXTitle("Run Number");
		_Golden_Run_Events_none2->SetYTitle("Events");
		_Golden_Run_Events_none2->Write();

		
		_Golden_Run_Charge_pos3->SetXTitle("Run Number");
		_Golden_Run_Charge_pos3->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_pos3->Write();

		_Golden_Run_Events_pos3->SetXTitle("Run Number");
		_Golden_Run_Events_pos3->SetYTitle("Events");
		_Golden_Run_Events_pos3->Write();

		_Golden_Run_Charge_neg3->SetXTitle("Run Number");
		_Golden_Run_Charge_neg3->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_neg3->Write();

		_Golden_Run_Events_neg3->SetXTitle("Run Number");
		_Golden_Run_Events_neg3->SetYTitle("Events");
		_Golden_Run_Events_neg3->Write();

		_Golden_Run_Charge_none3->SetXTitle("Run Number");
		_Golden_Run_Charge_none3->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_none3->Write();

		_Golden_Run_Events_none3->SetXTitle("Run Number");
		_Golden_Run_Events_none3->SetYTitle("Events");
		_Golden_Run_Events_none3->Write();

		_Golden_Run_Charge_pos4->SetXTitle("Run Number");
		_Golden_Run_Charge_pos4->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_pos4->Write();

		_Golden_Run_Events_pos4->SetXTitle("Run Number");
		_Golden_Run_Events_pos4->SetYTitle("Events");
		_Golden_Run_Events_pos4->Write();

		_Golden_Run_Charge_neg4->SetXTitle("Run Number");
		_Golden_Run_Charge_neg4->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_neg4->Write();

		_Golden_Run_Events_neg4->SetXTitle("Run Number");
		_Golden_Run_Events_neg4->SetYTitle("Events");
		_Golden_Run_Events_neg4->Write();

		_Golden_Run_Charge_none4->SetXTitle("Run Number");
		_Golden_Run_Charge_none4->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_none4->Write();

		_Golden_Run_Events_none4->SetXTitle("Run Number");
		_Golden_Run_Events_none4->SetYTitle("Events");
		_Golden_Run_Events_none4->Write();

		_Golden_Run_Charge_pos5->SetXTitle("Run Number");
		_Golden_Run_Charge_pos5->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_pos5->Write();

		_Golden_Run_Events_pos5->SetXTitle("Run Number");
		_Golden_Run_Events_pos5->SetYTitle("Events");
		_Golden_Run_Events_pos5->Write();

		_Golden_Run_Charge_neg5->SetXTitle("Run Number");
		_Golden_Run_Charge_neg5->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_neg5->Write();

		_Golden_Run_Events_neg5->SetXTitle("Run Number");
		_Golden_Run_Events_neg5->SetYTitle("Events");
		_Golden_Run_Events_neg5->Write();

		_Golden_Run_Charge_none5->SetXTitle("Run Number");
		_Golden_Run_Charge_none5->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_none5->Write();

		_Golden_Run_Events_none5->SetXTitle("Run Number");
		_Golden_Run_Events_none5->SetYTitle("Events");
		_Golden_Run_Events_none5->Write();

		_Golden_Run_Charge_pos6->SetXTitle("Run Number");
		_Golden_Run_Charge_pos6->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_pos6->Write();

		_Golden_Run_Events_pos6->SetXTitle("Run Number");
		_Golden_Run_Events_pos6->SetYTitle("Events");
		_Golden_Run_Events_pos6->Write();

		_Golden_Run_Charge_neg6->SetXTitle("Run Number");
		_Golden_Run_Charge_neg6->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_neg6->Write();

		_Golden_Run_Events_neg6->SetXTitle("Run Number");
		_Golden_Run_Events_neg6->SetYTitle("Events");
		_Golden_Run_Events_neg6->Write();

		_Golden_Run_Charge_none6->SetXTitle("Run Number");
		_Golden_Run_Charge_none6->SetYTitle("Faraday Cup Charge");
		_Golden_Run_Charge_none6->Write();

		_Golden_Run_Events_none6->SetXTitle("Run Number");
		_Golden_Run_Events_none6->SetYTitle("Events");
		_Golden_Run_Events_none6->Write();

	}
}