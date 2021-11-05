#include "histogram.hpp"


Histogram::Histogram(std::shared_ptr<Flags> flags_){
	//_RootOutputFile = fun::Name_File(output_file);
	//def = new TCanvas("def");
	Histogram::WQ2_Make(flags_);
	Histogram::Fid_Make(flags_);
	Histogram::SF_Make(flags_);
	Histogram::Delta_Make(flags_);
	Histogram::CC_Make(flags_);
	Histogram::Vertex_Make(flags_);
	Histogram::MM_Make(flags_);
	Histogram::CC_Eff_Make(flags_);
	Histogram::Friend_Make(flags_);
	//Histogram::Cross_Make(flags_);
	//Histogram::XY_Make(flags_);
	//Histogram::Fid_Det_Make(flags_);
	//Histogram::Charge_Make(flags_);
	//Histogram::Thrown_Make(flags_);
	//Histogram::WQ2_sf_Make(_envi);
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
	_RootOutputFile = fun::Name_File(flags_);
	def = new TCanvas("def");
	std::cout<< "Writing Plots" <<std::endl;
	_RootOutputFile->cd();
	std::cout<<"Writing WQ2 Histograms\n";
	Histogram::WQ2_Write(flags_);
	std::cout<<"Writing Fid Histograms\n";
	Histogram::Fid_Write(flags_);
	std::cout<<"Writing SF Histograms\n";
	Histogram::SF_Write(flags_);
	Histogram::Delta_Write(flags_);
	std::cout<<"Writing CC Histograms\n";
	Histogram::CC_Write(flags_);
	std::cout<<"Writing Vertex Histograms\n";
	Histogram::Vertex_Write(flags_);
	std::cout<<"Writing MM Histograms\n";
	Histogram::MM_Write(flags_);
	Histogram::CC_Eff_Write(flags_);
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
		/*
		for(int i=0; i<5; i++){
			std::cout<<space_dims[i] <<" ";
		}
		std::cout<<"\n";*/
		char hname[100];//For naming histograms
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
							//std::cout<<"Filled WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
							sprintf(hname,"W_Q2_%s_%s",_recon_[cart[0]],_weight_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
							_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
						}
					}
				}else{//Reconstructed Stuff
					if(_ecuts_[cart[4]]== _event_ && _top_[cart[2]]!=_mnone_ && _cut_[cart[3]]!=_no_cut_){//Event Selection
						//std::cout<<"\t\tMaking Event Selection WQ2\n";
						//std::cout<<"Filled WQ2 idx: "  <<_WQ2_hist.size() <<" "<<plot_4d.size() << " " <<plot_3d.size() << " " <<plot_2d.size() << " " <<plot_1d.size() <<"\n";
						sprintf(hname,"W_Q2_%s_%s_%s_%s",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]]);
						plot_1d.push_back(new TH2F(hname,hname,_wq2_xbin_,_wq2_xmin_,_wq2_xmax_,_wq2_ybin_,_wq2_ymin_,_wq2_ymax_));
						_WQ2_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]]=true;
						//fun::print_vector_idx(cart);
						//fun::print_vector_idx(WQ2_cut_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_weight_[cart[1]],_recon_[cart[0]],flags_));
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
			if(cart[0] == space_dims[0]-1 && plot_1d.size()>0){
				//std::cout<<"\t\tPush back 1d\n";
				plot_2d.push_back(plot_1d);
				//std::cout<<"Size of 1d: " <<plot_1d.size() <<" | Size of 2d: " <<plot_2d.size()<<"\n";
				plot_1d.clear();
			}
			if(cart[1] == space_dims[1]-1 && plot_2d.size()>0 && cart[0] == space_dims[0]-1){
				//std::cout<<"\t\tPush back 2d\n";
				plot_3d.push_back(plot_2d);
				//std::cout<<"Size of 2d: " <<plot_2d.size() <<" | Size of 3d: " <<plot_3d.size()<<"\n";
				plot_2d.clear();
			}
			if(cart[2] == space_dims[2]-1 && plot_3d.size()>0 && cart[1] == space_dims[1]-1 && cart[0] == space_dims[0]-1){
				//std::cout<<"\t\tPush back 3d\n";
				plot_4d.push_back(plot_3d);
				//std::cout<<"Size of 3d: " <<plot_3d.size() <<" | Size of 4d: " <<plot_4d.size()<<"\n";
				plot_3d.clear();
			}
			if(cart[3] == space_dims[3]-1 && plot_4d.size()>0 && cart[1] == space_dims[1]-1 && cart[0] == space_dims[0]-1 && cart[2]==space_dims[2]-1){
				//std::cout<<"\t\tPush back 4d\n";
				_WQ2_hist.push_back(plot_4d);
				//std::cout<<"Size of 4d: " <<plot_4d.size() <<" | Size of 5d: " <<_WQ2_hist.size()<<"\n";
				plot_4d.clear();
			}
			if(cart[4] == space_dims[4]-1 && cart[3] == space_dims[3]-1 && cart[1] == space_dims[1]-1 && cart[0] == space_dims[0]-1 && cart[2]==space_dims[2]-1){
				//std::cout<<"Size of 5d: " <<_WQ2_hist.size() <<"\n";
			}
		}
	}
}

std::vector<int> Histogram::WQ2_cut_idx(const char* ecut_, const char * cut_, const char* top_, const char * weight_, const char * recon_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	int ecut_idx = -1; 
	int cut_idx = -1; 
	//std::cout<<"input to WQ2 Cut idx: " <<ecut_ <<"\n\t" <<cut_ <<"\n\t" <<top_ <<"\n\t" <<weight_ <<"\n\t" <<recon_  <<"\n";
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
	if(flags_->Plot_WQ2()){
		//std::cout<<"Trying to fill WQ2: " <<ecut_ <<" " <<cut_ <<" " <<top_ <<" " <<thrown_ <<"\n";
		std::vector<int> idx;
		if(flags_->Flags::Sim()){//Weighted
			idx = Histogram::WQ2_cut_idx(ecut_,cut_,top_,_weighted_,thrown_,flags_);
			if(Histogram::OK_Idx(idx)){
				if(Made_WQ2_idx(ecut_,cut_,top_,_nweighted_,thrown_)){	
				//if(Histogram::Good_WQ2_idx(ecut_,cut_,top_,_nweighted_,_recon_[fun::truth_idx(thrown_)],flags_)){
					_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(W_,Q2_,weight_);
				}
			}
			idx.clear();
		}
		//if(Histogram::OK_Idx(idx)){
		if(Made_WQ2_idx(ecut_,cut_,top_,_nweighted_,thrown_)){
			idx = Histogram::WQ2_cut_idx(ecut_,cut_,top_,_nweighted_,thrown_,flags_);
			if(Histogram::OK_Idx(idx)){
				_WQ2_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(W_,Q2_,1.0);
			}
			idx.clear();
		}
		
		//std::cout<<"input to WQ2 Filling: " <<ecut_ <<"\n\t" <<cut_ <<"\n\t" <<top_ <<"\n\t" <<weight_ <<"\n\t" <<_recon_[fun::truth_idx(thrown_)]  <<"\n";
		//std::cout<<"Trying to fill WQ2:" <<ecut_ <<" " <<cut_ <<" " <<top_ <<" " <<_recon_[fun::truth_idx(thrown_)] <<"\n";
		
	}
}

void Histogram::WQ2_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Plot_WQ2()){
		std::cout<<"\tinside Writing WQ2\n";
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
	if(flags_->Plot_Fid(0) || flags_->Plot_Fid(1) || flags_->Plot_Fid(2) || flags_->Plot_Fid(3)){
		std::cout<<"Making Fid Histograms\n";
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
		std::vector<long> space_dims(6);
		space_dims[5] = std::distance(std::begin(_species_), std::end(_species_));//{ele,pro,pip,pim}
		space_dims[4] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));//Electron Cuts (electrons have more cuts than hadrons)
		space_dims[3]	= std::distance(std::begin(_sector_), std::end(_sector_));//Sectors
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));//Cut and anti cut
		space_dims[1] = std::distance(std::begin(_top_), std::end(_top_));//{Topologies  + all combined + None}
		space_dims[0] = std::distance(std::begin(_weight_), std::end(_weight_));//Not Weight vs. Weighted
		CartesianGenerator cart(space_dims);
		CartesianGenerator cart2(space_dims);
		char hname[100];//For naming histograms
		int idx = 0; 

		while(cart2.GetNextCombination()){
			//if(_cut_[cart[2]] != _no_cut_){
				made_1d.push_back(false);
			//}
			if(cart2[0]==space_dims[0]-1 && made_1d.size()>0){
				made_2d.push_back(made_1d);
				//std::cout<<"\t\t\t\t\tLength of made_1d: " <<made_1d.size() <<"\n";
				made_1d.clear();
			}
			if(cart2[1]==space_dims[1]-1 && cart2[0]==space_dims[0]-1 && made_2d.size()>0){
				made_3d.push_back(made_2d);
				//std::cout<<"\t\t\t\tLength of made_2d: " <<made_2d.size() <<"\n";
				made_2d.clear();
			}
			if(cart2[2]==space_dims[2]-1 && cart2[1]==space_dims[1]-1 && cart2[0]==space_dims[0]-1 && made_3d.size()>0){
				made_4d.push_back(made_3d);
				//std::cout<<"\t\t\tLength of made_3d: " <<made_3d.size() <<"\n";
				made_3d.clear();
			}
			if(cart2[3]==space_dims[3]-1 && cart2[2]==space_dims[2]-1 && cart2[1]==space_dims[1]-1 && cart2[0]==space_dims[0]-1 && made_4d.size()>0){
				made_5d.push_back(made_4d);
				//std::cout<<"\t\tLength of made_4d: " <<made_4d.size() <<"\n";
				made_4d.clear();
			}
			if(cart2[4]==space_dims[4]-1 && cart2[3]==space_dims[3]-1 && cart2[2]==space_dims[2]-1 && cart2[1]==space_dims[1]-1 && cart2[0]==space_dims[0]-1 && made_5d.size()>0){
				_Fid_made.push_back(made_5d);
				//std::cout<<"\t\t\tLength of made_3d: " <<made_3d.size() <<"\n";
				made_5d.clear();
			}
			if(cart2[5]==space_dims[5]-1 && cart2[4]==space_dims[4]-1 && cart2[3]==space_dims[3]-1 && cart2[2]==space_dims[2]-1 && cart2[1]==space_dims[1]-1 && cart2[0]==space_dims[0]-1){
				//std::cout<<"\tLength of made: " <<_Fid_made.size() <<"\n";
			}
		}

		while(cart.GetNextCombination()){
			/*std::cout<<"\t\tIndex: ";
			for(int i=0; i<5; i++){
				std::cout<<cart[i];
			}
			std::cout<<"\n";*/
			if(flags_->Plot_Fid(cart[5]) && (flags_->Flags::Sim() || _weight_[cart[0]]==_nweighted_)){
				if(_species_[cart[5]] == _ele_){
					//std::cout<<"\tFid for Electrons\n";
					if(_ecuts_[cart[4]] == _event_ && _top_[cart[1]]!=_mnone_ && _cut_[cart[2]]!=_no_cut_){
						//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n"; 
						sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]],_top_[cart[1]]);
						plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
						_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}else if(_ecuts_[cart[4]] != _event_ && _top_[cart[1]]==_mnone_){//No Event Selection
						if(_ecuts_[cart[4]]== _none_ && _cut_[cart[2]]==_no_cut_){// && _cut_[cart[2]]==_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
						}
						if(_ecuts_[cart[4]]== _sanity_ && _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
						}
						if(_ecuts_[cart[4]]== _fid_cut_ && flags_->Flags::Fid_Cut(0)&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
						}
						if(_ecuts_[cart[4]]== _sf_cut_ && flags_->Flags::SF_Cut()&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
						}
						if(_ecuts_[cart[4]]== _cc_cut_ && flags_->Flags::CC_Cut()&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
						}
						if(_ecuts_[cart[4]]== _ec_cut_ && flags_->Flags::EC_Cut()&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
						}
						if(_ecuts_[cart[4]]== _vertex_cut_ && flags_->Flags::Vertex_Cut()&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
						}
						if(_ecuts_[cart[4]]== _id_cut_ && flags_->Flags::ID_Cut()&& _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
						}
						if(_ecuts_[cart[4]]== _pid_ && _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
						}
					}
				}else{
					if(cart[4]<(std::distance(std::begin(_hcuts_), std::end(_hcuts_)))){//Make sure we aren't going over hcuts
						//std::cout<<"\tFid for Hadrons\n";
						if(_hcuts_[cart[4]] == _event_ && _top_[cart[1]]!=_mnone_ && _cut_[cart[2]]!=_no_cut_){
							//std::cout<<"\tFid Event for Hadrons\n";
							//std::cout<<"*Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"Fid_%s_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]],_top_[cart[1]]);
							plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
							_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
						}else if(_hcuts_[cart[4]] != _event_ && _top_[cart[1]]==_mnone_){//No Event Selection
							//std::cout<<"\tFid PID Hadrons\n";
							if(_hcuts_[cart[4]]== _none_ && _cut_[cart[2]]==_no_cut_ ){
							//	std::cout<<"\t\tFid None Hadrons\n";
								//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
							}
							if(_hcuts_[cart[4]]== _sanity_ && _cut_[cart[2]]!=_no_cut_){
								//std::cout<<"\t\tFid Sanity Hadrons\n";
								//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
							}
							if(_hcuts_[cart[4]]== _fid_cut_ && flags_->Flags::Fid_Cut(cart[5]) && _cut_[cart[2]]!=_no_cut_){
								//std::cout<<"\t\tFid Fid Hadrons\n";
								//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
							}
							if(_hcuts_[cart[4]]== _delta_cut_ && flags_->Flags::Delta_Cut(cart[5]) && _cut_[cart[2]]!=_no_cut_){
								//std::cout<<"\t\tFid delta Hadrons\n";
								//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
							}
							if(_hcuts_[cart[4]]== _pid_  && _cut_[cart[2]]!=_no_cut_){
								//std::cout<<"Fid Hist idx: " <<_Fid_hist.size() <<" " <<plot_5d.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"Fid_%s_%s_%s_%s_%s",_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_weight_[cart[0]]);
								plot_1d.push_back(new TH2F(hname,hname,_fid_xbin_,_fid_xmin_,_fid_xmax_,_fid_ybin_,_fid_ymin_,_fid_ymax_));
								_Fid_made[cart[5]][cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
							}
						}
					}
				}
			}
			if(cart[0] == space_dims[0]-1 && plot_1d.size()>0){
				plot_2d.push_back(plot_1d);
				plot_1d.clear();
			}
			if(cart[1] == space_dims[1]-1 && cart[0] == space_dims[0]-1 && plot_2d.size()>0){
				plot_3d.push_back(plot_2d);
				plot_2d.clear();
			}
			if(cart[2] == space_dims[2]-1 && cart[1] == space_dims[1]-1 && cart[0] == space_dims[0]-1 && plot_3d.size()>0){
				plot_4d.push_back(plot_3d);
				plot_3d.clear();
			}
			if(cart[3] == space_dims[3]-1 && cart[2] == space_dims[2]-1 && cart[1] == space_dims[1]-1 && cart[0] == space_dims[0]-1 && plot_4d.size()>0){
				plot_5d.push_back(plot_4d);
				plot_4d.clear();
			}
			if(cart[4] == space_dims[4]-1 && cart[3] == space_dims[3]-1 && cart[2] == space_dims[2]-1 && cart[1] == space_dims[1]-1 && cart[0] == space_dims[0]-1 && plot_5d.size()>0){
				_Fid_hist.push_back(plot_5d);
				plot_5d.clear();
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
		case 0:
			idx = idx_tmp;
		break;
		case 1:
			if(flags_->Flags::Fid_Cut(hadron_)){
				idx = idx_tmp; 
			}
		break;
		case 2:
			if(flags_->Flags::Delta_Cut(hadron_)){
				if(flags_->Flags::Fid_Cut(hadron_)){
					idx = idx_tmp; 
				}else{
					idx = idx_tmp-1; 
				}
			}
		break;
		case 3:
			if(flags_->Flags::Delta_Cut(hadron_)){
				if(flags_->Flags::Fid_Cut(hadron_)){
					idx = idx_tmp; 
				}else{
					idx = idx_tmp-1; 
				}
			}
		break;
		case 4:
			if(flags_->Flags::Delta_Cut(hadron_)){
				if(flags_->Flags::Fid_Cut(hadron_)){
					idx = idx_tmp; 
				}else{
					idx = idx_tmp-1; 
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
	


bool Histogram::Made_Fid_idx(const char* species_, const char* pcut_, const char* sector_, const char * cut_, const char * top_, const char * weight_){
	if(species_==_ele_){
		return _Fid_made[fun::species_idx(species_)][fun::ecut_idx(pcut_)][fun::sector_idx(sector_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::weight_idx(weight_)];
	}else{
		return _Fid_made[fun::species_idx(species_)][fun::hcut_idx(pcut_)][fun::sector_idx(sector_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::weight_idx(weight_)];
	}
}

std::vector<int> Histogram::Fid_cut_idx(const char* species_, const char* pcut_, const char* sector_, const char * cut_, const char * top_, const char * weight_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;// = {-1,-1,-1,-1,-1};
	int species_idx= -1; 
	int pcut_idx = -1;  
	int cut_idx = -1;
	for(int i=0; i<4; i++){
		if(species_ == _ele_){
			species_idx = 0; 
			pcut_idx = fun::ecut_idx(pcut_)+fun::ecut_offset(pcut_,flags_);
		}else if(species_ == _species_[i]){
			species_idx = i;
			pcut_idx = fun::hcut_idx(pcut_)+fun::hcut_offset(species_,pcut_,flags_);
		}
	}
	if(pcut_idx >= 0){
		//Species
		idx.push_back(species_idx);
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
			idx.push_back(fun::top_idx(top_));
		}else if(pcut_!=_event_ && top_==_mnone_){
			idx.push_back(0);
		}else{
			idx.push_back(-1);
		}
		//Weight
		if(flags_->Flags::Sim() || weight_==_nweighted_){
			idx.push_back(fun::weight_idx(weight_));
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
	}
	//std::cout<<"Fid idx";
	//fun::print_vector_idx(idx);
	return idx;
}

void Histogram::Fid_Fill(const char * species_, float theta_, float phi_, const char* pcut_, const char* sector_, const char *cut_, const char * top_, std::shared_ptr<Flags> flags_, float weight_){
	if(flags_->Plot_Fid(fun::species_idx(species_)) && !isnan(theta_) && !isnan(phi_)){
		//std::cout<<"Trying to fill Fid:" <<species_ <<" " <<pcut_ <<" "<<sector_ <<" " <<cut_ <<" " <<top_ <<" " <<weight_ <<"\n";
		int w_idx = -1;
		if(weight_!=1.0){
			w_idx = 1;
		}else{
			w_idx = 0;
		}
		std::vector<int> idx;
		if(flags_->Flags::Sim()){//Weighted Fiducial Plots
			idx = Histogram::Fid_cut_idx(species_,pcut_,sector_,cut_,top_,_weighted_,flags_);;
			if(Histogram::OK_Idx(idx)){
				if(Histogram::Made_Fid_idx(species_,pcut_,sector_,cut_,top_,_weighted_)){
					//fun::print_vector_idx(idx);
					_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(phi_,theta_,weight_);
				}
			}
			idx.clear();
			idx = Histogram::Fid_cut_idx(species_,pcut_,_sec_all_,cut_,top_,_weighted_,flags_);;
			if(Histogram::OK_Idx(idx)){
				if(Histogram::Made_Fid_idx(species_,pcut_,_sec_all_,cut_,top_,_weighted_)){
					//fun::print_vector_idx(idx);
					_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(phi_,theta_,weight_);
				}
			}
			idx.clear();
		}
		idx = Histogram::Fid_cut_idx(species_,pcut_,sector_,cut_,top_,_nweighted_,flags_);
		
		//if(pcut_==_event_){
			//std::cout<<"Trying to fill Fid: phi:" <<phi_ <<" theta:" <<theta_ <<" " <<species_ <<" " <<pcut_ <<" "<<sector_ <<" " <<cut_ <<" " <<top_ <<" " <<weight_ <<"\n";
		//}
		//fun::print_vector_idx(idx);
		if(Histogram::OK_Idx(idx)){
			if(Histogram::Made_Fid_idx(species_,pcut_,sector_,cut_,top_,_nweighted_)){
				//if(pcut_ == _event_){
				//	fun::print_vector_idx(idx);
				//	std::cout<<species_ <<" " <<pcut_ <<" " <<sector_ <<" " <<cut_ <<" " <<top_ <<" " <<_nweighted_ <<"\n";
				//}
				//fun::print_vector_idx(idx);
				_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(phi_,theta_,1.0);
			}
		}
		idx.clear();
		//Fill All Sectors
		idx = Histogram::Fid_cut_idx(species_,pcut_,_sec_all_,cut_,top_,_nweighted_,flags_);
		//fun::print_vector_idx(idx);
		if(Histogram::OK_Idx(idx)){
			if(Histogram::Made_Fid_idx(species_,pcut_,_sec_all_,cut_,top_,_nweighted_)){
				//if(pcut_ == _event_){
					//fun::print_vector_idx(idx);
				//}
				//fun::print_vector_idx(idx);
				_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Fill(phi_,theta_,1.0);
			}
		}
		idx.clear();
		//std::vector<int> idx = Histogram::Fid_cut_idx(species_,pcut_,cut_,top_,w_idx,flags_);
		
	}
}

void Histogram::Fid_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Plot_Fid(0) || flags_->Plot_Fid(1) || flags_->Plot_Fid(2) || flags_->Plot_Fid(3)){
		std::cout<<"Writing Fiducial Plots\n";
		char dir_name[100];
		TDirectory* dir_fid = _RootOutputFile->mkdir("Fiducial");
		dir_fid->cd();
		TDirectory* dir_fid_sub[4][2];
		for(int i=0; i<std::distance(std::begin(_species_), std::end(_species_)); i++){
			if(flags_->Flags::Plot_Fid(i)){
				sprintf(dir_name,"Fid %s PID Cuts",_species_[i]);
				dir_fid_sub[i][0] = dir_fid->mkdir(dir_name);
				sprintf(dir_name,"Fid Event %s Selection",_species_[i]);
				dir_fid_sub[i][1] = dir_fid->mkdir(dir_name);
			}
		}

		std::vector<long> space_dims(6);
		space_dims[5] = std::distance(std::begin(_species_), std::end(_species_));//{ele,pro,pip,pim}
		space_dims[4] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));//Electron Cuts (electrons have more cuts than hadrons)
		space_dims[3] = std::distance(std::begin(_sector_), std::end(_sector_));
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));//Cut and anti cut
		space_dims[1] = std::distance(std::begin(_top_), std::end(_top_));//{Topologies  + all combined + None}
		space_dims[0] = std::distance(std::begin(_weight_), std::end(_weight_));//Not Weight vs. Weighted
		CartesianGenerator cart(space_dims);

		std::vector<int> idx;

		while(cart.GetNextCombination()){
			if(flags_->Flags::Plot_Fid(cart[5])){
				if(_cut_[cart[2]]!=_no_cut_){
					if(_species_[cart[5]] == _ele_){
						idx = Histogram::Fid_cut_idx(_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[0]],flags_);
						//if(Histogram::OK_Idx(idx)){
						if(Histogram::OK_Idx(idx) && Histogram::Made_Fid_idx(_species_[cart[5]],_ecuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[0]])){
							//idx = Histogram::W2_cut_idx(_ecuts_[cart[5]],_cut_[cart[4]],_top_[cart[2]],cart[0],flags_,cart[1]);
							if(_ecuts_[cart[4]] == _event_){//Event Selection
								dir_fid_sub[cart[5]][1]->cd();
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("W (GeV)");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Q^{2} (GeV^{2}");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetOption("Colz");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
							}else{//PID Plots
								dir_fid_sub[cart[5]][0]->cd();
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("W (GeV)");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Q^{2} (GeV^{2}");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetOption("Colz");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
							}
						}
					}else{
						idx = Histogram::Fid_cut_idx(_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[0]],flags_);
						if(cart[4] < sizeof(_hcuts_) && Histogram::OK_Idx(idx) && Histogram::Made_Fid_idx(_species_[cart[5]],_hcuts_[cart[4]],_sector_[cart[3]],_cut_[cart[2]],_top_[cart[1]],_weight_[cart[0]])){
							if(_hcuts_[cart[4]]==_event_){//Event Selection
								dir_fid_sub[cart[5]][1]->cd();
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("W (GeV)");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Q^{2} (GeV^{2}");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetOption("Colz");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
							}else{//PID Plots
								dir_fid_sub[cart[5]][0]->cd();
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetXTitle("W (GeV)");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetYTitle("Q^{2} (GeV^{2}");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->SetOption("Colz");
								_Fid_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]][idx[5]]->Write();
							}
						}
					}
					idx.clear();
				}
			}
		}
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
		std::cout<<"Start Combo Stuff\n";
		std::vector<long> space_dims(5);
		space_dims[4] = std::distance(std::begin(_ecuts_),std::end(_ecuts_));//Ecuts
		space_dims[3] = std::distance(std::begin(_cut_),std::end(_cut_));//Cut applied?
		space_dims[2] = std::distance(std::begin(_top_),std::end(_top_));//Topology
		space_dims[1] = std::distance(std::begin(_sector_),std::end(_sector_));//sector
		space_dims[0] = Histogram::W_bins()+2;//W Dependence + exp range + full range
		char hname[100]; 
		CartesianGenerator cart2(space_dims);//For Made Histograms
		CartesianGenerator cart(space_dims);//For Histograms
		std::cout<<"Making 'Made' Array\n";
		//Make "Made" Vectors
		while(cart2.GetNextCombination()){
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
		}
		std::cout<<"Making Actual Histograms\n";
		//Make Actual Histograms
		while(cart.GetNextCombination()){
			if(_ecuts_[cart[4]]==_event_ && _top_[cart[2]]!=_mnone_){//Event Selection
				if(_cut_[cart[3]]!=_no_cut_){
					if(cart[0] < Histogram::W_bins()){//W Dependence
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"SF_%s_%s_%s_%s_W:%f-%f",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],Histogram::W_bot(cart[0]),Histogram::W_top(cart[0]));
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}else if(cart[0] == Histogram::W_bins()){//Experimental W range
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"SF_%s_%s_%s_%s_W:%s",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],"in_range");
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}else{//Full W range
					//	std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"SF_%s_%s_%s_%s_W:%s",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],"all");
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}
				}
			}else if(_ecuts_[cart[4]]!=_event_ && _top_[cart[2]]==_mnone_ && ((_cut_[cart[3]]!=_no_cut_ && _ecuts_[cart[4]]!=_none_) || ((_cut_[cart[3]]==_no_cut_ && _ecuts_[cart[4]]==_none_)))){//Electron ID Cuts
				if(fun::ecut_perform(_ecuts_[cart[4]],flags_)){
					if(cart[0] < Histogram::W_bins()){//W Dependence
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"SF_%s_%s_%s_%s_W:%f-%f",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],Histogram::W_bot(cart[0]),Histogram::W_top(cart[0]));
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}else if(cart[0] == Histogram::W_bins()){//Experimental W range
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"SF_%s_%s_%s_%s_W:%s",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],"in_range");
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
					}else{//Full W range
						//std::cout<<"SF Hist idx: " <<_SF_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n"; 
						sprintf(hname,"SF_%s_%s_%s_%s_W:%s",_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],"all");
						plot_1d.push_back(new TH2F(hname,hname,_sf_xbin_,_sf_xmin_,_sf_xmax_,_sf_ybin_,_sf_ymin_,_sf_ymax_));
						_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] = true;
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
		if(cut_==_no_cut_ && ecut_ == _none_ && top_==_mnone_){
			idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
		}else if(cut_ != _no_cut_ && ecut_!=_none_){
			idx.push_back(fun::ecut_idx(ecut_)+fun::ecut_offset(ecut_,flags_));
		}else{
			idx.push_back(-1);
		}
		if(cut_==_no_cut_){
			idx.push_back(0);
		}else{
			idx.push_back(fun::cut_idx(cut_));
		}
		if(ecut_==_event_ && top_ != _mnone_ ){
			idx.push_back(fun::top_idx(top_));
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
		//fun::print_vector_idx(idx);
		std::vector<int> idx2 = Histogram::SF_idx(ecut_,cut_,top_,_sec_all_,_W_var_,flags_,W_);
		if(Histogram::W_bin(W_) >=0){
			//std::cout<<"Looking at _SF_made idx: " <<fun::ecut_idx(ecut_) <<" " <<fun::cut_idx(cut_)<<" "<<fun::top_idx(top_)<<" "<<fun::sector_idx(sector_)<<" "<<Histogram::W_bin(W_)<<"\n";
			if(_SF_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::sector_idx(sector_)][Histogram::W_bin(W_)]){
				if(Histogram::OK_Idx(idx)){
					_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,sf_,weight_);
				}
				if(Histogram::OK_Idx(idx2)){
					_SF_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]][idx2[4]]->Fill(p_,sf_,weight_);
				}
			}
		}
		idx.clear();
		idx2.clear();
		//std::cout<<"Looking at _SF_made idx: " <<fun::ecut_idx(ecut_)<<" " <<fun::cut_idx(cut_)<<" "<<fun::top_idx(top_)<<" "<<fun::sector_idx(sector_)<<" "<<Histogram::W_bins()<<"\n";
		if(_SF_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::sector_idx(sector_)][Histogram::W_bins()] && cuts::in_range(W_)){
			idx = Histogram::SF_idx(ecut_,cut_,top_,sector_,_W_range_,flags_,W_);
			//fun::print_vector_idx(idx);
			idx2 = Histogram::SF_idx(ecut_,cut_,top_,_sec_all_,_W_range_,flags_,W_);
			if(Histogram::OK_Idx(idx)){
				_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,sf_,weight_);
			}
			if(Histogram::OK_Idx(idx2)){
				_SF_hist[idx2[0]][idx2[1]][idx2[2]][idx2[3]][idx2[4]]->Fill(p_,sf_,weight_);
			}
			idx.clear();
			idx2.clear();
		}
		//std::cout<<"Looking at _SF_made idx: " <<fun::ecut_idx(ecut_)<<" " <<fun::cut_idx(cut_)<<" "<<fun::top_idx(top_)<<" "<<fun::sector_idx(sector_)<<" "<<Histogram::W_bins()+1<<"\n";
		if(_SF_made[fun::ecut_idx(ecut_)][fun::cut_idx(cut_)][fun::top_idx(top_)][fun::sector_idx(sector_)][Histogram::W_bins()+1] && W_>0){
			idx = Histogram::SF_idx(ecut_,cut_,top_,sector_,_W_all_,flags_,W_);
			//fun::print_vector_idx(idx);
			idx2 = Histogram::SF_idx(ecut_,cut_,top_,_sec_all_,_W_all_,flags_,W_);
			if(Histogram::OK_Idx(idx)){
				_SF_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,sf_,weight_);
			}
			if(Histogram::OK_Idx(idx2)){
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
			if(_SF_made[cart[4]][cart[3]][cart[2]][cart[1]][cart[0]] && fun::ecut_perform(_ecuts_[cart[4]],flags_)){
				idx = Histogram::SF_idx(_ecuts_[cart[4]],_cut_[cart[3]],_top_[cart[2]],_sector_[cart[1]],_W_var_,flags_);
				idx[4] = (Histogram::W_bins()+1)-((Histogram::W_bins()+1)-cart[0]);
				//fun::print_vector_idx(idx);
				if(Histogram::OK_Idx(idx)){
					if(cart[0]<Histogram::W_bins()){
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

		std::vector<long> space_dims(3); //{ecuts,cut,top}
		space_dims[2] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		space_dims[1] = std::distance(std::begin(_cut_), std::end(_cut_));
		space_dims[0] = std::distance(std::begin(_top_), std::end(_top_));
		char hname[100];
		CartesianGenerator cart(space_dims);
		while(cart.GetNextCombination()){
			if(_ecuts_[cart[2]] == _none_ && _cut_[cart[1]] == _no_cut_ && _top_[cart[0]]==_mnone_){
				//std::cout<<"Making Vertex IDX: " <<_Vertex_hist.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
				sprintf(hname,"vertex_%s_%s",_none_,_no_cut_);
				plot_1d.push_back(new TH1F(hname,hname,_vertex_bin_,_vertex_min_,_vertex_max_));
			}else if(_ecuts_[cart[2]] != _none_ && _cut_[cart[1]] != _no_cut_){
				if(_ecuts_[cart[2]]==_event_ && _top_[cart[0]]!=_mnone_ && fun::top_perform(_top_[cart[0]],flags_)){
					//std::cout<<"Making Vertex IDX: " <<_Vertex_hist.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
					sprintf(hname,"vertex_%s_%s_%s",_ecuts_[cart[2]],_cut_[cart[1]],_top_[cart[0]]);
					plot_1d.push_back(new TH1F(hname,hname,_vertex_bin_,_vertex_min_,_vertex_max_));
				}else if(_ecuts_[cart[2]]!=_event_ && _top_[cart[0]]==_mnone_ && fun::ecut_perform(_ecuts_[cart[2]],flags_)){
					//std::cout<<"Making Vertex IDX: " <<_Vertex_hist.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
					sprintf(hname,"vertex_%s_%s_%s",_ecuts_[cart[2]],_cut_[cart[1]],_top_[cart[0]]);
					plot_1d.push_back(new TH1F(hname,hname,_vertex_bin_,_vertex_min_,_vertex_max_));
				}
			}
			if(cart[0] == space_dims[0]-1 ){
				if(plot_1d.size()>0){
					plot_2d.push_back(plot_1d);
					plot_1d.clear();
				}
				if(cart[1] == space_dims[1]-1){
					if(plot_2d.size()>0){
						_Vertex_hist.push_back(plot_2d);
						plot_2d.clear();
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::Vertex_idx(const char* ecut_, const char* cut_, const char* top_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
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
	return idx;
}

void Histogram::Vertex_Fill(float vz_, float weight_, const char* ecut_, const char* cut_, const char* top_, std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_Vertex()){
		std::vector<int> idx = Histogram::Vertex_idx(ecut_,cut_,top_,flags_);
		if(Histogram::OK_Idx(idx)){
			_Vertex_hist[idx[0]][idx[1]][idx[2]]->Fill(vz_, weight_);
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
		TDirectory* dir_vert_sub[std::distance(std::begin(_ecuts_), std::end(_ecuts_))][std::distance(std::begin(_cut_), std::end(_cut_))+1][std::distance(std::begin(_top_), std::end(_top_))+1];
		for(int i=0; i<std::distance(std::begin(_ecuts_), std::end(_ecuts_)); i++){
			if(fun::ecut_perform(_ecuts_[i],flags_)){
				sprintf(dir_name,"vertex_%s",_ecuts_[i]);
				dir_vert_sub[i][0][0] = dir_vert->mkdir(dir_name);
				if(_ecuts_[i] != _none_ ){
					for(int j=0; j<std::distance(std::begin(_cut_), std::end(_cut_)); j++){
						if(_cut_[j]!=_no_cut_){
							sprintf(dir_name,"vertex_%s_%s",_ecuts_[i],_cut_[j]);
							dir_vert_sub[i][j+1][0] = dir_vert_sub[i][0][0]->mkdir(dir_name);
							if(_ecuts_[i]==_event_){
								for(int k=0; k<std::distance(std::begin(_top_), std::end(_top_)); k++){
									if(_top_[k]!=_mnone_){
										sprintf(dir_name,"vertex_%s_%s_%s",_ecuts_[i],_cut_[j],_top_[k]);
										dir_vert_sub[i][j+1][k+1] = dir_vert_sub[i][j+1][0]->mkdir(dir_name);
									}
								}
							}
						}
					}
				}
			}
		}
		std::vector<long> space_dims(3); //{ecuts,cut,top}
		space_dims[2] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		space_dims[1] = std::distance(std::begin(_cut_), std::end(_cut_));
		space_dims[0] = std::distance(std::begin(_top_), std::end(_top_));
		CartesianGenerator cart(space_dims);
		std::vector<int> idx;
		while(cart.GetNextCombination()){
			if(_ecuts_[cart[2]]==_none_ && _cut_[cart[1]] == _no_cut_ && _top_[cart[0]]==_mnone_){
				idx = Histogram::Vertex_idx(_ecuts_[cart[2]],_cut_[cart[1]],_top_[cart[0]],flags_);
				if(Histogram::OK_Idx(idx)){
					dir_vert_sub[cart[2]][0][0]->cd();
					_Vertex_hist[idx[0]][idx[1]][idx[2]]->SetXTitle("Vz (nm)");
					_Vertex_hist[idx[0]][idx[1]][idx[2]]->SetYTitle("Yield");
					_Vertex_hist[idx[0]][idx[1]][idx[2]]->Write();
				}
				idx.clear();
			}else if(_ecuts_[cart[2]]!=_none_ && _cut_[cart[1]] != _no_cut_){
				if(_ecuts_[cart[2]]==_event_ && _top_[cart[0]]!=_mnone_){
					if(fun::top_perform(_top_[cart[0]],flags_)){
						idx = Histogram::Vertex_idx(_ecuts_[cart[2]],_cut_[cart[1]],_top_[cart[0]],flags_);
						if(Histogram::OK_Idx(idx)){
							dir_vert_sub[cart[2]][cart[1]+1][cart[0]+1]->cd();
							_Vertex_hist[idx[0]][idx[1]][idx[2]]->SetXTitle("Vz (nm)");
							_Vertex_hist[idx[0]][idx[1]][idx[2]]->SetYTitle("Yield");
							_Vertex_hist[idx[0]][idx[1]][idx[2]]->Write();
						}
						idx.clear();
					}
				}else if(_ecuts_[cart[2]]!=_event_ && _top_[cart[0]]==_mnone_){
					if(fun::ecut_perform(_ecuts_[cart[2]],flags_)){
						idx = Histogram::Vertex_idx(_ecuts_[cart[2]],_cut_[cart[1]],_top_[cart[0]],flags_);
						if(Histogram::OK_Idx(idx)){
							dir_vert_sub[cart[2]][cart[1]+1][0]->cd();
							_Vertex_hist[idx[0]][idx[1]][idx[2]]->SetXTitle("Vz (nm)");
							_Vertex_hist[idx[0]][idx[1]][idx[2]]->SetYTitle("Yield");
							_Vertex_hist[idx[0]][idx[1]][idx[2]]->Write();
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
		TH2F_ptr_4d plot_4d;
		TH2F_ptr_3d plot_3d;
		TH2F_ptr_2d plot_2d;
		TH2F_ptr_1d plot_1d;

		std::vector<long> space_dims(5);//{species,pcut,cut,top,w_dep}
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
		CartesianGenerator cart(space_dims);
		while(cart.GetNextCombination()){
			species_idx = cart[4];
			pcut_idx = cart[3];
			cut_idx = cart[2];
			top_idx = cart[1];
			w_idx = cart[0];
			if(_species_[species_idx]==_ele_){
				//std::cout<<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
				if(_ecuts_[pcut_idx]==_none_ && _cut_[cut_idx]==_no_cut_ && _top_[top_idx]==_mnone_){
					if(_W_dep_[w_idx]==_W_var_){
						for(int i=0; i<Histogram::W_bins(); i++){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_W:%f-%f",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],Histogram::W_bot(i),Histogram::W_top(i));
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}else if(_W_dep_[w_idx]==_W_range_){
						//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_);
						plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
					}
					else if(_W_dep_[w_idx]==_W_all_){
						//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_);
						plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
					}
				}else if(_ecuts_[pcut_idx]!=_none_ && _cut_[cut_idx]!=_no_cut_){
					if(_ecuts_[pcut_idx]==_event_ && _top_[top_idx]!=_mnone_ && fun::top_perform(_top_[top_idx],flags_)){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"delta_%s_%s_%s_%s_W:%f-%f",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],Histogram::W_bot(i),Histogram::W_top(i));
								plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
							}
						}else if(_W_dep_[w_idx]==_W_range_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
						else if(_W_dep_[w_idx]==_W_all_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}else if(_ecuts_[pcut_idx]!=_event_ && _top_[top_idx]==_mnone_ && fun::ecut_perform(_ecuts_[pcut_idx],flags_)){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"delta_%s_%s_%s_%s_W:%f-%f",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],Histogram::W_bot(i),Histogram::W_top(i));
								plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
							}
						}
						if(_W_dep_[w_idx]==_W_range_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}else if(_W_dep_[w_idx]==_W_all_){
						//	std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}
				}
			}else if(pcut_idx<std::distance(std::begin(_hcuts_), std::end(_hcuts_))){
				//std::cout<<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
				if(_hcuts_[pcut_idx]==_none_ && _cut_[cut_idx]==_no_cut_ && _top_[top_idx]==_mnone_){
					if(_W_dep_[w_idx]==_W_var_){
						for(int i=0; i<Histogram::W_bins(); i++){
						//	std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_W:%f-%f",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],Histogram::W_bot(i),Histogram::W_top(i));
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}else if(_W_dep_[w_idx]==_W_range_){
						//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_);
						plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
					}
					else if(_W_dep_[w_idx]==_W_all_){
						//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
						sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_);
						plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
					}
				}else if(_hcuts_[pcut_idx]!=_none_ && _cut_[cut_idx]!=_no_cut_){
					if(_hcuts_[pcut_idx]==_event_ && _top_[top_idx]!=_mnone_ && fun::top_perform(_top_[top_idx],flags_)){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
							//	std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"delta_%s_%s_%s_%s_W:%f-%f",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],Histogram::W_bot(i),Histogram::W_top(i));
								plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
							}
						}else if(_W_dep_[w_idx]==_W_range_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
						else if(_W_dep_[w_idx]==_W_all_){
						//	std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
					}else if(_hcuts_[pcut_idx]!=_event_ && _top_[top_idx]==_mnone_ && fun::hcut_perform(_species_[species_idx],_hcuts_[pcut_idx],flags_)){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
								sprintf(hname,"delta_%s_%s_%s_%s_W:%f-%f",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],Histogram::W_bot(i),Histogram::W_top(i));
								plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
							}
						}else if(_W_dep_[w_idx]==_W_range_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_);
							plot_1d.push_back(new TH2F(hname,hname,_delta_xbin_,_delta_xmin_,_delta_xmax_,_delta_ybin_,_delta_ymin_,_delta_ymax_));
						}
						else if(_W_dep_[w_idx]==_W_all_){
							//std::cout<<"Delta Histograms Idx: " <<_Delta_hist.size() <<" " <<plot_4d.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
							sprintf(hname,"delta_%s_%s_%s_%s_%s",_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_);
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
								_Delta_hist.push_back(plot_4d);
								plot_4d.clear();
							}
						}
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::Delta_idx(float W_, const char* species_, const char* pcut_, const char* cut_, const char* top_, const char* W_dep_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	idx.push_back(fun::species_idx(species_));
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
	return idx;
}


void Histogram::Delta_Fill(float p_, float dt_, float weight_, float W_, const char* species_, const char* pcut_, const char* cut_, const char* top_, std::shared_ptr<Flags> flags_ ){
	if(flags_->Flags::Plot_Delta(fun::species_idx(species_))){
		std::vector<int> idx;
		//std::cout<<"Filling Delta with p:" <<p_ <<" dt:" <<dt_ <<" " <<species_ <<" " <<pcut_ <<" " <<cut_  <<" " <<top_  <<" W:" <<W_ <<"\n";
		if(cuts::in_range(W_)){
			idx = Histogram::Delta_idx(W_,species_,pcut_,cut_,top_,_W_var_,flags_);
			//fun::print_vector_idx(idx); 
			if(Histogram::OK_Idx(idx)){
				_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,dt_,weight_);
			}
			idx.clear();
			idx = Histogram::Delta_idx(W_,species_,pcut_,cut_,top_,_W_range_,flags_);
			//fun::print_vector_idx(idx); 
			if(Histogram::OK_Idx(idx)){
				_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,dt_,weight_);
			}
			idx.clear();
		}
		idx = Histogram::Delta_idx(W_,species_,pcut_,cut_,top_,_W_all_,flags_);
		//fun::print_vector_idx(idx); 
		if(Histogram::OK_Idx(idx)){
			_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Fill(p_,dt_,weight_);
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
		TDirectory* dir_delta_sub[std::distance(std::begin(_species_), std::end(_species_))][std::distance(std::begin(_ecuts_), std::end(_ecuts_))+1][std::distance(std::begin(_cut_), std::end(_cut_))+1][std::distance(std::begin(_top_), std::end(_top_))+1][2+1];//Topology,cut,clean event, W dependence 
		//std::cout<<"\tMaking Directories\n";
		for(int i=0; i<std::distance(std::begin(_species_), std::end(_species_)); i++){
			if(flags_->Flags::Plot_Delta(i)){
				//std::cout<<"\tMaking Directory: " <<i<<" " <<0<<" " <<0<<" " <<0<<" " <<0 <<"\n";
				sprintf(dir_name,"delta_%s",_species_[i]);
				dir_delta_sub[i][0][0][0][0] = dir_delta->mkdir(dir_name);
				if(_species_[i]==_ele_){
					for(int j=0; j<std::distance(std::begin(_ecuts_), std::end(_ecuts_)); j++){
						if(_ecuts_[j]==_none_){
							sprintf(dir_name,"delta_%s_%s_%s",_species_[i],_ecuts_[j],_no_cut_);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0 <<"\n";
							dir_delta_sub[i][j+1][0][0][0] = dir_delta_sub[i][0][0][0][0]->mkdir(dir_name);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<1 <<"\n";
							sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_ecuts_[j],_no_cut_,_W_var_);
							dir_delta_sub[i][j+1][0][0][1] = dir_delta_sub[i][j+1][0][0][0]->mkdir(dir_name);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<2 <<"\n";
							sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_ecuts_[j],_no_cut_,_W_range_);
							dir_delta_sub[i][j+1][0][0][2] = dir_delta_sub[i][j+1][0][0][0]->mkdir(dir_name);
						}else if(fun::ecut_perform(_ecuts_[j],flags_)){
							sprintf(dir_name,"delta_%s_%s",_species_[i],_ecuts_[j]);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0 <<"\n";
							dir_delta_sub[i][j+1][0][0][0] = dir_delta_sub[i][0][0][0][0]->mkdir(dir_name);
							for(int k=0; k<std::distance(std::begin(_cut_), std::end(_cut_)); k++){
								if(_cut_[k]!=_no_cut_){
									//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<0<<" " <<0 <<"\n";
									sprintf(dir_name,"delta_%s_%s_%s",_species_[i],_ecuts_[j],_cut_[k]);
									dir_delta_sub[i][j+1][k+1][0][0] = dir_delta_sub[i][j+1][0][0][0]->mkdir(dir_name);
									for(int l=0; l<std::distance(std::begin(_top_), std::end(_top_)); l++){
										if(fun::top_perform(_top_[l],flags_)){
											if((_ecuts_[j]==_event_ && _top_[l]!=_mnone_) || (_ecuts_[j]!=_event_ && _top_[l]==_mnone_)){
												//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<0 <<"\n";
												sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_ecuts_[j],_cut_[k],_top_[l]);
												dir_delta_sub[i][j+1][k+1][l+1][0] = dir_delta_sub[i][j+1][k+1][0][0]->mkdir(dir_name);
												//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<1 <<"\n";
												sprintf(dir_name,"delta_%s_%s_%s_%s_%s",_species_[i],_ecuts_[j],_cut_[k],_top_[l],_W_var_);
												dir_delta_sub[i][j+1][k+1][l+1][1] = dir_delta_sub[i][j+1][k+1][l+1][0]->mkdir(dir_name);
												//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<2 <<"\n";
												sprintf(dir_name,"delta_%s_%s_%s_%s_%s",_species_[i],_ecuts_[j],_cut_[k],_top_[l],_W_range_);
												dir_delta_sub[i][j+1][k+1][l+1][2] = dir_delta_sub[i][j+1][k+1][l+1][0]->mkdir(dir_name);
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
						//	std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0 <<"\n";
							dir_delta_sub[i][j+1][0][0][0] = dir_delta_sub[i][0][0][0][0]->mkdir(dir_name);
							sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_hcuts_[j],_no_cut_,_W_var_);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<1 <<"\n";
							dir_delta_sub[i][j+1][0][0][1] = dir_delta_sub[i][j+1][0][0][0]->mkdir(dir_name);
							sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_hcuts_[j],_no_cut_,_W_range_);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0 <<" " <<0<<" " <<2 <<"\n";
							dir_delta_sub[i][j+1][0][0][2] = dir_delta_sub[i][j+1][0][0][0]->mkdir(dir_name);
						}else if(fun::hcut_perform(_species_[i],_hcuts_[j],flags_)){
							sprintf(dir_name,"delta_%s_%s",_species_[i],_hcuts_[j]);
							//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<0<<" " <<0<<" " <<0 <<"\n";
							dir_delta_sub[i][j+1][0][0][0] = dir_delta_sub[i][0][0][0][0]->mkdir(dir_name);
							for(int k=0; k<std::distance(std::begin(_cut_), std::end(_cut_)); k++){
								if(_cut_[k]!=_no_cut_){
									sprintf(dir_name,"delta_%s_%s_%s",_species_[i],_hcuts_[j],_cut_[k]);
									//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<0<<" " <<0 <<"\n";
									dir_delta_sub[i][j+1][k+1][0][0] = dir_delta_sub[i][j+1][0][0][0]->mkdir(dir_name);
									for(int l=0; l<std::distance(std::begin(_top_), std::end(_top_)); l++){
										if(fun::top_perform(_top_[l],flags_)){
											if((_hcuts_[j]==_event_ && _top_[l]!=_mnone_) || (_hcuts_[j]!=_event_ && _top_[l]==_mnone_)){
												sprintf(dir_name,"delta_%s_%s_%s_%s",_species_[i],_hcuts_[j],_cut_[k],_top_[l]);
												//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<0 <<"\n";
												dir_delta_sub[i][j+1][k+1][l+1][0] = dir_delta_sub[i][j+1][k+1][0][0]->mkdir(dir_name);
												sprintf(dir_name,"delta_%s_%s_%s_%s_%s",_species_[i],_hcuts_[j],_cut_[k],_top_[l],_W_var_);
												//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<1 <<"\n";
												dir_delta_sub[i][j+1][k+1][l+1][1] = dir_delta_sub[i][j+1][k+1][l+1][0]->mkdir(dir_name);
												sprintf(dir_name,"delta_%s_%s_%s_%s_%s",_species_[i],_hcuts_[j],_cut_[k],_top_[l],_W_range_);
												//std::cout<<"\tMaking Directory: " <<i<<" " <<j+1<<" " <<k+1<<" " <<l+1<<" " <<2 <<"\n";
												dir_delta_sub[i][j+1][k+1][l+1][2] = dir_delta_sub[i][j+1][k+1][l+1][0]->mkdir(dir_name);
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
		
		std::vector<long> space_dims(5);//{species,pcut,cut,top,w_dep}
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
		CartesianGenerator cart(space_dims);
		std::vector<int> idx;
		//std::cout<<"\tWriting Histograms\n";
		while(cart.GetNextCombination()){
			species_idx = cart[4];
			pcut_idx = cart[3];
			cut_idx = cart[2];
			top_idx = cart[1];
			w_idx = cart[0];
			//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
			if(flags_->Plot_Delta(species_idx)){
				if(_species_[species_idx]==_ele_){
					if(_ecuts_[pcut_idx]==_none_ && _cut_[cut_idx]==_no_cut_ && _top_[top_idx]==_mnone_){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								idx = Histogram::Delta_idx(Histogram::W_center(i),_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][0][0][1]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
								}
								idx.clear();
							}
						}else if(_W_dep_[w_idx]==_W_range_){
							idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
							if(Histogram::OK_Idx(idx)){
								//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
								dir_delta_sub[species_idx][pcut_idx+1][0][0][2]->cd();
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
							}
							idx.clear();
						}else if(_W_dep_[w_idx]==_W_all_){
							idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
							if(Histogram::OK_Idx(idx)){
								//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
								dir_delta_sub[species_idx][pcut_idx+1][0][0][2]->cd();
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
							}
							idx.clear();
						}
					}else if(_ecuts_[pcut_idx]!=_none_ && _cut_[cut_idx]!=_no_cut_){
						if(_ecuts_[pcut_idx]==_event_ && _top_[top_idx]!=_mnone_ && fun::top_perform(_top_[top_idx],flags_)){
							if(_W_dep_[w_idx]==_W_var_){
								for(int i=0; i<Histogram::W_bins(); i++){
									idx = Histogram::Delta_idx(Histogram::W_center(i),_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
									if(Histogram::OK_Idx(idx)){
										//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
										//fun::print_vector_idx(idx);
										//std::cout<<"\tPre directory: "<<species_idx<<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<1 <<"\n";
										dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][1]->cd();
										//std::cout<<"\tPost directory\n";
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
										//std::cout<<"\tPost Writing\n";
									}
									idx.clear();
								}
							}else if(_W_dep_[w_idx]==_W_range_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									//std::cout<<"\tPre directory: "<<species_idx<<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<2 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][2]->cd();
									//std::cout<<"\tPost directory\n";
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
									//std::cout<<"\tPost writing\n";
								}
								idx.clear();
							}else if(_W_dep_[w_idx]==_W_all_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									//std::cout<<"\tPre directory: "<<species_idx<<" " <<pcut_idx+1 <<" " <<cut_idx+1 <<" " <<top_idx+1 <<" " <<1 <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][2]->cd();
									//std::cout<<"\tPost directory\n";
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
									//std::cout<<"\tPost write\n";
								}
								idx.clear();
							}
						}else if(_ecuts_[pcut_idx]!=_event_ && _top_[top_idx]==_mnone_ && fun::ecut_perform(_ecuts_[pcut_idx],flags_)){
							if(_W_dep_[w_idx]==_W_var_){
								for(int i=0; i<Histogram::W_bins(); i++){
									idx = Histogram::Delta_idx(Histogram::W_center(i),_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
									if(Histogram::OK_Idx(idx)){
										//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
										dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][1]->cd();
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
									}
									idx.clear();
								}
							}else if(_W_dep_[w_idx]==_W_range_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
								}
								idx.clear();
							}else if(_W_dep_[w_idx]==_W_all_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_ecuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_ecuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
								}
								idx.clear();
							}
						}
					}
				}else if(pcut_idx<std::distance(std::begin(_hcuts_), std::end(_hcuts_))){
					if(_hcuts_[pcut_idx]==_none_ && _cut_[cut_idx]==_no_cut_ && _top_[top_idx]==_mnone_){
						if(_W_dep_[w_idx]==_W_var_){
							for(int i=0; i<Histogram::W_bins(); i++){
								idx = Histogram::Delta_idx(Histogram::W_center(i),_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][0][0][1]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
								}
								idx.clear();
							}
						}else if(_W_dep_[w_idx]==_W_range_){
							idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
							if(Histogram::OK_Idx(idx)){
								//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
								dir_delta_sub[species_idx][pcut_idx+1][0][0][2]->cd();
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
							}
							idx.clear();
						}else if(_W_dep_[w_idx]==_W_var_){
							idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
							if(Histogram::OK_Idx(idx)){
								//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
								dir_delta_sub[species_idx][pcut_idx+1][0][0][2]->cd();
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
								_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
							}
							idx.clear();
						}
					}else if(_hcuts_[pcut_idx]!=_none_ && _cut_[cut_idx]!=_no_cut_){
						if(_hcuts_[pcut_idx]==_event_ && _top_[top_idx]!=_mnone_ && fun::top_perform(_top_[top_idx],flags_)){
							if(_W_dep_[w_idx]==_W_var_){
								for(int i=0; i<Histogram::W_bins(); i++){
									idx = Histogram::Delta_idx(Histogram::W_center(i),_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
									if(Histogram::OK_Idx(idx)){
										//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
										//fun::print_vector_idx(idx);
										dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][1]->cd();
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
									}
									idx.clear();
								}
							}else if(_W_dep_[w_idx]==_W_range_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
								}
								idx.clear();
							}else if(_W_dep_[w_idx]==_W_all_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
								}
								idx.clear();
							}
						}else if(_hcuts_[pcut_idx]!=_event_ && _top_[top_idx]==_mnone_ && fun::hcut_perform(_species_[species_idx],_hcuts_[pcut_idx],flags_)){
							if(_W_dep_[w_idx]==_W_var_){
								for(int i=0; i<Histogram::W_bins(); i++){
									idx = Histogram::Delta_idx(Histogram::W_center(i),_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_var_,flags_);
									if(Histogram::OK_Idx(idx)){
										//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<" " <<i <<"\n";
										dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][1]->cd();
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
										_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
									}
									idx.clear();
								}
							}else if(_W_dep_[w_idx]==_W_range_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_range_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
								}
								idx.clear();
							}else if(_W_dep_[w_idx]==_W_all_){
								idx = Histogram::Delta_idx(Histogram::W_center(0),_species_[species_idx],_hcuts_[pcut_idx],_cut_[cut_idx],_top_[top_idx],_W_all_,flags_);
								if(Histogram::OK_Idx(idx)){
									//std::cout<<"Writing Here: " <<_species_[species_idx] <<" " <<_hcuts_[pcut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_W_dep_[w_idx] <<"\n";
									dir_delta_sub[species_idx][pcut_idx+1][cut_idx+1][top_idx+1][2]->cd();
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetXTitle("Momentum (GeV)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->SetYTitle("Delta T (ns)");
									_Delta_hist[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]]->Write();
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
		TH1F_ptr_3d plot_3d;
		TH1F_ptr_2d plot_2d;
		TH1F_ptr_1d plot_1d;

		std::vector<long> space_dims(5);
		space_dims[3] = std::distance(std::begin(_top_), std::end(_top_));//{ele,pro,pip,pim}
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));//Electron Cuts (electrons have more cuts than hadrons)
		space_dims[1] = std::distance(std::begin(_clean_event_),std::end(_clean_event_));//Clean event or not
		space_dims[0]	= (Histogram::W_bins()+2);//W Binning + All in range + All period
		//std::cout<<"W bins+2: " <<space_dims[0];
 		CartesianGenerator cart(space_dims);
		CartesianGenerator cart2(space_dims);
		char hname[100];//For naming histograms
		//std::cout<<"W Bins: " <<Histogram::W_bins() <<"\n";
		while(cart2.GetNextCombination()){
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
		}
		while(cart.GetNextCombination()){
			if(_top_[cart[3]]!=_mnone_ && _top_[cart[3]]!=_mall_ && _clean_event_[cart[1]]!=_isolated_){
				if(cart[0]< Histogram::W_bins()){
					//std::cout<<"MM idx: " <<_MM_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
					sprintf(hname,"MM_%s_%s_%s_W:%f-%f",_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],Histogram::W_bot(cart[0]),Histogram::W_top(cart[0]));
					plot_1d.push_back(new TH1F(hname,hname,_mm2_bin_[cart[3]],_mm2_min_[cart[3]],_mm2_max_[cart[3]]));
					_MM_made[cart[3]][cart[2]][cart[1]][cart[0]]=true;
				}else if(cart[0] == Histogram::W_bins()){
					//std::cout<<"MM idx: " <<_MM_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
					sprintf(hname,"MM_%s_%s_%s_W:%s",_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],"in_range");
					plot_1d.push_back(new TH1F(hname,hname,_mm2_bin_[cart[3]],_mm2_min_[cart[3]],_mm2_max_[cart[3]]));
					_MM_made[cart[3]][cart[2]][cart[1]][cart[0]]=true;
				}else if(cart[0] == Histogram::W_bins()+1){
					//std::cout<<"MM idx: " <<_MM_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
					sprintf(hname,"MM_%s_%s_%s_W:%s",_top_[cart[3]],_cut_[cart[2]],_clean_event_[cart[1]],"all");
					plot_1d.push_back(new TH1F(hname,hname,_mm2_bin_[cart[3]],_mm2_min_[cart[3]],_mm2_max_[cart[3]]));
					_MM_made[cart[3]][cart[2]][cart[1]][cart[0]]=true;
				}
			}
			if(cart[0] == space_dims[0]-1 && plot_1d.size()>0){
				plot_2d.push_back(plot_1d);
				plot_1d.clear();
			}
			if(cart[1] == space_dims[1]-1 && cart[0] == space_dims[0]-1 && plot_2d.size()>0){
				plot_3d.push_back(plot_2d);
				plot_2d.clear();
			}
			if(cart[2] == space_dims[2]-1 && cart[0] == space_dims[0]-1 && plot_3d.size()>0 && cart[1] == space_dims[1]-1){
				_MM_hist.push_back(plot_3d);
				plot_3d.clear();
			}
		}
	}
}

std::vector<int> Histogram::MM_idx(const char* top_, const char* cut_, const char * clean_, const char * W_dep_, float W_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx; 
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

void Histogram::MM_Fill(const char* top_, const char* cut_, const char * clean_,float MM_, float W_, float weight_, std::shared_ptr<Flags> flags_){
	//std::cout<<"\tTrying to fill MM top:" <<top_ <<" cut:" <<cut_ <<" clean:" <<clean_ <<" MM:" <<MM_ <<" W:" <<W_ <<"\n";
	std::vector<int> idx; 
	if(flags_->Flags::Plot_MM(fun::top_idx(top_))){
		if(top_ != _mnone_ && top_ != _mall_ && clean_!=_isolated_){
			//std::cout<<"W: " <<W_ <<"  => W Bin: " <<Histogram::W_bin(W_) <<"\n";
			if(Histogram::W_bin(W_)>=0){
				if(_MM_made[fun::top_idx(top_)][fun::cut_idx(cut_)][fun::clean_idx(clean_)][Histogram::W_bin(W_)]){//W Binning
					idx = Histogram::MM_idx(top_,cut_,clean_,_W_var_,W_,flags_);
					if(Histogram::OK_Idx(idx)){
						_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(MM_,weight_);
					}
					idx.clear();
				}
			}
			if(_MM_made[fun::top_idx(top_)][fun::cut_idx(cut_)][fun::clean_idx(clean_)][Histogram::W_bins()]){//W in range
				idx = Histogram::MM_idx(top_,cut_,clean_,_W_range_,W_,flags_);
				if(Histogram::OK_Idx(idx)){
					_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(MM_,weight_);
				}
				idx.clear();
			}
			if(_MM_made[fun::top_idx(top_)][fun::cut_idx(cut_)][fun::clean_idx(clean_)][Histogram::W_bins()+1]){//All W
				idx = Histogram::MM_idx(top_,cut_,clean_,_W_all_,W_,flags_);
				if(Histogram::OK_Idx(idx)){
					_MM_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(MM_,weight_);
				}
				idx.clear();
			}
		}
	}
}

void Histogram::MM_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_MM(0) || flags_->Flags::Plot_MM(1) || flags_->Flags::Plot_MM(2) || flags_->Flags::Plot_MM(3)){
		char dir_name[100];
		TDirectory* dir_mm = _RootOutputFile->mkdir("MM");
		dir_mm->cd();
		TDirectory* dir_mm_sub[4][std::distance(std::begin(_cut_), std::end(_cut_))+1][2+1][2+1];//Topology,cut,clean event, W dependence 
		for(int i=0; i<4; i++){
			sprintf(dir_name,"MM_%s",_top_[i]);
			dir_mm_sub[i][0][0][0] = dir_mm->mkdir(dir_name);
			for(int j=0; j<std::distance(std::begin(_cut_), std::end(_cut_)); j++){
				sprintf(dir_name,"MM_%s_%s",_top_[i],_cut_[j]);
				dir_mm_sub[i][1+j][0][0] = dir_mm_sub[i][0][0][0]->mkdir(dir_name);
				for(int k=0; k<(std::distance(std::begin(_clean_event_),std::end(_clean_event_))-1); k++){
					sprintf(dir_name,"MM_%s_%s_%s",_top_[i],_cut_[j],_clean_event_[k]);
					dir_mm_sub[i][j+1][k+1][0] = dir_mm_sub[i][j+1][0][0]->mkdir(dir_name); 
					sprintf(dir_name,"MM_%s_%s_W_dep",_top_[i],_cut_[j]);
					dir_mm_sub[i][j+1][k+1][1] = dir_mm_sub[i][j+1][k+1][0]->mkdir(dir_name);
					sprintf(dir_name,"MM_%s_%s_W_Range",_top_[i],_cut_[j]);
					dir_mm_sub[i][j+1][k+1][2] = dir_mm_sub[i][j+1][k+1][0]->mkdir(dir_name);
				}
				
			}
		}

		std::vector<long> space_dims(4);
		space_dims[3] = std::distance(std::begin(_top_), std::end(_top_));//{ele,pro,pip,pim}
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));//Electron Cuts (electrons have more cuts than hadrons)
		space_dims[1] = std::distance(std::begin(_clean_event_),std::end(_clean_event_));//Clean event or not
		space_dims[0]	= Histogram::W_bins()+2;//W Binning + All in range + All period 
		CartesianGenerator cart(space_dims);

		while(cart.GetNextCombination()){
			if(_top_[cart[3]]!=_mnone_ && _top_[cart[3]]!=_mall_ && _clean_event_[cart[1]]!=_isolated_){
				if(cart[0] >= Histogram::W_bins() && _MM_made[cart[3]][cart[2]][cart[1]][cart[0]]){
					dir_mm_sub[cart[3]][cart[2]+1][cart[1]+1][2]->cd();
					_MM_hist[cart[3]][cart[2]][cart[1]][cart[0]]->SetXTitle("MM (GeV^{2})");
					if(flags_->Flags::Sim()){
						_MM_hist[cart[3]][cart[2]][cart[1]][cart[0]]->SetXTitle("Weighted Yield");
					}else{
						_MM_hist[cart[3]][cart[2]][cart[1]][cart[0]]->SetXTitle("Yield");
					}
					_MM_hist[cart[3]][cart[2]][cart[1]][cart[0]]->Write();
				}else if(_MM_made[cart[3]][cart[2]][cart[1]][cart[0]]){
					dir_mm_sub[cart[3]][cart[2]+1][cart[1]+1][1]->cd();
					_MM_hist[cart[3]][cart[2]][cart[1]][cart[0]]->SetXTitle("MM (GeV^{2})");
					if(flags_->Flags::Sim()){
						_MM_hist[cart[3]][cart[2]][cart[1]][cart[0]]->SetXTitle("Weighted Yield");
					}else{
						_MM_hist[cart[3]][cart[2]][cart[1]][cart[0]]->SetXTitle("Yield");
					}
					_MM_hist[cart[3]][cart[2]][cart[1]][cart[0]]->Write();
				}
			}
		}
	}
}


 //*-------------------------------End MM Plot------------------------------*

 //*-------------------------------Start CC Eff Plot------------------------------*
void Histogram::CC_Eff_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_CC_Eff()){
		std::cout<<"Making CC Efficiency Plots\n";
		/*
		Histogram Splits
		- Electron Cuts
		- Cut Status
		- Topology
		*/
		TH2F_ptr_1d plot_1d_1;
		TH2F_ptr_2d plot_2d_1;
		TH2F_ptr_3d plot_3d_1;
		TH1F_ptr_1d plot_1d_2;
		TH1F_ptr_2d plot_2d_2;

		std::vector<long> space_dims(4);
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
			if(_ecuts_[ecut_idx] == _none_ && _cut_[cut_idx] == _no_cut_ && _top_[top_idx]==_none_){
			//	std::cout<<"CC Eff1 Idx: " <<_CC_Eff1_hist.size() <<" " <<plot_3d_1.size() <<" " <<plot_2d_1.size() <<" " <<plot_1d_1.size() <<"\n";
				sprintf(hname,"cc_eff_%s_%s_%s_%s",_ecuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx]);
				plot_1d_1.push_back(new TH2F(hname,hname,_cc_eff1_xbin_,_cc_eff1_xmin_,_cc_eff1_xmax_,_cc_eff1_ybin_,_cc_eff1_ymin_,_cc_eff1_ymax_));
				if(_sector_[sec_idx]==_sec_all_){
					//std::cout<<"\tCC Eff2 Idx: " <<_CC_Eff2_hist.size() <<" " <<plot_2d_2.size() <<" " <<plot_1d_2.size() <<"\n";
					sprintf(hname,"cc_eff_sec_yields_%s_%s_%s",_ecuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx]);
					plot_1d_2.push_back(new TH1F(hname,hname,_cc_eff2_xbin_,_cc_eff2_xmin_,_cc_eff2_xmax_));
				}
			}else if(_ecuts_[ecut_idx] != _none_ && _cut_[cut_idx] != _no_cut_){
				if(_ecuts_[ecut_idx] != _event_ && _top_[top_idx] == _mnone_){
					//std::cout<<"CC Eff1 Idx: " <<_CC_Eff1_hist.size() <<" " <<plot_3d_1.size() <<" " <<plot_2d_1.size() <<" " <<plot_1d_1.size() <<"\n";
					sprintf(hname,"cc_eff_%s_%s_%s_%s",_ecuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx]);
					plot_1d_1.push_back(new TH2F(hname,hname,_cc_eff1_xbin_,_cc_eff1_xmin_,_cc_eff1_xmax_,_cc_eff1_ybin_,_cc_eff1_ymin_,_cc_eff1_ymax_));
					if(_sector_[sec_idx]==_sec_all_){
						//std::cout<<"\tCC Eff2 Idx: " <<_CC_Eff2_hist.size() <<" " <<plot_2d_2.size() <<" " <<plot_1d_2.size() <<"\n";
						sprintf(hname,"cc_eff_sec_yields_%s_%s_%s",_ecuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx]);
						plot_1d_2.push_back(new TH1F(hname,hname,_cc_eff2_xbin_,_cc_eff2_xmin_,_cc_eff2_xmax_));
					}
				}else if(_ecuts_[ecut_idx] == _event_ && _top_[top_idx] != _mnone_){
					//std::cout<<"CC Eff1 Idx: " <<_CC_Eff1_hist.size() <<" " <<plot_3d_1.size() <<" " <<plot_2d_1.size() <<" " <<plot_1d_1.size() <<"\n";
					sprintf(hname,"cc_eff_%s_%s_%s_%s",_ecuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx]);
					plot_1d_1.push_back(new TH2F(hname,hname,_cc_eff1_xbin_,_cc_eff1_xmin_,_cc_eff1_xmax_,_cc_eff1_ybin_,_cc_eff1_ymin_,_cc_eff1_ymax_));
					if(_sector_[sec_idx]==_sec_all_){
						//std::cout<<"\tCC Eff2 Idx: " <<_CC_Eff2_hist.size() <<" " <<plot_2d_2.size() <<" " <<plot_1d_2.size() <<"\n";
						sprintf(hname,"cc_eff_sec_yields_%s_%s_%s",_ecuts_[ecut_idx],_cut_[cut_idx],_top_[top_idx]);
						plot_1d_2.push_back(new TH1F(hname,hname,_cc_eff2_xbin_,_cc_eff2_xmin_,_cc_eff2_xmax_));
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
							_CC_Eff1_hist.push_back(plot_3d_1);
							plot_3d_1.clear();
						}
					}
				}
				if(cart[2] == space_dims[2]-1 && plot_2d_2.size()>0 && _sector_[sec_idx]==_sec_all_){
					_CC_Eff2_hist.push_back(plot_2d_2);
						plot_2d_2.clear();
				}
			}
		}
	}
}

std::vector<int> Histogram::CC_Eff_idx(int which_, const char * ecut_, const char* cut_, const char* sector_, const char* top_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
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
	return idx;
}

void Histogram::CC_Eff_Fill(float p_, float theta_, float weight_, const char* ecut_, const char* cut_, const char* sector_, const char* top_, std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_CC_Eff()){
		std::vector<int> idx1 = Histogram::CC_Eff_idx(1,ecut_,cut_,sector_,top_,flags_);
		//fun::print_vector_idx(idx1);
		if(Histogram::OK_Idx(idx1)){
			_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->Fill(p_,theta_,weight_);
			idx1 = Histogram::CC_Eff_idx(1,ecut_,cut_,_sec_all_,top_,flags_);
			idx1.clear();
			if(Histogram::OK_Idx(idx1)){
				_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->Fill(p_,theta_,weight_);
			}
		}
		std::vector<int> idx2 = Histogram::CC_Eff_idx(2,ecut_,cut_,sector_,top_,flags_);
		//fun::print_vector_idx(idx2);
		if(Histogram::OK_Idx(idx2)){
			_CC_Eff2_hist[idx2[0]][idx2[1]][idx2[2]]->Fill(fun::sector_idx(sector_)+1,weight_);
		}
	}
}

void Histogram::CC_Eff_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Plot_CC_Eff()){
		std::cout<<"Writing CC Efficiency\n";
		char dirname[100];
		TDirectory* dir_cc_eff = _RootOutputFile->mkdir("CC Efficiency");
		TDirectory* dir_cc_eff_sub[std::distance(std::begin(_ecuts_), std::end(_ecuts_))][std::distance(std::begin(_cut_), std::end(_cut_))+1][std::distance(std::begin(_top_), std::end(_top_))+1];
		for(int i=0; i<std::distance(std::begin(_ecuts_), std::end(_ecuts_)); i++){
			if(fun::ecut_perform(_ecuts_[i],flags_)){
				//std::cout<<"making dir at " <<i <<" " <<0 <<" " <<0 <<" putting in dir_cc_eff\n";
				sprintf(dirname,"cc_eff_%s",_ecuts_[i]);
				dir_cc_eff_sub[i][0][0] = dir_cc_eff->mkdir(dirname);
				if(_ecuts_[i]==_event_){
					for(int j=0; j<std::distance(std::begin(_top_), std::end(_top_)); j++){
						if(_top_[j]!=_mnone_ && fun::top_perform(_top_[j],flags_)){
							//std::cout<<"\tmaking dir at " <<i <<" " <<0 <<" " <<j+1 <<" in " <<i <<" " <<0 <<" " <<0 <<"\n";
							sprintf(dirname,"cc_eff_%s_%s",_ecuts_[i],_top_[j]);
							dir_cc_eff_sub[i][0][j+1] = dir_cc_eff_sub[i][0][0]->mkdir(dirname);
							for(int k=0; k<std::distance(std::begin(_cut_), std::end(_cut_)); k++){
								if(_cut_[k]!=_no_cut_){
									//std::cout<<"\t\tmaking dir at " <<i <<" " <<k+1 <<" " <<j+1 <<" in " <<i <<" " <<0 <<" " <<j+1 <<"\n";
									sprintf(dirname,"cc_eff_%s_%s_%s",_ecuts_[i],_top_[j],_cut_[k]);
									dir_cc_eff_sub[i][k+1][j+1] = dir_cc_eff_sub[i][0][j+1]->mkdir(dirname);
								}
							}
						}
					}
				}else{
					for(int k=0; k<std::distance(std::begin(_cut_), std::end(_cut_)); k++){
						if(_cut_[k]!=_no_cut_ && _ecuts_[i]!=_none_){
							//std::cout<<"\tmaking dir at " <<i <<" " <<k+1 <<" " <<0 <<" in " <<i <<" " <<0 <<" " <<0 <<"\n";
							sprintf(dirname,"cc_eff_%s_%s",_ecuts_[i],_cut_[k]);
							dir_cc_eff_sub[i][k+1][0] = dir_cc_eff_sub[i][0][0]->mkdir(dirname);
						}
					}
				}
			}
		}
		std::vector<long> space_dims(4);
		space_dims[3] = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
		space_dims[2] = std::distance(std::begin(_cut_), std::end(_cut_));
		space_dims[1] = std::distance(std::begin(_sector_),std::end(_sector_));
		space_dims[0] = std::distance(std::begin(_top_), std::end(_top_));
		char hname[100];
		std::vector<int> idx1;
		std::vector<int> idx2;
		int ecut_idx;
		int cut_idx;
		int sec_idx;
		int top_idx;
		CartesianGenerator cart(space_dims);
		while(cart.GetNextCombination()){
			ecut_idx = cart[3];
			cut_idx = cart[2];
			sec_idx = cart[1];
			top_idx = cart[0];
			//std::cout<<"Looking at " <<_ecuts_[ecut_idx] <<" " <<_cut_[cut_idx] <<" " <<_top_[top_idx] <<" " <<_sector_[sec_idx] <<"\n";
			idx1 = Histogram::CC_Eff_idx(1,_ecuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx],flags_);
			idx2 = Histogram::CC_Eff_idx(2,_ecuts_[ecut_idx],_cut_[cut_idx],_sector_[sec_idx],_top_[top_idx],flags_);
			if(fun::ecut_perform(_ecuts_[ecut_idx],flags_)){
				if(_ecuts_[ecut_idx] == _none_ && _cut_[cut_idx] == _no_cut_ && _top_[top_idx]==_none_){
				//std::cout<<"Writing In Directory " <<fun::ecut_idx(_ecuts_[ecut_idx]) <<" 0 0\n"; 
				dir_cc_eff_sub[fun::ecut_idx(_ecuts_[ecut_idx])][0][0]->cd();
					if(Histogram::OK_Idx(idx1)){
						//fun::print_vector_idx(idx1);
						_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->SetXTitle("Momentum (GeV)");
						_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->SetYTitle("Theta (Degrees)");
						_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->Write();
					}
					if(_sector_[sec_idx]==_sec_all_){
						if(Histogram::OK_Idx(idx2)){
							//fun::print_vector_idx(idx2);
							_CC_Eff2_hist[idx2[0]][idx2[1]][idx2[2]]->SetXTitle("Sector");
							_CC_Eff2_hist[idx2[0]][idx2[1]][idx2[2]]->SetYTitle("Yield");
							_CC_Eff2_hist[idx2[0]][idx2[1]][idx2[2]]->Write();
						}
					}
				}else if(_ecuts_[ecut_idx] != _none_ && _cut_[cut_idx] != _no_cut_){
					if(_ecuts_[ecut_idx] != _event_ && _top_[top_idx] == _mnone_){
						//std::cout<<"Writing In Directory " <<fun::ecut_idx(_ecuts_[ecut_idx]) <<" " <<fun::cut_idx(_cut_[cut_idx])+1 <<" 0\n";
						dir_cc_eff_sub[fun::ecut_idx(_ecuts_[ecut_idx])][fun::cut_idx(_cut_[cut_idx])+1][0]->cd();
						if(Histogram::OK_Idx(idx1)){
							//fun::print_vector_idx(idx1);
							_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->SetXTitle("Momentum (GeV)");
							_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->SetYTitle("Theta (Degrees)");
							_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->Write();
						}
						if(_sector_[sec_idx]==_sec_all_){
							if(Histogram::OK_Idx(idx2)){
								//fun::print_vector_idx(idx2);
								_CC_Eff2_hist[idx2[0]][idx2[1]][idx2[2]]->SetXTitle("Sector");
								_CC_Eff2_hist[idx2[0]][idx2[1]][idx2[2]]->SetYTitle("Yield");
								_CC_Eff2_hist[idx2[0]][idx2[1]][idx2[2]]->Write();
							}
						}
					}else if(_ecuts_[ecut_idx] == _event_ && _top_[top_idx] != _mnone_){
						//std::cout<<"Writing In Directory " <<fun::ecut_idx(_ecuts_[ecut_idx]) <<" " <<fun::cut_idx(_cut_[cut_idx])+1 <<" " <<fun::top_idx(_top_[top_idx])+1 <<"\n";
						dir_cc_eff_sub[fun::ecut_idx(_ecuts_[ecut_idx])][fun::cut_idx(_cut_[cut_idx])+1][fun::top_idx(_top_[top_idx])+1]->cd();
						if(Histogram::OK_Idx(idx1)){
							//fun::print_vector_idx(idx1);
							_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->SetXTitle("Momentum (GeV)");
							_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->SetYTitle("Theta (Degrees)");
							_CC_Eff1_hist[idx1[0]][idx1[1]][idx1[2]][idx1[3]]->Write();
						}
						if(_sector_[sec_idx]==_sec_all_){
							if(Histogram::OK_Idx(idx2)){
								//fun::print_vector_idx(idx1);
								_CC_Eff2_hist[idx2[0]][idx2[1]][idx2[2]]->SetXTitle("Sector");
								_CC_Eff2_hist[idx2[0]][idx2[1]][idx2[2]]->SetYTitle("Yield");
								_CC_Eff2_hist[idx2[0]][idx2[1]][idx2[2]]->Write();
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

//*-------------------------------End CC Eff Plot------------------------------*

//*------------------------------- Start Friend Plot --------------------------*
std::vector<int> Histogram::Friend_Bin_Sizes(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Make_Friend()){
		std::vector<int> bin_sizes(7);
		bin_sizes[0] = Histogram::W_bins();
		bin_sizes[1] = std::distance(std::begin(_Q2_bins_), std::end(_Q2_bins_));
		bin_sizes[2] = _MM_bins_;
		bin_sizes[3] = _MM_bins_;
		bin_sizes[4] = _theta_bins_;
		bin_sizes[5] = _alpha_bins_;
		bin_sizes[6] = _phi_bins_;
		return bin_sizes;
	}
}

void Histogram::Friend_Make(std::shared_ptr<Flags> flags_){
	if(flags_->Make_Friend()){
		char hname[100];
		Int_t bins[7];
		std::vector<int> bins_vec = Histogram::Friend_Bin_Sizes(flags_);//A Bit redundant, but I believe it needs to be an Int array rather than a vector
		for(int k=0; k<7; k++){
			bins[k] = bins_vec[k];
		}
		for(int i=0; i<3; i++){//Topologies
			Double_t xmin[7] = {_W_min_,_Q2_min_,_MM_min_[i],_MM2_min_[i],_theta_min_,_alpha_min_,_phi_min_};
			Double_t xmax[7] = {_W_max_,_Q2_max_,_MM_max_[i],_MM2_max_[i],_theta_max_,_alpha_max_,_phi_max_};
			for(int j=0; j<5; j++){//Variable Sets
				sprintf(hname,"2#pi_off_proton_%s_%s",_var_names_[i],_top_[j]);
				_Friend[i][j] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
				_Friend[i][j]->GetAxis(1)->Set(5,_Q2_bins_);
				_Friend[i][j]->Sumw2();//Allow Weights and keep track of error
				if(flags_->Flags::Sim()){
					sprintf(hname,"Weight_2#pi_off_proton_%s_%s",_var_names_[i],_top_[j]);
					_Weight_Sum[i][j] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
					_Weight_Sum[i][j]->GetAxis(1)->Set(5,_Q2_bins_);
					_Weight_Sum[i][j]->Sumw2();
				}
			}
			if(flags_->Flags::Sim()){
				sprintf(hname,"Thrown_2#pi_off_proton_%s",_var_names_[i]);
				_Thrown[i] = new THnSparseD(hname,hname,7,bins,xmin,xmax);
				_Thrown[i]->GetAxis(1)->Set(5,_Q2_bins_);
				_Thrown[i]->Sumw2();//Allow Weights and keep track of error
			}
		}
		
	}
}


int Histogram::Friend_W_idx(float W_){
	return Histogram::W_bin(W_); 
}

int Histogram::Friend_Q2_idx(float Q2_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < std::distance(std::begin(_Q2_bins_), std::end(_Q2_bins_)); i++){//constants.hpp
    if(Q2_ < _Q2_bins_[i+1] && Q2_ >= _Q2_bins_[i]){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_MM_idx(float MM_, int var_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _MM_bins_; i++){//constants.hpp
    top = _MM_min_[var_] + (i+1)*((_MM_max_[var_]-_MM_min_[var_])/_MM_bins_);//constants.hpp
    bot = top - ((_MM_max_[var_]-_MM_min_[var_])/_MM_bins_); 
    if(MM_ < top && MM_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_MM2_idx(float MM_, int var_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _MM_bins_; i++){//constants.hpp
    top = _MM2_min_[var_] + (i+1)*((_MM2_max_[var_]-_MM2_min_[var_])/_MM_bins_);//constants.hpp
    bot = top - ((_MM2_max_[var_]-_MM2_min_[var_])/_MM_bins_); 
    if(MM_ < top && MM_ >= bot){
      j = i; 
    }
  }
  return j; 
}

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


std::vector<int>  Histogram::Friend_idx( float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_){
	std::vector<int> x(7); 
	int test = 0; 
	bool in_region = false; 
	//x[0] = top; //{pmiss,pipmiss,pimmiss,zeromiss,all}
	x[0] = Histogram::Friend_W_idx(W_);
	x[1] = Histogram::Friend_Q2_idx(Q2_);
	x[2] = Histogram::Friend_MM_idx(MM_,var_);
	x[3] = Histogram::Friend_MM2_idx(MM2_,var_);
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

void Histogram::Print_Friend_Bin(float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_, int var_){
	std::cout<<std::endl <<"--Printing Friend idx--" <<std::endl;
	std::cout<<"W: " <<W_ <<" in bin: " <<Histogram::Friend_W_idx(W_) <<std::endl;
	std::cout<<"Q2: " <<Q2_ <<" in bin: " <<Histogram::Friend_Q2_idx(Q2_) <<std::endl;
	std::cout<<"MM: " <<MM_ <<" in bin: " <<Histogram::Friend_MM_idx(MM_,var_) <<std::endl;
	std::cout<<"MM2: " <<MM2_ <<" in bin: " <<Histogram::Friend_MM2_idx(MM2_,var_) <<std::endl;
	std::cout<<"Theta: " <<theta_ <<" in bin: " <<Histogram::Friend_theta_idx(theta_) <<std::endl;
	std::cout<<"Alpha: " <<alpha_ <<" in bin: " <<Histogram::Friend_alpha_idx(alpha_) <<std::endl;
	std::cout<<"Phi: " <<phi_ <<" in bin: " <<Histogram::Friend_phi_idx(phi_) <<std::endl;
	//std::cout<<"Total Bin: " <<Histogram::Friend_idx(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_) <<std::endl;
}


void Histogram::Friend_Fill(const char* top_, float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_, bool thrown_, float weight_, std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Make_Friend() && fun::top_perform(top_,flags_)){
		//std::cout<<"W:" <<W_ <<" Q2:" <<Q2_ <<" MM:" <<MM_ <<" MM2:" <<MM2_ <<" theta:" <<theta_ <<" alpha:" <<alpha_ <<" phi:" <<phi_ <<" weight:" <<weight_ <<" thrown:" <<thrown_ <<"\n";
		if(!std::isnan(W_) && !std::isnan(Q2_) && !std::isnan(MM_) && !std::isnan(MM2_) && !std::isnan(theta_) && !std::isnan(alpha_) && !std::isnan(phi_) && !std::isnan(weight_)){
			if(Histogram::OK_Idx(Histogram::Friend_idx(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_))){
			//int *y = Friend_binning(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_);
			//std::cout<<std::endl <<"Y binning: " <<y;
				//Histogram::Print_Friend_Bin(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,var_);
				Double_t x[7] = { (double)W_, (double)Q2_, (double)MM_, (double)MM2_, (double)theta_, (double)alpha_, (double)phi_};
			//if(y[0]>=0){
				if(thrown_){
					_Thrown[var_]->Fill(x,weight_);
				}else{
					//std::cout<<"Filling Friend!\n";
					_Friend[var_][fun::top_idx(top_)]->Fill(x,weight_);
					//std::cout<<"Done filling friend\n";
					if(flags_->Flags::Sim()){
						_Weight_Sum[var_][fun::top_idx(top_)]->Fill(x,weight_*weight_);
					}
				}
				//std::cout<<std::endl <<"Filling Friend with " <<x <<" with weight " <<weight_;
			}
		}
	}
}

void Histogram::Friend_Write(std::shared_ptr<Flags> flags_){
	if(flags_->Flags::Make_Friend()){
		std::cout<<"Writing Friend\n";
		_SparseFile->cd();
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 5; j++){
				if(fun::top_perform(_top_[j],flags_)){
					_Friend[i][j]->Write();
					if(flags_->Flags::Sim()){
						_Weight_Sum[i][j]->Write();
					}
				}
			}
			if(flags_->Flags::Sim()){
				_Thrown[i]->Write();

			}
		}
	}
}

/*float Histogram::Friend_bin_reverse(int variable_, int bin_, int var_){//This only works for equally spaced bins. Will need rework when bins have varied sizes
	float val = NAN;
	float max = NAN;
	float min = NAN;
	int bins = -1; 
	switch(var_){
		/*case 0: 
			max = 4.5;
			min = -0.5;
			bins = _Friend_bins[var_]; 
		break;*/ //Top
		/*case 0:
			max = _W_max;
			min = _W_min;
			bins = _Friend_bins[var_] ; 
		break;//W
		case 1: 
			max = _Q2_max;
			min = _Q2_min;
			bins = _Friend_bins[var_]; 
		break;//Q2
		case 2: 
			max = _MM_max[channel_];
			min = _MM_min[channel_];
			bins = _Friend_bins[var_]; 
		break;//MM
		case 3: 
			max = _MM2_max[channel_];
			min = _MM2_min[channel_];
			bins = _Friend_bins[var_]; 
		break;//MM
		case 4: 
			max = _theta_max;
			min = _theta_min;
			bins = _Friend_bins[var_]; 
		break;//Theta
		case 5: 
			max = _alpha_max;
			min = _alpha_min;
			bins = _Friend_bins[var_]; 
		break;//Alpha
		case 6: 
			max = _phi_max;
			min = _phi_min;
			bins = _Friend_bins[var_]; 
		break;//phi
	}
	for(int i = 0; i < bins; i++){
		if(bin_ == i){
			val = ((max - min)/bins)*(i + 0.5) + min; 
		}
	}
	return val;
}*/





/*

//Sampling Fraction Cuts
void Histogram::Delta T_Make(std::shared_ptr<Environment> _envi){
	if(_envi->was_sf_plot()){
		char hname[100];
		std::vector<long> space_dims(5);
		space_dims[0] = 11; //Electron Cuts
		space_dims[1] = 30; //W binning
		space_dims[2] = 7;  //Sector
		space_dims[3] = 6;  //Topology
		space_dims[4] = 2; //cut v anti

		float top,bot; 

		CartesianGenerator cart(space_dims);

		while(cart.GetNextCombination()){
			if((cart[0] == 10 && cart[3] != 0) || (cart[0] != 10 && cart[3] == 0) && fun::hist_fitting(0,cart[0],cart[1],0,_envi->was_fit_type())){//Topology only matters for event selection cut
				if(cart[1] == 0 ){ //All W 
					sprintf(hname,"SF_%s_%s_%s_W:ALL_%s",eid_cut[cart[0]],cut_ver[cart[4]],sec_list[cart[2]],topologies[cart[3]]);
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]] = std::make_shared<TH2F>(hname,hname, SFxres, SFxmin, SFxmax, SFyres, SFymin, SFymax);
					//std::cout<<std::endl <<"Created plot: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3];

				}else{	//Specific W Bins
					top = Wbin_start + cart[1]*Wbin_res;
					bot = top - Wbin_res;
					sprintf(hname,"SF_%s_%s_%s_W:%f-%f_%s",eid_cut[cart[0]],cut_ver[cart[4]],sec_list[cart[2]],bot,top,topologies[cart[3]]);
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]] = std::make_shared<TH2F>(hname,hname, SFxres, SFxmin, SFxmax, SFyres, SFymin, SFymax);
					//std::cout<<std::endl <<"Created plot: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3];
				}
			}
		}
	}
}

void Histogram::SF_Fill(std::shared_ptr<Environment> _envi,int top, float p, float en, int cut, int cva, float W_, int sec, float weight_){
	if(_envi->was_sf_plot() && fun::hist_fitting(0,cut,Histogram::W_binning(W_),0,_envi->was_fit_type())){
		if(!std::isnan(p) && !std::isnan(en) && !std::isnan(W_)){
			SF_hist[cut][Histogram::W_binning(W_)][sec][top][cva]->Fill(p,en/p);
			SF_hist[cut][Histogram::W_binning(W_)][0][top][cva]->Fill(p,en/p);
			if(Histogram::W_binning(W_)!=0){
				if(_envi->was_sim()){
					SF_hist[cut][0][sec][top][cva]->Fill(p,en/p,weight_);
					SF_hist[cut][0][0][top][cva]->Fill(p,en/p,weight_);
				}else{
					SF_hist[cut][0][sec][top][cva]->Fill(p,en/p);
					SF_hist[cut][0][0][top][cva]->Fill(p,en/p);
				}
			}
		}
	}
}

void Histogram::SF_Write(std::shared_ptr<Environment> _envi){
	if(_envi->was_sf_plot()){
		std::cout<<"SF Plots: "; 
		TDirectory* dir_SF = _RootOutputFile->mkdir("SF Plots");
		dir_SF->cd();
		TDirectory* sf_dir[11];//Cut
		TDirectory* sf_dir_w[11];
		TDirectory* sf_dir_sec[11][7];
		TDirectory* sf_dir_top[11][5];
		//[8][6];// W binning, Sector, Topology
		char dir_name[100];
		for(int cut = 0 ; cut < 11; cut++){
			sprintf(dir_name,"SF_%s",eid_cut[cut]);
			sf_dir[cut] = dir_SF->mkdir(dir_name);
			//std::cout<<"Made Directory: " <<cut <<" 0 0 0" <<std::endl; 
			sprintf(dir_name,"SF_%s_%s",eid_cut[cut],W_dep_list[1]);
			sf_dir_w[cut] = sf_dir[cut]->mkdir(dir_name);
			//std::cout<<"Made Directory: " <<cut <<" " <<1 <<" 0 0" <<std::endl;
			for(int sec = 0; sec < 7 ; sec++){
				sprintf(dir_name,"SF_%s_%s",eid_cut[cut],sec_list[sec]);
				sf_dir_sec[cut][sec] = sf_dir[cut]->mkdir(dir_name);
				//std::cout<<"Made Directory: " <<cut <<" " <<0 <<" " <<sec <<" 0" <<std::endl;
			}
			for(int top = 1; top < 6; top++){
				sprintf(dir_name,"SF_%s_%s",eid_cut[cut],topologies[top]);
				sf_dir_top[cut][top] = sf_dir[cut]->mkdir(dir_name);
				//std::cout<<"Made Directory: " <<cut <<" " <<0 <<" " <<0 <<" "<<top <<std::endl;
			}
		}


		std::vector<long> space_dims(5);
		space_dims[0] = 11; //Electron Cuts
		space_dims[1] = 30; //W binning
		space_dims[2] = 7;  //Sector
		space_dims[3] = 6;  //Topology
		space_dims[4] = 2; //Cut v anti

		CartesianGenerator cart(space_dims);

		while(cart.GetNextCombination()){
			dir_SF->cd();
			//General Entry
			//std::cout<<"Curr Vals: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3]<<std::endl ;
			sf_dir[cart[0]]->cd();
			//std::cout<<"    Now In: " <<cart[0] /*<<" " <<"0" <<" " <<"0" <<" " <<"0"*///<<std::endl ; 
			//Want just cuts, so all Sectors, W, and use Combined Topology
/*			if(fun::hist_fitting(0,cart[0],cart[1],0,_envi->was_fit_type())){
				if(cart[2]==0 && cart[1]==0 && ((cart[0]!=10 && cart[3]==0)||(cart[0]==10 && cart[3]==5))){
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Sampling Fraction");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
				}
				//W binning
				if(cart[2]==0 && cart[1]!=0 && ((cart[0]!=10 && cart[3]==0)||(cart[0]==10 && cart[3]==5))){
					//std::cout<<"    Try In: " <<cart[0] <<" " <<"1" <<" " <<"0" <<" " <<"0"<<std::endl ;
					sf_dir_w[cart[0]]->cd();
					//std::cout<<"    Now In: " <<cart[0] <<" " <<"1" <<" " <<"0" <<" " <<"0"<<std::endl ; 
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Sampling Fraction");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
				}
				//Sector Binning
				if(cart[1]==0 && ((cart[0]!=10 && cart[3]==0)||(cart[0]==10 && cart[3]==5))){
					//std::cout<<"    Try In: " <<cart[0] <<" " <<"0" <<" " <<cart[2]+1 <<" " <<"0"<<std::endl ; 
					sf_dir_sec[cart[0]][cart[2]]->cd();
					//std::cout<<"    Now In: " <<cart[0] <<" " <<"0" <<" " <<cart[2]+1 <<" " <<"0"<<std::endl ; 
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Sampling Fraction");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
				}
				//Topology Binning
				if(cart[2]==0 && cart[1]==0 && (cart[0]==10 && cart[3]!=0)){
					//std::cout<<"    Now In: " <<cart[0] <<" " <<"0" <<" " <<"0" <<" " <<cart[3]<<std::endl ;
					sf_dir_top[cart[0]][cart[3]]->cd();
					//std::cout<<"    Now In: " <<cart[0] <<" " <<"0" <<" " <<"0" <<" " <<cart[3]<<std::endl ; 
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetXTitle("Momentum (GeV)");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetYTitle("Sampling Fraction");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->SetOption("COLZ");
					SF_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]]->Write();
				}
			}
		}
		std::cout<<"Done" <<std::endl;
	}
}


//Delta T Cuts
void Histogram::DT_Make(std::shared_ptr<Environment> _envi){
	char hname[100];
	std::vector<long> space_dims(6);
	space_dims[0] = 4;  //species
	space_dims[1] = 7;  //Cuts
	space_dims[2] = 30; //W Binning
	space_dims[3] = 7; //Sector
	space_dims[4] = 6; //topology
	space_dims[5] = 2; //cut v anti


	float bot,top;

	CartesianGenerator cart(space_dims);

	while(cart.GetNextCombination()){
		if(_envi->was_dt_plot(cart[0]) && fun::hist_fitting(cart[0],cart[1],cart[2],0,_envi->was_fit_type())){
			if((cart[1] == 6 && cart[4]!=0) || (cart[1]!=6 && cart[4] ==0)){
				if(cart[2] == 0){
					sprintf(hname,"%s_DeltaT_%s_%s_%s_W:ALL_%s",species[cart[0]],hid_cut[cart[1]],cut_ver[cart[5]],sec_list[cart[3]],topologies[cart[4]]);
					DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2F>(hname,hname, DTxres, DTxmin, DTxmax, DTyres, DTymin, DTymax);
					DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = true; 	
				}else if(cart[3]==0 && cart[4] == 0 && cart[1]!=6){//Looking at specific Cuts on fiducial and pre cut regimes 
					top = Wbin_start + cart[2]*Wbin_res;
					bot = top - Wbin_res;
					sprintf(hname,"%s_DeltaT_%s_%s_%s_W:%f-%f_%s",species[cart[0]],hid_cut[cart[1]],cut_ver[cart[5]],sec_list[cart[3]],bot,top,topologies[cart[4]]);
					DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH2F>(hname,hname, DTxres, DTxmin, DTxmax, DTyres, DTymin, DTymax);
					DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = true;
				}else{
					DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = false;
				}
			}
		}
	}
}

void Histogram::DT_Fill(std::shared_ptr<Environment> _envi,int top, int part, float p, float d, float t, float d0, float t0, int cut, int anti, float W_, int sec, float weight_){
	if(_envi->was_dt_plot(part) && fun::hist_fitting(part,cut,Histogram::W_binning(W_),0,_envi->was_fit_type())){
		if(!std::isnan(p) && !std::isnan(d) && !std::isnan(d) && !std::isnan(t) && !std::isnan(d0) && !std::isnan(t0) && !std::isnan(W_)){
			float dt = physics::delta_t(part, p, d, t, d0, t0);
			if(sec < 7){//If sector dependent
				if(top == 0 && cut!=6){//No event selection in the W variance. 
					if(DT_made_hist[part][cut][Histogram::W_binning(W_)][0][top][anti]){
						DT_hist[part][cut][Histogram::W_binning(W_)][0][top][anti]->Fill(p,dt,weight_);
					}else{
						std::cout<<"DT Would have segfaulted filling: " <<part <<" " <<cut <<" " <<Histogram::W_binning(W_) <<" " <<"0" <<" " <<top <<" " <<anti <<std::endl;
					}
				}
				if(DT_made_hist[part][cut][0][sec][top][anti]){
					DT_hist[part][cut][0][sec][top][anti]->Fill(p,dt,weight_);
				}else{
					std::cout<<"DT Would have segfaulted filling: " <<part <<" " <<cut <<" " <<0 <<" " <<sec <<" " <<top <<" " <<anti <<std::endl;
				}
				if(DT_made_hist[part][cut][0][0][top][anti]){
					DT_hist[part][cut][0][0][top][anti]->Fill(p,dt,weight_);
				}else{
					std::cout<<"DT Would have segfaulted filling: " <<part <<" " <<cut <<" " <<0 <<" " <<0 <<" " <<top <<" " <<anti <<std::endl;
				}
			}else{
				std::cout<<"DT would have had weird nonexistant sectors " <<std::endl;
			}
		}
	}
}
			
void Histogram::DT_Fill(std::shared_ptr<Environment> _envi,int top, int part, float p, float dt, int cut, int anti, float W_, int sec, float weight_){
	//std::cout<<"Filling DT Plot" <<std::endl;
	if(_envi->was_dt_plot(part) && fun::hist_fitting(part,cut,Histogram::W_binning(W_),0,_envi->was_fit_type())){
		if(!std::isnan(p) && !std::isnan(dt) && !std::isnan(W_)){
			if(sec < 7){
				if(cut!=6){//No event selection in the W variance. 
					//if(DT_made_hist[part][cut][Histogram::W_binning(W_)][sec][top][anti]){
					//	DT_hist[part][cut][Histogram::W_binning(W_)][sec][top][anti]->Fill(p,dt);
					//}else{
					//	std::cout<<"DT Would have segfaulted filling: " <<part <<" " <<cut <<" " <<Histogram::W_binning(W_) <<" " <<sec <<" " <<top <<" " <<anti <<std::endl;
					//}
					if(DT_made_hist[part][cut][Histogram::W_binning(W_)][0][top][anti]){
						DT_hist[part][cut][Histogram::W_binning(W_)][0][top][anti]->Fill(p,dt,weight_);
					}else{
						std::cout<<"DT Would have segfaulted filling: " <<part <<" " <<cut <<" " <<Histogram::W_binning(W_) <<" " <<sec <<" " <<top <<" " <<anti <<std::endl;
					}
				}
				if(DT_made_hist[part][cut][0][sec][top][anti]){
					DT_hist[part][cut][0][sec][top][anti]->Fill(p,dt,weight_);
				}else{
					std::cout<<"DT Would have segfaulted filling: " <<part <<" " <<cut <<" " <<0 <<" " <<sec <<" " <<top <<" " <<anti <<std::endl;
				}
				if(DT_made_hist[part][cut][0][0][top][anti]){
					DT_hist[part][cut][0][0][top][anti]->Fill(p,dt,weight_);
				}else{
					std::cout<<"DT Would have segfaulted filling: " <<part <<" " <<cut <<" " <<0 <<" " <<0 <<" " <<top <<" " <<anti <<std::endl;
				}
			}else{
				std::cout<<"DT would have had weird nonexistant sectors " <<std::endl;
			}
		}
	}
}

void Histogram::DT_Write(std::shared_ptr<Environment> _envi){
	bool do_dt_plots = false;
	for(int i = 0; i< 4; i++){
		if(_envi->was_dt_plot(i)){
			do_dt_plots = true;
		}
	}
	if(do_dt_plots){
		std::cout <<"Writing DT Plots: ";
		char dir_name[100]; 
		TDirectory* DT_plot = _RootOutputFile->mkdir("DT_plots");
		TDirectory* par_dt[4][9][2][8][6];//Particle, cut, W-dep, sector, topology
		//std::cout<<"Did I get here?" <<std::endl; 
		DT_plot->cd(); 
		for(int i = 0; i<4; i++){//Species
			if(_envi->was_dt_plot(i)){
				sprintf(dir_name,"%s_DT_plots",species[i]);
				par_dt[i][0][0][0][0]= DT_plot->mkdir(dir_name);
				DT_dir_made[i][0][0][0][0]=true;
				//std::cout<<"Made pointer:" <<i <<" 0 0 0 0" <<std::endl;
				for(int j = 1; j < 8; j++){//cut
					sprintf(dir_name,"%s_DT_%s",species[i],hid_cut[j-1]);
					par_dt[i][j][0][0][0] = par_dt[i][0][0][0][0]->mkdir(dir_name);
					DT_dir_made[i][j][0][0][0]=true;
					//std::cout<<"    Made Pointers" <<i <<j <<" 0 0 0"<<std::endl;
					//W Dependence
					sprintf(dir_name,"%s_DT_%s_%s",species[i],hid_cut[j-1],W_dep_list[1]);
					par_dt[i][j][1][0][0] = par_dt[i][j][0][0][0]->mkdir(dir_name);
					DT_dir_made[i][j][1][0][0]=true;
					//std::cout<<"    Made Pointers" <<i <<j <<" 1 0 0"<<std::endl;
					for(int k = 1; k < 8; k++){ //Sector 
						sprintf(dir_name,"%s_DT_%s_%s",species[i],hid_cut[j-1],sec_list[k-1]);
						par_dt[i][j][0][k][0] = par_dt[i][j][0][0][0]->mkdir(dir_name);
						DT_dir_made[i][j][0][k][0]=true;
						//std::cout<<"    Made Pointers" <<i <<j <<" 0 " <<k <<" 0"<<std::endl;
						//std::cout<<"Sector Pointer" <<std::endl;
						if(j == 7){//Event cut
							for(int l = 1; l < 6; l++){ //topology 
								sprintf(dir_name,"%s_DT_%s_%s",species[i],hid_cut[j-1],topologies[l]);
								par_dt[i][j][0][k][l] = par_dt[i][j][0][k][0]->mkdir(dir_name);
								DT_dir_made[i][j][0][k][l]=true;
								//std::cout<<"    Made Pointers" <<i <<j <<" 0 0 " <<k<<std::endl;
								//std::cout<<"Topology Pointer" <<std::endl;

							}
						}
					}
				}
			}
		}
		//std::cout<<"Made it through making directories" <<std::endl;

		std::vector<long> space_dims(6);
		space_dims[0] = 4;  //species
		space_dims[1] = 7;  //Cuts
		space_dims[2] = 30; //W Binning
		space_dims[3] = 7; //Sector
		space_dims[4] = 6; //topology
		space_dims[5] = 2; //Cut v anti

		CartesianGenerator cart(space_dims);
		while(cart.GetNextCombination()){
			if(_envi->was_dt_plot(cart[0]) && fun::hist_fitting(cart[0],cart[1],cart[2],0,_envi->was_fit_type())){
				if(DT_dir_made[0][0][0][0][0]){
					par_dt[cart[0]][0][0][0][0]->cd();//Main folder for particles
					if(DT_dir_made[cart[0]][cart[1]+1][0][0][0]){
						par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
						//std::cout <<"      Now Writing in " <<cart[0] <<" " <<cart[1]+1 <<" 0 0 0"<<std::endl;
						//All W, Sectors, and Combined Topology for Event selection, but still by Cut
						/*if(cart[2] ==0 &&  ((cart[1]!=6 && cart[4]==0)||(cart[1]==6 && cart[4]==5))){//got rid of cart[3]==0
							if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
								//std::cout<<"DT shouldn't have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
								DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("Momentum (GeV)");
								DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("Delta T (ns)");
								DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
								DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();	
							}else{
								std::cout<<"DT Would have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
							}		
						}*/
						
						//For W Range, but all sectors, combine topology for event selection, all sectors
		/*				if(_envi->was_W_dep_plot()){	
							if(cart[2]!=0 && cart[3] == 0 && ((cart[1]!=6 && cart[4]==0))){//||(cart[1]==6 && cart[4]==5))){ //Issue when trying to look at event selection for this. Not sure why, but getting rid of it solved it *shrug* 9/12/19
								//std::cout <<"          trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 1 0 0"<<std::endl;
								if(DT_dir_made[cart[0]][cart[1]+1][1][0][0]){
									par_dt[cart[0]][cart[1]+1][0][0][0]->cd();//Get into those CUTS
									par_dt[cart[0]][cart[1]+1][1][0][0]->cd();
									//std::cout <<"   We are writing" <<std::endl;
									if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
										//std::cout<<"DT shouldn't have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
										DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("Momentum (GeV)");
										DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("Delta T (ns)");
										DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
										DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();	
									}else{
										std::cout<<"DT Would have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
									}	
								}else{
									std::cout<<"cart values: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
									std::cout<<"DTWould have filled nonexistant dir on: " <<cart[0] <<" " <<cart[1]+1 <<" " <<1 <<" " <<0 <<" " <<0 <<" " <<0 <<std::endl;
								}
										
							}
						}
						//For Sector Range, but all W, combine topology for event selection
						if(cart[2] ==0 && ((cart[1]!=6 && cart[4]==0)||(cart[1]==6 && cart[4] > 0))){
							//std::cout <<"          Trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 0 " <<cart[3]+1 <<" 0"<<std::endl;
							if(DT_dir_made[cart[0]][cart[1]+1][0][cart[3]+1][cart[4]]){
								par_dt[cart[0]][cart[1]+1][0][cart[3]+1][cart[4]]->cd();//Get into those CUTS
								par_dt[cart[0]][cart[1]+1][0][cart[3]+1][cart[4]]->cd();
								//std::cout <<"  We are writing "<<std::endl;
								if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
									//std::cout<<"DT shouldn't have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
									DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("Momentum (GeV)");
									DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("Delta T (ns)");
									DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
									DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();	
								}else{
									std::cout<<"DT Would have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
								}
							}else{
								std::cout<<"cart values: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
								std::cout<<"DTWould have filled nonexistant dir on: " <<cart[0] <<" " <<cart[1]+1 <<" " <<0 <<" " <<cart[3]+1 <<" " <<0 <<" " <<0 <<std::endl;
							}		
						}
						//For topology Range, but all W, all sector
						/*if(cart[2] ==0 && cart[3] == 0 && cart[1] == 6 && cart[4]!=0){
							//std::cout <<"          Trying to Write in " <<cart[0] <<" " <<cart[1]+1 <<" 0 0 " <<cart[4]+1<<std::endl;
							if(DT_dir_made[cart[0]][cart[1]+1][0][cart[3]+1][cart[4]]){//Changed cart[4] + 1 to just cart[4] 5/5/20
								par_dt[cart[0]][cart[1]+1][0][cart[3]+1][cart[4]]->cd();//Get into those CUTS
								par_dt[cart[0]][cart[1]+1][0][cart[3]+1][cart[4]]->cd();
								//std::cout <<"  we are writing"<<std::endl;
								if(DT_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
									//std::cout<<"DT shouldn't have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
									DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetXTitle("Momentum (GeV)");
									DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetYTitle("Delta T (ns)");
									DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->SetOption("COLZ");
									DT_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]->Write();	
								}else{
									std::cout<<"DT Would have segfaulted Writing: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
								}
							}else{
								std::cout<<"cart values: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
								std::cout<<"DTWould have filled nonexistant dir on: " <<cart[0] <<" " <<cart[1]+1 <<" " <<cart[3]+1 <<" " <<0 <<" " <<0 <<" " <<0 <<std::endl;
							}		
						}*/
	/*				}else{
						std::cout<<"cart values: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
						std::cout<<"DTWould have filled nonexistant dir on: " <<cart[0] <<" " <<cart[1]+1 <<" " <<0 <<" " <<0 <<" " <<0 <<" " <<cart[4]+1 <<std::endl;
					}
								
				}else{
						std::cout<<"cart values: " <<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
						std::cout<<"DTWould have filled nonexistant dir on: " <<0 <<" " <<0 <<" " <<0 <<" " <<0 <<" " <<0 <<" " <<cart[4]+1 <<std::endl;
				}
			}
		}
		std::cout<<" Done" <<std::endl;
	}
}

//Min CC Cuts
void Histogram::CC_Make(std::shared_ptr<Environment> _envi){
	
	if(_envi->was_cc_plot() && !(_envi->was_sim())){
		char hname[100];
		std::vector<long> space_dims(6);
		space_dims[0] = 6;  //Sector
		space_dims[1] = 18; //Segment
		space_dims[2] = 11;  //Cut
		space_dims[3] = 4;  //Side of detector
		space_dims[4] = 6;  //Topology
		space_dims[5] = 2; //cut v anti

		CartesianGenerator cart(space_dims);

		while(cart.GetNextCombination()){
			if((cart[2] == 10 && cart[4]!=0) || (cart[2]!=10 && cart[4] ==0)){//Making sure event selection lines up with topologies
					sprintf(hname,"MinCC_%s_%s_%s_Seg%d_%s_%s",eid_cut[cart[2]],cut_ver[cart[5]],sec_list[cart[0]+1],cart[1]+1,CC_det_side[cart[3]],topologies[cart[4]]);
					CC_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]] = std::make_shared<TH1F>(hname,hname, MinCCres, MinCCmin, MinCCmax);
					CC_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]=true;//std::cout<<std::endl<<"Making CC plot: ";
					//std::cout<<" Made CC plot: "<<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
			}
			if(!CC_made_hist[cart[0]][cart[1]][cart[2]][cart[3]][cart[4]][cart[5]]){
				//std::cout<<"Didn't Make CC plot: "<<cart[0] <<" " <<cart[1] <<" " <<cart[2] <<" " <<cart[3] <<" " <<cart[4] <<" " <<cart[5] <<std::endl;
			}
		}
	}
}

void Histogram::CC_Fill(std::shared_ptr<Environment> _envi,int top, int sec, int segm, int nphe, int cut, int anti){
	//std::cout<<std::endl <<"top: " <<top <<" sec: " <<sec <<" segm: " <<segm <<" nphe: " <<nphe <<" cut: " <<cut <<" anti: " <<anti <<std::endl;
	//std::cout<<"filling plot";
	//std::cout<<std::endl <<sec <<" " <<detect::cc_segment(segm) <<" " <<cut <<" " <<detect::cc_lrc(segm) <<" " <<top <<" " <<anti <<std::endl;
	if(_envi->was_cc_plot() && !(_envi->was_sim())){
		if((cut == 10 && top!=0) || (cut!=10 && top ==0)){
			//if(((sec-1) < 6) && ((detect::cc_segment(segm)) <18)&& (cut < 11) && (detect::cc_lrc(segm)<4) && (top < 6) && (anti < 2)){
			if(CC_made_hist[sec-1][detect::cc_segment(segm)][cut][detect::cc_lrc(segm)][top][anti]){
				//CC_fill_hist[sec-1][detect::cc_segment(segm)][cut][detect::cc_lrc(segm)][top][anti]=true;
				//if(CC_fill_hist[sec-1][detect::cc_segment(segm)][cut][detect::cc_lrc(segm)][top][anti] && ){
					CC_hist[sec-1][detect::cc_segment(segm)][cut][detect::cc_lrc(segm)][top][anti]->Fill(nphe);
					CC_hist[sec-1][detect::cc_segment(segm)][cut][3][top][anti]->Fill(nphe);
				
				//}
				
			}else{
					std::cout<<std::endl <<"Would have segfaulted in CC filling";
					std::cout<<std::endl <<"for CC Plot: " <<sec-1 <<" " <<detect::cc_segment(segm) <<" " <<cut <<" " <<detect::cc_lrc(segm) <<" " <<top <<" " <<anti <<std::endl;
			}
		}else{
			//std::cout<<std::endl <<"You're trying to fill an event thing that doesn't correspond to the correct cut";
		}
	}

}
void Histogram::CC_Write(std::shared_ptr<Environment> _envi){
	if(_envi->was_cc_plot() && !(_envi->was_sim())){
		char dir_name[100]; 
		TDirectory* CC_plot = _RootOutputFile->mkdir("CC_plots");
		TDirectory* par_cc[6][12][5][6];
		for(int i = 0; i< 6; i++){//Sectors
			sprintf(dir_name,"CC_%s",sec_list[i+1]);
			par_cc[i][0][0][0] = CC_plot->mkdir(dir_name);//By sector
			//for(int j = 1; j < 19 ; j++){//Segments
				//sprintf(dir_name,"CC_%s_Segm%d",sec_list[i+1],j);
				//par_cc[i][j][0][0][0] = par_cc[i][0][0][0][0]->mkdir(dir_name);
				for(int k = 1; k < 5; k++){//Part of CC hit
					sprintf(dir_name,"CC_%s_%s",sec_list[i+1],CC_det_side[k-1]);
					par_cc[i][0][k][0] = par_cc[i][0][0][0]->mkdir(dir_name);
					for(int l = 1; l<12; l++){//EID Cuts
						sprintf(dir_name,"CC_%s_%s_%s",sec_list[i+1],CC_det_side[k-1],eid_cut[l-1]);
						par_cc[i][l][k][0] = par_cc[i][0][k][0]->mkdir(dir_name);
						if(l == 11){
							for(int m = 1; m<6; m++){//Topologies
								sprintf(dir_name,"CC_%s_%s_%s_%s",sec_list[i+1],CC_det_side[k-1],eid_cut[l-1],topologies[m]);
								par_cc[i][l][k][m] = par_cc[i][l][k][0]->mkdir(dir_name);
							}
						}
					}
				}
			//}
		}

		CC_plot->cd();
		for(int i = 0; i< 6; i++){//Sectors
			//std::cout<<" Trying to Enter: " <<i <<std::endl;
			par_cc[i][0][0][0]->cd();
			//std::cout<<" Did do it Enter: " <<i <<std::endl;
			for(int j = 0; j < 18 ; j++){//Segments
				//std::cout<<" Trying to Enter: " <<i <<" " <<j <<std::endl;
				//par_cc[i][0][0][0]->cd();
				//std::cout<<" Did do it Enter: " <<i <<" " <<j <<std::endl;
				for(int k = 1; k < 5; k++){//Part of CC hit
					//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" 0 " <<k <<std::endl;
					par_cc[i][0][k][0]->cd();
					//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" 0 " <<k <<std::endl;
					for(int l = 1; l<12; l++){//EID Cuts
						//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<std::endl;
						par_cc[i][l][k][0]->cd();
						//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<std::endl;
						for(int m = 0; m<6; m++){
							for(int n = 0; n<2; n++){//cut v anti-cut
								if(l == 11 && m!=0){//Event selected
									//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
									par_cc[i][l][k][m]->cd();
									//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
									CC_hist[i][j][l-1][k-1][m][n]->SetXTitle("num photoelectrons");
									CC_hist[i][j][l-1][k-1][m][n]->SetYTitle("Counts");
									CC_hist[i][j][l-1][k-1][m][n]->Write();
								}else if( l != 11 && m==0){//Not event selected
									//std::cout<<" Trying to Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
									par_cc[i][l][k][m]->cd();
									//std::cout<<" Did do it Enter: " <<i <<" " <<j <<" "<<l <<" " <<k <<" " <<m <<std::endl;
									CC_hist[i][j][l-1][k-1][m][n]->SetXTitle("num photoelectrons");
									CC_hist[i][j][l-1][k-1][m][n]->SetYTitle("Counts");
									CC_hist[i][j][l-1][k-1][m][n]->Write();
								}
							}
						}
					}
				}
			}
		}
	}

	

}
//Missing Mass Cuts
void Histogram::MM_Make(std::shared_ptr<Environment> _envi){
	bool do_mm_plots = false;
	for(int i = 0; i< 4; i++){
		if(_envi->was_MM_plot(i)){
			do_mm_plots = true;
		}
	}
	if(do_mm_plots){
		char hname[100];
		std::vector<long> space_dims(4);
		space_dims[0] = 4;  //Topology
		space_dims[1] = 3;  //Cuts {pre, cut, anti}
		space_dims[2] = 2; //squared vs linear
		space_dims[3] = 2; //Fitting vs. not fitting plots

		float MMmin = NAN; 
		float MMmax = NAN; 

		CartesianGenerator cart(space_dims);

		while(cart.GetNextCombination()){
			if(_envi->was_MM_plot(cart[0])){
				sprintf(hname,"%s_MM_%s_%s_%s",topologies[cart[0]+1],basic_cut[cart[1]],MM_sq[cart[2]],fit_q[cart[3]]);
				MMmin = MMxmin2[cart[0]];
				MMmax = MMxmax2[cart[0]];
				MM_hist[cart[0]][cart[1]][cart[2]][cart[3]] = std::make_shared<TH1F>(hname,hname,MMxres,MMmin,MMmax);
			}
		}
	}
}
void Histogram::MM_Fill(std::shared_ptr<Environment> _envi,int top, float mm, int cut, int square,bool fit, float weight_){
	if(_envi->was_MM_plot(top)){
		if(!std::isnan(mm)){
			MM_hist[top][cut][square][1]->Fill(mm,weight_);
			if(fit){
				MM_hist[top][cut][square][0]->Fill(mm,weight_);
			}
		}
	}
}
void Histogram::MM_Write(std::shared_ptr<Environment> _envi){
	bool do_mm_plots = false;
	for(int i = 0; i< 4; i++){
		if(_envi->was_MM_plot(i)){
			do_mm_plots = true;
		}
	}
	if(do_mm_plots){
		char dir_name[100];
		TDirectory * MM_plot = _RootOutputFile->mkdir("MM plots");
		MM_plot->cd();
		TDirectory * MM_dir[4][3][2];//top, cut, fit
		for(int k = 0; k<2; k++){
			for(int i = 0; i < 4; i++ ){
				if(_envi->was_MM_plot(i)){
					sprintf(dir_name,"MM %s %s",topologies[i+1],fit_q[k]);
					MM_dir[i][0][k] = MM_plot->mkdir(dir_name);
					for(int j = 1; j <3; j++){
						sprintf(dir_name,"MM %s %s %s",topologies[i+1],MM_sq[j-1],fit_q[k]);
						MM_dir[i][j][k] = MM_dir[i][0][k]->mkdir(dir_name);
					}
				}
			}
		}
		for(int p = 0; p<2 ; p++){
			for(int i = 0; i < 4; i++ ){//Topology without "All"
				if(_envi->was_MM_plot(i)){
					MM_dir[i][0][p]->cd();
					for(int j = 1; j <3; j++){//linear vs squared
						MM_dir[i][j][p]->cd();
						for( int k = 0 ; k<3; k++){//Cuts
							sprintf(dir_name,"MM %s (GeV %s)",MM_sq[j-1],MM_sq[j-1]);
							MM_hist[i][k][j-1][p]->SetXTitle(dir_name);
							MM_hist[i][k][j-1][p]->SetYTitle("Events");
							MM_hist[i][k][j-1][p]->Write();
						}
					}
				}
			}
		}
	}
	

}

void Histogram::XY_Make(std::shared_ptr<Environment> envi_){
	char hname[100];
	std::vector<long> space_dims(3);
	space_dims[0] = 3;  //Detectors {CC,SC,EC}
	space_dims[1] = 2;  //Simulation {Recon,Thrown}
	space_dims[2] = 4; //species

	CartesianGenerator cart(space_dims);
	if(envi_->was_xy_plot()){
		while(cart.GetNextCombination()){
			if(cart[1]==0){
				sprintf(hname,"XY_%s_%s",detectors[cart[0]+1],species[cart[2]]);
				XY_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,XYres,-XYmax[cart[0]],XYmax[cart[0]],XYres,-XYmax[cart[0]],XYmax[cart[0]]);
			}else if(envi_->was_sim() && cart[0]==0){
				sprintf(hname,"XY_%s_Thrown",species[cart[2]]);
				XY_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,XYres,-XYmax[cart[0]],XYmax[cart[0]],XYres,-XYmax[cart[0]],XYmax[cart[0]]);
			}
		}
	}
}

void Histogram::XY_Fill(std::shared_ptr<Environment> envi_, int species_, float x_, float y_, int detector_ , bool thrown_){
	//{1,2,3} -> {cc,sc,ed}
	if(envi_->was_xy_plot()){
		if(!std::isnan(x_) && !std::isnan(y_)){
			if(thrown_ && envi_->was_sim()){
				XY_hist[0][1][species_]->Fill(x_,y_);
			}else if(!thrown_){
				XY_hist[detector_][0][species_]->Fill(x_,y_);
			}
		}
	}
}

void Histogram::XY_Write(std::shared_ptr<Environment> envi_){
	if(envi_->was_xy_plot()){
		char dir_name[100];
		TDirectory * XY_plot = _RootOutputFile->mkdir("XY plots");
		XY_plot->cd();
		for(int l = 0; l < 4; l++){
			for(int k = 0; k<3; k++){
				XY_hist[k][0][l]->SetXTitle("X Position (mm)");//Check what these actual units are...
				XY_hist[k][0][l]->SetYTitle("Y Position (mm)");
				XY_hist[k][0][l]->Write();
			}
			if(envi_->was_sim()){
				XY_hist[0][1][l]->SetXTitle("X Position (mm)");//Check what these actual units are...
				XY_hist[0][1][l]->SetYTitle("Y Position (mm)");
				XY_hist[0][1][l]->Write();
			}
		}
	}
}

void Histogram::Fid_Det_Make(std::shared_ptr<Environment> envi_){
	char hname[100];
	std::vector<long> space_dims(3);//[3][4][7][2]
	space_dims[0] = 3;  //Detectors {CC,SC,EC}
	space_dims[1] = 4;  //Species
	space_dims[2] = 7; //All, Sector

	CartesianGenerator cart(space_dims);

	if(envi_->was_fid_det_plot()){
		while(cart.GetNextCombination()){
			if(cart[2]==0){
				if(envi_->was_sim()){
					sprintf(hname,"%sFid_%s_Sector:All_Sim",detectors[cart[0]+1],species[cart[1]]);
					Fid_Det_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
				}else{
					sprintf(hname,"%sFid_%s_Sector:All",detectors[cart[0]+1],species[cart[1]]);
					Fid_Det_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
				}
			}else{
				if(envi_->was_sim()){
					sprintf(hname,"%sFid_%s_Sector:%d_Sim",detectors[cart[0]+1],species[cart[1]],cart[2]);
					Fid_Det_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
				}else{
					sprintf(hname,"%sFid_%s_Sector:%d",detectors[cart[0]+1],species[cart[1]],cart[2]);
					Fid_Det_hist[cart[0]][cart[1]][cart[2]]= std::make_shared<TH2F>(hname,hname,FIDxres, FIDxmin, FIDxmax, FIDyres, FIDymin, FIDymax);
				}
			}
			
		}
	}
}

void Histogram::Fid_Det_Fill(std::shared_ptr<Environment> envi_, int species_, float theta_, float phi_, int sector_, int detector_){
	if(envi_->was_fid_det_plot()){
		if(!std::isnan(theta_) && !std::isnan(phi_)){
			Fid_Det_hist[detector_][species_][sector_]->Fill(phi_,theta_);
			Fid_Det_hist[detector_][species_][0]->Fill(phi_,theta_);
		}
	}
}

void Histogram::Fid_Det_Write(std::shared_ptr<Environment> envi_){
	char dir_name[100];
	if(envi_->was_fid_det_plot()){
		TDirectory * FidD_plot = _RootOutputFile->mkdir("Fid Detector plots");
		FidD_plot->cd();
		TDirectory * FidD_dir[4][4];
		for(int spe = 0; spe < 4; spe++){
			FidD_dir[spe][0] = FidD_plot->mkdir(species[spe]);
			FidD_dir[spe][0]->cd();
			for(int det = 0; det < 3; det++){
				sprintf(dir_name,"Det_%s_%s",detectors[det+1],species[spe]);
				FidD_dir[spe][det+1] = FidD_dir[spe][0]->mkdir(dir_name);
				FidD_dir[spe][det+1]->cd();
				for(int sec = 0; sec < 7; sec++){
				Fid_Det_hist[det][spe][sec]->SetXTitle("Phi (deg)");
				Fid_Det_hist[det][spe][sec]->SetYTitle("Theta (deg)");
				Fid_Det_hist[det][spe][sec]->Write();
				}
			}
		}
	}
}




void Histogram::Friend_Make(std::shared_ptr<Environment> _envi){
	//If we decide we want to fill this stuff
	if(_envi->was_Friend_plot()){
		for(int i=0; i<5; i++){
			for(int j=0; j<3; j++){
				char hname[100];
				sprintf(hname,"2#pi_off_proton_%s_%s",inter_had[j],topologies[i+1]);
				Double_t xmin[7] = {_W_min,_Q2_min,_MM_min[j],_MM2_min[j],_theta_min,_alpha_min,_phi_min};
				Double_t xmax[7] = {_W_max,_Q2_max,_MM_max[j],_MM2_max[j],_theta_max,_alpha_max,_phi_max};
				
				Friend[j][i] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin,xmax);
				//Friend[j][i] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,&{W_bins,_Q2_bins,MM1_bins,MM2_bins,theta_bins,alpha_bins,phi_bins});
				Friend[j][i]->GetAxis(1)->Set(5,_Q2_bins);
				//Friend[0]->Sumw2();//Must be called before error bars can be placed. 
				sprintf(hname,"2#pi_off_proton_#rho_%s",topologies[i+1]);
				Friend[1][i] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin2,xmax2);
				//Friend[1]->Sumw2();//Must be called before error bars can be placed. 
				sprintf(hname,"2#pi_off_proton_#Delta^{0}_%s",topologies[i+1]);
				Friend[2][i] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin3,xmax3);
				//Friend[2]->Sumw2();//Must be called before error bars can be placed. */ /*
			}
		}
	}
}

int Histogram::Friend_W_binning(float W_){
	int j = -1;
	float top, bot; 
	for(int i = 0; i < _Friend_bins[0]; i++){//constants.hpp
		top = _W_min + (i+1)*((_W_max-_W_min)/_Friend_bins[0]);//constants.hpp
	    bot = top - ((_W_max-_W_min)/_Friend_bins[0]); 
		if(W_ < top && W_ >= bot){
	    	j = i; 
	 	}
	}
	return j; 
}

int Histogram::Friend_Q2_binning(float Q2_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[1]; i++){//constants.hpp
    top = _Q2_min + (i+1)*((_Q2_max-_Q2_min)/_Friend_bins[1]);//constants.hpp
    bot = top - ((_Q2_max-_Q2_min)/_Friend_bins[1]); 
    if(Q2_ < top && Q2_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_MM_binning(float MM_, int chan){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[2]; i++){//constants.hpp
    top = _MM_min[chan] + (i+1)*((_MM_max[chan]-_MM_min[chan])/_Friend_bins[2]);//constants.hpp
    bot = top - ((_MM_max[chan]-_MM_min[chan])/_Friend_bins[2]); 
    if(MM_ < top && MM_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_MM2_binning(float MM_, int chan){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[3]; i++){//constants.hpp
    top = _MM2_min[chan] + (i+1)*((_MM2_max[chan]-_MM2_min[chan])/_Friend_bins[3]);//constants.hpp
    bot = top - ((_MM2_max[chan]-_MM2_min[chan])/_Friend_bins[3]); 
    if(MM_ < top && MM_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_theta_binning(float theta_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[4]; i++){//constants.hpp
    top = _theta_min + (i+1)*((_theta_max-_theta_min)/_Friend_bins[4]);//constants.hpp
    bot = top - ((_theta_max-_theta_min)/_Friend_bins[4]); 
    if(theta_ < top && theta_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_alpha_binning(float alpha_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[5]; i++){//constants.hpp
    top = _alpha_min + (i+1)*((_alpha_max-_alpha_min)/_Friend_bins[5]);//constants.hpp
    bot = top - ((_alpha_max-_alpha_min)/_Friend_bins[5]); 
    if(alpha_ < top && alpha_ >= bot){
      j = i; 
    }
  }
  return j; 
}

int Histogram::Friend_phi_binning(float phi_){
  int j = -1;
  float top, bot; 
  for(int i = 0; i < _Friend_bins[6]; i++){//constants.hpp
    top = _phi_min + (i+1)*((_phi_max-_phi_min)/_Friend_bins[6]);//constants.hpp
    bot = top - ((_phi_max-_phi_min)/_Friend_bins[6]); 
    if(phi_ < top && phi_ >= bot){
      j = i; 
    }
  }
  return j; 
}


int * Histogram::Friend_binning( float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int channel){
	int x[7]; 
	int test = 0; 
	bool in_region = false; 
	//x[0] = top; //{pmiss,pipmiss,pimmiss,zeromiss,all}
	x[0] = Friend_W_binning(W_);
	x[1] = Friend_Q2_binning(Q2_);
	x[2] = Friend_MM_binning(MM_,channel);
	x[3] = Friend_MM2_binning(MM2_,channel);
	x[4] = Friend_theta_binning(theta_);
	x[5] = Friend_alpha_binning(alpha_);
	x[6] = Friend_phi_binning(phi_);
	/*std::cout<<std::endl <<"Filling Friend in channel " <<channel <<std::endl; 
	std::cout<<"top = " <<top <<" bin: " <<x[0] <<std::endl;
	std::cout<<"W = " <<W_ <<" bin: " <<x[1] <<std::endl;
	std::cout<<"Q2 = " <<Q2_ <<" bin: " <<x[2] <<std::endl;
	std::cout<<"MM = " <<MM_ <<" bin: " <<x[3] <<std::endl;
	std::cout<<"theta = " <<theta_ <<" bin: " <<x[4] <<std::endl;
	std::cout<<"alpha = " <<alpha_ <<" bin: " <<x[5] <<std::endl;
	std::cout<<"phi = " <<phi_ <<" bin: " <<x[6] <<std::endl;
	*/   /*
	for(int i = 0; i < 7; i++){
		if(x[i] >= 0){
			test++;
		}
	}
	if(test < 7){
		x[0] = -1;
	}
	return x; 
}

void Histogram::Print_Friend_Bin(float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_, int chan_){
	std::cout<<std::endl <<"--Printing Friend Binning--" <<std::endl;
	std::cout<<"W: " <<W_ <<" in bin: " <<Friend_W_binning(W_) <<std::endl;
	std::cout<<"Q2: " <<Q2_ <<" in bin: " <<Friend_Q2_binning(Q2_) <<std::endl;
	std::cout<<"MM: " <<MM_ <<" in bin: " <<Friend_MM_binning(MM_,chan_) <<std::endl;
	std::cout<<"MM2: " <<MM2_ <<" in bin: " <<Friend_MM2_binning(MM2_,chan_) <<std::endl;
	std::cout<<"Theta: " <<theta_ <<" in bin: " <<Friend_theta_binning(theta_) <<std::endl;
	std::cout<<"Alpha: " <<alpha_ <<" in bin: " <<Friend_alpha_binning(alpha_) <<std::endl;
	std::cout<<"Phi: " <<phi_ <<" in bin: " <<Friend_phi_binning(phi_) <<std::endl;
	std::cout<<"Total Bin: " <<Friend_binning(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,0) <<std::endl;
}

bool Histogram::In_Friend_Bin(float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_, int chan_){
	bool pass = false;
	if(Friend_W_binning(W_) >= 0){
		if(Friend_Q2_binning(Q2_) >= 0){
			if(Friend_MM_binning(MM_,chan_) >= 0){
				if(Friend_MM2_binning(MM2_,chan_) >= 0){
					if(Friend_theta_binning(theta_) >= 0){
						if(Friend_alpha_binning(alpha_) >= 0){
							if(Friend_phi_binning(phi_) >= 0){
								pass = true;
							}
						}
					}
				}
			}
		}
	}
	return pass;
}

void Histogram::Friend_Fill(std::shared_ptr<Environment> _envi, int top_, float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int chan_, float weight_){
	if(_envi->was_Friend_plot()){
		//std::cout<<std::endl <<"W:" <<W_ <<" Q2:" <<Q2_ <<" MM:" <<MM_ <<" theta:" <<theta_ <<" alpha:" <<alpha_ <<" phi:" <<phi_ <<" weight:" <<weight_;
		if(!std::isnan(W_) && !std::isnan(Q2_) && !std::isnan(MM_) && !std::isnan(MM2_) && !std::isnan(theta_) && !std::isnan(alpha_) && !std::isnan(phi_) && !std::isnan(weight_)){
			if(Histogram::In_Friend_Bin(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,chan_)){
			//int *y = Friend_binning(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,chan_);
			//std::cout<<std::endl <<"Y binning: " <<y;
				//Histogram::Print_Friend_Bin(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,chan_);
				Double_t x[7] = { W_, Q2_, MM_, MM2_, theta_, alpha_, phi_};
			//if(y[0]>=0){
				Friend[chan_][top_]->Fill(x,weight_);
				//std::cout<<std::endl <<"Filling Friend with " <<x <<" with weight " <<weight_;
			}
		}
	}
}

void Histogram::Friend_Write(std::shared_ptr<Environment> _envi){
	/*if(_envi->Environment::was_Friend_plot()){
		char dir_name[100];
		TDirectory * Friend_plot = _RootOutputFile->mkdir("Friend plots");
		Friend_plot->cd();
		TH1D* MM_proj[3]; 
		TH1D* alpha_proj[3];
		TH1D* theta_proj[3];

		for(int i = 0; i< 3; i++){
			MM_proj[i] = Friend[i]->Projection(3);
			alpha_proj[i] = Friend[i]->Projection(5);
			theta_proj[i] = Friend[i]->Projection(4);
			Friend[i]->Write(); 
			MM_proj[i]->Draw();
			MM_proj[i]->Write();
			alpha_proj[i]->Draw();
			alpha_proj[i]->Write();
			theta_proj[i]->Draw();
			theta_proj[i]->Write();
		}
	}*/		/*
	if(_envi->Environment::was_Friend_plot()){
		SparseFile->cd();
		for(int i = 0; i< 3; i++){
			for(int j = 0; j<5; j++){
				Friend[i][j]->Write();
			}
		}
	}
}

float Histogram::Friend_bin_reverse(int var_, int bin_, int channel_){//This only works for equally spaced bins. Will need rework when bins have varied sizes
	float val = NAN;
	float max = NAN;
	float min = NAN;
	int bins = -1; 
	switch(var_){
		/*case 0: 
			max = 4.5;
			min = -0.5;
			bins = _Friend_bins[var_]; 
		break;*/ //Top
/*		case 0:
			max = _W_max;
			min = _W_min;
			bins = _Friend_bins[var_] ; 
		break;//W
		case 1: 
			max = _Q2_max;
			min = _Q2_min;
			bins = _Friend_bins[var_]; 
		break;//Q2
		case 2: 
			max = _MM_max[channel_];
			min = _MM_min[channel_];
			bins = _Friend_bins[var_]; 
		break;//MM
		case 3: 
			max = _MM2_max[channel_];
			min = _MM2_min[channel_];
			bins = _Friend_bins[var_]; 
		break;//MM
		case 4: 
			max = _theta_max;
			min = _theta_min;
			bins = _Friend_bins[var_]; 
		break;//Theta
		case 5: 
			max = _alpha_max;
			min = _alpha_min;
			bins = _Friend_bins[var_]; 
		break;//Alpha
		case 6: 
			max = _phi_max;
			min = _phi_min;
			bins = _Friend_bins[var_]; 
		break;//phi
	}
	for(int i = 0; i < bins; i++){
		if(bin_ == i){
			val = ((max - min)/bins)*(i + 0.5) + min; 
		}
	}
	return val;
}

void Histogram::Cross_Make(std::shared_ptr<Environment> envi_){
	char hname[100];
	
	/*
	0 - Pmiss Only
	1 - Pipmiss Only
	2 - Pimmiss Only
	3 - Zeromiss Only
	4 - Zeromiss + 3
	5 - Zeromiss + 2
	6 - Zeromiss + 1
	7 - Pmiss + Pipmiss
	8 - Pmiss + Pimmiss
	9 - Pipmiss + Pimmiss
	10 - No Zeromiss + 3
	11 - Multiples Pmiss
	12 - Multiples Pipmiss
	13 - Multiples Pimmiss
	14 - Multiples Zeromiss
	*/		/*
	if(envi_->was_cross_plot()){
		sprintf(hname,"Topology_Crossing_no_weights");
		Cross_hist[0] = std::make_shared<TH1F>(hname,hname, 15, -0.5, 14.5);
		sprintf(hname,"Topology_Crossing__weights");
		Cross_hist[1] = std::make_shared<TH1F>(hname,hname, 15, -0.5, 14.5);
	}
}

void Histogram::Cross_Fill(std::shared_ptr<Environment> envi_, int gevnt_[4], float weight_){
	bool top_pass[4] = {false,false,false,false};
	bool top_mult[4] = {false,false,false,false};
	int num_hadmiss = 0;
	int tot_evnt = 0; 
	for(int i = 0; i< 4; i++){
		if(gevnt_[i] > 0){
			top_pass[i] = true;
			if(i != 3){
				num_hadmiss += 1; 
			}
			if(gevnt_[i] > 1){
				top_mult[i] = true;
			}
		}
		tot_evnt += gevnt_[i];
	}
	if(envi_->was_cross_plot()){
		if(top_pass[0] && !top_pass[1] && !top_pass[2] && !top_pass[3] && !top_mult[0]){//Pmiss Only
			Cross_hist[0]->Fill(0.0);
			Cross_hist[1]->Fill(0.0,weight_);
		}else if(!top_pass[0] && top_pass[1] && !top_pass[2] && !top_pass[3] && !top_mult[1]){//Pipmiss only
			Cross_hist[0]->Fill(1.0);
			Cross_hist[1]->Fill(1.0,weight_);
		}else if(!top_pass[0] && !top_pass[1] && top_pass[2] && !top_pass[3] && !top_mult[2]){//Pimmiss only
			Cross_hist[0]->Fill(2.0);
			Cross_hist[1]->Fill(2.0,weight_);
		}else if(!top_pass[0] && !top_pass[1] && !top_pass[2] && top_pass[3] && !top_mult[3]){//Zeromiss only
			Cross_hist[0]->Fill(3.0);
			Cross_hist[1]->Fill(3.0,weight_);
		}else if(top_pass[0] && top_pass[1] && top_pass[2] && top_pass[3] && !top_mult[3]){//Zeromiss +3
			Cross_hist[0]->Fill(4.0);
			Cross_hist[1]->Fill(4.0,weight_);
		}else if((num_hadmiss == 2) && top_pass[3] && !top_mult[3]){//Zeromiss +2
			Cross_hist[0]->Fill(5.0);
			Cross_hist[1]->Fill(5.0,weight_);
		}else if((num_hadmiss == 1) && top_pass[3] && !top_mult[3]){//Zeromiss +1
			Cross_hist[0]->Fill(6.0);
			Cross_hist[1]->Fill(6.0,weight_);
		}else if(top_pass[0] && top_pass[1]){//P + Pip
			Cross_hist[0]->Fill(7.0);
			Cross_hist[1]->Fill(7.0,weight_);
		}else if(top_pass[0] && top_pass[2]){//P + Pip
			Cross_hist[0]->Fill(8.0);
			Cross_hist[1]->Fill(8.0,weight_);
		}else if(top_pass[1] && top_pass[2]){//P + Pip
			Cross_hist[0]->Fill(9.0);
			Cross_hist[1]->Fill(9.0,weight_);
		}else if(!top_pass[3] && (num_hadmiss == 3)){//P + Pip
			Cross_hist[0]->Fill(10.0);
			Cross_hist[1]->Fill(10.0,weight_);
		}else if(top_mult[0]){//P + Pip
			Cross_hist[0]->Fill(11.0);
			Cross_hist[1]->Fill(11.0,weight_);
		}else if(top_mult[1]){//P + Pip
			Cross_hist[0]->Fill(12.0);
			Cross_hist[1]->Fill(12.0,weight_);
		}else if(top_mult[2]){//P + Pip
			Cross_hist[0]->Fill(13.0);
			Cross_hist[1]->Fill(13.0,weight_);
		}else if(top_mult[3]){//P + Pip
			Cross_hist[0]->Fill(14.0);
			Cross_hist[1]->Fill(14.0,weight_);
		}
	}
}

void Histogram::Cross_Write(std::shared_ptr<Environment> envi_){
	//char dir_name[100];
	if(envi_->was_cross_plot()){
		std::cout<<"Cross Plots: ";
		TDirectory * Cross_plot = _RootOutputFile->mkdir("Cross Plots");
		Cross_plot->cd();
		Cross_hist[0]->SetXTitle("Event Topology Mixing");
		Cross_hist[0]->SetYTitle("Events");
		Cross_hist[0]->Write();
		Cross_hist[1]->SetXTitle("Event Topology Mixing");
		Cross_hist[1]->SetYTitle("Events");
		Cross_hist[1]->Write();
		std::cout<<"Done" <<std::endl;
	}
}

void Histogram::Thrown_Make(std::shared_ptr<Environment> _envi){
	//If we decide we want to fill this stuff
	if(_envi->was_Friend_plot() && _envi->was_sim()){
		char hname[100];
		sprintf(hname,"Thrown_2#pi_off_proton_#Delta^{++}");
		Double_t xmin1[7] = {_W_min,_Q2_min,_MM_min[0],_MM2_min[0],_theta_min,_alpha_min,_phi_min};
		Double_t xmax1[7] = {_W_max,_Q2_max,_MM_max[0],_MM2_max[0],_theta_max,_alpha_max,_phi_max};
		Double_t xmin2[7] = {_W_min,_Q2_min,_MM_min[1],_MM2_min[1],_theta_min,_alpha_min,_phi_min};
		Double_t xmax2[7] = {_W_max,_Q2_max,_MM_max[1],_MM2_max[1],_theta_max,_alpha_max,_phi_max};
		Double_t xmin3[7] = {_W_min,_Q2_min,_MM_min[2],_MM2_min[2],_theta_min,_alpha_min,_phi_min};
		Double_t xmax3[7] = {_W_max,_Q2_max,_MM_max[2],_MM2_max[2],_theta_max,_alpha_max,_phi_max};
		Thrown[0] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin1,xmax1);
		//Friend[0]->Sumw2();//Must be called before error bars can be placed. 
		sprintf(hname,"Thrown_2#pi_off_proton_#rho");
		Thrown[1] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin2,xmax2);
		//Friend[1]->Sumw2();//Must be called before error bars can be placed. 
		sprintf(hname,"Thrown_2#pi_off_proton_#Delta^{0}");
		Thrown[2] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin3,xmax3);
		//Friend[2]->Sumw2();//Must be called before error bars can be placed. 
	}
}

//Fill this with Thrown data after having already made the Friend
void Histogram::Thrown_Fill(std::shared_ptr<Environment> _envi, float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int chan_, float weight_){
	if(_envi->was_Friend_plot() && _envi->was_sim()){
		if(!std::isnan(W_) && !std::isnan(Q2_) && !std::isnan(MM_)&& !std::isnan(MM2_) && !std::isnan(theta_) && !std::isnan(alpha_) && !std::isnan(phi_) && !std::isnan(weight_)){
			if(Histogram::In_Friend_Bin(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,chan_)){
			//int *y = Friend_binning(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,chan_);
				Double_t x[7] = { W_, Q2_, MM_,MM2_, theta_, alpha_, phi_};
			//if(y[0]>=0){
				Thrown[chan_]->Fill(x,weight_);
			}
		}
	}
}

void Histogram::Thrown_Write(std::shared_ptr<Environment> _envi){
	//if(_envi->Environment::was_Friend_plot()){
		SparseFile->cd();
		for(int i = 0; i< 3; i++){
			Thrown[i]->Write();
		}
	//}
}


void Histogram::Acceptance_Make(std::shared_ptr<Environment> _envi){
	//If we decide we want to fill this stuff
	if(_envi->was_Friend_plot() && _envi->was_sim()){
		for(int i = 0; i< 5; i++){
			char hname[100];
			sprintf(hname,"2#pi_off_proton_#Delta^{++}_%s",topologies[i+1]);
			Double_t xmin1[7] = {_W_min,_Q2_min,_MM_min[0],_MM2_min[0],_theta_min,_alpha_min,_phi_min};
			Double_t xmax1[7] = {_W_max,_Q2_max,_MM_max[0],_MM2_max[0],_theta_max,_alpha_max,_phi_max};
			Double_t xmin2[7] = {_W_min,_Q2_min,_MM_min[1],_MM2_min[1],_theta_min,_alpha_min,_phi_min};
			Double_t xmax2[7] = {_W_max,_Q2_max,_MM_max[1],_MM2_max[1],_theta_max,_alpha_max,_phi_max};
			Double_t xmin3[7] = {_W_min,_Q2_min,_MM_min[2],_MM2_min[2],_theta_min,_alpha_min,_phi_min};
			Double_t xmax3[7] = {_W_max,_Q2_max,_MM_max[2],_MM2_max[2],_theta_max,_alpha_max,_phi_max};
			Acceptance[0][i] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin1,xmax1);
			//Friend[0]->Sumw2();//Must be called before error bars can be placed. 
			sprintf(hname,"2#pi_off_proton_#rho_%s",topologies[i+1]);
			Acceptance[1][i] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin2,xmax2);
			//Friend[1]->Sumw2();//Must be called before error bars can be placed. 
			sprintf(hname,"2#pi_off_proton_#Delta^{0}_%s",topologies[i+1]);
			Acceptance[2][i] = std::make_shared<THnSparseD>(hname,hname,7,_Friend_bins,xmin3,xmax3);
			//Friend[2]->Sumw2();//Must be called before error bars can be placed. 
		}
	}
}

//Fill this with Thrown data after having already made the Friend
void Histogram::Acceptance_Fill(std::shared_ptr<Environment> _envi, int top_, float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int chan_, float weight_){
	if(_envi->was_Friend_plot() && _envi->was_sim()){
		if(!std::isnan(W_) && !std::isnan(Q2_) && !std::isnan(MM_)&& !std::isnan(MM2_) && !std::isnan(theta_) && !std::isnan(alpha_) && !std::isnan(phi_) && !std::isnan(weight_)){
			int *y = Friend_binning(W_,Q2_,MM_,MM2_,theta_,alpha_,phi_,chan_);
			Double_t x[7] = { W_, Q2_, MM_,MM2_, theta_, alpha_, phi_};
			if(y[0]>=0){
				Acceptance[chan_][top_]->Fill(x,weight_);
			}
		}
	}
}


void Histogram::Charge_Make(std::shared_ptr<Environment> envi_){
	char hname[100];
	sprintf(hname,"Normalized Faraday Cup Yield");
	Charge_hist = std::make_shared<TH1F>(hname,hname, Charge_res, Charge_min, Charge_max);
}

void Histogram::Charge_Graph1(std::shared_ptr<Environment> envi_, std::vector<int> run_num_, std::vector<float> charge_, std::vector<float> run_size_){
	std::string file_name;
	  if(envi_->was_cluster()){
	    std::string path = fs::current_path().string();
	    //file_name = "/home/mclauchc/analysis/bin/analysis_runs/$name/$name.root";
	    file_name = "$boop/integrated_charge.root";
	    fun::replace(file_name,"$boop",path);
	    //replace(file_name, "$name", a_file_name);
	  }else{
	    file_name = "/Users/cmc/Desktop/analysis/analysis_clas6/bin/test/integrated_charge.root";
	    //file_name = "$name.root";
	    //replace(file_name, "$name", a_file_name);
	    //replace(file_name, "$name", a_file_name);
	  }
	OtherFile = std::make_shared<TFile>(file_name.c_str(),"RECREATE");
	//TFile OtherFile = TFile(file_name.c_str(),"RECREATE");
	def = new TCanvas("def1");
	OtherFile->cd();
	//OtherFile.cd();
	int size = run_num_.size();
	Float_t q_here[size];
	Float_t run_n[size]; 
	//Check to make sure they're the same size
	if(run_num_.size() == charge_.size() && charge_.size() == run_size_.size()){
		for(int i = 0; i<run_num_.size(); i++){
			q_here[i]=charge_[i]/run_size_[i];
			run_n[i]=run_num_[i];
			std::cout<<"\nRun: " <<run_n[i] << " Norm Q: " <<q_here[i];
		}

		TGraph* IntCharge = new TGraph(size,run_n,q_here);
		IntCharge->SetLineWidth(0);
		IntCharge->SetMarkerColor(2);
		IntCharge->SetMarkerStyle(21);
		IntCharge->GetXaxis()->SetTitle("Run Number");
		IntCharge->GetYaxis()->SetTitle("Normalized Charge (micro C/Event)");
		std::cout<<"Faraday Graph: ";
		TDirectory* Charge2_plot = OtherFile->mkdir("Normed Faraday Cup");
		//TDirectory* Charge2_plot = OtherFile.mkdir("Normed Faraday Cup");
		Charge2_plot->cd();
		//Charge2_plot.cd();
		IntCharge->SetTitle("Normalized Faraday Cup Charge per Run");
		//IntCharge->SetPointX("Run Number");
		//IntCharge->SetPointY("Faraday Cup Charge Normalized by Num Events");
		IntCharge->Write();
		//IntCharge.Write();
		std::cout<<"Done" <<std::endl;
	}else{
		std::cout<<"\nRun size: " <<run_num_.size() <<"\tCharge Size: " <<charge_.size();
	}
	
	
}
void Histogram::Charge_Graph2(std::shared_ptr<Environment> envi_, std::vector<std::vector<float>> run_num_, std::vector<std::vector<float>> charge_, std::vector<std::vector<float>> run_size_){
	std::string file_name;
	std::vector<float> run_num; 
	std::vector<float> charge; 
	std::vector<float> run_size;
	char hname[100];
	
	if(!envi_->was_sim()){
		sprintf(hname,"Normalized_Faraday_Cup_Yield");
		Find_Gold = std::make_shared<TH1F>(hname,hname, 80, 0.00002, 0.00006);


		  if(envi_->was_cluster()){
		    std::string path = fs::current_path().string();
		    //file_name = "/home/mclauchc/analysis/bin/analysis_runs/$name/$name.root";
		    file_name = "$boop/integrated_charge.root";
		    fun::replace(file_name,"$boop",path);
		    //replace(file_name, "$name", a_file_name);
		  }else{
		    file_name = "/Users/cmc/Desktop/analysis/analysis_clas6/bin/test/integrated_charge.root";
		    //file_name = "$name.root";
		    //replace(file_name, "$name", a_file_name);
		    //replace(file_name, "$name", a_file_name);
		  }
		OtherFile = std::make_shared<TFile>(file_name.c_str(),"RECREATE");
		//TFile OtherFile = TFile(file_name.c_str(),"RECREATE");
		def = new TCanvas("def1");
		OtherFile->cd();
		//OtherFile.cd();
		
		//Check to make sure they're the same size
		if(run_num_.size() == charge_.size() && charge_.size() == run_size_.size()){
			for(int j=0; j<run_num_.size(); j++){
				if(run_num_[j].size() == charge_[j].size() && charge_[j].size() == run_size_[j].size()){
					for(int k=0; k<run_num_[j].size(); k++){
						run_num.push_back(run_num_[j][k]);
						charge.push_back(charge_[j][k]);
						run_size.push_back(run_size_[j][k]);
					}
				}
			}
		}
		int size = run_num.size();
		Float_t q_here[size];
		Float_t run_n[size]; 
		if(run_num.size() == charge.size() && charge.size() == run_size.size()){
			for(int i = 0; i<run_num.size(); i++){
				q_here[i]=charge[i]/run_size[i];
				run_n[i]=run_num[i];
				//std::cout<<"\nRun: " <<run_n[i] << " Norm Q: " <<q_here[i];
				Find_Gold->Fill(q_here[i]);
			}

			TGraph* IntCharge = new TGraph(size,run_n,q_here);
			IntCharge->SetLineWidth(0);
			IntCharge->SetMarkerColor(2);
			IntCharge->SetMarkerStyle(21);
			IntCharge->GetXaxis()->SetTitle("Run Number");
			IntCharge->GetYaxis()->SetTitle("Normalized Charge (micro C/Event)");
			//std::cout<<"Faraday Graph: ";
			TDirectory* Charge2_plot = OtherFile->mkdir("Normed Faraday Cup");
			//TDirectory* Charge2_plot = OtherFile.mkdir("Normed Faraday Cup");
			Charge2_plot->cd();
			//Charge2_plot.cd();
			IntCharge->SetTitle("Normalized Faraday Cup Charge per Run");
			//IntCharge->SetPointX("Run Number");
			//IntCharge->SetPointY("Faraday Cup Charge Normalized by Num Events");
			Find_Gold->SetTitle("Normalized Faraday Cup Charge");
			IntCharge->Write();
			Find_Gold->Write();
			//IntCharge.Write();
			std::cout<<"Done" <<std::endl;
		}else{
			std::cout<<"\nRun size: " <<run_num_.size() <<"\tCharge Size: " <<charge_.size();
		}

		
		
	}
}

void Histogram::Charge_Fill(std::shared_ptr<Environment> envi_, int run_num_, float charge_, int events_){
	if(charge_ > 0){
		Charge_hist->Fill(run_num_,(float)events_/charge_);//treat the events/charge as a weight
	}else{
		Charge_hist->Fill(run_num_,0.0);//treat the events/charge as a weight
	}
	
}

void Histogram::Charge_Write(std::shared_ptr<Environment> envi_){
	std::cout<<"Faraday Plot: ";
	TDirectory* Charge_plot = _RootOutputFile->mkdir("Faraday Cup");
	Charge_plot->cd();
	Charge_hist->SetXTitle("Run Number");
	Charge_hist->SetYTitle("Faraday Cup Charge Normalized by Num Events");
	Charge_hist->Write();
	std::cout<<"Done" <<std::endl;
}


/*
void Histogram::WQ2_sf_Make(std::shared_ptr<Environment> envi_){
	char hname[100];
	sprintf(hname,"W_Q2_sf_elep<3.5GeV"); //constants.h and otherwise writing the specific cut to the right plot
	WQ2_hist_sf[0] = std::make_shared<TH2F>( hname, hname, WQxres, WQxmin, WQxmax, WQyres, WQymin, WQymax); // constants.h
	sprintf(hname,"W_Q2_sf_elep>3.5GeV"); //constants.h and otherwise writing the specific cut to the right plot
	WQ2_hist_sf[1] = std::make_shared<TH2F>( hname, hname, WQxres, WQxmin, WQxmax, WQyres, WQymin, WQymax); // constants.h
}
void Histogram::WQ2_sf_Fill(std::shared_ptr<Environment> envi_, float W_, float Q2_, float p_){
	if(p_ > 3.5){
		WQ2_hist_sf[1]->Fill(W_,Q2_);
	}else if(p_ < 3.5){
		WQ2_hist_sf[0]->Fill(W_,Q2_);
	}else{
		std::cout<<"what p is getting bad stuff? " <<p_ <<std::endl;
	}
}
void Histogram::WQ2_sf_Write(std::shared_ptr<Environment> envi_){
	char dir_name[100];
	TDirectory* dir_WQ2_sf = _RootOutputFile->mkdir("W vs. Q2 for simSF");
	dir_WQ2_sf->cd();
	for(int i = 0; i < 2; i++){
		WQ2_hist_sf[i]->SetXTitle("W (GeV)");
		WQ2_hist_sf[i]->SetYTitle("Q^{2} (GeV^{2}");
		WQ2_hist_sf[i]->SetOption("Colz");
		WQ2_hist_sf[i]->Write();
	}
}*/

/*
void Histogram::Event_Particle_Hist(std::shared_ptr<Environment> envi_, const Particle p1, float W_, int top_, int par_, bool pass_){
	std::cout<<"		Filling Particle Event" <<std::endl;
	int pass = -1; 
	if(pass_){
		pass = 0; 
	}else{
		pass = 1; 
	}
	//Electron
	if(pass != -1){
		if(par_ == 0 && _pid[0]){
			std::cout<<"			Electron cc_segm: " <<_cc_seg <<std::endl;
			Fid_Fill(envi_,top_+1,_theta,_phi,0,10,pass,W_,_p);
			SF_Fill(envi_,top_+1,_p,_etot,10,pass,W_,physics::get_sector(_phi));
			CC_Fill(envi_,top_+1,physics::get_sector(_phi),_cc_seg,_nphe,10,pass);
		}else if(par_ ==0){
			std::cout<<"			Electron issue, friend" <<std::endl;
		}
		if(par_ != 0 && _pid[par_]){
			Fid_Fill(envi_,top_+1,_theta,_phi,par_,6,pass,W_,_p);
			DT_Fill(envi_,top_+1,par_,_p,_dt[par_],6,pass,W_,physics::get_sector(_phi));
		}else if(par_ !=0){
			std::cout<<"			Hadron issue, friend" <<std::endl;
		}
	}
}
/*
void Histogram::Fill_EID(std::shared_ptr<Particle> par, float W_, float Q2_){
	bool san, fid, cc, sf;
	if(par->Particle::Is_sanity_pass(0)){
		Histogram::WQ2_Fill(0,1,W_,Q2_);
		Histogram::Fid_Fill(0,par->Particle::par_theta(),par->Particle::par_phi(),0,1,0,W_,par->Particle::par_p());
		Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),1,0,W_,sector[0]);
		Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),1,0);
		if(par->Particle::Is_fid_pass(0)){
			fid = true;
			Histogram::WQ2_Fill(0,2,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,2,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),2,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),2,0);
		}else{
			fid = false;
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,2,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),2,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),2,1);
		}
		if(par->Particle::Is_cc_pass()){
			cc = true;
			Histogram::WQ2_Fill(0,4,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,4,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),4,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),4,0);
		}else{
			cc = false;
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,4,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),4,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),4,1);
		}
		if(par->Is_sf_pass() && par->Is_min_cc_pass()){
			sf = true;
			Histogram::WQ2_Fill(0,3,W_,Q2_);
			Histogram::Fid_Fill(0,par->Particle::par_theta(),par->Particle::par_phi(),0,3,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),3,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),3,0);
		}else{
			sf = false;
			Histogram::Fid_Fill(0,par->Particle::par_theta(),par->Particle::par_phi(),0,3,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),3,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),3,1);
		}
		if(fid && sf){
			Histogram::WQ2_Fill(0,5,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,5,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),5,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),5,0);
		}else{
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,5,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),5,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),5,1);
		}
		if(fid && cc){
			Histogram::WQ2_Fill(0,6,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,6,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),6,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),6,0);
		}else{
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,6,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),6,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),6,1);
		}
		if(sf && cc){
			Histogram::WQ2_Fill(0,7,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,7,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),7,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),7,0);
		}else{
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,7,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),7,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),7,1);
		}
		if(sf && fid && cc){
			Histogram::WQ2_Fill(0,8,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,8,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),8,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),8,0);
			_elec = physics::Make_4Vector(par->Particle::par_p(),data->Branches::cx(0),data->Branches::cy(0),data->Branches::cz(0),me);
			good_electron++;
		}else{
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,8,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),8,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),8,1);
		}
		if(par->Particle::Bank()==ELECTRON){
			Histogram::WQ2_Fill(0,9,W_,Q2_);
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,9,0,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),9,0,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),9,0);
		}
		else{
			Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,9,1,W_,par->Particle::par_p());
			Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),9,1,W_,sector[0]);
			Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),9,1);
		}
	}else{
		san = false;
		Histogram::Fid_Fill(0,Particle::par_theta(),par->Particle::par_phi(),0,1,1,W_,par->Particle::par_p());
		Histogram::SF_Fill(0,par->Particle::par_p(),par->Particle::par_etot(),1,1,W_,sector[0]);
		Histogram::CC_Fill(0,par->Particle::par_cc_sect(),par->Particle::par_cc_segm(),par->Particle::par_nphe(),1,1);
	}
}

void Histogram::Fill_HID(std::shared_ptr<Particle> par){

}*/

