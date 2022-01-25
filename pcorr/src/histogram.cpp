#include "histogram.hpp"


Histogram::Histogram(std::shared_ptr<Flags> flags_){
	std::cout<<"Making Histograms\n";
	Histogram::ECorr_Angle_Make(flags_);
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
	Histogram::ECorr_Angle_Write(flags_);
	_RootOutputFile->Close();
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


int Histogram::Phi_Idx(float phi_){
	int idx=-1;
	float res = (_phi_max_-_phi_min_)/_phi_bins_;
	for(int i=0; i<_phi_bins_; i++){
		if(phi_ > (_phi_min_ + (float)i*res) || phi_<(_phi_min_ + ((float)i+1.0)*res)){
			idx = i;
		}
	}
	return idx;
}

float Histogram::Phi_Low(int bin_){
	if(bin_>=0 && bin_<_phi_bins_){
		float res = (_phi_max_-_phi_min_)/_phi_bins_;
		return _phi_min_ + (float)bin_*res;
	}else{
		return NAN;
	}
}
float Histogram::Phi_Top(int bin_){
	if(bin_>=0 && bin_<_phi_bins_){
		float res = (_phi_max_-_phi_min_)/_phi_bins_;
		return _phi_min_ + ((float)bin_+1.0)*res;
	}else{
		return NAN;
	}
}
float Histogram::Phi_Center(int bin_){
	if(bin_>=0 && bin_<_phi_bins_){
		float res = (_phi_max_-_phi_min_)/_phi_bins_;
		return _phi_min_ + ((float)bin_+0.5)*res;
	}else{
		return NAN;
	}
}

int Histogram::Theta_Idx(float theta_){
	int idx=-1;
	float res = (_theta_max_-_theta_min_)/_theta_bins_;
	for(int i=0; i<_theta_bins_; i++){
		if(theta_ > (_theta_min_ + (float)i*res) || theta_<(_theta_min_ + ((float)i+1.0)*res)){
			idx = i;
		}
	}
	return idx;
}
float Histogram::Theta_Low(int bin_){
	if(bin_>=0 && bin_<_theta_bins_){
		float res = (_theta_max_-_theta_min_)/_theta_bins_;
		return _theta_min_ + (float)bin_*res;
	}else{
		return NAN;
	}
}
float Histogram::Theta_Top(int bin_){
	if(bin_>=0 && bin_<_theta_bins_){
		float res = (_theta_max_-_theta_min_)/_theta_bins_;
		return _theta_min_ + ((float)bin_+1.0)*res;
	}else{
		return NAN;
	}
}

float Histogram::Theta_Center(int bin_){
	if(bin_>=0 && bin_<_theta_bins_){
		float res = (_theta_max_-_theta_min_)/_theta_bins_;
		return _theta_min_ + ((float)bin_+0.5)*res;
	}else{
		return NAN;
	}
}

//*------------------------------- Start Plot 1 Electron Angle Correction ---------------------------------*
void Histogram::ECorr_Angle_Make(std::shared_ptr<Flags> flags_){
	if(!flags_->Plot_E_PCorr()){ return;}
	std::cout<<"Making Electron Momentum Correction Delta Theta Plots\n";
	TH1F_ptr_3d plot_3d;
	TH1F_ptr_2d plot_2d;
	TH1F_ptr_1d plot_1d;


	std::vector<long> space_dims(4);
	space_dims[0] = _theta_bins_;//Theta bins
	space_dims[1] = _phi_bins_; //phi bins
	space_dims[2] = 6; //Sectors
	space_dims[3] = 3; //no corrections done, angle corrections done, full momentum corrections done
	CartesianGenerator cart(space_dims);
	char hname[100];
	int corr_idx = -1;
	int sec_idx = -1;
	int phi_idx = -1;
	int theta_idx = -1; 

	while(cart.GetNextCombination()){
		corr_idx=cart[3];
		sec_idx=cart[2];
		phi_idx=cart[1];
		theta_idx=cart[0];
		if((!flags_->Flags::E_Theta_Corr() && _ele_corr_[corr_idx]==_ele_angle_corr_) || (!flags_->Flags::E_PCorr() && _ele_corr_[corr_idx]==_ele_p_corr_)){
			//std::cout<<"\tFlag said not to make histograms for " <<corr_idx <<"\n";
		}else{
			sprintf(hname,"Delta_Theta_Sec:%s_Phi:%.2f-%.2f_Theta:%.2f-%.2f",_sector_[sec_idx],Phi_Low(phi_idx),Phi_Top(phi_idx),Theta_Low(theta_idx),Theta_Top(theta_idx));
			plot_1d.push_back(new TH1F(hname,hname,_delta_theta_res,_delta_theta_min,_delta_theta_max));
		}
		if(cart[0]==space_dims[0]-1){
			if(plot_1d.size()>0){
				plot_2d.push_back(plot_1d);
				plot_1d.clear();
			}
			if(cart[1]==space_dims[1]-1){
				if(plot_2d.size()>0){
					plot_3d.push_back(plot_2d);
					plot_2d.clear();
				}
				if(cart[2]==space_dims[2]-1){
					if(plot_3d.size()>0){
						_Delta_Theta_hist.push_back(plot_3d);
						plot_3d.clear();
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::ECorr_Angle_idx(float theta_, float phi_, int sector_, const char * corr_ , std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	if(flags_->Flags::E_Theta_Corr() && corr_ == _ele_angle_corr_){
		idx.push_back(1);
	}else if(flags_->Flags::E_PCorr() && corr_ == _ele_p_corr_){
		idx.push_back(2);
	}else if(corr_== _no_corr_){
		idx.push_back(0);
	}else{
		idx.push_back(-1);
	}
	idx.push_back(sector_-1);
	idx.push_back(Phi_Idx(phi_));
	idx.push_back(Theta_Idx(theta_));
	return idx;
}

void Histogram::ECorr_Angle_Fill(float delta_theta_, float theta_, float phi_, int sector_, const char* corr_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx = ECorr_Angle_idx(theta_,phi_,sector_,corr_,flags_);
	if(OK_Idx(idx)){
		_Delta_Theta_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(delta_theta_);
	}
}

void Histogram::ECorr_Angle_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_E_PCorr()){ return; }
	std::cout<<"Writing Electron Momentum Correction Plots\n";
	char dir_name[100];
	TDirectory* dir_e_corr = _RootOutputFile->mkdir("Electron PCorr");
	dir_e_corr->cd();
	TDirectory* dir_sub[6][_theta_bins_+1][_phi_bins_+1];//{sector,theta,phi}
	for(int i=0; i<6; i++){
		sprintf(dir_name,"Electron PCorr Sec:%s",i+1);
		dir_sub[i][0][0] = dir_e_corr->mkdir(dir_name);
		for(int j=0; j<_theta_bins_; j++){
			sprintf(dir_name,"Electron PCorr Sec:%s Theta:%.2f-%.2f",i+1,Theta_Low(j),Theta_Top(j));
			dir_sub[i][j+1][0] = dir_sub[i][0][0]->mkdir(dir_name);
			for(int k=0; k<_phi_bins_; k++){
				sprintf(dir_name,"Electron PCorr Sec:%s Theta:%.2f-%.2f Phi:%.2f-%.2f",i+1,Theta_Low(j),Theta_Top(j),Phi_Low(k),Phi_Top(k));
				dir_sub[i][j+1][k+1] = dir_sub[i][j+1][0]->mkdir(dir_name);
			} 
		}
	}
	std::vector<long> space_dims(4);
	space_dims[0] = _theta_bins_;//Theta bins
	space_dims[1] = _phi_bins_; //phi bins
	space_dims[2] = 6; //Sectors
	space_dims[3] = 3; //no corrections done, angle corrections done, full momentum corrections done
	CartesianGenerator cart(space_dims);
	char hname[100];
	int corr_idx = -1;
	int sec_idx = -1;
	int phi_idx = -1;
	int theta_idx = -1; 
	std::vector<int> idx;
	while(cart.GetNextCombination()){
		corr_idx=cart[3];
		sec_idx=cart[2];
		phi_idx=cart[1];
		theta_idx=cart[0];
		idx = ECorr_Angle_idx(Theta_Center(theta_idx),Phi_Center(phi_idx),sec_idx+1,_ele_corr_[corr_idx],flags_);
		if(OK_Idx(idx)){
			dir_sub[sec_idx][theta_idx+1][phi_idx+1]->cd();
			_Delta_Theta_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetXTitle("Delta Theta (deg)");
			_Delta_Theta_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetYTitle("Yield");
			_Delta_Theta_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Write();
		}
	}
}

//*------------------------------- End Plot 1 Electron Angle Correction ---------------------------------*
//*------------------------------- Start Elastic ---------------------------------*


//*------------------------------- End Elastic ---------------------------------*

