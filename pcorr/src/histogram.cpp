#include "histogram.hpp"


Histogram::Histogram(std::shared_ptr<Flags> flags_){
	std::cout<<"Making Histograms\n";
	Histogram::ECorr_Angle_Make(flags_);
	Histogram::Elastic_Make(flags_);
	Histogram::Check_Make(flags_);
	Histogram::Angular_Make(flags_);
	Histogram::E_PCorr_Make(flags_);
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
	Histogram::Elast_Write(flags_);
	Histogram::Check_Write(flags_);
	Histogram::Angular_Write(flags_);
	Histogram::E_PCorr_Write(flags_);
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
		if(phi_ > (_phi_min_ + (float)i*res) && phi_<(_phi_min_ + ((float)i+1.0)*res)){
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
		if(theta_ >= (_theta_min_ + (float)i*res) && theta_<(_theta_min_ + ((float)i+1.0)*res)){
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
			//std::cout<<"Making: " <<_Delta_Theta_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
			//std::cout<<"\tCorr idx: " <<corr_idx <<" sector idx:" <<sec_idx <<" phi_idx:" <<phi_idx <<" theta_idx:" <<theta_idx <<"\n";
			sprintf(hname,"Delta_Theta_%s_Sec:%s_Theta:%.2f-%.2f_Phi:%.2f-%.2f",_ele_corr_[corr_idx],_sector_[sec_idx],Theta_Low(theta_idx),Theta_Top(theta_idx),Phi_Low(phi_idx),Phi_Top(phi_idx));
			//std::cout<<"\tHistogram name: " <<hname <<"\n";
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
	}else if(flags_->Flags::E_Theta_Corr() && corr_ == _ele_angle_corr_){
		idx.push_back(1);
	}else if(flags_->Flags::E_PCorr() && corr_ == _ele_p_corr_){
		idx.push_back(2);
	}else{
		idx.push_back(-1);
	}
	idx.push_back(sector_-1);
	idx.push_back(Phi_Idx(phi_));
	idx.push_back(Theta_Idx(theta_));
	return idx;
}

void Histogram::ECorr_Angle_Fill(float delta_theta_, float theta_, float phi_, int sector_, const char* corr_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_E_PCorr()){ return; }
	std::vector<int> idx = ECorr_Angle_idx(theta_,phi_,sector_,corr_,flags_);
	if(OK_Idx(idx)){
		_Delta_Theta_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(delta_theta_);
		//std::cout<<"Tried to Fill |" <<delta_theta_ <<"| at | theta:" <<theta_ <<"__idx(" <<Theta_Idx(theta_) <<") " <<" phi:" <<phi_ <<"__idx(" <<Phi_Idx(phi_) <<") "<<" sector:" <<sector_ <<" corr_stat:" <<corr_ <<"\n";
		//std::cout<<"\tHist Index at: " <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<" " <<idx[3] <<"\n";
	}
}

void Histogram::ECorr_Angle_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_E_PCorr()){ return; }
	std::cout<<"Writing Electron Momentum Correction Plots\n";
	char dir_name[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_e_corr = _RootOutputFile->mkdir("Electron PCorr");
	//std::cout<<" Done\n";
	dir_e_corr->cd();
	//std::cout<<"Making Sub Directories: ";
	TDirectory* dir_sub[6][_theta_bins_+1][_phi_bins_+1];//{sector,theta,phi}
	//std::cout<<"Initialized\n";
	for(int i=0; i<6; i++){
		sprintf(dir_name,"Electron PCorr Sec:%d",i+1);
		//std::cout<<"\tMaking Sub Directory " <<i <<" " <<0 <<" " <<0 <<" : ";
		dir_sub[i][0][0] = dir_e_corr->mkdir(dir_name);
		//std::cout<<"Done\n";
		for(int j=0; j<_theta_bins_; j++){
			sprintf(dir_name,"Electron PCorr Sec:%d Theta:%.2f-%.2f",i+1,Theta_Low(j),Theta_Top(j));
			//std::cout<<"\t\tMaking Sub Directory " <<i <<" " <<j+1 <<" " <<0 <<" : ";
			dir_sub[i][j+1][0] = dir_sub[i][0][0]->mkdir(dir_name);
			//std::cout<<"Done\n";
			for(int k=0; k<_phi_bins_; k++){
				sprintf(dir_name,"Electron PCorr Sec:%d Theta:%.2f-%.2f Phi:%.2f-%.2f",i+1,Theta_Low(j),Theta_Top(j),Phi_Low(k),Phi_Top(k));
				//std::cout<<"\t\tMaking Final Sub Directory " <<i <<" " <<j+1 <<" " <<k+1 <<" : ";
				dir_sub[i][j+1][k+1] = dir_sub[i][j+1][0]->mkdir(dir_name);
				
				//std::cout<<"Done\n";
			} 
		}
	}
	std::cout<<"Directories Built\n";
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
			//std::cout<<"Trying to write at | theta:" <<Theta_Center(theta_idx) <<"idx(" <<theta_idx <<") " <<" phi:" <<Phi_Center(phi_idx) <<"idx(" <<phi_idx <<") "<<" sector:" <<sec_idx+1 <<" corr_stat:" <<_ele_corr_[corr_idx] <<"\n";
			//std::cout<<"\tGood Index at: " <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<" " <<idx[3] <<"\n";
			dir_sub[sec_idx][theta_idx+1][phi_idx+1]->cd();
			_Delta_Theta_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetXTitle("Delta Theta (deg)");
			_Delta_Theta_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetYTitle("Yield");
			_Delta_Theta_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Write();
		}
	}
}

//*------------------------------- End Plot 1 Electron Angle Correction ---------------------------------*
//*------------------------------- Start Elastic ---------------------------------*
void Histogram::Elastic_Make(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Elastic()){ return ;}

	TH1F_ptr_1d plot_1d;
	TH1F_ptr_2d plot_2d;


	std::vector<long> space_dims(3);
	space_dims[0] = 6; //Sectors
	space_dims[1] = 2; //Proton 35 deg cut
	space_dims[2] = 3; //Electron Corrections Performed
	CartesianGenerator cart(space_dims);
	char hname[100];

	while(cart.GetNextCombination()){
		if((!flags_->Flags::E_Theta_Corr() && _ele_corr_[cart[2]]==_ele_angle_corr_) || (!flags_->Flags::E_PCorr() && _ele_corr_[cart[2]]==_ele_p_corr_)){
			//
		}else{
			//std::cout<<"Trying to make: sector:" <<cart[0]+1 <<" ele_corr:" <<_ele_corr_[cart[2]] <<" pro_thres:" <<_proton_threshold_[cart[1]] <<"\n";
			//std::cout<<"Making Elastic Histogram: " <<_Elast_hist.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
			sprintf(hname,"Elastic_Peak_%s_Sector:%d_%s",_ele_corr_[cart[2]],cart[0]+1,_proton_threshold_[cart[1]]);
			plot_1d.push_back(new TH1F(hname,hname,_elast_res, _elast_min, _elast_max));
			if(cart[0] == space_dims[0]-1){
				if(plot_1d.size()>0){
					plot_2d.push_back(plot_1d);
					plot_1d.clear();
				}
				if(cart[1] == space_dims[1]-1){
					if(plot_2d.size()>0){
						_Elast_hist.push_back(plot_2d);
						plot_2d.clear();
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::Elastic_idx(int sector_, const char* corr_, const char* pro_thresh_, std::shared_ptr<Flags> flags_){
	
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
	if(pro_thresh_ == _no_pro_thresh_){
		idx.push_back(0);
	}else if(pro_thresh_ == _pro_thresh_){
		idx.push_back(1);
	}else{
		idx.push_back(-1);
	}
	idx.push_back(sector_-1);
	return idx;
}

void Histogram::Elastic_Fill(float W_, int sector_, const char* corr_, const char* pro_thresh_,std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Elastic()){ return ;}
	//std::cout<<"\tFilling Elastic Hist for " <<W_ <<" " <<sector_ <<" " <<corr_ <<" " <<pro_thresh_ <<"\n";
	std::vector<int> idx = Elastic_idx(sector_, corr_, pro_thresh_,flags_);
	if(OK_Idx(idx) && flags_->Flags::Plot_Elastic()){
		//std::cout<<"\t\tGood Index! Let's fill: " <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<"\n";
		_Elast_hist[idx[0]][idx[1]][idx[2]]->Fill(W_);
	}
}

void Histogram::Elast_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Elastic()){ return ;}
	std::cout<<"Writing Elastic Peak Plots\n";
	char dir_name[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_ela = _RootOutputFile->mkdir("Elastic Peak");
	//std::cout<<" Done\n";
	dir_ela->cd();
	//std::cout<<"Making Sub Directories: ";
	TDirectory* dir_sub[6][3];//{sector,proton_thesh}
	for(int i=0; i<6; i++){
		sprintf(dir_name,"Elastic Peak Sector:%d",i+1);
		dir_sub[i][0] = dir_ela->mkdir(dir_name);
		for(int j=0; j<2; j++){
			sprintf(dir_name,"Elastic Peak Sector:%d %s",i+1,_proton_threshold_[j]);
			dir_sub[i][j+1] = dir_sub[i][0]->mkdir(dir_name);
		}
	}
	std::vector<long> space_dims(3);
	space_dims[0] = 6; //Sectors
	space_dims[1] = 2; //Proton 35 deg cut
	space_dims[2] = 3; //Electron Corrections Performed
	CartesianGenerator cart(space_dims);
	char hname[100];
	std::vector<int> idx;

	while(cart.GetNextCombination()){
		if((!flags_->Flags::E_Theta_Corr() && _ele_corr_[cart[2]]==_ele_angle_corr_) || (!flags_->Flags::E_PCorr() && _ele_corr_[cart[2]]==_ele_p_corr_)){
			//
		}else{
			idx = Elastic_idx(cart[0]+1,_ele_corr_[cart[2]],_proton_threshold_[cart[1]],flags_);
			if(OK_Idx(idx)){
				dir_sub[cart[0]][cart[1]+1]->cd();
				_Elast_hist[idx[0]][idx[1]][idx[2]]->SetXTitle("W (GeV)");
				_Elast_hist[idx[0]][idx[1]][idx[2]]->SetYTitle("Yield");
				_Elast_hist[idx[0]][idx[1]][idx[2]]->Write();
			}
		}
	}
}

//*------------------------------- End Elastic ---------------------------------*
//*------------------------------- Start Check ---------------------------------*
void Histogram::Check_Make(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Check()){ return ;}

	TH2F_ptr_1d plot_1d;


	std::vector<long> space_dims(2);
	space_dims[0] = 6; //Sectors
	space_dims[1] = 4; //Types of Check Plots
	CartesianGenerator cart(space_dims);
	char hname[100];

	while(cart.GetNextCombination()){
		//std::cout<<"Trying to make: sector:" <<cart[0]+1 <<" ele_corr:" <<_ele_corr_[cart[2]] <<" pro_thres:" <<_proton_threshold_[cart[1]] <<"\n";
		//std::cout<<"Making Check Histogram: " <<_Check_hist.size() <<" " <<plot_1d.size() <<"\n";
		sprintf(hname,"%s_%d_Sector:%d",_check_names_[cart[1]],cart[0]+1);
		plot_1d.push_back(new TH2F(hname,hname,_check_xres[cart[1]], _check_xmin[cart[1]], _check_xmax[cart[1]],_check_yres[cart[1]], _check_ymin[cart[1]], _check_ymax[cart[1]]));
		if(cart[0] == space_dims[0]-1){
			if(plot_1d.size()>0){
				_Check_hist.push_back(plot_1d);
				plot_1d.clear();
			}
		}
	}
}

std::vector<int> Histogram::Check_idx(int sector_, const char* check_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	//std::cout<<"Trying to fill in " <<sector_ <<" " <<check_ <<"\n";
	if(!flags_->Plot_Check()){
		idx.push_back(-1);
		idx.push_back(-1);
		return idx;
	}
	for(int i=0; i<4; i++){
		if(check_ == _check_names_[i]){
			idx.push_back(i);
		}
	}
	if(idx.size()==0){
		idx.push_back(-1);
	}
	idx.push_back(sector_-1);
	//std::cout<<"idx: " <<idx[0] <<" " <<idx[1] <<"\n";
	return idx;
}

void Histogram::Check_Fill(float xval_, float yval_, int sector_, const char* check_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Check()){ return ;}
	//std::cout<<"\tFilling Check Hist for " <<xval_ <<" " <<yval_ <<" " <<sector_ <<" " <<check_ <<"\n";
	std::vector<int> idx = Check_idx(sector_, check_,flags_);
	if(OK_Idx(idx) && flags_->Flags::Plot_Check()){
		//std::cout<<"\t\tGood Index! Let's fill: " <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<"\n";
		_Check_hist[idx[0]][idx[1]]->Fill(xval_,yval_);
	}
}

void Histogram::Check_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Check()){ return ;}
	std::cout<<"Writing Check Plots\n";
	char dir_name[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_check = _RootOutputFile->mkdir("Check Plots");
	//std::cout<<" Done\n";
	dir_check->cd();
	//std::cout<<"Making Sub Directories: ";
	TDirectory* dir_sub[6][5];//{sector,proton_thesh}
	for(int i=0; i<6; i++){
		sprintf(dir_name,"Check Plots Sector:%d",i+1);
		dir_sub[i][0] = dir_check->mkdir(dir_name);
		for(int j=0; j<4; j++){
			sprintf(dir_name,"Check Plots Sector:%d %s",i+1,_check_names_[j]);
			dir_sub[i][j+1] = dir_sub[i][0]->mkdir(dir_name);
		}
	}
	std::vector<long> space_dims(2);
	space_dims[0] = 6; //Sectors
	space_dims[1] = 4; //Check Plots
	CartesianGenerator cart(space_dims);
	char hname[100];
	std::vector<int> idx;

	while(cart.GetNextCombination()){
		idx = Check_idx(cart[0]+1,_check_names_[cart[1]],flags_);
		if(OK_Idx(idx)){
			dir_sub[cart[0]][cart[1]+1]->cd();
			_Check_hist[idx[0]][idx[1]]->SetXTitle(_check_xnames[cart[1]]);
			_Check_hist[idx[0]][idx[1]]->SetYTitle(_check_ynames[cart[1]]);
			_Check_hist[idx[0]][idx[1]]->Write();
		}
	}
}

//*------------------------------- End Check ---------------------------------*
//*------------------------------- Start Angular ---------------------------------*
void Histogram::Angular_Make(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Fid(0)){ return ;}

	TH2F_ptr_1d plot_1d;
	TH2F_ptr_2d plot_2d;
	TH2F_ptr_3d plot_3d;


	std::vector<long> space_dims(4);
	space_dims[0] = 6; //Sectors
	space_dims[1] = 2; //Proton 35 deg cut
	space_dims[2] = 3; //Electron Corrections Performed
	space_dims[3] = 2; //Fiducial Cut
	CartesianGenerator cart(space_dims);
	char hname[100];

	while(cart.GetNextCombination()){
		if(!flags_->Flags::Plot_Fid(0)){
			//
		}else if(cart[3]!=1 || flags_->Flags::Fid_Cut(0)){
			//std::cout<<"Trying to make: sector:" <<cart[0]+1 <<" ele_corr:" <<_ele_corr_[cart[2]] <<" pro_thres:" <<_proton_threshold_[cart[1]] <<"\n";
			//std::cout<<"Making Angular Histogram: " <<_Elast_hist.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
			sprintf(hname,"Angular_Distribution_%s_Sector:%d_%s_%s",_ele_corr_[cart[2]],cart[0]+1,_proton_threshold_[cart[1]],_pcorr_ele_fid_[cart[3]]);
			plot_1d.push_back(new TH2F(hname,hname,_fid_xbin, _fid_xmin, _fid_xmax,_fid_ybin, _fid_ymin, _fid_ymax));
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
							_Angular_hist.push_back(plot_3d);
							plot_3d.clear();
						}
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::Angular_idx(int sector_, const char* corr_, const char* pro_thresh_, const char* fid_cut_, std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	//std::cout<<"Filling Angular: " <<sector_ <<" " <<corr_ <<" " <<pro_thresh_ <<" " <<fid_cut_ <<"\n";
	if(flags_->Fid_Cut(0) && fid_cut_ == _fid_cut_){
		idx.push_back(1);
	}else if(fid_cut_ == _no_cut_){
		idx.push_back(0);
	}else{
		idx.push_back(-1);
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
	if(pro_thresh_ == _no_pro_thresh_){
		idx.push_back(0);
	}else if(pro_thresh_ == _pro_thresh_){
		idx.push_back(1);
	}else{
		idx.push_back(-1);
	}
	idx.push_back(sector_-1);
	//fun::print_vector_idx(idx);
	return idx;
}

void Histogram::Angular_Fill(float theta_, float phi_, int sector_, const char* corr_, const char* pro_thresh_, const char* fid_cut_,std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Fid(0)){ return ;}
	//std::cout<<"\tFilling Angular Hist for " <<W_ <<" " <<sector_ <<" " <<corr_ <<" " <<pro_thresh_ <<"\n";
	std::vector<int> idx = Angular_idx(sector_, corr_, pro_thresh_,fid_cut_,flags_);
	if(OK_Idx(idx) && flags_->Flags::Plot_Fid(0)){
		//std::cout<<"\t\tGood Index! Let's fill: " <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<"\n";
		_Angular_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(phi_,theta_);
	}
}

void Histogram::Angular_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Fid(0)){ return ;}
	std::cout<<"Writing Angular Peak Plots\n";
	char dir_name[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_ang = _RootOutputFile->mkdir("Angular Plots");
	//std::cout<<" Done\n";
	dir_ang->cd();
	//std::cout<<"Making Sub Directories: ";
	TDirectory* dir_sub[6][3+1];//{sector,proton_thesh}
	for(int i=0; i<6; i++){
		sprintf(dir_name,"Angular Plots Sector:%d",i+1);
		dir_sub[i][0] = dir_ang->mkdir(dir_name);
		for(int j=0; j<3; j++){
			sprintf(dir_name,"Angular Plots Sector:%d %s",i+1,_ele_corr_[j]);
			dir_sub[i][j+1] = dir_sub[i][0]->mkdir(dir_name);
		}
	}
	std::vector<long> space_dims(4);
	space_dims[0] = 6; //Sectors
	space_dims[1] = 2; //Proton 35 deg cut
	space_dims[2] = 3; //Electron Corrections Performed
	space_dims[3] = 2; //Fiducial Cut
	CartesianGenerator cart(space_dims);
	char hname[100];
	std::vector<int> idx;

	while(cart.GetNextCombination()){
		if(!flags_->Flags::Plot_Fid(0)){
			//
		}else{
			idx = Angular_idx(cart[0]+1,_ele_corr_[cart[2]],_proton_threshold_[cart[1]],_pcorr_ele_fid_[cart[3]],flags_);
			if(OK_Idx(idx)){
				dir_sub[cart[0]][cart[2]+1]->cd();
				_Angular_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetXTitle("Phi (deg)");
				_Angular_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetYTitle("Theta (deg)");
				_Angular_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Write();
			}
		}
	}
}

//*------------------------------- End Elastic ---------------------------------*
//*------------------------------- Start Delta P ---------------------------------*
void Histogram::E_PCorr_Make(std::shared_ptr<Flags> flags_){
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
			//std::cout<<"Making: " <<_Delta_Theta_hist.size() <<" " <<plot_3d.size() <<" " <<plot_2d.size() <<" " <<plot_1d.size() <<"\n";
			//std::cout<<"\tCorr idx: " <<corr_idx <<" sector idx:" <<sec_idx <<" phi_idx:" <<phi_idx <<" theta_idx:" <<theta_idx <<"\n";
			sprintf(hname,"Delta_P_ele_%s_Sec:%s_Theta:%.2f-%.2f_Phi:%.2f-%.2f",_ele_corr_[corr_idx],_sector_[sec_idx],Theta_Low(theta_idx),Theta_Top(theta_idx),Phi_Low(phi_idx),Phi_Top(phi_idx));
			//std::cout<<"\tHistogram name: " <<hname <<"\n";
			plot_1d.push_back(new TH1F(hname,hname,_delta_p_res,_delta_p_min,_delta_p_max));
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
						_E_PCorr_hist.push_back(plot_3d);
						plot_3d.clear();
					}
				}
			}
		}
	}
}

std::vector<int> Histogram::E_PCorr_idx(float theta_, float phi_, int sector_, const char * corr_ , std::shared_ptr<Flags> flags_){
	std::vector<int> idx;
	if(flags_->Flags::E_Theta_Corr() && corr_ == _ele_angle_corr_){
		idx.push_back(1);
	}else if(flags_->Flags::E_PCorr() && corr_ == _ele_p_corr_){
		idx.push_back(2);
	}else if(corr_== _no_corr_){
		idx.push_back(0);
	}else if(flags_->Flags::E_Theta_Corr() && corr_ == _ele_angle_corr_){
		idx.push_back(1);
	}else if(flags_->Flags::E_PCorr() && corr_ == _ele_p_corr_){
		idx.push_back(2);
	}else{
		idx.push_back(-1);
	}
	idx.push_back(sector_-1);
	idx.push_back(Phi_Idx(phi_));
	idx.push_back(Theta_Idx(theta_));
	return idx;
}

void Histogram::E_PCorr_Fill(float delta_p_e_, float theta_, float phi_, int sector_, const char* corr_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_E_PCorr()){ return; }
	//std::cout<<"Filling E_PCorr Hist with " <<delta_p_e_ <<" theta: " <<theta_ <<" phi:" <<phi_ <<" sector:" <<sector_ <<" corr:" <<corr_ <<"\n";
	std::vector<int> idx = E_PCorr_idx(theta_,phi_,sector_,corr_,flags_);
	if(OK_Idx(idx)){
		_E_PCorr_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Fill(delta_p_e_);
		//std::cout<<"Tried to Fill |" <<delta_theta_ <<"| at | theta:" <<theta_ <<"__idx(" <<Theta_Idx(theta_) <<") " <<" phi:" <<phi_ <<"__idx(" <<Phi_Idx(phi_) <<") "<<" sector:" <<sector_ <<" corr_stat:" <<corr_ <<"\n";
		//std::cout<<"\tHist Index at: " <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<" " <<idx[3] <<"\n";
	}
}

void Histogram::E_PCorr_Write(std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_E_PCorr()){ return; }
	std::cout<<"Writing Electron Momentum Correction Plots\n";
	char dir_name[100];
	//std::cout<<"Making Large Directory:";
	TDirectory* dir_e_corr = _RootOutputFile->mkdir("Electron PCorr2");
	//std::cout<<" Done\n";
	dir_e_corr->cd();
	//std::cout<<"Making Sub Directories: ";
	TDirectory* dir_sub[6][_theta_bins_+1][_phi_bins_+1];//{sector,theta,phi}
	//std::cout<<"Initialized\n";
	for(int i=0; i<6; i++){
		sprintf(dir_name,"Electron PCorr2 Sec:%d",i+1);
		//std::cout<<"\tMaking Sub Directory " <<i <<" " <<0 <<" " <<0 <<" : ";
		dir_sub[i][0][0] = dir_e_corr->mkdir(dir_name);
		//std::cout<<"Done\n";
		for(int j=0; j<_theta_bins_; j++){
			sprintf(dir_name,"Electron PCorr2 Sec:%d Theta:%.2f-%.2f",i+1,Theta_Low(j),Theta_Top(j));
			//std::cout<<"\t\tMaking Sub Directory " <<i <<" " <<j+1 <<" " <<0 <<" : ";
			dir_sub[i][j+1][0] = dir_sub[i][0][0]->mkdir(dir_name);
			//std::cout<<"Done\n";
			for(int k=0; k<_phi_bins_; k++){
				sprintf(dir_name,"Electron PCorr2 Sec:%d Theta:%.2f-%.2f Phi:%.2f-%.2f",i+1,Theta_Low(j),Theta_Top(j),Phi_Low(k),Phi_Top(k));
				//std::cout<<"\t\tMaking Final Sub Directory " <<i <<" " <<j+1 <<" " <<k+1 <<" : ";
				dir_sub[i][j+1][k+1] = dir_sub[i][j+1][0]->mkdir(dir_name);
				
				//std::cout<<"Done\n";
			} 
		}
	}
	std::cout<<"Directories Built\n";
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
		idx = E_PCorr_idx(Theta_Center(theta_idx),Phi_Center(phi_idx),sec_idx+1,_ele_corr_[corr_idx],flags_);
		if(OK_Idx(idx)){
			//std::cout<<"Trying to write at | theta:" <<Theta_Center(theta_idx) <<"idx(" <<theta_idx <<") " <<" phi:" <<Phi_Center(phi_idx) <<"idx(" <<phi_idx <<") "<<" sector:" <<sec_idx+1 <<" corr_stat:" <<_ele_corr_[corr_idx] <<"\n";
			//std::cout<<"\tGood Index at: " <<idx[0] <<" " <<idx[1] <<" " <<idx[2] <<" " <<idx[3] <<"\n";
			dir_sub[sec_idx][theta_idx+1][phi_idx+1]->cd();
			_E_PCorr_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetXTitle("Delta P (GeV)");
			_E_PCorr_hist[idx[0]][idx[1]][idx[2]][idx[3]]->SetYTitle("Yield");
			_E_PCorr_hist[idx[0]][idx[1]][idx[2]][idx[3]]->Write();
		}
	}
}
//*------------------------------- End Delta P ---------------------------------*


































