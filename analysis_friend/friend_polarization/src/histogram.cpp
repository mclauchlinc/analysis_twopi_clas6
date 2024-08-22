#include "histogram.hpp"

Histogram::Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_){
    Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,holes_,flags_);
    Histogram::Rad_Corr();
    //Histogram::Sparse_7to5(flags_);
    Histogram::Polarization(flags_);
}

Histogram::Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_){
    Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,flags_);
    Histogram::Rad_Corr();
    //Histogram::Sparse_7to5(flags_);
    Histogram::Polarization(flags_);
}


void Histogram::Extract_5d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_){
    std::cout<<"Extract 5d Histograms\n";
	char hname[100];
    for(int i=0; i<_W_nbins_; i++){
        for(int j=0; j<_Q2_nbins_; j++){
            sprintf(hname,"Thrown_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _thrown_no_rad_5d[i][j] = (THnSparseD *)nr_sim_tree_->Get(hname);
            sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _sim_data_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _exp_data_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            _empty_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
            _acceptance_5d[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
            _acceptance_5d[i][j]->Divide(_thrown_5d[i][j]);
            sprintf(hname,"Acceptance_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _acceptance_5d[i][j]->SetNameTitle(hname,hname);
			for(long k=0; k<_acceptance_5d[i][j]->GetNbins(); k++){
				if(_acceptance_5d[i][j]->GetBinError(k)/_acceptance_5d[i][j]->GetBinContent(k) > _Acceptance_Rel_Error_Max_){
					_acceptance_5d[i][j]->SetBinContent(k,0.0);
				}
			}
			//std::cout<<"Making Yield\n";
            sprintf(hname,"N_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            //std::cout<<"Cloning Exp\n";
			_N_5d[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
            //std::cout<<"Renaming " <<i <<" " <<j <<"\n";
			_N_5d[i][j]->SetNameTitle(hname,hname);
            //std::cout<<"Empty Target Subtraction\n";
			_N_5d[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
            //std::cout<<"Divide by Acceptance\n";
			_N_5d[i][j]->Divide(_acceptance_5d[i][j]);
			//std::cout<<"Add Local Holes\n";
            sprintf(hname,"Localized_Holes_W:%.3f-%.3f_Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			//sprintf(hname,"Localized_Holes_50_W:%.3f-%.3f_Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _N_holes_5d[i][j] = (THnSparseD *)holes_->Get(hname);
	        _N_5d[i][j]->Add(_N_holes_5d[i][j]);
        }
    }
    std::cout<<"\nHistograms Extracted\nGetting Dimensions\n";
    for(int i = 0; i<_thrown_5d[0][0]->GetNdimensions(); i++){
		//_n_bins_7d.push_back(_thrown_7d->GetAxis(i)->GetNbins());
        _n_bins_5d.push_back(_thrown_5d[0][0]->GetAxis(i)->GetNbins());
        std::cout<<"Dimension " <<i <<" " <<_n_bins_5d[i] <<"\n";
	}
    std::cout<<"Dimensions extracted\n";
}

void Histogram::Extract_5d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_){
    std::cout<<"Extract 5d Histograms\n";
	char hname[100];
	sprintf(hname,"Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
	std::cout<<"Getting Thrown THnSparse " <<hname <<"\n";
    for(int i=0; i<_W_nbins_; i++){
        for(int j=0; j<_Q2_nbins_; j++){
            sprintf(hname,"Thrown_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            //std::cout<<"Getting Thrown THnSparses " <<hname <<"\n";
			_thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _thrown_no_rad_5d[i][j] = (THnSparseD *)nr_sim_tree_->Get(hname);
            sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            //std::cout<<"Getting Reconstructed THnSparses " <<hname <<"\n";
			_sim_data_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _exp_data_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            _empty_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
			//std::cout<<"Making Acceptance \n";
            _acceptance_5d[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
			//std::cout<<"\tDividing by Thrown\n";
            _acceptance_5d[i][j]->Divide(_thrown_5d[i][j]);
			sprintf(hname,"Acceptance_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _acceptance_5d[i][j]->SetNameTitle(hname,hname);
			for(long k=0; k<_acceptance_5d[i][j]->GetNbins(); k++){
				if(_acceptance_5d[i][j]->GetBinError(k)/_acceptance_5d[i][j]->GetBinContent(k) > _Acceptance_Rel_Error_Max_){
					_acceptance_5d[i][j]->SetBinContent(k,0.0);
				}
			}
			//std::cout<<"Cloning Exp\n";
            _N_5d[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
            //std::cout<<"Empty Target Sub\n" <<flags_.Flags::Qr() <<"\n";
			_N_5d[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
	        //std::cout<<"Divide by Acceptance\n";
			_N_5d[i][j]->Divide(_acceptance_5d[i][j]);
			/*if(flags_.Flags::Nonlocal_Holes()){
				std::cout<<"Making NonLocal Holes\n";
				_sim_holes_tmp_5d[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
				_sim_holes_tmp_5d[i][j]->Divide(_acceptance_5d[i][j]);
				_sim_holes_5d[i][j] = (THnSparseD*)_thrown_5d[i][j]->Clone();
				_sim_holes_5d[i][j]->Add(_sim_holes_tmp_5d[i][j],-1.0);
				for(long k=0; k<_sim_holes_5d[i][j]->GetNbins(); k++){
					_sim_holes_5d[i][j]->SetBinError(k,(_sim_holes_5d[i][j]->GetBinContent(k))/2.0);
				}
				_N_holes_5d[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_N_holes_5d[i][j]->Scale(fun::nSparseIntegral(_exp_data_5d[i][j])/fun::nSparseIntegral(_sim_data_5d[i][j]));
				_N_5d[i][j]->Add(_N_holes_5d[i][j],1.0);
			}*/
        }
    }
    for(int i = 0; i<_thrown_5d[0][0]->GetNdimensions(); i++){
		//_n_bins_7d.push_back(_thrown_7d->GetAxis(i)->GetNbins());
        _n_bins_5d.push_back(_thrown_5d[0][0]->GetAxis(i)->GetNbins());
	}
}

void Histogram::Rad_Corr(){
	std::cout<<"Calculating Radiative Effects\n";
    double rad_corr_top = 0.0;
    double rad_corr_bot = 0.0;
    _rad_corr = new TH2D("Radiative_Correction","Radiative_Correction",_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
    _rad_corr->GetYaxis()->Set(5,_Q2_bins_);
    
	for(int i=0; i<_W_nbins_; i++){
		for(int j=0; j<_Q2_nbins_; j++){
            //std::cout<<i <<" " <<j <<"\nthrown integral:" <<fun::Sparse_Integral(_thrown_no_rad_5d[i][j]) <<" no_rad integral:" <<fun::Sparse_Integral(_thrown_5d[i][j]) <<"\n";
            //_rad_corr->Fill(Histogram::W_mid(i),Histogram::Q2_mid(j),fun::Sparse_Integral(_thrown_no_rad_5d[i][j])/fun::Sparse_Integral(_thrown_5d[i][j]));
			_rad_corr->Fill(Histogram::W_mid(i),Histogram::Q2_mid(j),fun::Sparse_Integral(_thrown_5d[i][j])/fun::Sparse_Integral(_thrown_no_rad_5d[i][j]));
            //rad_corr_top+= fun::Sparse_Integral(_thrown_5d[i][j]);
            //rad_corr_bot+= fun::Sparse_Integral(_thrown_no_rad_5d[i][j]);
			rad_corr_top+= fun::Sparse_Integral(_thrown_no_rad_5d[i][j]);
            rad_corr_bot+= fun::Sparse_Integral(_thrown_5d[i][j]);
        }
    }
	//_rad_corr = _thrown_7d_no_rad->Projection(1,0);//First make it the projections
	//_rad_corr->SetNameTitle("Rad_Corr","Rad_Corr");
	//_rad_corr->Divide(_thrown_7d->Projection(1,0));//Then divide by the according projection
	//std::cout<<"Thrown 7d Integral: " <<fun::Sparse_Integral(_thrown_7d) <<" Thrown nr 7d Integral: " <<fun::Sparse_Integral(_thrown_7d_no_rad) <<"\n";
	//_rad_corr->Scale(fun::nSparseIntegral(_thrown_7d)/fun::nSparseIntegral(_thrown_7d_no_rad));//Another place for localized hole filling?
	//_rad_corr_mu = fun::Sparse_Integral(_thrown_7d)/fun::Sparse_Integral(_thrown_7d_no_rad);
	_rad_corr_mu = rad_corr_top/rad_corr_bot;
    _rad_corr->Scale(_rad_corr_mu);//Another place for localized hole filling?
	double_1d corr_tmp; 

	for(int i=0; i<_W_nbins_; i++){
		for(int j=0; j<_Q2_nbins_; j++){
			corr_tmp.push_back(_rad_corr->GetBinContent(i+1,j+1));
            std::cout<<"\t" <<i <<" " <<j <<" rad_corr: " <<_rad_corr->GetBinContent(i+1,j+1) <<"\n";
		}
		_rad_corr_array.push_back(corr_tmp);
		corr_tmp.clear();
	}
	
	/*_rad_corr_array = {	{ 0.775198,0.804581,0.811479,0.827047,0.829292 },
						{ 0.801247,0.821275,0.8326,0.844744,0.847486 },
						{ 0.821424,0.844058,0.85292,0.866613,0.866408 },
						{ 0.82837,0.848081,0.859096,0.868395,0.867181 },
						{ 0.861027,0.876247,0.890575,0.900434,0.896546 },
						{ 0.897684,0.911136,0.925244,0.93331,0.931846 },
						{ 0.911673,0.928591,0.948998,0.951515,0.947982 },
						{ 0.911908,0.92964,0.946975,0.95579,0.949468 },
						{ 0.927066,0.945327,0.961685,0.973025,0.972856 },
						{ 0.916798,0.936385,0.94859,0.955167,0.954566 },
						{ 0.923384,0.937913,0.952937,0.957886,0.95534 },
						{ 0.91669,0.924544,0.946102,0.94554,0.943438 },
						{ 0.919028,0.932221,0.956098,0.952331,0.950898 },
						{ 0.955865,0.968143,0.985961,0.991604,0.988503 },
						{ 0.9859,0.999213,1.02054,1.02307,1.01847 },
						{ 1.00696,1.01742,1.03747,1.04853,1.04314 },
						{ 1.01343,1.02234,1.04499,1.04835,1.04683 },
						{ 1.02198,1.04165,1.06371,1.05543,1.04857 },
						{ 1.03396,1.04435,1.06156,1.0582,1.05141 },
						{ 1.04005,1.05093,1.06718,1.06875,1.05725 },
						{ 1.04582,1.05348,1.07145,1.06954,1.05521 },
						{ 1.04437,1.05481,1.06588,1.07046,1.05948 },
						{ 1.05205,1.05808,1.08109,1.07789,1.06357 },
						{ 1.06119,1.07059,1.08582,1.08439,1.0787 },
						{ 1.07078,1.08555,1.10033,1.09523,1.0894 },
						{ 1.08771,1.09323,1.11688,1.10747,1.10258 },
						{ 1.0951,1.09919,1.11463,1.11389,1.10763 },
						{ 1.0888,1.10881,1.12096,1.11695,1.10986 },
						{ 1.09727,1.11428,1.12647,1.12352,1.11474 }};*/
}

void Histogram::Polarization(Flags flags_){
	std::cout<<"Polarization\n";
	//if(!flags_.Flags::Plot_Polarization()){
	//	std::cout<<"Not plotting polarization observables\n";
	//	return;
	//}
	char hname[100];
	char xlabel[100];
	char ylabel[100];
	//Get the 7dimensional bins ready
	TH1D_1d_star exp_ch_1d;
	std::cout<<"Making Output File\n";
	_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	std::cout<<"Made File:" <<flags_.Flags::Output_File().c_str() <<"\n";
	_RootOutputFile->cd();
	TDirectory* dir_P = _RootOutputFile->mkdir("Polarization Obs");
	TDirectory* dir_P1[4];
	TDirectory* dir_P2[4][_W_nbins_];
	TDirectory* dir_P3[4][_W_nbins_][_Q2_nbins_];
    TDirectory* dir_C1;
	char dirname[100];
	double denom = 1.0;
    //std::cout<<"\tMade some internal directories\n\t\tNow naming them\n";
	for(int k=0; k<4; k++){
		sprintf(dirname,"%s",_five_dim_[k]);
		dir_P1[k] = dir_P->mkdir(dirname);
		dir_P1[k]->cd();
		for(int j=0; j< _W_nbins_; j++){//W
			sprintf(dirname,"%s_W|%.3f-%.3f",_five_dim_[k],Histogram::W_low(j),Histogram::W_top(j));
			dir_P2[k][j] = dir_P1[k]->mkdir(dirname);
			for(int i=0; i<_Q2_nbins_; i++){
				sprintf(dirname,"%s_W|%.3f-%.3f_Q2|%.2f-%.2f",_five_dim_[k],Histogram::W_low(j),Histogram::W_top(j),Histogram::Q2_low(i),Histogram::Q2_top(i));
				dir_P3[k][j][i] = dir_P2[k][j]->mkdir(dirname);
			}
		}
	}
	for(int Wbin=0; Wbin<_W_nbins_; Wbin++){
		for(int Q2bin=0; Q2bin<_Q2_nbins_; Q2bin++){
			for(int Xij=0; Xij<4; Xij++){
				for(int Xijbin=0; Xijbin<_n_bins_5d[Xij]; Xijbin++){
					dir_P3[Xij][Wbin][Q2bin]->cd();
					sprintf(hname,"%s_2nd_order_diff_W:%.3f-%.3f_Q2:%.2f-%.2f_Xij:%.3f-%.3f_top:%s_var:%s",_five_dim_[Xij],Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin),_N_5d[Wbin][Q2bin]->GetAxis(Xij)->GetBinLowEdge(Xijbin),_N_5d[Wbin][Q2bin]->GetAxis(Xij)->GetBinUpEdge(Xijbin),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
					sprintf(xlabel,"%s %s",_five_dim_[4],_dim_units_[4]);
					sprintf(ylabel,"Diff CS (microbarns/(%s %s))",_dim_units_y_[Xij],_dim_units_y_[4]);
					if(Xij>1){
						denom *=(_N_5d[Wbin][Q2bin]->GetAxis(Xij)->GetBinUpEdge(Xijbin)-_N_5d[Wbin][Q2bin]->GetAxis(Xij)->GetBinLowEdge(Xijbin))*TMath::Pi()/180.0;//Divide by angle bins in radians
                	}else{
						denom *=(_N_5d[Wbin][Q2bin]->GetAxis(Xij)->GetBinUpEdge(Xijbin)-_N_5d[Wbin][Q2bin]->GetAxis(Xij)->GetBinLowEdge(Xijbin));//Phi in radians
					}
					denom *=(_N_5d[Wbin][Q2bin]->GetAxis(4)->GetBinUpEdge(2)-_N_5d[Wbin][Q2bin]->GetAxis(4)->GetBinLowEdge(2))*TMath::Pi()/180.0;
					_N_5d[Wbin][Q2bin]->GetAxis(Xij)->SetRange(Xijbin,Xijbin);
					exp_ch_1d.push_back(_N_5d[Wbin][Q2bin]->Projection(4,"E"));
					denom *= physics::Virtual_Photon_Flux((double)Histogram::W_mid(Wbin),(double)Histogram::Q2_mid(Q2bin),_beam_energy_[flags_.Flags::Run()]);
					denom *=flags_.Flags::L(flags_.Flags::Run());
					denom *= _rad_corr_array[Wbin][Q2bin];
					denom *=_W_res_;
			    	denom *=(_Q2_bins_[Q2bin+1]-_Q2_bins_[Q2bin]);
					exp_ch_1d[Xijbin]->Scale(1.0/denom);
					exp_ch_1d[Xijbin]->SetNameTitle(hname,hname);
					exp_ch_1d[Xijbin]->GetXaxis()->SetTitle(xlabel);
					exp_ch_1d[Xijbin]->GetYaxis()->SetTitle(ylabel);
					//std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<Xijbin <<"\n";
					exp_ch_1d[Xijbin]->Write();
					denom = 1.0;
					_N_5d[Wbin][Q2bin]->GetAxis(Xij)->SetRange();
				}
				exp_ch_1d.clear();
			}
		}
	}
	_RootOutputFile->Close();
	std::cout<<"\nCompleted Second Differential Cross Sections\n";
}

float Histogram::W_low(int i_){
    if(i_<0 || i_>=29){
        return NAN;
    }
    return _W_min_ + i_*_W_res_;
}
float Histogram::W_top(int i_){
    if(i_<0 || i_>=29){
        return NAN;
    }
    return _W_min_ + (i_+1)*_W_res_;
}
float Histogram::Q2_low(int i_){
    if(i_<0 || i_>=5){
        return NAN;
    }
    return _Q2_bins_[i_];
}
float Histogram::Q2_top(int i_){
    if(i_<0 || i_>=5){
        return NAN;
    }
    return _Q2_bins_[i_+1];
}

float Histogram::W_mid(int i_){
    if(i_<0 || i_>=29){
        return NAN;
    }
    return  (_W_min_ + (i_+0.5)*_W_res_);
}
float Histogram::Q2_mid(int i_){
    if(i_<0 || i_>=5){
        return NAN;
    }
    return (_Q2_bins_[i_] + _Q2_bins_[i_+1])/2.0;
}
/*
double Histogram::MM_max(int W_bin_, int var_set_){
	return _W_min_[var_set_]+_W_res_*(W_bin_+1)-_MM_offset[var_set_];
}

double Histogram::MM2_max(int W_bin_, int var_set_){
	return _W_min_[var_set_]+_W_res_*(W_bin_+1)-_MM2_offset[var_set_];
}*/

double Histogram::CosTheta(int theta_bin_){
    double theta_res = ((double)_theta_max_ - (double)_theta_min_)/(double)_theta_bins_;
    return abs(TMath::Cos((_theta_min_+(theta_res*theta_bin_))*TMath::Pi()/180.0)-TMath::Cos((_theta_min_+(theta_res*theta_bin_+1))*TMath::Pi()/180.0));
}

void Histogram::Make_Acceptance_Rel_Error(Flags flags_){
	char hname[100];
	std::cout<<"Making Relative Error Histograms\n";
	for(int a=0; a<_W_nbins_; a++){
		for(int b=0; b<_Q2_nbins_; b++){
			//std::cout<<"Starting with " <<a <<"," <<b <<"\n";
			//Int_t* coord = new Int_t[_acceptance_5d[i][j]->GetNbins()];
			_acceptance_rel_err_5d[a][b]=(THnSparseD*)_acceptance_5d[a][b]->Clone();
			sprintf(hname,"Acceptance Relative Error %s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(a),Histogram::W_top(a),Histogram::Q2_low(b),Histogram::Q2_top(b));
			_acceptance_rel_err_5d[a][b]->SetNameTitle(hname,hname);
			//std::cout<<"Making Rel Error Acceptance Hist for " <<a <<" " <<b <<"\n";
			//std::cout<<"  It is " <<_acceptance_5d[a][b]->GetNbins() <<" long\n";
			for(Long64_t i =0; i<_acceptance_5d[a][b]->GetNbins(); i++ ){
				//std::cout<<"\t" <<a <<"," <<b << " " <<i <<"/" <<_acceptance_5d[a][b]->GetNbins();
				if(_acceptance_5d[a][b]->GetBinContent(i)>0.0){
					//std::cout<<" acc val:" <<_acceptance_5d[a][b]->GetBinContent(i) <<" acc err:" <<_acceptance_5d[a][b]->GetBinError(i) <<"\n";
					_acceptance_rel_err_5d[a][b]->SetBinContent(i,_acceptance_5d[a][b]->GetBinError(i)/_acceptance_5d[a][b]->GetBinContent(i));
				}else{
					//std::cout<<" ** "  <<" acc val:" <<_acceptance_5d[a][b]->GetBinContent(i) <<" acc err:" <<_acceptance_5d[a][b]->GetBinError(i) <<"\n";
				}
			}
			//std::cout<<"Done with " <<a <<"," <<b <<"\n";
		}
	}
}


void Histogram::Acceptance_Rel_Error_Cut(){
	for(int a=0; a<_W_nbins_; a++){
		for(int b=0; b<_Q2_nbins_; b++){
			//Int_t* coord = new Int_t[_acceptance_5d[i][j]->GetNbins()];
			for(Long64_t i =0; i<_acceptance_5d[a][b]->GetNbins(); i++ ){
				if(_acceptance_rel_err_5d[a][b]->GetBinContent(i)>_Acceptance_Rel_Error_Max_){
					_acceptance_5d[a][b]->SetBinContent(i,0.0);
					_acceptance_5d[a][b]->SetBinError(i,0.0);
				}
			}
		}
	}
}