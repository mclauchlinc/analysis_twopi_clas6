#include "histogram.hpp"


Histogram::Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_){
    Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,flags_);
    //Histogram::Rad_Corr();
    //Histogram::Sparse_7to5(flags_);
    Histogram::Single_Diff(flags_);
}


void Histogram::Extract_5d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_){
    std::cout<<"Extract 5d Histograms\n";
	char hname[100];
    for(int i=0; i<_W_nbins_; i++){//W
        for(int j=0; j<_Q2_nbins_; j++){//Q2
			for(int k=0; k<3; k++){ //Var set
				sprintf(hname,"Thrown_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[k],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_thrown_5d[i][j][k] = (THnSparseD *)sim_tree_->Get(hname);
				_thrown_no_rad_5d[i][j][k] = (THnSparseD *)nr_sim_tree_->Get(hname);
				sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[k],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sim_data_5d[i][j][k] = (THnSparseD *)sim_tree_->Get(hname);
				_exp_data_5d[i][j][k] = (THnSparseD *)exp_tree_->Get(hname);
				_empty_5d[i][j][k] = (THnSparseD *)empty_tree_->Get(hname);
				/*
				_acceptance_5d[i][j][k] = (THnSparseD*)_sim_data_5d[i][j][k]->Clone();
				_acceptance_5d[i][j][k]->Divide(_thrown_5d[i][j][k]);
				sprintf(hname,"Acceptance_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[k],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_acceptance_5d[i][j][k]->SetNameTitle(hname,hname);
				sprintf(hname,"N_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[k],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_N_5d[i][j][k] = (THnSparseD*)_exp_data_5d[i][j][k]->Clone();
				_acceptance_5d[i][j][k]->SetNameTitle(hname,hname);
				_N_5d[i][j][k]->Add(_empty_5d[i][j][k],-flags_.Flags::Qr());//Empty target subtraction
				_N_5d[i][j][k]->Divide(_acceptance_5d[i][j][k]);
				*/
			}
        }
    }
    std::cout<<"\nHistograms Extracted\nGetting Dimensions\n";
    for(int i = 0; i<_thrown_5d[0][0][0]->GetNdimensions(); i++){
		//_n_bins_7d.push_back(_thrown_7d->GetAxis(i)->GetNbins());
        _n_bins_5d.push_back(_thrown_5d[0][0][0]->GetAxis(i)->GetNbins());
        std::cout<<"Dimension " <<i <<" " <<_n_bins_5d[i] <<"\n";
	}
    std::cout<<"Dimensions extracted\n";
    /*
    std::cout<<"Extract 7d Histograms\n";
	char hname[100];
	sprintf(hname,"Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
	std::cout<<"Getting Thrown THnSparse " <<hname <<"\n";
	_thrown_7d = (THnSparseD *)sim_tree_->Get(hname);
	std::cout<<"thrown dimensionality: " <<_thrown_7d->GetNdimensions() <<"\n";
	for(int i = 0; i<_thrown_7d->GetNdimensions(); i++){
		_n_bins_7d.push_back(_thrown_7d->GetAxis(i)->GetNbins());
	}
	//std::cout<<"printing thrown bin info\n";
	//Histogram::Print_Histogram_Bin_Info(_thrown_7d);
	Histogram::Extract_Bin_Info(flags_);
	std::cout<<"Getting Thrown THnSparse (no rad) " <<hname <<"\n";
	sprintf(hname,"Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
	_thrown_7d_no_rad = (THnSparseD *)nr_sim_tree_->Get(hname);//Not sure if I need the raw yield for this
	std::cout<<"\tdid it\n";
    sprintf(hname,"%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
    std::cout<<"Getting Exp THnSparse " <<hname <<"\n";
    _exp_data_7d = (THnSparseD *)exp_tree_->Get(hname);
    std::cout<<"Getting Exp Empty THnSparse " <<hname <<"\n";
    _empty_7d = (THnSparseD *)empty_tree_->Get(hname);
	std::cout<<"Getting Sim Recon THnSparse " <<hname <<"\n";
	sprintf(hname,"%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
	_sim_data_7d = (THnSparseD *)sim_tree_->Get(hname);
	std::cout<<"Calculating Acceptance\n";
	_acceptance_7d = (THnSparseD*)_sim_data_7d->Clone();
	_acceptance_7d->Divide(_thrown_7d);
	long num_are_one=0;
	std::cout<<"GetNbins() for acceptance_7d: " <<_acceptance_7d->GetNbins() <<"\n";
	for(long i=0; i<_acceptance_7d->GetNbins(); i++){
		if(_acceptance_7d->GetBinContent(i+1)==1.0){
			//_acceptance_7d->SetBinContent(i+1,0.0);
			num_are_one++;
		}
	}
	std::cout<<"num were one: " <<num_are_one <<"\n";
	_N=(THnSparseD*)_exp_data_7d->Clone();
	_N->Add(_empty_7d,-flags_.Flags::Qr());//Empty target subtraction
	_N->Divide(_acceptance_7d);//Acceptance Correction
	
	std::cout<<"Adding Holes\n";
	sprintf(hname,"Localized_Holes_50");
	_N_holes = (THnSparseD *)holes_->Get(hname);
	_N->Add(_N_holes);
	
	std::cout<<"First bin content of _N: " <<_N->GetBinContent(1) <<"\n";
    */
}

/*
void Histogram::Rad_Corr(){
	std::cout<<"Calculating Radiative Effects\n";
    double rad_corr_top[3] = {0.0,0.0,0.0};
    double rad_corr_bot[3] = {0.0,0.0,0.0};

	char hname[100];
    
 	for(int k=0; k<3; k++){
		sprintf(hname,"Radiative_Correction_%s",_var_set_[k]);
		_rad_corr[k] = new TH2D(hname,hname,_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
   		_rad_corr[k]->GetYaxis()->Set(5,_Q2_bins_);
		for(int i=0; i<_W_nbins_; i++){
			for(int j=0; j<_Q2_nbins_; j++){
				std::cout<<i <<" " <<j <<" " <<k <<"\nthrown integral:" <<fun::Sparse_Integral(_thrown_no_rad_5d[i][j][k]) <<" no_rad integral:" <<fun::Sparse_Integral(_thrown_5d[i][j][k]) <<"\n";
				_rad_corr[k]->Fill(Histogram::W_mid(i),Histogram::Q2_mid(j),fun::Sparse_Integral(_thrown_no_rad_5d[i][j][k])/fun::Sparse_Integral(_thrown_5d[i][j][k]));
				rad_corr_top[k]+= fun::Sparse_Integral(_thrown_5d[i][j][k]);
				rad_corr_bot[k]+= fun::Sparse_Integral(_thrown_no_rad_5d[i][j][k]);
			}
		}
		_rad_corr_mu[k] = rad_corr_top[k]/rad_corr_bot[k];
		_rad_corr[k]->Scale(_rad_corr_mu[k]);//Another place for localized hole filling?
		double_1d corr_tmp; 
		for(int i=0; i<_W_nbins_; i++){
			for(int j=0; j<_Q2_nbins_; j++){
				corr_tmp.push_back(_rad_corr[k]->GetBinContent(i+1,j+1));
				std::cout<<"\t" <<i <<" " <<j <<" rad_corr: " <<_rad_corr[k]->GetBinContent(i+1,j+1) <<"\n";
			}
			_rad_corr_array[k].push_back(corr_tmp);
			corr_tmp.clear();
		}
	}
}
*/

void Histogram::Single_Diff(Flags flags_){
    std::cout<<"Single Differential\n";
	if(!flags_.Flags::Plot_Single_Diff()){
		std::cout<<"\tNot Plotting Single Differential Cross Sections\n";
	 	return;
	}
	//TCanvas* def = new TCanvas("def1");
	//Convert the _exp_corr_5d and _exp_holes_5d to 3 dimensional sparse histograms for usage in plotting single differential cross sections
	char hname[100];
	char xlabel[100];
	char ylabel[100];
	int idx_2d[2]; 
	std::cout<<"Making Output File\n";
	_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	std::cout<<"Made File:" <<flags_.Flags::Output_File().c_str() <<"\n";
	_RootOutputFile->cd();
	//Get the 7dimensional bins ready
	TH1D_1d_star exp_ch_1d;
	TH1D_2d_star exp_ch_2d;
	std::cout<<"\tMaking Directory\n";
	TDirectory* dir_S = _RootOutputFile->mkdir("Single Differential CS");
    //std::cout<<"\t\tMade first directory\n";
	dir_S->cd();
    //_rad_corr->Write();
    //std::cout<<"\t\tWent int first directory\n";
	TDirectory* dir_S0[6];
	TDirectory* dir_S1[6][3];
	TDirectory* dir_S2[6][3][_Q2_nbins_];
    TDirectory* dir_C1;
	char dirname[100];
    sprintf(dirname,"Integrated Check");
    dir_C1 = dir_S->mkdir(dirname);
    //std::cout<<"\tMade some internal directories\n\t\tNow naming them\n";
	for(int l=0; l<6; l++){
		sprintf(dirname,"%s",_tree_types_[l]);
		dir_S0[l] = dir_S->mkdir(dirname);
		dir_S0[l]->cd();
		for(int k=0; k<3; k++){
			sprintf(dirname,"%s_%s",_tree_types_[l],_invar_mass_[k]);
			dir_S1[l][k] = dir_S0[l]->mkdir(dirname);
			for(int j=0; j< _Q2_nbins_; j++){//Q2
				sprintf(dirname,"%s_%s_Q2|%.2f-%.2f",_tree_types_[l],_invar_mass_[k],Histogram::Q2_low(j),Histogram::Q2_top(j));
				dir_S2[l][k][j] = dir_S1[l][k]->mkdir(dirname);
			}
		}
	}
    
	double denom = 1.0;
	int m=0;
	std::cout<<"\tDirectories Made\n\tWriting Histograms\n";
	for(int i=0; i<_W_nbins_; i++){//W
		for(int j=0; j< _Q2_nbins_; j++){//Q2
			for(int k=0; k<3; k++){//Which invariant mass
				if(k==0){
					m=2;
				}else{
					m=k-1;
				}
				for(int l=0; l<6; l++){
					dir_S2[l][k][j]->cd();
					//std::cout<<"\t\tWriting Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<k <<"\n";
					sprintf(hname,"%s_%s_diff_W:%.3f-%.3f_Q2:%.2f-%.2f_top:%s",_tree_types_[l],_invar_mass_[k],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),flags_.Flags::Top().c_str());
					sprintf(xlabel,"%s %s",_invar_mass_[k],"GeV");
					sprintf(ylabel,"Difference CS (microbarns/(%s))",_dim_units_y_[0]);
					//std::cout<<"\tpushing back projection\n";
					switch(l){
						case 0:
							exp_ch_1d.push_back(_exp_data_5d[i][j][k]->Projection(0,"E"));
							exp_ch_1d[l]->Add(_exp_data_5d[i][j][m]->Projection(1,"E"),-1.0);
						break;
						case 1:
							exp_ch_1d.push_back(_sim_data_5d[i][j][k]->Projection(0,"E"));
							exp_ch_1d[l]->Add(_sim_data_5d[i][j][m]->Projection(1,"E"),-1.0);
						break;
						case 2:
							exp_ch_1d.push_back(_empty_5d[i][j][k]->Projection(0,"E"));
							exp_ch_1d[l]->Add(_empty_5d[i][j][m]->Projection(1,"E"),-1.0);
						break;
						case 3:
							exp_ch_1d.push_back(_thrown_5d[i][j][k]->Projection(0,"E"));
							exp_ch_1d[l]->Add(_thrown_5d[i][j][m]->Projection(1,"E"),-1.0);
						break;
						case 4:
							exp_ch_1d.push_back(_thrown_no_rad_5d[i][j][k]->Projection(0,"E"));
							exp_ch_1d[l]->Add(_thrown_no_rad_5d[i][j][m]->Projection(1,"E"),-1.0);
						break;
						case 5:
							//exp_ch_1d.push_back(_N_5d[i][j][k]->Projection(0,"E"));
							//exp_ch_1d[l]->Add(_N_5d[i][j][m]->Projection(1,"E"),-1.0);
							exp_ch_1d.push_back(_thrown_no_rad_5d[i][j][k]->Projection(0,"E"));
							exp_ch_1d[l]->Add(_thrown_no_rad_5d[i][j][m]->Projection(1,"E"),-1.0);
						break;
					}
					//exp_ch_1d.push_back(_N_5d[i][j][k]->Projection(k,"E"));
					//exp_ch_1d[exp_ch_1d.size()-1].Add(_N_5d[i][j][k]->Projection(k,"E"))
					//std::cout<<"\tDenominator time\n\t\tvirtual photon flux\n";
					denom *= physics::Virtual_Photon_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]);
					//std::cout<<"Scaling Denominator post flux: " <<denom <<"\n";
					//std::cout<<"\t\tluminosity\n";
					denom *=flags_.Flags::L(0);
					//std::cout<<"Scaling Denominator post Luminosity: " <<denom <<"\n";
					//std::cout<<"\t\trad corr\n";
					//if(flags_.Flags::Rad_Corr()){
					//	denom *= _rad_corr_array[k][i][j];
					//	std::cout<<"Scaling Denominator post rad corr: " <<denom <<"\n";
					//}
					//std::cout<<"\t\tW bin\n";
					denom *=_W_res_;
					//std::cout<<"Scaling Denominator post W: " <<denom <<"\n";
					//std::cout<<"\t\tQ2 bin\n";
					denom *=(_Q2_bins_[j+1]-_Q2_bins_[j]);
					//std::cout<<"Scaling Denominator post Q2: " <<denom <<"\n";
					//if(k>1){
						//For Theta this will need to be undone and then modified for cosine theta
						//denom *=(_thrown_5d[i][j][k]->GetAxis(k)->GetBinUpEdge(2)-_thrown_5d[i][j][k]->GetAxis(k)->GetBinLowEdge(2))*TMath::Pi()/180.0;//Divide by angle bins in radians
						//denom *= Histogram::CosTheta()
						//std::cout<<"Scaling Denominator post Xij: " <<denom <<"\n";
					//}else{
						//std::cout<<"\t\txij bin\n";
						denom *=(_thrown_5d[i][j][k]->GetAxis(0)->GetBinUpEdge(2)-_thrown_5d[i][j][k]->GetAxis(0)->GetBinLowEdge(2));//Phi in radians
						//std::cout<<"Scaling Denominator post Xij: " <<denom <<"\n";
					//}
					//std::cout<<"\tScaling\n";
					//std::cout<<"Current Integral" <<exp_ch_1d[l]->Integral() <<"\n";
					//std::cout<<"Scaling Denominator: " <<denom <<"\n";
					//exp_ch_1d[k]->Scale(1.0/denom);
					//std::cout<<"Post Integral" <<exp_ch_1d[l]->Integral() <<"\n";

					exp_ch_1d[l]->SetNameTitle(hname,hname);
					exp_ch_1d[l]->GetXaxis()->SetTitle(xlabel);
					exp_ch_1d[l]->GetYaxis()->SetTitle(ylabel);
					//std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<k <<"\n";
					exp_ch_1d[l]->Write();
					denom = 1.0;
				}
				//exp_ch_2d.push_back(exp_ch_1d);
				exp_ch_1d.clear();
			}
		}
	}
    /*std::cout<<"moving on\n";
    TH1D* check_hist[5];
    //TH1D* check_hist2[5][_n_bins_7d[0]][_n_bins_7d[1]];
    dir_C1->cd();
    std::cout<<"\tCheck Directories Made\n\tWriting Histograms\n";
	for(int k=0; k<5; k++){
        sprintf(hname,"full_yield_%s_single_diff_top:%s_var:%s",_five_dim_[k],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
        sprintf(xlabel,"%s %s",_five_dim_[k],_dim_units_[k]);
        sprintf(ylabel,"Yield)",_dim_units_y_[k]);
        //std::cout<<"projection of full 7d\n";
        check_hist[k] = _exp_data_5d->Projection(k,"E");
        //check_hist[k] = _N->Projection(k,"E");
        check_hist[k]->SetNameTitle(hname,hname);
        check_hist[k]->GetXaxis()->SetTitle(xlabel);
        check_hist[k]->GetYaxis()->SetTitle(ylabel);
        //std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<k <<"\n";
        check_hist[k]->Write();
	}*/
	_RootOutputFile->Close();
	std::cout<<"\nCompleted Invariant mass differences\n";
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

double Histogram::MM_max(int W_bin_, int var_set_){
	return _W_min_+_W_res_*(W_bin_+1)-_MM_offset[var_set_];
}

double Histogram::MM2_max(int W_bin_, int var_set_){
	return _W_min_+_W_res_*(W_bin_+1)-_MM2_offset[var_set_];
}

double Histogram::CosTheta(int theta_bin_){
    double theta_res = ((double)_theta_max_ - (double)_theta_min_)/(double)_theta_bins_;
    return abs(TMath::Cos((_theta_min_+(theta_res*theta_bin_))*TMath::Pi()/180.0)-TMath::Cos((_theta_min_+(theta_res*theta_bin_+1))*TMath::Pi()/180.0));
}