#include "histogram.hpp"

Histogram::Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_){
    Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,holes_,flags_);
    Histogram::Rad_Corr();
    //Histogram::Sparse_7to5(flags_);
    Histogram::Single_Diff(flags_);
}

Histogram::Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_){
    Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,flags_);
    Histogram::Rad_Corr();
    //Histogram::Sparse_7to5(flags_);
    Histogram::Single_Diff(flags_);
}


void Histogram::Extract_5d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_){
    std::cout<<"Extract 5d Histograms\n";
	char hname[100];
	sprintf(hname,"Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
	std::cout<<"Getting Thrown THnSparse " <<hname <<"\n";
    for(int i=0; i<_W_nbins_; i++){
        for(int j=0; j<_Q2_nbins_; j++){
            sprintf(hname,"Thrown_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _thrown_no_rad_5d[i][j] = (THnSparseD *)nr_sim_tree_->Get(hname);
            sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _sim_data_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _exp_data_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            _empty_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            _acceptance_5d[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
            _acceptance_5d[i][j]->Divide(_thrown_5d[i][j]);
            _N_5d[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
            _N_5d[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
	        _N_5d[i][j]->Divide(_acceptance_5d[i][j]);
            sprintf(hname,"Localized_Holes_50_W:%.3f-%.3f_Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _N_holes_5d[i][j] = (THnSparseD *)holes_->Get(hname);
	        _N_5d[i][j]->Add(_N_holes_5d[i][j]);
        }
    }
    for(int i = 0; i<_thrown_5d[0][0]->GetNdimensions(); i++){
		//_n_bins_7d.push_back(_thrown_7d->GetAxis(i)->GetNbins());
        _n_bins_5d.push_back(_thrown_5d[0][0]->GetAxis(i)->GetNbins());
	}
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

void Histogram::Extract_5d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_){
    std::cout<<"Extract 5d Histograms\n";
	char hname[100];
	sprintf(hname,"Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
	std::cout<<"Getting Thrown THnSparse " <<hname <<"\n";
    for(int i=0; i<_W_nbins_; i++){
        for(int j=0; j<_Q2_nbins_; j++){
            sprintf(hname,"Thrown_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _thrown_no_rad_5d[i][j] = (THnSparseD *)nr_sim_tree_->Get(hname);
            sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _sim_data_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _exp_data_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            _empty_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            _acceptance_5d[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
            _acceptance_5d[i][j]->Divide(_thrown_5d[i][j]);
            _N_5d[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
            _N_5d[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
	        _N_5d[i][j]->Divide(_acceptance_5d[i][j]);
        }
    }
    for(int i = 0; i<_thrown_5d[0][0]->GetNdimensions(); i++){
		//_n_bins_7d.push_back(_thrown_7d->GetAxis(i)->GetNbins());
        _n_bins_5d.push_back(_thrown_5d[0][0]->GetAxis(i)->GetNbins());
	}
    /*
	_thrown_7d = (THnSparseD *)sim_tree_->Get(hname);
	std::cout<<"thrown dimensionality: " <<_thrown_7d->GetNdimensions() <<"\n";
	for(int i = 0; i<_thrown_5d->GetNdimensions(); i++){
		//_n_bins_7d.push_back(_thrown_7d->GetAxis(i)->GetNbins());
        _n_bins_5d.push_back(_thrown_7d->GetAxis(i)->GetNbins());
	}
	//std::cout<<"printing thrown bin info\n";
	//Histogram::Print_Histogram_Bin_Info(_thrown_7d);
	Histogram::Extract_Bin_Info(flags_);
	std::cout<<"Getting Thrown THnSparse (no rad) " <<hname <<"\n";
	sprintf(hname,"Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
	_thrown_no_rad = (THnSparseD *)nr_sim_tree_->Get(hname);//Not sure if I need the raw yield for this
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
	/*
    /*
	std::cout<<"Adding Holes\n";
	sprintf(hname,"Localized_Holes_50");
	_N_holes = (THnSparseD *)holes_->Get(hname);
	_N->Add(_N_holes);
	*/

	//std::cout<<"First bin content of _N: " <<_N->GetBinContent(1) <<"\n";

}
/*
void Histogram::Extract_Bin_Info(Flags flags_){
	std::cout<<"Extract Bin Info\n";
	int DIM = -1;
	std::vector<int> n_bins_5d;
	double_1d bin_lowside_7d_1d;
	double_1d bin_topside_7d_1d;
	double_1d bin_lowside_5d_1d;
	double_1d bin_topside_5d_1d;
	double_1d bin_low_7d_1d;
	double_1d bin_up_7d_1d;
	double_1d bin_mid_7d_1d;
	double_1d bin_size_7d_1d;
	double_1d bin_edges_7d_1d;
	double_1d bin_low_5d_1d;
	double_1d bin_up_5d_1d;
	double_1d bin_mid_5d_1d;
	double_1d bin_size_5d_1d;
	double_1d bin_edges_5d_1d;
	//Get the 7dimensional bins ready
	for(int bin_fun=0; bin_fun<5; bin_fun++){
		_n_bins_5d.push_back(_thrown_7d->GetAxis(bin_fun+2)->GetNbins());
	}
	DIM = _thrown_7d->GetNdimensions();
	for(int i = 0; i<DIM; i++){
		bin_lowside_7d_1d.push_back(_thrown_7d->GetAxis(i)->GetBinLowEdge(1));
		bin_topside_7d_1d.push_back(_thrown_7d->GetAxis(i)->GetBinUpEdge(_n_bins_7d[i]));
		for(int j=0; j<_n_bins_7d[i];j++){
			bin_low_7d_1d.push_back(_thrown_7d->GetAxis(i)->GetBinLowEdge(j+1));
			bin_up_7d_1d.push_back(_thrown_7d->GetAxis(i)->GetBinUpEdge(j+1));
			bin_mid_7d_1d.push_back((bin_up_7d_1d[j]+bin_low_7d_1d[j])/2.0);
			bin_size_7d_1d.push_back(bin_up_7d_1d[j]-bin_low_7d_1d[j]);
			bin_edges_7d_1d.push_back(_thrown_7d->GetAxis(i)->GetBinLowEdge(j+1));
			if(j==_n_bins_7d[i]-1){
				bin_edges_7d_1d.push_back(_thrown_7d->GetAxis(i)->GetBinUpEdge(j+1));
			}
		}
		_bin_low_7d.push_back(bin_low_7d_1d);
		bin_low_7d_1d.clear();
		_bin_up_7d.push_back(bin_up_7d_1d);
		bin_up_7d_1d.clear();
		_bin_mid_7d.push_back(bin_mid_7d_1d);
		bin_mid_7d_1d.clear();
		_bin_size_7d.push_back(bin_size_7d_1d);
		bin_size_7d_1d.clear();
		_bin_edges_7d.push_back(bin_edges_7d_1d);
		bin_edges_7d_1d.clear();
	}
	//Get the 5d bins ready
	for(int j=0; j<DIM-2; j++){
		bin_lowside_5d_1d.push_back(_thrown_7d->GetAxis(j+2)->GetBinLowEdge(1));
		bin_topside_5d_1d.push_back(_thrown_7d->GetAxis(j+2)->GetBinUpEdge(_n_bins_5d[j]));
		for(int k=0; k<_n_bins_5d[j];k++){
			bin_low_5d_1d.push_back(_thrown_7d->GetAxis(j+2)->GetBinLowEdge(k+1));
			bin_up_5d_1d.push_back(_thrown_7d->GetAxis(j+2)->GetBinUpEdge(k+1));
			bin_mid_5d_1d.push_back((bin_up_5d_1d[k]+bin_low_5d_1d[k])/2.0);
			bin_size_5d_1d.push_back(bin_up_5d_1d[k]-bin_low_5d_1d[k]);
			bin_edges_5d_1d.push_back(_thrown_7d->GetAxis(j+2)->GetBinUpEdge(k+1));
			if(k==_n_bins_5d[j]-1){
				bin_edges_5d_1d.push_back(_thrown_7d->GetAxis(j+2)->GetBinUpEdge(k+1));
			}
		}
		_bin_low_5d.push_back(bin_low_5d_1d);
		bin_low_5d_1d.clear();
		_bin_up_5d.push_back(bin_up_5d_1d);
		bin_up_5d_1d.clear();
		_bin_mid_5d.push_back(bin_mid_5d_1d);
		bin_mid_5d_1d.clear();
		_bin_size_5d.push_back(bin_size_5d_1d);
		bin_size_5d_1d.clear();
		_bin_edges_5d.push_back(bin_edges_5d_1d);
		bin_edges_5d_1d.clear();
	}
}*/
/*
void Histogram::Sparse_7to5(Flags flags_){
	std::cout<<"Sparse_7to5\n";
	std::cout<<"\tFilling Those Sparse Friends\n";
	int bin_5d[5];
	int bin_7d[7] = {1,1,1,1,1,1,1};
	int out_n = 0;
	std::cout<<"N bins: " <<_N->GetNbins() <<"\n";
	char hname[100];

	Sparse_1d_star N_1d;
	int num_bins=5;
	int bins[5]={2,3,4,5,6};

	for(int i=0; i<_W_nbins_; i++){
		for(int j=0; j<_Q2_nbins_; j++){
			_N->GetAxis(0)->SetRange(i+1,i+1);
			_N->GetAxis(1)->SetRange(j+1,j+1);
			N_1d.push_back((THnSparseD*)_N->Projection(num_bins,bins,"E"));
			sprintf(hname,"N_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			N_1d[j]->SetNameTitle(hname,hname);
		}
		//std::cout<<"\tPushing back all of W:" <<i <<"\n";
		//std::cout<<"\t\tPushing back N\n";
		_N_5d.push_back(N_1d);
		N_1d.clear();
	}
}*/

void Histogram::Rad_Corr(){
	std::cout<<"Calculating Radiative Effects\n";
    _rad_corr_mu = 0.0;
    for(int i=0; i<_W_nbins_; i++){
		for(int j=0; j<_Q2_nbins_; j++){
            _rad_corr->Fill(Histogram::W_mid(i),Histogram::Q2_mid(j),fun::Sparse_Integral(_thrown_no_rad_5d[i][j])/fun::Sparse_Integral(_thrown_5d[i][j]));
            _rad_corr_mu+= fun::Sparse_Integral(_thrown_5d[i][j])/fun::Sparse_Integral(_thrown_no_rad_5d[i][j]);
        }
    }
	//_rad_corr = _thrown_7d_no_rad->Projection(1,0);//First make it the projections
	//_rad_corr->SetNameTitle("Rad_Corr","Rad_Corr");
	//_rad_corr->Divide(_thrown_7d->Projection(1,0));//Then divide by the according projection
	//std::cout<<"Thrown 7d Integral: " <<fun::Sparse_Integral(_thrown_7d) <<" Thrown nr 7d Integral: " <<fun::Sparse_Integral(_thrown_7d_no_rad) <<"\n";
	//_rad_corr->Scale(fun::nSparseIntegral(_thrown_7d)/fun::nSparseIntegral(_thrown_7d_no_rad));//Another place for localized hole filling?
	//_rad_corr_mu = fun::Sparse_Integral(_thrown_7d)/fun::Sparse_Integral(_thrown_7d_no_rad);
	_rad_corr->Scale(_rad_corr_mu);//Another place for localized hole filling?
	double_1d corr_tmp; 
	for(int i=0; i<_W_nbins_; i++){
		for(int j=0; j<_Q2_nbins_; j++){
			corr_tmp.push_back(_rad_corr->GetBinContent(i+1,j+1));
		}
		_rad_corr_array.push_back(corr_tmp);
		corr_tmp.clear();
	}
}

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
	dir_S->cd();
	TDirectory* dir_S1[5];
	TDirectory* dir_S2[5][_n_bins_7d[1]];
    TDirectory* dir_C1;
	char dirname[100];
    sprintf(dirname,"Integrated Check");
    dir_C1 = dir_S->mkdir(dirname);
	for(int k=0; k<5; k++){
		sprintf(dirname,"%s",_five_dim_[k]);
		dir_S1[k] = dir_S->mkdir(dirname);
		dir_S1[k]->cd();
		for(int j=0; j< _Q2_nbins_; j++){//Q2
			sprintf(dirname,"%s_Q2|%.2f-%.2f",_five_dim_[k],Histogram::Q2_low(j),Histogram::Q2_top(j));
			dir_S2[k][j] = dir_S1[k]->mkdir(dirname);
		}
	}
	double denom = 1.0;
	std::cout<<"\tDirectories Made\n\tWriting Histograms\n";
	for(int i=0; i<_W_nbins_; i++){//W
		for(int j=0; j< _Q2_nbins_; j++){//Q2
			for(int k=0; k<5; k++){
				dir_S2[k][j]->cd();
				//std::cout<<"\t\tWriting Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<k <<"\n";
				sprintf(hname,"%s_single_diff_W:%.3f-%.3f_Q2:%.2f-%.2f_top:%s_var:%s",_five_dim_[k],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				sprintf(xlabel,"%s %s",_five_dim_[k],_dim_units_[k]);
				sprintf(ylabel,"Diff CS (microbarns/(%s))",_dim_units_y_[k]);
				exp_ch_1d.push_back(_N_5d[i][j]->Projection(k,"E"));
				denom *= physics::Virtual_Photon_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]);
				denom *=flags_.Flags::L(0);
				if(flags_.Flags::Rad_Corr()){
					denom *= _rad_corr_array[i][j];
				}
				denom *=_bin_size_7d[0][i];
				denom *=_bin_size_7d[1][j];
				if(k>1){
					denom *=_bin_size_5d[k][0]*TMath::Pi()/180.0;//Divide by angle bins in radians
				}else{
					denom *=_bin_size_5d[k][0];//Phi in radians
				}
				
				exp_ch_1d[k]->Scale(1.0/denom);
			
				exp_ch_1d[k]->SetNameTitle(hname,hname);
				exp_ch_1d[k]->GetXaxis()->SetTitle(xlabel);
				exp_ch_1d[k]->GetYaxis()->SetTitle(ylabel);
				//std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<k <<"\n";
				exp_ch_1d[k]->Write();
				denom = 1.0;
			}
			//exp_ch_2d.push_back(exp_ch_1d);
			exp_ch_1d.clear();
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
	std::cout<<"\nCompleted Single Differential Cross Sections\n";
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