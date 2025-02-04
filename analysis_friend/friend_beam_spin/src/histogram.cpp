#include "histogram.hpp"

Histogram::Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_){
    Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,holes_,flags_);
    Histogram::Rad_Corr();
    //Histogram::Sparse_7to5(flags_);
    //Histogram::Single_Diff(flags_);
	Histogram::Beam_Spin(flags_);
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
    for(int i=0; i<_W_nbins_; i++){
        for(int j=0; j<_Q2_nbins_; j++){
            sprintf(hname,"Thrown_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            //std::cout<<"\tExtracting Thrown Histograms: " <<hname <<"\n";
			_thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
           // _thrown_no_rad_5d[i][j] = (THnSparseD *)nr_sim_tree_->Get(hname);
            sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            //std::cout<<"\tExtracting Reconstructed Histograms: " <<hname <<"\n";
			_sim_data_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
           // _exp_data_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
           // _empty_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
			_empty_pos_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
			_empty_neg_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
			//sprintf(hname,"%s_%s_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
			sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			//std::cout<<"\tExtracting Reconstructed Pos Histograms: " <<hname <<"\n";
			_exp_data_pos_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
			_empty_pos_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
			
            //_empty_pos_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
			//sprintf(hname,"%s_%s_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
			sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			//std::cout<<"\tExtracting Reconstructed Neg Histograms: " <<hname <<"\n";
			_exp_data_neg_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            _empty_neg_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
			//_empty_neg_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
            _acceptance_5d[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
            _acceptance_5d[i][j]->Divide(_thrown_5d[i][j]);
            sprintf(hname,"Acceptance_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            //std::cout<<"\tMaking Acceptance Histogram: " <<hname <<"\n";
			_acceptance_5d[i][j]->SetNameTitle(hname,hname);
            for(long k=0; k<_acceptance_5d[i][j]->GetNbins(); k++){
				if(_acceptance_5d[i][j]->GetBinError(k)/_acceptance_5d[i][j]->GetBinContent(k) > 0.98){
					_acceptance_5d[i][j]->SetBinContent(k,0.0);
				}
			}
			sprintf(hname,"N_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            //std::cout<<"\tMaking Yield Histograms from clones: " <<hname <<"\n";
			//std::cout<<"\t\tnon-Helicity\n";
			_N_5d[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
			//std::cout<<"\t\tPos Helicity\n";
			_N_pos_5d[i][j] = (THnSparseD*)_exp_data_pos_5d[i][j]->Clone();
			//std::cout<<"\t\tNeg Helicity\n";
			_N_neg_5d[i][j] = (THnSparseD*)_exp_data_neg_5d[i][j]->Clone();
            //std::cout<<"\t\tNon-Hel Empty Subtraction\n";
			_N_5d[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
            _N_5d[i][j]->Divide(_acceptance_5d[i][j]);
			_N_5d[i][j]->SetNameTitle(hname,hname);
			//_N_pos_5d[i][j]->Add(_empty_pos_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
            _N_pos_5d[i][j]->Divide(_acceptance_5d[i][j]);
			_acceptance_5d[i][j]->SetNameTitle(hname,hname);
			//_N_neg_5d[i][j]->Add(_empty_neg_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
            _N_neg_5d[i][j]->Divide(_acceptance_5d[i][j]);
            //sprintf(hname,"Localized_Holes_50_W:%.3f-%.3f_Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            //_N_holes_5d[i][j] = (THnSparseD *)holes_->Get(hname);
	        //_N_5d[i][j]->Add(_N_holes_5d[i][j]);
			//_N_holes_5d[i][j] = (THnSparseD *)holes_->Get(hname);
			sprintf(hname,"Localized_Holes_50_W:%.3f-%.3f_Q2:%.2f-%.2f_pos",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			//std::cout<<"\tExtracting Localized Hole Histogram: " <<hname <<"\n";
			_N_holes_pos_5d[i][j] = (THnSparseD *)holes_->Get(hname);
			_N_pos_5d[i][j]->Add(_N_holes_pos_5d[i][j]);
			sprintf(hname,"N_pos_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			_N_pos_5d[i][j]->SetNameTitle(hname,hname);
			sprintf(hname,"Localized_Holes_50_W:%.3f-%.3f_Q2:%.2f-%.2f_neg",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			//std::cout<<"\tExtracting Localized Hole Histogram: " <<hname <<"\n";
			_N_holes_neg_5d[i][j] = (THnSparseD *)holes_->Get(hname);
			_N_neg_5d[i][j]->Add(_N_holes_neg_5d[i][j]);
			sprintf(hname,"N_neg_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			_N_neg_5d[i][j]->SetNameTitle(hname,hname);
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
            std::cout<<"Getting Thrown THnSparses " <<hname <<"\n";
			_thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _thrown_no_rad_5d[i][j] = (THnSparseD *)nr_sim_tree_->Get(hname);
            sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            std::cout<<"Getting Reconstructed THnSparses " <<hname <<"\n";
			_sim_data_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _exp_data_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            _empty_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
			//std::cout<<"Testing Integral methods for Experiment target filled\n";
			//std::cout<<"\tnSparseIntegral(): " <<fun::nSparseIntegral(_exp_data_5d[i][j]) <<"\n";
			//std::cout<<"\tSparse_Integral(): " <<fun::Sparse_Integral(_exp_data_5d[i][j]) <<"\n";
			//TH1D* tmp_hist = _exp_data_5d[i][j]->Projection(1);
			//std::cout<<"\tProj->Integral(): " <<tmp_hist->Integral() <<"\n";
			//std::cout<<"\t->Integral(): " <<_exp_data_5d[i][j]->Integral() <<"\n";
			//std::cout<<"Number of bins per axis\nSim\tExp\tEmp\tSim_no_rad\n";
			//for(int k=0; k<5; k++){
			//	std::cout<<_sim_data_5d[i][j]->GetAxis(k)->GetNbins() <<"\t" <<_exp_data_5d[i][j]->GetAxis(k)->GetNbins() <<"\t" <<_empty_5d[i][j]->GetAxis(k)->GetNbins() <<"\t" <<_thrown_no_rad_5d[i][j]->GetAxis(k)->GetNbins() <<"\n";
			//}
			std::cout<<"Making Acceptance \n";
            _acceptance_5d[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
			std::cout<<"\tDividing by Thrown\n";
            _acceptance_5d[i][j]->Divide(_thrown_5d[i][j]);
			for(long k=0; k<_acceptance_5d[i][j]->GetNbins(); k++){
				if(_acceptance_5d[i][j]->GetBinError(k)/_acceptance_5d[i][j]->GetBinContent(k) > 0.98){
					_acceptance_5d[i][j]->SetBinContent(k,0.0);
				}
			}
			std::cout<<"Making Yield\n";
            _N_5d[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
            _N_5d[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
	        _N_5d[i][j]->Divide(_acceptance_5d[i][j]);
			/*if(flags_.Flags::Nonlocal_Holes()){
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
	std::cout<<"Take 1\n";
	for (int i = 0; i < _W_nbins_; i++)
    {
        std::cout << " [" ;
        for (int j = 0; j < _Q2_nbins_; j++)
        {
            std::cout << _rad_corr->GetBinContent(i + 1, j + 1) << ", ";
        }
        std::cout<<"],"<<std::endl;
    }
	std::cout<<"Take 2\n";
	for (int i = 0; i < _W_nbins_; i++)
    {
        std::cout << " [" ;
        for (int j = 0; j < _Q2_nbins_; j++)
        {
            std::cout << _rad_corr_array[i][j] << ", ";
        }
        std::cout<<"],"<<std::endl;
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
    _rad_corr->Write();
    //std::cout<<"\t\tWent int first directory\n";
	TDirectory* dir_S1[5];
	TDirectory* dir_S2[5][_Q2_nbins_];
    TDirectory* dir_C1;
	char dirname[100];
    sprintf(dirname,"Integrated Check");
    dir_C1 = dir_S->mkdir(dirname);
    //std::cout<<"\tMade some internal directories\n\t\tNow naming them\n";
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
                //std::cout<<"\tpushing back projection\n";
				exp_ch_1d.push_back(_N_5d[i][j]->Projection(k,"E"));
				//std::cout<<"\tDenominator time\n\t\tvirtual photon flux\n";
                denom *= physics::Virtual_Photon_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]);
				//std::cout<<"Scaling Denominator post flux: " <<denom <<"\n";
				//std::cout<<"What it would have been with old flux: " <<denom/physics::ratio_Virtual_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]) <<"\n";
				//std::cout<<"\tThe ratio of difference is: " <<physics::ratio_Virtual_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]) <<"\n";
				//std::cout<<"Old virtual flux: " <<physics::old_Virtual_Photon_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]) <<" vs new: " <<physics::Virtual_Photon_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]) <<"\n";
                //std::cout<<"\t\tluminosity\n";
                denom *=flags_.Flags::L(0);
                //std::cout<<"Scaling Denominator post Luminosity: " <<denom <<"\n";
                //std::cout<<"\t\trad corr\n";
				if(flags_.Flags::Rad_Corr()){
					denom *= _rad_corr_array[i][j];
                    //std::cout<<"Scaling Denominator post rad corr: " <<denom <<"\n";
				}
                //std::cout<<"\t\tW bin\n";
				denom *=_W_res_;
                //std::cout<<"Scaling Denominator post W: " <<denom <<"\n";
                //std::cout<<"\t\tQ2 bin\n";
			    denom *=(_Q2_bins_[j+1]-_Q2_bins_[j]);
                //std::cout<<"Scaling Denominator post Q2: " <<denom <<"\n";
				if(k>1){
					//For Theta this will need to be undone and then modified for cosine theta
					denom *=(_thrown_5d[i][j]->GetAxis(k)->GetBinUpEdge(2)-_thrown_5d[i][j]->GetAxis(k)->GetBinLowEdge(2))*TMath::Pi()/180.0;//Divide by angle bins in radians
                    //denom *= Histogram::CosTheta()
                    //std::cout<<"Scaling Denominator post Xij: " <<denom <<"\n";
                }else{
                    //std::cout<<"\t\txij bin\n";
					denom *=(_thrown_5d[i][j]->GetAxis(k)->GetBinUpEdge(2)-_thrown_5d[i][j]->GetAxis(k)->GetBinLowEdge(2));//Phi in radians
                    //std::cout<<"Scaling Denominator post Xij: " <<denom <<"\n";
				}
				//std::cout<<"\tScaling\n";
                //std::cout<<"Current Integral" <<exp_ch_1d[k]->Integral() <<"\n";
                //std::cout<<"Scaling Denominator: " <<denom <<"\n";
				exp_ch_1d[k]->Scale(1.0/denom);
                //std::cout<<"Post Integral" <<exp_ch_1d[k]->Integral() <<"\n";
			
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

void Histogram::Beam_Spin(Flags flags_){
	//Will vary by W, Q2, and Phi 
	std::cout<<"Beam Spin Asymmetry\n";
	if(!flags_.Flags::Plot_Beam_Spin() || !flags_.Flags::Helicity()){
		std::cout<<"\tNot Plotting Beam Spin Asymmetry\n";
	 	return;
	}
	std::cout<<"Making Output File\n";
	_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	std::cout<<"Going into the rootfile: " <<_RootOutputFile->GetName() <<"\n";
	_RootOutputFile->cd();
	TCanvas* def = new TCanvas("def3");
	char hname[100];
	char dirname[100];
	sprintf(dirname,"ducsk_|%.2f-%.2f",Histogram::Q2_low(0),Histogram::Q2_top(0));
	std::cout<<"dirname: " <<dirname <<"\n";
	char xlabel[100];
	char ylabel[100];
	sprintf(xlabel,"Phi (deg)");
	sprintf(ylabel,"Asymmetry");
	TH1D_1d_star beam_spin_1d;
	TH1D_1d_star bot_spin_1d;
	TH1D_1d_star top_neg_1d;
	//TH3D * N_pos = _N_pos->THnSparse::Projection(0,1,6);
	//TH3D * N_neg = _N_neg->THnSparse::Projection(0,1,6);
	std::cout<<"\tMaking Directory\n";
	TDirectory* dir_B = _RootOutputFile->mkdir("Beam Spin");
	std::cout<<"\t\tMade main directory\n";
	//TDirectory* dir_B = _RootOutputFile.mkdir("Beam Spin");
	//dir_B->cd();
	std::cout<<"initialize sub directories \n";
	TDirectory* dir_B1[5];
	std::cout<<"Making sub directories\n";
	for(int j=0; j< 5; j++){//Q2
		sprintf(dirname,"Beam_Spin_Q2|%.2f-%.2f",Histogram::Q2_low(j),Histogram::Q2_top(j));
		std::cout<<"\t\tMaking Dir: " <<dirname <<"\n";
		dir_B1[j] = dir_B->mkdir(dirname);
	}
	int idx[3];
	std::cout<<"\tWriting Histograms\n";
	//THnSparseD* Beam_Spin_Top[_N_5d_pos.size()][_N_5d_pos[0].size()];
	//THnSparseD* Beam_Spin_Bot[_N_5d_pos.size()][_N_5d_pos[0].size()];
	//THnSparseD* Beam_Spin_R[_N_5d_pos.size()][_N_5d_pos[0].size()];
	//Convert the _exp_corr_5d and _exp_holes_5d to 4 dimensional sparse histograms for usage in plotting single differential cross sections
	for(int i=0; i<29; i++){//W
		for(int j=0; j< 5; j++){//Q2
			//std::cout<<"\t\tW:" <<i <<" Q2:" <<j <<"\r";
			sprintf(hname,"Beam_Spin_W|%.3f-%.3f_Q2|%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			//beam_spin_1d.push_back(new TH1D(hname,hname,_n_bins_5d[4],_bin_low_5d[4][0],_bin_up_5d[4][_n_bins_5d[4]-1]));
			beam_spin_1d.push_back(_N_pos_5d[i][j]->Projection(4,"E"));
			//beam_spin_1d.push_back(Beam_Spin_R[i][j]->Projection(4,"E"));
			top_neg_1d.push_back(_N_neg_5d[i][j]->Projection(4,"E"));
			beam_spin_1d[j]->SetNameTitle(hname,hname);
			beam_spin_1d[j]->Add(top_neg_1d[j],-1.0);
			sprintf(hname,"Top_neg_Beam_Spin_W|%.3f-%.3f_Q2|%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			top_neg_1d[j]->SetNameTitle(hname,hname);
			//beam_spin_1d[j]->Scale((1.0/_beam_pol_[0]));
			bot_spin_1d.push_back(_N_pos_5d[i][j]->Projection(4,"E"));
			bot_spin_1d[j]->Add(top_neg_1d[j],1.0);
			beam_spin_1d[j]->Divide(bot_spin_1d[j]);
			beam_spin_1d[j]->Scale((1.0/_beam_pol_[0]));
			//beam_spin_1d[j]->SetNameTitle(hname,hname);
			sprintf(hname,"Bot_Beam_Spin_W|%.3f-%.3f_Q2|%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			bot_spin_1d[j]->SetNameTitle(hname,hname);
			//std::cout<<"Checking Content for " <<hname <<"\n";
			//for(int k=0; k<_n_bins_5d[4]; k++){//Phi
			//	std::cout<<"\tPhi : " <<_bin_mid_7d[6][k] <<"| N_pos:" <<N_pos->GetBinContent(i,j,k) <<" | N_neg: " <<N_neg->GetBinContent(i,j,k) <<"\n";
			//	beam_spin_1d[beam_spin_1d.size()-1]->SetBinContent(k,(1.0/_beam_pol_[0])*(N_pos->GetBinContent(i,j,k)-N_neg->GetBinContent(i,j,k))/(N_pos->GetBinContent(i,j,k)+N_neg->GetBinContent(i,j,k)));//Beam polarization needs to be specified and is not currently accurate
			//}
			dir_B1[j]->cd();
			sprintf(hname,"Beam_Spin_W|%.3f-%.3f_Q2|%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			beam_spin_1d[j]->SetNameTitle(hname,hname);
			beam_spin_1d[j]->GetXaxis()->SetTitle(xlabel);
			beam_spin_1d[j]->GetYaxis()->SetTitle(ylabel);
			beam_spin_1d[j]->Write();
		}
		//_beam_spin_hist.push_back(beam_spin_1d);
		beam_spin_1d.clear();
		bot_spin_1d.clear();
		top_neg_1d.clear();
	}
	std::cout<<"\nFinished Beam Spin\n";
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
	return _W_min_+_W_res_*(W_bin_+1)-_MM_offset[var_set_];
}

double Histogram::MM2_max(int W_bin_, int var_set_){
	return _W_min_+_W_res_*(W_bin_+1)-_MM2_offset[var_set_];
}*/

double Histogram::CosTheta(int theta_bin_){
    double theta_res = ((double)_theta_max_ - (double)_theta_min_)/(double)_theta_bins_;
    return abs(TMath::Cos((_theta_min_+(theta_res*theta_bin_))*TMath::Pi()/180.0)-TMath::Cos((_theta_min_+(theta_res*theta_bin_+1))*TMath::Pi()/180.0));
}

