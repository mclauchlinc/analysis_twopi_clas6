#include "histogram.hpp"


Histogram::Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_){
    Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,flags_);
    //Histogram::Rad_Corr();
    //Histogram::Sparse_7to5(flags_);
    //Histogram::Single_Diff(flags_);
	Histogram::Localized_Holes(flags_,flags_.Flags::Min_Local_Dist(),20);
}


void Histogram::Extract_5d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_){
    std::cout<<"Extract 5d Histograms\n";
	char hname[100];
	sprintf(hname,"Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
	std::cout<<"Getting Thrown THnSparse " <<hname <<"\n";
	//int i = flags_.Flags::W_Bin();
	//int j = flags_.Flags::Q2_Bin();
    for(int i=0; i<_W_nbins_; i++){
        for(int j=0; j<_Q2_nbins_; j++){
            sprintf(hname,"Thrown_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _thrown_no_rad_5d[i][j] = (THnSparseD *)nr_sim_tree_->Get(hname);
            sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _sim_data_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _exp_data_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            //_empty_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            _acceptance_5d[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
            _acceptance_5d[i][j]->Divide(_thrown_5d[i][j]);
			std::cout<<"Acceptance Made, making relative error hist for it\n";
			Histogram::Fill_Acceptance_Rel_Hist();
			std::cout<<"Cutting Acceptance based on relative error cut\n";
			Histogram::Acceptance_Rel_Error_Cut()
			
            std::cout<<"Making the sync clone\n";
			_sim_exp_5d_sync[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
			sprintf(hname,"sync_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			_sim_exp_5d_sync[i][j]->SetNameTitle(hname,hname);
			_sim_exp_5d_sync[i][j]->Divide(_exp_data_5d[i][j]);
			if(_sim_exp_5d_sync[i][j]->GetNbins() >0){
				Int_t* coord = new Int_t[_sim_exp_5d_sync[i][j]->GetNdimensions()];
				for (Long64_t k = 0; k < _sim_exp_5d_sync[i][j]->GetNbins(); ++k) {
					if(_sim_exp_5d_sync[i][j]->GetBinContent(k,coord) > 0.0){
						//std::cout<<"\t\tInitial Bin content for bin " <<k <<" " <<_sim_exp_5d_sync[i][j]->GetBinContent(k,coord) <<" " <<_sim_exp_5d_sync[i][j]->GetBinError(k) <<"\n";
						_sim_exp_5d_sync[i][j]->SetBinContent(k,1.0);
						_sim_exp_5d_sync[i][j]->SetBinError(k,1.0);
						//std::cout<<"\tBin content for bin " <<k <<" " <<_sim_exp_5d_sync[i][j]->GetBinContent(k,coord) <<" " <<_sim_exp_5d_sync[i][j]->GetBinError(k) <<"\n";
					}
				}
			}
			std::cout<<"Made sync histogram fully for " <<i <<" " <<j <<"\n";
			_exp_data_5d_bop[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
			_exp_data_5d_bop[i][j]->Divide(_sim_exp_5d_sync[i][j]);
			sprintf(hname,"exp_bop_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			_exp_data_5d_bop[i][j]->SetNameTitle(hname,hname);
			_exp_data_5d_bop2[i][j] = (THnSparseD*)_exp_data_5d_bop[i][j]->Clone();
			_sim_data_5d_bop[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
			_sim_data_5d_bop[i][j]->Divide(_sim_exp_5d_sync[i][j]);
			sprintf(hname,"sim_bop_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			_sim_data_5d_bop[i][j]->SetNameTitle(hname,hname);
			_sim_data_5d_bop2[i][j] = (THnSparseD*)_sim_data_5d_bop[i][j]->Clone();

			_sim_holes_tmp_5d[i][j] = (THnSparseD*)_sim_data_5d_bop[i][j]->Clone();
			_sim_holes_tmp_5d[i][j]->Divide(_acceptance_5d[i][j]);
			std::cout<<"acceptance corr sim W:" <<i <<" Q2:" <<j <<" nbins:" <<_sim_holes_tmp_5d[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_sim_holes_tmp_5d[i][j]) <<"\n";
			std::cout<<"\tthrown n bins: " <<_thrown_5d[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_thrown_5d[i][j]) <<"\n";
			std::cout<<"\tsim n bins: " <<_sim_data_5d_bop[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_sim_data_5d_bop[i][j]) <<"\n";
			_sim_holes_5d[i][j] = (THnSparseD*)_thrown_5d[i][j]->Clone();
			_sim_holes_5d[i][j]->Add(_sim_holes_tmp_5d[i][j],-1.0);
			std::cout<<"\tSim holes nbins:" <<_sim_holes_5d[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_sim_holes_5d[i][j]) <<"\n";
			//_N_5d[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
            //_N_5d[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
	        //_N_5d[i][j]->Divide(_acceptance_5d[i][j]);
			std::cout<<"Checking for helicity\n";
			if(flags_.Flags::Helicity()){
				//std::cout<<"\t1\n";
				sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_exp_data_5d_pos[i][j] = (THnSparseD*)exp_tree_->Get(hname);
				//std::cout<<"\t2\n";
				_sim_exp_5d_sync_pos[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
				sprintf(hname,"sync_pos_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sim_exp_5d_sync_pos[i][j]->SetNameTitle(hname,hname);
				//std::cout<<"\t3\n";
				_sim_exp_5d_sync_pos[i][j]->Divide(_exp_data_5d_pos[i][j]);
				if(_sim_exp_5d_sync_pos[i][j]->GetNbins() >0){
					Int_t* coord_pos = new Int_t[_sim_exp_5d_sync_pos[i][j]->GetNdimensions()];
					for (Long64_t k = 0; k < _sim_exp_5d_sync_pos[i][j]->GetNbins(); ++k) {
						if(_sim_exp_5d_sync_pos[i][j]->GetBinContent(k,coord_pos) > 0.0){
							_sim_exp_5d_sync_pos[i][j]->SetBinContent(k,1.0);
							_sim_exp_5d_sync_pos[i][j]->SetBinError(k,1.0);
						}
					}
				
				}
				_exp_data_5d_bop_pos[i][j] = (THnSparseD*)_exp_data_5d_pos[i][j]->Clone();
				//std::cout<<"\t4\n";
				_exp_data_5d_bop_pos[i][j]->Divide(_sim_exp_5d_sync_pos[i][j]);
				sprintf(hname,"exp_bop_pos_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_exp_data_5d_bop_pos[i][j]->SetNameTitle(hname,hname);
				_exp_data_5d_bop2_pos[i][j] = (THnSparseD*)_exp_data_5d_bop_pos[i][j]->Clone();
				_sim_data_5d_bop_pos[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
				//std::cout<<"\t5\n";
				_sim_data_5d_bop_pos[i][j]->Divide(_sim_exp_5d_sync_pos[i][j]);
				sprintf(hname,"sim_bop_pos_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sim_data_5d_bop_pos[i][j]->SetNameTitle(hname,hname);
				_sim_data_5d_bop2_pos[i][j] = (THnSparseD*)_sim_data_5d_bop_pos[i][j]->Clone();

				_sim_holes_tmp_5d_pos[i][j] = (THnSparseD*)_sim_data_5d_bop_pos[i][j]->Clone();
				_sim_holes_tmp_5d_pos[i][j]->Divide(_acceptance_5d[i][j]);
				std::cout<<"acceptance corr sim W:" <<i <<" Q2:" <<j <<" nbins:" <<_sim_holes_tmp_5d[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_sim_holes_tmp_5d[i][j]) <<"\n";
				std::cout<<"\tthrown n bins: " <<_thrown_5d[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_thrown_5d[i][j]) <<"\n";
				std::cout<<"\tsim n bins: " <<_sim_data_5d_bop_pos[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_sim_data_5d_bop_pos[i][j]) <<"\n";
				_sim_holes_5d_pos[i][j] = (THnSparseD*)_thrown_5d[i][j]->Clone();
				_sim_holes_5d_pos[i][j]->Add(_sim_holes_tmp_5d_pos[i][j],-1.0);
				std::cout<<"\tSim holes nbins:" <<_sim_holes_5d_pos[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_sim_holes_5d_pos[i][j]) <<"\n";

				//Neg Time
				sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_exp_data_5d_neg[i][j] = (THnSparseD*)exp_tree_->Get(hname);
				_sim_exp_5d_sync_neg[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
				sprintf(hname,"sync_neg_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sim_exp_5d_sync_neg[i][j]->SetNameTitle(hname,hname);
				//std::cout<<"\t6\n";
				_sim_exp_5d_sync_neg[i][j]->Divide(_exp_data_5d_neg[i][j]);
				if(_sim_exp_5d_sync_neg[i][j]->GetNbins() >0){
					Int_t* coord_neg = new Int_t[_sim_exp_5d_sync_neg[i][j]->GetNdimensions()];
					for (Long64_t k = 0; k < _sim_exp_5d_sync_neg[i][j]->GetNbins(); ++k) {
						if(_sim_exp_5d_sync_neg[i][j]->GetBinContent(k,coord_neg) > 0.0){
							_sim_exp_5d_sync_neg[i][j]->SetBinContent(k,1.0);
							_sim_exp_5d_sync_neg[i][j]->SetBinError(k,1.0);
						}
					}
				}
				_exp_data_5d_bop_neg[i][j] = (THnSparseD*)_exp_data_5d_neg[i][j]->Clone();
				//std::cout<<"\t7\n";
				_exp_data_5d_bop_neg[i][j]->Divide(_sim_exp_5d_sync_neg[i][j]);
				sprintf(hname,"exp_bop_neg_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_exp_data_5d_bop_neg[i][j]->SetNameTitle(hname,hname);
				_exp_data_5d_bop2_neg[i][j] = (THnSparseD*)_exp_data_5d_bop_neg[i][j]->Clone();
				_sim_data_5d_bop_neg[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
				//std::cout<<"\t8\n";
				_sim_data_5d_bop_neg[i][j]->Divide(_sim_exp_5d_sync_neg[i][j]);
				sprintf(hname,"sim_bop_neg_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sim_data_5d_bop_neg[i][j]->SetNameTitle(hname,hname);
				_sim_data_5d_bop2_neg[i][j] = (THnSparseD*)_sim_data_5d_bop_neg[i][j]->Clone();

				_sim_holes_tmp_5d_neg[i][j] = (THnSparseD*)_sim_data_5d_bop_neg[i][j]->Clone();
				_sim_holes_tmp_5d_neg[i][j]->Divide(_acceptance_5d[i][j]);
				std::cout<<"acceptance corr sim W:" <<i <<" Q2:" <<j <<" nbins:" <<_sim_holes_tmp_5d[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_sim_holes_tmp_5d[i][j]) <<"\n";
				std::cout<<"\tthrown n bins: " <<_thrown_5d[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_thrown_5d[i][j]) <<"\n";
				std::cout<<"\tsim n bins: " <<_sim_data_5d_bop_neg[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_sim_data_5d_bop_neg[i][j]) <<"\n";
				_sim_holes_5d_neg[i][j] = (THnSparseD*)_thrown_5d[i][j]->Clone();
				_sim_holes_5d_neg[i][j]->Add(_sim_holes_tmp_5d_neg[i][j],-1.0);
				std::cout<<"\tSim holes nbins:" <<_sim_holes_5d_neg[i][j]->GetNbins() <<" and integral: " <<fun::nSparseIntegral(_sim_holes_5d_neg[i][j]) <<"\n";
			}
        }
    }
	//for(int k = 0; k<_thrown_5d[i][j]->GetNdimensions(); k++){
    for(int k = 0; k<_thrown_5d[0][0]->GetNdimensions(); k++){
		//_n_bins_7d.push_back(_thrown_7d->GetAxis(i)->GetNbins());
        //_n_bins_5d.push_back(_thrown_5d[i][j]->GetAxis(k)->GetNbins());
		_n_bins_5d.push_back(_thrown_5d[0][0]->GetAxis(k)->GetNbins());
	}
   

}


void Histogram::Localized_Holes(Flags flags_, int min_dist_, int max_dist_){
	if(!flags_.Flags::Plot_Localized_Holes()){
		return;
	}
	_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	_RootOutputFile->cd();
	std::vector<long> space_dims;
	char hname[100];
	std::cout<<"\n";
	for(int i=0; i<_n_bins_5d.size(); i++){
      	space_dims.push_back(_n_bins_5d[i]);
		std::cout<<"space_dims[" <<i <<"]:" <<space_dims[i] <<"\n";
	}
	_localized_hole_filling = false;
	CartesianGenerator cart(space_dims);
	int bin[_n_bins_5d.size()];
   	int bin2[_n_bins_5d.size()];
	bool look_further_all = false;
	bool look_further_pos = false;
	bool look_further_neg = false;
	int dist = 0; 
	int bin_low[5];
	int bin_top[5];
	int bin_low2[5];
	int bin_top2[5];
	std::vector<std::vector<int>> surr_bins; 
	long dist_dist[3][20] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}; 
	long total_ev[29][5];
	long curr_ev[29][5];
	int proj_bins[4] = {2,3,4,5};
	double loc_exp;
	double loc_sim;
	double loc_exp_pos;
	double loc_sim_pos;
	double loc_exp_neg;
	double loc_sim_neg;
	double loc_exp_err2;
	double loc_sim_err2;
	double loc_exp_err2_pos;
	double loc_sim_err2_pos;
	double loc_exp_err2_neg;
	double loc_sim_err2_neg; 
	int wq2_bin = 0;
	int binning[5];

	double rel_err = 0.0;
	double rel_err1 = 0.0;
	double rel_err2 = 0.0;
	double rel_err3 = 0.0;
	double rel_err1_pos = 0.0;
	double rel_err2_pos = 0.0;
	double rel_err3_pos = 0.0;
	double rel_err_pos = 0.0;
	double rel_err1_neg = 0.0;
	double rel_err2_neg = 0.0;
	double rel_err3_neg = 0.0;
	double rel_err_neg = 0.0;

	double sf_rel_err1;
	double sf_rel_err2;
	double sf_rel_err;
	double sf_rel_err1_pos;
	double sf_rel_err2_pos;
	double sf_rel_err_pos;
	double sf_rel_err1_neg;
	double sf_rel_err2_neg;
	double sf_rel_err_neg;


	//TH1D* tmp_exp;
	//TH1D* tmp_sim;
	//TH1D* tmp_exp_pos;
	//TH1D* tmp_sim_pos;
	//TH1D* tmp_exp_neg;
	//TH1D* tmp_sim_neg;
	
	//int i = flags_.Flags::W_Bin();
	//int j = flags_.Flags::Q2_Bin();
	//std::cout<<"Looking at W bin:" <<i <<" and Q2 bin:" <<j <<"\n";
	for(int i=0; i<_W_nbins_; i++){
		for(int j=0; j<_Q2_nbins_; j++){
			if(flags_.Flags::Helicity()){
				sprintf(hname,"pos scale factor distribution W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sf_hist_pos[i][j] = new TH1D(hname,hname,500,1.0,10000);
				sprintf(hname,"neg scale factor distribution W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sf_hist_neg[i][j] = new TH1D(hname,hname,500,1.0,10000);
				sprintf(hname,"pos localized holes relative error W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_err_hist_pos[i][j] = new TH1D(hname,hname,400,0.0,4.0);
				sprintf(hname,"neg localized holes relative error W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_err_hist_neg[i][j] = new TH1D(hname,hname,400,0.0,4.0);
				sprintf(hname,"pos localized radius W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_radius_hist_pos[i][j] = new TH1D(hname,hname,20,0.5,20.5);
				sprintf(hname,"neg localized radius W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_radius_hist_neg[i][j] = new TH1D(hname,hname,20,0.5,20.5);
			}else{
				sprintf(hname,"scale factor distribution W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sf_hist[i][j] = new TH1D(hname,hname,500,1.0,10000);
				sprintf(hname,"Acceptance Distribution W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_acceptance_dist[i][j] = new TH1D(hname,hname,200,0.0,0.2);
				sprintf(hname,"Relative Acceptance Distribution W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_rel_acceptance_dist[i][j] = new TH1D(hname,hname,100,0.0,1.5);
				sprintf(hname,"localized holes relative error W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_err_hist[i][j] = new TH1D(hname,hname,400,0.0,4.0);
				sprintf(hname,"localized radius W:%.3f-%.3f Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_radius_hist[i][j] = new TH1D(hname,hname,20,0.5,20.5);
				while(cart.GetNextCombination()){
					for(int q=0; q<5; q++){
						binning[q] = cart[q]+1;
					}
					if(_exp_data_5d[i][j]->GetBinContent(binning)){
						_acceptance_dist[i][j]->Fill(_acceptance_5d[i][j]->GetBinContent(binning));
					}
				}
				for (Long64_t l = 0; l < _acceptance_5d[i][j]->GetNbins(); ++l) {
					_rel_acceptance_dist[i][j]->Fill(_acceptance_5d[i][j]->GetBinError(l)/_acceptance_5d[i][j]->GetBinContent(l));
				}
				_acceptance_dist[i][j]->Write();
				_rel_acceptance_dist[i][j]->Write();
			}
		}
	}

	std::cout<<"*****************The helicity flag is " <<flags_.Flags::Helicity() <<"\n";
	for(int i=0; i<_W_nbins_; i++){
		for(int j=0; j<_Q2_nbins_; j++){
			for(int k=0; k<3; k++){
				for(int l=0; l<20; l++){
					dist_dist[k][l] = 0;
				}
			}
			total_ev[i][j] = _sim_holes_5d[i][j]->GetNbins();
			curr_ev[i][j] = 0;
			_scale_exp_5d[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
			_scale_exp_5d[i][j]->Scale(0.0);
			_scale_sim_5d[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
			_scale_sim_5d[i][j]->Scale(0.0);
			_scale_exp_5d_pos[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
			_scale_exp_5d_pos[i][j]->Scale(0.0);
			_scale_exp_5d_neg[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
			_scale_exp_5d_neg[i][j]->Scale(0.0);
			_scale_sim_5d_pos[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
			_scale_sim_5d_pos[i][j]->Scale(0.0);
			_scale_sim_5d_neg[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
			_scale_sim_5d_neg[i][j]->Scale(0.0);
			Int_t* coord = new Int_t[_sim_holes_5d[i][j]->GetNbins()];
			bool dubbin[5] = {false,false,false,false,false};
			while(cart.GetNextCombination()){
				dist = 0;
				for(int k = 0; k<_n_bins_5d.size(); k++){
					bin[k] = cart[k]+1;
				}
				if(flags_.Flags::Helicity()){
					if(_sim_holes_5d_pos[i][j]->GetBinContent(bin)>0.0){
						look_further_pos = true;
					}else{
						look_further_pos = false;
					}
					if(_sim_holes_5d_neg[i][j]->GetBinContent(bin)>0.0){
						look_further_neg = true;
					}else{
						look_further_neg = false;
					}
				}else{
					if(_sim_holes_5d[i][j]->GetBinContent(bin)>0.0){// && bin[0]>=3 && bin[0]<=16 && bin[1]>=3 && bin[1]<=16){
						look_further_all = true;
					}else{
						look_further_all = false;
				}
				}
				if(look_further_all || look_further_pos || look_further_neg){
					curr_ev[i][j] ++; 
					if((curr_ev[i][j]-1)%(total_ev[i][j]/100) == 0){
						std::cout<<"\r" <<"\t" <<(100*curr_ev[i][j]/total_ev[i][j]) <<"/100"  <<std::flush ;
					}
					if(flags_.Flags::Helicity()){
						look_further_all = false;
						//look_further_pos = true;
						//look_further_neg = true;
					}else{
						//look_further_all = true;
						look_further_pos = false;
						look_further_neg = false;
					}
					while(look_further_all || look_further_pos || look_further_neg){
						dist++;
						if(dist >= min_dist_){
							for(int n=0; n<5; n++){
								dubbin[n] = false;
							}
							for(int n=0; n<5; n++){
								if((bin[n]+dist) > _n_bins_5d[n]){
									bin_top[n] = _n_bins_5d[n];
								}else{
									bin_top[n] = (bin[n]+dist);
								}
								if((bin[n]-dist) < 1){
									bin_low[n] = 1;
								}else{
									bin_low[n] = (bin[n]-dist);
								}
								if(n >2){
									if((bin[n]+dist) > _n_bins_5d[n]){
										if((bin[n]-dist) < 1){
											bin_low2[n] = 1;
											bin_top2[n] = _n_bins_5d[n];
										}else{
											bin_low2[n] = 1;
											bin_top2[n] = (bin[n]+dist)-_n_bins_5d[n];
											dubbin[n] = true;
										}
									}else{
										if((bin[n]-dist) < 1){
											bin_low2[n] = (bin[n]-dist)+_n_bins_5d[n]-1;
											bin_top2[n] = _n_bins_5d[n];
											dubbin[n] = true;
										}else{
											bin_low2[n] = (bin[n]-dist);
											bin_top2[n] = (bin[n]+dist);
										}
									}
								}
							}
							if(dubbin[3] || dubbin[4]){
								for(int n=0; n<5; n++){
									if(!dubbin[n]){
										bin_low2[n] = bin_low[n];
										bin_top2[n] = bin_top[n];
									}
								}
							}
							for(int n=0; n<5; n++){
								//_sim_data_5d[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
								
								if(flags_.Flags::Helicity()){
									//_exp_data_5d_pos[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
									//_exp_data_5d_neg[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
									_exp_data_5d_bop_pos[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
									_exp_data_5d_bop_neg[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
									_sim_data_5d_bop_pos[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
									_sim_data_5d_bop_neg[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
								}else{
									//_exp_data_5d[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
									_sim_data_5d_bop[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
									_exp_data_5d_bop[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
								}
								if(dubbin[3] || dubbin[4]){
									
									if(flags_.Flags::Helicity()){
										//_exp_data_5d_pos[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
										//_exp_data_5d_neg[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
										_exp_data_5d_bop2_pos[i][j]->GetAxis(n)->SetRange(bin_low2[n],bin_top2[n]);
										_exp_data_5d_bop2_neg[i][j]->GetAxis(n)->SetRange(bin_low2[n],bin_top2[n]);
										_sim_data_5d_bop2_pos[i][j]->GetAxis(n)->SetRange(bin_low2[n],bin_top2[n]);
										_sim_data_5d_bop2_neg[i][j]->GetAxis(n)->SetRange(bin_low2[n],bin_top2[n]);
									}else{
										//_exp_data_5d[i][j]->GetAxis(n)->SetRange(bin_low[n],bin_top[n]);
										_exp_data_5d_bop2[i][j]->GetAxis(n)->SetRange(bin_low2[n],bin_top2[n]);
										_sim_data_5d_bop2[i][j]->GetAxis(n)->SetRange(bin_low2[n],bin_top2[n]);
									}
								}
								//std::cout<<"In  bin:" <<bin[n] <<" Low bin:" <<bin_low[n] <<" Top bin:" <<bin_top[n] <<"\n";
							}
							if(flags_.Flags::Helicity()){
								if(look_further_pos){
									loc_exp_pos = 0.0;
									loc_sim_pos = 0.0;
									loc_exp_err2_pos = 0.0;
									loc_sim_err2_pos = 0.0;
									//if(_exp_data_5d_pos[i][j]->GetNbins() > 0){
									if(_exp_data_5d_bop_pos[i][j]->GetNbins() > 0){
										//for( long v_pos_idx = 0; v_pos_idx< _exp_data_5d_pos[i][j]->GetNbins(); v_pos_idx++){
										for( long v_pos_idx = 0; v_pos_idx< _exp_data_5d_bop_pos[i][j]->GetNbins(); v_pos_idx++){
											//Double_t v_pos = _exp_data_5d_pos[i][j]->GetBinContent(v_pos_idx,coord);
											Double_t v_pos = _exp_data_5d_bop_pos[i][j]->GetBinContent(v_pos_idx,coord);
											bool pass_pos = true;
											for(int z = 0; z<5; z++){
												pass_pos &= coord[z] >= bin_low[z];
												pass_pos &= coord[z] <= bin_top[z];
											}
											if(pass_pos){
												if(v_pos >0.0){
													//Double_t w_pos = _sim_data_5d[i][j]->GetBinContent(coord);
													Double_t w_pos = _sim_data_5d_bop_pos[i][j]->GetBinContent(coord);
													if(w_pos > 0.0){
														loc_exp_pos += v_pos;
														loc_sim_pos += w_pos;
														//loc_exp_err2_pos += _exp_data_5d_pos[i][j]->GetBinError2(v_pos_idx);
														loc_exp_err2_pos += _exp_data_5d_bop_pos[i][j]->GetBinError2(v_pos_idx);
														//loc_sim_err2_pos += _sim_data_5d[i][j]->GetBinError2(_sim_data_5d[i][j]->GetBin(coord));
														loc_sim_err2_pos += _sim_data_5d_bop_pos[i][j]->GetBinError2(_sim_data_5d_bop_pos[i][j]->GetBin(coord));
													}
												}
											}
										}
										//rel_err1_pos = loc_exp_err2_pos*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_pos*loc_sim_pos);
										//rel_err2_pos = loc_exp_pos*loc_exp_pos*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim_pos*loc_sim_pos);
										//rel_err3_pos = loc_sim_err2_pos*loc_exp_pos*loc_exp_pos*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_pos*loc_sim_pos*loc_sim_pos*loc_sim_pos);
										//rel_err_pos = TMath::Sqrt(rel_err1_pos+rel_err2_pos+rel_err3_pos)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp_pos/loc_sim_pos);
											
									}
									if(dubbin){
										if(_exp_data_5d_bop2_pos[i][j]->GetNbins() > 0){
											//for( long v_pos_idx = 0; v_pos_idx< _exp_data_5d_pos[i][j]->GetNbins(); v_pos_idx++){
											for( long v_pos_idx = 0; v_pos_idx< _exp_data_5d_bop2_pos[i][j]->GetNbins(); v_pos_idx++){
												//Double_t v_pos = _exp_data_5d_pos[i][j]->GetBinContent(v_pos_idx,coord);
												Double_t v_pos = _exp_data_5d_bop2_pos[i][j]->GetBinContent(v_pos_idx,coord);
												bool pass_pos = true;
												for(int z = 0; z<5; z++){
													pass_pos &= coord[z] >= bin_low2[z];
													pass_pos &= coord[z] <= bin_top2[z];
												}
												if(pass_pos){
													if(v_pos >0.0){
														//Double_t w_pos = _sim_data_5d[i][j]->GetBinContent(coord);
														Double_t w_pos = _sim_data_5d_bop2_pos[i][j]->GetBinContent(coord);
														if(w_pos > 0.0){
															loc_exp_pos += v_pos;
															loc_sim_pos += w_pos;
															//loc_exp_err2_pos += _exp_data_5d_pos[i][j]->GetBinError2(v_pos_idx);
															loc_exp_err2_pos += _exp_data_5d_bop_pos[i][j]->GetBinError2(v_pos_idx);
															//loc_sim_err2_pos += _sim_data_5d[i][j]->GetBinError2(_sim_data_5d[i][j]->GetBin(coord));
															loc_sim_err2_pos += _sim_data_5d_bop2_pos[i][j]->GetBinError2(_sim_data_5d_bop2_pos[i][j]->GetBin(coord));
														}
													}
												}
											}
											//rel_err1_pos = loc_exp_err2_pos*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_pos*loc_sim_pos);
											//rel_err2_pos = loc_exp_pos*loc_exp_pos*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim_pos*loc_sim_pos);
											//rel_err3_pos = loc_sim_err2_pos*loc_exp_pos*loc_exp_pos*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_pos*loc_sim_pos*loc_sim_pos*loc_sim_pos);
											//rel_err_pos = TMath::Sqrt(rel_err1_pos+rel_err2_pos+rel_err3_pos)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp_pos/loc_sim_pos);
											
											//rel_err1_pos = loc_exp_err2_pos*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_pos*loc_sim_pos);
											//rel_err2_pos = loc_exp_pos*loc_exp_pos*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim_pos*loc_sim_pos);
											//rel_err3_pos = loc_sim_err2_pos*loc_exp_pos*loc_exp_pos*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_pos*loc_sim_pos*loc_sim_pos*loc_sim_pos);
											//rel_err_pos = TMath::Sqrt(rel_err1_pos+rel_err2_pos+rel_err3_pos)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp_pos/loc_sim_pos);

											//sf_rel_err1_pos = loc_exp_err2_pos/(loc_sim_pos*loc_sim_pos);
											//sf_rel_err2_pos = loc_sim_err2_pos*loc_exp_pos*loc_exp_pos/(loc_sim_pos*loc_sim_pos*loc_sim_pos*loc_sim_pos);
											//sf_rel_err_pos = TMath::Sqrt(sf_rel_err1_pos+sf_rel_err2_pos)/(loc_exp_pos/loc_sim_pos);
										}
									}
									if(_exp_data_5d_bop_pos[i][j]->GetNbins() > 0 || (dubbin && _exp_data_5d_bop2_pos[i][j]->GetNbins() > 0)){
										rel_err1_pos = loc_exp_err2_pos*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_pos*loc_sim_pos);
										rel_err2_pos = loc_exp_pos*loc_exp_pos*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim_pos*loc_sim_pos);
										rel_err3_pos = loc_sim_err2_pos*loc_exp_pos*loc_exp_pos*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_pos*loc_sim_pos*loc_sim_pos*loc_sim_pos);
										rel_err_pos = TMath::Sqrt(rel_err1_pos+rel_err2_pos+rel_err3_pos)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp_pos/loc_sim_pos);

										sf_rel_err1_pos = loc_exp_err2_pos/(loc_sim_pos*loc_sim_pos);
										sf_rel_err2_pos = loc_sim_err2_pos*loc_exp_pos*loc_exp_pos/(loc_sim_pos*loc_sim_pos*loc_sim_pos*loc_sim_pos);
										sf_rel_err_pos = TMath::Sqrt(sf_rel_err1_pos+sf_rel_err2_pos)/(loc_exp_pos/loc_sim_pos);
									}
									//loc_exp_pos = fun::nSparseIntegral(_exp_data_5d_pos[i][j]);
									//loc_sim_pos = fun::nSparseIntegral(_sim_data_5d[i][j]);
									//if(loc_exp_pos>0.0 && loc_sim_pos>0.0 && (rel_err_pos <= 0.5 || dist>=max_dist_)){
									if(loc_exp_pos>0.0 && loc_sim_pos>0.0 && (sf_rel_err_pos <= 0.2 || dist>=max_dist_)){	
										//std::cout<<"exp_pos: " <<loc_exp_pos <<"  sim_pos: " <<loc_sim_pos <<" for bin:" <<bin[0] <<" " <<bin[1] <<" " <<bin[2] <<" " <<bin[3] <<" " <<bin[4] <<"\n";
										//std::cout<<"\tIntegrating exp pos\n";
										_scale_exp_5d_pos[i][j]->SetBinContent(bin,loc_exp_pos);
										//std::cout<<"\tIntegrating sim pos\n";
										_scale_sim_5d_pos[i][j]->SetBinContent(bin,loc_sim_pos);
										look_further_pos=false;
										_scale_exp_5d_pos[i][j]->SetBinError2(_scale_exp_5d_pos[i][j]->GetBin(bin),loc_exp_err2_pos);
										_scale_sim_5d_pos[i][j]->SetBinError2(_scale_sim_5d_pos[i][j]->GetBin(bin),loc_sim_err2_pos);
										dist_dist[1][dist-1]++;
										_hole_radius_hist_pos[i][j]->Fill(dist);
									}
								}
								if(look_further_neg){
									loc_exp_neg = 0.0;
									loc_sim_neg = 0.0;
									loc_exp_err2_neg = 0.0;
									loc_sim_err2_neg = 0.0;
									//if(_exp_data_5d_neg[i][j]->GetNbins() > 0){
									if(_exp_data_5d_bop_neg[i][j]->GetNbins() > 0){
										//for( long v_neg_idx = 0; v_neg_idx< _exp_data_5d_neg[i][j]->GetNbins(); v_neg_idx++){
										for( long v_neg_idx = 0; v_neg_idx< _exp_data_5d_bop_neg[i][j]->GetNbins(); v_neg_idx++){
											//Double_t v_neg = _exp_data_5d_neg[i][j]->GetBinContent(v_neg_idx,coord);
											Double_t v_neg = _exp_data_5d_bop_neg[i][j]->GetBinContent(v_neg_idx,coord);
											bool pass_neg = true;
											for(int z = 0; z<5; z++){
												pass_neg &= coord[z] >= bin_low[z];
												pass_neg &= coord[z] <= bin_top[z];
											}
											if(pass_neg){
												if(v_neg >0.0){
													//Double_t w_neg = _sim_data_5d[i][j]->GetBinContent(coord);
													Double_t w_neg = _sim_data_5d_bop_neg[i][j]->GetBinContent(coord);
													if(w_neg > 0.0){
														loc_exp_neg += v_neg;
														loc_sim_neg += w_neg;
														//loc_exp_err2_neg += _exp_data_5d_neg[i][j]->GetBinError2(v_neg_idx);
														//loc_sim_err2_neg += _sim_data_5d[i][j]->GetBinError2(_sim_data_5d[i][j]->GetBin(coord));
														loc_exp_err2_neg += _exp_data_5d_bop_neg[i][j]->GetBinError2(v_neg_idx);
														loc_sim_err2_neg += _sim_data_5d_bop_neg[i][j]->GetBinError2(_sim_data_5d_bop_neg[i][j]->GetBin(coord));
													}
												}
											}
										}
										//rel_err1_neg = loc_exp_err2_neg*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_neg*loc_sim_neg);
										//rel_err2_neg = loc_exp_neg*loc_exp_neg*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim_neg*loc_sim_neg);
										//rel_err3_neg = loc_sim_err2_neg*loc_exp_neg*loc_exp_neg*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_neg*loc_sim_neg*loc_sim_neg*loc_sim_neg);
										//rel_err_neg = TMath::Sqrt(rel_err1_neg+rel_err2_neg+rel_err3_neg)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp_neg/loc_sim_neg);

										//sf_rel_err1_neg = loc_exp_err2_neg/(loc_sim_neg*loc_sim_neg);
										//sf_rel_err2_neg = loc_sim_err2_neg*loc_exp_neg*loc_exp_neg/(loc_sim_neg*loc_sim_neg*loc_sim_neg*loc_sim_neg);
										//sf_rel_err_neg = TMath::Sqrt(sf_rel_err1_neg+sf_rel_err2_neg)/(loc_exp_neg/loc_sim_neg);
									}
									if(dubbin){
										if(_exp_data_5d_bop2_neg[i][j]->GetNbins() > 0){
											//for( long v_neg_idx = 0; v_neg_idx< _exp_data_5d_neg[i][j]->GetNbins(); v_neg_idx++){
											for( long v_neg_idx = 0; v_neg_idx< _exp_data_5d_bop2_neg[i][j]->GetNbins(); v_neg_idx++){
												//Double_t v_neg = _exp_data_5d_neg[i][j]->GetBinContent(v_neg_idx,coord);
												Double_t v_neg = _exp_data_5d_bop2_neg[i][j]->GetBinContent(v_neg_idx,coord);
												bool pass_neg = true;
												for(int z = 0; z<5; z++){
													pass_neg &= coord[z] >= bin_low2[z];
													pass_neg &= coord[z] <= bin_top2[z];
												}
												if(pass_neg){
													if(v_neg >0.0){
														//Double_t w_neg = _sim_data_5d[i][j]->GetBinContent(coord);
														Double_t w_neg = _sim_data_5d_bop2_neg[i][j]->GetBinContent(coord);
														if(w_neg > 0.0){
															loc_exp_neg += v_neg;
															loc_sim_neg += w_neg;
															//loc_exp_err2_neg += _exp_data_5d_neg[i][j]->GetBinError2(v_neg_idx);
															//loc_sim_err2_neg += _sim_data_5d[i][j]->GetBinError2(_sim_data_5d[i][j]->GetBin(coord));
															loc_exp_err2_neg += _exp_data_5d_bop2_neg[i][j]->GetBinError2(v_neg_idx);
															loc_sim_err2_neg += _sim_data_5d_bop2_neg[i][j]->GetBinError2(_sim_data_5d_bop2_neg[i][j]->GetBin(coord));
														}
													}
												}
											}
											//rel_err1_neg = loc_exp_err2_neg*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_neg*loc_sim_neg);
											//rel_err2_neg = loc_exp_neg*loc_exp_neg*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim_neg*loc_sim_neg);
											//rel_err3_neg = loc_sim_err2_neg*loc_exp_neg*loc_exp_neg*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_neg*loc_sim_neg*loc_sim_neg*loc_sim_neg);
											//rel_err_neg = TMath::Sqrt(rel_err1_neg+rel_err2_neg+rel_err3_neg)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp_neg/loc_sim_neg);

											//sf_rel_err1_neg = loc_exp_err2_neg/(loc_sim_neg*loc_sim_neg);
											//sf_rel_err2_neg = loc_sim_err2_neg*loc_exp_neg*loc_exp_neg/(loc_sim_neg*loc_sim_neg*loc_sim_neg*loc_sim_neg);
											//sf_rel_err_neg = TMath::Sqrt(sf_rel_err1_neg+sf_rel_err2_neg)/(loc_exp_neg/loc_sim_neg);
										}
									}
									if(_exp_data_5d_bop_neg[i][j]->GetNbins() > 0 || (dubbin && _exp_data_5d_bop2_neg[i][j]->GetNbins() > 0)){
										rel_err1_neg = loc_exp_err2_neg*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_neg*loc_sim_neg);
										rel_err2_neg = loc_exp_neg*loc_exp_neg*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim_neg*loc_sim_neg);
										rel_err3_neg = loc_sim_err2_neg*loc_exp_neg*loc_exp_neg*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_neg*loc_sim_neg*loc_sim_neg*loc_sim_neg);
										rel_err_neg = TMath::Sqrt(rel_err1_neg+rel_err2_neg+rel_err3_neg)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp_neg/loc_sim_neg);

										sf_rel_err1_neg = loc_exp_err2_neg/(loc_sim_neg*loc_sim_neg);
										sf_rel_err2_neg = loc_sim_err2_neg*loc_exp_neg*loc_exp_neg/(loc_sim_neg*loc_sim_neg*loc_sim_neg*loc_sim_neg);
										sf_rel_err_neg = TMath::Sqrt(sf_rel_err1_neg+sf_rel_err2_neg)/(loc_exp_neg/loc_sim_neg);
									}
									//loc_exp_neg = fun::nSparseIntegral(_exp_data_5d_neg[i][j]);
									//loc_sim_neg = fun::nSparseIntegral(_sim_data_5d[i][j]);
									//if(loc_exp_neg>0.0 && loc_sim_neg>0.0 && (rel_err_neg <= 0.5 || dist>=max_dist_) ){
									if(loc_exp_neg>0.0 && loc_sim_neg>0.0 && (sf_rel_err_neg <= 0.2 || dist>=max_dist_) ){
										//std::cout<<"exp_neg: " <<loc_exp_neg <<"  sim_neg: " <<loc_sim_neg <<" for bin:" <<bin[0] <<" " <<bin[1] <<" " <<bin[2] <<" " <<bin[3] <<" " <<bin[4] <<"\n";
										//std::cout<<"\tIntegrating exp neg\n";
										_scale_exp_5d_neg[i][j]->SetBinContent(bin,loc_exp_neg);
										//std::cout<<"\tIntegrating sim neg\n";
										_scale_sim_5d_neg[i][j]->SetBinContent(bin,loc_sim_neg);
										look_further_neg=false;
										_scale_exp_5d_neg[i][j]->SetBinError2(_scale_exp_5d_neg[i][j]->GetBin(bin),loc_exp_err2_neg);
										_scale_sim_5d_neg[i][j]->SetBinError2(_scale_sim_5d_neg[i][j]->GetBin(bin),loc_sim_err2_neg);
										dist_dist[2][dist-1]++;
										_hole_radius_hist_neg[i][j]->Fill(dist);
										
									}
								}
							}else{
								if(look_further_all){
									loc_exp= 0.0;
									loc_sim = 0.0;
									loc_exp_err2= 0.0;
									loc_sim_err2 = 0.0;
									//if(_exp_data_5d[i][j]->GetNbins() > 0){
									if(_exp_data_5d_bop[i][j]->GetNbins() > 0){
										//for( long v_idx = 0; v_idx< _exp_data_5d[i][j]->GetNbins(); v_idx++){
										for( long v_idx = 0; v_idx< _exp_data_5d_bop[i][j]->GetNbins(); v_idx++){
											//Double_t v = _exp_data_5d[i][j]->GetBinContent(v_idx,coord);
											Double_t v = _exp_data_5d_bop[i][j]->GetBinContent(v_idx,coord);
											bool pass = true;
											for(int z = 0; z<5; z++){
												pass &= coord[z] >= bin_low[z];
												pass &= coord[z] <= bin_top[z];
											}
											if(pass){
												if(v >0.0){
													//Double_t w = _sim_data_5d[i][j]->GetBinContent(coord);
													Double_t w = _sim_data_5d_bop[i][j]->GetBinContent(coord);
													if(w > 0.0){
														loc_exp += v;
														loc_sim += w;
														//loc_exp_err2 += _exp_data_5d[i][j]->GetBinError2(v_idx);
														//loc_sim_err2 += _sim_data_5d[i][j]->GetBinError2(_sim_data_5d[i][j]->GetBin(coord));
														loc_exp_err2 += _exp_data_5d_bop[i][j]->GetBinError2(v_idx);
														loc_sim_err2 += _sim_data_5d_bop[i][j]->GetBinError2(_sim_data_5d_bop[i][j]->GetBin(coord));
													}
												}
											}//else{
											//	std::cout<<"improper coordinate: " <<coord[0] <<" " <<coord[1] <<" " <<coord[2] <<" " <<coord[3] <<" " <<coord[4] <<"\n";
											//}
										}
										//rel_err1 = loc_exp_err2*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim);
										//rel_err2 = loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim*loc_sim);
										//rel_err3 = loc_sim_err2*loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim*loc_sim*loc_sim);
										//rel_err = TMath::Sqrt(rel_err1+rel_err2+rel_err3)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp/loc_sim);
										//rel_err = TMath::Sqrt((loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim*loc_sim) ) + (loc_exp_err2*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim) ) + (loc_sim_err2*loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim*loc_sim*loc_sim) ))/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp/loc_sim);
										//std::cout<<"relative Error for dist " <<dist <<": 1 " <<rel_err1 <<" 2 " <<rel_err2 <<" 3 " <<rel_err3 <<" rel "<<rel_err <<" sf:" <<loc_exp/loc_sim <<" hole val:" <<_sim_holes_5d[i][j]->GetBinContent(bin) <<"\n";

										//sf_rel_err1 = loc_exp_err2/(loc_sim*loc_sim);
										//sf_rel_err2 = loc_sim_err2*loc_exp*loc_exp/(loc_sim*loc_sim*loc_sim*loc_sim);
										//sf_rel_err = TMath::Sqrt(sf_rel_err1+sf_rel_err2)/(loc_exp/loc_sim);
									}
									if(dubbin){
										if(_exp_data_5d_bop2[i][j]->GetNbins() > 0){
											//for( long v_idx = 0; v_idx< _exp_data_5d[i][j]->GetNbins(); v_idx++){
											for( long v_idx = 0; v_idx< _exp_data_5d_bop[i][j]->GetNbins(); v_idx++){
												//Double_t v = _exp_data_5d[i][j]->GetBinContent(v_idx,coord);
												Double_t v = _exp_data_5d_bop[i][j]->GetBinContent(v_idx,coord);
												bool pass = true;
												for(int z = 0; z<5; z++){
													pass &= coord[z] >= bin_low2[z];
													pass &= coord[z] <= bin_top2[z];
												}
												if(pass){
													if(v >0.0){
														//Double_t w = _sim_data_5d[i][j]->GetBinContent(coord);
														Double_t w = _sim_data_5d_bop2[i][j]->GetBinContent(coord);
														if(w > 0.0){
															loc_exp += v;
															loc_sim += w;
															//loc_exp_err2 += _exp_data_5d[i][j]->GetBinError2(v_idx);
															//loc_sim_err2 += _sim_data_5d[i][j]->GetBinError2(_sim_data_5d[i][j]->GetBin(coord));
															loc_exp_err2 += _exp_data_5d_bop2[i][j]->GetBinError2(v_idx);
															loc_sim_err2 += _sim_data_5d_bop2[i][j]->GetBinError2(_sim_data_5d_bop2[i][j]->GetBin(coord));
														}
													}
												}//else{
												//	std::cout<<"improper coordinate: " <<coord[0] <<" " <<coord[1] <<" " <<coord[2] <<" " <<coord[3] <<" " <<coord[4] <<"\n";
												//}
											}
											//rel_err1 = loc_exp_err2*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim);
											//rel_err2 = loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim*loc_sim);
											//rel_err3 = loc_sim_err2*loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim*loc_sim*loc_sim);
											//rel_err = TMath::Sqrt(rel_err1+rel_err2+rel_err3)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp/loc_sim);
											//rel_err = TMath::Sqrt((loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim*loc_sim) ) + (loc_exp_err2*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim) ) + (loc_sim_err2*loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim*loc_sim*loc_sim) ))/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp/loc_sim);
											//std::cout<<"relative Error for dist " <<dist <<": 1 " <<rel_err1 <<" 2 " <<rel_err2 <<" 3 " <<rel_err3 <<" rel "<<rel_err <<" sf:" <<loc_exp/loc_sim <<" hole val:" <<_sim_holes_5d[i][j]->GetBinContent(bin) <<"\n";

											//sf_rel_err1 = loc_exp_err2/(loc_sim*loc_sim);
											//sf_rel_err2 = loc_sim_err2*loc_exp*loc_exp/(loc_sim*loc_sim*loc_sim*loc_sim);
											//sf_rel_err = TMath::Sqrt(sf_rel_err1+sf_rel_err2)/(loc_exp/loc_sim);
										}
									}
									if((_exp_data_5d_bop[i][j]->GetNbins() > 0 || (dubbin && _exp_data_5d_bop2[i][j]->GetNbins() > 0))){
										rel_err1 = loc_exp_err2*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim);
										rel_err2 = loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim*loc_sim);
										rel_err3 = loc_sim_err2*loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim*loc_sim*loc_sim);
										rel_err = TMath::Sqrt(rel_err1+rel_err2+rel_err3)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp/loc_sim);
										//rel_err = TMath::Sqrt((loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim*loc_sim) ) + (loc_exp_err2*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim) ) + (loc_sim_err2*loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim*loc_sim*loc_sim) ))/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp/loc_sim);
										//std::cout<<"relative Error for dist " <<dist <<": 1 " <<rel_err1 <<" 2 " <<rel_err2 <<" 3 " <<rel_err3 <<" rel "<<rel_err <<" sf:" <<loc_exp/loc_sim <<" hole val:" <<_sim_holes_5d[i][j]->GetBinContent(bin) <<"\n";

										sf_rel_err1 = loc_exp_err2/(loc_sim*loc_sim);
										sf_rel_err2 = loc_sim_err2*loc_exp*loc_exp/(loc_sim*loc_sim*loc_sim*loc_sim);
										sf_rel_err = TMath::Sqrt(sf_rel_err1+sf_rel_err2)/(loc_exp/loc_sim);
									}
									//loc_exp = fun::nSparseIntegral(_exp_data_5d[i][j]);
									//loc_sim = fun::nSparseIntegral(_sim_data_5d[i][j]);
									//if(loc_exp>0.0 && loc_sim>0.0 && (rel_err <= 0.5 || dist>=max_dist_)){
									if(loc_exp>0.0 && loc_sim>0.0 && (sf_rel_err <= 0.2 || dist>=max_dist_)){
										//std::cout<<"dist: " <<dist <<"sf: " <<loc_exp/loc_sim <<" exp: " <<loc_exp <<"  sim: " <<loc_sim <<" exp_error: " <<loc_exp_err2 <<" sim_error: " <<loc_sim_err2 <<" sim holes: " <<_sim_holes_5d[i][j]->GetBinContent(_sim_holes_5d[i][j]->GetBin(bin)) <<" sim_holes err2: " <<_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin)) <<" for bin:" <<bin[0] <<" " <<bin[1] <<" " <<bin[2] <<" " <<bin[3] <<" " <<bin[4] <<"\n";
										//std::cout<<"Setting exp\n";
										_scale_exp_5d[i][j]->SetBinContent(bin,loc_exp);
										//std::cout<<"Setting Sim\n";
										_scale_sim_5d[i][j]->SetBinContent(bin,loc_sim);
										look_further_all=false;
										//std::cout<<"\t\tA Test: bin|" <<_scale_exp_7d->GetBinContent(bin) <<" vs. what should be there|" <<fun::Sparse_Integral(scale_exp) <<" after dist:" <<dist <<"\n";
										//std::cout<<"Setting exp bin error\n";
										dist_dist[0][dist-1]++;
										_hole_radius_hist[i][j]->Fill(dist);
										_scale_exp_5d[i][j]->SetBinError2(_scale_exp_5d[i][j]->GetBin(bin),loc_exp_err2);
										//std::cout<<"Setting sim bin error\n";
										_scale_sim_5d[i][j]->SetBinError2(_scale_sim_5d[i][j]->GetBin(bin),loc_sim_err2);
										
									}
								}
							}
							if(dist==max_dist_){
								look_further_all = false;
								look_further_pos = false;
								look_further_neg = false;
							}
						}
					}
				}
			}
			for(int n=0; n<5; n++){
				if(flags_.Flags::Helicity()){
					_exp_data_5d_pos[i][j]->GetAxis(n)->SetRange();
					_exp_data_5d_neg[i][j]->GetAxis(n)->SetRange();
				}else{
					_exp_data_5d[i][j]->GetAxis(n)->SetRange();
				}
				_sim_data_5d[i][j]->GetAxis(n)->SetRange();
			}
			wq2_bin++;
			std::cout<<"\nW:" <<Histogram::W_low(i) <<"-" <<Histogram::W_top(i) <<" Q2:" <<Histogram::Q2_low(j) <<"-" <<Histogram::Q2_top(j) <<" " <<wq2_bin <<"/" <<29*5 <<"\n";
			std::cout<<"\tThrown Bins: " <<_thrown_5d[i][j]->GetNbins() <<" Recon Bins: " <<_sim_data_5d[i][j]->GetNbins() <<" Sim Holes Bins: " <<_sim_holes_5d[i][j]->GetNbins() <<"\n";
			std::cout<<"\tShould have: " <<_thrown_5d[i][j]->GetNbins()-_sim_data_5d[i][j]->GetNbins() <<" bins localized\n";
			std::cout<<"\tScale Sim nbins = " <<_scale_sim_5d_pos[i][j]->GetNbins() <<" Scale exp nbins = " <<_scale_exp_5d_pos[i][j]->GetNbins() <<"\n"; 
			std::cout<<"\tScale Sim nbins = " <<_scale_sim_5d_neg[i][j]->GetNbins() <<" Scale exp nbins = " <<_scale_exp_5d_neg[i][j]->GetNbins() <<"\n"; 
			double exp_nonlocal = 0.0;
			double sim_nonlocal = 0.0;
			for( long v_idx = 0; v_idx< _exp_data_5d[i][j]->GetNbins(); v_idx++){
				Double_t v_nonlocal = _exp_data_5d[i][j]->GetBinContent(v_idx,coord);
				if(v_nonlocal >0.0){
					Double_t w_nonlocal = _sim_data_5d[i][j]->GetBinContent(coord);
					if(w_nonlocal > 0.0){
						exp_nonlocal += v_nonlocal;
						sim_nonlocal += w_nonlocal;
					}
				}
			}
			std::cout<<"\tloc exp: " <<exp_nonlocal <<" loc sim: " <<sim_nonlocal <<"\n";
			std::cout<<"\tNonLocal Scale would be = " <<exp_nonlocal/sim_nonlocal <<"\n";
			for(int n=0; n<5; n++){
				if(flags_.Flags::Helicity()){
					_exp_data_5d_pos[i][j]->GetAxis(n)->SetRange(1,_n_bins_5d[n]);
					_exp_data_5d_neg[i][j]->GetAxis(n)->SetRange(1,_n_bins_5d[n]);
				}else{
					_exp_data_5d[i][j]->GetAxis(n)->SetRange(1,_n_bins_5d[n]);
				}
				_sim_data_5d[i][j]->GetAxis(n)->SetRange(1,_n_bins_5d[n]);
			}
			double exp_nonlocal2 = 0.0;
			double sim_nonlocal2 = 0.0;
			for( long v_idx = 0; v_idx< _exp_data_5d[i][j]->GetNbins(); v_idx++){
				Double_t v_nonlocal2 = _exp_data_5d[i][j]->GetBinContent(v_idx,coord);
				if(v_nonlocal2 >0.0){
					Double_t w_nonlocal2 = _sim_data_5d[i][j]->GetBinContent(coord);
					if(w_nonlocal2 > 0.0){
						exp_nonlocal2 += v_nonlocal2;
						sim_nonlocal2 += w_nonlocal2;
					}
				}
			}
			std::cout<<"\tloc exp: " <<exp_nonlocal2 <<" loc sim: " <<sim_nonlocal2 <<"\n";
			std::cout<<"\tNonLocal Scale would be = " <<exp_nonlocal2/sim_nonlocal2 <<"\n";
			
			for(int k=0; k<3; k++){
				for(int l=0; l<20; l++){
					std::cout<<dist_dist[k][l] <<" ";
				}
				std::cout<<"\n";
			}
			if(flags_.Flags::Helicity()){
				//_scale_exp_7d_pos->Divide(_acceptance_7d);
				_scale_5d_pos[i][j]=(THnSparseD*)_scale_exp_5d_pos[i][j]->Clone();
				_scale_5d_pos[i][j]->Divide(_scale_sim_5d_pos[i][j]);
				std::cout<<"Filling _sf_hist_pos\n";
				for (Long64_t l = 0; l < _scale_5d_pos[i][j]->GetNbins(); ++l) {
					_sf_hist_pos[i][j]->Fill(_scale_5d_pos[i][j]->GetBinContent(l));
				}
				_N_holes_5d_pos[i][j]=(THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_N_holes_5d_pos[i][j]->Multiply(_scale_5d_pos[i][j]);
				sprintf(hname,"Localized_Holes_W:%.3f-%.3f_Q2:%.2f-%.2f_pos",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_N_holes_5d_pos[i][j]->SetNameTitle(hname,hname);
				_N_holes_5d_pos[i][j]->Write();
				_N_holes_fifty_5d_pos[i][j] = (THnSparseD*)_N_holes_5d_pos[i][j]->Clone();
				
				for(long m=0; m<_N_holes_fifty_5d_pos[i][j]->GetNbins(); m++){
					_N_holes_fifty_5d_pos[i][j]->SetBinError(m,_N_holes_fifty_5d_pos[i][j]->GetBinContent(m)/2.0);
					if(_N_holes_fifty_5d_pos[i][j]->GetBinContent(m)>0.0){
						_hole_err_hist_pos[i][j]->Fill(_N_holes_5d_pos[i][j]->GetBinError(m)/_N_holes_5d_pos[i][j]->GetBinContent(m));
					}
				}
				sprintf(hname,"Localized_Holes_50_W:%.3f-%.3f_Q2:%.2f-%.2f_pos",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_N_holes_fifty_5d_pos[i][j]->SetNameTitle(hname,hname);
				_N_holes_fifty_5d_pos[i][j]->Write();
				//_scale_exp_7d_neg->Divide(_acceptance_7d);
				_scale_5d_neg[i][j]=(THnSparseD*)_scale_exp_5d_neg[i][j]->Clone();
				_scale_5d_neg[i][j]->Divide(_scale_sim_5d_neg[i][j]);
				std::cout<<"Filling _sf_hist_neg\n";
				for (Long64_t l = 0; l < _scale_5d_neg[i][j]->GetNbins(); ++l) {
					_sf_hist_neg[i][j]->Fill(_scale_5d_neg[i][j]->GetBinContent(l));
				}
				_N_holes_5d_neg[i][j]=(THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_N_holes_5d_neg[i][j]->Multiply(_scale_5d_neg[i][j]);
				sprintf(hname,"Localized_Holes_W:%.3f-%.3f_Q2:%.2f-%.2f_neg",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_N_holes_5d_neg[i][j]->SetNameTitle(hname,hname);
				_N_holes_5d_neg[i][j]->Write();
				_N_holes_fifty_5d_neg[i][j] = (THnSparseD*)_N_holes_5d_neg[i][j]->Clone();
				for(long m=0; m<_N_holes_fifty_5d_neg[i][j]->GetNbins(); m++){
					_N_holes_fifty_5d_neg[i][j]->SetBinError(m,_N_holes_fifty_5d_neg[i][j]->GetBinContent(m)/2.0);
					if(_N_holes_fifty_5d_neg[i][j]->GetBinContent(m)>0.0){
						_hole_err_hist_neg[i][j]->Fill(_N_holes_5d_neg[i][j]->GetBinError(m)/_N_holes_5d_neg[i][j]->GetBinContent(m));
					}
				}
				sprintf(hname,"Localized_Holes_50_W:%.3f-%.3f_Q2:%.2f-%.2f_neg",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_N_holes_fifty_5d_neg[i][j]->SetNameTitle(hname,hname);
				_N_holes_fifty_5d_neg[i][j]->Write();
				_sf_hist_pos[i][j]->Write();
				_sf_hist_neg[i][j]->Write();
				_hole_err_hist_pos[i][j]->Write();
				_hole_err_hist_neg[i][j]->Write();
				_hole_radius_hist_pos[i][j]->Write();
				_hole_radius_hist_neg[i][j]->Write();
			}else{
				//std::cout<<"Dividing exp by sim scales\n";
				//_scale_exp_7d->Divide(_acceptance_7d);//They're both not acceptance corrected, so adding this would be pointless
				//_scale_7d=(THnSparseD*)_scale_sim_7d->Clone();
				//_scale_7d->Divide(_scale_exp_7d);
				_scale_5d[i][j]=(THnSparseD*)_scale_exp_5d[i][j]->Clone();
				_scale_5d[i][j]->Divide(_scale_sim_5d[i][j]);
				std::cout<<"Filling _sf_hist\n";
				for (Long64_t l = 0; l < _scale_5d[i][j]->GetNbins(); ++l) {
					_sf_hist[i][j]->Fill(_scale_5d[i][j]->GetBinContent(l));
				}
				_N_holes_5d[i][j]=(THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_N_holes_5d[i][j]->Multiply(_scale_5d[i][j]);
				sprintf(hname,"Localized_Holes_W:%.3f-%.3f_Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_N_holes_5d[i][j]->SetNameTitle(hname,hname);
				_N_holes_5d[i][j]->Write();
				_N_holes_fifty_5d[i][j] = (THnSparseD*)_N_holes_5d[i][j]->Clone();
				Int_t* coord_50 = new Int_t[_sim_holes_5d[i][j]->GetNbins()];
				for(long m=0; m<_N_holes_fifty_5d[i][j]->GetNbins(); m++){
					_N_holes_fifty_5d[i][j]->SetBinError(m,_N_holes_fifty_5d[i][j]->GetBinContent(m)/2.0);
					if(_N_holes_fifty_5d[i][j]->GetBinContent(m,coord_50)>0.0){
						_hole_err_hist[i][j]->Fill(_N_holes_5d[i][j]->GetBinError(m)/_N_holes_5d[i][j]->GetBinContent(m));
						/*if(coord_50[0]==17 && coord_50[1]==7 && coord_50[2]==8 && coord_50[3]==5 && coord_50[4]==10){
							std::cout<<"Fifty Holes Content: " <<_N_holes_fifty_5d[i][j]->GetBinContent(m) <<" and Error: " <<_N_holes_fifty_5d[i][j]->GetBinError(m) <<" Old Error: " <<_N_holes_5d[i][j]->GetBinError(m) <<"\n";
						}*/
					}
				}
				sprintf(hname,"Localized_Holes_50_W:%.3f-%.3f_Q2:%.2f-%.2f",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_N_holes_fifty_5d[i][j]->SetNameTitle(hname,hname);
				_N_holes_fifty_5d[i][j]->Write();
				_sf_hist[i][j]->Write();
				_hole_err_hist[i][j]->Write();
				_hole_radius_hist[i][j]->Write();
			}
		}
	}
	std::cout<<"Finished Making the Localized holes\n";
	//Assign Errors 
	_RootOutputFile->Close();
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

void Histogram::Make_Acceptance_Rel_Error(){
	char hname[100];
	for(int a=0; a<_W_nbins_; a++){
		for(int b=0; b<_Q2_nbins_; b++){
			//Int_t* coord = new Int_t[_acceptance_5d[i][j]->GetNbins()];
			_acceptance_rel_err_5d[a][b] = (THnSparse*) _acceptance_5d[a][b]->Clone();
			sprintf(hname,"Acceptance Relative Error %s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(a),Histogram::W_top(a),Histogram::Q2_low(b),Histogram::Q2_top(b));
			_acceptance_rel_err_5d[a][b]->SetNameTitle(hname,hname);
			for(Long64_t i =0; i<_acceptance_5d[a][b]->GetNbins(); i++ ){
				if(_acceptance_5d[a][b]->GetBinContent(i)>0.0){
					_acceptance_rel_err_5d[a][b]->SetBinContent(i,_acceptance_5d[a][b]->GetBinError(i)/_acceptance_5d[a][b]->GetBinContent(i));
				}
			}
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