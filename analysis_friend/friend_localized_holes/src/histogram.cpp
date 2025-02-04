#include "histogram.hpp"


Histogram::Histogram(TFile* exp_tree_, TFile* sim_tree_, Flags flags_){
    Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,flags_);
    //Histogram::Rad_Corr();
    //Histogram::Sparse_7to5(flags_);
    //Histogram::Single_Diff(flags_);
	//Histogram::Localized_Holes(flags_,flags_.Flags::Min_Local_Dist(),24);
	Histogram::Localized_Holes(flags_,23,24);
}


void Histogram::Extract_5d_Histograms(TFile *exp_tree_, TFile *sim_tree_, Flags flags_){
    std::cout<<"Extract 5d Histograms\n";
	char hname[100];
	char hname2[300];
	sprintf(hname,"Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
	std::cout<<"Getting Thrown THnSparse " <<hname <<"\n";

    for(int i=0; i<_W_nbins_; i++){
        for(int j=0; j<_Q2_nbins_; j++){
            //sprintf(hname,"Thrown_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            sprintf(hname,"Thrown_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			std::cout<<"reading " <<hname <<"\n";
			//sprintf(hname2,"Thrown_2#pi_off_proton_#Delta^{++}_W:1.400-1.425_Q2:2.00-2.40;1-2#pi_off_proton_#Delta^{0}_mall_MM_Delta0_W:2.100-2.125_Q1:4.20-5.00/%s",hname);
			_thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
			//_thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname2);
            //_thrown_no_rad_5d[i][j] = (THnSparseD *)nr_sim_tree_->Get(hname);
            //sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			sprintf(hname,"%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			std::cout<<"reading exp and sim " <<hname <<"\n";
			_sim_data_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            _exp_data_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
			std::cout<<"extracted exp and sim " <<hname <<"\n";
			std::cout<<"Making Acceptance\n";
			std::cout<<"\tCloning Sim\n";
            _acceptance_5d[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
			std::cout<<"\tDividing by thrown\n";
            _acceptance_5d[i][j]->Divide(_thrown_5d[i][j]);
			sprintf(hname,"acceptance_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			_acceptance_5d[i][j]->SetNameTitle(hname,hname);
			std::cout<<"Entering acceptance loop for AREC cuts\n";
			for(long k=0; k<_acceptance_5d[i][j]->GetNbins(); k++){
				if(_acceptance_5d[i][j]->GetBinContent(k) >0.0){
					if(_acceptance_5d[i][j]->GetBinError(k)/_acceptance_5d[i][j]->GetBinContent(k) > _Acceptance_Rel_Error_Max_[flags_.Flags::Acc_Rel_Error_Cut()]){
						_acceptance_5d[i][j]->SetBinContent(k,0.0);
						_acceptance_5d[i][j]->SetBinError(k,0.0);
					}
				}
			}
			_exp_data_5d_acc_corr[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
			_exp_data_5d_acc_corr[i][j]->Divide(_acceptance_5d[i][j]);
			sprintf(hname,"exp_acc_corr_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			_exp_data_5d_acc_corr[i][j]->SetNameTitle(hname,hname);
			_sim_data_5d_acc_corr[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
			_sim_data_5d_acc_corr[i][j]->Divide(_acceptance_5d[i][j]);
			sprintf(hname,"sim_acc_corr_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			_sim_data_5d_acc_corr[i][j]->SetNameTitle(hname,hname);
			_sim_holes_5d[i][j] = (THnSparseD*)_thrown_5d[i][j]->Clone();
			_sim_holes_5d[i][j]->Add(_sim_data_5d_acc_corr[i][j],-1.0);
			
			std::cout<<"Checking for helicity\n";
			if(flags_.Flags::Helicity()){
				//std::cout<<"\t1\n";
				//sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				sprintf(hname,"%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				std::cout<<"reading " <<hname <<"\n";
				_exp_data_5d_pos[i][j] = (THnSparseD*)exp_tree_->Get(hname);
				std::cout<<"How many bins are here? " <<_exp_data_5d_pos[i][j]->GetNbins() <<"\n";
				std::cout<<"cloning " <<hname <<"\n";
				_exp_data_5d_acc_corr_pos[i][j] = (THnSparseD*)_exp_data_5d_pos[i][j]->Clone();
				sprintf(hname,"Acc_Corr_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_exp_data_5d_acc_corr_pos[i][j]->SetNameTitle(hname,hname);
				std::cout<<"dividing by acceptance \n";
				_exp_data_5d_acc_corr_pos[i][j]->Divide(_acceptance_5d[i][j]);
				_sim_holes_5d_pos[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				
				
				//sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				sprintf(hname,"%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				std::cout<<"reading " <<hname <<"\n";
				_exp_data_5d_neg[i][j] = (THnSparseD*)exp_tree_->Get(hname);
				std::cout<<"cloning " <<hname <<"\n";
				_exp_data_5d_acc_corr_neg[i][j] = (THnSparseD*)_exp_data_5d_neg[i][j]->Clone();
				sprintf(hname,"Acc_Corr_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_exp_data_5d_acc_corr_neg[i][j]->SetNameTitle(hname,hname);
				std::cout<<"dividing by acceptance \n";
				_exp_data_5d_acc_corr_neg[i][j]->Divide(_acceptance_5d[i][j]);
				_sim_holes_5d_neg[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				
			}
        }
    }
    for(int k = 0; k<_thrown_5d[0][0]->GetNdimensions(); k++){
		_n_bins_5d.push_back(_thrown_5d[0][0]->GetAxis(k)->GetNbins());
	}
}


void Histogram::Localized_Holes(Flags flags_, int min_dist_, int max_dist_){
	if(!flags_.Flags::Plot_Localized_Holes()){
		return;
	}
	char aname[300];
	string bop = flags_.Flags::Output_File().c_str();
	bop.erase(bop.length()-5);
	sprintf(aname,"%s_glob.root",bop);
	//_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	_RootOutputFile = new TFile(aname,"RECREATE");
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
	long dist_dist[3][24] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}; 
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

	double int_exp_sync_to_sim;
	double int_sim;
	
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

	std::cout<<"Making acceptance output histograms\n";
	Double_t Q2_bin[6] = {_Q2_bins_[0],_Q2_bins_[1],_Q2_bins_[2],_Q2_bins_[3],_Q2_bins_[4],_Q2_bins_[5]};
	sprintf(hname,"Global Scale Factor");
	_global_sf = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
	sprintf(hname,"Global Scale Factor2");
	_global_sf2 = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
	sprintf(hname,"Fraction of Zero Acceptance Bins where Exp Exists");
	_zero_acc_w_exp_data = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
	sprintf(hname,"Fraction of Lost Exp Events due to Acceptance Limitations");
	_lost_exp_events = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
	sprintf(hname,"Fraction of Hole Bins to Exp Bins");
	_frac_bins_holes_to_exp = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
	sprintf(hname,"Fraction Hole Events to Exp Acceptance Corr Events");
	_frac_events_holes_to_exp = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);

	if(flags_.Flags::Helicity()){
		sprintf(hname,"Global Scale Factor Pos");
		_global_sf_pos = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
		sprintf(hname,"Global Scale Factor Neg");
		_global_sf_neg = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
		sprintf(hname,"Fraction of Zero Acceptance Bins where Pos Exp Exists");
		_zero_acc_w_exp_data_pos = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
		sprintf(hname,"Fraction of Zero Acceptance Bins where Neg Exp Exists");
		_zero_acc_w_exp_data_neg = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
		sprintf(hname,"Fraction of Lost Pos Exp Events due to Acceptance Limitations");
		_lost_exp_events_pos = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
		sprintf(hname,"Fraction of Lost Neg Exp Events due to Acceptance Limitations");
		_lost_exp_events_neg = new TH2D(hname,hname,29,_W_min_,_W_max_,5,Q2_bin);
	}
	std::cout<<"Beginning Acceptance output Histogram loop\n";
	for(int i=0; i<_W_nbins_; i++){
		for(int j=0; j<_Q2_nbins_; j++){
			if(flags_.Flags::Helicity()){
				sprintf(hname,"pos scale factor distribution W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sf_hist_pos[i][j] = new TH1D(hname,hname,200,1.0,500);
				//_sf_hist_pos[i][j] = new TH1D(hname,hname,200,0.0,100.0);
				sprintf(hname,"neg scale factor distribution W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sf_hist_neg[i][j] = new TH1D(hname,hname,200,1.0,500);
				//_sf_hist_neg[i][j] = new TH1D(hname,hname,200,0.0,100.0);
				sprintf(hname,"pos localized holes relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_err_hist_pos[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				sprintf(hname,"neg localized holes relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_err_hist_neg[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				sprintf(hname,"pos localized radius W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_radius_hist_pos[i][j] = new TH1D(hname,hname,24,0.5,24.5);
				sprintf(hname,"neg localized radius W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_radius_hist_neg[i][j] = new TH1D(hname,hname,24,0.5,24.5);
				sprintf(hname,"scale factor pos relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sf_rel_err_hist_pos[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				sprintf(hname,"scale factor neg relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sf_rel_err_hist_neg[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				float num_zeros_pos = 0.0;
				float num_total_pos = 0.0;
				double num_exp_events_pos = 0.0;
				double num_lost_exp_events_pos = 0.0;
				float num_zeros_neg = 0.0;
				float num_total_neg = 0.0;
				double num_exp_events_neg = 0.0;
				double num_lost_exp_events_neg = 0.0;
				//std::cout<<"made histograms for W:" <<W_mid(i) <<" Q2:" <<Q2_mid(j) <<"\n";
				while(cart.GetNextCombination()){
					for(int q=0; q<5; q++){
						binning[q] = cart[q]+1;
					}
					if(_exp_data_5d_pos[i][j]->GetBinContent(binning) > 0.0){
						num_total_pos += 1.0;
						num_exp_events_pos += _exp_data_5d_pos[i][j]->GetBinContent(binning);
						//_acceptance_dist_pos[i][j]->Fill(_acceptance_5d[i][j]->GetBinContent(binning));
						if(_acceptance_5d[i][j]->GetBinContent(binning) == 0.0){
							num_zeros_pos += 1.0;
							num_lost_exp_events_pos += _exp_data_5d_pos[i][j]->GetBinContent(binning);
						}
					}//else{
					//	_acceptance_dist[i][j]->Fill(_acceptance_5d[i][j]->GetBinContent(binning));
					//}
					if(_exp_data_5d_neg[i][j]->GetBinContent(binning) > 0.0){
						num_total_neg += 1.0;
						num_exp_events_neg += _exp_data_5d_neg[i][j]->GetBinContent(binning);
						//_acceptance_dist_neg[i][j]->Fill(_acceptance_5d[i][j]->GetBinContent(binning));
						if(_acceptance_5d[i][j]->GetBinContent(binning) == 0.0){
							num_zeros_neg += 1.0;
							num_lost_exp_events_neg += _exp_data_5d_neg[i][j]->GetBinContent(binning);
						}
					}
				}
				_lost_exp_events_pos->Fill(W_mid(i),Q2_mid(j),num_lost_exp_events_pos/num_exp_events_pos);
				_zero_acc_w_exp_data_pos->Fill(W_mid(i),Q2_mid(j),num_zeros_pos/num_total_pos);
				_lost_exp_events_neg->Fill(W_mid(i),Q2_mid(j),num_lost_exp_events_neg/num_exp_events_neg);
				_zero_acc_w_exp_data_neg->Fill(W_mid(i),Q2_mid(j),num_zeros_neg/num_total_neg);
				sprintf(hname,"pos global holes relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_global_hole_err_hist_pos[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				sprintf(hname,"neg global holes relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_global_hole_err_hist_neg[i][j] = new TH1D(hname,hname,400,0.0,1.6);
			}else{
				sprintf(hname,"scale factor distribution W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sf_hist[i][j] = new TH1D(hname,hname,120,1.0,240);

				//sprintf(hname,"old scale factor distribution W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				//_sf_hist_old[i][j] = new TH1D(hname,hname,100,0.0,240);
				//_sf_hist[i][j] = new TH1D(hname,hname,200,0.0,100.0);
				sprintf(hname,"Acceptance Distribution in Exp PS W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				//_acceptance_dist[i][j] = new TH1D(hname,hname,200,0.0,0.2);
				_acceptance_dist[i][j] = new TH1D(hname,hname,300,0.0,0.2);
				sprintf(hname,"Acceptance Distribution NonZero W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				//_acceptance_dist[i][j] = new TH1D(hname,hname,200,0.0,0.2);
				_acceptance_dist2[i][j] = new TH1D(hname,hname,300,0.0,0.2);
				sprintf(hname,"Acceptance Distribution Full W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				//_acceptance_dist[i][j] = new TH1D(hname,hname,200,0.0,0.2);
				_acceptance_dist3[i][j] = new TH1D(hname,hname,300,0.0,0.2);
				sprintf(hname,"Acceptance Distribution non-exp phase space W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				//_acceptance_dist[i][j] = new TH1D(hname,hname,200,0.0,0.2);
				_nacceptance_dist[i][j] = new TH1D(hname,hname,200,0.0,0.002);
				sprintf(hname,"Relative Acceptance Distribution W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_rel_acceptance_dist[i][j] = new TH1D(hname,hname,100,0.0,1.2);
				sprintf(hname,"localized holes relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_err_hist[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				//sprintf(hname,"old localized holes relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				//_hole_err_hist_old[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				sprintf(hname,"global holes relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_global_hole_err_hist[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				//sprintf(hname,"old global holes relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				//_global_hole_err_hist_old[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				sprintf(hname,"scale factor relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sf_rel_err_hist[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				//sprintf(hname,"old scale factor relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				//_sf_rel_err_hist_old[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				sprintf(hname,"localized radius W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_hole_radius_hist[i][j] = new TH1D(hname,hname,24,0.5,24.5);
				sprintf(hname,"exp acc corr relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_exp_acc_corr_rel_err[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				sprintf(hname,"sim acc corr relative error W:[%.3f GeV,%.3f GeV) ^{}Q^{2}:[%.2f ^{}GeV^{2}, %.2f ^{}GeV^{2})",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				_sim_acc_corr_rel_err[i][j] = new TH1D(hname,hname,400,0.0,1.6);
				float num_zeros = 0.0;
				float num_total = 0.0;
				double num_exp_events = 0.0;
				double num_lost_exp_events = 0.0;
				//std::cout<<"made histograms for W:" <<W_mid(i) <<" Q2:" <<Q2_mid(j) <<"\n";
				while(cart.GetNextCombination()){
					for(int q=0; q<5; q++){
						binning[q] = cart[q]+1;
					}
					if(_exp_data_5d[i][j]->GetBinContent(binning) > 0.0){
						num_total += 1.0;
						num_exp_events += _exp_data_5d[i][j]->GetBinContent(binning);
						_acceptance_dist[i][j]->Fill(_acceptance_5d[i][j]->GetBinContent(binning));
						if(_acceptance_5d[i][j]->GetBinContent(binning) == 0.0){
							num_zeros += 1.0;
							num_lost_exp_events += _exp_data_5d[i][j]->GetBinContent(binning);
						}
					}//else{
					//	_acceptance_dist[i][j]->Fill(_acceptance_5d[i][j]->GetBinContent(binning));
					//}
				}
				_lost_exp_events->Fill(W_mid(i),Q2_mid(j),num_lost_exp_events/num_exp_events);
				//std::cout<<"Filling zero acc histogram\n";
				_zero_acc_w_exp_data->Fill(W_mid(i),Q2_mid(j),num_zeros/num_total);
				double min_acc = 1.0;
				double max_acc = 0.0;
				//std::cout<<"filling rel error acceptance histograms\n";
				for (Long64_t l = 0; l < _acceptance_5d[i][j]->GetNbins(); ++l) {
					_acceptance_dist3[i][j]->Fill(_acceptance_5d[i][j]->GetBinContent(l));
					if(_acceptance_5d[i][j]->GetBinContent(l) >0.0){
						_rel_acceptance_dist[i][j]->Fill(_acceptance_5d[i][j]->GetBinError(l)/_acceptance_5d[i][j]->GetBinContent(l));
						_acceptance_dist2[i][j]->Fill(_acceptance_5d[i][j]->GetBinContent(l));
						if(_acceptance_5d[i][j]->GetBinContent(l)<min_acc ){
							min_acc = _acceptance_5d[i][j]->GetBinContent(l);
						}
						if(_acceptance_5d[i][j]->GetBinContent(l)>max_acc ){
							max_acc = _acceptance_5d[i][j]->GetBinContent(l);
						}
					}
				}
				//std::cout<<"\tW:" <<W_mid(i) <<" Q2:" <<Q2_mid(j) <<"\tmin_acc:" <<min_acc <<"\tmax_acc:" <<max_acc <<"\n"; 
				_acceptance_dist[i][j]->SetXTitle("Acceptance Value for Individual Bin");
				_acceptance_dist[i][j]->SetYTitle("Number of Bins");
				_acceptance_dist[i][j]->Write();
				_acceptance_dist2[i][j]->SetXTitle("Acceptance Value for Individual Bin");
				_acceptance_dist2[i][j]->SetYTitle("Number of Bins");
				_acceptance_dist2[i][j]->Write();
				_acceptance_dist3[i][j]->SetXTitle("Acceptance Value for Individual Bin");
				_acceptance_dist3[i][j]->SetYTitle("Number of Bins");
				_acceptance_dist3[i][j]->Write();
				_rel_acceptance_dist[i][j]->SetXTitle("Relative Acceptance Error for Individual Bin");
				_rel_acceptance_dist[i][j]->SetYTitle("Number of Bins");
				_rel_acceptance_dist[i][j]->Write();
			}
		}
	}
	std::cout<<"Made it out of acceptance output loop\n";
	if(flags_.Flags::Helicity()){
		_zero_acc_w_exp_data_pos->SetXTitle("W (GeV)");
		_zero_acc_w_exp_data_pos->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_zero_acc_w_exp_data_pos->Write();
		_zero_acc_w_exp_data_neg->SetXTitle("W (GeV)");
		_zero_acc_w_exp_data_neg->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_zero_acc_w_exp_data_neg->Write();
		_lost_exp_events_pos->SetXTitle("W (GeV)");
		_lost_exp_events_pos->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_lost_exp_events_pos->Write();
		_lost_exp_events_neg->SetXTitle("W (GeV)");
		_lost_exp_events_neg->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_lost_exp_events_neg->Write();
	}else{
		_zero_acc_w_exp_data->SetXTitle("W (GeV)");
		_zero_acc_w_exp_data->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_zero_acc_w_exp_data->Write();
		_lost_exp_events->SetXTitle("W (GeV)");
		_lost_exp_events->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_lost_exp_events->Write();
	}
	

	std::cout<<"*****************The helicity flag is " <<flags_.Flags::Helicity() <<"\n";
	for(int i=0; i<_W_nbins_; i++){
	//for(int i=0; i<10; i++){
		for(int j=0; j<_Q2_nbins_; j++){
			for(int k=0; k<3; k++){
				for(int l=0; l<24; l++){
					dist_dist[k][l] = 0;
				}
			}
			total_ev[i][j] = _sim_holes_5d[i][j]->GetNbins();
			curr_ev[i][j] = 0;
			if(flags_.Flags::Helicity()){
				_scale_exp_5d_pos[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_exp_5d_pos[i][j]->Scale(0.0);
				_scale_exp_5d_pos[i][j]->Sumw2();
				_scale_exp_5d_neg[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_exp_5d_neg[i][j]->Scale(0.0);
				_scale_exp_5d_neg[i][j]->Sumw2();
				_scale_sim_5d_pos[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_sim_5d_pos[i][j]->Scale(0.0);
				_scale_sim_5d_pos[i][j]->Sumw2();
				_scale_sim_5d_neg[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_sim_5d_neg[i][j]->Scale(0.0);
				_scale_sim_5d_neg[i][j]->Sumw2();
				_scale_exp_5d_global_pos[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_exp_5d_global_pos[i][j]->Scale(0.0);
				_scale_exp_5d_global_pos[i][j]->Sumw2();
				_scale_exp_5d_global_neg[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_exp_5d_global_neg[i][j]->Scale(0.0);
				_scale_exp_5d_global_neg[i][j]->Sumw2();
				_scale_sim_5d_global_pos[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_sim_5d_global_pos[i][j]->Scale(0.0);
				_scale_sim_5d_global_pos[i][j]->Sumw2();
				_scale_sim_5d_global_neg[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_sim_5d_global_neg[i][j]->Scale(0.0);
				_scale_sim_5d_global_neg[i][j]->Sumw2();
			}else{
				_scale_exp_5d[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_exp_5d[i][j]->Scale(0.0);
				_scale_exp_5d[i][j]->Sumw2();
				_scale_sim_5d[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_sim_5d[i][j]->Scale(0.0);
				_scale_sim_5d[i][j]->Sumw2();
				_scale_exp_5d_global[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_exp_5d_global[i][j]->Scale(0.0);
				_scale_exp_5d_global[i][j]->Sumw2();
				_scale_sim_5d_global[i][j] = (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_scale_sim_5d_global[i][j]->Scale(0.0);
				_scale_sim_5d_global[i][j]->Sumw2();
			}
			//std::cout<<"made scaling THnSparse histograms for W:" <<i <<" Q2:" <<j <<"\n";
			

			

			
			Int_t* coord = new Int_t[_sim_holes_5d[i][j]->GetNbins()];
			//std::cout<<"made initial THnSparse coordinate thing\n";
			bool dubbin[5] = {false,false,false,false,false};
			while(cart.GetNextCombination()){
				//std::cout<<"started cart loop\n";
				dist = 0;
				for(int k = 0; k<_n_bins_5d.size(); k++){
					bin[k] = cart[k]+1;
				}
				if(flags_.Flags::Helicity()){
					//std::cout<<"checking sim holes values for helicity\n";
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
					//std::cout<<"checked sim holes pos:" <<look_further_pos <<" neg:" <<look_further_neg <<"\n";
				}else{
					if(_sim_holes_5d[i][j]->GetBinContent(bin)>0.0){// && bin[0]>=3 && bin[0]<=16 && bin[1]>=3 && bin[1]<=16){
						look_further_all = true;
					}else{
						look_further_all = false;
					}
				}
				//std::cout<<"initial setting look_further flags\n";
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
					//std::cout<<"second setting look_further flags\n";
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
							//std::cout<<"set integration bins\n";
							if(flags_.Flags::Helicity()){
								if(look_further_pos){
									//std::cout<<"looking further in positive helicity\t dist:" <<dist <<"\n";
									loc_exp_pos = 0.0;
									loc_sim_pos = 0.0;
									loc_exp_err2_pos = 0.0;
									loc_sim_err2_pos = 0.0;
									//if(_exp_data_5d_pos[i][j]->GetNbins() > 0){
									if(_sim_data_5d_acc_corr[i][j]->GetNbins() > 0){
										for( long v_idx = 0; v_idx< _sim_data_5d_acc_corr[i][j]->GetNbins(); v_idx++){
											int passes = 0;
											Double_t v_pos = _sim_data_5d_acc_corr[i][j]->GetBinContent(v_idx,coord);
											if(v_pos>0.0){
												bool pass = true;
												bool pass_alpha = dubbin[3];
												bool pass_phi = dubbin[4];
												bool pass_both = (dubbin[3] && dubbin[4]);
												for(int z = 0; z<5; z++){
													pass &= coord[z] >= bin_low[z];
													pass &= coord[z] <= bin_top[z];
													if(dubbin[3]){
														if(z==3){
															pass_alpha &= coord[z] >= bin_low2[z];
															pass_alpha &= coord[z] <= bin_top2[z];
														}else{
															pass_alpha &= coord[z] >= bin_low[z];
															pass_alpha &= coord[z] <= bin_top[z];
														}
														if(dubbin[4]){
															if(z==4){
																pass_phi &= coord[z] >= bin_low2[z];
																pass_phi &= coord[z] <= bin_top2[z];
															}else{
																pass_phi &= coord[z] >= bin_low[z];
																pass_phi &= coord[z] <= bin_top[z];
															}
															if(z==3 || z==4){
																pass_both &= coord[z] >= bin_low2[z];
																pass_both &= coord[z] <= bin_top2[z];
															}else{
																pass_both &= coord[z] >= bin_low[z];
																pass_both &= coord[z] <= bin_top[z];
															}
														}
													}else if(dubbin[4]){
														pass_alpha = false;
														pass_both = false;
														if(z==4){
															pass_phi &= coord[z] >= bin_low2[z];
															pass_phi &= coord[z] <= bin_top2[z];
														}else{
															pass_phi &= coord[z] >= bin_low[z];
															pass_phi &= coord[z] <= bin_top[z];
														}
													}
												}
												
												if(pass){
													passes ++;
													Double_t w_pos = _exp_data_5d_acc_corr_pos[i][j]->GetBinContent(coord);
													loc_exp_pos += w_pos;
													loc_sim_pos += v_pos; 
													loc_exp_err2_pos += _exp_data_5d_acc_corr_pos[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2_pos += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}else if(pass_alpha){
													passes ++;
													Double_t w_pos = _exp_data_5d_acc_corr_pos[i][j]->GetBinContent(coord);
													loc_exp_pos += w_pos;
													loc_sim_pos += v_pos; 
													loc_exp_err2_pos += _exp_data_5d_acc_corr_pos[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2_pos += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}else if(pass_phi){
													passes ++;
													Double_t w_pos = _exp_data_5d_acc_corr_pos[i][j]->GetBinContent(coord);
													loc_exp_pos += w_pos;
													loc_sim_pos += v_pos; 
													loc_exp_err2_pos += _exp_data_5d_acc_corr_pos[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2_pos += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}else if(pass_both){
													passes ++;
													Double_t w_pos = _exp_data_5d_acc_corr_pos[i][j]->GetBinContent(coord);
													loc_exp_pos += w_pos;
													loc_sim_pos += v_pos; 
													loc_exp_err2_pos += _exp_data_5d_acc_corr_pos[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2_pos += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}
												if(passes >1){
													std::cout<<"# Passes:" <<passes <<"\n";
													std::cout<<"pass:"<<pass <<" alpha:" <<pass_alpha <<" phi:" <<pass_phi <<" both:" <<pass_both <<"\n";
												}
											}
										}
										rel_err1_pos = loc_exp_err2_pos*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_pos*loc_sim_pos);
										rel_err2_pos = loc_exp_pos*loc_exp_pos*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim_pos*loc_sim_pos);
										rel_err3_pos = loc_sim_err2_pos*loc_exp_pos*loc_exp_pos*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_pos*loc_sim_pos*loc_sim_pos*loc_sim_pos);
										rel_err_pos = TMath::Sqrt(rel_err1_pos+rel_err2_pos+rel_err3_pos)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp_pos/loc_sim_pos);

										sf_rel_err1_pos = loc_exp_err2_pos/(loc_sim_pos*loc_sim_pos);
										sf_rel_err2_pos = loc_sim_err2_pos*loc_exp_pos*loc_exp_pos/(loc_sim_pos*loc_sim_pos*loc_sim_pos*loc_sim_pos);
										sf_rel_err_pos = TMath::Sqrt(sf_rel_err1_pos+sf_rel_err2_pos)/(loc_exp_pos/loc_sim_pos);
									}
									if(loc_exp_pos>0.0 && loc_sim_pos>0.0 && (sf_rel_err_pos <= 0.2 || dist>=max_dist_)){	
										_scale_exp_5d_pos[i][j]->SetBinContent(bin,loc_exp_pos);
										_scale_sim_5d_pos[i][j]->SetBinContent(bin,loc_sim_pos);
										look_further_pos=false;
										_scale_exp_5d_pos[i][j]->SetBinError2(_scale_exp_5d_pos[i][j]->GetBin(bin),loc_exp_err2_pos);
										_scale_sim_5d_pos[i][j]->SetBinError2(_scale_sim_5d_pos[i][j]->GetBin(bin),loc_sim_err2_pos);
										dist_dist[1][dist-1]++;
										_hole_radius_hist_pos[i][j]->Fill(dist);
										_sf_rel_err_hist_pos[i][j]->Fill(sf_rel_err_pos);
									}
								}
								if(look_further_neg){
									loc_exp_neg = 0.0;
									loc_sim_neg = 0.0;
									loc_exp_err2_neg = 0.0;
									loc_sim_err2_neg = 0.0;
									//if(_exp_data_5d_neg[i][j]->GetNbins() > 0){
									if(_sim_data_5d_acc_corr[i][j]->GetNbins() > 0){
										for( long v_idx = 0; v_idx< _sim_data_5d_acc_corr[i][j]->GetNbins(); v_idx++){
											int passes = 0;
											Double_t v_neg = _sim_data_5d_acc_corr[i][j]->GetBinContent(v_idx,coord);
											if(v_neg>0.0){
												bool pass = true;
												bool pass_alpha = dubbin[3];
												bool pass_phi = dubbin[4];
												bool pass_both = (dubbin[3] && dubbin[4]);
												for(int z = 0; z<5; z++){
													pass &= coord[z] >= bin_low[z];
													pass &= coord[z] <= bin_top[z];
													if(dubbin[3]){
														if(z==3){
															pass_alpha &= coord[z] >= bin_low2[z];
															pass_alpha &= coord[z] <= bin_top2[z];
														}else{
															pass_alpha &= coord[z] >= bin_low[z];
															pass_alpha &= coord[z] <= bin_top[z];
														}
														if(dubbin[4]){
															if(z==4){
																pass_phi &= coord[z] >= bin_low2[z];
																pass_phi &= coord[z] <= bin_top2[z];
															}else{
																pass_phi &= coord[z] >= bin_low[z];
																pass_phi &= coord[z] <= bin_top[z];
															}
															if(z==3 || z==4){
																pass_both &= coord[z] >= bin_low2[z];
																pass_both &= coord[z] <= bin_top2[z];
															}else{
																pass_both &= coord[z] >= bin_low[z];
																pass_both &= coord[z] <= bin_top[z];
															}
														}
													}else if(dubbin[4]){
														pass_alpha = false;
														pass_both = false;
														if(z==4){
															pass_phi &= coord[z] >= bin_low2[z];
															pass_phi &= coord[z] <= bin_top2[z];
														}else{
															pass_phi &= coord[z] >= bin_low[z];
															pass_phi &= coord[z] <= bin_top[z];
														}
													}
												}
												
												if(pass){
													passes ++;
													Double_t w_neg = _exp_data_5d_acc_corr_neg[i][j]->GetBinContent(coord);
													loc_exp_neg += w_neg;
													loc_sim_neg += v_neg; 
													loc_exp_err2_neg += _exp_data_5d_acc_corr_neg[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2_neg += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}else if(pass_alpha){
													passes ++;
													Double_t w_neg = _exp_data_5d_acc_corr_neg[i][j]->GetBinContent(coord);
													loc_exp_neg += w_neg;
													loc_sim_neg += v_neg; 
													loc_exp_err2_neg += _exp_data_5d_acc_corr_neg[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2_neg += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}else if(pass_phi){
													passes ++;
													Double_t w_neg = _exp_data_5d_acc_corr_neg[i][j]->GetBinContent(coord);
													loc_exp_neg += w_neg;
													loc_sim_neg += v_neg; 
													loc_exp_err2_neg += _exp_data_5d_acc_corr_neg[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2_neg += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}else if(pass_both){
													passes ++;
													Double_t w_neg = _exp_data_5d_acc_corr_neg[i][j]->GetBinContent(coord);
													loc_exp_neg += w_neg;
													loc_sim_neg += v_neg; 
													loc_exp_err2_neg += _exp_data_5d_acc_corr_neg[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2_neg += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}
												if(passes >1){
													std::cout<<"# Passes:" <<passes <<"\n";
													std::cout<<"pass:"<<pass <<" alpha:" <<pass_alpha <<" phi:" <<pass_phi <<" both:" <<pass_both <<"\n";
												}
											}
										}
										rel_err1_neg = loc_exp_err2_neg*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_neg*loc_sim_neg);
										rel_err2_neg = loc_exp_neg*loc_exp_neg*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim_neg*loc_sim_neg);
										rel_err3_neg = loc_sim_err2_neg*loc_exp_neg*loc_exp_neg*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim_neg*loc_sim_neg*loc_sim_neg*loc_sim_neg);
										rel_err_neg = TMath::Sqrt(rel_err1_neg+rel_err2_neg+rel_err3_neg)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp_neg/loc_sim_neg);

										sf_rel_err1_neg = loc_exp_err2_neg/(loc_sim_neg*loc_sim_neg);
										sf_rel_err2_neg = loc_sim_err2_neg*loc_exp_neg*loc_exp_neg/(loc_sim_neg*loc_sim_neg*loc_sim_neg*loc_sim_neg);
										sf_rel_err_neg = TMath::Sqrt(sf_rel_err1_neg+sf_rel_err2_neg)/(loc_exp_neg/loc_sim_neg);
									}
									if(loc_exp_neg>0.0 && loc_sim_neg>0.0 && (sf_rel_err_neg <= 0.2 || dist>=max_dist_)){	
										_scale_exp_5d_neg[i][j]->SetBinContent(bin,loc_exp_neg);
										_scale_sim_5d_neg[i][j]->SetBinContent(bin,loc_sim_neg);
										look_further_neg=false;
										_scale_exp_5d_neg[i][j]->SetBinError2(_scale_exp_5d_neg[i][j]->GetBin(bin),loc_exp_err2_neg);
										_scale_sim_5d_neg[i][j]->SetBinError2(_scale_sim_5d_neg[i][j]->GetBin(bin),loc_sim_err2_neg);
										dist_dist[2][dist-1]++;
										_hole_radius_hist_neg[i][j]->Fill(dist);
										_sf_rel_err_hist_neg[i][j]->Fill(sf_rel_err_neg);
									}
								}
							}else{
								if(look_further_all){
									loc_exp= 0.0;
									loc_sim = 0.0;
									loc_exp_err2= 0.0;
									loc_sim_err2 = 0.0;

									if(_sim_data_5d_acc_corr[i][j]->GetNbins() > 0){
										for( long v_idx = 0; v_idx< _sim_data_5d_acc_corr[i][j]->GetNbins(); v_idx++){
											int passes = 0;
											Double_t v = _sim_data_5d_acc_corr[i][j]->GetBinContent(v_idx,coord);
											if(v>0.0){
												bool pass = true;
												bool pass_alpha = dubbin[3];
												bool pass_phi = dubbin[4];
												bool pass_both = (dubbin[3] && dubbin[4]);
												for(int z = 0; z<5; z++){
													pass &= coord[z] >= bin_low[z];
													pass &= coord[z] <= bin_top[z];
													if(dubbin[3]){
														if(z==3){
															pass_alpha &= coord[z] >= bin_low2[z];
															pass_alpha &= coord[z] <= bin_top2[z];
														}else{
															pass_alpha &= coord[z] >= bin_low[z];
															pass_alpha &= coord[z] <= bin_top[z];
														}
														if(dubbin[4]){
															if(z==4){
																pass_phi &= coord[z] >= bin_low2[z];
																pass_phi &= coord[z] <= bin_top2[z];
															}else{
																pass_phi &= coord[z] >= bin_low[z];
																pass_phi &= coord[z] <= bin_top[z];
															}
															if(z==3 || z==4){
																pass_both &= coord[z] >= bin_low2[z];
																pass_both &= coord[z] <= bin_top2[z];
															}else{
																pass_both &= coord[z] >= bin_low[z];
																pass_both &= coord[z] <= bin_top[z];
															}
														}
													}else if(dubbin[4]){
														pass_alpha = false;
														pass_both = false;
														if(z==4){
															pass_phi &= coord[z] >= bin_low2[z];
															pass_phi &= coord[z] <= bin_top2[z];
														}else{
															pass_phi &= coord[z] >= bin_low[z];
															pass_phi &= coord[z] <= bin_top[z];
														}
													}
												}
												if(pass){
													passes ++;
													Double_t w = _exp_data_5d_acc_corr[i][j]->GetBinContent(coord);
													loc_exp += w;
													loc_sim += v; 
													loc_exp_err2 += _exp_data_5d_acc_corr[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2 += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}else if(pass_alpha){
													passes ++;
													Double_t w = _exp_data_5d_acc_corr[i][j]->GetBinContent(coord);
													loc_exp += w;
													loc_sim += v; 
													loc_exp_err2 += _exp_data_5d_acc_corr[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2 += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}else if(pass_phi){
													passes ++;
													Double_t w = _exp_data_5d_acc_corr[i][j]->GetBinContent(coord);
													loc_exp += w;
													loc_sim += v; 
													loc_exp_err2 += _exp_data_5d_acc_corr[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2 += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}else if(pass_both){
													passes ++;
													Double_t w = _exp_data_5d_acc_corr[i][j]->GetBinContent(coord);
													loc_exp += w;
													loc_sim += v; 
													loc_exp_err2 += _exp_data_5d_acc_corr[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));
													loc_sim_err2 += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);
												}
												if(passes >1){
													std::cout<<"# Passes:" <<passes <<"\n";
													std::cout<<"pass:"<<pass <<" alpha:" <<pass_alpha <<" phi:" <<pass_phi <<" both:" <<pass_both <<"\n";
												}
											}
										}
									}
									if(_sim_data_5d_acc_corr[i][j]->GetNbins() > 0 ){
										rel_err1 = loc_exp_err2*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim);
										rel_err2 = loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinError2(_sim_holes_5d[i][j]->GetBin(bin))/(loc_sim*loc_sim);
										rel_err3 = loc_sim_err2*loc_exp*loc_exp*_sim_holes_5d[i][j]->GetBinContent(bin)*_sim_holes_5d[i][j]->GetBinContent(bin)/(loc_sim*loc_sim*loc_sim*loc_sim);
										rel_err = TMath::Sqrt(rel_err1+rel_err2+rel_err3)/(_sim_holes_5d[i][j]->GetBinContent(bin)*loc_exp/loc_sim);
										sf_rel_err1 = loc_exp_err2/(loc_sim*loc_sim);
										sf_rel_err2 = loc_sim_err2*loc_exp*loc_exp/(loc_sim*loc_sim*loc_sim*loc_sim);
										sf_rel_err = TMath::Sqrt(sf_rel_err1+sf_rel_err2)/(loc_exp/loc_sim);
									}
									if(loc_exp>0.0 && loc_sim>0.0 && ((sf_rel_err <= 0.2 && sf_rel_err > 0.0) || dist>=max_dist_)){
										_scale_exp_5d[i][j]->SetBinContent(bin,loc_exp);
										_scale_sim_5d[i][j]->SetBinContent(bin,loc_sim);
										look_further_all=false;
										dist_dist[0][dist-1]++;
										_hole_radius_hist[i][j]->Fill(dist);
										_scale_exp_5d[i][j]->SetBinError2(_scale_exp_5d[i][j]->GetBin(bin),loc_exp_err2);
										//std::cout<<"Setting sim bin error\n";
										_scale_sim_5d[i][j]->SetBinError2(_scale_sim_5d[i][j]->GetBin(bin),loc_sim_err2);
										_sf_rel_err_hist[i][j]->Fill(sf_rel_err);
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
					_exp_data_5d_acc_corr[i][j]->GetAxis(n)->SetRange();
				}
				_sim_data_5d[i][j]->GetAxis(n)->SetRange();
				_sim_data_5d_acc_corr[i][j]->GetAxis(n)->SetRange();
			}
			wq2_bin++;
			std::cout<<"\nW:" <<Histogram::W_low(i) <<"-" <<Histogram::W_top(i) <<" Q2:" <<Histogram::Q2_low(j) <<"-" <<Histogram::Q2_top(j) <<" " <<wq2_bin <<"/" <<29*5 <<"\n";
			//std::cout<<"\tThrown Bins: " <<_thrown_5d[i][j]->GetNbins() <<" Recon Bins: " <<_sim_data_5d[i][j]->GetNbins() <<" Sim Holes Bins: " <<_sim_holes_5d[i][j]->GetNbins() <<"\n";
			//std::cout<<"\tShould have: " <<_thrown_5d[i][j]->GetNbins()-_sim_data_5d[i][j]->GetNbins() <<" bins localized\n";
			if(flags_.Flags::Helicity()){
				std::cout<<"\tScale Sim pos nbins = " <<_scale_sim_5d_pos[i][j]->GetNbins() <<" Scale exp pos nbins = " <<_scale_exp_5d_pos[i][j]->GetNbins() <<"\n"; 
				std::cout<<"\tScale Sim neg nbins = " <<_scale_sim_5d_neg[i][j]->GetNbins() <<" Scale exp neg nbins = " <<_scale_exp_5d_neg[i][j]->GetNbins() <<"\n"; 
			}else{


				double exp_nonlocal = 0.0;
				double sim_nonlocal = 0.0;

				double exp_err_nonlocal = 0.0;
				double sim_err_nonlocal = 0.0;

				for( long v_idx = 0; v_idx< _sim_data_5d_acc_corr[i][j]->GetNbins(); v_idx++){
					Double_t v_nonlocal = _sim_data_5d_acc_corr[i][j]->GetBinContent(v_idx,coord);
					if(v_nonlocal >0.0){
						_sim_acc_corr_rel_err[i][j]->Fill(_sim_data_5d_acc_corr[i][j]->GetBinError(v_idx)/_sim_data_5d_acc_corr[i][j]->GetBinContent(v_idx));
						Double_t w_nonlocal = _exp_data_5d_acc_corr[i][j]->GetBinContent(coord);
						exp_nonlocal += w_nonlocal;
						sim_nonlocal += v_nonlocal;
						exp_err_nonlocal += _exp_data_5d_acc_corr[i][j]->GetBinError2(_exp_data_5d_acc_corr[i][j]->GetBin(coord));//w_nonlocal*w_nonlocal;
						sim_err_nonlocal += _sim_data_5d_acc_corr[i][j]->GetBinError2(v_idx);//v_nonlocal*v_nonlocal;
					}
				}
				_sim_acc_corr_rel_err[i][j]->Write();
				_exp_acc_corr_rel_err[i][j]->Write();

				exp_err_nonlocal = TMath::Sqrt(exp_err_nonlocal);
				sim_err_nonlocal = TMath::Sqrt(sim_err_nonlocal);
				_value_global_sf[i][j]= exp_nonlocal/sim_nonlocal;
				_value_global_sf_err[i][j] = TMath::Sqrt((exp_err_nonlocal/sim_nonlocal)*(exp_err_nonlocal/sim_nonlocal)  + (exp_nonlocal*sim_err_nonlocal/(sim_nonlocal*sim_nonlocal))*(exp_nonlocal*sim_err_nonlocal/(sim_nonlocal*sim_nonlocal)));
				
				std::cout<<"exp_global" <<exp_nonlocal <<" sim_global:" <<sim_nonlocal <<"\n";
				std::cout<<"exp_err_global" <<exp_err_nonlocal <<" sim_err_global:" <<sim_err_nonlocal <<"\n";
				std::cout<<"global: " <<_value_global_sf[i][j] <<" error:" <<_value_global_sf_err[i][j] <<"\n";
				//std::cout<<"Method 1:\n";
				std::cout<<"\tloc exp: " <<exp_nonlocal <<" loc sim: " <<sim_nonlocal <<"\n";
				std::cout<<"\tNonLocal Scale would be = " <<_value_global_sf[i][j] <<" with error:" <<_value_global_sf_err[i][j] <<"\n";
				_global_sf->Fill(W_mid(i),Q2_mid(j),_value_global_sf[i][j]);
				//_global_sf_old->Fill(W_mid(i),Q2_mid(j),_value_global_sf_old[i][j]);
			}
			if(flags_.Flags::Helicity()){
				double exp_nonlocal_pos = 0.0;
				double sim_nonlocal_pos = 0.0;
				double exp_err_nonlocal_pos = 0.0;
				double sim_err_nonlocal_pos = 0.0;

				double exp_nonlocal_neg = 0.0;
				double sim_nonlocal_neg = 0.0;
				double exp_err_nonlocal_neg = 0.0;
				double sim_err_nonlocal_neg = 0.0;
				for( long v_idx = 0; v_idx< _sim_data_5d_acc_corr[i][j]->GetNbins(); v_idx++){
					Double_t v_nonlocal = _sim_data_5d_acc_corr[i][j]->GetBinContent(v_idx,coord);
					if(v_nonlocal >0.0){
						Double_t w_nonlocal_pos = _exp_data_5d_acc_corr_pos[i][j]->GetBinContent(coord);
						exp_nonlocal_pos += w_nonlocal_pos;
						sim_nonlocal_pos += v_nonlocal;
						exp_err_nonlocal_pos += w_nonlocal_pos*w_nonlocal_pos;
						sim_err_nonlocal_pos += v_nonlocal*v_nonlocal;

						Double_t w_nonlocal_neg = _exp_data_5d_acc_corr_neg[i][j]->GetBinContent(coord);
						exp_nonlocal_neg += w_nonlocal_neg;
						sim_nonlocal_neg += v_nonlocal;
						exp_err_nonlocal_neg += w_nonlocal_neg*w_nonlocal_neg;
						sim_err_nonlocal_neg += v_nonlocal*v_nonlocal;

					}
				}
				exp_err_nonlocal_pos = TMath::Sqrt(exp_err_nonlocal_pos);
				sim_err_nonlocal_pos = TMath::Sqrt(sim_err_nonlocal_pos);
				_value_global_sf_pos[i][j] = exp_nonlocal_pos/sim_nonlocal_pos;
				_value_global_sf_err_pos[i][j] = TMath::Sqrt((exp_err_nonlocal_pos/sim_nonlocal_pos)*(exp_err_nonlocal_pos/sim_nonlocal_pos)  + (exp_nonlocal_pos*sim_err_nonlocal_pos/(sim_nonlocal_pos*sim_nonlocal_pos))*(exp_nonlocal_pos*sim_err_nonlocal_pos/(sim_nonlocal_pos*sim_nonlocal_pos)));
				
				exp_err_nonlocal_neg = TMath::Sqrt(exp_err_nonlocal_neg);
				sim_err_nonlocal_neg = TMath::Sqrt(sim_err_nonlocal_neg);
				_value_global_sf_neg[i][j] = exp_nonlocal_neg/sim_nonlocal_neg;
				_value_global_sf_err_neg[i][j] = TMath::Sqrt((exp_err_nonlocal_neg/sim_nonlocal_neg)*(exp_err_nonlocal_neg/sim_nonlocal_neg)  + (exp_nonlocal_neg*sim_err_nonlocal_neg/(sim_nonlocal_neg*sim_nonlocal_neg))*(exp_nonlocal_neg*sim_err_nonlocal_neg/(sim_nonlocal_neg*sim_nonlocal_neg)));
				
				std::cout<<"\tglobal Scale Pos would be = " <<_value_global_sf_pos[i][j] <<" with error:" <<_value_global_sf_err_pos[i][j] <<"\n";
				std::cout<<"\tglobal Scale Neg would be = " <<_value_global_sf_neg[i][j] <<" with error:" <<_value_global_sf_err_neg[i][j] <<"\n";
				_global_sf_pos->Fill(W_mid(i),Q2_mid(j),exp_nonlocal_pos/sim_nonlocal_pos);
				_global_sf_neg->Fill(W_mid(i),Q2_mid(j),exp_nonlocal_neg/sim_nonlocal_neg);
			}
			/*for(int n=0; n<5; n++){
				if(flags_.Flags::Helicity()){
					_exp_data_5d_pos[i][j]->GetAxis(n)->SetRange(1,_n_bins_5d[n]);
					_exp_data_5d_neg[i][j]->GetAxis(n)->SetRange(1,_n_bins_5d[n]);
				}else{
					_exp_data_5d[i][j]->GetAxis(n)->SetRange(1,_n_bins_5d[n]);
				}
				_sim_data_5d[i][j]->GetAxis(n)->SetRange(1,_n_bins_5d[n]);
			}*/
			
			for(int k=0; k<3; k++){
				for(int l=0; l<24; l++){
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
				sprintf(hname,"Localized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_pos",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				_N_holes_5d_pos[i][j]->SetNameTitle(hname,hname);
				_N_holes_5d_pos[i][j]->Write();
				_N_holes_global_5d_pos[i][j]=(THnSparseD*)_sim_holes_5d[i][j]->Clone();
				double tmp_val = 0.0;
				double tmp_err = 0.0;
				for( long g_idx = 0; g_idx< _N_holes_global_5d_pos[i][j]->GetNbins(); g_idx++){
					tmp_val = _N_holes_global_5d_pos[i][j]->GetBinContent(g_idx);
					tmp_err = _N_holes_global_5d_pos[i][j]->GetBinError(g_idx);
					_N_holes_global_5d_pos[i][j]->SetBinContent(g_idx,tmp_val * _value_global_sf_pos[i][j]);
					_N_holes_global_5d_pos[i][j]->SetBinError(g_idx,TMath::Sqrt((tmp_err * _value_global_sf_pos[i][j])*(tmp_err * _value_global_sf_pos[i][j]) + (tmp_val * _value_global_sf_err_pos[i][j])*(tmp_val * _value_global_sf_err_pos[i][j])));
					tmp_val = 0.0;
					tmp_err = 0.0;
				}
				sprintf(hname,"Globalized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_pos",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				_N_holes_global_5d_pos[i][j]->SetNameTitle(hname,hname);
				_N_holes_global_5d_pos[i][j]->Write();

				_scale_5d_neg[i][j]=(THnSparseD*)_scale_exp_5d_neg[i][j]->Clone();
				_scale_5d_neg[i][j]->Divide(_scale_sim_5d_neg[i][j]);
				std::cout<<"Filling _sf_hist_neg\n";
				for (Long64_t l = 0; l < _scale_5d_neg[i][j]->GetNbins(); ++l) {
					_sf_hist_neg[i][j]->Fill(_scale_5d_neg[i][j]->GetBinContent(l));
				}
				_N_holes_5d_neg[i][j]=(THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_N_holes_5d_neg[i][j]->Multiply(_scale_5d_neg[i][j]);
				sprintf(hname,"Localized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_neg",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				_N_holes_5d_neg[i][j]->SetNameTitle(hname,hname);
				_N_holes_5d_neg[i][j]->Write();
				_N_holes_global_5d_neg[i][j]=(THnSparseD*)_sim_holes_5d[i][j]->Clone();
				for( long g_idx = 0; g_idx< _N_holes_global_5d_neg[i][j]->GetNbins(); g_idx++){
					tmp_val = _N_holes_global_5d_neg[i][j]->GetBinContent(g_idx);
					tmp_err = _N_holes_global_5d_neg[i][j]->GetBinError(g_idx);
					_N_holes_global_5d_neg[i][j]->SetBinContent(g_idx,tmp_val * _value_global_sf_neg[i][j]);
					_N_holes_global_5d_neg[i][j]->SetBinError(g_idx,TMath::Sqrt((tmp_err * _value_global_sf_neg[i][j])*(tmp_err * _value_global_sf_neg[i][j]) + (tmp_val * _value_global_sf_err_neg[i][j])*(tmp_val * _value_global_sf_err_neg[i][j])));
					tmp_val = 0.0;
					tmp_err = 0.0;
				}
				//_N_holes_global_5d_neg[i][j]->Scale(_global_sf_neg->GetBinContent(i,j));
				sprintf(hname,"Globalized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_neg",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				_N_holes_global_5d_neg[i][j]->SetNameTitle(hname,hname);
				_N_holes_global_5d_neg[i][j]->Write();

				for(long m=0; m<_N_holes_5d_pos[i][j]->GetNbins(); m++){
					_hole_err_hist_pos[i][j]->Fill(_N_holes_5d_pos[i][j]->GetBinError(m)/_N_holes_5d_pos[i][j]->GetBinContent(m));
				}

				for(long m=0; m<_N_holes_global_5d_pos[i][j]->GetNbins(); m++){
					_global_hole_err_hist_pos[i][j]->Fill(_N_holes_global_5d_pos[i][j]->GetBinError(m)/_N_holes_global_5d_pos[i][j]->GetBinContent(m));
				}

				for(long m=0; m<_N_holes_5d_neg[i][j]->GetNbins(); m++){
					_hole_err_hist_neg[i][j]->Fill(_N_holes_5d_neg[i][j]->GetBinError(m)/_N_holes_5d_neg[i][j]->GetBinContent(m));
				}

				for(long m=0; m<_N_holes_global_5d_neg[i][j]->GetNbins(); m++){
					_global_hole_err_hist_neg[i][j]->Fill(_N_holes_global_5d_neg[i][j]->GetBinError(m)/_N_holes_global_5d_neg[i][j]->GetBinContent(m));
				}

				_sf_hist_pos[i][j]->Write();
				_sf_hist_neg[i][j]->Write();
				_hole_err_hist_pos[i][j]->Write();
				_hole_err_hist_neg[i][j]->Write();
				_sf_rel_err_hist_pos[i][j]->Write();
				_sf_rel_err_hist_neg[i][j]->Write();
				_hole_radius_hist_pos[i][j]->Write();
				_hole_radius_hist_neg[i][j]->Write();
				_global_hole_err_hist_pos[i][j]->Write();
				_global_hole_err_hist_neg[i][j]->Write();
				std::cout<<"Wrote Helicity Histograms\n";
			}else{

				_scale_5d[i][j]=(THnSparseD*)_scale_exp_5d[i][j]->Clone();
				_scale_5d[i][j]->Divide(_scale_sim_5d[i][j]);

				std::cout<<"Filling _sf_hist\n";
				for (Long64_t l = 0; l < _scale_5d[i][j]->GetNbins(); ++l) {
					_sf_hist[i][j]->Fill(_scale_5d[i][j]->GetBinContent(l));
				}

				_N_holes_5d[i][j]=(THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_N_holes_5d[i][j]->Multiply(_scale_5d[i][j]);
				sprintf(hname,"Localized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				_N_holes_5d[i][j]->SetNameTitle(hname,hname);
				_N_holes_5d[i][j]->Write();

				double tmp_val = 0.0;
				double tmp_err = 0.0;
				_N_holes_global_5d[i][j]=(THnSparseD*)_sim_holes_5d[i][j]->Clone();
				for( long g_idx = 0; g_idx< _N_holes_global_5d[i][j]->GetNbins(); g_idx++){
					tmp_val = _N_holes_global_5d[i][j]->GetBinContent(g_idx);
					tmp_err = _N_holes_global_5d[i][j]->GetBinError(g_idx);
					_N_holes_global_5d[i][j]->SetBinContent(g_idx,tmp_val * _value_global_sf[i][j]);
					_N_holes_global_5d[i][j]->SetBinError(g_idx,TMath::Sqrt((tmp_err * _value_global_sf[i][j])*(tmp_err * _value_global_sf[i][j]) + (tmp_val * _value_global_sf_err[i][j])*(tmp_val * _value_global_sf_err[i][j])));
					tmp_val = 0.0;
					tmp_err = 0.0;
				}
				sprintf(hname,"Globalized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				_N_holes_global_5d[i][j]->SetNameTitle(hname,hname);
				_N_holes_global_5d[i][j]->Write();

				tmp_val = 0.0;
				tmp_err = 0.0;

				for(long m=0; m<_N_holes_5d[i][j]->GetNbins(); m++){
					_hole_err_hist[i][j]->Fill(_N_holes_5d[i][j]->GetBinError(m)/_N_holes_5d[i][j]->GetBinContent(m));
				}

				for(long m=0; m<_N_holes_global_5d[i][j]->GetNbins(); m++){
					_global_hole_err_hist[i][j]->Fill(_N_holes_global_5d[i][j]->GetBinError(m)/_N_holes_global_5d[i][j]->GetBinContent(m));
				}

				_sf_hist[i][j]->SetXTitle("Localized Scale Factor");
				_sf_hist[i][j]->SetYTitle("Number of Bins");
				_sf_hist[i][j]->Write();
				_hole_err_hist[i][j]->SetXTitle("Relative Error for Localized Hole");
				_hole_err_hist[i][j]->SetYTitle("Number of Bins");
				_hole_err_hist[i][j]->Write();
				_global_hole_err_hist[i][j]->SetXTitle("Relative Error for Globalized Hole");
				_global_hole_err_hist[i][j]->SetYTitle("Number of Bins");
				_global_hole_err_hist[i][j]->Write();
				_sf_rel_err_hist[i][j]->SetXTitle("Relative Error for Localized Scale Factor");
				_sf_rel_err_hist[i][j]->SetYTitle("Number of Bins");
				_sf_rel_err_hist[i][j]->Write();
				_hole_radius_hist[i][j]->SetXTitle("Localized Integration Radius Used");
				_hole_radius_hist[i][j]->SetYTitle("Number of Bins");
				_hole_radius_hist[i][j]->Write();
				std::cout<<"Wrote Non-Helicity Histograms\n";

				std::cout<<"Wrote Non-JHelicity Hold Hist\n";

				long num_hole_bins = 0;
				long num_exp_bins = 0;
				double num_hole_evts = 0.0;
				double num_exp_evts = 0.0;
				for(long m=0; m<_N_holes_5d[i][j]->GetNbins(); m++){
					if(_N_holes_5d[i][j]->GetBinContent(m)>0.0){
						num_hole_bins++;
						num_hole_evts += _N_holes_5d[i][j]->GetBinContent(m);
					}
				}
				for(long m=0; m<_exp_data_5d_acc_corr[i][j]->GetNbins(); m++){
					if(_exp_data_5d_acc_corr[i][j]->GetBinContent(m)>0.0){
						num_exp_bins++;
						num_exp_evts += _exp_data_5d_acc_corr[i][j]->GetBinContent(m);
					}
				}
				_frac_bins_holes_to_exp->Fill(W_mid(i),Q2_mid(j),num_hole_bins/num_exp_bins);
				_frac_events_holes_to_exp->Fill(W_mid(i),Q2_mid(j),num_hole_evts/num_exp_evts);
			}
		}
	}

	

	if(flags_.Flags::Helicity()){
		_global_sf_pos->SetXTitle("W (GeV)");
		_global_sf_pos->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_global_sf_pos->Write();
		_global_sf_neg->SetXTitle("W (GeV)");
		_global_sf_neg->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_global_sf_neg->Write();
	}else{
		_global_sf->SetXTitle("W (GeV)");
		_global_sf->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_global_sf->Write();
		_frac_bins_holes_to_exp->SetXTitle("W (GeV)");
		_frac_bins_holes_to_exp->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_frac_bins_holes_to_exp->Write();
		_frac_events_holes_to_exp->SetXTitle("W (GeV)");
		_frac_events_holes_to_exp->SetYTitle("^{}Q^{2} (^{}GeV^{2})");
		_frac_events_holes_to_exp->Write();
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

void Histogram::Make_Acceptance_Rel_Error(Flags flags_){
	char hname[100];
	for(int a=0; a<_W_nbins_; a++){
		for(int b=0; b<_Q2_nbins_; b++){
			//Int_t* coord = new Int_t[_acceptance_5d[i][j]->GetNbins()];
			_acceptance_rel_err_5d[a][b] = (THnSparseD*) _acceptance_5d[a][b]->Clone();
			sprintf(hname,"Acceptance Relative Error %s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(a),Histogram::W_top(a),Histogram::Q2_low(b),Histogram::Q2_top(b));
			_acceptance_rel_err_5d[a][b]->SetNameTitle(hname,hname);
			for(Long64_t i =0; i<_acceptance_5d[a][b]->GetNbins(); i++ ){
				if(_acceptance_5d[a][b]->GetBinContent(i)>0.0){
					_acceptance_rel_err_5d[a][b]->SetBinContent(i,_acceptance_5d[a][b]->GetBinError(i)/_acceptance_5d[a][b]->GetBinContent(i));
				}
			}
		}
	}
}


void Histogram::Acceptance_Rel_Error_Cut(Flags flags_){
	for(int a=0; a<_W_nbins_; a++){
		for(int b=0; b<_Q2_nbins_; b++){
			//Int_t* coord = new Int_t[_acceptance_5d[i][j]->GetNbins()];
			for(Long64_t i =0; i<_acceptance_5d[a][b]->GetNbins(); i++ ){
				if(_acceptance_rel_err_5d[a][b]->GetBinContent(i)>_Acceptance_Rel_Error_Max_[flags_.Flags::Acc_Rel_Error_Cut()]){
					_acceptance_5d[a][b]->SetBinContent(i,0.0);
					_acceptance_5d[a][b]->SetBinError(i,0.0);
				}
			}
		}
	}
}