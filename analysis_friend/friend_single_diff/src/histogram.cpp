#include "histogram.hpp"

Histogram::Histogram(std::unique_ptr<TFile>& exp_tree_, std::unique_ptr<TFile>& sim_tree_, std::unique_ptr<TFile>& empty_tree_, std::unique_ptr<TFile>& holes_, Flags flags_){
//Histogram::Histogram(std::unique_ptr<TFile>& exp_tree_, std::unique_ptr<TFile>& sim_tree_, std::unique_ptr<TFile>& empty_tree_, std::unique_ptr<TFile>& nr_sim_tree_, std::unique_ptr<TFile>& holes_, Flags flags_){
	Histogram::Make_Acceptance_Rel_Hist(flags_);
	Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,holes_,flags_);
	//Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,holes_,flags_);
    //Histogram::Rad_Corr();
    //Histogram::Sparse_7to5(flags_);
	Histogram::Fill_Acceptance_Rel_Hist();
    Histogram::Single_Diff(flags_);
	Histogram::Polarization(flags_);
	Histogram::Beam_Spin(flags_);
}

Histogram::Histogram(std::unique_ptr<TFile>& exp_tree_, std::unique_ptr<TFile>& sim_tree_, std::unique_ptr<TFile>& empty_tree_, Flags flags_){
//Histogram::Histogram(std::unique_ptr<TFile>& exp_tree_, std::unique_ptr<TFile>& sim_tree_, std::unique_ptr<TFile>& empty_tree_, std::unique_ptr<TFile>& nr_sim_tree_, Flags flags_){
    Histogram::Make_Acceptance_Rel_Hist(flags_);
	Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,flags_);
	//Histogram::Extract_5d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,flags_);
    //Histogram::Rad_Corr();
    //Histogram::Sparse_7to5(flags_);
	Histogram::Fill_Acceptance_Rel_Hist();
    Histogram::Single_Diff(flags_);
}

void Histogram::Extract_5d_Histograms(std::unique_ptr<TFile>& exp_tree_, std::unique_ptr<TFile>& sim_tree_, std::unique_ptr<TFile>& empty_tree_, std::unique_ptr<TFile>& holes_, Flags flags_){
//void Histogram::Extract_5d_Histograms(std::unique_ptr<TFile>& exp_tree_, std::unique_ptr<TFile>& sim_tree_, std::unique_ptr<TFile>& empty_tree_, std::unique_ptr<TFile>& nr_sim_tree_, std::unique_ptr<TFile>& holes_, Flags flags_){
    std::cout<<"Extract 5d Histograms including kinematic holes\n";
	char hname[100];
    for(int i=0; i<_W_nbins_; i++){
	//for(int i=0; i<20; i++){
        for(int j=0; j<_Q2_nbins_; j++){
            //sprintf(hname,"Thrown_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			sprintf(hname,"Thrown_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            std::cout<<"extracting " <<hname <<"\n";
			_thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            //_thrown_no_rad_5d[i][j] = (THnSparseD *)nr_sim_tree_->Get(hname);
            //sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			sprintf(hname,"%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Read_Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            std::cout<<"extracting " <<hname <<"\n";
			std::cout<<"from sim";
			_sim_data_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
			std::cout<<" done\ncalculating acceptance";
			_acceptance_5d[i][j] = (THnSparseD*)_sim_data_5d[i][j]->Clone();
            _acceptance_5d[i][j]->Divide(_thrown_5d[i][j]);
            std::cout<<"initial division done\n";
			//sprintf(hname,"Acceptance_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
			sprintf(hname,"Acceptance_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            _acceptance_5d[i][j]->SetNameTitle(hname,hname);
            std::cout<<"AREC cut perform:";
			for(long k=0; k<_acceptance_5d[i][j]->GetNbins(); k++){
				if(_acceptance_5d[i][j]->GetBinContent(k) >0.0){
					if(_acceptance_5d[i][j]->GetBinError(k)/_acceptance_5d[i][j]->GetBinContent(k) > _Acceptance_Rel_Error_Max_[flags_.Flags::Acc_Rel_Error_Cut()]){
						_acceptance_5d[i][j]->SetBinContent(k,0.0);
						_acceptance_5d[i][j]->SetBinError(k,0.0);
					}
				}
			}
			std::cout<<" done\nGoing into exp files:";
			if(flags_.Flags::Plot_Beam_Spin()){
				//sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				sprintf(hname,"%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Read_Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				std::cout<<"extracting " <<hname <<"\n";
				std::cout<<"exp pos: ";
				_exp_data_5d_pos[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            	std::cout<<" done\nemp pos: ";
				_empty_5d_pos[i][j] = (THnSparseD *)empty_tree_->Get(hname);
				std::cout<<"done\n";
				//sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				sprintf(hname,"%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Read_Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				std::cout<<"extracting " <<hname <<"\n";
				_exp_data_5d_neg[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            	_empty_5d_neg[i][j] = (THnSparseD *)empty_tree_->Get(hname);
				//if(flags_.Flags::Plot_Local()){
				if(true){
					/*
					if(flags_.Flags::Plot_Local_Holes()){
						sprintf(hname,"Localized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_pos",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
						std::cout<<"extracting " <<hname <<"\n";
						_N_local_holes_5d_pos[i][j] = (THnSparseD *)holes_->Get(hname);
						std::cout<<"adding local holes to pos\n";
						sprintf(hname,"Localized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_neg",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
						_N_local_holes_5d_neg[i][j] = (THnSparseD *)holes_->Get(hname);
						std::cout<<"adding local holes to neg\n";
					}else{
					*/

					sprintf(hname,"N_local_holes_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					std::cout<<"making " <<hname <<"\n";
					_N_5d_local_holes_pos[i][j] = (THnSparseD*)_exp_data_5d_pos[i][j]->Clone();
					_N_5d_local_holes_pos[i][j]->SetNameTitle(hname,hname);
					_N_5d_local_holes_pos[i][j]->Add(_empty_5d_pos[i][j],-17.074);//Empty target subtraction
					_N_5d_local_holes_pos[i][j]->Divide(_acceptance_5d[i][j]);
					sprintf(hname,"Localized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_pos",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					std::cout<<"adding "<<i <<" " <<j <<" " <<hname <<"\n";
					_N_local_holes_5d_pos[i][j] = (THnSparseD *)holes_->Get(hname);
					_N_5d_local_holes_pos[i][j]->Add(_N_local_holes_5d_pos[i][j]);
					/*
						if((THnSparseD *)holes_->Get(hname)){
							_N_local_holes_5d_pos[i][j] = (THnSparseD *)holes_->Get(hname);
							_N_5d_local_holes_pos[i][j]->Add(_N_local_holes_5d_pos[i][j]);
						}else{
							std::cout<<"returned null pointer\n";
						}
					*/
						//_N_5d_local_holes_pos[i][j]->Add((THnSparseD *)holes_->Get(hname));

					sprintf(hname,"N_local_holes_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					std::cout<<"making " <<hname <<"\n";
					_N_5d_local_holes_neg[i][j] = (THnSparseD*)_exp_data_5d_neg[i][j]->Clone();
					_N_5d_local_holes_neg[i][j]->SetNameTitle(hname,hname);
					_N_5d_local_holes_neg[i][j]->Add(_empty_5d_neg[i][j],-17.468);//Empty target subtraction
					_N_5d_local_holes_neg[i][j]->Divide(_acceptance_5d[i][j]);
					sprintf(hname,"Localized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_neg",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					std::cout<<"adding " <<hname <<"\n";
					_N_local_holes_5d_neg[i][j] = (THnSparseD *)holes_->Get(hname);
					_N_5d_local_holes_neg[i][j]->Add(_N_local_holes_5d_neg[i][j]);
					/*
						if((THnSparseD *)holes_->Get(hname)){
							_N_local_holes_5d_neg[i][j] = (THnSparseD *)holes_->Get(hname);
							_N_5d_local_holes_neg[i][j]->Add(_N_local_holes_5d_neg[i][j]);
						}else{
							std::cout<<"returned null pointer\n";
						}
					*/
						//_N_5d_local_holes_neg[i][j]->Add((THnSparseD *)holes_->Get(hname));
					//}
					
					/*
					sprintf(hname,"N_local_holes_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					std::cout<<"Making " <<hname <<"\n";
					_N_5d_local_holes_neg[i][j] = (THnSparseD*)_exp_data_5d_neg[i][j]->Clone();
					_N_5d_local_holes_neg[i][j]->SetNameTitle(hname,hname);
					_N_5d_local_holes_neg[i][j]->Add(_empty_5d_neg[i][j],-flags_.Flags::Qr());//Empty target subtraction
					_N_5d_local_holes_neg[i][j]->Divide(_acceptance_5d[i][j]);

					sprintf(hname,"Localized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_neg",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					std::cout<<"extracting " <<hname <<"\n";
					_N_local_holes_5d_neg[i][j] = (THnSparseD *)holes_->Get(hname);
					std::cout<<"adding local holes to neg\n";
					_N_5d_local_holes_neg[i][j]->Add(_N_local_holes_5d_neg[i][j]);
					*/
				}

				//if(flags_.Flags::Plot_Global()){
				sprintf(hname,"N_global_holes_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"making " <<hname <<"\n";
				_N_5d_global_holes_pos[i][j] = (THnSparseD*)_exp_data_5d_pos[i][j]->Clone();
				_N_5d_global_holes_pos[i][j]->SetNameTitle(hname,hname);
				_N_5d_global_holes_pos[i][j]->Add(_empty_5d_pos[i][j],-17.074);//Empty target subtraction
				_N_5d_global_holes_pos[i][j]->Divide(_acceptance_5d[i][j]);

				sprintf(hname,"Globalized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_pos",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"extracting " <<hname <<"\n";
				_N_global_holes_5d_pos[i][j] = (THnSparseD *)holes_->Get(hname);
				std::cout<<"adding global holes to pos\n";
				_N_5d_global_holes_pos[i][j]->Add(_N_global_holes_5d_pos[i][j]);

				sprintf(hname,"N_global_holes_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_global_holes_neg[i][j] = (THnSparseD*)_exp_data_5d_neg[i][j]->Clone();
				_N_5d_global_holes_neg[i][j]->SetNameTitle(hname,hname);
				_N_5d_global_holes_neg[i][j]->Add(_empty_5d_neg[i][j],-17.468);//Empty target subtraction
				_N_5d_global_holes_neg[i][j]->Divide(_acceptance_5d[i][j]);

				sprintf(hname,"Globalized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_neg",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"extracting " <<hname <<"\n";
				_N_global_holes_5d_neg[i][j] = (THnSparseD *)holes_->Get(hname);
				std::cout<<"adding global holes to neg\n";
				_N_5d_global_holes_neg[i][j]->Add(_N_global_holes_5d_neg[i][j]);
				//}
				
				//if(flags_.Flags::Plot_Raw()){
				sprintf(hname,"N_raw_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"making " <<hname <<"\n";
				_N_5d_raw_yield_pos[i][j] = (THnSparseD*)_exp_data_5d_pos[i][j]->Clone();
				_N_5d_raw_yield_pos[i][j]->SetNameTitle(hname,hname);
				_N_5d_raw_yield_pos[i][j]->Add(_empty_5d_pos[i][j],-17.074);//Empty target subtraction
				_N_5d_raw_yield_pos[i][j]->Divide(_acceptance_5d[i][j]);

				sprintf(hname,"N_raw_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_raw_yield_neg[i][j] = (THnSparseD*)_exp_data_5d_neg[i][j]->Clone();
				_N_5d_raw_yield_neg[i][j]->SetNameTitle(hname,hname);
				_N_5d_raw_yield_neg[i][j]->Add(_empty_5d_neg[i][j],-17.468);//Empty target subtraction
				_N_5d_raw_yield_neg[i][j]->Divide(_acceptance_5d[i][j]);
				//}
			
				/*
				sprintf(hname,"Raw_Sum_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_raw_yield_sum[i][j] = (THnSparseD*)_N_5d_raw_yield_pos[i][j]->Clone();
				_N_5d_raw_yield_sum[i][j]->Add(_N_5d_raw_yield_neg[i][j],1.0);
				_N_5d_raw_yield_sum[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"Raw_Diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_raw_yield_diff[i][j] = (THnSparseD*)_N_5d_raw_yield_pos[i][j]->Clone();
				_N_5d_raw_yield_diff[i][j]->Add(_N_5d_raw_yield_neg[i][j],-1.0);
				_N_5d_raw_yield_diff[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"Raw_Ratio_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_raw_yield_ratio[i][j] = (THnSparseD*)_N_5d_raw_yield_diff[i][j]->Clone();
				_N_5d_raw_yield_ratio[i][j]->Divide(_N_5d_raw_yield_sum[i][j]);
				_N_5d_raw_yield_ratio[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"N_Local_Sum_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_local_holes_sum[i][j] = (THnSparseD*)_N_5d_local_holes_pos[i][j]->Clone();
				_N_5d_local_holes_sum[i][j]->Add(_N_5d_local_holes_neg[i][j],1.0);
				_N_5d_local_holes_sum[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"N_Local_Diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_local_holes_diff[i][j] = (THnSparseD*)_N_5d_local_holes_pos[i][j]->Clone();
				_N_5d_local_holes_diff[i][j]->Add(_N_5d_local_holes_neg[i][j],-1.0);
				_N_5d_local_holes_diff[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"N_Local_Ratio_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_local_holes_ratio[i][j] = (THnSparseD*)_N_5d_local_holes_diff[i][j]->Clone();
				_N_5d_local_holes_ratio[i][j]->Divide(_N_5d_local_holes_sum[i][j]);
				_N_5d_local_holes_ratio[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"N_Global_Sum_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_global_holes_sum[i][j] = (THnSparseD*)_N_5d_global_holes_pos[i][j]->Clone();
				_N_5d_global_holes_sum[i][j]->Add(_N_5d_global_holes_neg[i][j],1.0);
				_N_5d_global_holes_sum[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"N_Global_Diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_global_holes_diff[i][j] = (THnSparseD*)_N_5d_global_holes_pos[i][j]->Clone();
				_N_5d_global_holes_diff[i][j]->Add(_N_5d_global_holes_neg[i][j],-1.0);
				_N_5d_global_holes_diff[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"N_Global_Ratio_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_5d_global_holes_ratio[i][j] = (THnSparseD*)_N_5d_global_holes_diff[i][j]->Clone();
				_N_5d_global_holes_ratio[i][j]->Divide(_N_5d_global_holes_sum[i][j]);
				_N_5d_global_holes_ratio[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"Local_Holes_Sum_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_local_holes_5d_sum[i][j] = (THnSparseD*)_N_local_holes_5d_pos[i][j]->Clone();
				_N_local_holes_5d_sum[i][j]->Add(_N_local_holes_5d_neg[i][j],1.0);
				_N_local_holes_5d_sum[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"Local_Holes_Diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_local_holes_5d_diff[i][j] = (THnSparseD*)_N_local_holes_5d_pos[i][j]->Clone();
				_N_local_holes_5d_diff[i][j]->Add(_N_local_holes_5d_neg[i][j],-1.0);
				_N_local_holes_5d_diff[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"Local_Holes_Ratio_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_local_holes_5d_ratio[i][j] = (THnSparseD*)_N_local_holes_5d_diff[i][j]->Clone();
				_N_local_holes_5d_ratio[i][j]->Divide(_N_local_holes_5d_sum[i][j]);
				_N_local_holes_5d_ratio[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"Global_Holes_Sum_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_global_holes_5d_sum[i][j] = (THnSparseD*)_N_global_holes_5d_pos[i][j]->Clone();
				_N_global_holes_5d_sum[i][j]->Add(_N_global_holes_5d_neg[i][j],1.0);
				_N_global_holes_5d_sum[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"Global_Holes_Diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_global_holes_5d_diff[i][j] = (THnSparseD*)_N_global_holes_5d_pos[i][j]->Clone();
				_N_global_holes_5d_diff[i][j]->Add(_N_global_holes_5d_neg[i][j],-1.0);
				_N_global_holes_5d_diff[i][j]->SetNameTitle(hname,hname);

				sprintf(hname,"Global_Holes_Ratio_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"Making " <<hname <<"\n";
				_N_global_holes_5d_ratio[i][j] = (THnSparseD*)_N_global_holes_5d_diff[i][j]->Clone();
				_N_global_holes_5d_ratio[i][j]->Divide(_N_global_holes_5d_sum[i][j]);
				_N_global_holes_5d_ratio[i][j]->SetNameTitle(hname,hname);
				*/
			}else{
				//sprintf(hname,"%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				sprintf(hname,"%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Read_Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
				std::cout<<"extracting exp hist: " <<hname <<" =";
				_exp_data_5d[i][j] = (THnSparseD *)exp_tree_->Get(hname);
            	std::cout<<" done\nextracting emp hist: " <<hname <<" =";
				_empty_5d[i][j] = (THnSparseD *)empty_tree_->Get(hname);
				std::cout<<" done\n";
				std::cout<<flags_.Flags::Qr() <<"\n";
				sprintf(hname,"N_local_holes_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				_N_5d_local_holes[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
				_N_5d_local_holes[i][j]->SetNameTitle(hname,hname);
				_N_5d_local_holes[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
				_N_5d_local_holes[i][j]->Divide(_acceptance_5d[i][j]);

				sprintf(hname,"N_global_holes_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				_N_5d_global_holes[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
				_N_5d_global_holes[i][j]->SetNameTitle(hname,hname);
				_N_5d_global_holes[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
				_N_5d_global_holes[i][j]->Divide(_acceptance_5d[i][j]);

				sprintf(hname,"N_raw_%s_%s_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				_N_5d_raw_yield[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
				_N_5d_raw_yield[i][j]->SetNameTitle(hname,hname);
				_N_5d_raw_yield[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
				_N_5d_raw_yield[i][j]->Divide(_acceptance_5d[i][j]);

				sprintf(hname,"Localized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"extracting " <<hname <<"\n";
				_N_local_holes_5d[i][j] = (THnSparseD *)holes_->Get(hname);
				std::cout<<"adding local holes to total yield\n";
				_N_5d_local_holes[i][j]->Add(_N_local_holes_5d[i][j]);
				
				sprintf(hname,"Globalized_Holes_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_AREC:%s",Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				std::cout<<"extracting " <<hname <<"\n";
				_N_global_holes_5d[i][j] = (THnSparseD *)holes_->Get(hname);
				std::cout<<"adding global holes to total yield\n";
				_N_5d_global_holes[i][j]->Add(_N_global_holes_5d[i][j]);
			}

			
	        
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

//void Histogram::Extract_5d_Histograms(std::unique_ptr<TFile>& exp_tree_, std::unique_ptr<TFile>& sim_tree_, std::unique_ptr<TFile>& empty_tree_, std::unique_ptr<TFile>& nr_sim_tree_, Flags flags_){
void Histogram::Extract_5d_Histograms(std::unique_ptr<TFile>& exp_tree_, std::unique_ptr<TFile>& sim_tree_, std::unique_ptr<TFile>& empty_tree_, Flags flags_){
    std::cout<<"Extract 5d Histograms\n";
	char hname[100];
	sprintf(hname,"Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
	std::cout<<"Getting Thrown THnSparse " <<hname <<"\n";
    for(int i=0; i<_W_nbins_; i++){
        for(int j=0; j<_Q2_nbins_; j++){
            sprintf(hname,"Thrown_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j));
            std::cout<<"Getting Thrown THnSparses " <<hname <<"\n";
			_thrown_5d[i][j] = (THnSparseD *)sim_tree_->Get(hname);
            //_thrown_no_rad_5d[i][j] = (THnSparseD *)nr_sim_tree_->Get(hname);
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
				if(_acceptance_5d[i][j]->GetBinError(k)/_acceptance_5d[i][j]->GetBinContent(k) > 0.75){
					_acceptance_5d[i][j]->SetBinContent(k,0.0);
				}
			}
			std::cout<<"Making Yield\n";
            _N_5d[i][j] = (THnSparseD*)_exp_data_5d[i][j]->Clone();
            _N_5d[i][j]->Add(_empty_5d[i][j],-flags_.Flags::Qr());//Empty target subtraction
	        _N_5d[i][j]->Divide(_acceptance_5d[i][j]);
			/*
			if(flags_.Flags::Nonlocal_Holes()){
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
			}
			*/
        }
    }
    for(int i = 0; i<_thrown_5d[0][0]->GetNdimensions(); i++){
		//_n_bins_7d.push_back(_thrown_7d->GetAxis(i)->GetNbins());
        _n_bins_5d.push_back(_thrown_5d[0][0]->GetAxis(i)->GetNbins());
	}
}
/*
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
/*}
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
	char hname_raw[100];
	char hname_local[100];
	char hname_global[100];
	char hname_local_holes[100];
	char hname_global_holes[100];
	char xlabel[100];
	char ylabel[100];
	int idx_2d[2]; 
	std::cout<<"Making Output File\n";
	_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	std::cout<<"Made File:" <<flags_.Flags::Output_File().c_str() <<"\n";
	_RootOutputFile->cd();
	//Get the 7dimensional bins ready
	TH1D_1d_star exp_ch_1d_raw;
	TH1D_2d_star exp_ch_2d_raw;

	TH1D_1d_star exp_ch_1d_local;
	TH1D_2d_star exp_ch_2d_local;

	TH1D_1d_star exp_ch_1d_global;
	TH1D_2d_star exp_ch_2d_global;

	TH1D_1d_star exp_ch_1d_local_holes;
	TH1D_1d_star exp_ch_1d_global_holes;
	std::cout<<"\tMaking Directory\n";
	TDirectory* dir_S = _RootOutputFile->mkdir("Single Differential CS");
    //std::cout<<"\t\tMade first directory\n";
	dir_S->cd();
    //_rad_corr->Write();
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
			sprintf(dirname,"%s_Q2|[%.2f,%.2f)",_five_dim_[k],Histogram::Q2_low(j),Histogram::Q2_top(j));
			dir_S2[k][j] = dir_S1[k]->mkdir(dirname);
		}
	}
	sprintf(hname,"Filled-Target W-Q2 Yield");
	TH2D* _WQ2_filled_target = new TH2D(hname,hname,_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
	_WQ2_filled_target->GetYaxis()->Set(5,_Q2_bins_);
	sprintf(hname,"Empty-Target W-Q2 Yield");
	TH2D* _WQ2_empty_target = new TH2D(hname,hname,_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
	_WQ2_empty_target->GetYaxis()->Set(5,_Q2_bins_);
	sprintf(hname,"Sim Recon W-Q2 Yield");
	TH2D* _WQ2_sim_recon = new TH2D(hname,hname,_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
	_WQ2_sim_recon->GetYaxis()->Set(5,_Q2_bins_);
	sprintf(hname,"Sim Thrown W-Q2 Yield");
	TH2D* _WQ2_sim_thrown = new TH2D(hname,hname,_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
	_WQ2_sim_thrown->GetYaxis()->Set(5,_Q2_bins_);
	sprintf(hname,"Global Holes W-Q2 Yield");
	TH2D* _WQ2_global_holes = new TH2D(hname,hname,_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
	_WQ2_global_holes->GetYaxis()->Set(5,_Q2_bins_);
	sprintf(hname,"Local Holes W-Q2 Yield");
	TH2D* _WQ2_local_holes = new TH2D(hname,hname,_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
	_WQ2_local_holes->GetYaxis()->Set(5,_Q2_bins_);
	sprintf(hname,"Yield with Local Holes W-Q2 Yield");
	TH2D* _WQ2_yield_local = new TH2D(hname,hname,_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
	_WQ2_yield_local->GetYaxis()->Set(5,_Q2_bins_);
	sprintf(hname,"Yield with Global Holes W-Q2 Yield");
	TH2D* _WQ2_yield_global = new TH2D(hname,hname,_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
	_WQ2_yield_global->GetYaxis()->Set(5,_Q2_bins_);
	sprintf(hname,"Yield No Holes W-Q2 Yield");
	TH2D* _WQ2_yield_raw = new TH2D(hname,hname,_W_nbins_,_W_min_,_W_max_,_Q2_nbins_,_Q2_min_,_Q2_max_);
	_WQ2_yield_raw->GetYaxis()->Set(5,_Q2_bins_);


	for(int Wbin=0; Wbin< _W_nbins_; Wbin++){
		for(int Q2bin=0; Q2bin< _Q2_nbins_; Q2bin++){
			for(long i=0; i<_exp_data_5d[Wbin][Q2bin]->GetNbins(); i++){
				_WQ2_filled_target->Fill(Histogram::W_mid(Wbin),Histogram::Q2_mid(Q2bin),_exp_data_5d[Wbin][Q2bin]->GetBinContent(i));
			}
			for(long i=0; i<_empty_5d[Wbin][Q2bin]->GetNbins(); i++){
				_WQ2_empty_target->Fill(Histogram::W_mid(Wbin),Histogram::Q2_mid(Q2bin),_empty_5d[Wbin][Q2bin]->GetBinContent(i));
			}
			for(long i=0; i<_sim_data_5d[Wbin][Q2bin]->GetNbins(); i++){
				_WQ2_sim_recon->Fill(Histogram::W_mid(Wbin),Histogram::Q2_mid(Q2bin),_sim_data_5d[Wbin][Q2bin]->GetBinContent(i));
			}
			for(long i=0; i<_thrown_5d[Wbin][Q2bin]->GetNbins(); i++){
				_WQ2_sim_thrown->Fill(Histogram::W_mid(Wbin),Histogram::Q2_mid(Q2bin),_thrown_5d[Wbin][Q2bin]->GetBinContent(i));
			}
			for(long i=0; i<_N_5d_local_holes[Wbin][Q2bin]->GetNbins(); i++){
				_WQ2_yield_local->Fill(Histogram::W_mid(Wbin),Histogram::Q2_mid(Q2bin),_N_5d_local_holes[Wbin][Q2bin]->GetBinContent(i));
			}
			for(long i=0; i<_N_5d_global_holes[Wbin][Q2bin]->GetNbins(); i++){
				_WQ2_yield_global->Fill(Histogram::W_mid(Wbin),Histogram::Q2_mid(Q2bin),_N_5d_global_holes[Wbin][Q2bin]->GetBinContent(i));
			}
			for(long i=0; i<_N_5d_raw_yield[Wbin][Q2bin]->GetNbins(); i++){
				_WQ2_yield_raw->Fill(Histogram::W_mid(Wbin),Histogram::Q2_mid(Q2bin),_N_5d_raw_yield[Wbin][Q2bin]->GetBinContent(i));
			}
			for(long i=0; i<_N_local_holes_5d[Wbin][Q2bin]->GetNbins(); i++){
				_WQ2_local_holes->Fill(Histogram::W_mid(Wbin),Histogram::Q2_mid(Q2bin),_N_local_holes_5d[Wbin][Q2bin]->GetBinContent(i));
			}
			for(long i=0; i<_N_global_holes_5d[Wbin][Q2bin]->GetNbins(); i++){
				_WQ2_global_holes->Fill(Histogram::W_mid(Wbin),Histogram::Q2_mid(Q2bin),_N_global_holes_5d[Wbin][Q2bin]->GetBinContent(i));
			}
			
		}
	}
	dir_C1->cd();
	_WQ2_filled_target->SetXTitle("W (GeV)");
	_WQ2_filled_target->SetYTitle("Q2 (GeV^2)");
	_WQ2_filled_target->Write();

	_WQ2_empty_target->SetXTitle("W (GeV)");
	_WQ2_empty_target->SetYTitle("Q2 (GeV^2)");
	_WQ2_empty_target->Write();
	_WQ2_sim_recon->SetXTitle("W (GeV)");
	_WQ2_sim_recon->SetYTitle("Q2 (GeV^2)");
	_WQ2_sim_recon->Write();

	_WQ2_sim_thrown->SetXTitle("W (GeV)");
	_WQ2_sim_thrown->SetYTitle("Q2 (GeV^2)");
	_WQ2_sim_thrown->Write();

	_WQ2_yield_raw->SetXTitle("W (GeV)");
	_WQ2_yield_raw->SetYTitle("Q2 (GeV^2)");
	_WQ2_yield_raw->Write();

	_WQ2_yield_local->SetXTitle("W (GeV)");
	_WQ2_yield_local->SetYTitle("Q2 (GeV^2)");
	_WQ2_yield_local->Write();

	_WQ2_yield_global->SetXTitle("W (GeV)");
	_WQ2_yield_global->SetYTitle("Q2 (GeV^2)");
	_WQ2_yield_global->Write();

	_WQ2_local_holes->SetXTitle("W (GeV)");
	_WQ2_local_holes->SetYTitle("Q2 (GeV^2)");
	_WQ2_local_holes->Write();

	_WQ2_global_holes->SetXTitle("W (GeV)");
	_WQ2_global_holes->SetYTitle("Q2 (GeV^2)");
	_WQ2_global_holes->Write();


    
	double denom = 1.0;
	std::cout<<"\tDirectories Made\n\tWriting Histograms\n";
	for(int i=0; i<_W_nbins_; i++){//W
		for(int j=0; j< _Q2_nbins_; j++){//Q2
			for(int k=0; k<5; k++){
				dir_S2[k][j]->cd();
				//std::cout<<"\t\tWriting Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<k <<"\n";
				sprintf(hname_raw,"raw_%s_single_diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_top:%s_var:%s_AREC:%s",_five_dim_[k],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str(),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				sprintf(hname_local,"local_%s_single_diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_top:%s_var:%s_AREC:%s",_five_dim_[k],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str(),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				sprintf(hname_global,"global_%s_single_diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_top:%s_var:%s_AREC:%s",_five_dim_[k],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str(),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				sprintf(hname_local_holes,"local_holes_%s_single_diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_top:%s_var:%s_AREC:%s",_five_dim_[k],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str(),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				sprintf(hname_global_holes,"global_holes_%s_single_diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_top:%s_var:%s_AREC:%s",_five_dim_[k],Histogram::W_low(i),Histogram::W_top(i),Histogram::Q2_low(j),Histogram::Q2_top(j),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str(),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
				sprintf(xlabel,"%s %s",_five_dim_[k],_dim_units_[k]);
				sprintf(ylabel,"Diff CS Acceptance Corr Yield/(%s))",_dim_units_y_[k]);
                //std::cout<<"\tpushing back projection " <<_N_5d[i][j]->Projection(k,"E") <<"\n";
				exp_ch_1d_raw.push_back(_N_5d_raw_yield[i][j]->Projection(k,"E"));
				exp_ch_1d_local.push_back(_N_5d_local_holes[i][j]->Projection(k,"E"));
				exp_ch_1d_global.push_back(_N_5d_global_holes[i][j]->Projection(k,"E"));
				exp_ch_1d_local_holes.push_back(_N_local_holes_5d[i][j]->Projection(k,"E"));
				exp_ch_1d_global_holes.push_back(_N_global_holes_5d[i][j]->Projection(k,"E"));
				//std::cout<<"\tDenominator time\n\t\tvirtual photon flux\n";
                
				//denom *= physics::Virtual_Photon_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]);
				
				//std::cout<<"Scaling Denominator post flux: " <<denom <<"\n";
				//std::cout<<"What it would have been with old flux: " <<denom/physics::ratio_Virtual_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]) <<"\n";
				//std::cout<<"\tThe ratio of difference is: " <<physics::ratio_Virtual_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]) <<"\n";
				//std::cout<<"Old virtual flux: " <<physics::old_Virtual_Photon_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]) <<" vs new: " <<physics::Virtual_Photon_Flux((double)Histogram::W_mid(i),(double)Histogram::Q2_mid(j),_beam_energy_[0]) <<"\n";
                //std::cout<<"\t\tluminosity\n";
                
				//denom *=flags_.Flags::L(0);
                
				//std::cout<<"Scaling Denominator post Luminosity: " <<denom <<"\n";
                //std::cout<<"\t\trad corr\n";
				
				/*if(flags_.Flags::Rad_Corr()){
					denom *= _rad_corr_array[i][j];
                    //std::cout<<"Scaling Denominator post rad corr: " <<denom <<"\n";
				}*/

                //std::cout<<"\t\tW bin\n";
				
				//denom *=_W_res_;
                
				//std::cout<<"Scaling Denominator post W: " <<denom <<"\n";
                //std::cout<<"\t\tQ2 bin\n";
			    
				//denom *=(_Q2_bins_[j+1]-_Q2_bins_[j]);
                
				//std::cout<<"Scaling Denominator post Q2: " <<denom <<"\n";
				/*if(k>1){
					//For Theta this will need to be undone and then modified for cosine theta
					denom *=(_thrown_5d[i][j]->GetAxis(k)->GetBinUpEdge(2)-_thrown_5d[i][j]->GetAxis(k)->GetBinLowEdge(2))*TMath::Pi()/180.0;//Divide by angle bins in radians
                    //denom *= Histogram::CosTheta()
                    //std::cout<<"Scaling Denominator post Xij: " <<denom <<"\n";
                }else{
                    //std::cout<<"\t\txij bin\n";
					denom *=(_thrown_5d[i][j]->GetAxis(k)->GetBinUpEdge(2)-_thrown_5d[i][j]->GetAxis(k)->GetBinLowEdge(2));//Phi in radians
                    //std::cout<<"Scaling Denominator post Xij: " <<denom <<"\n";
				}*/

				//std::cout<<"\tScaling\n";
                //std::cout<<"Current Integral" <<exp_ch_1d[k]->Integral() <<"\n";
                //std::cout<<"Scaling Denominator: " <<denom <<"\n";
				
				//exp_ch_1d[k]->Scale(1.0/denom);
                
				//std::cout<<"Post Integral" <<exp_ch_1d[k]->Integral() <<"\n";
			
				exp_ch_1d_raw[k]->SetNameTitle(hname_raw,hname_raw);
				exp_ch_1d_raw[k]->GetXaxis()->SetTitle(xlabel);
				exp_ch_1d_raw[k]->GetYaxis()->SetTitle(ylabel);

				exp_ch_1d_local[k]->SetNameTitle(hname_local,hname_local);
				exp_ch_1d_local[k]->GetXaxis()->SetTitle(xlabel);
				exp_ch_1d_local[k]->GetYaxis()->SetTitle(ylabel);

				exp_ch_1d_global[k]->SetNameTitle(hname_global,hname_global);
				exp_ch_1d_global[k]->GetXaxis()->SetTitle(xlabel);
				exp_ch_1d_global[k]->GetYaxis()->SetTitle(ylabel);

				exp_ch_1d_local_holes[k]->SetNameTitle(hname_local_holes,hname_local_holes);
				exp_ch_1d_local_holes[k]->GetXaxis()->SetTitle(xlabel);
				exp_ch_1d_local_holes[k]->GetYaxis()->SetTitle(ylabel);

				exp_ch_1d_global_holes[k]->SetNameTitle(hname_global_holes,hname_global_holes);
				exp_ch_1d_global_holes[k]->GetXaxis()->SetTitle(xlabel);
				exp_ch_1d_global_holes[k]->GetYaxis()->SetTitle(ylabel);
				//std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<k <<"\n";
				exp_ch_1d_raw[k]->Write();
				exp_ch_1d_local[k]->Write();
				exp_ch_1d_global[k]->Write();
				exp_ch_1d_local_holes[k]->Write();
				exp_ch_1d_global_holes[k]->Write();
				//denom = 1.0;
			}
			//exp_ch_2d.push_back(exp_ch_1d);
			exp_ch_1d_raw.clear();
			exp_ch_1d_local.clear();
			exp_ch_1d_global.clear();
			exp_ch_1d_local_holes.clear();
			exp_ch_1d_global_holes.clear();
			_acceptance_rel_err_hist[i][j]->GetXaxis()->SetTitle("Relative Error in Acceptance");
			_acceptance_rel_err_hist[i][j]->GetYaxis()->SetTitle("Number of Bins");
			_acceptance_rel_err_hist[i][j]->Write();
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


void Histogram::Polarization(Flags flags_){
	std::cout<<"Polarization\n";
	if(!flags_.Flags::Plot_Polarization()){
		std::cout<<"Not plotting polarization observables\n";
		return;
	}
	char hname[500];
	char xlabel[100];
	char ylabel[100];
	//Get the 7dimensional bins ready
	TH1D_1d_star exp_ch_1d_raw;
	TH1D_1d_star exp_ch_1d_local;
	TH1D_1d_star exp_ch_1d_global;
	TH1D_1d_star exp_ch_1d_local_holes;
	TH1D_1d_star exp_ch_1d_global_holes;
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
			sprintf(dirname,"%s_W|[%.3f,%.3f)",_five_dim_[k],Histogram::W_low(j),Histogram::W_top(j));
			dir_P2[k][j] = dir_P1[k]->mkdir(dirname);
			for(int i=0; i<_Q2_nbins_; i++){
				sprintf(dirname,"%s_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",_five_dim_[k],Histogram::W_low(j),Histogram::W_top(j),Histogram::Q2_low(i),Histogram::Q2_top(i));
				dir_P3[k][j][i] = dir_P2[k][j]->mkdir(dirname);
			}
		}
	}
	std::cout<<"Looping to make the histograms\n";
	for(int Wbin=0; Wbin<_W_nbins_; Wbin++){
		for(int Q2bin=0; Q2bin<_Q2_nbins_; Q2bin++){
			for(int Xij=0; Xij<4; Xij++){
				//std::cout<<"entering directory:" <<Xij <<" " <<Wbin <<" " <<Q2bin <<"\n";
				dir_P3[Xij][Wbin][Q2bin]->cd();
				//std::cout<<"we did it! Now let's do histogram stuff\n";
				//std::cout<<"changing names of labels\n";
				sprintf(xlabel,"%s %s",_five_dim_[4],_dim_units_[4]);
				//sprintf(ylabel,"Diff CS (microbarns/(%s %s))",_dim_units_y_[Xij],_dim_units_y_[4]);
				sprintf(ylabel,"Acc Corr Yield");
				//std::cout<<"now looping through " <<_five_dim_[Xij] <<" bins\n";
				
				for(int Xijbin=0; Xijbin<_n_bins_5d[Xij]; Xijbin++){
					sprintf(hname,"raw_%s_2nd_order_diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_Xij:[%.3f,%.3f)_top:%s_var:%s_AREC:%s",_five_dim_[Xij],Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin),_N_5d_raw_yield[Wbin][Q2bin]->GetAxis(Xij)->GetBinLowEdge(Xijbin+1),_N_5d_raw_yield[Wbin][Q2bin]->GetAxis(Xij)->GetBinUpEdge(Xijbin+1),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str(),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					//std::cout<<"hname is now: " <<hname <<"\n";
					//std::cout<<"setting range for raw\n";
					_N_5d_raw_yield[Wbin][Q2bin]->GetAxis(Xij)->SetRange(Xijbin+1,Xijbin+1);
					//std::cout<<"pushing back raw\n";
					exp_ch_1d_raw.push_back(_N_5d_raw_yield[Wbin][Q2bin]->Projection(4,"E"));
					exp_ch_1d_raw[Xijbin]->SetNameTitle(hname,hname);
					exp_ch_1d_raw[Xijbin]->GetXaxis()->SetTitle(xlabel);
					exp_ch_1d_raw[Xijbin]->GetYaxis()->SetTitle(ylabel);
					//std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<Xijbin <<"\n";
					exp_ch_1d_raw[Xijbin]->Write();
					_N_5d_raw_yield[Wbin][Q2bin]->GetAxis(Xij)->SetRange();

					sprintf(hname,"local_%s_2nd_order_diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_Xij:[%.3f,%.3f)_top:%s_var:%s_AREC:%s",_five_dim_[Xij],Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin),_N_5d_local_holes[Wbin][Q2bin]->GetAxis(Xij)->GetBinLowEdge(Xijbin+1),_N_5d_local_holes[Wbin][Q2bin]->GetAxis(Xij)->GetBinUpEdge(Xijbin+1),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str(),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					//std::cout<<"setting range for local\n";
					_N_5d_local_holes[Wbin][Q2bin]->GetAxis(Xij)->SetRange(Xijbin+1,Xijbin+1);
					exp_ch_1d_local.push_back(_N_5d_local_holes[Wbin][Q2bin]->Projection(4,"E"));
					exp_ch_1d_local[Xijbin]->SetNameTitle(hname,hname);
					exp_ch_1d_local[Xijbin]->GetXaxis()->SetTitle(xlabel);
					exp_ch_1d_local[Xijbin]->GetYaxis()->SetTitle(ylabel);
					//std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<Xijbin <<"\n";
					exp_ch_1d_local[Xijbin]->Write();
					_N_5d_local_holes[Wbin][Q2bin]->GetAxis(Xij)->SetRange();

					sprintf(hname,"global_%s_2nd_order_diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_Xij:[%.3f,%.3f)_top:%s_var:%s_AREC:%s",_five_dim_[Xij],Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin),_N_5d_local_holes[Wbin][Q2bin]->GetAxis(Xij)->GetBinLowEdge(Xijbin+1),_N_5d_local_holes[Wbin][Q2bin]->GetAxis(Xij)->GetBinUpEdge(Xijbin+1),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str(),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					//std::cout<<"setting range for global\n";
					_N_5d_global_holes[Wbin][Q2bin]->GetAxis(Xij)->SetRange(Xijbin+1,Xijbin+1);
					exp_ch_1d_global.push_back(_N_5d_global_holes[Wbin][Q2bin]->Projection(4,"E"));
					exp_ch_1d_global[Xijbin]->SetNameTitle(hname,hname);
					exp_ch_1d_global[Xijbin]->GetXaxis()->SetTitle(xlabel);
					exp_ch_1d_global[Xijbin]->GetYaxis()->SetTitle(ylabel);
					//std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<Xijbin <<"\n";
					exp_ch_1d_global[Xijbin]->Write();
					_N_5d_global_holes[Wbin][Q2bin]->GetAxis(Xij)->SetRange();

					sprintf(hname,"local_holes_%s_2nd_order_diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_Xij:[%.3f,%.3f)_top:%s_var:%s_AREC:%s",_five_dim_[Xij],Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin),_N_5d_local_holes[Wbin][Q2bin]->GetAxis(Xij)->GetBinLowEdge(Xijbin+1),_N_5d_local_holes[Wbin][Q2bin]->GetAxis(Xij)->GetBinUpEdge(Xijbin+1),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str(),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					//std::cout<<"setting range for local holes\n";
					_N_local_holes_5d[Wbin][Q2bin]->GetAxis(Xij)->SetRange(Xijbin+1,Xijbin+1);
					exp_ch_1d_local_holes.push_back(_N_local_holes_5d[Wbin][Q2bin]->Projection(4,"E"));
					exp_ch_1d_local_holes[Xijbin]->SetNameTitle(hname,hname);
					exp_ch_1d_local_holes[Xijbin]->GetXaxis()->SetTitle(xlabel);
					exp_ch_1d_local_holes[Xijbin]->GetYaxis()->SetTitle(ylabel);
					//std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<Xijbin <<"\n";
					exp_ch_1d_local_holes[Xijbin]->Write();
					_N_local_holes_5d[Wbin][Q2bin]->GetAxis(Xij)->SetRange();

					sprintf(hname,"global_holes_%s_2nd_order_diff_W:[%.3f,%.3f)_Q2:[%.2f,%.2f)_Xij:[%.3f,%.3f)_top:%s_var:%s_AREC:%s",_five_dim_[Xij],Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin),_N_5d_local_holes[Wbin][Q2bin]->GetAxis(Xij)->GetBinLowEdge(Xijbin+1),_N_5d_local_holes[Wbin][Q2bin]->GetAxis(Xij)->GetBinUpEdge(Xijbin+1),flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str(),_cut_width_[flags_.Flags::Acc_Rel_Error_Cut()]);
					//std::cout<<"setting range for global holes\n";
					_N_global_holes_5d[Wbin][Q2bin]->GetAxis(Xij)->SetRange(Xijbin+1,Xijbin+1);
					exp_ch_1d_global_holes.push_back(_N_global_holes_5d[Wbin][Q2bin]->Projection(4,"E"));
					exp_ch_1d_global_holes[Xijbin]->SetNameTitle(hname,hname);
					exp_ch_1d_global_holes[Xijbin]->GetXaxis()->SetTitle(xlabel);
					exp_ch_1d_global_holes[Xijbin]->GetYaxis()->SetTitle(ylabel);
					//std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<Xijbin <<"\n";
					exp_ch_1d_global_holes[Xijbin]->Write();
					_N_global_holes_5d[Wbin][Q2bin]->GetAxis(Xij)->SetRange();
				}
				exp_ch_1d_raw.clear();
				exp_ch_1d_local.clear();
				exp_ch_1d_local_holes.clear();
				exp_ch_1d_global.clear();
				exp_ch_1d_global_holes.clear();
			}
		}
	}
	_RootOutputFile->Close();
	std::cout<<"\nCompleted Second Order Differential Cross Section Yields\n";
}

void Histogram::Beam_Spin(Flags flags_){
	std::cout<<"Beam_Spin\n";
	if(!flags_.Flags::Plot_Beam_Spin()){
		std::cout<<"Not plotting Beam Spin\n";
		return;
	}
	char hname[500];
	char xlabel[100];
	sprintf(xlabel,"Phi (deg)");
	char ylabel[100];
	sprintf(ylabel,"Beam Spin Asymmetry (no polarization% mod)");
	char dirname[100];
	/*
	TH1D_1d_star m1_raw_1d;
	TH1D_1d_star m1_local_1d;
	TH1D_1d_star m1_global_1d;
	TH1D_1d_star m1_local_holes_1d;
	TH1D_1d_star m1_global_holes_1d;  
	*/

	TH1D_1d_star raw_pos;
	TH1D_1d_star raw_neg;

	TH1D_1d_star local_pos;
	TH1D_1d_star local_neg;

	TH1D_1d_star global_pos;
	TH1D_1d_star global_neg;
	
	TH1D_1d_star local_holes_pos;
	TH1D_1d_star local_holes_neg;
	TH1D_1d_star global_holes_pos;
	TH1D_1d_star global_holes_neg;

	TH1D_1d_star m2_raw_1d_top;
	TH1D_1d_star m2_local_1d_top;
	TH1D_1d_star m2_global_1d_top;
	TH1D_1d_star m2_local_holes_1d_top;
	TH1D_1d_star m2_global_holes_1d_top;  

	TH1D_1d_star m2_raw_1d_bot;
	TH1D_1d_star m2_local_1d_bot;
	TH1D_1d_star m2_global_1d_bot;
	TH1D_1d_star m2_local_holes_1d_bot;
	TH1D_1d_star m2_global_holes_1d_bot;  
	
	TH1D_1d_star m2_raw_1d;
	TH1D_1d_star m2_local_1d;
	TH1D_1d_star m2_global_1d;
	TH1D_1d_star m2_local_holes_1d;
	TH1D_1d_star m2_global_holes_1d;  

	TH1D_1d_star raw_pos_1d;
	TH1D_1d_star local_pos_1d;
	TH1D_1d_star local_holes_pos_1d;
	TH1D_1d_star global_pos_1d;
	TH1D_1d_star global_holes_pos_1d;
	TH1D_1d_star raw_neg_1d;
	TH1D_1d_star local_neg_1d;
	TH1D_1d_star local_holes_neg_1d;
	TH1D_1d_star global_neg_1d;
	TH1D_1d_star global_holes_neg_1d;

	std::cout<<"Making Output File\n";
	_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	std::cout<<"Made File:" <<flags_.Flags::Output_File().c_str() <<"\n";
	_RootOutputFile->cd();
	TDirectory* dir_B = _RootOutputFile->mkdir("Beam Spin");
	TDirectory* dir_B1[5];
	for(int j=0; j< 5; j++){//Q2
		sprintf(dirname,"Beam_Spin_Q2|[%.2f,%.2f)",Histogram::Q2_low(j),Histogram::Q2_top(j));
		std::cout<<"\t\tMaking Dir: " <<dirname <<"\n";
		dir_B1[j] = dir_B->mkdir(dirname);
	}
	std::cout<<"Entering W-Q2 loop\n";
	for(int Wbin =0; Wbin<29; Wbin++){
		for(int Q2bin=0; Q2bin<5; Q2bin++){
			dir_B1[Q2bin]->cd();
			/*
			sprintf(hname,"raw_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_method1",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m1_raw_1d.push_back(_N_5d_raw_yield_ratio[Wbin][Q2bin]->Projection(4,"E"));
			m1_raw_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m1_raw_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m1_raw_1d[Q2bin]->SetNameTitle(hname,hname);
			m1_raw_1d[Q2bin]->Write();
			*/
			//if(flags_.Flags::Plot_Raw()){
			//std::cout<<" Raw Top Making\n";
			//std::cout<<"yield\n";
			raw_pos.push_back(_N_5d_raw_yield_pos[Wbin][Q2bin]->Projection(4,"E"));
			raw_neg.push_back(_N_5d_raw_yield_neg[Wbin][Q2bin]->Projection(4,"E"));
			sprintf(hname,"raw_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_top",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m2_raw_1d_top.push_back((TH1D*)raw_pos[Q2bin]->Clone());
			//m2_raw_1d_top.push_back(_N_5d_raw_yield_pos[Wbin][Q2bin]->Projection(4,"E"));
			//std::cout<<" Raw Top Change name\n";
			m2_raw_1d_top[Q2bin]->SetNameTitle(hname,hname);
			m2_raw_1d_bot.push_back((TH1D*)raw_pos[Q2bin]->Clone());
			//std::cout<<" Raw Top Subtract Neg\n";
			m2_raw_1d_top[Q2bin]->Add((TH1D*)raw_neg[Q2bin]->Clone(),-1.008);
			m2_raw_1d_top[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_raw_1d_top[Q2bin]->GetYaxis()->SetTitle(ylabel);
			//std::cout<<" Raw Bot Making\n";
			sprintf(hname,"raw_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_bot",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			//m2_raw_1d_bot.push_back(_N_5d_raw_yield_pos[Wbin][Q2bin]->Projection(4,"E"));
			
			m2_raw_1d_bot[Q2bin]->SetNameTitle(hname,hname);
			//std::cout<<" Raw bot Add Neg\n";
			m2_raw_1d_bot[Q2bin]->Add((TH1D*)raw_neg[Q2bin]->Clone(),1.008);
			m2_raw_1d_bot[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_raw_1d_bot[Q2bin]->GetYaxis()->SetTitle(ylabel);
			sprintf(hname,"raw_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m2_raw_1d.push_back(m2_raw_1d_top[Q2bin]);
			m2_raw_1d[Q2bin]->SetNameTitle(hname,hname);
			//std::cout<<" Raw Division\n";
			m2_raw_1d[Q2bin]->Divide(m2_raw_1d_bot[Q2bin]);
			m2_raw_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_raw_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_raw_1d[Q2bin]->Write();
			sprintf(hname,"raw_N_pos_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			//std::cout<<" Raw Pos yield\n";
			raw_pos_1d.push_back((TH1D*)raw_pos[Q2bin]->Clone());
			raw_pos_1d[Q2bin]->SetNameTitle(hname,hname);
			raw_pos_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			raw_pos_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			raw_pos_1d[Q2bin]->Write();
			sprintf(hname,"raw_N_neg_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			//std::cout<<" Raw Neg yield\n";
			raw_neg_1d.push_back((TH1D*)raw_neg[Q2bin]->Clone());
			raw_neg_1d[Q2bin]->SetNameTitle(hname,hname);
			raw_neg_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			raw_neg_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			raw_neg_1d[Q2bin]->Write();
			//}
			/*
			sprintf(hname,"local_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_method1",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m1_local_1d.push_back(_N_5d_local_holes_ratio[Wbin][Q2bin]->Projection(4,"E"));
			m1_local_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m1_local_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m1_local_1d[Q2bin]->SetNameTitle(hname,hname);
			m1_local_1d[Q2bin]->Write();
			*/
			//if(flags_.Flags::Plot_Local()){
			//std::cout<<"Made it to local stuff\n";
			//std::cout<<"yield + local holes\n";
			local_pos.push_back(_N_5d_local_holes_pos[Wbin][Q2bin]->Projection(4,"E"));
			local_neg.push_back(_N_5d_local_holes_neg[Wbin][Q2bin]->Projection(4,"E"));
			//std::cout<<"1\n";
			sprintf(hname,"local_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_top",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m2_local_1d_top.push_back((TH1D*)local_pos[Q2bin]->Clone());
			m2_local_1d_bot.push_back((TH1D*)local_pos[Q2bin]->Clone());
			m2_local_1d_top[Q2bin]->Add((TH1D*)local_neg[Q2bin]->Clone(),-1.008);
			m2_local_1d_top[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_local_1d_top[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_local_1d_top[Q2bin]->SetNameTitle(hname,hname);
			//std::cout<<"2\n";
			sprintf(hname,"local_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_bot",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			//m2_local_1d_bot.push_back(_N_5d_local_holes_pos[Wbin][Q2bin]->Projection(4,"E"));
			m2_local_1d_bot[Q2bin]->Add((TH1D*)local_neg[Q2bin]->Clone(),1.008);
			m2_local_1d_bot[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_local_1d_bot[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_local_1d_bot[Q2bin]->SetNameTitle(hname,hname);
			//std::cout<<"3\n";
			sprintf(hname,"local_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m2_local_1d.push_back(m2_local_1d_top[Q2bin]);
			m2_local_1d[Q2bin]->Divide(m2_local_1d_bot[Q2bin]);
			m2_local_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_local_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_local_1d[Q2bin]->SetNameTitle(hname,hname);
			m2_local_1d[Q2bin]->Write();
			//std::cout<<"4\n";
			sprintf(hname,"local_N_pos_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			local_pos_1d.push_back((TH1D*)local_pos[Q2bin]->Clone());
			local_pos_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			local_pos_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			local_pos_1d[Q2bin]->SetNameTitle(hname,hname);
			local_pos_1d[Q2bin]->Write();
			//std::cout<<"5\n";
			sprintf(hname,"local_N_neg_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			local_neg_1d.push_back((TH1D*)local_neg[Q2bin]->Clone());
			local_neg_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			local_neg_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			local_neg_1d[Q2bin]->SetNameTitle(hname,hname);
			local_neg_1d[Q2bin]->Write();
			//}

			/*
			sprintf(hname,"global_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_method1",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m1_global_1d.push_back(_N_5d_global_holes_ratio[Wbin][Q2bin]->Projection(4,"E"));
			m1_global_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m1_global_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m1_global_1d[Q2bin]->SetNameTitle(hname,hname);
			m1_global_1d[Q2bin]->Write();
			*/
			//if(flags_.Flags::Plot_Global()){
			//std::cout<<"yield + global holes\n";
			global_pos.push_back(_N_5d_global_holes_pos[Wbin][Q2bin]->Projection(4,"E"));
			global_neg.push_back(_N_5d_global_holes_neg[Wbin][Q2bin]->Projection(4,"E"));

			sprintf(hname,"global_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_top",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m2_global_1d_top.push_back((TH1D*)global_pos[Q2bin]->Clone());
			m2_global_1d_bot.push_back((TH1D*)global_pos[Q2bin]->Clone());
			m2_global_1d_top[Q2bin]->Add((TH1D*)global_neg[Q2bin]->Clone(),-1.008);
			m2_global_1d_top[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_global_1d_top[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_global_1d_top[Q2bin]->SetNameTitle(hname,hname);
			sprintf(hname,"global_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_bot",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			//m2_global_1d_bot.push_back(_N_5d_global_holes_pos[Wbin][Q2bin]->Projection(4,"E"));
			m2_global_1d_bot[Q2bin]->Add((TH1D*)global_neg[Q2bin]->Clone(),1.008);
			m2_global_1d_bot[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_global_1d_bot[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_global_1d_bot[Q2bin]->SetNameTitle(hname,hname);
			sprintf(hname,"global_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m2_global_1d.push_back(m2_global_1d_top[Q2bin]);
			m2_global_1d[Q2bin]->Divide(m2_global_1d_bot[Q2bin]);
			m2_global_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_global_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_global_1d[Q2bin]->SetNameTitle(hname,hname);
			m2_global_1d[Q2bin]->Write();
			sprintf(hname,"global_N_pos_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			global_pos_1d.push_back((TH1D*)global_pos[Q2bin]->Clone());
			global_pos_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			global_pos_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			global_pos_1d[Q2bin]->SetNameTitle(hname,hname);
			global_pos_1d[Q2bin]->Write();
			sprintf(hname,"global_N_neg_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			global_neg_1d.push_back((TH1D*)global_neg[Q2bin]->Clone());
			global_neg_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			global_neg_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			global_neg_1d[Q2bin]->SetNameTitle(hname,hname);
			global_neg_1d[Q2bin]->Write();
			//}
			
			/*
			sprintf(hname,"local_holes_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_method1",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m1_local_holes_1d.push_back(_N_local_holes_5d_ratio[Wbin][Q2bin]->Projection(4,"E"));
			m1_local_holes_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m1_local_holes_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m1_local_holes_1d[Q2bin]->SetNameTitle(hname,hname);
			m1_local_holes_1d[Q2bin]->Write();
			*/
			//if(flags_.Flags::Plot_Local_Holes()){
			//std::cout<<"local holes\n";
			local_holes_pos.push_back(_N_local_holes_5d_pos[Wbin][Q2bin]->Projection(4,"E"));
			local_holes_neg.push_back(_N_local_holes_5d_neg[Wbin][Q2bin]->Projection(4,"E"));
			sprintf(hname,"local_holes_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_top",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m2_local_holes_1d_top.push_back((TH1D*)local_holes_pos[Q2bin]->Clone());
			m2_local_holes_1d_bot.push_back((TH1D*)local_holes_pos[Q2bin]->Clone());
			m2_local_holes_1d_top[Q2bin]->Add((TH1D*)local_holes_neg[Q2bin]->Clone(),-1.008);
			m2_local_holes_1d_top[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_local_holes_1d_top[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_local_holes_1d_top[Q2bin]->SetNameTitle(hname,hname);
			m2_local_holes_1d_top[Q2bin]->Write();
			sprintf(hname,"local_holes_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_bot",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			//m2_local_holes_1d_bot.push_back(_N_local_holes_5d_pos[Wbin][Q2bin]->Projection(4,"E"));
			m2_local_holes_1d_bot[Q2bin]->Add((TH1D*)local_holes_neg[Q2bin]->Clone(),1.008);
			m2_local_holes_1d_bot[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_local_holes_1d_bot[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_local_holes_1d_bot[Q2bin]->SetNameTitle(hname,hname);
			m2_local_holes_1d_bot[Q2bin]->Write();
			sprintf(hname,"local_holes_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m2_local_holes_1d.push_back(m2_local_holes_1d_top[Q2bin]);
			m2_local_holes_1d[Q2bin]->Divide(m2_local_holes_1d_bot[Q2bin]);
			m2_local_holes_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_local_holes_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_local_holes_1d[Q2bin]->SetNameTitle(hname,hname);
			m2_local_holes_1d[Q2bin]->Write();
			sprintf(hname,"local_holes_N_pos_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			local_holes_pos_1d.push_back((TH1D*)local_holes_pos[Q2bin]->Clone());
			local_holes_pos_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			local_holes_pos_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			local_holes_pos_1d[Q2bin]->SetNameTitle(hname,hname);
			local_holes_pos_1d[Q2bin]->Write();
			sprintf(hname,"local_holes_N_neg_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			local_holes_neg_1d.push_back((TH1D*)local_holes_neg[Q2bin]->Clone());
			local_holes_neg_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			local_holes_neg_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			local_holes_neg_1d[Q2bin]->SetNameTitle(hname,hname);
			local_holes_neg_1d[Q2bin]->Write();
			//}

			/*
			sprintf(hname,"global_holes_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_method1",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m1_global_holes_1d.push_back(_N_global_holes_5d_ratio[Wbin][Q2bin]->Projection(4,"E"));
			m1_global_holes_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m1_global_holes_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m1_global_holes_1d[Q2bin]->SetNameTitle(hname,hname);
			m1_global_holes_1d[Q2bin]->Write();
			*/
			//if(flags_.Flags::Plot_Global_Holes()){
			//std::cout<<"global holes\n";
			global_holes_pos.push_back(_N_global_holes_5d_pos[Wbin][Q2bin]->Projection(4,"E"));
			global_holes_neg.push_back(_N_global_holes_5d_pos[Wbin][Q2bin]->Projection(4,"E"));
			sprintf(hname,"global_holes_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_top",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m2_global_holes_1d_top.push_back((TH1D*)global_holes_pos[Q2bin]->Clone());
			m2_global_holes_1d_bot.push_back((TH1D*)global_holes_pos[Q2bin]->Clone());
			m2_global_holes_1d_top[Q2bin]->Add((TH1D*)global_holes_neg[Q2bin]->Clone(),-1.008);
			m2_global_holes_1d_top[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_global_holes_1d_top[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_global_holes_1d_top[Q2bin]->SetNameTitle(hname,hname);
			//m2_global_holes_1d_top[Q2bin]->Write();
			sprintf(hname,"global_holes_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)_bot",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			//m2_global_holes_1d_bot.push_back(_N_global_holes_5d_pos[Wbin][Q2bin]->Projection(4,"E"));
			m2_global_holes_1d_bot[Q2bin]->Add((TH1D*)global_holes_neg[Q2bin]->Clone(),1.008);
			m2_global_holes_1d_bot[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_global_holes_1d_bot[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_global_holes_1d_bot[Q2bin]->SetNameTitle(hname,hname);
			//m2_global_holes_1d_bot[Q2bin]->Write();
			sprintf(hname,"global_holes_Beam_Spin_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			m2_global_holes_1d.push_back(m2_global_holes_1d_top[Q2bin]);
			m2_global_holes_1d[Q2bin]->Divide(m2_global_holes_1d_bot[Q2bin]);
			m2_global_holes_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			m2_global_holes_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			m2_global_holes_1d[Q2bin]->SetNameTitle(hname,hname);
			m2_global_holes_1d[Q2bin]->Write();
			sprintf(hname,"global_holes_N_pos_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			global_holes_pos_1d.push_back((TH1D*)global_holes_pos[Q2bin]->Clone());
			global_holes_pos_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			global_holes_pos_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			global_holes_pos_1d[Q2bin]->SetNameTitle(hname,hname);
			global_holes_pos_1d[Q2bin]->Write();
			sprintf(hname,"global_holes_N_neg_W|[%.3f,%.3f)_Q2|[%.2f,%.2f)",Histogram::W_low(Wbin),Histogram::W_top(Wbin),Histogram::Q2_low(Q2bin),Histogram::Q2_top(Q2bin));
			global_holes_neg_1d.push_back((TH1D*)global_holes_neg[Q2bin]->Clone());
			global_holes_neg_1d[Q2bin]->GetXaxis()->SetTitle(xlabel);
			global_holes_neg_1d[Q2bin]->GetYaxis()->SetTitle(ylabel);
			global_holes_neg_1d[Q2bin]->SetNameTitle(hname,hname);
			global_holes_neg_1d[Q2bin]->Write();
			//}
		}
		raw_pos.clear();
		raw_neg.clear();

		local_pos.clear();
		local_neg.clear();

		local_holes_pos.clear();
		local_holes_neg.clear();

		global_pos.clear();
		global_neg.clear();

		global_holes_pos.clear();
		global_holes_neg.clear();
		//m1_raw_1d.clear();
		m2_raw_1d.clear();
		//m1_local_1d.clear();
		m2_local_1d.clear();
		//m1_local_holes_1d.clear();
		m2_local_holes_1d.clear();
		//m1_global_1d.clear();
		m2_global_1d.clear();
		//m1_global_holes_1d.clear();
		m2_global_holes_1d.clear();
		raw_pos_1d.clear();
		local_pos_1d.clear();
		local_holes_pos_1d.clear();
		global_pos_1d.clear();
		global_holes_pos_1d.clear();
		raw_neg_1d.clear();
		local_neg_1d.clear();
		local_holes_neg_1d.clear();
		global_neg_1d.clear();
		global_holes_neg_1d.clear();
	}

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

void Histogram::Make_Acceptance_Rel_Hist(Flags flags_){
	char hname[100];

	for(int a=0; a<_W_nbins_; a++){
		for(int b=0; b<_Q2_nbins_; b++){
			sprintf(hname,"Acceptance_Relative_Error_%s_%s_W:%.3f-%.3f_Q2:%.2f-%.2f",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],Histogram::W_low(a),Histogram::W_top(a),Histogram::Q2_low(b),Histogram::Q2_top(b));
			_acceptance_rel_err_hist[a][b]= new TH1D(hname,hname,200,0.0,1.0);
		}
	}
}

void Histogram::Fill_Acceptance_Rel_Hist(){
	for(int a=0; a<_W_nbins_; a++){
		for(int b=0; b<_Q2_nbins_; b++){
			//Int_t* coord = new Int_t[_acceptance_5d[i][j]->GetNbins()];
			for(Long64_t i =0; i<_acceptance_5d[a][b]->GetNbins(); i++ ){
				if(_acceptance_5d[a][b]->GetBinContent(i)>0.0){
					_acceptance_rel_err_hist[a][b]->Fill(_acceptance_5d[a][b]->GetBinError(i)/_acceptance_5d[a][b]->GetBinContent(i));
				}
			}
		}
	}
}