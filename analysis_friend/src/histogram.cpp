#include "histogram.hpp"

Histogram::Histogram( TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_){
	//Make output file eventually RootOutputFile = fun::make_file
	std::cout<<"Name_File\n";
	_localized_hole_filling = flags_.Flags::Localized_Holes();
	//fun::Name_File(output_file_);
	Histogram::Name_Output(flags_);
	Histogram::Extract_7d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,flags_);//Gets the 7d histograms off the root file
	//Histogram::Extract_Bin_Info(flags_);//Extracts binning information about the 7d histograms just extracted //moved to inside Extract Histograms
	//Histogram::Skeleton_5D(flags_);//Creates empty vector arrays to map onto in the next step
	if(!flags_.Flags::Plot_Localized_Holes()){
		Histogram::Sparse_7to5(flags_);//Converts all the 7d histograms into usable 5d histograms for analysis
	}
	//Histogram::Single_Differential(flags_);//Single Differential Histograms
	//Histogram::Polarization_Observables(flags_);//Polarization Observables
	//Histogram::Beam_Spin(flags_);//Beam Spin Asymmetry
	//Histogram::Convert_to_Cross(flags_);
	//Histogram::Calc_Cross_Section(flags_);
	//Histogram::Acceptance_Errors(flags_);
	//Histogram::Dumb_Histograms(flags_);
	Histogram::Make_Histograms(flags_);
	//Histogram::Fill_Histograms(flags_);
	//Histogram::Write_Histograms(flags_);
}

Histogram::Histogram( TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_){
	//Make output file eventually RootOutputFile = fun::make_file
	std::cout<<"Name_File\n";
	_localized_hole_filling = flags_.Flags::Localized_Holes();
	//fun::Name_File(output_file_);
	Histogram::Name_Output(flags_);
	Histogram::Extract_7d_Histograms(exp_tree_,sim_tree_,empty_tree_,nr_sim_tree_,holes_,flags_);//Gets the 7d histograms off the root file
	//Histogram::Extract_Bin_Info(flags_);//Extracts binning information about the 7d histograms just extracted //moved to inside Extract Histograms
	//Histogram::Skeleton_5D(flags_);//Creates empty vector arrays to map onto in the next step
	Histogram::Sparse_7to5(flags_);//Converts all the 7d histograms into usable 5d histograms for analysis
	//Histogram::Single_Differential(flags_);//Single Differential Histograms
	//Histogram::Polarization_Observables(flags_);//Polarization Observables
	//Histogram::Beam_Spin(flags_);//Beam Spin Asymmetry
	//Histogram::Convert_to_Cross(flags_);
	//Histogram::Calc_Cross_Section(flags_);
	//Histogram::Acceptance_Errors(flags_);
	//Histogram::Dumb_Histograms(flags_);
	Histogram::Make_Histograms(flags_);
	//Histogram::Fill_Histograms(flags_);
	//Histogram::Write_Histograms(flags_);
}

void Histogram::Make_Histograms(Flags flags_){
	//_RootOutputFile=Histogram::Name_Output(flags_);
	
	std::cout<<"Making Histograms\n";
	//_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	Histogram::Calc_Error_R(flags_);
	Histogram::Single_Differential(flags_);//Single Differential Histograms
	Histogram::Polarization_Observables(flags_);//Polarization Observables
	Histogram::Beam_Spin(flags_);//Beam Spin Asymmetry
	Histogram::Make_WQ2(flags_);
	Histogram::Localized_Holes(flags_,1,8);//Given last test in May 2023 the max only needs to be 6
	//Histogram::Make_Acceptance(flags_);
	//Histogram::Make_Single_Diff(flags_);
	//Histogram::Make_Polarization(flags_);
	//Histogram::Make_Acceptance_Statistics(flags_); //Cannot perform here. Require rootfile level for sim recon data
	_RootOutputFile->Close();
	//Histogram::Make_Error_Hists(flags_);
}

void Histogram::Fill_Histograms(Flags flags_){
	//Histogram::Fill_Error_Hists(flags_);
}
	
void Histogram::Write_Histograms(Flags flags_){
	std::cout<<"Writing Histograms\n";
	Histogram::Acceptance_Histograms(flags_);
	Histogram::Hole_Histograms(flags_);
	//Histogram::Write_WQ2(flags_); //Already written when making
	//Histogram::Write_Acceptance(flags_); //Already written when making
	//Histogram::Write_Single_Diff(flags_);
	
	//Histogram::Write_Error_Hists(flags_);
	//Histogram::Write_Polarization(flags_);
	if(!flags_.Flags::Plot_Single_Diff() && !flags_.Flags::Plot_Pol()){	
		_RootOutputFile->Close();
	}
	std::cout<<"Output File: " <<flags_.Flags::Output_File() <<"\n";
	if(flags_.Flags::Plot_Single_Diff()){
		std::cout<<"Output File: " <<flags_.Flags::Output_File() <<"\n";
	}
	if(flags_.Flags::Plot_Pol()){
		std::cout<<"Output File: " <<flags_.Flags::Output_File() <<"\n";
	}
}

void Histogram::Name_Output(Flags flags_){
//std::shared_ptr<TFile> Histogram::Name_Output(Flags flags_){
	//_RootOutputFile=TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	
	//if(!flags_.Flags::Plot_Single_Diff() && !flags_.Flags::Plot_Pol()){	
		//_RootOutputFile=TFile::Open(flags_.Flags::Output_File().c_str(),"RECREATE");
		std::cout<<"Naming File1: " <<flags_.Flags::Output_File().c_str() <<"\n";
		//_RootOutputFile=TFile::Open(flags_.Flags::Output_File().c_str(),"RECREATE");
		//_RootOutputFile=std::make_shared<TFile>(flags_.Flags::Output_File().c_str(),"RECREATE");
	//}
	/*
	if(flags_.Flags::Plot_Single_Diff()){
		std::cout<<"Naming Single Diff File: " <<flags_.Flags::Output_File(1).c_str() <<"\n";
		_RootOutputFile=TFile::Open(flags_.Flags::Output_File(1).c_str(),"RECREATE");
	}
	if(flags_.Flags::Plot_Pol()){
		std::cout<<"Naming Polarization File: " <<flags_.Flags::Output_File(2).c_str() <<"\n";
		_RootOutputFile=TFile::Open(flags_.Flags::Output_File(2).c_str(),"RECREATE");
	}*/
}

void Histogram::Extract_7d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_){
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
	if(flags_.Flags::Helicity()){
		
		sprintf(hname,"%s_%s_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting Exp THnSparse Pos" <<hname <<"\n";
		_exp_data_7d_pos = (THnSparseD *)exp_tree_->Get(hname);
		_empty_7d_pos = (THnSparseD *)empty_tree_->Get(hname);
		std::cout<<"\tPositive Helicity Yields| \n\t\tfilled: " <<fun::nSparseIntegral(_exp_data_7d_pos) <<"\n\t\tempty: " <<fun::nSparseIntegral(_empty_7d_pos) <<"\n";
		sprintf(hname,"%s_%s_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting Exp THnSparse Neg" <<hname <<"\n";
		_exp_data_7d_neg = (THnSparseD *)exp_tree_->Get(hname);
		_empty_7d_neg = (THnSparseD *)empty_tree_->Get(hname);
		std::cout<<"\tNegative Helicity Yields| \n\t\tfilled: " <<fun::nSparseIntegral(_exp_data_7d_neg) <<"\n\t\tempty: " <<fun::nSparseIntegral(_empty_7d_neg) <<"\n";
		sprintf(hname,"%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting full Exp THnSparse " <<hname <<"\n";
		_exp_data_7d = (THnSparseD *)exp_tree_->Get(hname);//(THnSparseD *)_exp_data_7d_pos->Clone();
		_empty_7d = (THnSparseD *)empty_tree_->Get(hname);//(THnSparseD *)_empty_7d_pos->Clone();
		std::cout<<"\tTotal Helicity Yields| \n\t\tfilled: " <<fun::nSparseIntegral(_exp_data_7d) <<"\n\t\tempty: " <<fun::nSparseIntegral(_empty_7d) <<"\n";
	}else{
		sprintf(hname,"%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting Exp THnSparse " <<hname <<"\n";
		_exp_data_7d = (THnSparseD *)exp_tree_->Get(hname);
		std::cout<<"Getting Exp Empty THnSparse " <<hname <<"\n";
		_empty_7d = (THnSparseD *)empty_tree_->Get(hname);
	}
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
	if(flags_.Flags::Helicity()){
		_N_pos=(THnSparseD*)_exp_data_7d_pos->Clone();
		_N_pos->Add(_empty_7d_pos,-flags_.Flags::Qr());
		_N_pos->Divide(_acceptance_7d);
		_N_neg=(THnSparseD*)_exp_data_7d_neg->Clone();
		_N_neg->Add(_empty_7d_neg,-flags_.Flags::Qr());
		_N_neg->Divide(_acceptance_7d);
	}
	
	std::cout<<"Calculating Simulated 7d Holes\n";
	//Simulated holes will be opposite in yield from what it should be, due to order of operations
	std::cout<<"\tCloning simulation reconstructed\n";
	_sim_holes_tmp_7d=(THnSparseD*)_sim_data_7d->Clone();
	std::cout<<"\tDividing by acceptance\n";
	_sim_holes_tmp_7d->Divide(_acceptance_7d);
	std::cout<<"\tDifference from thrown\n";//This step takes a lot longer than I would have expected
	_sim_holes_7d=(THnSparseD*)_thrown_7d->Clone();
	_sim_holes_7d->Add(_sim_holes_tmp_7d,-1.0);
	std::cout<<"Making Experimental 7d Holes\n";
	
	//The world isn't ready for localized hole filling here yet 3/14/23
	if(_localized_hole_filling){
		if(flags_.Flags::Plot_Localized_Holes()){
			//Histogram::Localized_Holes(flags_,1,8);
		}else{
			//_N_holes=_exp_holes_7d=fun::Localized_Holes_5d_for_7d(_exp_data_7d,_sim_data_7d,_sim_holes_7d,_n_bins_7d);
			//_N->Add(_N_holes,1.0);
		}
	}else{
		_N_holes=(THnSparseD*)_sim_holes_7d->Clone();
		//This will be scaled properly at the 7to5 stage
	}
	
	std::cout<<"First bin content of _N: " <<_N->GetBinContent(1) <<"\n";
	if(flags_.Flags::Helicity()){
		if(_localized_hole_filling){
			if(flags_.Flags::Plot_Localized_Holes()){

			}else{
				_N_pos->Add(_N_holes_pos,1.0);
				_N_neg->Add(_N_holes_neg,1.0);
			}
			
		}else{
			_N_holes_pos=(THnSparseD*)_sim_holes_7d->Clone();
			_N_holes_neg=(THnSparseD*)_sim_holes_7d->Clone();
			//These will be scaled properly at the 7to5 stage
		}
	}
}

void Histogram::Extract_7d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_){
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
	if(flags_.Flags::Helicity()){
		
		sprintf(hname,"%s_%s_pos",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting Exp THnSparse Pos" <<hname <<"\n";
		_exp_data_7d_pos = (THnSparseD *)exp_tree_->Get(hname);
		_empty_7d_pos = (THnSparseD *)empty_tree_->Get(hname);
		std::cout<<"\tPositive Helicity Yields| \n\t\tfilled: " <<fun::nSparseIntegral(_exp_data_7d_pos) <<"\n\t\tempty: " <<fun::nSparseIntegral(_empty_7d_pos) <<"\n";
		sprintf(hname,"%s_%s_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting Exp THnSparse Neg" <<hname <<"\n";
		_exp_data_7d_neg = (THnSparseD *)exp_tree_->Get(hname);
		_empty_7d_neg = (THnSparseD *)empty_tree_->Get(hname);
		std::cout<<"\tNegative Helicity Yields| \n\t\tfilled: " <<fun::nSparseIntegral(_exp_data_7d_neg) <<"\n\t\tempty: " <<fun::nSparseIntegral(_empty_7d_neg) <<"\n";
		sprintf(hname,"%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting full Exp THnSparse " <<hname <<"\n";
		_exp_data_7d = (THnSparseD *)exp_tree_->Get(hname);//(THnSparseD *)_exp_data_7d_pos->Clone();
		_empty_7d = (THnSparseD *)empty_tree_->Get(hname);//(THnSparseD *)_empty_7d_pos->Clone();
		std::cout<<"\tTotal Helicity Yields| \n\t\tfilled: " <<fun::nSparseIntegral(_exp_data_7d) <<"\n\t\tempty: " <<fun::nSparseIntegral(_empty_7d) <<"\n";
	}else{
		sprintf(hname,"%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting Exp THnSparse " <<hname <<"\n";
		_exp_data_7d = (THnSparseD *)exp_tree_->Get(hname);
		std::cout<<"Getting Exp Empty THnSparse " <<hname <<"\n";
		_empty_7d = (THnSparseD *)empty_tree_->Get(hname);
	}
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
	if(flags_.Flags::Helicity()){
		_N_pos=(THnSparseD*)_exp_data_7d_pos->Clone();
		_N_pos->Add(_empty_7d_pos,-flags_.Flags::Qr());
		_N_pos->Divide(_acceptance_7d);
		_N_neg=(THnSparseD*)_exp_data_7d_neg->Clone();
		_N_neg->Add(_empty_7d_neg,-flags_.Flags::Qr());
		_N_neg->Divide(_acceptance_7d);
	}
	std::cout<<"Adding Holes\n";
	sprintf(hname,"Localized_Holes_50");
	_N_holes = (THnSparseD *)holes_->Get(hname);
	_N->Add(_N_holes);
	
	std::cout<<"First bin content of _N: " <<_N->GetBinContent(1) <<"\n";
	if(flags_.Flags::Helicity()){
		sprintf(hname,"Localized_Holes_50_pos");
		_N_holes_pos = (THnSparseD *)holes_->Get(hname);
		_N_pos->Add(_N_holes_pos);
		sprintf(hname,"Localized_Holes_50_neg");
		_N_holes_neg = (THnSparseD *)holes_->Get(hname);
		_N_neg->Add(_N_holes_neg);
	}
}

void Histogram::Acceptance(){
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
}

void Histogram::Rad_Corr(){
	std::cout<<"Calculating Radiative Effects\n";
	_rad_corr = _thrown_7d_no_rad->Projection(1,0);//First make it the projections
	_rad_corr->SetNameTitle("Rad_Corr","Rad_Corr");
	_rad_corr->Divide(_thrown_7d->Projection(1,0));//Then divide by the according projection
	std::cout<<"Thrown 7d Integral: " <<fun::Sparse_Integral(_thrown_7d) <<" Thrown nr 7d Integral: " <<fun::Sparse_Integral(_thrown_7d_no_rad) <<"\n";
	//_rad_corr->Scale(fun::nSparseIntegral(_thrown_7d)/fun::nSparseIntegral(_thrown_7d_no_rad));//Another place for localized hole filling?
	_rad_corr_mu = fun::Sparse_Integral(_thrown_7d)/fun::Sparse_Integral(_thrown_7d_no_rad);
	_rad_corr->Scale(_rad_corr_mu);//Another place for localized hole filling?
	double_1d corr_tmp; 
	for(int i=0; i<_n_bins_7d[0]; i++){
		for(int j=0; j<_n_bins_7d[1]; j++){
			corr_tmp.push_back(_rad_corr->GetBinContent(i+1,j+1));
		}
		_rad_corr_array.push_back(corr_tmp);
		corr_tmp.clear();
	}
}

void Histogram::Make_N_7d(Flags flags_){
	if(flags_.Flags::Plot_Localized_Holes()){
		return;
	}
	if(flags_.Flags::Run()==2){//This needs a lot to be done a lot of other places

		_N=(THnSparseD*)_exp_data_7d->Clone();
		_N->Add(_empty_7d,-flags_.Flags::Qr());//Empty target subtraction
		_N->Divide(_acceptance_7d);//Acceptance Correction
		if(flags_.Flags::Helicity()){
			_N_pos=(THnSparseD*)_exp_data_7d_pos->Clone();
			_N_pos->Add(_empty_7d_pos,-flags_.Flags::Qr());
			_N_pos->Divide(_acceptance_7d);
			_N_neg=(THnSparseD*)_exp_data_7d_neg->Clone();
			_N_neg->Add(_empty_7d_neg,-flags_.Flags::Qr());
			_N_neg->Divide(_acceptance_7d);
		}
		_N2=(THnSparseD*)_exp_data2_7d->Clone();
		_N2->Add(_empty2_7d,-flags_.Flags::Qr());//Empty target subtraction
		_N2->Divide(_acceptance2_7d);//Acceptance Correction
		if(flags_.Flags::Helicity()){
			_N2_pos=(THnSparseD*)_exp_data2_7d_pos->Clone();
			_N2_pos->Add(_empty2_7d_pos,-flags_.Flags::Qr());
			_N2_pos->Divide(_acceptance2_7d);
			_N2_neg=(THnSparseD*)_exp_data2_7d_neg->Clone();
			_N2_neg->Add(_empty2_7d_neg,-flags_.Flags::Qr());
			_N2_neg->Divide(_acceptance2_7d);
		}
	}else{
		_N=(THnSparseD*)_exp_data_7d->Clone();
		_N->Add(_empty_7d,-flags_.Flags::Qr());//Empty target subtraction
		_N->Divide(_acceptance_7d);//Acceptance Correction
		if(flags_.Flags::Helicity()){
			_N_pos=(THnSparseD*)_exp_data_7d_pos->Clone();
			_N_pos->Add(_empty_7d_pos,-flags_.Flags::Qr());
			_N_pos->Divide(_acceptance_7d);
			_N_neg=(THnSparseD*)_exp_data_7d_neg->Clone();
			_N_neg->Add(_empty_7d_neg,-flags_.Flags::Qr());
			_N_neg->Divide(_acceptance_7d);
		}
		if(flags_.Flags::Localized_Holes()){
			Histogram::Localized_Holes(flags_,1,8);
			//_N_holes=_exp_holes_7d=fun::Localized_Holes_5d_for_7d(_exp_data_7d,_sim_data_7d,_sim_holes_7d,_n_bins_7d);
			_N->Add(_N_holes,1.0);
			if(flags_.Flags::Helicity()){
				_N_pos->Add(_N_holes_pos,1.0);
				_N_neg->Add(_N_holes_neg,1.0);
			}
		}else{
			std::cout<<"Calculating Simulated 7d Holes\n";
			//Simulated holes will be opposite in yield from what it should be, due to order of operations
			std::cout<<"\tCloning simulation reconstructed\n";
			_sim_holes_tmp_7d=(THnSparseD*)_sim_data_7d->Clone();
			std::cout<<"\tDividing by acceptance\n";
			_sim_holes_tmp_7d->Divide(_acceptance_7d);
			std::cout<<"\tDifference from thrown\n";//This step takes a lot longer than I would have expected
			_sim_holes_7d=(THnSparseD*)_thrown_7d->Clone();
			_sim_holes_7d->Add(_sim_holes_tmp_7d,-1.0);
			std::cout<<"Making Experimental 7d Holes\n";
			_N_holes=(THnSparseD*)_sim_holes_7d->Clone();//This will be scaled in 7to5
		}
	}
	
}

void Histogram::Localized_Holes(Flags flags_, int min_dist_ = 1, int max_dist_ = -1){
	if(!flags_.Flags::Plot_Localized_Holes()){
		return;
	}
	_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	_RootOutputFile->cd();
	_scale_exp_7d;
	_scale_sim_7d;
	_scale_7d;
	_sim_holes_7d;
	_N_holes;
	std::vector<long> space_dims;
	for(int i=0; i<_n_bins_7d.size(); i++){
		//space_dims.push_back(_n_bins_7d[1]);
      	space_dims.push_back(_n_bins_7d[i]);
	}
	_localized_hole_filling = false;
	CartesianGenerator cart(space_dims);
   	int bin[_n_bins_7d.size()];
   	int bin2[_n_bins_7d.size()];
	bool look_further_all = false;
	bool look_further_pos = false;
	bool look_further_neg = false;
	int dist = 0; 
	int bin_low[5];
	int bin_top[5];
   	std::vector<std::vector<int>> surr_bins; 
	_scale_exp_7d=(THnSparseD*)_sim_holes_7d->Clone();
	_scale_sim_7d=(THnSparseD*)_sim_holes_7d->Clone();
	_scale_exp_7d_pos=(THnSparseD*)_sim_holes_7d->Clone();
	_scale_exp_7d_neg=(THnSparseD*)_sim_holes_7d->Clone();
	_scale_sim_7d_pos=(THnSparseD*)_sim_holes_7d->Clone();
	_scale_sim_7d_neg=(THnSparseD*)_sim_holes_7d->Clone();
	int proj_bins[4] = {2,3,4,5};

	long dist_dist[3][14] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0}}; 

	THnSparseD* scale_exp;
	THnSparseD* scale_sim;
	THnSparseD* scale_exp_pos;
	//THnSparseD* scale_sim_pos;
	THnSparseD* scale_exp_neg;

	double loc_exp;
	double loc_sim;
	double loc_exp_pos;
	double loc_sim_pos;
	double loc_exp_neg;
	double loc_sim_neg;

	long total_ev = _sim_holes_7d->GetNbins();
	long curr_ev = 0;
	//THnSparseD* scale_sim_neg;
   	while(cart.GetNextCombination()){
      	dist = 0;
      	for(int j = 0; j<_n_bins_7d.size(); j++){
    		bin[j] = cart[j]+1;
			_sim_data_7d->GetAxis(j)->SetRange();
			if(flags_.Flags::Helicity()){
				_exp_data_7d_pos->GetAxis(j)->SetRange();
				_exp_data_7d_neg->GetAxis(j)->SetRange();
			}else{
				_exp_data_7d->GetAxis(j)->SetRange();
			}
      	}
		if(_sim_holes_7d->GetBinContent(bin)>0.0){
			curr_ev ++; 
			if((curr_ev-1)%(total_ev/1000) == 0){
				std::cout<<"\r" <<"\t" <<(1000*curr_ev/total_ev) <<"/1000"  <<std::flush ;
			}
			if(curr_ev%(total_ev/200) == 0){
				std::cout<<"\n";
				for(int i=0; i<3; i++){
					for(int j=0; j<14; j++){
						std::cout<<dist_dist[i][j] <<" ";
					}
					std::cout<<"\n";
				}
				std::cout<<"\n";
			}
		
			if(flags_.Flags::Helicity()){
				look_further_pos = true;
				look_further_neg = true;
				_scale_exp_7d_pos->SetBinContent(bin,0.0);
				_scale_sim_7d_pos->SetBinContent(bin,0.0);
				_scale_exp_7d_neg->SetBinContent(bin,0.0);
				_scale_sim_7d_neg->SetBinContent(bin,0.0);
			}else{
				look_further_all = true;
				_scale_exp_7d->SetBinContent(bin,0.0);
				_scale_sim_7d->SetBinContent(bin,0.0);
			}
			while(look_further_all || look_further_all || look_further_all){
				//std::cout<<"\nLooking further\n";
				dist++;
				if(dist >= min_dist_){
					//std::cout<<"Setting ranges\n";
					
					_sim_data_7d->GetAxis(0)->SetRange(bin[0],bin[0]);
					_sim_data_7d->GetAxis(1)->SetRange(bin[1],bin[1]);
					if(flags_.Flags::Helicity()){
						_exp_data_7d_pos->GetAxis(0)->SetRange(bin[0],bin[0]);
						_exp_data_7d_pos->GetAxis(1)->SetRange(bin[1],bin[1]);
						_exp_data_7d_neg->GetAxis(0)->SetRange(bin[0],bin[0]);
						_exp_data_7d_neg->GetAxis(1)->SetRange(bin[1],bin[1]);
					}else{
						_exp_data_7d->GetAxis(0)->SetRange(bin[0],bin[0]);
						_exp_data_7d->GetAxis(1)->SetRange(bin[1],bin[1]);
					}
					for(int i=2; i<7; i++){
						if((bin[i]+dist) > _n_bins_7d[i]){
							bin_top[i-2] = _n_bins_7d[i];
						}else{
							bin_top[i-2] = (bin[i]+dist);
						}
						if((bin[i]-dist) < 1){
							bin_low[i-2] = 1;
						}else{
							bin_low[i-2] = (bin[i]-dist);
						}
					}
					for(int i=0; i<5; i++){
						_sim_data_7d->GetAxis(i+2)->SetRange(bin_low[i],bin_top[i]);
						if(flags_.Flags::Helicity()){
							_exp_data_7d_pos->GetAxis(i+2)->SetRange(bin_low[i],bin_top[i]);
							_exp_data_7d_neg->GetAxis(i+2)->SetRange(bin_low[i],bin_top[i]);
						}else{
							_exp_data_7d->GetAxis(i+2)->SetRange(bin_low[i],bin_top[i]);
						}
					}
					//std::cout<<"Making Projections\nscale_exp\n";
					
					//std::cout<<"scale_sim\n";
					scale_sim = (THnSparseD*)_sim_data_7d->Projection(4,proj_bins,"E");
					if(flags_.Flags::Helicity()){
						//std::cout<<"scale_pos\n";
						scale_exp_pos = (THnSparseD*)_exp_data_7d_pos->Projection(4,proj_bins,"E");
						//std::cout<<"scale_neg\n";
						scale_exp_neg = (THnSparseD*)_exp_data_7d_neg->Projection(4,proj_bins,"E");
					}else{
						scale_exp = (THnSparseD*)_exp_data_7d->Projection(4,proj_bins,"E");
					}
					//std::cout<<"Integral of post sim_7d " <<fun::Sparse_Integral(_sim_data_7d) <<"\n";
					//std::cout<<"Integral of post proj sim_7d " <<fun::Sparse_Integral(scale_sim) <<"\n";
					//std::cout<<"Num bins of post sim_7d " <<_sim_data_7d->GetNbins() <<"\n";
					//std::cout<<"Num bins of post proj sim_7d " <<scale_sim->GetNbins() <<"\n";
					
					//std::cout<<"continue to helicity\n";
					if(flags_.Flags::Helicity()){
						if(look_further_pos){
							//std::cout<<"Look further pos\n";
							loc_exp_pos=fun::Sparse_Integral(scale_exp_pos);
							loc_sim_pos=fun::Sparse_Integral(scale_sim);
							if(loc_exp_pos>0.0){
								//std::cout<<"Integrating exp pos\n";
								_scale_exp_7d_pos->SetBinContent(bin,loc_exp_pos);
							}
							if(loc_sim_pos>0.0){
								//std::cout<<"Integrating sim pos\n";
								_scale_sim_7d_pos->SetBinContent(bin,loc_sim_pos);
							}
						}
						if(look_further_neg){
							//std::cout<<"Look further neg\n";
							loc_exp_neg=fun::Sparse_Integral(scale_exp_neg);
							//std::cout<<"Look further sim neg\n";
							loc_sim_neg=fun::Sparse_Integral(scale_sim);
							if(loc_exp_neg>0.0){
								//std::cout<<"Integrating exp neg\n";
								_scale_exp_7d_neg->SetBinContent(bin,loc_exp_neg);
							}
							if(loc_sim_neg>0.0){
								//std::cout<<"Integrating sim neg\n";
								_scale_sim_7d_neg->SetBinContent(bin,loc_sim_neg);
							}
						}
					}else{
						if(look_further_all){
							//std::cout<<"Integrating exp\n";
							loc_exp = fun::Sparse_Integral(scale_exp);
							//std::cout<<"Integrating sim\n";
							loc_sim = fun::Sparse_Integral(scale_sim);
							//std::cout<<"next\n";
							//std::cout<<"\tLocalized Hole Filling at ";
							//for(int k=0; k<_n_bins_7d.size(); k++ ){
								//std::cout<<bin[k] <<" ";
							//}
							//std::cout<<"with dist=" <<dist <<"\n\t\tscale_exp: " <<fun::Sparse_Integral(scale_exp);
							//std::cout<<"with Integral function: " <<_exp_data_7d->ComputeIntegral() <<"\n";
							if(loc_exp>0.0){
								//std::cout<<"Setting exp\n";
								_scale_exp_7d->SetBinContent(bin,loc_exp);
							}
							//std::cout<<" and what bin is now: " <<_scale_exp_7d->GetBinContent(bin) <<"\n";
							//std::cout<<"\n\t\tscale_sim: " <<fun::Sparse_Integral(scale_sim) <<"\n";
							if(loc_sim>0.0){
								//std::cout<<"Setting Sim\n";
								_scale_sim_7d->SetBinContent(bin,loc_sim);
							}
						}
					}
					if(flags_.Flags::Helicity()){
						if(loc_exp_pos>0.0 && loc_sim_pos>0.0){
							//std::cout<<"non-zero result for both positive. Filling scales\n";
							look_further_pos=false;
							_scale_exp_7d_pos->SetBinError2(_scale_exp_7d_pos->GetBin(bin),fun::Sparse_Integral_Error2(scale_exp_pos));
							_scale_sim_7d_pos->SetBinError2(_scale_sim_7d_pos->GetBin(bin),fun::Sparse_Integral_Error2(scale_sim));
							dist_dist[1][dist-1]++;
						}
						if(loc_exp_neg>0.0 && loc_sim_neg>0.0){
							//std::cout<<"non-zero result for both negative. Filling scales\n";
							look_further_neg=false;
							_scale_exp_7d_neg->SetBinError2(_scale_exp_7d_neg->GetBin(bin),fun::Sparse_Integral_Error2(scale_exp_neg));
							_scale_sim_7d_neg->SetBinError2(_scale_sim_7d_neg->GetBin(bin),fun::Sparse_Integral_Error2(scale_sim));
							dist_dist[2][dist-1]++;
						}
					}else{
						if(loc_exp>0.0 && loc_sim>0.0){
							//std::cout<<"Non-zero result for both. Filling scales\n";
							look_further_all=false;
							//std::cout<<"\t\tA Test: bin|" <<_scale_exp_7d->GetBinContent(bin) <<" vs. what should be there|" <<fun::Sparse_Integral(scale_exp) <<" after dist:" <<dist <<"\n";
							//std::cout<<"Setting exp bin error\n";
							dist_dist[0][dist-1]++;
							_scale_exp_7d->SetBinError2(_scale_exp_7d->GetBin(bin),fun::Sparse_Integral_Error2(scale_exp));
							//std::cout<<"Setting sim bin error\n";
							_scale_sim_7d->SetBinError2(_scale_sim_7d->GetBin(bin),fun::Sparse_Integral_Error2(scale_sim));
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
	std::cout<<"Finished Making the Localized holes\n";
	//Assign Errors 
	for(int i=0; i<3; i++){
		std::cout<<"\n";
		for(int j=0; j<14; j++){
			std::cout<<dist_dist[i][j] <<" ";
		}
		std::cout<<"\n\n";
	}
	std::cout<<"Named File: " <<flags_.Flags::Output_File().c_str() <<"\n";
	char hname[100];
	if(flags_.Flags::Helicity()){
		//_scale_exp_7d_pos->Divide(_acceptance_7d);
		_scale_7d_pos=(THnSparseD*)_scale_exp_7d_pos->Clone();
		_scale_7d_pos->Divide(_scale_sim_7d_pos);
		_N_holes_pos=(THnSparseD*)_sim_holes_7d->Clone();
		_N_holes_pos->Multiply(_scale_7d_pos);
		sprintf(hname,"Localized_Holes_pos");
		_N_holes_pos->SetNameTitle(hname,hname);
		_N_holes_pos->Write();
		THnSparseD* _N_holes_fifty_pos = (THnSparseD*)_N_holes_pos->Clone();
		for(long i=0; i<_N_holes_fifty_pos->GetNbins(); i++){
			_N_holes_fifty_pos->SetBinError(i,_N_holes_fifty_pos->GetBinContent(i)/2.0);
		}
		sprintf(hname,"Localized_Holes_50_pos");
		_N_holes_fifty_pos->SetNameTitle(hname,hname);
		_N_holes_fifty_pos->Write();
		//_scale_exp_7d_neg->Divide(_acceptance_7d);
		_scale_7d_neg=(THnSparseD*)_scale_exp_7d_neg->Clone();
		_scale_7d_neg->Divide(_scale_sim_7d_neg);
		_N_holes_neg=(THnSparseD*)_sim_holes_7d->Clone();
		_N_holes_neg->Multiply(_scale_7d_neg);
		sprintf(hname,"Localized_Holes_neg");
		_N_holes_neg->SetNameTitle(hname,hname);
		_N_holes_neg->Write();
		THnSparseD* _N_holes_fifty_neg = (THnSparseD*)_N_holes_neg->Clone();
		for(long i=0; i<_N_holes_fifty_neg->GetNbins(); i++){
			_N_holes_fifty_neg->SetBinError(i,_N_holes_fifty_neg->GetBinContent(i)/2.0);
		}
		sprintf(hname,"Localized_Holes_50_neg");
		_N_holes_fifty_neg->SetNameTitle(hname,hname);
		_N_holes_fifty_neg->Write();
	}else{
		std::cout<<"Dividing exp by sim scales\n";
		//_scale_exp_7d->Divide(_acceptance_7d);//They're both not acceptance corrected, so adding this would be pointless
		//_scale_7d=(THnSparseD*)_scale_sim_7d->Clone();
		//_scale_7d->Divide(_scale_exp_7d);
		_scale_7d=(THnSparseD*)_scale_exp_7d->Clone();
		_scale_7d->Divide(_scale_sim_7d);
		std::cout<<"Multiplying simulated holes by scale factor\n";
		_N_holes=(THnSparseD*)_sim_holes_7d->Clone();
		//_N_holes->Divide(_scale_7d);
		_N_holes->Multiply(_scale_7d);
		sprintf(hname,"Localized_Holes");
		_N_holes->SetNameTitle(hname,hname);
		std::cout<<"Writing to RootFile\n";
		_N_holes->Write();
		THnSparseD* _N_holes_fifty = (THnSparseD*)_N_holes->Clone();
		for(long i=0; i<_N_holes_fifty->GetNbins(); i++){
			_N_holes_fifty->SetBinError(i,_N_holes_fifty->GetBinContent(i)/2.0);
		}
		sprintf(hname,"Localized_Holes_50");
		_N_holes_fifty->SetNameTitle(hname,hname);
		_N_holes_fifty->Write();
	}
	_RootOutputFile->Close();
}

//which->{exp,sim,thrown}
void	Histogram::Sparse_7to5(Flags flags_){
	std::cout<<"Sparse_7to5\n";
	std::cout<<"\tFilling Those Sparse Friends\n";
	int bin_5d[5];
	int bin_7d[7] = {1,1,1,1,1,1,1};
	int out_n = 0;
	int out_n_pos = 0; 
	int out_n_neg=0; 
	int out_exp = 0; 
	int out_sim=0;
	int out_thrown=0; 
	int out_empty=0;
	int out_accept=0; 
	std::cout<<"N bins: " <<_N->GetNbins() <<"\n";
	if(flags_.Flags::Helicity()){
		std::cout<<"N_pos bins: " <<_N_pos->GetNbins() <<"\n";
		std::cout<<"N_neg bins: " <<_N_neg->GetNbins() <<"\n";
	}
	std::cout<<"Exp bins: " <<_exp_data_7d->GetNbins() <<"\n";
	std::cout<<"Sim bins: " <<_sim_data_7d->GetNbins() <<"\n";
	std::cout<<"Empty bins: " <<_empty_7d->GetNbins() <<"\n";
	std::cout<<"Thrown bins: " <<_thrown_7d->GetNbins() <<"\n";
	std::cout<<"Acceptance bins: " <<_acceptance_7d->GetNbins() <<"\n";
	std::cout<<"Holes bins: " <<_N_holes->GetNbins() <<"\n";
	char hname[100];

	Sparse_1d_star N_1d;
	Sparse_1d_star N_1d_pos;
	Sparse_1d_star N_1d_neg;
	Sparse_1d_star exp_1d;
	Sparse_1d_star sim_1d;
	Sparse_1d_star empty_1d;
	Sparse_1d_star accept_1d;
	Sparse_1d_star thrown_1d;
	Sparse_1d_star thrown_no_rad_1d;
	double_1d scale_factor_1d;
	Sparse_1d_star sim_holes_1d;


	int num_bins=5;
	int bins[5]={2,3,4,5,6};

	for(int i=0; i<_n_bins_7d[0]; i++){
		for(int j=0; j<_n_bins_7d[1]; j++){
			_N->GetAxis(0)->SetRange(i+1,i+1);
			_N->GetAxis(1)->SetRange(j+1,j+1);
			N_1d.push_back((THnSparseD*)_N->Projection(num_bins,bins,"E"));
			sprintf(hname,"N_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			N_1d[j]->SetNameTitle(hname,hname);
			if(flags_.Flags::Helicity()){
				_N_pos->GetAxis(0)->SetRange(i+1,i+1);
				_N_pos->GetAxis(1)->SetRange(j+1,j+1);
				N_1d_pos.push_back((THnSparseD*)_N_pos->Projection(num_bins,bins,"E"));
				sprintf(hname,"N_pos_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				N_1d_pos[j]->SetNameTitle(hname,hname);
				_N_neg->GetAxis(0)->SetRange(i+1,i+1);
				_N_neg->GetAxis(1)->SetRange(j+1,j+1);
				N_1d_neg.push_back((THnSparseD*)_N_neg->Projection(num_bins,bins,"E"));
				sprintf(hname,"N_neg_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				N_1d_neg[j]->SetNameTitle(hname,hname);
			}
			_exp_data_7d->GetAxis(0)->SetRange(i+1,i+1);
			_exp_data_7d->GetAxis(1)->SetRange(j+1,j+1);
			exp_1d.push_back((THnSparseD*)_exp_data_7d->Projection(num_bins,bins,"E"));
			sprintf(hname,"exp_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			exp_1d[j]->SetNameTitle(hname,hname);
			_sim_data_7d->GetAxis(0)->SetRange(i+1,i+1);
			_sim_data_7d->GetAxis(1)->SetRange(j+1,j+1);
			sim_1d.push_back((THnSparseD*)_sim_data_7d->Projection(num_bins,bins,"E"));
			sprintf(hname,"sim_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			sim_1d[j]->SetNameTitle(hname,hname);
			_empty_7d->GetAxis(0)->SetRange(i+1,i+1);
			_empty_7d->GetAxis(1)->SetRange(j+1,j+1);
			empty_1d.push_back((THnSparseD*)_empty_7d->Projection(num_bins,bins,"E"));
			sprintf(hname,"empty_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			empty_1d[j]->SetNameTitle(hname,hname);
			_thrown_7d->GetAxis(0)->SetRange(i+1,i+1);
			_thrown_7d->GetAxis(1)->SetRange(j+1,j+1);
			thrown_1d.push_back((THnSparseD*)_thrown_7d->Projection(num_bins,bins,"E"));
			sprintf(hname,"thrown_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			thrown_1d[j]->SetNameTitle(hname,hname);
			_thrown_7d_no_rad->GetAxis(0)->SetRange(i+1,i+1);
			_thrown_7d_no_rad->GetAxis(1)->SetRange(j+1,j+1);
			thrown_no_rad_1d.push_back((THnSparseD*)_thrown_7d_no_rad->Projection(num_bins,bins,"E"));
			sprintf(hname,"thrown_no_rad_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			thrown_no_rad_1d[j]->SetNameTitle(hname,hname);
			_acceptance_7d->GetAxis(0)->SetRange(i+1,i+1);
			_acceptance_7d->GetAxis(1)->SetRange(j+1,j+1);
			accept_1d.push_back((THnSparseD*)_acceptance_7d->Projection(num_bins,bins,"E"));
			sprintf(hname,"acceptance_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			accept_1d[j]->SetNameTitle(hname,hname);
			if(!_localized_hole_filling){
				//std::cout<<"\tPushing Back Scale Factor for W:" <<i <<" Q2:" <<j <<"\n";
				scale_factor_1d.push_back(fun::Sparse_Integral(exp_1d[j])/fun::Sparse_Integral(sim_1d[j]));
				//std::cout<<"\t\tPushing Back Sim Holes " <<"\n";
				_sim_holes_7d->GetAxis(0)->SetRange(i+1,i+1);
				_sim_holes_7d->GetAxis(1)->SetRange(j+1,j+1);
				sim_holes_1d.push_back((THnSparseD*)_sim_holes_7d->Projection(num_bins,bins,"E"));
				sprintf(hname,"sim_holes_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				sim_holes_1d[j]->SetNameTitle(hname,hname);
				//std::cout<<"\t\tAdding Scaled Holes to N" <<"\n";
				N_1d[j]->Add(sim_holes_1d[j],scale_factor_1d[j]);
			}
		}
		//std::cout<<"\tPushing back all of W:" <<i <<"\n";
		//std::cout<<"\t\tPushing back N\n";
		_N_5d.push_back(N_1d);
		N_1d.clear();
		if(flags_.Flags::Helicity()){
			//std::cout<<"\t\tPushing back N_pos\n";
			_N_5d_pos.push_back(N_1d_pos);
			N_1d_pos.clear();
			//std::cout<<"\t\tPushing back N_neg\n";
			_N_5d_neg.push_back(N_1d_neg);
			N_1d_neg.clear();
		}
		//std::cout<<"\t\tPushing back thrown\n";
		_thrown_5d.push_back(thrown_1d);
		thrown_1d.clear();
		//std::cout<<"\t\tPushing back thrown no rad\n";
		_thrown_5d_no_rad.push_back(thrown_no_rad_1d);
		thrown_no_rad_1d.clear();
		//std::cout<<"\t\tPushing back acceptance\n";
		_acceptance_5d.push_back(accept_1d);
		accept_1d.clear();
		//std::cout<<"\t\tPushing back empty\n";
		_empty_5d.push_back(empty_1d);
		empty_1d.clear();
		//std::cout<<"\t\tPushing back exp\n";
		_exp_data_5d.push_back(exp_1d);
		exp_1d.clear();
		//std::cout<<"\t\tPushing back sim\n";
		_sim_data_5d.push_back(sim_1d);
		sim_1d.clear();
		if(!_localized_hole_filling){
			//std::cout<<"\t\tPushing back scale factor\n";
			_scale_factor_5d.push_back(scale_factor_1d);
			scale_factor_1d.clear();
			//std::cout<<"\t\tPushing back sim holes\n";
			_sim_holes_5d.push_back(sim_holes_1d);
			sim_holes_1d.clear();
		}
	}
}
//For Single Differential bins
void Histogram::Single_Differential(Flags flags_){
	std::cout<<"Single Differential\n";
	if(!flags_.Flags::Plot_Single_Diff()){
		std::cout<<"\tNot Plotting Single Differential Cross Sections\n";
	 	return;
	}
	//_RootOutputFile->cd();
	//TCanvas* def = new TCanvas("def1");
	//Convert the _exp_corr_5d and _exp_holes_5d to 3 dimensional sparse histograms for usage in plotting single differential cross sections
	char hname[100];
	char xlabel[100];
	char ylabel[100];
	int idx_2d[2]; 
	std::cout<<"Making Output File\n";
	_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	std::cout<<flags_.Flags::Output_File().c_str() <<"\n";
	//_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	_RootOutputFile->cd();
	//Get the 7dimensional bins ready
	TH1D_1d_star exp_ch_1d;
	TH1D_2d_star exp_ch_2d;
	std::cout<<"\tMaking Directory\n";
	TDirectory* dir_S = _RootOutputFile->mkdir("Single Differential CS");
	//TDirectory* dir_S = _RootOutputFile.mkdir("Single Differential CS");
	dir_S->cd();
	TDirectory* dir_S1[5];
	TDirectory* dir_S2[5][_n_bins_7d[1]];
	char dirname[100];
	for(int k=0; k<5; k++){
		sprintf(dirname,"%s",_five_dim_[k]);
		dir_S1[k] = dir_S->mkdir(dirname);
		dir_S1[k]->cd();
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			sprintf(dirname,"%s_Q2|%.2f-%.2f",_five_dim_[k],_bin_low_7d[1][j],_bin_up_7d[1][j]);
			dir_S2[k][j] = dir_S1[k]->mkdir(dirname);
		}
	}
	double denom = 1.0;
	std::cout<<"\tDirectories Made\n\tWriting Histograms\n";
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			for(int k=0; k<5; k++){
				dir_S2[k][j]->cd();
				//std::cout<<"\t\tWriting Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<k <<"\n";
				sprintf(hname,"%s_single_diff_W:%.3f-%.3f_Q2:%.2f-%.2f_top:%s_var:%s",_five_dim_[k],_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				sprintf(xlabel,"%s %s",_five_dim_[k],_dim_units_[k]);
				sprintf(ylabel,"Diff CS (microbarns/(%s))",_dim_units_y_[k]);
				//std::cout<<"\t\t\tNamed stuff. Let's Project!\n";
				exp_ch_1d.push_back(_N_5d[i][j]->Projection(k,"E"));
				//std::cout<<"\t\t\tProjected. Let's Flux:\n";
				if(!flags_.Flags::Flux_Included()){//Scaling by the virtual photon flux at the center of the W Q2 bin if not already scaled
					//denom = physics::Virtual_Photon_Flux(_bin_mid_7d[0][i],_bin_mid_7d[1][j],_beam_energy_[0]);
					denom *= physics::Virtual_Photon_Flux(_bin_mid_7d[0][i],_bin_mid_7d[1][j],_beam_energy_[0]);
					//exp_ch_1d[k]->Scale(1.0/physics::Virtual_Photon_Flux(_bin_mid_7d[0][i],_bin_mid_7d[1][j],_beam_energy_[0]));
					//std::cout<<"scale denom post photon flux: " <<denom <<"\n";
				}
				denom *=flags_.Flags::L(0);
				//exp_ch_1d[k]->Scale(1.0/flags_.L(0)); //Only implemented for e16 so far
				//std::cout<<"scale denom post Luminosity: " <<denom <<"\n";
				
				//std::cout<<"\t\t\tFluxed stuff. Let's Rad Corr!\n";
				//std::cout<<"The Rad corr for this one is: " <<_rad_corr_array[i][j] <<"\n";
				if(flags_.Flags::Rad_Corr()){
					denom *= _rad_corr_array[i][j];
					//exp_ch_1d[k]->Scale(1.0/_rad_corr_array[i][j]);
					//std::cout<<"scale denom post rad corr: " <<denom <<"\n";
				}
				//std::cout<<"\t\t\tRad Corred. Let's Bin Normalize:\n";	
				denom *=_bin_size_7d[0][i];
				denom *=_bin_size_7d[1][j];
				if(k>1){
					denom *=_bin_size_5d[k][0]*TMath::Pi()/180.0;//Divide by angle bins in radians
				}else{
					denom *=_bin_size_5d[k][0];//Phi in radians
				}
				
				//std::cout<<"scale denom post bin divide: " <<denom <<"\n";
				exp_ch_1d[k]->Scale(1.0/denom);
				//exp_ch_1d[k]->Scale(1.0/(_bin_size_7d[0][i]*_bin_size_7d[1][j]*_bin_size_5d[k][0]));//Scaling by the size of the W Q^2 bins
				//std::cout<<"\t\t\tBin Normalized!\n";
				//exp_ch_1d[k]->Scale(1.0/denom);
				
				//for(int l=0; l<_n_bins_5d[k]; l++){//Scaling by the size of the given variable bin size
				//	double prev_val = exp_ch_1d[k]->GetBinContent(l+1);
				//	exp_ch_1d[k]->SetBinContent(l+1,prev_val/_bin_size_5d[k][l]);
				//}
				exp_ch_1d[k]->SetNameTitle(hname,hname);
				exp_ch_1d[k]->GetXaxis()->SetTitle(xlabel);
				exp_ch_1d[k]->GetYaxis()->SetTitle(ylabel);
				std::cout<<"Writing Histogram for W:" <<i <<" Q2:" <<j <<" Xij:" <<k <<"\n";
				exp_ch_1d[k]->Write();
				denom = 1.0;
			}
			//exp_ch_2d.push_back(exp_ch_1d);
			exp_ch_1d.clear();
		}
		//_single_diff_hist.push_back(exp_ch_2d);
		//exp_ch_2d.clear();
	}
	//_RootOutputFile->Close();
	std::cout<<"\nComplted Single Differential Cross Sections\n";
}

//For Polarization Observables
void Histogram::Polarization_Observables(Flags flags_){
	std::cout<<"Polarization Observables\n";
	if(!flags_.Flags::Plot_Pol()){
		std::cout<<"\tNot Plotting Polarization Cross Sections\n";
	 	return;
	}
	//Histogram::Calc_Error_R();
	std::cout<<"looking into the rootoutput file\n";
	_RootOutputFile->cd();
	char hname[100];
	char xlabel[100];
	char ylabel[100];
	//_RootOutputFile = new TFile(flags_.Flags::Output_File().c_str(),"RECREATE");
	TH1D_1d_star exp_ch_1d;
	TH1D_2d_star exp_ch_2d;
	TH1D_3d_star exp_ch_3d;
	TH1D* pol_obs;
	std::cout<<"\tMaking Directory\n";
	TDirectory* dir_P = _RootOutputFile->mkdir("Polarization CS");
	//TDirectory* dir_P = _RootOutputFile.mkdir("Polarization CS");
	dir_P->cd();
	TDirectory* dir_P1[4];
	TDirectory* dir_P2[5][_n_bins_7d[1]];
	TDirectory* dir_P3[5][_n_bins_7d[1]][_n_bins_5d[0]];
	TDirectory* dir_P4[5][_n_bins_7d[1]][_n_bins_5d[0]];
	char dirname[100];
	//char hname[100];
	for(int k=0; k<4; k++){
		if(flags_.Flags::Plot_Polarization(k)){
			sprintf(dirname,"%s",_five_dim_[k]);
			//std::cout<<"\tMaking TDirectory: " <<dirname <<"\n";
			dir_P1[k] = dir_P->mkdir(dirname);
			dir_P1[k]->cd();
			for(int j=0; j<_n_bins_7d[1]; j++){
				sprintf(dirname,"%s_Q2|%.2f-%.2f",_five_dim_[k],_bin_low_7d[1][j],_bin_up_7d[1][j]);
				dir_P2[k][j] = dir_P1[k]->mkdir(dirname);
				dir_P2[k][j]->cd();
				for(int m=0; m<_n_bins_5d[k]; m++){
					if(flags_.Flags::Error()){
						sprintf(dirname,"Error_%s_Q2|%.3f-%.3f_%s|%.3f-%.3f",_five_dim_[k],_bin_low_7d[1][j],_bin_up_7d[1][j],_five_dim_[k],_bin_low_5d[k][m],_bin_up_5d[k][m]);
						//std::cout<<"\t\tMaking TDirectory: " <<dirname <<"\n";
						dir_P4[k][j][m] = dir_P2[k][j]->mkdir(dirname);
						dir_P4[k][j][m]->cd();
					}else{
						sprintf(dirname,"%s_Q2|%.3f-%.3f_%s|%.3f-%.3f",_five_dim_[k],_bin_low_7d[1][j],_bin_up_7d[1][j],_five_dim_[k],_bin_low_5d[k][m],_bin_up_5d[k][m]);
						dir_P3[k][j][m] = dir_P2[k][j]->mkdir(dirname);
						dir_P3[k][j][m]->cd();
					}
				}
			}
		}
	}
	double denom=1.0;
	std::cout<<"\tWriting the Histograms\n";
	//Convert the _exp_corr_5d and _exp_holes_5d to 4 dimensional sparse histograms for usage in plotting single differential cross sections
	TH1D_1d_star exp_1d;
	TH1D_1d_star empty_1d;
	TH1D_1d_star thrown_1d;
	TH1D_1d_star sim_1d;
	TH1D_1d_star accept_1d;
	TH1D_1d_star err_hist_1d;
	double err1, err2, err3, err4, err5, err6, err7;
	double sf;
	double bin_error;

	double accept, exp, empty, sim, thrown;
	TH1D_1d_star err_hist_1d2;
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			for(int k=0; k<4; k++){//Xij
				for(int l=0; l<_n_bins_5d[k]; l++){//Bins of Xij
					if(flags_.Flags::Plot_Polarization(k)){
						if(!flags_.Flags::Error()){
							std::cout<<"\t\tAt W:" <<i <<" \tQ2:" <<j <<" \tXij:" <<k <<" \tbin:" <<l <<" \r";
							sprintf(xlabel,"Phi (deg)");
							sprintf(ylabel,"Diff CS (microbarns/(%s %s))",_dim_units_y_[k],_dim_units_y_[4]);
							_N_5d[i][j]->GetAxis(k)->SetRange(l+1,l+1);
							exp_ch_1d.push_back(_N_5d[i][j]->Projection(4,"E"));
							sprintf(hname,"CS_2d_%s:%.3f-%.3f_W:%.3f-%.3f_Q2:%.2f-%.2f_top:%s_var:%s",_five_dim_[k],_bin_low_5d[k][l],_bin_up_5d[k][l],_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
							denom *= physics::Virtual_Photon_Flux(_bin_mid_7d[0][i],_bin_mid_7d[1][j],_beam_energy_[0]);
							denom *=flags_.L(0);//Right now just e16
							if(flags_.Flags::Rad_Corr()){
								denom *=_rad_corr_array[i][j];
							}
							denom *= _bin_size_7d[0][i];
							denom *= _bin_size_7d[1][j];
							if(k>1){
								denom *=_bin_size_5d[k][l]*TMath::Pi()/180.0;
							}else{
								denom *= _bin_size_5d[k][l];
							}
							denom *=_bin_size_5d[4][0]*TMath::Pi()/180.0;
							dir_P3[k][j][l]->cd();
							exp_ch_1d[l]->Scale(1.0/denom);
							exp_ch_1d[l]->SetNameTitle(hname,hname);
							exp_ch_1d[l]->GetXaxis()->SetTitle(xlabel);
							exp_ch_1d[l]->GetYaxis()->SetTitle(ylabel);
							exp_ch_1d[l]->Write();
						}else{
							denom = 1.0;
							//std::cout<<"\t\tAt W:" <<i <<" \tQ2:" <<j <<" \tXij:" <<k <<" \tbin:" <<l <<" \t\n";
							sprintf(xlabel,"Phi (deg)");
							sprintf(ylabel,"Diff CS (microbarns/(%s %s))",_dim_units_y_[k],_dim_units_y_[4]);
							_N_5d[i][j]->GetAxis(k)->SetRange(l+1,l+1);
							exp_ch_1d.push_back(_N_5d[i][j]->Projection(4,"E"));
							sprintf(hname,"CS_2d_%s:%.3f-%.3f_W:%.3f-%.3f_Q2:%.2f-%.2f_top:%s_var:%s",_five_dim_[k],_bin_low_5d[k][l],_bin_up_5d[k][l],_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
							_exp_data_5d[i][j]->GetAxis(k)->SetRange(l+1,l+1);
							exp_1d.push_back(_exp_data_5d[i][j]->Projection(4,"E"));
							_empty_5d[i][j]->GetAxis(k)->SetRange(l+1,l+1);
							empty_1d.push_back(_empty_5d[i][j]->Projection(4,"E"));
							_thrown_5d[i][j]->GetAxis(k)->SetRange(l+1,l+1);
							thrown_1d.push_back(_thrown_5d[i][j]->Projection(4,"E"));
							_sim_data_5d[i][j]->GetAxis(k)->SetRange(l+1,l+1);
							sim_1d.push_back(_sim_data_5d[i][j]->Projection(4,"E"));
							_acceptance_5d[i][j]->GetAxis(k)->SetRange(l+1,l+1);
							accept_1d.push_back(_acceptance_5d[i][j]->Projection(4,"E"));
							sprintf(hname,"Error_CS_2d_%s:%.3f-%.3f_W:%.3f-%.3f_Q2:%.2f-%.2f_top:%s_var:%s",_five_dim_[k],_bin_low_5d[k][l],_bin_up_5d[k][l],_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
							//TH1D* pol_err = new TH1D(hname,hname,_n_bins_5d[4],_bin_low_5d[4][0],_bin_up_5d[4][_n_bins_5d[4]-1]);
							err_hist_1d.push_back(new TH1D(hname,hname,_n_bins_5d[4],_bin_low_5d[4][0],_bin_up_5d[4][_n_bins_5d[4]-1]));
							for(int m=0; m<_n_bins_5d[4]; m++){
								if(exp_ch_1d[l]->GetBinContent(m+1)>0.0){
									if(_localized_hole_filling){
										std::cout<<"localized hole filling stuff\n";
									}else{
										//sf=fun::nSparseIntegral(_exp_data_7d)/fun::nSparseIntegral(_sim_data_7d);
										sf = _scale_factor_5d[i][j];
									}
									accept = accept_1d[l]->GetBinContent(m+1);
									thrown = thrown_1d[l]->GetBinContent(m+1);
									sim = sim_1d[l]->GetBinContent(m+1);
									exp = exp_1d[l]->GetBinContent(m+1);
									empty = empty_1d[l]->GetBinContent(m+1);
									if(accept>0.0){
										accept = accept_1d[l]->GetBinContent(m+1);
										err1 = 1.0/(accept);//Nx
										if(flags_.Flags::Qr()>0.0 && !isnan(flags_.Flags::Qr())){
											err2 = flags_.Flags::Qr()/(accept);//Ne
										}else{
											err2 = 0.0;
										}
										err4 = thrown - (sim/accept);//sf
										err6 = sf/accept;//Nr
									}else{
										err1 = 0.0;
										err2 = 0.0;
										err4 = 0.0;
										err6 = 0.0;
									}
									err3 = 1.0/_rad_corr_array[i][j];//R
									err5 = sf;//Nt
									if(sim && exp>0.0 && empty>0.0 && flags_.Flags::Qr()>0.0 && !isnan(flags_.Flags::Qr())){
										err7 = (sim-exp+flags_.Flags::Qr()*empty);//A
									}else{
										if(sim>0.0){
											if(exp>0.0){
												err7 = (sim-exp);//A
											}else if(empty>0.0){
												err7 = (sim+flags_.Flags::Qr()*empty);//A
											}else{
												err7 = sim;//A
											}
										}else if(exp>0.0){
											if(empty>0.0){
												err7 = (-exp+flags_.Flags::Qr()*empty);//A
											}else{
												err7 = (-exp);//A
											}
										}else if(empty>0.0){
											err7 = (flags_.Flags::Qr()*empty);//A
										}else{
											err7 = 0.0;
										}
										
									}
									if(exp>0.0){
										err1 *= exp;//Nx
									}else{
										err1 *= 0.0;//Nx
									}
									if(empty>0.0 ){
										err2 *= empty_1d[l]->GetBinError(m+1);//Ne
									}else{
										err2 *= 0.0;//Ne
									}
									err3 *= _rad_error[i][j];//R
									//sf error is if it's calculated from the 7 dimensional space, which is not correct 5/8/23
									err4 *= TMath::Sqrt((fun::Sparse_Integral_Error(_exp_data_5d[i][j])*fun::Sparse_Integral_Error(_exp_data_5d[i][j])/(fun::Sparse_Integral(_sim_data_5d[i][j])*fun::Sparse_Integral(_sim_data_5d[i][j]))+(fun::Sparse_Integral(_exp_data_5d[i][j])*fun::Sparse_Integral(_exp_data_5d[i][j])*fun::Sparse_Integral_Error(_sim_data_5d[i][j])*fun::Sparse_Integral_Error(_sim_data_5d[i][j])/(fun::Sparse_Integral(_sim_data_5d[i][j])*fun::Sparse_Integral(_sim_data_5d[i][j])*fun::Sparse_Integral(_sim_data_5d[i][j])*fun::Sparse_Integral(_sim_data_5d[i][j])))));//sf
									//err4 *= Histogram::Error_SF();//Will need a fun function for uncertainty in sf
									if(sim>0.0 ){
										err6 *= sim_1d[l]->GetBinError(m+1);//Nr
									}else{
										err6 *= 0.0;//Nr
									}
									if(thrown>0.0){
										err5 *= thrown_1d[l]->GetBinError(m+1);//Nt
										err7 *= TMath::Sqrt((sim_1d[l]->GetBinError(m+1)*sim_1d[l]->GetBinError(m+1)/(thrown*thrown)+(sim*sim*thrown_1d[l]->GetBinError(m+1)*thrown_1d[l]->GetBinError(m+1)/(thrown*thrown*thrown*thrown))));//A
									}else{
										err5 *= 0.0;//Nt
										err7 *= 0.0;//A
									}
									bin_error = TMath::Sqrt((err1*err1)+(err2*err2)+(err3*err3)+(err4*err4)+(err5*err5)+(err6*err6)+(err7*err7));
									//err7 *= accept_1d[l]->GetBinError(m+1);//Will likely need a function for this
									/*std::cout<<"\t\tAgain we are at bin:" <<m <<" with ";
									std::cout<<"err1:" <<err1 <<" ";
									std::cout<<"err2:" <<err2 <<" ";
									std::cout<<"err3:" <<err3 <<" ";
									std::cout<<"err4:" <<err4 <<" ";
									std::cout<<"err5:" <<err5 <<" ";
									std::cout<<"err6:" <<err6 <<" ";
									std::cout<<"err7:" <<err7 <<"\n";*/
									//pol_err->SetBinContent(m+1,TMath::Sqrt((err1*err1)+(err2*err2)+(err3*err3)+(err4*err4)+(err5*err5)+(err6*err6)+(err7*err7)));
									err_hist_1d[l]->SetBinContent(m+1,bin_error);
									//std::cout<<"\t\tFull Error: " <<TMath::Sqrt((err1*err1)+(err2*err2)+(err3*err3)+(err4*err4)+(err5*err5)+(err6*err6)+(err7*err7)) <<" check bin:" <<pol_err->GetBinContent(m+1) <<"\n";
									//std::cout<<"\tExp/Emp/Sim/Thr/Acc: " <<(exp_1d[l]->GetBinContent(m+1)>0.0) <<(empty_1d[l]->GetBinContent(m+1)>0.0) <<(sim_1d[l]->GetBinContent(m+1)>0.0) <<(thrown_1d[l]->GetBinContent(m+1)>0.0) <<(accept_1d[l]->GetBinContent(m+1)>0.0) <<" | bin content:" <<exp_ch_1d[l]->GetBinContent(m+1) <<"\n";
									//std::cout<<"\t\tErr1: " <<err1*err1/(bin_error*bin_error) <<"| Err2: " <<err2*err2/(bin_error*bin_error) <<"| Err3: " <<err3*err3/(bin_error*bin_error)  <<"| Err4: " <<err4*err4/(bin_error*bin_error)  <<"| Err5: " <<err5*err5/(bin_error*bin_error)  <<"| Err6: " <<err6*err6/(bin_error*bin_error)  <<"| Err7: " <<err7*err7/(bin_error*bin_error) <<"\n";
									//std::cout<<"\t\tCalculated: " <<err_hist_1d[l]->GetBinContent(m+1) <<" for N: " <<exp_ch_1d[l]->GetBinError(m+1) <<"\n";
									err1 = 0.0;
									err2 = 0.0;
									err3 = 0.0;
									err4 = 0.0;
									err5 = 0.0;
									err6 = 0.0;
									err7 = 0.0;
								}else{
									err_hist_1d[l]->SetBinContent(m+1,0.0);
								}
							}
							denom *= physics::Virtual_Photon_Flux(_bin_mid_7d[0][i],_bin_mid_7d[1][j],_beam_energy_[0]);
							denom *=flags_.L(0);//Right now just e16
							if(flags_.Flags::Rad_Corr()){
								denom *=_rad_corr_array[i][j];
							}
							denom *= _bin_size_7d[0][i];
							denom *= _bin_size_7d[1][j];
							if(k>1){
								denom *=_bin_size_5d[k][l]*TMath::Pi()/180.0;
							}else{
								denom *= _bin_size_5d[k][l];
							}
							denom *=_bin_size_5d[4][0]*TMath::Pi()/180.0;
							dir_P4[k][j][l]->cd();
							//pol_err->Scale(1.0/denom);
							//pol_err->SetNameTitle(hname,hname);
							//pol_err->GetXaxis()->SetTitle(xlabel);
							//pol_err->GetYaxis()->SetTitle(ylabel);
							//pol_err->Write();
							err_hist_1d[l]->Scale(1.0/denom);
							err_hist_1d[l]->SetNameTitle(hname,hname);
							err_hist_1d[l]->GetXaxis()->SetTitle(xlabel);
							err_hist_1d[l]->GetYaxis()->SetTitle(ylabel);
							err_hist_1d[l]->Write();
						}
					}
					denom=1.0;
					_N_5d[i][j]->GetAxis(k)->SetRange();//Set the range back to normal
					if(flags_.Flags::Error()){
						_exp_data_5d[i][j]->GetAxis(k)->SetRange();
						_sim_data_5d[i][j]->GetAxis(k)->SetRange();
						_empty_5d[i][j]->GetAxis(k)->SetRange();
						_thrown_5d[i][j]->GetAxis(k)->SetRange();
						_acceptance_5d[i][j]->GetAxis(k)->SetRange();
					}
				}
				//exp_ch_2d.push_back(exp_ch_1d);
				exp_ch_1d.clear();
				if(flags_.Flags::Error()){
					exp_1d.clear();
					empty_1d.clear();
					sim_1d.clear();
					thrown_1d.clear();
					accept_1d.clear();
					err_hist_1d.clear();
				}
			}
			//exp_ch_3d.push_back(exp_ch_2d);
			//exp_ch_2d.clear();
		}
		//_polarization_hist.push_back(exp_ch_3d);
		//exp_ch_3d.clear();
	}
	//_RootOutputFile3->Close();
	std::cout<<"\nFinished Polarization\n";
}

void Histogram::Beam_Spin(Flags flags_){
	//Will vary by W, Q2, and Phi 
	std::cout<<"Beam Spin Asymmetry\n";
	if(!flags_.Flags::Plot_Beam_Spin() || !flags_.Flags::Helicity()){
		std::cout<<"\tNot Plotting Beam Spin Asymmetry\n";
	 	return;
	}
	_RootOutputFile->cd();
	TCanvas* def = new TCanvas("def3");
	char hname[100];
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
	//TDirectory* dir_B = _RootOutputFile.mkdir("Beam Spin");
	dir_B->cd();
	TDirectory* dir_B1[_n_bins_7d[1]];
	char dirname[100];	
	for(int j=0; j< _n_bins_7d[1]; j++){//Q2
		sprintf(dirname,"Beam_Spin_Q2|%.2f-%.2f",_bin_low_7d[1][j],_bin_up_7d[1][j]);
		dir_B1[j] = dir_B->mkdir(dirname);
	}
	int idx[3];
	std::cout<<"\tWriting Histograms\n";
	//THnSparseD* Beam_Spin_Top[_N_5d_pos.size()][_N_5d_pos[0].size()];
	//THnSparseD* Beam_Spin_Bot[_N_5d_pos.size()][_N_5d_pos[0].size()];
	//THnSparseD* Beam_Spin_R[_N_5d_pos.size()][_N_5d_pos[0].size()];
	//Convert the _exp_corr_5d and _exp_holes_5d to 4 dimensional sparse histograms for usage in plotting single differential cross sections
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			//std::cout<<"\t\tW:" <<i <<" Q2:" <<j <<"\r";
			sprintf(hname,"Beam_Spin_W|%.3f-%.3f_Q2|%.2f-%.2f",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j]);
			//beam_spin_1d.push_back(new TH1D(hname,hname,_n_bins_5d[4],_bin_low_5d[4][0],_bin_up_5d[4][_n_bins_5d[4]-1]));
			beam_spin_1d.push_back(_N_5d_pos[i][j]->Projection(4,"E"));
			//beam_spin_1d.push_back(Beam_Spin_R[i][j]->Projection(4,"E"));
			top_neg_1d.push_back(_N_5d_neg[i][j]->Projection(4,"E"));
			beam_spin_1d[j]->SetNameTitle(hname,hname);
			beam_spin_1d[j]->Add(top_neg_1d[j],-1.0);
			sprintf(hname,"Top_neg_Beam_Spin_W|%.3f-%.3f_Q2|%.2f-%.2f",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j]);
			top_neg_1d[j]->SetNameTitle(hname,hname);
			//beam_spin_1d[j]->Scale((1.0/_beam_pol_[0]));
			bot_spin_1d.push_back(_N_5d_pos[i][j]->Projection(4,"E"));
			bot_spin_1d[j]->Add(top_neg_1d[j],1.0);
			beam_spin_1d[j]->Divide(bot_spin_1d[j]);
			beam_spin_1d[j]->Scale((1.0/_beam_pol_[0]));
			//beam_spin_1d[j]->SetNameTitle(hname,hname);
			sprintf(hname,"Bot_Beam_Spin_W|%.3f-%.3f_Q2|%.2f-%.2f",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j]);
			bot_spin_1d[j]->SetNameTitle(hname,hname);
			//std::cout<<"Checking Content for " <<hname <<"\n";
			//for(int k=0; k<_n_bins_5d[4]; k++){//Phi
			//	std::cout<<"\tPhi : " <<_bin_mid_7d[6][k] <<"| N_pos:" <<N_pos->GetBinContent(i,j,k) <<" | N_neg: " <<N_neg->GetBinContent(i,j,k) <<"\n";
			//	beam_spin_1d[beam_spin_1d.size()-1]->SetBinContent(k,(1.0/_beam_pol_[0])*(N_pos->GetBinContent(i,j,k)-N_neg->GetBinContent(i,j,k))/(N_pos->GetBinContent(i,j,k)+N_neg->GetBinContent(i,j,k)));//Beam polarization needs to be specified and is not currently accurate
			//}
			dir_B1[j]->cd();
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

void Histogram::Make_WQ2(Flags flags_){
	std::cout<<"Making W Q2 histograms\n";
	if(!flags_.Flags::Plot_WQ2()){
		std::cout<<"\tNot making WQ2 Histograms\n";
		return; 
	}
	std::cout<<"Making Directory for W Q2 distributions\n";
	TCanvas* def = new TCanvas("def4");
	char hname[100];
	_RootOutputFile->cd();
	TDirectory* dir_W = _RootOutputFile->mkdir("W Q2 Distributions");
	dir_W->cd();	
	TH2D* exp_wq2 = _exp_data_7d->Projection(1,0,"E");
	sprintf(hname,"Exp_Data_W_vs_Q2");
	exp_wq2->SetNameTitle(hname,hname);
	exp_wq2->SetXTitle("W (GeV)");
	exp_wq2->SetYTitle("Q2 (GeV^2)");
	exp_wq2->Write();
	TH2D* sim_wq2 = _sim_data_7d->Projection(1,0,"E");
	sprintf(hname,"Sim_Data_W_vs_Q2");
	sim_wq2->SetNameTitle(hname,hname);
	sim_wq2->SetXTitle("W (GeV)");
	sim_wq2->SetYTitle("Q2 (GeV^2)");
	sim_wq2->Write();
	TH2D* thr_wq2 = _thrown_7d->Projection(1,0,"E");
	sprintf(hname,"Thrown_Data_W_vs_Q2");
	thr_wq2->SetNameTitle(hname,hname);
	thr_wq2->SetXTitle("W (GeV)");
	thr_wq2->SetYTitle("Q2 (GeV^2)");
	thr_wq2->Write();
	TH2D* tot_wq2 = _N->Projection(1,0,"E");
	sprintf(hname,"N_Data_W_vs_Q2");
	tot_wq2->SetNameTitle(hname,hname);
	tot_wq2->SetXTitle("W (GeV)");
	tot_wq2->SetYTitle("Q2 (GeV^2)");
	tot_wq2->Write();
	if(flags_.Flags::Helicity()){
		TH2D* exp_pos_wq2 = _exp_data_7d_pos->Projection(1,0,"E");
		sprintf(hname,"Exp_Data_Pos_W_vs_Q2");
		exp_pos_wq2->SetNameTitle(hname,hname);
		exp_pos_wq2->SetXTitle("W (GeV)");
		exp_pos_wq2->SetYTitle("Q2 (GeV^2)");
		exp_pos_wq2->Write();
		TH2D* exp_neg_wq2 = _exp_data_7d_neg->Projection(1,0,"E");
		sprintf(hname,"Exp_Data_Neg_W_vs_Q2");
		exp_neg_wq2->SetNameTitle(hname,hname);
		exp_neg_wq2->SetXTitle("W (GeV)");
		exp_neg_wq2->SetYTitle("Q2 (GeV^2)");
		exp_neg_wq2->Write();
		TH2D* empty_wq2_pos = _empty_7d_pos->Projection(1,0,"E");
		sprintf(hname,"Empty_Pos_W_vs_Q2");
		empty_wq2_pos->SetNameTitle(hname,hname);
		empty_wq2_pos->SetXTitle("W (GeV)");
		empty_wq2_pos->SetYTitle("Q2 (GeV^2)");
		empty_wq2_pos->Write();
		TH2D* empty_wq2_neg = _empty_7d_neg->Projection(1,0,"E");
		sprintf(hname,"Empty_Neg_W_vs_Q2");
		empty_wq2_neg->SetNameTitle(hname,hname);
		empty_wq2_neg->SetXTitle("W (GeV)");
		empty_wq2_neg->SetYTitle("Q2 (GeV^2)");
		empty_wq2_neg->Write();
		TH2D* tot_pos_wq2 = _N_pos->Projection(1,0,"E");
		sprintf(hname,"N_Data_Pos_W_vs_Q2");
		tot_pos_wq2->SetNameTitle(hname,hname);
		tot_pos_wq2->SetXTitle("W (GeV)");
		tot_pos_wq2->SetYTitle("Q2 (GeV^2)");
		tot_pos_wq2->Write();
		TH2D* tot_neg_wq2 = _N_neg->Projection(1,0,"E");
		sprintf(hname,"N_Data_Neg_W_vs_Q2");
		tot_neg_wq2->SetNameTitle(hname,hname);
		tot_neg_wq2->SetXTitle("W (GeV)");
		tot_neg_wq2->SetYTitle("Q2 (GeV^2)");
		tot_neg_wq2->Write();
	}
	TH2D* acc_wq2 = _acceptance_7d->Projection(1,0,"E");
	sprintf(hname,"Acceptance_W_vs_Q2");
	acc_wq2->SetNameTitle(hname,hname);
	acc_wq2->SetXTitle("W (GeV)");
	acc_wq2->SetYTitle("Q2 (GeV^2)");
	acc_wq2->Write();
	TH2D* empty_wq2 = _empty_7d->Projection(1,0,"E");
	sprintf(hname,"Empty_W_vs_Q2");
	empty_wq2->SetNameTitle(hname,hname);
	empty_wq2->SetXTitle("W (GeV)");
	empty_wq2->SetYTitle("Q2 (GeV^2)");
	empty_wq2->Write();
	sprintf(hname,"Rad_Corr_W_vs_Q2");
	_rad_corr->SetNameTitle(hname,hname);
	_rad_corr->SetXTitle("W (GeV)");
	_rad_corr->SetYTitle("Q2 (GeV^2)");
	_rad_corr->Write();
	sprintf(hname,"Sim_Holes_W_vs_Q2");
	TH2D* sim_hole_wq2= _sim_holes_7d->Projection(1,0,"E");
	sim_hole_wq2->SetNameTitle(hname,hname);
	sim_hole_wq2->SetXTitle("W (GeV)");
	sim_hole_wq2->SetYTitle("Q2 (GeV^2)");
	sim_hole_wq2->Write();
	sprintf(hname,"Holes_W_vs_Q2");
	TH2D* exp_hole_wq2= _N_holes->Projection(1,0,"E");
	exp_hole_wq2->SetNameTitle(hname,hname);
	exp_hole_wq2->SetXTitle("W (GeV)");
	exp_hole_wq2->SetYTitle("Q2 (GeV^2)");
	exp_hole_wq2->Write();
	sprintf(hname,"Holes_Pos_W_vs_Q2");
	TH2D* exp_hole_pos_wq2= _N_holes_pos->Projection(1,0,"E");
	exp_hole_pos_wq2->SetNameTitle(hname,hname);
	exp_hole_pos_wq2->SetXTitle("W (GeV)");
	exp_hole_pos_wq2->SetYTitle("Q2 (GeV^2)");
	exp_hole_pos_wq2->Write();
	sprintf(hname,"Holes_Neg_W_vs_Q2");
	TH2D* exp_hole_neg_wq2= _N_holes_neg->Projection(1,0,"E");
	exp_hole_neg_wq2->SetNameTitle(hname,hname);
	exp_hole_neg_wq2->SetXTitle("W (GeV)");
	exp_hole_neg_wq2->SetYTitle("Q2 (GeV^2)");
	exp_hole_neg_wq2->Write();
}

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
}

//Make the skeletons of the histograms to subsequently fill
void Histogram::Skeleton_5D(Flags flags_){
	std::cout<<"Skeleton_5D\n";
	char hname[100];
	int DIM = _exp_data_7d->GetNdimensions();
	int bin_sizes_5d[5];
	double bin_low_5d[5];
	double bin_up_5d[5];
	Sparse_1d_star exp_1d;
	Sparse_1d_star sim_1d;
	Sparse_1d_star empty_1d;
	Sparse_1d_star accept_1d;
	Sparse_1d_star thr_1d;
	Sparse_1d_star n_1d;
	Sparse_1d_star n_1d_pos;
	Sparse_1d_star n_1d_neg;
	Sparse_1d_star hole_exp_1d;
	Sparse_1d_star hole_exp_1d_pos;
	Sparse_1d_star hole_exp_1d_neg;
	Sparse_1d_star hole_sim_1d;
	Sparse_1d_star hole_sim_1d_pos;
	Sparse_1d_star hole_sim_1d_neg;
	for(int boop=0; boop<5; boop++){//Loop over all bins
				_n_bins_5d[boop] = _n_bins_7d[boop+2];
				bin_sizes_5d[boop] = _n_bins_5d[boop];//Get sizes for individual bins
				bin_low_5d[boop] = _bin_low_5d[boop][0];//Get low end for each bin
				bin_up_5d[boop] = _bin_up_5d[boop][_bin_up_5d[boop].size()-1];//Get high end for each bin
	}
	std::cout<<"\n# W bins: " <<_n_bins_7d[0] <<" # Wbins " <<_empty_7d->GetAxis(0)->GetNbins() <<"\n\n";
	for(int k=0; k<_n_bins_7d[0]; k++){//W	
		for(int l=0; l< _n_bins_7d[1]; l++){//Q2
			sprintf(hname,"exp_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			exp_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			sprintf(hname,"empty_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			empty_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			sprintf(hname,"sim_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			sim_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			if(_localized_hole_filling){
				sprintf(hname,"loc_hole_exp_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				hole_exp_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
				sprintf(hname,"loc_hole_exp_pos_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				hole_exp_1d_pos.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
				sprintf(hname,"loc_hole_exp_neg_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				hole_exp_1d_neg.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
				sprintf(hname,"loc_hole_sim_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				hole_sim_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
				sprintf(hname,"loc_hole_sim_pos_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				hole_sim_1d_pos.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
				sprintf(hname,"loc_hole_sim_neg_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				hole_sim_1d_neg.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			}
			sprintf(hname,"acceptance_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			accept_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			sprintf(hname,"thrown_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Var_Set().c_str());
			thr_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			sprintf(hname,"total_yield_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Var_Set().c_str());
			n_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			sprintf(hname,"pos_total_yield_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Var_Set().c_str());
			n_1d_pos.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			sprintf(hname,"neg_total_yield_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Var_Set().c_str());
			n_1d_neg.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
		}//Q2
		_exp_data_5d.push_back(exp_1d);
		_empty_5d.push_back(empty_1d);
		_sim_data_5d.push_back(sim_1d);
		_acceptance_5d.push_back(accept_1d);
		_thrown_5d.push_back(thr_1d);
		_N_5d.push_back(n_1d);
		_N_5d_pos.push_back(n_1d_pos);
		_N_5d_neg.push_back(n_1d_neg);
		if(_localized_hole_filling){
			_scale_exp_5d.push_back(hole_exp_1d);
			_scale_exp_5d_pos.push_back(hole_exp_1d_pos);
			_scale_exp_5d_neg.push_back(hole_exp_1d_neg);
			_scale_sim_5d.push_back(hole_sim_1d);
			_scale_sim_5d_pos.push_back(hole_sim_1d_pos);
			_scale_sim_5d_neg.push_back(hole_sim_1d_neg);
			hole_exp_1d.clear();
			hole_exp_1d_pos.clear();
			hole_exp_1d_neg.clear();
			hole_sim_1d.clear();
			hole_sim_1d_pos.clear();
			hole_sim_1d_neg.clear();
		}
		exp_1d.clear();
		empty_1d.clear();
		sim_1d.clear();
		accept_1d.clear();
		thr_1d.clear();
		n_1d.clear();
		n_1d_pos.clear();
		n_1d_neg.clear();
	}//W
	//for(int k=0; k<_n_bins_7d[0]; k++){//W
		//for(int l=0; l< _n_bins_7d[1]; l++){//Q2
			//_exp_data_5d[k][l]->Sumw2();
			//_sim_data_5d[k][l]->Sumw2();
			//_empty_5d[k][l]->Sumw2();
			//_thrown_5d[k][l]->Sumw2();
			//_acceptance_5d[k][l]->Sumw2();
			//_N_5d[k][l]->Sumw2();
		//}
	//}
	std::cout<<"End Skeleton_5D\n";
}

float Histogram::Bin_Size( int variable, int bin_7d){
	return _bin_up_7d[variable][bin_7d] - _bin_low_7d[variable][bin_7d];
}

float Histogram::W_Bin_Size( int bin_7d){
	return Histogram::Bin_Size(0,bin_7d);
}

float Histogram::Q2_Bin_Size( int bin_7d){
	return Histogram::Bin_Size(1,bin_7d);
}

float Histogram::X_Bin_Size( int variable, int bin_7d){
	return Histogram::Bin_Size(variable+2,bin_7d);
}

float Histogram::Phi_Bin_Size( int bin_7d){
	return Histogram::Bin_Size(6,bin_7d);
}

double Histogram::Bin_Center_7d( int dim, int bin){
	return (_bin_up_7d[dim][bin] + _bin_low_7d[dim][bin])/2.0;
}

void Histogram::XandPhi_BinHistograms(Flags flags_){
	std::cout<<"\tX and Phi Bin Histograms\n";
	char hname[100];
	//Make the histograms
	
	
	//std::cout<<"\n\tLowest phi: " <<_bin_low_5d[4][0] <<" | highest phi: " <<_bin_up_5d[4][_n_bins_5d[4]-1];
	_phi_bin_sizes = new TH1D(hname,hname,_n_bins_5d[4],_bin_low_5d[4][0],_bin_up_5d[4][_n_bins_5d[4]-1]);
	for(int j=0; j<_n_bins_5d[4]; j++){
		_phi_bin_sizes->Fill(j+1,Histogram::Phi_Bin_Size(j));
	}
	TDirectory* dir_phi_bins = _RootOutputFile->mkdir("Phi_Bins");
	dir_phi_bins->cd();
	_phi_bin_sizes->SetXTitle("Bin");
	_phi_bin_sizes->SetXTitle("Bin Size");
	_phi_bin_sizes->Write();
	TDirectory* dir_x_bins = _RootOutputFile->mkdir("X_Bins");
	dir_x_bins->cd();
	
	for(int x=0; x<4; x++){
		//std::cout<<"X: " <<x <<"\n\tLowest X: " <<_bin_low_5d[x][0] <<" | highest x: " <<_bin_up_5d[x][_n_bins_5d[x]-1] <<"\n";
		sprintf(hname,"%s_bin_sizes",_five_dim_[x]);
		_X_bin_sizes.push_back(new TH1D(hname,hname,_n_bins_5d[x],_bin_low_5d[x][0],_bin_up_5d[x][_n_bins_5d[x]-1]));
		for(int i=0; i<_n_bins_5d[x]; i++){
			_X_bin_sizes[x]->Fill(i+1,Histogram::X_Bin_Size(x,i));
		}
		_X_bin_sizes[x]->SetNameTitle(hname,hname);
		_X_bin_sizes[x]->SetXTitle("Bin");
		_X_bin_sizes[x]->SetXTitle("Bin Size");
		_X_bin_sizes[x]->Write();
	}
}

void Histogram::Acceptance_Histograms(Flags flags_){
	std::vector<long> space_dims;
	for(int i=0; i<_n_bins_7d.size(); i++){
      space_dims.push_back(_n_bins_7d[i]);
	}
	CartesianGenerator cart(space_dims);
   	int bin[_n_bins_7d.size()];
	char hname[100];
	long n_bins_accept = _acceptance_7d->GetNbins();
	std::cout<<"Number of bins filled in acceptance_7d " <<n_bins_accept <<"\n";
	std::cout<<"Number of bins filled in thrown " <<_thrown_7d->GetNbins() <<"\n";
	std::cout<<"Number of bins filled in recon " <<_sim_data_7d->GetNbins() <<"\n";
	sprintf(hname,"Acceptance_Value_Distribution0");
	TH1D* accept_vals = new TH1D(hname,hname,400,0.0,4.0);
	sprintf(hname,"Acceptance_Value_Distribution1");
	TH1D* accept_vals1 = new TH1D(hname,hname,400,0.0,4.0);
	sprintf(hname,"Acceptance_Value_Distribution2");
	TH1D* accept_vals2 = new TH1D(hname,hname,400,0.0,4.0);
	sprintf(hname,"Acceptance_Value_Distribution3");
	TH1D* accept_vals3 = new TH1D(hname,hname,400,0.0,4.0);
	sprintf(hname,"Acceptance_Value_Distribution4");
	TH1D* accept_vals4 = new TH1D(hname,hname,400,0.0,4.0);
	sprintf(hname,"Acceptance_Value_Distribution5");
	TH1D* accept_vals5 = new TH1D(hname,hname,400,0.0,4.0);
	//TH1D* accept_vals = new TH1D(hname,hname,,0.0,1.0);
	TH1D_2d_star accept_vals6;
	TH1D_2d_star accept_vals7;
	TH1D_1d_star accept_vals61;
	TH1D_1d_star accept_vals71;
	for(int i=0; i<_n_bins_7d[0]; i++){
		for(int j=0; j<_n_bins_7d[1]; j++){
			sprintf(hname,"Acceptance_Dis_W:%.3f-%.3f_Q2:%.3f-%.3f",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j]);
			accept_vals61.push_back(new TH1D(hname,hname,400,0.0,4.0));
			sprintf(hname,"Acceptance_Dis2_W:%.3f-%.3f_Q2:%.3f-%.3f",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j]);
			accept_vals71.push_back(new TH1D(hname,hname,400,0.0,4.0));
		}
		accept_vals6.push_back(accept_vals61);
		accept_vals61.clear();
		accept_vals7.push_back(accept_vals71);
		accept_vals71.clear();
	}
	double duck = NAN;
	long num_z = 0; 
	long num_z2 = 0; 
   	while(cart.GetNextCombination()){
      	for(int j = 0; j<_n_bins_7d.size(); j++){
    		bin[j] = (int)cart[j]+1;
      	}
		num_z2++;
		if(_acceptance_7d->GetBinContent(bin)>0.0){
			num_z++;
			accept_vals->Fill(_acceptance_7d->GetBinContent(bin));
		}
	}
	std::cout<<"total poss bins: " <<num_z2 <<" and one's looped over: " <<num_z <<"\n";
	long num_bins_nonzero = 0; 
	Sparse_2d_star acceptance_5d_new;
	Sparse_1d_star acceptance_5d_new1;
	int proj_bins[5] = {2,3,4,5,6};
	int bin_sizes_5d[5];
	double bin_low_5d[5];
	double bin_up_5d[5];
	for(int boop=0; boop<5; boop++){//Loop over all bins
				bin_sizes_5d[boop] = _n_bins_5d[boop];//Get sizes for individual bins
				bin_low_5d[boop] = _bin_low_5d[boop][0];//Get low end for each bin
				bin_up_5d[boop] = _bin_up_5d[boop][_bin_up_5d[boop].size()-1];//Get high end for each bin
	}
	for(int i=0; i<_n_bins_7d[0]; i++){
		for(int j=0; j<_n_bins_7d[1]; j++){
			sprintf(hname,"acceptance_5d_W:%.3f-%.3f_Q2:%.3f-%.3f",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j]);
			acceptance_5d_new1.push_back(new THnSparseD(hname,hname,5,bin_sizes_5d,bin_low_5d,bin_up_5d));
		}
		acceptance_5d_new.push_back(acceptance_5d_new1);
		acceptance_5d_new1.clear();
	}
	std::cout<<"Made skeleton of acceptance 5d \n";
	int bin_5d[5];
	for(long i = 0; i<n_bins_accept; i++){
		if(_acceptance_7d->GetBinContent(i,bin)>0.0){
			//std::cout<<"Filling 7d bin: ";
			//for(int k=0; k<7; k++){
			//	std::cout<<bin[k] <<" ";
			//}
			//std::cout<<"with: " <<_acceptance_7d->GetBinContent(i) <<" and error2: " <<_acceptance_7d->GetBinError2(i) <<"\nCheck 5d bin: ";
			for(int j=0; j<5; j++){
				bin_5d[j] = bin[j+2];
				//std::cout<<bin_5d[j] <<" ";
			}
			//std::cout<<"\n";
			acceptance_5d_new[bin[0]-1][bin[1]-1]->SetBinContent(bin_5d,_acceptance_7d->GetBinContent(i));
			acceptance_5d_new[bin[0]-1][bin[1]-1]->SetBinError2(acceptance_5d_new[bin[0]-1][bin[1]-1]->GetBin(bin_5d),_acceptance_7d->GetBinError2(i));
			//std::cout<<"Filled 5d bin with: " <<acceptance_5d_new[bin[0]-1][bin[1]-1]->GetBinContent(bin_5d) <<" and error2: " <<acceptance_5d_new[bin[0]-1][bin[1]-1]->GetBinError2(acceptance_5d_new[bin[0]-1][bin[1]-1]->GetBin(bin_5d))<<"\n";
		}
	}
	for(long i = 0; i<n_bins_accept; i++){
		if(_acceptance_7d->GetBinContent(i)>0.0){
			num_bins_nonzero++;
			accept_vals1->Fill(_acceptance_7d->GetBinContent(i));
		}
		duck=_acceptance_7d->GetBinContent(i,bin);
		if(duck>0.0){
			accept_vals2->Fill(_acceptance_7d->GetBinContent(i));
		}
		if(_acceptance_7d->GetBinContent(bin)>0.0){
			accept_vals3->Fill(_acceptance_7d->GetBinContent(bin));
		}
		if(_acceptance_7d->GetBinContent(i,bin)>0.0){
			for(int j=0; j<5; j++){
				bin_5d[j] = bin[j+2];
			}
			accept_vals5->Fill(acceptance_5d_new[bin[0]-1][bin[1]-1]->GetBinContent(bin_5d));
			accept_vals6[bin[0]-1][bin[1]-1]->Fill(acceptance_5d_new[bin[0]-1][bin[1]-1]->GetBinContent(bin_5d));
		}
	}
	std::cout<<"Final Loop\n";
	for(int i=0; i<_n_bins_7d[0]; i++){
		for(int j=0; j<_n_bins_7d[1]; j++){
			std::cout<<"How many bins for new" <<i <<" " <<j <<" -> " <<acceptance_5d_new[i][j]->GetNbins() <<"\n";
			std::cout<<"How many bins for old" <<i <<" " <<j <<" -> " <<_acceptance_5d[i][j]->GetNbins() <<"\n";
			for(long k=0; k<_acceptance_5d[i][j]->GetNbins(); k++){
				if(_acceptance_5d[i][j]->GetBinContent(k)>0.0){
					accept_vals4->Fill(_acceptance_5d[i][j]->GetBinContent(k));
				}
			}
			for(long k=0; k<acceptance_5d_new[i][j]->GetNbins(); k++){
				if(acceptance_5d_new[i][j]->GetBinContent(k)>0.0){
					accept_vals7[i][j]->Fill(acceptance_5d_new[i][j]->GetBinContent(k));
					accept_vals->Fill(acceptance_5d_new[i][j]->GetBinContent(k));
				}
			}
		}
	}
	std::cout<<"Nonzero acceptance bins: " <<num_bins_nonzero <<"\n";
	TDirectory* dir_acc = _RootOutputFile->mkdir("Acceptance");
	dir_acc->cd();
	accept_vals->SetXTitle("Acceptance Value");
	accept_vals->SetYTitle("Num Bins");
	accept_vals->Write();
	accept_vals1->SetXTitle("Acceptance Value");
	accept_vals1->SetYTitle("Num Bins");
	accept_vals1->Write();
	accept_vals2->SetXTitle("Acceptance Value");
	accept_vals2->SetYTitle("Num Bins");
	accept_vals2->Write();
	accept_vals3->SetXTitle("Acceptance Value");
	accept_vals3->SetYTitle("Num Bins");
	accept_vals3->Write();
	accept_vals4->SetXTitle("Acceptance Value");
	accept_vals4->SetYTitle("Num Bins");
	accept_vals4->Write();
	accept_vals5->SetXTitle("Acceptance Value");
	accept_vals5->SetYTitle("Num Bins");
	accept_vals5->Write();
	for(int i=0; i<_n_bins_7d[0]; i++){
		for(int j=0; j<_n_bins_7d[1]; j++){
			accept_vals6[i][j]->SetXTitle("Acceptance Value");
			accept_vals6[i][j]->SetYTitle("Num Bins");
			accept_vals6[i][j]->Write();
			accept_vals7[i][j]->SetXTitle("Acceptance Value");
			accept_vals7[i][j]->SetYTitle("Num Bins");
			accept_vals7[i][j]->Write();
		}
	}
}

void Histogram::Hole_Histograms(Flags flags_){
	std::vector<long> space_dims;
	for(int i=0; i<_n_bins_7d.size(); i++){
      space_dims.push_back(_n_bins_7d[i]);
	}
	CartesianGenerator cart(space_dims);
   	int bin[_n_bins_7d.size()];
	char hname[100];
	sprintf(hname,"Hole Fraction");
	TH1D* hole_frac = new TH1D(hname,hname,300,0.0,1.0);
   	while(cart.GetNextCombination()){
      	for(int j = 0; j<_n_bins_7d.size(); j++){
    		bin[j] = cart[j]+1;
      	}
		if(_N->GetBinContent(bin)>0.0){
			hole_frac->Fill(_N_holes->GetBinContent(bin)/_N->GetBinContent(bin));

		}
	}
	TDirectory* dir_acc = _RootOutputFile->mkdir("Holes");
	dir_acc->cd();
	hole_frac->SetXTitle("Hole Fraction");
	hole_frac->SetXTitle("Num Bins");
	hole_frac->Write();
}

void Histogram::Print_Histogram_Bin_Info(THnSparseD* hist_){
	std::cout<<"Printing Histogram Bin Info\n\tTotal Dimensions: " <<hist_->GetNdimensions() <<"\n";
	for(int i=0; i<hist_->GetNdimensions(); i++){
		std::cout<<"\tDimension " <<i <<"\n";
		std::cout<<"\t\tNum Bins: " <<hist_->GetAxis(i)->GetNbins() <<"\n";
		for(int j=0; j<hist_->GetAxis(i)->GetNbins(); j++){
			std::cout<<"\t\t\tBin " <<j <<": " <<hist_->GetAxis(i)->GetBinLowEdge(j+1) <<" - " <<hist_->GetAxis(i)->GetBinUpEdge(j+1) <<"\n";
		}
	}
}

void Histogram::Calc_Error_R(Flags flags_){
	if(!flags_.Flags::Error()){
		return;
	}
	std::cout<<"Calculating Statistical Error for Radiative Correction\n";
	TH2D* Nt = _thrown_7d->Projection(1,0,"E");
	Nt->SetNameTitle("thrown_nt","thrown_nt");
	//double Nt_err = Nt->GetBinError(Q2_bin_+1,W_bin_+1);
	TH2D* Ntn = _thrown_7d_no_rad->Projection(1,0,"E");
	Ntn->SetNameTitle("thrown_ntn","thrown_ntn");
	//double Ntn_err = Ntn->GetBinError(Q2_bin_+1,W_bin_+1);
	//TH2D* Nt = _thrown_5d[W_bin_][Q2_bin_]->Projection(1,0);
	double int_Nt_err = fun::Sparse_Integral_Error(_thrown_7d);
	double int_Ntn_err = fun::Sparse_Integral_Error(_thrown_7d_no_rad);
	std::cout<<"int_Nt_err: " <<int_Nt_err <<"  int_Ntn_err: " <<int_Ntn_err <<"\n";
	double t_int = fun::Sparse_Integral(_thrown_7d);
	double tn_int = fun::Sparse_Integral(_thrown_7d_no_rad);
	std::cout<<"int_Nt: " <<t_int <<"  int_Ntn: " <<tn_int <<"\n";
	double mu_err = TMath::Sqrt((int_Nt_err/tn_int)*(int_Nt_err/tn_int) + (int_Ntn_err*t_int/(tn_int*tn_int))*(int_Ntn_err*t_int/(tn_int*tn_int)));
	std::cout<<"mu: " <<_rad_corr_mu <<" mu_err: " <<mu_err <<"\n"; 
	double_1d r_err_1d;
	double Nt_err, Ntn_err; 
	double err1, err2, err3;
	
	for(int i=0; i<_n_bins_7d[0]; i++){
		for(int j=0; j<_n_bins_7d[1]; j++){
			//std::cout<<"W:" <<i <<" Q2:" <<j <<"\n";
			//Nt_err = Nt->GetBinError(i+1,j+1);
			Nt_err = fun::Sparse_Integral_Error(_thrown_5d[i][j]);
			Ntn_err = fun::Sparse_Integral_Error(_thrown_5d_no_rad[i][j]);
			err1 = (_rad_corr_mu*Ntn_err/fun::Sparse_Integral(_thrown_5d[i][j]));
			err2 = (_rad_corr_mu*fun::Sparse_Integral(_thrown_5d_no_rad[i][j])*Nt_err/(fun::Sparse_Integral(_thrown_5d[i][j])*fun::Sparse_Integral(_thrown_5d[i][j])));
			err3 = fun::Sparse_Integral(_thrown_5d_no_rad[i][j])*mu_err/fun::Sparse_Integral(_thrown_5d[i][j]);
			r_err_1d.push_back(TMath::Sqrt(err1*err1+err2*err2+err3*err3));
		}
		_rad_error.push_back(r_err_1d);
		r_err_1d.clear();
	}
}

/*void Histogram::Calc_Error_A(Flags flags_){
	if(!flags_.Flags::Error()){
		return;
	}
	Sparse_1d_star err_1d;
	for(int i=0; i<_n_bins_7d[0]; i++){
		for(int j=0; j<_n_bins_7d[1]; j++){
			err_1d.push_back((THnSparse*) )
		}
	}
}

void Histogram::Calc_Error_Single_Diff(Flags flags_){
	if(!flags_.Flags::Error()){
		return;
	}

}
void Histogram::Calc_Error_Polarization(Flags flags_){
	if(!flags_.Flags::Error()){
		return;
	}

}
void Histogram::Calc_Error_Beam_Spin(Flags flags_){
	if(!flags_.Flags::Error()){
		return;
	}

}*/