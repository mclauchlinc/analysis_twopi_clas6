#include "histogram.hpp"

Histogram::Histogram(const std::string& output_file, TFile *exp_tree, TFile *sim_tree, Flags flags_){
	//Make output file eventually RootOutputFile = fun::make_file
	std::cout<<"Name_File\n";
	_RootOutputFile = fun::Name_File(output_file);
	Histogram::Extract_7d_Histograms(exp_tree,sim_tree,flags_);//Gets the 7d histograms off the root file
	Histogram::Extract_Bin_Info(flags_);//Extracts binning information about the 7d histograms just extracted
	Histogram::Skeleton_5D(flags_);//Creates empty vector arrays to map onto in the next step
	Histogram::Sparse_7to5(flags_);//Converts all the 7d histograms into usable 5d histograms for analysis
	Histogram::Sparse_5to3(flags_);
	Histogram::Sparse_5to4(flags_);
	Histogram::Convert_to_Cross(flags_);
	Histogram::Calc_Cross_Section(flags_);
	Histogram::Acceptance_Errors(flags_);
	Histogram::Make_Histograms(flags_);
	Histogram::Fill_Histograms(flags_);
	Histogram::Write_Histograms(flags_);
}

void Histogram::Make_Histograms(Flags flags_){
	_RootOutputFile=Histogram::Name_Output(flags_);
	TCanvas* def = new TCanvas("def");
	std::cout<<"Making Histograms\n";
	Histogram::Make_WQ2(flags_);
	Histogram::Make_Acceptance(flags_);
	Histogram::Make_Single_Diff(flags_);
	Histogram::Make_Polarization(flags_);
	//Histogram::Make_Acceptance_Statistics(flags_); //Cannot perform here. Require rootfile level for sim recon data
	Histogram::Make_Error_Hists(flags_);
}

void Histogram::Fill_Histograms(Flags flags_){
	Histogram::Fill_Error_Hists(flags_);
}
	
void Histogram::Write_Histograms(Flags flags_){
	std::cout<<"Writing Histograms\n";
	//Histogram::Write_WQ2(flags_); //Already written when making
	//Histogram::Write_Acceptance(flags_); //Already written when making
	//Histogram::Write_Single_Diff(flags_);
	Histogram::Write_Error_Hists(flags_);
	//Histogram::Write_Polarization(flags_);
	_RootOutputFile->Close();
	std::cout<<"Output File: " <<flags_.Flags::Output_File() <<"\n";
}

std::shared_ptr<TFile> Histogram::Name_Output(Flags flags_){
	return std::make_shared<TFile>(flags_.Flags::Output_File().c_str(),"RECREATE");
}

void Histogram::Extract_7d_Histograms(TFile *exp_tree, TFile *sim_tree, Flags flags_){
	std::cout<<"Extract 7d Histograms\n";
	char hname[100];
	//static const char * _sparse_names_[] = {"2#pi_off_proton_#Delta^{++}","2#pi_off_proton_#rho","2#pi_off_proton_#Delta^{0}"};
	//static const char * topo[] = {"Pmiss","PIPmiss","PIMmiss","Zeromiss","ALLmiss"};
	if(flags_.Flags::Flux_Included()){
		sprintf(hname,"Scaled_Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
		std::cout<<"Getting Thrown THnSparse " <<hname <<"\n";
		_thrown_7d = (THnSparseD *)sim_tree->Get(hname);
		if(flags_.Flags::Helicity()){
			sprintf(hname,"Scaled_%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);//_top_[flags_.Flags::Top_idx()]);
			std::cout<<"Getting Exp THnSparse Pos" <<hname <<"\n";
			_exp_data_7d_pos = (THnSparseD *)exp_tree->Get(hname);
			sprintf(hname,"Scaled_%s_%s_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);//_top_[flags_.Flags::Top_idx()]);
			std::cout<<"Getting Exp THnSparse Neg" <<hname <<"\n";
			_exp_data_7d_neg = (THnSparseD *)exp_tree->Get(hname);
			std::cout<<"Getting Exp THnSparse " <<hname <<"\n";
			_exp_data_7d = Histogram::Add_Sparse(_exp_data_7d_pos,_exp_data_7d_neg);
		}else{
			sprintf(hname,"Scaled_%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);//_top_[flags_.Flags::Top_idx()]);
			std::cout<<"Getting Exp THnSparse " <<hname <<"\n";
			_exp_data_7d = (THnSparseD *)exp_tree->Get(hname);
		}
		sprintf(hname,"Scaled_%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting Sim Recon THnSparse " <<hname <<"\n";
		_sim_data_7d = (THnSparseD *)sim_tree->Get(hname);
		sprintf(hname,"Weight_%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting Weight THnSparse " <<hname <<"\n";
	}else{
		sprintf(hname,"Thrown_%s",_sparse_names_[flags_.Flags::Var_idx()]);
		std::cout<<"Getting Thrown THnSparse " <<hname <<"\n";
		_thrown_7d = (THnSparseD *)sim_tree->Get(hname);
		if(flags_.Flags::Helicity()){
			sprintf(hname,"%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);//_top_[flags_.Flags::Top_idx()]);
			std::cout<<"Getting Exp THnSparse Pos" <<hname <<"\n";
			_exp_data_7d_pos = (THnSparseD *)exp_tree->Get(hname);
			sprintf(hname,"%s_%s_neg",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);//_top_[flags_.Flags::Top_idx()]);
			std::cout<<"Getting Exp THnSparse Neg" <<hname <<"\n";
			_exp_data_7d_neg = (THnSparseD *)exp_tree->Get(hname);
			std::cout<<"Getting Exp THnSparse " <<hname <<"\n";
			_exp_data_7d = Histogram::Add_Sparse(_exp_data_7d_pos,_exp_data_7d_neg);
		}else{
			sprintf(hname,"%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);//_top_[flags_.Flags::Top_idx()]);
			std::cout<<"Getting Exp THnSparse " <<hname <<"\n";
			_exp_data_7d = (THnSparseD *)exp_tree->Get(hname);
		}
		sprintf(hname,"%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting Sim Recon THnSparse " <<hname <<"\n";
		_sim_data_7d = (THnSparseD *)sim_tree->Get(hname);
		sprintf(hname,"Weight_%s_%s",_sparse_names_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
		std::cout<<"Getting Weight THnSparse " <<hname <<"\n";
	}
	_sim_weight_sq_7d = (THnSparseD *)sim_tree->Get(hname);
	_acceptance_7d = (THnSparseD*)_sim_data_7d->Clone();
	_acceptance_7d->Divide(_thrown_7d);
	if(fun::nSparseIntegral(_sim_data_7d)>0.0){
		std::cout<<"\tSim Data 7d integral is non-zero\n";
		_scale_factor_7d=fun::nSparseIntegral(_exp_data_7d)/fun::nSparseIntegral(_sim_data_7d);
		std::cout<<"\tScale factor 7d set to " <<_scale_factor_7d <<"\n";
	}else{
		std::cout<<"\tSim Data 7d integral is zero\n";
		_scale_factor_7d=0.0;
		std::cout<<"\tScale factor 7d set to " <<_scale_factor_7d <<"\n";
	}
	_exp_corr_7d=(THnSparseD*)_exp_data_7d->Clone();
	_exp_corr_7d->Divide(_acceptance_7d);
	//std::cout<<"\nexp_corr integral: " <<fun::nSparseIntegral(_exp_data_7d);
	//std::cout<<"Getting Sim Corr "<<_sparse_names_[i] <<" " <<_top_[j] <<"\n";
	_sim_corr_7d=(THnSparseD*)_sim_data_7d->Clone();
	_sim_corr_7d->Divide(_acceptance_7d);
}



void Histogram::Sparse_Add_7d(THnSparseD &h0,THnSparseD* h1, THnSparseD* h2, int sign){
	std::cout<<"Sparse Add 7d\n";
	//THnSparseD output = THnSparseD(hname,hname,h1.GetNdimensions(),fun::Vector_Array(_n_bins_7d),fun::Vector_Array(_bin_lowside_7d.at(var_set)),fun::Vector_Array(_bin_topside_7d.at(var_set)));
	int dim1 = h1->GetNdimensions();
	//int dim1 = h1->GetNdimensions();
	int dim2 = h2->GetNdimensions();
	//int dim0 = h0->GetNdimensions();
	int dim0 = h0.GetNdimensions();
	int bin_7d[7];
	if(dim1 == dim2 && dim1==dim0 && dim0==7){
		//sprintf(hname,"exp_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[l][0][i],_bin_up_7d[l][0][i],_bin_low_7d[l][1][j],_bin_up_7d[l][1][j],topo[k],var_set[l]);
		//output = THnSparseD(hname,hname,dim1,fun::Vector_Array(_n_bins_7d),fun::Vector_Array(_bin_lowside_7d.at(var_set)),fun::Vector_Array(_bin_topside_7d.at(var_set)));
		for(int idx =0; idx<dim1; idx++){
			//h0->GetAxis(idx)->Set(_n_bins_7d[idx],fun::Vector_Array(_bin_edges_7d.at(var_set).at(idx)));
			h0.GetAxis(idx)->Set(_n_bins_7d[idx],fun::Vector_Array(_bin_edges_7d[idx]));
		}
		for(int i=0; i<_n_bins_7d[0]; i++){//W
			bin_7d[0]=i;
			for(int j=0; j< _n_bins_7d[1]; j++){//Q2
				bin_7d[1]=j;
				for(int k = 0; k< _n_bins_7d[2]; k++){ //MM1
					bin_7d[2]=k;
					for(int l = 0; l< _n_bins_7d[3]; l++){ //MM2
						bin_7d[3]=l;
						for(int n = 0; n< _n_bins_7d[4]; n++){//theta
							bin_7d[4]=n;
							for(int m = 0; m< _n_bins_7d[5]; m++){ //alpha
								bin_7d[5]=m;
								for(int o = 0; o< _n_bins_7d[6]; o++){//Phi
									bin_7d[6]=o;
									//Only fill the sparse histogram where values are non-zero
									if(h1->THnSparse::GetBinContent(bin_7d) > 0 && h2->THnSparse::GetBinContent(bin_7d) > 0){//If the bin doesn't exist it will yield 0, but not create a bin
										//h0->THnSparse::AddBinContent(bin_7d,h1->THnSparse::GetBinContent(bin_7d)+sign*h2->THnSparse::GetBinContent(bin_7d));
										h0.THnSparse::AddBinContent(bin_7d,h1->THnSparse::GetBinContent(bin_7d)+sign*h2->THnSparse::GetBinContent(bin_7d));
									}else if(h1->THnSparse::GetBinContent(bin_7d) > 0 ){//If the bin doesn't exist it will yield 0, but not create a bin
										//h0->THnSparse::AddBinContent(bin_7d,h1->THnSparse::GetBinContent(bin_7d));
										h0.THnSparse::AddBinContent(bin_7d,h1->THnSparse::GetBinContent(bin_7d));
									}else if(h2->THnSparse::GetBinContent(bin_7d) > 0){//If the bin doesn't exist it will yield 0, but not create a bin
										//h0->THnSparse::AddBinContent(bin_7d,sign*h2->THnSparse::GetBinContent(bin_7d));
										h0.THnSparse::AddBinContent(bin_7d,sign*h2->THnSparse::GetBinContent(bin_7d));
									}
								}
							}
						}
					}
				}
			}
		}
	}else{
		std::cout<<"histograms do not share the same dimensions" <<std::endl;
	}
}

void Histogram::Sparse_Add_5d(THnSparseD* &h0, THnSparseD* h1, THnSparseD* h2, int sign){
	std::cout<<"Sparse Add 5d\n";
	//THnSparseD output = THnSparseD(hname,hname,h1.GetNdimensions(),fun::Vector_Array(_n_bins_5d),fun::Vector_Array(_bin_lowside_5d.at(var_set)),fun::Vector_Array(_bin_topside_5d.at(var_set)));
	std::cout<<"\tAdding 5d Hist\n";
	int dim1 = h1->GetNdimensions();
	int dim2 = h2->GetNdimensions();
	int dim0 = h0->GetNdimensions();
	//int dim0 = h0.GetNdimensions();
	int bin_5d[5];
	if(dim1 == dim2 && dim1==dim0 && dim0==5){
		//output = THnSparseD(hname,hname,dim1,fun::Vector_Array(_n_bins_5d),fun::Vector_Array(_bin_lowside_5d.at(var_set)),fun::Vector_Array(_bin_topside_5d.at(var_set)));
		//not sure why that's here
		//for(int idx =0; idx<dim1; idx++){
		//	h0->GetAxis(idx)->Set(_n_bins_5d[idx],fun::Vector_Array(_bin_edges_5d.at(var_set).at(idx)));
		//}
		for(int k = 0; k< _n_bins_5d[0]; k++){
			bin_5d[0]=k;
			for(int l = 0; l< _n_bins_5d[1]; l++){
				bin_5d[1]=l;
				for(int n = 0; n< _n_bins_5d[2]; n++){
					bin_5d[2]=n;
					for(int m = 0; m< _n_bins_5d[3]; m++){
						bin_5d[3]=m;
						for(int o = 0; o< _n_bins_5d[4]; o++){
							bin_5d[4]=o;
							//Only fill the sparse histogram where values are non-zero
							if(h1->THnSparse::GetBinContent(bin_5d) > 0 && h2->THnSparse::GetBinContent(bin_5d) > 0){//If the bin doesn't exist it will yield 0, but not create a bin
								h0->THnSparse::AddBinContent(bin_5d,h1->THnSparse::GetBinContent(bin_5d)+sign*h2->THnSparse::GetBinContent(bin_5d));
								//h0.THnSparse::AddBinContent(bin_5d,h1->THnSparse::GetBinContent(bin_5d)+sign*h2->THnSparse::GetBinContent(bin_5d));
							}else if(h1->THnSparse::GetBinContent(bin_5d) > 0 ){//If the bin doesn't exist it will yield 0, but not create a bin
								h0->THnSparse::AddBinContent(bin_5d,h1->THnSparse::GetBinContent(bin_5d));
								//h0.THnSparse::AddBinContent(bin_5d,h1->THnSparse::GetBinContent(bin_5d));
							}else if(h2->THnSparse::GetBinContent(bin_5d) > 0){//If the bin doesn't exist it will yield 0, but not create a bin
								h0->THnSparse::AddBinContent(bin_5d,sign*h2->THnSparse::GetBinContent(bin_5d));
								//h0.THnSparse::AddBinContent(bin_5d,sign*h2->THnSparse::GetBinContent(bin_5d));
							}
						}
					}
				}
			}
		}
	}else{
		std::cout<<"histograms do not share the same dimensions" <<std::endl;
	}
}

//which->{exp,sim,thrown}
void	Histogram::Sparse_7to5(Flags flags_){
	std::cout<<"Sparse_7to5\n";
	std::cout<<"\tFilling Those Sparse Friends\n";
	//Get the 7dimensional bins ready
	int bin_5d[5];
	int bin_7d[7];
	//std::cout<<"\nPart 1";
	//long n_combo = _n_bins_7d[0]*_n_bins_7d[1]*_n_bins_7d[2]*_n_bins_7d[3]*_n_bins_7d[4]*_n_bins_7d[5]*_n_bins_7d[6];
	//std::cout<<"Bins to iterate: " <<n_combo <<"\nPer:";

	long idx = 0;
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		bin_7d[0]=i;
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			bin_7d[1]=j;
			for(int k = 0; k< _n_bins_7d[2]; k++){//MM1
				bin_5d[0]=k;
				bin_7d[2]=k;
				for(int l = 0; l< _n_bins_7d[3]; l++){//MM2
					bin_5d[1]=l;
					bin_7d[3]=l;
					for(int m = 0; m< _n_bins_7d[4]; m++){//Theta
						bin_5d[2]=m;
						bin_7d[4]=m;
						for(int n = 0; n< _n_bins_7d[5]; n++){//Alpha
							bin_5d[3]=n;
							bin_7d[5]=n;
							for(int o = 0; o< _n_bins_7d[6]; o++){//Phi
								//if(((o+1)*(n+1)*(m+1)*(l+1)*(k+1)*(j+1)*(i+1)*-1)%(n_combo/100)==0){
								//	std::cout<<"\r" <<"\t" <<(100*((o+1)*(n+1)*(m+1)*(l+1)*(k+1)*(j+1)*(i+1)*-1)/n_combo) <<" %" <<std::flush;
								//}
								bin_5d[4]=o;
								bin_7d[6]=o;
								//Only fill the sparse histogram where values are non-zero
								//Exp 7d-5d
								//std::cout<<"\nMaking an Experimental Histogram: ";
								
								if(_exp_data_7d->THnSparse::GetBinContent(bin_7d) > 0.0){//If the bin doesn't exist it will yield 0, but not create a bin
									_exp_data_5d[i][j]->THnSparse::AddBinContent(bin_5d,_exp_data_7d->THnSparse::GetBinContent(bin_7d));
									//std::cout<<"Success at "<<bin_5d;
								}//Need this for scale factor_5d

								
								//Sim Recon 7d-5d
								//std::cout<<"\nMaking a Reconstructed Histogram: ";
								if(_sim_data_7d->THnSparse::GetBinContent(bin_7d) > 0.0){//If the bin doesn't exist it will yield 0, but not create a bin
									_sim_data_5d[i][j]->THnSparse::AddBinContent(bin_5d,_sim_data_7d->THnSparse::GetBinContent(bin_7d));
									//std::cout<<"Sucess at " <<bin_5d;
								}

								//Thrown 7d-5d
								//std::cout<<"\nMaking a Thrown Histogram: ";
								if(_thrown_7d->THnSparse::GetBinContent(bin_7d) > 0.0){//If the bin doesn't exist it will yield 0, but not create a bin
									_thrown_5d[i][j]->THnSparse::AddBinContent(bin_5d,_thrown_7d->THnSparse::GetBinContent(bin_7d));
									//std::cout<<"Sucess at " <<bin_5d;
								}
								//Acceptance 7d-5d
								//std::cout<<"\nMaking an Acceptance Histogram: ";
								

								if(_acceptance_7d->THnSparse::GetBinContent(bin_7d) > 0.0){
									_acceptance_5d[i][j]->THnSparse::AddBinContent(bin_5d,_acceptance_7d->THnSparse::GetBinContent(bin_7d));
									//std::cout<<"Sucess at " <<bin_5d;
								}
								//Acceptance Corrected Yields
								if(_exp_corr_7d->THnSparse::GetBinContent(bin_7d) > 0.0){
									_exp_corr_5d[i][j]->THnSparse::AddBinContent(bin_5d,_exp_corr_7d->THnSparse::GetBinContent(bin_7d));
									//_exp_corr_holes_5d[i][j]->THnSparse::AddBinContent(bin_5d,_exp_corr_7d->THnSparse::GetBinContent(bin_7d));
								}
								if(_sim_corr_7d->THnSparse::GetBinContent(bin_7d) > 0.0){
									_sim_corr_5d[i][j]->THnSparse::AddBinContent(bin_5d,_sim_corr_7d->THnSparse::GetBinContent(bin_7d));
								}
								/*idx= idx +1;
								if(idx%(n_combo/100)==0){
									std::cout<<"\r" <<"\t" <<(100*idx/n_combo) <<" %   "  <<std::flush;
								}*/
							}
						}
					}
				}
			}
			_scale_factor_5d[i][j]=fun::nSparseIntegral(_exp_data_5d[i][j])/fun::nSparseIntegral(_sim_data_5d[i][j]);
			//std::cout<<"Exp Hist 5D Yields| var:" <<var <<"top:" <<top <<"W:" <<i <<" Q2:" <<j <<" => " <<fun::nSparseIntegral(_exp_data_5d[top][i][j]) <<"\n";
			//std::cout<<"Sim Hist 5D Yields| var:" <<var <<"top:" <<top <<"W:" <<i <<" Q2:" <<j <<" => " <<fun::nSparseIntegral(_sim_data_5d[top][i][j]) <<"\n";
			//if(top==0){
				//std::cout<<"Thr Hist 5D Yields| var:" <<var <<"W:" <<i <<" Q2:" <<j <<" => " <<fun::nSparseIntegral(_thrown_5d[i][j]) <<"\n";
			//}
			//if(fun::nSparseIntegral(_acceptance[i][j])>0.0 && fun::nSparseIntegral(_thrown_5d[i][j])>0.0){
			//	_acceptance[i][j]->Divide(_thrown_5d[i][j]);
			//}
			if(fun::nSparseIntegral(_acceptance_5d[i][j]) > 0.0 && fun::nSparseIntegral(_sim_corr_5d[i][j]) >0.0){
				//std::cout<<"\nhave Acceptance and Sim data. ";
				//std::cout<<"\n\t|Adding Sim data to corr hist |";
				//_sim_corr_5d[i][j]= (THnSparseD*)_sim_data_5d[i][j]->Clone();
				//std::cout<<"\n\tDividing by acceptance |";
				//_sim_corr_5d[i][j]->Divide(_acceptance[i][j]);
				//std::cout<<"\n\tMade Sim Acceptance Corrected |";
				for(int k = 0; k< _n_bins_7d[2]; k++){
					bin_5d[0]=k;
					for(int l = 0; l< _n_bins_7d[3]; l++){
						bin_5d[1]=l;
						for(int n = 0; n< _n_bins_7d[4]; n++){
							bin_5d[2]=n;
							for(int m = 0; m< _n_bins_7d[5]; m++){
								bin_5d[3]=m;
								for(int o = 0; o< _n_bins_7d[6]; o++){
									bin_5d[4]=o;
					//std::cout<<"\n\tMaking Sim Holes: Adding Thrown |";
									
									if(_thrown_5d[i][j]->THnSparse::GetBinContent(bin_5d)>0.0){
										/*std::cout<<"\n5d bin: {";
										for(int boor=0; boor<5; boor++){
											std::cout<<bin_5d[boor];
											if(boor!=4){
												std::cout<<", ";
											}else{
												std::cout<<"}";
											}
										}
										std::cout<<"\n\t thrown bin: " <<_thrown_5d[i][j]->THnSparse::GetBinContent(bin_5d);*/
										_sim_holes_5d[i][j]->THnSparse::AddBinContent(bin_5d,_thrown_5d[i][j]->THnSparse::GetBinContent(bin_5d));
									}
					//std::cout<<"\n\tSubtracting Sim Corr |";
									if(_sim_corr_5d[i][j]->THnSparse::GetBinContent(bin_5d) >0.0){
										/*std::cout<<"\n5d bin: {";
										for(int boor=0; boor<5; boor++){
											std::cout<<bin_5d[boor];
											if(boor!=4){
												std::cout<<", ";
											}else{
												std::cout<<"}";
											}
										}
										std::cout<<"\n\t -sim corr bin: " <<-1*(_sim_corr_5d[i][j]->GetBinContent(bin_5d));*/
										_sim_holes_5d[i][j]->THnSparse::AddBinContent(bin_5d,-1*(_sim_corr_5d[i][j]->THnSparse::GetBinContent(bin_5d)));
									}
								}
							}
						}
					}
				}
				//std::cout<<"\nSucess at " <<bin_5d;
			}
			//Scale Factor
			//std::cout<<"\nMaking a Scale Factor: ";
			if(fun::nSparseIntegral(_exp_data_5d[i][j]) >0.0 && fun::nSparseIntegral(_sim_data_5d[i][j]) >0.0){
				_scale_factor_5d[i][j]=fun::nSparseIntegral(_exp_data_5d[i][j])/fun::nSparseIntegral(_sim_data_5d[i][j]);
				//std::cout<<"Success at " <<bin_5d; 
			}else{
				_scale_factor_5d[i][j]=0.0;
				//std::cout<<"Zeroed at " <<bin_5d;
			}
			//Experimental correction 
			//std::cout<<"\nMaking Exp Corr Histogram: ";
			/*if(fun::nSparseIntegral(_acceptance[i][j]) > 0 && fun::nSparseIntegral(_sim_data_5d[i][j])>0.0){
				//_exp_corr_5d[i][j]->THnSparse::AddBinContent(bin_5d,_exp_data_5d[i][j]->THnSparse::GetBinContent(bin_5d));
				_exp_corr_5d[i][j]= (THnSparseD*)_exp_data_5d[i][j]->Clone();
				_exp_corr_5d[i][j]->Divide(_acceptance[i][j]);
				//std::cout<<"Sucess at " <<bin_5d;
			}*/
			//Experimental Hole estimate
			//std::cout<<"\nMaking Exp Holes: ";
			//std::cout<<"\nvar:" <<var <<" top:" <<top <<" W:" <<i <<" Q2:" <<j;
			//std::cout<<"\nSim_holes integral: " <<fun::nSparseIntegral(_sim_holes_5d[i][j]);
			//std::cout<<"\nScale Factor: " <<_scale_factor_5d[i][j];
			if(fun::nSparseIntegral(_sim_holes_5d[i][j])>0.0 && _scale_factor_5d[i][j]>0.0){
				//_exp_holes_5d[i][j]->THnSparse::AddBinContent(bin_5d,_sim_holes_5d[i][j]->THnSparse::GetBinContent(bin_5d));
				//std::cout<<"\nvar:" <<var <<" top:" <<top <<" W:" <<i <<" Q2:" <<j <<"| Int sim_holes:" <<fun::nSparseIntegral(_sim_holes_5d[i][j]);
				_exp_holes_5d[i][j]= (THnSparseD*)_sim_holes_5d[i][j]->Clone();
				_exp_holes_5d[i][j]->Scale(_scale_factor_5d[i][j]);
				_exp_corr_holes_5d[i][j] = (THnSparseD*)_exp_corr_5d[i][j]->Clone();
				//std::cout<<"\n\tExp Corr integral: " <<fun::nSparseIntegral(_exp_corr_5d[i][j]);
				//std::cout<<"\n\tExp holes integral: " <<fun::nSparseIntegral(_exp_holes_5d[i][j]);
				_exp_corr_holes_5d[i][j]->Add(_exp_holes_5d[i][j]);
				//std::cout<<"\n\tExp Corr +holes integral: " <<fun::nSparseIntegral(_exp_corr_holes_5d[i][j]);
				//std::cout<<"\n\tSuccess at " <<bin_5d;
			}
		}
	}
	std::cout<<"For Krishna: W:" <<_bin_low_7d[0][12] <<"-" <<_bin_up_7d[0][13] <<" Q2:" <<_bin_low_7d[1][2] <<"-" <<_bin_up_7d[1][2] <<" gives integral: " <<(fun::nSparseIntegral(_exp_data_5d[12][3])+fun::nSparseIntegral(_exp_data_5d[13][3])) <<"\n";
}
//For Single Differential bins
void Histogram::Sparse_5to3(Flags flags_){
	std::cout<<"Sparse 5 to 3\n";
	if(!flags_.Flags::Plot_Single_Diff()){
		std::cout<<"\tNot Plotting Single Differential Cross Sections\n";
	 	return;
	}
	//Convert the _exp_corr_5d and _exp_holes_5d to 3 dimensional sparse histograms for usage in plotting single differential cross sections
	//std::cout<<"\nPart 1";
	char hname[100];
	char xlabel[100];
	char ylabel[100];

	//Get the 7dimensional bins ready
	TH1D_1d_star exp_ch_1d;
	TH1D_2d_star exp_ch_2d;
	TH1D * exp_corr;
	TH1D * exp_holes;
	THnSparseD * exp_corr_holes_5d;
	//std::cout<<"\nPart 2";
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			for(int k=0; k<4; k++){
				//std::cout<<"\nPart 3 " <<k;
				sprintf(hname,"%s_yield_acc_corr_holes_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_five_dim_[k],_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
				sprintf(xlabel,"%s",_five_dim_[k]);
				sprintf(ylabel,"Yield");
				exp_ch_1d.push_back(_exp_corr_holes_5d[i][j]->Projection(k));
				//exp_ch_1d.push_back(_thrown_5d[i][j]->Projection(k));
				if(!flags_.Flags::Flux_Included()){//Scaling by the virtual photon flux at the center of the W Q2 bin if not already scaled
					exp_ch_1d[k]->Scale(1.0/physics::Virtual_Photon_Flux(_bin_mid_7d[0][i],_bin_mid_7d[1][j],_beam_energy_[0]));
				}
				//exp_ch_1d[k].Scale(1.0/flags_.Luminosity(0)); //Not yet implemented 12/6/22
				//exp_ch_1d[k].Scale(1.0/_Rad_Corr[i][j])//Not yet implemented 12/6/22
				exp_ch_1d[k]->Scale(1.0/(_bin_size_7d[0][i]*_bin_size_7d[1][j]));//Scaling by the size of the W Q^2 bins
				for(int l=0; l<_n_bins_5d[k]; l++){//Scaling by the size of the given variable bin size
					double prev_val = exp_ch_1d[k]->GetBinContent(l);
					exp_ch_1d[k]->SetBinContent(l,prev_val/_bin_size_5d[k][l]);
				}
				exp_ch_1d[k]->SetNameTitle(hname,hname);
				exp_ch_1d[k]->GetXaxis()->SetTitle(xlabel);
				exp_ch_1d[k]->GetYaxis()->SetTitle(ylabel);
			}
			//std::cout<<"\nPart 8 Q2:" <<j;
			exp_ch_2d.push_back(exp_ch_1d);
			exp_ch_1d.clear();
		}
		//std::cout<<"\nPart 9 W:" <<i;
		_exp_corr_holes_3d.push_back(exp_ch_2d);
		exp_ch_2d.clear();
	}
}

//For Polarization Observables
void Histogram::Sparse_5to4(Flags flags_){
	std::cout<<"Sparse 5 to 4\n";
	if(!flags_.Flags::Plot_Pol()){
		std::cout<<"\tNot Plotting Polarization Cross Sections\n";
	 	return;
	}
	char hname[100];
	char xlabel[100];
	char ylabel[100];
	TH1D_1d_star exp_ch_1d;
	TH1D_2d_star exp_ch_2d;
	TH1D_3d_star exp_ch_3d;
	TH2D * exp_corr_2deg;
	TH2D * exp_holes_2deg;
	TH1D * exp_corr;
	TH1D * exp_holes;
	Int_t proj;
	Int_t proj_bins[2];
	THnSparseD * exp_corr_holes_5d;
	TH2D* exp_corr_holes_2d;
	//Convert the _exp_corr_5d and _exp_holes_5d to 4 dimensional sparse histograms for usage in plotting single differential cross sections
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			for(int k=0; k<4; k++){//Xij
				//std::cout<<"For W: " <<i <<" Q2: " <<j <<" and dim: " <<_five_dim_[k] <<" the 5d integral is: " <<fun::nSparseIntegral(_exp_corr_holes_5d[i][j]) <<"\n";
				//std::cout<<"\tFor " <<_five_dim_[k] <<" the 2d-proj integral is " <<_exp_corr_holes_5d[i][j]->Projection(k,4)->Integral() <<"\n";
				//std::cout<<"\tNumber of bins in x " <<_exp_corr_holes_5d[i][j]->Projection(k,4)->GetNbinsX() <<" and in y " <<_exp_corr_holes_5d[i][j]->Projection(k,4)->GetNbinsY() <<"\n";
				for(int l=0; l<_n_bins_5d[k]; l++){//Bins of Xij
					//exp_corr_holes_5d = (THnSparseD*)_exp_corr_holes_5d[i][j]->Clone()->;
					//exp_corr_holes_2d = _exp_corr_holes_5d[i][j]->Projection(4,k);
					
					//sprintf(hname,"%s+%f",_five_dim_[k],l);
					//std::cout<<"\t\tFor " <<_five_dim_[k] <<" bin " <<l <<" the integral is " <<_exp_corr_holes_5d[i][j]->Projection(4,k)->ProjectionX("",l,l)->Integral() <<"\n";
					sprintf(xlabel,"Phi (deg)");
					sprintf(ylabel,"Yield");
					//exp_corr_holes_2d->SetNameTitle(hname,hname);
					//exp_ch_1d.push_back(_exp_corr_holes_5d[i][j]->Projection(4,k)->ProjectionX(_five_dim_[k],l,l));
					exp_ch_1d.push_back(_exp_corr_holes_5d[i][j]->Projection(k,4)->ProjectionX("",l,l,"e"));
					sprintf(hname,"%s:%.3f-%.3f_yield_corr_holes_W:%.3f-%.3f_Q2:%.2f-%.2f_top:%s_var:%s",_five_dim_[k],_bin_low_5d[k][l],_bin_up_5d[k][l],_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
					//std::cout<<"\t\t\tNaming histogram " <<hname <<"\n";
					if(!flags_.Flags::Flux_Included()){
						exp_ch_1d[l]->Scale(1.0/physics::Virtual_Photon_Flux(_bin_mid_7d[0][i],_bin_mid_7d[1][j],_beam_energy_[0]));//Needs way to set beam energy depending on which run 12/6/22
					}
					//exp_ch_1d[k]->Scale(1.0/flags_.Luminosity(0)); //Not yet implemented 12/6/22
					//exp_ch_1d[k]->Scale(1.0/_Rad_Corr[i][j])//Not yet implemented 12/6/22
					exp_ch_1d[l]->Scale(1.0/(_bin_size_7d[0][i]*_bin_size_7d[1][j]));//Scaling by the size of the W Q^2 bins
					exp_ch_1d[l]->Scale(1.0/_bin_size_5d[k][l]);//Scaling by the size of the given variable bin
					for(int m=0; m<_n_bins_5d[k]; m++){//Scaling by the size of the given phi bin
						double prev_val = exp_ch_1d[l]->GetBinContent(m);
						exp_ch_1d[l]->SetBinContent(m,prev_val/_bin_size_5d[4][m]);
					}
					exp_ch_1d[l]->SetNameTitle(hname,hname);
					exp_ch_1d[l]->GetXaxis()->SetTitle(xlabel);
					exp_ch_1d[l]->GetYaxis()->SetTitle(ylabel);
					//_exp_corr_4d[i][j][k][l]=Add(_exp_corr_4d[i][j][k][l],_exp_holes_4d[i][j][k][l]);
				}
				exp_ch_2d.push_back(exp_ch_1d);
				exp_ch_1d.clear();
			}
			exp_ch_3d.push_back(exp_ch_2d);
			exp_ch_2d.clear();
		}
		_exp_corr_holes_4d.push_back(exp_ch_3d);
		exp_ch_3d.clear();
	}
}

void Histogram::Beam_Spin(Flags flags_){
	//Will vary by W, Q2, and Phi 
	std::cout<<"Beam Spin Asymmetry\n";
	if(!flags_.Flags::Plot_Beam_Spin()){
		std::cout<<"\tNot Plotting Beam Spin Asymmetry\n";
	 	return;
	}
	char hname[100];
	char xlabel[100];
	char ylabel[100];
	TH1D_1d_star exp_ch_1d;
	TH1D_2d_star exp_ch_2d;
	TH1D_3d_star exp_ch_3d;
	TH2D * exp_corr_2deg;
	TH2D * exp_holes_2deg;
	TH1D * exp_corr;
	TH1D * exp_holes;
	Int_t proj;
	Int_t proj_bins[2];
	THnSparseD * exp_corr_holes_5d;
	TH2D* exp_corr_holes_2d;
	//Convert the _exp_corr_5d and _exp_holes_5d to 4 dimensional sparse histograms for usage in plotting single differential cross sections
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			sprintf(hname,"Beam_Spin_W|%f-%f_Q2|%f-%f",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j]);
			exp_ch_1d.push_back(new TH1D(hname,hname,_n_bins_5d[4],_bin_low_5d[4][0],_bin_up_5d[4][_n_bins_5d[4]-1]));
			for(int k=0; k<_n_bins_5d[4]; k++){//Phi
				//Needs to have the proper separation for postiive and negative helicities
				//exp_ch_1d->SetBinContent((1.0/_beam_pol_[0])*());//Beam polarization needs to be specified and is not currently accurate
			}
		}
	}
}

void Histogram::Extract_Bin_Info(Flags flags_){
	std::cout<<"Extract Bin Info\n";
	int DIM = -1;
	std::vector<int> n_bins_5d;
	/*std::vector<double> bin_low_set;
	std::vector<double> bin_up_set;
	std::vector<double> bin_mid_set;
	std::vector<double> bin_size_set;

	std::vector<std::vector<double>> bin_low_set_1;
	std::vector<std::vector<double>> bin_up_set_1;
	std::vector<std::vector<double>> bin_mid_set_1;
	std::vector<std::vector<double>> bin_size_set_1;

	std::vector<std::vector<std::vector<double>>> bin_low_set_2;
	std::vector<std::vector<std::vector<double>>> bin_up_set_2;
	std::vector<std::vector<std::vector<double>>> bin_mid_set_2;
	std::vector<std::vector<std::vector<double>>> bin_size_set_2;*/
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
		_n_bins_5d.push_back(_exp_data_7d->GetAxis(bin_fun+2)->GetNbins());;
	}
	//_n_bins_5d.push_back(n_bins_5d);
	//n_bins_5d.clear();
	DIM = _exp_data_7d->GetNdimensions();
	for(int i = 0; i<DIM; i++){
		_n_bins_7d.push_back(_exp_data_7d->GetAxis(i)->GetNbins());
		bin_lowside_7d_1d.push_back(_exp_data_7d->GetAxis(i)->GetBinLowEdge(1));
		bin_topside_7d_1d.push_back(_exp_data_7d->GetAxis(i)->GetBinUpEdge(_n_bins_7d[i]));
		for(int j=0; j<_n_bins_7d[i];j++){
			bin_low_7d_1d.push_back(_exp_data_7d->GetAxis(i)->GetBinLowEdge(j+1));
			bin_up_7d_1d.push_back(_exp_data_7d->GetAxis(i)->GetBinUpEdge(j+1));
			bin_mid_7d_1d.push_back((bin_up_7d_1d[j]+bin_low_7d_1d[j])/2.0);
			bin_size_7d_1d.push_back(bin_up_7d_1d[j]-bin_low_7d_1d[j]);
			bin_edges_7d_1d.push_back(_exp_data_7d->GetAxis(i)->GetBinLowEdge(j+1));
			if(j==_n_bins_7d[i]-1){
				bin_edges_7d_1d.push_back(_exp_data_7d->GetAxis(i)->GetBinUpEdge(j+1));
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
		//if(var==0 && top==0){
		//	_n_bins_5d.push_back(_exp_data_7d->GetAxis(j+2)->GetNbins());
			//std::cout<<"\nVar set: " <<var <<" num_bins: " <<_exp_data_7d->GetAxis(j+2)->GetNbins();
		//}
		bin_lowside_5d_1d.push_back(_exp_data_7d->GetAxis(j+2)->GetBinLowEdge(1));
		bin_topside_5d_1d.push_back(_exp_data_7d->GetAxis(j+2)->GetBinUpEdge(_n_bins_5d[j]));
		for(int k=0; k<_n_bins_5d[j];k++){
			bin_low_5d_1d.push_back(_exp_data_7d->GetAxis(j+2)->GetBinLowEdge(k+1));
			bin_up_5d_1d.push_back(_exp_data_7d->GetAxis(j+2)->GetBinUpEdge(k+1));
			bin_mid_5d_1d.push_back((bin_up_5d_1d[k]+bin_low_5d_1d[k])/2.0);
			bin_size_5d_1d.push_back(bin_up_5d_1d[k]-bin_low_5d_1d[k]);
			bin_edges_5d_1d.push_back(_exp_data_7d->GetAxis(j+2)->GetBinUpEdge(k+1));
			if(k==_n_bins_5d[j]-1){
				bin_edges_5d_1d.push_back(_exp_data_7d->GetAxis(j+2)->GetBinUpEdge(k+1));
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
	//std::cout<<"\nPart 2";
	int bin_sizes_5d[5];
	double bin_low_5d[5];
	double bin_up_5d[5];

	Sparse_3d_star exp_3d;
	Sparse_2d_star exp_2d;
	Sparse_1d_star exp_1d;
	Sparse_3d_star sim_3d;
	Sparse_2d_star sim_2d;
	Sparse_1d_star sim_1d;
	Sparse_3d_star exp_hole_3d;
	Sparse_2d_star exp_hole_2d;
	Sparse_1d_star exp_hole_1d;
	Sparse_1d_star exp_corr_hole_1d;
	Sparse_2d_star exp_corr_hole_2d;
	Sparse_3d_star exp_corr_hole_3d;
	Sparse_3d_star sim_hole_3d;
	Sparse_2d_star sim_hole_2d;
	Sparse_1d_star sim_hole_1d;
	Sparse_3d_star exp_corr_3d;
	Sparse_2d_star exp_corr_2d;
	Sparse_1d_star exp_corr_1d;
	Sparse_3d_star sim_corr_3d;
	Sparse_2d_star sim_corr_2d;
	Sparse_1d_star sim_corr_1d;
	Sparse_3d_star cross_3d;
	Sparse_2d_star cross_2d;
	Sparse_1d_star cross_1d;
	Sparse_3d_star accept_3d;
	Sparse_2d_star accept_2d;
	Sparse_1d_star accept_1d;
	Sparse_2d_star thr_2d;
	Sparse_1d_star thr_1d;
	Sparse_2d_star cs_5d_2d;
	Sparse_1d_star cs_5d_1d;

	double_3d scale_3d;
	double_2d scale_2d;
	double_1d scale_1d;

	/*Sparse_3d exp_set;
	Sparse_2d exp_set_1;
	Sparse_1d exp_set_2;
	Sparse_2d thrown_set;
	Sparse_1d thrown_set_1;*/
	
	for(int k=0; k<_n_bins_7d[0]; k++){//W	
		for(int l=0; l< _n_bins_7d[1]; l++){//Q2
			for(int boop=0; boop<5; boop++){//Loop over all bins
				bin_sizes_5d[boop] = _n_bins_5d[boop];//Get sizes for individual bins
				bin_low_5d[boop] = _bin_low_5d[boop][0];//Get low end for each bin
				bin_up_5d[boop] = _bin_up_5d[boop][_bin_up_5d[boop].size()-1];//Get high end for each bin
			}
			//exp
			sprintf(hname,"exp_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			//_exp_data_5d[j][k].push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			exp_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			//sim
			sprintf(hname,"sim_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			//_sim_data_5d[j][k].push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			sim_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			
			
			sprintf(hname,"sim_holes_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			//_sim_holes_5d[j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins,bin_low_5d,bin_up_5d));
			sim_hole_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));

			sprintf(hname,"exp_holes_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			//std::cout<<"\nexpHoles hname: " <<hname;
			//_exp_holes_5d[j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins,bin_low_5d,bin_up_5d));
			exp_hole_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			sprintf(hname,"exp_corr_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			//_exp_corr_5d[j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins,bin_low_5d,bin_up_5d));
			exp_corr_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));

			sprintf(hname,"exp_corr_hole_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			//_exp_corr_5d[j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins,bin_low_5d,bin_up_5d));
			exp_corr_hole_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));

			sprintf(hname,"sim_corr_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			//_sim_corr_5d[j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins,bin_low_5d,bin_up_5d));
			sim_corr_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));

			sprintf(hname,"cross_section_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			//_cross_section_5d[j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins,bin_low_5d,bin_up_5d));
			cross_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));

			sprintf(hname,"acceptance_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
			//_acceptance[j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins,bin_low_5d,bin_up_5d));
			accept_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			
			sprintf(hname,"thrown_5d_W:%.3f-%.3f_Q2:%.3f-%.3f_var:%s",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Var_Set().c_str());
			//_thrown_5d[l].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins,bin_low_5d,bin_up_5d));
			thr_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
			
			scale_1d.push_back(NAN);	
			//curr_hist=curr_hist+1;
			//percent=100.0*(curr_hist/num_hist);
			//std::cout<<"Skeleton_5D " <<percent <<"\r";
			//std::cout<<"\r" <<"\t" <<"Skeleton_5D " <<(100*curr_hist/num_hist) <<" %"  <<std::flush;
		}//Q2
		
		_exp_data_5d.push_back(exp_1d);
		_sim_data_5d.push_back(sim_1d);
		_sim_holes_5d.push_back(sim_hole_1d);
		_exp_holes_5d.push_back(exp_hole_1d);
		_exp_corr_holes_5d.push_back(exp_corr_hole_1d);
		_sim_corr_5d.push_back(sim_corr_1d);
		_exp_corr_5d.push_back(exp_corr_1d);
		_cross_section_5d.push_back(cross_1d);
		_acceptance_5d.push_back(accept_1d);
		_scale_factor_5d.push_back(scale_1d);
		_thrown_5d.push_back(thr_1d);
		
		exp_1d.clear();
		sim_1d.clear();
		sim_hole_1d.clear();
		exp_hole_1d.clear();
		exp_corr_1d.clear();
		exp_corr_hole_1d.clear();
		sim_corr_1d.clear();
		cross_1d.clear();
		accept_1d.clear();
		scale_1d.clear();
		thr_1d.clear();
		cs_5d_1d.clear();
		
		/*std::vector<THnSparseD> exp_set;
		std::vector<THnSparseD> sim_set;
		_exp_data_5d.push_back(exp_set);
		_sim_data_5d.push_back(sim_set);
		exp_set.clear();//std::vector::erase(exp_set.begin(),exp_set.begin()+exp_set.size());
		sim_set.clear();//std::vector::erase(sim_set.begin(),sim_set.begin()+sim_set.size());
		if(k==0){
			//thrown_set_1.push_back(thrown_set);
			//thrown_set.std::vector::erase(thrown_set.begin(),thrown_set.size());
		}*/
	}//W
	std::cout<<"End Skeleton_5D\n";
	//std::cout<<"\nCheck: " <<fun::nSparseIntegral(_exp_data_5d[0][0][0][0]) <<"\n";
}

void	Histogram::Calc_Acceptance(Flags flags_){
	std::cout<<"Calce Acceptance\n";
	/*std::vector<THnSparseD> accept_set;
	std::vector<THnSparseD> exp_corr_set;
	std::vector<THnSparseD> sim_corr_set;
	std::vector<double> n_thrown_set;
	std::vector<double> n_exp_set;
	std::vector<double> n_sim_set;
	std::vector<double> scale_set;*/
	std::cout<<"\nPart 1\n";
	Sparse_3d_star nsparse_set_3d;
	Sparse_2d_star nsparse_set_2d;
	Sparse_1d_star nsparse_set_1d;
	double_3d double_set_3d;
	double_2d double_set_2d;
	double_1d double_set_1d;

	Sparse_3d_star accept_3d;
	Sparse_3d_star exp_corr_3d;
	Sparse_3d_star sim_corr_3d;
	Sparse_2d_star accept_2d;
	Sparse_2d_star exp_corr_2d;
	Sparse_2d_star sim_corr_2d;
	Sparse_1d_star accept_1d;
	Sparse_1d_star exp_corr_1d;
	Sparse_1d_star sim_corr_1d;
	double_3d n_exp_3d;
	double_2d n_exp_2d;
	double_1d n_exp_1d;
	double_3d n_sim_3d;
	double_2d n_sim_2d;
	double_1d n_sim_1d;
	double_2d n_thr_2d;
	double_1d n_thr_1d;
	double_3d scale_3d;
	double_2d scale_2d;
	double_1d scale_1d;
	//THnSparseD empty = new THnSparseD("dummmy","dummy",5,{1,2,3,4,5},{1.0,1.0,1.0,1.0,1.0},{5.0,5.0,5.0,5.0,5.0}); 
	for(int k=0; k<_n_bins_7d[0]; k++){
		//std::cout<<"Part 4 idx: "<<i <<" " <<j <<" " <<k <<"\n";
		/*_acceptance[i][j].push_back(nsparse_set_1d);
		_exp_corr_5d[i][j].push_back(nsparse_set_1d);
		_sim_corr_5d[i][j].push_back(nsparse_set_1d);
		_n_exp_corr[i][j].push_back(double_set_1d);
		_n_sim_corr[i][j].push_back(double_set_1d);
		_scale_factor_5d[i][j].push_back(double_set_1d);
		if(j==0){
			_n_thrown[i].push_back(double_set_1d);
			_n_thrown[i][k].clear();
		}
		_acceptance[i][j][k].clear();
		_exp_corr_5d[i][j][k].clear();
		_sim_corr_5d[i][j][k].clear();
		_n_exp_corr[i][j][k].clear();
		_n_sim_corr[i][j][k].clear();
		_scale_factor_5d[i][j][k].clear();*/

		//std::cout<<"--Start Size of Vector Arrays--\nAccept: " <<_acceptance[i][j][k].size() <<"\n";
		for(int l=0; l<_n_bins_7d[1]; l++){
			//std::cout<<"Part 5 idx: "<<i <<" " <<j <<" " <<k <<" " <<l <<"\n";
			//std::cout<<"acceptance calc\n";
			if(fun::nSparseIntegral(_sim_data_5d[k][l]) > 0.0 && fun::nSparseIntegral(_thrown_5d[k][l]) > 0.0 ){//.at(i).at(j).at(k).at(l))
				std::cout<<"Sim Yield: " <<fun::nSparseIntegral(_sim_data_5d[k][l]) <<"\nThrown Yield: " <<fun::nSparseIntegral(_thrown_5d[k][l]) <<"\n";
				//std::cout<<"*-*idx_l = " <<l <<" | idx_size = " <<_acceptance[k].size() <<"\n";
				//_acceptance[k].push_back(_sim_data_5d[k][l]);
				accept_1d.push_back(_sim_data_5d[k][l]);
				//std::cout<<"*** Proper Fill *** idx_l = " <<l <<" | idx_size = " <<_acceptance[k].size() <<"\n";
				std::cout<<"*** Proper Fill *** idx_l = " <<l <<" | idx_size = " <<accept_1d.size() <<"\n";
				//std::cout<<"Acceptance Pre-Div Yield: " <<fun::nSparseIntegral(_acceptance[k][l]) <<"\n";
				//std::cout<<"Acceptance Pre-Div Yield: " <<fun::nSparseIntegral(accept_1d[l]) <<"\n";
				//_acceptance[k][l]->Divide(_thrown_5d[k][l]);
				accept_1d[l]->Divide(_thrown_5d[k][l]);
				//std::cout<<"Acceptance Divided Yield: " <<fun::nSparseIntegral(_acceptance[k][l]) <<"\n";
				std::cout<<"Acceptance Divided Yield: " <<fun::nSparseIntegral(accept_1d[l]) <<"\n";
			}else{
			
				//_acceptance[k].push_back(_sim_data_5d[k][l]);
				accept_1d.push_back(_sim_data_5d[k][l]);
				//std::cout<<"*** Empty Fill *** idx_l = " <<l <<" | idx_size = " <<_acceptance[k].size() <<"\n";
				std::cout<<"*** Empty Fill *** idx_l = " <<l <<" | idx_size = " <<accept_1d.size() <<"\n";
			}

			//std::cout<<"acceptance corr calc\n";
			//if(fun::nSparseIntegral(_acceptance[k][l]) >0.0){
			if(fun::nSparseIntegral(accept_1d[l]) >0.0){
				//std::cout<<"exp corr calc\n";
				if(fun::nSparseIntegral(_exp_data_5d[k][l]) > 0.0 ){
					//std::cout<<"Exp Yield: " <<fun::nSparseIntegral(_exp_data_5d[k][l]) <<"\n";
					//_exp_corr_5d[k].push_back(_exp_data_5d[k][l]);
					exp_corr_1d.push_back(_exp_data_5d[k][l]);
					//std::cout<<"*** Empty EXP Corr Fill *** idx_l = " <<l <<" | idx_size = " <<_exp_corr_5d[k].size() <<"\n";
					//std::cout<<"*** Proper EXP Corr Fill *** idx_l = " <<l <<" | idx_size = " <<exp_corr_1d.size() <<"\n";
					//_exp_corr_5d[k][l]->Divide(_acceptance[k][l]);
					//_n_exp_corr[k].push_back(fun::nSparseIntegral(_exp_corr_5d[k][l]));
					exp_corr_1d[l]->Divide(accept_1d[l]);
					n_exp_1d.push_back(fun::nSparseIntegral(exp_corr_1d[l]));
				}else{
					//_n_exp_corr[k].push_back(0.0);
					//_exp_corr_5d[k].push_back(_exp_data_5d[k][l]);
					n_exp_1d.push_back(0.0);
					exp_corr_1d.push_back(_exp_data_5d[k][l]);
				}
				//std::cout<<"sim corr calc\n";
				if(fun::nSparseIntegral(_sim_data_5d[k][l])>0.0){
					//_sim_corr_5d[k].push_back(_sim_data_5d[k][l]);
					//_sim_corr_5d[k][l]->Divide(_acceptance[k][l]);
					//_n_exp_corr[k].push_back(fun::nSparseIntegral(_sim_corr_5d[k][l]));
					sim_corr_1d.push_back(_sim_data_5d[k][l]);
					sim_corr_1d[l]->Divide(accept_1d[l]);
					n_sim_1d.push_back(fun::nSparseIntegral(sim_corr_1d[l]));
				}else{
					//_n_sim_corr[k].push_back(0.0);
					//_sim_corr_5d[k].push_back(_sim_data_5d[k][l]);
					n_sim_1d.push_back(0.0);
					sim_corr_1d.push_back(_sim_data_5d[k][l]);
				}
			}else{
				//_n_exp_corr[k].push_back(0.0);
				//_n_sim_corr[k].push_back(0.0);
				//_exp_corr_5d[k].push_back(_exp_data_5d[k][l]);
				//_sim_corr_5d[k].push_back(_sim_data_5d[k][l]);
				n_exp_1d.push_back(0.0);
				n_sim_1d.push_back(0.0);
				exp_corr_1d.push_back(_exp_data_5d[k][l]);
				sim_corr_1d.push_back(_sim_data_5d[k][l]);
			}
			//std::cout<<"Corr Yields\n" <<"exp: " <<_n_exp_corr[k][l] <<"| sim: " <<_n_sim_corr[k][l];
			std::cout<<"Corr Yields\n" <<"exp: " <<n_exp_1d[l] <<"| sim: " <<n_sim_1d[l];
			//if(_n_sim_corr[k][l]>0.0){
			if(n_sim_1d[l]>0.0){
				//_scale_factor_5d[k].push_back(_n_exp_corr[k][l]/_n_sim_corr[k][l]);
				scale_1d.push_back(n_exp_1d[l]/n_sim_1d[l]);

			}else{
				//_scale_factor_5d[k].push_back(0.0);
				scale_1d.push_back(0.0);
			}
			//std::cout<<" scale_factor: " <<_scale_factor_5d[k][l] <<"\n";
			std::cout<<" scale_factor: " <<scale_1d[l] <<"\n";
			/*accept_set.push_back(THnSparse::Divide(_sim_data_5d[k][l],_thrown_5d[k][l]));
			accept_set.push_back(_sim_data_5d[k][l]);
			accept_set.at(l).THnSparse::Divide(_thrown_5d[k][l]);
			exp_corr_set.push_back(THnSparse::Divide(_exp_data_5d[k][l],accept_set[l]));
			sim_corr_set.push_back(THnSparse::Divide(_sim_data_5d[k][l],accept_set[l]));
			n_exp_set.push_back(exp_corr_set->THnSparse::ComputeIntegral());
			n_sim_set.push_back(sim_corr_set->THnSparse::ComputeIntegral());
			scale_set.push_back(n_exp_set[l]/n_sim_set[l]);*/
			
			if(fun::nSparseIntegral(_thrown_5d[k][l]) > 0.0){
				//_n_thrown[k].push_back(fun::nSparseIntegral(_thrown_5d[k][l]));
				n_thr_1d.push_back(fun::nSparseIntegral(_thrown_5d[k][l]));
			}else{
				//_n_thrown[k].push_back(0.0);
				n_thr_1d.push_back(0.0);
			}
			
		}//l loop
		_acceptance_5d.push_back(accept_1d);
		_exp_corr_5d.push_back(exp_corr_1d);
		_sim_corr_5d.push_back(sim_corr_1d);
		_n_exp_corr.push_back(n_exp_1d);
		_n_sim_corr.push_back(n_sim_1d);
		_scale_factor_5d.push_back(scale_1d);
		_n_thrown.push_back(n_thr_1d);
		n_thr_1d.clear();
		accept_1d.clear();
		exp_corr_1d.clear();
		sim_corr_1d.clear();
		n_exp_1d.clear();
		n_sim_1d.clear();
		scale_1d.clear();
		
	}//k loop
}

void	Histogram::Calc_Holes_Sim(Flags flags_){
	std::cout<<"Calc Holes Sim\n";
	std::vector<THnSparseD> sim_hole_set;
	for(int k=0; k<_n_bins_7d[0]; k++){
		for(int l=0; l<_n_bins_7d[1]; l++){
			//sim_hole_set.push_back(_thrown_5d.at(i).at(k).at(l)-_sim_corr_5d.at(i).at(j).at(k).at(l));
			Sparse_Add_5d(_sim_holes_5d[k][l], _thrown_5d[k][l], _sim_corr_5d[k][l], -1);
		}
		//_sim_holes_5d[k].push_back(sim_hole_set);
		//sim_hole_set.std::vector::erase(sim_hole_set.begin(),sim_hole_set.size());
	}
}

void	Histogram::Calc_Holes_Exp(Flags flags_){
	std::cout<<"Calc Holes Exp\n";
	Sparse_1d_star exp_hole_1d;
	for(int k=0; k<_n_bins_7d[0]; k++){
		for(int l=0; l<_n_bins_7d[1]; l++){
			exp_hole_1d.push_back(_sim_holes_5d[k][l]);
			//_exp_holes_5d[k][l]=_sim_holes_5d[k][l];
			//exp_hole_set.push_back(_sim_holes_5d[k][l]);
			//exp_hole_set[l].Scale(_scale_factor_5d[k][l]);
			exp_hole_1d[l]->Scale(_scale_factor_5d[k][l]);
			//_exp_holes_5d[k][l]->Scale(_scale_factor_5d[k][l]);
		}
		//_exp_holes_5d.at(i).at(j).at(k).push_back(exp_hole_set);
		//exp_hole_set.std::vector::erase(exp_hole_set.begin(),exp_hole_set.size());
		_exp_holes_5d.push_back(exp_hole_1d);
		exp_hole_1d.clear();
	}
}

void	Histogram::Calc_Holes(Flags flags_){
	Histogram::Calc_Holes_Sim(flags_);
	Histogram::Calc_Holes_Exp(flags_);
}

void Histogram::Calc_Cross_Section(Flags flags_){
	std::cout<<"Calculate 7 dimensional Cross Section";
	//Get the 7dimensional bins ready
	int bin_5d[5];
	int bin_7d[7];

	long idx = 0;
	std::cout<<"\tcloning\n";
	TH2D_1d_star cs_2d_1d;
	TH2D_2d_star cs_2d_2d;
	//_cross_section_7d=(THnSparseD*)_exp_corr_7d->THnSparse::Clone();
	//_exp_corr_holes_5d[i][j] = (THnSparseD*)_exp_corr_5d[i][j]->Clone();
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		bin_7d[0]=i;
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			bin_7d[1]=j;
			for(int k = 0; k< _n_bins_7d[2]; k++){//MM1
				//_cross_section_5d[i][j]=(THnSparseD*)_exp_data_5d[i][j]->Clone();
				bin_5d[0]=k;
				bin_7d[2]=k;
				for(int l = 0; l< _n_bins_7d[3]; l++){//MM2
					bin_5d[1]=l;
					bin_7d[3]=l;
					for(int m = 0; m< _n_bins_7d[4]; m++){//Theta
						bin_5d[2]=m;
						bin_7d[4]=m;
						for(int n = 0; n< _n_bins_7d[5]; n++){//Alpha
							bin_5d[3]=n;
							bin_7d[5]=n;
							for(int o = 0; o< _n_bins_7d[6]; o++){//Phi
								bin_5d[4]=o;
								//std::cout<<"\t\tAdding bin content: " <<bin_7d <<" | " <<bin_5d <<"\n";
								bin_7d[6]=o;
								//_cross_section_7d->THnSparse::AddBinContent(bin_7d,_exp_corr_7d->THnSparse::GetBinContent(bin_7d));//Should replace this and empty piece with empty subtracted then acceptance divided
								//_cross_section_7d->THnSparse::AddBinContent(bin_7d,flags_.Flags::Qr()*_empty_7d->THnSparse::GetBinContent(bin_7d));
								//_cross_section_7d->THnSparse::AddBinContent(bin_7d,_scale_factor_7d*_exp_holes_7d->THnSparse::GetBinContent(bin_7d));
								
								_cross_section_5d[i][j]->THnSparse::AddBinContent(bin_5d,_exp_data_5d[i][j]->THnSparse::GetBinContent(bin_5d));//Should replace this and empty piece with empty subtracted then acceptance divided
								//_cross_section_7d->THnSparse::AddBinContent(bin_7d,flags_.Flags::Qr()*_empty_7d->THnSparse::GetBinContent(bin_7d));
								_cross_section_5d[i][j]->THnSparse::AddBinContent(bin_5d,_scale_factor_5d[i][j]*_exp_holes_5d[i][j]->THnSparse::GetBinContent(bin_5d));
							}
						}
					}
				}
			}
			for(int p=0; i<4; i++){
				cs_2d_1d.push_back(_cross_section_5d[i][j]->THnSparse::Projection(p,4));
			}
			cs_2d_2d.push_back(cs_2d_1d);
			cs_2d_1d.clear();
		}
		_cross_section_2d.push_back(cs_2d_2d);
		cs_2d_2d.clear();
	}
}

void	Histogram::Make_Single_Diff(Flags flags_){
	std::cout<<"Making Single Differential";
	if(!flags_.Flags::Plot_Single_Diff()){
		std::cout<<" | Not plotting single differential cross sections\n";
		return;
	}
	char place[100];
	//Making the histograms
	TH1D_1d_star single_1d;
	TH1D_2d_star single_2d;
	for(int i=0; i<_n_bins_7d[0]; i++){
		for(int j=0; j<_n_bins_7d[1]; j++){
			for(int k=0; k<5; k++){
				sprintf(place,"Projection in W(%.3f-%.3f) Q2(%.3f-%.3f)",_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j]);
				std::cout<<"\t" <<place <<" with integral of: " <<fun::nSparseIntegral(_cross_section_5d[i][j]) <<"\n";
				single_1d.push_back(_cross_section_5d[i][j]->Projection(k));
			}
			single_2d.push_back(single_1d);
			single_1d.clear();
		}
		_single_diff_hist.push_back(single_2d);
		single_2d.clear();
	}


	std::cout<<"\tMaking Directory\n";
	TDirectory* dir_S = _RootOutputFile->mkdir("Single Differential CS");
	dir_S->cd();
	//gDirectory->pwd();
	TDirectory* dir_S1[5];
	TDirectory* dir_S2[5][_n_bins_7d[1]];
	char dirname[100];
	char hname[100];
	for(int k=0; k<5; k++){
		sprintf(dirname,"%s",_five_dim_[k]);
		dir_S1[k] = dir_S->mkdir(dirname);
		dir_S1[k]->cd();
		for(int j=0; j<_n_bins_7d[1]; j++){
			sprintf(dirname,"%s_Q2|%.3f-%.3f",_five_dim_[k],_bin_low_7d[1][j],_bin_up_7d[1][j]);
			dir_S2[k][j] = dir_S1[k]->mkdir(dirname);
			dir_S2[k][j]->cd();
			for(int i=0; i<_n_bins_7d[0]; i++){
				//std::cout<<"\t\tTrying to make projection in final directory\n\t\t\tk= " <<k <<" j= " <<j <<" i= " <<i <<"\n";
				//gDirectory->pwd();
				//_single_diff_hist[i][j][k]=_cross_section_5d[i][j]->Projection(k);//Oh, I haven't built out the single_diff at all
				sprintf(hname,"Single_Differential_Cross_Section_%s_W=(%.3f-%.3f)_Q2=(%.3f-%.3f)",_five_dim_[k],_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j]);
				//std::cout<<"\t\t\t" <<hname <<"\n";
				_single_diff_hist[i][j][k]->SetTitle(hname);
				//std::cout<<"\t\t\tAssigned title: " <<hname <<"\n"; 
				_single_diff_hist[i][j][k]->SetXTitle(_five_dim_[k]);
				//std::cout<<"\t\t\tLabeled X axis: " <<_five_dim_[k] <<"\n"; 
				_single_diff_hist[i][j][k]->SetYTitle("Differential Cross Section (microbarns/unit)");
				//std::cout<<"\t\t\tLabeled Y axis: " <<"Differential Cross Section (microbarns/unit)" <<"\n"; 
				_single_diff_hist[i][j][k]->Write();
				//std::cout<<"\t\t\tHistogram Written\n";
			}
		}
	}
}

void	Histogram::Make_Polarization(Flags flags_){
	//In no way complete 10/4/22
	std::cout<<"Make Polarization Histograms";
	if(!flags_.Flags::Plot_Pol()){
		std::cout<<" | Not plotting Polarization cross sections\n";
		return;
	}
	TDirectory* dir_P = _RootOutputFile->mkdir("Polarization CS");
	dir_P->cd();
	TDirectory* dir_P1[4];
	TDirectory* dir_P2[5][_n_bins_7d[1]];
	char dirname[100];
	char hname[100];
	for(int k=0; k<4; k++){
		sprintf(dirname,"%s",_five_dim_[k]);
		dir_P1[k] = dir_P->mkdir(dirname);
		dir_P1[k]->cd();
		for(int j=0; j<_n_bins_7d[1]; j++){
			sprintf(dirname,"%s_Q2|%.3f-%.3f",_five_dim_[k],_bin_low_7d[1][j],_bin_up_7d[1][j]);
			dir_P2[k][j] = dir_P->mkdir(dirname);
			dir_P2[k][j]->cd();
			for(int i=0; i<_n_bins_7d[0]; i++){
				_single_diff_hist[i][j][k]=_cross_section_5d[i][j]->Projection(k);
				sprintf(hname,"Polarization_Cross_Section_%s_W|%.3f-%.3f_Q2|%.3f-%.3f",_five_dim_[k],_bin_low_7d[0][i],_bin_up_7d[0][i],_bin_low_7d[1][j],_bin_up_7d[1][j]);
				_single_diff_hist[i][j][k]->SetTitle(hname);
				_single_diff_hist[i][j][k]->SetXTitle(_five_dim_[k]);
				_single_diff_hist[i][j][k]->SetXTitle("Polarization Cross Section (microbarns/unit)");
				_single_diff_hist[i][j][k]->Write();
			}
		}
	}
}

void	Histogram::Make_Integrated(Flags flags_){
	std::cout<<"Calc";
}

void Histogram::Make_WQ2(Flags flags_){
	std::cout<<"\tMaking WQ2 Plots\n";
	if(!flags_.Flags::Plot_WQ2()){
		std::cout<<"\t\tNot making WQ2 Plots\n";
		return;
		std::cout<<"\t\tProperly exited...?\n";
	}
	char hname[100];
	sprintf(hname,"Exp_WQ2_Top:%s_Var:%s",flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
	_exp_hist_wq2 = _exp_data_7d->Projection(1,0);
	std::cout<<"\t\tExperimental WQ2 total yield: " <<_exp_hist_wq2->Integral() <<"\n";
	_exp_hist_wq2->SetNameTitle(hname,hname);
	sprintf(hname,"Sim_WQ2_Top:%s_Var:%s",flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
	_sim_hist_wq2 = _sim_data_7d->Projection(1,0);
	std::cout<<"\t\tSimulated WQ2 total yield: " <<_sim_hist_wq2->Integral() <<"\n";
	_sim_hist_wq2->SetNameTitle(hname,hname);
	sprintf(hname,"Exp_Corr_WQ2_Top:%s_Var:%s",flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
	_exp_corr_hist_wq2 = _exp_corr_7d->Projection(1,0);
	std::cout<<"\t\tExperimental acceptance corrected WQ2 total yield: " <<_exp_corr_hist_wq2->Integral() <<"\n";
	_exp_corr_hist_wq2->SetNameTitle(hname,hname);
	sprintf(hname,"Sim_Corr_WQ2_Top:%s_Var:%s",flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
	_sim_corr_hist_wq2 = _sim_corr_7d->Projection(1,0);
	std::cout<<"\t\tSim acceptance corrected WQ2 total yield: " <<_sim_corr_hist_wq2->Integral() <<"\n";
	_sim_corr_hist_wq2->SetNameTitle(hname,hname);
	sprintf(hname,"Thr_WQ2_Var:%s",flags_.Flags::Var_Set().c_str());
	_thr_hist_wq2 = _thrown_7d->Projection(1,0);
	std::cout<<"\t\tThrown WQ2 total yield: " <<_thr_hist_wq2->Integral() <<"\n";
	_thr_hist_wq2->SetNameTitle(hname,hname);
	std::cout<<"\tWriting W vs Q2 Plots\n";
	if(!flags_.Plot_WQ2()){
		std::cout<<"\tNot Plotting WQ2\n";
	 	return;
	 	std::cout<<"\t\tProperly exited...?\n";
	}else{
		std::cout<<"\t\tMaking Directory\n";
		TDirectory* dir_WQ2 = _RootOutputFile->mkdir("WQ2 Plots");
		dir_WQ2->cd();
		std::cout<<"\t\tWriting in Directory\n";
		std::cout<<"\t\tWriting exp hist\n";
		_exp_hist_wq2->SetXTitle("W (GeV)");
		_exp_hist_wq2->SetYTitle("Q^{2} (GeV^{2}");
		_exp_hist_wq2->SetOption("Colz");
		_exp_hist_wq2->Write();
		std::cout<<"\t\tWriting sim hist\n";
		_sim_hist_wq2->SetXTitle("W (GeV)");
		_sim_hist_wq2->SetYTitle("Q^{2} (GeV^{2}");
		_sim_hist_wq2->SetOption("Colz");
		_sim_hist_wq2->Write();
		std::cout<<"\t\tWriting exp acceptance corrected hist\n";
		_exp_corr_hist_wq2->SetXTitle("W (GeV)");
		_exp_corr_hist_wq2->SetYTitle("Q^{2} (GeV^{2}");
		_exp_corr_hist_wq2->SetOption("Colz");
		_exp_corr_hist_wq2->Write();
		std::cout<<"\t\tWriting sim acceptance corrected hist\n";
		_sim_corr_hist_wq2->SetXTitle("W (GeV)");
		_sim_corr_hist_wq2->SetYTitle("Q^{2} (GeV^{2}");
		_sim_corr_hist_wq2->SetOption("Colz");
		_sim_corr_hist_wq2->Write();
	}
}

//Plots showing the number of empty reconstructed simulation bins where there is data from experiment
//This doesn't work the way I had very much hoped it would work. Cannot make this histogram as intended here
/*
void Histogram::Make_Acceptance_Statistics(Flags flags_){
	std::cout<<"\tMake Acceptance Statistic Plot? ";
	if(flags_.Flags::Plot_Acceptance()){
		std::cout<<"Yes\n";
		char hname[100];
		sprintf(hname,"Determination of Proper Simulation Statistics");
		_acc_zero_exp = new TH1D(hname,hname,200,0.0,1.0);
		int bin_7d[7];
		int sim_bins_total = _sim_data_7d->THnSparse::GetNbins();
		int sim_bin = 0;
		float frac = (float)sim_bin/(float)sim_bins_total;
		std::cout<<"Total Reconstructed Bins: " <<sim_bins_total <<"\n";
		for(int i=0; i<_n_bins_7d[0]; i++){//W
			bin_7d[0]=i;
			for(int j=0; j< _n_bins_7d[1]; j++){//Q2
				bin_7d[1]=j;
				for(int k = 0; k< _n_bins_7d[2]; k++){//MM1
					bin_7d[2]=k;
					for(int l = 0; l< _n_bins_7d[3]; l++){//MM2
						bin_7d[3]=l;
						for(int m = 0; m< _n_bins_7d[4]; m++){//Theta
							bin_7d[4]=m;
							for(int n = 0; n< _n_bins_7d[5]; n++){//Alpha
								bin_7d[5]=n;
								for(int o = 0; o< _n_bins_7d[6]; o++){//Phi
									bin_7d[6]=o;
									if(_sim_data_7d->THnSparse::GetBinContent(bin_7d)>0.0){
										sim_bin++;
										frac = (float)sim_bin/(float)sim_bins_total;
									}
									if(_exp_data_7d->THnSparse::GetBinContent(bin_7d) > 0.0){
										if(_sim_data_7d->THnSparse::GetBinContent(bin_7d)==0.0){
											_acc_zero_exp->Fill(frac);
										}		
									}
								}
							}
						}
					}
				}
			}
		}
		TDirectory* dir_acc_zero = _RootOutputFile->mkdir("Acceptance Statistics");
		dir_acc_zero->cd();
		_acc_zero_exp->SetXTitle("Fraction of Simulation Performed");
		_acc_zero_exp->SetYTitle("Number of empty Rec Bins in Filled Exp bins");
		_acc_zero_exp->Write();
	}
}*/


//These are the 1 dimensional Acceptance plots that vary by W and Q2
void Histogram::Make_Acceptance(Flags flags_){
	if(flags_.Flags::Plot_Acceptance()){
		std::cout<<"Making Acceptance Plots\n";
		std::cout<<"\tMaking Acceptance Directory\n";
		TDirectory* dir_A = _RootOutputFile->mkdir("Acceptance Plots");
		dir_A->cd();
		std::cout<<"\t\tMaking sub directories\n";
		TDirectory* dir_A_1[_n_bins_7d[1]];
		TDirectory* dir_A_2[_n_bins_7d[1]][_n_bins_7d[0]];
		//TH1D_5d _accept_hist_1;//(var,top,w,q2,{MM1,MM2,theta,alpha,phi})
		//TH1D_3d _accept_hist_2;//(var,top,{w,q2})
		//std::cout<<"\nMaking WQ2 Plots";
		char hname[100];
		char dir_name[100];
		TH1D_1d_star accept_hist_01;
		TH1D_2d_star accept_hist_02;
		//TH1D* accept_hist_00;
		for(int k=0; k<_n_bins_7d[0]; k++){
			/*std::cout<<"\t\tIterating by Q2. Bin:" <<l <<"\n";
			sprintf(dir_name,"Acceptance_Q2:%.3f-%.3f",_bin_low_7d[1][l],_bin_up_7d[1][l]);
			std::cout<<"\t\tdir_name: " <<dir_name <<"\n";
			dir_A_1[l] = dir_A->mkdir(dir_name);*/
			//_accept_hist_1.push_back(accept_hist_03);
			for(int l=0; l<_n_bins_7d[1]; l++){
				//std::cout<<"\t\t\tacceptance integral: " <<fun::nSparseIntegral(_acceptance_5d[k][l]) <<" for W Q2 bins: " <<k <<" " <<l <<"\n";
				//std::cout<<"\t\t\tQ2 bin:" <<l <<" W bin:" <<k <<"\n";
				//sprintf(dir_name,"Acceptance_W:%.3f-%.3f_Q2:%.3f-%.3f",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l]);
				//std::cout<<"\t\t\tdirectory name: " <<dir_name <<"\n";		
				//_accept_hist_1[k].push_back(accept_hist_04);
				//std::cout<<"\t\t\t\t5d Acceptance isn't empty\n";
				//dir_A_2[l][k] = dir_A_1[l]->mkdir(dir_name);
				//dir_A_2[l][k]->cd();
				//std::cout<<"\t\t\t\tEntered Directory for saving\n";
				for(int m=0; m<5; m++){
					//std::cout<<"\nstep 1";
					//sprintf(hname,"acceptance_%s_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_set:%s_",_five_dim_[m],_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
					//std::cout<<"\tGetting hist: " <<hname <<"\n";
					accept_hist_01.push_back(_acceptance_5d[k][l]->Projection(m));
					//accept_hist_00 = _acceptance_5d[k][l]->Projection(m);
					//std::cout<<"\t\tWriting 1D Acceptance for W:" <<k <<" Q2:" <<l <<" X:" <<m <<"\n";
					//accept_hist_00->SetTitle(hname);
					//accept_hist_00->SetYTitle("Integrated Acceptance");
					//accept_hist_00->SetXTitle(_five_dim_[m]);
					//accept_hist_00->Write();
				}
				//std::cout<<"\nstep 2";
				accept_hist_02.push_back(accept_hist_01);
				accept_hist_01.clear();
			}
			//std::cout<<"\nstep 3";
			_accept_hist_1.push_back(accept_hist_02);
			accept_hist_02.clear();
		}
		_accept_hist_2[0]=(TH1D*)_acceptance_7d->Projection(0)->Clone();
		_accept_hist_2[1]=(TH1D*)_acceptance_7d->Projection(1)->Clone();
		std::cout<<"\tWriting W Q2 projections of acceptance\n";
		dir_A->cd();
		_accept_hist_2[0]->SetXTitle("W (GeV)");
		_accept_hist_2[0]->SetYTitle("Integrated Acceptance");
		_accept_hist_2[0]->Write();
		_accept_hist_2[1]->SetXTitle("Q^{2} (GeV^{2})");
		_accept_hist_2[1]->SetYTitle("Integrated Acceptance");
		_accept_hist_2[1]->Write();

		//TString* path; 
		//char path[100];
		for(int l=0; l<_n_bins_7d[1]; l++){
			std::cout<<"\t\tIterating by Q2. Bin:" <<l <<"\n";
			sprintf(dir_name,"Acceptance_Q2:%.3f-%.3f",_bin_low_7d[1][l],_bin_up_7d[1][l]);
			std::cout<<"\t\t" <<dir_name <<"\n";
			dir_A_1[l] = dir_A->mkdir(dir_name);
			dir_A_1[l]->cd();
			gDirectory->pwd();
			for(int k =0; k<_n_bins_7d[0]; k++){
				if(fun::nSparseIntegral(_acceptance_5d[k][l])>0.0){
					std::cout<<"\t\t\tacceptance integral: " <<fun::nSparseIntegral(_acceptance_5d[k][l]) <<" for W Q2 bins: " <<k <<" " <<l <<"\n";
					std::cout<<"\t\t\tQ2 bin:" <<l <<" W bin:" <<k <<"\n";
					sprintf(dir_name,"Acceptance_W:%.3f-%.3f_Q2:%.3f-%.3f",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l]);
					std::cout<<"\t\t\tMaking dir:" <<dir_name <<"\n";
					dir_A_2[l][k] = dir_A_1[l]->mkdir(dir_name);
					//std::cout<<"\t\tdirectory contents " <<dir_A_1[l]->Print() <<"\n";
					dir_A_2[l][k]->cd();
					dir_A_2[l][k]->Print();
					//std::cout<<"\t\t pwd: " <<dir_A_2[l][k]->TDirectory::pwd() <<"\n"; 
					gDirectory->pwd();
					//gDirectory->GetDirectory();
					//std::cout<<"\t\t\t\tlook at full path: " <<path <<"\n";
					std::cout<<"\t\t\tEntered Directory for saving " <<dir_name <<"\n";
					for(int m=0; m<5; m++){
						std::cout<<"\t\t\t\tWriting 1D Acceptance for W:" <<k <<" Q2:" <<l <<" X:" <<m <<"\n";
						sprintf(hname,"Acceptance_%s_W:%.3f-%.3f_Q2:%.3f-%.3f_top:%s_set:%s_",_five_dim_[m],_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l],flags_.Flags::Top().c_str(),flags_.Flags::Var_Set().c_str());
						_accept_hist_1[k][l][m]->SetTitle(hname);
						_accept_hist_1[k][l][m]->SetYTitle("Integrated Acceptance");
						_accept_hist_1[k][l][m]->SetXTitle(_five_dim_[m]);
						_accept_hist_1[k][l][m]->Write();
					}
				}
			}
		}
	}
}
	//void Write_5d_Yield();
	//void Write_5d_Cross_Section();
	//void Write_5d_Holes();
	//void Write_7d_Holes();
	//void Write_Single_Diff();
	//void Write_Polarization();
	//void Write_Integrated();
void Histogram::Convert_to_Cross(Flags flags_){
	std::cout<<"Convert to Cross\n";
	std::cout<<"\tMaking X and Phi histograms with their bin widths\n";
	Histogram::XandPhi_BinHistograms(flags_);
	Histogram::Convert_Single_Diff_to_Cross(flags_);
	Histogram::Convert_Polarization_to_Cross(flags_);
}

void Histogram::Convert_Single_Diff_to_Cross(Flags flags_){
	std::cout<<"\tConvert Single Diff to Cross\n";
	if(!flags_.Flags::Plot_Single_Diff()){
		std::cout<<"\tNot Plotting Single Differentials\n";
	 	return;
	 	std::cout<<"\t\tProperly exited...?\n";
	}
	for(int k=0; k<_n_bins_7d[0]; k++){//W
		//_accept_hist_1.push_back(accept_hist_03);
		for(int l=0; l<_n_bins_7d[1]; l++){//Q2
			for(int m=0; m<4; m++){//Xij
				_exp_corr_holes_3d[k][l][m]->Scale(1./physics::Luminosity(_Q_tot_));
				_exp_corr_holes_3d[k][l][m]->Scale(1./physics::Virtual_Photon_Flux(Histogram::Bin_Center_7d(0,k),Histogram::Bin_Center_7d(0,k),_beam_energy_[flags_.Flags::Run()]));
				//_exp_corr_holes_3d[k][l][m]->Scale(1./physics::Radiative_corr()); //We don't have this correction yet
				_exp_corr_holes_3d[k][l][m]->Scale(1./Histogram::W_Bin_Size(k));
				_exp_corr_holes_3d[k][l][m]->Scale(1./Histogram::Q2_Bin_Size(l));
				//if(i==j && j==k && k==l && l==0){
					//std::cout<<"\t\texp corr holes number of x bins:" <<_exp_corr_holes_3d[k][l][m]->GetNbinsX() <<"\n";
					//std::cout<<"\t\texp corr holes low edge:" <<_exp_corr_holes_3d[k][l][m]->GetXaxis()->GetBinLowEdge(0)<<"\n";
					//std::cout<<"\t\texp corr holes top edge:" <<_exp_corr_holes_3d[k][l][m]->GetXaxis()->GetBinUpEdge(_exp_corr_holes_3d[k][l][m]->GetNbinsX()-1)<<"\n";
					//std::cout<<"\t\tX bin size num x bins:" <<_X_bin_sizes[i][m]->GetNbinsX() <<"\n";
					//std::cout<<"\t\tX bin low edge:" <<_X_bin_sizes[i][m]->GetXaxis()->GetBinLowEdge(0) <<"\n";
					//std::cout<<"\t\tX bin top edge:" <<_X_bin_sizes[i][m]->GetXaxis()->GetBinUpEdge(_X_bin_sizes[i][m]->GetNbinsX()-1) <<"\n\n";
				//}
				if(_exp_corr_holes_3d[k][l][m]->GetNbinsX() == _X_bin_sizes[m]->GetNbinsX()){
					_exp_corr_holes_3d[k][l][m]->Divide(_X_bin_sizes[m]);
				}else{
					std::cout<<"\t\t\tBin sizes did not match: exp= " <<_exp_corr_holes_3d[k][l][m]->GetNbinsX() <<" X= "<<_X_bin_sizes[m]->GetNbinsX() <<"\n";
					std::cout<<"\t\t\tLow Ends: exp= " <<_exp_corr_holes_3d[k][l][m]->GetXaxis()->GetBinLowEdge(0) <<" X= "<<_X_bin_sizes[m]->GetXaxis()->GetBinLowEdge(0) <<"\n";
					std::cout<<"\t\t\tHigh Ends: exp= " <<_exp_corr_holes_3d[k][l][m]->GetXaxis()->GetBinUpEdge(_exp_corr_holes_3d[k][l][m]->GetNbinsX()-1) <<" X= "<<_X_bin_sizes[m]->GetXaxis()->GetBinUpEdge(_X_bin_sizes[m]->GetNbinsX()-1) <<"\n";
				}
			}
		}
	}
}	

void Histogram::Convert_Polarization_to_Cross(Flags flags_){
	std::cout<<"\tConvert Polarization to Cross\n";
	if(!flags_.Flags::Plot_Pol()){
		std::cout<<"\t\tNot Plotting Polarization Cross Sections\n";
	 	return;
	 	std::cout<<"\t\tProperly exited...?\n";
	}
	for(int k=0; k<_n_bins_7d[0]; k++){//W
		//_accept_hist_1.push_back(accept_hist_03);
		for(int l=0; l<_n_bins_7d[1]; l++){//Q2
			for(int m=0; m<4; m++){//Xij
				for(int n=0; n<_n_bins_5d[m]; n++){
					_exp_corr_holes_4d[k][l][m][n]->Scale(1./physics::Luminosity(_Q_tot_));
					_exp_corr_holes_4d[k][l][m][n]->Scale(1./physics::Virtual_Photon_Flux(Histogram::Bin_Center_7d(0,k),Histogram::Bin_Center_7d(0,k),_beam_energy_[flags_.Flags::Run()]));
					//_exp_corr_holes_4d[k][l][m][n]->Scale(1./physics::Radiative_corr()); //We don't have this correction yet
					_exp_corr_holes_4d[k][l][m][n]->Scale(1./Histogram::W_Bin_Size(k));
					_exp_corr_holes_4d[k][l][m][n]->Scale(1./Histogram::Q2_Bin_Size(l));
					_exp_corr_holes_4d[k][l][m][n]->Scale(1./Histogram::X_Bin_Size(m,n));
					if(_exp_corr_holes_4d[k][l][m][n]->GetNbinsX() == _phi_bin_sizes->GetNbinsX()){
						_exp_corr_holes_4d[k][l][m][n]->Divide(_phi_bin_sizes);
					}else{
						std::cout<<"\t\t\t\tIndex: " <<" " <<k <<" " <<l <<" " <<m <<" " <<n <<"\n";
						std::cout<<"\t\t\tBin sizes did not match: exp= " <<_exp_corr_holes_4d[k][l][m][n]->GetNbinsX() <<" X= "<<_phi_bin_sizes->GetNbinsX() <<"\n";
						std::cout<<"\t\t\tLow Ends: exp= " <<_exp_corr_holes_4d[k][l][m][n]->GetXaxis()->GetBinLowEdge(0) <<" X= "<<_phi_bin_sizes->GetXaxis()->GetBinLowEdge(0) <<"\n";
						std::cout<<"\t\t\tHigh Ends: exp= " <<_exp_corr_holes_4d[k][l][m][n]->GetXaxis()->GetBinUpEdge(_exp_corr_holes_4d[k][l][m][n]->GetNbinsX()-1) <<" X= "<<_phi_bin_sizes->GetXaxis()->GetBinUpEdge(_X_bin_sizes[m]->GetNbinsX()-1) <<"\n";
					}
				}
			}
		}
	}
}

//It appears That I cannot access any of these histograms at this point for some reason, so I'm moving the writing to where I make them. Will address this issue at a later date 9/23/22
void Histogram::Write_WQ2(Flags flags_){
	std::cout<<"Writing W vs Q2 Plots\n";
	if(!flags_.Plot_WQ2()){
		std::cout<<"\tNot Plotting WQ2\n";
	 	return;
	 	std::cout<<"\t\tProperly exited...?\n";
	}
	std::cout<<"\tMaking Directory\n";
	TDirectory* dir_WQ2 = _RootOutputFile->mkdir("WQ2 Plots");
	dir_WQ2->cd();
	std::cout<<"\tWriting in Directory\n";
	std::cout<<"\t\t\tChecking exp hist yield: " <<_exp_hist_wq2->Integral() <<"\n";
	std::cout<<"\t\tWriting exp hist\n";
	_exp_hist_wq2->SetXTitle("W (GeV)");
	_exp_hist_wq2->SetYTitle("Q^{2} (GeV^{2}");
	_exp_hist_wq2->SetOption("Colz");
	_exp_hist_wq2->Write();
	std::cout<<"\t\tWriting sim hist\n";
	_sim_hist_wq2->Write();
	std::cout<<"\t\tWriting exp acceptance corrected hist\n";
	_exp_corr_hist_wq2->Write();
	std::cout<<"\t\tWriting sim acceptance corrected hist\n";
	_sim_corr_hist_wq2->Write();
}

void Histogram::Write_Acceptance(Flags flags_){
	std::cout<<"Writing Acceptance Plots\n";
	if(!flags_.Flags::Plot_Acceptance()){
		std::cout<<"\tNot Plotting Acceptance\n";
		return;
		std::cout<<"\t\tProperly exited...?\n";
	}
	std::cout<<"\tMaking directory\n";
	TDirectory* dir_A = _RootOutputFile->mkdir("Acceptance Plots");
	dir_A->cd();
	std::cout<<"\t\tPart 1\n";
	TDirectory* dir_A_1[_n_bins_7d[1]];
	TDirectory* dir_A_2[_n_bins_7d[1]][_n_bins_7d[0]];
	char dir_name[100];
	for(int l=0; l<_n_bins_7d[1]; l++){
		std::cout<<"\t\tPart 4\n";
		sprintf(dir_name,"Acceptance_Q2:%.3f-%.3f",_bin_low_7d[1][l],_bin_up_7d[1][l]);
		std::cout<<"\t\t\t" <<dir_name <<"\n";
		dir_A_1[l] = dir_A->mkdir(dir_name);
		for(int k =0; k<_n_bins_7d[0]; k++){
			if(fun::nSparseIntegral(_acceptance_5d[k][l])>0.0){
				std::cout<<"acceptance integral: " <<fun::nSparseIntegral(_acceptance_5d[k][l]) <<" for W Q2 bins: " <<k <<" " <<l <<"\n";
				std::cout<<"\t\tPart 5\n";
				//std::cout<<"\t\t index:" <<i <<" " <<j <<" " <<k <<" " <<l <<"\n";
				sprintf(dir_name,"Acceptance_W:%.3f-%.3f_Q2:%.3f-%.3f",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l]);
				std::cout<<"\t\t\t" <<dir_name <<"\n";
				dir_A_2[l][k] = dir_A_1[k]->mkdir(dir_name);
				if(fun::nSparseIntegral(_acceptance_5d[k][l])>0.0){
					std::cout<<"\t\tPart 7e\n";
					dir_A_2[l][k]->cd();
					std::cout<<"\t\tPart 7\n";
					for(int m=0; m<5; m++){
						std::cout<<"\t\tPart 8\n";
						_accept_hist_1[k][l][m]->Write();
					}
				}
			}
		}
	}
	std::cout<<"Writing W Q2 projections of acceptance\n";
	dir_A->cd();
	_accept_hist_2[0]->SetXTitle("W (GeV)");
	_accept_hist_2[0]->SetYTitle("Integrated Acceptance");
	_accept_hist_2[0]->Write();
	_accept_hist_2[1]->SetXTitle("Q^{2} (GeV^{2})");
	_accept_hist_2[1]->SetYTitle("Integrated Acceptance");
	_accept_hist_2[1]->Write();
	//std::cout<<"\nPart 7a";
	//std::cout<<"\nPart 9";
	/*if(fun::nSparseIntegral(_acceptance_7d)>0.0){
		//std::cout<<"\nPart 10";
		dir_A_2->cd();
		
	}*/
	std::cout<<": Done\n";
}

void Histogram::Write_Single_Diff(Flags flags_){
	std::cout<<"Writing Single Differential Plots\n";
	if(!flags_.Flags::Plot_Single_Diff()){
		std::cout<<"\tNot Plotting Single Differentials\n";
	 	return;
	 	std::cout<<"\t\tProperly exited...?\n";
	}
	TDirectory* dir_sd = _RootOutputFile->mkdir("Single Diff Plots");
	dir_sd->cd();
	TDirectory* dir_sd_1[4];
	TDirectory* dir_sd_2[4][_n_bins_7d[1]];
	char dir_name[100];
	for(int k=0; k<4; k++){//Xij
		sprintf(dir_name,"Single_Diff_X:%s",_five_dim_[k]);
		dir_sd_1[k] = dir_sd->mkdir(dir_name);
		for(int l=0; l<_n_bins_7d[1]; l++){//Q2
			sprintf(dir_name,"Single_Diff_X:%s_Q2:%f-%f",_five_dim_[k],_bin_low_7d[1][l],_bin_up_7d[1][l]);
			dir_sd_2[k][l] = dir_sd_1[k]->mkdir(dir_name);
			dir_sd_2[k][l]->cd();
			for(int m=0; m<_n_bins_7d[0]; m++){
				_exp_corr_holes_3d[m][l][k]->Write();
			}
		}
	}
}

void Histogram::Write_Polarization(Flags flags_){
	std::cout<<"Writing Polarization Plots\n";
	if(!flags_.Flags::Plot_Pol()){
		std::cout<<"\tNot Plotting Polarization Cross Sections\n";
	 	return;
	 	std::cout<<"\t\tProperly exited...?\n";
	}
	TDirectory* dir_pol = _RootOutputFile->mkdir("Polarization Plots");
	dir_pol->cd();
	TDirectory* dir_pol_1[_n_bins_7d[0]];
	TDirectory* dir_pol_2[_n_bins_7d[0]][_n_bins_7d[1]];
	TDirectory* dir_pol_3[_n_bins_7d[0]][_n_bins_7d[1]][4];
	char dir_name[100];
	std::cout<<"\tWriting Polarization Dir and Hist\n";
	
	for(int k=0; k<_n_bins_7d[0]; k++){//W
		sprintf(dir_name,"Polarization_W:%f-%f",_bin_low_7d[0][k],_bin_up_7d[0][k]);
		dir_pol_1[k] = dir_pol->mkdir(dir_name);
		for(int l=0; l<_n_bins_7d[1]; l++){//Q2
			sprintf(dir_name,"Polarization_W:%f-%f_Q2:%f-%f",_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l]);
			dir_pol_2[k][l] = dir_pol_1[k]->mkdir(dir_name);
			for(int m=0; m<4; m++){//Xij
				sprintf(dir_name,"Polarization_X:%s_W:%f-%f_Q2:%f-%f",_five_dim_[m],_bin_low_7d[0][k],_bin_up_7d[0][k],_bin_low_7d[1][l],_bin_up_7d[1][l]);
				dir_pol_3[k][l][m] = dir_pol_2[k][l]->mkdir(dir_name);
				dir_pol_3[k][l][m]->cd();
				for(int n=0; n<_n_bins_5d[m]; n++){//Xij bins
					_exp_corr_holes_4d[k][l][m][n]->Write();
				}
			}
		}
	}
	std::cout<<"Done Writing Polarization Hists\n";
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
		_phi_bin_sizes->Fill(j,Histogram::Phi_Bin_Size(j));
	}
	for(int x=0; x<4; x++){
		//std::cout<<"X: " <<x <<"\n\tLowest X: " <<_bin_low_5d[x][0] <<" | highest x: " <<_bin_up_5d[x][_n_bins_5d[x]-1] <<"\n";
		sprintf(hname,"%s_bin_sizes",_five_dim_[x]);
		_X_bin_sizes.push_back(new TH1D(hname,hname,_n_bins_5d[x],_bin_low_5d[x][0],_bin_up_5d[x][_n_bins_5d[x]-1]));
		for(int i=0; i<_n_bins_5d[x]; i++){
			_X_bin_sizes[x]->Fill(i,Histogram::X_Bin_Size(x,i));
		}
	}
}


void Histogram::Acceptance_Errors(Flags flags_){
	std::cout<<"Acceptance Errors\n";
	int bin_7d[7];
	//_acceptance_eff_7d=(THnSparseD*)_acceptance_7d[var][top]->Clone();
	_acceptance_err_7d[0]=(THnSparseD*)_acceptance_7d->Clone();//Unweighted
	_acceptance_err_7d[1]=(THnSparseD*)_acceptance_7d->Clone();//Weighted
	_acceptance_eff_7d = (THnSparseD*)_acceptance_7d->Clone();
	//std::cout<<"\tAcceptance Cloned into error histograms\n";
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		bin_7d[0]=i;
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			bin_7d[1]=j;
			for(int k = 0; k< _n_bins_7d[2]; k++){ //MM1`
				bin_7d[2]=k;
				for(int l = 0; l< _n_bins_7d[3]; l++){ //MM2
					bin_7d[3]=l;
					for(int n = 0; n< _n_bins_7d[4]; n++){//theta
						bin_7d[4]=n;
						for(int m = 0; m< _n_bins_7d[5]; m++){ //alpha
							bin_7d[5]=m;
							for(int o = 0; o< _n_bins_7d[6]; o++){//Phi
								bin_7d[6]=o;
								if(_thrown_7d->THnSparse::GetBinContent(bin_7d) > 0.0 && _sim_data_7d->GetBinContent(bin_7d) > 0.0){
									//std::cout<<"\t\tFilling Errors and Efficiencies at: " <<i <<" " <<j <<" " <<k <<" " <<l <<" " <<m <<" " <<n <<" " <<o <<"\n";
									//"Unweighted" Error
									//std::cout<<"\t\t\tAcceptance at bin: " <<_acceptance_7d->THnSparse::GetBinContent(bin_7d) <<"\n";
									_acceptance_err_7d[0]->SetBinContent(bin_7d,physics::Error(_thrown_7d->THnSparse::GetBinContent(bin_7d),_sim_data_7d->GetBinContent(bin_7d)));
									//Weighted Error
									_acceptance_err_7d[1]->SetBinContent(bin_7d,physics::Error(_thrown_7d->THnSparse::GetBinContent(bin_7d),_sim_data_7d->GetBinContent(bin_7d),_sim_weight_sq_7d->GetBinContent(bin_7d)));
									//Efficiency
									_acceptance_eff_7d->SetBinContent(bin_7d,(_sim_data_7d->GetBinContent(bin_7d)/_thrown_7d->THnSparse::GetBinContent(bin_7d)));
									_acceptance_eff_7d->SetBinError(bin_7d,_acceptance_err_7d[1]->GetBinContent(bin_7d));//Give it the weighted error
								}
							}
						}
					}
				}
			}
		}
	}
}

void Histogram::Make_Error_Hists(Flags flags_){
	char hname[100];
	if(!flags_.Flags::Plot_Eff()){
		std::cout<<"\tNot Plotting Efficiency\n";
	 	return;
	}
	if(!flags_.Flags::Plot_Err()){
		std::cout<<"\tNot Plotting Errors\n";
	 	return;
	}

	TH1D_1d_star rel_err_w_1d;
	TH1D_1d_star rel_err_nw_1d;
	TH1D_1d_star yield_1d;
	TH1D_1d_star zero_thr_1d;
	
	sprintf(hname,"Fraction_Acc_Zero_wExp");
	_acc_zero_exp=new TH1D(hname,hname,101,-0.005,1.005);
	sprintf(hname,"Thrown_Exp_Ratio_%s_%s",_var_set_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()]);
	_thr_exp_ratio = new TH2D(hname,hname,_n_bins_7d[0],_acceptance_7d->GetAxis(0)->GetBinLowEdge(0),_acceptance_7d->GetAxis(0)->GetBinUpEdge(_n_bins_7d[0]-1),_n_bins_7d[1],_acceptance_7d->GetAxis(1)->GetBinLowEdge(0),_acceptance_7d->GetAxis(1)->GetBinUpEdge(_n_bins_7d[1]-1));
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			sprintf(hname,"Acc_Yield_%s_%s_W:%f-%f_Q:%f_%F",_var_set_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d->GetAxis(1)->GetBinUpEdge(j));
			yield_1d.push_back(new TH1D(hname,hname,1001,-0.001,0.2));
			sprintf(hname,"Fraction_Acc_Zero_wThr_%s_%s_W:%f-%f_Q:%f_%F",_var_set_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d->GetAxis(1)->GetBinUpEdge(j));
			zero_thr_1d.push_back(new TH1D(hname,hname,1001,-0.001,0.2));
			sprintf(hname,"Relative_Acc_Error_Weighted_%s_%s_W:%f-%f_Q:%f_%F",_var_set_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d->GetAxis(1)->GetBinUpEdge(j));
			rel_err_w_1d.push_back(new TH1D(hname,hname,1401,-0.005,1.405));
			sprintf(hname,"Relative_Acc_Error_Unweighted_%s_%s_W:%f-%f_Q:%f_%F",_var_set_[flags_.Flags::Var_idx()],_top_[flags_.Flags::Top_idx()],_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d->GetAxis(1)->GetBinUpEdge(j));
			rel_err_nw_1d.push_back(new TH1D(hname,hname,1401,-0.005,1.405));
			
		}
		if(yield_1d.size()>0){
			_acc_yield.push_back(yield_1d);
			yield_1d.clear();
		}
		if(zero_thr_1d.size()>0){
			_acc_zero_thr.push_back(zero_thr_1d);
			zero_thr_1d.clear();
		}
		if(rel_err_w_1d.size()>0){
			_acc_rel_error_weighted.push_back(rel_err_w_1d);
			rel_err_w_1d.clear();
		}
		if(rel_err_nw_1d.size()>0){
			_acc_rel_error_unweighted.push_back(rel_err_nw_1d);
			rel_err_nw_1d.clear();
		}
	}
}

void Histogram::Fill_Error_Hists(Flags flags_){
	std::cout<<"Filling Acceptance Error and Efficiency Histograms\n";
	if(!flags_.Flags::Plot_Eff()){
		std::cout<<"\tNot Plotting Efficiency\n";
	 	return;
	}
	if(!flags_.Flags::Plot_Err()){
		std::cout<<"\tNot Plotting Errors\n";
	 	return;
	}
	int bin_7d[7];
	int bin_5d[5];
	for(int i=0; i<_n_bins_7d[0]; i++){//W
		bin_7d[0]=i;
		for(int j=0; j< _n_bins_7d[1]; j++){//Q2
			bin_7d[1]=j;
			_thr_exp_ratio->Fill(Histogram::Bin_Center_7d(0,i),Histogram::Bin_Center_7d(0,i),fun::nSparseIntegral(_thrown_5d[i][j])/fun::nSparseIntegral(_exp_data_5d[i][j]));
			for(int k = 0; k< _n_bins_7d[2]; k++){ //MM1
				bin_7d[2]=k;
				bin_5d[0]=k;
				for(int l = 0; l< _n_bins_7d[3]; l++){ //MM2
					bin_7d[3]=l;
					bin_5d[1]=l;
					for(int n = 0; n< _n_bins_7d[4]; n++){//theta
						bin_7d[4]=n;
						bin_5d[2]=n;
						for(int m = 0; m< _n_bins_7d[5]; m++){ //alpha
							bin_7d[5]=m;
							bin_5d[3]=m;
							for(int o = 0; o< _n_bins_7d[6]; o++){//Phi
								bin_7d[6]=o;
								bin_5d[4]=o;
								if(_exp_data_7d->GetBinContent(bin_7d)>0.0){
									_acc_zero_exp->Fill(_acceptance_7d->GetBinContent(bin_7d));
									_acc_yield[i][j]->Fill(_acceptance_5d[i][j]->GetBinContent(bin_5d));
								}
								if(_thrown_7d->GetBinContent(bin_7d)>0.0){
									_acc_zero_thr[i][j]->Fill(_acceptance_7d->GetBinContent(bin_7d));
									_acc_rel_error_weighted[i][j]->Fill(_acceptance_err_7d[1]->GetBinContent(bin_7d)/_acceptance_7d->GetBinContent(bin_7d));
									_acc_rel_error_unweighted[i][j]->Fill(_acceptance_err_7d[0]->GetBinContent(bin_7d)/_acceptance_7d->GetBinContent(bin_7d));
								}
							}
						}
					}
				}
			}
		}
	}
	
}

void Histogram::Write_Error_Hists(Flags flags_){
	char dirname[100];
	std::cout<<"Writing Error Histogram\n";
	if(!flags_.Flags::Plot_Eff()){
		std::cout<<"\tNot Plotting Efficiency\n";
	 	return;
	 	std::cout<<"\t\tProperly exited...?\n";
	}
	if(!flags_.Flags::Plot_Err()){
		std::cout<<"\tNot Plotting Errors\n";
	 	return;
	 	std::cout<<"\t\tProperly exited...?\n";
	}
	TDirectory* dir_err = _RootOutputFile->mkdir("Acceptance Error Plots");
	TDirectory* dir_err_zero = dir_err->mkdir("Zero Acceptance");
	TDirectory* dir_err_zero_exp = dir_err_zero->mkdir("Zero Acceptance Exp Non-Zero");
	TDirectory* dir_err_zero_thr = dir_err_zero->mkdir("Zero Acceptance Thrown Non-Zero");
	TDirectory* dir_err_zero_thr_sub[_n_bins_7d[0]][_n_bins_7d[1]+1];
	TDirectory* dir_err_yield = dir_err->mkdir("Acceptance Yield Distribution");
	TDirectory* dir_err_yield_sub[_n_bins_7d[0]][_n_bins_7d[1]+1];
	TDirectory* dir_err_ratio = dir_err->mkdir("Thrown to Exp Ratio");
	TDirectory* dir_err_rel = dir_err->mkdir("Relative Error");
	TDirectory* dir_err_rel_sub[_n_bins_7d[0]][_n_bins_7d[1]+1];
	
	dir_err_ratio->cd();
	_thr_exp_ratio->SetXTitle("W (GeV)");
	_thr_exp_ratio->SetYTitle("Q2 (GeV^2)");
	_thr_exp_ratio->Write();
	dir_err_zero_exp->cd();
	_acc_zero_exp->SetXTitle("Fraction of Sim Performed");
	_acc_zero_exp->SetYTitle("Fraction of Bins of Zero");
	_acc_zero_exp->Write();
	for(int i=0; i<_n_bins_7d[0]; i++){
		sprintf(dirname,"Zero Acceptance Thr Non-Zero W:%f-%f",_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i));
		dir_err_zero_thr_sub[i][0] = dir_err_zero_thr->mkdir(dirname);
		sprintf(dirname,"Acceptance Yield Distribution W:%f-%f",_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i));
		dir_err_yield_sub[i][0] = dir_err_yield->mkdir(dirname);
		sprintf(dirname,"Relative Error W:%f-%f",_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i));
		dir_err_rel_sub[i][0] = dir_err_rel->mkdir(dirname);
		for(int j=0; j<_n_bins_7d[1]; j++){
			sprintf(dirname,"Zero Acceptance Thr Non-Zero W:%f-%f Q2:%f-%f",_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d->GetAxis(1)->GetBinUpEdge(j));
			dir_err_zero_thr_sub[i][j+1] = dir_err_zero_thr_sub[i][0]->mkdir(dirname);
			dir_err_zero_thr_sub[i][j+1]->cd();
			_acc_zero_thr[i][j]->SetXTitle("Fraction of Sim Performed");
			_acc_zero_thr[i][j]->SetYTitle("Fraction of Bins of Zero");
			_acc_zero_thr[i][j]->Write();
			sprintf(dirname,"Acceptance Yield Distribution W:%f-%f Q2:%f-%f",_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d->GetAxis(1)->GetBinUpEdge(j));
			dir_err_yield_sub[i][j+1] = dir_err_yield_sub[i][0]->mkdir(dirname);
			dir_err_yield_sub[i][j+1]->cd();
			_acc_yield[i][j]->SetXTitle("Acceptance Yield (5d)");
			_acc_yield[i][j]->SetYTitle("Number of Bins");
			_acc_yield[i][j]->Write();
			sprintf(dirname,"Relative Error W:%f-%f Q2:%f-%f",_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d->GetAxis(1)->GetBinUpEdge(j));
			dir_err_rel_sub[i][j+1] = dir_err_rel_sub[i][0]->mkdir(dirname);
			dir_err_rel_sub[i][j+1]->cd();
			_acc_rel_error_weighted[i][j]->SetXTitle("Relative Weighted Error");//{top,W,Q2}
			_acc_rel_error_weighted[i][j]->SetYTitle("Number of Bins");
			_acc_rel_error_weighted[i][j]->Write();
			_acc_rel_error_unweighted[i][j]->SetXTitle("Relative Unweighted Error");//{var,top,W,Q2}
			_acc_rel_error_unweighted[i][j]->SetYTitle("Number of Bins");
			_acc_rel_error_unweighted[i][j]->Write();
		}
	}
}

THnSparseD* Histogram::Add_Sparse(THnSparse * h1_, THnSparse * h2_){
	int dim1 = h1_->GetNdimensions();
	//int dim1 = h1->GetNdimensions();
	int dim2 = h2_->GetNdimensions();
	
	if(dim1 == dim2 && dim1==7){
		THnSparseD* new_h = (THnSparseD*)h1_->Clone();
		int bin_7d[7];
		//sprintf(hname,"exp_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[l][0][i],_bin_up_7d[l][0][i],_bin_low_7d[l][1][j],_bin_up_7d[l][1][j],topo[k],var_set[l]);
		//output = THnSparseD(hname,hname,dim1,fun::Vector_Array(_n_bins_7d),fun::Vector_Array(_bin_lowside_7d.at(var_set)),fun::Vector_Array(_bin_topside_7d.at(var_set)));
		for(int i=0; i<_n_bins_7d[0]; i++){//W
			bin_7d[0]=i;
			for(int j=0; j< _n_bins_7d[1]; j++){//Q2
				bin_7d[1]=j;
				for(int k = 0; k< _n_bins_7d[2]; k++){ //MM1
					bin_7d[2]=k;
					for(int l = 0; l< _n_bins_7d[3]; l++){ //MM2
						bin_7d[3]=l;
						for(int n = 0; n< _n_bins_7d[4]; n++){//theta
							bin_7d[4]=n;
							for(int m = 0; m< _n_bins_7d[5]; m++){ //alpha
								bin_7d[5]=m;
								for(int o = 0; o< _n_bins_7d[6]; o++){//Phi
									bin_7d[6]=o;
									//Only fill the sparse histogram where values are non-zero
									if(h2_->THnSparse::GetBinContent(bin_7d) > 0){//If the bin doesn't exist it will yield 0, but not create a bin
										new_h->THnSparse::AddBinContent(bin_7d,h2_->THnSparse::GetBinContent(bin_7d));
									}
								}
							}
						}
					}
				}
			}
		}
		return new_h;
	}else if(dim1 == dim2 && dim1==5){
		THnSparseD* new_h = (THnSparseD*)h1_->Clone();
		int bin_5d[5];
		//sprintf(hname,"exp_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[l][0][i],_bin_up_7d[l][0][i],_bin_low_7d[l][1][j],_bin_up_7d[l][1][j],topo[k],var_set[l]);
		//output = THnSparseD(hname,hname,dim1,fun::Vector_Array(_n_bins_7d),fun::Vector_Array(_bin_lowside_7d.at(var_set)),fun::Vector_Array(_bin_topside_7d.at(var_set)));
		for(int k = 0; k< _n_bins_7d[2]; k++){ //MM1
			bin_5d[0]=k;
			for(int l = 0; l< _n_bins_7d[3]; l++){ //MM2
				bin_5d[1]=l;
				for(int n = 0; n< _n_bins_7d[4]; n++){//theta
					bin_5d[2]=n;
					for(int m = 0; m< _n_bins_7d[5]; m++){ //alpha
						bin_5d[3]=m;
						for(int o = 0; o< _n_bins_7d[6]; o++){//Phi
							bin_5d[4]=o;
							//Only fill the sparse histogram where values are non-zero
							if(h2_->THnSparse::GetBinContent(bin_5d) > 0){//If the bin doesn't exist it will yield 0, but not create a bin
								new_h->THnSparse::AddBinContent(bin_5d,h2_->THnSparse::GetBinContent(bin_5d));
							}
						}
					}
				}
			}
		}
		return new_h;
	}
}	

