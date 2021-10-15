#include "histogram.hpp"

Histogram::Histogram(const std::string& output_file, TFile *exp_tree, TFile *sim_tree){
	//Make output file eventually RootOutputFile = fun::make_file
	std::cout<<"\nName_File";
	_RootOutputFile = fun::Name_File(output_file);
	Histogram::Extract_7d_Histograms(exp_tree,sim_tree);//Gets the 7d histograms off the root file
	Histogram::Extract_Bin_Info();//Extracts binning information about the 7d histograms just extracted
	Histogram::Skeleton_5D();//Creates empty vector arrays to map onto in the next step
	Histogram::Sparse_7to5();//Converts all the 7d histograms into usable 5d histograms for analysis
	Histogram::Sparse_5to3();
	Histogram::Sparse_5to4();
	Histogram::Convert_to_Cross();
	Histogram::Acceptance_Errors();
	Histogram::Make_Histograms();
	Histogram::Fill_Error_Hists();
	Histogram::Write_Histograms();
}

void Histogram::Make_Histograms(){
	std::cout<<"Making Histograms\n";
	Histogram::Make_WQ2();
	Histogram::Make_Acceptance();
	Histogram::Make_Error_Hists();
}

void Histogram::Fill_Histograms(){

}
	
void Histogram::Write_Histograms(){
	//Histogram::Write_WQ2();
	//Histogram::Write_Acceptance();
	Histogram::Write_Single_Diff();
	Histogram::Write_Polarization();
	Histogram::Write_Error_Hists();

}

void Histogram::Extract_7d_Histograms(TFile *exp_tree, TFile *sim_tree){
	std::cout<<"Extract 7d Histograms\n";
	char hname[100];
	//static const char * _sparse_names_[] = {"2#pi_off_proton_#Delta^{++}","2#pi_off_proton_#rho","2#pi_off_proton_#Delta^{0}"};
	//static const char * topo[] = {"Pmiss","PIPmiss","PIMmiss","Zeromiss","ALLmiss"};

	for(int i =0; i<_n_var_sets; i++){
		//std::cout<<"Getting Thrown "<<_sparse_names_[i] <<"\n";
		sprintf(hname,"Thrown_%s",_sparse_names_[i]);
		_thrown_7d[i] = (THnSparseD *)sim_tree->Get(hname);
		//std::cout<<"\tThrown integral: " <<fun::nSparseIntegral(_thrown_7d[i]) <<"\n";
		for(int j = 0; j<_n_topology; j++){
			//std::cout<<std::endl <<topo[j];
			//std::cout<<"Getting Exp "<<_sparse_names_[i] <<" " <<_topo_[j] <<"\n";
			sprintf(hname,"%s_%s",_sparse_names_[i],_topo_[j]);//_top_[j]);
			//std::cout<<"\nTHnSparse Name: " <<hname <<std::endl;
			_exp_data_7d[i][j] = (THnSparseD *)exp_tree->Get(hname);
			//std::cout<<"\tExp_data integral: " <<fun::nSparseIntegral(_exp_data_7d[i][j]) <<"\n";
			//std::cout<<"Getting Sim "<<_sparse_names_[i] <<" " <<_top_[j] <<"\n";
			sprintf(hname,"%s_%s",_sparse_names_[i],_top_[j]);
			_sim_data_7d[i][j] = (THnSparseD *)sim_tree->Get(hname);
			sprintf(hname,"Weight_%s_%s",_sparse_names_[i],_top_[j]);
			_sim_weight_sq_7d[i][j] = (THnSparseD *)sim_tree->Get(hname);
			//std::cout<<"\tSim_data integral: " <<fun::nSparseIntegral(_sim_data_7d[i][j]) <<"\n";
			//std::cout<<"Getting Acceptance "<<_sparse_names_[i] <<" " <<_top_[j] <<"\n";
			_acceptance_7d[i][j] = (THnSparseD*)_sim_data_7d[i][j]->Clone();
			_acceptance_7d[i][j]->Divide(_thrown_7d[i]);
			//std::cout<<"\tAcceptance integral: " <<fun::nSparseIntegral(_acceptance_7d[i][j]) <<"\n";
			//std::cout<<"Getting 7d Scale Factor "<<_sparse_names_[i] <<" " <<_top_[j] <<"\n";
			if(fun::nSparseIntegral(_sim_data_7d[i][j])>0.0){
				std::cout<<"\tSim Data 7d integral is non-zero\n";
				_scale_factor_7d[i][j]=fun::nSparseIntegral(_exp_data_7d[i][j])/fun::nSparseIntegral(_sim_data_7d[i][j]);
			}else{
				std::cout<<"\tSim Data 7d integral is zero\n";
				_scale_factor_7d[i][j]=0.0;
			}
			//std::cout<<"\nScale factor: " <<_scale_factor_7d[i][j];
			//std::cout<<"Getting Exp Corr "<<_sparse_names_[i] <<" " <<_top_[j] <<"\n";
			_exp_corr_7d[i][j]=(THnSparseD*)_exp_data_7d[i][j]->Clone();
			_exp_corr_7d[i][j]->Divide(_acceptance_7d[i][j]);
			//std::cout<<"\nexp_corr integral: " <<fun::nSparseIntegral(_exp_data_7d[i][j]);
			//std::cout<<"Getting Sim Corr "<<_sparse_names_[i] <<" " <<_top_[j] <<"\n";
			_sim_corr_7d[i][j]=(THnSparseD*)_sim_data_7d[i][j]->Clone();
			_sim_corr_7d[i][j]->Divide(_acceptance_7d[i][j]);
			//std::cout<<"\nsim_corr integral: " <<fun::nSparseIntegral(_exp_data_7d[i][j]);
			//Histogram::Sparse_Add_7d(_sim_holes_7d[i][j],_thrown_7d[i],sim_corr_7d[i][j],-1.0,i);
			//_exp_holes_7d[i][j]= (THnSparseD*)_sim_holes_7d[i][j]->Clone();
			//_exp_data_7d[i][j]->Scale(_scale_factor_7d[i][j]);
			//std::cout<<"Exp Hist Yield 7D: var:" <<var_set[i] <<" top:" <<topo[j] <<" => "	<<fun::nSparseIntegral(_exp_data_7d[i][j]) <<"| numbins:" <<_exp_data_7d[i][j]->GetNbins() <<"\n";
			//std::cout<<"Sim Hist Yield 7D: var:" <<var_set[i] <<" top:" <<topo[j] <<" => "	<<fun::nSparseIntegral(_sim_data_7d[i][j]) <<"| numbins:" <<_sim_data_7d[i][j]->GetNbins()<<"\n";	
		}
	}
}



void Histogram::Sparse_Add_7d(THnSparseD &h0,THnSparseD* h1, THnSparseD* h2, int sign, int var){
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
			h0.GetAxis(idx)->Set(_n_bins_7d[idx],fun::Vector_Array(_bin_edges_7d.at(var).at(idx)));
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

void Histogram::Sparse_Add_5d(THnSparseD* &h0, THnSparseD* h1, THnSparseD* h2, int sign, int var){
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
		for(int k = 0; k< _n_bins_5d[var][0]; k++){
			bin_5d[0]=k;
			for(int l = 0; l< _n_bins_5d[var][1]; l++){
				bin_5d[1]=l;
				for(int n = 0; n< _n_bins_5d[var][2]; n++){
					bin_5d[2]=n;
					for(int m = 0; m< _n_bins_5d[var][3]; m++){
						bin_5d[3]=m;
						for(int o = 0; o< _n_bins_5d[var][4]; o++){
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
void	Histogram::Sparse_7to5(){
	std::cout<<"Sparse_7to5\n";
	std::cout<<"\tFilling Those Sparse Friends\n";
	//Get the 7dimensional bins ready
	int bin_5d[5];
	int bin_7d[7];
	//std::cout<<"\nPart 1";
	long n_combo = _n_var_sets*_n_topology*_n_bins_7d[0]*_n_bins_7d[1]*_n_bins_7d[2]*_n_bins_7d[3]*_n_bins_7d[4]*_n_bins_7d[5]*_n_bins_7d[6];
	std::cout<<"Bins to iterate: " <<n_combo <<"\nPer:";

	long idx = 0;
	for(int var = 0; var<_n_var_sets; var++){		
		//std::cout<<" var set" <<var;
		for(int top=0; top<_n_topology; top++){
			for(int i=0; i<_n_bins_7d[0]; i++){//W
				bin_7d[0]=i;
				for(int j=0; j< _n_bins_7d[1]; j++){//Q2
					bin_7d[1]=j;
					for(int k = 0; k< _n_bins_7d[2]; k++){
						bin_5d[0]=k;
						bin_7d[2]=k;
						for(int l = 0; l< _n_bins_7d[3]; l++){
							bin_5d[1]=l;
							bin_7d[3]=l;
							for(int m = 0; m< _n_bins_7d[4]; m++){
								bin_5d[2]=m;
								bin_7d[4]=m;
								for(int n = 0; n< _n_bins_7d[5]; n++){
									bin_5d[3]=n;
									bin_7d[5]=n;
									for(int o = 0; o< _n_bins_7d[6]; o++){
										if(((o+1)*(n+1)*(m+1)*(l+1)*(k+1)*(j+1)*(i+1)*(top+1)*(var+1)-1)%(n_combo/100)==0){
											std::cout<<"\r" <<"\t" <<(100*((o+1)*(n+1)*(m+1)*(l+1)*(k+1)*(j+1)*(i+1)*(top+1)*(var+1)-1)/n_combo) <<" %" <<std::flush;
										}
										bin_5d[4]=o;
										bin_7d[6]=o;
										//Only fill the sparse histogram where values are non-zero
										//Exp 7d-5d
										//std::cout<<"\nMaking an Experimental Histogram: ";
										
										if(_exp_data_7d[var][top]->THnSparse::GetBinContent(bin_7d) > 0.0){//If the bin doesn't exist it will yield 0, but not create a bin
											_exp_data_5d[var][top][i][j]->THnSparse::AddBinContent(bin_5d,_exp_data_7d[var][top]->THnSparse::GetBinContent(bin_7d));
											//std::cout<<"Success at "<<bin_5d;
										}//Need this for scale factor_5d

										
										//Sim Recon 7d-5d
										//std::cout<<"\nMaking a Reconstructed Histogram: ";
										if(_sim_data_7d[var][top]->THnSparse::GetBinContent(bin_7d) > 0.0){//If the bin doesn't exist it will yield 0, but not create a bin
											_sim_data_5d[var][top][i][j]->THnSparse::AddBinContent(bin_5d,_sim_data_7d[var][top]->THnSparse::GetBinContent(bin_7d));
											//std::cout<<"Sucess at " <<bin_5d;
										}

										//Thrown 7d-5d
										//std::cout<<"\nMaking a Thrown Histogram: ";
										if(_thrown_7d[var]->THnSparse::GetBinContent(bin_7d) > 0.0 && top==0){//If the bin doesn't exist it will yield 0, but not create a bin
											_thrown_5d[var][i][j]->THnSparse::AddBinContent(bin_5d,_thrown_7d[var]->THnSparse::GetBinContent(bin_7d));
											//std::cout<<"Sucess at " <<bin_5d;
										}
										//Acceptance 7d-5d
										//std::cout<<"\nMaking an Acceptance Histogram: ";
										

										if(_acceptance_7d[var][top]->THnSparse::GetBinContent(bin_7d) > 0.0){
											_acceptance_5d[var][top][i][j]->THnSparse::AddBinContent(bin_5d,_acceptance_7d[var][top]->THnSparse::GetBinContent(bin_7d));
											//std::cout<<"Sucess at " <<bin_5d;
										}
										//Acceptance Corrected Yields
										if(_exp_corr_7d[var][top]->THnSparse::GetBinContent(bin_7d) > 0.0){
											_exp_corr_5d[var][top][i][j]->THnSparse::AddBinContent(bin_5d,_exp_corr_7d[var][top]->THnSparse::GetBinContent(bin_7d));
											//_exp_corr_holes_5d[var][top][i][j]->THnSparse::AddBinContent(bin_5d,_exp_corr_7d[var][top]->THnSparse::GetBinContent(bin_7d));
										}
										if(_sim_corr_7d[var][top]->THnSparse::GetBinContent(bin_7d) > 0.0){
											_sim_corr_5d[var][top][i][j]->THnSparse::AddBinContent(bin_5d,_sim_corr_7d[var][top]->THnSparse::GetBinContent(bin_7d));
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
					_scale_factor_5d[var][top][i][j]=fun::nSparseIntegral(_exp_data_5d[var][top][i][j])/fun::nSparseIntegral(_sim_data_5d[var][top][i][j]);
					//std::cout<<"Exp Hist 5D Yields| var:" <<var <<"top:" <<top <<"W:" <<i <<" Q2:" <<j <<" => " <<fun::nSparseIntegral(_exp_data_5d[top][var][i][j]) <<"\n";
					//std::cout<<"Sim Hist 5D Yields| var:" <<var <<"top:" <<top <<"W:" <<i <<" Q2:" <<j <<" => " <<fun::nSparseIntegral(_sim_data_5d[top][var][i][j]) <<"\n";
					//if(top==0){
						//std::cout<<"Thr Hist 5D Yields| var:" <<var <<"W:" <<i <<" Q2:" <<j <<" => " <<fun::nSparseIntegral(_thrown_5d[var][i][j]) <<"\n";
					//}
					//if(fun::nSparseIntegral(_acceptance[var][top][i][j])>0.0 && fun::nSparseIntegral(_thrown_5d[var][i][j])>0.0){
					//	_acceptance[var][top][i][j]->Divide(_thrown_5d[var][i][j]);
					//}
					if(fun::nSparseIntegral(_acceptance_5d[var][top][i][j]) > 0.0 && fun::nSparseIntegral(_sim_corr_5d[var][top][i][j]) >0.0){
						//std::cout<<"\nhave Acceptance and Sim data. ";
						//std::cout<<"\n\t|Adding Sim data to corr hist |";
						//_sim_corr_5d[var][top][i][j]= (THnSparseD*)_sim_data_5d[var][top][i][j]->Clone();
						//std::cout<<"\n\tDividing by acceptance |";
						//_sim_corr_5d[var][top][i][j]->Divide(_acceptance[var][top][i][j]);
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
											
											if(_thrown_5d[var][i][j]->THnSparse::GetBinContent(bin_5d)>0.0){
												/*std::cout<<"\n5d bin: {";
												for(int boor=0; boor<5; boor++){
													std::cout<<bin_5d[boor];
													if(boor!=4){
														std::cout<<", ";
													}else{
														std::cout<<"}";
													}
												}
												std::cout<<"\n\t thrown bin: " <<_thrown_5d[var][i][j]->THnSparse::GetBinContent(bin_5d);*/
												_sim_holes_5d[var][top][i][j]->THnSparse::AddBinContent(bin_5d,_thrown_5d[var][i][j]->THnSparse::GetBinContent(bin_5d));
											}
							//std::cout<<"\n\tSubtracting Sim Corr |";
											if(_sim_corr_5d[var][top][i][j]->THnSparse::GetBinContent(bin_5d) >0.0){
												/*std::cout<<"\n5d bin: {";
												for(int boor=0; boor<5; boor++){
													std::cout<<bin_5d[boor];
													if(boor!=4){
														std::cout<<", ";
													}else{
														std::cout<<"}";
													}
												}
												std::cout<<"\n\t -sim corr bin: " <<-1*(_sim_corr_5d[var][top][i][j]->GetBinContent(bin_5d));*/
												_sim_holes_5d[var][top][i][j]->THnSparse::AddBinContent(bin_5d,-1*(_sim_corr_5d[var][top][i][j]->THnSparse::GetBinContent(bin_5d)));
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
					if(fun::nSparseIntegral(_exp_data_5d[var][top][i][j]) >0.0 && fun::nSparseIntegral(_sim_data_5d[var][top][i][j]) >0.0){
						_scale_factor_5d[var][top][i][j]=fun::nSparseIntegral(_exp_data_5d[var][top][i][j])/fun::nSparseIntegral(_sim_data_5d[var][top][i][j]);
						//std::cout<<"Success at " <<bin_5d; 
					}else{
						_scale_factor_5d[var][top][i][j]=0.0;
						//std::cout<<"Zeroed at " <<bin_5d;
					}
					//Experimental correction 
					//std::cout<<"\nMaking Exp Corr Histogram: ";
					/*if(fun::nSparseIntegral(_acceptance[var][top][i][j]) > 0 && fun::nSparseIntegral(_sim_data_5d[var][top][i][j])>0.0){
						//_exp_corr_5d[var][top][i][j]->THnSparse::AddBinContent(bin_5d,_exp_data_5d[var][top][i][j]->THnSparse::GetBinContent(bin_5d));
						_exp_corr_5d[var][top][i][j]= (THnSparseD*)_exp_data_5d[var][top][i][j]->Clone();
						_exp_corr_5d[var][top][i][j]->Divide(_acceptance[var][top][i][j]);
						//std::cout<<"Sucess at " <<bin_5d;
					}*/
					//Experimental Hole estimate
					//std::cout<<"\nMaking Exp Holes: ";
					//std::cout<<"\nvar:" <<var <<" top:" <<top <<" W:" <<i <<" Q2:" <<j;
					//std::cout<<"\nSim_holes integral: " <<fun::nSparseIntegral(_sim_holes_5d[var][top][i][j]);
					//std::cout<<"\nScale Factor: " <<_scale_factor_5d[var][top][i][j];
					if(fun::nSparseIntegral(_sim_holes_5d[var][top][i][j])>0.0 && _scale_factor_5d[var][top][i][j]>0.0){
						//_exp_holes_5d[var][top][i][j]->THnSparse::AddBinContent(bin_5d,_sim_holes_5d[var][top][i][j]->THnSparse::GetBinContent(bin_5d));
						//std::cout<<"\nvar:" <<var <<" top:" <<top <<" W:" <<i <<" Q2:" <<j <<"| Int sim_holes:" <<fun::nSparseIntegral(_sim_holes_5d[var][top][i][j]);
						_exp_holes_5d[var][top][i][j]= (THnSparseD*)_sim_holes_5d[var][top][i][j]->Clone();
						_exp_holes_5d[var][top][i][j]->Scale(_scale_factor_5d[var][top][i][j]);
						_exp_corr_holes_5d[var][top][i][j] = (THnSparseD*)_exp_corr_5d[var][top][i][j]->Clone();
						//std::cout<<"\n\tExp Corr integral: " <<fun::nSparseIntegral(_exp_corr_5d[var][top][i][j]);
						//std::cout<<"\n\tExp holes integral: " <<fun::nSparseIntegral(_exp_holes_5d[var][top][i][j]);
						_exp_corr_holes_5d[var][top][i][j]->Add(_exp_holes_5d[var][top][i][j]);
						//std::cout<<"\n\tExp Corr +holes integral: " <<fun::nSparseIntegral(_exp_corr_holes_5d[var][top][i][j]);
						//std::cout<<"\n\tSuccess at " <<bin_5d;
					}
				}
			}
		}
	}
	std::cout<<"\n";
}
//For Single Differential bins
void Histogram::Sparse_5to3(){
	std::cout<<"Sparse 5 to 3\n";
	//Convert the _exp_corr_5d and _exp_holes_5d to 3 dimensional sparse histograms for usage in plotting single differential cross sections
	//std::cout<<"\nPart 1";
	char hname[100];
	char xlabel[100];
	char ylabel[100];

	//Get the 7dimensional bins ready
	TH1D_1d_star exp_ch_1d;
	TH1D_2d_star exp_ch_2d;
	TH1D_3d_star exp_ch_3d;
	TH1D_4d_star exp_ch_4d;
	TH1D * exp_corr;
	TH1D * exp_holes;
	THnSparseD * exp_corr_holes_5d;
	//std::cout<<"\nPart 2";
	for(int var = 0; var<_n_var_sets; var++){
		for(int top=0; top<_n_topology; top++){
			for(int i=0; i<_n_bins_7d[0]; i++){//W
				for(int j=0; j< _n_bins_7d[1]; j++){//Q2
					for(int k=0; k<4; k++){
						//std::cout<<"\nPart 3 " <<k;
						
						sprintf(hname,"%s_yield_acc_corr_holes_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_five_dim_[k],_bin_low_7d[var][0][i],_bin_up_7d[var][0][i],_bin_low_7d[var][1][j],_bin_up_7d[var][1][j],_top_[top],_var_set_[var]);
						sprintf(xlabel,"%s",_five_dim_[k]);
						sprintf(ylabel,"Yield");
						exp_ch_1d.push_back(_exp_corr_holes_5d[var][top][i][j]->Projection(k));
						//exp_ch_1d.push_back(_thrown_5d[var][i][j]->Projection(k));
						exp_ch_1d[k]->SetNameTitle(hname,hname);
						exp_ch_1d[k]->GetXaxis()->SetTitle(xlabel);
						exp_ch_1d[k]->GetYaxis()->SetTitle(ylabel);
					}
					//std::cout<<"\nPart 8 Q2:" <<j;
					exp_ch_2d.push_back(exp_ch_1d);
					exp_ch_1d.clear();
				}
				//std::cout<<"\nPart 9 W:" <<i;
				exp_ch_3d.push_back(exp_ch_2d);
				exp_ch_2d.clear();
			}
			//std::cout<<"\nPart 10 top:" <<top;
			exp_ch_4d.push_back(exp_ch_3d);
			exp_ch_3d.clear();
		}
		//std::cout<<"\nPart 11 var:" <<var;
		_exp_corr_holes_3d.push_back(exp_ch_4d);
		exp_ch_4d.clear();
	}
}

//For Polarization Observables
void Histogram::Sparse_5to4(){
	std::cout<<"Sparse 5 to 4\n";
	char hname[100];
	char xlabel[100];
	char ylabel[100];
	TH1D_1d_star exp_ch_1d;
	TH1D_2d_star exp_ch_2d;
	TH1D_3d_star exp_ch_3d;
	TH1D_4d_star exp_ch_4d;
	TH1D_5d_star exp_ch_5d;
	TH2D * exp_corr_2deg;
	TH2D * exp_holes_2deg;
	TH1D * exp_corr;
	TH1D * exp_holes;
	Int_t proj;
	Int_t proj_bins[2];
	THnSparseD * exp_corr_holes_5d;
	TH2D* exp_corr_holes_2d;
	//Convert the _exp_corr_5d and _exp_holes_5d to 4 dimensional sparse histograms for usage in plotting single differential cross sections
	for(int var = 0; var<_n_var_sets; var++){
		for(int top=0; top<_n_topology; top++){
			for(int i=0; i<_n_bins_7d[0]; i++){//W
				for(int j=0; j< _n_bins_7d[1]; j++){//Q2
					for(int k=0; k<4; k++){//Xij
						for(int l=0; l<_n_bins_5d[var][k]; l++){//Bins of Xij
							//exp_corr_holes_5d = (THnSparseD*)_exp_corr_holes_5d[var][top][i][j]->Clone()->;
							exp_corr_holes_2d = _exp_corr_holes_5d[var][top][i][j]->Projection(4,k);
							sprintf(hname,"%s+%f",_five_dim_[k],l);
							
							sprintf(xlabel,"Phi (deg)");
							sprintf(ylabel,"Yield");
							exp_corr_holes_2d->SetNameTitle(hname,hname);
							//exp_ch_1d.push_back(_exp_corr_holes_5d[var][top][i][j]->Projection(4,k)->ProjectionX(_five_dim_[k],l,l));
							exp_ch_1d.push_back(exp_corr_holes_2d->ProjectionX(_five_dim_[k],l,l+1));
							sprintf(hname,"%s:%f-%f_yield_acc_corr_holes_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_five_dim_[k],_bin_low_5d[var][k][l],_bin_up_5d[var][k][l],_bin_low_7d[var][0][i],_bin_up_7d[var][0][i],_bin_low_7d[var][1][j],_bin_up_7d[var][1][j],_top_[top],_var_set_[var]);
							exp_ch_1d[l]->SetNameTitle(hname,hname);
							exp_ch_1d[l]->GetXaxis()->SetTitle(xlabel);
							exp_ch_1d[l]->GetYaxis()->SetTitle(ylabel);
							//_exp_corr_4d[var][top][i][j][k][l]=Add(_exp_corr_4d[var][top][i][j][k][l],_exp_holes_4d[var][top][i][j][k][l]);
						}
						exp_ch_2d.push_back(exp_ch_1d);
						exp_ch_1d.clear();
					}
					exp_ch_3d.push_back(exp_ch_2d);
					exp_ch_2d.clear();
				}
				exp_ch_4d.push_back(exp_ch_3d);
				exp_ch_3d.clear();
			}
			exp_ch_5d.push_back(exp_ch_4d);
			exp_ch_4d.clear();
		}
		_exp_corr_holes_4d.push_back(exp_ch_5d);
		exp_ch_5d.clear();
	}
}

void Histogram::Extract_Bin_Info(){
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
	double_2d bin_low_7d_2d;
	double_2d bin_up_7d_2d;
	double_2d bin_mid_7d_2d;
	double_2d bin_size_7d_2d;
	double_2d bin_edges_7d_2d;
	double_1d bin_lowside_7d_1d;
	double_1d bin_topside_7d_1d;
	double_2d bin_low_5d_2d;
	double_2d bin_up_5d_2d;
	double_2d bin_mid_5d_2d;
	double_2d bin_size_5d_2d;
	double_2d bin_edges_5d_2d;
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
	for(int var = 0; var<_n_var_sets; var++){
		int top = 0;
		for(int bin_fun=0; bin_fun<5; bin_fun++){
			n_bins_5d.push_back(_exp_data_7d[var][0]->GetAxis(bin_fun+2)->GetNbins());;
		}
		_n_bins_5d.push_back(n_bins_5d);
		n_bins_5d.clear();
			DIM = _exp_data_7d[var][top]->GetNdimensions();
			for(int i = 0; i<DIM; i++){
				if(top ==0){//Just need to do this once in the loop 
					_n_bins_7d.push_back(_exp_data_7d[var][top]->GetAxis(i)->GetNbins());
				}
				bin_lowside_7d_1d.push_back(_exp_data_7d[var][top]->GetAxis(i)->GetBinLowEdge(1));
				bin_topside_7d_1d.push_back(_exp_data_7d[var][top]->GetAxis(i)->GetBinUpEdge(_n_bins_7d[i]));
				for(int j=0; j<_n_bins_7d[i];j++){
					bin_low_7d_1d.push_back(_exp_data_7d[var][top]->GetAxis(i)->GetBinLowEdge(j+1));
					bin_up_7d_1d.push_back(_exp_data_7d[var][top]->GetAxis(i)->GetBinUpEdge(j+1));
					bin_mid_7d_1d.push_back((bin_up_7d_1d[j]+bin_low_7d_1d[j])/2.0);
					bin_size_7d_1d.push_back(bin_up_7d_1d[j]-bin_low_7d_1d[j]);
					bin_edges_7d_1d.push_back(_exp_data_7d[var][top]->GetAxis(i)->GetBinLowEdge(j+1));
					if(j==_n_bins_7d[i]-1){
						bin_edges_7d_1d.push_back(_exp_data_7d[var][top]->GetAxis(i)->GetBinUpEdge(j+1));
					}
				}
				bin_low_7d_2d.push_back(bin_low_7d_1d);
				bin_low_7d_1d.clear();
				bin_up_7d_2d.push_back(bin_up_7d_1d);
				bin_up_7d_1d.clear();
				bin_mid_7d_2d.push_back(bin_mid_7d_1d);
				bin_mid_7d_1d.clear();
				bin_size_7d_2d.push_back(bin_size_7d_1d);
				bin_size_7d_1d.clear();
				bin_edges_7d_2d.push_back(bin_edges_7d_1d);
				bin_edges_7d_1d.clear();
			}
			_bin_low_7d.push_back(bin_low_7d_2d);
			bin_low_7d_2d.clear();
			_bin_up_7d.push_back(bin_up_7d_2d);
			bin_up_7d_2d.clear();
			_bin_mid_7d.push_back(bin_mid_7d_2d);
			bin_mid_7d_2d.clear();
			_bin_size_7d.push_back(bin_size_7d_2d);
			bin_size_7d_2d.clear();
			_bin_edges_7d.push_back(bin_edges_7d_2d);
			bin_edges_7d_2d.clear();
			_bin_lowside_7d.push_back(bin_lowside_7d_1d);
			bin_lowside_7d_1d.clear();
			_bin_topside_7d.push_back(bin_topside_7d_1d);
			bin_topside_7d_1d.clear();
			//Get the 5d bins ready
			for(int j=0; j<DIM-2; j++){
				//if(var==0 && top==0){
				//	_n_bins_5d.push_back(_exp_data_7d[var][top]->GetAxis(j+2)->GetNbins());
					//std::cout<<"\nVar set: " <<var <<" num_bins: " <<_exp_data_7d[var][top]->GetAxis(j+2)->GetNbins();
				//}
				bin_lowside_5d_1d.push_back(_exp_data_7d[var][top]->GetAxis(j+2)->GetBinLowEdge(1));
				bin_topside_5d_1d.push_back(_exp_data_7d[var][top]->GetAxis(j+2)->GetBinUpEdge(_n_bins_5d[var][j]));
				for(int k=0; k<_n_bins_5d[var][j];k++){
					bin_low_5d_1d.push_back(_exp_data_7d[var][top]->GetAxis(j+2)->GetBinLowEdge(k+1));
					bin_up_5d_1d.push_back(_exp_data_7d[var][top]->GetAxis(j+2)->GetBinUpEdge(k+1));
					bin_mid_5d_1d.push_back((bin_up_5d_1d[k]+bin_low_5d_1d[k])/2.0);
					bin_size_5d_1d.push_back(bin_up_5d_1d[k]-bin_low_5d_1d[k]);
					bin_edges_5d_1d.push_back(_exp_data_7d[var][top]->GetAxis(j+2)->GetBinUpEdge(k+1));
					if(k==_n_bins_5d[var][j]-1){
						bin_edges_5d_1d.push_back(_exp_data_7d[var][top]->GetAxis(j+2)->GetBinUpEdge(k+1));
					}
				}
				bin_low_5d_2d.push_back(bin_low_5d_1d);
				bin_low_5d_1d.clear();
				bin_up_5d_2d.push_back(bin_up_5d_1d);
				bin_up_5d_1d.clear();
				bin_mid_5d_2d.push_back(bin_mid_5d_1d);
				bin_mid_5d_1d.clear();
				bin_size_5d_2d.push_back(bin_size_5d_1d);
				bin_size_5d_1d.clear();
				bin_edges_5d_2d.push_back(bin_edges_5d_1d);
				bin_edges_5d_1d.clear();
			}
			_bin_low_5d.push_back(bin_low_5d_2d);
			bin_low_5d_2d.clear();
			_bin_up_5d.push_back(bin_up_5d_2d);
			bin_up_5d_2d.clear();
			_bin_mid_5d.push_back(bin_mid_5d_2d);
			bin_mid_5d_2d.clear();
			_bin_size_5d.push_back(bin_size_5d_2d);
			bin_size_5d_2d.clear();
			_bin_edges_5d.push_back(bin_edges_5d_2d);
			bin_edges_5d_2d.clear();
			_bin_lowside_5d.push_back(bin_lowside_5d_1d);
			bin_lowside_5d_1d.clear();
			_bin_topside_5d.push_back(bin_topside_5d_1d);
			bin_topside_5d_1d.clear();
	}

}

//Make the skeletons of the histograms to subsequently fill
void Histogram::Skeleton_5D(){
	std::cout<<"Skeleton_5D\n";
	char hname[100];
	int DIM = _exp_data_7d[0][0]->GetNdimensions();
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

	double_3d scale_3d;
	double_2d scale_2d;
	double_1d scale_1d;

	/*Sparse_3d exp_set;
	Sparse_2d exp_set_1;
	Sparse_1d exp_set_2;
	Sparse_2d thrown_set;
	Sparse_1d thrown_set_1;*/
	for(int i = 0; i<_n_var_sets; i++){//var sets
		for(int j=0; j<_n_topology; j++){//topology
			for(int k=0; k<_n_bins_7d[0]; k++){//W	
				for(int l=0; l< _n_bins_7d[1]; l++){//Q2
					for(int boop=0; boop<5; boop++){//Loop over all bins
						bin_sizes_5d[boop] = _n_bins_5d[i][boop];//Get sizes for individual bins
						bin_low_5d[boop] = _bin_low_5d[i][boop][0];//Get low end for each bin
						bin_up_5d[boop] = _bin_up_5d[i][boop][_bin_up_5d[i][boop].size()-1];//Get high end for each bin
					}
					//exp
					sprintf(hname,"exp_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_top_[j],_var_set_[i]);
					//_exp_data_5d[i][j][k].push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
					exp_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
					//sim
					sprintf(hname,"sim_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_top_[j],_var_set_[i]);
					//_sim_data_5d[i][j][k].push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
					sim_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
					
					
					sprintf(hname,"sim_holes_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_top_[j],_var_set_[i]);
					//_sim_holes_5d[i][j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins[i],bin_low_5d,bin_up_5d));
					sim_hole_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));

					sprintf(hname,"exp_holes_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_top_[j],_var_set_[i]);
					//std::cout<<"\nexpHoles hname: " <<hname;
					//_exp_holes_5d[i][j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins[i],bin_low_5d,bin_up_5d));
					exp_hole_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
					sprintf(hname,"exp_corr_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_top_[j],_var_set_[i]);
					//_exp_corr_5d[i][j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins[i],bin_low_5d,bin_up_5d));
					exp_corr_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));

					sprintf(hname,"exp_corr_hole_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_top_[j],_var_set_[i]);
					//_exp_corr_5d[i][j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins[i],bin_low_5d,bin_up_5d));
					exp_corr_hole_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));

					sprintf(hname,"sim_corr_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_top_[j],_var_set_[i]);
					//_sim_corr_5d[i][j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins[i],bin_low_5d,bin_up_5d));
					sim_corr_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));

					sprintf(hname,"xs_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_top_[j],_var_set_[i]);
					//_cross_section_5d[i][j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins[i],bin_low_5d,bin_up_5d));
					cross_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));

					sprintf(hname,"acceptance_5d_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_top_[j],_var_set_[i]);
					//_acceptance[i][j][k].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins[i],bin_low_5d,bin_up_5d));
					accept_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
					
					if(j ==0){
						sprintf(hname,"thrown_5d_W:%f-%f_Q2:%f-%f_var:%s",_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_var_set_[i]);
						//_thrown_5d[l][i].push_back(new THnSparseD(hname,hname,DIM-2,_n_bins[i],bin_low_5d,bin_up_5d));
						thr_1d.push_back(new THnSparseD(hname,hname,DIM-2,bin_sizes_5d,bin_low_5d,bin_up_5d));
					}
					scale_1d.push_back(NAN);	
					//curr_hist=curr_hist+1;
					//percent=100.0*(curr_hist/num_hist);
					//std::cout<<"Skeleton_5D " <<percent <<"\r";
					//std::cout<<"\r" <<"\t" <<"Skeleton_5D " <<(100*curr_hist/num_hist) <<" %"  <<std::flush;
				}//Q2
				
				exp_2d.push_back(exp_1d);
				sim_2d.push_back(sim_1d);
				sim_hole_2d.push_back(sim_hole_1d);
				exp_hole_2d.push_back(exp_hole_1d);
				exp_corr_hole_2d.push_back(exp_corr_hole_1d);
				sim_corr_2d.push_back(sim_corr_1d);
				exp_corr_2d.push_back(exp_corr_1d);
				cross_2d.push_back(cross_1d);
				accept_2d.push_back(accept_1d);
				scale_2d.push_back(scale_1d);
				if(j==0){
					thr_2d.push_back(thr_1d);
					thr_1d.clear();
				}
				
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
			exp_3d.push_back(exp_2d);
			sim_3d.push_back(sim_2d);
			sim_hole_3d.push_back(sim_hole_2d);
			exp_hole_3d.push_back(exp_hole_2d);
			sim_corr_3d.push_back(sim_corr_2d);
			exp_corr_3d.push_back(exp_corr_2d);
			exp_corr_hole_3d.push_back(exp_corr_hole_2d);
			cross_3d.push_back(cross_2d);
			accept_3d.push_back(accept_2d);
			scale_3d.push_back(scale_2d);
			exp_2d.clear();
			sim_2d.clear();
			sim_hole_2d.clear();
			exp_hole_2d.clear();
			exp_corr_hole_2d.clear();
			exp_corr_2d.clear();
			sim_corr_2d.clear();
			cross_2d.clear();
			accept_2d.clear();
			scale_2d.clear();
			//_thrown_5d.push_back(thr_2d););
		}//Top loop
		_exp_data_5d.push_back(exp_3d);
		_sim_data_5d.push_back(sim_3d);
		_sim_holes_5d.push_back(sim_hole_3d);
		_exp_holes_5d.push_back(exp_hole_3d);
		_exp_corr_holes_5d.push_back(exp_corr_hole_3d);
		_exp_corr_5d.push_back(exp_corr_3d);
		_sim_corr_5d.push_back(sim_corr_3d);
		_cross_section_5d.push_back(cross_3d);
		_acceptance_5d.push_back(accept_3d);
		_thrown_5d.push_back(thr_2d);
		_scale_factor_5d.push_back(scale_3d);
		exp_3d.clear();
		sim_3d.clear();
		sim_hole_3d.clear();
		exp_hole_3d.clear();
		exp_corr_hole_3d.clear();
		exp_corr_3d.clear();
		sim_corr_3d.clear();
		cross_3d.clear();
		accept_3d.clear();
		thr_2d.clear();
		scale_3d.clear();
		/*_exp_data_5d.push_back(exp_set_2);
		_sim_data_5d.push_back(sim_set_2);
		exp_set_2.std::vector::erase(exp_set_2.begin(),exp_set_2.begin()+exp_set_2.size());
		sim_set_2.std::vector::erase(sim_set_2.begin(),sim_set_2.begin()+sim_set_2.size());*/
	}//var set loop
	//std::cout<<"\n";
	std::cout<<"End Skeleton_5D\n";
	//std::cout<<"\nCheck: " <<fun::nSparseIntegral(_exp_data_5d[0][0][0][0]) <<"\n";
}

void	Histogram::Calc_Acceptance(){
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
	for(int i =0; i<_n_var_sets; i++){
		//std::cout<<"Part 2 idx: "<<i <<"\n";
		/*_acceptance.push_back(nsparse_set_3d);
		_exp_corr_5d.push_back(nsparse_set_3d);
		_sim_corr_5d.push_back(nsparse_set_3d);
		_n_exp_corr.push_back(double_set_3d);
		_n_sim_corr.push_back(double_set_3d);
		_scale_factor_5d.push_back(double_set_3d);
		_n_thrown.push_back(double_set_2d);*/
		for(int j = 0; j<_n_topology; j++){
			//std::cout<<"Part 3 idx: "<<i <<" " <<j <<"\n";
			/*_acceptance[i].push_back(nsparse_set_2d);
			_exp_corr_5d[i].push_back(nsparse_set_2d);
			_sim_corr_5d[i].push_back(nsparse_set_2d);
			_n_exp_corr[i].push_back(double_set_2d);
			_n_sim_corr[i].push_back(double_set_2d);
			_scale_factor_5d[i].push_back(double_set_2d);*/
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
					if(fun::nSparseIntegral(_sim_data_5d[i][j][k][l]) > 0.0 && fun::nSparseIntegral(_thrown_5d[i][k][l]) > 0.0 ){//.at(i).at(j).at(k).at(l))
						std::cout<<"Sim Yield: " <<fun::nSparseIntegral(_sim_data_5d[i][j][k][l]) <<"\nThrown Yield: " <<fun::nSparseIntegral(_thrown_5d[i][k][l]) <<"\n";
						//std::cout<<"*-*idx_l = " <<l <<" | idx_size = " <<_acceptance[i][j][k].size() <<"\n";
						//_acceptance[i][j][k].push_back(_sim_data_5d[i][j][k][l]);
						accept_1d.push_back(_sim_data_5d[i][j][k][l]);
						//std::cout<<"*** Proper Fill *** idx_l = " <<l <<" | idx_size = " <<_acceptance[i][j][k].size() <<"\n";
						std::cout<<"*** Proper Fill *** idx_l = " <<l <<" | idx_size = " <<accept_1d.size() <<"\n";
						//std::cout<<"Acceptance Pre-Div Yield: " <<fun::nSparseIntegral(_acceptance[i][j][k][l]) <<"\n";
						//std::cout<<"Acceptance Pre-Div Yield: " <<fun::nSparseIntegral(accept_1d[l]) <<"\n";
						//_acceptance[i][j][k][l]->Divide(_thrown_5d[i][k][l]);
						accept_1d[l]->Divide(_thrown_5d[i][k][l]);
						//std::cout<<"Acceptance Divided Yield: " <<fun::nSparseIntegral(_acceptance[i][j][k][l]) <<"\n";
						std::cout<<"Acceptance Divided Yield: " <<fun::nSparseIntegral(accept_1d[l]) <<"\n";
					}else{
					
						//_acceptance[i][j][k].push_back(_sim_data_5d[i][j][k][l]);
						accept_1d.push_back(_sim_data_5d[i][j][k][l]);
						//std::cout<<"*** Empty Fill *** idx_l = " <<l <<" | idx_size = " <<_acceptance[i][j][k].size() <<"\n";
						std::cout<<"*** Empty Fill *** idx_l = " <<l <<" | idx_size = " <<accept_1d.size() <<"\n";
					}

					//std::cout<<"acceptance corr calc\n";
					//if(fun::nSparseIntegral(_acceptance[i][j][k][l]) >0.0){
					if(fun::nSparseIntegral(accept_1d[l]) >0.0){
						//std::cout<<"exp corr calc\n";
						if(fun::nSparseIntegral(_exp_data_5d[i][j][k][l]) > 0.0 ){
							//std::cout<<"Exp Yield: " <<fun::nSparseIntegral(_exp_data_5d[i][j][k][l]) <<"\n";
							//_exp_corr_5d[i][j][k].push_back(_exp_data_5d[i][j][k][l]);
							exp_corr_1d.push_back(_exp_data_5d[i][j][k][l]);
							//std::cout<<"*** Empty EXP Corr Fill *** idx_l = " <<l <<" | idx_size = " <<_exp_corr_5d[i][j][k].size() <<"\n";
							//std::cout<<"*** Proper EXP Corr Fill *** idx_l = " <<l <<" | idx_size = " <<exp_corr_1d.size() <<"\n";
							//_exp_corr_5d[i][j][k][l]->Divide(_acceptance[i][j][k][l]);
							//_n_exp_corr[i][j][k].push_back(fun::nSparseIntegral(_exp_corr_5d[i][j][k][l]));
							exp_corr_1d[l]->Divide(accept_1d[l]);
							n_exp_1d.push_back(fun::nSparseIntegral(exp_corr_1d[l]));
						}else{
							//_n_exp_corr[i][j][k].push_back(0.0);
							//_exp_corr_5d[i][j][k].push_back(_exp_data_5d[i][j][k][l]);
							n_exp_1d.push_back(0.0);
							exp_corr_1d.push_back(_exp_data_5d[i][j][k][l]);
						}
						//std::cout<<"sim corr calc\n";
						if(fun::nSparseIntegral(_sim_data_5d[i][j][k][l])>0.0){
							//_sim_corr_5d[i][j][k].push_back(_sim_data_5d[i][j][k][l]);
							//_sim_corr_5d[i][j][k][l]->Divide(_acceptance[i][j][k][l]);
							//_n_exp_corr[i][j][k].push_back(fun::nSparseIntegral(_sim_corr_5d[i][j][k][l]));
							sim_corr_1d.push_back(_sim_data_5d[i][j][k][l]);
							sim_corr_1d[l]->Divide(accept_1d[l]);
							n_sim_1d.push_back(fun::nSparseIntegral(sim_corr_1d[l]));
						}else{
							//_n_sim_corr[i][j][k].push_back(0.0);
							//_sim_corr_5d[i][j][k].push_back(_sim_data_5d[i][j][k][l]);
							n_sim_1d.push_back(0.0);
							sim_corr_1d.push_back(_sim_data_5d[i][j][k][l]);
						}
					}else{
						//_n_exp_corr[i][j][k].push_back(0.0);
						//_n_sim_corr[i][j][k].push_back(0.0);
						//_exp_corr_5d[i][j][k].push_back(_exp_data_5d[i][j][k][l]);
						//_sim_corr_5d[i][j][k].push_back(_sim_data_5d[i][j][k][l]);
						n_exp_1d.push_back(0.0);
						n_sim_1d.push_back(0.0);
						exp_corr_1d.push_back(_exp_data_5d[i][j][k][l]);
						sim_corr_1d.push_back(_sim_data_5d[i][j][k][l]);
					}
					//std::cout<<"Corr Yields\n" <<"exp: " <<_n_exp_corr[i][j][k][l] <<"| sim: " <<_n_sim_corr[i][j][k][l];
					std::cout<<"Corr Yields\n" <<"exp: " <<n_exp_1d[l] <<"| sim: " <<n_sim_1d[l];
					//if(_n_sim_corr[i][j][k][l]>0.0){
					if(n_sim_1d[l]>0.0){
						//_scale_factor_5d[i][j][k].push_back(_n_exp_corr[i][j][k][l]/_n_sim_corr[i][j][k][l]);
						scale_1d.push_back(n_exp_1d[l]/n_sim_1d[l]);

					}else{
						//_scale_factor_5d[i][j][k].push_back(0.0);
						scale_1d.push_back(0.0);
					}
					//std::cout<<" scale_factor: " <<_scale_factor_5d[i][j][k][l] <<"\n";
					std::cout<<" scale_factor: " <<scale_1d[l] <<"\n";
					/*accept_set.push_back(THnSparse::Divide(_sim_data_5d[i][j][k][l],_thrown_5d[i][k][l]));
					accept_set.push_back(_sim_data_5d[i][j][k][l]);
					accept_set.at(l).THnSparse::Divide(_thrown_5d[i][k][l]);
					exp_corr_set.push_back(THnSparse::Divide(_exp_data_5d[i][j][k][l],accept_set[l]));
					sim_corr_set.push_back(THnSparse::Divide(_sim_data_5d[i][j][k][l],accept_set[l]));
					n_exp_set.push_back(exp_corr_set->THnSparse::ComputeIntegral());
					n_sim_set.push_back(sim_corr_set->THnSparse::ComputeIntegral());
					scale_set.push_back(n_exp_set[l]/n_sim_set[l]);*/
					if(j==0){
						if(fun::nSparseIntegral(_thrown_5d[i][k][l]) > 0.0){
							//_n_thrown[i][k].push_back(fun::nSparseIntegral(_thrown_5d[i][k][l]));
							n_thr_1d.push_back(fun::nSparseIntegral(_thrown_5d[i][k][l]));
						}else{
							//_n_thrown[i][k].push_back(0.0);
							n_thr_1d.push_back(0.0);
						}
						
					}
					
				}//l loop
				accept_2d.push_back(accept_1d);
				exp_corr_2d.push_back(exp_corr_1d);
				sim_corr_2d.push_back(sim_corr_1d);
				n_exp_2d.push_back(n_exp_1d);
				n_sim_2d.push_back(n_sim_1d);
				scale_2d.push_back(scale_1d);
				if(j==0){
					n_thr_2d.push_back(n_thr_1d);
					n_thr_1d.clear();
				}
				accept_1d.clear();
				exp_corr_1d.clear();
				sim_corr_1d.clear();
				n_exp_1d.clear();
				n_sim_1d.clear();
				scale_1d.clear();
				
			}//k loop
			accept_3d.push_back(accept_2d);
			exp_corr_3d.push_back(exp_corr_2d);
			sim_corr_3d.push_back(sim_corr_2d);
			n_exp_3d.push_back(n_exp_2d);
			n_sim_3d.push_back(n_sim_2d);
			scale_3d.push_back(scale_2d);
			accept_2d.clear();
			exp_corr_2d.clear();
			sim_corr_2d.clear();
			n_exp_2d.clear();
			n_sim_2d.clear();
			scale_2d.clear();
		}//j loop
		_acceptance_5d.push_back(accept_3d);
		_exp_corr_5d.push_back(exp_corr_3d);
		_sim_corr_5d.push_back(sim_corr_3d);
		_n_exp_corr.push_back(n_exp_3d);
		_n_sim_corr.push_back(n_sim_3d);
		_scale_factor_5d.push_back(scale_3d);
		_n_thrown.push_back(n_thr_2d);
		n_thr_2d.clear();
		accept_3d.clear();
		exp_corr_3d.clear();
		sim_corr_3d.clear();
		n_exp_3d.clear();
		n_sim_3d.clear();
		scale_3d.clear();
	}//i loop
	for(int i=0; i<_n_var_sets; i++){
		for(int j=0; j<_n_topology; j++){
			for(int k=0; k<_n_bins_7d[0]; k++){
				for(int l=0; l<_n_bins_7d[1]; l++){
					std::cout<<"\n Acceptance in Bin: " <<i <<" " <<j <<" " <<k <<" " <<l <<" _| " <<fun::nSparseIntegral(_acceptance_5d[i][j][k][l]) <<"\n";
				}
			}
		}
	}
}

void	Histogram::Calc_Holes_Sim(){
	std::cout<<"Calc Holes Sim\n";
	std::vector<THnSparseD> sim_hole_set;
	for(int i =0; i<_n_var_sets; i++){
		for(int j = 0; j<_n_topology; j++){
			for(int k=0; k<_n_bins_7d[0]; k++){
				for(int l=0; l<_n_bins_7d[1]; l++){
					//sim_hole_set.push_back(_thrown_5d.at(i).at(k).at(l)-_sim_corr_5d.at(i).at(j).at(k).at(l));
					Sparse_Add_5d(_sim_holes_5d[i][j][k][l], _thrown_5d[i][k][l], _sim_corr_5d[i][j][k][l], -1, i);
				}
				//_sim_holes_5d[i][j][k].push_back(sim_hole_set);
				//sim_hole_set.std::vector::erase(sim_hole_set.begin(),sim_hole_set.size());
			}
		}
	}
}

void	Histogram::Calc_Holes_Exp(){
	std::cout<<"Calc Holes Exp\n";
	Sparse_3d_star exp_hole_3d;
	Sparse_2d_star exp_hole_2d;
	Sparse_1d_star exp_hole_1d;
	for(int i =0; i<_n_var_sets; i++){
		for(int j = 0; j<_n_topology; j++){
			for(int k=0; k<_n_bins_7d[0]; k++){
				for(int l=0; l<_n_bins_7d[1]; l++){
					exp_hole_1d.push_back(_sim_holes_5d[i][j][k][l]);
					//_exp_holes_5d[i][j][k][l]=_sim_holes_5d[i][j][k][l];
					//exp_hole_set.push_back(_sim_holes_5d[i][j][k][l]);
					//exp_hole_set[l].Scale(_scale_factor_5d[i][j][k][l]);
					exp_hole_1d[l]->Scale(_scale_factor_5d[i][j][k][l]);
					//_exp_holes_5d[i][j][k][l]->Scale(_scale_factor_5d[i][j][k][l]);
				}
				//_exp_holes_5d.at(i).at(j).at(k).push_back(exp_hole_set);
				//exp_hole_set.std::vector::erase(exp_hole_set.begin(),exp_hole_set.size());
				exp_hole_2d.push_back(exp_hole_1d);
				exp_hole_1d.clear();
			}
			exp_hole_3d.push_back(exp_hole_2d);
			exp_hole_2d.clear();
		}
		_exp_holes_5d.push_back(exp_hole_3d);
		exp_hole_3d.clear();
	}
}

void	Histogram::Calc_Holes(){
	Histogram::Calc_Holes_Sim();
	Histogram::Calc_Holes_Exp();
}

void Histogram::Calc_Cross_Section(){
	std::cout<<"Calc";

}

void	Histogram::Make_Single_Diff(){
	std::cout<<"Calc";
}

void	Histogram::Make_Polarization(){
	std::cout<<"Calc";
}

void	Histogram::Make_Integrated(){
	std::cout<<"Calc";
}

void Histogram::Make_WQ2(){
	std::cout<<"\nMaking WQ2 Plots";
	char hname[100];
	for(int i =0; i<_n_var_sets; i++){
		for(int j=0; j<_n_topology; j++){
			sprintf(hname,"Exp_WQ2_Top:%s_Var:%s",_top_[j],_var_set_[i]);
			_exp_hist_wq2[i][j]= _exp_data_7d[i][j]->Projection(1,0);
			_exp_hist_wq2[i][j]->SetNameTitle(hname,hname);
			sprintf(hname,"Sim_WQ2_Top:%s_Var:%s",_top_[j],_var_set_[i]);
			_sim_hist_wq2[i][j]= _sim_data_7d[i][j]->Projection(1,0);
			_sim_hist_wq2[i][j]->SetNameTitle(hname,hname);
			sprintf(hname,"Exp_Corr_WQ2_Top:%s_Var:%s",_top_[j],_var_set_[i]);
			_exp_corr_hist_wq2[i][j]= _exp_corr_7d[i][j]->Projection(1,0);
			_exp_corr_hist_wq2[i][j]->SetNameTitle(hname,hname);
			sprintf(hname,"Sim_Corr_WQ2_Top:%s_Var:%s",_top_[j],_var_set_[i]);
			_sim_corr_hist_wq2[i][j]= _sim_corr_7d[i][j]->Projection(1,0);
			_sim_corr_hist_wq2[i][j]->SetNameTitle(hname,hname);
		}
		sprintf(hname,"Thr_WQ2_Var:%s",_var_set_[i]);
		_thr_hist_wq2[i]= _thrown_7d[i]->Projection(1,0);
		_thr_hist_wq2[i]->SetNameTitle(hname,hname);
	}
}

void Histogram::Make_Acceptance(){
	std::cout<<"\nMaking Acceptance Plots";
	//TH1D_5d _accept_hist_1;//(var,top,w,q2,{MM1,MM2,theta,alpha,phi})
	//TH1D_3d _accept_hist_2;//(var,top,{w,q2})
	//std::cout<<"\nMaking WQ2 Plots";
	char hname[100];
	TH1D_1d_star accept_hist_01;
	TH1D_2d_star accept_hist_02;
	TH1D_3d_star accept_hist_03;
	TH1D_4d_star accept_hist_04;
	for(int i =0; i<_n_var_sets; i++){
		//_accept_hist_1.push_back(accept_hist_01);
		for(int j=0; j<_n_topology; j++){
			//_accept_hist_1[i].push_back(accept_hist_02);
			for(int k=0; k<_n_bins_7d[0]; k++){
				//_accept_hist_1[i][j].push_back(accept_hist_03);
				for(int l=0; l<_n_bins_7d[1]; l++){
					//_accept_hist_1[i][j][k].push_back(accept_hist_04);
					for(int m=0; m<5; m++){
						//std::cout<<"\nstep 1";
						sprintf(hname,"acceptance_%s_W:%f-%f_Q2:%f-%f_top:%s_var:%s",_five_dim_[m],_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l],_top_[j],_var_set_[i]);
						accept_hist_01.push_back(_acceptance_5d[i][j][k][l]->Projection(m));
					}
					//std::cout<<"\nstep 2";
					accept_hist_02.push_back(accept_hist_01);
					accept_hist_01.clear();
				}
				//std::cout<<"\nstep 3";
				accept_hist_03.push_back(accept_hist_02);
				accept_hist_02.clear();
			}
			_accept_hist_2[i][j][0]=(TH1D*)_acceptance_7d[i][j]->Projection(0)->Clone();
			_accept_hist_2[i][j][1]=(TH1D*)_acceptance_7d[i][j]->Projection(1)->Clone();
			//std::cout<<"\nstep 4";
			accept_hist_04.push_back(accept_hist_03);
			accept_hist_03.clear();
		}
		//std::cout<<"\nstep 5";
		_accept_hist_1.push_back(accept_hist_04);
		accept_hist_04.clear();
	}
}
	//void Write_5d_Yield();
	//void Write_5d_Cross_Section();
	//void Write_5d_Holes();
	//void Write_7d_Holes();
	//void Write_Single_Diff();
	//void Write_Polarization();
	//void Write_Integrated();
void Histogram::Convert_to_Cross(){
	std::cout<<"Convert to Cross\n";
	std::cout<<"\tMaking X and Phi histograms with their bin widths\n";
	Histogram::XandPhi_BinHistograms();
	Histogram::Convert_Single_Diff_to_Cross();
	Histogram::Convert_Polarization_to_Cross();
}

void Histogram::Convert_Single_Diff_to_Cross(){
	std::cout<<"\tConvert Single Diff to Cross\n";
	for(int i =0; i<_n_var_sets; i++){//Variable Sets
		//_accept_hist_1.push_back(accept_hist_01);
		for(int j=0; j<_n_topology; j++){//Topology
			//_accept_hist_1[i].push_back(accept_hist_02);
			for(int k=0; k<_n_bins_7d[0]; k++){//W
				//_accept_hist_1[i][j].push_back(accept_hist_03);
				for(int l=0; l<_n_bins_7d[1]; l++){//Q2
					for(int m=0; m<4; m++){//Xij
						_exp_corr_holes_3d[i][j][k][l][m]->Scale(1./physics::Luminosity(_Q_tot_));
						_exp_corr_holes_3d[i][j][k][l][m]->Scale(1./physics::Virtual_Photon_Flux(Histogram::Bin_Center_7d(i,0,k),Histogram::Bin_Center_7d(i,0,k),_energy_e16_));
						//_exp_corr_holes_3d[i][j][k][l][m]->Scale(1./physics::Radiative_corr()); //We don't have this correction yet
						_exp_corr_holes_3d[i][j][k][l][m]->Scale(1./Histogram::W_Bin_Size(i,k));
						_exp_corr_holes_3d[i][j][k][l][m]->Scale(1./Histogram::Q2_Bin_Size(i,l));
						//if(i==j && j==k && k==l && l==0){
							//std::cout<<"\t\texp corr holes number of x bins:" <<_exp_corr_holes_3d[i][j][k][l][m]->GetNbinsX() <<"\n";
							//std::cout<<"\t\texp corr holes low edge:" <<_exp_corr_holes_3d[i][j][k][l][m]->GetXaxis()->GetBinLowEdge(0)<<"\n";
							//std::cout<<"\t\texp corr holes top edge:" <<_exp_corr_holes_3d[i][j][k][l][m]->GetXaxis()->GetBinUpEdge(_exp_corr_holes_3d[i][j][k][l][m]->GetNbinsX()-1)<<"\n";
							//std::cout<<"\t\tX bin size num x bins:" <<_X_bin_sizes[i][m]->GetNbinsX() <<"\n";
							//std::cout<<"\t\tX bin low edge:" <<_X_bin_sizes[i][m]->GetXaxis()->GetBinLowEdge(0) <<"\n";
							//std::cout<<"\t\tX bin top edge:" <<_X_bin_sizes[i][m]->GetXaxis()->GetBinUpEdge(_X_bin_sizes[i][m]->GetNbinsX()-1) <<"\n\n";
						//}
						if(_exp_corr_holes_3d[i][j][k][l][m]->GetNbinsX() == _X_bin_sizes[i][m]->GetNbinsX()){
							_exp_corr_holes_3d[i][j][k][l][m]->Divide(_X_bin_sizes[i][m]);
						}else{
							std::cout<<"\t\t\tBin sizes did not match: exp= " <<_exp_corr_holes_3d[i][j][k][l][m]->GetNbinsX() <<" X= "<<_X_bin_sizes[i][m]->GetNbinsX() <<"\n";
							std::cout<<"\t\t\tLow Ends: exp= " <<_exp_corr_holes_3d[i][j][k][l][m]->GetXaxis()->GetBinLowEdge(0) <<" X= "<<_X_bin_sizes[i][m]->GetXaxis()->GetBinLowEdge(0) <<"\n";
							std::cout<<"\t\t\tHigh Ends: exp= " <<_exp_corr_holes_3d[i][j][k][l][m]->GetXaxis()->GetBinUpEdge(_exp_corr_holes_3d[i][j][k][l][m]->GetNbinsX()-1) <<" X= "<<_X_bin_sizes[i][m]->GetXaxis()->GetBinUpEdge(_X_bin_sizes[i][m]->GetNbinsX()-1) <<"\n";
						}
					}
				}
			}
		}
	}
}	

void Histogram::Convert_Polarization_to_Cross(){
	std::cout<<"\tConvert Polarization to Cross\n";
	for(int i =0; i<_n_var_sets; i++){
		//_accept_hist_1.push_back(accept_hist_01);
		for(int j=0; j<_n_topology; j++){
			//_accept_hist_1[i].push_back(accept_hist_02);
			for(int k=0; k<_n_bins_7d[0]; k++){//W
				//_accept_hist_1[i][j].push_back(accept_hist_03);
				for(int l=0; l<_n_bins_7d[1]; l++){//Q2
					for(int m=0; m<4; m++){//Xij
						for(int n=0; n<_n_bins_5d[i][m]; n++){
							_exp_corr_holes_4d[i][j][k][l][m][n]->Scale(1./physics::Luminosity(_Q_tot_));
							_exp_corr_holes_4d[i][j][k][l][m][n]->Scale(1./physics::Virtual_Photon_Flux(Histogram::Bin_Center_7d(i,0,k),Histogram::Bin_Center_7d(i,0,k),_energy_e16_));
							//_exp_corr_holes_4d[i][j][k][l][m][n]->Scale(1./physics::Radiative_corr()); //We don't have this correction yet
							_exp_corr_holes_4d[i][j][k][l][m][n]->Scale(1./Histogram::W_Bin_Size(i,k));
							_exp_corr_holes_4d[i][j][k][l][m][n]->Scale(1./Histogram::Q2_Bin_Size(i,l));
							_exp_corr_holes_4d[i][j][k][l][m][n]->Scale(1./Histogram::X_Bin_Size(i,m,n));
							if(_exp_corr_holes_4d[i][j][k][l][m][n]->GetNbinsX() == _phi_bin_sizes[i]->GetNbinsX()){
								_exp_corr_holes_4d[i][j][k][l][m][n]->Divide(_phi_bin_sizes[i]);
							}else{
								std::cout<<"\t\t\t\tIndex: " <<i <<" " <<j <<" " <<k <<" " <<l <<" " <<m <<" " <<n <<"\n";
								std::cout<<"\t\t\tBin sizes did not match: exp= " <<_exp_corr_holes_4d[i][j][k][l][m][n]->GetNbinsX() <<" X= "<<_phi_bin_sizes[i]->GetNbinsX() <<"\n";
								std::cout<<"\t\t\tLow Ends: exp= " <<_exp_corr_holes_4d[i][j][k][l][m][n]->GetXaxis()->GetBinLowEdge(0) <<" X= "<<_phi_bin_sizes[i]->GetXaxis()->GetBinLowEdge(0) <<"\n";
								std::cout<<"\t\t\tHigh Ends: exp= " <<_exp_corr_holes_4d[i][j][k][l][m][n]->GetXaxis()->GetBinUpEdge(_exp_corr_holes_4d[i][j][k][l][m][n]->GetNbinsX()-1) <<" X= "<<_phi_bin_sizes[i]->GetXaxis()->GetBinUpEdge(_X_bin_sizes[i][m]->GetNbinsX()-1) <<"\n";
							}
						}
					}
				}
			}
		}
	}
}

void Histogram::Write_WQ2(){
	std::cout<<"Writing W vs Q2 Plots\n";
	TDirectory* dir_WQ2 = _RootOutputFile->mkdir("WQ2 Plots");
	dir_WQ2->cd();
	TDirectory* dir_wq2_1[_n_var_sets];
	TDirectory* dir_wq2_2[_n_var_sets][_n_topology];
	char dir_name[100];
	for(int i=0; i<_n_var_sets ; i++){
		sprintf(dir_name,"WQ2_Var:%s",_var_set_[i]);
		dir_wq2_1[i] = dir_WQ2->mkdir(dir_name);
		for(int j=0; j<_n_topology; j++){
			sprintf(dir_name,"WQ2_Top:%s_Var:%s",_top_[j],_var_set_[i]);
			dir_wq2_2[i][j] = dir_wq2_1[i]->mkdir(dir_name);
		}
	}

	for(int i=0; i<_n_var_sets ; i++){
		dir_wq2_1[i]->cd();
		_thr_hist_wq2[i]->Write();
		for(int j=0; j<_n_topology; j++){
			dir_wq2_2[i][j]->cd();
			_exp_hist_wq2[i][j]->Write();
			_sim_hist_wq2[i][j]->Write();
			_exp_corr_hist_wq2[i][j]->Write();
			_sim_corr_hist_wq2[i][j]->Write();
		}
	}
}

void Histogram::Write_Acceptance(){
	std::cout<<"Writing Acceptance Plots\n";
	TDirectory* dir_A = _RootOutputFile->mkdir("Acceptance Plots");
	dir_A->cd();
	//std::cout<<"\nPart 1";
	TDirectory* dir_A_1[_n_var_sets];
	TDirectory* dir_A_2[_n_var_sets][_n_topology];
	TDirectory* dir_A_3[_n_var_sets][_n_topology][_n_bins_7d[0]];
	TDirectory* dir_A_4[_n_var_sets][_n_topology][_n_bins_7d[0]][_n_bins_7d[1]];
	TDirectory* dir_A_3b[_n_var_sets][_n_topology];
	char dir_name[100];
	for(int i=0; i<_n_var_sets ; i++){
		//std::cout<<"\nPart 2";
		sprintf(dir_name,"Acceptance_Var:%s",_var_set_[i]);
		dir_A_1[i] = dir_A->mkdir(dir_name);
		for(int j=0; j<_n_topology; j++){
			//std::cout<<"\nPart 3";
			sprintf(dir_name,"Acceptance_Top:%s_Var:%s",_top_[j],_var_set_[i]);
			dir_A_2[i][j] = dir_A_1[i]->mkdir(dir_name);
			for(int l=0; l<_n_bins_7d[1]; l++){
				//std::cout<<"\nPart 4";
				sprintf(dir_name,"Acceptance_Top:%s_Var:%s_Q2:%f-%f",_top_[j],_var_set_[i],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l]);
				dir_A_3[i][j][l] = dir_A_2[i][j]->mkdir(dir_name);
				for(int k =0; k<_n_bins_7d[0]; k++){
					if(fun::nSparseIntegral(_acceptance_5d[i][j][k][l])>0.0){
						//std::cout<<"\nPart 5";
						//std::cout<<"\t index:" <<i <<" " <<j <<" " <<k <<" " <<l <<"\n";
						sprintf(dir_name,"Acceptance_Top:%s_Var:%s_W:%f-%f_Q2:%f-%f",_top_[j],_var_set_[i],_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l]);
						dir_A_4[i][j][l][k] = dir_A_3[i][j][k]->mkdir(dir_name);
					}
				}
			}
			//std::cout<<"\nPart 6";
			//sprintf(dir_name,"Acceptance_Top:%s_Var:%s_W&Q2",topo[j],var_set[i]);
			//dir_A_3b[i][j] = dir_A_2[i][j]->mkdir(dir_name);
		}
	}
	//std::cout<<"\nPart 7a";
	for(int i=0; i<_n_var_sets ; i++){
		//std::cout<<"\nPart 7b";
		dir_A_1[i]->cd();
		for(int j=0; j<_n_topology; j++){
			//std::cout<<"\nPart 7c";
			dir_A_2[i][j]->cd();
			for(int k=0; k<_n_bins_7d[1]; k++){
				//std::cout<<"\nPart 7d";
				dir_A_3[i][j][k]->cd();
				for(int l =0; l<_n_bins_7d[0]; l++){
					if(fun::nSparseIntegral(_acceptance_5d[i][j][k][l])>0.0){
						//std::cout<<"\nPart 7e";
						dir_A_4[i][j][k][l]->cd();
						//std::cout<<"\nPart 7";
					
						for(int m=0; m<5; m++){
							//std::cout<<"\nPart 8";
							_accept_hist_1[i][j][l][k][m]->Write();
						}
					}
				}
			}
			//std::cout<<"\nPart 9";
			if(fun::nSparseIntegral(_acceptance_7d[i][j])>0.0){
				//std::cout<<"\nPart 10";
				dir_A_2[i][j]->cd();
				_accept_hist_2[i][j][0]->Write();
				_accept_hist_2[i][j][1]->Write();
			}
		}
	}
	std::cout<<": Done\n";
}

void Histogram::Write_Single_Diff(){
	std::cout<<"Writing Single Differential Plots\n";
	TDirectory* dir_sd = _RootOutputFile->mkdir("Single Diff Plots");
	dir_sd->cd();
	TDirectory* dir_sd_1[_n_var_sets];
	TDirectory* dir_sd_2[_n_var_sets][_n_topology];
	TDirectory* dir_sd_3[_n_var_sets][_n_topology][4];
	TDirectory* dir_sd_4[_n_var_sets][_n_topology][4][_n_bins_7d[1]];
	char dir_name[100];
	for(int i=0; i<_n_var_sets ; i++){//Variable Set
		sprintf(dir_name,"Single_Diff_Var:%s",_var_set_[i]);
		dir_sd_1[i] = dir_sd->mkdir(dir_name);
		for(int j=0; j<_n_topology; j++){//Topology
			sprintf(dir_name,"Single_Diff_Var:%s_Top:%s",_var_set_[i],_top_[j]);
			dir_sd_2[i][j] = dir_sd_1[i]->mkdir(dir_name);
			for(int k=0; k<4; k++){//Xij
				sprintf(dir_name,"Single_Diff_Var:%s_Top:%s_X:%s",_var_set_[i],_top_[j],_five_dim_[k]);
				dir_sd_3[i][j][k] = dir_sd_2[i][j]->mkdir(dir_name);
				for(int l=0; l<_n_bins_7d[1]; l++){//Q2
					sprintf(dir_name,"Single_Diff_Var:%s_Top:%s_X:%s_Q2:%f-%f",_var_set_[i],_top_[j],_five_dim_[k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l]);
					dir_sd_4[i][j][k][l] = dir_sd_3[i][j][k]->mkdir(dir_name);
				}
			}
			
		}
	}
	for(int i=0; i<_n_var_sets ; i++){
		dir_sd_1[i]->cd();
		for(int j=0; j<_n_topology; j++){
			dir_sd_2[i][j]->cd();
			for(int k=0; k<4; k++){
				dir_sd_3[i][j][k]->cd();
				for(int l=0; l<_n_bins_7d[1]; l++){
					dir_sd_4[i][j][k][l]->cd();
					for(int m=0; m<_n_bins_7d[0]; m++){
						_exp_corr_holes_3d[i][j][m][l][k]->Write();
					}
				}
			}
			
		}
	}
}

void Histogram::Write_Polarization(){
	std::cout<<"Writing Polarization Plots\n";
	TDirectory* dir_pol = _RootOutputFile->mkdir("Polarization Plots");
	dir_pol->cd();
	TDirectory* dir_pol_1[_n_var_sets];
	TDirectory* dir_pol_2[_n_var_sets][_n_topology];
	TDirectory* dir_pol_3[_n_var_sets][_n_topology][_n_bins_7d[0]];
	TDirectory* dir_pol_4[_n_var_sets][_n_topology][_n_bins_7d[0]][_n_bins_7d[1]];
	TDirectory* dir_pol_5[_n_var_sets][_n_topology][_n_bins_7d[0]][_n_bins_7d[1]][4];
	char dir_name[100];
	std::cout<<"\tNesting Polarization Directories\n";
	for(int i=0; i<_n_var_sets ; i++){//Variable Set
		sprintf(dir_name,"Polarization_Var:%s",_var_set_[i]);
		dir_pol_1[i] = dir_pol->mkdir(dir_name);
		for(int j=0; j<_n_topology; j++){//Topology
			sprintf(dir_name,"Polarization_Var:%s_Top:%s",_var_set_[i],_top_[j]);
			dir_pol_2[i][j] = dir_pol_1[i]->mkdir(dir_name);
			for(int k=0; k<_n_bins_7d[0]; k++){//W
				sprintf(dir_name,"Polarization_Var:%s_Top:%s_W:%f-%f",_var_set_[i],_top_[j],_bin_low_7d[i][0][k],_bin_up_7d[i][0][k]);
				dir_pol_3[i][j][k] = dir_pol_2[i][j]->mkdir(dir_name);
				for(int l=0; l<_n_bins_7d[1]; l++){//Q2
					sprintf(dir_name,"Polarization_Var:%s_Top:%s_W:%f-%f_Q2:%f-%f",_var_set_[i],_top_[j],_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l]);
					dir_pol_4[i][j][k][l] = dir_pol_3[i][j][k]->mkdir(dir_name);
					for(int m=0; m<4; m++){//Xij
						sprintf(dir_name,"Polarization_X:%s_Var:%s_Top:%s_W:%f-%f_Q2:%f-%f",_five_dim_[m],_var_set_[i],_top_[j],_bin_low_7d[i][0][k],_bin_up_7d[i][0][k],_bin_low_7d[i][1][l],_bin_up_7d[i][1][l]);
						dir_pol_5[i][j][k][l][m] = dir_pol_4[i][j][k][l]->mkdir(dir_name);
					}
				}
			}
			
		}
	}
	std::cout<<"\tWriting Polarization Histograms\n";
	for(int i=0; i<_n_var_sets ; i++){//var
		dir_pol_1[i]->cd();
		for(int j=0; j<_n_topology; j++){//top
			dir_pol_2[i][j]->cd();
			for(int k=0; k<_n_bins_7d[0]; k++){//W
				dir_pol_3[i][j][k]->cd();
				for(int l=0; l<_n_bins_7d[1]; l++){//Q2
					dir_pol_4[i][j][k][l]->cd();
					for(int m=0; m<4; m++){//Xij
						dir_pol_5[i][j][k][l][m]->cd();
						for(int n=0; n<_n_bins_5d[i][m]; n++){//Xij bins
							_exp_corr_holes_4d[i][j][k][l][m][n]->Write();
						}
					}
				}
			}
			
		}
	}
	std::cout<<"Done Writing Polarization Hists\n";
}

float Histogram::Bin_Size(int var_set, int variable, int bin_7d){
	return _bin_up_7d[var_set][variable][bin_7d] - _bin_low_7d[var_set][variable][bin_7d];
}

float Histogram::W_Bin_Size(int var_set, int bin_7d){
	return Histogram::Bin_Size(var_set,0,bin_7d);
}

float Histogram::Q2_Bin_Size(int var_set, int bin_7d){
	return Histogram::Bin_Size(var_set,1,bin_7d);
}

float Histogram::X_Bin_Size(int var_set, int variable, int bin_7d){
	return Histogram::Bin_Size(var_set,variable+2,bin_7d);
}

float Histogram::Phi_Bin_Size(int var_set, int bin_7d){
	return Histogram::Bin_Size(var_set,6,bin_7d);
}

double Histogram::Bin_Center_7d(int var_set, int dim, int bin){
	return (_bin_up_7d[var_set][dim][bin] + _bin_low_7d[var_set][dim][bin])/2.0;
}

void Histogram::XandPhi_BinHistograms(){
	std::cout<<"\tX and Phi Bin Histograms\n";
	char hname[100];

	TH1D_1d_star xbinsizes;
	TH1D* phi_bin_size;
	//Make the histograms
	
	for(int var =0; var<_n_var_sets; var++){
		sprintf(hname,"phi_bin_sizes_var: %s",_var_set_[var]);
		std::cout<<"var: " <<var <<"\n\tLowest phi: " <<_bin_low_5d[var][4][0] <<" | highest phi: " <<_bin_up_5d[var][4][_n_bins_5d[var][4]-1];
		_phi_bin_sizes.push_back(new TH1D(hname,hname,_n_bins_5d[var][4],_bin_low_5d[var][4][0],_bin_up_5d[var][4][_n_bins_5d[var][4]-1]));
		for(int j=0; j<_n_bins_5d[var][4]; j++){
			_phi_bin_sizes[var]->Fill(j,Histogram::Phi_Bin_Size(var,j));
		}
		for(int x=0; x<4; x++){
			std::cout<<"var: " <<var <<" X: " <<x <<"\n\tLowest X: " <<_bin_low_5d[var][x][0] <<" | highest x: " <<_bin_up_5d[var][x][_n_bins_5d[var][x]-1];
			sprintf(hname,"%s_bin_sizes_%s",_five_dim_[x],_var_set_[var]);
			xbinsizes.push_back(new TH1D(hname,hname,_n_bins_5d[var][x],_bin_low_5d[var][x][0],_bin_up_5d[var][x][_n_bins_5d[var][x]-1]));
			for(int i=0; i<_n_bins_5d[var][x]; i++){
				xbinsizes[x]->Fill(i,Histogram::X_Bin_Size(var,x,i));
			}
		}
		_X_bin_sizes.push_back(xbinsizes);
		xbinsizes.clear();
	}
}


void Histogram::Acceptance_Errors(){
	int bin_7d[7];
	for(int var=0; var<_n_var_sets; var++){
		for(int top=0; top<_n_topology; top++){
			//_acceptance_eff_7d[var][top]=(THnSparseD*)_acceptance_7d[var][top]->Clone();
			_acceptance_err_7d[var][top][0]=(THnSparseD*)_acceptance_7d[var][top]->Clone();//Unweighted
			_acceptance_err_7d[var][top][1]=(THnSparseD*)_acceptance_7d[var][top]->Clone();//Weighted
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
										//"Unweighted" Error
										_acceptance_err_7d[var][top][0]->SetBinContent(bin_7d,physics::Error(_thrown_7d[var]->THnSparse::GetBinContent(bin_7d),_sim_data_7d[var][top]->GetBinContent(bin_7d)));
										//Weighted Error
										_acceptance_err_7d[var][top][1]->SetBinContent(bin_7d,physics::Error(_thrown_7d[var]->THnSparse::GetBinContent(bin_7d),_sim_data_7d[var][top]->GetBinContent(bin_7d),_sim_weight_sq_7d[var][top]->GetBinContent(bin_7d)));
										//Efficiency
										//_acceptance_eff_7d[var][top]->SetBinContent(bin_7d,_sim_data_7d[var][top]->GetBinContent(bin_7d)/_thrown_7d[var]->THnSparse::GetBinContent(bin_7d));
										//_acceptance_eff_7d[var][top]->SetBinError(bin_7d,_acceptance_err_7d[var][top][1]);//Give it the weighted error
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void Histogram::Make_Error_Hists(){
	char hname[100];

	//TH1D_1d_star _zero_1d;

	TH2D_1d_star _ratio_1d;

	TH1D_1d_star _rel_err_w_1d;
	TH1D_2d_star _rel_err_w_2d;
	TH1D_3d_star _rel_err_w_3d;

	TH1D_1d_star _rel_err_nw_1d;
	TH1D_2d_star _rel_err_nw_2d;
	TH1D_3d_star _rel_err_nw_3d;

	TH1D_3d_star _yield_3d;
	TH1D_2d_star _yield_2d;
	TH1D_1d_star _yield_1d;

	//TH1D_3d_star _zero_thr_3d;
	//TH1D_2d_star _zero_thr_2d;
	//TH1D_1d_star _zero_thr_1d;
	for(int var=0; var<_n_var_sets; var++){
		for(int top=0; top<_n_topology; top++){
			//sprintf(hname,"Fraction_Acc_Zero_%s_%s",_var_set_[var],_top_[top]);
			//_zero_1d.push_back(new TH1D(hname,hname,101,-0.005,1.005));
			sprintf(hname,"Thrown_Exp_Ratio_%s_%s",_var_set_[var],_top_[top]);
			_ratio_1d.push_back(new TH2D(hname,hname,_n_bins_7d[0],_acceptance_7d[var][top]->GetAxis(0)->GetBinLowEdge(0),_acceptance_7d[var][top]->GetAxis(0)->GetBinUpEdge(_n_bins_7d[0]-1),_n_bins_7d[1],_acceptance_7d[var][top]->GetAxis(1)->GetBinLowEdge(0),_acceptance_7d[var][top]->GetAxis(1)->GetBinUpEdge(_n_bins_7d[1]-1)));
			for(int i=0; i<_n_bins_7d[0]; i++){//W
				for(int j=0; j< _n_bins_7d[1]; j++){//Q2
					sprintf(hname,"Acc_Yield_%s_%s_W:%f-%f_Q:%f_%F",_var_set_[var],_top_[top],_acceptance_7d[var][top]->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d[var][top]->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d[var][top]->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d[var][top]->GetAxis(1)->GetBinUpEdge(j));
					_yield_1d.push_back(new TH1D(hname,hname,1001,-0.001,0.2));
					//sprintf(hname,"Fraction_Acc_Zero_%s_%s_W:%f-%f_Q:%f_%F",_var_set_[var],_top_[top],_acceptance_7d->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d->GetAxis(1)->GetBinUpEdge(j));
					//_zero_thr_1d.push_back(new TH1D(hname,hname,1001,-0.001,0.2));
					sprintf(hname,"Relative_Acc_Error_Weighted_%s_%s_W:%f-%f_Q:%f_%F",_var_set_[var],_top_[top],_acceptance_7d[var][top]->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d[var][top]->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d[var][top]->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d[var][top]->GetAxis(1)->GetBinUpEdge(j));
					_rel_err_w_1d.push_back(new TH1D(hname,hname,1401,-0.005,1.405));
					sprintf(hname,"Relative_Acc_Error_Unweighted_%s_%s_W:%f-%f_Q:%f_%F",_var_set_[var],_top_[top],_acceptance_7d[var][top]->GetAxis(0)->GetBinLowEdge(i),_acceptance_7d[var][top]->GetAxis(0)->GetBinUpEdge(i),_acceptance_7d[var][top]->GetAxis(1)->GetBinLowEdge(j),_acceptance_7d[var][top]->GetAxis(1)->GetBinUpEdge(j));
					_rel_err_nw_1d.push_back(new TH1D(hname,hname,1401,-0.005,1.405));
					if(j==_n_bins_7d[1]-1){
						if(_yield_1d.size()>0){
							_yield_2d.push_back(_yield_1d);
							_yield_1d.clear();
						}
						/*if(_zero_thr_1d.size()>0){
							_zero_thr_2d.push_back(_zero_thr_1d);
							_zero_thr_1d.clear();
						}*/
						if(_rel_err_nw_1d.size()>0){
							_rel_err_nw_2d.push_back(_rel_err_w_1d);
							_rel_err_nw_1d.clear();
						}
					}
				}
				if(i==_n_bins_7d[0]-1){
					if(_yield_2d.size()>0){
						_yield_3d.push_back(_yield_2d);
						_yield_2d.clear();
					}
					/*if(_zero_thr_2d.size()>0){
						_zero_thr_3d.push_back(_zero_thr_2d);
						_zero_thr_2d.clear();
					}*/
					if(_rel_err_w_2d.size()>0){
						_rel_err_w_3d.push_back(_rel_err_w_2d);
						_rel_err_w_2d.clear();
					}
					if(_rel_err_nw_2d.size()>0){
						_rel_err_nw_3d.push_back(_rel_err_nw_2d);
						_rel_err_nw_2d.clear();
					}
				}
			}
			if(top==_n_topology-1){
				if(_yield_3d.size()>0){
					_acc_yield.push_back(_yield_3d);
					_yield_3d.clear();
				}
				/*if(_zero_thr_3d.size()>0){
					_acc_zero_thr.push_back(_zero_thr_3d);
					_zero_thr_3d.clear();
				}*/
				if(_rel_err_w_3d.size()>0){
					_acc_rel_error_weighted.push_back(_rel_err_w_3d);
					_rel_err_w_3d.clear();
				}
				if(_rel_err_nw_3d.size()>0){
					_acc_rel_error_unweighted.push_back(_rel_err_nw_3d);
					_rel_err_nw_3d.clear();
				}
				/*if(_zero_1d.size()>0){
					_acc_zero_exp.push_back(_zero_1d);
					_zero_1d.clear();
				}*/
				if(_ratio_1d.size()>0){
					_thr_exp_ratio.push_back(_ratio_1d);
					_ratio_1d.clear();
				}
			}
		}
	}
}

void Histogram::Fill_Error_Hists(){
	std::cout<<"Filling Acceptance Error and Efficiency Histograms\n";
	int bin_7d[7];
	int bin_5d[5];
	for(int var=0; var<_n_var_sets; var++){
		for(int top=0; top<_n_topology; top++){
			for(int i=0; i<_n_bins_7d[0]; i++){//W
				bin_7d[0]=i;
				for(int j=0; j< _n_bins_7d[1]; j++){//Q2
					bin_7d[1]=j;
					_thr_exp_ratio[var][top]->Fill(Histogram::Bin_Center_7d(var,0,i),Histogram::Bin_Center_7d(var,0,i),fun::nSparseIntegral(_thrown_5d[var][i][j])/fun::nSparseIntegral(_exp_data_5d[var][top][i][j]));
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
										if(_acceptance_7d[var][top]->GetBinContent(bin_7d)==0){
											if(_exp_data_7d[var][top]->GetBinContent(bin_7d)>0.0){
												//_acc_zero_exp[var][top]->Fill(1.0);
												_acc_yield[var][top][i][j]->Fill(_acceptance_5d[var][top][i][j]->GetBinContent(bin_5d));
											}
											if(_thrown_7d[var][top].GetBinContent(bin_7d)>0.0){
												//_acc_zero_thr[var][top][i][j]->Fill(1.0);
												_acc_rel_error_weighted[var][top][i][j]->Fill(_acceptance_err_7d[var][top][1]->GetBinContent(bin_7d)/_acceptance_7d[var][top]->GetBinContent(bin_7d));
												_acc_rel_error_unweighted[var][top][i][j]->Fill(_acceptance_err_7d[var][top][0]->GetBinContent(bin_7d)/_acceptance_7d[var][top]->GetBinContent(bin_7d));
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
}

void Histogram::Write_Error_Hists(){
	char dirname[100];
	TDirectory* dir_err = _RootOutputFile->mkdir("Acceptance Error Plots");
	//TDirectory* dir_err_zero = dir_err->mkdir("Zero Acceptance");
	//TDirectory* dir_err_zero_exp = dir_err_zero->mkdir("Zero Acceptance Exp Non-Zero");
	//TDirectory* dir_err_zero_exp_sub[_n_var_sets][_n_topology+1];
	//TDirectory* dir_err_zero_thr = dir_err_zero->mkdir("Zero Acceptance Thrown Non-Zero");
	//TDirectory* dir_err_zero_thr_sub[_n_var_sets][_n_topology+1][_n_bins_7d[0]+1][_n_bins_7d[1]+1];
	TDirectory* dir_err_yield = dir_err->mkdir("Acceptance Yield Distribution");
	TDirectory* dir_err_yield_sub[_n_var_sets][_n_topology+1][_n_bins_7d[0]+1][_n_bins_7d[1]+1];
	TDirectory* dir_err_ratio = dir_err->mkdir("Thrown to Exp Ratio");
	TDirectory* dir_err_ratio_sub[_n_var_sets][_n_topology+1];
	TDirectory* dir_err_rel = dir_err->mkdir("Relative Error");
	TDirectory* dir_err_rel_sub[_n_var_sets][_n_topology+1][_n_bins_7d[0]+1][_n_bins_7d[1]+1];
	for(int var=0; var<_n_var_sets; var++){
		//sprintf(dirname,"Zero Acceptance Exp Non-Zero %s",_var_set_[var]);
		//dir_err_zero_exp_sub[var][0] = dir_err_zero_exp->mkdir(dirname);
		//sprintf(dirname,"Zero Acceptance Thr Non-Zero %s",_var_set_[var]);
		//dir_err_zero_thr_sub[var][0][0][0] = dir_err_zero_thr->mkdir(dirname);
		sprintf(dirname,"Acceptance Yield Distribution %s",_var_set_[var]);
		dir_err_yield_sub[var][0][0][0] = dir_err_yield->mkdir(dirname);
		sprintf(dirname,"Thrown to Exp Ratio %s",_var_set_[var]);
		dir_err_ratio_sub[var][0] = dir_err_yield->mkdir(dirname);
		sprintf(dirname,"Relative Error %s",_var_set_[var]);
		dir_err_rel_sub[var][0][0][0] = dir_err_rel->mkdir(dirname);
		for(int top=0; top<_n_topology; top++){
			/*sprintf(dirname,"Zero Acceptance Exp Non-Zero %s",_var_set_[var]);
			dir_err_zero_exp_sub[var][top+1] = dir_err_zero_exp_sub[var][0]->mkdir(dirname);
			dir_err_zero_exp_sub[var][top+1]->cd();
			_acc_zero_exp[var][top]->SetXTitle("Fraction of Sim Performed");
			_acc_zero_exp[var][top]->SetYTitle("Fraction of Bins of Zero");
			_acc_zero_exp[var][top]->Write();
			sprintf(dirname,"Zero Acceptance Thr Non-Zero %s",_var_set_[var]);
			dir_err_zero_thr_sub[var][top+1][0][0] = dir_err_zero_thr_sub[var][0][0][0]->mkdir(dirname);*/
			sprintf(dirname,"Acceptance Yield Distribution %s",_var_set_[var]);
			dir_err_yield_sub[var][top+1][0][0] = dir_err_yield_sub[var][0][0][0]->mkdir(dirname);
			sprintf(dirname,"Thrown to Exp Ratio %s",_var_set_[var]);
			dir_err_ratio_sub[var][top+1] = dir_err_ratio_sub[var][0]->mkdir(dirname);
			dir_err_ratio_sub[var][top+1]->cd();
			_thr_exp_ratio[var][top]->SetXTitle("W (GeV)");
			_thr_exp_ratio[var][top]->SetYTitle("Q2 (GeV^2)");
			_thr_exp_ratio[var][top]->Write();
			sprintf(dirname,"Relative Error %s",_var_set_[var]);
			dir_err_rel_sub[var][top+1][0][0] = dir_err_rel_sub[var][0][0][0]->mkdir(dirname);
			for(int i=0; i<_n_bins_7d[0]; i++){
				//sprintf(dirname,"Zero Acceptance Thr Non-Zero %s",_var_set_[var]);
				//dir_err_zero_thr_sub[var][top+1][i+1][0] = dir_err_zero_thr_sub[var][top+1][0][0]->mkdir(dirname);
				sprintf(dirname,"Acceptance Yield Distribution %s",_var_set_[var]);
				dir_err_yield_sub[var][top+1][i+1][0] = dir_err_yield_sub[var][top+1][0][0]->mkdir(dirname);
				sprintf(dirname,"Relative Error %s",_var_set_[var]);
				dir_err_rel_sub[var][top+1][i+1][0] = dir_err_rel_sub[var][top+1][0][0]->mkdir(dirname);
				for(int j=0; j<_n_bins_7d[1]; j++){
					/*sprintf(dirname,"Zero Acceptance Thr Non-Zero %s",_var_set_[var]);
					dir_err_zero_thr_sub[var][top+1][i+1][j+1] = dir_err_zero_thr_sub[var][top+1][i+1][0]->mkdir(dirname);
					dir_err_zero_thr_sub[var][top+1][i+1][j+1]->cd();
					_acc_zero_thr[var][top][i][j]->SetXTitle("Fraction of Sim Performed");
					_acc_zero_thr[var][top][i][j]->SetYTitle("Fraction of Bins of Zero");
					_acc_zero_thr[var][top][i][j]->Write();*/
					sprintf(dirname,"Acceptance Yield Distribution %s",_var_set_[var]);
					dir_err_yield_sub[var][top+1][i+1][j+1] = dir_err_yield_sub[var][top+1][i+1][0]->mkdir(dirname);
					dir_err_yield_sub[var][top+1][i+1][j+1]->cd();
					_acc_yield[var][top][i][j]->SetXTitle("Acceptance Yield (5d)");
					_acc_yield[var][top][i][j]->SetYTitle("Number of Bins");
					_acc_yield[var][top][i][j]->Write();
					sprintf(dirname,"Relative Error %s",_var_set_[var]);
					dir_err_rel_sub[var][top+1][i+1][j+1] = dir_err_rel_sub[var][top+1][i+1][0]->mkdir(dirname);
					dir_err_rel_sub[var][top+1][i+1][j+1]->cd();
					_acc_rel_error_weighted[var][top][i][j]->SetXTitle("Relative Weighted Error");//{var,top,W,Q2}
					_acc_rel_error_weighted[var][top][i][j]->SetYTitle("Number of Bins");
					_acc_rel_error_weighted[var][top][i][j]->Write();
					_acc_rel_error_unweighted[var][top][i][j]->SetXTitle("Relative Unweighted Error");//{var,top,W,Q2}
					_acc_rel_error_unweighted[var][top][i][j]->SetYTitle("Number of Bins");
					_acc_rel_error_unweighted[var][top][i][j]->Write();
				}
			}
		}
	}
}

























