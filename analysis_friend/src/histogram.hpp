#ifndef HISTOGRAM_H_GUARD
#define HISTOGRAM_H_GUARD

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "constants.hpp"
#include "functions.hpp"
//#include "CartesianGenerator.hpp"
//#include "physics.hpp"
#include "THnSparse.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include "physics.hpp"
#include "TLatex.h"
#include "flags.hpp"

using Sparse_5d = std::vector<std::vector<std::vector<std::vector<std::vector<THnSparseD>>>>>;
using Sparse_4d = std::vector<std::vector<std::vector<std::vector<THnSparseD>>>>;
using Sparse_3d = std::vector<std::vector<std::vector<THnSparseD>>>;
using Sparse_2d = std::vector<std::vector<THnSparseD>>;
using Sparse_1d =std::vector<THnSparseD>;

using Sparse_5d_star = std::vector<std::vector<std::vector<std::vector<std::vector<THnSparseD*>>>>>;
using Sparse_4d_star = std::vector<std::vector<std::vector<std::vector<THnSparseD*>>>>;
using Sparse_3d_star = std::vector<std::vector<std::vector<THnSparseD*>>>;
using Sparse_2d_star = std::vector<std::vector<THnSparseD*>>;
using Sparse_1d_star =std::vector<THnSparseD*>;
//using float_3d = std::vector<std::vector<std::vector<float>>>;
using double_4d = std::vector<std::vector<std::vector<std::vector<double>>>>;
using double_3d = std::vector<std::vector<std::vector<double>>>;
using double_2d = std::vector<std::vector<double>>;
using double_1d = std::vector<double>;

using TH2D_5d = std::vector<std::vector<std::vector<std::vector<std::vector<TH2D>>>>>;
using TH2D_4d = std::vector<std::vector<std::vector<std::vector<TH2D>>>>;
using TH2D_3d = std::vector<std::vector<std::vector<TH2D>>>;
using TH2D_2d = std::vector<std::vector<TH2D>>;
using TH2D_1d = std::vector<TH2D>;

using TH2D_5d_star = std::vector<std::vector<std::vector<std::vector<std::vector<TH2D*>>>>>;
using TH2D_4d_star = std::vector<std::vector<std::vector<std::vector<TH2D*>>>>;
using TH2D_3d_star = std::vector<std::vector<std::vector<TH2D*>>>;
using TH2D_2d_star = std::vector<std::vector<TH2D*>>;
using TH2D_1d_star = std::vector<TH2D*>;

using TH1D_6d = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<TH1D>>>>>>;
using TH1D_5d = std::vector<std::vector<std::vector<std::vector<std::vector<TH1D>>>>>;
using TH1D_4d = std::vector<std::vector<std::vector<std::vector<TH1D>>>>;
using TH1D_3d = std::vector<std::vector<std::vector<TH1D>>>;
using TH1D_2d = std::vector<std::vector<TH1D>>;
using TH1D_1d = std::vector<TH1D>;

using TH1F_6d = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<TH1F>>>>>>;
using TH1F_5d = std::vector<std::vector<std::vector<std::vector<std::vector<TH1F>>>>>;
using TH1F_4d = std::vector<std::vector<std::vector<std::vector<TH1F>>>>;
using TH1F_3d = std::vector<std::vector<std::vector<TH1F>>>;
using TH1F_2d = std::vector<std::vector<TH1F>>;
using TH1F_1d = std::vector<TH1F>;

using TH1D_6d_star = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<TH1D*>>>>>>;
using TH1D_5d_star = std::vector<std::vector<std::vector<std::vector<std::vector<TH1D*>>>>>;
using TH1D_4d_star = std::vector<std::vector<std::vector<std::vector<TH1D*>>>>;
using TH1D_3d_star = std::vector<std::vector<std::vector<TH1D*>>>;
using TH1D_2d_star = std::vector<std::vector<TH1D*>>;
using TH1D_1d_star = std::vector<TH1D*>;

using TH1F_6d_star = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*>>>>>>;
using TH1F_5d_star = std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*>>>>>;
using TH1F_4d_star = std::vector<std::vector<std::vector<std::vector<TH1F*>>>>;
using TH1F_3d_star = std::vector<std::vector<std::vector<TH1F*>>>;
using TH1F_2d_star = std::vector<std::vector<TH1F*>>;
using TH1F_1d_star = std::vector<TH1F*>;

using TH2F_6d_star = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<TH2F*>>>>>>;
using TH2F_5d_star = std::vector<std::vector<std::vector<std::vector<std::vector<TH2F*>>>>>;
using TH2F_4d_star = std::vector<std::vector<std::vector<std::vector<TH2F*>>>>;
using TH2F_3d_star = std::vector<std::vector<std::vector<TH2F*>>>;
using TH2F_2d_star = std::vector<std::vector<TH2F*>>;
using TH2F_1d_star = std::vector<TH2F*>;

static const float _mp_ = 0.93828;	//Mass of proton in GeV
static const float _mpi_ = 0.1395; //Mass of Pion in GeV

//Copied from the Event Selection 
static int _n_var_ = 3; //Number of variable sets
static float _W_min_ = 1.4;
static float _W_max_ = 2.125;
static float _W_res_ = 0.025;
static float _Q2_min_ = 2.0;
static float _Q2_max_ = 5.0;
static float _Q2_bins_[6] = {2.0,2.4,3.0,3.5,4.2,5.0};
static float _MM_min_[3] = {_W_min_-_mpi_,_W_min_-_mp_,_W_min_-_mpi_};//{1.1,0.3,1.1};//Changed 6/19/23
static float _MM_max_[3] = {_W_max_-_mpi_,_W_max_-_mp_,_W_max_-_mpi_};//{2.0,1.1,2.0};//Changed 6/19/23
static int _MM_bins_ = 14;
static float _MM2_min_[3] = {_W_min_-_mp_,_W_min_-_mpi_,_W_min_-_mpi_};//{0.3,1.1,1.1};//Changed 6/19/23
static float _MM2_max_[3] = {_W_max_-_mp_,_W_max_-_mpi_,_W_max_-_mpi_};
static int _theta_bins_ = 10;
static float _theta_min_ = 0.0;
static float _theta_max_ = 180.0;
static int _alpha_bins_ = 10;
static float _alpha_min_ = 0.0;
static float _alpha_max_ = 360.;
static int _phi_bins_ = 10;
static float _phi_min_ = 0.0; 
static float _phi_max_ = 360.0;
static char * _bin_names_[] = {"W","Q2","MM1","MM2","theta","alpha","phi"};


class Histogram{
protected:
	//RootFiles
	TFile* _RootOutputFile;
	//TFile* _RootOutputFile1;//Beam Spin and others
	//TFile* _RootOutputFile2;//Single Diff 
	//TFile* _RootOutputFile3;//Polarization

	bool _localized_hole_filling=false;

	//Sparse Binning
	std::vector<int> _n_bins_7d;//Number of bins in a given dimension {W,Q2,MM1,MM2,theta,alpha,phi}
	std::vector<int> _n_bins_5d; //Number of bins in a given dimension {MM1,MM2,theta,alpha,phi}
	double_2d _bin_low_7d;//Low edge of a bin  {variable,bin # inside variable}
	double_2d _bin_up_7d;//Top Edge of a bin {variable,bin # inside variable}
	double_2d _bin_lowside_7d; //low values of the bins {variable,bin # inside variable}
	double_2d _bin_topside_7d; //top values of the bins {variable,bin # inside variable}
	double_2d _bin_mid_7d;//Middle value of a bin {variable,bin # inside variable}
	double_2d _bin_size_7d;//Width of a bin {variable,bin # inside variable}
	double_2d _bin_edges_7d;//Edges of bins {variable,bin # inside variable}
	double_2d _bin_low_5d;//Low edge of a bin {variable,bin # inside variable}
	double_2d _bin_up_5d;//Top Edge of a bin {variable,bin # inside variable}
	double_2d _bin_mid_5d;//Middle value of a bin {variable,bin # inside variable}
	double_2d _bin_size_5d;//Width of a bin {variable,bin # inside variable}
	double_2d _bin_edges_5d;//Edges of a bin {variable,bin # inside variable}
	double_2d _bin_lowside_5d;// {variable,bin # inside variable}
	double_2d _bin_topside_5d;// {variable,bin # inside variable}


	//Sparse Histograms
	//Experimental Reconstruction
		//Target Filled
			//CC Efficiency Applied
	THnSparseD *_exp_data_7d;
	THnSparseD *_exp_data_7d_pos;//Positive Helicity
	THnSparseD *_exp_data_7d_neg;//Negative Helicity
	THnSparseD *_exp_data2_7d;
	THnSparseD *_exp_data2_7d_pos;//Positive Helicity
	THnSparseD *_exp_data2_7d_neg;//Negative Helicity
	Sparse_2d_star _exp_data_5d;
	Sparse_2d_star _exp_data_5d_pos;//Positive Helicity
	Sparse_2d_star _exp_data_5d_neg;//Negative Helicity
			//Raw Yields
	//THnSparseD *_exp_data_7d_y;
	//THnSparseD *_exp_data_7d_pos_y;//Positive Helicity
	//THnSparseD *_exp_data_7d_neg_y;//Negative Helicity
	//Sparse_2d_star _exp_data_5d_y;
	//Sparse_2d_star _exp_data_5d_pos_y;//Positive Helicity
	//Sparse_2d_star _exp_data_5d_neg_y;//Negative Helicity
			//Summed CC Efficiencies by bin
	//THnSparseD *_exp_weight_7d;
	//THnSparseD *_exp_weight_7d_pos;//Positive Helicity
	//THnSparseD *_exp_weight_7d_neg;//Negative Helicity
	//Sparse_2d_star _exp_weight_5d;
	//Sparse_2d_star _exp_weight_5d_pos;//Positive Helicity
	//Sparse_2d_star _exp_weight_5d_neg;//Negative Helicity
		//Target Empty
			//CC Efficiency Applied
	THnSparseD *_empty_7d;
	THnSparseD *_empty_7d_pos;//Positive Helicity
	THnSparseD *_empty_7d_neg;//Negative Helicity
	THnSparseD *_empty2_7d;
	THnSparseD *_empty2_7d_pos;//Positive Helicity
	THnSparseD *_empty2_7d_neg;//Negative Helicity
	Sparse_2d_star _empty_5d;
	Sparse_2d_star _empty_5d_pos;//Positive Helicity
	Sparse_2d_star _empty_5d_neg;//Negative Helicity
			//Raw Yields
	//THnSparseD *_empty_7d_y;
	//THnSparseD *_empty_7d_pos_y;//Positive Helicity
	//THnSparseD *_empty_7d_neg_y;//Negative Helicity
	//Sparse_2d_star _empty_5d_y;
	//Sparse_2d_star _empty_5d_pos_y;//Positive Helicity
	//Sparse_2d_star _empty_5d_neg_y;//Negative Helicity
			//Summed CC Efficiencies by bin
	//THnSparseD *_empty_weight_7d;
	//THnSparseD *_empty_weight_7d_pos;//Positive Helicity
	//THnSparseD *_empty_weight_7d_neg;//Negative Helicity
	//Sparse_2d_star _empty_weight_5d;
	//Sparse_2d_star _empty_weight_5d_pos;//Positive Helicity
	//Sparse_2d_star _empty_weight_5d_neg;//Negative Helicity
	//Simulation
		//Thrown
	THnSparseD *_thrown_7d;//Weighted
	THnSparseD *_thrown2_7d;//Weighted
	//THnSparseD *_thrown_7d_y;//Unweighted
	THnSparseD *_thrown_7d_no_rad;//No radiative effects, weighted
	THnSparseD *_thrown2_7d_no_rad;//No radiative effects, weighted
	//THnSparseD *_thrown_weight_7d;//Summed weight per bin
	Sparse_2d_star _thrown_5d;//Weighted
	//Sparse_2d_star _thrown_5d_y;//Unweighted
	Sparse_2d_star _thrown_5d_no_rad;//No radiative effects, weighted
	//Sparse_2d_star _thrown_weight_5d;//Summed weight per bin
		//Reconstructed
	THnSparseD *_sim_data_7d;//Weighted
	THnSparseD *_sim_data2_7d;//Weighted
	//THnSparseD *_sim_data_7d_y;//Unweighted
	THnSparseD * _sim_holes_7d;//Simulation Holes
	THnSparseD * _sim_holes_tmp_7d;
	Sparse_2d_star _sim_holes_5d;
	//THnSparseD *_sim_weight_7d;//summed weights for each 7d bin for simulation
	Sparse_2d_star _sim_data_5d;//Weighted
	//Sparse_2d_star _sim_data_5d_y;//Unweighted
	//Sparse_2d_star _sim_weight_5d;//summed weights for each 7d bin for simulation

	THnSparseD *_acceptance_7d;//Acceptance
	THnSparseD *_acceptance2_7d;//Acceptance
	Sparse_2d_star _acceptance_5d;//Acceptance
	Sparse_2d_star _acceptance_5d_error;
	THnSparseD * _N;//Experimental Weighted and Corrected Yield
	THnSparseD * _N_pos;//Experimental Weighted and Corrected Yield
	THnSparseD * _N_neg;//Experimental Weighted and Corrected Yield
	Sparse_2d_star _N_5d;//Experimental Weighted and Corrected Yield
	Sparse_2d_star _N_5d_pos;//Experimental Weighted and Corrected Yield
	Sparse_2d_star _N_5d_neg;//Experimental Weighted and Corrected Yield
	THnSparseD * _N_holes;//Experimental Holes
	THnSparseD * _N_holes_pos;//Experimental Holes
	THnSparseD * _N_holes_neg;//Experimental Holes
	THnSparseD * _N2;//Experimental Weighted and Corrected Yield
	THnSparseD * _N2_pos;//Experimental Weighted and Corrected Yield
	THnSparseD * _N2_neg;//Experimental Weighted and Corrected Yield
	Sparse_2d_star _N2_5d;//Experimental Weighted and Corrected Yield
	Sparse_2d_star _N2_5d_pos;//Experimental Weighted and Corrected Yield
	Sparse_2d_star _N2_5d_neg;//Experimental Weighted and Corrected Yield
	THnSparseD * _N2_holes;//Experimental Holes
	THnSparseD * _N2_holes_pos;//Experimental Holes
	THnSparseD * _N2_holes_neg;//Experimental Holes


	THnSparseD * _scale_7d; 
	THnSparseD * _scale_7d_pos; 
	THnSparseD * _scale_7d_neg; 
	THnSparseD * _scale_exp_7d;//The Localized Scale Factor
	THnSparseD * _scale_sim_7d;//The Localized Scale Factor
	Sparse_2d_star _scale_exp_5d;//Localized Scale Factor 3d 
	Sparse_2d_star _scale_sim_5d;//Localized Scale Factor 3d 
	THnSparseD * _scale_exp_7d_pos;//The Localized Scale Factor
	THnSparseD * _scale_sim_7d_pos;//The Localized Scale Factor
	Sparse_2d_star _scale_exp_5d_pos;//Localized Scale Factor 3d 
	Sparse_2d_star _scale_sim_5d_pos;//Localized Scale Factor 3d 
	THnSparseD * _scale_exp_7d_neg;//The Localized Scale Factor
	THnSparseD * _scale_sim_7d_neg;//The Localized Scale Factor
	Sparse_2d_star _scale_exp_5d_neg;//Localized Scale Factor 3d 
	Sparse_2d_star _scale_sim_5d_neg;//Localized Scale Factor 3d 
	double_2d _scale_factor_5d;//5d Integrated Scale Factor

	TH1D_1d_star _X_bin_sizes; //Size of individual bins for non-phi variables {MM1,MM2,theta,alpha}
	TH1D* _phi_bin_sizes;//Width of phi bins 

	TH1F_2d_star _sim_eff_1d;//Simulation efficiency {W,Q2}
	TH1F_2d_star _sim_eff_ratio_1d;//Sim efficiency ratio (){W,Q2}
	TH2F* _sim_eff_wq2;//Sim efficiency as a function of w and Q2 

	TH1D* _acc_zero_exp;//In bins where exp_recon!=0 how many bins of the acceptance are zero
	TH1D_2d_star _acc_yield;//Distribution of acceptance yield {W,Q2}
	TH1D_2d_star _acc_zero_thr;//In bins where thrown!=0 how many bins of the acceptance are zero {W,Q2}
	TH2D* _thr_exp_ratio;//Ratio of thrown to exp_recon as a function of W and Q2
	TH1D_2d_star _acc_rel_error_weighted;//Distribution of Relative error of the Acceptance weighted {W,Q2}
	TH1D_2d_star _acc_rel_error_unweighted;//Distribution of Relative error of the Acceptance unweighted {W,Q2}

	TH2D* _rad_corr;//Radiative Corrections
	double_2d _rad_corr_array;
	double _rad_corr_mu;//integral thrown over integral thrown no-radiative effects
	double_2d _rad_error;//Statistical Error for radiative corrections
	

	//Topology Yields
	double_2d _n_exp_corr;
	double_2d _n_sim_corr;
	double_2d _n_thrown;
	//Sparse Binning

	//W Q2 Histograms
	TH2D* _exp_hist_wq2;
	TH2D* _sim_hist_wq2;
	TH2D* _exp_corr_hist_wq2;
	TH2D* _sim_corr_hist_wq2;
	TH2D* _thr_hist_wq2;

	//Single Differential Histograms
	TH1D_3d_star _single_diff_hist;//{w,q2,X} X->{MM1,MM2,theta,alpha,phi}

	//Polarization Histograms
	TH1D_4d_star _polarization_hist;//{w,q2,X} X->{MM1,MM2,theta,alpha}

	//Acceptance Histograms
	TH1D_3d_star _accept_hist_1;//{w,q2,X} X->{MM1,MM2,theta,alpha,phi}
	TH1D* _accept_hist_2[2];//{w,q2}

	//Beam Spin Histograms
	TH1D_2d_star _beam_spin_hist;

	/*
	//5D THnSparse
	Sparse_2d_star _exp_5d;
	Sparse_2d_star _exp_5d_pos;
	Sparse_2d_star _exp_5d_neg;
	Sparse_2d_star _sim_5d;
	Sparse_2d_star _empty_5d;
	Sparse_2d_star _empty_5d_pos;
	Sparse_2d_star _empty_5d_neg;
	Sparse_2d_star _thrown_5d;
	Sparse_2d_star _norad_5d;
	Sparse_2d_star _acceptance_5d;
	Sparse_2d_star _sim_holes_5d;
	Sparse_2d_star _exp_holes_5d;
	Sparse_2d_star _exp_holes_5d_pos;
	Sparse_2d_star _exp_holes_5d_neg;
	Sparse_2d_star _scale_sim_5d;
	Sparse_2d_star _scale_exp_5d;
	Sparse_2d_star _scale_sim_5d_pos;
	Sparse_2d_star _scale_exp_5d_pos;
	Sparse_2d_star _scale_sim_5d_neg;
	Sparse_2d_star _scale_exp_5d_neg;
	Sparse_2d_star _N_5d;
	Sparse_2d_star _N_5d_pos;
	Sparse_2d_star _N_5d_neg; 

	//1D Phi Histograms
	TH1D_4d_star _exp_1d;
	TH1D_4d_star _exp_1d_pos;
	TH1D_4d_star _exp_1d_neg;
	TH1D_4d_star _sim_1d;
	TH1D_4d_star _empty_1d;
	TH1D_4d_star _empty_1d_pos;
	TH1D_4d_star _empty_1d_neg;
	TH1D_4d_star _thrown_1d;
	TH1D_4d_star _acceptance_1d;
	TH1D_4d_star _sim_holes_1d;
	TH1D_4d_star _exp_holes_1d;
	TH1D_4d_star _exp_holes_1d_pos;
	TH1D_4d_star _exp_holes_1d_neg;
	TH1D_4d_star _scale_sim_1d;
	TH1D_4d_star _scale_exp_1d;
	TH1D_4d_star _N_1d;
	TH1D_4d_star _N_1d_pos;
	TH1D_4d_star _N_1d_neg;
	*/
	
public:
	Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_);
	Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_);
	void Make_Histograms(Flags flags_);
	void Fill_Histograms(Flags flags_);
	void Write_Histograms(Flags flags_);
	void Name_Output(Flags flags_);
	//std::shared_ptr<TFile> Name_Output(Flags flags_);
	//void Extract_Event_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_);
	//void Calculate_5d_Holes(Flags flags_);
	//void Get_N(Flags flags_);
	//void Single_Differential(Flags flags_);
	//void Polarization_Observables(Flags flags_);//Making Polarization Histograms
	//void Beam_Spin(Flags flags_);

	void Extract_7d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_);//(TFile* exp_tree, TFile* sim_tree, Flags flags_);
	void Extract_7d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_);
	void Acceptance();
	void Rad_Corr();
	void Make_N_7d(Flags flags_);
	void Localized_Holes(Flags flags_, int min_dist_, int max_dist_);
	void Sparse_7to5(Flags flags_);//Convert 7d histograms to 5d histograms
	void Single_Differential(Flags flags_);//Making Single Differential Histograms
	void Polarization_Observables(Flags flags_);//Making Polarization Histograms
	void Beam_Spin(Flags flags_);
	void Make_WQ2(Flags flags_);
	void Dumb_Histograms(Flags flags_);
	void Extract_Bin_Info(Flags flags_);//Extract binning information for 7 and 5d histograms
	void Skeleton_5D(Flags flags_);//Create Empty 5D THnSparse to fill
	
	float Bin_Size( int variable, int bin_7d);
	float W_Bin_Size( int bin_7d);
	float Q2_Bin_Size( int bin_7d);
	float X_Bin_Size( int variable, int bin_7d);
	float Phi_Bin_Size( int bin_7d);
	double Bin_Center_7d( int dim, int bin);
	void XandPhi_BinHistograms(Flags flags_);
	void Acceptance_Errors(Flags flags_);
	void Make_Error_Hists(Flags flags_);
	void Fill_Error_Hists(Flags flags_);
	void Write_Error_Hists(Flags flags_);
	void Acceptance_Histograms(Flags flags_);
	void Hole_Histograms(Flags flags_);
	/*int W_bins();
	int W_bin(float W_);
	float W_bot(int bin_);
	float W_top(int bin_);
	float Q2_bot(int bin_);
	float Q2_top(int bin_);
	int Q2_bin(float Q2_);
	int Xij_bin(int var_, int var_set_, float Xij_);
	float W_center(int bin_);
	float Q2_center(int bin_);
	float MM_center(int bin_);
	float MM2_center(int bin_);
	float Theta_center(int bin_);
	float Alpha_center(int bin_);
	float Phi_center(int bin_);
	double W_bin_size(int bin_);
	double Q2_bin_size(int bin_);
	double Xij_size(int Xij_idx_, int bin_);
	double Xij_bot(int Xij_idx_, int bin_);
	double Xij_top(int Xij_idx_, int bin_);
	*/
	void Print_Histogram_Bin_Info(THnSparseD* hist_);
	void Calc_Error_R(Flags flags_);
	//void Calc_Error_A(Flags flags_);
	//void Calc_Error_Single_Diff(Flags flags_);
	//void Calc_Error_Polarization(Flags flags_);
	//void Calc_Error_Beam_Spin(Flags flags_);
};

#endif
