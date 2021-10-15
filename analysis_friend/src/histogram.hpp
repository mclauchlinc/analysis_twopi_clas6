#ifndef HISTOGRAM_H_GUARD
#define HISTOGRAM_H_GUARD

#include "TH1.h"
#include "TH2.h"
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


class Histogram{
protected:
	//RootFiles
	std::shared_ptr<TFile> _RootOutputFile;

	static const int _n_var_sets = 3;
	static const int _n_topology = 5;

	//Sparse Binning
	std::vector<int> _n_bins_7d;//Number of bins in a given dimension //Difficult to define inside a class like this
	std::vector<std::vector<int>> _n_bins_5d; //var set, then 5dim variables
	double_3d _bin_low_7d;//Low edge of a bin  [var_set,variable,bin # inside variable]
	double_3d _bin_up_7d;//Top Edge of a bin 
	double_2d _bin_lowside_7d; //low values of the bins
	double_2d _bin_topside_7d; //top values of the bins
	double_3d _bin_mid_7d;//Middle value of a bin
	double_3d _bin_size_7d;//Width of a bin
	double_3d _bin_edges_7d;
	double_3d _bin_low_5d;//Low edge of a bin [var,bins]
	double_3d _bin_up_5d;//Top Edge of a bin
	double_3d _bin_mid_5d;//Middle value of a bin
	double_3d _bin_size_5d;//Width of a bin
	double_3d _bin_edges_5d;//Edges of a bin
	double_2d _bin_lowside_5d;
	double_2d _bin_topside_5d;


	//Sparse Histograms
	THnSparseD *_exp_data_7d[_n_var_sets][_n_topology];
	THnSparseD *_sim_data_7d[_n_var_sets][_n_topology];
	THnSparseD *_sim_weight_sq_7d[_n_var_sets][_n_topology];
	Sparse_4d_star _exp_data_5d;//[_n_var_sets][_n_topology];//would have binning for W and then Q2
	Sparse_4d_star _sim_data_5d;//[_n_var_sets][_n_topology];//would have binning for W and then Q2
	THnSparseD *_thrown_7d[_n_var_sets];//would have binning for W and then Q2
	Sparse_3d_star _thrown_5d;//[_n_var_sets];//would have binning for W and then Q2
	Sparse_4d_star _sim_holes_5d;//[_n_var_sets][_n_topology];//would have binning for W and then Q2
	Sparse_4d_star _exp_holes_5d;//[_n_var_sets][_n_topology];//would have binning for W and then Q2
	Sparse_4d_star _exp_corr_5d;//[_n_var_sets][_n_topology];//would have binning for W and then Q2
	Sparse_4d_star _sim_corr_5d;//[_n_var_sets][_n_topology];//would have binning for W and then Q2//needed due to holes present in reconstructed not in thrown
	Sparse_4d_star _exp_corr_holes_5d;
	double_4d _scale_factor_5d;//[_n_var_sets][_n_topology]; //would have binning for W and then Q2= NAN;
	double _scale_factor_7d[_n_var_sets][_n_topology];
	THnSparseD *_exp_holes_7d[_n_var_sets][_n_topology];
	THnSparseD *_sim_holes_7d[_n_var_sets][_n_topology];
	THnSparseD *_exp_corr_7d[_n_var_sets][_n_topology];
	THnSparseD *_sim_corr_7d[_n_var_sets][_n_topology];
	THnSparseD *_cross_section_7d[_n_var_sets][_n_topology];
	Sparse_4d_star _cross_section_5d;//[_n_var_sets][_n_topology];
	Sparse_4d_star _acceptance_5d;//[_n_var_sets][_n_topology];
	THnSparseD *_acceptance_7d[_n_var_sets][_n_topology];
	TH1D_5d_star _exp_corr_holes_3d; //For single differential {var,top,W,Q2,Xij}[Xij bins]
	TH1D_6d_star _exp_corr_holes_4d; //For Polarization Observables {var,top,W,Q2,Xij,Xij_bin}[phi bins]
	THnSparseD *_acceptance_eff_7d[_n_var_sets][_n_topology];
	THnSparseD *_acceptance_err_7d[_n_var_sets][_n_topology][2];//Last one is Weighted or not

	TH1D_2d_star _X_bin_sizes; //var set, Xij
	TH1D_1d_star _phi_bin_sizes;

	TH1F_3d_star _sim_eff_1d;//{top,W,Q2}
	TH1F_3d_star _sim_eff_ratio_1d;//{top,W,Q2}
	TH2F_1d_star _sim_eff_wq2;//{top}

	TH1D_2d_star _acc_zero_exp;//{var,top}
	TH1D_4d_star _acc_yield;//{var,top,W,Q2}
	TH1D_4d_star _acc_zero_thr;//{var,top,W,Q2}
	TH2D_2d_star _thr_exp_ratio;//{var,top}
	TH1D_4d_star _acc_rel_error_weighted;//{var,top,W,Q2}
	TH1D_4d_star _acc_rel_error_unweighted;//{var,top,W,Q2}



	//Topology Yields
	double_4d _n_exp_corr;
	double_4d _n_sim_corr;
	double_3d _n_thrown;
	//Sparse Binning

	//W Q2 Histograms
	TH2D *_exp_hist_wq2[_n_var_sets][_n_topology];
	TH2D *_sim_hist_wq2[_n_var_sets][_n_topology];
	TH2D *_exp_corr_hist_wq2[_n_var_sets][_n_topology];
	TH2D *_sim_corr_hist_wq2[_n_var_sets][_n_topology];
	TH2D *_thr_hist_wq2[_n_topology];

	//Acceptance Histograms
	TH1D_5d_star _accept_hist_1;//(var,top,w,q2,{MM1,MM2,theta,alpha,phi})
	TH1D *_accept_hist_2[_n_var_sets][_n_topology][2];//(var,top,{w,q2})
	
public:
	Histogram(const std::string& output_file, TFile* exp_tree, TFile* sim_tree, Flags flags_);
	//Histogram(const std::string& output_file, TFile* exp_tree, TFile* sim_tree, TFile* empty_tree, Flags flags_);
	//Histogram(const std::string& output_file, TFile* exp_tree, TFile* sim_tree, TFile* exp_tree2, TFile* sim_tree2, Flags flags_);
	//Histogram(const std::string& output_file, TFile* exp_tree, TFile* sim_tree, TFile* empty_tree, TFile* exp_tree2, TFile* sim_tree2, TFile* empty_tree, Flags flags_);
	void Make_Histograms();
	void Fill_Histograms();
	void Write_Histograms();
	void Extract_7d_Histograms(TFile* exp_tree, TFile* sim_tree);
	void Sparse_Add_7d(THnSparseD &h0, THnSparseD* h1, THnSparseD* h2, int sign, int var);//Add/Subtract Sparse Histograms
	void Sparse_Add_5d(THnSparseD* &h0, THnSparseD* h1, THnSparseD* h2, int sign, int var);//Add/Subtract Sparse Histograms
	void Sparse_7to5(Flags flags_);//Convert 7d histograms to 5d histograms
	void Sparse_5to3(Flags flags_);//For Single Differential bins
	void Sparse_5to4(Flags flags_);//For Polarization Observables
	void Extract_Bin_Info();//Extract binning information for 7 and 5d histograms
	void Skeleton_5D();//Create Empty 5D THnSparse to fill
	void Calc_Acceptance();//Calculate acceptance from 5d histograms
	void Calc_Holes_Sim();//Calculate Estimated Holes in Simulation
	void Calc_Holes_Exp();//Calculate Estimated Holes in Experiment
	void Calc_Holes();//Calculate Estimated Holes for Sim and Exp
	void Calc_Cross_Section();//Calculate Cross Section
	void Make_Single_Diff();//Make Single Differential Cross Section Histograms
	void Make_Polarization();//Make Histograms to extract polarization observables
	void Make_Integrated();//Make Integrated Cross Section Histograms
	void Make_WQ2();//Make WQ2 histograms to show binnning
	void Make_Acceptance();
	//void Write_5d_Yield();
	//void Write_5d_Cross_Section();
	//void Write_5d_Holes();
	//void Write_7d_Holes();
	//void Write_Single_Diff();
	//void Write_Polarization();
	//void Write_Integrated();
	void Convert_to_Cross();
	void Convert_Single_Diff_to_Cross();
	void Convert_Polarization_to_Cross();
	void Write_WQ2();
	void Write_Acceptance();
	void Write_Single_Diff();
	void Write_Polarization();
	float Bin_Size(int var_set, int variable, int bin_7d);
	float W_Bin_Size(int var_set, int bin_7d);
	float Q2_Bin_Size(int var_set, int bin_7d);
	float X_Bin_Size(int var_set, int variable, int bin_7d);
	float Phi_Bin_Size(int var_set, int bin_7d);
	double Bin_Center_7d(int var_set, int dim, int bin);
	void XandPhi_BinHistograms();
	void Acceptance_Errors();
	void Make_Error_Hists();
	void Fill_Error_Hists();
	void Write_Error_Hists();
};

#endif
