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


class Histogram{
protected:
	//RootFiles
	std::shared_ptr<TFile> _RootOutputFile;

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
	THnSparseD *_exp_data_7d_pos;//Experimental reconstruction 7d with postitive helicity
	THnSparseD *_exp_data_7d_neg;//Experimental reconstruction 7d with negative helicity
	THnSparseD *_exp_data_7d;//Experimental Reconstruction 7 dimension
	THnSparseD *_sim_data_7d;//Simulated Reconstruction
	THnSparseD *_thrown_7d_no_rad;//Simulated thrown with no radiative effects
	THnSparseD *_sim_weight_sq_7d;//summed square of weights for each 7d bin for simulation
	THnSparseD *_empty_7d;//Empty Reconstruction
	THnSparseD *_empty_7d_pos;//Empty Reconstruction positive helicity
	THnSparseD *_empty_7d_neg;//Empty Reconstruction negative helicity
	
	Sparse_2d_star _exp_data_5d;//Exp Recon 5 Dimension {W,Q2}
	Sparse_2d_star _empty_5d;//Exp Empty Recon 5d {W,Q2}
	Sparse_2d_star _sim_data_5d;//Sim Recon 5 Dimension {W,Q2}
	THnSparseD *_thrown_7d;//Thrown Simulation
	Sparse_2d_star _thrown_5d;//Thrown Sim 5d {W,Q2}
	Sparse_2d_star _sim_holes_5d;//Sim Recon Holes {W,Q2}
	Sparse_2d_star _exp_holes_5d;//Exp Recon Holes {W,Q2}
	Sparse_2d_star _exp_corr_5d;//Exp Acceptance Corrected {W,Q2}
	Sparse_2d_star _empty_corr_5d;//Exp Acceptance Corrected {W,Q2}
	Sparse_2d_star _sim_corr_5d;//Sim Recon Acceptance Corrected  {W,Q2}
	Sparse_2d_star _exp_corr_holes_5d;//Experimental Accept Corrected and Hole Filled {W,Q2}
	double_2d _scale_factor_5d;//5d scale factor of simulation to experimental yields {W,Q2}
	double _scale_factor_7d;//7d scale factor of simulation to experimental yields 
	THnSparseD *_exp_holes_7d;//Exp holes
	THnSparseD *_sim_holes_7d;//Sim Holes
	THnSparseD *_exp_corr_7d;//Exp accept corrected
	THnSparseD *_empty_corr_7d;//Exp accept corrected
	THnSparseD *_sim_corr_7d;//Sim accept corrected
	THnSparseD *_cross_section_7d;//Cross Section 
	Sparse_2d_star _cross_section_5d;//;
	Sparse_2d_star _acceptance_5d;//;
	TH2D_3d_star _cross_section_2d;//
	THnSparseD *_acceptance_7d;
	TH1D_3d_star _exp_corr_holes_3d; //For single differential {W,Q2,Xij}[Xij bins]
	TH1D_4d_star _exp_corr_holes_4d; //For Polarization Observables {W,Q2,Xij,Xij_bin}[phi bins]
	THnSparseD *_acceptance_eff_7d; // Acceptance Efficiency? 
	THnSparseD *_acceptance_err_7d[2];//Last one is Weighted or not
	 

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

	//std::vector<std::vector<double>> _rad_corr;//Radiative Correction
	//std::vector<std::vector<double>> _rad_corr2;//Radiative Correction2

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
	TH1D_3d_star _polarization_hist;//{w,q2,X} X->{MM1,MM2,theta,alpha}

	//Acceptance Histograms
	TH1D_3d_star _accept_hist_1;//{w,q2,X} X->{MM1,MM2,theta,alpha,phi}
	TH1D* _accept_hist_2[2];//{w,q2}

	//Beam Spin Histograms
	TH1D_2d_star _beam_spin_hist;
	
public:
	Histogram(const std::string& output_file_, TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_);
	//Histogram(const std::string& output_file, TFile* exp_tree, TFile* sim_tree, TFile* empty_tree, Flags flags_);
	//Histogram(const std::string& output_file, TFile* exp_tree, TFile* sim_tree, TFile* exp_tree2, TFile* sim_tree2, Flags flags_);
	//Histogram(const std::string& output_file, TFile* exp_tree, TFile* sim_tree, TFile* empty_tree, TFile* exp_tree2, TFile* sim_tree2, TFile* empty_tree, Flags flags_);
	void Make_Histograms(Flags flags_);
	void Fill_Histograms(Flags flags_);
	void Write_Histograms(Flags flags_);
	std::shared_ptr<TFile> Name_Output(Flags flags_);
	void Extract_7d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_);//(TFile* exp_tree, TFile* sim_tree, Flags flags_);
	void Sparse_Add_7d(THnSparseD &h0, THnSparseD* h1, THnSparseD* h2, int sign);//Add/Subtract Sparse Histograms
	void Sparse_Add_5d(THnSparseD* &h0, THnSparseD* h1, THnSparseD* h2, int sign);//Add/Subtract Sparse Histograms
	void Sparse_7to5(Flags flags_);//Convert 7d histograms to 5d histograms
	void Sparse_5to3(Flags flags_);//Making Single Differential Histograms
	void Sparse_5to4(Flags flags_);//Making Polarization Histograms
	void Extract_Bin_Info(Flags flags_);//Extract binning information for 7 and 5d histograms
	void Skeleton_5D(Flags flags_);//Create Empty 5D THnSparse to fill
	void Calc_Acceptance(Flags flags_);//Calculate acceptance from 5d histograms
	void Calc_Holes_Sim(Flags flags_);//Calculate Estimated Holes in Simulation
	void Calc_Holes_Exp(Flags flags_);//Calculate Estimated Holes in Experiment
	void Calc_Holes(Flags flags_);//Calculate Estimated Holes for Sim and Exp
	void Calc_Cross_Section(Flags flags_);//Calculate Cross Section
	void Make_Single_Diff(Flags flags_);//Make Single Differential Cross Section Histograms
	void Make_Polarization(Flags flags_);//Make Histograms to extract polarization observables
	void Beam_Spin(Flags flags_);
	void Make_Integrated(Flags flags_);//Make Integrated Cross Section Histograms
	void Make_WQ2(Flags flags_);//Make WQ2 histograms to show binnning
	//void Make_Acceptance_Statistics(Flags flags_);//Make Histograms for determining proper Acceptance Statistics //Cannot make this here in the form the data has already been placed. Need on rootfile level for sim recon
	void Make_Acceptance(Flags flags_);
	//void Write_5d_Yield();
	//void Write_5d_Cross_Section();
	//void Write_5d_Holes();
	//void Write_7d_Holes();
	//void Write_Single_Diff();
	//void Write_Polarization();
	//void Write_Integrated();
	void Convert_to_Cross(Flags flags_);
	void Convert_Single_Diff_to_Cross(Flags flags_);
	void Convert_Polarization_to_Cross(Flags flags_);
	void Write_WQ2(Flags flags_);
	void Write_Acceptance(Flags flags_);
	void Write_Single_Diff(Flags flags_);
	void Write_Polarization(Flags flags_);
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
	THnSparseD* Add_Sparse(THnSparse * h1_, THnSparse * h2_);
	//void Radiative_Correction(Flags flags_);
	//*-------------------------------Start Hole Identification-----------------------*
	int Hole_Size(THnSparse* hist_); 
	
};

#endif
