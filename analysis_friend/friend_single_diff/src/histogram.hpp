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

//static const float _mp_ = 0.93828;	//Mass of proton in GeV
//static const float _mpi_ = 0.1395; //Mass of Pion in GeV

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
    THnSparseD *_sim_data_7d;//Weighted
	THnSparseD *_empty_7d;
	THnSparseD *_thrown_7d;//Weighted
	THnSparseD *_thrown_7d_no_rad;//No radiative effects, weighted
	THnSparseD * _sim_holes_7d;//Simulation Holes
	THnSparseD * _sim_holes_tmp_7d;
    THnSparseD * _N_holes;//Experimental Holes
	THnSparseD *_acceptance_7d;//Acceptance
	THnSparseD * _N;//Experimental Weighted and Corrected Yield
	Sparse_2d_star _N_5d;//Experimental Weighted and Corrected Yield

	TH1D_1d_star _X_bin_sizes; //Size of individual bins for non-phi variables {MM1,MM2,theta,alpha}
	TH1D* _phi_bin_sizes;//Width of phi bins 

	TH2D* _rad_corr;//Radiative Corrections
	double_2d _rad_corr_array;
	double _rad_corr_mu;//integral thrown over integral thrown no-radiative effects
	double_2d _rad_error;//Statistical Error for radiative corrections
	

	//Topology Yields
	double_2d _n_exp_corr;
	double_2d _n_sim_corr;
	double_2d _n_thrown;

	//Single Differential Histograms
	TH1D_3d_star _single_diff_hist;//{w,q2,X} X->{MM1,MM2,theta,alpha,phi}

	
public:
	Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_);
    Histogram(TFile* exp_tree_, TFile* sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_);
	void Extract_7d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, TFile *holes_, Flags flags_);
    void Extract_7d_Histograms(TFile *exp_tree_, TFile *sim_tree_, TFile *empty_tree_, TFile *nr_sim_tree_, Flags flags_);
	void Extract_Bin_Info(Flags flags_);//Extract binning information for 7 and 5d histograms
	void Sparse_7to5(Flags flags_);//Convert 7d histograms to 5d histograms
    void Rad_Corr();
	void Single_Diff(Flags flags_);//Making Single Differential Histograms
	

	
};

#endif

