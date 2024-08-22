#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "constants.hpp"
#include "functions.hpp"
#include "CartesianGenerator.hpp"
#include "physics.hpp"
#include "detectors.hpp"
#include "THnSparse.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include "TGraph.h"
#include <iterator>
#include "cuts.hpp"
#include <cmath>
#include "TThread.h"
#include "TLatex.h"
#include <mutex>//For pausing the multithreading for the filling of THnSparse histograms 
//#include <unistd.h>//to allow us to use chdir
//#include "TImage.h"
//#include "particle.hpp"
//#include "variables.h"
//#include "CartesianGenerator.hh"

//Binning Info
//W Q2
static float _wq2_xmin_ = -0.01;//GeV
static float _wq2_xmax_ = 3.99; //GeV
static int _wq2_xbin_ = 200;//GeV
static float _wq2_ymin_ = -0.01;//GeV^2
static float _wq2_ymax_ = 8.99; //GeV^2
static int _wq2_ybin_ = 200;//GeV^2

//PID Plots
//Fiducial 
static float _fid_xmin_ = -30.0;//Degrees
static float _fid_xmax_ = 30.0; //Degrees
static int _fid_xbin_ = 400;//Degrees
static float _fid_ymin_ = 0.0;//Degrees
static float _fid_ymax_ = 180.0; //Degrees
static int _fid_ybin_ = 300;//Degrees
static float _geo_fid_xmin_[3] = {-250.0,-250.0,-200.0};//Made for detector centering {-600.0,0.0,-400.0};//mm
static float _geo_fid_xmax_[3] = {250.0,250.0,200.0};//Made for detector centering {600.0,400.0,400.0}; //mm
static int _geo_fid_xbin_[3] = {200,200,200};//cm
static float _geo_fid_ymin_[3] = {50.0,50.0,50.0};//Made for detector centering{-600.0,-250.0,-400.0};//mm
static float _geo_fid_ymax_[3] = {500.0,500.0,500.0};//Made for detector centering{600.0,250.0,400.0}; //mm
static int _geo_fid_ybin_[3] = {200,250,250};//cm
static int _geo_fid_cc_segments_ = 18;
static int _geo_fid_cc_sides_ = 3;
static int _geo_fid_sc_paddles_ = 48;
static float _geo_fid2_xmin_[3] = {-600.0,-500.0,-500.0};
static float _geo_fid2_xmax_[3] = {600.0,500.0,500.0};
static int _geo_fid2_xbin_[3] = {600,600,600};
static float _geo_fid2_ymin_[3] = {-600.0,-500.0,-500.0};
static float _geo_fid2_ymax_[3] = {600.0,500.0,500.0};
static int _geo_fid2_ybin_[3] = {600,600,600};
//Delta T
static float _delta_xmin_ = 0.0;//GeV
static float _delta_xmax_ = 6.0; //GeV
static int _delta_xbin_ = 1000;//GeV
static float _delta_ymin_ = -6.0;//ns
static float _delta_ymax_ = 6.0; //ns
static int _delta_ybin_ = 400;//ns
//SF
static float _sf_xmin_ = 0.0;//GeV
static float _sf_xmax_ = 6.0; //GeV
static int _sf_xbin_ = 500;//GeV
static float _sf_ymin_ = 0;//unitless
static float _sf_ymax_ = 1.0; //unitless
static int _sf_ybin_ = 300;//unitless
//CC
//nphe
static float _cc_xmin_ = -0.5;//Segment
static float _cc_xmax_ = 500.5; //Segment
static int _cc_xbin_ = 501;//Segment
//EC
static float _ec_xmin_ = 0.0;//GeV
static float _ec_xmax_ = 6.0; //GeV
static int _ec_xbin_ = 100;//GeV
//MM
static float _mm_min_[4] = {0.836,-0.158,-0.158,-0.045};
static float _mm_max_[4] = {1.2,0.273,0.273,0.045};
static int _mm_bin_[4] = {200,200,200,200};
//MM2
static float _mm2_min_[4] = {0.7,-0.025,-0.025,-0.002};
static float _mm2_max_[4] = {1.5,0.075,0.075,0.002};
static int _mm2_bin_[4] = {100,100,100,100};
//EC
static float _ec_min_ = 0.0;
static float _ec_max_ = 6.0;
static int _ec_bin_ = 200; 

//Detector Plots
//Vertex
static float _vertex_min_[2] = {-10.0,-35.0};
static float _vertex_max_[2] = {10.0,2.0};
static int _vertex_bin_[2] = {300,400};
//CC Eff
	//Momentum vs. Theta
static float _cc_eff1_xmin_ = 0.0;//GeV
static float _cc_eff1_xmax_ = 5.0;//GeV
static int _cc_eff1_xbin_ = 300;
static float _cc_eff1_ymin_ = 0.0;//Degrees
static float _cc_eff1_ymax_ = 120.0;//Degrees
static int _cc_eff1_ybin_ = 300;
	//Sector and Segment yields
static float _cc_eff2_xmin_ = 0.5;//Sector
static float _cc_eff2_xmax_ = 6.5; //Sector
static int _cc_eff2_xbin_ = 6;//Sector

//SC Eff
static float _sc_eff_xmin_ = -1.5;//Sector
static float _sc_eff_xmax_ = 58.5; //Sector
static int _sc_eff_xbin_ = 60;//Sector

//For Momentum Binning
static float _p_bin_min_ = 0.5;//GeV
static float _p_bin_max_ = 5.0;//GeV
static int _p_bin_bins_= 25;//Steps

//Detector Parsing
static int _n_cc_seg_ = 18;
static int _n_cc_lrc_ = 3; 
static int _n_sec_ = 6; 


//THnSparse Binning
static int _n_var_ = 3; //Number of variable sets
static float _W_min_ = 1.4;
static float _W_max_ = 2.125;
static float _W_res_ = 0.025;
static float _Q2_min_ = 2.0;
static float _Q2_max_ = 5.0;
static float _Q2_bins_[6] = {2.0,2.4,3.0,3.5,4.2,5.0};
//Changed to these arrays to match Arjun's binning 
static float _MM_min_[3] = {1.07784,0.27914,1.07784};//{1.077835,0.309141,1.077835};//changed 7-12-23//{1.1,0.3,1.1};//Changed 5/5/23 //Changed 8-28-23 to match Arjun
//static float _MM_max_[3] = {_W_max_-_mpi_,_W_max_-_mp_,_W_max_-_mpi_};//{2.0,1.1,2.0};//Changed 5/5/23
static float _MM_max_[3][29] = {	{1.28543, 1.43543, 1.43543, 1.43543, 1.43543, 1.43543, 1.43543, 1.58543, 1.58543, 1.58543, 1.58543, 1.58543, 1.58543, 1.71043, 1.71043, 1.71043, 1.71043, 1.71043, 1.86043, 1.86043, 1.86043, 1.86043, 1.86043, 1.86043, 1.98543, 1.98543, 1.98543, 1.98543, 1.98543}, 
									{0.486728, 0.636728, 0.636728, 0.636728, 0.636728, 0.636728, 0.636728, 0.786728, 0.786728, 0.786728, 0.786728, 0.786728, 0.786728, 0.911728, 0.911728, 0.911728, 0.911728, 0.911728, 1.06173, 1.06173, 1.06173, 1.06173, 1.06173, 1.06173, 1.18673, 1.18673, 1.18673, 1.18673, 1.18673}, 
									{1.28543, 1.43543, 1.43543, 1.43543, 1.43543, 1.43543, 1.43543, 1.58543, 1.58543, 1.58543, 1.58543, 1.58543, 1.58543, 1.71043, 1.71043, 1.71043, 1.71043, 1.71043, 1.86043, 1.86043, 1.86043, 1.86043, 1.86043, 1.86043, 1.98543, 1.98543, 1.98543, 1.98543, 1.98543}};
static int _MM_bins_ = 14;//8/28/23 to Arjun's actual binning;
static int _MM_wider_ = 5; 
static float _MM2_min_[3] = {_MM_min_[1],_MM_min_[2],_MM_min_[0]};//changed 7-12-23//{0.3,1.1,1.1};//Changed 5/5/23//Changed 8-28-23 to match Arjun
//static float _MM2_max_[3][29] = {_MM_max_[1],_MM_max_[2],_MM_max_[0]};//{1.1,2.0,2.0};//Changed 5/5/23
static float _MM2_max_[3][29] = {	{0.486728, 0.636728, 0.636728, 0.636728, 0.636728, 0.636728, 0.636728, 0.786728, 0.786728, 0.786728, 0.786728, 0.786728, 0.786728, 0.911728, 0.911728, 0.911728, 0.911728, 0.911728, 1.06173, 1.06173, 1.06173, 1.06173, 1.06173, 1.06173, 1.18673, 1.18673, 1.18673, 1.18673, 1.18673},
									{1.28543, 1.43543, 1.43543, 1.43543, 1.43543, 1.43543, 1.43543, 1.58543, 1.58543, 1.58543, 1.58543, 1.58543, 1.58543, 1.71043, 1.71043, 1.71043, 1.71043, 1.71043, 1.86043, 1.86043, 1.86043, 1.86043, 1.86043, 1.86043, 1.98543, 1.98543, 1.98543, 1.98543, 1.98543},  
									{1.28543, 1.43543, 1.43543, 1.43543, 1.43543, 1.43543, 1.43543, 1.58543, 1.58543, 1.58543, 1.58543, 1.58543, 1.58543, 1.71043, 1.71043, 1.71043, 1.71043, 1.71043, 1.86043, 1.86043, 1.86043, 1.86043, 1.86043, 1.86043, 1.98543, 1.98543, 1.98543, 1.98543, 1.98543}};
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
static float _ppi_offset_ = 0.012565; //changed 7-12-23//GeV offset from threshold MM value of missing pion
static float _pipi_offset_ = 0.0324895; //changed 7-12-23//GeV offset from threshold MM value of missing proton
static float _MM_offset_[3] = {_mpi_+_ppi_offset_,_mp_+_pipi_offset_,_mpi_+_ppi_offset_};
static float _MM2_offset_[3] = {_MM_offset_[1],_MM_offset_[2],_MM_offset_[0]};

static int _top_check_bins_ = 16;
static int _top_check_min_ = 0;
static int _top_check_max_ = 16;


//static float _bin_low_ = {_W_min_,_Q2_min_,_MM}


//using TH2F_ptr = std::make_shared<TH2F*>;
//using TH1F_ptr = std::make_shared<TH1F*>;
//using THn_ptr = std::make_shared<THnSparseD*>;
//using TGraph_ptr = std::make_shared<TGraph*>;
//using TDir_ptr = std::make_shared<TDirectory*>;


using Bool_1d = std::vector<bool>;
using Bool_2d = std::vector<std::vector<bool>>;
using Bool_3d = std::vector<std::vector<std::vector<bool>>>;
using Bool_4d = std::vector<std::vector<std::vector<std::vector<bool>>>>;
using Bool_5d = std::vector<std::vector<std::vector<std::vector<std::vector<bool>>>>>;
using Bool_6d = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<bool>>>>>>;



using TH2F_ptr = std::shared_ptr<TH2F>;
using TH2F_ptr_1d = std::vector<TH2F*>;
using TH2F_ptr_2d = std::vector<std::vector<TH2F*>>;
using TH2F_ptr_3d = std::vector<std::vector<std::vector<TH2F*>>>;
using TH2F_ptr_4d = std::vector<std::vector<std::vector<std::vector<TH2F*>>>>;
using TH2F_ptr_5d = std::vector<std::vector<std::vector<std::vector<std::vector<TH2F*>>>>>;
using TH2F_ptr_6d = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<TH2F*>>>>>>;
using TH2F_ptr_7d = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<TH2F*>>>>>>>;

using TH1F_ptr = std::shared_ptr<TH1F>;
using TH1F_ptr_1d = std::vector<TH1F*>;
using TH1F_ptr_2d = std::vector<std::vector<TH1F*>>;
using TH1F_ptr_3d = std::vector<std::vector<std::vector<TH1F*>>>;
using TH1F_ptr_4d = std::vector<std::vector<std::vector<std::vector<TH1F*>>>>;
using TH1F_ptr_5d = std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*>>>>>;
using TH1F_ptr_6d = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*>>>>>>;

using TH1D_ptr = std::shared_ptr<TH1D>;
using TH1D_ptr_1d = std::vector<TH1D*>;
using TH1D_ptr_2d = std::vector<std::vector<TH1D*>>;
using TH1D_ptr_3d = std::vector<std::vector<std::vector<TH1D*>>>;
using TH1D_ptr_4d = std::vector<std::vector<std::vector<std::vector<TH1D*>>>>;
using TH1D_ptr_5d = std::vector<std::vector<std::vector<std::vector<std::vector<TH1D*>>>>>;
using TH1D_ptr_6d = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<TH1D*>>>>>>;
using TH1D_ptr_7d = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<TH1D*>>>>>>>;

using TDir_ptr_1d = std::vector<TDirectory*>;
using TDir_ptr_2d = std::vector<std::vector<TDirectory*>>;
using TDir_ptr_3d = std::vector<std::vector<std::vector<TDirectory*>>>;
using TDir_ptr_4d = std::vector<std::vector<std::vector<std::vector<TDirectory*>>>>;
using TDir_ptr_5d = std::vector<std::vector<std::vector<std::vector<std::vector<TDirectory*>>>>>;

using Sparse_ptr = std::shared_ptr<THnSparseD*>;
using Sparse_ptr_1d = std::vector<THnSparseD*>;
using Sparse_ptr_2d = std::vector<std::vector<THnSparse*>>;

//Canvas information
//const Double_t WQ2_cw = 1200;
//const Double_t WQ2_ch = 600;
//using TH1I_ptr = std::shared_ptr<TH1I>;
//Lengths of different arrays
const static int _len_ecuts_ = std::distance(std::begin(_ecuts_), std::end(_ecuts_));
const static int _len_top_ = std::distance(std::begin(_top_), std::end(_top_));
const static int _len_cut_ = std::distance(std::begin(_cut_), std::end(_cut_));
const static int _len_recon_ = std::distance(std::begin(_recon_), std::end(_recon_));
const static int _len_weight_ = std::distance(std::begin(_weight_), std::end(_weight_));
const static int _len_species_ = std::distance(std::begin(_species_), std::end(_species_));

const static int _wider_ = 5;//The additional bins on either side of Friend Invariant Masses


class Histogram {
protected:
	
	std::shared_ptr<TFile> _RootOutputFile;
	std::shared_ptr<TFile> _SparseFile;
	std::shared_ptr<TFile> _HistImageFile;

	TCanvas* def;

	//Plot Formation Constants
	//W Q2
	 int WQxres = 200;
	 int WQyres = 200;
	 double WQxmin = -0.01;
	 double WQymin = -0.01;
	 double WQxmax = 3.99;
	 double WQymax = 8.99;

	//Electron Sampling Fraction
	 int SFxres = 500;
	 int SFyres = 300;
	 double SFxmin = 0.0;
	 double SFymin = 0.0;
	 double SFxmax = 6.0;
	 double SFymax = 1.0;
	//Minimum Electron Energy
	 int MEExres = 300;
	 int MEEyres = 300;
	 double MEExmin = 0.0;
	 double MEEymin = 0.0;
	 double MEExmax = 6.0;
	 double MEEymax = 1.0;
	//Fiducial Cuts
	 int FIDxres = 400;
	 int FIDyres = 300;
	 double FIDxmin = -30.0;
	 double FIDymin = 0.0;
	 double FIDxmax = 30.0;
	 double FIDymax = 180.0;
	 int GEOFIDxres = 200;
	 int GEOFIDyres = 200;
	 double GEOFIDxmin = -30.0;
	 double GEOFIDymin = 0.0;
	 double GEOFIDxmax = 30.0;
	 double GEOFIDymax = 180.0;
	//Delta_t
	 int DTxres = 1000;//300;
	 int DTyres = 400;
	 double DTxmin = 0.0;
	 double DTymin = -10.0;//-4.0;
	 double DTxmax = 4.5;//7.0;
	 double DTymax = 10.0;//4.0;
	//Missing Mass
	 int MMxres = 400;
	 double MMxmin = -0.2;
	 double MMxmax = 3.0;
	//Alpha
	 int alphaxres = 100;
	 double alphaxmin = 0.0;
	 double alphaxmax = 3.2;
	//Binning for Analysis
	 double Wmin = 1;
	 double Wres = 0.5;
	 double Wmax = 3;
	 double Q2min = 1.5;
	 double Q2max = 5.0;
	 double Q2res = 0.5;
	 //CC Min
	 double MinCCmin = -0.5;
	 double MinCCmax = 501.5;
	 int MinCCres = 502;
	//Electron Angle Corrections
	double min_theta_pecorr = 13.0;
	double max_theta_pecorr = 26.0;
	double res_theta_pecorr = 1.0;
	double min_phi_pecorr = -25.0;
	double max_phi_pecorr = 25.0;
	double res_phi_pecorr = 1.0;
	float min_delta_theta = -0.4;
	float max_delta_theta = 0.4;
	int bins_delta_theta = 400;
	float min_delta_p = 0.9;
	float max_delta_p = 1.1;
	int bins_delta_p = 200;


	//binning
	 float Wbin_res = 0.025;//The width of a W bin //30 steps
	 float Wbin_start = 1.4;//The starting of W bins

	 float Q2bin_res = 0.5;//6 steps
	 float Q2bin_start = 2.0; 

	double Yth_start = 9.0;
	double Yth_res = 18.0;
	double Yal_start = 18.0;//10.0; 0-36, 36-72, 72-108,  
	double Yal_res = 36.0;//40.0;
	double YM_start[3] = {1.1,1.1,0.3};
	double YM_res[3] = {0.06,0.06,0.06};

	long _ntop[4] = {0,0,0,0};//Number of selected events for each topology
	long _nptop[4] = {0,0,0,0};//Number of "missed" potential events for each topology

	long _n_wq2[3][29][5][4] = {{{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}},{{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}},{{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}},{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}}};

	long _nbad_angles[3][4][3] = {{{0,0,0},{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,0},{0,0,0}}};//topology,angle

	long _event_npass[3][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0}};

	long _ctop[4] = {0,0,0,0};//clean topologies
	long _itop[4] = {0,0,0,0};//isolated topologies
	long _cnitop[4] = {0,0,0,0};//clean but not isolated topologies
	long _caitop[4] = {0,0,0,0};//clean and isolated topologies

	long _nevnts = 0;
	long _n_attempted_evnts = 0;
	long _nvar[3] = {0,0,0};
	long _bvar[3] = {0,0,0};


	 //Making the Histograms

	TH2F_ptr_5d _WQ2_hist;
	TH2F_ptr_1d _WQ2_yield_hist;
	Bool_5d _WQ2_made;
	//TH2F_ptr_5d _Made_WQ2_hist;//[11][6][2][2];//electron cuts, topologies (including pre), Recon vs. thrown, weight (for data this should always be "Recon")
	TH2F_ptr_7d _Fid_hist;
	Bool_6d _Fid_made;
	TH2F_ptr_6d _Geo_Fid_hist;
	TH2F_ptr_6d _Geo_Fid2_hist;

	//CC Histograms
	Bool_6d _CC_made;
	TH1F_ptr_6d _CC_hist;

	//Kinematic Efficencies
	TH1F_ptr_4d _Kinematic_Eff2_hist;
	TH2F_ptr_5d _Kinematic_Eff1_hist;
	TH2F_ptr_3d _Kinematic_Eff_SC_Seg_hist;

	//SC Histograms
	TH1F_ptr_4d _SC_Eff_hist;

	//MM Histograms
	TH1F_ptr_5d _MM_hist;
	TH1F_ptr_5d _MM_lin_hist;
	Bool_4d _MM_made;

	TH2F_ptr_5d _SF_hist;
	Bool_5d _SF_made;

	//Vertex Histograms
	TH1F_ptr_4d _Vertex_hist;

	//EC Histograms
	TH1F_ptr _EC_hist;

	//Delta Histograms
	TH2F_ptr_6d _Delta_hist;

	//Electron Momentum Correction
	TH1F_ptr_3d _Pecorr_Angle_hist;
	TH1F_ptr_3d _Pecorr_Mag_hist;
	TH2F_ptr_1d _Pecorr_Angle_Dist_hist;
	TH2F_ptr_2d _Pecorr_Angle_Dist2_hist;

	//Elastic Peak
	TH1F* _Elastic_Peak_hist[6][3];//[sector][no corr, e_theta corr, all e corr]

	//Proton Energy Loss Correction
	TH1F_ptr_3d _Proton_ELoss_hist;

	//THnSparse Friend Histograms
	THnSparseD* _Friend[3][29][5];//Only doing mall for memory issues 6-28-23[5];//[30][5];//{Variable sets,topologies} => top->{pro,pip,pim,zero,all}
	THnSparseD* _Friend1[3][29][5];//Only doing mall for memory issues 6-28-23[5];//[30][5];//For positive Helicity {Variable sets,topologies} => top->{pro,pip,pim,zero,all}
	THnSparseD* _Friend2[3][29][5];//Only doing mall for memory issues 6-28-23[5];//[30][5];//For negative Helicity {Variable sets,topologies} => top->{pro,pip,pim,zero,all}
	THnSparseD* _Thrown[3][29][5];//[30][5];//Variable Sets
	THnSparseD* _Duo[3][29][5];//Systematics for multiple events counted
	THnSparseD* _Duo1[3][29][5];//Systematics for multiple events counted
	THnSparseD* _Duo2[3][29][5];//Systematics for multiple events counted
	
	//TH1D* _W_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	//TH1D* _Q2_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _MM1_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _MM2_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _Theta_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _Alpha_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _Phi_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];

	//TH1D* _W_Dist_thr[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	//TH1D* _Q2_Dist_thr[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _MM1_Dist_thr[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _MM2_Dist_thr[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _Theta_Dist_thr[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _Alpha_Dist_thr[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _Phi_Dist_thr[3][29][5];//Only doing mall for memory issues 6-28-23[5];

	TH1D* _Top_Check_hist;
	TH1D* _Top_Check_hist2[0];

	TH1F_ptr_3d _PCorr_Check_hist;

	TH1D_ptr_6d _Bin_Center_hist;//[W][Q2][Var Set][Xij][Bin Xij][W,Q2,Xij projection]
	TH1D_ptr_6d _Bin_Center2_hist1;//[W][Q2][Var Set][Xij][Bin Xij][Phi Bin][W,Q2,Xij,Phi projection]
	TH1D_ptr_6d _Bin_Center2_hist2;//[W][Q2][Var Set][Xij][Bin Xij][Phi Bin][W,Q2,Xij,Phi projection]
	TH1D_ptr_6d _Bin_Center2_hist3;//[W][Q2][Var Set][Xij][Bin Xij][Phi Bin][W,Q2,Xij,Phi projection]
	TH1D_ptr_6d _Bin_Center2_hist4;//[W][Q2][Var Set][Xij][Bin Xij][Phi Bin][W,Q2,Xij,Phi projection]


public:
	Histogram(std::shared_ptr<Flags> flags_);
	bool OK_Idx(std::vector<int> idx_);
	void Write(std::shared_ptr<Flags> flags_);
	//void Print(const std::string& output_dir_, std::shared_ptr<Flags> flags_);
	//W Qsquared plots
	int P_bin(float p_);
	float P_Min(int p_bin_);
	float P_Max(int p_bin_);
	float P_center(int p_bin_);
	int W_bins();
	int W_bin(float W_);
	float W_bot(int bin_);
	float W_top(int bin_);
	float W_center(int bin_);
	int Q2_bin(float Q2_);
	float Q2_bot(int bin_);
	float Q2_top(int bin_);
	int W_binning(float W_);
	int p_binning(float p_);
	char Part_cut(int species, int cut);
	int ecut_idx(char * ecut_, std::shared_ptr<Flags> flags_);
	int ecut_idx(const char * ecut_, std::shared_ptr<Flags> flags_);
	 //*--------------------------------Begin W Q2 Plot----------------------------*
	void WQ2_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> WQ2_cut_idx(const char* ecut_, const char * cut_, const char * top_, const char * weight_, const char * recon_, std::shared_ptr<Flags> flags_);
	bool Made_WQ2_idx(const char* ecut_, const char * cut_, const char* top_, const char * weight_, const char * recon_);
	void WQ2_Fill(float W_, float Q2_,const char* ecut_,const  char *cut_, const char * top_, const char * thrown_, std::shared_ptr<Flags> flags_, float weight_);
	void WQ2_Write(std::shared_ptr<Flags> flags_);
	 //*--------------------------------End W Q2 Plot----------------------------*

	 //*-------------------------------Start Fid Plot----------------------------*

	void Fid_Make(std::shared_ptr<Flags> flags_);
	int hcut_idx(char * hcut_, int hadron_, std::shared_ptr<Flags> flags_);
	int hcut_idx(const char * hcut_, int hadron_, std::shared_ptr<Flags> flags_);
	//bool Made_Fid_idx(const char* species_, const char* pcut_,const char * sector_, const char * cut_, const char * top_, const char * weight_);
	std::vector<int> Fid_cut_idx(const char* species_, const char* pcut_,const char * sector_, const char * cut_, const char * top_, const char * weight_,const char * p_dep_, float p_, std::shared_ptr<Flags> flags_);
	void Fid_Fill(const char * species_, float theta_, float phi_,const char* pcut_,const char * sector_,const char *cut_, const char* top_, float p_, std::shared_ptr<Flags> flags_, float weight_);
	void Fid_Write(std::shared_ptr<Flags> flags_);
	 //*--------------------------------End Fid Plot----------------------------*

	//*--------------------------------End SF Plot----------------------------*
	void SF_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> SF_idx(const char* ecut_, const char* cut_, const char* top_, const char* sector_, const char * W_dep_, std::shared_ptr<Flags> flags_, float W_=NAN);
	void SF_Fill(float p_, float sf_, float W_, const char* ecut_, const char* cut_, const char* top_, const char * sector_, float weight_, std::shared_ptr<Flags> flags_);
	std::vector<int> SF_dir_idx(const char* ecut_, const char* cut_, const char* top_, const char* sector_, const char* W_dep_, std::shared_ptr<Flags> flags_);
	void SF_Write(std::shared_ptr<Flags> flags_);
	//*--------------------------------End SF Plot----------------------------*
	//*-------------------------------Start CC Plot----------------------------*
	void CC_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> CC_idx(const char * ecut_, const char * cut_, const char* top_, const char * sector_, const char* side_, int seg_, std::shared_ptr<Flags> flags_);
	void CC_Fill(int nphe, int seg_, const char * ecut_, const char* cut_, const char* top_, const char * sector_, const char* side_, std::shared_ptr<Flags> flags_);
	void CC_Write(std::shared_ptr<Flags> flags_);
 	//*-------------------------------End CC Plot-----------------------------*

 	//*-------------------------------Start EC Plot----------------------------*
 	//*-------------------------------End EC Plot------------------------------*

 	//*-------------------------------Start Vertex Plot----------------------------*
 	void Vertex_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> Vertex_idx(const char* ecut_, const char* cut_, const char* top_, const char* sector_, std::shared_ptr<Flags> flags_);
	void Vertex_Fill(float vz_, float weight_, const char* ecut_, const char* cut_, const char* top_, const char* sector_, std::shared_ptr<Flags> flags_);
	void Vertex_Write(std::shared_ptr<Flags> flags_);
 	//*-------------------------------End Vertex Plot------------------------------*

 	//*-------------------------------Start Delta T Plot----------------------------*
 	void Delta_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> Delta_idx(float W_, const char* sector_,const char* species_, const char* pcut_, const char* cut_, const char* top_, const char* W_dep_, std::shared_ptr<Flags> flags_);
	void Delta_Fill(float p_, float dt_, float weight_, float W_, const char* species_, const char* pcut_, const char* cut_, const char* top_, const char* sector_,std::shared_ptr<Flags> flags_ );
	void Delta_Write(std::shared_ptr<Flags> flags_);
 	//*-------------------------------End Delta T Plot------------------------------*
	//*-------------------------------Start MM Plot----------------------------*
	void MM_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> MM_idx(const char* top_, const char* cut_, const char * clean_, const char * sector_, const char * W_dep_, float W_, std::shared_ptr<Flags> flags_);
	void MM_Fill(const char* top_, const char* cut_, const char * clean_, const char * sector_, float MM_, float W_, float weight_, std::shared_ptr<Flags> flags_);
	void MM_Write(std::shared_ptr<Flags> flags_);
	//*-------------------------------End MM Plot----------------------------*
	//*-------------------------------Start Kinematic Efficiency Plot----------------------------*
	void Kinematic_Eff_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> Kinematic_Eff_idx(int which_, const char* species_, const char * ecut_, const char* cut_, const char* sector_, const char* top_, std::shared_ptr<Flags> flags_);
	void Kinematic_Eff_Fill(float p_, float theta_, float weight_, const char* species_, const char* ecut_, const char* cut_, const char* sector_, const char* top_, std::shared_ptr<Flags> flags_);
	void Kinematic_Eff_Write(std::shared_ptr<Flags> flags_);
	//*-------------------------------End Kinematic Efficiency Plot----------------------------*
	//*-------------------------------Start SC Eff Plot------------------------------*
	void SC_Eff_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> SC_Eff_idx(const char * species_, const char * sector_, const char * pcut_, const char* cut_, std::shared_ptr<Flags> flags_);
	void SC_Eff_Fill(int sc_pd_, float weight_, const char* species_, const char* pcut_, const char* cut_, const char * sector_, std::shared_ptr<Flags> flags_);
	void SC_Eff_Write(std::shared_ptr<Flags> flags_);
	//*-------------------------------End SC Eff Plot------------------------------*
	//*-------------------------------Start Friend Plot----------------------------*
	std::vector<int> Friend_Bin_Sizes(std::shared_ptr<Flags> flags_);
	void Friend_Make(std::shared_ptr<Flags> flags_);
	int Friend_W_idx(float W_);
	int Friend_Q2_idx(float Q2_);
	//int Friend_MM_idx(float MM_, int var_);
	//int Friend_MM2_idx(float MM_, int var_);
	int Friend_theta_idx(float theta_);
	int Friend_alpha_idx(float alpha_);
	int Friend_phi_idx(float phi_);
	int Friend_MM_idx(float MM_, int var_, int W_bin_);
	int Friend_MM2_idx(float MM_, int var_, int W_bin_);
	std::vector<int>  Friend_idx( float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_);
	void Print_Friend_Bin(float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_, int var_);
	bool Good_Friend_Idx(float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_);
	void Friend_Fill(const char* top_, float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_, bool thrown_, float weight_, int helicity_, float cc_eff_, int top_passed_, std::shared_ptr<Flags> flags_);
	void Duo_Fill(const char* top_, float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_, bool thrown_, float weight_, int helicity_, float cc_eff_, int top_passed_, std::shared_ptr<Flags> flags_);
	void Friend_Write(std::shared_ptr<Flags> flags_);
	//*-------------------------------End Friend Plot----------------------------*
	//*-------------------------------Start PCorr Check Plot----------------------------*
	void PCorr_Check_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> PCorr_Check_idx(int sector_, const char* top_, const char* corr_, std::shared_ptr<Flags> flags_);
	void PCorr_Check_Fill(float MM2_, int sector_, const char* top_, const char* corr_, std::shared_ptr<Flags> flags_);
	void PCorr_Check_Write(std::shared_ptr<Flags> flags_);
	//*-------------------------------End PCorr Check Plot----------------------------*
	//*------------------------------Start Bin Centering Corrections------------------*
	double Xij_Bin_Min(int bin_, int Xij_, int W_bin_, int var_set_);
	double Xij_Bin_Max(int bin_, int Xij_, int W_bin_, int var_set_);
	void Bin_Centering_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> Bin_Centering_idx(float W_, float Q2_, int var_set_, int Xij_, double Xij_val_, int variable_, std::shared_ptr<Flags> flags_);
	void Bin_Centering_Fill(float W_, float Q2_, int var_set_, int Xij_, double Xij_val_, int variable_, double weight_, std::shared_ptr<Flags> flags_);
	void Bin_Centering_Write(std::shared_ptr<Flags> flags_);
	//*------------------------------End Bin Centering Corrections------------------*
	void Clean_Top_Increment(int i_);
	void Isolated_Top_Increment(int i_);
	long Clean_Top(int i_);
	long Isolated_Top(int i_);
	void Clean_Not_Isolated_Top_Increment(int i_);
	long Clean_Not_Isolated_Top(int i_);
	void Clean_and_Isolated_Top_Increment(int i_);
	long Clean_and_Isolated_Top(int i_);
	void Top_Increment(int i_);	
	void Top_Pot_Increment(int i_);	
	long NTop(int i_);
	long NPTop(int i_);
	long N_Yield(int var_, int Wbin_, int Q2bin_, int top_);
	long Bad_Angles(int var_, int top_, int angle_);
	long Bad_MM(int var_, int top_, int par_);
	long Ev_No_Pass(int var_, int top_);
	//*------------------------------Start Electron Angle Corrections------------------*
	void Ele_Angle_Corr_Make(std::shared_ptr<Flags> flags_);
	int Ele_Angle_Corr_Theta_idx(float theta_);
	int Ele_Angle_Corr_Phi_idx(float phi_);
	std::vector<int> Ele_Angle_Corr_idx(int sector_, float theta_, float phi_, std::shared_ptr<Flags> flags_);
	void Ele_Angle_Corr_Fill(int sector_, float theta_, float phi_, float W_, float thetap_, std::shared_ptr<Flags> flags_);
	void Ele_Angle_Corr_Write(std::shared_ptr<Flags> flags_);
	//*------------------------------End Electron Angle Corrections------------------*
	//*------------------------------Start Electron Momentum Magnitude Corrections------------------*
	//Please note that the theta being put in should be corrected
	void Ele_Mag_Corr_Make(std::shared_ptr<Flags> flags_);
	int Ele_Mag_Corr_Theta_idx(float theta_);
	int Ele_Mag_Corr_Phi_idx(float phi_);
	std::vector<int> Ele_Mag_Corr_idx(int sector_, float theta_, float phi_, std::shared_ptr<Flags> flags_);
	void Ele_Mag_Corr_Fill(float pe_, int sector_, float theta_, float phi_, float W_, std::shared_ptr<Flags> flags_);
	void Ele_Pro_Angle_Dist_Fill(const char* sector_, float theta_,  float thetap_, const char* top_, std::shared_ptr<Flags> flags_);
	void Ele_Mag_Corr_Write(std::shared_ptr<Flags> flags_);
	//*------------------------------End Electron Momentum Magnitude Corrections------------------*
	//*------------------------------Start Proton Energy Loss Corrections------------------*

	//*------------------------------End Proton Energy Loss Corrections------------------*
	//*------------------------------Start Elastic Peak------------------*
	void Elastic_Peak_Make(std::shared_ptr<Flags> flags_);
	void Elastic_Peak_Fill(float W_, int corr_, int sector_, std::shared_ptr<Flags> flags_);
	void Elastic_Peak_Write(std::shared_ptr<Flags> flags_);
	//*------------------------------End Elastic Peak------------------*
	//*------------------------------Start Detector Geometric Cut Histograms------------------*
	int Geo_Fid_Paddle_Idx(int det_, int seg_idx_);
	int Geo_Fid_Side_Idx(int det_, int seg_idx_);
	int Geo_Fid_Seg_Idx(int det_, int pad_idx_, int side_idx_);
	void Geo_Fid_Make(std::shared_ptr<Flags> flags_);
	//int Geo_Fid_Theta_idx(float theta_);
	//int Geo_Fid_Phi_idx(float phi_);
	std::vector<int> Geo_Fid_idx(const char* species_,const char* detector_, const char* sec_, const char* cut_, const char* pcut_, int det_seg_, int det_side_, std::shared_ptr<Flags> flags_);
	//void Geo_Fid_Fill(int sector_, float theta_, float phi_, float W_, float thetap_, std::shared_ptr<Flags> flags_);
	void Geo_Fid_Fill(float x_, float y_, float weight_, const char* species_, const char* detector_, const char* sec_ ,const char* cut_, const char* pcut_, int det_seg_, int det_side_, std::shared_ptr<Flags> flags_);
	void Geo_Fid_Write(std::shared_ptr<Flags> flags_);
	//*------------------------------End Detector Geometric Cut Histograms------------------*
	//*------------------------------Start Bin Centering Corrections for Polarization------------------*
	//double Xij_Bin_Min(int bin_, int Xij_, int W_bin_, int var_set_);
	//double Xij_Bin_Max(int bin_, int Xij_, int W_bin_, int var_set_);
	void Bin_Centering2_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> Bin_Centering2_idx(float W_, float Q2_, int var_set_, int Xij_, double Xij_val_, double phi_val_, int variable_, std::shared_ptr<Flags> flags_);
	void Bin_Centering2_Fill(float W_, float Q2_, int var_set_, int Xij_, double Xij_val_, double phi_val_, int variable_, double weight_, std::shared_ptr<Flags> flags_);
	void Bin_Centering2_Write(std::shared_ptr<Flags> flags_);
	//*------------------------------End Bin Centering Corrections for Polarization------------------*
	void Make_Top_Check(std::shared_ptr<Flags> flags_);
	void Fill_Top_Check(bool mpro_pass_, bool mpip_pass_, bool mpim_pass_, bool mzero_pass_, int npro_, int npip_, int npim_, int nzero_, std::shared_ptr<Flags> flags_);
	void Write_Top_Check(std::shared_ptr<Flags> flags_);

	//*----------- Event Idx ------------*
	//int Friend_Event_Idx(Event event_, int i_, int var_);
	int Friend_Event_Idx1(float W_);
	int Friend_Event_Idx2(float Q2_);
	int Friend_Event_Idx3(float MM_, float W_, int var_);
	int Friend_Event_Idx4(float MM2_, float W_, int var_);
	int Friend_Event_Idx5(float theta_);
	int Friend_Event_Idx6(float alpha_);
	int Friend_Event_Idx7(float phi_);

	long Number_of_Events();
	long Number_of_Attempted_Events();
	long Number_of_Filled_Var(int i_);
	long Number_of_Unfilled_Var(int i_);

	//*--------------------SC Segment Separated Kinematic Efficiency Plots ---------------------*
	void Kin_Eff_SC_Seg_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> Kin_Eff_SC_Seg_idx(const char* species_, const char* sector_, int seg_, std::shared_ptr<Flags> flags_);
	void Kin_Eff_SC_Seg_Fill(float p_, float theta_, float weight_, const char* species_, const char* sector_, int seg_, std::shared_ptr<Flags> flags_);
	void Kin_Eff_SC_Seg_Write(std::shared_ptr<Flags> flags_);
};	


 


#endif