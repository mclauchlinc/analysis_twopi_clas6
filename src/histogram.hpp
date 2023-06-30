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
static float _mm_min_[4] = {0.7,0.0,0.0,-0.2};
static float _mm_max_[4] = {1.5,0.3,0.3,0.2};
static int _mm_bin_[4] = {200,200,200,200};
//MM2
static float _mm2_min_[4] = {0.7,0.0,0.0,-0.02};
static float _mm2_max_[4] = {1.5,0.09,0.09,0.02};
static int _mm2_bin_[4] = {100,100,100,100};
//EC
static float _ec_min_ = 0.0;
static float _ec_max_ = 6.0;
static int _ec_bin_ = 200; 

//Detector Plots
//Vertex
static float _vertex_min_ = -12.0;
static float _vertex_max_ = 5.0;
static int _vertex_bin_ = 200;
//CC Eff
	//Momentum vs. Theta
static float _cc_eff1_xmin_ = 0.0;//GeV
static float _cc_eff1_xmax_ = 5.0;//GeV
static int _cc_eff1_xbin_ = 300;
static float _cc_eff1_ymin_ = 0.0;//Degrees
static float _cc_eff1_ymax_ = 80.0;//Degrees
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
static float _MM_min_[3] = {1.1000,0.2780,1.1000};//{1.1,0.3,1.1};//Changed 5/5/23
static float _MM_max_[3] = {_W_max_-_mpi_,_W_max_-_mp_,_W_max_-_mpi_};//{2.0,1.1,2.0};//Changed 5/5/23
static int _MM_bins_ = 14;
static float _MM2_min_[3] = {0.2780,1.1000,1.1000};//{0.3,1.1,1.1};//Changed 5/5/23
static float _MM2_max_[3] = {_W_max_-_mp_,_W_max_-_mpi_,_W_max_-_mpi_};//{1.1,2.0,2.0};//Changed 5/5/23
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

	 //Making the Histograms

	TH2F_ptr_5d _WQ2_hist;
	Bool_5d _WQ2_made;
	//TH2F_ptr_5d _Made_WQ2_hist;//[11][6][2][2];//electron cuts, topologies (including pre), Recon vs. thrown, weight (for data this should always be "Recon")
	TH2F_ptr_7d _Fid_hist;
	Bool_6d _Fid_made;

	//CC Histograms
	Bool_6d _CC_made;
	TH1F_ptr_6d _CC_hist;

	//Kinematic Efficencies
	TH1F_ptr_4d _Kinematic_Eff2_hist;
	TH2F_ptr_5d _Kinematic_Eff1_hist;

	//SC Histograms
	TH1F_ptr_4d _SC_Eff_hist;

	//MM Histograms
	TH1F_ptr_5d _MM_hist;
	Bool_4d _MM_made;

	TH2F_ptr_5d _SF_hist;
	Bool_5d _SF_made;

	//Vertex Histograms
	TH1F_ptr_4d _Vertex_hist;

	//EC Histograms
	TH1F_ptr _EC_hist;

	//Delta Histograms
	TH2F_ptr_6d _Delta_hist;

	//THnSparse Friend Histograms
	THnSparseD* _Friend[3][29][5];//Only doing mall for memory issues 6-28-23[5];//[30][5];//{Variable sets,topologies} => top->{pro,pip,pim,zero,all}
	THnSparseD* _Friend1[3][29][5];//Only doing mall for memory issues 6-28-23[5];//[30][5];//For positive Helicity {Variable sets,topologies} => top->{pro,pip,pim,zero,all}
	THnSparseD* _Friend2[3][29][5];//Only doing mall for memory issues 6-28-23[5];//[30][5];//For negative Helicity {Variable sets,topologies} => top->{pro,pip,pim,zero,all}
	THnSparseD* _Thrown[3][29][5];//[30][5];//Variable Sets
	
	//TH1D* _W_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	//TH1D* _Q2_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _MM1_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _MM2_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _Theta_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _Alpha_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];
	TH1D* _Phi_Dist[3][29][5];//Only doing mall for memory issues 6-28-23[5];


	TH1F_ptr_3d _PCorr_Check_hist;


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
	bool Made_Fid_idx(const char* species_, const char* pcut_,const char * sector_, const char * cut_, const char * top_, const char * weight_);
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
	int Friend_MM_idx(float MM_, int var_);
	int Friend_MM2_idx(float MM_, int var_);
	int Friend_theta_idx(float theta_);
	int Friend_alpha_idx(float alpha_);
	int Friend_phi_idx(float phi_);
	std::vector<int>  Friend_idx( float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_);
	void Print_Friend_Bin(float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_, int var_);
	void Friend_Fill(const char* top_, float W_, float Q2_, float MM_, float MM2_, float theta_, float alpha_, float phi_ , int var_, bool thrown_, float weight_, int helicity_, float cc_eff_, std::shared_ptr<Flags> flags_);
	void Friend_Write(std::shared_ptr<Flags> flags_);
	//*-------------------------------End Friend Plot----------------------------*
	//*-------------------------------Start PCorr Check Plot----------------------------*
	void PCorr_Check_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> PCorr_Check_idx(int sector_, const char* top_, const char* corr_, std::shared_ptr<Flags> flags_);
	void PCorr_Check_Fill(float MM2_, int sector_, const char* top_, const char* corr_, std::shared_ptr<Flags> flags_);
	void PCorr_Check_Write(std::shared_ptr<Flags> flags_);
	//*-------------------------------End PCorr Check Plot----------------------------*
	
};





#endif