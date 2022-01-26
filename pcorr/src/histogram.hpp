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
#include "THnSparse.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include "TGraph.h"
#include <iterator>
#include "cuts.hpp"
#include <cmath>
#include "TLatex.h"
#include <mutex>//For pausing the multithreading for the filling of THnSparse histograms 
//#include <unistd.h>//to allow us to use chdir
//#include "TImage.h"
//#include "particle.hpp"
//#include "variables.h"
//#include "CartesianGenerator.hh"

//Binning Info
static float _theta_min_ = 10.0;
static float _theta_max_ = 30.0;
static int _theta_bins_ = 40;

static float _phi_min_ = -25.0;
static float _phi_max_ = 25.0;
static int _phi_bins_ = 100;

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

using TDir_ptr_1d = std::vector<TDirectory*>;
using TDir_ptr_2d = std::vector<std::vector<TDirectory*>>;
using TDir_ptr_3d = std::vector<std::vector<std::vector<TDirectory*>>>;
using TDir_ptr_4d = std::vector<std::vector<std::vector<std::vector<TDirectory*>>>>;
using TDir_ptr_5d = std::vector<std::vector<std::vector<std::vector<std::vector<TDirectory*>>>>>;

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
	std::shared_ptr<TFile> _HistImageFile;

	TCanvas* def;
	float _delta_theta_max = 0.2;
	float _delta_theta_min = -0.2;
	int _delta_theta_res = 200;

	float _elast_max = 1.4;
	float _elast_min = 0.5;
	int _elast_res = 300;

	TH1F_ptr_4d _Delta_Theta_hist;
	TH1F_ptr_3d _Elast_hist;

public:
	Histogram(std::shared_ptr<Flags> flags_);
	bool OK_Idx(std::vector<int> idx_);
	void Write(std::shared_ptr<Flags> flags_);
	//void Print(const std::string& output_dir_, std::shared_ptr<Flags> flags_);
	//W Qsquared plots
	int Phi_Idx(float phi_);
	float Phi_Low(int bin_);
	float Phi_Top(int bin_);
	float Phi_Center(int bin_);
	int Theta_Idx(float theta_);
	float Theta_Low(int bin_);
	float Theta_Top(int bin_);
	float Theta_Center(int bin_);
	//*------------------------------- Start Plot 1 Electron Angle Correction ---------------------------------*
	void ECorr_Angle_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> ECorr_Angle_idx(float theta_, float phi_, int sector_, const char * corr_ , std::shared_ptr<Flags> flags_);
	void ECorr_Angle_Fill(float delta_theta_, float theta_, float phi_, int sector_, const char* corr_, std::shared_ptr<Flags> flags_);
	void ECorr_Angle_Write(std::shared_ptr<Flags> flags_);
	//*------------------------------- End Plot 1 Electron Angle Correction ---------------------------------*
	void Elastic_Make(std::shared_ptr<Flags> flags_);
	std::vector<int> Elastic_idx(int sector_, const char* corr_, const char* pro_thresh_, std::shared_ptr<Flags> flags_);
	void Elastic_Fill(float W_, int sector_, const char* corr_, const char* pro_thresh_, std::shared_ptr<Flags> flags_);
	void Elast_Write(std::shared_ptr<Flags> flags_);
};





#endif