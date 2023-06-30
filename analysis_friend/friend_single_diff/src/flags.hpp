#ifndef FLAGS_HPP
#define FLAGS_HPP

#include <stdio.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <stdlib.h>
#include "constants.hpp"

/*
These flags serve to tell the analysis a number of conditions 
necessary to perform the analysis, such as the run group, whether 
experimental or simulated data, and which cuts to perform. 
Libraries that use this include:
- event_analysis.hpp
- pid.hpp
- particle.hpp
- main.hpp

*/
static std::string _give_top_ = "-t=";//Determined which file is extracted
static std::string _real_top_ = "-r=";//Tells program which topology is truly being looked at

//Topologies are listed in Constants.hpp

static std::string _histograms_ = "-h="; 
static std::string _all_ = "all";
static std::string _polarization_ = "pol";
static std::string _single_diff_ = "single_diff";
static std::string _beam_spin_ = "beam_spin";
static std::string _eff_ = "eff";
static std::string _err_ = "err";
static std::string _wq2_ = "wq2";
static std::string _acceptance_ = "accept";
static std::string _pol_mm1_ = "pol_mm1";
static std::string _pol_mm2_ = "pol_mm2";
static std::string _pol_theta_ = "pol_theta";
static std::string _pol_alpha_ = "pol_alpha";
static std::string _error_ = "error";
static std::string _local_holes_ = "local_holes";
static std::string _plot_local_holes_ = "plot_local_holes";


static std::string _output_name_ = "-name=";
static std::string _output_name1_ = "-name1=";
static std::string _output_name2_ = "-name2=";
static std::string _output_name3_ = "-name3=";
static std::string _sim_file_ = "-sim=";
static std::string _exp_file_ = "-exp=";
//static std::string _weight_file_ = "-weight=";
static std::string _empty_file_ = "-empty=";
static std::string _sim_no_rad_file_ = "-no_rad=";
static std::string _sim_file2_ = "-sim2=";
static std::string _exp_file2_ = "-exp2=";
//static std::string _weight_file2_ = "-weight2=";
static std::string _empty_file2_ = "-empty2=";

static std::string _sim_no_rad_file2_ = "-no_rad2=";
static std::string _loc_holes1_ = "-holes=";
static std::string _loc_holes2_ = "-holes2=";
static std::string _no_loc_holes_ = "-nolocalholes";
//static std::string _exp_pos_file_ = "expp=";
//static std::string _exp_pos_file2_ = "expp2=";
//static std::string _exp_neg_file_ = "expn=";
//static std::string _exp_neg_file2_ = "expn2=";

//Other Info about the runf
static std::string _info_ = "-i=";
static std::string _run_e16_ = "e16";
static std::string _run_e1f_ = "e1f";
static std::string _both_ = "both";
static std::string _var_pim_ = "var_pim";
static std::string _var_pro_ = "var_pro";
static std::string _var_pip_ = "var_pip";
static std::string _flux_included_ = "flux_weighted";
static std::string _helicity_ = "helicity";

//Luminosity
static std::string _lumen_ = "-l=";
static std::string _lumen2_ = "-l2=";
static std::string _charge1_= "-q1=";
static std::string _charge2_= "-q2=";

//Charge Ratio
static std::string _qratio_ = "-qr=";
static std::string _qratio2_ = "-qr2=";

static std::string _image_name_ = "-image=";








class Flags {
private:
	//Run Information
	int _run_group = -1; //{0,1}->{e1-6,e1f}
	bool _has_empty = false;
	bool _has_exp = false;
	bool _has_sim = false;
	bool _has_weight = false;

	//File Locations/Names
	std::string _sim_loc = "";//Simulation file for
	std::string _sim_no_rad_loc = "";//Simulation file with no radiative corrections
	std::string _exp_loc = "";//Experimental data file, ideally weighted with CC efficiencies and virtual photon flux
	std::string _empty_loc = "";//Experimental empty target file
	std::string _weight_loc = "";//I don't remember the use for this, but I believe it comes into play for error analysis
	std::string _sim_loc2 = "";
	std::string _sim_no_rad_loc2 = "";
	std::string _exp_loc2 = "";
	std::string _exp_pos_loc = "";//Experimental Data for postitive helicity
	std::string _exp_neg_loc = "";//Experimental Data for negative helicity
	std::string _exp_pos_loc2 = "";
	std::string _exp_neg_loc2 = "";
	std::string _empty_loc2 = "";
	std::string _weight_loc2 = "";
	std::string _output_name1 = "";
	std::string _output_name2 = "";
	std::string _output_name3 = "";
	std::string _output_name="";
	std::string _image_name = "";
	std::string _localized_holes_name = "";
	std::string _localized_holes_name2 = "";

	//Procedure Flags
	bool _helicity = false;
	bool _rad_corr = false;
	bool _flux_inc = false;

	//Plots
	bool _plot_all = false;
	bool _plot_polarization = false;
	bool _plot_single_diff = false;
	bool _plot_beam_spin = false;
	bool _plot_eff = false;
	bool _plot_err = false;
	bool _plot_wq2 = false;
	bool _plot_acceptance = false;
	bool _plot_pol_mm1 = false;
	bool _plot_pol_mm2 = false;
	bool _plot_pol_theta = false;
	bool _plot_pol_alpha = false;
	bool _plot_error = false;
	bool _localized_holes = false;
	bool _plot_localized_holes = false;
	bool _has_localized_holes = false;

	std::string _real_top = "";
	std::string _input_top = "";

	bool _make_image = false;

	std::string _var_set = "";
	int _var_idx = -1; 

	float _luminosity = NAN;//Integrated luminosity for the given run
	float _luminosity2 = NAN;
	double _charge1=NAN;
	double _charge2=NAN;


	float _Qr = NAN; //Integrated Charge Ratio between filled and empty target 
	float _Qr2 = NAN;

	
public:
	Flags();
	void Help();
	void Read_Flags(int argc_, char** argv_);
	void ID_Flag(std::string str_);
	void Plot_Flag(std::string str_);
	void Info_Flag(std::string str_);
	int Run();
	std::string Run(std::string str_);
	//Plotting
	bool Plot_Pol();
	bool Plot_Single_Diff();
	bool Plot_Beam_Spin();
	bool Plot_Eff();
	bool Plot_Err();
	bool Plot_WQ2();
	bool Plot_Acceptance();
	//File Names
	std::string Sim_File();
	std::string Exp_File();
	std::string Empty_File();
	std::string Weight_File();
	std::string Sim_File2();
	std::string Exp_File2();
	std::string Empty_File2();
	std::string Weight_File2();
	std::string Output_File();
	std::string Image_File();
	std::string Exp_Hel();
	std::string Exp_Hel2();
	std::string Sim_No_Rad();
	std::string Sim_No_Rad2();
	std::string Holes_File();
	std::string Holes_File2();
	//Have Things
	bool Has_Weight();
	bool Has_Sim();
	bool Has_Exp();
	bool Has_Empty();
	bool Make_Image();
	int Var_idx();
	std::string Var_Set();
	std::string Top();
	int Top_idx();
	bool Helicity();
	bool Rad_Corr();
	bool Flux_Included();
	float Qr();//Charge Ratio for exp target vs. no target
	float Qr2();
	double L(int i_);
	bool Plot_Polarization(int i_);
	bool Error();
	bool Localized_Holes();
	bool Plot_Localized_Holes();
	bool Has_Localized_Holes();
};


#endif