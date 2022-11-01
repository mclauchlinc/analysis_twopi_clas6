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

static std::string _output_name_ = "-name=";
static std::string _sim_file_ = "-sim=";
static std::string _exp_file_ = "-exp=";
static std::string _weight_file_ = "-weight=";
static std::string _empty_file_ = "-empty=";
static std::string _sim_file2_ = "-sim2=";
static std::string _exp_file2_ = "-exp2=";
static std::string _weight_file2_ = "-weight2=";
static std::string _empty_file2_ = "-empty2=";

//Other Info about the run
static std::string _info_ = "-i=";
static std::string _run_e16_ = "e16";
static std::string _run_e1f_ = "e1f";
static std::string _both_ = "both";
static std::string _var_pim_ = "var_pim";
static std::string _var_pro_ = "var_pro";
static std::string _var_pip_ = "var_pip";

//Luminosity
static std::string _lumen_ = "-l=";
static std::string _lumen2_ = "-l2=";

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
	std::string _sim_loc = "";
	std::string _exp_loc = "";
	std::string _empty_loc = "";
	std::string _weight_loc = "";
	std::string _sim_loc2 = "";
	std::string _exp_loc2 = "";
	std::string _empty_loc2 = "";
	std::string _weight_loc2 = "";
	std::string _output_name = "";
	std::string _image_name = "";

	//Plots
	bool _plot_all = false;
	bool _plot_polarization = false;
	bool _plot_single_diff = false;
	bool _plot_beam_spin = false;
	bool _plot_eff = false;
	bool _plot_err = false;
	bool _plot_wq2 = false;
	bool _plot_acceptance = false;

	std::string _real_top = "";
	std::string _input_top = "";

	bool _make_image = false;

	std::string _var_set = "";
	int _var_idx = -1; 

	float _luminosity = NAN;
	float _luminosity2 = NAN;


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

};


#endif