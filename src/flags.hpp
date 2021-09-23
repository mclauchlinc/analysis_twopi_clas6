#ifndef FLAGS_HPP
#define FLAGS_HPP

#include <stdio.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <stdlib.h>

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


static std::string _all_f_ = "all";
//Target Flags
static std::string _target_ = "-t=";
static std::string _run_group_[2] = {"e16","e1f"}; //Run group
static std::string _filled_ = "filled"; //Was the target filled
static std::string _sim_ = "sim"; //Is this simulated data
static std::string _location_ = "-loc="; //path to text file containing all file locations
static std::string _helicity_ = "hel"; //Are we separating out beam helicity?
//PID Cuts Performed
static std::string _cut_f_ = "-cut=";
static std::string _ncut_ = "-ncut=";
static std::string _fid_cut_f_[4] = {"fid_ele","fid_pro","fid_pip","fid_pim"}; //{ele,pro,pip,pim}
static std::string _delta_cut_f_[4] = {"dt_ele","dt_pro","dt_pip","dt_pim"};; //{ele,pro,pip,pim}
static std::string _sf_cut_f_ = "sf";
static std::string _cc_cut_f_ = "cc";
static std::string _ec_cut_f_ = "ec";
static std::string _beta_cut_f_[4] = {"beta_ele","beta_pro","beta_pip","beta_pim"};; //{ele,pro,pip,pim}
//Other Cuts
static std::string _vertex_cut_f_ = "vertex";
//Event Selection
static std::string _mm_cut_f_[4] = {"mm_pro,mm_pip,mm_pim,mm_zero"};
//Efficiency Cuts
static std::string _eff_cut_ = "-eff=";
static std::string _neff_cut_ = "-neff=";
static std::string _geo_cut_ = "-geo=";
static std::string _ngeo_cut_ = "-ngeo=";
static std::string _eff_cc_ = "cc";
static std::string _eff_dc_ = "dc";
static std::string _eff_ec_ = "ec";
static std::string _eff_sc_ = "sc";
//Corrections
static std::string _corr_ = "-corr=";
static std::string _ncorr_ = "-ncorr=";
static std::string _p_corr_ = "p_corr";//Momentum
static std::string _e_corr_ = "e_corr";//Energy
static std::string _v_corr_ = "v_corr";//Vertex
//Histograms
static std::string _plot_ = "-plot=";
static std::string _nplot_ = "-nplot=";
static std::string _plot_wq2_ = "wq2";
static std::string _plot_fid_[4] = {"fid_ele","fid_pro","fid_pip","fid_pim"}; //{ele,pro,pip,pim}
static std::string _plot_sf_ = "sf";
static std::string _plot_cc_ = "cc";
static std::string _plot_dt_[4] = {"dt_ele","dt_pro","dt_pip","dt_pim"}; //{ele,pro,pip,pim}
static std::string _plot_cc_eff_ = "cc_eff";
static std::string _plot_sc_eff_ = "sc_eff";
static std::string _plot_ec_eff_ = "ec_eff";
static std::string _plot_dc_eff_ = "dc_eff";
static std::string _plot_cc_geo_ = "cc_geo";
static std::string _plot_sc_geo_ = "sc_geo";
static std::string _plot_ec_geo_ = "ec_geo";
static std::string _plot_mm_[4] = {"mm_pro","mm_pip","mm_pim","mm_zero"}; //{pro,pip,pim,zero}
static std::string _plot_beta_[4] = {"beta_ele","beta_pro","beta_pip","beta_pim"}; //{ele,pro,pip,pim}
static std::string _plot_vertex_ = "vertex";
//Histogram Separation
//Will have space for this
//THnSparse
static std::string _make_friend_ = "-friend";
//Other Flags
static std::string _output_name_ = "-name=";
static std::string _num_files_ = "-n=";



class Flags {
private:
	//Data Information
	int _run_group = -1; //{0,1}->{e1-6,e1f}
	bool _filled = false; //Was the target filled
	bool _sim = false; //Is this simulated data
	std::string _location = ""; //
	bool _helicity = false; //Are we separating out beam helicity?
	int _num_files = -1;
	//PID Cuts Performed
	bool _cut_all = false;
	bool _fid_cut[4] = {false,false,false,false}; //{ele,pro,pip,pim}
	bool _delta_cut[4] = {false,false,false,false}; //{ele,pro,pip,pim}
	bool _sf_cut = false;
	bool _cc_cut = false;
	bool _ec_cut = false;
	bool _beta_cut[4] = {false,false,false,false}; //{ele,pro,pip,pim}
	//Other Cuts
	bool _vertex_cut = false;
	//Event Selection
	bool _mm_cut[4] = {false,false,false,false}; //{pro,pip,pim,zero}
	//Efficiency Cuts
	bool _cc_eff_cut = false;
	bool _ec_eff_cut = false;
	bool _sc_eff_cut = false;
	bool _dc_eff_cut = false;
	bool _eff_all = false;
	bool _cc_geo_cut = false;
	bool _ec_geo_cut = false;
	bool _sc_geo_cut = false;
	bool _geo_all = true;
	//Corrections
	bool _corr_all = false;
	bool _p_corr = false;//Momentum
	bool _e_corr = false;//Energy
	bool _v_corr = false;//Vertex
	//Histograms
	bool _plot_all = false;
	bool _plot_wq2 = false;
	bool _plot_fid[4] = {false,false,false,false}; //{ele,pro,pip,pim}
	bool _plot_sf = false;
	bool _plot_cc = false;
	bool _plot_dt[4] = {false,false,false,false}; //{ele,pro,pip,pim}
	bool _plot_cc_geo = false;
	bool _plot_sc_geo = false;
	bool _plot_ec_geo = false;
	bool _plot_mm[4] = {false,false,false,false}; //{pro,pip,pim,zero}
	bool _plot_beta[4] = {false,false,false,false}; //{ele,pro,pip,pim}
	bool _plot_cc_eff = false;
	bool _plot_sc_eff = false;
	bool _plot_ec_eff = false;
	bool _plot_dc_eff = false;
	bool _plot_vertex = false;

	//Histogram Separation
		//Will have space for this
	//THnSparse
	bool _make_friend = false;
	std::string _output_name = "";
public:
	Flags();
	void Help();
	void Read_Flags(int argc, char** argv);
	void ID_Flag(std::string str);
	void Cut_Flag(std::string str, bool include);
	void Run_Flag(std::string str);
	void Plot_Flag(std::string str, bool include);
	void Eff_Flag(std::string str, bool include);
	void Geo_Flag(std::string str, bool include);
	void Corr_Flag(std::string str, bool include);
	//Run Info
	int Run();
	bool Sim();
	bool Fill();
	bool Helicity();
	std::string Files();
	int Num_Files();
	//PID Cuts 
	bool Fid_Cut(int particle);
	bool Delta_Cut(int particle);
	bool SF_Cut();
	bool CC_Cut();
	bool EC_Cut();
	bool Beta_Cut(int particle);
	//Other Cuts
	bool Vertex_Cut();
	//Event Selection 
	bool MM_Cut(int topology);
	//Efficiency Cuts
	bool CC_Eff();
	bool EC_Eff();
	bool SC_Eff();
	bool DC_Eff();
	bool CC_Geo();
	bool EC_Geo();
	bool SC_Geo();
	//Corrections
	bool P_Corr();
	bool Delta_Corr();
	bool E_Corr();
	bool Vertex_Corr();
	//Histograms
	bool Plot_WQ2();
	bool Plot_Fid(int particle);
	bool Plot_SF();
	bool Plot_CC();
	bool Plot_Delta(int particle);
	bool Plot_CC_Geo();
	bool Plot_EC_Geo();
	bool Plot_SC_Geo();
	bool Plot_MM(int topology);
	bool Plot_Beta(int particle);
	bool Plot_Vertex();
	//Portions of Histograms
	bool Ele_Cut(const char * ecut_);

	//Histogram Separation

	//THnSparse
	bool Make_Friend();

	//Output File
	std::string Output_Name();

};


#endif