#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD

#include <unordered_map>
#include "TMath.h"
#include "TLorentzVector.h"

//static int num_mixed_p_pip = 0; 

static const int MAX_PARTS = 20; 

static const float c_special = 29.9792458; //speed of light in cm/ns
static const float c_convert = 10000000; //Convert c_special to m/s

//Beam energies in GeV
static const float energy_e16 = 5.754;//5.7696;5.754
//3_14 was 5.754, but led to bad results: updated to 5.759 after reading Paremuzyan's thesis 
//Check 4.794 GeV for e16
static const float energy_e1f = 5.499;

//Masses of relevant particles
static const	float me = 0.0005109989; //mass of electron in GeV
static const	float mp = 0.93828;	//Mass of proton in GeV
static const	float mpi = 0.1395;	//Mass of pion in GeV

static const std::string x_e16_local = "exp_e16_local";
static const std::string x_e1f_local = "exp_e1f_local";
static const std::string empty_e16_local = "empty_e16_local";
static const std::string empty_e1f_local = "empty_e1f_local";
static const std::string x_e16_cluster = "exp_e16_cluster";
static const std::string x_e1f_cluster = "exp_e1f_cluster";
static const std::string empty_e16_cluster = "empty_e16_cluster";
static const std::string empty_e1f_cluster = "empty_e1f_cluster";

static const std::string path_x_e16_local = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/NickSkim_e16_PlateIN.txt";
static const std::string path_x_e1f_local = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/NickSkim_e16_PlateIN.txt";
static const std::string path_empty_e16_local = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/NickSkim_e16_PlateIN.txt";
static const std::string path_empty_e1f_local = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/NickSkim_e16_PlateIN.txt";
static const std::string path_x_e16_cluster = "/home/mclauchc/analysis/Path_Files/cluster_e16_path.txt";
static const std::string path_x_e1f_cluster = "/home/mclauchc/analysis/Path_Files/cluster_e1f_file_paths.txt";
static const std::string path_empty_e16_cluster = "/home/mclauchc/analysis/Path_Files/empty_cluster_e16_file_paths.txt";
static const std::string path_empty_e1f_cluster = "/home/mclauchc/analysis/Path_Files/empty_cluster_e1f_sim_files.txt";

static std::unordered_map<std::string, std::string> filepath_map = 	{	{x_e16_local,path_x_e16_local},
																		{x_e1f_local,path_x_e1f_local},
																		{empty_e16_local,path_empty_e16_local},
																		{empty_e1f_local,path_empty_e1f_local},
																		{x_e16_cluster,path_x_e16_cluster},
																		{x_e1f_cluster,path_x_e1f_cluster},
																		{empty_e16_cluster,path_empty_e16_cluster},
																		{empty_e1f_cluster,path_empty_e1f_cluster}
																	};

static std::unordered_map<std::string, bool> fileloc_map = 	{	{x_e16_local,true},
																		{x_e1f_local,true},
																		{empty_e16_local,true},
																		{empty_e1f_local,true},
																		{x_e16_cluster,false},
																		{x_e1f_cluster,false},
																		{empty_e16_cluster,false},
																		{empty_e1f_cluster,false}
																	};

static std::unordered_map<std::string, int> histbound_map = 	{	{x_e16_local,0},
																		{x_e1f_local,1},
																		{empty_e16_local,2},
																		{empty_e1f_local,3},
																		{x_e16_cluster,0},
																		{x_e1f_cluster,1},
																		{empty_e16_cluster,2},
																		{empty_e1f_cluster,3}
																	};

static std::unordered_map<std::string, int> string_cut1_map = 	{	{x_e16_local,51},
																		{x_e1f_local,51},//need to update
																		{empty_e16_local,51},//need to update
																		{empty_e1f_local,51},//need to update
																		{x_e16_cluster,48},
																		{x_e1f_cluster,34},//need to update
																		{empty_e16_cluster,48},//need to update
																		{empty_e1f_cluster,48}//need to update
																	};
static std::unordered_map<std::string, int> string_cut2_map = 	{	{x_e16_local,64},
																		{x_e1f_local,64},//need to update
																		{empty_e16_local,64},//need to update
																		{empty_e1f_local,64},//need to update
																		{x_e16_cluster,61},
																		{x_e1f_cluster,47},//need to update
																		{empty_e16_cluster,61},//need to update
																		{empty_e1f_cluster,61}//need to update
																	};



#endif
