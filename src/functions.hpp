#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "constants.hpp"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "branches.hpp"
//#include "environment.hpp"
#include "histogram.hpp"
#include <sys/stat.h>
#include <unistd.h>
#include <experimental/filesystem>
#include "flags.hpp"
#include <iterator>

namespace fs = std::experimental::filesystem;

static std::string _file_name_ = "/Users/cmc/Desktop/analysis/analysis_clas6/bin/$name/$name.root";



namespace fun {
bool IsPathExist(const std::string &s);

bool replace(std::string& str, const std::string& from, const std::string& to);

std::shared_ptr<TFile> Name_File(std::string file_name_);
std::shared_ptr<TFile> Name_File(std::string file_name_, std::shared_ptr<Flags> flag_);
//std::shared_ptr<TFile> Name_File(std::string file_name_, bool cluster=false);
std::shared_ptr<TFile> Name_Image_File(std::string file_name_);
std::shared_ptr<TFile> Name_Tree_File(std::string file_name_, bool thrown_ = false, bool cluster = false);

std::vector<std::string> read_file_list(std::string path, int thread_num);

void removeTree(std::string file_name);

void loadChain(std::shared_ptr<TChain> chain_, std::string file_, int thread_id_, int max_);

char* appendCharToCharArray(char* array, char a);

bool no_pro_pip_match(int idx1, int idx2[20]);

bool hist_fitting(int species_, int cut_, int Wbin_, int pbin_, int fit_);

int Make_Dir(std::string a_dir_name);

std::string get_current_dir();

std::shared_ptr<TFile> Name_Sparse(std::string a_file_name, bool cluster = false);

int extract_run_number(std::string file_name, bool cluster);

float extract_run_number_float(std::string file_name, bool cluster);

//bool extract_hex(int hex_, int pow_, int row_, int max_pow_);

//void print(auto s, int indent=0);

bool good_idx(int idx_);
void print_vector_idx(std::vector<int> vec_);
void print_vector_idx(std::vector<long> vec_);
int top_idx(const char* top_);
int ecut_idx(const char* ecut_);
int hcut_idx(const char* hcut_);
int weight_idx(const char* hcut_);
int recon_idx(const char* recon_);
int cut_idx(const char* cut_);
int species_idx(const char* species_);
int sector_idx(const char* sector_);
int ecut_offset(const char * ecut_, std::shared_ptr<Flags> flags_);
bool ecut_perform(const char* ecut_, std::shared_ptr<Flags> flags_);
int hcut_offset(const char * species_, const char * hcut_, std::shared_ptr<Flags> flags_);
bool hcut_perform(const char * species_,const char* ecut_, std::shared_ptr<Flags> flags_);
bool pcut_perform(const char * species_, const char* pcut_, std::shared_ptr<Flags> flags_);
int pcut_offset(const char * species_, const char * pcut_, std::shared_ptr<Flags> flags_);
int clean_idx(const char * clean_);
int sim_idx(bool sim_);
bool top_perform(const char* top_, std::shared_ptr<Flags> flags_);
int truth_idx(bool pass_);
int top_offset(const char * top_, std::shared_ptr<Flags> flags_);
//int array_size(char* array_[]);
//int array_size(const char* array_[]);
}

#endif