#ifndef FUNCTIONS_H_GUARD
#define FUNCTIONS_H_GUARD

#include "constants.hpp"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include "branches.hpp"
#include "histogram.hpp"
#include <sys/stat.h>
#include <unistd.h>



namespace fun {
	bool IsPathExist(const std::string &s);
	bool replace(std::string& str, const std::string& from, const std::string& to);
	std::shared_ptr<TFile> Name_File(std::string file_name_);
	std::vector<std::string> read_file_list(std::string path_);
	void loadChain(std::shared_ptr<TChain> c, std::string file, int max);
	char* appendCharToCharArray(char* array, char a);
	int extract_run_number(std::string file_name, std::string data_set_);
	float extract_run_number_float(std::string file_name, std::string data_set_);
	int run_number(std::string file_name_, std::string front_);
	int run_segment(std::string file_name_, std::string front_, std::string mid_);
	float run_num_seg(int run_num_, int run_seg_);
	int run_num_idx(int run_num_,std::vector<int> run_nums_);
}

#endif