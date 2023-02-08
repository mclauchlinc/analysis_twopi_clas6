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
#include "histogram.hpp"
#include <sys/stat.h>
#include <unistd.h>
#include <experimental/filesystem>
#include "THnSparse.h"
#include "CartesianGenerator.hpp"

namespace fs = std::experimental::filesystem;


namespace fun {
double* Vector_Array(std::vector<double> vec);
int* Vector_Array(std::vector<int> vec);
bool replace(std::string& str, const std::string& from, const std::string& to);
std::shared_ptr<TFile> Name_File(std::string a_file_name);
double nSparseIntegral(THnSparseD* nhist_);
THnSparseD* Add_THnSparse(THnSparseD* hist1_, THnSparseD* hist2_, int sign_, std::vector<int> num_bins_);
THnSparseD* Localized_Holes(THnSparseD* exp_hist_, THnSparseD* sim_hist_, THnSparseD* sim_hole_hist_, std::vector<int> num_bins_);
std::vector<std::vector<int>> Surrounding_Bin(int* bin_, int dist_, std::vector<int> num_bins_);

};

#endif