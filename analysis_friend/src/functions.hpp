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

namespace fs = std::experimental::filesystem;


namespace fun {
double* Vector_Array(std::vector<double> vec);
int* Vector_Array(std::vector<int> vec);
bool replace(std::string& str, const std::string& from, const std::string& to);
std::shared_ptr<TFile> Name_File(std::string a_file_name);
double nSparseIntegral(THnSparseD* nhist_);

};

#endif