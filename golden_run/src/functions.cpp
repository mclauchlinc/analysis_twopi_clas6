#include "functions.hpp"

bool fun::IsPathExist(const std::string &s){
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}

bool fun::replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::shared_ptr<TFile> fun::Name_File(std::string file_name_)
{
  std::cout<<"Creating Output File: " <<file_name_ <<"\n";
	return std::make_shared<TFile>(file_name_.c_str(),"RECREATE");
}


std::vector<std::string> fun::read_file_list(std::string path_){
  std::ifstream infile(path_.c_str()); // in file stream
  std::vector<std::string> result;
  std::string line;
  int t = 0;
  while(getline(infile,line)) { //getline sees if there is a line available
    result.push_back(line);//Gets the current line
    t++;
  }
  return result;
}


void fun::loadChain(std::shared_ptr<TChain> c, std::string file, int max)
{
  std::vector<std::string> filelist = fun::read_file_list(file);//read_file_list(file); //creates a vector of file names
  //If not specified will take in all the files in the text file
  int test = filelist.size();
  if(max > test)
  {
    std::cout<< "You tried to add too many files. This has been corrected" <<std::endl <<"Remember that you may only add " <<test <<" files" <<std::endl;
  }
  if(max == -1 || max > test) {//In case one tries to add too many files
    max = filelist.size();
  }
  //If specified then it will take in that number of files 
  for(int i = 0; i < max; i++) {
    c->AddFile(filelist[i].c_str());
  }
}

char* fun::appendCharToCharArray(char* array, char a)
{
    size_t len = strlen(array);

    char* ret = new char[len+2];

    strcpy(ret, array);    
    ret[len] = a;
    ret[len+1] = '\0';

    return ret;
}

int fun::extract_run_number(std::string file_name, std::string data_set_){//Right now only optimized for experimental e16, might need to change numbers
  int result = 0;
  std::stringstream boop(file_name.std::string::substr(string_cut1_map[data_set_],5));
  boop >> result;
  return result;
}

float fun::extract_run_number_float(std::string file_name, std::string data_set_){
  float result;
  std::string inter_1 = file_name;
  std::string inter_2 = file_name;
  std::string intermediary1 = inter_1.std::string::substr(string_cut1_map[data_set_],5);
  std::string intermediary2 = inter_2.std::string::substr(string_cut2_map[data_set_],2);
  float par1 = 100.0*std::stof(intermediary1);
  float par2 = std::stof(intermediary2);
  //std::cout<<"\nRun: "<<par1 <<" | file: " <<par2 <<" | cut1:" <<string_cut1_map[data_set_] <<" | cut2: " <<; 
  result = (par1+ par2);
  return result;
}

int fun::run_number(std::string file_name_, std::string front_){
  int result = 0;
  std::stringstream name(file_name_.std::string::substr(front_.size(),front_.size()+5));
  name >> result;
  return result;
}

int fun::run_segment(std::string file_name_, std::string front_, std::string mid_){
  int result = 0;
  std::stringstream name(file_name_.std::string::substr(front_.size()+5+mid_.size(),front_.size()+5+mid_.size()+2));
  name >> result;
  return result;
}

float fun::run_num_seg(int run_num_, int run_seg_){
  return (float)run_num_ + (float)run_seg_/100.0;
}

int fun::run_num_idx(int run_num_,std::vector<int> run_nums_){
  for(int i=0; i<run_nums_.size(); i++){
    if(run_nums_[i] == run_num_){
      return i;
    }
  }
  return -1;
}

