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

std::shared_ptr<TFile> fun::Name_File(std::shared_ptr<Flags> flags_){
  return std::make_shared<TFile>(flags_->Flags::Output_Name().c_str(),"RECREATE");
}

std::shared_ptr<TFile> fun::Name_Image(std::shared_ptr<Flags> flags_){
  if(flags_->Flags::Make_Image()){
    std::cout<<std::endl <<"Named Image File: " <<flags_->Flags::Image_Name().c_str();
    return std::make_shared<TFile>(flags_->Flags::Image_Name().c_str(),"RECREATE");
  }
}

std::shared_ptr<TFile> fun::Name_Sparse(std::shared_ptr<Flags> flags_){
  if(flags_->Flags::Make_Friend()){
    std::cout<<"Named Sparse File: " <<flags_->Flags::Friend_Name() <<"\n";
    return std::make_shared<TFile>(flags_->Flags::Friend_Name().c_str(),"RECREATE");
  }
}

std::vector<std::string> fun::read_file_list(std::string path, int thread_num){
  std::ifstream infile(path.c_str()); // in file stream
  std::vector<std::string> result;
  std::string line;
  int t = 0;
  while(getline(infile,line)) { //getline sees if there is a line available
    if(thread_num == (t%_NUM_THREADS_)){
      result.push_back(line);//Gets the current line
    }
    t++;
  }
  return result;
}
void fun::removeTree(std::string file_name){
  TFile *file=new TFile((file_name).c_str(),"update");
  std::string object_to_remove="h10;1";
  //the object can be a tree, a histogram, etc, in this case "test1" is a TTree
  //notice the ";1" which means cycle 1; to remove all cycles do ";*"
  //if your object is not at the top directory, but in a directory in the .root file, called foo
  // you do first
  //file->cd("foo");
  //then continue with the Delete command which is only applied to the current gDirectory
  gDirectory->Delete(object_to_remove.c_str());
  file->Close();
}

void fun::loadChain(std::shared_ptr<TChain> chain_, std::string file_, int thread_id_, int max_){
  std::vector<std::string> filelist = fun::read_file_list(file_,thread_id_);//read_file_list(file); //creates a vector of file names
  //If not specified will take in all the files in the text file
  int test = filelist.size();
  if(max_ > test)
  {
    std::cout<< "You tried to add too many files. This has been corrected" <<std::endl <<"Remember that you may only add " <<test <<" files" <<std::endl;
  }
  if(max_ == -1 || max_ > test) {//In case one tries to add too many files
    max_ = filelist.size();
  }
  //If specified then it will take in that number of files 
  for(int i = 0; i < max_; i++) {
    //if(run_type ==3 || run_type == 4){//With some of the larger sim files this seems to be an issue where there are multiple trees in the sim files..?
    //  fun::removeTree(filelist[i]);
    //}
    chain_->AddFile(filelist[i].c_str());
  }
}

/*char* fun::appendCharToCharArray(char* array, char a)
{
    size_t len = strlen(array);

    char* ret = new char[len+2];

    strcpy(ret, array);    
    ret[len] = a;
    ret[len+1] = '\0';

    return ret;
}

bool fun::no_pro_pip_match(int idx1, int idx2[20]){//Designed to check to see if there was a double identification on a pip/proton
  bool pass = true;
  for(int i=0; i< 20; i++){
    if(idx1 == idx2[i]){
      pass = false;
    }
  }
  return pass; 
}


int fun::Make_Dir(std::string a_dir_name){
  std::string dir_name = "$name";
  if(fun::IsPathExist(a_dir_name)){
    return 0;
  }else{
    replace(dir_name, "$name", a_dir_name.c_str());
    return mkdir(dir_name.c_str(),0777);
  }
}
*/

int fun::extract_run_number(std::string file_name, bool cluster){
  int result = 0;
  if(cluster){
    std::stringstream boop(file_name.std::string::substr(48,5));//Changed to 74 
    boop >> result;
  }else{
    std::stringstream boop(file_name.std::string::substr(51,5));
    boop >> result;
  }
  return result;
}

float fun::extract_run_number_float(std::string file_name, bool cluster){
  float result;
  std::string inter_1 = file_name;
  std::string inter_2 = file_name;
  if(cluster){
    std::string intermediary1 = inter_1.std::string::substr(48,5);
    std::string intermediary2 = inter_2.std::string::substr(61,2);
    float par1 = 100.0*std::stof(intermediary1);
    float par2 = std::stof(intermediary2);
    result = (par1+ par2);
  }else{
    std::string intermediary1 = inter_1.std::string::substr(51,5);
    std::string intermediary2 = inter_2.std::string::substr(64,2);
    float par1 = 100.0*std::stof(intermediary1);
    float par2 = std::stof(intermediary2);
    result = (par1+ par2);
  }
  return result;
}

void fun::print_vector_idx(std::vector<int> vec_){
  std::cout<<"\tIndex: ";
  for(int i=0; i<vec_.size(); i++){
    std::cout<<vec_[i] <<" ";
  }
  std::cout<<"\n";
}

void fun::print_vector_idx(std::vector<long> vec_){
  std::cout<<"\tIndex: ";
  for(int i=0; i<vec_.size(); i++){
    std::cout<<vec_[i] <<" ";
  }
  std::cout<<"\n";
}

int fun::top_idx(const char* top_){
  int top_idx = -1;
  for(int i = 0; i<std::distance(std::begin(_top_), std::end(_top_)); i++){
    if(_top_[i] == top_){
      top_idx = i;
    }
  }
  return top_idx;
}

int fun::ecut_idx(const char* ecut_){
  int ecut_idx = -1;
  for(int i = 0; i<std::distance(std::begin(_ecuts_),std::end(_ecuts_)); i++){
    if(_ecuts_[i] == ecut_){
      ecut_idx = i;
    }
  }
  return ecut_idx;
}
int fun::hcut_idx(const char* hcut_){
  int hcut_idx = -1;
  for(int i = 0; i<std::distance(std::begin(_hcuts_), std::end(_hcuts_)); i++){
    if(_hcuts_[i] == hcut_){
      hcut_idx = i;
    }
  }
  return hcut_idx;
}

int fun::weight_idx(const char* weight_){
  int weight_idx = -1;
  for(int i = 0; i<std::distance(std::begin(_weight_), std::end(_weight_)); i++){
    if(_weight_[i] == weight_){
      weight_idx = i;
    }
  }
  return weight_idx;
}
int fun::recon_idx(const char* recon_){
  int recon_idx = -1;
  for(int i = 0; i<std::distance(std::begin(_recon_), std::end(_recon_)); i++){
    if(_recon_[i] == recon_){
      recon_idx = i;
    }
  }
  return recon_idx;
}
int fun::cut_idx(const char* cut_){
  int cut_idx = -1;
  for(int i = 0; i<std::distance(std::begin(_cut_), std::end(_cut_)); i++){
    if(_cut_[i] == cut_){
      cut_idx = i;
    }
  }
  return cut_idx;
}

int fun::species_idx(const char* species_){
  int species_idx = -1;
  for(int i = 0; i<std::distance(std::begin(_species_), std::end(_species_)); i++){
    if(_species_[i] == species_){
      species_idx = i;
    }
  }
  return species_idx;
}

int fun::sector_idx(const char* sector_){
  int sector_idx = -1;
  for(int i = 0; i<std::distance(std::begin(_sector_), std::end(_sector_)); i++){
    if(_sector_[i] == sector_){
      sector_idx = i;
    }
  }
  return sector_idx;
}

int fun::ecut_offset(const char * ecut_, std::shared_ptr<Flags> flags_){
  int offset = 0; 
  if(fun::ecut_perform(ecut_,flags_)){
    int ecut_idx= fun::ecut_idx(ecut_);
    for(int i=0; i<ecut_idx; i++){
      if(!fun::ecut_perform(_ecuts_[i],flags_)){
        offset -= 1; 
      }
    }
  }
  return offset; 
}

bool fun::ecut_perform(const char* ecut_, std::shared_ptr<Flags> flags_){
  bool pass = false;
  if(ecut_ == _none_ || ecut_==_sanity_ || ecut_==_pid_){
    pass = true;
  }else{
    if(ecut_==_fid_cut_ && flags_->Flags::Fid_Cut(0)){
      pass = true;
    }else if(ecut_==_sf_cut_ && flags_->Flags::SF_Cut()){
      pass = true;
    }else if(ecut_==_cc_cut_ && flags_->Flags::CC_Cut()){
      pass = true;
    }else if(ecut_==_ec_cut_ && flags_->Flags::EC_Cut()){
      pass = true;
    }else if(ecut_==_vertex_cut_ && flags_->Flags::Vertex_Cut()){
      pass = true;
    }else if(ecut_==_delta_cut_ && flags_->Flags::Delta_Cut(0)){
      pass = true;
    }else if(ecut_==_event_){
      pass = true;
    }
  }
  return pass; 
}
int fun::hcut_offset(const char * species_, const char * hcut_, std::shared_ptr<Flags> flags_){
  int offset = 0; 
  if(fun::hcut_perform(species_,hcut_,flags_)){
    int hcut_idx= fun::hcut_idx(hcut_);
    for(int i=0; i<hcut_idx; i++){
      if(!fun::hcut_perform(species_,_hcuts_[i],flags_)){
        offset -= 1; 
      }
    }
  }
  return offset; 
}

bool fun::hcut_perform(const char * species_,const char* hcut_, std::shared_ptr<Flags> flags_){
  bool pass = false;
  if(species_ != _ele_){
    if(hcut_ == _none_ || hcut_==_sanity_ || hcut_==_pid_){
      pass = true;
    }else{
      if(hcut_==_fid_cut_ && flags_->Flags::Fid_Cut(fun::species_idx(species_))){
        pass = true;
      }else if(hcut_==_delta_cut_ && flags_->Flags::Delta_Cut(fun::species_idx(species_))){
        pass = true;
      }else if(hcut_==_event_ ){
        pass = true;
      }
    }
  }
  return pass; 
}

bool fun::pcut_perform(const char * species_, const char* pcut_, std::shared_ptr<Flags> flags_){
  bool pass = false;
  if(species_==_ele_){
    pass = fun::ecut_perform(pcut_,flags_);
  }else{
    pass = fun::hcut_perform(species_,pcut_,flags_);
  }
  return pass; 
}

int fun::pcut_offset(const char * species_, const char * pcut_, std::shared_ptr<Flags> flags_){
  int offset = 0;
  if(species_==_ele_){
    offset = fun::ecut_offset(pcut_,flags_);
  }else{
    offset = fun::hcut_offset(species_,pcut_,flags_);
  }
  return offset; 
}

int fun::clean_idx(const char * clean_){
  int clean_idx = -1;
  for(int i = 0; i<std::distance(std::begin(_clean_event_), std::end(_clean_event_)); i++){
    if(_clean_event_[i] == clean_){
      clean_idx = i;
    }
  }
  return clean_idx;
}

int fun::sim_idx(bool sim_){
  int sim_idx = -1;
  if(sim_){
    sim_idx = 1;
  }else{
    sim_idx =0;
  }
  return sim_idx;
}

bool fun::top_perform(const char* top_, std::shared_ptr<Flags> flags_){
  bool pass = false;
  if(top_ == _mall_ || top_==_mnone_){
    pass = true;
  }else if(flags_->Flags::MM_Cut(fun::top_idx(top_))){
    pass = true;
  }
  return pass; 
}

int fun::truth_idx(bool pass_){
  if(pass_){
    return 1;
  }else{
    return 0;
  }
}

int fun::top_offset(const char * top_, std::shared_ptr<Flags> flags_){
  int top_idx = fun::top_idx(top_);
  int idx = 0; 
  if(top_ == _mnone_){
    return -top_idx;
  }else{
    for(int i=0; i<top_idx; i++){
      if(!flags_->Flags::MM_Cut(i)){
        idx -= 1; 
      }
    }
    return idx; 
  }
}

int fun::cc_side_idx(const char * side_){
  int side_idx = -1; 
  if(side_ == _left_){
    side_idx = 0;
  }else if(side_ == _coinc_){
    side_idx = 1;
  }else if(side_ == _right_){
    side_idx = 2;
  }
  return side_idx;
}

/*
int fun::array_size(char * array_[]){
  return std::distance(std::begin(array_), std::end(array_));
}
int fun::array_size(const char * array_[]){
  return std::distance(std::begin(array_), std::end(array_));
}*/
