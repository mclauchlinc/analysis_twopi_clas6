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

std::vector<std::string> fun::read_file_list(std::string path, int thread_num, int max_, std::shared_ptr<Flags> flags_){
  std::ifstream infile(path.c_str()); // in file stream
  std::vector<std::string> result;
  std::vector<std::string> filelist;
  std::string line;
  /*while(getline(infile,line)) { //getline sees if there is a line available
    filelist.push_back(line);//Gets the current line
  }
  int test = filelist.size();
  std::cout<<"****filelist size is " <<test <<"\n";
  if(max_ > test)
  {
    std::cout<< "You tried to add too many files. This has been corrected" <<std::endl <<"Remember that you may only add " <<test <<" files" <<std::endl;
    max_ = filelist.size();
  }*/
  std::cout<<"max files : " <<max_ <<"\n";
  int t = 0;
  bool limit = true;
  while(getline(infile,line) && limit) { //getline sees if there is a line available
    if(thread_num == (t%(flags_->Flags::Num_Cores()))){//_NUM_THREADS_)){
      result.push_back(line);//Gets the current line
    }
    t++;
    if(t >= max_ && max_>0){
      limit = false;
    }
  }
  //if(max_ == -1 || max_ > test) {//In case one tries to add too many files
    
  //}
  /*if(max_ == -1){
    std::cout<<"max is -1 so we'll do it the old way\n";
    std::cout<<"Num cores:" <<flags_->Flags::Num_Cores() <<"\n";
    int t = 0;
    bool limit = true;
    while(getline(infile,line) && limit) { //getline sees if there is a line available
      if(thread_num == (t%(flags_->Flags::Num_Cores()))){//_NUM_THREADS_)){
        result.push_back(line);//Gets the current line
      }
      t++;
      if(t >= max_){
        limit = false;
      }
    }
  }else{
    std::cout<<"Max is not -1 so we'll do size of " <<max_ <<" cool \n";
    for(int i=0; i<max_; i++){
      if(thread_num == (i%flags_->Flags::Num_Cores())){//_NUM_THREADS_)){
        result.push_back(line);//Gets the current line
      }
    }
  }*/
  
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

void fun::loadChain(std::shared_ptr<TChain> chain_, std::string file_, int thread_id_, int max_, std::shared_ptr<Flags> flags_){
  std::cout<<"Loading up the vector\n";
  std::vector<std::string> filelist = fun::read_file_list(file_,thread_id_,max_,flags_);//read_file_list(file); //creates a vector of file names
  std::cout<<"Vector Loaded\n" <<"\tit has a length of " <<filelist.size() <<"\n";
  //If not specified will take in all the files in the text file
  /*int test = filelist.size();
  if(max_ > test)
  {
    std::cout<< "You tried to add too many files. This has been corrected" <<std::endl <<"Remember that you may only add " <<test <<" files" <<std::endl;
  }
  if(max_ == -1 || max_ > test) {//In case one tries to add too many files
    max_ = filelist.size();
  }*/
  //If specified then it will take in that number of files 
  for(int i = 0; i < filelist.size(); i++) {
    //if(run_type ==3 || run_type == 4){//With some of the larger sim files this seems to be an issue where there are multiple trees in the sim files..?
    //  fun::removeTree(filelist[i]);
    //}
    std::cout<<"\r\tAdding file number " <<i <<std::flush;
    //std::cout<<"\tAdding file number " <<i <<"\n";
    chain_->AddFile(filelist[i].c_str());
  }
  std::cout<<"\nChain loaded\n";
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

int fun::extract_run_number(std::string file_name_, std::shared_ptr<Flags> flags_){
  int result = 0;
  //std::cout<<"File name: " <<file_name_ <<std::endl;
  //std::cout<<"Base: " <<flags_->Flags::Base() <<"\n";
  std::stringstream name(file_name_.std::string::substr(flags_->Flags::Base().size(),flags_->Flags::Base().size()+5));
  name >> result;
  return result;
}

int fun::extract_run_number_sim(std::string file_name_, std::shared_ptr<Flags> flags_){
  int result = 0;
  //std::cout<<"File name: " <<file_name_ <<std::endl;
  //std::cout<<"Base: " <<flags_->Flags::Base() <<"\n";
  //std::cout<<"\tbase is: " <<flags_->Flags::Base() <<"\n";
  //std::cout<<"\tsize is: " <<flags_->Flags::Base().size() <<"\n";
  std::string tmp_string = "";
  for(int i=0; i<20; i++){
    std::stringstream nametmp(file_name_.std::string::substr(flags_->Flags::Base().size()+i,1));
    //std::cout<<"Interval: " <<flags_->Flags::Base().size()+i <<" to " <<flags_->Flags::Base().size()+(1+i) <<"\n";
    //std::cout<<"interval check " <<i <<" is " <<nametmp.str() <<"\n";
    if(nametmp.str() == "_"){
      std::stringstream name(file_name_.std::string::substr(flags_->Flags::Base().size(),i));
      name >> result;
      return result;
    } 
  }
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

int fun::species_offset(const char* species_, const char* pcut_, std::shared_ptr<Flags> flags_){
  //For delta, beta, and fid
  int offset=0;
  if(pcut_==_fid_cut_){
    for(int i=0; i<fun::species_idx(species_); i++){
      if(!flags_->Flags::Plot_Fid(i)){
        offset-=1;
      }
    }
  }
  if(pcut_==_delta_cut_){
    for(int i=0; i<fun::species_idx(species_); i++){
      if(!flags_->Flags::Plot_Delta(i)){
        offset-=1;
      }
    }
  }
  if(pcut_==_geo_sc_cut_){
    for(int i=0; i<fun::species_idx(species_); i++){
      if(!flags_->Flags::Plot_Fid_Geo(i,1)){
        offset-=1;
      }
    }
  }
  if(pcut_==_geo_ec_cut_){
    for(int i=0; i<fun::species_idx(species_); i++){
      if(!flags_->Flags::Plot_Fid_Geo(i,2)){
        offset-=1;
      }
    }
  }
  if(pcut_==_kin_eff_cut_){
    for(int i=0; i<fun::species_idx(species_); i++){
      if(!flags_->Flags::Plot_Kin_Eff(i)){
        offset-=1;
      }
    }
  }
  /*if(pcut_==_beta_cut_){
    for(int i=0; i<fun::species_idx(species_);i++){
      if(!flags_->Flags::Plot_Beta(i)){
        offset-=1;
      }
    }
  }*/
  return offset;
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
    }else if(ecut_ == _sc_eff_cut_ && flags_->Flags::SC_Eff()){
      pass = true;
    }else if(ecut_==_vertex_cut_ && flags_->Flags::Vertex_Cut()){
      pass = true;
    }else if(ecut_ == _id_cut_ && flags_->Flags::ID_Cut()){
      pass = true;
    }else if(ecut_ == _geo_cc_cut_ && flags_->Flags::Geo_Cut(0,0)){
      pass = true;
    }else if(ecut_ == _geo_sc_cut_ && flags_->Flags::Geo_Cut(0,1)){
      pass = true;
    }else if(ecut_ == _geo_ec_cut_ && flags_->Flags::Geo_Cut(0,2)){
      pass = true;
    }else if(ecut_ == _kin_eff_cut_ && flags_->Flags::Kin_Eff_Cut(0)){
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
    if(hcut_ == _none_ || hcut_==_sanity_ || hcut_==_pid_ ){
      pass = true;
    }else{
      if(hcut_==_fid_cut_ && flags_->Flags::Fid_Cut(fun::species_idx(species_))){
        pass = true;
      }else if(hcut_==_delta_cut_ && flags_->Flags::Delta_Cut(fun::species_idx(species_))){
        pass = true;
      }else if(hcut_ == _sc_eff_cut_ && flags_->Flags::SC_Eff()){
        pass = true;
      }else if(hcut_ == _id_cut_ && flags_->Flags::ID_Cut()){
        pass = true;
      }else if(hcut_ == _geo_sc_cut_ && flags_->Flags::Geo_Cut(fun::species_idx(species_),1)){
        pass = true;
      }else if(hcut_ == _geo_ec_cut_ && flags_->Flags::Geo_Cut(fun::species_idx(species_),2)){
        pass = true;
      }else if(hcut_ == _kin_eff_cut_ && flags_->Flags::Kin_Eff_Cut(fun::species_idx(species_))){
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

bool fun::is_empty(int run_num_, std::shared_ptr<Flags> flags_){
  switch(flags_->Run()){
    case 0:
      for(int i=0; i<std::distance(std::begin(_empty_e16_), std::end(_empty_e16_)); i++){
        if(run_num_ == _empty_e16_[i]){
          return true;
        }
      }
    break;
    case 1:
      for(int i=0; i<std::distance(std::begin(_empty_e1f_), std::end(_empty_e1f_)); i++){
        if(run_num_ == _empty_e16_[i]){
          return true;
        }
      }
    break;
    default:
      std::cout<<"<in is_empty> Improper Run Group Given " <<flags_->Run() <<"\n";
    break;
  }
  return false;
}

bool fun::is_full(int run_num_, std::shared_ptr<Flags> flags_){
  std::cout<<"Checking run num: " <<run_num_ <<" for run: " <<run_groups[flags_->Run()] <<"\n";
    switch(flags_->Run()){
    case 0:
      if(run_num_ >= _e16_run_bounds_[0] && run_num_ <= _e16_run_bounds_[1]){
        for(int i=0; i<std::distance(std::begin(_empty_e16_), std::end(_empty_e16_)); i++){
          if(run_num_ == _empty_e16_[i]){
            return false;
          }
        }
        return true;
      }else{
        std::cout<<"Wasn't within the range of e16\n";
      }
    break;
    case 1:
      std::cout<<"is full e1f bounds: " <<_e1f_run_bounds_[0] <<"-" <<_e1f_run_bounds_[2] <<"\n";
      if(run_num_ >= _e1f_run_bounds_[0] && run_num_ <= _e1f_run_bounds_[1]){
        for(int i=0; i<std::distance(std::begin(_empty_e1f_), std::end(_empty_e1f_)); i++){
          if(run_num_ == _empty_e16_[i]){
            return false;
          }
        }
        return true;
      }else{
        std::cout<<"Wasn't within the range of e1f\n";
      }
    break;
    default:
      std::cout<<"<in is_full> Improper Run Group Given " <<flags_->Run() <<"\n";
    break;
  }
  return false;
}

bool fun::correct_run_num(int run_num_, std::shared_ptr<Flags> flags_){
  std::cout<<"*Checking run number*\n";
  if(flags_->Flags::Fill()){//Was the target filled?
    std::cout<<"\tTarget Filled\n";
    return fun::is_full(run_num_,flags_);
  }else{
    std::cout<<"\tTarget Empty\n";
    return fun::is_empty(run_num_,flags_);
  }
}

int fun::real_helicity(int hel_, int run_num_, std::shared_ptr<Flags> flags_){
  int hel = hel_*_plate_sign_[flags_->Run()]; 
  switch(flags_->Run()){
    case 0:
      for(int i=0; i<std::distance(std::begin(_plate_swap_e16_), std::end(_plate_swap_e16_)); i++){
        if(run_num_ >= _plate_swap_e16_[i] ){
          hel = hel*-1;
        }
      }
      return hel; 
    break;
    case 1:
      for(int i=0; i<std::distance(std::begin(_plate_swap_e1f_), std::end(_plate_swap_e1f_)); i++){
        if(run_num_ >= _plate_swap_e1f_[i]){
          hel = hel*-1;
        }
      }
      return hel; 
    break;
    default:
      std::cout<<"Entered wrong run <fun::real_helicity>\n";
      return 0; 
    break;
  }
}

bool fun::is_num_in_list(int num_, const int list_[]){
  bool present = false;
  std::cout<<"size of array: " <<(sizeof(list_)/sizeof(list_[0])) <<"\n";
  
  return present;
}

bool fun::correct_run(int run_num_, std::shared_ptr<Flags> flags_){
  if(flags_->Sim()){ return true;}
  bool correct = false;
  if(flags_->Run()==_e16_){//e16
    if(flags_->Fill()){//
      /*for(int i=0; i<(sizeof(_beam_e16_)/sizeof(_beam_e16_[0])); i++){
        if(run_num_ == _beam_e16_[i]){
          return true;
        }
      }*/
      for(int i=0; i<(sizeof(_good_beam_e16_)/sizeof(_good_beam_e16_[0])); i++){
        if(run_num_ == _good_beam_e16_[i]){
          return true;
        }
      }
    }else{
      /*for(int i=0; i<(sizeof(_empty_e16_)/sizeof(_empty_e16_[0])); i++){
        if(run_num_ == _empty_e16_[i]){
          return true;
        }
      }*/
      for(int i=0; i<(sizeof(_good_empty_e16_)/sizeof(_good_empty_e16_[0])); i++){
        if(run_num_ == _good_empty_e16_[i]){
          return true;
        }
      }
    }
  }else if(flags_->Run()==_e1f_){//e1f
    //std::cout<<"We are doing e1f!\n";
    if(flags_->Fill()){//
      for(int i=0; i<(sizeof(_beam_e1f_)/sizeof(_beam_e1f_[0])); i++){
        if(run_num_ == _beam_e1f_[i]){
          return true;
        }
      }
    }else{
      for(int i=0; i<(sizeof(_empty_e1f_)/sizeof(_empty_e1f_[0])); i++){
        if(run_num_ == _empty_e1f_[i]){
          return true;
        }
      }
    }
  }
  return correct;
}

  float fun::poly_4(float x_, float a_, float b_, float c_, float d_, float e_){
    return a_*pow(x_,4) + b_*pow(x_,3) + c_*pow(x_,2) + d_*x_ + e_;
  }

  float fun::poly_3(float x_, float a_, float b_, float c_, float d_){
    return a_*pow(x_,3) + b_*pow(x_,2) + c_*x_ + d_;
  }

  float fun::poly_2(float x_, float a_, float b_, float c_){
    return a_*pow(x_,2) + b_*x_,2 + c_;
  }

  bool fun::check_sec(const char * sec_, float x_, float y_){
    return sec_ == _sector_[physics::get_sector(physics::get_phi(x_,y_))-1];
  }

/*
int fun::array_size(char * array_[]){
  return std::distance(std::begin(array_), std::end(array_));
}
int fun::array_size(const char * array_[]){
  return std::distance(std::begin(array_), std::end(array_));
}*/

int fun::ele_cut_width(const char* cut_, std::shared_ptr<Flags> flags_){
  if(cut_ == _fid_cut_){
    return flags_->Flags::Fid_Cut_Width(0);
  }else if(cut_ == _delta_cut_){
    return flags_->Flags::Delta_Cut_Width(0);
  }else if(cut_ == _geo_cc_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(0,0);
  }else if(cut_ == _geo_sc_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(0,1);
  }else if(cut_ == _geo_ec_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(0,2);
  }else if(cut_ == _vertex_cut_){
    return flags_->Flags::Vertex_Cut_Width();
  }else if(cut_ == _kin_eff_cut_){
    return flags_->Flags::Kin_Eff_Cut_Width(0);
  }else if(cut_ == _cc_cut_){
    return flags_->Flags::Min_CC_Cut_Width();
  }else if(cut_ == _sf_cut_){
    return flags_->Flags::SF_Cut_Width();
  }else{
    return 1;
  }
}
int fun::pro_cut_width(const char* cut_, std::shared_ptr<Flags> flags_){
  if(cut_ == _fid_cut_){
    return flags_->Flags::Fid_Cut_Width(1);
  }else if(cut_ == _delta_cut_){
    return flags_->Flags::Delta_Cut_Width(1);
  }else if(cut_ == _geo_cc_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(1,0);
  }else if(cut_ == _geo_sc_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(1,1);
  }else if(cut_ == _geo_ec_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(1,2);
  }else if(cut_ == _kin_eff_cut_){
    return flags_->Flags::Kin_Eff_Cut_Width(1);
  }else{
    return 1;
  }
}
int fun::pip_cut_width(const char* cut_, std::shared_ptr<Flags> flags_){
  if(cut_ == _fid_cut_){
    return flags_->Flags::Fid_Cut_Width(2);
  }else if(cut_ == _delta_cut_){
    return flags_->Flags::Delta_Cut_Width(2);
  }else if(cut_ == _geo_cc_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(2,0);
  }else if(cut_ == _geo_sc_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(2,1);
  }else if(cut_ == _geo_ec_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(2,2);
  }else if(cut_ == _kin_eff_cut_){
    return flags_->Flags::Kin_Eff_Cut_Width(2);
  }else{
    return 1;
  }
}
int fun::pim_cut_width(const char* cut_, std::shared_ptr<Flags> flags_){
  if(cut_ == _fid_cut_){
    return flags_->Flags::Fid_Cut_Width(3);
  }else if(cut_ == _delta_cut_){
    return flags_->Flags::Delta_Cut_Width(3);
  }else if(cut_ == _geo_cc_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(3,0);
  }else if(cut_ == _geo_sc_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(3,1);
  }else if(cut_ == _geo_ec_cut_){
    return flags_->Flags::Fid_Geo_Cut_Width(3,2);
  }else if(cut_ == _kin_eff_cut_){
    return flags_->Flags::Kin_Eff_Cut_Width(3);
  }else{
    return 1;
  }
}
int fun::cut_width(const char* species_, const char* cut_, std::shared_ptr<Flags> flags_){
  if(species_ == _ele_){
    return fun::ele_cut_width(cut_,flags_);
  }else if(species_ == _pro_){
    return fun::pro_cut_width(cut_,flags_);
  }else if(species_ == _pip_){
    return fun::pip_cut_width(cut_,flags_);
  }else if(species_ == _pim_){
    return fun::pim_cut_width(cut_,flags_);
  }else{
    return 1;
  }
}

int fun::geo_det_idx(const char* detector_){
  if(detector_ == _cc_){
    return 0;
  }else if(detector_ == _sc_){
    return 1;
  }else if(detector_ == _ec_){
    return 2;
  }else{
    return -1;
  }
}

bool fun::vector_in_vector_of_vectors(std::vector<std::vector<int>> vov_, std::vector<int> vec_){
  if(vov_.size()==0){return false;}
  for(int i=0; i<vov_.size(); i++){
    if(vov_[i] == vec_){
      return true;
    }
  }
  return false;
}

bool fun::idx_in_vector_of_idx(std::vector<int> vec_, int idx_){
  if(vec_.size()==0){return false;}
  for(int i=0; i<vec_.size(); i++){
    if(vec_[i] == idx_){
      return true;
    }
  }
  return false;
}

bool fun::Half_Wave(int run_num_, std::shared_ptr<Flags> flags_){
	bool plate_in = false;
	int above = 0;
	int below = 0;
	bool at = false; 
	if(flags_->Run()==_e16_){
		for(int i=0; i<sizeof(_plate_swap_e16_); i++){
			if(run_num_ > _plate_swap_e16_[i]){
				below +=1; 
			}else if(run_num_ == _plate_swap_e16_[i]){
				at = true; 
			}else if(run_num_ < _plate_swap_e16_[i]){
				above += 1; 
			}
		}
	}else if(flags_->Run()==_e1f_){
		for(int i=0; i<sizeof(_plate_swap_e1f_); i++){
			if(run_num_ > _plate_swap_e1f_[i]){
				below +=1; 
			}else if(run_num_ == _plate_swap_e1f_[i]){
				at = true; 
			}else if(run_num_ < _plate_swap_e1f_[i]){
				above += 1; 
			}
		}
	}
	if(at){
		if(below%2 == 0){
			plate_in = true;
		}else{
			plate_in = false;
		}
	}else{
		if(below%2 == 0){
			plate_in = false;
		}else{
			plate_in = true;
		}
	}
	return plate_in;
}


//Correct Helicity accoridng to half wave plate status
float fun::Corr_Helicity(float helicity_, int run_num_, std::shared_ptr<Flags> flags_){
	
	int eh = 0; 
	if(helicity_ >= 1000) eh = 1; 
	if(helicity_ <= -1000) eh = -1; 
	if(helicity_ < 1000 && helicity_ > -1000) eh = 0; 
	//if(plate_stat == 0 ) eh = 1; 
	if(fun::Half_Wave(run_num_,flags_)){
		return eh*-1;
	}else{
		return  eh;
	}
}