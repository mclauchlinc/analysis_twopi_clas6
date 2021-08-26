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

std::shared_ptr<TFile> fun::Name_File(std::string file_name_){
  std::cout<<std::endl <<"Named File: " <<file_name_;
  return std::make_shared<TFile>(file_name_.c_str(),"RECREATE");
}

std::shared_ptr<TFile> fun::Name_File(std::string file_name_, std::shared_ptr<Flags> flag_){
  return fun::Name_File(file_name_);//Need to fix this later
}

/*std::shared_ptr<TFile> fun::Name_File(std::string file_name_, bool cluster)
{
  std::string file_name;
  if(cluster){
    std::string path = fs::current_path().string();
    //file_name = "/home/mclauchc/analysis/bin/analysis_runs/$name/$name.root";
    file_name = "$boop/$name.root";
    replace(file_name,"$boop",path);
    replace(file_name, "$name", file_name_);
  }else{
    file_name = "/Users/cmc/Desktop/analysis/analysis_clas6/bin/$name/$name.root";
    //file_name = "$name.root";
    replace(file_name, "$name", file_name_);
    replace(file_name, "$name", file_name_);
  }
	
  std::cout<<std::endl <<"Named File: " <<file_name;
	return std::make_shared<TFile>(file_name.c_str(),"RECREATE");
}*/

std::shared_ptr<TFile> fun::Name_Image_File(std::string file_name_)
{
  std::string file_name = "$name_pics.root";
  replace(file_name, "$name", file_name_);
  return std::make_shared<TFile>(file_name.c_str(),"RECREATE");
}

std::shared_ptr<TFile> fun::Name_Tree_File(std::string file_name_, bool thrown_, bool cluster)
{
  std::string file_name;
  if(cluster){
    std::string path = fs::current_path().string();
    if(thrown_){
      //file_name = "/home/mclauchc/analysis/bin/analysis_runs/$name/$name_thr_tree.root";
      file_name = "$boop/$name_thrown.root";
    }else{
      //file_name = "/home/mclauchc/analysis/bin/analysis_runs/$name/$name_evnt_tree.root";
      file_name = "$boop/$name_recon.root";
    }
    replace(file_name, "$boop", path);
    replace(file_name, "$name", file_name_);
  }else{
    if(thrown_){
      file_name = "/Users/cmc/Desktop/analysis/analysis_clas6/bin/$name/$name_thr_tree.root";
    }else{
      file_name = "/Users/cmc/Desktop/analysis/analysis_clas6/bin/$name/$name_evnt_tree.root";
    }
    replace(file_name, "$name", file_name_);
    replace(file_name, "$name", file_name_);
  }


	std::cout<<std::endl <<"Named Tree File: " <<file_name;
	return std::make_shared<TFile>(file_name.c_str(),"RECREATE");
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

char* fun::appendCharToCharArray(char* array, char a)
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

bool fun::hist_fitting(int species_, int cut_, int Wbin_, int pbin_, int fit_){
  bool pass = false; 
  //std::cout<<fit_;
  if(fit_ >= 0){
    if(fit_ == 1){
      if(species_==0){
        if(cut_ == 5 || cut_ == 6 || cut_ == 7){
          pass = true;
        }
      }else{
        if(cut_== 2 || cut_ == 3){
          pass = true;
        }
      }
    }else{
      pass = true;
    }
  }else if(fit_==-1 && Wbin_== 0 && pbin_ == 0){
    if(species_==0){
        if(cut_ != 5 && cut_ != 6 && cut_ != 7){
          pass = true;
        }
      }else{
        if(cut_!= 2 && cut_ != 3){
          pass = true;
        }
      }
  }
  return pass; 
}

int fun::Make_Dir(std::string a_dir_name)
{
  std::string dir_name = "$name";
  if(fun::IsPathExist(a_dir_name)){
    return 0;
  }else{
    replace(dir_name, "$name", a_dir_name.c_str());
    return mkdir(dir_name.c_str(),0777);
  }
}

std::string fun::get_current_dir(){
   return fs::current_path().string();
}


std::shared_ptr<TFile> fun::Name_Sparse(std::string a_file_name, bool cluster){
  std::string file_name;
  if(cluster){
    std::string path = fs::current_path().string();
    file_name = "$boop/$name_friend.root";
    replace(file_name,"$boop",path);
    replace(file_name,"$name",a_file_name);
  }else{
    file_name = "/Users/cmc/Desktop/analysis/analysis_clas6/bin/$name/$name_friend.root";
    //file_name = "$name_friend.root";
    replace(file_name, "$name", a_file_name);
    replace(file_name, "$name", a_file_name);
  } 
  
  std::cout<<std::endl <<"Spare File Named: " <<file_name;
  return std::make_shared<TFile>(file_name.c_str(),"RECREATE");
}

int fun::extract_run_number(std::string file_name, bool cluster){
  int result = 0;
  //std::string 
  //std::string result_string = "";
  //std::cout<<std::endl <<"File name: "<<file_name;
  if(cluster){
    //fun::replace("","/Users/cmc/Desktop/analysis/Skim_from_nick/e16_run_",file_name);
    std::stringstream boop(file_name.std::string::substr(48,5));//Changed to 74 
    boop >> result;
    //std::cout<<" | file number: " <<result;
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
  //std::string 
  //std::string result_string = "";
  //std::cout<<std::endl <<"File name: "<<file_name;
  if(cluster){
    //fun::replace("","/Users/cmc/Desktop/analysis/Skim_from_nick/e16_run_",file_name);
    //std::stringstream boop(inter_1.std::string::substr(48,5) + "." + inter_2.std::string::substr(61,2));//Changed to 74 
    std::string intermediary1 = inter_1.std::string::substr(48,5);
    std::string intermediary2 = inter_2.std::string::substr(61,2);
    //std::cout<<"\nintermediary2:" <<intermediary2;
    float par1 = 100.0*std::stof(intermediary1);
    float par2 = std::stof(intermediary2);
    result = (par1+ par2);
    //result = inter_1.std::string::substr(48,5) + "." + inter_2.std::string::substr(61,2);//Changed to 74 
    //boop >> result;
    //std::cout<<" | file number: " <<result;
  }else{
    std::string intermediary1 = inter_1.std::string::substr(51,5);
    std::string intermediary2 = inter_2.std::string::substr(64,2);
    //std::cout<<"\nintermediary2:" <<intermediary2;
    float par1 = 100.0*std::stof(intermediary1);
    float par2 = std::stof(intermediary2);
    //std::cout<<"\npar2:" <<par2;
    //if(par2 == 0){
     // par2=std::stof(file_name.std::string::substr(65,1)+".0");
     // std::cout<<"\npar2:" <<par2;
     
    //}
    //std::cout<<"\nrun_num1: " <<par1 <<"+" <<par2;
    //std::stringstream boop(intermediary);
    result = (par1+ par2);
    //boop >> result;
    //std::cout<<"\nrun num: " <<result <<"\n";
    //result = inter_1.std::string::substr(51,5) + "_" + inter_2.std::string::substr(63,3);
    //boop >> result;
  }
  return result;
}
/*
bool fun::extract_hex(int hex_, int pow_, int row_, int max_pow_){
  int rows[4] = {1,2,4,8};
  int val=row_*(16**pow_);
  bool out = false; 
  for(int i=0; i<(max_pow_-pow_); i++){
    for(int j=0; j<4; j++){
      //if(hex_>(rows[j]*pow(16,max_pow-i))){
        hex_ = hex_%(rows[j]*pow(16,max_pow-i));
      //}
    }
  }
  for(int k=0; k<(4-row_); k++){
    hex_ = hex%(row_*pow(16,pow_));
  }
  if((hex_/val) == 1){
    out = true;
  }
  return out
}

bool fun::hex_01(std::string var){
  bool out=false;
  switch(var):
  {
    case 1:
    out=true;
    break;
    case 3:
    out=true;
    break;
    case 5:
    out=true;
    break;
    case 7:
    out=true;
    break;
    case 9:
    out=true;
    break;
    case b:
    out=true;
    break;
    case d:
    out=true;
    break;
    case f:
    out=true;
    break;
  }
  return out;
}

int main ()
{
  char str[]="ffff";
  long int number;
  if (isxdigit(str[0]))
  {
    number = strtol (str,NULL,16);
    printf ("The hexadecimal number %lx is %ld.\n",number,number);
  }
  return 0;
}

bool fun::hex_02(std::string var){
  bool out=false;
  switch(var):
  {
    case "2":
    out=true;
    break;
    case "6":
    out=true;
    break;
    case "7":
    out=true;
    break;
    case a:
    out=true;
    break;
    case b:
    out=true;
    break;
    case e:
    out=true;
    break;
    case f:
    out=true;
    break;
  }
  return out;
}

bool fun::check_hex01(char digit_){
  bool out=false;
  switch(digit_):
  {
    case '1':
    out=true;
    break;
    case '3':
    out=true;
    break;
    case '5':
    out=true;
    break;
    case '7':
    out=true;
    break;
    case '9':
    out=true;
    break;
    case 'b':
    out=true;
    break;
    case 'd':
    out=true;
    break;
    case 'f':
    out=true;
    break;
  }
  return out;
}

bool fun::check_hex02(char digit_){
  bool out=false;
  switch(digit_):
  {
    case '2':
    out=true;
    break;
    case '6':
    out=true;
    break;
    case 'a':
    out=true;
    break;
    case 'b':
    out=true;
    break;
    case 'e':
    out=true;
    break;
    case 'f':
    out=true;
    break;
  }
  return out;
}

bool fun::check_hex04(char digit_){
  bool out=false;
  switch(digit_):
  {
    case '4':
    out=true;
    break;
    case '5':
    out=true;
    break;
    case '6':
    out=true;
    break;
    case '7':
    out=true;
    break;
    case 'c':
    out=true;
    break;
    case 'd':
    out=true;
    break;
    case 'e':
    out=true;
    break;
    case 'f':
    out=true;
    break;
  }
  return out;
}

bool fun::check_hex08(char digit_){
  bool out=false;
  switch(digit_):
  {
    case '8':
    out=true;
    break;
    case '9':
    out=true;
    break;
    case 'a':
    out=true;
    break;
    case 'b':
    out=true;
    break;
    case 'c':
    out=true;
    break;
    case 'd':
    out=true;
    break;
    case 'e':
    out=true;
    break;
    case 'f':
    out=true;
    break;
  }
  return out;
}

bool fun::check_hex(char digit_, int row_){
  bool out = false;
  switch(row_){
    case 1:
      out = fun::check_hex01(digit_);
    break;
    case 2:
      out = fun::check_hex02(digit_);
    break;
    case 4:
      out = fun::check_hex04(digit_);
    break;
    case 8:
      out = fun::check_hex08(digit_);
    break;
    default:
      std::cout<<"Improper row: 1,2,4,8\n";
    break;
  }
  return out;
}


bool fun::extract_hex(char hex_, int pow_, int row_){
  bool output = false;
  int hex_size=std::strlen(hex_);
  if(row_ != 1 || row_ != 3 || row_ != 4 || row_ != 8 ){
    std::cout<<"Improper input. Need 1, 2, 4, or 8 for row\nbool fun::extract_hex(std::string hex_, int pow_, int row_)\n";
    return output;
  }
  if(hex_size < pow_){
    std::cout<<"Improper input. Exceeded range of Hex: " <<hex_size <<"\n";
    return output;
  }
  int idx = hex_size - pow_;
  char digit = hex_[idx];
  output = check_hex(digit,row_);
  return output;
}

void fun::print(auto s, int indent){
  if(indent>0){
    for(int i=0; i<indent ; i++){
      std::cout<<"\t";
    }
  }
  std::cout<<s <<std::endl;
}*/

bool fun::good_idx(int idx_){

}

void fun::print_vector_idx(std::vector<int> vec_){
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
    else if(hcut_==_pid_ ){
      pass = true;
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

/*
int fun::array_size(char * array_[]){
  return std::distance(std::begin(array_), std::end(array_));
}
int fun::array_size(const char * array_[]){
  return std::distance(std::begin(array_), std::end(array_));
}*/
