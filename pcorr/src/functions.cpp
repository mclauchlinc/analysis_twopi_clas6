#include "functions.hpp"


std::shared_ptr<TFile> fun::Name_File(std::shared_ptr<Flags> flags_){
  return std::make_shared<TFile>(flags_->Flags::Output_Name().c_str(),"RECREATE");
}

float fun::theta_calc(float theta_p_, float beam_energy_, bool radians_){
	float theta = theta_p_*_radian_;
	float denom = (beam_energy_ + _mp_)*TMath::Tan(theta);
	if(radians_){
		return 2*TMath::ATan(_mp_/denom);
	}else{
		return 2*TMath::ATan(_mp_/denom)*_degree_;
	}
}

float fun::delta_theta(float theta_e_, float theta_p_, float beam_energy_){
	return fun::theta_calc(theta_p_,beam_energy_,false) - theta_e_;
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

float fun::theta(std::shared_ptr<Branches> data_, int idx_, bool radians_){
	if(radians_){
		return TMath::ACos(data_->cz(idx_));
	}
	return TMath::ACos(data_->cz(idx_))*_degree_;
}

int fun::get_sector(float phi_){
	int sector;
	if(phi_>-30 && phi_ <=30){
		sector = 1;
	}else if(phi_>30 && phi_<=90){
		sector = 2;
	}else if(phi_>90 && phi_ <=150){
		sector = 3;
	}else if(phi_>150 || phi_<=-150){
		sector = 4;
	}else if(phi_>-150 && phi_<=-90){
		sector = 5;
	}else if(phi_>-90 && phi_<=-30){
		sector = 6;
	}//got rid of pointless "else" statement 3/10/2017
	return sector;
}

int fun::get_sector(std::shared_ptr<Branches> data_, int idx_){
	return get_sector(fun::phi(data_,idx_,false,false));
}


float fun::phi_center( float phi_){
	double phi_corr;
	int sector = get_sector(phi_);
	switch(sector){
		case 1:
			phi_corr = phi_;
		break;
		case 2:
			phi_corr = phi_-60.0;
		break;
		case 3:
			phi_corr = phi_-120.0;
		break;
		case 4:
			if(phi_<=-150.0){
				phi_corr = phi_+180.0;
			}else if(phi_>=150.0){
				phi_corr = phi_-180.0;
			}
		break;
		case 5:
			phi_corr = phi_+120.0;
		break;
		case 6:
			phi_corr = phi_+60.0;
		break;
		default:
			std::cout<<"Improper phi. Cannot determine sector\n";
		break;
	}
	//Not working, but elegant. Adjust for elegance later
	/*phi0 = phi0 +180;
	int mod6 = ((int)phi0+210)/60;
	phi_corr = phi0 - (double)mod6*60.0;
	*/
	return phi_corr;
}


float fun::phi(std::shared_ptr<Branches> data_, int idx_, bool center_, bool radians_){
	float phi = TMath::ATan2(data_->cy(idx_),data_->cx(idx_));
	//std::cout<<"\tPhi start (deg): " <<phi*_degree_ <<"\n";
	if(center_){
		if(radians_){
			return phi_center(phi*_degree_)/_degree_;
		}else{
			//std::cout<<"\tPhi Centerd (deg): " <<phi_center(phi*_degree_) <<"\n";
			return phi_center(phi*_degree_);
		}
	}else{
		if(radians_){
			return phi;
		}else{
			return phi*_degree_;
		}
	}
}

float fun::p_calc(float theta_e_, float beam_energy_){
	float denom = 1+ (2*beam_energy_ * pow(TMath::Sin(theta_e_*_radian_/2),2))/_mp_;
	return beam_energy_/denom;
}

float fun::delta_p_e(float p_e_, float theta_e_, float beam_energy_){
	return fun::p_calc(theta_e_, beam_energy_) - p_e_;
}

TLorentzVector fun::Set_k_mu(int run_){
	TLorentzVector k_mu;
	switch(run_){
		case 0:
		k_mu = _k_mu_e16_; //Constants.hpp
		break;
		case 1:
		k_mu = _k_mu_e1f_; //Constants.hpp
		break;
		default:
		k_mu = _k_mu_e16_;
		break;
	}
	return k_mu; 
}

float fun::W(TLorentzVector k_mu_prime_, int run_){
	TLorentzVector k_mu = fun::Set_k_mu(run_);
	TLorentzVector q_mu = k_mu - k_mu_prime_;
	return (_p_mu_ + q_mu).Mag();
}

TLorentzVector fun::Make_4Vector(float p, float theta, float phi, float m){
	TLorentzVector k_mu;
	TVector3 k_mu_3(p*TMath::Sin(theta*TMath::Pi()/180.0)*TMath::Cos(phi*TMath::Pi()/180.0), p*TMath::Sin(theta*TMath::Pi()/180.0)*TMath::Sin(phi*TMath::Pi()/180.0), p*TMath::Cos(theta*TMath::Pi()/180.0));
	k_mu.SetVectM(k_mu_3,m);
	return k_mu;
}
void fun::loadChain(std::shared_ptr<TChain> chain_, std::string file_, int thread_id_, int max_, std::shared_ptr<Flags> flags_){
  std::vector<std::string> filelist = fun::read_file_list(file_,thread_id_,flags_);//read_file_list(file); //creates a vector of file names
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
std::vector<std::string> fun::read_file_list(std::string path, int thread_num, std::shared_ptr<Flags> flags_){
  std::ifstream infile(path.c_str()); // in file stream
  std::vector<std::string> result;
  std::string line;
  int t = 0;
  while(getline(infile,line)) { //getline sees if there is a line available
    if(thread_num == (t%flags_->Flags::Num_Cores())){//_NUM_THREADS_)){
      result.push_back(line);//Gets the current line
    }
    t++;
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