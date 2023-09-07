#include "corrections.hpp"

float corr::power_10(float num_, int power_of_ten_){
	//std::cout<<"----Power of 10-----\n";
	//std::cout<<num_ <<"x10^" <<power_of_ten_ <<"\n";
	if(power_of_ten_==0){
		return num_; 
	}else if(power_of_ten_ > 0.0){
		for(int i=0; i<power_of_ten_; i++){
			num_ = num_*10.0;
		}
	}else if(power_of_ten_ < 0.0){
		for(int i=0; i<abs(power_of_ten_); i++){
			num_ = num_/10.0;
		}
	}
	//std::cout<<"gets us: " <<num_ <<"\n";
	return num_;
}

float corr::power(float num_, int power_){
	float out = num_;
	if(power_==0){
		return num_;
	}else if(power_>0){
		for(int i=0; i<power_; i++){
			out = out*num_;
		}
	}else if(power_<0){
		for(int i=0; i>power_; i--){
			out = out/num_;
		}
	}
	return out;
}

float corr::p_corr_e(float p_e_, float theta_e_, float phi_e_, int run_, bool centered_, int sector_idx_){//Will need exponential stuff
	float phi_const[4] = {0.0,0.0,0.0,0.0};
	if(centered_){
		for(int i=0; i<4; i++){
			for(int j=0; j<3; j++){
				//phi_const[i] += corr::power_10(_p_e_part[run_][sector_][i][j],_angle_e_expon[run_][sector_][i][j])*pow(theta_e_,2-j);
				//phi_const[i] += corr::power_10(_p_e_part[run_][sector_][i][j],_angle_e_expon[run_][sector_][i][j])*corr::power(theta_e_,2-j);
				phi_const[i] += _p_e_part[run_][sector_idx_][i][j]*corr::power(theta_e_,2-j);
			}
		}
		return fun::poly_3(phi_e_,phi_const[0],phi_const[1],phi_const[2],phi_const[3]);
	}
	int sector = physics::get_sector(phi_e_);
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++){
			//phi_const[i] += corr::power_10(_p_e_part[run_][sector][i][j],_angle_e_expon[run_][sector][i][j])*pow(theta_e_,2-j);
			//phi_const[i] += corr::power_10(_p_e_part[run_][sector][i][j],_angle_e_expon[run_][sector][i][j])*corr::power(theta_e_,2-j);
			phi_const[i] += _p_e_part[run_][sector-1][i][j]*corr::power(theta_e_,2-j);
		}
	}
	float phi = physics::phi_center(phi_e_);
	return p_e_*fun::poly_3(phi,phi_const[0],phi_const[1],phi_const[2],phi_const[3]);
}

float corr::theta_e_corr(float theta_e_, float phi_e_, int run_, bool centered_, int sector_){
	//std::cout<<"Correcting Theta:" <<theta_e_ <<"\n";
	if(theta_e_ >= 35.0 ){
		return theta_e_;
	}
	float phi_const[5] = {0.0,0.0,0.0,0.0,0.0};
	if(centered_){
		for(int i=0; i<5; i++){
			for(int j=0; j<3; j++){
				//std::cout<<"\t\tTheta Parameters: " <<_angle_e_part[run_][sector_][i][j] <<" x10^" <<_angle_e_expon[run_][sector_][i][j] <<" = " <<corr::power_10(_angle_e_part[run_][sector_][i][j],_angle_e_expon[run_][sector_][i][j]) <<"\n";
				//std::cout<<"\t\tTheta to power: " <<corr::power(theta_e_,2-j) <<"\n";
				//std::cout<<"\t\tPhi Constant added: " <<corr::power_10(_angle_e_part[run_][sector_][i][j],_angle_e_expon[run_][sector_][i][j])*corr::power(theta_e_,2-j) <<"\n";
				//phi_const[i] += corr::power_10(_angle_e_part[run_][sector_][i][j],_angle_e_expon[run_][sector_][i][j])*pow(theta_e_,2-j);
				//phi_const[i] += corr::power_10(_angle_e_part[run_][sector_][i][j],_angle_e_expon[run_][sector_][i][j])*corr::power(theta_e_,2-j);
				phi_const[i] += _angle_e_part[run_][sector_][i][j]*corr::power(theta_e_,2-j);
				//std::cout<<"\t\t\tPhi Constant " <<i <<" :" <<phi_const[i] <<"\n";
			}
		}
		//std::cout<<"\tPoly: " <<fun::poly_4(phi_e_,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]) <<" c0:" <<phi_const[0] <<" c1:" <<phi_const[1] <<" c2:" <<phi_const[2] <<" c3:" <<phi_const[3] <<" c4:" <<phi_const[4] <<"\n"; 
		//std::cout<<"\tCorrected Theta: " <<theta_e_-fun::poly_4(phi_e_,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]) <<"\n";
		return theta_e_-fun::poly_4(phi_e_,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]);
	}
	int sector = physics::get_sector(phi_e_);
	for(int i=0; i<5; i++){
		for(int j=0; j<3; j++){
			//std::cout<<"\t\tTheta Parameters: " <<_angle_e_part[run_][sector][i][j] <<" x10^" <<_angle_e_expon[run_][sector][i][j] <<" = " <<corr::power_10(_angle_e_part[run_][sector][i][j],_angle_e_expon[run_][sector][i][j]) <<"\n";
			//std::cout<<"\t\tTheta to power: " <<corr::power(theta_e_,2-j) <<"\n";
			//std::cout<<"\t\tPhi Constant added: " <<corr::power_10(_angle_e_part[run_][sector][i][j],_angle_e_expon[run_][sector][i][j])*corr::power(theta_e_,2-j) <<"\n";
			//phi_const[i] += corr::power_10(_angle_e_part[run_][sector][i][j],_angle_e_expon[run_][sector][i][j])*pow(theta_e_,2-j);
			//phi_const[i] += corr::power_10(_angle_e_part[run_][sector][i][j],_angle_e_expon[run_][sector][i][j])*corr::power(theta_e_,2-j);
			phi_const[i] += _angle_e_part[run_][sector][i][j]*corr::power(theta_e_,2-j);
			//std::cout<<"\t\t\tPhi Constant " <<i <<" :" <<phi_const[i] <<"\n";
		}
	}
	float phi = physics::phi_center(phi_e_);
	//std::cout<<"\tPoly: " <<fun::poly_4(phi,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]) <<" c0:" <<phi_const[0] <<" c1:" <<phi_const[1] <<" c2:" <<phi_const[2] <<" c3:" <<phi_const[3] <<" c4:" <<phi_const[4] <<"\n"; 
	//std::cout<<"\tCorrected Theta: " <<theta_e_-fun::poly_4(phi,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]) <<"\n";
	return theta_e_-fun::poly_4(phi,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]);
}