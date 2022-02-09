#include "corrections.hpp"

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

float corr::p_corr_e(float p_e_, float theta_e_, float phi_e_, int run_, bool centered_, int sector_){
	float phi_const[4];
	if(centered_){
		for(int i=0; i<4; i++){
			for(int j=0; j<3; j++){
				phi_const[i] += _p_e_par[run_][sector_][i][j]*power(theta_e_,2-j);
			}
		}
		return fun::poly_3(phi_e_,phi_const[0],phi_const[1],phi_const[2],phi_const[3]);
	}
	int sector = fun::get_sector(phi_e_);
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++){
			phi_const[i] += _p_e_par[run_][sector][i][j]*power(theta_e_,2-j);
		}
	}
	float phi = fun::phi_center(phi_e_);
	return p_e_*fun::poly_3(phi,phi_const[0],phi_const[1],phi_const[2],phi_const[3]);
}

float corr::theta_e_corr(float theta_e_, float phi_e_, int run_, bool centered_, int sector_){
	//std::cout<<"Correcting Theta: " <<theta_e_ <<"\n"; 
	float phi_const[5]={0.0,0.0,0.0,0.0,0.0};
	if(centered_){
		for(int i=0; i<5; i++){
			for(int j=0; j<3; j++){
				//std::cout<<"Adding phi power (" <<i <<") theta power of " <<2-j <<" with angle par: " <<_angle_e_par[run_][sector_][i][j] <<" and result: " <<_angle_e_par[run_][sector_][i][j]*power(theta_e_,2-j) <<"\n";
				phi_const[i] += _angle_e_par[run_][sector_][i][j]*power(theta_e_,2-j);
				//std::cout<<"\tCurrent phi_const: " <<phi_const[i] <<"\n";
			}
		}
		//std::cout<<"\tPoly: " <<fun::poly_4(phi_e_,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]) <<" c0:" <<phi_const[0] <<" c1:" <<phi_const[1] <<" c2:" <<phi_const[2] <<" c3:" <<phi_const[3] <<" c4:" <<phi_const[4] <<"\n"; 
		//std::cout<<"\tCorrected Theta: " <<theta_e_-fun::poly_4(phi_e_,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]) <<"\n";
		return theta_e_-fun::poly_4(phi_e_,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]);
	}
	int sector = fun::get_sector(phi_e_);
	for(int i=0; i<5; i++){
		for(int j=0; j<3; j++){
			phi_const[i] += _angle_e_par[run_][sector][i][j]*pow(theta_e_,2-j);
		}
	}
	float phi = fun::phi_center(phi_e_);
	//std::cout<<"\tPoly: " <<fun::poly_4(phi_e_,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]) <<" c0:" <<phi_const[0] <<" c1:" <<phi_const[1] <<" c2:" <<phi_const[2] <<" c3:" <<phi_const[3] <<" c4:" <<phi_const[4] <<"\n";
	//std::cout<<"\tCorrected Theta: " <<theta_e_-fun::poly_4(phi,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]) <<"\n";
	return theta_e_-fun::poly_4(phi,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]);
}