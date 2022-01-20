#include "corrections.hpp"

float corr:p_corr_e(float p_e_, float theta_e_, float phi_e_, int run_, bool centered_, int sector_){
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

float corr:theta_e_corr(float theta_e_, float phi_e_, int run_, bool centered_, int sector_){
	float phi_const[4];
	if(centered_){
		for(int i=0; i<5; i++){
			for(int j=0; j<3; j++){
				phi_const[i] += _angle_e_par[run_][sector_][i][j]*power(theta_e_,2-j);
			}
		}
		return fun::poly_4(phi_e_,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]);
	}
	int sector = fun::get_sector(phi_e_);
	for(int i=0; i<5; i++){
		for(int j=0; j<3; j++){
			phi_const[i] += _angle_e_par[run_][sector][i][j]*power(theta_e_,2-j);
		}
	}
	float phi = fun::phi_center(phi_e_);
	return theta_e_-fun::poly_4(phi,phi_const[0],phi_const[1],phi_const[2],phi_const[3],phi_const[4]);
}