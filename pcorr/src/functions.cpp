#include "functions.hpp"

float fun::theta_calc(float theta_p_, float beam_energy_, bool radians_){
	float theta = theta_p_*_radian_;
	float denom = (beam_energy_ + _mp_)*Tan(theta);
	if(radians_){
		return 2*ATan(_mp_/denom);
	}else{
		return 2*ATan(_mp_/denom)*_degree_;
	}
}

float delta_theta(float theta_e_, float theta_p_, float beam_energy_){
	return fun::theta_calc(theta_p_,beam_energy_) - theta_e_;
}

float fun::poly_4(float x_, float a_, float b_, float c_, float d_, float e_){
	return a_*power(x_,4) + b_*power(x_,3) + c_*power(x_,2) + d_*x_ + e_;
}

float fun::poly_3(float x_, float a_, float b_, float c_, float d_){
	return a_*power(x_,3) + b_*power(x_,2) + c_*x_ + d_;
}

float fun::poly_2(float x_, float a_, float b_, float c_){
	return a_*power(x_,2) + b_*x_,2 + c_;
}

float fun::theta(std::shared_ptr<Branches> data_, int idx_, bool radians_){
	if(radians_){
		return TMath::ACos(cz_)*_degree_;
	}
	return TMath::ACos(cz_)*_degree_;
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
	return get_sector(phi(data_),idx_);
}


float fun::phi_center( float phi_)
{
	double phi_corr;
	int sector = get_sector(phi_);
	switch(sector_){
		case 1:
			phi_corr = phi_;
		break;
		case 2:
			phi_corr = phi_-60;
		break;
		case 3:
			phi_corr = phi_-120;
		break;
		case 4:
			if(phi_<=-150){
				phi_corr = phi_+180;
			}else if(phi_>=150){
				phi_corr = phi_-180;
			}
		break;
		case 5:
			phi_corr = phi_+120;
		break;
		case 6:
			phi_corr = phi_+60;
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


float fun::phi(std::shared_ptr<Branches> data_, int idx_, bool center_, bool radians_){
	float phi = TMath::ATan2(cy_,cx_);
	if(center_){
		phi = phi_center(phi*_degree_);
	}
	if(radians_){
		return phi;
	}
	return phi*_degree_;
}

float fun::p_calc(float theta_e_, float beam_energy_){
	float denom = 1+ (2*beam_energy_ * power(TMath::Sin(theta_e_*_radian_/2),2))/_mp_;
	return beam_energy_/denom;
}

float fun::delta_p_e(float p_e_, float theta_e_, float beam_energy_){
	return fun::p_calc(theta_e_, beam_energy_) - p_e_;
}

float fun::W(TLorentzVector k_mu_prime_, int run_){
	TLorentzVector k_mu = physics::Set_k_mu(run_);
	TLorentzVector q_mu = k_mu - k_mu_prime_;
	return (_p_mu_ + q_mu).Mag();
}

TLorentzVector fun::Make_4Vector(float p, float theta, float phi, float m){
	TLorentzVector k_mu;
	TVector3 k_mu_3(p*TMath::Sin(theta*TMath::Pi()/180.0)*TMath::Cos(phi*TMath::Pi()/180.0), p*TMath::Sin(theta*TMath::Pi()/180.0)*TMath::Sin(phi*TMath::Pi()/180.0), p*TMath::Cos(theta*TMath::Pi()/180.0));
	k_mu.SetVectM(k_mu_3,m);
	return k_mu;
}