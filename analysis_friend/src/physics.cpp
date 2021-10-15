#include "physics.hpp"

double physics::Virtual_Photon_Flux(double W_, double Q2_, double E_){
	double num1 = _alpha_em_ * W_ * (pow(W_,2.) - pow(_mp_,2.));
	double den1 = 4*TMath::Pi()*pow(E_,2.)*pow(_mp_,2.)*Q2_;
	double num2 = Q2_+pow(W_,2.) - pow(_mp_,2.);
	double den2 = 2*_mp_;
	double num3 = 2*(Q2_ + pow((num2/den2),2.));
	double den3 = 4*E_*(E_ - (num2/den2)) - Q2_;
	double bigden = 1-(1/(1+(num3/den3)));
	return num1/(den1*bigden);
}

double physics::Luminosity(double Qtot_, double corr_factor_){
	return corr_factor_*Qtot_*_length_target_*_density_target_*_avo_n_/(_qe_*_MH_);
}

double physics::Error(double N_gen_, double N_rec_, double weight_sum_){
	return sqrt((((N_gen_-2*N_rec_)/pow(N_gen_,3.0))*weight_sum_)+((pow(N_rec_,2.0)/pow(N_gen_,4.0))*weight_sum_));
}

double physics::Error(double N_gen_, double N_rec_){
	return sqrt((N_rec_*(N_gen_-N_rec_))/pow(N_gen_,3.0));
}

//double physics::Radiative_Corr();