#include "physics.hpp"

double physics::Virtual_Photon_Flux(double W_, double Q2_, double E_){
	/*double num1 = _alpha_em_ * W_ * (pow(W_,2.) - pow(_mp_,2.));
	double den1 = 4*TMath::Pi()*pow(E_,2.)*pow(_mp_,2.)*Q2_;
	double num2 = Q2_+pow(W_,2.) - pow(_mp_,2.);
	double den2 = 2*_mp_;
	double num3 = 2*(Q2_ + pow((num2/den2),2.));
	double den3 = 4*E_*(E_ - (num2/den2)) - Q2_;
	double bigden = 1-(1/(1+(num3/den3)));
	return num1/(den1*bigden);*/
	double omega = (W_*W_+Q2_-_mp_*_mp_)/(2*_mp_);
	double lead = _alpha_em_/(4 *TMath::Pi()*E_*E_*_mp_*_mp_);
	double num1 = 2*(Q2_+omega*omega);
	double den1 = 4*E_*(E_-omega)-Q2_;
	double epsilon = 1.0/(1+(num1/den1));
	return lead*W_*(W_*W_-_mp_*_mp_)/((1-epsilon)*Q2_);
}

double physics::Luminosity(double Qtot_, double corr_factor_){
	return corr_factor_*Qtot_*_length_target_*_density_target_*_avo_n_/(_qe_*_MH_);
}

double physics::Error(double N_gen_, double N_rec_, double weight_sum_){
	//std::cout<<"N_gen = " <<N_gen_ <<" N_rec_ = " <<N_rec_ <<" weight_sum = " <<weight_sum_ <<"\n";
	if(N_gen_ > 0.0 && N_rec_ > 0.0 && weight_sum_ > 0.0){
		return sqrt((((N_gen_-2*N_rec_)/pow(N_gen_,3.0))*weight_sum_)+((pow(N_rec_,2.0)/pow(N_gen_,4.0))*weight_sum_));
	}else{
		return 0.0;
	}
}

double physics::Error(double N_gen_, double N_rec_){
	//std::cout<<"N_gen = " <<N_gen_ <<" N_rec_ = " <<N_rec_ <<"\n";
	if(N_gen_ > 0.0 && N_rec_ > 0.0){
		return sqrt((N_rec_*(N_gen_-N_rec_))/pow(N_gen_,3.0));
	}else{
		return 0.0;
	}
}

double physics::old_Virtual_Photon_Flux(double W_, double Q2_, double E_){
	double num1 = _alpha_em_ * W_ * (pow(W_,2.) - pow(_mp_,2.));
	double den1 = 4*TMath::Pi()*pow(E_,2.)*pow(_mp_,2.)*Q2_;
	double num2 = Q2_+pow(W_,2.) - pow(_mp_,2.);
	double den2 = 2*_mp_;
	double num3 = 2*(Q2_ + pow((num2/den2),2.));
	double den3 = 4*E_*(E_ - (num2/den2)) - Q2_;
	double bigden = 1-(1/(1+(num3/den3)));
	return num1/(den1*bigden);
}

double physics::ratio_Virtual_Flux(double W_, double Q2_, double E_){
	return old_Virtual_Photon_Flux(W_,Q2_,E_)/Virtual_Photon_Flux(W_,Q2_,E_);
}

//double physics::Radiative_Corr();