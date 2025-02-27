#include "detectors.hpp"

float detect::x_det_center(float x_,float y_,int sector_){
	float rot_angle = _detector_center_angles_[sector_-1]*TMath::Pi()/180.0;
	return TMath::Cos(rot_angle)*x_-TMath::Sin(rot_angle)*y_;
}
float detect::y_det_center(float x_, float y_,int sector_){
	float rot_angle = _detector_center_angles_[sector_-1]*TMath::Pi()/180.0;
	return TMath::Sin(rot_angle)*x_+TMath::Cos(rot_angle)*y_;
}
float detect::phi_det_center(float phi_,int sector_){

}
float detect::theta_det_center(float theta_,int sector_){

}

int detect::cc_segment(int cc_segm){
	int seg = -1;

	if(cc_segm > 0 && cc_segm <200){
		seg = ((cc_segm)/10)-1;
		//std::cout<<"cc_seg: " <<cc_segm <<" -> " <<seg <<"             " <<std::endl;
	}
	if(cc_segm > 1000 && cc_segm <1200){
		seg = ((cc_segm-1000)/10)-1;
		//std::cout<<"cc_seg: " <<cc_segm <<" -> " <<seg<<"             " <<std::endl;
	}
	if(cc_segm > 2000 && cc_segm <2200){
		seg = ((cc_segm-2000)/10)-1;
		//std::cout<<"cc_seg: " <<cc_segm <<" -> " <<seg<<"             " <<std::endl;
	}
	return seg;
}

//Defining left, right, or coincidence hits in the CC
int detect::cc_lrc(int cc_segm){
	int po = -1;
	if(cc_segm > 0 && cc_segm <200){
		po = 0;//right
	}
	if(cc_segm > 1000 && cc_segm <1200){
		po = 1;//coincident
	}
	if(cc_segm > 2000 && cc_segm <2200){
		po = 2;//left
	}
	
	return po;
}

TVector3 detect::cc_sector_center(TVector3 v_){
	float _phi_ = physics::get_phi(v_[0],v_[1]);
	int _sector_ = physics::get_sector(_phi_);
	float _x_ = NAN;
	float _y_ = NAN;
	float _z_ = v_[3];
	_x_ = physics::CCX_Rotate(v_[0],v_[1],_sector_);
	_y_ = physics::CCY_Rotate(v_[0],v_[1],_sector_);
	TVector3 _v_ = {_x_,_y_,_z_}; 
	return _v_;
}
TVector3 detect::cc_sector_center(float x_, float y_, float z_){
	float _phi_ = physics::get_phi(x_,y_);
	int _sector_ = physics::get_sector(_phi_);
	float _x_ = NAN;
	float _y_ = NAN;
	float _z_ = z_;
	_x_ = physics::CCX_Rotate(x_,y_,_sector_);
	_y_ = physics::CCY_Rotate(x_,y_,_sector_);
	TVector3 _v_ = {_x_,_y_,_z_}; 
	return _v_;
}

float detect::cc_theta(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_){
	TVector3 n(cx_sc_, cy_sc_, cz_sc_);
	TVector3 p0(x_sc_, y_sc_, z_sc_);
	TVector3 s(_Acc_,_Bcc_,_Ccc_);
	//TVector3 s(physics::X_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),physics::Y_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),Ccc);//Get Projective plane for particular secgtor we're working in
	float tmag = TMath::Abs((physics::Dot_Product(s,p0)+_Dcc_)/physics::Dot_Product(s,n)); 
	TVector3 t = tmag*n;
	TVector3 _p_ = p0 + t; 
	return TMath::ACos(_p_[2]/physics::Vec3_Mag(_p_))*(180.0/TMath::Pi());
}

float detect::cc_theta(std::shared_ptr<Branches> data_, int idx_){
	return detect::cc_theta(data_->dc_cxsc(idx_), data_->dc_cysc(idx_), data_->dc_czsc(idx_), data_->dc_xsc(idx_), data_->dc_ysc(idx_), data_->dc_zsc(idx_));
}

float detect::cc_phi(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_){
	TVector3 n(cx_sc_, cy_sc_, cz_sc_);
	TVector3 p0(x_sc_, y_sc_, z_sc_);
	TVector3 s(_Acc_,_Bcc_,_Ccc_);
	//TVector3 s(physics::X_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),physics::Y_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),Ccc);//Get Projective plane for particular secgtor we're working in
	float tmag = TMath::Abs((physics::Dot_Product(s,p0)+_Dcc_)/physics::Dot_Product(s,n)); 
	TVector3 t = tmag*n;
	TVector3 _p_ = p0 + t; 
	return TMath::ATan2(_p_[1],_p_[0])*(180.0/TMath::Pi());
}

float detect::cc_phi(std::shared_ptr<Branches> data_, int idx_){
	return detect::cc_phi(data_->dc_cxsc(idx_), data_->dc_cysc(idx_), data_->dc_czsc(idx_), data_->dc_xsc(idx_), data_->dc_ysc(idx_), data_->dc_zsc(idx_));
}

float detect::cc_x(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_, int sec_){
	TVector3 n(cx_sc_, cy_sc_, cz_sc_);
	TVector3 p0(x_sc_, y_sc_, z_sc_);
	TVector3 s(_Acc_,_Bcc_,_Ccc_);
	//TVector3 s(physics::X_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),physics::Y_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),Ccc);//Get Projective plane for particular secgtor we're working in
	float tmag = TMath::Abs((physics::Dot_Product(s,p0)+_Dcc_)/physics::Dot_Product(s,n)); 
	TVector3 t = tmag*n;
	TVector3 _p_ = p0 + t; 
	return physics::X_Rotate(_p_[0],_p_[1],sec_);
}

float detect::cc_y(float cx_sc_, float cy_sc_, float cz_sc_, float x_sc_, float y_sc_, float z_sc_, int sec_){
	TVector3 n(cx_sc_, cy_sc_, cz_sc_);
	TVector3 p0(x_sc_, y_sc_, z_sc_);
	TVector3 s(_Acc_,_Bcc_,_Ccc_);
	//TVector3 s(physics::X_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),physics::Y_Rotate(Acc,Bcc,physics::get_sector(physics::get_phi(cx_sc_,cy_sc_))),Ccc);//Get Projective plane for particular secgtor we're working in
	float tmag = TMath::Abs((physics::Dot_Product(s,p0)+_Dcc_)/physics::Dot_Product(s,n)); 
	TVector3 t = tmag*n;
	TVector3 _p_ = p0 + t; 
	return physics::Y_Rotate(_p_[0],_p_[1],sec_);
}

float detect::cc_x(std::shared_ptr<Branches> data_, int idx_){
	return detect::cc_x(data_->dc_cxsc(idx_), data_->dc_cysc(idx_), data_->dc_czsc(idx_), data_->dc_xsc(idx_), data_->dc_ysc(idx_), data_->dc_zsc(idx_), data_->sc_sect(idx_));
}

float detect::cc_y(std::shared_ptr<Branches> data_, int idx_){
	return detect::cc_y(data_->dc_cxsc(idx_), data_->dc_cysc(idx_), data_->dc_czsc(idx_), data_->dc_xsc(idx_), data_->dc_ysc(idx_), data_->dc_zsc(idx_), data_->sc_sect(idx_)); 
}

float detect::cc_eff(int run_, int sec_, int seg_, int pos_){
	return detect::cc_eff_vals[run_][sec_-1][pos_][seg_][1];//Right now just mid cut
}


//SC functions
float detect::sc_theta(std::shared_ptr<Branches> data_, int idx_){
	return physics::get_theta(data_->Branches::dc_czsc(idx_)); 
}

float detect::sc_phi(std::shared_ptr<Branches> data_, int idx_){
	return physics::get_phi(data_->dc_xsc(idx_),data_->dc_ysc(idx_));
}


//EC functions
float detect::ec_theta(std::shared_ptr<Branches> data_, int idx_){
	TVector3 _p_(data_->Branches::ech_x(idx_),data_->Branches::ech_y(idx_),data_->Branches::ech_z(idx_));
	float pmag = physics::Vec3_Mag(_p_);
	return physics::get_theta(data_->Branches::ech_z(idx_)/pmag); 
}

float detect::ec_phi(std::shared_ptr<Branches> data_, int idx_){
	float _phi_ = physics::get_phi(data_->ech_x(idx_),data_->ech_y(idx_));
	return physics::phi_center(_phi_);
}




