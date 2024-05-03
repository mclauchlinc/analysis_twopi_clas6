#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "physics.hpp"
#include "branches.hpp"
#include "cuts.hpp"
//#include "event_class.hpp"
#include "histogram.hpp"

#include "physics.hpp"
#include "pid.hpp"
#include "flags.hpp"
#include "corrections.hpp"
#include "detectors.hpp"
/*
	This class is for individual particles in an event
	It's filled with all the relevant kinematic, detector, and cut info
	to establish what the particle is and allow for easier 
	management of the data. 

	For plotting, a Particle is fed to a Plotting Function in
	plotting.hpp
*/


class Particle{
private:
	int _idx = -1; 
	int _run = -1; //{0,1} -> {e16,e1f}
	bool _sim = false;
	bool _thrown = false;
	float _weight = NAN;
	float _cc_eff = NAN;

	float _p = NAN;//In lab frame
	int _q = 0; 
	float _dt[4] = {NAN,NAN,NAN,NAN};
	float _theta = NAN;//In lab frame
	float _phi = NAN;//In lab frame
	float _vz = NAN;//In lab frame
	float _vx = NAN;//In lab frame
	float _vy = NAN;//In lab frame
	float _sf = NAN;//In lab frame
	float _etot = NAN;//Energy deposited in EC
	int _cc_seg = -1; //Segment of CC hit
	int _cc_sect = -1; //Sector of CC hit
	int _cc_lrc = -1; //Side of CC hit {0,1,2} -> {left, coic, right}
	int _nphe = -1;//number photo electrons generated in CC


	bool _sanity_pass[4] = {false,false,false,false}; //Pass Sanity Check {ele,pro,pip,pim}
	bool _min_ec_pass = false;//Pass Min EC Check
	bool _fid_pass[4] = {false,false,false,false}; //Pass Fiducial check {ele,pro,pip,pim}
	bool _sf_pass = false;//Pass Sampling Fraction
	bool _cc_pass = false; //Pass Min CC 
	bool _dt_pass[4] = {false,false,false,false};//Pass Delta t {ele,pro,pip,pim}
	bool _beta_pass[4] = {false,false,false,false};//Pass Beta {ele,pro,pip,pim}
	bool _id_pass[4] = {false,false,false,false};//Pass ID bank {ele,pro,pip,pim}
	bool _sp_dt_pass[3] = {false,false, false};
	bool _vertex_pass = false;//Did it fall within the determined vertex region?
	bool _p_corr = false; //Performed Momentum Correction
	int _id_crisis = 0;//{0,1,2} = {none, pro/pip, pim/e}
	//Efficiency Passes
	bool _sc_eff_pass[4] = {false,false,false,false};
	std::vector<bool> _pid;// = {false,false,false,false};//Identified {ele,pro,pip,pim}
	bool _ided = false; //Has gone through PID
	bool _event[4] = {false,false,false,false};

	float _x[3] = {NAN,NAN,NAN};//cc,sc,ec 
	float _y[3] = {NAN,NAN,NAN};//cc,sc,ec
	float _dtheta[3] = {NAN,NAN,NAN};
	float _dphi[3] = {NAN,NAN,NAN};

	int _sc_pd = -1;//sc_paddle


public:
	Particle(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_, bool thrown_ = false);
	void PID_Thrown(int idx_, std::shared_ptr<Branches> data_);
	void PID_Recon(int idx_, std::shared_ptr<Branches> data_, std::shared_ptr<Flags> flags_);
	
	//Pass checks for all relevant pieces
	bool Pass_Sanity(int i);
	bool Pass_ec();
	bool Pass_fid(int i);
	bool Pass_sf();
	bool Pass_cc();
	bool Pass_dt(int i);
	bool Pass_id(int i);
	bool Pass_pid(int i);
	bool Pass_vertex();
	bool Pass_SC_Eff(int i);
	bool Corr_p();
	int ID_crisis();
	bool IDed();
	bool Is_Sim();
	bool Is_Thrown();
	bool Is_Elec();
	bool Is_Pro();
	bool Is_Pip();
	bool Is_Pim();
	float Get_p();
	float Get_theta();
	float Get_phi();
	int Get_run();
	int Get_idx();
	float Get_Weight();
	int Get_q();
	int Get_sc_pd();

	float W();
	float Q2();

	int Sector();

	float Get_sf();
	float Get_etot();
	int Get_cc_seg();
	int Get_cc_lrc();
	int Get_nphe();
	float Get_vz();
	float Get_vx();
	float Get_vy();
	float Get_delta(int par_);
	float Get_x(int det_);
	float Get_y(int det_);

	TLorentzVector Get_4Vec(int i);//Which assumed mass (needed due to dual id of proton and pip)

	//void Fill_Par_Event(std::shared_ptr<Environment> envi_, std::shared_ptr<Histogram> hist_, float W_, int top_, int par_, bool pass_);

	//void Check_Particle();
	
};


#endif