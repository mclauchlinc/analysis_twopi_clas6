#ifndef EVENT_HPP
#define EVENT_HPP

#include "histogram.hpp"
#include "particle.hpp"
#include "cuts.hpp"
#include "physics.hpp"
#include "constants.hpp"
#include "branches.hpp"


class Event{
private:
	//Topology
	int _run = -1; //{0,1} -> {e16,e1f}
	bool _top[4] = {false,false,false,false}; //Which topology was looked at {pmiss,pipmiss,pimmiss,zeromiss}
	bool _pass[4] = {false,false,false,false}; //Did it pass this topology's MM cut? 
	int _pass_top = -1; //Final assignment of Event Topology {zero,pim,pip,pro}
	float _weight = 1.0; 
	float _ev_weight = NAN;
	bool _sim = false; 
	bool _thrown = false;
	bool _filled_correctly = false; //Was this event filled correctly?
	int _hel = 0;
	bool _COM = false;
	bool _full_event = false;
	int _set = -1; 
	float _MM = -99.9;//For exclusive topology do not want 0.0 coming back unnecessarily 
	float _MM2 = -99.9;//For exclusive topology do not want 0.0 coming back unnecessarily 
	float _error = NAN;

	bool _part[4] = {false,false,false,false}; //Which particles are present?

	//Particle Attributes {ele',pro',pip',pim'}
	float _p_lab[4] = {NAN,NAN,NAN,NAN};//Particle momentum in lab frame
	float _theta_lab[4] = {NAN,NAN,NAN,NAN};//Particle theta  in lab frame
	float _phi_lab[4] = {NAN,NAN,NAN,NAN}; //Particle phi in lab frame
	float _p[4] = {NAN,NAN,NAN,NAN};//Particle momentum in COM frame
	float _theta[4] = {NAN,NAN,NAN,NAN};//Particle theta in COM frame
	float _phi[4] = {NAN,NAN,NAN,NAN}; //Particle phi in COM frame
	float _sf = NAN;//Sampling Fraction of Electron;
	float _etot = NAN; //Energy deposited in EC by electron
	int _nphe = -99; //Number photo electrons put in CC
	int _cc_seg = -99; //CC segment
	int _cc_lrc = -99; //CC side {left,coinc,right}
	float _cc_eff = NAN;//CC Efficiency
	float _dt[4] = {NAN,NAN,NAN,NAN};
	int _sc_pd[4] = {-1,-1,-1,-1};
	
	//Virtual Photon Flux
	float _virtual_photon_flux = NAN;

	float _vz = NAN;
	float _vx = NAN;
	float _vy = NAN;

	TLorentzVector _k1_lab = {NAN,NAN,NAN,NAN};
	//TLorentzVector _k1_lab_e16 = physics::Make_4Vector(_energy_e16_,0.0,0.0,1.0,_me_);//{NAN,NAN,NAN,NAN};
	//TLorentzVector _k1_lab_e1f = physics::Make_4Vector(_energy_e1f_,0.0,0.0,1.0,_me_);//{NAN,NAN,NAN,NAN};
	//TLorentzVector _k1_lab[] = {_k1_lab_e16,_k1_lab_e1f};
	TLorentzVector _k1 = {NAN,NAN,NAN,NAN};//beam
	TLorentzVector _p1 = {NAN,NAN,NAN,NAN};//target proton
	TLorentzVector _vec_lab[4] = {{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN}};//4 vectors in lab frame
	TLorentzVector _vec[4] = {{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN},{NAN,NAN,NAN,NAN}};//4 vectors in COM frame

	//Event Binning
	float _W = NAN; 
	float _Q2 = NAN; 
	float _MMb[3] = {NAN,NAN,NAN};//Combined missing mass for given expression of data 
	float _MM2b[3] = {NAN,NAN,NAN};//Combined missing mass for given expression of data 
	float _thetab[3] = {NAN,NAN,NAN};
	float _alphab[3] = {NAN,NAN,NAN};
	float _phib[3] = {NAN,NAN,NAN};


public:
	Event(int top_, Particle p0_, Particle p1_, Particle p2_, std::shared_ptr<Flags> flags_, float weight_, int hel_, bool thrown_ = false);
	Event(int top_, Particle p0_, Particle p1_, Particle p2_, Particle p3_, std::shared_ptr<Flags> flags_, float weight_, int hel_, bool thrown_ = false);
	bool Pass_Top(int i);
	bool Pass();
	int Top();
	bool Check_Particles(int top_, Particle p0_, Particle p1_, Particle p2_, std::shared_ptr<Flags> flags_);
	bool Check_Particles(int top_, Particle p0_, Particle p1_, Particle p2_, Particle p3_, std::shared_ptr<Flags> flags_);
	void Extract_Particles(int top_, Particle p0_, Particle p1_, Particle p2_, std::shared_ptr<Flags> flags_);
	void Extract_Particles(int top_, Particle p0_, Particle p1_, Particle p2_, Particle p3_, std::shared_ptr<Flags> flags_);
	void Calc_Error();
	float Error();
	float Weight();
	float MM();
	float MM2();
	float Theta(int particle_);
	float Phi(int particle_);
	float W();
	float Q2();
	float SF();
	float P(int particle_, bool COM_=false);
	float Delta(int particle_);
	float Vz();
	float Vx();
	float Vy();
	int CC_seg();
	int CC_side();
	float CC_eff();
	int nphe();
	float MMb(int i);
	float MM2b(int i);
	float Thetab(int i);
	float Phib(int i);
	float Alphab(int i);
	void Missing_Hadron();
	void COM();
	void Vars();
	int Helicity();
	float Get_Px(int particle_, bool COM_);
	float Get_Py(int particle_, bool COM_);
	float Get_Pz(int particle_, bool COM_);
	float Get_P0(int particle_, bool COM_);
	float Get_Beam_Comp(int component_, bool COM_);
	float Get_Target_Comp(int component_, bool COM_);
	int Get_PID(int particle_);
	void Get_Angles(bool COM_=false);
	bool Was_COM();
	int Run();
	int Sector(int particle_);
	float Virtual_Photon_Flux();
	int SC_pd(int i);
	void Check_Event(bool thrown_, std::shared_ptr<Flags> flags_);
};



#endif