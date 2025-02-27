#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD

#include <unordered_map>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLatex.h"


static const double _alpha_em_=1./137.035999084; //electromag coupling constant
static const float _mp_ = 0.93828;//mass of proton
static const float _mpi_ = 0.1395;

//Beam energies in GeV
static const float _energy_e16_ = 5.754;
static const float _energy_e1f_ = 5.499;
static const float _beam_energy_[] = {_energy_e16_,_energy_e1f_};

//Beam Polarizations
static const float _pol_e16_ = 0.70;//Not accurate and needs to be updated 12/6/22
static const float _pol_e1f_ = 0.70;//Not accurate and needs to be updated 12/6/22
static const float _beam_pol_[] = {_pol_e16_,_pol_e1f_};

//Cut Values
static const double _Acceptance_Rel_Error_Max_ = 0.75;

//Target info for luminosity
static const double _length_target_ = 5.; //cm
static const double _density_target_ = 0.073; //g/cm^3 
static const double _avo_n_ = 6.022*pow(10.0,-23);// mol^-1 avogadro's #
static const double _qe_ = 1.602*pow(10.,-19.); //Coulombs charge of electron
static const double _MH_ = 1.007; //g/mol molar mass of hydrogen
static const double _Q_tot_ = 21.32*pow(10.,-3.); //Coulombs from Arjun's Thesis. need to calculate myself 

static const char * _sparse_names_[] = {"2#pi_off_proton_#Delta^{++}","2#pi_off_proton_#rho","2#pi_off_proton_#Delta^{0}"};
static const  char * _mpro_ = "mpro";
static const  char * _mpip_ = "mpip";
static const  char * _mpim_ = "mpim";
static const  char * _mzero_ = "mzero";
static const  char * _mall_ = "mall";
static const  char* _top_[] = {_mpro_,_mpip_,_mpim_,_mzero_,_mall_};
static const char * _topo_[] = {"Pmiss","PIPmiss","PIMmiss","Zeromiss","ALLmiss"};
static const char * _var_set_[] = {"PIM","Pro","PIP"};
static const bool _cluster_loc_[] = {false,true};
static const char * _five_dim_[] = {"MM1","MM2","Theta","Alpha","Phi"};
static const char * _cross_five_dim_[] = {"MM1","MM2","Theta)","Alpha","Phi"};
static const char * _dim_units_[] = {"GeV","GeV","Angle (deg)","Angle (deg)","Angle (deg)"};
static const char * _dim_units_y_[] = {"GeV","GeV","rad","rad","rad"};

//const int * some_bins[] = {{0,4},{1,4},{2,4},{3,4}};

#endif