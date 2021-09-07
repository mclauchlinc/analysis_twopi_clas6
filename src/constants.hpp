#ifndef CONST inline ANTS_HPP
#define CONST inline ANTS_HPP

#include <unordered_map>
#include "TMath.h"
#include "TLorentzVector.h"

//static int num_mixed_p_pip = 0; 

static const inline  int MAX_PARTS = 20; 
static const inline  int _NUM_THREADS_ = 6;//4;

static const inline  float _c_special_ = 29.9792458; //speed of light in cm/ns
static const inline  float _c_convert_ = 10000000; //Convert c_special to m/s

//Beam energies in GeV
static const inline  float _energy_e16_ = 5.754;//5.7696;5.754
//3_14 was 5.754, but led to bad results: updated to 5.759 after reading Paremuzyan's thesis 
//Check 4.794 GeV for e16
static const inline  float _energy_e1f_ = 5.499;

static const inline float _beam_energy_[] = {_energy_e16_,_energy_e1f_};

static const inline  float _fine_structure_ = 1.0/137.0; 

//Masses of relevant particles
static const inline 	float _me_ = 0.0005109989; //mass of electron in GeV
static const inline 	float _mp_ = 0.93828;	//Mass of proton in GeV
static const inline 	float _mpi_ = 0.1395;	//Mass of pion in GeV

static const inline  float _degree_ = 180.0/TMath::Pi();


//Initial Four Vectors

const inline  TLorentzVector _k_mu_e16_(0.0,0.0,sqrt(_energy_e16_*_energy_e16_-_me_*_me_),_energy_e16_);
const inline  TLorentzVector _k_mu_e1f_(0.0,0.0,sqrt(_energy_e1f_*_energy_e1f_-_me_*_me_),_energy_e1f_);
const inline  TLorentzVector _p_mu_(0.0,0.0,0.0,_mp_);

//e1-6 Luminosity Values
static const inline  float lt_e16 = 5.0; //Target length in cm
static const inline  float Dt_e16 = 0.073; //Density of target in g/cm^3
static const inline  float NA = 6.022 * pow(10.0,23); //Avogadro's number
static const inline  float qe = 1.602 * pow(10.0,-19); // fundamental Coulomb charge 
static const inline  float Mt_e16 = 1.007; //Molar mass of target in g/mole
static const inline  float Qt_e16 = 21.32*pow(10.0,-3); //Total charge incident on target from Arjun in Coulombs //Will need to verify myself
static const inline  float L_e16 = Qt_e16*lt_e16*Dt_e16*NA/(qe*Mt_e16);
static const inline  float cm2_to_mbarn = pow(10.0,-27);

//Particle ID Numbers
//anti particles are negative
static const inline  	int _PROTON_ = 2212;
static const inline 	int _ELECTRON_ = 11;
static const inline 	int _PION_ = 211;
static const inline 	int _PION_0_ = 111;

// Beam spot from Arjun's Fortran vertex_e1f.f
static const inline  float _x_beam_ = 0.090;
static const inline  float _y_beam_ = -0.345;

//Fiducial Cut Parameters (electrons have e, hadrons have h)
//Arjun's cut_fid_e1f.f

static const inline  float pim_bot_MM[11] = {0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.083,0.109,0.109};
static const inline  float pim_top_MM[11] = {0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.34,0.415,0.415};


static const inline  float _XYmax_[3] = {545.0,545.0,545.0};
static const inline  int _XYres_ = 350;

//CC Projection Plane
static const inline  float _Acc_ = -0.000785;
static const inline  float _Bcc_ = 0.0;
static const inline  float _Ccc_ = -0.00168;
static const inline  float _Dcc_ = 1.0;


//binning
static const inline  float Wbin_res = 0.025;//The width of a W bin //30 steps
static const inline  float Wbin_start = 1.4;//The starting of W bins

static const inline  float Q2bin_res = 0.5;//6 steps
static const inline  float Q2bin_start = 2.0; 

static const inline  float pbin_res = 0.18;//Range: 0-5.0
static const inline  float pbin_start = 0.5;

//For cross sections
static const inline  Int_t xc_bins[7]= {29,5,14,14,10,10,10};


//Paths for file names
static const inline  std::string path1 = "/home/mclauchlinc/Desktop/analysis/nick_convert_e16.txt";
static const inline  std::string path2 = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/e16_10_18_17_ntple.txt";
static const inline  std::string path3 = "/Users/cmc/Desktop/analysis/Nick_skim_e16.txt";
static const inline  std::string path4 = "/Users/cmc/Desktop/analysis/arjun_sim.txt";
static const inline  std::string path3p = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/NickSkim_e16_PlateIN.txt";
static const inline  std::string path3n = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/NickSkim_e16_PlateOUT.txt";
static const inline  std::string paths1 = "/Users/cmc/Desktop/analysis/simulation/sim_e16_group1.txt";
static const inline  std::string paths2 = "/Users/cmc/Desktop/analysis/simulation/sim_e16_pre_gpp.txt";
static const inline  std::string paths3n = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gsim_no_gpp_nmcdata.txt";
static const inline  std::string paths3y = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gsim_no_gpp_mcdata.txt";
static const inline  std::string paths3 = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/sim_e16.txt";//"/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/sim_e16_fix.txt";
static const inline  std::string paths32 = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/sim_e16_2.txt";
static const inline  std::string paths_ar = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/arjun_sim_e16.txt";
static const inline  std::string path_ts = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/test_sim_e16.txt";
static const inline  std::string path_ts_gpp = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gpp_test_sim_e16.txt";
static const inline  std::string path_ts_gsim = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gsim_test_sim_e16.txt";
static const inline  std::string path_ts_gsim_1 = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gsim_test_e16_1.txt";
static const inline  std::string paths_ar_a1c = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/arjun_recooked_gsim_a1c.txt";
static const inline  std::string paths_ar_ana = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/arjun_recooked_gsim_ana.txt";
static const inline  std::string path_ce16 = "/home/mclauchc/analysis/Path_Files/cluster_e16_path.txt";
static const inline  std::string path_ce16s = "/home/mclauchc/analysis/Path_Files/cluster_e16_sim_files.txt";
static const inline  std::string path_ce1f = "/home/mclauchc/analysis/Path_Files/cluster_e1f_file_paths.txt";
static const inline  std::string path_ce1fs = "/home/mclauchc/analysis/Path_Files/cluster_e1f_sim_files.txt";

static const inline  int _plate_swap_e16_[] = {1000,1002};
static const inline  int _plate_swap_e1f_[] = {1000,1002};

static const inline  int _e16_ = 0;
static const inline  int _e1f_ = 1; 

//Cut Names
static const inline  char * _none_ = "none";
static const inline  char * _sanity_ = "sanity";
static const inline  char * _fid_cut_ = "fid";
static const inline  char * _delta_cut_ = "delta";
static const inline  char * _sf_cut_ = "sf";
static const inline  char * _cc_cut_ = "min_cc";
static const inline  char * _ec_cut_ = "min_ec";
static const inline  char * _vertex_cut_ = "vertex";
static const inline  char * _beta_cut_ = "beta";
static const inline  char * _pid_ = "pid";
static const inline  char * _event_ = "event";
static const inline  char * _ele_ = "ele";
static const inline  char * _pro_ = "pro";
static const inline  char * _pip_ = "pip";
static const inline  char * _pim_ = "pim";
//Event Indexes
static const inline  char * _mpro_ = "mpro";
static const inline  char * _mpip_ = "mpip";
static const inline  char * _mpim_ = "mpim";
static const inline  char * _mzero_ = "mzero";
static const inline  char * _mall_ = "mall";
static const inline  char * _mnone_ = "mnone";

static const inline  char * _nweighted_ = "no-weight";
static const inline  char * _weighted_ = "weight";

static const inline  char * _cut_applied_ = "cut";
static const inline  char * _anti_cut_ = "anti";
static const inline  char * _no_cut_ = "no-cut";

static const inline  char * _nthrown_ = "recon";
static const inline  char * _thrown_ = "thrown";

static const inline  char * _dirty_ = "dirty";
static const inline  char * _clean_ = "clean";
static const inline  char * _isolated_ = "isolated";

static const inline  char * _W_var_ = "W_Dep";
static const inline  char * _W_range_ = "in_range";
static const inline  char * _W_all_ = "all";
static const inline  char * _W_dep_[] = {_W_var_,_W_range_,_W_all_};


static const inline  char * _true_ = "true";
static const inline  char * _false_ = "false";
static const inline  char * _truth_[] = {_false_,_true_};

//Parsing Out PID Cuts
static const inline  char* _ecuts_[] = {_none_,_sanity_, _fid_cut_, _sf_cut_, _cc_cut_, _ec_cut_, _vertex_cut_,_pid_,_event_};//Be sure to add statements for flags if modified
static const inline  char* _hcuts_[] = {_none_,_sanity_, _fid_cut_, _delta_cut_,_pid_,_event_};//Be sure to add statements for flags if modified
static const inline  char* _top_[] = {_mpro_,_mpip_,_mpim_,_mzero_,_mall_,_mnone_};
static const inline  char* _cut_[] = {_cut_applied_,_anti_cut_,_no_cut_};
static const inline  char* _species_[] = {_ele_,_pro_,_pip_,_pim_};
static const inline  char* _recon_[] = {_nthrown_,_thrown_};
static const inline  char* _weight_[] = {_nweighted_,_weighted_};
static const inline  char* _clean_event_[] = {_dirty_,_clean_,_isolated_};

static const inline  char* _sec1_ = "sec1";
static const inline  char* _sec2_ = "sec2";
static const inline  char* _sec3_ = "sec3";
static const inline  char* _sec4_ = "sec4";
static const inline  char* _sec5_ = "sec5";
static const inline  char* _sec6_ = "sec6";
static const inline  char* _sec_all_ = "all";
static const inline  char* _sector_[] = {_sec1_,_sec2_,_sec3_,_sec4_,_sec5_,_sec6_,_sec_all_};

static const inline  char* _left_ = "left";
static const inline  char* _right_ = "right";
static const inline  char* _coinc_ = "coinc";
static const inline  char* _cc_sides_[] = {_left_,_right_,_coinc_};



static const inline  char * species[] = {_ele_,_pro_,_pip_,_pim_};//4
static const inline  char * eid_cut[] = {"pre","sanity","fid","sf","min_cc","fid+sf","fid+cc","sf+cc","eid","bank","event"}; //11"min_cc","min_ec","eid","bank","event"};//"fid+sf","fid+cc","sf+cc","eid","bank","event"}; //11
static const inline  char * cut_ver[] = {"cut","anticut"};
static const inline  char * hid_cut[] = {"pre","sanity","fid","dt","hid","bank","event"};//,"pWQ2"}; //7
static const inline  char * topologies[] = {"None","Pmiss","PIPmiss","PIMmiss","Zeromiss","ALLmiss"}; //6
static const inline  char * sec_list[] = {"all_sectors","sec1","sec2","sec3","sec4","sec5","sec6"}; //7`
static const inline  char * W_dep_list[] = {"No_W_Dep","W_Dep"};
static const inline  char * CC_det_side[] = {"Left","Coince","Right","All"};
static const inline  char * basic_cut[] = {"pre","cut","anti"};
static const inline  char * MM_sq[] = {"linear","squared"};
static const inline  char * par_cut[4][11] = {{"pre","sanity","fid","sf","min_cc","fid+sf","fid+cc","sf+cc","eid","bank","event"},{"pre","sanity","fid","dt","hid","bank","event"},{"pre","sanity","fid","dt","hid","bank","event"},{"pre","sanity","fid","dt","hid","bank","event"}};
static const inline  char * fit_q[] = {"4fit","4show"};
static const inline  char * throw_stat[] = {"recon","thrown"};//For reconst inline ructed vs. thrown events
static const inline  char * w_stat[] = {"nweight","weight"};
static const inline  char * detectors[] = {"dc","cc","sc","ec"};
static const inline  char * inter_had[] = {"#Delta^{++}","#rho","#Delta^{0}"};

#endif
