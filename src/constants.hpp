#ifndef CONST ANTS_HPP
#define CONST ANTS_HPP

#include <unordered_map>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLatex.h"

//static int num_mixed_p_pip = 0; 

static const  int MAX_PARTS = 20; 
static const  int _NUM_THREADS_ = 6;//20;//4;

static const  float _c_special_ = 29.9792458; //speed of light in cm/ns
static const  float _c_convert_ = 10000000; //Convert c_special to m/s

//Beam energies in GeV
static const  float _energy_e16_ = 5.754;//5.7696;5.754
//3_14 was 5.754, but led to bad results: updated to 5.759 after reading Paremuzyan's thesis 
//Check 4.794 GeV for e16
static const  float _energy_e1f_ = 5.499;

static const float _beam_energy_[] = {_energy_e16_,_energy_e1f_};

static const  float _fine_structure_ = 1.0/137.0; 

//Masses of relevant particles
static const 	float _me_ = 0.0005109989; //mass of electron in GeV
static const 	float _mp_ = 0.93828;	//Mass of proton in GeV
static const 	float _mpi_ = 0.1395;	//Mass of pion in GeV

static const  float _degree_ = 180.0/TMath::Pi();


//Initial Four Vectors

const  TLorentzVector _k_mu_e16_(0.0,0.0,sqrt(_energy_e16_*_energy_e16_-_me_*_me_),_energy_e16_);
const  TLorentzVector _k_mu_e1f_(0.0,0.0,sqrt(_energy_e1f_*_energy_e1f_-_me_*_me_),_energy_e1f_);
const  TLorentzVector _p_mu_(0.0,0.0,0.0,_mp_);

//e1-6 Luminosity Values
static const  float lt_e16 = 5.0; //Target length in cm
static const  float Dt_e16 = 0.073; //Density of target in g/cm^3
static const  float NA = 6.022 * pow(10.0,23); //Avogadro's number
static const  float qe = 1.602 * pow(10.0,-19); // fundamental Coulomb charge 
static const  float Mt_e16 = 1.007; //Molar mass of target in g/mole
static const  float Qt_e16 = 21.32*pow(10.0,-3); //Total charge incident on target from Arjun in Coulombs //Will need to verify myself
static const  float L_e16 = Qt_e16*lt_e16*Dt_e16*NA/(qe*Mt_e16);
static const  float cm2_to_mbarn = pow(10.0,-27);

//Particle ID Numbers
//anti particles are negative
static const  	int _PROTON_ = 2212;
static const 	int _ELECTRON_ = 11;
static const 	int _PION_ = 211;
static const 	int _PION_0_ = 111;

// Beam spot from Arjun's Fortran vertex_e1f.f
static const  float _x_beam_ = 0.090;
static const  float _y_beam_ = -0.345;

//Fiducial Cut Parameters (electrons have e, hadrons have h)
//Arjun's cut_fid_e1f.f

static const  float pim_bot_MM[11] = {0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.0781216,0.083,0.109,0.109};
static const  float pim_top_MM[11] = {0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.2506164,0.34,0.415,0.415};


static const  float _XYmax_[3] = {545.0,545.0,545.0};
static const  int _XYres_ = 350;

//CC Projection Plane
static const  float _Acc_ = -0.000785;
static const  float _Bcc_ = 0.0;
static const  float _Ccc_ = -0.00168;
static const  float _Dcc_ = 1.0;


//binning
static const  float Wbin_res = 0.025;//The width of a W bin //30 steps
static const  float Wbin_start = 1.4;//The starting of W bins

static const  float Q2bin_res = 0.5;//6 steps
static const  float Q2bin_start = 2.0; 

static const  float pbin_res = 0.18;//Range: 0-5.0
static const  float pbin_start = 0.5;

//For cross sections
static const  Int_t xc_bins[7]= {29,5,14,14,10,10,10};


//Paths for file names
static const  std::string path1 = "/home/mclauchlinc/Desktop/analysis/nick_convert_e16.txt";
static const  std::string path2 = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/e16_10_18_17_ntple.txt";
static const  std::string path3 = "/Users/cmc/Desktop/analysis/Nick_skim_e16.txt";
static const  std::string path4 = "/Users/cmc/Desktop/analysis/arjun_sim.txt";
static const  std::string path3p = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/NickSkim_e16_PlateIN.txt";
static const  std::string path3n = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/NickSkim_e16_PlateOUT.txt";
static const  std::string paths1 = "/Users/cmc/Desktop/analysis/simulation/sim_e16_group1.txt";
static const  std::string paths2 = "/Users/cmc/Desktop/analysis/simulation/sim_e16_pre_gpp.txt";
static const  std::string paths3n = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gsim_no_gpp_nmcdata.txt";
static const  std::string paths3y = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gsim_no_gpp_mcdata.txt";
static const  std::string paths3 = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/sim_e16.txt";//"/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/sim_e16_fix.txt";
static const  std::string paths32 = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/sim_e16_2.txt";
static const  std::string paths_ar = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/arjun_sim_e16.txt";
static const  std::string path_ts = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/test_sim_e16.txt";
static const  std::string path_ts_gpp = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gpp_test_sim_e16.txt";
static const  std::string path_ts_gsim = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gsim_test_sim_e16.txt";
static const  std::string path_ts_gsim_1 = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/gsim_test_e16_1.txt";
static const  std::string paths_ar_a1c = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/arjun_recooked_gsim_a1c.txt";
static const  std::string paths_ar_ana = "/Users/cmc/Desktop/analysis/analysis_clas6/Path_Files/arjun_recooked_gsim_ana.txt";
static const  std::string path_ce16 = "/home/mclauchc/analysis/Path_Files/cluster_e16_path.txt";
static const  std::string path_ce16s = "/home/mclauchc/analysis/Path_Files/cluster_e16_sim_files.txt";
static const  std::string path_ce1f = "/home/mclauchc/analysis/Path_Files/cluster_e1f_file_paths.txt";
static const  std::string path_ce1fs = "/home/mclauchc/analysis/Path_Files/cluster_e1f_sim_files.txt";

static const  int _plate_sign_[] = {1,1};//1 is positive, -1 is negative
static const  int _plate_swap_e16_[] = {30703,30916,31143,31255};//Starts "out"
static const  int _plate_swap_e1f_[] = {38092,38114,38131,38132,38137,38143,38194,38199,38200,38203,38204,38207,38265,38290,38300,38548,38681};//"Starts "pos"
//static const  int _empty_e16_[] = {30825,30962,31104,31128,31252,31254,31300,31344};//From Arjun analysis
static const  int _empty_e16_[] = {30745,30746, 30824, 30825, 30867, 30961, 30962, 31029, 31104, 31252,31253,31254, 31298, 31299, 31300, 31344, 31382, 31393};//Beam and empty from e16_runs.csv
static const  int _other_e16_[] = {30909, 30910, 30767, 30775, 30776, 30932, 30933, 31076, 31343};//Empty, but bad B field, Moller, or other from e16_runs.csv
static const  int _beam_e16_[] = {30540,30541,30542,30549,30550,30551,30565,30568,30570,30571,30572,30579,30580,30582,30583,30584,30586,30587,30588,30589,30590,30591,30592,30593,30594,30595,30599,30600,30601,30602,30603,30604,30605,30606,30610,30612,30613,30614,30615,30616,30617,30618,30619,30620,30621,30622,30623,30624,30625,30626,30627,30628,30629,30631,30632,30633,30636,30638,30640,30678,30680,30681,30682,30684,30686,30687,30693,30694,30695,30698,30698,30699,30700,30701,30702,30703,30704,30705,30707,30708,30710,30711,30712,30713,30718,30719,30720,30722,30723,30724,30725,30726,30727,30728,30729,30730,30731,30732,30734,30735,30736,30737,30739,30740,30741,30742,30743,30745,30746,30747,30748,30749,30750,30753,30755,30756,30757,30758,30759,30760,30761,30762,30763,30766,30767,30768,30769,30773,30774,30775,30776,30777,30778,30779,30780,30781,30783,30784,30785,30786,30787,30789,30790,30791,30792,30794,30795,30799,30800,30801,30802,30803,30804,30805,30806,30807,30808,30809,30810,30811,30812,30813,30814,30815,30816,30817,30818,30819,30820,30821,30824,30825,30826,30827,30828,30829,30830,30831,30832,30833,30834,30835,30836,30837,30838,30839,30840,30841,30842,30843,30847,30848,30849,30850,30851,30852,30853,30854,30855,30856,30857,30858,30859,30860,30861,30862,30863,30864,30865,30866,30867,30909,30910,30912,30913,30914,30915,30916,30917,30918,30919,30920,30921,30922,30923,30924,30925,30926,30927,30928,30929,30930,30931,30932,30933,30934,30935,30940,30941,30942,30943,30944,30945,30946,30947,30948,30949,30950,30951,30952,30953,30954,30955,30956,30957,30958,30960,30961,30962,30963,30964,30965,30966,30967,30968,30969,30970,30971,30972,30973,30974,30975,30976,30977,30978,30985,30990,30991,30992,30993,30994,30995,30996,30997,30998,30999,31000,31001,31002,31003,31004,31005,31006,31007,31008,31009,31011,31012,31013,31014,31015,31016,31017,31018,31019,31029,31030,31031,31032,31033,31034,31035,31036,31037,31038,31039,31040,31041,31042,31043,31044,31045,31046,31047,31048,31049,31050,31051,31052,31053,31056,31057,31058,31059,31060,31061,31062,31063,31064,31065,31066,31068,31069,31070,31071,31072,31073,31074,31075,31076,31079,31080,31081,31082,31083,31084,31085,31086,31087,31088,31089,31090,31091,31092,31093,31094,31096,31097,31098,31099,31100,31101,31102,31103,31104,31108,31109,31110,31111,31112,31113,31114,31115,31116,31117,31118,31119,31120,31121,31122,31123,31124,31125,31141,31142,31143,31144,31145,31146,31147,31158,31159,31160,31162,31163,31164,31165,31166,31167,31168,31169,31170,31171,31172,31173,31174,31175,31176,31177,31178,31179,31180,31181,31182,31183,31184,31185,31186,31188,31189,31190,31191,31192,31193,31194,31195,31196,31197,31199,31200,31201,31203,31206,31207,31208,31209,31210,31211,31212,31213,31215,31216,31217,31218,31219,31220,31221,31222,31223,31224,31225,31226,31227,31228,31229,31231,31232,31233,31234,31235,31236,31237,31238,31239,31240,31241,31242,31243,31244,31245,31246,31247,31248,31249,31250,31252,31253,31254,31255,31256,31257,31258,31259,31260,31261,31262,31263,31264,31265,31266,31267,31268,31269,31270,31275,31276,31277,31278,31279,31280,31281,31282,31283,31284,31285,31286,31287,31288,31289,31290,31291,31292,31293,31294,31295,31296,31298,31299,31300,31301,31302,31303,31304,31305,31306,31307,31308,31309,31310,31311,31312,31313,31314,31315,31316,31317,31320,31321,31322,31323,31324,31325,31326,31327,31328,31329,31330,31331,31332,31333,31334,31335,31336,31337,31338,31339,31340,31341,31342,31343,31344,31345,31346,31347,31348,31349,31350,31351,31352,31353,31354,31355,31356,31357,31358,31359,31360,31361,31362,31363,31364,31365,31366,31367,31368,31369,31370,31371,31372,31373,31374,31375,31376,31377,31378,31379,31381,31382,31392,31393,31394,31395,31396,31397,31399,31400,31401,31402,31403,31404,31405,31406,31407,31408,31409,31410,31411,31412,31413,31414,31415,31416,31417,31418,31419,31422,31423,31424,31425,31426,31431,31432,31433,31434,31435,31436,31437,31438,31439,31440,31441,31442,31443,31444,31445,31446,31447,31448,31449,31450,31451,31452,31453,31454,31455,31456,31457,31458,31459,31460,31461,31462,31463,31464,31465,31466,31467,31468,31469,31470,31471,31472,31473,31474,31475,31476,31477,31478,31479,31480,31481,31482,31483,31484};//From e16_runs.csv
static const  int _empty_e1f_[] = {37854,37856,38111,38112,38587,38588,38590};
static const  int _other_e1f_[] = {10000,20000};
static const  int _beam_e1f_[] = {30000,40000};
static const  int _bad_helicity_e1f_[] = {38131};
static const  int _e16_run_bounds_[] = {30540,31484};
static const  int _e1f_run_bounds_[] = {37658,38751};

static const  int _e16_ = 0;
static const  int _e1f_ = 1; 

//Cut Names
static const  char * _none_ = "none";
static const  char * _sanity_ = "sanity";
static const  char * _fid_cut_ = "fid";
static const  char * _delta_cut_ = "delta";
static const  char * _sf_cut_ = "sf";
static const  char * _cc_cut_ = "min_cc";
static const  char * _ec_cut_ = "min_ec";
static const  char * _vertex_cut_ = "vertex";
static const  char * _beta_cut_ = "beta";
static const  char * _id_cut_ = "id";
static const  char * _pid_ = "pid";
static const  char * _event_ = "event";
static const  char * _ele_ = "ele";
static const  char * _pro_ = "pro";
static const  char * _pip_ = "pip";
static const  char * _pim_ = "pim";
//Event Indexes
static const  char * _mpro_ = "mpro";
static const  char * _mpip_ = "mpip";
static const  char * _mpim_ = "mpim";
static const  char * _mzero_ = "mzero";
static const  char * _mall_ = "mall";
static const  char * _mnone_ = "mnone";

static const  char * _nweighted_ = "no-weight";
static const  char * _weighted_ = "weight";

static const  char * _cut_applied_ = "cut";
static const  char * _anti_cut_ = "anti";
static const  char * _no_cut_ = "no-cut";

static const  char * _nthrown_ = "recon";
static const  char * _thrown_ = "thrown";

static const  char * _dirty_ = "dirty";
static const  char * _clean_ = "clean";
static const  char * _isolated_ = "isolated";

static const  char * _W_var_ = "W_Dep";
static const  char * _W_range_ = "in_range";
static const  char * _W_all_ = "all";
static const  char * _W_dep_[] = {_W_var_,_W_range_,_W_all_};

static const  char * _p_dep_ = "p_dep";
static const  char * _no_p_dep_ = "no_p_dep";

static const  char * _p_look_[] = {_no_p_dep_,_p_dep_};


static const  char * _true_ = "true";
static const  char * _false_ = "false";
static const  char * _truth_[] = {_false_,_true_};

//Parsing Out PID Cuts
static const  char* _ecuts_[] = {_none_,_sanity_, _fid_cut_, _sf_cut_, _cc_cut_, _ec_cut_, _vertex_cut_,_id_cut_,_pid_,_event_};//Be sure to add statements for flags if modified
static const  char* _hcuts_[] = {_none_,_sanity_, _fid_cut_, _delta_cut_,_id_cut_,_pid_,_event_};//Be sure to add statements for flags if modified
static const  char* _top_[] = {_mpro_,_mpip_,_mpim_,_mzero_,_mall_,_mnone_};
static const  char* _cut_[] = {_cut_applied_,_anti_cut_,_no_cut_};
static const  char* _species_[] = {_ele_,_pro_,_pip_,_pim_};
static const  char* _recon_[] = {_nthrown_,_thrown_};
static const  char* _weight_[] = {_nweighted_,_weighted_};
static const  char* _clean_event_[] = {_dirty_,_clean_,_isolated_};

static const  char* _sec1_ = "sec1";
static const  char* _sec2_ = "sec2";
static const  char* _sec3_ = "sec3";
static const  char* _sec4_ = "sec4";
static const  char* _sec5_ = "sec5";
static const  char* _sec6_ = "sec6";
static const  char* _sec_all_ = "all";
static const  char* _sector_[] = {_sec1_,_sec2_,_sec3_,_sec4_,_sec5_,_sec6_,_sec_all_};

static const  char* _left_ = "left";
static const  char* _right_ = "right";
static const  char* _coinc_ = "coinc";
static const  char* _cc_sides_[] = {_left_,_coinc_,_right_};

static const  char* _W_ = "W";
static const  char* _Q2_ = "Q2";
static const  char* _mpro_mpip_ = "MM_Delta++";
static const  char* _mpro_mpim_ = "MM_Delta0";
static const  char* _mpip_mpim_ = "MM_Rho";
static const  char* _MM1_[] = {_mpro_mpip_,_mpip_mpim_,_mpro_mpim_};
static const  char* _MM2_[] = {_mpip_mpim_,_mpro_mpim_,_mpro_mpip_};
static const  char* _theta_ = "Theta";
static const  char* _alpha_ = "Alpha";
static const  char* _phi_ = "Phi";
static const  char* _friend_pars_[] = {_W_,_Q2_,_MM1_[0],_MM2_[0],_theta_,_alpha_,_phi_};
static const  char* _var_names_[] =  {"#Delta^{++}","#rho","#Delta^{0}"};

static const  char* _ele_angle_corr_ = "e_theta_corr";
static const  char* _ele_p_corr_ = "e_p_corr";
static const  char* _no_corr_ = "no_corr";
static const  char* _ele_corr_[] = {_no_corr_,_ele_angle_corr_,_ele_p_corr_};



static const  char * species[] = {_ele_,_pro_,_pip_,_pim_};//4
static const  char * eid_cut[] = {"pre","sanity","fid","sf","min_cc","fid+sf","fid+cc","sf+cc","eid","bank","event"}; //11"min_cc","min_ec","eid","bank","event"};//"fid+sf","fid+cc","sf+cc","eid","bank","event"}; //11
static const  char * cut_ver[] = {"cut","anticut"};
static const  char * hid_cut[] = {"pre","sanity","fid","dt","hid","bank","event"};//,"pWQ2"}; //7
static const  char * topologies[] = {"None","Pmiss","PIPmiss","PIMmiss","Zeromiss","ALLmiss"}; //6
static const  char * sec_list[] = {"all_sectors","sec1","sec2","sec3","sec4","sec5","sec6"}; //7`
static const  char * W_dep_list[] = {"No_W_Dep","W_Dep"};
static const  char * CC_det_side[] = {"Left","Coince","Right","All"};
static const  char * basic_cut[] = {"pre","cut","anti"};
static const  char * MM_sq[] = {"linear","squared"};
static const  char * par_cut[4][11] = {{"pre","sanity","fid","sf","min_cc","fid+sf","fid+cc","sf+cc","eid","bank","event"},{"pre","sanity","fid","dt","hid","bank","event"},{"pre","sanity","fid","dt","hid","bank","event"},{"pre","sanity","fid","dt","hid","bank","event"}};
static const  char * fit_q[] = {"4fit","4show"};
static const  char * throw_stat[] = {"recon","thrown"};//For reconst ructed vs. thrown events
static const  char * w_stat[] = {"nweight","weight"};
static const  char * detectors[] = {"dc","cc","sc","ec"};
static const  char * inter_had[] = {"#Delta^{++}","#rho","#Delta^{0}"};

#endif
