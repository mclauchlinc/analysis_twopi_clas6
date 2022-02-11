#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include "TFile.h"
#include "TH1.h"
//#include "event_class.hpp"
#include "histogram.hpp" 
#include "constants.hpp"
#include "functions.hpp"
#include "branches.hpp"
#include <unistd.h>
#include "corrections.hpp"
#include "cuts.hpp"

std::string comp; //Variable for choosing which data set will be used
char* output_name;//Variable for the output file name. This is reassigned through input parameters
int file_num = -1;//The initial assignment for the number of files in the program. -1 equates to all of them
char envi_name[100];//name for the environment file
char output_name2[100];
char* curr_file_name;
bool cluster = false;
std::vector<float> fara_q_charge;
std::vector<int> fara_q_run;
std::vector<float> fara_q_size;
std::vector<std::vector<float>> fara_q_run_f;
std::vector<std::vector<float>> fara_q_charge_f;
std::vector<std::vector<float>> fara_q_size_f;
int run_n = 0;
//int old_run_num[_NUM_THREADS_];
//int run_number[_NUM_THREADS_];

int _cut_par=0;
int _hist_par=0;
int _other_par=0;
std::string _file_list;

float _theta_cut_ = 35.0;


size_t run(std::shared_ptr<TChain> chain_, std::shared_ptr<Histogram> hists_, int thread_id_, std::shared_ptr<Flags> flags_){//, int &num_ppip){
	//Number of events in this thread
	int num_events = (int) chain_->GetEntries();
	//Print out information about the thread
	std::cout<<"Thread " <<thread_id_ <<": " <<num_events <<" Events\n";
	
	//Make a data object which all the branches can be accessed from
	auto data = std::make_shared<Branches>(chain_,flags_->Flags::Sim());

	//int run_num = fun::extract_run_number(chain_->GetFile()->GetName(),flags_); //Not Finished so using temporary run number
	for(size_t curr_event = 0; curr_event < num_events; curr_event++){
		//Get singular event
		chain_->GetEntry(curr_event);
		//Update on Progress through Analysis
		if(thread_id_ == 0 && curr_event%(num_events/100) == 0){
			//curr_file_name = 
			std::cout<<"\r" <<"\t" <<(100*curr_event/num_events) <<" %"  <<std::flush ;//<<"|| File: " <<chain_->GetFile()->GetName() <<std::flush;//;
		}
		if(data->gpart()>=2 && data->q(1)>0.0){
			//std::cout<<"We are indeed plotting pcorr. Let's make some W's\n";
			float p_e = data->p(0);

			float theta_p = fun::theta(data,1);
			float theta_e = fun::theta(data,0);
			float phi_e = fun::phi(data,0);
			float phi_e_center = fun::phi(data,0,true);
			//std::cout<<"Electron Phi: " <<phi_e <<" | centered: " <<phi_e_center <<" gives bin: " <<hists_->Histogram::Phi_Idx(phi_e_center) <<" which says range is: " <<hists_->Histogram::Phi_Low(hists_->Histogram::Phi_Idx(phi_e_center)) <<"_" <<hists_->Histogram::Phi_Top(hists_->Histogram::Phi_Idx(phi_e_center)) <<"\n";
			std::cout<<"Electron Theta: " <<theta_e <<" gives bin: " <<hists_->Histogram::Theta_Idx(theta_e) <<" which says range is: " <<hists_->Histogram::Theta_Low(hists_->Histogram::Theta_Idx(theta_e)) <<"_" <<hists_->Histogram::Theta_Top(hists_->Histogram::Theta_Idx(theta_e)) <<"\n";
			int sector_e = fun::get_sector(data, 0);
			float theta_e_corr = corr::theta_e_corr(theta_e,phi_e_center,flags_->Run(),true,sector_e);
			TLorentzVector k_mu = fun::Make_4Vector(p_e, theta_e, phi_e, _me_);
			float W = fun::W(k_mu,flags_->Run());
			TLorentzVector k_mu_acorr = fun::Make_4Vector(p_e, theta_e_corr, phi_e, _me_);
			float W_acorr = fun::W(k_mu_acorr,flags_->Run());
			//TLorentzVector k_mu_pcorr = fun::Make_4Vector(corr::p_corr_e(data->p(0),corr::theta_e_corr(fun::theta(data,0),fun::phi(data,0,true),flags_->Run(),true,),fun::phi(data,0,true),flags_->Run()), corr::theta_e_corr(fun::theta(data,0),fun::phi(data,0,true),flags_->Run()), fun::phi(data,0), _me_);
			//float W_pcorr = fun::W(k_mu_pcorr,flags_->Run());
			if(flags_->Plot_Check()){
				if(cuts::W_cut(W,flags_->Run())){
					hists_->Histogram::Check_Fill(fun::theta(data,1),fun::theta(data,0),fun::get_sector(data, 0),_check_1_,flags_);
					hists_->Histogram::Check_Fill(fun::theta_calc(fun::theta(data,1), _beam_energy_[flags_->Run()]),fun::theta(data,0),fun::get_sector(data, 0),_check_2_,flags_);
					hists_->Histogram::Check_Fill(fun::delta_theta(fun::theta(data,0), fun::theta(data,1), _beam_energy_[flags_->Run()]),fun::theta(data,0),fun::get_sector(data, 0),_check_3_,flags_);
					hists_->Histogram::Check_Fill(fun::delta_theta(fun::theta(data,0), fun::theta(data,1), _beam_energy_[flags_->Run()]),fun::phi(data,0,true),fun::get_sector(data, 0),_check_4_,flags_);
				}
			}
			//Elastic Histograms
			if(flags_->Plot_Elastic()){
				//std::cout<<"Filling Elastic Histograms\n";
				if(flags_->Flags::Fid_Cut(0)){
					hists_->Histogram::Angular_Fill(theta_e, phi_e_center, sector_e,_no_corr_, _no_pro_thresh_, _no_cut_,flags_);
					if(theta_p>=_theta_cut_){
						hists_->Histogram::Angular_Fill(theta_e, phi_e_center, fun::get_sector(data, 0), _no_corr_, _pro_thresh_, _no_cut_,flags_);
					}
					if(cuts::fid_e(data->p(0),theta_e,phi_e_center,flags_)){
						hists_->Histogram::Elastic_Fill(W,sector_e,_no_corr_,_no_pro_thresh_,flags_);
						hists_->Histogram::Angular_Fill(theta_e, phi_e_center, sector_e,_no_corr_, _no_pro_thresh_, _fid_cut_,flags_);
						if(theta_p>=_theta_cut_){
							hists_->Histogram::Elastic_Fill(W,sector_e,_no_corr_,_pro_thresh_,flags_);
							hists_->Histogram::Angular_Fill(theta_e, phi_e_center, sector_e, _no_corr_, _pro_thresh_, _fid_cut_,flags_);
						}
					}
					if(flags_->Flags::E_Theta_Corr() && cuts::fid_e(data->p(0),theta_e_corr,phi_e_center,flags_)){
						//TLorentzVector k_mu_acorr = fun::Make_4Vector(data->p(0), theta_e_corr, fun::phi(data,0), _me_);
						//float W_acorr = fun::W(k_mu_acorr,flags_->Run());
						hists_->Histogram::Elastic_Fill(W_acorr,sector_e,_ele_angle_corr_,_no_pro_thresh_,flags_);
						hists_->Histogram::Angular_Fill(theta_e_corr, phi_e_center, sector_e,_ele_angle_corr_, _no_pro_thresh_, _no_cut_,flags_);
						if(theta_p>=_theta_cut_){
							hists_->Histogram::Elastic_Fill(W_acorr,sector_e,_ele_angle_corr_,_pro_thresh_,flags_);
							hists_->Histogram::Angular_Fill(theta_e_corr, phi_e_center, sector_e,_ele_angle_corr_, _pro_thresh_, _no_cut_,flags_);
						}
					}
					if(flags_->Flags::E_PCorr() && cuts::fid_e(corr::p_corr_e(data->p(0),theta_e_corr,phi_e_center,flags_->Run()),theta_e_corr,phi_e_center,flags_)){
						//TLorentzVector k_mu_pcorr = fun::Make_4Vector(corr::p_corr_e(data->p(0), theta_e_corr, fun::phi(data,0), _me_);
						//float W_pcorr = fun::W(k_mu_pcorr,flags_->Run());
						//hists_->Histogram::Elastic_Fill(W_pcorr,sector_e,_ele_p_corr_,_no_pro_thresh_,flags_);
						//hists_->Histogram::Angular_Fill(theta_e_corr, phi_e_center, sector_e,_ele_p_corr_, _no_pro_thresh_, _no_cut_,flags_);
						//if(theta_p>=_theta_cut_){
						//	hists_->Histogram::Elastic_Fill(W_pcorr,sector_e,_ele_p_corr_,_pro_thresh_,flags_);
						//	hists_->Histogram::Angular_Fill(theta_e_corr, phi_e_center, sector_e,_ele_p_corr_, _pro_thresh_, _no_cut_,flags_);
						//}
					}
				}else{
					hists_->Histogram::Elastic_Fill(W,sector_e,_no_corr_,_no_pro_thresh_,flags_);
					hists_->Histogram::Angular_Fill(fun::theta(data,0), phi_e_center, sector_e,_no_corr_, _no_pro_thresh_, _no_cut_,flags_);
					if(fun::theta(data,1)>=_theta_cut_){
						hists_->Histogram::Elastic_Fill(W,sector_e,_no_corr_,_pro_thresh_,flags_);
						hists_->Histogram::Angular_Fill(fun::theta(data,0), phi_e_center, sector_e, _no_corr_, _pro_thresh_, _no_cut_,flags_);
					}
					if(flags_->Flags::E_Theta_Corr()){
						//TLorentzVector k_mu_acorr = fun::Make_4Vector(data->p(0), theta_e_corr, fun::phi(data,0), _me_);
						//float W_acorr = fun::W(k_mu_acorr,flags_->Run());
						hists_->Histogram::Elastic_Fill(W_acorr,sector_e,_ele_angle_corr_,_no_pro_thresh_,flags_);
						hists_->Histogram::Angular_Fill(theta_e_corr, phi_e_center, sector_e,_ele_angle_corr_, _no_pro_thresh_, _no_cut_,flags_);
						if(theta_p>=_theta_cut_){
							hists_->Histogram::Elastic_Fill(W_acorr,sector_e,_ele_angle_corr_,_pro_thresh_,flags_);
							hists_->Histogram::Angular_Fill(theta_e_corr, phi_e_center, sector_e,_ele_angle_corr_, _pro_thresh_, _no_cut_,flags_);
						}
					}
					if(flags_->Flags::E_PCorr() && cuts::fid_e(corr::p_corr_e(data->p(0),theta_e_corr,phi_e_center,flags_->Run()),theta_e_corr,phi_e_center,flags_)){
						TLorentzVector k_mu_pcorr = fun::Make_4Vector(corr::p_corr_e(data->p(0),theta_e_corr,phi_e_center,flags_->Run()), theta_e_corr, fun::phi(data,0), _me_);
						float W_pcorr = fun::W(k_mu_pcorr,flags_->Run());
						hists_->Histogram::Elastic_Fill(W_pcorr,sector_e,_ele_p_corr_,_no_pro_thresh_,flags_);
						hists_->Histogram::Angular_Fill(theta_e_corr, phi_e_center, sector_e,_ele_p_corr_, _no_pro_thresh_, _no_cut_,flags_);
						if(theta_p>=_theta_cut_){
							hists_->Histogram::Elastic_Fill(W_pcorr,sector_e,_ele_p_corr_,_pro_thresh_,flags_);
							hists_->Histogram::Angular_Fill(theta_e_corr, phi_e_center, sector_e,_ele_p_corr_, _pro_thresh_, _no_cut_,flags_);
						}
					}
				}
			}
			//Delta Theta Histograms
			if(flags_->Plot_E_PCorr()){
				if(flags_->Flags::Fid_Cut(0)){
					if(cuts::fid_e(p_e,theta_e,phi_e_center,flags_)){
					//std::cout<<"Filling Delta Theta Histograms\n";
						if(cuts::W_cut(W,flags_->Run())){
							if(theta_p>=_theta_cut_){	
								float delta_theta = fun::delta_theta(theta_e, theta_p, _beam_energy_[flags_->Run()]);
								hists_->Histogram::ECorr_Angle_Fill(delta_theta, theta_e,phi_e_center,fun::get_sector(data, 0), _no_corr_, flags_);
							}
						}
					}
					if(flags_->Flags::E_Theta_Corr() && cuts::fid_e(p_e,theta_e_corr,phi_e_center,flags_)){	
						if(cuts::W_cut(W_acorr,flags_->Run())){
							if(theta_p>=_theta_cut_){
								float delta_theta_acorr = fun::delta_theta(theta_e_corr, theta_p, _beam_energy_[flags_->Run()]);
								hists_->Histogram::ECorr_Angle_Fill(delta_theta_acorr, theta_e_corr,phi_e_center,fun::get_sector(data, 0), _ele_angle_corr_, flags_);
							}
						}
					}
					if(flags_->Flags::E_PCorr()){
						TLorentzVector k_mu_pcorr = fun::Make_4Vector(corr::p_corr_e(p_e,theta_e_corr,phi_e_center,flags_->Run()), theta_e_corr, fun::phi(data,0), _me_);
						float W_pcorr = fun::W(k_mu_pcorr,flags_->Run());
						if(cuts::W_cut(W_pcorr,flags_->Run())){
							if(theta_p>=_theta_cut_){
								
								float delta_theta_pcorr = fun::delta_theta(theta_e_corr, theta_p, _beam_energy_[flags_->Run()]);
								hists_->Histogram::ECorr_Angle_Fill(delta_theta_pcorr, theta_e,phi_e_center,sector_e, _ele_p_corr_, flags_);
							}
						}
					}
				}else{
				//std::cout<<"Filling Delta Theta Histograms\n";
					if(cuts::W_cut(W,flags_->Run())){
						if(theta_p>=_theta_cut_){	
							float delta_theta = fun::delta_theta(theta_e, theta_p, _beam_energy_[flags_->Run()]);
							hists_->Histogram::ECorr_Angle_Fill(delta_theta, theta_e,phi_e_center,sector_e, _no_corr_, flags_);
						}
					}
					if(flags_->Flags::E_Theta_Corr()){	
						//TLorentzVector k_mu_acorr = fun::Make_4Vector(p_e, theta_e_corr, fun::phi(data,0), _me_);
						//float W_acorr = fun::W(k_mu_acorr,flags_->Run());
						if(cuts::W_cut(W_acorr,flags_->Run())){
							if(theta_e_corr>=_theta_cut_){
								float delta_theta_acorr = fun::delta_theta(theta_e_corr, theta_p, _beam_energy_[flags_->Run()]);
								hists_->Histogram::ECorr_Angle_Fill(delta_theta_acorr, theta_e_corr,phi_e_center,sector_e, _ele_angle_corr_, flags_);
							}
						}
					}
					if(flags_->Flags::E_PCorr()){
						TLorentzVector k_mu_pcorr = fun::Make_4Vector(corr::p_corr_e(data->p(0),theta_e_corr,phi_e_center,flags_->Run()), theta_e_corr,phi_e, _me_);
						float W_pcorr = fun::W(k_mu_pcorr,flags_->Run());
						if(cuts::W_cut(W_pcorr,flags_->Run())){
							if(theta_e_corr>=_theta_cut_){
								
								float delta_theta_pcorr = fun::delta_theta(theta_e_corr, theta_p, _beam_energy_[flags_->Run()]);
								hists_->Histogram::ECorr_Angle_Fill(delta_theta_pcorr, theta_e_corr,phi_e_center,sector_e, _ele_p_corr_, flags_);
							}
						}
					}
				}
			}
			//if(flags_->Flags::W_Plot()){
			//	hists_->Histogram::W_Fill(W,fun::get_sector(data,0));
			//}
		}
	}
}


size_t run_files( std::shared_ptr<Histogram> hists_, int thread_id_, int max_, std::shared_ptr<Flags> flags_){//, int &num_ppip){std::vector<std::string> files_,
	//Called once per thread
	//Make a new chain to process for this thread
	auto chain = std::make_shared<TChain>("h10");
	//Add every file to the chain
	std::cout<<"Loading Chain\n";
	fun::loadChain(chain, flags_->Flags::Files(), thread_id_, flags_->Flags::Num_Files(),flags_);
	//for(auto in:files_) chain->Add(in.c_str());
	//Run the function over each thread
	return run(chain,hists_,thread_id_,flags_);//,num_ppip);
}

#endif