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
		if(data->gpart()==2 && flags_->Plot_E_PCorr()){
			TLorentzVector k_mu = fun::Make_4Vector(data->p(0), fun::theta(data,0), fun::phi(data,0), _me_);
			float W = fun::W(k_mu,flags_->Run());
			if(cuts::W_cut(W,flags_->Run())){
				if(fun::theta(data,1)>=_theta_cut_){	
					float delta_theta = fun::delta_theta(fun::theta(data,0), fun::theta(data,1), _beam_energy_[flags_->Run()]);
					hists_->Histogram::ECorr_Angle_Fill(delta_theta, fun::theta(data,0),fun::phi(data,0,true),fun::get_sector(data, 0), _no_corr_, flags_);
				}
				if(flags_->Flags::E_Theta_Corr()){
					TLorentzVector k_mu_acorr = fun::Make_4Vector(data->p(0), corr::theta_e_corr(fun::theta(data,0),fun::phi(data,0,true),flags_->Run()), fun::phi(data,0), _me_);
					float W_acorr = fun::W(k_mu_acorr,flags_->Run());
					if(fun::theta(data,1)>=_theta_cut_){
						float delta_theta_acorr = fun::delta_theta(corr::theta_e_corr(fun::theta(data,0),fun::phi(data,0,true),flags_->Run()), fun::theta(data,1), _beam_energy_[flags_->Run()]);
						hists_->Histogram::ECorr_Angle_Fill(delta_theta_acorr, corr::theta_e_corr(fun::theta(data,0),fun::phi(data,0,true),flags_->Run()),fun::phi(data,0,true),fun::get_sector(data, 0), _ele_angle_corr_, flags_);
					}
				}
				if(flags_->Flags::E_PCorr()){
					TLorentzVector k_mu_pcorr = fun::Make_4Vector(corr::p_corr_e(data->p(0),corr::theta_e_corr(fun::theta(data,0),fun::phi(data,0,true),flags_->Run()),fun::phi(data,0,true),flags_->Run()), corr::theta_e_corr(fun::theta(data,0),fun::phi(data,0,true),flags_->Run()), fun::phi(data,0), _me_);
					float W_acorr = fun::W(k_mu_pcorr,flags_->Run());
					if(fun::theta(data,1)>=_theta_cut_){
						float delta_theta_pcorr = fun::delta_theta(corr::theta_e_corr(fun::theta(data,0),fun::phi(data,0,true),flags_->Run()), fun::theta(data,1), _beam_energy_[flags_->Run()]);
						hists_->Histogram::ECorr_Angle_Fill(delta_theta_pcorr, fun::theta(data,0),fun::phi(data,0,true),fun::get_sector(data, 0), _ele_angle_corr_, flags_);
					}
				}
			}
			//if(flags_->Flags::W_Plot()){
			//	hists_->Histogram::W_Fill(W,fun::get_sector(data,0));
			//}
		}

		//Particle ID, Event Selection, and Histogram Filling
		//auto analysis = std::make_shared<Analysis>(data,hists_, thread_id_, run_num, flags_);
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