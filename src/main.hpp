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
#include "forest.hpp"
#include "event_analysis.hpp"
#include <unistd.h>

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

//double q_prev = NAN;
//double q_tot = NAN;


size_t run(std::shared_ptr<TChain> chain_, std::shared_ptr<Histogram> hists_, int thread_id_, std::shared_ptr<double> qtot_, std::shared_ptr<Flags> flags_){//, int &num_ppip){
	//Number of events in this thread
	std::cout<<"Time to run!\n\tCalculating number of events for this thread\n";
	int num_events = (long) chain_->GetEntries();
	//Print out information about the thread
	std::cout<<"Thread " <<thread_id_ <<": " <<num_events <<" Events\n";
	double q_prev=0.0;
	double q_tot=0.0;
	//Make a data object which all the branches can be accessed from
	auto data = std::make_shared<Branches>(chain_,flags_->Flags::Sim());
	float pe = 0.0;
	float W = 0.0;
	float theta_e = 0.0;
	float theta_p = 0.0;
	float phi_e = 0.0;
	float phi_e_center =0.0;
	int sector = 0;

	int run_num = 0;//fun::extract_run_number(chain_->GetFile()->GetName(),flags_); //Not Finished so using temporary run number
	for(size_t curr_event = 0; curr_event < num_events; curr_event++){
		//Get singular event
		chain_->GetEntry(curr_event);
		run_num = fun::extract_run_number(chain_->GetFile()->GetName(),flags_);
		//std::cout<<"Run Number is: " <<run_num <<" and it: ";
		if(fun::correct_run(run_num,flags_)){
			//std::cout<<"passed\n";
			//std::cout<<"Run Number is: " <<run_num <<" and it: passed";
			//Update on Progress through Analysis
			//if((thread_id_ == 0 || flags_->Flags::Make_Friend()) && curr_event%(num_events/100) == 0){
			if(thread_id_ == 0 && curr_event%(num_events/100) == 0){
				//curr_file_name = 
				std::cout<<"\r" <<"\t" <<(100*curr_event/num_events) <<" %"  <<std::flush ;//<<"|| File: " <<chain_->GetFile()->GetName() <<std::flush;//;
			}
			if(data->Branches::q_l()>0 && !flags_->Flags::Sim()){
				if(data->Branches::q_l()>q_prev && data->Branches::q_l()>0){
					//qtot_ = std::make_shared<double>((double)data->Branches::q_l() - q_prev) + qtot_; 
					q_tot += data->Branches::q_l() - q_prev;
				}
				q_prev = data->Branches::q_l();
			}
			//Particle ID, Event Selection, and Histogram Filling
			auto analysis = std::make_shared<Analysis>(data,hists_, thread_id_, run_num, flags_);
			if(flags_->Flags::Plot_Electron_Angle_Corr() && !flags_->Flags::Sim() && data->Branches::gpart()>=2){
				phi_e = physics::get_phi(0,data);
				pe = data->Branches::p(0);
				sector = physics::get_sector(phi_e);
				theta_p = physics::get_theta(1,data);
				phi_e_center = physics::phi_center(phi_e);
				if(flags_->Flags::E_Theta_Corr()){
					theta_e = corr::theta_e_corr(physics::get_theta(0,data),phi_e_center,flags_->Flags::Run(),true,sector-1);
					//std::cout<<"Correcting Electron Theta " <<physics::get_theta(0,data) <<" -> " <<theta_e <<"\n";
					W = physics::W(physics::Make_4Vector(true,pe,theta_e,phi_e_center,_me_),flags_->Flags::Run());
				}else{
					theta_e = physics::get_theta(0,data);
					W = physics::W(physics::Make_4Vector(true,pe,theta_e,phi_e,_me_),flags_->Flags::Run());
				}
				hists_->Histogram::Ele_Angle_Corr_Fill(sector,theta_e,phi_e_center,W,theta_p,flags_);
			}
			if(flags_->Flags::Plot_Electron_Mag_Corr() && !flags_->Flags::Sim() && data->Branches::gpart()>=2){
				if(flags_->Flags::E_Theta_Corr()){
					theta_e = corr::theta_e_corr(physics::get_theta(0,data),phi_e,flags_->Flags::Run(),false,sector-1);
					phi_e = physics::get_phi(0,data);
					sector = physics::get_sector(phi_e);
					phi_e_center = physics::phi_center(phi_e);
					if(flags_->Flags::E_PCorr()){
						pe = corr::p_corr_e(data->Branches::p(0),theta_e,phi_e_center,flags_->Flags::Run(),true,sector-1);
						W = physics::W(physics::Make_4Vector(true,pe,theta_e,phi_e_center,_me_),flags_->Flags::Run());
					}else{
						pe = data->Branches::p(0);
						W = physics::W(physics::Make_4Vector(true,pe,theta_e,phi_e,_me_),flags_->Flags::Run());
					}
					hists_->Histogram::Ele_Mag_Corr_Fill(pe,sector,theta_e, phi_e_center,W,flags_);
				}
			}
		}else{
			//std::cout<<"failed\n";
		}
	}
	if(!flags_->Flags::Sim()){
		std::cout<<"\n\n***For thread " <<thread_id_ <<" the integrated charge is: " <<q_tot <<"***\n\n";
		//It seems like the q_tot is adding up between all the threads, but separately too? So whatever the last thread to truly finish is ends up having the true q_tot
	}
}


size_t run_files( std::shared_ptr<Histogram> hists_, int thread_id_, int max_, std::shared_ptr<double> qtot_, std::shared_ptr<Flags> flags_){//, int &num_ppip){std::vector<std::string> files_,
	//Called once per thread
	//Make a new chain to process for this thread
	auto chain = std::make_shared<TChain>("h10");
	//Add every file to the chain
	std::cout<<"Loading Chain\n";
	fun::loadChain(chain, flags_->Flags::Files(), thread_id_, flags_->Flags::Num_Files(),flags_);
	//for(auto in:files_) chain->Add(in.c_str());
	//Run the function over each thread
	return run(chain,hists_,thread_id_,qtot_,flags_);//,num_ppip);
}

#endif