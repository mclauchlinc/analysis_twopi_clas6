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
//std::shared_ptr<float> q_inc[15] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};//If this is less than the number of threads being attempted, it will fail
int _cut_par=0;
int _hist_par=0;
int _other_par=0;
std::string _file_list;

//double q_prev = NAN;
//std::shared_ptr<float> q_tot = NAN;


size_t run(std::shared_ptr<TChain> chain_, std::shared_ptr<Histogram> hists_, int thread_id_, std::shared_ptr<float> qtot_, std::shared_ptr<Flags> flags_){//, int &num_ppip){
	//Number of events in this thread
	std::cout<<"Time to run!\n\tCalculating number of events for this thread\n";
	int num_events = (long) chain_->GetEntries();
	//Print out information about the thread
	std::cout<<"Thread " <<thread_id_ <<": " <<num_events <<" Events\n";
	double q_prev=0.0;
	double q_tot=0.0;

	double q_prev_pos = 0.0;
	double q_tot_pos = 0.0;
	double q_prev_neg = 0.0;
	double q_tot_neg = 0.0;

	

	float curr_hel = 0.0;
	float prev_hel = 0.0;

	float delta_q2 = 0.0;
	float q_prev2 = 0.0;
	float hel_2 = 0.0;
	float hel_2_prev = 0.0;

	long number_of_events = 0;
	long number_of_recon = 0;
	//q_inc[thread_idx_] = 0.0;
	//if(std::distance(std::begin(q_inc), std::end(q_inc)) < (thread_id_+1)){
	//	std::cout<<"\n**ERROR TOO MANY THREADS FOR Q_INC**\n";
	//}
	//Make a data object which all the branches can be accessed from
	auto data = std::make_shared<Branches>(chain_,flags_->Flags::Sim());
	float pe = 0.0;
	float W = 0.0;
	float W0 = 0.0;
	float theta_e = 0.0;
	float theta_p = 0.0;
	float phi_e = 0.0;
	float phi_e_center =0.0;
	int sector = 0;

	int hel_idx = -1;
	float weight = 0.0;
	float Q2 = 0.0;


	float w0 = 0.0;
	float old_max_q = 0.0;
	float curr_max_q = 0.0;
	float run_sum_q = 0.0;
	int prev_run = 0;
	int run_evts = 0;

	std::string prev_filename = "";
	std::string curr_filename = "";


	int run_num = 0;//fun::extract_run_number(chain_->GetFile()->GetName(),flags_); //Not Finished so using temporary run number
	int prev_run_num = 0;
	float delta_q = 0.0;
	std::cout<<"Begin dive into Thread:" <<thread_id_ <<"\n";
	for(size_t curr_event = 0; curr_event < num_events; curr_event++){
		//Get singular event
		chain_->GetEntry(curr_event);
		if(!flags_->Sim()){
			run_num = fun::extract_run_number(chain_->GetFile()->GetName(),flags_);
			
			//std::cout<<"Run group: " <<flags_->Run() <<" Run Number is: " <<run_num <<" and it: ";
			if(flags_->Flags::Golden_Run()){
				if(flags_->Flags::Helicity()){
					hel_2 = fun::Corr_Helicity(data->Branches::evntclas2(),run_num,flags_);
				}else{
					hel_2 = 0.0;
				}
				if(curr_event == 0){
					prev_run = run_num;
					q_prev2 = data->Branches::q_l();
					delta_q2 = 0.0;
					hists_->Histogram::Golden_Indiv_Fill(run_num,delta_q2,hel_2);
				}else{
					curr_filename = chain_->GetFile()->GetName();
					if(prev_filename == curr_filename && data->Branches::q_l() >= q_prev2 && q_prev2>0.0){
					//if(prev_run == run_num && data->Branches::q_l() >= q_prev2){//&& curr_event != num_events-1 ){
						delta_q2 = data->Branches::q_l() - q_prev2; 
						hists_->Histogram::Golden_Indiv_Fill(run_num,delta_q2,hel_2_prev);
					}
					if(prev_filename == curr_filename && data->Branches::q_l() >= q_prev2){
						delta_q2 = data->Branches::q_l() - q_prev2; 
						hists_->Histogram::Golden_Indiv_Fill2(run_num,delta_q2,hel_2_prev);
					}
					if(prev_filename == curr_filename ){
						delta_q2 = data->Branches::q_l() - q_prev2; 
						hists_->Histogram::Golden_Indiv_Fill3(run_num,delta_q2,hel_2_prev);
					}
					if(prev_run == run_num  && data->Branches::q_l() >= q_prev2 && q_prev2>0.0){
					//if(prev_run == run_num && data->Branches::q_l() >= q_prev2){//&& curr_event != num_events-1 ){
						delta_q2 = data->Branches::q_l() - q_prev2; 
						hists_->Histogram::Golden_Indiv_Fill4(run_num,delta_q2,hel_2_prev);
					}
					if(prev_run == run_num  && data->Branches::q_l() >= q_prev2){
						delta_q2 = data->Branches::q_l() - q_prev2; 
						hists_->Histogram::Golden_Indiv_Fill5(run_num,delta_q2,hel_2_prev);
					}
					if(prev_run == run_num ){
						delta_q2 = data->Branches::q_l() - q_prev2; 
						hists_->Histogram::Golden_Indiv_Fill6(run_num,delta_q2,hel_2_prev);
					}
				}
				q_prev2 = data->Branches::q_l();
				hel_2_prev = hel_2;
			
				prev_run = run_num;
				prev_filename = chain_->GetFile()->GetName();
			}
			

		}
		/*
		if(flags_->Sim()){
			run_num = fun::extract_run_number_sim(chain_->GetFile()->GetName(),flags_);
			//std::cout<<"run_num sim: " <<run_num <<"\n";
			//if(flags_->Flags::Sim()){
			//	if(data->Branches::vz(0) < 0.0 && data->Branches::vz(0) > -8.0){
			//		hists_->Histogram::Sim_Vertex_Fill(run_num,flags_);
			//	}
			//}
		}*/
		

		
		if(fun::correct_run(run_num,flags_)){
			//std::cout<<" event is valid for what we are doing\n";
			//std::cout<<"Run Number is: " <<run_num <<" and it: passed";
			//Update on Progress through Analysis
			//if((thread_id_ == 0 || flags_->Flags::Make_Friend()) && curr_event%(num_events/100) == 0){
			if(thread_id_ == 0 && curr_event%(num_events/100) == 0){
			//if(thread_id_ == 0 && curr_event%(num_events/10000) == 0){
				//curr_file_name = 
				std::cout<<"\r" <<"\t" <<(100*curr_event/num_events) <<" %"  <<std::flush ;//<<"|| File: " <<chain_->GetFile()->GetName() <<std::flush;//;
				//std::cout<<"\r" <<"\t" <<(100.0*(float)(curr_event/num_events)) <<" %"  <<std::flush ;//<<"|| File: " <<chain_->GetFile()->GetName() <<std::flush;//;
				//std::cout<<"\r" <<"\t" <<(100*curr_event/num_events) <<"/10000"  <<std::flush ;//<<"|| File: " <<chain_->GetFile()->GetName() <<std::flush;//;
			}
			if(data->Branches::q_l()>0 && !flags_->Flags::Sim()){
				if(data->Branches::q_l()>q_prev && data->Branches::q_l()>0.0 && q_prev>0.0){
					if(run_num == prev_run_num || prev_run_num == 0){
						//delta_q = data->Branches::q_l() - q_prev;
						//qtot_ = std::make_shared<double>((double)data->Branches::q_l() - q_prev) + qtot_; 
						//*q_inc[thread_id_] += delta_q;
						q_tot += (data->Branches::q_l() - q_prev);
						//TThread::Lock();
						//*qtot_ += delta_q;
						//TThread::UnLock();
						//qtot_ += q_inc[thread_id_];
					}
				}
				q_prev = data->Branches::q_l();
			}
			//Particle ID, Event Selection, and Histogram Filling
			//std::cout<<"Going into a given event\n";
			number_of_events++;
			
			auto analysis = std::make_shared<Analysis>(data,hists_, thread_id_, run_num, flags_);
			if(flags_->Flags::Helicity() && !flags_->Flags::Sim()){
				curr_hel = analysis->Analysis::Corr_Helicity(data->Branches::evntclas2(),run_num,flags_);
				if(curr_hel == 1.0){
					if(data->Branches::q_l()>q_prev_pos && data->Branches::q_l()>0.0 && q_prev_pos>0.0){
						if(run_num == prev_run_num || prev_run_num == 0){
							if(prev_hel == -1.0){
								q_tot_neg += data->Branches::q_l() - q_prev_neg;
							}else if(prev_hel == 1.0){
								q_tot_pos += data->Branches::q_l() - q_prev_pos;
							}
						}
					}
					
				}else if(curr_hel == -1.0){
					if(data->Branches::q_l()>q_prev_neg && data->Branches::q_l()>0.0 && q_prev_neg>0.0){
						if(run_num == prev_run_num || prev_run_num == 0){
							if(prev_hel == 1.0){
								q_tot_pos += data->Branches::q_l() - q_prev_pos;
							}else if(prev_hel == -1.0){
								q_tot_neg += data->Branches::q_l() - q_prev_neg;
							}
						}
					}
				}
				q_prev_pos = data->Branches::q_l();
				q_prev_neg = data->Branches::q_l();
				prev_hel = curr_hel; 
			}
			
			
			if(analysis->Analysis::Valid_Event()){
				number_of_recon++;
			}
			if((flags_->Flags::Plot_Electron_Angle_Corr() || flags_->Flags::Plot_Elastic()) && !flags_->Flags::Sim() && data->Branches::gpart()>=2){
				phi_e = physics::get_phi(0,data);
				pe = data->Branches::p(0);
				sector = physics::get_sector(phi_e);
				theta_p = physics::get_theta(1,data);
				phi_e_center = physics::phi_center(phi_e);
				theta_e = physics::get_theta(0,data);
				W0 = physics::W(physics::Make_4Vector(true,pe,theta_e,phi_e,_me_),flags_->Flags::Run());
				hists_->Histogram::Elastic_Peak_Fill(W0,0,sector,flags_);
				if(W0 >= 0.7 && W0 < 1.05){
					hists_->Histogram::PCorr_Elast_Theta_Fill(theta_e, theta_p, sector, flags_);
				}
				//hists_->Histogram::Elastic_Peak_Fill(physics::W(physics::Make_4Vector(true,pe,physics::get_theta(0,data),phi_e,_me_),flags_->Flags::Run()),0,sector,flags_);
				if(flags_->Flags::E_Theta_Corr()){
					theta_e = corr::theta_e_corr(physics::get_theta(0,data),phi_e_center,flags_->Flags::Run(),true,sector-1);
					//std::cout<<"Correcting Electron Theta " <<physics::get_theta(0,data) <<" -> " <<theta_e <<"\n";
					W = physics::W(physics::Make_4Vector(true,pe,theta_e,phi_e,_me_),flags_->Flags::Run());
					hists_->Histogram::Elastic_Peak_Fill(W,1,sector,flags_);
				}else{
					theta_e = physics::get_theta(0,data);
					W = physics::W(physics::Make_4Vector(true,pe,theta_e,phi_e,_me_),flags_->Flags::Run());
				}
				hists_->Histogram::Ele_Angle_Corr_Fill(sector,theta_e,phi_e_center,W,theta_p,flags_);
			}
			//std::cout<<"Checking plot P mag corr\n";
			if(flags_->Flags::Plot_Electron_Mag_Corr() && !flags_->Flags::Sim() && data->Branches::gpart()>=2){
				//std::cout<<"Passed plot p mag corr \tChecking e theta corr corr\n";
				if(flags_->Flags::E_Theta_Corr()){
					//std::cout<<"passed e theta corr\n";
					//theta_e = physics::get_theta(0,data);
					phi_e = physics::get_phi(0,data);
					sector = physics::get_sector(phi_e);
					phi_e_center = physics::phi_center(phi_e);
					hists_->Histogram::Elastic_Peak_Fill(physics::W(physics::Make_4Vector(true,pe,physics::get_theta(0,data),phi_e,_me_),flags_->Flags::Run()),0,sector,flags_);
					theta_e = corr::theta_e_corr(physics::get_theta(0,data),phi_e_center,flags_->Flags::Run(),true,sector-1);
					if(flags_->Flags::E_PCorr()){
						pe = corr::p_corr_e(data->Branches::p(0),theta_e,phi_e_center,flags_->Flags::Run(),true,sector-1);
						W = physics::W(physics::Make_4Vector(true,pe,theta_e,phi_e,_me_),flags_->Flags::Run());
						hists_->Histogram::Elastic_Peak_Fill(W,2,sector,flags_);
					}else{
						pe = data->Branches::p(0);
						//std::cout<<"\npe:" <<pe <<" theta:" <<theta_e <<" phi:" <<phi_e ;
						//physics::Print_4Vec(physics::Make_4Vector(true,pe,theta_e,phi_e,_me_));
						W = physics::W(physics::Make_4Vector(true,pe,theta_e,phi_e,_me_),flags_->Flags::Run());
						//hists_->Histogram::Elastic_Peak_Fill(W,1,sector,flags_);
					}
					hists_->Histogram::Ele_Mag_Corr_Fill(pe,sector,theta_e, phi_e_center,W,flags_);
				}
			}
		}else{
			//if(run_num != prev_run_num){
			//	std::cout<<"\n\tfailed for run " <<run_num <<"\n";
			//}
		}
		prev_run_num = run_num;
	}
	if(!flags_->Flags::Sim()){
		std::cout<<"\n\n***For thread " <<thread_id_ <<" the integrated charge is: " <<q_tot <<"***\n\n";
		if(flags_->Flags::Helicity()){
			std::cout<<"\n\n***For thread " <<thread_id_ <<" the pos hel integrated charge is: " <<q_tot_pos <<"***\n\n";
			std::cout<<"\n\n***For thread " <<thread_id_ <<" the neg hel integrated charge is: " <<q_tot_neg <<"***\n\n";
		}
		//std::cout<<"\n\n***For thread " <<thread_id_ <<" the total integrated charge now is: " <<qtot_ <<"***\n\n";
		//std::cout<<"\n\n***For thread " <<thread_id_ <<" the integrated charge is: " <<*q_inc[thread_id_] <<"***\n\n";
		//It seems like the q_tot is adding up between all the threads, but separately too? So whatever the last thread to truly finish is ends up having the true q_tot
	}
	std::cout<<"\t**thread " <<thread_id_ <<" should have had " << num_events <<" events, but actually had" <<number_of_events <<" attempted events with " <<number_of_recon <<" events reconstructed\n";
}


size_t run_files( std::shared_ptr<Histogram> hists_, int thread_id_, int max_, std::shared_ptr<float> qtot_, std::shared_ptr<Flags> flags_){//, int &num_ppip){std::vector<std::string> files_,
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