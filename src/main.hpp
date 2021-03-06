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




size_t run(std::shared_ptr<TChain> chain_, std::shared_ptr<Histogram> hists_, int thread_id_, std::shared_ptr<Flags> flags_){//, int &num_ppip){
	//Number of events in this thread
	int num_events = (int) chain_->GetEntries();
	//Print out information about the thread
	std::cout<<"Thread " <<thread_id_ <<": " <<num_events <<" Events\n";
	
	//Make a data object which all the branches can be accessed from
	auto data = std::make_shared<Branches>(chain_,flags_->Flags::Sim());

	int run_num = 0;fun::extract_run_number(chain_->GetFile()->GetName(),flags_); //Not Finished so using temporary run number
	for(size_t curr_event = 0; curr_event < num_events; curr_event++){
		//Get singular event
		chain_->GetEntry(curr_event);
		run_num = fun::extract_run_number(chain_->GetFile()->GetName(),flags_);
		//std::cout<<"Run Number is: " <<run_num <<" and it: ";
		if(fun::correct_run(run_num,flags_)){
			//std::cout<<"passed\n";
			//Update on Progress through Analysis
			if(thread_id_ == 0 && curr_event%(num_events/100) == 0){
				//curr_file_name = 
				std::cout<<"\r" <<"\t" <<(100*curr_event/num_events) <<" %"  <<std::flush ;//<<"|| File: " <<chain_->GetFile()->GetName() <<std::flush;//;
			}
			//Particle ID, Event Selection, and Histogram Filling
			auto analysis = std::make_shared<Analysis>(data,hists_, thread_id_, run_num, flags_);
		}else{
			//std::cout<<"failed\n";
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