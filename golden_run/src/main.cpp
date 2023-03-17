#include "main.hpp"
#include <future>
#include <thread>
#include "TROOT.h"

//./golden_run $dataset_$where #files
//dataset -> {exp_e16,exp_e1f,empty_e16,empty_e1f}
//where->{local_cluster}

//Main body
int main(int argc, char **argv){
	std::cout<<"\nBegin\n";
	//start timer
	auto start = std::chrono::high_resolution_clock::now();
	if(argc!=9){
		std::cout<<"Did not enter correct number of arguments\n";
		return 0;
	}
	std::cout<<"Reading in Parameters\n";
	std::string file_location = argv[1];
	int file_num = std::atoi(argv[2]);
	std::string output_name = argv[3];
	std::string remove_front = argv[4];
	std::string remove_mid = argv[5];
	float low_bound = atof(argv[6]);
	float top_bound = atof(argv[7]);
	std::string run_group = argv[8];

	auto chain = std::make_shared<TChain>("h10");
	//Add every file to the chain
	std::cout<<"Loading Data from " <<file_location <<"\n";
	fun::loadChain(chain,file_location,file_num);
	
	//Make output directory
	std::string out_dir;

	
	size_t num_of_events = (int) chain->GetEntries();
	std::cout<<"Loaded " <<num_of_events <<" events\n";
	//Get the Data
	auto data = std::make_shared<Branches>(chain);

	//Make Histograms
	std::cout<<"Making Histograms\n";
	auto hist = Histogram(run_group);

	std::vector<int> run_nums;
	float old_max_q;//Faraday Cup charge from previous event
	float curr_max_q;//Faraday Cup charge from current event
	std::string curr_file_name;//File name for current event
	std::string old_file_name;//File name for previous event
	float sum_q;//Summed Faraday Cup Charge
	int run_seg;//Segment within Run
	int run_num;//Run number
	float sum_q_segment = 0.0;//Faraday Summed Charge on segment
	float sum_q_run = 0.0;//Faraday Summed Charge on run
	int run_seg_size = 0;//Number events on segment
	int run_num_size = 0;//Number events on Run

	float max_q_seg = 0.0;
	float max_event_seg = 0.0;
	float max_q_run = 0.0;
	float max_event_run = 0.0;
	float min_q_seg = NAN;
	float min_event_seg = NAN;
	float min_q_run = NAN;
	float min_event_run = NAN;


	std::cout<<"Begin Charge Extraction\n";
	chain->GetEntry(0);
	old_max_q = data->Branches::q_l();
	old_file_name = chain->GetFile()->GetName();
	//std::cout<<"File: " <<old_file_name <<"\n";
	for(size_t curr_event = 0; curr_event < num_of_events; curr_event++){
			//Get singular event
			chain->GetEntry(curr_event);
			curr_file_name = chain->GetFile()->GetName();
			if(curr_file_name == old_file_name){//If we're in the same file
				if(old_max_q==0.0){//First q_l of every file is zero, and will have an artificial jump otherwise
					old_max_q = data->Branches::q_l();
				}
				curr_max_q = data->Branches::q_l();
				if(curr_max_q > old_max_q){
					sum_q += curr_max_q - old_max_q;
					sum_q_segment += curr_max_q - old_max_q;
					sum_q_run += curr_max_q - old_max_q;
					//std::cout<<"\t" <<"run:" <<run_num <<" seg:" <<run_seg <<" sum_q added "<<curr_max_q-old_max_q <<" to be: " <<sum_q <<" seg:" <<sum_q_segment <<" run:" <<sum_q_run <<"\n";
					//std::cout<<"\t\tOld " <<old_max_q <<" and New: " <<curr_max_q <<"\n";
				}
				run_num = fun::run_number(chain->GetFile()->GetName(),remove_front);
				if(fun::run_num_idx(run_num,run_nums)==-1){
					run_nums.push_back(run_num);
				}
				run_seg = fun::run_segment(chain->GetFile()->GetName(),remove_front,remove_mid);
				run_seg_size+=1;
				run_num_size+=1;
				if(curr_event == num_of_events){
					hist.Fill_Norm_Seg(run_num,run_seg,sum_q_segment/(float)run_seg_size,run_nums);
					hist.Fill_Norm_Run(run_num,sum_q_run/(float)run_num_size);
					hist.Fill_Charge_v_Event_Seg(sum_q_segment,(float)run_seg_size);
					hist.Fill_Charge_v_Event_Run(sum_q_run,(float)run_num_size);
				}
				old_max_q = data->Branches::q_l();//Prepping charge for next event
				old_file_name = chain->GetFile()->GetName();//Prepping file name for next event
			}else{//Start fresh with a new file
				if(sum_q_segment > max_q_seg){
					max_q_seg = sum_q_segment;
				}
				if( (float)run_seg_size > max_event_seg){
					max_event_seg = run_seg_size;
				}
				if(sum_q_segment < min_q_seg && sum_q_segment>0.0){
					min_q_seg = sum_q_segment;
				}
				if( (float)run_seg_size < min_event_seg && (float)run_seg_size>0.0){
					min_event_seg = run_seg_size;
				}
				hist.Fill_Norm_Seg(run_num,run_seg,sum_q_segment/(float)run_seg_size,run_nums);//Normalized Charge Yield for individual files
				hist.Fill_Charge_v_Event_Seg(sum_q_segment,(float)run_seg_size);
				sum_q_segment = 0.0;
				run_seg_size = 1;
				if(run_num == fun::run_number(chain->GetFile()->GetName(),remove_front)){//Still in same run, but different file
					run_num_size+=1;
				}else{//Different Run
					if(sum_q_run > max_q_run){
						max_q_run = sum_q_run;
					}
					if( (float)run_num_size > max_event_run){
						max_event_run = run_num_size;
					}
					if(sum_q_run < min_q_run && sum_q_run>0.0){
						min_q_run = sum_q_run;
					}
					if( (float)run_num_size < min_event_run && (float)run_num_size>0.0){
						min_event_run = run_num_size;
					}
					hist.Fill_Norm_Run(run_num,sum_q_run/(float)run_num_size);//Normalized Charge Yield for Individual Runs
					hist.Fill_Charge_v_Event_Run(sum_q_run,(float)run_num_size);
					run_num_size=1;
					sum_q_run = 0.0;
				}
				old_max_q = data->Branches::q_l();//Prepping charge for next event
				old_file_name = chain->GetFile()->GetName();//Prepping file name for next event
				run_num = fun::run_number(chain->GetFile()->GetName(),remove_front);//Getting Run number
				run_seg = fun::run_segment(chain->GetFile()->GetName(),remove_front,remove_mid);//Getting run segment
				if(fun::run_num_idx(run_num,run_nums)==-1){
					run_nums.push_back(run_num);
				}
				//std::cout<<"\tNew File starting run:" <<run_num << " seg:" <<run_seg <<" charge: " <<old_max_q <<"\n";
			}
			if(curr_event%(num_of_events/100) == 0){
				std::cout<<"\r" <<"\t" <<(100*curr_event/num_of_events) <<" %"  <<std::flush;//<<"|| File: " <<chain->GetFile()->GetName() <<std::flush;//;
			}
	}
	std::cout<<"\nTotal Charge: " <<sum_q <<"\n";
	std::cout<<"For Individual Files\n";
	std::cout<<"\tMin Charge: " <<min_q_seg <<" Max Charge: " <<max_q_seg <<"\n";
	std::cout<<"\tMin Event: " <<min_event_seg <<" Max Event: " <<max_event_seg <<"\n";
	std::cout<<"For Individual Runs\n";
	std::cout<<"\tMin Charge: " <<min_q_run <<" Max Charge: " <<max_q_run <<"\n";
	std::cout<<"\tMin Event: " <<min_event_run <<" Max Event: " <<max_event_run <<"\n";
	hist.Write(run_nums,output_name);
	std::cout<<"Complete\n";
	//Timer and Efficiency Counters
	std::cout.imbue(std::locale(""));//Puts commas in appropriately? Apparently?
	std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now()-start);
	std::cout<< elapsed_full.count() << " Sec" <<std::endl;//Total elapsed time
	//std::cout<< events/elapsed_full.count() <<" Hz" <<std::endl; 
	return 0; 

}