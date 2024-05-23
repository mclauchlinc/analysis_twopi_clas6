#include "main.hpp"
#include <future>
#include <thread>
#include "TROOT.h"



//Main body
int main(int argc, char **argv){

	//Apparently one needs this to ensure R0t doesn't break with multiple threads
	ROOT::EnableThreadSafety();

	//start timer
	std::cout<<"Starting Clock\n";
	auto start = std::chrono::high_resolution_clock::now();
	//std::cout<<"The res for DTx is: " <<DTxres <<" and the rest for DTy is: " <<DTyres <<std::endl;
	
	
	//std::vector<std::vector<std::string>> infilenames2(_NUM_THREADS_);//For two lists of plate stuff

	std::cout<<"Reading Flags\n";
	auto flags = std::make_shared<Flags>();//argc,argv);
	//std::cout<<"Output Flags before filling\n";
	//flags->Flags::Print_Flags();
	flags->Flags::Read_Flags(argc,argv);
	std::cout<<"Output Flags after filling\n";
	flags->Flags::Print_Flags();
	//std::vector<std::vector<std::string>> infilenames(flags->Flags::Num_Cores());

	int num_files = flags->Flags::Num_Files();
	std::cout<<"Looking to have " <<flags->Flags::Num_Cores() <<" cores being used\n";
	//for(int i=0; i<flags->Flags::Num_Cores(); i++){
	//	infilenames[i]  = fun::read_file_list(flags->Flags::Files(),i,flags);
	//}
	
	// Make a set of threads (Futures are special threads which return a value)
  	std::future<size_t> threads[flags->Flags::Num_Cores()];

	//Define variable for total number of events
	size_t events = 0; 

	std::shared_ptr<double> q_tot;

	for(int i=0; i<std::distance(std::begin(_ecuts_), std::end(_ecuts_)); i++){
		std::cout<<"ecut:" <<_ecuts_[i] <<" cut:" <<fun::ecut_perform(_ecuts_[i],flags) <<" mod_idx:" <<fun::ecut_idx(_ecuts_[i])+fun::ecut_offset(_ecuts_[i],flags) <<"\n";
	}

	//Make histograms objects as a shared pointer that all threads will have
	std::cout<<"Making Histogram\n";
	auto hists = std::make_shared<Histogram>(flags);//Check on this

	//Make relevant TTrees and Event Rootfile
	//std::cout<<"Making TTree\n";
	//auto forest = std::make_shared<Forest>(flags); 
	for(int i=0; i<5; i++){
		std::cout<<"power i:" <<i <<" = " <<corr::power(2.0, i) <<"\n";
	}

	std::future<bool> fut;


	std::cout<<"Am i making bin centering plots? " <<flags->Flags::Plot_Bin_Centering() <<"\n";
	//For each thread
	std::cout<<"Running Multithreaded\n";
	for(int i = 0; i<flags->Flags::Num_Cores(); i++){
		//Set the thread to run asynchronously
		//The function running is the first argument
		//The functions arguments are all remaining arguments
		threads[i] = std::async(run_files, hists, i, file_num, q_tot,flags);//, num_mixed_p_pip[i]);infilenames.at(i)
		/*
		I keep running into issues with putting all the files onto one chain
		It says "illegal instruction" when trying to GetEntries() for the chain
		if all the simulation files are on one chain. 
		Further, THnSparse don't like multiple threads writing to them at the same
		time. So to resolve this, I'm still having multithreaded, but now waiting
		for each thread to finish first so only one thread is writing to the 
		THnSparse at a time. This will be very slow
		*/
		//if(flags->Flags::Make_Friend()){
		//	threads[i].wait(); 
		//}
	}
	for(int j = 0; j<flags->Flags::Num_Cores(); j++){//_NUM_THREADS_
		threads[j].wait(); //Wait until all threads are complete to move on
	}

	//For each thread to see how many events each thread successfully analyized
	for(int i = 0; i<flags->Flags::Num_Cores(); i++){//_NUM_THREADS_
		events += threads[i].get();
	}
	std::cout<<"\tClean Event Dist:\n";
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<_top_[i] <<":" <<hists->Histogram::Clean_Top(i);
	}
	std::cout<<"\n\tIsolated Event Dist:\n";
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<_top_[i] <<":" <<hists->Histogram::Isolated_Top(i);
	}
	std::cout<<"\n\tClean not Isolated Event Dist:\n";
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<_top_[i] <<":" <<hists->Histogram::Clean_Not_Isolated_Top(i);
	}
	std::cout<<"\n\tClean and Isolated Event Dist:\n";
	for(int i=0; i<4; i++){
		std::cout<<"\t" <<_top_[i] <<":" <<hists->Histogram::Clean_and_Isolated_Top(i);
	}
	//std::cout<<std::endl <<"Total Number of Files: " <<envi->Environment::was_num_file() <<std::endl; 
	std::cout<<"\n***Integrated Charge was " <<q_tot <<"**\n";

	long total_passed_events[3][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0}};
	long total_bad_angles[3][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0}};
	for(int i=0; i<4; i++){
		std::cout<<topologies[i+1] <<" measured:" <<hists->Histogram::NTop(i) <<" missed:" <<hists->Histogram::NPTop(i) <<"\n";
		for(int l=0; l<3; l++){
			for(int j=0; j<29; j++){
				for(int k=0; k<5; k++){
					total_passed_events[l][i]+= hists->Histogram::N_Yield(l,j,k,i);
				}
			}
			std::cout<<"\tin defined range for var " <<l <<": " <<total_passed_events[l][i] <<"\n";
			std::cout<<"\tBad Theta Angles in var " <<l <<":" <<hists->Bad_Angles(l,i,0) <<"\n";
			std::cout<<"\tBad Alpha Angles in var " <<l <<":" <<hists->Bad_Angles(l,i,1) <<"\n";
			std::cout<<"\tBad Phi Angles in var " <<l <<":" <<hists->Bad_Angles(l,i,2) <<"\n";
			std::cout<<"*Events not passed in var " <<l <<":" <<hists->Ev_No_Pass(l,i) <<"\n";
		}
		
	}
	std::cout<<"\nWriting Histograms\n";
	hists->Histogram::Write(flags);
	//hists->Histogram::Print(output_name,envi);

	std::cout<<"\n\nOutput File:\n\t" <<flags->Flags::Output_Name() <<"\n\n";

	std::cout<<std::endl;
	//Timer and Efficiency Counters
	std::cout.imbue(std::locale(""));//Puts commas in appropriately? Apparently?
	std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now()-start);
	std::cout<< elapsed_full.count() << " Sec" <<std::endl;//Total elapsed time
	std::cout<< events/elapsed_full.count() <<" Hz" <<std::endl; 
	return 0; 

}