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

	//for threading

	std::vector<std::vector<std::string>> infilenames(_NUM_THREADS_);
	//std::vector<std::vector<std::string>> infilenames2(_NUM_THREADS_);//For two lists of plate stuff

	std::cout<<"Reading Flags\n";
	auto flags = std::make_shared<Flags>();//argc,argv);
	flags->Flags::Read_Flags(argc,argv);

	std::cout<<"Test for PIM Fid plotting "<<flags->Plot_Fid(3) <<"\n";
	std::cout<<"Test for perform  "<<fun::ecut_perform(_event_,flags) <<"\n";


	int num_files = flags->Flags::Num_Files();
	for(int i=0; i<_NUM_THREADS_; i++){
		infilenames[i]  = fun::read_file_list(flags->Flags::Files(),i);
	}
	
	// Make a set of threads (Futures are special threads which return a value)
  	std::future<size_t> threads[_NUM_THREADS_];

	//Define variable for total number of events
	size_t events = 0; 



	//Make histograms objects as a shared pointer that all threads will have
	std::cout<<"Making Histogram\n";
	auto hists = std::make_shared<Histogram>(flags);//Check on this

	//Make relevant TTrees and Event Rootfile
	std::cout<<"Making TTree\n";
	auto forest = std::make_shared<Forest>(flags); 

	std::future<bool> fut;

	//For each thread
	std::cout<<"Running Multithreaded\n";
	for(int i = 0; i<_NUM_THREADS_; i++){
		//Set the thread to run asynchronously
		//The function running is the first argument
		//The functions arguments are all remaining arguments
		threads[i] = std::async(run_files, infilenames.at(i), hists, forest, i, file_num, flags);//, num_mixed_p_pip[i]);
	}
	for(int j = 0; j<_NUM_THREADS_; j++){
		threads[j].wait(); //Wait until all threads are complete to move on
	}

	//For each thread to see how many events each thread successfully analyized
	for(int i = 0; i<_NUM_THREADS_; i++){
		events += threads[i].get();
	}
	//std::cout<<std::endl <<"Total Number of Files: " <<envi->Environment::was_num_file() <<std::endl; 
	
	std::cout<<"\nWriting Histograms\n";
	hists->Histogram::Write(flags);
	//hists->Histogram::Print(output_name,envi);

	std::cout<<std::endl;
	//Timer and Efficiency Counters
	std::cout.imbue(std::locale(""));//Puts commas in appropriately? Apparently?
	std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now()-start);
	std::cout<< elapsed_full.count() << " Sec" <<std::endl;//Total elapsed time
	std::cout<< events/elapsed_full.count() <<" Hz" <<std::endl; 
	return 0; 

}