#include "main.hpp"


int main(int argc, char **argv){
	std::cout<<"Started Program\n";
	//start the clock
	auto start = std::chrono::high_resolution_clock::now();

	/*
	Inputs should go
	./read_friend <experimental root file> <simulated root file w/ Radiative Effects> <simulated root file w/o Radiative Effects>
	*/
	if(argc != 4){//6){//When I have things for radiative corrections we can add this part
		std::cout<<"Hey, you didn't enter all the things required";
	}else{
		std::cout<<"\nExperiment File: " <<argv[1];
		std::cout<<"\nSimulated File: " <<argv[2];
		//std::cout<<"Non Rad Sim File:" <<argv[3];
		std::cout<<"\nOutput File Name: " <<argv[3];//[5]
	}

	//Read in the Various THnSparse Histograms
	TFile *exp_file = new TFile(argv[1],"READ");
	TFile *sim_file = new TFile(argv[2],"READ");

	std::string output_filename = argv[3];//[5]

	auto hist = Histogram(output_filename,exp_file,sim_file);
	
}
