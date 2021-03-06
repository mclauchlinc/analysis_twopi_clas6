#include "main.hpp"


int main(int argc, char **argv){
	//start the clock
	auto start = std::chrono::high_resolution_clock::now();

	auto flags = Flags();//argc,argv);
	flags.Flags::Read_Flags(argc,argv);
	std::cout<<"Test Plot Single Diff " <<flags.Flags::Plot_Single_Diff() <<"\n";

	TFile *exp_file;
	TFile *sim_file;
	TFile *empty_file;
	TFile *exp_file2;
	TFile *sim_file2;
	TFile *empty_file2;
	if(flags.Flags::Run("both") == _both_ ){//Looking at both e16 and e1f results
		exp_file = new TFile(flags.Flags::Exp_File().c_str());
		sim_file = new TFile(flags.Flags::Sim_File().c_str());
		exp_file2 = new TFile(flags.Flags::Exp_File2().c_str());
		sim_file2 = new TFile(flags.Flags::Sim_File2().c_str());
		if(flags.Flags::Has_Empty()){
			empty_file = new TFile(flags.Flags::Empty_File().c_str());
			empty_file2 = new TFile(flags.Flags::Empty_File2().c_str());
			//auto hist Histogram(flags.Flags::Output_File(),exp_file,sim_file,empty_file,exp_file2,sim_file2,empty_file2,flags);
		}else{
			//auto hist Histogram(flags.Flags::Output_File(),exp_file,sim_file,exp_file2,sim_file2,flags);
		}
	}else if(flags.Flags::Has_Exp() && flags.Flags::Has_Sim()){//Looking at just e16 or e1f
		exp_file = new TFile(flags.Flags::Exp_File().c_str());
		sim_file = new TFile(flags.Flags::Sim_File().c_str());
		if(flags.Flags::Has_Empty()){
			empty_file = new TFile(flags.Flags::Empty_File().c_str());
			//auto hist Histogram(flags.Flags::Output_File(),exp_file,sim_file,empty_file,flags);
		}else{
			auto hist = Histogram(flags.Flags::Output_File(),exp_file,sim_file,flags);
		}
	}
}
