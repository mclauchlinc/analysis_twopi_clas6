#include "main.hpp"


int main(int argc, char **argv){
	//start the clock
	auto start = std::chrono::high_resolution_clock::now();

	auto flags = Flags();//argc,argv);
	flags.Flags::Read_Flags(argc,argv);

	exp_file = new TFile(flags.Flags::Exp_Hel().c_str());
    sim_file = new TFile(flags.Flags::Sim_File().c_str());
    sim_no_rad_file = new TFile(flags.Flags::Sim_No_Rad().c_str());
    empty_file = new TFile(flags.Flags::Empty_File().c_str());
    auto hist = Histogram(exp_file,sim_file,empty_file,sim_no_rad_file,flags);
    

}
