#include "main.hpp"


int main(int argc, char **argv){
	//start the clock
	auto start = std::chrono::high_resolution_clock::now();

	auto flags = Flags();//argc,argv);
	flags.Flags::Read_Flags(argc,argv);

	exp_file = new TFile(flags.Flags::Exp_Hel().c_str());
    sim_file = new TFile(flags.Flags::Sim_File().c_str());
    //sim_no_rad_file = new TFile(flags.Flags::Sim_No_Rad().c_str());
    //empty_file = new TFile(flags.Flags::Empty_File().c_str());
    auto hist = Histogram(exp_file,sim_file,flags);
    

	std::cout<<"\nFlags:\n";
	std::cout<<"Sim File: " <<flags.Flags::Sim_File() <<"\n";
	std::cout<<"Exp File: " <<flags.Flags::Exp_File() <<"\n";
	std::cout<<"Top : " <<flags.Flags::Top() <<"\n";
	std::cout<<"Top idx: " <<flags.Flags::Top_idx() <<"\n";
	//std::cout<<"Read Top : " <<flags.Flags::Read_Top() <<"\n";
	//std::cout<<"Read Top idx: " <<flags.Flags::Read_Top_idx() <<"\n";
	std::cout<<"Var: " <<flags.Flags::Var_Set() <<"\n";
	std::cout<<"Var idx: " <<flags.Flags::Var_idx() <<"\n";
	std::cout<<"Helicity: " <<flags.Flags::Helicity() <<"\n";
	std::cout<<"AREC: " <<flags.Flags::Acc_Rel_Error_Cut() <<"\n";
	std::cout<<"Min Dist: " <<flags.Flags::Min_Local_Dist() <<"\n";
	std::cout<<"\nOutput File: " <<flags.Flags::Output_File() <<"\n\n";

}
