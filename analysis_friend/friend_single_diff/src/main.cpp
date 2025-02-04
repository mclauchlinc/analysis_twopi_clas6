#include "main.hpp"


int main(int argc, char **argv){
	//start the clock
	auto start = std::chrono::high_resolution_clock::now();

	auto flags = Flags();//argc,argv);
	flags.Flags::Read_Flags(argc,argv);


    std::unique_ptr<TFile> exp_file(new TFile(flags.Flags::Exp_Hel().c_str()));
    std::unique_ptr<TFile> sim_file(new TFile(flags.Flags::Sim_File().c_str()));
    //std::unique_ptr<TFile> sim_no_rad_file(new TFile(flags.Flags::Sim_No_Rad().c_str()));
    std::unique_ptr<TFile> empty_file(new TFile(flags.Flags::Empty_File().c_str()));
    /*
	auto exp_file = std::unique_ptr<TFile>(flags.Flags::Exp_Hel().c_str());
    auto sim_file = std::unique_ptr<TFile>(flags.Flags::Sim_File().c_str());
    auto sim_no_rad_file = std::unique_ptr<TFile>(flags.Flags::Sim_No_Rad().c_str());
    auto empty_file = std::unique_ptr<TFile>(flags.Flags::Empty_File().c_str());
    */
    if(flags.Flags::Has_Localized_Holes()){
        //auto hole_file = std::unique_ptr<TFile>(flags.Flags::Holes_File().c_str());
        std::unique_ptr<TFile> hole_file(new TFile(flags.Flags::Holes_File().c_str()));
        auto hist = Histogram(exp_file,sim_file,empty_file,hole_file,flags);
        //auto hist = Histogram(exp_file,sim_file,empty_file,sim_no_rad_file,hole_file,flags);
    }else{
        auto hist = Histogram(exp_file,sim_file,empty_file,flags);
        //auto hist = Histogram(exp_file,sim_file,empty_file,sim_no_rad_file,flags);
    }
    
    std::cout<<"Output File: " <<flags.Flags::Output_File() <<"\n";

}
