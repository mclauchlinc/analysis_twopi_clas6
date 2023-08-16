#include "main.hpp"


int main(int argc, char **argv){
	//start the clock
	auto start = std::chrono::high_resolution_clock::now();

	auto flags = Flags();//argc,argv);
	flags.Flags::Read_Flags(argc,argv);

    auto exp_file = new TFile(flags.Flags::Exp_Hel().c_str());
    auto sim_file = new TFile(flags.Flags::Sim_File().c_str());
    auto sim_no_rad_file = new TFile(flags.Flags::Sim_No_Rad().c_str());
    auto empty_file = new TFile(flags.Flags::Empty_File().c_str());

    auto hist = Histogram(exp_file,sim_file,empty_file,sim_no_rad_file,flags);
    /*
	exp_pim_file = new TFile(flags.Flags::Exp_Hel("PIM").c_str());
    sim_pim_file = new TFile(flags.Flags::Sim_File("PIM").c_str());
    sim_no_rad_pim_file = new TFile(flags.Flags::Sim_No_Rad("PIM").c_str());
    empty_pim_file = new TFile(flags.Flags::Empty_File("PIM").c_str());
    exp_pro_file = new TFile(flags.Flags::Exp_Hel("Pro").c_str());
    sim_pro_file = new TFile(flags.Flags::Sim_File("Pro").c_str());
    sim_no_rad_pro_file = new TFile(flags.Flags::Sim_No_Rad("Pro").c_str());
    empty_pro_file = new TFile(flags.Flags::Empty_File("Pro").c_str());
    exp_pip_file = new TFile(flags.Flags::Exp_Hel("PIP").c_str());
    sim_pip_file = new TFile(flags.Flags::Sim_File("PIP").c_str());
    sim_no_rad_pip_file = new TFile(flags.Flags::Sim_No_Rad("PIP").c_str());
    empty_pip_file = new TFile(flags.Flags::Empty_File("PIP").c_str());
    */
    /*
    if(flags.Flags::Has_Localized_Holes()){
        hole_file = new TFile(flags.Flags::Holes_File().c_str());
        auto hist = Histogram(exp_file,sim_file,empty_file,sim_no_rad_file,hole_file,flags);
    }else{
        auto hist = Histogram(exp_file,sim_file,empty_file,sim_no_rad_file,flags);
    }
    */

}
