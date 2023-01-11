#include "main.hpp"


int main(int argc, char **argv){
	//start the clock
	auto start = std::chrono::high_resolution_clock::now();

	auto flags = Flags();//argc,argv);
	flags.Flags::Read_Flags(argc,argv);
	std::cout<<"Test Plot Single Diff " <<flags.Flags::Plot_Single_Diff() <<"\n";

	

	std::cout<<"Run both? " <<flags.Flags::Run("both") <<"\n";
	std::cout<<"Run both? " <<flags.Flags::Run() <<"\n";
	if(flags.Flags::Run("both") == _both_ ){//Looking at both e16 and e1f results
		if(flags.Flags::Helicity()){
			std::cout<<"Working with both e16 and e1f\n";
			std::cout<<"With Helicity Tracking";
			exp_pos_file = new TFile(flags.Flags::Exp_Hel(1).c_str());
			exp_neg_file = new TFile(flags.Flags::Exp_Hel(-1).c_str());
			exp_pos_file2 = new TFile(flags.Flags::Exp_Hel(1).c_str());
			exp_neg_file2 = new TFile(flags.Flags::Exp_Hel(-1).c_str());
			sim_file = new TFile(flags.Flags::Sim_File().c_str());
			sim_file2 = new TFile(flags.Flags::Sim_File2().c_str());
			if(flags.Flags::Rad_Corr()){
				sim_no_rad_file = new TFile(flags.Flags::Sim_No_Rad().c_str());
				sim_no_rad_file2 = new TFile(flags.Flags::Sim_No_Rad2().c_str());
			}
		}else{
			std::cout<<"Working with both e16 and e1f\n";
			exp_file = new TFile(flags.Flags::Exp_File().c_str());
			sim_file = new TFile(flags.Flags::Sim_File().c_str());
			exp_file2 = new TFile(flags.Flags::Exp_File2().c_str());
			sim_file2 = new TFile(flags.Flags::Sim_File2().c_str());
			if(flags.Flags::Rad_Corr()){
				sim_no_rad_file = new TFile(flags.Flags::Sim_No_Rad().c_str());
				sim_no_rad_file2 = new TFile(flags.Flags::Sim_No_Rad2().c_str());
			}
		}
		if(flags.Flags::Has_Empty()){
			empty_file = new TFile(flags.Flags::Empty_File().c_str());
			empty_file2 = new TFile(flags.Flags::Empty_File2().c_str());
			//auto hist Histogram(flags.Flags::Output_File(),exp_file,sim_file,empty_file,exp_file2,sim_file2,empty_file2,flags);
		}else{
			//auto hist Histogram(flags.Flags::Output_File(),exp_file,sim_file,exp_file2,sim_file2,flags);
		}
		if(flags.Flags::Has_Weight()){
			weight_file = new TFile(flags.Flags::Weight_File().c_str());
			weight_file2 = new TFile(flags.Flags::Weight_File2().c_str());
		}
		//auto hist = Histogram(flags.Flags::Output_File(),exp_file,sim_file,flags);
	}else if(flags.Flags::Has_Exp() && flags.Flags::Has_Sim()){//Looking at just e16 or e1f
		std::cout<<"Working with just " <<flags.Flags::Run("type") <<"\n";
		if(flags.Flags::Helicity()){
			std::cout<<"Helcity stuff!\n";
			exp_file = new TFile(flags.Flags::Exp_Hel(1).c_str());
			//exp_neg_file = new TFile(flags.Flags::Exp_Hel(-1).c_str());
			sim_file = new TFile(flags.Flags::Sim_File().c_str());
			if(flags.Flags::Rad_Corr()){
				sim_no_rad_file = new TFile(flags.Flags::Sim_No_Rad().c_str());
			}
			if(flags.Flags::Has_Empty()){
				empty_file = new TFile(flags.Flags::Empty_File().c_str());
				std::cout<<"Running single run with empty and no rad\n";
				std::cout<<"empty file: " <<flags.Flags::Empty_File().c_str() <<"\nno rad file: " <<flags.Flags::Sim_No_Rad().c_str() <<"\n";
				auto hist = Histogram(flags.Flags::Output_File(),exp_file,sim_file,empty_file,sim_no_rad_file,flags);
			}else{
				//auto hist = Histogram(flags.Flags::Output_File(),exp_file,sim_file,flags);
			}
		}else{
			std::cout<<"No Helcity stuff!\n";
			exp_file = new TFile(flags.Flags::Exp_File().c_str());
			sim_file = new TFile(flags.Flags::Sim_File().c_str());
			if(flags.Flags::Rad_Corr()){
				sim_no_rad_file = new TFile(flags.Flags::Sim_No_Rad().c_str());
				//sim_no_rad_file2 = new TFile(flags.Flags::Sim_No_Rad2());
			}
			//auto hist = Histogram(flags.Flags::Output_File(),exp_file,sim_file,flags);
		}
		
		
		//Haven't implemented empty file or cc_acceptance integration yet
		/*if(flags.Flags::Has_Weight()){
			weight_file = new TFile(flags.Flags::Weight_File().c_str());
		}
		if(flags.Flags::Has_Empty()){
			empty_file = new TFile(flags.Flags::Empty_File().c_str());
			//auto hist Histogram(flags.Flags::Output_File(),exp_file,sim_file,empty_file,flags);
		}else if(flags.Flags::Has_Weight()){
			auto hist = Histogram(flags.Flags::Output_File(),exp_file,sim_file,flags);
		}else{
			auto hist = Histogram(flags.Flags::Output_File(),exp_file,sim_file,flags);
		}*/
		
	}
}
