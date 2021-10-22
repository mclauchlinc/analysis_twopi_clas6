#include "histogram.hpp"


Histogram::Histogram(std::string run_group_){
	if(run_group_ == "e16"){
		_run = 0;
	}else if(run_group_ == "e1f"){
		_run = 1;
	}else{
		_run = -1;
	}
	Histogram::Make();
}

//Histogram::~Histogram() { this->Write(); }

void Histogram::Make(){
	Histogram::Make_Norm_Seg();
	Histogram::Make_Norm_Run();
}

void Histogram::Write(std::vector<int> runs_,std::string output_file_){
	std::cout<<"Writing Histograms to " <<output_file_ <<"\n";
	_RootOutputFile = fun::Name_File(output_file_);
	Norm_Seg_Write(runs_);
	Norm_Run_Write();
}

void Histogram::Make_Norm_Seg(){
	std::cout<<"\tMake Normalized Charge Hists by Segment\n";
	char hname[100];
	sprintf(hname,"Normalized Faraday Charge Distribution by File");
	_Norm_Seg_Int = new TH1F(hname,hname,200,_Norm_Seg_Int_Bounds[_run][0],_Norm_Seg_Int_Bounds[_run][1]);
}

void Histogram::Make_Norm_Run(){
	std::cout<<"\tMake Normalized Charge Hists by Run\n";
	char hname[100];
	sprintf(hname,"Normalized Faraday Charge Distribution by Run");
	_Norm_Run_Int = new TH1F(hname,hname,200,_Norm_Run_Int_Bounds[_run][0],_Norm_Run_Int_Bounds[_run][1]);
	sprintf(hname,"Normalized Faraday Charge by Run");
	_Norm_Run = new TH1I(hname,hname,_Norm_Run_Max[_run]-_Norm_Run_Min[_run],_Norm_Run_Min[_run],_Norm_Run_Max[_run]);
}


void Histogram::Fill_Norm_Seg(int run_num_, int run_seg_, float norm_q_, std::vector<int> runs_ ){
	//std::cout<<"Trying to fill Norm Segment for run:" <<run_num_ <<" at seg:" <<run_seg_ <<"\n";
	char hname[100];
	if(fun::run_num_idx(run_num_,runs_) == -1 || _Norm_Seg.size()!=runs_.size()){
		//std::cout<<"\tNot registered run number yet, so filling index: ";
		sprintf(hname,"Normalized Faraday Charge for Run:%d",run_num_);
		//std::cout<<"\ttangent: " <<_Norm_Seg_Bounds[1]-_Norm_Seg_Bounds[0] <<"\n";
		_Norm_Seg.push_back(new TH1I(hname,hname,_Norm_Seg_Bounds[1]-_Norm_Seg_Bounds[0],_Norm_Seg_Bounds[0],_Norm_Seg_Bounds[1]));
		//std::cout<<_Norm_Seg.size()-1 <<"\n";
		_Norm_Seg[_Norm_Seg.size()-1]->Fill(run_seg_,norm_q_);
	}else{
		//std::cout<<"\tHist size " <<_Norm_Seg.size() <<" Filling index " <<fun::run_num_idx(run_num_,runs_) <<"\n";
		_Norm_Seg[fun::run_num_idx(run_num_,runs_)]->Fill(run_seg_,norm_q_);
	}
	//std::cout<<"Filling Norm Seg Int\n";
	_Norm_Seg_Int->Fill(norm_q_);
}

void Histogram::Fill_Norm_Run(int run_num_, float norm_q_){
	//std::cout<<"Filling Norm Run: " <<run_num_ <<"  with norm_q:" <<norm_q_ <<"\n";
	_Norm_Run->Fill(run_num_,norm_q_);
	_Norm_Run_Int->Fill(norm_q_);
}

void Histogram::Write(const std::string& output_file, std::vector<int> runs_){
	_RootOutputFile = fun::Name_File(output_file);
	std::cout<< "Writing Plots:\n";
	_RootOutputFile->cd();
	Histogram::Norm_Seg_Write(runs_);
	Histogram::Norm_Run_Write();
	std::cout<<"\tHistograms Done!\n";
}


void Histogram::Norm_Seg_Write(std::vector<int> runs_){
	std::cout<<"\tWriting Segment Histograms\n";
	TDirectory* seg_dir = _RootOutputFile->mkdir("Normalized Charge By File");
	seg_dir->cd();
	_Norm_Seg_Int->SetXTitle("Normalized Charge");
	_Norm_Seg_Int->SetYTitle("Number of Files");
	_Norm_Seg_Int->Write();
	TDirectory* seg_dir_sub = seg_dir->mkdir("Normalized Charge By File By Run");
	seg_dir_sub->cd();
	for(int i=0; i<runs_.size(); i++){
		_Norm_Seg[i]->SetXTitle("Segment");
		_Norm_Seg[i]->SetYTitle("Normalized Charge");
		_Norm_Seg[i]->Write();
	}
}

void Histogram::Norm_Run_Write(){
	std::cout<<"\tWriting Run Histogram\n";
	TDirectory* run_dir = _RootOutputFile->mkdir("Normalized Charge by Run");
	run_dir->cd();
	_Norm_Run->SetXTitle("Run Number");
	_Norm_Run->SetYTitle("Normalized Charge");
	_Norm_Run->Write();
	_Norm_Run_Int->SetYTitle("Number of Runs");
	_Norm_Run_Int->SetXTitle("Normalized Charge");
	_Norm_Run_Int->Write();
}
