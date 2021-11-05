#ifndef HISTOGRAM_H_GUARD
#define HISTOGRAM_H_GUARD

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "constants.hpp"
#include "functions.hpp"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include "TGraph.h"
//#include <unistd.h>//to allow us to use chdir
//#include "TImage.h"
//#include "particle.hpp"
//#include "variables.h"
//#include "CartesianGenerator.hh"


class Histogram {
protected:
	std::shared_ptr<TFile> _RootOutputFile;

	int _run = -1;

	int _Norm_Run_Max[2] = {31500,51000};//Need to settle e1f numbers
	int _Norm_Run_Min[2] = {30500,50000};//Need to settle e1f numbers
	int _Norm_Seg_Bounds[2] = {0,20};
	float _Norm_Seg_Int_Bounds[2][2] = {{0.000020,0.000060},{0.000020,0.000060}};
	float _Norm_Run_Int_Bounds[2][2] = {{0.000020,0.000060},{0.000020,0.000060}};
	float _Norm_Seg_Int_Inv_Bounds[2][2] = {{20000.0,30000.},{20000.0,30000.}};
	float _Norm_Run_Int_Inv_Bounds[2][2] = {{20000.0,30000.},{20000.0,30000.}};
	float _Charge_Seg_Bounds[2] = {0.050,60.0};
	float _Event_Seg_Bounds[2] = {400.0,110000.0};
	float _Charge_Run_Bounds[2] = {0.5,350.0};
	float _Event_Run_Bounds[2] = {1500.0,1500000.0};

	std::vector<TH1I*> _Norm_Seg;//Normalized Integrated Charge by Run and Segment
	TH1I* _Norm_Run;//Normalized Integrated Charge by Run
	TH1F* _Norm_Seg_Int;//Normalized Integrated Charge Distribution across all Runs and Segments
	TH1F* _Norm_Run_Int;//Normalized Integrated Charge Distribution across all Runs
	std::vector<TH1I*> _Norm_Seg_Inv;//Normalized Integrated Charge by Run and Segment
	TH1I* _Norm_Run_Inv;//Normalized Integrated Charge by Run
	TH1F* _Norm_Seg_Int_Inv;//Normalized Integrated Charge Distribution across all Runs and Segments
	TH1F* _Norm_Run_Int_Inv;//Normalized Integrated Charge Distribution across all Runs
	TH2F* _Charge_v_Events_Seg;
	TH2F* _Charge_v_Events_Run;


public:
	Histogram(std::string run_group_);
	void Make();
	void Write(std::vector<int> runs_,std::string output_file_);
	void Make_Charge_v_Event();
	void Make_Norm_Seg();
	void Make_Norm_Run();
	void Fill_Norm_Seg(int run_num_, int run_seg_, float norm_q_, std::vector<int> runs_ );
	void Fill_Norm_Run(int run_num_, float norm_q_);
	void Fill_Charge_v_Event_Seg(float charge_, float event_);
	void Fill_Charge_v_Event_Run(float charge_, float event_);
	void Write(const std::string& output_file, std::vector<int> runs_);
	void Norm_Seg_Write(std::vector<int> runs_);
	void Norm_Run_Write();

	
};





#endif