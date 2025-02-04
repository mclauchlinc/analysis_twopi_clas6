#ifndef EVENT_ANALYSIS_HPP
#define EVENT_ANALYSIS_HPP

#include "TLorentzVector.h"
#include <iostream>
#include "histogram.hpp"
#include "constants.hpp"
#include "branches.hpp"
#include "cuts.hpp"
#include "physics.hpp"
#include "forest.hpp"
//#include "ntuple.hpp"
#include "functions.hpp"
#include "particle.hpp"
#include "event.hpp"
#include "plot.hpp"

/*
This class is designed to perform the bulk of the analysis
- Particle Identification
- Event Selection
- Histogram Filling
- THnSparse Filling
Basically setting everyting up for pushing out into the world 
*/

class Analysis{
private:
	int _gevts = 0; //Number of good events
	int _npart = 0;//Number of particles measured
	int _gtop[4] = {0,0,0,0}; //Number of good events for each topology
	int _ntop[4] = {0,0,0,0}; //Number of candidate events for each topology
	int _gpart[4] = {0,0,0,0}; //Number of identified particles of each species
	std::vector<int> _pro_idx;
	std::vector<int> _pip_idx;
	std::vector<int> _pim_idx;
	std::vector<std::vector<int>> _idx_mpro;
	std::vector<std::vector<int>> _idx_mpip;
	std::vector<std::vector<int>> _idx_mpim;
	std::vector<std::vector<int>> _idx_mzero;

	long _number_of_events = 0;
	long _number_of_recon = 0;

	long _number_of_events_pos = 0;
	long _number_of_recon_pos = 0;

	long _number_of_events_neg = 0;
	long _number_of_recon_neg = 0;

	int _hel = 0;

	std::vector<Particle> _tParticle; //Just 4
	std::vector<Particle> _rParticle;	//Reconstructed Particles
	std::vector<Event> _rEvent;	//Reconstructed Events
	std::vector<Event> _tEvent; //Just 1
	std::vector<Event> _gEvent;	//Good events in a given topology
	std::vector<Event> _iEvent; //Just 1

	std::vector<int> _gevt_idx;
	std::vector<int> _gevt_idx_mpro;
	std::vector<int> _gevt_idx_mpip;
	std::vector<int> _gevt_idx_mpim;
	std::vector<int> _gevt_idx_mzero;

	int _gevents = 0; 

	int _set = -1; //{0,1} -> {e16,e1f}
	bool _sim = false;
	bool _recon = false;

	float _weight = 1.0;//with event considerations
	float _pweight = 1.0;//Pre-event considerations

	int _top_passed = -1;



public:
	Analysis(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, int thread_id_, int run_num_, std::shared_ptr<Flags> flags_);
	//(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Environment> envi_, int run_type_, std::shared_ptr<forest> a_forest_, int thread_id_, int run_num_);
	//Number of potential events given identified particles
	void Num_top(std::shared_ptr<Flags> flags_);
	//Particle ID
	void PID(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	//Event ID
	void Event_ID(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	//Isolate a Chosen Event for Analysis
	void Isolate_Event(std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	//Isolating a Chosen Event from a given Topology
	void Isolate_Top(int top_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	//Getting the index of an event given its topology and top index
	int Event_idx(int top_, int top_idx);
	int gEvent_idx(int top_, int top_idx_);
	//Number of good events
	//int Gevts();
	//Number of possible events in given topology
	//int Ntop(int i);
	//Determine whether half wave plate is in or out
	bool Half_Wave(int run_num_, std::shared_ptr<Flags> flags_);
	//Correct Helicity accoridng to half wave plate status
	float Corr_Helicity(float helicity_, int run_num_, std::shared_ptr<Flags> flags_);
	void Plot_Particles(std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	void Plot_Events(std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	bool Valid_Event();



};


#endif