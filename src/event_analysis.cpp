#include "event_analysis.hpp"
		  
//Analysis::Analysis(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Environment> envi_, int run_type_,std::shared_ptr<Forest> a_Forest_, int thread_id_, int run_num_){
Analysis::Analysis(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, int thread_id_, int run_num_, std::shared_ptr<Flags> flags_){
	if(thread_id_==0){
		//std::cout<<"Starting Event Analysis\n";
	}
	//int npart = 0; //# of Detected Particles
	if(flags_->Flags::Sim()){
		//std::cout<<"looking at Thrown particles\n";
		_npart = data_->Branches::npart();
		_hel = 1; 
		_weight = data_->Branches::weight();
		//Thrown Particle
		int thr_idx[4] = {-1,-1,-1,-1};
		for(int i=0; i<4; i++){
			_tParticle.push_back(Particle(i,data_,flags_,true));
			plot::plot_thrown(_tParticle[i],hist_,flags_);
			if(_tParticle[i].Particle::Is_Elec()){
				thr_idx[0] = i;
			}else if(_tParticle[i].Particle::Is_Pro()){
				thr_idx[1] = i;
			}else if(_tParticle[i].Particle::Is_Pip()){
				thr_idx[2] = i;
			}else if(_tParticle[i].Particle::Is_Pim()){
				thr_idx[3] = i;
			}
		}
		//Thrown Event
		
		_tEvent.push_back(Event(fun::top_idx(_mzero_),_tParticle[thr_idx[0]],_tParticle[thr_idx[1]],_tParticle[thr_idx[2]],_tParticle[thr_idx[3]],flags_,_weight,1,true));
		/*if(!_tEvent[0].Event::W() >= _W_min_ || !_tEvent[0].Event::W() < _W_max_){
			std::cout<<"Thrown W:" <<_tEvent[0].Event::W() <<"\n";
		}
		if(!_tEvent[0].Event::Q2() >= _Q2_bins_[0] || !_tEvent[0].Event::Q2() < _Q2_bins_[5]){
			std::cout<<"Thrown Q2:" <<_tEvent[0].Event::Q2() <<"\n";
		}*/
		/*double thr_W = _tEvent[0].Event::W();
		double thr_Q2 = _tEvent[0].Event::Q2();
		bool plot_thrown = true;
		plot_thrown &= (thr_W >= _W_min_);
		plot_thrown &= (thr_W >= _W_max_);
		plot_thrown &= (thr_Q2 >= _Q2_min_);
		plot_thrown &= (thr_Q2 < _Q2_max_);

		if(plot_thrown){	
			plot::plot_event(_tEvent[0],hist_,flags_,true);
		}*/
	}else{
		_npart = data_->Branches::gpart();
		_weight = 1.0;
		if(flags_->Flags::Helicity()){
			_hel = Analysis::Corr_Helicity(data_->Branches::evntclas2(),run_num_,flags_);//Need to pair with library of files where half-wave plate is in or out
		}
	}
	//std::cout<<"size of thrown vector " <<_tEvent.size() <<"\n";
	//Particle ID
	if(thread_id_==0){
		//std::cout<<"There are "<<_npart <<" good particles to ID\nBegin PID\n";
	}
	if(_npart>0){
		Analysis::PID(data_, hist_, flags_);//Identify and Plot Particles
		Analysis::Plot_Particles(hist_,flags_);
	}
	//Event ID
	if(thread_id_==0){
		//std::cout<<"Begin Event ID\n";
	}
	if(_npart>2){
		Analysis::Event_ID(data_, hist_, flags_);//Identify and Plot Events (and Event PID)
		if(thread_id_==0){
			//std::cout<<"Isolating Event\n";
		}
		Analysis::Isolate_Event(hist_, flags_);//Isolate Definitive Event and Plot Accordingly	
		Analysis::Plot_Events(hist_,flags_);
	}else if(flags_->Flags::Sim()){
		Plot_Events(hist_,flags_);
	}
}

void Analysis::Num_top(std::shared_ptr<Flags> flags_){
	//Pro Missing
	if(_rParticle.size()>0){
		if(_rParticle[0].Particle::Is_Elec()){
			//Proton Missing
			if(flags_->Flags::MM_Cut(0)){
				if(_pip_idx.size() > 0 && _pim_idx.size() > 0){
					std::vector<int> pro_tmp;
					_ntop[0] = _pip_idx.size() * _pim_idx.size();
					for(int i=0; i<_pip_idx.size(); i++){
						for(int j=0; j<_pim_idx.size(); j++){
							pro_tmp.push_back(_pip_idx[i]);
							pro_tmp.push_back(_pim_idx[j]);
							_idx_mpro.push_back(pro_tmp);
							pro_tmp.clear();
						}
					}
				}
			}
			//Pip Missing
			if(flags_->Flags::MM_Cut(1)){
				if(_pro_idx.size() > 0 && _pim_idx.size() > 0){
					std::vector<int> pip_tmp;
					_ntop[1] = _pim_idx.size() * _pro_idx.size();
					for(int i=0; i<_pro_idx.size(); i++){
						for(int j=0; j<_pim_idx.size(); j++){
							pip_tmp.push_back(_pro_idx[i]);
							pip_tmp.push_back(_pim_idx[j]);
							_idx_mpip.push_back(pip_tmp);
							pip_tmp.clear();
						}
					}
				}
			}
			//Pim Missing
			if(flags_->Flags::MM_Cut(2)){
				if(_pro_idx.size() > 0 && _pip_idx.size() > 0){
					std::vector<int> pim_tmp;
					for(int i=0; i<_pip_idx.size(); i++){
						for(int j=0; j<_pro_idx.size(); j++){
							if(_pro_idx[j] != _pip_idx[i]){
								_ntop[2]+=1;
								pim_tmp.push_back(_pro_idx[j]);
								pim_tmp.push_back(_pip_idx[i]);
								_idx_mpim.push_back(pim_tmp);
								pim_tmp.clear();
							}
						}
					}
				}
			}
			//Zero Missing
			if(flags_->Flags::MM_Cut(3)){
				if(_pro_idx.size() > 0 && _pip_idx.size() > 0 && _pim_idx.size() > 0){
					//std::cout<<"Looking at a Potential Exclusive topology\t";
					std::vector<int> zero_tmp;
					for(int i=0; i<_pip_idx.size(); i++){
						for(int j=0; j<_pro_idx.size(); j++){
							if(_pro_idx[j] != _pip_idx[i]){
								_ntop[3]+=1;
								for(int k=0; k<_pim_idx.size(); k++){
									zero_tmp.push_back(_pro_idx[j]);
									zero_tmp.push_back(_pip_idx[i]);
									zero_tmp.push_back(_pim_idx[k]);
									_idx_mzero.push_back(zero_tmp);
									zero_tmp.clear();
								}
							}
						}
					}
					_ntop[3] = _ntop[3]*_pim_idx.size();
					//std::cout<<"Found: " <<_ntop[3] <<" possibilities\n";
				}
			}
			//if(_ntop[0] + _ntop[1] +  _ntop[2] + _ntop[3]>=1 ){
		//		std::cout<<"\tNum Possible Topoogies: " <<_ntop[0] <<_ntop[1] <<_ntop[2] <<_ntop[3] <<"\n";
		//	}
		}
	}
}

void Analysis::PID(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_,std::shared_ptr<Flags> flags_){
	//std::cout<<"\t****IDing Possible Particle Events******\n";
	for(int j=0; j<_npart; j++){
		_rParticle.push_back(Particle(j,data_,flags_));//Fill and ID Particle
		plot::plot_pid(_rParticle[j],hist_,flags_);//Plot on relevant PID Plots
		if(_rParticle[0].Particle::Is_Elec()){
			for(int k=1; k<4; k++){//Figuring out which particles we have
				if(_rParticle[j].Particle::Pass_pid(k)){
					switch(k){
						case 1:
							_pro_idx.push_back(j);
						break;
						case 2:
							_pip_idx.push_back(j);
						break;
						case 3:
							_pim_idx.push_back(j);
						break;
						default:
							std::cout<<"\tNot a valid index for PID\n";
						break;
					}
				}
			}
		}
	}
	//std::cout<<"\t#Protons:" <<_pro_idx.size() <<"\n";
	//std::cout<<"\t#PIP:" <<_pip_idx.size() <<"\n";
	//std::cout<<"\t#PIM:" <<_pim_idx.size() <<"\n";
	//if(_pim_idx.size()>0){
	//	std::cout<<"\tnum pim poss:" <<_pim_idx.size() <<"\n";
	//}
	if(_pro_idx.size()*_pip_idx.size()>0 || _pro_idx.size()*_pim_idx.size()>0 || _pip_idx.size()*_pim_idx.size()>0 ){
		//std::cout<<"\tEvent PID\t# particles IDed: " <<_pro_idx.size() <<" " <<_pip_idx.size() <<" " <<_pim_idx.size() <<"\n";
	}
}

void Analysis::Event_ID(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	Analysis::Num_top(flags_);
	
	for(int i=0; i<4; i++){//Topologies
		if(_ntop[i]>0){//If there are any possible events for a given topology
			for(int j=0; j<_ntop[i]; j++){//loop through those possibilities 
				switch(i){
					case 0:
						//std::cout<<"\t\tidx: " <<i <<" pip idx:" <<_idx_mpro[j][0] <<" pim idx:" <<_idx_mpro[j][1] <<"\n";
						_rEvent.push_back(Event(i,_rParticle[0],_rParticle[_idx_mpro[j][0]],_rParticle[_idx_mpro[j][1]],flags_,_weight,_hel));
						//_ntop[0]+=1;
					break;
					case 1:
						//std::cout<<"\t\tidx: " <<i <<" pro idx:" <<_idx_mpip[j][0] <<" pim idx:" <<_idx_mpip[j][1] <<"\n";
						_rEvent.push_back(Event(i,_rParticle[0],_rParticle[_idx_mpip[j][0]],_rParticle[_idx_mpip[j][1]],flags_,_weight,_hel));
						//_ntop[1]+=1;
					break;
					case 2:
						//std::cout<<"\t\tidx: " <<i <<" pro idx:" <<_idx_mpim[j][0] <<" pip idx:" <<_idx_mpim[j][1] <<"\n";
						_rEvent.push_back(Event(i,_rParticle[0],_rParticle[_idx_mpim[j][0]],_rParticle[_idx_mpim[j][1]],flags_,_weight,_hel));
						//_ntop[2]+=1;
					break;
					case 3:
						//std::cout<<"\t\tidx: " <<i <<" pro idx:" <<_idx_mzero[j][0] <<" pip idx:" <<_idx_mzero[j][1] <<" pim idx:" <<_idx_mzero[j][2]<<"\n";
						_rEvent.push_back(Event(i,_rParticle[0],_rParticle[_idx_mzero[j][0]],_rParticle[_idx_mzero[j][1]],_rParticle[_idx_mzero[j][2]],flags_,_weight,_hel));
						//_ntop[3]+=1;
					break;
					default:
						std::cout<<"\tIncorrect Topology Search\n";
					break;
				}
			}
		}
	}
	
	for(int idx=0; idx<_rEvent.size(); idx++){
		if(_rEvent[idx].Event::Pass() && (_rEvent[idx].Event::W() >= _W_min_ && _rEvent[idx].Event::W() < _W_max_) && (_rEvent[idx].Event::Q2() >= _Q2_bins_[0] && _rEvent[idx].Event::Q2() < _Q2_bins_[5])){
			if(_rEvent[idx].Event::No_Nan()){
				_gevts+=1;
				_gevt_idx.push_back(idx);
				switch(_rEvent[idx].Event::Top()){
					case 0:
						_gevt_idx_mpro.push_back(_gEvent.size());
					break;
					case 1:
						_gevt_idx_mpip.push_back(_gEvent.size());
					break;
					case 2:
						_gevt_idx_mpim.push_back(_gEvent.size());
					break;
					case 3:
						_gevt_idx_mzero.push_back(_gEvent.size());
					break;
					default:
						std::cout<<"Invalid Top for Good Event\n";
					break;
				}
				_gEvent.push_back(_rEvent[idx]);
				_gtop[_rEvent[idx].Event::Top()]+=1;
			}
		}
	}
	//if(_rEvent.size()>0){
//		std::cout<<"_rEvent size :" <<_rEvent.size() <<"\n";
//
//		std::cout<<"_gEvent size :" <<_gEvent.size() <<"\n";
//	}
	
	if(_gEvent.size() != _gevts){
		std::cout<<"\nGood Event Size:" <<_gEvent.size() <<" _gevts:" <<_gevts <<"\n";
	}
	/*
	for(int j=0; j<4; j++){//Topology
		if(_ntop[j]>0){//Number of possible in that topology
			for(int k=0; k<_ntop[j]; k++){//Loop through the number in that topology
				//plot::plot_event(_rEvent[idx],hist_,flags_);
				if(_rEvent[idx].Event::Pass_Top(j)){
					_gEvent.push_back(_rEvent[idx]);
					_gevts+=1;
					_gevt_idx.push_back(idx);
					_gtop[j]+=1;
				}
				idx++;
			}
		}
	}*/
	//std::cout<<"\t\tNum good events: " <<_gtop[0] <<_gtop[1] <<_gtop[2] <<_gtop[3] <<"\n";
		
}

void Analysis::Isolate_Event(std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Isolated()){return;}
	if(_gevts>0){
		if(_gevts>1){
			//std::cout<<"Isolating Events: "<<_gtop[0] <<" " <<_gtop[1] <<" " <<_gtop[2] <<" " <<_gtop[3] <<"\n";
			//std::cout<<"Isolating Events with " <<_gevts <<" good events => " <<_gtop[0] <<_gtop[1] <<_gtop[2] <<_gtop[3]  <<"\n";
		}else{
			//std::cout<<"\tonly one good event\n";
		}
		if(_gtop[3]>0){
		//if(_gtop[3]>0 && _gtop[0]>0 && _gtop[1]>0 && _gtop[2]>0){	
			_top_passed = 3;
			hist_->Histogram::Top_Increment(3);
			for(int i=0; i<3; i++){
				if(_gtop[i]>0){
					hist_->Histogram::Top_Pot_Increment(i);
				}
			}
			if(_gtop[3]>1){
				Analysis::Isolate_Top(3,hist_,flags_);
			}else{
				_iEvent.push_back(_gEvent[Analysis::gEvent_idx(3,0)]);
				if(_gEvent[Analysis::gEvent_idx(3,0)].Event::Pass()){
					//hist_->Histogram::Isolated_Top_Increment(3);
					//std::cout<<"isolated to Zero Missing\n";
				}else{
					std::cout<<"We have a weird problem with the event passage\n";
				}
			}
		}else{
			if(_gtop[2]>0){
				//std::cout<<"Good Pim Missing\n";
				_top_passed = 2;
				hist_->Histogram::Top_Increment(2);
				for(int i=0; i<2; i++){
					if(_gtop[i]>0){
						hist_->Histogram::Top_Pot_Increment(i);
					}
				}
				if(_gtop[2]>1){
					Analysis::Isolate_Top(2,hist_,flags_);
				}else{
					_iEvent.push_back(_gEvent[Analysis::gEvent_idx(2,0)]);
					//hist_->Histogram::Isolated_Top_Increment(2);
					//std::cout<<"isolated to PIM Missing\n";
				}
				//std::cout<<"\tisolated event size:" <<_iEvent.size() <<"\n";
			}else{
				if(_gtop[1]>0){
					_top_passed = 1;
					hist_->Histogram::Top_Increment(1);
					if(_gtop[0]>0){
						hist_->Histogram::Top_Pot_Increment(0);
					}
					if(_gtop[1]>1){
						Analysis::Isolate_Top(1,hist_,flags_);
					}else{
						_iEvent.push_back(_gEvent[Analysis::gEvent_idx(1,0)]);
						//hist_->Histogram::Isolated_Top_Increment(1);
						//std::cout<<"isolated to PIP Missing\n";
					}
				}else{
					if(_gtop[0]>0){
						hist_->Histogram::Top_Increment(0);
						_top_passed = 0;
						if(_gtop[0]>1){
							Analysis::Isolate_Top(0,hist_,flags_);
						}else{
							_iEvent.push_back(_gEvent[Analysis::gEvent_idx(0,0)]);
							//hist_->Histogram::Isolated_Top_Increment(0);
							//std::cout<<"isolated to Proton Missing\n";
							
						}
					}
				}
			}
		}
	}
}

void Analysis::Isolate_Top(int top_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	if(!flags_->Flags::Plot_Isolated()){return;}
	//std::cout<<"Isolating Topology: "<<_top_[top_] <<" " <<_gtop[top_] <<"\n";
	int idx = -1;
	float curr_err = 100.0;
	if(_gtop[top_]>0){
		curr_err = _gEvent[Analysis::gEvent_idx(top_,0)].Event::Error();
		idx = Analysis::gEvent_idx(top_,0);
		for(int i=0; i<_gtop[top_]; i++){
			//std::cout<<"\t"<<i-1<<" event index: " <<Analysis::gEvent_idx(top_,i-1) <<" error: "<<_gEvent[Analysis::gEvent_idx(top_,i-1)].Event::Error() <<"\n";
			//std::cout<<"\t"<<i<<" event index: " <<Analysis::gEvent_idx(top_,i) <<" error: "<<_gEvent[Analysis::gEvent_idx(top_,i)].Event::Error() <<"\n";
			if(_gEvent[Analysis::gEvent_idx(top_,i)].Event::Error() < curr_err){
				idx = Analysis::gEvent_idx(top_,i);
				curr_err = _gEvent[Analysis::gEvent_idx(top_,i)].Event::Error();
				//std::cout<<"\t\tIndex: " <<idx <<"\n";
			}
			/*else if(i==_gtop[top_]-1){
				idx = Analysis::gEvent_idx(top_,i);
				//std::cout<<"\t\tIndex: " <<idx <<"\n";
			}
			if(_gEvent[Analysis::gEvent_idx(top_,i-1)].Event::Error() < _gEvent[Analysis::gEvent_idx(top_,i)].Event::Error()){
				idx = Analysis::gEvent_idx(top_,i-1);
				//std::cout<<"\t\tIndex: " <<idx <<"\n";
			}else if(i==_gtop[top_]-1){
				idx = Analysis::gEvent_idx(top_,i);
				//std::cout<<"\t\tIndex: " <<idx <<"\n";
			}*/
		}
		if(idx>=0){
			_iEvent.push_back(_gEvent[idx]);
			//std::cout<<"isolated to " <<_top_[top_] <<" Index " <<idx <<"\n";
			//hist_->Histogram::Isolated_Top_Increment(top_);
		}else{
			std::cout<<"Bad Index for Good Event\n";
		}
	}else{
		std::cout<<"Did not have any good events for given topology: " <<_top_[top_] <<"\n"; 
	}
}

int Analysis::Event_idx(int top_, int top_idx_){
	int idx = 0;
	for(int i=0; i<top_; i++){
		idx+=_ntop[i];
	}
	return idx+top_idx_;
}

int Analysis::gEvent_idx(int top_, int top_idx_){
	int idx = -1;
	switch(top_){
		case 0:
			idx = _gevt_idx_mpro[top_idx_];
		break;
		case 1:
			idx = _gevt_idx_mpip[top_idx_];
		break;
		case 2:
			idx = _gevt_idx_mpim[top_idx_];
		break;
		case 3:
			idx = _gevt_idx_mzero[top_idx_];
		break;
		default:
			std::cout<<"Invalid Top for Good Event Index\n";
		break;
	}
	return idx;
}

//Determine whether half wave plate is in or out
bool Analysis::Half_Wave(int run_num_, std::shared_ptr<Flags> flags_){
	bool plate_in = false;
	int above = 0;
	int below = 0;
	bool at = false; 
	if(flags_->Run()==_e16_){
		for(int i=0; i<sizeof(_plate_swap_e16_); i++){
			if(run_num_ > _plate_swap_e16_[i]){
				below +=1; 
			}else if(run_num_ == _plate_swap_e16_[i]){
				at = true; 
			}else if(run_num_ < _plate_swap_e16_[i]){
				above += 1; 
			}
		}
	}else if(flags_->Run()==_e1f_){
		for(int i=0; i<sizeof(_plate_swap_e1f_); i++){
			if(run_num_ > _plate_swap_e1f_[i]){
				below +=1; 
			}else if(run_num_ == _plate_swap_e1f_[i]){
				at = true; 
			}else if(run_num_ < _plate_swap_e1f_[i]){
				above += 1; 
			}
		}
	}
	if(at){
		if(below%2 == 0){
			plate_in = true;
		}else{
			plate_in = false;
		}
	}else{
		if(below%2 == 0){
			plate_in = false;
		}else{
			plate_in = true;
		}
	}
	return plate_in;
}


//Correct Helicity accoridng to half wave plate status
float Analysis::Corr_Helicity(float helicity_, int run_num_, std::shared_ptr<Flags> flags_){
	
	int eh = 0; 
	if(helicity_ >= 1000) eh = 1; 
	if(helicity_ <= -1000) eh = -1; 
	if(helicity_ < 1000 && helicity_ > -1000) eh = 0; 
	//if(plate_stat == 0 ) eh = 1; 
	if(Analysis::Half_Wave(run_num_,flags_)){
		return eh*-1;
	}else{
		return  eh;
	}
}

void Analysis::Plot_Particles(std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//std::cout<<"\t****Plotting Particles******\n";
	for(int i=0; i<_rParticle.size(); i++){
		if(_rParticle[i].Get_q()<0){
			if(_rParticle[i].Get_idx()==0){
				plot::plot_pid(_rParticle[i],hist_,flags_);
			}else if(_pim_idx.size()==1){
				plot::plot_pid(_rParticle[i],hist_,flags_);
			}
		}else if(_rParticle[i].Get_q()>0){
			if(_pro_idx.size()==1){
				plot::plot_pid(_rParticle[i],hist_,flags_);
			}
			if(_pip_idx.size()==1){
				plot::plot_pid(_rParticle[i],hist_,flags_);
			}
		}
		//plot::plot_pid(_rParticle[i],hist_,flags_);
	}
}

void Analysis::Plot_Events(std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//std::cout<<"\t\tPlotting Events\n";
	//if(_rEvent.size() >0){
//		std::cout<<"_rEvent size:" <<_rEvent.size() <<" _gEvent size:" <<_gEvent.size() <<" _iEvent size:" <<_iEvent.size() <<"\n";
//	}

	if(_rEvent.size() == 1){
		//if((_rEvent[0].Event::W() >= _W_min_ && _rEvent[0].Event::W() < _W_max_) && (_rEvent[0].Event::Q2() >= _Q2_bins_[0] && _rEvent[0].Event::Q2() < _Q2_bins_[5])){
		if(_rEvent[0].Event::No_Nan()){
			plot::plot_clean_event(_rEvent[0],hist_,flags_);
			plot::plot_event(_rEvent[0],hist_,flags_,false);
		}
		//}
		//std::cout<<"_rEvent W:" <<_rEvent[0].Event::W() <<" Q2:" <<_rEvent[0].Event::Q2() <<"\n";
	//if(_gEvent.size() == 1){
		//std::cout<<"\nPlot Clean event\n";
		//if((_rEvent[0].Event::W() >= _W_min_ && _rEvent[0].Event::W() < _W_max_) && (_rEvent[0].Event::Q2() >= _Q2_bins_[0] && _rEvent[0].Event::Q2() < _Q2_bins_[5])){
		//if((_gEvent[0].Event::W() >= _W_min_ && _gEvent[0].Event::W() < _W_max_) && (_gEvent[0].Event::Q2() >= _Q2_bins_[0] && _gEvent[0].Event::Q2() < _Q2_bins_[5])){
			//if(!_iEvent.size()>0){
			//	//hist_->Histogram::Clean_Not_Isolated_Top_Increment(_rEvent[0].Event::Top());
			//	hist_->Histogram::Clean_Not_Isolated_Top_Increment(_gEvent[0].Event::Top());
			//}else{
			//	//hist_->Histogram::Clean_and_Isolated_Top_Increment(_rEvent[0].Event::Top());
			//	hist_->Histogram::Clean_and_Isolated_Top_Increment(_gEvent[0].Event::Top());
			//}
		//}
		
		//plot::plot_clean_event(_gEvent[0],hist_,flags_);
		//hist_->Histogram::Clean_Top_Increment(_rEvent[0].Event::Top());
		//hist_->Histogram::Clean_Top_Increment(_gEvent[0].Event::Top());
		
		//plot::plot_event(_gEvent[0],hist_,flags_,false);
	}else if(_rEvent.size()>1){
	//}else if(_gEvent.size()>1){
		for(int i=0; i<_rEvent.size(); i++){
		//for(int i=0; i<_gEvent.size(); i++){
			//if((_rEvent[i].Event::W() >= _W_min_ && _rEvent[i].Event::W() < _W_max_) && (_rEvent[i].Event::Q2() >= _Q2_bins_[i] && _rEvent[0].Event::Q2() < _Q2_bins_[5])){
				//std::cout<<"_rEvent W:" <<_rEvent[i].Event::W() <<" Q2:" <<_rEvent[i].Event::Q2() <<"\n";
			//if(_rEvent[i].Event::No_Nan()){
			plot::plot_event(_rEvent[i],hist_,flags_,false);
			if(_top_[_rEvent[i].Event::Top()]==_mzero_ && _ntop[3]==1){// && _ntop[0]>0 && _ntop[1]>0 && _ntop[2]>0){
				plot::plot_clean_event(_rEvent[i],hist_,flags_);
			}else if(_top_[_rEvent[i].Event::Top()]==_mpim_ && _ntop[2]==1){
				plot::plot_clean_event(_rEvent[i],hist_,flags_);
			}else if(_top_[_rEvent[i].Event::Top()]==_mpip_ &&_ntop[1] == 1){
				plot::plot_clean_event(_rEvent[i],hist_,flags_);
			}else if(_top_[_rEvent[i].Event::Top()]==_mpro_ & _ntop[0] == 1){
				plot::plot_clean_event(_rEvent[i],hist_,flags_);
			}
			//}
			//}
			//std::cout<<"\nPlot dirty events\n";
			
			//plot::plot_event(_gEvent[i],hist_,flags_,false);
			
			//if(_top_[_gEvent[i].Event::Top()]==_mzero_){
				//if(_ntop[3]==1 && _ntop[0] && _ntop[1] && _ntop[2]){
					//if(!_iEvent.size()>0){
						//hist_->Histogram::Clean_Not_Isolated_Top_Increment(_rEvent[0].Event::Top());
					//	hist_->Histogram::Clean_Not_Isolated_Top_Increment(_gEvent[0].Event::Top());
					//}else{
						//hist_->Histogram::Clean_and_Isolated_Top_Increment(_rEvent[0].Event::Top());
					//	hist_->Histogram::Clean_and_Isolated_Top_Increment(_gEvent[0].Event::Top());
					//}
				
					//plot::plot_clean_event(_gEvent[i],hist_,flags_);
					//hist_->Histogram::Clean_Top_Increment(3);
				//}
			
		}
	}
	if(_gEvent.size() == 1){
		hist_->Histogram::Clean_Top_Increment(_gEvent[0].Event::Top());
		if(!_iEvent.size()>0){
			//hist_->Histogram::Clean_Not_Isolated_Top_Increment(_rEvent[0].Event::Top());
			hist_->Histogram::Clean_Not_Isolated_Top_Increment(_gEvent[0].Event::Top());
		}else{
			//hist_->Histogram::Clean_and_Isolated_Top_Increment(_rEvent[0].Event::Top());
			hist_->Histogram::Clean_and_Isolated_Top_Increment(_gEvent[0].Event::Top());
		}
	}else{
		bool mpro_pass = false;
		bool mpip_pass = false;
		bool mpim_pass = false;
		bool mzero_pass = false;
		int mpro_n = 0;
		int mpip_n = 0;
		int mpim_n = 0;
		int mzero_n = 0;
		for(int j=0; j<_gEvent.size(); j++){
			if(_gEvent[j].Event::Top()==0){
				mpro_pass = true;
				mpro_n += 1;
			}
			if(_gEvent[j].Event::Top()==1){
				mpip_pass = true;
				mpip_n += 1;
			}
			if(_gEvent[j].Event::Top()==2){
				mpim_pass = true;
				mpim_n += 1;
			}
			if(_gEvent[j].Event::Top()==3){
				mzero_pass = true;
				mzero_n += 1;
			}
			hist_->Histogram::Fill_Top_Check(mpro_pass,mpip_pass,mpim_pass,mzero_pass,mpro_n,mpip_n,mpim_n,mzero_n,flags_);
		}
		if(flags_->Flags::Plot_Isolated()){
			if(_gtop[3]==1){
				for(int i=0; i<_gEvent.size(); i++){
					if(_gEvent[i].Event::Top()==3 && _gtop[0]>0 && _gtop[1]>0 && _gtop[2]>0){
						hist_->Histogram::Clean_Top_Increment(_gEvent[i].Event::Top());
						if(!_iEvent.size()>0){
							//hist_->Histogram::Clean_Not_Isolated_Top_Increment(_rEvent[0].Event::Top());
							hist_->Histogram::Clean_Not_Isolated_Top_Increment(_gEvent[i].Event::Top());
						}else{
							//hist_->Histogram::Clean_and_Isolated_Top_Increment(_rEvent[0].Event::Top());
							hist_->Histogram::Clean_and_Isolated_Top_Increment(_gEvent[i].Event::Top());
						}
					}
				}
			}
		}else{
			bool any_difference = false;
			int num_forms = 1;
			std::vector<int> great_mixed_idx1;
			std::vector<int> great_mixed_idx2;
			std::vector<int> great_mixed_idx3;
			std::vector<std::vector<int>> good_mixed_friend_idx1;
			std::vector<std::vector<int>> good_mixed_friend_idx2;
			std::vector<std::vector<int>> good_mixed_friend_idx3;
			std::vector<std::vector<int>> great_mixed_friend_idx1;
			std::vector<std::vector<int>> great_mixed_friend_idx2;
			std::vector<std::vector<int>> great_mixed_friend_idx3;
			if(_gEvent.size()>=1){
				std::cout<<"---new event set----\n";
				for(int i=0; i< _gEvent.size(); i++){
					good_mixed_friend_idx1.push_back(hist_->Histogram::Friend_idx(_gEvent[i].Event::W(),_gEvent[i].Event::Q2(),_gEvent[i].Event::MMb(0),_gEvent[i].Event::MM2b(0),_gEvent[i].Event::Thetab(0),_gEvent[i].Event::Alphab(0),_gEvent[i].Event::Phib(0),0));
					good_mixed_friend_idx1.push_back(hist_->Histogram::Friend_idx(_gEvent[i].Event::W(),_gEvent[i].Event::Q2(),_gEvent[i].Event::MMb(1),_gEvent[i].Event::MM2b(1),_gEvent[i].Event::Thetab(1),_gEvent[i].Event::Alphab(1),_gEvent[i].Event::Phib(1),1));
					good_mixed_friend_idx1.push_back(hist_->Histogram::Friend_idx(_gEvent[i].Event::W(),_gEvent[i].Event::Q2(),_gEvent[i].Event::MMb(2),_gEvent[i].Event::MM2b(2),_gEvent[i].Event::Thetab(2),_gEvent[i].Event::Alphab(2),_gEvent[i].Event::Phib(2),2));
					std::cout<<"\tEvent: " <<i <<"\n";
					std::cout<<"\tTopology:" <<_top_[_gEvent[i].Event::Top()] <<"\n";
					std::cout<<"\t" <<hist_->Histogram::Friend_Event_Idx1(_gEvent[i].Event::W()) <<"\n";
					std::cout<<"\t" <<hist_->Histogram::Friend_Event_Idx2(_gEvent[i].Event::Q2()) <<"\n";
					std::cout<<"\t";
					for(int j=0; j<3; j++){
						std::cout<<hist_->Histogram::Friend_Event_Idx3(_gEvent[i].Event::MMb(j),_gEvent[i].Event::W(),j) <<" ";
						if(hist_->Histogram::Friend_Event_Idx3(_gEvent[i].Event::MMb(j),_gEvent[i].Event::W(),j)==-1){
							std::cout<<" MM:" <<_gEvent[i].Event::MMb(j) <<" ";
						}
					}
					std::cout<<"\n\t";
					for(int j=0; j<3; j++){
						std::cout<<hist_->Histogram::Friend_Event_Idx4(_gEvent[i].Event::MM2b(j),_gEvent[i].Event::W(),j)<<" ";
						if(hist_->Histogram::Friend_Event_Idx4(_gEvent[i].Event::MM2b(j),_gEvent[i].Event::W(),j)==-1){
							std::cout<<" MM2:" <<_gEvent[i].Event::MM2b(j) <<" ";
						}
					}
					std::cout<<"\n\t";
					for(int j=0; j<3; j++){
						std::cout<<hist_->Histogram::Friend_Event_Idx5(_gEvent[i].Event::Thetab(j)) <<" ";
					}
					std::cout<<"\n\t";
					for(int j=0; j<3; j++){
						std::cout<<hist_->Histogram::Friend_Event_Idx6(_gEvent[i].Event::Alphab(j))<<" ";
					}
					std::cout<<"\n\t";
					for(int j=0; j<3; j++){
						std::cout<<hist_->Histogram::Friend_Event_Idx7(_gEvent[i].Event::Phib(j)) <<" ";
					}
					std::cout<<"\n\t" <<_gEvent[i].Event::Weight() <<"\n";
				}
				bool boppin1 = true;
				bool boppin2 = true;
				for(int i=0; i< _gEvent.size(); i++){
					for(int j=0; j< _gEvent.size(); j++){
						if(i < j){
							if( good_mixed_friend_idx1[i]==good_mixed_friend_idx1[j]){
								if(great_mixed_friend_idx1.size()==0){
									good_mixed_friend_idx1.push_back(good_mixed_friend_idx1[i]);
									great_mixed_idx1.push_back(i);
								}
							}else{
								if(great_mixed_friend_idx1.size()==0){
									good_mixed_friend_idx1.push_back(good_mixed_friend_idx1[i]);
									great_mixed_idx1.push_back(i);
								}else{
									boppin1 = true;
									boppin2 = true;
									for(int k=0; k<great_mixed_friend_idx1.size(); k++){
										if(great_mixed_friend_idx1[k] == good_mixed_friend_idx1[i]){
											boppin1 = false;
										}
										if(great_mixed_friend_idx1[k] == good_mixed_friend_idx1[j]){
											boppin2 = false;
										}
									}
									if(boppin1){
										great_mixed_friend_idx1.push_back(good_mixed_friend_idx1[i]);
										great_mixed_idx1.push_back(i);
									}
									if(boppin2){
										great_mixed_friend_idx1.push_back(good_mixed_friend_idx1[j]);
										great_mixed_idx1.push_back(j);
									}
								}
							}
							if( good_mixed_friend_idx2[i]==good_mixed_friend_idx2[j]){
								if(great_mixed_friend_idx2.size()==0){
									good_mixed_friend_idx2.push_back(good_mixed_friend_idx2[i]);
									great_mixed_idx2.push_back(i);
								}
							}else{
								if(great_mixed_friend_idx2.size()==0){
									good_mixed_friend_idx2.push_back(good_mixed_friend_idx2[i]);
									great_mixed_idx2.push_back(i);
								}else{
									boppin1 = true;
									boppin2 = true;
									for(int k=0; k<great_mixed_friend_idx2.size(); k++){
										if(great_mixed_friend_idx2[k] == good_mixed_friend_idx2[i]){
											boppin1 = false;
										}
										if(great_mixed_friend_idx2[k] == good_mixed_friend_idx2[j]){
											boppin2 = false;
										}
									}
									if(boppin1){
										great_mixed_friend_idx2.push_back(good_mixed_friend_idx2[i]);
										great_mixed_idx2.push_back(i);
									}
									if(boppin2){
										great_mixed_friend_idx2.push_back(good_mixed_friend_idx2[j]);
										great_mixed_idx2.push_back(j);
									}
								}
							}
							if( good_mixed_friend_idx3[i]==good_mixed_friend_idx3[j]){
								if(great_mixed_friend_idx3.size()==0){
									good_mixed_friend_idx3.push_back(good_mixed_friend_idx3[i]);
									great_mixed_idx3.push_back(i);
								}
							}else{
								if(great_mixed_friend_idx3.size()==0){
									good_mixed_friend_idx3.push_back(good_mixed_friend_idx3[i]);
									great_mixed_idx3.push_back(i);
								}else{
									boppin1 = true;
									boppin2 = true;
									for(int k=0; k<great_mixed_friend_idx3.size(); k++){
										if(great_mixed_friend_idx3[k] == good_mixed_friend_idx3[i]){
											boppin1 = false;
										}
										if(great_mixed_friend_idx3[k] == good_mixed_friend_idx3[j]){
											boppin2 = false;
										}
									}
									if(boppin1){
										great_mixed_friend_idx3.push_back(good_mixed_friend_idx3[i]);
										great_mixed_idx3.push_back(i);
									}
									if(boppin2){
										great_mixed_friend_idx3.push_back(good_mixed_friend_idx3[j]);
										great_mixed_idx3.push_back(j);
									}
								}
							}
						}
					}
				}
			}
			if(_gEvent.size()>1){
				for(int i=0; i<_gEvent.size()-1; i++){
					
					for(int j=0; j<3; j++){
						if(hist_->Histogram::Friend_Event_Idx1(_gEvent[i].Event::W()) != hist_->Histogram::Friend_Event_Idx1(_gEvent[i+1].Event::W())){
							any_difference = true;
							std::cout<<"\tdiff in W| " <<_gEvent[i].Event::W() <<" vs. " <<_gEvent[i+1].Event::W() <<"  i:" <<i <<" j:" <<i+1 <<"\n";
						}
						if(hist_->Histogram::Friend_Event_Idx2(_gEvent[i].Event::Q2()) != hist_->Histogram::Friend_Event_Idx2(_gEvent[i+1].Event::Q2())){
							any_difference = true;
							std::cout<<"\tdiff in Q2| " <<_gEvent[i].Event::Q2() <<" vs. " <<_gEvent[i+1].Event::Q2() <<"  i:" <<i <<" j:" <<i+1 <<"\n";
						}
						if(hist_->Histogram::Friend_Event_Idx3(_gEvent[i].Event::MMb(j),_gEvent[i].Event::W(),j) != hist_->Histogram::Friend_Event_Idx3(_gEvent[i+1].Event::MMb(j),_gEvent[i+1].Event::W(),j)){
							any_difference = true;
							std::cout<<"\tdiff in MM1| " <<_gEvent[i].Event::MMb(j) <<" vs. " <<_gEvent[i+1].Event::MMb(j) <<"  i:" <<i <<" j:" <<i+1 <<"\n";
						}
						if(hist_->Histogram::Friend_Event_Idx4(_gEvent[i].Event::MM2b(j),_gEvent[i].Event::W(),j) != hist_->Histogram::Friend_Event_Idx4(_gEvent[i+1].Event::MM2b(j),_gEvent[i+1].Event::W(),j)){
							any_difference = true;
							std::cout<<"\tdiff in MM12| " <<_gEvent[i].Event::MM2b(j) <<" vs. " <<_gEvent[i+1].Event::MM2b(j) <<"  i:" <<i <<" j:" <<i+1 <<"\n";
						}
						if(hist_->Histogram::Friend_Event_Idx5(_gEvent[i].Event::Thetab(j)) != hist_->Histogram::Friend_Event_Idx5(_gEvent[i+1].Event::Thetab(j))){
							any_difference = true;
							std::cout<<"\tdiff in theta| " <<_gEvent[i].Event::Thetab(j) <<" vs. " <<_gEvent[i+1].Event::Thetab(j) <<"  i:" <<i <<" j:" <<i+1 <<"\n";
						}
						if(hist_->Histogram::Friend_Event_Idx6(_gEvent[i].Event::Alphab(j)) != hist_->Histogram::Friend_Event_Idx6(_gEvent[i+1].Event::Alphab(j))){
							any_difference = true;
							std::cout<<"\tdiff in alpha| " <<_gEvent[i].Event::Alphab(j) <<" vs. " <<_gEvent[i+1].Event::Alphab(j) <<"  i:" <<i <<" j:" <<i+1 <<"\n";
						}	
						if(hist_->Histogram::Friend_Event_Idx7(_gEvent[i].Event::Phib(j)) != hist_->Histogram::Friend_Event_Idx7(_gEvent[i+1].Event::Phib(j))){
							any_difference = true;
							std::cout<<"\tdiff in phi| " <<_gEvent[i].Event::Phib(j) <<" vs. " <<_gEvent[i+1].Event::Phib(j) <<"  i:" <<i <<" j:" <<i+1 <<"\n";
						}
						if(_gEvent[i].Event::Weight() != _gEvent[i+1].Event::Weight()){
							any_difference = true;
							std::cout<<"\tdiff in weight| " <<_gEvent[i].Event::Weight() <<" vs. " <<_gEvent[i+1].Event::Weight() <<"  i:" <<i <<" j:" <<i+1 <<"\n";
						}
					}
				}
				if(any_difference && great_mixed_friend_idx1.size()==0){
					std::cout<<"A difference, but not reflected in idx1!\n";
				}
				for(int i=0; i<_gEvent.size(); i++){
					for(int j=0; j<4; j++){
						if(_gtop[i] == 1 && _gEvent[i].Event::Top()==i){
							hist_->Histogram::Clean_Not_Isolated_Top_Increment(_gEvent[i].Event::Top());
						}
					}
					if(any_difference){
						fun::print_vector_idx(great_mixed_idx1);
						fun::print_vector_idx(great_mixed_idx2);
						fun::print_vector_idx(great_mixed_idx3);
						for(int j=0; j<great_mixed_friend_idx1.size(); j++){
							plot::plot_mixed_events(_gEvent[great_mixed_idx1[j]],hist_,flags_, _gEvent[great_mixed_idx1[j]].Event::Top(),great_mixed_friend_idx1.size());
						}
						//plot::plot_mixed_events(_gEvent[i],hist_,flags_, _gEvent[i].Event::Top(),_gEvent.size());
					}else if(i==0){
						plot::plot_mixed_events(_gEvent[i],hist_,flags_, _gEvent[i].Event::Top(),1);
					}
				}
			}else if(_gEvent.size()==1){
				hist_->Histogram::Clean_Not_Isolated_Top_Increment(_gEvent[0].Event::Top());
			}
		}
	}

	
	if(flags_->Flags::Sim() && _tEvent[0].Event::W()>= _W_min_ && _tEvent[0].Event::W() < _W_max_ && _tEvent[0].Event::Q2()>= _Q2_min_ && _tEvent[0].Event::Q2() < _Q2_max_){
		//std::cout<<"Thrown W:" <<_tEvent[0].Event::W() <<" Q2:" <<_tEvent[0].Event::Q2() <<"\n";
		plot::plot_event(_tEvent[0],hist_,flags_,true);
	}
	//if(_gEvent.size()>0){
//		std::cout<<" _iEvent.size()= " <<_iEvent.size() <<" _top_passed:"<<_top_passed <<"\n";
//	}
	if(_iEvent.size()>0){
		
		if(_iEvent.size()==1){
			hist_->Histogram::Isolated_Top_Increment(_top_passed);
			plot::plot_isolated_event(_iEvent[0],hist_,flags_, _top_passed);
		}else{
			std::cout<<"size: " <<_iEvent.size() <<"  Didn't isolate event: there's more than one here!\n";
		}
	}
}

