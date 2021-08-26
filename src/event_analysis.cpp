#include "event_analysis.hpp"
		  
//Analysis::Analysis(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Environment> envi_, int run_type_,std::shared_ptr<Forest> a_Forest_, int thread_id_, int run_num_){
Analysis::Analysis(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Forest> forest_, int thread_id_, int run_num_, std::shared_ptr<Flags> flags_){
	if(thread_id_==0){
		//std::cout<<"Starting Event Analysis\n";
	}
	//int npart = 0; //# of Detected Particles
	if(flags_->Flags::Sim()){
		_npart = data_->Branches::npart();
		_hel = 1; 
		_weight = data_->Branches::weight();
		//Thrown Particle
		for(int i=0; i<4; i++){
			_tParticle.push_back(Particle(i,data_,flags_,true));
		}
		//Thrown Event
		_tEvent.push_back(Event(3,_tParticle[0],_tParticle[1],_tParticle[2],_tParticle[3],flags_,_weight,1,true));
		plot::plot_event(_tEvent[0],hist_,flags_,true);
	}else{
		_npart = data_->Branches::gpart();
		_weight = 1.0;
		if(flags_->Flags::Helicity()){
			_hel = Analysis::Corr_Helicity(data_->Branches::evntclas2(),run_num_,flags_);//Need to pair with library of files where half-wave plate is in or out
		}
	}
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
		//Analysis::Isolate_Event(hist_, flags_);//Isolate Definitive Event and Plot Accordingly
		Analysis::Plot_Events(hist_,flags_);
	}
}

void Analysis::Num_top(){
	//Pro Missing
	if(_rParticle.size()>0){
		if(_rParticle[0].Particle::Is_Elec()){
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
			//Pip Missing
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
			//Pim Missing
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
			//Zero Missing
			if(_pro_idx.size() > 0 && _pip_idx.size() > 0 && _pim_idx.size() > 0){
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
			}
			//std::cout<<"\tNum Possible Topoogies: " <<_ntop[0] <<_ntop[1] <<_ntop[2] <<_ntop[3] <<"\n";
		}
	}
}

void Analysis::PID(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_,std::shared_ptr<Flags> flags_){
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
}

void Analysis::Event_ID(std::shared_ptr<Branches> data_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	Analysis::Num_top();
	for(int i=0; i<4; i++){
		if(_ntop[i]>0){
			for(int j=0; j<_ntop[i]; j++){
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
	int idx = 0; 
	for(int j=0; j<4; j++){
		if(_ntop[j]>0){
			for(int k=0; k<_ntop[j]; k++){
				if(_rEvent[idx].Event::Pass()){
					_gEvent.push_back(_rEvent[idx]);
					plot::plot_event(_rEvent[idx],hist_,flags_);
					_gevts+=1;
					_gevt_idx.push_back(idx);
					_gtop[j]+=1;
				}
				idx++;
			}
		}
	}
	//std::cout<<"\t\tNum good events: " <<_gtop[0] <<_gtop[1] <<_gtop[2] <<_gtop[3] <<"\n";
		
}

void Analysis::Isolate_Event(std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	if(_gevts>0){
		//std::cout<<"Isolating Events with " <<_gevts <<" good events => " <<_gtop[0] <<_gtop[1] <<_gtop[2] <<_gtop[3]  <<"\n";
		if(_gtop[3]>0){
			if(_gtop[3]>1){
				Analysis::Isolate_Top(3,hist_,flags_);
			}else{
				_iEvent.push_back(_rEvent[Analysis::Event_idx(3,0)]);
			}
		}else{
			if(_gtop[2]>0){
				if(_gtop[2]>1){
					Analysis::Isolate_Top(2,hist_,flags_);
				}else{
					_iEvent.push_back(_rEvent[Analysis::Event_idx(2,0)]);
				}
			}else{
				if(_gtop[1]>0){
					if(_gtop[2]>1){
						Analysis::Isolate_Top(1,hist_,flags_);
					}else{
						_iEvent.push_back(_rEvent[Analysis::Event_idx(1,0)]);
					}
				}else{
					if(_gtop[0]>0){
						if(_gtop[2]>1){
							Analysis::Isolate_Top(0,hist_,flags_);
						}else{
							_iEvent.push_back(_rEvent[Analysis::Event_idx(0,0)]);
						}
					}
				}
			}
		}
		if(_iEvent.size()==1){
			plot::plot_isolated_event(_iEvent[0],hist_,flags_);
		}else{
			std::cout<<"size: " <<_iEvent.size() <<"  Didn't isolate event: there's more than one here!\n";
		}
	}
}

void Analysis::Isolate_Top(int top_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	int idx = -1;
	if(_gtop[top_]>0){
		for(int i=1; i<_gtop[top_]; i++){
			if(_gEvent[Analysis::Event_idx(top_,i-1)].Event::Error() < _gEvent[Analysis::Event_idx(top_,i)].Event::Error()){
				idx = Analysis::Event_idx(top_,i-1);
			}
		}
		_iEvent.push_back(_gEvent[idx]);
	}
	
}

int Analysis::Event_idx(int top_, int top_idx_){
	int idx = 0;
	for(int i=0; i<top_; i++){
		idx+=_ntop[i];
	}
	return idx+top_idx_;
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
	if(Analysis::Half_Wave(run_num_,flags_)){
		return helicity_*-1;
	}else{
		return  helicity_;
	}
}

void Analysis::Plot_Particles(std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//std::cout<<"\tPlotting Particles\n";
	for(int i=0; i<_rParticle.size(); i++){
		plot::plot_pid(_rParticle[i],hist_,flags_);
	}
}

void Analysis::Plot_Events(std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_){
	//std::cout<<"\t\tPlotting Events\n";
	if(_rEvent.size() == 1){
		//std::cout<<"\nPlot Clean event\n";
		plot::plot_clean_event(_rEvent[0],hist_,flags_);
	}else{
		for(int i=0; i<_rEvent.size(); i++){
			//std::cout<<"\nPlot dirty events\n";
			plot::plot_event(_rEvent[i],hist_,flags_);
		}
	}
	if(flags_->Flags::Sim()){
		plot::plot_event(_tEvent[0],hist_,flags_,true);
	}
}

