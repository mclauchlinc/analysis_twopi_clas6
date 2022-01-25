#include "cuts.hpp"

bool cuts::W_cut(float W_, int run_){
	if(W_ > _W_cut_param_[0] && W_ < _W_cut_param_[1]){
		return true;
	}
	return false;
}