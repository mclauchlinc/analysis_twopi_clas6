#ifndef PLOT_HPP
#define PLOT_HPP

#include "histogram.hpp"
#include "particle.hpp"
#include "event.hpp"
#include "flags.hpp"
#include "forest.hpp"

namespace plot{
	//Plotting Particle ID
	void plot_pid(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	void plot_ele(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	void plot_pro(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	void plot_pip(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	void plot_pim(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	//Individual Plotting
	void plot_thrown(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	void plot_no_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_);
	void plot_sanity_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_);
	void plot_vertex_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_);
	void plot_fid_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_);
	void plot_delta_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_);
	void plot_sf_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_);
	void plot_cc_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_);
	void plot_ec_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_);
	//void plot_beta_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_);
	void plot_pid_cut(Particle particle_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, int par_);
	//Event Plotting
	void plot_event(Event event_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_, bool thrown_ = false);
	void plot_clean_event(Event event_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
	void plot_isolated_event(Event event_, std::shared_ptr<Histogram> hist_, std::shared_ptr<Flags> flags_);
}

#endif