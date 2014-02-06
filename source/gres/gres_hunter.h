/* GRESHunter class
 *
 * 
*/

#ifndef INC_GRESHUNTER
#define INC_GRESHUNTER

#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include <TH1.h> 

#include "../bumphunter/bumphunter.h"
#include "../event/event.h"


/* NAMESPACE */
namespace analysis
{

	class greshunter
	{

	public:
		
		/* con & destructor */
		greshunter() = default;
		
		/* copy & assignment */
		
		/* properties */
		void add_background(const std::vector<event*> &events, double weight);
		void set_signal(const std::vector<event*> &events, double weight);
		
		void set_topology(const std::vector<int> & top);
		void set_top_cuts(const std::vector<double> & cuts);
		
		void set_folder(std::string f);
		
		/* bump hunting */
		void run();
		void run(const std::vector<int> & comb, double mass_cut);

	private:
	
		std::vector< std::vector<event*> > events_bkg;
		std::vector<double> weight_bkg;
	
		std::vector<event*> events_sig;
		double weight_sig;
		
		std::string gres_folder;
		
		std::vector<int> topology;
		std::vector<double> top_cuts;
		
	};

/* NAMESPACE */
}

#endif
