/* Invariant mass cut classes
 *
 * Provides a set of cuts related to invariant mass by overloading the 
 * cut base class.
*/

#ifndef INC_CUTS_MASS
#define INC_CUTS_MASS

#include <iostream>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "../event/event.h"
#include "../cuts/cuts.h"


/* NAMESPACE */
namespace analysis
{

	/* invariant jet combination mass cut class */
	class cut_gres : public cut
	{
	
	public:
		cut_gres();

		void init();
		
		bool operator() (const event *ev);
		
		bool passed_top21(const event *ev);
		bool passed_top31(const event *ev);
		bool passed_top22(const event *ev);

		std::string name() const;

		void set_topology(const std::vector<int> & top);
		void set_top_cuts(const std::vector<double> & cuts);

		void calc_nrjets();
		void remove_singleres();
		void calc_nrcombinations();
		void calc_combinations();
		void calc_pairings();

	private:
		bool is_valid_pairing(const std::vector<std::vector<unsigned int> > &pairing) const;

	private:
		std::vector<int> topology;
		std::vector<double> top_cuts;

		unsigned int nr_jets;
		unsigned int nr_combinations;
		std::vector<int> reduced_topology;
		std::vector<double> reduced_cuts;
		std::vector<std::vector<std::vector<unsigned int> > > resonance_combinations;
		std::vector<std::vector<std::vector<unsigned int> > > resonance_pairings;

	};

/* NAMESPACE */
}

#endif
