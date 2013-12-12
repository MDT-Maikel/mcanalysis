/* Invariant mass cut classes
 *
 * Provides a set of cuts related to invariant mass by overloading the 
 * cut base class.
*/

#ifndef INC_CUTS_ADVANCED
#define INC_CUTS_ADVANCED

#include <iostream>
#include <string>
#include <vector>

#include "boost/lexical_cast.hpp"

#include "../event/event.h"
#include "cuts.h"


/* NAMESPACE */
namespace analysis
{

	/* invariant mass cut class */
	class cut_mass : public cut
	{
	
	public:
		cut_mass(double minm = 0, double maxm = 0);
	
		virtual void set_min(double min);
		virtual void set_max(double max);

	private:
		double min_mass;
		double max_mass;

	};

/* NAMESPACE */
}

#endif
