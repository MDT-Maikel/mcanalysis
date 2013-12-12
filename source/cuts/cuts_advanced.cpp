/* Invariant mass cut classes
 *
 * Provides a set of cuts related to invariant mass by overloading the 
 * cut base class.
*/

#include "cuts_advanced.h"


/* NAMESPACE */
namespace analysis
{

	/* invariant mass cut class */

	cut_mass::cut_mass(double minm, double maxm)
	{
		min_mass = minm;
		max_mass = maxm;
	}

	void cut_mass::set_min(double min)
	{
		min_mass = min;
	}
	
	void cut_mass::set_max(double max)
	{
		max_mass = max;
	}

/* NAMESPACE */
}
