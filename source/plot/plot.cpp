/* Plot class
 *
 *
*/

#include "plot.h"


/* NAMESPACE */
namespace analysis
{

	/* plot base class */

	plot::plot()
	{
		// standard histogram options		
		hist.set_title(name() + ".ps");
		hist.set_leg_title("Legend");
		hist.set_bins(100);
		hist.set_logy(true);	
	}

	plot::~plot()
	{
	}

	void plot::add_sample(const std::vector<event*> &events, const std::string &name)
	{
		// implementation: add data to histogram and set sample name
	}

	void plot::run()
	{
		// draw histogram
		hist.normalize();
		hist.draw();
	}	

	std::string plot::name() const
	{
		return "WARNING: using base plot class";
	}

/* NAMESPACE */
}
