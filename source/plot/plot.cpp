/* Plot class
 *
 *
*/

#include "plot.h"


/* NAMESPACE */
namespace analysis
{

	/* plot base class */

	plot::plot(std::string new_name)
	{
		// standard histogram options		
		hist.set_title(name() + ".png");
		hist.set_leg_title("Legend");
		hist.set_bins(100);
		hist.set_logy(true);
		// name defaults to empty string
		plot_name = new_name;
	}

	plot::~plot()
	{
	}

	void plot::add_sample(const std::vector<event*> &events, const std::string &name, double weight)
	{
		// implementation: add data to histogram and set sample name
	}

	void plot::run()
	{
		// set histogram title based on name
		hist.set_title("plot_" + name());
		// draw histograms
		hist.draw();
	}

	/* plot properties */

	std::string plot::name() const
	{
		if (plot_name != "")
			return plot_name;		
		return "WARNING: using base plot class";
	}

	void plot::set_name(std::string n)
	{
		plot_name = n;
	}

	void plot::set_logy(bool on)
	{
		hist.set_logy(on);
	}

	void plot::set_stacked(bool on)
	{
		hist.set_stacked(on);
	}

	void plot::set_normalized(bool on)
	{
		hist.set_normalized(on);
	}

/* NAMESPACE */
}
