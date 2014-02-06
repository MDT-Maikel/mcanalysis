/* Plot class
 *
 *
*/

#include "plot.h"


/* NAMESPACE */
namespace analysis
{

	/* plotting class */

	plot::plot(std::string name, std::string folder)
	{
		// name defaults to empty string
		plot_name = name;
		plot_folder = folder;

		// standard histogram options		
		hist.set_title(plot_folder + plot_name);
		hist.set_leg_title("Legend");
		hist.set_bins(100);
		hist.set_logy(true);	
	}

	plot::~plot()
	{
	}
	
	/* plot data */
	
	void plot::add_sample(const std::vector<event*> &events, plot_implementation *plot_imp, const std::string &name, double weight)
	{
		std::vector<double> result;
		for (unsigned int j = 0; j < events.size(); j++)
		{
			double res = (*plot_imp)(events[j]);
			result.push_back(res);
		}
		hist.add_sample(result, name, weight);	
	}

	void plot::add_sample(const std::vector<event*> &events, double(*func)(const event*), const std::string &name, double weight)
	{
		std::vector<double> result;
		for (unsigned int j = 0; j < events.size(); j++)
		{
			double res = func(events[j]);
			result.push_back(res);
		}
		hist.add_sample(result, name, weight);	
	}

	void plot::run()
	{
		// set histogram title based on name
		hist.set_title(plot_folder + plot_name);
		// draw histograms
		hist.draw();
	}

	/* plot properties */

	void plot::set_name(std::string n)
	{
		plot_name = n;
	}

	void plot::set_folder(std::string f)
	{
		plot_folder = f;
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
	
	void plot::set_bins(double nbins, double min, double max)
	{
		hist.set_bins(nbins);
		hist.set_range(min, max);		
	}
	
	/* standard plot definitions */
	
	double plot_pt(const event *ev)
	{
		if (ev->size() == 0)
			return 0;
		return (*ev)[0]->pt();		
	}

	double plot_met(const event *ev)
	{
		return ev->met();	
	}
	
	double plot_ht(const event *ev)
	{
		return ev->ht(ptype_all, 40, 2.8);	
	}

	double plot_mass(const event *ev)
	{
		return ev->mass();
	}
	
	double plot_mt2(const event *ev)
	{
		return ev->mt2();	
	}

/* NAMESPACE */
}
