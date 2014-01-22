/* 2D Plot class
 *
 *
*/

#include "plot2d.h"


/* NAMESPACE */
namespace analysis
{

	/* plot base class */

	plot2d::plot2d(std::string new_name)
	{
		// name defaults to empty string
		plot_name = new_name;
		plot_folder = "";

		// standard histogram options		
		hist.set_title(plot_folder + "plot_" + name());
		hist.set_palette();
		hist.set_x_bins(100);
		hist.set_y_bins(100);
	}

	plot2d::~plot2d()
	{
	}

	void plot2d::add_sample(const std::vector<event*> &events, const std::string &name, double weight)
	{
		// implementation: add data to histogram and set sample name
	}

	void plot2d::run()
	{
		// set histogram title based on name
		hist.set_title(plot_folder + "plot_" + name());
		// draw histograms
		hist.draw();
	}

	/* plot properties */

	std::string plot2d::name() const
	{
		if (plot_name != "")
			return plot_name;		
		return "WARNING: using base plot class";
	}

	void plot2d::set_name(std::string n)
	{
		plot_name = n;
	}

	void plot2d::set_folder(std::string f)
	{
		plot_folder = f;
	}

	void plot2d::set_x_bins(double nbins, double min, double max)
	{
		hist.set_x_bins(nbins);
		hist.set_x_range(min, max);		
	}
	
	void plot2d::set_y_bins(double nbins, double min, double max)
	{
		hist.set_y_bins(nbins);
		hist.set_y_range(min, max);		
	}
	
	/* plot met class */

	plot2d_ht_met::plot2d_ht_met() : plot2d()
	{
		hist.set_x_label("ht");
		hist.set_y_label("met");
		hist.set_title("plot_" + name());
	}

	void plot2d_ht_met::add_sample(const std::vector<event*> &events, const std::string &name, double weight)
	{
		std::vector< std::vector<double> > ht_met;
		for (unsigned int j = 0; j < events.size(); j++)
		{
			event *ev = events[j];				
			ht_met.push_back({ev->ht(particle::type_jet, 40, 2.8), ev->met()});
		}
		hist.add_sample_xy(ht_met);
	}	

	std::string plot2d_ht_met::name() const
	{
		if (plot_name != "")
			return plot_name;		
		return "ht_met";
	}

/* NAMESPACE */
}
