/* 2D Plot class
 *
 *
*/

#include "plot2d.h"


/* NAMESPACE */
namespace analysis
{

	/* plot base class */

	plot2d::plot2d(std::string name, std::string folder)
	{
		// name defaults to empty string
		plot_name = name;
		plot_folder = folder;

		// standard histogram options		
		hist.set_title(plot_folder + plot_name);
		hist.set_palette();
		hist.set_x_bins(50);
		hist.set_y_bins(50);
	}

	plot2d::~plot2d()
	{
	}

	/* plot data */

	void plot2d::add_sample(const std::vector<event*> &events, const std::string &name, double weight)
	{
		std::vector< std::vector<double> > result;
		for (unsigned int j = 0; j < events.size(); j++)
		{
			event *ev = events[j];
			// use ftor's to get the plot data
			if (!!x_ftor && !!y_ftor)
				result.push_back({(*x_ftor)(ev), (*y_ftor)(ev)});
		}
		hist.add_sample_xy(result);	
	}

	void plot2d::run()
	{
		// set histogram title based on name
		hist.set_title(plot_folder + plot_name);
		// draw histograms
		hist.draw();
	}

	/* plot properties */

	void plot2d::set_name(std::string n)
	{
		plot_name = n;
	}

	void plot2d::set_folder(std::string f)
	{
		plot_folder = f;
	}
	
	void plot2d::set_functors(plot_default *x, plot_default *y)
	{
		x_ftor = x;
		y_ftor = y;		
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
	
/* NAMESPACE */
}
