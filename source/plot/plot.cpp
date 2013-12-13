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
		hist.set_title(name() + ".ps");
		hist.set_leg_title("Legend");		
	}

	plot::~plot()
	{
	}

	void plot::run(std::vector<event*> &events)
	{

	}		

	std::string plot::name() const
	{
		return "WARNING: using base plot class";
	}

	/* plot pt class */

	plot_pt::plot_pt(int t, unsigned int n) : plot()
	{
		type = t;
		number = n;
		hist.set_x_label("pt");
		hist.set_y_label("# events");
		hist.set_leg_title("Legend");
		hist.set_title(name() + ".ps");
	}

	void plot_pt::run(std::vector<event*> &events)
	{
		hist.set_bins(100);
		//hist.set_range(0, 2000);		

		// construct the pt data from the events
		for (unsigned int i = 1; i <= number; i++)
		{
			std::vector<double> pt;
			for (unsigned int j = 0; j < events.size(); j++)
			{
				event *ev = events[j];				
				particle *p = ev->get(type, i);
				if (p)
					pt.push_back(p->pt());
			}
			hist.add_sample(pt, "particle " + boost::lexical_cast<std::string>(i));
		}

		hist.normalize();
		hist.draw();
	}	

	std::string plot_pt::name() const
	{
		return "pt_plot";
	}

/* NAMESPACE */
}
