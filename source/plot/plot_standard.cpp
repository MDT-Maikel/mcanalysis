/* Standard plot classes
 *
 *
*/

#include "plot_standard.h"


/* NAMESPACE */
namespace analysis
{

	/* plot pt class */

	plot_pt::plot_pt(int t, unsigned int n) : plot()
	{
		type = t;
		number = n;
		hist.set_x_label("pt");
		hist.set_y_label("# events");
		hist.set_leg_title("Legend");
		hist.set_title("plot_" + name() + ".ps");
	}

	void plot_pt::add_sample(const std::vector<event*> &events, const std::string &name)
	{
		std::vector<double> pt;
		for (unsigned int j = 0; j < events.size(); j++)
		{
			event *ev = events[j];				
			particle *p = ev->get(type, number);
			if (p)
				pt.push_back(p->pt());
		}
		hist.add_sample(pt, name);
	}	

	std::string plot_pt::name() const
	{
		return "pt particle " + boost::lexical_cast<std::string>(number) + "";
	}

	/* plot met class */

	plot_met::plot_met() : plot()
	{
		hist.set_x_label("met");
		hist.set_y_label("# events");
		hist.set_leg_title("Legend");
		hist.set_title("plot_" + name() + ".ps");
	}

	void plot_met::add_sample(const std::vector<event*> &events, const std::string &name)
	{
		std::vector<double> met;
		for (unsigned int j = 0; j < events.size(); j++)
		{
			event *ev = events[j];				
			met.push_back(ev->met());
		}
		hist.add_sample(met, name);
	}	

	std::string plot_met::name() const
	{
		return "met";
	}

	/* plot ht class */

	plot_ht::plot_ht() : plot()
	{
		hist.set_x_label("ht");
		hist.set_y_label("# events");
		hist.set_leg_title("Legend");
		hist.set_title("plot_" + name() + ".ps");
	}

	void plot_ht::add_sample(const std::vector<event*> &events, const std::string &name)
	{
		std::vector<double> ht;
		for (unsigned int j = 0; j < events.size(); j++)
		{
			event *ev = events[j];				
			ht.push_back(ev->ht(particle::type_jet, 40, 2.8));
		}
		hist.add_sample(ht, name);
	}	

	std::string plot_ht::name() const
	{
		return "ht";
	}

/* NAMESPACE */
}
