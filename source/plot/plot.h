/* Plot class
 *
 *
*/

#ifndef INC_PLOT
#define INC_PLOT

#include <string>
#include <vector>

#include "boost/lexical_cast.hpp"

#include "../event/event.h"
#include "../histogram/histogram.h"


/* NAMESPACE */
namespace analysis
{
	
	/* cut base class */
	class plot_default
	{
		
	public:
		// constructor can be default
		plot_default() = default;
		// virtual destructor needed for being a overloadable class
		virtual ~plot_default() {};
		// all cuts are based on this class and should overload this
		// operator, and return true for events that pass the cuts
		virtual double operator() (const event *ev) { return 0; }
	
	};

	/* plotting class */
	class plot
	{

	public:
	
		/* con & destructor */
		plot(std::string name = "", std::string folder = "");
		~plot();

		/* plot data */
		void add_sample(const std::vector<event*> &events, plot_default *plot_imp, const std::string &name = "", double weight = 1);
		void run(); 
		
		/* plot properties */
		void set_name(std::string n);
		void set_folder(std::string f);

		void set_x_label(std::string x);
		void set_y_label(std::string y);
		void set_leg_title(std::string t);

		void set_logy(bool on);
		void set_stacked(bool on);
		void set_normalized(bool on);
		
		void set_bins(double nbins, double min, double max);

	private:

		histogram hist;
		std::string plot_name;
		std::string plot_folder;

	};
	
/* NAMESPACE */
}

/* CUT IMPLEMENTATIONS */
#include "plot_default.h"

#endif
