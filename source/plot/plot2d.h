/* 2D Plot class
 *
 *
*/

#ifndef INC_PLOT2D
#define INC_PLOT2D

#include <string>
#include <vector>

#include "boost/lexical_cast.hpp"

#include "../event/event.h"
#include "../histogram/histogram2D.h"
#include "../plot/plot.h"


/* NAMESPACE */
namespace analysis
{

	/* plot base class */
	class plot2d
	{

	public:
	
		/* con & destructor */
		plot2d(std::string name = "", std::string folder = "");
		~plot2d();

		/* plot data */
		void add_sample(const std::vector<event*> &events, const std::string &name = "", double weight = 1);
		void run(); 
		
		/* plot properties */
		std::string name() const;
		void set_name(std::string n);
		void set_folder(std::string f);
		
		void set_functors(plot_default *x, plot_default *y);
		
		void set_x_bins(double nbins, double min, double max);
		void set_y_bins(double nbins, double min, double max);

	private:

		histogram2D hist;
		std::string plot_name;
		std::string plot_folder;
		
		plot_default *x_ftor;
		plot_default *y_ftor;

	};

/* NAMESPACE */
}

#endif
