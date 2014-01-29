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
#include "../plot/plot.h"
#include "../histogram/histogram2D.h"


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
		
		void set_functions(double(*x)(const event*), double(*y)(const event*));
		
		void set_x_bins(double nbins, double min, double max);
		void set_y_bins(double nbins, double min, double max);

	private:

		histogram2D hist;
		std::string plot_name;
		std::string plot_folder;
		
		double (*x_func)(const event*);
		double (*y_func)(const event*);

	};
	


/* NAMESPACE */
}

#endif
