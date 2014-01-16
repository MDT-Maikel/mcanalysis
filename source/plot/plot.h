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

	/* plot base class */
	class plot
	{

	public:
	
		/* con & destructor */
		plot(std::string new_name = "");
		~plot();

		/* plot properties */
		virtual void add_sample(const std::vector<event*> &events, const std::string &name, double weight = 1);
		void run(); 
		
		/* plot properties */
		virtual std::string name() const;
		void set_name(std::string n);
		void set_folder(std::string f);

		void set_logy(bool on);
		void set_stacked(bool on);
		void set_normalized(bool on);
		
		void set_bins(double nbins, double min, double max);

	public:

		histogram hist;
		std::string plot_name;
		std::string plot_folder;

	private:

	};

/* NAMESPACE */
}

#endif
