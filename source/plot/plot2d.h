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


/* NAMESPACE */
namespace analysis
{

	/* plot base class */
	class plot2d
	{

	public:
	
		/* con & destructor */
		plot2d(std::string new_name = "");
		~plot2d();

		/* plot properties */
		virtual void add_sample(const std::vector<event*> &events, const std::string &name, double weight = 1);
		void run(); 
		
		/* plot properties */
		virtual std::string name() const;
		void set_name(std::string n);
		void set_folder(std::string f);
		
		void set_x_bins(double nbins, double min, double max);
		void set_y_bins(double nbins, double min, double max);

	public:

		histogram2D hist;
		std::string plot_name;
		std::string plot_folder;

	private:

	};
	
	/* plot ht vs met class */
	class plot2d_ht_met : public plot2d
	{
		
	public:
	
		plot2d_ht_met();

		void add_sample(const std::vector<event*> &events, const std::string &name, double weight = 1); 

		/* plot properties */
		std::string name() const; 		
	};

/* NAMESPACE */
}

#endif
