/* Standard plot classes
 *
 *
*/

#ifndef INC_PLOT_STANDARD
#define INC_PLOT_STANDARD

#include <string>
#include <vector>

#include "boost/lexical_cast.hpp"

#include "plot.h"
#include "../event/event.h"
#include "../histogram/histogram.h"


/* NAMESPACE */
namespace analysis
{

	/* plot pt class */
	class plot_pt : public plot
	{
		
	public:
		plot_pt(int t, unsigned int n);

		void add_sample(const std::vector<event*> &events, const std::string &name); 

		/* plot properties */
		std::string name() const; 

	private:
		int type;
		unsigned int number;

	};

	/* plot met class */
	class plot_met : public plot
	{
		
	public:
		plot_met();

		void add_sample(const std::vector<event*> &events, const std::string &name); 

		/* plot properties */
		std::string name() const; 

	private:

	};

	/* plot ht class */
	class plot_ht : public plot
	{
		
	public:
		plot_ht();

		void add_sample(const std::vector<event*> &events, const std::string &name); 

		/* plot properties */
		std::string name() const; 

	private:

	};

	/* plot mass class */
	class plot_mass : public plot
	{
		
	public:
		plot_mass();

		void add_sample(const std::vector<event*> &events, const std::string &name); 

		/* plot properties */
		std::string name() const; 
		void set_type(int t);
		void set_comb(const std::vector<int> &p);

	private:
		int type;
		std::vector<int> comb;
	};

/* NAMESPACE */
}

#endif
