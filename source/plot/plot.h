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
		plot();
		~plot();
		
		virtual void run(std::vector<event*> &events); 
		
		/* plot properties */
		virtual std::string name() const;

	public:
		histogram hist;

	private:

	};

	/* plot pt class */
	class plot_pt : public plot
	{
		
	public:
		plot_pt(int t, unsigned int n);

		void run(std::vector<event*> &events); 

		/* plot properties */
		std::string name() const; 

	private:
		int type;
		unsigned int number;

	};

/* NAMESPACE */
}

#endif
