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

		virtual void add_sample(const std::vector<event*> &events, const std::string &name);
		
		void run(); 
		
		/* plot properties */
		virtual std::string name() const;

	public:
		histogram hist;

	private:

	};

/* NAMESPACE */
}

#endif
