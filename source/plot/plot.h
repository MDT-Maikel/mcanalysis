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
		plot();
		~plot();

		/* plot properties */
		virtual void add_sample(const std::vector<event*> &events, const std::string &name);
		void run(); 
		
		/* plot properties */
		virtual std::string name() const;

		void set_logy(bool on);
		void set_stacked(bool on);
		void set_normalized(bool on);

	public:

		histogram hist;

	private:

	};

/* NAMESPACE */
}

#endif
