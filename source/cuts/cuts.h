/* Cut & Count classes
 *
 * Provides a class which cuts & counts in a list of events.
 * This class is based on a simple virtual cut class which 
 * can be overloaded by specialised cuts.
*/

#ifndef INC_CUTS
#define INC_CUTS

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../event/event.h"


/* NAMESPACE */
namespace analysis
{

	/* cut base class */
	class cut
	{
		
	public:
		cut();
		virtual ~cut();

		virtual bool passed(const event *ev);

		double efficiency() const;
		virtual std::string name() const;

		// methods to track the efficiency
		void increase_total();
		void increase_passed();
		void clear();

		unsigned int get_total() const;
		unsigned int get_passed() const;

	private:
		// number of events that have been handled by the cut
		unsigned int total;
		// number of events that have passed the cut
		unsigned int pass;

	};

	/* cut & count class */
	class cuts 
	{

	public:
		cuts();
		
		void add_cut(cut *add);
		void apply(std::vector<event*> &events);
		void clear();
	
		void write(std::ostream& os) const;
		void write(std::ofstream& ofs) const;	

	private:
		std::vector<cut*> cut_list;
		unsigned int total;
		unsigned int pass;

	};

/* NAMESPACE */
}

#endif
