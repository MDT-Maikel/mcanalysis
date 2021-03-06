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
		// constructor can be default
		cut() = default;
		// virtual destructor needed for being a overloadable class
		virtual ~cut() {};
		// all cuts are based on this class and should overload this
		// operator, and return true for events that pass the cuts
		virtual bool operator() (const event *ev) { return false; }
	
	};
	
	/* cut & count class */
	class cuts 
	{

	public:
		cuts();
		
		void add_cut(cut *add, std::string n = "");
		void apply(std::vector<event*> &events);
		const std::vector<event*> reduce(const std::vector<event*> &events) const;
		double efficiency() const;
		double efficiency(unsigned int p, unsigned int t) const;
		void clear();
	
		void write(std::ostream& os) const;
		void write(std::ofstream& ofs) const;	

	private:
		std::vector<cut*> list_cuts;
		std::vector<std::string> list_names;
		std::vector<unsigned int> list_total;
		std::vector<unsigned int> list_pass;
		unsigned int total;
		unsigned int pass;

	};
	
/* NAMESPACE */
}

/* CUT IMPLEMENTATIONS */
#include "cuts_default.h"

#endif
