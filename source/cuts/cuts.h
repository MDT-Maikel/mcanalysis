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
		cut() = default;
		virtual ~cut() {};
		virtual bool operator() (const event *ev) { return false; }
	
	};
	
	/* cut & count class */
	class cuts 
	{

	public:
		cuts();
		
		void add_cut(cut *add, std::string n = "");
		void apply(std::vector<event*> &events);
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

/* TEMPLATE IMPLEMENTATIONS */
#include "cuts_imp.h"

#endif
