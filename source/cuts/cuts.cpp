/* Cut & Count classes
 *
 * Provides a class which cuts & counts in a list of events.
 * This class is based on a simple virtual cut class which 
 * can be overloaded by specialised cuts.
*/

#include "cuts.h"


/* NAMESPACE */
namespace analysis
{

	/* cut & count class */

	cuts::cuts()
	{
		total = 0;
		pass = 0;
	}

	void cuts::add_cut(cut *add, std::string n)
	{
		list_cuts.push_back(add);;
		list_names.push_back(n);
		list_total.push_back(0);
		list_pass.push_back(0);
	}

	void cuts::apply(std::vector<event*> &events) 
	{
		// store the total number of events	
		total = events.size();		
	
		// loop over all cuts and events
		for (unsigned int i = 0; i < list_cuts.size(); i++)
		{
			// get the cut
			cut *apply_cut = list_cuts[i];

			// loop over all events
			for (int index = events.size() - 1; index >= 0; index--)
			{	
				event *ev = events[index];
				// test the cut on the event and determine whether it passed				
				if ((*apply_cut)(ev))
				{
					list_pass[i]++;
				}
				else
				{
					// delete also the pointer to the event
					delete ev;
					events.erase(events.begin() + index);
				}
				list_total[i]++;
			}
		}

		// store the passed number of events
		pass = events.size();
	}

	double cuts::efficiency() const
	{
		if (total == 0)
			return 0.0;
		return static_cast<double>(pass) / total;
	}
	
	double cuts::efficiency(unsigned int p, unsigned int t) const
	{
		if (t == 0)
			return 0.0;
		return static_cast<double>(p) / t;
	}

	void cuts::clear()
	{
		total = 0;
		pass = 0;
		// also clear the different totals in the lists
		for (unsigned int i = 0; i < list_cuts.size(); i++)
		{
			list_pass[i] = 0;
			list_total[i] = 0;
		}
	}

	void cuts::write(std::ostream& os) const
	{
		os << "Efficiencies for each of the different cuts:" << std::endl;		
		for (unsigned int i = 0; i < list_cuts.size(); i++)
		{
			unsigned int p = list_pass[i];
			unsigned int t = list_total[i];
			os << "cut: " << list_names[i] << " -> efficiency: " << 100 * efficiency(p, t) << "%";
			os << " (" << p << "/" << t << ")" << std::endl;
		}
		os << "total efficiency: " << 100 * efficiency() << "%";
		os << " (" << pass << "/" << total << ")" << std::endl;
	}

	void cuts::write(std::ofstream& ofs) const
	{
		ofs << "Efficiencies of different cuts:" << std::endl;		
		for (unsigned int i = 0; i < list_cuts.size(); i++)
		{
			unsigned int p = list_pass[i];
			unsigned int t = list_total[i];
			ofs << "cut: " << list_names[i] << " -> efficiency: " << 100 * efficiency(p, t) << "%";
			ofs << " (" << p << "/" << t << ")" << std::endl;
		}
		ofs << "total efficiency: " << 100 * efficiency() << "%";
		ofs << " (" << pass << "/" << total << ")" << std::endl;
	}

/* NAMESPACE */
}
