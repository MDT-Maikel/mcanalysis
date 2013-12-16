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

	/* cut base class */

	cut::cut()
	{
		total = 0;
		pass = 0;
	}

	cut::~cut()
	{
	}

	bool cut::passed(const event *ev)
	{
		std::cout << "WARNING: called virtual cut base class." << std::endl; 
		total++;
		pass++;
		return true; 
	}

	double cut::efficiency() const
	{
		if (total == 0)
			return 0.0;
		return static_cast<double>(pass) / total;
	}

	std::string cut::name() const
	{
		return "WARNING: using base cut class";
	}

	void cut::increase_total()
	{
		total++;
	}

	void cut::increase_passed()
	{
		pass++;
	}

	void cut::clear()
	{
		total = 0;
		pass = 0;
	}

	unsigned int cut::get_total() const
	{
		return total;
	}

	unsigned int cut::get_passed() const
	{
		return pass;
	}

	/* cut & count class */

	cuts::cuts()
	{
		total = 0;
		pass = 0;
	}

	void cuts::add_cut(cut *add)
	{
		cut_list.push_back(add);
	}

	void cuts::apply(std::vector<event*> &events)
	{
		for (int index = events.size() - 1; index >= 0; index--)
		{				
			total++;			

			bool passed = true;
			event *ev = events[index];

			for (unsigned int i = 0; i < cut_list.size(); i++)
			{
				cut *apply_cut = cut_list[i];
				if (!apply_cut->passed(ev))
				{
					passed = false;
					break;
				}
			}
			if (!passed)
			{
				// delete also the pointer to the event.
				delete events[index];
				events.erase(events.begin() + index);
			}
			else
			{
				pass++;
			}
		}
	}

	void cuts::clear()
	{
		total = 0;
		pass = 0;
		// also clear the different cuts in the list
		for (unsigned int i = 0; i < cut_list.size(); i++)
			cut_list[i]->clear();		
	}

	void cuts::write(std::ostream& os) const
	{
		os << "Efficiencies for each of the different cuts:" << std::endl;		
		for (unsigned int i = 0; i < cut_list.size(); i++)
		{
			cut *print_cut = cut_list[i];
			os << "cut: " << cut_list[i]->name() << " -> efficiency: " << 100.0 * print_cut->efficiency() << "%";
			os << " (" << print_cut->get_passed() << "/" << print_cut->get_total() << ")" << std::endl;
		}
		double efficiency = 0.0;
		if (total != 0)
			efficiency = 100.0 * static_cast<double>(pass) / total;
		os << "total efficiency: " << efficiency << "%";
		os << " (" << pass << "/" << total << ")" << std::endl;
	}

	void cuts::write(std::ofstream& ofs) const
	{
		ofs << "Efficiencies of different cuts:" << std::endl;		
		for (unsigned int i = 0; i < cut_list.size(); i++)
		{
			cut *print_cut = cut_list[i];
			ofs << "cut: " << cut_list[i]->name() << " -> efficiency: " << 100.0 * print_cut->efficiency() << "%";
			ofs << " (" << print_cut->get_passed() << "/" << print_cut->get_total() << ")" << std::endl;
		}
		double efficiency = 0.0;
		if (total != 0)
			efficiency = 100.0 * static_cast<double>(pass) / total;
		ofs << "total efficiency: " << efficiency << "%";
		ofs << " (" << pass << "/" << total << ")" << std::endl;
	}

/* NAMESPACE */
}
