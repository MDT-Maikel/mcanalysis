/* Standard cut classes
 *
 * Provides a set of standard cuts by overloading the 
 * cut base class.
*/

#include "cuts_standard.h"


/* NAMESPACE */
namespace analysis
{

	/* pt cut class */

	cut_pt::cut_pt(int t, unsigned int n, double pt, double eta)
	{
		type = t;
		number = n;
		ptcut = pt;
		etamax = eta;
	}

	void cut_pt::init()
	{
	}

	bool cut_pt::passed(const event *ev) 
	{
		// increase total number of events which went through this cut.	
		increase_total();
		//std::cout << "called pt cut" << std::endl;
		const particle *p = ev->get(type, number, etamax); 
		if (!p)
			return false;
		//std::cout << "type: " << p->type() << ", pt: " << p->pt() << std::endl;
		bool has_passed = p->pt() > ptcut;
		if (has_passed)
			increase_passed();
		return has_passed; 
	}

	std::string cut_pt::name() const
	{
		return "for jet " + std::to_string(number) + ": pt > " + boost::lexical_cast<std::string>(ptcut) + " & eta < " + boost::lexical_cast<std::string>(etamax);
	}

	void cut_pt::set_type(int to_type)
	{
		type = to_type;
	}

	void cut_pt::set_number(unsigned int to_number)
	{
		number = to_number;
	}

	void cut_pt::set_pt(double pt)
	{
		ptcut = pt;
	}

	void cut_pt::set_eta(double eta)
	{
		etamax = eta;
	}

	/* MET cut class */

	cut_met::cut_met(double met)
	{
		metcut = met;
	}

	void cut_met::init()
	{
	}

	bool cut_met::passed(const event *ev)
	{
		// increase total number of events which went through this cut.	
		increase_total();
		//std::cout << "called MET cut: " << ev->met() << std::endl;
		bool has_passed = ev->met() > metcut;
		if (has_passed)
			increase_passed();
		return has_passed;
	}

	std::string cut_met::name() const
	{
		return "met > " + boost::lexical_cast<std::string>(metcut);
	}

	void cut_met::set_met(double met)
	{
		metcut = met;
	}

	/* HT cut class */

	cut_ht::cut_ht(int t, double pt, double eta, double ht)
	{
		type = t;
		min_pt = pt;
		max_eta = eta;
		htcut = ht;
	}

	void cut_ht::init()
	{
	}

	bool cut_ht::passed(const event *ev)
	{
		// increase total number of events which went through this cut.	
		increase_total();
		//std::cout << "called HT cut: " << ev->ht(type, min_pt, max_eta) << std::endl;
		bool has_passed = ev->ht(type, min_pt, max_eta) > htcut;
		if (has_passed)
			increase_passed();
		return has_passed;
	}

	std::string cut_ht::name() const
	{
		return "ht > " + boost::lexical_cast<std::string>(htcut);
	}

	void cut_ht::set_ht(double ht)
	{
		htcut = ht;
	}

	/* veto cut class */

	cut_veto::cut_veto(int t, double pt, double eta)
	{
		type = t;
		min_pt = pt;
		max_eta = eta;
	}

	void cut_veto::init()
	{
	}

	bool cut_veto::passed(const event *ev)
	{
		// increase total number of events which went through this cut.	
		increase_total();
		unsigned int index = 1;
		bool has_passed = true;
		const particle *p;
		while ((p = ev->get(type, index)))
		{
			if (p->pt() > min_pt && std::abs(p->eta()) < max_eta) 
			{
				has_passed = false;
				break;
			}
			index++;
		}
		if (has_passed)
			increase_passed();
		return has_passed;
	}

	std::string cut_veto::name() const
	{
		return "veto(" + boost::lexical_cast<std::string>(type) + "): pt > " + boost::lexical_cast<std::string>(min_pt) + " & eta < " + boost::lexical_cast<std::string>(max_eta);
	}

/* NAMESPACE */
}
