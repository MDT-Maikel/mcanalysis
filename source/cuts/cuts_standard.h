/* Standard cut classes
 *
 * Provides a set of standard cuts by overloading the 
 * cut base class.
*/

#ifndef INC_CUTS_STANDARD
#define INC_CUTS_STANDARD

#include <cmath>
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>

#include "../event/event.h"
#include "cuts.h"


/* NAMESPACE */
namespace analysis
{

	/* pt cut class */
	class cut_pt : public cut
	{
	
	public:
		cut_pt(int t = 0, unsigned int n = 0, double pt = 0, double eta = 5);
	
		bool passed(const event *ev);
		
		std::string name() const;

		void set_type(int to_type);
		void set_number(unsigned int to_number);
		void set_pt(double pt);
		void set_eta(double eta);

	private:
		int type;
		unsigned int number;
		double ptcut;
		double etamax;

	};

	/* MET cut class */
	class cut_met : public cut
	{

	public:
		cut_met(double met = 0);

		bool passed(const event *ev);

		std::string name() const;
		
		void set_met(double met);

	private:
		double metcut;

	};

	/* HT cut class */
	class cut_ht : public cut
	{

	public:
		cut_ht(int t = 0, double pt = 0, double eta = 0, double ht = 0);
	
		bool passed(const event *ev);

		std::string name() const;

		void set_ht(double ht);

	private:
		int type;
		double min_pt;
		double max_eta;
		double htcut;		

	};

	/* veto cut class */
	class cut_veto : public cut
	{

	public:
		cut_veto(int t = 0, double pt = 0, double eta = 0);

		bool passed(const event *ev);

		std::string name() const;

	private: 
		int type;
		double min_pt;
		double max_eta;

	};

/* NAMESPACE */
}

#endif
