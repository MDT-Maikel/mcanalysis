/* Event class
 *
 * 
*/

#ifndef INC_EVENT
#define INC_EVENT

#include <algorithm>
#include <cmath>
#include <vector>

#include "../particle/particle.h"
#include "../particle/lhco.h"
#include "../particle/lhe.h"
#include "mt2_bisect.h"


/* NAMESPACE */
namespace analysis
{

	class event
	{

	public:
		
		/* con & destructor */
		event() = default;
		~event();
		
		/* copy & assignment */
		event(const event&) = default;
		event & operator = (const event&) = default;

		/* member operations */
		particle* operator[] (int n);
		const particle* operator[] (int n) const;
		unsigned int size() const;
		void resize(unsigned int n);
		void push_back(particle *p);
		void erase(int n);
		void clear();

		/* member access */
		particle* get(unsigned int type, unsigned int number) const;
		particle* get(unsigned int type, unsigned int number, double max_eta) const;

		/* kinematics */
		double met() const;
		double ht(unsigned int type, double min_pt, double max_eta) const;
		double mass() const;
		double mass(unsigned int type, const std::vector<int> &comb) const;
		double mt2(double mn = 0) const;
		
		/* utility */
		void sort_pt();

		/* input & output */
		void write(std::ostream& os) const;
		void read(std::istream& is);
		void write(std::ofstream& ofs) const;
		void read(std::ifstream& ifs, particle *type);

	private:

		/* event data members */
		std::vector<particle*> particles;

	};
	
	/* utility functions */
	double mass(std::vector<particle*> particles);

/* NAMESPACE */
}

#endif
