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
		unsigned int size();
		void resize(unsigned int n);
		void push_back(particle *p);
		void clear();

		/* input & output */
		void write(std::ostream& os) const;
		void read(std::istream& is);
		void write(std::ofstream& ofs) const;
		void read(std::ifstream& ifs, particle *type);

	private:

		/* event data members */
		std::vector<particle*> particles;

	};

/* NAMESPACE */
}

#endif
