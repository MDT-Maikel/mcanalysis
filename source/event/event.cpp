/* Event class
 *
 * 
*/

#include "event.h" 

/* NAMESPACE */
namespace analysis
{

	/* con & destructor */

	event::~event() 
	{
		for (unsigned int i = 0; i < particles.size(); i++) 
			delete particles[i];
		particles.clear();
	}

	/* member operations */

	particle* event::operator[] (int n) 
	{
		return particles[n]; 
	}

	const particle* event::operator[] (int n) const 
	{
		return particles[n]; 
	}

	unsigned int event::size() 
	{ 
		return particles.size(); 
	}

	void event::resize(unsigned int n) 
	{ 
		particles.resize(n); 
	}

	void event::push_back(particle *p) 
	{ 
		particles.push_back(p); 
	}

	void event::clear() 
	{ 
		particles.clear(); 
	}

	/* input & output */

	void event::write(std::ostream& os) const
	{
		os << "test" << std::endl;
		for (unsigned int i = 0; i < particles.size(); i++)
			particles[i]->write(os);
	}

	void event::read(std::istream& is)
	{


	}

	void event::write(std::ofstream& ofs) const
	{
		ofs << "test" << std::endl;
		for (unsigned int i = 0; i < particles.size(); i++)
			particles[i]->write(ofs);
	}

	void event::read(std::ifstream& ifs, particle *type)
	{
		particles.clear();		
		// assume that ifstream is at the start of an event
		//int temp;
		while (ifs)
		{
			//ifs >> temp;
			//std::cout << temp << std::endl;			
			lhe *p = new lhe;
			p->read(ifs);
			particles.push_back(p);
		}
		std::cout << particles.size() << std::endl;
	}





/* NAMESPACE */
}
