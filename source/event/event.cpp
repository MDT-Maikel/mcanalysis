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
		{
			delete particles[i];
			particles[i] = nullptr;
		}
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

	unsigned int event::size() const
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

	void event::erase(int n)
	{
		particles.erase(particles.begin() + n);
	}

	void event::clear() 
	{ 
		particles.clear(); 
	}

	/* member access */

	particle* event::get(int type, unsigned int number) const
	{
		// TODO: safety checks for number		
		unsigned int count = 0;		
		for (unsigned int index = 0; index < size(); index++)
		{
			if ((particles[index])->type() == type)
			{
				count++;
				if (count == number)
					return particles[index];
			}
		}
		return nullptr;
	}

	particle* event::get(int type, unsigned int number, double max_eta) const
	{
		// TODO: safety checks for number	
		unsigned int count = 0;		
		for (unsigned int index = 0; index < size(); index++)
		{
			if (particles[index]->type() == type && std::abs(particles[index]->eta()) < max_eta)
			{
				count++;
				if (count == number)
					return particles[index];
			}
		}
		return nullptr;
	}

	/* kinematics */

	double event::met() const
	{
		double met = 0;
		for (unsigned int index = 0; index < size(); index++)
		{
			particle *p = particles[index];
			if (p->type() == particle::type_met)
				met += p->pt();
		}
		return met;
	}

	double event::ht(int type, double min_pt, double max_eta) const
	{
		double ht = 0;
		for (unsigned int index = 0; index < size(); index++)
		{
			particle *p = particles[index];
			if (p->type() == type && p->pt() > min_pt && std::abs(p->eta()) < max_eta)
				ht += p->pt();
		}
		return ht;
	}

	double event::mass() const
	{
		double pe = 0.0; double px = 0.0; double py = 0.0; double pz = 0.0;		
		for (unsigned int index = 0; index < size(); index++)
		{
			particle *p = particles[index];
			pe += p->pe();
			px += p->px();
			py += p->py();
			pz += p->pz();
		}
		double inv_mass = std::pow(pe, 2.0) - std::pow(px, 2.0) - std::pow (py, 2.0) - std::pow(pz, 2.0);
		return std::sqrt(std::max(inv_mass, 0.0));
	}

	/* input & output */

	void event::write(std::ostream& os) const
	{
		for (unsigned int i = 0; i < particles.size(); i++)
			particles[i]->write(os);
	}

	void event::read(std::istream& is)
	{


	}

	void event::write(std::ofstream& ofs) const
	{
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
