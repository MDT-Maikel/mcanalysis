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
		// sort particles according to pt
		sort_pt();
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

	particle* event::get(unsigned int type, unsigned int number) const
	{
		// TODO: safety checks for number		
		unsigned int count = 0;		
		for (unsigned int index = 0; index < size(); index++)
		{
			if (particles[index]->type() & type)
			{
				count++;
				if (count == number)
					return particles[index];
			}
		}
		return nullptr;
	}

	particle* event::get(unsigned int type, unsigned int number, double max_eta) const
	{
		// TODO: safety checks for number	
		unsigned int count = 0;		
		for (unsigned int index = 0; index < size(); index++)
		{
			if ((particles[index]->type() & type) && std::abs(particles[index]->eta()) < max_eta)
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
		double px_inv = 0;
		double py_inv = 0;
		for (unsigned int index = 0; index < size(); index++)
		{
			particle *p = particles[index];
			if (p->is_final() && (p->type() & ptype_met))
			{
				px_inv += p->px();
				py_inv += p->py();				
			}
		}
		return std::sqrt(px_inv * px_inv + py_inv * py_inv);
	}

	// returns the ht of all particles satisfying the conditions
	double event::ht(unsigned int type, double min_pt, double max_eta) const
	{
		double ht = 0;
		for (unsigned int index = 0; index < size(); index++)
		{
			particle *p = particles[index];
			if (p->type() & type)
				if (p->is_final() && p->pt() > min_pt && std::abs(p->eta()) < max_eta)
					ht += p->pt();
		}
		return ht;
	}

	// returns the invariant mass of all final state particles in the event
	double event::mass() const
	{
		double pe = 0.0; double px = 0.0; double py = 0.0; double pz = 0.0;		
		for (unsigned int index = 0; index < size(); index++)
		{
			particle *p = particles[index];
			if (p->is_final())
			{
				pe += p->pe();
				px += p->px();
				py += p->py();
				pz += p->pz();
			}
		}
		double inv_mass = std::pow(pe, 2.0) - std::pow(px, 2.0) - std::pow (py, 2.0) - std::pow(pz, 2.0);
		return std::sqrt(std::max(inv_mass, 0.0));
	}

	// returns the invariant mass of the combination of particle of the requested type
	double event::mass(unsigned int type, const std::vector<int> &comb) const
	{
		double pe = 0.0; double px = 0.0; double py = 0.0; double pz = 0.0;
		unsigned int count = 0;	
		for (unsigned int index = 0; index < size(); index++)
		{
			particle *p = particles[index];
			if ((p->type() & type) && p->is_final())
			{
				count++;
				if (std::find(comb.begin(), comb.end(), count) != comb.end())
				{
					pe += p->pe();
					px += p->px(); 
					py += p->py();
					pz += p->pz();
				}
			}
		}
		double inv_mass = std::pow(pe, 2.0) - std::pow(px, 2.0) - std::pow (py, 2.0) - std::pow(pz, 2.0);
		return std::sqrt(std::max(inv_mass, 0.0));
	}
	
	// returns the mt2 for the event
	double event::mt2(double mn) const
	{
		// identify leptons within the event
		std::vector< particle* > leptons;
		for (unsigned int index = 0; index < size(); index++)
		{
			particle *p = particles[index];
			if (p->is_final() && (p->type() & ptype_lepton))
				leptons.push_back(p);
		}
		// TODO: generalise?
		// at least two leptons has to be present in order to calculate mT2
		if (leptons.size() < 2)
			return 0;

		// extract MET components
		double px_inv = 0;
		double py_inv = 0;
		for (unsigned int index = 0; index < size(); index++)
		{
			particle *p = particles[index];
			if (p->is_final() && p->type() & ptype_met)
			{
				px_inv += p->px();
				py_inv += p->py();				
			}
		}

		// calculate mT2 using the two leading leptons, MET and mn as input values
		double pa[3]    = { 0, leptons[0]->px(), leptons[0]->py() };
		double pb[3]    = { 0, leptons[1]->px(), leptons[1]->py() };
		double pmiss[3] = { 0, px_inv, py_inv };

		mt2_bisect::mt2 mt2_event;
		mt2_event.set_momenta(pa,pb,pmiss);
		mt2_event.set_mn(mn);
		double mt2_value = mt2_event.get_mt2();
		
		return mt2_value;

	}
	
	/* utility */
	
	// sorts the current particles in the event by its pt, highest pt first
	void event::sort_pt()
	{
		std::sort(particles.begin(), particles.end(), compare_pt);
	}
	
	void event::sort_type()
	{
		std::sort(particles.begin(), particles.end(), compare_type);
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
	
	/* utility functions */
	
	double mass(std::vector<particle*> particles)
	{
		double pe = 0.0; double px = 0.0; double py = 0.0; double pz = 0.0;		
		for (unsigned int index = 0; index < particles.size(); index++)
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
	
	void delete_events(std::vector<event*> & events)
	{
		for (unsigned int i = 0; i < events.size(); ++i)
		{
			delete events[i];			
		}
		events.clear();				
	}

/* NAMESPACE */
}
