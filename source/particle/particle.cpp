/* Particle base class
 *
 * Provides a virtual base class for all particle/object formats.
*/

#include "particle.h"


/* NAMESPACE */
namespace analysis
{

	/* con & destructor */
	
	particle::particle() {}
	particle::~particle() {}

	/* properties: type */

	int particle::id() const { return 0; }
	unsigned int particle::type() const { return ptype_none; }
	
	/* properties: state */

	void particle::set_final(bool is_final) { /* to be overloaded */ }
	bool particle::is_final() const { return true; }
	
	/* properties: quantum numbers */
	
	double particle::charge() const { return 0; }

	/* properties: kinematics */

	double particle::px() const { return 0; }
	double particle::py() const { return 0; }
	double particle::pz() const { return 0; }
	double particle::pe() const { return 0; }

	double particle::pt() const { return 0; }
	double particle::eta() const { return 0; }
	double particle::phi() const { return 0; }
	double particle::mass() const { return 0; }

	double particle::y() const { return 0; }

	/* input & output */

	void particle::write(std::ostream& os) const {}
	void particle::read(std::ostream& is) {}
	void particle::write(std::ofstream& ofs) const {}
	void particle::read(std::ifstream& ifs) {}
	void particle::write(ogzstream& ogzs) const {}
	void particle::read(igzstream& igzs) {}
	
	/* utility functions */
	
	bool compare_pt(particle *p1, particle *p2)
	{
		return p1->pt() > p2->pt();		
	}
	
	bool compare_type(particle *p1, particle *p2)
	{
		return p1->type() < p2->type();		
	}
	
	std::string ptype_to_string(unsigned int type)
	{
		// basic types
		if (type == ptype_none)
			return "none";
		if (type == ptype_all)
			return "all";
		
		// other types might be concatenated 
		std::string stype = "";
		
		if (type & ptype_lepton)
		{
			stype += stype.size() > 0 ? "|l" : "l";
		}
		else
		{
			if (type & ptype_electron)
				stype += stype.size() > 0 ? "|e" : "e";
			if (type & ptype_muon)
				stype += stype.size() > 0 ? "|mu" : "mu";
		}		
		if (type & ptype_tau)
			stype += stype.size() > 0 ? "|tau" : "tau";
		if (type & ptype_jet)
			stype += stype.size() > 0 ? "|jet" : "jet";
		if (type & ptype_met)
			stype += stype.size() > 0 ? "|met" : "met";	
		
		return stype;
	}

/* NAMESPACE */
}
