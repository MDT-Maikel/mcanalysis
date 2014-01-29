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

	/* properties */

	int particle::id() const { return 0; }
	int particle::type() const { return 0; }
	void particle::set_final(bool is_final) { /* to be overloaded */ };
	bool particle::is_final() const { return true; }

	/* kinematics */

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

/* NAMESPACE */
}
