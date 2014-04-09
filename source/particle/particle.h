/* Particle base class
 *
 * Provides a virtual base class for all particle/object formats.
*/

#ifndef INC_PARTICLE
#define INC_PARTICLE

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../../deps/gzstream/gzstream.h"


/* NAMESPACE */
namespace analysis
{
	
	// constants for the particle types
	const int
		pid_nutbar     = -16,
		pid_taubar     = -15,
		pid_numbar     = -14,
		pid_muonbar    = -13,
		pid_nuebar     = -12,
		pid_positron   = -11,
		pid_topbar     = -6,
		pid_bottombar  = -5,
		pid_charmbar   = -4,
		pid_strangebar = -3,
		pid_upbar      = -2,
		pid_downbar    = -1,
		pid_none       = 0,
		pid_down       = 1,
		pid_up         = 2,
		pid_strange    = 3,
		pid_charm      = 4,
		pid_bottom     = 5,
		pid_top        = 6,
		pid_electron   = 11,
		pid_nue        = 12,
		pid_muon       = 13,
		pid_num        = 14,
		pid_tau        = 15,
		pid_nut        = 16,
		pid_gluon      = 21,
		pid_photon     = 22,
		pid_z          = 23,
		pid_w          = 24,
		pid_h          = 25;
	
	// constants for the particle types
	const unsigned int 
		ptype_none      = 0,
		ptype_photon    = 1,
		ptype_electron  = 1 << 1,
		ptype_muon      = 1 << 2,
		ptype_tau       = 1 << 3,
		ptype_jet       = 1 << 4,
		ptype_met       = 1 << 6,
		ptype_lepton    = ptype_electron | ptype_muon,
		ptype_leptonall = ptype_lepton | ptype_tau,
		ptype_all       = ~ptype_none;

	class particle
	{
		
	public:

		/* con & destructor */
		particle();
		virtual ~particle();

		/* copy & assignment */
		particle(const particle&) = default;
    	particle& operator = (const particle&) = default;

		/* properties: type */
		virtual int id() const;
		virtual unsigned int type() const;
		
		/* properties: state */
		virtual void set_final(bool is_final);
		virtual bool is_final() const;
		
		/* properties: quantum numbers */
		virtual double charge() const;
		virtual double bjet() const;

		/* properties: kinematics */
		virtual double px() const;
		virtual double py() const;
		virtual double pz() const;
		virtual double pe() const;

		virtual double pt() const;
		virtual double eta() const;
		virtual double phi() const;
		virtual double mass() const;

		virtual double y() const;

		/* input & output */
		virtual void write(std::ostream& os) const;
		virtual void read(std::ostream& is);
		virtual void write(std::ofstream& ofs) const;
		virtual void read(std::ifstream& ifs);
		virtual void write(ogzstream& ogzs) const;
		virtual void read(igzstream& igzs);

	};
	
	/* utility functions */
	
	bool compare_pt(particle *p1, particle *p2);
	bool compare_type(particle *p1, particle *p2);
	std::string ptype_to_string(unsigned int type);
	
	double delta_eta(const particle *p1, const particle *p2);
	double delta_phi(const particle *p1, const particle *p2);
	double delta_r(const particle *p1, const particle *p2);

/* NAMESPACE */
}

#endif
