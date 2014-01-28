/* Particle base class
 *
 * Provides a virtual base class for all particle/object formats.
*/

#ifndef INC_PARTICLE
#define INC_PARTICLE

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "../../deps/gzstream/gzstream.h"


/* NAMESPACE */
namespace analysis
{

	class particle
	{
	public:

		/* con & destructor */
		particle();
		virtual ~particle();

		/* copy & assignment */
		particle(const particle&) = default;
    	particle& operator = (const particle&) = default;

	public:

		/* constants */
		// particle ids: from the pdg codes
		enum particle_id 
		{
			nuebar = -12, positron,	topbar = -6, bottombar, charmbar, strangebar, upbar, downbar, gluon, up, down, strange, charm, bottom, top, electron = 11, nue
		};
		enum particle_type
		{
			type_unknown = -1, type_photon = 0, type_electron = 1, type_muon = 2, type_tau = 3, type_jet = 4, type_met = 6
		};

	public:

		/* properties */
		virtual int id() const;
		virtual int type() const;
		virtual bool is_final() const;

		/* kinematics */
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

/* NAMESPACE */
}

#endif
