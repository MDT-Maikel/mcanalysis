/* LHE particle class
 *
 *
*/

#ifndef INC_LHE
#define INC_LHE

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include <boost/lexical_cast.hpp>

#include "../../deps/gzstream/gzstream.h"

#include "particle.h"


/* NAMESPACE */
namespace analysis
{

	class lhe : public particle
	{
	
	public:

		/* con & destructor */
		lhe() = default;
		lhe(double px, double py, double pz, double pe, double mass);
		//virtual ~lhe() = default;	

		/* copy & assignment */

	public:

		/* properties */
		int id() const;
		int type() const;
		bool is_final() const;

		/* kinematics */
		double px() const;
		double py() const;
		double pz() const;
		double pe() const;

		double pt() const;
		double eta() const;
		double phi() const;
		double mass() const;
		double y() const;

		/* input & output */
		void write(std::ostream& os) const;
		void read(std::istream& is);
		void write(std::ofstream& ofs) const;
		void read(std::ifstream& ifs);
		void write(ogzstream& ogzs) const;
		void read(igzstream& igzs);

	private:

		/* lhe data members */
		int p_id;
		int p_inout;
		int p_mother1;
		int p_mother2;
		int p_color1;
		int p_color2;
		double p_px;
		double p_py;
		double p_pz;
		double p_pe;
		double p_mass;
		double p_btag;
		double p_hel;

	};

/* NAMESPACE */
}

#endif
