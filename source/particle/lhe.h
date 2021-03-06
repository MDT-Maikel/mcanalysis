/* LHE particle class
 *
 *
*/

#ifndef INC_LHE
#define INC_LHE

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
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
		lhe(double px, double py, double pz, double pe, double mass = 0);
		//virtual ~lhe() = default;	

		/* copy & assignment */
		// lhe(const lhe&) = default;
    	// lhe& operator = (const lhe&) = default;
    	virtual lhe* clone() const;

		/* properties: type */
		int id() const;
		unsigned int type() const;
		
		/* properties: state */
		void set_final(bool is_final);
		bool is_final() const;
		
		/* properties: quantum numbers */
		double charge() const;
		double bjet() const;

		/* properties: kinematics */
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
