/* LHCO particle class
 *
 *
*/

#ifndef INC_LHCO
#define INC_LHCO

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../../deps/gzstream/gzstream.h"

#include "particle.h"


/* NAMESPACE */
namespace analysis
{

	class lhco : public particle
	{
	
	public:

		/* con & destructor */
		lhco(unsigned int type = ptype_none, double eta = 0, double phi = 0, double pt = 0, double jmass = 0, double ntrk = 0, double btag = 0, double hadem = 0, double dum1 = 0, double dum2 = 0);
		//~lhco() = default;

		/* copy & assignment */
		lhco(const lhco&) = default;
    	lhco& operator = (const lhco&) = default;

		/* properties: type */
		int id() const;
		unsigned int type() const;
		
		/* properties: state */
		void set_final(bool is_final);
		bool is_final() const;
		
		/* properties: quantum numbers */
		double charge() const;

		/* properties: kinematics */
		double pt() const;
		double eta() const;
		double phi() const;
		double mass() const;

		double px() const;
		double py() const;
		double pz() const;
		double pe() const;
		
		double y() const;

		/* input & output */
		int type_to_int(unsigned int type) const;
		void write(std::ostream& os) const;
		void read(std::istream& is);
		void write(std::ofstream& ofs) const;
		void read(std::ifstream& ifs);
		void write(ogzstream& ogzs) const;
		void read(igzstream& igzs);

	private:

		/* lhco data members */
		unsigned int p_type;
		double p_eta;
		double p_phi;
		double p_pt;
		double p_jmass;
		double p_ntrk;
		double p_btag;
		double p_hadem;
		double p_dum1;
		double p_dum2;

	};

/* NAMESPACE */
}

#endif
