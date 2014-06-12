/* LHCO particle class
 *
 *
*/

#include "lhco.h"


/* NAMESPACE */
namespace analysis
{

	/* con & destructor */

	lhco::lhco(unsigned int type, double eta, double phi, double pt, double jmass, double ntrk, double btag, double hadem, double dum1, double dum2)
	{
		p_type = type;
		p_eta = eta;
		p_phi = phi;
		p_pt = pt;
		p_jmass = jmass;
		p_ntrk = ntrk;
		p_btag = btag;
		p_hadem = hadem;
		p_dum1 = dum1;
		p_dum2 = dum2;
	}
	
	/* copy & assignment */
			
	lhco* lhco::clone() const
	{
		return new lhco(*this);
	}

	/* properties: type */

	int lhco::id() const { return 0; }
	unsigned int lhco::type() const { return p_type; }
	
	/* properties: state */
	
	void lhco::set_final(bool is_final) { /* to be overloaded */ };
	bool lhco::is_final() const { return true; }
	
	/* properties: quantum numbers */
	
	double lhco::charge() const
	{
		// TODO: for other types and tau -1 vs -3
		if (p_type & ptype_leptonall)
			return p_ntrk;
		return 0.0;
	}

	double lhco::bjet() const
	{
		if (p_type & ptype_jet)
			return p_btag;
		return 0.0;
	}

	/* properties: kinematics */

	double lhco::pt() const { return p_pt; }
	double lhco::eta() const { return p_eta; }
	double lhco::phi() const { return p_phi; }
	double lhco::mass() const { return p_jmass; }

	double lhco::px() const { return pt() * std::cos(phi()); }
	double lhco::py() const { return pt() * std::sin(phi()); }
	double lhco::pz() const { return pt() * std::sinh(eta()); }
	double lhco::pe() const { return std::sqrt(pow(pt() * std::cosh(eta()), 2) + pow(mass(), 2)); }
	
	double lhco::y() const { return 0.5 * std::log((pe() + pz()) / (pe() - pz())); }

	/* input & output */

	int lhco::type_to_int(unsigned int type) const
	{
		switch (type)
		{
			case ptype_photon:
				return 0;
			case ptype_electron:
				return 1;
			case ptype_muon:
				return 2;
			case ptype_tau:
				return 3;
			case ptype_jet:
				return 4;
			case ptype_met:
				return 6;
			default:
				return -1;		
		}
		return -1;
	}

	void lhco::write(std::ostream& os) const
	{
		os << type_to_int(p_type) << "\t" << std::fixed; 
		os << std::setprecision(3) << p_eta << "\t" << p_phi << "\t"; 
		os << std::setprecision(2) << p_pt << "\t" << p_jmass << "\t"; 
		os << std::setprecision(1) << p_ntrk << "\t" << p_btag << "\t"; 
		os << std::setprecision(2) << p_hadem << "\t"; 
		os << std::setprecision(1) << p_dum1 << "\t" << p_dum2 << std::endl; 
	}

	void lhco::read(std::istream& is)
	{
		unsigned int type;
		is >> type >> p_eta >> p_phi >> p_pt >> p_jmass >> p_ntrk >> p_btag >> p_hadem >> p_dum1 >> p_dum2;
		p_type = 1 << type;
	}

	void lhco::write(std::ofstream& ofs) const
	{
		ofs << type_to_int(p_type) << "\t" << std::fixed; 
		ofs << std::setprecision(3) << p_eta << "\t" << p_phi << "\t";
		ofs << std::setprecision(2) << p_pt << "\t" << p_jmass << "\t";
		ofs << std::setprecision(1) << p_ntrk << "\t" << p_btag << "\t";
		ofs << std::setprecision(2) << p_hadem << "\t";
		ofs << std::setprecision(1) << p_dum1 << "\t" << p_dum2 << std::endl;
	}

	void lhco::read(std::ifstream& ifs)
	{
		unsigned int type;
		ifs >> type >> p_eta >> p_phi >> p_pt >> p_jmass >> p_ntrk >> p_btag >> p_hadem >> p_dum1 >> p_dum2;
		p_type = 1 << type;
	}

	void lhco::write(ogzstream& ogzs) const
	{
		ogzs << type_to_int(p_type) << "\t" << std::fixed;
		ogzs << std::setprecision(3) << p_eta << "\t" << p_phi << "\t";
		ogzs << std::setprecision(2) << p_pt << "\t" << p_jmass << "\t";
		ogzs << std::setprecision(1) << p_ntrk << "\t" << p_btag << "\t";
		ogzs << std::setprecision(2) << p_hadem << "\t";
		ogzs << std::setprecision(1) << p_dum1 << "\t" << p_dum2 << std::endl;
	}

	void lhco::read(igzstream& igzs)
	{
		unsigned int type;
		igzs >> type >> p_eta >> p_phi >> p_pt >> p_jmass >> p_ntrk >> p_btag >> p_hadem >> p_dum1 >> p_dum2;
		p_type = 1 << type;
	}

/* NAMESPACE */
}
