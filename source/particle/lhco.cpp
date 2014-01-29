/* LHCO particle class
 *
 *
*/

#include "lhco.h"


/* NAMESPACE */
namespace analysis
{

	/* con & destructor */

	lhco::lhco(int type, double eta, double phi, double pt, double jmass, double ntrk, double btag, double hadem, double dum1, double dum2)
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

	/* properties */

	int lhco::id() const { return 0; }
	int lhco::type() const { return p_type; }
	void lhco::set_final(bool is_final) { /* to be overloaded */ };
	bool lhco::is_final() const { return true; }

	/* kinematics */

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

	void lhco::write(std::ostream& os) const
	{
		os << p_type << "\t" << std::fixed; 
		os << std::setprecision(3) << p_eta << "\t" << p_phi << "\t"; 
		os << std::setprecision(2) << p_pt << "\t" << p_jmass << "\t"; 
		os << std::setprecision(1) << p_ntrk << "\t" << p_btag << "\t"; 
		os << std::setprecision(2) << p_hadem << "\t"; 
		os << std::setprecision(1) << p_dum1 << "\t" << p_dum2 << std::endl; 
	}

	void lhco::read(std::istream& is)
	{
		is >> p_type >> p_eta >> p_phi >> p_pt >> p_jmass >> p_ntrk >> p_btag >> p_hadem >> p_dum1 >> p_dum2;
	}

	void lhco::write(std::ofstream& ofs) const
	{
		ofs << p_type << "\t" << std::fixed; 
		ofs << std::setprecision(3) << p_eta << "\t" << p_phi << "\t";
		ofs << std::setprecision(2) << p_pt << "\t" << p_jmass << "\t";
		ofs << std::setprecision(1) << p_ntrk << "\t" << p_btag << "\t";
		ofs << std::setprecision(2) << p_hadem << "\t";
		ofs << std::setprecision(1) << p_dum1 << "\t" << p_dum2 << std::endl;
	}

	void lhco::read(std::ifstream& ifs)
	{
		ifs >> p_type >> p_eta >> p_phi >> p_pt >> p_jmass >> p_ntrk >> p_btag >> p_hadem >> p_dum1 >> p_dum2;
	}

	void lhco::write(ogzstream& ogzs) const
	{
		ogzs << p_type << "\t" << std::fixed;
		ogzs << std::setprecision(3) << p_eta << "\t" << p_phi << "\t";
		ogzs << std::setprecision(2) << p_pt << "\t" << p_jmass << "\t";
		ogzs << std::setprecision(1) << p_ntrk << "\t" << p_btag << "\t";
		ogzs << std::setprecision(2) << p_hadem << "\t";
		ogzs << std::setprecision(1) << p_dum1 << "\t" << p_dum2 << std::endl;
	}

	void lhco::read(igzstream& igzs)
	{
		igzs >> p_type >> p_eta >> p_phi >> p_pt >> p_jmass >> p_ntrk >> p_btag >> p_hadem >> p_dum1 >> p_dum2;
	}

/* NAMESPACE */
}
