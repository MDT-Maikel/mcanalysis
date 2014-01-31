/* LHE particle class
 *
 *
*/

#include "lhe.h"


/* NAMESPACE */
namespace analysis
{
	
	/* con & destructor */
	
	lhe::lhe(double px, double py, double pz, double pe, double mass)
	{
		p_id = 0;
		p_px = px;
		p_py = py;
		p_pz = pz;
		p_pe = pe;
		p_mass = mass;
	}

	/* properties: type */

	int lhe::id() const { return p_id; }
	
	unsigned int lhe::type() const 
	{ 
		switch (p_id)
		{		
		// photon
		case pid_photon:
			return ptype_photon;
		// electron		
		case pid_electron: case pid_positron:
			return ptype_electron;
		// muon
		case pid_muon: case pid_muonbar:
			return ptype_muon;
		// tau
		case pid_tau: case pid_taubar:
			return ptype_tau;
		// jet
		case pid_gluon:
		case pid_down: case pid_up: case pid_strange: case pid_charm: case pid_bottom: case pid_top:
		case pid_downbar: case pid_upbar: case pid_strangebar: case pid_charmbar: case pid_bottombar: case pid_topbar:
			return ptype_jet;	
		// missing energy
		case pid_nue: case pid_num: case pid_nut: case pid_nuebar: case pid_numbar: case pid_nutbar:
		case 8880022:
			return ptype_met;
		// unknown
		default:
			return ptype_none;	
		}		
		// unknown
		return ptype_none; 
	}
	
	/* properties: state */
	
	void lhe::set_final(bool is_final)
	{
		p_inout = is_final ? 1 : 0;
	}
	
	bool lhe::is_final() const 
	{ 
		return p_inout == 1; 
	}
	
	/* properties: quantum numbers */
	
	double lhe::charge() const
	{
		switch (p_id)
		{	
		case pid_electron: case pid_muon: case pid_tau:
			return -1.0;
		case pid_positron: case pid_muonbar: case pid_taubar:
			return 1.0;
		case pid_down: case pid_strange: case pid_bottom:
			return -1.0 / 3.0;
		case pid_up: case pid_charm: case pid_top:
			return 2.0 / 3.0;
		case pid_downbar: case pid_strangebar: case pid_bottombar:
			return 1.0 / 3.0;
		case pid_upbar: case pid_charmbar: case pid_topbar:
			return -2.0 / 3.0;
		case pid_w:
			return 1.0; // TODO: differentiate between Wplus and Wminus
		case pid_nue: case pid_num: case pid_nut: case pid_nuebar: case pid_numbar: case pid_nutbar:
		case pid_gluon: case pid_photon: case pid_z:
			return 0.0;				
		default:
			return 0.0;		
		}
		return 0.0;
	}

	/* properties: kinematics */

	double lhe::px() const { return p_px; }
	double lhe::py() const { return p_py; }
	double lhe::pz() const { return p_pz; }
	double lhe::pe() const { return p_pe; }

	double lhe::pt() const { return std::sqrt(p_px * p_px + p_py * p_py); }
	double lhe::eta() const { return -std::log(std::abs(std::tan(std::acos(p_pz / std::sqrt(p_px * p_px + p_py * p_py + p_pz * p_pz)) / 2))); }
	double lhe::phi() const 
	{ 
		if (p_px < 0 && p_py >= 0)
			return 4 * std::atan(1) + std::atan(p_py / p_px);
		if (p_px < 0 && p_py < 0)
			return std::atan(p_py / p_px) - 4 * std::atan(1);
		return std::atan(p_py / p_px);
	}
	double lhe::mass() const { return p_mass; }
	
	double lhe::y() const { return 0.5 * std::log((pe() + pz()) / (pe() - pz())); }

	/* input & output */

	void lhe::write(std::ostream& os) const
	{
		os << p_id << "\t" << p_inout << "\t" << p_mother1 << "\t" << p_mother2 << "\t" << p_color1 << "\t" << p_color2 << "\t" << std::fixed;
		os << std::scientific << std::setprecision(7) << p_px << "\t" << p_py << "\t" << p_pz << "\t" << p_pe << "\t" << p_mass << "\t";
		os << std::fixed << std::setprecision(1) << p_btag << "\t" << p_hel << std::endl;
	}

	void lhe::read(std::istream& is)
	{
		std::string str_px, str_py, str_pz, str_pe, str_mass;
		is >> p_id >> p_inout >> p_mother1 >> p_mother2 >> p_color1 >> p_color2 >> str_px >> str_py >> str_pz >> str_pe >> str_mass >> p_btag >> p_hel;
		p_px = boost::lexical_cast<double>(str_px);
		p_py = boost::lexical_cast<double>(str_py);
		p_pz = boost::lexical_cast<double>(str_pz);
		p_pe = boost::lexical_cast<double>(str_pe);
		p_mass = boost::lexical_cast<double>(str_mass);	
	}

	void lhe::write(std::ofstream& ofs) const
	{
		ofs << p_id << "\t" << p_inout << "\t" << p_mother1 << "\t" << p_mother2 << "\t" << p_color1 << "\t" << p_color2 << "\t" << std::fixed;
		ofs << std::scientific << std::setprecision(7) << p_px << "\t" << p_py << "\t" << p_pz << "\t" << p_pe << "\t" << p_mass << "\t";
		ofs << std::fixed << std::setprecision(1) << p_btag << "\t" << p_hel << std::endl;
	}

	void lhe::read(std::ifstream& ifs) 
	{
		std::string str_px, str_py, str_pz, str_pe, str_mass;
		ifs >> p_id >> p_inout >> p_mother1 >> p_mother2 >> p_color1 >> p_color2 >> str_px >> str_py >> str_pz >> str_pe >> str_mass >> p_btag >> p_hel;
		p_px = boost::lexical_cast<double>(str_px);
		p_py = boost::lexical_cast<double>(str_py);
		p_pz = boost::lexical_cast<double>(str_pz);
		p_pe = boost::lexical_cast<double>(str_pe);
		p_mass = boost::lexical_cast<double>(str_mass);	
	}

	void lhe::write(ogzstream& ogzs) const
	{
		ogzs << p_id << "\t" << p_inout << "\t" << p_mother1 << "\t" << p_mother2 << "\t" << p_color1 << "\t" << p_color2 << "\t" << std::fixed;
		ogzs << std::scientific << std::setprecision(7) << p_px << "\t" << p_py << "\t" << p_pz << "\t" << p_pe << "\t" << p_mass << "\t";
		ogzs << std::fixed << std::setprecision(1) << p_btag << "\t" << p_hel << std::endl;
	}

	void lhe::read(igzstream& igzs)
	{
		std::string str_px, str_py, str_pz, str_pe, str_mass;
		igzs >> p_id >> p_inout >> p_mother1 >> p_mother2 >> p_color1 >> p_color2 >> str_px >> str_py >> str_pz >> str_pe >> str_mass >> p_btag >> p_hel;
		p_px = boost::lexical_cast<double>(str_px);
		p_py = boost::lexical_cast<double>(str_py);
		p_pz = boost::lexical_cast<double>(str_pz);
		p_pe = boost::lexical_cast<double>(str_pe);
		p_mass = boost::lexical_cast<double>(str_mass);		
	}

/* NAMESPACE */
}
