/* LHE particle class
 *
 *
*/

#include "lhe.h"


/* NAMESPACE */
namespace analysis
{

	/* properties */

	int lhe::id() const { return p_id; }
	int lhe::type() const { return 0; }

	/* kinematics */

	double lhe::px() const { return p_px; }
	double lhe::py() const { return p_py; }
	double lhe::pz() const { return p_pz; }
	double lhe::pe() const { return p_pe; }

	double lhe::pt() const { return p_px * p_px + p_py * p_py; }
	double lhe::eta() const { return -std::log(std::abs(std::tan(std::acos(p_pz / std::sqrt(p_px * p_px + p_py * p_py + p_pz * p_pz)) / 2))); }
	double lhe::phi() const 
	{ 
		if (p_px >= 0 && p_py >= 0)
			return std::atan(p_py / p_px);
		if (p_px >= 0 && p_py < 0)
			return std::atan(p_py / p_px) + 8 * std::atan(1);
		return std::atan(p_py / p_px) + 4 * std::atan(1);
	}
	double lhe::mass() const { return p_mass; }
	double lhe::y() const { return 0.5 * std::log((pe() + pz()) / (pe() - pz())); }

	/* input & output */

	void lhe::write(std::ostream& os) const
	{
		os << p_id << "\t" << p_inout << "\t" << p_mother1 << "\t" << p_mother2 << "\t" << p_color1 << "\t" << p_color2 << "\t" << std::fixed;
		os << std::setprecision(7) << p_px << "\t" << p_py << "\t" << p_pz << "\t" << p_pe << "\t" << p_mass << "\t";
		os << std::setprecision(1) << p_btag << "\t" << p_hel << std::endl;
	}

	void lhe::read(std::istream& is)
	{
		is >> p_id >> p_inout >> p_mother1 >> p_mother2 >> p_color1 >> p_color2 >> p_px >> p_py >> p_pz >> p_pe >> p_mass >> p_btag >> p_hel;
	}

	void lhe::write(std::ofstream& ofs) const
	{
		ofs << p_id << "\t" << p_inout << "\t" << p_mother1 << "\t" << p_mother2 << "\t" << p_color1 << "\t" << p_color2 << "\t" << std::fixed;
		ofs << std::setprecision(7) << p_px << "\t" << p_py << "\t" << p_pz << "\t" << p_pe << "\t" << p_mass << "\t";
		ofs << std::setprecision(1) << p_btag << "\t" << p_hel << std::endl;
	}

	void lhe::read(std::ifstream& ifs) 
	{
		ifs >> p_id >> p_inout >> p_mother1 >> p_mother2 >> p_color1 >> p_color2 >> p_px >> p_py >> p_pz >> p_pe >> p_mass >> p_btag >> p_hel;
	}

	void lhe::write(ogzstream& ogzs) const
	{
		ogzs << p_id << "\t" << p_inout << "\t" << p_mother1 << "\t" << p_mother2 << "\t" << p_color1 << "\t" << p_color2 << "\t" << std::fixed;
		ogzs << std::setprecision(7) << p_px << "\t" << p_py << "\t" << p_pz << "\t" << p_pe << "\t" << p_mass << "\t";
		ogzs << std::setprecision(1) << p_btag << "\t" << p_hel << std::endl;
	}

	void lhe::read(igzstream& igzs)
	{
		igzs >> p_id >> p_inout >> p_mother1 >> p_mother2 >> p_color1 >> p_color2 >> p_px >> p_py >> p_pz >> p_pe >> p_mass >> p_btag >> p_hel;
	}

/* NAMESPACE */
}
