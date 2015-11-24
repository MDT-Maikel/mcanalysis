/* CMS paired dijet analysis
 *
 * Performs the CMS dijet analysis based on CMS-PAS-EXO-11-016,
 * these limits can be compared to the model-independent limits
 * from figure 3 in their conference note. This code calculates
 * the efficiency for the given events file.
 * 
*/

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "cuts/cuts.h"
#include "event/event.h"
#include "particle/particle.h"
#include "utility/utility.h"



using namespace std;
using namespace analysis;
using namespace boost;
using namespace boost::filesystem;

// cms 4 jet cut
class cut_cms_4jets: public cut
{
public:
	cut_cms_4jets() {}

	bool operator() (const event *ev) 
	{ 
		// get the fourth jet
		const particle *j4 = ev->get(ptype_jet, 4, 2.5);
		
		// reject event if jets not found
		if (!j4)
			return false;
		
		// reject event if pt <= 110 GeV
		if (j4->pt() <= 110)
			return false;

		// event passed
		return true;
	}
};

// cms dijet pairing
class cut_cms_dijetpair: public cut
{
public:
	cut_cms_dijetpair() {}

	bool operator() (const event *ev) 
	{ 
		// get the four jets
		const particle *j1 = ev->get(ptype_jet, 1, 2.5);
		const particle *j2 = ev->get(ptype_jet, 2, 2.5);
		const particle *j3 = ev->get(ptype_jet, 3, 2.5);
		const particle *j4 = ev->get(ptype_jet, 4, 2.5);

		bool deltar_pairing1 = min(delta_r(j1, j2), delta_r(j3, j4)) >= 0.7;
		bool deltar_pairing2 = min(delta_r(j1, j3), delta_r(j2, j4)) >= 0.7;
		bool deltar_pairing3 = min(delta_r(j1, j4), delta_r(j2, j3)) >= 0.7;

		double deltam_pairing1 = abs(mass({j1, j2}) - mass({j3, j4}));
		double deltam_pairing2 = abs(mass({j1, j3}) - mass({j2, j4}));
		double deltam_pairing3 = abs(mass({j1, j4}) - mass({j2, j3}));

		double mavg_pairing1 = 0.5 * (mass({j1, j2}) + mass({j3, j4}));
		double mavg_pairing2 = 0.5 * (mass({j1, j3}) + mass({j2, j4}));
		double mavg_pairing3 = 0.5 * (mass({j1, j4}) + mass({j2, j3}));

		double delta_pairing1 = min(j1->pt() + j2->pt(), j3->pt() + j4->pt()) - mavg_pairing1;
		double delta_pairing2 = min(j1->pt() + j3->pt(), j2->pt() + j4->pt()) - mavg_pairing2;
		double delta_pairing3 = min(j1->pt() + j4->pt(), j2->pt() + j3->pt()) - mavg_pairing3;

		if (deltar_pairing1 && deltam_pairing1 < 0.15 * mavg_pairing1 && (deltam_pairing1 < deltam_pairing2 || !deltar_pairing2) && (deltam_pairing1 < deltam_pairing3 || !deltar_pairing3))
		{
			if (delta_pairing1 > 25.0)
				return true;
		}
		else if (deltar_pairing2 && deltam_pairing2 < 0.15 * mavg_pairing2 && (deltam_pairing2 < deltam_pairing1 || !deltar_pairing1) && (deltam_pairing2 < deltam_pairing3 || !deltar_pairing3))
		{
			if (delta_pairing2 > 25.0)
				return true;
		}
		else if (deltar_pairing3 && deltam_pairing3 < 0.15 * mavg_pairing3 && (deltam_pairing3 < deltam_pairing1 || !deltar_pairing1) && (deltam_pairing3 < deltam_pairing2 || !deltar_pairing2))
		{
			if (delta_pairing3 > 25.0)
				return true;
		}
		// event not passed
		return false;
	}
};

// main program with one argument representing the input file 
int main(int argc, char* argv[])
{
	// make sure there is one arguments
	if (argc != 2) 
	{
		cout << "specify at least one argument: <input.mc.gz>" << endl;
		return EXIT_FAILURE;
  	}

	// convert the arguments to the input file and output directory
	string input_file = argv[1];
	
	// load the events dependent on whether they are .lhe.gz or .lhco.gz
	vector<event*> events; 	
	read_lhco(events, input_file);
	
	// initialize cutlist and add all cms dijet cuts
	cuts dijet;
	cut_cms_4jets *dijet_cut = new cut_cms_4jets();
	dijet.add_cut(dijet_cut, "cms 4jet cut");
	cut_cms_dijetpair *dijet_cut_pair = new cut_cms_dijetpair();
	dijet.add_cut(dijet_cut_pair, "cms dijet pair cut");
	dijet.apply(events);
	//dijet.write(cout);
	double acceptance = dijet.efficiency();
	//cout << "acceptance: " << setprecision(6) << acceptance << endl;
	
	// delete all the cut pointers
	delete dijet_cut;
	delete dijet_cut_pair;
		
	// clear remaining event pointers
	delete_events(events);
	
	// just output the acceptance
	cout << acceptance << endl;
	
	// finished the analysis
	return EXIT_SUCCESS;	
}
