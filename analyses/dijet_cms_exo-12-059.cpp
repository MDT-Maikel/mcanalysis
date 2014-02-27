/* CMS dijet analysis
 *
 * Performs the CMS dijet analysis based on CMS-PAS-EXO-12-059,
 * these limits can be compared to the model-independent limits
 * from figure 4 in their conference note. This code calculates
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

// atlas dijet cut
class cut_cms_dijet: public cut
{
public:
	cut_cms_dijet() {}

	bool operator() (const event *ev) 
	{ 
		// get first two jets
		const particle *j1 = ev->get(ptype_jet, 1);
		const particle *j2 = ev->get(ptype_jet, 2);
		
		// reject event if jets not found
		if (!j1 || !j2)
			return false;
		
		// reject event if pt <= 30 GeV or |eta| >= 2.5 for either jet
		if (j1->pt() <= 30 || j2->pt() <= 30)
			return false;
		if (abs(j1->eta()) >= 2.5 || abs(j2->eta()) >= 2.5)
			return false;
			
		// init wide jets as four vectors
		double wpx1 = j1->px(), wpy1 = j1->py(), wpz1 = j1->pz(), wpe1 = j1->pe(); 
		double wpx2 = j2->px(), wpy2 = j2->py(), wpz2 = j2->pz(), wpe2 = j2->pe(); 
		// construct wide jets by adding all other jets (pt > 30 GeV and |eta| < 2.5)
		unsigned int i = 1;
		const particle *jet;
		do
		{
			// find all jets with pt > 30 GeV and |eta| < 2.5
			jet = ev->get(ptype_jet, i, 2.5);
			// increase index
			i++;
			
			// check jet
			if (!jet)
				continue;			
			if (jet->pt() <= 30)
				continue;
			// don't accept the jet if it is the first or the second
			if (jet == j1 || jet == j2)
				continue;

			// calculate delta R1 and delta R2
			double dr1 = delta_r(jet, j1);
			double dr2 = delta_r(jet, j2);
			
			// add the jet to one of the wide jets
			if (dr1 < 1.1 && dr1 < dr2)
			{
				wpx1 += jet->px();
				wpy1 += jet->py();
				wpz1 += jet->pz();
				wpe1 += jet->pe();
			}
			if (dr2 < 1.1 && dr2 < dr1)
			{
				wpx2 += jet->px();
				wpy2 += jet->py();
				wpz2 += jet->pz();
				wpe2 += jet->pe();
			}
			
		} while (jet);
		
		const lhe *wj1 = new lhe(wpx1, wpy1, wpz1, wpe1);
		const lhe *wj2 = new lhe(wpx2, wpy2, wpz2, wpe2);

		// reject event if wide jets have |delta eta| >= 1.3
		if (delta_eta(wj1, wj2) >= 1.3)
			return false;
		
		// reject event if wide jets have |eta| >= 2.5
		if (abs(wj1->eta()) >= 2.5 || abs(wj2->eta()) >= 2.5)
			return false;
		
		// reject event if mjj <= 890 GeV 
		if (mass({wj1, wj2}) <= 890)
			return false;
			
		// trigger: reject event if HT < 650 GeV and mjj < 750 GeV
		// TODO
		
		// event passed
		return true;
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
	cut_cms_dijet *dijet_cut = new cut_cms_dijet();
	dijet.add_cut(dijet_cut, "cms dijet cuts");
	dijet.apply(events);
	//dijet.write(cout);
	double acceptance = dijet.efficiency();
	//cout << "acceptance: " << setprecision(6) << acceptance << endl;
	
	// delete all the cut pointers
	delete dijet_cut;
		
	// clear remaining event pointers
	delete_events(events);
	
	// just put out the acceptance
	cout << acceptance << endl;
	
	// finished the analysis
	return EXIT_SUCCESS;	
}
