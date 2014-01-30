/* Cut Tests
 *
 * Test cuts classes.
 * 
*/

#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>

#include <boost/filesystem.hpp>

#include "../../source/event/event.h"
#include "../../source/utility/utility.h"
#include "../../source/cuts/cuts.h"

using namespace std;
using namespace boost::filesystem;
using namespace analysis;


// main program
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;
	
	// load the lhco and lhe events for this test
	vector<event*> events_lhco;
	vector<event*> events_lhe;
	read_lhco(events_lhco, "input/test_cuts_events.lhco.gz");
	read_lhe(events_lhe, "input/test_cuts_events.lhe.gz");
	
	// initiate general cut class
	cuts test_cuts;
	
	// create pt cut on leading jet
	cut_pt *pt1 = new cut_pt(200, particle::type_jet, 1, 2.5);
	test_cuts.add_cut(pt1, "pt(j1) > 200 GeV");
	
	// create pt cut on second jet
	cut_pt *pt2 = new cut_pt(200, particle::type_jet, 2, 2.5);
	test_cuts.add_cut(pt2, "pt(j2) > 200 GeV");
	
	// create met cut
	cut_met *met = new cut_met(100);
	test_cuts.add_cut(met, "met > 100 GeV");
	
	// create ht cut
	cut_ht *ht = new cut_ht(400, particle::type_jet, 20, 5.0);
	test_cuts.add_cut(ht, "ht > 400 GeV");
	
	// create veto cut
	cut_veto *veto = new cut_veto(particle::type_electron, 20, 2.5);
	test_cuts.add_cut(veto, "electron veto");
	
	// run the cuts on the LHCO sample
	test_cuts.apply(events_lhco);
	test_cuts.write(cout);
	double eff_lhco = test_cuts.efficiency();
	test_cuts.clear();
	
	// run the cuts on the LHE sample
	test_cuts.apply(events_lhe);
	test_cuts.write(cout);
	double eff_lhe = test_cuts.efficiency();
	test_cuts.clear();
	
	// determine success
	bool test_cuts_passed = eff_lhe / eff_lhco > 0.8 && eff_lhe / eff_lhco < 1.2;
		
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Cuts test: completed in " << duration << " seconds." << endl;
	cout << "Cuts for LHCO and LHE have " << (test_cuts_passed ? "passed!" : "failed!") << endl;
	cout << "=====================================================================" << endl;
	
	// return whether tests passed
	if (test_cuts_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;	
}
