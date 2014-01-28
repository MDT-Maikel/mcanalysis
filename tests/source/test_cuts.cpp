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
#include "../../source/cuts/cuts_standard.h"

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
	
	// initiate general cut class and the cuts
	double pt1_cut = 200;
	double pt2_cut = 100;
	double met_cut = 100;
	double ht_cut = 400;
	cuts test_cuts;
	
	// create pt cut on leading jet
	cut_pt *pt1 = new cut_pt(particle::type_jet, 1, pt1_cut, 2.5);
	test_cuts.add_cut(pt1);
	
	// create pt cut on second jet
	cut_pt *pt2 = new cut_pt(particle::type_jet, 2, pt2_cut, 2.5);
	test_cuts.add_cut(pt2);
	
	// create met cut
	cut_met *met = new cut_met(met_cut);
	test_cuts.add_cut(met);
	
	// create ht cut
	cut_ht *ht = new cut_ht(particle::type_jet, 30, 2.5, ht_cut);
	test_cuts.add_cut(ht);
	
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
	bool test_cuts_passed = eff_lhe / eff_lhco > 0.9 && eff_lhe / eff_lhco < 1.1;
		
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
