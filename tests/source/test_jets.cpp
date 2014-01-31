/* FastJet Tests
 *
 * Tests our jet analysis suite based on FastJet
 * 
*/

#include <iostream>
#include <iomanip>  
#include <cmath>
#include <ctime>

#include "Pythia8/Pythia.h"
#include "Pythia8/FastJet3.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JHTopTagger.hh"

#include "../../source/jet_analysis/jet_analysis.h"

using namespace std;
using namespace analysis;
using namespace Pythia8;
using namespace boost::filesystem;


// main program
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;

	// initialise jet_analysis class
	jet_analysis test;
	int nr_events = 100;
	test.set_nEvents(nr_events);
	test.undo_BDRSTagging();

	// Set lhe input files
	test.add_lhe("input/test_fastjet_events.lhe");

	// initialisation (possibly specifying the TopTagger)
	// test.initialise(); // notice: default TopTagger is JHTopTagging, not needed to be specified
	fastjet::PseudoJet (jet_analysis::*TopTagger)(const fastjet::PseudoJet &) = &jet_analysis::HEPTopTagging;
	test.initialise(TopTagger);

	// apply jet pt cut
	cuts cut_list;
	cut_pt *pt; pt = new cut_pt(200, particle::type_jet, 1, 5.0);
	cut_list.add_cut(pt, "pt(j) > 200 GeV");
	test.reduce_sample(cut_list);

	// extract tgging information
	int ntops = 1;
	double eff = test.require_top_tagged(ntops);
	cout << setprecision(2) << endl << "Efficiency of tagging requirement: " << eff << endl;

	// determine success
	bool test_jets_passed = true;

	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Jet analysis test: completed in " << duration << " seconds." << endl;
	cout << "Jet tests have " << (test_jets_passed ? "passed!" : "failed!") << endl;
	cout << "=====================================================================" << endl;
	
	// clear remaining pointers
	delete pt;

	// return whether tests passed
	if (test_jets_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}
