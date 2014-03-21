/* FastJet Merging Tests
 *
 * Tests our jet analysis suite and merging procedure based on FastJet and Pythia8
 * 
*/

#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>

#include "Pythia8/Pythia.h"
#include "Pythia8/FastJet3.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JHTopTagger.hh"

#include "jet_analysis/jet_analysis.h"

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
	test.set_nEvents(1000);
	test.undo_TopTagging();
	test.undo_BDRSTagging();
	// test.set_fast_showering();

	// Set lhe input file and merging (njets merged)
	test.import_lhe("../../files/tests/input/w_production");
	test.set_merging_process("pp>LEPTONS,NEUTRINOS");
	test.set_merging_njmax(2);
	test.set_merging_scale(15);

	// initialisation (without TopTagger specification)
	test.initialise();

	// apply jet pt cut
	cuts cut_list;
	cut_pt *pt; pt = new cut_pt(20, ptype_jet, 1, 5.0);
	cut_list.add_cut(pt, "pt(j) > 20 GeV");
	test.reduce_sample(cut_list);

	// determine success
	bool test_jets_passed = true;

	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Merging test: completed in " << duration << " seconds." << endl;
	cout << "Merging test have " << (test_jets_passed ? "passed!" : "failed!") << endl;
	cout << "=====================================================================" << endl;
	
	// clear remaining pointers
	delete pt;

	// return whether tests passed
	if (test_jets_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}
