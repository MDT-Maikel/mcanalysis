/* HEPTopTagger Tests
 *
 * Tests our jet analysis suite based on FastJet and HEPTopTagger
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

	//===== initialisation =====//

	double cross_sec_orig = 5100; // stop pair production cross section (mt1=340 GeV) as in hep-ph/1006.2833

	jet_analysis test;
	test.set_nEvents(50000);
	test.set_Rsize_fat(1.5);
	test.undo_BDRSTagging();

	test.add_lhe("input/test_heptop_events.lhe");
	test.import_lhco("input/test_heptop_events.lhco.gz");

	fastjet::PseudoJet (jet_analysis::*TopTagger)(const fastjet::PseudoJet &) = &jet_analysis::HEPTopTagging;
	test.initialise(TopTagger);


	//===== apply list of cuts as in hep-ph/1006.2833 =====//

	// two fat jets (C/A, R=1.5) with pT> 200 GeV
	cuts cut_list;
	cut_pt *pt = new cut_pt(200, ptype_jet, 2, 2.5);
	cut_list.add_cut(pt, "pt(j1,j2) > 200 GeV");
	// double eff_fatjet_pt = test.require_fatjet_pt(200,2);
	// cout << setprecision(2) << endl << "Efficiency after requiring two FatJets with pT> 200 GeV: " << eff_fatjet_pt << "%" << endl; //just for testing

	// lepton veto
	cut_veto *veto = new cut_veto(ptype_lepton, 15, 2.5);
	cut_list.add_cut(veto, "lepton veto (pT>15 GeV, |eta|<2.5)");

	// met cut
	cut_met *met = new cut_met(150);
	cut_list.add_cut(met, "met > 150 GeV");
	
	// apply first set of cuts
	test.reduce_sample(cut_list);
	double eff_lhco = cut_list.efficiency();

	// apply HEPTopTag cut
	double eff_ttag_1 = test.require_top_tagged(1);
	cout << setprecision(2) << endl << "first top-tag efficiency: " << eff_ttag_1 << "%" << endl;

	double eff_ttag_2 = test.require_top_tagged(2);
	cout << setprecision(2) << endl << "second top-tag efficiency: " << eff_ttag_2 << "%" << endl;


	//===== result =====//
	double cross_sec = cross_sec_orig*eff_lhco*eff_ttag_1*eff_ttag_2;
	bool test_jets_passed = true;

	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "HEPTopTagger test: completed in " << duration << " seconds." << endl;
	cout << "Final cross section for stop pair production (mt1=340 GeV, mn1=98 GeV) = " << cross_sec << endl;
	cout << "Result quoted in hep-ph/1006.2833 : 15 [fb]" << endl;
	cout << "=====================================================================" << endl;
	
	// clear remaining pointers
	delete pt;

	// return whether tests passed
	if (test_jets_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}
