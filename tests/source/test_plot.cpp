/* Plot Tests
 *
 * Test plot and plot2d class.
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
#include "../../source/plot/plot.h"
#include "../../source/plot/plot_standard.h"
#include "../../source/plot/plot2d.h"

using namespace std;
using namespace boost::filesystem;
using namespace analysis;


// main program
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;
	
	// remove possible existing output files
	remove("output/plot_test_ht.png");
	remove("output/plot_test_met.png");
	remove("output/plot_test_2d_lhco.png");
	remove("output/plot_test_2d_lhe.png");
	
	// load the lhco and lhe events for this test
	vector<event*> events_lhco;
	vector<event*> events_lhe;
	read_lhco(events_lhco, "input/test_plot_events.lhco.gz");
	read_lhe(events_lhe, "input/test_plot_events.lhe.gz");
	
	// make a 1D ht plot
	plot_ht ht;
	ht.set_name("test_ht");
	ht.set_folder("output/");
	ht.add_sample(events_lhco, "LHCO");
	ht.add_sample(events_lhe, "LHE");
	ht.run();
	bool test_plot_ht_passed = is_regular_file("output/plot_test_ht.png");
	
	// make a 1D met plot
	plot_met met;
	met.set_name("test_met");
	met.set_folder("output/");
	met.add_sample(events_lhco, "LHCO");
	met.add_sample(events_lhe, "LHE");
	met.run();
	bool test_plot_met_passed = is_regular_file("output/plot_test_met.png");
	
	// make a 2D ht vs met plot for LHCO
	plot2d_ht_met lhco_2d;
	lhco_2d.set_name("test_2d_lhco");
	lhco_2d.set_folder("output/");
	lhco_2d.set_x_bins(25, 0, 1200);
	lhco_2d.set_y_bins(25, 0, 800);
	lhco_2d.add_sample(events_lhco, "events");
	lhco_2d.run();
	bool test_plot_lhco_2d_passed = is_regular_file("output/plot_test_2d_lhco.png");
	
	// make a 2D ht vs met plot for LHCO
	plot2d_ht_met lhe_2d;
	lhe_2d.set_name("test_2d_lhe");
	lhe_2d.set_folder("output/");
	lhe_2d.set_x_bins(25, 0, 1200);
	lhe_2d.set_y_bins(25, 0, 800);
	lhe_2d.add_sample(events_lhe, "events");
	lhe_2d.run();
	bool test_plot_lhe_2d_passed = is_regular_file("output/plot_test_2d_lhe.png");
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Plot test: completed in " << duration << " seconds." << endl;
	cout << "HT plot creation has " << (test_plot_ht_passed ? "succeeded!" : "failed!") << endl; 
	cout << "Check: reasonable HT distribution for LHCO and LHE." << endl;
	cout << "MET plot creation has " << (test_plot_met_passed ? "succeeded!" : "failed!") << endl; 
	cout << "Check: reasonable MET distribution for LHCO and LHE." << endl;
	cout << "2D HT vs MET plot (LHCO) creation has " << (test_plot_lhco_2d_passed ? "succeeded!" : "failed!") << endl;
	cout << "Check: previous distributions combined in 2D." << endl;
	cout << "2D HT vs MET plot (LHE) creation has " << (test_plot_lhe_2d_passed ? "succeeded!" : "failed!") << endl;
	cout << "Check: previous distributions combined in 2D." << endl;
	cout << "=====================================================================" << endl;
	
	// return whether tests passed
	if (test_plot_ht_passed && test_plot_met_passed && test_plot_lhco_2d_passed && test_plot_lhe_2d_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;	
}
