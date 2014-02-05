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

#include "event/event.h"
#include "utility/utility.h"
#include "plot/plot.h"
#include "plot/plot2d.h"

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
	remove("output/test_plot_ht.png");
	remove("output/test_plot_met.png");
	remove("output/test_plot_2d_lhco.png");
	remove("output/test_plot_2d_lhe.png");
	
	// load the lhco and lhe events for this test
	vector<event*> events_lhco;
	vector<event*> events_lhe;
	read_lhco(events_lhco, "input/test_plot_events.lhco.gz");
	read_lhe(events_lhe, "input/test_plot_events.lhe.gz");
	
	// make a 1D ht plot
	plot ht("test_plot_ht", "output/");
	ht.add_sample(events_lhco, plot_ht, "LHCO");
	ht.add_sample(events_lhe, plot_ht, "LHE");
	ht.run();
	bool test_plot_ht_passed = is_regular_file("output/test_plot_ht.png");
	
	// make a 1D met plot
	plot met("test_plot_met", "output/");
	met.add_sample(events_lhco, plot_met, "LHCO");
	met.add_sample(events_lhe, plot_met, "LHE");
	met.run();
	bool test_plot_met_passed = is_regular_file("output/test_plot_met.png");
	
	// make a 2D ht vs met plot for LHCO
	plot2d lhco_2d("test_plot_2d_lhco", "output/");
	lhco_2d.set_functions(plot_ht, plot_met);
	lhco_2d.set_x_bins(25, 0, 1200);
	lhco_2d.set_y_bins(25, 0, 800);
	lhco_2d.add_sample(events_lhco, "LHCO");
	lhco_2d.run();
	bool test_plot_lhco_2d_passed = is_regular_file("output/test_plot_2d_lhco.png");
	
	// make a 2D ht vs met plot for LHCO
	plot2d lhe_2d("test_plot_2d_lhe", "output/");
	lhe_2d.set_functions(plot_ht, plot_met);
	lhe_2d.set_x_bins(25, 0, 1200);
	lhe_2d.set_y_bins(25, 0, 800);
	lhe_2d.add_sample(events_lhe, "LHE");
	lhe_2d.run();
	bool test_plot_lhe_2d_passed = is_regular_file("output/test_plot_2d_lhe.png");
	
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
