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
	remove("output/plot_test_htmet.png");
	
	// load the lhco events for this test
	vector<event*> events;
	read_lhco(events, "input/test_plot_events.lhco.gz");
	
	// make a 1D ht plot
	plot_ht ht;
	ht.set_name("test_ht");
	ht.set_folder("output/");
	ht.add_sample(events, "events");
	ht.run();
	bool test_plot_ht_passed = is_regular_file("output/plot_test_ht.png");
	
	// make a 1D met plot
	plot_met met;
	met.set_name("test_met");
	met.set_folder("output/");
	met.add_sample(events, "events");
	met.run();
	bool test_plot_met_passed = is_regular_file("output/plot_test_met.png");
	
	// make a 2D ht vs met plot
	plot2d_ht_met htmet;
	htmet.set_name("test_htmet");
	htmet.set_folder("output/");
	htmet.set_x_bins(25, 0, 1200);
	htmet.set_y_bins(25, 0, 800);
	htmet.add_sample(events, "events");
	htmet.run();
	bool test_plot_htmet_passed = is_regular_file("output/plot_test_htmet.png");
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Plot test: completed in " << duration << " seconds." << endl;
	cout << "HT plot creation has " << (test_plot_ht_passed ? "succeeded!" : "failed!") << endl; 
	cout << "Check: reasonable HT distribution." << endl;
	cout << "MET plot creation has " << (test_plot_met_passed ? "succeeded!" : "failed!") << endl; 
	cout << "Check: reasonable MET distribution." << endl;
	cout << "2D HT vs MET plot creation has " << (test_plot_htmet_passed ? "succeeded!" : "failed!") << endl;
	cout << "Check: previous distribution combined in 2D." << endl;
	cout << "=====================================================================" << endl;
	
	// return whether tests passed
	if (test_plot_ht_passed && test_plot_met_passed && test_plot_htmet_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;	
}
