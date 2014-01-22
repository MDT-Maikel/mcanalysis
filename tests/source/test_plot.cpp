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
	clock_t clock_old;
	clock_old = clock();
	double duration;
	
	// keep track of success
	bool test_plot_passed = true;
		
	// load the lhco events for this test
	vector<event*> events;
	read_lhco(events, "input/test_plot_events.lhco.gz");
	
	// make a 1D ht plot
	plot_ht ht;
	ht.set_name("test_ht");
	ht.set_folder("output/");
	ht.add_sample(events, "events");
	ht.run();
	test_plot_passed = test_plot_passed && is_regular_file("output/plot_test_ht.png");
	
	// make a 1D met plot
	plot_met met;
	met.set_name("test_met");
	met.set_folder("output/");
	met.add_sample(events, "events");
	met.run();
	test_plot_passed = test_plot_passed && is_regular_file("output/plot_test_met.png");
	
	// make a 2D ht vs met plot
	plot2d_ht_met ht_met;
	ht_met.set_name("test_ht_met");
	ht_met.set_folder("output/");
	ht_met.set_x_bins(25, 0, 1200);
	ht_met.set_y_bins(25, 0, 800);
	ht_met.add_sample(events, "events");
	ht_met.run();
	test_plot_passed = test_plot_passed && is_regular_file("output/plot_test_ht_met.png");
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Plot test: completed in " << duration << " seconds." << endl;
	cout << "Plotting 1D and 2D plots has " << (test_plot_passed ? "succeeded!" : "failed!") << endl;
	cout << "=====================================================================" << endl;
	
	// return whether tests passed
	if (test_plot_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;	
}
