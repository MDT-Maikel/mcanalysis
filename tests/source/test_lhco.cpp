/* LHCO Tests
 *
 * Test LHCO class by reading and writing events.
 * 
*/

#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>

#include <boost/filesystem.hpp>

#include "../../source/event/event.h"
#include "../../source/particle/particle.h"
#include "../../source/particle/lhco.h"
#include "../../source/utility/utility.h"

using namespace std;
using namespace boost::filesystem;
using namespace analysis;


// main program
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;

	// check if file exists, if so read the lhco events 
	bool test_reading_passed = true;
	vector<event*> events;
	if (is_regular_file("input/test_lhco_events.lhco.gz"))
	{
		read_lhco(events, "input/test_lhco_events.lhco.gz");
	}
	else
	{
		cout << "File (input/test_lhco_events.lhco.gz) for reading lhco events not available." << endl;
		test_reading_passed = false;
	}
	
	// check approximate transverse momentum balance of each event
	bool test_momentum_passed = true;
	if (test_reading_passed)
	{
		for (unsigned int i = 0; i < events.size(); i++)
		{
			event *ev = events[i];
			double px_sum = 0;
			double py_sum = 0;
			for (unsigned int j = 0; j < ev->size(); j++)
			{
				px_sum += (*ev)[j]->px();
				py_sum += (*ev)[j]->py();
			}
			// transverse momentum imbalance may at most be 200 GeV
			if (sqrt(px_sum * px_sum + py_sum * py_sum) > 200)
			{
				cout << "momentum imbalance found (px, py): (" << px_sum << ", " << py_sum << ")" << endl;
				test_momentum_passed = false;
			}
		}	
	}
		
	// if reading has succeeded write the events to file
	bool test_writing_passed = false;
	if (test_reading_passed)
	{
		write_lhco(events, "output/test_lhco_events.lhco.gz");
		
		// test whether the writing has succeeded
		if (is_regular_file("output/test_lhco_events.lhco.gz"))
		{
			// check file size of the newly created lhco file 
			// has to be within 90 to 110 percent of old file
			unsigned int in_size = file_size("input/test_lhco_events.lhco.gz");
			unsigned int out_size = file_size("output/test_lhco_events.lhco.gz");
			if (out_size < 11 * in_size / 10 && out_size > 9 * in_size / 10)
				test_writing_passed = true;
				
			// delete the created file as it is not needed for any checks
			remove("output/test_lhco_events.lhco.gz");
		}
	}
			
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "LHCO test: completed in " << duration << " seconds." << endl;
	cout << "LHCO reading has " << (test_reading_passed ? "succeeded!" : "failed!") << endl;
	cout << "LHCO writing has " << (test_writing_passed ? "succeeded!" : "failed!") << endl;
	cout << "LHCO momentum balance has " << (test_momentum_passed ? "succeeded!" : "failed!") << endl;
	cout << "=====================================================================" << endl;
	
	// return whether tests passed
	if (test_reading_passed && test_writing_passed && test_momentum_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}
