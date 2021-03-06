/* LHE Tests
 *
 * Test LHE class by reading and writing events.
 * 
*/

#include <cmath>
#include <ctime>
#include <iostream>
#include <random>
#include <vector>

#include <boost/filesystem.hpp>

#include "event/event.h"
#include "particle/particle.h"
#include "particle/lhe.h"
#include "utility/utility.h"

using namespace std;
using namespace boost::filesystem;
using namespace analysis;


// main program
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;
	
	// keep track of success
	bool test_reading_passed = true;
		
	// check if file exists, if so read the lhco events 
	vector<event*> events;
	if (is_regular_file("../../files/tests/input/test_lhe_events.lhe.gz"))
	{
		read_lhe(events, "../../files/tests/input/test_lhe_events.lhe.gz");
	}
	else
	{
		cout << "File (../../files/tests/input/test_lhe_events.lhe.gz) for reading lhe events not available." << endl;
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
				// only take the particles in the final state
				if ((*ev)[j]->is_final())
				{				
					px_sum += (*ev)[j]->px();
					py_sum += (*ev)[j]->py();
				}
			}
			// transverse momentum imbalance may at most be 1 GeV
			if (sqrt(px_sum * px_sum + py_sum * py_sum) > 1)
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
		write_lhe(events, "../../files/tests/output/test_lhe_events.lhe.gz");
		
		// test whether the writing has succeeded		
		if (is_regular_file("../../files/tests/output/test_lhe_events.lhe.gz"))
		{
			// check file size of the newly created lhco file 
			// has to be within 70 to 110 percent of old file
			unsigned int in_size = file_size("../../files/tests/input/test_lhe_events.lhe.gz");
			unsigned int out_size = file_size("../../files/tests/output/test_lhe_events.lhe.gz");
			if (out_size < 11 * in_size / 10 && out_size > 7 * in_size / 10)
				test_writing_passed = true;
				
			// delete the created file as it is not needed for any checks
			remove("../../files/tests/output/test_lhe_events.lhe.gz");
		}
	}
			
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "LHE test: completed in " << duration << " seconds." << endl;
	cout << "LHE reading has " << (test_reading_passed ? "succeeded!" : "failed!") << endl;
	cout << "LHE writing has " << (test_writing_passed ? "succeeded!" : "failed!") << endl;
	cout << "LHE momentum balance has " << (test_momentum_passed ? "succeeded!" : "failed!") << endl;
	cout << "=====================================================================" << endl;
	
	// clear remaining event pointers
	delete_events(events);
	
	// return whether tests passed
	if (test_reading_passed && test_writing_passed && test_momentum_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}
