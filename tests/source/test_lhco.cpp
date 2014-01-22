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
	clock_t clock_old;
	clock_old = clock();
	double duration;
	
	// keep track of success
	bool test_lhco_passed = true;
		
	// check if file exists, if so read the lhco events 
	vector<event*> events;
	if (is_regular_file("input/test_lhco_events.lhco.gz"))
	{
		read_lhco(events, "input/test_lhco_events.lhco.gz");
	}
	else
	{
		cout << "File (input/test_lhco_events.lhco.gz) for reading lhco events not available." << endl;
		test_lhco_passed = false;
	}
		
	// if reading has succeeded write the events to file
	if (test_lhco_passed)
	{
		write_lhco(events, "output/test_lhco_events.lhco.gz");
		
		// test whether the writing has succeeded
		test_lhco_passed = false;
		if (is_regular_file("output/test_lhco_events.lhco.gz"))
		{
			// check file size of the newly created lhco file 
			// has to be within 90 to 110 percent of old file
			unsigned int in_size = file_size("input/test_lhco_events.lhco.gz");
			unsigned int out_size = file_size("output/test_lhco_events.lhco.gz");
			if (out_size < 11 * in_size / 10 && out_size > 9 * in_size / 10)
				test_lhco_passed = true;
				
			// delete the created file as it is not needed for any checks
			remove("output/test_lhco_events.lhco.gz");
		}
	}
			
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "LHCO test: completed in " << duration << " seconds." << endl;
	cout << "LHCO reading and writing has " << (test_lhco_passed ? "succeeded!" : "failed!") << endl;
	cout << "=====================================================================" << endl;
	
	// return whether tests passed
	if (test_lhco_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}
