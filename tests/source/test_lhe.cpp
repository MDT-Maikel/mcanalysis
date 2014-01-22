/* LHE Tests
 *
 * Test LHE class by reading and writing events.
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
#include "../../source/particle/lhe.h"
#include "../../source/utility/utility.h"

using namespace std;
using namespace boost::filesystem;
using namespace analysis;


// main program
int main(int argc, const char* argv[])
{
	// TODO: the real implementation depends on read_lhe() 
	// and write_lhe() to be implemented
	cout << "LHE test still needs to be implemented." << endl;
	return EXIT_FAILURE;	
	
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;
	
	// keep track of success
	bool test_lhe_passed = true;
		
	// check if file exists, if so read the lhco events 
	vector<event*> events;
	if (is_regular_file("input/test_lhe_events.lhe.gz"))
	{
		//read_lhe(events, "input/test_lhe_events.lhe.gz");
	}
	else
	{
		cout << "File (input/test_lhe_events.lhe.gz) for reading lhe events not available." << endl;
		test_lhe_passed = false;
	}
		
	// if reading has succeeded write the events to file
	if (test_lhe_passed)
	{
		//write_lhe(events, "output/test_lhe_events.lhe.gz");
		
		// test whether the writing has succeeded
		test_lhe_passed = false;
		if (is_regular_file("output/test_lhco_events.lhe.gz"))
		{
			// check file size of the newly created lhe file 
			// has to be within 90 to 110 percent of old file
			unsigned int in_size = file_size("input/test_lhe_events.lhe.gz");
			unsigned int out_size = file_size("output/test_lhe_events.lhe.gz");
			if (out_size < 11 * in_size / 10 && out_size > 9 * in_size / 10)
				test_lhe_passed = true;
				
			// delete the created file as it is not needed for any checks
			remove("output/test_lhe_events.lhe.gz");
		}
	}
			
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "LHE test: completed in " << duration << " seconds." << endl;
	cout << "LHE reading and writing has " << (test_lhe_passed ? "succeeded!" : "failed!") << endl;
	cout << "=====================================================================" << endl;
	
	// return whether tests passed
	if (test_lhe_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}
