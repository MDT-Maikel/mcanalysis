/* Unit Tests
 *
 * Run all tests in the unit tests folder and print whether they all succeed.
 * 
*/

#include <cstdlib>
#include <iostream>

using namespace std;


// main program
int main(int argc, const char* argv[])
{	
	// check whether system is available to execute other tests
	if (!system(nullptr)) 
		exit(EXIT_FAILURE);
		
	// log start
	cout << "=====================================================================" << endl;
	cout << "=-------------------------------------------------------------------=" << endl;
	cout << "=====================================================================" << endl;	
	
	// keep track of whether all tests have passed
	bool all_tests_passed = true;	
	
	// run: test_event
	bool test_event_passed = true;
	int result_event = system("./test_event") / 256;
	if (result_event == EXIT_FAILURE)
	{
		all_tests_passed = false;
		test_event_passed = false;
	}
	
	// run: test_histogram
	bool test_histogram_passed = true;
	int result_histogram = system("./test_histogram") / 256;
	if (result_histogram == EXIT_FAILURE)
	{
		all_tests_passed = false;
		test_histogram_passed = false;
	}
	
	// run: test_bumphunter	
	bool test_bumphunter_passed = true;
	int result_bumphunter = system("./test_bumphunter") / 256;
	if (result_bumphunter == EXIT_FAILURE)
	{
		all_tests_passed = false;
		test_bumphunter_passed = false;
	}
	
	// log results of all tests
	cout << "=====================================================================" << endl;
	cout << "=-------------------------------------------------------------------=" << endl;
	cout << "=====================================================================" << endl;
	cout << "All tests completed. Results per test:" << endl;
	cout << "Event test has " << (test_event_passed ? "passed!" : "failed!") << endl;
	cout << "Histogram test has " << (test_histogram_passed ? "passed!" : "failed!") << endl;	
	cout << "BumpHunter test has " << (test_bumphunter_passed ? "passed!" : "failed!") << endl;	
	cout << "=====================================================================" << endl;
	if (all_tests_passed)
	{
		cout << "All tests have passed!" << endl;
	}
	else
	{
		cout << "Not all tests have passed!" << endl;
		cout << "See above for details." << endl;
	}
	cout << "=====================================================================" << endl;
	cout << "=-------------------------------------------------------------------=" << endl;
	cout << "=====================================================================" << endl;
	
	// done and return
	return EXIT_SUCCESS;
}
