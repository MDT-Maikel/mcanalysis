/* Unit Tests
 *
 * Run all tests in the unit tests folder and print whether they all succeed.
 * 
 * First run the tests for the different particle classes (lhco, lhe) and then
 * for the event class based on them. After that test the histogram classes and
 * the plotting classes which are based on those. Finally test more advanced
 * classes like the bump hunter class. 
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
		
	// log start with instructions
	cout << "=====================================================================" << endl;
	cout << "=-------------------------------------------------------------------=" << endl;
	cout << "=====================================================================" << endl;
	cout << "= Starting different unit tests, all tests should pass. When all    =" << endl;
	cout << "= tests are done, a summary is given which tests have failed. Look  =" << endl;
	cout << "= at the specific test logging to see what has gone wrong within    =" << endl;
	cout << "= the test and fix the source accordingly before pushing!           =" << endl;
	cout << "=====================================================================" << endl;
	cout << "=-------------------------------------------------------------------=" << endl;
	cout << "=====================================================================" << endl;
	cout << endl << endl << endl;
	
	// run: test_lhco
	cout << "=====================================================================" << endl;
	cout << "= TETS: LHCO                                                        =" << endl;
	cout << "=====================================================================" << endl;
	bool test_lhco_passed = true;
	int result_lhco = system("./test_lhco") / 256;
	if (result_lhco == EXIT_FAILURE)
		test_lhco_passed = false;
	cout << endl << endl << endl;
			
	// run: test_lhe
	cout << "=====================================================================" << endl;
	cout << "= TEST: LHE                                                         =" << endl;
	cout << "=====================================================================" << endl;
	bool test_lhe_passed = true;
	int result_lhe = system("./test_lhe") / 256;
	if (result_lhe == EXIT_FAILURE)
		test_lhe_passed = false;
	cout << endl << endl << endl;
	
	// run: test_event
	cout << "=====================================================================" << endl;
	cout << "= TEST: EVENT                                                       =" << endl;
	cout << "=====================================================================" << endl;
	bool test_event_passed = true;
	int result_event = system("./test_event") / 256;
	if (result_event == EXIT_FAILURE)
		test_event_passed = false;
	cout << endl << endl << endl;
	
	// run: test_cuts
	cout << "=====================================================================" << endl;
	cout << "= TEST: CUTS                                                        =" << endl;
	cout << "=====================================================================" << endl;
	bool test_cuts_passed = true;
	int result_cuts = system("./test_cuts") / 256;
	if (result_cuts == EXIT_FAILURE)
		test_cuts_passed = false;
	cout << endl << endl << endl;
	
	// run: test_histogram
	cout << "=====================================================================" << endl;
	cout << "= TEST: HISTOGRAM                                                   =" << endl;
	cout << "=====================================================================" << endl;
	bool test_histogram_passed = true;
	int result_histogram = system("./test_histogram") / 256;
	if (result_histogram == EXIT_FAILURE)
		test_histogram_passed = false;
	cout << endl << endl << endl;
	
	// run: test_plot
	cout << "=====================================================================" << endl;
	cout << "= TEST: PLOT                                                        =" << endl;
	cout << "=====================================================================" << endl;
	bool test_plot_passed = true;
	int result_plot = system("./test_plot") / 256;
	if (result_plot == EXIT_FAILURE)
		test_plot_passed = false;
	cout << endl << endl << endl;
	
	// run: test_bumphunter	
	cout << "=====================================================================" << endl;
	cout << "= TEST: BUMPHUNTER                                                  =" << endl;
	cout << "=====================================================================" << endl;
	bool test_bumphunter_passed = true;
	int result_bumphunter = system("./test_bumphunter") / 256;
	if (result_bumphunter == EXIT_FAILURE)
		test_bumphunter_passed = false;
	cout << endl << endl << endl;
	
	// run: test_jets
	cout << "=====================================================================" << endl;
	cout << "= TEST: JET ANALYSIS                                                =" << endl;
	cout << "=====================================================================" << endl;
	bool test_jets_passed = true;
	int result_jets = system("./test_jets") / 256;
	if (result_jets == EXIT_FAILURE)
		test_jets_passed = false;
	cout << endl << endl << endl;
	
	// determine success of all test
	bool all_tests_passed = test_lhco_passed && test_lhe_passed && test_event_passed && test_cuts_passed
		&& test_histogram_passed && test_plot_passed && test_bumphunter_passed && test_jets_passed;
	
	// log results of all tests
	cout << "=====================================================================" << endl;
	cout << "=-------------------------------------------------------------------=" << endl;
	cout << "=====================================================================" << endl;
	cout << "= All tests completed. Results per test:                            =" << endl;
	cout << "= LHCO test has " << (test_lhco_passed ? "passed" : "failed");
	cout << "!                                             =" << endl;
	cout << "= LHE test has " << (test_lhe_passed ? "passed" : "failed");
	cout << "!                                              =" << endl;
	cout << "= Event test has " << (test_event_passed ? "passed" : "failed");
	cout << "!                                            =" << endl;
	cout << "= Cuts test has " << (test_cuts_passed ? "passed" : "failed");
	cout << "!                                             =" << endl;
	cout << "= Histogram test has " << (test_histogram_passed ? "passed" : "failed");
	cout << "!                                        =" << endl;
	cout << "= Plot test has " << (test_plot_passed ? "passed" : "failed");
	cout << "!                                             =" << endl;
	cout << "= BumpHunter test has " << (test_bumphunter_passed ? "passed" : "failed");
	cout << "!                                       =" << endl;
	cout << "= Jet analysis test has " << (test_jets_passed ? "passed" : "failed");
	cout << "!                                     =" << endl;
	cout << "=====================================================================" << endl;
	if (all_tests_passed)
	{
		cout << "= All tests have passed!                                            =" << endl;
	}
	else
	{
		cout << "= At least one of the tests failed!                                 =" << endl;
		cout << "= See above results for details.                                    =" << endl;
	}
	cout << "=====================================================================" << endl;
	cout << "=-------------------------------------------------------------------=" << endl;
	cout << "=====================================================================" << endl;
	
	// done and return
	return EXIT_SUCCESS;
}
