/* Calculate Value
 *
 * Used by the shell scripts to calculate values of dependent parameters for the parameter card.
 * Will take at least one argument, namely a string which identifies the parameter. But also 
 * more independent parameters needed for the calculation. 
 * 
 * Compile this script and put the binary in the /scripts directory.
 * 
*/

#include <iostream> 
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

#include "boost/lexical_cast.hpp"


// functions to calculate values
double calc_value1();
double calc_value2();


// main program
int main(int argc, const char* argv[])
{
	// make sure there is at least one argument
	if (argc <= 1) 
	{
		cout << "specify at least one argument: <id> <par1> <par2> ... <parN>" << endl;
		return EXIT_FAILURE;
  	}
	
	// get the identifier
	string id = argv[1];
	
	// store the parameters in a vector
	vector<double> pars;
	for (unsigned int i = 2; i < argc; ++i)
		pars.push_back(boost::lexical_cast<double>(argv[i]));
	
	// determine the value of the given identifier
	double value = 0; // standard zero
	
	if (id == "<placeholder 1>") // replace me with a real identifier
		value = calc_value1(); 
	else if (id == "<placeholder 2>") // replace me with a real identifier
		value = calc_value2();

	// output the value in scientific form
	cout << setprecision(6) << scientific;
	cout << new_val << endl;

	return EXIT_SUCCESS;
}

// replace me with a real calculation
double calc_value1()
{
	return 1.0;	
}

// replace me with a real calculation
double calc_value2()
{
	return 2.0;
}
