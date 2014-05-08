/* Calculate Value
 *
 * Used by the shell scripts to calculate values of dependent parameters for the parameter card.
 * Will take at least one argument, namely a string which identifies the parameter. But also 
 * more independent parameters needed for the calculation. 
 * 
 * Compile this script and put the binary in the /scripts directory.
 * 
*/

#include <cmath>
#include <iomanip>
#include <iostream> 
#include <string>
#include <vector>

#include "boost/lexical_cast.hpp"


// functions to calculate values
double calc_MQ(std::vector<double> &pars);
double calc_gstar(std::vector<double> &pars);
double calc_width(std::vector<double> &pars);

// main program
int main(int argc, const char* argv[])
{
	// make sure there are three arguments
	if (argc <= 3) 
	{
		std::cout << "specify three arguments: <id -> MQ,gstar or WTP> <MQ> <gstar>" << std::endl;
		return EXIT_FAILURE;
  	}
	
	// get the identifier
	std::string id = argv[1];
	
	// store the parameters in a vector
	std::vector<double> pars;
	for (unsigned int i = 2; i < argc; ++i)
		pars.push_back(boost::lexical_cast<double>(argv[i]));
	
	// determine the value of the given identifier
	double value = 0; // standard zero
	
	if (id == "MQ")
		value = calc_MQ(pars); 
	else if (id == "gstar")
		value = calc_gstar(pars);
	else if (id == "WTP")
		value = calc_width(pars);

	// output the value in scientific form
	std::cout << std::setprecision(6) << std::scientific;
	std::cout << value << std::endl;

	return EXIT_SUCCESS;
}

// top partner mass
double calc_MQ(std::vector<double> &pars)
{
	return pars[0];	
}

// coupling gstar
double calc_gstar(std::vector<double> &pars)
{
	return pars[1];
}

// top partner width
double calc_width(std::vector<double> &pars)
{
	double MQ = pars[0];
	double gstar = pars[1];
	double width = pow(MQ,-3) * pow(gstar,2) * (
				1.69638*pow(10,-7)*pow(MQ,2) * (13959.0 + pow(MQ,2)) * pow( (pow(MQ,2)-pow(297,2)) * (pow(MQ,2)-pow(47,2)), 0.5) + 
				( 161.277 - 0.00834385*pow(MQ,2) + 1.64078*pow(10,-7)*pow(MQ,4) ) * pow( 4.52363*pow(10,8) - 75798.4*pow(MQ,2) + pow(MQ,4), 0.5) + 
				( -26.601 + 0.00207649*pow(MQ,2) + 3.28157*pow(10,-7)*pow(MQ,4) ) * pow( 4.03204*pow(10,7) - 12788.0*pow(MQ,2) + pow(MQ,4), 0.5)
				);
	return width;
}
