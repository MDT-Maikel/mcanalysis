/* Merge MC
 *
 * merges different Monte Carlo files into one file
 * 
*/

#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "utility/utility.h"
#include "event/event.h"


using namespace std;
using namespace analysis;
using namespace boost::filesystem;


// main program with two arguments representing the input file 
// and the output folder for the plots
int main(int argc, char* argv[])
{
	// make sure there are two arguments
	if (argc != 3) 
	{
		cout << "specify two arguments: <input_dir/> <output_file>" << endl;
		return EXIT_FAILURE;
  	}

	// convert the arguments to the input file and output directory
	string input_dir = argv[1];
	string output_file = argv[2];
	
	// get the files to load from the folder
	vector<path> files(get_files(input_dir, ".lhco.gz", true));
	
	// print the files to load
	for (unsigned int i = 0; i < files.size(); i++)
		cout << files[i].string() << endl;
	
	// load the files into an event vector
	vector<event*> events;
	for (unsigned int i = 0; i < files.size(); i++)
		read_lhco(events, files[i]);

	// write the events into the output file
	write_lhco(events, output_file);

	// clear remaining event pointers
	delete_events(events);

	// finished the merging
	return EXIT_SUCCESS;
}
