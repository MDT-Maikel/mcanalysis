/* Merge MC
 *
 * merges different Monte Carlo files into one file
 * 
*/

#include <iostream>
#include <string>
#include <vector>

#include <getopt.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "utility/utility.h"
#include "event/event.h"

using namespace std;
using namespace boost;
using namespace boost::filesystem;
using namespace analysis;


// utility functions
void read_options(int &argc, char* argv[], bool &exit_program, string &file_type, bool &recursive, string &input_dir, string &output_file);
void print_help();
void print_version();

// main program with two arguments representing the input file 
// and the output folder for the plots
int main(int argc, char* argv[])
{
	// read command line options
	bool exit_program = false;
	string file_type = "lhco";
	bool recursive = false;
	string input_dir;
	string output_file;
	read_options(argc, argv, exit_program, file_type, recursive, input_dir, output_file);
	
	// exit program if requested
	if (exit_program)
		return EXIT_SUCCESS;
		
	// make sure the output file's directory exists
	string output_dir = output_file;
	while (output_dir.back() != '/' && output_dir.back() != '\\' && output_dir.size() > 0)
		output_dir.erase(output_dir.end() - 1);
	if (output_dir.size() > 0 && !is_directory(output_dir))
		create_directory(output_dir);
	
	// get the files to load from the folder
	vector<path> files(get_files(input_dir, "." + file_type + ".gz", recursive));
	
	// load the files into an event vector and print while loading
	vector<event*> events;
	cout << "Merging files:" << endl;
	for (unsigned int i = 0; i < files.size(); i++)
	{
		cout << "  " << files[i].string() << endl;
		read_events(events, files[i]);
	}

	// write the events into the output file and print output dir
	cout << "Done, now writing to:" << endl;
	cout << "  " << output_file << endl;
	if (file_type == "lhco")
		write_lhco(events, output_file);
	if (file_type == "lhe")
		write_lhe(events, output_file);	

	// clear remaining event pointers
	delete_events(events);

	// finished the merging
	return EXIT_SUCCESS;
}

// reads in the command line options
void read_options(int &argc, char* argv[], bool &exit_program, string &file_type, bool &recursive, string &input_dir, string &output_file)
{
	// values will be set by getopt
	extern char *optarg; 
	extern int optind;
	
	// program options
	const struct option longopts[] =
	{
		{"help",       no_argument,       0, 'h'},
		{"version",    no_argument,       0, 'v'},
		{"filetype",   required_argument, 0, 'f'},
		{"recursive",  no_argument,       0, 'r'},
		{0,            0,                 0, 0  },
	};
	
	// read in the options using getopt_long
	int index;  
	int arg = 0;
	while (arg != -1)
	{
		arg = getopt_long(argc, argv, "hvf:r", longopts, &index);

		switch (arg)
		{
		// check for --help (-h) first and print
		case 'h':
			print_help();
			exit_program = true;
			return;
		// check for --version (-v) second and print
		case 'v':
			print_version();
			exit_program = true;
			return;
		// check for --filetype=XXX (-f XXX)
		case 'f':
			file_type = lexical_cast<string>(optarg);
			break;
		// check for --recursive (-r)
		case 'r':
			recursive = true;
			break;			
		// default
		default:
			/* EMPTY */;
		}
	}
	
	// retrieve the input & output file strings, otherwise print warnings
	if (argc - optind < 2)
	{
		std::cout << "Warning: did not specify either input directory or output file." << std::endl;
		print_help();
		exit_program = true;	
	}
	else if (argc - optind > 2)
	{
		std::cout << "Warning: specified to many arguments." << std::endl;
		print_help();
		exit_program = true;	
	}
	else
	{	
		input_dir = argv[optind];
		output_file = argv[optind + 1];
	}	
}

// prints the help output to the screen
void print_help()
{
	std::cout << "Usage: mergemc [OPTION]... [INPUT DIR] [OUTPUT FILE]" << std::endl;
	std::cout << "Merges Monte Carlo event samples in LHCO and LHE format" << std::endl;
	std::cout << "from INPUT DIR and writes the resulting file to OUTPUT FILE." << std::endl;
	std::cout << std::endl;
	std::cout << "The following options are available:" << std::endl;
	std::cout << "  -h, --help        display this help and exit" << std::endl;
	std::cout << "  -v, --version     output version information and exit" << std::endl;
	std::cout << "  -f, --filetype    MC file type: lhco (default) or lhe" << std::endl;
	std::cout << "  -r, --recursive   also processes files in sub folders" << std::endl;
}

// prints the version output to the screen
void print_version()
{
	std::cout << "mergemc version 1.0" << std::endl;
}
