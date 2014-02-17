/* Generate HepMC
 * 
 * Takes the output of a usual event generator in the .lhe format then showers it with
 * Pythia8 using standard settings and then creates an output using the HepMC format.
 * 
*/

#include <iostream>
#include <string>

#include <getopt.h>

#include <boost/lexical_cast.hpp>

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using namespace std;
using namespace boost;
using namespace Pythia8; 


// utility functions
void read_options(int &argc, char* argv[], bool &exit_program, bool &pythia_fast, unsigned int &max_events, string &input_file, string &output_file);
void print_help();
void print_version();

// main program 
int main(int argc, char* argv[])
{
	// read command line options
	bool exit_program = false;
	bool pythia_fast = false;
	unsigned int max_events = 1000000;
	string input_file;
	string output_file;
	read_options(argc, argv, exit_program, pythia_fast, max_events, input_file, output_file);
	
	// exit program if requested
	if (exit_program)
		return EXIT_SUCCESS;
		
	// make sure the output file's directory exists
	// TODO
	
	// initialize Pythia with speed-up settings
	Pythia pythia;
	pythia.readString("Print:quiet=ON");
	if (pythia_fast)
	{
		pythia.readString("PartonLevel:MPI=OFF");
		pythia.readString("PartonLevel:Remnants=OFF");
		pythia.readString("Check:Event=OFF");
		pythia.readString("HadronLevel:all=OFF");
	}
	
	// load LHE input file 
	pythia.init(input_file);	

	// interface for conversion from Pythia8::Event to HepMC event
	HepMC::Pythia8ToHepMC ToHepMC;

	// specify file where HepMC events will be stored
	HepMC::IO_GenEvent output_hepmc(output_file, std::ios::out);

	// begin event loop: maximum a million events
	int event_cnt = 0;
	while (event_cnt < max_events) 
	{
		// check whether the showering succeeded and the loop
		// has not arrived at the end of the file yet
		if (!pythia.next()) 
		{
			if (pythia.info.atEndOfFile())
				break;
			else
				continue;
		}
		
		// increase the event counter and log progress
		event_cnt++;
		if (event_cnt % 1000 == 0)
			cout << "processed " << event_cnt << " events" << std::endl;
			
		// construct new empty HepMC event and fill it
		// units will be as chosen for HepMC build; but can be changed
		// by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)  
		HepMC::GenEvent* hepmc_evt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(pythia, hepmc_evt);

		// write the HepMC event to file and delete
		output_hepmc << hepmc_evt;
		delete hepmc_evt;
	}
	
	// succeeded with the generation and storage
	return EXIT_SUCCESS;
}

// reads in the command line options
void read_options(int &argc, char* argv[], bool &exit_program, bool &pythia_fast, unsigned int &max_events, string &input_file, string &output_file)
{
	// values will be set by getopt
	extern char *optarg; 
	extern int optind;
	
	// program options
	const struct option longopts[] =
	{
		{"help",       no_argument,       0, 'h'},
		{"version",    no_argument,       0, 'v'},
		{"fast",       no_argument,       0, 'f'},
		{"maxevents",  required_argument, 0, 'm'},
		{0,            0,                 0, 0  },
	};
	
	// read in the options using getopt_long
	int index;  
	int arg = 0;
	while (arg != -1)
	{
		arg = getopt_long(argc, argv, "hvf", longopts, &index);

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
		// check for --fast (-f)
		case 'f':
			pythia_fast = true;
			break;
		// check for --maxevents=XXX (-m XXX)
		case 'm':
			max_events = lexical_cast<unsigned int>(optarg);
			break;			
		// default
		default:
			/* EMPTY */;
		}
	}
	
	// retrieve the input & output file strings, otherwise print warnings
	if (argc - optind < 2)
	{
		std::cout << "Warning: did not specify either input or output file." << std::endl;
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
		input_file = argv[optind];
		output_file = argv[optind + 1];
	}	
}

// prints the help output to the screen
void print_help()
{
	std::cout << "Usage: gen_hepmc [OPTION]... [INPUT FILE] [OUTPUT_FILE]" << std::endl;
	std::cout << "Runs Pythia8 on the INPUT FILE which must be in .lhe format and" << std::endl;
	std::cout << "writes the results to OUTPUT FILE in the .hepmc format." << std::endl;
	std::cout << std::endl;
	std::cout << "The following options are available:" << std::endl;
	std::cout << "  -h, --help        display this help and exit" << std::endl;
	std::cout << "  -v, --version     output version information and exit" << std::endl;
	std::cout << "  -f, --fast        turn off Pythia8 advanced options like" << std::endl;
	std::cout << "                      MPI, remnants, hadronlevel, event check" << std::endl;
	std::cout << "  -m, --maxevents   specify a maximum number of events to process" << std::endl;
}

// prints the version output to the screen
void print_version()
{
	std::cout << "gen_hepmc version 1.0" << std::endl;
}
