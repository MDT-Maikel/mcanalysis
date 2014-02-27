/* Generate HepMC
 * 
 * Takes the output of a usual event generator in the .lhe format then showers it with
 * Pythia8 using standard settings and then creates an output using the HepMC format.
 * 
*/

#include <iostream>
#include <cstdlib>
#include <string>

#include <getopt.h>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using namespace std;
using namespace boost;
using namespace boost::filesystem;
using namespace Pythia8; 


// utility functions
void read_options(int &argc, char* argv[], bool &exit_program, bool &pythia_fast, unsigned int &max_events, 
	bool &merging, string &process, int &njmax, int &njadd, double &scale, string &input_file, string &output_file);
void print_help();
void print_version();

// main program 
int main(int argc, char* argv[])
{
	//============= Read command line options =============//
	bool exit_program = false;
	bool pythia_fast = false;

	unsigned int max_events = 1000000;

	bool merging = false;
	string merge_process;
	int merge_njmax;
	int merge_njadd;
	double merge_scale;

	string input_file;
	string output_file;

	read_options(argc, argv, exit_program, pythia_fast, max_events, merging, merge_process, merge_njmax, merge_njadd, merge_scale, input_file, output_file);
	
	// exit program if requested
	if (exit_program)
		return EXIT_SUCCESS;

	//============= Basic setup =============//
	// Make sure the output file's directory exists
	string output_dir = output_file;
	while (output_dir.back() != '/' && output_dir.back() != '\\' && output_dir.size() > 0)
		output_dir.erase(output_dir.end() - 1);
	if (!is_directory(output_dir))
		create_directory(output_dir);

	// Interface for conversion from Pythia8::Event to HepMC event
	HepMC::Pythia8ToHepMC ToHepMC;

	// Specify file where HepMC events will be stored
	HepMC::IO_GenEvent output_hepmc(output_file, std::ios::out);
	
	// Pythia settings and merging variables
	Pythia pythia;
	pythia.settings.flag("Print:quiet", true);
	if ( pythia_fast )
	{
		pythia.settings.flag("PartonLevel:MPI", false);
		pythia.settings.flag("PartonLevel:Remnants", false);
		pythia.settings.flag("Check:Event", false);
		pythia.settings.flag("HadronLevel:all", false);
	}
	int njetmerged, njetcounterLO;

	// Merging setup
	if ( merging )
	{
		// Check merging settings
		if ( (merge_process == "") || (merge_njmax == 0) || (merge_scale == 0.) )
		{
			cout << "Error while initialising Pythia (Merging settings)" << endl;
			exit (EXIT_FAILURE);
		}
		pythia.settings.flag("Merging:doKTMerging", true);
		pythia.settings.word("Merging:Process", merge_process);
		pythia.settings.mode("Merging:nJetMax", merge_njmax);
		pythia.settings.parm("Merging:TMS", merge_scale);
		njetmerged = merge_njmax;
		njetcounterLO = merge_njadd;
	}
	else
		njetcounterLO = 0; // only 0-jet sample

	//============= Cross section estimation procedure =============//
	// Save estimates in vectors
	vector< double > xsecLO;
	vector< double > nAcceptLO;
	bool fsr, isr, mpi, had;

	if ( merging )
	{
		cout << "\n\n ================ Start cross section estimation ================" << endl << endl;
		
		// Switch off all showering and MPI when extimating the cross section after the merging scale cut
		fsr = pythia.flag("PartonLevel:FSR");
		isr = pythia.flag("PartonLevel:ISR");
		mpi = pythia.flag("PartonLevel:MPI");
		had = pythia.flag("HadronLevel:all");
		pythia.settings.flag("PartonLevel:FSR", false);
		pythia.settings.flag("PartonLevel:ISR", false);
		pythia.settings.flag("PartonLevel:MPI", false);
		pythia.settings.flag("HadronLevel:all", false);

		pythia.settings.flag("Merging:doXSectionEstimate", true);

		while(njetcounterLO >= 0) 
		{
			// Set appropriate LHE file name
			string lhe_file;
			if ( njetcounterLO == 0 )
				lhe_file = input_file + ".lhe";
			else
				lhe_file = input_file + "_j" + lexical_cast<string>(njetcounterLO) + ".lhe";

			// LHE initialisation
			pythia.settings.mode("Merging:nRequested", njetmerged);
			pythia.settings.word("Beams:LHEF", lhe_file);
			if ( !pythia.init(lhe_file) )
			{
				cout << "Error while initialising Pythia (pythia.init)" << endl;
				exit (EXIT_FAILURE);
			}

			// Start event loop
			for( int iEvent = 0; iEvent < max_events; ++iEvent )
			{
				// Generate event
				if ( !pythia.next() )
				{
					if( pythia.info.atEndOfFile() ) 
						break;
					else 
						continue;
				}
			} // End of event loop

			// Store cross section
			xsecLO.push_back(pythia.info.sigmaGen());
			nAcceptLO.push_back(pythia.info.nAccepted());

			// Restart with ME of a reduced the number of jets
			if( njetcounterLO > 0 )
			{
				njetcounterLO--;
				njetmerged--;
			}
			else
				break;

		} // End while( njetcounterLO>=0 )

		// Reset values
		njetmerged = merge_njmax;
		njetcounterLO = merge_njadd;
	}


	//============= Event generation and matching =============//
		cout << "\n\n ================ Start Event analysis ================" << endl;

  		// Cross section and error variables
		double sigmaTotal  = 0.;
		double errorTotal  = 0.;
		int sizeLO = static_cast<int>(xsecLO.size()), iNow;

  		// Loop over different LHE files with additional external jets
		while(njetcounterLO >= 0)
		{
			// Set appropriate LHE file name
			string lhe_file;
			if ( njetcounterLO == 0 )
				lhe_file = input_file + ".lhe";
			else
				lhe_file = input_file + "_j" + lexical_cast<string>(njetcounterLO) + ".lhe";

			// Additional merging settings: LHE input and total jet to be merged
			if ( merging )
			{
				// Possibly switch showering and multiple interaction back on
				pythia.settings.flag("Merging:doXSectionEstimate", false);
				pythia.settings.flag("PartonLevel:FSR", fsr);
  				pythia.settings.flag("PartonLevel:ISR", isr);
  				pythia.settings.flag("PartonLevel:MPI", mpi);
  				pythia.settings.flag("HadronLevel:all", had);

				pythia.settings.mode("Merging:nRequested", njetmerged);
    			pythia.settings.word("Beams:LHEF", lhe_file);
    			iNow = sizeLO-1-njetcounterLO;
    			cout << "\n\nStart analysis of " << njetcounterLO << " jets sample" << endl;
			}

			// LHE initialisation
			if ( !pythia.init(lhe_file) )
			{
				cout << "Error while initialising Pythia (pythia.init)" << endl;
				exit (EXIT_FAILURE);
			}


			//============= Event loop =============//
			for (int iEvent = 0; iEvent < max_events; ++iEvent)
			{
				// Generate event
				if ( !pythia.next() )
				{
					if( pythia.info.atEndOfFile() ) 
						break;
					else 
						continue;
				}

				// Get event weight(s) of merging procedure
				double weight = 1.0, evtweight;
				if ( merging )
				{
					weight = pythia.info.mergingWeight();
					evtweight = pythia.info.weight();
					weight *= evtweight;
				}

				// Do not consider zero-weight events in merging procedure
				if ( weight == 0. ) 
					continue;

				// Increase the event counter and log progress
				if ( (iEvent+1)%100 == 0 && njetcounterLO == 0 )
					cout << "processed " << iEvent+1 << " events" << "\r" << flush;

				// Construct new empty HepMC event and fill it
				// units will be as chosen for HepMC build; but can be changed
				// by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)  
				HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();

				// Get correct cross section from previous estimate and set event weigth
				double normhepmc;
				if ( merging )
				{
					normhepmc = xsecLO[iNow] / nAcceptLO[iNow];
					sigmaTotal += weight*normhepmc;
					errorTotal += Pythia8::pow2(weight*normhepmc);
					hepmcevt->weights().push_back(weight*normhepmc);
				}
				
				// Fill HepMC event
				ToHepMC.fill_next_event( pythia, hepmcevt );

				// Report cross section to hepmc
				if ( merging )
				{
					HepMC::GenCrossSection xsec;
					xsec.set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
					hepmcevt->set_cross_section( xsec );
				}

				// Write the HepMC event to file and delete
				output_hepmc << hepmcevt;
				delete hepmcevt;

			} // End of event loop


		//============= Restart with ME of a reduced the number of jets =============//
		if( njetcounterLO > 0 )
		{
			njetcounterLO--;
			njetmerged--;
		}
		else
			break;

	} // End while( njetcounterLO>=0 )
	
	// succeeded with the generation and storage
	return EXIT_SUCCESS;
}

// reads in the command line options
void read_options(int &argc, char* argv[], bool &exit_program, bool &pythia_fast, unsigned int &max_events, 
	bool &merging, string &process, int &njmax, int &njadd, double &scale, string &input_file, string &output_file)
{
	// values will be set by getopt
	extern char *optarg; 
	extern int optind;
	
	// program options
	const struct option longopts[] =
	{
		{"help",        	no_argument,       0, 'h'},
		{"version",     	no_argument,       0, 'v'},
		{"fast",        	no_argument,       0, 'f'},
		{"maxevents",   	required_argument, 0, 'm'},
		{"merge_process", 	required_argument, 0, 'n'},
		{"merge_njmax", 	required_argument, 0, 'o'},
		{"merge_njadd",  	required_argument, 0, 'p'},
		{"merge_scale", 	required_argument, 0, 'q'},
		{0,             	0,                 0, 0  },
	};
	
	// read in the options using getopt_long
	int index;  
	int arg = 0;
	while (arg != -1)
	{
		arg = getopt_long(argc, argv, "hvfm:n:o:p:q:", longopts, &index);

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

		// check for --merge_process=XXX (-n XXX)
		case 'n':
			process = optarg;
			merging = true;
			break;

		// check for --merge_njmax=XXX (-o XXX)
		case 'o':
			njmax = lexical_cast<int>(optarg);
			njadd = njmax;
			break;

		// check for --merge_njadd=XXX (-p XXX)
		case 'p':
			njadd = lexical_cast<int>(optarg);
			break;	

		// check for --merge_scale=XXX (-q XXX)
		case 'q':
			scale = lexical_cast<double>(optarg);
			break;	

		// default
		default:
			/* EMPTY */;
		}
	}
	
	// retrieve the input & output file strings, otherwise print warnings
	if (argc - optind < 2)
	{
		cout << "Warning: did not specify either input or output file." << endl;
		print_help();
		exit_program = true;	
	}
	else if (argc - optind > 2)
	{
		cout << "Warning: specified to many arguments." << endl;
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
	cout << "Usage: gen_hepmc [OPTION]... [INPUT FILE] [OUTPUT_FILE]" << endl;
	cout << "Runs Pythia8 on the INPUT FILE which must be in .lhe format and" << endl;
	cout << "writes the results to OUTPUT FILE in the .hepmc format." << endl;
	cout << endl;
	cout << "The following options are available:" << endl;
	cout << "  -h, --help         	display this help and exit" << endl;
	cout << "  -v, --version      	output version information and exit" << endl;
	cout << "  -f, --fast         	turn off Pythia8 advanced options like" << endl;
	cout << "                       MPI, remnants, hadronlevel, event check" << endl;
	cout << "  -m, --maxevents    	specify a maximum number of events to process" << endl;
	cout << "  -n, --merge_process  specify the merging hard process" << endl;
	cout << "  -o, --merge_njmax  	specify the max number of jets to be merged" << endl;
	cout << "  -p, --merge_njadd  	specify the additional number of jets besides the hard process" << endl;
	cout << "  -q, --merge_scale   	specify the merging scale" << endl;
}

// prints the version output to the screen
void print_version()
{
	cout << "gen_hepmc version 1.0" << endl;
}
