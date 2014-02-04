/* Generate HepMC
 * 
 * Takes the output of a usual event generator in the .lhe format then showers it with
 * Pythia8 using standard settings and then creates an output using the HepMC format.
 * 
*/

#include <iostream>

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using namespace std;
using namespace Pythia8; 


// main program with two arguments representing the in and output files
int main(int argc, const char* argv[])
{
	// make sure there are two arguments
	if (argc != 3) 
	{
		cout << "specify two arguments: <input.lhe> <output.hepmc>" << endl;
		return EXIT_FAILURE;
  	}

	// convert the arguments to the in and output strings
	string input_file = argv[1];
	string output_file = argv[2];

	// initialize Pythia with speed-up settings
	Pythia pythia;
	pythia.readString("Print:quiet=ON");
	pythia.readString("PartonLevel:MPI=OFF");
	pythia.readString("PartonLevel:Remnants=OFF");
	pythia.readString("Check:Event=OFF");
	pythia.readString("HadronLevel:all=OFF");
	
	// load LHE input file 
	pythia.init(input_file);	

	// interface for conversion from Pythia8::Event to HepMC event
	HepMC::Pythia8ToHepMC ToHepMC;

	// specify file where HepMC events will be stored
	HepMC::IO_GenEvent output_hepmc(output_file, std::ios::out);

	// begin event loop: maximum a million events
	int event_cnt = 0;
	while (event_cnt < 1000000) 
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
		HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(pythia, hepmcevt);

		// write the HepMC event to file and delete
		output_hepmc << hepmcevt;
		delete hepmcevt;
	}
	
	// succeeded with the generation and storage
	return EXIT_SUCCESS;
}
