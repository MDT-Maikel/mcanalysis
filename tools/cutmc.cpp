/* Cut MC
 *
 * Performs a whole bunch of different cuts, specified by the user for the given input file.
 * 
*/

#include <iostream>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "event/event.h"
#include "utility/utility.h"
#include "cuts/cuts.h"


using namespace std;
using namespace analysis;
using namespace boost;
using namespace boost::filesystem;


// necessary function prototypes
void read_events(vector<event*> & events, const string & input_file);
void read_cuts(vector<cut*> & cutlist, vector<string> & namelist);
void read_cuts_atlas_2013_091(vector<cut*> & cutlist, vector<string> & namelist, double jet_pt = 80, unsigned int nr_jets = 6, unsigned int nr_bjets = 0);
void perform_cut_pt(vector<event*> & events);


// main program with two arguments representing the input file 
// and the output file
int main(int argc, char* argv[])
{
	// make sure there are two arguments
	if (argc != 2) 
	{
		cout << "specify at least one argument: <input.mc.gz> <config file> <output_dir/>" << endl;
		return EXIT_FAILURE;
  	}

	// convert the arguments to the input file and output directory
	string input_file = argv[1];
	//string output_dir = argv[2];
		
	// create output directory if it does not exist yet
	//if (!is_directory(output_dir))
	//	create_directory(output_dir);
	
	// load the events dependent on whether they are .lhe.gz or .lhco.gz
	vector<event*> events; 	
	read_events(events, input_file);
	
	// initialize all the cuts in a std vector
	cuts cutmc;
	vector<cut*> cutlist;
	vector<string> namelist;
	
	// read the cuts from the config file
	//read_cuts(cutlist, namelist);
	read_cuts_atlas_2013_091(cutlist, namelist);
	read_cuts_atlas_2013_091(cutlist, namelist, 100);
	read_cuts_atlas_2013_091(cutlist, namelist, 120);
	read_cuts_atlas_2013_091(cutlist, namelist, 140);
	read_cuts_atlas_2013_091(cutlist, namelist, 160);
	
	// add all the cuts and apply & print them
	for (unsigned int i = 0; i < cutlist.size(); ++i)
		cutmc.add_cut(cutlist[i], namelist[i]);	
	cutmc.apply(events);
	cutmc.write(cout);
	
	// delete all the cut pointers
	for (unsigned int i = 0; i < cutlist.size(); ++i)
		delete cutlist[i];
	cutlist.clear();
		
	// clear remaining event pointers
	delete_events(events);
	
	// finished the plotting
	return EXIT_SUCCESS;	
}

// read the events dependent on the file type
void read_events(vector<event*> & events, const string & input_file)
{
	// determine whether the input file is *.lhe.gz then load events
	std::regex lhe_match("(.*)(lhe.gz)");
	if (std::regex_match(input_file, lhe_match))
	{
		read_lhe(events, input_file);
		return;
	}
	
	// determine whether the input file is *.lhco.gz then load events
	std::regex lhco_match("(.*)(lhco.gz)");
	if (std::regex_match(input_file, lhco_match))
	{
		read_lhco(events, input_file);
		return;
	}
	
	// if reached here, the process failed and events remain unchanged
	return;
}

void read_cuts(vector<cut*> & cutlist, vector<string> & namelist)
{
	
	
		

	
}


/* Implementation of the ATLAS-CONF-2013-091 analysis
 * jet definition: anti-kt R=0.4, |eta| < 2.8, pT > 40 GeV
 * b-jet definition: same as jet with |eta| < 2.5 
 * 
 * variable cuts:
 * pT of the jets: from 80 to 220 GeV in steps of 20 GeV
 * number of jets: from 6 to 7 jets
 * number of btagged jets: from 0 to 2 b-tagged jets; TODO 
 * 
*/
void read_cuts_atlas_2013_091(vector<cut*> & cutlist, vector<string> & namelist, double jet_pt, unsigned int nr_jets, unsigned int nr_bjets)
{
	// create pt cut on the nth jet with |eta| < 2.8
	cut_pt *pt_jetn = new cut_pt(jet_pt, ptype_jet, nr_jets, 2.8);
	string name_jetn = "pt(j" + lexical_cast<string>(nr_jets) + ") > " + lexical_cast<string>(jet_pt) + " GeV";
	
	// push back cuts and names into the lists
	cutlist.push_back(pt_jetn);
	namelist.push_back(name_jetn);	
}

void perform_cut_pt(vector<event*> & events)
{
	
	
	
}

