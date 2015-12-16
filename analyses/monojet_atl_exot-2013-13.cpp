/* ATLAS monojet + met analysis
 *
 * Performs the CMS dijet analysis based on ATL-EXOT-2013-13,
 * these limits can be compared to the model-independent limits
 * from table 6 in their conference note. This code calculates
 * the efficiency for the given events file.
 * 
*/

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "cuts/cuts.h"
#include "event/event.h"
#include "particle/particle.h"
#include "utility/utility.h"



using namespace std;
using namespace analysis;
using namespace boost;
using namespace boost::filesystem;


// atl pt over met cut
class cut_atl_jetmet: public cut
{
public:
	cut_atl_jetmet() {}

	bool operator() (const event *ev) 
	{ 
		// get the leading jet and the met
		const particle *j1 = ev->get(ptype_jet, 1, 2.0);
		// reject event if jets not found
		if (!j1)
			return false;
		// event passed
		return j1->pt() > ev->met() / 2.0;
	}
};

// atl delta phi cut
class cut_atl_delta_phi: public cut
{
public:
	cut_atl_delta_phi() {}

	bool operator() (const event *ev) 
	{ 
		const particle *j1 = ev->get(ptype_jet, 1, 2.0);
		const particle *met = ev->get(ptype_met, 1);
		if (!j1 || !met)
			return false;
		return delta_phi(j1, met) > 1.0;
	}
};

// main program with one argument representing the input file 
int main(int argc, char* argv[])
{
	// make sure there is one arguments
	if (argc != 2) 
	{
		cout << "specify at least one argument: <input.mc.gz>" << endl;
		return EXIT_FAILURE;
  	}

	// convert the arguments to the input file and output directory
	string input_file = argv[1];
	
	// load the events dependent on whether they are .lhe.gz or .lhco.gz
	vector<event*> events; 	
	read_lhco(events, input_file);
	
	// initialize cutlist and add all cms dijet cuts
	cut_pt *jet_cut = new cut_pt(120.0, ptype_jet, 1, 2.0);
	cut_atl_jetmet *jetmet_cut = new cut_atl_jetmet();
	cut_atl_delta_phi *delphi_cut = new cut_atl_delta_phi();

	
	cuts atl_monojet;
	atl_monojet.add_cut(jet_cut, "jet cut");
	atl_monojet.add_cut(jetmet_cut, "pt/met cut");
	atl_monojet.add_cut(delphi_cut, "delta phi cut");
	atl_monojet.apply(events);
	double acc = atl_monojet.efficiency();

	vector<double> met_cuts = {150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 700.0};
	vector<double> acceptances = {acc};
	for (int i = 0; i < met_cuts.size(); i++)
	{
		cuts atl_monojet_met;
		cut_met *met_cut = new cut_met(met_cuts[i]);
		atl_monojet_met.add_cut(met_cut, "met cut");
		atl_monojet_met.apply(events);
		acceptances.push_back(acceptances[i] * atl_monojet_met.efficiency());
		delete met_cut;
	}
	
	// delete all the cut pointers
	delete jet_cut, jetmet_cut, delphi_cut;
		
	// clear remaining event pointers
	delete_events(events);
	
	// just output the acceptance
	for (int i = 1; i < acceptances.size(); i++)
		cout << acceptances[i] << " ";
	cout << endl;
	
	// finished the analysis
	return EXIT_SUCCESS;	
}

