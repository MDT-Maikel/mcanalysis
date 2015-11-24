/* ATLAS jet mass analysis
 *
 * Performs the CMS dijet analysis based on ATL_SUSY-2013-07,
 * these limits can be compared to the model-independent limits
 * from table VII in their conference note. This code calculates
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


// atl 4 jet cut
class cut_atl_4jets: public cut
{
public:
	cut_atl_4jets() {}

	bool operator() (const event *ev) 
	{ 
		// get the fourth jet
		const particle *j4 = ev->get(ptype_jet, 4, 2.5);
		
		// reject event if jets not found
		if (!j4)
			return false;
		
		// reject event if pt <= 100 GeV
		if (j4->pt() <= 100)
			return false;

		// event passed
		return true;
	}
};

// atl delta eta cut
class cut_atl_delta_eta: public cut
{
public:
	cut_atl_delta_eta() {}

	bool operator() (const event *ev) 
	{ 
		const particle *j1 = ev->get(ptype_jet, 1, 2.5);
		const particle *j2 = ev->get(ptype_jet, 2, 2.5);

		return delta_eta(j1, j2) < 0.7;
	}
};

// atl pt jet3 cut
class cut_atl_pt_jet3: public cut
{
public:
	cut_atl_pt_jet3(double pt) : pt_cut(pt) {}

	bool operator() (const event *ev) 
	{ 
		const particle *j3 = ev->get(ptype_jet, 3, 2.5);

		return j3->pt() > pt_cut;
	}
};

// atl jet mass cut
class cut_atl_jet_mass: public cut
{
public:
	cut_atl_jet_mass(double m) : jet_mass(m) {}

	bool operator() (const event *ev) 
	{ 
		const particle *j1 = ev->get(ptype_jet, 1, 2.5);
		const particle *j2 = ev->get(ptype_jet, 2, 2.5);
		const particle *j3 = ev->get(ptype_jet, 3, 2.5);
		const particle *j4 = ev->get(ptype_jet, 4, 2.5);

		double jet_mass = j1->mass() + j2->mass() + j3->mass() + j4->mass();
	
		return jet_mass < jet_mass;
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
	cuts jet_mass;
	cut_atl_4jets *fourjet_cut = new cut_atl_4jets();
	jet_mass.add_cut(fourjet_cut, "atl 4jet cut");
	cut_atl_delta_eta *eta_cut = new cut_atl_delta_eta();
	jet_mass.add_cut(eta_cut, "atl 4jet cut");
	cut_atl_pt_jet3 *jetpt_cut = new cut_atl_pt_jet3(250.0);
	jet_mass.add_cut(jetpt_cut, "atl 4jet cut");
	cut_atl_jet_mass *jetmass_cut = new cut_atl_jet_mass(625.0);
	jet_mass.add_cut(jetmass_cut, "atl jet mass cut");
	jet_mass.apply(events);
	//dijet.write(cout);
	double acceptance = jet_mass.efficiency();
	//cout << "acceptance: " << setprecision(6) << acceptance << endl;
	
	// delete all the cut pointers
	delete fourjet_cut;
	delete eta_cut;
	delete jetpt_cut;
	delete jetmass_cut;
		
	// clear remaining event pointers
	delete_events(events);
	
	// just output the acceptance
	cout << acceptance << endl;
	
	// finished the analysis
	return EXIT_SUCCESS;	
}

