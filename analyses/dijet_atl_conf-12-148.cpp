/* ATLAS dijet analysis
 *
 * Performs the ATLAS dijet analysis based on ATLAS-CONF-2012-148,
 * these limits can be compared to the model-independent limits in
 * table 1 of their conference note. This calculates the efficiency
 * for the given event file.
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

// atlas dijet cut
class cut_atlas_dijet: public cut
{
public:
	cut_atlas_dijet() {}

	bool operator() (const event *ev) 
	{ 
		// get first two jets
		const particle *j1 = ev->get(ptype_jet, 1);
		const particle *j2 = ev->get(ptype_jet, 2);
		
		// reject event if jets not found
		if (!j1 || !j2)
			return false;
		
		// reject event if |y| > 2.8 for either jet
		if (abs(j1->y()) > 2.8 || abs(j2->y()) > 2.8)
			return false;
			
		// reject event if |y*| = |y1-y2|/2 >= 0.6
		if (abs(j1->y() - j2->y()) / 2 >= 0.6)
			return false;
		
		// reject event if m(jj) <= 1000
		if (mass({j1, j2}) <= 1000)
			return false;
			
		// reject event if pt(j) < 150 GeV for either jet
		if (j1->pt() < 150 || j2->pt() < 150)
			return false;
		
		// event passed
		return true;
	}
};

// atlas dijet cut
class cut_atlas_dijet_mass: public cut
{
public:
	cut_atlas_dijet_mass(double m) : res_mass(m) {}
	
	bool operator() (const event *ev) 
	{
		// get first two jets
		const particle *j1 = ev->get(ptype_jet, 1);
		const particle *j2 = ev->get(ptype_jet, 2);
		
		// reject event if jets not found
		if (!j1 || !j2)
			return false;
			
		// reject event if not 0.8 M < m(jj) < 1.2 M
		if (mass({j1, j2}) < 0.8 * res_mass || mass({j1, j2}) > 1.2 * res_mass)
			return false;	
		
		return true;
	}
private:
	double res_mass;	
};


// main program with two arguments representing the input file and resonance mass
int main(int argc, char* argv[])
{
	// make sure there are two arguments
	if (argc != 3) 
	{
		cout << "specify at least one argument: <input.mc.gz> <mass>" << endl;
		return EXIT_FAILURE;
  	}

	// convert the arguments to the input file and output directory
	string input_file = argv[1];
	double resonance_mass = lexical_cast<double>(argv[2]);
	
	// load the events dependent on whether they are .lhe.gz or .lhco.gz
	vector<event*> events; 	
	read_lhco(events, input_file);
	
	// initialize cutlist and add all atlas dijet cuts
	cuts dijet;
	cut_atlas_dijet *dijet_cut = new cut_atlas_dijet();
	dijet.add_cut(dijet_cut, "atlas dijet basic cuts");
	cut_atlas_dijet_mass *dijet_cut_mass = new cut_atlas_dijet_mass(resonance_mass);
	dijet.add_cut(dijet_cut_mass, "atlas dijet mass range");
	dijet.apply(events);
	//atl_dijet.write(cout);
	double mod_acceptance = dijet.efficiency();
	//cout << "modified acceptance: " << setprecision(6) << mod_acceptance << endl;
	
	// calculate average mass
	unsigned int total = 0;
	double average_mass = 0;
	while (total < events.size())
	{
		const event *ev = events[total];
		average_mass += mass({ev->get(ptype_jet, 1), ev->get(ptype_jet, 2)});
		total++;
	}
	average_mass /= total;
	//cout << "average mass: " << average_mass << " GeV" << endl;
	
	// delete all the cut pointers
	delete dijet_cut;
	delete dijet_cut_mass;
		
	// clear remaining event pointers
	delete_events(events);
	
	// just put out the mod acc \t avg mass
	cout << mod_acceptance << "\t" << average_mass << endl;
	
	// finished the analysis
	return EXIT_SUCCESS;	
}
