/* ATLAS F_chi analysis
 *
 * Performs the ATLAS dijet analysis based on arXiv:1210.1718v3.
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
#include "plot/plot.h"
#include "utility/utility.h"



using namespace std;
using namespace analysis;
using namespace boost;
using namespace boost::filesystem;

// atlas dijet cut
class cut_atlas_fchi: public cut
{
public:
	cut_atlas_fchi() {}

	bool operator() (const event *ev) 
	{ 
		// get the two leading jets
		const particle *j1 = ev->get(ptype_jet, 1);
		const particle *j2 = ev->get(ptype_jet, 2);
		
		// reject event if jets not found
		if (!j1 || !j2)
			return false;
		
		// reject event if |y| >= 4.4 for either jet
		if (abs(j1->y()) >= 4.4 || abs(j2->y()) >= 4.4)
			return false;
			
		// reject event if |y*| = |y1-y2|/2 >= 1.7
		if (abs(j1->y() - j2->y()) / 2 >= 1.7)
			return false;
			
		// reject event if |yB| = |y1+y2|/2 >= 1.1
		if (abs(j1->y() + j2->y()) / 2 >= 1.1)
			return false;
		
		// reject event if m(jj) <= 800
		if (mass({j1, j2}) <= 800)
			return false;
		
		// event passed
		return true;
	}
};

// plot fchi
class plot_fchi : public plot_default
{
public:
	plot_fchi() {}

	double operator() (const event *ev)
	{
		const particle *j1 = ev->get(ptype_jet, 1);
		const particle *j2 = ev->get(ptype_jet, 2);
		
		return exp(abs(j1->y() - j2->y()));
	}
};


// main program with two arguments representing the input file and resonance mass
int main(int argc, char* argv[])
{
	// make sure there are two arguments
	if (argc != 2) 
	{
		cout << "specify at least one argument: <input.mc.gz>" << endl;
		return EXIT_FAILURE;
  	}

	// convert the arguments to the input file
	string input_file = argv[1];
	
	// load the events dependent on whether they are .lhe.gz or .lhco.gz
	vector<event*> events; 	
	read_lhco(events, input_file);
	
	// initialize cutlist and add all atlas fchi cuts
	cuts fchi;
	cut_atlas_fchi *fchi_cut = new cut_atlas_fchi();
	fchi.add_cut(fchi_cut, "atlas fchi basic cuts");
	fchi.apply(events);
	double acceptance = fchi.efficiency();
	
	// plot fchi
	plot p_fchi("plot_fchi", "");
	plot_fchi *f_fchi = new plot_fchi();
	p_fchi.add_sample(events, f_fchi, "signal", 1);
	p_fchi.run();	
	
	// delete all the cut pointers
	delete fchi_cut;
	delete f_fchi;
		
	// clear remaining event pointers
	delete_events(events);
	
	// just put out the mod acc \t avg mass
	cout << acceptance << endl;
	
	// finished the analysis
	return EXIT_SUCCESS;	
}
