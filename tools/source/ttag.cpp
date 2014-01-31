/* Top tagging 
 *
 * 
*/

#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>

#include <boost/filesystem.hpp>

#include "../../source/event/event.h"
#include "../../source/utility/utility.h"
#include "../../source/cuts/cuts.h"
#include "../../source/plot/plot.h"
#include "../../source/plot/plot2d.h"
#include "../../source/jet_analysis/jet_analysis.h"

using namespace std;
using namespace boost::filesystem;
using namespace Pythia8;
using namespace analysis;


// four lepton cut 
class cut_4lepton : public cut
{

public:
	cut_4lepton(double pt, double eta) : pt_cut(pt), eta_max(eta) {};

	// event passes if two pairs of opposite sign leptons are found
	bool operator() (const event *ev) 
	{ 
		// get the four leptons with highest pt
		const particle *p = ev->get(ptype_lepton, 4, eta_max);
		if (!p || p->pt() < pt_cut)
			return false;
		
		// check whether they are two pairs of opposite sign
		double charge = 0;
		for (unsigned int i = 1; i <= 4; i++)
		{
			particle *p = ev->get(ptype_lepton, i, eta_max);
			charge += p->charge();			
		}
		if (charge == 0) // PRECISION??
			return true;
		
		// not passed
		return false;
	};

private:
	double pt_cut;
	double eta_max;

};

// four lepton mass plot
double plot_leptonmass(const event *ev)
{
	// we know that there are four leading leptons with charge sum zero
	// combine the oppositely charged ones into a mass and sum the two mass pairs
	std::vector<particle*> lepplus;
	std::vector<particle*> lepminus;
	for (unsigned int i = 1; i <= 4; i++)
	{
		particle *p = ev->get(ptype_lepton, i, 5.0);
		if (p->charge() == 1.0)
			lepplus.push_back(p);
		else
			lepminus.push_back(p);
	}
	
	return mass({lepplus[0], lepminus[0]}) + mass({lepplus[1], lepminus[1]});
}

// main program
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;

	// initialise jet_analysis class
	jet_analysis ttag;
	int nEvents = 50000;
	ttag.set_nEvents(nEvents);
	ttag.undo_BDRSTagging();
	
	// load, shower and cluster the events
	ttag.add_lhe("input/thth_tztz/mass_1000.lhe");
	fastjet::PseudoJet (jet_analysis::*TopTagger)(const fastjet::PseudoJet &) = &jet_analysis::HEPTopTagging;
	ttag.initialise(TopTagger);
	
	// load the lhco and lhe events
	vector<event*> ttag_lhe;
	//vector<event*> ttag_lhco;
	read_lhe(ttag_lhe, "input/thth_tztz/mass_1000.lhe.gz");
	//read_lhe(ttag_lhco, "input/thth_tztz/mass_1000.lhco.gz");
	
	// initiate general cut class and specific cuts
	cuts ttag_cuts;
	cut_4lepton *fourlepton = new cut_4lepton(10, 5.0);
	ttag_cuts.add_cut(fourlepton, "4 leptons");
	
	// reduce the ttag sample requiring four leptons 
	// and 2 tops to be HEP top tagged
	ttag.reduce_sample(ttag_cuts);
	double eff = 0;
	int ntops = 2;
	//eff = btag.require_top_tagged(ntops);
	cout << setprecision(2);
	cout << endl << "Efficiency of tagging requirement: " << eff << endl;
	cout << "size: " << ttag.map_lhco_taggedJets.size() << endl;
	
	// get event sample from ttag
	std::vector<event*> ttag_tagged = ttag.events();
	
	// run the cuts both samples
	ttag_cuts.apply(ttag_lhe);
	ttag_cuts.write(cout);
	double eff_lhco = ttag_cuts.efficiency();
	ttag_cuts.clear();
	//ttag_cuts.apply(ttag_lhco);
	//ttag_cuts.write(cout);
	//double eff_lhe = ttag_cuts.efficiency();
	//ttag_cuts.clear();
	
	// plot lepton masses
	plot lmass("test_plot_leptonmass", "output/");
	lmass.add_sample(ttag_lhe, plot_leptonmass, "LHE");
	//lmass.add_sample(ttag_lhco, plot_leptonmass, "LHCO");
	lmass.add_sample(ttag_tagged, plot_leptonmass, "Tagged");
	lmass.run();
	
	// clear remaining pointers
	delete fourlepton;
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Tag program completed in " << duration << " seconds." << endl;
	cout << "=====================================================================" << endl;
}

