/* ttzz tagging 
 *
 * 
*/

#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>

#include <boost/filesystem.hpp>

#include "fastjet/PseudoJet.hh"

#include "particle/lhco.h"
#include "particle/particle.h"
#include "event/event.h"
#include "utility/utility.h"
#include "cuts/cuts.h"
#include "plot/plot.h"
#include "plot/plot2d.h"
#include "jet_analysis/jet_analysis.h"

using namespace std;
using namespace boost::filesystem;
using namespace Pythia8;
using namespace fastjet;
using namespace analysis;


// four lepton cut 
class cut_4lepton : public cut
{
public:
	cut_4lepton(double pt, double eta) : pt_cut(pt), eta_max(eta) {}

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
	}
private:
	double pt_cut;
	double eta_max;
};

// four lepton mass cut 
class cut_4lepton_mass : public cut
{
public:
	cut_4lepton_mass(double mass, double eta) : mass_range(mass), eta_max(eta) {}

	// event passes if two pairs of opposite sign leptons are found
	bool operator() (const event *ev) 
	{
		double mz = 91.1876;
		std::vector<const particle*> lepplus;
		std::vector<const particle*> lepminus;
		for (unsigned int i = 1; i <= 4; i++)
		{
			const particle *p = ev->get(ptype_lepton, i, eta_max);
			if (p->charge() == 1.0)
				lepplus.push_back(p);
			else
				lepminus.push_back(p);
		}
		// test first combination
		double mass1 = mass({lepplus[0], lepminus[0]});
		double mass2 = mass({lepplus[1], lepminus[1]});
		if (mass1 > mz - mass_range && mass1 < mz + mass_range &&
			mass2 > mz - mass_range && mass2 < mz + mass_range)
			return true;
		// test second combination	
		mass1 = mass({lepplus[0], lepminus[1]});
		mass2 = mass({lepplus[1], lepminus[0]});
		if (mass1 > mz - mass_range && mass1 < mz + mass_range &&
			mass2 > mz - mass_range && mass2 < mz + mass_range)
			return true;
		// not passed
		return false;
	}
private:
	double mass_range;
	double eta_max;
};

// four lepton mass plot
class plot_leptonmass : public plot_default
{
public:
	plot_leptonmass() {}

	double operator() (const event *ev)
	{
		// we know that there are four leading leptons with charge sum zero
		// combine the oppositely charged ones into a mass and sum the two mass pairs
		std::vector<const particle*> lepplus;
		std::vector<const particle*> lepminus;
		for (unsigned int i = 1; i <= 4; i++)
		{
			const particle *p = ev->get(ptype_lepton, i, 2.5);
			if (p->charge() == 1.0)
				lepplus.push_back(p);
			else
				lepminus.push_back(p);
		}
		
		return mass({lepplus[0], lepminus[0]}) + mass({lepplus[1], lepminus[1]});
	}
};

// four lepton mass plot
class plot_partnermass : public plot_default
{
public:
	plot_partnermass() {}

	double operator() (const event *ev)
	{
		double mz = 91.1876;
		double mass_range = 10;
		// we know that there are four leading leptons with charge sum zero
		// combine the oppositely charged ones into a mass and sum the two mass pairs
		std::vector<const particle*> lepplus;
		std::vector<const particle*> lepminus;
		for (unsigned int i = 1; i <= 4; i++)
		{
			const particle *p = ev->get(ptype_lepton, i, 2.5);
			if (p->charge() == 1.0)
				lepplus.push_back(p);
			else
				lepminus.push_back(p);
		}
		vector<const particle*> pair1;
		vector<const particle*> pair2;		
		
		double mass1 = mass({lepplus[0], lepminus[0]});
		double mass2 = mass({lepplus[1], lepminus[1]});
		if (mass1 > mz - mass_range && mass1 < mz + mass_range &&
			mass2 > mz - mass_range && mass2 < mz + mass_range)
		{
			pair1 = {lepplus[0], lepminus[0]};
			pair2 = {lepplus[1], lepminus[1]};				
		}
		else
		{
			pair1 = {lepplus[0], lepminus[1]};
			pair2 = {lepplus[1], lepminus[0]};	
		}
				
		pair1.push_back(ev->get(ptype_jet, 1));
		pair2.push_back(ev->get(ptype_jet, 2));
		
		return (mass(pair1) + mass(pair2)) / 2;
	}
};

// main program
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;

	// initialise jet_analysis class and do settings
	jet_analysis ttag;
	int nEvents = 50000;
	ttag.set_nEvents(nEvents);
	ttag.undo_BDRSTagging();
	ttag.set_Rsize_fat(1.5);
	ttag.set_fast_showering();
		
	// load, shower and cluster the events
	ttag.add_lhe("../../files/tools/input_ttag/input_ttag.lhe");
	fastjet::PseudoJet (jet_analysis::*TopTagger)(const fastjet::PseudoJet &) = &jet_analysis::HEPTopTagging;
	ttag.import_lhco("../../files/tools/input_ttag/input_ttag.lhco.gz");
	ttag.initialise(TopTagger);
	
	// initiate general cut class and specific cuts
	cuts ttag_cuts;
	cut_4lepton *fourlepton = new cut_4lepton(10, 2.5);
	ttag_cuts.add_cut(fourlepton, "4 leptons");
	
	// reduce the ttag sample requiring four leptons 
	// and 2 tops to be HEP top tagged
	ttag.reduce_sample(ttag_cuts);
	double eff = 0;
	int ntops = 2;
	eff = ttag.require_top_tagged(ntops);
	cout << setprecision(2);
	cout << endl << "Efficiency of tagging requirement: " << 100 * eff << " %" << endl;
	cout << "size: " << ttag.map_lhco_taggedJets.size() << endl;
	
	// get event sample from ttag
	//std::vector<event*> ttag_tagged = ttag.events();
	
	// plot lepton masses
	//plot lmass("test_plot_leptonmass", "../../files/tools/output_ttag/");
	//lmass.add_sample(ttag_tagged, new plot_leptonmass, "Tagged");
	//lmass.run();
	
	// add four lepton mass cut and apply to sample
	cut_4lepton_mass *fourleptonmass = new cut_4lepton_mass(10, 2.5);
	ttag_cuts.add_cut(fourleptonmass, "4 lepton mass");
	ttag.reduce_sample(ttag_cuts);
	
	// combine fat tops with the leptons
	vector<event*> ttag_lep = ttag.events();
	vector< vector<PseudoJet> > ttag_fatjets = ttag.fatjets();
	vector<event*> ttag_events;
	for (unsigned int i = 0; i < ttag_lep.size(); ++i)
	{
		event *newev = new event();
		// push back leptons TODO, break after four leptons.
		event *ev = ttag_lep[i];
		for (unsigned int j = 0; j < ev->size(); ++j)
		{
			if ((*ev)[j]->type() & ptype_lepton)
			{
				int p_type = ptype_lepton;
				double	p_eta 	= (*ev)[j]->eta(), 
						p_phi 	= (*ev)[j]->phi(), 
						p_pt 	= (*ev)[j]->pt(), 
						p_m 	= (*ev)[j]->mass(),
						p_charge = (*ev)[j]->charge();
				lhco *p = new lhco(p_type, p_eta, p_phi, p_pt, p_m, p_charge);
				newev->push_back(p);		
			}
		}		
		// push back tops		
		vector<PseudoJet> fatjets = ttag_fatjets[i];
		for (unsigned int j = 0; j < 2; ++j)
		{
			int p_type = ptype_jet;
			double	p_eta 	= fatjets[j].eta(), 
					p_phi 	= fatjets[j].phi(), 
					p_pt 	= fatjets[j].pt(), 
					p_m 	= fatjets[j].m();
			lhco *p = new lhco(p_type, p_eta, p_phi, p_pt, p_m);
			newev->push_back(p);
		}
		ttag_events.push_back(newev);
	}
	
	// plot lepton masses
	plot pmass("test_plot_partnermass", "../../files/tools/output_ttag/");
	pmass.add_sample(ttag_events, new plot_partnermass, "Tagged");
	pmass.run();	
	
	// clear remaining pointers
	delete fourlepton;
	delete fourleptonmass;
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Tag program completed in " << duration << " seconds." << endl;
	cout << "=====================================================================" << endl;
	
	// finished the program
	return EXIT_SUCCESS;
}

