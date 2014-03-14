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
using namespace boost;
using namespace boost::filesystem;
using namespace Pythia8;
using namespace fastjet;
using namespace analysis;


// ====================================================== // 
//	          Auxiliary functions and classes 		   	  //
// ====================================================== //

// identify lepton candidates to reconstruct the Z boson
vector< const particle* > identify_candidate_leptons(const vector< const particle* > & leptons)
{
	// identify leading-pT lepton
	const particle *leading_l = leptons[0];
	unsigned int leading_index = 0;
	for (unsigned int i = 0; i < leptons.size(); ++i)
	{
		if ( leptons[i]->pt() > leading_l->pt() )
		{
			leading_l = leptons[i];
			leading_index = i;
		}	
	}

	// identify lepton closest in delta_r
	const particle *closest_l = new particle;
	double delta_r_min = 10000;
	for (unsigned int i = 0; i < leptons.size(); ++i)
	{
		if ( i!=leading_index && delta_r(leptons[leading_index],leptons[i]) < delta_r_min )
		{
			delta_r_min = delta_r(leptons[leading_index],leptons[i]);
			closest_l = leptons[i];
		}
	}

	vector< const particle* > test_leptons;
	test_leptons.push_back(leading_l);
	test_leptons.push_back(closest_l);

	// return lepton candidates
	return test_leptons;
};


// identify top candidate in the same emisphere of lepton candidates
PseudoJet identify_candidate_top(const vector< PseudoJet > & fatjets, vector< const particle* > & leptons)
{
	// identify top-tagged jets
	vector< PseudoJet > topjets;
	for (unsigned int i = 0; i < fatjets.size(); ++i)
	{
		if ( fatjets[i].user_info<TagInfo>().top_tag() )
			topjets.push_back(fatjets[i]);
	}

	// identify top candidate closest in delta_r wrt lepton candidates
	PseudoJet top_candidate;
	double delta_r_min = 10000;
	for (unsigned int i = 0; i < topjets.size(); ++i)
	{
		double deltaEta1 = topjets[i].eta() - leptons[0]->eta();
		double deltaPhi1 = topjets[i].phi() - leptons[0]->phi();
		double deltaR1   = sqrt( pow(deltaEta1,2.0) + pow(deltaPhi1,2.0) );

		double deltaEta2 = topjets[i].eta() - leptons[1]->eta();
		double deltaPhi2 = topjets[i].phi() - leptons[1]->phi();
		double deltaR2   = sqrt( pow(deltaEta2,2.0) + pow(deltaPhi2,2.0) );

		double deltaRtest = (deltaR1+deltaR2)/2;

		if ( deltaRtest < delta_r_min )
		{
			delta_r_min = deltaRtest;
			top_candidate = topjets[i];
		}
	}

	// return top candidate
	return top_candidate;
};


// basic cut: at least two opposite sign leptons need to be present, with invariant mass near the Z boson
class cut_2osl : public cut
{
public:
	cut_2osl(double pt, double eta) : pt_min(pt), eta_max(eta) {}

	bool operator() (const event *ev) 
	{ 
		// extract all visible leptons
		vector< const particle* > leptons;
		for (unsigned int i = 0; i < ev->size(); ++i)
		{
			if ( (*ev)[i]->type() & ptype_lepton && (*ev)[i]->pt() > pt_min && abs((*ev)[i]->eta()) < eta_max )
				leptons.push_back((*ev)[i]);
		}

		// check whether at least two visible leptons are present
		if( leptons.size() < 2 )
			return false;

		// identify lepton candidates
		vector< const particle* > test_leptons = identify_candidate_leptons(leptons);
		const particle *leading_l = new particle;
		const particle *closest_l = new particle;
		leading_l = test_leptons[0];
		closest_l = test_leptons[1];
		
		// check whether they are opposite sign leptons
		double charge = 0.;
		charge += leading_l->charge();
		charge += closest_l->charge();
		if ( charge != 0. )
			return false;

		// check whether their invariant mass is close to the Z boson mass
		double mz = 91.1876;
		double mass_range = 10.;
		double inv_mass = mass({leading_l, closest_l});
		if ( !(inv_mass > mz - mass_range && inv_mass < mz + mass_range) )
			return false;
		
		// cut passed
		return true;
	}
private:
	double pt_min;
	double eta_max;
};

// basic cut: upper bound on missing transverse energy of the event
class cut_no_met : public cut
{		
public:
	cut_no_met(double met) : met_cut(met) {};
	
	bool operator() (const event *ev) 
	{ 
		return ev->met() < met_cut;
	};		
private:
	double met_cut;		
};


// Top partner mass plot
class plot_THmass : public plot_default
{
public:
	plot_THmass() {}

	double operator() (const event *ev)
	{
		vector<const particle*> constituents;
		constituents.push_back( ev->get(ptype_lepton, 1) );
		constituents.push_back( ev->get(ptype_lepton, 2) );
		constituents.push_back( ev->get(ptype_jet, 1) );

		return mass(constituents);
	}
};


// ================================== // 
//	          Main program 		   	  //
// ================================== //
int main(int argc, const char* argv[])
{
	//===== initialisation =====//
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;

	// folder definitions
	// string mc_folder = "/data/btag/mc/";
	// string event_folder = "/data/btag/events/";
	// string plot_folder = "/data/btag/plots/";
	string mc_folder = "/Users/marco/Documents/University/PhD/Works/MC_data/MadGraph/lht_ttag/";
	string event_folder = "/Users/marco/Documents/University/PhD/Works/LHT_tagged_top/events/";
	string plot_folder = "/Users/marco/Documents/University/PhD/Works/LHT_tagged_top/plots/";

	//===== basic cuts definition =====//
	cuts basic_cuts;
	cut_2osl *osl = new cut_2osl(10, 2.5);
	basic_cuts.add_cut(osl, "2 opposite sign leptons");
	// cut_no_met *met = new cut_no_met(30);
	// basic_cuts.add_cut(met, "missing energy veto (>30 GeV)");
	// cut_ht *ht = new cut_ht(300,ptype_jet);
	// basic_cuts.add_cut(ht, "ht>300 GeV");

	//===== signal analysis =====//
	// signal info
	int partner_mass = 1000;
	string signal_proc = "thth_tztz_2l";

	// jet_analysis settings
	jet_analysis thth_tztz;
	thth_tztz.set_nEvents(50000);
	thth_tztz.undo_BDRSTagging();
	thth_tztz.set_Rsize_fat(1.5);

	// load files, shower and cluster the events
	thth_tztz.import_lhe(mc_folder + "signal/" + signal_proc + "/Events/run_mass" + lexical_cast<string>(partner_mass));
	thth_tztz.import_lhco(event_folder + "signal/" + signal_proc + "/run_mass" + lexical_cast<string>(partner_mass) + ".lhco.gz");
	// thth_tztz.import_lhe("../../files/tools/input/thth_tztz/mass_1000");
	// thth_tztz.import_lhco("../../files/tools/input/thth_tztz/mass_1000.lhco.gz");
	PseudoJet (jet_analysis::*TopTagger)(const PseudoJet &) = &jet_analysis::HEPTopTagging;
	thth_tztz.initialise(TopTagger);

	// apply cuts and extract efficiencies
	double eff_basic_signal = thth_tztz.reduce_sample(basic_cuts); // require: 2 osl which reconstruct a Z, met>30 GeV veto, ht>300 GeV
	double eff_fatjpt_signal = thth_tztz.require_fatjet_pt(200,2); // require at least 2 fatjets with pT>200 GeV
	double eff_ttag_signal = thth_tztz.require_top_tagged(2); // require 2 fatjets to be HEP Top-Tagged

	unsigned int remaining_events = thth_tztz.map_lhco_taggedJets.size();
	cout << "\n2 Fatjet pT>200 GeV Efficiency: " << setprecision(4) << 100 * eff_fatjpt_signal << " %" << endl;
	cout << "\nTop Tagging Efficiency: " << setprecision(4) << 100 * eff_ttag_signal << " %" << endl;
	cout << "\nRemaining events: " << remaining_events << endl;
	cout << "\nSignal Efficiency: " << setprecision(4) << 100 * eff_basic_signal * eff_fatjpt_signal * eff_ttag_signal << " %" << endl;
	
	// identify Top partner constituents
	vector< event* > signal_lhco = thth_tztz.events();
	vector< vector< PseudoJet > > signal_fatjets = thth_tztz.fatjets();
	vector< event* > signal_reconstructed;
	for (unsigned int i = 0; i < signal_lhco.size(); ++i) // loop over events
	{
		event *newev = new event();

		// extract lepton candidates
		event *ev = signal_lhco[i];
		vector< const particle* > leptons;
		for (unsigned int j = 0; j < ev->size(); ++j)
		{
			if ( (*ev)[j]->type() & ptype_lepton && (*ev)[j]->pt() > 10. && abs((*ev)[j]->eta()) < 2.5 )
				leptons.push_back((*ev)[j]);
		}		
		vector< const particle* > test_leptons = identify_candidate_leptons(leptons);
		for (unsigned int j = 0; j < test_leptons.size(); ++j)
		{
			int p_type = ptype_lepton;
			double	p_eta 	 = test_leptons[j]->eta(), 
					p_phi 	 = test_leptons[j]->phi(), 
					p_pt 	 = test_leptons[j]->pt(), 
					p_m 	 = test_leptons[j]->mass(),
					p_charge = test_leptons[j]->charge();
			lhco *p = new lhco(p_type, p_eta, p_phi, p_pt, p_m, p_charge);
			newev->push_back(p);
		}

		// extract top candidate		
		vector< PseudoJet > fatjets = signal_fatjets[i];
		PseudoJet top_candidate = identify_candidate_top(fatjets, test_leptons);
		int p_type = ptype_jet;
		double	p_eta 	= top_candidate.eta(), 
				p_phi 	= top_candidate.phi(), 
				p_pt 	= top_candidate.pt(), 
				p_m 	= top_candidate.m();
		lhco *p = new lhco(p_type, p_eta, p_phi, p_pt, p_m);
		newev->push_back(p);

		// merge leptons and tagged top into a single event pointer
		signal_reconstructed.push_back(newev);
	}
	
	//===== plot top partner mass =====//
	plot mass("plot_THmass" + lexical_cast<string>(partner_mass), plot_folder);
	plot_THmass *THmass = new plot_THmass();
	mass.add_sample(signal_reconstructed, THmass, "Signal");
	mass.run();	
	
	//===== finalise program =====//
	// clear remaining pointers
	delete osl;
	// delete met;
	// delete ht;
	delete THmass;
	
	// clear remaining event pointers
	delete_events(signal_reconstructed);
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Program completed in " << duration << " seconds." << endl;
	cout << "=====================================================================" << endl;
	
	// finished the program
	return EXIT_SUCCESS;
}

