/* thth->tztz tagging 
 *
 * 
*/

#include <cmath>
#include <ctime>
#include <iostream>
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

// function prototypes
bool load_settings_input(const string &settings_file, string &input_sig_lhe, string &input_sig_lhco, double &sig_xsec, int &nr_events);
bool load_settings_output(const string &settings_file, string &output_lhco, string &output_xsec);
bool load_settings_merging(const string &settings_file, bool &pythia_fast, bool &merging_on, string &merging_process, int &merging_njets, double &merging_scale);
PseudoJet identify_candidate_top(const vector<PseudoJet> & fatjets, vector<const particle*> & leptons);
vector<const particle*> identify_candidate_leptons(const vector<const particle*> & leptons);
void get_top_partner_constituents(const vector<event*> &signal_lhco, const vector<vector<PseudoJet> > &signal_fatjets, vector<event*> &signal_reconstructed);

// basic cut: at least two opposite sign leptons need to be present, with invariant mass near the Z boson
class cut_2osl : public cut
{
public:
	cut_2osl(double pt, double eta) : pt_min(pt), eta_max(eta) {}

	bool operator() (const event *ev) 
	{ 
		// extract all visible leptons
		vector<const particle*> leptons;
		for (unsigned int i = 0; i < ev->size(); ++i)
		{
			if ((*ev)[i]->type() & ptype_lepton && (*ev)[i]->pt() > pt_min && abs((*ev)[i]->eta()) < eta_max)
				leptons.push_back((*ev)[i]);
		}

		// check whether at least two visible leptons are present
		if (leptons.size() < 2)
			return false;

		// identify lepton candidates
		vector<const particle*> test_leptons = identify_candidate_leptons(leptons);
		const particle *leading_l;
		const particle *closest_l;
		leading_l = test_leptons[0];
		closest_l = test_leptons[1];
		
		if (closest_l == nullptr)
			return false;
		
		// check whether they are opposite sign leptons
		/*double charge = 0.;
		charge += leading_l->charge();
		charge += closest_l->charge();
		if ( charge != 0. )
			return false;*/

		// check whether their invariant mass is close to the Z boson mass
		double mz = 91.1876;
		double mass_range = 10.;
		double inv_mass = mass({leading_l, closest_l});
		if (inv_mass < mz - mass_range || inv_mass > mz + mass_range)
			return false;
		
		// cut passed
		return true;
	}
private:
	double pt_min;
	double eta_max;
};

// main program: may have one argument
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;
	
	// take single argument specifying settings file if available
	string settings_file = "tztag_basic.cmnd";
	if (argc >= 2) 
		settings_file = argv[1];

	// read the input settings from command file
	string input_lhe, input_lhco;
	double input_xsec = 0;
	int nr_events = 10000;
	if (!load_settings_input(settings_file, input_lhe, input_lhco, input_xsec, nr_events))
		return EXIT_FAILURE;
	// read the output settings from command file
	string output_lhco, output_xsec;
	if (!load_settings_output(settings_file, output_lhco, output_xsec))
		return EXIT_FAILURE;		
	// read the merging settings from command file
	bool pythia_fast = false, merging_on = false;
	string merging_process;
	int merging_njets = 2;
	double merging_scale = 0;
	if (!load_settings_merging(settings_file, pythia_fast, merging_on, merging_process, merging_njets, merging_scale))
		return EXIT_FAILURE;

	// basic cuts definition
	cuts basic_cuts;
	cut_2osl *osl = new cut_2osl(10, 2.5);
	basic_cuts.add_cut(osl, "2 opposite sign leptons");
	// cut_ht *ht = new cut_ht(300,ptype_jet);
	// basic_cuts.add_cut(ht, "ht>300 GeV");

	// jet_analysis settings
	jet_analysis thth_tztz;
	thth_tztz.set_nEvents(nr_events);
	thth_tztz.undo_BDRSTagging();
	thth_tztz.set_Rsize_fat(1.5);
	if (pythia_fast)
		thth_tztz.set_fast_showering();

	// load files, shower and cluster the events
	thth_tztz.import_lhe(input_lhe);
	thth_tztz.import_lhco(input_lhco);
	if (merging_on)
	{
		thth_tztz.set_merging_process(merging_process);
		thth_tztz.set_merging_njmax(merging_njets);
		thth_tztz.set_merging_scale(merging_scale);
		thth_tztz.AllowPythiaDecay();
	}	
	PseudoJet (jet_analysis::*TopTagger)(const PseudoJet &) = &jet_analysis::HEPTopTagging;
	thth_tztz.initialise(TopTagger);

	// apply cuts and extract efficiencies
	double eff_basic_signal = thth_tztz.reduce_sample(basic_cuts); // require: 2 osl which reconstruct a Z, ht>300 GeV
	double eff_fatjpt_signal = thth_tztz.require_fatjet_pt(200, 2); // require at least 2 fatjets with pT>200 GeV
	double eff_ttag_signal = thth_tztz.require_top_tagged(2); // require 2 fatjets to be HEP Top-Tagged

	unsigned int remaining_events = thth_tztz.map_lhco_taggedJets.size();
	cout << "\n2 Fatjet pT>200 GeV Efficiency: " << setprecision(4) << 100 * eff_fatjpt_signal << " %" << endl;
	cout << "\nTop Tagging Efficiency: " << setprecision(4) << 100 * eff_ttag_signal << " %" << endl;
	cout << "\nRemaining events: " << remaining_events << endl;
	cout << "\nSignal Efficiency: " << setprecision(4) << 100 * eff_basic_signal * eff_fatjpt_signal * eff_ttag_signal << " %" << endl;
	
	// identify top partner constituents
	vector<event*> signal_reconstructed;
	get_top_partner_constituents(thth_tztz.events(), thth_tztz.fatjets(), signal_reconstructed);
	
	// calculate cross section and store
	ofstream ofs_txt;
	ofs_txt.open(output_xsec.c_str());
	ofs_txt << input_xsec << "\t" << eff_basic_signal * eff_fatjpt_signal * eff_ttag_signal << "\t" << input_xsec * eff_basic_signal * eff_fatjpt_signal * eff_ttag_signal << endl;
	ofs_txt.close();
	
	// store lhco events
	write_lhco(signal_reconstructed, output_lhco);
	
	// clear remaining pointers
	delete osl;
	// delete ht;
	
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

bool load_settings_input(const string &settings_file, string &input_sig_lhe, string &input_sig_lhco, double &sig_xsec, int &nr_events)
{
	// read input settings
	input_sig_lhe = read_settings<string>(settings_file, static_cast<string>("INPUT_LHE"));
	input_sig_lhco = read_settings<string>(settings_file, static_cast<string>("INPUT_LHCO"));
	sig_xsec = read_settings<double>(settings_file, static_cast<string>("INPUT_XSEC"));
	nr_events = read_settings<int>(settings_file, static_cast<string>("NR_EVENTS"));
	
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded input settings from " << settings_file << ":" << endl;
	cout << "Input input LHE file: " << input_sig_lhe << endl;
	cout << "Input input LHCO file: " << input_sig_lhco << endl;
	cout << "Input cross section: " << sig_xsec << endl;
	cout << "Input nr of events: " << nr_events << endl;
	cout << "################################################################################" << endl;
	return true;
}

bool load_settings_output(const string &settings_file, string &output_lhco, string &output_xsec)
{
	// read output settings
	output_lhco = read_settings<string>(settings_file, static_cast<string>("OUTPUT_LHCO"));
	output_xsec = read_settings<string>(settings_file, static_cast<string>("OUTPUT_XSEC"));
		
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded output settings from " << settings_file << ":" << endl;
	cout << "Results output lhco: " << output_lhco << endl;
	cout << "Results output xsec: " << output_xsec << endl;
	cout << "################################################################################" << endl;
	return true;	
}

bool load_settings_merging(const string &settings_file, bool &pythia_fast, bool &merging_on, string &merging_process, int &merging_njets, double &merging_scale)
{
	// read merging settings
	pythia_fast = read_settings<bool>(settings_file, static_cast<string>("PYTHIA_FAST"));
	merging_on = read_settings<bool>(settings_file, static_cast<string>("MERGING_ON"));
	merging_process = read_settings<string>(settings_file, static_cast<string>("MERGING_PROCESS"));
	merging_njets = read_settings<int>(settings_file, static_cast<string>("MERGING_NJETS"));
	merging_scale = read_settings<double>(settings_file, static_cast<string>("MERGING_SCALE"));
	
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded merging settings from " << settings_file << ":" << endl;
	cout << "Pythia fast: " << pythia_fast << endl;
	cout << "Merging on: " << merging_on << endl;
	cout << "Merging process: " << merging_process << endl;
	cout << "Merging nr. jets: " << merging_njets << endl;
	cout << "Merging scale: " << merging_scale << endl;
	cout << "################################################################################" << endl;
	return true;
}

// identify lepton candidates to reconstruct the Z boson
vector<const particle*> identify_candidate_leptons(const vector<const particle*> & leptons)
{
	// identify leading-pT lepton
	const particle *leading_l = leptons[0];
	unsigned int leading_index = 0;
	for (unsigned int i = 0; i < leptons.size(); ++i)
	{
		if (leptons[i]->pt() > leading_l->pt())
		{
			leading_l = leptons[i];
			leading_index = i;
		}	
	}

	// identify lepton closest in delta_r
	const particle *closest_l = nullptr;
	double delta_r_min = 10000;
	for (unsigned int i = 0; i < leptons.size(); ++i)
	{
		if (i == leading_index)
			continue;

		if (leptons[i]->type() != leading_l->type())
			continue;
			
		if (leptons[i]->charge() + leading_l->charge() != 0)
			continue;		
		
		if (delta_r(leptons[leading_index], leptons[i]) < delta_r_min)
		{
			delta_r_min = delta_r(leptons[leading_index], leptons[i]);
			closest_l = leptons[i];
		}
	}

	// return lepton candidates
	vector<const particle*> test_leptons;
	test_leptons.push_back(leading_l);
	test_leptons.push_back(closest_l);
	return test_leptons;
}

// identify top candidate in the same emisphere of lepton candidates
PseudoJet identify_candidate_top(const vector<PseudoJet> & fatjets, vector<const particle*> & leptons)
{
	// identify top-tagged jets
	vector< PseudoJet > topjets;
	for (unsigned int i = 0; i < fatjets.size(); ++i)
	{
		if (fatjets[i].user_info<TagInfo>().top_tag())
			topjets.push_back(fatjets[i]);
	}

	// identify top candidate closest in delta_r wrt lepton candidates
	PseudoJet top_candidate;
	double delta_r_min = 10000;
	for (unsigned int i = 0; i < topjets.size(); ++i)
	{
		double deltaEta1 = topjets[i].eta() - leptons[0]->eta();
		double deltaPhi1 = abs(topjets[i].phi() - leptons[0]->phi());
		deltaPhi1 = min(deltaPhi1, 8 * atan(1) - deltaPhi1);
		double deltaR1   = sqrt( pow(deltaEta1,2.0) + pow(deltaPhi1,2.0) );

		double deltaEta2 = topjets[i].eta() - leptons[1]->eta();
		double deltaPhi2 = abs(topjets[i].phi() - leptons[1]->phi());
		deltaPhi2 = min(deltaPhi2, 8 * atan(1) - deltaPhi2);
		double deltaR2   = sqrt(pow(deltaEta2, 2.0) + pow(deltaPhi2, 2.0));

		double deltaRtest = (deltaR1+deltaR2)/2;

		if (deltaRtest < delta_r_min)
		{
			delta_r_min = deltaRtest;
			top_candidate = topjets[i];
		}
	}

	// return top candidate
	return top_candidate;
}


void get_top_partner_constituents(const vector<event*> &signal_lhco, const vector<vector<PseudoJet> > &signal_fatjets, vector<event*> &signal_reconstructed)
{
	// loop over events
	for (unsigned int i = 0; i < signal_lhco.size(); ++i)
	{
		event *newev = new event();

		// extract lepton candidates
		event *ev = signal_lhco[i];
		vector< const particle* > leptons;
		for (unsigned int j = 0; j < ev->size(); ++j)
		{
			if ((*ev)[j]->type() & ptype_lepton && (*ev)[j]->pt() > 10. && abs((*ev)[j]->eta()) < 2.5)
				leptons.push_back((*ev)[j]);
		}		
		vector<const particle*> test_leptons = identify_candidate_leptons(leptons);
		for (unsigned int j = 0; j < test_leptons.size(); ++j)
		{
			int p_type = ptype_electron;
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
		
		// add met to the event for lhco format sanity
		lhco *met = new lhco(ptype_met, 0, 0, 0, 0);
		newev->push_back(met);

		// merge leptons and tagged top into a single event pointer
		signal_reconstructed.push_back(newev);
	}
	return;
}
