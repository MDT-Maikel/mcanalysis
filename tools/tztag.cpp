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
bool load_settings_general(string &output_folder);
bool load_settings_signal(string &input_sig_lhe, string &input_sig_lhco, double &sig_xsec);
bool load_settings_analysis(int &partner_mass, string &process);
PseudoJet identify_candidate_top(const vector<PseudoJet> & fatjets, vector<const particle*> & leptons);
vector<const particle*> identify_candidate_leptons(const vector<const particle*> & leptons);

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

// lepton pair mass plot
class plot_leptonmass : public plot_default
{
public:
	plot_leptonmass() {}

	double operator() (const event *ev)
	{
		return mass({ev->get(ptype_lepton, 1), ev->get(ptype_lepton, 2)});
	}
};

// top partner mass plot
class plot_thmass : public plot_default
{
public:
	plot_thmass() {}

	double operator() (const event *ev)
	{
		return mass({ev->get(ptype_lepton, 1), ev->get(ptype_lepton, 2), ev->get(ptype_jet, 1)});
	}
};

// main program 
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;

	// read the general settings from tztag.cmnd file
	string output_folder;
	if (!load_settings_general(output_folder))
		return EXIT_FAILURE;
	// read the background settings from tztag.cmnd file
	// read the signal settings from tztag.cmnd file
	string input_sig_lhe, input_sig_lhco;
	double sig_xsec = 0;
	if (!load_settings_signal(input_sig_lhe, input_sig_lhco, sig_xsec))
		return EXIT_FAILURE;
	// read the analysis settings from tztag.cmnd file
	int partner_mass = 1000;
	string signal_proc = "thth_tztz_2l";
	if (!load_settings_analysis(partner_mass, signal_proc))
		return EXIT_FAILURE;

	// basic cuts definition
	cuts basic_cuts;
	cut_2osl *osl = new cut_2osl(10, 2.5);
	basic_cuts.add_cut(osl, "2 opposite sign leptons");
	// cut_no_met *met = new cut_no_met(30);
	// basic_cuts.add_cut(met, "missing energy veto (>30 GeV)");
	// cut_ht *ht = new cut_ht(300,ptype_jet);
	// basic_cuts.add_cut(ht, "ht>300 GeV");

	// jet_analysis settings
	jet_analysis thth_tztz;
	thth_tztz.set_nEvents(5000);
	thth_tztz.undo_BDRSTagging();
	thth_tztz.set_Rsize_fat(1.5);

	// load files, shower and cluster the events
	thth_tztz.import_lhe(input_sig_lhe);
	thth_tztz.import_lhco(input_sig_lhco);
	PseudoJet (jet_analysis::*TopTagger)(const PseudoJet &) = &jet_analysis::HEPTopTagging;
	thth_tztz.initialise(TopTagger);

	// apply cuts and extract efficiencies
	double eff_basic_signal = thth_tztz.reduce_sample(basic_cuts); // require: 2 osl which reconstruct a Z, met>30 GeV veto, ht>300 GeV
	double eff_fatjpt_signal = thth_tztz.require_fatjet_pt(200, 2); // require at least 2 fatjets with pT>200 GeV
	double eff_ttag_signal = thth_tztz.require_top_tagged(2); // require 2 fatjets to be HEP Top-Tagged

	unsigned int remaining_events = thth_tztz.map_lhco_taggedJets.size();
	cout << "\n2 Fatjet pT>200 GeV Efficiency: " << setprecision(4) << 100 * eff_fatjpt_signal << " %" << endl;
	cout << "\nTop Tagging Efficiency: " << setprecision(4) << 100 * eff_ttag_signal << " %" << endl;
	cout << "\nRemaining events: " << remaining_events << endl;
	cout << "\nSignal Efficiency: " << setprecision(4) << 100 * eff_basic_signal * eff_fatjpt_signal * eff_ttag_signal << " %" << endl;
	
	// identify top partner constituents
	vector<event*> signal_lhco = thth_tztz.events();
	vector<vector<PseudoJet> > signal_fatjets = thth_tztz.fatjets();
	vector<event*> signal_reconstructed;
	for (unsigned int i = 0; i < signal_lhco.size(); ++i) // loop over events
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
	
	// plot lepton pair mass
	plot lmass("plot_leptonmass", output_folder);
	plot_leptonmass *lepton_mass = new plot_leptonmass();
	lmass.add_sample(signal_reconstructed, lepton_mass, "Signal");
	lmass.run();	
	
	// plot top partner mass
	plot pmass("plot_thmass" + lexical_cast<string>(partner_mass), output_folder);
	plot_thmass *th_mass = new plot_thmass();
	pmass.add_sample(signal_reconstructed, th_mass, "Signal");
	pmass.run();	
	
	// clear remaining pointers
	delete osl;
	// delete met;
	// delete ht;
	delete lepton_mass;
	delete th_mass;
	
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

bool load_settings_general(string &output_folder)
{
	// load the settings from the tztag.cmnd, also initialize loading variables
	string settings_file = "tztag.cmnd";
	
	// read general settings
	output_folder = read_settings<string>(settings_file, static_cast<string>("OUTPUT_FOLDER"));
	
	
	// create output directory if it does not exist yet
	if (!is_directory(output_folder))
		create_directory(output_folder);
		
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded general settings from tztag.cmnd:" << endl;
	cout << "Results output folder: " << output_folder << endl;
	cout << "################################################################################" << endl;
	return true;	
}

bool load_settings_signal(string &input_sig_lhe, string &input_sig_lhco, double &sig_xsec)
{
	// load the settings from the tztag.cmnd, also initialize loading variables
	string settings_file = "tztag.cmnd";
	
	// read signal settings
	input_sig_lhe = read_settings<string>(settings_file, static_cast<string>("INPUT_SIGNAL_LHE"));
	input_sig_lhco = read_settings<string>(settings_file, static_cast<string>("INPUT_SIGNAL_LHCO"));
	sig_xsec = read_settings<double>(settings_file, static_cast<string>("SIGNAL_XSEC"));
	
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded signal settings from tztag.cmnd:" << endl;
	cout << "Signal input LHE file: " << input_sig_lhe << endl;
	cout << "Signal input LHCO file: " << input_sig_lhco << endl;
	cout << "Signal cross section: " << sig_xsec << endl;	
	cout << "################################################################################" << endl;
	return true;
}

bool load_settings_analysis(int &partner_mass, string &process)
{
	// load the settings from the tztag.cmnd, also initialize loading variables
	string settings_file = "tztag.cmnd";
	
	// read analysis settings
	partner_mass = read_settings<int>(settings_file, static_cast<string>("PARTNER_MASS"));
	process = read_settings<string>(settings_file, static_cast<string>("PROCESS_NAME"));
	
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded analysis settings from tztag.cmnd:" << endl;
	cout << "Partner mass (GeV): " << partner_mass << endl;
	cout << "Process name: " << process << endl;
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
		double deltaR2   = sqrt( pow(deltaEta2,2.0) + pow(deltaPhi2,2.0) );

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
