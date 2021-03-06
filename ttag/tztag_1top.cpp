/* thth->tztz 1 top tagging 
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
vector<const particle*> identify_candidate_leptons(const vector<const particle*> & leptons);
PseudoJet identify_candidate_top(const vector<PseudoJet> & fatjets, vector<const particle*> & leptons);
double cut_ptT(jet_analysis &analysis, double pt_min);
double cut_etaT(jet_analysis &analysis, double eta_max);
double cut_deltaRtb(jet_analysis &analysis, double delta_r);
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
		vector<const particle*> l_candidates = identify_candidate_leptons(leptons);
		const particle *first_l = l_candidates[0];
		const particle *second_l = l_candidates[1];
		
		if (first_l == nullptr || second_l == nullptr)
			return false;
		
		// cut passed
		return true;
	}
private:
	double pt_min;
	double eta_max;
};

// basic cut: pt of the reconstructed Z boson
class cut_ptZ : public cut
{
public:
	cut_ptZ(double pt) : pt_min(pt) {}

	bool operator() (const event *ev) 
	{ 
		// extract all visible leptons
		vector<const particle*> leptons;
		for (unsigned int i = 0; i < ev->size(); ++i)
		{
			if ((*ev)[i]->type() & ptype_lepton && (*ev)[i]->pt() > 25. && abs((*ev)[i]->eta()) < 2.5)
				leptons.push_back((*ev)[i]);
		}

		// identify lepton candidates (notice that the osl leptons reconstructing the Z need to be present)
		vector<const particle*> l_candidates = identify_candidate_leptons(leptons);
		const particle *first_l = l_candidates[0];
		const particle *second_l = l_candidates[1];
		
		// reconstruct Z-boson transverse momentum
		double px_Z = first_l->px() + second_l->px();
		double py_Z = first_l->py() + second_l->py();
		double pt_Z = sqrt(pow(px_Z, 2.0) + pow(py_Z, 2.0));

		// check if pT(Z) > pt_min
		if (pt_Z < pt_min)
			return false;
		
		// cut passed
		return true;
	}
private:
	double pt_min;
};

// basic cut: eta of the reconstructed Z boson
class cut_etaZ : public cut
{
public:
	cut_etaZ(double eta) : eta_max(eta) {}

	bool operator() (const event *ev) 
	{ 
		// extract all visible leptons
		vector<const particle*> leptons;
		for (unsigned int i = 0; i < ev->size(); ++i)
		{
			if ((*ev)[i]->type() & ptype_lepton && (*ev)[i]->pt() > 25. && abs((*ev)[i]->eta()) < 2.5)
				leptons.push_back((*ev)[i]);
		}

		// identify lepton candidates (notice that the osl leptons reconstructing the Z need to be present)
		vector<const particle*> l_candidates = identify_candidate_leptons(leptons);
		const particle *first_l = l_candidates[0];
		const particle *second_l = l_candidates[1];
		
		// reconstruct Z-boson 4-vector
		const double px_Z = first_l->px() + second_l->px();
		const double py_Z = first_l->py() + second_l->py();
		const double pz_Z = first_l->pz() + second_l->pz();
		const double pe_Z = first_l->pe() + second_l->pe();
		PseudoJet Zboson(px_Z, py_Z, pz_Z, pe_Z);

		// extract eta of Z
		double eta_Z = abs(Zboson.eta());

		// check if eta(Z) < eta_max
		if (eta_Z > eta_max)
			return false;
		
		// cut passed
		return true;
	}
private:
	double eta_max;
};

// basic cut: at least n visible jets
class cut_njet : public cut
{
public:
	cut_njet(int n, double pt, double eta) : n_j(n), pt_min(pt), eta_max(eta) {}

	bool operator() (const event *ev)
	{ 
		// extract all visible jets
		vector<const particle*> jets;
		for (unsigned int i = 0; i < ev->size(); ++i)
		{
			if ((*ev)[i]->type() & ptype_jet && (*ev)[i]->pt() > pt_min && abs((*ev)[i]->eta()) < eta_max)
				jets.push_back((*ev)[i]);
		}

		// check whether at least n visible jets are present
		if (jets.size() < n_j)
			return false;
		
		// cut passed
		return true;
	}
private:
	int n_j;
	double pt_min;
	double eta_max;
};

// basic cut: at least n visible b-jets
class cut_bjet : public cut
{
public:
	cut_bjet(int n, double pt, double eta) : n_b(n), pt_min(pt), eta_max(eta) {}

	bool operator() (const event *ev)
	{ 
		// extract all visible b-jets
		vector<const particle*> bjets;
		for (unsigned int i = 0; i < ev->size(); ++i)
		{
			if ((*ev)[i]->bjet() != 0.0 && (*ev)[i]->pt() > pt_min && abs((*ev)[i]->eta()) < eta_max)
				bjets.push_back((*ev)[i]);
		}

		// check whether at least n visible b-jets are present
		if (bjets.size() < n_b)
			return false;
		
		// cut passed
		return true;
	}
private:
	int n_b;
	double pt_min;
	double eta_max;
};

// main program: may have one argument
int main(int argc, const char* argv[])
{
	// initialise timing procedure
	clock_t clock_old = clock();
	double duration;
	
	// take single argument specifying settings file if available
	string settings_file = "tztag_1top.cmnd";
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

	// jet_analysis initialisation (settings, load files, shower and cluster the events)
	jet_analysis thth_tztz;
	thth_tztz.set_nEvents(nr_events);
	thth_tztz.undo_BDRSTagging();
	thth_tztz.set_Rsize_fat(1.5);
	PseudoJet (jet_analysis::*TopTagger)(const PseudoJet &) = &jet_analysis::HEPTopTagging;
	if (pythia_fast)
		thth_tztz.set_fast_showering();

	thth_tztz.import_lhe(input_lhe);
	thth_tztz.import_lhco(input_lhco);
	if (merging_on)
	{
		thth_tztz.set_merging_process(merging_process);
		thth_tztz.set_merging_njmax(merging_njets);
		thth_tztz.set_merging_scale(merging_scale);
		thth_tztz.AllowPythiaDecay();
	}
	thth_tztz.initialise(TopTagger);

	// basic cuts definition
	cuts basic_cuts;
	cut_2osl *osl = new cut_2osl(25, 2.5);
	basic_cuts.add_cut(osl, "2 opposite sign leptons within R=1.0 cone");
	cut_ptZ *ptZ = new cut_ptZ(225);
	basic_cuts.add_cut(ptZ, "pT(Z)>225 GeV");
	cut_etaZ *etaZ = new cut_etaZ(2.3);
	basic_cuts.add_cut(etaZ, "eta(Z)<2.3");
	cut_ht *ht = new cut_ht(400, ptype_jet, 30, 3.0);
	basic_cuts.add_cut(ht, "ht>400 GeV");
	cut_njet *njet = new cut_njet(0, 30, 3.0);
	basic_cuts.add_cut(njet, "at least 0 jets pt>30 GeV");
	cut_bjet *bjet = new cut_bjet(1, 40, 2.5);
	basic_cuts.add_cut(bjet, "at least 1 b-jet pt>40 GeV");

	// apply cuts and extract efficiencies
	double eff_basic = thth_tztz.reduce_sample(basic_cuts); // require: 2 osl which reconstruct a Z, pT(Z), eta(Z), HT, nj, nb, pT(b) cuts
	double eff_fatjpt = thth_tztz.require_fatjet_pt(200, 1); // require at least 1 fatjet with pT>200 GeV
	double eff_ttag = thth_tztz.require_top_tagged(1); // require at least 1 fatjet to be HEP Top-Tagged
	double eff_ptT = cut_ptT(thth_tztz, 250); // require pT(t)>250 GeV
	double eff_deltaRtb = cut_deltaRtb(thth_tztz, 0.8); // require deltaR(t, b)<0.8 for at least one b-jet

	unsigned int remaining_events = thth_tztz.map_lhco_taggedJets.size();
	cout << "\nFatjet pT>200 GeV Efficiency: " << setprecision(4) << 100 * eff_fatjpt << " %" << endl;
	cout << "\nTop Tagging Efficiency: " << setprecision(4) << 100 * eff_ttag << " %" << endl;
	cout << "\npT(t)>250 GeV Efficiency: " << setprecision(4) << 100 * eff_ptT << " %" << endl;
	cout << "\nDeltaR(t,b)<0.8 Efficiency: " << setprecision(4) << 100 * eff_deltaRtb << " %" << endl;
	cout << "\nRemaining events: " << remaining_events << endl;
	cout << "\nSignal Efficiency: " << setprecision(4) << 100 * eff_basic * eff_fatjpt * eff_ttag * eff_ptT * eff_deltaRtb << " %" << endl;
	
	// identify top partner constituents
	vector<event*> signal_reconstructed;
	get_top_partner_constituents(thth_tztz.events(), thth_tztz.fatjets(), signal_reconstructed);
	
	// calculate cross section and store
	ofstream ofs_txt;
	ofs_txt.open(output_xsec.c_str());
	ofs_txt << input_xsec << "\t" << eff_basic * eff_fatjpt * eff_ttag * eff_ptT * eff_deltaRtb << "\t" << input_xsec * eff_basic * eff_fatjpt * eff_ttag * eff_ptT * eff_deltaRtb << endl;
	ofs_txt.close();
	
	// store lhco events
	write_lhco(signal_reconstructed, output_lhco);
	
	// clear remaining pointers
	delete osl;
	delete ptZ;
	delete etaZ;
	delete ht;
	delete njet;
	delete bjet;
	
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
	cout << "Input LHE file: " << input_sig_lhe << endl;
	cout << "Input LHCO file: " << input_sig_lhco << endl;
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
	cout << "Results output LHCO: " << output_lhco << endl;
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
	const particle *candidate1 = nullptr;
	const particle *candidate2 = nullptr;

	double delta_r_max = 1.0;
	double mz = 91.1876;
	double mass_range = 10.;
	double mass_diff_min = 1000.;

	for (unsigned int i = 0; i < leptons.size()-1; ++i)
	{
		const particle *first_l = leptons[i];

		for (unsigned int j = i+1; j < leptons.size(); ++j)
		{
			const particle *second_l = leptons[j];

			// same flavour required
			if (second_l->type() != first_l->type())
				continue;
			
			// opposite charge required
			if (second_l->charge() + first_l->charge() != 0)
				continue;

			// delta_r < delta_r_max required
			if (delta_r(first_l, second_l) > delta_r_max)
				continue;

			// inv. mass required to be within 10 GeV from mZ
			double inv_mass = mass({first_l, second_l});
			if (abs(inv_mass - mz) > mass_range)
				continue;

			// select osl leptons with inv. mass closest to mZ
			if (abs(inv_mass - mz) < mass_diff_min)
			{
				candidate1 = first_l;
				candidate2 = second_l;
				mass_diff_min = abs(inv_mass - mz);
			}
				 
		} // end second loop over leptons
		
	} // end first loop over leptons

	// return lepton candidates
	vector<const particle*> l_candidates;
	l_candidates.push_back(candidate1);
	l_candidates.push_back(candidate2);
	return l_candidates;
}

// identify top candidate most back-to-back wrt the reconstructed Z
PseudoJet identify_candidate_top(const vector<PseudoJet> & fatjets, vector<const particle*> & leptons)
{
	// reconstruct ideal back-to-back top wrt Z-boson
	const double px_Z = leptons[0]->px() + leptons[1]->px();
	const double py_Z = leptons[0]->py() + leptons[1]->py();
	const double pz_Z = leptons[0]->pz() + leptons[1]->pz();
	const double pe_Z = leptons[0]->pe() + leptons[1]->pe();
	PseudoJet Ideal_T(-px_Z, -py_Z, -pz_Z, pe_Z);

	// identify top-tagged jets
	vector< PseudoJet > topjets;
	for (unsigned int i = 0; i < fatjets.size(); ++i)
	{
		if (fatjets[i].user_info<TagInfo>().top_tag())
			topjets.push_back(fatjets[i]);
	}

	// identify top candidate
	PseudoJet top_candidate;
	double delta_r_min = 1000.;
	for (unsigned int i = 0; i < topjets.size(); ++i)
	{
		// deltaR between top candidate and back-to-back ideal direction
		double deltaEta = topjets[i].eta() - Ideal_T.eta();
		double deltaPhi = abs(topjets[i].phi() - Ideal_T.phi());
		deltaPhi = min(deltaPhi, 8 * atan(1) - deltaPhi);
		double deltaR = sqrt( pow(deltaEta,2.0) + pow(deltaPhi,2.0) );

		// minimize deltaR
		if (deltaR < delta_r_min)
		{
			delta_r_min = deltaR;
			top_candidate = topjets[i];
		}
	}

	// return top candidate
	return top_candidate;
}

// cut: require pT(t)>pT_min
double cut_ptT(jet_analysis &analysis, double pt_min)
{
	map< event *, vector< PseudoJet > >::iterator it;
	map< event *, vector< PseudoJet > > original_map = analysis.map_lhco_taggedJets;
	map< event *, vector< PseudoJet > > reduced_map;
	int map_size = 0, passed = 0;
	double eff;

	// loop over events
	for (it=original_map.begin(); it!=original_map.end(); ++it)
	{
		map_size++;

		// extract lepton candidates
		event *ev = it->first;
		vector< const particle* > leptons;
		for (unsigned int j = 0; j < ev->size(); ++j)
		{
			if ((*ev)[j]->type() & ptype_lepton && (*ev)[j]->pt() > 25. && abs((*ev)[j]->eta()) < 2.5)
				leptons.push_back((*ev)[j]);
		}		
		vector<const particle*> l_candidates = identify_candidate_leptons(leptons);

		// identify top candidate		
		vector< PseudoJet > fatjets = it->second;
		PseudoJet top_candidate = identify_candidate_top(fatjets, l_candidates);

		// evaluate pT(t): the cut is passed if pT(t)>pT_min
		double pt_T = top_candidate.pt();
		if (pt_T > pt_min)
		{
			passed++;
			reduced_map.insert( make_pair(it->first,it->second) );
		}

	}

	// calculate efficiency
	eff = (map_size == 0 ? 0.0 : static_cast<double>(passed) / map_size);

	// resize the map_lhco_taggedJets
	analysis.map_lhco_taggedJets = reduced_map;

	return eff;
}

// cut: require eta(t)<eta_max
double cut_etaT(jet_analysis &analysis, double eta_max)
{
	map< event *, vector< PseudoJet > >::iterator it;
	map< event *, vector< PseudoJet > > original_map = analysis.map_lhco_taggedJets;
	map< event *, vector< PseudoJet > > reduced_map;
	int map_size = 0, passed = 0;
	double eff;

	// loop over events
	for (it=original_map.begin(); it!=original_map.end(); ++it)
	{
		map_size++;

		// extract lepton candidates
		event *ev = it->first;
		vector< const particle* > leptons;
		for (unsigned int j = 0; j < ev->size(); ++j)
		{
			if ((*ev)[j]->type() & ptype_lepton && (*ev)[j]->pt() > 25. && abs((*ev)[j]->eta()) < 2.5)
				leptons.push_back((*ev)[j]);
		}		
		vector<const particle*> l_candidates = identify_candidate_leptons(leptons);

		// identify top candidate		
		vector< PseudoJet > fatjets = it->second;
		PseudoJet top_candidate = identify_candidate_top(fatjets, l_candidates);

		// evaluate eta(t): the cut is passed if eta(t)<eta_max
		double eta_T = abs(top_candidate.eta());
		if (eta_T < eta_max)
		{
			passed++;
			reduced_map.insert( make_pair(it->first,it->second) );
		}

	}

	// calculate efficiency
	eff = (map_size == 0 ? 0.0 : static_cast<double>(passed) / map_size);

	// resize the map_lhco_taggedJets
	analysis.map_lhco_taggedJets = reduced_map;

	return eff;
}

// cut: require deltaR(t, b)<delta_r for at least one b-jet
double cut_deltaRtb(jet_analysis &analysis, double delta_r)
{
	map< event *, vector< PseudoJet > >::iterator it;
	map< event *, vector< PseudoJet > > original_map = analysis.map_lhco_taggedJets;
	map< event *, vector< PseudoJet > > reduced_map;
	int map_size = 0, passed = 0;
	double eff;

	// loop over events
	for (it=original_map.begin(); it!=original_map.end(); ++it)
	{
		map_size++;

		// extract lepton candidates
		event *ev = it->first;
		vector< const particle* > leptons;
		for (unsigned int j = 0; j < ev->size(); ++j)
		{
			if ((*ev)[j]->type() & ptype_lepton && (*ev)[j]->pt() > 25. && abs((*ev)[j]->eta()) < 2.5)
				leptons.push_back((*ev)[j]);
		}		
		vector<const particle*> l_candidates = identify_candidate_leptons(leptons);

		// identify top candidate		
		vector< PseudoJet > fatjets = it->second;
		PseudoJet top_candidate = identify_candidate_top(fatjets, l_candidates);
		double	top_eta = top_candidate.eta(), 
				top_phi = top_candidate.phi();

		// extract all visible b-jets
		vector<const particle*> bjets;
		for (unsigned int j = 0; j < ev->size(); ++j)
		{
			if ((*ev)[j]->bjet() != 0.0 && (*ev)[j]->pt() > 30. && abs((*ev)[j]->eta()) < 2.5)
				bjets.push_back((*ev)[j]);
		}

		// evaluate deltaR(t, b): the cut is passed if at least one b-jet is within a cone of delta_r wrt the reconstructed top
		for (unsigned int j = 0; j < bjets.size(); ++j)
		{
			double	b_eta = bjets[j]->eta(), 
					b_phi = bjets[j]->phi();

			double deltaEta = top_eta - b_eta;
			double deltaPhi = abs(top_phi - b_phi);
			deltaPhi = min(deltaPhi, 8 * atan(1) - deltaPhi);
			double deltaR = sqrt( pow(deltaEta,2.0) + pow(deltaPhi,2.0) );

			if (deltaR < delta_r)
			{
				passed++;
				reduced_map.insert( make_pair(it->first,it->second) );
				break;
			}		
		}

	}

	// calculate efficiency
	eff = (map_size == 0 ? 0.0 : static_cast<double>(passed) / map_size);

	// resize the map_lhco_taggedJets
	analysis.map_lhco_taggedJets = reduced_map;

	return eff;
}

// reconstruct top partner consistituents
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
			if ((*ev)[j]->type() & ptype_lepton && (*ev)[j]->pt() > 25. && abs((*ev)[j]->eta()) < 2.5)
				leptons.push_back((*ev)[j]);
		}		
		vector<const particle*> l_candidates = identify_candidate_leptons(leptons);
		for (unsigned int j = 0; j < l_candidates.size(); ++j)
		{
			int p_type = ptype_electron;
			double	p_eta 	 = l_candidates[j]->eta(), 
					p_phi 	 = l_candidates[j]->phi(), 
					p_pt 	 = l_candidates[j]->pt(), 
					p_m 	 = l_candidates[j]->mass(),
					p_charge = l_candidates[j]->charge();
			lhco *p = new lhco(p_type, p_eta, p_phi, p_pt, p_m, p_charge);
			newev->push_back(p);
		}

		// extract top candidate		
		vector< PseudoJet > fatjets = signal_fatjets[i];
		PseudoJet top_candidate = identify_candidate_top(fatjets, l_candidates);
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
