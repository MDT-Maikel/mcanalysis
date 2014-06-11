/* thth->tztz cutmap 
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
bool load_settings(const string &settings_file, string &input_sig_lhco, string &output_file_cutmap);
vector<const particle*> identify_candidate_leptons(const vector<const particle*> & leptons, const double RLL_max);

// basic cut: at least two opposite sign leptons need to be present, with invariant mass near the Z boson
class cut_2osl : public cut
{
public:
	cut_2osl(double pt, double eta, double RLL) : pt_min(pt), eta_max(eta), RLL_max(RLL) {}

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
		vector<const particle*> l_candidates = identify_candidate_leptons(leptons, RLL_max);
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
	double RLL_max;
};

// basic cut: pt of the reconstructed Z boson
class cut_ptZ : public cut
{
public:
	cut_ptZ(double pt, double RLL) : pt_min(pt), RLL_max(RLL) {}

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
		vector<const particle*> l_candidates = identify_candidate_leptons(leptons, RLL_max);
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
	double RLL_max;
};

// basic cut: eta of the reconstructed Z boson
class cut_etaZ : public cut
{
public:
	cut_etaZ(double eta, double RLL) : eta_max(eta), RLL_max(RLL) {}

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
		vector<const particle*> l_candidates = identify_candidate_leptons(leptons, RLL_max);
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
	double RLL_max;
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
	string settings_file = "tztag_cutmap.cmnd";
	if (argc >= 2) 
		settings_file = argv[1];

	// read the input settings from command file
	string input_lhco, output_cutmap;
	if (!load_settings(settings_file, input_lhco, output_cutmap))
		return EXIT_FAILURE;
	ofstream cutmap_table(output_cutmap.c_str());

	// load lhco events
	vector<event*> events;
	read_lhco(events, input_lhco);

	// cut values, efficiency variables, cutted samples
	int ptZ_cut, ht_cut, nj_cut, ptB_cut; 
	double RLL_cut, etaZ_cut;
	double tot_eff, RLL_eff, ptZ_eff, etaZ_eff, ht_eff, nj_eff, ptB_eff;
	vector<event*> RLL_data, ptZ_data, etaZ_data, ht_data, nj_data, ptB_data;

	// cutmap loggin variables
	int cutmap_done = 0, cutmap_total = 7 * 7 * 8 * 6 * 4 * 6, cutmap_logging = 7 * 7; 

	// Delat_R(LL) cut: 7 steps
	RLL_data = events;
	RLL_eff = 1;
	for (RLL_cut = 1.6; RLL_cut <= 0.8; RLL_cut -= 0.2)
  	{
		if (RLL_data.size() > 0)
		{
			cuts RLL_cuts;
			cut_2osl *osl = new cut_2osl(25, 2.5, RLL_cut);
			RLL_cuts.add_cut(osl);
			RLL_cuts.apply(RLL_data);
			RLL_eff *= RLL_cuts.efficiency();

			RLL_cuts.clear();
			delete osl;
		}

		// pT(Z) cut: 7 steps.
		ptZ_data = RLL_data;
		ptZ_eff = 1;
		for (ptZ_cut = 200; ptZ_cut <= 350; ptZ_cut += 25)
		{
			if (ptZ_data.size() > 0)
			{	
				cuts ptZ_cuts;
				cut_ptZ *ptZ = new cut_ptZ(ptZ_cut, RLL_cut);
				ptZ_cuts.add_cut(ptZ);
				ptZ_cuts.apply(ptZ_data);
				ptZ_eff *= ptZ_cuts.efficiency();

				ptZ_cuts.clear();
				delete ptZ;
			}				

			// eta(Z) cut: 8 steps.
			etaZ_data = ptZ_data;
			etaZ_eff = 1;
			for (etaZ_cut = 2.5; etaZ_cut >= 1.1; etaZ_cut -= 0.2)
			{
				if (etaZ_data.size() > 0)
				{
					cuts etaZ_cuts;
					cut_etaZ *etaZ = new cut_etaZ(etaZ_cut, RLL_cut);
					etaZ_cuts.add_cut(etaZ);
					etaZ_cuts.apply(etaZ_data);
					etaZ_eff *= etaZ_cuts.efficiency();

					etaZ_cuts.clear();
					delete etaZ;
				}

				// HT cut: 6 steps.
				ht_data = etaZ_data;
				ht_eff = 1;
				for (ht_cut = 400; ht_cut <= 900; ht_cut += 100)
				{
					if (ht_data.size() > 0)
					{
						cuts ht_cuts;
						cut_ht *ht = new cut_ht(ht_cut, ptype_jet, 30, 3.0);
						ht_cuts.add_cut(ht);
						ht_cuts.apply(ht_data);
						ht_eff *= ht_cuts.efficiency();

						ht_cuts.clear();
						delete ht;
					}
	
					// n_jets cut: 4 steps.
					nj_data = ht_data;
					nj_eff = 1;
					for (nj_cut = 0; nj_cut <=6; nj_cut += 2)
					{
						if (nj_data.size() > 0)
						{	

							cuts nj_cuts;
							cut_njet *njet = new cut_njet(nj_cut, 30, 3.0);
							nj_cuts.add_cut(njet);
							nj_cuts.apply(nj_data);
							nj_eff *= nj_cuts.efficiency();

							nj_cuts.clear();
							delete njet;		
						}

						// pT(B) cut: 6 steps.
						ptB_data = nj_data;
						ptB_eff = 1;
						for (ptB_cut = 40; ptB_cut <=140; ptB_cut += 20)
						{
							if (ptB_data.size() > 0)
							{	

								cuts ptB_cuts;
								cut_bjet *bjet = new cut_bjet(1, ptB_cut, 2.8);
								ptB_cuts.add_cut(bjet);
								ptB_cuts.apply(ptB_data);
								ptB_eff *= ptB_cuts.efficiency();

								ptB_cuts.clear();
								delete bjet;	
							}

							// Calculate combined efficiency and store it in the cutmap table.
							tot_eff = RLL_eff * ptZ_eff * etaZ_eff * ht_eff * nj_eff * ptB_eff;
							cutmap_table << RLL_cut << "\t" << RLL_eff << "\t" << ptZ_cut << "\t" << ptZ_eff << "\t" << etaZ_cut << "\t" << etaZ_eff << "\t";
							cutmap_table << ht_cut << "\t" << ht_eff << "\t" << nj_cut << "\t" << nj_eff << "\t" << ptB_cut << "\t" << ptB_eff << "\t" << tot_eff << endl;

							// Print progress in steps of some total.
							cutmap_done++;
							if (cutmap_done % (cutmap_total / cutmap_logging) == 0)
								cout << "Cutmap in progress: done " << cutmap_done << "/" << cutmap_total << " steps" <<endl;

						}
					}
				}
			}			
		}
	}

	// Close the write-to text stream.
	cutmap_table.close();
	
	// clear remaining event pointers
	delete_events(events);
	delete_events(RLL_data);
	delete_events(ptZ_data);
	delete_events(etaZ_data);
	delete_events(ht_data);
	delete_events(nj_data);
	delete_events(ptB_data);
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Program completed in " << duration << " seconds." << endl;
	cout << "=====================================================================" << endl;
	
	// finished the program
	return EXIT_SUCCESS;
}

// load settings
bool load_settings(const string &settings_file, string &input_sig_lhco, string &output_file_cutmap)
{
	// read settings
	input_sig_lhco = read_settings<string>(settings_file, static_cast<string>("INPUT_LHCO"));
	output_file_cutmap = read_settings<string>(settings_file, static_cast<string>("OUTPUT_CUTMAP"));
	
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded settings from " << settings_file << ":" << endl;
	cout << "Input LHCO file: " << input_sig_lhco << endl;
	cout << "Output Cutmap: " << output_file_cutmap << endl;
	cout << "################################################################################" << endl;
	return true;
}

// identify lepton candidates to reconstruct the Z boson
vector<const particle*> identify_candidate_leptons(const vector<const particle*> & leptons, const double RLL_max)
{
	const particle *candidate1 = nullptr;
	const particle *candidate2 = nullptr;

	double delta_r_max = RLL_max;
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
