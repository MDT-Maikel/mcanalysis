/* GRES Basic
 *
 * Command line program that applies the basic cuts for the GRES Hunter algorithm 
 * to a lhco file without header. The resulting lhco file is then stored in a new
 * file. 
 *
 * This algorithm implements the following cuts (in this order):
 *  - Number of jets
 *  - pT of X leading jets
 *  - Lepton veto
 *  - Invariant mass combinations
 *
 * ./gres_basic <input file name> <output file name>
 * 
*/

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include <TH1.h> 

#include "utility/utility.h"
#include "event/event.h"
#include "cuts/cuts.h"
#include "plot/plot.h"
#include "bumphunter/bumphunter.h"
#include "gres/gres_cuts.h"
#include "gres/gres_hunter.h"


using namespace std;
using namespace analysis;
using namespace boost;
using namespace boost::filesystem;
using namespace boost::math;


// function prototypes
bool load_settings(vector<string> &input_files_bkg, vector<unsigned int> &input_count_bkg, vector<string> &output_files_bkg, vector<string> &names_bkg, vector<double> &sigmas_bkg, string &input_file_sig, string &output_file_sig, string &name_sig, double &sigma_sig, string &output_folder, string &output_txt, int &nr_jets, double &nr_jets_ptcut, double &nr_jets_etacut, vector<int> &topology, vector<double> &topology_cuts, double &luminosity);

void load_events(vector< vector<event*> > &events_bkgs, vector<event*> &events_sig, const vector<string> &input_files_bkg, const vector<unsigned int> &input_count_bkg, const string &input_file_sig);

void cut_basic(vector< vector<event*> > &events_bkgs, const vector<double> &sigmas_bkg, vector<double> &efficiencies_bkg, vector<event*> &events_sig, const double &sigma_sig, double &efficiency_sig, ofstream &ofs_txt, const double &luminosity, const int &nr_jets, const double &nr_jets_ptcut, const double &nr_jets_etacut, const vector<int> &topology, const vector<double> &topology_cuts);

void plot_invmass(const vector< vector<event*> > &events_bkgs, const vector<string> &names_bkg, const vector<double> &sigmas_bkg, const vector<double> &efficiencies_bkg, const vector<event*> &events_sig, const string &name_sig, const double &sigma_sig, const double &efficiency_sig, const string &output_folder, const double &luminosity, const vector<int> &topology, const vector<double> &topology_cuts);

void bump_hunter(const vector< vector<event*> > &events_bkgs, const vector<string> &names_bkg, const vector<double> &sigmas_bkg, const vector<double> &efficiencies_bkg, const vector<event*> &events_sig, const string &name_sig, const double &sigma_sig, const double &efficiency_sig, const string &output_folder, const double &luminosity, const vector<int> &topology, const vector<double> &topology_cuts);


class plot_masscomb : public plot_implementation
{
public:
	plot_masscomb(int t, const std::vector<int> & c) : type(t), comb(c) {}

	double operator() (const event *ev) 
	{ 
		ev->mass(type, comb);
	}

private:
	int type;
	std::vector<int> comb;
};

// main program
int main(int argc, const char* argv[])
{
	// background variables
	vector<string> input_files_bkg, output_files_bkg, names_bkg;
	vector<unsigned int> input_count_bkg;
	vector<double> sigmas_bkg, efficiencies_bkg;
	// signal variables
	string input_file_sig, output_file_sig, name_sig;
	double sigma_sig, efficiency_sig;
	// output variables
	string output_folder, output_txt;
	// default settings for the basic analysis, overloaded by settings loading
	int nr_jets = 3; double nr_jets_ptcut = 100, nr_jets_etacut = 2.8;
	// default settings for the topology and topology cuts
	// an topology mass cut on a 1-resonance corresponds to a jet pT cut
	// otherwise it corresponds to an invariant mass cut
	vector<int> topology; topology.push_back(2); topology.push_back(1);
	vector<double> topology_cuts; topology_cuts.push_back(800); topology_cuts.push_back(500);
	double luminosity = 1;
	
	// load the settings from the gres_basic_card.txt file, stop program if failed
	if (!load_settings(input_files_bkg, input_count_bkg, output_files_bkg, names_bkg, sigmas_bkg, input_file_sig, output_file_sig, name_sig, sigma_sig, output_folder, output_txt, nr_jets, nr_jets_ptcut, nr_jets_etacut, topology, topology_cuts, luminosity))
		return 1;

	// load the events
	vector< vector<event*> > events_bkgs; 
	vector<event*> events_sig;
	load_events(events_bkgs, events_sig, input_files_bkg, input_count_bkg, input_file_sig);

	// open a text stream for efficiency storing
	ofstream ofs_txt;
	ofs_txt.open((output_folder + output_txt).c_str());

	// perform the basic cuts on background and signal
	cut_basic(events_bkgs, sigmas_bkg, efficiencies_bkg, events_sig, sigma_sig, efficiency_sig, ofs_txt, luminosity, nr_jets, nr_jets_ptcut, nr_jets_etacut, topology, topology_cuts);

	// plot the invariant mass combination of signal compared to background
	plot_invmass(events_bkgs, names_bkg, sigmas_bkg, efficiencies_bkg, events_sig, name_sig, sigma_sig, efficiency_sig, output_folder, luminosity, topology, topology_cuts);

	// perform bumphunter algorithm based on the topology
	bump_hunter(events_bkgs, names_bkg, sigmas_bkg, efficiencies_bkg, events_sig, name_sig, sigma_sig, efficiency_sig, output_folder, luminosity, topology, topology_cuts);

	// write the remaining events to file
	for (unsigned int i = 0; i < events_bkgs.size(); ++i)
	{
		write_lhco(events_bkgs[i], output_folder + output_files_bkg[i]);
	}
	write_lhco(events_sig, output_folder + output_file_sig);

	// close the write-to text stream
	ofs_txt.close();
}


// loads the settings for the gres_basic from the gres_basic_card.txt
bool load_settings(vector<string> &input_files_bkg, vector<unsigned int> &input_count_bkg, vector<string> &output_files_bkg, vector<string> &names_bkg, vector<double> &sigmas_bkg, string &input_file_sig, string &output_file_sig, string &name_sig, double &sigma_sig, string &output_folder, string &output_txt, int &nr_jets, double &nr_jets_ptcut, double &nr_jets_etacut, vector<int> &topology, vector<double> &topology_cuts, double &luminosity)
{
	// load the settings from the gres_basic_card.txt, also initialize loading variables
	string settings_file = "../../files/gres/gres_basic_card.txt";

	// variables to store the settings
	stringstream topology_stream, top_cuts_stream;

	// read background settings
	input_files_bkg.clear();
	input_files_bkg = read_settings_list<string>(settings_file, static_cast<string>("BKG_INPUT_FILES"));
	input_count_bkg.clear();
	input_count_bkg = read_settings_list<unsigned int>(settings_file, static_cast<string>("BKG_INPUT_COUNT"));
	output_files_bkg.clear();
	output_files_bkg = read_settings_list<string>(settings_file, static_cast<string>("BKG_OUTPUT_FILES"));
	names_bkg.clear();
	names_bkg = read_settings_list<string>(settings_file, static_cast<string>("BKG_NAMES"));
	sigmas_bkg.clear();
	sigmas_bkg =  read_settings_list<double>(settings_file, static_cast<string>("BKG_SIGMAS"));

	// read signal settings
	input_file_sig = read_settings<string>(settings_file, static_cast<string>("SIG_INPUT_FILE"));
	output_file_sig = read_settings<string>(settings_file, static_cast<string>("SIG_OUTPUT_FILE"));
	name_sig = read_settings<string>(settings_file, static_cast<string>("SIG_NAME"));
	sigma_sig = read_settings<double>(settings_file, static_cast<string>("SIG_SIGMA"));

	// read general settings
	output_folder = read_settings<string>(settings_file, static_cast<string>("OUTPUT_FOLDER"));
	output_txt = read_settings<string>(settings_file, static_cast<string>("OUTPUT_TXT"));
	luminosity = read_settings<double>(settings_file, static_cast<string>("LUMINOSITY"));

	// read jets settings
	nr_jets = read_settings<int>(settings_file, static_cast<string>("NUMBER_JETS"));
	nr_jets_ptcut = read_settings<double>(settings_file, static_cast<string>("NUMBER_JETS_PT"));
	nr_jets_etacut = read_settings<double>(settings_file, static_cast<string>("NUMBER_JETS_ETA"));

	// read topology settings
	topology.clear();
	topology = read_settings_list<int>(settings_file, static_cast<string>("TOPOLOGY"));
	topology_cuts.clear();
	topology_cuts = read_settings_list<double>(settings_file, static_cast<string>("TOPOLOGY_CUTS"));

	// check the consistency of the settings
	if (input_files_bkg.size() != output_files_bkg.size())
	{
		cout << "ERROR: Entries BKG_INPUT_FILES and BKG_OUTPUT_FILES in gres_basic_card.txt have different sizes." << endl;
		return false;
	}
	if (topology.size() != topology_cuts.size())
	{
		cout << "ERROR: Entries TOPOLOGY and TOPOLOGY_CUTS in gres_basic_card.txt differ in size." << endl;
		return false;
	}
	int topology_content = 0;
	for (unsigned int i = 0; i < topology.size(); i++)
		topology_content += topology[i];
	if (topology_content != nr_jets)
	{
		cout << "ERROR: Resonance combination specified in TOPOLOGY in gres_basic_card.txt does not match NUMBER_JETS." << endl;
		//return false;
	}

	// create output directory if it does not exist yet
	if (!is_directory(output_folder))
		create_directory(output_folder);

	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded gres_basic_card.txt and start analysis with following settings:" << endl;
	cout << "LHCO background input files: {";
	for (unsigned int i = 0; i < input_files_bkg.size(); ++i)
	{	
		cout << input_files_bkg[i];
		if (i != input_files_bkg.size() - 1)
			cout << ", ";
	}
	cout << "}" << endl;
	cout << "LHCO background output files: {";
	for (unsigned int i = 0; i < output_files_bkg.size(); ++i)
	{	
		cout << output_files_bkg[i];
		if (i != output_files_bkg.size() - 1)
			cout << ", ";
	}
	cout << "}" << endl;
	cout << "Background cross sections: {";
	for (unsigned int i = 0; i < sigmas_bkg.size(); ++i)
	{	
		cout << sigmas_bkg[i];
		if (i != sigmas_bkg.size() - 1)
			cout << ", ";
	}
	cout << "} pb" << endl;
	cout << "LHCO signal input file: " << input_file_sig << endl;
	cout << "LHCO signal output file: " << output_file_sig << endl;
	cout << "Signal cross section: " << sigma_sig << " pb" << endl;
	cout << "Results output folder: " << output_folder << endl;
	cout << "Efficiency output file: " << output_txt << endl;
	cout << "Do analysis for " << nr_jets << " jets, with pT > " << nr_jets_ptcut << " GeV and eta < " << nr_jets_etacut << endl;
	cout << "Considering the following topology: {";
	for (unsigned int i = 0; i < topology.size(); ++i)
	{	
		cout << topology[i];
		if (i != topology.size() - 1)
			cout << ", ";
	}
	cout << "} resonance structure." << endl;
	cout << "Performing the following topology cuts: {";
	for (unsigned int i = 0; i < topology_cuts.size(); i++)
	{	
		cout << topology_cuts[i];
		if (i != topology_cuts.size() - 1)
			cout << ", ";
	}
	cout << "} GeV." << endl;	
	cout << "################################################################################" << endl;
	return true;
}

void load_events(vector< vector<event*> > &events_bkgs, vector<event*> &events_sig, const vector<string> &input_files_bkg, const vector<unsigned int> &input_count_bkg, const string &input_file_sig)
{
	// logging
	unsigned int nr_samples = 0;
	cout << "Loading event samples ..." << endl;

	// load backgrounds	
	for (unsigned int i = 0; i < input_files_bkg.size(); ++i)
	{
		path input_file_bkg(input_files_bkg[i]);
		vector<path> files;
		if (is_directory(input_file_bkg))
			files = get_files(input_file_bkg.string(), ".lhco.gz", true);
		else
			files.push_back(input_file_bkg);

		vector<event*> events_bkg;
		unsigned int count_bkg = input_count_bkg[i];
		for (unsigned int j = 0; j < min(static_cast<unsigned int>(files.size()), count_bkg); ++j)
		{
			// log which file is being read			
			cout << "... reading " << files[j] << " events ..." << "\r" << flush;
			nr_samples++;

			read_lhco(events_bkg, files[j]);
		}
		events_bkgs.push_back(events_bkg);
	}

	// load signal
	cout << "... reading sample " << input_file_sig << " ..." << "\r" << flush;
	nr_samples++;
	read_lhco(events_sig, input_file_sig);
	
	// finished, log this
	cout << "... finished loading " << nr_samples << " event samples                     " << endl;
	cout << "################################################################################" << endl;
}

void cut_basic(vector< vector<event*> > &events_bkgs, const vector<double> &sigmas_bkg, vector<double> &efficiencies_bkg, vector<event*> &events_sig, const double &sigma_sig, double &efficiency_sig, ofstream &ofs_txt, const double &luminosity, const int &nr_jets, const double &nr_jets_ptcut, const double &nr_jets_etacut, const vector<int> &topology, const vector<double> &topology_cuts)
{
	// set up the different cuts
	cuts basic_cuts;

	// cut: amount of jets with pt and eta
	cut_pt *cut_nrjets; 
	cut_nrjets = new cut_pt(nr_jets_ptcut, ptype_jet, nr_jets, nr_jets_etacut);
	basic_cuts.add_cut(cut_nrjets, "nr jets");

	// cut: specific pt cuts on the leading jets
	vector<double> pt_cuts;
	for (unsigned int i = 0; i < topology.size(); i++)
		if (topology[i] == 1)
			pt_cuts.push_back(topology_cuts[i]);

	vector<cut_pt*> jet_cuts;
	for (unsigned int i = 0; i < pt_cuts.size(); i++)
	{
		cut_pt *cut_leadjets; 
		cut_leadjets = new cut_pt(pt_cuts[i], ptype_jet, i + 1, 2.8);
		jet_cuts.push_back(cut_leadjets);
		basic_cuts.add_cut(cut_leadjets, "pt jet n");
	}

	// cut: lepton vetos
	cut_veto *veto_el; veto_el = new cut_veto(ptype_electron, 10, 2.5);
	cut_veto *veto_mu; veto_mu = new cut_veto(ptype_muon, 10, 2.5);
	cut_veto *veto_tau; veto_tau = new cut_veto(ptype_tau, 20, 2.5);
	basic_cuts.add_cut(veto_el, "electron veto");
	basic_cuts.add_cut(veto_mu, "muon veto");
	basic_cuts.add_cut(veto_tau, "tau veto");

	// cut: invariant mass cuts
	cut_gres *cut_mass;
	cut_mass = new cut_gres();
	cut_mass->set_topology(topology);
	cut_mass->set_top_cuts(topology_cuts);
	cut_mass->init();
	basic_cuts.add_cut(cut_mass, "mass cut");

	// apply the basic cuts to the background and print the results
	for (unsigned int i = 0; i < events_bkgs.size(); ++i)
	{
		basic_cuts.clear();		
		basic_cuts.apply(events_bkgs[i]);
		basic_cuts.write(cout);
		ofs_txt << "Analyzing background " << i << ":" << endl;
		basic_cuts.write(ofs_txt);
		efficiencies_bkg.push_back(basic_cuts.efficiency());
		cout << "Corresponding number of events: " << luminosity * efficiencies_bkg[i] * sigmas_bkg[i] << endl;
	}
	
	// apply the basic cuts to the signal and print the results
	basic_cuts.clear();
	basic_cuts.apply(events_sig);
	basic_cuts.write(cout);
	ofs_txt << "Analyzing signal:" << endl;
	basic_cuts.write(ofs_txt);
	efficiency_sig = basic_cuts.efficiency();
	cout << "Corresponding number of events: " << luminosity * efficiency_sig * sigma_sig << endl;

	// delete cuts
	delete cut_nrjets;
	for (unsigned int i = 0; i < jet_cuts.size(); i++)
		delete jet_cuts[i];
	jet_cuts.clear();
	delete veto_el; delete veto_mu; delete veto_tau;
	delete cut_mass;
}

void plot_invmass(const vector< vector<event*> > &events_bkgs, const vector<string> &names_bkg, const vector<double> &sigmas_bkg, const vector<double> &efficiencies_bkg, const vector<event*> &events_sig, const string &name_sig, const double &sigma_sig, const double &efficiency_sig, const string &output_folder, const double &luminosity, const vector<int> &topology, const vector<double> &topology_cuts)
{
	// an invariant mass spectrum for both background and signal and for each of the resonance and invariant mass combinations is needed
	// however resonances of the same size don't need to be done twice, therefore first reduce the topology list to only distinctly sized resonances
	unsigned int nr_jets = 0;
	for (unsigned int i = 0; i < topology.size(); i++)
		nr_jets += topology[i];

	vector<int> red_topology;
	for (unsigned int i = 0; i < topology.size(); i++)
	{
		int resonance = topology[i];
		bool add = (resonance > 1);
		for (unsigned int j = 0; j < red_topology.size(); j++)
			if (red_topology[j] == resonance)
				add = false;
		if (add)
			red_topology.push_back(resonance);	
	}

	// log plotting for which topology
	cout << "################################################################################" << endl;
	cout << "Plotting invariant mass combinations for topology: {";
	for (unsigned int i = 0; i < red_topology.size(); i++)
	{	
		cout << red_topology[i];
		if (i != red_topology.size() - 1)
			cout << ", ";
	}
	cout << "} resonance structure." << endl;

	// produce all the needed invariant mass spectra, these consist of all combinations nr_jets over resonance size
	for (unsigned int i = 0; i < red_topology.size(); i++)
	{
		int resonance = red_topology[i];
		unsigned int nr_combs = static_cast<unsigned int>(binomial_coefficient<double>(nr_jets, resonance)); // TODO: Is this failsafe?
		vector<vector<int> > combinations;
		vector<int> start_comb;
		for (unsigned int j = resonance; j > 0; j--)
			start_comb.push_back(j);
		// produce all combinations which lead to the resonance
		for (unsigned int j = 0; j < nr_combs; j++)
		{
			combinations.push_back(start_comb);
			start_comb[0]++;
			for (unsigned int k = 0; k < start_comb.size() - 1; k++)
			{
				if (static_cast<unsigned int>(start_comb[k]) > nr_jets - k)
				{
					start_comb[k + 1]++;
					for (int l = k; l >= 0; l--)
						start_comb[l] = start_comb[l + 1] + 1;
				}
			}
		}

		// make the invariant mass spectrum plot for each of the combinations
		for (unsigned int j = 0; j < combinations.size(); j++)
		{
			// binlist name
			vector<int> comb = combinations[j];
			string file_name = "njet" + lexical_cast<string>(nr_jets) + "_comb";
			for (unsigned int k = 0; k < comb.size(); k++)
	 			file_name += lexical_cast<string>(comb[k]);	
	 	
	 		// mass plot	
	 		plot_masscomb *ftor_mass = new plot_masscomb(ptype_jet, comb); 
	 		plot pmass(file_name, output_folder);
	 		pmass.set_normalized(true);
			pmass.set_bins(50, 0, 3000);
			
			// add bkg and signal samples
			for (unsigned int k = 0; k < events_bkgs.size(); ++k)
				pmass.add_sample(events_bkgs[k], ftor_mass, names_bkg[k], luminosity * efficiencies_bkg[k] * sigmas_bkg[k]);
			pmass.add_sample(events_sig, ftor_mass, name_sig, luminosity * efficiency_sig * sigma_sig);

			// make normal, stacked and stacked/no-log-scale plot
			pmass.run();
			pmass.set_name(file_name + "_stacked");
			pmass.set_stacked(true);
			pmass.run();
			pmass.set_name(file_name + "_nolog");
			pmass.set_logy(false);
			pmass.run();
		}
	}
	cout << "################################################################################" << endl;
}

void bump_hunter(const vector< vector<event*> > &events_bkgs, const vector<string> &names_bkg, const vector<double> &sigmas_bkg, const vector<double> &efficiencies_bkg, const vector<event*> &events_sig, const string &name_sig, const double &sigma_sig, const double &efficiency_sig, const string &output_folder, const double &luminosity, const vector<int> &topology, const vector<double> &topology_cuts)
{
	cout << "################################################################################" << endl;
	
	// create a bumphunter analysis
	greshunter hunt;
	hunt.set_folder(output_folder);
	hunt.set_topology(topology);
	hunt.set_top_cuts(topology_cuts);
	
	// add event samples for background and signal
	for (unsigned int i = 0; i < events_bkgs.size(); i++)
	{
		vector<event*> events = events_bkgs[i];
		hunt.add_background(events, luminosity * efficiencies_bkg[i] * sigmas_bkg[i]);
	}
	hunt.set_signal(events_sig, luminosity * efficiency_sig * sigma_sig);
	
	// run the bumphunter analysis
	hunt.run();
	
	cout << "################################################################################" << endl;
	
	
	
	/*
	vector< vector<double> > mass_bkg;
	vector<double> mass_sig;


	// get background mass
	for (unsigned int i = 0; i < events_bkgs.size(); i++)
	{
		vector<event*> events = events_bkgs[i];
		vector<double> mass;
		for (unsigned int j = 0; j < events.size(); j++)
		{
			event *ev = events[j];
			mass.push_back(ev->mass(particle::type_jet, {2,3}));
		}
		mass_bkg.push_back(mass);	
	}

	// get signal mass
	for (unsigned int i = 0; i < events_sig.size(); i++)
	{
		event *ev = events_sig[i];
		mass_sig.push_back(ev->mass(particle::type_jet, {2,3}));
	}

	// histograms
	TH1* hist_bkg;
	hist_bkg = new TH1F("", "bkg", 100, 0, 4000);
	TH1* hist_sig;
	hist_sig = new TH1F("", "sig", 100, 0, 4000);

	// fill the histograms
	for (unsigned int i = 0; i < mass_bkg.size(); i++)
	{
		vector<double> mass = mass_bkg[i];
		for (unsigned int j = 0; j < mass_sig.size(); j++)
		{
			hist_bkg->Fill(mass[j], luminosity * efficiencies_bkg[i] * sigmas_bkg[i]);
			hist_sig->Fill(mass[j], luminosity * efficiencies_bkg[i] * sigmas_bkg[i]);
		}
	}
	for (unsigned int i = 0; i < mass_sig.size(); i++)
			hist_sig->Fill(mass_sig[i], luminosity * efficiency_sig * sigma_sig);
			
	// start bumphunter
	TBumpHunter b(hist_bkg, hist_sig);
	b.SetNPseudoExperiments(1000);
	b.SetBinModel(1);
	b.SetSearchRegion(800,3000);
	b.Run();*/

}


