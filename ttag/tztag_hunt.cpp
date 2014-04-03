/* thth->tztz hunt
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
bool load_settings_mcinput(const string &settings_file, vector<string> &bkg_lhco, vector<double> &bkg_xsec, string &sig_lhco, double &sig_xsec);
bool load_settings_general(const string &settings_file, string &output_folder, double &luminosity);

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

// deltaR(Z, tagged top)
class plot_deltarZt : public plot_default
{
public:
	plot_deltarZt() {}

	double operator() (const event *ev)
	{
		// extract eta and phi of leptons and top-jet
		double etaL1 = (ev->get(ptype_lepton, 1))->eta();
		double etaL2 = (ev->get(ptype_lepton, 2))->eta();
		double etaJ = (ev->get(ptype_jet, 1))->eta();

		double phiL1 = (ev->get(ptype_lepton, 1))->phi();
		double phiL2 = (ev->get(ptype_lepton, 2))->phi();
		double phiJ = (ev->get(ptype_jet, 1))->phi();

		// evaluate deltaR(l1,j)
		double deltaEta1 = etaJ - etaL1;
		double deltaPhi1 = abs(phiJ - phiL1);
		deltaPhi1 = min(deltaPhi1, 8 * atan(1) - deltaPhi1);
		double deltaR1   = sqrt( pow(deltaEta1,2.0) + pow(deltaPhi1,2.0) );

		// evaluate deltaR(l2,j)
		double deltaEta2 = etaJ - etaL2;
		double deltaPhi2 = abs(phiJ - phiL2);
		deltaPhi2 = min(deltaPhi2, 8 * atan(1) - deltaPhi2);
		double deltaR2   = sqrt(pow(deltaEta2, 2.0) + pow(deltaPhi2, 2.0));

		// return average deltaR
		double deltaR = (deltaR1+deltaR2)/2;
		return deltaR;
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

// main program: may have one argument
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;
	
	// take single argument specifying settings file if available
	string settings_file = "tztag_hunt_1top.cmnd";
	if (argc >= 2) 
		settings_file = argv[1];
		
	// read the mcinput settings from command file	
	vector<string> bkg_lhco;
	vector<double> bkg_xsec;
	string sig_lhco;
	double sig_xsec;
	if (!load_settings_mcinput(settings_file, bkg_lhco, bkg_xsec, sig_lhco, sig_xsec))
		return EXIT_FAILURE;
	// read the general settings from command file
	double luminosity = 1;
	string output_folder = "output/";
	if (!load_settings_general(settings_file, output_folder, luminosity))
		return EXIT_FAILURE;
	luminosity *= 1000;
	
	// make sure the output file's directory exists
	if (!is_directory(output_folder))
		create_directory(output_folder);
		
	// load the events
	vector<vector<event*> > bkg_evts;
	vector<event*> sig_evts;
	for (unsigned int i = 0; i < bkg_lhco.size(); ++i)
	{
		vector<event*> evts;
		read_lhco(evts, bkg_lhco[i]);
		bkg_evts.push_back(evts);
	}	
	read_lhco(sig_evts, sig_lhco);
	
	// plot lepton pair mass
	plot lmass("plot_leptonmass", output_folder);
	plot_leptonmass *lepton_mass = new plot_leptonmass();
	for (unsigned int i = 0; i < bkg_evts.size(); ++i)
	{
		double weight = bkg_xsec[i] * luminosity / bkg_evts[i].size();
		lmass.add_sample(bkg_evts[i], lepton_mass, "bkg" + lexical_cast<string>(i), weight);
	}
	lmass.add_sample(sig_evts, lepton_mass, "signal", sig_xsec * luminosity / sig_evts.size());
	lmass.run();	

	// plot deltaR(Z, tagged top)
	plot drZt("plot_deltarZt", output_folder);
	drZt.set_normalized(true);
	plot_deltarZt *deltarZt = new plot_deltarZt();
	for (unsigned int i = 0; i < bkg_evts.size(); ++i)
	{
		drZt.add_sample(bkg_evts[i], deltarZt, "bkg" + lexical_cast<string>(i));
	}
	drZt.add_sample(sig_evts, deltarZt, "signal");
	drZt.run();
	
	// plot top partner mass
	plot pmass("plot_thmass", output_folder);
	pmass.set_bins(50, 200, 2200);
	plot_thmass *th_mass = new plot_thmass();
	for (unsigned int i = 0; i < bkg_evts.size(); ++i)
	{
		double weight = bkg_xsec[i] * luminosity / bkg_evts[i].size();
		pmass.add_sample(bkg_evts[i], th_mass, "bkg" + lexical_cast<string>(i), weight);
	}
	pmass.add_sample(sig_evts, th_mass, "signal", sig_xsec * luminosity / sig_evts.size());
	pmass.run();	
	
	// clear remaining pointers
	delete lepton_mass;
	delete th_mass;
	
	// clear remaining event pointers
	for (unsigned int i = 0; i < bkg_evts.size(); ++i)
		delete_events(bkg_evts[i]);
	delete_events(sig_evts);
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Program completed in " << duration << " seconds." << endl;
	cout << "=====================================================================" << endl;
	
	// finished the program
	return EXIT_SUCCESS;
}

bool load_settings_mcinput(const string &settings_file, vector<string> &bkg_lhco, vector<double> &bkg_xsec, string &sig_lhco, double &sig_xsec)
{
	// read input settings
	bkg_lhco.clear(); 
	bkg_xsec.clear();
	bkg_lhco = read_settings_list<string>(settings_file, static_cast<string>("BKG_LHCO"));
	bkg_xsec = read_settings_list<double>(settings_file, static_cast<string>("BKG_XSEC"));
	sig_lhco = read_settings<string>(settings_file, static_cast<string>("SIG_LHCO"));
	sig_xsec = read_settings<double>(settings_file, static_cast<string>("SIG_XSEC"));
	
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded input settings from " << settings_file << ":" << endl;
	cout << "Background LHCO: {";
	for (unsigned int i = 0; i < bkg_lhco.size(); ++i)
	{	
		cout << bkg_lhco[i];
		if (i != bkg_lhco.size() - 1)
			cout << ", ";
	}
	cout << "}" << endl;
	cout << "Background xsec: {";
	for (unsigned int i = 0; i < bkg_xsec.size(); ++i)
	{	
		cout << bkg_xsec[i];
		if (i != bkg_xsec.size() - 1)
			cout << ", ";
	}
	cout << "} pb" << endl;
	cout << "Signal LHCO: " << sig_lhco << endl;
	cout << "Signal xsec: " << sig_xsec << endl;
	cout << "################################################################################" << endl;
	return true;
}

bool load_settings_general(const string &settings_file, string &output_folder, double &luminosity)
{
	// read input settings
	output_folder = read_settings<string>(settings_file, static_cast<string>("OUTPUT_FOLDER"));
	luminosity = read_settings<double>(settings_file, static_cast<string>("LUMINOSITY"));
	
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded input settings from " << settings_file << ":" << endl;
	cout << "Output folder: " << output_folder << endl;
	cout << "Luminosity (ifb): " << luminosity << endl;
	cout << "################################################################################" << endl;
	return true;
}

