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

#include <TF1.h>
#include <TH1.h>

#include "particle/lhco.h"
#include "particle/particle.h"
#include "event/event.h"
#include "utility/utility.h"
#include "cuts/cuts.h"
#include "plot/plot.h"
#include "plot/plot2d.h"
#include "jet_analysis/jet_analysis.h"
#include "bumphunter/bumphunter.h"

using namespace std;
using namespace boost;
using namespace boost::filesystem;
using namespace Pythia8;
using namespace fastjet;
using namespace analysis;

// function prototypes
bool load_settings_mcinput(const string &settings_file, vector<string> &bkg_lhco, vector<double> &bkg_xsec, vector<double> &bkg_kfac, double &bkg_scal, string &sig_lhco, double &sig_xsec, double &sig_kfac, double &sig_scal);
bool load_settings_general(const string &settings_file, string &output_folder, double &luminosity);
double get_thmass(const event *ev);

// main program: may have one argument
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;
	
	// take single argument specifying settings file if available
	string settings_file = "tztag_hunt.cmnd";
	if (argc >= 2) 
		settings_file = argv[1];
		
	// read the mcinput settings from command file	
	vector<string> bkg_lhco;
	vector<double> bkg_xsec;
	vector<double> bkg_kfac;
	double bkg_scal;
	string sig_lhco;
	double sig_xsec;
	double sig_kfac;
	double sig_scal;
	if (!load_settings_mcinput(settings_file, bkg_lhco, bkg_xsec, bkg_kfac, bkg_scal, sig_lhco, sig_xsec, sig_kfac, sig_scal))
		return EXIT_FAILURE;
	// read the general settings from command file
	string output_folder = "output/";
	double luminosity = 1;
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

	// define histrograms
	int nbins = 20;
	double hmin = 500.;
	double hmax = 1500.;

	TH1F* hist_b = new TH1F("","bkg",nbins,hmin,hmax);
	TH1F* hist_sb = new TH1F("","sig+bkg",nbins,hmin,hmax);

	// fill the histograms hist_b and hist_sb with background samples
	for (unsigned int i = 0; i < bkg_evts.size(); ++i)
	{
		vector<event*> events = bkg_evts[i];
		double xsec = bkg_xsec[i]*bkg_kfac[i]*bkg_scal;
		double weight = xsec * luminosity / events.size();
		
		for (unsigned int j = 0; j < events.size(); j++)
		{
			double val = get_thmass(events[j]);
			hist_b->Fill(val, weight);
			hist_sb->Fill(val, weight);
		}
	}

	// fill the histogram hist_sb with signal sample
	vector<event*> events = sig_evts;
	double xsec = sig_xsec*sig_kfac*sig_scal;
	double weight = xsec * luminosity / events.size();
	
	for (unsigned int j = 0; j < events.size(); j++)
	{
		double val = get_thmass(events[j]);
		hist_sb->Fill(val, weight);
	}

	// // for testing
	// for (int i = 1; i <= 15; ++i)
	// {
	// 	double bin = hist_sb->GetBinCenter(i);
	// 	double sb_val = hist_sb->GetBinContent(i);
	// 	double b_val = hist_b->GetBinContent(i);
	// 	cout << "hist_sb bin " << bin << " content: " << sb_val << endl;
	// 	cout << "hist_b bin " << bin << " content: " << b_val << endl << endl;
	// }
	// gStyle->SetOptStat(0);
	// TCanvas* c1 = new TCanvas("","");
	// c1->SetTicks(1,1);
	// c1->SetBottomMargin(0.2);
	// c1->SetLeftMargin(0.2);
	// c1->SetLogy(1);
	// hist_sb->Draw("");
	// c1->Print("output/test.png");
	
	// create a bumphunter instance and do basic settings
	bumphunter hunt(hist_b, hist_sb);
	hunt.set_folder(output_folder);
	hunt.set_name("tztag_hunt_poisson");
	hunt.SetNPseudoExperiments(1.0e+8);
	hunt.SetBinModel(bumphunter::BUMP_POISSON);
	hunt.SetSearchRegion(800, 1500);
	hunt.SetMinWindowSize(1);
	hunt.SetMaxWindowSize(3);
	hunt.SetWindowStepSize(1);
	
	// run bumphunter analysis
	hunt.run();
	double sigma_poisson = hunt.get_global_sigma();
	
	// run bumphunter analysis a second time with poisson-gamma bin model
	hunt.SetBinModel(bumphunter::BUMP_POISSON_GAMMA);
	hunt.set_name("tztag_hunt_poissongamma");
	hunt.run();
	double sigma_poissongamma = hunt.get_global_sigma();

	// clear remaining pointers
	delete hist_b;
	delete hist_sb;
	
	// clear remaining event pointers
	for (unsigned int i = 0; i < bkg_evts.size(); ++i)
		delete_events(bkg_evts[i]);
	delete_events(sig_evts);
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Program completed in " << duration << " seconds." << endl;
	cout << "Poisson bin model significance : " << max(0.0,sigma_poisson) << endl;
	cout << "PoissonGamma bin model significance : " << max(0.0,sigma_poissongamma) << endl;
	cout << "=====================================================================" << endl;
	
	// finished the program
	return EXIT_SUCCESS;
}

// extract top partner mass from event
double get_thmass(const event *ev)
{
	return mass({ev->get(ptype_lepton, 1), ev->get(ptype_lepton, 2), ev->get(ptype_jet, 1)});
}

// load settings functions
bool load_settings_mcinput(const string &settings_file, vector<string> &bkg_lhco, vector<double> &bkg_xsec, vector<double> &bkg_kfac, double &bkg_scal, string &sig_lhco, double &sig_xsec, double &sig_kfac, double &sig_scal)
{
	// read input settings
	bkg_lhco.clear(); 
	bkg_xsec.clear();
	bkg_lhco = read_settings_list<string>(settings_file, static_cast<string>("BKG_LHCO"));
	bkg_xsec = read_settings_list<double>(settings_file, static_cast<string>("BKG_XSEC"));
	bkg_kfac = read_settings_list<double>(settings_file, static_cast<string>("BKG_KFAC"));
	bkg_scal = read_settings<double>(settings_file, static_cast<string>("BKG_SCAL"));
	sig_lhco = read_settings<string>(settings_file, static_cast<string>("SIG_LHCO"));
	sig_xsec = read_settings<double>(settings_file, static_cast<string>("SIG_XSEC"));
	sig_kfac = read_settings<double>(settings_file, static_cast<string>("SIG_KFAC"));
	sig_scal = read_settings<double>(settings_file, static_cast<string>("SIG_SCAL"));
	
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
	cout << "Background K-factors: {";
	for (unsigned int i = 0; i < bkg_kfac.size(); ++i)
	{	
		cout << bkg_kfac[i];
		if (i != bkg_kfac.size() - 1)
			cout << ", ";
	}
	cout << "}" << endl;
	cout << "Background rescaling factor: " << bkg_scal << endl;
	cout << "Signal LHCO: " << sig_lhco << endl;
	cout << "Signal xsec: " << sig_xsec << " pb" << endl;
	cout << "Signal K-factor: " << sig_kfac << endl;
	cout << "Signal rescaling factor: " << sig_scal << endl;
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

