/* thth->tztz cutmap findmax
 *
 * 
*/

#include <cstdio>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include "utility/utility.h"

using namespace std;
using namespace boost;
using namespace boost::filesystem;
using namespace analysis;

// function prototypes
bool load_settings(const string &settings_file, int &lumi, int &ns_min);
bool load_files(const string &settings_file, string &file_ttz, string &file_zjj, string &file_tz, string &file_tbz, string &file_signal);
bool load_xsec(const string &settings_file, double &xsec_ttz, double &xsec_zjj, double &xsec_tz, double &xsec_tbz, double &xsec_signal);

// auxiliary class CutFlow, where the cut values, S/B ratio and remaining signal events are stored
class CutFlow
{
public:
	int ptZ, ht, nj, ptB;
	double RLL, etaZ, ns, eff;

	CutFlow() {};
};

// main program: may have one argument
int main(int argc, const char* argv[])
{
	// take single argument specifying settings file if available
	string settings_file = "tztag_cutmap_findmax.cmnd";
	if (argc >= 2) 
		settings_file = argv[1];

	// read the input settings from command file
	int lumi, ns_min;
	string file_ttz, file_zjj, file_tz, file_tbz, file_signal;
	double xsec_ttz, xsec_zjj, xsec_tz, xsec_tbz, xsec_signal;
	if (!load_settings(settings_file, lumi, ns_min))
		return EXIT_FAILURE;
	if (!load_files(settings_file, file_ttz, file_zjj, file_tz, file_tbz, file_signal))
		return EXIT_FAILURE;
	if (!load_xsec(settings_file, xsec_ttz, xsec_zjj, xsec_tz, xsec_tbz, xsec_signal))
		return EXIT_FAILURE;

	// load the cutmaps into streams
	ifstream ifs_ttz, ifs_zjj, ifs_tz, ifs_tbz, ifs_signal;
	ifs_ttz.open(file_ttz.c_str(), ios::in);
	ifs_zjj.open(file_zjj.c_str(), ios::in);
	ifs_tz.open(file_tz.c_str(), ios::in);
	ifs_tbz.open(file_tbz.c_str(), ios::in);
	ifs_signal.open(file_signal.c_str(), ios::in);

	// loop over the signal cutmap assuming all cutmaps have the same length
	vector<CutFlow> cutmap;
	double RLL, ptZ, etaZ, ht, nj, ptB, dump;
	double eff_ttz, eff_zjj, eff_tz, eff_tbz, eff_signal;
	double ns, nb, efficiency;

	while (ifs_signal)
	{
		// read the efficiency data per background type and signal
		ifs_ttz >> RLL >> dump >> ptZ >> dump >> etaZ >> dump >> ht >> dump >> nj >> dump >> ptB >> dump >> eff_ttz;
		ifs_zjj >> RLL >> dump >> ptZ >> dump >> etaZ >> dump >> ht >> dump >> nj >> dump >> ptB >> dump >> eff_zjj;
		ifs_tz >> RLL >> dump >> ptZ >> dump >> etaZ >> dump >> ht >> dump >> nj >> dump >> ptB >> dump >> eff_tz;
		ifs_tbz >> RLL >> dump >> ptZ >> dump >> etaZ >> dump >> ht >> dump >> nj >> dump >> ptB >> dump >> eff_tbz;
		ifs_signal >> RLL >> dump >> ptZ >> dump >> etaZ >> dump >> ht >> dump >> nj >> dump >> ptB >> dump >> eff_signal;

		// multiply efficiencies with cross section and luminosity
		eff_ttz *= xsec_ttz*lumi;
		eff_zjj *= xsec_zjj*lumi;
		eff_tz *= xsec_tz*lumi;
		eff_tbz *= xsec_tbz*lumi;
		
		// Determine relative signal over background and store into a list
		nb = eff_ttz + eff_zjj + eff_tz + eff_tbz;
		ns = eff_signal*xsec_signal*lumi;
		if (nb == 0)
			efficiency = 0;
		else
    		efficiency = ns / nb;

  		CutFlow flow;
  		flow.RLL = RLL; flow.ptZ = ptZ; flow.etaZ = etaZ; flow.ht = ht; flow.nj = nj; flow.ptB = ptB; flow.ns = ns; flow.eff = efficiency;
  		cutmap.push_back(flow);
	}
	

	// Find maximum value point in cutmap, with at least ns_min remaining signal events
	double eff_max = 0;
	int i_max;

	for (int i = 0; i < cutmap.size(); i++)
	{
		if (cutmap[i].eff >= eff_max && cutmap[i].ns >= ns_min)
		{
			eff_max = cutmap[i].eff;
			i_max = i;
		}
	}

	// Print the results of the optimization to the screen.
	cout << "S/B maximized (" << cutmap[i_max].eff << ") for:" << endl;
	cout << "\t DeltaR(L,L) cut = " << cutmap[i_max].RLL << endl; 
	cout << "\t ptZ cut = " << cutmap[i_max].ptZ << " GeV" << endl;
	cout << "\t etaZ cut = " << cutmap[i_max].etaZ << endl;
	cout << "\t ht cut = " << cutmap[i_max].ht << " GeV" << endl;
	cout << "\t nj cut = " << cutmap[i_max].nj << endl;
	cout << "\t ptB cut = " << cutmap[i_max].ptB << " GeV" << endl;
	cout << "with remaining " << cutmap[i_max].ns << " signal events." << endl;

	// finished the program
	return EXIT_SUCCESS;
}

// load settings
bool load_settings(const string &settings_file, int &lumi, int &ns_min)
{
	// read settings
	lumi = read_settings<int>(settings_file, static_cast<string>("LUMINOSITY"));
	lumi *= 1000;
	ns_min = read_settings<int>(settings_file, static_cast<string>("NS_MIN"));
	
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded settings from " << settings_file << ":" << endl;
	cout << "Input luminosity [fb-1]: " << lumi/1000 << endl;
	cout << "Input ns_min: " << ns_min << endl;
	cout << "################################################################################" << endl;
	return true;
}

// load files
bool load_files(const string &settings_file, string &file_ttz, string &file_zjj, string &file_tz, string &file_tbz, string &file_signal)
{
	// read settings
	file_ttz = read_settings<string>(settings_file, static_cast<string>("TTZ_FILE"));
	file_zjj = read_settings<string>(settings_file, static_cast<string>("ZJJ_FILE"));
	file_tz = read_settings<string>(settings_file, static_cast<string>("TZ_FILE"));
	file_tbz = read_settings<string>(settings_file, static_cast<string>("TBZ_FILE"));
	file_signal = read_settings<string>(settings_file, static_cast<string>("SIGNAL_FILE"));
	
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded settings from " << settings_file << ":" << endl;
	cout << "Input file ttz: " << file_ttz << endl;
	cout << "Input file zjj: " << file_zjj << endl;
	cout << "Input file tz: " << file_tz << endl;
	cout << "Input file tbz: " << file_tbz << endl;
	cout << "Input file signal: " << file_signal << endl;
	cout << "################################################################################" << endl;
	return true;
}

// load cross sections
bool load_xsec(const string &settings_file, double &xsec_ttz, double &xsec_zjj, double &xsec_tz, double &xsec_tbz, double &xsec_signal)
{
	// read settings
	xsec_ttz = read_settings<double>(settings_file, static_cast<string>("TTZ_XSEC"));
	xsec_zjj = read_settings<double>(settings_file, static_cast<string>("ZJJ_XSEC"));
	xsec_tz = read_settings<double>(settings_file, static_cast<string>("TZ_XSEC"));
	xsec_tbz = read_settings<double>(settings_file, static_cast<string>("TBZ_XSEC"));
	xsec_signal = read_settings<double>(settings_file, static_cast<string>("SIGNAL_XSEC"));
	
	// display the loaded settings if no errors occur
	cout << "################################################################################" << endl;
	cout << "Loaded settings from " << settings_file << ":" << endl;
	cout << "Input xsec ttz [pb]: " << xsec_ttz << endl;
	cout << "Input xsec zjj [pb]: " << xsec_zjj << endl;
	cout << "Input xsec tz [pb]: " << xsec_tz << endl;
	cout << "Input xsec tbz [pb]: " << xsec_tbz << endl;
	cout << "Input xsec signal [pb]: " << xsec_signal << endl;
	cout << "################################################################################" << endl;
	return true;
}