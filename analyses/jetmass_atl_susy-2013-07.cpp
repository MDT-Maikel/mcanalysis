/* ATLAS jet mass analysis
 *
 * Performs the CMS dijet analysis based on ATL_SUSY-2013-07,
 * these limits can be compared to the model-independent limits
 * from table VII in their conference note. This code calculates
 * the efficiency for the given events file.
 * 
*/

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "cuts/cuts.h"
#include "event/event.h"
#include "particle/particle.h"
#include "utility/utility.h"



using namespace std;
using namespace analysis;
using namespace boost;
using namespace boost::filesystem;


// atl 4 jet cut
class cut_atl_4jets: public cut
{
public:
	cut_atl_4jets() {}

	bool operator() (const event *ev) 
	{ 
		// get the fourth jet
		const particle *j4 = ev->get(ptype_jet, 4, 2.5);
		
		// reject event if jets not found
		if (!j4)
			return false;
		
		// reject event if pt <= 100 GeV
		if (j4->pt() <= 100)
			return false;

		// event passed
		return true;
	}
};

// atl delta eta cut
class cut_atl_delta_eta: public cut
{
public:
	cut_atl_delta_eta() {}

	bool operator() (const event *ev) 
	{ 
		const particle *j1 = ev->get(ptype_jet, 1, 2.5);
		const particle *j2 = ev->get(ptype_jet, 2, 2.5);

		return delta_eta(j1, j2) < 0.7;
	}
};

// atl pt jet3 cut
class cut_atl_pt_jet3: public cut
{
public:
	cut_atl_pt_jet3(double pt) : pt_cut(pt) {}

	bool operator() (const event *ev) 
	{ 
		const particle *j3 = ev->get(ptype_jet, 3, 2.5);

		return j3->pt() > pt_cut;
	}
private:
	double pt_cut;
};

// atl jet mass cut
class cut_atl_jet_mass: public cut
{
public:
	cut_atl_jet_mass(double mmin, double mmax) : jet_mass_min(mmin), jet_mass_max(mmax) {}

	bool operator() (const event *ev) 
	{ 
		const particle *j1 = ev->get(ptype_jet, 1, 2.5);
		const particle *j2 = ev->get(ptype_jet, 2, 2.5);
		const particle *j3 = ev->get(ptype_jet, 3, 2.5);
		const particle *j4 = ev->get(ptype_jet, 4, 2.5);

		double jet_mass = j1->mass() + j2->mass() + j3->mass() + j4->mass();
	
		return jet_mass > jet_mass_min && jet_mass < jet_mass_max;
	}
private:
	double jet_mass_min;
	double jet_mass_max;
};

// main program with one argument representing the input file 
int main(int argc, char* argv[])
{
	// make sure there is one arguments
	if (argc != 2) 
	{
		cout << "specify at least one argument: <input.mc.gz>" << endl;
		return EXIT_FAILURE;
  	}

	// convert the arguments to the input file and output directory
	string input_file = argv[1];
	
	// load the events dependent on whether they are .lhe.gz or .lhco.gz
	vector<event*> events; 	
	read_lhco(events, input_file);
	
	// initialize cutlist and add all cms dijet cuts

	cut_atl_4jets *fourjet_cut = new cut_atl_4jets();
	cut_atl_delta_eta *eta_cut = new cut_atl_delta_eta();
	cut_atl_pt_jet3 *jetpt_cut_100 = new cut_atl_pt_jet3(100.0);	
	cut_atl_pt_jet3 *jetpt_cut_250 = new cut_atl_pt_jet3(250.0);
	cut_atl_jet_mass *jetmass_cut_full = new cut_atl_jet_mass(625.0, 10000.0);
	cut_atl_jet_mass *jetmass_cut_1 = new cut_atl_jet_mass(350.0, 400.0);
	cut_atl_jet_mass *jetmass_cut_2 = new cut_atl_jet_mass(400.0, 450.0);
	cut_atl_jet_mass *jetmass_cut_3 = new cut_atl_jet_mass(450.0, 525.0);
	cut_atl_jet_mass *jetmass_cut_4 = new cut_atl_jet_mass(525.0, 725.0);
	cut_atl_jet_mass *jetmass_cut_5 = new cut_atl_jet_mass(725.0, 10000.0);

	cuts sr1;
	sr1.add_cut(fourjet_cut, "atl 4jet cut");
	sr1.add_cut(eta_cut, "atl 4jet cut");
	sr1.add_cut(jetpt_cut_250, "atl 4jet cut");
	sr1.add_cut(jetmass_cut_full, "atl jet mass cut");
	vector<event*> events_sr1 = copy_events(events);
	sr1.apply(events_sr1);
	double acc_sr1 = sr1.efficiency();
	delete_events(events_sr1);

	cuts sr100_1;
	sr100_1.add_cut(fourjet_cut, "atl 4jet cut");
	sr100_1.add_cut(eta_cut, "atl 4jet cut");
	sr100_1.add_cut(jetpt_cut_100, "atl 4jet cut");
	sr100_1.add_cut(jetmass_cut_1, "atl jet mass cut");
	vector<event*> events_sr100_1 = copy_events(events);
	sr100_1.apply(events_sr100_1);
	double acc_sr100_1 = sr100_1.efficiency();
	delete_events(events_sr100_1);

	cuts sr100_2;
	sr100_2.add_cut(fourjet_cut, "atl 4jet cut");
	sr100_2.add_cut(eta_cut, "atl 4jet cut");
	sr100_2.add_cut(jetpt_cut_100, "atl 4jet cut");
	sr100_2.add_cut(jetmass_cut_2, "atl jet mass cut");
	vector<event*> events_sr100_2 = copy_events(events);
	sr100_2.apply(events_sr100_2);
	double acc_sr100_2 = sr100_2.efficiency();
	delete_events(events_sr100_2);

	cuts sr100_3;
	sr100_3.add_cut(fourjet_cut, "atl 4jet cut");
	sr100_3.add_cut(eta_cut, "atl 4jet cut");
	sr100_3.add_cut(jetpt_cut_100, "atl 4jet cut");
	sr100_3.add_cut(jetmass_cut_3, "atl jet mass cut");
	vector<event*> events_sr100_3 = copy_events(events);
	sr100_3.apply(events_sr100_3);
	double acc_sr100_3 = sr100_3.efficiency();
	delete_events(events_sr100_3);

	cuts sr100_4;
	sr100_4.add_cut(fourjet_cut, "atl 4jet cut");
	sr100_4.add_cut(eta_cut, "atl 4jet cut");
	sr100_4.add_cut(jetpt_cut_100, "atl 4jet cut");
	sr100_4.add_cut(jetmass_cut_4, "atl jet mass cut");
	vector<event*> events_sr100_4 = copy_events(events);
	sr100_4.apply(events_sr100_4);
	double acc_sr100_4 = sr100_4.efficiency();
	delete_events(events_sr100_4);

	cuts sr100_5;
	sr100_5.add_cut(fourjet_cut, "atl 4jet cut");
	sr100_5.add_cut(eta_cut, "atl 4jet cut");
	sr100_5.add_cut(jetpt_cut_100, "atl 4jet cut");
	sr100_5.add_cut(jetmass_cut_5, "atl jet mass cut");
	vector<event*> events_sr100_5 = copy_events(events);
	sr100_5.apply(events_sr100_5);
	double acc_sr100_5 = sr100_5.efficiency();
	delete_events(events_sr100_5);

	cuts sr250_1;
	sr250_1.add_cut(fourjet_cut, "atl 4jet cut");
	sr250_1.add_cut(eta_cut, "atl 4jet cut");
	sr250_1.add_cut(jetpt_cut_250, "atl 4jet cut");
	sr250_1.add_cut(jetmass_cut_1, "atl jet mass cut");
	vector<event*> events_sr250_1 = copy_events(events);
	sr250_1.apply(events_sr250_1);
	double acc_sr250_1 = sr250_1.efficiency();
	delete_events(events_sr250_1);

	cuts sr250_2;
	sr250_2.add_cut(fourjet_cut, "atl 4jet cut");
	sr250_2.add_cut(eta_cut, "atl 4jet cut");
	sr250_2.add_cut(jetpt_cut_250, "atl 4jet cut");
	sr250_2.add_cut(jetmass_cut_2, "atl jet mass cut");
	vector<event*> events_sr250_2 = copy_events(events);
	sr250_2.apply(events_sr250_2);
	double acc_sr250_2 = sr250_2.efficiency();
	delete_events(events_sr250_2);

	cuts sr250_3;
	sr250_3.add_cut(fourjet_cut, "atl 4jet cut");
	sr250_3.add_cut(eta_cut, "atl 4jet cut");
	sr250_3.add_cut(jetpt_cut_250, "atl 4jet cut");
	sr250_3.add_cut(jetmass_cut_3, "atl jet mass cut");
	vector<event*> events_sr250_3 = copy_events(events);
	sr250_3.apply(events_sr250_3);
	double acc_sr250_3 = sr250_3.efficiency();
	delete_events(events_sr250_3);

	cuts sr250_4;
	sr250_4.add_cut(fourjet_cut, "atl 4jet cut");
	sr250_4.add_cut(eta_cut, "atl 4jet cut");
	sr250_4.add_cut(jetpt_cut_250, "atl 4jet cut");
	sr250_4.add_cut(jetmass_cut_4, "atl jet mass cut");
	vector<event*> events_sr250_4 = copy_events(events);
	sr250_4.apply(events_sr250_4);
	double acc_sr250_4 = sr250_4.efficiency();
	delete_events(events_sr250_4);

	cuts sr250_5;
	sr250_5.add_cut(fourjet_cut, "atl 4jet cut");
	sr250_5.add_cut(eta_cut, "atl 4jet cut");
	sr250_5.add_cut(jetpt_cut_250, "atl 4jet cut");
	sr250_5.add_cut(jetmass_cut_5, "atl jet mass cut");
	vector<event*> events_sr250_5 = copy_events(events);
	sr250_5.apply(events_sr250_5);
	double acc_sr250_5 = sr250_5.efficiency();
	delete_events(events_sr250_5);

	// delete all the cut pointers
	delete fourjet_cut, eta_cut, jetpt_cut_100, jetpt_cut_250, jetmass_cut_full, jetmass_cut_1, jetmass_cut_2, jetmass_cut_3, jetmass_cut_4, jetmass_cut_5;
		
	// clear remaining event pointers
	delete_events(events);
	
	// just output the acceptance
	cout << acc_sr1 << " " << acc_sr100_1 << " " << acc_sr100_2 << " " << acc_sr100_3 << " " << acc_sr100_4 << " " << acc_sr100_5 << " " << acc_sr250_1 << " " << acc_sr250_2 << " " << acc_sr250_3 << " " << acc_sr250_4 << " " << acc_sr250_5 << endl;
	
	// finished the analysis
	return EXIT_SUCCESS;	
}

