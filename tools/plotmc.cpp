/* Plot MC
 *
 * Plots a whole bunch of standard plots, specified by the user for the given input file.
 * 
*/

#include <iostream>

#include <boost/filesystem.hpp>

#include "cuts/cuts.h"
#include "event/event.h"
#include "plot/plot.h"
#include "utility/utility.h"


using namespace std;
using namespace boost;
using namespace boost::filesystem;
using namespace analysis;


// necessary function prototypes
void read_events(vector<event*> & events, const string & input_file);
void plot_particle_pt(const vector<event*> & events, const string & output_dir, unsigned int type, unsigned int number);
void plot_particle_eta(const vector<event*> & events, const string & output_dir, unsigned int type, unsigned int number);
void plot_particle_phi(const vector<event*> & events, const string & output_dir, unsigned int type, unsigned int number);
void plot_event_met(const vector<event*> & events, const string & output_dir);


// main program with two arguments representing the input file 
// and the output folder for the plots
int main(int argc, char* argv[])
{
	// make sure there are two arguments
	if (argc != 3) 
	{
		cout << "specify two arguments: <input.mc.gz> <output_dir/>" << endl;
		return EXIT_FAILURE;
  	}

	// convert the arguments to the input file and output directory
	string input_file = argv[1];
	string output_dir = argv[2];
		
	// create output directory if it does not exist yet
	if (!is_directory(output_dir))
		create_directory(output_dir);
	
	// load the events dependent on whether they are .lhe.gz or .lhco.gz
	vector<event*> events; 	
	read_events(events, input_file);
	
	// plot pt, eta and phi of the n leading jets
	plot_particle_pt(events, output_dir, ptype_jet, 2);
	plot_particle_eta(events, output_dir, ptype_jet, 3);
	plot_particle_phi(events, output_dir, ptype_electron | ptype_muon | ptype_tau, 1);
	
	// plot general event variables
	plot_event_met(events, output_dir);
	
	// clear remaining event pointers
	delete_events(events);
	
	// finished the plotting
	return EXIT_SUCCESS;
}

// read the events dependent on the file type
void read_events(vector<event*> & events, const string & input_file)
{
	// determine whether the input file is *.lhe.gz then load events
	std::regex lhe_match("(.*)(lhe.gz)");
	if (std::regex_match(input_file, lhe_match))
	{
		read_lhe(events, input_file);
		return;
	}
	
	// determine whether the input file is *.lhco.gz then load events
	std::regex lhco_match("(.*)(lhco.gz)");
	if (std::regex_match(input_file, lhco_match))
	{
		read_lhco(events, input_file);
		return;
	}
	
	// if reached here, the process failed and events remain unchanged
	return;
}


void plot_particle_pt(const vector<event*> & events, const string & output_dir, unsigned int type, unsigned int number)
{
	// cut away events that do not have that specific particle
	cuts cp;	
	cut_particle *ft_p = new cut_particle(type, number);
	cp.add_cut(ft_p);
	const vector<event*> red_events(cp.reduce(events));

	// log details of what is being plotted
	cout << "plotting pt of particle: " << ptype_to_string(type) << " and number: " << number << " for " << red_events.size() << " events" << endl;
	
	// plot the pt of the nth particle
	plot_pt *ft_pt = new plot_pt(type, number);
	string plot_name = "pt_" + ptype_to_string(type) + "_" + boost::lexical_cast<std::string>(number);
	plot pt(plot_name, output_dir);
	pt.add_sample(red_events, ft_pt, "sample");
	pt.run();
}

void plot_particle_eta(const vector<event*> & events, const string & output_dir, unsigned int type, unsigned int number)
{
	// cut away events that do not have that specific particle
	cuts cp;	
	cut_particle *ft_p = new cut_particle(type, number);
	cp.add_cut(ft_p);
	const vector<event*> red_events(cp.reduce(events));

	// log details of what is being plotted
	cout << "plotting pt of particle: " << ptype_to_string(type) << " and number: " << number << " for " << red_events.size() << " events" << endl;
	
	// plot the pt of the nth particle
	plot_eta *ft_eta = new plot_eta(type, number);
	string plot_name = "eta_" + ptype_to_string(type) + "_" + boost::lexical_cast<std::string>(number);
	plot eta(plot_name, output_dir);
	eta.add_sample(red_events, ft_eta, "sample");
	eta.run();
}

void plot_particle_phi(const vector<event*> & events, const string & output_dir, unsigned int type, unsigned int number)
{
	// cut away events that do not have that specific particle
	cuts cp;	
	cut_particle *ft_p = new cut_particle(type, number);
	cp.add_cut(ft_p);
	const vector<event*> red_events(cp.reduce(events));

	// log details of what is being plotted
	cout << "plotting pt of particle: " << ptype_to_string(type) << " and number: " << number << " for " << red_events.size() << " events" << endl;
	
	// plot the pt of the nth particle
	plot_phi *ft_phi = new plot_phi(type, number);
	string plot_name = "phi_" + ptype_to_string(type) + "_" + boost::lexical_cast<std::string>(number);
	plot phi(plot_name, output_dir);
	phi.add_sample(red_events, ft_phi, "sample");
	phi.run();
}

void plot_event_met(const vector<event*> & events, const string & output_dir)
{
	// log details of what is being plotted
	cout << "plotting met for " << events.size() << " events" << endl;
	
	// plot the met of the event
	plot_met *ft_met = new plot_met();
	plot met("met", output_dir);
	met.add_sample(events, ft_met, "sample");
	met.run();
}

