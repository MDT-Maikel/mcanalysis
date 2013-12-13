/* Plot MC
 *
 * 
*/

#include <iostream>

#include "boost/filesystem.hpp"

#include "../../source/utility/utility.h"
#include "../../source/event/event.h"
#include "../../source/cuts/cuts.h"
#include "../../source/cuts/cuts_standard.h"
#include "../../source/cuts/cuts_advanced.h"
#include "../../source/histogram/histogram.h"
#include "../../source/plot/plot.h"
#include "../../source/plot/plot_standard.h"


using namespace std;
using namespace analysis;
using namespace boost::filesystem;


int main(int argc, char* argv[])
{
	
	//vector<path> files(get_files("../files/big", ".lhco.gz", true));
	//for (unsigned int i = 0; i < files.size(); i++)
	//	cout << files[i].string() << endl;
	
	// load events
	vector<event*> qcd_2j;
	read_lhco(qcd_2j, "../files/various/qcd_2j.lhco.gz");
	vector<event*> qcd_3j;
	read_lhco(qcd_3j, "../files/various/qcd_3j.lhco.gz");
	vector<event*> qcd_4j;
	read_lhco(qcd_4j, "../files/various/qcd_4j.lhco.gz");

	// plotting
	plot_pt pt1(particle::type_jet, 1);
	pt1.add_sample(qcd_2j, "qcd_2j");
	pt1.add_sample(qcd_3j, "qcd_3j");
	pt1.add_sample(qcd_4j, "qcd_4j");
	pt1.run();

	plot_pt pt2(particle::type_jet, 2);
	pt2.add_sample(qcd_2j, "qcd_2j");
	pt2.add_sample(qcd_3j, "qcd_3j");
	pt2.add_sample(qcd_4j, "qcd_4j");
	pt2.run();

	plot_met met;
	met.add_sample(qcd_2j, "qcd_2j");
	met.add_sample(qcd_3j, "qcd_3j");
	met.add_sample(qcd_4j, "qcd_4j");
	met.run();

	plot_ht ht;
	ht.add_sample(qcd_2j, "qcd_2j");
	ht.add_sample(qcd_3j, "qcd_3j");
	ht.add_sample(qcd_4j, "qcd_4j");
	ht.run();

	// delete events
	for (unsigned int i = 0; i < qcd_2j.size(); i++)
		delete qcd_2j[i];
	qcd_2j.clear();
	for (unsigned int i = 0; i < qcd_3j.size(); i++)
		delete qcd_3j[i];
	qcd_3j.clear();
	for (unsigned int i = 0; i < qcd_4j.size(); i++)
		delete qcd_4j[i];
	qcd_4j.clear();

	return 0;
}
