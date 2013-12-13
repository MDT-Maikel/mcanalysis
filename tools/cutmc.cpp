/* Cut MC
 *
 * 
*/

#include <iostream>

#include "boost/filesystem.hpp"

#include "../source/utility/utility.h"
#include "../source/event/event.h"
#include "../source/cuts/cuts.h"
#include "../source/cuts/cuts_standard.h"
#include "../source/cuts/cuts_advanced.h"
#include "../source/histogram/histogram.h"


using namespace std;
using namespace analysis;
using namespace boost::filesystem;


int main(int argc, char* argv[])
{
	
	vector<path> files(get_files("../files/sized", ".lhco.gz", true));
	for (unsigned int i = 0; i < files.size(); i++)
		cout << files[i].string() << endl;
	
	vector<event*> events;
	for (unsigned int i = 0; i < files.size(); i++)
		read_lhco(events, files[i]);

	// electron veto
	cut_veto *veto; veto = new cut_veto(particle::type_electron, 40, 3.14);
	cut_veto *vetomu; vetomu = new cut_veto(particle::type_muon, 40, 3.14);	

	// pt cut
	cut_pt *pt; pt = new cut_pt;
	pt->set_type(particle::type_jet);
	pt->set_number(1);
	pt->set_pt(130);

	// MET cut
	cut_met *met; met = new cut_met;
	met->set_met(20);

	// HT cut
	cut_ht *ht; ht = new cut_ht(particle::type_jet, 10, 3.14, 210);

	cuts test;
	test.add_cut(veto);
	test.add_cut(vetomu);
	test.add_cut(pt);
	test.add_cut(met);
	test.add_cut(ht);
	test.apply(events);

	test.write(cout);

	write_lhco(events, "./cutted.lhco.gz");

	// delete cuts
	delete pt; delete met; delete ht;


	// plot met
	vector<double> metlist;
	for (unsigned int i = 0; i < events.size(); i++)
		metlist.push_back(events[i]->met());

	vector<double> htlist;
	for (unsigned int i = 0; i < events.size(); i++)
		htlist.push_back(events[i]->ht(particle::type_jet, 40, 2.5));

	histogram tast;

	tast.set_title("met.ps");
	tast.set_bins(100);
	tast.set_range(0,1000);
	tast.set_x_label("met");
	tast.set_y_label("# events");
	tast.set_leg_title("Legend");		
	tast.add_sample(metlist);
	tast.add_sample(htlist);
	tast.normalize();
	tast.draw();


	// delete events
	for (unsigned int i = 0; i < events.size(); i++)
		delete events[i];
	events.clear();

	return 0;
}
