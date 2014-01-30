/* Merge MC
 *
 * 
*/

#include <iostream>

#include "boost/filesystem.hpp"

#include "../../source/utility/utility.h"
#include "../../source/event/event.h"


using namespace std;
using namespace analysis;
using namespace boost::filesystem;


int main(int argc, char* argv[])
{
	
	vector<path> files(get_files("files/", ".lhco.gz", true));
	
	for (unsigned int i = 0; i < files.size(); i++)
		cout << files[i].string() << endl;
	
	vector<event*> events;
	for (unsigned int i = 0; i < files.size(); i++)
		read_lhco(events, files[i]);

	write_lhco(events, "./merged.lhco.gz");

	// delete events
	for (unsigned int i = 0; i < events.size(); i++)
		delete events[i];
	events.clear();

	return 0;
}
