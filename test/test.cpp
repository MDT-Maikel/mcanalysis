//  std includes
#include <iostream> // std input-output stream
#include <vector> // vector type
#include <string> // string type
#include <sstream> // stringstream type
#include <fstream> // ifstream type
#include <cmath> // math library

// local includes
#include "../source/histogram/histogram.h"

using namespace std;
using namespace analysis;

int main(int argc, const char* argv[]) 
{

	histogram test;

	//=== Set options ===//
	test.set_ps_title("hello world.ps");
	test.set_hist_bins(10);
	test.set_hist_range(0,10);
	test.set_x_label("x label");
	test.set_y_label("y label");
	test.set_leg_title("Legend");

	vector<double> a = {1,2,3,4,5,1,2,3,4,5,5};
	vector<double> b = {1,2,3,3,3,2,2,3,4,5,5};			

	//=== Add samples ===//
	test.add_sample(a);
	test.add_sample(b);
	//test.add_sample("sample3");
	//test.add_sample("sample4");

	//=== Draw histograms ===//
	test.draw();

	return 0;
}
