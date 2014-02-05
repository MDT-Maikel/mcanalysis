/* Histogram Tests
 *
 * Test all the different histogram classes.
 * 
*/

#include <iostream> 
#include <cmath>
#include <ctime>
#include <random>
#include <vector> 

#include <boost/filesystem.hpp>

#include "histogram/histogram.h"
#include "histogram/histogram2D.h"

using namespace std;
using namespace boost::filesystem;
using namespace analysis;


// main program
int main(int argc, const char* argv[]) 
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;
	
	// remove possible existing output files
	remove("output/test_histogram_1d.png");
	remove("output/test_histogram_2d_xy.png");
	remove("output/test_histogram_2d_xyz.png");
	
	
	/* 1D histogram */
		
	// test one dimensional histogram: counting x 
	histogram hist_1d;
	hist_1d.set_title("output/test_histogram_1d");
	hist_1d.set_bins(100);
	hist_1d.set_range(0, 100);
	hist_1d.set_x_label("x");
	hist_1d.set_y_label("y");
	hist_1d.set_stacked(true);
	
	// create three three samples
	vector<double> exp, gauss, gamma;
	random_device rd;
	exponential_distribution<double> dist_exp(0.05);
	normal_distribution<double> dist_gauss(20.0, 5.0);
	gamma_distribution<double> dist_gamma(3.0, 2.0);
	for (int i = 0; i < 50000; i++)
	{
		exp.push_back(dist_exp(rd));
		gauss.push_back(dist_gauss(rd));
		gamma.push_back(50.0 + dist_gamma(rd));		
	}

	// add samples to histogram and draw
	hist_1d.add_sample(exp, "exp", 10);
	hist_1d.add_sample(gauss, "gauss");
	hist_1d.add_sample(gamma, "gamma");
	hist_1d.draw();
	bool hist_1d_success = is_regular_file("output/test_histogram_1d.png");


	/* 2D histogram: XY pairs */

	// test two dimensional histogram: counting x vs y
	unsigned int nr_bins = 10;
	histogram2D hist_2d_xy;
	hist_2d_xy.set_title("output/test_histogram_2d_xy");
	hist_2d_xy.set_x_bins(nr_bins);
	hist_2d_xy.set_y_bins(nr_bins);
	hist_2d_xy.set_x_range(0, 100);
	hist_2d_xy.set_y_range(0, 100);
	hist_2d_xy.set_x_label("x");
	hist_2d_xy.set_y_label("y");
	hist_2d_xy.set_palette();
	
	// create the xy sample using a gaussian for both x and y
	normal_distribution<double> dist_gauss_xy(50.0, 25.0);
	vector< vector<double> > sample_xy;
	for (unsigned int i = 0; i < 500000; i++)
		sample_xy.push_back({dist_gauss_xy(rd), dist_gauss_xy(rd)});
		
	// draw the histogram	
	hist_2d_xy.add_sample_xy(sample_xy);
	hist_2d_xy.draw();
	bool hist_2d_xy_success = is_regular_file("output/test_histogram_2d_xy.png");
		
		
	/* 2D histogram: XYZ pairs */
	
	// test two-dimensional histogram: plotting (x, y, z) pairs
	histogram2D hist_2d_xyz;
	hist_2d_xyz.set_x_bins(2 * nr_bins);
	hist_2d_xyz.set_y_bins(2 * nr_bins);
	hist_2d_xyz.set_x_range(0, 100);
	hist_2d_xyz.set_y_range(0, 100);
	hist_2d_xyz.set_title("output/test_histogram_2d_xyz");
	hist_2d_xyz.set_x_label("x");
	hist_2d_xyz.set_y_label("y");
	hist_2d_xyz.set_palette();
	
	// create the xyz sample by using the same sample as for the xy case
	vector< vector<double> > sample_xyz;
	for (unsigned int i = 0; i < nr_bins; i++)
	{
		vector<double> zero_vec(nr_bins, 0.0);
		sample_xyz.push_back(zero_vec);
	}
	for (unsigned int i = 0; i < sample_xy.size(); i++)
	{
		double x = sample_xy[i][0];
		double y = sample_xy[i][1];
		
		unsigned int posx = static_cast<unsigned int>(round(nr_bins * x / 100));
		unsigned int posy = static_cast<unsigned int>(round(nr_bins * y / 100));
		
		if (posx >= 0 && posx < nr_bins && posy >= 0 && posy < nr_bins)
			sample_xyz[posx][posy]++;
	}
	
	// convert binned data to list
	vector< vector<double> > list_xyz;
	for (unsigned int i = 0; i < sample_xyz.size(); i++)
	{
		vector<double> xyz = sample_xyz[i];
		for (unsigned int j = 0; j < xyz.size(); j++)
		{
			double x = static_cast<double>(i) * 100.0 / nr_bins + 50.0 / nr_bins;
			double y = static_cast<double>(j) * 100.0 / nr_bins + 50.0 / nr_bins;
			list_xyz.push_back({x, y, sample_xyz[i][j]});
		}
	}
	
	// draw the histogram
	hist_2d_xyz.add_sample_xyz(list_xyz);
	hist_2d_xyz.draw();
	bool hist_2d_xyz_success = is_regular_file("output/test_histogram_2d_xyz.png");
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Histogram test: completed in " << duration << " seconds." << endl;
	cout << "1D Histogram creation has " << (hist_1d_success ? "succeeded!" : "failed!") << endl; 
	cout << "Check: exponential falling background with gaussian and gamma peaks." << endl;
	cout << "2D Histogram XY creation has " << (hist_2d_xy_success ? "succeeded!" : "failed!") << endl;
	cout << "Check: 2D gaussian peak around (50, 50)." << endl;
	cout << "2D Histogram XYZ creation has " << (hist_2d_xyz_success ? "succeeded!" : "failed!") << endl;
	cout << "Check: same 2D gaussian peak as XY with double resolution." << endl;
	cout << "=====================================================================" << endl;
		
	// return whether tests passed
	if (hist_1d_success && hist_2d_xy_success && hist_2d_xyz_success)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}
