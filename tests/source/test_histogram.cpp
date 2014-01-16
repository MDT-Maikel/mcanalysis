/* Histogram Tests
 *
*/

#include <iostream> 
#include <random>
#include <vector> 

#include "../../source/histogram/histogram.h"
#include "../../source/histogram/histogram2D.h"

using namespace std;
using namespace analysis;


// main program
int main(int argc, const char* argv[]) 
{
	// test one dimensional histogram: counting x 
	histogram hist_1d;
	hist_1d.set_title("output/test_histogram_1D");
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


	// test two dimensional histogram: counting x vs y
	histogram2D hist_2d_xy;
	hist_2d_xy.set_title("output/test_histogram_2D_XY");
	hist_2d_xy.set_x_bins(100);
	hist_2d_xy.set_y_bins(100);
	hist_2d_xy.set_x_range(0, 100);
	hist_2d_xy.set_y_range(0, 100);
	hist_2d_xy.set_x_label("x");
	hist_2d_xy.set_y_label("y");
	
	// create the xy sample using an exp for x and a gaussian for y
	vector< vector<double> > sample_xy;
	for (int i = 0; i < 500000; i++)
		sample_xy.push_back({dist_exp(rd), dist_gauss(rd)});
	
	hist_2d_xy.add_sample_xy(sample_xy);
	hist_2d_xy.draw();
		
	
	// test two-dimensional histogram: plotting (x, y, z) pairs
	histogram2D hist_2d_xyz;
	hist_2d_xyz.set_x_bins(10);
	hist_2d_xyz.set_y_bins(10);
	hist_2d_xyz.set_x_range(1, 6);
	hist_2d_xyz.set_y_range(1, 10);
	hist_2d_xyz.set_title("output/test_histogram_2D_XYZ");
	hist_2d_xyz.set_hist_title("histogram");
	hist_2d_xyz.set_x_label("x label");
	hist_2d_xyz.set_y_label("y label");
	hist_2d_xyz.set_palette();

	// sample to histogram and draw
	vector< vector<double> > sample_xyz = {
		{1, 1, 9}, {1, 2, 3}, {1, 3, 4}, {1, 4, 5}, {1, 5, 6}, {1, 6, 7}, {1, 7, 8}, {1, 9, 10},
		{2, 1, 2}, {2, 2, 3}, {2, 3, 4}, {2, 4, 5}, {2, 5, 6}, {2, 6, 7}, {2, 7, 8}, {2, 9, 10}, 
		{3, 1, 2}, {3, 2, 3}, {3, 3, 4}, {3, 4, 5}, {3, 5, 6}, {3, 6, 7}, {3, 7, 8}, {3, 9, 10}, 
		{4, 1, 2}, {4, 2, 3}, {4, 3, 4}, {4, 4, 5}, {4, 5, 6}, {4, 6, 7}, {4, 7, 8}, {4, 9, 10}, 
		{5, 1, 2}, {5, 2, 3}, {5, 3, 4}, {5, 4, 5}, {5, 5, 6}, {5, 6, 7}, {5, 7, 8}, {5, 9, 10}		
	};
	hist_2d_xyz.add_sample_xyz(sample_xyz);
	hist_2d_xyz.draw();
}
