/* BumpHunter Tests
 * 
 * Test the different features and methods of the bumphunter class.
 *
*/

#include <ctime>
#include <iostream>
#include <random>
#include <string>

#include <boost/filesystem.hpp>

#include <TF1.h>
#include <TH1.h> 

#include "../../source/bumphunter/bumphunter.h"

using namespace std;
using namespace boost::filesystem;
using namespace analysis;


// main program
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old;
	clock_old = clock();
	double duration;
	
	// remove possible existing output files
	remove("output/test_bumphunter_pred_data.png");
	remove("output/test_bumphunter_poisson.png");
	remove("output/test_bumphunter_gaussian.png");
	remove("output/test_bumphunter_poissongamma.png");
	
	// initiate root with a random seed		
	std::random_device rd;
	gRandom->SetSeed(rd());
	
	// parameters
	double gauss_height = 60;
	double gauss_width = 20;
	double gauss_mean = 500;
	double exp_const = 1500;
	double exp_decay = 0.003;
	int nr_bins = 20;
	
	// create data: falling exponent
	TF1 *f_pred = new TF1("fpred", "[0]*exp([1]*x)", 0, 1000);
	f_pred->SetParameters(exp_const, -exp_decay);
		
	// create prediction: same exponent + gaussian resonance
	TF1 *f_data = new TF1("fdata", "[0]*exp([1]*x) + gaus(2)", 0, 1000);
	f_data->SetParameters(exp_const, -exp_decay, gauss_height, gauss_mean, gauss_width);

	// fill histograms with data and prediction
	TH1F* hist_pred = new TH1F("", "pred", nr_bins, 0, 1000);
	hist_pred->FillRandom("fpred", 100000);
	TH1F* hist_data = new TH1F("", "data", nr_bins, 0, 1000);
	hist_data->FillRandom("fdata", 100000);
	
	// print the prediction and data histograms
	TCanvas* canvas = new TCanvas("test_bumphunter", "");
	hist_pred->Draw("same");
	hist_data->Draw("same");
	canvas->Print("output/test_bumphunter_pred_data.png");
	
	// create a bumphunter instance and do basic settings
	bumphunter hunt(hist_pred, hist_data);
	hunt.set_folder("output/");
	hunt.set_name("test_bumphunter_poisson");
	hunt.SetNPseudoExperiments(5000);
	hunt.SetBinModel(bumphunter::BUMP_POISSON);
	hunt.SetSearchRegion(100, 1000);
	hunt.SetMinWindowSize(2);
	hunt.SetMaxWindowSize(4);
	hunt.SetWindowStepSize(1);
	
	// run bumphunter analysis
	hunt.run();
	double sigma_poisson = hunt.get_global_sigma();
	
	// run bumphunter analysis a second time with gaussian bin model
	hunt.SetBinModel(bumphunter::BUMP_GAUSSIAN);
	hunt.set_name("test_bumphunter_gaussian");
	hunt.run();
	double sigma_gaussian = hunt.get_global_sigma();
	
	// run bumphunter analysis a second time with poisson-gamma bin model
	hunt.SetBinModel(bumphunter::BUMP_POISSON_GAMMA);
	hunt.set_name("test_bumphunter_poissongamma");
	hunt.run();
	double sigma_poissongamma = hunt.get_global_sigma();
	
	// determine success conditions
	bool bumphunter_success = is_regular_file("output/test_bumphunter_poisson.png");
	bumphunter_success = bumphunter_success && is_regular_file("output/test_bumphunter_gaussian.png");
	bumphunter_success = bumphunter_success && is_regular_file("output/test_bumphunter_poissongamma.png");
	bumphunter_success = bumphunter_success && sigma_poisson >= sigma_gaussian;
	bumphunter_success = bumphunter_success && sigma_poisson >= sigma_poissongamma;
	
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "BumpHunter test: hunt completed in " << duration << " seconds." << endl;
	cout << "BumpHunter algorithm has " << (bumphunter_success ? "succeeded!" : "failed!") << endl;
	cout << "Poisson bin model significance: " << sigma_poisson << " sigma" << endl;
	cout << "Gaussian bin model significance: " << sigma_gaussian << " sigma" << endl;
	cout << "PoissonGamma bin model significance: " << sigma_poissongamma << " sigma" << endl;
	cout << "CHECK: what exactly?" << endl;
	cout << "=====================================================================" << endl;
		
	// clear remaining pointers
	delete canvas;
	delete f_pred;
	delete f_data;
	delete hist_pred;
	delete hist_data;
	
	// return whether tests passed
	if (bumphunter_success)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}
