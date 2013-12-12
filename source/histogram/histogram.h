/* Histogram class
 *
 * 
*/

#ifndef INC_HISTOGRAM
#define INC_HISTOGRAM

#include <iostream> // std input-output stream
#include <vector> // vector type
#include <string> // string type
#include <sstream> // stringstream type
#include <fstream> // ifstream type
#include <cmath> // math library

#include <TCanvas.h> //basic canvas class
#include <TStyle.h> // style class
#include <TLegend.h> // draw legend class
#include <TH1.h> // 1D histogram class
// #include <TF1.h> // draw 1D functions class
// #include <TH2.h> // 2D histogram class
// #include <TGraph.h> // graphics class
// #include <TNtuple.h> // N-tuple class
// #include <TVector3.h> // three-vector class
// #include <TLorentzVector.h> // Lorentz vector class
// #include <TRandom.h> // random number generation class
// #include <TRandom3.h> // extended random number generation class
// #include <TROOT.h> // entry point to the ROOT system
// #include <TChain.h> // access to collection of files containg TTree objects
// #include <TFile.h> // a ROOT file is a suite of consecutive data records (TKey's)
// #include <TTree.h> // a TTree object consists of a list of independent branches (TBranch)
// #include <TPaveText.h> // a PaveText is a Pave (see TPave) with text, lines or/and boxes inside.
// #include <Rtypes.h> // basic types used by ROOT


namespace analysis {

	class histogram
	{
	private:
		TCanvas* c1;
		std::vector< std::vector<double> > samples;
		double nbins;
		double hmin; 
		double hmax;
		std::string ps_title;
		std::string leg_title;
		std::string x_label;
		std::string y_label;
		bool is_normalised;

	public:
		//=== Class constructor ===//
		histogram();

		//=== Class destructor ===//
		~histogram();

		//=== Set histogram Options ===//
		void set_ps_title(std::string ps);
		void set_hist_bins(double nb);
		void set_hist_range(double min, double max);
		void set_x_label(std::string x);
		void set_y_label(std::string y);
		void set_leg_title(std::string title);
		void add_sample(std::vector<double> sample);
		void normalise();

		//=== Draw histogram ===//
		void draw();
	};


/* NAMESPACE */
}

#endif
