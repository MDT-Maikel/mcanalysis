/* 2D Histogram class
 *
 * 
*/

#ifndef INC_HISTOGRAM2D
#define INC_HISTOGRAM2D

#include <iostream> // std input-output stream
#include <vector> // vector type
#include <string> // string type
#include <sstream> // stringstream type
#include <fstream> // ifstream type
#include <cmath> // math library
#include <map> // map type
#include <list> // list type
#include <boost/python.hpp> // python converter

#include <TCanvas.h> //basic canvas class
#include <TStyle.h> // style class
#include <TLegend.h> // draw legend class
#include <TH2F.h> // 2D histogram class
#include <TColor.h> // modify color class


 namespace analysis {

	class histogram2D
	{
	private:
		int xbins;
		int ybins;
		double xmin; 
		double xmax;
		double ymin; 
		double ymax;
		std::string ps_title;
		std::string hist_title;
		std::string x_label;
		std::string y_label;
		std::string infile;

	public:
		//=== Class constructor ===//
		histogram2D();

		//=== Class destructor ===//
		// ~histogram2D();

		//=== Set histogram Options ===//
		void set_x_bins(int nbx);
		void set_y_bins(int nby);
		void set_x_range(double x_min, double x_max);
		void set_y_range(double y_min, double y_max);
		void set_ps_title(std::string ps);
		void set_hist_title(std::string title);
		void set_x_label(std::string x);
		void set_y_label(std::string y);
		void add_sample(std::string sample);
		void set_palette(std::string name = "", const Int_t NCont = 999);

		//=== Interpolating functions ===//
		double Interpolate2D(std::list< std::list<double> > values, std::list<double> xy);

		//=== Draw histogram ===//
		void draw();
	};


/* NAMESPACE */
}

#endif
