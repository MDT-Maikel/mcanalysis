/* Histogram class
 *
 * 
*/

#ifndef INC_HISTOGRAM
#define INC_HISTOGRAM

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string> 
#include <vector> 

#include <TCanvas.h>
#include <TColor.h>
#include <THStack.h>
#include <TH1.h> 
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h> 


#include "boost/lexical_cast.hpp"


namespace analysis {

	class histogram
	{

	public:

		/* con & destructor */
		histogram();
		~histogram();

		/* histogram options */
		void set_title(std::string ps);
		void set_bins(double nb);
		void set_range(double min, double max);
		void set_x_label(std::string x);
		void set_y_label(std::string y);
		void set_leg_title(std::string title);
		void set_logy(bool on);
		void set_normalized(bool on);
		void set_stacked(bool on);

		/* histogram data */
		void add_sample(const std::vector<double> & sample, const std::string & name = "", const std::string & line_color = "standard");

		/* histogram drawing */
		void draw();

	private:

		std::vector< std::vector<double> > sample_list;
		std::vector<std::string> sample_colors;
		std::vector<std::string> sample_names;

		/* histogram options */
		double nbins;
		double hmin; 
		double hmax;
		bool auto_range;
		std::string ps_title;
		std::string leg_title;
		std::string x_label;
		std::string y_label;
		bool is_normalised;
		bool has_logy;
		bool is_stacked;

	};

/* NAMESPACE */
}

#endif
