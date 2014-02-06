/* 2D Histogram class
 *
 * 
*/

#ifndef INC_HISTOGRAM2D
#define INC_HISTOGRAM2D

#include <iostream> 
#include <vector> 
#include <string>
#include <fstream>
#include <cmath>
#include <map>
#include <list>

#include <boost/python.hpp>

#include <TCanvas.h> 
#include <TStyle.h>
#include <TLegend.h>
#include <TH2F.h>
#include <TColor.h>


 namespace analysis 
 {

	class histogram2D
	{

	public:
	
		/* con & destructor */
		histogram2D();
		~histogram2D();

		/* histogram options */
		void set_x_bins(int nbx);
		void set_y_bins(int nby);
		void set_x_range(double x_min, double x_max);
		void set_y_range(double y_min, double y_max);
		void set_title(std::string ps);
		void set_x_label(std::string x);
		void set_y_label(std::string y);
		void set_palette(std::string name = "", const Int_t ncont = 999);
		
		/* histogram data */
		void add_sample_xyz(const std::vector< std::vector<double> > & list_xyz); 
		void add_sample_xy(const std::vector< std::vector<double> > & list_xy);

		/* interpolation */
		double Interpolate2D(const std::list< std::list<double> > & values, const std::list<double> & xy);

		/* histogram drawing */
		void draw();
		
	private:
	
		/* histogram data */
		std::vector< std::vector<double> > sample_xyz;
		std::vector< std::vector< std::vector<double> > > sample_xy;
		
		/* histogram options */
		int xbins;
		int ybins;
		double xmin; 
		double xmax;
		double ymin; 
		double ymax;
		std::string ps_title;
		std::string x_label;
		std::string y_label;
		
	};


/* NAMESPACE */
}

#endif
