/* 2D Histogram class
 *
 * 
*/

#include "histogram2D.h"


namespace analysis 
{

	/* con & destructor */
	
	histogram2D::histogram2D() 
	{
		
		xbins = 100;
		ybins = 100;
		xmin = 0;
		xmax = 100;
		ymin = 0;
		ymax = 100;

		ps_title = "new 2D histogram";
    	x_label = ""; 
    	y_label = "";
	}

	histogram2D::~histogram2D()
	{
	}	

	/* histogram options */
	
	void histogram2D::set_x_bins(int nbx)
	{
		xbins = nbx;
	}

	void histogram2D::set_y_bins(int nby)
	{
		ybins = nby;
	}

	void histogram2D::set_x_range(double x_min, double x_max)
	{
		xmin = x_min;
		xmax = x_max;
	}

	void histogram2D::set_y_range(double y_min, double y_max)
	{
		ymin = y_min;
		ymax = y_max;
	}

	void histogram2D::set_title(std::string ps)
	{
		ps_title = ps;
	}

	void histogram2D::set_x_label(std::string x)
	{
		x_label = x;
	}

	void histogram2D::set_y_label(std::string y)
	{
		y_label = y;
	}

	void histogram2D::set_palette(std::string name, const Int_t ncont)
	{
    	const Int_t NRGBs = 5;

    	// default palette, looks cool
		Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
		Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
		Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
		Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };

    	if (name == "gray" || name == "grayscale")
    	{
    		stops[0] = 0.00; red[0] = 1.00; green[0] = 1.00; blue[0] = 1.00;
    		stops[1] = 0.34; red[1] = 0.84; green[1] = 0.84; blue[1] = 0.84;
    		stops[2] = 0.61; red[2] = 0.61; green[2] = 0.61; blue[2] = 0.61;
    		stops[3] = 0.84; red[3] = 0.34; green[3] = 0.34; blue[3] = 0.34;
    		stops[4] = 1.00; red[4] = 1.00; green[4] = 0.00; blue[4] = 0.00;
    	}
    	
    	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, ncont);
    	gStyle->SetNumberContours(ncont);
	}
	
	/* histogram data */
	void histogram2D::add_sample_xyz(const std::vector< std::vector<double> > & list_xyz)
	{
		sample_xyz = list_xyz;		
	}
	
	void histogram2D::add_sample_xy(const std::vector< std::vector<double> > & list_xy)
	{
		sample_xy.push_back(list_xy);
	}

	/* interpolation */
	
	double histogram2D::Interpolate2D(const std::list< std::list<double> > & values, const std::list<double> & xy)
	{
		try
    	{
			Py_Initialize();

	    	PyRun_SimpleString("import sys, string, random");
	    	PyRun_SimpleString("sys.path.append('../../source/histogram')");
	    	
	    	boost::python::list l_values;
	  		typename std::list< std::list<double> >::const_iterator it_outer;
	  		typename std::list<double>::const_iterator it_inner;

	  		for (it_outer = values.begin(); it_outer != values.end(); ++it_outer)
	  		{
	  			boost::python::list l_inner;
	  			for (it_inner = (*it_outer).begin(); it_inner != (*it_outer).end(); ++it_inner)
	  			{
	  				l_inner.append(*it_inner);
	  			}
	  			l_values.append(l_inner);
	  		}  

	  		boost::python::list l_xy;
	  		for (it_inner = xy.begin(); it_inner != xy.end(); ++it_inner)
	  		{
	  			l_xy.append(*it_inner);
	  		}

	    	boost::python::object module_ = boost::python::import("interpolate");
	    	boost::python::object result_ = module_.attr("Interpolate2D")(l_values, l_xy);
	    	double result = boost::python::extract<double>(result_);

	    	Py_Finalize();

	    	return result;
    	}
    	catch (boost::python::error_already_set const&)
	   	{
	        PyErr_Print();
	        return 1.0;
	   	}
	}

	/* histogram drawing */
	
	void histogram2D::draw()
	{	
		// define histogram with options
		Int_t font = 42;
		gStyle->SetOptStat(0);
		gStyle->SetTextFont(font);
		gStyle->SetFrameLineWidth(3);
		gStyle->SetHistLineWidth(3);
		gStyle->SetLineWidth(3);
		gStyle->SetTitleXOffset (1.);
		gStyle->SetTitleYOffset (1.);
		gStyle->SetLabelSize(0.05, "XYZ");
		gStyle->SetLabelFont(font, "XYZ");
		gStyle->SetTitleSize(0.06, "XYZ");
		gStyle->SetTitleFont(font, "XYZ");
		gStyle->SetTitleX(0.15);

		// create a canvas with options
		TCanvas* canvas = new TCanvas("", "");
		canvas->SetTicks(1, 1);
		canvas->SetBottomMargin(0.15);
		canvas->SetLeftMargin(0.15);
		canvas->SetRightMargin(0.25);
		canvas->SetLogx(0);
		canvas->SetLogy(0);
		canvas->SetLogz(0);

		// set axes' labels
		std::string labels = ";" + x_label + ";" + y_label;

  		// construct 2D histogram
	    TH2D* hist2D = new TH2D("", labels.c_str(), xbins, xmin, xmax, ybins, ymin, ymax);
	    
	    // draw 2D histogram dependent on which kind of sample has been specified
	    if (sample_xy.size() > 0) // draw for xy sample
	    {
			// fill 2D histogram with xy list values
			for (unsigned int i = 0; i < sample_xy.size(); i++)
			{
				std::vector< std::vector<double> > sample = sample_xy[i];
				for (unsigned int j = 0; j < sample.size(); j++)
				{
					std::vector<double> xy = sample[j];
					hist2D->Fill(xy[0], xy[1]);			
				}
			}
		}
	    else if (sample_xyz.size() > 0) // draw xyz sample
	    {
			// convert list_xyz to map
			std::map<std::list<double>, double> map_values;
			std::list< std::list<double> > list_values; 
			for (unsigned int i = 0; i < sample_xyz.size(); i++)
			{
				std::vector<double> xyz = sample_xyz[i];
				std::list<double> coord = {xyz[0], xyz[1]};
				map_values.insert(std::make_pair(coord, xyz[2]));
				list_values.push_back({xyz[0], xyz[1], xyz[2]});		
			}
			
			// fill 2D histogram with xyz values
			double value;
			double xbin_size = (xmax - xmin) / xbins;
			double ybin_size = (ymax - ymin) / ybins;

			for (double x = (xmin + xbin_size / 2); x < xmax; x += xbin_size)
			{
				for (double y = (ymin + ybin_size / 2); y < ymax; y += ybin_size)
				{	
					std::list<double> coord = {x, y};
					std::map<std::list<double>, double>::iterator it;
					it = map_values.find(coord);

					if (it != map_values.end())
						value = map_values[coord];
					else    				
						value = Interpolate2D(list_values, coord);

					hist2D->Fill(x, y, value);
				}
			}
		}
		else // no sample specified
		{
			std::cout << "WARNING: Trying to draw 2D histogram but no sample specified." << std::endl;
			delete hist2D;
			delete canvas;
			return;
		}

		// draw histogram
  		hist2D->Draw("colz");

  		// export histogram
		canvas->Print((ps_title + ".png").c_str());
		
		// delete pointers to canvas and histogram
		delete hist2D;
		delete canvas;
	}

/* NAMESPACE */
}
