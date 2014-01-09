/* 2D Histogram class
 *
 * 
*/

#include "histogram2D.h"

namespace analysis {

	// ========================================== // 
	//	       	 	Class constructor	   		  //
	// ========================================== //
	histogram2D::histogram2D() 
	{
		
		xbins = 100;
		ybins = 100;
		xmin = 0;
		xmax = 100;
		ymin = 0;
		ymax = 100;

		ps_title = "new 2D histogram";
    	hist_title = "";
    	x_label = ""; 
    	y_label = "";
    	infile = "";
	}

	// histogram2D::~histogram2D()
	// {
	// }

	// ========================================== // 
	//	        Set Histogram Options 		   	  //
	// ========================================== //
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

	void histogram2D::set_ps_title(std::string ps)
	{
		ps_title = ps;
	}

	void histogram2D::set_hist_title(std::string title)
	{
		hist_title = title;
	}

	void histogram2D::set_x_label(std::string x)
	{
		x_label = x;
	}

	void histogram2D::set_y_label(std::string y)
	{
		y_label = y;
	}

	void histogram2D::add_sample(std::string sample)
	{
		infile = sample;
	}

	void histogram2D::set_palette(std::string name, const Int_t NCont)
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
    	
    	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    	gStyle->SetNumberContours(NCont);
	}


	// ========================================== // 
	//	        Interpolating functions		   	  //
	// ========================================== //
	double histogram2D::Interpolate2D(std::list< std::list<double> > values, std::list<double> xy)
	{
		namespace py = boost::python;
		
		try
    	{
			Py_Initialize();

	    	PyRun_SimpleString("import sys, string, random");
	    	PyRun_SimpleString("sys.path.append('../source/histogram2D')");
	    	// PyRun_SimpleString("sys.path.append(\".\")");

	    	py::list l_values;
	  		typename std::list< std::list<double> >::const_iterator it_outer;
	  		typename std::list<double>::const_iterator it_inner;

	  		for (it_outer = values.begin(); it_outer != values.end(); ++it_outer)
	  		{
	  			py::list l_inner;
	  			for (it_inner = (*it_outer).begin(); it_inner != (*it_outer).end(); ++it_inner)
	  			{
	  				l_inner.append(*it_inner);
	  			}
	  			l_values.append(l_inner);
	  		}  

	  		py::list l_xy;
	  		for (it_inner = xy.begin(); it_inner != xy.end(); ++it_inner)
	  		{
	  			l_xy.append(*it_inner);
	  		}

	    	py::object module_ = py::import("interpolate");
	    	py::object result_ = module_.attr("Interpolate2D")(l_values,l_xy);
	    	double result = py::extract<double>(result_);

	    	Py_Finalize();

	    	return result;
    	}
    	catch (py::error_already_set const&)
	   	{
	        PyErr_Print();
	        return 1.0;
	   	}
	}


	// ========================================== // 
	//	        	Draw Histogram   			  //
	// ========================================== // 
	void histogram2D::draw()
	{	
		// ========================================== // 
		//	        	Define Histogram  	 		  //
		// ========================================== // 

		//=== Options ===//
		Int_t iFont=42;
		//gROOT->SetStyle("Plain");
		gStyle->SetOptStat(0);
		gStyle->SetTextFont(iFont);
		gStyle->SetFrameLineWidth(3);
		gStyle->SetHistLineWidth(3);
		gStyle->SetLineWidth(3);
		gStyle->SetTitleXOffset (1.);
		gStyle->SetTitleYOffset (1.);
		gStyle->SetLabelSize(0.05,"XYZ");
		gStyle->SetLabelFont(iFont, "XYZ");
		gStyle->SetTitleSize(0.06,"XYZ");
		gStyle->SetTitleFont(iFont, "xyz");
		gStyle->SetTitleX(0.15);

		TCanvas* c1;
		c1 = new TCanvas("","");

		c1->SetTicks(1,1);
		c1->SetBottomMargin(0.15);
		c1->SetLeftMargin(0.15);
		c1->SetRightMargin(0.25);
		c1->SetLogx(0);
		c1->SetLogy(0);
		c1->SetLogz(0);

		//=== Set labels ===//
		std::stringstream labels;
	    labels << hist_title << ";" << x_label << ";" << y_label;

  		//=== Construction of 2D histogram ===//
	    TH2D* hist2D;
	    hist2D = new TH2D("",labels.str().c_str(),xbins,xmin,xmax,ybins,ymin,ymax);

	    /* read values */
	    double val1, val2, val3;
	    std::list<double> coord;
	    std::map<std::list<double>, double> map_values;
	    std::list<double> list_tmp;
	    std::list< std::list<double> > list_values; 

	    ifstream ifs;
  		ifs.open(infile.c_str(), std::ios::in);

  		while(true){
  			ifs >> val1;
    		if( ifs.eof() ) break;
    		ifs >> val2 >> val3;
  			
			/* map construction */
  			coord.push_back(val1);
  			coord.push_back(val2);
  			map_values.insert(std::make_pair(coord,val3));
  			coord.clear();

  			/* list construction */
  			list_tmp.push_back(val1);
  			list_tmp.push_back(val2);
  			list_tmp.push_back(val3);
  			list_values.push_back(list_tmp);
  			list_tmp.clear();
  		}

  		/* fill 2D histogram */
  		double value;
  		double xbin_size = (xmax-xmin)/(xbins);
  		double ybin_size = (ymax-ymin)/(ybins);

  		// std::cout <<"x binsize:"<< xbin_size << std::endl;
  		// std::cout <<"y binsize:"<< ybin_size << std::endl;

	    for (double x=(xmin+xbin_size/2); x<xmax; x+=xbin_size)
	    {
	    	for (double y=(ymin+ybin_size/2); y<ymax; y+=ybin_size)
	    	{	
	    		coord.push_back(x);
  				coord.push_back(y);
	    		std::map<std::list<double>, double>::iterator it;
	    		it = map_values.find(coord);

	    		if(it!=map_values.end())
	    		{
	    			value = map_values[coord];
	    		}
    			else
    			{
    				// std::cout << "interpolation... ";
    				value = Interpolate2D(list_values,coord);
    			}

    			// std::cout << "(x,y,value)= " << x <<","<< y <<","<< value << std::endl;

        		hist2D->Fill(x,y,value);
        		coord.clear();
	    	}
	    }

		// ========================================== // 
		//		        Draw and Export	     		  //
		// ========================================== //

		//=== Draw histograms ===//
  		hist2D->Draw("colz");

  		//=== Export .ps ===//
		c1->Print(ps_title.c_str());
	}

	// delete c1;

/* NAMESPACE */
}
