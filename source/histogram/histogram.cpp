/* Histogram class
 *
 * 
*/

#include "histogram.h"

namespace analysis {

	/* con & destructor */

	histogram::histogram() 
	{
		c1 = new TCanvas("", "");
		nbins = 100;
		hmin = 0;
		hmax = 100;
		auto_range = true;

		ps_title = "new histogram";
    	leg_title = "";
    	x_label = ""; 
    	y_label = "";
		is_normalised = false;
	}

	histogram::~histogram()
	{
		delete c1;
	}

	/* histogram options */
	
	void histogram::set_title(std::string ps)
	{
		ps_title = ps;
	}

	void histogram::set_bins(double nb)
	{
		nbins = nb;
	}

	void histogram::set_range(double min, double max)
	{
		hmin = min;
		hmax = max;
		auto_range = false;
	}

	void histogram::set_x_label(std::string x)
	{
		x_label = x;
	}

	void histogram::set_y_label(std::string y)
	{
		y_label = y;
	}

	void histogram::set_leg_title(std::string title)
	{
		leg_title = title;
	}

	void histogram::add_sample(const std::vector<double> & sample, const std::string & name, const std::string & line_color)
	{
		sample_list.push_back(sample);
		sample_colors.push_back(line_color);
		sample_names.push_back(name);
	}

	void histogram::normalize()
	{
		is_normalised = true;
	}

	/* draw histogram */

	void histogram::draw()
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
		gStyle->SetLineWidth(2);
		gStyle->SetTitleXOffset (1.);
		gStyle->SetTitleYOffset (1.);
		gStyle->SetLabelSize(0.05,"XYZ");
		gStyle->SetLabelFont(iFont, "XYZ");
		gStyle->SetTitleSize(0.055,"XYZ");
		gStyle->SetTitleFont(iFont, "xyz");
		gStyle->SetEndErrorSize(0);

		c1->SetTicks(1,1);
		c1->SetBottomMargin(0.2);
		c1->SetLeftMargin(0.2);
		c1->SetLogy(0);

		/* set colors and line shapes */
	    Color_t col[4] = {kBlack, kRed, kBlue, kGreen};
	    int line_shape[4] = {1, 7, 3, 4};

		/* set axes' labels */
		std::stringstream labels;
	    labels << ";" << x_label << ";" << y_label;

	    /* variables for setting the axes scales */
  		double max = -100;
  		double ymax;
	
  		/* construction of histogram collection */
  		const int nsamples = sample_list.size();
	    TH1F* hist[nsamples];

	    for (int iprc = 0; iprc < nsamples; ++iprc)
	    {
			std::vector<double> list = sample_list[iprc];

	    	// find max and min entry for auto range option.
			if (auto_range)
			{
				auto result = std::minmax_element(list.begin(), list.end());
				double min = *result.first;
				double max = *result.second;
				if (min < hmin) hmin = min;
				if (max > hmax) hmax = max;
			}

			// construct one histogram per sample
			hist[iprc] = new TH1F("", labels.str().c_str(), nbins, hmin, hmax);

			// read the sample and fill the histogram
			for (unsigned int i = 0; i < list.size(); i++)
			{
				hist[iprc]->Fill(list[i]);
			}

			if (is_normalised)
			{
				Double_t norm = 1;
				Double_t scale = norm / hist[iprc]->Integral();
				hist[iprc]->Scale(scale);
			}

	    	// for(int i=0; i<nbins; i++) hist[iprc]->SetBinContent(i, 1./(i+1));

	    	hist[iprc]->SetLineColor(col[iprc]);
    		hist[iprc]->SetLineWidth(3);
    		hist[iprc]->SetLineStyle(line_shape[iprc]);
    		ymax = hist[iprc]->GetMaximum();
    		if(ymax > max) max = ymax;
	    }

	    //=== Set the y-axis scale ===//
  		hist[0]->SetMaximum(1.1 * max);


		/* make a legend */

		// make legend at position
		double x1 = 0.65;
		double x2 = 0.89;
		double y1 = 0.87;
		double y2 = 0.64;
		TLegend *leg = new TLegend(x1, y1, x2, y2, leg_title.c_str());

		for (int iprc = 0; iprc < nsamples; iprc++)
    		leg->AddEntry(hist[iprc], sample_names[iprc].c_str());

  		//=== Set Legend cosmetic ===//
		leg->SetBorderSize(0);
		leg->SetLineWidth(0.0);
		leg->SetMargin(0.3);
		leg->SetTextSize(0.045);
		leg->SetTextFont(42);
		leg->SetFillColor(0);


		// ========================================== // 
		//		        Draw and Export	     		  //
		// ========================================== //

		//=== Draw histograms ===//
  		for (int iprc = 0; iprc < nsamples; ++iprc)
  		{
  			const char* option = "";
  			if (iprc > 0)
  			{
  				option = "same";
  			}
  			hist[iprc]->Draw(option);
  		}

  		//=== Draw legend ===//
  		leg->Draw("same");

  		//=== Export .ps ===//
		c1->Print(ps_title.c_str());
	}

/* NAMESPACE */
}
