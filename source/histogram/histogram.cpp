/* Histogram class
 *
 * 
*/

#include "histogram.h"


namespace analysis 
{

	/* con & destructor */

	histogram::histogram() 
	{
		nbins = 100;
		hmin = 0;
		hmax = 100;
		auto_range = true;

		ps_title = "new histogram";
    	leg_title = "";
    	x_label = ""; 
    	y_label = "";
		is_normalised = false;
		is_stacked = false;
		has_logy = false;
	}

	histogram::~histogram()
	{
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

	void histogram::set_logy(bool on)
	{
		has_logy = on;
	}

	void histogram::set_normalized(bool on)
	{
		is_normalised = on;
	}

	void histogram::set_stacked(bool on)
	{
		is_stacked = on;
	}

	/* histogram data */

	void histogram::add_sample(const std::vector<double> & sample, const std::string & name, double weight)
	{
		sample_list.push_back(sample);
		sample_names.push_back(name);
		sample_weights.push_back(weight);
	}

	/* histogram drawing */

	void histogram::draw()
	{	
		// create a canvas with a random name		
		std::random_device rd;
		TCanvas* canvas = new TCanvas(boost::lexical_cast<std::string>(rd()).c_str(), "");

		// options 
		Int_t font = 42;
		gStyle->SetOptStat(0);
		gStyle->SetTextFont(font);
		gStyle->SetFrameLineWidth(3);
		gStyle->SetHistLineWidth(3);
		gStyle->SetLineWidth(2);
		gStyle->SetTitleXOffset (1.);
		gStyle->SetTitleYOffset (1.);
		gStyle->SetLabelSize(0.05,"XYZ");
		gStyle->SetLabelFont(font, "XYZ");
		gStyle->SetTitleSize(0.055,"XYZ");
		gStyle->SetTitleFont(font, "XYZ");
		gStyle->SetEndErrorSize(0);

		canvas->SetTicks(1,1);
		canvas->SetBottomMargin(0.2);
		canvas->SetLeftMargin(0.2);
		canvas->SetLogy(0);
		if (has_logy)
			canvas->SetLogy(1);

		// set colors and line shapes
		Color_t kTransBlack = 1700, kTransRed = 1701, kTransBlue = 1702, kTransGreen = 1703, kTransYellow = 1704, kTransPink = 1705;
 		std::vector<Color_t> colors = {kBlack, kRed, kBlue, kGreen, kYellow, kPink};
		std::vector<Color_t> colors_trans = {kTransBlack, kTransRed, kTransBlue, kTransGreen, kTransYellow, kTransPink};
		for (unsigned int i = 0; i < colors.size(); ++i)
		{
			TColor *c = gROOT->GetColor(colors_trans[i]);
			if (!c)
				c = new TColor(colors_trans[i], 0.0, 0.0, 0.0);
				
			TColor *copy = gROOT->GetColor(colors[i]);
			Float_t r, g, b;
			copy->GetRGB(r, g, b);
			c->SetRGB(r, g, b);
			c->SetAlpha(0.5);
		}

		// set axes' labels 
		std::string labels = ";" + x_label + ";" + y_label;

	    // variables for setting the axes scales 
  		double max = -100;
  		double ymax;
	
  		// construction of histogram collection
  		const int nsamples = sample_list.size();
	    TH1D* hist[nsamples];

	    for (int iprc = 0; iprc < nsamples; ++iprc)
	    {
			std::vector<double> list = sample_list[iprc];

	    	// find max and min entry for auto range option.
			if (auto_range && !list.empty())
			{
				auto result = std::minmax_element(list.begin(), list.end());
				double min = *result.first;
				double max = *result.second;
				if (min < hmin) hmin = min;
				if (max > hmax) hmax = max;
			}

			// construct one histogram per sample
			hist[iprc] = new TH1D("", labels.c_str(), nbins, hmin, hmax);

			// read the sample and fill the histogram
			for (unsigned int i = 0; i < list.size(); i++)
				hist[iprc]->Fill(list[i], sample_weights[iprc]);

			if (is_normalised)
			{
				Double_t norm = sample_weights[iprc];
				Double_t scale = norm / hist[iprc]->Integral();
				hist[iprc]->Scale(scale);
			}

	    	hist[iprc]->SetLineColor(colors[iprc]);
    		hist[iprc]->SetLineWidth(2);
    		ymax = hist[iprc]->GetMaximum();
    		if(ymax > max) max = ymax;
	    }

	    // set the y-axis scale
		if (has_logy)
			max *= 3;
		else 
			max *= 1.1;
  		hist[0]->SetMaximum(max);

		// make legend at position
		double x1 = 0.65;
		double x2 = 0.89;
		double y1 = 0.87;
		double y2 = 0.64;
		TLegend *legend = new TLegend(x1, y1, x2, y2, leg_title.c_str());

		for (int iprc = 0; iprc < nsamples; iprc++)
    		legend->AddEntry(hist[iprc], sample_names[iprc].c_str());

  		// set legend cosmetics 
		legend->SetBorderSize(0);
		legend->SetLineWidth(0.0);
		legend->SetMargin(0.3);
		legend->SetTextSize(0.045);
		legend->SetTextFont(42);
		legend->SetFillColor(0);

		// draw histograms: stacked or not
		THStack stack(boost::lexical_cast<std::string>(rd()).c_str(), labels.c_str());
  		for (int iprc = 0; iprc < nsamples; ++iprc)
  		{
			if (is_stacked)
				hist[iprc]->SetFillColor(colors_trans[iprc]);

			stack.Add(hist[iprc]);
  		}
		if (is_stacked)
			stack.Draw();
		else
			stack.Draw("nostack");

  		// draw legend 
  		legend->Draw("same");

  		// export .png
		canvas->Print((ps_title + ".png").c_str());

		// delete pointers to canvas, histograms and legend
		for (int iprc = 0; iprc < nsamples; ++iprc)
  		{
			delete hist[iprc];
  		}
		delete legend;
		delete canvas;
	}

/* NAMESPACE */
}
