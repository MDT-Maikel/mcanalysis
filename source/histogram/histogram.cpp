/* Histogram class
 *
 * Plots a simple 1D histogram with options.
*/

#include "histogram.h"


/* NAMESPACE */
namespace analysis 
{

	/* con & destructor */

	histogram::histogram() 
	{
		nbins = 100;
		// set hmin and hmax to unrealistic values
		// which are then automatically set later
		hmin = 10000;
		hmax = -10000;
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
		if (min >= max)
			return;
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
		Color_t kTransBlack = 1700, kTransRed = 1701, kTransBlue = 1702, kTransGreen = 1703, kTransMagenta = 1704, kTransCyan = 1705;
 		std::vector<Color_t> colors = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan};
		std::vector<Color_t> colors_trans = {kTransBlack, kTransRed, kTransBlue, kTransGreen, kTransMagenta, kTransCyan};
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
	    std::vector<TH1D*> hist;

	    for (int i = 0; i < sample_list.size(); ++i)
	    {
			std::vector<double> list = sample_list[i];

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
			hist.push_back(new TH1D("", labels.c_str(), nbins, hmin, hmax));

			// read the sample and fill the histogram
			for (unsigned int j = 0; j < list.size(); j++)
				hist[i]->Fill(list[j], sample_weights[i]);

			if (is_normalised)
			{
				Double_t norm = sample_weights[i];
				Double_t scale = norm / hist[i]->Integral();
				hist[i]->Scale(scale);
			}

	    	hist[i]->SetLineColor(colors[i]);
    		hist[i]->SetLineWidth(2);
    		ymax = hist[i]->GetMaximum();
    		if (ymax > max) 
				max = ymax;
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

		for (int i = 0; i < sample_list.size(); ++i)
    		legend->AddEntry(hist[i], sample_names[i].c_str());

  		// set legend cosmetics 
		legend->SetBorderSize(0);
		legend->SetLineWidth(0.0);
		legend->SetMargin(0.3);
		legend->SetTextSize(0.045);
		legend->SetTextFont(42);
		legend->SetFillColor(0);

		// draw histograms: stacked or not
		THStack stack(boost::lexical_cast<std::string>(rd()).c_str(), labels.c_str());
		for (int i = 0; i < sample_list.size(); ++i)
  		{
			if (is_stacked)
				hist[i]->SetFillColor(colors_trans[i]);

			stack.Add(hist[i]);
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
		for (int i = 0; i < hist.size(); ++i)
  		{
			delete hist[i];
  		}
		delete legend;
		delete canvas;
	}

/* NAMESPACE */
}
