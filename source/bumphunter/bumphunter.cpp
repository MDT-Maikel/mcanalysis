/* BumpHunter class
 *
 * implementation of the BumpHunter algorithm based on: 
 * 	C++code: https://twiki.cern.ch/twiki/bin/view/Main/BumpHunter by W. Buttinger (2012)
 * 	article: http://arxiv.org/pdf/1101.0390v2.pdf by G. Choudalakis (2011)
 * 
*/

#include "bumphunter.h" 


/* NAMESPACE */
namespace analysis
{
	
	/* con & destructor */
	
	bumphunter::bumphunter(TH1* mc, TH1* data) : TObject() 
	{ 
		init();
		SetChannelDistribution(mc, data);
	}
	
	bumphunter::bumphunter() : TObject() 
	{ 
		init(); 
	};
	
	bumphunter::~bumphunter() 
	{
	}

	/* initialization */

	void bumphunter::init() 
	{
		// general settings
		hunt_name = "";
		hunt_folder = "";
		
		// adding data
		current_channel = 0;
		
		// hunt properties: general 
		fBinModel = BUMP_POISSON_GAMMA; 
		fTestStatisticType = BUMPHUNTER;
		fnPseudo = 0;	
		fMultiChannelBumpOverlapFactor = 1.;
		
		/* hunt properties: bump window */
		fMinWindowSize = -1;
		fMaxWindowSize = -1;
		fWindowSizeStep = -1;

		// hunt properties: results 
		local_pvalue = 1;
		global_pvalue = 1;	
		
		// statistics 
		pRand.SetSeed(123345);
		//r.SetSeed(431155);
	}
	
	/* general settings */
	
	void bumphunter::set_name(std::string n)
	{
		hunt_name = n;
	}

	void bumphunter::set_folder(std::string f)
	{
		hunt_folder = f;
	}
	
	/* data adding */
	
	//specify the background distribution (e.g. the poisson means and the uncertainty on these means)
	void bumphunter::SetChannelDistribution(TH1* hist, TH1* data) 
	{ 
		//fHistBack = hist; fHistData=data;
		if (channel_list.size() == 0)
			channel_list.resize(1);
		channel_list[0].hist_bkg = hist; 
		channel_list[0].hist_data = data;
		if (channel_list[0].hist_pseudodata) 
			delete channel_list[0].hist_pseudodata;
		channel_list[0].hist_pseudodata = (TH1*)data->Clone(Form("pseudo_%s",data->GetName())); 
		channel_list[0].hist_pseudodata->SetDirectory(0);
		channel_list[0].hist_pseudodata->Reset();
		SetChannel(0);
		EvaluateSearchPattern(); 
	}

	void bumphunter::AddChannelDistribution(TH1* hist, TH1* data) 
	{ 
		channel t;
		t.hist_bkg = hist; 
		t.hist_data = hist; 
		t.hist_pseudodata = (TH1*) data->Clone(Form("pseudo_%s", data->GetName())); 
		t.hist_pseudodata->SetDirectory(0);
		t.hist_pseudodata->Reset();
		channel_list.push_back(t);
		SetChannel(channel_list.size()-1);
		EvaluateSearchPattern();
	}
	
	void bumphunter::SetChannel(Int_t in)
	{ 
		current_channel = in;
	}
	
	Int_t bumphunter::GetNChannels() 
	{
		return channel_list.size(); 
	}
	
	/* hunt properties: general */
	
	void bumphunter::SetBinModel(Int_t model) 
	{ 
		if (model != 0 && model != 1 && model != 2 && model != 3)
			return;
		fBinModel = model;
	}
	
	void bumphunter::SetTestStatisticType(Int_t t) 
	{
		fTestStatisticType = t;
	}

	void bumphunter::SetNPseudoExperiments(Int_t n) 
	{
		fnPseudo = n;
	}
			
	void bumphunter::SetBumpOverlapFactor(Double_t in) 
	{ 
		fMultiChannelBumpOverlapFactor = in; 
	}
	
	/* hunt properties: search region */
	
	void bumphunter::SetSearchRegion(Double_t low, Double_t high)
	{
		channel_list[current_channel].search_min = low; 
		channel_list[current_channel].search_max = high; 
		EvaluateSearchPattern(); 
	}
	
	Double_t bumphunter::GetSearchLowEdge() 
	{
		channel& curr_channel = channel_list[current_channel];
		if (!curr_channel.hist_bkg) 
			return 0;
		if (curr_channel.search_min != curr_channel.search_min)
		{
			// look for first non-zero mc prediction bin
			int i = 1;
			while (curr_channel.hist_bkg->GetBinContent(i) == 0 && curr_channel.hist_bkg->GetBinError(i) == 0 && i <= curr_channel.hist_bkg->GetNbinsX())
				i++;
			if (i > curr_channel.hist_bkg->GetNbinsX()) 
			{
				Error("GetSearchLowEdge", "All bins are empty!?"); 
				return 0;
			}
			return curr_channel.hist_bkg->GetBinLowEdge(i);
		}
		return curr_channel.hist_bkg->GetBinLowEdge(curr_channel.hist_bkg->FindFixBin(curr_channel.search_min));
	}

	Double_t bumphunter::GetSearchHighEdge() 
	{
		channel& curr_channel = channel_list[current_channel];
		if (!curr_channel.hist_bkg) 
			return 0;
		if (curr_channel.search_max != curr_channel.search_max) 
		{
			// look for first zero mc prediction bin after low edge
			int i = curr_channel.hist_bkg->FindFixBin(GetSearchLowEdge());
			while (curr_channel.hist_bkg->GetBinContent(i) != 0 && i <= curr_channel.hist_bkg->GetNbinsX())
				i++;
			if (i == curr_channel.hist_bkg->FindFixBin(GetSearchLowEdge())) 
			{
				Error("GetSearchHighEdge", "All bins are empty!?"); 
				return 0;
			}
			return curr_channel.hist_bkg->GetBinLowEdge(i);
		}
		return curr_channel.hist_bkg->GetBinLowEdge(curr_channel.hist_bkg->FindFixBin(curr_channel.search_max) + 1);
	}
	
	/* hunt properties: bump window */
	
	void bumphunter::SetMinWindowSize(Int_t in) 
	{
		fMinWindowSize = in;
		EvaluateSearchPattern();
	}
	
	void bumphunter::SetMaxWindowSize(Int_t in) 
	{
		fMaxWindowSize = in;
		EvaluateSearchPattern();
	}
	
	void bumphunter::SetWindowStepSize(Int_t in) 
	{
		fWindowSizeStep = in;
		EvaluateSearchPattern();
	}
	
	Double_t bumphunter::GetMinWindowSize() 
	{
		channel& curr_channel = channel_list[current_channel];
		if (!curr_channel.hist_bkg)
			return 0;
		// default is twice the average bin width
		if (fMinWindowSize < 0) 
			return 1;
		return fMinWindowSize;
	}

	Double_t bumphunter::GetMaxWindowSize() 
	{
		channel& curr_channel = channel_list[current_channel];
		if (!curr_channel.hist_bkg)
			return 0;
		// default is half the search window
		if (fMaxWindowSize < 0)
			return (GetSearchHighEdge() - GetSearchLowEdge()) /2.;
		return fMaxWindowSize;
	}

	Double_t bumphunter::GetWindowSizeStep() 
	{
		channel& curr_channel = channel_list[current_channel];
		if (!curr_channel.hist_bkg) 
			return 0;
		// default is the average bin width in the search region, or 1 bin if using bins for the window size
		if (fWindowSizeStep < 0) 
			return 1;
		return fWindowSizeStep;
	}
	
	Double_t bumphunter::GetBumpLowEdge() 
	{
		return channel_list[current_channel].bump_window_min;
	}
	
	Double_t bumphunter::GetBumpHighEdge() 
	{ 
		return channel_list[current_channel].bump_window_max; 
	}

	/* hunt properties: results */
	
	double bumphunter::get_global_pvalue()
	{
		return global_pvalue; 
	}

	double bumphunter::get_local_pvalue()
	{
		return local_pvalue;
	}
	
	double bumphunter::get_global_sigma()
	{
		return GetZValue(global_pvalue); 
	}

	double bumphunter::get_local_sigma()
	{
		return GetZValue(local_pvalue);
	}
	
	/* statistics */
	
	Double_t bumphunter::GetPoissonPValue(Double_t nobs, Double_t e) 
	{
		if (nobs > e) // excess
		{ 
			double p = ROOT::Math::inc_gamma_c(nobs, e);
			if (p == 1. && nobs > 100.) // excess very excessive, and have a large nobs so no big problem summing one extra term
				return ROOT::Math::inc_gamma(nobs, e);
			else 
				return 1 - p;
		} 
		return ROOT::Math::inc_gamma_c(nobs + 1, e); // deficit
	}

	Double_t bumphunter::GetGaussianPValue(Double_t nobs, Double_t e, Double_t sigma) 
	{
		if (nobs > e) 
			return ROOT::Math::gaussian_cdf_c(nobs, sigma, e);
		return ROOT::Math::gaussian_cdf(nobs, sigma, e);
	}

	// return the p-value (and associated zValue = number of sigma) of seeing nobs "events" under the model a poisson distribution 
	// with the value of the mean parameter being itself uncertain. Assumption is that the mean parameter is distributed as a gamma
	// density, with the expectation (mean) of the gamma function = e and the variance of the gamma density equal to err^2
	Double_t bumphunter::GetPoissonConvGammaPValue(Double_t nobs, Double_t e, Double_t err) 
	{
		//parameters (a,b) of gamma density can be written in terms of the expectation and variance:
		//Gammma(x, a,b);   - note this is not equal to gamma(x) or gamma(a,b), which are different functions
		double b = e / (err * err); // = E/V
		double a = e * b; // = E^2/V

		double pval = 0.;
		//decide if we should ignore systematics or not 
		if (a > 100 * nobs) 
		{
			//stat error is big enough to ignore the syst error
			//considering only stat error means the p-value is given by:
			// (nData>nMC): pval = sum(x = nData->inf, Poisson(x,nMC)) = 1 - sum(x = 0->nData-1,Poisson(x,nMC))
			// (nData<=nMC): pval = sum(x = 0->nData, Poisson(x,nMC))
			// But sum(x = 0->nData-1,Poisson(x,nMC)) = gamma(nData,nMC)/gamma(nData); <---- see maths websites
			// so we have:
			// (nData>nMC): pval = 1 - gamma(nData,nMC)/gamma(nData);
			// (nData<=nMC): pval = gamma(nData+1,nMC)/gamma(nData);
			// .....And ROOT provides gamma(a,b)/gamma(a) = ROOT::Math::inc_gamma_c(a,b)
			//pval = (nObs<=E) ? ROOT::Math::inc_gamma_c(nObs+1,E) : (1. - ROOT::Math::inc_gamma_c(nObs,E));
			pval = GetPoissonPValue(nobs, e);
		} 
		else 
		{
			//use recursive formula to solve:
			// (nData>nMC): pval = 1 - sum(x=0->nData-1, Integral(y=0->inf, Poisson(x,y)Gamma(y,a,b) dy))
			// (nData<=nMC): pval = sum(x=0->nData, Integral(y=0->inf, Poisson(x,y)Gamma(y,a,b) dy))
			//i.e. integrating out the unknown parameter y
			// Recursive formula: P(n;A,B) = P(n-1;A,B) (A+n-1)/(n*(1+B))
			unsigned stop = nobs;
			if (nobs > e) 
				--stop;
			double sum = 0;
			if (a > 100) 
			{
				// NB: must work in log-scale otherwise troubles!
				double log_prob = a * log(b / (1 + b));
				sum = exp(log_prob); // P(n=0)
				for (unsigned u = 1; u <= stop; ++u)
				{
					log_prob += log((a + u - 1) / (u * (1 + b)));
					sum += exp(log_prob);
				}
			} 
			else 
			{
				double p0 = pow(b / (1 + b), a); // P(0;A,B)
				double pLast = p0;
				sum = p0;
				for (unsigned k = 1; k <= stop; ++k) 
				{
					double p = pLast * (a + k - 1) / (k * (1 + b));
					sum += p;
					pLast = p;
				}
			}
			pval = (nobs > e) ? 1 - sum : sum; 
		} 
		//bool overflow(false);
		//zValue = GetZValue(pValue,overflow) //large z-values correspond to small p-values.... i.e. significant diff
		//zValue = (nObs<E) ? -1.*zValue : zValue; //flip the z-values of deficits
		return pval;
	}

	Double_t bumphunter::GetZValue(Double_t pvalue, bool& overflow) 
	{
		if (pvalue < 0.000000001) 
		{
			overflow = true;
			return 6.; 
		}
		return sqrt(2.) * TMath::ErfInverse(1. - 2. * pvalue);
	}
	
	Double_t bumphunter::GetZValue(Double_t pvalue) 
	{
		if (pvalue < 0.000000001) 
			return 6.; 
		return sqrt(2.) * TMath::ErfInverse(1. - 2. * pvalue);
	}
	
	/* bump hunting */
	
	Int_t bumphunter::run() 
	{
		// create canvas for plotting the results of this run
		TCanvas* fBumpHunterResultCanvas;
		fBumpHunterResultCanvas = new TCanvas("bhCanvas", "BumpHunter Results", 500, 600);
		fBumpHunterResultCanvas->Divide(1,3);
		fBumpHunterResultCanvas->cd(1);

		// get current channel and set the bump window to zero
		channel& curr_channel = channel_list[current_channel];
		curr_channel.bump_window_min = 0.;
		curr_channel.bump_window_max = 0.;

		// evaluate bumphunter for data, automatically defaults to single channel 
		Double_t tobs = EvaluateMultiChannelTestStatistic(false);
		local_pvalue = exp(-tobs);
	
		if (curr_channel.bump_window_min == 0 && curr_channel.bump_window_max == 0) 
		{
			Error("Run","No bumps found?"); 
			return -1;
		}
		
		// loop over the channels to show most significant bumps 
		for (unsigned int i = 0; i < channel_list.size(); i++) 
		{
			bool overflow = false;
			double zval = GetZValue(channel_list[i].bump_pvalue, overflow);
			if(overflow) 
			{
				Info("Run","Channel #%d: Observed most significant bump in [%f,%f] (> %f local sigma)", i, channel_list[i].bump_window_min, channel_list[i].bump_window_max, zval);
			} 
			else 
			{
				Info("Run","Channel #%d: Observed most significant bump in [%f,%f] (%f local sigma)", i, channel_list[i].bump_window_min, channel_list[i].bump_window_max, zval);
			}
		}
		
		if (channel_list.size() > 1) 
		{
			bool overflow = false;
			double zval = GetZValue(local_pvalue, overflow);
			if(overflow) 
			{
				Info("Run","MultiChannel p-value (with OverlapFactor=%f): %f (> %f sigma)", fMultiChannelBumpOverlapFactor, local_pvalue, zval);
			}
			else
			{
				Info("Run","MultiChannel p-value (with OverlapFactor=%f): %f (%f sigma)", fMultiChannelBumpOverlapFactor, local_pvalue, zval);
			}
		}
		
		Double_t b1 = curr_channel.bump_window_min; 
		Double_t b2 = curr_channel.bump_window_max;
		
		TH1D* fBumpHunterStatisticPDF = new TH1D("bumpHunterPDF", "BumpHunter Test Statistic PDF", 100, 0, tobs * 1.5);
		fBumpHunterStatisticPDF->SetDirectory(0);
		fBumpHunterStatisticPDF->GetXaxis()->SetTitle("BumpHunter Statistic");
		fBumpHunterStatisticPDF->GetYaxis()->SetTitle("Probability Density");

		Int_t nDone = 0; //Int_t nGreater=0;

		global_pvalue = 1.0;

		std::vector<Double_t> trial;
		std::vector<Double_t> pValues;
		std::vector<Double_t> pValuesLow; 
		std::vector<Double_t> pValuesHigh;

		TH1D* nGreater = new TH1D("nGreater", "", 1, 0, 1);
		TH1D* nTotal = new TH1D("nTotal", "", 1, 0, 1);

	
	
		// 
   		Int_t nPseudo = fnPseudo;
		if (nPseudo == 0) // if no nPseudo set, then use the pvalue to estimate the number of pseudo needed
		{
			nPseudo = (local_pvalue < 0.0000001) ? 1000000 : 1. / local_pvalue; 
			if (nPseudo < 1000) 
				nPseudo = 1000;
		}

		Info("Run","Performing %d Pseudo-experiments....",nPseudo);
		std::cout << "|0%--------------------------------------------100%|" << std::endl;
		std::cout << "|" << std::flush;
		Int_t tickPoint = nPseudo / 50;
		Int_t graphPoint = nPseudo / 1000 + 1;
		while (/*!converged && */nDone < nPseudo) 
		{
			if ((nDone % tickPoint) == 0)
				std::cout << "*" << std::flush;
			Double_t tPseudo = EvaluateMultiChannelTestStatistic(true);
			fBumpHunterStatisticPDF->Fill(tPseudo);
			if (tPseudo > tobs)
				nGreater->Fill(0.5); //nGreater++;
			nDone++;
			nTotal->Fill(0.5);
			global_pvalue = (double) nGreater->GetEntries() / (double) nDone;
			if (nDone == nPseudo || nGreater->GetEntries() == 1 || (nDone % graphPoint) == 0) 
			{
				pValues.push_back(global_pvalue);
				trial.push_back(nDone);
				TGraphAsymmErrors f; 
				f.BayesDivide(nGreater,nTotal); 
				pValuesLow.push_back(f.GetErrorYlow(0));
				pValuesHigh.push_back(f.GetErrorYhigh(0)); // original: pValuesLow.push_back(f.GetErrorYhigh(0))
			}
		}
		std::cout << "|" << std::endl;

		TGraphAsymmErrors* fBumpHunterStatisticConvergenceGraph = new TGraphAsymmErrors(trial.size(), &trial[0], &pValues[0], &pValuesLow[0], &pValuesHigh[0]);

		delete nGreater; 
		delete nTotal;

		curr_channel.bump_window_min = b1;
		curr_channel.bump_window_max = b2;

		// draw the background and data with bump window onto the canvas
		curr_channel.hist_bkg->Draw();
		curr_channel.hist_bkg->GetXaxis()->SetRangeUser(GetSearchLowEdge(),GetSearchHighEdge());
		curr_channel.hist_data->Draw("same");
		TLine *bumpLineLeft = new TLine(curr_channel.bump_window_min, 0, curr_channel.bump_window_min, curr_channel.hist_bkg->GetMaximum());
		TLine *bumpLineRight = new TLine(curr_channel.bump_window_max, 0, curr_channel.bump_window_max, curr_channel.hist_bkg->GetMaximum());
		bumpLineLeft->SetLineColor(kRed);
		bumpLineLeft->SetLineStyle(2);
		bumpLineLeft->SetLineWidth(2);
		bumpLineLeft->Draw("same");
		bumpLineRight->SetLineColor(kRed);
		bumpLineRight->SetLineStyle(2);
		bumpLineRight->SetLineWidth(2);
		bumpLineRight->Draw("same");

		fBumpHunterResultCanvas->cd(2);
		fBumpHunterStatisticPDF->Scale(1. / fBumpHunterStatisticPDF->GetSumOfWeights());
		fBumpHunterStatisticPDF->Draw();
		TLine *l1 = new TLine(tobs, 0, tobs, fBumpHunterStatisticPDF->GetMaximum() / 2.);
		l1->SetLineColor(kBlue);
		l1->SetLineWidth(2);
		l1->Draw("same");
		fBumpHunterResultCanvas->cd(3);
		fBumpHunterStatisticConvergenceGraph->Draw("ALP");
		fBumpHunterStatisticConvergenceGraph->GetHistogram()->GetXaxis()->SetTitle("# Pseudo-experiments");
		fBumpHunterStatisticConvergenceGraph->GetHistogram()->GetYaxis()->SetTitle("P-Value");
		TLine *l = new TLine(fBumpHunterStatisticConvergenceGraph->GetHistogram()->GetXaxis()->GetXmin(),pValues[pValues.size()-1],fBumpHunterStatisticConvergenceGraph->GetHistogram()->GetXaxis()->GetXmax(),pValues[pValues.size()-1]);
		l->SetLineColor(kRed);
		l->Draw("same");
		TLatex n;
		n.SetNDC();n.SetTextFont(43);n.SetTextSize(14);n.SetTextColor(kRed);
		std::stringstream m;m.precision(3);
		bool overflow = false;
		
		if (global_pvalue > 0.) 
		{
			 double zVal = GetZValue(global_pvalue, overflow);
			 m << zVal << "#sigma (local ";
			 overflow=false;
			 zVal = GetZValue(local_pvalue, overflow);
			 if(overflow) m << "> ";
			 m << zVal << "#sigma)";
		}
		else
		{
			 double zVal = GetZValue(1./double(nPseudo),overflow);
			 m << " > " << zVal << "#sigma (local ";
			 overflow=false;
			 zVal = GetZValue(local_pvalue, overflow);
			 if(overflow) m << "> ";
			 m << zVal << "#sigma)";
		}
		n.DrawLatex(0.5, 0.6, m.str().c_str());

		if (global_pvalue > 0.) 
			Info("Run","global p = %f (%f sigma)",global_pvalue,GetZValue(global_pvalue,overflow));
		else //none of the pseudodata was bigger than our pvalue, so estimate the global p-value as 1/nPseudo 
		{
			global_pvalue = 1. / double(nPseudo);
			Info("Run","global p < %f (%f sigma)",1./double(nPseudo),GetZValue(1./double(nPseudo),overflow));			
		}
  
		// export the final plot to a png file
		fBumpHunterResultCanvas->Print((hunt_folder + hunt_name + ".png").c_str());
		delete fBumpHunterResultCanvas;
		delete fBumpHunterStatisticPDF;
		delete fBumpHunterStatisticConvergenceGraph;
		return 0;
	}
	

	/* TO SORT */

	//decide what sequence to use for the various windows for the current settings
	void bumphunter::EvaluateSearchPattern() 
	{
		channel& curr_channel = channel_list[current_channel];
		curr_channel.central_windows.clear();

		TH1* fHistBack = curr_channel.hist_bkg;


		Int_t startBin = fHistBack->FindFixBin(GetSearchLowEdge()); //first bin to use in search
		Int_t stopBin = fHistBack->FindFixBin(GetSearchHighEdge())-1; //last bin to use in search

		if (stopBin < startBin) 			
			return;

		Double_t windowSizeStepSize = GetWindowSizeStep();
		Double_t currentWindowSize = GetMinWindowSize();
		Double_t maxWsize = GetMaxWindowSize();

		bool hasGoodWindows = true;
		while (hasGoodWindows)
		{
			hasGoodWindows = false;
			Int_t firstBin = startBin;
			Int_t lastBin = firstBin + static_cast<int>(currentWindowSize - 0.5);
			while (lastBin <= stopBin) 
			{
				// check the current bin isn't too wide 
				if ((lastBin-firstBin+1) <= maxWsize) 
				{
					std::pair<Int_t,Int_t> p; 
					p.first = firstBin; 
					p.second=lastBin; 
					curr_channel.central_windows.push_back(p);
					hasGoodWindows = true;
				}
				// shift the bump window by half the current window size, but at least one
				firstBin += std::max(currentWindowSize / 2, 1.0);
				lastBin = firstBin + static_cast<int>(currentWindowSize - 0.5);
			}
			currentWindowSize += windowSizeStepSize;
		}
	}

	void bumphunter::PrintSearchPattern() 
	{
	   std::vector<std::pair<Int_t,Int_t> >& myWindows = channel_list[current_channel].central_windows;
	   for (std::vector<std::pair<Int_t,Int_t> >::iterator it = myWindows.begin(); it != myWindows.end(); ++it)
	   {
		  Info("PrintSearchPattern","%d to %d [%f,%f]", it->first, it->second, channel_list[current_channel].hist_bkg->GetBinLowEdge(it->first), channel_list[current_channel].hist_bkg->GetBinLowEdge(it->second + 1));
	   }
	}

	// TODO: try to improve this algorithm
	Double_t bumphunter::EvaluateTestStatistic(TH1* data, Bool_t printOut) 
	{
		//start at low edge
		//calculate pvalue in the given window
		//shift the window across max(first-bin-width in window,window/2)
		//recalculate pvalue
		//keep going until right edge of window is beyond high edge
		//increase window size by (highEdge-lowEdge)/(nBins in range from lowEdge to highEdge) - for equal bin widths this is just one bin
		//keep repeating all this until window size > maxWindowSize

		Double_t minPValue = 1.;

		channel& curr_channel = channel_list[current_channel];
		TH1* fHistBack = curr_channel.hist_bkg;

		std::vector<std::pair<Int_t,Int_t> >& myWindows = curr_channel.central_windows;

		for(std::vector<std::pair<Int_t,Int_t> >::iterator it = myWindows.begin(); it != myWindows.end(); ++it)
		{
			Double_t localPValue = 1.;        
			Double_t nObs = data->Integral(it->first, it->second);
			Double_t errExp(0.);
			Double_t nExp(0.);
			if (fBinModel == 3) 
			{
				//need to treat errors as correlated, so add up all the errors in the range 
				for (int i=it->first; i <= it->second; i++)
				{
					nExp += fHistBack->GetBinContent(i); 
					errExp += fHistBack->GetBinError(i);
				}
			}
			else
			{
				nExp = fHistBack->IntegralAndError(it->first,it->second,errExp);
			}

			if (fTestStatisticType == BUMPHUNTER && nObs < nExp) 
				localPValue = 1.; // dips are considered insignificant if only bumphunting
			else if (fTestStatisticType == DIPHUNTER && nObs > nExp) 
				localPValue = 1.; // bumps are considered insigificant if diphunting
			else if (fBinModel == 0) 
				localPValue = GetPoissonConvGammaPValue(nObs, nExp, errExp);
			else if (fBinModel == 1) 
				localPValue = GetPoissonPValue(nObs, nExp);
			else if (fBinModel == 2) 
				localPValue = GetGaussianPValue(nObs, nExp, errExp);
			else if (fBinModel == 3) 
				localPValue = GetPoissonConvGammaPValue(nObs, nExp, errExp);
				
			if (printOut) 
				Info("EvaluateTestStatistic", "search region: [%f,%f] gave pvalue of %f (mc=%f, data=%f)", fHistBack->GetBinLowEdge(it->first), fHistBack->GetBinLowEdge(it->second + 1), localPValue, nExp, nObs);
			if (localPValue < minPValue) 
			{
				if(printOut) 
					Info("EvaluateTestStatistic","smallest so far");
				minPValue = localPValue;
				curr_channel.bump_window_min = fHistBack->GetBinLowEdge(it->first);
				curr_channel.bump_window_max = fHistBack->GetBinLowEdge(it->second) + fHistBack->GetBinWidth(it->second);
				curr_channel.bump_pvalue = localPValue;
			}
		}

		return -log(minPValue);
	}

	Double_t bumphunter::EvaluateMultiChannelTestStatistic(Bool_t generatePseudo, Bool_t printOut) 
	{
		if (channel_list.size() == 1) 
			return EvaluateTestStatistic(((generatePseudo)? GenerateToyMC() : channel_list[0].hist_data));
		// loop over channels, requiring the the worst bump overlap to still be better than the overlapfactor 
		std::vector<Double_t> bumpOverlaps;
		Double_t commonWindowLow = -DBL_MAX; 
		Double_t commonWindowHigh = DBL_MAX;
		Double_t out = 0.;

		for (unsigned int i = 0; i < channel_list.size(); i++)
		{
			SetChannel(i);
			channel& currChannel = channel_list[current_channel];
			// can just add the nll of the pvalues to make the total test statistic
			out += ((generatePseudo) ? EvaluateTestStatistic(GenerateToyMC()) :  EvaluateTestStatistic(currChannel.hist_data));
			bumpOverlaps.push_back(1.);
			// update the common window
			if (commonWindowLow < currChannel.bump_window_min)
				commonWindowLow = currChannel.bump_window_min;
			if (commonWindowHigh < currChannel.bump_window_max)
				commonWindowHigh = currChannel.bump_window_max;
			// calculate overlap factors for all bump windows so far 
			for (unsigned int j = 0; j < bumpOverlaps.size(); j++)
			{
				Double_t low = (commonWindowLow<channel_list[j].bump_window_min) ? channel_list[j].bump_window_min : commonWindowLow;
				Double_t high = (commonWindowHigh>channel_list[j].bump_window_max) ? channel_list[j].bump_window_max : commonWindowHigh;
				bumpOverlaps[j] = (high - low) / (channel_list[j].bump_window_max - channel_list[j].bump_window_min);
				if (printOut) 
					Info("EvaluateMultiChannelTestStatistic","bump overlap (%d) = %f", j, bumpOverlaps[j]);
				if (bumpOverlaps[j] < fMultiChannelBumpOverlapFactor) 
					return 0; //one bump not overlapping enough, so just return 0;
			}
		}
		return out;
	}

	TH1* bumphunter::GenerateToyMC() 
	{
	   return GeneratePseudoData(channel_list[current_channel].hist_bkg);
	}

	TH1* bumphunter::GeneratePseudoData(TH1* bkg) 
	{
		if (!bkg) 
		{
			Error("GeneratePseudoData", "No background distribution given"); 
			return 0;
		}

		TH1* pseudodata = channel_list[current_channel].hist_pseudodata;

		pseudodata->Reset();

		//Double_t gausFactor = ROOT::Math::gaussian_cdf(pRand.Gaus(0,1),1); //where in the quantile spectrum to be for correlated uncertainties
		Double_t gausFactor = pRand.Uniform(); // choose quantile to use uniformly so that we fairly sample the possible gamma means

		Int_t startBin = bkg->FindFixBin(GetSearchLowEdge()); //first bin to use in search
		Int_t stopBin = bkg->FindFixBin(GetSearchHighEdge()) - 1; //last bin to use in search

		//loop over bins in the search region and fill with pseudo data
		for (Int_t j = startBin; j<= stopBin; j++)
		{
			Double_t mean = bkg->GetBinContent(j);
			Double_t err = bkg->GetBinError(j);

			if (mean == 0 && err == 0)
				continue; // leaves bin as 0 entries

			// define gamma parameters (a,b)
			Double_t b = mean / (err * err); // = E/V
			Double_t a = mean * b; // = E^2/V
			Double_t randomMean(0.); //used in binModel=0

			// which model are using 
			switch (fBinModel) 
			{
			case BUMP_POISSON_GAMMA: //poisson convoluted with a gamma distribution for the mean parameter
				//pick a random value for the mean
				randomMean = r.Gamma(a, 1./b);
				pseudodata->SetBinContent(j, pRand.PoissonD(randomMean));
				break;
			case BUMP_POISSON: //poisson with no uncertainty on the mean (i.e. the bin error is meaningless 
				pseudodata->SetBinContent(j, pRand.PoissonD(mean));
				break;
			case BUMP_GAUSSIAN: //gaussian with mean and sigma as given 
				pseudodata->SetBinContent(j, pRand.Gaus(mean, err));
				//alternative: pseudodata->SetBinContent(j, std::normal_distribution<double>(mean, err)(rd));
			case BUMP_POISSON_GAUSSIAN: //poisson convoluted with a gamma on the mean, but all the errors are correlated
				//use the gausfactor to decide how far from the 0.5 (mean) position to stray 
				randomMean = ROOT::Math::gamma_quantile(gausFactor, a, 1. / b);
				pseudodata->SetBinContent(j, pRand.PoissonD(randomMean));
				break;
			}
		}
		return pseudodata;
	}

/* NAMESPACE */
}
