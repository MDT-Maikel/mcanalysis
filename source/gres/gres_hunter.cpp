/* GRESHunter class
 *
 * 
*/

#include "gres_hunter.h" 


/* NAMESPACE */
namespace analysis
{

	/* properties */
	
	void greshunter::add_background(const std::vector<event*> &events, double weight)
	{
		events_bkg.push_back(events);
		weight_bkg.push_back(weight);		
	}
	
	void greshunter::set_signal(const std::vector<event*> &events, double weight)
	{
		events_sig = events;
		weight_sig = weight;		
	}
	
	void greshunter::set_topology(const std::vector<int> & top)
	{
		topology = top;
	}

	void greshunter::set_top_cuts(const std::vector<double> & cuts)
	{
		top_cuts = cuts;
	}
	
	void greshunter::set_folder(std::string f)
	{
		gres_folder = f;
	}
	
	/* bump hunting */
	
	void greshunter::run()
	{
		unsigned int nr_jets = 0;
		for (unsigned int i = 0; i < topology.size(); i++)
			nr_jets += topology[i];
		
		// only look at resonances of different sizes
		std::vector<int> red_topology;
		std::vector<double> red_cuts;
		for (unsigned int i = 0; i < topology.size(); i++)
		{
			int resonance = topology[i];
			bool add = (resonance > 1);
			for (unsigned int j = 0; j < red_topology.size(); j++)
				if (red_topology[j] == resonance)
					add = false;
			if (add)
			{
				red_topology.push_back(resonance);
				red_cuts.push_back(top_cuts[i]);	
			}
		}
		
		std::cout << "GRES Hunter for topology: {";
		for (unsigned int i = 0; i < red_topology.size(); i++)
		{	
			std::cout << red_topology[i];
			if (i != red_topology.size() - 1)
				std::cout << ", ";
		}
		std::cout << "} resonance structure." << std::endl;
		
		// produce all the needed invariant mass spectra, these consist of all combinations nr_jets over resonance size
		for (unsigned int i = 0; i < red_topology.size(); i++)
		{
			int resonance = red_topology[i];
			unsigned int nr_combs = static_cast<unsigned int>(boost::math::binomial_coefficient<double>(nr_jets, resonance)); // TODO: Is this failsafe?
			std::vector<std::vector<int> > combinations;
			std::vector<int> start_comb;
			for (unsigned int j = resonance; j > 0; j--)
				start_comb.push_back(j);
			// produce all combinations which lead to the resonance
			for (unsigned int j = 0; j < nr_combs; j++)
			{
				combinations.push_back(start_comb);
				start_comb[0]++;
				for (unsigned int k = 0; k < start_comb.size() - 1; k++)
				{
					if (static_cast<unsigned int>(start_comb[k]) > nr_jets - k)
					{
						start_comb[k + 1]++;
						for (int l = k; l >= 0; l--)
							start_comb[l] = start_comb[l + 1] + 1;
					}
				}
			}

			// make the invariant mass spectrum plot for each of the combinations
			for (unsigned int j = 0; j < combinations.size(); j++)
			{
				std::vector<int> comb = combinations[j];
				run(comb, red_cuts[i]);				
			}
		}		
	}
	
	void greshunter::run(const std::vector<int> & comb, double mass_cut)
	{
		// turn comb into string and print to screen
		std::cout << "RUN BUMPHUNTER FOR: {";
		for (unsigned int i = 0; i < comb.size(); i++)
		{	
			std::cout << comb[i];
			if (i != comb.size() - 1)
				std::cout << ", ";
		}
		std::cout << "}." << std::endl;		
		
		// convert background events into mass data	
		std::vector< std::vector<double> > mass_bkg;
		for (unsigned int i = 0; i < events_bkg.size(); i++)
		{
			std::vector<event*> events = events_bkg[i];
			std::vector<double> mass;
			for (unsigned int j = 0; j < events.size(); j++)
			{
				event *ev = events[j];
				mass.push_back(ev->mass(ptype_jet, comb));
			}
			mass_bkg.push_back(mass);	
		}

		// convert signal events into mass data
		std::vector<double> mass_sig;
		for (unsigned int i = 0; i < events_sig.size(); i++)
		{
			event *ev = events_sig[i];
			mass_sig.push_back(ev->mass(ptype_jet, comb));
		}
		
		
		// create and fill the background and signal histograms and normalize them
		std::vector<TH1F*> hist_bkgs;
		for (unsigned int i = 0; i < mass_bkg.size(); i++)
		{
			TH1F* hist;
			hist = new TH1F("", "bkg", 50, 0, 3000);
			std::vector<double> mass = mass_bkg[i];
			for (unsigned int j = 0; j < mass.size(); j++)
				hist->Fill(mass[j], weight_bkg[i]);
			hist->Scale(weight_bkg[i] / hist->Integral());	
			hist_bkgs.push_back(hist);
		}
		
		TH1F* hist_sig;
		hist_sig = new TH1F("", "sig", 50, 0, 3000);
		for (unsigned int i = 0; i < mass_sig.size(); i++)
			hist_sig->Fill(mass_sig[i], weight_sig);
		hist_sig->Scale(weight_sig / hist_sig->Integral());
		
		// make prediction and data histograms
		TH1F* hist_pred = new TH1F("", "pred", 50, 0, 3000);
		TH1F* hist_data = new TH1F("", "data", 50, 0, 3000);
		hist_data->Add(hist_sig);
		for (unsigned int i = 0; i < hist_bkgs.size(); i++)
		{
			hist_pred->Add(hist_bkgs[i]);
			hist_data->Add(hist_bkgs[i]);
		}
		
		// ERROR: plot histograms
		//std::random_device rd;
		//TCanvas* canvas = new TCanvas(boost::lexical_cast<std::string>(rd()).c_str(), "");
		//hist_pred->Draw("");
		//canvas->Print("TEST_PRED.png");
		//TCanvas* canvas2 = new TCanvas(boost::lexical_cast<std::string>(rd()).c_str(), "");
		//hist_data->Draw("");
		//canvas2->Print("TEST_DATA.png");
		
		// set up bumphunter analysis
		bumphunter hunt(hist_pred, hist_data);
		hunt.set_folder(gres_folder);
		std::string comb_string = "";
		for (unsigned int i = 0; i < comb.size(); i++)
			comb_string += boost::lexical_cast<std::string>(comb[i]);
		hunt.set_name("invmass_" + comb_string);
		hunt.SetNPseudoExperiments(1000);
		hunt.SetBinModel(1);
		hunt.SetSearchRegion(mass_cut, 3000);
	
		// run bumphunter analysis
		hunt.run();
	}


/* NAMESPACE */
}
