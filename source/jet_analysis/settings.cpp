/* Jet analysis class: Settings and input
 *
 * 
*/

#include "jet_analysis.h"


/* NAMESPACE */
namespace analysis 
{

	// ========================================== // 
	//	      Set parameters and input files  	  //
	// ========================================== //

	void jet_analysis::add_lhe(const std::string & name)
	{
		lhe_input.push_back(name);
	}

	void jet_analysis::import_lhco(const std::string & name)
	{
		importedLHCO = true;
		lhco_input = name;
	}

	void jet_analysis::set_nEvents(const int & events)
	{
		nEvent = events;
	}

	// ===== TODO: update using p_type ===== //
	void jet_analysis::set_Isolation(const std::string & type, const double & eta, const double & pt, const double & Rcone, const double & EtCone)
	{
		if (type == "electron")
		{
			electronMaxEta = eta;
			electronMinPt = pt;
			deltaR_IsolatedLepton = Rcone;
		}
		else if (type == "muon")
		{
			muonMaxEta = eta;
			muonMinPt = pt;
			deltaR_IsolatedLepton = Rcone;
			sumEtInCone_IsolatedMuon = EtCone;
		}
		else if (type == "photon")
		{
			photonMaxEta = eta;
			photonMinPt = pt;
			deltaR_IsolatedPhoton = Rcone;
			sumEtInCone_IsolatedPhoton = EtCone;
		}
		else
			std::cout << "Wrong type assignment." << std::endl;
	}

	void jet_analysis::set_Rsize_fat(const double & R)
	{
		Rsize_fat = R;
	}

	void jet_analysis::set_Rsize_skinny(const double & R)
	{
		Rsize_skinny = R;
	}

	void jet_analysis::set_JHTopTagger(const double & delta_p, const double & delta_r)
	{
		DoTopTagger = true;
		JHTopTagger_delta_p = delta_p;
		JHTopTagger_delta_r = delta_r;
	}

	void jet_analysis::set_BDRS_w_range(const double & w_min, const double & w_max)
	{
		DoBDRS = true;
		BDRS_w_min = w_min;
		BDRS_w_max = w_max;
	}

	void jet_analysis::set_BDRS_higgs_range(const double & higgs_min, const double & higgs_max)
	{
		DoBDRS = true;
		BDRS_higgs_min = higgs_min;
		BDRS_higgs_max = higgs_max;
	}

	void jet_analysis::undo_TopTagging()
	{
		DoTopTagger = false;
	}

	void jet_analysis::undo_BDRSTagging()
	{
		DoBDRS = false;
	}


	// ========================================== // 
	//	            Merging settings	 	   	  //
	// ========================================== //

	void jet_analysis::set_merging(const double & ms, const int & njmax, const std::string & process)
	{
		DoMerging = true;
		MergingTMS = ms;
		MergingNJetMax = njmax;
		MergingProcess = process;
	}


	// ========================================== // 
	//	         Extract lhco from map 		   	  //
	// ========================================== //
	std::vector< event * > jet_analysis::events()
	{
		// map_lhco_taggedJets iterator definition
		std::map< event *, std::vector <fastjet::PseudoJet> >::iterator it;

		// extract vector of event pointers "events"
		std::vector< event * > events;
		for (it=map_lhco_taggedJets.begin(); it!=map_lhco_taggedJets.end(); ++it)
			events.push_back(it->first);
			
		return events;
	}


/* NAMESPACE */
}