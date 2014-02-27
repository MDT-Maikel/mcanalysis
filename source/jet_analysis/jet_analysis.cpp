/* Jet analysis class: Con- & destructor and external functions/classes
 *
 * 
*/

#include "jet_analysis.h"


/* NAMESPACE */
namespace analysis 
{

	// ========================================== // 
	//	       	 Class con- & destructor 		  //
	// ========================================== //

	jet_analysis::jet_analysis() 
	{
		// Basic settings and flags
		nEvent = 1000;
		firstEvent = true;
		printTopTagDetails = true;
		printBDRSDetails = true;
		fast_showering = false;
		importedLHE = false;
		importedLHCO = false;

		// Detector range and Isolation parameters
		MaxEta = 4.9;
		jetMinPt = 20.;
		electronMaxEta = 2.47;
		electronMinPt = 20.;
		muonMaxEta = 2.4;
		muonMinPt = 10.;
		deltaR_IsolatedLepton = 0.2; 
		sumEtInCone_IsolatedMuon = 1.8; 
		photonMaxEta = 2.37;
		photonMinPt = 20.;
		deltaR_IsolatedPhoton = 0.2; 
		sumEtInCone_IsolatedPhoton = 2.6;

		// Merging procedure flags
		DoMerging = false;
		Process = false;
		NJetMax = false;
		Scale = false;

		// Jet clustering parameters
		algorithm_fat = fastjet::cambridge_algorithm;
		Rsize_fat = 1.4;
		algorithm_skinny = fastjet::antikt_algorithm;
		Rsize_skinny = 0.4;

		// JHTopTagger parameters
		DoTopTagger = true;
		JHTopTagger_delta_p = 0.05; 
		JHTopTagger_delta_r = 0.19;

		//BDRS parameters
		DoBDRS = true;
		BDRS_w_min = 65.;
		BDRS_w_max = 95.;
		BDRS_higgs_min = 100.;
		BDRS_higgs_max = 130.;
	}

	jet_analysis::~jet_analysis()
	{
	}


	// ========================================== // 
	//	          Tagging information 		   	  //
	// ========================================== //
	// default con- & destructor
	TagInfo::TagInfo(const bool & top_tag, const bool & w_tag, const bool & h_tag)
	{
		_top_tag = top_tag;
		_w_tag = w_tag;
		_h_tag = h_tag;
	}

	TagInfo::~TagInfo()
	{
	}

	// access to tagging information
	bool TagInfo::top_tag() const
		{ return _top_tag; }

	bool TagInfo::w_tag() const
		{ return _w_tag; }

	bool TagInfo::h_tag() const
		{ return _h_tag; }
		

	// ========================================== // 
	//	         b-flavour information 		   	  //
	// ========================================== //
	// default con- & destructor
	FlavourInfo::FlavourInfo(const bool & b)
	{
		_b = b;
	}

	FlavourInfo::~FlavourInfo()
	{
	}

	// access to flavour information
	bool FlavourInfo::b_type() const
		{ return _b; }


	// ========================================== // 
	//	      b-hadron pdg identification 	   	  //
	// ========================================== //
	bool isBHadron(int pdg)
	{
		while (pdg!=0)
		{
			int tmp = pdg%10;
			if (tmp == 5)
				return true;
			else
				pdg = int((pdg-tmp)/10);
		}

		return false;
	}

	// ========================================== // 
	//	            Jet << operator 		   	  //
	// ========================================== //
	ostream & operator<<(ostream & ostr, const fastjet::PseudoJet & jet) 
	{
	  ostr << "pt, y, phi =" << setprecision(2)
	       << " " << setw(9) << jet.perp() 
	       << " " << setw(9)  <<  jet.rap()  
	       << " " << setw(9)  <<  jet.phi()
	       << ", mass = " << setw(9) << jet.m();
	  return ostr;
	}

/* NAMESPACE */
}
