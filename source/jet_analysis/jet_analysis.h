/* Jet analysis class
 *
 * 
*/

#ifndef INC_JET_ANALYSIS
#define INC_JET_ANALYSIS

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string> 
#include <vector> 

#include "boost/filesystem.hpp"
#include "boost/lexical_cast.hpp"

#include "Pythia8/Pythia.h"
#include "Pythia8/FastJet3.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JHTopTagger.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Filter.hh"

#include "../particle/lhco.h"
#include "../event/event.h"
#include "../utility/utility.h"
#include "../cuts/cuts.h"
#include "HEPTopTagger.h"


/* NAMESPACE */
namespace analysis
{

	/* jet analysis */
	class jet_analysis
	{

	public:
		/* class con- & destructor: jet_analysis.cpp */
		jet_analysis();
		~jet_analysis();

		/* set parameters and input files: settings.cpp */
		void import_lhe(const std::string & name);
		void set_fast_showering();
		void import_lhco(const std::string & name);
		void set_nEvents(const int & events);
		void set_Isolation(const std::string & type, const double & eta, const double & pt, const double & Rcone, const double & EtCone = 1.8);
		void set_Rsize_fat(const double & R);
		void set_Rsize_skinny(const double & R);
		void set_skinnyjet_algorithm(const std::string & type);
		void set_JHTopTagger(const double & delta_p, const double & delta_r);
		void set_BDRS_w_range(const double & w_min, const double & w_max);
		void set_BDRS_higgs_range(const double & higgs_min, const double & higgs_max);
		void undo_TopTagging();
		void undo_BDRSTagging();

		/* merging settings: settings.cpp */
		void set_merging_process(const std::string & process);
		void set_merging_njmax(const int & njet);
		void set_merging_njadditional(const int & njet);
		void set_merging_scale(const double & scale);
		bool MergingSettings();
		void AllowPythiaDecay();

		/* isolation functions: isolation.cpp */
		bool isolatedElectron(const int & j, const Pythia8::Event & particles);
		bool isolatedMuon(const int & j, const Pythia8::Event & particles);
		bool isolatedPhoton(const int & j, const Pythia8::Event & particles);
		bool JetElectronOverlapping(const fastjet::PseudoJet & jet, const std::vector< fastjet::PseudoJet > & leptons);

		/* tag and cut functions: tagandcut.cpp */
		fastjet::PseudoJet JHTopTagging(const fastjet::PseudoJet & jet);
		fastjet::PseudoJet HEPTopTagging(const fastjet::PseudoJet & jet);
		fastjet::PseudoJet BDRSTagging(const fastjet::PseudoJet & jet);
		double reduce_sample(cuts cut_list);
		double require_fatjet_pt(const double & ptcut, const int & n = 1);
		double require_top_tagged(const int & n);
		double require_higgs_tagged(const int & n);
		double require_w_tagged(const int & n);
		double require_t_or_w_tagged(const int & n);
		
		/* extract lhco or fatjets from map: tagandcut.cpp */
		std::vector< event* > events();
		std::vector< std::vector< fastjet::PseudoJet > > fatjets();

		/* initialisation function: initialisation.cpp */
		void initialise(fastjet::PseudoJet (jet_analysis::*TopTagger)(const fastjet::PseudoJet &) = &jet_analysis::JHTopTagging);
	
	public: // TODO: private
		
		/* basic settings, input and flags */
		std::string lhe_input;
		std::string lhco_input;
		int nEvent;
		bool firstEvent;
		bool printTopTagDetails;
		bool printBDRSDetails;
		bool fast_showering;
		bool importedLHE;
		bool importedLHCO;

		/* detector range and isolation parameters */
		double MaxEta;
		double jetMinPt;
		double electronMaxEta;
		double electronMinPt;
		double muonMaxEta;
		double muonMinPt;
		double deltaR_IsolatedLepton; 
		double sumEtInCone_IsolatedMuon;
		double photonMaxEta;
		double photonMinPt;
		double deltaR_IsolatedPhoton; 
		double sumEtInCone_IsolatedPhoton; 

		/* merging flags and parameters */
		bool DoMerging;
		bool Process;
		bool NJetMax;
		bool Scale;
		bool PythiaDecay;
		std::string MergingProcess;
		int MergingNJetMax;
		double MergingScale;

		/* jet clustering parameters */
		fastjet::JetAlgorithm algorithm_fat;
		double Rsize_fat;
		fastjet::JetAlgorithm algorithm_skinny;
		double Rsize_skinny;
		
		/* JHTopTagger parameters */
		bool DoTopTagger;
		double JHTopTagger_delta_p; 
		double JHTopTagger_delta_r;

		/* BDRS parameters */
		bool DoBDRS;
		double BDRS_w_min;
		double BDRS_w_max;
		double BDRS_higgs_min;
		double BDRS_higgs_max;

		/* map (lhco, taggedJets) */
		std::map< event *, std::vector< fastjet::PseudoJet > > map_lhco_taggedJets;

	};

	/* tagging information */
	class TagInfo : public fastjet::PseudoJet::UserInfoBase
	{

	public:
		/* default con- & destructor */
		TagInfo(
			const bool & top_tag = false, 
			const bool & w_tag = false, 
			const bool & h_tag = false
			);
		~TagInfo();

		/* access to tagging information */
		bool top_tag() const;
		bool w_tag() const;
		bool h_tag() const;
		
			
	protected:
		bool _top_tag; // top tagging information
		bool _w_tag;   // w tagging information
		bool _h_tag;   // higgs tagging information

	};

	/* b-flavour information */
	class FlavourInfo : public fastjet::PseudoJet::UserInfoBase
	{

	public:
		/* default con- & destructor */
		FlavourInfo( const bool & b = false );
		~FlavourInfo();

		/* access to flavour information */
		bool b_type() const;
		
	protected:
		bool _b; // b-quark jet

	};

	/* b-hadron pdg identification */
	bool isBHadron(int pdg);

	/* jet << operator */
	ostream & operator<<(ostream &, const fastjet::PseudoJet &);

/* NAMESPACE */
}

#endif
