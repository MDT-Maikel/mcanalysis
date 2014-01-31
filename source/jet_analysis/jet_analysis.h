/* Jet analysis class
 *
 * 
*/

#ifndef INC_JET_ANALYSIS
#define INC_JET_ANALYSIS

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

	class jet_analysis
	{
	// ===== TODO: private? ===== //
	public:
		// lhe input files, number of events to be analysed by Pythia and print options
		std::vector< std::string > lhe_input;
		int nEvent;
		bool firstEvent;
		bool printTopTagDetails;
		bool printBDRSDetails;

		// Detector range and Isolation parameters
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

		// Merging procedure
		bool DoMerging;
		double MergingTMS;
		int MergingNJetMax;
		std::string MergingProcess;

		// Jet clustering parameters
		double Rsize_fat;
		double Rsize_skinny;

		// JHTopTagger parameters
		bool DoTopTagger;
		double JHTopTagger_delta_p; 
		double JHTopTagger_delta_r;

		// BDRS parameters
		bool DoBDRS;
		double BDRS_w_min;
		double BDRS_w_max;
		double BDRS_higgs_min;
		double BDRS_higgs_max;

		// map (lhco with skinny jets, taggedJets)
		std::map< event *, std::vector <fastjet::PseudoJet> > map_lhco_taggedJets;

	public:
		//=== Class con- & destructor ===//
		jet_analysis();
		~jet_analysis();

		//=== Set parameters and lhe input ===//
		void add_lhe(const std::string & name);
		void set_nEvents(const int & events);
		void set_Isolation(const std::string & type, const double & eta, const double & pt, const double & Rcone, const double & EtCone = 1.8);
		void set_Rsize_fat(const double & R);
		void set_Rsize_skinny(const double & R);
		void set_JHTopTagger(const double & delta_p, const double & delta_r);
		void set_BDRS_w_range(const double & w_min, const double & w_max);
		void set_BDRS_higgs_range(const double & higgs_min, const double & higgs_max);
		void undo_TopTagging();
		void undo_BDRSTagging();

		//=== Merging procedure ===//
		void set_merging(const double & ms, const int & njmax, const std::string & process);

		//=== Isolation functions ===//
		bool isolatedElectron(const int & j, const Pythia8::Event & particles);
		bool isolatedMuon(const int & j, const Pythia8::Event & particles);
		bool isolatedPhoton(const int & j, const Pythia8::Event & particles);
		bool JetElectronOverlapping(const fastjet::PseudoJet & jet, const std::vector< fastjet::PseudoJet > & leptons);

		//=== Tagging functions ===//
		fastjet::PseudoJet JHTopTagging(const fastjet::PseudoJet & jet);
		fastjet::PseudoJet HEPTopTagging(const fastjet::PseudoJet & jet);
		fastjet::PseudoJet BDRSTagging(const fastjet::PseudoJet & jet);

		//=== Analysis core function ===//
		void initialise(fastjet::PseudoJet (jet_analysis::*TopTagger)(const fastjet::PseudoJet &) = &jet_analysis::JHTopTagging);
		std::vector< event* > events();

		//=== Apply cut functions ===//
		void reduce_sample(cuts cut_list);
		double require_top_tagged(const int & n);
		double require_higgs_tagged(const int & n);
		double require_w_tagged(const int & n);
		double require_t_or_w_tagged(const int & n);

	};


	//=== Tagging information ===//
	class TagInfo : public fastjet::PseudoJet::UserInfoBase
	{
	protected:
		bool _top_tag;	// top tagging information
		bool _w_tag;	// w tagging information
		bool _h_tag;	// higgs tagging information

	public:
		// default con- & destructor
		TagInfo(
			const bool & top_tag = false, 
			const bool & w_tag = false, 
			const bool & h_tag = false
			);
		~TagInfo();

		// access to tagging information
		bool top_tag() const;
		bool w_tag() const;
		bool h_tag() const;

	};


	//=== b-flavour information ===//
	class FlavourInfo : public fastjet::PseudoJet::UserInfoBase
	{
	protected:
		bool _b;	// b-quark/jet

	public:
		// default con- & destructor
		FlavourInfo( const bool & b = false );
		~FlavourInfo();

		// access to flavour information
		bool b_type() const;

	};


	//=== b-hadron pdg identification ===//
	bool isBHadron(int pdg);

	//=== Jet << operator ===//
	ostream & operator<<(ostream &, const fastjet::PseudoJet &);

/* NAMESPACE */
}

#endif
