/* jet_analysis class
 *
 * 
*/

#ifndef INC_JET_ANALYSIS
#define INC_JET_ANALYSIS

#include <cmath>
#include <fstream>
#include <iostream>
#include <string> 
#include <vector> 
#include <map>

#include "Pythia8/Pythia.h"
#include "Pythia8/FastJet3.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
// #include "fastjet/CompositeJetStructure.hh.hh"
#include "fastjet/tools/JHTopTagger.hh"
#include <fastjet/tools/MassDropTagger.hh>
#include <fastjet/tools/Filter.hh>

#include "boost/filesystem.hpp"
#include "../particle/lhco.h"
#include "../event/event.h"
#include "../utility/utility.h"
#include "../cuts/cuts.h"
#include "../cuts/cuts_standard.h"
#include "../cuts/cuts_advanced.h"


using namespace Pythia8;

namespace analysis {

	class jet_analysis
	{

	public: //TODO: private?

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

		// lhe input & number of events to be analysed by Phytia
		std::string lhe_name;
		int nEvent;
		bool firstEvent;
		bool printTopTagDetails;
		bool printBRDSDetails;

		// Jet clustering parameters
		double Rsize_fat;
		double Rsize_skinny;

		// JHTopTagger parameters
		bool DoTopTagger;
		double JHTopTagger_delta_p; 
		double JHTopTagger_delta_r;

		// map (lhco with skinny jets, taggedJets)
		std::map< event *, std::vector <fastjet::PseudoJet> > map_lhco_fatJets;

	public:

		//=== Class con- & destructor ===//
		jet_analysis();
		~jet_analysis();

		//=== Set parameters and lhe input ===//
		void set_Isolation(const std::string & type, const double & eta, const double & pt, const double & Rcone, const double & EtCone = 1.8);
		void set_nEvent(const int & events);
		void set_Rsize_fat(const double & R);
		void set_Rsize_skinny(const double & R);
		void set_JHTopTagger(const double & delta_p, const double & delta_r);
		void add_lhe(const std::string & name);

		//=== Isolation functions ===//
		bool isolatedElectron(const int & j, const Event & particles);
		bool isolatedMuon(const int & j, const Event & particles);
		bool isolatedPhoton(const int & j, const Event & particles);
		bool JetElectronOverlapping(const fastjet::PseudoJet & jet, const std::vector< fastjet::PseudoJet > & leptons);

		//=== Tagging functions ===//
		fastjet::PseudoJet JHTopTagging(const fastjet::PseudoJet & jet);
		fastjet::PseudoJet BRDSTagging(const fastjet::PseudoJet & jet);

		//=== Apply cut function ===//
		void reduce_sample(cuts cut_list);
		
		//=== Analysis core function ===//
		void initialise(fastjet::PseudoJet (jet_analysis::*TopTagger)(const fastjet::PseudoJet &));

		// template <typename cut_type>
		// void reduce_sample(cut_type cut_def);

	};


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


	//=== Jet << operator ===//
	ostream & operator<<(ostream &, const fastjet::PseudoJet &);

/* NAMESPACE */
}

/* 	TEMPLATE IMPLEMENTATIONS */
// #include "jet_analysis.tpp"

#endif
