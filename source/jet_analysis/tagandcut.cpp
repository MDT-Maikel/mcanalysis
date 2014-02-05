/* Jet analysis class: Tag and cut functions
 *
 * 
*/

#include "jet_analysis.h"


/* NAMESPACE */
namespace analysis 
{

	// ========================================== // 
	//	       	 Tagging functions	 		   	  //
	// ========================================== //

	fastjet::PseudoJet jet_analysis::JHTopTagging(const fastjet::PseudoJet & jet)
	{
		// Definition of JH Top-tagger algorithm
		fastjet::JHTopTagger top_tagger(JHTopTagger_delta_p, JHTopTagger_delta_r);
		top_tagger.set_top_selector(fastjet::SelectorMassRange(150,200));
		top_tagger.set_W_selector  (fastjet::SelectorMassRange( 65, 95));
		fastjet::PseudoJet tagged = top_tagger(jet);

		// Print Top-tagging details
		if (printTopTagDetails)
			cout << "\nRan the following top tagger:\n" << top_tagger.description() << endl;

		// Set TagInfo
		if ( tagged!=0 ) // return top-tagged jet
		{
			bool top_tag = true;
			bool w_tag = false;
			bool h_tag = false;
			tagged.set_user_info(new TagInfo(top_tag, w_tag, h_tag));
		}

		return tagged;
	}

	fastjet::PseudoJet jet_analysis::HEPTopTagging(const fastjet::PseudoJet & jet)
	{
		// Definition of HEP Top-tagger algorithm
		double topmass=172.3;
		double wmass=80.4;
		HEPTopTagger top_tagger(*(jet.associated_cs()),jet,topmass,wmass);
		top_tagger.set_top_range(150.,200.);
		top_tagger.run_tagger();		

		// Print Top-tagging details
		if (printTopTagDetails)
		{
			std::cout << "\nRan HEPTopTagger with the following parameters:" << std::endl; 
			top_tagger.get_setting();
			std::cout << std::endl;
			// top_tagger.get_info();
			// cout << endl;
		}

		// Set TagInfo
		if ( top_tagger.is_masscut_passed() ) // return top-tagged jet
		{
			fastjet::PseudoJet top = top_tagger.top_candidate();
			// fastjet::PseudoJet b = top_tagger.top_subjets().at(0);
			// fastjet::PseudoJet W1 = top_tagger.top_subjets().at(1);
			// fastjet::PseudoJet W2 = top_tagger.top_subjets().at(2);
			bool top_tag = true;
			bool w_tag = false;
			bool h_tag = false;
			top.set_user_info(new TagInfo(top_tag, w_tag, h_tag));

			return top;
		}
		else // return void jet
		{
			fastjet::PseudoJet voidJet;
			return voidJet;
		}

	}

	fastjet::PseudoJet jet_analysis::BDRSTagging(const fastjet::PseudoJet & jet)
	{
		/* from example 12 of FastJet */

		// Definition of Mass drop tagger with \mu=0.667 and ycut=0.09
		fastjet::MassDropTagger md_tagger(0.667, 0.09);
		fastjet::PseudoJet tagged = md_tagger(jet);

		// Print BDRS-tagging details
		if (printBDRSDetails)
			cout << "\nRan the following BDRS tagger:\n" << md_tagger.description() << endl << endl;

		if ( tagged!=0 )
		{
			// Filter the tagged jet, to remove UE & pileup contamination
			fastjet::PseudoJet parent1 = tagged.pieces()[0];
	  		fastjet::PseudoJet parent2 = tagged.pieces()[1];
	  		double   Rjj = parent1.delta_R(parent2);
	  		double   Rfilt = min(Rjj/2, 0.3); // somewhat arbitrary choice
	  		unsigned nfilt = 3;               // number of pieces we'll take

	  		fastjet::JetAlgorithm			cambridge = fastjet::cambridge_algorithm;
	  		fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;
			fastjet::Strategy               strategy = fastjet::Best;
			fastjet::Filter filter(
				fastjet::JetDefinition(cambridge, Rfilt, recombScheme, strategy),
				fastjet::SelectorNHardest(nfilt)
				);
	  		fastjet::PseudoJet filtered = filter(tagged);
	  		std::vector< fastjet::PseudoJet > filtered_pieces = filtered.pieces();

	  		// Set TagInfo
	  		fastjet::PseudoJet candidate;
	  		candidate = join(filtered_pieces[0], filtered_pieces[1]);
	  		double invMass = candidate.m();
			if ( BDRS_w_min < invMass && invMass < BDRS_w_max ) // return W-tagged jet
			{
				bool top_tag = false;
				bool w_tag = true;
				bool h_tag = false;
				filtered.set_user_info(new TagInfo(top_tag, w_tag, h_tag));
				return filtered;
			}
			else if ( BDRS_higgs_min < invMass && invMass < BDRS_higgs_max ) // return h-tagged jet
			{
				bool top_tag = false;
				bool w_tag = false;
				bool h_tag = true;
				filtered.set_user_info(new TagInfo(top_tag, w_tag, h_tag));
				return filtered;
			}
			else // return void jet
			{
				fastjet::PseudoJet voidJet;
				return voidJet;
			}
		}

		return tagged;
	}


	// ========================================== // 
	//	          Apply cut function 		   	  //
	// ========================================== //

	double jet_analysis::reduce_sample(cuts cut_list)
	{	
		// map_lhco_taggedJets iterator definition
		std::map< event *, std::vector <fastjet::PseudoJet> >::iterator it;

		// extract vector of event pointers "events"
		std::vector< event * > events;
		for (it=map_lhco_taggedJets.begin(); it!=map_lhco_taggedJets.end(); ++it)
			events.push_back(it->first);

		// perform the cuts defined in cut_list on "events"
		cut_list.apply(events);
		cut_list.write(cout);
		double eff = cut_list.efficiency();

		// retrieve from map_lhco_taggedJets only the elements which passed the cuts
		std::map< event *, std::vector <fastjet::PseudoJet> > reduced_map;
		for (unsigned int i = 0; i < events.size(); ++i)
		{
			it = map_lhco_taggedJets.find(events[i]);
			if ( it!=map_lhco_taggedJets.end() )
				reduced_map.insert( 
					std::make_pair( events[i], 
									map_lhco_taggedJets[events[i]] ) 
					); // end insert
		}

		// resize the map_lhco_taggedJets
		map_lhco_taggedJets = reduced_map;

		// return total efficiency
		return eff;
	}

	double jet_analysis::require_fatjet_pt(const double & ptcut, const int & n)
	{
		std::map< event *, std::vector <fastjet::PseudoJet> >::iterator it;
		std::map< event *, std::vector <fastjet::PseudoJet> > reduced_map;
		int map_size = 0, passed = 0;
		double eff;

		for (it=map_lhco_taggedJets.begin(); it!=map_lhco_taggedJets.end(); ++it)
		{
			map_size++;
			int jetcount = 0;

			if ((it->second).size()!=0)
			{
				for (unsigned int i = 0; i < (it->second).size(); ++i)
				{
					fastjet::PseudoJet fatjet;
					fatjet = (it->second)[i];

					if (fatjet.pt() > ptcut)
						jetcount++;

				} // for (it->second).size() -> calculates how many FatJets satisfy pT requirement within the event

				if (jetcount >= n)
				{
					passed++;
					reduced_map.insert( std::make_pair(it->first,it->second) );
				}

			} // if ((it->second).size()!=0) -> look among events with reconstructed FatJets

		} // for (map_lhco_taggedJets)

		// calculate efficiency
		eff = (double) passed/map_size;

		// resize the map_lhco_taggedJets
		map_lhco_taggedJets = reduced_map;

		return eff;
	}

	double jet_analysis::require_top_tagged(const int & n)
	{
		std::map< event *, std::vector <fastjet::PseudoJet> >::iterator it;
		std::map< event *, std::vector <fastjet::PseudoJet> > reduced_map;
		int map_size = 0, passed = 0;
		double eff;

		for (it=map_lhco_taggedJets.begin(); it!=map_lhco_taggedJets.end(); ++it)
		{
			map_size++;
			int toptag = 0;

			if ((it->second).size()!=0)
			{
				for (unsigned int i = 0; i < (it->second).size(); ++i)
				{
					fastjet::PseudoJet tagged;
					tagged = (it->second)[i];

					if (tagged.user_info<TagInfo>().top_tag())
					{
						// cout << "Found tagged top" << endl;
						// cout << "top candidate:     " << tagged << endl;
						// cout << "|_ W   candidate:  " << tagged.structure_of<fastjet::JHTopTagger>().W() << endl;
						// cout << "|  |_  W subjet 1: " << tagged.structure_of<fastjet::JHTopTagger>().W1() << endl;
						// cout << "|  |_  W subjet 2: " << tagged.structure_of<fastjet::JHTopTagger>().W2() << endl;
						// cout << "|  cos(theta_W) =  " << tagged.structure_of<fastjet::JHTopTagger>().cos_theta_W() << endl;
						// cout << "|_ non-W subjet:   " << tagged.structure_of<fastjet::JHTopTagger>().non_W() << endl;
						toptag++;
					}

				} // for (it->second).size() -> calculates how many tops are tagged within the event

				if (toptag >= n)
				{
					passed++;
					reduced_map.insert( std::make_pair(it->first,it->second) );
				}

			} // if ((it->second).size()!=0) -> look among events with tagged jets

		} // for (map_lhco_taggedJets)

		// calculate efficiency
		eff = (double) passed/map_size;

		// resize the map_lhco_taggedJets
		map_lhco_taggedJets = reduced_map;

		return eff;
	}

	double jet_analysis::require_higgs_tagged(const int & n)
	{
		std::map< event *, std::vector <fastjet::PseudoJet> >::iterator it;
		std::map< event *, std::vector <fastjet::PseudoJet> > reduced_map;
		int map_size = 0, passed = 0;
		double eff;

		for (it=map_lhco_taggedJets.begin(); it!=map_lhco_taggedJets.end(); ++it)
		{
			map_size++;
			int higgstag = 0;

			if ((it->second).size()!=0)
			{
				for (unsigned int i = 0; i < (it->second).size(); ++i)
				{
					fastjet::PseudoJet tagged;
					tagged = (it->second)[i];

					if (tagged.user_info<TagInfo>().h_tag())
					{
						// cout << "Found tagged higgs" << endl;
						// std::vector< fastjet::PseudoJet > tagged_pieces = tagged.pieces();
						// fastjet::PseudoJet candidate;
						// candidate = join(tagged_pieces[0], tagged_pieces[1]);
						// double invMass = candidate.m();

						// cout << "Filtered pieces are " << endl;
						// for (unsigned i = 0; i < 3 && i < tagged_pieces.size(); i++) 
						// {
						// cout << " " << tagged_pieces[i] << endl;
						// }
						// cout << "Filtered total is " << endl;
						// cout << " " << tagged << endl;
						// cout << "Invariant mass of the two leading pieces: " << invMass << endl;
						higgstag++;
					}

				} // for (it->second).size() -> calculates how many Higgses are tagged within the event

				if (higgstag >= n)
				{
					passed++;
					reduced_map.insert( std::make_pair(it->first,it->second) );
				}

			} // if ((it->second).size()!=0) -> look among events with tagged jets

		} // for (map_lhco_taggedJets)

		// calculate efficiency
		eff = (double) passed/map_size;
		
		// resize the map_lhco_taggedJets
		map_lhco_taggedJets = reduced_map;

		return eff;
	}

	double jet_analysis::require_w_tagged(const int & n)
	{
		std::map< event *, std::vector <fastjet::PseudoJet> >::iterator it;
		std::map< event *, std::vector <fastjet::PseudoJet> > reduced_map;
		int map_size = 0, passed = 0;
		double eff;

		for (it=map_lhco_taggedJets.begin(); it!=map_lhco_taggedJets.end(); ++it)
		{
			map_size++;
			int wtag = 0;

			if ((it->second).size()!=0)
			{
				for (unsigned int i = 0; i < (it->second).size(); ++i)
				{
					fastjet::PseudoJet tagged;
					tagged = (it->second)[i];

					if (tagged.user_info<TagInfo>().w_tag())
					{
						// cout << "Found tagged W" << endl;
						// std::vector< fastjet::PseudoJet > tagged_pieces = tagged.pieces();
						// fastjet::PseudoJet candidate;
						// candidate = join(tagged_pieces[0], tagged_pieces[1]);
						// double invMass = candidate.m();

						// cout << "Filtered pieces are " << endl;
						// for (unsigned i = 0; i < 3 && i < tagged_pieces.size(); i++) 
						// {
						// cout << " " << tagged_pieces[i] << endl;
						// }
						// cout << "Filtered total is " << endl;
						// cout << " " << tagged << endl;
						// cout << "Invariant mass of the two leading pieces: " << invMass << endl;
						wtag++;
					}

				} // for (it->second).size() -> calculates how many W are tagged within the event

				if (wtag >= n)
				{
					passed++;
					reduced_map.insert( std::make_pair(it->first,it->second) );
				}

			} // if ((it->second).size()!=0) -> look among events with tagged jets

		} // for (map_lhco_taggedJets)

		// calculate efficiency
		eff = (double) passed/map_size;
		
		// resize the map_lhco_taggedJets
		map_lhco_taggedJets = reduced_map;

		return eff;
	}

	double jet_analysis::require_t_or_w_tagged(const int & n)
	{
		std::map< event *, std::vector <fastjet::PseudoJet> >::iterator it;
		std::map< event *, std::vector <fastjet::PseudoJet> > reduced_map;
		int map_size = 0, passed = 0;
		double eff;

		for (it=map_lhco_taggedJets.begin(); it!=map_lhco_taggedJets.end(); ++it)
		{
			map_size++;
			int twtag = 0;

			if ((it->second).size()!=0)
			{
				for (unsigned int i = 0; i < (it->second).size(); ++i)
				{
					fastjet::PseudoJet tagged;
					tagged = (it->second)[i];

					if ( tagged.user_info<TagInfo>().top_tag() || tagged.user_info<TagInfo>().w_tag() )
						twtag++;

				} // for (it->second).size() -> calculates how many W are tagged within the event

				if (twtag >= n)
				{
					passed++;
					reduced_map.insert( std::make_pair(it->first,it->second) );
				}

			} // if ((it->second).size()!=0) -> look among events with tagged jets

		} // for (map_lhco_taggedJets)

		// calculate efficiency
		eff = (double) passed/map_size;
		
		// resize the map_lhco_taggedJets
		map_lhco_taggedJets = reduced_map;

		return eff;
	}


/* NAMESPACE */
}