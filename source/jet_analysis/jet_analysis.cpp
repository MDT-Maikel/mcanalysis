/* jet_analysis class
 *
 * 
*/

#include "jet_analysis.h"

namespace analysis {

	// ========================================== // 
	//	       	 Class con- & destructor 		  //
	// ========================================== //

	jet_analysis::jet_analysis() 
	{
		// lhe input files, number of events to be analysed by Pythia and print options
		nEvent = 1000;
		firstEvent = true;
		printTopTagDetails = true;
		printBDRSDetails = true;

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

		// Merging procedure
		DoMerging = false;

		// Jet clustering parameters
		Rsize_fat = 1.4;
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
	//	       Set parameters and lhe input   	  //
	// ========================================== //

	void jet_analysis::add_lhe(const std::string & name)
	{
		lhe_input.push_back(name);
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
	//	            Merging procedure	 	   	  //
	// ========================================== //

	void jet_analysis::set_merging(const double & ms, const int & njmax, const std::string & process)
	{
		DoMerging = true;
		MergingTMS = ms;
		MergingNJetMax = njmax;
		MergingProcess = process;
	}


	// ========================================== // 
	//	            Isolation functions 	   	  //
	// ========================================== //

	bool jet_analysis::isolatedElectron(const int & j, const Event & particles) 
	{
		long id = particles[j].idAbs();

		// Enter the loop if particles[j] is an electron. Otherwise return false
		if ( id == 11 ) 
		{
			fastjet::PseudoJet electron = particles[j];

			// Check if electron is within range. Otherwise return false
			if ( abs(electron.pt()) > electronMinPt && 
			     abs(electron.eta()) < electronMaxEta ) 
			{
				// Check all other particles in event: add up their transverse
				// energies in a cone of radius deltaR_IsolatedLepton
				double sumEtInCone = 0.;
				for (int i=0; i<particles.size(); i++) 
				{
					// Only consider final state particles of the event record
					if (particles[i].isFinal())
					{
						fastjet::PseudoJet x = particles[i];

						// Only count visible particles that are not our electron in question
						if (i != j && 
						  particles[i].idAbs() != 12 && 
						  particles[i].idAbs() != 14 && 
						  particles[i].idAbs() != 16 && 
						  particles[i].idAbs() != 1000022 &&
						  particles[i].idAbs() != 8880022 &&
						  abs(particles[i].eta()) < MaxEta) 
						{
							double deltaEta = x.eta() - electron.eta();
							double deltaPhi = x.phi() - electron.phi();
							double deltaR   = sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
							if (deltaR < deltaR_IsolatedLepton)
								sumEtInCone += abs(x.pt());

						} // end if(visible and not considered electron)

					} // end if(final state)

				} // end for(particles.size())

				// isolation criterium: pT within cone less than 10% of pT(e)
				if (sumEtInCone < 0.1*abs(electron.pt())) 
					return true;

			} // end if(pTmin, etamax)

		} // end if(electron)

		return false;
	}

	bool jet_analysis::isolatedMuon(const int & j, const Event & particles) 
	{
		long id = particles[j].idAbs();

		// Enter the loop if particles[j] is a muon. Otherwise return false
		if ( id == 13 ) 
		{
			fastjet::PseudoJet muon = particles[j];

			// Check if muon is within range. Otherwise return false
			if ( abs(muon.pt()) > muonMinPt && 
			     abs(muon.eta()) < muonMaxEta ) 
			{
				// Check all other particles in event: add up their transverse
				// energies in a cone of radius deltaR_IsolatedLepton
				double sumEtInCone = 0.;
				for (int i=0; i<particles.size(); i++) 
				{
					// Only consider final state particles of the event record
					if (particles[i].isFinal())
					{
						fastjet::PseudoJet x = particles[i];

						// Only count visible particles that are not our muon in question
						if (i != j && 
						  particles[i].idAbs() != 12 && 
						  particles[i].idAbs() != 14 && 
						  particles[i].idAbs() != 16 && 
						  particles[i].idAbs() != 1000022 &&
						  particles[i].idAbs() != 8880022 &&
						  abs(particles[i].eta()) < MaxEta) 
						{
							double deltaEta = x.eta() - muon.eta();
							double deltaPhi = x.phi() - muon.phi();
							double deltaR   = sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
							if (deltaR < deltaR_IsolatedLepton)
								sumEtInCone += abs(x.pt());

						} // end if(visible and not considered muon)

					} // end if(final state)

				} // end for(particles.size())

				// isolation criterium: pT within cone less than sumEtInCone_IsolatedMuon GeV
				if (sumEtInCone < sumEtInCone_IsolatedMuon) 
					return true;

			} // end if(pTmin, etamax)

		} // end if(muon)

		return false;
	}

	bool jet_analysis::isolatedPhoton(const int & j, const Event & particles) 
	{
		long id = particles[j].idAbs();

		// Enter the loop if particles[j] is a photon. Otherwise return false
		if ( id == 22 ) 
		{
			fastjet::PseudoJet photon = particles[j];

			// Check if photon is within range. Otherwise return false
			if ( abs(photon.pt()) > photonMinPt && 
			     abs(photon.eta()) < photonMaxEta ) 
			{
				// Check all other particles in event: add up their transverse
				// energies in a cone of radius deltaR_IsolatedPhoton
				double sumEtInCone = 0.;
				for (int i=0; i<particles.size(); i++) 
				{
					// Only consider final state particles of the event record
					if (particles[i].isFinal())
					{
						fastjet::PseudoJet x = particles[i];

						// Only count visible particles that are not our photon in question
						if (i != j && 
						  particles[i].idAbs() != 12 && 
						  particles[i].idAbs() != 14 && 
						  particles[i].idAbs() != 16 && 
						  particles[i].idAbs() != 1000022 &&
						  particles[i].idAbs() != 8880022 &&
						  abs(particles[i].eta()) < MaxEta) 
						{
							double deltaEta = x.eta() - photon.eta();
							double deltaPhi = x.phi() - photon.phi();
							double deltaR   = sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
							if (deltaR < deltaR_IsolatedPhoton)
								sumEtInCone += abs(x.pt());

						} // end if(visible and not considered photon)

					} // end if(final state)

				} // end for(particles.size())

				// isolation criterium: pT within cone less than sumEtInCone_IsolatedPhoton GeV
				if (sumEtInCone < sumEtInCone_IsolatedPhoton) 
					return true;

			} // end if(pTmin, etamax)

		} // end if(photon)

		return false;
	}

	bool jet_analysis::JetElectronOverlapping(const fastjet::PseudoJet & jet, const std::vector< fastjet::PseudoJet > & leptons) 
	{
		if (leptons.size() == 0)
			return false;

		for (unsigned int i = 0; i < leptons.size(); ++i)
		{
			double deltaR = jet.delta_R(leptons[i]);
			if ( leptons[i].user_info<Pythia8::Particle>().idAbs() == 11 && deltaR < 0.2 )
				return true;
		}

		return false;
	}


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
			cout << "\nRan HEPTopTagger with the following parameters:" << endl; 
			top_tagger.get_setting();
			cout << endl;
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

	void jet_analysis::reduce_sample(cuts cut_list)
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


	// ========================================== // 
	//	         Analysis core function 	   	  //
	// ========================================== //

	void jet_analysis::initialise(fastjet::PseudoJet (jet_analysis::*TopTagger)(const fastjet::PseudoJet &))
	{
		//============= Initialisation =============//
		Pythia pythia;

		// Set the merging procedure if required
		if (DoMerging)
		{	
			pythia.readString("Merging:doKTMerging = on");

			std::string tms_string = "Merging:TMS = " + boost::lexical_cast<std::string>(MergingTMS);
			pythia.settings.readString(tms_string.c_str());

			std::string nj_string = "Merging:nJetMax = " + boost::lexical_cast<std::string>(MergingNJetMax);
			pythia.settings.readString(nj_string.c_str());

			std::string prc_string = "Merging:Process = " + MergingProcess;
			pythia.settings.readString(prc_string.c_str());
		}		

		// ===== TODO: modify for merging procedure ===== //
		// Initialise lhe input files
		if ( lhe_input.size() != 0 )
			pythia.init(lhe_input[0].c_str());
		else
			cout << "Error: no input file specified." << endl;

		// Select FastJet algorithm and parameters
		std::vector <fastjet::PseudoJet>  isolLeptons;
		std::vector <fastjet::PseudoJet>  isolPhotons;
		std::vector <fastjet::PseudoJet>  fjInputs;
	
		fastjet::RecombinationScheme      recombScheme = fastjet::E_scheme;
		fastjet::Strategy                 strategy = fastjet::Best;

		fastjet::JetDefinition            *fatJetDef;
		fastjet::JetAlgorithm             cambridge = fastjet::cambridge_algorithm;
		fatJetDef = new fastjet::JetDefinition(cambridge, Rsize_fat, recombScheme, strategy);

		fastjet::JetDefinition            *skinnyJetDef;
		fastjet::JetAlgorithm             antikt = fastjet::antikt_algorithm;
		skinnyJetDef = new fastjet::JetDefinition(antikt, Rsize_skinny, recombScheme, strategy);


		//============= Event loop =============//
		for (int iEvent = 0; iEvent < nEvent; ++iEvent)
		{
			// Generate event
			if (!pythia.next())
			{
				if( pythia.info.atEndOfFile() ) break;
				else continue;
			}

			// Reset FastJet inputs
			isolLeptons.clear();
			isolPhotons.clear();
			fjInputs.clear();

			// Define Etmiss vector
			double pxmiss = 0.0;
			double pymiss = 0.0;

			// Index vector to take trace of outgoing b-quarks of the hardest process
			std::vector< int > index;


			//============= Collect Isolated particles and FastJet input within each event=============//
			for (int i = 0; i < pythia.event.size(); ++i) 
			{ 
				// Identify the index of the outgoing b-quarks of the hardest process
				long pdg = pythia.event[i].idAbs();
				int status = pythia.event[i].status();
				if (pdg == 5 && status == -23)
					index.push_back(i);

				// Final state only
				if (!pythia.event[i].isFinal()) 
					continue;

				// Don't include particles that go down the beampipe
				if ( pythia.event[i].eta() > MaxEta )
					continue;

				// Don't include invisible particles (neutrinos neutralinos, heavy photon)
				long id = pythia.event[i].idAbs();
				if ( id == 12 || id == 14 || id == 16 || id == 1000022 || id == 8880022 ) 
					continue;

				// Take trace of Etmiss vector components from visible particles
				pxmiss += pythia.event[i].px();
				pymiss += pythia.event[i].py();

				// Collect isolated leptons
				if ( isolatedElectron(i, pythia.event) || 
					 isolatedMuon(i, pythia.event) ) 
				{
					isolLeptons.push_back(pythia.event[i]);
					continue;
				}

				// Collect isolated photons
				if ( isolatedPhoton(i, pythia.event) ) 
				{
					isolPhotons.push_back(pythia.event[i]);
					continue;
				}

				// Collect FastJet input
				fjInputs.push_back(pythia.event[i]);

			} // End of particle loop

			// Print Warning
			if ( fjInputs.size() == 0 ) 
			{
				cout << "Event no. " << iEvent+1 << " with no final state particles." << endl;
				continue;
			}

			// Sort isolLeptons, isolPhotons by pT
			isolLeptons = sorted_by_pt( isolLeptons );
			isolPhotons = sorted_by_pt( isolPhotons );


			//============= Run FastJet algorithm on fatJets and skinnyJets =============//
			fastjet::ClusterSequence *CSfatJets; 
			CSfatJets = new fastjet::ClusterSequence(fjInputs, *fatJetDef);

			fastjet::ClusterSequence *CSskinnyJets;
			CSskinnyJets = new fastjet::ClusterSequence(fjInputs, *skinnyJetDef);

			// Print FastJet details
			if (firstEvent)
			{
				cout << "\nSkinny Jets clustering:\n" 
					 << skinnyJetDef->description() << ". ";
				cout << "Strategy adopted by FastJet was "
				     << CSskinnyJets->strategy_string() << "." << endl;

				cout << "\nFat Jets clustering:\n" 
					 << fatJetDef->description() << ". ";
				cout << "Strategy adopted by FastJet was " 
					 << CSfatJets->strategy_string() << "." << endl;
			}

			// Extract inclusive jets sorted by pT (note minimum pT veto)
			std::vector< fastjet::PseudoJet >  skinnyJets;
			std::vector< fastjet::PseudoJet >  fatJets;
			fatJets = sorted_by_pt( CSfatJets->inclusive_jets(jetMinPt) );
			skinnyJets = sorted_by_pt( CSskinnyJets->inclusive_jets(jetMinPt) );

			// ===== TODO: implement b-flavour information ===== //
			// // Determine b-tag information of skinnyJets
			// for (unsigned int i = 0; i < skinnyJets.size(); ++i)
			// {
			// 	std::vector< fastjet::PseudoJet > pieces = skinnyJets[i].constituents();
			// 	for (unsigned int j = 0; j < pieces.size(); ++j)
			// 	{
			// 		bool isB = false;
			// 		for (unsigned int k = 0; k < index.size(); ++k)
			// 		{
			// 	        int iAncestor = index[k];
			// 	        int iPos = pieces[j].user_info<Pythia8::Particle>().mother1();

			// 	        bool isMatched = pythia.event.isAncestor(iPos,iAncestor);

			// 	        // if a b-quark ancestor is found: set b-tag information, remove b-index and jump at the next skinnyJet (isB = true)
			// 	        if(isMatched)
			// 	        {
			// 	        	skinnyJets[i].set_user_info(new FlavourInfo(true));
			// 	        	isB = true;
			// 				index.erase(index.begin()+k);
			// 				break;
			// 	        }
			// 		} // end index loop

			// 		if (isB)
			// 			break;

			// 		// if no b-quark is found among the ancestors of the skinnyJet constituents, then set no-b information	
			// 		if (j == (pieces.size()-1))
			// 			skinnyJets[i].set_user_info(new FlavourInfo(false));

			// 	} // end pieces loop

			// } // end skinnyJets loop

			// // Check b-quark history of final state particles
			// for (unsigned int i = 0; i < fjInputs.size(); ++i)
			// {
			// 	for (unsigned int k = 0; k < index.size(); ++k)
			// 	{
			//         int iAncestor = index[k];
			//         int iPos = fjInputs[i].user_info<Pythia8::Particle>().mother1();

			//         bool isMatched = pythia.event.isAncestor(iPos,iAncestor);

			//         if(isMatched)
			//         {
			//         	if(firstEvent)
			//         		cout << "found ancestor, with id: " << iAncestor << endl;;
			//         }
			// 	} // end index loop

			// } // end fjInputs loop

			// Temporary null-assignment of b-tag information
			for (unsigned int i = 0; i < skinnyJets.size(); ++i)
				skinnyJets[i].set_user_info(new FlavourInfo(false));


			//============= Store clustered events in lhco format =============//
			event *ev = new event;

			// Translate isolLeptons into lhco format and push back into the list of event pointers
			for (unsigned int i = 0; i < isolLeptons.size(); ++i)
			{
				int p_type;
				int 	p_id 	= isolLeptons[i].user_info<Pythia8::Particle>().idAbs();
				double	p_eta 	= isolLeptons[i].eta(), 
						p_phi 	= isolLeptons[i].phi(), 
						p_pt 	= isolLeptons[i].pt(),
						p_ch	= isolLeptons[i].user_info<Pythia8::Particle>().charge();

				if (p_id == 11)
					p_type = particle::type_electron;
				else
					p_type = particle::type_muon;

				lhco *p = new lhco(p_type, p_eta, p_phi, p_pt, 0.0, p_ch);
				ev->push_back(p);
			}

			// Translate isolPhotons into lhco format and push back into the list of event pointers
			for (unsigned int i = 0; i < isolPhotons.size(); ++i)
			{
				int p_type = particle::type_photon;
				double	p_eta 	= isolPhotons[i].eta(), 
						p_phi 	= isolPhotons[i].phi(), 
						p_pt 	= isolPhotons[i].pt();

				lhco *p = new lhco(p_type, p_eta, p_phi, p_pt);
				ev->push_back(p);
			}

			// Translate skinnyJets into lhco format and push back into the list of event pointers
			for (unsigned int i = 0; i < skinnyJets.size(); ++i)
			{
				// no Jet-Electron overlapping within R=0.2 cone from the electron
				if (!JetElectronOverlapping(skinnyJets[i],isolLeptons))
				{
					int p_type = particle::type_jet;
					double	p_eta 	= skinnyJets[i].eta(), 
							p_phi 	= skinnyJets[i].phi(), 
							p_pt 	= skinnyJets[i].pt(), 
							p_m 	= skinnyJets[i].m(),
							p_bjet	= skinnyJets[i].user_info<FlavourInfo>().b_type();

					lhco *p = new lhco(p_type, p_eta, p_phi, p_pt, p_m, 0.0, p_bjet);
					ev->push_back(p);
				}
			}

			// Translate Etmiss into lhco format and push back into the list of event pointers
			int p_type = particle::type_met;
			double Etmiss = sqrt( pxmiss*pxmiss + pymiss*pymiss );
			double Etmiss_phi = atan2(pymiss,pxmiss);
			lhco *p = new lhco(p_type, 0.0, Etmiss_phi, Etmiss);
			ev->push_back(p);

			// First event only
			if (firstEvent)
			{
				cout << "\nFirst event details:" << endl;
  				ev->write(cout);
  			}


			//============= Tagging analysis =============//
			std::vector <fastjet::PseudoJet>  taggedJets;
  			for (unsigned int i = 0; i < fatJets.size(); ++i)
  			{
  				// Top tagging
  				if (DoTopTagger)
  				{
	  				fastjet::PseudoJet TopTagged;
	  				TopTagged=(this->*TopTagger)(fatJets[i]);
	  				if ( printTopTagDetails )
	  					printTopTagDetails = false;
	  				if ( TopTagged!=0 )
	  				{
	  					taggedJets.push_back(TopTagged);
	  					continue;
	  				}
  				}
  				
  				// BDRS tagging: only on non-top tagged fat jets
  				if (DoBDRS)
  				{
					fastjet::PseudoJet BDRSTagged;
					BDRSTagged = BDRSTagging(fatJets[i]);
					if ( printBDRSDetails )
						printBDRSDetails = false;
					if ( BDRSTagged!=0 )
						taggedJets.push_back(BDRSTagged);
				}

  			}


  			//============= Fill the (lhco, taggedJets) map =============//
			map_lhco_taggedJets.insert(std::make_pair(ev,taggedJets));

			// Print evolution
			if ((iEvent+1)%100 == 0)
				cout << "Event no.: " << iEvent+1 << "\r" << flush;

			// Exit First event loop
			if (firstEvent)
		   		firstEvent = false;

			delete CSfatJets;
			delete CSskinnyJets;

		} // End of event loop

		// // Print statistics
		// pythia.statistics();

		cout << "\n" << endl;		

		delete fatJetDef;
		delete skinnyJetDef;

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
