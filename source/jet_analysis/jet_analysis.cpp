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

		// lhe input & number of events to be analysed by Phytia
		lhe_name = "";
		nEvent = 1000;
		firstEvent = true;
		printTopTagDetails = true;
		printBRDSDetails = true;

		// Jet clustering parameters
		Rsize_fat = 1.4;
		Rsize_skinny = 0.4;

		// JHTopTagger parameters
		DoTopTagger = true;
		JHTopTagger_delta_p = 0.05; 
		JHTopTagger_delta_r = 0.19;
	}

	jet_analysis::~jet_analysis()
	{
	}


	// ========================================== // 
	//	       Set parameters and lhe input   	  //
	// ========================================== //

	void jet_analysis::set_Isolation(const std::string & type, const double & eta, const double & pt, const double & Rcone, const double & EtCone) //TODO: update using p_type
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
		{
			std::cout << "Wrong type assignment." << std::endl;
		}
	}

	void jet_analysis::set_nEvent(const int & events)
	{
		nEvent = events;
	}

	void jet_analysis::add_lhe(const std::string & name)
	{
		lhe_name = name;
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
		JHTopTagger_delta_p = delta_p;
		JHTopTagger_delta_r = delta_r;
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
						  particles[i].idAbs() != 1000022 
						  && abs(particles[i].eta()) < MaxEta) 
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
						  particles[i].idAbs() != 1000022 
						  && abs(particles[i].eta()) < MaxEta) 
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
						  particles[i].idAbs() != 1000022 
						  && abs(particles[i].eta()) < MaxEta) 
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
		{
			return false;
		}

		for (unsigned int i = 0; i < leptons.size(); ++i)
		{
			double deltaR   = jet.delta_R(leptons[i]);
			if (/*leptons[i].idAbs() == 11 &&*/ deltaR < 0.2) //TODO: add pdg information
			{
				return true;
			}
		}

		return false;
	}


	// ========================================== // 
	//	       	 Tagging functions	 		   	  //
	// ========================================== //

	fastjet::PseudoJet jet_analysis::JHTopTagging(const fastjet::PseudoJet & jet)
	{
		// Tagging algorithm
		fastjet::JHTopTagger top_tagger(JHTopTagger_delta_p, JHTopTagger_delta_r);
		top_tagger.set_top_selector(fastjet::SelectorMassRange(150,200));
		top_tagger.set_W_selector  (fastjet::SelectorMassRange( 65, 95));
		fastjet::PseudoJet tagged = top_tagger(jet);

		if (printTopTagDetails)
		{
			// Print Top-tagging details
			cout << "\nRan the following top tagger:\n" 
		   		 << top_tagger.description() << endl;
		}

		// Set TagInfo
		if ( tagged!=0 )
		{
			bool top_tag = true;
			bool w_tag = false;
			bool h_tag = false;
			tagged.set_user_info(new TagInfo(top_tag, w_tag, h_tag));
		}

		return tagged;
	}

	fastjet::PseudoJet jet_analysis::BRDSTagging(const fastjet::PseudoJet & jet)
	{
		// from example 12 of FastJet

		// Definition of Mass drop tagger with \mu=0.667 and ycut=0.09
		fastjet::MassDropTagger md_tagger(0.667, 0.09);
		fastjet::PseudoJet tagged = md_tagger(jet);

		if (printBRDSDetails)
		{
			// Print BRDS-tagging details
			cout << "\nRan the following BRDS tagger:\n" 
		   		 << md_tagger.description() << endl << endl;
		}

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
			if ( 65. < invMass && invMass < 95. )
			{
				bool top_tag = false;
				bool w_tag = true;
				bool h_tag = false;
				filtered.set_user_info(new TagInfo(top_tag, w_tag, h_tag));
			}
			else if ( 100. < invMass && invMass < 130. )
			{
				bool top_tag = false;
				bool w_tag = false;
				bool h_tag = true;
				filtered.set_user_info(new TagInfo(top_tag, w_tag, h_tag));
			}
			else
			{
				bool top_tag = false;
				bool w_tag = false;
				bool h_tag = false;
				filtered.set_user_info(new TagInfo(top_tag, w_tag, h_tag));
			}

			return filtered;
		}

		return tagged;
	}


	// ========================================== // 
	//	          Apply cut function 		   	  //
	// ========================================== //

	void jet_analysis::reduce_sample(cuts cut_list)
	{	
		// map_lhco_fatJets iterator definition
		std::map< event *, std::vector <fastjet::PseudoJet> >::iterator it;

		// extract vector of event pointers "events"
		std::vector< event * > events;
		for (it=map_lhco_fatJets.begin(); it!=map_lhco_fatJets.end(); ++it)
		{
			events.push_back(it->first);
		}

		// perform the cuts defined in cut_list on "events"
		cut_list.apply(events);
		cut_list.write(cout);

		// retrieve from map_lhco_fatJets only the elements which passed the cuts
		std::map< event *, std::vector <fastjet::PseudoJet> > reduced_map;
		for (unsigned int i = 0; i < events.size(); ++i)
		{
			it = map_lhco_fatJets.find(events[i]);
			if ( it!=map_lhco_fatJets.end() )
			{
				reduced_map.insert( 
					std::make_pair( events[i], 
									map_lhco_fatJets[events[i]] ) 
					); // end insert
			}
		}

		// resize the map_lhco_fatJets
		map_lhco_fatJets = reduced_map;
	}


	// ========================================== // 
	//	         Analysis core function 	   	  //
	// ========================================== //

	void jet_analysis::initialise(fastjet::PseudoJet (jet_analysis::*TopTagger)(const fastjet::PseudoJet &))
	{
		// Initialise Pythia on LHE file
		Pythia pythia;
		// pythia.init("./w2emveb_light_lhc8_0.lhe");
		pythia.init(lhe_name.c_str());

		// Fastjet analysis: select algorithm and parameters
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

		// Event loop
		for (int iEvent = 0; iEvent < nEvent; ++iEvent)
		{
			// Generate event
			if (!pythia.next()) continue;

			// Reset Fastjet inputs
			isolLeptons.clear();
			isolPhotons.clear();
			fjInputs.clear();

			// Loop over event record to decide what to pass to FastJet
			for (int i = 0; i < pythia.event.size(); ++i) 
			{ 
				// Final state only
				if (!pythia.event[i].isFinal()) 
					continue;

				// Don't include particles that go down the beampipe
				if ( pythia.event[i].eta() > MaxEta )
					continue;

				// Don't include neutrinos or neutralinos
				long id = pythia.event[i].idAbs();
				if ( id == 12 || id == 14 || id == 16 || id == 1000022 ) 
					continue;

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

				// Collect Fastjet input, with TagInfo
				fastjet::PseudoJet jetCandidate;
				jetCandidate = pythia.event[i];
				// jetCandidate.set_user_info(new TagInfo()); //TODO: insert b-tag info
				fjInputs.push_back(jetCandidate);
			}

			// Print Warning
			if ( fjInputs.size() == 0 ) 
				{
					cout << "Event no. " << iEvent+1 << " with no final state particles." << endl;
					continue;
				}

			// Sort isolLeptons, isolPhotons by pT
			isolLeptons = sorted_by_pt( isolLeptons );
			isolPhotons = sorted_by_pt( isolPhotons );

			// Run Fastjet algorithm for fatJets and skinnyJets
			fastjet::ClusterSequence *CSfatJets; 
			CSfatJets = new fastjet::ClusterSequence(fjInputs, *fatJetDef);

			fastjet::ClusterSequence *CSskinnyJets;
			CSskinnyJets = new fastjet::ClusterSequence(fjInputs, *skinnyJetDef);

			// First event only
			if (firstEvent)
			{
				// Print FastJet details
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
			std::vector <fastjet::PseudoJet>  skinnyJets;
			std::vector <fastjet::PseudoJet>  fatJets;
			fatJets = sorted_by_pt( CSfatJets->inclusive_jets(jetMinPt) );
			skinnyJets = sorted_by_pt( CSskinnyJets->inclusive_jets(jetMinPt) );


			//=== Store clustered events ===//
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
				{
					p_type = particle::type_electron;
				}
				else
				{
					p_type = particle::type_muon;
				}

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
							p_m 	= skinnyJets[i].m();

					lhco *p = new lhco(p_type, p_eta, p_phi, p_pt, p_m);
					ev->push_back(p);
				}
			}

			// Translate Etmiss into lhco format and push back into the list of event pointers
			int p_type = particle::type_met;
			double Etmiss, pxmiss = 0.0, pymiss = 0.0;
			for (unsigned int i = 0; i < skinnyJets.size(); i++) 
			{
				pxmiss += skinnyJets[i].px();        
				pymiss += skinnyJets[i].py();
			}
			for (unsigned int i = 0; i < isolLeptons.size(); i++) 
			{
				pxmiss += isolLeptons[i].px();
				pymiss += isolLeptons[i].py();
			}
			for (unsigned int i = 0; i < isolPhotons.size(); i++) 
			{
				pxmiss += isolPhotons[i].px();
				pymiss += isolPhotons[i].py();
			}
			Etmiss = sqrt( pxmiss*pxmiss + pymiss*pymiss );
			double Etmiss_phi = atan2(pymiss,pxmiss);
			lhco *p = new lhco(p_type, 0.0, Etmiss_phi, Etmiss);
			ev->push_back(p);

			// First event only
			if (firstEvent)
			{
				cout << "\nFirst event details:" << endl;
  				ev->write(cout);
  			}

			// Tagging analysis
			std::vector <fastjet::PseudoJet>  taggedJets;
  			for (unsigned int i = 0; i < fatJets.size(); ++i)
  			{
  				// Top tagging
  				fastjet::PseudoJet TopTagged;
  				TopTagged=(this->*TopTagger)(fatJets[i]);
  				if ( printTopTagDetails )
  					printTopTagDetails = false;
  				if ( TopTagged!=0 )
  				{
  					taggedJets.push_back(TopTagged);
  					continue;
  				}
  				
  				// BRDS tagging: only on non-top tagged fat jets
				fastjet::PseudoJet BRDSTagged;
				BRDSTagged = BRDSTagging(fatJets[i]);
				if ( printBRDSDetails )
					printBRDSDetails = false;
				if ( BRDSTagged!=0 )
					taggedJets.push_back(BRDSTagged);

  			}

			// Fill the (lhco, taggedJets) map
			map_lhco_fatJets.insert(std::make_pair(ev,taggedJets));

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
	//	            TagInfo class 			   	  //
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
