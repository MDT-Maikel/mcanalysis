/* Jet analysis class: Initialisation function
 *
 * 
*/

#include "jet_analysis.h"


/* NAMESPACE */
namespace analysis 
{

	// ========================================== // 
	//	         Initialisation function 	   	  //
	// ========================================== //

	void jet_analysis::initialise(fastjet::PseudoJet (jet_analysis::*TopTagger)(const fastjet::PseudoJet &))
	{
		//============= Initialisation =============//
		Pythia8::Pythia pythia;
		pythia.readString("Print:quiet=ON");
		if (fast_showering)
		{
			pythia.readString("PartonLevel:MPI=OFF");
			pythia.readString("PartonLevel:Remnants=OFF");
			pythia.readString("Check:Event=OFF");
			pythia.readString("HadronLevel:all=OFF");
		}

		// Set the merging procedure if required
		if (DoMerging)
		{	
			pythia.readString("Merging:doKTMerging = on");
			pythia.settings.readString(("Merging:TMS = " + boost::lexical_cast<std::string>(MergingTMS)).c_str());
			pythia.settings.readString(("Merging:nJetMax = " + boost::lexical_cast<std::string>(MergingNJetMax)).c_str());
			pythia.settings.readString(("Merging:Process = " + MergingProcess).c_str());
		}		

		// ===== TODO: modify for merging procedure ===== //
		// Initialise lhe input files
		if ( lhe_input.size() != 0 )
			pythia.init(lhe_input[0].c_str());
		else
			cout << "Error: no input file specified." << endl;

		// Initialise lhco input files
		vector<event*> events_lhco;
		if ( importedLHCO )
			read_lhco(events_lhco, lhco_input);

		// Select FastJet algorithm and parameters
		std::vector <fastjet::PseudoJet>  isolLeptons;
		std::vector <fastjet::PseudoJet>  isolPhotons;
		std::vector <fastjet::PseudoJet>  fjInputs;
	
		fastjet::RecombinationScheme      recombScheme = fastjet::E_scheme;
		fastjet::Strategy                 strategy = fastjet::Best;

		fastjet::JetDefinition            *fatJetDef;
		fatJetDef = new fastjet::JetDefinition(algorithm_fat, Rsize_fat, recombScheme, strategy);

		fastjet::JetDefinition            *skinnyJetDef;
		skinnyJetDef = new fastjet::JetDefinition(algorithm_skinny, Rsize_skinny, recombScheme, strategy);


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
				if (isolatedElectron(i, pythia.event) || isolatedMuon(i, pythia.event)) 
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


			//============= Run FastJet algorithm on fatJets and (if required) on skinnyJets =============//
			fastjet::ClusterSequence *CSfatJets; 
			CSfatJets = new fastjet::ClusterSequence(fjInputs, *fatJetDef);

			fastjet::ClusterSequence *CSskinnyJets;
			CSskinnyJets = new fastjet::ClusterSequence(fjInputs, *skinnyJetDef);

			// Print FastJet details
			if (firstEvent)
			{
				if ( !importedLHCO )
				{
					cout << "\nSkinny Jets clustering:\n" 
						 << skinnyJetDef->description() << ". ";
					cout << "Strategy adopted by FastJet was "
					     << CSskinnyJets->strategy_string() << "." << endl;
				}

				cout << "\nFat Jets clustering:\n" 
					 << fatJetDef->description() << ". ";
				cout << "Strategy adopted by FastJet was " 
					 << CSfatJets->strategy_string() << "." << endl;
			}

			// Extract inclusive jets sorted by pT (note minimum pT veto)
			std::vector< fastjet::PseudoJet >  fatJets;
			fatJets = fastjet::sorted_by_pt( CSfatJets->inclusive_jets(jetMinPt) );

			std::vector< fastjet::PseudoJet >  skinnyJets;
			if ( !importedLHCO )
				skinnyJets = fastjet::sorted_by_pt( CSskinnyJets->inclusive_jets(jetMinPt) );

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


			//============= If requested, store clustered events in lhco format, otherwise read from input lhco =============//
			event *ev = new event;

			if ( !importedLHCO )
			{
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
						p_type = ptype_electron;
					else
						p_type = ptype_muon;

					lhco *p = new lhco(p_type, p_eta, p_phi, p_pt, 0.0, p_ch);
					ev->push_back(p);
				}

				// Translate isolPhotons into lhco format and push back into the list of event pointers
				for (unsigned int i = 0; i < isolPhotons.size(); ++i)
				{
					int p_type = ptype_photon;
					double	p_eta 	= isolPhotons[i].eta(), 
							p_phi 	= isolPhotons[i].phi(), 
							p_pt 	= isolPhotons[i].pt();

					lhco *p = new lhco(p_type, p_eta, p_phi, p_pt);
					ev->push_back(p);
				}

				// Temporary null-assignment of b-tag information
				for (unsigned int i = 0; i < skinnyJets.size(); ++i)
					skinnyJets[i].set_user_info(new FlavourInfo(false));

				// Translate skinnyJets into lhco format and push back into the list of event pointers
				for (unsigned int i = 0; i < skinnyJets.size(); ++i)
				{
					// no Jet-Electron overlapping within R=0.2 cone from the electron
					if (!JetElectronOverlapping(skinnyJets[i],isolLeptons))
					{
						int p_type = ptype_jet;
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
				int p_type = ptype_met;
				double Etmiss = sqrt( pxmiss*pxmiss + pymiss*pymiss );
				double Etmiss_phi = atan2(pymiss,pxmiss);
				lhco *p = new lhco(p_type, 0.0, Etmiss_phi, Etmiss);
				ev->push_back(p);

  			}
  			else
  				ev = events_lhco[iEvent];

  			// Print details of first event
			if (firstEvent)
			{
				cout << endl << "First event details:" << endl;
  				ev->write(cout);
  				cout << endl;
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
					{
						taggedJets.push_back(BDRSTagged);
						continue;
					}
				}

				// push_back remaining non-tagged fatjets
				fastjet::PseudoJet fatjet = fatJets[i];
				fatjet.set_user_info(new TagInfo());
				taggedJets.push_back(fatjet);

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


/* NAMESPACE */
}