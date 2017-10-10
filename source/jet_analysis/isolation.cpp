/* Jet analysis class: Isolation functions
 *
 * 
*/

#include "jet_analysis.h"


/* NAMESPACE */
namespace analysis 
{

	/* isolation functions */

	bool jet_analysis::isolatedElectron(const int & j, const Pythia8::Event & particles) 
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
				// energies in a cone of radius deltaR_IsoEl
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
							abs(particles[i].eta()) < MaxEta &&
							particles[i].pT() > pTminTrack_IsoEl
							) 
						{
							double deltaEta = x.eta() - electron.eta();
							double deltaPhi = abs(x.phi() - electron.phi());
							deltaPhi = std::min(deltaPhi, 8 * atan(1) - deltaPhi);
							double deltaR = sqrt( pow(deltaEta,2.0) + pow(deltaPhi,2.0) );

							if (deltaR < deltaR_IsoEl)
								sumEtInCone += abs(x.pt());

						} // end if(visible and not considered electron)

					} // end if(final state)

				} // end for(particles.size())

				// isolation criterium: pT within cone less than pTfracMax_IsoEl% of pT(e)
				if (sumEtInCone < pTfracMax_IsoEl*abs(electron.pt())) 
					return true;

			} // end if(pTmin, etamax)

		} // end if(electron)

		return false;
	}

	bool jet_analysis::isolatedMuon(const int & j, const Pythia8::Event & particles) 
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
				// energies in a cone of radius deltaR_IsoMuon
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
							abs(particles[i].eta()) < MaxEta &&
							particles[i].pT() > pTminTrack_IsoMuon
							) 
						{
							double deltaEta = x.eta() - muon.eta();
							double deltaPhi = abs(x.phi() - muon.phi());
							deltaPhi = std::min(deltaPhi, 8 * atan(1) - deltaPhi);
							double deltaR = sqrt( pow(deltaEta,2.0) + pow(deltaPhi,2.0) );

							if (deltaR < deltaR_IsoMuon)
								sumEtInCone += abs(x.pt());

						} // end if(visible and not considered muon)

					} // end if(final state)

				} // end for(particles.size())

				// isolation criterium: isolation criterium: pT within cone less than pTfracMax_IsoMuon% of pT(mu)
				if (sumEtInCone < pTfracMax_IsoMuon*abs(muon.pt())) 
					return true;

			} // end if(pTmin, etamax)

		} // end if(muon)

		return false;
	}

	bool jet_analysis::isolatedPhoton(const int & j, const Pythia8::Event & particles) 
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
				// energies in a cone of radius deltaR_IsoGamma
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
							abs(particles[i].eta()) < MaxEta &&
							particles[i].pT() > pTminTrack_IsoGamma
							) 
						{
							double deltaEta = x.eta() - photon.eta();
							double deltaPhi = abs(x.phi() - photon.phi());
							deltaPhi = std::min(deltaPhi, 8 * atan(1) - deltaPhi);
							double deltaR = sqrt( pow(deltaEta,2.0) + pow(deltaPhi,2.0) );

							if (deltaR < deltaR_IsoGamma)
								sumEtInCone += abs(x.pt());

						} // end if(visible and not considered photon)

					} // end if(final state)

				} // end for(particles.size())

				// isolation criterium: pT within cone less than pTfracMax_IsoGamma% of pT(gamma)
				if (sumEtInCone < pTfracMax_IsoGamma*abs(photon.pt())) 
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

/* NAMESPACE */
}
