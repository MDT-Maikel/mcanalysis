/* Invariant mass cut classes
 *
 * Provides a set of cuts related to invariant mass by overloading the 
 * cut base class.
*/

#include "gres_cuts.h"


/* NAMESPACE */
namespace analysis
{

	/* invariant jet combination mass cut class */

	cut_gres::cut_gres()
	{
		topology = {1, 2};		
	}

	void cut_gres::init()
	{
		// calculate the number of jets
		calc_nrjets();
		// remove single resonances from topology list
		remove_singleres();
		// calculate number of combinations
		calc_nrcombinations();
		// calculate all combinations
		calc_combinations();
		// calculate all pairing of resonance combinations
		calc_pairings();

		// LOGGING: topology
		std::cout << "CUT_GRES: reduced topology: {";
		for (unsigned int i = 0; i < reduced_topology.size(); i++)
		{	
			std::cout << reduced_topology[i];
			if (i != reduced_topology.size() - 1)
				std::cout << ", ";
		}
		std::cout << "} with Njet = " << nr_jets << std::endl;
	
		// LOGGING: combinations
		for (unsigned int i = 0; i < resonance_combinations.size(); i++)
		{
			std::vector<std::vector<unsigned int> > combinations = resonance_combinations[i];
			std::cout << "CUT_GRES: for resonance with Njets = " << reduced_topology[i] << " we have " << combinations.size() << " combinations: " << std::endl << "\t";		
			for (unsigned int j = 0; j < combinations.size(); j++)
			{
				std::vector<unsigned int> comb = combinations[j];
				std::cout << "{";
				for (unsigned int k = 0; k < comb.size(); k++)
				{	
					std::cout << comb[k];
					if (k != comb.size() - 1)
						std::cout << ", ";
				}
				std::cout << "}; ";
			}
			std::cout << std::endl;
		}

		// LOGGING: number of combinations and pairings
		std::cout << "CUT_GRES: number of combinations = " << nr_combinations << " with the following pairings: " << std::endl;
		for (unsigned int i = 0; i < resonance_pairings.size(); i++)
		{
			std::cout << "[";
			std::vector<std::vector<unsigned int> > pair = resonance_pairings[i];
			for (unsigned int j = 0; j < pair.size(); j++)
			{
				std::cout << "{";
				std::vector<unsigned int> res = pair[j];
				for (unsigned int k = 0; k < res.size(); k++)
				{	
					std::cout << res[k];
					if (k != res.size() - 1)
						std::cout << ", ";
				}
				std::cout << "}";

				if (j != pair.size() - 1)
					std::cout << "; ";

			}
			std::cout << "]" << std::endl;
		}		

	}

	/*bool cut_gres::passed(const event *ev) 
	{
		if (nr_jets == 3 && reduced_topology.size() == 1 && reduced_topology[0] == 2)
			return passed_top21(ev);
		if (nr_jets == 4 && reduced_topology.size() == 1 && reduced_topology[0] == 3)
			return passed_top31(ev);
		if (nr_jets == 4 && reduced_topology.size() == 2 && reduced_topology[0] == 2 && reduced_topology[1] == 2)
			return passed_top22(ev);
	
		return false;
	}*/

	bool cut_gres::operator() (const event *ev)
	{
		// event did not pass by default
		bool has_passed = false;

		// loop over all resonance pairings
		for (unsigned int i = 0; i < nr_combinations; ++i)
		{
			std::vector< std::vector<unsigned int> > pair = resonance_pairings[i];
			bool pair_passed = true;			
	
			// loop over all resonances in the pair
			for (unsigned int j = 0; j < pair.size(); j++)
			{
				std::vector<unsigned int> comb = pair[j];
				
				// construct the resonance event
				event *res = new event;
				for	(unsigned int k = 0; k < comb.size(); k++)
					res->push_back(ev->get(ptype_jet, comb[k], 2.8));

				// check whether the pair passed the cuts
				if (res->mass() < reduced_cuts[j])
					pair_passed = false;
	
			}
			
			// all pairs passed the cuts
			if (pair_passed)
			{
				has_passed = true;	
				break;
			}
		}

		return has_passed;		
	}

	// special case for topology {2,1}
	bool cut_gres::passed_top21(const event *ev)
	{
		bool has_passed = false;
		std::vector< std::vector<int> > combinations = {{1, 2}, {1, 3}, {2, 3}};
		for (unsigned int i = 0; i < combinations.size(); i++)
		{
			particle *p1 = ev->get(ptype_jet, combinations[i][0], 2.8);
			particle *p2 = ev->get(ptype_jet, combinations[i][1], 2.8);
			event *jets = new event;
			jets->push_back(p1); jets->push_back(p2);
			
			if (jets->mass() > reduced_cuts[0])
			{
				has_passed = true;
				break;
			}
		}
		return has_passed;		
	}

	// special case for topology {3,1}
	bool cut_gres::passed_top31(const event *ev)
	{
		bool has_passed = false;		
		std::vector< std::vector<int> > combinations = {{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}};
		for (unsigned int i = 0; i < combinations.size(); i++)
		{
			particle *p1 = ev->get(ptype_jet, combinations[i][0], 2.8);
			particle *p2 = ev->get(ptype_jet, combinations[i][1], 2.8);
			particle *p3 = ev->get(ptype_jet, combinations[i][2], 2.8);
			event *jets = new event;
			jets->push_back(p1); jets->push_back(p2); jets->push_back(p3);
			
			if (jets->mass() > reduced_cuts[0])
			{
				has_passed = true;
				break;
			}
		}
		return has_passed;
	}

	// special case for topology {2,2}
	bool cut_gres::passed_top22(const event *ev)
	{
		bool has_passed = false;
		std::vector< std::vector<int> > comb1 = {{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}};
		std::vector< std::vector<int> > comb2 = {{3, 4}, {2, 4}, {2, 3}, {1, 4}, {1, 3}, {1, 2}};
		for (unsigned int i = 0; i < comb1.size(); i++)
		{
			event *res1 = new event;
			particle *p1 = ev->get(ptype_jet, comb1[i][0], 2.8);
			particle *p2 = ev->get(ptype_jet, comb1[i][1], 2.8);
			res1->push_back(p1); res1->push_back(p2);

			event *res2 = new event;
			particle *p3 = ev->get(ptype_jet, comb2[i][0], 2.8);
			particle *p4 = ev->get(ptype_jet, comb2[i][1], 2.8);
			res2->push_back(p3); res2->push_back(p4);
			
			if (res1->mass() > reduced_cuts[0] && res2->mass() > reduced_cuts[1])
			{
				has_passed = true;
				break;
			}
		}
		return has_passed;
	}

	std::string cut_gres::name() const
	{
		return "invariant mass combination";
	}

	void cut_gres::set_topology(const std::vector<int> & top)
	{
		topology = top;
		reduced_topology = top;
	}

	void cut_gres::set_top_cuts(const std::vector<double> & cuts)
	{
		top_cuts = cuts;
		reduced_cuts = cuts;
	}

	void cut_gres::calc_nrjets()
	{
		nr_jets = 0;
		for (unsigned int i = 0; i < topology.size(); i++)
			nr_jets += topology[i];
	}

	void cut_gres::remove_singleres()
	{
		reduced_topology = topology;
		reduced_cuts = top_cuts;
		// remove single resonances from topology list
		for (int i = reduced_topology.size() - 1; i >= 0; --i)
		{
			if (reduced_topology[i] == 1)
			{
				reduced_topology.erase(reduced_topology.begin() + i);
				reduced_cuts.erase(reduced_cuts.begin() + i);
			}
		}
	}

	void cut_gres::calc_nrcombinations()
	{
		nr_combinations = 1;
		int count = nr_jets;
		for (unsigned int i = 0; i < reduced_topology.size(); ++i)
		{
			nr_combinations *= static_cast<unsigned int>(boost::math::binomial_coefficient<double>(count, reduced_topology[i]));
			count -= reduced_topology[i];
		}
	}

	void cut_gres::calc_combinations()
	{
		resonance_combinations.clear();		

		// make the list of jet combinations which is possible for each of the resonances in the topology
		for (unsigned int i = 0; i < reduced_topology.size(); i++)
		{
			unsigned int resonance = reduced_topology[i];
			unsigned int nr_combs = static_cast<unsigned int>(boost::math::binomial_coefficient<double>(nr_jets, resonance)); // TODO: Is this failsafe?
			std::vector<std::vector<unsigned int> > combinations;
			std::vector<unsigned int> start_comb;
			for (int j = resonance; j > 0; j--)
				start_comb.push_back(j);
			// produce all combinations which lead to the resonance
			for (unsigned int j = 0; j < nr_combs; j++)
			{
				combinations.push_back(start_comb);
				start_comb[0]++;
				for (unsigned int k = 0; k < start_comb.size() - 1; k++)
				{
					if (start_comb[k] > nr_jets - k)
					{
						start_comb[k + 1]++;
						for (int l = k; l >= 0; l--)
							start_comb[l] = start_comb[l + 1] + 1;
					}
				}
			}
			resonance_combinations.push_back(combinations);
		}
	}

	void cut_gres::calc_pairings()
	{
		resonance_pairings.clear();		

		// standard ones for the basic three
		/*if (nr_jets == 3 && reduced_topology.size() == 1 && reduced_topology[0] == 2)
		{
			resonance_pairings = {{{1, 2}}, {{1, 3}}, {{2, 3}}};
			return;
		}
		if (nr_jets == 4 && reduced_topology.size() == 1 && reduced_topology[0] == 3)
		{
			resonance_pairings = {{{1, 2, 3}}, {{1, 2, 4}}, {{1, 3, 4}}, {{2, 3, 4}}};
			return;
		}
		if (nr_jets == 4 && reduced_topology.size() == 2 && reduced_topology[0] == 2 && reduced_topology[1] == 2)
		{
			resonance_pairings = {{{1, 2}, {3, 4}}, {{1, 3}, {2, 4}}, {{1, 4}, {2, 3}}, {{2, 3}, {1, 4}}, {{2, 4}, {1, 3}}, {{3, 4}, {1, 2}}};
			return;
		}*/
	
		// a list with zero for each resonance
		std::vector<unsigned int> counter(reduced_topology.size(), 0);
	
		// loop over number of topology combinations
		// and check all combinations of resonance assignments
		unsigned int i = 0;
		while (i < nr_combinations)
		{
			// construct resonance pairing
			std::vector<std::vector<unsigned int> > pairing;
			for (unsigned int j = 0; j < reduced_topology.size(); j++)
				pairing.push_back(resonance_combinations[j][counter[j]]);

			// increase pairing counters
			for (int k = reduced_topology.size() - 1; k >= 0; k--)
			{
				counter[k]++;
				if (counter[k] >= resonance_combinations[k].size())
					counter[k] = 0;
				else
					break;
			}

			// test pairing and put into list if valid
			if (is_valid_pairing(pairing))
			{
				resonance_pairings.push_back(pairing);
				i++;
			}			
		}

	}

	// a pairing is valid if each jet only appears once
	bool cut_gres::is_valid_pairing(const std::vector<std::vector<unsigned int> > &pairing) const
	{
		// make a list of jets.
		std::vector<unsigned int> jets;
		for (unsigned int i = 0; i < pairing.size(); i++)
		{
			std::vector<unsigned int> pair = pairing[i];
			for (unsigned int j = 0; j < pair.size(); j++)
				jets.push_back(pair[j]);
		}
		
		bool is_valid = true;
		for (unsigned int i = 0; i < jets.size(); i++)
		{
			unsigned int compare = jets[i];
			for (unsigned int j = i + 1; j < jets.size(); j++)
			{
			 	if (jets[j] == compare)
				{	
					is_valid = false;
					break;
				}
			}
		}	

		return is_valid;
	}

/* NAMESPACE */
}
