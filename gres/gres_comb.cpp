/* GRES Combinatorics
 *
 * Command line program that calculates the combinatorics of any generalised search with Njet and Rmax.
 *
 * ./gres_comb <Njet> <Rmax>
 *
*/

#include <sstream>
#include <iostream>
#include <vector>

using namespace std;


int main(int argc, char* argv[])
{
	// Check for two arguments.
	if (argc != 3) 
	{
		cout << "Specify two arguments: <Njet> <Rmax>" << endl;
		return 0;
  	}

	// Get Njet and Rmax.
	string njet_in = argv[1];
	string rmax_in = argv[2];
	long njet, rmax;
	istringstream(njet_in) >> njet;
	istringstream(rmax_in) >> rmax;

	// Make sure Njet and Rmax are chosen reasonably.
	njet = max((long) 2, njet);
	rmax = max((long) 2, rmax); rmax = min(njet, rmax);
	cout << "Calculating number of combinations for Njet = " << njet << " and Rmax = " << rmax << "." << endl;

	// Initialize combinations vector with standard combination.	
	vector< vector<int> > combinations = {{2}};
	
	// Continue to add combinations until you run out of new ones.
	for (unsigned int size = 3; size <= njet; size++)
	{
		vector< vector<int> > new_combinations;
		// Build new combinations on the last cycle of new combinations.
		for (unsigned int i = 0; i < combinations.size(); i++)
		{
			vector<int> comb(combinations[i]);
			
			// No new combination based on this combination if already contains Njet jets.
			int comb_size = 0;
			for (unsigned int j = 0; j < comb.size(); j++)
				comb_size += comb[j];
			if (comb_size >= njet)
				continue;

			// Make a new combination by adding a trailing one if first element equals 2 and every other element equals 1.
			bool is_oneortwo = (comb[0] == 2);
			for (unsigned int j = 1; j < comb.size(); j++)
				if (comb[j] > 1)
					is_oneortwo = false;
			if (is_oneortwo)
			{
				vector<int> new_comb(comb);
				new_comb.push_back(1);
				new_combinations.push_back(new_comb);
			}
			
			// Make new combinations by adding 1 to each of its elements.
			// But only add one if it does not exceed the element to the left.
			// And only add one if the element to the right is one.
			for (unsigned int j = 0; j < comb.size(); j++)
			{
				if (comb[j] >= rmax)
					continue;
				if (j != 0 && comb[j] >= comb[j - 1])
					continue;
				if (j != comb.size() - 1 && comb[j + 1] != 1)
					continue;
				vector<int> new_comb(comb);
				new_comb[j]++;
				new_combinations.push_back(new_comb);
			}			
		}

		// Stop looking if no new combinations have been found.
		if (new_combinations.size() == 0)
			break;
	
		// Replace the existing combinations with the new combinations.
		unsigned int csize = combinations.size();
		for (unsigned int i = 0; i < new_combinations.size(); i++)
		{
			if (i < csize)
				combinations[i] = new_combinations[i];
			else
				combinations.push_back(new_combinations[i]);
		}
	}

	// Remove all double combinations.
	cout << "Removing all double combinations, currently " << combinations.size() << " combinations." << endl;
	for (unsigned int i = 0; i < combinations.size() - 1; i++)
	{
		vector<int> comb(combinations[i]);
		for (unsigned int j = combinations.size() - 1; j > i; j--)
		{
			// Remove every combination which is equal to comb.
			vector<int> compare(combinations[j]);
			if (comb.size() != compare.size())
				continue;
			bool is_equal = true;
			for (unsigned int k = 0; k < comb.size(); k++)
			{
				if (comb[k] != compare[k])
				{
					is_equal = false;
					break;
				}
			}
			if (is_equal)
				combinations.erase(combinations.begin() + j);
		}
	}

	// Print the result of the calculations.
	int total_combinations = combinations.size();
	cout << "The number combinations for Njet = " << njet << " and Rmax = " << rmax << " equals " << total_combinations << "." << endl;
	/*cout << "These combinations are:" << endl;
	for (int i = 0; i < total_combinations; i++)
	{ 
		vector<int> comb(combinations[i]);
		cout << "{";
		for (int j = 0; j < comb.size(); j++)
		{
			cout << comb[j];
			if (j < comb.size() - 1)
				cout << ", ";
		}
		cout << "}; ";	
	}
	cout << endl;*/
}


