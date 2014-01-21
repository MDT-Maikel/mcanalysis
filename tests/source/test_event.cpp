/* Event Tests
 *
 * Test event class as well as the different particle classes it has been build upon.
 * 
*/

#include <iostream> 
#include <cmath>
#include <ctime>
#include <random>
#include <vector> 

#include "../../source/event/event.h"
#include "../../source/particle/particle.h"
#include "../../source/particle/lhco.h"
#include "../../source/particle/lhe.h"

using namespace std;
using namespace analysis;


// main program
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old;
	clock_old = clock();
	double duration;
	
	// set test precision
	double test_precision = 0.0001;
	
	// initiate random device
	std::random_device rd;
	
	// create an event of lhco particles
	event ev_lhco;
	for (unsigned int i = 0; i < 10000; i++)
	{
		double eta = normal_distribution<double>(0.0, 1.5)(rd);
		double phi = uniform_real_distribution<double>(-3.14, 3.14)(rd);
		double pt = uniform_real_distribution<double>(0.0, 1000.0)(rd);
		double mass = uniform_real_distribution<double>(0.0, 10.0)(rd);
		lhco* p = new lhco(0, eta, phi, pt, mass);
		ev_lhco.push_back(p);
	}

	// convert into an event of lhe particles
	event ev_lhe;
	for (unsigned int i = 0; i < ev_lhco.size(); i++)
	{
		double px = ev_lhco[i]->px();
		double py = ev_lhco[i]->py();
		double pz = ev_lhco[i]->pz();
		double pe = ev_lhco[i]->pe();
		double mass = ev_lhco[i]->mass();
		lhe* p = new lhe(px, py, pz, pe, mass);
		ev_lhe.push_back(p);
	}
	
	// test for each particle if kinematics are equivalent between lhco and lhe
	// test: pt, eta, phi and mass
	bool test_lhco_lhe_passed = true;
	for (unsigned int i = 0; i < ev_lhco.size(); i++)
	{
		// test pt
		if (!(fabs(ev_lhco[i]->pt() - ev_lhe[i]->pt()) < test_precision))
		{
			cout << std::setprecision(3) << "pt: " << ev_lhco[i]->pt() << " != " << ev_lhe[i]->pt() << endl;
			test_lhco_lhe_passed = false;
			break;			
		}
		// test eta
		if (!(fabs(ev_lhco[i]->eta() - ev_lhe[i]->eta()) < test_precision))
		{
			cout << std::setprecision(3) << "eta: " << ev_lhco[i]->eta() << " != " << ev_lhe[i]->eta() << endl;
			test_lhco_lhe_passed = false;
			break;			
		}
		// test phi
		if (!(fabs(ev_lhco[i]->phi() - ev_lhe[i]->phi()) < test_precision))
		{
			cout << std::setprecision(3) << "phi: " << ev_lhco[i]->phi() << " != " << ev_lhe[i]->phi() << endl;
			cout << std::setprecision(3) << "phi diff: " << ev_lhco[i]->phi() - ev_lhe[i]->phi() << endl;
			test_lhco_lhe_passed = false;
			break;			
		}		
		// test lhe lorentz invariance
		double px = ev_lhe[i]->px();
		double py = ev_lhe[i]->py();
		double pz = ev_lhe[i]->pz();
		double pe = ev_lhe[i]->pe();
		double mass = ev_lhe[i]->mass();		
		if (!(fabs(pe * pe - px * px - py * py - pz * pz - mass * mass) < test_precision))
		{
			cout << "lhe lorentz invariance failed" << endl;
			cout << std::setprecision(3) << "E^2 - p^2 != m^2: " << pe * pe - px * px - py * py - pz * pz << " != " << mass * mass << endl;
			test_lhco_lhe_passed = false;
			break;			
		}		
	}
	
	// test event functions between both formats
	bool test_event_passed = true;
	// test the event.mass() function
	if (!fabs(ev_lhco.mass() - ev_lhe.mass()) < test_precision)
	{
		cout << std::setprecision(3) << "event.mass(): " << ev_lhco.mass() << " != " << ev_lhe.mass() << endl;
		test_event_passed = false;
	}
	// test the event.ht() function
	if (!fabs(ev_lhco.ht(0, 5.0, 5.0) - ev_lhe.ht(0, 5.0, 5.0)) < test_precision)
	{
		cout << std::setprecision(3) << "event.ht(0, 5.0, 5.0): " << ev_lhco.ht(0, 5.0, 5.0) << " != " << ev_lhe.ht(0, 5.0, 5.0) << endl;
		test_event_passed = false;
	}
		
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	clock_old = clock();
	cout << "=====================================================================" << endl;
	cout << "Event & particle test: completed in " << duration << " seconds." << endl;
	cout << "Kinematics checks between lhco and lhe classes have ";
	if (test_lhco_lhe_passed)
		cout << "passed!" << endl;
	else
		cout << "failed!" << endl;
	cout << "Event function checks between lhco and lhe classes have ";
	if (test_event_passed)
		cout << "passed!" << endl;
	else
		cout << "failed!" << endl;
	cout << "=====================================================================" << endl;
}
