/* Event Tests
 *
 * Test event class as well as the different particle classes it has been build upon.
 * 
*/

#include <cmath>
#include <ctime>
#include <iostream> 
#include <random>
#include <vector> 

#include "event/event.h"
#include "particle/lhco.h"
#include "particle/lhe.h"
#include "particle/particle.h"

using namespace std;
using namespace analysis;


// main program
int main(int argc, const char* argv[])
{
	// initiate timing procedure
	clock_t clock_old = clock();
	double duration;
		
	// set test precision
	double test_precision = 0.00001;
	
	// initiate random device
	std::random_device rd;
	
	// create an event of lhco particles
	event *ev_lhco = new event;
	for (unsigned int i = 0; i < 1000; i++)
	{
		double eta = normal_distribution<double>(0.0, 1.5)(rd);
		double phi = uniform_real_distribution<double>(-3.14, 3.14)(rd);
		double pt = uniform_real_distribution<double>(0.0, 1000.0)(rd);
		double mass = uniform_real_distribution<double>(0.0, 10.0)(rd);
		lhco* p = new lhco(ptype_none, eta, phi, pt, mass);
		ev_lhco->push_back(p);
	}

	// convert into an event of lhe particles
	event *ev_lhe = new event;
	for (unsigned int i = 0; i < ev_lhco->size(); i++)
	{
		particle *p_lhco = (*ev_lhco)[i];
		double px = p_lhco->px();
		double py = p_lhco->py();
		double pz = p_lhco->pz();
		double pe = p_lhco->pe();
		double mass = p_lhco->mass();
		lhe* p = new lhe(px, py, pz, pe, mass);
		p->set_final(true);
		ev_lhe->push_back(p);
	}
	
	// test for each particle if kinematics are equivalent between lhco and lhe
	// test: pt, eta, phi and mass
	bool test_lhco_lhe_passed = true;
	for (unsigned int i = 0; i < ev_lhco->size(); i++)
	{
		particle *p_lhco = (*ev_lhco)[i];
		particle *p_lhe = (*ev_lhe)[i];		
		// test pt
		if (!(fabs((p_lhco->pt() - p_lhe->pt()) / p_lhco->pt()) < test_precision))
		{
			cout << std::setprecision(3) << "pt: " << p_lhco->pt() << " != " << p_lhe->pt() << " " << i << endl;
			test_lhco_lhe_passed = false;
			//break;			
		}
		// test eta
		if (!(fabs((p_lhco->eta() - p_lhe->eta()) / p_lhco->eta()) < test_precision))
		{
			cout << std::setprecision(3) << "eta: " << p_lhco->eta() << " != " << p_lhe->eta() << endl;
			test_lhco_lhe_passed = false;
			//break;			
		}
		// test phi
		if (!(fabs((p_lhco->phi() - p_lhe->phi()) / p_lhco->phi()) < test_precision))
		{
			cout << std::setprecision(3) << "phi: " << p_lhco->phi() << " != " << p_lhe->phi() << endl;
			cout << std::setprecision(3) << "phi diff: " << p_lhco->phi() - p_lhe->phi() << endl;
			test_lhco_lhe_passed = false;
			//break;			
		}		
		// test lhe lorentz invariance
		double px = p_lhe->px();
		double py = p_lhe->py();
		double pz = p_lhe->pz();
		double pe = p_lhe->pe();
		double mass = p_lhe->mass();		
		if (!(fabs((pe * pe - px * px - py * py - pz * pz - mass * mass) / (pe * pe)) < test_precision))
		{
			cout << "lhe lorentz invariance failed" << endl;
			cout << std::setprecision(3) << "E^2 - p^2 != m^2: " << pe * pe - px * px - py * py - pz * pz << " != " << mass * mass << endl;
			test_lhco_lhe_passed = false;
			//break;			
		}		
	}
	
	// test event functions between both formats
	bool test_event_passed = true;
	// test the event.mass() function
	if (!fabs((ev_lhco->mass() - ev_lhe->mass()) / ev_lhco->mass()) < test_precision)
	{
		cout << std::setprecision(6) << "event.mass(): " << ev_lhco->mass() << " != " << ev_lhe->mass() << endl;
		test_event_passed = false;
	}
	// test the event.ht() function
	if (fabs((ev_lhco->ht(ptype_all, 5.0, 5.0) - ev_lhe->ht(ptype_all, 5.0, 5.0)) / ev_lhco->ht(ptype_all, 5.0, 5.0)) > test_precision)
	{
		cout << std::setprecision(6) << "event.ht(-1, 5.0, 5.0): " << ev_lhco->ht(ptype_all, 5.0, 5.0) << " != " << ev_lhe->ht(ptype_all, 5.0, 5.0) << endl;
		test_event_passed = false;
	}
		
	// log results
	duration = (clock() - clock_old) / static_cast<double>(CLOCKS_PER_SEC);
	cout << "=====================================================================" << endl;
	cout << "Event & particle test: completed in " << duration << " seconds." << endl;
	cout << "Kinematics checks between lhco and lhe classes have " << (test_lhco_lhe_passed ? "passed!" : "failed!") << endl;
	cout << "Event function checks between lhco and lhe classes have " << (test_event_passed ? "passed!" : "failed!") << endl;
	cout << "=====================================================================" << endl;
	
	// clear remaining event pointers
	delete ev_lhco;
	delete ev_lhe;
	
	// return whether tests passed
	if (test_lhco_lhe_passed && test_event_passed)
		return EXIT_SUCCESS;
	return EXIT_FAILURE;
}
