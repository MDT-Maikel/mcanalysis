/* MT2 Bisect
 *
 * Based on an implementation by:
 * 		Hsin-Chia Cheng, Zhenyu Han    
 *		Dec 11, 2008, v1.01a 
 *
*/

/*******************************************************************************
  Usage: 

  1. Define an object of type "mt2":
     
     mt2_bisect::mt2 mt2_event;
 
  2. Set momenta and the mass of the invisible particle, mn:
 
     mt2_event.set_momenta( pa, pb, pmiss );
     mt2_event.set_mass( mn );
 
     where array pa[0..2], pb[0..2], pmiss[0..2] contains (mass,px,py) 
     for the visible particles and the missing momentum. pmiss[0] is not used. 
     All quantities are given in double.    

  3. Use mt2::get_mt2() to obtain the value of mt2:

     double mt2_value = mt2_event.get_mt2();       
          
*******************************************************************************/ 

#ifndef INC_MT2_BISECT
#define INC_MT2_BISECT

#include <iostream>
#include <math.h>


/* NAMESPACE */
namespace analysis
{
namespace mt2_bisect
{

	/*The user can change the desired precision below, the larger one of the following two definitions is used. Relative precision less than 0.00001 is not guaranteed to be achievable--use with caution*/ 

	#define RELATIVE_PRECISION 0.00001 //defined as precision = RELATIVE_PRECISION * scale, where scale = max{Ea, Eb}
	#define ABSOLUTE_PRECISION 0.0     //absolute precision for mt2, unused by default

	//Reserved for expert
	#define MIN_MASS  0.1   //if ma<MINMASS and mb<MINMASS, use massless code
	#define ZERO_MASS 0.000 //give massless particles a small mass
	#define SCANSTEP 0.1

	class mt2
	{  
	public:

		mt2();
		void   mt2_bisect();
		void   mt2_massless();
		void   set_momenta(double *pa0, double *pb0, double* pmiss0);
		void   set_mn(double mn);
		double get_mt2();
		void   print();
		int    nevt;

	private:  

		bool   solved;
		bool   momenta_set;
		double mt2_b;

		int nsols(double Dsq);
		int nsols_massless(double Dsq);
		inline int signchange_n(long double t1, long double t2, long double t3, long double t4, long double t5);
		inline int signchange_p(long double t1, long double t2, long double t3, long double t4, long double t5);
		int scan_high(double &Deltasq_high);
		int find_high(double &Deltasq_high);
		
		//data members
		double pax, pay, ma, Ea;
		double pmissx, pmissy;
		double pbx, pby, mb, Eb;
		double mn, mn_unscale;

		//auxiliary definitions
		double masq, Easq;
		double mbsq, Ebsq;
		double pmissxsq, pmissysq;
		double mnsq;

		//auxiliary coefficients
		double a1, b1, c1, a2, b2, c2, d1, e1, f1, d2, e2, f2;
		double d11, e11, f12, f10, d21, d20, e21, e20, f22, f21, f20;

		double scale;
		double precision;
		  
	};

/* NAMESPACE */
}
}

#endif
