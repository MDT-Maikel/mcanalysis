/* Cut & Count classes
 *
 * Provides a class which cuts & counts in a list of events.
 * This class is based on a simple virtual cut class which 
 * can be overloaded by specialised cuts.
*/


/* NAMESPACE */
namespace analysis
{
	
	/* cut: particle */
	class cut_particle: public cut
	{
	public:
		cut_particle(unsigned int t, unsigned int n) : type(t), number(n) {};
	
		bool operator() (const event *ev) 
		{ 
			const particle *p = ev->get(type, number); 
			if (!p)
				return false;
			return true;
		};
	private:
		unsigned int type;
		unsigned int number;
	};
	
	/* cut: pt */
	class cut_pt : public cut
	{	
	public:
		cut_pt(double pt, unsigned int t, unsigned int n, double eta) : pt_cut(pt), type(t), number(n), eta_max(eta) {};
	
		bool operator() (const event *ev) 
		{ 
			const particle *p = ev->get(type, number, eta_max); 
			if (!p)
				return false;
			return p->pt() > pt_cut;
		};
	private:
		double pt_cut;
		unsigned int type;
		unsigned int number;
		double eta_max;
	};
	
	/* cut: met */
	class cut_met : public cut
	{		
	public:
		cut_met(double met) : met_cut(met) {};
		
		bool operator() (const event *ev) 
		{ 
			return ev->met() > met_cut;
		};		
	private:
		double met_cut;		
	};
	
	/* cut: ht */
	class cut_ht : public cut
	{		
	public:
		cut_ht(double ht, unsigned int t, double pt, double eta) : ht_cut(ht), type(t), min_pt(pt), max_eta(eta) {};
		
		bool operator() (const event *ev) 
		{ 
			return ev->ht(type, min_pt, max_eta) > ht_cut;
		};		
	private:
		double ht_cut;
		unsigned int type;
		double min_pt;
		double max_eta;		
	};
	
	/* cut: veto */
	class cut_veto: public cut
	{		
	public:
		cut_veto(unsigned int t, double pt, double eta) : type(t), min_pt(pt), max_eta(eta) {};
		
		bool operator() (const event *ev) 
		{ 
			unsigned int index = 1;
			bool has_passed = true;
			const particle *p;
			while ((p = ev->get(type, index)))
			{
				if (p->pt() > min_pt && std::abs(p->eta()) < max_eta) 
				{
					has_passed = false;
					break;
				}
				index++;
			}
			return has_passed;
		};		
	private: 
		unsigned int type;
		double min_pt;
		double max_eta;			
	};
	
	/* cut: mt2 */
	class cut_mt2: public cut
	{		
	public:
		cut_mt2(double mt2, double mn) : mt2_cut(mt2), mt2_mn(mn) {};
		
		bool operator() (const event *ev) 
		{ 
			return ev->mt2(mt2_mn) > mt2_cut;
		};		
	private: 
		double mt2_cut;
		double mt2_mn;		
	};
	
/* NAMESPACE */
}
