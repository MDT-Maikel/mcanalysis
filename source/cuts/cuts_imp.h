/* Cut & Count classes
 *
 * Provides a class which cuts & counts in a list of events.
 * This class is based on a simple virtual cut class which 
 * can be overloaded by specialised cuts.
*/


/* NAMESPACE */
namespace analysis
{
	
	/* cut: pt */
	class cut_pt : public cut
	{
	
	public:
		cut_pt(double pt, int t, unsigned int n, double eta) : pt_cut(pt), type(t), number(n), eta_max(eta) {};
	
		bool operator() (const event *ev) 
		{ 
			const particle *p = ev->get(type, number, eta_max); 
			if (!p)
				return false;
			return p->pt() > pt_cut;
		};

	private:
		double pt_cut;
	
		int type;
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
		cut_ht(double ht, int t, double pt, double eta) : ht_cut(ht), type(t), min_pt(pt), max_eta(eta) {};
		
		bool operator() (const event *ev) 
		{ 
			return ev->ht(type, min_pt, max_eta) > ht_cut;
		};
		
	private:
		double ht_cut;
	
		int type;
		double min_pt;
		double max_eta;
		
	};
	
	/* cut: veto */
	class cut_veto: public cut
	{
		
	public:
		cut_veto(int t, double pt, double eta) : type(t), min_pt(pt), max_eta(eta) {};
		
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
		int type;
		double min_pt;
		double max_eta;			
		
	};
	
/* NAMESPACE */
}
