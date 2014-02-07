/* Plot class
 *
 *
*/


/* NAMESPACE */
namespace analysis
{

	/* plot: pt */
	class plot_pt : public plot_default
	{
	public:
		plot_pt(unsigned int t, unsigned int n, double eta = 5.0) : type(t), number(n), eta_max(eta) {}

		double operator() (const event *ev) 
		{ 
			particle *p = ev->get(type, number, eta_max);
			if (p)
				return p->pt();
			return 0;
		}

	private:
		unsigned int type;
		unsigned int number;
		double eta_max;
	};
	
	/* plot: eta */
	class plot_eta : public plot_default
	{
	public:
		plot_eta(unsigned int t, unsigned int n) : type(t), number(n) {}

		double operator() (const event *ev) 
		{ 
			particle *p = ev->get(type, number);
			if (p)
				return p->eta();
			return 0;
		}

	private:
		unsigned int type;
		unsigned int number;
	};
	
	/* plot: phi */
	class plot_phi : public plot_default
	{
	public:
		plot_phi(unsigned int t, unsigned int n) : type(t), number(n) {}

		double operator() (const event *ev) 
		{ 
			particle *p = ev->get(type, number);
			if (p)
				return p->phi();
			return 0;
		}

	private:
		unsigned int type;
		unsigned int number;
	};
	
	/* cut: delta pt */
	class plot_deltapt : public plot_default
	{
	public:
		plot_deltapt(unsigned int t1, unsigned int n1, unsigned int t2, unsigned int n2, double eta = 5.0) :
			type1(t1), number1(n1), type2(t2), number2(n2), eta_max(eta) {}

		double operator() (const event *ev) 
		{ 
			particle *p1 = ev->get(type1, number1, eta_max);
			particle *p2 = ev->get(type2, number2, eta_max);
			if (p1 && p2)
				return p1->pt() - p2->pt();
			return 0;
		}

	private:
		unsigned int type1;
		unsigned int type2;
		unsigned int number1;
		unsigned int number2;
		double eta_max;
	};
	
	/* plot: met */
	class plot_met : public plot_default
	{
	public:
		plot_met() {}

		double operator() (const event *ev) 
		{ 
			return ev->met();
		}
	};
	
	/* plot: ht */
	class plot_ht : public plot_default
	{
	public:
		plot_ht(unsigned int t, double pt = 0.0, double eta = 5.0) : type(t), min_pt(pt), max_eta(eta) {};
		
		double operator() (const event *ev) 
		{ 
			return ev->ht(type, min_pt, max_eta);
		};

	private:
		unsigned int type;
		double min_pt;
		double max_eta;
	};

	/* plot: mass */
	class plot_mass : public plot_default
	{
	public:
		plot_mass(int t, const std::vector<int> & c) : type(t), comb(c) {}

		double operator() (const event *ev) 
		{ 
			ev->mass(type, comb);
		}

	private:
		int type;
		std::vector<int> comb;
	};

	/* plot: mt2 */
	class plot_mt2 : public plot_default
	{
	public:
		plot_mt2() {}

		double operator() (const event *ev) 
		{ 
			return ev->mt2();
		}
	};

/* NAMESPACE */
}
