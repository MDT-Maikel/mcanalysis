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
