/* BumpHunter class
 *
 * implementation of the BumpHunter algorithm based on: 
 * 	C++code: https://twiki.cern.ch/twiki/bin/view/Main/BumpHunter by W. Buttinger (2012)
 * 	article: http://arxiv.org/pdf/1101.0390v2.pdf by G. Choudalakis (2011)
 * 
*/

#ifndef INC_BUMPHUNTER
#define INC_BUMPHUNTER

#include <algorithm>
#include <cfloat>
#include <iostream>
#include <sstream>
#include <vector>

#include <TH1D.h>
#include <TMath.h>
#include <Math/SpecFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>
#include <TRandom3.h>
#include <Math/GSLRndmEngines.h>
#include <Math/Random.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TList.h>
#include <TLatex.h>
#include <TLine.h>
#include <Math/QuantFuncMathCore.h>


/* NAMESPACE */
namespace analysis
{

	class bumphunter : public TObject 
	{
		
	public:
	
		/* hunt properties: data types */
		enum TestStatisticType {BUMPHUNTER, DIPHUNTER, DISCREPANCYHUNTER, LOGLIKELIHOOD};
		enum bin_model {BUMP_POISSON_GAMMA = 0, BUMP_POISSON = 1, BUMP_GAUSSIAN = 2, BUMP_POISSON_GAUSSIAN = 3}; 

		/* channel: background, data, pseudodata and bump window */
		struct channel 
		{ 
			TH1* hist_bkg;
			TH1* hist_data;
			TH1* hist_pseudodata;
			Double_t search_min; // start of the bump hunt
			Double_t search_max; // end of the bump hunt
			Double_t bump_window_min; // minimal bump window size
			Double_t bump_window_max; // maximal bump window size
			Double_t bump_pvalue; // pvalue of the bump
			std::vector<std::pair<Int_t, Int_t> > central_windows;
			channel() : hist_bkg(0), hist_data(0), hist_pseudodata(0), search_min(0./0.), search_max(0./0.), bump_window_min(0), bump_window_max(0), bump_pvalue(1) {}
		};

	public:
	
		/* con & destructor */
		bumphunter(TH1* mc, TH1* data);
		bumphunter();
		~bumphunter();
		
	private:
	
		/* initialization */
		void init();
	
	public:
	
		/* general settings */
		void set_name(std::string n);
		void set_folder(std::string f);
	
		/* adding data */
		void SetChannelDistribution(TH1* hist, TH1* data);
		void AddChannelDistribution(TH1* hist, TH1* data);
		void SetChannel(Int_t in);
		Int_t GetNChannels();
				
		/* hunt properties: general */
		void SetBinModel(Int_t model);
		void SetTestStatisticType(Int_t t);
		void SetNPseudoExperiments(Int_t n);
		void SetBumpOverlapFactor(Double_t in);
		
		/* hunt properties: search region */
		void SetSearchRegion(Double_t low, Double_t high); 
		Double_t GetSearchLowEdge();
		Double_t GetSearchHighEdge();
		
		/* hunt properties: bump window */
		void SetMinWindowSize(Int_t in);
		void SetMaxWindowSize(Int_t in);
		void SetWindowStepSize(Int_t in);
		Double_t GetMinWindowSize();
		Double_t GetMaxWindowSize();
		Double_t GetWindowSizeStep();
		
		/* hunt properties: bump details */
		Double_t GetBumpLowEdge();
		Double_t GetBumpHighEdge();
	
		/* hunt properties: results */
		double get_global_pvalue();
		double get_local_pvalue();

		/* statistics */
		static Double_t GetPoissonPValue(Double_t nobs, Double_t e);
		static Double_t GetGaussianPValue(Double_t nobs, Double_t e, Double_t sigma);
		static Double_t GetPoissonConvGammaPValue(Double_t nobs, Double_t e, Double_t err);
		static Double_t GetZValue(Double_t pvalue, bool& overflow);
		
		/* bump hunting */
		Int_t run(); // do bumphunting. 0 = success

		/* TO SORT */
		// update the current pseudo-experiment histogram with a new set of values
		TH1* GenerateToyMC();

		TH1* GeneratePseudoData(TH1* mc); 

		// evaluates the bumphunter test statistic for the given data
		Double_t EvaluateTestStatistic(TH1* data, Bool_t printOut = false);

		// evaluate for multichannel
		Double_t EvaluateMultiChannelTestStatistic(Bool_t generatePseudo, Bool_t printOut = false);

		void EvaluateSearchPattern();
		void PrintSearchPattern();
		
	private:
	
		/* general settings */
		std::string hunt_name;
		std::string hunt_folder;
		
		/* adding data */
		std::vector<channel> channel_list;
		unsigned int current_channel;
		
		/* hunt properties: general */
		Int_t fBinModel;
		Int_t fTestStatisticType;
		Int_t fnPseudo;			
		Double_t fMultiChannelBumpOverlapFactor; //how much of the bump's regions must overlap. Default is 1 (perfect overlap)
		
		/* hunt properties: bump window */
		Double_t fMinWindowSize;
		Double_t fMaxWindowSize; 
		Double_t fWindowSizeStep;

		/* hunt properties: results */
		Double_t local_pvalue; //smallest p-value found in the actual data.. i.e. the bump's local pValue
		Double_t global_pvalue;
		
		/* statistics */
		TRandom3 pRand;
		ROOT::Math::Random<ROOT::Math::GSLRngMT> r; //used for gamma distribution random sampling
		
	};

/* NAMESPACE */
}

#endif
