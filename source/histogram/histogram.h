#ifndef __ANALYSIS_histogram_H__
#define __ANALYSIS_histogram_H__ 1

//  std includes
#include <iostream> // std input-output stream
#include <vector> // vector type
#include <string> // string type
#include <sstream> // stringstream type
#include <fstream> // ifstream type
#include <cmath> // math library

//  root includes
#include <TCanvas.h> //basic canvas class
#include <TStyle.h> // style class
#include <TLegend.h> // draw legend class
#include <TH1.h> // 1D histogram class
// #include <TF1.h> // draw 1D functions class
// #include <TH2.h> // 2D histogram class
// #include <TGraph.h> // graphics class
// #include <TNtuple.h> // N-tuple class
// #include <TVector3.h> // three-vector class
// #include <TLorentzVector.h> // Lorentz vector class
// #include <TRandom.h> // random number generation class
// #include <TRandom3.h> // extended random number generation class
// #include <TROOT.h> // entry point to the ROOT system
// #include <TChain.h> // access to collection of files containg TTree objects
// #include <TFile.h> // a ROOT file is a suite of consecutive data records (TKey's)
// #include <TTree.h> // a TTree object consists of a list of independent branches (TBranch)
// #include <TPaveText.h> // a PaveText is a Pave (see TPave) with text, lines or/and boxes inside.
// #include <Rtypes.h> // basic types used by ROOT

using namespace std;

namespace analysis {

  class histogram
  {
  private:
    TCanvas* c1;
    vector<string> samples;
    double nbins;
    double hmin; 
    double hmax;
    string ps_title;
    string leg_title;
    string x_label;
    string y_label;

  public:
    //=== Class constructor ===//
    histogram();

    //=== Class destructor ===//
    // ~histogram() = default;

    //=== Set histogram Options ===//
    void set_ps_title(string ps);
    void set_hist_bins(double nb);
    void set_hist_range(double min, double max);
    void set_x_label(string x);
    void set_y_label(string y);
    void set_leg_title(string title);
    void add_sample(string file_name);
    
    //=== Draw histogram ===//
    void draw();
  };

}

#endif