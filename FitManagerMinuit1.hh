#ifndef FITMANAGER_MINUIT1_HH
#define FITMANAGER_MINUIT1_HH

#include "TMinuit.hh"
extern PdfBase* pdfPointer;
extern int numPars;

void FitFun(int &npar, double *gin, double &fun, double *fp, int iflag);

class FitManager {
public:
  FitManager (PdfBase* dat);
  FitManager (PdfBase* dat,bool hesse);
  ~FitManager ();
  void setMaxCalls (double mxc) {overrideCallLimit = mxc;}
  void setupMinuit ();
  void runMigrad ();
  void fit ();
  TMinuit* getMinuitObject () {return minuit;}
  void getMinuitValues () const;
  TMinuit* minuit;
  void setHesse(bool hesse){runhesse = true}
private:
  double overrideCallLimit;
  bool runhesse;
};

#endif
