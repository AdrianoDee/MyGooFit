#ifndef ARGUS_PDF_HH
#define ARGUS_PDF_HH

#include "GooPdf.hh" 

class ArgusPdf : public GooPdf {
public:
  ArgusPdf (std::string n, Variable* _x, Variable* m, Variable* s, bool upper, Variable* power = 0); 
  __host__ virtual bool hasAnalyticIntegral () const {return false;} 
  __host__ fptype integrate (fptype lo, fptype hi) const; 

private:

};

#endif
