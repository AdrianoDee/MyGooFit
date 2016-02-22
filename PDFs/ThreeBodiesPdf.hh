#ifndef THREE_PDF_HH
#define THREE_PDF_HH

#include "GooPdf.hh" 

class ThreePdf : public GooPdf {
public:
  
    ThreePdf (std::string n, Variable* _x);//,Variable* nBk);
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return false;} 



private:

};

#endif

