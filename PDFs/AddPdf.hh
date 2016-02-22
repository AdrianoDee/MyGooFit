#ifndef ADD_PDF_HH
#define ADD_PDF_HH

#include "GooPdf.hh" 

class AddPdf : public GooPdf {
public:

  AddPdf (std::string n, std::vector<Variable*> weights, std::vector<PdfBase*> comps);
  AddPdf (std::string n, std::vector<Variable*> weights, std::vector<PdfBase*> comps,unsigned int integr);
  AddPdf (std::string n, Variable* frac1, PdfBase* func1, PdfBase* func2);
  AddPdf (std::string n, Variable* frac1, PdfBase* func1, PdfBase* func2,unsigned int integr);
  __host__ virtual fptype normalise () const;
  __host__ virtual bool hasAnalyticIntegral () const {return false;}

protected:
  __host__ virtual double sumOfNll (int numVars) const;

private:
  bool extended; 
};

#endif
