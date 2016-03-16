#ifndef INDEX_PDF_HH
#define INDEX_PDF_HH

#include "GooPdf.hh" 

class IndexPdf : public GooPdf {
public:
  IndexPdf (std::string n,vector<Variable*>& b, vector<GooPdf*>& t);
  // Map function m must be custom written to correspond to order of function list t. 
  __host__ fptype normalise () const;
private:

};

#endif
