#ifndef EXPGAUS_PDF_HH
#define EXPGAUS_PDF_HH

#include "GooPdf.hh" 

class ExpGausPdf : public GooPdf {
public:
  ExpGausPdf (std::string n, Variable* _x, Variable* m, Variable* s, Variable* t); 

private:

};

#endif
