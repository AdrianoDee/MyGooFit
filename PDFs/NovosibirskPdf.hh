#ifndef NOVOSIBIRSK_PDF_HH
#define NOVOSIBIRSK_PDF_HH

#include "GooPdf.hh" 

class NovosibirskPdf : public GooPdf {
public:
  NovosibirskPdf (std::string n, Variable* _x, Variable* m, Variable* s, Variable* t); 

private:

};

#endif
