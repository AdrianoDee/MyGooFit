#ifndef VOIGTIAN_THRESH_PDF_HH
#define VOIGTIAN_THRESH_PDF_HH

#include "GooPdf.hh" 

class VoigtianThreshPdf : public GooPdf {
public:
  VoigtianThreshPdf (std::string n, Variable* _x, Variable* m, Variable* s,Variable* w);

private:

};

#endif
