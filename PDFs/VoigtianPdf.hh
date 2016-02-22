#ifndef VOIGTIAN_PDF_HH
#define VOIGTIAN_PDF_HH

#include "GooPdf.hh" 

class VoigtianPdf : public GooPdf {
public:
  VoigtianPdf (std::string n, Variable* _x, Variable* m, Variable* s,Variable* w);
  VoigtianPdf (std::string n, Variable* _x, Variable* m, Variable* s,Variable* w,Variable* o,Variable* lowerCut = 0, Variable* upperCut = 0);
private:
    Variable* offset;
    Variable* lowerThreshold;
    Variable* upperThreshold;
};

#endif
