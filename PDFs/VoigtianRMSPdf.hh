#ifndef VOIGTIANRMS_PDF_HH
#define VOIGTIANRMS_PDF_HH

#include "GooPdf.hh" 

class VoigtianRMSPdf : public GooPdf {
public:
  VoigtianRMSPdf (std::string n, Variable* _x, Variable* m, Variable* w);

private:

};

#endif
