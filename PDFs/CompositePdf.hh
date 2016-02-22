#ifndef COMPOSITE_PDF_HH
#define COMPOSITE_PDF_HH

#include "GooPdf.hh" 

// Composites of arbitrary functions, ie f(x) = h(g(x)) 
// for any h and g. In principle we should allow multi-
// dimensional compositing, eg f(x, y) = i(g(x, y), h(x, y)).
// Not implemented yet. 

class CompositePdf : public GooPdf {
public:
  CompositePdf (std::string n, PdfBase* core, PdfBase* shell); // Where 'core' corresponds to 'g' and 'shell' to 'h'. 
  __host__ virtual fptype normalise () const;

private:

};

#endif
