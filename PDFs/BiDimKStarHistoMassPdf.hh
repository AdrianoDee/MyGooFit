#ifndef BIDIMKSTARMASS_PDF_HH
#define BIDIMKSTARMASS_PDF_HH

#include "GooPdf.hh"
#include "BinnedDataSet.hh"

class BiDIm : public GooPdf {
public:
  BiDimKStarHistoMassPdf (std::string n,
			  BinnedDataSet* x,
			  std::vector<Variable*> obses);
  //__host__ virtual fptype normalise () const;

private:
  thrust::device_vector<fptype>* dev_base_histogram;
  fptype totalEvents;
  fptype* host_constants;
  int numVars;
};

#endif
