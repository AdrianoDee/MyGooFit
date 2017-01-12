#ifndef BIDIMHIST_PDF_HH
#define BIDIMHIST_PDF_HH

#include "GooPdf.hh"
#include "BinnedDataSet.hh"

class BiDimHistoPdf : public GooPdf {
public:
  BiDimHistoPdf (std::string n,
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
