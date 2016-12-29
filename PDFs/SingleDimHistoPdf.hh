#ifndef SINGLEDIM_PDF_HH
#define SINGLEDIM_PDF_HH

#include "GooPdf.hh"
#include "SingleDimHistoPdf.hh"

class SingleDimHistoPdf : public GooPdf {
public:
  SingleDimHistoPdf (std::string n,
			  BinnedDataSet* x,
        std::vector<Variable*> obses,unsigned int interOrder);
  //__host__ virtual fptype normalise () const;

private:
  thrust::device_vector<fptype>* dev_base_histogram;
  fptype totalEvents;
  fptype* host_constants;
  int numVars;
};

#endif
