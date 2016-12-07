#ifndef THRUST_PDF_FUNCTOR_HH
#define THRUST_PDF_FUNCTOR_HH
#include <stdlib.h>
#include <math.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/iterator/zip_iterator.h>
#include "thrust/iterator/constant_iterator.h"
#include <thrust/transform_reduce.h>
#include <thrust/transform.h>
#include <thrust/fill.h>
#include <thrust/host_vector.h>
#include <cmath>
#include <cassert>
#include <set>

#include "PdfBase.hh"

EXEC_TARGET int dev_powi (int base, int exp); // Implemented in SmoothHistogramPdf.

#define CALLS_TO_PRINT 10
typedef fptype (*device_function_ptr) (fptype*, fptype*, unsigned int*);            // Pass event, parameters, index into parameters.
typedef fptype (*device_metric_ptr) (fptype, fptype*, unsigned int);

extern void* host_fcn_ptr;
extern void* host_int_ptr;
extern void* host_pnt_ptr;

class MetricTaker;

class GooPdf : public PdfBase {
public:

  GooPdf (Variable* x, std::string n);
  __host__ virtual double calculateNLL () const;
  __host__ void evaluateAtPoints (std::vector<fptype>& points) const;
  __host__ void evaluateAtPoints (Variable* var, std::vector<fptype>& res);
  __host__ void evaluateAtPointsInt (Variable* var, std::vector<fptype>& res);
  __host__ virtual fptype normalise () const;
  __host__ virtual fptype normaliseInt () const;
  __host__ virtual fptype integrate (fptype lo, fptype hi) const {return 0;}
  __host__ virtual fptype doIntegral ();
  __host__ virtual bool hasAnalyticIntegral () const {return false;}
  __host__ fptype getValue ();
  __host__ void getCompProbsAtDataPoints (std::vector<std::vector<fptype> >& values);
  __host__ void getCompProbsAtDataPoints (std::vector<std::vector<fptype> >& values,unsigned int evts);
  __host__ void getCompProbsAtDataPointsUn (std::vector<std::vector<fptype> >& values);
  __host__ void initialise (std::vector<unsigned int> pindices, void* dev_functionPtr = host_fcn_ptr);
  __host__ void initialiseInt (std::vector<unsigned int> pindices, void* dev_functionPtr = host_fcn_ptr,  void* dev_integralPtr = host_int_ptr, void* dev_pointPtr = host_pnt_ptr);
  __host__ void scan (Variable* var, std::vector<fptype>& values);
  __host__ virtual void setFitControl (FitControl* const fc, bool takeOwnerShip = true);
  __host__ virtual void setMetrics ();
  __host__ virtual void setMetricsInt();
  __host__ void setParameterConstantness (bool constant = true);

  __host__ virtual void transformGrid (fptype* host_output);
  static __host__ int findFunctionIdx (void* dev_functionPtr);
  __host__ void setDebugMask (int mask, bool setSpecific = true) const;

protected:
  __host__ virtual double sumOfNll (int numVars) const;
  MetricTaker* logger;
private:

};

class MetricTaker : public thrust::unary_function<thrust::tuple<int, fptype*, int>, fptype> {
public:

  MetricTaker (PdfBase* dat, void* dev_functionPtr);
  MetricTaker (PdfBase* dat, void* dev_functionPtr,unsigned int flagInt);
  MetricTaker (int fIdx, int pIdx);
  EXEC_TARGET fptype operator () (thrust::tuple<int, fptype*, int> t) const;           // Event number, dev_event_array (pass this way for nvcc reasons), event size
  EXEC_TARGET fptype operator () (thrust::tuple<int, int, fptype*> t) const;           // Event number, event size, normalisation ranges (for binned stuff, eg integration)
  EXEC_TARGET fptype operator () (thrust::tuple<fptype*,int, int> t) const;       // Bin integral

private:

  unsigned int metricIndex; // Function-pointer index of processing function, eg logarithm, chi-square, other metric.
  unsigned int functionIdx; // Function-pointer index of actual PDF
  unsigned int functionIntegralIdx; // Function-pointer index of bin integral PDF
  unsigned int functionPointIdx;
  unsigned int parameters;
  unsigned int intFlag;


};



#endif
