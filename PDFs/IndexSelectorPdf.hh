#ifndef STAIRCASE_THRUST_FUNCTOR_HH
#define STAIRCASE_THRUST_FUNCTOR_HH

#include "GooPdf.hh" 
#include <vector>

class StaircasePdf : public GooPdf {
public:
  StaircasePdf (std::string n, Variable* _x, const std::vector<Variable*> &x0list);
  //__host__ fptype integrate (fptype lo, fptype hi) const;
  //__host__ virtual bool hasAnalyticIntegral () const {return true;}
private:
};

#endif
