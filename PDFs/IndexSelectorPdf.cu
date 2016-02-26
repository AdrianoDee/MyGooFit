#include "IndexSelectorPdf.hh"

__device__ fptype device_IndexSelector (fptype* evt, fptype* p, unsigned int* indices)
{
  fptype x = evt[indices[2 + indices[0]]];
  unsigned int index_int = FLOOR(0.5 + x);
  unsigned int ret = indices[1];
  // indices[1] should be the number of step points we have
  for(unsigned int i = 0; i < indices[1]; i++)
  {
    if(index_int < p[indices[2 + i]])
    {
      ret = i;
      break;
    }
  }

  //printf("Staircase returning %u based on %u\n", ret, x_int);
  return ret;
}

__device__ device_function_ptr ptr_to_Index = device_IndexSelector;
device_function_ptr ptr_to_Index = device_IndexSelector;

__host__ IndexSelectorPdf::IndexSelectorPdf(std::string n, Variable* _index, const std::vector<Variable*> &fIndexList)
  : GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(fIndexList.size());
  for(std::vector<Variable*>::const_iterator index0 = fIndexList.begin(); index0 != fIndexList.end(); x0++)
    pindices.push_back(registerParameter(*index0));

  GET_FUNCTION_ADDR(ptr_to_Index);
  initialise(pindices);
  std::cout << "IndexSelectorPdf::IndexSelectorPdf(" << n << ", ...)" << std::endl;
}

//__host__ fptype StaircasePdf::integrate (fptype lo, fptype hi) const {
//  unsigned int* indices = host_indices+parameters;

  // the part where the function is zero is

  //fptype x0 = host_params[indices[1]];
  //return (hi - x0);
//}