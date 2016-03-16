#include "IndexPdf.hh"

EXEC_TARGET fptype device_Index (fptype* evt, fptype* p, unsigned int* indices) {
  // Structure : nP index 1 inde 2 ... index n functionIndex1 parameterIndex1 ....  functionIndex1 parameterIndex1
  // nP = 3*no.of.indeces
  // Find mapping between event variables and function to evaluate
  fptype x = evt[indices[2 + indices[0]]];
  unsigned int compareIndex;

  for (size_t i = 1; i <= indeces[0]/3; i++) {
    /* code */
    if(p[indeces[i]] == x){
      compareIndex = i;
      break;
    }
  }
  //unsigned int mapFunction = (int) FLOOR(0.5 + x
    // This is an index into the IndexPdf's list of functions
  //int targetFunction = (int) FLOOR(0.5 + (*(reinterpret_cast<device_function_ptr>(device_function_table[mapFunction])))(evt, p, paramIndices + indices[2]));
  int targetFunction = compareIndex;

  targetFunction *= 2; // Because there are two pieces of information about each function
  targetFunction += indeces[0]; // Because first function information begins at index 3

  //fptype ret = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[targetFunction]])))(evt, p, paramIndices + indices[targetFunction + 1]);
  fptype ret = callFunction(evt+compareIndex, indices[targetFunction], indices[targetFunction + 1]);
  ret *= normalisationFactors[indices[targetFunction + 1]];
  //if (gpuDebug & 1)
  //if ((gpuDebug & 1) && (0 == BLOCKIDX) && (0 == THREADIDX))
  //printf("[%i, %i] Mapped: %i (%f %f %f %f) %f\n", BLOCKIDX, THREADIDX, targetFunction, evt[0], evt[1], evt[2], evt[3], ret);
  return ret;
}

MEM_DEVICE device_function_ptr ptr_to_Mapped = device_Mapped;

__host__ IndexPdf::IndexPdf (std::string n,vector<Variable*>& b, vector<GooPdf*>& t)
  : GooPdf(0, n)
{
  //components.push_back(m);
  std::vector<unsigned int> pindices;
  //pindices.push_back(m->getFunctionIndex());
  //pindices.push_back(m->getParameterIndex());
  if (b.size() != t.size()) abortWithCudaPrintFlush(__FILE__, __LINE__, getName() + " Index Size =! Pdf Size", this);

  for (vector<Variable*>::iterator v = b.begin(); v != b.end(); ++v) {
    //components.push_back(*f);
    pindices.push_back(registerParameter(v));
  }

  pindices.push_back(registerParameter(mean))
  std::set<int> functionIndicesUsed;
  for (vector<GooPdf*>::iterator f = t.begin(); f != t.end(); ++f) {
    components.push_back(*f);
    pindices.push_back((*f)->getFunctionIndex());
    pindices.push_back((*f)->getParameterIndex());
    functionIndicesUsed.insert((*f)->getFunctionIndex());
  }
  if (functionIndicesUsed.size() > 1) {
    std::cout << "Warning: More than one function type given to IndexPdf "
	      << getName()
	      << " constructor. This may slow execution by causing sequential evaluations.\n";
  }

  getObservables(observables);
  GET_FUNCTION_ADDR(ptr_to_Mapped);
  initialise(pindices);
}

__host__ fptype IndexPdf::normalise () const {
  //std::cout << "Normalising IndexPdf " << getName() << std::endl;
  fptype ret = 0;
  for (unsigned int i = 1; i < components.size(); ++i) { // No need to normalise mapping function.
    fptype curr = components[i]->normalise();
    ret += curr;
  }
  host_normalisation[parameters] = 1.0;
  return ret;
}
