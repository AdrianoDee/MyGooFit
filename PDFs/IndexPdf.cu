#include "IndexPdf.hh"

EXEC_TARGET fptype device_Index (fptype* evt, fptype* p, unsigned int* indices) {
  // Structure : nP index 1 inde 2 ... index n functionIndex1 parameterIndex1 ....  functionIndex1 parameterIndex1
  // nP = 3*no.of.indeces
  // Find mapping between event variables and function to
  for (size_t i = 0; i < 20; i++) {
    printf("indices[%d] = %.3f \n",i,indices[i]);
  }
  fptype x = evt[indices[2 + indices[0]]];
  unsigned int compareIndex;
  int debug = 0;
  printf("here %d \n",debug);debug++;
  for (size_t i = 1; i <= indices[0]/3; i++) {
    /* code */
    printf("x = %.2f Parameter %d =  %.2f \n",x,i,p[indices[i]]);debug++;
    if(p[indices[i]] == x){
      compareIndex = i;
      printf("Inside -> x = %.2f Parameter %d =  %.2f \n",x,i,p[indices[i]]);debug++;
      break;
    }
  }
  printf("here %d \n",debug);debug++;
  //unsigned int mapFunction = (int) FLOOR(0.5 + x
    // This is an index into the IndexPdf's list of functions
  //int targetFunction = (int) FLOOR(0.5 + (*(reinterpret_cast<device_function_ptr>(device_function_table[mapFunction])))(evt, p, paramIndices + indices[2]));
  int targetFunction = compareIndex;
  printf("here %d \n",debug);debug++;
  targetFunction *= 2; // Because there are two pieces of information about each function
  targetFunction += indices[0]; // Because first function information begins at index 3
  printf("here %d \n",debug);debug++;
  //fptype ret = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[targetFunction]])))(evt, p, paramIndices + indices[targetFunction + 1]);
  fptype ret = callFunction(evt+compareIndex, indices[targetFunction], indices[targetFunction + 1]);
  printf("x = %.2f compareIndex = %u indices[0] = %.2f targetFunction = %d pdfX = %.2f \n",x,compareIndex,indices[0],targetFunction,evt[indices[2 + indices[targetFunction]]+compareIndex+targetFunction]);
  ret *= normalisationFactors[indices[targetFunction + 1]];
  //if (gpuDebug & 1)
  //if ((gpuDebug & 1) && (0 == BLOCKIDX) && (0 == THREADIDX))
  //printf("[%i, %i] Mapped: %i (%f %f %f %f) %f\n", BLOCKIDX, THREADIDX, targetFunction, evt[0], evt[1], evt[2], evt[3], ret);
  return ret;
}

MEM_DEVICE device_function_ptr ptr_to_Index = device_Index;

__host__ IndexPdf::IndexPdf (std::string n,vector<Variable*>& vars,vector<Variable*>& b, vector<GooPdf*>& t)
  : GooPdf(0, n)
{
  for (unsigned int i = 0; i < vars.size(); ++i) {
    registerObservable(vars[i]);
  }
  //components.push_back(m);
  std::vector<unsigned int> pindices;
  //pindices.push_back(m->getFunctionIndex());
  //pindices.push_back(m->getParameterIndex());
  if (b.size() != t.size()) abortWithCudaPrintFlush(__FILE__, __LINE__, getName() + " Index Size =! Pdf Size", this);

  for (vector<Variable*>::iterator v = b.begin(); v != b.end(); ++v) {
    //components.push_back(*f);
    pindices.push_back(registerParameter(*v));
  }

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

  //getObservables(observables);
  GET_FUNCTION_ADDR(ptr_to_Index);
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
