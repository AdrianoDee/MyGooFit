#include "IndexPdf.hh"

EXEC_TARGET fptype device_Index (fptype* evt, fptype* p, unsigned int* indices) {
  // Structure : nP index 1 inde 2 ... index n functionIndex1 functionPointIndex1 parameterIndex1 ....  functionIndexn functionPointIndexn parameterIndexn
  // nP = 3*no.of.indeces
  // Find mapping between event variables and function to
  //printf("indices[0] = %.3f indices[1] = %.3f indices[2] = %.3f indices[2] = %.3f indices[3] = %.3f indices[4] = %.3f indices[5] = %.3f indices[6] = %.3f indices[7] = %.3f indices[8] = %.3f\n",indices[0],indices[1],indices[2],indices[3],indices[4],indices[5],indices[6],indices[7],indices[8]);
  //printf("p[0] = %.3f p[1] = %.3f p[2] = %.3f p[2] = %.3f p[3] = %.3f p[4] = %.3f p[5] = %.3f p[6] = %.3f p[7] = %.3f p[8] = %.3f\n",p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8]);

  unsigned int numParams = indices[0];
  unsigned int numIndex = numParams/3;
  fptype index = evt[indices[2 + indices[0]]];
  //std::vector <fptype> x;
  unsigned int functionIndex = 1;
  unsigned int compareIndex = 0;
  //printf(" %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e \n",p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],p[17],p[18],p[19]);
  //printf(" %u = %u = %u = %u = %u = %u = %u = %u = %u = %u = %u = %u = %u = %u = %u = %u = %u = %u = %u = %u \n %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e = %.2e \n %.2f - %.2f - %.2f ",indices[0],indices[1],indices[2],indices[3],indices[4],indices[5],indices[6],indices[7],indices[8],indices[9],indices[10],indices[11],indices[12],indices[13],indices[14],indices[15],p[indices[1]],p[indices[2]],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],p[15],p[16],p[17],evt[indices[2 + indices[0]]],evt[indices[2 + indices[0]]+1],evt[indices[2 + indices[0]]+2]);
  for (int i = 1; i <= numIndex; i++) {
    /* code */
    //printf("x = %.2f Parameter %d =  %.2f \n",x,i,p[indices[i]]);debug++;
    if(p[indices[i]] == index){
      compareIndex = i;
      //printf("Inside -> x = %.2f Parameter %d =  %.2f \n",x,i,p[indices[i]]);debug++;
    }
  }

  if (compareIndex == 0) return 1.0;
  //printf("here %d \n",debug);debug++;
  //unsigned int mapFunction = (int) FLOOR(0.5 + x
  //This is an index into the IndexPdf's list of functions
  //int targetFunction = (int) FLOOR(0.5 + (*(reinterpret_cast<device_function_ptr>(device_function_table[mapFunction])))(evt, p, paramIndices + indices[2]));
  functionIndex = numIndex+3*(compareIndex-1)+2;
  //unsigned int numObservables = indices[indices[0]+1];
  //printf("p1 = %.2f p2 = %.2f indices[0] = %d - numObservables = %u - numIndex = %u index = %.2f z = %.2f x = %.2f y = %.2f compareIndex = %u point[0] = %.2f functionIndex = %u\n",p[indices[1]],p[indices[2]],indices[0],numObservables,numIndex,index,evt[indices[2 + indices[0]]+1],evt[indices[2 + indices[0]]+2],evt[indices[2 + indices[0]]+3],compareIndex,(evt+compareIndex+1)[0],functionIndex);
  //printf("here %d \n",debug);debug++;
  //Because there are two pieces of information about each function
  //Because first function information begins at index 3
  //printf("here %d \n",debug);debug++;
  //fptype ret = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[targetFunction]])))(evt, p, paramIndices + indices[targetFunction + 1]);
  fptype ret = callFunction((evt+compareIndex+1), indices[functionIndex], indices[functionIndex + 1]);
  //fptype ret = 1.0;
  //printf("x = %.2f compareIndex = %u indices[0] = %.2f targetFunction = %d pdfX = %.2f \n",x,compareIndex,indices[0],targetFunction,evt[indices[2 + indices[targetFunction]]+compareIndex+targetFunction]);
  ret *= normalisationFactors[indices[functionIndex]];
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

  assert((t.size() == b.size()));

  for (vector<Variable*>::iterator v = b.begin(); v != b.end(); ++v) {
    //components.push_back(*f);
    pindices.push_back(registerParameter(*v));
  }

  std::set<int> functionIndicesUsed;
  for (vector<GooPdf*>::iterator f = t.begin(); f != t.end(); ++f) {
    components.push_back(*f);
    pindices.push_back((*f)->getFunctionIndex());
    pindices.push_back((*f)->getFunctionPntIndex());
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
