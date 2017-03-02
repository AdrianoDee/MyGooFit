#include "BiDimHistoPdf.hh"
#define JUMP 22222222
MEM_CONSTANT fptype* dev_base_bidimhisto[100]; // Multiple histograms for the case of multiple PDFs



   EXEC_TARGET fptype interpolateArrays (fptype* xArray, fptype* yArray, int intOrder,fptype xvalue)
   {

      printf("Bin histo pdf 3.0 = %.3f %d \n",xvalue,intOrder);


      fptype den,dif,dift,ho,hp,w,y,dy;
      fptype coeffC[20], coeffD[20];

      dif = fabs(xvalue-xArray[0]) ;

      int ns=1;

      for(int intexInter =1 ; intexInter<=intOrder+1 ; ++intexInter)
      {
        dift=fabs(xvalue-xArray[intexInter-1]);
        if (dift<dif)
        {
           ns = intexInter;
           dif = dift ;
        }

        coeffC[intexInter] = yArray[intexInter-1];
        coeffD[intexInter] = yArray[intexInter-1];

        // printf("Bin histo pdf 3 = %.3f %d %.3f %.3f %.3f %.3f \n",xvalue,intexInter,dift,dif,xArray[intexInter-1],coeffC[intexInter],coeffD[intexInter]);


      }

      y=yArray[--ns] ;

      for(int m=1 ; m<intOrder+1; m++)
      {
        for(int intexInter=1 ; intexInter<=intOrder+1-m ; intexInter++)
        {
          ho=xArray[intexInter-1]-xvalue ;
          hp=xArray[intexInter-1+m]-xvalue ;
          w=coeffC[intexInter+1]-coeffD[intexInter] ;
          den=ho-hp ;
          if (den==0.)
          {
            return 0. ;
          }
          den = w/den ;
          coeffD[intexInter]=hp*den ;
          coeffC[intexInter]=ho*den;
          }
          dy = (2*ns)<(intOrder+1-m) ? coeffC[ns+1] : coeffD[ns--] ;
          y += dy ;

          // printf("Bin histo pdf 4 = %.3f %.3f %.3f %.3f %.3f %.3f\n",xvalue,ho,hp,w,den,dy);


        }

        return y;

   }

     EXEC_TARGET fptype device_BiDimHistoPdf (fptype* evt, fptype* p, unsigned int* indices) {
       // Structure is
       // nP totalHis./tograms (idx1 limit1 step1 bins1) (idx2 limit2 step2 bins2) nO o1 o2
       // where limit and step are indices into functorConstants.

       int numVars = (indices[0] - 1) / 4;
       int globalBin = 0;
       int previous = 1;
       int myHistogramIndex = indices[1];
       fptype binDistances[10]; // Ten dimensions should be more than enough!
       // Distance from bin center in units of bin width in each dimension.
       //int holdObs;

       //fptype one,two;
       unsigned int observablesSeen = 0;
       for (int i = 0; i < numVars; ++i) {
         fptype currVariable = 0;

         currVariable = evt[indices[indices[0] + 2 + observablesSeen++]];

         int lowerBoundIdx   = 3 + 4*i;
         fptype lowerBound   = functorConstants[indices[lowerBoundIdx + 0]];
         fptype step         = functorConstants[indices[lowerBoundIdx + 1]];

         currVariable   -= lowerBound;
         currVariable   /= step;

         int localNumBins = indices[4*(i+1) + 1];
         int localBin    = (int) FLOOR(currVariable);
         binDistances[i] = currVariable - localBin - fptype(0.5);
         globalBin      += previous * localBin;
         previous       *= localNumBins;
       }

       fptype* myHistogram = dev_base_bidimhisto[myHistogramIndex];
       fptype ret = 0;

       fptype totalWeight = 0;
       int totalBins = dev_powi(3, numVars);
       for (int i = 0; i < totalBins; ++i) {
         int currBin = globalBin;
         int localPrevious = 1;
         int trackingBin = globalBin;
         bool offSomeAxis = false;
         fptype currentWeight = 0;

         for (int v = 0; v < numVars; ++v) {
           int localNumBins = indices[4*(v+1) + 1];
           int offset = ((i / dev_powi(3, v)) % 3) - 1;

           currBin += offset * localPrevious;
           localPrevious *= localNumBins;

           int currVarBin = trackingBin % localNumBins;
           trackingBin /= localNumBins;
           if (currVarBin + offset < 0) offSomeAxis = true;
           if (currVarBin + offset >= localNumBins) offSomeAxis = true;

           fptype currDist = binDistances[v];
           currDist -= offset;
           currentWeight += currDist*currDist;

         }

         // Only interpolate the four closest boxes (in two dimensions; more in three dimensions).
         currentWeight = currentWeight > 0 ? (currentWeight <= SQRT((fptype) numVars) ? 1 / SQRT(currentWeight) : 0) : 0;
         fptype currentEntry = offSomeAxis ? 0 : myHistogram[currBin];
         ret += currentWeight * currentEntry;
         totalWeight += currentWeight;

       }

       ret /= totalWeight;
       return ret;
     }

MEM_DEVICE device_function_ptr ptr_to_BiDimHistogram = device_BiDimHistoPdf;

__host__ BiDimHistoPdf::BiDimHistoPdf (std::string n,
							 BinnedDataSet* x,
							 std::vector<Variable*> obses)
  : GooPdf(0, n)
  , numVars(x->numVariables())
{
  // if(numVars>2) abortWithCudaPrintFlush(__FILE__, __LINE__, "Only the first two variables will be taken into account !\n");

  int numConstants = 2*numVars;
  registerConstants(numConstants);
  static unsigned int totalHistograms = 0;
  host_constants = new fptype[numConstants];
  totalEvents = 0;

  std::cout<<"Num histograms : "<<totalHistograms<<std::endl;

  std::vector<unsigned int> pindices;
  pindices.push_back(totalHistograms);

  int varIndex = 0;
  for (varConstIt var = x->varsBegin(); var != x->varsEnd(); ++var) {
    if (std::find(obses.begin(), obses.end(), *var) != obses.end()) {
      registerObservable(*var);
      pindices.push_back(JUMP);
    }
    else {
      abortWithCudaPrintFlush(__FILE__, __LINE__, "The BinnedDataSet provided variables are different from p.d.f. observables \n");
    }

    pindices.push_back(cIndex + 2*varIndex + 0 ); //cIndex is no. of constants index
    pindices.push_back(cIndex + 2*varIndex + 1 );
    pindices.push_back((*var)->numbins);

    // NB, do not put cIndex here, it is accounted for by the offset in MEMCPY_TO_SYMBOL below.
    host_constants[2*varIndex + 0 ] = (*var)->lowerlimit;
    host_constants[2*varIndex + 1 ] = ((*var)->upperlimit - (*var)->lowerlimit) / (*var)->numbins;
    varIndex++;
  }

  unsigned int numbins = x->getNumBins();
  thrust::host_vector<fptype> host_histogram;
  for (unsigned int i = 0; i < numbins; ++i) {
    fptype curr = x->getBinContent(i);
    host_histogram.push_back(curr);
    totalEvents += curr;
  }


  // printf("NumVars Declaration = %d \n",numVars);

  MEMCPY_TO_SYMBOL(functorConstants, host_constants, numConstants*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);
  // printf("NumVars Declaration = %d \n",numVars);
  dev_base_histogram = new thrust::device_vector<fptype>(host_histogram);
  static fptype* dev_address[1];
  dev_address[0] = (&((*dev_base_histogram)[0])).get();
  MEMCPY_TO_SYMBOL(dev_base_bidimhisto, dev_address, sizeof(fptype*), totalHistograms*sizeof(fptype*), cudaMemcpyHostToDevice);
  // printf("NumVars Declaration = %d \n",numVars);
  GET_FUNCTION_ADDR(ptr_to_BiDimHistogram);
  initialise(pindices);

  totalHistograms++;

  // printf("NumVars Declaration = %d \n",numVars);
}
