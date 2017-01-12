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


   EXEC_TARGET fptype lastInter(int intOrder,int nS,fptype* xArray,fptype* yArray,fptype* cC, fptype* cD,fptype xV)
   {

      fptype ho = 0.0;
      fptype hp = 0.0;
      fptype w = 0.0;

      fptype den = 0.0;

      fptype y = yArray[--nS];

      fptype dy = 0.0;

      //
      for(int m=1 ; m<intOrder+1; m++)
      {
        for(int intexInter=1 ; intexInter<=intOrder+1-m ; intexInter++)
        {
          ho=xArray[intexInter-1]-xV ;
          hp=xArray[intexInter-1+m]-xV ;
          w=cC[intexInter+1]-cD[intexInter] ;
          den=ho-hp ;
          if (den==0.)
          {
            return 0. ;
          }
          den = w/den ;
          cD[intexInter]=hp*den ;
          cC[intexInter]=ho*den;
        }

           dy = (2*nS)<(intOrder+1-m) ? cC[nS+1] : cD[nS--] ;
           y += dy ;

        // printf("Bin histo pdf 4 = %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",xval,ho,hp,w,den,dy,y);


       }


   }

   EXEC_TARGET fptype interSingleDimensionMulti (fptype yval, int otherBin, int localNumBins,int otherNumBins, fptype step, fptype lowerBound, fptype xval, int intOrder,fptype* histogram)
   {

     fptype den = 0.0 ,dif = 0.0 ,dift = 0.0 ,ho = 0.0 ,hp = 0.0 ,w = 0.0 ,y = 0.0 ,dy = 0.0 ;
     fptype coeffC[30], coeffD[30];

     int intexInter = 1;
     int ns=1;
     int m = 1;
     int index = 1;

     fptype xarr[30];
     fptype yarr[30];

     for (size_t i = 0; i < 30; i++) {
       coeffC[i] = 0.0;
       coeffD[i] = 0.0;

       xarr[i] = 0.0;
       yarr[i] = 0.0;

     }

     int localBin    = (int) FLOOR((xval-lowerBound)/step); // Int_t fbinC = dim.getBin(*binning) ;

     fptype binCenter = (fptype)localBin*step+lowerBound-0.5*step;
     fptype upperBound   = lowerBound + step*localNumBins;

     int binOffset = (xval<binCenter)? 1 : 0;
     int fbinLo  = localBin - intOrder/2 - binOffset;//Int_t fbinLo = fbinC-intOrder/2 - ((xval<binning->binCenter(fbinC))?1:0) ;

     int gB = fbinLo + otherNumBins*otherBin;

    //  // printf("Bin histo Multi pdf 1 = %.3f %.3f %d %.3f %d %.3f %.3f %.3f %d %d %d %d %.3f\n",xval,yval,localBin,binCenter,fbinLo,lowerBound,step,upperBound,intOrder,otherBin,otherNumBins,gB,histogram[gB]);
      // printf("Bin histo Multi pdf 1 = %.3f %.3f %d %.3f %d %.3f %.3f %.3f %d %d %d %d %.3f\n",xval,yval,localBin,binCenter,fbinLo,lowerBound,step,upperBound,intOrder,otherBin,otherNumBins,gB,xval);

     for (index=fbinLo ; index<=intOrder+fbinLo ; ++index)
     {
       int ibin = 0;
       int globalBin = 0;

       if (index>=0 && index<localNumBins) {
         ibin = index;
         xarr[index-fbinLo] = lowerBound+ibin*step+step*0.5;
         globalBin = ibin + localNumBins*otherBin;
         yarr[index-fbinLo] = histogram[globalBin];

         // printf("Bin histo Multi pdf 2 = %.3f %.3f %d %d %d %d %d %d %d %.2f %.3f %.3f \n",xval,yval,localBin, otherBin, otherNumBins,globalBin, index,ibin,localNumBins,(fptype) sizeof(histogram),xarr[index-fbinLo],histogram[0]);
       } else if (index>=localNumBins) {
        //  ibin = 2*localNumBins-index-1 ;
         xarr[index-fbinLo] = upperBound+(1e-10)*(index-localNumBins+1);
         yarr[index-fbinLo] = 0.0 ;
       } else {
         ibin = -index - 1 ;
         xarr[index-fbinLo] = lowerBound-ibin*(1e-10);
         yarr[index-fbinLo] = 0.0 ;
       }
     }

     __syncthreads();


     dif = fabs(xval-xarr[0]) ;


     for(int intexInter =1 ; intexInter<=intOrder+1 ; ++intexInter)
     {
       // printf("Bin histo pdf 3.1 = %.3f %d %.3f %d %d\n",xval,intOrder,dif,ns,intexInter);
       dift=fabs(xval-xarr[intexInter-1]);
       // printf("Bin histo pdf 3.2 = %.3f %d %.3f %d %.3f %d\n",xval,intOrder,dif,ns,dift,intexInter);

       if (dift<dif)
       {
          ns = intexInter;
          dif = dift ;
       }

      //  // printf("Bin histo pdf 3.3 = %.3f %d %.3f %d %d\n",xval,intOrder,dif,ns,intexInter);

        coeffC[intexInter] = yarr[intexInter-1];
        coeffD[intexInter] = yarr[intexInter-1];

        // // printf("Bin histo pdf 3.4 = %.3f %d %.3f %.3f %.3f %.3f \n",xval,intexInter,dift,dif,xarr[intexInter-1],coeffC[intexInter],coeffD[intexInter]);


     }

     fptype ret = lastInter(intOrder,ns,xarr,yarr,coeffC,coeffD,xval);

     return ret;

   }

   EXEC_TARGET fptype device_BiDimHistosPdf (fptype* evt, fptype* p, unsigned int* indices) {

     int numVars = 2;// (indices[0] - 1) / 3;


     int globalBin = 0;
     int previousNofBins = 1;
     int myHistogramIndex = indices[1];
     int interpolationOrder = indices[2];

     // printf("Indices %d %d %d %d %d %d %d %d %d \n",indices[0],indices[1],indices[2],indices[3],indices[4],indices[5],indices[6],indices[7],indices[8]);

     fptype* myHistogram = dev_base_bidimhisto[myHistogramIndex];

     if(numVars==2)
     {
       // printf("NumVars = 2\n");
       fptype var[2],lowerBound[2],step[2],upperBound[2],binCenter[2];
       int bins[2],localBin[2],binOffset[2];

       for (int i = 0; i < numVars; ++i) {

         var[i] = evt[indices[indices[0] + 2] + i];
         bins[i] = indices[3*(i+1) + 1 + 1];
         int lowerBoundIdx   = 2 + 3*i + 1;
         lowerBound[i]   = functorConstants[indices[lowerBoundIdx + 0]];
         step[i]         = functorConstants[indices[lowerBoundIdx + 1]];
         upperBound[i]   = lowerBound[i] + step[i]*bins[i];
         localBin[i]     = (int) FLOOR((var[i]-lowerBound[i])/step[i]);
         binCenter[i]    = (fptype)localBin[i]*step[i]+lowerBound[i]-0.5*step[i];
         binOffset[i]    = (var[i]<binCenter[i])? 1 : 0;

       }

       // printf("Bin histo Multi pdf var0 = %.3f %d %.3f %.3f %.3f %d %.3f %.3f %d \n",var[0],bins[0],lowerBound[0],step[0],upperBound[0],localBin[0],binCenter[0],binOffset[0],indices[2]);
       // printf("Bin histo Multi pdf var1 = %.3f %d %.3f %.3f %.3f %d %.3f %.3f %d \n",var[1],bins[1],lowerBound[1],step[1],upperBound[1],localBin[1],binCenter[1],binOffset[1],indices[2]);

       int ybinLo = localBin[1]-indices[2]/2 - binOffset[1];
       int ybinLoUp = ybinLo+1;
       int iBin = ybinLo;

       fptype yarr[2];
       fptype xarr[2];

       if (ybinLo>=0 && ybinLo<bins[1]) xarr[0] = lowerBound[1] + ybinLo*step[1]+step[1]*0.5;
       else if(ybinLo>bins[1])
       {
         iBin = 2*bins[1]-ybinLo-1;
         xarr[0] = upperBound[1] - lowerBound[1] + iBin*step[1]+step[1]*0.5;
       }
       else if(ybinLo<0)
       {
         iBin = -ybinLo -1;
         xarr[0] = lowerBound[1] + lowerBound[1] + iBin*step[1]+step[1]*0.5;
       }

       yarr[0] = interSingleDimensionMulti(xarr[0],iBin, bins[0], bins[1], step[0], lowerBound[0], var[0], indices[2], myHistogram);

       iBin = ybinLoUp;

       if (ybinLoUp>=0 && ybinLoUp<bins[1]) xarr[1] = lowerBound[1] + iBin*step[1]+step[1]*0.5;
       else if(ybinLoUp>bins[1])
       {
         iBin = 2*bins[1]-ybinLo-1;
         xarr[1] = upperBound[1] - lowerBound[1] + iBin*step[1]+step[1]*0.5;
       }
       else if(ybinLoUp<0)
       {
         iBin = -ybinLo -1;
         xarr[1] = lowerBound[1] + lowerBound[1] + iBin*step[1]+step[1]*0.5;
       }

       yarr[1] = interSingleDimensionMulti(xarr[1],iBin, bins[0], bins[1], step[0], lowerBound[0], var[0], indices[2], myHistogram);

       //  // printf("Interpolation : %d %d %.3f %.3f \n",iBin,yIndex,xarr[yIndex-ybinLo],yarr[yIndex-ybinLo]);



        //fptype ret = interpolateArrays(xarr,yarr,interpolationOrder+1,var[1]);

        fptype ret = 1.0;

        return ret;


     }
     else
      return 0.0;

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

       fptype holdcurrVariable;

       //fptype one,two;
       unsigned int observablesSeen = 0;
       for (int i = 0; i < numVars; ++i) {
         fptype currVariable = 0;
         unsigned int varIndex = indices[2 + 4*i];
         currVariable = evt[indices[indices[0] + 2 + observablesSeen++]];

         int lowerBoundIdx   = 3 + 4*i;
         fptype lowerBound   = functorConstants[indices[lowerBoundIdx + 0]];
         fptype step         = functorConstants[indices[lowerBoundIdx + 1]];

         holdcurrVariable = currVariable;

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
  if(numVars>2) abortWithCudaPrintFlush(__FILE__, __LINE__, "Only the first two variables will be taken into account !\n");

  int numConstants = 2*numVars;
  registerConstants(numConstants);
  static unsigned int totalHistograms = 0;
  host_constants = new fptype[numConstants];
  totalEvents = 0;

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
