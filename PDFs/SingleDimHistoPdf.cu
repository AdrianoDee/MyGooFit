#include "SingleDimHistoPdf.hh"

MEM_CONSTANT fptype* dev_base_sidimhisto[100]; // Multiple histograms for the case of multiple PDFs



  //  EXEC_TARGET fptype interpolateArrays (fptype* xArray, fptype* yArray, int intOrder,fptype xvalue)
  //  {
   //
  //     printf("Bin histo pdf 3.0 = %.3f %d \n",xvalue,intOrder);
   //
   //
  //     fptype den,dif,dift,ho,hp,w,y,dy;
  //     fptype coeffC[20], coeffD[20];
   //
  //     dif = fabs(xvalue-xArray[0]) ;
   //
  //     int ns=1;
   //
  //     for(int intexInter =1 ; intexInter<=intOrder+1 ; ++intexInter)
  //     {
  //       dift=fabs(xvalue-xArray[intexInter-1]);
  //       if (dift<dif)
  //       {
  //          ns = intexInter;
  //          dif = dift ;
  //       }
   //
  //       coeffC[intexInter] = yArray[intexInter-1];
  //       coeffD[intexInter] = yArray[intexInter-1];
   //
  //       printf("Bin histo pdf 3 = %.3f %d %.3f %.3f %.3f %.3f \n",xvalue,intexInter,dift,dif,xArray[intexInter-1],coeffC[intexInter],coeffD[intexInter]);
   //
   //
  //     }
   //
  //     y=yArray[--ns] ;
   //
  //     for(int m=1 ; m<intOrder+1; m++)
  //     {
  //       for(int intexInter=1 ; intexInter<=intOrder+1-m ; intexInter++)
  //       {
  //         ho=xArray[intexInter-1]-xvalue ;
  //         hp=xArray[intexInter-1+m]-xvalue ;
  //         w=coeffC[intexInter+1]-coeffD[intexInter] ;
  //         den=ho-hp ;
  //         if (den==0.)
  //         {
  //           return 0. ;
  //         }
  //         den = w/den ;
  //         coeffD[intexInter]=hp*den ;
  //         coeffC[intexInter]=ho*den;
  //         }
  //         dy = (2*ns)<(intOrder+1-m) ? coeffC[ns+1] : coeffD[ns--] ;
  //         y += dy ;
   //
  //         printf("Bin histo pdf 4 = %.3f %.3f %.3f %.3f %.3f %.3f\n",xvalue,ho,hp,w,den,dy);
   //
   //
  //       }
   //
  //       return y;
   //
  //  }


  //  EXEC_TARGET fptype interSingleDimension (int localNumBins, fptype step, fptype lowerBound, fptype xval, int intOrder,fptype* histogram)
  //  {
   //
  //    int localBin    = (int) FLOOR((xval-lowerBound)/step); // Int_t fbinC = dim.getBin(*binning) ;
   //
  //    fptype binCenter = (fptype)localBin*step+lowerBound-0.5*step;
  //    fptype upperBound   = lowerBound + step*localNumBins;
   //
  //    int binOffset = (xval<binCenter)? 1 : 0;
  //    int fbinLo  = localBin - intOrder/2 - binOffset;//Int_t fbinLo = fbinC-intOrder/2 - ((xval<binning->binCenter(fbinC))?1:0) ;
   //
  //    fptype xarr[20];
  //    fptype yarr[20];
   //
  //    printf("Bin histo pdf 1 = %.3f %d %.3f %d %.3f %.3f %.3f %d \n",xval,localBin,binCenter,fbinLo,lowerBound,step,upperBound,intOrder);
   //
  //    for (int index=fbinLo ; index<=intOrder+fbinLo ; ++index)
  //    {
  //      int ibin ;
  //      if (index>=0 && index<localNumBins) {
  //        ibin = index;
  //        xarr[index-fbinLo] = lowerBound+ibin*step-step*0.5;
  //        yarr[index-fbinLo] = histogram[ibin];
  //        printf("Bin histo pdf 2 = %.3f %d %d %d %d %.3f %.3f \n",xval,localBin,index,ibin,localNumBins,xarr[index-fbinLo],histogram[ibin]);
  //      } else if (index>=localNumBins) {
  //       //  ibin = 2*localNumBins-index-1 ;
  //        printf("Over binning 2 \n");
  //        xarr[index-fbinLo] = upperBound+(1e-10)*(index-localNumBins+1);
  //        yarr[index-fbinLo] = 0.0 ;
  //      } else {
  //        printf("Under binning 2 \n");
  //        ibin = -index - 1 ;
  //        xarr[index-fbinLo] = lowerBound-ibin*(1e-10);
  //        yarr[index-fbinLo] = 0.0 ;
  //      }
  //    }
   //
  //    fptype ret = interpolateArrays(xarr,yarr,intOrder+1,xval);
   //
  //    return ret;
   //
  //  }

  //  EXEC_TARGET fptype interSingleDimensionMulti (int otherBin, int localNumBins,int otherNumBins, fptype step, fptype lowerBound, fptype xval, int intOrder,fptype* histogram)
  //  {
   //
  //    int localBin    = (int) FLOOR((xval-lowerBound)/step); // Int_t fbinC = dim.getBin(*binning) ;
   //
  //    fptype binCenter = (fptype)localBin*step+lowerBound-0.5*step;
  //    fptype upperBound   = lowerBound + step*localNumBins;
   //
  //    int binOffset = (xval<binCenter)? 1 : 0;
  //    int fbinLo  = localBin - intOrder/2 - binOffset;//Int_t fbinLo = fbinC-intOrder/2 - ((xval<binning->binCenter(fbinC))?1:0) ;
   //
  //    fptype xarr[20];
  //    fptype yarr[20];
   //
  //    printf("Bin histo Multi pdf 1 = %.3f %d %.3f %d %.3f %.3f %.3f %d \n",xval,localBin,binCenter,fbinLo,lowerBound,step,upperBound,intOrder);
   //
  //    for (int index=fbinLo ; index<=intOrder+fbinLo ; ++index)
  //    {
  //      int ibin ;
  //      int globalBin;
   //
  //      if (index>=0 && index<localNumBins) {
  //        ibin = index;
  //        xarr[index-fbinLo] = lowerBound+ibin*step-step*0.5;
  //        globalBin = otherBin + otherNumBins*ibin;
  //        yarr[index-fbinLo] = histogram[globalBin];
  //        printf("Bin histo Multi pdf 2 = %.3f %d %d %d %d %.3f %.3f \n",xval,localBin,index,ibin,localNumBins,xarr[index-fbinLo],histogram[ibin]);
  //      } else if (index>=localNumBins) {
  //       //  ibin = 2*localNumBins-index-1 ;
  //        xarr[index-fbinLo] = upperBound+(1e-10)*(index-localNumBins+1);
  //        yarr[index-fbinLo] = 0.0 ;
  //      } else {
  //        ibin = -index - 1 ;
  //        xarr[index-fbinLo] = lowerBound-ibin*(1e-10);
  //        yarr[index-fbinLo] = 0.0 ;
  //      }
  //    }
   //
  //    fptype ret = interpolateArrays(xarr,yarr,intOrder+1,xval);
   //
  //    return ret;
   //
  //  }

   EXEC_TARGET fptype device_SiDimHistoPdf (fptype* evt, fptype* p, unsigned int* indices) {
     // Structure is
     // nP totalHistograms interPolationOrder (limit1 step1 bins1) (limit2 step2 bins2) nO o1 o2
     // where limit and step are indices into functorConstants.

     int numVars = (indices[0] - 1) / 3;

     printf("NumVars = %d \n",numVars);

     int globalBin = 0;
     int previousNofBins = 1;
     int myHistogramIndex = indices[1];
     int interpolationOrder = indices[2];

     fptype* myHistogram = dev_base_sidimhisto[myHistogramIndex];

     if(numVars==1)
     {
       printf("NumVars = 1\n");

       int i = 0;

       int localNumBins = indices[3*(i+1) + 1 + 1];
       int lowerBoundIdx   = 2 + 3*i + 1;
       fptype lowerBound   = functorConstants[indices[lowerBoundIdx + 0]];
       fptype step         = functorConstants[indices[lowerBoundIdx + 1]];
       fptype upperBound   = lowerBound + step*localNumBins;

       fptype currVariable = evt[indices[indices[0] + 2]];

       if(currVariable<lowerBound || currVariable >upperBound) return 0.0;

       fptype ret = interSingleDimension(localNumBins, step, lowerBound, currVariable, interpolationOrder, myHistogram);

       return ret;

     }
     else
      return 0.0;

     }

MEM_DEVICE device_function_ptr ptr_to_SiDimHistogram = device_SiDimHistoPdf;

__host__ SingleDimHistoPdf::SingleDimHistoPdf (std::string n,
							 BinnedDataSet* x,
							 std::vector<Variable*> obses, unsigned int interOrder)
  : GooPdf(0, n)
  , numVars(x->numVariables())
{
  if(numVars>2) abortWithCudaPrintFlush(__FILE__, __LINE__, "Only the first two variables will be taken into account !\n");
  if(interOrder>20) abortWithCudaPrintFlush(__FILE__, __LINE__, "Interpolation order must be smaller than 20! \n");

  printf("NumVars Declaration = %d \n",numVars);

  int numConstants = 2*numVars+1;
  registerConstants(numConstants);
  static unsigned int totalHistograms = 0;
  host_constants = new fptype[numConstants];
  totalEvents = 0;

  std::vector<unsigned int> pindices;
  pindices.push_back(totalHistograms);

  host_constants[0] = interOrder;

  pindices.push_back(interOrder);

  int varIndex = 0;
  for (varConstIt var = x->varsBegin(); var != x->varsEnd(); ++var) {
    if (std::find(obses.begin(), obses.end(), *var) != obses.end()) {
      registerObservable(*var);
    }
    else {
      abortWithCudaPrintFlush(__FILE__, __LINE__, "The BinnedDataSet provided variables are different from p.d.f. observables \n");
    }

    pindices.push_back(cIndex + 2*varIndex + 0 + 1); //cIndex is no. of constants index
    pindices.push_back(cIndex + 2*varIndex + 1 + 1);
    pindices.push_back((*var)->numbins);

    // NB, do not put cIndex here, it is accounted for by the offset in MEMCPY_TO_SYMBOL below.
    host_constants[2*varIndex + 0 + 1] = (*var)->lowerlimit;
    host_constants[2*varIndex + 1 + 1] = ((*var)->upperlimit - (*var)->lowerlimit) / (*var)->numbins;
    varIndex++;
  }

  unsigned int numbins = x->getNumBins();
  thrust::host_vector<fptype> host_histogram;
  for (unsigned int i = 0; i < numbins; ++i) {
    fptype curr = x->getBinContent(i);
    host_histogram.push_back(curr);
    totalEvents += curr;
  }
  MEMCPY_TO_SYMBOL(functorConstants, host_constants, numConstants*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);

  dev_base_histogram = new thrust::device_vector<fptype>(host_histogram);
  static fptype* dev_address[1];
  dev_address[0] = (&((*dev_base_histogram)[0])).get();
  MEMCPY_TO_SYMBOL(dev_base_sidimhisto, dev_address, sizeof(fptype*), totalHistograms*sizeof(fptype*), cudaMemcpyHostToDevice);
  GET_FUNCTION_ADDR(ptr_to_SiDimHistogram);
  initialise(pindices);

  totalHistograms++;
}
