#include "QuadDimHistoKStarPdf.hh"
#include "MatrixPdf.hh"
#define JUMPH 22222222

MEM_CONSTANT fptype* dev_base_quadkstarhisto[5]; // Multiple histograms for the case of multiple PDFs

// EXEC_TARGET fptype cosTheta_FromMasses(const fptype sameSideM2, const fptype oppositeSideM2, const fptype psi_nSM2, const fptype motherM2, const fptype refM2, const fptype otherM2) {
//
//   fptype num = (sameSideM2/2)*(motherM2 + refM2 - oppositeSideM2) - (1./4.)*(motherM2 - psi_nSM2 + sameSideM2)*(sameSideM2 - otherM2 + refM2) ;
//   fptype denom2 = ((1./4.)*POW(motherM2 - psi_nSM2 + sameSideM2,2) - sameSideM2*motherM2) * ((1./4.)*POW(sameSideM2 - otherM2 + refM2,2) - sameSideM2*refM2) ;
//
//   return (num / SQRT(denom2)) ;
//
// }

EXEC_TARGET fptype device_QuadDimHistoKStarPdf (fptype* evt, fptype* p, unsigned int* indices) {
  // Structure is
  // nP totalHis./tograms (idx1 limit1 step1 bins1) (idx2 limit2 step2 bins2) nO o1 o2
  // where limit and step are indices into functorConstants.

  int numVars = 4;//(indices[0] - 1) / 4;
  int globalBin = 0;
  int previous = 1;
  int myHistogramIndex = indices[1];
  fptype binDistances[10]; // Ten dimensions should be more than enough!
  // Distance from bin center in units of bin width in each dimension.
  //int holdObs;

  fptype psi_nS = 1.0;

  fptype MPsi_nS = MJpsi;
  fptype mkp = evt[indices[2 + indices[0]]];
  fptype cJ = evt[indices[2 + indices[0]]+2];
  fptype mPsiP = evt[indices[2 + indices[0]]+1];
  fptype phi = evt[indices[2 + indices[0]]+3];

  fptype mKP2 = mkp*mkp;
  fptype mPsiP2 = mPsiP*mPsiP;
  fptype MPsi_nS2 = MPsi_nS*MPsi_nS;

  fptype cKs = cosTheta_FromMasses(mKP2, mPsiP2, MPsi_nS2, MBd2, MKaon2, MPion2);

  //printf("Quad : %.2f %.2f  %.2f  %.2f  %.2f  %.2f \n - mpk = %.2f (%.2f - %.2f) mPsi = %.2f cJ = %.2f cKs = %.2f phi = %.2f \n",evt[indices[2 + indices[0]]],evt[indices[2 + indices[0]]+1],evt[indices[2 + indices[0]]+2],evt[indices[2 + indices[0]]+3],evt[indices[2 + indices[0]]+4],evt[indices[2 + indices[0]]+5],mkp,MBd - MPsi_nS,MKaon + MPion,mPsiP,cJ,cKs,phi);

  if ((mkp < MKaon + MPion) || (mkp > MBd - MPsi_nS) || (mPsiP < MPsi_nS + MPion) || (mPsiP > MBd - MKaon)) {
//    printf("Returning 0: point out of the Dalitz borders!\n");
    return 0.; }
  if (FABS(cKs) > 1.0 || FABS(phi) > devPi || FABS(cJ) > 1.0)
    {
//	printf("Returning 0: point out of the angles borders!\n");
	return 0.;
	}

  // if (!(dalitz_contour_dev(mkp,mPsiP,false,psi_nS)))
  // {
  //   printf("mpk = %.2f (%.2f - %.2f) cJ = %.2f cKs = %.2f phi = %.2f \n",mkp,MBd - MPsi_nS,MKaon + MPion,cJ,cKs,phi);
  //   return 0.;
  // }

  //fptype one,two;
  unsigned int observablesSeen = 0;
  for (int i = 0; i < numVars; ++i) {

    fptype currVariable = 0;
    //unsigned int varIndex = indices[2 + 4*i];
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

  fptype* myHistogram = dev_base_quadkstarhisto[myHistogramIndex];
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

MEM_DEVICE device_function_ptr ptr_to_QuadDimHistoKStar = device_QuadDimHistoKStarPdf;

__host__ QuadDimHistoKStarPdf::QuadDimHistoKStarPdf (std::string n,
				       BinnedDataSet* x,
				       std::vector<Variable*> obses)
  : GooPdf(0, n)
  , numVars(x->numVariables())
{
   if(numVars!=4) abortWithCudaPrintFlush(__FILE__, __LINE__, "Only FOUR variables please !\n");

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
      pindices.push_back(JUMPH);
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
  MEMCPY_TO_SYMBOL(dev_base_quadkstarhisto, dev_address, sizeof(fptype*), totalHistograms*sizeof(fptype*), cudaMemcpyHostToDevice);
  // printf("NumVars Declaration = %d \n",numVars);
  GET_FUNCTION_ADDR(ptr_to_QuadDimHistoKStar);
  initialise(pindices);

  totalHistograms++;

  // printf("NumVars Declaration = %d \n",numVars);
}
