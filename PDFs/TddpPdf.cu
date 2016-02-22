#include "TddpPdf.hh"
#include <complex>
using std::complex; 

const int resonanceOffset = 8; // Offset of the first resonance into the parameter index array 
// Offset is number of parameters, constant index, indices for tau, xmix, and ymix, index
// of resolution function, and finally number of resonances (not calculable from nP
// because we don't know what the efficiency and time resolution might need). Efficiency 
// and time-resolution parameters are after the resonance information. 
const unsigned int SPECIAL_RESOLUTION_FLAG = 999999999; 

// The function of this array is to hold all the cached waves; specific 
// waves are recalculated when the corresponding resonance mass or width 
// changes. Note that in a multithread environment each thread needs its
// own cache, hence the '10'. Ten threads should be enough for anyone! 
MEM_DEVICE WaveHolder* cWaves[10]; 

EXEC_TARGET bool inDalitz (fptype m12, fptype m13, fptype bigM, fptype dm1, fptype dm2, fptype dm3) {
  if (m12 < POW(dm1 + dm2, 2)) return false; // This m12 cannot exist, it's less than the square of the (1,2) particle mass.
  if (m12 > POW(bigM - dm3, 2)) return false;   // This doesn't work either, there's no room for an at-rest 3 daughter. 
  
  // Calculate energies of 1 and 3 particles in m12 rest frame. 
  fptype e1star = 0.5 * (m12 - dm2*dm2 + dm1*dm1) / SQRT(m12); 
  fptype e3star = 0.5 * (bigM*bigM - m12 - dm3*dm3) / SQRT(m12); 

  // Bounds for m13 at this value of m12.
  fptype minimum = POW(e1star + e3star, 2) - POW(SQRT(e1star*e1star - dm1*dm1) + SQRT(e3star*e3star - dm3*dm3), 2);
  if (m13 < minimum) return false;
  fptype maximum = POW(e1star + e3star, 2) - POW(SQRT(e1star*e1star - dm1*dm1) - SQRT(e3star*e3star - dm3*dm3), 2);
  if (m13 > maximum) return false;

  return true; 
}

EXEC_TARGET inline int parIndexFromResIndex (int resIndex) {
  return resonanceOffset + resIndex*resonanceSize; 
}

EXEC_TARGET devcomplex<fptype> getResonanceAmplitude (fptype m12, fptype m13, fptype m23, 
						     unsigned int functionIdx, unsigned int pIndex) {
  resonance_function_ptr func = reinterpret_cast<resonance_function_ptr>(device_function_table[functionIdx]);
  return (*func)(m12, m13, m23, paramIndices + pIndex); 
}

EXEC_TARGET ThreeComplex device_Tddp_calcIntegrals (fptype m12, fptype m13, int res_i, int res_j, fptype* p, unsigned int* indices) {
  // For calculating Dalitz-plot integrals. What's needed is the products 
  // AiAj*, AiBj*, and BiBj*, where 
  // Ai = BW_i(x, y) + BW_i(y, x)
  // and Bi reverses the sign of the second BW. 
  // This function returns the above values at a single point. 
  // NB: Multiplication by efficiency is done by the calling function. 
  // Note that this function expects
  // to be called on a normalisation grid, not on 
  // observed points, that's why it doesn't use 
  // cWaves. No need to cache the values at individual
  // grid points - we only care about totals. 

  fptype motherMass = functorConstants[indices[1] + 0]; 
  fptype daug1Mass  = functorConstants[indices[1] + 1]; 
  fptype daug2Mass  = functorConstants[indices[1] + 2]; 
  fptype daug3Mass  = functorConstants[indices[1] + 3];  

  ThreeComplex ret; 
  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
  fptype m23 = motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - m12 - m13; 

  int parameter_i = parIndexFromResIndex(res_i);
  int parameter_j = parIndexFromResIndex(res_j);

  //fptype amp_real             = p[indices[parameter_i+0]];
  //fptype amp_imag             = p[indices[parameter_i+1]];
  unsigned int functn_i = indices[parameter_i+2];
  unsigned int params_i = indices[parameter_i+3];
  devcomplex<fptype> ai = getResonanceAmplitude(m12, m13, m23, functn_i, params_i);
  devcomplex<fptype> bi = getResonanceAmplitude(m13, m12, m23, functn_i, params_i);

  unsigned int functn_j = indices[parameter_j+2];
  unsigned int params_j = indices[parameter_j+3];
  devcomplex<fptype> aj = conj(getResonanceAmplitude(m12, m13, m23, functn_j, params_j));
  devcomplex<fptype> bj = conj(getResonanceAmplitude(m13, m12, m23, functn_j, params_j)); 

  ret = ThreeComplex((ai*aj).real, (ai*aj).imag, (ai*bj).real, (ai*bj).imag, (bi*bj).real, (bi*bj).imag);
  return ret; 
}

EXEC_TARGET fptype device_Tddp (fptype* evt, fptype* p, unsigned int* indices) {
  fptype motherMass = functorConstants[indices[1] + 0]; 
  fptype daug1Mass  = functorConstants[indices[1] + 1]; 
  fptype daug2Mass  = functorConstants[indices[1] + 2]; 
  fptype daug3Mass  = functorConstants[indices[1] + 3]; 

  fptype m12 = evt[indices[4 + indices[0]]]; 
  fptype m13 = evt[indices[5 + indices[0]]];

  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return 0; 
  int evtNum = (int) FLOOR(0.5 + evt[indices[6 + indices[0]]]); 

  devcomplex<fptype> sumWavesA(0, 0);
  devcomplex<fptype> sumWavesB(0, 0); 
  devcomplex<fptype> sumRateAA(0, 0); 
  devcomplex<fptype> sumRateAB(0, 0); 
  devcomplex<fptype> sumRateBB(0, 0); 

  unsigned int numResonances = indices[6]; 
  unsigned int cacheToUse    = indices[7]; 

  for (int i = 0; i < numResonances; ++i) {
    int paramIndex  = parIndexFromResIndex(i);
    fptype amp_real = p[indices[paramIndex+0]];
    fptype amp_imag = p[indices[paramIndex+1]];

    devcomplex<fptype> matrixelement(thrust::get<0>(cWaves[cacheToUse][evtNum*numResonances + i]),
				     thrust::get<1>(cWaves[cacheToUse][evtNum*numResonances + i])); 
    matrixelement.multiply(amp_real, amp_imag); 
    sumWavesA += matrixelement; 

#ifdef DEBUGSUMRATES
    if (25 > evtNum) {
      devcomplex<fptype> waveA_i(thrust::get<0>(cWaves[cacheToUse][evtNum*numResonances + i]),
				 thrust::get<1>(cWaves[cacheToUse][evtNum*numResonances + i])); 
      devcomplex<fptype> waveB_i(thrust::get<2>(cWaves[cacheToUse][evtNum*numResonances + i]),
				 thrust::get<3>(cWaves[cacheToUse][evtNum*numResonances + i])); 

      for (int j = 0; j < numResonances; ++j) {
	int paramIndex_j  = parIndexFromResIndex(j);
	fptype amp_real_j = p[indices[paramIndex_j+0]];
	fptype amp_imag_j = p[indices[paramIndex_j+1]];

	devcomplex<fptype> waveA_j(thrust::get<0>(cWaves[cacheToUse][evtNum*numResonances + j]),
				   thrust::get<1>(cWaves[cacheToUse][evtNum*numResonances + j])); 

	devcomplex<fptype> waveB_j(thrust::get<2>(cWaves[cacheToUse][evtNum*numResonances + j]),
				   thrust::get<3>(cWaves[cacheToUse][evtNum*numResonances + j])); 
	devcomplex<fptype> amps(amp_real, -amp_imag);
	amps.multiply(amp_real_j, amp_imag_j); 

	devcomplex<fptype> rateAA = conj(waveA_i)*waveA_j*amps;
	devcomplex<fptype> rateAB = conj(waveA_i)*waveB_j*amps;
	devcomplex<fptype> rateBB = conj(waveB_i)*waveB_j*amps;
	
	sumRateAA += rateAA;
	sumRateAB += rateAB;
	sumRateBB += rateBB;
      }
      waveA_i.multiply(amp_real, amp_imag);
      waveB_i.multiply(amp_real, amp_imag);
    }
#endif

    matrixelement = devcomplex<fptype>(thrust::get<2>(cWaves[cacheToUse][evtNum*numResonances + i]),
				       thrust::get<3>(cWaves[cacheToUse][evtNum*numResonances + i])); 
    matrixelement.multiply(amp_real, amp_imag); 
    sumWavesB += matrixelement; 
  } 

  fptype _tau     = p[indices[2]];
  fptype _xmixing = p[indices[3]];
  fptype _ymixing = p[indices[4]];
  
  fptype _time    = evt[indices[2 + indices[0]]];
  fptype _sigma   = evt[indices[3 + indices[0]]];

  //if ((gpuDebug & 1) && (0 == BLOCKIDX) && (0 == THREADIDX)) 
  //if (0 == evtNum) printf("TDDP: (%f, %f) (%f, %f)\n", sumWavesA.real, sumWavesA.imag, sumWavesB.real, sumWavesB.imag);
  //printf("TDDP: %f %f %f %f | %f %f %i\n", m12, m13, _time, _sigma, _xmixing, _tau, evtNum); 


  /*
  fptype ret = 0; 
  ret += (norm2(sumWavesA) + norm2(sumWavesB))*COSH(_ymixing * _time);
  ret += (norm2(sumWavesA) - norm2(sumWavesB))*COS (_xmixing * _time);
  sumWavesA *= conj(sumWavesB); 
  ret -= 2*sumWavesA.real * SINH(_ymixing * _time);
  ret -= 2*sumWavesA.imag * SIN (_xmixing * _time); // Notice sign difference wrt to Mikhail's code, because I have AB* and he has A*B. 
  ret *= EXP(-_time); 
  */

  fptype term1 = norm2(sumWavesA) + norm2(sumWavesB);
  fptype term2 = norm2(sumWavesA) - norm2(sumWavesB);
  sumWavesA *= conj(sumWavesB); 
  //printf("(%i, %i) TDDP: %f %f %f %f %f %f %f\n", BLOCKIDX, THREADIDX, term1, term2, sumWavesA.real, sumWavesA.imag, m12, m13, _tau);

  // Cannot use callFunction on resolution function. 
  int effFunctionIdx = parIndexFromResIndex(numResonances); 
  int resFunctionIdx = indices[5]; 
  int resFunctionPar = 2 + effFunctionIdx; 
  fptype ret = 0; 
  int md0_offset = 0; 
  if (resFunctionIdx == SPECIAL_RESOLUTION_FLAG) {
    // In this case there are multiple resolution functions, they are stored after the efficiency function,
    // and which one we use depends on the measured mother-particle mass. 
    md0_offset = 1; 
    fptype massd0 = evt[indices[7 + indices[0]]]; 
    fptype minMass = functorConstants[indices[1] + 6];
    fptype md0Step = functorConstants[indices[1] + 7];
    int res_to_use = (massd0 <= minMass) ? 0 : (int) FLOOR((massd0 - minMass) / md0Step); 
    int maxFcn     = indices[2+effFunctionIdx]; 
    if (res_to_use > maxFcn) res_to_use = maxFcn; 

    // Now calculate index of resolution function.
    // At the end of the array are indices efficiency_function, efficiency_parameters, maxFcn, res_function_1, res_function_1_nP, par1, par2 ... res_function_2, res_function_2_nP, ... 
    res_to_use = 3 + effFunctionIdx + res_to_use * (2 + indices[effFunctionIdx + 4]); 
    // NB this assumes all resolution functions have the same number of parameters. The number
    // of parameters in the first resolution function is stored in effFunctionIdx+3; add one to 
    // account for the index of the resolution function itself in the device function table, one
    // to account for the number-of-parameters index, and this is the space taken up by each 
    // resolution function. Multiply by res_to_use to get the number of spaces to skip to get to 
    // the one we want. 

    resFunctionIdx = indices[res_to_use]; 
    resFunctionPar = res_to_use + 1; 
  }
  
  ret = (*(reinterpret_cast<device_resfunction_ptr>(device_function_table[resFunctionIdx])))(term1, term2, sumWavesA.real, sumWavesA.imag,
											     _tau, _time, _xmixing, _ymixing, _sigma, 
											     p, indices + resFunctionPar); 
  
  // For the reversed (mistagged) fraction, we make the 
  // interchange A <-> B. So term1 stays the same, 
  // term2 changes sign, and AB* becomes BA*. 
  // Efficiency remains the same for the mistagged part,
  // because it depends on the momenta of the pi+ and pi-,
  // which don't change even though we tagged a D0 as D0bar. 
  
  fptype mistag = functorConstants[indices[1] + 5]; 
  if (mistag > 0) { // This should be either true or false for all events, so no branch is caused.
    // See header file for explanation of 'mistag' variable - it is actually the probability
    // of having the correct sign, given that we have a correctly reconstructed D meson. 
    mistag = evt[indices[md0_offset + 7 + indices[0]]]; 
    ret *= mistag; 
    ret += (1 - mistag) * (*(reinterpret_cast<device_resfunction_ptr>(device_function_table[resFunctionIdx])))(term1, -term2, sumWavesA.real, -sumWavesA.imag,
													   _tau, _time, _xmixing, _ymixing, _sigma, 
													   p, &(indices[resFunctionPar])); 
  }
   
  fptype eff = callFunction(evt, indices[effFunctionIdx], indices[effFunctionIdx + 1]); 
  //internalDebug = 0; 
  ret *= eff;

  //if ((gpuDebug & 1) && (4000 > evtNum)) {
    //printf("FULLPRINT: %i: %f %f\n", evtNum, ret*normalisationFactors[(indices - paramIndices)], eff);
    //printf("FULLPRINT: %i: %f %f (%f %f %f %f)\n", evtNum, ret, eff, m12, m13, _time, _sigma); 
  //} 



  //if (evtNum < 50) {
  //if ((gpuDebug & 1) && (0 == THREADIDX)) {
  //if ((gpuDebug & 1) && (180 == evtNum)) {
  //if ((0 == THREADIDX) && (0 == BLOCKIDX) && (gpuDebug & 1)) {
  //sumRateAA *= eff;
  //sumRateAB *= eff;
  //sumRateBB *= eff;
  //printf("Signal1 %i (%f, %f) (%f, %f) (%f, %f)\n", evtNum, sumRateAA.real, sumRateAA.imag, sumRateAB.real, sumRateAB.imag, sumRateBB.real, sumRateBB.imag); 
  //printf("Signal1 %i (%f, %f) (%f, %f) (%f, %f %f) %f\n", evtNum, term1, term2, sumWavesA.real, sumWavesA.imag, _tau, _xmixing, _ymixing, ret);
  //fptype m23 = motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - m12 - m13; 
  //printf("Signal2 %i (%f, %f, %f) %f, %f | %f %f %f\n", evtNum, m12, m13, m23, _time, _sigma, eff, ret, normalisationFactors[(indices - paramIndices)]);
  //printf("Signal3 %f %f %f %f %f %f %f %f\n", cudaArray[indices[effFunctionIdx+1]+1], cudaArray[indices[effFunctionIdx+1]+2], cudaArray[indices[effFunctionIdx+1]+3], 
  //cudaArray[indices[effFunctionIdx+1]+4], cudaArray[indices[effFunctionIdx+1]+5], cudaArray[indices[effFunctionIdx+1]+6], 
  //cudaArray[indices[effFunctionIdx+1]+7], cudaArray[indices[effFunctionIdx+1]+8]); 
  //}

  //printf("(%i, %i) TDDP: %f %f %f %f %f %i %f\n", BLOCKIDX, THREADIDX, _time, _sigma, m12, m13, term1, evtNum, ret);
  //if ((gpuDebug & 1) && (isnan(ret)))
  //printf("(%i, %i) TDDP: %f %f %f %f %i %i %f\n", BLOCKIDX, THREADIDX, _time, _sigma, m12, m13, evtNum, indices[6 + indices[0]], evt[indices[6 + indices[0]]]);
  //if ((gpuDebug & 1) && (isnan(ret)))
  //printf("(%i, %i) TDDP: %f %f %f %f %f %f %f\n", BLOCKIDX, THREADIDX, term1, term2, sumWavesA.real, sumWavesA.imag, _xmixing, _ymixing, _tau);

  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_Tddp = device_Tddp; 

__host__ TddpPdf::TddpPdf (std::string n, Variable* _dtime, Variable* _sigmat, Variable* m12, Variable* m13, Variable* eventNumber, DecayInfo* decay, MixingTimeResolution* r, GooPdf* efficiency, Variable* mistag) 
  : GooPdf(_dtime, n) 
  , decayInfo(decay)
  , _m12(m12)
  , _m13(m13)
  , dalitzNormRange(0)
  , cachedWaves(0) 
  , integrals(0)
  , resolution(r) 
  , forceRedoIntegrals(true)
  , totalEventSize(5) // Default 5 = m12, m13, time, sigma_t, evtNum 
  , cacheToUse(0) 
  , integrators(0)
  , calculators(0) 
{
  // NB, _dtime already registered!
  registerObservable(_sigmat);
  registerObservable(_m12);
  registerObservable(_m13);
  registerObservable(eventNumber); 

  fptype decayConstants[6];
  decayConstants[5] = 0; 
  
  if (mistag) {
    registerObservable(mistag);
    totalEventSize = 6; 
    decayConstants[5] = 1; // Flags existence of mistag
  } 

  std::vector<unsigned int> pindices;
  pindices.push_back(registerConstants(6)); 
  decayConstants[0] = decayInfo->motherMass;
  decayConstants[1] = decayInfo->daug1Mass;
  decayConstants[2] = decayInfo->daug2Mass;
  decayConstants[3] = decayInfo->daug3Mass;
  decayConstants[4] = decayInfo->meson_radius;
  MEMCPY_TO_SYMBOL(functorConstants, decayConstants, 6*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);  
  
  pindices.push_back(registerParameter(decayInfo->_tau));
  pindices.push_back(registerParameter(decayInfo->_xmixing));
  pindices.push_back(registerParameter(decayInfo->_ymixing));
  assert(resolution->getDeviceFunction() >= 0); 
  pindices.push_back((unsigned int) resolution->getDeviceFunction()); 
  pindices.push_back(decayInfo->resonances.size()); 

  static int cacheCount = 0; 
  cacheToUse = cacheCount++; 
  pindices.push_back(cacheToUse); 

  for (std::vector<ResonancePdf*>::iterator res = decayInfo->resonances.begin(); res != decayInfo->resonances.end(); ++res) {
    pindices.push_back(registerParameter((*res)->amp_real));
    pindices.push_back(registerParameter((*res)->amp_imag));
    pindices.push_back((*res)->getFunctionIndex());
    pindices.push_back((*res)->getParameterIndex());
    (*res)->setConstantIndex(cIndex); 
    components.push_back(*res);
  }

  pindices.push_back(efficiency->getFunctionIndex()); 
  pindices.push_back(efficiency->getParameterIndex());
  components.push_back(efficiency); 

  resolution->createParameters(pindices, this); 
  GET_FUNCTION_ADDR(ptr_to_Tddp);
  initialise(pindices);

  redoIntegral = new bool[decayInfo->resonances.size()];
  cachedMasses = new fptype[decayInfo->resonances.size()];
  cachedWidths = new fptype[decayInfo->resonances.size()];
  integrals    = new ThreeComplex**[decayInfo->resonances.size()];
  integrators  = new SpecialDalitzIntegrator**[decayInfo->resonances.size()];
  calculators  = new SpecialWaveCalculator*[decayInfo->resonances.size()];

  for (int i = 0; i < decayInfo->resonances.size(); ++i) {
    redoIntegral[i] = true;
    cachedMasses[i] = -1;
    cachedWidths[i] = -1; 
    integrators[i]  = new SpecialDalitzIntegrator*[decayInfo->resonances.size()];
    calculators[i]  = new SpecialWaveCalculator(parameters, i); 
    integrals[i]    = new ThreeComplex*[decayInfo->resonances.size()];
    
    for (int j = 0; j < decayInfo->resonances.size(); ++j) {
      integrals[i][j]   = new ThreeComplex(0, 0, 0, 0, 0, 0);
      integrators[i][j] = new SpecialDalitzIntegrator(parameters, i, j); 
    }
  }

  addSpecialMask(PdfBase::ForceSeparateNorm); 
}

__host__ TddpPdf::TddpPdf (std::string n, Variable* _dtime, Variable* _sigmat, Variable* m12, Variable* m13, Variable* eventNumber, DecayInfo* decay, vector<MixingTimeResolution*>& r, GooPdf* efficiency, Variable* md0, Variable* mistag) 
  : GooPdf(_dtime, n) 
  , decayInfo(decay)
  , _m12(m12)
  , _m13(m13)
  , dalitzNormRange(0)
  , cachedWaves(0) 
  , integrals(0)
  , resolution(r[0])  // Only used for normalisation, which only depends on x and y - it doesn't matter which one we use. 
  , forceRedoIntegrals(true)
  , totalEventSize(6) // This case adds the D0 mass by default. 
  , cacheToUse(0) 
  , integrators(0)
  , calculators(0) 
{
  // NB, _dtime already registered!
  registerObservable(_sigmat);
  registerObservable(_m12);
  registerObservable(_m13);
  registerObservable(eventNumber); 
  registerObservable(md0); 

  fptype decayConstants[8];
  decayConstants[5] = 0; 
  decayConstants[6] = md0->lowerlimit;
  decayConstants[7] = (md0->upperlimit - md0->lowerlimit) / r.size();
  
  if (mistag) {
    registerObservable(mistag);
    totalEventSize++; 
    decayConstants[5] = 1; // Flags existence of mistag
  } 

  std::vector<unsigned int> pindices;
  pindices.push_back(registerConstants(8)); 
  decayConstants[0] = decayInfo->motherMass;
  decayConstants[1] = decayInfo->daug1Mass;
  decayConstants[2] = decayInfo->daug2Mass;
  decayConstants[3] = decayInfo->daug3Mass;
  decayConstants[4] = decayInfo->meson_radius;
  MEMCPY_TO_SYMBOL(functorConstants, decayConstants, 8*sizeof(fptype), cIndex*sizeof(fptype), cudaMemcpyHostToDevice);  
  
  pindices.push_back(registerParameter(decayInfo->_tau));
  pindices.push_back(registerParameter(decayInfo->_xmixing));
  pindices.push_back(registerParameter(decayInfo->_ymixing));
  pindices.push_back(SPECIAL_RESOLUTION_FLAG); // Flag existence of multiple resolution functions.
  pindices.push_back(decayInfo->resonances.size()); 

  static int cacheCount = 0; 
  cacheToUse = cacheCount++; 
  pindices.push_back(cacheToUse); 

  for (std::vector<ResonancePdf*>::iterator res = decayInfo->resonances.begin(); res != decayInfo->resonances.end(); ++res) {
    pindices.push_back(registerParameter((*res)->amp_real));
    pindices.push_back(registerParameter((*res)->amp_imag));
    pindices.push_back((*res)->getFunctionIndex());
    pindices.push_back((*res)->getParameterIndex());
    (*res)->setConstantIndex(cIndex); 
    components.push_back(*res);
  }

  pindices.push_back(efficiency->getFunctionIndex());
  pindices.push_back(efficiency->getParameterIndex());
  components.push_back(efficiency); 

  pindices.push_back(r.size() - 1); // Highest index, not number of functions. 
  for (int i = 0; i < r.size(); ++i) {
    assert(r[i]->getDeviceFunction() >= 0); 
    pindices.push_back((unsigned int) r[i]->getDeviceFunction()); 
    r[i]->createParameters(pindices, this); 
  }

  GET_FUNCTION_ADDR(ptr_to_Tddp);
  initialise(pindices);

  redoIntegral = new bool[decayInfo->resonances.size()];
  cachedMasses = new fptype[decayInfo->resonances.size()];
  cachedWidths = new fptype[decayInfo->resonances.size()];
  integrals    = new ThreeComplex**[decayInfo->resonances.size()];
  integrators  = new SpecialDalitzIntegrator**[decayInfo->resonances.size()];
  calculators  = new SpecialWaveCalculator*[decayInfo->resonances.size()];

  for (int i = 0; i < decayInfo->resonances.size(); ++i) {
    redoIntegral[i] = true;
    cachedMasses[i] = -1;
    cachedWidths[i] = -1; 
    integrators[i]  = new SpecialDalitzIntegrator*[decayInfo->resonances.size()];
    calculators[i]  = new SpecialWaveCalculator(parameters, i); 
    integrals[i]    = new ThreeComplex*[decayInfo->resonances.size()];
    
    for (int j = 0; j < decayInfo->resonances.size(); ++j) {
      integrals[i][j]   = new ThreeComplex(0, 0, 0, 0, 0, 0);
      integrators[i][j] = new SpecialDalitzIntegrator(parameters, i, j); 
    }
  }

  addSpecialMask(PdfBase::ForceSeparateNorm); 
}

__host__ void TddpPdf::setDataSize (unsigned int dataSize, unsigned int evtSize) {
  // Default 5 is m12, m13, time, sigma_t, evtNum
  totalEventSize = evtSize;
  assert(totalEventSize >= 5); 

  if (cachedWaves) {
    delete cachedWaves;
  }

  numEntries = dataSize; 
  cachedWaves = new thrust::device_vector<WaveHolder>(dataSize*decayInfo->resonances.size());
  void* dummy = thrust::raw_pointer_cast(cachedWaves->data()); 
  MEMCPY_TO_SYMBOL(cWaves, &dummy, sizeof(WaveHolder*), cacheToUse*sizeof(WaveHolder*), cudaMemcpyHostToDevice); 
  setForceIntegrals(); 
}

__host__ fptype TddpPdf::normalise () const {
  recursiveSetNormalisation(1); // Not going to normalise efficiency, 
  // so set normalisation factor to 1 so it doesn't get multiplied by zero. 
  // Copy at this time to ensure that the SpecialWaveCalculators, which need the efficiency, 
  // don't get zeroes through multiplying by the normFactor. 
  MEMCPY_TO_SYMBOL(normalisationFactors, host_normalisation, totalParams*sizeof(fptype), 0, cudaMemcpyHostToDevice); 
  //std::cout << "TDDP normalisation " << getName() << std::endl;

  int totalBins = _m12->numbins * _m13->numbins;
  if (!dalitzNormRange) {
    gooMalloc((void**) &dalitzNormRange, 6*sizeof(fptype));
  
    fptype* host_norms = new fptype[6];
    host_norms[0] = _m12->lowerlimit;
    host_norms[1] = _m12->upperlimit;
    host_norms[2] = _m12->numbins;
    host_norms[3] = _m13->lowerlimit;
    host_norms[4] = _m13->upperlimit;
    host_norms[5] = _m13->numbins;
    MEMCPY(dalitzNormRange, host_norms, 6*sizeof(fptype), cudaMemcpyHostToDevice);
    delete[] host_norms; 
  }

  for (unsigned int i = 0; i < decayInfo->resonances.size(); ++i) {
    redoIntegral[i] = forceRedoIntegrals; 
    if (!(decayInfo->resonances[i]->parametersChanged())) continue;
    redoIntegral[i] = true; 
    decayInfo->resonances[i]->storeParameters();
  }
  forceRedoIntegrals = false; 

  // Only do this bit if masses or widths have changed.  
  thrust::constant_iterator<fptype*> arrayAddress(dalitzNormRange); 
  thrust::counting_iterator<int> binIndex(0); 

  // NB, SpecialWaveCalculator assumes that fit is unbinned! 
  // And it needs to know the total event size, not just observables
  // for this particular PDF component. 
  thrust::constant_iterator<fptype*> dataArray(dev_event_array); 
  thrust::constant_iterator<int> eventSize(totalEventSize);
  thrust::counting_iterator<int> eventIndex(0); 

  static int normCall = 0; 
  normCall++; 

  for (int i = 0; i < decayInfo->resonances.size(); ++i) {
    if (redoIntegral[i]) {
      
      thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, dataArray, eventSize)),
			thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
			strided_range<thrust::device_vector<WaveHolder>::iterator>(cachedWaves->begin() + i, cachedWaves->end(), decayInfo->resonances.size()).begin(), 
			*(calculators[i]));
      //std::cout << "Integral for resonance " << i << " " << numEntries << " " << totalEventSize << std::endl; 
    }
    
    // Possibly this can be done more efficiently by exploiting symmetry? 
    for (int j = 0; j < decayInfo->resonances.size(); ++j) {
      if ((!redoIntegral[i]) && (!redoIntegral[j])) continue; 
      ThreeComplex dummy(0, 0, 0, 0, 0, 0);
      SpecialComplexSum complexSum; 
      (*(integrals[i][j])) = thrust::transform_reduce(thrust::make_zip_iterator(thrust::make_tuple(binIndex, arrayAddress)),
						      thrust::make_zip_iterator(thrust::make_tuple(binIndex + totalBins, arrayAddress)),
						      *(integrators[i][j]), 
						      dummy, 
						      complexSum); 
      /*
      std::cout << "With resonance " << j << ": " 
		<< thrust::get<0>(*(integrals[i][j])) << " " 
		<< thrust::get<1>(*(integrals[i][j])) << " " 
		<< thrust::get<2>(*(integrals[i][j])) << " " 
		<< thrust::get<3>(*(integrals[i][j])) << " " 
		<< thrust::get<4>(*(integrals[i][j])) << " " 
		<< thrust::get<5>(*(integrals[i][j])) << std::endl; 
      */
    }
  }      

  // End of time-consuming integrals. 

  complex<fptype> integralA_2(0, 0);
  complex<fptype> integralB_2(0, 0);
  complex<fptype> integralABs(0, 0);
  for (unsigned int i = 0; i < decayInfo->resonances.size(); ++i) {
    int param_i = parameters + resonanceOffset + resonanceSize*i; 
    complex<fptype> amplitude_i(host_params[host_indices[param_i]], host_params[host_indices[param_i + 1]]);

    for (unsigned int j = 0; j < decayInfo->resonances.size(); ++j) {
      int param_j = parameters + resonanceOffset + resonanceSize*j; 
      complex<fptype> amplitude_j(host_params[host_indices[param_j]], -host_params[host_indices[param_j + 1]]); // Notice complex conjugation

      integralA_2 += (amplitude_i * amplitude_j * complex<fptype>(thrust::get<0>(*(integrals[i][j])), thrust::get<1>(*(integrals[i][j])))); 
      integralABs += (amplitude_i * amplitude_j * complex<fptype>(thrust::get<2>(*(integrals[i][j])), thrust::get<3>(*(integrals[i][j]))));
      integralB_2 += (amplitude_i * amplitude_j * complex<fptype>(thrust::get<4>(*(integrals[i][j])), thrust::get<5>(*(integrals[i][j]))));

      /*
      if (cpuDebug & 1) {
	int idx = i * decayInfo->resonances.size() + j;     
	if (0 == host_callnumber) std::cout << "Integral contribution " << i << ", " << j << " " << idx << " : "
					    << amplitude_i << " "
					    << amplitude_j << " (" 
					    << real(amplitude_i * amplitude_j * complex<fptype>(thrust::get<0>(*(integrals[i][j])), thrust::get<1>(*(integrals[i][j])))) << ", "
					    << imag(amplitude_i * amplitude_j * complex<fptype>(thrust::get<0>(*(integrals[i][j])), thrust::get<1>(*(integrals[i][j])))) << ") ("
					    << real(amplitude_i * amplitude_j * complex<fptype>(thrust::get<2>(*(integrals[i][j])), thrust::get<3>(*(integrals[i][j])))) << ", "
					    << imag(amplitude_i * amplitude_j * complex<fptype>(thrust::get<2>(*(integrals[i][j])), thrust::get<3>(*(integrals[i][j])))) << ") ("
					    << real(amplitude_i * amplitude_j * complex<fptype>(thrust::get<4>(*(integrals[i][j])), thrust::get<5>(*(integrals[i][j])))) << ", "
					    << imag(amplitude_i * amplitude_j * complex<fptype>(thrust::get<4>(*(integrals[i][j])), thrust::get<5>(*(integrals[i][j])))) << ") "
					    << thrust::get<0>(*(integrals[i][j])) << ", "
					    << thrust::get<1>(*(integrals[i][j])) << ") ("
					    << thrust::get<2>(*(integrals[i][j])) << ", "
					    << thrust::get<3>(*(integrals[i][j])) << ") ("
					    << thrust::get<4>(*(integrals[i][j])) << ", "
					    << thrust::get<5>(*(integrals[i][j])) << ") ("
					    << real(integralA_2) << ", " << imag(integralA_2) << ") "
					    << std::endl; 
      }
      */

    }
  }

  double dalitzIntegralOne = real(integralA_2); // Notice that this is already the abs2, so it's real by construction; but the compiler doesn't know that. 
  double dalitzIntegralTwo = real(integralB_2);
  double dalitzIntegralThr = real(integralABs);
  double dalitzIntegralFou = imag(integralABs); 
  
  fptype tau     = host_params[host_indices[parameters + 2]];
  fptype xmixing = host_params[host_indices[parameters + 3]];
  fptype ymixing = host_params[host_indices[parameters + 4]];

  fptype ret = resolution->normalisation(dalitzIntegralOne, dalitzIntegralTwo, dalitzIntegralThr, dalitzIntegralFou, tau, xmixing, ymixing); 

  double binSizeFactor = 1;
  binSizeFactor *= ((_m12->upperlimit - _m12->lowerlimit) / _m12->numbins);
  binSizeFactor *= ((_m13->upperlimit - _m13->lowerlimit) / _m13->numbins);
  ret *= binSizeFactor;

#ifdef CUDAPRINT 
  std::cout << "Tddp dalitz integrals: " 
	    << dalitzIntegralOne * binSizeFactor << " "
	    << dalitzIntegralTwo * binSizeFactor << " "
	    << dalitzIntegralThr * binSizeFactor << " "
	    << dalitzIntegralFou * binSizeFactor << " | "
	    << binSizeFactor << " " 
	    << tau << " " << xmixing << " " << ymixing << " "
	    << ret << " " 
	    << std::endl; 
#endif
    
  host_normalisation[parameters] = 1.0/ret;
  //std::cout << "End of TDDP normalisation: " << ret << " " << host_normalisation[parameters] << " " << binSizeFactor << std::endl; 
  return (fptype) ret; 
}
//#endif 

SpecialDalitzIntegrator::SpecialDalitzIntegrator (int pIdx, unsigned int ri, unsigned int rj) 
  : resonance_i(ri)
  , resonance_j(rj)
  , parameters(pIdx) 
{}

EXEC_TARGET ThreeComplex SpecialDalitzIntegrator::operator () (thrust::tuple<int, fptype*> t) const {
  // Bin index, base address [lower, upper, numbins] 
  // Notice that this is basically MetricTaker::operator (binned) with the special-case knowledge
  // that event size is two, and that the function to call is dev_Tddp_calcIntegrals.

  int globalBinNumber  = thrust::get<0>(t);
  fptype lowerBoundM12 = thrust::get<1>(t)[0];
  fptype upperBoundM12 = thrust::get<1>(t)[1];  
  int numBinsM12       = (int) FLOOR(thrust::get<1>(t)[2] + 0.5); 
  int binNumberM12     = globalBinNumber % numBinsM12;
  fptype binCenterM12  = upperBoundM12 - lowerBoundM12;
  binCenterM12        /= numBinsM12;
  binCenterM12        *= (binNumberM12 + 0.5); 
  binCenterM12        += lowerBoundM12; 

  globalBinNumber     /= numBinsM12; 
  fptype lowerBoundM13 = thrust::get<1>(t)[3];
  fptype upperBoundM13 = thrust::get<1>(t)[4];  
  int numBinsM13       = (int) FLOOR(thrust::get<1>(t)[5] + 0.5); 
  fptype binCenterM13  = upperBoundM13 - lowerBoundM13;
  binCenterM13        /= numBinsM13;
  binCenterM13        *= (globalBinNumber + 0.5); 
  binCenterM13        += lowerBoundM13; 

  //if (0 == THREADIDX) cuPrintf("%i %i %i %f %f operator\n", thrust::get<0>(t), thrust::get<0>(t) % numBinsM12, globalBinNumber, binCenterM12, binCenterM13);
  unsigned int* indices = paramIndices + parameters;   
  ThreeComplex ret = device_Tddp_calcIntegrals(binCenterM12, binCenterM13, resonance_i, resonance_j, cudaArray, indices); 

  fptype fakeEvt[10]; // Need room for many observables in case m12 or m13 were assigned a high index in an event-weighted fit. 
  fakeEvt[indices[indices[0] + 2 + 2]] = binCenterM12;
  fakeEvt[indices[indices[0] + 2 + 3]] = binCenterM13;
  unsigned int numResonances = indices[6]; 
  int effFunctionIdx = parIndexFromResIndex(numResonances); 
  //if (thrust::get<0>(t) == 19840) {internalDebug1 = BLOCKIDX; internalDebug2 = THREADIDX;}
  //fptype eff = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[effFunctionIdx]])))(fakeEvt, cudaArray, paramIndices + indices[effFunctionIdx + 1]);
  fptype eff = callFunction(fakeEvt, indices[effFunctionIdx], indices[effFunctionIdx + 1]); 
  //if (thrust::get<0>(t) == 19840) {
  //internalDebug1 = -1; 
  //internalDebug2 = -1;
    //printf("Efficiency: %i %f %f %f %i\n", thrust::get<0>(t), binCenterM12, binCenterM13, eff, effFunctionIdx);
    //printf("Efficiency: %f %f %f %f %f %i %i\n", fakeEvt[0], fakeEvt[1], fakeEvt[2], fakeEvt[3], fakeEvt[4], indices[indices[0] + 2 + 2], indices[indices[0] + 2 + 3]); 
  //}

  // Multiplication by eff, not sqrt(eff), is correct:
  // These complex numbers will not be squared when they
  // go into the integrals. They've been squared already,
  // as it were. 
  thrust::get<0>(ret) *= eff;
  thrust::get<1>(ret) *= eff;
  thrust::get<2>(ret) *= eff;
  thrust::get<3>(ret) *= eff;
  thrust::get<4>(ret) *= eff;
  thrust::get<5>(ret) *= eff;
  return ret; 
}

SpecialWaveCalculator::SpecialWaveCalculator (int pIdx, unsigned int res_idx) 
  : resonance_i(res_idx)
  , parameters(pIdx)
{}

EXEC_TARGET WaveHolder SpecialWaveCalculator::operator () (thrust::tuple<int, fptype*, int> t) const {
  // Calculates the BW values for a specific resonance. 
  // The 'A' wave stores the value at each point, the 'B' 
  // at the opposite (reversed) point. 

  WaveHolder ret;

  int evtNum = thrust::get<0>(t); 
  fptype* evt = thrust::get<1>(t) + (evtNum * thrust::get<2>(t)); 

  unsigned int* indices = paramIndices + parameters;   // Jump to TDDP position within parameters array
  fptype m12 = evt[indices[4 + indices[0]]]; 
  fptype m13 = evt[indices[5 + indices[0]]];

  fptype motherMass = functorConstants[indices[1] + 0]; 
  fptype daug1Mass  = functorConstants[indices[1] + 1]; 
  fptype daug2Mass  = functorConstants[indices[1] + 2]; 
  fptype daug3Mass  = functorConstants[indices[1] + 3];  

  if (!inDalitz(m12, m13, motherMass, daug1Mass, daug2Mass, daug3Mass)) return ret;
  fptype m23 = motherMass*motherMass + daug1Mass*daug1Mass + daug2Mass*daug2Mass + daug3Mass*daug3Mass - m12 - m13; 

  int parameter_i = parIndexFromResIndex(resonance_i); // Find position of this resonance relative to TDDP start 
  unsigned int functn_i = indices[parameter_i+2];
  unsigned int params_i = indices[parameter_i+3];

  devcomplex<fptype> ai = getResonanceAmplitude(m12, m13, m23, functn_i, params_i);
  devcomplex<fptype> bi = getResonanceAmplitude(m13, m12, m23, functn_i, params_i);

  //printf("Amplitudes %f, %f => (%f %f) (%f %f)\n", m12, m13, ai.real, ai.imag, bi.real, bi.imag);

  thrust::get<0>(ret) = ai.real;
  thrust::get<1>(ret) = ai.imag;
  thrust::get<2>(ret) = bi.real;
  thrust::get<3>(ret) = bi.imag; 

  return ret; 
}

