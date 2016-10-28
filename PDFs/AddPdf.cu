#include "AddPdf.hh"

//K5
 MEM_CONSTANT fptype xK5[3] = {0.9258200997725514615666,0.5773502691896257645091,0.00};
 MEM_CONSTANT fptype wK5[3] = {0.197979797979797979798,0.4909090909090909090909,0.6222222222222222222222};
 
//K7 
 MEM_CONSTANT fptype xK7[7]= {0.9604912687080202834235,0.7745966692414833770359,0.4342437493468025580021,0.0,-0.9604912687080202834235,-0.7745966692414833770359,-0.4342437493468025580021};
 MEM_CONSTANT fptype wK7[7]= {0.1046562260264672651938,0.268488089868333440729,0.401397414775962222905,0.450916538658474142345,0.1046562260264672651938,0.268488089868333440729,0.401397414775962222905};

//K33
 MEM_CONSTANT fptype xK33[33] = {0.998239274145444514183282371262429,0.989400934991649932596154173450333,0.971505950969392594304109858259183,0.944575023073232576077988415534608,0.909157667012342948622557332882809,0.865631202387831743880467897712393,0.814240287062444468094572577588398,0.755404408355003033895101194847442,0.689741106681762303876199858016583,0.617876244402643748446671764048791,0.5404076763521397436254283440515912,0.458016777657227386342419442983578,0.3714837808784162886840279470795229,0.2816035507792589132304605014604961,0.1891685790180837263147120866349424,0.0950125098376374401853193354249581,0.00,-0.998239274145444514183282371262429,-0.989400934991649932596154173450333,-0.971505950969392594304109858259183,-0.944575023073232576077988415534608,-0.909157667012342948622557332882809,-0.865631202387831743880467897712393,-0.814240287062444468094572577588398,-0.755404408355003033895101194847442,-0.689741106681762303876199858016583,-0.617876244402643748446671764048791,-0.5404076763521397436254283440515912,-0.458016777657227386342419442983578,-0.3714837808784162886840279470795229,-0.2816035507792589132304605014604961,-0.1891685790180837263147120866349424,-0.0950125098376374401853193354249581};
 MEM_CONSTANT fptype wK33[33] = {0.004742777049247317906344087722423,0.01325793068809115724543298539375,0.0224988594400494440294746719427368,0.03126054364738052823966115413873805,0.039512951202421964003078606929551,0.0475062159764070172349863273944227,0.0552056330954221723031200013936293,0.062358806011834855470397014742265,0.0688629951915312430013066966746289,0.074769823885599557533866452931194,0.080053941263719286445983306344537,0.084595803792590637585459382901613,0.08833750257911273141357403046803,0.09129203282819166227222659416544712,0.0934386740609212304781471900796438,0.0947284012472300413267339687994708,0.0951542160804983070204150559558403,0.004742777049247317906344087722423,0.01325793068809115724543298539375,0.0224988594400494440294746719427368,0.03126054364738052823966115413873805,0.039512951202421964003078606929551,0.0475062159764070172349863273944227,0.0552056330954221723031200013936293,0.062358806011834855470397014742265,0.0688629951915312430013066966746289,0.074769823885599557533866452931194,0.080053941263719286445983306344537,0.084595803792590637585459382901613,0.08833750257911273141357403046803,0.09129203282819166227222659416544712,0.0934386740609212304781471900796438,0.0947284012472300413267339687994708};

 //Legendre integration
 //fptype xLegs[16] = {0.095012509837637,0.281603550779259,0.458016777657227,0.617876244402644,0.755404408355003,0.865631202387832,0.944575023073233 ,0.989400934991650,-0.095012509837637,-0.281603550779259,-0.458016777657227,-0.617876244402644,-0.755404408355003,-0.865631202387832,-0.944575023073233,-0.989400934991650};
 //fptype wLegs[16] = {0.189450610455069,0.182603415044924,0.169156519395003,0.149595988816577,0.124628971255534,0.095158511682493,0.062253523938648,0.027152459411754,0.189450610455069,0.182603415044924,0.169156519395003,0.149595988816577,0.124628971255534,0.095158511682493,0.062253523938648,0.027152459411754};
 
//Legendre integration
 MEM_CONSTANT fptype xK16[16] = {0.095012509837637,0.281603550779259,0.458016777657227,0.617876244402644,0.755404408355003,0.865631202387832,0.944575023073233 ,0.989400934991650,-0.095012509837637,-0.281603550779259,-0.458016777657227,-0.617876244402644,-0.755404408355003,-0.865631202387832,-0.944575023073233,-0.989400934991650};
 MEM_CONSTANT fptype wK16[16] = {0.189450610455069,0.182603415044924,0.169156519395003,0.149595988816577,0.124628971255534,0.095158511682493,0.062253523938648,0.027152459411754,0.189450610455069,0.182603415044924,0.169156519395003,0.149595988816577,0.124628971255534,0.095158511682493,0.062253523938648,0.027152459411754};
 


EXEC_TARGET fptype device_AddPdfs_Bin (fptype* evt, fptype* p, unsigned int* indices) {
  int numParameters = indices[0];
  fptype totalWeight = 0;
  fptype xmid = evt[indices[2 + indices[0]]];;
  //fptype xnex = evt[indices[2 + indices[0]]+3]; // 3 = 1 observable (x) + binvalue + binvolume
  fptype step = evt[indices[2 + indices[0]]+2];
  //ONLY for binned datasets
  fptype binLow = xmid-step/2.0;
  fptype binHig = xmid+step/2.0;
  fptype lH = (binHig-binLow)/2.0;

  //printf("low = %f -- high = %f \n",binLow,binHig);
  
  fptype fX[4][16]={0};
  
 fptype* xKro;
 fptype* wKro;
 int kronMax = 0;
	
  if(xmid<=1.05||xmid>=1.65){
	  xKro = xK16;
	  wKro = wK16;
	  kronMax =16;
	  }else{
		xKro = xK7;
		wKro = wK7;
		kronMax = 7;
	  }
	  
	 //printf("KronMaax = %d\n",kronMax);
 
 fptype Xs[16];
 
 for(int b =0;b<kronMax;b++){ Xs[b] = xmid + lH*xKro[b];};

  //STEPS FOR INTEGRATION
  //fptype Xs[5];
  //printf("----INTEGRAL X Steps----\n");
  /*for (int k=0;k<=n;k++){

    Xs[k]=(-2+k)*step+xmid;

  }*/

  //printf("----FUNCTIONS 1----\n");
 // f
  int counter = 0;
  for (int i = 1; i < numParameters-4; i += 4) {
    totalWeight += p[indices[i+3]];
    for(int h=0;h<kronMax;h++){

    fptype curr= callFunction(Xs+h, indices[i+1], indices[i+2]);
    fptype weight = p[indices[i+3]];
    fX[counter][h] = weight * curr * normalisationFactors[indices[i+2]];

    //if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
    //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    //printf("Add comp %i: %f * %f * %f = %f (%f)\n", i, weight, curr, normalisationFactors[indices[i+1]], weight*curr*normalisationFactors[indices[i+1]], ret);
    }
	++counter;
  }
  // numParameters does not count itself. So the array structure for two functions is
  // nP | F P w | F P
  // in which nP = 5. Therefore the parameter index for the last function pointer is nP, and the function index is nP-1.
  //fptype last = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[numParameters-1]])))(evt, p, paramIndices + indices[numParameters]);
  //printf("----FUNCTIONS 2----\n");
  for(int h=0;h<kronMax;h++){
  fptype last = callFunction(Xs+h, indices[numParameters - 1], indices[numParameters]);
  fX[counter][h] = (1 - totalWeight) * last * normalisationFactors[indices[numParameters]];
  }

  fptype fXtot[16];

  //printf("----FUNCTIONS tot----\n");
  for(int l=0;l<kronMax;l++){
    fXtot[l]=fX[0][l]+fX[1][l];
  }
  //INTEGRAL
  fptype integralF = 0;

   for (int k=0;k<kronMax;k++){
        integralF += wKro[k]*fXtot[k];
		//printf("|| Integral F = %f || \n",integralF);
   }

	fptype integralTot = integralF*lH;
   
    //printf("|| xmid = %f legeH = %f legeK = %f fX3 = %f fX4 =%f fX5 = %f fX6 = %f w3 = %f w4 = %f w5 = %f w6 = %f integral = %f||\n",xmid,legeH,legeK,fXtot[2],fXtot[3],fXtot[4],fXtot[5],wLegs[2],wLegs[3],wLegs[4],wLegs[5],integralTot);
    
	return  integralTot;
   //printf(".");
  //if ((THREADIDX < 50) && (isnan(ret))) printf("NaN final component %f %f\n", last, totalWeight);

  //if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
  //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))

}

EXEC_TARGET fptype device_AddPdfs (fptype* evt, fptype* p, unsigned int* indices) {

  int numParameters = indices[0];
  fptype ret = 0;
  fptype totalWeight = 0;
  for (int i = 1; i < numParameters-3; i += 3) {
    totalWeight += p[indices[i+2]];
    fptype curr = callFunction(evt, indices[i], indices[i+1]);
    fptype weight = p[indices[i+2]];
    ret += weight * curr * normalisationFactors[indices[i+1]];

    //if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
    //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    //printf("Add comp %i: %f * %f * %f = %f (%f)\n", i, weight, curr, normalisationFactors[indices[i+1]], weight*curr*normalisationFactors[indices[i+1]], ret);

  }

 // printf("device_AddPdfs");
  // numParameters does not count itself. So the array structure for two functions is
  // nP | F P w | F P
  // in which nP = 5. Therefore the parameter index for the last function pointer is nP, and the function index is nP-1.
  //fptype last = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[numParameters-1]])))(evt, p, paramIndices + indices[numParameters]);
  fptype last = callFunction(evt, indices[numParameters - 1], indices[numParameters]);
  ret += (1 - totalWeight) * last * normalisationFactors[indices[numParameters]];

  //if ((THREADIDX < 50) && (isnan(ret))) printf("NaN final component %f %f\n", last, totalWeight);

  //if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
  //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
  //printf("Add final: %f * %f * %f = %f (%f)\n", (1 - totalWeight), last, normalisationFactors[indices[numParameters]], (1 - totalWeight) *last* normalisationFactors[indices[numParameters]], ret);

  return ret;
}

EXEC_TARGET fptype device_AddPdfs_Bin_Mid (fptype* evt, fptype* p, unsigned int* indices) {

   //FOR INTEGRAL CONSTRUCTOR
   //nP | f fp P W | f fp P W | f fp P
   // 0 | 1  2 3 4 | 5  6 7 8 | 9 10 11

  int numParameters = indices[0];
  fptype ret = 0;
  fptype totalWeight = 0;
  for (int i = 1; i < numParameters-4; i += 4) {
    totalWeight += p[indices[i+3]];
    fptype curr = callFunction(evt, indices[i], indices[i+2]);
    fptype weight = p[indices[i+3]];
    ret += weight * curr * normalisationFactors[indices[i+2]];

    //if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
    //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    //printf("Add comp %i: %f * %f * %f = %f (%f)\n", i, weight, curr, normalisationFactors[indices[i+1]], weight*curr*normalisationFactors[indices[i+1]], ret);

  }

  // numParameters does not count itself. So the array structure for two functions is
  // nP | F P w | F P
  // in which nP = 5. Therefore the parameter index for the last function pointer is nP, and the function index is nP-1.
  //fptype last = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[numParameters-1]])))(evt, p, paramIndices + indices[numParameters]);
  fptype last = callFunction(evt, indices[numParameters - 1], indices[numParameters]);
  ret += (1 - totalWeight) * last * normalisationFactors[indices[numParameters]];

  //fptype pdf1 = ret;
  //fptype pdf2 =  (1 - totalWeight) * last * normalisationFactors[indices[numParameters]];
  //printf("pdf1 = %f pdf2 = %f \n",pdf1,pdf2);

  //if ((THREADIDX < 50) && (isnan(ret))) printf("NaN final component %f %f\n", last, totalWeight);

  //if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
  //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
  //printf("Add final: %f * %f * %f = %f (%f)\n", (1 - totalWeight), last, normalisationFactors[indices[numParameters]], (1 - totalWeight) *last* normalisationFactors[indices[numParameters]], ret);

  return ret;
}

EXEC_TARGET fptype device_AddPdfsExt (fptype* evt, fptype* p, unsigned int* indices) {
  // numParameters does not count itself. So the array structure for two functions is
  // nP | F P w | F P w
  // in which nP = 6.
  int numParameters = indices[0];
  fptype ret = 0;
  fptype totalWeight = 0;
  for (int i = 1; i < numParameters; i += 3) {
    //fptype curr = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[i]])))(evt, p, paramIndices + indices[i+1]);
    fptype curr = callFunction(evt, indices[i], indices[i+1]);
    fptype weight = p[indices[i+2]];
    ret += weight * curr * normalisationFactors[indices[i+1]];

    totalWeight += weight;
    //if ((gpuDebug & 1) && (THREADIDX == 0) && (0 == BLOCKIDX))
    //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    //printf("AddExt: %i %E %f %f %f %f %f %f\n", i, curr, weight, ret, totalWeight, normalisationFactors[indices[i+1]], evt[0], evt[8]);
  }
  ret /= totalWeight;
  //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
  //if ((gpuDebug & 1) && (THREADIDX == 0) && (0 == BLOCKIDX))
  //printf("AddExt result: %f\n", ret);

  return ret;
}

//ONLY FOR BINNED DATASET -> BIN INTEGRAL!
EXEC_TARGET fptype device_AddPdfsExt_Bin (fptype* evt, fptype* p, unsigned int* indices) {
  //printf("INTEGRATION EXT\n");
  fptype fX [3][16]={0};//MATRIX of FUNCTIONS -> ROW = COMPONENT COLUMN = BIN STEP
  // | f1(x0) f1(x1) f1(x2) f1(x3) f1(x4) ....|
  // | f2(x0) f2(x1) f2(x2) f2(x3) f2(x4) ....|
  // | f3(x0) f3(x1) f3(x2) f3(x3) f3(x4) ....|
  // etc.
  // | f(x0)  f(x1)  f(x2)  f(x3)  f(x4) | LAST ROW tot function f(x0) = w1f1(x0) +  w2f2(x0) + ..... + wnfn(x0)
  int numParameters = indices[0];
  fptype totalWeight = 0;

  //ONLY for binned datasets
  fptype xmid = evt[indices[2 + indices[0]]];;
  //fptype xnex = evt[indices[2 + indices[0]]+3]; // 3 = 1 observable (x) + binvalue + binvolume
  fptype step = evt[indices[2 + indices[0]]+2];
  
   fptype binLow = xmid-step/2.0;
	fptype binHig = xmid+step/2.0;
	fptype legeH = (binHig-binLow)/2.0;
	fptype legeK = (binHig+binLow)/2.0;
  
   //Legendre integration
	fptype xLegs[16] = {0.095012509837637,0.281603550779259,0.458016777657227,0.617876244402644,0.755404408355003,0.865631202387832,0.944575023073233 ,0.989400934991650,-0.095012509837637,-0.281603550779259,-0.458016777657227,-0.617876244402644,-0.755404408355003,-0.865631202387832,-0.944575023073233,-0.989400934991650};
	fptype wLegs[16] = {0.189450610455069,0.182603415044924,0.169156519395003,0.149595988816577,0.124628971255534,0.095158511682493,0.062253523938648,0.027152459411754,0.189450610455069,0.182603415044924,0.169156519395003,0.149595988816577,0.124628971255534,0.095158511682493,0.062253523938648,0.027152459411754};
	fptype Xs[16];
	
	for(int b =0;b<16;b++){ Xs[b] = legeK + legeH*xLegs[b];};

  int compCnt = 0;  //TO COUNT THE NUMBER OF COMPONENTS
  //np | F Fp P W | F Fp P W |
  //0  | 1  2 3 4 | 5  6 7 8 |

  for (int i = 1; i < numParameters-1; i +=4) {
    fptype weight = p[indices[i+3]];
    totalWeight += p[indices[i+3]];
    //CYCLE ON THE X STEPS
    for(int h=0;h<=16;h++){
    //printf("Xs+%d = %f",h,*(Xs+h));
    fptype curr= callFunction(Xs+h, indices[i+1], indices[i+2]);
    fX[compCnt][h] = weight * curr * normalisationFactors[indices[i+2]];
    //double test= weight * curr * normalisationFactors[indices[i+2]];
    //printf("curr = %f \n",curr);
    //printf("weight = %f \n",weight);
    //printf("norm = %f \n",normalisationFactors[indices[i+2]]);
    //printf("test = %f\n",test);
    //if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
    //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    //printf("Add comp %i: %f * %f * %f = %f (%f)\n", i, weight, curr, normalisationFactors[indices[i+1]], weight*curr*normalisationFactors[indices[i+1]], ret);
    }

    compCnt++;
  }
  // numParameters does not count itself. So the array structure for two functions is
  // nP | F P w | F P
  // in which nP = 5. Therefore the parameter index for the last function pointer is nP, and the function index is nP-1.
  //fptype last = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[numParameters-1]])))(evt, p, paramIndices + indices[numParameters]);
 // printf("----FUNCTIONS 2----\n");
 // for(int k=0;k<20;k++) fX[k]=0;

 // printf("----FUNCTIONS tot----\n");
 
 fptype fXtot[16]={0};
 
 //TOTAL FUNCTION EVALUATION
  for(int p=0;p<=16;p++){
  for(int l=0;l<compCnt;l++){

    fXtot[p] += fX[l][p];

  //  printf("f[%d] = %f \n",l,fX[l]);
  }}
  
  //INTEGRAL
  fptype integralF = 0;
  
  for (int k=0;k<16;k++){
        integralF += wLegs[k]*fXtot[k];
		//printf("|| Integral F = %f || \n",integralF);
   }
   
   fptype integralTot = integralF*legeH;
   
   //printf("|| xmid = %f legeH = %f legeK = %f fX3 = %f fX4 =%f fX5 = %f fX6 = %f w3 = %f w4 = %f w5 = %f w6 = %f integral = %f||\n",xmid,legeH,legeK,fXtot[2],fXtot[3],fXtot[4],fXtot[5],wLegs[2],wLegs[3],wLegs[4],wLegs[5],integralTot);
    
	return  integralTot;


  //if ((THREADIDX < 50) && (isnan(ret))) printf("NaN final component %f %f\n", last, totalWeight);

  //if ((gpuDebug & 1) && (0 == THREADIDX) && (0 == BLOCKIDX))
  //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
  //printf("Add final: %f * %f * %f = %f (%f)\n", (1 - totalWeight), last, normalisationFactors[indices[numParameters]], (1 - totalWeight) *last* normalisationFactors[indices[numParameters]], ret);

}

//ONLY FOR BINNED DATASET -> BIN INTEGRAL!
EXEC_TARGET fptype device_AddPdfsExt_Bin_Mid (fptype* evt, fptype* p, unsigned int* indices)  {
  // numParameters does not count itself. So the array structure for two functions is
  // nP | F Fp P w | F Fp P w
  // in which nP = 6.
  int numParameters = indices[0];
  fptype ret = 0;
  fptype totalWeight = 0;
  for (int i = 1; i < numParameters; i += 4) {
    //fptype curr = (*(reinterpret_cast<device_function_ptr>(device_function_table[indices[i]])))(evt, p, paramIndices + indices[i+1]);
    fptype curr = callFunction(evt, indices[i], indices[i+1]);
    fptype weight = p[indices[i+2]];
    ret += weight * curr * normalisationFactors[indices[i+1]];

    totalWeight += weight;
    //if ((gpuDebug & 1) && (THREADIDX == 0) && (0 == BLOCKIDX))
    //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
    //printf("AddExt: %i %E %f %f %f %f %f %f\n", i, curr, weight, ret, totalWeight, normalisationFactors[indices[i+1]], evt[0], evt[8]);
  }
  ret /= totalWeight;
  //if ((1 > (int) floor(0.5 + evt[8])) && (gpuDebug & 1) && (paramIndices + debugParamIndex == indices))
  //if ((gpuDebug & 1) && (THREADIDX == 0) && (0 == BLOCKIDX))
  //printf("AddExt result: %f\n", ret);

  return ret;
}

MEM_DEVICE device_function_ptr ptr_to_AddPdfs = device_AddPdfs;
MEM_DEVICE device_function_ptr ptr_to_AddPdfs_Bin_Mid = device_AddPdfs_Bin_Mid;
MEM_DEVICE device_function_ptr ptr_to_AddPdfs_Bin = device_AddPdfs_Bin;
MEM_DEVICE device_function_ptr ptr_to_AddPdfsExt_Bin_Mid = device_AddPdfsExt_Bin_Mid;
MEM_DEVICE device_function_ptr ptr_to_AddPdfsExt = device_AddPdfsExt;
MEM_DEVICE device_function_ptr ptr_to_AddPdfsExt_Bin = device_AddPdfsExt_Bin;

AddPdf::AddPdf (std::string n, std::vector<Variable*> weights, std::vector<PdfBase*> comps)
  : GooPdf(0, n)
  , extended(true)
{

  assert((weights.size() == comps.size()) || (weights.size() + 1 == comps.size()));

  // Indices stores (function index)(function parameter index)(weight index) triplet for each component.
  // Last component has no weight index unless function is extended.
  for (std::vector<PdfBase*>::iterator p = comps.begin(); p != comps.end(); ++p) {
    components.push_back(*p);
    assert(components.back());
  }

  getObservables(observables);

  std::vector<unsigned int> pindices;
  for (unsigned int w = 0; w < weights.size(); ++w) {
    assert(components[w]);
    pindices.push_back(components[w]->getFunctionIndex());
    pindices.push_back(components[w]->getParameterIndex());
    pindices.push_back(registerParameter(weights[w]));
  }
  assert(components.back());
  if (weights.size() < components.size()) {
    pindices.push_back(components.back()->getFunctionIndex());
    pindices.push_back(components.back()->getParameterIndex());
    extended = false;
  }


  if (extended) GET_FUNCTION_ADDR(ptr_to_AddPdfsExt);
  else GET_FUNCTION_ADDR(ptr_to_AddPdfs);

  initialise(pindices);
}


AddPdf::AddPdf (std::string n, std::vector<Variable*> weights, std::vector<PdfBase*> comps,unsigned int integr)
  : GooPdf(0, n)
  , extended(true)
{

  assert((weights.size() == comps.size()) || (weights.size() + 1 == comps.size()));

  // Indices stores (function index)(function parameter index)(weight index) triplet for each component.
  // Last component has no weight index unless function is extended.
  for (std::vector<PdfBase*>::iterator p = comps.begin(); p != comps.end(); ++p) {
    components.push_back(*p);
    assert(components.back());
  }

  getObservables(observables);

  std::vector<unsigned int> pindices;
  for (unsigned int w = 0; w < weights.size(); ++w) {
    assert(components[w]);
    pindices.push_back(components[w]->getFunctionIndex());
    pindices.push_back(components[w]->getFunctionPntIndex());
    pindices.push_back(components[w]->getParameterIndex());
    pindices.push_back(registerParameter(weights[w]));
  }
  assert(components.back());
  if (weights.size() < components.size()) {
    pindices.push_back(components.back()->getFunctionIndex());
    pindices.push_back(components.back()->getFunctionPntIndex());
    pindices.push_back(components.back()->getParameterIndex());
    extended = false;
  }


  if (extended){
        GET_FUNCTION_ADDR(ptr_to_AddPdfsExt_Bin_Mid);
        GET_INTEGRAL_ADDR(ptr_to_AddPdfsExt_Bin);
        GET_ATPOINTS_ADDR(ptr_to_AddPdfsExt_Bin_Mid);
        }
  else {GET_FUNCTION_ADDR(ptr_to_AddPdfs_Bin_Mid);
        GET_INTEGRAL_ADDR(ptr_to_AddPdfs_Bin);
        GET_ATPOINTS_ADDR(ptr_to_AddPdfs_Bin_Mid);}

  initialiseInt(pindices);
}



AddPdf::AddPdf (std::string n, Variable* frac1, PdfBase* func1, PdfBase* func2)
  : GooPdf(0, n)
  , extended(false)
{
  // Special-case constructor for common case of adding two functions.
  components.push_back(func1);
  components.push_back(func2);
  getObservables(observables);

  std::vector<unsigned int> pindices;
  pindices.push_back(func1->getFunctionIndex());
  pindices.push_back(func1->getParameterIndex());

  pindices.push_back(registerParameter(frac1));

  pindices.push_back(func2->getFunctionIndex());
  pindices.push_back(func2->getParameterIndex());

  GET_FUNCTION_ADDR(ptr_to_AddPdfs);

  initialise(pindices);
}

AddPdf::AddPdf (std::string n, Variable* frac1, PdfBase* func1, PdfBase* func2, unsigned int integr)
  : GooPdf(0, n)
  , extended(false)
{
  // Special-case constructor for common case of adding two functions.
  components.push_back(func1);
  components.push_back(func2);
  getObservables(observables);

  std::vector<unsigned int> pindices;
  pindices.push_back(func1->getFunctionIndex());
  pindices.push_back(func1->getFunctionPntIndex());
  pindices.push_back(func1->getParameterIndex());

  pindices.push_back(registerParameter(frac1));

  pindices.push_back(func2->getFunctionIndex());
  pindices.push_back(func2->getFunctionPntIndex());
  pindices.push_back(func2->getParameterIndex());

  GET_FUNCTION_ADDR(ptr_to_AddPdfs_Bin_Mid);
  GET_INTEGRAL_ADDR(ptr_to_AddPdfs_Bin);
  GET_ATPOINTS_ADDR(ptr_to_AddPdfs_Bin_Mid);

  initialiseInt(pindices);
}

__host__ fptype AddPdf::normalise () const {
  //if (cpuDebug & 1) std::cout << "Normalising AddPdf " << getName() << std::endl;

  fptype ret = 0;
  fptype totalWeight = 0;
  for (unsigned int i = 0; i < components.size()-1; ++i) {
    fptype weight = host_params[host_indices[parameters + 3*(i+1)]];
    totalWeight += weight;
    fptype curr = components[i]->normalise();
    ret += curr*weight;
  }
  fptype last = components.back()->normalise();
  if (extended) {
    fptype lastWeight = host_params[host_indices[parameters + 3*components.size()]];
    totalWeight += lastWeight;
    ret += last * lastWeight;
    ret /= totalWeight;
  }
  else {
    ret += (1 - totalWeight) * last;
  }
  host_normalisation[parameters] = 1.0;

  if (getSpecialMask() & PdfBase::ForceCommonNorm) {
    // Want to normalise this as
    // (f1 A + (1-f1) B) / int (f1 A + (1-f1) B)
    // instead of default
    // (f1 A / int A) + ((1-f1) B / int B).

    for (unsigned int i = 0; i < components.size(); ++i) {
      host_normalisation[components[i]->getParameterIndex()] = (1.0 / ret);
    }
  }

  //if (cpuDebug & 1) std::cout << getName() << " integral returning " << ret << std::endl;
  return ret;
}

__host__ double AddPdf::sumOfNll (int numVars) const {
  static thrust::plus<double> cudaPlus;
  thrust::constant_iterator<int> eventSize(numVars);
  thrust::constant_iterator<fptype*> arrayAddress(dev_event_array);
  double dummy = 0;
  thrust::counting_iterator<int> eventIndex(0);
  double ret = thrust::transform_reduce(thrust::make_zip_iterator(thrust::make_tuple(eventIndex, arrayAddress, eventSize)),
					thrust::make_zip_iterator(thrust::make_tuple(eventIndex + numEntries, arrayAddress, eventSize)),
					*logger, dummy, cudaPlus);

  if (extended) {
    fptype expEvents = 0;
    //std::cout << "Weights:";
    for (unsigned int i = 0; i < components.size(); ++i) {
      expEvents += host_params[host_indices[parameters + 3*(i+1)]];
      //std::cout << " " << host_params[host_indices[parameters + 3*(i+1)]];
    }
    // Log-likelihood of numEvents with expectation of exp is (-exp + numEvents*ln(exp) - ln(numEvents!)).
    // The last is constant, so we drop it; and then multiply by minus one to get the negative log-likelihood.
    ret += (expEvents - numEvents*log(expEvents));
    //std::cout << " " << expEvents << " " << numEvents << " " << (expEvents - numEvents*log(expEvents)) << std::endl;
  }

  return ret;
}
