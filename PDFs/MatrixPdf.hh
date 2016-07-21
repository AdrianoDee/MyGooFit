#ifndef MATRIX_PDF_HH
#define MATRIX_PDF_HH

#include "TComplex.h"
#include "TMath.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <iostream>

#include "devcomplex.hh"

#define KSTARSIZE 7
#define SKIP -12

/*
#define ZERO 0
#define M1 1
#define P1 2
#define M1H 3

#define K892_1 8921
#define H892_1_0 89210
#define H892_1_p1 89211
#define H892_1_m1 89212

#define K800_0 8000
#define H800_0_0 8000

#define K1410_1 14101
#define H1410_1_0 141010
#define H1410_1_p1 141011
#define H1410_1_m1 141012

#define K1430_0 14300
#define H1430_0_0 143000

#define K1430_2 14302
#define H1430_2_0 143020
#define H1430_2_p1 143020
#define H1430_2_m1 143020*/


MEM_CONSTANT fptype MLb = 5.61951;
MEM_CONSTANT fptype MBd = 5.27961;
MEM_CONSTANT fptype MPsi2S = 3.686109;
MEM_CONSTANT fptype MJpsi = 3.096916;
MEM_CONSTANT fptype MProton = 0.938272046;
MEM_CONSTANT fptype MKaon = 0.493677;
MEM_CONSTANT fptype MPion = 0.13957018;


 MEM_CONSTANT fptype  MLb2 = 31.5788926401;
 MEM_CONSTANT fptype  MLb4 = 997.226460374962;
 MEM_CONSTANT fptype  MBd2 = 27.8742817521;
 MEM_CONSTANT fptype  MBd4 = 776.975583195455;
 MEM_CONSTANT fptype  MPsi2S2 = 13.587399559881;
 MEM_CONSTANT fptype  MPsi2S4 = 184.6174267998544;
 MEM_CONSTANT fptype  MJpsi2 = 9.590888711055999;
 MEM_CONSTANT fptype  MJpsi4 = 91.98514626786141;
 MEM_CONSTANT fptype  MJpsi4mTwoMJpsi2MLb2pMLb4 = 483.4723167836545;
 MEM_CONSTANT fptype  MPsi2S4mTwoMPsi2S2MBd2pMBd4 = 204.1150027743444;
 MEM_CONSTANT fptype  MJpsi4mTwoMJpsi2MBd2pMBd4 = 334.2824610932961;
 MEM_CONSTANT fptype  TwoMJpsi2pTwoMLb2 = 82.33956270231201;
 MEM_CONSTANT fptype  TwoMPsi2S2pTwoMBd2 = 82.92336262396199;
 MEM_CONSTANT fptype  TwoMJpsi2pTwoMBd2 = 74.930340926312;
 MEM_CONSTANT fptype  InvTwoMLb = 0.08897572920058866;
 MEM_CONSTANT fptype  InvTwoMBd = 0.09470396487619351;
 MEM_CONSTANT fptype  MProton2 = 0.8803544323050262;
 MEM_CONSTANT fptype  MProton4 = 0.7750239264791049;
 MEM_CONSTANT fptype  MKaon2 = 0.243716980329;
 MEM_CONSTANT fptype  MKaon4 = 0.05939796650068616;
 MEM_CONSTANT fptype  MPion2 = 0.0194798351452324;
 MEM_CONSTANT fptype  MPion4 = 0.0003794639772854313;
 MEM_CONSTANT fptype  MKaon4mTwoMKaon2MProton2pMProton4 = 0.405307245258527;
 MEM_CONSTANT fptype  MKaon4mTwoMKaon2MPion2pMPion4 = 0.05028229728016605;
 MEM_CONSTANT fptype  TwoMKaon2pTwoMProton2 = 2.248142825268052;
 MEM_CONSTANT fptype  TwoMKaon2pTwoMPion2 = 0.5263936309484647;


// K*
const fptype M892 = 0.89581 ; const fptype G892 = 0.0474; // From PDG charged only K*(892)
const fptype M1410 = 1.414; const fptype G1410 = 0.232; // K*1410
const fptype M800 = 0.931; const fptype G800 = 0.578; // K*800
const fptype M1430_0 = 1.425; const fptype G1430_0 = 0.270; // K*1430_0
const fptype M1430_2 = 1.4324; const fptype G1430_2 = 0.109; // K*1430_2

const fptype K892_1 = 892.1;
const fptype K800_0 = 800.0;
const fptype K1410_1 = 1410.1;
const fptype K1430_0 = 1430.0;
const fptype K1430_2 = 1430.2;

class MatrixPdf : public GooPdf {
public:
  MatrixPdf(std::string n, Variable* _x, Variable* _cJ, Variable* _cKs, Variable* _phi,
  	const std::vector<Variable*>& _amplitudeGooVars,
    const std::host_vector<fptype>& _KStarVector,
  	fptype* _psi_nS,
  	fptype* _dRadB0, fptype* _dRadKs);

  __host__ fptype integrate (fptype lo, fptype hi) const;
  __host__ virtual bool hasAnalyticIntegral () const {return true;}

  /*
  EXEC_TARGET devcomplex<fptype> matrixElement(fptype mkp, fptype* p,unsigned int* indices,fptype helDmu);

  EXEC_TARGET devcomplex<fptype> RFunction(fptype mkp,fptype RMass, fptype RGamma, fptype MomMass, int LminMom, int LminR, fptype DB0, fptype DKs);
  EXEC_TARGET devcomplex<fptype> AngularTerm(fptype* p,unsigned int* indices, fptype spinR, fptype helJ, fptype helDmu,int iKStar);
  EXEC_TARGET fptype BlattWeisskopf(int Lmin, fptype q, fptype q0, fptype D);
  EXEC_TARGET fptype BWGamma(fptype mkp,fptype RMass, fptype RGamma, int Lmin, fptype D);
  EXEC_TARGET devcomplex<fptype> BW(fptype mkp,fptype RMass, fptype RGamma, int Lmin, fptype D);
  EXEC_TARGET devcomplex<fptype> H(fptype* p,unsigned int* indices, fptype helJ,int iKStar);
  EXEC_TARGET fptype Pmom(fptype mkp);
  EXEC_TARGET fptype Qmom(fptype mkp);
  EXEC_TARGET fptype PhiPHSP(fptype mkp);
  EXEC_TARGET fptype ME2();

  EXEC_TARGET fptype Wignerd_R(fptype spinR, fptype helJ, fptype cKs);
  EXEC_TARGET devcomplex<fptype> WignerD_J(fptype helJ, fptype helDmu, fptype angle,fptype cJ);*/


private:
  //HOST SIDE
  std::host_vector< fptype > KStarVector;
  //map< TString,RooRealProxy* > amplitudeVarProxy_map;
  //std::map<std::string,Variable*> amplitudeVars_map;
  fptype* psi_nS;
  fptype* dRadB0;
  fptype* dRadKs;
  int numberOfKStar;

};

#endif
