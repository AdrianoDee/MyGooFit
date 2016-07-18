#ifndef MATRIX_PDF_HH
#define MATRIX_PDF_HH

#include "TComplex.h"
#include "TMath.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "devcomplex.hh"

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

/*
MEM_CONSTANT fptype MLb = 5.61951;
MEM_CONSTANT fptype MBd = 5.27961;
MEM_CONSTANT fptype MPsi2S = 3.686109;
MEM_CONSTANT fptype MJpsi = 3.096916;
MEM_CONSTANT fptype MProton = 0.938272046;
MEM_CONSTANT fptype MKaon = 0.493677;
MEM_CONSTANT fptype MPion = 0.13957018;*/

/*
MEM_CONSTANT fptype MLb2 = MLb*MLb;
MEM_CONSTANT fptype MLb4 = MLb2*MLb2;
MEM_CONSTANT fptype MBd2 = MBd*MBd;
MEM_CONSTANT fptype MBd4 = MBd2*MBd2;
MEM_CONSTANT fptype MPsi2S2 = MPsi2S*MPsi2S;
MEM_CONSTANT fptype MPsi2S4 = MPsi2S2*MPsi2S2;
MEM_CONSTANT fptype MJpsi2 = MJpsi*MJpsi;
MEM_CONSTANT fptype MJpsi4 = MJpsi2*MJpsi2;
MEM_CONSTANT fptype MJpsi4mTwoMJpsi2MLb2pMLb4 = MJpsi4 - 2.*MJpsi2*MLb2 + MLb4;
MEM_CONSTANT fptype MPsi2S4mTwoMPsi2S2MBd2pMBd4 = MPsi2S4 - 2.*MPsi2S2*MBd2 + MBd4;
MEM_CONSTANT fptype MJpsi4mTwoMJpsi2MBd2pMBd4 = MJpsi4 - 2.*MJpsi2*MBd2 + MBd4;
MEM_CONSTANT fptype TwoMJpsi2pTwoMLb2 = 2.*(MJpsi2 + MLb2);
MEM_CONSTANT fptype TwoMPsi2S2pTwoMBd2 = 2.*(MPsi2S2 + MBd2);
MEM_CONSTANT fptype TwoMJpsi2pTwoMBd2 = 2.*(MJpsi2 + MBd2);
MEM_CONSTANT fptype InvTwoMLb = 1./(2.*MLb);
MEM_CONSTANT fptype InvTwoMBd = 1./(2.*MBd);

MEM_CONSTANT fptype MProton2 = MProton*MProton;
MEM_CONSTANT fptype MProton4 = MProton2*MProton2;
MEM_CONSTANT fptype MKaon2 = MKaon*MKaon;
MEM_CONSTANT fptype MKaon4 = MKaon2*MKaon2;
MEM_CONSTANT fptype MPion2 = MPion*MPion;
MEM_CONSTANT fptype MPion4 = MPion2*MPion2;
MEM_CONSTANT fptype MKaon4mTwoMKaon2MProton2pMProton4 = MKaon4 - 2.*MKaon2*MProton2 + MProton4;
MEM_CONSTANT fptype MKaon4mTwoMKaon2MPion2pMPion4 = MKaon4 - 2.*MKaon2*MPion2 + MPion4;
MEM_CONSTANT fptype TwoMKaon2pTwoMProton2 = 2.*(MKaon2 + MProton2);
MEM_CONSTANT fptype TwoMKaon2pTwoMPion2 = 2.*(MKaon2 + MPion2);*/


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
    const std::vector< std::vector<fptype> >& _KstarDotSpin,
  	const vector<std::string>& _varNames,
  	const std::vector<Variable*>& _amplitudeGooVars,
  	const std::string& _psi_nS,
  	const fptype& _dRadB0, const fptype& _dRadKs);

  __host__ fptype integrate (fptype lo, fptype hi) const;
  __host__ virtual bool hasAnalyticIntegral () const {return true;}


private:
  //HOST SIDE
  std::vector< std::vector<fptype> > KstarDotSpin;
  vector< std::string > varNames;
  std::vector<Variable*> amplitudeGooVars;
  //map< TString,RooRealProxy* > amplitudeVarProxy_map;
  //std::map<std::string,Variable*> amplitudeVars_map;
  std::string psi_nS;
  fptype dRadB0, dRadKs;
  //DEVICE SIDE
  thrust::device_vector<thrust::device_vector<fptype> > d_KstarDotSpin;

};

#endif