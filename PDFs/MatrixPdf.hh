#ifndef MATRIX_PDF_HH
#define MATRIX_PDF_HH

#include "TComplex.h"
#include "TMath.h"


const fptype MLb = 5.61951;
const fptype MBd = 5.27961;
const fptype MPsi2S = 3.686109;
const fptype MJpsi = 3.096916;
const fptype MProton = 0.938272046;
const fptype MKaon = 0.493677;
const fptype MPion = 0.13957018;


const fptype MLb2 = MLb*MLb;
const fptype MLb4 = MLb2*MLb2;
const fptype MBd2 = MBd*MBd;
const fptype MBd4 = MBd2*MBd2;
const fptype MPsi2S2 = MPsi2S*MPsi2S;
const fptype MPsi2S4 = MPsi2S2*MPsi2S2;
const fptype MJpsi2 = MJpsi*MJpsi;
const fptype MJpsi4 = MJpsi2*MJpsi2;
const fptype MJpsi4mTwoMJpsi2MLb2pMLb4 = MJpsi4 - 2.*MJpsi2*MLb2 + MLb4;
const fptype MPsi2S4mTwoMPsi2S2MBd2pMBd4 = MPsi2S4 - 2.*MPsi2S2*MBd2 + MBd4;
const fptype MJpsi4mTwoMJpsi2MBd2pMBd4 = MJpsi4 - 2.*MJpsi2*MBd2 + MBd4;
const fptype TwoMJpsi2pTwoMLb2 = 2.*(MJpsi2 + MLb2);
const fptype TwoMPsi2S2pTwoMBd2 = 2.*(MPsi2S2 + MBd2);
const fptype TwoMJpsi2pTwoMBd2 = 2.*(MJpsi2 + MBd2);
const fptype InvTwoMLb = 1./(2.*MLb);
const fptype InvTwoMBd = 1./(2.*MBd);

const fptype MProton2 = MProton*MProton;
const fptype MProton4 = MProton2*MProton2;
const fptype MKaon2 = MKaon*MKaon;
const fptype MKaon4 = MKaon2*MKaon2;
const fptype MPion2 = MPion*MPion;
const fptype MPion4 = MPion2*MPion2;
const fptype MKaon4mTwoMKaon2MProton2pMProton4 = MKaon4 - 2.*MKaon2*MProton2 + MProton4;
const fptype MKaon4mTwoMKaon2MPion2pMPion4 = MKaon4 - 2.*MKaon2*MPion2 + MPion4;
const fptype TwoMKaon2pTwoMProton2 = 2.*(MKaon2 + MProton2);
const fptype TwoMKaon2pTwoMPion2 = 2.*(MKaon2 + MPion2);

// Lambda*
const fptype M1600 = 1.600 ;
const fptype G1600 = 0.150;
const fptype M1670 = 1.670;
const fptype G1670 = 0.035;
// K*
/*
const fptype M892 = 0.89581 ; // From PDG charged only K*(892)
const fptype G892 = 0.0474; // From PDG charged only K*(892)
const fptype M1410 = 1.414;
const fptype G1410 = 0.232;
const fptype M1430 = 1.425;
const fptype G1430 = 0.270;
*/

//const fptype dRad = 3.0;

class MatrixPdf : public GooPdf {
public:
  MatrixPdf (std::string n, Variable* _x);
  __host__ fptype integrate (fptype lo, fptype hi) const;
  __host__ virtual bool hasAnalyticIntegral () const {return true;}


private:

};

#endif
