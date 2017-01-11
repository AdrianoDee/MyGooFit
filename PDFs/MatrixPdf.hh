#ifndef MATRIX_PDF_HH
#define MATRIX_PDF_HH

#include "GooPdf.hh"
#include "BinnedDataSet.hh"

#include "/lustrehome/cristella/work/Z_analysis/AA_fit/constants.h"
/*
//const fptype MBd = 5.27962;      // from PDG: http://pdglive.lbl.gov/Particle.action?node=S042&init=
const fptype MBd = 5.2794;         // from EvtGen: https://github.com/cms-sw/cmssw/blob/86088d61d757bad0e55addda14859c5ea6108d84/GeneratorInterface/ExternalDecays/data/evt.pdl#L79
//const fptype MPsi2S = 3.686097; // from PDG: http://pdglive.lbl.gov/Particle.action?node=M071&init=
const fptype MPsi2S = 3.68596; // from EvtGen https://github.com/cms-sw/cmssw/blob/86088d61d757bad0e55addda14859c5ea6108d84/GeneratorInterface/ExternalDecays/data/evt.pdl#L123
//const fptype MJpsi = 3.0969; // from PDG: http://pdglive.lbl.gov/Particle.action?node=M070&init=
const fptype MJpsi = 3.09687; // from EvtGen: https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/GeneratorInterface/ExternalDecays/data/evt.pdl#L37
//const fptype MKaon = 0.493677;   // from PDG: http://pdglive.lbl.gov/Particle.action?node=S010&init=
const fptype MKaon = 0.493677;     // from EvtGen: https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/GeneratorInterface/ExternalDecays/data/evt.pdl#L13
//const fptype MPion = 0.13957018; // from PDG: http://pdglive.lbl.gov/Particle.action?node=S008&init=
const fptype MPion = 0.13957;      // from EvtGen: https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/GeneratorInterface/ExternalDecays/data/evt.pdl#L9

const fptype MBd2 = MBd*MBd;
const fptype MBd4 = MBd2*MBd2;
const fptype MPsi2S2 = MPsi2S*MPsi2S;
const fptype MPsi2S4 = MPsi2S2*MPsi2S2;
const fptype MJpsi2 = MJpsi*MJpsi;
const fptype MJpsi4 = MJpsi2*MJpsi2;
const fptype MKaon2 = MKaon*MKaon;
const fptype MKaon4 = MKaon2*MKaon2;
const fptype MPion2 = MPion*MPion;
const fptype MPion4 = MPion2*MPion2;
*/
const fptype MPsi2S4mTwoMPsi2S2MBd2pMBd4 = MPsi2S4 - 2.*MPsi2S2*MBd2 + MBd4;
const fptype MJpsi4mTwoMJpsi2MBd2pMBd4 = MJpsi4 - 2.*MJpsi2*MBd2 + MBd4;
const fptype TwoMPsi2S2pTwoMBd2 = 2.*(MPsi2S2 + MBd2);
const fptype TwoMJpsi2pTwoMBd2 = 2.*(MJpsi2 + MBd2);
const fptype InvTwoMBd = 1./(2.*MBd);

const fptype MKaon4mTwoMKaon2MPion2pMPion4 = MKaon4 - 2.*MKaon2*MPion2 + MPion4;
const fptype TwoMKaon2pTwoMPion2 = 2.*(MKaon2 + MPion2);

#define KSTARSIZE 7
#define SKIP -12

#define M1HEL 123
#define P1HEL 234
#define ZEROHEL 0101

#define ZEROSPIN 0001
#define ONESPIN 1001
#define TWOSPIN 2001
#define THREESPIN 3001

#define PSIONE 1000
#define PSITWO 2000

class MatrixPdf : public GooPdf {
public:
  /*MatrixPdf(std::string n, Variable* _x, Variable* _cJ, Variable* _cKs, Variable* _phi,
    std::vector<Variable*>& _KParameters,
    Variable* _psi_nS, Variable* _dRadB0, Variable* _dRadKs);*/
  //MatrixPdf(std::string n, Variable* _x, Variable* _cJ, Variable* _cKs, Variable* _phi, Variable* _Mass,Variable* _Gamma,Variable* _Spin,Variable* _a,Variable* _b, Variable* _psi_nS, Variable* _dRadB0, Variable* _dRadKs);
  MatrixPdf(std::string n, Variable* _x, Variable* _mJP, Variable* _cJ, Variable* _phi, Variable* _B0beauty,
	    std::vector<Variable*> _Masses,std::vector<Variable*> _Gamma,std::vector<Variable*> _Spin,std::vector<Variable*> _a,std::vector<Variable*> _b,
	    Variable* _psi_nS, Variable* _dRadB0, Variable* _dRadKs);
  MatrixPdf(std::string n, Variable* _x, Variable* _mJP,Variable* _cJ, Variable* _phi,
	    std::vector<Variable*> _Masses,std::vector<Variable*> _Gamma,std::vector<Variable*> _Spin,std::vector<Variable*> _a,std::vector<Variable*> _b,
	    Variable* _psi_nS, Variable* _dRadB0, Variable* _dRadKs, BinnedDataSet* x);
  __host__ virtual bool hasAnalyticIntegral () const {return false;}

  //__host__ fptype integrate (fptype lo, fptype hi) const;

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
  EXEC_TARGET devcomplex<fptype> WignerD_J(fptype helJ, fptype helDmu, fptype angle,fptype cJ);
  */


private:
  //HOST SIDE
  Variable* psi_nS;
  Variable* dRadB0;
  Variable* dRadKs;

  thrust::device_vector<fptype>* dev_base_matrix_histo;
  thrust::device_vector<fptype>* dev_lowerlimits;
  thrust::device_vector<fptype>* dev_upperlimits;
  thrust::device_vector<fptype>* dev_steps;
  thrust::device_vector<fptype>* dev_bins;

  fptype totalEvents;


};

#endif
