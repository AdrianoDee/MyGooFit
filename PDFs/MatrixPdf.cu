/*****************************************************************************
 * Project: GooFit                                                           *
 *                                                                           *
 * This code was autogenerated by                             *
 *                                                                           *
 * A simple AA PDF class by Ivan Heredia de la Cruz on 4/25/16.              *
 *****************************************************************************/


#include <math.h>
#include "TMath.h"

#include "MatrixPdf.hh"
#include "devcomplex.hh"

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


EXEC_TARGET fptype BlattWeisskopf(int Lmin, fptype q, fptype q0, fptype D)
{
    fptype Dq = D*q;
    fptype Dq0 = D*q0;
    fptype Dq2 = Dq*Dq;
    fptype Dq02 = Dq0*Dq0;
    fptype Dq4 = Dq2*Dq2;
    fptype Dq04 = Dq02*Dq02;
    fptype Dq6 = Dq4*Dq2;
    fptype Dq06 = Dq04*Dq02;
    fptype rootterm = -1;

    if (Lmin==0)
      return 1.;
    else if (Lmin==1)
      rootterm = (1. + Dq02) / (1. + Dq2) ;
    else if (Lmin==2)
      rootterm = (9. + 3*Dq02 +Dq04) / (9. + 3*Dq2 +Dq4) ;
    else if (Lmin==3)
      rootterm = (225. + 45*Dq02 + 6*Dq04 + Dq06) / (225. + 45*Dq2 + 6*Dq4 + Dq6) ;
    else { printf("WARNING! Lmin = %d is not implemented for BlattWeisskopf functions at the moment. Returning 1 -> \"AngularTerm\" = 0\n",Lmin);
      return 1.;
    }

    if (rootterm > 0)
        return SQRT( rootterm );
    else { printf("WARNING! In \"BlattWeisskopf\" function: rootterm <= 0 -> returning 1\n");
        return 1.;
    }
}


EXEC_TARGET fptype Pmom(fptype mkp,fptype psiN)
{

    fptype  MPsi2S4mTwoMPsi2S2MBd2pMBd4 = 204.1150027743444;
    fptype  MJpsi4mTwoMJpsi2MBd2pMBd4 = 334.2824610932961;

    fptype  TwoMPsi2S2pTwoMBd2 = 82.92336262396199;
    fptype  TwoMJpsi2pTwoMBd2 = 74.930340926312;

    fptype  InvTwoMBd = 0.09470396487619351;

    fptype mkp2 = mkp*mkp;
    fptype rootterm = 0;

    if (psiN==1.0)
      rootterm = MJpsi4mTwoMJpsi2MBd2pMBd4 + mkp2*(mkp2 - TwoMJpsi2pTwoMBd2);
    else if (psiN==2.0)
      rootterm = MPsi2S4mTwoMPsi2S2MBd2pMBd4 + mkp2*(mkp2 - TwoMPsi2S2pTwoMBd2);
    else
      //cout <<"psi_nS = " <<psi_nS <<" not allowed in \"Pmom\" function at the moment. Keeping rootterm at 0" <<endl;
      printf("psi_nS = %.2f not allowed in \"Pmom\" function at the moment. Keeping rootterm at 0\n",psiN);

    //cout <<"\nrootterm for psi_nS = " <<psi_nS <<" and mkp = " <<mkp <<": " <<rootterm <<endl;
    if (rootterm > 0)
        return InvTwoMBd * SQRT(rootterm);
    else { //cout <<"WARNING! In \"Pmom\" function: rootterm (" <<rootterm <<") <= 0 for mkp = " <<mkp <<" and psi_nS = " <<psi_nS <<" -> returning 0" <<endl;
           printf("WARNING! In \"Pmom\" function: rootterm (%.2f) <= 0 for mkp = %.2f and psi_nS = %.2f  -> returning 0 \n",rootterm,mkp,psiN);
           return 0.;
    }

}


EXEC_TARGET fptype Qmom(fptype mkp)
{

      fptype mkp2 = mkp*mkp;


     fptype  MKaon4mTwoMKaon2MPion2pMPion4 = 0.05028229728016605;

     fptype  TwoMKaon2pTwoMPion2 = 0.5263936309484647;

    fptype rootterm = MKaon4mTwoMKaon2MPion2pMPion4 + mkp2*(mkp2 - TwoMKaon2pTwoMPion2) ;
    if (rootterm > 0)
        return 0.5*SQRT(rootterm)/mkp;
    else { //cout <<"WARNING! In \"Qmom\" function: rootterm <= 0 for mkp = " <<mkp <<" -> returning 0" <<endl;
            printf("WARNING! In \"Qmom\" function: rootterm (%.2f) <= 0 for mkp = %.2f  -> returning 0 \n",rootterm,mkp);

        return 0.;
    }

}

EXEC_TARGET fptype BWGamma(fptype mkp,fptype RMass, fptype RGamma, int Lmin, fptype D)
{
    fptype QmKP = Qmom(mkp);
    fptype QRMass = Qmom(RMass);
    int expoterm = 2*Lmin + 1 ;

    fptype BWG = ( RGamma * RMass * POW(QmKP/QRMass,expoterm) * POW(BlattWeisskopf(Lmin, QmKP, QRMass, D),2) ) / mkp;
    //cout <<"BWGamma for RMass = " <<RMass <<": " <<BWG <<endl;
    return BWG ;

}

EXEC_TARGET devcomplex<fptype> BW(fptype mkp,fptype RMass, fptype RGamma, int Lmin, fptype D)
{

    fptype num1term = RMass*RMass - mkp*mkp ;
    fptype num2term = RMass * BWGamma(mkp,RMass, RGamma, Lmin, D) ;
    fptype denoterm = num1term*num1term + num2term*num2term ;

    devcomplex<fptype> bw (num1term / denoterm, num2term / denoterm);
    //cout <<"BW for RMass = " <<RMass <<": " <<bw <<endl;
    return bw ;

}

EXEC_TARGET devcomplex<fptype> H(fptype* p,unsigned int* indices, fptype helJ,int iKStar)
{

  devcomplex<fptype> result(0.0,0.0);
  devcomplex<fptype> imUnit(0.0,1.0);

  int whichOfThree = 0;

  if(helJ==P1HEL) whichOfThree = 1;
  if(helJ==M1HEL) whichOfThree = 2;

  //fptype a = p[indices[3+(iKStar+whichOfThree)*6+4]];
  //fptype b = p[indices[3+(iKStar+whichOfThree)*6+5]];

  fptype a = p[indices[7]];
  fptype b = p[indices[8]];

  result = exp(imUnit*b);

  return a * result ;

}

EXEC_TARGET fptype Wignerd_R(fptype spinR, fptype helJ, fptype cKs)
{

  if (spinR==0.0)
    return 1. ;
  else if (spinR==1.0)
    if (helJ==M1HEL)
      return +(SIN(ACOS(cKs)) / root2) ;
    else if (helJ==ZEROHEL)
      return cKs ;
    else if (helJ==P1HEL)
      return -(SIN(ACOS(cKs)) / root2) ;
    else {
      printf("Wignerd_R Spin 1.0 PRINFT TO BE CONFIGURED returning 0\n");
      //cout <<"helJ = " <<helJ <<" is not allowed for spinR-" <<spinR <<" Wigner d^{spinR}_{helJ,0} functions. Returning 0" <<endl;
      return 0 ;
    }
  else if (spinR==2.0)
    if (helJ==M1HEL)
      return +(SIN(2*ACOS(cKs)) * SQRT(3./8.)) ;
    else if (helJ==ZEROHEL)
      return +(3*POW(cKs,2) -1)/2. ;
    else if (helJ==P1HEL)
      return -(SIN(2*ACOS(cKs)) * SQRT(3./8.)) ;
    else {
      printf("Wignerd_R Spin 2.0 PRINFT TO BE CONFIGURED returning 0\n");
      //cout <<"helJ = " <<helJ <<" is not allowed for spinR-" <<spinR <<" Wigner d^{spinR}_{helJ,0} functions. Returning 0" <<endl;
      return 0 ;
    }
  else if (spinR==3.0)
    if (helJ==M1HEL)
      return +(SIN(ACOS(cKs))*(5*COS(2*ACOS(cKs)) + 3.) * SQRT(3.)/8.) ;
    else if (helJ==ZEROHEL)
      return +(5*POW(cKs,3) - 3*cKs)/2. ;
    else if (helJ==P1HEL)
      return -(SIN(ACOS(cKs))*(5*COS(2*ACOS(cKs)) + 3.) * SQRT(3.)/8.) ;
    else {
      printf("Wignerd_R Spin 3.0 PRINFT TO BE CONFIGURED returning 0\n");
      //cout <<"helJ = " <<helJ <<" is not allowed for spinR-" <<spinR <<" Wigner d^{spinR}_{helJ,0} functions. Returning 0" <<endl;
      return 0 ;
    }
  else {
    printf("PRINFT TO BE CONFIGURED returning 0\n");
    //cout <<"spinR = " <<spinR <<" is not implemented for Wigner d^{spinR}_{helJ,0} functions at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl;
    return 0 ;
  }
}

EXEC_TARGET devcomplex<fptype> WignerD_J(fptype helJ, fptype helDmu, fptype angle,fptype cJ)
{

  devcomplex<fptype> imUnit(0.0,1.0);
  devcomplex<fptype> nullo(0.0,0.0);

  if (helJ==M1HEL) {
    if (helDmu==M1HEL)
      return ((+1. + cJ)*exp(imUnit*angle))*.5;
    else if (helDmu==P1HEL)
      return (-1.0)*((-1. + cJ)*exp(imUnit*angle))*.5;
    else {
      printf("PRINFT TO BE CONFIGURED returning 0\n");
      //cout <<"helDmu = " <<helDmu <<" not allowed in \"WignerD_J\" functions for helJ = " <<helJ <<" at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
      return nullo; }
  } else if (helJ==ZEROHEL) {
    if (helDmu==M1HEL)
      return devcomplex<fptype>((-1.0)*(SQRT(1. - POW(cJ,2))/root2),0.0);
    else if (helDmu==P1HEL)
      return devcomplex<fptype>((SQRT(1. - POW(cJ,2))/root2),0.0);
    else {
      printf("PRINFT TO BE CONFIGURED returning 0\n");
      //cout <<"helDmu = " <<helDmu <<" not allowed in \"WignerD_J\" functions for helJ = " <<helJ <<" at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
      return nullo; }
  } else if(helJ==P1HEL) {
    if (helDmu==M1HEL)
      return (-1.0)*(-1. + cJ)/(2.*exp(imUnit*angle));
    else if (helDmu==P1HEL)
      return (+1. + cJ)/(2.*exp(imUnit*angle));
    else {
      printf("WignerD_J HelDem PRINFT TO BE CONFIGURED returning 0\n");
      //cout <<"helDmu = " <<helDmu <<" not allowed in \"WignerD_J\" functions for helJ = " <<helJ <<" at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
      return nullo; }
  } else {
    printf("PRINFT TO BE CONFIGURED returning 0\n");
    //cout <<"helJ = " <<helJ <<" not allowed in \"WignerD_J\" functions at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
    return nullo;
  }

}


EXEC_TARGET devcomplex<fptype> AngularTerm(fptype cJ, fptype cKs, fptype phi, fptype* p,unsigned int* indices,fptype spinR, fptype helJ, fptype helDmu,int iKStar)
{
  devcomplex<fptype> result;

  result = H(p,indices,helJ,iKStar) * Wignerd_R(spinR, helJ,cKs) * conj( WignerD_J(helJ, helDmu, phi,cJ) ) ;
  //cout <<"\nAngularTerm for K* " <<R <<" and helDmu = " <<helDmu <<" and helJ = " <<helJ <<" is made of Wignerd_R(spinR, helJ) * cWignerD_J(helJ, helDmu, phi) = " <<Wignerd_R(spinR, helJ) <<" * " <<cWignerD_J( WignerD_J(helJ, helDmu, phi) ) <<endl;
  //cout <<"It is multiplied by H(R,helJ) = H(" <<R <<"," <<helJ <<") = " <<H(R,helJ) <<endl;
  printf("====== AngularTerm = (%.3f , %.3f) cJ = %.3f cKs = %.3f phi = %.3f \n",cJ,cKs,phi);

  return result;
}

EXEC_TARGET devcomplex<fptype> RFunction(fptype mkp,fptype RMass, fptype RGamma, fptype MomMass, int LminMom, int LminR, fptype DB0, fptype DKs,fptype psi_nS)
{
    fptype PmKP = Pmom(mkp,psi_nS);
    fptype PRMass = Pmom(RMass,psi_nS);
    fptype QmKP = Qmom(mkp);
    fptype QRMass = Qmom(RMass);

    fptype blatt1 = BlattWeisskopf(LminMom, PmKP, PRMass, DB0);
    fptype blatt2 = BlattWeisskopf(LminR, QmKP, QRMass, DKs);
    devcomplex<fptype> bw = BW(mkp,RMass, RGamma, LminR, DKs);

    fptype pow1 = POW(PmKP/MomMass,LminMom);
    fptype pow2 = POW(QmKP/mkp,LminR);

    //TComplex RFunc = BlattWeisskopf(LminMom, PmKP, PRMass, D) * POW(PmKP/MomMass,LminMom) * BW(RMass, RGamma, LminR, D) * BlattWeisskopf(LminR, QmKP, QRMass, D) * POW(QmKP/RMass,LminR);
    devcomplex<fptype> RFunc = blatt1 * pow1 * bw * blatt2 * pow2;
    //cout <<"BlattWeisskopf(LminR, QmKP, QRMass, D) for RMass " <<RMass <<" = " <<BlattWeisskopf(LminR, QmKP, QRMass, D) <<endl;
    //cout <<"BlattWeisskopf(LminMom, PmKP, PRMass, D) for RMass " <<RMass <<" = " <<BlattWeisskopf(LminMom, PmKP, PRMass, D) <<endl;
    //cout <<"BlattWeisskopf(LminMom, PmKP, PRMass, D) * BlattWeisskopf(LminR, QmKP, QRMass, D) for RMass " <<RMass <<" = " <<BlattWeisskopf(LminMom, PmKP, PRMass, D) * BlattWeisskopf(LminR, QmKP, QRMass, D) <<endl;
    //cout <<"POW(QmKP/RMass,LminR) for RMass " <<RMass <<" = " <<POW(QmKP/mKP,LminR) <<endl;
    //cout <<"POW(PmKP/MomMass,LminMom) * POW(QmKP/RMass,LminR) for RMass " <<RMass <<" = " <<(POW(PmKP/MomMass,LminMom) * POW(QmKP/RMass,LminR)) <<endl;
    //cout <<"\nRFunction for RMass " <<RMass <<" = " <<RFunc <<"\n\n" <<endl;

    printf("\n RFunction (%.3f ,%.3f) for RMass = %.3f Mkp = %.3f  \n PmKP  = %.3f PRMass = %.3f    \n QmKP  = %.3f QRMass = %.3f    \n QmKP = %.3f LminMom = %.3f    \n LminMom = %.3f LminR = %.3f    \n DB0 = %.3f DKs = %.3f    \n BlattWeisskopf LminM = %.3f    \n BlattWeisskopf LminR = %.3f    \n Power1  = %.3f    \n Power2  = %.3f    \n BW  = (%.3f ,%.3f)",RFunc.real,RFunc.imag,RMass,mkp,PmKP,PRMass,QmKP,QRMass,LminMom,LminR,DB0,DKs,blatt1,blatt2,pow1,pow2,bw.real,bw.imag);

    //printf("======= RFunction (%.2f ,%.2f)\n  RMass = %.3f Mass KPi = %.3f ",RFunc.real,RFunc.imag,RMass,mkp);

    return RFunc ;
}

EXEC_TARGET devcomplex<fptype> matrixElement(fptype mkp, fptype cJ, fptype cKs, fptype phi, fptype* p,unsigned int* indices,fptype helDmu)
{

  int numParams = indices[0];
  fptype MBd = 5.27961;
  //int numberOfKStar = indices[0]/6;


  fptype psi_nS = p[indices[1]];
  fptype dRadB0 = p[indices[2]];
  fptype dRadKs = p[indices[3]];

  devcomplex<fptype> matrixElement (0.0,0.0);

  int iKStar = 0;

  // K+ and pi- have 0 spin -> second last argument of K* RFunction is = spin(K*)

  /*for (int iKStar=0; iKStar<numParams-3; iKStar += 6) {

    fptype Mass = p[indices[3+1+iKStar]];
    fptype Spin = p[indices[3+2+iKStar]];
    fptype Gamma = p[indices[3+3+iKStar]];

    printf("Mass = %f Gamma = %f Spin = %f \n",Mass,Gamma,Spin);

    devcomplex<fptype> matrixElement_R(0.0,0.0);
    if (Spin==0.0) { // for spin0 K*, third last argument = spin(psi_nS) = spin.Atoi() + 1 = 1
      matrixElement_R = RFunction(mkp,Mass,Gamma, MBd, Spin+1, Spin, dRadB0, dRadKs,psi_nS) *
	               AngularTerm(cJ,cKs,phi,p,indices,Spin, ZEROHEL, helDmu,iKStar) ;
    } else
              {
                matrixElement_R = RFunction(mkp,Mass,Gamma, MBd, Spin-1, Spin,dRadB0,dRadKs,psi_nS) *
	               ( AngularTerm(cJ,cKs,phi,p,indices,Spin, M1HEL, helDmu,iKStar) + AngularTerm(cJ,cKs,phi,p,indices, Spin, ZEROHEL, helDmu,iKStar) + AngularTerm(cJ,cKs,phi,p,indices,Spin, P1HEL, helDmu,iKStar) ) ;
                   iKStar +=6*2;
               }
    //cout <<"\nAngularTerm.Rho() for " <<R <<" = " <<(AngularTerm(R, spin, "0", helDmu)).Rho() <<endl;
    //cout <<"matrixElement for (R,helDmu) = (" <<R <<"," <<helDmu <<") = H(R,helJ) * RFunction * AngularTerm = " <<matrixElement_R <<endl;
    matrixElement += matrixElement_R;
    //cout <<"matrixElement_R.Rho2() for (R,helDmu) = (" <<R <<"," <<helDmu <<") = " <<matrixElement_R.Rho2() <<"\n\n" <<endl;
  }*/

  fptype Mass = p[indices[4]];
  fptype Gamma = p[indices[5]];
  fptype Spin = p[indices[6]];

  //printf("Mass = %f Gamma = %f Spin = %f psi_nS = %f dRadB0 = %f dRadKs = %f \n",Mass,Gamma,Spin,psi_nS,dRadB0,dRadKs);

  devcomplex<fptype> matrixElement_R(0.0,0.0);
  if (Spin==0.0) { // for spin0 K*, third last argument = spin(psi_nS) = spin.Atoi() + 1 = 1
    matrixElement_R = RFunction(mkp,Mass,Gamma, MBd, Spin+1, Spin, dRadB0, dRadKs,psi_nS) *
               AngularTerm(cJ,cKs,phi,p,indices,Spin, ZEROHEL, helDmu,iKStar) ;
  } else
            {
              matrixElement_R = RFunction(mkp,Mass,Gamma, MBd, Spin-1, Spin,dRadB0,dRadKs,psi_nS) *
               ( AngularTerm(cJ,cKs,phi,p,indices,Spin, M1HEL, helDmu,iKStar) + AngularTerm(cJ,cKs,phi,p,indices, Spin, ZEROHEL, helDmu,iKStar) + AngularTerm(cJ,cKs,phi,p,indices,Spin, P1HEL, helDmu,iKStar) ) ;
                 iKStar +=6*2;
             }
  //cout <<"\nAngularTerm.Rho() for " <<R <<" = " <<(AngularTerm(R, spin, "0", helDmu)).Rho() <<endl;
  //cout <<"matrixElement for (R,helDmu) = (" <<R <<"," <<helDmu <<") = H(R,helJ) * RFunction * AngularTerm = " <<matrixElement_R <<endl;
  matrixElement += matrixElement_R;
  //cout <<"matrixElement_R.Rho2() for (R,helDmu) = (" <<R <<"," <<helDmu <<") = " <<matrixElement_R.Rho2() <<"\n\n" <<endl;

  printf("======= Matrix Element HEL = %.2f \n  Mass KPi = %.3f cJ = %.3f  cKs = %.3f phi = %.3f \n mE = ( %.3f , %.3f )",helDmu,mkp,cJ,cKs,phi,matrixElement.real,matrixElement.imag);

  return matrixElement;

}

EXEC_TARGET fptype ME2(fptype mkp, fptype cJ, fptype cKs, fptype phi, fptype* p,unsigned int* indices)
{
  //cout <<"\nME(\"m1\") + ME(\"p1\") = " <<ME("m1") <<" + " <<ME("p1") <<endl;
  //cout <<"ME(\"m1\").Rho2() + ME(\"p1\").Rho2() = " <<ME("m1").Rho2() <<" + " <<ME("p1").Rho2() <<endl;
  return matrixElement(mkp,cJ,cKs,phi,p,indices,M1HEL).abs2() + matrixElement(mkp,cJ,cKs,phi,p,indices,P1HEL).abs2() ;
}

EXEC_TARGET fptype PhiPHSP(fptype mkp,fptype psiN)
{
    //printf("=======Phase Space \n");
    fptype p = Pmom(mkp,psiN);
    fptype q = Qmom(mkp);
    fptype phsp = p*q;
    //printf(" Mass KPi = %.3f Phase space = %.3f\n",mkp,phsp);
    //printf("==================");

    return  Pmom(mkp,psiN) * Qmom(mkp);
}

EXEC_TARGET fptype device_Matrix (fptype* evt, fptype* p, unsigned int* indices) {

  int numParams = indices[0];

  fptype mkp = evt[indices[2 + indices[0]]];
  fptype cJ = evt[indices[2 + indices[0]]+1];
  fptype cKs = evt[indices[2 + indices[0]]+2];
  fptype phi = evt[indices[2 + indices[0]]+3];

  fptype psi_nS = p[indices[1]];
  fptype dRadB0 = p[indices[2]];
  fptype dRadKs = p[indices[3]];
  fptype MKaon = 0.493677; fptype MPion = 0.13957018;
  fptype MBd = 5.27961;

  //printf("%.2f %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f ",psi_nS,dRadB0,dRadKs,p[indices[4]],p[indices[5]],p[indices[6]],p[indices[7]],p[indices[8]]);



  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
   fptype MPsi_nS = 0.;

   if (psi_nS==1.0)
     MPsi_nS = 3.096916;
   else if (psi_nS==2.0)
     MPsi_nS = 3.686109;
   else
      printf("PRINFT TO BE CONFIGURED = 0 mpk = %.2f cJ = %.2f cKs = %.2f phi = %.2f psi_nS = %f \n",mkp,cJ,cKs,phi,psi_nS);
      //printf("mpk = %.2f (%.2f - %.2f) cJ = %.2f cKs = %.2f phi = %.2f \n",mkp,MBd - MPsi_nS,MKaon + MPion,cJ,cKs,phi);

  if ((mkp < MKaon + MPion) || (mkp > MBd - MPsi_nS)){
    printf("Out of the borders \n");
    return 0.;}
  else{
      //printf("Device Matrix mkp = %.2f cJ = %.2f cKs = %.2f phi = %.2f \n",mkp,cJ,cKs,phi);
      return ME2(mkp,cJ,cKs,phi,p,indices) * PhiPHSP(mkp,psi_nS);}


}

/*
EXEC_TARGET fptype device_Matrix_Point (fptype* point, fptype* p, unsigned int* indices) {

  return 0;

}

EXEC_TARGET fptype device_Matrix_Bin (fptype* point, fptype* p, unsigned int* indices) {

  return 0;

}*/


MEM_DEVICE device_function_ptr ptr_to_Matrix = device_Matrix;
//MEM_DEVICE device_function_ptr ptr_to_Matrix_Point = device_Matrix_Point;
//MEM_DEVICE device_function_ptr ptr_to_Matrix_Bin = device_Matrix_Bin;
/*
__host__ MatrixPdf::MatrixPdf (std::string n, Variable* _x, Variable* _cJ, Variable* _cKs, Variable* _phi,
  std::vector<Variable*>& _KParameters,
  Variable* _psi_nS, Variable* _dRadB0, Variable* _dRadKs)*/
__host__ MatrixPdf::MatrixPdf (std::string n, Variable* _x, Variable* _cJ, Variable* _cKs, Variable* _phi,
    Variable* _Mass,Variable* _Gamma,Variable* _Spin,Variable* _a,Variable* _b,
    Variable* _psi_nS, Variable* _dRadB0, Variable* _dRadKs)
  : GooPdf(0, n),
  psi_nS(_psi_nS),dRadB0(_dRadB0),dRadKs(_dRadKs)
{

  registerObservable(_x);
  registerObservable(_cJ);
  registerObservable(_cKs);
  registerObservable(_phi);

  std::vector<unsigned int> pindices;

  pindices.push_back(registerParameter(_psi_nS));  // p[indices[1]]
  pindices.push_back(registerParameter(_dRadB0));  // p[indices[2]]
  pindices.push_back(registerParameter(_dRadKs));  // p[indices[3]]

  pindices.push_back(registerParameter(_Mass));  // p[indices[4]]
  pindices.push_back(registerParameter(_Gamma));  // p[indices[5]]
  pindices.push_back(registerParameter(_Spin));  // p[indices[6]]
  pindices.push_back(registerParameter(_a));  // p[indices[7]]
  pindices.push_back(registerParameter(_b));  // p[indices[8]]

  /*for (int j = 0 ; j < (int)_KParameters.size(); j++) {
    pindices.push_back(registerParameter(_KParameters[j]));
  }*/

  /*
  for (int j = 0 ; j < (int)_amplitudeGooVars.size(); j++) {
    pindices.push_back(registerParameter(_amplitudeGooVars[j]));
  }*/
  /*
  gooMalloc((void**) d_psi_nS,sizeof(fptype));
  gooMalloc((void**) d_dRadB0,sizeof(fptype));
  gooMalloc((void**) d_dRadKs,sizeof(fptype));

  MEMCPY(d_psi_nS,psi_nS,sizeof(fptype),cudaMemcpyHostToDevice);
  MEMCPY(d_dRadB0,dRadB0,sizeof(fptype),cudaMemcpyHostToDevice);
  MEMCPY(d_dRadKs,dRadKs,sizeof(fptype),cudaMemcpyHostToDevice);

  gooMalloc((void**) d_numberOfKStar,sizeof(int));
  MEMCPY(d_numberOfKStar,&numberOfKStar,sizeof(int),cudaMemcpyHostToDevice);

  gooMalloc((void**) d_KStarVector,sizeof(fptype)*(int)KStarVector.size());
  MEMCPY(d_KStarVector,&KStarVector[0],sizeof(int)*KStarVector.size(),cudaMemcpyHostToDevice);
  */
  GET_FUNCTION_ADDR(ptr_to_Matrix);
  //GET_INTEGRAL_ADDR(ptr_to_Matrix_Bin);
  //GET_ATPOINTS_ADDR(ptr_to_Matrix_Point);
  initialise(pindices);
}
