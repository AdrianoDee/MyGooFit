/*****************************************************************************
 * Project: GooFit                                                           *
 *                                                                           *
 * This code was autogenerated by                             *
 *                                                                           *
 * A simple AA PDF class by Ivan Heredia de la Cruz on 4/25/16.              *
 *****************************************************************************/


#include <math.h>
//#include "TMath.h"

#include "MatrixPdf.hh"
#include "BinnedDataSet.hh"
#include "devcomplex.hh"

//#define MDEBUGGING 1

EXEC_TARGET fptype cosTheta_FromMasses(const fptype sameSideM2, const fptype oppositeSideM2, const fptype psi_nSM2, const fptype motherM2, const fptype refM2, const fptype otherM2) {

  fptype num = (sameSideM2/2)*(motherM2 + refM2 - oppositeSideM2) - (1./4.)*(motherM2 - psi_nSM2 + sameSideM2)*(sameSideM2 - otherM2 + refM2) ;
  fptype denom2 = ((1./4.)*POW(motherM2 - psi_nSM2 + sameSideM2,2) - sameSideM2*motherM2) * ((1./4.)*POW(sameSideM2 - otherM2 + refM2,2) - sameSideM2*refM2) ;

  return (num / SQRT(denom2)) ;

}

EXEC_TARGET bool dalitz_contour_dev(const fptype mKP, const fptype mPsiP, const bool massSquared, const fptype psi_nS) {

  fptype MPsi_nS = 0;

  printf("Npsi problems? %.2f %d\n",psi_nS,psi_nS);

  if (psi_nS == 1.0)
  {
    printf("One? %.2f %d\n",psi_nS,psi_nS);
    MPsi_nS = MJpsi;
  }
  else if (psi_nS == 2.0)
    {
      printf("Two? %.2f %d\n",psi_nS,psi_nS);
      MPsi_nS = MPsi2S;
    }
  else
  {
    printf("Mmm nspsi problem %.2f %d\n",psi_nS,psi_nS);
    return false;
  }
  printf("No mpsipi problems \n");
  fptype mKP_1 = mKP;
  fptype mPsiP_1 = mPsiP;

  if (massSquared) {
    mKP_1 = SQRT( mKP );
    mPsiP_1 = SQRT( mPsiP );
  }

  if ((mKP_1 < MKaon + MPion) || (mKP_1 > MBd - MPsi_nS) || (mPsiP_1 < MPsi_nS + MPion) || (mPsiP_1 > MBd - MKaon))
    return false;
  else { // Dalitz border from PDG KINEMATICS 43.4.3.1.
    fptype E_P = (mPsiP_1*mPsiP_1 - MJpsi2 + MPion2)/(2*mPsiP_1) ;
    fptype E_K = (MBd2 - mPsiP_1*mPsiP_1 - MKaon2)/(2*mPsiP_1) ;
    fptype E_PpE_K_2 = POW((E_P + E_K),2.);
    fptype sqrt_E_P2mMP2 = SQRT(E_P*E_P - MPion2);
    fptype sqrt_E_K2mMK2 = SQRT(E_K*E_K - MKaon2);
    fptype mKP2_min = E_PpE_K_2 - POW(sqrt_E_P2mMP2 + sqrt_E_K2mMK2,2.);
    fptype mKP2_max = E_PpE_K_2 - POW(sqrt_E_P2mMP2 - sqrt_E_K2mMK2,2.);
    if ((mKP_1*mKP_1 < mKP2_min) || (mKP_1*mKP_1 > mKP2_max))
      return true;
  }

  return false ;

}

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


EXEC_TARGET fptype Pmom(fptype mkp,const fptype psiN)
{

    fptype mkp2 = mkp*mkp;
    fptype rootterm = 0;

    //printf("=======Pmom mpk = %.2f psiN = %.2f \n",mkp,psiN);

    //if (psiN==1.0)
    if (true)
      rootterm = MJpsi4mTwoMJpsi2MBd2pMBd4 + mkp2*(mkp2 - TwoMJpsi2pTwoMBd2);
    else if (psiN==2.0)
      rootterm = MPsi2S4mTwoMPsi2S2MBd2pMBd4 + mkp2*(mkp2 - TwoMPsi2S2pTwoMBd2);
    else
      {printf("psi_nS = %.2f not allowed in \"Pmom\" function at the moment. Keeping rootterm at 0\n",psiN);}

    if (rootterm > 0){
        return InvTwoMBd * SQRT(rootterm);}
    else { //cout <<"WARNING! In \"Pmom\" function: rootterm (" <<rootterm <<") <= 0 for mkp = " <<mkp <<" and psi_nS = " <<psi_nS <<" -> returning 0" <<endl;
           printf("WARNING! In \"Pmom\" function: rootterm (%.2f) <= 0 for mkp = %.2f and psi_nS = %.2f  -> returning 0 \n",rootterm,mkp,psiN);
           return 0.;
    }

}


EXEC_TARGET fptype Qmom(fptype mkp)
{

  fptype mkp2 = mkp*mkp;

  fptype rootterm = MKaon4mTwoMKaon2MPion2pMPion4 + mkp2*(mkp2 - TwoMKaon2pTwoMPion2) ;
  if (rootterm >= 0)
    return 0.5*SQRT(rootterm)/mkp;
  else {
    printf("WARNING! In \"Qmom\" function: rootterm (%f) < 0 for mkp = %f  -> returning 0 \n", rootterm, mkp);
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

EXEC_TARGET devcomplex<fptype> H(fptype* p,unsigned int* indices, fptype helJ,int lA)
{

  devcomplex<fptype> result(0.0,0.0);
  devcomplex<fptype> imUnit(0.0,1.0);

  int whichOfThree = 0;


  if(helJ==P1HEL) whichOfThree = 1;
  if(helJ==M1HEL) whichOfThree = 2;

  // psi_ns dRadB0 dRadKs m1 g1 s1 m2 g2 s2 . . . mn  gn  sn   a1          b1             a2 b2 . . . an                  bn                    //
  //Indeces
  // 2      3      4      5  6  7  8  9  10        4+n 5+n 6+n  4+nKstars*3 4+nKstars*3+1              4+nKstars*3+(n-1)*2 4+nKstars*3+(n-1)*2+1

  int noOfMasses = indices[1];
  fptype a = p[indices[5+noOfMasses*3+lA*2+whichOfThree*2]];
  fptype b = p[indices[5+noOfMasses*3+lA*2+whichOfThree*2+1]];

  #ifdef MDEBUGGING
  if(helJ==ZEROHEL) printf("Which of Three : %d Index : %d  a = %.3f  b = %.3f for helJ = 0 (%.2f) noOfMasses = %d lA = %d \n ",whichOfThree,5+noOfMasses*3+lA+whichOfThree*2,a,b,helJ,noOfMasses,lA);
  if(helJ==M1HEL) printf("Which of Three : %d Index : %d  a = %.3f  b = %.3f for helJ = M1 (%.2f) noOfMasses = %d lA = %d \n",whichOfThree,5+noOfMasses*3+lA+whichOfThree*2,a,b,helJ,noOfMasses,lA);
  if(helJ==P1HEL) printf("Which of Three : %d Index : %d  a = %.3f  b = %.3f for helJ = P1 (%.2f) noOfMasses = %d lA = %d \n",whichOfThree,5+noOfMasses*3+lA+whichOfThree*2,a,b,helJ,noOfMasses,lA);
  #endif

  //fptype a = p[indices[7]];
  //fptype b = p[indices[8]];



  result = exp(imUnit*b);

  result *= a;

  #ifdef MDEBUGGING
  // if(helJ==ZEROHEL) printf("H = (%.3f,%.3f)  a = %.3f  b = %.3f for helJ = 0 (%.2f) \n",result.real,result.imag,a,b,helJ);
  // if(helJ==M1HEL) printf("H = (%.3f,%.3f)  a = %.3f  b = %.3f for helJ = M1 (%.2f) \n",result.real,result.imag,a,b,helJ);
  // if(helJ==P1HEL) printf("H = (%.3f,%.3f)  a = %.3f  b = %.3f for helJ = P1(%.2f) \n",result.real,result.imag,a,b,helJ);
  #endif

  return result ;

}

EXEC_TARGET fptype Wignerd_R(fptype spinR, fptype helJ, fptype cKs)
{

  if (spinR==0.0)
    return 1. ;
  else if (spinR==1.0)
    if (helJ==M1HEL){
      fptype result = +(SIN(ACOS(cKs)) / root2);
      //printf("helj = m1 cKs = %.3f -> Wignerd_R = %.3f \n",cKs,result);
      return result;
      }
    else if (helJ==ZEROHEL){
      //printf("helj = 0 cKs = %.3f -> Wignerd_R = %.3f \n",cKs,cKs);
      return cKs ;}
    else if (helJ==P1HEL){
      fptype result = -(SIN(ACOS(cKs)) / root2);
      //printf("helj = p1 cKs = %.3f -> Wignerd_R = %.3f \n",cKs,result);
      return result;
    }
    else {
      printf("Wignerd_R Spin 1.0 not Allowed  returning 0\n");
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
      printf("Wignerd_R Spin 2.0 not Allowed  returning 0\n");
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
      printf("Wignerd_R Spin 3.0 not Allowed returning 0\n");
      //cout <<"helJ = " <<helJ <<" is not allowed for spinR-" <<spinR <<" Wigner d^{spinR}_{helJ,0} functions. Returning 0" <<endl;
      return 0 ;
    }
  else {
    printf("Wignerd_R Spin %.2f not Allowed  returning 0\n",spinR);
    //cout <<"spinR = " <<spinR <<" is not implemented for Wigner d^{spinR}_{helJ,0} functions at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl;
    return 0 ;
  }
}

EXEC_TARGET devcomplex<fptype> WignerD_J(fptype helJ, fptype helDmu, fptype angle,fptype cJ)
{
  devcomplex<fptype> imUnit(0.0,1.0);
  devcomplex<fptype> result(0.0,0.0);

  //printf("imUinit = (%.3f , %.3f) nullo = (%.3f , %.3f) \n",imUnit.real,imUnit.imag,nullo.real,nullo.imag);

  if (helJ==M1HEL) {
    if (helDmu==M1HEL){
      result = ((+1. + cJ)*exp(imUnit*angle))*.5;
      //printf("helj = m1 (%.3f) helDmu = (%.3f)  m1 cJ = %.3f angle(phi) = %.3f -> Wignerd_D = (%.3f,%.3f) \n",helJ,helDmu,cJ,angle,result.real,result.imag);
    }
    else if (helDmu==P1HEL){
      result = (-1.0)*((-1. + cJ)*exp(imUnit*angle))*.5;
      //printf("helj = m1 helDmu = p1 cJ = %.3f angle(phi) = %.3f -> Wignerd_D = (%.3f,%.3f) \n",cJ,angle,result.real,result.imag);
    }
    else {
      printf("WignerD_J :Not Allowed helDmu = %.2f with helJ = M1HEL returning 0\n",helDmu);
      //cout <<"helDmu = " <<helDmu <<" not allowed in \"WignerD_J\" functions for helJ = " <<helJ <<" at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
    }
  }
  else if (helJ==ZEROHEL) {
    //printf("Helj = zerohel \n");
    if (helDmu==M1HEL){
      result = devcomplex<fptype>((-1.0)*(SQRT(1. - POW(cJ,2))/root2),0.0);
      //printf("helj = 0 helDmu = m1 cJ = %.3f angle(phi) = %.3f -> Wignerd_D = (%.3f,%.3f) \n",cJ,angle,result.real,result.imag);
    }
    else if (helDmu==P1HEL){
      result = devcomplex<fptype>((SQRT(1. - POW(cJ,2))/root2),0.0);
      //printf("helj = 0 helDmu = p1 cJ = %.3f angle(phi) = %.3f -> Wignerd_D = (%.3f,%.3f) \n",cJ,angle,result.real,result.imag);
    }
    else {
      printf("WignerD_J :Not Allowed helDmu = %.2f with helJ = 0 returning 0\n",helDmu);
      //cout <<"helDmu = " <<helDmu <<" not allowed in \"WignerD_J\" functions for helJ = " <<helJ <<" at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
      }
  }
  else if(helJ==P1HEL) {
    if (helDmu==M1HEL)
      {
        result = (-1.0)*(-1. + cJ)/(2.*exp(imUnit*angle));
        //printf("helj = p1 helDmu = m1 cJ = %.3f angle(phi) = %.3f -> Wignerd_D = (%.3f,%.3f) \n",cJ,angle,result.real,result.imag);
      }
    else if (helDmu==P1HEL)
      {
        result = (+1. + cJ)/(2.*exp(imUnit*angle));
        //printf("helj = p1 helDmu = m1 cJ = %.3f angle(phi) = %.3f -> Wignerd_D = (%.3f,%.3f) \n",cJ,angle,result.real,result.imag);
      }
    else {
      printf("WignerD_J :Not Allowed helDmu = %.2f with helJ = P1 returning 0\n",helDmu);
      //cout <<"helDmu = " <<helDmu <<" not allowed in \"WignerD_J\" functions for helJ = " <<helJ <<" at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
      }
  } else {
    printf("WignerD_J :Not Allowed helDmu = %.2f with helJ = %.2f returning 0\n",helDmu,helJ);
    //cout <<"helJ = " <<helJ <<" not allowed in \"WignerD_J\" functions at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
  }

  return result;
}


EXEC_TARGET devcomplex<fptype> AngularTerm(fptype cJ, fptype cKs, fptype phi, fptype* p,unsigned int* indices,fptype spinR, fptype helJ, fptype helDmu,int iKStar)
{
  devcomplex<fptype> result;

  devcomplex<fptype> HH = H(p,indices,helJ,iKStar);

  fptype WR = Wignerd_R(spinR, helJ,cKs);
  devcomplex<fptype> cWR = conj( WignerD_J(helJ, helDmu, phi,cJ) );

  result = HH * WR * cWR ;

  #ifdef MDEBUGGING
  //printf("====== AngularTerm = (%.3f , %.3f) cJ = %.3f cKs = %.3f phi = %.3f \n ====== HH = (%.3f , %.3f) \n ====== WR = %.3f \n -- cKs = %.3f \n ====== cWR = (%.3f , %.3f) -- phi = %.3f  --- cJ = %.3f \n spin = %.3f \n",result.real,result.imag,cJ,cKs,phi,HH.real,HH.imag,WR,cKs,cWR.real,cWR.imag,phi,cJ,spinR);
  #endif

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

  //printf("\n RFunction (%.3f ,%.3f) for RMass = %.3f Mkp = %.3f  \n PmKP  = %.3f PRMass = %.3f    \n QmKP  = %.3f QRMass = %.3f  \n LminMom = %.3f LminR = %.3f    \n DB0 = %.3f DKs = %.3f    \n BlattWeisskopf LminM = %.3f    \n BlattWeisskopf LminR = %.3f    \n Power1  = %.3f    \n Power2  = %.3f    \n BW  = (%.3f ,%.3f)",RFunc.real,RFunc.imag,RMass,mkp,PmKP,PRMass,QmKP,QRMass,LminMom,LminR,DB0,DKs,blatt1,blatt2,pow1,pow2,bw.real,bw.imag);

  return RFunc ;
}

EXEC_TARGET devcomplex<fptype> matrixElement(fptype mkp, fptype cJ, fptype cKs, fptype phi, fptype* p,unsigned int* indices,fptype helDmu)
{
  unsigned int nOfKstar = indices[1];
  //int numberOfKStar = indices[0]/6;


  fptype psi_nS = p[indices[2]];
  fptype dRadB0 = p[indices[3]];
  fptype dRadKs = p[indices[4]];
  //int nKStars = (int) p[indices[4]];

  devcomplex<fptype> matrixElement (0.0,0.0);

  // K+ and pi- have 0 spin -> second last argument of K* RFunction is = spin(K*)
  //int numParams = indices[0];
  //printf("psi_nS = %f dRadB0 = %f dRadKs = %f nKStars = %d numparm = %d \n",psi_nS,dRadB0,dRadKs,nOfKstar,numParams);

  int lastAmplitude = 0;

  for (int iKStar=0; iKStar<nOfKstar; ++iKStar) {

    fptype Mass = p[indices[4+1+iKStar*3]];
    fptype Gamma = p[indices[4+2+iKStar*3]];
    fptype Spin = p[indices[4+3+iKStar*3]];

    // 2) psi_ns 3) DB0 4) DBK 5) M1 6) G1 7) S1 8) A1 9)B1 10)M2
    devcomplex<fptype> matrixElement_R(0.0,0.0);

    #ifdef MDEBUGGING
    printf("Mass = %f Gamma = %f Spin = %f psi_nS = %f dRadB0 = %f dRadKs = %f \n",Mass,Gamma,Spin,psi_nS,dRadB0,dRadKs);
    #endif

    if (Spin==0.0) { // for spin0 K*, third last argument = spin(psi_nS) = spin.Atoi() + 1 = 1
      matrixElement_R = RFunction(mkp,Mass,Gamma, MBd, Spin+1, Spin, dRadB0, dRadKs,psi_nS) *
      AngularTerm(cJ,cKs,phi,p,indices,Spin, ZEROHEL, helDmu,lastAmplitude) ;
      ++lastAmplitude;
    }
       else
              {
                matrixElement_R = RFunction(mkp,Mass,Gamma, MBd, Spin-1, Spin,dRadB0,dRadKs,psi_nS) *
	               ( AngularTerm(cJ,cKs,phi,p,indices,Spin, M1HEL, helDmu,lastAmplitude) + AngularTerm(cJ,cKs,phi,p,indices, Spin, ZEROHEL, helDmu,lastAmplitude) + AngularTerm(cJ,cKs,phi,p,indices,Spin, P1HEL, helDmu,lastAmplitude) ) ;
                   lastAmplitude +=3;
               }
    //cout <<"\nAngularTerm.Rho() for " <<R <<" = " <<(AngularTerm(R, spin, "0", helDmu)).Rho() <<endl;
    //cout <<"matrixElement for (R,helDmu) = (" <<R <<"," <<helDmu <<") = H(R,helJ) * RFunction * AngularTerm = " <<matrixElement_R <<endl;
    matrixElement += matrixElement_R;
    //cout <<"matrixElement_R.Rho2() for (R,helDmu) = (" <<R <<"," <<helDmu <<") = " <<matrixElement_R.Rho2() <<"\n\n" <<endl;
  }

  //printf("======= Matrix Element HEL = %.2f \n  Mass KPi = %.3f cJ = %.3f  cKs = %.3f phi = %.3f \n mE = ( %.3f , %.3f )",helDmu,mkp,cJ,cKs,phi,matrixElement.real,matrixElement.imag);

  return matrixElement;

}

EXEC_TARGET fptype ME2(fptype mkp, fptype cJ, fptype cKs, fptype phi, fptype* p,unsigned int* indices)
{
  //cout <<"\nME(\"m1\") + ME(\"p1\") = " <<ME("m1") <<" + " <<ME("p1") <<endl;
  //cout <<"ME(\"m1\").Rho2() + ME(\"p1\").Rho2() = " <<ME("m1").Rho2() <<" + " <<ME("p1").Rho2() <<endl;

  fptype finalDevice1 = matrixElement(mkp,cJ,cKs,phi,p,indices,M1HEL).abs2();
  fptype finalDevice2 = matrixElement(mkp,cJ,cKs,phi,p,indices,P1HEL).abs2();

  //printf("Matrix Element HelDMu m1 = %.3f Matrix Element HelDMu p1 = %.3f ( mkp = %.3f - cJ = %.3f cKs = %.3f phi = %.3f)\n",finalDevice1,finalDevice2,mkp,cJ,cKs,phi);

  return finalDevice1 + finalDevice2 ;
}

EXEC_TARGET fptype PhiPHSP(fptype mkp,fptype psiN)
{
    //printf("=======Phase Space mpk = %.2f psiN = %.2f \n",mkp,psiN);
    const fptype psin = psiN;
    fptype p = Pmom(mkp,psin);
    fptype q = Qmom(mkp);
    //fptype phsp = p*q;
    //printf(" Mass KPi = %.3f Phase space = %.3f\n",mkp,phsp);
    //printf("==================");

    return Pmom(mkp,psiN) * Qmom(mkp);
}



EXEC_TARGET fptype device_Matrix_B0 (fptype* evt, fptype* p, unsigned int* indices) {

  #ifdef MDEBUGGING
  printf("Zero paramater set %d %d %d %d %.2f %.2f %.2f %.2f  %.2f \n",indices[0],indices[1],indices[2],indices[3],p[indices[0]],p[indices[1]],p[indices[2]],p[indices[3]],p[indices[4]]);
  printf("First K paramater set  %.2f %.2f  %.2f  %.2f  %.2f \n",p[indices[5]],p[indices[6]],p[indices[7]],p[indices[8]],p[indices[9]]);
  printf("Second K paramater set  %.2f %.2f  %.2f  %.2f  %.2f \n",p[indices[10]],p[indices[11]],p[indices[12]],p[indices[13]],p[indices[14]]);
  printf("Third K paramater set  %.2f %.2f  %.2f  %.2f  %.2f \n",p[indices[15]],p[indices[16]],p[indices[17]],p[indices[18]],p[indices[19]]);
  #endif

  fptype mkp = evt[indices[2 + indices[0]]];
  fptype mPsiP = evt[indices[2 + indices[0]]+1];
  fptype cJ = evt[indices[2 + indices[0]]+2];
  //fptype cKs = evt[indices[2 + indices[0]]+2];
  fptype phi = evt[indices[2 + indices[0]]+3];
  fptype b0Flag = evt[indices[2 + indices[0]]+4];

  if (b0Flag<0.0)
    phi *= -1.0;

  fptype psi_nS = p[indices[2]];

  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
  fptype MPsi_nS = 0.;
  if (psi_nS==1.0)
    MPsi_nS = MJpsi;
  else if (psi_nS==2.0)
    MPsi_nS = MPsi2S;
  else {
    printf("\nMatrix P.d.f not configured for psi_nS = %.0f",psi_nS);
    //printf("mpk = %.2f (%.2f - %.2f) cJ = %.2f cKs = %.2f phi = %.2f \n",mkp,MBd - MPsi_nS,MKaon + MPion,cJ,cKs,phi);
  }

  fptype mKP2 = mkp*mkp;
  fptype mPsiP2 = mPsiP*mPsiP;
  fptype MPsi_nS2 = MPsi_nS*MPsi_nS;

  if ((mkp < MKaon + MPion) || (mkp > MBd - MPsi_nS) || (mPsiP < MPsi_nS + MPion) || (mPsiP > MBd - MKaon)) {
    //printf("Returning 0: point out of the Dalitz borders!\n");
    return 0.; }

  fptype cKs = cosTheta_FromMasses(mKP2, mPsiP2, MPsi_nS2, MBd2, MKaon2, MPion2);

  //fptype dRadB0 = p[indices[3]];
  //fptype dRadKs = p[indices[4]];
  //printf("Hei mpk = %.2f cJ = %.2f cKs = %.2f phi = %.2f psi_nS = %.2f mPSi = %.2f \n",mkp,cJ,cKs,phi,psi_nS,mPsiP);

  if (FABS(cKs) > 1.0 || FABS(phi) > devPi || FABS(cJ) > 1.0) {
    //printf("\nReturning 0 : ckS > 1 or < -1 ");
    return 0.; }
  else {
    fptype MEME = ME2(mkp,cJ,cKs,phi,p,indices);
    fptype phiPhase = PhiPHSP(mkp,psi_nS);
    fptype result = MEME * phiPhase;
    //printf("Device Matrix = %.3f MEME = %.3f phiPhase = %.3f with mkp = %.2f cJ = %.2f cKs = %.2f phi = %.2f mPSIP = %.2f \n",result,MEME,phiPhase,mkp,cJ,cKs,phi,mPsiP);
    return result;
  }
}

EXEC_TARGET fptype device_Matrix(fptype* evt, fptype* p, unsigned int* indices) {

  #ifdef MDEBUGGING
  printf("Zero paramater set %d %d %d %d %.2f %.2f %.2f %.2f  %.2f \n",indices[0],indices[1],indices[2],indices[3],p[indices[0]],p[indices[1]],p[indices[2]],p[indices[3]],p[indices[4]]);
  printf("First K paramater set  %.2f %.2f  %.2f  %.2f  %.2f \n",p[indices[5]],p[indices[6]],p[indices[7]],p[indices[8]],p[indices[9]]);
  printf("Second K paramater set  %.2f %.2f  %.2f  %.2f  %.2f \n",p[indices[10]],p[indices[11]],p[indices[12]],p[indices[13]],p[indices[14]]);
  printf("Third K paramater set  %.2f %.2f  %.2f  %.2f  %.2f \n",p[indices[15]],p[indices[16]],p[indices[17]],p[indices[18]],p[indices[19]]);
  #endif

  fptype mkp = evt[indices[2 + indices[0]]];
  fptype mPsiP = evt[indices[2 + indices[0]]+1];
  fptype cJ = evt[indices[2 + indices[0]]+2];
  //fptype cKs = evt[indices[2 + indices[0]]+2];
  fptype phi = evt[indices[2 + indices[0]]+3];

  fptype psi_nS = p[indices[2]];

  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
  fptype MPsi_nS = 0.;
  if (psi_nS==1.0)
    MPsi_nS = MJpsi;
  else if (psi_nS==2.0)
    MPsi_nS = MPsi2S;
  else {
    printf("\nMatrix P.d.f not configured for psi_nS = %.0f",psi_nS);
    //printf("mpk = %.2f (%.2f - %.2f) cJ = %.2f cKs = %.2f phi = %.2f \n",mkp,MBd - MPsi_nS,MKaon + MPion,cJ,cKs,phi);
  }

  fptype mKP2 = mkp*mkp;
  fptype mPsiP2 = mPsiP*mPsiP;
  fptype MPsi_nS2 = MPsi_nS*MPsi_nS;

  if ((mkp < MKaon + MPion) || (mkp > MBd - MPsi_nS) || (mPsiP < MPsi_nS + MPion) || (mPsiP > MBd - MKaon)) {
    //printf("Returning 0: point out of the Dalitz borders!\n");
    return 0.; }

  fptype cKs = cosTheta_FromMasses(mKP2, mPsiP2, MPsi_nS2, MBd2, MKaon2, MPion2);

  //fptype dRadB0 = p[indices[3]];
  //fptype dRadKs = p[indices[4]];
  //printf("Hei mpk = %.2f cJ = %.2f cKs = %.2f phi = %.2f psi_nS = %.2f mPSi = %.2f \n",mkp,cJ,cKs,phi,psi_nS,mPsiP);

  if (FABS(cKs) > 1.0 || FABS(phi) > devPi || FABS(cJ) > 1.0) {
    //printf("\nReturning 0 : ckS > 1 or < -1 ");
    return 0.; }
  else {
    fptype MEME = ME2(mkp,cJ,cKs,phi,p,indices);
    fptype phiPhase = PhiPHSP(mkp,psi_nS);
    fptype result = MEME * phiPhase;
    //printf("Device Matrix = %.3f MEME = %.3f phiPhase = %.3f with mkp = %.2f cJ = %.2f cKs = %.2f phi = %.2f mPSIP = %.2f \n",result,MEME,phiPhase,mkp,cJ,cKs,phi,mPsiP);
    return result;
  }

}

EXEC_TARGET fptype device_Matrix_Bar(fptype* evt, fptype* p, unsigned int* indices) {

  #ifdef MDEBUGGING
  printf("Zero paramater set %d %d %d %d %.2f %.2f %.2f %.2f  %.2f \n",indices[0],indices[1],indices[2],indices[3],p[indices[0]],p[indices[1]],p[indices[2]],p[indices[3]],p[indices[4]]);
  printf("First K paramater set  %.2f %.2f  %.2f  %.2f  %.2f \n",p[indices[5]],p[indices[6]],p[indices[7]],p[indices[8]],p[indices[9]]);
  printf("Second K paramater set  %.2f %.2f  %.2f  %.2f  %.2f \n",p[indices[10]],p[indices[11]],p[indices[12]],p[indices[13]],p[indices[14]]);
  printf("Third K paramater set  %.2f %.2f  %.2f  %.2f  %.2f \n",p[indices[15]],p[indices[16]],p[indices[17]],p[indices[18]],p[indices[19]]);
  #endif

  fptype mkp = evt[indices[2 + indices[0]]];
  fptype mPsiP = evt[indices[2 + indices[0]]+1];
  fptype cJ = evt[indices[2 + indices[0]]+2];
  //fptype cKs = evt[indices[2 + indices[0]]+2];
  fptype phi = -evt[indices[2 + indices[0]]+3];

  fptype psi_nS = p[indices[2]];

  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
  fptype MPsi_nS = 0.;
  if (psi_nS==1.0)
    MPsi_nS = MJpsi;
  else if (psi_nS==2.0)
    MPsi_nS = MPsi2S;
  else {
    printf("\nMatrix P.d.f not configured for psi_nS = %.0f",psi_nS);
    //printf("mpk = %.2f (%.2f - %.2f) cJ = %.2f cKs = %.2f phi = %.2f \n",mkp,MBd - MPsi_nS,MKaon + MPion,cJ,cKs,phi);
  }

  fptype mKP2 = mkp*mkp;
  fptype mPsiP2 = mPsiP*mPsiP;
  fptype MPsi_nS2 = MPsi_nS*MPsi_nS;

  if ((mkp < MKaon + MPion) || (mkp > MBd - MPsi_nS) || (mPsiP < MPsi_nS + MPion) || (mPsiP > MBd - MKaon)) {
    //printf("Returning 0: point out of the Dalitz borders!\n");
    return 0.; }

  fptype cKs = cosTheta_FromMasses(mKP2, mPsiP2, MPsi_nS2, MBd2, MKaon2, MPion2);

  //fptype dRadB0 = p[indices[3]];
  //fptype dRadKs = p[indices[4]];
  //printf("Hei mpk = %.2f cJ = %.2f cKs = %.2f phi = %.2f psi_nS = %.2f mPSi = %.2f \n",mkp,cJ,cKs,phi,psi_nS,mPsiP);

  if (FABS(cKs) > 1.0 || FABS(phi) > devPi || FABS(cJ) > 1.0) {

    //printf("\nReturning 0 : ckS > 1 or < -1 ");
    return 0.; }
  else {
    fptype MEME = ME2(mkp,cJ,cKs,phi,p,indices);
    fptype phiPhase = PhiPHSP(mkp,psi_nS);
    fptype result = MEME * phiPhase;
    //printf("Device Matrix = %.3f MEME = %.3f phiPhase = %.3f with mkp = %.2f cJ = %.2f cKs = %.2f phi = %.2f mPSIP = %.2f \n",result,MEME,phiPhase,mkp,cJ,cKs,phi,mPsiP);
    return result;
  }

}


MEM_DEVICE device_function_ptr ptr_to_Matrix_Bar = device_Matrix_Bar;
MEM_DEVICE device_function_ptr ptr_to_Matrix_B0 = device_Matrix_B0;
MEM_DEVICE device_function_ptr ptr_to_Matrix = device_Matrix;

//MEM_DEVICE device_function_ptr ptr_to_Matrix_Point = device_Matrix_Point;
//MEM_DEVICE device_function_ptr ptr_to_Matrix_Bin = device_Matrix_Bin;

__host__ MatrixPdf::MatrixPdf(std::string n, Variable* _mkp, Variable* _mJP,Variable* _cJ, Variable* _phi,
        std::vector<Variable*> _Masses,std::vector<Variable*> _Gammas,std::vector<Variable*> _Spins,std::vector<Variable*> _a,std::vector<Variable*> _b,
        Variable* _psi_nS, Variable* _dRadB0, Variable* _dRadKs)
  : GooPdf(0, n),
  psi_nS(_psi_nS),dRadB0(_dRadB0),dRadKs(_dRadKs)
{

  unsigned int noOfKStars = 0;
  unsigned int noOfMasses = (int) _Masses.size();

  for (int j = 0 ; j < _Masses.size(); j++) {

    if(_Spins[j]->value > 0.0)
      noOfKStars += 3;
    else
      ++noOfKStars;
  }

  printf("Number of K* \t\t\t = %d\n", noOfKStars);
  printf("Number of masses \t\t = %d\n", noOfMasses);
  printf("Amplitudes vector size \t\t = %d \n",_a.size());

  if(noOfKStars != (int) _a.size())
      abortWithCudaPrintFlush(__FILE__, __LINE__, "No. of kStars different from no. of amplitudes and phases provided \n");

  registerObservable(_mkp);
  registerObservable(_mJP);
  registerObservable(_cJ);
  registerObservable(_phi);

  std::vector<unsigned int> pindices;

  pindices.push_back(noOfMasses);
  pindices.push_back(registerParameter(_psi_nS));  // p[indices[2]]
  pindices.push_back(registerParameter(_dRadB0));  // p[indices[3]]
  pindices.push_back(registerParameter(_dRadKs));  // p[indices[4]]

  //Parameter vector
  // psi_ns dRadB0 dRadKs m1 g1 s1 m2 g2 s2 . . . mn  gn  sn   a1          b1             a2 b2 . . . an                  bn                    //
  //Indeces
  // 2      3      4      5  6  7  8  9  10        4+n 5+n 6+n  4+nKstars*3 4+nKstars*3+1              4+nKstars*3+(n-1)*2 4+nKstars*3+(n-1)*2+1
  for (int j = 0 ; j < _Masses.size(); j++) {

    pindices.push_back(registerParameter(_Masses[j]));  // p[indices[5]]
    pindices.push_back(registerParameter(_Gammas[j]));  // p[indices[6]]
    pindices.push_back(registerParameter(_Spins[j]));
  }

  for (int j = 0 ; j < (int)_a.size(); j++) {
    //pindices.push_back(registerParameter(_helj[j]));
    pindices.push_back(registerParameter(_a[j]));
    pindices.push_back(registerParameter(_b[j]));

  }

  GET_FUNCTION_ADDR(ptr_to_Matrix);
  //GET_INTEGRAL_ADDR(ptr_to_Matrix_Bin);
  //GET_ATPOINTS_ADDR(ptr_to_Matrix_Point);
  initialise(pindices);
}

__host__ MatrixPdf::MatrixPdf(std::string n, std::vector<Variable*> _Masses,std::vector<Variable*> _Gammas,std::vector<Variable*> _Spins,std::vector<Variable*> _a,std::vector<Variable*> _b,
                  Variable* _psi_nS, Variable* _dRadB0, Variable* _dRadKs,Variable* _mkp, Variable* _mJP,Variable* _cJ, Variable* _phi)
  : GooPdf(0, n),
  psi_nS(_psi_nS),dRadB0(_dRadB0),dRadKs(_dRadKs)
{

  unsigned int noOfKStars = 0;
  unsigned int noOfMasses = (int) _Masses.size();

  for (int j = 0 ; j < _Masses.size(); j++) {

    if(_Spins[j]->value > 0.0)
      noOfKStars += 3;
    else
      ++noOfKStars;
  }

  printf("Number of K* \t\t\t = %d\n", noOfKStars);
  printf("Number of masses \t\t = %d\n", noOfMasses);
  printf("Amplitudes vector size \t\t = %d \n",_a.size());

  if(noOfKStars != (int) _a.size())
      abortWithCudaPrintFlush(__FILE__, __LINE__, "No. of kStars different from no. of amplitudes and phases provided \n");

  registerObservable(_mkp);
  registerObservable(_mJP);
  registerObservable(_cJ);
  registerObservable(_phi);

  std::vector<unsigned int> pindices;

  pindices.push_back(noOfMasses);
  pindices.push_back(registerParameter(_psi_nS));  // p[indices[2]]
  pindices.push_back(registerParameter(_dRadB0));  // p[indices[3]]
  pindices.push_back(registerParameter(_dRadKs));  // p[indices[4]]

  //Parameter vector
  // psi_ns dRadB0 dRadKs m1 g1 s1 m2 g2 s2 . . . mn  gn  sn   a1          b1             a2 b2 . . . an                  bn                    //
  //Indeces
  // 2      3      4      5  6  7  8  9  10        4+n 5+n 6+n  4+nKstars*3 4+nKstars*3+1              4+nKstars*3+(n-1)*2 4+nKstars*3+(n-1)*2+1
  for (int j = 0 ; j < _Masses.size(); j++) {

    pindices.push_back(registerParameter(_Masses[j]));  // p[indices[5]]
    pindices.push_back(registerParameter(_Gammas[j]));  // p[indices[6]]
    pindices.push_back(registerParameter(_Spins[j]));
  }

  for (int j = 0 ; j < (int)_a.size(); j++) {
    //pindices.push_back(registerParameter(_helj[j]));
    pindices.push_back(registerParameter(_a[j]));
    pindices.push_back(registerParameter(_b[j]));

  }

  GET_FUNCTION_ADDR(ptr_to_Matrix_Bar);
  //GET_INTEGRAL_ADDR(ptr_to_Matrix_Bin);
  //GET_ATPOINTS_ADDR(ptr_to_Matrix_Point);
  initialise(pindices);
}

__host__ MatrixPdf::MatrixPdf(std::string n, Variable* _mkp, Variable* _mJP,Variable* _cJ, Variable* _phi,Variable* _B0beauty,
        std::vector<Variable*> _Masses,std::vector<Variable*> _Gammas,std::vector<Variable*> _Spins,std::vector<Variable*> _a,std::vector<Variable*> _b,
        Variable* _psi_nS, Variable* _dRadB0, Variable* _dRadKs)
  : GooPdf(0, n),
  psi_nS(_psi_nS),dRadB0(_dRadB0),dRadKs(_dRadKs)
{

  unsigned int noOfKStars = 0;
  unsigned int noOfMasses = (int) _Masses.size();

  for (int j = 0 ; j < _Masses.size(); j++) {

    if(_Spins[j]->value>0.0)
      noOfKStars += 3;
    else
      ++noOfKStars;
  }

  printf("Number of K* \t\t\t = %d\n", noOfKStars);
  printf("Number of masses \t\t = %d\n", noOfMasses);
  printf("Amplitudes vector size \t\t = %d \n",_a.size());

  if(noOfKStars != (int) _a.size()) abortWithCudaPrintFlush(__FILE__, __LINE__, "No. of kStars different from no. of amplitudes and phases provided \n");

  registerObservable(_mkp);
  registerObservable(_mJP);
  registerObservable(_cJ);
  registerObservable(_phi);
  registerObservable(_B0beauty);

  std::vector<unsigned int> pindices;

  pindices.push_back(noOfMasses);
  pindices.push_back(registerParameter(_psi_nS));  // p[indices[2]]
  pindices.push_back(registerParameter(_dRadB0));  // p[indices[3]]
  pindices.push_back(registerParameter(_dRadKs));  // p[indices[4]]

  //Parameter vector
  // psi_ns dRadB0 dRadKs m1 g1 s1 m2 g2 s2 . . . mn  gn  sn   a1          b1             a2 b2 . . . an                  bn                    //
  //Indeces
  // 2      3      4      5  6  7  8  9  10        4+n 5+n 6+n  4+nKstars*3 4+nKstars*3+1              4+nKstars*3+(n-1)*2 4+nKstars*3+(n-1)*2+1
  for (int j = 0 ; j < _Masses.size(); j++) {

    pindices.push_back(registerParameter(_Masses[j]));  // p[indices[5]]
    pindices.push_back(registerParameter(_Gammas[j]));  // p[indices[6]]
    pindices.push_back(registerParameter(_Spins[j]));
  }

  for (int j = 0 ; j < (int)_a.size(); j++) {
    //pindices.push_back(registerParameter(_helj[j]));
    pindices.push_back(registerParameter(_a[j]));
    pindices.push_back(registerParameter(_b[j]));

  }



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
  GET_FUNCTION_ADDR(ptr_to_Matrix_B0);
  //GET_INTEGRAL_ADDR(ptr_to_Matrix_Bin);
  //GET_ATPOINTS_ADDR(ptr_to_Matrix_Point);
  initialise(pindices);
}
