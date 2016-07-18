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

EXEC_TARGET devcomplex<fptype> matrixElement(string helDmu) const
{
  /*
  // K+ and pi- have 0 spin -> second last argument of K* RFunction is = spin(K*)
  return
    RFunction(M892, G892, MBd, 0, 1, dRad) * ( AngularTerm("K*(892)", "1", "m1", helDmu) +
						      AngularTerm("K*(892)", "1", "0", helDmu) +
						      AngularTerm("K*(892)", "1", "p1", helDmu) )
    // + ...
    ;
  // any other K* should be added above
  */
  devcomplex<fptype> matrixElement (0.0,0.0);
  // K+ and pi- have 0 spin -> second last argument of K* RFunction is = spin(K*)
  for (int iKstar_S=0; iKstar_S<(int)Kstar_spin.size(); ++iKstar_S) {
    TString R = Kstar_spin[iKstar_S].first ;
    TString spin = R(Kstar_spin[iKstar_S].first.Length() -1) ;
    TString mass = R(0, Kstar_spin[iKstar_S].first.Length() -2) ;
    TComplex matrixElement_R = 0.;
    if (spin.EqualTo("0")) { // for spin0 K*, third last argument = spin(psi_nS) = spin.Atoi() + 1 = 1
      matrixElement_R = RFunction(Kstar_spin[iKstar_S].second.first, Kstar_spin[iKstar_S].second.second, MBd, spin.Atoi()+1, spin.Atoi(), dRadB0, dRadKs) *
	               AngularTerm(R, spin, "0", helDmu) ;
    } else { // for non-0 spin K*, third last argument = spin(K*) - spin(psi_nS) = spin.Atoi() - 1
      matrixElement_R = RFunction(Kstar_spin[iKstar_S].second.first, Kstar_spin[iKstar_S].second.second, MBd, spin.Atoi()-1, spin.Atoi(), dRadB0, dRadKs) *
	               ( AngularTerm(R, spin, "m1", helDmu) + AngularTerm(R, spin, "0", helDmu) + AngularTerm(R, spin, "p1", helDmu) ) ;
    }
    //cout <<"\nAngularTerm.Rho() for " <<R <<" = " <<(AngularTerm(R, spin, "0", helDmu)).Rho() <<endl;
    //cout <<"matrixElement for (R,helDmu) = (" <<R <<"," <<helDmu <<") = H(R,helJ) * RFunction * AngularTerm = " <<matrixElement_R <<endl;
    matrixElement += matrixElement_R;
    //cout <<"matrixElement_R.Rho2() for (R,helDmu) = (" <<R <<"," <<helDmu <<") = " <<matrixElement_R.Rho2() <<"\n\n" <<endl;
  }
  return matrixElement ;

}

__host__ MatrixPdf::MatrixPdf (std::string n, Variable* _x,
  Variable* _cJ, Variable* _cKs, Variable* _phi,
  const std::vector< std::vector<fptype>> _KstarDotSpin,
  const vector<std::string> _varNames,
  const std::vector<Variable*> _amplitudeGooVars,
  const std::string _psi_nS,
  const fptype& _dRadB0, const fptype& _dRadKs)
  : GooPdf(_x, n),KstarDotSpin(_KstarDotSpin),
  varNames(_varNames),amplitudeGooVars(_amplitudeGooVars),
  psi_nS(_psi_nS),dRadB0(_dRadB0),dRadKs(_dRadKs)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(_cJ));
  pindices.push_back(registerParameter(_cKs));
  pindices.push_back(registerParameter(_phi));

  for (vector<Variable*>::iterator  = _amplitudeGooVars.begin(); v != _amplitudeGooVars.end(); ++v) {
    pindices.push_back(registerParameter(*v));
  }

  GET_FUNCTION_ADDR(ptr_to_Matrix);
  GET_INTEGRAL_ADDR(ptr_to_Matrix_Bin);
  GET_ATPOINTS_ADDR(ptr_to_Matrix_Point);
  initialiseInt(pindices);
}

MEM_DEVICE device_function_ptr ptr_to_Matrix = device_Matrix;
MEM_DEVICE device_function_ptr ptr_to_Matrix_Point = device_Matrix_Point;
MEM_DEVICE device_function_ptr ptr_to_Matrix_Bin = device_Matrix_Bin;


/*
EXEC_TARGET devcomplex<fptype> WignerD_J(int helJ,int helDmu, fptype angle) const
{

  devcomplex<fptype> imJ(0.,1.);

  if (helJ==M1) {
    if (helDmu==M1)
      return +((+1. + cJ)*exp(imJ*angle))/2.;
    else if (helDmu==P1)
      return -((-1. + cJ)*exp(imJ*angle))/2.;
    else {
      cout <<"helDmu = " <<helDmu <<" not allowed in \"WignerD_J\" functions for helJ = " <<helJ <<" at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
      return 0; }
  } else if (helJ==ZERO) {
    if (helDmu==M1)
      return -(pow(1. - pow(cJ,2))/TMath::Sqrt2());
    else if (helDmu==P1)
      return +(pow(1. - pow(cJ,2))/TMath::Sqrt2());
    else {
      cout <<"helDmu = " <<helDmu <<" not allowed in \"WignerD_J\" functions for helJ = " <<helJ <<" at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
      return 0; }
  } else if(helJ==P1) {
    if (helDmu==M1)
      return -(-1. + cJ)/(2.*exp(imJ*angle));
    else if (helDmu==P1)
      return +(+1. + cJ)/(2.*exp(imJ*angle));
    else {
      cout <<"helDmu = " <<helDmu <<" not allowed in \"WignerD_J\" functions for helJ = " <<helJ <<" at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
      return 0; }
  } else {
    cout <<"helJ = " <<helJ <<" not allowed in \"WignerD_J\" functions at the moment. Returning 0 -> \"AngularTerm\" = 0" <<endl ;
    return 0;
  }

}

// H term in slide 11 second last line for Lambda*(1600)
EXEC_TARGET fptype HLs1600(std::string help) const
{

    if(help==M1H)
        return -1.;
    else if(help==P1H)
        return 1.;
    else { cout <<"WARNING! In \"HLs1600\" function: help = " <<help <<" -> returning 0" <<endl;
        return 0.;
    }
}

// H term in slide 11 second last line for Lambda*(1670)
EXEC_TARGET fptype HLs1670(std::string help) const
{

    if(help==M1H)
        return 1.;
    else if(help==P1H)
        return 1.;
    else { cout <<"WARNING! In \"HLs1670\" function: help = " <<help <<" -> returning 0" <<endl;
        return 0.;
    }

}

EXEC_TARGET devcomplex<fptype> ME( std::string helDmu ) const
{
  /*
  // K+ and pi- have 0 spin -> second last argument of K* RFunction is = spin(K*)
  return
    RFunction(M892, G892, MBd, 0, 1, dRad) * ( AngularTerm("K*(892)", "1", M1, helDmu) +
						      AngularTerm("K*(892)", "1", ZERO, helDmu) +
						      AngularTerm("K*(892)", "1", "p1", helDmu) )
    // + ...
    ;
  // any other K* should be added above


  devcomplex<fptype> matrixElement(0.,0.);
  // K+ and pi- have 0 spin -> second last argument of K* RFunction is = spin(K*)
  for (Int_t iKstar_S=0; iKstar_S<(Int_t)Kstar_spin.size(); ++iKstar_S) {
    TString R = Kstar_spin[iKstar_S].first ;
    TString spin = R(Kstar_spin[iKstar_S].first.Length() -1) ;
    TString mass = R(0, Kstar_spin[iKstar_S].first.Length() -2) ;
    devcomplex<fptype> matrixElement_R = 0.;
    if (spin.EqualTo(ZERO)) { // for spin0 K*, third last argument = spin(psi_nS) = spin.Atoi() + 1 = 1
      matrixElement_R = RFunction(Kstar_spin[iKstar_S].second.first, Kstar_spin[iKstar_S].second.second, MBd, spin.Atoi()+1, spin.Atoi(), dRadB0, dRadKs) *
	               AngularTerm(R, spin, ZERO, helDmu) ;
    } else { // for non-0 spin K*, third last argument = spin(K*) - spin(psi_nS) = spin.Atoi() - 1
      matrixElement_R = RFunction(Kstar_spin[iKstar_S].second.first, Kstar_spin[iKstar_S].second.second, MBd, spin.Atoi()-1, spin.Atoi(), dRadB0, dRadKs) *
	               ( AngularTerm(R, spin, M1, helDmu) + AngularTerm(R, spin, ZERO, helDmu) + AngularTerm(R, spin, "p1", helDmu) ) ;
    }
    //cout <<"\nAngularTerm.Rho() for " <<R <<" = " <<(AngularTerm(R, spin, ZERO, helDmu)).Rho() <<endl;
    //cout <<"matrixElement for (R,helDmu) = (" <<R <<"," <<helDmu <<") = H(R,helJ) * RFunction * AngularTerm = " <<matrixElement_R <<endl;
    matrixElement += matrixElement_R;
    //cout <<"matrixElement_R.Rho2() for (R,helDmu) = (" <<R <<"," <<helDmu <<") = " <<matrixElement_R.Rho2() <<"\n\n" <<endl;
  }
  return matrixElement ;

}

EXEC_TARGET fptype ME2() const
{
  //cout <<"\nME(\"m1\") + ME(\"p1\") = " <<ME(M1) <<" + " <<ME("p1") <<endl;
  //cout <<"ME(\"m1\").Rho2() + ME(\"p1\").Rho2() = " <<ME(M1).Rho2() <<" + " <<ME("p1").Rho2() <<endl;
  return ME(M1).Rho2() + ME(P1).Rho2() ;
}

//TComplex myPDF::PDF() const
EXEC_TARGET fptype PDF() const
{
  //cout <<"\nME2() = " <<ME2() <<endl;
  return ME2() * PhiPHSP(mKP); // missing * efficiency(from reconstructed PHSP MC)

}

*/


EXEC_TARGET fptype device_Matrix (fptype* point, fptype* p, unsigned int* indices) {

  fptype x = evt[indices[2 + indices[0]]];

  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE
   fptype MPsi_nS = 0.;
   if (psi_nS==1)
     MPsi_nS = 3.096916;
   else if (psi_nS==2)
     MPsi_nS = 3.686109;
   else
     cout <<"psi_nS = " <<psi_nS <<" not allowed in the \"evaluate\" function at the moment. Keeping MPsi_nS to 0" <<endl;

  if ((mKP < MKaon + MPion) || (mKP > MBd - MPsi_nS))
    return 0.;
  else
    return PDF();
  return 0;

}

Double_t myPDF::ME2() const
{
  //cout <<"\nME(\"m1\") + ME(\"p1\") = " <<ME("m1") <<" + " <<ME("p1") <<endl;
  //cout <<"ME(\"m1\").Rho2() + ME(\"p1\").Rho2() = " <<ME("m1").Rho2() <<" + " <<ME("p1").Rho2() <<endl;
  return ME("m1").Rho2() + ME("p1").Rho2() ;
}

//TComplex myPDF::PDF() const
Double_t myPDF::PDF() const
{
  //cout <<"\nME2() = " <<ME2() <<endl;
  return ME2() * PhiPHSP(mKP); // missing * efficiency(from reconstructed PHSP MC)

}



EXEC_TARGET fptype device_Matrix_Point (fptype* point, fptype* p, unsigned int* indices) {

  return 0;

}

EXEC_TARGET fptype device_Matrix_Bin (fptype* point, fptype* p, unsigned int* indices) {

  return 0;

}
