#include "VoigtianPdf.hh"
#include <limits>
#include "Faddeeva.hh"
#include "devcomplex.hh"

#define M_2PI 6.28318530717958
//#define ROOT2 1.41421356

// tables for Pade approximation
MEM_CONSTANT fptype C[7] = { 65536.0, -2885792.0, 69973904.0, -791494704.0,
			     8962513560.0, -32794651890.0, 175685635125.0 };
MEM_CONSTANT fptype D[7] = { 192192.0, 8648640.0, 183783600.0, 2329725600.0,
			     18332414100.0, 84329104860.0, 175685635125.0 };


//#define UNROLL_LOOP 1

#ifndef UNROLL_LOOP
MEM_CONSTANT fptype n1[12] = { 0.25, 1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0, 20.25, 25.0, 30.25, 36.0 };
MEM_CONSTANT fptype e1[12] = { 0.7788007830714049,    0.3678794411714423,
			       1.053992245618643e-1,  1.831563888873418e-2,
			       1.930454136227709e-3,  1.234098040866795e-4,
			       4.785117392129009e-6,  1.125351747192591e-7,
			       1.605228055185612e-9,  1.388794386496402e-11,
			       7.287724095819692e-14, 2.319522830243569e-16 };

// table 2: coefficients for h = 0.53
MEM_CONSTANT fptype n2[12] = { 0.2809, 1.1236, 2.5281, 4.4944, 7.0225, 10.1124,
			 13.7641, 17.9776, 22.7529, 28.09, 33.9889, 40.4496 };
MEM_CONSTANT fptype e2[12] = { 0.7551038420890235,    0.3251072991205958,
			       7.981051630007964e-2,  1.117138143353082e-2,
			       0.891593719995219e-3,  4.057331392320188e-5,
			       1.052755021528803e-6,  1.557498087816203e-8,
			       1.313835773243312e-10, 6.319285885175346e-13,
			       1.733038792213266e-15, 2.709954036083074e-18 };

EXEC_TARGET devcomplex<fptype> device_Faddeeva_2 (const devcomplex<fptype>& z) {
  fptype *n,*e,t,u,r,s,d,f,g,h;
  devcomplex<fptype> c,d2,v;
  int i;

  s = norm2(z); // NB: norm2 is correct although CPU version calls the function 'norm'.
  if (s < 1e-7) {
    // use Pade approximation
    devcomplex<fptype> zz = z*z;
    v  = exp(zz); // Note lower-case! This is our own already-templated exp function for devcomplex, no need for float/double define.
    c  = C[0];
    d2 = D[0];
    for (i = 1; i <= 6; i++) {
      c  = c  * zz + C[i];
      d2 = d2 * zz + D[i];
    }
    return fptype(1.0) / v + devcomplex<fptype>(0.0,M_2_SQRTPI) * c/d2 * z * v;
  }

  // use trapezoid rule
  // select default table 1
  n = n1;
  e = e1;
  r = M_1_PI * 0.5;

  // if z is too close to a pole select table 2
  if (FABS(z.imag) < 0.01 && FABS(z.real) < 6.01) {
    // h = modf(2*FABS(z.real),&g);
    // Equivalent to above. Do this way because nvcc only knows about double version of modf.
    h = FABS(z.real)*2;
    g = FLOOR(h);
    h -= g;

    if (h < 0.02 || h > 0.98) {
      n = n2;
      e = e2;
      r = M_1_PI * 0.53;
    }
  }

  d = (z.imag - z.real) * (z.imag + z.real);
  f = 4 * z.real * z.real * z.imag * z.imag;

  g = h = 0.0;
  for (i = 0; i < 12; i++) {
    t = d + n[i];
    u = e[i] / (t * t + f);
    g += (s + n[i]) * u;
    h += (s - n[i]) * u;
  }
  u = 1 / s;

  c = r * devcomplex<fptype>(z.imag * (u + 2.0 * g),
			     z.real * (u + 2.0 * h) );

  if (z.imag < M_2PI) {
    s = 2.0 / r;
    t = s * z.real;
    u = s * z.imag;
    s = SIN(t);
    h = COS(t);
    f = EXP(- u) - h;
    g = 2.0 * EXP(d-u) / (s * s + f * f);
    u = 2.0 * z.real * z.imag;
    h = COS(u);
    t = SIN(u);
    c += g * devcomplex<fptype>( (h * f - t * s), -(h * s + t * f));
  }
  return c;
}

#else
EXEC_TARGET devcomplex<fptype> device_Faddeeva_2 (const devcomplex<fptype>& z) {
  fptype u,s,d,f,g,h;
  devcomplex<fptype> c,d2,v;

  s = norm2(z); // NB: norm2 is correct although CPU version calls the function 'norm'.
  if (s < 1e-7) {
    // use Pade approximation
    devcomplex<fptype> zz = z*z;
    v  = exp(zz); // Note lower-case! This is our own already-templated exp function for devcomplex, no need for float/double define.
    c  = C[0];
    d2 = D[0];
    for (int i = 1; i < 7; ++i) {
      c  = c  * zz + C[i];
      d2 = d2 * zz + D[i];
    }
    return fptype(1.0) / v + devcomplex<fptype>(0.0,M_2_SQRTPI) * c/d2 * z * v;
  }

  // use trapezoid rule
  fptype r = M_1_PI * 0.5;
  bool useDefault = true;

  // if z is too close to a pole select table 2
  if (FABS(z.imag) < 0.01 && FABS(z.real) < 6.01) {
    // h = modf(2*FABS(z.real),&g);
    // Equivalent to above. Do this way because nvcc only knows about double version of modf.
    h = FABS(z.real)*2;
    g = FLOOR(h);
    h -= g;

    if (h < 0.02 || h > 0.98) {
      useDefault = false;
      r = M_1_PI * 0.53;
    }
  }

  d = (z.imag - z.real) * (z.imag + z.real);
  f = 4 * z.real * z.real * z.imag * z.imag;

  g = h = 0.0;
  fptype currentN = (useDefault ? 0.25 : 0.2809);
  fptype currentE = (useDefault ? 0.7788007830714049 : 0.7551038420890235);
  fptype t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 1.0 : 1.1236);
  currentE = (useDefault ? 0.3678794411714423 : 0.3251072991205958);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 2.25 : 2.5281);
  currentE = (useDefault ? 1.053992245618643e-1  : 7.981051630007964e-2);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 4.0 : 4.4944);
  currentE = (useDefault ? 1.930454136227709e-3  : 0.891593719995219e-3);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 6.25 : 7.0225);
  currentE = (useDefault ? 4.785117392129009e-6  : 1.052755021528803e-6);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 9.0 : 10.1124);
  currentE = (useDefault ? 1.605228055185612e-9  : 1.313835773243312e-10);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 12.25 : 13.7641);
  currentE = (useDefault ? 7.287724095819692e-14 : 1.733038792213266e-15);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 16.0 : 17.9776);
  currentE = (useDefault ? 1.831563888873418e-2  : 1.117138143353082e-2);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 20.25 : 22.7529);
  currentE = (useDefault ? 1.234098040866795e-4  : 4.057331392320188e-5);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 25.0 : 28.09);
  currentE = (useDefault ? 1.125351747192591e-7  : 1.557498087816203e-8);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 30.25 : 33.9889);
  currentE = (useDefault ? 1.388794386496402e-11 : 6.319285885175346e-13);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  currentN = (useDefault ? 36.0 : 40.4496);
  currentE = (useDefault ? 2.319522830243569e-16 : 2.709954036083074e-18);
  t = d + currentN;
  u = currentE / (t*t + f);
  g += (s + currentN)*u;
  h += (s - currentN)*u;

  u = 1 / s;
  c = r * devcomplex<fptype>(z.imag * (u + 2.0 * g),
			     z.real * (u + 2.0 * h) );

  if (z.imag < M_2PI) {
    s = 2.0 / r;
    t = s * z.real;
    u = s * z.imag;
    s = SIN(t);
    h = COS(t);
    f = EXP(- u) - h;
    g = 2.0 * EXP(d-u) / (s * s + f * f);
    u = 2.0 * z.real * z.imag;
    h = COS(u);
    t = SIN(u);
    c += g * devcomplex<fptype>( (h * f - t * s), -(h * s + t * f));
  }

  return c;
}
#endif
EXEC_TARGET fptype device_Voigtian_Point (fptype* point, fptype* p, unsigned int* indices) {
  fptype x = point[0];

  fptype m = p[indices[1]];
  fptype w = p[indices[2]];
  fptype s = p[indices[3]];

#ifdef CUDAPRINT
  //if ((0 == THREADIDX) && (0 == BLOCKIDX))
  //cuPrintf("Values %f %i %i %i %f %f %f %i %i\n", x, indices[1], indices[2], indices[3], m, w, s, indices, callnumber);
  if (callnumber < 1) cuPrintf("Voigtian Values %f %i %i %i %f %f %f %i\n", x, indices[1], indices[2], indices[3], m, w, s, callnumber);

#endif

  // return constant for zero width and sigma
  if ((0==s) && (0==w)) return 1;
  //assert(s > 0); // No device-side assert?!
  //assert(w > 0);

  fptype arg = x - m;

  // Breit-Wigner for zero sigma
  if (0==s) return (1/(arg*arg+0.25*w*w));

  fptype coef = -0.5/(s*s);
  // Gauss for zero width
  if (0==w) return EXP(coef*arg*arg);

  // actual Voigtian for non-trivial width and sigma
  //fptype c = 1./(ROOT2*s);
  fptype c = 0.707106781187; // 1/root(2)
  c /= s;
  fptype a = 0.5*c*w;
  fptype u = c*arg;
  devcomplex<fptype> z(u,a) ;
  devcomplex<fptype> v = device_Faddeeva_2(z);

#define rsqrtPi 0.5641895835477563
  return c*rsqrtPi*v.real;
}

EXEC_TARGET fptype device_Voigtian (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x = evt[0];

  fptype m = p[indices[1]];
  fptype w = p[indices[2]];
  fptype s = p[indices[3]];

#ifdef CUDAPRINT
  //if ((0 == THREADIDX) && (0 == BLOCKIDX))
  //cuPrintf("Values %f %i %i %i %f %f %f %i %i\n", x, indices[1], indices[2], indices[3], m, w, s, indices, callnumber);
  if (callnumber < 1) cuPrintf("Voigtian Values %f %i %i %i %f %f %f %i\n", x, indices[1], indices[2], indices[3], m, w, s, callnumber);

#endif

  // return constant for zero width and sigma
  if ((0==s) && (0==w)) return 1;
  //assert(s > 0); // No device-side assert?!
  //assert(w > 0);

  fptype arg = x - m;

  // Breit-Wigner for zero sigma
  if (0==s) return (1/(arg*arg+0.25*w*w));

  fptype coef = -0.5/(s*s);
  // Gauss for zero width
  if (0==w) return EXP(coef*arg*arg);

  // actual Voigtian for non-trivial width and sigma
  //fptype c = 1./(ROOT2*s);
  fptype c = 0.707106781187; // 1/root(2)
  c /= s;
  fptype a = 0.5*c*w;
  fptype u = c*arg;
  devcomplex<fptype> z(u,a) ;
  devcomplex<fptype> v = device_Faddeeva_2(z);

#define rsqrtPi 0.5641895835477563
  return c*rsqrtPi*v.real;
}

EXEC_TARGET fptype device_VoigtianOffset (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x = evt[0];

  fptype m = p[indices[1]];
  fptype w = p[indices[2]];
  fptype s = p[indices[3]];
  fptype o = p[indices[4]];
  fptype lo = p[indices[5]];
  fptype up = p[indices[6]];
  
  //Thresholds
  if(x> up|| x<lo) return 0.0;
  
  //Offset
  x -= o;
  
#ifdef CUDAPRINT
  //if ((0 == THREADIDX) && (0 == BLOCKIDX))
  //cuPrintf("Values %f %i %i %i %f %f %f %i %i\n", x, indices[1], indices[2], indices[3], m, w, s, indices, callnumber);
  if (callnumber < 1) cuPrintf("Voigtian Values %f %i %i %i %f %f %f %i\n", x, indices[1], indices[2], indices[3], m, w, s, callnumber);

#endif

  // return constant for zero width and sigma
  if ((0==s) && (0==w)) return 1;
  //assert(s > 0); // No device-side assert?!
  //assert(w > 0);

  fptype arg = x - m;

  // Breit-Wigner for zero sigma
  if (0==s) return (1/(arg*arg+0.25*w*w));

  fptype coef = -0.5/(s*s);
  // Gauss for zero width
  if (0==w) return EXP(coef*arg*arg);

  // actual Voigtian for non-trivial width and sigma
  //fptype c = 1./(ROOT2*s);
  fptype c = 0.707106781187; // 1/root(2)
  c /= s;
  fptype a = 0.5*c*w;
  fptype u = c*arg;
  devcomplex<fptype> z(u,a) ;
  devcomplex<fptype> v = device_Faddeeva_2(z);

#define rsqrtPi 0.5641895835477563
  return c*rsqrtPi*v.real;
}

EXEC_TARGET fptype device_Voigtian_Bin (fptype* evt, fptype* p, unsigned int* indices) {

   //ONLY for binned datasets
  fptype xmid = evt[indices[2 + indices[0]]];
  fptype xnex = evt[indices[2 + indices[0]]+3]; // 3 = 1 observable (x) + binvalue + binvolume
  fptype bin = xnex-xmid;
  fptype step = bin; //Integration steps
  step /= 4;


  //STEPS FOR INTEGRATION
  fptype Xs[5];

  for(int l=0;l<=4;l++){
  Xs[l] = xmid +(-2+l)*step;
  }

  //NEWTON-COTES COEFFICIENTS
  fptype cnw = 2.0;
  cnw /=45.0;
  fptype b0 = 7.0;
  fptype b1 = 32.0;
  fptype b2 = 12.0;
  fptype b3 = 32.0;
  fptype b4 = 7.0;

  fptype a0 = b0*cnw;
  fptype a1 = b1*cnw;
  fptype a2 = b2*cnw;
  fptype a3 = b3*cnw;
  fptype a4 = b4*cnw;

  //PARAMETERS
  fptype m = p[indices[1]];
  fptype w = p[indices[2]];
  fptype s = p[indices[3]];

  fptype fx[5] = {1.0,1.0,1.0,1.0,1.0};

  for(int k=0;k<=4;k++){
  // return constant for zero width and sigma
  if ((0==s) && (0==w)) fx[k]=1;
  //assert(s > 0); // No device-side assert?!
  //assert(w > 0);

  fptype arg = Xs[k] - m;

  // Breit-Wigner for zero sigma
  if (0==s) fx[k]=(1/(arg*arg+0.25*w*w));

  fptype coef = -0.5/(s*s);
  // Gauss for zero width
  if (0==w) fx[k]=EXP(coef*arg*arg);

  // actual Voigtian for non-trivial width and sigma
  //fptype c = 1./(ROOT2*s);
  fptype c = 0.707106781187; // 1/root(2)
  c /= s;
  fptype a = 0.5*c*w;
  fptype u = c*arg;
  devcomplex<fptype> z(u,a) ;
  devcomplex<fptype> v = device_Faddeeva_2(z);

#define rsqrtPi 0.5641895835477563
  fx[k] = c*rsqrtPi*v.real;
  }

//INTEGRAL
  fptype integralF = a0*fx[0] + a1*fx[1] + a2*fx[2]  + a3*fx[3]  + a4*fx[4] ;
  integralF *= step;

/*
  printf("device_Vogitian %f %f\n", xmid, xnex);
  printf("device_Vogitian Bin NEWTON \n");
  printf("device_Vogitian Newton %f %f %f %f %f %f %f\n",cnw, a0, a1, a2, a3, a4,step);
  printf("device_Vogitian Bin VALUES \n");
  printf("device_Vogitian Values %f %f %f (xmid) %f %f %f %f %f \n", Xs[0],Xs[1],Xs[2],Xs[3],Xs[4], m, s,w);
  printf("device_Vogitian Bin FUNCTIONS \n");
  printf("device_Vogitian Functions %f %f %f f(xmid) %f %f \n", fx[0], fx[1],fx[2],fx[3],fx[4]);
  printf("device_Vogitian Bin INTEGRAL \n");
  printf("device_Vogitian Integral = (%f + %f + %f + %f + %f)*(%f)^-1 = %f \n",a0*fx[0], a1*fx[1], a2*fx[2], a3*fx[3], a4*fx[4], integralF);
*/

  return integralF;

}


MEM_DEVICE device_function_ptr ptr_to_Voigtian = device_Voigtian;
MEM_DEVICE device_function_ptr ptr_to_VoigtianOffset = device_VoigtianOffset;
MEM_DEVICE device_function_ptr ptr_to_Voigtian_Bin = device_Voigtian_Bin;
MEM_DEVICE device_function_ptr ptr_to_Voigtian_Point = device_Voigtian_Point;


__host__ VoigtianPdf::VoigtianPdf (std::string n, Variable* _x, Variable* m, Variable* s, Variable* w)
: GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(m));
  pindices.push_back(registerParameter(s));
  pindices.push_back(registerParameter(w));
  GET_FUNCTION_ADDR(ptr_to_Voigtian);
  GET_INTEGRAL_ADDR(ptr_to_Voigtian_Bin);
  GET_ATPOINTS_ADDR(ptr_to_Voigtian_Point);
  initialiseInt(pindices);
}

__host__ VoigtianPdf::VoigtianPdf (std::string n, Variable* _x, Variable* m, Variable* s,Variable* w,Variable* o, Variable* lowerCut, Variable* upperCut)
: GooPdf(_x, n)
,offset(o)
,lowerThreshold(lowerCut)
,upperThreshold(upperCut)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(m));
  pindices.push_back(registerParameter(s));
  pindices.push_back(registerParameter(w));
  pindices.push_back(registerParameter(o));
  pindices.push_back(registerParameter(lowerCut));
  pindices.push_back(registerParameter(upperCut));
  GET_FUNCTION_ADDR(ptr_to_VoigtianOffset);
  initialise(pindices);
}

