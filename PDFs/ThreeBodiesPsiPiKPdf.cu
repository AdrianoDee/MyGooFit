#include "ThreeBodiesPsiPiKPdf.hh"
#include "TRandom.h"

#define

fptype AREA[256] = {0, 0.00237088, 0.00649312, 0.0117822, 0.0180134, 0.0250543, 0.0328146, 0.0412273, 0.0502399, 0.05981, 0.0699022, 0.0804861, 0.0915358, 0.103028, 0.114943, 0.127261, 0.139967, 0.153046, 0.166483, 0.180265, 0.194382, 0.208821, 0.223573, 0.238629, 0.253978, 0.269612, 0.285524, 0.301705, 0.318149, 0.334848, 0.351796, 0.368986, 0.386412, 0.404069, 0.421951, 0.440052, 0.458367, 0.476891, 0.49562, 0.514548, 0.533671, 0.552985, 0.572485, 0.592168, 0.612028, 0.632063, 0.652268, 0.672639, 0.693174, 0.713867, 0.734717, 0.755719, 0.77687, 0.798168, 0.819608, 0.841187, 0.862903, 0.884753, 0.906733, 0.928841, 0.951073, 0.973428, 0.995902, 1.01849, 1.0412, 1.06401, 1.08694, 1.10997, 1.1331, 1.15634, 1.17967, 1.2031, 1.22663, 1.25024, 1.27395, 1.29774, 1.32162, 1.34558, 1.36962, 1.39374, 1.41794, 1.4422, 1.46655, 1.49096, 1.51543, 1.53998, 1.56458, 1.58925, 1.61398, 1.63877, 1.66361, 1.6885, 1.71345, 1.73845, 1.76349, 1.78858, 1.81371, 1.83889, 1.86411, 1.88936, 1.91465, 1.93998, 1.96534, 1.99073, 2.01615, 2.0416, 2.06707, 2.09257, 2.11809, 2.14363, 2.16919, 2.19477, 2.22036, 2.24597, 2.27159, 2.29721, 2.32285, 2.3485, 2.37414, 2.3998, 2.42545, 2.4511, 2.47676, 2.50241, 2.52805, 2.55369, 2.57932, 2.60494, 2.63054, 2.65614, 2.68172, 2.70728, 2.73283, 2.75835, 2.78385, 2.80933, 2.83479, 2.86022, 2.88562, 2.91099, 2.93633, 2.96163, 2.9869, 3.01214, 3.03733, 3.06249, 3.08761, 3.11268, 3.1377, 3.16268, 3.18762, 3.2125, 3.23733, 3.2621, 3.28683, 3.31149, 3.3361, 3.36064, 3.38513, 3.40955, 3.4339, 3.45819, 3.4824, 3.50655, 3.53063, 3.55463, 3.57855, 3.60239, 3.62616, 3.64984, 3.67344, 3.69696, 3.72039, 3.74372, 3.76697, 3.79012, 3.81318, 3.83614, 3.859, 3.88176, 3.90442, 3.92697, 3.94942, 3.97175, 3.99398, 4.01608, 4.03808, 4.05995, 4.08171, 4.10334, 4.12485, 4.14623, 4.16748, 4.1886, 4.20959, 4.23044, 4.25114, 4.27171, 4.29213, 4.31241, 4.33254, 4.35251, 4.37233, 4.39199, 4.41148, 4.43082, 4.44999, 4.46898, 4.48781, 4.50645, 4.52492, 4.5432, 4.5613, 4.57921, 4.59692, 4.61443, 4.63174, 4.64884, 4.66574, 4.68242, 4.69888, 4.71511, 4.73112, 4.74689, 4.76242, 4.77771, 4.79275, 4.80753, 4.82204, 4.83629, 4.85026, 4.86394, 4.87733, 4.89042, 4.9032, 4.91566, 4.92779, 4.93957, 4.951, 4.96206, 4.97274, 4.98302, 4.99288, 5.0023, 5.01126, 5.01972, 5.02766, 5.03504, 5.04179, 5.04787, 5.05317, 5.05755, 5.06075, 5.06184, 5.06184, 5.06184};

fptype PUNTI[256] = {1.019, 1.02165, 1.0243, 1.02694, 1.02959, 1.03224, 1.03489, 1.03753, 1.04018, 1.04283, 1.04548, 1.04813, 1.05077, 1.05342, 1.05607, 1.05872, 1.06137, 1.06401, 1.06666, 1.06931, 1.07196, 1.0746, 1.07725, 1.0799, 1.08255, 1.0852, 1.08784, 1.09049, 1.09314, 1.09579, 1.09844, 1.10108, 1.10373, 1.10638, 1.10903, 1.11167, 1.11432, 1.11697, 1.11962, 1.12227, 1.12491, 1.12756, 1.13021, 1.13286, 1.13551, 1.13815, 1.1408, 1.14345, 1.1461, 1.14874, 1.15139, 1.15404, 1.15669, 1.15934, 1.16198, 1.16463, 1.16728, 1.16993, 1.17257, 1.17522, 1.17787, 1.18052, 1.18317, 1.18581, 1.18846, 1.19111, 1.19376, 1.19641, 1.19905, 1.2017, 1.20435, 1.207, 1.20964, 1.21229, 1.21494, 1.21759, 1.22024, 1.22288, 1.22553, 1.22818, 1.23083, 1.23348, 1.23612, 1.23877, 1.24142, 1.24407, 1.24671, 1.24936, 1.25201, 1.25466, 1.25731, 1.25995, 1.2626, 1.26525, 1.2679, 1.27055, 1.27319, 1.27584, 1.27849, 1.28114, 1.28378, 1.28643, 1.28908, 1.29173, 1.29438, 1.29702, 1.29967, 1.30232, 1.30497, 1.30761, 1.31026, 1.31291, 1.31556, 1.31821, 1.32085, 1.3235, 1.32615, 1.3288, 1.33145, 1.33409, 1.33674, 1.33939, 1.34204, 1.34468, 1.34733, 1.34998, 1.35263, 1.35528, 1.35792, 1.36057, 1.36322, 1.36587, 1.36852, 1.37116, 1.37381, 1.37646, 1.37911, 1.38175, 1.3844, 1.38705, 1.3897, 1.39235, 1.39499, 1.39764, 1.40029, 1.40294, 1.40559, 1.40823, 1.41088, 1.41353, 1.41618, 1.41882, 1.42147, 1.42412, 1.42677, 1.42942, 1.43206, 1.43471, 1.43736, 1.44001, 1.44265, 1.4453, 1.44795, 1.4506, 1.45325, 1.45589, 1.45854, 1.46119, 1.46384, 1.46649, 1.46913, 1.47178, 1.47443, 1.47708, 1.47972, 1.48237, 1.48502, 1.48767, 1.49032, 1.49296, 1.49561, 1.49826, 1.50091, 1.50356, 1.5062, 1.50885, 1.5115, 1.51415, 1.51679, 1.51944, 1.52209, 1.52474, 1.52739, 1.53003, 1.53268, 1.53533, 1.53798, 1.54063, 1.54327, 1.54592, 1.54857, 1.55122, 1.55386, 1.55651, 1.55916, 1.56181, 1.56446, 1.5671, 1.56975, 1.5724, 1.57505, 1.57769, 1.58034, 1.58299, 1.58564, 1.58829, 1.59093, 1.59358, 1.59623, 1.59888, 1.60153, 1.60417, 1.60682, 1.60947, 1.61212, 1.61476, 1.61741, 1.62006, 1.62271, 1.62536, 1.628, 1.63065, 1.6333, 1.63595, 1.6386, 1.64124, 1.64389, 1.64654, 1.64919, 1.65183, 1.65448, 1.65713, 1.65978, 1.66243, 1.66507, 1.66772, 1.67037, 1.67302, 1.67567, 1.67831, 1.68096, 1.68361, 1.68626, 1.6889, 1.69155, 1.6942};

EXEC_TARGET fptype device_ThreeBodiesPsiPiK (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x = evt[indices[2 + indices[0]]];

  fptype mP = p[indices[1]];
  fptype m1 = p[indices[2]];
  fptype m2 = p[indices[3]];
  fptype m3 = p[indices[4]];

  fptype ret = isnan(sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x)) ? 0 : (sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x));


return ret;
}

EXEC_TARGET fptype device_ThreeBodiesPsiPiK_Point (fptype* point, fptype* p, unsigned int* indices) {
  fptype x = point[0];

  fptype mP = p[indices[1]];
  fptype m1 = p[indices[2]];
  fptype m2 = p[indices[3]];
  fptype m3 = p[indices[4]];

  fptype ret = isnan(sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x)) ? 0 : (sqrt(pow(x,4) + pow(m1,4) + pow(m2,4) - 2*pow(x,2)*pow(m1,2) - 2*pow(x,2)*pow(m2,2) - 2*pow(m1,2)*pow(m2,2)) * sqrt(pow(mP,4) + pow(x,4) + pow(m3,4) - 2*pow(mP,2)*pow(x,2) - 2*pow(mP,2)*pow(m3,2) - 2*pow(x,2)*pow(m3,2) ) / (x));

  return ret;
}

/*
EXEC_TARGET fptype device_Three_Bin (fptype* evt, fptype* p, unsigned int* indices) {
  fptype xmid = evt[indices[2 + indices[0]]];
  fptype xnex = evt[indices[2 + indices[0]]+3];

  fptype bin = xnex-xmid;
  fptype step = bin; //Integration steps
  int n = 4;
  step /= n;

  //STEPS FOR INTEGRATION
  fptype Xs[5];

  for (int k=0;k<=n;k++){

    Xs[k]=(-2+k)*step+xmid;

  }

  //NEWTON-COTES COEFFICIENTS
  fptype cnw = 2.0;
  cnw /=45.0;
  fptype b[5] = {7.0,32.0,12.0,32.0,7.0};

  fptype a[5];

  for (int k=0;k<=n;k++){
        a[k] = b[k]*cnw;
  }

  //FUNCTION IN Xi
  fptype fX[5];
  for (int k=0;k<=n;k++){
    fptype x;
    x=Xs[k];

    fX[k] = (x<1.015)||(x>1.7)? 0 : isnan((sqrt( pow(x+3.0967,4) + pow(3.0967,4) + pow(1.01946,4) - 2*pow(x+3.0967,2)*pow(3.0967,2) - 2*pow(3.0967,2)*pow(1.01946,2) - 2*pow(x+3.0967,2)*pow(1.01946,2) ) * sqrt( pow(5.279,4) + pow(x+3.0967,4) + pow(0.493677,4) - 2*pow(5.279,2)*pow(x+3.0967,2) - 2*pow(5.279,2)*pow(0.493677,2) - 2*pow(x+3.0967,2)*pow(0.493677,2) ) / (x+3.0967)))? 0 : (sqrt( pow(x+3.0967,4) + pow(3.0967,4) + pow(1.01946,4) - 2*pow(x+3.0967,2)*pow(3.0967,2) - 2*pow(3.0967,2)*pow(1.01946,2) - 2*pow(x+3.0967,2)*pow(1.01946,2) ) * sqrt( pow(5.279,4) + pow(x+3.0967,4) + pow(0.493677,4) - 2*pow(5.279,2)*pow(x+3.0967,2) - 2*pow(5.279,2)*pow(0.493677,2) - 2*pow(x+3.0967,2)*pow(0.493677,2) ) / (x+3.0967));

  }

  //INTEGRAL
  fptype integralF = 0;

   for (int k=0;k<=n;k++){
        integralF += a[k]*fX[k];
   }
    integralF *= step;

return integralF;

}
*/

MEM_DEVICE device_function_ptr ptr_to_ThreeBodiesPsiPiK = device_Three;
//MEM_DEVICE device_function_ptr ptr_to_ThreeBodiesPsiPiK_Bin = device_Three_Bin;
MEM_DEVICE device_function_ptr ptr_to_ThreeBodiesPsiPiK_Point = device_ThreeBodiesPsiPiK_Point;


__host__ ThreeBodiesPsiPiK::ThreeBodiesPsiPiK (std::string n, Variable* _x)
  : GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;

  GET_FUNCTION_ADDR(ptr_to_ThreeBodiesPsiPiK);
  //GET_INTEGRAL_ADDR(ptr_to_Three_Bin);
  GET_ATPOINTS_ADDR(ptr_to_ThreeBodiesPsiPiK_Point);
  initialiseInt(pindices);
}

__host__ fptype ThreeBodiesPsiPiK::integrate (fptype lo, fptype hi) const {

  int loind=0, hiind=0;
  fptype max =2.0;
  fptype min =1.0;

  if(lo<=1.0){
  if(hi>=2) return 5.06186;
  else {
  while(hi<max){
  max=PUNTI[255-hiind];
  hiind++;
  }
  return AREA[255-hiind];
  }}

  if(hi>=2.0){
  if(lo<=1) return 5.06186;
  else {
  while(lo>min){
  min=PUNTI[loind];
  loind++;
  }
  return AREA[loind];
  }}

  while(lo>min){
  min=PUNTI[loind];
  loind++;
  }

  while(hi<max){
  max=PUNTI[255-hiind];
  hiind++;
  }

  return AREA[255-hiind]-AREA[loind];
}

__global__ void fillHisto (TH1F* histo,unsigned int nevents,fptype fmax){

/*
    fptype roll,func,x;
    fptype xmin = histo->GetBinLowEdge(1);
    fptype xmax = histo->GetBinLowEdge(histo->GetNbinsX()+1);
    long int ms; struct timeval tp;

    for (int j = 0; j < nevents; ++j) {

		    gettimeofday(&tp,NULL);
        ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
        TRandom3 ranGen(ms);

  			x = donram.Uniform(xmax-xmins)+xmin;
  			func = background(xvar->value);
  			roll = donram.Uniform(fmax);

  			if (roll > func) {
  				--j;
  				continue; }

  			if ((x < xmin) || (x > xmax)) {
  				--j;
  				continue;}

  			histo->Fill(x);

  		}

*/
}
