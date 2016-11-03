#ifdef OMP_ON
#include <omp.h>
#endif

PdfBase* pdfPointer;
FitManager* currGlue = 0;
int numPars = 0;
vector<Variable*> vars;
int fitMan = 0;

void specialTddpPrint (double fun);

FitManager::FitManager (PdfBase* dat)
  : minuit(0)
  , overrideCallLimit(-1),runhesse(false),runminos(false)
{
  pdfPointer = dat;
  currGlue = this;
}

FitManager::FitManager (PdfBase* dat,bool hesse,bool minos)
  : minuit(0)
  , overrideCallLimit(-1),runhesse(hesse),runminos(minos)
{
  pdfPointer = dat;
  currGlue = this;
}

FitManager::~FitManager () {
  if (minuit) delete minuit;
}

void FitManager::setupMinuit () {
  //printf("SetUp Minuit %d \n",fitMan); fitMan++;
  vars.clear();
  pdfPointer->getParameters(vars);
  //printf("SetUp Minuit %d \n",fitMan); fitMan++;
  numPars = vars.size();
  if (minuit) delete minuit;
  minuit = new TMinuit(numPars);
  int maxIndex = 0;
  int counter = 0;
  for (std::vector<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    minuit->DefineParameter(counter, (*i)->name.c_str(), (*i)->value, (*i)->error, (*i)->lowerlimit, (*i)->upperlimit);
    if ((*i)->fixed) minuit->FixParameter(counter);
    counter++;
    if (maxIndex < (*i)->getIndex()) maxIndex = (*i)->getIndex();
  }
  //printf("SetUp Minuit %d \n",fitMan); fitMan++;
  numPars = maxIndex+1;
  pdfPointer->copyParams();
  minuit->SetFCN(FitFun);
  //minuit->SetPrintLevel(-1);
}

void FitManager::fit () {

   setupMinuit();
   runMigrad();

   if(runhesse) runHesse();
   if(runminos) runMinos();

}

void FitManager::fitOrderd (std::vector< std::string > algos) {

   setupMinuit();

   std::cout<<" ===================== Setting Up Minimisation Algos "<<st::endl;
   for (size_t i = 0; i < algos.size(); i++) {
     std::cout<<" Running : "<<std::endl;
     if(algos[i] == "MIGRAD") std::cout<<"\t\t - Migrad"<<std::endl;
     else if(algos[i] == "MIGRAD") std::cout<<"\t\t - Migrad"<<std::endl;
     else if(algos[i] == "MIGRAD") std::cout<<"\t\t - Migrad"<<std::endl;
     else {
       std::cout<<" INVALID ALGO INPUT ================================= "<<std::endl;
       std::cout<<" Options: \"MIGRAD\" - \"HESSE\" - \"MINOS\" "std::endl;
       exit;
     }
   }

   for (size_t i = 0; i < algos.size(); i++) {
     if(algos[i] == "MIGRAD") runMigrad();
     else if(algos[i] == "MIGRAD") runHesse();
     else if(algos[i] == "MIGRAD") runMinos();
     else {
       std::cout<<" INVALID ALGO INPUT ================================= "<<std::endl;
       std::cout<<" Options: \"MIGRAD\" - \"HESSE\" - \"MINOS\" "std::endl;
     }
   }

}

void FitManager::runMigrad () {
  //printf("Run Migrad %d \n",fitMan); fitMan++;
  assert(minuit);
  host_callnumber = 0;
  if (0 < overrideCallLimit) {
    std::cout << "Calling MIGRAD with call limit " << overrideCallLimit << std::endl;
    double plist[1];
    plist[0] = overrideCallLimit;
    int err = 0;
    minuit->mnexcm("MIGRAD", plist, 1, err);
  }
  else{

//   double arglist[10];
//   arglist [0]=2;
//   int ierflg=0;
//   minuit->mnexcm("SET STR",arglist,1,ierflg);

   minuit->Migrad();

  }

}

void FitManager::runHesse () {
  //printf("Run Migrad %d \n",fitMan); fitMan++;
  assert(minuit);
  host_callnumber = 0;
  if (0 < overrideCallLimit) {
    std::cout << "Calling HESSE with call limit " << overrideCallLimit << std::endl;
    double plist[1];
    plist[0] = overrideCallLimit;
    int err = 0;
    minuit->mnexcm("HESSE", plist, 1, err);
  }
  else{
    int err;
    double tmp[1];
    tmp[0] = 0;
    minuit->mnexcm("HESSE", tmp, 1, err);

  }

}

void FitManager::runMinos () {
  //printf("Run Minos %d \n",fitMan); fitMan++;
  assert(minuit);
  host_callnumber = 0;
  if (0 < overrideCallLimit) {
    std::cout << "Calling HESSE with call limit " << overrideCallLimit << std::endl;
    double plist[1];
    plist[0] = overrideCallLimit;
    int err = 0;
    minuit->mnexcm("MINOS", plist, 1, err);
  }
  else{
    int err;
    double tmp[1];
    tmp[0] = 0;
    minuit->mnexcm("MINOS", tmp, 1, err);

  }

}

void FitManager::getMinuitValues () const {
  int counter = 0;
  for (std::vector<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    minuit->GetParameter(counter++, (*i)->value, (*i)->error);
  }
}

void FitFun(int &npar, double *gin, double &fun, double *fp, int iflag) {
  //printf("FitFun %d \n",fitMan); fitMan++;
  vector<double> pars;
  // Notice that npar is number of variable parameters, not total.
  pars.resize(numPars);
  //printf("FitFun %d \n",fitMan); fitMan++;
  int counter = 0;

  for (std::vector<Variable*>::iterator i = vars.begin(); i != vars.end(); ++i) {
    if (isnan(fp[counter])) cout << "Variable " << (*i)->name << " " << (*i)->index << " is NaN\n";

    pars[(*i)->getIndex()] = fp[counter++] + (*i)->blind;

 }
  //printf("FitFun %d \n",fitMan); fitMan++;
  pdfPointer->copyParams(pars);
  //printf("FitFun %d \n",fitMan); fitMan++;
  fun = pdfPointer->calculateNLL();
  //printf("FitFun %d \n",fitMan); fitMan++;
  host_callnumber++;

#ifdef PRINTCALLS
  specialTddpPrint(fun);
#endif
}

#ifdef PRINTCALLS
void specialTddpPrint (double fun) {
  // Stupid amplitude-fit debugging method.
  cout << "Function call " << host_callnumber << ": " << fun << "\n";
  currGlue->getMinuitValues();
  int varCount = 1;
  for (std::vector<Variable*>::iterator v = vars.begin(); v != vars.end(); ++v) {
    if (!(*v)) cout << "Null!" << endl;
    if ((*v)->fixed) continue;

    const fptype _mD0 = 1.86484;
    const fptype _mD02 = _mD0 *_mD0;
    const fptype _mD02inv = 1./_mD02;
    double stupidSpecialModifier = 1; // Mikhail interprets some of the weights differently.
    if (((*v)->name == "f0_980_amp_real") ||
	((*v)->name == "f0_980_amp_imag") ||
	((*v)->name == "f0_1370_amp_real") ||
	((*v)->name == "f0_1370_amp_imag") ||
	((*v)->name == "f0_1500_amp_real") ||
	((*v)->name == "f0_1500_amp_imag") ||
	((*v)->name == "f0_1710_amp_real") ||
	((*v)->name == "f0_1710_amp_imag") ||
	((*v)->name == "f0_600_amp_real") ||
	((*v)->name == "f0_600_amp_imag")) stupidSpecialModifier = -_mD02;
    else if (((*v)->name == "f2_1270_amp_real") ||
	     ((*v)->name == "f2_1270_amp_imag")) stupidSpecialModifier = -_mD02inv;
    else if (((*v)->name == "nonr_amp_real") ||
	     ((*v)->name == "nonr_amp_imag")) stupidSpecialModifier = -1;

    cout.width(20);
    cout << (*v)->name;
    cout.setf(ios_base::right,ios_base::adjustfield);
    cout.width(3);
    cout << varCount++;
    cout.setf(ios_base::right,ios_base::adjustfield); cout.precision(8);
    cout << "  ";         cout.width(12);
    cout << (*v)->value / stupidSpecialModifier;
    cout.setf(ios_base::right,ios_base::adjustfield); cout.precision(8);
    cout << "  ";         cout.width(12);
    cout << (*v)->error;
    cout << endl;
  }

  cout << endl;
}
#endif
