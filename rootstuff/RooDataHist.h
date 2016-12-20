/*****************************************************************************
    2  * Project: RooFit                                                           *
    3  * Package: RooFitCore                                                       *
    4  * @(#)root/roofitcore:$Id$
    5  * Authors:                                                                  *
    6  *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
    7  *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
    8  *                                                                           *
    9  * Copyright (c) 2000-2005, Regents of the University of California          *
   10  *                          and Stanford University. All rights reserved.    *
   11  *                                                                           *
   12  * Redistribution and use in source and binary forms,                        *
   13  * with or without modification, are permitted according to the terms        *
   14  * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
   15  *****************************************************************************/
   16 
   17 /**
   18 \file RooDataHist.cxx
   19 \class RooDataHist
   20 \ingroup Roofitcore
   21 
   22 RooDataSet is a container class to hold N-dimensional binned data. Each bins central
   23 coordinates in N-dimensional space are represented by a RooArgSet of RooRealVar, RooCategory
   24 or RooStringVar objects, thus data can be binned in real and/or discrete dimensions
   25 **/
   26 
   27 #include "RooFit.h"
   28 #include "Riostream.h"
   29 
   30 #include "TH1.h"
   31 #include "TH1.h"
   32 #include "TDirectory.h"
   33 #include "TMath.h"
   34 #include "RooMsgService.h"
   35 #include "RooDataHist.h"
   36 #include "RooDataHistSliceIter.h"
   37 #include "RooAbsLValue.h"
   38 #include "RooArgList.h"
   39 #include "RooRealVar.h"
   40 #include "RooMath.h"
   41 #include "RooBinning.h"
   42 #include "RooPlot.h"
   43 #include "RooHistError.h"
   44 #include "RooCategory.h"
   45 #include "RooCmdConfig.h"
   46 #include "RooLinkedListIter.h"
   47 #include "RooTreeDataStore.h"
   48 #include "RooVectorDataStore.h"
   49 #include "TTree.h"
   50 #include "RooTrace.h"
   51 #include "RooTreeData.h"
   52 
   53 using namespace std ;
   54 
   55 ClassImp(RooDataHist)
   56 ;
   57 
   58 
   59 
   60 ////////////////////////////////////////////////////////////////////////////////
   61 /// Default constructor
   62 
   63 RooDataHist::RooDataHist() : _pbinvCacheMgr(0,10)
   64 {
   65   _arrSize = 0 ;
   66   _wgt = 0 ;
   67   _errLo = 0 ;
   68   _errHi = 0 ;
   69   _sumw2 = 0 ;
   70   _binv = 0 ;
   71   _pbinv = 0 ;
   72   _curWeight = 0 ;
   73   _curIndex = -1 ;
   74   _realIter = _realVars.createIterator() ;
   75   _binValid = 0 ;
   76   _curSumW2 = 0 ;
   77   _curVolume = 1 ;
   78   _curWgtErrHi = 0 ;
   79   _curWgtErrLo = 0 ;
   80   _cache_sum_valid = 0 ;
   81   TRACE_CREATE
   82 }
   83 
   84 
   85 
   86 ////////////////////////////////////////////////////////////////////////////////
   87 /// Constructor of an empty data hist from a RooArgSet defining the dimensions
   88 /// of the data space. The range and number of bins in each dimensions are taken
   89 /// from getMin()getMax(),getBins() of each RooAbsArg representing that
   90 /// dimension.
   91 ///
   92 /// For real dimensions, the fit range and number of bins can be set independently
   93 /// of the plot range and number of bins, but it is advisable to keep the
   94 /// ratio of the plot bin width and the fit bin width an integer value.
   95 /// For category dimensions, the fit ranges always comprises all defined states
   96 /// and each state is always has its individual bin
   97 ///
   98 /// To effective achive binning of real dimensions with variable bins sizes,
   99 /// construct a RooThresholdCategory of the real dimension to be binned variably.
  100 /// Set the thresholds at the desired bin boundaries, and construct the
  101 /// data hist as function of the threshold category instead of the real variable.
  102 
  103 RooDataHist::RooDataHist(const char *name, const char *title, const RooArgSet& vars, const char* binningName) :
  104   RooAbsData(name,title,vars), _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
  105 {
  106   // Initialize datastore
  107   _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) :
  108                                          ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;
  109 
  110   initialize(binningName) ;
  111 
  112   _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  113 
  114   appendToDir(this,kTRUE) ;
  115   TRACE_CREATE
  116 }
  117 
  118 
  119 
  120 ////////////////////////////////////////////////////////////////////////////////
  121 /// Constructor of a data hist from an existing data collection (binned or unbinned)
  122 /// The RooArgSet 'vars' defines the dimensions of the histogram.
  123 /// The range and number of bins in each dimensions are taken
  124 /// from getMin()getMax(),getBins() of each RooAbsArg representing that
  125 /// dimension.
  126 ///
  127 /// For real dimensions, the fit range and number of bins can be set independently
  128 /// of the plot range and number of bins, but it is advisable to keep the
  129 /// ratio of the plot bin width and the fit bin width an integer value.
  130 /// For category dimensions, the fit ranges always comprises all defined states
  131 /// and each state is always has its individual bin
  132 ///
  133 /// To effective achive binning of real dimensions with variable bins sizes,
  134 /// construct a RooThresholdCategory of the real dimension to be binned variably.
  135 /// Set the thresholds at the desired bin boundaries, and construct the
  136 /// data hist as function of the threshold category instead of the real variable.
  137 ///
  138 /// If the constructed data hist has less dimensions that in source data collection,
  139 /// all missing dimensions will be projected.
  140 
  141 RooDataHist::RooDataHist(const char *name, const char *title, const RooArgSet& vars, const RooAbsData& data, Double_t wgt) :
  142   RooAbsData(name,title,vars), _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
  143 {
  144   // Initialize datastore
  145   _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) :
  146                                          ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;
  147 
  148   initialize() ;
  149   _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  150 
  151   add(data,(const RooFormulaVar*)0,wgt) ;
  152   appendToDir(this,kTRUE) ;
  153   TRACE_CREATE
  154 }
  155 
  156 
  157 
  158 ////////////////////////////////////////////////////////////////////////////////
  159 /// Constructor of a data hist from a map of TH1,TH2 or TH3 that are collated into a x+1 dimensional
  160 /// RooDataHist where the added dimension is a category that labels the input source as defined
  161 /// in the histMap argument. The state names used in histMap must correspond to predefined states
  162 /// 'indexCat'
  163 ///
  164 /// The RooArgList 'vars' defines the dimensions of the histogram.
  165 /// The ranges and number of bins are taken from the input histogram and must be the same in all histograms
  166 
  167 RooDataHist::RooDataHist(const char *name, const char *title, const RooArgList& vars, RooCategory& indexCat,
  168           map<string,TH1*> histMap, Double_t wgt) :
  169   RooAbsData(name,title,RooArgSet(vars,&indexCat)),
  170   _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
  171 {
  172   // Initialize datastore
  173   _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) :
  174                                          ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;
  175 
  176   importTH1Set(vars, indexCat, histMap, wgt, kFALSE) ;
  177 
  178   _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  179   TRACE_CREATE
  180 }
  181 
  182 
  183 
  184 ////////////////////////////////////////////////////////////////////////////////
  185 /// Constructor of a data hist from a map of RooDataHists that are collated into a x+1 dimensional
  186 /// RooDataHist where the added dimension is a category that labels the input source as defined
  187 /// in the histMap argument. The state names used in histMap must correspond to predefined states
  188 /// 'indexCat'
  189 ///
  190 /// The RooArgList 'vars' defines the dimensions of the histogram.
  191 /// The ranges and number of bins are taken from the input histogram and must be the same in all histograms
  192 
  193 RooDataHist::RooDataHist(const char *name, const char *title, const RooArgList& vars, RooCategory& indexCat,
  194           map<string,RooDataHist*> dhistMap, Double_t wgt) :
  195   RooAbsData(name,title,RooArgSet(vars,&indexCat)),
  196   _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
  197 {
  198   // Initialize datastore
  199   _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) :
  200                                          ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;
  201 
  202   importDHistSet(vars, indexCat, dhistMap, wgt) ;
  203 
  204   _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  205   TRACE_CREATE
  206 }
  207 
  208 
  209 
  210 ////////////////////////////////////////////////////////////////////////////////
  211 /// Constructor of a data hist from an TH1,TH2 or TH3
  212 /// The RooArgSet 'vars' defines the dimensions of the histogram. The ranges
  213 /// and number of bins are taken from the input histogram, and the corresponding
  214 /// values are set accordingly on the arguments in 'vars'
  215 
  216 RooDataHist::RooDataHist(const char *name, const char *title, const RooArgList& vars, const TH1* hist, Double_t wgt) :
  217   RooAbsData(name,title,vars), _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
  218 {
  219   // Initialize datastore
  220   _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) :
  221                                          ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;
  222 
  223   // Check consistency in number of dimensions
  224   if (vars.getSize() != hist->GetDimension()) {
  225     coutE(InputArguments) << "RooDataHist::ctor(" << GetName() << ") ERROR: dimension of input histogram must match "
  226            << "number of dimension variables" << endl ;
  227     assert(0) ;
  228   }
  229 
  230   importTH1(vars,*hist,wgt, kFALSE) ;
  231 
  232   _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  233   TRACE_CREATE
  234 }
  235 
  236 
  237 
  238 ////////////////////////////////////////////////////////////////////////////////
  239 /// Constructor of a binned dataset from a RooArgSet defining the dimensions
  240 /// of the data space. The range and number of bins in each dimensions are taken
  241 /// from getMin() getMax(),getBins() of each RooAbsArg representing that
  242 /// dimension.
  243 ///
  244 /// This constructor takes the following optional arguments
  245 ///
  246 /// Import(TH1&, Bool_t impDens) -- Import contents of the given TH1/2/3 into this binned dataset. The
  247 ///                                 ranges and binning of the binned dataset are automatically adjusted to
  248 ///                                 match those of the imported histogram.
  249 ///
  250 ///                                 Please note: for TH1& with unequal binning _only_,
  251 ///                                 you should decide if you want to import the absolute bin content,
  252 ///                                 or the bin content expressed as density. The latter is default and will
  253 ///                                 result in the same histogram as the original TH1. For certain type of
  254 ///                                 bin contents (containing efficiencies, asymmetries, or ratio is general)
  255 ///                                 you should import the absolute value and set impDens to kFALSE
  256 ///
  257 ///
  258 /// Weight(Double_t)          -- Apply given weight factor when importing histograms
  259 ///
  260 /// Index(RooCategory&)       -- Prepare import of multiple TH1/1/2/3 into a N+1 dimensional RooDataHist
  261 ///                              where the extra discrete dimension labels the source of the imported histogram
  262 ///                              If the index category defines states for which no histogram is be imported
  263 ///                              the corresponding bins will be left empty.
  264 ///
  265 /// Import(const char*, TH1&) -- Import a THx to be associated with the given state name of the index category
  266 ///                              specified in Index(). If the given state name is not yet defined in the index
  267 ///                              category it will be added on the fly. The import command can be specified
  268 ///                              multiple times.
  269 /// Import(map<string,TH1*>&) -- As above, but allows specification of many imports in a single operation
  270 ///
  271 
  272 RooDataHist::RooDataHist(const char *name, const char *title, const RooArgList& vars, const RooCmdArg& arg1, const RooCmdArg& arg2, const RooCmdArg& arg3,
  273           const RooCmdArg& arg4,const RooCmdArg& arg5,const RooCmdArg& arg6,const RooCmdArg& arg7,const RooCmdArg& arg8) :
  274   RooAbsData(name,title,RooArgSet(vars,(RooAbsArg*)RooCmdConfig::decodeObjOnTheFly("RooDataHist::RooDataHist", "IndexCat",0,0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8))),
  275   _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
  276 {
  277   // Initialize datastore
  278   _dstore = (defaultStorageType==Tree) ? ((RooAbsDataStore*) new RooTreeDataStore(name,title,_vars)) :
  279                                          ((RooAbsDataStore*) new RooVectorDataStore(name,title,_vars)) ;
  280 
  281   // Define configuration for this method
  282   RooCmdConfig pc(Form("RooDataHist::ctor(%s)",GetName())) ;
  283   pc.defineObject("impHist","ImportHisto",0) ;
  284   pc.defineInt("impDens","ImportHisto",0) ;
  285   pc.defineObject("indexCat","IndexCat",0) ;
  286   pc.defineObject("impSliceHist","ImportHistoSlice",0,0,kTRUE) ; // array
  287   pc.defineString("impSliceState","ImportHistoSlice",0,"",kTRUE) ; // array
  288   pc.defineObject("impSliceDHist","ImportDataHistSlice",0,0,kTRUE) ; // array
  289   pc.defineString("impSliceDState","ImportDataHistSlice",0,"",kTRUE) ; // array
  290   pc.defineDouble("weight","Weight",0,1) ;
  291   pc.defineObject("dummy1","ImportDataHistSliceMany",0) ;
  292   pc.defineObject("dummy2","ImportHistoSliceMany",0) ;
  293   pc.defineMutex("ImportHisto","ImportHistoSlice","ImportDataHistSlice") ;
  294   pc.defineDependency("ImportHistoSlice","IndexCat") ;
  295   pc.defineDependency("ImportDataHistSlice","IndexCat") ;
  296 
  297   RooLinkedList l ;
  298   l.Add((TObject*)&arg1) ;  l.Add((TObject*)&arg2) ;
  299   l.Add((TObject*)&arg3) ;  l.Add((TObject*)&arg4) ;
  300   l.Add((TObject*)&arg5) ;  l.Add((TObject*)&arg6) ;
  301   l.Add((TObject*)&arg7) ;  l.Add((TObject*)&arg8) ;
  302 
  303   // Process & check varargs
  304   pc.process(l) ;
  305   if (!pc.ok(kTRUE)) {
  306     assert(0) ;
  307     return ;
  308   }
  309 
  310   TH1* impHist = static_cast<TH1*>(pc.getObject("impHist")) ;
  311   Bool_t impDens = pc.getInt("impDens") ;
  312   Double_t initWgt = pc.getDouble("weight") ;
  313   const char* impSliceNames = pc.getString("impSliceState","",kTRUE) ;
  314   const RooLinkedList& impSliceHistos = pc.getObjectList("impSliceHist") ;
  315   RooCategory* indexCat = static_cast<RooCategory*>(pc.getObject("indexCat")) ;
  316   const char* impSliceDNames = pc.getString("impSliceDState","",kTRUE) ;
  317   const RooLinkedList& impSliceDHistos = pc.getObjectList("impSliceDHist") ;
  318 
  319 
  320   if (impHist) {
  321 
  322     // Initialize importing contents from TH1
  323     importTH1(vars,*impHist,initWgt, impDens) ;
  324 
  325   } else if (indexCat) {
  326 
  327 
  328     if (impSliceHistos.GetSize()>0) {
  329 
  330       // Initialize importing mapped set of TH1s
  331       map<string,TH1*> hmap ;
  332       char tmp[1024] ;
  333       strlcpy(tmp,impSliceNames,1024) ;
  334       char* token = strtok(tmp,",") ;
  335       TIterator* hiter = impSliceHistos.MakeIterator() ;
  336       while(token) {
  337    hmap[token] = (TH1*) hiter->Next() ;
  338    token = strtok(0,",") ;
  339       }
  340       importTH1Set(vars,*indexCat,hmap,initWgt,kFALSE) ;
  341     } else {
  342 
  343       // Initialize importing mapped set of RooDataHists
  344       map<string,RooDataHist*> dmap ;
  345       char tmp[1024] ;
  346       strlcpy(tmp,impSliceDNames,1024) ;
  347       char* token = strtok(tmp,",") ;
  348       TIterator* hiter = impSliceDHistos.MakeIterator() ;
  349       while(token) {
  350    dmap[token] = (RooDataHist*) hiter->Next() ;
  351    token = strtok(0,",") ;
  352       }
  353       importDHistSet(vars,*indexCat,dmap,initWgt) ;
  354 
  355     }
  356 
  357 
  358   } else {
  359 
  360     // Initialize empty
  361     initialize() ;
  362     appendToDir(this,kTRUE) ;
  363 
  364   }
  365 
  366   _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  367   TRACE_CREATE
  368 
  369 }
  370 
  371 
  372 
  373 
  374 ////////////////////////////////////////////////////////////////////////////////
  375 /// Import data from given TH1/2/3 into this RooDataHist
  376 
  377 void RooDataHist::importTH1(const RooArgList& vars, const TH1& histo, Double_t wgt, Bool_t doDensityCorrection)
  378 {
  379   // Adjust binning of internal observables to match that of input THx
  380   Int_t offset[3] ;
  381   adjustBinning(vars,histo,offset) ;
  382 
  383   // Initialize internal data structure
  384   initialize() ;
  385   appendToDir(this,kTRUE) ;
  386 
  387   // Define x,y,z as 1st, 2nd and 3rd observable
  388   RooRealVar* xvar = (RooRealVar*) _vars.find(vars.at(0)->GetName()) ;
  389   RooRealVar* yvar = (RooRealVar*) (vars.at(1) ? _vars.find(vars.at(1)->GetName()) : 0 ) ;
  390   RooRealVar* zvar = (RooRealVar*) (vars.at(2) ? _vars.find(vars.at(2)->GetName()) : 0 ) ;
  391 
  392   // Transfer contents
  393   Int_t xmin(0),ymin(0),zmin(0) ;
  394   RooArgSet vset(*xvar) ;
  395   Double_t volume = xvar->getMax()-xvar->getMin() ;
  396   xmin = offset[0] ;
  397   if (yvar) {
  398     vset.add(*yvar) ;
  399     ymin = offset[1] ;
  400     volume *= (yvar->getMax()-yvar->getMin()) ;
  401   }
  402   if (zvar) {
  403     vset.add(*zvar) ;
  404     zmin = offset[2] ;
  405     volume *= (zvar->getMax()-zvar->getMin()) ;
  406   }
  407   //Double_t avgBV = volume / numEntries() ;
  408 //   cout << "average bin volume = " << avgBV << endl ;
  409 
  410   Int_t ix(0),iy(0),iz(0) ;
  411   for (ix=0 ; ix < xvar->getBins() ; ix++) {
  412     xvar->setBin(ix) ;
  413     if (yvar) {
  414       for (iy=0 ; iy < yvar->getBins() ; iy++) {
  415    yvar->setBin(iy) ;
  416    if (zvar) {
  417      for (iz=0 ; iz < zvar->getBins() ; iz++) {
  418        zvar->setBin(iz) ;
  419        Double_t bv = doDensityCorrection ? binVolume(vset) : 1;
  420        add(vset,bv*histo.GetBinContent(ix+1+xmin,iy+1+ymin,iz+1+zmin)*wgt,bv*TMath::Power(histo.GetBinError(ix+1+xmin,iy+1+ymin,iz+1+zmin)*wgt,2)) ;
  421      }
  422    } else {
  423      Double_t bv = doDensityCorrection ? binVolume(vset) : 1;
  424      add(vset,bv*histo.GetBinContent(ix+1+xmin,iy+1+ymin)*wgt,bv*TMath::Power(histo.GetBinError(ix+1+xmin,iy+1+ymin)*wgt,2)) ;
  425    }
  426       }
  427     } else {
  428       Double_t bv = doDensityCorrection ? binVolume(vset) : 1 ;
  429       add(vset,bv*histo.GetBinContent(ix+1+xmin)*wgt,bv*TMath::Power(histo.GetBinError(ix+1+xmin)*wgt,2)) ;
  430     }
  431   }
  432 
  433 }
  434 
  435 
  436 
  437 
  438 
  439 ////////////////////////////////////////////////////////////////////////////////
  440 /// Import data from given set of TH1/2/3 into this RooDataHist. The category indexCat labels the sources
  441 /// in the constructed RooDataHist. The stl map provides the mapping between the indexCat state labels
  442 /// and the import source
  443 
  444 void RooDataHist::importTH1Set(const RooArgList& vars, RooCategory& indexCat, map<string,TH1*> hmap, Double_t wgt, Bool_t doDensityCorrection)
  445 {
  446   RooCategory* icat = (RooCategory*) _vars.find(indexCat.GetName()) ;
  447 
  448   TH1* histo(0) ;
  449   Bool_t init(kFALSE) ;
  450   for (map<string,TH1*>::iterator hiter = hmap.begin() ; hiter!=hmap.end() ; ++hiter) {
  451     // Store pointer to first histogram from which binning specification will be taken
  452     if (!histo) {
  453       histo = hiter->second ;
  454     }
  455     // Define state labels in index category (both in provided indexCat and in internal copy in dataset)
  456     if (!indexCat.lookupType(hiter->first.c_str())) {
  457       indexCat.defineType(hiter->first.c_str()) ;
  458       coutI(InputArguments) << "RooDataHist::importTH1Set(" << GetName() << ") defining state \"" << hiter->first << "\" in index category " << indexCat.GetName() << endl ;
  459     }
  460     if (!icat->lookupType(hiter->first.c_str())) {
  461       icat->defineType(hiter->first.c_str()) ;
  462     }
  463   }
  464 
  465   // Check consistency in number of dimensions
  466   if (histo && (vars.getSize() != histo->GetDimension())) {
  467     coutE(InputArguments) << "RooDataHist::ctor(" << GetName() << ") ERROR: dimension of input histogram must match "
  468            << "number of continuous variables" << endl ;
  469     assert(0) ;
  470   }
  471 
  472   // Copy bins and ranges from THx to dimension observables
  473   Int_t offset[3] ;
  474   adjustBinning(vars,*histo,offset) ;
  475 
  476   // Initialize internal data structure
  477   if (!init) {
  478     initialize() ;
  479     appendToDir(this,kTRUE) ;
  480     init = kTRUE ;
  481   }
  482 
  483   // Define x,y,z as 1st, 2nd and 3rd observable
  484   RooRealVar* xvar = (RooRealVar*) _vars.find(vars.at(0)->GetName()) ;
  485   RooRealVar* yvar = (RooRealVar*) (vars.at(1) ? _vars.find(vars.at(1)->GetName()) : 0 ) ;
  486   RooRealVar* zvar = (RooRealVar*) (vars.at(2) ? _vars.find(vars.at(2)->GetName()) : 0 ) ;
  487 
  488   // Transfer contents
  489   Int_t xmin(0),ymin(0),zmin(0) ;
  490   RooArgSet vset(*xvar) ;
  491   Double_t volume = xvar->getMax()-xvar->getMin() ;
  492   xmin = offset[0] ;
  493   if (yvar) {
  494     vset.add(*yvar) ;
  495     ymin = offset[1] ;
  496     volume *= (yvar->getMax()-yvar->getMin()) ;
  497   }
  498   if (zvar) {
  499     vset.add(*zvar) ;
  500     zmin = offset[2] ;
  501     volume *= (zvar->getMax()-zvar->getMin()) ;
  502   }
  503   Double_t avgBV = volume / numEntries() ;
  504 
  505   Int_t ic(0),ix(0),iy(0),iz(0) ;
  506   for (ic=0 ; ic < icat->numBins(0) ; ic++) {
  507     icat->setBin(ic) ;
  508     histo = hmap[icat->getLabel()] ;
  509     for (ix=0 ; ix < xvar->getBins() ; ix++) {
  510       xvar->setBin(ix) ;
  511       if (yvar) {
  512    for (iy=0 ; iy < yvar->getBins() ; iy++) {
  513      yvar->setBin(iy) ;
  514      if (zvar) {
  515        for (iz=0 ; iz < zvar->getBins() ; iz++) {
  516          zvar->setBin(iz) ;
  517          Double_t bv = doDensityCorrection ? binVolume(vset)/avgBV : 1;
  518          add(vset,bv*histo->GetBinContent(ix+1+xmin,iy+1+ymin,iz+1+zmin)*wgt,bv*TMath::Power(histo->GetBinError(ix+1+xmin,iy+1+ymin,iz+1+zmin)*wgt,2)) ;
  519        }
  520      } else {
  521        Double_t bv = doDensityCorrection ? binVolume(vset)/avgBV : 1;
  522        add(vset,bv*histo->GetBinContent(ix+1+xmin,iy+1+ymin)*wgt,bv*TMath::Power(histo->GetBinError(ix+1+xmin,iy+1+ymin)*wgt,2)) ;
  523      }
  524    }
  525       } else {
  526    Double_t bv = doDensityCorrection ? binVolume(vset)/avgBV : 1;
  527    add(vset,bv*histo->GetBinContent(ix+1+xmin)*wgt,bv*TMath::Power(histo->GetBinError(ix+1+xmin)*wgt,2)) ;
  528       }
  529     }
  530   }
  531 
  532 }
  533 
  534 
  535 
  536 ////////////////////////////////////////////////////////////////////////////////
  537 /// Import data from given set of TH1/2/3 into this RooDataHist. The category indexCat labels the sources
  538 /// in the constructed RooDataHist. The stl map provides the mapping between the indexCat state labels
  539 /// and the import source
  540 
  541 void RooDataHist::importDHistSet(const RooArgList& /*vars*/, RooCategory& indexCat, std::map<std::string,RooDataHist*> dmap, Double_t initWgt)
  542 {
  543   RooCategory* icat = (RooCategory*) _vars.find(indexCat.GetName()) ;
  544 
  545   for (map<string,RooDataHist*>::iterator diter = dmap.begin() ; diter!=dmap.end() ; ++diter) {
  546 
  547     // Define state labels in index category (both in provided indexCat and in internal copy in dataset)
  548     if (!indexCat.lookupType(diter->first.c_str())) {
  549       indexCat.defineType(diter->first.c_str()) ;
  550       coutI(InputArguments) << "RooDataHist::importDHistSet(" << GetName() << ") defining state \"" << diter->first << "\" in index category " << indexCat.GetName() << endl ;
  551     }
  552     if (!icat->lookupType(diter->first.c_str())) {
  553       icat->defineType(diter->first.c_str()) ;
  554     }
  555   }
  556 
  557   initialize() ;
  558   appendToDir(this,kTRUE) ;
  559 
  560 
  561   for (map<string,RooDataHist*>::iterator diter = dmap.begin() ; diter!=dmap.end() ; ++diter) {
  562 
  563     RooDataHist* dhist = diter->second ;
  564 
  565     icat->setLabel(diter->first.c_str()) ;
  566 
  567     // Transfer contents
  568     for (Int_t i=0 ; i<dhist->numEntries() ; i++) {
  569       _vars = *dhist->get(i) ;
  570       add(_vars,dhist->weight()*initWgt, pow(dhist->weightError(SumW2),2) ) ;
  571     }
  572 
  573   }
  574 }
  575 
  576 ////////////////////////////////////////////////////////////////////////////////
  577 /// Adjust binning specification on first and optionally second and third
  578 /// observable to binning in given reference TH1. Used by constructors
  579 /// that import data from an external TH1
  580 
  581 void RooDataHist::adjustBinning(const RooArgList& vars, const TH1& href, Int_t* offset)
  582 {
  583   // X
  584   RooRealVar* xvar = (RooRealVar*) _vars.find(*vars.at(0)) ;
  585   if (!dynamic_cast<RooRealVar*>(xvar)) {
  586     coutE(InputArguments) << "RooDataHist::adjustBinning(" << GetName() << ") ERROR: dimension " << xvar->GetName() << " must be real" << endl ;
  587     assert(0) ;
  588   }
  589 
  590   Double_t xlo = ((RooRealVar*)vars.at(0))->getMin() ;
  591   Double_t xhi = ((RooRealVar*)vars.at(0))->getMax() ;
  592   Int_t xmin(0) ;
  593   if (href.GetXaxis()->GetXbins()->GetArray()) {
  594 
  595     RooBinning xbins(href.GetNbinsX(),href.GetXaxis()->GetXbins()->GetArray()) ;
  596 
  597     Double_t tolerance = 1e-6*xbins.averageBinWidth() ;
  598 
  599     // Adjust xlo/xhi to nearest boundary
  600     Double_t xloAdj = xbins.binLow(xbins.binNumber(xlo+tolerance)) ;
  601     Double_t xhiAdj = xbins.binHigh(xbins.binNumber(xhi-tolerance)) ;
  602     xbins.setRange(xloAdj,xhiAdj) ;
  603 
  604     ((RooRealVar*)vars.at(0))->setBinning(xbins) ;
  605     if (fabs(xloAdj-xlo)>tolerance||fabs(xhiAdj-xhi)<tolerance) {
  606       coutI(DataHandling) << "RooDataHist::adjustBinning(" << GetName() << "): fit range of variable " << xvar->GetName() << " expanded to nearest bin boundaries: ["
  607            << xlo << "," << xhi << "] --> [" << xloAdj << "," << xhiAdj << "]" << endl ;
  608     }
  609 
  610     xvar->setBinning(xbins) ;
  611     xmin = xbins.rawBinNumber(xloAdj+tolerance) ;
  612     if (offset) {
  613       offset[0] = xmin ;
  614     }
  615 
  616   } else {
  617 
  618     RooBinning xbins(href.GetXaxis()->GetXmin(),href.GetXaxis()->GetXmax()) ;
  619     xbins.addUniform(href.GetNbinsX(),href.GetXaxis()->GetXmin(),href.GetXaxis()->GetXmax()) ;
  620 
  621     Double_t tolerance = 1e-6*xbins.averageBinWidth() ;
  622 
  623     // Adjust xlo/xhi to nearest boundary
  624     Double_t xloAdj = xbins.binLow(xbins.binNumber(xlo+tolerance)) ;
  625     Double_t xhiAdj = xbins.binHigh(xbins.binNumber(xhi-tolerance)) ;
  626     xbins.setRange(xloAdj,xhiAdj) ;
  627     ((RooRealVar*)vars.at(0))->setRange(xloAdj,xhiAdj) ;
  628     if (fabs(xloAdj-xlo)>tolerance||fabs(xhiAdj-xhi)<tolerance) {
  629       coutI(DataHandling) << "RooDataHist::adjustBinning(" << GetName() << "): fit range of variable " << xvar->GetName() << " expanded to nearest bin boundaries: ["
  630            << xlo << "," << xhi << "] --> [" << xloAdj << "," << xhiAdj << "]" << endl ;
  631     }
  632 
  633 
  634     RooUniformBinning xbins2(xloAdj,xhiAdj,xbins.numBins()) ;
  635     xvar->setBinning(xbins2) ;
  636     xmin = xbins.rawBinNumber(xloAdj+tolerance) ;
  637     if (offset) {
  638       offset[0] = xmin ;
  639     }
  640   }
  641 
  642 
  643 
  644   // Y
  645   RooRealVar* yvar = (RooRealVar*) (vars.at(1) ? _vars.find(*vars.at(1)) : 0 ) ;
  646   Int_t ymin(0) ;
  647   if (yvar) {
  648     Double_t ylo = ((RooRealVar*)vars.at(1))->getMin() ;
  649     Double_t yhi = ((RooRealVar*)vars.at(1))->getMax() ;
  650 
  651     if (!dynamic_cast<RooRealVar*>(yvar)) {
  652       coutE(InputArguments) << "RooDataHist::adjustBinning(" << GetName() << ") ERROR: dimension " << yvar->GetName() << " must be real" << endl ;
  653       assert(0) ;
  654     }
  655 
  656     if (href.GetYaxis()->GetXbins()->GetArray()) {
  657 
  658       RooBinning ybins(href.GetNbinsY(),href.GetYaxis()->GetXbins()->GetArray()) ;
  659 
  660       Double_t tolerance = 1e-6*ybins.averageBinWidth() ;
  661 
  662       // Adjust ylo/yhi to nearest boundary
  663       Double_t yloAdj = ybins.binLow(ybins.binNumber(ylo+tolerance)) ;
  664       Double_t yhiAdj = ybins.binHigh(ybins.binNumber(yhi-tolerance)) ;
  665       ybins.setRange(yloAdj,yhiAdj) ;
  666       ((RooRealVar*)vars.at(1))->setBinning(ybins) ;
  667       if (fabs(yloAdj-ylo)>tolerance||fabs(yhiAdj-yhi)<tolerance) {
  668    coutI(DataHandling) << "RooDataHist::adjustBinning(" << GetName() << "): fit range of variable " << yvar->GetName() << " expanded to nearest bin boundaries: ["
  669              << ylo << "," << yhi << "] --> [" << yloAdj << "," << yhiAdj << "]" << endl ;
  670       }
  671 
  672       yvar->setBinning(ybins) ;
  673       ymin = ybins.rawBinNumber(yloAdj+tolerance) ;
  674       if (offset) {
  675    offset[1] = ymin ;
  676       }
  677 
  678     } else {
  679 
  680       RooBinning ybins(href.GetYaxis()->GetXmin(),href.GetYaxis()->GetXmax()) ;
  681       ybins.addUniform(href.GetNbinsY(),href.GetYaxis()->GetXmin(),href.GetYaxis()->GetXmax()) ;
  682 
  683       Double_t tolerance = 1e-6*ybins.averageBinWidth() ;
  684 
  685       // Adjust ylo/yhi to nearest boundary
  686       Double_t yloAdj = ybins.binLow(ybins.binNumber(ylo+tolerance)) ;
  687       Double_t yhiAdj = ybins.binHigh(ybins.binNumber(yhi-tolerance)) ;
  688       ybins.setRange(yloAdj,yhiAdj) ;
  689       ((RooRealVar*)vars.at(1))->setRange(yloAdj,yhiAdj) ;
  690       if (fabs(yloAdj-ylo)>tolerance||fabs(yhiAdj-yhi)<tolerance) {
  691    coutI(DataHandling) << "RooDataHist::adjustBinning(" << GetName() << "): fit range of variable " << yvar->GetName() << " expanded to nearest bin boundaries: ["
  692              << ylo << "," << yhi << "] --> [" << yloAdj << "," << yhiAdj << "]" << endl ;
  693       }
  694 
  695       RooUniformBinning ybins2(yloAdj,yhiAdj,ybins.numBins()) ;
  696       yvar->setBinning(ybins2) ;
  697       ymin = ybins.rawBinNumber(yloAdj+tolerance) ;
  698       if (offset) {
  699    offset[1] = ymin ;
  700       }
  701 
  702     }
  703   }
  704 
  705   // Z
  706   RooRealVar* zvar = (RooRealVar*) (vars.at(2) ? _vars.find(*vars.at(2)) : 0 ) ;
  707   Int_t zmin(0) ;
  708   if (zvar) {
  709     Double_t zlo = ((RooRealVar*)vars.at(2))->getMin() ;
  710     Double_t zhi = ((RooRealVar*)vars.at(2))->getMax() ;
  711 
  712     if (!dynamic_cast<RooRealVar*>(zvar)) {
  713       coutE(InputArguments) << "RooDataHist::adjustBinning(" << GetName() << ") ERROR: dimension " << zvar->GetName() << " must be real" << endl ;
  714       assert(0) ;
  715     }
  716 
  717     if (href.GetZaxis()->GetXbins()->GetArray()) {
  718 
  719       RooBinning zbins(href.GetNbinsZ(),href.GetZaxis()->GetXbins()->GetArray()) ;
  720 
  721       Double_t tolerance = 1e-6*zbins.averageBinWidth() ;
  722 
  723       // Adjust zlo/zhi to nearest boundary
  724       Double_t zloAdj = zbins.binLow(zbins.binNumber(zlo+tolerance)) ;
  725       Double_t zhiAdj = zbins.binHigh(zbins.binNumber(zhi-tolerance)) ;
  726       zbins.setRange(zloAdj,zhiAdj) ;
  727       ((RooRealVar*)vars.at(2))->setBinning(zbins) ;
  728       if (fabs(zloAdj-zlo)>tolerance||fabs(zhiAdj-zhi)<tolerance) {
  729    coutI(DataHandling) << "RooDataHist::adjustBinning(" << GetName() << "): fit range of variable " << zvar->GetName() << " expanded to nearest bin boundaries: ["
  730              << zlo << "," << zhi << "] --> [" << zloAdj << "," << zhiAdj << "]" << endl ;
  731       }
  732 
  733       zvar->setBinning(zbins) ;
  734       zmin = zbins.rawBinNumber(zloAdj+tolerance) ;
  735       if (offset) {
  736    offset[2] = zmin ;
  737       }
  738 
  739     } else {
  740 
  741       RooBinning zbins(href.GetZaxis()->GetXmin(),href.GetZaxis()->GetXmax()) ;
  742       zbins.addUniform(href.GetNbinsZ(),href.GetZaxis()->GetXmin(),href.GetZaxis()->GetXmax()) ;
  743 
  744       Double_t tolerance = 1e-6*zbins.averageBinWidth() ;
  745 
  746       // Adjust zlo/zhi to nearest boundary
  747       Double_t zloAdj = zbins.binLow(zbins.binNumber(zlo+tolerance)) ;
  748       Double_t zhiAdj = zbins.binHigh(zbins.binNumber(zhi-tolerance)) ;
  749       zbins.setRange(zloAdj,zhiAdj) ;
  750       ((RooRealVar*)vars.at(2))->setRange(zloAdj,zhiAdj) ;
  751       if (fabs(zloAdj-zlo)>tolerance||fabs(zhiAdj-zhi)<tolerance) {
  752    coutI(DataHandling) << "RooDataHist::adjustBinning(" << GetName() << "): fit range of variable " << zvar->GetName() << " expanded to nearest bin boundaries: ["
  753              << zlo << "," << zhi << "] --> [" << zloAdj << "," << zhiAdj << "]" << endl ;
  754       }
  755 
  756       RooUniformBinning zbins2(zloAdj,zhiAdj,zbins.numBins()) ;
  757       zvar->setBinning(zbins2) ;
  758       zmin = zbins.rawBinNumber(zloAdj+tolerance) ;
  759       if (offset) {
  760    offset[2] = zmin ;
  761       }
  762     }
  763   }
  764 
  765 }
  766 
  767 
  768 
  769 
  770 
  771 ////////////////////////////////////////////////////////////////////////////////
  772 /// Initialization procedure: allocate weights array, calculate
  773 /// multipliers needed for N-space to 1-dim array jump table,
  774 /// and fill the internal tree with all bin center coordinates
  775 
  776 void RooDataHist::initialize(const char* binningName, Bool_t fillTree)
  777 {
  778 
  779   // Save real dimensions of dataset separately
  780   RooAbsArg* real ;
  781   _iterator->Reset() ;
  782   while((real=(RooAbsArg*)_iterator->Next())) {
  783     if (dynamic_cast<RooAbsReal*>(real)) _realVars.add(*real);
  784   }
  785   _realIter = _realVars.createIterator() ;
  786 
  787   // Fill array of LValue pointers to variables
  788   _iterator->Reset();
  789   RooAbsArg* rvarg;
  790   while((rvarg=(RooAbsArg*)_iterator->Next())) {
  791     if (binningName) {
  792       RooRealVar* rrv = dynamic_cast<RooRealVar*>(rvarg);
  793       if (rrv) {
  794    rrv->setBinning(rrv->getBinning(binningName));
  795       }
  796     }
  797     // coverity[FORWARD_NULL]
  798     _lvvars.push_back(dynamic_cast<RooAbsLValue*>(rvarg));
  799     // coverity[FORWARD_NULL]
  800     const RooAbsBinning* binning = dynamic_cast<RooAbsLValue*>(rvarg)->getBinningPtr(0);
  801     _lvbins.push_back(binning ? binning->clone() : 0);
  802   }
  803 
  804 
  805   // Allocate coefficients array
  806   _idxMult.resize(_vars.getSize()) ;
  807 
  808   _arrSize = 1 ;
  809   _iterator->Reset() ;
  810   RooAbsLValue* arg ;
  811   Int_t n(0), i ;
  812   while((arg=dynamic_cast<RooAbsLValue*>(_iterator->Next()))) {
  813 
  814     // Calculate sub-index multipliers for master index
  815     for (i=0 ; i<n ; i++) {
  816       _idxMult[i] *= arg->numBins() ;
  817     }
  818     _idxMult[n++] = 1 ;
  819 
  820     // Calculate dimension of weight array
  821     _arrSize *= arg->numBins() ;
  822   }
  823 
  824   // Allocate and initialize weight array if necessary
  825   if (!_wgt) {
  826     _wgt = new Double_t[_arrSize] ;
  827     _errLo = new Double_t[_arrSize] ;
  828     _errHi = new Double_t[_arrSize] ;
  829     _sumw2 = new Double_t[_arrSize] ;
  830     _binv = new Double_t[_arrSize] ;
  831 
  832     // Refill array pointers in data store when reading
  833     // from Streamer
  834     if (fillTree==kFALSE) {
  835       _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  836     }
  837 
  838     for (i=0 ; i<_arrSize ; i++) {
  839       _wgt[i] = 0 ;
  840       _errLo[i] = -1 ;
  841       _errHi[i] = -1 ;
  842       _sumw2[i] = 0 ;
  843     }
  844   }
  845 
  846   if (!fillTree) return ;
  847 
  848   // Fill TTree with bin center coordinates
  849   // Calculate plot bins of components from master index
  850 
  851   Int_t ibin ;
  852   for (ibin=0 ; ibin<_arrSize ; ibin++) {
  853     _iterator->Reset() ;
  854     RooAbsLValue* arg2 ;
  855     Int_t j(0), idx(0), tmp(ibin) ;
  856     Double_t theBinVolume(1) ;
  857     while((arg2=dynamic_cast<RooAbsLValue*>(_iterator->Next()))) {
  858       idx  = tmp / _idxMult[j] ;
  859       tmp -= idx*_idxMult[j++] ;
  860       RooAbsLValue* arglv = dynamic_cast<RooAbsLValue*>(arg2) ;
  861       arglv->setBin(idx) ;
  862       theBinVolume *= arglv->getBinWidth(idx) ;
  863 //       cout << "init: bin width at idx=" << idx << " = " << arglv->getBinWidth(idx) << " binv[" << idx << "] = " << theBinVolume << endl ;
  864     }
  865     _binv[ibin] = theBinVolume ;
  866 //     cout << "binv[" << ibin << "] = " << theBinVolume << endl ;
  867     fill() ;
  868   }
  869 
  870 
  871 }
  872 
  873 
  874 ////////////////////////////////////////////////////////////////////////////////
  875 
  876 void RooDataHist::checkBinBounds() const
  877 {
  878   if (!_binbounds.empty()) return;
  879   for (std::vector<const RooAbsBinning*>::const_iterator it = _lvbins.begin();
  880       _lvbins.end() != it; ++it) {
  881     _binbounds.push_back(std::vector<Double_t>());
  882     if (*it) {
  883       std::vector<Double_t>& bounds = _binbounds.back();
  884       bounds.reserve(2 * (*it)->numBins());
  885       for (Int_t i = 0; i < (*it)->numBins(); ++i) {
  886    bounds.push_back((*it)->binLow(i));
  887    bounds.push_back((*it)->binHigh(i));
  888       }
  889     }
  890   }
  891 }
  892 
  893 ////////////////////////////////////////////////////////////////////////////////
  894 /// Copy constructor
  895 
  896 RooDataHist::RooDataHist(const RooDataHist& other, const char* newname) :
  897   RooAbsData(other,newname), RooDirItem(), _idxMult(other._idxMult), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(other._pbinvCacheMgr,0), _cache_sum_valid(0)
  898 {
  899   Int_t i ;
  900 
  901   // Allocate and initialize weight array
  902   _arrSize = other._arrSize ;
  903   _wgt = new Double_t[_arrSize] ;
  904   _errLo = new Double_t[_arrSize] ;
  905   _errHi = new Double_t[_arrSize] ;
  906   _binv = new Double_t[_arrSize] ;
  907   _sumw2 = new Double_t[_arrSize] ;
  908   for (i=0 ; i<_arrSize ; i++) {
  909     _wgt[i] = other._wgt[i] ;
  910     _errLo[i] = other._errLo[i] ;
  911     _errHi[i] = other._errHi[i] ;
  912     _sumw2[i] = other._sumw2[i] ;
  913     _binv[i] = other._binv[i] ;
  914   }
  915 
  916   // Save real dimensions of dataset separately
  917   RooAbsArg* arg ;
  918   _iterator->Reset() ;
  919   while((arg=(RooAbsArg*)_iterator->Next())) {
  920     if (dynamic_cast<RooAbsReal*>(arg)) _realVars.add(*arg) ;
  921   }
  922   _realIter = _realVars.createIterator() ;
  923 
  924   // Fill array of LValue pointers to variables
  925   _iterator->Reset() ;
  926   RooAbsArg* rvarg ;
  927   while((rvarg=(RooAbsArg*)_iterator->Next())) {
  928     // coverity[FORWARD_NULL]
  929     _lvvars.push_back(dynamic_cast<RooAbsLValue*>(rvarg)) ;
  930     // coverity[FORWARD_NULL]
  931     const RooAbsBinning* binning = dynamic_cast<RooAbsLValue*>(rvarg)->getBinningPtr(0) ;
  932     _lvbins.push_back(binning ? binning->clone() : 0) ;
  933   }
  934 
  935   _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  936 
  937  appendToDir(this,kTRUE) ;
  938 }
  939 
  940 
  941 
  942 ////////////////////////////////////////////////////////////////////////////////
  943 /// Constructor of a data hist from (part of) an existing data hist. The dimensions
  944 /// of the data set are defined by the 'vars' RooArgSet, which can be identical
  945 /// to 'dset' dimensions, or a subset thereof. Reduced dimensions will be projected
  946 /// in the output data hist. The optional 'cutVar' formula variable can used to
  947 /// select the subset of bins to be copied.
  948 ///
  949 /// For most uses the RooAbsData::reduce() wrapper function, which uses this constructor,
  950 /// is the most convenient way to create a subset of an existing data
  951 
  952 RooDataHist::RooDataHist(const char* name, const char* title, RooDataHist* h, const RooArgSet& varSubset,
  953           const RooFormulaVar* cutVar, const char* cutRange, Int_t nStart, Int_t nStop, Bool_t copyCache) :
  954   RooAbsData(name,title,varSubset),
  955   _wgt(0), _binValid(0), _curWeight(0), _curVolume(1), _pbinv(0), _pbinvCacheMgr(0,10), _cache_sum_valid(0)
  956 {
  957   // Initialize datastore
  958   _dstore = new RooTreeDataStore(name,title,*h->_dstore,_vars,cutVar,cutRange,nStart,nStop,copyCache) ;
  959 
  960   initialize(0,kFALSE) ;
  961 
  962   _dstore->setExternalWeightArray(_wgt,_errLo,_errHi,_sumw2) ;
  963 
  964   // Copy weight array etc
  965   Int_t i ;
  966   for (i=0 ; i<_arrSize ; i++) {
  967     _wgt[i] = h->_wgt[i] ;
  968     _errLo[i] = h->_errLo[i] ;
  969     _errHi[i] = h->_errHi[i] ;
  970     _sumw2[i] = h->_sumw2[i] ;
  971     _binv[i] = h->_binv[i] ;
  972   }
  973 
  974 
  975   appendToDir(this,kTRUE) ;
  976   TRACE_CREATE
  977 }
  978 
  979 
  980 ////////////////////////////////////////////////////////////////////////////////
  981 /// Construct a clone of this dataset that contains only the cached variables
  982 
  983 RooAbsData* RooDataHist::cacheClone(const RooAbsArg* newCacheOwner, const RooArgSet* newCacheVars, const char* newName)
  984 {
  985   checkInit() ;
  986 
  987   RooDataHist* dhist = new RooDataHist(newName?newName:GetName(),GetTitle(),this,*get(),0,0,0,2000000000,kTRUE) ;
  988 
  989   RooArgSet* selCacheVars = (RooArgSet*) newCacheVars->selectCommon(dhist->_cachedVars) ;
  990   dhist->attachCache(newCacheOwner, *selCacheVars) ;
  991   delete selCacheVars ;
  992 
  993   return dhist ;
  994 }
  995 
  996 
  997 
  998 ////////////////////////////////////////////////////////////////////////////////
  999 /// Implementation of RooAbsData virtual method that drives the RooAbsData::reduce() methods
 1000 
 1001 RooAbsData* RooDataHist::reduceEng(const RooArgSet& varSubset, const RooFormulaVar* cutVar, const char* cutRange,
 1002                Int_t nStart, Int_t nStop, Bool_t /*copyCache*/)
 1003 {
 1004   checkInit() ;
 1005   RooArgSet* myVarSubset = (RooArgSet*) _vars.selectCommon(varSubset) ;
 1006   RooDataHist *rdh = new RooDataHist(GetName(), GetTitle(), *myVarSubset) ;
 1007   delete myVarSubset ;
 1008 
 1009   RooFormulaVar* cloneVar = 0;
 1010   RooArgSet* tmp(0) ;
 1011   if (cutVar) {
 1012     // Deep clone cutVar and attach clone to this dataset
 1013     tmp = (RooArgSet*) RooArgSet(*cutVar).snapshot() ;
 1014     if (!tmp) {
 1015       coutE(DataHandling) << "RooDataHist::reduceEng(" << GetName() << ") Couldn't deep-clone cut variable, abort," << endl ;
 1016       return 0 ;
 1017     }
 1018     cloneVar = (RooFormulaVar*) tmp->find(*cutVar) ;
 1019     cloneVar->attachDataSet(*this) ;
 1020   }
 1021 
 1022   Int_t i ;
 1023   Double_t lo,hi ;
 1024   Int_t nevt = nStop < numEntries() ? nStop : numEntries() ;
 1025   TIterator* vIter = get()->createIterator() ;
 1026   for (i=nStart ; i<nevt ; i++) {
 1027     const RooArgSet* row = get(i) ;
 1028 
 1029     Bool_t doSelect(kTRUE) ;
 1030     if (cutRange) {
 1031       RooAbsArg* arg ;
 1032       vIter->Reset() ;
 1033       while((arg=(RooAbsArg*)vIter->Next())) {
 1034    if (!arg->inRange(cutRange)) {
 1035      doSelect = kFALSE ;
 1036      break ;
 1037    }
 1038       }
 1039     }
 1040     if (!doSelect) continue ;
 1041 
 1042     if (!cloneVar || cloneVar->getVal()) {
 1043       weightError(lo,hi,SumW2) ;
 1044       rdh->add(*row,weight(),lo*lo) ;
 1045     }
 1046   }
 1047   delete vIter ;
 1048 
 1049   if (cloneVar) {
 1050     delete tmp ;
 1051   }
 1052 
 1053     return rdh ;
 1054   }
 1055 
 1056 
 1057 
 1058 ////////////////////////////////////////////////////////////////////////////////
 1059 /// Destructor
 1060 
 1061 RooDataHist::~RooDataHist()
 1062 {
 1063   if (_wgt) delete[] _wgt ;
 1064   if (_errLo) delete[] _errLo ;
 1065   if (_errHi) delete[] _errHi ;
 1066   if (_sumw2) delete[] _sumw2 ;
 1067   if (_binv) delete[] _binv ;
 1068   if (_realIter) delete _realIter ;
 1069   if (_binValid) delete[] _binValid ;
 1070   vector<const RooAbsBinning*>::iterator iter = _lvbins.begin() ;
 1071   while(iter!=_lvbins.end()) {
 1072     delete *iter ;
 1073     iter++ ;
 1074   }
 1075 
 1076    removeFromDir(this) ;
 1077   TRACE_DESTROY
 1078 }
 1079 
 1080 
 1081 
 1082 
 1083 ////////////////////////////////////////////////////////////////////////////////
 1084 
 1085 Int_t RooDataHist::getIndex(const RooArgSet& coord, Bool_t fast)
 1086 {
 1087   checkInit() ;
 1088   if (fast) {
 1089     _vars.assignFast(coord,kFALSE) ;
 1090   } else {
 1091     _vars.assignValueOnly(coord) ;
 1092   }
 1093   return calcTreeIndex() ;
 1094 }
 1095 
 1096 
 1097 
 1098 
 1099 ////////////////////////////////////////////////////////////////////////////////
 1100 /// Calculate the index for the weights array corresponding to
 1101 /// to the bin enclosing the current coordinates of the internal argset
 1102 
 1103 Int_t RooDataHist::calcTreeIndex() const
 1104 {
 1105   Int_t masterIdx(0), i(0) ;
 1106   vector<RooAbsLValue*>::const_iterator iter = _lvvars.begin() ;
 1107   vector<const RooAbsBinning*>::const_iterator biter = _lvbins.begin() ;
 1108   for (;iter!=_lvvars.end() ; ++iter) {
 1109     const RooAbsBinning* binning = (*biter) ;
 1110     masterIdx += _idxMult[i++]*(*iter)->getBin(binning) ;
 1111     biter++ ;
 1112   }
 1113   return masterIdx ;
 1114 }
 1115 
 1116 
 1117 
 1118 ////////////////////////////////////////////////////////////////////////////////
 1119 /// Debug stuff, should go...
 1120 
 1121 void RooDataHist::dump2()
 1122 {
 1123   cout << "_arrSize = " << _arrSize << endl ;
 1124   for (Int_t i=0 ; i<_arrSize ; i++) {
 1125     cout << "wgt[" << i << "] = " << _wgt[i] << "sumw2[" << i << "] = " << _sumw2[i] << " vol[" << i << "] = " << _binv[i] << endl ;
 1126   }
 1127 }
 1128 
 1129 
 1130 
 1131 ////////////////////////////////////////////////////////////////////////////////
 1132 /// Back end function to plotting functionality. Plot RooDataHist on given
 1133 /// frame in mode specified by plot options 'o'. The main purpose of
 1134 /// this function is to match the specified binning on 'o' to the
 1135 /// internal binning of the plot observable in this RooDataHist.
 1136 
 1137 RooPlot *RooDataHist::plotOn(RooPlot *frame, PlotOpt o) const
 1138 {
 1139   checkInit() ;
 1140   if (o.bins) return RooAbsData::plotOn(frame,o) ;
 1141 
 1142   if(0 == frame) {
 1143     coutE(InputArguments) << ClassName() << "::" << GetName() << ":plotOn: frame is null" << endl;
 1144     return 0;
 1145   }
 1146   RooAbsRealLValue *var= (RooAbsRealLValue*) frame->getPlotVar();
 1147   if(0 == var) {
 1148     coutE(InputArguments) << ClassName() << "::" << GetName()
 1149     << ":plotOn: frame does not specify a plot variable" << endl;
 1150     return 0;
 1151   }
 1152 
 1153   RooRealVar* dataVar = (RooRealVar*) _vars.find(*var) ;
 1154   if (!dataVar) {
 1155     coutE(InputArguments) << ClassName() << "::" << GetName()
 1156     << ":plotOn: dataset doesn't contain plot frame variable" << endl;
 1157     return 0;
 1158   }
 1159 
 1160   o.bins = &dataVar->getBinning() ;
 1161   o.correctForBinWidth = kFALSE ;
 1162   return RooAbsData::plotOn(frame,o) ;
 1163 }
 1164 
 1165 
 1166 
 1167 
 1168 ////////////////////////////////////////////////////////////////////////////////
 1169 
 1170 Double_t RooDataHist::weightSquared() const {
 1171   return _curSumW2 ;
 1172 }
 1173 
 1174 
 1175 
 1176 ////////////////////////////////////////////////////////////////////////////////
 1177 /// Return the weight at given coordinates with optional
 1178 /// interpolation. If intOrder is zero, the weight
 1179 /// for the bin enclosing the coordinates
 1180 /// contained in 'bin' is returned. For higher values,
 1181 /// the result is interpolated in the real dimensions
 1182 /// of the dataset
 1183 
  Double_t RooDataHist::weight(const RooArgSet& bin, Int_t intOrder, Bool_t correctForBinSize, Bool_t cdfBoundaries, Bool_t oneSafe)
  {
    //cout << "RooDataHist::weight(" << bin << "," << intOrder << "," << correctForBinSize << "," << cdfBoundaries << "," << oneSafe << ")" << endl ;
  
      checkInit() ;
   
     // Handle illegal intOrder values
     if (intOrder<0) {
       coutE(InputArguments) << "RooDataHist::weight(" << GetName() << ") ERROR: interpolation order must be positive" << endl ;
       return 0 ;
     }
   
     // Handle no-interpolation case
     if (intOrder==0) {
       _vars.assignValueOnly(bin,oneSafe) ;
       Int_t idx = calcTreeIndex() ;
       //cout << "intOrder 0, idx = " << idx << endl ;
       if (correctForBinSize) {
         //calculatePartialBinVolume(*get()) ;
         //cout << "binw[" << idx << "] = " << _wgt[idx] <<  " / " << _binv[idx] << endl ;
         return _wgt[idx] / _binv[idx] ;
       } else {
         //cout << "binw[" << idx << "] = " << _wgt[idx] << endl ;
         return _wgt[idx] ;
       }
     }
   
     // Handle all interpolation cases
     _vars.assignValueOnly(bin) ;
   
     Double_t wInt(0) ;
     if (_realVars.getSize()==1) {
 1216 
 1217     // 1-dimensional interpolation
 1218     RooFIter realIter = _realVars.fwdIterator() ;
 1219     RooRealVar* real=(RooRealVar*)realIter.next() ;
 1220     const RooAbsBinning* binning = real->getBinningPtr(0) ;
 1221     wInt = interpolateDim(*real,binning,((RooAbsReal*)bin.find(*real))->getVal(), intOrder, correctForBinSize, cdfBoundaries) ;
 1222 
     } else if (_realVars.getSize()==2) {
   
       // 2-dimensional interpolation
       RooFIter realIter = _realVars.fwdIterator() ;
       RooRealVar* realX=(RooRealVar*)realIter.next() ;
       RooRealVar* realY=(RooRealVar*)realIter.next() ;
       Double_t xval = ((RooAbsReal*)bin.find(*realX))->getVal() ;
       Double_t yval = ((RooAbsReal*)bin.find(*realY))->getVal() ;
   
       Int_t ybinC = realY->getBin() ;
       Int_t ybinLo = ybinC-intOrder/2 - ((yval<realY->getBinning().binCenter(ybinC))?1:0) ;
       Int_t ybinM = realY->numBins() ;
   
       Int_t i ;
       Double_t yarr[10] ;
       Double_t xarr[10] ;
       const RooAbsBinning* binning = realX->getBinningPtr(0) ;
       for (i=ybinLo ; i<=intOrder+ybinLo ; i++) {
         Int_t ibin ;
         if (i>=0 && i<ybinM) {
      // In range
      ibin = i ;
      realY->setBin(ibin) ;
      xarr[i-ybinLo] = realY->getVal() ;
         } else if (i>=ybinM) {
      // Overflow: mirror
      ibin = 2*ybinM-i-1 ;
      realY->setBin(ibin) ;
      xarr[i-ybinLo] = 2*realY->getMax()-realY->getVal() ;
         } else {
      // Underflow: mirror
      ibin = -i -1;
      realY->setBin(ibin) ;
      xarr[i-ybinLo] = 2*realY->getMin()-realY->getVal() ;
         }
         yarr[i-ybinLo] = interpolateDim(*realX,binning,xval,intOrder,correctForBinSize,kFALSE) ;
       }
 1260 
 1261     if (gDebug>7) {
 1262       cout << "RooDataHist interpolating data is" << endl ;
 1263       cout << "xarr = " ;
 1264       for (int q=0; q<=intOrder ; q++) cout << xarr[q] << " " ;
 1265       cout << " yarr = " ;
 1266       for (int q=0; q<=intOrder ; q++) cout << yarr[q] << " " ;
 1267       cout << endl ;
 1268     }
 1269     wInt = RooMath::interpolate(xarr,yarr,intOrder+1,yval) ;
 1270 
 1271   } else {
 1272 
 1273     // Higher dimensional scenarios not yet implemented
 1274     coutE(InputArguments) << "RooDataHist::weight(" << GetName() << ") interpolation in "
 1275     << _realVars.getSize() << " dimensions not yet implemented" << endl ;
 1276     return weight(bin,0) ;
 1277 
 1278   }
 1279 
 1280   // Cut off negative values
 1281 //   if (wInt<=0) {
 1282 //     wInt=0 ;
 1283 //   }
 1284 
 1285   //cout << "RooDataHist wInt = " << wInt << endl ;
 1286   return wInt ;
 1287 }
 1288 
 1289 
 1290 
 1291 
 1292 ////////////////////////////////////////////////////////////////////////////////
 1293 /// Return the error on current weight
 1294 
 1295 void RooDataHist::weightError(Double_t& lo, Double_t& hi, ErrorType etype) const
 1296 {
 1297   checkInit() ;
 1298 
 1299   switch (etype) {
 1300 
 1301   case Auto:
 1302     throw string(Form("RooDataHist::weightError(%s) error type Auto not allowed here",GetName())) ;
 1303     break ;
 1304 
 1305   case Expected:
 1306     throw string(Form("RooDataHist::weightError(%s) error type Expected not allowed here",GetName())) ;
 1307     break ;
 1308 
 1309   case Poisson:
 1310     if (_curWgtErrLo>=0) {
 1311       // Weight is preset or precalculated
 1312       lo = _curWgtErrLo ;
 1313       hi = _curWgtErrHi ;
 1314       return ;
 1315     }
 1316 
 1317     // Calculate poisson errors
 1318     Double_t ym,yp ;
 1319     RooHistError::instance().getPoissonInterval(Int_t(weight()+0.5),ym,yp,1) ;
 1320     _curWgtErrLo = weight()-ym ;
 1321     _curWgtErrHi = yp-weight() ;
 1322     _errLo[_curIndex] = _curWgtErrLo ;
 1323     _errHi[_curIndex] = _curWgtErrHi ;
 1324     lo = _curWgtErrLo ;
 1325     hi = _curWgtErrHi ;
 1326     return ;
 1327 
 1328   case SumW2:
 1329     lo = sqrt(_curSumW2) ;
 1330     hi = sqrt(_curSumW2) ;
 1331     return ;
 1332 
 1333   case None:
 1334     lo = 0 ;
 1335     hi = 0 ;
 1336     return ;
 1337   }
 1338 }
 1339 
 1340 
 1341 // wve adjust for variable bin sizes
 1342 
 1343 ////////////////////////////////////////////////////////////////////////////////
 1344 /// Perform boundary safe 'intOrder'-th interpolation of weights in dimension 'dim'
 1345 /// at current value 'xval'
 1346 
 1347 Double_t RooDataHist::interpolateDim(RooRealVar& dim, const RooAbsBinning* binning, Double_t xval, Int_t intOrder, Bool_t correctForBinSize, Bool_t cdfBoundaries)
 1348 {
 1349   // Fill workspace arrays spanning interpolation area
 1350   Int_t fbinC = dim.getBin(*binning) ;
 1351   Int_t fbinLo = fbinC-intOrder/2 - ((xval<binning->binCenter(fbinC))?1:0) ;
 1352   Int_t fbinM = dim.numBins(*binning) ;
 1353 
 1354 
 1355   Int_t i ;
 1356   Double_t yarr[10] ;
 1357   Double_t xarr[10] ;
 1358   for (i=fbinLo ; i<=intOrder+fbinLo ; i++) {
 1359     Int_t ibin ;
 1360     if (i>=0 && i<fbinM) {
 1361       // In range
 1362       ibin = i ;
 1363       dim.setBinFast(ibin,*binning) ;
 1364       //cout << "INRANGE: dim.getVal(ibin=" << ibin << ") = " << dim.getVal() << endl ;
 1365       xarr[i-fbinLo] = dim.getVal() ;
 1366       Int_t idx = calcTreeIndex() ;
 1367       yarr[i-fbinLo] = _wgt[idx] ;
 1368       if (correctForBinSize) yarr[i-fbinLo] /=  _binv[idx] ;
 1369     } else if (i>=fbinM) {
 1370       // Overflow: mirror
 1371       ibin = 2*fbinM-i-1 ;
 1372       dim.setBinFast(ibin,*binning) ;
 1373       //cout << "OVERFLOW: dim.getVal(ibin=" << ibin << ") = " << dim.getVal() << endl ;
 1374       if (cdfBoundaries) {
 1375    xarr[i-fbinLo] = dim.getMax()+1e-10*(i-fbinM+1) ;
 1376    yarr[i-fbinLo] = 1.0 ;
 1377       } else {
 1378    Int_t idx = calcTreeIndex() ;
 1379    xarr[i-fbinLo] = 2*dim.getMax()-dim.getVal() ;
 1380    yarr[i-fbinLo] = _wgt[idx] ;
 1381    if (correctForBinSize) yarr[i-fbinLo] /=  _binv[idx] ;
 1382       }
 1383     } else {
 1384       // Underflow: mirror
 1385       ibin = -i - 1 ;
 1386       dim.setBinFast(ibin,*binning) ;
 1387       //cout << "UNDERFLOW: dim.getVal(ibin=" << ibin << ") = " << dim.getVal() << endl ;
 1388       if (cdfBoundaries) {
 1389    xarr[i-fbinLo] = dim.getMin()-ibin*(1e-10) ; ;
 1390    yarr[i-fbinLo] = 0.0 ;
 1391       } else {
 1392    Int_t idx = calcTreeIndex() ;
 1393    xarr[i-fbinLo] = 2*dim.getMin()-dim.getVal() ;
 1394    yarr[i-fbinLo] = _wgt[idx] ;
 1395    if (correctForBinSize) yarr[i-fbinLo] /=  _binv[idx] ;
 1396       }
 1397     }
 1398     //cout << "ibin = " << ibin << endl ;
 1399   }
 1400 //   for (int k=0 ; k<=intOrder ; k++) {
 1401 //     cout << "k=" << k << " x = " << xarr[k] << " y = " << yarr[k] << endl ;
 1402 //   }
 1403   dim.setBinFast(fbinC,*binning) ;
 1404   Double_t ret = RooMath::interpolate(xarr,yarr,intOrder+1,xval) ;
 1405   return ret ;
 1406 }
 1407 
 1408 
 1409 
 1410 
 1411 ////////////////////////////////////////////////////////////////////////////////
 1412 /// Increment the weight of the bin enclosing the coordinates given
 1413 /// by 'row' by the specified amount. Add the sum of weights squared
 1414 /// for the bin by 'sumw2' rather than wgt^2
 1415 
 1416 void RooDataHist::add(const RooArgSet& row, Double_t wgt, Double_t sumw2)
 1417 {
 1418   checkInit() ;
 1419 
 1420 //   cout << "RooDataHist::add() accepted coordinate is " << endl ;
 1421 //   _vars.Print("v") ;
 1422 
 1423   _vars = row ;
 1424   Int_t idx = calcTreeIndex() ;
 1425   _wgt[idx] += wgt ;
 1426   _sumw2[idx] += (sumw2>0?sumw2:wgt*wgt) ;
 1427   _errLo[idx] = -1 ;
 1428   _errHi[idx] = -1 ;
 1429 
 1430   _cache_sum_valid = kFALSE ;
 1431 }
 1432 
 1433 
 1434 
 1435 ////////////////////////////////////////////////////////////////////////////////
 1436 /// Increment the weight of the bin enclosing the coordinates
 1437 /// given by 'row' by the specified amount. Associate errors
 1438 /// [wgtErrLo,wgtErrHi] with the event weight on this bin.
 1439 
 1440 void RooDataHist::set(const RooArgSet& row, Double_t wgt, Double_t wgtErrLo, Double_t wgtErrHi)
 1441 {
 1442   checkInit() ;
 1443 
 1444   _vars = row ;
 1445   Int_t idx = calcTreeIndex() ;
 1446   _wgt[idx] = wgt ;
 1447   _errLo[idx] = wgtErrLo ;
 1448   _errHi[idx] = wgtErrHi ;
 1449 
 1450   _cache_sum_valid = kFALSE ;
 1451 }
 1452 
 1453 
 1454 
 1455 ////////////////////////////////////////////////////////////////////////////////
 1456 /// Increment the weight of the bin enclosing the coordinates
 1457 /// given by 'row' by the specified amount. Associate errors
 1458 /// [wgtErrLo,wgtErrHi] with the event weight on this bin.
 1459 
 1460 void RooDataHist::set(Double_t wgt, Double_t wgtErr)
 1461 {
 1462   checkInit() ;
 1463 
 1464   if (_curIndex<0) {
 1465     _curIndex = calcTreeIndex() ;
 1466   }
 1467 
 1468   _wgt[_curIndex] = wgt ;
 1469   _errLo[_curIndex] = wgtErr ;
 1470   _errHi[_curIndex] = wgtErr ;
 1471   _sumw2[_curIndex] = wgtErr*wgtErr ;
 1472 
 1473   _cache_sum_valid = kFALSE ;
 1474 }
 1475 
 1476 
 1477 
 1478 ////////////////////////////////////////////////////////////////////////////////
 1479 /// Increment the weight of the bin enclosing the coordinates
 1480 /// given by 'row' by the specified amount. Associate errors
 1481 /// [wgtErrLo,wgtErrHi] with the event weight on this bin.
 1482 
 1483 void RooDataHist::set(const RooArgSet& row, Double_t wgt, Double_t wgtErr)
 1484 {
 1485   checkInit() ;
 1486 
 1487   _vars = row ;
 1488   Int_t idx = calcTreeIndex() ;
 1489   _wgt[idx] = wgt ;
 1490   _errLo[idx] = wgtErr ;
 1491   _errHi[idx] = wgtErr ;
 1492   _sumw2[idx] = wgtErr*wgtErr ;
 1493 
 1494   _cache_sum_valid = kFALSE ;
 1495 }
 1496 
 1497 
 1498 
 1499 ////////////////////////////////////////////////////////////////////////////////
 1500 /// Add all data points contained in 'dset' to this data set with given weight.
 1501 /// Optional cut string expression selects the data points to be added and can
 1502 /// reference any variable contained in this data set
 1503 
 1504 void RooDataHist::add(const RooAbsData& dset, const char* cut, Double_t wgt)
 1505 {
 1506   RooFormulaVar cutVar("select",cut,*dset.get()) ;
 1507   add(dset,&cutVar,wgt) ;
 1508 }
 1509 
 1510 
 1511 
 1512 ////////////////////////////////////////////////////////////////////////////////
 1513 /// Add all data points contained in 'dset' to this data set with given weight.
 1514 /// Optional RooFormulaVar pointer selects the data points to be added.
 1515 
 1516 void RooDataHist::add(const RooAbsData& dset, const RooFormulaVar* cutVar, Double_t wgt)
 1517 {
 1518   checkInit() ;
 1519 
 1520   RooFormulaVar* cloneVar = 0;
 1521   RooArgSet* tmp(0) ;
 1522   if (cutVar) {
 1523     // Deep clone cutVar and attach clone to this dataset
 1524     tmp = (RooArgSet*) RooArgSet(*cutVar).snapshot() ;
 1525     if (!tmp) {
 1526       coutE(DataHandling) << "RooDataHist::add(" << GetName() << ") Couldn't deep-clone cut variable, abort," << endl ;
 1527       return ;
 1528     }
 1529 
 1530     cloneVar = (RooFormulaVar*) tmp->find(*cutVar) ;
 1531     cloneVar->attachDataSet(dset) ;
 1532   }
 1533 
 1534 
 1535   Int_t i ;
 1536   for (i=0 ; i<dset.numEntries() ; i++) {
 1537     const RooArgSet* row = dset.get(i) ;
 1538     if (!cloneVar || cloneVar->getVal()) {
 1539        add(*row,wgt*dset.weight(), wgt*wgt*dset.weightSquared()) ;
 1540     }
 1541   }
 1542 
 1543   if (cloneVar) {
 1544     delete tmp ;
 1545   }
 1546 
 1547   _cache_sum_valid = kFALSE ;
 1548 }
 1549 
 1550 
 1551 
 1552 ////////////////////////////////////////////////////////////////////////////////
 1553 /// Return the sum of the weights of all hist bins.
 1554 ///
 1555 /// If correctForBinSize is specified, the sum of weights
 1556 /// is multiplied by the N-dimensional bin volume,
 1557 /// making the return value the integral over the function
 1558 /// represented by this histogram
 1559 
 1560 Double_t RooDataHist::sum(Bool_t correctForBinSize, Bool_t inverseBinCor) const
 1561 {
 1562   checkInit() ;
 1563 
 1564   // Check if result was cached
 1565   Int_t cache_code = 1 + (correctForBinSize?1:0) + ((correctForBinSize&&inverseBinCor)?1:0) ;
 1566   if (_cache_sum_valid==cache_code) {
 1567     return _cache_sum ;
 1568   }
 1569 
 1570   Int_t i ;
 1571   Double_t total(0), carry(0);
 1572   for (i=0 ; i<_arrSize ; i++) {
 1573 
 1574     Double_t theBinVolume = correctForBinSize ? (inverseBinCor ? 1/_binv[i] : _binv[i]) : 1.0 ;
 1575     // cout << "total += " << _wgt[i] << "*" << theBinVolume << endl ;
 1576     Double_t y = _wgt[i]*theBinVolume - carry;
 1577     Double_t t = total + y;
 1578     carry = (t - total) - y;
 1579     total = t;
 1580   }
 1581 
 1582   // Store result in cache
 1583   _cache_sum_valid=cache_code ;
 1584   _cache_sum = total ;
 1585 
 1586   return total ;
 1587 }
 1588 
 1589 
 1590 
 1591 ////////////////////////////////////////////////////////////////////////////////
 1592 /// Return the sum of the weights of a multi-dimensional slice of the histogram
 1593 /// by summing only over the dimensions specified in sumSet.
 1594 ///
 1595 /// The coordinates of all other dimensions are fixed to those given in sliceSet
 1596 ///
 1597 /// If correctForBinSize is specified, the sum of weights
 1598 /// is multiplied by the M-dimensional bin volume, (M = N(sumSet)),
 1599 /// making the return value the integral over the function
 1600 /// represented by this histogram
 1601 
 1602 Double_t RooDataHist::sum(const RooArgSet& sumSet, const RooArgSet& sliceSet, Bool_t correctForBinSize, Bool_t inverseBinCor)
 1603 {
 1604   checkInit() ;
 1605 
 1606   RooArgSet varSave ;
 1607   varSave.addClone(_vars) ;
 1608 
 1609   RooArgSet* sliceOnlySet = new RooArgSet(sliceSet) ;
 1610   sliceOnlySet->remove(sumSet,kTRUE,kTRUE) ;
 1611 
 1612   _vars = *sliceOnlySet ;
 1613   calculatePartialBinVolume(*sliceOnlySet) ;
 1614   delete sliceOnlySet ;
 1615 
 1616   TIterator* ssIter = sumSet.createIterator() ;
 1617 
 1618   // Calculate mask and refence plot bins for non-iterating variables
 1619   RooAbsArg* arg ;
 1620   Bool_t* mask = new Bool_t[_vars.getSize()] ;
 1621   Int_t*  refBin = new Int_t[_vars.getSize()] ;
 1622 
 1623   Int_t i(0) ;
 1624   _iterator->Reset() ;
 1625   while((arg=(RooAbsArg*)_iterator->Next())) {
 1626     if (sumSet.find(*arg)) {
 1627       mask[i] = kFALSE ;
 1628     } else {
 1629       mask[i] = kTRUE ;
 1630       // coverity[FORWARD_NULL]
 1631       refBin[i] = (dynamic_cast<RooAbsLValue*>(arg))->getBin() ;
 1632     }
 1633     i++ ;
 1634   }
 1635 
 1636   // Loop over entire data set, skipping masked entries
 1637   Double_t total(0), carry(0);
 1638   Int_t ibin ;
 1639   for (ibin=0 ; ibin<_arrSize ; ibin++) {
 1640 
 1641     Int_t idx(0), tmp(ibin), ivar(0) ;
 1642     Bool_t skip(kFALSE) ;
 1643 
 1644     // Check if this bin belongs in selected slice
 1645     _iterator->Reset() ;
 1646     // coverity[UNUSED_VALUE]
 1647     while((!skip && (arg=(RooAbsArg*)_iterator->Next()))) {
 1648       idx  = tmp / _idxMult[ivar] ;
 1649       tmp -= idx*_idxMult[ivar] ;
 1650       if (mask[ivar] && idx!=refBin[ivar]) skip=kTRUE ;
 1651       ivar++ ;
 1652     }
 1653 
 1654     if (!skip) {
 1655       Double_t theBinVolume = correctForBinSize ? (inverseBinCor ? 1/(*_pbinv)[i] : (*_pbinv)[i] ) : 1.0 ;
 1656 //       cout << "adding bin[" << ibin << "] to sum wgt = " << _wgt[ibin] << " binv = " << theBinVolume << endl ;
 1657       Double_t y = _wgt[ibin]*theBinVolume - carry;
 1658       Double_t t = total + y;
 1659       carry = (t - total) - y;
 1660       total = t;
 1661     }
 1662   }
 1663   delete ssIter ;
 1664 
 1665   delete[] mask ;
 1666   delete[] refBin ;
 1667 
 1668   _vars = varSave ;
 1669 
 1670   return total ;
 1671 }
 1672 
 1673 ////////////////////////////////////////////////////////////////////////////////
 1674 /// Return the sum of the weights of a multi-dimensional slice of the histogram
 1675 /// by summing only over the dimensions specified in sumSet.
 1676 ///
 1677 /// The coordinates of all other dimensions are fixed to those given in sliceSet
 1678 ///
 1679 /// If correctForBinSize is specified, the sum of weights
 1680 /// is multiplied by the M-dimensional bin volume, (M = N(sumSet)),
 1681 /// or the fraction of it that falls inside the range rangeName,
 1682 /// making the return value the integral over the function
 1683 /// represented by this histogram
 1684 ///
 1685 /// If correctForBinSize is not specified, the weights are multiplied by the
 1686 /// fraction of the bin volume that falls inside the range, i.e. a factor or
 1687 /// binVolumeInRange/totalBinVolume.
 1688 
 1689 Double_t RooDataHist::sum(const RooArgSet& sumSet, const RooArgSet& sliceSet,
 1690    Bool_t correctForBinSize, Bool_t inverseBinCor,
 1691    const std::map<const RooAbsArg*, std::pair<Double_t, Double_t> >& ranges)
 1692 {
 1693   checkInit();
 1694   checkBinBounds();
 1695   RooArgSet varSave;
 1696   varSave.addClone(_vars);
 1697   {
 1698     RooArgSet sliceOnlySet(sliceSet);
 1699     sliceOnlySet.remove(sumSet,kTRUE,kTRUE);
 1700     _vars = sliceOnlySet;
 1701   }
 1702 
 1703   // Calculate mask and refence plot bins for non-iterating variables,
 1704   // and get ranges for iterating variables
 1705   std::vector<bool> mask(_vars.getSize());
 1706   std::vector<Int_t> refBin(_vars.getSize());
 1707   std::vector<Double_t> rangeLo(_vars.getSize(), -std::numeric_limits<Double_t>::infinity());
 1708   std::vector<Double_t> rangeHi(_vars.getSize(), +std::numeric_limits<Double_t>::infinity());
 1709 
 1710   _iterator->Reset();
 1711   RooAbsArg* arg;
 1712   for (Int_t i = 0; (arg=(RooAbsArg*)_iterator->Next()); ++i) {
 1713     RooAbsArg* sumsetv = sumSet.find(*arg);
 1714     RooAbsArg* slicesetv = sliceSet.find(*arg);
 1715     mask[i] = !sumsetv;
 1716     if (mask[i]) {
 1717       // coverity[FORWARD_NULL]
 1718       refBin[i] = (dynamic_cast<RooAbsLValue*>(arg))->getBin();
 1719     }
 1720     std::map<const RooAbsArg*, std::pair<Double_t, Double_t> >::const_iterator
 1721    it = ranges.find(sumsetv ? sumsetv : slicesetv);
 1722     if (ranges.end() != it) {
 1723       rangeLo[i] = it->second.first;
 1724       rangeHi[i] = it->second.second;
 1725     }
 1726   }
 1727 
 1728   // Loop over entire data set, skipping masked entries
 1729   Double_t total(0), carry(0);
 1730   for (Int_t ibin = 0; ibin < _arrSize; ++ibin) {
 1731     // Check if this bin belongs in selected slice
 1732     _iterator->Reset();
 1733     Bool_t skip(kFALSE);
 1734     // coverity[UNUSED_VALUE]
 1735     for (Int_t ivar = 0, tmp = ibin;
 1736    (!skip && (arg=(RooAbsArg*)_iterator->Next())); ++ivar) {
 1737       const Int_t idx = tmp / _idxMult[ivar];
 1738       tmp -= idx*_idxMult[ivar];
 1739       if (mask[ivar] && idx!=refBin[ivar]) skip=kTRUE;
 1740     }
 1741     if (skip) continue;
 1742     _iterator->Reset();
 1743     // work out bin volume
 1744     Double_t theBinVolume = 1.;
 1745     for (Int_t ivar = 0, tmp = ibin;
 1746    (arg=(RooAbsArg*)_iterator->Next()); ++ivar) {
 1747       const Int_t idx = tmp / _idxMult[ivar];
 1748       tmp -= idx*_idxMult[ivar];
 1749       if (_binbounds[ivar].empty()) continue;
 1750       const Double_t binLo = _binbounds[ivar][2 * idx];
 1751       const Double_t binHi = _binbounds[ivar][2 * idx + 1];
 1752       if (binHi < rangeLo[ivar] || binLo > rangeHi[ivar]) {
 1753    // bin is outside of allowed range - effective bin volume is zero
 1754    theBinVolume = 0.;
 1755    break;
 1756       }
 1757       theBinVolume *=
 1758    (std::min(rangeHi[ivar], binHi) - std::max(rangeLo[ivar], binLo));
 1759     }
 1760     const Double_t corrPartial = theBinVolume / _binv[ibin];
 1761     if (0. == corrPartial) continue;
 1762     const Double_t corr = correctForBinSize ? (inverseBinCor ? 1. / _binv[ibin] : _binv[ibin] ) : 1.0;
 1763     //cout << "adding bin[" << ibin << "] to sum wgt = " << _wgt[ibin] << " binv = " << theBinVolume << " _binv[" << ibin << "] " << _binv[ibin] << endl;
 1764 
 1765     const Double_t y = _wgt[ibin] * corr * corrPartial - carry;
 1766     const Double_t t = total + y;
 1767     carry = (t - total) - y;
 1768     total = t;
 1769   }
 1770 
 1771   _vars = varSave;
 1772 
 1773   return total;
 1774 }
 1775 
 1776 
 1777 
 1778 ////////////////////////////////////////////////////////////////////////////////
 1779 /// Fill the transient cache with partial bin volumes with up-to-date
 1780 /// values for the partial volume specified by observables 'dimSet'
 1781 
 1782 void RooDataHist::calculatePartialBinVolume(const RooArgSet& dimSet) const
 1783 {
 1784   // Allocate cache if not yet existing
 1785   vector<Double_t> *pbinv = _pbinvCacheMgr.getObj(&dimSet) ;
 1786   if (pbinv) {
 1787     _pbinv = pbinv ;
 1788     return ;
 1789   }
 1790 
 1791   pbinv = new vector<Double_t>(_arrSize) ;
 1792 
 1793   // Calculate plot bins of components from master index
 1794   Bool_t* selDim = new Bool_t[_vars.getSize()] ;
 1795   _iterator->Reset() ;
 1796   RooAbsArg* v ;
 1797   Int_t i(0) ;
 1798   while((v=(RooAbsArg*)_iterator->Next())) {
 1799     selDim[i++] = dimSet.find(*v) ? kTRUE : kFALSE ;
 1800   }
 1801 
 1802   // Recalculate partial bin volume cache
 1803   Int_t ibin ;
 1804   for (ibin=0 ; ibin<_arrSize ; ibin++) {
 1805     _iterator->Reset() ;
 1806     RooAbsLValue* arg ;
 1807     Int_t j(0), idx(0), tmp(ibin) ;
 1808     Double_t theBinVolume(1) ;
 1809     while((arg=dynamic_cast<RooAbsLValue*>(_iterator->Next()))) {
 1810       idx  = tmp / _idxMult[j] ;
 1811       tmp -= idx*_idxMult[j++] ;
 1812       if (selDim[j-1]) {
 1813    RooAbsLValue* arglv = dynamic_cast<RooAbsLValue*>(arg) ;
 1814    theBinVolume *= arglv->getBinWidth(idx) ;
 1815       }
 1816     }
 1817     (*pbinv)[ibin] = theBinVolume ;
 1818   }
 1819 
 1820   delete[] selDim ;
 1821 
 1822   // Put in cache (which takes ownership)
 1823   _pbinvCacheMgr.setObj(&dimSet,pbinv) ;
 1824 
 1825   // Publicize the array
 1826   _pbinv = pbinv ;
 1827 }
 1828 
 1829 
 1830 
 1831 ////////////////////////////////////////////////////////////////////////////////
 1832 /// Return the number of bins
 1833 
 1834 Int_t RooDataHist::numEntries() const
 1835 {
 1836   return RooAbsData::numEntries() ;
 1837 }
 1838 
 1839 
 1840 
 1841 ////////////////////////////////////////////////////////////////////////////////
 1842 
 1843 Double_t RooDataHist::sumEntries() const
 1844 {
 1845   Int_t i ;
 1846   Double_t n(0), carry(0);
 1847   for (i=0 ; i<_arrSize ; i++) {
 1848     if (!_binValid || _binValid[i]) {
 1849       Double_t y = _wgt[i] - carry;
 1850       Double_t t = n + y;
 1851       carry = (t - n) - y;
 1852       n = t;
 1853     }
 1854   }
 1855   return n ;
 1856 }
 1857 
 1858 
 1859 
 1860 ////////////////////////////////////////////////////////////////////////////////
 1861 /// Return the sum of weights in all entries matching cutSpec (if specified)
 1862 /// and in named range cutRange (if specified)
 1863 /// Return the
 1864 
 1865 Double_t RooDataHist::sumEntries(const char* cutSpec, const char* cutRange) const
 1866 {
 1867   checkInit() ;
 1868 
 1869   if (cutSpec==0 && cutRange==0) {
 1870     return sumEntries();
 1871   } else {
 1872 
 1873     // Setup RooFormulaVar for cutSpec if it is present
 1874     RooFormula* select = 0 ;
 1875     if (cutSpec) {
 1876       select = new RooFormula("select",cutSpec,*get()) ;
 1877     }
 1878 
 1879     // Otherwise sum the weights in the event
 1880     Double_t sumw(0), carry(0);
 1881     Int_t i ;
 1882     for (i=0 ; i<numEntries() ; i++) {
 1883       get(i) ;
 1884       if (select && select->eval()==0.) continue ;
 1885       if (cutRange && !_vars.allInRange(cutRange)) continue ;
 1886 
 1887       if (!_binValid || _binValid[i]) {
 1888    Double_t y = weight() - carry;
 1889    Double_t t = sumw + y;
 1890    carry = (t - sumw) - y;
 1891    sumw = t;
 1892       }
 1893     }
 1894 
 1895     if (select) delete select ;
 1896 
 1897     return sumw ;
 1898   }
 1899 }
 1900 
 1901 
 1902 
 1903 ////////////////////////////////////////////////////////////////////////////////
 1904 /// Reset all bin weights to zero
 1905 
 1906 void RooDataHist::reset()
 1907 {
 1908   // WVE DO NOT CALL RooTreeData::reset() for binned
 1909   // datasets as this will delete the bin definitions
 1910 
 1911   Int_t i ;
 1912   for (i=0 ; i<_arrSize ; i++) {
 1913     _wgt[i] = 0. ;
 1914     _errLo[i] = -1 ;
 1915     _errHi[i] = -1 ;
 1916   }
 1917   _curWeight = 0 ;
 1918   _curWgtErrLo = -1 ;
 1919   _curWgtErrHi = -1 ;
 1920   _curVolume = 1 ;
 1921 
 1922   _cache_sum_valid = kFALSE ;
 1923 
 1924 }
 1925 
 1926 
 1927 
 1928 ////////////////////////////////////////////////////////////////////////////////
 1929 /// Return an argset with the bin center coordinates for
 1930 /// bin sequential number 'masterIdx'. For iterative use.
 1931 
 1932 const RooArgSet* RooDataHist::get(Int_t masterIdx) const
 1933 {
 1934   checkInit() ;
 1935   _curWeight = _wgt[masterIdx] ;
 1936   _curWgtErrLo = _errLo[masterIdx] ;
 1937   _curWgtErrHi = _errHi[masterIdx] ;
 1938   _curSumW2 = _sumw2[masterIdx] ;
 1939   _curVolume = _binv[masterIdx] ;
 1940   _curIndex  = masterIdx ;
 1941   return RooAbsData::get(masterIdx) ;
 1942 }
 1943 
 1944 
 1945 
 1946 ////////////////////////////////////////////////////////////////////////////////
 1947 /// Return a RooArgSet with center coordinates of the bin
 1948 /// enclosing the point 'coord'
 1949 
 1950 const RooArgSet* RooDataHist::get(const RooArgSet& coord) const
 1951 {
 1952   ((RooDataHist*)this)->_vars = coord ;
 1953   return get(calcTreeIndex()) ;
 1954 }
 1955 
 1956 
 1957 
 1958 ////////////////////////////////////////////////////////////////////////////////
 1959 /// Return the volume of the bin enclosing coordinates 'coord'
 1960 
 1961 Double_t RooDataHist::binVolume(const RooArgSet& coord)
 1962 {
 1963   checkInit() ;
 1964   ((RooDataHist*)this)->_vars = coord ;
 1965   return _binv[calcTreeIndex()] ;
 1966 }
 1967 
 1968 
 1969 ////////////////////////////////////////////////////////////////////////////////
 1970 /// Set all the event weight of all bins to the specified value
 1971 
 1972 void RooDataHist::setAllWeights(Double_t value)
 1973 {
 1974   for (Int_t i=0 ; i<_arrSize ; i++) {
 1975     _wgt[i] = value ;
 1976   }
 1977 
 1978   _cache_sum_valid = kFALSE ;
 1979 }
 1980 
 1981 
 1982 
 1983 ////////////////////////////////////////////////////////////////////////////////
 1984 /// Create an iterator over all bins in a slice defined by the subset of observables
 1985 /// listed in sliceArg. The position of the slice is given by otherArgs
 1986 
 1987 TIterator* RooDataHist::sliceIterator(RooAbsArg& sliceArg, const RooArgSet& otherArgs)
 1988 {
 1989   // Update to current position
 1990   _vars = otherArgs ;
 1991   _curIndex = calcTreeIndex() ;
 1992 
 1993   RooAbsArg* intArg = _vars.find(sliceArg) ;
 1994   if (!intArg) {
 1995     coutE(InputArguments) << "RooDataHist::sliceIterator() variable " << sliceArg.GetName() << " is not part of this RooDataHist" << endl ;
 1996     return 0 ;
 1997   }
 1998   return new RooDataHistSliceIter(*this,*intArg) ;
 1999 }
 2000 
 2001 
 2002 ////////////////////////////////////////////////////////////////////////////////
 2003 /// Change the name of the RooDataHist
 2004 
 2005 void RooDataHist::SetName(const char *name)
 2006 {
 2007   if (_dir) _dir->GetList()->Remove(this);
 2008   TNamed::SetName(name) ;
 2009   if (_dir) _dir->GetList()->Add(this);
 2010 }
 2011 
 2012 
 2013 ////////////////////////////////////////////////////////////////////////////////
 2014 /// Change the title of this RooDataHist
 2015 
 2016 void RooDataHist::SetNameTitle(const char *name, const char* title)
 2017 {
 2018   if (_dir) _dir->GetList()->Remove(this);
 2019   TNamed::SetNameTitle(name,title) ;
 2020   if (_dir) _dir->GetList()->Add(this);
 2021 }
 2022 
 2023 
 2024 ////////////////////////////////////////////////////////////////////////////////
 2025 /// Print value of the dataset, i.e. the sum of weights contained in the dataset
 2026 
 2027 void RooDataHist::printValue(ostream& os) const
 2028 {
 2029   os << numEntries() << " bins (" << sumEntries() << " weights)" ;
 2030 }
 2031 
 2032 
 2033 
 2034 
 2035 ////////////////////////////////////////////////////////////////////////////////
 2036 /// Print argument of dataset, i.e. the observable names
 2037 
 2038 void RooDataHist::printArgs(ostream& os) const
 2039 {
 2040   os << "[" ;
 2041   _iterator->Reset() ;
 2042   RooAbsArg* arg ;
 2043   Bool_t first(kTRUE) ;
 2044   while((arg=(RooAbsArg*)_iterator->Next())) {
 2045     if (first) {
 2046       first=kFALSE ;
 2047     } else {
 2048       os << "," ;
 2049     }
 2050     os << arg->GetName() ;
 2051   }
 2052   os << "]" ;
 2053 }
 2054 
 2055 
 2056 
 2057 ////////////////////////////////////////////////////////////////////////////////
 2058 /// Cache the datahist entries with bin centers that are inside/outside the
 2059 /// current observable definitio
 2060 
 2061 void RooDataHist::cacheValidEntries()
 2062 {
 2063   checkInit() ;
 2064 
 2065   if (!_binValid) {
 2066     _binValid = new Bool_t[_arrSize] ;
 2067   }
 2068   TIterator* iter = _vars.createIterator() ;
 2069   RooAbsArg* arg ;
 2070   for (Int_t i=0 ; i<_arrSize ; i++) {
 2071     get(i) ;
 2072     _binValid[i] = kTRUE ;
 2073     iter->Reset() ;
 2074     while((arg=(RooAbsArg*)iter->Next())) {
 2075       // coverity[CHECKED_RETURN]
 2076       _binValid[i] &= arg->inRange(0) ;
 2077     }
 2078   }
 2079   delete iter ;
 2080 
 2081 }
 2082 
 2083 
 2084 ////////////////////////////////////////////////////////////////////////////////
 2085 /// Return true if currently loaded coordinate is considered valid within
 2086 /// the current range definitions of all observables
 2087 
 2088 Bool_t RooDataHist::valid() const
 2089 {
 2090   // If caching is enabled, use the precached result
 2091   if (_binValid) {
 2092     return _binValid[_curIndex] ;
 2093   }
 2094 
 2095   return kTRUE ;
 2096 }
 2097 
 2098 
 2099 
 2100 ////////////////////////////////////////////////////////////////////////////////
 2101 /// Returns true if datasets contains entries with a non-integer weight
 2102 
 2103 Bool_t RooDataHist::isNonPoissonWeighted() const
 2104 {
 2105   for (int i=0 ; i<numEntries() ; i++) {
 2106     if (fabs(_wgt[i]-Int_t(_wgt[i]))>1e-10) return kTRUE ;
 2107   }
 2108   return kFALSE ;
 2109 }
 2110 
 2111 
 2112 
 2113 
 2114 ////////////////////////////////////////////////////////////////////////////////
 2115 /// Print the details on the dataset contents
 2116 
 2117 void RooDataHist::printMultiline(ostream& os, Int_t content, Bool_t verbose, TString indent) const
 2118 {
 2119   RooAbsData::printMultiline(os,content,verbose,indent) ;
 2120 
 2121   os << indent << "Binned Dataset " << GetName() << " (" << GetTitle() << ")" << endl ;
 2122   os << indent << "  Contains " << numEntries() << " bins with a total weight of " << sumEntries() << endl;
 2123 
 2124   if (!verbose) {
 2125     os << indent << "  Observables " << _vars << endl ;
 2126   } else {
 2127     os << indent << "  Observables: " ;
 2128     _vars.printStream(os,kName|kValue|kExtras|kTitle,kVerbose,indent+"  ") ;
 2129   }
 2130 
 2131   if(verbose) {
 2132     if (_cachedVars.getSize()>0) {
 2133       os << indent << "  Caches " << _cachedVars << endl ;
 2134     }
 2135   }
 2136 }
 2137 
 2138 
 2139 
 2140 ////////////////////////////////////////////////////////////////////////////////
 2141 /// Stream an object of class RooDataHist.
 2142 
  void RooDataHist::Streamer(TBuffer &R__b)
  {
 2145    if (R__b.IsReading()) {
 2146 
 2147      UInt_t R__s, R__c;
 2148      Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
 2149 
 2150       if (R__v>2) {
 2151 
 2152    R__b.ReadClassBuffer(RooDataHist::Class(),this,R__v,R__s,R__c);
 2153    initialize(0,kFALSE) ;
 2154 
 2155       } else {
 2156 
 2157    // Legacy dataset conversion happens here. Legacy RooDataHist inherits from RooTreeData
 2158    // which in turn inherits from RooAbsData. Manually stream RooTreeData contents on
 2159    // file here and convert it into a RooTreeDataStore which is installed in the
 2160    // new-style RooAbsData base class
 2161 
 2162    // --- This is the contents of the streamer code of RooTreeData version 2 ---
 2163    UInt_t R__s1, R__c1;
 2164    Version_t R__v1 = R__b.ReadVersion(&R__s1, &R__c1); if (R__v1) { }
 2165 
 2166    RooAbsData::Streamer(R__b);
 2167    TTree* X_tree(0) ; R__b >> X_tree;
 2168    RooArgSet X_truth ; X_truth.Streamer(R__b);
 2169    TString X_blindString ; X_blindString.Streamer(R__b);
 2170    R__b.CheckByteCount(R__s1, R__c1, RooTreeData::Class());
 2171    // --- End of RooTreeData-v1 streamer
 2172 
 2173    // Construct RooTreeDataStore from X_tree and complete initialization of new-style RooAbsData
 2174    _dstore = new RooTreeDataStore(X_tree,_vars) ;
 2175    _dstore->SetName(GetName()) ;
 2176    _dstore->SetTitle(GetTitle()) ;
 2177    _dstore->checkInit() ;
 2178 
 2179    RooDirItem::Streamer(R__b);
 2180    R__b >> _arrSize;
 2181    delete [] _wgt;
 2182    _wgt = new Double_t[_arrSize];
 2183    R__b.ReadFastArray(_wgt,_arrSize);
 2184    delete [] _errLo;
 2185    _errLo = new Double_t[_arrSize];
 2186    R__b.ReadFastArray(_errLo,_arrSize);
 2187    delete [] _errHi;
 2188    _errHi = new Double_t[_arrSize];
 2189    R__b.ReadFastArray(_errHi,_arrSize);
 2190    delete [] _sumw2;
 2191    _sumw2 = new Double_t[_arrSize];
 2192    R__b.ReadFastArray(_sumw2,_arrSize);
 2193    delete [] _binv;
 2194    _binv = new Double_t[_arrSize];
 2195    R__b.ReadFastArray(_binv,_arrSize);
 2196    _realVars.Streamer(R__b);
 2197    R__b >> _curWeight;
 2198    R__b >> _curWgtErrLo;
 2199    R__b >> _curWgtErrHi;
 2200    R__b >> _curSumW2;
 2201    R__b >> _curVolume;
 2202    R__b >> _curIndex;
 2203    R__b.CheckByteCount(R__s, R__c, RooDataHist::IsA());
 2204 
 2205       }
 2206 
 2207    } else {
 2208 
 2209       R__b.WriteClassBuffer(RooDataHist::Class(),this);
 2210    }
 2211 }
 2212 
