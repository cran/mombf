#include <R.h>
#include <Rinternals.h>
#include "cstat.h"
#include "marginals.h"

//Global variables defined for minimization/integration routines
struct marginalPars f2opt_pars, f2int_pars;


//********************************************************************************************
// GENERAL ALGEBRA
//********************************************************************************************

//multiply symmetric A[1..fi][1..fi] * x[sel[0]..sel[fi-1]]
//Note: A is indexed at 1. x and sel are indexed at 0. ans is indexed at 1.
void Asym_xsel(double **A, int fi, double *x, int *sel, double *ans) {
  int _i, _j;
  for (_i=1;_i<=fi;_i++) {
    for (_j=_i, ans[_i]=0; _j<=fi; _j++) { ans[_i]+= A[_i][_j] * x[sel[_j-1]]; }
    for (_j= 1; _j<_i; _j++) { ans[_i]+= A[_j][_i] * x[sel[_j-1]]; }
  } 
}

//Add constant ct to diagonal elements in XtX[sel,sel]. XtX[0..p-1][0..p-1] is formatted as a vector indexed at 0, V[sel[0]..sel[nsel-1]][sel[0]..sel[nsel-1]] as a matrix indexed at 1, sel is indexed at 0
//Note: Only diagonal & upper-diagonal elements in V are set.
void addct2XtX(double *ct, double *XtX, int *sel, int *nsel, int *p, double **V) {
  int i,j;
  for (i=1;i<=(*nsel);i++) { V[i][i]= XtX[sel[i-1]*(*p)+sel[i-1]] + (*ct); }
  for (i=1;i<=(*nsel);i++) {
    for (j=i+1;j<=(*nsel);j++) {
      V[i][j]= XtX[sel[j-1]*(*p) + sel[i-1]];
    }
  }
}


//********************************************************************************************
// GENERAL MARGINAL DENSITY CALCULATION ROUTINES
//********************************************************************************************

void set_marginalPars(struct marginalPars *pars, int *n,int *p,double *y,double *sumy2,double *x,double *XtX,double *ytX,int *method,int *B,double *alpha,double *lambda,double *phi,double *tau,int *r,double *prDeltap,double *parprDeltap, int *logscale) {
  (*pars).n= n;
  (*pars).p= p;
  (*pars).y= y;
  (*pars).sumy2= sumy2;
  (*pars).x= x;
  (*pars).XtX= XtX;
  (*pars).ytX= ytX;
  (*pars).method= method;
  (*pars).B= B;
  (*pars).alpha= alpha;
  (*pars).lambda= lambda;
  (*pars).phi= phi;
  (*pars).tau= tau;
  (*pars).r= r;
  (*pars).prDeltap= prDeltap;
  (*pars).parprDeltap= parprDeltap;
  (*pars).logscale= logscale;
}

void set_f2opt_pars(double *m, double **S, double *XtX, double *ytX, double *phi, double *tau, int *r, int *n, int *p, int *sel, int *nsel) {
  f2opt_pars.m= m;
  f2opt_pars.S= S;
  f2opt_pars.XtX= XtX;
  f2opt_pars.ytX= ytX;
  f2opt_pars.phi= phi;
  f2opt_pars.tau= tau;
  f2opt_pars.r= r;
  f2opt_pars.n= n;
  f2opt_pars.p= p;
  f2opt_pars.sel= sel;
  f2opt_pars.nsel= nsel;
}

void set_f2int_pars(double *XtX, double *ytX, double *tau, int *n, int *p, int *sel, int *nsel, double *y, double *sumy2, int *method, int *B, double *alpha, double *lambda, int *logscale) {
  f2int_pars.XtX= XtX;
  f2int_pars.ytX= ytX;
  f2int_pars.tau= tau;
  f2int_pars.n= n;
  f2int_pars.p= p;
  f2int_pars.sel= sel;
  f2int_pars.nsel= nsel;
  f2int_pars.y= y;
  f2int_pars.sumy2= sumy2;
  f2int_pars.method= method;
  f2int_pars.B= B;
  f2int_pars.alpha= alpha;
  f2int_pars.lambda= lambda;
  f2int_pars.logscale= logscale;
}


//********************************************************************************************
// MODEL SELECTION ROUTINES
//********************************************************************************************

//modelSelectionC: Gibbs sampler for model selection in linear regression for several choices of prior distribution
//Input parameters
// - knownphi: is residual variance phi known?
// - priorCoef: 0 for product MOM, 1 for product iMOM, 2 for product eMOM
// - priorDelta: 0 for uniform, 1 for binomial, 2 for binomial with beta hyper-prior for success prob
// - niter: number of Gibbs iterations
// - ndeltaini: length of deltaini
// - deltaini: vector with indexes of covariates initially in the model (both deltaini and its indexes must be indexed at 0)
// - verbose: set verbose==1 to print iteration progress every 10% of the iterations
// - pars: struct of type marginalPars containing parameters needed to evaluate the marginal density of the data & prior on model space
//Output
// - postSample: matrix with niter rows and p columns with posterior samples for covariate inclusion/exclusion (formatted as a vector in column order)
// - postOther: matrix with niter rows and nOther columns with posterior samples for other model parameters
// - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
// - postMode: model with highest posterior probability amongst all those visited
// - postModeProb: unnormalized posterior prob of posterior mode (log scale)
// - postProb: unnormalized posterior prob of each visited model (log scale)

SEXP modelSelectionCI(SEXP SpostSample, SEXP SpostOther, SEXP Smargpp, SEXP SpostMode, SEXP SpostModeProb, SEXP SpostProb, SEXP Sknownphi, SEXP SpriorCoef, SEXP Sniter, SEXP Sthinning, SEXP Sburnin, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose) {
  int logscale=1;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars, INTEGER(Sn), INTEGER(Sp), REAL(Sy), REAL(Ssumy2), REAL(Sx), REAL(SXtX), REAL(SytX), INTEGER(Smethod), INTEGER(SB), REAL(Salpha),REAL(Slambda), REAL(Sphi), REAL(Stau), INTEGER(Sr), REAL(SprDeltap), REAL(SparprDeltap), &logscale);
  modelSelectionC(INTEGER(SpostSample), REAL(SpostOther), REAL(Smargpp), INTEGER(SpostMode), REAL(SpostModeProb), REAL(SpostProb), INTEGER(Sknownphi), INTEGER(SpriorCoef), INTEGER(SpriorDelta), INTEGER(Sniter), INTEGER(Sthinning), INTEGER(Sburnin), INTEGER(Sndeltaini), INTEGER(Sdeltaini), INTEGER(Sverbose), &pars);

  PROTECT(ans = allocVector(REALSXP, 1));
  *REAL(ans)= 1.0;
  UNPROTECT(1);
  return ans;
}

void modelSelectionC(int *postSample, double *postOther, double *margpp, int *postMode, double *postModeProb, double *postProb, int *knownphi, int *prCoef, int *prDelta, int *niter, int *thinning, int *burnin, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars) {
  int i, j, k, *sel, *selnew, *selaux, nsel, nselnew, niter10, niterthin, savecnt, ilow, iupper;
  double currentJ, newM, newP, newJ, ppnew, u;
  pt2margFun marginalFunction=NULL, priorFunction=NULL; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);

  marginalFunction= set_marginalFunction(prCoef, knownphi);
  priorFunction= set_priorFunction(prDelta);
 
  sel= ivector(0,*(*pars).p); selnew= ivector(0,*(*pars).p);

  //Initialize
  if (*verbose ==1) Rprintf("Running Gibbs sampler");
  niterthin= floor((*niter - *burnin +.0)/(*thinning+.0));
  if (*niter >10) { niter10= *niter/10; } else { niter10= 1; }
  for (j=0; j< *(*pars).p; j++) { margpp[j]= 0; }
  nsel= *ndeltaini;
  for (j=0; j< nsel; j++) { sel[j]= deltaini[j]; postSample[deltaini[j]*niterthin]= postMode[deltaini[j]]= 1; }
  if ((*prDelta)==2) { postOther[0]= *(*pars).prDeltap; }
  currentJ= marginalFunction(sel,&nsel,pars) + priorFunction(sel,&nsel,pars);
  postProb[0]= *postModeProb= currentJ;
  if (*burnin >0) { ilow=-(*burnin); savecnt=0; iupper= *niter - *burnin +1; } else { ilow=1; savecnt=1; iupper= *niter; } //if no burnin, start at i==1 & save initial value

  //Iterate
  for (i=ilow; i< iupper; i++) {
    for (j=0; j< *(*pars).p; j++) {
      sel2selnew(j,sel,&nsel,selnew,&nselnew); //copy sel into selnew, adding/removing jth element
      newM= marginalFunction(selnew,&nselnew,pars);
      newP= priorFunction(selnew,&nselnew,pars);
      newJ= newM + newP;
      if (newJ > *postModeProb) {   //update posterior mode
        *postModeProb= newJ;
        for (k=0; k< *(*pars).p; k++) { postMode[k]= 0; }
        for (k=0; k< nselnew; k++) { postMode[selnew[k]]= 1; } 
      }
      ppnew= 1.0/(1.0+exp(currentJ-newJ));
      if (i>=0) { if (nselnew>nsel) { margpp[j]+= ppnew; } else { margpp[j]+= (1-ppnew); } }
      u= runif();
      if (u<ppnew) {  //update model indicator
        selaux= sel; sel=selnew; selnew=selaux; nsel=nselnew; currentJ= newJ;
      }
    }  //end j for
    if ((i>0) && ((i%(*thinning))==0)) {
      if ((*prDelta)==2) { 
        *(*pars).prDeltap= rbetaC(nsel + (*pars).parprDeltap[0], *(*pars).p - nsel + (*pars).parprDeltap[1]);
        postOther[savecnt]= *(*pars).prDeltap;
      }
      for (j=0; j<nsel; j++) { postSample[sel[j]*niterthin+savecnt]= 1; }
      postProb[savecnt]= currentJ;
      savecnt++;
    }
    if ((*verbose==1) && ((i%niter10)==0)) { Rprintf("."); }
  }
  if (iupper>ilow) { for (j=0; j< *(*pars).p; j++) { margpp[j] /= (iupper-imax_xy(0,ilow)+.0); } } //from sum to average
  if (*verbose==1) printf(" Done.\n");

  free_ivector(sel,0,*(*pars).p); free_ivector(selnew,0,*(*pars).p);
}

//greedyVarSelC: greedy search for posterior mode in variable selection.
//               Similar to Gibbs sampling, except that deterministic updates are made iff there is an increase in post model prob
//               The scheme proceeds until no variable is included/excluded or niter iterations are reached
// Input arguments: same as in modelSelectionC.

SEXP greedyVarSelCI(SEXP SpostMode, SEXP SothersMode, SEXP SpostModeProb, SEXP Sknownphi, SEXP SpriorCoef, SEXP Sniter, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose) {
  int logscale=1;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars, INTEGER(Sn), INTEGER(Sp), REAL(Sy), REAL(Ssumy2), REAL(Sx), REAL(SXtX), REAL(SytX), INTEGER(Smethod), INTEGER(SB), REAL(Salpha),REAL(Slambda), REAL(Sphi), REAL(Stau), INTEGER(Sr), REAL(SprDeltap), REAL(SparprDeltap), &logscale);
  greedyVarSelC(INTEGER(SpostMode),REAL(SothersMode),REAL(SpostModeProb),INTEGER(Sknownphi),INTEGER(SpriorCoef),INTEGER(SpriorDelta),INTEGER(Sniter),INTEGER(Sndeltaini),INTEGER(Sdeltaini),INTEGER(Sverbose),&pars);
  PROTECT(ans = allocVector(REALSXP, 1));
  *REAL(ans)= 1.0;
  UNPROTECT(1);
  return ans;

}

void greedyVarSelC(int *postMode, double *othersMode, double *postModeProb, int *knownphi, int *prCoef, int *prDelta, int *niter, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars) {
  int i, j, *sel, *selnew, *selaux, nsel, nselnew, nchanges;
  double newJ, a, b;
  pt2margFun marginalFunction=NULL, priorFunction=NULL; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);

  marginalFunction= set_marginalFunction(prCoef, knownphi);
  priorFunction= set_priorFunction(prDelta);
  sel= ivector(0,*(*pars).p); selnew= ivector(0,*(*pars).p);

  //Initialize
  if (*verbose==1) Rprintf("Greedy searching posterior mode... ");
  for (j=0, nsel=*ndeltaini; j< nsel; j++) { sel[j]= deltaini[j]; postMode[deltaini[j]]= 1; }
  *postModeProb= marginalFunction(sel,&nsel,pars) + priorFunction(sel,&nsel,pars);

  //Iterate
  for (i=0, nchanges=1; (i< *niter) && (nchanges>0); i++) {
    for (j=0, nchanges=0; j< *(*pars).p; j++) {
      sel2selnew(j,sel,&nsel,selnew,&nselnew); //copy sel into selnew, adding/removing jth element
      newJ= marginalFunction(selnew,&nselnew,pars) + priorFunction(selnew,&nselnew,pars);
      if (newJ > *postModeProb) {
        *postModeProb= newJ;  //update post mode prob
        if (postMode[j]==0) { postMode[j]= 1; } else { postMode[j]= 0; }  //update post mode
        //for (k=0; k< *(*pars).p; k++) { postMode[k]= 0; }
        //for (k=0; k< nselnew; k++) { postMode[selnew[k]]= 1; } 
        selaux= sel; sel=selnew; selnew=selaux; nsel=nselnew; //update model indicator
        nchanges++;
      }
    } //end j for
    if ((*prDelta)==2) { 
      a= nsel + (*pars).parprDeltap[0]; b= *(*pars).p - nsel + (*pars).parprDeltap[1];
      if ((a>1) && (b>1)) { *(*pars).prDeltap= (a-1)/(a+b-2); } else { *(*pars).prDeltap= a/(a+b); }
    }
  }
  othersMode[0]= *(*pars).prDeltap;
  if (*verbose==1) Rprintf("Done.\n");

  free_ivector(sel,0,*(*pars).p); free_ivector(selnew,0,*(*pars).p);
}


pt2margFun set_marginalFunction(int *prCoef, int *knownphi) {
  //Returns pointer to function to compute the marginal density of the data for a given model indicator
  // - prCoef: 0 for product MOM, 1 for product iMOM, 2 for product eMOM
  // - knownphi: 1 if residual variance phi is know, 0 otherwise. 
  // Note: if phi known, when actually calling the returned pt2margFun, phi must be set in the parameter of type struct marginalPars *
  pt2margFun ans=NULL;
  if (*prCoef==0) {
    if (*knownphi==1) { ans= pmomMarginalKC; } else { ans= pmomMarginalUC; }
  } else if (*prCoef==1) {
    if (*knownphi==1) { ans= pimomMarginalKC; } else { ans= pimomMarginalUC; }
  } else if (*prCoef==2) {
    if (*knownphi==1) { ans= pemomMarginalKC; } else { ans= pemomMarginalUC; }
  }
  return ans;
}

pt2margFun set_priorFunction(int *prDelta) {
  //Returns pointer to function to compute the prior probability of a model indicator
  // - prDelta: 0 for uniform, 1 for binomial, 2 for beta-binomial
  pt2margFun ans=NULL;
  if (*prDelta==0) { ans= unifPrior; } else if ((*prDelta==1) || (*prDelta==2)) { ans= binomPrior; }
  return ans;
}

void sel2selnew(int newelem,int *sel,int *nsel,int *selnew,int *nselnew) {
//Copy sel into selnew. 
// - If j in sel, don't copy it in selnew and set nselnew=nsel-1. 
// - If j not in sel, add it to selnew and set nselnew=nsel+1.
  int i, found;
  for (i=0, found=0; (i< *nsel) && (found==0); i++) { selnew[i]= sel[i]; found= (newelem==sel[i]); }
  if (found==0) { //add newelem
    selnew[i]= newelem; (*nselnew)= (*nsel)+1;
  } else {  //remove new elem
    for (i=i; i< *nsel; i++) { selnew[i-1]= sel[i]; }
    (*nselnew)= (*nsel)-1;
  }
}


//********************************************************************************************
// PRIORS ON MODEL SPACE (always return on log scale)
//********************************************************************************************

double unifPrior(int *sel, int *nsel, struct marginalPars *pars) { return 0.0; }

double binomPrior(int *sel, int *nsel, struct marginalPars *pars) {
  //nsel ~ Binom(p,prDeltap)
  return dbinomial(*nsel,*(*pars).p,*(*pars).prDeltap,1);
} 


//*************************************************************************************
// PRODUCT MOM ROUTINES
//*************************************************************************************

double f2opt_mom(double *th) {
  return fmomNegC_non0(th+1,f2opt_pars.m+1,f2opt_pars.S,f2opt_pars.phi,f2opt_pars.tau,f2opt_pars.r,f2opt_pars.n,f2opt_pars.nsel);
}

//Note: th and m are assumed to be indexed at 0; S indexed at 1
double fmomNegC_non0(double *th, double *m, double **S, double *phi, double *tau, int *r, int *n, int *nsel) {
  int i;
  double ans, sumlogth, *z;
  z= dvector(0,*nsel);
  for (i=0, sumlogth=0; i<(*nsel); i++) { sumlogth+= log(th[i]*th[i]); z[i]= th[i]-m[i]; }
  ans= .5*quadratic_xtAx(z-1,S,1,*nsel)/(*phi) - (*r +.0)*sumlogth;
  //ans= .5*quadratic_xtAselx(z, XtXplusct, p, nsel, sel)/(*phi) - (*r +.0)*sumlogth;
  free_dvector(z,0,*nsel);
  return ans;
}

void fppmomNegC_non0(double **ans, double *th, double **S, double *phi, double *tau, int *r, int *n, int *nsel) {
  int i, j;
  for (i=1; i<=(*nsel); i++) { ans[i][i]= S[i][i]/(*phi) + 2.0*(*r)/(th[i]*th[i]); }
  for (i=1; i<=(*nsel); i++) { for (j=i+1; j<=(*nsel); j++) { ans[i][j]= ans[j][i]= S[i][j]/(*phi); } }
}

void momIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *n, int *nsel, double *m, double **S, double *detS, double *phi, double *tau, int *r, int *logscale) {
  int i, emptyint, iter, maxit=50;
  double emptydouble, ftol= 1.0e-2, **Vopt, detVopt, **dirth;

  Vopt= dmatrix(1,*nsel,1,*nsel); dirth= dmatrix(1,*nsel,1,*nsel);
  set_f2opt_pars(m,S,&emptydouble,&emptydouble,phi,tau,r,n,nsel,&emptyint,nsel);

  //Minimization
  for (i=1; i<=(*nsel); i++) { thopt[i]= m[i]; }  //init
  ddiag(dirth,1,*nsel);
  minimize(thopt, dirth, *nsel, ftol, &iter, fopt, f2opt_mom, maxit);

  //Laplace approx
  fppmomNegC_non0(Vopt,thopt,S,phi,tau,r,n,nsel);
  invdet_posdef(Vopt,*nsel,Voptinv,&detVopt);

  (*ILaplace)= -(*fopt) + .5*(log(*detS)-log(detVopt)- (*nsel)*log(*phi)) ;

  if ((*logscale)!=1) { (*ILaplace)= exp(*ILaplace); }
  free_dmatrix(Vopt,1,*nsel,1,*nsel); free_dmatrix(dirth,1,*nsel,1,*nsel);
}


//Monter Carlo evaluation of E(prod(z^(2*r))), where z ~ N(m,Sinv)
double MC_mom_normal(double *m,double **Sinv,int *r,int *nsel, int *B) {
  int i;
  double **cholSinv, *thsim, ans, normfac;

  thsim= dvector(1,*nsel);
  cholSinv= dmatrix(1,*nsel,1,*nsel);
  choldc(Sinv,*nsel,cholSinv); //compute cholesky decomposition
  normfac= rsumlogsq(m,r,nsel);
  for (i=0, ans=0; i<(*B); i++) {
    rmvnormC(thsim,*nsel,m,cholSinv);
    ans+= exp(rsumlogsq(thsim,r,nsel) - normfac); 
  }
  ans= log(ans/(*B +.0)) + normfac;

  free_dvector(thsim,1,*nsel);
  free_dmatrix(cholSinv,1,*nsel,1,*nsel);
  return ans;
}

//Monter Carlo evaluation of E(prod(z^(2*r))), where z ~ T_nu(m,Sinv)
double MC_mom_T(double *m,double **Sinv,int *nu,int *r,int *nsel, int *B) {
  int i;
  double **cholSinv, *thsim, ans, normfac;

  thsim= dvector(1,*nsel);
  cholSinv= dmatrix(1,*nsel,1,*nsel);
  choldc(Sinv,*nsel,cholSinv); //compute cholesky decomposition
  normfac= rsumlogsq(m,r,nsel);
  for (i=0, ans=0; i<(*B); i++) {
    rmvtC(thsim,*nsel,m,cholSinv,*nu);
    ans+= exp(rsumlogsq(thsim,r,nsel) - normfac); 
  }
  ans= log(ans/(*B +.0)) + normfac;

  free_dvector(thsim,1,*nsel);
  free_dmatrix(cholSinv,1,*nsel,1,*nsel);
  return ans;
}


// PRODUCT MOMENT MARGINALS
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - r: MOM power parameter
// - method==0 for Laplace; method==1 for Monte Carlo; method==2 for plug-in
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale

SEXP pmomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale) {
  struct marginalPars pars;
  double *rans, emptydouble=0;
  SEXP ans;

  set_marginalPars(&pars,INTEGER(Sn),INTEGER(Sp),REAL(Sy),REAL(Ssumy2),&emptydouble,REAL(SXtX),REAL(SytX),INTEGER(Smethod),INTEGER(SB),&emptydouble,&emptydouble,REAL(Sphi),REAL(Stau),INTEGER(Sr),&emptydouble,&emptydouble,INTEGER(Slogscale));
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pmomMarginalKC(INTEGER(Ssel),INTEGER(Snsel),&pars);
  UNPROTECT(1);
  return ans;
}

// Function to compute r * sum(log(th^2))
double rsumlogsq(double *th, int *r, int *nsel) {
  int i; double ans;
  for (i=1, ans=0; i<=(*nsel); i++) { ans+= log(th[i]*th[i]); }
  ans*= (*r);
  return(ans);
}

double pmomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  int i,j;
  double *m, s, **S, **Sinv, detS, num, den, logtau= log(*(*pars).tau), tauinv= 1.0/(*(*pars).tau), logphi= log(*(*pars).phi), ans, *thopt, **Voptinv, fopt;

  if (*nsel ==0) {
    m= dvector(1,1);
    m[1]=0; s= sqrt(*(*pars).phi);
    ans= dnormC_jvec((*pars).y,*(*pars).n,m[1],s,1);
    free_dvector(m,1,1);
  } else {
    m= dvector(1,*nsel);
    S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&tauinv,(*pars).XtX,sel,nsel,(*pars).p,S);
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);

    num= -.5*(*(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel))/(*(*pars).phi);
    den= .5*((*(*pars).n +.0)*(LOG_M_2PI+logphi) + log(detS) + (*nsel)*logtau) + (*nsel)*(*(*pars).r)*(logtau+logphi+ldoublefact(2*(*(*pars).r)-1));
    if (*(*pars).method ==0) { //Laplace
      thopt= dvector(1,*nsel); Voptinv= dmatrix(1,*nsel,1,*nsel);
      momIntegralApproxC(&ans,thopt,Voptinv,&fopt,(*pars).n,nsel,m,S,&detS,(*pars).phi,(*pars).tau,(*pars).r,(*pars).logscale);
      free_dvector(thopt,1,*nsel); free_dmatrix(Voptinv,1,*nsel,1,*nsel);
    } else if (*(*pars).method ==1) { //MC
      for (i=1; i<=(*nsel); i++) { Sinv[i][i]= (*(*pars).phi)*Sinv[i][i]; for (j=i+1; j<=(*nsel); j++) { Sinv[i][j]=Sinv[j][i]= (*(*pars).phi)*Sinv[i][j]; } }
      ans= MC_mom_normal(m,Sinv,(*pars).r,nsel,(*pars).B);
    } else if (*(*pars).method ==2) { //Plug-in
      ans= rsumlogsq(m,(*pars).r,nsel);
    }
    ans+= num - den;
    free_dvector(m,1,*nsel);
    free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;  
}


SEXP pmomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda) {
  double *rans, emptydouble=0;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars,INTEGER(Sn),INTEGER(Sp),REAL(Sy),REAL(Ssumy2),REAL(Sx),REAL(SXtX),REAL(SytX),INTEGER(Smethod),INTEGER(SB),REAL(Salpha),REAL(Slambda),&emptydouble,REAL(Stau),INTEGER(Sr),&emptydouble,&emptydouble,INTEGER(Slogscale));
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pmomMarginalUC(INTEGER(Ssel), INTEGER(Snsel), &pars);
  UNPROTECT(1);
  return ans;
}


double pmomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  int i, j, nu;
  double num, den, ans, term1, *m, **S, **Sinv, detS, *thopt, **Voptinv, fopt, phiadj, tauinv= 1.0/(*(*pars).tau), nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda);
  if (*nsel ==0) {
    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*LOG_M_PI + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);
  } else {
    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&tauinv,(*pars).XtX,sel,nsel,(*pars).p,S);
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    nuhalf= (*(*pars).r)*(*nsel) + .5*(*(*pars).n + *(*pars).alpha);
    nu= 2*nuhalf;

    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2) - log(*(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel)));
    den= (*nsel)*ldoublefact(2*(*(*pars).r)-1.0) + .5*(*(*pars).n * LOG_M_2PI + log(detS)) + (*nsel)*(.5 + *(*pars).r)*log(*(*pars).tau) + gamln(&alphahalf);
    if (*(*pars).method ==0) { //Laplace
      thopt= dvector(1,*nsel); Voptinv= dmatrix(1,*nsel,1,*nsel);
      phiadj= (nu+.0)/(nu-2.0);
      momIntegralApproxC(&ans,thopt,Voptinv,&fopt,(*pars).n,nsel,m,S,&detS,&phiadj,(*pars).tau,(*pars).r,(*pars).logscale);
      free_dvector(thopt,1,*nsel); free_dmatrix(Voptinv,1,*nsel,1,*nsel);
    } else if (*(*pars).method ==1) {  //MC
      term1= (*(*pars).lambda + *(*pars).sumy2 - quadratic_xseltAxsel((*pars).ytX,Sinv,1,nsel,sel))/(nu+.0);
      for (i=1; i<= *nsel; i++) { for (j=i; j<= *nsel; j++) { Sinv[i][j]= Sinv[j][i]= Sinv[i][j]*term1; } } //Vinv matrix
      ans= MC_mom_T(m,Sinv,&nu,(*pars).r,nsel,(*pars).B);
    } else if (*(*pars).method ==2) {  //Plug-in
      ans= rsumlogsq(m,(*pars).r,nsel);
    }
    ans+= num - den;
    free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}


//********************************************************************************************
// PRODUCT IMOM ROUTINES
//********************************************************************************************

//fimomNeg: minus log integrand of the function needed to compute product iMOM marginal density of the data conditional under a given linear model
//
// fimomNegC
// Input
// - th: theta value at which to evaluate the function (includes coef for all covariates, even those excluded in the current model)
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - n: sample size (length of y)
// - p: number of columns in XtX
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// Output: scalar evaluating the function for a single value of theta
double fimomNegC(double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i;
  double ans, ytXth, sumlogth, suminvth, th2;
  for (i=0, ytXth=0, sumlogth=0, suminvth=0; i<(*nsel); i++) {
    ytXth+= ytX[sel[i]] * th[sel[i]];
    th2= th[sel[i]] * th[sel[i]];
    suminvth+= 1/th2;
    sumlogth+= log(th2);
  }
  ans= .5*(quadratic_xseltAselxsel(th, XtX, p, nsel, sel) - 2*ytXth)/(*phi) + (*tau)*(*phi)*suminvth + sumlogth;
  return ans;
}


double f2opt_imom(double *th) {
  double ans;
  ans= fimomNegC_non0(th+1,f2opt_pars.XtX,f2opt_pars.ytX,f2opt_pars.phi,f2opt_pars.tau,f2opt_pars.n,f2opt_pars.p,f2opt_pars.sel,f2opt_pars.nsel);
  return(ans);
}

double fimomNegC_non0(double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
//same as fimomNegC but loops over all elem in th (i.e. th has length *nsel and contains non-zero elements only). th is indexed at 0.
  int i;
  double ans, ytXth, sumlogth, suminvth, th2;
  for (i=0, ytXth=0, sumlogth=0, suminvth=0; i<(*nsel); i++) {
    ytXth+= ytX[sel[i]] * th[i];
    th2= th[i] * th[i];
    suminvth+= 1/th2;
    sumlogth+= log(th2);
  }
  ans= .5*(quadratic_xtAselx(th, XtX, p, nsel, sel) - 2*ytXth)/(*phi) + (*tau)*(*phi)*suminvth + sumlogth;
  return ans;
}


// fimomNegC_mat
// - ans: output vector (length nrow) evaluating the function for several theta values
// - nrow: number of rows in the th matrix
// - th: matrix (nrow rows, nsel columns) with non-zero th elements. Given as a vector in row order.
// - rest as for fimomNegC
/*void fimomNegC_non0_mat(double *ans, int *nrow, double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i, j;
  double ytXth, sumlogth, suminvth, th2;
  for (j=0; j<(*nrow); j++) {
    for (i=0, ytXth=0, sumlogth=0, suminvth=0; i<(*nsel); i++) {
      ytXth+= ytX[sel[i]] * th[j*(*nsel)+i];
      th2= th[j*(*nsel)+i] * th[j*(*nsel)+i];
      suminvth+= 1/th2;
      sumlogth+= log(th2);
    }
    ans[j]= .5*(quadratic_xtAselx(th + j*(*nsel), XtX, p, nsel, sel) - 2*ytXth)/(*phi) + (*tau)*(*phi)*suminvth + sumlogth;
  }
}
*/


//Hessian of fimomNegC
// - ans: hessian matrix evaluated at th (indexed at 1, i.e. ans[1:(*nsel)][1:(*nsel)])
// - th: th[1:(*nsel)] indicates point at which to evaluate the hessian.
// - Other arguments as in fimomNegC_non0
void fppimomNegC_non0(double **ans, double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i, j;
  double th2;

  for (i=1; i<=(*nsel); i++) {
    th2= th[i]*th[i];
    ans[i][i]= XtX[sel[i-1]*(*p)+sel[i-1]]/(*phi) + 6.0*(*tau)*(*phi)/(th2*th2) - 1.0/th2;
  }
  for (i=1; i<=(*nsel); i++) {
    for (j=i+1; j<=(*nsel); j++) {
      ans[i][j]= ans[j][i]= XtX[sel[i-1]*(*p)+sel[j-1]]/(*phi);
    }
  }
}


void imomIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *sel, int *nsel, int *n, int *p, double *XtX, double *ytX, double *phi, double *tau, int *logscale) {
  int iter, maxit=50, emptyint;
  double **V, **Vinv, ftol= 1.0e-2, **dirth, **Vopt, detVopt, emptydouble=0, **emptymatrix;

  V= dmatrix(1,*nsel,1,*nsel); Vinv= dmatrix(1,*nsel,1,*nsel); Vopt= dmatrix(1,*nsel,1,*nsel); dirth= dmatrix(1,*nsel,1,*nsel);
  emptymatrix= dmatrix(1,1,1,1);
  //Initialize
  addct2XtX(tau,XtX,sel,nsel,p,V); //add tau to XtX diagonal, store in V
  inv_posdef_upper(V,*nsel,Vinv);
  Asym_xsel(Vinv,*nsel,ytX,sel,thopt);  //product Vinv * selected elements in ytX
  ddiag(dirth,1,*nsel);
  set_f2opt_pars(&emptydouble,emptymatrix,XtX,ytX,phi,tau,&emptyint,n,p,sel,nsel);

  //Minimization
  minimize(thopt, dirth, *nsel, ftol, &iter, fopt, f2opt_imom, maxit);

  //Laplace approx
  fppimomNegC_non0(Vopt,thopt,XtX,ytX,phi,tau,n,p,sel,nsel);
  invdet_posdef(Vopt,*nsel,Voptinv,&detVopt);
  (*ILaplace)= -(*fopt) - 0.5*log(detVopt);

  free_dmatrix(V,1,*nsel,1,*nsel); free_dmatrix(Vinv,1,*nsel,1,*nsel); free_dmatrix(Vopt,1,*nsel,1,*nsel); free_dmatrix(dirth,1,*nsel,1,*nsel);
  free_dmatrix(emptymatrix,1,1,1,1);
  if ((*logscale)!=1) { (*ILaplace)= exp(*ILaplace); }
}


//Product iMOM marginal density for known phi
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - method: method to approximate the marginal. method==0 for Laplace approximation (may underestimate the true value), method==1 for Importance Sampling Monte Carlo 
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale

SEXP pimomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale) {
  int *sel=INTEGER(Ssel), *nsel=INTEGER(Snsel), *n=INTEGER(Sn), *p=INTEGER(Sp), *method=INTEGER(Smethod), *B=INTEGER(SB), *logscale=INTEGER(Slogscale), r=1;
  double *y=REAL(Sy), *sumy2=REAL(Ssumy2), *XtX=REAL(SXtX), *ytX=REAL(SytX), *phi=REAL(Sphi), *tau=REAL(Stau), *rans, emptydouble=0;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars,n,p,y,sumy2,&emptydouble,XtX,ytX,method,B,&emptydouble,&emptydouble,phi,tau,&r,&emptydouble,&emptydouble,logscale);
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pimomMarginalKC(sel, nsel, &pars);
  UNPROTECT(1);
  return ans;
}


double pimomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  int one=1;
  double k, ans, m, s, ILaplace, *thopt, **Voptinv, fopt;
  thopt= dvector(1,*nsel);
  Voptinv= dmatrix(1,*nsel,1,*nsel);
  if ((*nsel)==0) {
    m= 0;
    s= sqrt(*(*pars).phi);
    ans= dnormC_jvec((*pars).y,*(*pars).n,m,s,1);
  } else {
    imomIntegralApproxC(&ILaplace,thopt,Voptinv,&fopt,sel,nsel,(*pars).n,(*pars).p,(*pars).XtX,(*pars).ytX,(*pars).phi,(*pars).tau,&one);
    k= .5*((*nsel)*log(*(*pars).tau) - (*(*pars).sumy2)/(*(*pars).phi) - (*(*pars).n +.0) *LOG_M_2PI - (*(*pars).n - *nsel)*log(*(*pars).phi) - (*nsel)*LOG_M_PI);
    if ((*(*pars).method)==0) {
      ans= k + ILaplace;
    } else {
      ans= k + IS_imom(thopt,Voptinv,sel,nsel,(*pars).n,(*pars).p,(*pars).XtX,(*pars).ytX,(*pars).phi,(*pars).tau,(*pars).B);
    }
  }
  if ((*(*pars).logscale)!=1) { ans= exp(ans); }
  free_dvector(thopt,1,*nsel);
  free_dmatrix(Voptinv,1,*nsel,1,*nsel);
  return(ans);
}


//Evaluation of iMOM integral via Importance Sampling (result is returned in log-scale)
double IS_imom(double *thopt, double **Voptinv, int *sel, int *nsel, int *n, int *p, double *XtX, double *ytX, double *phi, double *tau, int *B) {
  int i,j;
  double *sdprop, **Vprop, *sopt, **cholVprop, **cholVpropinv, detVpropinv, *mprop, *thsim, *logr, maxlogr, ans;

  sdprop= dvector(1,*nsel); sopt= dvector(1,*nsel);
  mprop= dvector(1,*nsel); thsim= dvector(1, *nsel);
  logr= dvector(0,999);
  Vprop= dmatrix(1,*nsel,1,*nsel); cholVprop= dmatrix(1,*nsel,1,*nsel); cholVpropinv= dmatrix(1,*nsel,1,*nsel);

  for (i=1; i<=(*nsel); i++) {
    mprop[i]= 0;
    sopt[i]= sqrt(Voptinv[i][i]);
    sdprop[i]= .5*fabs(thopt[i] + 2*dsign(thopt[i])*sopt[i]);
  }
  for (i=1; i<=(*nsel); i++) {
    for (j=i; j<=(*nsel); j++) {
      Vprop[i][j]= Vprop[j][i]= sdprop[i]*sdprop[j]*Voptinv[i][j]/(sopt[i]*sopt[j]);
    }
  }
  choldc(Vprop,*nsel,cholVprop);
  choldc_inv(Vprop,*nsel,cholVpropinv);
  detVpropinv= choldc_det(cholVpropinv, *nsel);
  rmvtC(thsim, *nsel, mprop, cholVprop, 1);
  maxlogr= logr[0]= -fimomNegC_non0(thsim+1,XtX,ytX,phi,tau,n,p,sel,nsel) - dmvtC(thsim,*nsel,mprop,cholVpropinv,detVpropinv,1,1);
  for (i=1;i<1000;i++) {
    rmvtC(thsim, *nsel, mprop, cholVprop, 1);
    logr[i]= -fimomNegC_non0(thsim+1,XtX,ytX,phi,tau,n,p,sel,nsel) - dmvtC(thsim,*nsel,mprop,cholVpropinv,detVpropinv,1,1);
    if (logr[i]>maxlogr) { maxlogr= logr[i]; }
  }
  for (i=0, ans=0; i<1000; i++) { ans+= exp(logr[i]-maxlogr+500); }
  for (i=1000;i<(*B);i++) {
    rmvtC(thsim, *nsel, mprop, cholVprop, 1);
    ans+= exp(-fimomNegC_non0(thsim+1,XtX,ytX,phi,tau,n,p,sel,nsel) - dmvtC(thsim,*nsel,mprop,cholVpropinv,detVpropinv,1,1) -maxlogr+500);
  }
  ans= log(ans/(.0+ (*B))) + maxlogr-500;

  free_dvector(sdprop,1,*nsel); free_dvector(sopt,1,*nsel);
  free_dvector(mprop, 1,*nsel); free_dvector(thsim, 1, *nsel);
  free_dvector(logr,0,999);
  free_dmatrix(Vprop,1,*nsel,1,*nsel); free_dmatrix(cholVprop,1,*nsel,1,*nsel); free_dmatrix(cholVpropinv,1,*nsel,1,*nsel);
  return(ans);
}


//Product iMOM marginal density for unknown phi
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - x: design matrix (includes all covariates, even those excluded under the current model)
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - tau: prior dispersion parameter
// - method: method to approximate the integral for known phi. Integral wrt phi is performed numerically (Romberg). 0 for Laplace approx which may underestimate true value, 1 for exact evaluation which can be very computationally expensive. 2 ('Hybrid') uses Laplace approx for all phi values, and then estimates the Laplace error by doing a single exact evaluation for a value of phi close to the posterior mode
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale
// - alpha, lambda: prior for phi (residual variance) is Inverse Gamma (.5*alpha,.5*lambda)
double f2int_imom(double phi) {
  double ans, *inputphi= f2int_pars.phi;
  f2int_pars.phi= &phi;
  ans= pimomMarginalKC(f2int_pars.sel,f2int_pars.nsel,&f2int_pars) * dinvgammaC(phi,.5*(*f2int_pars.alpha),.5*(*f2int_pars.lambda));
  f2int_pars.phi= inputphi;
  //ans= pimomMarginalKC(f2int_pars.sel,f2int_pars.nsel,f2int_pars.n,f2int_pars.p,f2int_pars.y,f2int_pars.sumy2,f2int_pars.XtX,f2int_pars.ytX,&phi,f2int_pars.tau,&f2int_pars.method,f2int_pars.B,&logscale) * dinvgammaC(phi,.5*(*f2int_pars.alpha),.5*(*f2int_pars.lambda));
  return(ans);
}

SEXP pimomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda) {
  int *sel=INTEGER(Ssel), *nsel=INTEGER(Snsel), *n=INTEGER(Sn), *p=INTEGER(Sp), *method=INTEGER(Smethod), *B=INTEGER(SB), *logscale=INTEGER(Slogscale), r=1;
  double *y=REAL(Sy), *sumy2=REAL(Ssumy2), *x=REAL(Sx), *XtX=REAL(SXtX), *ytX=REAL(SytX), *tau=REAL(Stau), *alpha=REAL(Salpha), *lambda=REAL(Slambda), *rans, emptydouble=0;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars,n,p,y,sumy2,x,XtX,ytX,method,B,alpha,lambda,&emptydouble,tau,&r,&emptydouble,&emptydouble,logscale);
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pimomMarginalUC(sel, nsel, &pars);
  UNPROTECT(1);
  return ans;
}


double pimomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  int i, j, zero=0;
  double ans, er, sumer2, **V, **Vinv, *thest, ypred, phiest, intmc, intlapl, adj, *inputphi;
  //double *alpha, double *lambda
  set_f2int_pars((*pars).XtX,(*pars).ytX,(*pars).tau,(*pars).n,(*pars).p,sel,nsel,(*pars).y,(*pars).sumy2,(*pars).method,(*pars).B,(*pars).alpha,(*pars).lambda,&zero);
  if ((*(*pars).method)==2) {
    V= dmatrix(1,*nsel,1,*nsel); Vinv= dmatrix(1,*nsel,1,*nsel);
    thest= dvector(1,*nsel);

    addct2XtX((*pars).tau,(*pars).XtX,sel,nsel,(*pars).p,V); //add tau to diagonal elem of XtX
    inv_posdef_upper(V,*nsel,Vinv);
    Asym_xsel(Vinv,*nsel,(*pars).ytX,sel,thest);
    for (i=0, sumer2=0; i<(*(*pars).n); i++) {
      for (j=1, ypred=0; j<=(*nsel); j++) { ypred += (*pars).x[i + (*(*pars).n)*sel[j-1]] * thest[j]; }
      er= (*pars).y[i] - ypred;
      sumer2+= er*er;
    }
    phiest= (sumer2 + (*(*pars).lambda))/(*(*pars).alpha + *(*pars).n);
    inputphi= (*pars).phi; (*pars).phi= &phiest; 
    (*(*pars).method)= 1; //IS evaluation of marginal for phi=phiest
    intmc= pimomMarginalKC(sel, nsel, pars); 
    (*(*pars).method)= 0; //Laplace approx for phi=phiest
    intlapl= pimomMarginalKC(sel, nsel, pars); 
    (*pars).phi= inputphi; (*(*pars).method)= 2;  //reset input values for phi, method
    if (intlapl==0) { intmc+= 1.0e-300; intlapl+= 1.0e-300; } //avoid numerical zero
    adj= intmc/intlapl;
    f2int_pars.method= &zero;  //set method to eval marginal for known phi to Laplace approx

    free_dmatrix(V,1,*nsel,1,*nsel); free_dmatrix(Vinv,1,*nsel,1,*nsel);
    free_dvector(thest,1,*nsel);
  } else {
    adj= 1.0;
  }

  ans= adj * (qromo(f2int_imom,0.0,100,midpnt) + qromo(f2int_imom,100,1.0e30,midinf));
  if ((*(*pars).logscale)==1) ans= log(ans);
  return(ans);
}

/*
double pimomMarginalUC(int *sel, int *nsel, int *n, int *p, double *y, double *sumy2, double *x, double *XtX, double *ytX, double *tau, int *method, int *B, int *logscale, double *alpha, double *lambda) {
  int i, j, zero=0, one=1;
  double ans, er, sumer2, **V, **Vinv, *thest, ypred, phiest, intmc, intlapl, adj;

  set_f2int_pars(XtX,ytX,tau,n,p,sel,nsel,y,sumy2,method,B,alpha,lambda);
  if ((*method)==2) {
    V= dmatrix(1,*nsel,1,*nsel); Vinv= dmatrix(1,*nsel,1,*nsel);
    thest= dvector(1,*nsel);

    addct2XtX(tau,XtX,sel,nsel,p,V); //add tau to diagonal elem of XtX
    inv_posdef_upper(V,*nsel,Vinv);
    Asym_xsel(Vinv,*nsel,ytX,sel,thest);
    for (i=0, sumer2=0; i<(*n); i++) {
      for (j=1, ypred=0; j<=(*nsel); j++) { ypred += x[i + (*n)*sel[j-1]] * thest[j]; }
      er= y[i] - ypred;
      sumer2+= er*er;
    }
    phiest= (sumer2 + (*lambda))/(*alpha + *n);
    intmc= pimomMarginalKC(sel,nsel,n,p,y,sumy2,XtX,ytX,&phiest,tau,&one,B,&zero); //IS evaluation of marginal for phi=phiest
    intlapl= pimomMarginalKC(sel,nsel,n,p,y,sumy2,XtX,ytX,&phiest,tau,&zero,B,&zero); //Laplace approx for phi=phiest
    if (intlapl==0) { intmc+= 1.0e-300; intlapl+= 1.0e-300; } //avoid numerical zero
    adj= intmc/intlapl;
    f2int_pars.method= &zero;  //set method to eval marginal for known phi to Laplace approx

    free_dmatrix(V,1,*nsel,1,*nsel); free_dmatrix(Vinv,1,*nsel,1,*nsel);
    free_dvector(thest,1,*nsel);
  } else {
    adj= 1.0;
  }

  ans= adj * (qromo(f2int_imom,0.0,100,midpnt) + qromo(f2int_imom,100,1.0e30,midinf));
  if ((*logscale)==1) ans= log(ans);
  return(ans);
}
*/



//*************************************************************************************
// Product eMOM routines
//*************************************************************************************

double pemomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  return 0.0;
}

double pemomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  return 0.0;
}
