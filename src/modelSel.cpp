#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "cstat.h"
#include "modelSel.h"
#include "do_mombf.h"
#include "modselIntegrals.h"
#include <vector>
#include "Polynomial.h"

//Global variables defined for minimization/integration routines
struct marginalPars f2opt_pars, f2int_pars;


//*************************************************************************************
//SETTING PRIOR & MARGINALS
//*************************************************************************************

pt2margFun set_marginalFunction(int *prCoef, int *knownphi, int *family) {
  //Returns pointer to function to compute the marginal density of the data for a given model indicator
  // - prCoef: 0 for product MOM, 1 for product iMOM, 2 for product eMOM, 3 for Zellner's prior
  // - knownphi: 1 if residual variance phi is known, 0 otherwise. knownphi==1 currently only allowed for Normal residuals
  // - family: distribution of residuals. 1 for Normal, 2 for two-piece Normal; 3 for Laplace,; 4 for two-piece Laplace
  // Note: if phi known, when actually calling the returned pt2margFun, phi must be set in the parameter of type struct marginalPars *
  pt2margFun ans=NULL;
  if ((*family)==1) {  //Normal residuals
    if (*prCoef==0) {
      if (*knownphi==1) { ans= pmomMarginalKC; } else { ans= pmomMarginalUC; }
    } else if (*prCoef==1) {
      if (*knownphi==1) { ans= pimomMarginalKC; } else { ans= pimomMarginalUC; }
    } else if (*prCoef==2) {
      if (*knownphi==1) { ans= pemomMarginalKC; } else { ans= pemomMarginalUC; }
    } else if (*prCoef==3) {
      if (*knownphi==1) { ans= zellnerMarginalKC; } else { ans= zellnerMarginalUC; }
    }
  } else if ((*family)==2) { //Two-piece Normal residuals
    if (*prCoef==0) {
      ans= pmomMargSkewNormU;
    } else if (*prCoef==1) {
      ans= pimomMargSkewNormU;
    } else if (*prCoef==2) {
      ans= pemomMargSkewNormU;
    } else if (*prCoef==3) {
      Rprintf("Zellner prior with two-piece Normal residuals not currently implemented");
    }
  } else {
    Rf_error("This error distribution is not available");
  }
  return ans;
}


pt2margFun set_priorFunction(int *prDelta) {
  //Returns pointer to function to compute the prior probability of a model indicator
  // - prDelta: 0 for uniform, 1 for binomial, 2 for beta-binomial
  pt2margFun ans=NULL;
  if (*prDelta==0) { ans= unifPrior; } else if (*prDelta==1) { ans= binomPrior; } else if (*prDelta==2) { ans= betabinPrior; }
  return ans;
}

pt2modavgPrior set_priorFunction_modavg(int *priorModel) {
  //Returns pointer to function to compute the prior probability of a model indicator
  // - priorModel: 0 for uniform, 1 for binomial, 2 for beta-binomial
  pt2modavgPrior ans=NULL;
  if (*priorModel==0) { ans= unifPrior_modavg; } else if (*priorModel==1) { ans= binomPrior_modavg; } else if (*priorModel==2) { ans= betabinPrior_modavg; }
  return ans;
}


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

//multiply symmetric A[1..ncolA][1..ncolA] (formatted as vector) with x[1..nsel].
//Use only selected elems in A and all elems in x[1..nsel]
void Asel_x(double *A, int ncolA, double *x, int nsel, int *sel, double *ans) {
  int _i, _j;
  for (_i=1;_i<= nsel;_i++) {
    for (_j=1, ans[_i]=0; _j<= nsel; _j++) { ans[_i]+= A[sel[_j]*ncolA+sel[_i]] * x[_j]; }
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


//*************************************************************************************
// MODEL AVERAGING ROUTINES
//*************************************************************************************

void set_modavgPars(struct modavgPars *pars, int *n, int *p1, int *p2, int *isbinary, int *ybinary, double *y, double *sumy2, double *x1, double *x2, double *XtX, double *ytX, double *cholS2, double *S2inv, double *cholS2inv, double *colsumx1sq, double *alpha, double *lambda, int *priorCoef, int *r, double *tau1, double *tau2, int *priorTau1, double *atau1, double *btau1, int *priorModel, double *prModelpar) {
  (*pars).n= n;
  (*pars).p1= p1;
  (*pars).p2= p2;
  (*pars).isbinary= isbinary;
  (*pars).ybinary= ybinary;
  (*pars).y= y;
  (*pars).sumy2= sumy2;
  (*pars).x1= x1;
  (*pars).x2= x2;
  (*pars).XtX= XtX;
  (*pars).ytX= ytX;
  (*pars).cholS2= cholS2;
  (*pars).S2inv= S2inv;
  (*pars).cholS2inv= cholS2inv;
  (*pars).colsumx1sq= colsumx1sq;
  (*pars).alpha= alpha;
  (*pars).lambda= lambda;
  (*pars).priorCoef= priorCoef;
  (*pars).r= r;
  (*pars).tau1= tau1;
  (*pars).tau2= tau2;
  (*pars).priorTau1= priorTau1;
  (*pars).atau1= atau1;
  (*pars).btau1= btau1;
  (*pars).priorModel= priorModel;
  (*pars).prModelpar= prModelpar;
}


//MH within Gibbs scheme to sample from the joint posterior of (theta,delta) under a product MOM prior and a linear regression model
//Input
// - pars: data, pre-computed quantities and prior parameters. See definition of struct modavgPars for details
// - niter: number of MCMC iterations
// - thinning: only 1 out of each thinning iterations is kept
// - burnin: burning
// - niniModel: number of variables initially in the model
// - iniModel: initial model (vector of length p1)
// - iniCoef1: initial coefficient values for variables being selected (length p1)
// - iniCoef2: initial coefficient values for variables always in the model (length p2)
// - iniPhi: initial residual variance value
// - iniOthers: initial values for other variables. Currently this indicates initial tau value, and is only use if (*pars).priorTau1 !=0.
// - verbose: set verbose==1 to print iteration progress every 10% of the iterations
//Output
// - postModel: MCMC saves for variable inclusion indicators
// - margpp: marginal posterior prob for inclusion of each covariate (uses MH acceptance prob, which is more precise than simply averaging inclusion indicators)
// - postCoef1: MCMC saves for regression coefficients of variables being selected
// - postCoef2: MCMC saves for regression coefficients of variables which are always in the model
// - postPhi: MCMC saves for residual variance
// - postOther: MCMC saves for other parameters. Currently saves tau values (pMOM prior precision) when (*pars).priorTau1 != 0.

SEXP pmomLM_I(SEXP postModel, SEXP margpp, SEXP postCoef1, SEXP postCoef2, SEXP postPhi, SEXP postOther, SEXP niter, SEXP thinning, SEXP burnin, SEXP niniModel, SEXP iniModel, SEXP iniCoef1, SEXP iniCoef2, SEXP iniPhi, SEXP iniOthers, SEXP verbose, SEXP n, SEXP p1, SEXP p2, SEXP isbinary, SEXP ybinary, SEXP y, SEXP sumy2, SEXP x1, SEXP x2, SEXP XtX, SEXP ytX, SEXP cholS2, SEXP S2inv, SEXP cholS2inv, SEXP colsumx1sq, SEXP alpha, SEXP lambda, SEXP priorCoef, SEXP r, SEXP tau1, SEXP tau2, SEXP priorTau1, SEXP atau1, SEXP btau1, SEXP priorModel, SEXP prModelpar) {
  struct modavgPars pars;
  SEXP ans;

  set_modavgPars(&pars,INTEGER(n),INTEGER(p1),INTEGER(p2),INTEGER(isbinary),INTEGER(ybinary),REAL(y),REAL(sumy2),REAL(x1),REAL(x2),REAL(XtX),REAL(ytX),REAL(cholS2),REAL(S2inv),REAL(cholS2inv),REAL(colsumx1sq),REAL(alpha),REAL(lambda),INTEGER(priorCoef),INTEGER(r),REAL(tau1),REAL(tau2),INTEGER(priorTau1),REAL(atau1),REAL(btau1),INTEGER(priorModel),REAL(prModelpar));
  pmomLM(INTEGER(postModel), REAL(margpp), REAL(postCoef1), REAL(postCoef2), REAL(postPhi), REAL(postOther), &pars, INTEGER(niter), INTEGER(thinning), INTEGER(burnin), INTEGER(niniModel), INTEGER(iniModel), REAL(iniCoef1), REAL(iniCoef2), REAL(iniPhi), REAL(iniOthers), INTEGER(verbose));

  PROTECT(ans = allocVector(REALSXP, 1));
  *REAL(ans)= 1.0;
  UNPROTECT(1);
  return ans;
}


void pmomLM(int *postModel, double *margpp, double *postCoef1, double *postCoef2, double *postPhi, double *postOther, struct modavgPars *pars, int *niter, int *thinning, int *burnin, int *niniModel, int *iniModel, double *iniCoef1, double *iniCoef2, double *iniPhi, double *iniOthers, int *verbose) {
  int i, j, k, ilow, iupper, savecnt, niterthin, niter10, nsel= *niniModel, *curModel, newdelta, n=*(*pars).n, p1=*(*pars).p1, p2=*(*pars).p2, psn, resupdate, isbinary=*(*pars).isbinary;
  double *res, *partialres, sumres2, sumpartialres2, newcoef, *curCoef1, *curCoef2, curPhi=1.0, *linpred1, *linpred2, pinclude, *temp;
  if (*verbose) Rprintf("Running MCMC");
  niterthin= (int) floor((*niter - *burnin +.0)/(*thinning +.0));
  if (*niter >10) { niter10= *niter/10; } else { niter10= 1; }
  if (*burnin >0) { ilow= - *burnin; savecnt=0; iupper= *niter - *burnin +1; } else { ilow=0; savecnt=1; iupper= *niter; }
  //Initialize
  curModel= ivector(0,p1); curCoef1= dvector(0,p1); curCoef2= dvector(0,p2); linpred1= dvector(0,n); linpred2= dvector(0,n);
  res= dvector(0, n); partialres= dvector(0,n);
  for (i=0; i<p1; i++) { margpp[i]= 0; curModel[i]= postModel[i*niterthin]= iniModel[i]; curCoef1[i]= postCoef1[i*niterthin]= iniCoef1[i]; }
  for (i=0; i<p2; i++) { curCoef2[i]= postCoef2[i*niterthin]= iniCoef2[i]; }
  if (isbinary) { curPhi= postPhi[0]= 1.0; } else { curPhi= postPhi[0]= *iniPhi; }
  postOther[0]= iniOthers[0];
  Avecx((*pars).x1, curCoef1, linpred1, 0, n-1, 0, p1-1);
  Avecx((*pars).x2, curCoef2, linpred2, 0, n-1, 0, p2-1);
  if (isbinary) sample_latentProbit((*pars).y,res,&sumres2,(*pars).ybinary,linpred1,linpred2,pars);
  for (i=0, sumres2=0; i<n; i++) { res[i]= partialres[i]= (*pars).y[i] - linpred1[i] - linpred2[i]; sumres2+= res[i]*res[i]; }
  sumpartialres2= sumres2;
  //MCMC iterations
  for (i=ilow; i< iupper; i++) {
    //Sample (curCoef1,curModel)
    for (j=0; j< *(*pars).p1; j++) {
      if (curModel[j]) {
	for (k=0, sumpartialres2=0; k<n; k++) { partialres[k]= res[k] + curCoef1[j] * ((*pars).x1[n*j+k]); sumpartialres2+= partialres[k]*partialres[k]; }
      }
      MHTheta1pmom(&newdelta, &newcoef, &pinclude, &resupdate, res, partialres, &sumres2, &sumpartialres2, j, &nsel, curModel, curCoef1, &curPhi, pars);
      if (newdelta > curModel[j]) { nsel++; } else if (newdelta < curModel[j]) { nsel--; }
      curModel[j]= newdelta; curCoef1[j]= newcoef;
      if (i>=0) margpp[j]+= pinclude;
      if (resupdate) { temp= partialres; partialres= res; res=temp; }
    }
    //Sample curCoef2
    for (k=0; k<n; k++) res[k]+= linpred2[k];
    simTheta2(curCoef2, res, &curPhi, pars);
    Avecx((*pars).x2, curCoef2, linpred2, 0, n-1, 0, p2-1);
    for (k=0, sumres2=0; k<n; k++) { res[k]-= linpred2[k]; sumres2+= res[k]*res[k]; }
    //Sample phi
    if (isbinary==0) { curPhi= simPhipmom(&nsel, curModel, curCoef1, curCoef2, &sumres2, pars); }
    //Sample tau
    if (*(*pars).priorTau1 !=0) { *(*pars).tau1= simTaupmom(&nsel, curModel, curCoef1, &curPhi, pars); }
    //Sample latent variables (only for probit model)
    if (isbinary) {
      Avecx((*pars).x1, curCoef1, linpred1, 0, n-1, 0, p1-1);  //update linpred1 (linpred2 already updated)
      sample_latentProbit((*pars).y,res,&sumres2,(*pars).ybinary,linpred1,linpred2,pars);
    }
    //Save values
    if ((i>0) && (i % (*thinning))==0) {
      for (j=0; j<p1; j++) {
        psn= niterthin*j+savecnt;
        postModel[psn]= curModel[j];
        postCoef1[psn]= curCoef1[j];
      }
      for (j=0; j<p2; j++) postCoef2[niterthin*j+savecnt]= curCoef2[j];
      postPhi[savecnt]= curPhi;
      if (*(*pars).priorTau1 !=0) postOther[savecnt]= *(*pars).tau1;
      savecnt++;
    }
    if ((*verbose ==1) && (i%niter10)==0) Rprintf(".");
  } //end MCMC for
  if (iupper>ilow) { for (j=0; j< p1; j++) { margpp[j] /= (iupper-imax_xy(0,ilow)+.0); } } //from sum to average
  if (*verbose ==1) Rprintf("Done.\n");
  free_ivector(curModel,0,p1); free_dvector(curCoef1,0,p1); free_dvector(curCoef2,0,p2); free_dvector(linpred1,0,n); free_dvector(linpred2,0,n);
  free_dvector(res,0,n); free_dvector(partialres,0,n);
}


//Sample from the posterior of latent variables in probit model given the regression coefficients (i.e. the linear predictor)
//Input:
// - ybinary: response variable (1: success; 0: failure)
// - linpred1: linear predictor for current regression coefficients associated to variables under selection
// - linpred2: linear predictor associated to adjustment variables
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output:
// - y: sampled values for the latent variables
// - res: sampled residuals, i.e. y - linpred1 - linpred2
// - sumres2: sum(res^2)
// - (*pars).ytX: updated t(y) %*% x1
// - (*pars).sumy2: updated sum(y^2)
void sample_latentProbit(double *y, double *res, double *sumres2, int *ybinary, double *linpred1, double *linpred2, struct modavgPars *pars) {
  int i;
  double linpred, plinpred, u;
  for (i=0, *sumres2=0, *(*pars).sumy2=0; i< *(*pars).n; i++) {
    linpred= linpred1[i] + linpred2[i];
    plinpred= pnormC(-linpred,0,1);
    if (ybinary[i]) {
      u= plinpred + (1.0-plinpred) * runif();  //u ~ Unif(plinpred,1)
    } else {
      u= plinpred * runif(); //u ~ Unif(0,plinpred)
    }
    res[i]= qnormC(u,0,1);
    (*sumres2)+= res[i]*res[i];
    y[i]= linpred + res[i];
    (*(*pars).sumy2)+= y[i]*y[i];
  }
  Atvecx((*pars).x1,y,(*pars).ytX,0,*(*pars).p1 -1,0,*(*pars).n -1); //update ytX=Xty
}

//Univariate MH update of (coef,delta)
//Input:
// - res: vector with residuals. Only used if curModel[j]==0
// - partialres: vector with partial residuals removing the variable from the model. Only used if curModel[j]==1
// - sumres2: sum(res*res)
// - sumpartialres2: sum(partialres*partialres)
// - j: index of the variable for which move is to be proposed
// - curModel: vector of 0's and 1's indicating which variables are currently in the model
// - curCoef1: current values of the regression coefficients for variables undergoing selection
// - curPhi: current value of the residual variance
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output:
// - newdelta: new value
// - newcoef: new coefficient
// - pinclude: probability of including the variable in the model
// - res: on input, vector with residuals given the current coefficient value. On output, residuals given the updated coefficient value.
void MHTheta1pmom(int *newdelta, double *newcoef, double *pinclude, int *resupdate, double *res, double *partialres, double *sumres2, double *sumpartialres2, int j, int *nsel, int *curModel, double *curCoef1, double *curPhi, struct modavgPars *pars) {
  int n= *(*pars).n, logscale=1, nsel0, nsel1, deltaprop, nu, i;
  double m1, *xj, m0, logbf, logpratio, thetaprop, m, S, propPars[6], lhood, lprior, lprop, lambda=0.0, num, den, sqrtPhi=sqrt(*curPhi);
  pt2modavgPrior priorFunction= NULL;
  *resupdate= 0;
  xj= (*pars).x1+j*n; //pointer to variable j in x1
  priorFunction= set_priorFunction_modavg((*pars).priorModel);
  //Propose delta
  if (curModel[j]) {
    m1= pmomMargKuniv(partialres, xj, sumpartialres2, (*pars).colsumx1sq+j, &n, curPhi, (*pars).tau1, (*pars).r, &logscale);
    m0= dnormC_jvec(partialres, *(*pars).n, 0, sqrtPhi, 1);
    nsel0= *nsel -1; nsel1= *nsel;
  } else {
    m1= pmomMargKuniv(res, xj, sumres2, (*pars).colsumx1sq+j, &n, curPhi, (*pars).tau1, (*pars).r, &logscale);
    m0= dnormC_jvec(res, *(*pars).n, 0, sqrtPhi, 1);
    nsel0= *nsel; nsel1= *nsel +1;
  }
  logbf= m0-m1;
  logpratio= priorFunction(curModel, &nsel0, pars) - priorFunction(curModel, &nsel1, pars); //we use curModel in both cases as priorFunction currently only depends on nb vars
  *pinclude= 1.0/(1.0+exp(logbf+logpratio));
  if (runif() < *pinclude) { deltaprop=1; } else { deltaprop=0; }
  //Propose coef
  nu= (int) sqrt((double) n);
  if ((curModel[j]==0) && (deltaprop==0)) {  //proposal is to keep variable out of the model
    *newdelta=0; *newcoef=0;
  } else {
    S= (*pars).colsumx1sq[j] + 1.0/(*(*pars).tau1);
    if (curModel[j]) {
      for (i=0, m=0; i<n; i++) m+= xj[i]*partialres[i];
      m= m/S;
      proposalpmom(propPars, &m, &S, curPhi, (*pars).r, (*pars).tau1, &n, partialres, xj, &m1, &nu);
    } else {
      for (i=0, m=0; i<n; i++) m+= xj[i]*res[i];
      m= m/S;
      proposalpmom(propPars, &m, &S, curPhi, (*pars).r, (*pars).tau1, &n, res, xj, &m1, &nu);
    }
    if (curModel[j] && deltaprop) {  //proposal is to keep variable in the model
      thetaprop= rtmixC(propPars, propPars+2, propPars+4, nu, 2);
      for (i=0, lhood=0; i<n; i++) {
        partialres[i]-= thetaprop*xj[i];
        lhood+= dnormC(partialres[i],0,sqrtPhi,1) - dnormC(res[i],0,sqrtPhi,1);
      }
      lprior= dmom(thetaprop,0,*(*pars).tau1,*curPhi,*(*pars).r,1) - dmom(curCoef1[j],0,*(*pars).tau1,*curPhi,*(*pars).r,1);
      lprop= dtmixC(curCoef1[j],propPars,propPars+2,propPars+4,nu,2,1) - dtmixC(thetaprop,propPars,propPars+2,propPars+4,nu,2,1);
      lambda= exp(lhood+lprior+lprop);
    } else if ((curModel[j]==0) && deltaprop) { //proposal is to add variable to the model
      thetaprop= rtmixC(propPars, propPars+2, propPars+4, nu, 2);
      for (i=0, num=0; i<n; i++) {
        partialres[i]= res[i] - thetaprop*xj[i];
        num+= dnormC(partialres[i],0,sqrtPhi,1);
      }
      num+= dmom(thetaprop,0,*(*pars).tau1,*curPhi,*(*pars).r,1);
      den= dtmixC(thetaprop,propPars,propPars+2,propPars+4,nu,2,1) + m1;
      lambda= exp(num-den);
    } else {    //(curModel[j] && (deltaprop==0)), i.e. proposal is to drop variable from the model
      thetaprop=0;
      num= dtmixC(curCoef1[j],propPars,propPars+2,propPars+4,nu,2,1) + m1;
      for (i=0, den=0; i<n; i++) { den+= dnormC(res[i],0,sqrtPhi,1); }
      den+= dmom(curCoef1[j],0,*(*pars).tau1,*curPhi,*(*pars).r,1);
      lambda= exp(num-den);
    }
    if (runif()<lambda) {
      *newdelta=deltaprop; *newcoef= thetaprop;
      *resupdate= 1;  //signal that res and partialres have to be interchanged after exiting the function
      for (i=0, *sumres2=0; i< n; i++) (*sumres2)+= partialres[i]*partialres[i];
    } else {
      *newdelta=curModel[j]; *newcoef= curCoef1[j];
    }
  }
}

//Find parameters for univariate pmom proposal distribution (2 component mixture of T distributions)
//   Posterior: N(e; xj*theta; phi*I) * pmom(theta; phi, r, tau) / m1
//   Mixture: w1 * T_nu(theta;mu1,sigma21) + (1-w1) * T_nu(theta;mu2,sigma22)  (sigma21, sigma22 denote variances)
//Input:
// - m, S: posterior location & scale parameters
// - phi: residual variance
// - r: product MOM power parameter
// - tau1: product MOM prior dispersion parameter
// - e: response variable
// - xj: predictor
// - m1: normalization constant
// - nu: desired degrees of freedom
//Output: means in propPars[0:1], SD in propPars[2:3], weights in propPars[4:5]
void proposalpmom(double *propPars, double *m, double *S, double *phi, int *r, double *tau1, int *n, double *e, double *xj, double *m1, int *nu) {
  int i;
  double eps, fmode, sqrtPhi=sqrt(*phi), temp, temp2, doubler=2*(*r), ct2;
  //Find modes
  eps= sqrt((*m)*(*m) + 8.0*(*r)*(*phi)/(*S));
  propPars[0]= .5*(*m - eps); propPars[1]= .5*(*m + eps);
  //Find density at the mode
  for (i=0, fmode=0; i< *n; i++) fmode+= dnormC(e[i],propPars[1]*xj[i],sqrtPhi,1);
  fmode+= dmom(propPars[1],0,*tau1,*phi,*r,1) - *m1;
  fmode= exp(fmode);
  //Proposal variances
  temp= (*S)/(*phi);
  propPars[2]= sqrt(1.0/(temp + doubler/(propPars[0]*propPars[0]))); propPars[3]= sqrt(1.0/(temp + doubler/(propPars[1]*propPars[1])));
  temp2= .5*(*nu); temp= temp2+.5;
  //Weights
  ct2= exp(gamln(&temp) - .5*log((double) (*nu)) - gamln(&temp2) - .5*log(M_PI*propPars[3]*propPars[3]));
  propPars[4]= max_xy(0,(fmode-ct2)/(dnormC(propPars[1],propPars[0],propPars[2],0) - ct2));
  propPars[5]= 1-propPars[4];
}

//Expectation of prod_j th_j^{2*power} under multivariate Normal/T with mean m, covariance S, dimension n, degrees of freedom dof (dof=-1 for Normal)
SEXP eprod_I(SEXP m, SEXP S, SEXP n, SEXP power, SEXP dof) {
  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, 1));
  //int i,j,nn=INTEGER(n)[0]; double **sigma;
  //sigma= dmatrix(1,nn,1,nn);
  //for (i=1; i<=nn; i++) {
  //  for (j=1; j<i; j++) { sigma[i][j]= sigma[j][i]= REAL(S)[(j-1)*nn+(i-1)]; }
  //  sigma[i][i]= REAL(S)[(i-1)*nn+(i-1)];
  //}
  //*REAL(ans)= mvtexpect(REAL(m)-1, sigma, nn, INTEGER(power)[0], REAL(dof)[0]);
  //free_dmatrix(sigma, 1,nn,1,nn);
  *REAL(ans)= mvtexpect_vec(REAL(m), REAL(S), INTEGER(n)[0], INTEGER(power)[0], REAL(dof)[0]);
  UNPROTECT(1);
  return ans;
}



//Univariate marginal density under a product MOM prior (known variance case)
// integral N(y; x*theta, phi*I) * (theta^2/(tau*phi))^r * N(theta; 0; tau*phi) / (2r-1)!! d theta
// - y: response variable (must be a vector)
// - x: design matrix (must be a vector)
// - sumy2: sum(y*y)
// - n: length of y
// - phi: residual variance
// - tau: prior variance parameter
// - logscale: if set to 1 the log of the integral is returned
double pmomMargKuniv(double *y, double *x, double *sumy2, double *sumxsq, int *n, double *phi, double *tau, int *r, int *logscale) {
  int i; double ans, m, s, I, doubler=2.0*(*r);
  s= *sumxsq + 1.0/(*tau);
  for (i=0, m=0; i< *n; i++) { m+= y[i]*x[i]; }
  m/= s;
  I= log(mnorm(doubler,m,sqrt(*phi/s)));
  ans= I -.5*(*sumy2 - s*m*m)/(*phi) - .5*(*n)*log(2*M_PI*(*phi)) - .5*(log(s)+log(*tau)) - ldoublefact(doubler-1) - (*r)*log((*tau)*(*phi));
  if (*logscale ==0) ans=exp(ans);
  return(ans);
}

//Sample from conditional posterior of coefficients for variables not undergoing selection
//Input:
// - partialres: partial residuals obtained by removing adjustment variables from the model
// - phi: current value of the residual variance
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output:
// - theta2: sample from conditional posterior of theta2 given the data and all other parameters
void simTheta2(double *theta2, double *partialres, double *phi, struct modavgPars *pars) {
  int i, j; double *tmp, *m, **cholS, sqrtPhi=sqrt(*phi);
  //m= S2inv * t(x2) * partialres
  tmp= dvector(0,*(*pars).p2); m= dvector(0,*(*pars).p2); cholS= dmatrix(1,*(*pars).p2,1,*(*pars).p2);
  Atvecx((*pars).x2, partialres, tmp, 0, *(*pars).p2 -1, 0, *(*pars).n -1);
  Avecx((*pars).S2inv, tmp, m, 0, *(*pars).p2, 0, *(*pars).p2);
  //S= S2inv * phi
  for (i=0; i< *(*pars).p2; i++) { for (j=0; j< *(*pars).p2; j++) { cholS[i+1][j+1]= sqrtPhi * (*pars).cholS2inv[i+j*(*(*pars).p2)]; } }
  //Generate theta2 ~ N(m,S)
  rmvnormC(theta2-1,*(*pars).p2,m-1,cholS);
  free_dvector(tmp,0,*(*pars).p2); free_dvector(m,0,*(*pars).p2); free_dmatrix(cholS,1,*(*pars).p2,1,*(*pars).p2);
}


//Sample from conditional posterior of residual variance phi under a product MOM prior
//Input:
// - curModel: vector of 0's and 1's indicating which variables are currently in the model
// - curCoef1: current values of the regression coefficients for variables undergoing selection
// - curCoef2: current values of the regression coefficients for variables not undergoing selection
// - ssr: residual sum of squares for current coefficient values
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output: random draw from conditional posterior of phi given the data and all other parameteres
double simPhipmom(int *nsel, int *curModel, double *curCoef1, double *curCoef2, double *ssr, struct modavgPars *pars) {
  int i; double a, b, sumth1, sumth2;
  a= *(*pars).alpha + *(*pars).n + (2*(*(*pars).r)+1)*(*nsel) + *(*pars).p2;
  for (i=0, sumth1=0; i< *(*pars).p1; i++) { if (curModel[i]==1) sumth1+= curCoef1[i]*curCoef1[i]; }
  for (i=0, sumth2=0; i< *(*pars).p2; i++) { sumth2+= curCoef2[i]*curCoef2[i]; }
  b= *(*pars).lambda + sumth1/(*(*pars).tau1) + sumth2/(*(*pars).tau2) + *ssr;
  return(1.0/rgammaC(.5*a,.5*b));
}

//Sample from conditional posterior of tau given all other parameters
//Input:
// - nsel: number of variables undergoing selection currently in the model
// - curModel: vector of 0's and 1's indicating which variables are currently in the model
// - curCoef1: current values of the regression coefficients for variables undergoing selection
// - curPhi: current value for the residual variance
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output: random draw from conditional posterior of tau given the data and all other parameters
double simTaupmom(int *nsel, int *curModel, double *curCoef1, double *curPhi, struct modavgPars *pars) {
  int i; double a, b, sumth1;
  a= *(*pars).atau1 + (2*(*(*pars).r)+1)*(*nsel);
  for (i=0, sumth1=0; i< *(*pars).p1; i++) { if (curModel[i]==1) sumth1+= curCoef1[i]*curCoef1[i]; }
  b= *(*pars).btau1 + sumth1/(*curPhi);
  return(1.0/rgammaC(.5*a,.5*b));
}


//********************************************************************************************
// GENERAL MARGINAL DENSITY CALCULATION ROUTINES
//********************************************************************************************

void set_marginalPars(struct marginalPars *pars, int *n,int *p,double *y,double *sumy2,double *x,double *XtX,double *ytX,int *method,int *B,double *alpha,double *lambda,double *phi,double *tau,double *taualpha, int *r,double *prDeltap,double *parprDeltap, int *logscale, double *offset) {
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
  (*pars).taualpha= taualpha;
  (*pars).r= r;
  (*pars).prDeltap= prDeltap;
  (*pars).parprDeltap= parprDeltap;
  (*pars).logscale= logscale;
  (*pars).offset= offset;
}

void set_f2opt_pars(double *m, double **S, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *phi, double *tau, int *r, int *n, int *p, int *sel, int *nsel) {
  f2opt_pars.m= m;
  f2opt_pars.S= S;
  f2opt_pars.sumy2= sumy2;
  f2opt_pars.XtX= XtX;
  f2opt_pars.ytX= ytX;
  f2opt_pars.alpha= alpha;
  f2opt_pars.lambda= lambda;
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

//modelSelectionGibbs: Gibbs sampler for model selection in linear regression for several choices of prior distribution
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

SEXP modelSelectionCI(SEXP SpostSample, SEXP SpostOther, SEXP Smargpp, SEXP SpostMode, SEXP SpostModeProb, SEXP SpostProb, SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP Sniter, SEXP Sthinning, SEXP Sburnin, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staualpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose) {
  int logscale=1;
  double offset=0;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars, INTEGER(Sn), INTEGER(Sp), REAL(Sy), REAL(Ssumy2), REAL(Sx), REAL(SXtX), REAL(SytX), INTEGER(Smethod), INTEGER(SB), REAL(Salpha),REAL(Slambda), REAL(Sphi), REAL(Stau), REAL(Staualpha), INTEGER(Sr), REAL(SprDeltap), REAL(SparprDeltap), &logscale, &offset);
  modelSelectionGibbs2(INTEGER(SpostSample), REAL(SpostOther), REAL(Smargpp), INTEGER(SpostMode), REAL(SpostModeProb), REAL(SpostProb), INTEGER(Sknownphi), INTEGER(Sfamily), INTEGER(SpriorCoef), INTEGER(SpriorDelta), INTEGER(Sniter), INTEGER(Sthinning), INTEGER(Sburnin), INTEGER(Sndeltaini), INTEGER(Sdeltaini), INTEGER(Sverbose), &pars);

  PROTECT(ans = allocVector(REALSXP, 1));
  *REAL(ans)= 1.0;
  UNPROTECT(1);
  return ans;
}



void modelSelectionGibbs(int *postSample, double *postOther, double *margpp, int *postMode, double *postModeProb, double *postProb, int *knownphi, int *family, int *prCoef, int *prDelta, int *niter, int *thinning, int *burnin, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars) {
  int i, j, k, *sel, *selnew, *selaux, nsel, nselnew, niter10, niterthin, savecnt, ilow, iupper;
  double currentJ, newM, newP, newJ=0.0, ppnew, u;
  pt2margFun marginalFunction=NULL, priorFunction=NULL; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);

  marginalFunction= set_marginalFunction(prCoef, knownphi, family);
  priorFunction= set_priorFunction(prDelta);

  sel= ivector(0,*(*pars).p); selnew= ivector(0,*(*pars).p);

  //Initialize
  if (*verbose ==1) Rprintf("Running Gibbs sampler");
  niterthin= (int) floor((*niter - *burnin +.0)/(*thinning+.0));
  if (*niter >10) { niter10= *niter/10; } else { niter10= 1; }
  for (j=0; j< *(*pars).p; j++) { margpp[j]= 0; }
  nsel= *ndeltaini;
  for (j=0; j< nsel; j++) { sel[j]= deltaini[j]; postMode[deltaini[j]]= 1; }
  if ((*prDelta)==2) { postOther[0]= *(*pars).prDeltap; }
  currentJ= marginalFunction(sel,&nsel,pars) + priorFunction(sel,&nsel,pars);
  postProb[0]= *postModeProb= currentJ;
  if (*burnin >0) { ilow=-(*burnin); savecnt=0; iupper= *niter - *burnin +1; } else { ilow=1; savecnt=1; iupper= *niter; } //if no burnin, start at i==1 & save initial value

  //Iterate
  for (i=ilow; i< iupper; i++) {
    for (j=0; j< *(*pars).p; j++) {
      sel2selnew(j,sel,&nsel,selnew,&nselnew); //copy sel into selnew, adding/removing jth element
      if (nselnew <= (*(*pars).p)) {
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
      } else {
 	ppnew= -0.01;
      }
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
  if (*verbose==1) Rprintf(" Done.\n");

  free_ivector(sel,0,*(*pars).p); free_ivector(selnew,0,*(*pars).p);
}

//Same as modelSelectionGibbs, but log(integrated likelihood) + log(prior) are managed through an object of class modselIntegrals
void modelSelectionGibbs2(int *postSample, double *postOther, double *margpp, int *postMode, double *postModeProb, double *postProb, int *knownphi, int *family, int *prCoef, int *prDelta, int *niter, int *thinning, int *burnin, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars) {
  int i, j, k, *sel, *selnew, *selaux, nsel, nselnew, niter10, niterthin, savecnt, ilow, iupper;
  double currentJ, newJ=0.0, ppnew, u;
  pt2margFun marginalFunction=NULL, priorFunction=NULL; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);

  marginalFunction= set_marginalFunction(prCoef, knownphi, family);
  priorFunction= set_priorFunction(prDelta);

  modselIntegrals *integrals= new modselIntegrals(marginalFunction, priorFunction, *(*pars).p);

  sel= ivector(0,*(*pars).p); selnew= ivector(0,*(*pars).p);

  //Initialize
  if (*verbose ==1) Rprintf("Running Gibbs sampler");
  niterthin= (int) floor((*niter - *burnin +.0)/(*thinning+.0));
  if (*niter >10) { niter10= *niter/10; } else { niter10= 1; }
  for (j=0; j< *(*pars).p; j++) { margpp[j]= 0; }
  nsel= *ndeltaini;
  for (j=0; j< nsel; j++) { sel[j]= deltaini[j]; postMode[deltaini[j]]= 1; }
  if ((*prDelta)==2) { postOther[0]= *(*pars).prDeltap; }
  currentJ= integrals->getJoint(sel,&nsel,pars);
  postProb[0]= *postModeProb= currentJ;
  if (*burnin >0) { ilow=-(*burnin); savecnt=0; iupper= *niter - *burnin +1; } else { ilow=1; savecnt=1; iupper= *niter; } //if no burnin, start at i==1 & save initial value

  //Iterate
  for (i=ilow; i< iupper; i++) {
    for (j=0; j< *(*pars).p; j++) {
      sel2selnew(j,sel,&nsel,selnew,&nselnew); //copy sel into selnew, adding/removing jth element
      if (nselnew <= (*(*pars).p)) {
        newJ= integrals->getJoint(selnew,&nselnew,pars);
        if (newJ > *postModeProb) {   //update posterior mode
          *postModeProb= newJ;
          for (k=0; k< *(*pars).p; k++) { postMode[k]= 0; }
          for (k=0; k< nselnew; k++) { postMode[selnew[k]]= 1; }
        }
        ppnew= 1.0/(1.0+exp(currentJ-newJ));
        if (i>=0) { if (nselnew>nsel) { margpp[j]+= ppnew; } else { margpp[j]+= (1-ppnew); } }
      } else {
	ppnew= -0.01;
      }
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
  if (*verbose==1) Rprintf(" Done.\n");

  free_ivector(sel,0,*(*pars).p); free_ivector(selnew,0,*(*pars).p);
  delete integrals;
}

//greedyVarSelC: greedy search for posterior mode in variable selection.
//               Similar to Gibbs sampling, except that deterministic updates are made iff there is an increase in post model prob
//               The scheme proceeds until no variable is included/excluded or niter iterations are reached
// Input arguments: same as in modelSelectionC.
SEXP greedyVarSelCI(SEXP SpostMode, SEXP SpostModeProb, SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP Sniter, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staualpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose) {
  int logscale=1;
  double offset=0;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars, INTEGER(Sn), INTEGER(Sp), REAL(Sy), REAL(Ssumy2), REAL(Sx), REAL(SXtX), REAL(SytX), INTEGER(Smethod), INTEGER(SB), REAL(Salpha),REAL(Slambda), REAL(Sphi), REAL(Stau), REAL(Staualpha), INTEGER(Sr), REAL(SprDeltap), REAL(SparprDeltap), &logscale, &offset);
  greedyVarSelC(INTEGER(SpostMode),REAL(SpostModeProb),INTEGER(Sknownphi),INTEGER(Sfamily),INTEGER(SpriorCoef),INTEGER(SpriorDelta),INTEGER(Sniter),INTEGER(Sndeltaini),INTEGER(Sdeltaini),INTEGER(Sverbose),&pars);
  PROTECT(ans = allocVector(REALSXP, 1));
  *REAL(ans)= 1.0;
  UNPROTECT(1);
  return ans;

}

void greedyVarSelC(int *postMode, double *postModeProb, int *knownphi, int *family, int *prCoef, int *prDelta, int *niter, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars) {
  int i, j, *sel, *selnew, *selaux, nsel, nselnew, nchanges;
  double newJ;
  pt2margFun marginalFunction=NULL, priorFunction=NULL; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);

  marginalFunction= set_marginalFunction(prCoef, knownphi, family);
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
        selaux= sel; sel=selnew; selnew=selaux; nsel=nselnew; //update model indicator
        nchanges++;
      }
    } //end j for
  }
  if (*verbose==1) Rprintf("Done.\n");

  free_ivector(sel,0,*(*pars).p); free_ivector(selnew,0,*(*pars).p);
}



void sel2selnew(int newelem,int *sel,int *nsel,int *selnew,int *nselnew) {
//Copy sel into selnew.
// - If j in sel, don't copy it in selnew and set nselnew=nsel-1.
// - If j not in sel, add it to selnew and set nselnew=nsel+1.
  int i, ii, found;
  for (i=0, found=0; (i< *nsel) && (found==0); i++) { selnew[i]= sel[i]; found= (newelem==sel[i]); }
  if (found==0) { //add newelem
    selnew[i]= newelem; (*nselnew)= (*nsel)+1;
  } else {  //remove new elem
    for (ii=i; ii< *nsel; ii++) { selnew[ii-1]= sel[ii]; }
    (*nselnew)= (*nsel)-1;
  }
}


//********************************************************************************************
// PRIORS ON MODEL SPACE (always return on log scale)
//********************************************************************************************

double unifPrior(int *sel, int *nsel, struct marginalPars *pars) { return 0.0; }
double unifPrior_modavg(int *sel, int *nsel, struct modavgPars *pars) { return 0.0; }

//nsel ~ Binom(p,prDeltap)
double binomPrior(int *sel, int *nsel, struct marginalPars *pars) {
  return dbinomial(*nsel,*(*pars).p,*(*pars).prDeltap,1);
}
double binomPrior_modavg(int *sel, int *nsel, struct modavgPars *pars) {
  return dbinomial(*nsel,*(*pars).p1,(*pars).prModelpar[0],1);
}

//nsel ~ Beta-Binomial(prModelpar[0],prModelPar[1])
double betabinPrior(int *sel, int *nsel, struct marginalPars *pars) {
  return bbPrior(*nsel,*(*pars).p,(*pars).parprDeltap[0],(*pars).parprDeltap[1],1);
}
double betabinPrior_modavg(int *sel, int *nsel, struct modavgPars *pars) {
  return bbPrior(*nsel,*(*pars).p1,(*pars).prModelpar[0],(*pars).prModelpar[1],1);
}




//*************************************************************************************
// LEAST SQUARES
//*************************************************************************************


void leastsquares(double *theta, double *phi, double *ypred, double *y, double *x, double *XtX, double *ytX, int *n, int *sel, int *nsel) {
  //Least squares estimate for y= x[,sel] %*% theta + e, where e ~ N(0,phi)   (phi is the variance)
  //Input
  // - y: observed response
  // - x: predictors in vector format
  // - n: length(y)
  // - sel: variables in x to be included in the model are given by sel[0], sel[1], etc.
  // - nsel: length(sel)
  //Output
  // - theta: least squares estimate for theta
  // - phi: MLE for residual variance (i.e. SSR/n). If <1.0e-10 then phi=1.0e-10 is returned
  // - ypred: predicted y, i.e. x[,sel] %*% theta where theta is the least squares estimate
  int i;
  double zero=0, **S, **Sinv, detS, e;

  S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
  (*phi)= 0;

  if ((*nsel)>0) {
    //Least squares
    addct2XtX(&zero,XtX,sel,nsel,nsel,S);
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,ytX,sel,theta);

    //MLE for residual variance
    Aselvecx(x, theta+1, ypred, 0, (*n) -1, sel, nsel);
    for (i=0; i< (*n); i++) { e= y[i]-ypred[i]; (*phi) += e*e; }

  } else {

    for (i=0; i< (*n); i++) { (*phi) += y[i]*y[i]; }

  }

  (*phi)= (*phi) / (*n);
  if ((*phi) < 1.0e-10) { (*phi)= 1.0e-10; }

  free_dmatrix(S, 1,*nsel,1,*nsel); free_dmatrix(Sinv, 1,*nsel,1,*nsel);
}



//*************************************************************************************
// TWO-PIECE NORMAL ROUTINES
//*************************************************************************************

double pmomMargSkewNormU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=1;
  return nlpMargSkewNorm(sel, nsel, pars, &prior);
}

double pimomMargSkewNormU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=2;
  return nlpMargSkewNorm(sel, nsel, pars, &prior);
}

double pemomMargSkewNormU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=3;
  return nlpMargSkewNorm(sel, nsel, pars, &prior);
}


double nlpMargSkewNorm(int *sel, int *nsel, struct marginalPars *pars, int *prior) {
//Integrated likelihood for linear regression model y= Xtheta + e where e ~ two-piece Normal(0,vartheta,alpha) (vartheta prop to variance, alpha gives asymmetry)
// and priors theta ~ pMOM/piMOM/peMOM(0,tau*vartheta), vartheta ~ IG(alpha/2,lambda/2), atanh(alpha) ~ pMOM(0,taualpha)
// Input
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - pars: parameters needed to compute marginal likelihood
// - prior: prior==1 for pMOM, prior==2 for piMOM, prior==3 for peMOM
// Output: integrated likelihood

  bool initmle=true, posdef;
  int maxit= 20, p= (*nsel)+2, n= (*((*pars).n));
  double ans, *thmode, fmode, **hess, **cholhess, d, det, *ypred;

  thmode= dvector(1,p); hess= dmatrix(1, p, 1, p); ypred=dvector(0,n-1);

  postmodeSkewNorm(thmode, &fmode, hess, sel, nsel, (*pars).n, (*pars).y, (*pars).x, (*pars).XtX, (*pars).ytX, &maxit, (*pars).tau, (*pars).taualpha, (*pars).alpha, (*pars).lambda, &initmle, prior);

  if ((*(*pars).method ==0) | (*(*pars).method ==0)) { //Laplace or MC

    cholhess= dmatrix(1,p,1,p);
    d= p + 2.0;
    choldc(hess,p,cholhess,&posdef);
    if (!posdef) {
      int i;
      double lmin=0, *vals;
      vals= dvector(1,p);
      eigenvals(hess,p,vals);
      for (i=1; i<=p; i++) if (vals[i]<lmin) lmin= vals[i];
      lmin = -lmin + .01;
      for (i=1; i<=p; i++) hess[i][i] += lmin;
      choldc(hess,p,cholhess,&posdef);
      free_dvector(vals,1,p);
    }
    det= choldc_det(cholhess, p);

    if (*(*pars).method ==0) { //Laplace

      ans= -fmode + 0.5 * d * LOG_M_2PI - 0.5*log(det);

    } else if (*(*pars).method ==0) { //Monte Carlo

      int i, j, nu=3;
      double *thsim, **cholV, **cholVinv, ctnu= (nu+2.0)/(nu+.0), detVinv, term1, term2;

      thsim= dvector(1, p); cholV= dmatrix(1,p,1,p); cholVinv= dmatrix(1,p,1,p);

      thmode[p+1]= log(thmode[p+1]); thmode[p+2]= atanh(thmode[p+2]);
      cholS_inv(cholhess, p, cholV);
      for (i=1; i<=p; i++) {
	for (j=i; j<=p; j++) {
	  cholV[i][j]= sqrt(ctnu) * cholV[i][j];
	  cholVinv[i][j]= sqrt(ctnu) * hess[i][j];
	}
      }
      detVinv= exp(log(det) + p * log(ctnu));

      ans= 0;
      for (i=1; i<= (*(*pars).B); i++) {
	rmvtC(thsim, p, thmode, cholV, nu);
	fnegSkewnorm(&term1,ypred,thsim,sel,nsel,(*pars).n,(*pars).y,(*pars).x,(*pars).XtX,(*pars).tau,(*pars).taualpha,(*pars).alpha,(*pars).lambda,prior,true);
        term2= -dmvtC(thsim, p, thmode, cholVinv, detVinv, nu, 1);
	ans += exp(-term1 + fmode + term2);
      }
      ans= log(ans / ((*(*pars).B)+.0)) - fmode;

      free_dvector(thsim, 1,p); free_dmatrix(cholV, 1,p,1,p); free_dmatrix(cholVinv, 1,p,1,p);
    }

    free_dmatrix(cholhess, 1,p,1,p);

  } else {
    Rf_error("Method must be 'Laplace' or 'MC'");
  }
  if (*((*pars).logscale) == 0) ans= exp(ans);

  free_dvector(thmode, 0,p); free_dmatrix(hess, 1,p,1,p); free_dvector(ypred,0,n-1);
  return(ans);

}



void postmodeSkewNorm(double *thmode, double *fmode, double **hess, int *sel, int *nsel, int *n, double *y, double *x, double *XtX, double *ytX, int *maxit, double *tau, double *taualpha, double *alpha, double *lambda, bool *initmle, int *prior) {
//Posterior mode for two-piece normal under pMOM, piMOM or peMOM prior on (theta,atanh(alpha)) and vartheta ~ IG(alpha/2,lambda/2)
  // Input
  // - y: observed response
  // - x: design matrix
  // - maxit: maximum number of iterations
  // - tau: dispersion parameter for MOM prior on theta
  // - tau.alpha: dispersion parameter for MOM prior on atanh(alpha)
  // - alpha, lambda: prior on vartheta ~ IG(alpha/2,lambda/2)
  // - init: init=='mle' to initialize at MLE, else initialize at least squares
  // Ouput
  // - thmode: posterior mode for (theta,vartheta,alpha) (i.e. in original parameterization)
  // - fmode: minus log-joint evaluated at thmode
  // - hess: hessian evaluated at thmode

  bool posdef;
  int i, ii, j, p=(*nsel)+2;
  double err, damp, *g, **H, **Hinv, *delta, lmin=0, *vals, fnew, *thnew, *ypred;

  ypred= dvector(0,*n -1);

  if (*initmle) {  //Initialize at MLE

    mleSkewnorm(thmode, y, x, maxit);

  } else {  //Initialize at least-squares for theta; set (phi,alpha)= argmax likelihood for given theta

    double s1=0, s2=0, pows1, pows2;

    leastsquares(thmode, thmode+p+1, ypred, y, x, XtX, ytX, n, sel, nsel);

    for (i=0; i<(*n); i++) {
      if (y[i]<=ypred[i]) { s1+= pow(y[i]-ypred[i], 2.0); } else { s2+= pow(y[i]-ypred[i], 2.0); }
    }

    pows1= pow(s1, 1.0/3.0); pows2= pow(s2, 1.0/3.0);
    thmode[p+2]= (pows1 - pows2)/(pows1 + pows2);  //estimate for alpha
    thmode[p+1]= (0.25/((*n)+.0)) * pow(pows1 + pows2, 3.0);  //estimate for phi

  }

  thmode[p-1]= log(thmode[p-1]); //phi
  thmode[p]= atanh(thmode[p]);   //alpha (Note: atanh(z)= 0.5*(log(1+z)-log(1-z)))

  g= dvector(1,p); delta= dvector(1,p); thnew= dvector(1,p);
  H= dmatrix(1,p,1,p); Hinv= dmatrix(1,p,1,p);

  i=1; err=1; damp=0.01;

  fnegSkewnorm(fmode,ypred,thmode,sel,nsel,n,y,x,XtX,tau,taualpha,alpha,lambda,prior,true);

  while ((err>0.0001) & (i<(*maxit))) {

    fpnegSkewnorm(g,thmode,ypred,sel,nsel,n,y,x,tau,taualpha,alpha,lambda,prior); //gradient
    fppnegSkewnorm(H,thmode,ypred,sel,nsel,n,y,x,tau,taualpha,alpha,lambda,prior); //Hessian

    choldc_inv(H,p,Hinv,&posdef);

    if (posdef) {
      Ax(Hinv,g,delta,1,p,1,p);
    } else {
      //Ensure H is posdef
      vals= dvector(1,p);
      eigenvals(H,p,vals);
      for (j=1; j<=p; j++) if (vals[j]<lmin) lmin= vals[j];
      lmin = -lmin + .01;
      for (j=1; j<=p; i++) H[j][j] += lmin;
      choldc_inv(H,p,Hinv,&posdef);
      Ax(Hinv,g,delta,1,p,1,p);
      free_dvector(vals,1,p);
    }

    for (j=1; j<=p; j++) { thnew[j]= thmode[j] - delta[j]; }
    fnegSkewnorm(&fnew,ypred,thnew,sel,nsel,n,y,x,XtX,tau,taualpha,alpha,lambda,prior,true);

    //If Newton update fails, use Levenberg-Marquardt (LMA)
    ii= 1;
    while ((fnew > (*fmode)) & (ii<5)) {
      for (j=1; j<=p; j++) H[j][j] *= (1.0+damp);
      choldc_inv(H,p,Hinv,&posdef);
      Ax(Hinv,g,delta,1,p,1,p);
      for (j=1; j<=p; j++) { thnew[j]= thmode[j] - delta[j]; }
      fnegSkewnorm(&fnew,ypred,thnew,sel,nsel,n,y,x,XtX,tau,taualpha,alpha,lambda,prior,true);
      ii++;
    }

    //If new value improves target function, update thmode, fmode
    if (fnew<(*fmode)) {
      err= 0;
      for (j=1; j<=p; j++) {
	err= max_xy(err,fabs(delta[j]));
	thmode[j]= thnew[j];
      }
      (*fmode)= fnew;
      i++;
    } else {
      i= (*maxit);
    }
  }

  thmode[p-1]= exp(thmode[p-1]);
  thmode[p]= tanh(thmode[p]); //Note: tanh(z)= -1 + 2/(1+exp(-2*z))

  for (i=1; i<=p; i++) {
    hess[i][i]= H[i][i];
    for (j=1; j<i; j++) { hess[i][j]= hess[j][i]= H[i][j]; }
  }

  free_dvector(ypred, 0,*n -1); free_dvector(g,1,p); free_dvector(delta,1,p); free_dvector(thnew,1,p);
  free_dmatrix(H,1,p,1,p); free_dmatrix(Hinv,1,p,1,p);
}




void fnegSkewnorm(double *ans, double *ypred, double *th, int *sel, int *nsel, int *n, double *y, double *x, double *XtX, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, bool logscale) {
//Negative log-joint for two-piece Normal under MOM/eMOM/iMOM prior on coef and IG on variance
// Input
// - th[1..nsel+2]: (theta, log(vartheta), atanh(alpha)) where theta=regression coef, vartheta \propto variance and alpha=asymmetry parameter in [-1,1]
// - Other parameters as in postmodeSkewNorm
// Output: value of minus the log-joint evaluated at th
  double scale, alpha;

  scale= exp(th[*nsel +1]);
  alpha= tanh(th[*nsel +2]);
  loglSkewnorm(ans, ypred, th, nsel, sel, n, &scale, &alpha, y, x, XtX);
  (*ans)= -(*ans);

  if ((*prior)==1) {

    if ((*nsel)>0) {
      (*ans) += -dmomvec(th+1,*nsel,0.0,*tau,scale,1,1) - dmom(alpha,0.0,*taualpha,1.0,1,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    } else {
      (*ans) += -dmom(alpha,0.0,*taualpha,1.0,1,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    }

  } else if ((*prior)==2) {

    if ((*nsel)>0) {
      (*ans) += -dimomvec(th+1,*nsel,0.0,*tau,scale,1) - dimom(alpha,0.0,*taualpha,1.0,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    } else {
      (*ans) += -dimom(alpha,0.0,*taualpha,1.0,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    }

  } else if ((*prior)==3) {

    if ((*nsel)>0) {
      (*ans) += -demomvec(th+1,*nsel,*tau,scale,1) - demom(alpha,*taualpha,1.0,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    } else {
      (*ans) += -demom(alpha,*taualpha,1.0,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    }

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

  if (!logscale) { (*ans)= exp(*ans); }
}



void fpnegSkewnorm(double *g, double *th, double *ypred, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior) {
 //Gradient of fnegSkewnorm
  int i, one=1, nselplus1= (*nsel)+1;
  double *gprior, zero=0;

  gprior= dvector(1,(*nsel)+2);

  loglnegGradSkewNorm(g,th,nsel,sel,n,y,ypred,x);

  if ((*prior)==1) {

    dmomiggrad(gprior+1,&nselplus1,th+1,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) { g[i] -= gprior[i]; }

    dmomgrad(gprior+(*nsel)+2,&one,th+(*nsel)+2,&zero,taualpha);
    g[(*nsel)+2] -= gprior[(*nsel)+2];

  } else if ((*prior)==2) {

    dimomiggrad(gprior+1,&nselplus1,th+1,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) { g[i] -= gprior[i]; }

    dimomgrad(gprior+(*nsel)+2,&one,th+(*nsel)+2,&zero,taualpha);
    g[(*nsel)+2] -= gprior[(*nsel)+2];

  } else if ((*prior)==3) {

    dimomiggrad(gprior+1,&nselplus1,th+1,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) { g[i] -= gprior[i]; }

    dimomgrad(gprior+(*nsel)+2,&one,th+(*nsel)+2,&zero,taualpha);
    g[(*nsel)+2] -= gprior[(*nsel)+2];

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

  free_dvector(gprior,1,(*nsel)+2);
}


void fppnegSkewnorm(double **H, double *th, double *ypred, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior) {
 //Hessian of fnegSkewnorm

  int i, j, one=1, nselplus1= (*nsel)+1;
  double **Hprior, hprioralpha, zero=0;

  Hprior= dmatrix(1,nselplus1,1,nselplus1);

  loglnegHessSkewNorm(H,th,nsel,sel,n,y,ypred,x);

  if ((*prior)==1) {

    dmomighess(Hprior,&nselplus1,th+1,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) {
      H[i][i] -= Hprior[i][i];
      for (j=1; j<i; j++) {
	H[i][j]= H[j][i]= H[i][j] - Hprior[i][j];
      }
    }

    dmomhess(&hprioralpha,&one,th+(*nsel)+2,&zero,taualpha);
    H[(*nsel)+2][(*nsel)+2] -= hprioralpha;

  } else if ((*prior)==2) {

    dimomighess(Hprior,&nselplus1,th+1,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) {
      H[i][i] -= Hprior[i][i];
      for (j=1; j<i; j++) {
	H[i][j]= H[j][i]= H[i][j] - Hprior[i][j];
      }
    }

    dimomhess(&hprioralpha,&one,th+(*nsel)+2,&zero,taualpha);
    H[(*nsel)+2][(*nsel)+2] -= hprioralpha;

  } else if ((*prior)==3) {

    demomighess(Hprior,&nselplus1,th+1,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) {
      H[i][i] -= Hprior[i][i];
      for (j=1; j<i; j++) {
	H[i][j]= H[j][i]= H[i][j] - Hprior[i][j];
      }
    }

    demomhess(&hprioralpha,&one,th+(*nsel)+2,&zero,taualpha);
    H[(*nsel)+2][(*nsel)+2] -= hprioralpha;

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

  free_dmatrix(Hprior,1,nselplus1,1,nselplus1);
}



void loglSkewnorm(double *ans, double *ypred, double *th, int *nsel, int *sel, int *n, double *scale, double *alpha, double *y, double *x, double *XtX) {
  //Log-likelihood function of a linear model with two-piece normal errors evaluated at th=(theta,scale,alpha)
  int i;
  double w1, w2;

  w1= 0.5 / (pow(1.0 + (*alpha),2) * (*scale));
  w2= 0.5 / (pow(1.0 - (*alpha),2) * (*scale));
  (*ans)= -0.5*(*n)*(LOG_M_2PI + log(*scale));

  if ((*nsel)>0) {

    ypred= dvector(0,*n -1);
    Aselvecx(x, th+1, ypred, 0, (*n) -1, sel, nsel); //ypred= x %*% th

    for (i=0; i<(*n); i++) {

      if (y[i]<ypred[i]) { (*ans) -= w1 * pow(y[i]-ypred[i],2); } else { (*ans) -= w2 * pow(y[i]-ypred[i],2); }

    }

    free_dvector(ypred,0,*n -1);

  } else {

    for (i=0; i<(*n); i++) {

      if (y[i]<0) { (*ans) -= w1 * pow(y[i],2); } else { (*ans) -= w2 * pow(y[i],2); }

    }

  }

}


void loglnegGradSkewNorm(double *g, double *th, int *nsel, int *sel, int *n, double *y, double *ypred, double *x) {
  //Gradient of minus the log-likelihood function of a linear model with two-piece Normal errors
  int i;
  double sigma, alpha, alphat, w1, w2, ws1, ws2, *y0, *Wy0, y0Wy0=0, y0Wsy0=0;

  Wy0= dvector(0,*n -1);
  sigma= exp(th[*nsel +1]);
  alphat= tanh(th[*nsel +2]);
  alpha= th[*nsel+2];

  w1= 1.0 / (pow(1.0 + alphat,2));
  w2= 1.0 / (pow(1.0 - alphat,2));
  ws1= -2.0 / (pow(cosh(alpha),2) * pow(1.0+alphat,3));
  ws2=  2.0 / (pow(cosh(alpha),2) * pow(1.0-alphat,3));

  if ((*nsel)>0) {

    y0= dvector(0,*n -1);

    for (i=0; i<(*n); i++) {
      y0[i]= (y[i]-ypred[i]);

      if (y[i]<ypred[i]) {
	Wy0[i]= w1 * y0[i];
	y0Wsy0+= pow(y0[i],2)*ws1;
      } else {
	Wy0[i]= w2 * y0[i];
	y0Wsy0+= pow(y0[i],2)*ws2;
      }
      y0Wy0+= y0[i] * Wy0[i];

    }

    Atselvecx(x, Wy0, g+1, 0, (*n)-1, sel, nsel);  //g[1:nsel]= t(x) %*% Wy0
    for (i=1; i<=(*nsel); i++) g[i]= -g[i]/sigma;
    free_dvector(y0,0,*n -1);

  } else {

    for (i=0; i<(*n); i++) {
      if (y[i]<0) {
	Wy0[i]= w1 * y[i];
	y0Wsy0+= pow(y[i],2)*ws1;
      } else {
	Wy0[i]= w2 * y[i];
	y0Wsy0+= pow(y[i],2)*ws2;
      }
      y0Wy0+= y[i] * Wy0[i];
    }

  }

  g[*nsel +1]= 0.5*(*n) - 0.5*y0Wy0/sigma;

  g[*nsel +2]= 0.5*y0Wsy0/sigma;

  free_dvector(Wy0,0,*n -1);
}



void loglnegHessSkewNorm(double **H, double *th, int *nsel, int *sel, int *n, double *y, double *ypred, double *x) {
  //Hessian of minus the log-likelihood function of a linear model with two-piece Normal errors
  //NOTE: log.lik.hess function in R, but we need to change the sign of the output
  int i, j, k, idxi, idxj;
  double sigma, alphat, alpha, w, w1, w2, ws1, ws2, wss1, wss2, *y0, *Wy0, *Wsy0, y0Wy0=0, y0Wsy0=0, y0Wssy0=0;

  Wy0= dvector(0,*n -1); Wsy0= dvector(0,*n -1);
  sigma= exp(th[*nsel +1]);
  alphat= tanh(th[*nsel +2]);
  alpha= th[*nsel+2];

  w1= 1.0 / (pow(1.0 + alphat,2));
  w2= 1.0 / (pow(1.0 - alphat,2));
  ws1= -2.0 / (pow(cosh(alpha),2) * pow(1.0+alphat,3));
  ws2=  2.0 / (pow(cosh(alpha),2) * pow(1.0-alphat,3));
  wss1= 2.0 * exp(-2.0*alpha) + 4.0*exp(-4.0*alpha);
  wss2= 2.0 * exp(2.0*alpha) + 4.0*exp(4.0*alpha);

  if ((*nsel)>0) {

    y0= dvector(0,*n -1);

    for (i=0; i<(*n); i++) {
      y0[i]= (y[i]-ypred[i]);

      if (y[i]<ypred[i]) {
	Wy0[i]= w1 * y0[i];
	Wsy0[i]= ws1 * y0[i];
	y0Wsy0+= pow(y0[i],2)*ws1;
	y0Wssy0+= pow(y0[i],2)*wss1;
      } else {
	Wy0[i]= w2 * y0[i];
	Wsy0[i]= ws2 * y0[i];
	y0Wsy0+= pow(y0[i],2)*ws2;
	y0Wssy0+= pow(y0[i],2)*wss2;
      }
      y0Wy0+= y0[i] * Wy0[i];

    }

    free_dvector(y0,0,*n -1);

    //Compute H[1:*nsel,1:*nsel] <- t(X0)%*%W%*%X0/sigma
    for (i=1; i<=(*nsel); i++) {
      idxi= (*n)*i;
      for (j=i; j<=(*nsel); j++) {
	idxj= (*n)*j;
	H[i][j]= 0;
	for (k=0; k<(*n); k++) {
	  if (y[i]<ypred[i]) { w= w1; } else { w= w2; }
	  H[i][j] += x[k+ idxi] * x[k +idxj] * w;  //x[k][i] * x[k][j] * w
	}
      }
      for (j=1; j<i; j++) { H[i][j]= H[j][i]; }
    }
    for (i=1; i<=(*nsel); i++) { for (j=1; j<=(*nsel); j++) { H[i][j]= H[i][j]/sigma; H[j][i]= H[i][j]; } }

    //Compute H[1:*nsel,*nsel+1] <- t(X0)%*%Wy0/sigma
    Atselvecx(x, Wy0, (double *)H +(*nsel)*(*nsel), 0, (*n)-1, sel, nsel);
    for (i=1; i<=(*nsel); i++) { H[i][*nsel +1]= H[i][*nsel +1]/sigma; H[*nsel +1][i]= H[i][*nsel +1]; }

    //Compute H[1:*nsel,*nsel+2] <- -t(X0)%*%Wsy0/sigma
    Atselvecx(x, Wsy0, (double *)H +(*nsel)*(*nsel +1), 0, (*n)-1, sel, nsel);
    for (i=1; i<=(*nsel); i++) { H[i][*nsel +2]= -H[i][*nsel +2]/sigma; H[*nsel +2][i]= H[i][*nsel +2]; }

  } else {

    for (i=0; i<(*n); i++) {
      if (y[i]<0) {
	Wy0[i]= w1 * y[i];
	Wsy0[i]= ws1 * y[i];
	y0Wsy0+= pow(y[i],2)*ws1;
	y0Wssy0+= pow(y[i],2)*wss1;
      } else {
	Wy0[i]= w2 * y[i];
	Wsy0[i]= ws2 * y[i];
	y0Wsy0+= pow(y[i],2)*ws2;
	y0Wssy0+= pow(y[i],2)*wss2;
      }
      y0Wy0+= y[i] * Wy0[i];
    }

    //Compute H[1:*nsel,*nsel+1] <- t(X0)%*%Wy0/sigma
    Atselvecx(x, Wy0, (double *)H +(*nsel)*(*nsel), 0, (*n)-1, sel, nsel);
    for (i=1; i<=(*nsel); i++) { H[i][*nsel +1]= H[i][*nsel +1]/sigma; H[*nsel +1][i]= H[i][*nsel +1]; }

    //Compute H[1:*nsel,*nsel+2] <- -t(X0)%*%Wsy0/sigma
    Atselvecx(x, Wsy0, (double *)H +(*nsel)*(*nsel +1), 0, (*n)-1, sel, nsel);
    for (i=1; i<=(*nsel); i++) { H[i][*nsel +2]= -H[i][*nsel +2]/sigma; H[*nsel +2][i]= H[i][*nsel +2]; }

  }

  H[*nsel +1][*nsel +1]= 0.5 * y0Wy0 / sigma;
  H[*nsel +2][*nsel +2]= 0.5 * y0Wssy0 / sigma;
  H[*nsel +1][*nsel +2]= H[*nsel +2][*nsel +1] = -0.5 * y0Wsy0 / sigma;

  free_dvector(Wy0, 0,*n -1); free_dvector(Wsy0, 0,*n -1);
}


//TO DO: mleSkewnorm NEEDS TO BE ADAPTED FROM R CODE
void mleSkewnorm(double *thmode, double *y, double *x, int *maxit) {
  thmode[1]= 0.0;
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
  int i, emptyint, iter, maxit=100;
  double emptydouble, ftol= 1.0e-5, **Vopt, detVopt, **dirth;

  Vopt= dmatrix(1,*nsel,1,*nsel); dirth= dmatrix(1,*nsel,1,*nsel);
  set_f2opt_pars(m,S,&emptydouble,&emptydouble,&emptydouble,&emptydouble,&emptydouble,phi,tau,r,n,nsel,&emptyint,nsel);

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
  bool posdef;
  int i;
  double **cholSinv, *thsim, ans, normfac;

  thsim= dvector(1,*nsel);
  cholSinv= dmatrix(1,*nsel,1,*nsel);
  choldc(Sinv,*nsel,cholSinv,&posdef); //compute cholesky decomposition
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
  bool posdef;
  int i;
  double **cholSinv, *thsim, ans, normfac;

  thsim= dvector(1,*nsel);
  cholSinv= dmatrix(1,*nsel,1,*nsel);
  choldc(Sinv,*nsel,cholSinv,&posdef); //compute cholesky decomposition
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
// - method==0 for Laplace; method==1 for Monte Carlo; method==2 for plug-in; method== -1 for exact calculation if p<20
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale

SEXP pmomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale) {
  struct marginalPars pars;
  double *rans, emptydouble=0, offset=0, *taualpha=NULL;
  SEXP ans;

  set_marginalPars(&pars,INTEGER(Sn),INTEGER(Sp),REAL(Sy),REAL(Ssumy2),&emptydouble,REAL(SXtX),REAL(SytX),INTEGER(Smethod),INTEGER(SB),&emptydouble,&emptydouble,REAL(Sphi),REAL(Stau),taualpha,INTEGER(Sr),&emptydouble,&emptydouble,INTEGER(Slogscale),&offset);
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
  double *m, s, **S, **Sinv, detS, num, den, logtau= log(*(*pars).tau), tauinv= 1.0/(*(*pars).tau), logphi= log(*(*pars).phi), ans=0.0, *thopt, **Voptinv, fopt;

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

    if ((*(*pars).method ==0) | ((*(*pars).method == -1) & ((*nsel)>10)))  { //Laplace

      thopt= dvector(1,*nsel); Voptinv= dmatrix(1,*nsel,1,*nsel);
      momIntegralApproxC(&ans,thopt,Voptinv,&fopt,(*pars).n,nsel,m,S,&detS,(*pars).phi,(*pars).tau,(*pars).r,(*pars).logscale);
      free_dvector(thopt,1,*nsel); free_dmatrix(Voptinv,1,*nsel,1,*nsel);

    } else if (*(*pars).method ==1) { //MC

      for (i=1; i<=(*nsel); i++) { Sinv[i][i]= (*(*pars).phi)*Sinv[i][i]; for (j=i+1; j<=(*nsel); j++) { Sinv[i][j]=Sinv[j][i]= (*(*pars).phi)*Sinv[i][j]; } }
      ans= MC_mom_normal(m,Sinv,(*pars).r,nsel,(*pars).B);

    } else if (*(*pars).method ==2) { //Plug-in

      ans= rsumlogsq(m,(*pars).r,nsel);

    } else if ((*(*pars).method == -1) & ((*nsel)<=10)) { //Exact

      Voptinv= dmatrix(1,*nsel,1,*nsel);
      for (i=1; i<= *nsel; i++) for (j=i; j<= *nsel; j++) Voptinv[i][j]= Voptinv[j][i]= Sinv[i][j] * (*(*pars).phi);
      ans= log(mvtexpect(m, Voptinv, *nsel, 2, -1));
      free_dmatrix(Voptinv,1,*nsel,1,*nsel);

    }
    ans+= num - den;
    free_dvector(m,1,*nsel);
    free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}


SEXP pmomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda) {
  double *rans, emptydouble=0, offset=0, *taualpha=NULL;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars,INTEGER(Sn),INTEGER(Sp),REAL(Sy),REAL(Ssumy2),REAL(Sx),REAL(SXtX),REAL(SytX),INTEGER(Smethod),INTEGER(SB),REAL(Salpha),REAL(Slambda),&emptydouble,REAL(Stau),taualpha,INTEGER(Sr),&emptydouble,&emptydouble,INTEGER(Slogscale),&offset);
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pmomMarginalUC(INTEGER(Ssel), INTEGER(Snsel), &pars);
  UNPROTECT(1);
  return ans;
}


double pmomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  int i, j, nu;
  double num, den, ans=0.0, term1, *m, **S, **Sinv, detS, *thopt, **Voptinv, fopt, phiadj, tauinv= 1.0/(*(*pars).tau), nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda), ss;
  if (*nsel ==0) {
    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*(LOG_M_PI) + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);
  } else {
    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&tauinv,(*pars).XtX,sel,nsel,(*pars).p,S);
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    nuhalf= (*(*pars).r)*(*nsel) + .5*(*(*pars).n + *(*pars).alpha);
    nu= (int) (2.0*nuhalf);

    ss= *(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel);
    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2.0) - log(ss));
    den= (*nsel)*ldoublefact(2*(*(*pars).r)-1.0) + .5*(*(*pars).n * LOG_M_2PI + log(detS)) + (*nsel)*(.5 + *(*pars).r)*log(*(*pars).tau) + gamln(&alphahalf);

    if ((*(*pars).method ==0) | ((*(*pars).method == -1) & ((*nsel)>10)))  { //Laplace

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

    } else if ((*(*pars).method == -1) & ((*nsel)<=10)) { //Exact

      Voptinv= dmatrix(1,*nsel,1,*nsel);
      for (i=1; i<= *nsel; i++) for (j=i; j<= *nsel; j++) Voptinv[i][j]= Voptinv[j][i]= Sinv[i][j] * ss / (nu+.0);
      ans= log(mvtexpect(m, Voptinv, *nsel, 2, nu));
      free_dmatrix(Voptinv,1,*nsel,1,*nsel);

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


//Hessian of fimomNegC
// - ans: hessian matrix evaluated at th (indexed at 1, i.e. ans[1:(*nsel)][1:(*nsel)])
// - th: th[1:(*nsel)] indicates point at which to evaluate the hessian.
// - Other arguments as in fimomNegC_non0
void fppimomNegC_non0(double **ans, double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i, j;
  double th2;

  for (i=1; i<=(*nsel); i++) {
    th2= th[i]*th[i];
    ans[i][i]= XtX[sel[i-1]*(*p)+sel[i-1]]/(*phi) + 6.0*(*tau)*(*phi)/(th2*th2) - 2.0/th2;
  }
  for (i=1; i<=(*nsel); i++) {
    for (j=i+1; j<=(*nsel); j++) {
      ans[i][j]= ans[j][i]= XtX[sel[i-1]*(*p)+sel[j-1]]/(*phi);
    }
  }
}


void imomModeK(double *th, PolynomialRootFinder::RootStatus_T *status, double *XtX, double *ytX, double *phi, double *tau, int *sel, int *nsel, int *p) {
  //piMOM mode when phi is known using gradient algorithm
  // - th: contains initial estimate at input and mode at output
  // - status: indicates if root finding has been successful
  bool found=false;
  int i, j, niter=0, root_count;
  double err= 1.0, *coef, *real_vector, *imag_vector;
  Polynomial poly;

  coef= dvector(0,4);
  real_vector= dvector(0,4);
  imag_vector= dvector(0,4);

  coef[0]= 2.0 * (*tau) * (*phi);
  coef[1]= 0.0;
  coef[2]= -2;
  while ((err > 1.0e-5) & (niter<50)) {
    err= 0;
    for (i=1; i<=(*nsel); i++) {
      coef[3]= ytX[sel[i-1]];
      for (j=1; j<i; j++) { coef[3]-= XtX[sel[i-1]*(*p)+sel[j-1]] * th[j]; }
      for (j=i+1; j<=(*nsel); j++) { coef[3]-= XtX[sel[i-1]*(*p)+sel[j-1]] * th[j]; }
      coef[3]= coef[3]/(*phi);
      coef[4]= -XtX[sel[i-1]*(*p)+sel[i-1]]/(*phi);
      poly.SetCoefficients(coef, 4);
      (*status)= poly.FindRoots(real_vector,imag_vector,&root_count);

      j=0; found= false;
      while ((!found) & (j<=4)) {
	if (fabs(imag_vector[j])<1.0e-5) {
	  if (((real_vector[j]>0) & (th[i]>0)) | ((real_vector[j]<0) & (th[i]<0))) {
	    err= max_xy(err,fabs(th[i] - real_vector[j]));
	    th[i]= real_vector[j];
	    found= true;
	  }
	}
	j++;
      }

    }
    niter++;
  }

  free_dvector(coef,0,4); free_dvector(real_vector,0,4); free_dvector(imag_vector,0,4);
}


void imomIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *sel, int *nsel, int *n, int *p, double *XtX, double *ytX, double *phi, double *tau, int *logscale, int *hessian) {
  int iter, maxit=100, emptyint;
  double **V, **Vinv, ftol= 1.0e-5, **dirth, **Vopt, detVopt, emptydouble=0, **emptymatrix;
  PolynomialRootFinder::RootStatus_T status;

  V= dmatrix(1,*nsel,1,*nsel); Vinv= dmatrix(1,*nsel,1,*nsel); Vopt= dmatrix(1,*nsel,1,*nsel); dirth= dmatrix(1,*nsel,1,*nsel);
  emptymatrix= dmatrix(1,1,1,1);
  //Initialize
  addct2XtX(tau,XtX,sel,nsel,p,V); //add tau to XtX diagonal, store in V
  inv_posdef_upper(V,*nsel,Vinv);
  Asym_xsel(Vinv,*nsel,ytX,sel,thopt);  //product Vinv * selected elements in ytX
  //Minimization
  imomModeK(thopt,&status,XtX,ytX,phi,tau,sel,nsel,p);
  set_f2opt_pars(&emptydouble,emptymatrix,&emptydouble,XtX,ytX,&emptydouble,&emptydouble,phi,tau,&emptyint,n,p,sel,nsel);
  if (status == PolynomialRootFinder::SUCCESS) {
    (*fopt)= f2opt_imom(thopt);
  } else {
    ddiag(dirth,1,*nsel);
    minimize(thopt, dirth, *nsel, ftol, &iter, fopt, f2opt_imom, maxit);
  }

  if (*hessian == 1) {
    //Laplace approx
    fppimomNegC_non0(Vopt,thopt,XtX,ytX,phi,tau,n,p,sel,nsel);
    invdet_posdef(Vopt,*nsel,Voptinv,&detVopt);
    (*ILaplace)= -(*fopt) - 0.5*log(detVopt);
  } else {
    (*ILaplace)= -(*fopt) - 0.5*(*nsel)*log(*n +.0);  //BIC-type approximation
  }

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
  double *y=REAL(Sy), *sumy2=REAL(Ssumy2), *XtX=REAL(SXtX), *ytX=REAL(SytX), *phi=REAL(Sphi), *tau=REAL(Stau), *rans, emptydouble=0, offset=0, *taualpha=NULL;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars,n,p,y,sumy2,&emptydouble,XtX,ytX,method,B,&emptydouble,&emptydouble,phi,tau,taualpha,&r,&emptydouble,&emptydouble,logscale,&offset);
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pimomMarginalKC(sel, nsel, &pars);
  UNPROTECT(1);
  return ans;
}


double pimomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  int one=1, hessian;
  double k, ans, m, s, ILaplace, *thopt, **Voptinv, fopt;
  thopt= dvector(1,*nsel);
  Voptinv= dmatrix(1,*nsel,1,*nsel);
  if ((*nsel)==0) {
    m= 0;
    s= sqrt(*(*pars).phi);
    ans= dnormC_jvec((*pars).y,*(*pars).n,m,s,1);
  } else {
    if (*(*pars).method == 2) { hessian=0; } else { hessian=1; }
    imomIntegralApproxC(&ILaplace,thopt,Voptinv,&fopt,sel,nsel,(*pars).n,(*pars).p,(*pars).XtX,(*pars).ytX,(*pars).phi,(*pars).tau,&one,&hessian);
    k= .5*((*nsel)*log(*(*pars).tau) - (*(*pars).sumy2)/(*(*pars).phi) - (*(*pars).n +.0) *LOG_M_2PI - (*(*pars).n - *nsel)*log(*(*pars).phi) - (*nsel)*LOG_M_PI);
    if (((*(*pars).method)==0) || ((*(*pars).method)==2)) {
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
  bool posdef;
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
  choldc(Vprop,*nsel,cholVprop,&posdef);
  choldc_inv(Vprop,*nsel,cholVpropinv,&posdef);
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




double f2opt_imomU(double *th) {
  //last element in th corresponds to eta=log(tau), i.e. log residual variance
  double ans;
  ans= fimomUNegC_non0(th+1,f2opt_pars.sumy2,f2opt_pars.XtX,f2opt_pars.ytX,f2opt_pars.alpha,f2opt_pars.lambda,f2opt_pars.tau,f2opt_pars.n,f2opt_pars.p,f2opt_pars.sel,f2opt_pars.nsel);
  return(ans);
}

double fimomUNegC_non0(double *th, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *n, int *p, int *sel, int *nsel) {
//loops over all elem in th (i.e. th has length *nsel+1 and contains non-zero elements only). th is indexed at 0.
//Note: last element in th corresponds to eta=log(tau), i.e. log residual variance
  int i;
  double ans, ytXth, sumlogth, suminvth, th2, eta, phi;
  eta= th[*nsel]; phi= exp(eta);
  for (i=0, ytXth=0, sumlogth=0, suminvth=0; i<(*nsel); i++) {
    ytXth+= ytX[sel[i]] * th[i];
    th2= th[i] * th[i];
    suminvth+= 1/th2;
    sumlogth+= log(th2);
  }
  ans= .5*(*lambda + *sumy2 - 2*ytXth + quadratic_xtAselx(th, XtX, p, nsel, sel))/phi + (*tau)*phi*suminvth + sumlogth + .5*eta*(*n - *nsel + *alpha);
  return ans;
}


//Hessian of fimomUNegC
// - ans: hessian matrix evaluated at th (indexed at 1, i.e. ans[1:(*nsel)+1][1:(*nsel)+1])
// - th: th[1:(*nsel)+1] indicates point at which to evaluate the hessian.
// - Other arguments as in fimomNegC_non0
void fppimomUNegC_non0(double **ans, double *th, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i, j;
  double th2, eta, phi, suminvth, ytXth, *XtXth;

  XtXth= dvector(1,*nsel);
  eta= th[*nsel +1]; phi= exp(eta);
  Asel_x(XtX,*p,th,*nsel,sel-1,XtXth);
  for (i=1, ytXth=0, suminvth=0; i<=(*nsel); i++) {
    th2= th[i]*th[i];
    ans[i][i]= XtX[sel[i-1]*(*p)+sel[i-1]]/phi + 6.0*(*tau)*phi/(th2*th2) - 2.0/th2;
    ans[i][*nsel+1]= ans[*nsel+1][i]= -2.0*(*tau)*phi/(th2*th[i]) - (XtXth[i]-ytX[sel[i-1]])/phi;
    ytXth+= ytX[sel[i-1]] * th[i];
    suminvth+= 1/(th[i]*th[i]);
  }
  for (i=1; i<=(*nsel); i++) {
    for (j=i+1; j<=(*nsel); j++) {
      ans[i][j]= ans[j][i]= XtX[sel[i-1]*(*p)+sel[j-1]]/phi;
    }
  }
  ans[*nsel+1][*nsel+1]= .5*(*lambda + *sumy2 - 2*ytXth + quadratic_xtAselx(th+1, XtX, p, nsel, sel))/phi + (*tau)*phi*suminvth;
  free_dvector(XtXth,1,*nsel);
}

void imomModeU(double *th, PolynomialRootFinder::RootStatus_T *status, double *sumy2, double *XtX, double *ytX, double *tau, double *alpha, double *lambda, int *sel, int *nsel, int *n, int *p) {
  //piMOM mode when phi is unknown using gradient algorithm
  // - th: contains initial estimate (theta,log(phi)) at input and mode at output
  // - status: indicates if root finding has been successful
  bool found=false;
  int i, j, niter=0, root_count;
  double err= 1.0, *coef, *real_vector, *imag_vector, phi, phinew, a, b, b2, c, d, suminvth2, *XtXth;
  Polynomial poly;

  phi= exp(th[*nsel +1]);
  b= (*n -(*nsel) + *alpha)/2.0;
  b2= b*b;

  coef= dvector(0,4);
  real_vector= dvector(0,4);
  imag_vector= dvector(0,4);
  XtXth= dvector(1,*nsel);

  coef[1]= 0.0;
  coef[2]= -2;
  while ((err > 1.0e-5) & (niter<50)) {
    coef[0]= 2.0 * (*tau) * phi;
    suminvth2= 0.0;
    err= 0;
    //Update th
    for (i=1; i<=(*nsel); i++) {
      coef[3]= ytX[sel[i-1]];
      for (j=1; j<i; j++) { coef[3]-= XtX[sel[i-1]*(*p)+sel[j-1]] * th[j]; }
      for (j=i+1; j<=(*nsel); j++) { coef[3]-= XtX[sel[i-1]*(*p)+sel[j-1]] * th[j]; }
      coef[3]= coef[3]/phi;
      coef[4]= -XtX[sel[i-1]*(*p)+sel[i-1]]/phi;
      poly.SetCoefficients(coef, 4);
      (*status)= poly.FindRoots(real_vector,imag_vector,&root_count);

      j=0; found= false;
      while ((!found) & (j<=4)) {
	if (fabs(imag_vector[j])<1.0e-5) {
	  if (((real_vector[j]>0) & (th[i]>0)) | ((real_vector[j]<0) & (th[i]<0))) {
	    err= max_xy(err, fabs(th[i] - real_vector[j]));
	    th[i]= real_vector[j];
	    suminvth2 += 1.0/(th[i]*th[i]);
	    found= true;
	  }
	}
	j++;
      }
    }

    //Update phi
    a= (*tau) * suminvth2;
    c= 0;
    Asel_x(XtX,*p,th,*nsel,sel-1,XtXth);
    for (i=1; i<=(*nsel); i++) { c += -2.0*ytX[sel[i-1]]*th[i] + th[i]*XtXth[i]; }
    c= -.5*(*lambda + *sumy2 + c);
    d= sqrt(b2 - 4.0*a*c);

    if (-b > d) { phinew= (-b-d)/(2.0*a); } else { phinew= (-b+d)/(2.0*a); }
    err= max_xy(err,fabs(phi-phinew));
    phi= phinew;

    niter++;
  }

  th[*nsel +1]= log(phi);

  free_dvector(coef,0,4); free_dvector(real_vector,0,4); free_dvector(imag_vector,0,4);
  free_dvector(XtXth,1,*nsel);
}


void imomUIntegralApproxC(double *ILaplace, double *thopt, int *sel, int *nsel, int *n, int *p, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *logscale, int *hessian) {
  int iter, maxit=100, emptyint;
  double ftol= 1.0e-10, **dirth, **Vopt, **Voptinv, detVopt, emptydouble=0, **emptymatrix, fopt;
  PolynomialRootFinder::RootStatus_T status;

  Vopt= dmatrix(1,*nsel +1,1,*nsel +1); Voptinv= dmatrix(1,*nsel +1,1,*nsel +1);
  dirth= dmatrix(1,*nsel +1,1,*nsel +1);
  emptymatrix= dmatrix(1,1,1,1);
  //Initialize
  set_f2opt_pars(&emptydouble,emptymatrix,sumy2,XtX,ytX,alpha,lambda,&emptydouble,tau,&emptyint,n,p,sel,nsel);

  //Minimization
  imomModeU(thopt,&status,sumy2,XtX,ytX,tau,alpha,lambda,sel,nsel,n,p);
  set_f2opt_pars(&emptydouble,emptymatrix,sumy2,XtX,ytX,alpha,lambda,&emptydouble,tau,&emptyint,n,p,sel,nsel);
  if (status == PolynomialRootFinder::SUCCESS) {
    fopt= f2opt_imomU(thopt);
  } else {
    ddiag(dirth,1,*nsel +1);
    minimize(thopt, dirth, *nsel +1, ftol, &iter, &fopt, f2opt_imomU, maxit);
  }

  if (*hessian ==1) {
    //Laplace approx
    fppimomUNegC_non0(Vopt,thopt,sumy2,XtX,ytX,alpha,lambda,tau,n,p,sel,nsel);
    invdet_posdef(Vopt,*nsel +1,Voptinv,&detVopt);
    (*ILaplace)= -fopt - 0.5*log(detVopt) + .5*(*nsel)*log(2.0*(*tau));
  } else {
    (*ILaplace)= -fopt - 0.5*(*nsel)*log(*n +.0) + .5*(*nsel)*log(2.0*(*tau));  //BIC-type approximation
  }

  free_dmatrix(Vopt,1,*nsel +1,1,*nsel +1); free_dmatrix(Voptinv,1,*nsel +1,1,*nsel +1);
  free_dmatrix(dirth,1,*nsel +1,1,*nsel+1); free_dmatrix(emptymatrix,1,1,1,1);
  if ((*logscale)!=1) { (*ILaplace)= exp(*ILaplace); }
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
// - method: method to approximate the integral for known phi. 0 for Laplace approx which may underestimate true value, 1 for exact evaluation which can be very computationally expensive. 2 ('Hybrid') integrates wrt phi numerically (Romberg) and wrt theta via Laplace approx (Laplace error is adjusted based on exact evaluation for single phi value)
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale
// - alpha, lambda: prior for phi (residual variance) is Inverse Gamma (.5*alpha,.5*lambda)
double f2int_imom(double phi) {
  int one=1, *inputlog= f2int_pars.logscale;
  double ans, *inputphi= f2int_pars.phi;
  f2int_pars.phi= &phi;
  f2int_pars.logscale= &one;
  ans= exp(pimomMarginalKC(f2int_pars.sel,f2int_pars.nsel,&f2int_pars) + dinvgammaC(phi,.5*(*f2int_pars.alpha),.5*(*f2int_pars.lambda),1) - *(f2int_pars.offset));
  f2int_pars.phi= inputphi;
  f2int_pars.logscale= inputlog;
  return(ans);
}


SEXP pimomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda) {
  int *sel=INTEGER(Ssel), *nsel=INTEGER(Snsel), *n=INTEGER(Sn), *p=INTEGER(Sp), *method=INTEGER(Smethod), *B=INTEGER(SB), *logscale=INTEGER(Slogscale), r=1;
  double *y=REAL(Sy), *sumy2=REAL(Ssumy2), *x=REAL(Sx), *XtX=REAL(SXtX), *ytX=REAL(SytX), *tau=REAL(Stau), *alpha=REAL(Salpha), *lambda=REAL(Slambda), *rans, emptydouble=0, offset=0, *taualpha=NULL;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars,n,p,y,sumy2,x,XtX,ytX,method,B,alpha,lambda,&emptydouble,tau,taualpha,&r,&emptydouble,&emptydouble,logscale,&offset);
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pimomMarginalUC(sel, nsel, &pars);
  UNPROTECT(1);
  return ans;
}


double pimomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  int i, j, zero=0, one=1, *inputlog, hessian;
  double ans=0, er, sumer2, **V, **Vinv, *thest, ypred, phiest, intmc, intlapl, *inputphi, num, den, term1, alphahalf=.5*(*(*pars).alpha);

  if (*nsel ==0) {
    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*LOG_M_PI + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);
    if ((*(*pars).logscale)!=1) ans= exp(ans);
  } else {
    V= dmatrix(1,*nsel,1,*nsel);
    Vinv= dmatrix(1,*nsel,1,*nsel);
    thest= dvector(1,*nsel+1);

    addct2XtX((*pars).tau,(*pars).XtX,sel,nsel,(*pars).p,V); //add tau to diagonal elem of XtX
    inv_posdef_upper(V,*nsel,Vinv);
    Asym_xsel(Vinv,*nsel,(*pars).ytX,sel,thest);
    for (i=0, sumer2=0; i<(*(*pars).n); i++) {
      for (j=1, ypred=0; j<=(*nsel); j++) { ypred += (*pars).x[i + (*(*pars).n)*sel[j-1]] * thest[j]; }
      er= (*pars).y[i] - ypred;
      sumer2+= er*er;
    }
    phiest= (sumer2 + (*(*pars).lambda))/(*(*pars).alpha + *(*pars).n);
    if (((*(*pars).method)==0) || ((*(*pars).method)==2)) {  //Laplace or Plug-in

      if (*(*pars).method == 2) { hessian=0; } else { hessian=1; }
      thest[*nsel +1]= log(phiest);
      imomUIntegralApproxC(&ans,thest,sel,nsel,(*pars).n,(*pars).p,(*pars).sumy2,(*pars).XtX,(*pars).ytX,(*pars).alpha,(*pars).lambda,(*pars).tau,&one,&hessian);
      ans= ans + alphahalf*log(.5*(*(*pars).lambda)) - .5*(*(*pars).n)*LOG_M_2PI - gamln(&alphahalf);
      if ((*(*pars).logscale)!=1) ans= exp(ans);

    } else if ((*(*pars).method)==1) {  //MC for each fixed phi + univariate integration
      set_f2int_pars((*pars).XtX,(*pars).ytX,(*pars).tau,(*pars).n,(*pars).p,sel,nsel,(*pars).y,(*pars).sumy2,(*pars).method,(*pars).B,(*pars).alpha,(*pars).lambda,&zero);
      inputphi= (*pars).phi; (*pars).phi= &phiest;
      (*(*pars).method)= 0; inputlog= (*pars).logscale; (*pars).logscale= &one; //Laplace approx for phi=phiest
      intlapl= pimomMarginalKC(sel, nsel, pars);
      (*pars).phi= inputphi; (*(*pars).method)= 1; (*pars).logscale= inputlog;  //reset input values for phi, method
      f2int_pars.offset= &intlapl; //f2int_imom returns result divided by exp(intlapl) to avoid numerical overflow
      ans= intlapl + log(qromo(f2int_imom,0.0,100,midpnt) + qromo(f2int_imom,100,1.0e30,midinf));
      if ((*(*pars).logscale)==0) ans= exp(ans);

    } else if ((*(*pars).method)==3) {  //Hybrid Laplace - MC - Univariate integration
      set_f2int_pars((*pars).XtX,(*pars).ytX,(*pars).tau,(*pars).n,(*pars).p,sel,nsel,(*pars).y,(*pars).sumy2,(*pars).method,(*pars).B,(*pars).alpha,(*pars).lambda,&zero);
      inputphi= (*pars).phi; (*pars).phi= &phiest;
      (*(*pars).method)= 1; //IS evaluation of marginal for phi=phiest
      intmc= pimomMarginalKC(sel, nsel, pars);
      (*(*pars).method)= 0; //Laplace approx for phi=phiest
      intlapl= pimomMarginalKC(sel, nsel, pars);
      (*pars).phi= inputphi; (*(*pars).method)= 2;  //reset input values for phi, method
      if (intlapl==0) { intmc+= 1.0e-300; intlapl+= 1.0e-300; } //avoid numerical zero
      f2int_pars.method= &zero;  //set method to eval marginal for known phi to Laplace approx
      f2int_pars.offset= &intlapl; //f2int_imom returns result divided by exp(intlapl) to avoid numerical overflow
      ans= intmc + log(qromo(f2int_imom,0.0,100,midpnt) + qromo(f2int_imom,100,1.0e30,midinf)); //adjusment is intmc - intlapl, but intlapl is the offset so needs to added back in
      if ((*(*pars).logscale)==0) ans= exp(ans);
    }
    free_dmatrix(V,1,*nsel,1,*nsel);
    free_dmatrix(Vinv,1,*nsel,1,*nsel);
    free_dvector(thest,1,*nsel+1);
  }
  return(ans);
}



//*************************************************************************************
// Product eMOM routines
//*************************************************************************************

double pemomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  return 0.0;
}

double pemomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  return 0.0;
}



//*************************************************************************************
// Zellner's prior routines
//*************************************************************************************

// Marginal likelihood for linear models under Zellner's prior
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
// - logscale: if set to 1 result is returned in log scale

SEXP zellnerMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Slogscale) {
  struct marginalPars pars;
  int emptyint=0;
  double *rans, emptydouble=0, offset=0, *taualpha=NULL;
  SEXP ans;

  set_marginalPars(&pars,INTEGER(Sn),INTEGER(Sp),REAL(Sy),REAL(Ssumy2),&emptydouble,REAL(SXtX),REAL(SytX),&emptyint,&emptyint,&emptydouble,&emptydouble,REAL(Sphi),REAL(Stau),taualpha,&emptyint,&emptydouble,&emptydouble,INTEGER(Slogscale),&offset);
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= zellnerMarginalKC(INTEGER(Ssel),INTEGER(Snsel),&pars);
  UNPROTECT(1);
  return ans;
}


double zellnerMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  int i,j;
  double *m, s, **S, **Sinv, detS, num, den, adj, tau= *(*pars).tau, logphi= log(*(*pars).phi), ans=0.0, zero=0;

  if (*nsel ==0) {

    m= dvector(1,1);
    m[1]=0; s= sqrt(*(*pars).phi);
    ans= dnormC_jvec((*pars).y,*(*pars).n,m[1],s,1);
    free_dvector(m,1,1);

  } else {

    m= dvector(1,*nsel);
    S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);  //copy XtX into S
    adj= (tau+1)/tau;
    for (i=1; i<=(*nsel); i++) {
      S[i][i]= S[i][i] * adj;
      for (j=i+1; j<=(*nsel); j++) { S[i][j]= S[i][j] * adj; }
    }
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);

    num= -.5*(*(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel))/(*(*pars).phi);
    den= .5*((*(*pars).n +.0)*(LOG_M_2PI+logphi) + (*nsel)*log(tau+1.0));
    //den= .5*((*(*pars).n +.0)*(LOG_M_2PI+logphi) + log(detS) + (*nsel)*logtau);
    ans= num - den;

    free_dvector(m,1,*nsel);
    free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}


SEXP zellnerMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Slogscale, SEXP Salpha, SEXP Slambda) {
  int emptyint=0;
  double *rans, emptydouble=0, offset=0, *taualpha=NULL;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars,INTEGER(Sn),INTEGER(Sp),REAL(Sy),REAL(Ssumy2),REAL(Sx),REAL(SXtX),REAL(SytX),&emptyint,&emptyint,REAL(Salpha),REAL(Slambda),&emptydouble,REAL(Stau),taualpha,&emptyint,&emptydouble,&emptydouble,INTEGER(Slogscale),&offset);
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= zellnerMarginalUC(INTEGER(Ssel), INTEGER(Snsel), &pars);
  UNPROTECT(1);
  return ans;
}

double zellnerMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  int i, j;
  double num, den, ans=0.0, term1, *m, **S, **Sinv, detS, adj, tau= *(*pars).tau, nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda), ss, zero=0;
  if (*nsel ==0) {

    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*(LOG_M_PI) + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);

  } else {

    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);  //copy XtX onto S
    adj= (tau+1)/tau;
    for (i=1; i<=(*nsel); i++) {
      S[i][i]= S[i][i] * adj;
      for (j=i+1; j<=(*nsel); j++) { S[i][j]= S[i][j] * adj; }
    }
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    nuhalf= .5*(*(*pars).n + *(*pars).alpha);

    ss= *(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel);
    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2.0) - log(ss));
    den= .5*(*(*pars).n * LOG_M_2PI) + .5 * (*nsel) *log((*(*pars).tau)+1.0) + gamln(&alphahalf);
    //den= .5*(*(*pars).n * LOG_M_2PI + log(detS)) + .5 * (*nsel) *log(*(*pars).tau) + gamln(&alphahalf);
    ans= num - den;

    free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);

  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}
