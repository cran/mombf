/***********************************************************
 Basic, input-output and matrix manipulation

 Authors. Peter Mueller, Stephen Morris, David Rossell
          (some routines obtained from other sources)
 Last modified. 06 2008
***********************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cstat.h"

static const char interface_c_sccs_id[] = "%W%";
static const char mess_c_sccs_id[] = "%W%";
static const char nrutil_c_sccs_id[] = "%W%"; 
static const char vector_c_sccs_id[] = "%W%";
static const char css_c_sccs_id[] = "@(#)$Workfile: rand.c$ $Revision: 5$";
static long is1=123456789, is2=981963; 
static int set=0;

FILE	*ifile, *fopen(), *ofile;
int nv = 0;

long Xm1,Xm2,Xa1,Xa2,Xcg1[32],Xcg2[32],Xa1w,Xa2w,Xig1[32],Xig2[32],Xlg1[32],
    Xlg2[32],Xa1vw,Xa2vw;
long Xqanti[32];


/************************************************************************
                         STATISTICAL FUNCTIONS
************************************************************************/

/* Sample mean of elements 0 through lim of vector x */
double meani(int *x, int lim) 
{
    int i;
    double value;
    for(i=0,value=0; i<=lim; i++) { value += x[i]; }
    value *= 1.0/(lim+1.0);
    return value;
}

double wmeani(int *x, int lim, double *w)
{
    int i;
    double value, wtot;
    for(i=0,value=0,wtot=0; i<=lim; i++) { value += w[i]*x[i]; wtot += w[i]; }
    value *= 1.0/wtot;
    return value;
}

/* Sample mean of elements 0 through lim (both included) of vector x */
double meanx(double *x, int lim) 
{
    int i;
    double value;
    for(i=0,value=0; i<=lim; i++) value += x[i];
    value *= 1.0/(lim+1.0);
    return value;
}

double wmeanx(double *x, int lim, double *w)
{
    int i;
    double value, wtot;
    for(i=0,value=0,wtot=0; i<=lim; i++) {
        value += w[i]*x[i]; wtot += w[i];
    }
    value *= 1.0/wtot;
    return value;
}


/* Sample variance of elements 0 through lim of vector x, unb=1 returns unbiased version */
double vari(int *x, int lim, int unb) 
{
    int i;
    double value;
    for(i=0,value=0; i<=lim; i++) { value += pow(x[i],2) / (1.0+lim); }
    value -= pow(meani(x,lim),2);
    if (unb==1 && lim>0) value *= (1.0+lim)/(0.0+lim);
    return value;
}

double wvari(int *x, int lim, double *w) 
{
    int i;
    double value, wtot;
    for(i=0,value=0,wtot=0; i<=lim; i++) { value += w[i]*pow(x[i],2); wtot += w[i]; }
    value = value/wtot - pow(wmeani(x,lim,w),2);
    return value;
}


/* Sample variance of elements 0 through lim of vector x, unb=1 returns unbiased version */
double varx(double *x, int lim, int unb) 
{
    int i;
    double value;
    for(i=0,value=0; i<=lim; i++) { value += pow(x[i],2) / (1.0+lim); }
    value -= pow(meanx(x,lim),2);
    if (unb==1 && lim>0) value *= (1.0+lim)/(lim+0.0);
    return value;
}

double wvarx(double *x, int lim, double *w) 
{
    int i;
    double value, wtot;
    for(i=0,value=0,wtot=0; i<=lim; i++) { value += w[i]*pow(x[i],2); wtot += w[i]; }
    value = value/wtot - pow(wmeanx(x,lim,w),2);
    return value;
}


double cv(double *x, int ini, int fi) {
  //Compute coefficient of variation of x[ini..fi]
  int i; double m, s, ans;

  for (i=ini, m=s=0; i<=fi; i++) {
    m+= x[i];
    s+= x[i]*x[i];
  }
  m= m/(1.0+fi-ini); 
  s= s/(.0+fi-ini) - m*m*(1.0+fi-ini)/(.0+fi-ini);
  ans= sqrt(s)/m;
  return(ans);
}


double cvinv(double *x, int ini, int fi) {
  //Compute coefficient of variation of x[ini..fi]
  int i; double m, s, ans;

  for (i=ini, m=s=0; i<=fi; i++) {
    m+= 1.0/x[i];
    s+= 1.0/(x[i]*x[i]);
  }
  m= m/(1.0+fi-ini); 
  s= s/(.0+fi-ini) - m*m*(1.0+fi-ini)/(.0+fi-ini);
  ans= sqrt(s)/m;
  return(ans);
}


/* Column means */
void colMeans(double *m, double *x, int nrow, int ncol) {
  //x is assumed to be in row order, that is x[0], x[1] ... x[ncol-1] are elem in 1st row
  int i,j;
  for (j=0; j<ncol; j++) { m[j]= 0; }
  for (i=0; i<nrow; i++) {
    for (j=0; j<ncol; j++) {
      m[j] += x[i*ncol+j];
    }
  }
  for (j=0; j<ncol; j++) { m[j]= m[j]/(nrow+.0); }
}

void colVar(double *v, double *x, int nrow, int ncol) {
  //x is assumed to be in row order, that is x[0], x[1] ... x[ncol-1] are elem in 1st row
  int i,j; double *m, *m2;

  m= dvector(0,ncol-1); m2= dvector(0,ncol-1);
  for (j=0; j<ncol; j++) { m[j]= m2[j]= 0; }
  for (i=0; i<nrow; i++) {
    for (j=0; j<ncol; j++) {
      m[j]+= x[i*ncol+j];
      m2[j] += x[i*ncol+j]*x[i*ncol+j];
    }
  }
  for (j=0; j<ncol; j++) { m[j]= m[j]/(.0+nrow); v[j]= m2[j]/(nrow-1.0) - m[j]*m[j]*(nrow+.0)/(nrow-1.0); }
  free_dvector(m,0,ncol-1); free_dvector(m2,0,ncol-1);

}


void colCV(double *cv, double *x, int nrow, int ncol) {
  //x is assumed to be in row order, that is x[0], x[1] ... x[ncol-1] are elem in 1st row
  int i,j; double *m, *s;
  m= dvector(0,ncol); s= dvector(0,ncol);
  for (j=0; j<ncol; j++) { m[j]= s[j]= 0; }
  for (i=0; i<nrow; i++) {
    for (j=0; j<ncol; j++) {
      m[j] += x[i*ncol+j];
      s[j] += x[i*ncol+j]*x[i*ncol+j];
    }
  }
  for (j=0; j<ncol; j++) { 
    m[j]= m[j]/(nrow+.0); 
    s[j]= s[j]/(nrow-1.0) - m[j]*m[j]*(nrow+.0)/(nrow-1.0);
    cv[j]= sqrt(s[j])/m[j];
  }
  free_dvector(m,0,ncol); free_dvector(s,0,ncol);
}


void colCVinv(double *cv, double *x, int nrow, int ncol) {
  //x is assumed to be in row order, that is x[0], x[1] ... x[ncol-1] are elem in 1st row
  int i,j; double *m, *s;
  m= dvector(0,ncol); s= dvector(0,ncol);
  for (j=0; j<ncol; j++) { m[j]= s[j]= 0; }
  for (i=0; i<nrow; i++) {
    for (j=0; j<ncol; j++) {
      m[j] += 1.0/x[i*ncol+j];
      s[j] += 1.0/(x[i*ncol+j]*x[i*ncol+j]);
    }
  }
  for (j=0; j<ncol; j++) { 
    m[j]= m[j]/(nrow+.0); 
    s[j]= s[j]/(nrow-1.0) - m[j]*m[j]*(nrow+.0)/(nrow-1.0);
    cv[j]= sqrt(s[j])/m[j];
  }
  free_dvector(m,0,ncol); free_dvector(s,0,ncol);
}




/************************************************************************
                         BASIC BAYESIAN MODELS
************************************************************************/

/********************************************
 *         normal_normal
 ********************************************/
void nn_bayes(double *mpo, double **Spo, double **Spo_inv, int p, double r1, double *mpr, double **Spr_inv, double r2, double *y, double **Slik_inv)
/* prior: N(x; mpr, r1*Spr)
   likl:  N(y; x, r2*Slik)
   p: dimensionality
   returns:  post N(x; mpo,Spo)
   Spo = (1/r1*Spr_inv + 1/r2*Slik_inv)^-1
   mpo = Spo*(1/r1*Spr*mpr + 1/r2*Slik*y)
   NOTE: input vectors and matrices must start at position 1, not 0
*/
{ double *z;

  z = dvector(1,p);

  rA_plus_sB(1.0/r1, Spr_inv, 1.0/r2, Slik_inv, Spo_inv,1,p,1,p);
  inv_posdef(Spo_inv,p,Spo); 
  rAx_plus_sBy(1.0/r1, Spr_inv, mpr, 1.0/r2, Slik_inv, y, z,1,p,1,p);
  Ax(Spo,z,mpo,1,p,1,p);
  
  free_dvector(z,1,p);
}

void nn_bayes_rand(double *theta, int p, double r1, double **Spr_inv, double *mpr, double r2, double **Slik_inv, double *y) {
/* same as nn_bayes, but returns a draw only 
   returns:
   theta: draw from the posterior N(theta; mpo, Spo):
*/
  double *z, **S, **S_inv, *m, **cholS;

  /* allocate memory */
  z = dvector(0,p-1);
  m = dvector(0,p-1);
  S = dmatrix(0,p-1,0,p-1);
  S_inv = dmatrix(0,p-1,0,p-1);
  cholS= dmatrix(0,p-1,0,p-1);

  rA_plus_sB(1.0/r1, Spr_inv, 1.0/r2, Slik_inv, S_inv,1,p,1,p);
  inv_posdef(S_inv,p,S);
  rAx_plus_sBy(1.0/r1, Spr_inv, mpr, 1.0/r2, Slik_inv, y, z,1,p,1,p);
  Ax(S,z,m,1,p,1,p);

  choldc(S,p,cholS);
  rmvnormC(theta,p,m,cholS);

  free_dvector(z,0,p-1);
  free_dvector(m,0,p-1);
  free_dmatrix(S,0,p-1,0,p-1);
  free_dmatrix(S_inv,0,p-1,0,p-1);
  free_dmatrix(cholS,0,p-1,0,p-1);

}

double nn_integral(double *x, double *rx, double **Vxinv, double *detVx, double *mpr, double *rpr, double **Vprinv, double *detVpr, int *p, int *logscale) {
  // Compute normal-normal integral
  //   N(x;beta,rx*Vx) * N(beta;mpr,rpr*Vpr) with respect to beta
  //
  // - Vxinv, detVx: inverse and determinant of Vx
  // - Vprinv, detVpr: inverse and determinant of Vpr
  // - p: dimensionality
  // - logscale: if not 0, result is returned in log scale

  double ans, detsum, **Vsum, **Vsuminv, **cholVsum, *m;

  m= dvector(1,*p);
  Vsum= dmatrix(1,*p,1,*p); Vsuminv= dmatrix(1,*p,1,*p); cholVsum= dmatrix(1,*p,1,*p);

  rA_plus_sB(1.0/(*rx),Vxinv,1.0/(*rpr),Vprinv,Vsuminv,1,*p,1,*p);
  choldc_inv(Vsuminv,*p,cholVsum);
  detsum= choldc_det(cholVsum,*p);
  inv_posdef_chol(cholVsum,*p,Vsum);
  rAx_plus_sBy(1.0/(*rx),Vxinv,x,1.0/(*rpr),Vprinv,mpr,m,1,*p,1,*p);
  ans= xtAy(m,Vsum,m,1,*p) - xtAy(x,Vxinv,x,1,*p) - xtAy(mpr,Vprinv,mpr,1,*p);

  ans= .5*ans - 0.5*((*p+.0)*LOG_M_2PI + log(*detVx) + log(*detVpr) - log(detsum));
  if (*logscale != 0) ans= exp(ans);

  free_dvector(m,1,*p); 
  free_dmatrix(Vsum,1,*p,1,*p); free_dmatrix(Vsuminv,1,*p,1,*p); free_dmatrix(cholVsum,1,*p,1,*p);

  return(ans);
}


void lm (double *b, double **XtX, double **invXtX, double *Xty, double *s, double *ypred, double *y, double **X, int *n, int *p, int *useXtX) {
 //Fits classical multiple linear regression
 /* Input
      - y: response variable y[1..n]
      - X: design matrix X[1..n][1..p]
      - n: number of observations
      - p: number of covariates
      - useXtX: if set to 0 the inverse of X'X is computed, otherwise the supplied value is used
    Ouput
      - b: least-squares estimate for regression coefficients
      - XtX, invXtX: X'X and its inverse (if useXtX==0 they're ouput param, otherwise they're input)
      - Xty: vector X'y (if useXtX==0 it's an output param, otherwise it's input)
      - s: residual variance (dividing by n-p)
      - ypred: predicted values i.e. X'b
 */
  int i;

  if (*n<*p) errorC("lm","Linear model with more variables than observations",0);

  if (*useXtX==0) {
    AtB(X,1,*n,1,*p,X,1,*n,1,*p,XtX);
    inv_posdef(XtX,*p,invXtX);
    Atx(X,y,Xty,1,*n,1,*p); //X'y
  }

  Ax(invXtX,Xty,b,1,*p,1,*p); //least squares estimate
  Ax(X,b,ypred,1,*n,1,*p); //predicted values

  for (*s= 0, i=1; i<=(*n); i++) { (*s) += (y[i]-ypred[i])*(y[i]-ypred[i]); }
  (*s)= (*s)/(*n- *p);

}

void lmbayes (double *bpost, double *spost, double *b, double **Vb, double *a_s, double *b_s, double **XtX, double **invXtX, double *Xty, int *B, double *y, double **X, int *n, int *p, int *useXtX, double *mpr, double **Spr_inv, double *tauprior, double *nu0, double *s0) {
 //Bayesian conjugate multiple linear regression
   // y ~ N(X'beta,sigma^2)
   // beta ~ N(mpr,sigma^2*Spr) (if tauprior<=0)
   // beta ~ N(mpr,tauprior*sigma^2*(X'X)^{-1}) (if tauprior>0) 
   //        e.g. tauprior==n gives unit information prior
   // sigma^2 ~ IG(.5*nu0,.5*s0)  (nu0: prior sample size; s0: prior sum of squares)
 /* Input
      - B: number of posterior samples to draw (can be 0)
      - y: response variable y[1..n]
      - X: design matrix X[1..n][1..p]
      - n: number of observations
      - p: number of covariates
      - useXtX: if set to 0 the inverse of X'X is computed, otherwise the supplied value is used
      - mpr, Spr_inv, tauprior: prior parameters for beta
      - nu0, s0: prior for sigma^2 is IG(.5*nu0,.5*s0)
    Ouput
      - bpost: matrix (B rows, p cols) with samples from the posterior of beta. Starts at bpost[1].
      - spost: vector (B rows) with samples from the posterior of sigma^2. Starts at spost[1].
      - b, Vb: posterior for regression coef is N(b,sigma^2*Vb)
      - a_s, b_s: posterior for sigma^2 is IG(a_s,b_s)
      - XtX, invXtX: X'X and its inverse (if useXtX==0 they're ouput param, otherwise they're input)
      - Xty: vector X'y (if useXtX==0 it's an output param, otherwise it's input)
 */
  int i, j, one=1;
  double *b_ls, s_ls, *ypred, **Vb_inv, *zeroes, **cholVb;

  if (*useXtX==0) {
    AtB(X,1,*n,1,*p,X,1,*n,1,*p,XtX);
    inv_posdef(XtX,*p,invXtX);
    Atx(X,y,Xty,1,*n,1,*p); //X'y
  }

  b_ls= dvector(1,*p); ypred= dvector(1,*n);
  lm(b_ls,XtX,invXtX,Xty,&s_ls,ypred,y,X,n,p,&one);  //least-squares fit

  *a_s= .5*(*nu0 + *n); *b_s= .5*(*s0 + (*n- *p)*s_ls); //posterior for sigma^2

  Vb_inv= dmatrix(1,*p,1,*p);   //posterior for beta
  if (*tauprior > 0) {
    nn_bayes(b,Vb,Vb_inv,*p,*tauprior,mpr,XtX,1.0,b_ls,XtX);
  } else {
    nn_bayes(b,Vb,Vb_inv,*p,1.0,mpr,Spr_inv,1.0,b_ls,XtX);
  }

  if (*B>0) {             //posterior samples
    cholVb= dmatrix(1,*p,1,*p);
    choldc(Vb,*p,cholVb); //cholesky decomp of posterior covar for beta
    zeroes= dvector(1,*p);
    for (i=1; i<=(*p); i++) { zeroes[i]= 0; }
    for (i=1; i<=(*B); i++) {
      spost[i]= 1.0/rgammaC(*a_s,*b_s);
      rmvnormC(bpost+(i-1)*(*p),*p,zeroes,cholVb);
      for (j=1; j<=(*p); j++) { bpost[(i-1)*(*p)+j]= bpost[(i-1)*(*p)+j]*sqrt(spost[i])+b[j]; }
    }
    free_dvector(zeroes,1,*p);
    free_dmatrix(cholVb,1,*p,1,*p);
  }

  free_dvector(b_ls,1,*p); free_dvector(ypred,1,*n); free_dmatrix(Vb_inv,1,*p,1,*p);

}


void lmbayes_knownvar (double *bpost, double *b, double **Vb, double **XtX, double **invXtX, double *Xty, double *sigma, int *B, double *y, double **X, int *n, int *p, int *useXtX, double *mpr, double **Spr_inv, double *tauprior) {
 //Bayesian conjugate multiple linear regression with known var
   // y ~ N(X'beta,sigma^2)
   // beta ~ N(mpr,sigma^2*Spr) (if tauprior<=0)
   // beta ~ N(mpr,tauprior*sigma^2*(X'X)^{-1}) (if tauprior>0) 
   //        e.g. tauprior==n gives unit information prior
 /* Input
      - sigma: residual standard deviation 
      - B: number of posterior samples to draw (can be 0)
      - y: response variable y[1..n]
      - X: design matrix X[1..n][1..p]
      - n: number of observations
      - p: number of covariates
      - useXtX: if set to 0 the inverse of X'X is computed, otherwise the supplied value is used
      - mpr, Spr_inv, tauprior: prior parameters for beta
    Ouput
      - bpost: matrix (B rows, p cols) with samples from the posterior of beta. Starts at bpost[1].
      - b, Vb: posterior for regression coef is N(b,sigma^2*Vb)
      - XtX, invXtX: X'X and its inverse (if useXtX==0 they're ouput param, otherwise they're input)
      - Xty: vector X'y (if useXtX==0 it's an output param, otherwise it's input)
 */
  int i, j, one=1;
  double *b_ls, s_ls, *ypred, **Vb_inv, *zeroes, **cholVb;

  if (*useXtX==0) {
    AtB(X,1,*n,1,*p,X,1,*n,1,*p,XtX);
    inv_posdef(XtX,*p,invXtX);
    Atx(X,y,Xty,1,*n,1,*p); //X'y
  }

  b_ls= dvector(1,*p); ypred= dvector(1,*n);
  lm(b_ls,XtX,invXtX,Xty,&s_ls,ypred,y,X,n,p,&one);  //least-squares fit

  Vb_inv= dmatrix(1,*p,1,*p);   //posterior for beta
  if (*tauprior > 0) {
    nn_bayes(b,Vb,Vb_inv,*p,*tauprior,mpr,XtX,1.0,b_ls,XtX);
  } else {
    nn_bayes(b,Vb,Vb_inv,*p,1.0,mpr,Spr_inv,1.0,b_ls,XtX);
  }

  if (*B>0) {             //posterior samples
    cholVb= dmatrix(1,*p,1,*p);
    choldc(Vb,*p,cholVb); //cholesky decomp of posterior covar for beta
    zeroes= dvector(1,*p);
    for (i=1; i<=(*p); i++) { zeroes[i]= 0; }
    for (i=1; i<=(*B); i++) {
      rmvnormC(bpost+(i-1)*(*p),*p,zeroes,cholVb);
      for (j=1; j<=(*p); j++) { bpost[(i-1)*(*p)+j]= bpost[(i-1)*(*p)+j]*(*sigma)+b[j]; }
    }
    free_dvector(zeroes,1,*p);
    free_dmatrix(cholVb,1,*p,1,*p);
  }

  free_dvector(b_ls,1,*p); free_dvector(ypred,1,*n); free_dmatrix(Vb_inv,1,*p,1,*p);

}


/************************************************************************
                         INPUT/OUTPUT FUNCTIONS
************************************************************************/

/* open file for input */
FILE *openIn(char *name)
{
  if ((ifile=fopen(name,"r"))==NULL){
    fserror("openIn","open file",name);
  }
  return(ifile);
}

/* open file for output */

FILE *openOut(char *name)
{
  if ((ofile=fopen(name,"w"))==NULL){
    fserror("openOut","open file",name);
  }
  return ofile;
}


/* --------------------   read in   ------------------------- */

//void scanFloat(char *txt, float *f)
//{
//  fscanf(ifile,txt);
//  if (fscanf(ifile," %f ",f) != 1) {
//    fserror ("scanFloat","read float",txt);
//  }
//}
// 
//void scanDouble(char *txt,double *f)
//{
//  fscanf(ifile,txt);
//  if (fscanf(ifile," %lf ",f) != 1) {
//    fserror ("scanDouble","read double",txt);
//  }
//}
// 
//void fscanDouble(FILE *ifile, char *txt, double *f)
//{
//  fscanf(ifile,txt);
//  if (fscanf(ifile," %lf ",f) != 1) {
//    fserror ("fscanDouble","read double",txt);
//  }
//}
// 
//void scanInt(char *txt, int *n)
//{
//  int 	s;
// 
//  fscanf(ifile,txt);
//  if ((s = fscanf(ifile," %d ",n)) != 1) {
//    fserror ("scanInt","read int",txt);
//  }
//}
// 
//void fscanInt(FILE *ifile, char *txt, int *n)
//{
// 	int 	s;
// 
//  fscanf(ifile,txt);
//  if ((s = fscanf(ifile," %d ",n)) != 1) {
//    fserror ("fscanInt","read int",txt);
//  }
//}
// 
// 
//void scanLong(char *txt, long *n)
//{
//        int     s;
// 
//  fscanf(ifile,txt);
//  if ((s = fscanf(ifile," %ld ",n)) != 1) {
//    fserror ("scanLong","read long",txt);
//  }
//}
// 
///* --------------------   read arrays   ------------------------- */
// 
//void scanFloatArray(char *txt, float *x, int n)
//{
// 	scanArray(txt,x,n);
//}
// 
//void scanArray(char *txt, float *x, int n)
//{
// 	int	i; 
// 
//  fscanf(ifile,txt);
//  for(i=0;i<n;i++){
//    if (fscanf(ifile," %f ",&x[i]) != 1) {
//      fserror ("scanArray","read float array",txt);
//    }
//  }
//}
// 
//void scanDoubleArray(char *txt, double *x, int n)
//{
//  int	i;
// 
//  fscanf(ifile,txt);
//  for(i=0;i<n;i++){
//    if (fscanf(ifile," %lg ",&x[i]) != 1) {
//      fserror ("scanDoubleArray",
// 	       "read double array",txt);
//    }
//  }
//}
// 
//void fscanDoubleArray(FILE *in, double *x, int n)
//{
//  int	i;
// 
//  for(i=0;i<n;i++){
//    if (fscanf(in," %lg ",&x[i]) != 1) {
//      /* printf("i=%d\n",i); */
//      fserror("fscanDoubleArray","read double array","");
//    }
//  }
//}
// 
//void scanString(char *txt, char *s, int n)
//{
//  fgets(s,n,ifile);
//}
//  
// 
//void fscanString(FILE *ifile, char *txt, char *s, int n)
//{
//  fgets(s,n,ifile);
//}
//  
//void scanDoubleMatrix(char *txt,double **x,int r,int c)
//{
//  int	i,j;
//  
//  fscanf(ifile,txt);
//  for(i=0;i<r;i++)
//    for(j=0;j<c;j++){
//      if (fscanf(ifile," %lg ",&x[i][j]) != 1) {
// 	fserror ("scanDoubleMatrix","read double matrix",txt);
//      }
//    }
//}
// 
//void fscanDoubleMatrix(FILE *ifile, double **x,int r,int c)
//{
//  int	i,j;
//  
//  for(i=0;i<r;i++)
//    for(j=0;j<c;j++){
//      if (fscanf(ifile," %lg ",&x[i][j]) != 1) {
// 	exit(1);
//      }
//    }
//}
// 
//void scanIntArray(char *txt, int *x, int n)
//{
//  int	i;
// 
//  fscanf(ifile,txt);
//  for(i=0;i<n;i++){
//    if (fscanf(ifile," %d ",&x[i]) != 1) {
//      fserror ("scanIntArray","read int array",txt);
//    }
//  }
//}
// 
//void fscanIntArray(FILE *ifile, int *x, int n)
//{
//  int	i;
// 
//  for(i=0;i<n;i++){
//    if (fscanf(ifile," %d ",&x[i]) != 1) {
//      fserror("fscanIntArray","read int array","");
//    }
//  }
//}


/* ------------------------  write scalars  ------------------------ */
void writeInt(int i)
{
  int s;
  s=fprintf(ofile,"%d\n",i);
  if(s<0)
    fserror("writeInt","write int","");
}

void writeLong(long i)
{
  int s;
  s=fprintf(ofile,"%ld\n",i);
  if (s<0)
    fserror("writeLong","write long","");
    
}


void writeFloat(float x)
{
  int s;
  s=fprintf(ofile,"%f\n",x);
  if (s<0)
    fserror("writeFloat","write float","");
  
}

void writeDouble(double x)
{
  int s;
  s=fprintf(ofile,"%5.3e\n",x);
  if (s<0)
    fserror("writeDouble","write double","");
}


/* -----------------------  write arrays   --------------------- */

void fwriteDoubleArray(FILE *f, double *x, int rows, int cols)
{
  int	i,j,s1,s2;
  
  s1 = 0;
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      if(j%10 == 9)
	fprintf(f,"\n\t");
      s1=fprintf(f,"%5.3e ",x[i*cols+j]);
      if (s1<0) break;
    }
    s2=fprintf(f,"\n");
    if ((s2<0)|(s1<0))
      fserror("fwriteDoubleArray","write double array","");
  }
}

void fwriteIntArray(FILE *f, int *x, int rows, int cols)
{
  int	i,j,s1,s2;
  
  s1 = 0;
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      if(j%10 == 9)
	fprintf(f,"\n\t");
      s1=fprintf(f,"%d\t",x[i*cols+j]);
      if (s1<0) break;
    }
    s2=fprintf(f,"\n");
    if ((s2<0)|(s1<0))
      fserror("fwriteIntArray","write int array","");
  }
}


void writeIntArray(int *x, int rows, int cols)
{
  fwriteIntArray(ofile,x,rows,cols);
}

void fwriteIntMatrix(FILE *f, int **x, int rows, int cols)
{
  int	i,j,s1;
  
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      if(j%10 == 9)
	fprintf(f,"\n\t");
      s1=fprintf(f,"%d\t",x[i][j]);
      if (s1<0) 
	fserror("fwriteIntMatrix","write int matrix","");
    }
    fprintf(f,"\n");
  }
}

void writeIntMatrix(int **x, int rows, int cols)
{
  fwriteIntMatrix(ofile,x,rows,cols);
}

void writeDoubleArray(double *x,int rows,int cols)
{
  fwriteDoubleArray(ofile,x,rows,cols);
}

void fwriteDoubleMatrix2(FILE *f, double **x, int rows, int cols)
{
  int	i, j, s;
  
  for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
      if(j%10 == 9){
	fprintf(f,"\n\t");
      }
      s=fprintf(f,"%5.3e ",x[i][j]);
      if (s<0)
	fserror("fwriteDoubleMatrix2","write double matrix","");
    }
    fprintf(f,"\n");
  }
}

void writeDoubleMatrix2(double **x, int rows, int cols)
{
  fwriteDoubleMatrix2(ofile,x,rows,cols);
}

void writeDoubleMatrix(double **x, int rows, int cols)
{
  int	i,j, c,s;
  
  for(i=0;i<rows;i++){
    for(j=0, c=0;j<cols;j++){
      if(++c > 10){
	fprintf(ofile,"\n\t");
	c = 0;
      }
      s=fprintf(ofile,"%5.3e ",x[i][j]);
      if (s<0)
	fserror("fwriteDoubleMatrix","write double matrix","");
    }
    fprintf(ofile,"\n");
  }
}

void writeFloatArray(float *x,int rows,int cols)
{
 writeArray(x,rows,cols);
}

void writeArray(float *x,int rows,int cols)
{
  int	
    i,j, c,s;
  
  for(i=0;i<rows;i++){
    for(j=0,c=0;j<cols;j++){
      if (c++>9){
	fprintf(ofile,"\n\t");
	c = 0;
      }
      s=fprintf(ofile,"%5.3e ",x[i*cols+j]);
      if (s<0)
	fserror("fwriteDoubleMatrix","write double matrix","");
    }
    fprintf(ofile,"\n");
  }
}

/***************************   error   *************************** */ 

void fserror(char *proc, char *act, char *what){
  /* writes error message and aborts */
  
  fprintf(stderr, "\n ** Error ");
  if (proc[0]!='\0') /* not empty */
    fprintf(stderr, " in function '%s', ", proc);
  if (act[0]!='\0') /* not empty */
    fprintf(stderr, " trying to %s", act);
  if (what[0]!='\0') /* not empty */
    fprintf(stderr, " '%s'", what);
  fprintf(stderr, "\n ** .. exiting program");
  fprintf(stderr, " (from a function in 'interface.c').\n");
  
  exit(1);
}

/******************************************************************************************
                                      DEBUG MESSAGES ETC.
******************************************************************************************/

void errorC(char *module, char *mess, int nr)              
{
   	printf("\n *** ERROR # %d in %s***\n %s",nr,module, mess);
	printf(  " exiting program \n");
	exit(1);
}

void err_msg(char *fct, char *txt, int n1, int n2, int n3)
{ /* print error message */
  printf("\n\n *** Error in %s \n", fct);
  printf(txt,n1,n2,n3); /* n1,n2 allows to include numbers in txt */
  printf("\n");
  exit(1);
}



/******************************************************************************************
                          MEMORY ALLOCATION
******************************************************************************************/

float *vector(int nl,int nh) 
{ 
        float *v; 
 
        v=(float *)calloc((unsigned) (nh-nl+1),sizeof(float)); 
        if (!v) nrerror("vector","allocate a float vector",""); 
        return v-nl; 
} 

double  *dvector(int nl,int nh) 
{ 
        double  *v; 
 
	nv += (nh-nl+1); 
        v=(double  *)calloc((unsigned) (nh-nl+1),sizeof(double)); 
        if (!v)  
	  nrerror("dvector","allocate a double vector",""); 
        return v-nl; 
} 

double  **dmatrix(int nrl,int nrh,int ncl,int nch) 
{ 
        int i; 
        double  **m; 
	 
	nv += (nrh-nrl+1)*(nch-ncl+1); 
        m=(double  **)calloc((unsigned) (nrh-nrl+1), 
				   sizeof(double  *)); 
        if (!m)  
	  nrerror("dmatrix","allocate a double matrix (1st dim)",""); 
        m -= nrl; 
 
        for(i=nrl;i<=nrh;i++) { 
                m[i]=(double  *)calloc((unsigned) (nch-ncl+1), 
					     sizeof(double)); 
                if (!m[i])  
	  	  nrerror("dmatrix","allocate a double matrix (2nd dim)",""); 
                m[i] -= ncl; 
        } 
        return m; 
} 

/*
double ***darray_3(int lo, int hi) 
{ 
  double ***m; 
 
  nv += (hi-lo+1); 
  m=(double ***)calloc((unsigned) (hi-lo+1),sizeof(double **)); 
  if (!m)  
    nrerror("darray_3","allocate a 3dim double array ",""); 
  m -= lo; 
  return m; 
} 


double ***darray3(int n, int p, int q) 
{ 
  double ***a; 
  int i; 
  a = darray_3(0,n); 
  for(i=0;i<n;i++) 
    a[i] = dmatrix(0,p,0,q); 
  return(a); 
} 
*/

double ***darray3(int n1,int n2,int n3)
/***********************************************************************
  allocates space for a 3 index array 0..n1-1, 0...n2-1, 0...n3-1
***********************************************************************/
{
  double ***a;
  int  i, j;

  a = (double ***) malloc(n1 * sizeof(double **));
  if(a == NULL) nrerror("darray3","allocate a 3dim double array ","");

  a[0] = (double **) malloc(n1 * n2 * sizeof(double *));
  if(a[0] == NULL) nrerror("darray3","allocate a 3dim double array ","");
  for(i=1;i<n1;i++) a[i] = a[i-1] + n2;

  a[0][0] = (double *) malloc(n1 * n2 * n3 * sizeof(double));
  if(a[0][0] == NULL) nrerror("darray3","allocate a 3dim double array ","");
  for(i=0;i<n1;i++) 
    for(j=0;j<n2;j++) 
      a[i][j] = a[0][0] + n2*n3*i + j*n3;
    
  return a;
}

int  *ivector(int nl,int nh) 
{ 
        int  *v; 
 
	nv += (nh-nl+1); 
        v=(int  *)calloc((unsigned) (nh-nl+1),sizeof(int)); 
        if (!v) nrerror("ivector","allocate an int vector",""); 
        return v-nl; 
} 

int  **imatrix(int nrl,int nrh,int ncl,int nch) 
{ 
        int i, **m; 
 
	nv += (nrh-nrl+1)*(nch-ncl+1); 
        m=(int  **)calloc((unsigned) (nrh-nrl+1),sizeof(int  *)); 
        if (!m)  
	  nrerror("imatrix","allocate a int matrix (1st dim).",""); 
        m -= nrl; 
 
        for(i=nrl;i<=nrh;i++) { 
                m[i]=(int  *)calloc((unsigned) (nch-ncl+1),sizeof(int)); 
                if (!m[i])  
	  	  nrerror("imatrix","allocate a int matrix (2nd dim).",""); 
                m[i] -= ncl; 
        } 
        return m; 
} 

/*
int ***iarray_3(int lo, int hi) 
{ 
  int ***m; 
 
  nv += (hi-lo+1); 
  m=(int ***)calloc((unsigned) (hi-lo+1),sizeof(int **)); 
  if (!m)  
    nrerror("iarray_3","allocate a 3dim int array ",""); 
  m -= lo; 
  return m; 
}

int ***iarray3(int p1, int p2, int p3) 
{ 
  int ***m, i; 
 
  m = iarray_3(0,p1); 
  for (i=0;i<p1;i++) 
    m[i] = imatrix(0,p2,0,p3); 
  return m; 
} 
*/


int ***iarray3(int n1,int n2,int n3)
/***********************************************************************
  allocates space for a 3 index array 0..n1-1, 0...n2-1, 0...n3-1
***********************************************************************/
{
  int ***a, i, j;

  a = (int ***) malloc(n1 * sizeof(int **));
  if(a == NULL)  nrerror("iarray3","allocate a 3dim double array ","");

  a[0] = (int **) malloc(n1 * n2 * sizeof(int *));
  if(a[0] == NULL)  nrerror("iarray3","allocate a 3dim double array ","");
  for(i=1;i<n1;i++) a[i] = a[i-1] + n2;

  a[0][0] = (int *) malloc(n1 * n2 * n3 * sizeof(int));
  if(a[0][0] == NULL)  nrerror("iarray3","allocate a 3dim double array ","");
  for(i=0;i<n1;i++) 
    for(j=0;j<n2;j++) 
      a[i][j] = a[0][0] + n2*n3*i + j*n3;
    
  return a;
}


void free_vector(float  *v,int nl,int nh) 
{ 
        if( (v+nl) != NULL ) free((char  *) (v+nl)); 
	nv -= (nh-nl+1); 
} 

void free_dvector(double  *v,int nl,int nh) 
{ 
        if( (v+nl) != NULL ) free((char  *) (v+nl)); 
	nv -= (nh-nl+1); 
} 

void free_ivector(int  *v,int nl,int nh) 
{ 
        if( (v+nl) != NULL ) free((char  *) (v+nl)); 
	nv -= (nh-nl+1); 
} 

void free_dmatrix(double  **m,int nrl,int nrh,int ncl,int nch) 
{ 
        int i; 
 
        for(i=nrh;i>=nrl;i--) {if( (m[i]+ncl) != NULL )  
				 free((char  *) (m[i]+ncl));} 
        if( (m+nrl) != NULL ) free((char  *) (m+nrl)); 
        nv -= (nch-ncl+1)*(nrh-nrl+1); 
} 
 
void free_imatrix(int  **m,int nrl,int nrh,int ncl,int nch) 
{ 
        int i; 
 
        for(i=nrh;i>=nrl;i--) {if( (m[i]+ncl) != NULL )  
				 free((char  *) (m[i]+ncl));} 
        if( (m+nrl) != NULL ) free((char  *) (m+nrl)); 
        nv -= (nch-ncl+1)*(nrh-nrl+1); 
} 

void free_darray3(double ***a, int n1, int n2, int n3) {
  free((char*) (a[0][0]));
  free((char*) (a[0]));
  free((char*) (a));
}

void free_iarray3(int ***a, int n1, int n2, int n3) {
        free((char*) (a[0][0]));
        free((char*) (a[0]));
        free((char*) (a));
}


void nrerror(char *proc, char *act, char *what) 
{ 
  void exit(); 
 
  fprintf(stderr, "\n ** Error "); 
  if (proc[0]!='\0') /* not empty */ 
    fprintf(stderr, " in function '%s', ", proc); 
  if (act[0]!='\0') /* not empty */ 
    fprintf(stderr, " trying to %s", act); 
  if (what[0]!='\0') /* not empty */ 
    fprintf(stderr, " '%s',", what); 
  else 
    fprintf(stderr, ", "); 
  fprintf(stderr, "\n ** .. exiting program.\n"); 
  fprintf(stderr, "\n ** (a function in 'nrutil.c').\n"); 
  exit(1); 
} 


/************************************************************************
                          MATHEMATICAL FUNCTIONS
************************************************************************/

double ldoublefact(double x) {
  //log-double factorial
  //Returns x*(x-2)*...*k, where k is the first term <2. Corresponds to double factorial for integer x
  int _i; double ans=0;
  for (_i=x; _i>=2; _i-=2) { ans+= log(_i); }
  return(ans);
}

double gamln(double *a)
/*
-----------------------------------------------------------------------
            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
-----------------------------------------------------------------------
     WRITTEN BY ALFRED H. MORRIS
          NAVAL SURFACE WARFARE CENTER
          DAHLGREN, VIRGINIA
--------------------------
     D = 0.5*(LN(2*PI) - 1)
--------------------------
*/
{
static double c0 = .833333333333333e-01;
static double c1 = -.277777777760991e-02;
static double c2 = .793650666825390e-03;
static double c3 = -.595202931351870e-03;
static double c4 = .837308034031215e-03;
static double c5 = -.165322962780713e-02;
static double d = .418938533204673e0;
static double gamln,t,w;
static int i,n;
static double T1;
/*
     ..
     .. Executable Statements ..
*/
    if(*a > 0.8e0) goto S10;
    gamln = gamln1(a)-log(*a);
    return gamln;
S10:
    if(*a > 2.25e0) goto S20;
    t = *a-0.5e0-0.5e0;
    gamln = gamln1(&t);
    return gamln;
S20:
    if(*a >= 10.0e0) goto S40;
    n = (long)(*a - 1.25e0);
    t = *a;
    w = 1.0e0;
    for(i=1; i<=n; i++) {
        t -= 1.0e0;
        w = t*w;
    }
    T1 = t-1.0e0;
    gamln = gamln1(&T1)+log(w);
    return gamln;
S40:
    t = pow(1.0e0/ *a,2.0);
    w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/ *a;
    gamln = d+w+(*a-0.5e0)*(log(*a)-1.0e0);
    return gamln;
}

double gamln1(double *a)
/*
-----------------------------------------------------------------------
     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
-----------------------------------------------------------------------
*/
{
static double p0 = .577215664901533e+00;
static double p1 = .844203922187225e+00;
static double p2 = -.168860593646662e+00;
static double p3 = -.780427615533591e+00;
static double p4 = -.402055799310489e+00;
static double p5 = -.673562214325671e-01;
static double p6 = -.271935708322958e-02;
static double q1 = .288743195473681e+01;
static double q2 = .312755088914843e+01;
static double q3 = .156875193295039e+01;
static double q4 = .361951990101499e+00;
static double q5 = .325038868253937e-01;
static double q6 = .667465618796164e-03;
static double r0 = .422784335098467e+00;
static double r1 = .848044614534529e+00;
static double r2 = .565221050691933e+00;
static double r3 = .156513060486551e+00;
static double r4 = .170502484022650e-01;
static double r5 = .497958207639485e-03;
static double s1 = .124313399877507e+01;
static double s2 = .548042109832463e+00;
static double s3 = .101552187439830e+00;
static double s4 = .713309612391000e-02;
static double s5 = .116165475989616e-03;
static double gamln1,w,x;
/*
     ..
     .. Executable Statements ..
*/
    if(*a >= 0.6e0) goto S10;
    w = ((((((p6**a+p5)**a+p4)**a+p3)**a+p2)**a+p1)**a+p0)/((((((q6**a+q5)**a+
      q4)**a+q3)**a+q2)**a+q1)**a+1.0e0);
    gamln1 = -(*a*w);
    return gamln1;
S10:
    x = *a-0.5e0-0.5e0;
    w = (((((r5*x+r4)*x+r3)*x+r2)*x+r1)*x+r0)/(((((s5*x+s4)*x+s3)*x+s2)*x+s1)*x
      +1.0e0);
    gamln1 = x*w;
    return gamln1;
}



// digamma(x) returns the derivative of the natural log of the gamma function, e.g. gamma'(x)/gamma(x)
//  x must be positive and finite

double digamma(double x) {

  double stirling[] = {
    -8.333333333333333e-02, 8.333333333333333e-03,
    -3.968253968253968e-03, 4.166666666666667e-03,
    -7.575757575757576e-03, 2.109279609279609e-02,
    -8.333333333333334e-02, 4.432598039215686e-01,
    -3.053954330270120e+00, 2.645621212121212e+01,
    -2.814601449275362e+02, 3.607510546398047e+03 
  };

  long i;
  double lower= 1.0e-8, upper= 19.5, euler_one= .422784335098467139393488, ans, x_inv, x_pow;


  if (x<=0) errorC("digamma",  "digamma must be given a positive argument", 1);;

  if (x<lower) {
    ans = -1.0 / x - 1.0/(1.0+x) + euler_one;
    return(ans);
  }

  ans= 0.0;
  while (x<upper) {
    ans= ans - 1.0/x;
    x= x + 1.0;
  }

  x_inv= 1.0/x;
  ans= ans + log(x) - 0.5*x_inv;

  x_inv= x_inv*x_inv;
  x_pow= x_inv;

  for (i=0; i<12; i++) {
    ans= ans + stirling[i] * x_pow;
    x_pow= x_pow * x_inv;
  }
  return(ans);

}


/* Bernoulli numbers of even order from 2 to 60 */
static double bernou[30] = {1.0/6.0, -1.0/30.0, 1.0/42.0, -1.0/30.0, 5.0/66.0, -691.0/2730.0, 7.0/6.0, -3617.0/510.0, 43867.0/798.0,
-174611.0/330.0, 854513.0/138.0, -236364091.0/2730.0, 8553103.0/6.0, -23749461029.0/870.0, 8615841276005.0/14322.0,
-7709321041217.0/510.0, 2577687858367.0/6.0, -1.371165521e13, 4.883323190e14, -1.929657934e16,
8.416930476e17, -4.033807185e19, 2.115074864e21, -1.208662652e23, 7.500866746e24, -5.038778101e26,
3.652877648e28, -2.849876930e30, 2.386542750e32, -2.139994926e34};

double trigamma(double x) { 

if (x>1.0e-5) {                         //fast approx to trigamma
  return(1/(x*x) + 1/((x+1)*(x+1)) + 1/((x+2)*(x+2)) + 1/(x+3) + .5/((x+3)*(x+3)) + 1/(6.0*pow(x+3,3)));
} else {
  return(polygamma(x,1,.0001,100,5,1)); //slower computation
}
 
} 

double polygamma(double x, long n, double low, double high, long terms, double nfact)
//     long n, terms;
//     double x, low, high, nfact;
{
/* polygamma function of a real positive x of order n (n==1 is trigamma etc.).
 * setting low=0.0001, high=100 and terms=5 usually gives good results
 * nfact is n! e.g. 1 for trigamma etc.
*  no checks are made here on the suitability of arguments */
long i;
double asign, ans = 0.0, nd = (double) n, nexp, ser = 0.0;
double t0, x2_inv;

 asign = (n % 2) ? 1.0 : -1.0;
 if(x < low) { return(asign * nfact / nd * pow(x, - nd) * (1.0 + nd * .5 / x)); }
 nexp = - nd - 1.0;
 while(x < high) {
   ans = ans + asign * nfact * pow(x, nexp);
   x = x + 1.0;
 }
 t0 = nfact / nd * pow(x, - nd);
 ser = t0 * ( 1.0 + nd * .5 / x);
 x2_inv = pow(x, -2.0);
 for(i=0; i < terms; i++) {
   if(n ==1) {
     t0 = t0 * x2_inv;
   } else {
     t0 = (2.0 * i + nd + 3.0) / (2.0 * i + 4.0) * (2.0 * i + nd + 2.0) / (2.0 * i + 3.0) * t0 * x2_inv;
   }
   ser = ser + bernou[i] * t0;
 }
 ans = ans + asign * ser;
 return(ans);
}


/* log of Beta function */
double lnbeta(double a, double b) {
  double c;
  c= a+b;
  return(gamln(&a)+gamln(&b)-gamln(&c));
}


double betacf(double a, double b, double x) {
//Used by pbetaC: Evaluates continued fraction for incomplete beta function by modified Lentz's
//method (x5.2).
  int m,m2, MAXIT=100;
  double aa,c,d,del,h,qab,qam,qap, EPS=3.0e-7, FPMIN=1.0e-30;
  qab=a+b; //These q's will be used in factors that occur in the coe cients (6.4.6).
  qap=a+1.0; 
  qam=a-1.0;
  c=1.0; //First step of Lentz's method.
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d; //One step (the even one) of the recurrence.
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d; //Next step of the recurrence (the odd one).
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break; //Are we done?
  }
  if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf","","");
  return(h);
}

double lnchoose(int n, int k) {
  double a= 1.0+n-k, b= 1.0+k;
  return lnbeta(a,b) - log(1.0+n);
}

double choose(int n, int k) { return exp(lnchoose(n,k)); }

double logit(double x) { return(log(x/(1-x))); }

double ilogit(double x) { return(1.0/(1.0+exp(-x))); }

double dsign(double x) {
//Returns 1.0 if x>=0, -1.0 if x<0
  double ans;
  if (x>=0) { 
    ans= 1.0;
  } else {
    ans= -1.0;
  }
  return ans;
}

double isign(int x) {
//Returns 1.0 if x>0, 0 if x==0, -1.0 if x<0
  double ans;
  if (x==0) {
    ans= 0.0;
  } else if (x>0) {
    ans= 1.0;
  } else {
    ans= -1.0;
  }
  return ans;
}


/************************************************************************
                            VECTOR ALGEBRA
************************************************************************/

void grid (double x0, double xn, int n, double *x) 
{ 
  int i; 
  double dx,xi; 
  dx = (xn-x0)/(n-1.0); 
  for(i=0,xi=x0;i<n;i++,xi+=dx) 
    x[i]=xi; 
} 

void rA(double r,double **A, double **B, int rowini, int rowfi, int colini, int colfi) {
  //Multiply matrix A[1..p][1..q] by scalar r, store results in matrix B
  int _i, _j;
  for(_i=rowini; _i<=(rowfi); _i++){			 
    for(_j=colini; _j<=(colfi); _j++)			 
      B[_i][_j] = r* A[_i][_j];	 
  } 
} 

void A_plus_B(double **A, double **B, double **C, int rowini, int rowfi, int colini, int colfi) {
 //Sum matrix A[rowini..rowfi][colini..colfi] + B[rowini..rowfi][colini..colfi], store results in C
  int _i, _j;
  for (_i=rowini; _i<=rowfi; _i++) { for (_j=colini; _j<=colfi; _j++) { C[_i][_j]= A[_i][_j] + B[_i][_j]; } }
}


void  rA_plus_sB(double r, double **A, double s, double **B, double **C, int rowini, int rowfi, int colini, int colfi) {
  //Multiply scalar r times matrix A, add scalar s times matrix B, store results in C
  int _i, _j;
  for(_i=rowini;_i<=(rowfi);_i++)  
    for(_j=colini;_j<=(colfi);_j++)			 
      C[_i][_j]= r*A[_i][_j]+s*B[_i][_j];  
} 

void rAx_plus_sBy(double r, double **A, double *x, double s, double **B, double *y, double *z, int rowini, int rowfi, int colini, int colfi) {
  //Scalar*matrix*vector + scalar*matrix*vector
  int _i, _j; 
  for(_i=rowini;_i<=rowfi;_i++) 
    for(z[_i]=0,_j=colini; _j<=rowfi; _j++) 
      z[_i] += r*A[_i][_j]*x[_j] + s*B[_i][_j]*y[_j]; 
} 

void Ax_plus_y(double **A, double *x, double *y, double *z, int ini, int fi) { 
  //Multiply matrix A[ini..fi][ini..fi] by vector x[ini..fi] and add vector y[ini..fi]
  //Store result in vector z
  int _i,_j;
  for(_i=ini;_i<=fi;_i++) 
    for(z[_i]=y[_i],_j=ini; _j<=fi; _j++) 
      z[_i] += A[_i][_j]*x[_j]; 
} 

void xA(double *x,double **A,double *z, int ini, int fi) { 
  int _i, _j;
  for(_i=ini;_i<=(fi);_i++){				 
    for(z[_i]=0,_j=ini; _j<=(fi); _j++)		 
      z[_i]+=A[_j][_i]*x[_j];	 
  } 
}
 
void Ax(double **A,double *x,double *z, int rowini, int rowfi, int colini, int colfi) { 
  int _i, _j;
  for(_i=rowini;_i<=rowfi;_i++){				 
    for(z[_i]=0,_j=colini; _j<=colfi; _j++)		 
      z[_i]+=A[_i][_j]*x[_j];	 
  } 
} 

double xtAy (double *x, double **A, double *y, int ini, int fi) { 
  int _i, _j; double z; 
  for(z=0,_i=ini;_i<=fi;_i++) 
    for(_j=ini; _j<=fi; _j++) 
      z += A[_i][_j]*x[_j]*y[_i]; 
  return(z); 
}
 
double quadratic_xtAx(double *x, double **A, int ini, int fi) {
 //t(vector)*matrix*vector for quadratic forms (A must be symmetric)
 //Note: this routine is faster than xtAy for symmetric A (saves 25%-50% operations)
  int _i, _j; double z;
  for (z=0,_i=ini; _i<=fi; _i++) {
    z+= A[_i][_i]*x[_i]*x[_i];
    for (_j=_i+1; _j<=fi; _j++) {
      z+= 2*A[_i][_j]*x[_i]*x[_j];
    }
  }
  return(z);
}

double quadratic_xseltAselxsel(double *x, double *A, int *ncolA, int *nsel, int *sel) {
 //t(x[sel])*A[sel,sel]*x[sel] for quadratic forms (A must be symmetric and given as a vector)
 // - ncolA: number of columns in A
 // - nsel: length of vector sel
 // - sel: vector with indexes for positions in x and (rows,columns) in A to be used in the operation
 //Note: this routine is faster than xtAy for symmetric A (saves 25%-50% operations)
  int _i, _j; double z;
  for (z=0,_i=0; _i<=(*nsel)-1; _i++) {
    z+= A[sel[_i]*(*ncolA)+sel[_i]]*x[sel[_i]]*x[sel[_i]];
    for (_j=_i+1; _j<=(*nsel)-1; _j++) {
      z+= 2*A[sel[_i]*(*ncolA)+sel[_j]]*x[sel[_i]]*x[sel[_j]];
    }
  }
  return(z);
}

double quadratic_xtAselx(double *x, double *A, int *ncolA, int *nsel, int *sel) {
 //t(x)*A[sel,sel]*x for quadratic forms (A must be symmetric and given as a vector)
 // - ncolA: number of columns in A
 // - nsel: length of vector sel
 // - sel: vector with indexes for (rows,columns) in A to be used in the operation
 //Note: this routine is faster than xtAy for symmetric A (saves 25%-50% operations)
  int _i, _j; double z;
  for (z=0,_i=0; _i<=(*nsel)-1; _i++) {
    z+= A[sel[_i]*(*ncolA)+sel[_i]]*x[_i]*x[_i];
    for (_j=_i+1; _j<=(*nsel)-1; _j++) {
      z+= 2*A[sel[_i]*(*ncolA)+sel[_j]]*x[_i]*x[_j];
    }
  }
  return(z);
}


double quadratic_xseltAxsel(double *x, double **A, int ini, int *nsel, int *sel) {
 //t(x[sel])*A*x[sel] for quadratic forms (A must be symmetric and given as a matrix)
 // - ini: element in A are indexed from ini to *nsel - ini
 // - nsel: length of vector sel
 // - sel: vector with indexes for positions in x and (rows,columns) in A to be used in the operation
 //Note: this routine is faster than xtAy for symmetric A (saves 25%-50% operations)
  int _i, _j, rowA; double z;
  for (z=0,_i=0; _i<=(*nsel)-1; _i++) {
    rowA= ini + _i;
    z+= A[rowA][rowA]*x[sel[_i]]*x[sel[_i]];
    for (_j=_i+1; _j<=(*nsel)-1; _j++) {
      z+= 2*A[rowA][ini+ _j]*x[sel[_i]]*x[sel[_j]];
    }
  }
  return(z);
}


void Atx(double **A,double *x,double *z, int rowini, int rowfi, int colini, int colfi) {
  int _i, _j; 
  for(_i=colini;_i<=colfi;_i++){				 
    for(z[_i]=0,_j=rowini; _j<=rowfi; _j++)		 
      z[_i]+=A[_j][_i]*x[_j];	 
  } 
} 

void AtB(double **A, int rowiniA, int rowfiA, int coliniA, int colfiA, double **B, int rowiniB, int rowfiB, int coliniB, int colfiB, double **C) { 
  int _i, _j, _k;
  if ((rowfiA-rowiniA) != (rowfiB-rowiniB)) errorC("AtB","dimensions don't match",1); 
  for(_i=coliniA;_i<=colfiA;_i++)			 
    for(_j=coliniB;_j<=colfiB;_j++)			 
      for(C[_i][_j]=0,_k=rowiniA;_k<=rowfiA;_k++)	 
	C[_i][_j]+=A[_k][_i]*B[_k][_j]; 
} 

void a_plus_b(double *a, double *b, double *c, int ini, int fi) {
  int _i;
  for (_i=ini;_i<=fi;_i++) { c[_i]= a[_i]+b[_i]; }
}

void a_prod_b(double *a, double *b, double *c, int ini, int fi) {
  int _i;
  for (_i=ini;_i<=fi;_i++) { c[_i]= a[_i]*b[_i]; }
}

void a_prod_b_sel(double *a, double *b, double *c, int *lengtha, int *nsel, int *sel) {
  int _i;
  for (_i=0;_i<=(*nsel-1);_i++) { c[sel[_i]]= a[sel[_i]]*b[sel[_i]]; }
}

void a_zero(double *a, int p) {
  int _i;
  for (_i=0;_i<p;_i++) a[_i]= 0.0;
}

void R_zero(double **A, int p,int q) 
{
  int _i, _j;
  for(_i=0;_i<p;_i++) 
    for(_j=0;_j<q;_j++) 
      A[_i][_j] = 0.0; 
} 

void ddiag(double **A, int ini, int fi) {
//Diagonal matrix
  int _i, _j;
  for (_i=ini; _i<= fi; _i++) {
    for (_j=ini; _j<= fi; _j++) {
      if (_i==_j) A[_i][_j]= 1; else A[_i][_j]= 0;
    }
  }
}

int iabs(int x) {
  return (x>0) ? x : -x;
}

int imax_xy(int x, int y) 
{ 
  return (x>y) ? x : y; 
} 

int imin_xy(int x, int y) 
{ 
  return (x<y) ? x : y; 
} 

double max_xy(double x, double y) 
{ 
  return (x>y) ? x : y; 
} 

double min_xy(double x, double y) 
{ 
  return (x<y) ? x : y; 
} 


void minvec(double *x, int ini, int fi, double *xmin, int *minpos) {
//Minimum of vector x[ini..fi] is returned in xmin. Position at which min occurs is returned in minpos
  int _i;
  *xmin= x[ini]; *minpos= ini;
  for (_i=ini+1;_i<=fi;_i++) { if (x[_i]<(*xmin)) { *xmin= x[_i]; *minpos= _i; } }
}

void maxvec(double *x, int ini, int fi, double *xmax, int *maxpos) {
//Maximum of vector x[ini..fi] is returned in xmax. Position at which max occurs is returned in maxpos
  int _i;
  *xmax= x[ini]; *maxpos= ini;
  for (_i=ini+1;_i<=fi;_i++) { if (x[_i]>(*xmax)) { *xmax= x[_i]; *maxpos= _i; } }
}


void choldc(double **a, int n, double **aout) {
/*Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
decomposition, A = L * L' . On input, only the upper triangle of a need be given; 
 The Cholesky factor L is returned in the lower triangle of aout (upper-diag elem are set to 0) */
  int i,j,k;
  double sum, *p;

  for (i=1;i<=n;i++) { for (j=i;j<=n;j++) { aout[i][j]= a[i][j]; } }  //copy a into aout
  p= dvector(1,n);
  for (i=1;i<=n;i++) {
    for (j=i;j<=n;j++) {
      for (sum=aout[i][j],k=i-1;k>=1;k--) sum -= aout[i][k]*aout[j][k];
      if (i == j) {
	if (sum <= 0.0) nrerror("choldc failed","","matrix is not positive definite");
	aout[i][i]=sqrt(sum);
      } else aout[j][i]=sum/aout[i][i];
    }
  }
  free_dvector(p,1,n);
  for (i=1;i<=n;i++) { for (j=i+1;j<=n;j++) { aout[i][j]= 0; } }  //set upper-diagonal elem to 0
}

void choldc_inv(double **a, int n, double **aout) {
  /*Given a positive-definite symmetric matrix a[1..n][1..n], this routine computes the inverse
   of its Cholesky matrix. That is, if A=L * L' it returns the inverse of L
   (note that inv(A)= inv(L)' * inv(L)) */
  int i,j,k;
  double sum;

  choldc(a,n,aout);
  for (i=1;i<=n;i++) {
    aout[i][i]=1.0/aout[i][i];
    for (j=i+1;j<=n;j++) {
      sum=0.0;
      for (k=i;k<j;k++) sum -= aout[j][k]*aout[k][i];
      aout[j][i]=sum/aout[j][j];
    }
  }
}

double choldc_det(double **chols, int n) {
  //Find determinant of the matrix having chols as its Cholesky decomposition
  //Example of usage
  //  choldc(S,n,cholS);
  //  det= choldc_det(cholS,n);
  //
  //Another example
  //  choldc_inv(S,n,cholSinv);
  //  det= 1.0/choldc_det(cholSinv,n);
  int i; double det;
  for (det=1,i=1; i<=n; i++) { det *= chols[i][i]*chols[i][i]; }
  return(det);
}


void inv_posdef(double **a, int n, double **aout) {
  /* Inverse of a symmetric, positive definite matrix a[1..n][1..n] using Cholesky decomposition
     Result is returned in aout */
  int i,j,k;
  double **b, sum;

  b= dmatrix(1,n,1,n);
  choldc_inv(a,n,b);
  for (i=1; i<=n; i++) {
    for (j=i; j<=n; j++) {
      for (k=1,sum=0;k<=n;k++) { sum+= b[k][i]*b[k][j]; }
      aout[i][j]= sum;
    }
  }
  for (i=2; i<=n; i++) { for (j=1;j<i;j++) { aout[i][j]= aout[j][i]; }}
  free_dmatrix(b,1,n,1,n);

}

void inv_posdef_upper(double **a, int n, double **aout) {
  /* Inverse of a symmetric, positive definite matrix a[1..n][1..n] using Cholesky decomposition
     Result is returned in aout 
     Function does the same as inv_posdef, except that here only upper triangular elements are returned*/
  int i,j,k;
  double **b, sum;

  b= dmatrix(1,n,1,n);
  choldc_inv(a,n,b);
  for (i=1; i<=n; i++) {
    for (j=i; j<=n; j++) {
      for (k=1,sum=0;k<=n;k++) { sum+= b[k][i]*b[k][j]; }
      aout[i][j]= sum;
    }
  }
  free_dmatrix(b,1,n,1,n);

}


void invdet_posdef(double **a, int n, double **aout, double *det_a) {
  /* Inverse and determinant of a positive definite matrix a[1..n][1..n] using Cholesky decomposition
     Inverse is returned in aout, determinant in det_a */
  int i,j,k;
  double **b, sum;

  b= dmatrix(1,n,1,n);
  choldc_inv(a,n,b);
  for (i=1, *det_a=1; i<=n; i++) { (*det_a)*= 1/(b[i][i]*b[i][i]); }
  for (i=1; i<=n; i++) {
    for (j=i; j<=n; j++) {
      for (k=1,sum=0;k<=n;k++) { sum+= b[k][i]*b[k][j]; }
      aout[i][j]= sum;
    }
  }
  for (i=2; i<=n; i++) { for (j=1;j<i;j++) { aout[i][j]= aout[j][i]; }}
  free_dmatrix(b,1,n,1,n);

}


void inv_posdef_chol(double **invchol, int n, double **aout) {
  /* Inverse of a positive definite matrix with inverse of Cholesky decomposition stored in invchol
     Result is returned in aout */
  // Example of usage
  //   choldc_inv(a,n,invchol);
  //   inv_posdef_chol(invchol,n,ainv);
  int i,j,k;
  double sum;

  for (i=1; i<=n; i++) {
    for (j=i; j<=n; j++) {
      for (k=1,sum=0;k<=n;k++) { sum+= invchol[k][i]*invchol[k][j]; }
      aout[i][j]= sum;
    }
  }
  for (i=2; i<=n; i++) { for (j=1;j<i;j++) { aout[i][j]= aout[j][i]; }}

}


/* LU decomposition, Inverse and determinant of a non-singular matrix */
void ludc(double **a, int n, int *indx, double *d) {
/* Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
indx[1..n] is an output vector that records the row permutation e ected by the partial
pivoting; d is output as  1 depending on whether the number of row interchanges was even
or odd, respectively. This routine is used in combination with lu_solve to solve linear equations
or invert a matrix. */

//  int i,imax,j,k; //initialized imax to 1 to avoid warning when compiling
  int i,imax=1,j,k;
  double big,dum,sum,temp,TINY=1.0e-20;
  double *vv; //vv stores the implicit scaling of each row.
  vv=dvector(1,n);
  *d=1.0; //No row interchanges yet.
  for (i=1;i<=n;i++) { //Loop over rows to get the implicit scaling information
    big= 0.0;
    for (j=1;j<=n;j++) if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp","",""); //No nonzero largest element.
    vv[i]=1.0/big; //Save the scaling.
  }
  for (j=1;j<=n;j++) { //This is the loop over columns of Crout's method.
    for (i=1;i<j;i++) { //This is equation (2.3.12) except for i = j.
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0; //Initialize for the search for largest pivot element.
    for (i=j;i<=n;i++) { //This is i = j of equation (2.3.12) and i = j+1 : ::N of equation (2.3.13).
      sum=a[i][j]; 
      for (k=1;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) { //Is gure of merit for the pivot better than the best so far?
	big=dum;
	imax=i;
      }
    }
    if (j != imax) { //Do we need to interchange rows?
      for (k=1;k<=n;k++) { ///Yes, do so...
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d); //...and change the parity of d.
      vv[imax]=vv[j]; //Also interchange the scale factor.
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    /*If the pivot element is zero the matrix is singular (at least to the precision of the
      algorithm). For some applications on singular matrices, it is desirable to substitute
      TINY for zero. */
    if (j != n) { //Now,  nally, divide by the pivot element.
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  } //Go back for the next column in the reduction.
  free_dvector(vv,1,n);
}


void lu_solve(double **a, int n, int *indx, double b[]) {
/*Solves the set of n linear equations A X = B. Here a[1..n][1..n] is input, not as the matrix
A but rather as its LU decomposition, determined by the routine ludc. indx[1..n] is input
as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
B, and returns with the solution vector X. a, n, and indx are not modi ed by this routine
and can be left in place for successive calls with di erent right-hand sides b. This routine takes
into account the possibility that b will begin with many zero elements, so it is efficient for use
in matrix inversion. */
/* Usage:
     ludc(a,n,indx,&d);
     lu_solve(a,n,indx,b);  //answer in b, original matrix A has been destroyed
     lu_solve(a,n,indx,c);  //now we solve AX=c, since A hasn't changed we only call lu_solve
 */

  int i,ii=0,ip,j;
  double sum;
  for (i=1;i<=n;i++) { 
    /*When ii is set to a positive value, it will become the
      index of the  rst nonvanishing element of b. Wenow
      do the forward substitution, equation (2.3.6). The
      only new wrinkle is to unscramble the permutation as we go. */
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i; //Nonzero element encountered, from now we do sums in the loop above.
    b[i]=sum; 
  }
  for (i=n;i>=1;i--) { //Now we do the backsubstitution, equation (2.3.7).
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i]; //Store a component of the solution vector X.
  } //All done!
}

void lu_inverse(double **a, int n, double **aout) {
/* Inverse of a non-singular matrix a. Result is stored in aout, original matrix a is destroyed */
double d,*col;
int i,j,*indx;

 indx= ivector(1,n); col= dvector(1,n);
 ludc(a,n,indx,&d); //Decompose the matrix just once.
 for(j=1;j<=n;j++) { //Find inverse by columns.
   for(i=1;i<=n;i++) col[i]=0.0;
   col[j]=1.0;
   lu_solve(a,n,indx,col);
   for(i=1;i<=n;i++) aout[i][j]=col[i];
 }
 free_ivector(indx,1,n); free_dvector(col,1,n);

}

double lu_det(double **a, int n) {
/*Determinant of a matrix a. Original matrix a is destroyed and its LU decomposition is returned */
  double d;
  int j,*indx;

  indx= ivector(1,n);
  ludc(a,n,indx,&d); //This returns d as +/-1.
  for(j=1;j<=n;j++) d *= a[j][j];
  free_ivector(indx,1,n);
  return(d);
}



int dcompare (const void *a, const void *b) {
const double *da = (const double *) a;
const double *db = (const double *) b;
     
return (*da > *db) - (*da < *db);
}

void dvecsort(double *v, int size) {
qsort (v, size, sizeof (double), dcompare);
}


void dindexsort(double *x, int *index, int ilo, int ihi, int incr) {
/* Sorts vector of doubles x by rearranging values in index with quicksort algorithm e.g. x[index[ilo]], x[index[ilo+1]]... x[index[ihi]] is ordered */
/* Input
   - x: vector of doubles that we want to order from position ilo to ihi
   - index: vector of integers indexing the values of x
   - ilo: first element of x we want to order
   - ihi: last element of x we want to order
   - incr: for incr==1 x is returned in increasing order; incr==-1 in decreasing order
   Output: vector index rearranged so that x[index[lo]], x[index[lo+1]]... x[index[ihi]] is ordered
*/

int pivot;              // pivot value for partitioning array
int ulo, uhi;           // indices at ends of unpartitioned region
int tempEntry;          // temporary entry used for swapping
int sortlo, sortup;      // indicate if sub-vectors are sorted so no further subdivision is needed

if (ilo >= ihi) { return; }

sortlo= sortup= 1;
pivot = (ilo + ihi)/2;                                                // Select a pivot value
ulo = ilo; uhi = ihi;                                                 // Initialize ends of unpartitioned region
// While the unpartitioned region is not empty, try to reduce its size.
while (ulo < uhi) {
  if ((x[index[uhi]]*incr) > (x[index[pivot]]*incr)) {                // Here, we can reduce the size of the unpartitioned region and try again.
    if ((uhi<ihi) && ((x[index[uhi]]*incr)>(x[index[uhi+1]]*incr))) sortup= 0;      // Check if upper subvector is ordered
    uhi--;
    if ((uhi==pivot) && (ulo<pivot))  { tempEntry= index[pivot]; index[pivot]= index[pivot-1]; index[pivot-1]= tempEntry; pivot--; } 
  } else {                                                            // Here, x[index[uhi]] <= x[index[pivot]], so swap entries at indices ulo and uhi.
    tempEntry = index[ulo]; index[ulo] = index[uhi]; index[uhi] = tempEntry;
    if (pivot==ulo) pivot= uhi;
    if ((ulo>ilo) && ((x[index[ulo]]*incr)<(x[index[ulo-1]]*incr))) sortlo= 0;
    ulo++;                                                            // Reduce the size of the unpartitioned region
    if ((ulo==pivot) && (uhi>(pivot+1)))  { tempEntry= index[pivot]; index[pivot]= index[pivot+1]; index[pivot+1]= tempEntry; pivot++; }
  }
}

// Entries from ilo to pivot - 1 are < or > pivot and from pivot+1 to ihi are > or < pivot. The two regions can be sorted recursively.
if ((sortlo==0) && (ilo<(pivot-1))) dindexsort(x, index, ilo, pivot - 1, incr);
if ((sortup==0) && (ihi>(pivot+1))) dindexsort(x, index, pivot + 1, ihi, incr);

}


void iindexsort(int *x, int *index, int ilo, int ihi, int incr) {
/* Sorts vector of integers x by rearranging values in index with quicksort algorithm e.g. x[index[ilo]], x[index[ilo+1]]... x[index[ihi]] is ordered */
/* Input
   - x: vector of doubles that we want to order from position ilo to ihi
   - index: vector of integers indexing the values of x
   - ilo: first element of x we want to order
   - ihi: last element of x we want to order
   - incr: for incr==1 x is returned in increasing order; incr==-1 in decreasing order
   Output: vector index rearranged so that x[index[lo]], x[index[lo+1]]... x[index[ihi]] is ordered
*/

int pivot;              // pivot value for partitioning array
int ulo, uhi;           // indices at ends of unpartitioned region
int tempEntry;          // temporary entry used for swapping
int sortlo, sortup;      // indicate if sub-vectors are sorted so no further subdivision is needed

if (ilo >= ihi) { return; }

sortlo= sortup= 1;
pivot = (ilo + ihi)/2;                                                // Select a pivot value
ulo = ilo; uhi = ihi;                                                 // Initialize ends of unpartitioned region
// While the unpartitioned region is not empty, try to reduce its size.
while (ulo < uhi) {
  if ((x[index[uhi]]*incr) > (x[index[pivot]]*incr)) {                // Here, we can reduce the size of the unpartitioned region and try again.
    if ((uhi<ihi) && ((x[index[uhi]]*incr)>(x[index[uhi+1]]*incr))) sortup= 0;      // Check if upper subvector is ordered
    uhi--;
    if ((uhi==pivot) && (ulo<pivot))  { tempEntry= index[pivot]; index[pivot]= index[pivot-1]; index[pivot-1]= tempEntry; pivot--; } 
  } else {                                                            // Here, x[index[uhi]] <= x[index[pivot]], so swap entries at indices ulo and uhi.
    tempEntry = index[ulo]; index[ulo] = index[uhi]; index[uhi] = tempEntry;
    if (pivot==ulo) pivot= uhi;
    if ((ulo>ilo) && ((x[index[ulo]]*incr)<(x[index[ulo-1]]*incr))) sortlo= 0;
    ulo++;                                                            // Reduce the size of the unpartitioned region
    if ((ulo==pivot) && (uhi>(pivot+1)))  { tempEntry= index[pivot]; index[pivot]= index[pivot+1]; index[pivot+1]= tempEntry; pivot++; }
  }
}

// Entries from ilo to pivot - 1 are < or > pivot and from pivot+1 to ihi are > or < pivot. The two regions can be sorted recursively.
if ((sortlo==0) && (ilo<(pivot-1))) iindexsort(x, index, ilo, pivot - 1, incr);
if ((sortup==0) && (ihi>(pivot+1))) iindexsort(x, index, pivot + 1, ihi, incr);

}


/**************************************************************/
/* Random sampling                                            */
/**************************************************************/

void samplei_wr(int *x, int popsize, int n) { 
//Sample of size n without replacement from vector x of length popsize
//Result is returned in the first n elements of x
  int i, r, temp;
  for (i=0; i<n; i++) {
    r= i + (popsize - i - 1)*runif();
    temp= x[i]; x[i]= x[r]; x[r]= temp;
  }
}

void sampled_wr(double *x, int popsize, int n) { //same for vector of doubles
//Sample of size n without replacement from vector x of length popsize
//Result is returned in the first n elements of x
  int i, r;
  double temp;
  for (i=0; i<n; i++) {
    r= i + (popsize - i - 1)*runif();
    temp= x[i]; x[i]= x[r]; x[r]= temp;
  }
}


/************************************************************************
                       RANDOM VARIATE GENERATION
************************************************************************/

/* call setall(is1,is2) */ 
void setseed(long is1, long is2) 
{ 
  set=1;   
  setall(is1,is2); 
} 

double runif() 
{ 
  double  
    x; 
 
  if (set==0){ 
    setall(is1,is2); 
    set=1; 
  } 
   
  /* assign to double x for conversion */ 
  x = genunf(0.0,1.0); 
  return(x); 
} 

double dunifC(double x, double a, double b) {
  //Density of a Unif(a,b)
  return ((x>a) && (x<b)) ? (1.0/(b-a)) : 0; 
}


int runifdisc(int min, int max) {
//Returns integer value between min and max (both included)
  return(min + runif()*(max+1-min));
}

int rdisc(double *probs, int nvals) {
/* Random deviates from a discrete distribution with values 0,1...nvals-1 and probabilities probs[0],probs[1]...probs[nvals-1] */
/* Returns 0 with probability probs[0], 1 with probability probs[1] etc. */
  int i;
  double u, pcum;
  u= runif();
  pcum= probs[0]; i=1;
  while((pcum<u) && (i<nvals)) { pcum += probs[i]; i++; }
  return(i-1);
}

double rbetaC(double alpha, double beta) 
{ 
      double gamdev(); 
      double x, y; 
       
      x = gamdev(alpha);               /* X ~ gamma(alpha) */ 
      y = gamdev(beta);                /* Y ~ gamma(beta)  */ 
 
      return (x/(x+y));     /* X/(X+y) ~ ebta(alpha, beta) */ 
} 

/* CDF of a Beta distribution */
double pbetaC(double x, double a, double b) {
  double bt, c;
  if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai","","");
  if (x == 0.0 || x == 1.0) bt=0.0;
  else { //Factors in front of the continued fraction.
    c= a+b;
    bt=exp(gamln(&c)-gamln(&a)-gamln(&b)+a*log(x)+b*log(1.0-x));
  }
  if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a; //Use continued fraction directly.
  else return(1.0-bt*betacf(b,a,1.0-x)/b); //Use continued fraction after symmetry transformation.
}

void rdirichlet(double *w, double *alpha, int *p) 
/* Draws from Dirichlet with parameter alpha. 
   The value is saved in w, and p is the dimensionality of the parameter 
*/ 
{ 
  double s, a,b, W; 
  int j; 
 
  for(s=0,j=0;j<*p;j++) s += alpha[j]; 
  for(j=0,W=1.0,b=s;j<*p-1;j++){ 
    a = alpha[j]; 
    b -= alpha[j]; 
    w[j] = rbetaC(a,b)*W; 
    W -= w[j]; 
  } 
  w[*p-1] = W; 
  if (W < 0) printf("\n\n **** non-pos dirich gen!! **\n"); 
} 


double ddirichlet(double *w, double *alpha, int *p) {
 //Evaluates Dirichlet density at w (alpha: params; p: dimensionality of w and alpha)
  int i; double ans, sumalpha;

  for (i=0, ans=0, sumalpha=0; i< *p; i++) {
    ans+= (alpha[i]-1)*log(w[i]) - gamln(alpha+i);
    sumalpha+= alpha[i];
  }
  ans+= gamln(&sumalpha);
  return(exp(ans));
}

double gamdev(double alpha) 
{ 
  double  
    value; 
  double 
    a; 
  double  
    gengam(double,double); 
   
  a = alpha; /* type conversion */ 
  value = gengam(1.0,a); 
  return value; 
} 

/* *************************************************' 
   normal cdf and inv cdf 
 ************************************************* */ 
double	pnormC(double y, double m, double s) 
/* returns cdf of normal N(m,s^2) at x */ 
{	 
  double  /* used to all be float... */
    p, q, mean,sd,bound,x,z; 
  double  
    cdf; 
  int  
    status,which; 
 
  /* primitive type conversion */ 
  x = y;  
  mean = m; 
  sd = s; 
  which = 1; 
  z = (x-mean)/sd; 
 
  if (z < -5.0) 
    p = 2.86e-7F; 
  else if (z > 5.0) 
    p = 0.9999997F; 
  else 
    cdfnor(&which,&p, &q,     &x,     &mean,  &sd,    &status,    &bound); 
   
  cdf = p; /* another primitive type conversion */ 
  return cdf; 
}


double dnormC(double y, double m, double s, int logscale) {
/* Density of univariate Normal(m,s^2) evaluated at y. log==1 returns in log-scale */

  if (logscale==1) return(-log(SQ_M_PI_2)-log(s)-.5*(y-m)*(y-m)/(s*s));
  else return(exp(-.5*(y-m)*(y-m)/(s*s))/(SQ_M_PI_2 * s));

}

double dnormC_jvec(double *y, int n, double m, double s, int logscale) { 
//Joint density of y[0]...y[n-1] under Normal(m,s^2), i.e. returns scalar
  int _i; double ans;
  for (_i=0, ans=0; _i<n; _i++) {
    ans+= dnormC(y[_i],m,s,1);
  }
  if (logscale!=1) ans= exp(ans);
  return(ans);
}

double dmvnormC(double *y, int n, double *mu, double **cholsinv, double det, int logscale) { 
/* Density of multivariate Normal evaluated at y[1]...y[n]. mu is the mean. chols and
   are the Cholesky decomposition and the determinant of the inverse covariance matrix.
   Example of usage: 
     choldc_inv(s,n,cholsinv); 
     det= choldc_det(cholsinv,n);
     dmvnormC(y,n,mu,cholsinv,det,0); */
  int i;  
  double *z,*z2, res;

  //Find (y-mu)' * cholsinv' * cholsinv * (y-mu)
  z= dvector(1,n); z2= dvector(1,n);
  for (i=1; i<=n; i++) { z[i]= y[i]-mu[i]; }
  Ax(cholsinv,z,z2,1,n,1,n);
  for (res=0, i=1; i<=n; i++) { res += z2[i]*z2[i]; }
  free_dvector(z,1,n); free_dvector(z2,1,n);

  if(logscale==1) return(-n*log(SQ_M_PI_2) + .5*log(det) - .5*res);
  else return(exp(-n*log(SQ_M_PI_2) + .5*log(det) - .5*res));

}


double	qnormC (double cdf, double m, double s) 
/* returns inv cdf of normal N(m,s^2) at p */ 
{	 
  double /* used to all be floats.. */
    p,q,mean,sd,bound,x; 
  double  
    y; 
  int  
    status,which; 
 
  if( (cdf < 0.0) | (cdf > 1.0) ){ 
    errorC("qnormC",  "Tried inverse cdf with p<0 or p>1", 1); 
  }  
  /* par check */ 
  if (cdf <= 2.86e-07) 
    y = -5.0*s+m; 
  else if (cdf >= 0.9999997) 
    y = 5.0*s+m; 
  else{ 
    /* primitive type conversion */ 
    p = cdf; 
    q = 1.0-p;
    mean = m; 
    sd = s; 
    which = 2; 
     
    cdfnor(&which,&p,&q, &x,&mean,&sd,&status,&bound); 
   
    y = x; /* another primitive type conversion */ 
  } 
  return y; 
} 


/* returns draw from binomial(n,p) */  
int rbinomial(int n, double p) 
{ 	 
  int i,x; 
  for(i=0,x=0;i<n;i++) 
    x += (runif() < p) ? 1 : 0; 
  return x; 
} 

double dbinomial(int x, int n, double p, int logscale) {
  double ans;
  ans= lnchoose(n,x) + (x+.0)*log(p) + (n-x+.0)*log(1-p);
  if (logscale==1) return(ans); else return(exp(ans));
}


/* returns draw from multinomial with cell prob pr  
  ------------------------------------------   
  Value: 
    x:   vector of indices indicating draws 
         x in [0,..n_cell-1] 
  ------------------------------------------   
  Args: 
    n_draw: number of draws 
    n_cell: number of cells 
    pr:     n_cell vector of cell probs 
            (not necessarily standardized) 
    x:      n_draw vector of indices 
            (OUTPUT) 
*/ 
void rmultinomial(int n_draw, int n_cell, double *pr, int *x) 
{ 
  double *cum_p, uj; 
  int i,j; 
 
  cum_p = dvector(0,n_cell); 
 
  for(i=1,cum_p[0]=pr[0];i<n_cell;i++)     
    cum_p[i] = cum_p[i-1]+pr[i]; 
  for(j=0;j<n_draw;j++){ 
    uj = runif()*cum_p[n_cell-1]; 
    for(i=0; ((uj > cum_p[i]) & (i<n_cell)); i++); 
    x[j] = i; 
  } 
 
  free_dvector(cum_p,0,n_cell); 
} 


// Draw from univariate Normal(mu,s^2)
double rnormC(double mu, double s) {

static int iset=0;
static double gset;
double fac,rsq,v1,v2;

  if (iset == 0) { //We don't have an extra deviate handy, so
    do {
      v1=2.0*runif()-1.0; //pick two uniform numbers in the square extending from
      v2=2.0*runif()-1.0; //-1 to +1 in each direction,
      rsq=v1*v1+v2*v2;       //see if they are in the unit circle,
    } while (rsq >= 1.0 || rsq == 0.0); //and if they are not, try again.
    fac=sqrt(-2.0*log(rsq)/rsq);
    //Now make the Box-Muller transformation to get two normal deviates. Return one and
    //save the other for next time.
    gset=v1*fac;
    iset=1; //Set flag.
    return v2*fac*s + mu;
  } else { //We have an extra deviate handy,
    iset=0; //so unset the flag, and return it.
    return gset*s + mu; 
  }
}


/* Draw from a univariate truncated Normal giving truncation points */
double rnorm_trunc(double ltrunc, double rtrunc, double m, double s) {
  // ltrunc: left truncation point; rtrunc: right truncation point, m: mean; s: SD
  double lprob, rprob;
  lprob= pnormC(ltrunc,m,s); rprob= pnormC(rtrunc,m,s);
  return(rnorm_trunc_prob(lprob,rprob,m,s));
}

/* Draw from a univariate truncated Normal giving truncation probabilities */
double rnorm_trunc_prob(double lprob, double rprob, double m, double s) {
  // lprob, rprob: prob to the left of the truncation points; m: mean; s: SD
  //e.g. lprob=.05, rprob=.99 means we're truncating the lower 5% and the upper 1%
  double u;
  if (lprob>=rprob) nrerror("rnorm_trunc_prob","left truncation probability is larger than right truncation probability","");
  u= lprob + runif()*(rprob-lprob);  //generate uniform between lprob, rprob
  u= qnormC(u,m,s);
  return(u);

}


// Draw from multivar Normal with n dimensions
void rmvnormC(double *y, int n, double *mu, double **chols) {
/* Result is stored in y[1..n]. mu is the location parameter, chols is the Cholesky decomposition
   of the covariance matrix. That is, the covariance is s and s=chols*chols'
   Note: both y and mu should have length n, and s should be an n*n matrix. The routine doesn't
   check it 
   Example: choldc(s,n,chols); //compute cholesky decomposition
            rmvnormC(y,n,mu,chols); //generate random variate */

  int i;
  double *z;

  z= dvector(1,n);
  for (i=1;i<=n;i++) { z[i]= rnormC(0,1); } //generate n independent draws from a univariate Normal
  Ax_plus_y(chols,z,mu,y,1,n);             //compute mu + chols*z
  free_dvector(z,1,n);

}


double dtC(double y, double mu, double s, int nu) { 
//Density of t with nu df, location mu and scale s^2
  double normk, t1, t2;

  t2= .5*nu; t1= t2 + .5;

  normk= exp(gamln(&t1)-gamln(&t2))/(sqrt(nu*M_PI)*s);
  return(normk*pow(1+(y-mu)*(y-mu)/(s*s*(nu+.0)),-t1));

}

double dmvtC(double *y, int n, double *mu, double **cholsinv, double det, int nu, int logscale) {
/* Density of multivariate T with nu df and n dimensions. mu: location, cholsinv, det: cholesky
   decomp and determinant of inverse cov matrix, logscale: set to 1 to return density in log-scale
   Example of usage: 
     choldc_inv(s,n,cholsinv); 
     det= choldc_det(cholsinv,n);
     dmvtC(y,n,mu,cholsinv,det,nu,0); */
  int i;  
  double *z,*z2, res, t1, t2, normk;

  //Find (y-mu)' * cholsinv' * cholsinv * (y-mu)
  z= dvector(1,n); z2= dvector(1,n);
  for (i=1; i<=n; i++) { z[i]= y[i]-mu[i]; }
  Ax(cholsinv,z,z2,1,n,1,n);
  for (res=0, i=1; i<=n; i++) { res += z2[i]*z2[i]; }
  free_dvector(z,1,n); free_dvector(z2,1,n);

  t2= .5*nu; t1= t2+.5*(n+.0);
  normk= gamln(&t1)-gamln(&t2)-.5*(n+.0)*(log(nu+.0)+log(M_PI))+.5*log(det);

  if (logscale==1) return(normk -t1*log(1+res/(nu+.0)));
  else return(exp(normk) * pow(1+res/(nu+.0),-t1));

}

/* Draw from a univariate standard t with nu degrees of freedom */
double rtC(int nu) {
  double x, z;

  z= rnormC(0,1);
  x= gengam(.5,nu/2.0);  //draw from chi-square with nu degrees of freedom
  return(z*sqrt(nu/x));

}

/* Draw from a univariate truncated t with nu degrees of freedom giving truncation points */
double rt_trunc(int nu, double ltrunc, double rtrunc) {
  // nu: degrees of freedom; ltrunc: left truncation point; rtrunc: right truncation point
  double lprob, rprob;
  lprob= ptC(ltrunc,nu); rprob= ptC(rtrunc,nu);
  return(rt_trunc_prob(nu,lprob,rprob));
}

/* Draw from a univariate truncated t with nu degrees of freedom giving truncation probabilities */
double rt_trunc_prob(int nu, double lprob, double rprob) {
  //nu: degrees of freedom; lprob, rprob: prob to the left of the truncation points
  //e.g. lprob=.05, rprob=.99 means we're truncating the lower 5% and the upper 1%
  double u;
  if (lprob>=rprob) nrerror("rt_trunc_prob","left truncation probability is larger than right truncation probability","");
  u= lprob + runif()*(rprob-lprob);  //generate uniform between lprob, rprob
  return(qtC(u,nu));

}

/* Find quantiles of a t-distribution with nu degrees of freedom */
double qtC(double p, int nu) {
  /* p: probability; nu: degrees of freedom; lower_tail==1 means p indicates left tail area */
 /* * @author Sundar Dorai-Raj
  * * See the GNU General Public License for more details at http://www.gnu.org * * */
  // Algorithm 396: Student's t-quantiles by G.W. Hill CACM 13(10), 619-620, October 1970
  int neg;
  double ndf, eps, P, q, prob, a, b, c, d, y, x;

  ndf= nu + 0.0;
  if(p<=0 || p>=1 || ndf<1) return(-1);
  eps=1e-12;
  if(p > 0.5) {
    neg = 0;
    P = 2.0 * (1.0 - p);
  }
  else {
    neg = 1;
    P = 2.0 * p;
  }

  if(abs(ndf - 2.0) < eps) {   /* df ~= 2 */
    q=sqrt(2.0 / (P * (2.0 - P)) - 2.0);
  }
  else if (ndf < 1 + eps) {   /* df ~= 1 */
    prob = P * M_PI_2;
    q = cos(prob)/sin(prob);
  }
  else {      /*-- usual case;  including, e.g.,  df = 1.1 */
    a = 1.0 / (ndf - 0.5);
    b = 48.0 / (a * a);
    c = ((20700.0 * a / b - 98.0) * a - 16.0) * a + 96.36;
    d = ((94.5 / (b + c) - 3.0) / b + 1.0) * sqrt(a * M_PI_2) * ndf;
    y = pow(d * P, 2.0 / ndf);
    if (y > 0.05 + a) {  /* Asymptotic inverse expansion about normal */
      x = qnormC(0.5*P,0.0,1.0);
      y = x * x;
      if (ndf < 5)
	c += 0.3 * (ndf - 4.5) * (x + 0.6);
      c = (((0.05 * d * x - 5.0) * x - 7.0) * x - 2.0) * x + b + c;
      y = (((((0.4 * y + 6.3) * y + 36.0) * y + 94.5) / c - y - 3.0) / b + 1.0) * x;
      y = a * y * y;
      if (y > 0.002) y = exp(y) - 1.0;
      else { /* Taylor of    e^y -1 : */
	y = (0.5 * y + 1.0) * y;
      }
    } else {
      y = ((1.0 / (((ndf + 6.0) / (ndf * y) - 0.089 * d - 0.822) * (ndf + 2.0) * 3.0) + 0.5 / (ndf + 4.0))* y - 1.0) * (ndf + 1.0) / (ndf + 2.0) + 1.0 / y;
    }
    q = sqrt(ndf * y);
  }
  if (neg==1) q = -q;
  return(q);
}

/* CDF of a t-Student distribution */
double ptC(double x, int nu) {

  if (x>0) { return(1-0.5*pbetaC((nu+0.0)/(x*x+nu),0.5*nu,0.5)); }
  else if (x<0) { return(0.5*pbetaC((nu+0.0)/(x*x+nu),0.5*nu,0.5)); }
  else { return(0.5); }

}


// Draw from multivar T with n dimensions and nu degrees of freedom
void rmvtC(double *y, int n, double *mu, double **chols, int nu) {
/* Result is stored in y[1..n]. mu is the location parameter, chols is the Cholesky decomposition
   of the covariance matrix. That is, the covariance is s*nu/(nu-2) and s=chols*chols'
   and nu are the degrees of freedom 
   Note: both y and mu should have length n, and s should be an n*n matrix. The routine doesn't
   check it 
   Example: choldc(s,n,chols); //compute cholesky decomposition
            rmvtC(y,n,mu,chols,nu); //generate random variate */

  int i;
  double x, *z;

  x= sqrt(nu/gengam(.5,nu/2.0));  //draw from chi-square with nu degrees of freedom
  z= dvector(1,n);
  for (i=1;i<=n;i++) { z[i]= x*rnormC(0,1); } //multiple n indep normal draws by the common chi-square
  Ax_plus_y(chols,z,mu,y,1,n);          //compute mu + chols*z
  free_dvector(z,1,n);

}


double rgammaC(double a, double b) {
  //Generate from a Gamma(a,b) (a is shape; b location; mean= a/b)
  return(gengam(b,a));
}

double dgammaC(double x, double a, double b) {
  //Density of a Gamma(a,b) (a is shape; b location; mean= a/b)
  if (x!=0) { 
    return(exp(a*log(b)-gamln(&a)+(a-1)*log(x)-x*b)); 
  } else {
    if (a==1) return(b); else return(0);
  }
}

double dinvgammaC(double x, double a, double b) {
 //Density of an Inverse Gamma(a,b) (a: shape; b: location; mean of 1/x= a/b)
  if (x!=0) { 
    return(exp(a*log(b)-gamln(&a)-(a+1)*log(x)-b/x)); 
  } else {
    return(0);
  }
}


/************************************************************************
                       MORE RANDOM VARIATE STUFF
************************************************************************/

/************************************************************************
FIFDINT:
Truncates a double precision number to an integer and returns the
value in a double.
************************************************************************/
double fifdint(double a)
/* a     -     number to be truncated */
{
  long temp;
  temp = (long)(a);
  return (double)(temp);
}


void cdfnor(int *which,double *p,double *q,double *x,double *mean,
	    double *sd,int *status,double *bound)
/**********************************************************************

      void cdfnor(int *which,double *p,double *q,double *x,double *mean,
            double *sd,int *status,double *bound)

               Cumulative Distribution Function
               NORmal distribution


                              Function


     Calculates any one parameter of the normal
     distribution given values for the others.


                              Arguments


     WHICH  --> Integer indicating  which of the  next  parameter
     values is to be calculated using values  of the others.
     Legal range: 1..4
               iwhich = 1 : Calculate P and Q from X,MEAN and SD
               iwhich = 2 : Calculate X from P,Q,MEAN and SD
               iwhich = 3 : Calculate MEAN from P,Q,X and SD
               iwhich = 4 : Calculate SD from P,Q,X and MEAN

     P <--> The integral from -infinity to X of the normal density.
            Input range: (0,1].

     Q <--> 1-P.
            Input range: (0, 1].
            P + Q = 1.0.

     X < --> Upper limit of integration of the normal-density.
             Input range: ( -infinity, +infinity)

     MEAN <--> The mean of the normal density.
               Input range: (-infinity, +infinity)

     SD <--> Standard Deviation of the normal density.
             Input range: (0, +infinity).

     STATUS <-- 0 if calculation completed correctly
               -I if input parameter number I is out of range
                1 if answer appears to be lower than lowest
                  search bound
                2 if answer appears to be higher than greatest
                  search bound
                3 if P + Q .ne. 1

     BOUND <-- Undefined if STATUS is 0

               Bound exceeded by parameter number I if STATUS
               is negative.

               Lower search bound if STATUS is 1.

               Upper search bound if STATUS is 2.


                              Method




     A slightly modified version of ANORM from

     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
     Package of Special Function Routines and Test Drivers"
     acm Transactions on Mathematical Software. 19, 22-32.

     is used to calulate the  cumulative standard normal distribution.

     The rational functions from pages  90-95  of Kennedy and Gentle,
     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
     starting values to Newton's Iterations which compute the inverse
     standard normal.  Therefore no  searches  are necessary for  any
     parameter.

     For X < -15, the asymptotic expansion for the normal is used  as
     the starting value in finding the inverse standard normal.
     This is formula 26.2.12 of Abramowitz and Stegun.


                              Note


      The normal density is proportional to
      exp( - 0.5 * (( X - MEAN)/SD)**2)

**********************************************************************/
{
static int K1 = 1;
static double z,pq;
/*
     ..
     .. Executable Statements ..
*/
/*
     Check arguments
*/
    *status = 0;
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
/*
     P
*/
    if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p <= 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
/*
     Q
*/
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 1) goto S150;
/*
     P + Q
*/
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*spmpar(&K1))) goto S140;
    if(!(pq < 0.0e0)) goto S120;
    *bound = 0.0e0;
    goto S130;
S120:
    *bound = 1.0e0;
S130:
    *status = 3;
    return;
S150:
S140:
    if(*which == 4) goto S170;
/*
     SD
*/
    if(!(*sd <= 0.0e0)) goto S160;
    *bound = 0.0e0;
    *status = -6;
    return;
S170:
S160:
/*
     Calculate ANSWERS
*/
    if(1 == *which) {
/*
     Computing P
*/
        z = (*x-*mean)/ *sd;
        cumnor(&z,p,q);
    }
    else if(2 == *which) {
/*
     Computing X
*/
        z = dinvnr(p,q);
        *x = *sd*z+*mean;
    }
    else if(3 == *which) {
/*
     Computing the MEAN
*/
        z = dinvnr(p,q);
        *mean = *x-*sd*z;
    }
    else if(4 == *which) {
/*
     Computing SD
*/
        z = dinvnr(p,q);
        *sd = (*x-*mean)/z;
    }
    return;
}

double spmpar(int *i)
/*
-----------------------------------------------------------------------
 
     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
 
        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
 
        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
 
        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
 
-----------------------------------------------------------------------
     WRITTEN BY
        ALFRED H. MORRIS, JR.
        NAVAL SURFACE WARFARE CENTER
        DAHLGREN VIRGINIA
-----------------------------------------------------------------------
-----------------------------------------------------------------------
     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
-----------------------------------------------------------------------
*/
{
static int K1 = 4;
static int K2 = 8;
static int K3 = 9;
static int K4 = 10;
static double spmpar,b,binv,bm1,one,w,z;
static int emax,emin,ibeta,m;
/*
     ..
     .. Executable Statements ..
*/
    if(*i > 1) goto S10;
    b = ipmpar(&K1);
    m = ipmpar(&K2);
    spmpar = pow(b,(double)(1-m));
    return spmpar;
S10:
    if(*i > 2) goto S20;
    b = ipmpar(&K1);
    emin = ipmpar(&K3);
    one = 1.0;
    binv = one/b;
    w = pow(b,(double)(emin+2));
    spmpar = w*binv*binv*binv;
    return spmpar;
S20:
    ibeta = ipmpar(&K1);
    m = ipmpar(&K2);
    emax = ipmpar(&K4);
    b = ibeta;
    bm1 = ibeta-1;
    one = 1.0;
    z = pow(b,(double)(m-1));
    w = ((z-one)*b+bm1)/(b*z);
    z = pow(b,(double)(emax-2));
    spmpar = w*z*b*b;
    return spmpar;
}


void cumnor(double *arg,double *result,double *ccum)
/*
**********************************************************************
 
     void cumnor(double *arg,double *result,double *ccum)
 
 
                              Function
 
 
     Computes the cumulative  of    the  normal   distribution,   i.e.,
     the integral from -infinity to x of
          (1/sqrt(2*pi)) exp(-u*u/2) du
 
     X --> Upper limit of integration.
                                        X is DOUBLE PRECISION
 
     RESULT <-- Cumulative normal distribution.
                                        RESULT is DOUBLE PRECISION
 
     CCUM <-- Compliment of Cumulative normal distribution.
                                        CCUM is DOUBLE PRECISION
 
     Renaming of function ANORM from:

     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
     Package of Special Function Routines and Test Drivers"
     acm Transactions on Mathematical Software. 19, 22-32.

     with slight modifications to return ccum and to deal with
     machine constants.
 
**********************************************************************
  Original Comments:
------------------------------------------------------------------
 
 This function evaluates the normal distribution function:
 
                              / x
                     1       |       -t*t/2
          P(x) = ----------- |      e       dt
                 sqrt(2 pi)  |
                             /-oo
 
   The main computation evaluates near-minimax approximations
   derived from those in "Rational Chebyshev approximations for
   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
   This transportable program uses rational functions that
   theoretically approximate the normal distribution function to
   at least 18 significant decimal digits.  The accuracy achieved
   depends on the arithmetic system, the compiler, the intrinsic
   functions, and proper selection of the machine-dependent
   constants.
 
*******************************************************************
*******************************************************************
 
 Explanation of machine-dependent constants.
 
   MIN   = smallest machine representable number.
 
   EPS   = argument below which anorm(x) may be represented by
           0.5  and above which  x*x  will not underflow.
           A conservative value is the largest machine number X
           such that   1.0 + X = 1.0   to machine precision.
*******************************************************************
*******************************************************************
 
 Error returns
 
  The program returns  ANORM = 0     for  ARG .LE. XLOW.
 
 
 Intrinsic functions required are:
 
     ABS, AINT, EXP
 
 
  Author: W. J. Cody
          Mathematics and Computer Science Division
          Argonne National Laboratory
          Argonne, IL 60439
 
  Latest modification: March 15, 1992
 
------------------------------------------------------------------
*/
{
static double a[5] = {
    2.2352520354606839287e00,1.6102823106855587881e02,1.0676894854603709582e03,
    1.8154981253343561249e04,6.5682337918207449113e-2
};
static double b[4] = {
    4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
    4.5507789335026729956e04
};
static double c[9] = {
    3.9894151208813466764e-1,8.8831497943883759412e00,9.3506656132177855979e01,
    5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
    1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
};
static double d[8] = {
    2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
    6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
    3.8912003286093271411e04,1.9685429676859990727e04
};
static double half = 0.5e0;
static double p[6] = {
    2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
    1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
};
static double one = 1.0e0;
static double q[5] = {
    1.28426009614491121e00,4.68238212480865118e-1,6.59881378689285515e-2,
    3.78239633202758244e-3,7.29751555083966205e-5
};
static double sixten = 1.60e0;
static double sqrpi = 3.9894228040143267794e-1;
static double thrsh = 0.66291e0;
static double root32 = 5.656854248e0;
static double zero = 0.0e0;
static int K1 = 1;
static int K2 = 2;
static int i;
static double del,eps,temp,x,xden,xnum,y,xsq,min;
/*
------------------------------------------------------------------
  Machine dependent constants
------------------------------------------------------------------
*/
    eps = spmpar(&K1)*0.5e0;
    min = spmpar(&K2);
    x = *arg;
    y = fabs(x);
    if(y <= thrsh) {
/*
------------------------------------------------------------------
  Evaluate  anorm  for  |X| <= 0.66291
------------------------------------------------------------------
*/
        xsq = zero;
        if(y > eps) xsq = x*x;
        xnum = a[4]*xsq;
        xden = xsq;
        for(i=0; i<3; i++) {
            xnum = (xnum+a[i])*xsq;
            xden = (xden+b[i])*xsq;
        }
        *result = x*(xnum+a[3])/(xden+b[3]);
        temp = *result;
        *result = half+temp;
        *ccum = half-temp;
    }
/*
------------------------------------------------------------------
  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
------------------------------------------------------------------
*/
    else if(y <= root32) {
        xnum = c[8]*y;
        xden = y;
        for(i=0; i<7; i++) {
            xnum = (xnum+c[i])*y;
            xden = (xden+d[i])*y;
        }
        *result = (xnum+c[7])/(xden+d[7]);
        xsq = fifdint(y*sixten)/sixten;
        del = (y-xsq)*(y+xsq);
        *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
        *ccum = one-*result;
        if(x > zero) {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
/*
------------------------------------------------------------------
  Evaluate  anorm  for |X| > sqrt(32)
------------------------------------------------------------------
*/
    else  {
        *result = zero;
        xsq = one/(x*x);
        xnum = p[5]*xsq;
        xden = xsq;
        for(i=0; i<4; i++) {
            xnum = (xnum+p[i])*xsq;
            xden = (xden+q[i])*xsq;
        }
        *result = xsq*(xnum+p[4])/(xden+q[4]);
        *result = (sqrpi-*result)/y;
        xsq = fifdint(x*sixten)/sixten;
        del = (x-xsq)*(x+xsq);
        *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
        *ccum = one-*result;
        if(x > zero) {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
    if(*result < min) *result = 0.0e0;
/*
------------------------------------------------------------------
  Fix up for negative argument, erf, etc.
------------------------------------------------------------------
----------Last card of ANORM ----------
*/
    if(*ccum < min) *ccum = 0.0e0;
}

double dinvnr(double *p,double *q)
/*
**********************************************************************
 
     double dinvnr(double *p,double *q)
     Double precision NoRmal distribution INVerse
 
 
                              Function
 
 
     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
 
 
                              Arguments
 
 
     P --> The probability whose normal deviate is sought.
                    P is DOUBLE PRECISION
 
     Q --> 1-P
                    P is DOUBLE PRECISION
 
 
                              Method
 
 
     The  rational   function   on  page 95    of Kennedy  and  Gentle,
     Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
     value for the Newton method of finding roots.
 
 
                              Note
 
 
     If P or Q .lt. machine EPS returns +/- DINVNR(EPS)
 
**********************************************************************
*/
{
#define maxit 100
#define eps 1.0e-13
#define r2pi 0.3989422804014326e0
#define nhalf -0.5e0
#define dennor(x) (r2pi*exp(nhalf*(x)*(x)))
static double dinvnr,strtx,xcur,cum,ccum,pp,dx;
static int i;
static unsigned long qporq;
/*
     ..
     .. Executable Statements ..
*/
/*
     FIND MINIMUM OF P AND Q
*/
    qporq = *p <= *q;
    if(!qporq) goto S10;
    pp = *p;
    goto S20;
S10:
    pp = *q;
S20:
/*
     INITIALIZATION STEP
*/
    strtx = stvaln(&pp);
    xcur = strtx;
/*
     NEWTON INTERATIONS
*/
    for(i=1; i<=maxit; i++) {
        cumnor(&xcur,&cum,&ccum);
        dx = (cum-pp)/dennor(xcur);
        xcur -= dx;
        if(fabs(dx/xcur) < eps) goto S40;
    }
    dinvnr = strtx;
/*
     IF WE GET HERE, NEWTON HAS FAILED
*/
    if(!qporq) dinvnr = -dinvnr;
    return dinvnr;
S40:
/*
     IF WE GET HERE, NEWTON HAS SUCCEDED
*/
    dinvnr = xcur;
    if(!qporq) dinvnr = -dinvnr;
    return dinvnr;
#undef maxit
#undef eps
#undef r2pi
#undef nhalf
#undef dennor
}


double stvaln(double *p)
/*
**********************************************************************
 
     double stvaln(double *p)
                    STarting VALue for Neton-Raphon
                calculation of Normal distribution Inverse
 
 
                              Function
 
 
     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
 
 
                              Arguments
 
 
     P --> The probability whose normal deviate is sought.
                    P is DOUBLE PRECISION
 
 
                              Method
 
 
     The  rational   function   on  page 95    of Kennedy  and  Gentle,
     Statistical Computing, Marcel Dekker, NY , 1980.
 
**********************************************************************
*/
{
static double xden[5] = {
    0.993484626060e-1,0.588581570495e0,0.531103462366e0,0.103537752850e0,
    0.38560700634e-2
};
static double xnum[5] = {
    -0.322232431088e0,-1.000000000000e0,-0.342242088547e0,-0.204231210245e-1,
    -0.453642210148e-4
};
static int K1 = 5;
static double stvaln,sign,y,z;
/*
     ..
     .. Executable Statements ..
*/
    if(!(*p <= 0.5e0)) goto S10;
    sign = -1.0e0;
    z = *p;
    goto S20;
S10:
    sign = 1.0e0;
    z = 1.0e0-*p;
S20:
    y = sqrt(-(2.0e0*log(z)));
    stvaln = y+devlpl(xnum,&K1,&y)/devlpl(xden,&K1,&y);
    stvaln = sign*stvaln;
    return stvaln;
}

double devlpl(double a[],int *n,double *x)
/*
**********************************************************************
 
     double devlpl(double a[],int *n,double *x)
              Double precision EVALuate a PoLynomial at X
 
 
                              Function
 
 
     returns
          A(1) + A(2)*X + ... + A(N)*X**(N-1)
 
 
                              Arguments
 
 
     A --> Array of coefficients of the polynomial.
                                        A is DOUBLE PRECISION(N)
 
     N --> Length of A, also degree of polynomial - 1.
                                        N is INTEGER
 
     X --> Point at which the polynomial is to be evaluated.
                                        X is DOUBLE PRECISION
 
**********************************************************************
*/
{
static double devlpl,term;
static int i;
/*
     ..
     .. Executable Statements ..
*/
    term = a[*n-1];
    for(i= *n-1-1; i>=0; i--) term = a[i]+term**x;
    devlpl = term;
    return devlpl;
}


int ipmpar(int*);
/*
-----------------------------------------------------------------------
 
     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...
 
  INTEGERS.
 
     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM
 
               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
 
               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
 
     IPMPAR(1) = A, THE BASE.
 
     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.
 
     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
 
  FLOATING-POINT NUMBERS.
 
     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
     NONZERO NUMBERS ARE REPRESENTED IN THE FORM
 
               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
 
               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
 
     IPMPAR(4) = B, THE BASE.
 
  SINGLE-PRECISION
 
     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.
 
     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
 
     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
 
  DOUBLE-PRECISION
 
     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.
 
     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
 
     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
 
-----------------------------------------------------------------------
 
     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED REMOVE
     THE COMMENT DELIMITORS FROM THE DEFINITIONS DIRECTLY BELOW THE NAME
     OF THE MACHINE
 
-----------------------------------------------------------------------
 
     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.
 
-----------------------------------------------------------------------
     .. Scalar Arguments ..
*/
int ipmpar(int *i)
{
static int imach[11];
static int ipmpar;
/*     MACHINE CONSTANTS FOR AMDAHL MACHINES. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 16;
   imach[5] = 6;
   imach[6] = -64;
   imach[7] = 63;
   imach[8] = 14;
   imach[9] = -64;
   imach[10] = 63;
*/
/*     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
       PC 7300, AND AT&T 6300. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM. */
/*
   imach[1] = 2;
   imach[2] = 33;
   imach[3] = 8589934591;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -256;
   imach[7] = 255;
   imach[8] = 60;
   imach[9] = -256;
   imach[10] = 255;
*/
/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM. */
/*
   imach[1] = 2;
   imach[2] = 39;
   imach[3] = 549755813887;
   imach[4] = 8;
   imach[5] = 13;
   imach[6] = -50;
   imach[7] = 76;
   imach[8] = 26;
   imach[9] = -50;
   imach[10] = 76;
*/
/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS. */
/*
   imach[1] = 2;
   imach[2] = 39;
   imach[3] = 549755813887;
   imach[4] = 8;
   imach[5] = 13;
   imach[6] = -50;
   imach[7] = 76;
   imach[8] = 26;
   imach[9] = -32754;
   imach[10] = 32780;
*/
/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
       60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
       ARITHMETIC (NOS OPERATING SYSTEM). */
/*
   imach[1] = 2;
   imach[2] = 48;
   imach[3] = 281474976710655;
   imach[4] = 2;
   imach[5] = 48;
   imach[6] = -974;
   imach[7] = 1070;
   imach[8] = 95;
   imach[9] = -926;
   imach[10] = 1070;
*/
/*     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
       ARITHMETIC (NOS/VE OPERATING SYSTEM). */
/*
   imach[1] = 2;
   imach[2] = 63;
   imach[3] = 9223372036854775807;
   imach[4] = 2;
   imach[5] = 48;
   imach[6] = -4096;
   imach[7] = 4095;
   imach[8] = 96;
   imach[9] = -4096;
   imach[10] = 4095;
*/
/*     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3. */
/*
   imach[1] = 2;
   imach[2] = 63;
   imach[3] = 9223372036854775807;
   imach[4] = 2;
   imach[5] = 47;
   imach[6] = -8189;
   imach[7] = 8190;
   imach[8] = 94;
   imach[9] = -8099;
   imach[10] = 8190;
*/
/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200. */
/*
   imach[1] = 2;
   imach[2] = 15;
   imach[3] = 32767;
   imach[4] = 16;
   imach[5] = 6;
   imach[6] = -64;
   imach[7] = 63;
   imach[8] = 14;
   imach[9] = -64;
   imach[10] = 63;
*/
/*     MACHINE CONSTANTS FOR THE HARRIS 220. */
/*
   imach[1] = 2;
   imach[2] = 23;
   imach[3] = 8388607;
   imach[4] = 2;
   imach[5] = 23;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 38;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000
       AND DPS 8/70 SERIES. */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 63;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HP 2100
       3 WORD DOUBLE PRECISION OPTION WITH FTN4 */
/*
   imach[1] = 2;
   imach[2] = 15;
   imach[3] = 32767;
   imach[4] = 2;
   imach[5] = 23;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 39;
   imach[9] = -128;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HP 2100
       4 WORD DOUBLE PRECISION OPTION WITH FTN4 */
/*
   imach[1] = 2;
   imach[2] = 15;
   imach[3] = 32767;
   imach[4] = 2;
   imach[5] = 23;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 55;
   imach[9] = -128;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HP 9000. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -126;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
       THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
       5/7/9 AND THE SEL SYSTEMS 85/86. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 16;
   imach[5] = 6;
   imach[6] = -64;
   imach[7] = 63;
   imach[8] = 14;
   imach[9] = -64;
   imach[10] = 63;
*/
/*     MACHINE CONSTANTS FOR THE IBM PC. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
       MACFORTRAN II. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 56;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR). */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 54;
   imach[9] = -101;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR). */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 62;
   imach[9] = -128;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
       32-BIT INTEGER ARITHMETIC. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 56;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
       SERIES (MIPS R3000 PROCESSOR). */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
       3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
       PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300). */

   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 60;
   imach[9] = -1024;
   imach[10] = 1023;
*/
/*     MACHINE CONSTANTS FOR THE VAX 11/780. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 56;
   imach[9] = -127;
   imach[10] = 127;
*/
    ipmpar = imach[*i];
    return ipmpar;
}


/************************************************************************
                        EVEN MORE STUFF
************************************************************************/

double genunf(double low,double high)
/*
**********************************************************************
     double genunf(double low,double high)
               GeNerate Uniform Real between LOW and HIGH
                              Function
     Generates a real uniformly distributed between LOW and HIGH.
                              Arguments
     low --> Low bound (exclusive) on real value to be generated
     high --> High bound (exclusive) on real value to be generated
**********************************************************************
*/
{
static double genunf;

    if(!(low > high)) goto S10;
    printf("LOW > HIGH in GENUNF: LOW %16.6E HIGH: %16.6E\n",low,high);
    puts("Abort");
    exit(1);
S10:
    genunf = low+(high-low)*ranf();
    return genunf;
}


double gengam(double a,double r)
/*
**********************************************************************
     double gengam(double a,double r)
           GENerates random deviates from GAMma distribution
                              Function
     Generates random deviates from the gamma distribution whose
     density is
          (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
                              Arguments
     a --> Location parameter of Gamma distribution
     r --> Shape parameter of Gamma distribution
     CAREFUL: order of parameters is reversed wrt usual notation. mean is r/a; var is r/a^2
                              Method
     Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.
     For details see:
               (Case R >= 1.0)
               Ahrens, J.H. and Dieter, U.
               Generating Gamma Variates by a
               Modified Rejection Technique.
               Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
     Algorithm GD
               (Case 0.0 <= R <= 1.0)
               Ahrens, J.H. and Dieter, U.
               Computer Methods for Sampling from Gamma,
               Beta, Poisson and Binomial Distributions.
               Computing, 12 (1974), 223-246/
     Adapted algorithm GS.
**********************************************************************
*/
{
static double gengam;

    gengam = sgamma(r);
    gengam /= a;
    return gengam;
}


double sgamma(double a)
/*
**********************************************************************
                                                                      
                                                                      
     (STANDARD-)  G A M M A  DISTRIBUTION                             
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
               PARAMETER  A >= 1.0  !                                 
                                                                      
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               GENERATING GAMMA VARIATES BY A                         
               MODIFIED REJECTION TECHNIQUE.                          
               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  
                                                                      
     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     
                                 (STRAIGHTFORWARD IMPLEMENTATION)     
                                                                      
     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
     SUNIF.  The argument IR thus goes away.                          
                                                                      
**********************************************************************
                                                                      
               PARAMETER  0.0 < A < 1.0  !                            
                                                                      
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               COMPUTER METHODS FOR SAMPLING FROM GAMMA,              
               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              
               COMPUTING, 12 (1974), 223 - 246.                       
                                                                      
     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    
                                                                      
**********************************************************************
     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
*/
{
extern double fsign( double num, double sign );
static double q1 = 4.166669E-2;
static double q2 = 2.083148E-2;
static double q3 = 8.01191E-3;
static double q4 = 1.44121E-3;
static double q5 = -7.388E-5;
static double q6 = 2.4511E-4;
static double q7 = 2.424E-4;
static double a1 = 0.3333333;
static double a2 = -0.250003;
static double a3 = 0.2000062;
static double a4 = -0.1662921;
static double a5 = 0.1423657;
static double a6 = -0.1367177;
static double a7 = 0.1233795;
static double e1 = 1.0;
static double e2 = 0.4999897;
static double e3 = 0.166829;
static double e4 = 4.07753E-2;
static double e5 = 1.0293E-2;
static double aa = 0.0;
static double aaa = 0.0;
static double sqrt32 = 5.656854;
static double sgamma,s2,s,d,t,x,u,r,q0,b,si,c,v,q,e,w,p;
    if(a == aa) goto S10;
    if(a < 1.0) goto S120;
/*
     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
*/
    aa = a;
    s2 = a-0.5;
    s = sqrt(s2);
    d = sqrt32-12.0*s;
S10:
/*
     STEP  2:  T=STANDARD NORMAL DEVIATE,
               X=(S,1/2)-NORMAL DEVIATE.
               IMMEDIATE ACCEPTANCE (I)
*/
    t = snorm();
    x = s+0.5*t;
    sgamma = x*x;
    if(t >= 0.0) return sgamma;
/*
     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
*/
    u = ranf();
    if(d*u <= t*t*t) return sgamma;
/*
     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
*/
    if(a == aaa) goto S40;
    aaa = a;
    r = 1.0/ a;
    q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
/*
               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
*/
    if(a <= 3.686) goto S30;
    if(a <= 13.022) goto S20;
/*
               CASE 3:  A .GT. 13.022
*/
    b = 1.77;
    si = 0.75;
    c = 0.1515/s;
    goto S40;
S20:
/*
               CASE 2:  3.686 .LT. A .LE. 13.022
*/
    b = 1.654+7.6E-3*s2;
    si = 1.68/s+0.275;
    c = 6.2E-2/s+2.4E-2;
    goto S40;
S30:
/*
               CASE 1:  A .LE. 3.686
*/
    b = 0.463+s-0.178*s2;
    si = 1.235;
    c = 0.195/s-7.9E-2+1.6E-2*s;
S40:
/*
     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
*/
    if(x <= 0.0) goto S70;
/*
     STEP  6:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S50;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S60;
S50:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S60:
/*
     STEP  7:  QUOTIENT ACCEPTANCE (Q)
*/
    if(log(1.0-u) <= q) return sgamma;
S70:
/*
     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
               U= 0,1 -UNIFORM DEVIATE
               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
*/
    e = sexpo();
    u = ranf();
    u += (u-1.0);
    t = b+fsign(si*e,u);
/*
     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
*/
    if(t < -0.7187449) goto S70;
/*
     STEP 10:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S80;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S90;
S80:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S90:
/*
     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
*/
    if(q <= 0.0) goto S70;
    if(q <= 0.5) goto S100;
    w = exp(q)-1.0;
    goto S110;
S100:
    w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
S110:
/*
               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
*/
    if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
    x = s+0.5*t;
    sgamma = x*x;
    return sgamma;
S120:
/*
     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
*/
    aa = 0.0;
    b = 1.0+0.3678794*a;
S130:
    p = b*ranf();
    if(p >= 1.0) goto S140;
    sgamma = exp(log(p)/ a);
    if(sexpo() < sgamma) goto S130;
    return sgamma;
S140:
    sgamma = -log((b-p)/ a);
    if(sexpo() < (1.0-a)*log(sgamma)) goto S130;
    return sgamma;
}


double snorm(void)
/*
**********************************************************************
                                                                      
                                                                      
     (STANDARD-)  N O R M A L  DISTRIBUTION                           
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             
               SAMPLING FROM THE NORMAL DISTRIBUTION.                 
               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          
                                                                      
     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  
     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  
                                                                      
     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
     SUNIF.  The argument IR thus goes away.                          
                                                                      
**********************************************************************
     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/
{
static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
};
static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
};
static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
};
static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
};
static long i;
static double snorm,u,s,ustar,aa,w,y,tt;
    u = ranf();
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (long) (u);
    if(i == 32) i = 31;
    if(i == 0) goto S100;
/*
                                START CENTER
*/
    ustar = u-(double)i;
    aa = *(a+i-1);
S40:
    if(ustar <= *(t+i-1)) goto S60;
    w = (ustar-*(t+i-1))**(h+i-1);
S50:
/*
                                EXIT   (BOTH CASES)
*/
    y = aa+w;
    snorm = y;
    if(s == 1.0) snorm = -y;
    return snorm;
S60:
/*
                                CENTER CONTINUED
*/
    u = ranf();
    w = u*(*(a+i)-aa);
    tt = (0.5*w+aa)*w;
    goto S80;
S70:
    tt = u;
    ustar = ranf();
S80:
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S70;
    ustar = ranf();
    goto S40;
S100:
/*
                                START TAIL
*/
    i = 6;
    aa = *(a+31);
    goto S120;
S110:
    aa += *(d+i-1);
    i += 1;
S120:
    u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
    w = u**(d+i-1);
    tt = (0.5*w+aa)*w;
    goto S160;
S150:
    tt = u;
S160:
    ustar = ranf();
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S150;
    u = ranf();
    goto S140;
}

double fsign( double num, double sign )
/* Transfers sign of argument sign to argument num */
{
if ( ( sign > 0.0 && num < 0.0) || ( sign < 0.0 && num > 0.0) )
    return -num;
else return num;
}



double sexpo(void)
/*
**********************************************************************
                                                                      
                                                                      
     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                
                                                                      
                                                                      
**********************************************************************
**********************************************************************
                                                                      
     FOR DETAILS SEE:                                                 
                                                                      
               AHRENS, J.H. AND DIETER, U.                            
               COMPUTER METHODS FOR SAMPLING FROM THE                 
               EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  
               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               
                                                                      
     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       
                                                                      
     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
     SUNIF.  The argument IR thus goes away.                          
                                                                      
**********************************************************************
     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
*/
{
static double q[8] = {
    0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,1.0
};
static long i;
static double sexpo,a,u,ustar,umin;
static double *q1 = q;
    a = 0.0;
    u = ranf();
    goto S30;
S20:
    a += *q1;
S30:
    u += u;
    if(u <= 1.0) goto S20;
    u -= 1.0;
    if(u > *q1) goto S60;
    sexpo = a+u;
    return sexpo;
S60:
    i = 1;
    ustar = ranf();
    umin = ustar;
S70:
    ustar = ranf();
    if(ustar < umin) umin = ustar;
    i += 1;
    if(u > *(q+i-1)) goto S70;
    sexpo = a+umin**q1;
    return sexpo;
}

long mltmod(long a,long s,long m)
/*
**********************************************************************
     long mltmod(long a,long s,long m)
                    Returns (A*S) MOD M
     This is a transcription from Pascal to Fortran of routine
     MULtMod_Decompos from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     a, s, m  -->
**********************************************************************
*/
{
#define h 32768L
static long mltmod,a0,a1,k,p,q,qh,rh;
/*
     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
      machine. On a different machine recompute H
*/
    if(!(a <= 0 || a >= m || s <= 0 || s >= m)) goto S10;
    puts(" a, m, s out of order in mltmod - ABORT!");
    printf(" a = %12ld s = %12ld m = %12ld\n",a,s,m);
    puts(" mltmod requires: 0 < a < m; 0 < s < m");
    exit(1);
S10:
    if(!(a < h)) goto S20;
    a0 = a;
    p = 0;
    goto S120;
S20:
    a1 = a/h;
    a0 = a-h*a1;
    qh = m/h;
    rh = m-h*qh;
    if(!(a1 >= h)) goto S50;
    a1 -= h;
    k = s/qh;
    p = h*(s-k*qh)-k*rh;
S30:
    if(!(p < 0)) goto S40;
    p += m;
    goto S30;
S40:
    goto S60;
S50:
    p = 0;
S60:
/*
     P = (A2*S*H)MOD M
*/
    if(!(a1 != 0)) goto S90;
    q = m/a1;
    k = s/q;
    p -= (k*(m-a1*q));
    if(p > 0) p -= m;
    p += (a1*(s-k*q));
S70:
    if(!(p < 0)) goto S80;
    p += m;
    goto S70;
S90:
S80:
    k = p/qh;
/*
     P = ((A2*H + A1)*S)MOD M
*/
    p = h*(p-k*qh)-k*rh;
S100:
    if(!(p < 0)) goto S110;
    p += m;
    goto S100;
S120:
S110:
    if(!(a0 != 0)) goto S150;
/*
     P = ((A2*H + A1)*H*S)MOD M
*/
    q = m/a0;
    k = s/q;
    p -= (k*(m-a0*q));
    if(p > 0) p -= m;
    p += (a0*(s-k*q));
S130:
    if(!(p < 0)) goto S140;
    p += m;
    goto S130;
S150:
S140:
    mltmod = p;
    return mltmod;
#undef h
}

double ranf(void)
/*
**********************************************************************
     double ranf(void)
                RANDom number generator as a Function
     Returns a random doubleing point number from a uniform distribution
     over 0 - 1 (endpoints of this interval are not returned) using the
     current generator
     This is a transcription from Pascal to Fortran of routine
     Uniform_01 from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
{
static double ranf;
/*
     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
      and is currently 2147483563. If M1 changes, change this also.
*/
    ranf = ignlgi()*4.656613057E-10;
    return ranf;
}

void gscgn(long getset,long *g)
/*
**********************************************************************
     void gscgn(long getset,long *g)
                         Get/Set GeNerator
     Gets or returns in G the number of the current generator
                              Arguments
     getset --> 0 Get
                1 Set
     g <-- Number of the current random number generator (1..32)
**********************************************************************
*/
{
#define numg 32L
static long curntg = 1;
    if(getset == 0) *g = curntg;
    else  {
        if(*g < 0 || *g > numg) {
            puts(" Generator number out of range in GSCGN");
            exit(0);
        }
        curntg = *g;
    }
#undef numg
}

void gsrgs(long getset,long *qvalue)
/*
**********************************************************************
     void gsrgs(long getset,long *qvalue)
               Get/Set Random Generators Set
     Gets or sets whether random generators set (initialized).
     Initially (data statement) state is not set
     If getset is 1 state is set to qvalue
     If getset is 0 state returned in qvalue
**********************************************************************
*/
{
static long qinit = 0;

    if(getset == 0) *qvalue = qinit;
    else qinit = *qvalue;
}

void gssst(long getset,long *qset)
/*
**********************************************************************
     void gssst(long getset,long *qset)
          Get or Set whether Seed is Set
     Initialize to Seed not Set
     If getset is 1 sets state to Seed Set
     If getset is 0 returns T in qset if Seed Set
     Else returns F in qset
**********************************************************************
*/
{
static long qstate = 0;
    if(getset != 0) qstate = 1;
    else  *qset = qstate;
}

void setall(long iseed1,long iseed2)
/*
**********************************************************************
     void setall(long iseed1,long iseed2)
               SET ALL random number generators
     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
     initial seeds of the other generators are set accordingly, and
     all generators states are set to these seeds.
     This is a transcription from Pascal to Fortran of routine
     Set_Initial_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First of two integer seeds
     iseed2 -> Second of two integer seeds
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern void gssst(long getset,long *qset);
extern void gscgn(long getset,long *g);
extern long Xm1,Xm2,Xa1vw,Xa2vw,Xig1[],Xig2[];
static long T1;
static long g,ocgn;
static long qrgnin;
    T1 = 1;
/*
     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
      HAS BEEN CALLED.
*/
    gssst(1,&T1);
    gscgn(0L,&ocgn);
/*
     Initialize Common Block if Necessary
*/
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    *Xig1 = iseed1; 
    *Xig2 = iseed2;
    initgn(-1L);
    for(g=2; g<=numg; g++) {
        *(Xig1+g-1) = mltmod(Xa1vw,*(Xig1+g-2),Xm1);
        *(Xig2+g-1) = mltmod(Xa2vw,*(Xig2+g-2),Xm2);
        gscgn(1L,&g);
        initgn(-1L);
    }
    gscgn(1L,&ocgn);
#undef numg
}

void initgn(long isdtyp)
/*
**********************************************************************
     void initgn(long isdtyp)
          INIT-ialize current G-e-N-erator
     Reinitializes the state of the current generator
     This is a transcription from Pascal to Fortran of routine
     Init_Generator from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     isdtyp -> The state to which the generator is to be set
          isdtyp = -1  => sets the seeds to their initial value
          isdtyp =  0  => sets the seeds to the first value of
                          the current block
          isdtyp =  1  => sets the seeds to the first value of
                          the next block
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern void gscgn(long getset,long *g);
extern long Xm1,Xm2,Xa1w,Xa2w,Xig1[],Xig2[],Xlg1[],Xlg2[],Xcg1[],Xcg2[];
static long g;
static long qrgnin;
/*
     Abort unless random number generator initialized
*/
    gsrgs(0L,&qrgnin);
    if(qrgnin) goto S10;
    printf(
      " INITGN called before random number generator  initialized -- abort!\n");
    exit(1);
S10:
    gscgn(0L,&g);
    if(-1 != isdtyp) goto S20;
    *(Xlg1+g-1) = *(Xig1+g-1);
    *(Xlg2+g-1) = *(Xig2+g-1);
    goto S50;
S20:
    if(0 != isdtyp) goto S30;
    goto S50;
S30:
/*
     do nothing
*/
    if(1 != isdtyp) goto S40;
    *(Xlg1+g-1) = mltmod(Xa1w,*(Xlg1+g-1),Xm1);
    *(Xlg2+g-1) = mltmod(Xa2w,*(Xlg2+g-1),Xm2);
    goto S50;
S40:
    printf("isdtyp not in range in INITGN");
    exit(1);
S50:
    *(Xcg1+g-1) = *(Xlg1+g-1);
    *(Xcg2+g-1) = *(Xlg2+g-1);
#undef numg
}


long ignlgi(void)
/*
**********************************************************************
     long ignlgi(void)
               GeNerate LarGe Integer
     Returns a random integer following a uniform distribution over
     (1, 2147483562) using the current generator.
     This is a transcription from Pascal to Fortran of routine
     Random from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern void gssst(long getset,long *qset);
extern void gscgn(long getset,long *g);
extern long Xm1,Xm2,Xa1,Xa2,Xcg1[],Xcg2[];
extern long Xqanti[];
static long ignlgi,curntg,k,s1,s2,z;
static long qqssd,qrgnin;
/*
     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
     THIS ROUTINE  2) A CALL TO SETALL.
*/
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    gssst(0,&qqssd);
    if(!qqssd) setall(1234567890L,123456789L);
/*
     Get Current Generator
*/
    gscgn(0L,&curntg);
    s1 = *(Xcg1+curntg-1);
    s2 = *(Xcg2+curntg-1);
    k = s1/53668L;
    s1 = Xa1*(s1-k*53668L)-k*12211;
    if(s1 < 0) s1 += Xm1;
    k = s2/52774L;
    s2 = Xa2*(s2-k*52774L)-k*3791;
    if(s2 < 0) s2 += Xm2;
    *(Xcg1+curntg-1) = s1;
    *(Xcg2+curntg-1) = s2;
    z = s1-s2;
    if(z < 1) z += (Xm1-1);
    if(*(Xqanti+curntg-1)) z = Xm1-z;
    ignlgi = z;
    return ignlgi;
#undef numg
}

void inrgcm(void)
/*
**********************************************************************
     void inrgcm(void)
          INitialize Random number Generator CoMmon
                              Function
     Initializes common area  for random number  generator.  This saves
     the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
     assuring that the routine is loaded with the other routines.
**********************************************************************
*/
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern long Xm1,Xm2,Xa1,Xa2,Xa1w,Xa2w,Xa1vw,Xa2vw;
extern long Xqanti[];
static long T1;
static long i;
/*
     V=20;                            W=30;
     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
    An efficient way to precompute a**(2*j) MOD m is to start with
    a and square it j times modulo m using the function MLTMOD.
*/
    Xm1 = 2147483563L;
    Xm2 = 2147483399L;
    Xa1 = 40014L;
    Xa2 = 40692L;
    Xa1w = 1033780774L;
    Xa2w = 1494757890L;
    Xa1vw = 2082007225L;
    Xa2vw = 784306273L;
    for(i=0; i<numg; i++) *(Xqanti+i) = 0;
    T1 = 1;
/*
     Tell the world that common has been initialized
*/
    gsrgs(1L,&T1);
#undef numg
}




/************************************************************************
                      INTEGRATION
************************************************************************/

double midpnt(double (*func)(double), double a, double b, int n)
/* This routine computes the nth stage of refinement of an extended midpoint rule. func is input
   as a pointer to the function to be integrated between limits a and b. When called with n=1, 
   the routine return the crudest estimate. Subsequent calls with n=2,3... (in that order) will
   improve the accuracy of s by adding (2/3)*3^(n-1) additional interior points.
   s should not be modified between sequential calls
*/
{

  #define FUNC(x) ((*func)(x))

  double x, tnm, sum, del, ddel;
  static double s;
  int it, j;

  if (n==1) {
    return(s=(b-a)*FUNC(0.5*(a+b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC(x);
      x += ddel;
      sum += FUNC(x);
      x += del;
    }
    return(s= (s+(b-a)*sum/tnm)/3.0);
  }

}

double midinf(double (*funk)(double), double aa, double bb, int n)
/* This routine is an exact replacement for midpnt i.e. returns the nth stage of refinement of
   the integral of funk from aa to bb, except that the function is evaluated at evenly spaced
   points in 1/x rather than in x. This allows the upper limit bb to be as large and positive as
   the computer allows, or the lower limit aa to be as large and negative, but not both.
   aa and bb must have the same sign, and they cannot be equal to zero.
*/
{

  #define FUNK(x) ((*funk)(1.0/(x))/((x)*(x)))
  double x,tnm,sum,del,ddel,b,a;
  static double s;
  int it,j;

  b= 1.0/aa;  //change the limits of integration
  a= 1.0/bb;

  if (n==1) {  //from here on the routine is identical to midpnt
    return(s=(b-a)*FUNK(0.5*(a+b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNK(x);
      x += ddel;
      sum += FUNK(x);
      x += del;
    }
    return(s= (s+(b-a)*sum/tnm)/3.0);
  }

}

double qromo(double (*func)(double), double a, double b, double (*choose)(double(*)(double), double, double, int)) {
/* Romberg integration on an open interval. Returns the integral of the function func from a to b,
   using any specified integrating function choose and Romberg's method. Normally choose will be an
   open formula like midpnt or midinf, not evaluating the function at the endpoints. It is assumed that choose
   triples the number of steps on each call, and that its error series contains only even powers of
   the number of steps.
   Integration is performed by Romberg's method of order 2K (K=2 is Simpson's rule)

   USAGE EXAMPLES
   answer= qromo(f,0.0,2.0,midpnt); //integrate f from 0 to 2
   answer= qromo(f,0.0,2.0,midpnt) + qromo(f,2.0,1.0e30,midinf) //integrate f from 0 to 1.0e30 (cutpoint should be chosen in the tail of fx)

 */

#define EPS 1.0e-6
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

  void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
  int j;
  double ss, dss, h[JMAXP+1], s[JMAXP];

  h[1]= 1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=(*choose)(func,a,b,j);
    if (j >= K) {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) <= EPS*fabs(ss)) return ss;
    }
    h[j+1]= h[j]/9.0;  //this is where the assumption of step tripling and even error series is used
  }
  nrerror("qromo","integrate a function","");
  return(0.0); //never get here

}


/************************************************************************
          INTERPOLATION, EXTRAPOLATION AND SPLINES
************************************************************************/

void polint (double xa[], double ya[], int n, double x, double *y, double *dy)
/* Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and
   and error estimate dy. If P(x) is the polynomial of degree N-1 such that P(xa[i])=ya[i] i=1..n
   then the returned value y=P(x).
*/
{

  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c, *d;

  dif= fabs(x-xa[1]);
  c= dvector(1,n);
  d= dvector(1,n);

  for (i=1; i<=n; i++) { //here we find the index ns of the closest table entry
    if ((dift=fabs(x-xa[i])) < dif) { ns=i; dif= dift; }
    c[i]= ya[i];
    d[i]= ya[i];
  }
  *y= ya[ns--];          //this is the initial approximation to y
  for (m=1; m<n; m++) {  //for each column of the tableau
    for (i=1; i<=n-m; i++) {  //we loop over the current c's and d's and update them
      ho= xa[i]-x;
      hp=xa[i+m]-x;
      w= c[i+1] - d[i];
      if ( (den=ho-hp) == 0.0) nrerror("polint","increment in x axis in 0 units (two input x values are identical)","");
      den= w/den;
      d[i]= hp*den;
      c[i]= ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    //after each column in the tableau is completed, we decide which correction (c or d) to add to our accumulating value of y
    //i.e. which path to take through the tableau: forking up or down. We do this in such a way as to take the most "straigh line"
    //route to its apex, updating ns accordingly to keep track of where we are. This route keeps the partial approximations
    //centered on the target x. The last dy added is thus the error indication.
  }

  free_dvector(d,1,n);
  free_dvector(c,1,n);

}


double bspline_singlex(double x, int j, int degree, double *knots) {
/*Returns the jth B-spline basis evaluated at single value x
  - x: value at which to evaluate the B-spline basis
  - j: basis
  - degree: degree of the B-spline (0: piecewise constant, 1:linear etc.)
  - knots: sequence of knots
*/
  double ans;
  if (degree==0) {
    if (knots[j]<=x && x<knots[j+1]) ans= 1.0; else ans= 0.0;
  } else {
    ans= bspline_singlex(x,j,degree-1,knots)*(x-knots[j])/(knots[j+degree]-knots[j]) + bspline_singlex(x,j+1,degree-1,knots)*(knots[j+degree+1]-x)/(knots[j+degree+1]-knots[j+1]);
  }
  return(ans);
}

void bspline(double **W, double *x, int *nx, int *degree, double *knots, int *nknots) {
 //B-spline basis eval at vector of values x. Normalized to sum to 1 at any x value.
 /* Input
      - x: vector of values at which to evaluate the B-spline basis
      - nx: length of x
      - degree: degree of the spline (0: piecewise constant, 1: linear etc.)
      - knots: vector with positions of the knots
      - nknots: length of knots
    Output
      - W: matrix with nx rows and nknots-degree-1 columns containing the B-spline basis
   */
  int i,j;
  if (*nknots<(*degree+2)) {
    printf("error: number of knots must be >= degree + 2");
  } else {
    for (i=0; i<(*nx); i++) {
      for (j=0; j<(*nknots - *degree -1); j++) {
        W[i][j]= bspline_singlex(x[i],j,*degree,knots);
      }
    }
  }
}

void bspline_vec(double *W, double *x, int *nx, int *degree, double *knots, int *nknots) {
  //same as routine bspline but returs a vector so that it can be called from R
  int i,j;
  double **Wtemp;

Wtemp= dmatrix(0,*nx,0,*nknots- *degree -1);
bspline(Wtemp,x,nx,degree,knots,nknots);
for (i=0; i<(*nx); i++) { for (j=0; j<(*nknots - *degree -1); j++) { W[i*(*nknots - *degree -1)+j]= Wtemp[i][j]; } }
free_dmatrix(Wtemp,0,*nx,0,*nknots- *degree -1);
}

void mspline(double **W, double *x, int *nx, int *degree, double *knots, int *nknots) {
  //M-spline basis eval at vector of values x. Normalized to integrate to 1 wrt x
  int i,j;
  if (*nknots<(*degree+2)) {
    printf("error: number of knots must be >= degree + 2");
  } else {
    for (i=0; i<(*nx); i++) {
      for (j=0; j<(*nknots - *degree -1); j++) {
        W[i][j]= bspline_singlex(x[i],j,*degree,knots)*(*degree+1.0)/(knots[j+ *degree +1]-knots[j]);
      }
    }
  }
}

void mspline_vec(double *W, double *x, int *nx, int *degree, double *knots, int *nknots) {
  //same as routine mspline but returs a vector so that it can be called from R
  int i,j;
  double **Wtemp;

Wtemp= dmatrix(0,*nx,0,*nknots- *degree -1);
mspline(Wtemp,x,nx,degree,knots,nknots);
for (i=0; i<(*nx); i++) { for (j=0; j<(*nknots - *degree -1); j++) { W[i*(*nknots - *degree -1)+j]= Wtemp[i][j]; } }
free_dmatrix(Wtemp,0,*nx,0,*nknots- *degree -1);
}


/************************************************************************
            FUNCTION OPTIMIZATION 
************************************************************************/

double univmin(double ax, double bx, double cx, double (*f)(double), double eps, double *xmin, int itmax) {
/* PURPOSE: UNIVARIATE OPTIMIZATION WITHOUT DERIVATIVE INFORMATION
   Given a function f and a bracketing triplet of abscissas ax, bx, cx (such that bx is between ax & cx
   and f(bx) is less than f(ax) and f(cx)), this routine isolates the minimum to a fractional precision
   of about eps using Brent's method. The abscissa of the minimum is returned as xmin and the value
   of the function at the minimum is the returned value.
 */
#define CGOLD .3819660  //golden ratio
#define ZEPS 1.0e-10   //protects against trying to achieve fractional accuracy when min is exactly zero
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b)>=0.0 ? fabs(a) : -fabs(a))

  int iter;
  double a,b,d=1,etemp,fu,fv,fw,fx,p,q,r,eps1,eps2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx);     //a,b must be in ascending order but input abscissas need not be
  b=(ax > cx ? ax : cx);

  x=w=v=bx;                  //initializations
  fw=fv=fx=(*f)(x);

  for (iter=1; iter<=itmax; iter++) {      //main program loop
    xm=0.5*(a+b);
    eps2= 2.0*(eps1=eps*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (eps2-0.5*(b-a))) {  //test for done here
      *xmin= x;
      return fx;
    }
    if (fabs(e) > eps1) {                  //construct a trial parabolic fit
      r= (x-w)*(fx-fv);
      q= (x-v)*(fx-fw);
      p= (x-v)*q - (x-w)*r;
      q= 2.0*(q-r);
      if (q>0.0) p= -p;
      q= fabs(q);
      etemp= e;
      e= d;
      //determine acceptability of parabolic fit & take the golden section step into the larger of the two segments
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) { 
	d= CGOLD*(e=(x >= xm ? a-x : b-x));
      } else {
	d=p/q;                            //take the parabolic step
	u= x+d;
	if (u-a < eps2 || b-u < eps2) d=SIGN(eps1,xm-x);
      }
    } else {
      d= CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u= (fabs(d) >= eps1 ? x+d : x+SIGN(eps1,d));
    fu= (*f)(u);                         //this is the one function evaluation per iteration
    if (fu <= fx) {                      //now decide what to do with our function evaluation
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u);                     //housekeeping
      SHFT(fv,fw,fx,fu);
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    }                                   //done with housekeeping. Back for another iteration
  }

  *xmin= x;                             //only get here if iteration limit is reached
  return fx;

}

#define MOV3(a,b,c,d,e,f) (a)=(d);(b)=(e);(c)=(f);
double dunivmin(double ax,double bx,double cx,double (*f)(double),double (*df)(double),double eps,double *xmin,int itmax) {
/* Given a function f and its derivative function df, and given a bracketing triplet of abscissas ax,
   bx, cx [such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)],
   this routine isolates the minimum to a fractional precision of about eps using a modification of
   Brent's method that uses derivatives. The abscissa of the minimum is returned as xmin, and
the minimum function value is returned as dunivmin, the returned function value.
*/

  #define ZEPS 1.0e-10   //protects against trying to achieve fractional accuracy when min is exactly zero
  int iter,ok1,ok2; //Will be used as flags for whether proposed steps are acceptable or not.
  double a,b,d=1,d1,d2,du,dv,dw,dx,e=0.0; 
  double fu,fv,fw,fx,olde,eps1,eps2,u,u1,u2,v,w,x,xm;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  dw=dv=dx=(*df)(x); 
  for (iter=1;iter<=itmax;iter++) {
    xm=0.5*(a+b);
    eps1=eps*fabs(x)+ZEPS;
    eps2=2.0*eps1;
    if (fabs(x-xm) <= (eps2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > eps1) {
      d1=2.0*(b-a); //Initialize these d's to an out-of-bracket value
      d2=d1;
      if (dw != dx) d1=(w-x)*dx/(dx-dw); //Secant method with one point.
      if (dv != dx) d2=(v-x)*dx/(dx-dv); //And the other.
      //Which of these two estimates of d shall we take? We will insist that they be within
      //the bracket, and on the side pointed to by the derivative at x:
      u1=x+d1;
      u2=x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde=e; //Movement on the step before last.
      e=d;
      if (ok1 || ok2) { //Take only an acceptable d, and if both are acceptable, then take the smallest one.
	if (ok1 && ok2)
	  d=(fabs(d1) < fabs(d2) ? d1 : d2);
	else if (ok1)
	  d=d1;
	else
	  d=d2;
	if (fabs(d) <= fabs(0.5*olde)) {
	  u=x+d;
	  if (u-a < eps2 || b-u < eps2)
	    d=SIGN(eps1,xm-x);
	} else { //Bisect, not golden section.
	  d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	  //Decide which segment by the sign of the derivative.
	}
      } else {
	d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
    } else {
      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
    }
    if (fabs(d) >= eps1) {
      u=x+d;
      fu=(*f)(u);
    } else {
      u=x+SIGN(eps1,d);
      fu=(*f)(u);
      if (fu > fx) { //If the minimum step in the downhill direction takes us uphill, then we are done.
	*xmin=x;
	return fx;
      }
    }
    du=(*df)(u); //Now all the housekeeping, sigh.
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      MOV3(v,fv,dv, w,fw,dw);
      MOV3(w,fw,dw, x,fx,dx);
      MOV3(x,fx,dx, u,fu,du);
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	MOV3(v,fv,dv, w,fw,dw);
	MOV3(w,fw,dw, u,fu,du);
      } else if (fu < fv || v == x || v == w) {
	MOV3(v,fv,dv, u,fu,du);
      }
    }
  }

  *xmin= x;                             //only get here if iteration limit is reached
  return fx;

}



void minimize(double th[],double **dirini,int n,double eps,int *iter,double *fret,double (*f)(double []),int itmax) {
/* Multivariate function minimization. */
/* Input/Output
   - th[1..n]: initial parameter values.
   - dirini[1..n][1..n]: matrix with initial directions, typically the n unit vectors
   - n: length of th
   - eps: relative tolerance to achieve convergence
   - f: function to minimize (must take a vector of doubles as input)
   Output
   - th[1..n]: optimal parameter values
   - dirini[1..n][1..n]: directions used in the last iteration
   - iter: number of iterations used
   - fret: value of the function at the optimum
*/

  int i,ibig,j,converged=0;
  double del,fth,fthtt,t,*tht,*thtt,*dirinit;

  tht=dvector(1,n);
  thtt=dvector(1,n);
  dirinit=dvector(1,n);
  //initial parameter and function values
  *fret=(*f)(th);
  for (j=1;j<=n;j++) tht[j]=th[j];
  for (*iter=1;(*iter <itmax) && (converged==0);++(*iter)) {
    fth=(*fret);
    ibig=0;
    del=0.0;
    for (i=1;i<=n;i++) {  //minimize along all directions
      for (j=1;j<=n;j++) dirinit[j]=dirini[j][i];
      fthtt=(*fret);
      dirmin(th,dirinit,n,fret,f,itmax,eps);
      if (fabs(fthtt-(*fret)) > del) {  //
	del=fabs(fthtt-(*fret));
	ibig=i;
      }
    }
    for (j=1;j<=n;j++) { //extrapolated point, average direction and function value at extrapolated point
      thtt[j]=2.0*th[j]-tht[j];
      dirinit[j]=th[j]-tht[j];
      tht[j]=th[j];
    }
    fthtt=(*f)(thtt);
    if (fthtt < fth) {
      t=2.0*(fth-2.0*(*fret)+fthtt)*sqrt(fth-(*fret)-del)-del*sqrt(fth-fthtt);
      if (t < 0.0) {
	dirmin(th,dirinit,n,fret,f,itmax,eps);
	for (j=1;j<=n;j++) {
	  dirini[j][ibig]=dirini[j][n];
	  dirini[j][n]=dirinit[j];
	}
      }
    }
    if (2.0*fabs(fth-(*fret)) <= eps*(fabs(fth)+fabs(*fret))) converged=1;
  }

  free_dvector(dirinit,1,n);
  free_dvector(thtt,1,n);
  free_dvector(tht,1,n);
}


#define TINY 1.0e-25
int ncom; //Global variables communicate with f1dim.
double *pcom,*xicom,(*nrfunc)(double []);

void dirmin(double p[], double xi[], int n, double *fret, double (*func)(double []), int itmax, double dirminEPS) {
/* Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
   resets p to where the function func(p) takes on a minimum along the direction xi from p,
   and replaces xi by the actual vector displacement that p was moved. Also returns as fret
   the value of func at the returned location p. This is actually all accomplished by calling the
   routines mnbrak and univmin.
*/

double univmin(double ax, double bx, double cx,double (*f)(double), double eps, double *xmin, int itmax);
double f1dim(double x);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
double *fc, double (*func)(double));
int j;
double xx,xmin,fx,fb,fa,bx,ax;

 ncom=n; //Define the global variables.
 pcom=dvector(1,n);
 xicom=dvector(1,n);
 nrfunc=func;
 for (j=1;j<=n;j++) {
   pcom[j]=p[j];
   xicom[j]=xi[j];
 }
 ax=0.0; //Initial guess for brackets.
 xx=1.0;
 mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
 *fret=univmin(ax,xx,bx,f1dim,dirminEPS,&xmin,itmax);
 for (j=1;j<=n;j++) { //Construct the vector results to return.
   xi[j] *= xmin;
   p[j] += xi[j];
 }
 free_dvector(xicom,1,n);
 free_dvector(pcom,1,n);
}

double f1dim(double x) {
//Must accompany dirmin.
  int j;
  double f,*xt;
  xt=dvector(1,ncom);
  for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  free_dvector(xt,1,ncom);
  return f;
}

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,double (*func)(double)) {
/* Given a function func, and given distinct initial points ax and bx, this routine searches in
   the downhill direction (defined by the function as evaluated at the initial points) and returns
   new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
   values at the three points, fa, fb, and fc.
*/
#define GOLD 1.618034  //default ratio by which successive intervals are magnified;
#define GLIMIT 100.0   //maximum magnification allowed for a parabolic-fit step.

double ulim,u,r,q,fu,dum;

 *fa=(*func)(*ax);
 *fb=(*func)(*bx);
 if (*fb > *fa) { //Switch roles of a and b so that we can go downhill in the direction from a to b.
   SHFT(dum,*ax,*bx,dum);
   SHFT(dum,*fb,*fa,dum);
 }
 *cx=(*bx)+GOLD*(*bx-*ax); //First guess for c.
 *fc=(*func)(*cx);
 while (*fb > *fc) {       //Keep returning here until we bracket.
   r=(*bx-*ax)*(*fb-*fc);  //Compute u by parabolic extrapolation from a, b, c. TINY prevents any division by zero.
   q=(*bx-*cx)*(*fb-*fa);
   u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
   ulim=(*bx)+GLIMIT*(*cx-*bx); //We won't go farther than this. Test various possibilities:
   if ((*bx-u)*(u-*cx) > 0.0) { //Parabolic u is between b and c: try it.
     fu=(*func)(u);
     if (fu < *fc) {              //Got a minimum between b and c.
       *ax=(*bx);
       *bx=u;
       *fa=(*fb);
       *fb=fu;
       return;
     } else if (fu > *fb) {       //Got a minimum between between a and u.
       *cx=u;
       *fc=fu;
       return;
     }
     u=(*cx)+GOLD*(*cx-*bx); //Parabolic fit was no use. Use default magnification.
     fu=(*func)(u);
   } else if ((*cx-u)*(u-ulim) > 0.0) { //Parabolic fit is between c and its allowed limit
     fu=(*func)(u);
     if (fu < *fc) {
       SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
       SHFT(*fb,*fc,fu,(*func)(u));
     }
   } else if ((u-ulim)*(ulim-*cx) >= 0.0) { //Limit parabolic u to maximum allowed value
     u=ulim;
     fu=(*func)(u);
   } else {                                 //Reject parabolic u, use default magnification
     u=(*cx)+GOLD*(*cx-*bx);
     fu=(*func)(u);
   }
   SHFT(*ax,*bx,*cx,u); //Eliminate oldest point and continue.
   SHFT(*fa,*fb,*fc,fu);
 }
}
