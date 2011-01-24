#include <R.h>
#include <Rinternals.h>

//Define structures
struct marginalPars {
  int *sel;
  int *nsel;
  int *n;
  int *p;
  double *y;
  double *sumy2;
  double *x;
  double *XtX;
  double *ytX;
  double *m;  //Sinv * Xty   (needed by mom and emom)
  double **S;  //XtX + I/tau  (needed by mom and emom)
  int *method;
  int *B;
  double *alpha;    //prior for residual variance is IG(.5*alpha,.5*lambda)
  double *lambda;
  double *phi;      //residual variance
  double *tau;      //coefficients prior dispersion parameter
  int *r;           //MOM power parameter for prior on coefficients
  double *prDeltap; //For Binomial prior on model space, prDeltap is the prob of success
  int *logscale;
};

typedef double(*pt2margFun)(int *, int *, struct marginalPars *);  //pointer to function to compute marginal densities & prior prob

/*
struct optpars { double *XtX; double *ytX; double *phi; double *tau; int *n; int *p; int *sel; int *nsel; };
struct phiintpars { double *XtX; double *ytX; double *tau; int *n; int *p; int *sel; int *nsel; double *y; double *sumy2; int method;  int *B; double *alpha; double *lambda; };
*/

//General Algebra
void Asym_xsel(double **A, int fi, double *x, int *sel, double *ans);  //multiply symmetric A[1..fi][1..fi] * x[sel[0]..sel[fi-1]]; Return in ans[1..fi]
void addct2XtX(double *ct, double *XtX, int *sel, int *nsel, int *p, double **V); //add constant to diagonal elem of XtX

//General marginal density calculation routines
void set_marginalPars(struct marginalPars *pars, int *n,int *p,double *y,double *sumy2,double *x,double *XtX,double *ytX,int *method,int *B,double *alpha,double *lambda,double *phi,double *tau,int *r,double *prDeltap, int *logscale);
void set_f2opt_pars(double *m, double **S, double *XtX, double *ytX, double *phi, double *tau, int *r, int *n, int *p, int *sel, int *nsel);
void set_f2int_pars(double *XtX, double *ytX, double *tau, int *n, int *p, int *sel, int *nsel, double *y, double *sumy2, int *method, int *B, double *alpha, double *lambda, int *logscale);

//*************************************************************************************
// Model Selection Routines
//*************************************************************************************

void modelSelectionC(int *postSample, double *margpp, int *postMode, double *postModeProb, double *postProb, int *knownphi, int *prCoef, int *prDelta, int *niter, int *thinning, int *burnin, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars);
void greedyVarSelC(int *postMode, double *postModeProb, int *knownphi, int *prCoef, int *prDelta, int *niter, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars);
pt2margFun set_marginalFunction(int *prCoef, int *knownphi);
pt2margFun set_priorFunction(int *prDelta);
void sel2selnew(int newelem,int *sel,int *nsel,int *selnew,int *nselnew);

// Priors on Model Space (always return on log scale)
double unifPrior(int *sel, int *nsel, struct marginalPars *pars);
double binomPrior(int *sel, int *nsel, struct marginalPars *pars);


//*************************************************************************************
// Product MOM routines
//*************************************************************************************

double f2opt_mom(double *th);
double fmomNegC_non0(double *th, double *m, double **S, double *phi, double *tau, int *r, int *n, int *nsel);
void fppmomNegC_non0(double **ans, double *th, double **S, double *phi, double *tau, int *r, int *n, int *nsel);
void momIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *n, int *nsel, double *m, double **S, double *detS, double *phi, double *tau, int *r, int *logscale);

double rsumlogsq(double *th, int *r, int *nsel);  //compute r*sum(log(th^2))
double pmomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double MC_mom(double *m,double **Sinv,int *r,int *nsel, int *B);  //MC evaluation of E(prod(th^2r)) for th ~ N(m,Sinv)
double MC_mom_T(double *m,double **Sinv,int *nu,int *r,int *nsel, int *B); //MC evaluation of E(prod(th^2r)) for th ~ T_nu(m,Sinv)

//*************************************************************************************
// Product iMOM routines
//*************************************************************************************

double f2opt_imom(double *th);
double fimomNegC(double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel);
double fimomNegC_non0(double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel);
//void fimomNegC_non0_mat(double *ans, int *nrow, double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel);
void fppimomNegC_non0(double **ans, double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel);
void imomIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *sel, int *nsel, int *n, int *p, double *XtX, double *ytX, double *phi, double *tau, int *logscale);

SEXP pimomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale);
double pimomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double IS_imom(double *thopt, double **Voptinv, int *sel, int *nsel, int *n, int *p, double *XtX, double *ytX, double *phi, double *tau, int *B);

double f2int_imom(double phi);
SEXP pimomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda);
double pimomMarginalUC(int *sel, int *nsel, struct marginalPars *pars);


//*************************************************************************************
// Product MOM routines
//*************************************************************************************

double pmomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double pmomMarginalUC(int *sel, int *nsel, struct marginalPars *pars);


//*************************************************************************************
// Product eMOM routines
//*************************************************************************************

double pemomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double pemomMarginalUC(int *sel, int *nsel, struct marginalPars *pars);
