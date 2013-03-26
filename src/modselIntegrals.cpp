#include "modselIntegrals.h"
#include "cstat.h"
using namespace std;

modselIntegrals::modselIntegrals(pt2margFun marfun, pt2margFun priorfun, int nvars) {
  int i;

  this->maxVars= nvars;
  this->marginalFunction= marfun;
  this->priorFunction= priorfun;

  this->zerochar = (char *) calloc(nvars, sizeof(char));
  for (i=0; i<nvars; i++) this->zerochar[i]= '0';
  //this->zerochar= charvector(0,maxVars-1);
}

modselIntegrals::~modselIntegrals() {

  free((char  *) this->zerochar);
  //free_charvector(this->zerochar, 0, maxVars-1);

}

//Return log(marginal likelihood) + log(prior). Uses logjointSaved if available, else adds result to logjointSaved
// Input: 
//
//   - sel: integer vector [0..maxVars-1], where 0's and 1's indicate covariates out/in the model (respectively)
//   - nsel: number of covariates in the model (i.e. sum(sel))
//   - pars: struct of type marginalPars containing parameters needed to evaluate the marginal density of the data & prior on model space
//
// Output: evaluates log joint. It returns previously saved results in logjointSaved if available, else it performs the computation and saves the result in logjointSaved
double modselIntegrals::getJoint(int *sel, int *nsel, struct marginalPars *pars) {
  int i;
  double ans;

  for (i=0; i< *nsel; i++) zerochar[sel[i]]= '1';
  std::string s (zerochar);
  
  if (logjointSaved.count(s) > 0) { 
    ans= logjointSaved[s];
  } else {
    ans= marginalFunction(sel,nsel,pars) + priorFunction(sel,nsel,pars);
    logjointSaved[s]= ans;
  }

  for (i=0; i<= *nsel; i++) this->zerochar[sel[i]]= '0';

  return ans;
} 
