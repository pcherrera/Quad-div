#include "mex.h"
#include<math.h>

#include "mexHelpers/funSigmaTau.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) {
  
 
  double * a = mxGetPr(prhs[0]);
  double * b = mxGetPr(prhs[1]);
  double * c = mxGetPr(prhs[2]);
  double * d = mxGetPr(prhs[3]);
  
  int p = (int) mxGetScalar(prhs[4]);
  int q = (int) mxGetScalar(prhs[5]);

  /* number of elements */
  int nE = mxGetM(prhs[0]);
  
  int dimP = (p+1)*(p+2);
  int dimQ = (q+1)*(q+2);
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(dimP*dimQ*nE,1,mxREAL);
  
  /* Get pointer to Matrix */
  double * M = mxGetPr(plhs[0]);
  
  /* Define a function pointer */
  void (*compFun)(double*,double,double,double,double) ;
  
  /* Depending on p and q choose the right function */
  if(p==0) {
      if(q==1)
          compFun = &compTsigmaTauP0Q1;
      else if(q==2)
          compFun = &compTsigmaTauP0Q2;
      else if(q==3)
          compFun = &compTsigmaTauP0Q3;
      else if(q==4)
          compFun = &compTsigmaTauP0Q4;
      else if(q==5)
          compFun = &compTsigmaTauP0Q5;
  }                    
  else if(p==1) {
      if(q==1)
          compFun = &compTsigmaTauP1Q1;
      else if(q==2)
          compFun = &compTsigmaTauP1Q2;
      else if(q==3)
          compFun = &compTsigmaTauP1Q3;
      else if(q==4)
          compFun = &compTsigmaTauP1Q4;
      else if(q==5)
          compFun = &compTsigmaTauP1Q5;
  }
  else if(p==2) {
      if(q==1)
          compFun = &compTsigmaTauP2Q1;
      else if(q==2)
          compFun = &compTsigmaTauP2Q2;
      else if(q==3)
          compFun = &compTsigmaTauP2Q3;
      else if(q==4)
          compFun = &compTsigmaTauP2Q4;
      else if(q==5)
          compFun = &compTsigmaTauP2Q5;
  }
  else if(p==3) {
      if(q==1)
          compFun = &compTsigmaTauP3Q1;
      else if(q==2)
          compFun = &compTsigmaTauP3Q2;
      else if(q==3)
          compFun = &compTsigmaTauP3Q3;
      else if(q==4)
          compFun = &compTsigmaTauP3Q4;
      else if(q==5)
          compFun = &compTsigmaTauP3Q5;
  }
      
  
  //mexPrintf("p = %d, q = %d, dimP = %d, dimQ = %d\n",p,q,dimP,dimQ);
  
  /* for loop */
  int i = 0;
 
  #pragma omp parallel for
  for(i=0;i<nE; ++i) {
      //compMtmpP0Q2(M+i*dimP*dimQ,a[i],b[i],c[i],d[i]);
      compFun(M+i*dimP*dimQ,a[i],b[i],c[i],d[i]);
  }
 
}