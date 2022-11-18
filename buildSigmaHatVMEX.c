#include "mex.h"
#include<math.h>

#include "mexHelpers/funSigmaHatV.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) {
  
 
  double * e1 = mxGetPr(prhs[0]);
  double * e2 = mxGetPr(prhs[1]);
  double * e3 = mxGetPr(prhs[2]);
  double * o1 = mxGetPr(prhs[3]);
  double * o2 = mxGetPr(prhs[4]);
  double * o3 = mxGetPr(prhs[5]);
  
  int p = (int) mxGetScalar(prhs[6]);
  int q = (int) mxGetScalar(prhs[7]);

  /* number of elements */
  int nE = mxGetM(prhs[0]);
  
  int dimP = (p+1)*3;
  int dimQ = (q+1)*(q+2)/2;
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(dimP*dimQ*nE,1,mxREAL);
  
  /* Get pointer to Matrix */
  double * M = mxGetPr(plhs[0]);
  
  /* Define a function pointer */
  void (*compFun)(double*,double,double,double,double,double,double) ;
  
  /* Depending on p and q choose the right function */
  if(p==0) {
      if(q==1)
          compFun = &compTSigmaHatVP0Q1;
      else if(q==2)
          compFun = &compTSigmaHatVP0Q2;
      else if(q==3)
          compFun = &compTSigmaHatVP0Q3;
      else if(q==4)
          compFun = &compTSigmaHatVP0Q4;
      else if(q==5)
          compFun = &compTSigmaHatVP0Q5;
  }                    
  else if(p==1) {
      if(q==1)
          compFun = &compTSigmaHatVP1Q1;
      else if(q==2)
          compFun = &compTSigmaHatVP1Q2;
      else if(q==3)
          compFun = &compTSigmaHatVP1Q3;
      else if(q==4)
          compFun = &compTSigmaHatVP1Q4;
      else if(q==5)
          compFun = &compTSigmaHatVP1Q5;
  }
  else if(p==2) {
      if(q==1)
          compFun = &compTSigmaHatVP2Q1;
      else if(q==2)
          compFun = &compTSigmaHatVP2Q2;
      else if(q==3)
          compFun = &compTSigmaHatVP2Q3;
      else if(q==4)
          compFun = &compTSigmaHatVP2Q4;
      else if(q==5)
          compFun = &compTSigmaHatVP2Q5;
  }
  else if(p==3) {
      if(q==1)
          compFun = &compTSigmaHatVP3Q1;
      else if(q==2)
          compFun = &compTSigmaHatVP3Q2;
      else if(q==3)
          compFun = &compTSigmaHatVP3Q3;
      else if(q==4)
          compFun = &compTSigmaHatVP3Q4;
      else if(q==5)
          compFun = &compTSigmaHatVP3Q5;
  }
      
  
  //mexPrintf("p = %d, q = %d, dimP = %d, dimQ = %d\n",p,q,dimP,dimQ);
  
  /* for loop */
  int i = 0;
 
  #pragma omp parallel for
  for(i=0;i<nE; ++i) {
      //compMtmpP0Q2(M+i*dimP*dimQ,a[i],b[i],c[i],d[i]);
      compFun(M+i*dimP*dimQ,e1[i],e2[i],e3[i],o1[i],o2[i],o3[i]);
  }
 
}