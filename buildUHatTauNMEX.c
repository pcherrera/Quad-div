#include "mex.h"
#include<math.h>

#include "mexHelpers/funUhatTauN.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ) {
  
 
  double * a1 = mxGetPr(prhs[0]);
  double * a2 = mxGetPr(prhs[1]);
  double * b1 = mxGetPr(prhs[2]);
  double * b2 = mxGetPr(prhs[3]);
  double * c1 = mxGetPr(prhs[4]);
  double * c2 = mxGetPr(prhs[5]);
  double * o1 = mxGetPr(prhs[6]);
  double * o2 = mxGetPr(prhs[7]);
  double * o3 = mxGetPr(prhs[8]);
  
  int p = (int) mxGetScalar(prhs[9]);
  int q = (int) mxGetScalar(prhs[10]);

  /* number of elements */
  int nE = mxGetM(prhs[0]);
  
  int dimP = p*3;
  int dimQ = (q+1)*(q+2);
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(dimP*dimQ*nE,1,mxREAL);
  
  /* Get pointer to Matrix */
  double * M = mxGetPr(plhs[0]);
  
  /* Define a function pointer */
  void (*compFun)(double*,double,double,double,double,double,double,double,double,double) ;
  
  /* Depending on p and q choose the right function */
  if(p==1) {
      if(q==1)
          compFun = &compTUHatTauNP1Q1;
      else if(q==2)
          compFun = &compTUHatTauNP1Q2;
      else if(q==3)
          compFun = &compTUHatTauNP1Q3;
      else if(q==4)
          compFun = &compTUHatTauNP1Q4;
      else if(q==5)
          compFun = &compTUHatTauNP1Q5;
  }                    
  else if(p==2) {
      if(q==1)
          compFun = &compTUHatTauNP2Q1;
      else if(q==2)
          compFun = &compTUHatTauNP2Q2;
      else if(q==3)
          compFun = &compTUHatTauNP2Q3;
      else if(q==4)
          compFun = &compTUHatTauNP2Q4;
      else if(q==5)
          compFun = &compTUHatTauNP2Q5;
  }
  else if(p==3) {
      if(q==1)
          compFun = &compTUHatTauNP3Q1;
      else if(q==2)
          compFun = &compTUHatTauNP3Q2;
      else if(q==3)
          compFun = &compTUHatTauNP3Q3;
      else if(q==4)
          compFun = &compTUHatTauNP3Q4;
      else if(q==5)
          compFun = &compTUHatTauNP3Q5;
  }
  else if(p==4) {
      if(q==1)
          compFun = &compTUHatTauNP4Q1;
      else if(q==2)
          compFun = &compTUHatTauNP4Q2;
      else if(q==3)
          compFun = &compTUHatTauNP4Q3;
      else if(q==4)
          compFun = &compTUHatTauNP4Q4;
      else if(q==5)
          compFun = &compTUHatTauNP4Q5;
  }
      
  
  //mexPrintf("p = %d, q = %d, dimP = %d, dimQ = %d\n",p,q,dimP,dimQ);
  
  /* for loop */
  int i = 0;
 
  #pragma omp parallel for
  for(i=0;i<nE; ++i) {
      //compMtmpP0Q2(M+i*dimP*dimQ,a[i],b[i],c[i],d[i]);
      compFun(M+i*dimP*dimQ,a1[i],a2[i],b1[i],b2[i],c1[i],c2[i],o1[i],o2[i],o3[i]);
  }
 
}