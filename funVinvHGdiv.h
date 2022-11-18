#include "mex.h"
#include<math.h>
#include "invSymMat.h"
// #define phi0 2.449489742783178
// #define phi1 -3.162277660168380
// #define phi2 -0.935414346693485
inline void compVinvHDGdivQ0(double * M, double a, double b, double c)
inline void compVinvHDGdivQ1(double * M, double a, double b, double c)
// Automated generated file
#include "funSigmaGdivTau.h"
void compTsigmaGdivTauQ0(double * M, double a, double b, double c, double d) {
M[0]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/2 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/2
M[1]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/2 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/2
M[2]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/2 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/2
M[3]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/2 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/2
}
void compTsigmaGdivTauQ1(double * M, double a, double b, double c, double d) {
M[0]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/2 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/2
M[1]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/2 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/2
M[2]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[3]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[4]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[5]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[6]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/2 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/2
M[7]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/2 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/2
M[8]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[9]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[10]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[11]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[12]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/2 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/2
M[13]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/2 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/2
M[14]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[15]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[16]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[17]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[18]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/2 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/2
M[19]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/2 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/2
M[20]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[21]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[22]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[23]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[24]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/2 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/2
M[25]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/2 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/2
M[26]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[27]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[28]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[29]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[30]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/2 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/2
M[31]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/2 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/2
M[32]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[33]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
M[34]=(a*(a*c11 + c*c12) + c*(a*c12 + c*c22))^2/12 + (a*(b*c11 + c12*d) + c*(b*c12 + c22*d))^2/12
M[35]=(b*(a*c11 + c*c12) + d*(a*c12 + c*c22))^2/12 + (b*(b*c11 + c12*d) + d*(b*c12 + c22*d))^2/12
}
