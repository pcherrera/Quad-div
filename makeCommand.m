%% Make commands:

%% Bilinear forms
%mex buildSigmaGradVMEX.c CFLAGS="\$CFLAGS -std=c99 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims
mex buildSigmaTauMEX.c CFLAGS="\$CFLAGS -std=c99 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims
mex buildSigmaGdivTauMEX.c CFLAGS="\$CFLAGS -std=c99 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims
%mex buildUDivTauMEX.c CFLAGS="\$CFLAGS -std=c99 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims
%mex buildSigmaHatVMEX.c CFLAGS="\$CFLAGS -std=c99 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims
%mex buildUHatTauNMEX.c CFLAGS="\$CFLAGS -std=c99 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims
%mex buildUVMEX.c CFLAGS="\$CFLAGS -std=c99 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -largeArrayDims

%% (inverse) Riesz matrices
%mex buildVinvH1MEX.c -largeArrayDims -lmwlapack CFLAGS="\$CFLAGS -std=c99 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
%mex buildVinvHDIVMEX.c -largeArrayDims -lmwlapack CFLAGS="\$CFLAGS -std=c99 -Wall -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"