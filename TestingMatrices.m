clear all
close all
%% Mesh 
%addpath Data/
%addpath Data1/
%load Meshes/quadMesh.mat
h=1;
%coordinates = h*coordinates;
coordinates = [0,0;h,0;0,h;h,h];
elements = [2,3,1;3,2,4];
edges = [1 2;1 4;2 4;3 2;4 3];
boundary = [1 2;2 3;3 4;1 4];
element2edges = [1 2 3;3 4 5];

%% Example 
U = @(x) [x(:,1).^4,x(:,2).^4];
U_1 = @(x) x(:,1).^4;
U_2 = @(x) x(:,2).^4;
divU = @(x) 4*x(:,1).^3 + 4*x(:,2).^3;
gdivU = @(x) [12*x(:,1).^2,12*x(:,2).^2];
gdivU_1 = @(x) 12*x(:,1).^2;
gdivU_2 = @(x) 12*x(:,2).^2;
f = @(x) U(x)+0*[24+0*x(:,1),24+0*x(:,2)];
u0 = @(x) U(x);

%% ansatz space
p = 0;
pU1 = p;
pU2 = p;

% trace Variables
pU1hat = 1;
qU1hat = 0;
pU2hat = 1;
qU2hat = 0;

% test space
qV1 = 3;
qV2 = 3;

%% Set adaptivity parameter \theta
% 0<theta<1 is adaptive refinement using 
% theta = 1 is uniform mesh refinement
theta = 1;

% store some values (number of elements, vertices, etc.)
nE = size(elements,1);
nC = size(coordinates,1);
nB = size(boundary,1);
nEd = size(edges,1);
%nBEd = size(boundary2edges,1);

tmp = edges(element2edges,:);
idx1 = tmp(1:nE,1) == elements(:,1);
idx2 = tmp(nE+1:2*nE,1) == elements(:,2);
idx3 = tmp(2*nE+1:3*nE,1) == elements(:,3);
orientation = -ones(nE,3);
orientation([idx1 idx2 idx3]) = 1;

fprintf('Step %d, nE = %d, nC = %d, nB = %d, nEd = %d\n',1,nE,nC,nB,nEd);

    
% dimension of discretization spaces
dimU1 = nE*(pU1+1)*(pU1+2);
dimU2 = nE*(pU2+1)*(pU2+2);
dimU = dimU1+dimU2;
    
dimU1qhat = nEd*(qU1hat+1);
dimU1phat = nC+nEd*(pU1hat-1);
dimU2qhat = nEd*(qU2hat+1);
dimU2phat = nC+nEd*(pU2hat-1);
dimUhat = dimU1phat+dimU1qhat+dimU2phat+dimU2qhat;
    
dimV1 = nE*(qV1+1)*(qV1+2);
dimV2 = nE*(qV2+1)*(qV2+2);
dimV = dimV1+dimV2;
    
% degrees of freedom (approximately)
dofDPG = dimU + dimUhat;

%fprintf('j = %d\n    Building B ...\n',j);
B = buildBMatLap(coordinates,elements,edges,element2edges,...
pU1,pU2,pU1hat,qU1hat,pU2hat,qU2hat,qV1,qV2);
D= full(B);
%fprintf('    Building V ...\n');
Vinv = buildVinvMatLap(coordinates,elements,qV1,qV2);
E = full(Vinv);    
%% *** right-hand side
bVec = buildRHSLap(coordinates,elements,qV1,qV2,f);
    
    %% boundary dof's 
%      jdx =[];
%      if(pU1hat>=1)
%         %B(:,dimU+dimSigma+nC(j)-nB(j)+1:dimU+dimSigma+nC(j)) = [];
%         jdx = [jdx;(dimU+nC-nB+1:dimU+nC)'];
%      end
%      if(pU1hat>=2)
%         %B(:,dimU+dimSigma+nC(j)-nB(j)+1:dimU+dimSigma+nC(j)) = [];
%         jdx = [jdx;(dimU+nC+nEd-nB+1:dimU+nC+nEd)'];
%      end
%     
%      kdx = [];
%      if (qU1hat>=0)
%          kdx = [kdx;(dimU+dimU1phat+nEd-nBEd+1:dimU+dimU1phat+nEd)'];
%      end
%      if (qU1hat>=1)
%          kdx = [kdx;(dimU+dimU1phat+2*nEd-nBEd+1:dimU+dimU1phat+2*nEd)']; 
%      end
%      if (qU1hat>=2) 
%          kdx = [kdx; (dimU+dimU1phat+3*nEd-nBEd+1:dimU+dimU1phat+3*nEd)'];
%      end        

    
    
%     %% Define system matrix: S = B'V^{-1} B where V is Riesz-matrix
%     b = B'*(Vinv*bVec);
%     S = B'*(Vinv*B);
%     
%     %*** determine indices corresponding to the dof's
%     dof = size(S,1);
%     freeonodes = setdiff( (1:dof)',[jdx;kdx]);
% 
%     %*** initialise solution vector x
%     %fprintf('    Solve ...\n');
%     x = zeros(dof,1);
%     
%     %*** Dirichlet boundary values
%     x(jdx) = computeL2projBOUC1(coordinates(nC-nB+1:end,:),boundary-nC+nB,divU,pU1hat);
%     x(kdx) = computeL2projNeumannPp(coordinates(nC-nB+1:end,:),boundary-nC+nB,u0,qU1hat);
%     %x(jdx) = 0;
%     %x(kdx) = 0;
%     b = b-S*x;
%     
%     x(freeonodes) = S(freeonodes,freeonodes)\b(freeonodes);
% 
%     nEle = nE;
% 
%     %*** DPG error estimator
%     err = sqrt( (bVec-B*x)'*Vinv*(bVec-B*x) );
% 
%     %*** L2 error of u1
%     U1H = [x(1:2:dimU1) x(2:2:dimU1)];
%     tmpErr1 = computeErrL2(coordinates,elements,U1H,U,pU1);
%     errU1 = sqrt( sum( sum(tmpErr1,2) ) );
%     
%         
%     %*** L2 error of u2
%     U2H = [x(1+dimU1:2:dimU1+dimU2) x(2+dimU1:2:dimU1+dimU2)];
%     tmpErr2 = computeErrL2(coordinates,elements,U2H,gdivU,pU2);
%     errU2 = sqrt( sum( sum(tmpErr2,2) ) );
%     
%     % *** L2 error of proj U1
%     %U1P = computeL2projToPp_Vec(coordinates,elements,U_1,U_2,pU1);
%     %tmpL2NormProjU1 = computeErrL2(coordinates,elements,U1P,U,pU1);
%     %errU1P(j) = sqrt( sum( sum(tmpL2NormProjU1,2) ) );
%     
%     % *** L2 error of proj U2
%     %U2P = computeL2projToPp_Vec(coordinates,elements,gdivU_1,gdivU_2,pU2);
%     %tmpL2NormProjU2 = computeErrL2(coordinates,elements,U2P,gdivU,pU2);
%     %errU2P(j) = sqrt( sum( sum(tmpL2NormProjU2,2) ) );
%         
%     % Relative Error 
%     %tmpL2NormU1 = computeErrL2(coordinates,elements,0,U,pU1);
%     %normU1(j) = sqrt(sum(sum(tmpL2NormU1,2)));
%     %rel_errU1(j) = (errU1(j)/normU1(j)); %Relative Error for U1 dpg solution
%     %rel_errU1P(j) = (errU1P(j)/normU1(j)); %Relative Error for U1 projection
%     
%     %*** Functions
%     xU1 = [x(1:2:dimU1) x(2:2:dimU1)];   
%     xU2 = [x(dimU1+1:2:dimU1+dimU2) x(dimU1+2:2:dimU1+dimU2)];   
%     
% 
% 
%     %*** error quantity (= DPG-error estimator)
%     est = (bVec-B*x).*(Vinv*(bVec-B*x));
% 
% 
%     % H(Gdiv) parts for V1
%     est1 = est(1:dimV1);
%     est1 = sum( reshape(est1,(qV1+1)*(qV1+2),nE),1)';
%     % H(Gdiv) parts for V2
%     est2 = est(dimV1+1:dimV1+dimV2);
%     est2 = sum( reshape(est2,(qV2+1)*(qV2+2),nE),1)';
%    
%     % overall volume-elementwise estimator:
%     estVol = est1 +est2 ;


