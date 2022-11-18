function [B,I,J] = buildBMatLap(coordinates,elements,edges,element2edges,...
    pU1,pU2,pU1hat,qU1hat,pU2hat,qU2hat,qV1,qV2)
%% builds bilinear form with 
% Field Variables: u1, u2. 
% Trace Variables: u1.n, divu1, u2.n, divu2. 
% Test Variables: v1, v2. 
% 
% Bilinear form (u1,grad div v2 + v1)+(u2,grad div v1 - v2) ...
% + < div u1,v2.n> - <u1.n,div v2> + < div u2, v1.n> - <u2.n,div v1>

nE = size(elements,1);
nEd = size(edges,1);
nC = size(coordinates,1);

tmp = edges(element2edges,:);
idx1 = tmp(1:nE,1) == elements(:,1);
idx2 = tmp(nE+1:2*nE,1) == elements(:,2);
idx3 = tmp(2*nE+1:3*nE,1) == elements(:,3);
orientation = -ones(nE,3);
orientation([idx1 idx2 idx3]) = -1;

%% 
B = [];
I = [];
J = [];

dimU1 = nE*(pU1+1)*(pU1+2);
dimU2 = nE*(pU2+1)*(pU2+2);
dimqU1hat = nEd*(qU1hat+1);
dimpU1hat = nC + nEd*(pU1hat-1);
dimqU2hat = nEd*(qU2hat+1);
dimpU2hat = nC + nEd*(pU2hat-1);

dimV1 = nE*(qV1+1)*(qV1+2);
dimV2 = nE*(qV2+1)*(qV2+2);

dimU = dimU1+ dimpU1hat + dimqU1hat + dimU2 + dimpU2hat + dimqU2hat;
dimV = dimV1 + dimV2;

%% (u1,v1)
[Btmp,Itmp,Jtmp] = buildBLFsigmaTau(coordinates,elements,pU1,qV1,[]);

I = [I;Itmp];
J = [J;Jtmp];
B = [B;Btmp];

%% (u1,grad div(v2))
[Btmp,Itmp,Jtmp] = buildBLFsigmaGdivTau(coordinates,elements,pU1,qV2);

I = [I;Itmp+dimV1];
J = [J;Jtmp];
B = [B;Btmp];

%% (u2,grad div(v1))
[Btmp,Itmp,Jtmp] = buildBLFsigmaGdivTau(coordinates,elements,pU2,qV1);

I = [I;Itmp];
J = [J;Jtmp+dimU1];
B = [B;Btmp];

%% -(u2,v2)
[Btmp,Itmp,Jtmp] = buildBLFsigmaTau(coordinates,elements,pU2,qV2,[]);

I = [I;Itmp+dimV1];
J = [J;Jtmp+dimU1];
B = [B;-Btmp];

%% (div hat u1,[v2.n])
[Btmp,Itmp,Jtmp] = buildBLFuHatTauN(coordinates,elements,edges,element2edges,orientation,pU1hat,qV2);

I = [I;Itmp+dimV1];
J = [J;Jtmp+dimU1+dimU2+dimqU1hat];
B = [B;Btmp];

%% (u1.n,[div v2])
[Btmp,Itmp,Jtmp] = buildBLFSigmaHatDivTau(coordinates,edges,element2edges,orientation,qU1hat,qV2);

I = [I;Itmp+dimV1];
J = [J;Jtmp+dimU1+dimU2];
B = [B;-Btmp];

%% (div hat u2,v1.n)
[Btmp,Itmp,Jtmp] = buildBLFuHatTauN(coordinates,elements,edges,element2edges,orientation,pU2hat,qV1);

I = [I;Itmp];
J = [J;Jtmp+dimU1+dimU2+dimqU1hat+dimpU1hat+dimqU1hat];
B = [B;Btmp];

%% (u2.n,[div v1])
[Btmp,Itmp,Jtmp] = buildBLFSigmaHatDivTau(coordinates,edges,element2edges,orientation,qU2hat,qV1);

I = [I;Itmp];
J = [J;Jtmp+dimU1+dimU2+dimqU1hat+dimpU1hat];
B = [B;-Btmp];

%% OUTPUT
if(nargout==1)
    B = sparse(I,J,B,dimV,dimU);
end