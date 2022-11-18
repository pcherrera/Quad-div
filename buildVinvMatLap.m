function [V,I,J] = buildVinvMatLap(coordinates,elements,qV1,qV2)

nE = size(elements,1);

%% 
V = [];
I = [];
J = [];

dimV1 = nE*(qV1+1)*(qV1+2);
dimV2 = nE*(qV2+1)*(qV2+2);
dimV = dimV1+dimV2;

%% v1 inner products (broken H(Gdiv))
[Vtmp,Itmp,Jtmp] = buildVinvMatHGDIVexact(coordinates,elements,qV1,[],[]);

I = [I;Itmp];
J = [J;Jtmp];
V = [V;Vtmp];

%% v2 inner products (broken H(Gdiv))
[Vtmp,Itmp,Jtmp] = buildVinvMatHGDIVexact(coordinates,elements,qV2,[],[]);

I = [I;Itmp+dimV1];
J = [J;Jtmp+dimV1];
V = [V;Vtmp];

%% OUTPUT
%if(nargout==1)
    V = sparse(I,J,V,dimV,dimV);
%end