function [B I J] = buildBLFuHatTauN(coordinates,elements,edges,element2edges,orientation,p,q)

dimQ = (q+1)*(q+2);
dimP = p*3;

nE = size(elements,1);
nEd = size(edges,1);
nC = size(coordinates,1);

% edge coefficients
d21 = coordinates(elements(:,2),:)-coordinates(elements(:,1),:);
d32 = coordinates(elements(:,3),:)-coordinates(elements(:,2),:);
d13 = coordinates(elements(:,1),:)-coordinates(elements(:,3),:);
a1 = d21(:,1); a2 = d21(:,2);
b1 = d32(:,1); b2 = d32(:,2);
c1 = d13(:,1); c2 = d13(:,2);
clear d21 d32 d13

% % orientation
% if(p==1)
%     orientP = ones(size(orientation));
% elseif(p==2)
%     orientP = [ones(size(orientation)) ones(size(orientation))];
% elseif(p==3)
%     orientP = [ones(size(orientation)) ones(size(orientation)) orientation];
% end
% 
% %% get the L2 matrix
% uHatTauF = getFblfUHATTAUN(p,q);
% 
% %% go over all elements
% B = zeros(dimP*dimQ*nE,1);
% 
% I = zeros(dimP*dimQ*nE,1);
% J = zeros(dimP*dimQ*nE,1);
% 
% for jdx=1:nE
%     
%     tmpO = repmat(orientP(jdx,:),dimQ,1).*uHatTauF(a1(jdx),a2(jdx),b1(jdx),b2(jdx),c1(jdx),c2(jdx));
%     B( (jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = tmpO(:);
%     
%     Itmp = repmat((1:dimQ)',1,dimP) + (jdx-1)*dimQ;
%     I( (jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = Itmp(:);
%     
%     if(p==1)
%         Jtmp = repmat(elements(jdx,:),dimQ,1);
%     elseif(p==2) 
%         Jtmp = [repmat(elements(jdx,:),dimQ,1) repmat(element2edges(jdx,:)+nC,dimQ,1)];
%     elseif(p==3)
%         Jtmp = [repmat(elements(jdx,:),dimQ,1) repmat(element2edges(jdx,:)+nC,dimQ,1) repmat(element2edges(jdx,:)+nC+nEd,dimQ,1)];
%     end
%     J( (jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = Jtmp(:);
% end

Itmp = repmat((1:dimQ)',1,dimP);
Itmp = Itmp(:);
I = repmat(Itmp,nE,1)+reshape(repmat(  (1:dimQ:dimQ*nE)-1, dimQ*dimP,1),dimQ*dimP*nE,1);
clear Itmp

Jtmp = [repmat(elements(:,1)',dimQ,1);repmat(elements(:,2)',dimQ,1);repmat(elements(:,3)',dimQ,1)];
if(p==2)
    Jtmp = [Jtmp; repmat(element2edges(:,1)'+nC,dimQ,1);repmat(element2edges(:,2)'+nC,dimQ,1);repmat(element2edges(:,3)'+nC,dimQ,1)];
elseif(p==3)
    Jtmp = [Jtmp; repmat(element2edges(:,1)'+nC,dimQ,1);repmat(element2edges(:,2)'+nC,dimQ,1);repmat(element2edges(:,3)'+nC,dimQ,1); ...
        repmat(element2edges(:,1)'+nC+nEd,dimQ,1);repmat(element2edges(:,2)'+nC+nEd,dimQ,1);repmat(element2edges(:,3)'+nC+nEd,dimQ,1)];
elseif(p==4)
    Jtmp = [Jtmp; repmat(element2edges(:,1)'+nC,dimQ,1);repmat(element2edges(:,2)'+nC,dimQ,1);repmat(element2edges(:,3)'+nC,dimQ,1); ...
        repmat(element2edges(:,1)'+nC+nEd,dimQ,1);repmat(element2edges(:,2)'+nC+nEd,dimQ,1);repmat(element2edges(:,3)'+nC+nEd,dimQ,1); ...
        repmat(element2edges(:,1)'+nC+2*nEd,dimQ,1);repmat(element2edges(:,2)'+nC+2*nEd,dimQ,1);repmat(element2edges(:,3)'+nC+2*nEd,dimQ,1)];
end
J = Jtmp(:);
clear Jtmp

B = buildUHatTauNMEX(a1,a2,b1,b2,c1,c2,orientation(:,1),orientation(:,2),orientation(:,3),p,q);


if(nargout==1)
    B = sparse(I,J,B,dimQ*nE,nC+(p-1)*nEd);
end