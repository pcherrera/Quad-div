function [B,I,J] = buildBLFsigmaGdivTau(coordinates,elements,p,q)

dimP = (p+2)*(p+1);
dimQ = (q+1)*(q+2);

d21 = coordinates(elements(:,2),:)-coordinates(elements(:,1),:);
d31 = coordinates(elements(:,3),:)-coordinates(elements(:,1),:);

area2 = (d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1));

nE = size(elements,1);

%if(isempty(A))
%    A = [ones(nE,1) zeros(nE,1) zeros(nE,1) ones(nE,1)];
%end



% %% get the L2 matrix
%blfMat = compSigmaGdivTau(d21(:,1),d21(:,2),d31(:,1),d31(:,2),p,q);
 
%% go over all elements
B = zeros(dimP*dimQ*nE,1);
    for jdx=1:nE 
        tmp = area2(jdx)*compSigmaGdivTau(d21(jdx,1),d21(jdx,2),d31(jdx,1),d31(jdx,2),p,q);
        B((jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = tmp;
    end

% mex function
%B = buildSigmaGdivTauMEX(area2.*A(:,1),area2.*A(:,2),area2.*A(:,3),area2.*A(:,4),p,q);

Itmp = repmat((1:dimQ)',1,dimP);
Itmp = Itmp(:);
Jtmp = repmat(1:dimP,dimQ,1);
Jtmp = Jtmp(:);

I = repmat(Itmp,nE,1)+reshape(repmat(  (1:dimQ:dimQ*nE)-1, dimQ*dimP,1),dimQ*dimP*nE,1);
J = repmat(Jtmp,nE,1)+reshape(repmat(  (1:dimP:dimP*nE)-1, dimQ*dimP,1),dimQ*dimP*nE,1);

if(nargout==1)
    B = sparse(I,J,B,dimQ*nE,dimP*nE);
end