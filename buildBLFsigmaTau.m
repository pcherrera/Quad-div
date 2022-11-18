function [B I J] = buildBLFsigmaTau(coordinates,elements,p,q,A)

dimP = (p+1)*(p+2);
dimQ = (q+1)*(q+2);

d21 = coordinates(elements(:,2),:)-coordinates(elements(:,1),:);
d31 = coordinates(elements(:,3),:)-coordinates(elements(:,1),:);

area2 = (d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1));

nE = size(elements,1);

if(isempty(A))
    A = [ones(nE,1) zeros(nE,1) zeros(nE,1) ones(nE,1)];
end

% %% get the L2 matrix
% blfMat = getFblfSIGMATAU(p,q);
% 
% %% go over all elements
% B = zeros(dimP*dimQ*nE,1);
% for jdx=1:nE
%     
%     tmp = area2(jdx)*blfMat(A(jdx,1),A(jdx,2),A(jdx,3),A(jdx,4));
%     
%     B( (jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = tmp(:);
% end

% mex function
B = buildSigmaTauMEX(area2.*A(:,1),area2.*A(:,2),area2.*A(:,3),area2.*A(:,4),p,q);

Itmp = repmat((1:dimQ)',1,dimP);
Itmp = Itmp(:);
Jtmp = repmat(1:dimP,dimQ,1);
Jtmp = Jtmp(:);

I = repmat(Itmp,nE,1)+reshape(repmat(  (1:dimQ:dimQ*nE)-1, dimQ*dimP,1),dimQ*dimP*nE,1);
J = repmat(Jtmp,nE,1)+reshape(repmat(  (1:dimP:dimP*nE)-1, dimQ*dimP,1),dimQ*dimP*nE,1);

if(nargout==1)
    B = sparse(I,J,B,dimQ*nE,dimP*nE);
end