function [B I J] = buildBLFsigmaHatV(coordinates,edges,element2edges,orientation,p,q)

dimQ = (q+1)*(q+2)/2;
dimP = (p+1)*3;

ce1 = coordinates(edges(:,1),:);
ce2 = coordinates(edges(:,2),:);
edgL = sqrt(sum( (ce1-ce2).^2,2));



edgLe2e = edgL(element2edges);

nE = size(element2edges,1);
nEd = size(edges,1);

% %%
% if(p==0)
%     orientP = orientation;
% elseif(p==1)
%     orientP = [orientation ones(size(orientation))];
% elseif(p==2)
%     orientP = [orientation ones(size(orientation)) orientation];
% end
% %% get the L2 matrix
% blfMat = getSigmaHatVmat(p,q);
% 
% %% go over all elements
% B = zeros(dimP*dimQ*nE,1);
% 
% I = zeros(dimP*dimQ*nE,1);
% J = zeros(dimP*dimQ*nE,1);
% 
% for jdx=1:nE
%     
%     tmpO = repmat(edgLe2e(jdx,:),dimQ,p+1).*repmat(orientP(jdx,:),dimQ,1).*blfMat;
%     B( (jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = tmpO(:);
%     
%     Itmp = repmat((1:dimQ)',1,dimP) + (jdx-1)*dimQ;
%     I( (jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = Itmp(:);
%     
%     if(p==0)
%         Jtmp = repmat(element2edges(jdx,:),dimQ,1);
%     elseif(p==1)
%         Jtmp = repmat( [element2edges(jdx,:) element2edges(jdx,:)+nEd],dimQ,1);
%     elseif(p==2)
%         Jtmp = repmat( [element2edges(jdx,:) element2edges(jdx,:)+nEd element2edges(jdx,:)+2*nEd],dimQ,1);
%     end
%     J( (jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = Jtmp(:);
% end

Itmp = repmat((1:dimQ)',1,dimP);
Itmp = Itmp(:);
I = repmat(Itmp,nE,1)+reshape(repmat(  (1:dimQ:dimQ*nE)-1, dimQ*dimP,1),dimQ*dimP*nE,1);
clear Itmp

Jtmp = [repmat(element2edges(:,1)',dimQ,1);repmat(element2edges(:,2)',dimQ,1);repmat(element2edges(:,3)',dimQ,1)];
if(p==1)
    Jtmp = [Jtmp;Jtmp+nEd];
elseif(p==2)
    Jtmp = [Jtmp;Jtmp+nEd;Jtmp+2*nEd];
elseif(p==3)
    Jtmp = [Jtmp;Jtmp+nEd;Jtmp+2*nEd;Jtmp+3*nEd];
end
J = Jtmp(:);
clear Jtmp

B = buildSigmaHatVMEX(edgLe2e(:,1),edgLe2e(:,2),edgLe2e(:,3),orientation(:,1),orientation(:,2),orientation(:,3),p,q);


if(nargout==1)
    B = sparse(I,J,B,dimQ*nE,(p+1)*nEd);
end