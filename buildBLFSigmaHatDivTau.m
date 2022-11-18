function [B,I,J] = buildBLFSigmaHatDivTau(coordinates,edges,element2edges,orientation,p,q)
    dimP = 3*(p+1);
    dimQ = (q+1)*(q+2);

    ce1 = coordinates(edges(:,1),:);
    ce2 = coordinates(edges(:,2),:);
    edgL = sqrt(sum((ce1-ce2).^2,2));


    edgLe2e = edgL(element2edges);

    nE = size(element2edges,1);
    nEd = size(edges,1);
    
    b21 = coordinates(edges(element2edges(:,1),1),:)-coordinates(edges(element2edges(:,1),2),:);
    b31 = coordinates(edges(element2edges(:,3),1),:)-coordinates(edges(element2edges(:,3),2),:);
    area2 = (b21(:,1).*b31(:,2)-b21(:,2).*b31(:,1));

% %%
    if(p==0)
        orientP = orientation;
    elseif(p==1)
        orientP = [orientation ones(size(orientation))];
    elseif(p==2)
        orientP = [orientation ones(size(orientation)) orientation];
    end
%% get the L2 matrix
% blfMat = getSigmaHatVmat(p,q);
% 
% %% go over all elements
    B = zeros(dimP*dimQ*nE,1);
% 
    I = zeros(dimP*dimQ*nE,1);
    J = zeros(dimP*dimQ*nE,1);
% 
for jdx=1:nE
     
    tmpO = (1/area2(jdx))*compSigmaHatDivTau(edgLe2e(jdx,1),edgLe2e(jdx,2),edgLe2e(jdx,3),orientation(jdx,1),orientation(jdx,2),orientation(jdx,3),b31(jdx,1),b31(jdx,2),b21(jdx,1),b21(jdx,2),p,q)';
    B( (jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = tmpO;
    
    Itmp = repmat((1:dimQ)',1,dimP) + (jdx-1)*dimQ;
    I( (jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = Itmp(:);
    
    if(p==0)
        Jtmp = repmat(element2edges(jdx,:),dimQ,1);
    elseif(p==1)
        Jtmp = repmat( [element2edges(jdx,:) element2edges(jdx,:)+nEd],dimQ,1);
    elseif(p==2)
        Jtmp = repmat( [element2edges(jdx,:) element2edges(jdx,:)+nEd element2edges(jdx,:)+2*nEd],dimQ,1);
    end
    J( (jdx-1)*dimP*dimQ+1:jdx*dimP*dimQ) = Jtmp(:);
end

if(nargout==1)
    B = sparse(I,J,B,dimQ*nE,(p+1)*nEd);
end

end
