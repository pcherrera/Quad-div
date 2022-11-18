function [Vinv,I,J] = buildVinvMatHGDIVexact(coordinates,elements,p,A,C)
%% buildVinvMatHGDIVexact computes
%   the inverse of the Galerkinmatrix of the Riesz operator associated to
%   the inner product
%
%       (v,w)_V = (C(x)v,w) + (A(x)grad(div v), grad(div w))
%
%   on each elementwise. The resulting matrix is a block-diagonal matrix,
%   where the first block corresponds to the first element given by the
%   matrices (coordinates,elements), the second block to the second, etc.
%   The output is either the matrix itself or the triple (Vinv,I,J), such
%   that sparse(I,J,Vinv) is the assembled matrix.
%
%   p \in\{1,2,3,4} is the polynomial degree, the basis then has dimension
%       (p+1)*(p+2)
%   here we use a hierarchical basis of Gauss-Lobatto type
%
%   If C = [], then it is assumed that C = Identity matriz, if b is a vector 
%   with length(C) == nE (number of elements) then C is thought of being
%   elementwise constant
%   
%   Also if A = [], then A = Identity matrix
%   Or A(j,:) = [a11 a12 a21 a22] elementwise (note that a12==a21!!! is too
%   much ;-)
%
%      

%% Some useful def.
    nE = size(elements,1);
    if(p>5 || p<+0)
        error('Polynomial degree %d not allowed!',p);
    end

    dimPp = (p+1)*(p+2);

    d21 = coordinates(elements(:,2),:)-coordinates(elements(:,1),:);
    d31 = coordinates(elements(:,3),:)-coordinates(elements(:,1),:);
    area2 = (d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1));

    if(isempty(C))
        C = [ones(nE,1) zeros(nE,1) zeros(nE,1) ones(nE,1)];
    end

    if(isempty(A))
        A = [ones(nE,1) zeros(nE,1) zeros(nE,1) ones(nE,1)];
    end

%% Get Matrix and go over all elements
    Minv = zeros(dimPp^2*nE,1);
    for jdx=1:nE
        d211 = d21(jdx,1);
        d212 = d21(jdx,2);
        d311 = d31(jdx,1);
        d312 = d31(jdx,2);
        if( abs(A(jdx,2)-A(jdx,3))>=1e-14)
            error('Input matrix A not symmetric?')
        end
     
        %Mtmp = compVinvHGdiv(A(jdx,1),A(jdx,2),A(jdx,4),C(jdx,1),C(jdx,2),C(jdx,4),d211,d212,d311,d312,p);
        Mtmp = compVinvHGdiv(d211,d212,d311,d312,p);  
        Minv(1+dimPp^2*(jdx-1):dimPp^2*jdx) = Mtmp(:);
    end
%% indices
    Itmp = repmat((1:dimPp)',1,dimPp);
    Itmp = Itmp(:);
    Jtmp = repmat(1:dimPp,dimPp,1);
    Jtmp = Jtmp(:);

    I = repmat(Itmp,nE,1)+reshape(repmat(  (1:dimPp:dimPp*nE)-1, dimPp^2,1),dimPp^2*nE,1);
    J = repmat(Jtmp,nE,1)+reshape(repmat(  (1:dimPp:dimPp*nE)-1, dimPp^2,1),dimPp^2*nE,1);

    if(nargout==1)
        Vinv = sparse(I,J,Minv,dimPp*nE,dimPp*nE);
    elseif (nargout==3)
        Vinv = Minv;
    end
end
