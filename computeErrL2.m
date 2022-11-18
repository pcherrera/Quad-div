function err = computeErrL2(coordinates,elements,uh,u,p)
% uh is elementwise polynomial function
%
% We use the same basis (of course) as in the DPG matrices

%*** determine dimensions of uh
dim = size(uh,2);

dimU = (p+1)*(p+2)/2;

%***
err = zeros(size(elements,1),dim);

%*** quadrature rule on reference triangle Tref = conv{(0,0),(1,0),(0,1)}
tmp_pos = [6-sqrt(15) ; 9+2*sqrt(15) ; 6+sqrt(15) ; 9-2*sqrt(15) ; 7]/21;
quad_vertices = tmp_pos([1 1 ; 2 1 ; 1 2 ; 3 4 ; 3 3 ; 4 3 ; 5 5]);
tmp_wts = [155-sqrt(15) 155+sqrt(15) 270]/2400;
quad_weights = tmp_wts([1 1 1 2 2 2 3]);

% xR = [-0.333333333333333 -0.059715871789770 -0.059715871789770 -0.880568256420460 -0.797426985353088 -0.797426985353088 0.594853970706174];
% yR = [-0.333333333333333 -0.059715871789770 -0.880568256420460 -0.059715871789770 -0.797426985353088 0.594853970706174 -0.797426985353088];
% wR = [0.450000000000000 0.264788305577012 0.264788305577012 0.264788305577012 0.251878361089654 0.251878361089654 0.251878361089654]';
% w = 0.25*wR;
% x = 0.5*(xR+1);
% y = 0.5*(yR+1);

 x = [-0.333333333333333
  -0.479308067841920
  -0.479308067841920
  -0.041383864316160
  -0.869739794195568
  -0.869739794195568
  0.739479588391136
  -0.374269007990252
  0.276888377139620
  -0.902619369149368
  -0.374269007990252
  0.276888377139620
  -0.902619369149368];

y = [-0.333333333333333
  -0.479308067841920
  -0.041383864316160
  -0.479308067841920
  -0.869739794195568
  0.739479588391136
  -0.869739794195568
  0.276888377139620
  -0.902619369149368
  -0.374269007990252
  -0.902619369149368
  -0.374269007990252
  0.276888377139620];

w = [-0.299140088935364
 0.351230514866416
 0.351230514866416
 0.351230514866416
 0.106694471217676
 0.106694471217676
 0.106694471217676
 0.154227521780514
 0.154227521780514
 0.154227521780514
 0.154227521780514
 0.154227521780514
 0.154227521780514];

% %% 12-point rule
% x = [-0.501426509658180
%   -0.501426509658180
%   0.002853019316358
%   -0.873821971016996
%   -0.873821971016996
%   0.747643942033992
%   -0.379295097932432
%   0.273004998242798
%   -0.893709900310366
%   -0.379295097932432
%   0.273004998242798
%   -0.893709900310366];
% 
% y =  [-0.501426509658180
%   0.002853019316358
%   -0.501426509658180
%   -0.873821971016996
%   0.747643942033992
%   -0.873821971016996
%   0.273004998242798
%   -0.893709900310366
%   -0.379295097932432
%   -0.893709900310366
%   -0.379295097932432
%   0.273004998242798];
% 
% w = [0.233572551452758
%  0.233572551452758
%  0.233572551452758
%  0.101689812740414
%  0.101689812740414
%  0.101689812740414
%  0.165702151236748
%  0.165702151236748
%  0.165702151236748
%  0.165702151236748
%  0.165702151236748
%  0.165702151236748];

%% transform
x = (x+1)/2;
y = (y+1)/2;
quad_weights = w'/4;

quad_vertices = [x y];




%*** the remaining code is independent of the chosen quadrature rule
nT = size(elements,1);
nQ = size(quad_vertices,1);

%*** first vertices of triangles and corresponding edge vectors
c1  = coordinates(elements(:,1),:);
d21 = coordinates(elements(:,2),:) - c1;
d31 = coordinates(elements(:,3),:) - c1;

%*** compute vector of triangle areas 2*area(T)
area2 = d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1);

%*** build matrix of quadrature vertices by use of affine transformation of Tref
jacobian1 = reshape(repmat([d21(:,1);d31(:,1)],1,nQ)',nT*nQ,2);
jacobian2 = reshape(repmat([d21(:,2);d31(:,2)],1,nQ)',nT*nQ,2);
zref = repmat(quad_vertices,nT,1);
sx = [ sum(zref.*jacobian1,2) sum(zref.*jacobian2,2) ] ...
    + reshape(repmat(c1(:),1,nQ)',nT*nQ,2);

%*** evaluate volume force f at all quadrature vertices z
usx = u(sx);

if(p>=3)
x = quad_vertices(:,1);
y = quad_vertices(:,2);
eta{2} = x; eta{3} = y;

eta{7} = (1-x-y).*eta{2}.*(eta{2}-(1-x-y));
eta{8} = eta{2}.*eta{3}.*(eta{3}-eta{2});
eta{9} = eta{3}.*(1-x-y).*((1-x-y)-eta{3});
% Bubble function degree 3
eta{10} = (1-x-y).*eta{2}.*eta{3};
end

if(p>=4)
% Edge functions, degree 4
eta{11} = (1-x-y).*eta{2}.*( (eta{2}-(1-x-y)).^2-1 );
eta{12} = eta{2}.*eta{3}.*( (eta{3}-eta{2}).^2-1 );
eta{13} = eta{3}.*(1-x-y).*( ((1-x-y)-eta{3}).^2-1 );

% Bubble functions up to degree 4
eta{14} = (1-x-y).*eta{2}.*eta{3}.*((1-x-y)-eta{3});
eta{15} = (1-x-y).*eta{2}.*eta{3}.*(eta{2}-(1-x-y));
end

for j=1:dim
    uEval = reshape(usx(:,j),nQ,nT);
    
    if(p==0)
        err(:,j) = area2.*(quad_weights*(uEval - repmat(uh(:,j)',nQ,1)).^2)';
    elseif(p==1)
        uhP1 = repmat(uh(1:dimU:dimU*nT,j)',nQ,1) + quad_vertices(:,1)*uh(2:dimU:dimU*nT,j)' + quad_vertices(:,2)*uh(3:dimU:dimU*nT,j)';
        err(:,j) = area2.*(quad_weights*(uEval - uhP1).^2)';
    elseif(p==2)
        uhP1 = repmat(uh(1:dimU:dimU*nT,j)',nQ,1) + quad_vertices(:,1)*uh(2:dimU:dimU*nT,j)' + quad_vertices(:,2)*uh(3:dimU:dimU*nT,j)' ...
           + ( (1-quad_vertices(:,1)-quad_vertices(:,2)).*quad_vertices(:,1) )*uh(4:dimU:dimU*nT,j)' ...
           + ( quad_vertices(:,1).*quad_vertices(:,2) )*uh(5:dimU:dimU*nT,j)' ...
           + ( (1-quad_vertices(:,1)-quad_vertices(:,2)).*quad_vertices(:,2) )*uh(6:dimU:dimU*nT,j)';
        err(:,j) = area2.*(quad_weights*(uEval - uhP1).^2)';
    elseif(p==3)
        uhP1 = repmat(uh(1:dimU:dimU*nT,j)',nQ,1) + quad_vertices(:,1)*uh(2:dimU:dimU*nT,j)' + quad_vertices(:,2)*uh(3:dimU:dimU*nT,j)' ...
           + ( (1-quad_vertices(:,1)-quad_vertices(:,2)).*quad_vertices(:,1) )*uh(4:dimU:dimU*nT,j)' ...
           + ( quad_vertices(:,1).*quad_vertices(:,2) )*uh(5:dimU:dimU*nT,j)' ...
           + ( (1-quad_vertices(:,1)-quad_vertices(:,2)).*quad_vertices(:,2) )*uh(6:dimU:dimU*nT,j)' ...
           + eta{7}*uh(7:dimU:dimU*nT,j)' + eta{8}*uh(8:dimU:dimU*nT,j)' + eta{9}*uh(9:dimU:dimU*nT,j)' + eta{10}*uh(10:dimU:dimU*nT,j)';
        err(:,j) = area2.*(quad_weights*(uEval - uhP1).^2)';
    elseif(p==4)
        uhP1 = repmat(uh(1:dimU:dimU*nT,j)',nQ,1) + quad_vertices(:,1)*uh(2:dimU:dimU*nT,j)' + quad_vertices(:,2)*uh(3:dimU:dimU*nT,j)' ...
           + ( (1-quad_vertices(:,1)-quad_vertices(:,2)).*quad_vertices(:,1) )*uh(4:dimU:dimU*nT,j)' ...
           + ( quad_vertices(:,1).*quad_vertices(:,2) )*uh(5:dimU:dimU*nT,j)' ...
           + ( (1-quad_vertices(:,1)-quad_vertices(:,2)).*quad_vertices(:,2) )*uh(6:dimU:dimU*nT,j)' ...
           + eta{7}*uh(7:dimU:dimU*nT,j)' + eta{8}*uh(8:dimU:dimU*nT,j)' + eta{9}*uh(9:dimU:dimU*nT,j)' + eta{10}*uh(10:dimU:dimU*nT,j)' ...
           + eta{11}*uh(11:dimU:dimU*nT,j)' + eta{12}*uh(12:dimU:dimU*nT,j)' + eta{13}*uh(13:dimU:dimU*nT,j)' + eta{14}*uh(14:dimU:dimU*nT,j)' + eta{15}*uh(15:dimU:dimU*nT,j)';
        err(:,j) = area2.*(quad_weights*(uEval - uhP1).^2)';        
    end
end
    

err = abs(err);