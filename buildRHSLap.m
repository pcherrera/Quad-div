function bVec = buildRHSLap(coordinates,elements,qV1,qV2,f)

nE = size(elements,1);

%*** rhs vector
dimV1 = (qV1+1)*(qV1+2);
dimV2 = (qV2+1)*(qV2+2);
dimV = dimV1+dimV2;

bVec = zeros(dimV*nE,1);


%% compute <f,v>_\Omega
%*** by using 7-point gauss-quadrature on reference triangle
%w = [0.1125 0.0661970764 0.0661970764 0.0661970764 0.0629695903 0.0629695903 0.0629695903]';
%x = [0.333333333 0.470142064 0.059715872 0.470142064 0.101286507 0.797426985 0.101286507]';
%y = [0.333333333 0.470142064 0.470142064 0.059715872 0.101286507 0.101286507 0.797426985]';
%gauss_order = length(w);

%*** quadrature rule on reference triangle Tref = conv{(0,0),(1,0),(0,1)}
% tmp_pos = [6-sqrt(15) ; 9+2*sqrt(15) ; 6+sqrt(15) ; 9-2*sqrt(15) ; 7]/21;
% quad_vertices = tmp_pos([1 1 ; 2 1 ; 1 2 ; 3 4 ; 3 3 ; 4 3 ; 5 5]);
% x = quad_vertices(:,1);
% y = quad_vertices(:,2);
% tmp_wts = [155-sqrt(15) 155+sqrt(15) 270]/2400;
% w = tmp_wts([1 1 1 2 2 2 3])';

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



gauss_order = length(w);

%% transform
x = (x+1)/2;
y = (y+1)/2;
w = w/4;

% define basis functions on reference triangle
% ordered by, i=1,2,3 correspond to nodes in array elements
%             i=4,5,6 correspont to edges of a triangle, where 
%                [elements(j,1) elements(j,2)] is first edge (i=4), etc...
% and evaluate at gaussian quadrature nodes
% eta_1 = 1-x-y;
% eta_2 = x;
% eta_3 = y;
% eta_4 = (1-x-y).*x;
% eta_5 = x.*y;
% eta_6 = (1-x-y).*y;

%
eta{2} = 1-x-y;
eta{1} = x;
eta{3} = y;

% Define kernel functions for lobatto shape function
% phi0 = 2*sqrt(3/2);
% phi1 = -2*sqrt(5/2);
% phi2 = -1/2*sqrt(7/2);
phi0 = 1; phi1 = 1; phi2 = 1; phi3 = 1;

% Edge functions, degree 2
eta{4} = eta{1}.*eta{2}*phi0; %edge 1
eta{5} = eta{2}.*eta{3}*phi0; %edge 2
eta{6} = eta{3}.*eta{1}*phi0; %edge 3

% Edge functions, degree 3
eta{7} = eta{1}.*eta{2}.*(eta{2}-eta{1})*phi1;
eta{8} = eta{2}.*eta{3}.*(eta{3}-eta{2})*phi1;
eta{9} = eta{3}.*eta{1}.*(eta{1}-eta{3})*phi1;

% Bubble function degree 3
eta{10} = eta{1}.*eta{2}.*eta{3}*phi0^2;

% Edge functions, degree 4
eta{11} = eta{1}.*eta{2}.*( (eta{2}-eta{1}).^2-1 )*phi2;
eta{12} = eta{2}.*eta{3}.*( (eta{3}-eta{2}).^2-1 )*phi2;
eta{13} = eta{3}.*eta{1}.*( (eta{1}-eta{3}).^2-1 )*phi2;

% Bubble functions up to degree 4
eta{14} = eta{1}.*eta{2}.*eta{3}.*(eta{1}-eta{3})*phi0*phi1;
eta{15} = eta{1}.*eta{2}.*eta{3}.*(eta{2}-eta{1})*phi0*phi1;

% Edge functions, degree 5
eta{16} = eta{1}.*eta{2}.*( 7*(eta{2}-eta{1}).^2-3 ).*(eta{2}-eta{1});
eta{17} = eta{2}.*eta{3}.*( 7*(eta{3}-eta{2}).^2-3 ).*(eta{3}-eta{2});
eta{18} = eta{3}.*eta{1}.*( 7*(eta{1}-eta{3}).^2-3 ).*(eta{1}-eta{3});

% Bubble functions up to degree 5
eta{19} = eta{1}.*eta{2}.^3.*eta{3};
eta{20} = eta{1}.^2.*eta{2}.^2.*eta{3};
eta{21} = eta{1}.^3.*eta{2}.*eta{3};

%Build Basis functions for P(T)^2
k=0;
for j=1:2:42
    k=k+1;
    etaVec{j} = [eta{k} zeros(size(x))];
    etaVec{j+1} = [zeros(size(x)) eta{k}];
end


%*** First vertex of elements and corresponding edge vectors 
c1 = coordinates(elements(:,1),:);
d21 = coordinates(elements(:,2),:)-c1;
d31 = coordinates(elements(:,3),:)-c1;
%*** Vector of element areas 2*|T|
area2 = (d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1));
% 
% %*** evalutate f at Transformed nodes
nodes = [d21(:,1)*x'+d31(:,1)*y'; d21(:,2)*x'+d31(:,2)*y']+repmat(c1(:),1,gauss_order);
nodes = reshape(nodes',nE*gauss_order,2);
%f_T = reshape(feval(f,nodes)',gauss_order,nE)';
f_T =  feval(f,nodes);
f_T1 = -reshape(f_T(:,1)',gauss_order,nE)';
f_T2 = -reshape(f_T(:,2)',gauss_order,nE)';

% %*** integrate
% erg = [repmat(eta_1',nE,1).*f_T; repmat(eta_2',nE,1).*f_T; repmat(eta_3',nE,1).*f_T; ...
%     repmat(eta_4',nE,1).*f_T; repmat(eta_5',nE,1).*f_T; repmat(eta_6',nE,1).*f_T];
erg = [];
for k=1:dimV1
    erg = [erg; repmat(etaVec{k}(:,1)',nE,1).*f_T1 + repmat(etaVec{k}(:,2)',nE,1).*f_T2];
end
erg = ( repmat(area2,1,dimV1).*reshape(erg*w,nE,dimV1) )';
bVec(1:dimV1*nE) = erg(:);

% for dog = 1:dimV1
%     for E = 1:nE
%         T_area2 = area2(E);
%         for q = 1:gauss_order
%             f_at_node = f([nodes(E,q), nodes(nE+E,q)]);
%             bVec(nE*(dog-1)+E) = bVec(nE*(dog-1)+E) + T_area2*w(q)*(f_at_node(1,1)*etaVec{dog}(q,1)+f_at_node(1,2)*etaVec{dog}(q,2));
%         end 
%     end
% end
