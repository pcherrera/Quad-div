function proj = computeL2projBOUC1(coordinates,elements,g,pU)
%% compute L2 proj. to C1 elements on boundary,
%
%   p=1,2,3,4 possible
%
% local basis is given by the functions
%   1-t, t, t(1-t), t(1-t)(2t-1), t(1-t)( (2t-1)^2-1 )
%
% ordering of dofs: dofs according to coordinates field, then to elements
% field (first #ELE corresponds to p=2, etc.)

nE = size(elements,1);
nC = size(coordinates,1);

area = sqrt( sum((coordinates(elements(:,1),:)-...
    coordinates(elements(:,2),:)).^2,2) );

dimU = nC+(pU-1)*nE;

% <\eta_i, g>
% mit gauss-quad.
% w = [0.236926885056189 0.478628670499366 0.568888888888889 0.478628670499366 0.236926885056189];
% xi = [-0.906179845938664 -0.538469310105683 0 0.538469310105683 0.906179845938664];
xi = [-0.9894009349916499325961542, ...
    -0.9445750230732325760779884, ...
    -0.8656312023878317438804679, ...
    -0.7554044083550030338951012, ...
    -0.6178762444026437484466718, ...
    -0.4580167776572273863424194, ...
    -0.2816035507792589132304605, ...
    -0.0950125098376374401853193, ...
    0.0950125098376374401853193, ...
    0.2816035507792589132304605, ...
    0.4580167776572273863424194, ...
    0.6178762444026437484466718, ...
    0.7554044083550030338951012, ...
    0.8656312023878317438804679, ...
    0.9445750230732325760779884, ...
    0.9894009349916499325961542];

w = [0.0271524594117540948517806; ...
    0.0622535239386478928628438; ...
    0.0951585116824927848099251; ...
    0.1246289712555338720524763; ...
    0.1495959888165767320815017; ...
    0.1691565193950025381893121; ...
    0.1826034150449235888667637; ...
    0.1894506104550684962853967; ...
    0.1894506104550684962853967; ...
    0.1826034150449235888667637; ...
    0.1691565193950025381893121; ...
    0.1495959888165767320815017; ...
    0.1246289712555338720524763; ...
    0.0951585116824927848099251; ...
    0.0622535239386478928628438; ...
    0.0271524594117540948517806]';

%*** transform to [0,1]
xi = 0.5*(xi+1);
w = w/2;

gauss_order = length(w);

%eta_g = zeros(dimU,1);
g_eval_xi = zeros(nE,gauss_order);

a = coordinates(elements(:,1),:);
b = coordinates(elements(:,2),:);

%*** basis functions
eta{1} = 1-xi;
eta{2} = xi;

if(pU>=2)
    eta{3} = xi.*(1-xi);
end
if(pU>=3)
    eta{4} = xi.*(1-xi).*(2*xi-1);
end
if(pU>=4)
    eta{5} = xi.*(1-xi).*( (2*xi-1).^2-1);
end



% Evaluate g at quadrature points
if(nargin(g)==3)
    for i=1:gauss_order
        g_eval_xi(:,i) = g( (b-a)*xi(i) + a ,a,b );
    end
else % assume g = g(x)
    for i=1:gauss_order
        g_eval_xi(:,i) = g( (b-a)*xi(i) + a );
    end
end

for j=1:pU+1
    eta{j} = repmat(eta{j},nE,1);
    gEta{j} = area.*( (eta{j}.*g_eval_xi)*(w'));
end

%% RHS
eta_g = accumarray(elements(:),[gEta{1};gEta{2}],[nC 1]);

for j=2:pU
    eta_g = [eta_g;gEta{j+1}];
end

%% Sys matrix
L2 = [   1/3,   1/6,  1/12, -1/60, -1/15; ...
   1/6,   1/3,  1/12,  1/60, -1/15; ...
  1/12,  1/12,  1/30,     0, -1/35; ...
 -1/60,  1/60,     0, 1/210,     0; ...
 -1/15, -1/15, -1/35,     0, 8/315];

L2 = L2(1:pU+1,1:pU+1);

L = repmat(L2(:),nE,1).*reshape(repmat(area',(pU+1)^2,1),(pU+1)^2*nE,1);

I = zeros((pU+1)^2*nE,1);
J = I;


for j=1:nE
    Itmp = repmat([elements(j,1:2) nC+(0:nE:(pU-1)*nE-1)+j]',1,pU+1);
    Jtmp = Itmp';
    
    I( (j-1)*(pU+1)^2+1:j*(pU+1)^2) = Itmp(:);
    J( (j-1)*(pU+1)^2+1:j*(pU+1)^2) = Jtmp(:);
end

L = sparse(I,J,L,nC+(pU-1)*nE,nC+(pU-1)*nE);

% compute L^2-proj
proj = L\eta_g;