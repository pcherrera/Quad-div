function proj = computeL2projNeumannPp(coordinates,elements,g,pU)
%% compute L2 proj. to C1 elements on boundary,
%
%   p=0,1,2,3 possible
%
% local basis is given by the functions
%   1, 2t-1, 4t(1-t), 4t(1-t)(2t-1),
%
% ordering of dofs: dofs according to coordinates field, then to elements
% field (first #ELE corresponds to p=2, etc.)

nE = size(elements,1);
nC = size(coordinates,1);

area = sqrt( sum((coordinates(elements(:,1),:)-...
    coordinates(elements(:,2),:)).^2,2) );

dimU = (pU+1)*nE;

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
g_eval_xi = zeros(nE*gauss_order,2);

a = coordinates(elements(:,1),:);
b = coordinates(elements(:,2),:);
%normal vector
z = (a-b); 
n = [z(:,2),-z(:,1)]./sqrt((z(:,1).^2+z(:,2).^2)); 
%*** basis functions
if (pU>=0)
    eta{1} = 1 ;
end 
if (pU>=1)
    eta{2} = 2*xi-1;
end
if(pU>=2)
    eta{3} = 4*xi.*(1-xi);
end
if(pU>=3)
    eta{4} = 4*xi.*(1-xi).*(2*xi-1);
end


for i=1:gauss_order
    g_eval_xi(nE*(i-1)+1:nE*i,(1:2)) = g( (b-a)*xi(i) + a );
end
g_eval_xi1 = g_eval_xi(:,1);
g_eval_xi2 = g_eval_xi(:,2);
g_eval_xi1 = reshape(g_eval_xi1,nE,gauss_order);
g_eval_xi2 = reshape(g_eval_xi2,nE,gauss_order);
n1  = repmat(n(:,1),1,gauss_order);
n2  = repmat(n(:,2),1,gauss_order);
g_eval_xi = g_eval_xi1.*n1+g_eval_xi2.*n2;

eta{1} = ones(nE,gauss_order);
gEta{1} = area.*( (eta{1}.*g_eval_xi)*(w'));

if (pU>=1)
    for j=2:pU+1
        eta{j} = repmat(eta{j},nE,1);
        gEta{j} = area.*( (eta{j}.*g_eval_xi)*(w'));
    end
end

%% RHS
eta_g = gEta{1};

for j=2:pU+1
    eta_g = [eta_g;gEta{j}];
end

%% Sys matrix
L2 = [   1,     0,  2/3,        0; ...
         0,   1/3,    0,     2/15; ...
       2/3,     0, 8/15,        0; ...
         0,  2/15,    0,    8/105; ...
 ];

L2 = L2(1:pU+1,1:pU+1);

L = repmat(L2(:),nE,1).*reshape(repmat(area',(pU+1)^2,1),(pU+1)^2*nE,1);

I = zeros((pU+1)^2*nE,1);
J = I;



for j=1:nE
    Itmp = repmat([j (nE:nE:(pU+1)*nE-1)+j]',1,pU+1);
    Jtmp = Itmp';
    
    I( (j-1)*(pU+1)^2+1:j*(pU+1)^2) = Itmp(:);
    J( (j-1)*(pU+1)^2+1:j*(pU+1)^2) = Jtmp(:);
end

L = sparse(I,J,L,dimU,dimU);

% compute L^2-proj
proj = L\eta_g;