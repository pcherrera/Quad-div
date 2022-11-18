function fh = computeL2projToPp(coordinates,elements,f,p)

%*** We use the same hierarchical basis as for the ansatz functions, that
%is
%
% eta1 = 1, eta2 = x, eta3 = y,
% eta4 = (1-x-y)*x eta5 = x*y, eta6 = y*(1-x-y)
% % Edge functions, degree 3
% eta{7} = eta{1}*eta{2}*(eta{2}-eta{1});
% eta{8} = eta{2}*eta{3}*(eta{3}-eta{2});
% eta{9} = eta{3}*eta{1}*(eta{1}-eta{3});
% 
% % Bubble function degree 3
% eta{10} = eta{1}*eta{2}*eta{3};


%*** get inverse of L2 mass matrix on reference element
if(p==0)
    L2inv = 2;
elseif(p==1)
    L2inv = [18   -24   -24; -24 48 24; -24 24 48];
elseif(p==2)
    L2inv = [72      -60         -60        -180           0        -180; ...
            -60      120          60           0        -180         180; ...
            -60       60         120         180        -180           0; ...
            -180       0         180        1080         180         180; ...
            0       -180        -180         180        1080         180; ...
            -180     180           0         180         180        1080];
elseif(p==3)
    L2inv = [  200,  -220, -220,  -420,     0,  -420,   560,    0, -560,  1120; ...
             -220,   440,  220,     0,  -420,   420, -1120,  560,  560,     0; ...
             -220,   220,  440,   420,  -420,     0,  -560, -560, 1120,     0; ...
             -420,     0,  420,  2520,  1260,  1260,     0, -840,  840, -8400; ...
                0,  -420, -420,  1260,  2520,  1260,   840,    0, -840, -8400; ...
             -420,   420,    0,  1260,  1260,  2520,  -840,  840,    0, -8400; ...
              560, -1120, -560,     0,   840,  -840,  5600, -280, -280,     0; ...
                0,   560, -560,  -840,     0,   840,  -280, 5600, -280,     0; ...
             -560,   560, 1120,   840,  -840,     0,  -280, -280, 5600,     0; ...
             1120,     0,    0, -8400, -8400, -8400,     0,    0,    0, 58800];
elseif(p==4)
    L2inv = [   450,  -420,   -420,  -2520,      0,  -2520,   1260,     0,  -1260,   1470,   -1575,       0,   -1575,      3675,     -3675; ...
              -420,   840,    420,      0,  -2520,   2520,  -2520,  1260,   1260,      0,       0,   -1575,    1575,         0,     11025; ...
              -420,   420,    840,   2520,  -2520,      0,  -1260, -1260,   2520,      0,    1575,   -1575,       0,    -11025,         0; ...
             -2520,     0,   2520,  40320,   7560,   7560,      0, -5040,   5040, -12600,   31500,    3150,    3150,    -31500,    -15750; ...
                 0, -2520,  -2520,   7560,  40320,   7560,   5040,     0,  -5040, -12600,    3150,   31500,    3150,     15750,    -15750; ...
             -2520,  2520,      0,   7560,   7560,  40320,  -5040,  5040,      0, -12600,    3150,    3150,   31500,     15750,     31500; ...
              1260, -2520,  -1260,      0,   5040,  -5040,  12600, -3780,  -3780,      0,       0,    1575,   -1575,         0,    -61425; ...
                 0,  1260,  -1260,  -5040,      0,   5040,  -3780, 12600,  -3780,      0,   -1575,       0,    1575,     61425,     61425; ...
             -1260,  1260,   2520,   5040,  -5040,      0,  -3780, -3780,  12600,      0,    1575,   -1575,       0,    -61425,         0; ...
              1470,     0,      0, -12600, -12600, -12600,      0,     0,      0,  59850,   -3150,   -3150,   -3150,         0,         0; ...
             -1575,     0,   1575,  31500,   3150,   3150,      0, -1575,   1575,  -3150, 55125/2,  1575/4,  1575/4,   -1575/2,   -1575/4; ...
                 0, -1575,  -1575,   3150,  31500,   3150,   1575,     0,  -1575,  -3150,  1575/4, 55125/2,  1575/4,    1575/4,   -1575/4; ...
             -1575,  1575,      0,   3150,   3150,  31500,  -1575,  1575,      0,  -3150,  1575/4,  1575/4, 55125/2,    1575/4,    1575/2; ...
              3675,     0, -11025, -31500,  15750,  15750,      0, 61425, -61425,      0, -1575/2,  1575/4,  1575/4, 1472625/2, 1472625/4; ...
             -3675, 11025,      0, -15750, -15750,  31500, -61425, 61425,      0,      0, -1575/4, -1575/4,  1575/2, 1472625/4, 1472625/2];
 
else
    error('Pol. degree %d not supported',p)
end

%*** quadrature rule on reference triangle Tref = conv{(0,0),(1,0),(0,1)}
tmp_pos = [6-sqrt(15) ; 9+2*sqrt(15) ; 6+sqrt(15) ; 9-2*sqrt(15) ; 7]/21;
quad_vertices = tmp_pos([1 1 ; 2 1 ; 1 2 ; 3 4 ; 3 3 ; 4 3 ; 5 5]);
tmp_wts = [155-sqrt(15) 155+sqrt(15) 270]/2400;
w = tmp_wts([1 1 1 2 2 2 3])';

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



%% transform
x = (x+1)/2;
y = (y+1)/2;
w = w/4;

quad_vertices = [x y];

%*** the remaining code is independent of the chosen quadrature rule
nE = size(elements,1);
nQ = size(quad_vertices,1);

%*** first vertices of triangles and corresponding edge vectors
c1  = coordinates(elements(:,1),:);
d21 = coordinates(elements(:,2),:) - c1;
d31 = coordinates(elements(:,3),:) - c1;

%*** basis functions evaluated at quadrature nodes
x = quad_vertices(:,1); y = quad_vertices(:,2);

eta{1} = ones(size(x));
eta{2} = x;
eta{3} = y;

% Edge functions, degree 2
eta{4} = (1-x-y).*x; %edge 1
eta{5} = x.*y; %edge 2
eta{6} = y.*(1-x-y); %edge 3

if(p>=3)
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

%*** evalutate f at Transformed nodes
nodes = [d21(:,1)*x'+d31(:,1)*y'; d21(:,2)*x'+d31(:,2)*y']+repmat(c1(:),1,nQ);
nodes = reshape(nodes',nE*nQ,2);
f_T = reshape(feval(f,nodes)',nQ,nE)';

%*** dimension of basis
dimV = (p+1)*(p+2)/2;

%*** integrate
% erg = [repmat(eta_1',nE,1).*f_T; repmat(eta_2',nE,1).*f_T; repmat(eta_3',nE,1).*f_T; ...
%     repmat(eta_4',nE,1).*f_T; repmat(eta_5',nE,1).*f_T; repmat(eta_6',nE,1).*f_T];
erg = [];
for k=1:dimV
    erg = [erg; repmat(eta{k}',nE,1).*f_T];
end
erg = L2inv*reshape(erg*w,nE,dimV)';
fh = erg(:);


    

    