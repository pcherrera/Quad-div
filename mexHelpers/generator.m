function generator(filename)

syms x y a b c d
vars = [x,y];

p = 0;
q = 3;

eta{1} = 1;
eta{2} = x;
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

fid = fopen(filename, 'w');
fprintf(fid,'// Automated generated file\n');
fprintf(fid,'#ifndef _FUN_SIGMAGDIVTAU_H_\n');
fprintf(fid,'#define _FUN_SIGMAGDIVTAU_H_\n');
for i=0:p
    for j=1:q+1
        fprintf(fid,'inline void compTsigmaGdivTauP%dQ%d(double * M, double a, double b, double c, double d)\n',i,j);
    end
end

fprintf(fid,'// Automated generated file\n');
fprintf(fid,'#include "funSigmaGdivTau.h"\n');
for i=0:p
    for j=1:q+1
        fprintf(fid,'void compTsigmaGradVP%dQ%d(double * M, double a, double b, double c, double d) {\n',i,j);
        dimI = (i+1)*(i+2);
        dimJ = (j+1)*(j+2);
        for r=0:dimI*dimJ-1 
            r1 = rem(r,dimI)+1;
            r2 = rem(r,dimJ)+1;
            fprintf(fid,'M[%d]=%s\n',r,int(int(dot(etaVec{r1},gradient(divergence(etaVec{r2},vars),vars)),y,0,1-x),x,0,1));            
            %fprintf('M[%d]=%s\n',r,int(int(dot(etaVec{r1},etaVec{r2}),y,0,1-x),x,0,1));            
        end
        fprintf(fid,'}\n');
    end
end 
fclose(fid);