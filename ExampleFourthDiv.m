%% Fourth Order Div problem
% ultra-weak DPG formulation
% 
% MODEL PROBLEM is
% % grad div(u2) + u1  = f    in Omega
% % grad div(u1) - u2  = 0    in Omega
% %             div u1 = 0    on Gamma
% %               n.u1 = 0    on Gamma

%%
clear all 
close all

%addpath Data/
% addpath Data1/
load Meshes/quadMesh.mat
h=1;
coordinates = h*coordinates;
%% Example 1: Smooth solution on square domain (0,1)^2
%
% U = @(x)[x(:,1).^2.*(x(:,1)-1).^2.*x(:,2).^2.*(x(:,2)-1).^2 ... 
%     x(:,1).^2.*(x(:,1)-1).^2.*x(:,2).^2.*(x(:,2)-1).^2];
% U_1 = @(x) x(:,1).^2.*(x(:,1)-1).^2.*x(:,2).^2.*(x(:,2)-1).^2;
% U_2 = @(x) x(:,1).^2.*(x(:,1)-1).^2.*x(:,2).^2.*(x(:,2)-1).^2;
% 
% divU=@(x) 2.*x(:,1).*x(:,2).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ... 
%     + 2.*x(:,1).^2.*x(:,2).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ... 
%     + x(:,1).^2.*x(:,2).^2.*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ...
%     + x(:,1).^2.*x(:,2).^2.*(2.*x(:,2) - 2).*(x(:,1) - 1).^2;
% 
% gdivU = @(x) [2.*x(:,1).^2.*x(:,2).^2.*(x(:,2) - 1).^2 ... 
%     + 2.*x(:,2).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ... 
%     + 4.*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ...
%     + 2.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 ... 
%     + 2.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ...
%     + x(:,1).^2.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) ...
%     + 4.*x(:,1).*x(:,2).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ...
%      2.*x(:,1).^2.*x(:,2).^2.*(x(:,1) - 1).^2 + 2.*x(:,1).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ... 
%     + 2.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 2.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%     + 4.*x(:,1).^2.*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + x(:,1).^2.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) ...
%     + 4.*x(:,1).*x(:,2).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2];
% 
% gdivU_1 = @(x) 2.*x(:,1).^2.*x(:,2).^2.*(x(:,2) - 1).^2 ... 
%     + 2.*x(:,2).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ... 
%     + 4.*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ...
%     + 2.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 ... 
%     + 2.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ...
%     + x(:,1).^2.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) ...
%     + 4.*x(:,1).*x(:,2).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2;
% 
% gdivU_2 = @(x) 2.*x(:,1).^2.*x(:,2).^2.*(x(:,1) - 1).^2 + 2.*x(:,1).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ... 
%     + 2.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 2.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%     + 4.*x(:,1).^2.*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + x(:,1).^2.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) ...
%     + 4.*x(:,1).*x(:,2).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2;
% 
% 
% dgdivU = @(x) 2.*x(:,1).^2.*x(:,2).^2.*(2.*x(:,1) - 2) + 4.*x(:,1).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ... 
%     + 2.*x(:,1).^2.*x(:,2).^2.*(2.*x(:,2) - 2) + 4.*x(:,2).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ... 
%     + 2.*x(:,1).^2.*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + 6.*x(:,1).^2.*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 ... 
%     + 6.*x(:,2).^2.*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + 2.*x(:,2).^2.*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 ...
%     + 4.*x(:,1).*x(:,2).^2.*(x(:,1) - 1).^2 + 12.*x(:,1).^2.*x(:,2).*(x(:,1) - 1).^2 + 12.*x(:,1).*x(:,2).^2.*(x(:,2) - 1).^2 ...
%     + 4.*x(:,1).^2.*x(:,2).*(x(:,2) - 1).^2 + 4.*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) ...
%     + 4.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 8.*x(:,1).*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ...
%     + 8.*x(:,1).*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2;
% 
% gdgdivU =@(x) [4.*x(:,1).^2.*x(:,2).^2 + 4.*x(:,1).^2.*(x(:,2) - 1).^2 + 4.*x(:,2).^2.*(x(:,1) - 1).^2 + 24.*x(:,2).^2.*(x(:,2) - 1).^2 ... 
%     + 4.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + 8.*x(:,1).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%     + 12.*x(:,1).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 12.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%     + 8.*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 24.*x(:,1).*x(:,2).*(x(:,1) - 1).^2 + ...
%     24.*x(:,1).*x(:,2).*(x(:,2) - 1).^2 + 6.*x(:,1).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + ...
%     +6.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 8.*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2) ... 
%     + 12.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2) + 12.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2) ... 
%     + 8.*x(:,1).^2.*x(:,2).*(2.*x(:,2) - 2) + 16.*x(:,1).*x(:,2).*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) ...
%     , 4.*x(:,1).^2.*x(:,2).^2 + 24.*x(:,1).^2.*(x(:,1) - 1).^2 + 4.*x(:,1).^2.*(x(:,2) - 1).^2 + 4.*x(:,2).^2.*(x(:,1) - 1).^2 ... 
%     + 4.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + 8.*x(:,1).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%     + 12.*x(:,1).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 12.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%     + 8.*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 24.*x(:,1).*x(:,2).*(x(:,1) - 1).^2 + 24.*x(:,1).*x(:,2).*(x(:,2) - 1).^2 ...
%     + 6.*x(:,1).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 6.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) ...
%     + 8.*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2) + 12.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2) + 12.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2) ... 
%     + 8.*x(:,1).^2.*x(:,2).*(2.*x(:,2) - 2) + 16.*x(:,1).*x(:,2).*(2.*x(:,1) - 2).*(2.*x(:,2) - 2)];
% 
% gdgdivU_1 = @(x) 4.*x(:,1).^2.*x(:,2).^2 + 4.*x(:,1).^2.*(x(:,2) - 1).^2 + 4.*x(:,2).^2.*(x(:,1) - 1).^2 + 24.*x(:,2).^2.*(x(:,2) - 1).^2 ... 
%     + 4.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + 8.*x(:,1).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%     + 12.*x(:,1).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 12.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%     + 8.*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 24.*x(:,1).*x(:,2).*(x(:,1) - 1).^2 + ...
%     24.*x(:,1).*x(:,2).*(x(:,2) - 1).^2 + 6.*x(:,1).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + ...
%     +6.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 8.*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2) ... 
%     + 12.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2) + 12.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2) ... 
%     + 8.*x(:,1).^2.*x(:,2).*(2.*x(:,2) - 2) + 16.*x(:,1).*x(:,2).*(2.*x(:,1) - 2).*(2.*x(:,2) - 2);
% gdgdivU_2 = @(x) 4.*x(:,1).^2.*x(:,2).^2 + 24.*x(:,1).^2.*(x(:,1) - 1).^2 + 4.*x(:,1).^2.*(x(:,2) - 1).^2 + 4.*x(:,2).^2.*(x(:,1) - 1).^2 ... 
%     + 4.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + 8.*x(:,1).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%     + 12.*x(:,1).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 12.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%     + 8.*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 24.*x(:,1).*x(:,2).*(x(:,1) - 1).^2 + 24.*x(:,1).*x(:,2).*(x(:,2) - 1).^2 ...
%     + 6.*x(:,1).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 6.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) ...
%     + 8.*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2) + 12.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2) + 12.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2) ... 
%     + 8.*x(:,1).^2.*x(:,2).*(2.*x(:,2) - 2) + 16.*x(:,1).*x(:,2).*(2.*x(:,1) - 2).*(2.*x(:,2) - 2);
% 
% f = @(x) (gdgdivU(x)+U(x));
% u0 = @(x) U(x); 



%% Example 2: Smooth solution on square domain (0,1)..^2
%   U = @(x) [ x(:,1).^2.*x(:,2).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2, sin(pi.*x(:,1)).^2.*sin(pi.*x(:,2)).^2 ];
%   U_1 = @(x) x(:,1).^2.*x(:,2).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2;
%   U_2 = @(x) sin(pi.*x(:,1)).^2.*sin(pi.*x(:,2)).^2;
%   divU = @(x) 2.*x(:,1).*x(:,2).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + x(:,1).^2.*x(:,2).^2.*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + 2.*pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1)).^2.*sin(pi.*x(:,2)); 
%   gdivU = @(x) [2.*x(:,1).^2.*x(:,2).^2.*(x(:,2) - 1).^2 ...
%       + 2*x(:,2).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ...
%       + 4*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ...
%       + 4*pi^2.*cos(pi.*x(:,1)).*cos(pi.*x(:,2)).*sin(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
%            2*pi.^2.*cos(pi.*x(:,2)).^2.*sin(pi.*x(:,1)).^2 ... 
%           - 2*pi.^2.*sin(pi.*x(:,1)).^2.*sin(pi.*x(:,2)).^2 ...
%           + 2*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 ...
%           + 2*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 ... 
%           + x(:,1).^2.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) ...
%           + 4.*x(:,1).*x(:,2).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2];
%   
%  gdiv U in components
%   gdivU_1 = @(x) 2.*x(:,1).^2.*x(:,2).^2.*(x(:,2) - 1).^2 + 2.*x(:,2).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + 4.*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + 4.*pi.^2.*cos(pi.*x(:,1)).*cos(pi.*x(:,2)).*sin(pi.*x(:,1)).*sin(pi.*x(:,2));
%   gdivU_2 = @(x) 2.*pi.^2.*cos(pi.*x(:,2)).^2.*sin(pi.*x(:,1)).^2 - 2.*pi.^2.*sin(pi.*x(:,1)).^2.*sin(pi.*x(:,2)).^2 + 2.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 2.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + x(:,1).^2.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 4.*x(:,1).*x(:,2).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 ;
%  
%        
%    
%   dgdivU = @(x) 2.*x(:,1).^2.*x(:,2).^2.*(2.*x(:,1) - 2) + 4.*x(:,1).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + 2.*x(:,1).^2.*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + 6.*x(:,2).^2.*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + 4.*x(:,1).*x(:,2).^2.*(x(:,1) - 1).^2 + 12.*x(:,1).*x(:,2).^2.*(x(:,2) - 1).^2 + 4.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 4.*pi.^3.*cos(pi.*x(:,1)).^2.*cos(pi.*x(:,2)).*sin(pi.*x(:,2)) - 12.*pi.^3.*cos(pi.*x(:,2)).*sin(pi.*x(:,1)).^2.*sin(pi.*x(:,2)) + 8.*x(:,1).*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2;
%   
%   gdgdivU = @(x)  [4.*x(:,1).^2.*x(:,2).^2 + 4.*x(:,1).^2.*(x(:,2) - 1).^2 + 4.*x(:,2).^2.*(x(:,1) - 1).^2 + 24.*x(:,2).^2.*(x(:,2) - 1).^2 + 4.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + 8.*x(:,1).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + 8.*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 8.*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2) + 8.*x(:,1).^2.*x(:,2).*(2.*x(:,2) - 2) + 16.*x(:,1).*x(:,2).*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) - 32.*pi.^4.*cos(pi.*x(:,1)).*cos(pi.*x(:,2)).*sin(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
%                12.*pi.^4.*sin(pi.*x(:,1)).^2.*sin(pi.*x(:,2)).^2 + 12.*x(:,1).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 12.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + 24.*x(:,1).*x(:,2).*(x(:,1) - 1).^2 + 24.*x(:,1).*x(:,2).*(x(:,2) - 1).^2 + 6.*x(:,1).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 6.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 4.*pi.^4.*cos(pi.*x(:,1)).^2.*cos(pi.*x(:,2)).^2 + 12.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2) + 12.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2) - 4.*pi.^4.*cos(pi.*x(:,1)).^2.*sin(pi.*x(:,2)).^2 - 12.*pi.^4.*cos(pi.*x(:,2)).^2.*sin(pi.*x(:,1)).^2];
%   f = @(x) (gdgdivU(x) + U(x));
% 
%   grad div grad div in components
%   gdgdivU_1 = @(x)  4.*x(:,1).^2.*x(:,2).^2 + 4.*x(:,1).^2.*(x(:,2) - 1).^2 + 4.*x(:,2).^2.*(x(:,1) - 1).^2 + 24.*x(:,2).^2.*(x(:,2) - 1).^2 + 4.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + 8.*x(:,1).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + 8.*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 8.*x(:,1).*x(:,2).^2.*(2.*x(:,1) - 2) + 8.*x(:,1).^2.*x(:,2).*(2.*x(:,2) - 2) + 16.*x(:,1).*x(:,2).*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) - 32.*pi.^4.*cos(pi.*x(:,1)).*cos(pi.*x(:,2)).*sin(pi.*x(:,1)).*sin(pi.*x(:,2)); 
%   gdgdivU_2 = @(x)  12.*pi.^4.*sin(pi.*x(:,1)).^2.*sin(pi.*x(:,2)).^2 + 12.*x(:,1).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2 + 12.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 + 24.*x(:,1).*x(:,2).*(x(:,1) - 1).^2 + 24.*x(:,1).*x(:,2).*(x(:,2) - 1).^2 + 6.*x(:,1).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 6.*x(:,2).^2.*(2.*x(:,1) - 2).*(2.*x(:,2) - 2) + 4.*pi.^4.*cos(pi.*x(:,1)).^2.*cos(pi.*x(:,2)).^2 + 12.*x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2) + 12.*x(:,1).*x(:,2).^2.*(2.*x(:,2) - 2) - 4.*pi.^4.*cos(pi.*x(:,1)).^2.*sin(pi.*x(:,2)).^2 - 12.*pi.^4.*cos(pi.*x(:,2)).^2.*sin(pi.*x(:,1)).^2;
%   
%   
%   u0 = @(x) U(x);
   

%% Example 3 
% U = @(x) [ x(:,1).*x(:,2).*(x(:,1) - 1) + x(:,1).*(x(:,1) - 1).*(x(:,2) - 1), - x(:,1).*x(:,2).*(x(:,2) - 1) - x(:,2).*(x(:,1) - 1).*(x(:,2) - 1)]; 
% U_1 = @(x) x(:,1).*x(:,2).*(x(:,1) - 1) + x(:,1).*(x(:,1) - 1).*(x(:,2) - 1);
% U_2 = @(x) - x(:,1).*x(:,2).*(x(:,2) - 1) - x(:,2).*(x(:,1) - 1).*(x(:,2) - 1);
% divU = @(x) 0*x(:,1)+0*x(:,2);
% gdivU = @(x) 0*[x(:,1),x(:,2)];
% gdivU_1 = @(x) 0*x(:,1)+0*x(:,2);
% gdivU_2 = @(x) 0*x(:,1)+0*x(:,2);
% dgdivU = @(x) 0*x(:,1)+0*x(:,2);
% f = @(x) (U(x));
% u0 = @(x) (U(x));   

%% Example 4: Singular solution on Lshape domain
% load Meshes/LshapeMesh.mat 
% 
% U = @(x) gradULSHAPE(x);
% U_1 = @(x) gradU_1LSHAPE(x);
% U_2 = @(x) gradU_2LSHAPE(x);
% divU = @(x) 0*x(:,1)+0*x(:,2);
% gdivU = @(x) 0*[x(:,1),x(:,2)];
% gdivU_1 = @(x) 0*x(:,1)+0*x(:,2);
% gdivU_2 = @(x) 0*x(:,1)+0*x(:,2);
% dgdivU = @(x) 0*x(:,1)+0*x(:,2);
% f = @(x) U(x);
% u0 = @(x) gradULSHAPE(x);
% zero = @(x) 0*x(:,1)+0*x(:,2);
% zero_vec = @(x) 0;

%% Example 5 
% curlv = @(x) [x(:,1).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + x(:,1).^2.*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2, ...
%  - x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 - 2.*x(:,1).*x(:,2).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2];
% U = @(x) curlv(x);
% U_1 = @(x) x(:,1).^2.*(x(:,1) - 1).^2.*(x(:,2) - 1).^2 + x(:,1).^2.*x(:,2).*(2.*x(:,2) - 2).*(x(:,1) - 1).^2;
% U_2 = @(x) - x(:,1).^2.*x(:,2).*(2.*x(:,1) - 2).*(x(:,2) - 1).^2 - 2.*x(:,1).*x(:,2).*(x(:,1) - 1).^2.*(x(:,2) - 1).^2;
% divU = @(x) 0*x(:,1)+0*x(:,2);
% gdivU = @(x) 0*[x(:,1),x(:,2)];
% gdivU_1 = @(x) 0*x(:,1)+0*x(:,2);
% gdivU_2 = @(x) 0*x(:,1)+0*x(:,2);
% dgdivU = @(x) 0*x(:,1)+0*x(:,2);
% f = @(x) U(x);
% u0 = @(x) U(x);
%% Example 6 
% U = @(x) [x(:,1) x(:,2)];
% U_1 = @(x) x(:,1);
% U_2 = @(x) x(:,2);
% divU = @(x) 2 + 0*x(:,1)+0*x(:,2);
% gdivU = @(x) 0*[x(:,1),x(:,2)];
% gdivU_1 = @(x) 0*x(:,1)+0*x(:,2);
% gdivU_2 = @(x) 0*x(:,1)+0*x(:,2);
% dgdivU = @(x) 0*x(:,1)+0*x(:,2);
% f = @(x) U(x);
% u0 = @(x) U(x);
%% Example 7 Lshape u=curlv

% load Meshes/LshapeMesh.mat 
% U = @(x) curlVLSHAPE(x);
% U_1 = @(x) curlV_1LSHAPE(x);
% U_2 = @(x) curlV_2LSHAPE(x);
% divU = @(x) 0*x(:,1)+0*x(:,2);
% gdivU = @(x) 0*[x(:,1),x(:,2)];
% gdivU_1 = @(x) 0*x(:,1)+0*x(:,2);
% gdivU_2 = @(x) 0*x(:,1)+0*x(:,2);
% dgdivU = @(x) 0*x(:,1)+0*x(:,2);
% f = @(x) U(x);
% u0 = @(x) curlVLSHAPE(x);
% zero = @(x) 0*x(:,1)+0*x(:,2);
% zero_vec = @(x) 0;
%% Example 8
U = @(x) [x(:,1).^4,x(:,2).^4];
U_1 = @(x) x(:,1).^4;
U_2 = @(x) x(:,2).^4;
divU = @(x) 4*x(:,1).^3 + 4*x(:,2).^3;
gdivU = @(x) [12*x(:,1).^2,12*x(:,2).^2];
gdivU_1 = @(x) 12*x(:,1).^2;
gdivU_2 = @(x) 12*x(:,2).^2;
f = @(x) U(x)+[24,24];
u0 = @(x) U(x);

%% ansatz space
p = 0;
pU1 = p;
pU2 = p;

% trace Variables
pU1hat = 1;
qU1hat = 1;
pU2hat = 1;
qU2hat = 1;

% test space
qV1 = 3;
qV2 = 3;


%%
%% Other parameters
% marked elements, max number of iterations, and max number of elements in
% mesh 
% algorithm terminates after the mesh contains >= nEmax elements
markedE = [];
maxL = 100;
nEmax = 4000;

%% Set adaptivity parameter \theta
% 0<theta<1 is adaptive refinement using 
% theta = 1 is uniform mesh refinement
theta = 1;
tic;
for j=1:maxL
    % Refine mesh
    [coordinates,elements,boundary] = refineMeshFB(coordinates,elements,boundary,markedE,[]);
     
    [coordinates,elements,boundary,edges,element2edges,...
        boundary2edges] = prepareFBMesh(coordinates,elements,boundary);

    % store some values (number of elements, vertices, etc.)
    nE(j) = size(elements,1);
    nC(j) = size(coordinates,1);
    nB(j) = size(boundary,1);
    nEd(j) = size(edges,1);
    nBEd(j) = size(boundary2edges,1);
    fprintf('Step %d, nE = %d, nC = %d, nB = %d, nEd = %d\n',j,nE(j),nC(j),nB(j),nEd(j));

    
    % dimension of discretization spaces
    dimU1 = nE(j)*(pU1+1)*(pU1+2);
    dimU2 = nE(j)*(pU2+1)*(pU2+2);
    dimU = dimU1+dimU2;
    
    dimU1qhat = nEd(j)*(qU1hat+1);
    dimU1phat = nC(j)+nEd(j)*(pU1hat-1);
    dimU2qhat = nEd(j)*(qU2hat+1);
    dimU2phat = nC(j)+nEd(j)*(pU2hat-1);
    dimUhat = dimU1phat+dimU1qhat+dimU2phat+dimU2qhat;
    
    dimV1 = nE(j)*(qV1+1)*(qV1+2);
    dimV2 = nE(j)*(qV2+1)*(qV2+2);
    dimV = dimV1+dimV2;
    
    % degrees of freedom (approximately)
    dofDPG(j) = dimU + dimUhat;

    %fprintf('j = %d\n    Building B ...\n',j);
    B = buildBMatLap(coordinates,elements,edges,element2edges,...
    pU1,pU2,pU1hat,qU1hat,pU2hat,qU2hat,qV1,qV2);
    %fprintf('    Building V ...\n');
    Vinv = buildVinvMatLap(coordinates,elements,qV1,qV2);
    
    %*** right-hand side
    bVec = buildRHSLap(coordinates,elements,qV1,qV2,f);
    
%     bVec = zeros(dimV,1);
%     x1=coordinates(elements(:,1),:);

%     x2=coordinates(elements(:,2),:);
%     x3=coordinates(elements(:,3),:);
%     baricenter = (x1+x2+x3)/3;
%     f_baricenter = feval(f,baricenter);
%     f_baricenter = f_baricenter';
%     [sigtauTmp,Itmp,Jtmp] = buildBLFsigmaTau(coordinates,elements,0,2,[]);
%     Sigtau = sparse(Itmp,Jtmp,sigtauTmp,dimV1,dimU1);
%     bVec(1:dimV1) = -Sigtau*f_baricenter(:);
    
    %% boundary dof's 
     jdx =[];
     if(pU1hat>=1)
        %B(:,dimU+dimSigma+nC(j)-nB(j)+1:dimU+dimSigma+nC(j)) = [];
        jdx = [jdx;(dimU+nC(j)-nB(j)+1:dimU+nC(j))'];
     end
     if(pU1hat>=2)
        %B(:,dimU+dimSigma+nC(j)-nB(j)+1:dimU+dimSigma+nC(j)) = [];
        jdx = [jdx;(dimU+nC(j)+nEd(j)-nB(j)+1:dimU+nC(j)+nEd(j))'];
     end
    
     kdx = [];
     if (qU1hat>=0)
         kdx = [kdx;(dimU+dimU1phat+nEd(j)-nBEd(j)+1:dimU+dimU1phat+nEd(j))'];
     end
     if (qU1hat>=1)
         kdx = [kdx;(dimU+dimU1phat+2*nEd(j)-nBEd(j)+1:dimU+dimU1phat+2*nEd(j))']; 
     end
     if (qU1hat>=2) 
         kdx = [kdx; (dimU+dimU1phat+3*nEd(j)-nBEd(j)+1:dimU+dimU1phat+3*nEd(j))'];
     end        

    
    
    %% Define system matrix: S = B'V^{-1} B where V is Riesz-matrix
    b = B'*(Vinv*bVec);
    S = B'*(Vinv*B);

    %*** determine indices corresponding to the dof's
    dof = size(S,1);
    freeonodes = setdiff( (1:dof)',[jdx;kdx]);

    %*** initialise solution vector x
    %fprintf('    Solve ...\n');
    x = zeros(dof,1);
    
    %*** Dirichlet boundary values
    x(jdx) = computeL2projBOUC1(coordinates(nC(j)-nB(j)+1:end,:),boundary-nC(j)+nB(j),divU,pU1hat);
    x(kdx) = computeL2projNeumannPp(coordinates(nC(j)-nB(j)+1:end,:),boundary-nC(j)+nB(j),u0,qU1hat);
    %x(jdx) = 0;
    %x(kdx) = 0;
    b = b-S*x;
    
    x(freeonodes) = S(freeonodes,freeonodes)\b(freeonodes);

    nEle(j) = nE(j);

    %*** DPG error estimator
    err(j) = sqrt( (bVec-B*x)'*Vinv*(bVec-B*x) );

    %*** L2 error of u1
    U1H = [x(1:2:dimU1) x(2:2:dimU1)];
    tmpErr1 = computeErrL2(coordinates,elements,U1H,U,pU1);
    errU1(j) = sqrt( sum( sum(tmpErr1,2) ) );
    
        
    %*** L2 error of u2
    U2H = [x(1+dimU1:2:dimU1+dimU2) x(2+dimU1:2:dimU1+dimU2)];
    tmpErr2 = computeErrL2(coordinates,elements,U2H,gdivU,pU2);
    errU2(j) = sqrt( sum( sum(tmpErr2,2) ) );
    
    % *** L2 error of proj U1
    %U1P = computeL2projToPp_Vec(coordinates,elements,U_1,U_2,pU1);
    %tmpL2NormProjU1 = computeErrL2(coordinates,elements,U1P,U,pU1);
    %errU1P(j) = sqrt( sum( sum(tmpL2NormProjU1,2) ) );
    
    % *** L2 error of proj U2
    %U2P = computeL2projToPp_Vec(coordinates,elements,gdivU_1,gdivU_2,pU2);
    %tmpL2NormProjU2 = computeErrL2(coordinates,elements,U2P,gdivU,pU2);
    %errU2P(j) = sqrt( sum( sum(tmpL2NormProjU2,2) ) );
        
    % Relative Error 
    %tmpL2NormU1 = computeErrL2(coordinates,elements,0,U,pU1);
    %normU1(j) = sqrt(sum(sum(tmpL2NormU1,2)));
    %rel_errU1(j) = (errU1(j)/normU1(j)); %Relative Error for U1 dpg solution
    %rel_errU1P(j) = (errU1P(j)/normU1(j)); %Relative Error for U1 projection
    
    %*** Functions
    xU1 = [x(1:2:dimU1) x(2:2:dimU1)];   
    xU2 = [x(dimU1+1:2:dimU1+dimU2) x(dimU1+2:2:dimU1+dimU2)];   
    


    %*** error quantity (= DPG-error estimator)
    est = (bVec-B*x).*(Vinv*(bVec-B*x));


    % H(Gdiv) parts for V1
    est1 = est(1:dimV1);
    est1 = sum( reshape(est1,(qV1+1)*(qV1+2),nE(j)),1)';
    % H(Gdiv) parts for V2
    est2 = est(dimV1+1:dimV1+dimV2);
    est2 = sum( reshape(est2,(qV2+1)*(qV2+2),nE(j)),1)';
   
    % overall volume-elementwise estimator:
    estVol = est1 +est2 ;

    % mark elements
    markedE = markElements(theta,estVol);
    if(theta==1)
        markedE = (1:nE(j))';
    end

    if(nE(j)>=nEmax)
        break;
    end   
end
toc;
%% Plot convergence rates
%close all
minP=min([pU1;pU2;pU1hat;qU1hat;pU2hat;qU2hat])+1;
figure
%loglog(dofDPG,err,'--x',dofDPG,errU1,'--o',dofDPG,errU2,'--d',dofDPG,errU3,'--^',dofDPG,errU4,'--s',...
%    dofDPG,err(1)*(dofDPG/dofDPG(1)).^(-minP/2),'--k',dofDPG,errU1(1)*(dofDPG/dofDPG(1)).^(-(pU1+1)/2),'--k',...
%    dofDPG,errU2(1)*(dofDPG/dofDPG(1)).^(-(pU2+1)/2),'--k',dofDPG,errU3(1)*(dofDPG/dofDPG(1)).^(-(pU3+1)/2),'--k',...
%    dofDPG,errU4(1)*(dofDPG/dofDPG(1)).^(-(pU4+1)/2),'--k','LineWidth',2.5)

loglog(dofDPG,err,'--x',dofDPG,errU1,'--o',dofDPG,errU2,'-d',dofDPG,err(1)*(dofDPG/dofDPG(1)).^(-minP/2),'--k'...
,'LineWidth',2.5);


legend('$\eta$','\boldmath$\|u_1-u_{1,h}\|$','\boldmath$\|u_2-u_{2,h}\|$','$O(h)$','interpreter','latex','Location','southwest');
xlabel('Degrees of Freedom','interpreter','latex')
%xtitle(['Error degree  p = ', num2str(p)]);




%% plot solution component
% for now it only works for polynomial degrees pU = 0 and pU = 1 
%   figure
%   trisurf(elements,coordinates(:,1),coordinates(:,2),x(dimU+dimU1hat+(1:nC(j))))
%   if(p<=2)    
%     figure;  
%     trisurf(elements,coordinates(:,1),coordinates(:,2),U_1(coordinates));
%     title('Exact solution U1(:,1)');
%     
%     figure;  
%     trisurf(elements,coordinates(:,1),coordinates(:,2),U_2(coordinates));
%     title('Exact solution U2(:,1)');
%     
%     figure;
%     myP0plot(coordinates,elements,x(1:2:dimU1),pU1);
%     title('Approximated solution U1(:,1) ');
    
%     figure;
%     myP0plot(coordinates,elements,x(2:2:dimU1),pU1);
%     title('Approximated solution U2(:,1) ');
%     
%    end 