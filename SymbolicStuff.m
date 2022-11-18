% Symbolic stuff 
clc
clear all
%% Start Variables 

syms a1 a2 b1 b2 c1 c2 x y u v a b c d real
syms e1 e2 e3 o1 o2 o3 b11 b12 b21 b22 x y real

eta{1} = x;
eta{3} = y;
eta{2} = 1-x-y;


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

for j =1:9 % degrees p =1 ,2 ,3
    etaHat{j} = eta {j};
end
for j =1:3 % degree p =4
    etaHat{9+j} = eta{10+j};
end

nu{1} = 1;
nu{2} = x;
nu{3} = y;

%B = [b1-a1,c1-a1;b2-a2,c2-a2];
B = [b11,b12;b21,b22];
C = [a b; c d];
Binv = [ b22/(b11*b22 - b12*b21), -b12/(b11*b22 - b12*b21); -b21/(b11*b22 - b12*b21),  b11/(b11*b22 - b12*b21)];
k=0;
for j=1:2:42
    k=k+1;
    ETA{j} = [eta{k} 0]';
    ETA{j+1} = [0 eta{k}]';
end

for j=4:15
    nu{j}=eta{j};
end

for j =1:13
    NU{(j-1)*2+1} = [nu{j} 0]';
    NU{(j-1)*2+2} = [0 nu{j}]';
end

syms t
%sigmaHat{1}=1;
%sigmaHat{2}=2*t-1;
%sigmaHat{3}=4*(1-t)*t;
%sigmaHat{4}=4*(2*t-1)*t*(1-t);

sigmaHat{1}=1;
sigmaHat{2}=t;
sigmaHat{3}=t^2;
sigmaHat{4}=t^3;



r=0;
s=0;
E = [e1,e2,e3]; 
O = [o1,o2,o3];


%% (u,v)
r=0;
p = 0;
q = 3;
dimP = (p+1)*(q+2);
dimQ = (q+1)*(q+2);
for i=1:dimP
    for j=1:dimQ
        r=r+1;
        udv = dot(NU{i},C*ETA{j});
        R = int(int(udv,y,0,1-x),x,0,1);
        R = simplify(R);
        fprintf('M(1,%d)= %s; \n',r,R);
    end
end
%% <u,v.n>
clc;
r=0;
p = 1;
q = 3;
dimP = 3*p;
dimQ = (q+1)*(q+2);
for i=1:dimP
    for j=1:dimQ
        u = matlabFunction(etaHat{i},'Vars',{x,y});
        v1 = matlabFunction(ETA{j}(1,:),'Vars',{x,y});
        v2 = matlabFunction(ETA{j}(2,:),'Vars',{x,y});
        R1 = int(u(t,0)*(v1(t,0)*(a2-b2)-v2(t,0)*(a1-b1)),0,1);
        R2 = int(u(t,1-t)*(v1(t,1-t)*(c2-a2)-v2(t,1-t)*(c1-a1)),0,1);
        R3 = int(u(0,t)*(v1(0,t)*(b2-c2)-v2(0,t)*(b1-c1)),0,1);
        r=r+1;
        fprintf('M(1,%d)= %s; \n',r,R1+R2+R3);
    end
end

%% (grad div u, grad div u) + (u,u) 
clc;
r=0;
q=4;
dimQ = (q+1)*(q+2);
for i=1:dimQ
    for j=1:dimQ
        r=r+1;
        dudv = det(B)*dot(B'\gradient(divergence(B\ETA{i},[x,y]),[x,y]),...
        B'\gradient(divergence(B\ETA{j},[x,y]),[x,y])); 
        uv = det(B)*dot(ETA{i},ETA{j});
        R = int(int(uv+dudv,y,0,1-x),x,0,1); 
        fprintf('M(1,%d)=%s;\n',r,R);
    end
end    
%% (u,grad div v) 
clc
r=0; 
p=0;
q=2;
dimP = (p+1)*(p+2);
dimQ = (q+1)*(q+2);
for i = 1:dimP
    for j = 1:dimQ
        r=r+1;
        udv = dot(NU{i},Binv'*gradient(divergence(Binv*ETA{j},[x,y]),[x,y]));
        R = int(int(udv,y,0,1-x),x,0,1);
        R = simplify(R);
        fprintf('M(1,%d)=%s;\n',r,R)
    end
end
%% < u.n,div v >
clc
r=0;
p=0;
q=4;
dimP = 3*(p+1);
dimQ = (q+1)*(q+2);
for i=1:dimP
    for j = 1:dimQ
        r=r+1;
        u = sigmaHat{fix((i-1)/3)+1};
        dv = matlabFunction(divergence(adjoint(B)*ETA{j},[x,y]),'Vars',{x,y,b11,b12,b21,b22});
        %dv = matlabFunction(divergence(ETA{j},[x,y]),'Vars',{x,y});
        
        if (rem(i-1,3)+1 ==  1)
           R = O(rem(i-1,3)+1)*E(rem(i-1,3)+1)*int(u*dv(t,0,b11,b12,b21,b22),t,0,1);
        end
        if (rem(i-1,3)+1 == 2) 
           R = O(rem(i-1,3)+1)*E(rem(i-1,3)+1)*int(u*dv(0,t,b11,b12,b21,b22),t,0,1);
        end
        if (rem(i-1,3)+1 == 3)
           R = O(rem(i-1,3)+1)*E(rem(i-1,3)+1)*int(u*dv(t,1-t,b11,b12,b21,b22),t,0,1);
        end
        %if (rem(i-1,3)+1 ==  1)
        %   R = O(rem(i-1,3)+1)*E(rem(i-1,3)+1)*int(u*dv(t,0),t,0,1);
        %elseif (rem(i-1,3)+1 == 2) 
        %   R = O(rem(i-1,3)+1)*E(rem(i-1,3)+1)*int(u*dv(0,-t),t,0,1);
        %else 
        %   R = O(rem(i-1,3)+1)*E(rem(i-1,3)+1)*int(u*dv(t,1-t),t,0,1);
        %end
        fprintf('M(1,%d)= %s; \n',r,R);
    end
end
%% <u.n,v> 
r=0;
p=1;
q=3;
dimP = 3*(p+1);
dimQ = (q+1)*(q+2)/2;
for i=1:dimP
    for j=1:dimQ
        r=r+1;
        u = sigmaHat{fix((i-1)/3)+1};
        dv = matlabFunction(eta{j},'Vars',{x,y});
        if (rem(i-1,3)+1 ==  1)
           R = O(rem(i-1,3)+1)*E(rem(i-1,3)+1)*int(u*dv(t,0),t,0,1);
        end
        if (rem(i-1,3)+1 == 2) 
           R = O(rem(i-1,3)+1)*E(rem(i-1,3)+1)*int(u*dv(0,t),t,0,1);
        end
        if (rem(i-1,3)+1 == 3)
           R = O(rem(i-1,3)+1)*E(rem(i-1,3)+1)*int(u*dv(t,1-t),t,0,1);
        end
        
        fprintf('M(1,%d)= %s;\n',r,R);
    end
end
%% Scalling 
clc
syms h
r=0;
p=0;
q=4;
dimP = 3*(p+1);
dimQ = (q+1)*(q+2)/2;

fprintf('funciones test para T1\n');
for i=1:dimQ
    R1 = matlabFunction(eta{i},'Vars',{x,y});
    S = R1(x/h,y/h);
    fprintf('%s \n',S);
end
fprintf('funciones Test para T2\n');
for i=1:dimQ
    R1 = matlabFunction(eta{i},'Vars',{x,y});
    S = simplify(R1((x-h)/h,(y-h)/h));
    fprintf('%s \n',S);
end

