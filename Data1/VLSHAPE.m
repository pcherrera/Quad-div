function out = VLSHAPE(x)

alpha = 2/3;
[t,r]=cart2pol(x(:,1),x(:,2));

%*** evaluate Dirichlet data
out= r.^(alpha).*cos(alpha.*t);

