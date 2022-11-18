function grad = curlVLSHAPE(x)

grad = zeros(size(x,1),2);

alpha = 2/3;

%*** compute gradient gradu(x)
[t,r]=cart2pol(x(:,1),x(:,2));
grad(:,1)= sin(t).*(alpha*r.^(alpha-1).*cos(alpha*t))+(1./r).*cos(t).*(-alpha.*sin(alpha*t).*r.^alpha);
grad(:,2)= -cos(t).*(alpha*r.^(alpha-1).*cos(alpha*t))+(1./r).*sin(t).*(-alpha.*sin(alpha*t).*r.^alpha);
