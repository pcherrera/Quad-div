function grad=curlV_2LSHAPE(x)

grad = zeros(size(x,1),1);
alpha = 2/3;

%*** compute gradient gradu(x)
[t,r]=cart2pol(x(:,1),x(:,2));
grad(:,1)=-cos(t).*(alpha*r.^(alpha-1).*cos(alpha*t))+(1./r).*sin(t).*(-alpha.*sin(alpha*t).*r.^alpha);

